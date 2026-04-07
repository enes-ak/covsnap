"""Self-contained interactive HTML report generator for covsnap.

Produces a single HTML file with all CSS and JavaScript inline.
No external dependencies required -- the report works when opened
from any location on disk.
"""

from __future__ import annotations

import html
import json
from datetime import datetime, timezone
from typing import Any

from covsnap import ANNOTATION_VERSION, BUILD, __version__
from covsnap.metrics import TargetResult
from covsnap.report import ClassifyParams, ReportContext, _rationale, _status_description


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def write_html_report(path: str, ctx: ReportContext) -> None:
    """Write an interactive HTML coverage report to *path*."""
    if ctx.run_date is None:
        ctx.run_date = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    doc = _build_html(ctx)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(doc)


# ---------------------------------------------------------------------------
# Internal builder
# ---------------------------------------------------------------------------

_E = html.escape  # shorthand


def _format_date(raw: str) -> str:
    """Convert ISO date like '2026-04-07T11:48:18Z' to '07/04/2026 11:48:18'."""
    for fmt in ("%Y-%m-%dT%H:%M:%SZ", "%Y-%m-%dT%H:%M:%S %Z", "%Y-%m-%d %H:%M:%S %Z", "%Y-%m-%d %H:%M:%S"):
        try:
            dt = datetime.strptime(raw, fmt)
            return dt.strftime("%d/%m/%Y %H:%M:%S")
        except ValueError:
            continue
    return raw


def _fmt_region(contig: str, start: int, end: int) -> str:
    return f"{contig}:{start + 1:,}-{end:,}"


def _status_css_class(status: str) -> str:
    return {
        "PASS": "pass",
        "LOW_COVERAGE": "low",
        "DROP_OUT": "drop",
        "UNEVEN": "uneven",
        "LOW_EXON": "lexon",
    }.get(status, "low")


# Severity ranking for "worst status" computation
_STATUS_SEVERITY = {
    "PASS": 0,
    "LOW_EXON": 1,
    "UNEVEN": 2,
    "LOW_COVERAGE": 3,
    "DROP_OUT": 4,
}


def _build_html(ctx: ReportContext) -> str:
    """Assemble the full HTML document."""
    parts: list[str] = []
    w = parts.append

    # --- Compute summary statistics ---
    total = len(ctx.results)
    n_pass = sum(1 for r in ctx.results if r.coverage_status == "PASS")
    pass_rate = (n_pass / total * 100) if total > 0 else 0.0

    total_bp = sum(r.length_bp for r in ctx.results)
    if total_bp > 0:
        overall_mean = sum(r.mean_depth * r.length_bp for r in ctx.results) / total_bp
        overall_pct20 = sum(r.pct_thresholds.get(20, 0.0) * r.length_bp for r in ctx.results) / total_bp
    else:
        overall_mean = 0.0
        overall_pct20 = 0.0

    worst_status = "PASS"
    for r in ctx.results:
        if _STATUS_SEVERITY.get(r.coverage_status, 0) > _STATUS_SEVERITY.get(worst_status, 0):
            worst_status = r.coverage_status

    # Count low exons
    n_low_exons = 0
    if ctx.exon_statuses:
        for gene, statuses in ctx.exon_statuses.items():
            n_low_exons += sum(1 for s in statuses if s == "LOW_EXON")

    # --- Exon bar chart data (sorted by exon_number ascending) ---
    exon_bar_data: list[dict[str, Any]] = []
    if ctx.exon_results:
        for gene, exons in ctx.exon_results.items():
            meta_list = (ctx.exon_metadata or {}).get(gene, [])
            status_list = (ctx.exon_statuses or {}).get(gene, [])
            for i, er in enumerate(exons):
                meta = meta_list[i] if i < len(meta_list) else {}
                status = status_list[i] if i < len(status_list) else "OK"
                exon_bar_data.append({
                    "gene": gene,
                    "exon": meta.get("exon_number", i + 1),
                    "exon_id": meta.get("exon_id", er.target_id),
                    "pct20": round(er.pct_thresholds.get(20, 0.0), 2),
                    "mean": round(er.mean_depth, 2),
                    "status": status,
                })
        # Sort by exon number ascending
        exon_bar_data.sort(key=lambda d: (d["gene"], d["exon"]))

    sorted_thresholds = sorted(ctx.thresholds)

    # ================================================================
    # HTML document
    # ================================================================
    w("<!DOCTYPE html>")
    w('<html lang="en">')
    w("<head>")
    w('<meta charset="utf-8">')
    w('<meta name="viewport" content="width=device-width, initial-scale=1">')
    w(f"<title>covsnap Report — {_E(ctx.sample_name)}</title>")
    w("<style>")
    w(_CSS)
    w("</style>")
    w("</head>")
    w("<body>")

    # -- Header --
    w('<header class="header">')
    w('<div class="container">')
    w('<h1>covsnap Coverage Report</h1>')
    w('<div class="meta-grid">')
    _meta_item(w, "Tool version", f"covsnap {__version__}")
    _meta_item(w, "Date", _E(_format_date(ctx.run_date or "")))
    _meta_item(w, "Engine", f"{_E(ctx.engine_used)} {_E(ctx.engine_version)}")
    _meta_item(w, "Annotation", f"GENCODE v44 ({BUILD})")
    _meta_item(w, "Input", f"<code>{_E(ctx.bam_path)}</code>")
    _meta_item(w, "Sample", _E(ctx.sample_name))
    _meta_item(w, "Contig style", f"{_E(ctx.contig_style)}-prefixed (auto-detected)")
    if ctx.reference:
        _meta_item(w, "Reference", f"<code>{_E(ctx.reference)}</code>")
    if ctx.bed_path:
        _meta_item(w, "Input BED", f"<code>{_E(ctx.bed_path)}</code>")
    w("</div>")  # meta-grid
    w("</div>")  # container
    w("</header>")

    w('<main class="container">')

    # -- Summary cards --
    w('<section class="summary-cards">')
    _summary_card(w, "Overall Mean Depth", f"{overall_mean:.1f}x", "primary")
    _summary_card(w, "Overall %&ge;20x", f"{overall_pct20:.1f}%", "primary")
    _summary_card(w, "Targets Analyzed", str(total), "secondary")
    _summary_card(w, "PASS Rate", f"{pass_rate:.0f}%",
                  "pass" if pass_rate >= 90 else ("low" if pass_rate >= 50 else "drop"))
    _summary_card(w, "Worst Status", worst_status, _status_css_class(worst_status))
    if ctx.exon_results:
        _summary_card(w, "Low Exons", str(n_low_exons),
                      "pass" if n_low_exons == 0 else "drop")
    w("</section>")

    # -- BED Guardrails --
    if ctx.bed_guardrail_message:
        w('<div class="alert alert-warning">')
        w(f"<strong>BED Guardrail Warning:</strong> {_E(ctx.bed_guardrail_message)}")
        w("</div>")

    # -- Exon Coverage Bar Chart (only if exon data exists) --
    if exon_bar_data:
        w('<section class="section">')
        w("<h2>Exon Coverage Overview</h2>")
        w('<p class="section-desc">Horizontal bars show %&ge;20x coverage per exon. '
          'Color indicates quality: teal (adequate) to amber (warning) to red (poor).</p>')
        w('<div class="exon-chart-container" id="exon-chart-container"></div>')
        w("</section>")

    # Tooltip div (shared)
    w('<div id="tooltip" class="tooltip"></div>')

    # -- Target summary table --
    w('<section class="section">')
    w("<h2>Target Summary</h2>")
    w('<div class="filter-bar">')
    w('<input type="text" id="filter-text" placeholder="Search targets..." class="filter-input">')
    w('<select id="filter-status" class="filter-select">')
    w('<option value="">All statuses</option>')
    for s in ["PASS", "LOW_COVERAGE", "DROP_OUT", "UNEVEN", "LOW_EXON"]:
        w(f'<option value="{s}">{s}</option>')
    w("</select>")
    w("</div>")

    w('<div class="table-wrap">')
    w('<table class="data-table" id="target-table">')
    w("<thead><tr>")
    headers = ["Target", "Region (1-based)", "Length", "Mean Depth", "Median", "%&ge;20x", "%Zero", "Status"]
    for i, h in enumerate(headers):
        w(f'<th data-col="{i}" class="sortable">{h} <span class="sort-arrow"></span></th>')
    w("</tr></thead>")
    w("<tbody>")
    for r in ctx.results:
        pct20 = r.pct_thresholds.get(20, 0.0)
        w("<tr>")
        w(f"<td>{_E(r.target_id)}</td>")
        w(f"<td>{_E(_fmt_region(r.contig, r.start, r.end))}</td>")
        w(f'<td data-sort="{r.length_bp}">{r.length_bp:,} bp</td>')
        w(f'<td data-sort="{r.mean_depth:.2f}">{r.mean_depth:.2f}</td>')
        w(f'<td data-sort="{r.median_depth:.0f}">{r.median_depth:.0f}</td>')
        w(f'<td data-sort="{pct20:.2f}">{pct20:.2f}%</td>')
        w(f'<td data-sort="{r.pct_zero:.2f}">{r.pct_zero:.2f}%</td>')
        w(f'<td><span class="badge badge-{_status_css_class(r.coverage_status)}">{_E(r.coverage_status)}</span></td>')
        w("</tr>")
    w("</tbody>")
    w("</table>")
    w("</div>")
    w("</section>")

    # -- Detailed per-target sections (accordion) --
    w('<section class="section">')
    w("<h2>Detailed Results</h2>")
    for r in ctx.results:
        is_problem = r.coverage_status != "PASS"
        open_cls = " open" if is_problem else ""
        w(f'<div class="accordion{open_cls}" id="detail-{_E(r.target_id)}">')
        w(f'<div class="accordion-header" onclick="this.parentElement.classList.toggle(\'open\')">')
        w(f'<span class="accordion-title">{_E(r.target_id)}</span>')
        w(f'<span class="badge badge-{_status_css_class(r.coverage_status)}">{_E(r.coverage_status)}</span>')
        w(f'<span class="accordion-summary">{_E(_status_description(r.coverage_status))}</span>')
        w('<span class="accordion-chevron"></span>')
        w("</div>")
        w('<div class="accordion-body">')
        w('<table class="metric-table">')
        w("<tbody>")
        _metric_row(w, "Region", _E(_fmt_region(r.contig, r.start, r.end)))
        _metric_row(w, "Length", f"{r.length_bp:,} bp")
        _metric_row(w, "Mean depth", f"{r.mean_depth:.2f}x")
        _metric_row(w, "Median depth", f"{r.median_depth:.0f}x")
        _metric_row(w, "Min / Max depth", f"{r.min_depth} / {r.max_depth}")
        _metric_row(w, "Std deviation", f"{r.stdev_depth:.2f}")
        _metric_row(w, "% at 0x", f"{r.pct_zero:.2f}%")
        for t in sorted_thresholds:
            _metric_row(w, f"% &ge; {t}x", f"{r.pct_thresholds.get(t, 0.0):.2f}%")
        _metric_row(w, "Low-cov blocks", f"{r.n_lowcov_blocks} blocks, {r.lowcov_total_bp:,} bp")
        w("</tbody>")
        w("</table>")
        w(f'<div class="rationale"><strong>Classification rationale:</strong> {_E(_rationale(r, ctx.classify_params))}</div>')
        w("</div>")  # accordion-body
        w("</div>")  # accordion
    w("</section>")

    # -- Gene-level results (region mode) --
    if ctx.gene_results:
        w('<section class="section">')
        w("<h2>Genes in Region</h2>")
        w(f"<p>{len(ctx.gene_results)} gene(s) found overlapping the queried region.</p>")
        w('<div class="table-wrap">')
        w('<table class="data-table" id="geneTable">')
        w("<thead><tr>")
        w("<th>Gene</th><th>Region</th><th>Length</th><th>Mean</th><th>Median</th>")
        w("<th>%&ge;20x</th><th>%Zero</th><th>Status</th>")
        w("</tr></thead><tbody>")
        for i, gr in enumerate(ctx.gene_results):
            meta = ctx.gene_metadata[i] if ctx.gene_metadata and i < len(ctx.gene_metadata) else {}
            pct20 = gr.pct_thresholds.get(20, 0.0)
            gene_type = meta.get("gene_type", "")
            strand = meta.get("strand", "")
            w("<tr>")
            w(f"<td><strong>{_E(gr.target_id)}</strong>")
            if gene_type or strand:
                w(f' <small class="gene-annotation">({_E(gene_type)}, {_E(strand)})</small>')
            w("</td>")
            w(f"<td>{_E(_fmt_region(gr.contig, gr.start, gr.end))}</td>")
            w(f"<td>{gr.length_bp:,} bp</td>")
            w(f"<td>{gr.mean_depth:.2f}</td>")
            w(f"<td>{gr.median_depth:.0f}</td>")
            w(f"<td>{pct20:.2f}%</td>")
            w(f"<td>{gr.pct_zero:.2f}%</td>")
            w(f'<td><span class="badge badge-{_status_css_class(gr.coverage_status)}">{_E(gr.coverage_status)}</span></td>')
            w("</tr>")
        w("</tbody></table>")
        w("</div>")
        w("</section>")

    # -- Exon-level results table --
    if ctx.exon_results:
        w('<section class="section">')
        w("<h2>Exon-Level Results</h2>")
        for gene, exon_list in ctx.exon_results.items():
            w(f"<h3>{_E(gene)}</h3>")
            meta_list = (ctx.exon_metadata or {}).get(gene, [])
            status_list = (ctx.exon_statuses or {}).get(gene, [])

            # Build sortable list by exon_number
            exon_rows = []
            for i, er in enumerate(exon_list):
                meta = meta_list[i] if i < len(meta_list) else {}
                exon_num = meta.get("exon_number", i + 1)
                exon_id = meta.get("exon_id", er.target_id)
                status = status_list[i] if i < len(status_list) else "OK"
                exon_rows.append((exon_num, exon_id, er, status))
            exon_rows.sort(key=lambda x: (x[0] if isinstance(x[0], int) else 0))

            w('<div class="table-wrap">')
            w('<table class="data-table exon-table">')
            w("<thead><tr><th>Exon</th><th>Exon ID</th><th>Region</th><th>Length</th>"
              "<th>Mean</th><th>%&ge;20x</th><th>%Zero</th><th>Flag</th></tr></thead>")
            w("<tbody>")
            for exon_num, exon_id, er, status in exon_rows:
                pct20 = er.pct_thresholds.get(20, 0.0)
                row_cls = ' class="exon-problem"' if status == "LOW_EXON" else ""
                w(f"<tr{row_cls}>")
                w(f"<td>{_E(str(exon_num))}</td>")
                w(f"<td>{_E(str(exon_id))}</td>")
                w(f"<td>{_E(_fmt_region(er.contig, er.start, er.end))}</td>")
                w(f"<td>{er.length_bp:,} bp</td>")
                w(f"<td>{er.mean_depth:.2f}</td>")
                w(f"<td>{pct20:.2f}%</td>")
                w(f"<td>{er.pct_zero:.2f}%</td>")
                w(f'<td><span class="badge badge-{"lexon" if status == "LOW_EXON" else "pass"}">{_E(status)}</span></td>')
                w("</tr>")
            w("</tbody></table>")
            w("</div>")
        w("</section>")

    # -- Low-coverage blocks --
    if ctx.emit_lowcov:
        w('<section class="section">')
        w("<h2>Low-Coverage Blocks</h2>")
        w(f"<p><em>Threshold: depth &lt; {ctx.lowcov_threshold}, min length: {ctx.lowcov_min_len} bp</em></p>")
        any_blocks = any(r.lowcov_blocks for r in ctx.results)
        if any_blocks:
            w('<div class="table-wrap">')
            w('<table class="data-table">')
            w("<thead><tr><th>Target</th><th>Block (1-based)</th><th>Length</th><th>Mean Depth</th></tr></thead>")
            w("<tbody>")
            for r in ctx.results:
                for b in r.lowcov_blocks:
                    w("<tr>")
                    w(f"<td>{_E(r.target_id)}</td>")
                    w(f"<td>{_E(_fmt_region(r.contig, b.start, b.end))}</td>")
                    w(f"<td>{b.length:,} bp</td>")
                    w(f"<td>{b.mean_depth:.1f}x</td>")
                    w("</tr>")
            w("</tbody></table>")
            w("</div>")
        else:
            w("<p>No low-coverage blocks detected.</p>")
        w("</section>")

    # -- Classification heuristics reference --
    p = ctx.classify_params
    w('<section class="section">')
    w("<h2>Classification Heuristics Reference</h2>")
    w('<div class="table-wrap">')
    w('<table class="data-table">')
    w("<thead><tr><th>Status</th><th>Rule</th></tr></thead>")
    w("<tbody>")
    w(f'<tr><td><span class="badge badge-pass">PASS</span></td>'
      f"<td>pct_ge_20 &ge; {p.pass_pct_ge_20}% AND pct_zero &le; {p.pass_max_pct_zero}%</td></tr>")
    w(f'<tr><td><span class="badge badge-low">LOW_COVERAGE</span></td>'
      f"<td>pct_ge_20 &lt; {p.pass_pct_ge_20}% (and not DROP_OUT)</td></tr>")
    w(f'<tr><td><span class="badge badge-drop">DROP_OUT</span></td>'
      f"<td>pct_zero &gt; {p.dropout_pct_zero}% OR zero-block &ge; {p.dropout_zero_block_bp} bp</td></tr>")
    w(f'<tr><td><span class="badge badge-uneven">UNEVEN</span></td>'
      f"<td>mean_depth &gt; 20 AND CV &gt; {p.uneven_cv}</td></tr>")
    w(f'<tr><td><span class="badge badge-lexon">LOW_EXON</span></td>'
      f"<td>(exon mode) any exon pct_ge_20 &lt; {p.exon_pct_ge_20}% or pct_zero &gt; {p.exon_max_pct_zero}%</td></tr>")
    w("</tbody></table>")
    w("</div>")
    w("<p>Evaluation order: DROP_OUT &rarr; UNEVEN &rarr; LOW_EXON &rarr; LOW_COVERAGE &rarr; PASS.</p>")
    w("</section>")

    # -- Glossary / Quick Reference --
    w('<section class="section glossary">')
    w("<h2>Glossary / Quick Reference</h2>")
    w('<dl class="glossary-list">')
    glossary = [
        ("PASS", "Coverage is adequate. &ge;95% of bases have &ge;20x depth."),
        ("LOW_COVERAGE", "Insufficient depth. &lt;95% of bases reach 20x."),
        ("DROP_OUT", "Complete loss of coverage in significant portions (&gt;5% zero-depth or &ge;500bp zero block)."),
        ("UNEVEN", "High mean depth but highly variable (CV &gt;1.0). May indicate amplification bias."),
        ("LOW_EXON", "One or more individual exons have critically low coverage."),
        ("Mean Depth", "Average number of reads covering each base position."),
        ("Median Depth", "Middle value of per-base depths (more robust than mean)."),
        ("%&ge;20x", "Percentage of bases covered by at least 20 reads (clinical adequacy threshold)."),
        ("%Zero", "Percentage of bases with no coverage at all."),
        ("CV (Coefficient of Variation)", "Standard deviation divided by mean. Higher = more uneven."),
    ]
    for term, defn in glossary:
        w(f"<dt>{term}</dt>")
        w(f"<dd>{defn}</dd>")
    w("</dl>")
    w("</section>")

    w("</main>")

    # -- Footer --
    w('<footer class="footer">')
    w(f'<div class="container">Generated by covsnap {_E(__version__)}</div>')
    w("</footer>")

    # -- JavaScript --
    w("<script>")
    w(f"const EXON_BAR_DATA = {json.dumps(exon_bar_data)};")
    w(f"const THRESHOLDS = {json.dumps(sorted_thresholds)};")
    w(_JS)
    w("</script>")

    w("</body>")
    w("</html>")
    return "\n".join(parts)


# ---------------------------------------------------------------------------
# HTML helpers
# ---------------------------------------------------------------------------


def _meta_item(w, label: str, value: str) -> None:
    w(f'<div class="meta-item"><span class="meta-label">{label}</span>'
      f'<span class="meta-value">{value}</span></div>')


def _summary_card(w, title: str, value: str, color: str) -> None:
    w(f'<div class="summary-card sc-{color}">'
      f'<div class="sc-value">{value}</div>'
      f'<div class="sc-title">{title}</div></div>')


def _metric_row(w, label: str, value: str) -> None:
    w(f"<tr><td class=\"metric-label\">{label}</td><td class=\"metric-value\">{value}</td></tr>")


# ---------------------------------------------------------------------------
# CSS (inline)
# ---------------------------------------------------------------------------

_CSS = """\
:root {
  --bg: #F7FAFA;
  --card: #FFFFFF;
  --primary: #0D9488;
  --primary-light: #CCFBF1;
  --primary-dark: #0F766E;
  --secondary: #64748B;
  --text: #0F172A;
  --text-light: #475569;
  --border: #E2E8F0;
  --border-light: #F1F5F9;
  --status-pass: #0D9488;
  --status-pass-bg: #CCFBF1;
  --status-low: #D97706;
  --status-low-bg: #FEF3C7;
  --status-drop: #DC2626;
  --status-drop-bg: #FEE2E2;
  --status-uneven: #2563EB;
  --status-uneven-bg: #DBEAFE;
  --status-lexon: #7C3AED;
  --status-lexon-bg: #EDE9FE;
}
*, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
body {
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
  background: var(--bg);
  color: var(--text);
  font-size: 15px;
  line-height: 1.65;
}
.container { max-width: 1200px; margin: 0 auto; padding: 0 24px; }

/* Header */
.header {
  background: linear-gradient(135deg, var(--primary), var(--primary-dark));
  color: #fff;
  padding: 32px 0;
}
.header h1 { font-size: 1.75rem; margin-bottom: 16px; font-weight: 700; }
.meta-grid { display: flex; flex-wrap: wrap; gap: 12px 28px; }
.meta-item { display: flex; flex-direction: column; }
.meta-label { font-size: 0.7rem; text-transform: uppercase; letter-spacing: 0.06em; opacity: 0.8; }
.meta-value { font-size: 0.9rem; font-weight: 500; }
.meta-value code { background: rgba(255,255,255,0.15); padding: 1px 6px; border-radius: 4px; font-size: 0.85em; }

/* Summary cards */
.summary-cards {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(170px, 1fr));
  gap: 16px;
  margin: 28px 0;
}
.summary-card {
  background: var(--card);
  border-radius: 10px;
  padding: 18px 20px;
  border: 1px solid var(--border);
  border-left: 4px solid var(--secondary);
  box-shadow: 0 1px 3px rgba(0,0,0,0.05);
}
.sc-value { font-size: 1.5rem; font-weight: 700; color: var(--text); }
.sc-title { font-size: 0.78rem; color: var(--text-light); margin-top: 2px; text-transform: uppercase; letter-spacing: 0.04em; }
.sc-primary  { border-left-color: var(--primary); }
.sc-secondary { border-left-color: var(--secondary); }
.sc-pass     { border-left-color: var(--status-pass); }
.sc-low      { border-left-color: var(--status-low); }
.sc-drop     { border-left-color: var(--status-drop); }
.sc-uneven   { border-left-color: var(--status-uneven); }
.sc-lexon    { border-left-color: var(--status-lexon); }

/* Sections */
.section { margin: 36px 0; }
.section h2 {
  font-size: 1.3rem;
  margin-bottom: 14px;
  padding-bottom: 8px;
  border-bottom: 2px solid var(--border);
  color: var(--text);
}
.section h3 { font-size: 1.1rem; margin: 18px 0 10px; color: var(--text); }
.section-desc { color: var(--text-light); font-size: 0.9rem; margin-bottom: 16px; }

/* Alert */
.alert { padding: 12px 16px; border-radius: 8px; margin-bottom: 16px; font-size: 0.9rem; }
.alert-warning { background: var(--status-low-bg); border: 1px solid var(--status-low); color: #92400E; }

/* Tables */
.table-wrap { overflow-x: auto; }
.data-table {
  width: 100%;
  border-collapse: collapse;
  background: var(--card);
  border-radius: 8px;
  overflow: hidden;
  box-shadow: 0 1px 3px rgba(0,0,0,0.05);
}
.data-table th, .data-table td {
  padding: 11px 16px;
  text-align: left;
  border-bottom: 1px solid var(--border);
  font-size: 0.875rem;
}
.data-table thead th {
  background: var(--border-light);
  font-weight: 600;
  cursor: pointer;
  user-select: none;
  white-space: nowrap;
  position: sticky;
  top: 0;
  z-index: 2;
}
.data-table th.sortable:hover { background: var(--border); }
.sort-arrow { font-size: 0.7em; opacity: 0.4; }
.sort-arrow.asc::after { content: " \\25B2"; opacity: 1; }
.sort-arrow.desc::after { content: " \\25BC"; opacity: 1; }
.data-table tbody tr:nth-child(even) { background: var(--border-light); }
.data-table tbody tr:hover { background: #E0F2F1; }

/* Exon table - problem row */
.exon-problem { border-left: 3px solid var(--status-lexon) !important; }

/* Badges */
.badge {
  display: inline-block;
  padding: 2px 10px;
  border-radius: 12px;
  font-size: 0.75rem;
  font-weight: 600;
  letter-spacing: 0.02em;
}
.badge-pass   { background: var(--status-pass-bg);   color: #115E59; }
.badge-low    { background: var(--status-low-bg);    color: #92400E; }
.badge-drop   { background: var(--status-drop-bg);   color: #991B1B; }
.badge-uneven { background: var(--status-uneven-bg); color: #1E40AF; }
.badge-lexon  { background: var(--status-lexon-bg);  color: #5B21B6; }

/* Filter bar */
.filter-bar { display: flex; gap: 12px; margin-bottom: 14px; flex-wrap: wrap; }
.filter-input, .filter-select {
  padding: 8px 14px;
  border: 1px solid var(--border);
  border-radius: 6px;
  font-size: 0.875rem;
  background: var(--card);
  transition: border-color 0.15s;
}
.filter-input:focus, .filter-select:focus { outline: none; border-color: var(--primary); }
.filter-input { flex: 1; min-width: 200px; }
.filter-select { min-width: 160px; }

/* Accordion detail cards */
.accordion {
  background: var(--card);
  border: 1px solid var(--border);
  border-radius: 10px;
  margin-bottom: 12px;
  box-shadow: 0 1px 3px rgba(0,0,0,0.04);
  overflow: hidden;
}
.accordion-header {
  display: flex;
  align-items: center;
  gap: 12px;
  padding: 14px 20px;
  cursor: pointer;
  user-select: none;
  transition: background 0.15s;
}
.accordion-header:hover { background: var(--border-light); }
.accordion-title { font-weight: 600; font-size: 1rem; }
.accordion-summary { flex: 1; color: var(--text-light); font-size: 0.85rem; }
.accordion-chevron {
  width: 20px;
  height: 20px;
  flex-shrink: 0;
  position: relative;
}
.accordion-chevron::before {
  content: "";
  display: block;
  width: 8px; height: 8px;
  border-right: 2px solid var(--secondary);
  border-bottom: 2px solid var(--secondary);
  transform: rotate(-45deg);
  position: absolute;
  top: 4px; left: 6px;
  transition: transform 0.2s;
}
.accordion.open .accordion-chevron::before { transform: rotate(45deg); top: 2px; }
.accordion-body {
  max-height: 0;
  overflow: hidden;
  transition: max-height 0.3s ease;
}
.accordion.open .accordion-body {
  max-height: 2000px;
}
.accordion-body > * { padding: 0 20px; }
.accordion-body > *:last-child { padding-bottom: 20px; }

/* Metric table (2-col inside accordion) */
.metric-table { width: 100%; border-collapse: collapse; margin: 4px 0 12px; }
.metric-table td { padding: 6px 16px; font-size: 0.875rem; border-bottom: 1px dotted var(--border); }
.metric-table .metric-label { color: var(--text-light); width: 45%; }
.metric-table .metric-value { font-weight: 600; }

/* Rationale */
.rationale {
  margin: 10px 0 0;
  padding: 10px 14px;
  font-size: 0.85rem;
  color: var(--text-light);
  background: var(--border-light);
  border-radius: 6px;
  border-left: 3px solid var(--border);
}

/* Gene annotation */
.gene-annotation { color: var(--secondary); }

/* Exon bar chart */
.exon-chart-container {
  background: var(--card);
  border: 1px solid var(--border);
  border-radius: 10px;
  padding: 20px;
  box-shadow: 0 1px 3px rgba(0,0,0,0.05);
  overflow-x: auto;
}
.exon-bar-row {
  display: flex;
  align-items: center;
  margin-bottom: 4px;
  font-size: 0.82rem;
}
.exon-bar-label {
  width: 140px;
  flex-shrink: 0;
  text-align: right;
  padding-right: 12px;
  color: var(--text-light);
  font-weight: 500;
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
}
.exon-bar-track {
  flex: 1;
  height: 22px;
  background: var(--border-light);
  border-radius: 4px;
  position: relative;
  min-width: 100px;
}
.exon-bar-fill {
  height: 100%;
  border-radius: 4px;
  transition: width 0.3s;
  min-width: 1px;
}
.exon-bar-pct {
  width: 55px;
  flex-shrink: 0;
  text-align: right;
  padding-left: 8px;
  font-weight: 600;
  font-size: 0.82rem;
}

/* Glossary */
.glossary { margin-bottom: 32px; }
.glossary-list { margin: 0; padding: 0; }
.glossary-list dt {
  font-weight: 600;
  color: var(--text);
  margin-top: 12px;
  font-size: 0.9rem;
}
.glossary-list dd {
  color: var(--text-light);
  margin-left: 0;
  padding-left: 16px;
  font-size: 0.88rem;
  border-left: 2px solid var(--border);
  margin-top: 2px;
}

/* Tooltip */
.tooltip {
  position: fixed;
  background: var(--text);
  color: #fff;
  padding: 8px 12px;
  border-radius: 6px;
  font-size: 0.8rem;
  pointer-events: none;
  opacity: 0;
  transition: opacity 0.15s;
  z-index: 100;
  max-width: 280px;
  line-height: 1.5;
}

/* Footer */
.footer {
  text-align: center;
  padding: 24px 0;
  margin-top: 40px;
  border-top: 1px solid var(--border);
  color: var(--secondary);
  font-size: 0.85rem;
}

/* Responsive */
@media (max-width: 700px) {
  .summary-cards { grid-template-columns: repeat(2, 1fr); }
  .meta-grid { gap: 8px 16px; }
  .exon-bar-label { width: 90px; font-size: 0.75rem; }
  .accordion-summary { display: none; }
}

/* Print */
@media print {
  .header { background: var(--primary) !important; -webkit-print-color-adjust: exact; print-color-adjust: exact; }
  .accordion-body { max-height: none !important; }
  .accordion-chevron { display: none; }
  .filter-bar { display: none; }
  .tooltip { display: none; }
  body { font-size: 12px; }
  .summary-card, .accordion, .data-table { break-inside: avoid; }
}
"""

# ---------------------------------------------------------------------------
# JavaScript (inline)
# ---------------------------------------------------------------------------

_JS = r"""
(function() {
  "use strict";

  // -- Tooltip --
  var tip = document.getElementById("tooltip");
  function showTip(evt, html) {
    tip.innerHTML = html;
    tip.style.opacity = "1";
    moveTip(evt);
  }
  function moveTip(evt) {
    tip.style.left = (evt.clientX + 12) + "px";
    tip.style.top  = (evt.clientY + 12) + "px";
  }
  function hideTip() { tip.style.opacity = "0"; }

  // -- Exon horizontal bar chart --
  (function() {
    var container = document.getElementById("exon-chart-container");
    if (!container || EXON_BAR_DATA.length === 0) return;

    EXON_BAR_DATA.forEach(function(d) {
      var row = document.createElement("div");
      row.className = "exon-bar-row";

      var label = document.createElement("div");
      label.className = "exon-bar-label";
      label.textContent = d.gene + " E" + d.exon;
      label.title = d.exon_id;
      row.appendChild(label);

      var track = document.createElement("div");
      track.className = "exon-bar-track";

      var fill = document.createElement("div");
      fill.className = "exon-bar-fill";
      var pct = Math.min(d.pct20, 100);
      fill.style.width = pct + "%";

      // Color: teal (good) -> amber (warning) -> red (poor)
      var color;
      if (pct >= 95) {
        color = "#0D9488";      // teal
      } else if (pct >= 80) {
        color = "#0EA5E9";      // sky — transition
      } else if (pct >= 50) {
        color = "#D97706";      // amber
      } else {
        color = "#DC2626";      // red
      }
      fill.style.background = color;

      track.appendChild(fill);
      row.appendChild(track);

      var pctLabel = document.createElement("div");
      pctLabel.className = "exon-bar-pct";
      pctLabel.textContent = d.pct20.toFixed(1) + "%";
      pctLabel.style.color = color;
      row.appendChild(pctLabel);

      // Tooltip on hover
      row.addEventListener("mousemove", function(e) {
        showTip(e, "<b>" + d.gene + " exon " + d.exon + "</b><br>%\u226520x: " + d.pct20.toFixed(1) + "%<br>Mean: " + d.mean.toFixed(1) + "x");
      });
      row.addEventListener("mouseleave", hideTip);

      container.appendChild(row);
    });
  })();

  // -- Sortable table --
  (function() {
    var table = document.getElementById("target-table");
    if (!table) return;
    var headers = table.querySelectorAll("th.sortable");
    var currentCol = -1, ascending = true;

    headers.forEach(function(th) {
      th.addEventListener("click", function() {
        var col = parseInt(th.dataset.col, 10);
        if (col === currentCol) { ascending = !ascending; }
        else { currentCol = col; ascending = true; }

        headers.forEach(function(h) {
          h.querySelector(".sort-arrow").className = "sort-arrow";
        });
        th.querySelector(".sort-arrow").className = "sort-arrow " + (ascending ? "asc" : "desc");

        var tbody = table.querySelector("tbody");
        var rows = Array.from(tbody.querySelectorAll("tr"));
        rows.sort(function(a, b) {
          var cellA = a.children[col];
          var cellB = b.children[col];
          var va = cellA.dataset.sort !== undefined ? cellA.dataset.sort : cellA.textContent;
          var vb = cellB.dataset.sort !== undefined ? cellB.dataset.sort : cellB.textContent;
          var na = parseFloat(va), nb = parseFloat(vb);
          if (!isNaN(na) && !isNaN(nb)) {
            return ascending ? na - nb : nb - na;
          }
          va = va.toLowerCase(); vb = vb.toLowerCase();
          if (va < vb) return ascending ? -1 : 1;
          if (va > vb) return ascending ? 1 : -1;
          return 0;
        });
        rows.forEach(function(r) { tbody.appendChild(r); });
      });
    });
  })();

  // -- Filter --
  (function() {
    var table = document.getElementById("target-table");
    var textInput = document.getElementById("filter-text");
    var statusSelect = document.getElementById("filter-status");
    if (!table || !textInput || !statusSelect) return;

    function applyFilter() {
      var text = textInput.value.toLowerCase();
      var status = statusSelect.value;
      var rows = table.querySelectorAll("tbody tr");
      rows.forEach(function(row) {
        var cells = row.querySelectorAll("td");
        var rowText = row.textContent.toLowerCase();
        var rowStatus = cells[cells.length - 1] ? cells[cells.length - 1].textContent.trim() : "";
        var matchText = !text || rowText.indexOf(text) !== -1;
        var matchStatus = !status || rowStatus === status;
        row.style.display = (matchText && matchStatus) ? "" : "none";
      });
    }
    textInput.addEventListener("input", applyFilter);
    statusSelect.addEventListener("change", applyFilter);
  })();

})();
"""
