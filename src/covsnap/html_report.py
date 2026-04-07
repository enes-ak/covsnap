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


def _status_color(status: str) -> str:
    return {
        "PASS": "var(--success)",
        "LOW_COVERAGE": "var(--warning)",
        "DROP_OUT": "var(--danger)",
        "UNEVEN": "var(--primary)",
        "LOW_EXON": "var(--secondary)",
    }.get(status, "var(--secondary)")


def _build_html(ctx: ReportContext) -> str:
    """Assemble the full HTML document."""
    parts: list[str] = []
    w = parts.append

    # --- Status counts ---
    status_counts: dict[str, int] = {}
    for r in ctx.results:
        status_counts[r.coverage_status] = status_counts.get(r.coverage_status, 0) + 1
    total = len(ctx.results)

    # --- Serialize target data for JS charts ---
    chart_data = []
    for r in ctx.results:
        chart_data.append({
            "id": r.target_id,
            "mean": round(r.mean_depth, 2),
            "status": r.coverage_status,
            "pct_thresholds": {str(k): round(v, 2) for k, v in sorted(r.pct_thresholds.items())},
        })

    # --- Exon heatmap data ---
    exon_heatmap_data: list[dict[str, Any]] = []
    if ctx.exon_results:
        for gene, exons in ctx.exon_results.items():
            meta_list = (ctx.exon_metadata or {}).get(gene, [])
            for i, er in enumerate(exons):
                meta = meta_list[i] if i < len(meta_list) else {}
                exon_heatmap_data.append({
                    "gene": gene,
                    "exon": meta.get("exon_number", i + 1),
                    "exon_id": meta.get("exon_id", er.target_id),
                    "pct20": round(er.pct_thresholds.get(20, 0.0), 2),
                    "mean": round(er.mean_depth, 2),
                })

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

    # ── Header ──
    w('<header class="header">')
    w('<div class="container">')
    w(f'<h1>covsnap Coverage Report</h1>')
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

    # ── Summary cards ──
    w('<section class="cards">')
    _card(w, "Total Targets", str(total), "primary")
    _card(w, "PASS", str(status_counts.get("PASS", 0)), "success")
    _card(w, "LOW_COVERAGE", str(status_counts.get("LOW_COVERAGE", 0)), "warning")
    _card(w, "DROP_OUT", str(status_counts.get("DROP_OUT", 0)), "danger")
    _card(w, "UNEVEN", str(status_counts.get("UNEVEN", 0)), "primary")
    _card(w, "LOW_EXON", str(status_counts.get("LOW_EXON", 0)), "secondary")
    w("</section>")

    # ── BED Guardrails ──
    if ctx.bed_guardrail_message:
        w('<div class="alert alert-warning">')
        w(f"<strong>BED Guardrail Warning:</strong> {_E(ctx.bed_guardrail_message)}")
        w("</div>")

    # ── Charts ──
    w('<section class="section">')
    w("<h2>Coverage Overview</h2>")
    w('<div class="chart-row">')
    w('<div class="chart-container"><h3>Per-Target Mean Depth</h3><svg id="chart-depth"></svg></div>')
    w('<div class="chart-container"><h3>Threshold Coverage Distribution</h3><svg id="chart-thresholds"></svg></div>')
    w("</div>")
    if exon_heatmap_data:
        w('<div class="chart-container full-width"><h3>Exon-Level Coverage Heatmap (%&ge;20x)</h3><svg id="chart-exon-heatmap"></svg></div>')
    w("</section>")

    # Tooltip div (shared)
    w('<div id="tooltip" class="tooltip"></div>')

    # ── Target summary table ──
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
        w(f"<td data-sort=\"{r.length_bp}\">{r.length_bp:,} bp</td>")
        w(f'<td data-sort="{r.mean_depth:.2f}">{r.mean_depth:.2f}</td>')
        w(f'<td data-sort="{r.median_depth:.0f}">{r.median_depth:.0f}</td>')
        w(f'<td data-sort="{pct20:.2f}">{pct20:.2f}%</td>')
        w(f'<td data-sort="{r.pct_zero:.2f}">{r.pct_zero:.2f}%</td>')
        w(f'<td><span class="badge badge-{_E(r.coverage_status)}">{_E(r.coverage_status)}</span></td>')
        w("</tr>")
    w("</tbody>")
    w("</table>")
    w("</section>")

    # ── Detailed per-target sections ──
    w('<section class="section">')
    w("<h2>Detailed Results</h2>")
    for r in ctx.results:
        pct20 = r.pct_thresholds.get(20, 0.0)
        w(f'<div class="detail-card" id="detail-{_E(r.target_id)}">')
        w(f'<h3>{_E(r.target_id)} &mdash; <span class="badge badge-{_E(r.coverage_status)}">{_E(r.coverage_status)}</span></h3>')
        w(f"<p>{_E(_status_description(r.coverage_status))}</p>")
        w('<div class="metric-grid">')
        _metric(w, "Mean depth", f"{r.mean_depth:.2f}x")
        _metric(w, "Median depth", f"{r.median_depth:.0f}x")
        _metric(w, "Min / Max", f"{r.min_depth} / {r.max_depth}")
        _metric(w, "Std deviation", f"{r.stdev_depth:.2f}")
        _metric(w, "% at 0x", f"{r.pct_zero:.2f}%")
        for t in sorted_thresholds:
            _metric(w, f"% &ge; {t}x", f"{r.pct_thresholds.get(t, 0.0):.2f}%")
        _metric(w, "Low-cov blocks", f"{r.n_lowcov_blocks} blocks, {r.lowcov_total_bp:,} bp")
        w("</div>")
        w(f'<p class="rationale"><strong>Classification rationale:</strong> {_E(_rationale(r, ctx.classify_params))}</p>')
        w("</div>")
    w("</section>")

    # ── Exon-level results ──
    if ctx.exon_results:
        w('<section class="section">')
        w("<h2>Exon-Level Results</h2>")
        for gene, exon_list in ctx.exon_results.items():
            w(f"<h3>{_E(gene)}</h3>")
            meta_list = (ctx.exon_metadata or {}).get(gene, [])
            status_list = (ctx.exon_statuses or {}).get(gene, [])
            w('<table class="data-table">')
            w("<thead><tr><th>Exon</th><th>Exon ID</th><th>Region</th><th>Length</th><th>Mean</th><th>%&ge;20x</th><th>%Zero</th><th>Flag</th></tr></thead>")
            w("<tbody>")
            for i, er in enumerate(exon_list):
                meta = meta_list[i] if i < len(meta_list) else {}
                exon_id = meta.get("exon_id", er.target_id)
                exon_num = meta.get("exon_number", "?")
                status = status_list[i] if i < len(status_list) else "OK"
                pct20 = er.pct_thresholds.get(20, 0.0)
                row_cls = ' class="exon-low"' if status == "LOW_EXON" else ""
                w(f"<tr{row_cls}>")
                w(f"<td>{_E(str(exon_num))}</td>")
                w(f"<td>{_E(str(exon_id))}</td>")
                w(f"<td>{_E(_fmt_region(er.contig, er.start, er.end))}</td>")
                w(f"<td>{er.length_bp:,} bp</td>")
                w(f"<td>{er.mean_depth:.2f}</td>")
                w(f"<td>{pct20:.2f}%</td>")
                w(f"<td>{er.pct_zero:.2f}%</td>")
                w(f"<td><span class=\"badge badge-{'LOW_EXON' if status == 'LOW_EXON' else 'PASS'}\">{_E(status)}</span></td>")
                w("</tr>")
            w("</tbody></table>")
        w("</section>")

    # ── Low-coverage blocks ──
    if ctx.emit_lowcov:
        w('<section class="section">')
        w("<h2>Low-Coverage Blocks</h2>")
        w(f"<p><em>Threshold: depth &lt; {ctx.lowcov_threshold}, min length: {ctx.lowcov_min_len} bp</em></p>")
        any_blocks = any(r.lowcov_blocks for r in ctx.results)
        if any_blocks:
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
        else:
            w("<p>No low-coverage blocks detected.</p>")
        w("</section>")

    # ── Classification heuristics reference ──
    p = ctx.classify_params
    w('<section class="section">')
    w("<h2>Classification Heuristics Reference</h2>")
    w('<table class="data-table">')
    w("<thead><tr><th>Status</th><th>Rule</th></tr></thead>")
    w("<tbody>")
    w(f"<tr><td><span class=\"badge badge-PASS\">PASS</span></td><td>pct_ge_20 &ge; {p.pass_pct_ge_20}% AND pct_zero &le; {p.pass_max_pct_zero}%</td></tr>")
    w(f"<tr><td><span class=\"badge badge-LOW_COVERAGE\">LOW_COVERAGE</span></td><td>pct_ge_20 &lt; {p.pass_pct_ge_20}% (and not DROP_OUT)</td></tr>")
    w(f"<tr><td><span class=\"badge badge-DROP_OUT\">DROP_OUT</span></td><td>pct_zero &gt; {p.dropout_pct_zero}% OR zero-block &ge; {p.dropout_zero_block_bp} bp</td></tr>")
    w(f"<tr><td><span class=\"badge badge-UNEVEN\">UNEVEN</span></td><td>mean_depth &gt; 20 AND CV &gt; {p.uneven_cv}</td></tr>")
    w(f"<tr><td><span class=\"badge badge-LOW_EXON\">LOW_EXON</span></td><td>(exon mode) any exon pct_ge_20 &lt; {p.exon_pct_ge_20}% or pct_zero &gt; {p.exon_max_pct_zero}%</td></tr>")
    w("</tbody></table>")
    w("<p>Evaluation order: DROP_OUT &rarr; UNEVEN &rarr; LOW_EXON &rarr; LOW_COVERAGE &rarr; PASS.</p>")
    w("</section>")

    w("</main>")

    # ── Footer ──
    w('<footer class="footer">')
    w(f'<div class="container">Generated by covsnap {_E(__version__)}</div>')
    w("</footer>")

    # ── JavaScript ──
    w("<script>")
    w(f"const CHART_DATA = {json.dumps(chart_data)};")
    w(f"const THRESHOLDS = {json.dumps(sorted_thresholds)};")
    w(f"const EXON_HEATMAP = {json.dumps(exon_heatmap_data)};")
    w(_JS)
    w("</script>")

    w("</body>")
    w("</html>")
    return "\n".join(parts)


# ---------------------------------------------------------------------------
# HTML helpers
# ---------------------------------------------------------------------------


def _meta_item(w, label: str, value: str) -> None:
    w(f'<div class="meta-item"><span class="meta-label">{label}</span><span class="meta-value">{value}</span></div>')


def _card(w, title: str, value: str, color: str) -> None:
    w(f'<div class="card card-{color}"><div class="card-value">{value}</div><div class="card-title">{title}</div></div>')


def _metric(w, label: str, value: str) -> None:
    w(f'<div class="metric"><span class="metric-label">{label}</span><span class="metric-value">{value}</span></div>')


# ---------------------------------------------------------------------------
# CSS (inline)
# ---------------------------------------------------------------------------

_CSS = """\
:root {
  --bg: #F8FAFC;
  --card: #FFFFFF;
  --primary: #2563EB;
  --secondary: #64748B;
  --accent: #14B8A6;
  --success: #22C55E;
  --warning: #F59E0B;
  --danger: #EF4444;
  --text: #0F172A;
  --border: #E2E8F0;
}
*, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
body {
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
  background: var(--bg);
  color: var(--text);
  line-height: 1.6;
}
.container { max-width: 1200px; margin: 0 auto; padding: 0 24px; }

/* Header */
.header {
  background: linear-gradient(135deg, var(--primary), #1E40AF);
  color: #fff;
  padding: 32px 0;
}
.header h1 { font-size: 1.75rem; margin-bottom: 16px; }
.meta-grid { display: flex; flex-wrap: wrap; gap: 12px 24px; }
.meta-item { display: flex; flex-direction: column; }
.meta-label { font-size: 0.7rem; text-transform: uppercase; letter-spacing: 0.05em; opacity: 0.8; }
.meta-value { font-size: 0.9rem; font-weight: 500; }
.meta-value code { background: rgba(255,255,255,0.15); padding: 1px 6px; border-radius: 4px; font-size: 0.85em; }

/* Cards */
.cards { display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 16px; margin: 24px 0; }
.card {
  background: var(--card);
  border-radius: 10px;
  padding: 20px;
  text-align: center;
  border: 1px solid var(--border);
  box-shadow: 0 1px 3px rgba(0,0,0,0.06);
}
.card-value { font-size: 2rem; font-weight: 700; }
.card-title { font-size: 0.8rem; color: var(--secondary); margin-top: 4px; text-transform: uppercase; letter-spacing: 0.04em; }
.card-primary .card-value { color: var(--primary); }
.card-success .card-value { color: var(--success); }
.card-warning .card-value { color: var(--warning); }
.card-danger .card-value  { color: var(--danger); }
.card-secondary .card-value { color: var(--secondary); }

/* Sections */
.section { margin: 32px 0; }
.section h2 { font-size: 1.35rem; margin-bottom: 16px; border-bottom: 2px solid var(--border); padding-bottom: 8px; }
.section h3 { font-size: 1.1rem; margin: 16px 0 8px; }

/* Alert */
.alert { padding: 12px 16px; border-radius: 8px; margin-bottom: 16px; }
.alert-warning { background: #FEF3C7; border: 1px solid var(--warning); color: #92400E; }

/* Tables */
.data-table { width: 100%; border-collapse: collapse; background: var(--card); border-radius: 8px; overflow: hidden; box-shadow: 0 1px 3px rgba(0,0,0,0.06); }
.data-table th, .data-table td { padding: 10px 14px; text-align: left; border-bottom: 1px solid var(--border); font-size: 0.875rem; }
.data-table th { background: #F1F5F9; font-weight: 600; cursor: pointer; user-select: none; white-space: nowrap; }
.data-table th.sortable:hover { background: #E2E8F0; }
.sort-arrow { font-size: 0.7em; opacity: 0.4; }
.sort-arrow.asc::after { content: " \\25B2"; opacity: 1; }
.sort-arrow.desc::after { content: " \\25BC"; opacity: 1; }
.data-table tbody tr:hover { background: #F8FAFC; }
.exon-low { background: #FEF2F2 !important; }

/* Badges */
.badge { display: inline-block; padding: 2px 10px; border-radius: 12px; font-size: 0.75rem; font-weight: 600; letter-spacing: 0.02em; }
.badge-PASS { background: #DCFCE7; color: #166534; }
.badge-LOW_COVERAGE { background: #FEF9C3; color: #854D0E; }
.badge-DROP_OUT { background: #FEE2E2; color: #991B1B; }
.badge-UNEVEN { background: #DBEAFE; color: #1E40AF; }
.badge-LOW_EXON { background: #F1F5F9; color: #475569; }

/* Filter bar */
.filter-bar { display: flex; gap: 12px; margin-bottom: 12px; flex-wrap: wrap; }
.filter-input, .filter-select {
  padding: 8px 12px;
  border: 1px solid var(--border);
  border-radius: 6px;
  font-size: 0.875rem;
  background: var(--card);
}
.filter-input { flex: 1; min-width: 200px; }
.filter-select { min-width: 160px; }

/* Detail cards */
.detail-card {
  background: var(--card);
  border: 1px solid var(--border);
  border-radius: 10px;
  padding: 20px;
  margin-bottom: 16px;
  box-shadow: 0 1px 3px rgba(0,0,0,0.06);
}
.detail-card h3 { margin-bottom: 8px; }
.metric-grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(180px, 1fr)); gap: 8px 16px; margin: 12px 0; }
.metric { display: flex; justify-content: space-between; padding: 4px 0; border-bottom: 1px dotted var(--border); }
.metric-label { color: var(--secondary); font-size: 0.85rem; }
.metric-value { font-weight: 600; font-size: 0.85rem; }
.rationale { margin-top: 8px; font-size: 0.85rem; color: var(--secondary); }

/* Charts */
.chart-row { display: grid; grid-template-columns: 1fr 1fr; gap: 24px; }
@media (max-width: 800px) { .chart-row { grid-template-columns: 1fr; } }
.chart-container {
  background: var(--card);
  border: 1px solid var(--border);
  border-radius: 10px;
  padding: 16px;
  box-shadow: 0 1px 3px rgba(0,0,0,0.06);
}
.chart-container.full-width { margin-top: 24px; }
.chart-container h3 { font-size: 0.95rem; margin-bottom: 12px; }
.chart-container svg { width: 100%; height: 320px; }
.chart-container.full-width svg { height: 200px; }

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
"""

# ---------------------------------------------------------------------------
# JavaScript (inline)
# ---------------------------------------------------------------------------

_JS = r"""
(function() {
  "use strict";

  // ── Tooltip ──
  const tip = document.getElementById("tooltip");
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

  // ── Status colors ──
  const SC = {
    PASS: "#22C55E", LOW_COVERAGE: "#F59E0B", DROP_OUT: "#EF4444",
    UNEVEN: "#2563EB", LOW_EXON: "#64748B"
  };

  // ── Chart A: Per-target mean depth bar chart ──
  (function() {
    const svg = document.getElementById("chart-depth");
    if (!svg || CHART_DATA.length === 0) return;
    const W = 100, H = 100;
    svg.setAttribute("viewBox", "0 0 " + W + " " + H);
    svg.setAttribute("preserveAspectRatio", "none");
    const n = CHART_DATA.length;
    const maxVal = Math.max(...CHART_DATA.map(d => d.mean), 1);
    const barW = (W - 10) / n;
    const oX = 8;

    // Y-axis label
    for (let i = 0; i <= 4; i++) {
      const y = H - 8 - (H - 16) * i / 4;
      const val = Math.round(maxVal * i / 4);
      const t = document.createElementNS("http://www.w3.org/2000/svg", "text");
      t.setAttribute("x", oX - 1);
      t.setAttribute("y", y + 1);
      t.setAttribute("text-anchor", "end");
      t.setAttribute("font-size", "3");
      t.setAttribute("fill", "#64748B");
      t.textContent = val;
      svg.appendChild(t);
      const line = document.createElementNS("http://www.w3.org/2000/svg", "line");
      line.setAttribute("x1", oX);
      line.setAttribute("x2", W);
      line.setAttribute("y1", y);
      line.setAttribute("y2", y);
      line.setAttribute("stroke", "#E2E8F0");
      line.setAttribute("stroke-width", "0.2");
      svg.appendChild(line);
    }

    CHART_DATA.forEach(function(d, i) {
      const barH = (d.mean / maxVal) * (H - 16);
      const x = oX + i * barW + barW * 0.1;
      const y = H - 8 - barH;
      const r = document.createElementNS("http://www.w3.org/2000/svg", "rect");
      r.setAttribute("x", x);
      r.setAttribute("y", y);
      r.setAttribute("width", barW * 0.8);
      r.setAttribute("height", barH);
      r.setAttribute("fill", SC[d.status] || "#64748B");
      r.setAttribute("rx", "0.5");
      r.addEventListener("mousemove", function(e) {
        showTip(e, "<b>" + d.id + "</b><br>Mean: " + d.mean.toFixed(2) + "x<br>Status: " + d.status);
      });
      r.addEventListener("mouseleave", hideTip);
      svg.appendChild(r);

      // Label
      if (n <= 30) {
        const t = document.createElementNS("http://www.w3.org/2000/svg", "text");
        t.setAttribute("x", x + barW * 0.4);
        t.setAttribute("y", H - 3);
        t.setAttribute("text-anchor", "middle");
        t.setAttribute("font-size", Math.min(3, 80 / n).toFixed(1));
        t.setAttribute("fill", "#64748B");
        t.textContent = d.id.length > 10 ? d.id.slice(0, 9) + "…" : d.id;
        svg.appendChild(t);
      }
    });
  })();

  // ── Chart B: Threshold coverage stacked bar ──
  (function() {
    const svg = document.getElementById("chart-thresholds");
    if (!svg || CHART_DATA.length === 0) return;
    const W = 100, H = 100;
    svg.setAttribute("viewBox", "0 0 " + W + " " + H);
    svg.setAttribute("preserveAspectRatio", "none");
    const n = CHART_DATA.length;
    const barW = (W - 10) / n;
    const oX = 8;
    const colors = ["#DBEAFE", "#93C5FD", "#60A5FA", "#3B82F6", "#2563EB", "#1D4ED8", "#1E3A8A"];

    // Y-axis
    for (let i = 0; i <= 4; i++) {
      const y = H - 8 - (H - 16) * i / 4;
      const t = document.createElementNS("http://www.w3.org/2000/svg", "text");
      t.setAttribute("x", oX - 1);
      t.setAttribute("y", y + 1);
      t.setAttribute("text-anchor", "end");
      t.setAttribute("font-size", "3");
      t.setAttribute("fill", "#64748B");
      t.textContent = (25 * i) + "%";
      svg.appendChild(t);
    }

    CHART_DATA.forEach(function(d, i) {
      const x = oX + i * barW + barW * 0.1;
      const w = barW * 0.8;
      // Stack thresholds from highest to lowest
      const keys = THRESHOLDS.slice().reverse();
      keys.forEach(function(thr, ti) {
        const pct = d.pct_thresholds[String(thr)] || 0;
        const barH = (pct / 100) * (H - 16);
        const y = H - 8 - barH;
        const r = document.createElementNS("http://www.w3.org/2000/svg", "rect");
        r.setAttribute("x", x);
        r.setAttribute("y", y);
        r.setAttribute("width", w);
        r.setAttribute("height", barH);
        r.setAttribute("fill", colors[Math.min(ti, colors.length - 1)]);
        r.setAttribute("rx", "0.3");
        r.addEventListener("mousemove", function(e) {
          showTip(e, "<b>" + d.id + "</b><br>≥" + thr + "x: " + pct.toFixed(1) + "%");
        });
        r.addEventListener("mouseleave", hideTip);
        svg.appendChild(r);
      });

      if (n <= 30) {
        const t = document.createElementNS("http://www.w3.org/2000/svg", "text");
        t.setAttribute("x", x + w / 2);
        t.setAttribute("y", H - 3);
        t.setAttribute("text-anchor", "middle");
        t.setAttribute("font-size", Math.min(3, 80 / n).toFixed(1));
        t.setAttribute("fill", "#64748B");
        t.textContent = d.id.length > 10 ? d.id.slice(0, 9) + "…" : d.id;
        svg.appendChild(t);
      }
    });
  })();

  // ── Chart C: Exon heatmap ──
  (function() {
    const svg = document.getElementById("chart-exon-heatmap");
    if (!svg || EXON_HEATMAP.length === 0) return;
    const W = 100, H = 100;
    svg.setAttribute("viewBox", "0 0 " + W + " " + H);
    const n = EXON_HEATMAP.length;
    const cellW = Math.min((W - 6) / n, 12);
    const cellH = 20;
    const oX = 3;
    const oY = (H - cellH) / 2;

    EXON_HEATMAP.forEach(function(d, i) {
      const pct = d.pct20;
      // Color gradient: red (0%) -> yellow (50%) -> green (100%)
      var r, g, b;
      if (pct < 50) {
        r = 239; g = Math.round(68 + (pct / 50) * 180); b = 68;
      } else {
        r = Math.round(239 - ((pct - 50) / 50) * 205); g = 197; b = Math.round(68 + ((pct - 50) / 50) * 26);
      }
      const fill = "rgb(" + r + "," + g + "," + b + ")";
      const rect = document.createElementNS("http://www.w3.org/2000/svg", "rect");
      rect.setAttribute("x", oX + i * cellW);
      rect.setAttribute("y", oY);
      rect.setAttribute("width", cellW - 0.5);
      rect.setAttribute("height", cellH);
      rect.setAttribute("fill", fill);
      rect.setAttribute("rx", "1");
      rect.addEventListener("mousemove", function(e) {
        showTip(e, "<b>" + d.gene + " exon " + d.exon + "</b><br>%≥20x: " + pct.toFixed(1) + "%<br>Mean: " + d.mean.toFixed(1) + "x");
      });
      rect.addEventListener("mouseleave", hideTip);
      svg.appendChild(rect);

      // Label below
      if (n <= 40) {
        const t = document.createElementNS("http://www.w3.org/2000/svg", "text");
        t.setAttribute("x", oX + i * cellW + (cellW - 0.5) / 2);
        t.setAttribute("y", oY + cellH + 6);
        t.setAttribute("text-anchor", "middle");
        t.setAttribute("font-size", Math.min(3, 60 / n).toFixed(1));
        t.setAttribute("fill", "#64748B");
        t.textContent = d.exon;
        svg.appendChild(t);
      }
    });

    // Legend
    var legendX = oX + n * cellW + 4;
    if (legendX > W - 20) legendX = W - 20;
    [0, 50, 100].forEach(function(v, li) {
      var cr, cg, cb;
      if (v < 50) { cr = 239; cg = 68; cb = 68; }
      else if (v === 50) { cr = 239; cg = 248; cb = 68; }
      else { cr = 34; cg = 197; cb = 94; }
      const rect = document.createElementNS("http://www.w3.org/2000/svg", "rect");
      rect.setAttribute("x", legendX);
      rect.setAttribute("y", oY + li * 7);
      rect.setAttribute("width", 5);
      rect.setAttribute("height", 5);
      rect.setAttribute("fill", "rgb(" + cr + "," + cg + "," + cb + ")");
      svg.appendChild(rect);
      const t = document.createElementNS("http://www.w3.org/2000/svg", "text");
      t.setAttribute("x", legendX + 7);
      t.setAttribute("y", oY + li * 7 + 4);
      t.setAttribute("font-size", "3");
      t.setAttribute("fill", "#64748B");
      t.textContent = v + "%";
      svg.appendChild(t);
    });
  })();

  // ── Sortable table ──
  (function() {
    const table = document.getElementById("target-table");
    if (!table) return;
    const headers = table.querySelectorAll("th.sortable");
    let currentCol = -1, ascending = true;

    headers.forEach(function(th) {
      th.addEventListener("click", function() {
        const col = parseInt(th.dataset.col, 10);
        if (col === currentCol) { ascending = !ascending; }
        else { currentCol = col; ascending = true; }

        // Update arrows
        headers.forEach(function(h) {
          h.querySelector(".sort-arrow").className = "sort-arrow";
        });
        th.querySelector(".sort-arrow").className = "sort-arrow " + (ascending ? "asc" : "desc");

        const tbody = table.querySelector("tbody");
        const rows = Array.from(tbody.querySelectorAll("tr"));
        rows.sort(function(a, b) {
          const cellA = a.children[col];
          const cellB = b.children[col];
          let va = cellA.dataset.sort !== undefined ? cellA.dataset.sort : cellA.textContent;
          let vb = cellB.dataset.sort !== undefined ? cellB.dataset.sort : cellB.textContent;
          const na = parseFloat(va), nb = parseFloat(vb);
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

  // ── Filter ──
  (function() {
    const table = document.getElementById("target-table");
    const textInput = document.getElementById("filter-text");
    const statusSelect = document.getElementById("filter-status");
    if (!table || !textInput || !statusSelect) return;

    function applyFilter() {
      const text = textInput.value.toLowerCase();
      const status = statusSelect.value;
      const rows = table.querySelectorAll("tbody tr");
      rows.forEach(function(row) {
        const cells = row.querySelectorAll("td");
        const rowText = row.textContent.toLowerCase();
        const rowStatus = cells[cells.length - 1] ? cells[cells.length - 1].textContent.trim() : "";
        const matchText = !text || rowText.indexOf(text) !== -1;
        const matchStatus = !status || rowStatus === status;
        row.style.display = (matchText && matchStatus) ? "" : "none";
      });
    }
    textInput.addEventListener("input", applyFilter);
    statusSelect.addEventListener("change", applyFilter);
  })();

})();
"""
