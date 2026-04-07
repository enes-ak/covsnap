# mosdepth Single-Pass + HTML Report Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Eliminate mosdepth's redundant second pass for ~50% speedup, and add a self-contained interactive HTML report output.

**Architecture:** Feature 1 simplifies `_run_mosdepth` to a single mosdepth invocation (per-base only), computing all metrics — mean, thresholds, min/max/median/stdev, lowcov — from `TargetAccumulator` which already supports this. Feature 2 adds `html_report.py` with inline CSS/JS/SVG charts, wired via `--html-out` CLI flag.

**Tech Stack:** Python 3.9+, pysam, mosdepth, vanilla HTML/CSS/JS (no external libraries)

---

## File Structure

| Action | File | Responsibility |
|--------|------|----------------|
| Modify | `src/covsnap/engines.py` | Remove double mosdepth pass, simplify to single per-base run |
| Create | `src/covsnap/html_report.py` | Self-contained HTML report generator |
| Modify | `src/covsnap/cli.py` | Add `--html-out` flag, wire `write_html_report` |
| Modify | `tests/test_engines.py` | Verify single-pass mosdepth still produces correct results |
| Create | `tests/test_html_report.py` | Test HTML report generation |
| Modify | `tests/test_cli_smoke.py` | Add CLI smoke test for `--html-out` |

---

## Task 1: Simplify mosdepth to single pass

**Files:**
- Modify: `src/covsnap/engines.py:288-518`
- Test: `tests/test_engines.py`

The key insight: `TargetAccumulator` already computes thresholds (`_threshold_counts`), mean (`_sum/_n`), min/max, median, stdev, and lowcov from streaming per-base data. The first mosdepth pass (thresholds + regions) is redundant.

- [ ] **Step 1: Write a test that verifies single-pass produces correct metrics**

In `tests/test_engines.py`, add a test inside `TestMosdepthEngine` that checks all metric fields are populated correctly from a single mosdepth run:

```python
def test_single_pass_metrics_complete(self, synthetic_bam, tmp_output_dir):
    """Verify mosdepth produces all metrics (thresholds, mean, etc.) from single pass."""
    results = compute_depth(
        bam_path=synthetic_bam,
        regions=[("chr17", 43044294, 43050000, "test_region")],
        engine="mosdepth",
        thresholds=[1, 5, 10, 20, 30, 50, 100],
        threads=2,
        tmp_dir=str(tmp_output_dir),
    )
    assert len(results) == 1
    r = results[0]
    # All threshold percentages must be present
    for t in [1, 5, 10, 20, 30, 50, 100]:
        assert t in r.pct_thresholds, f"Missing threshold {t}"
        assert 0 <= r.pct_thresholds[t] <= 100
    # Mean must be computed (not zero placeholder for region with reads)
    assert isinstance(r.mean_depth, float)
    # Min/max/median/stdev must be set
    assert r.min_depth >= 0
    assert r.max_depth >= r.min_depth
    assert isinstance(r.median_depth, float)
    assert isinstance(r.stdev_depth, float)
    assert r.stdev_depth >= 0
```

- [ ] **Step 2: Run the test to verify it passes with current double-pass code**

Run: `cd /home/enes/Desktop/denemeler/covinspector && python -m pytest tests/test_engines.py::TestMosdepthEngine::test_single_pass_metrics_complete -v`

Expected: PASS (current code already produces these metrics, just inefficiently)

- [ ] **Step 3: Rewrite `_run_mosdepth` to single pass**

Replace `_run_mosdepth` in `src/covsnap/engines.py` (lines 288-367). The new version runs mosdepth once with per-base output, then streams the per-base file through `_stream_mosdepth_perbase` (which already exists and already computes everything via `TargetAccumulator`):

```python
def _run_mosdepth(
    bam_path: str,
    targets: list[_TargetSpec],
    thresholds: list[int],
    lowcov_threshold: int,
    lowcov_min_len: int,
    threads: int,
    tmp_dir: str,
    reference: Optional[str],
) -> list[TargetResult]:
    """Run mosdepth once and compute all metrics from per-base output."""
    bed_path = _write_tmp_bed(targets, tmp_dir)
    prefix = os.path.join(tmp_dir, "covsnap_mosdepth")

    cmd = [
        "mosdepth",
        "--by", bed_path,
        "--threads", str(threads),
        prefix,
        bam_path,
    ]
    if reference:
        cmd.extend(["--fasta", reference])

    logger.info("Running mosdepth: %s", " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"mosdepth failed (exit {proc.returncode}): {proc.stderr}")

    perbase_file = f"{prefix}.per-base.bed.gz"
    return _stream_mosdepth_perbase(
        perbase_file, targets, thresholds, lowcov_threshold, lowcov_min_len,
    )
```

- [ ] **Step 4: Remove dead functions `_parse_mosdepth_thresholds` and `_parse_mosdepth_regions`**

Delete the following functions from `src/covsnap/engines.py`:
- `_parse_mosdepth_thresholds` (lines 370-413)
- `_parse_mosdepth_regions` (lines 418-438)

- [ ] **Step 5: Run all engine tests**

Run: `cd /home/enes/Desktop/denemeler/covinspector && python -m pytest tests/test_engines.py -v`

Expected: All tests PASS (mosdepth tests may skip if mosdepth not installed)

- [ ] **Step 6: Run full test suite to check nothing else broke**

Run: `cd /home/enes/Desktop/denemeler/covinspector && python -m pytest tests/ -v`

Expected: All tests PASS

- [ ] **Step 7: Commit**

```bash
git add src/covsnap/engines.py tests/test_engines.py
git commit -m "perf: simplify mosdepth to single pass

Remove redundant first mosdepth invocation (thresholds + regions).
All metrics are now computed from a single per-base pass via
TargetAccumulator, which already supports streaming threshold
computation. This eliminates ~50% of mosdepth runtime."
```

---

## Task 2: Create HTML report module — layout and CSS

**Files:**
- Create: `src/covsnap/html_report.py`
- Test: `tests/test_html_report.py`

This task builds the HTML skeleton, CSS styles, and the main `write_html_report` function that accepts a `ReportContext` (same as the markdown report).

- [ ] **Step 1: Write a test for basic HTML report generation**

Create `tests/test_html_report.py`:

```python
"""Tests for HTML report generation."""

from covsnap.html_report import write_html_report
from covsnap.metrics import TargetResult
from covsnap.report import ClassifyParams, ReportContext


def _make_result(**overrides) -> TargetResult:
    defaults = dict(
        target_id="TEST_GENE",
        contig="chr1",
        start=0,
        end=1000,
        length_bp=1000,
        mean_depth=100.0,
        median_depth=95.0,
        min_depth=5,
        max_depth=300,
        stdev_depth=30.0,
        pct_zero=0.5,
        pct_thresholds={1: 99.5, 5: 99.0, 10: 98.5, 20: 97.0, 30: 95.0, 50: 80.0, 100: 40.0},
        n_lowcov_blocks=0,
        lowcov_total_bp=0,
        lowcov_blocks=[],
        coverage_status="PASS",
    )
    defaults.update(overrides)
    return TargetResult(**defaults)


def _make_context(results=None, **overrides) -> ReportContext:
    if results is None:
        results = [_make_result()]
    defaults = dict(
        results=results,
        engine_used="samtools",
        engine_version="1.19",
        bam_path="/data/sample.bam",
        sample_name="SAMPLE_01",
        contig_style="chr",
        thresholds=[1, 5, 10, 20, 30, 50, 100],
        run_date="2026-04-07T00:00:00Z",
    )
    defaults.update(overrides)
    return ReportContext(**defaults)


class TestHTMLReportBasic:
    def test_produces_valid_html(self, tmp_path):
        path = str(tmp_path / "report.html")
        ctx = _make_context()
        write_html_report(path, ctx)

        with open(path) as f:
            html = f.read()
        assert html.startswith("<!DOCTYPE html>")
        assert "</html>" in html
        assert "<style>" in html
        assert "</style>" in html

    def test_contains_metadata(self, tmp_path):
        path = str(tmp_path / "report.html")
        ctx = _make_context()
        write_html_report(path, ctx)

        with open(path) as f:
            html = f.read()
        assert "SAMPLE_01" in html
        assert "samtools" in html
        assert "covsnap" in html.lower()

    def test_contains_target_data(self, tmp_path):
        path = str(tmp_path / "report.html")
        ctx = _make_context()
        write_html_report(path, ctx)

        with open(path) as f:
            html = f.read()
        assert "TEST_GENE" in html
        assert "PASS" in html

    def test_self_contained_no_external_links(self, tmp_path):
        path = str(tmp_path / "report.html")
        ctx = _make_context()
        write_html_report(path, ctx)

        with open(path) as f:
            html = f.read()
        # No external CSS/JS references
        assert 'href="http' not in html
        assert 'src="http' not in html
        assert '<link rel="stylesheet"' not in html

    def test_multiple_targets(self, tmp_path):
        path = str(tmp_path / "report.html")
        results = [
            _make_result(target_id="BRCA1", coverage_status="PASS"),
            _make_result(target_id="TP53", coverage_status="LOW_COVERAGE",
                         pct_thresholds={1: 99.0, 5: 95.0, 10: 90.0, 20: 80.0, 30: 70.0, 50: 50.0, 100: 10.0}),
        ]
        ctx = _make_context(results=results)
        write_html_report(path, ctx)

        with open(path) as f:
            html = f.read()
        assert "BRCA1" in html
        assert "TP53" in html
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd /home/enes/Desktop/denemeler/covinspector && python -m pytest tests/test_html_report.py -v`

Expected: FAIL — `ModuleNotFoundError: No module named 'covsnap.html_report'`

- [ ] **Step 3: Create the HTML report module**

Create `src/covsnap/html_report.py` with the full implementation. The module must:

1. Accept a `ReportContext` (same as markdown report)
2. Generate self-contained HTML with inline `<style>` and `<script>` tags
3. Use the specified color palette
4. Include: header metadata, summary cards, target summary table (sortable), detailed per-target sections
5. No external dependencies

```python
"""Self-contained interactive HTML report for covsnap."""

from __future__ import annotations

import html as html_mod
import json
from datetime import datetime, timezone
from typing import Any, Optional

from covsnap import ANNOTATION_VERSION, BUILD, __version__
from covsnap.metrics import TargetResult
from covsnap.report import ReportContext, _rationale, _status_description


def write_html_report(path: str, ctx: ReportContext) -> None:
    """Write an interactive HTML coverage report."""
    if ctx.run_date is None:
        ctx.run_date = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    targets_json = _build_targets_json(ctx)
    exon_json = _build_exon_json(ctx)

    parts: list[str] = []
    parts.append(_HTML_HEAD)
    parts.append(_build_css())
    parts.append("</style></head><body>")
    parts.append(_build_header(ctx))
    parts.append(_build_summary_cards(ctx))
    parts.append(_build_target_table(ctx))
    parts.append(_build_detail_sections(ctx))
    if ctx.exon_results:
        parts.append(_build_exon_sections(ctx))
    if ctx.emit_lowcov:
        parts.append(_build_lowcov_section(ctx))
    parts.append(_build_heuristics_reference(ctx))
    parts.append(_build_charts_container())
    parts.append(f"<script>const TARGETS={targets_json};const EXONS={exon_json};</script>")
    parts.append(_build_js())
    parts.append("</body></html>")

    with open(path, "w") as fh:
        fh.write("\n".join(parts))


# ---------------------------------------------------------------------------
# HTML skeleton
# ---------------------------------------------------------------------------

_HTML_HEAD = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>covsnap Coverage Report</title>
<style>"""


def _build_css() -> str:
    return """
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
  --shadow: 0 1px 3px rgba(0,0,0,0.08), 0 1px 2px rgba(0,0,0,0.06);
  --radius: 8px;
}

*, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }

body {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
  background: var(--bg);
  color: var(--text);
  line-height: 1.6;
  padding: 2rem;
  max-width: 1400px;
  margin: 0 auto;
}

h1 { font-size: 1.75rem; font-weight: 700; margin-bottom: 0.25rem; }
h2 { font-size: 1.25rem; font-weight: 600; margin: 2rem 0 1rem; color: var(--primary); border-bottom: 2px solid var(--border); padding-bottom: 0.5rem; }
h3 { font-size: 1.1rem; font-weight: 600; margin: 1rem 0 0.5rem; }

.header { margin-bottom: 2rem; }
.header p { color: var(--secondary); font-size: 0.9rem; }

.meta-grid {
  display: grid;
  grid-template-columns: repeat(auto-fill, minmax(250px, 1fr));
  gap: 0.5rem;
  background: var(--card);
  padding: 1rem;
  border-radius: var(--radius);
  box-shadow: var(--shadow);
  margin-top: 1rem;
  font-size: 0.875rem;
}
.meta-grid .label { color: var(--secondary); font-weight: 500; }
.meta-grid .value { font-weight: 600; }

/* Summary cards */
.cards { display: grid; grid-template-columns: repeat(auto-fill, minmax(160px, 1fr)); gap: 1rem; margin: 1.5rem 0; }
.card {
  background: var(--card);
  border-radius: var(--radius);
  padding: 1.25rem;
  box-shadow: var(--shadow);
  text-align: center;
  border-top: 3px solid var(--border);
}
.card .card-value { font-size: 2rem; font-weight: 700; }
.card .card-label { font-size: 0.8rem; color: var(--secondary); text-transform: uppercase; letter-spacing: 0.5px; margin-top: 0.25rem; }
.card.pass { border-top-color: var(--success); }
.card.pass .card-value { color: var(--success); }
.card.fail { border-top-color: var(--danger); }
.card.fail .card-value { color: var(--danger); }
.card.warn { border-top-color: var(--warning); }
.card.warn .card-value { color: var(--warning); }
.card.info { border-top-color: var(--primary); }
.card.info .card-value { color: var(--primary); }
.card.accent { border-top-color: var(--accent); }
.card.accent .card-value { color: var(--accent); }

/* Tables */
.table-wrap {
  background: var(--card);
  border-radius: var(--radius);
  box-shadow: var(--shadow);
  overflow-x: auto;
  margin: 1rem 0;
}
table { width: 100%; border-collapse: collapse; font-size: 0.85rem; }
thead { background: var(--bg); }
th {
  padding: 0.75rem 1rem;
  text-align: left;
  font-weight: 600;
  color: var(--secondary);
  white-space: nowrap;
  cursor: pointer;
  user-select: none;
  position: relative;
}
th:hover { color: var(--primary); }
th .sort-arrow { font-size: 0.7rem; margin-left: 4px; opacity: 0.4; }
th .sort-arrow.active { opacity: 1; color: var(--primary); }
td { padding: 0.6rem 1rem; border-top: 1px solid var(--border); }
tr:hover td { background: #F1F5F9; }

/* Filter */
.filter-bar {
  display: flex;
  gap: 0.75rem;
  align-items: center;
  margin: 1rem 0;
  flex-wrap: wrap;
}
.filter-bar input[type="text"] {
  padding: 0.5rem 0.75rem;
  border: 1px solid var(--border);
  border-radius: var(--radius);
  font-size: 0.85rem;
  width: 220px;
}
.filter-bar select {
  padding: 0.5rem 0.75rem;
  border: 1px solid var(--border);
  border-radius: var(--radius);
  font-size: 0.85rem;
  background: var(--card);
}

/* Status badges */
.badge {
  display: inline-block;
  padding: 0.2rem 0.6rem;
  border-radius: 9999px;
  font-size: 0.75rem;
  font-weight: 600;
  letter-spacing: 0.3px;
}
.badge-PASS { background: #DCFCE7; color: #166534; }
.badge-LOW_COVERAGE { background: #FEF9C3; color: #854D0E; }
.badge-DROP_OUT { background: #FEE2E2; color: #991B1B; }
.badge-UNEVEN { background: #E0E7FF; color: #3730A3; }
.badge-LOW_EXON { background: #FFE4E6; color: #9F1239; }
.badge-OK { background: #DCFCE7; color: #166534; }

/* Detail panels */
.detail-panel {
  background: var(--card);
  border-radius: var(--radius);
  box-shadow: var(--shadow);
  padding: 1.25rem;
  margin: 1rem 0;
}
.detail-panel .metric-grid {
  display: grid;
  grid-template-columns: repeat(auto-fill, minmax(180px, 1fr));
  gap: 0.75rem;
}
.metric-item { padding: 0.5rem; background: var(--bg); border-radius: 6px; }
.metric-item .metric-label { font-size: 0.75rem; color: var(--secondary); }
.metric-item .metric-value { font-size: 1.1rem; font-weight: 600; }
.rationale { margin-top: 0.75rem; padding: 0.75rem; background: var(--bg); border-radius: 6px; font-size: 0.85rem; color: var(--secondary); border-left: 3px solid var(--primary); }

/* Charts */
.chart-section { margin: 2rem 0; }
.chart-container {
  background: var(--card);
  border-radius: var(--radius);
  box-shadow: var(--shadow);
  padding: 1.5rem;
  margin: 1rem 0;
}
svg text { font-family: inherit; }
.chart-bar { transition: opacity 0.2s; }
.chart-bar:hover { opacity: 0.8; }
.chart-tooltip {
  position: absolute;
  background: var(--text);
  color: white;
  padding: 0.4rem 0.6rem;
  border-radius: 4px;
  font-size: 0.75rem;
  pointer-events: none;
  opacity: 0;
  transition: opacity 0.15s;
  z-index: 100;
}

/* Heatmap */
.heatmap-cell {
  stroke: var(--card);
  stroke-width: 1;
}

/* Heuristics table */
.heuristics { font-size: 0.85rem; }
.heuristics td:first-child { font-weight: 600; white-space: nowrap; }

/* Footer */
.footer { margin-top: 3rem; padding-top: 1rem; border-top: 1px solid var(--border); font-size: 0.8rem; color: var(--secondary); text-align: center; }

@media print {
  body { padding: 0; }
  .filter-bar { display: none; }
  .card, .table-wrap, .detail-panel, .chart-container { box-shadow: none; border: 1px solid var(--border); }
}
"""


# ---------------------------------------------------------------------------
# Section builders
# ---------------------------------------------------------------------------

def _esc(text: Any) -> str:
    return html_mod.escape(str(text))


def _fmt_region(contig: str, start: int, end: int) -> str:
    return f"{contig}:{start + 1:,}-{end:,}"


def _build_header(ctx: ReportContext) -> str:
    meta_items = [
        ("Tool version", f"covsnap {__version__}"),
        ("Date", ctx.run_date),
        ("Engine", f"{ctx.engine_used} {ctx.engine_version}"),
        ("Annotation", f"GENCODE v44 ({BUILD})"),
        ("Input", ctx.bam_path),
        ("Sample", ctx.sample_name),
        ("Contig style", f"{ctx.contig_style}-prefixed"),
    ]
    if ctx.reference:
        meta_items.append(("Reference", ctx.reference))
    if ctx.bed_path:
        meta_items.append(("Input BED", ctx.bed_path))

    grid = "".join(
        f'<div><span class="label">{_esc(label)}</span><br><span class="value">{_esc(val)}</span></div>'
        for label, val in meta_items
    )
    return f"""
<div class="header">
  <h1>covsnap Coverage Report</h1>
  <p>Coverage inspector for targeted sequencing QC</p>
  <div class="meta-grid">{grid}</div>
</div>"""


def _build_summary_cards(ctx: ReportContext) -> str:
    total = len(ctx.results)
    counts: dict[str, int] = {}
    for r in ctx.results:
        counts[r.coverage_status] = counts.get(r.coverage_status, 0) + 1

    n_pass = counts.get("PASS", 0)
    n_low = counts.get("LOW_COVERAGE", 0)
    n_drop = counts.get("DROP_OUT", 0)
    n_uneven = counts.get("UNEVEN", 0)
    n_lexon = counts.get("LOW_EXON", 0)

    cards = [
        ("info", str(total), "Total Targets"),
        ("pass", str(n_pass), "PASS"),
        ("warn", str(n_low), "Low Coverage"),
        ("fail", str(n_drop), "Drop Out"),
        ("accent", str(n_uneven), "Uneven"),
    ]
    if n_lexon > 0:
        cards.append(("fail", str(n_lexon), "Low Exon"))

    items = "".join(
        f'<div class="card {cls}"><div class="card-value">{val}</div><div class="card-label">{label}</div></div>'
        for cls, val, label in cards
    )
    return f'<div class="cards">{items}</div>'


def _build_target_table(ctx: ReportContext) -> str:
    sorted_thresholds = sorted(ctx.thresholds)

    ths = "".join(
        f'<th data-sort="num">{_esc(col)} <span class="sort-arrow">&#9650;</span></th>'
        for col in ["Mean", "Median", "Min", "Max", "StdDev", "%Zero"]
        + [f"%>={t}x" for t in sorted_thresholds]
    )

    rows_html = ""
    for r in ctx.results:
        region = _fmt_region(r.contig, r.start, r.end)
        badge = f'<span class="badge badge-{_esc(r.coverage_status)}">{_esc(r.coverage_status)}</span>'
        tds = "".join(f"<td>{v}</td>" for v in [
            f"{r.mean_depth:.2f}", f"{r.median_depth:.0f}", str(r.min_depth), str(r.max_depth),
            f"{r.stdev_depth:.2f}", f"{r.pct_zero:.2f}%",
        ] + [f"{r.pct_thresholds.get(t, 0.0):.2f}%" for t in sorted_thresholds])

        rows_html += f"""<tr>
  <td><strong>{_esc(r.target_id)}</strong></td>
  <td>{_esc(region)}</td>
  <td>{r.length_bp:,}</td>
  {tds}
  <td>{badge}</td>
</tr>"""

    return f"""
<h2>Target Summary</h2>
<div class="filter-bar">
  <input type="text" id="tableFilter" placeholder="Filter targets..." oninput="filterTable()">
  <select id="statusFilter" onchange="filterTable()">
    <option value="">All statuses</option>
    <option value="PASS">PASS</option>
    <option value="LOW_COVERAGE">LOW_COVERAGE</option>
    <option value="DROP_OUT">DROP_OUT</option>
    <option value="UNEVEN">UNEVEN</option>
    <option value="LOW_EXON">LOW_EXON</option>
  </select>
</div>
<div class="table-wrap">
<table id="targetTable">
  <thead><tr>
    <th data-sort="text">Target <span class="sort-arrow">&#9650;</span></th>
    <th data-sort="text">Region <span class="sort-arrow">&#9650;</span></th>
    <th data-sort="num">Length <span class="sort-arrow">&#9650;</span></th>
    {ths}
    <th data-sort="text">Status <span class="sort-arrow">&#9650;</span></th>
  </tr></thead>
  <tbody>{rows_html}</tbody>
</table>
</div>"""


def _build_detail_sections(ctx: ReportContext) -> str:
    sorted_thresholds = sorted(ctx.thresholds)
    sections = ['<h2>Detailed Results</h2>']

    for r in ctx.results:
        badge = f'<span class="badge badge-{_esc(r.coverage_status)}">{_esc(r.coverage_status)}</span>'
        region = _fmt_region(r.contig, r.start, r.end)

        metrics = [
            ("Mean Depth", f"{r.mean_depth:.2f}x"),
            ("Median Depth", f"{r.median_depth:.0f}x"),
            ("Min / Max", f"{r.min_depth} / {r.max_depth}"),
            ("Std Deviation", f"{r.stdev_depth:.2f}"),
            ("% at 0x", f"{r.pct_zero:.2f}%"),
        ]
        for t in sorted_thresholds:
            metrics.append((f"% >= {t}x", f"{r.pct_thresholds.get(t, 0.0):.2f}%"))
        metrics.append(("Low-cov Blocks", f"{r.n_lowcov_blocks} ({r.lowcov_total_bp:,} bp)"))

        grid = "".join(
            f'<div class="metric-item"><div class="metric-label">{_esc(label)}</div><div class="metric-value">{_esc(val)}</div></div>'
            for label, val in metrics
        )

        desc = _status_description(r.coverage_status)
        rationale = _rationale(r, ctx.classify_params)

        sections.append(f"""
<div class="detail-panel">
  <h3>{_esc(r.target_id)} — {badge}</h3>
  <p style="color:var(--secondary);font-size:0.85rem;margin:0.5rem 0">{_esc(region)} &middot; {r.length_bp:,} bp &middot; {_esc(desc)}</p>
  <div class="metric-grid">{grid}</div>
  <div class="rationale"><strong>Classification rationale:</strong> {_esc(rationale)}</div>
</div>""")

    return "\n".join(sections)


def _build_exon_sections(ctx: ReportContext) -> str:
    if not ctx.exon_results:
        return ""

    sections = []
    for gene_name, exon_list in ctx.exon_results.items():
        metadata_list = (ctx.exon_metadata or {}).get(gene_name, [])
        status_list = (ctx.exon_statuses or {}).get(gene_name, [])

        rows = ""
        for i, er in enumerate(exon_list):
            meta = metadata_list[i] if i < len(metadata_list) else {}
            exon_id = meta.get("exon_id", er.target_id)
            exon_num = meta.get("exon_number", "?")
            status = status_list[i] if i < len(status_list) else "OK"
            pct_20 = er.pct_thresholds.get(20, 0.0)
            badge = f'<span class="badge badge-{_esc(status)}">{_esc(status)}</span>'

            rows += f"""<tr>
  <td>{_esc(exon_num)}</td>
  <td>{_esc(exon_id)}</td>
  <td>{_esc(_fmt_region(er.contig, er.start, er.end))}</td>
  <td>{er.length_bp:,}</td>
  <td>{er.mean_depth:.2f}</td>
  <td>{er.median_depth:.0f}</td>
  <td>{er.pct_zero:.2f}%</td>
  <td>{pct_20:.2f}%</td>
  <td>{badge}</td>
</tr>"""

        sections.append(f"""
<h2>Exon-Level Results ({_esc(gene_name)})</h2>
<div class="table-wrap">
<table>
  <thead><tr>
    <th>Exon</th><th>Exon ID</th><th>Region</th><th>Length</th>
    <th>Mean</th><th>Median</th><th>%Zero</th><th>%>=20x</th><th>Flag</th>
  </tr></thead>
  <tbody>{rows}</tbody>
</table>
</div>
<div id="heatmap-{_esc(gene_name)}" class="chart-container"><h3>Exon Coverage Heatmap</h3><div class="heatmap-target" data-gene="{_esc(gene_name)}"></div></div>
""")

    return "\n".join(sections)


def _build_lowcov_section(ctx: ReportContext) -> str:
    any_blocks = any(r.lowcov_blocks for r in ctx.results)
    header = f"""
<h2>Low-Coverage Blocks</h2>
<p style="color:var(--secondary);font-size:0.85rem">Threshold: depth &lt; {ctx.lowcov_threshold}, min length: {ctx.lowcov_min_len} bp</p>"""

    if not any_blocks:
        return header + '<p style="margin-top:1rem">No low-coverage blocks detected.</p>'

    rows = ""
    for r in ctx.results:
        for b in r.lowcov_blocks:
            rows += f"""<tr>
  <td>{_esc(r.target_id)}</td>
  <td>{_esc(_fmt_region(r.contig, b.start, b.end))}</td>
  <td>{b.length:,}</td>
  <td>{b.mean_depth:.1f}x</td>
</tr>"""

    return header + f"""
<div class="table-wrap">
<table>
  <thead><tr><th>Target</th><th>Block</th><th>Length</th><th>Mean Depth</th></tr></thead>
  <tbody>{rows}</tbody>
</table>
</div>"""


def _build_heuristics_reference(ctx: ReportContext) -> str:
    p = ctx.classify_params
    return f"""
<h2>Classification Heuristics Reference</h2>
<div class="table-wrap">
<table class="heuristics">
  <thead><tr><th>Status</th><th>Rule</th><th>Evaluation Order</th></tr></thead>
  <tbody>
    <tr><td><span class="badge badge-DROP_OUT">DROP_OUT</span></td><td>pct_zero &gt; {p.dropout_pct_zero}% OR zero-block &ge; {p.dropout_zero_block_bp} bp</td><td>1st</td></tr>
    <tr><td><span class="badge badge-UNEVEN">UNEVEN</span></td><td>mean_depth &gt; 20 AND CV &gt; {p.uneven_cv}</td><td>2nd</td></tr>
    <tr><td><span class="badge badge-LOW_EXON">LOW_EXON</span></td><td>(exon mode) any exon pct_ge_20 &lt; {p.exon_pct_ge_20}% or pct_zero &gt; {p.exon_max_pct_zero}%</td><td>3rd</td></tr>
    <tr><td><span class="badge badge-LOW_COVERAGE">LOW_COVERAGE</span></td><td>pct_ge_20 &lt; {p.pass_pct_ge_20}%</td><td>4th</td></tr>
    <tr><td><span class="badge badge-PASS">PASS</span></td><td>pct_ge_20 &ge; {p.pass_pct_ge_20}% AND pct_zero &le; {p.pass_max_pct_zero}%</td><td>5th</td></tr>
  </tbody>
</table>
</div>
<div class="footer">Generated by covsnap {__version__}</div>"""


def _build_charts_container() -> str:
    return """
<h2>Charts</h2>
<div class="chart-section">
  <div class="chart-container"><h3>Mean Depth by Target</h3><div id="chart-depth"></div></div>
  <div class="chart-container"><h3>Threshold Coverage Distribution</h3><div id="chart-thresholds"></div></div>
</div>
<div class="chart-tooltip" id="tooltip"></div>"""


# ---------------------------------------------------------------------------
# JSON data for charts
# ---------------------------------------------------------------------------

def _build_targets_json(ctx: ReportContext) -> str:
    sorted_thresholds = sorted(ctx.thresholds)
    data = []
    for r in ctx.results:
        entry = {
            "id": r.target_id,
            "mean": r.mean_depth,
            "median": r.median_depth,
            "min": r.min_depth,
            "max": r.max_depth,
            "stdev": r.stdev_depth,
            "pct_zero": r.pct_zero,
            "status": r.coverage_status,
            "thresholds": {str(t): r.pct_thresholds.get(t, 0.0) for t in sorted_thresholds},
        }
        data.append(entry)
    return json.dumps(data)


def _build_exon_json(ctx: ReportContext) -> str:
    if not ctx.exon_results:
        return "{}"
    data: dict[str, list[dict]] = {}
    for gene_name, exon_list in ctx.exon_results.items():
        metadata_list = (ctx.exon_metadata or {}).get(gene_name, [])
        status_list = (ctx.exon_statuses or {}).get(gene_name, [])
        exons = []
        for i, er in enumerate(exon_list):
            meta = metadata_list[i] if i < len(metadata_list) else {}
            status = status_list[i] if i < len(status_list) else "OK"
            exons.append({
                "num": meta.get("exon_number", i + 1),
                "id": meta.get("exon_id", er.target_id),
                "mean": er.mean_depth,
                "pct_ge_20": er.pct_thresholds.get(20, 0.0),
                "pct_zero": er.pct_zero,
                "status": status,
                "length": er.length_bp,
            })
        data[gene_name] = exons
    return json.dumps(data)


# ---------------------------------------------------------------------------
# JavaScript (inline)
# ---------------------------------------------------------------------------

def _build_js() -> str:
    return """<script>
(function() {
  'use strict';

  // --- Table sorting ---
  const table = document.getElementById('targetTable');
  if (table) {
    const headers = table.querySelectorAll('th');
    let sortCol = -1, sortAsc = true;

    headers.forEach((th, idx) => {
      th.addEventListener('click', () => {
        if (sortCol === idx) { sortAsc = !sortAsc; }
        else { sortCol = idx; sortAsc = true; }
        sortTable(idx, sortAsc, th.dataset.sort || 'text');
        headers.forEach(h => h.querySelector('.sort-arrow').classList.remove('active'));
        th.querySelector('.sort-arrow').classList.add('active');
        th.querySelector('.sort-arrow').innerHTML = sortAsc ? '&#9650;' : '&#9660;';
      });
    });

    function sortTable(col, asc, type) {
      const tbody = table.querySelector('tbody');
      const rows = Array.from(tbody.querySelectorAll('tr'));
      rows.sort((a, b) => {
        let va = a.cells[col].textContent.trim().replace(/[,%]/g, '');
        let vb = b.cells[col].textContent.trim().replace(/[,%]/g, '');
        if (type === 'num') {
          va = parseFloat(va) || 0;
          vb = parseFloat(vb) || 0;
          return asc ? va - vb : vb - va;
        }
        return asc ? va.localeCompare(vb) : vb.localeCompare(va);
      });
      rows.forEach(r => tbody.appendChild(r));
    }
  }

  // --- Table filtering ---
  window.filterTable = function() {
    const text = document.getElementById('tableFilter').value.toLowerCase();
    const status = document.getElementById('statusFilter').value;
    const rows = document.querySelectorAll('#targetTable tbody tr');
    rows.forEach(r => {
      const matchText = !text || r.textContent.toLowerCase().includes(text);
      const matchStatus = !status || r.textContent.includes(status);
      r.style.display = (matchText && matchStatus) ? '' : 'none';
    });
  };

  // --- Tooltip helper ---
  const tooltip = document.getElementById('tooltip');
  function showTip(evt, html) {
    tooltip.innerHTML = html;
    tooltip.style.opacity = '1';
    tooltip.style.left = (evt.pageX + 10) + 'px';
    tooltip.style.top = (evt.pageY - 30) + 'px';
  }
  function hideTip() { tooltip.style.opacity = '0'; }

  // --- SVG helper ---
  function svg(tag, attrs, children) {
    const el = document.createElementNS('http://www.w3.org/2000/svg', tag);
    for (const [k, v] of Object.entries(attrs || {})) el.setAttribute(k, v);
    if (typeof children === 'string') el.textContent = children;
    else if (children) children.forEach(c => { if (c) el.appendChild(c); });
    return el;
  }

  // --- Colors ---
  const statusColor = {PASS:'#22C55E', LOW_COVERAGE:'#F59E0B', DROP_OUT:'#EF4444', UNEVEN:'#818CF8', LOW_EXON:'#FB7185'};
  const thresholdColors = ['#DBEAFE','#BFDBFE','#93C5FD','#60A5FA','#3B82F6','#2563EB','#1D4ED8'];

  // --- Chart 1: Mean Depth Bar ---
  function drawDepthChart() {
    const container = document.getElementById('chart-depth');
    if (!container || !TARGETS.length) return;
    const W = Math.max(container.clientWidth - 20, 400);
    const barW = Math.min(60, (W - 80) / TARGETS.length - 4);
    const H = 300;
    const maxVal = Math.max(...TARGETS.map(t => t.mean), 10);
    const padL = 55, padB = 80, padT = 20;
    const chartH = H - padB - padT;
    const chartW = TARGETS.length * (barW + 4);

    const s = svg('svg', {width: Math.max(W, chartW + padL + 20), height: H, style:'display:block;margin:0 auto'});

    // Y axis
    for (let i = 0; i <= 4; i++) {
      const val = (maxVal / 4 * i).toFixed(0);
      const y = padT + chartH - (chartH * i / 4);
      s.appendChild(svg('line', {x1:padL, y1:y, x2:padL+chartW, y2:y, stroke:'#E2E8F0', 'stroke-width':1}));
      s.appendChild(svg('text', {x:padL-8, y:y+4, 'text-anchor':'end', fill:'#64748B', 'font-size':'11'}, val + 'x'));
    }

    TARGETS.forEach((t, i) => {
      const x = padL + i * (barW + 4) + 2;
      const h = (t.mean / maxVal) * chartH;
      const y = padT + chartH - h;
      const bar = svg('rect', {x:x, y:y, width:barW, height:h, rx:2, fill: statusColor[t.status] || '#3B82F6', class:'chart-bar'});
      bar.addEventListener('mousemove', e => showTip(e, '<b>'+t.id+'</b><br>Mean: '+t.mean.toFixed(2)+'x<br>Status: '+t.status));
      bar.addEventListener('mouseleave', hideTip);
      s.appendChild(bar);

      const label = svg('text', {x: x + barW/2, y: H - padB + 14, 'text-anchor':'end', 'font-size':'10', fill:'#64748B',
        transform: 'rotate(-45,' + (x + barW/2) + ',' + (H - padB + 14) + ')'}, t.id);
      s.appendChild(label);
    });

    container.appendChild(s);
  }

  // --- Chart 2: Threshold Stacked Bar ---
  function drawThresholdChart() {
    const container = document.getElementById('chart-thresholds');
    if (!container || !TARGETS.length) return;
    const W = Math.max(container.clientWidth - 20, 400);
    const barW = Math.min(60, (W - 80) / TARGETS.length - 4);
    const H = 300;
    const padL = 55, padB = 80, padT = 20;
    const chartH = H - padB - padT;
    const chartW = TARGETS.length * (barW + 4);

    const keys = Object.keys(TARGETS[0].thresholds).sort((a,b) => +a - +b);
    const s = svg('svg', {width: Math.max(W, chartW + padL + 20), height: H, style:'display:block;margin:0 auto'});

    // Y axis (0-100%)
    for (let i = 0; i <= 4; i++) {
      const y = padT + chartH - (chartH * i / 4);
      s.appendChild(svg('line', {x1:padL, y1:y, x2:padL+chartW, y2:y, stroke:'#E2E8F0', 'stroke-width':1}));
      s.appendChild(svg('text', {x:padL-8, y:y+4, 'text-anchor':'end', fill:'#64748B', 'font-size':'11'}, (25*i)+'%'));
    }

    TARGETS.forEach((t, i) => {
      const x = padL + i * (barW + 4) + 2;
      // Draw stacked segments from highest threshold down
      const reversedKeys = [...keys].reverse();
      reversedKeys.forEach((k, ki) => {
        const pct = t.thresholds[k] || 0;
        const h = (pct / 100) * chartH;
        const y = padT + chartH - h;
        const color = thresholdColors[Math.min(keys.indexOf(k), thresholdColors.length - 1)];
        const bar = svg('rect', {x:x, y:y, width:barW, height:h, fill:color, class:'chart-bar', opacity: 0.9});
        bar.addEventListener('mousemove', e => showTip(e, '<b>'+t.id+'</b><br>>=' + k + 'x: ' + pct.toFixed(1) + '%'));
        bar.addEventListener('mouseleave', hideTip);
        s.appendChild(bar);
      });

      const label = svg('text', {x: x + barW/2, y: H - padB + 14, 'text-anchor':'end', 'font-size':'10', fill:'#64748B',
        transform: 'rotate(-45,' + (x + barW/2) + ',' + (H - padB + 14) + ')'}, t.id);
      s.appendChild(label);
    });

    // Legend
    const legendY = padT;
    keys.forEach((k, i) => {
      const lx = padL + chartW + 10;
      const ly = legendY + i * 18;
      s.appendChild(svg('rect', {x:lx, y:ly, width:12, height:12, rx:2, fill:thresholdColors[Math.min(i, thresholdColors.length-1)]}));
      s.appendChild(svg('text', {x:lx+16, y:ly+10, 'font-size':'10', fill:'#64748B'}, '>=' + k + 'x'));
    });

    container.appendChild(s);
  }

  // --- Chart 3: Exon Heatmap ---
  function drawExonHeatmaps() {
    if (!EXONS || !Object.keys(EXONS).length) return;
    for (const [gene, exons] of Object.entries(EXONS)) {
      const target = document.querySelector('.heatmap-target[data-gene="'+gene+'"]');
      if (!target || !exons.length) continue;

      const cellW = Math.min(50, Math.max(20, (target.clientWidth - 80) / exons.length));
      const cellH = 40;
      const padL = 10, padT = 30, padB = 60;
      const W = padL + exons.length * cellW + 20;
      const H = padT + cellH + padB;

      const s = svg('svg', {width: W, height: H, style:'display:block'});

      exons.forEach((ex, i) => {
        const x = padL + i * cellW;
        // Color based on pct_ge_20: green (100) -> yellow (90) -> red (0)
        const pct = ex.pct_ge_20;
        let fill;
        if (pct >= 95) fill = '#22C55E';
        else if (pct >= 90) fill = '#86EFAC';
        else if (pct >= 80) fill = '#FDE047';
        else if (pct >= 60) fill = '#FB923C';
        else fill = '#EF4444';

        const rect = svg('rect', {x:x, y:padT, width:cellW-1, height:cellH, rx:3, fill:fill, class:'heatmap-cell'});
        rect.addEventListener('mousemove', e => showTip(e,
          '<b>Exon '+ex.num+'</b> ('+ex.id+')<br>Mean: '+ex.mean.toFixed(1)+'x<br>%>=20x: '+ex.pct_ge_20.toFixed(1)+'%<br>Status: '+ex.status));
        rect.addEventListener('mouseleave', hideTip);
        s.appendChild(rect);

        // Exon number label
        s.appendChild(svg('text', {x: x + (cellW-1)/2, y: padT + cellH/2 + 4, 'text-anchor':'middle',
          'font-size':'10', fill:'white', 'font-weight':'600'}, String(ex.num)));
        // Bottom label
        s.appendChild(svg('text', {x: x + (cellW-1)/2, y: padT + cellH + 16, 'text-anchor':'middle',
          'font-size':'9', fill:'#64748B'}, ex.pct_ge_20.toFixed(0) + '%'));
      });

      // Color legend
      const legendItems = [
        {label:'>=95%', color:'#22C55E'}, {label:'90-95%', color:'#86EFAC'},
        {label:'80-90%', color:'#FDE047'}, {label:'60-80%', color:'#FB923C'}, {label:'<60%', color:'#EF4444'}
      ];
      legendItems.forEach((li, idx) => {
        const lx = padL + idx * 70;
        const ly = padT + cellH + 32;
        s.appendChild(svg('rect', {x:lx, y:ly, width:10, height:10, rx:2, fill:li.color}));
        s.appendChild(svg('text', {x:lx+14, y:ly+9, 'font-size':'9', fill:'#64748B'}, li.label));
      });

      target.appendChild(s);
    }
  }

  // Init
  drawDepthChart();
  drawThresholdChart();
  drawExonHeatmaps();
})();
</script>"""
```

- [ ] **Step 4: Run the tests**

Run: `cd /home/enes/Desktop/denemeler/covinspector && python -m pytest tests/test_html_report.py -v`

Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add src/covsnap/html_report.py tests/test_html_report.py
git commit -m "feat: add self-contained interactive HTML report module

Generates a single HTML file with all CSS/JS inline.
Features: summary cards, sortable/filterable tables,
mean depth bar chart, threshold stacked bars, exon heatmap.
Light theme with professional color palette."
```

---

## Task 3: Wire HTML report into CLI

**Files:**
- Modify: `src/covsnap/cli.py:32,127-146,554-635`
- Test: `tests/test_cli_smoke.py`

- [ ] **Step 1: Write a CLI smoke test for --html-out**

Add to `tests/test_cli_smoke.py`:

```python
class TestCLIHtmlOutput:
    @pytest.fixture(autouse=True)
    def _samtools(self):
        _requires_samtools()

    def test_html_output_produced(self, synthetic_bam, tmp_output_dir):
        raw_out = str(tmp_output_dir / "raw.tsv")
        report_out = str(tmp_output_dir / "report.md")
        html_out = str(tmp_output_dir / "report.html")

        main([
            synthetic_bam,
            "BRCA1",
            "--raw-out", raw_out,
            "--report-out", report_out,
            "--html-out", html_out,
            "--engine", "samtools",
        ])

        assert os.path.isfile(html_out)
        with open(html_out) as f:
            html = f.read()
        assert "<!DOCTYPE html>" in html
        assert "BRCA1" in html
        assert "SYNTH_SAMPLE" in html

    def test_html_output_with_bed(self, synthetic_bam, sample_bed, tmp_output_dir):
        raw_out = str(tmp_output_dir / "raw.tsv")
        report_out = str(tmp_output_dir / "report.md")
        html_out = str(tmp_output_dir / "report.html")

        main([
            synthetic_bam,
            "--bed", sample_bed,
            "--raw-out", raw_out,
            "--report-out", report_out,
            "--html-out", html_out,
            "--engine", "samtools",
        ])

        assert os.path.isfile(html_out)
        with open(html_out) as f:
            html = f.read()
        assert "region_1" in html
        assert "region_2" in html
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd /home/enes/Desktop/denemeler/covinspector && python -m pytest tests/test_cli_smoke.py::TestCLIHtmlOutput -v`

Expected: FAIL — unrecognized argument --html-out

- [ ] **Step 3: Add --html-out argument and wire the writer**

In `src/covsnap/cli.py`:

1. Add import at line 33 (after the report import):
```python
from covsnap.html_report import write_html_report
```

2. Add `--html-out` argument in `build_parser()`, inside the "Output options" group (after `--json-out`):
```python
    out.add_argument(
        "--html-out",
        metavar="FILE",
        default=None,
        help="Optional interactive HTML report (self-contained, no external dependencies).",
    )
```

3. Add HTML output writing after the JSON section (after line 603 in the current file, inside `main()`):
```python
    # HTML report (optional)
    if args.html_out:
        write_html_report(args.html_out, report_ctx)
        logger.info("HTML report written to %s", args.html_out)
```

4. Update the summary line to mention HTML when used (modify the stderr summary block):
```python
    # ── Summary to stderr ──
    if not args.quiet:
        n_pass = sum(1 for r in results if r.coverage_status == "PASS")
        n_total = len(results)
        out_parts = [f"Raw: {args.raw_out}", f"Report: {args.report_out}"]
        if args.html_out:
            out_parts.append(f"HTML: {args.html_out}")
        print(
            f"[covsnap] Done. {n_pass}/{n_total} targets PASS. "
            + "  ".join(out_parts),
            file=sys.stderr,
        )
```

- [ ] **Step 4: Run the HTML CLI tests**

Run: `cd /home/enes/Desktop/denemeler/covinspector && python -m pytest tests/test_cli_smoke.py::TestCLIHtmlOutput -v`

Expected: PASS

- [ ] **Step 5: Run full test suite**

Run: `cd /home/enes/Desktop/denemeler/covinspector && python -m pytest tests/ -v`

Expected: All tests PASS

- [ ] **Step 6: Commit**

```bash
git add src/covsnap/cli.py tests/test_cli_smoke.py
git commit -m "feat: wire --html-out flag into CLI

Adds --html-out argument to generate self-contained HTML reports.
HTML output is optional and works alongside existing TSV/MD/JSON outputs."
```
