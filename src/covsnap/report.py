"""Coverage classification and human-readable Markdown report generation.

Classification heuristics (evaluated in order, first match wins):

    DROP_OUT     — pct_zero > dropout_pct_zero OR zero-block >= 500 bp
    UNEVEN       — mean_depth > 20 AND cv > uneven_cv
    LOW_EXON     — (exon mode) any exon pct_ge_20 < exon_pct_ge_20 or pct_zero > exon_max_pct_zero
    LOW_COVERAGE — pct_ge_20 < pass_pct_ge_20
    PASS         — pct_ge_20 >= pass_pct_ge_20 AND pct_zero <= pass_max_pct_zero
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Optional

from covsnap import ANNOTATION_VERSION, BUILD, __version__
from covsnap.metrics import TargetResult


# ---------------------------------------------------------------------------
# Classification parameters
# ---------------------------------------------------------------------------


@dataclass
class ClassifyParams:
    """Thresholds for coverage classification."""

    pass_pct_ge_20: float = 95.0
    pass_max_pct_zero: float = 1.0
    dropout_pct_zero: float = 5.0
    dropout_zero_block_bp: int = 500
    uneven_cv: float = 1.0
    exon_pct_ge_20: float = 90.0
    exon_max_pct_zero: float = 5.0


# ---------------------------------------------------------------------------
# Classifier
# ---------------------------------------------------------------------------


def classify_target(
    result: TargetResult,
    params: ClassifyParams,
    exon_results: Optional[list[TargetResult]] = None,
) -> str:
    """Assign a coverage_status to a TargetResult.

    Mutates ``result.coverage_status`` in-place and returns the status string.
    """
    pct_zero = result.pct_zero
    pct_ge_20 = result.pct_thresholds.get(20, 0.0)
    mean = result.mean_depth
    stdev = result.stdev_depth

    # 1. DROP_OUT
    if pct_zero > params.dropout_pct_zero:
        result.coverage_status = "DROP_OUT"
        return "DROP_OUT"

    # Check for long zero-depth blocks
    for block in result.lowcov_blocks:
        if block.length >= params.dropout_zero_block_bp:
            # Check if the block is actually zero-depth (not just low)
            if block.mean_depth < 0.5:
                result.coverage_status = "DROP_OUT"
                return "DROP_OUT"

    # 2. UNEVEN
    if mean > 20:
        cv = stdev / mean if mean > 0 else 0.0
        if cv > params.uneven_cv:
            result.coverage_status = "UNEVEN"
            return "UNEVEN"

    # 3. LOW_EXON (only if exon results provided)
    if exon_results:
        for er in exon_results:
            exon_pct_ge_20 = er.pct_thresholds.get(20, 0.0)
            if exon_pct_ge_20 < params.exon_pct_ge_20 or er.pct_zero > params.exon_max_pct_zero:
                result.coverage_status = "LOW_EXON"
                return "LOW_EXON"

    # 4. LOW_COVERAGE
    if pct_ge_20 < params.pass_pct_ge_20:
        result.coverage_status = "LOW_COVERAGE"
        return "LOW_COVERAGE"

    # 5. PASS
    if pct_ge_20 >= params.pass_pct_ge_20 and pct_zero <= params.pass_max_pct_zero:
        result.coverage_status = "PASS"
        return "PASS"

    # Edge case: pct_ge_20 is sufficient but pct_zero is between
    # pass_max_pct_zero and dropout_pct_zero
    result.coverage_status = "LOW_COVERAGE"
    return "LOW_COVERAGE"


def classify_exon(
    result: TargetResult,
    params: ClassifyParams,
) -> str:
    """Classify a single exon. Returns 'OK' or 'LOW_EXON'."""
    pct_ge_20 = result.pct_thresholds.get(20, 0.0)
    if pct_ge_20 < params.exon_pct_ge_20 or result.pct_zero > params.exon_max_pct_zero:
        return "LOW_EXON"
    return "OK"


# ---------------------------------------------------------------------------
# Report context
# ---------------------------------------------------------------------------


@dataclass
class ReportContext:
    """All data needed to render the report."""

    results: list[TargetResult]
    engine_used: str
    engine_version: str
    bam_path: str
    sample_name: str
    contig_style: str
    reference: Optional[str] = None
    bed_path: Optional[str] = None
    bed_guardrail_message: str = ""
    exon_results: Optional[dict[str, list[TargetResult]]] = None
    exon_metadata: Optional[dict[str, list[dict[str, Any]]]] = None
    exon_statuses: Optional[dict[str, list[str]]] = None
    classify_params: ClassifyParams = field(default_factory=ClassifyParams)
    lowcov_threshold: int = 10
    lowcov_min_len: int = 50
    emit_lowcov: bool = False
    lowcov_bed_path: str = ""
    thresholds: list[int] = field(default_factory=lambda: [1, 5, 10, 20, 30, 50, 100])
    run_date: Optional[str] = None


# ---------------------------------------------------------------------------
# Markdown report writer
# ---------------------------------------------------------------------------


def write_report(path: str, ctx: ReportContext) -> None:
    """Write the interpreted coverage report as Markdown."""
    if ctx.run_date is None:
        ctx.run_date = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    lines: list[str] = []
    _w = lines.append

    # ── Header ──
    _w("# covsnap Coverage Report")
    _w("")
    _w("| Field | Value |")
    _w("|-------|-------|")
    _w(f"| **Tool version** | {__version__} |")
    _w(f"| **Date** | {ctx.run_date} |")
    _w(f"| **Engine** | {ctx.engine_used} {ctx.engine_version} |")
    _w(f"| **Annotation** | GENCODE v44 ({BUILD}) |")
    _w(f"| **Input** | `{ctx.bam_path}` |")
    _w(f"| **Sample name** | {ctx.sample_name} |")
    _w(f"| **Contig style** | {ctx.contig_style}-prefixed (auto-detected) |")
    if ctx.reference:
        _w(f"| **Reference** | `{ctx.reference}` |")
    if ctx.bed_path:
        _w(f"| **Input BED** | `{ctx.bed_path}` |")
    _w("")

    # ── Coordinate convention ──
    _w("## Coordinate Convention")
    _w("")
    _w("- Input region strings are interpreted as **1-based inclusive** (samtools convention).")
    _w("- Input BED files use **0-based half-open** (BED standard).")
    _w("- All output BED coordinates and raw TSV start/end columns use **0-based half-open**.")
    _w("- Report display positions are **1-based inclusive** for readability.")
    _w("")
    _w("---")
    _w("")

    # ── BED Guardrails ──
    _w("## BED Guardrails")
    _w("")
    if ctx.bed_guardrail_message:
        _w(f"> **WARNING:** {ctx.bed_guardrail_message}")
    else:
        _w("No BED guardrails were triggered.")
    _w("")
    _w("---")
    _w("")

    # ── Target Summary Table ──
    _w("## Target Summary")
    _w("")
    _w("| Target | Region (1-based) | Length | Mean Depth | Median | %\u226520x | %Zero | Status |")
    _w("|--------|-------------------|--------|------------|--------|-------|-------|--------|")

    for r in ctx.results:
        region_str = _format_region_1based(r.contig, r.start, r.end)
        length_str = f"{r.length_bp:,} bp"
        pct_20 = r.pct_thresholds.get(20, 0.0)
        status_fmt = f"**{r.coverage_status}**"
        _w(
            f"| {r.target_id} "
            f"| {region_str} "
            f"| {length_str} "
            f"| {r.mean_depth:.2f} "
            f"| {r.median_depth:.0f} "
            f"| {pct_20:.2f}% "
            f"| {r.pct_zero:.2f}% "
            f"| {status_fmt} |"
        )

    _w("")
    _w("---")
    _w("")

    # ── Detailed per-target results ──
    _w("## Detailed Results")
    _w("")

    sorted_thresholds = sorted(ctx.thresholds)

    for r in ctx.results:
        _w(f"### {r.target_id} \u2014 {r.coverage_status}")
        _w("")
        _w(_status_description(r.coverage_status))
        _w("")
        _w("| Metric | Value |")
        _w("|--------|-------|")
        _w(f"| Mean depth | {r.mean_depth:.2f}x |")
        _w(f"| Median depth | {r.median_depth:.0f}x |")
        _w(f"| Min / Max | {r.min_depth} / {r.max_depth} |")
        _w(f"| Std deviation | {r.stdev_depth:.2f} |")
        _w(f"| % at 0x | {r.pct_zero:.2f}% |")
        for t in sorted_thresholds:
            _w(f"| % \u2265 {t}x | {r.pct_thresholds.get(t, 0.0):.2f}% |")
        _w(f"| Low-coverage blocks | {r.n_lowcov_blocks} blocks, {r.lowcov_total_bp:,} bp total |")
        _w("")

        # Classification rationale
        pct_20 = r.pct_thresholds.get(20, 0.0)
        _w(f"**Classification rationale:** {_rationale(r, ctx.classify_params)}")
        _w("")

    # ── Exon-level results ──
    if ctx.exon_results:
        for gene_name, exon_list in ctx.exon_results.items():
            _w(f"## Exon-Level Results ({gene_name})")
            _w("")
            _w("| Exon | Exon ID | Region (1-based) | Length | Mean | %\u226520x | %Zero | Flag |")
            _w("|------|---------|-------------------|--------|------|-------|-------|------|")

            metadata_list = (ctx.exon_metadata or {}).get(gene_name, [])
            status_list = (ctx.exon_statuses or {}).get(gene_name, [])

            for i, er in enumerate(exon_list):
                meta = metadata_list[i] if i < len(metadata_list) else {}
                exon_id = meta.get("exon_id", er.target_id)
                exon_num = meta.get("exon_number", "?")
                status = status_list[i] if i < len(status_list) else "OK"
                region_str = _format_region_1based(er.contig, er.start, er.end)
                pct_20 = er.pct_thresholds.get(20, 0.0)

                bold = "**" if status == "LOW_EXON" else ""
                _w(
                    f"| {bold}{exon_num}{bold} "
                    f"| {bold}{exon_id}{bold} "
                    f"| {bold}{region_str}{bold} "
                    f"| {bold}{er.length_bp:,} bp{bold} "
                    f"| {bold}{er.mean_depth:.2f}{bold} "
                    f"| {bold}{pct_20:.2f}%{bold} "
                    f"| {bold}{er.pct_zero:.2f}%{bold} "
                    f"| {bold}{status}{bold} |"
                )
            _w("")

    # ── Low-coverage blocks ──
    if ctx.emit_lowcov:
        _w("## Low-Coverage Blocks")
        _w("")
        _w(f"*(Threshold: depth < {ctx.lowcov_threshold}, min length: {ctx.lowcov_min_len} bp)*")
        _w("")

        any_blocks = any(r.lowcov_blocks for r in ctx.results)
        if any_blocks:
            _w("| Target | Block (1-based) | Length | Mean Depth |")
            _w("|--------|-----------------|--------|------------|")
            for r in ctx.results:
                for b in r.lowcov_blocks:
                    block_str = _format_region_1based(r.contig, b.start, b.end)
                    _w(
                        f"| {r.target_id} | {block_str} "
                        f"| {b.length:,} bp | {b.mean_depth:.1f}x |"
                    )
            _w("")
            if ctx.lowcov_bed_path:
                _w(f"Full low-coverage BED written to: `{ctx.lowcov_bed_path}`")
        else:
            _w("No low-coverage blocks detected.")
        _w("")

    _w("---")
    _w("")

    # ── Classification summary ──
    _w("## Classification Summary")
    _w("")
    status_counts: dict[str, int] = {}
    for r in ctx.results:
        status_counts[r.coverage_status] = status_counts.get(r.coverage_status, 0) + 1
    total = len(ctx.results)

    _w("| Status | Count | Percentage |")
    _w("|--------|-------|------------|")
    for status in ["PASS", "LOW_COVERAGE", "LOW_EXON", "DROP_OUT", "UNEVEN"]:
        count = status_counts.get(status, 0)
        if count > 0 or status in ("PASS", "LOW_COVERAGE"):
            pct = count / total * 100 if total > 0 else 0
            _w(f"| {status} | {count} | {pct:.1f}% |")
    _w("")

    _w("---")
    _w("")

    # ── Heuristics reference ──
    _w("## Classification Heuristics Reference")
    _w("")
    _w("| Status | Rule |")
    _w("|--------|------|")
    p = ctx.classify_params
    _w(f"| PASS | pct_ge_20 \u2265 {p.pass_pct_ge_20}% AND pct_zero \u2264 {p.pass_max_pct_zero}% |")
    _w(f"| LOW_COVERAGE | pct_ge_20 < {p.pass_pct_ge_20}% (and not DROP_OUT) |")
    _w(f"| DROP_OUT | pct_zero > {p.dropout_pct_zero}% OR zero-block \u2265 {p.dropout_zero_block_bp} bp |")
    _w(f"| UNEVEN | mean_depth > 20 AND coefficient_of_variation > {p.uneven_cv} |")
    _w(f"| LOW_EXON | (exon mode) any exon pct_ge_20 < {p.exon_pct_ge_20}% or pct_zero > {p.exon_max_pct_zero}% |")
    _w("")
    _w("Evaluation order: DROP_OUT \u2192 UNEVEN \u2192 LOW_EXON \u2192 LOW_COVERAGE \u2192 PASS.")
    _w("")
    _w("---")
    _w("")
    _w(f"*Generated by covsnap {__version__}*")
    _w("")

    with open(path, "w") as fh:
        fh.write("\n".join(lines))


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _format_region_1based(contig: str, start_0: int, end_0: int) -> str:
    """Format a 0-based half-open region as 1-based inclusive for display."""
    return f"{contig}:{start_0 + 1:,}-{end_0:,}"


def _status_description(status: str) -> str:
    """Human-readable one-liner for each status."""
    descs = {
        "PASS": "Coverage is adequate across the target region.",
        "LOW_COVERAGE": "Coverage is below the 20x threshold across a significant fraction of the target.",
        "DROP_OUT": "Complete loss of coverage detected in substantial portions of the target.",
        "UNEVEN": "High average coverage but highly variable distribution across the target.",
        "LOW_EXON": "One or more exons have critically low coverage.",
    }
    return descs.get(status, "")


def _rationale(r: TargetResult, params: ClassifyParams) -> str:
    """Generate a human-readable rationale string for the classification."""
    pct_20 = r.pct_thresholds.get(20, 0.0)
    status = r.coverage_status

    if status == "PASS":
        return (
            f"pct_ge_20 ({pct_20:.2f}%) \u2265 {params.pass_pct_ge_20}% "
            f"AND pct_zero ({r.pct_zero:.2f}%) \u2264 {params.pass_max_pct_zero}% \u2192 PASS"
        )
    if status == "LOW_COVERAGE":
        return (
            f"pct_ge_20 ({pct_20:.2f}%) < {params.pass_pct_ge_20}% \u2192 LOW_COVERAGE"
        )
    if status == "DROP_OUT":
        if r.pct_zero > params.dropout_pct_zero:
            return (
                f"pct_zero ({r.pct_zero:.2f}%) > {params.dropout_pct_zero}% \u2192 DROP_OUT"
            )
        return f"Long zero-depth block detected (\u2265 {params.dropout_zero_block_bp} bp) \u2192 DROP_OUT"
    if status == "UNEVEN":
        cv = r.stdev_depth / r.mean_depth if r.mean_depth > 0 else 0
        return (
            f"mean_depth ({r.mean_depth:.2f}) > 20 AND "
            f"CV ({cv:.2f}) > {params.uneven_cv} \u2192 UNEVEN"
        )
    if status == "LOW_EXON":
        return "One or more exons below threshold \u2192 LOW_EXON"
    return status
