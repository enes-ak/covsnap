"""Coverage classification and report context.

Classification heuristics (evaluated in order, first match wins):

    DROP_OUT     — pct_zero > dropout_pct_zero OR zero-block >= 500 bp
    UNEVEN       — mean_depth > 20 AND cv > uneven_cv
    LOW_EXON     — (exon mode) any exon pct_ge_20 < exon_pct_ge_20 or pct_zero > exon_max_pct_zero
    LOW_COVERAGE — pct_ge_20 < pass_pct_ge_20
    PASS         — pct_ge_20 >= pass_pct_ge_20 AND pct_zero <= pass_max_pct_zero
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional

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
    # Gene-level results for region mode (region → overlapping genes)
    gene_results: Optional[list[TargetResult]] = None
    gene_metadata: Optional[list[dict[str, Any]]] = None


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


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
