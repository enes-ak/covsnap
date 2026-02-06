"""Output writers for raw TSV, JSON, exon TSV, and low-coverage BED."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from typing import Any, Optional

from covsnap import ANNOTATION_VERSION, BUILD, __version__
from covsnap.metrics import TargetResult


# ---------------------------------------------------------------------------
# Raw TSV
# ---------------------------------------------------------------------------

# Column order for the raw TSV
RAW_COLUMNS = [
    "target_id",
    "contig",
    "start",
    "end",
    "length_bp",
    "mean_depth",
    "median_depth",
    "min_depth",
    "max_depth",
    "stdev_depth",
    "pct_zero",
    # pct_ge_X columns are inserted dynamically
    "n_lowcov_blocks",
    "lowcov_total_bp",
    "engine_used",
    "bam_path",
    "sample_name",
    "build",
    "annotation_version",
]


def write_raw_tsv(
    path: str,
    results: list[TargetResult],
    thresholds: list[int],
    engine_used: str,
    bam_path: str,
    sample_name: str,
    run_date: Optional[str] = None,
) -> None:
    """Write per-target metrics to a raw TSV file.

    Parameters
    ----------
    path : str
        Output file path.
    results : list of TargetResult
        One result per target.
    thresholds : list of int
        Depth thresholds (must match what was used during computation).
    engine_used, bam_path, sample_name : str
        Metadata to include in each row and the header comment.
    run_date : str or None
        ISO-8601 date string; auto-generated if None.
    """
    if run_date is None:
        run_date = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    sorted_thresholds = sorted(thresholds)
    pct_columns = [f"pct_ge_{t}" for t in sorted_thresholds]

    # Build full header
    header_parts = list(RAW_COLUMNS[:11])  # up to and including pct_zero
    header_parts.extend(pct_columns)
    header_parts.extend(RAW_COLUMNS[11:])  # from n_lowcov_blocks onward

    with open(path, "w") as fh:
        # Metadata comment line
        fh.write(
            f"#covsnap_version={__version__}\t"
            f"annotation={ANNOTATION_VERSION}\t"
            f"build={BUILD}\t"
            f"engine={engine_used}\t"
            f"date={run_date}\n"
        )
        # Header
        fh.write("\t".join(header_parts) + "\n")

        # Data rows
        for r in results:
            row = [
                r.target_id,
                r.contig,
                str(r.start),
                str(r.end),
                str(r.length_bp),
                f"{r.mean_depth:.2f}",
                f"{r.median_depth:.1f}",
                str(r.min_depth),
                str(r.max_depth),
                f"{r.stdev_depth:.2f}",
                f"{r.pct_zero:.2f}",
            ]
            # Threshold percentages
            for t in sorted_thresholds:
                row.append(f"{r.pct_thresholds.get(t, 0.0):.2f}")
            # Remaining columns
            row.extend([
                str(r.n_lowcov_blocks),
                str(r.lowcov_total_bp),
                engine_used,
                bam_path,
                sample_name,
                BUILD,
                ANNOTATION_VERSION,
            ])
            fh.write("\t".join(row) + "\n")


# ---------------------------------------------------------------------------
# JSON output
# ---------------------------------------------------------------------------


def write_json(
    path: str,
    results: list[TargetResult],
    thresholds: list[int],
    engine_used: str,
    bam_path: str,
    sample_name: str,
    run_date: Optional[str] = None,
) -> None:
    """Write per-target metrics as JSON."""
    if run_date is None:
        run_date = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    sorted_thresholds = sorted(thresholds)

    targets_list: list[dict[str, Any]] = []
    for r in results:
        entry: dict[str, Any] = {
            "target_id": r.target_id,
            "contig": r.contig,
            "start": r.start,
            "end": r.end,
            "length_bp": r.length_bp,
            "mean_depth": r.mean_depth,
            "median_depth": r.median_depth,
            "min_depth": r.min_depth,
            "max_depth": r.max_depth,
            "stdev_depth": r.stdev_depth,
            "pct_zero": r.pct_zero,
        }
        for t in sorted_thresholds:
            entry[f"pct_ge_{t}"] = r.pct_thresholds.get(t, 0.0)
        entry["n_lowcov_blocks"] = r.n_lowcov_blocks
        entry["lowcov_total_bp"] = r.lowcov_total_bp
        entry["coverage_status"] = r.coverage_status
        entry["lowcov_blocks"] = [
            {"start": b.start, "end": b.end, "mean_depth": round(b.mean_depth, 2)}
            for b in r.lowcov_blocks
        ]
        targets_list.append(entry)

    doc = {
        "covsnap_version": __version__,
        "annotation": ANNOTATION_VERSION,
        "build": BUILD,
        "engine": engine_used,
        "date": run_date,
        "bam_path": bam_path,
        "sample_name": sample_name,
        "targets": targets_list,
    }

    with open(path, "w") as fh:
        json.dump(doc, fh, indent=2)
        fh.write("\n")


# ---------------------------------------------------------------------------
# Exon TSV
# ---------------------------------------------------------------------------

EXON_COLUMNS = [
    "target_id",
    "exon_id",
    "exon_number",
    "contig",
    "start",
    "end",
    "length_bp",
    "mean_depth",
    "median_depth",
    "pct_zero",
    "pct_ge_20",
    "pct_ge_30",
]


def write_exon_tsv(
    path: str,
    exon_results: list[TargetResult],
    gene_name: str,
    exon_metadata: list[dict[str, Any]],
    run_date: Optional[str] = None,
) -> None:
    """Write exon-level metrics TSV.

    *exon_metadata* provides the exon_id, exon_number, etc. for each exon,
    in the same order as *exon_results*.
    """
    if run_date is None:
        run_date = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    with open(path, "w") as fh:
        fh.write(
            f"#covsnap_version={__version__}\t"
            f"annotation={ANNOTATION_VERSION}\t"
            f"build={BUILD}\n"
        )
        fh.write("\t".join(EXON_COLUMNS) + "\n")

        for r, meta in zip(exon_results, exon_metadata):
            row = [
                gene_name,
                meta.get("exon_id", r.target_id),
                str(meta.get("exon_number", 0)),
                r.contig,
                str(r.start),
                str(r.end),
                str(r.length_bp),
                f"{r.mean_depth:.2f}",
                f"{r.median_depth:.1f}",
                f"{r.pct_zero:.2f}",
                f"{r.pct_thresholds.get(20, 0.0):.2f}",
                f"{r.pct_thresholds.get(30, 0.0):.2f}",
            ]
            fh.write("\t".join(row) + "\n")


# ---------------------------------------------------------------------------
# Low-coverage BED
# ---------------------------------------------------------------------------


def write_lowcov_bed(
    path: str,
    results: list[TargetResult],
    lowcov_threshold: int,
    lowcov_min_len: int,
) -> None:
    """Write low-coverage blocks as a BED file."""
    with open(path, "w") as fh:
        fh.write(
            f'#track name="covsnap_lowcov" '
            f'description="Low-coverage blocks (depth < {lowcov_threshold}, '
            f'min {lowcov_min_len}bp)"\n'
        )
        for r in results:
            for b in r.lowcov_blocks:
                fh.write(
                    f"{r.contig}\t{b.start}\t{b.end}\t"
                    f"{r.target_id}\tmean_depth={b.mean_depth:.1f}\n"
                )
