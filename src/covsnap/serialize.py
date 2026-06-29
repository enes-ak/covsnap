"""Machine-readable serializers for covsnap reports.

Mirrors html_report.write_html_report: every writer consumes a ReportContext
(the single source of truth) and performs no metric computation.
"""

from __future__ import annotations

import json

from covsnap import __version__
from covsnap.metrics import TargetResult
from covsnap.report import ReportContext


def _target_to_dict(r: TargetResult) -> dict:
    return {
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
        "pct_thresholds": {str(t): v for t, v in r.pct_thresholds.items()},
        "coverage_status": r.coverage_status,
        "n_lowcov_blocks": r.n_lowcov_blocks,
        "lowcov_total_bp": r.lowcov_total_bp,
        "lowcov_blocks": [
            {
                "contig": r.contig,
                "start": b.start,
                "end": b.end,
                "length": b.length,
                "mean_depth": round(b.mean_depth, 2),
            }
            for b in r.lowcov_blocks
        ],
        "histogram": {str(d): c for d, c in r.histogram.items()},
    }


def report_to_dict(ctx: ReportContext) -> dict:
    """Build the canonical native-JSON structure from a ReportContext."""
    targets = []
    for r in ctx.results:
        td = _target_to_dict(r)
        exon_list = (ctx.exon_results or {}).get(r.target_id)
        if exon_list:
            statuses = (ctx.exon_statuses or {}).get(r.target_id, [])
            metas = (ctx.exon_metadata or {}).get(r.target_id, [])
            exons = []
            for i, er in enumerate(exon_list):
                ed = _target_to_dict(er)
                meta = metas[i] if i < len(metas) else {}
                ed["exon_number"] = meta.get("exon_number", i + 1)
                ed["status"] = statuses[i] if i < len(statuses) else ""
                exons.append(ed)
            td["exons"] = exons
        targets.append(td)

    status_counts: dict[str, int] = {}
    for r in ctx.results:
        status_counts[r.coverage_status] = status_counts.get(r.coverage_status, 0) + 1

    genes = [_target_to_dict(g) for g in (ctx.gene_results or [])]

    p = ctx.classify_params
    return {
        "covsnap_version": __version__,
        "schema_version": "1.0",
        "run": {
            "sample_name": ctx.sample_name,
            "bam_path": ctx.bam_path,
            "reference": ctx.reference,
            "engine_used": ctx.engine_used,
            "engine_version": ctx.engine_version,
            "contig_style": ctx.contig_style,
            "run_date": ctx.run_date,
            "bed_path": ctx.bed_path,
            "exon_only": ctx.exon_only,
            "thresholds": list(ctx.thresholds),
            "classify_params": {
                "pass_pct_ge_20": p.pass_pct_ge_20,
                "pass_max_pct_zero": p.pass_max_pct_zero,
                "dropout_pct_zero": p.dropout_pct_zero,
                "dropout_zero_block_bp": p.dropout_zero_block_bp,
                "uneven_cv": p.uneven_cv,
                "exon_pct_ge_20": p.exon_pct_ge_20,
                "exon_max_pct_zero": p.exon_max_pct_zero,
            },
        },
        "summary": {
            "n_targets": len(ctx.results),
            "n_pass": sum(1 for r in ctx.results if r.coverage_status == "PASS"),
            "status_counts": status_counts,
        },
        "targets": targets,
        "genes": genes,
    }


def write_json_report(path: str, ctx: ReportContext) -> None:
    """Write the full native JSON report to *path*."""
    with open(path, "w") as f:
        json.dump(report_to_dict(ctx), f, indent=2)


# Severity order for worst_status aggregation (higher = more severe).
_STATUS_SEVERITY = {
    "PASS": 0,
    "LOW_COVERAGE": 1,
    "LOW_EXON": 2,
    "UNEVEN": 3,
    "DROP_OUT": 4,
}


def write_multiqc_report(path: str, ctx: ReportContext) -> None:
    """Write a MultiQC custom-content file (one sample-level summary row)."""
    results = ctx.results
    n_targets = len(results)
    n_pass = sum(1 for r in results if r.coverage_status == "PASS")
    total_len = sum(r.length_bp for r in results) or 1

    mean_depth = sum(r.mean_depth * r.length_bp for r in results) / total_len
    row = {
        "n_targets": n_targets,
        "n_pass": n_pass,
        "pct_targets_pass": round(n_pass / n_targets * 100, 2) if n_targets else 0.0,
        "mean_depth": round(mean_depth, 2),
    }
    if 20 in ctx.thresholds:
        pct20 = sum(r.pct_thresholds.get(20, 0.0) * r.length_bp for r in results) / total_len
        row["pct_ge_20"] = round(pct20, 2)
    row["worst_status"] = (
        max((r.coverage_status for r in results), key=lambda s: _STATUS_SEVERITY.get(s, 0))
        if results
        else ""
    )

    payload = {
        "id": "covsnap",
        "section_name": "covsnap Coverage QC",
        "description": "Per-sample targeted-sequencing coverage QC summary from covsnap.",
        "plot_type": "table",
        "pconfig": {
            "id": "covsnap_table",
            "title": "covsnap: Coverage QC",
            "namespace": "covsnap",
        },
        "data": {ctx.sample_name: row},
    }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)


def write_tsv_report(path: str, ctx: ReportContext) -> None:
    """Write a flat, target-level TSV (one row per target; exons excluded)."""
    thresholds = sorted(ctx.thresholds)
    columns = [
        "sample_name", "target_id", "contig", "start", "end", "length_bp",
        "mean_depth", "median_depth", "min_depth", "max_depth", "stdev_depth", "pct_zero",
    ]
    columns += [f"pct_ge_{t}" for t in thresholds]
    columns += ["coverage_status", "n_lowcov_blocks", "lowcov_total_bp"]

    with open(path, "w", encoding="utf-8") as f:
        f.write("\t".join(columns) + "\n")
        for r in ctx.results:
            row = [
                ctx.sample_name, r.target_id, r.contig, str(r.start), str(r.end),
                str(r.length_bp), str(r.mean_depth), str(r.median_depth),
                str(r.min_depth), str(r.max_depth), str(r.stdev_depth), str(r.pct_zero),
            ]
            row += [str(r.pct_thresholds.get(t, 0.0)) for t in thresholds]
            row += [r.coverage_status, str(r.n_lowcov_blocks), str(r.lowcov_total_bp)]
            f.write("\t".join(row) + "\n")
