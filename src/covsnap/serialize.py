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
