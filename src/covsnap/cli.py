"""CLI entry point for covsnap.

Orchestrates: argument parsing → input validation → annotation lookup →
engine execution → classification → HTML report writing.
"""

from __future__ import annotations

import argparse
import logging
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from typing import Any, Optional

from covsnap import ANNOTATION_VERSION, BUILD, __version__
from covsnap.annotation import (
    detect_contig_style,
    get_sample_name,
    is_region_string,
    lookup_exons,
    lookup_gene,
    lookup_genes_in_region,
    parse_region,
    suggest_genes,
    translate_contig,
)
from covsnap.bed import enforce_limits
from covsnap.engines import compute_depth, ensure_index, select_engine
from covsnap.html_report import write_html_report
from covsnap.metrics import TargetResult
from covsnap.report import (
    ClassifyParams,
    ReportContext,
    classify_exon,
    classify_target,
)

logger = logging.getLogger("covsnap")


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="covsnap",
        description=(
            "covsnap — Coverage inspector for targeted sequencing QC (hg38 only).\n\n"
            "Computes per-target and optionally per-exon depth metrics from BAM/CRAM\n"
            "files, producing an interactive HTML report with PASS/FAIL classifications."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  covsnap sample.bam BRCA1                              # gene mode\n"
            "  covsnap sample.bam BRCA1 --exons                      # gene + exon breakdown\n"
            "  covsnap sample.bam chr17:43044295-43125482             # region mode (auto-finds genes & exons)\n"
            "  covsnap sample.bam --bed targets.bed                   # BED mode\n"
            "  covsnap sample.bam BRCA1 -o my_report.html             # custom output path\n"
            "  covsnap sample.cram BRCA1 --reference hg38.fa          # CRAM input\n"
            "  covsnap sample.bam BRCA1 --emit-lowcov                 # include low-coverage blocks\n"
            "  covsnap sample.bam BRCA1 --engine samtools --threads 8 # engine selection\n"
        ),
    )

    parser.add_argument("--version", action="version", version=f"covsnap {__version__}")

    # ── Positional ──
    parser.add_argument(
        "alignment",
        help="Path to BAM or CRAM file.",
    )
    parser.add_argument(
        "target",
        nargs="?",
        default=None,
        help=(
            "Gene symbol (e.g. BRCA1) or genomic region "
            "(e.g. chr17:43044295-43125482). Mutually exclusive with --bed."
        ),
    )

    # ── Input options ──
    inp = parser.add_argument_group("Input options")
    inp.add_argument(
        "--bed",
        metavar="BED",
        default=None,
        help="Path to a BED file of target regions. Mutually exclusive with positional target.",
    )
    inp.add_argument(
        "--exons",
        action="store_true",
        help="Enable exon-level statistics (gene mode only).",
    )
    inp.add_argument(
        "--reference",
        metavar="FASTA",
        default=None,
        help="Reference FASTA for CRAM decoding.",
    )
    inp.add_argument(
        "--no-index",
        action="store_true",
        help="Do not create .bai/.csi/.crai index if missing.",
    )

    # ── Engine options ──
    eng = parser.add_argument_group("Engine options")
    eng.add_argument(
        "--engine",
        choices=["auto", "mosdepth", "samtools"],
        default="auto",
        help="Depth computation engine (default: auto).",
    )
    eng.add_argument(
        "--threads",
        type=int,
        default=4,
        metavar="N",
        help="Parallel workers for samtools / threads for mosdepth (default: 4).",
    )

    # ── Output options ──
    out = parser.add_argument_group("Output options")
    out.add_argument(
        "-o", "--output",
        metavar="FILE",
        default="covsnap.report.html",
        help="HTML report output path (default: covsnap.report.html).",
    )

    # ── Low-coverage options ──
    lc = parser.add_argument_group("Low-coverage options")
    lc.add_argument(
        "--emit-lowcov",
        action="store_true",
        help="Include low-coverage blocks in the report.",
    )
    lc.add_argument(
        "--lowcov-threshold",
        type=int,
        default=10,
        help="Depth below which a position is 'low coverage' (default: 10).",
    )
    lc.add_argument(
        "--lowcov-min-len",
        type=int,
        default=50,
        help="Minimum contiguous length (bp) for a low-coverage block (default: 50).",
    )

    # ── BED guardrails ──
    bg = parser.add_argument_group("BED guardrail options")
    bg.add_argument(
        "--max-targets",
        type=int,
        default=2000,
        help="Maximum number of BED target intervals (default: 2000).",
    )
    bg.add_argument(
        "--max-total-bp",
        type=int,
        default=50_000_000,
        help="Maximum total target bases (default: 50000000).",
    )
    bg.add_argument(
        "--max-bed-bytes",
        type=_parse_size,
        default=50 * 1024 * 1024,
        metavar="BYTES",
        help="Maximum BED file size in bytes (default: 50MB). Accepts suffixes: K, M, G.",
    )
    bg.add_argument(
        "--on-large-bed",
        choices=["error", "warn_and_clip", "warn_and_sample"],
        default="warn_and_clip",
        help="Action when BED limits are exceeded (default: warn_and_clip).",
    )
    bg.add_argument(
        "--large-bed-seed",
        type=int,
        default=42,
        help="Random seed for warn_and_sample mode (default: 42).",
    )

    # ── Coverage thresholds ──
    ct = parser.add_argument_group("Coverage thresholds")
    ct.add_argument(
        "--pct-thresholds",
        default="1,5,10,20,30,50,100",
        help="Comma-separated depth thresholds for pct_ge_X columns (default: 1,5,10,20,30,50,100).",
    )

    # ── Classification tuning ──
    cl = parser.add_argument_group("Classification tuning")
    cl.add_argument("--pass-pct-ge-20", type=float, default=95.0,
                    help="Minimum pct_ge_20 for PASS (default: 95.0).")
    cl.add_argument("--pass-max-pct-zero", type=float, default=1.0,
                    help="Maximum pct_zero for PASS (default: 1.0).")
    cl.add_argument("--dropout-pct-zero", type=float, default=5.0,
                    help="pct_zero above which DROP_OUT is called (default: 5.0).")
    cl.add_argument("--uneven-cv", type=float, default=1.0,
                    help="CV above which UNEVEN is called (default: 1.0).")
    cl.add_argument("--exon-pct-ge-20", type=float, default=90.0,
                    help="Minimum exon pct_ge_20 for LOW_EXON (default: 90.0).")
    cl.add_argument("--exon-max-pct-zero", type=float, default=5.0,
                    help="Maximum exon pct_zero for LOW_EXON (default: 5.0).")

    # ── General ──
    gen = parser.add_argument_group("General")
    gen.add_argument(
        "-v", "--verbose",
        action="count",
        default=0,
        help="Increase logging verbosity (can repeat: -vv).",
    )
    gen.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress all non-error output to stderr.",
    )

    return parser


def _parse_size(value: str) -> int:
    """Parse a size string like '50M' or '1G' into bytes."""
    value = value.strip().upper()
    multipliers = {"K": 1024, "M": 1024**2, "G": 1024**3}
    if value[-1] in multipliers:
        return int(float(value[:-1]) * multipliers[value[-1]])
    return int(value)


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def _validate_args(args: argparse.Namespace) -> None:
    """Validate CLI arguments; exit with a clear error on failure."""
    # Mutual exclusion: target vs --bed
    if args.target is not None and args.bed is not None:
        _error("Provide either a positional target (gene/region) or --bed, not both.", code=1)

    if args.target is None and args.bed is None:
        _error(
            "Must provide a target: gene symbol, region string (positional), or --bed file.",
            code=1,
        )

    # Alignment file exists
    if not os.path.isfile(args.alignment):
        _error(f"Alignment file not found: {args.alignment}", code=1)

    # BED file exists
    if args.bed is not None and not os.path.isfile(args.bed):
        _error(f"BED file not found: {args.bed}", code=1)

    # CRAM reference check
    if args.alignment.endswith(".cram"):
        if args.reference is None and not os.environ.get("REF_PATH") and not os.environ.get("REF_CACHE"):
            _error(
                "CRAM input requires --reference <fasta> or REF_PATH/REF_CACHE environment variable.",
                code=5,
            )

    # Threads
    if args.threads < 1:
        _error("--threads must be >= 1.", code=1)

    # Parse thresholds
    try:
        thresholds = [int(t.strip()) for t in args.pct_thresholds.split(",")]
        if any(t < 0 for t in thresholds):
            raise ValueError("negative threshold")
        args._parsed_thresholds = thresholds
    except ValueError:
        _error("--pct-thresholds must be comma-separated positive integers.", code=1)

    # Exons mode requires gene target (not region or bed)
    if args.exons and args.bed is not None:
        _error("--exons is only supported in gene mode (positional gene target), not with --bed.", code=1)


def _error(message: str, code: int = 1) -> None:
    """Print an error to stderr and exit."""
    print(f"[covsnap] ERROR: {message}", file=sys.stderr)
    sys.exit(code)


# ---------------------------------------------------------------------------
# Engine version detection
# ---------------------------------------------------------------------------


def _get_engine_version(engine: str) -> str:
    """Get version string for the active engine."""
    try:
        if engine == "mosdepth":
            proc = subprocess.run(
                ["mosdepth", "--version"],
                capture_output=True, text=True,
            )
            # mosdepth outputs: "mosdepth X.Y.Z"
            return proc.stdout.strip().split()[-1] if proc.stdout.strip() else ""
        elif engine == "samtools":
            proc = subprocess.run(
                ["samtools", "--version"],
                capture_output=True, text=True,
            )
            # First line: "samtools X.Y"
            first_line = proc.stdout.strip().split("\n")[0] if proc.stdout.strip() else ""
            return first_line.split()[-1] if first_line else ""
    except Exception:
        pass
    return ""


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main(argv: Optional[list[str]] = None) -> None:
    """CLI entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)

    # ── Logging setup ──
    if args.quiet:
        log_level = logging.ERROR
    elif args.verbose >= 2:
        log_level = logging.DEBUG
    elif args.verbose == 1:
        log_level = logging.INFO
    else:
        log_level = logging.WARNING

    logging.basicConfig(
        level=log_level,
        format="[covsnap] %(levelname)s: %(message)s",
        stream=sys.stderr,
    )

    # ── Validate ──
    _validate_args(args)
    thresholds: list[int] = args._parsed_thresholds  # type: ignore[attr-defined]

    run_date = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    # ── Detect contig style ──
    try:
        contig_style = detect_contig_style(args.alignment)
    except ValueError as exc:
        _error(str(exc), code=1)
        return  # unreachable, but helps type checker

    logger.info("Detected contig style: %s", contig_style)

    # ── Select engine ──
    try:
        engine = select_engine(args.engine)
    except RuntimeError as exc:
        _error(str(exc), code=2)
        return

    engine_version = _get_engine_version(engine)
    logger.info("Using engine: %s %s", engine, engine_version)

    # ── Ensure index ──
    try:
        ensure_index(args.alignment, no_index=args.no_index)
    except RuntimeError as exc:
        _error(str(exc), code=1)
        return

    # ── Get sample name ──
    sample_name = get_sample_name(args.alignment)

    # ── Resolve targets ──
    # regions: list of (contig, start, end, name) in BAM contig style
    regions: list[tuple[str, int, int, str]] = []
    gene_name_resolved: Optional[str] = None
    bed_guardrail_message = ""

    if args.bed is not None:
        # BED mode
        result = enforce_limits(
            args.bed,
            max_targets=args.max_targets,
            max_total_bp=args.max_total_bp,
            max_bed_bytes=args.max_bed_bytes,
            on_large_bed=args.on_large_bed,
            seed=args.large_bed_seed,
        )
        bed_guardrail_message = result.message
        for iv in result.intervals:
            # BED is already 0-based half-open
            bam_contig = translate_contig(iv.contig, contig_style)
            regions.append((bam_contig, iv.start, iv.end, iv.name))

    elif args.target is not None:
        if is_region_string(args.target):
            # Region mode
            try:
                contig, start, end = parse_region(args.target)
            except ValueError as exc:
                _error(str(exc), code=1)
                return
            bam_contig = translate_contig(contig, contig_style)
            region_name = f"{contig}:{start + 1}-{end}"
            regions.append((bam_contig, start, end, region_name))
        elif ":" in args.target:
            # Contains ':' but didn't match region regex — malformed region
            _error(
                f"'{args.target}' looks like a region string but does not match "
                "the expected format 'contig:start-end' (e.g. chr17:43044295-43125482).",
                code=1,
            )
            return
        else:
            # Gene mode
            gene_info = lookup_gene(args.target)
            if gene_info is None:
                suggestions = suggest_genes(args.target)
                suggestion_str = ", ".join(suggestions) if suggestions else "(no close matches)"
                _error(
                    f"Gene '{args.target}' not found in hg38 gene index ({ANNOTATION_VERSION}). "
                    f"Did you mean: {suggestion_str}?",
                    code=3,
                )
                return
            gene_name_resolved = gene_info["gene_name"]
            bam_contig = translate_contig(gene_info["contig"], contig_style)
            regions.append((
                bam_contig,
                gene_info["start"],
                gene_info["end"],
                gene_name_resolved,
            ))

    if not regions:
        _error("No target regions resolved. Nothing to do.", code=1)
        return

    logger.info("Processing %d target region(s)", len(regions))

    # ── Compute depth metrics ──
    try:
        results = compute_depth(
            bam_path=args.alignment,
            regions=regions,
            engine=engine,
            thresholds=thresholds,
            lowcov_threshold=args.lowcov_threshold,
            lowcov_min_len=args.lowcov_min_len,
            threads=args.threads,
            reference=args.reference,
        )
    except RuntimeError as exc:
        _error(str(exc), code=2)
        return

    # Attach shared metadata
    for r in results:
        r.engine_used = engine
        r.bam_path = args.alignment
        r.sample_name = sample_name

    # ── Exon-level analysis (gene mode only) ──
    exon_results_map: Optional[dict[str, list[TargetResult]]] = None
    exon_metadata_map: Optional[dict[str, list[dict[str, Any]]]] = None
    exon_status_map: Optional[dict[str, list[str]]] = None

    if args.exons and gene_name_resolved:
        exon_records = lookup_exons(gene_name_resolved)
        if exon_records:
            exon_regions = [
                (
                    translate_contig(e["contig"], contig_style),
                    e["start"],
                    e["end"],
                    e.get("exon_id", f"exon_{e['exon_number']}"),
                )
                for e in exon_records
            ]
            try:
                exon_depth_results = compute_depth(
                    bam_path=args.alignment,
                    regions=exon_regions,
                    engine=engine,
                    thresholds=thresholds,
                    lowcov_threshold=args.lowcov_threshold,
                    lowcov_min_len=args.lowcov_min_len,
                    threads=args.threads,
                    reference=args.reference,
                )
                for er in exon_depth_results:
                    er.engine_used = engine
                    er.bam_path = args.alignment
                    er.sample_name = sample_name

                exon_results_map = {gene_name_resolved: exon_depth_results}
                exon_metadata_map = {gene_name_resolved: exon_records}

                # Classify exons
                params = _build_classify_params(args)
                exon_statuses = [classify_exon(er, params) for er in exon_depth_results]
                exon_status_map = {gene_name_resolved: exon_statuses}
            except RuntimeError as exc:
                logger.warning("Exon-level analysis failed: %s", exc)
        else:
            logger.warning(
                "No exon records found for %s. Exon index may not be installed.",
                gene_name_resolved,
            )

    # ── Region → gene/exon overlay (region mode) ──
    gene_results_list: Optional[list[TargetResult]] = None
    gene_metadata_list: Optional[list[dict[str, Any]]] = None

    is_region_mode = (args.target is not None and is_region_string(args.target)
                      and gene_name_resolved is None)

    if is_region_mode and len(regions) == 1:
        reg_contig, reg_start, reg_end, _ = regions[0]
        # Find overlapping genes using chr-style contig for annotation lookup
        chr_contig = translate_contig(reg_contig, "chr")
        overlapping_genes = lookup_genes_in_region(chr_contig, reg_start, reg_end)

        if overlapping_genes:
            logger.info("Found %d gene(s) in region", len(overlapping_genes))

            # Compute gene-level coverage
            gene_regions = [
                (
                    translate_contig(g["contig"], contig_style),
                    max(g["start"], reg_start),  # clip to input region
                    min(g["end"], reg_end),
                    g["gene_name"],
                )
                for g in overlapping_genes
            ]
            try:
                gene_depth_results = compute_depth(
                    bam_path=args.alignment,
                    regions=gene_regions,
                    engine=engine,
                    thresholds=thresholds,
                    lowcov_threshold=args.lowcov_threshold,
                    lowcov_min_len=args.lowcov_min_len,
                    threads=args.threads,
                    reference=args.reference,
                )
                for gr in gene_depth_results:
                    gr.engine_used = engine
                    gr.bam_path = args.alignment
                    gr.sample_name = sample_name

                gene_results_list = gene_depth_results
                gene_metadata_list = overlapping_genes

                # Per-gene exon analysis (parallel)
                exon_results_map = exon_results_map or {}
                exon_metadata_map = exon_metadata_map or {}
                exon_status_map = exon_status_map or {}
                params = _build_classify_params(args)

                # Prepare exon jobs
                exon_jobs: list[tuple[str, list[tuple[str, int, int, str]], list[dict[str, Any]]]] = []
                for g_info in overlapping_genes:
                    gname = g_info["gene_name"]
                    exon_records = lookup_exons(gname)
                    if not exon_records:
                        continue
                    clipped_exon_regions = []
                    clipped_exon_records = []
                    for e in exon_records:
                        e_start = max(e["start"], reg_start)
                        e_end = min(e["end"], reg_end)
                        if e_start < e_end:
                            clipped_exon_regions.append((
                                translate_contig(e["contig"], contig_style),
                                e_start, e_end,
                                e.get("exon_id", f"exon_{e['exon_number']}"),
                            ))
                            clipped_exon_records.append(e)
                    if clipped_exon_regions:
                        exon_jobs.append((gname, clipped_exon_regions, clipped_exon_records))

                if exon_jobs:
                    n_workers = min(args.threads, len(exon_jobs))
                    with ThreadPoolExecutor(max_workers=n_workers) as pool:
                        futures = {}
                        for gname, exon_regs, exon_recs in exon_jobs:
                            future = pool.submit(
                                compute_depth,
                                bam_path=args.alignment,
                                regions=exon_regs,
                                engine=engine,
                                thresholds=thresholds,
                                lowcov_threshold=args.lowcov_threshold,
                                lowcov_min_len=args.lowcov_min_len,
                                threads=1,  # each worker gets 1 thread
                                reference=args.reference,
                            )
                            futures[future] = (gname, exon_recs)

                        for future in as_completed(futures):
                            gname, exon_recs = futures[future]
                            try:
                                exon_depth = future.result()
                                for er in exon_depth:
                                    er.engine_used = engine
                                    er.bam_path = args.alignment
                                    er.sample_name = sample_name
                                exon_results_map[gname] = exon_depth
                                exon_metadata_map[gname] = exon_recs
                                exon_status_map[gname] = [classify_exon(er, params) for er in exon_depth]
                            except RuntimeError as exc:
                                logger.warning("Exon analysis failed for %s: %s", gname, exc)
            except RuntimeError as exc:
                logger.warning("Gene-level analysis failed: %s", exc)
        else:
            logger.info("No genes found in region %s:%d-%d", chr_contig, reg_start, reg_end)

    # ── Classify targets ──
    params = _build_classify_params(args)
    for r in results:
        exon_list = None
        if exon_results_map and r.target_id in exon_results_map:
            exon_list = exon_results_map[r.target_id]
        classify_target(r, params, exon_results=exon_list)

    # Classify genes if present
    if gene_results_list:
        for gr in gene_results_list:
            gene_exons = exon_results_map.get(gr.target_id) if exon_results_map else None
            classify_target(gr, params, exon_results=gene_exons)

    # ── Write HTML report ──
    report_ctx = ReportContext(
        results=results,
        engine_used=engine,
        engine_version=engine_version,
        bam_path=args.alignment,
        sample_name=sample_name,
        contig_style=contig_style,
        reference=args.reference,
        bed_path=args.bed,
        bed_guardrail_message=bed_guardrail_message,
        exon_results=exon_results_map if exon_results_map else None,
        exon_metadata=exon_metadata_map if exon_metadata_map else None,
        exon_statuses=exon_status_map if exon_status_map else None,
        classify_params=params,
        lowcov_threshold=args.lowcov_threshold,
        lowcov_min_len=args.lowcov_min_len,
        emit_lowcov=args.emit_lowcov,
        thresholds=thresholds,
        run_date=run_date,
        gene_results=gene_results_list,
        gene_metadata=gene_metadata_list,
    )
    write_html_report(args.output, report_ctx)
    logger.info("HTML report written to %s", args.output)

    # ── Summary to stderr ──
    if not args.quiet:
        n_pass = sum(1 for r in results if r.coverage_status == "PASS")
        n_total = len(results)
        print(
            f"[covsnap] Done. {n_pass}/{n_total} targets PASS. "
            f"Report: {args.output}",
            file=sys.stderr,
        )


def _build_classify_params(args: argparse.Namespace) -> ClassifyParams:
    """Build ClassifyParams from CLI arguments."""
    return ClassifyParams(
        pass_pct_ge_20=args.pass_pct_ge_20,
        pass_max_pct_zero=args.pass_max_pct_zero,
        dropout_pct_zero=args.dropout_pct_zero,
        uneven_cv=args.uneven_cv,
        exon_pct_ge_20=args.exon_pct_ge_20,
        exon_max_pct_zero=args.exon_max_pct_zero,
    )


if __name__ == "__main__":
    main()
