"""Depth computation engines: samtools depth and mosdepth.

Both engines stream per-base depth data and feed it into
:class:`~covsnap.metrics.TargetAccumulator` instances,
so RAM usage is O(targets) regardless of region size.
"""

from __future__ import annotations

import gzip
import logging
import math
import os
import shutil
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Optional

from covsnap.metrics import LowCovBlock, TargetAccumulator, TargetResult

logger = logging.getLogger(__name__)


def select_engine(engine: str) -> str:
    """Resolve an engine name to a concrete engine.

    Parameters
    ----------
    engine : str
        One of ``"auto"``, ``"mosdepth"``, ``"samtools"``, ``"pysam"``.

    Returns
    -------
    str
        The selected engine name.

    Raises
    ------
    RuntimeError
        If the requested engine (or any engine in auto mode) is not found.
    """
    if engine == "pysam":
        return "pysam"

    if engine == "mosdepth":
        if shutil.which("mosdepth"):
            return "mosdepth"
        raise RuntimeError("mosdepth not found on PATH. Install it or use --engine samtools.")

    if engine == "samtools":
        if shutil.which("samtools"):
            return "samtools"
        raise RuntimeError("samtools not found on PATH.")

    # auto mode: pysam first (fastest, no subprocess), then mosdepth, then samtools
    logger.info("Auto-selected engine: pysam")
    return "pysam"


# ---------------------------------------------------------------------------
# Target specification helpers
# ---------------------------------------------------------------------------


class _TargetSpec:
    """Lightweight container for a region to analyse."""

    __slots__ = ("target_id", "contig", "start", "end")

    def __init__(self, target_id: str, contig: str, start: int, end: int):
        self.target_id = target_id
        self.contig = contig
        self.start = start
        self.end = end

    @property
    def length(self) -> int:
        return self.end - self.start


def _make_target_specs(
    regions: list[tuple[str, int, int, str]],
) -> list[_TargetSpec]:
    """Convert (contig, start, end, name) tuples to _TargetSpec list."""
    return [_TargetSpec(target_id=name, contig=contig, start=start, end=end) for contig, start, end, name in regions]


def _write_tmp_bed(targets: list[_TargetSpec], tmp_dir: str) -> str:
    """Write targets to a temporary BED file, return path."""
    bed_path = os.path.join(tmp_dir, "covsnap_targets.bed")
    with open(bed_path, "w") as fh:
        for t in targets:
            fh.write(f"{t.contig}\t{t.start}\t{t.end}\t{t.target_id}\n")
    return bed_path


# ---------------------------------------------------------------------------
# Public compute API
# ---------------------------------------------------------------------------


def compute_depth(
    bam_path: str,
    regions: list[tuple[str, int, int, str]],
    engine: str,
    thresholds: list[int],
    lowcov_threshold: int = 10,
    lowcov_min_len: int = 50,
    threads: int = 4,
    tmp_dir: Optional[str] = None,
    reference: Optional[str] = None,
) -> list[TargetResult]:
    """Compute depth metrics for the given regions.

    Parameters
    ----------
    bam_path : str
        Path to BAM/CRAM file (must be indexed).
    regions : list of (contig, start, end, name)
        Target regions in 0-based half-open coordinates.
    engine : str
        ``"samtools"`` or ``"mosdepth"``.
    thresholds : list of int
        Depth thresholds for pct_ge_X columns.
    lowcov_threshold : int
        Depth below which a position is "low coverage".
    lowcov_min_len : int
        Minimum contiguous length for a low-cov block.
    threads : int
        Number of parallel workers for samtools, threads for mosdepth.
    tmp_dir : str or None
        Temporary directory (auto-created if None).
    reference : str or None
        Reference FASTA for CRAM files.

    Returns
    -------
    list of TargetResult
        One result per region, in input order.
    """
    specs = _make_target_specs(regions)
    cleanup_tmp = tmp_dir is None
    if tmp_dir is None:
        tmp_dir = tempfile.mkdtemp(prefix="covsnap_")

    try:
        if engine == "pysam":
            return _run_pysam(
                bam_path,
                specs,
                thresholds,
                lowcov_threshold,
                lowcov_min_len,
                reference,
            )
        elif engine == "samtools":
            return _run_samtools_parallel(
                bam_path,
                specs,
                thresholds,
                lowcov_threshold,
                lowcov_min_len,
                threads,
                tmp_dir,
                reference,
            )
        elif engine == "mosdepth":
            return _run_mosdepth(
                bam_path,
                specs,
                thresholds,
                lowcov_threshold,
                lowcov_min_len,
                threads,
                tmp_dir,
                reference,
            )
        else:
            raise ValueError(f"Unknown engine: {engine}")
    finally:
        if cleanup_tmp:
            _cleanup_dir(tmp_dir)


def _cleanup_dir(path: str) -> None:
    """Remove a temporary directory and its contents."""
    try:
        for f in os.listdir(path):
            os.remove(os.path.join(path, f))
        os.rmdir(path)
    except OSError:
        pass


# ---------------------------------------------------------------------------
# samtools depth engine — parallel
# ---------------------------------------------------------------------------


def _run_samtools_parallel(
    bam_path: str,
    targets: list[_TargetSpec],
    thresholds: list[int],
    lowcov_threshold: int,
    lowcov_min_len: int,
    threads: int,
    tmp_dir: str,
    reference: Optional[str],
) -> list[TargetResult]:
    """Run samtools depth in parallel chunks for multiple regions.

    When there are multiple targets and threads > 1, splits targets
    into chunks and runs a separate samtools depth process per chunk.
    Each chunk handles a subset of non-overlapping regions, and results
    are merged back in original order.
    """
    n_targets = len(targets)
    workers = min(threads, n_targets)

    if workers <= 1 or n_targets <= 1:
        return _run_samtools_chunk(
            bam_path,
            targets,
            thresholds,
            lowcov_threshold,
            lowcov_min_len,
            tmp_dir,
            reference,
        )

    # Split targets into chunks, preserving original indices
    indexed_targets = list(enumerate(targets))
    chunk_size = math.ceil(n_targets / workers)
    chunks = [indexed_targets[i : i + chunk_size] for i in range(0, n_targets, chunk_size)]

    logger.info("Parallel samtools: %d chunks across %d workers", len(chunks), workers)

    results = [None] * n_targets  # type: ignore[list-item]

    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {}
        for ci, chunk in enumerate(chunks):
            chunk_specs = [t for _, t in chunk]
            chunk_indices = [i for i, _ in chunk]
            chunk_tmp = os.path.join(tmp_dir, f"chunk_{ci}")
            os.makedirs(chunk_tmp, exist_ok=True)
            # Serialize specs for subprocess
            chunk_regions = [(s.contig, s.start, s.end, s.target_id) for s in chunk_specs]
            future = pool.submit(
                _samtools_chunk_worker,
                bam_path,
                chunk_regions,
                thresholds,
                lowcov_threshold,
                lowcov_min_len,
                chunk_tmp,
                reference,
            )
            futures[future] = chunk_indices

        for future in as_completed(futures):
            chunk_indices = futures[future]
            chunk_results = future.result()
            for idx, result in zip(chunk_indices, chunk_results):
                results[idx] = result

    return results  # type: ignore[return-value]


def _samtools_chunk_worker(
    bam_path: str,
    regions: list[tuple[str, int, int, str]],
    thresholds: list[int],
    lowcov_threshold: int,
    lowcov_min_len: int,
    tmp_dir: str,
    reference: Optional[str],
) -> list[TargetResult]:
    """Worker function for parallel samtools execution."""
    specs = _make_target_specs(regions)
    return _run_samtools_chunk(
        bam_path,
        specs,
        thresholds,
        lowcov_threshold,
        lowcov_min_len,
        tmp_dir,
        reference,
    )


def _run_samtools_chunk(
    bam_path: str,
    targets: list[_TargetSpec],
    thresholds: list[int],
    lowcov_threshold: int,
    lowcov_min_len: int,
    tmp_dir: str,
    reference: Optional[str],
) -> list[TargetResult]:
    """Stream ``samtools depth -a`` output and accumulate metrics for one chunk."""
    # Sort targets by contig and start for cursor-based parsing
    sorted_targets = sorted(targets, key=lambda t: (t.contig, t.start))
    # Keep original order mapping for output
    original_order = {id(t): i for i, t in enumerate(targets)}

    bed_path = _write_tmp_bed(sorted_targets, tmp_dir)

    cmd = ["samtools", "depth", "-a", "-b", bed_path, bam_path]
    if reference:
        cmd.extend(["--reference", reference])

    logger.info("Running: %s", " ".join(cmd))
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,  # line-buffered
    )

    # Create accumulators for each target
    accumulators = [
        TargetAccumulator(
            target_id=t.target_id,
            contig=t.contig,
            start=t.start,
            end=t.end,
            thresholds=thresholds,
            lowcov_threshold=lowcov_threshold,
            lowcov_min_len=lowcov_min_len,
        )
        for t in sorted_targets
    ]

    cursor = 0
    n_targets = len(sorted_targets)

    assert proc.stdout is not None
    for line in proc.stdout:
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        contig = parts[0]
        pos_1based = int(parts[1])
        depth = int(parts[2])
        pos = pos_1based - 1  # convert to 0-based

        # Advance cursor past completed targets
        while cursor < n_targets:
            t = sorted_targets[cursor]
            if contig == t.contig and pos < t.end:
                break
            if contig != t.contig or pos >= t.end:
                cursor += 1
                continue
            break

        if cursor >= n_targets:
            break

        t = sorted_targets[cursor]
        if contig == t.contig and t.start <= pos < t.end:
            accumulators[cursor].add_single(depth, pos)

    proc.wait()
    stderr_output = proc.stderr.read() if proc.stderr else ""
    if proc.returncode != 0:
        logger.warning("samtools depth exited with code %d: %s", proc.returncode, stderr_output)

    # Finalize and reorder to original input order
    results_sorted = [acc.finalize() for acc in accumulators]

    # Map back to original order
    results = [TargetResult] * len(targets)  # type: ignore[list-item]
    for st, res in zip(sorted_targets, results_sorted):
        orig_idx = original_order[id(st)]
        results[orig_idx] = res

    return results  # type: ignore[return-value]


# ---------------------------------------------------------------------------
# mosdepth engine
# ---------------------------------------------------------------------------


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
        "--by",
        bed_path,
        "--threads",
        str(threads),
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
        perbase_file,
        targets,
        thresholds,
        lowcov_threshold,
        lowcov_min_len,
    )


def _stream_mosdepth_perbase(
    path: str,
    targets: list[_TargetSpec],
    thresholds: list[int],
    lowcov_threshold: int,
    lowcov_min_len: int,
) -> list[TargetResult]:
    """Stream mosdepth per-base.bed.gz to compute min/max/median/stdev/lowcov.

    The per-base file is run-length encoded: each line is
    ``chrom  start  end  depth`` representing a contiguous block of
    bases with the same depth.  We feed these blocks into
    TargetAccumulators.
    """
    sorted_targets = sorted(enumerate(targets), key=lambda x: (x[1].contig, x[1].start))

    accumulators = [
        TargetAccumulator(
            target_id=t.target_id,
            contig=t.contig,
            start=t.start,
            end=t.end,
            thresholds=thresholds,
            lowcov_threshold=lowcov_threshold,
            lowcov_min_len=lowcov_min_len,
        )
        for _, t in sorted_targets
    ]

    if not os.path.exists(path):
        logger.warning("mosdepth per-base file not found: %s", path)
        return [acc.finalize() for acc in accumulators]

    cursor = 0
    n_targets = len(sorted_targets)

    with gzip.open(path, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            contig = parts[0]
            block_start = int(parts[1])
            block_end = int(parts[2])
            depth = int(float(parts[3]))
            block_end - block_start

            # Find overlapping targets
            while cursor < n_targets:
                _, t = sorted_targets[cursor]
                if contig == t.contig and block_start < t.end:
                    break
                if contig != t.contig or block_start >= t.end:
                    cursor += 1
                    continue
                break

            if cursor >= n_targets:
                break

            # Check overlaps with current and subsequent targets
            for k in range(cursor, n_targets):
                _, t = sorted_targets[k]
                if t.contig != contig or t.start >= block_end:
                    break
                # Clip to target
                clip_start = max(block_start, t.start)
                clip_end = min(block_end, t.end)
                if clip_start < clip_end:
                    accumulators[k].add_block(depth, clip_end - clip_start, clip_start)

    # Finalize and reorder to original order
    finalized = [acc.finalize() for acc in accumulators]
    results = [None] * len(targets)  # type: ignore[list-item]
    for (orig_idx, _), res in zip(sorted_targets, finalized):
        results[orig_idx] = res

    return results  # type: ignore[return-value]


# ---------------------------------------------------------------------------
# pysam engine — direct API, no subprocess
# ---------------------------------------------------------------------------


def _run_pysam(
    bam_path: str,
    targets: list[_TargetSpec],
    thresholds: list[int],
    lowcov_threshold: int,
    lowcov_min_len: int,
    reference: Optional[str],
) -> list[TargetResult]:
    """Compute depth via pysam count_coverage (no subprocess needed)."""
    import numpy as np
    import pysam

    sorted_thresholds = sorted(thresholds)

    open_kwargs: dict = {"mode": "rb"}
    if reference:
        open_kwargs["reference_filename"] = reference

    try:
        af = pysam.AlignmentFile(bam_path, **open_kwargs)
    except Exception as exc:
        raise RuntimeError(f"pysam failed to open {bam_path}: {exc}") from exc

    results: list[TargetResult] = []
    try:
        for t in targets:
            try:
                cov = af.count_coverage(
                    t.contig,
                    t.start,
                    t.end,
                    quality_threshold=0,
                    read_callback="all",
                )
            except ValueError as exc:
                raise RuntimeError(f"pysam count_coverage failed for {t.contig}:{t.start}-{t.end}: {exc}") from exc

            depths = (
                np.array(cov[0], dtype=np.int32)
                + np.array(cov[1], dtype=np.int32)
                + np.array(cov[2], dtype=np.int32)
                + np.array(cov[3], dtype=np.int32)
            )

            length = len(depths)
            if length == 0:
                results.append(_empty_result(t, sorted_thresholds))
                continue

            mean_depth = float(depths.mean())
            stdev_depth = float(depths.std())
            min_depth = int(depths.min())
            max_depth = int(depths.max())
            median_depth = float(np.median(depths))
            zero_count = int(np.sum(depths == 0))
            pct_zero = zero_count / length * 100

            pct_thresholds = {}
            for th in sorted_thresholds:
                pct_thresholds[th] = round(float(np.sum(depths >= th)) / length * 100, 2)

            # Histogram for merge compatibility
            unique, counts = np.unique(depths, return_counts=True)
            histogram = dict(zip(unique.tolist(), counts.tolist()))

            # Low-coverage blocks
            lowcov_blocks = _detect_lowcov_blocks(
                depths,
                t.start,
                lowcov_threshold,
                lowcov_min_len,
            )
            lowcov_total_bp = sum(b.length for b in lowcov_blocks)

            results.append(
                TargetResult(
                    target_id=t.target_id,
                    contig=t.contig,
                    start=t.start,
                    end=t.end,
                    length_bp=length,
                    mean_depth=round(mean_depth, 2),
                    median_depth=round(median_depth, 1),
                    min_depth=min_depth,
                    max_depth=max_depth,
                    stdev_depth=round(stdev_depth, 2),
                    pct_zero=round(pct_zero, 2),
                    pct_thresholds=pct_thresholds,
                    n_lowcov_blocks=len(lowcov_blocks),
                    lowcov_total_bp=lowcov_total_bp,
                    lowcov_blocks=lowcov_blocks,
                    histogram=histogram,
                )
            )
    finally:
        af.close()

    return results


def _empty_result(t: _TargetSpec, thresholds: list[int]) -> TargetResult:
    """Return a zero-filled TargetResult for an empty region."""
    return TargetResult(
        target_id=t.target_id,
        contig=t.contig,
        start=t.start,
        end=t.end,
        length_bp=0,
        mean_depth=0,
        median_depth=0,
        min_depth=0,
        max_depth=0,
        stdev_depth=0,
        pct_zero=0,
        pct_thresholds={th: 0.0 for th in thresholds},
        n_lowcov_blocks=0,
        lowcov_total_bp=0,
        lowcov_blocks=[],
    )


def _detect_lowcov_blocks(
    depths,  # numpy array
    region_start: int,
    lowcov_threshold: int,
    lowcov_min_len: int,
) -> list[LowCovBlock]:
    """Detect contiguous low-coverage blocks from a numpy depth array."""
    import numpy as np

    is_low = depths < lowcov_threshold
    if not np.any(is_low):
        return []

    blocks: list[LowCovBlock] = []
    # Find transitions
    padded = np.concatenate([[False], is_low, [False]])
    diff = np.diff(padded.astype(np.int8))
    starts = np.where(diff == 1)[0]
    ends = np.where(diff == -1)[0]

    for s, e in zip(starts, ends):
        block_len = e - s
        if block_len >= lowcov_min_len:
            block_depths = depths[s:e]
            blocks.append(
                LowCovBlock(
                    start=region_start + int(s),
                    end=region_start + int(e),
                    depth_sum=int(block_depths.sum()),
                    length=block_len,
                )
            )
    return blocks


def ensure_index(bam_path: str, no_index: bool = False) -> None:
    """Ensure the BAM/CRAM file has an index.

    Creates one with ``samtools index`` if missing and *no_index* is False.
    Raises RuntimeError if the index is missing and *no_index* is True.
    """
    # Check for existing index
    for suffix in (".bai", ".csi", ".crai"):
        if os.path.exists(bam_path + suffix):
            return
    # Also check path without extension + suffix (e.g. sample.bam → sample.bai)
    base, _ = os.path.splitext(bam_path)
    for suffix in (".bai", ".csi", ".crai"):
        if os.path.exists(base + suffix):
            return

    if no_index:
        raise RuntimeError(
            f"Index not found for '{bam_path}' and --no-index prevents creation. "
            f"Create it manually: samtools index {bam_path}"
        )

    if not shutil.which("samtools"):
        raise RuntimeError(f"Index not found for '{bam_path}' and samtools is not on PATH to create one.")

    logger.info("Creating index for %s ...", bam_path)
    proc = subprocess.run(
        ["samtools", "index", bam_path],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"samtools index failed: {proc.stderr}")
    logger.info("Index created for %s", bam_path)
