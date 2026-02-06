"""Depth computation engines: samtools depth and mosdepth.

Both engines stream per-base depth data and feed it into
:class:`~covsnap.metrics.TargetAccumulator` instances,
so RAM usage is O(targets) regardless of region size.
"""

from __future__ import annotations

import gzip
import logging
import os
import shutil
import subprocess
import tempfile
from typing import Optional

from covsnap.metrics import TargetAccumulator, TargetResult

logger = logging.getLogger(__name__)


def select_engine(engine: str) -> str:
    """Resolve an engine name to a concrete engine.

    Parameters
    ----------
    engine : str
        One of ``"auto"``, ``"mosdepth"``, ``"samtools"``.

    Returns
    -------
    str
        The selected engine name (``"mosdepth"`` or ``"samtools"``).

    Raises
    ------
    RuntimeError
        If the requested engine (or any engine in auto mode) is not found.
    """
    if engine == "mosdepth":
        if shutil.which("mosdepth"):
            return "mosdepth"
        raise RuntimeError(
            "mosdepth not found on PATH. Install it or use --engine samtools."
        )

    if engine == "samtools":
        if shutil.which("samtools"):
            return "samtools"
        raise RuntimeError("samtools not found on PATH.")

    # auto mode
    if shutil.which("mosdepth"):
        logger.info("Auto-selected engine: mosdepth")
        return "mosdepth"
    if shutil.which("samtools"):
        logger.info("Auto-selected engine: samtools")
        return "samtools"

    raise RuntimeError(
        "Neither mosdepth nor samtools found on PATH. "
        "Install at least one: conda install -c bioconda samtools"
    )


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
    return [
        _TargetSpec(target_id=name, contig=contig, start=start, end=end)
        for contig, start, end, name in regions
    ]


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
        Threads for mosdepth.
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
        if engine == "samtools":
            return _run_samtools(
                bam_path, specs, thresholds, lowcov_threshold, lowcov_min_len,
                tmp_dir, reference,
            )
        elif engine == "mosdepth":
            return _run_mosdepth(
                bam_path, specs, thresholds, lowcov_threshold, lowcov_min_len,
                threads, tmp_dir, reference,
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
# samtools depth engine
# ---------------------------------------------------------------------------


def _run_samtools(
    bam_path: str,
    targets: list[_TargetSpec],
    thresholds: list[int],
    lowcov_threshold: int,
    lowcov_min_len: int,
    tmp_dir: str,
    reference: Optional[str],
) -> list[TargetResult]:
    """Stream ``samtools depth -a`` output and accumulate metrics."""
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

    # Build an index: (contig, start, end) → accumulator index
    # for fast lookup. Since targets are sorted and non-overlapping,
    # we use a cursor.
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
    """Run mosdepth and parse its outputs for full metrics."""
    bed_path = _write_tmp_bed(targets, tmp_dir)
    prefix = os.path.join(tmp_dir, "covsnap_mosdepth")

    threshold_str = ",".join(str(t) for t in sorted(thresholds))

    cmd = [
        "mosdepth",
        "--by", bed_path,
        "--thresholds", threshold_str,
        "--threads", str(threads),
        "--no-per-base",  # We'll do a separate per-base pass
        prefix,
        bam_path,
    ]
    if reference:
        cmd.extend(["--fasta", reference])

    logger.info("Running mosdepth (thresholds): %s", " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"mosdepth failed (exit {proc.returncode}): {proc.stderr}")

    # Parse thresholds output for pct_ge_X
    thresholds_file = f"{prefix}.thresholds.bed.gz"
    threshold_data = _parse_mosdepth_thresholds(thresholds_file, targets, thresholds)

    # Parse regions output for mean depth
    regions_file = f"{prefix}.regions.bed.gz"
    mean_depths = _parse_mosdepth_regions(regions_file, targets)

    # Now run mosdepth again WITH per-base to get min/max/median/stdev/lowcov
    # (removing --no-per-base)
    cmd_perbase = [
        "mosdepth",
        "--by", bed_path,
        "--threads", str(threads),
        prefix + "_pb",
        bam_path,
    ]
    if reference:
        cmd_perbase.extend(["--fasta", reference])

    logger.info("Running mosdepth (per-base): %s", " ".join(cmd_perbase))
    proc2 = subprocess.run(cmd_perbase, capture_output=True, text=True)
    if proc2.returncode != 0:
        raise RuntimeError(f"mosdepth per-base failed (exit {proc2.returncode}): {proc2.stderr}")

    perbase_file = f"{prefix}_pb.per-base.bed.gz"
    perbase_results = _stream_mosdepth_perbase(
        perbase_file, targets, thresholds, lowcov_threshold, lowcov_min_len,
    )

    # Merge: use thresholds from mosdepth output (more efficient),
    # but min/max/median/stdev/lowcov from per-base streaming
    results: list[TargetResult] = []
    for i, t in enumerate(targets):
        pb = perbase_results[i]
        # Override pct_thresholds with mosdepth's exact counts
        if i < len(threshold_data):
            pb.pct_thresholds = threshold_data[i]
            pb.pct_zero = round(
                (1.0 - threshold_data[i].get(1, 0) / 100.0) * 100, 2
            ) if 1 in threshold_data[i] else pb.pct_zero
        # Override mean_depth with mosdepth's value
        if i < len(mean_depths):
            pb.mean_depth = mean_depths[i]
        results.append(pb)

    return results


def _parse_mosdepth_thresholds(
    path: str,
    targets: list[_TargetSpec],
    thresholds: list[int],
) -> list[dict[int, float]]:
    """Parse mosdepth thresholds.bed.gz for per-target threshold percentages."""
    results: list[dict[int, float]] = []
    if not os.path.exists(path):
        logger.warning("mosdepth thresholds file not found: %s", path)
        return results

    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            # parts: chrom, start, end, [name,] T1_count, T2_count, ...
            # mosdepth thresholds output: each value is the number of bases >= threshold
            contig = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            length = end - start

            # Values start at index 4 if there's a name column, else index 3
            # Detect by checking if we have enough numeric columns
            val_start = 3
            if len(parts) > 3 and not parts[3].replace(".", "").isdigit():
                val_start = 4

            pct_map: dict[int, float] = {}
            for j, t in enumerate(sorted(thresholds)):
                idx = val_start + j
                if idx < len(parts):
                    try:
                        count = float(parts[idx])
                        pct_map[t] = round(count / length * 100, 2) if length > 0 else 0.0
                    except ValueError:
                        pct_map[t] = 0.0
                else:
                    pct_map[t] = 0.0
            results.append(pct_map)

    return results


def _parse_mosdepth_regions(
    path: str,
    targets: list[_TargetSpec],
) -> list[float]:
    """Parse mosdepth regions.bed.gz for mean depths."""
    means: list[float] = []
    if not os.path.exists(path):
        logger.warning("mosdepth regions file not found: %s", path)
        return means

    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            # Last column is mean depth
            try:
                means.append(round(float(parts[-1]), 2))
            except (ValueError, IndexError):
                means.append(0.0)

    return means


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
            block_count = block_end - block_start

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
        raise RuntimeError(
            f"Index not found for '{bam_path}' and samtools is not on PATH to create one."
        )

    logger.info("Creating index for %s ...", bam_path)
    proc = subprocess.run(
        ["samtools", "index", bam_path],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"samtools index failed: {proc.stderr}")
    logger.info("Index created for %s", bam_path)
