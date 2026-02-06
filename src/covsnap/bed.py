"""BED file streaming parser with guardrails for large files.

All BED processing uses streaming iteration — no file is ever loaded
entirely into memory.  Guardrails enforce limits on target count,
total base-pairs, and file size, with configurable responses
(error, clip, or sample).
"""

from __future__ import annotations

import logging
import os
import random
import sys
from dataclasses import dataclass, field
from typing import Iterator, Optional

logger = logging.getLogger(__name__)


@dataclass
class BedInterval:
    """A single BED interval."""

    contig: str
    start: int  # 0-based
    end: int  # half-open
    name: str = ""

    @property
    def length_bp(self) -> int:
        return self.end - self.start


@dataclass
class BedLimitStats:
    """Statistics from a pre-scan of a BED file."""

    n_targets: int = 0
    total_bp: int = 0
    file_bytes: int = 0
    exceeds_targets: bool = False
    exceeds_bp: bool = False
    exceeds_bytes: bool = False

    @property
    def any_exceeded(self) -> bool:
        return self.exceeds_targets or self.exceeds_bp or self.exceeds_bytes


@dataclass
class BedEnforceResult:
    """Result of enforcing BED limits."""

    intervals: list[BedInterval] = field(default_factory=list)
    clipped: bool = False
    sampled: bool = False
    original_stats: Optional[BedLimitStats] = None
    kept_count: int = 0
    kept_bp: int = 0
    message: str = ""


def stream_bed_intervals(bed_path: str) -> Iterator[BedInterval]:
    """Yield BED intervals one at a time from a file.

    Skips comment lines (starting with ``#``) and blank lines.
    Accepts BED3+ format (contig, start, end, optional name, ...).
    """
    with open(bed_path) as fh:
        for line_no, line in enumerate(fh, 1):
            line = line.rstrip("\n\r")
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                logger.warning("BED line %d: fewer than 3 columns, skipping: %s", line_no, line)
                continue
            try:
                contig = parts[0]
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                logger.warning("BED line %d: non-integer coordinates, skipping: %s", line_no, line)
                continue
            if start >= end:
                logger.warning("BED line %d: start >= end (%d >= %d), skipping", line_no, start, end)
                continue
            name = parts[3] if len(parts) >= 4 else f"{contig}:{start}-{end}"
            yield BedInterval(contig=contig, start=start, end=end, name=name)


def check_bed_limits(
    bed_path: str,
    max_targets: int = 2000,
    max_total_bp: int = 50_000_000,
    max_bed_bytes: int = 50 * 1024 * 1024,
) -> BedLimitStats:
    """Single-pass scan of a BED file to check limits.

    Does NOT store intervals in memory — only accumulates counts.
    """
    stats = BedLimitStats()
    stats.file_bytes = os.path.getsize(bed_path)
    stats.exceeds_bytes = stats.file_bytes > max_bed_bytes

    for iv in stream_bed_intervals(bed_path):
        stats.n_targets += 1
        stats.total_bp += iv.length_bp

    stats.exceeds_targets = stats.n_targets > max_targets
    stats.exceeds_bp = stats.total_bp > max_total_bp
    return stats


def enforce_limits(
    bed_path: str,
    max_targets: int = 2000,
    max_total_bp: int = 50_000_000,
    max_bed_bytes: int = 50 * 1024 * 1024,
    on_large_bed: str = "warn_and_clip",
    seed: int = 42,
) -> BedEnforceResult:
    """Check BED limits and enforce the configured policy.

    Parameters
    ----------
    bed_path : str
        Path to the BED file.
    max_targets, max_total_bp, max_bed_bytes :
        Configurable limits.
    on_large_bed : str
        One of ``"error"``, ``"warn_and_clip"``, ``"warn_and_sample"``.
    seed : int
        Random seed for ``warn_and_sample`` mode.

    Returns
    -------
    BedEnforceResult
        Contains the (possibly reduced) interval list and metadata
        about what action was taken.

    Raises
    ------
    SystemExit
        If ``on_large_bed == "error"`` and any limit is exceeded.
    """
    stats = check_bed_limits(bed_path, max_targets, max_total_bp, max_bed_bytes)

    if not stats.any_exceeded:
        # All within limits — return all intervals
        intervals = list(stream_bed_intervals(bed_path))
        return BedEnforceResult(
            intervals=intervals,
            clipped=False,
            sampled=False,
            original_stats=stats,
            kept_count=len(intervals),
            kept_bp=sum(iv.length_bp for iv in intervals),
            message="All BED limits satisfied.",
        )

    # Limits exceeded
    limit_msgs: list[str] = []
    if stats.exceeds_targets:
        limit_msgs.append(f"target count ({stats.n_targets}) exceeds --max-targets ({max_targets})")
    if stats.exceeds_bp:
        limit_msgs.append(f"total bp ({stats.total_bp:,}) exceeds --max-total-bp ({max_total_bp:,})")
    if stats.exceeds_bytes:
        limit_msgs.append(
            f"file size ({stats.file_bytes:,} bytes) exceeds --max-bed-bytes ({max_bed_bytes:,})"
        )
    detail = "; ".join(limit_msgs)

    if on_large_bed == "error":
        logger.error("BED guardrail: %s", detail)
        sys.exit(4)

    if on_large_bed == "warn_and_clip":
        return _clip_bed(bed_path, max_targets, max_total_bp, stats, detail)

    if on_large_bed == "warn_and_sample":
        return _sample_bed(bed_path, max_targets, max_total_bp, stats, detail, seed)

    raise ValueError(f"Unknown --on-large-bed value: {on_large_bed}")


def _clip_bed(
    bed_path: str,
    max_targets: int,
    max_total_bp: int,
    stats: BedLimitStats,
    detail: str,
) -> BedEnforceResult:
    """Keep the first N targets until both limits are satisfied."""
    intervals: list[BedInterval] = []
    running_bp = 0

    for iv in stream_bed_intervals(bed_path):
        if len(intervals) >= max_targets:
            break
        if running_bp + iv.length_bp > max_total_bp:
            break
        intervals.append(iv)
        running_bp += iv.length_bp

    msg = (
        f"BED guardrails triggered (warn_and_clip): {detail}. "
        f"Retained first {len(intervals)} of {stats.n_targets} targets "
        f"({running_bp:,} bp of {stats.total_bp:,} bp)."
    )
    logger.warning("%s", msg)

    return BedEnforceResult(
        intervals=intervals,
        clipped=True,
        sampled=False,
        original_stats=stats,
        kept_count=len(intervals),
        kept_bp=running_bp,
        message=msg,
    )


def _sample_bed(
    bed_path: str,
    max_targets: int,
    max_total_bp: int,
    stats: BedLimitStats,
    detail: str,
    seed: int,
) -> BedEnforceResult:
    """Deterministic reservoir sampling of targets.

    Uses reservoir sampling (Algorithm R) with the given seed so that
    the same seed always selects the same targets.  After sampling,
    intervals are sorted by genomic position and trimmed if total bp
    still exceeds the limit.
    """
    rng = random.Random(seed)
    reservoir: list[BedInterval] = []

    for i, iv in enumerate(stream_bed_intervals(bed_path)):
        if i < max_targets:
            reservoir.append(iv)
        else:
            j = rng.randint(0, i)
            if j < max_targets:
                reservoir[j] = iv

    # Sort by genomic position for consistent output
    reservoir.sort(key=lambda iv: (iv.contig, iv.start))

    # Trim by total bp if needed
    final: list[BedInterval] = []
    running_bp = 0
    for iv in reservoir:
        if running_bp + iv.length_bp > max_total_bp:
            break
        final.append(iv)
        running_bp += iv.length_bp

    msg = (
        f"BED guardrails triggered (warn_and_sample): {detail}. "
        f"Sampled {len(final)} of {stats.n_targets} targets "
        f"({running_bp:,} bp of {stats.total_bp:,} bp) with seed={seed}."
    )
    logger.warning("%s", msg)

    return BedEnforceResult(
        intervals=final,
        sampled=True,
        clipped=False,
        original_stats=stats,
        kept_count=len(final),
        kept_bp=running_bp,
        message=msg,
    )
