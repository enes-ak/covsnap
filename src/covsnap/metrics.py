"""Streaming metrics accumulator for per-base depth data.

Computes mean, median, min, max, stdev, threshold percentages, and
low-coverage blocks without storing the full per-base depth array.

Uses:
- Welford's online algorithm for mean and variance (batch-capable).
- An exact-count histogram (depth → count) for median estimation.
- Streaming block detection for low-coverage regions.
"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any


@dataclass
class LowCovBlock:
    """A contiguous block of low-coverage positions."""

    start: int  # 0-based
    end: int  # half-open
    depth_sum: int = 0
    length: int = 0

    @property
    def mean_depth(self) -> float:
        return self.depth_sum / self.length if self.length > 0 else 0.0


@dataclass
class TargetResult:
    """Computed coverage metrics for a single target region."""

    target_id: str
    contig: str
    start: int  # 0-based
    end: int  # half-open
    length_bp: int
    mean_depth: float
    median_depth: float
    min_depth: int
    max_depth: int
    stdev_depth: float
    pct_zero: float
    pct_thresholds: dict[int, float]  # threshold → percentage
    n_lowcov_blocks: int
    lowcov_total_bp: int
    lowcov_blocks: list[LowCovBlock]
    engine_used: str = ""
    bam_path: str = ""
    sample_name: str = "NA"
    coverage_status: str = ""  # filled by classifier


class TargetAccumulator:
    """Streaming accumulator for depth values over a target region.

    Supports two modes of input:

    - ``add_single(depth, pos)`` — for samtools depth (one base at a time).
    - ``add_block(depth, count, block_start)`` — for mosdepth per-base BED
      (run-length encoded blocks of identical depth).

    Call ``finalize()`` to produce a ``TargetResult``.
    """

    def __init__(
        self,
        target_id: str,
        contig: str,
        start: int,
        end: int,
        thresholds: list[int],
        lowcov_threshold: int = 10,
        lowcov_min_len: int = 50,
    ) -> None:
        self.target_id = target_id
        self.contig = contig
        self.start = start
        self.end = end
        self.expected_length = end - start
        self.thresholds = sorted(thresholds)
        self.lowcov_threshold = lowcov_threshold
        self.lowcov_min_len = lowcov_min_len

        # Welford's accumulators
        self._n: int = 0
        self._sum: int = 0
        self._min: int = 2**31
        self._max: int = 0
        self._welford_mean: float = 0.0
        self._welford_m2: float = 0.0

        # Histogram: depth → base count (exact for median computation)
        self._histogram: dict[int, int] = defaultdict(int)

        # Threshold counts: how many bases >= each threshold
        self._threshold_counts: dict[int, int] = {t: 0 for t in self.thresholds}
        self._zero_count: int = 0

        # Low-coverage block tracking
        self._in_lowcov: bool = False
        self._lowcov_block_start: int = 0
        self._lowcov_block_depth_sum: int = 0
        self._lowcov_block_len: int = 0
        self._lowcov_blocks: list[LowCovBlock] = []

        # Track the last position we've seen (for gap detection)
        self._last_pos: int = start - 1

    def add_single(self, depth: int, pos: int) -> None:
        """Add a single base position with the given depth.

        *pos* is 0-based.  Positions must arrive in order.
        """
        # Handle any gap between last_pos and this pos (missing = 0 depth)
        gap = pos - self._last_pos - 1
        if gap > 0:
            self._add_bulk(0, gap, self._last_pos + 1)

        self._add_bulk(depth, 1, pos)
        self._last_pos = pos

    def add_block(self, depth: int, count: int, block_start: int) -> None:
        """Add *count* consecutive positions with the same *depth*.

        Used for mosdepth's run-length-encoded per-base BED output.
        *block_start* is the 0-based start position of this block.
        """
        # Handle gap before this block
        gap = block_start - self._last_pos - 1
        if gap > 0:
            self._add_bulk(0, gap, self._last_pos + 1)

        # Clip block to target boundaries
        block_end = block_start + count
        clipped_start = max(block_start, self.start)
        clipped_end = min(block_end, self.end)
        if clipped_start >= clipped_end:
            self._last_pos = max(self._last_pos, block_end - 1)
            return
        clipped_count = clipped_end - clipped_start

        self._add_bulk(depth, clipped_count, clipped_start)
        self._last_pos = clipped_end - 1

    def _add_bulk(self, depth: int, count: int, block_start: int) -> None:
        """Core accumulation: add *count* bases with identical *depth*."""
        if count <= 0:
            return

        # --- Welford's batch update for identical values ---
        # Given existing (n, mean, M2) and k new values all equal to x:
        #   n_new = n + k
        #   delta = x - mean
        #   mean_new = mean + delta * k / n_new
        #   M2_new = M2 + delta^2 * n * k / n_new
        old_n = self._n
        new_n = old_n + count
        delta = depth - self._welford_mean
        self._welford_mean += delta * count / new_n
        self._welford_m2 += delta * delta * old_n * count / new_n

        self._n = new_n
        self._sum += depth * count

        # Min/Max
        if depth < self._min:
            self._min = depth
        if depth > self._max:
            self._max = depth

        # Histogram
        self._histogram[depth] += count

        # Threshold counts
        for t in self.thresholds:
            if depth >= t:
                self._threshold_counts[t] += count

        if depth == 0:
            self._zero_count += count

        # Low-coverage block tracking
        block_end = block_start + count
        if depth < self.lowcov_threshold:
            if not self._in_lowcov:
                # Start new low-cov block
                self._in_lowcov = True
                self._lowcov_block_start = block_start
                self._lowcov_block_depth_sum = depth * count
                self._lowcov_block_len = count
            else:
                # Extend existing block
                self._lowcov_block_depth_sum += depth * count
                self._lowcov_block_len += count
        else:
            if self._in_lowcov:
                self._close_lowcov_block()

    def _close_lowcov_block(self) -> None:
        """Close the current low-coverage block if it meets min length."""
        if self._in_lowcov and self._lowcov_block_len >= self.lowcov_min_len:
            self._lowcov_blocks.append(
                LowCovBlock(
                    start=self._lowcov_block_start,
                    end=self._lowcov_block_start + self._lowcov_block_len,
                    depth_sum=self._lowcov_block_depth_sum,
                    length=self._lowcov_block_len,
                )
            )
        self._in_lowcov = False
        self._lowcov_block_len = 0
        self._lowcov_block_depth_sum = 0

    def finalize(self) -> TargetResult:
        """Finalise accumulation and return computed metrics.

        Fills any trailing gap (positions not reported up to target end)
        with zero depth, then computes all derived metrics.
        """
        # Fill trailing gap to target end
        trailing = self.end - self._last_pos - 1
        if trailing > 0:
            self._add_bulk(0, trailing, self._last_pos + 1)

        # Close any open low-cov block
        self._close_lowcov_block()

        n = self._n
        length = self.expected_length

        # If we never saw any positions, fill with zeros
        if n == 0:
            n = length
            self._zero_count = length
            self._histogram[0] = length
            self._min = 0
            self._max = 0

        # Mean
        mean_depth = self._sum / n if n > 0 else 0.0

        # Stdev (population stdev)
        variance = self._welford_m2 / n if n > 1 else 0.0
        stdev = math.sqrt(variance)

        # Median from histogram
        median_depth = self._compute_median()

        # Percentages
        pct_zero = (self._zero_count / n * 100) if n > 0 else 0.0
        pct_thresholds = {}
        for t in self.thresholds:
            pct_thresholds[t] = (self._threshold_counts[t] / n * 100) if n > 0 else 0.0

        # Low-cov totals
        lowcov_total_bp = sum(b.length for b in self._lowcov_blocks)

        return TargetResult(
            target_id=self.target_id,
            contig=self.contig,
            start=self.start,
            end=self.end,
            length_bp=length,
            mean_depth=round(mean_depth, 2),
            median_depth=round(median_depth, 1),
            min_depth=self._min if self._min != 2**31 else 0,
            max_depth=self._max,
            stdev_depth=round(stdev, 2),
            pct_zero=round(pct_zero, 2),
            pct_thresholds={t: round(v, 2) for t, v in pct_thresholds.items()},
            n_lowcov_blocks=len(self._lowcov_blocks),
            lowcov_total_bp=lowcov_total_bp,
            lowcov_blocks=self._lowcov_blocks,
        )

    def _compute_median(self) -> float:
        """Compute median depth from the exact histogram."""
        if not self._histogram:
            return 0.0
        total = sum(self._histogram.values())
        if total == 0:
            return 0.0

        mid = total / 2.0
        cumulative = 0
        prev_depth = 0

        for depth in sorted(self._histogram.keys()):
            cumulative += self._histogram[depth]
            if cumulative >= mid:
                if total % 2 == 1:
                    return float(depth)
                # Even count: average of the two middle values
                if cumulative - self._histogram[depth] < mid:
                    # Both middle values are in this bin
                    return float(depth)
                # One middle value is in the previous bin
                return (prev_depth + depth) / 2.0
            prev_depth = depth

        return 0.0


def merge_accumulator_results(
    results: list[TargetResult],
    extra_fields: dict[str, Any],
) -> list[TargetResult]:
    """Attach shared metadata fields (engine_used, bam_path, etc.) to results."""
    for r in results:
        for k, v in extra_fields.items():
            if hasattr(r, k):
                setattr(r, k, v)
    return results
