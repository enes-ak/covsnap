"""Tests for the streaming metrics accumulator."""

import math

from covsnap.metrics import TargetAccumulator


class TestTargetAccumulator:
    """Test the TargetAccumulator with known depth patterns."""

    def _make_acc(self, start=0, end=100, thresholds=None):
        if thresholds is None:
            thresholds = [1, 5, 10, 20, 30, 50, 100]
        return TargetAccumulator(
            target_id="test",
            contig="chr1",
            start=start,
            end=end,
            thresholds=thresholds,
            lowcov_threshold=10,
            lowcov_min_len=5,
        )

    def test_uniform_depth(self):
        """All positions have the same depth."""
        acc = self._make_acc(start=0, end=10)
        for pos in range(10):
            acc.add_single(50, pos)
        result = acc.finalize()

        assert result.length_bp == 10
        assert result.mean_depth == 50.0
        assert result.median_depth == 50.0
        assert result.min_depth == 50
        assert result.max_depth == 50
        assert result.stdev_depth == 0.0
        assert result.pct_zero == 0.0
        assert result.pct_thresholds[1] == 100.0
        assert result.pct_thresholds[50] == 100.0
        assert result.pct_thresholds[100] == 0.0

    def test_all_zeros(self):
        """All positions have zero depth."""
        acc = self._make_acc(start=0, end=10)
        for pos in range(10):
            acc.add_single(0, pos)
        result = acc.finalize()

        assert result.mean_depth == 0.0
        assert result.median_depth == 0.0
        assert result.min_depth == 0
        assert result.max_depth == 0
        assert result.pct_zero == 100.0
        assert result.pct_thresholds[1] == 0.0

    def test_mixed_depths(self):
        """Mix of zero and non-zero depths."""
        acc = self._make_acc(start=0, end=10)
        depths = [0, 0, 5, 10, 20, 30, 50, 100, 100, 100]
        for pos, d in enumerate(depths):
            acc.add_single(d, pos)
        result = acc.finalize()

        assert result.length_bp == 10
        assert result.min_depth == 0
        assert result.max_depth == 100
        assert result.pct_zero == 20.0  # 2 out of 10
        expected_mean = sum(depths) / len(depths)
        assert abs(result.mean_depth - expected_mean) < 0.01

    def test_gap_handling(self):
        """Positions with gaps are filled with zero depth."""
        acc = self._make_acc(start=0, end=10)
        # Only report positions 0, 1, 8, 9 — gap from 2-7
        acc.add_single(50, 0)
        acc.add_single(50, 1)
        acc.add_single(50, 8)
        acc.add_single(50, 9)
        result = acc.finalize()

        assert result.length_bp == 10
        # 4 positions at 50, 6 at 0
        assert result.pct_zero == 60.0
        expected_mean = (50 * 4) / 10
        assert abs(result.mean_depth - expected_mean) < 0.01

    def test_block_mode(self):
        """add_block accumulates correctly for run-length-encoded input."""
        acc = self._make_acc(start=0, end=100)
        acc.add_block(depth=30, count=50, block_start=0)
        acc.add_block(depth=60, count=50, block_start=50)
        result = acc.finalize()

        assert result.length_bp == 100
        assert result.min_depth == 30
        assert result.max_depth == 60
        expected_mean = (30 * 50 + 60 * 50) / 100
        assert abs(result.mean_depth - expected_mean) < 0.01
        # Median: 50 values of 30, 50 values of 60 → median between 30 and 60
        assert result.median_depth == 30.0 or result.median_depth == 45.0

    def test_lowcov_block_detection(self):
        """Low-coverage blocks are detected correctly."""
        acc = self._make_acc(start=0, end=20, thresholds=[1, 10, 20])
        # First 10 positions: depth 0 (below threshold 10, block length 10 >= min 5)
        for pos in range(10):
            acc.add_single(0, pos)
        # Next 10 positions: depth 50 (above threshold)
        for pos in range(10, 20):
            acc.add_single(50, pos)
        result = acc.finalize()

        assert result.n_lowcov_blocks == 1
        assert result.lowcov_blocks[0].start == 0
        assert result.lowcov_blocks[0].end == 10
        assert result.lowcov_blocks[0].length == 10
        assert result.lowcov_total_bp == 10

    def test_lowcov_block_too_short(self):
        """Low-coverage blocks shorter than min_len are not reported."""
        acc = self._make_acc(start=0, end=10)
        # 3 positions at 0 (below min_len=5), then 7 at 50
        for pos in range(3):
            acc.add_single(0, pos)
        for pos in range(3, 10):
            acc.add_single(50, pos)
        result = acc.finalize()

        assert result.n_lowcov_blocks == 0

    def test_welford_stdev(self):
        """Verify Welford's algorithm produces correct stdev."""
        acc = self._make_acc(start=0, end=4)
        depths = [10, 20, 30, 40]
        for pos, d in enumerate(depths):
            acc.add_single(d, pos)
        result = acc.finalize()

        # Population stdev of [10, 20, 30, 40]
        mean = 25.0
        variance = sum((d - mean) ** 2 for d in depths) / len(depths)
        expected_stdev = math.sqrt(variance)
        assert abs(result.stdev_depth - expected_stdev) < 0.01

    def test_trailing_gap_filled(self):
        """Positions beyond last reported pos up to target end are filled with 0."""
        acc = self._make_acc(start=0, end=10)
        # Only report first 5 positions
        for pos in range(5):
            acc.add_single(20, pos)
        result = acc.finalize()

        # 5 at 20, 5 at 0
        assert result.pct_zero == 50.0
        expected_mean = (20 * 5) / 10
        assert abs(result.mean_depth - expected_mean) < 0.01

    def test_no_data_gives_zeros(self):
        """Target with no depth data produces all-zero metrics."""
        acc = self._make_acc(start=0, end=10)
        result = acc.finalize()

        assert result.mean_depth == 0.0
        assert result.median_depth == 0.0
        assert result.min_depth == 0
        assert result.max_depth == 0
        assert result.pct_zero == 100.0

    def test_block_clipping(self):
        """Blocks partially outside target are clipped to target boundaries."""
        acc = self._make_acc(start=10, end=20)
        # Block from 5 to 25 — should be clipped to 10-20
        acc.add_block(depth=30, count=20, block_start=5)
        result = acc.finalize()

        assert result.length_bp == 10
        # All 10 positions within target should have depth 30
        assert result.mean_depth == 30.0
        assert result.pct_thresholds[20] == 100.0
