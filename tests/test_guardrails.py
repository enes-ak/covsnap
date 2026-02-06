"""Tests for BED file guardrails (streaming, limits, clipping, sampling)."""

import pytest

from covsnap.bed import (
    BedInterval,
    check_bed_limits,
    enforce_limits,
    stream_bed_intervals,
)


class TestBedStreaming:
    """Verify BED files are processed in streaming fashion."""

    def test_bed_stats_no_storage(self, large_bed):
        """check_bed_limits scans without storing intervals."""
        stats = check_bed_limits(large_bed)
        assert stats.n_targets == 3000
        assert stats.total_bp == 3000 * 500
        assert stats.file_bytes > 0

    def test_streaming_iteration(self, sample_bed):
        """stream_bed_intervals yields intervals one at a time."""
        intervals = list(stream_bed_intervals(sample_bed))
        assert len(intervals) == 2
        assert intervals[0].contig == "chr17"
        assert intervals[0].start == 43044294
        assert intervals[0].end == 43050000
        assert intervals[0].name == "region_1"

    def test_skips_comments_and_blanks(self, tmp_path):
        bed = tmp_path / "comments.bed"
        bed.write_text(
            "# comment line\n"
            "\n"
            "track name=test\n"
            "chr1\t100\t200\tgood_interval\n"
            "browser position chr1:1-1000\n"
            "chr1\t300\t400\tanother_good\n"
        )
        intervals = list(stream_bed_intervals(str(bed)))
        assert len(intervals) == 2

    def test_skips_invalid_lines(self, tmp_path):
        bed = tmp_path / "bad.bed"
        bed.write_text(
            "chr1\t100\t200\tok\n"
            "chr1\tnot_a_number\t200\n"
            "chr1\t500\t300\n"  # start >= end
            "chr1\t600\t700\tok2\n"
        )
        intervals = list(stream_bed_intervals(str(bed)))
        assert len(intervals) == 2  # only the valid ones


class TestBedLimits:
    """Test limit checking and enforcement."""

    def test_within_limits(self, sample_bed):
        result = enforce_limits(
            sample_bed,
            max_targets=2000,
            max_total_bp=50_000_000,
            max_bed_bytes=50 * 1024 * 1024,
            on_large_bed="error",
        )
        assert result.clipped is False
        assert result.sampled is False
        assert len(result.intervals) == 2

    def test_exceeds_max_targets_error(self, large_bed):
        with pytest.raises(SystemExit) as exc_info:
            enforce_limits(
                large_bed,
                max_targets=2000,
                max_total_bp=50_000_000,
                max_bed_bytes=50 * 1024 * 1024,
                on_large_bed="error",
            )
        assert exc_info.value.code == 4

    def test_clip_mode_respects_limit(self, large_bed):
        result = enforce_limits(
            large_bed,
            max_targets=100,
            max_total_bp=50_000_000,
            max_bed_bytes=50 * 1024 * 1024,
            on_large_bed="warn_and_clip",
        )
        assert result.clipped is True
        assert len(result.intervals) == 100
        assert result.intervals[0].name == "target_0000"
        assert result.intervals[99].name == "target_0099"

    def test_sample_mode_deterministic(self, large_bed):
        result1 = enforce_limits(
            large_bed,
            max_targets=50,
            max_total_bp=50_000_000,
            max_bed_bytes=50 * 1024 * 1024,
            on_large_bed="warn_and_sample",
            seed=42,
        )
        result2 = enforce_limits(
            large_bed,
            max_targets=50,
            max_total_bp=50_000_000,
            max_bed_bytes=50 * 1024 * 1024,
            on_large_bed="warn_and_sample",
            seed=42,
        )
        names1 = [iv.name for iv in result1.intervals]
        names2 = [iv.name for iv in result2.intervals]
        assert names1 == names2
        assert len(result1.intervals) == 50

    def test_sample_mode_different_seeds(self, large_bed):
        result1 = enforce_limits(
            large_bed,
            max_targets=50,
            max_total_bp=50_000_000,
            max_bed_bytes=50 * 1024 * 1024,
            on_large_bed="warn_and_sample",
            seed=42,
        )
        result2 = enforce_limits(
            large_bed,
            max_targets=50,
            max_total_bp=50_000_000,
            max_bed_bytes=50 * 1024 * 1024,
            on_large_bed="warn_and_sample",
            seed=99,
        )
        names1 = [iv.name for iv in result1.intervals]
        names2 = [iv.name for iv in result2.intervals]
        assert names1 != names2


class TestBedTotalBpLimit:
    """Test the max_total_bp limit independently."""

    def test_clip_by_total_bp(self, tmp_path):
        bed = tmp_path / "bp_test.bed"
        lines = []
        for i in range(10):
            start = i * 20000
            end = start + 10000
            lines.append(f"chr1\t{start}\t{end}\tt_{i}")
        bed.write_text("\n".join(lines) + "\n")

        result = enforce_limits(
            str(bed),
            max_targets=1000,
            max_total_bp=50000,
            max_bed_bytes=50 * 1024 * 1024,
            on_large_bed="warn_and_clip",
        )
        assert result.clipped is True
        assert len(result.intervals) == 5
        total_bp = sum(iv.length_bp for iv in result.intervals)
        assert total_bp == 50000
