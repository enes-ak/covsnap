"""Tests for depth computation engines."""

import shutil

import pytest

from covsnap.engines import compute_depth, select_engine


class TestEngineSelection:
    """Test engine auto-selection logic."""

    def test_auto_selects_something(self):
        engine = select_engine("auto")
        assert engine in ("mosdepth", "samtools")

    def test_explicit_samtools(self):
        if not shutil.which("samtools"):
            pytest.skip("samtools not installed")
        engine = select_engine("samtools")
        assert engine == "samtools"

    def test_missing_engine_raises(self):
        if shutil.which("mosdepth"):
            pytest.skip("mosdepth is installed, cannot test missing-engine path")
        with pytest.raises(RuntimeError, match="not found"):
            select_engine("mosdepth")


class TestSamtoolsEngine:
    """Test the samtools depth streaming engine."""

    @pytest.fixture(autouse=True)
    def _require_samtools(self):
        if not shutil.which("samtools"):
            pytest.skip("samtools not installed")

    def test_single_region(self, synthetic_bam):
        results = compute_depth(
            bam_path=synthetic_bam,
            regions=[("chr17", 43044294, 43050000, "test_region")],
            engine="samtools",
            thresholds=[1, 5, 10, 20, 30, 50, 100],
            lowcov_threshold=10,
            lowcov_min_len=50,
        )
        assert len(results) == 1
        r = results[0]
        assert r.contig == "chr17"
        assert r.start == 43044294
        assert r.end == 43050000
        assert r.length_bp == 5706
        assert isinstance(r.mean_depth, float)
        assert r.mean_depth >= 0
        assert r.min_depth >= 0
        assert r.max_depth >= r.min_depth
        assert 0 <= r.pct_zero <= 100
        assert 0 <= r.pct_thresholds.get(20, 0) <= 100

    def test_result_has_no_raw_array(self, synthetic_bam):
        """Engine returns aggregated stats, not raw per-base arrays."""
        results = compute_depth(
            bam_path=synthetic_bam,
            regions=[("chr17", 43044294, 43125482, "BRCA1")],
            engine="samtools",
            thresholds=[1, 5, 10, 20, 30, 50, 100],
        )
        r = results[0]
        assert isinstance(r.mean_depth, (int, float))
        assert isinstance(r.median_depth, (int, float))
        assert not hasattr(r, "per_base_depths")

    def test_multiple_regions(self, synthetic_bam):
        regions = [
            ("chr17", 43044294, 43050000, "region_1"),
            ("chr17", 43090000, 43095000, "region_2"),
        ]
        results = compute_depth(
            bam_path=synthetic_bam,
            regions=regions,
            engine="samtools",
            thresholds=[1, 10, 20],
        )
        assert len(results) == 2
        assert results[0].target_id == "region_1"
        assert results[1].target_id == "region_2"

    def test_region_outside_reads(self, synthetic_bam):
        """Region with no reads produces all-zero metrics."""
        results = compute_depth(
            bam_path=synthetic_bam,
            regions=[("chr17", 1000000, 1001000, "empty_region")],
            engine="samtools",
            thresholds=[1, 10, 20],
        )
        assert len(results) == 1
        r = results[0]
        assert r.mean_depth == 0.0
        assert r.pct_zero == 100.0


class TestMosdepthEngine:
    """Test mosdepth engine (skipped if mosdepth not installed)."""

    @pytest.fixture(autouse=True)
    def _require_mosdepth(self):
        if not shutil.which("mosdepth"):
            pytest.skip("mosdepth not installed")

    def test_single_region(self, synthetic_bam, tmp_output_dir):
        results = compute_depth(
            bam_path=synthetic_bam,
            regions=[("chr17", 43044294, 43050000, "test_region")],
            engine="mosdepth",
            thresholds=[1, 5, 10, 20, 30, 50, 100],
            threads=2,
            tmp_dir=str(tmp_output_dir),
        )
        assert len(results) == 1
        r = results[0]
        assert r.mean_depth >= 0
        assert r.length_bp == 5706
