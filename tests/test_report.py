"""Tests for classification logic and report generation."""

import os

from covsnap.metrics import LowCovBlock, TargetResult
from covsnap.report import ClassifyParams, classify_target


def _make_result(**overrides) -> TargetResult:
    """Create a TargetResult with sensible defaults, overriding as needed."""
    defaults = dict(
        target_id="TEST",
        contig="chr1",
        start=0,
        end=1000,
        length_bp=1000,
        mean_depth=100.0,
        median_depth=95.0,
        min_depth=5,
        max_depth=300,
        stdev_depth=30.0,
        pct_zero=0.5,
        pct_thresholds={1: 99.5, 5: 99.0, 10: 98.5, 20: 97.0, 30: 95.0, 50: 80.0, 100: 40.0},
        n_lowcov_blocks=0,
        lowcov_total_bp=0,
        lowcov_blocks=[],
    )
    defaults.update(overrides)
    return TargetResult(**defaults)


class TestClassification:
    """Test coverage classification heuristics."""

    def test_pass(self):
        r = _make_result(pct_zero=0.5, pct_thresholds={20: 97.0})
        params = ClassifyParams()
        status = classify_target(r, params)
        assert status == "PASS"
        assert r.coverage_status == "PASS"

    def test_low_coverage(self):
        r = _make_result(pct_zero=0.5, pct_thresholds={20: 80.0})
        params = ClassifyParams()
        status = classify_target(r, params)
        assert status == "LOW_COVERAGE"

    def test_dropout_by_pct_zero(self):
        r = _make_result(pct_zero=10.0, pct_thresholds={20: 70.0})
        params = ClassifyParams()
        status = classify_target(r, params)
        assert status == "DROP_OUT"

    def test_dropout_by_zero_block(self):
        block = LowCovBlock(start=100, end=700, depth_sum=0, length=600)
        r = _make_result(
            pct_zero=2.0,
            pct_thresholds={20: 90.0},
            lowcov_blocks=[block],
        )
        params = ClassifyParams()
        status = classify_target(r, params)
        assert status == "DROP_OUT"

    def test_uneven(self):
        # High mean but high CV
        r = _make_result(
            mean_depth=50.0,
            stdev_depth=60.0,  # CV = 60/50 = 1.2 > 1.0
            pct_zero=0.5,
            pct_thresholds={20: 97.0},
        )
        params = ClassifyParams()
        status = classify_target(r, params)
        assert status == "UNEVEN"

    def test_uneven_not_triggered_low_mean(self):
        # Low mean: UNEVEN check requires mean > 20
        r = _make_result(
            mean_depth=10.0,
            stdev_depth=20.0,  # CV = 2.0 but mean < 20
            pct_zero=0.5,
            pct_thresholds={20: 50.0},
        )
        params = ClassifyParams()
        status = classify_target(r, params)
        assert status == "LOW_COVERAGE"  # not UNEVEN

    def test_low_exon(self):
        # Main target is good, but one exon is bad
        r = _make_result(pct_zero=0.5, pct_thresholds={20: 97.0})
        bad_exon = _make_result(pct_zero=0.0, pct_thresholds={20: 85.0})
        params = ClassifyParams()
        status = classify_target(r, params, exon_results=[bad_exon])
        assert status == "LOW_EXON"

    def test_classification_order_dropout_before_uneven(self):
        """DROP_OUT takes priority over UNEVEN."""
        r = _make_result(
            mean_depth=50.0,
            stdev_depth=60.0,
            pct_zero=10.0,  # triggers DROP_OUT
            pct_thresholds={20: 97.0},
        )
        params = ClassifyParams()
        status = classify_target(r, params)
        assert status == "DROP_OUT"

    def test_custom_thresholds(self):
        """Custom classification thresholds work."""
        r = _make_result(pct_zero=0.5, pct_thresholds={20: 92.0})
        # Default: pass_pct_ge_20=95 → LOW_COVERAGE
        assert classify_target(r, ClassifyParams()) == "LOW_COVERAGE"
        # Custom: pass_pct_ge_20=90 → PASS
        r.coverage_status = ""  # reset
        assert classify_target(r, ClassifyParams(pass_pct_ge_20=90.0)) == "PASS"
