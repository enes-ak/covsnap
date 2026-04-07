"""Tests for HTML report generation."""

from covsnap.html_report import write_html_report
from covsnap.metrics import TargetResult
from covsnap.report import ClassifyParams, ReportContext


def _make_result(**overrides) -> TargetResult:
    defaults = dict(
        target_id="TEST_GENE",
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
        coverage_status="PASS",
    )
    defaults.update(overrides)
    return TargetResult(**defaults)


def _make_context(results=None, **overrides) -> ReportContext:
    if results is None:
        results = [_make_result()]
    defaults = dict(
        results=results,
        engine_used="samtools",
        engine_version="1.19",
        bam_path="/data/sample.bam",
        sample_name="SAMPLE_01",
        contig_style="chr",
        thresholds=[1, 5, 10, 20, 30, 50, 100],
        run_date="2026-04-07T00:00:00Z",
    )
    defaults.update(overrides)
    return ReportContext(**defaults)


class TestHTMLReportBasic:
    def test_produces_valid_html(self, tmp_path):
        path = str(tmp_path / "report.html")
        ctx = _make_context()
        write_html_report(path, ctx)

        with open(path) as f:
            html = f.read()
        assert html.startswith("<!DOCTYPE html>")
        assert "</html>" in html
        assert "<style>" in html
        assert "</style>" in html

    def test_contains_metadata(self, tmp_path):
        path = str(tmp_path / "report.html")
        ctx = _make_context()
        write_html_report(path, ctx)

        with open(path) as f:
            html = f.read()
        assert "SAMPLE_01" in html
        assert "samtools" in html
        assert "covsnap" in html.lower()

    def test_contains_target_data(self, tmp_path):
        path = str(tmp_path / "report.html")
        ctx = _make_context()
        write_html_report(path, ctx)

        with open(path) as f:
            html = f.read()
        assert "TEST_GENE" in html
        assert "PASS" in html

    def test_self_contained_no_external_links(self, tmp_path):
        path = str(tmp_path / "report.html")
        ctx = _make_context()
        write_html_report(path, ctx)

        with open(path) as f:
            html = f.read()
        assert 'href="http' not in html
        assert 'src="http' not in html
        assert '<link rel="stylesheet"' not in html

    def test_multiple_targets(self, tmp_path):
        path = str(tmp_path / "report.html")
        results = [
            _make_result(target_id="BRCA1", coverage_status="PASS"),
            _make_result(target_id="TP53", coverage_status="LOW_COVERAGE",
                         pct_thresholds={1: 99.0, 5: 95.0, 10: 90.0, 20: 80.0, 30: 70.0, 50: 50.0, 100: 10.0}),
        ]
        ctx = _make_context(results=results)
        write_html_report(path, ctx)

        with open(path) as f:
            html = f.read()
        assert "BRCA1" in html
        assert "TP53" in html
