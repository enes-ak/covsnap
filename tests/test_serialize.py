"""Tests for machine-readable serializers (JSON / TSV / MultiQC)."""

import json

from covsnap.metrics import LowCovBlock, TargetResult
from covsnap.report import ClassifyParams, ReportContext
from covsnap.serialize import (
    report_to_dict,
    write_json_report,
)


def _make_result(**overrides) -> TargetResult:
    defaults = dict(
        target_id="BRCA1",
        contig="chr17",
        start=43044294,
        end=43125483,
        length_bp=81189,
        mean_depth=142.3,
        median_depth=140.0,
        min_depth=0,
        max_depth=310,
        stdev_depth=38.1,
        pct_zero=0.12,
        pct_thresholds={1: 99.9, 5: 99.8, 10: 99.5, 20: 98.6, 30: 97.0, 50: 92.0, 100: 71.0},
        n_lowcov_blocks=0,
        lowcov_total_bp=0,
        lowcov_blocks=[],
        histogram={0: 96, 140: 12000},
        coverage_status="PASS",
    )
    defaults.update(overrides)
    return TargetResult(**defaults)


def _make_context(results=None, **overrides) -> ReportContext:
    if results is None:
        results = [_make_result()]
    defaults = dict(
        results=results,
        engine_used="mosdepth",
        engine_version="0.3.8",
        bam_path="/data/sample.bam",
        sample_name="SAMPLE_01",
        contig_style="chr",
        thresholds=[1, 5, 10, 20, 30, 50, 100],
        run_date="2026-06-29",
        classify_params=ClassifyParams(),
    )
    defaults.update(overrides)
    return ReportContext(**defaults)


class TestReportToDict:
    def test_top_level_keys(self):
        d = report_to_dict(_make_context())
        assert set(d) >= {
            "covsnap_version", "schema_version", "run", "summary", "targets", "genes"
        }
        assert d["schema_version"] == "1.0"

    def test_run_section_has_all_classify_params(self):
        d = report_to_dict(_make_context())
        cp = d["run"]["classify_params"]
        assert set(cp) == {
            "pass_pct_ge_20", "pass_max_pct_zero", "dropout_pct_zero",
            "dropout_zero_block_bp", "uneven_cv", "exon_pct_ge_20", "exon_max_pct_zero",
        }
        assert cp["pass_pct_ge_20"] == 95.0

    def test_summary_counts(self):
        results = [_make_result(coverage_status="PASS"),
                   _make_result(target_id="TP53", coverage_status="LOW_COVERAGE")]
        d = report_to_dict(_make_context(results=results))
        assert d["summary"]["n_targets"] == 2
        assert d["summary"]["n_pass"] == 1
        assert d["summary"]["status_counts"]["LOW_COVERAGE"] == 1

    def test_threshold_and_histogram_keys_are_strings(self):
        d = report_to_dict(_make_context())
        t = d["targets"][0]
        assert all(isinstance(k, str) for k in t["pct_thresholds"])
        assert t["pct_thresholds"]["20"] == 98.6
        assert all(isinstance(k, str) for k in t["histogram"])

    def test_lowcov_blocks_serialized(self):
        r = _make_result(
            n_lowcov_blocks=1,
            lowcov_total_bp=60,
            lowcov_blocks=[LowCovBlock(start=43100000, end=43100060, depth_sum=192, length=60)],
        )
        d = report_to_dict(_make_context(results=[r]))
        blk = d["targets"][0]["lowcov_blocks"][0]
        assert blk == {
            "contig": "chr17", "start": 43100000, "end": 43100060,
            "length": 60, "mean_depth": 3.2,
        }

    def test_exons_nested_when_present(self):
        gene = _make_result(target_id="BRCA1")
        exon = _make_result(target_id="BRCA1_exon1", length_bp=213)
        ctx = _make_context(
            results=[gene],
            exon_results={"BRCA1": [exon]},
            exon_statuses={"BRCA1": ["OK"]},
            exon_metadata={"BRCA1": [{"exon_number": 1, "exon_id": "ENSE0001"}]},
        )
        d = report_to_dict(ctx)
        exons = d["targets"][0]["exons"]
        assert len(exons) == 1
        assert exons[0]["exon_number"] == 1
        assert exons[0]["status"] == "OK"

    def test_exons_absent_when_no_exon_data(self):
        d = report_to_dict(_make_context())
        assert "exons" not in d["targets"][0]

    def test_genes_populated_only_in_region_mode(self):
        assert report_to_dict(_make_context())["genes"] == []
        gene = _make_result(target_id="BRCA1")
        d = report_to_dict(_make_context(gene_results=[gene]))
        assert len(d["genes"]) == 1
        assert d["genes"][0]["target_id"] == "BRCA1"


class TestJsonWriter:
    def test_write_json_round_trips(self, tmp_path):
        path = str(tmp_path / "out.json")
        write_json_report(path, _make_context())
        with open(path) as f:
            loaded = json.load(f)
        assert loaded["run"]["sample_name"] == "SAMPLE_01"
        assert loaded["targets"][0]["target_id"] == "BRCA1"
