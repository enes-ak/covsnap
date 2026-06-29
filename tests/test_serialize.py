"""Tests for machine-readable serializers (JSON / TSV / MultiQC)."""

import json

from covsnap.metrics import LowCovBlock, TargetResult
from covsnap.report import ClassifyParams, ReportContext
from covsnap.serialize import (
    report_to_dict,
    write_json_report,
    write_tsv_report,
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
        assert set(d) >= {"covsnap_version", "schema_version", "run", "summary", "targets", "genes"}
        assert d["schema_version"] == "1.0"

    def test_run_section_has_all_classify_params(self):
        d = report_to_dict(_make_context())
        cp = d["run"]["classify_params"]
        assert set(cp) == {
            "pass_pct_ge_20",
            "pass_max_pct_zero",
            "dropout_pct_zero",
            "dropout_zero_block_bp",
            "uneven_cv",
            "exon_pct_ge_20",
            "exon_max_pct_zero",
        }
        assert cp["pass_pct_ge_20"] == 95.0

    def test_summary_counts(self):
        results = [_make_result(coverage_status="PASS"), _make_result(target_id="TP53", coverage_status="LOW_COVERAGE")]
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
            "contig": "chr17",
            "start": 43100000,
            "end": 43100060,
            "length": 60,
            "mean_depth": 3.2,
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


class TestTsvWriter:
    def _read(self, path):
        with open(path) as f:
            return [line.rstrip("\n").split("\t") for line in f]

    def test_header_and_row_count(self, tmp_path):
        results = [_make_result(), _make_result(target_id="TP53")]
        path = str(tmp_path / "out.tsv")
        write_tsv_report(path, _make_context(results=results))
        rows = self._read(path)
        assert len(rows) == 3  # header + 2 targets
        assert rows[0][0] == "sample_name"
        assert rows[0][1] == "target_id"

    def test_threshold_columns_present(self, tmp_path):
        path = str(tmp_path / "out.tsv")
        write_tsv_report(path, _make_context())
        header = self._read(path)[0]
        assert "pct_ge_20" in header
        assert "pct_ge_100" in header
        # ascending order of thresholds
        assert header.index("pct_ge_1") < header.index("pct_ge_20")

    def test_row_values(self, tmp_path):
        path = str(tmp_path / "out.tsv")
        write_tsv_report(path, _make_context())
        rows = self._read(path)
        header, row = rows[0], rows[1]
        rec = dict(zip(header, row))
        assert rec["sample_name"] == "SAMPLE_01"
        assert rec["target_id"] == "BRCA1"
        assert rec["pct_ge_20"] == "98.6"
        assert rec["coverage_status"] == "PASS"
        # constant column count across all rows
        assert all(len(r) == len(header) for r in rows)


class TestMultiqcWriter:
    def _load(self, path):
        with open(path) as f:
            return json.load(f)

    def test_required_keys(self, tmp_path):
        from covsnap.serialize import write_multiqc_report

        path = str(tmp_path / "out_mqc.json")
        write_multiqc_report(path, _make_context())
        d = self._load(path)
        assert d["id"] == "covsnap"
        assert d["plot_type"] == "table"
        assert "SAMPLE_01" in d["data"]

    def test_aggregation_weighted_and_counts(self, tmp_path):
        from covsnap.serialize import write_multiqc_report

        # Two targets, different lengths and pct_ge_20, one PASS one LOW_COVERAGE.
        r1 = _make_result(
            target_id="A",
            length_bp=100,
            mean_depth=100.0,
            pct_thresholds={20: 100.0},
            coverage_status="PASS",
        )
        r2 = _make_result(
            target_id="B",
            length_bp=300,
            mean_depth=20.0,
            pct_thresholds={20: 50.0},
            coverage_status="LOW_COVERAGE",
        )
        path = str(tmp_path / "out_mqc.json")
        write_multiqc_report(path, _make_context(results=[r1, r2], thresholds=[20]))
        row = self._load(path)["data"]["SAMPLE_01"]
        assert row["n_targets"] == 2
        assert row["n_pass"] == 1
        assert row["pct_targets_pass"] == 50.0
        # length-weighted mean_depth = (100*100 + 20*300) / 400 = 40.0
        assert row["mean_depth"] == 40.0
        # length-weighted pct_ge_20 = (100*100 + 50*300) / 400 = 62.5
        assert row["pct_ge_20"] == 62.5
        assert row["worst_status"] == "LOW_COVERAGE"

    def test_worst_status_severity_order(self, tmp_path):
        from covsnap.serialize import write_multiqc_report

        results = [
            _make_result(target_id="A", coverage_status="PASS"),
            _make_result(target_id="B", coverage_status="DROP_OUT"),
            _make_result(target_id="C", coverage_status="UNEVEN"),
        ]
        path = str(tmp_path / "out_mqc.json")
        write_multiqc_report(path, _make_context(results=results))
        assert self._load(path)["data"]["SAMPLE_01"]["worst_status"] == "DROP_OUT"

    def test_pct_ge_20_omitted_when_threshold_absent(self, tmp_path):
        from covsnap.serialize import write_multiqc_report

        path = str(tmp_path / "out_mqc.json")
        write_multiqc_report(path, _make_context(thresholds=[1, 10, 30]))
        assert "pct_ge_20" not in self._load(path)["data"]["SAMPLE_01"]


class TestFormatParsing:
    def test_default_single(self):
        from covsnap.serialize import parse_formats

        assert parse_formats("html") == ["html"]

    def test_multi_and_dedupe(self):
        from covsnap.serialize import parse_formats

        assert parse_formats("html, json , json,tsv") == ["html", "json", "tsv"]

    def test_case_insensitive(self):
        from covsnap.serialize import parse_formats

        assert parse_formats("JSON,MultiQC") == ["json", "multiqc"]

    def test_invalid_token_raises(self):
        import pytest

        from covsnap.serialize import parse_formats

        with pytest.raises(ValueError):
            parse_formats("html,pdf")

    def test_empty_raises(self):
        import pytest

        from covsnap.serialize import parse_formats

        with pytest.raises(ValueError):
            parse_formats("")


class TestPathDerivation:
    def test_strips_html_extension(self):
        from covsnap.serialize import derive_output_paths

        paths = derive_output_paths("covsnap.report.html", ["html", "json", "tsv", "multiqc"])
        assert paths["html"] == "covsnap.report.html"
        assert paths["json"] == "covsnap.report.json"
        assert paths["tsv"] == "covsnap.report.tsv"
        assert paths["multiqc"] == "covsnap.report_mqc.json"

    def test_no_html_extension_used_as_stem(self):
        from covsnap.serialize import derive_output_paths

        paths = derive_output_paths("myreport", ["html", "multiqc"])
        assert paths["html"] == "myreport.html"
        assert paths["multiqc"] == "myreport_mqc.json"

    def test_only_requested_formats(self):
        from covsnap.serialize import derive_output_paths

        paths = derive_output_paths("r.html", ["json"])
        assert paths == {"json": "r.json"}
