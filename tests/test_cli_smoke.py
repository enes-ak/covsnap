"""CLI smoke tests for covsnap.

These tests verify that the CLI entry point works end-to-end,
produces expected output files, and contains required columns/sections.
All tests use synthetic BAM data (no real sequencing data or network needed).
"""

import csv
import json
import os
import shutil

import pytest

from covsnap.cli import main


def _requires_samtools():
    if not shutil.which("samtools"):
        pytest.skip("samtools not installed")


class TestCLIVersion:
    def test_version_flag(self, capsys):
        with pytest.raises(SystemExit) as exc_info:
            main(["--version"])
        assert exc_info.value.code == 0
        captured = capsys.readouterr()
        assert "covsnap" in captured.out.lower() or "0.1.0" in captured.out


class TestCLIValidation:
    def test_no_target_errors(self, synthetic_bam):
        with pytest.raises(SystemExit) as exc_info:
            main([synthetic_bam])
        assert exc_info.value.code == 1

    def test_mutual_exclusion_target_and_bed(self, synthetic_bam, sample_bed):
        with pytest.raises(SystemExit) as exc_info:
            main([synthetic_bam, "BRCA1", "--bed", sample_bed])
        assert exc_info.value.code == 1

    def test_unknown_gene_error(self, synthetic_bam):
        _requires_samtools()
        with pytest.raises(SystemExit) as exc_info:
            main([synthetic_bam, "BRAC1"])
        assert exc_info.value.code == 3

    def test_invalid_region_format(self, synthetic_bam):
        with pytest.raises(SystemExit) as exc_info:
            main([synthetic_bam, "chr17:abc-def"])
        assert exc_info.value.code == 1

    def test_missing_alignment_file(self, tmp_path):
        with pytest.raises(SystemExit) as exc_info:
            main([str(tmp_path / "nonexistent.bam"), "BRCA1"])
        assert exc_info.value.code == 1

    def test_missing_bed_file(self, synthetic_bam):
        with pytest.raises(SystemExit) as exc_info:
            main([synthetic_bam, "--bed", "/nonexistent/targets.bed"])
        assert exc_info.value.code == 1


class TestCLIGeneMode:
    @pytest.fixture(autouse=True)
    def _samtools(self):
        _requires_samtools()

    def test_gene_mode_produces_outputs(self, synthetic_bam, tmp_output_dir):
        raw_out = str(tmp_output_dir / "raw.tsv")
        report_out = str(tmp_output_dir / "report.md")

        main([
            synthetic_bam,
            "BRCA1",
            "--raw-out", raw_out,
            "--report-out", report_out,
            "--engine", "samtools",
        ])
        assert os.path.isfile(raw_out)
        assert os.path.isfile(report_out)

    def test_raw_tsv_has_required_columns(self, synthetic_bam, tmp_output_dir):
        raw_out = str(tmp_output_dir / "raw.tsv")
        report_out = str(tmp_output_dir / "report.md")

        main([
            synthetic_bam,
            "BRCA1",
            "--raw-out", raw_out,
            "--report-out", report_out,
            "--engine", "samtools",
        ])

        with open(raw_out) as f:
            lines = [line for line in f if not line.startswith("#covsnap")]
            reader = csv.DictReader(lines, delimiter="\t")
            rows = list(reader)

            required_cols = [
                "target_id", "contig", "start", "end", "length_bp",
                "mean_depth", "median_depth", "min_depth", "max_depth",
                "pct_zero", "pct_ge_1", "pct_ge_5", "pct_ge_10",
                "pct_ge_20", "pct_ge_30", "pct_ge_50", "pct_ge_100",
                "n_lowcov_blocks", "lowcov_total_bp",
                "engine_used", "bam_path", "sample_name", "build",
                "annotation_version",
            ]
            for col in required_cols:
                assert col in reader.fieldnames, f"Missing column: {col}"

            assert len(rows) >= 1
            assert rows[0]["target_id"] == "BRCA1"
            assert rows[0]["build"] == "hg38"
            assert rows[0]["annotation_version"] == "gencode_v44"
            assert rows[0]["engine_used"] == "samtools"
            assert rows[0]["sample_name"] == "SYNTH_SAMPLE"

    def test_report_has_required_sections(self, synthetic_bam, tmp_output_dir):
        raw_out = str(tmp_output_dir / "raw.tsv")
        report_out = str(tmp_output_dir / "report.md")

        main([
            synthetic_bam,
            "BRCA1",
            "--raw-out", raw_out,
            "--report-out", report_out,
            "--engine", "samtools",
        ])

        with open(report_out) as f:
            content = f.read()
            assert "# covsnap Coverage Report" in content
            assert "BRCA1" in content
            assert "hg38" in content
            assert "GENCODE v44" in content or "gencode_v44" in content
            assert "## Target Summary" in content
            assert "## Classification Heuristics Reference" in content


class TestCLIRegionMode:
    @pytest.fixture(autouse=True)
    def _samtools(self):
        _requires_samtools()

    def test_region_mode(self, synthetic_bam, tmp_output_dir):
        raw_out = str(tmp_output_dir / "raw_region.tsv")
        report_out = str(tmp_output_dir / "report_region.md")

        main([
            synthetic_bam,
            "chr17:43044295-43125482",
            "--raw-out", raw_out,
            "--report-out", report_out,
            "--engine", "samtools",
        ])
        assert os.path.isfile(raw_out)

        with open(raw_out) as f:
            lines = [line for line in f if not line.startswith("#covsnap")]
            reader = csv.DictReader(lines, delimiter="\t")
            rows = list(reader)
            assert len(rows) == 1
            assert rows[0]["contig"] == "chr17"


class TestCLIBedMode:
    @pytest.fixture(autouse=True)
    def _samtools(self):
        _requires_samtools()

    def test_bed_mode(self, synthetic_bam, sample_bed, tmp_output_dir):
        raw_out = str(tmp_output_dir / "raw_bed.tsv")
        report_out = str(tmp_output_dir / "report_bed.md")

        main([
            synthetic_bam,
            "--bed", sample_bed,
            "--raw-out", raw_out,
            "--report-out", report_out,
            "--engine", "samtools",
        ])

        with open(raw_out) as f:
            lines = [line for line in f if not line.startswith("#covsnap")]
            reader = csv.DictReader(lines, delimiter="\t")
            rows = list(reader)
            assert len(rows) == 2


class TestCLIJsonOutput:
    @pytest.fixture(autouse=True)
    def _samtools(self):
        _requires_samtools()

    def test_json_output(self, synthetic_bam, tmp_output_dir):
        raw_out = str(tmp_output_dir / "raw.tsv")
        report_out = str(tmp_output_dir / "report.md")
        json_out = str(tmp_output_dir / "raw.json")

        main([
            synthetic_bam,
            "BRCA1",
            "--raw-out", raw_out,
            "--report-out", report_out,
            "--json-out", json_out,
            "--engine", "samtools",
        ])

        assert os.path.isfile(json_out)
        with open(json_out) as f:
            data = json.load(f)
        assert "targets" in data
        assert data["build"] == "hg38"
        assert len(data["targets"]) >= 1
        assert data["targets"][0]["target_id"] == "BRCA1"


class TestCLILowcov:
    @pytest.fixture(autouse=True)
    def _samtools(self):
        _requires_samtools()

    def test_lowcov_bed_output(self, synthetic_bam, tmp_output_dir):
        raw_out = str(tmp_output_dir / "raw.tsv")
        report_out = str(tmp_output_dir / "report.md")
        lowcov_bed = str(tmp_output_dir / "lowcov.bed")

        main([
            synthetic_bam,
            "BRCA1",
            "--raw-out", raw_out,
            "--report-out", report_out,
            "--emit-lowcov",
            "--lowcov-bed", lowcov_bed,
            "--engine", "samtools",
        ])

        assert os.path.isfile(lowcov_bed)


class TestCLIExonMode:
    @pytest.fixture(autouse=True)
    def _samtools(self):
        _requires_samtools()

    def test_exons_flag_produces_exon_tsv(self, synthetic_bam, tmp_output_dir):
        from covsnap.annotation import _has_full_index

        if not _has_full_index():
            pytest.skip("Full gene/exon index not built")

        raw_out = str(tmp_output_dir / "raw.tsv")
        report_out = str(tmp_output_dir / "report.md")
        exon_out = str(tmp_output_dir / "exons.tsv")

        main([
            synthetic_bam,
            "BRCA1",
            "--exons",
            "--raw-out", raw_out,
            "--report-out", report_out,
            "--exon-out", exon_out,
            "--engine", "samtools",
        ])

        assert os.path.isfile(exon_out)

        # Verify exon TSV structure
        with open(exon_out) as f:
            lines = [line for line in f if not line.startswith("#covsnap")]
            reader = csv.DictReader(lines, delimiter="\t")
            rows = list(reader)

            # BRCA1 MANE Select has 23 exons
            assert len(rows) == 23

            required = [
                "target_id", "exon_id", "exon_number", "contig",
                "start", "end", "length_bp", "mean_depth",
                "median_depth", "pct_zero", "pct_ge_20", "pct_ge_30",
            ]
            for col in required:
                assert col in reader.fieldnames, f"Missing exon column: {col}"

            for row in rows:
                assert row["target_id"] == "BRCA1"
                assert row["contig"] == "chr17"
                assert row["exon_id"].startswith("ENSE")

    def test_exons_report_contains_exon_table(self, synthetic_bam, tmp_output_dir):
        from covsnap.annotation import _has_full_index

        if not _has_full_index():
            pytest.skip("Full gene/exon index not built")

        raw_out = str(tmp_output_dir / "raw.tsv")
        report_out = str(tmp_output_dir / "report.md")
        exon_out = str(tmp_output_dir / "exons.tsv")

        main([
            synthetic_bam,
            "BRCA1",
            "--exons",
            "--raw-out", raw_out,
            "--report-out", report_out,
            "--exon-out", exon_out,
            "--engine", "samtools",
        ])

        with open(report_out) as f:
            content = f.read()
            assert "## Exon-Level Results (BRCA1)" in content
            assert "ENSE" in content  # exon IDs present

    def test_exons_with_bed_is_rejected(self, synthetic_bam, sample_bed):
        """--exons is only valid in gene mode, not with --bed."""
        with pytest.raises(SystemExit) as exc_info:
            main([
                synthetic_bam,
                "--bed", sample_bed,
                "--exons",
            ])
        assert exc_info.value.code == 1


class TestCLIBedGuardrails:
    @pytest.fixture(autouse=True)
    def _samtools(self):
        _requires_samtools()

    def test_large_bed_error_mode(self, synthetic_bam, large_bed):
        with pytest.raises(SystemExit) as exc_info:
            main([
                synthetic_bam,
                "--bed", large_bed,
                "--on-large-bed", "error",
                "--max-targets", "2000",
            ])
        assert exc_info.value.code == 4
