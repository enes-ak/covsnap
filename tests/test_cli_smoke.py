"""CLI smoke tests for covsnap.

These tests verify that the CLI entry point works end-to-end,
produces expected HTML output, and handles errors properly.
All tests use synthetic BAM data (no real sequencing data or network needed).
"""

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

    def test_gene_mode_produces_html(self, synthetic_bam, tmp_output_dir):
        html_out = str(tmp_output_dir / "report.html")
        main([
            synthetic_bam,
            "BRCA1",
            "-o", html_out,
            "--engine", "samtools",
        ])
        assert os.path.isfile(html_out)
        with open(html_out) as f:
            html = f.read()
        assert "<!DOCTYPE html>" in html
        assert "BRCA1" in html
        assert "SYNTH_SAMPLE" in html

    def test_html_has_required_content(self, synthetic_bam, tmp_output_dir):
        html_out = str(tmp_output_dir / "report.html")
        main([
            synthetic_bam,
            "BRCA1",
            "-o", html_out,
            "--engine", "samtools",
        ])
        with open(html_out) as f:
            html = f.read()
        assert "covsnap" in html.lower()
        assert "hg38" in html
        assert "samtools" in html
        # Self-contained: no external links
        assert 'href="http' not in html
        assert 'src="http' not in html

    def test_default_output_name(self, synthetic_bam, tmp_output_dir, monkeypatch):
        monkeypatch.chdir(tmp_output_dir)
        main([
            synthetic_bam,
            "BRCA1",
            "--engine", "samtools",
        ])
        assert os.path.isfile(str(tmp_output_dir / "covsnap.report.html"))


class TestCLIRegionMode:
    @pytest.fixture(autouse=True)
    def _samtools(self):
        _requires_samtools()

    def test_region_mode(self, synthetic_bam, tmp_output_dir):
        html_out = str(tmp_output_dir / "report.html")
        main([
            synthetic_bam,
            "chr17:43044295-43125482",
            "-o", html_out,
            "--engine", "samtools",
        ])
        assert os.path.isfile(html_out)
        with open(html_out) as f:
            html = f.read()
        assert "chr17" in html


class TestCLIBedMode:
    @pytest.fixture(autouse=True)
    def _samtools(self):
        _requires_samtools()

    def test_bed_mode(self, synthetic_bam, sample_bed, tmp_output_dir):
        html_out = str(tmp_output_dir / "report.html")
        main([
            synthetic_bam,
            "--bed", sample_bed,
            "-o", html_out,
            "--engine", "samtools",
        ])
        assert os.path.isfile(html_out)
        with open(html_out) as f:
            html = f.read()
        assert "region_1" in html
        assert "region_2" in html


class TestCLIExonMode:
    @pytest.fixture(autouse=True)
    def _samtools(self):
        _requires_samtools()

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
