"""Tests for interactive mode prompt logic."""

import argparse
from unittest.mock import patch

import pytest

from covsnap.cli import main
from covsnap.interactive import collect_inputs


class TestCollectInputsGeneMode:
    """Test that collect_inputs returns a correct Namespace for gene mode."""

    @patch("covsnap.interactive.questionary")
    def test_gene_mode_basic(self, mock_q):
        """Minimal gene mode: BAM + gene, no exons, defaults."""
        answers = iter([
            "/data/sample.bam",   # alignment path
            "Gene symbol",        # analysis mode
            "BRCA1",              # gene name
            False,                # exon detail
            "auto",               # engine
            "covsnap.report.html",  # output path (accept default)
            False,                # no advanced settings
        ])

        def make_ask():
            val = next(answers)
            mock_obj = type("Mock", (), {"ask": lambda self: val})()
            return mock_obj

        mock_q.path.side_effect = lambda *a, **kw: make_ask()
        mock_q.select.side_effect = lambda *a, **kw: make_ask()
        mock_q.text.side_effect = lambda *a, **kw: make_ask()
        mock_q.confirm.side_effect = lambda *a, **kw: make_ask()

        ns = collect_inputs()

        assert ns.alignment == "/data/sample.bam"
        assert ns.target == "BRCA1"
        assert ns.bed is None
        assert ns.exons is False
        assert ns.engine == "auto"
        assert ns.output == "covsnap.report.html"


class TestCollectInputsRegionMode:
    @patch("covsnap.interactive.questionary")
    def test_region_mode(self, mock_q):
        answers = iter([
            "/data/sample.bam",
            "Genomic region",
            "chr17:43044295-43125482",
            "auto",
            "covsnap.report.html",
            False,
        ])

        def make_ask():
            val = next(answers)
            return type("Mock", (), {"ask": lambda self: val})()

        mock_q.path.side_effect = lambda *a, **kw: make_ask()
        mock_q.select.side_effect = lambda *a, **kw: make_ask()
        mock_q.text.side_effect = lambda *a, **kw: make_ask()
        mock_q.confirm.side_effect = lambda *a, **kw: make_ask()

        ns = collect_inputs()

        assert ns.target == "chr17:43044295-43125482"
        assert ns.bed is None
        assert ns.exons is False


class TestCollectInputsBEDMode:
    @patch("covsnap.interactive.questionary")
    def test_bed_mode(self, mock_q):
        answers = iter([
            "/data/sample.bam",
            "BED file",
            "/data/targets.bed",
            "samtools",
            "my_report.html",
            False,
        ])

        def make_ask():
            val = next(answers)
            return type("Mock", (), {"ask": lambda self: val})()

        mock_q.path.side_effect = lambda *a, **kw: make_ask()
        mock_q.select.side_effect = lambda *a, **kw: make_ask()
        mock_q.text.side_effect = lambda *a, **kw: make_ask()
        mock_q.confirm.side_effect = lambda *a, **kw: make_ask()

        ns = collect_inputs()

        assert ns.target is None
        assert ns.bed == "/data/targets.bed"


class TestCollectInputsAdvanced:
    @patch("covsnap.interactive.questionary")
    def test_advanced_settings(self, mock_q):
        answers = iter([
            "/data/sample.bam",
            "Gene symbol",
            "TP53",
            True,               # exons
            "samtools",
            "tp53.html",
            True,               # advanced settings
            True,               # low-coverage
            "20",               # lowcov threshold
            "100",              # lowcov min len
            "98.0",             # pass_pct_ge_20
            "0.5",              # pass_max_pct_zero
            "3.0",              # dropout_pct_zero
            "0.8",              # uneven_cv
            "8",                # threads
        ])

        def make_ask():
            val = next(answers)
            return type("Mock", (), {"ask": lambda self: val})()

        mock_q.path.side_effect = lambda *a, **kw: make_ask()
        mock_q.select.side_effect = lambda *a, **kw: make_ask()
        mock_q.text.side_effect = lambda *a, **kw: make_ask()
        mock_q.confirm.side_effect = lambda *a, **kw: make_ask()

        ns = collect_inputs()

        assert ns.exons is True
        assert ns.emit_lowcov is True
        assert ns.lowcov_threshold == 20
        assert ns.lowcov_min_len == 100
        assert ns.pass_pct_ge_20 == 98.0
        assert ns.pass_max_pct_zero == 0.5
        assert ns.dropout_pct_zero == 3.0
        assert ns.uneven_cv == 0.8
        assert ns.threads == 8


class TestCollectInputsCRAM:
    @patch("covsnap.interactive.questionary")
    def test_cram_asks_reference(self, mock_q):
        answers = iter([
            "/data/sample.cram",
            "/data/hg38.fa",        # reference (asked because .cram)
            "Gene symbol",
            "BRCA1",
            False,
            "auto",
            "covsnap.report.html",
            False,
        ])

        def make_ask():
            val = next(answers)
            return type("Mock", (), {"ask": lambda self: val})()

        mock_q.path.side_effect = lambda *a, **kw: make_ask()
        mock_q.select.side_effect = lambda *a, **kw: make_ask()
        mock_q.text.side_effect = lambda *a, **kw: make_ask()
        mock_q.confirm.side_effect = lambda *a, **kw: make_ask()

        ns = collect_inputs()

        assert ns.alignment == "/data/sample.cram"
        assert ns.reference == "/data/hg38.fa"


class TestUserCancellation:
    @patch("covsnap.interactive.questionary")
    def test_ctrl_c_returns_none(self, mock_q):
        mock_obj = type("Mock", (), {"ask": lambda self: None})()
        mock_q.path.return_value = mock_obj

        result = collect_inputs()
        assert result is None


class TestCLINoArgsTriggersInteractive:
    @patch("covsnap.cli.collect_inputs")
    def test_no_args_calls_interactive(self, mock_collect):
        """covsnap with no args should call collect_inputs."""
        mock_collect.return_value = None  # simulate Ctrl-C
        with pytest.raises(SystemExit) as exc_info:
            main([])
        mock_collect.assert_called_once()
        assert exc_info.value.code == 0
