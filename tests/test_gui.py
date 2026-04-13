"""Tests for covsnap GUI module."""

import argparse
import os
from unittest.mock import patch

import pytest

from covsnap.gui import _DEFAULTS


# Skip all Tk-dependent tests when no display is available
_has_display = bool(os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"))


class TestDefaults:
    """Verify _DEFAULTS dict matches CLI defaults."""

    def test_defaults_has_required_keys(self):
        required = [
            "bed", "exons", "reference", "no_index", "engine", "output",
            "threads", "emit_lowcov", "lowcov_threshold", "lowcov_min_len",
            "pct_thresholds", "pass_pct_ge_20", "pass_max_pct_zero",
            "dropout_pct_zero", "uneven_cv", "verbose", "quiet",
        ]
        for key in required:
            assert key in _DEFAULTS, f"Missing default: {key}"

    def test_default_values(self):
        assert _DEFAULTS["engine"] == "auto"
        assert _DEFAULTS["output"] == "covsnap.report.html"
        assert _DEFAULTS["threads"] == 4
        assert _DEFAULTS["pass_pct_ge_20"] == 95.0


@pytest.mark.skipif(not _has_display, reason="No display available")
class TestGUIWithDisplay:
    """Tests that require a real Tk display."""

    def test_gui_creates_and_cancels(self):
        from covsnap.gui import CovSnapGUI
        gui = CovSnapGUI()
        gui._on_cancel()
        assert gui.result is None

    def test_run_gene_mode(self):
        from covsnap.gui import CovSnapGUI
        gui = CovSnapGUI()
        gui.alignment_var.set("/data/sample.bam")
        gui.mode_var.set("gene")
        gui.target_var.set("BRCA1")
        gui._on_run()
        assert gui.result is not None
        assert gui.result.alignment == "/data/sample.bam"
        assert gui.result.target == "BRCA1"
        assert gui.result.bed is None

    def test_run_bed_mode(self):
        from covsnap.gui import CovSnapGUI
        gui = CovSnapGUI()
        gui.alignment_var.set("/data/sample.bam")
        gui.mode_var.set("bed")
        gui.bed_var.set("/data/targets.bed")
        gui._on_run()
        assert gui.result is not None
        assert gui.result.target is None
        assert gui.result.bed == "/data/targets.bed"

    def test_run_region_mode(self):
        from covsnap.gui import CovSnapGUI
        gui = CovSnapGUI()
        gui.alignment_var.set("/data/sample.bam")
        gui.mode_var.set("region")
        gui.target_var.set("chr17:43044295-43125482")
        gui._on_run()
        assert gui.result is not None
        assert gui.result.target == "chr17:43044295-43125482"
        assert gui.result.exons is False

    def test_run_empty_alignment_shows_error(self):
        from covsnap.gui import CovSnapGUI
        gui = CovSnapGUI()
        gui.alignment_var.set("")
        with patch("covsnap.gui.messagebox.showerror") as mock_err:
            gui._on_run()
            mock_err.assert_called_once()
        assert gui.result is None

    def test_advanced_settings_applied(self):
        from covsnap.gui import CovSnapGUI
        gui = CovSnapGUI()
        gui.alignment_var.set("/data/sample.bam")
        gui.mode_var.set("gene")
        gui.target_var.set("TP53")
        gui.threads_var.set("8")
        gui.pass_pct_var.set("98.0")
        gui._on_run()
        assert gui.result is not None
        assert gui.result.threads == 8
        assert gui.result.pass_pct_ge_20 == 98.0


class TestCLINoArgsTriggersGUI:
    @patch("covsnap.cli._run_interactive")
    def test_no_args_calls_interactive(self, mock_interactive):
        from covsnap.cli import main
        mock_interactive.side_effect = SystemExit(0)
        with pytest.raises(SystemExit) as exc_info:
            main([])
        mock_interactive.assert_called_once()
        assert exc_info.value.code == 0
