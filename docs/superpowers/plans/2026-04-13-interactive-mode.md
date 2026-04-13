# Interactive Mode Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an interactive TUI that launches when `covsnap` is run without arguments, guiding non-CLI users through input selection with questionary prompts.

**Architecture:** A new `interactive.py` module collects user inputs via questionary prompts and returns an `argparse.Namespace` matching the existing CLI contract. `cli.py:main()` detects no-args invocation and delegates to interactive mode, then continues with the same validation and execution pipeline.

**Tech Stack:** questionary (MIT, built on prompt_toolkit), Python 3.9+

---

## File Structure

| File | Action | Responsibility |
|---|---|---|
| `pyproject.toml` | Modify | Add `questionary` dependency |
| `src/covsnap/interactive.py` | Create | All interactive prompt logic |
| `src/covsnap/cli.py` | Modify | Hook interactive mode on no-args |
| `tests/test_interactive.py` | Create | Unit tests for interactive module |

---

### Task 1: Add questionary dependency

**Files:**
- Modify: `pyproject.toml:28-31`

- [ ] **Step 1: Add questionary to dependencies**

In `pyproject.toml`, add `questionary` to the `dependencies` list:

```toml
dependencies = [
    "pysam>=0.22",
    "numpy>=1.24",
    "questionary>=2.0",
]
```

- [ ] **Step 2: Install and verify**

Run: `pip install -e ".[dev]" && python -c "import questionary; print(questionary.__version__)"`
Expected: prints version >= 2.0

- [ ] **Step 3: Commit**

```bash
git add pyproject.toml
git commit -m "deps: add questionary for interactive mode"
```

---

### Task 2: Create interactive module — core prompt flow

**Files:**
- Create: `src/covsnap/interactive.py`
- Test: `tests/test_interactive.py`

- [ ] **Step 1: Write failing tests**

Create `tests/test_interactive.py`:

```python
"""Tests for interactive mode prompt logic."""

import argparse
from unittest.mock import patch

import pytest

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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_interactive.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'covsnap.interactive'`

- [ ] **Step 3: Implement `collect_inputs()`**

Create `src/covsnap/interactive.py`:

```python
"""Interactive mode for covsnap.

Presents a guided prompt flow using questionary, returning an
argparse.Namespace compatible with the CLI pipeline.
"""

from __future__ import annotations

import argparse
from typing import Optional

import questionary

from covsnap import __version__

# Defaults matching cli.py build_parser()
_DEFAULTS = dict(
    bed=None,
    exons=False,
    reference=None,
    no_index=False,
    engine="auto",
    output="covsnap.report.html",
    threads=4,
    emit_lowcov=False,
    lowcov_threshold=10,
    lowcov_min_len=50,
    max_targets=2000,
    max_total_bp=50_000_000,
    max_bed_bytes=50 * 1024 * 1024,
    on_large_bed="warn_and_clip",
    large_bed_seed=42,
    pct_thresholds="1,5,10,20,30,50,100",
    pass_pct_ge_20=95.0,
    pass_max_pct_zero=1.0,
    dropout_pct_zero=5.0,
    uneven_cv=1.0,
    exon_pct_ge_20=90.0,
    exon_max_pct_zero=5.0,
    verbose=0,
    quiet=False,
)


def _ask(prompt_obj):
    """Run a questionary prompt; return None on Ctrl-C / empty."""
    result = prompt_obj.ask()
    return result


def collect_inputs() -> Optional[argparse.Namespace]:
    """Run the interactive prompt flow.

    Returns an argparse.Namespace matching cli.py expectations,
    or None if the user cancels (Ctrl-C).
    """
    print(f"\n  covsnap v{__version__} — Interactive Mode\n")

    # ── 1. Alignment file ──
    alignment = _ask(questionary.path(
        "BAM/CRAM file path:",
        only_directories=False,
    ))
    if alignment is None:
        return None

    # ── 1b. Reference for CRAM ──
    reference = None
    if alignment.endswith(".cram"):
        reference = _ask(questionary.path(
            "Reference FASTA path (required for CRAM):",
            only_directories=False,
        ))
        if reference is None:
            return None

    # ── 2. Analysis mode ──
    mode = _ask(questionary.select(
        "Analysis mode:",
        choices=["Gene symbol", "Genomic region", "BED file"],
    ))
    if mode is None:
        return None

    # ── 3. Target input ──
    target = None
    bed = None

    if mode == "Gene symbol":
        target = _ask(questionary.text(
            "Gene symbol (e.g. BRCA1, TP53):",
        ))
        if target is None:
            return None

        # ── 4. Exon detail (gene mode only) ──
        exons = _ask(questionary.confirm(
            "Include exon-level detail?",
            default=False,
        ))
        if exons is None:
            return None

    elif mode == "Genomic region":
        target = _ask(questionary.text(
            "Genomic region (e.g. chr17:43044295-43125482):",
        ))
        if target is None:
            return None
        exons = False

    else:  # BED file
        bed = _ask(questionary.path(
            "BED file path:",
            only_directories=False,
        ))
        if bed is None:
            return None
        exons = False

    # ── 5. Engine ──
    engine = _ask(questionary.select(
        "Depth engine:",
        choices=["auto", "samtools", "mosdepth"],
        default="auto",
    ))
    if engine is None:
        return None

    # ── 6. Output path ──
    output = _ask(questionary.text(
        "Output HTML path:",
        default="covsnap.report.html",
    ))
    if output is None:
        return None

    # ── 7. Advanced settings ──
    advanced = _ask(questionary.confirm(
        "Configure advanced settings? (not required — sensible defaults are used)",
        default=False,
    ))
    if advanced is None:
        return None

    # Start with defaults
    settings = dict(_DEFAULTS)
    settings.update(
        alignment=alignment,
        target=target,
        bed=bed,
        exons=exons,
        reference=reference,
        engine=engine,
        output=output,
    )

    if advanced:
        settings.update(_collect_advanced(mode))

    return argparse.Namespace(**settings)


def _collect_advanced(mode: str) -> dict:
    """Collect advanced settings; returns dict of overrides."""
    overrides: dict = {}

    # ── Low-coverage ──
    emit_lowcov = _ask(questionary.confirm(
        "Include low-coverage blocks in report?",
        default=False,
    ))
    if emit_lowcov is None:
        return overrides
    overrides["emit_lowcov"] = emit_lowcov

    if emit_lowcov:
        val = _ask(questionary.text(
            "Low-coverage depth threshold:",
            default="10",
        ))
        if val is not None:
            overrides["lowcov_threshold"] = int(val)

        val = _ask(questionary.text(
            "Minimum low-coverage block length (bp):",
            default="50",
        ))
        if val is not None:
            overrides["lowcov_min_len"] = int(val)

    # ── Classification thresholds ──
    val = _ask(questionary.text(
        "PASS minimum pct_ge_20:",
        default="95.0",
    ))
    if val is not None:
        overrides["pass_pct_ge_20"] = float(val)

    val = _ask(questionary.text(
        "PASS maximum pct_zero:",
        default="1.0",
    ))
    if val is not None:
        overrides["pass_max_pct_zero"] = float(val)

    val = _ask(questionary.text(
        "DROP_OUT pct_zero threshold:",
        default="5.0",
    ))
    if val is not None:
        overrides["dropout_pct_zero"] = float(val)

    val = _ask(questionary.text(
        "UNEVEN CV threshold:",
        default="1.0",
    ))
    if val is not None:
        overrides["uneven_cv"] = float(val)

    # ── Threads ──
    val = _ask(questionary.text(
        "Number of threads:",
        default="4",
    ))
    if val is not None:
        overrides["threads"] = int(val)

    return overrides
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_interactive.py -v`
Expected: all 6 tests PASS

- [ ] **Step 5: Commit**

```bash
git add src/covsnap/interactive.py tests/test_interactive.py
git commit -m "feat: add interactive mode prompt flow"
```

---

### Task 3: Hook interactive mode into CLI

**Files:**
- Modify: `src/covsnap/cli.py:332-335`
- Test: `tests/test_interactive.py` (add integration test)

- [ ] **Step 1: Write failing test**

Add to `tests/test_interactive.py`:

```python
class TestCLINoArgsTriggersInteractive:
    @patch("covsnap.cli.collect_inputs")
    def test_no_args_calls_interactive(self, mock_collect):
        """covsnap with no args should call collect_inputs."""
        mock_collect.return_value = None  # simulate Ctrl-C
        with pytest.raises(SystemExit) as exc_info:
            main([])
        mock_collect.assert_called_once()
        assert exc_info.value.code == 0
```

Add necessary import at top:

```python
from covsnap.cli import main
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_interactive.py::TestCLINoArgsTriggersInteractive -v`
Expected: FAIL — `main([])` currently calls `_validate_args` which errors with "Must provide a target"

- [ ] **Step 3: Modify `cli.py:main()` to detect no-args**

In `src/covsnap/cli.py`, modify the `main()` function. Replace the first few lines:

```python
def main(argv: Optional[list[str]] = None) -> None:
    """CLI entry point."""
    parser = build_parser()

    # No arguments → interactive mode
    if argv is not None and len(argv) == 0:
        _run_interactive()
        return
    if argv is None:
        import sys as _sys
        if len(_sys.argv) == 1:
            _run_interactive()
            return

    args = parser.parse_args(argv)
```

Add the `_run_interactive` helper at the bottom of cli.py (before `if __name__`):

```python
def _run_interactive() -> None:
    """Launch interactive mode and run the pipeline with collected inputs."""
    try:
        from covsnap.interactive import collect_inputs
    except ImportError:
        print(
            "[covsnap] ERROR: Interactive mode requires 'questionary'. "
            "Install it with: pip install questionary",
            file=sys.stderr,
        )
        sys.exit(1)

    args = collect_inputs()
    if args is None:
        # User cancelled (Ctrl-C)
        print("\n[covsnap] Cancelled.", file=sys.stderr)
        sys.exit(0)

    # Fill in fields that _validate_args expects
    if not hasattr(args, "pct_thresholds"):
        args.pct_thresholds = "1,5,10,20,30,50,100"

    # Re-enter normal pipeline from validation onward
    # Set up logging (interactive → default WARNING)
    logging.basicConfig(
        level=logging.WARNING,
        format="[covsnap] %(levelname)s: %(message)s",
        stream=sys.stderr,
    )

    _validate_args(args)
    _run_pipeline(args)
```

Extract the pipeline logic from `main()` into `_run_pipeline(args)` so both CLI and interactive mode share the same code path. This means taking everything after `_validate_args(args)` in the current `main()` (lines 355–683) and moving it into:

```python
def _run_pipeline(args: argparse.Namespace) -> None:
    """Run the analysis pipeline with validated args."""
    thresholds: list[int] = args._parsed_thresholds
    # ... (rest of current main() from line 357 onward, unchanged)
```

And update `main()` to call `_run_pipeline(args)` after validation.

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_interactive.py -v`
Expected: all 7 tests PASS

- [ ] **Step 5: Run full test suite**

Run: `pytest tests/ -v --tb=short`
Expected: all existing tests still pass (no regression)

- [ ] **Step 6: Commit**

```bash
git add src/covsnap/cli.py tests/test_interactive.py
git commit -m "feat: hook interactive mode into CLI on no-args invocation"
```

---

### Task 4: Add summary confirmation before execution

**Files:**
- Modify: `src/covsnap/interactive.py`
- Test: `tests/test_interactive.py`

- [ ] **Step 1: Write failing test**

Add to `tests/test_interactive.py`:

```python
class TestSummaryDisplay:
    @patch("covsnap.interactive.questionary")
    def test_summary_shown_before_run(self, mock_q, capsys):
        """collect_inputs should print a summary before final confirmation."""
        answers = iter([
            "/data/sample.bam",
            "Gene symbol",
            "BRCA1",
            False,
            "auto",
            "covsnap.report.html",
            False,              # no advanced
            True,               # confirm run
        ])

        def make_ask():
            val = next(answers)
            return type("Mock", (), {"ask": lambda self: val})()

        mock_q.path.side_effect = lambda *a, **kw: make_ask()
        mock_q.select.side_effect = lambda *a, **kw: make_ask()
        mock_q.text.side_effect = lambda *a, **kw: make_ask()
        mock_q.confirm.side_effect = lambda *a, **kw: make_ask()

        ns = collect_inputs()

        captured = capsys.readouterr()
        assert "sample.bam" in captured.out
        assert "BRCA1" in captured.out
        assert ns is not None
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_interactive.py::TestSummaryDisplay -v`
Expected: FAIL — summary not printed yet, and current flow doesn't have a final confirm step

- [ ] **Step 3: Add summary display and confirm to `collect_inputs()`**

Add a `_print_summary()` function and a final confirmation prompt at the end of `collect_inputs()`, just before the `return` statement:

```python
def _print_summary(settings: dict) -> None:
    """Print a summary of collected settings."""
    print("\n" + "=" * 50)
    print("  Run Summary")
    print("=" * 50)
    print(f"  Alignment:  {settings['alignment']}")
    if settings.get('target'):
        print(f"  Target:     {settings['target']}")
    if settings.get('bed'):
        print(f"  BED file:   {settings['bed']}")
    if settings.get('reference'):
        print(f"  Reference:  {settings['reference']}")
    print(f"  Exons:      {'yes' if settings.get('exons') else 'no'}")
    print(f"  Engine:     {settings['engine']}")
    print(f"  Output:     {settings['output']}")
    if settings.get('emit_lowcov'):
        print(f"  Low-cov:    threshold={settings['lowcov_threshold']}, min_len={settings['lowcov_min_len']}")
    if settings.get('threads', 4) != 4:
        print(f"  Threads:    {settings['threads']}")
    print("=" * 50 + "\n")
```

In `collect_inputs()`, after building `settings` dict and before `return`, add:

```python
    _print_summary(settings)

    confirm = _ask(questionary.confirm("Run analysis?", default=True))
    if not confirm:
        return None

    return argparse.Namespace(**settings)
```

- [ ] **Step 4: Update earlier tests to include the final confirmation answer**

Each existing test's `answers` iterator needs one more `True` at the end (for the "Run analysis?" confirm). Update all test answer lists — append `True` as the last element.

- [ ] **Step 5: Run all interactive tests**

Run: `pytest tests/test_interactive.py -v`
Expected: all tests PASS

- [ ] **Step 6: Commit**

```bash
git add src/covsnap/interactive.py tests/test_interactive.py
git commit -m "feat: add summary confirmation before analysis runs"
```

---

### Task 5: Manual integration test and final polish

**Files:**
- Modify: `src/covsnap/interactive.py` (if needed)

- [ ] **Step 1: Run full test suite**

Run: `pytest tests/ -v --tb=short`
Expected: all tests pass, no regressions

- [ ] **Step 2: Manual smoke test**

Run: `covsnap` (no arguments) in terminal. Walk through the interactive flow with a synthetic or real BAM file. Verify:
- Tab-completion works for file paths
- CRAM triggers reference prompt
- Gene/Region/BED modes work
- Advanced settings are clearly marked optional
- Summary is displayed before confirmation
- Analysis runs and produces HTML report

- [ ] **Step 3: Commit any fixes**

```bash
git add -u
git commit -m "fix: interactive mode polish from manual testing"
```

- [ ] **Step 4: Final commit**

```bash
git add -A
git commit -m "feat: interactive mode for covsnap (questionary-based TUI)"
```
