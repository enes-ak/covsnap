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

    _print_summary(settings)

    confirm = _ask(questionary.confirm("Run analysis?", default=True))
    if not confirm:
        return None

    return argparse.Namespace(**settings)


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
