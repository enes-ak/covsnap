"""Tkinter GUI for covsnap.

Provides a simple graphical interface for users who prefer not to use
the command line. Launches when ``covsnap`` is invoked without arguments.
Returns an ``argparse.Namespace`` compatible with the CLI pipeline,
or *None* if the user closes the window.
"""

from __future__ import annotations

import argparse
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from typing import Optional

from covsnap import __version__

# ── Defaults (must stay in sync with cli.py build_parser) ──────────────
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


class CovSnapGUI:
    """Main application window."""

    def __init__(self) -> None:
        self.result: Optional[argparse.Namespace] = None

        self.root = tk.Tk()
        self.root.title(f"covsnap v{__version__}")
        self.root.resizable(False, False)

        # ── Variables ──
        r = self.root
        self.alignment_var = tk.StringVar(master=r)
        self.mode_var = tk.StringVar(master=r, value="gene")
        self.target_var = tk.StringVar(master=r)
        self.bed_var = tk.StringVar(master=r)
        self.reference_var = tk.StringVar(master=r)
        self.exons_var = tk.BooleanVar(master=r, value=False)
        self.engine_var = tk.StringVar(master=r, value="auto")
        self.output_var = tk.StringVar(master=r, value="covsnap.report.html")
        # Advanced
        self.emit_lowcov_var = tk.BooleanVar(master=r, value=False)
        self.lowcov_threshold_var = tk.StringVar(master=r, value="10")
        self.lowcov_min_len_var = tk.StringVar(master=r, value="50")
        self.pass_pct_var = tk.StringVar(master=r, value="95.0")
        self.pass_zero_var = tk.StringVar(master=r, value="1.0")
        self.dropout_var = tk.StringVar(master=r, value="5.0")
        self.uneven_var = tk.StringVar(master=r, value="1.0")
        self.threads_var = tk.StringVar(master=r, value="4")

        self._build_ui()
        self._on_mode_change()

    # ── UI Construction ────────────────────────────────────────────────

    def _build_ui(self) -> None:
        pad = dict(padx=8, pady=3)
        root = self.root

        # Title
        title = ttk.Label(root, text=f"covsnap v{__version__}", font=("", 13, "bold"))
        title.grid(row=0, column=0, columnspan=3, pady=(12, 4))
        subtitle = ttk.Label(root, text="Coverage inspector for targeted sequencing QC")
        subtitle.grid(row=1, column=0, columnspan=3, pady=(0, 10))

        # ── Alignment ──
        row = 2
        ttk.Label(root, text="BAM / CRAM file:").grid(row=row, column=0, sticky="e", **pad)
        ttk.Entry(root, textvariable=self.alignment_var, width=45).grid(row=row, column=1, **pad)
        ttk.Button(root, text="Browse...", command=self._browse_alignment).grid(row=row, column=2, **pad)

        # ── Reference (CRAM) ──
        row = 3
        ttk.Label(root, text="Reference FASTA:").grid(row=row, column=0, sticky="e", **pad)
        ttk.Entry(root, textvariable=self.reference_var, width=45).grid(row=row, column=1, **pad)
        ttk.Button(root, text="Browse...", command=self._browse_reference).grid(row=row, column=2, **pad)
        self.ref_hint = ttk.Label(root, text="(only needed for CRAM)", foreground="gray")
        self.ref_hint.grid(row=4, column=1, sticky="w", padx=8)

        # ── Separator ──
        ttk.Separator(root, orient="horizontal").grid(row=5, column=0, columnspan=3, sticky="ew", pady=8)

        # ── Mode selection ──
        row = 6
        ttk.Label(root, text="Analysis mode:").grid(row=row, column=0, sticky="e", **pad)
        mode_frame = ttk.Frame(root)
        mode_frame.grid(row=row, column=1, sticky="w", **pad)
        for text, val in [("Gene symbol", "gene"), ("Genomic region", "region"), ("BED file", "bed")]:
            ttk.Radiobutton(mode_frame, text=text, variable=self.mode_var, value=val,
                            command=self._on_mode_change).pack(side="left", padx=(0, 12))

        # ── Target (gene / region) ──
        row = 7
        self.target_label = ttk.Label(root, text="Gene symbol:")
        self.target_label.grid(row=row, column=0, sticky="e", **pad)
        self.target_entry = ttk.Entry(root, textvariable=self.target_var, width=45)
        self.target_entry.grid(row=row, column=1, **pad)
        self.target_hint = ttk.Label(root, text="e.g. BRCA1, TP53", foreground="gray")
        self.target_hint.grid(row=8, column=1, sticky="w", padx=8)

        # ── BED file ──
        row = 9
        self.bed_label = ttk.Label(root, text="BED file:")
        self.bed_label.grid(row=row, column=0, sticky="e", **pad)
        self.bed_entry = ttk.Entry(root, textvariable=self.bed_var, width=45)
        self.bed_entry.grid(row=row, column=1, **pad)
        self.bed_button = ttk.Button(root, text="Browse...", command=self._browse_bed)
        self.bed_button.grid(row=row, column=2, **pad)

        # ── Exons checkbox ──
        row = 10
        self.exons_check = ttk.Checkbutton(root, text="Include exon-level detail",
                                            variable=self.exons_var)
        self.exons_check.grid(row=row, column=1, sticky="w", **pad)

        # ── Separator ──
        ttk.Separator(root, orient="horizontal").grid(row=11, column=0, columnspan=3, sticky="ew", pady=8)

        # ── Engine ──
        row = 12
        ttk.Label(root, text="Engine:").grid(row=row, column=0, sticky="e", **pad)
        engine_combo = ttk.Combobox(root, textvariable=self.engine_var,
                                     values=["auto", "samtools", "mosdepth"],
                                     state="readonly", width=15)
        engine_combo.grid(row=row, column=1, sticky="w", **pad)

        # ── Output ──
        row = 13
        ttk.Label(root, text="Output file:").grid(row=row, column=0, sticky="e", **pad)
        ttk.Entry(root, textvariable=self.output_var, width=45).grid(row=row, column=1, **pad)
        ttk.Button(root, text="Browse...", command=self._browse_output).grid(row=row, column=2, **pad)

        # ── Advanced settings (collapsible) ──
        ttk.Separator(root, orient="horizontal").grid(row=14, column=0, columnspan=3, sticky="ew", pady=8)

        self.advanced_open = tk.BooleanVar(master=root, value=False)
        self.advanced_toggle = ttk.Button(root, text="Advanced Settings (optional) \u25b6",
                                           command=self._toggle_advanced)
        self.advanced_toggle.grid(row=15, column=0, columnspan=3, pady=(0, 4))

        self.advanced_frame = ttk.LabelFrame(root, text="Advanced Settings", padding=8)
        # Not gridded initially — shown/hidden by _toggle_advanced

        self._build_advanced(self.advanced_frame)

        # ── Run / Cancel ──
        btn_frame = ttk.Frame(root)
        btn_frame.grid(row=17, column=0, columnspan=3, pady=(8, 12))
        ttk.Button(btn_frame, text="Run Analysis", command=self._on_run).pack(side="left", padx=8)
        ttk.Button(btn_frame, text="Cancel", command=self._on_cancel).pack(side="left", padx=8)

    def _build_advanced(self, frame: ttk.LabelFrame) -> None:
        pad = dict(padx=6, pady=2)

        # Low-coverage
        ttk.Checkbutton(frame, text="Emit low-coverage blocks",
                         variable=self.emit_lowcov_var).grid(row=0, column=0, columnspan=2, sticky="w", **pad)

        ttk.Label(frame, text="Low-cov threshold:").grid(row=1, column=0, sticky="e", **pad)
        ttk.Entry(frame, textvariable=self.lowcov_threshold_var, width=10).grid(row=1, column=1, sticky="w", **pad)

        ttk.Label(frame, text="Low-cov min length:").grid(row=2, column=0, sticky="e", **pad)
        ttk.Entry(frame, textvariable=self.lowcov_min_len_var, width=10).grid(row=2, column=1, sticky="w", **pad)

        ttk.Separator(frame, orient="horizontal").grid(row=3, column=0, columnspan=2, sticky="ew", pady=6)

        # Classification thresholds
        ttk.Label(frame, text="PASS min pct_ge_20:").grid(row=4, column=0, sticky="e", **pad)
        ttk.Entry(frame, textvariable=self.pass_pct_var, width=10).grid(row=4, column=1, sticky="w", **pad)

        ttk.Label(frame, text="PASS max pct_zero:").grid(row=5, column=0, sticky="e", **pad)
        ttk.Entry(frame, textvariable=self.pass_zero_var, width=10).grid(row=5, column=1, sticky="w", **pad)

        ttk.Label(frame, text="DROP_OUT pct_zero:").grid(row=6, column=0, sticky="e", **pad)
        ttk.Entry(frame, textvariable=self.dropout_var, width=10).grid(row=6, column=1, sticky="w", **pad)

        ttk.Label(frame, text="UNEVEN CV:").grid(row=7, column=0, sticky="e", **pad)
        ttk.Entry(frame, textvariable=self.uneven_var, width=10).grid(row=7, column=1, sticky="w", **pad)

        ttk.Separator(frame, orient="horizontal").grid(row=8, column=0, columnspan=2, sticky="ew", pady=6)

        ttk.Label(frame, text="Threads:").grid(row=9, column=0, sticky="e", **pad)
        ttk.Entry(frame, textvariable=self.threads_var, width=10).grid(row=9, column=1, sticky="w", **pad)

    # ── Mode switching ─────────────────────────────────────────────────

    def _on_mode_change(self) -> None:
        mode = self.mode_var.get()
        if mode == "gene":
            self.target_label.config(text="Gene symbol:")
            self.target_hint.config(text="e.g. BRCA1, TP53")
            self._show_target(True)
            self._show_bed(False)
            self.exons_check.grid()
        elif mode == "region":
            self.target_label.config(text="Genomic region:")
            self.target_hint.config(text="e.g. chr17:43044295-43125482")
            self._show_target(True)
            self._show_bed(False)
            self.exons_check.grid_remove()
            self.exons_var.set(False)
        else:  # bed
            self._show_target(False)
            self._show_bed(True)
            self.exons_check.grid_remove()
            self.exons_var.set(False)

    def _show_target(self, show: bool) -> None:
        if show:
            self.target_label.grid()
            self.target_entry.grid()
            self.target_hint.grid()
        else:
            self.target_label.grid_remove()
            self.target_entry.grid_remove()
            self.target_hint.grid_remove()

    def _show_bed(self, show: bool) -> None:
        if show:
            self.bed_label.grid()
            self.bed_entry.grid()
            self.bed_button.grid()
        else:
            self.bed_label.grid_remove()
            self.bed_entry.grid_remove()
            self.bed_button.grid_remove()

    def _toggle_advanced(self) -> None:
        if self.advanced_open.get():
            self.advanced_frame.grid_remove()
            self.advanced_open.set(False)
            self.advanced_toggle.config(text="Advanced Settings (optional) \u25b6")
        else:
            self.advanced_frame.grid(row=16, column=0, columnspan=3, padx=12, pady=(0, 4), sticky="ew")
            self.advanced_open.set(True)
            self.advanced_toggle.config(text="Advanced Settings (optional) \u25bc")

    # ── File dialogs ───────────────────────────────────────────────────

    def _browse_alignment(self) -> None:
        path = filedialog.askopenfilename(
            title="Select BAM/CRAM file",
            filetypes=[("BAM files", "*.bam"), ("CRAM files", "*.cram"), ("All files", "*.*")],
        )
        if path:
            self.alignment_var.set(path)

    def _browse_reference(self) -> None:
        path = filedialog.askopenfilename(
            title="Select Reference FASTA",
            filetypes=[("FASTA files", "*.fa *.fasta *.fa.gz"), ("All files", "*.*")],
        )
        if path:
            self.reference_var.set(path)

    def _browse_bed(self) -> None:
        path = filedialog.askopenfilename(
            title="Select BED file",
            filetypes=[("BED files", "*.bed"), ("All files", "*.*")],
        )
        if path:
            self.bed_var.set(path)

    def _browse_output(self) -> None:
        path = filedialog.asksaveasfilename(
            title="Save HTML report as",
            defaultextension=".html",
            filetypes=[("HTML files", "*.html"), ("All files", "*.*")],
        )
        if path:
            self.output_var.set(path)

    # ── Validation & Run ───────────────────────────────────────────────

    def _on_run(self) -> None:
        alignment = self.alignment_var.get().strip()
        if not alignment:
            messagebox.showerror("Error", "Please select a BAM/CRAM file.")
            return

        mode = self.mode_var.get()
        target = None
        bed = None

        if mode in ("gene", "region"):
            target = self.target_var.get().strip()
            if not target:
                label = "gene symbol" if mode == "gene" else "genomic region"
                messagebox.showerror("Error", f"Please enter a {label}.")
                return
        else:
            bed = self.bed_var.get().strip()
            if not bed:
                messagebox.showerror("Error", "Please select a BED file.")
                return

        reference = self.reference_var.get().strip() or None

        try:
            threads = int(self.threads_var.get())
            lowcov_threshold = int(self.lowcov_threshold_var.get())
            lowcov_min_len = int(self.lowcov_min_len_var.get())
            pass_pct = float(self.pass_pct_var.get())
            pass_zero = float(self.pass_zero_var.get())
            dropout = float(self.dropout_var.get())
            uneven = float(self.uneven_var.get())
        except ValueError:
            messagebox.showerror("Error", "Invalid number in advanced settings.")
            return

        settings = dict(_DEFAULTS)
        settings.update(
            alignment=alignment,
            target=target,
            bed=bed,
            exons=self.exons_var.get(),
            reference=reference,
            engine=self.engine_var.get(),
            output=self.output_var.get().strip() or "covsnap.report.html",
            threads=threads,
            emit_lowcov=self.emit_lowcov_var.get(),
            lowcov_threshold=lowcov_threshold,
            lowcov_min_len=lowcov_min_len,
            pass_pct_ge_20=pass_pct,
            pass_max_pct_zero=pass_zero,
            dropout_pct_zero=dropout,
            uneven_cv=uneven,
        )

        self.result = argparse.Namespace(**settings)
        self.root.destroy()

    def _on_cancel(self) -> None:
        self.result = None
        self.root.destroy()

    # ── Public API ─────────────────────────────────────────────────────

    def run(self) -> Optional[argparse.Namespace]:
        """Show the window and block until the user acts. Returns Namespace or None."""
        self.root.mainloop()
        return self.result


def run_gui() -> Optional[argparse.Namespace]:
    """Entry point: create and run the GUI, return collected args or None."""
    app = CovSnapGUI()
    return app.run()
