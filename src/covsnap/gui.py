"""Tkinter GUI for covsnap.

Provides a simple graphical interface for users who prefer not to use
the command line. Launches when ``covsnap`` is invoked without arguments.
Returns an ``argparse.Namespace`` compatible with the CLI pipeline,
or *None* if the user closes the window.
"""

from __future__ import annotations

import argparse
import os
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from typing import Callable, Optional

from covsnap import __version__

# ── Defaults (must stay in sync with cli.py build_parser) ──────────────
_DEFAULTS = dict(
    bed=None,
    exons=False,
    exon_only=False,
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

    def __init__(self, pipeline_fn: Optional[Callable] = None) -> None:
        self.result: Optional[argparse.Namespace] = None
        self.pipeline_fn = pipeline_fn

        self.root = tk.Tk()
        self.root.title(f"covsnap v{__version__}")
        self.root.resizable(False, False)

        # ── Window icon ──
        self._set_icon()

        # ── Variables ──
        r = self.root
        self.alignment_var = tk.StringVar(master=r)
        self.mode_var = tk.StringVar(master=r, value="gene")
        self.target_var = tk.StringVar(master=r)
        self.bed_var = tk.StringVar(master=r)
        self.reference_var = tk.StringVar(master=r)
        self.exons_var = tk.BooleanVar(master=r, value=False)
        self.exon_only_var = tk.BooleanVar(master=r, value=False)
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

    # ── Icon ──────────────────────────────────────────────────────────

    def _set_icon(self) -> None:
        """Set window icon from the bundled logo."""
        try:
            logo_path = os.path.join(
                os.path.dirname(__file__), "data", "covsnap_logo.png",
            )
            if os.path.exists(logo_path):
                icon = tk.PhotoImage(file=logo_path)
                self.root.iconphoto(True, icon)
                self._icon_ref = icon  # prevent garbage collection
        except Exception:
            pass  # icon is cosmetic, never fail on it

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
        self.target_hint = ttk.Label(root, text="e.g. BRCA1 or BRCA1,TP53,ETFDH", foreground="gray")
        self.target_hint.grid(row=8, column=1, sticky="w", padx=8)

        # ── BED file ──
        row = 9
        self.bed_label = ttk.Label(root, text="BED file:")
        self.bed_label.grid(row=row, column=0, sticky="e", **pad)
        self.bed_entry = ttk.Entry(root, textvariable=self.bed_var, width=45)
        self.bed_entry.grid(row=row, column=1, **pad)
        self.bed_button = ttk.Button(root, text="Browse...", command=self._browse_bed)
        self.bed_button.grid(row=row, column=2, **pad)

        # ── Exons checkboxes ──
        row = 10
        self.exons_check = ttk.Checkbutton(root, text="Show exon-level detail",
                                            variable=self.exons_var)
        self.exons_check.grid(row=row, column=1, sticky="w", **pad)

        row = 11
        self.exon_only_check = ttk.Checkbutton(
            root, text="Exon-only metrics (exclude introns)",
            variable=self.exon_only_var,
        )
        self.exon_only_check.grid(row=row, column=1, sticky="w", **pad)

        # ── Separator ──
        ttk.Separator(root, orient="horizontal").grid(row=12, column=0, columnspan=3, sticky="ew", pady=8)

        # ── Engine ──
        row = 13
        ttk.Label(root, text="Engine:").grid(row=row, column=0, sticky="e", **pad)
        engine_combo = ttk.Combobox(root, textvariable=self.engine_var,
                                     values=["auto", "samtools", "mosdepth"],
                                     state="readonly", width=15)
        engine_combo.grid(row=row, column=1, sticky="w", **pad)

        # ── Output ──
        row = 14
        ttk.Label(root, text="Output file:").grid(row=row, column=0, sticky="e", **pad)
        ttk.Entry(root, textvariable=self.output_var, width=45).grid(row=row, column=1, **pad)
        ttk.Button(root, text="Browse...", command=self._browse_output).grid(row=row, column=2, **pad)

        # ── Advanced settings (collapsible) ──
        ttk.Separator(root, orient="horizontal").grid(row=15, column=0, columnspan=3, sticky="ew", pady=8)

        self.advanced_open = tk.BooleanVar(master=root, value=False)
        self.advanced_toggle = ttk.Button(root, text="Advanced Settings (optional) \u25b6",
                                           command=self._toggle_advanced)
        self.advanced_toggle.grid(row=16, column=0, columnspan=3, pady=(0, 4))

        self.advanced_frame = ttk.LabelFrame(root, text="Advanced Settings", padding=8)
        # Not gridded initially — shown/hidden by _toggle_advanced

        self._build_advanced(self.advanced_frame)

        # ── Status label ──
        self.status_var = tk.StringVar(master=root)
        self.status_label = ttk.Label(root, textvariable=self.status_var, foreground="gray")
        self.status_label.grid(row=19, column=0, columnspan=3, pady=(4, 0))

        # ── Run / Cancel ──
        btn_frame = ttk.Frame(root)
        btn_frame.grid(row=20, column=0, columnspan=3, pady=(4, 12))
        self.run_button = ttk.Button(btn_frame, text="Run Analysis", command=self._on_run)
        self.run_button.pack(side="left", padx=8)
        self.cancel_button = ttk.Button(btn_frame, text="Cancel", command=self._on_cancel)
        self.cancel_button.pack(side="left", padx=8)

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
            self.exon_only_check.grid()
        elif mode == "region":
            self.target_label.config(text="Genomic region:")
            self.target_hint.config(text="e.g. chr17:43044295-43125482")
            self._show_target(True)
            self._show_bed(False)
            self.exons_check.grid_remove()
            self.exon_only_check.grid_remove()
            self.exons_var.set(False)
            self.exon_only_var.set(False)
        else:  # bed
            self._show_target(False)
            self._show_bed(True)
            self.exons_check.grid_remove()
            self.exon_only_check.grid_remove()
            self.exons_var.set(False)
            self.exon_only_var.set(False)

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
            exon_only=self.exon_only_var.get(),
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

        if self.pipeline_fn is None:
            # No pipeline provided (e.g. testing) — just close
            self.root.destroy()
            return

        # Disable controls and show status
        self.run_button.config(state="disabled")
        self.status_var.set("Analyzing... please wait.")
        self.status_label.config(foreground="#0D9488")
        self.root.update_idletasks()

        # Run pipeline in background thread
        thread = threading.Thread(target=self._run_pipeline_thread, daemon=True)
        thread.start()

    def _run_pipeline_thread(self) -> None:
        """Run pipeline in a background thread; schedule UI update when done."""
        try:
            output_path = self.pipeline_fn(self.result)
            self.root.after(0, self._on_pipeline_done, output_path, None)
        except SystemExit as exc:
            # _error() prints to stderr then calls sys.exit — capture the
            # last stderr line so the GUI can show the real message.
            import io, contextlib
            msg = getattr(exc, '_covsnap_message', None)
            if msg is None:
                msg = f"Process exited with code {exc.code}. Check terminal for details."
            self.root.after(0, self._on_pipeline_done, None, msg)
        except Exception as exc:
            self.root.after(0, self._on_pipeline_done, None, str(exc))

    def _on_pipeline_done(self, output_path: Optional[str], error: Optional[str]) -> None:
        """Handle pipeline completion on the main thread."""
        if error:
            self.status_var.set("")
            self.run_button.config(state="normal")
            messagebox.showerror("Error", f"Analysis failed:\n{error}")
            return

        abs_path = os.path.abspath(output_path)
        self.status_var.set(f"Done! Report: {abs_path}")
        self.status_label.config(foreground="#0D9488")
        self.run_button.config(state="normal")
        messagebox.showinfo("Analysis Complete", f"Report saved to:\n{abs_path}")

    def _on_cancel(self) -> None:
        self.result = None
        self.root.destroy()

    # ── Public API ─────────────────────────────────────────────────────

    def run(self) -> Optional[argparse.Namespace]:
        """Show the window and block until the user acts. Returns Namespace or None."""
        self.root.mainloop()
        return self.result


def run_gui(pipeline_fn: Optional[Callable] = None) -> Optional[argparse.Namespace]:
    """Entry point: create and run the GUI.

    If *pipeline_fn* is provided, the pipeline runs inside the GUI window
    (with status feedback). Otherwise the GUI just collects inputs and returns.
    """
    app = CovSnapGUI(pipeline_fn=pipeline_fn)
    return app.run()
