# covsnap

[![Bioconda](https://img.shields.io/conda/vn/bioconda/covsnap.svg)](https://anaconda.org/bioconda/covsnap)
[![PyPI](https://img.shields.io/pypi/v/covsnap.svg)](https://pypi.org/project/covsnap/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18732742.svg)](https://doi.org/10.5281/zenodo.18732742)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**Coverage inspector for targeted sequencing QC (hg38)**

covsnap computes per-target (and optionally per-exon) depth-of-coverage metrics from BAM/CRAM files aligned to the human reference genome **hg38**. It produces a self-contained interactive HTML report with automated PASS/FAIL classification heuristics — designed for clinical and research sequencing QC workflows.

---

## Screenshots

### Interactive HTML Report

<p align="center">
  <img src="docs/screenshots/html-report-full.png" alt="covsnap HTML report showing summary cards, exon coverage chart, and PASS/FAIL classification" width="800">
</p>

### Graphical Interface (GUI)

<p align="center">
  <img src="docs/screenshots/gui.png" alt="covsnap Tkinter GUI with file pickers, mode selection, and analysis options" width="500">
</p>

Run `covsnap` with no arguments to launch the GUI.

---

## Key Features

- **Graphical interface** — Run `covsnap` with no arguments to launch a Tkinter GUI with file pickers, mode selection, and progress feedback. Works on Linux, macOS, and Windows.
- **Gene-aware analysis** — Look up genes by symbol (e.g. `BRCA1`) or analyze multiple genes at once with a comma-separated list (e.g. `BRCA1,TP53,ETFDH`). Ships with a built-in dictionary of ~60 clinically relevant genes and an optional full GENCODE v44 tabix index covering 62,700+ genes.
- **Exon-level resolution** — Per-exon depth metrics via the `--exons` flag using MANE Select transcripts from GENCODE v44.
- **Region and BED modes** — Accepts genomic coordinates (`chr17:43044295-43125482`) or a BED file of arbitrary target intervals. Region mode auto-discovers overlapping genes and exons.
- **Interactive HTML report** — Single self-contained HTML file with summary cards, exon bar charts with smooth color gradients, accordion details, glossary, and PASS/FAIL classifications.
- **Streaming architecture** — O(1) memory per target using Welford's online algorithm for mean/variance and histogram-based exact median. No per-base depth arrays are ever held in memory.
- **Parallel execution** — Concurrent samtools and region/exon analysis for faster results.
- **Dual engine support** — Prefers [mosdepth](https://github.com/brentp/mosdepth) when available; falls back to `samtools depth`.
- **Contig auto-detection** — Transparently handles both `chr`-prefixed (UCSC) and non-prefixed (Ensembl/1000G) BAM contig naming.
- **Gene alias resolution** — Common aliases like `HER2 -> ERBB2` and `P53 -> TP53` are resolved automatically, with fuzzy suggestions for typos.
- **BED guardrails** — Configurable limits on target count, total bases, and file size to prevent accidental whole-exome/whole-genome runs.
- **Classification heuristics** — Automated PASS / DROP_OUT / UNEVEN / LOW_EXON / LOW_COVERAGE calls with tunable thresholds.

---

## Installation

### From Bioconda (recommended)

```bash
conda install -c bioconda covsnap
```

### From PyPI

```bash
pip install covsnap
```

### From source

```bash
git clone https://github.com/enes-ak/covsnap.git
cd covsnap
pip install .
```

### With development/test dependencies

```bash
pip install ".[dev]"
```

### Runtime requirements

| Dependency | Version | Required? |
|---|---|---|
| Python | >= 3.9 | Yes |
| pysam | >= 0.22 | Yes |
| numpy | >= 1.24 | Yes |
| samtools | any recent | Yes (engine) |
| mosdepth | >= 0.3 | Optional (preferred engine) |

> **Note:** At least one of `samtools` or `mosdepth` must be on your `$PATH`. When `--engine auto` (the default), covsnap prefers mosdepth and falls back to samtools.

---

## Quick Start

### Graphical interface

Run covsnap with no arguments to launch the GUI:

```bash
covsnap
```

A window opens where you can select your BAM file, choose analysis mode, configure options, and run the analysis — all without typing commands.

### Gene mode

Analyze coverage for a gene by name:

```bash
covsnap sample.bam BRCA1
```

This produces `covsnap.report.html` — an interactive HTML report with coverage metrics and PASS/FAIL classification.

### Multiple genes

Analyze several genes in a single run with a comma-separated list:

```bash
covsnap sample.bam BRCA1,TP53,ETFDH --exons
```

### With exon-level detail

```bash
covsnap sample.bam BRCA1 --exons
```

### Region mode

Specify an explicit genomic region (1-based inclusive coordinates). Overlapping genes and exons are auto-discovered:

```bash
covsnap sample.bam chr17:43044295-43125482
```

### BED mode

Use a BED file of target intervals:

```bash
covsnap sample.bam --bed targets.bed
```

### Custom output path

```bash
covsnap sample.bam BRCA1 -o my_report.html
```

### CRAM files

```bash
covsnap sample.cram BRCA1 --reference hg38.fa
```

---

## HTML Report

covsnap produces a single self-contained HTML file (no external dependencies) containing:

- **Summary cards** — key metrics at a glance (mean depth, coverage breadth, classification)
- **Exon bar chart** — per-exon coverage with smooth HSL color gradient (red → amber → teal)
- **Accordion details** — expandable per-target and per-exon metrics
- **Low-coverage blocks** — contiguous regions below threshold (when `--emit-lowcov` is used)
- **Glossary** — definitions of all metrics and classification terms

---

## Classification Heuristics

Each target is classified using ordered heuristics (first match wins):

| Status | Condition |
|---|---|
| **DROP_OUT** | `pct_zero > 5%` OR any zero-coverage block >= 500 bp |
| **UNEVEN** | `mean_depth > 20` AND coefficient of variation > 1.0 |
| **LOW_EXON** | Any exon with `pct_ge_20 < 90%` or `pct_zero > 5%` (exon mode only) |
| **LOW_COVERAGE** | `pct_ge_20 < 95%` |
| **PASS** | `pct_ge_20 >= 95%` AND `pct_zero <= 1%` |

All thresholds are tunable via CLI flags:

```bash
covsnap sample.bam BRCA1 \
    --pass-pct-ge-20 98.0 \
    --pass-max-pct-zero 0.5 \
    --dropout-pct-zero 3.0 \
    --uneven-cv 0.8
```

---

## BED Guardrails

When using `--bed`, covsnap enforces limits to prevent accidental whole-exome/whole-genome processing:

| Parameter | Default | Flag |
|---|---|---|
| Max target intervals | 2,000 | `--max-targets` |
| Max total base pairs | 50 Mb | `--max-total-bp` |
| Max BED file size | 50 MB | `--max-bed-bytes` |

When limits are exceeded, the behavior is controlled by `--on-large-bed`:

| Mode | Behavior |
|---|---|
| `error` | Exit with code 4 |
| `warn_and_clip` (default) | Keep the first N targets that fit within limits |
| `warn_and_sample` | Reservoir sample N targets (deterministic with `--large-bed-seed`) |

---

## Building the Full Gene Index

The package ships with a built-in dictionary of ~60 clinically relevant genes. For access to the full GENCODE v44 catalog (62,700+ genes, 201,000+ MANE Select exons), build the tabix index:

```bash
# Download GENCODE v44 GTF (requires ~1.5 GB)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz

# Build the index
python scripts/build_gene_index.py gencode.v44.annotation.gtf.gz

# Files are written to src/covsnap/data/
```

This creates:
- `hg38_genes.tsv.gz` + `.tbi` — Gene-level tabix index
- `hg38_exons.bed.gz` + `.tbi` — Exon-level tabix index (MANE Select only)
- `hg38_gene_aliases.json.gz` — Gene alias mapping

After building, reinstall the package to include the index files:

```bash
pip install .
```

---

## Full CLI Reference

```
covsnap [-h] [--version] [--bed BED] [--exons] [--reference FASTA]
             [--no-index] [--engine {auto,mosdepth,samtools}]
             [--threads N] [-o FILE] [--emit-lowcov]
             [--lowcov-threshold N] [--lowcov-min-len N]
             [--max-targets N] [--max-total-bp N] [--max-bed-bytes BYTES]
             [--on-large-bed {error,warn_and_clip,warn_and_sample}]
             [--large-bed-seed N] [--pct-thresholds LIST]
             [--pass-pct-ge-20 F] [--pass-max-pct-zero F]
             [--dropout-pct-zero F] [--uneven-cv F]
             [--exon-pct-ge-20 F] [--exon-max-pct-zero F]
             [-v] [--quiet]
             alignment [target]
```

### Positional arguments

| Argument | Description |
|---|---|
| `alignment` | Path to BAM or CRAM file |
| `target` | Gene symbol, comma-separated gene list, or genomic region. Mutually exclusive with `--bed` |

### Commonly used options

| Flag | Description | Default |
|---|---|---|
| `--bed BED` | BED file of target intervals | — |
| `--exons` | Enable exon-level statistics (gene mode only) | off |
| `--reference FASTA` | Reference FASTA for CRAM decoding | — |
| `--engine` | Depth engine: `auto`, `mosdepth`, `samtools` | `auto` |
| `--threads N` | Parallel workers for samtools / threads for mosdepth | 4 |
| `-o FILE` / `--output FILE` | HTML report output path | `covsnap.report.html` |
| `--emit-lowcov` | Include low-coverage blocks in the report | off |
| `-v` / `--verbose` | Increase verbosity (repeatable) | — |
| `--quiet` | Suppress non-error output | off |

---

## Coordinate Convention

All output coordinates use **0-based half-open** intervals, consistent with BED format:

```
# A 100 bp region starting at position 1000
contig    start    end      length_bp
chr17     999      1099     100
```

User-facing region input accepts **1-based inclusive** coordinates (e.g. `chr17:1000-1099`), which are internally converted.

---

## Examples

### Gene mode with custom output

```bash
covsnap sample.bam BRCA1 -o results/brca1.html
```

### Multiple genes with exon breakdown

```bash
covsnap sample.bam BRCA1,TP53,ETFDH --exons -o panel_report.html
```

### Multi-gene panel via BED

```bash
covsnap sample.bam --bed panel_targets.bed -o panel_report.html
```

### Exon-level analysis with low-coverage output

```bash
covsnap sample.bam BRCA1 --exons --emit-lowcov --lowcov-threshold 20
```

### Strict BED guardrails

```bash
covsnap sample.bam --bed wes_targets.bed \
    --on-large-bed error \
    --max-targets 500 \
    --max-total-bp 10000000
```

### Using samtools explicitly with more threads

```bash
covsnap sample.bam TP53 --engine samtools --threads 8
```

---

## Exit Codes

| Code | Meaning |
|---|---|
| 0 | Success |
| 1 | Invalid arguments or input validation failure |
| 2 | Engine error (samtools/mosdepth failure) |
| 3 | Unknown gene name (with fuzzy suggestions printed to stderr) |
| 4 | BED guardrail limits exceeded (when `--on-large-bed error`) |
| 5 | CRAM reference not provided (missing `--reference` and no `REF_PATH`/`REF_CACHE`) |

---

## Running Tests

```bash
pip install ".[test]"
pytest
```

The test suite uses synthetic BAM files generated on the fly (no real sequencing data needed). Tests requiring the full GENCODE index or mosdepth are automatically skipped if unavailable.

---

## Project Structure

```
covsnap/
├── src/covsnap/
│   ├── __init__.py          # Version, build, annotation constants
│   ├── cli.py               # CLI entry point and orchestration
│   ├── annotation.py        # Gene lookup, contig detection, region parsing
│   ├── bed.py               # Streaming BED parser with guardrails
│   ├── metrics.py           # TargetAccumulator (Welford + histogram)
│   ├── engines.py           # samtools / mosdepth depth computation
│   ├── gui.py               # Tkinter graphical interface
│   ├── html_report.py       # Self-contained interactive HTML report
│   ├── report.py            # Classification heuristics
│   └── data/                # Gene/exon tabix indexes (GENCODE v44)
├── tests/                   # Comprehensive test suite
├── scripts/
│   └── build_gene_index.py  # GENCODE GTF → tabix index builder
├── recipes/conda/           # Bioconda-compatible recipe
└── pyproject.toml
```

---

## License

MIT License. See [LICENSE](LICENSE) for details.
