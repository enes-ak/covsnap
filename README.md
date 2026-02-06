# covsnap

**Coverage inspector for targeted sequencing QC (hg38)**

covsnap computes per-target (and optionally per-exon) depth-of-coverage metrics from BAM/CRAM files aligned to the human reference genome **hg38**. It produces a machine-readable TSV alongside a human-readable Markdown report with automated PASS/FAIL classification heuristics — designed for clinical and research sequencing QC workflows.

---

## Key Features

- **Gene-aware analysis** — Look up genes by symbol (e.g. `BRCA1`, `TP53`). Ships with a built-in dictionary of ~60 clinically relevant genes and an optional full GENCODE v44 tabix index covering 62,700+ genes.
- **Exon-level resolution** — Per-exon depth metrics via the `--exons` flag using MANE Select transcripts from GENCODE v44.
- **Region and BED modes** — Accepts genomic coordinates (`chr17:43044295-43125482`) or a BED file of arbitrary target intervals.
- **Streaming architecture** — O(1) memory per target using Welford's online algorithm for mean/variance and histogram-based exact median. No per-base depth arrays are ever held in memory.
- **Dual engine support** — Prefers [mosdepth](https://github.com/brentp/mosdepth) when available; falls back to `samtools depth`.
- **Contig auto-detection** — Transparently handles both `chr`-prefixed (UCSC) and non-prefixed (Ensembl/1000G) BAM contig naming.
- **Gene alias resolution** — Common aliases like `HER2 -> ERBB2` and `P53 -> TP53` are resolved automatically, with fuzzy suggestions for typos.
- **BED guardrails** — Configurable limits on target count, total bases, and file size to prevent accidental whole-exome/whole-genome runs.
- **Classification heuristics** — Automated PASS / DROP_OUT / UNEVEN / LOW_EXON / LOW_COVERAGE calls with tunable thresholds.
- **Multiple output formats** — Raw TSV (25 columns), JSON, exon TSV, low-coverage BED, and interpreted Markdown report.

---

## Installation

### From source (pip)

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

### Gene mode

Analyze coverage for a gene by name:

```bash
covsnap sample.bam BRCA1
```

This produces two files in the current directory:
- `covsnap.raw.tsv` — machine-readable metrics
- `covsnap.report.md` — human-readable interpreted report

### Region mode

Specify an explicit genomic region (1-based inclusive coordinates):

```bash
covsnap sample.bam chr17:43044295-43125482
```

### BED mode

Use a BED file of target intervals:

```bash
covsnap sample.bam --bed targets.bed
```

### CRAM files

```bash
covsnap sample.cram BRCA1 --reference hg38.fa
```

### With exon-level detail

```bash
covsnap sample.bam BRCA1 --exons
```

This additionally writes `covsnap.exons.tsv` with per-exon metrics for each MANE Select exon.

---

## Output Files

### Raw TSV (`covsnap.raw.tsv`)

Tab-separated file with a metadata header line and 25 columns:

| Column | Description |
|---|---|
| `target_id` | Gene symbol or BED interval label |
| `contig` | Chromosome (always hg38 style) |
| `start` | 0-based start coordinate |
| `end` | Half-open end coordinate |
| `length_bp` | Target length in base pairs |
| `mean_depth` | Mean read depth |
| `median_depth` | Exact median depth (histogram-based) |
| `min_depth` | Minimum depth at any position |
| `max_depth` | Maximum depth at any position |
| `stdev_depth` | Standard deviation (Welford's algorithm) |
| `pct_zero` | Percentage of bases with zero coverage |
| `pct_ge_1` | Percentage of bases with depth >= 1 |
| `pct_ge_5` | Percentage of bases with depth >= 5 |
| `pct_ge_10` | Percentage of bases with depth >= 10 |
| `pct_ge_20` | Percentage of bases with depth >= 20 |
| `pct_ge_30` | Percentage of bases with depth >= 30 |
| `pct_ge_50` | Percentage of bases with depth >= 50 |
| `pct_ge_100` | Percentage of bases with depth >= 100 |
| `n_lowcov_blocks` | Number of contiguous low-coverage blocks |
| `lowcov_total_bp` | Total bases in low-coverage blocks |
| `engine_used` | `samtools` or `mosdepth` |
| `bam_path` | Path to the input alignment file |
| `sample_name` | Sample name from BAM `@RG SM` tag |
| `build` | Always `hg38` |
| `annotation_version` | Always `gencode_v44` |

The threshold columns (`pct_ge_X`) are customizable with `--pct-thresholds`.

### Exon TSV (`covsnap.exons.tsv`)

Written when `--exons` is used. Contains 12 columns per exon:

| Column | Description |
|---|---|
| `target_id` | Parent gene symbol |
| `exon_id` | GENCODE exon stable ID (e.g. `ENSE00003896691.1`) |
| `exon_number` | Exon number within the MANE Select transcript |
| `contig` | Chromosome |
| `start` | 0-based start |
| `end` | Half-open end |
| `length_bp` | Exon length |
| `mean_depth` | Mean depth across the exon |
| `median_depth` | Exact median depth |
| `pct_zero` | Percentage of zero-coverage bases |
| `pct_ge_20` | Percentage >= 20x |
| `pct_ge_30` | Percentage >= 30x |

### JSON output (`--json-out`)

Optional structured JSON with all metrics, suitable for programmatic consumption and database loading.

### Low-coverage BED (`--emit-lowcov`)

BED file of contiguous low-coverage blocks. Controlled by:
- `--lowcov-threshold` (default: 10) — depth below which a position is "low-coverage"
- `--lowcov-min-len` (default: 50) — minimum block length to report

### Markdown report (`covsnap.report.md`)

Human-readable report including:
- Run metadata (sample, build, engine, date)
- Per-target summary table with classification
- Detailed per-target metrics
- Exon-level breakdown (when `--exons` is used)
- Low-coverage block listing
- Classification heuristics reference

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
             [--threads N] [--raw-out FILE] [--report-out FILE]
             [--json-out FILE] [--exon-out FILE] [--emit-lowcov]
             [--lowcov-bed FILE] [--lowcov-threshold N] [--lowcov-min-len N]
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
| `target` | Gene symbol or genomic region. Mutually exclusive with `--bed` |

### Commonly used options

| Flag | Description | Default |
|---|---|---|
| `--bed BED` | BED file of target intervals | — |
| `--exons` | Enable exon-level statistics (gene mode only) | off |
| `--reference FASTA` | Reference FASTA for CRAM decoding | — |
| `--engine` | Depth engine: `auto`, `mosdepth`, `samtools` | `auto` |
| `--threads N` | Threads for mosdepth | 4 |
| `--raw-out FILE` | Raw TSV output path | `covsnap.raw.tsv` |
| `--report-out FILE` | Report markdown path | `covsnap.report.md` |
| `--json-out FILE` | JSON output path | — |
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

### Basic gene analysis with custom output paths

```bash
covsnap sample.bam BRCA1 \
    --raw-out results/brca1.tsv \
    --report-out results/brca1.report.md
```

### Multi-gene panel via BED

```bash
covsnap sample.bam --bed panel_targets.bed \
    --raw-out panel_results.tsv \
    --report-out panel_report.md \
    --json-out panel_results.json
```

### Exon-level analysis with low-coverage output

```bash
covsnap sample.bam BRCA1 \
    --exons \
    --exon-out brca1_exons.tsv \
    --emit-lowcov \
    --lowcov-bed brca1_lowcov.bed \
    --lowcov-threshold 20
```

### Strict BED guardrails

```bash
covsnap sample.bam --bed wes_targets.bed \
    --on-large-bed error \
    --max-targets 500 \
    --max-total-bp 10000000
```

### Using samtools explicitly

```bash
covsnap sample.bam TP53 --engine samtools
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
│   ├── output.py            # TSV, JSON, BED writers
│   ├── report.py            # Classification heuristics + Markdown report
│   └── data/                # Gene/exon tabix indexes (GENCODE v44)
├── tests/                   # Comprehensive test suite (~90 tests)
├── scripts/
│   └── build_gene_index.py  # GENCODE GTF → tabix index builder
├── examples/                # Sample output files
├── recipes/conda/           # Bioconda-compatible recipe
└── pyproject.toml
```

---

## License

MIT License. See [LICENSE](LICENSE) for details.
