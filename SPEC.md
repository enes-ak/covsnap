# covsnap — CLI Specification

**Version:** 0.1.0
**Build:** hg38-only
**Annotation:** GENCODE v44 (Ensembl 110), frozen at build time
**License:** MIT

---

## 1. Synopsis

```
covsnap <BAM/CRAM> <TARGET>
covsnap <BAM/CRAM> --bed <targets.bed>
```

Where `<TARGET>` is one of:
- A **gene symbol** (e.g., `BRCA1`)
- A **genomic region** string (e.g., `chr17:43044295-43125482`)

---

## 2. Coordinate Convention

| Context | Convention | Example |
|---------|-----------|---------|
| Input region string (`chr:start-end`) | **1-based inclusive** (samtools style) | `chr17:43044295-43125482` |
| Input BED file | **0-based half-open** (BED standard) | `chr17\t43044294\t43125482` |
| Packaged gene/exon index | **0-based half-open** | internal |
| Output raw TSV `start`/`end` columns | **0-based half-open** | consistent with BED |
| Output low-coverage BED | **0-based half-open** | standard BED |
| Report markdown positions | **1-based inclusive** (human-readable) | for display only |

When a user supplies a region string like `chr17:43044295-43125482`, the tool internally converts to 0-based half-open (`chr17, 43044294, 43125482`) for all computation, and stores 0-based half-open in raw output. The report displays 1-based inclusive for readability.

---

## 3. Packaged Data

The tool ships with **no internet requirement at runtime**. All annotation data is embedded in the package under `data/`:

### 3.1 Gene Index

**File:** `data/hg38_genes.tsv.gz` + `data/hg38_genes.tsv.gz.tbi`

Derived from GENCODE v44 comprehensive gene annotation (`gencode.v44.annotation.gtf.gz`), filtered to `gene` features only.

**Columns (tab-separated, bgzipped, tabix-indexed on contig/start/end):**

```
#contig  start   end     gene_name  gene_id           strand  gene_type
chr17    43044294  43125482  BRCA1  ENSG00000012048.24  -     protein_coding
```

- `contig`: chr-prefixed (canonical hg38 style)
- `start`/`end`: 0-based half-open
- `gene_name`: HGNC symbol (primary), case-preserved
- `gene_id`: Ensembl gene ID with version
- `strand`: `+` or `-`
- `gene_type`: protein_coding, lncRNA, etc.

**Generation command (build-time only, not at runtime):**

```bash
# Performed once during package build; output shipped with package
python scripts/build_gene_index.py \
  --gtf gencode.v44.annotation.gtf.gz \
  --out data/hg38_genes.tsv.gz
tabix -s1 -b2 -e3 data/hg38_genes.tsv.gz
```

A second lookup dict `data/hg38_gene_aliases.json.gz` maps common aliases and case-folded names to canonical gene_name for fuzzy matching.

### 3.2 Exon Index (optional, for `--exons` mode)

**File:** `data/hg38_exons.bed.gz` + `data/hg38_exons.bed.gz.tbi`

Derived from the same GENCODE v44, extracting `exon` features for MANE Select or canonical transcripts (one transcript per gene).

**Columns:**

```
#contig  start   end     exon_id                gene_name  exon_number  transcript_id
chr17    43044294  43045802  ENSE00003896691.1  BRCA1      24           ENST00000357654.9
```

Tabix-indexed on contig/start/end.

---

## 4. Contig Style Detection

The tool **auto-detects** the contig naming style from the BAM/CRAM header:

1. Read `@SQ` lines from the alignment header.
2. If any `SN` field starts with `chr` → **chr-style** detected.
3. If `SN` fields are bare numbers (`1`, `2`, ..., `X`, `Y`, `MT`) → **nochr-style** detected.
4. Packaged indices use chr-style.

**Behavior:**
- If BAM is chr-style: use packaged indices directly. No translation needed.
- If BAM is nochr-style: the tool prepends `chr` when querying the gene/exon index, and strips `chr` when constructing mosdepth/samtools regions. The contig translation is transparent and bidirectional.
- If detection is ambiguous (mixed or unrecognized): emit a clear actionable error:
  ```
  ERROR: Could not determine contig naming style from BAM header.
  Found contigs: ['scaffold_1', 'scaffold_2', ...]
  Expected either 'chr1'-style or '1'-style contigs for hg38.
  ```

---

## 5. Full CLI Help Text

```
usage: covsnap [-h] [--version]
                    [--bed BED]
                    [--exons]
                    [--engine {auto,mosdepth,samtools}]
                    [--reference FASTA]
                    [--no-index]
                    [--threads N]
                    [--raw-out FILE]
                    [--report-out FILE]
                    [--json-out FILE]
                    [--exon-out FILE]
                    [--emit-lowcov]
                    [--lowcov-bed FILE]
                    [--lowcov-threshold INT]
                    [--lowcov-min-len INT]
                    [--max-targets INT]
                    [--max-total-bp INT]
                    [--max-bed-bytes BYTES]
                    [--on-large-bed {error,warn_and_clip,warn_and_sample}]
                    [--large-bed-seed INT]
                    [--pct-thresholds LIST]
                    alignment [target]

covsnap — Coverage inspector for targeted sequencing QC (hg38 only).

Computes per-target and optionally per-exon depth metrics from BAM/CRAM
files, producing a machine-readable raw TSV and a human-readable
interpreted report with PASS/FAIL classifications.

Positional arguments:
  alignment             Path to BAM or CRAM file.
  target                Gene symbol (e.g. BRCA1) or genomic region
                        (e.g. chr17:43044295-43125482).
                        Mutually exclusive with --bed.

Input options:
  --bed BED             Path to a BED file of target regions.
                        Mutually exclusive with positional target.
  --exons               Enable exon-level statistics (gene mode only).
                        Requires the packaged hg38_exons.bed.gz index.
  --reference FASTA     Reference FASTA for CRAM decoding.
                        Required if input is CRAM and REF_PATH/REF_CACHE
                        environment variables are not set.
  --no-index            Do not create .bai/.csi/.crai index if missing.
                        Error if index is required but absent.

Engine options:
  --engine {auto,mosdepth,samtools}
                        Depth computation engine (default: auto).
                        auto   — use mosdepth if on PATH, else samtools.
                        mosdepth — require mosdepth; error if not found.
                        samtools — use 'samtools depth' streaming mode.
  --threads N           Threads for mosdepth (default: 4). Ignored for
                        samtools engine.

Output options:
  --raw-out FILE        Raw metrics TSV output path
                        (default: covsnap.raw.tsv).
  --report-out FILE     Interpreted report markdown output path
                        (default: covsnap.report.md).
  --json-out FILE       Optional raw metrics in JSON format.
  --exon-out FILE       Exon-level metrics TSV (default:
                        covsnap.exons.tsv). Only written if --exons.

Low-coverage output:
  --emit-lowcov         Enable writing low-coverage BED blocks.
  --lowcov-bed FILE     Output BED for low-coverage blocks
                        (default: covsnap.lowcov.bed).
  --lowcov-threshold INT
                        Depth below which a position is "low coverage"
                        (default: 10).
  --lowcov-min-len INT  Minimum contiguous length (bp) for a low-coverage
                        block to be reported (default: 50).

BED guardrail options:
  --max-targets INT     Maximum number of BED target intervals
                        (default: 2000).
  --max-total-bp INT    Maximum total target bases (default: 50000000).
  --max-bed-bytes BYTES Maximum BED file size in bytes (default: 52428800,
                        i.e. 50 MB). Accepts suffixes: K, M, G.
  --on-large-bed {error,warn_and_clip,warn_and_sample}
                        Action when BED limits are exceeded
                        (default: warn_and_clip).
                        error          — abort with a clear message.
                        warn_and_clip  — keep first N targets up to limits;
                                         warn in report.
                        warn_and_sample — deterministic random sampling of
                                          N targets; warn in report.
  --large-bed-seed INT  Random seed for warn_and_sample mode
                        (default: 42).

Coverage thresholds:
  --pct-thresholds LIST Comma-separated list of depth thresholds for
                        pct_ge_X columns (default: 1,5,10,20,30,50,100).
                        Custom thresholds replace defaults entirely.

Classification tuning:
  --pass-pct-ge-20 FLOAT
                        Minimum pct_ge_20 for PASS (default: 95.0).
  --pass-max-pct-zero FLOAT
                        Maximum pct_zero for PASS (default: 1.0).
  --dropout-pct-zero FLOAT
                        pct_zero above which DROP_OUT is called
                        (default: 5.0).
  --uneven-cv FLOAT     Coefficient of variation above which UNEVEN is
                        called (default: 1.0). Computed as stdev/mean.
  --exon-pct-ge-20 FLOAT
                        Minimum exon pct_ge_20 for LOW_EXON
                        (default: 90.0).
  --exon-max-pct-zero FLOAT
                        Maximum exon pct_zero for LOW_EXON
                        (default: 5.0).

General:
  -h, --help            Show this help message and exit.
  --version             Show program version, annotation version,
                        and build info.
  --verbose             Increase logging verbosity (can repeat: -vv).
  --quiet               Suppress all non-error output to stderr.
```

---

## 6. Argument Validation Rules

| Rule | Condition | Error message |
|------|-----------|---------------|
| Mutual exclusion | Both positional `target` and `--bed` supplied | `ERROR: Provide either a positional target (gene/region) or --bed, not both.` |
| No target | Neither positional `target` nor `--bed` supplied | `ERROR: Must provide a target: gene symbol, region string, or --bed file.` |
| CRAM without ref | Input is `.cram` and `--reference` not given and `$REF_PATH`/`$REF_CACHE` unset | `ERROR: CRAM input requires --reference <fasta> or REF_PATH/REF_CACHE env variable.` |
| BAM/CRAM exists | File does not exist or is not readable | `ERROR: Alignment file not found: <path>` |
| BED exists | `--bed` file not found | `ERROR: BED file not found: <path>` |
| Region format | Positional target matches neither gene lookup nor `contig:start-end` regex | `ERROR: '<target>' is not a known gene name and does not match region format 'contig:start-end'. Did you mean: <suggestions>?` |
| Region bounds | `start >= end` in parsed region | `ERROR: Invalid region: start (X) must be less than end (Y).` |
| Gene not found | Gene not in packaged index | `ERROR: Gene 'BRAC1' not found in hg38 gene index (GENCODE v44). Did you mean: BRCA1, BRCA2?` (case-insensitive fuzzy match) |
| Engine not found | `--engine mosdepth` but `mosdepth` not on PATH | `ERROR: mosdepth not found on PATH. Install it or use --engine samtools.` |
| Index missing | No `.bai`/`.csi`/`.crai` and `--no-index` is set | `ERROR: Index not found for <file> and --no-index prevents creation. Run 'samtools index <file>' manually.` |
| Threads | `--threads` < 1 | `ERROR: --threads must be >= 1.` |
| Thresholds | Non-integer or negative values in `--pct-thresholds` | `ERROR: --pct-thresholds must be comma-separated positive integers.` |

---

## 7. Engine Selection Logic

```
if --engine == "mosdepth":
    require mosdepth on PATH → error if absent
elif --engine == "samtools":
    use samtools depth
elif --engine == "auto":        # default
    if shutil.which("mosdepth"):
        use mosdepth
    elif shutil.which("samtools"):
        use samtools depth
    else:
        error: "Neither mosdepth nor samtools found on PATH."
```

### 7.1 mosdepth mode

- Construct a temporary BED of target region(s).
- Run: `mosdepth --by <tmp.bed> --thresholds <T1>,<T2>,... --threads N <prefix> <bam>`
- Parse `<prefix>.thresholds.bed.gz` and `<prefix>.regions.bed.gz` for per-target stats.
- For low-coverage blocks: parse `<prefix>.per-base.bed.gz` in streaming mode (never load entirely).

### 7.2 samtools depth mode

- Run: `samtools depth -a -b <tmp.bed> <bam>` (streams to stdout).
- Parse line-by-line (contig, pos, depth) — one line per base, streamed.
- Accumulate running stats per target (histogram bins, min, max, count) in O(1) RAM per target.
- Use a histogram with bins [0, 1, 5, 10, 20, 30, 50, 100, max] to compute percentile thresholds without storing per-base array.

---

## 8. Streaming & Memory Guardrails

### 8.1 Per-base depth: never loaded entirely

Both engines stream per-base data. The tool maintains per-target accumulators:

```python
class TargetAccumulator:
    count: int = 0          # positions seen
    depth_sum: int = 0      # for mean
    min_depth: int = MAX_INT
    max_depth: int = 0
    histogram: dict[int, int]  # bin → count (fixed bins)
    # For low-cov block detection:
    in_lowcov: bool = False
    lowcov_start: int = 0
    lowcov_blocks: list[tuple[int, int]]  # (start, end) pairs
    # For variance (Welford's online algorithm):
    M2: float = 0.0
    running_mean: float = 0.0
```

Memory per target: ~200 bytes + lowcov_blocks list. For 2000 targets: ~400 KB + lowcov blocks.

### 8.2 BED file streaming

BED files are **never read entirely into memory**. The file is iterated line-by-line with a streaming parser:

1. **Pre-check phase** (single pass, no storage):
   - Stat file size → compare to `--max-bed-bytes`.
   - Count lines and accumulate total bp → compare to `--max-targets` and `--max-total-bp`.
   - This pass does NOT store intervals.

2. **If limits exceeded**, apply `--on-large-bed` policy:
   - `error`: abort immediately with message.
   - `warn_and_clip`: re-iterate, yield first N intervals until both `--max-targets` and `--max-total-bp` are satisfied, then stop.
   - `warn_and_sample`: use reservoir sampling (seed = `--large-bed-seed`) during a single pass, selecting up to `--max-targets` intervals with total bp ≤ `--max-total-bp`.

3. **Sort/merge** (if needed): pipe through `bedtools sort | bedtools merge` using `subprocess.Popen` pipes — no temp file needed for the sorted output (streaming).

4. **Second pass** (execution): iterate selected intervals, dispatch to engine.

### 8.3 Limits summary

| Limit | Default | Flag |
|-------|---------|------|
| Max target intervals | 2,000 | `--max-targets` |
| Max total bp | 50,000,000 | `--max-total-bp` |
| Max BED file size | 50 MB (52,428,800 bytes) | `--max-bed-bytes` |

---

## 9. Classification Heuristics

Applied per target in the report. Evaluated in order; first match wins:

| Status | Condition | Meaning |
|--------|-----------|---------|
| `DROP_OUT` | `pct_zero > 5.0` **OR** any zero-depth block ≥ 500 bp | Complete loss of coverage in substantial portions |
| `UNEVEN` | `mean_depth > 20` **AND** `cv > 1.0` (cv = stdev / mean) | High average but erratic distribution |
| `LOW_EXON` | Exon mode enabled **AND** any exon has `pct_ge_20 < 90` or `pct_zero > 5` | Specific exons underperforming |
| `LOW_COVERAGE` | `pct_ge_20 < 95.0` **AND** `pct_zero <= 1.0` | Generally insufficient depth |
| `PASS` | `pct_ge_20 >= 95.0` **AND** `pct_zero <= 1.0` | Adequate coverage |

All thresholds are configurable via CLI flags (see Section 5).

---

## 10. Output Schemas

### 10.1 Raw Metrics TSV (`covsnap.raw.tsv`)

**Header line** (prefixed with `#`):

```
#covsnap_version=0.1.0	annotation=gencode_v44	build=hg38	engine=mosdepth	date=2026-02-06T14:30:00Z
```

**Columns:**

| # | Column | Type | Description |
|---|--------|------|-------------|
| 1 | `target_id` | str | Gene name, region string, or BED interval label |
| 2 | `contig` | str | Chromosome (matches BAM contig style) |
| 3 | `start` | int | 0-based start |
| 4 | `end` | int | 0-based half-open end |
| 5 | `length_bp` | int | `end - start` |
| 6 | `mean_depth` | float | Mean per-base depth (2 decimal places) |
| 7 | `median_depth` | float | Median per-base depth (estimated from histogram) |
| 8 | `min_depth` | int | Minimum observed depth |
| 9 | `max_depth` | int | Maximum observed depth |
| 10 | `stdev_depth` | float | Standard deviation (Welford's, 2 decimal places) |
| 11 | `pct_zero` | float | % bases with depth == 0 |
| 12 | `pct_ge_1` | float | % bases with depth >= 1 |
| 13 | `pct_ge_5` | float | % bases with depth >= 5 |
| 14 | `pct_ge_10` | float | % bases with depth >= 10 |
| 15 | `pct_ge_20` | float | % bases with depth >= 20 |
| 16 | `pct_ge_30` | float | % bases with depth >= 30 |
| 17 | `pct_ge_50` | float | % bases with depth >= 50 |
| 18 | `pct_ge_100` | float | % bases with depth >= 100 |
| 19 | `n_lowcov_blocks` | int | Number of low-coverage contiguous blocks |
| 20 | `lowcov_total_bp` | int | Total bases in low-coverage blocks |
| 21 | `engine_used` | str | `mosdepth` or `samtools` |
| 22 | `bam_path` | str | Input alignment file path |
| 23 | `sample_name` | str | From `@RG SM:` tag, or `NA` |
| 24 | `build` | str | Always `hg38` |
| 25 | `annotation_version` | str | `gencode_v44` |

All `pct_*` values are in range [0.0, 100.0] with 2 decimal places.

### 10.2 Exon Metrics TSV (`covsnap.exons.tsv`)

Written only when `--exons` is supplied and target is a gene.

| # | Column | Type | Description |
|---|--------|------|-------------|
| 1 | `target_id` | str | Parent gene name |
| 2 | `exon_id` | str | Ensembl exon ID |
| 3 | `exon_number` | int | Exon number in transcript |
| 4 | `contig` | str | Chromosome |
| 5 | `start` | int | 0-based start |
| 6 | `end` | int | 0-based half-open end |
| 7 | `length_bp` | int | `end - start` |
| 8 | `mean_depth` | float | Mean per-base depth |
| 9 | `median_depth` | float | Median per-base depth |
| 10 | `pct_zero` | float | % bases with depth == 0 |
| 11 | `pct_ge_20` | float | % bases with depth >= 20 |
| 12 | `pct_ge_30` | float | % bases with depth >= 30 |

### 10.3 Low-Coverage BED (`covsnap.lowcov.bed`)

Written only when `--emit-lowcov` is supplied. Standard BED3+ format:

```
#track name="covsnap_lowcov" description="Low-coverage blocks (depth < 10, min 50bp)"
chr17	43067600	43067750	BRCA1	mean_depth=2.3
chr17	43091400	43091520	BRCA1	mean_depth=0.0
```

### 10.4 JSON Output (`--json-out`)

Same data as raw TSV, structured as:

```json
{
  "covsnap_version": "0.1.0",
  "annotation": "gencode_v44",
  "build": "hg38",
  "engine": "mosdepth",
  "date": "2026-02-06T14:30:00Z",
  "targets": [
    {
      "target_id": "BRCA1",
      "contig": "chr17",
      "start": 43044294,
      "end": 43125482,
      ...
    }
  ]
}
```

---

## 11. Indexing Behavior

If the BAM/CRAM file lacks an index (`.bai`, `.csi`, or `.crai`):

1. **Default:** Automatically create the index using `samtools index <file>`.
   - Print to stderr: `INFO: Creating index for <file>...`
2. **`--no-index`:** Do NOT create an index. If the engine requires one (mosdepth always does, samtools depth with `-b` does), emit:
   ```
   ERROR: Index required but not found for '<file>', and --no-index is set.
   Create it manually: samtools index <file>
   ```

---

## 12. Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success (all targets processed) |
| 1 | General error (invalid arguments, file not found, etc.) |
| 2 | Engine not found or engine error |
| 3 | Gene/region not found in annotation |
| 4 | BED guardrail limit exceeded with `--on-large-bed error` |
| 5 | CRAM reference missing |
| 10 | Partial success (some targets failed; report notes which) |

---

## 13. Logging

- `stderr` only (stdout is never used for messages).
- Default: `WARNING` level.
- `--verbose` (`-v`): `INFO` level.
- `-vv`: `DEBUG` level.
- `--quiet`: `ERROR` only.
- Log lines are prefixed: `[covsnap] LEVEL: message`
