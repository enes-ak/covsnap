---
title: 'covsnap: A single-command coverage inspector for targeted sequencing quality control'
tags:
  - Python
  - bioinformatics
  - coverage analysis
  - targeted sequencing
  - quality control
  - next-generation sequencing
authors:
  - name: Enes Ak
    orcid: 0000-0001-5931-730X
    affiliation: 1
affiliations:
  - name: Independent developer
    index: 1
date: 13 April 2026
bibliography: paper.bib
---

# Summary

Targeted sequencing panels are widely used in clinical diagnostics and research to interrogate specific genomic regions at high depth. Ensuring adequate and uniform coverage across all targeted genes and exons is a critical quality control (QC) step before downstream variant calling. However, existing tools for coverage assessment typically require multi-step workflows: computing raw depth statistics with one tool, parsing and aggregating the output with custom scripts, and manually interpreting the results against laboratory-defined thresholds. This fragmented process creates friction for both bioinformaticians building pipelines and bench scientists reviewing sequencing runs.

**covsnap** is a Python command-line tool and graphical application that compresses this workflow into a single command. Given a BAM or CRAM file and one or more gene symbols, a genomic region, or a BED file, covsnap computes per-target and per-exon depth-of-coverage metrics and produces a self-contained interactive HTML report with automated PASS/FAIL classifications.

# Statement of need

Coverage QC for targeted panels typically involves combining `samtools depth` or `mosdepth` [@Pedersen2018] with ad hoc scripts to aggregate per-base depth into per-target statistics, compute coverage breadth at clinically relevant thresholds (e.g., percentage of bases with depth $\geq$ 20x), and flag problematic targets. Tools such as Picard's `CollectHsMetrics` [@Picard] address hybrid-selection metrics but require pre-prepared BED/interval files and produce tabular output that needs further interpretation. `Qualimap` [@GarciAlcalde2012] provides general BAM QC but lacks gene-aware targeted panel analysis with automated classification.

covsnap fills this gap by providing:

- **Gene-aware analysis without external annotation files.** Users specify gene symbols (e.g., `BRCA1,TP53,ETFDH`), and covsnap resolves coordinates from a bundled GENCODE v44 index covering 62,700+ genes. No internet access or GTF downloads are required.
- **Automated classification heuristics.** Each target receives a status (PASS, DROP_OUT, UNEVEN, LOW_EXON, or LOW_COVERAGE) based on tunable thresholds, replacing manual interpretation.
- **Interactive HTML reports.** A single self-contained HTML file with summary cards, exon-level bar charts, expandable detail sections, and a glossary — suitable for direct review or archival in a LIMS.
- **A graphical interface.** Running `covsnap` with no arguments launches a cross-platform Tkinter GUI with native file pickers and progress feedback, making the tool accessible to users without command-line experience.

# Implementation

covsnap is implemented in Python and supports Python 3.9+. Its architecture is organized around four concerns: annotation lookup, depth computation, classification, and reporting.

**Annotation.** A bundled gene dictionary of approximately 60 clinically relevant genes enables immediate use. An optional full GENCODE v44 tabix index (built from a GTF with a provided script) extends coverage to 62,700+ genes and 201,000+ MANE Select exons. Gene alias resolution (e.g., HER2 $\rightarrow$ ERBB2) and fuzzy matching for typos are included. Contig naming style (UCSC `chr`-prefixed vs. Ensembl numeric) is auto-detected from BAM headers.

**Depth computation.** covsnap delegates per-base depth calculation to `mosdepth` (preferred) or `samtools depth`, selected automatically or by the user. Depth streams are consumed in a single pass using Welford's online algorithm [@Welford1962] for mean and variance and a histogram-based exact median, achieving O(1) memory per target regardless of region size. Multi-target and multi-gene analyses run in parallel using Python's `ThreadPoolExecutor`.

**Classification.** Ordered heuristics assign each target a status based on coverage breadth ($\%\geq$20x), zero-coverage fraction, coefficient of variation, and per-exon metrics. All thresholds are user-configurable.

**Reporting.** The HTML report is generated as a single file with embedded CSS and JavaScript — no external dependencies. It includes summary cards, per-exon coverage bar charts with HSL-interpolated color gradients, collapsible detail sections, and a glossary of metrics and classification terms.

**BED guardrails** prevent accidental whole-exome or whole-genome analysis by enforcing configurable limits on target count, total bases, and file size, with options to clip or sample when limits are exceeded.

The tool is packaged with `setuptools`, distributed via PyPI and Bioconda, and registered on bio.tools. A comprehensive test suite of 101 tests uses synthetic BAM files generated with pysam, requiring no real sequencing data.

# Availability

covsnap is available under the MIT license at [https://github.com/enes-ak/covsnap](https://github.com/enes-ak/covsnap). It can be installed via `pip install covsnap` or `conda install -c bioconda covsnap`.

# References
