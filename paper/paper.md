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

Coverage QC for targeted panels typically involves combining `samtools depth` or `mosdepth` [@Pedersen2018] with ad hoc scripts to aggregate per-base depth into per-target statistics, compute coverage breadth at clinically relevant thresholds (e.g., percentage of bases with depth $\geq$ 20x), and flag problematic targets. This multi-step process is error-prone, difficult to standardize, and inaccessible to laboratory personnel without bioinformatics training.

covsnap addresses these issues by providing gene-aware analysis without external annotation files, automated classification heuristics with tunable thresholds, interactive HTML reports suitable for direct review or LIMS archival, and a cross-platform graphical interface for users without command-line experience. The target audience includes clinical laboratory scientists validating panel runs, bioinformaticians building sequencing QC pipelines, and researchers performing targeted resequencing studies.

# State of the field

Several tools address aspects of sequencing coverage analysis but leave gaps that covsnap fills.

`mosdepth` [@Pedersen2018] is a fast, widely used tool for computing depth statistics from BAM/CRAM files. However, it produces raw per-base or per-region depth output that requires downstream aggregation and interpretation. It has no built-in gene awareness, classification logic, or reporting capability.

Picard's `CollectHsMetrics` [@Picard] computes hybrid-selection metrics including on-target percentage, fold enrichment, and coverage uniformity. It requires pre-prepared BED and interval list files, produces tabular output without automated pass/fail interpretation, and does not support gene-symbol-based queries.

`Qualimap` [@GarciAlcalde2012] provides comprehensive BAM QC including coverage histograms and GC-bias plots. While useful for general alignment quality assessment, it lacks gene-aware targeted panel analysis and does not provide per-target or per-exon classification heuristics.

`sambamba depth` offers fast multi-threaded depth computation but, like mosdepth, provides raw metrics without classification or reporting.

covsnap differentiates itself by bundling annotation lookup (GENCODE v44, 62,700+ genes), streaming depth aggregation, classification, and interactive HTML reporting into a single tool that requires only a gene symbol and BAM file as input.

# Software design

covsnap's architecture is organized around four loosely coupled concerns, each implemented as a separate module with well-defined interfaces.

**Annotation** (`annotation.py`). A bundled gene dictionary of approximately 60 clinically relevant genes enables immediate use. An optional full GENCODE v44 tabix index extends coverage to 62,700+ genes and 201,000+ MANE Select exons. Gene alias resolution (e.g., HER2 $\rightarrow$ ERBB2) and fuzzy matching for typos are included. Contig naming style (UCSC `chr`-prefixed vs. Ensembl numeric) is auto-detected from BAM headers.

**Depth computation** (`engines.py`). covsnap delegates per-base depth calculation to `mosdepth` (preferred) or `samtools depth`, selected automatically or by the user. Depth streams are consumed in a single pass using Welford's online algorithm [@Welford1962] for mean and variance and a histogram-based exact median, achieving O(1) memory per target regardless of region size. This streaming design was chosen over in-memory depth arrays to support large panels without excessive memory consumption. Multi-target and multi-gene analyses run in parallel using Python's `ThreadPoolExecutor`.

**Classification** (`report.py`). Ordered heuristics assign each target a status (PASS, DROP_OUT, UNEVEN, LOW_EXON, or LOW_COVERAGE) based on coverage breadth ($\%\geq$20x), zero-coverage fraction, coefficient of variation, and per-exon metrics. All thresholds are user-configurable via CLI flags, enabling adaptation to different laboratory protocols and quality standards.

**Reporting** (`html_report.py`). The HTML report is generated as a single file with embedded CSS and JavaScript — no external dependencies. It includes summary cards, per-exon coverage bar charts with HSL-interpolated color gradients, collapsible detail sections, and a glossary of metrics and classification terms. The self-contained design ensures reports can be attached to emails, archived in a LIMS, or opened on any device with a web browser.

**BED guardrails** (`bed.py`) prevent accidental whole-exome or whole-genome analysis by enforcing configurable limits on target count, total bases, and file size, with options to clip or sample when limits are exceeded.

The tool is packaged with `setuptools`, distributed via PyPI and Bioconda, and registered on bio.tools. A comprehensive test suite of 101 tests uses synthetic BAM files generated with pysam, requiring no real sequencing data.

# Research impact statement

covsnap was developed to address a recurring need encountered during targeted sequencing QC in both clinical and research settings. Since its initial release in February 2026, it has been downloaded over 85 times from Bioconda and over 100 times from PyPI (as of April 2026). The tool is registered on bio.tools (https://bio.tools/covsnap), indexed in Bioconda, and available on PyPI, making it discoverable through standard bioinformatics package channels. Its single-command design and built-in gene index lower the barrier to adoption for laboratories that lack dedicated bioinformatics support, and its interactive HTML reports provide an auditable QC artifact suitable for clinical sequencing workflows.

# AI usage disclosure

Generative AI tools (Claude, Anthropic) were used during the development of covsnap to assist with code implementation, test writing, and documentation drafting, including this manuscript. All AI-generated code was reviewed, tested (101 automated tests), and validated by the author. The software architecture, design decisions, and scientific content were directed by the author.

# Acknowledgements

The author thanks the Bioconda and PyPI communities for providing packaging infrastructure, and the developers of mosdepth, samtools, and pysam whose tools covsnap builds upon.

# References
