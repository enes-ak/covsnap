# Contributing to covsnap

Thanks for your interest in contributing. This document describes how to
report issues, seek support, and contribute code.

## Reporting issues and seeking support

Please use the [GitHub issue tracker](https://github.com/enes-ak/covsnap/issues)
for all of the following:

- **Bug reports** — include `covsnap --version`, your Python version and OS,
  the full command you ran, and the error output. A minimal reproducer (a small
  BAM file or a public sample name) makes the report actionable.
- **Feature requests** — describe the use case and, if possible, an example of
  the command you wish existed.
- **Usage questions** — feel free to open an issue; there is no separate
  support channel.

Search existing open and closed issues before opening a new one.

## Contributing code

### Development setup

```bash
git clone https://github.com/enes-ak/covsnap
cd covsnap
pip install -e ".[dev]"
```

This installs covsnap in editable mode along with `ruff`, `pytest`, and
`pytest-timeout`.

Runtime dependencies not installable via pip:

- `samtools` (required) — available via conda or your system package manager.
- `mosdepth` (optional, preferred engine) — `conda install -c bioconda mosdepth`.

### Running tests

```bash
pytest
```

The suite uses synthetic BAM files generated with `pysam`; no real sequencing
data is needed. Tests that require `samtools` or `mosdepth` are skipped if
those tools are not on `$PATH`. Expected runtime is under two minutes.

### Code style

covsnap uses `ruff` for linting and formatting:

```bash
ruff check src/ tests/
ruff format src/ tests/
```

CI runs both checks; pull requests that fail linting will not be merged.

### Pull request process

1. Fork the repository and create a topic branch from `master`.
2. Make your changes, adding tests that cover new behavior.
3. Run `pytest` and `ruff` locally — CI runs both across Python 3.9–3.12.
4. Open a pull request against `master` describing the change and its
   motivation. Reference any related issues.
5. A maintainer will review; respond to feedback by pushing additional
   commits to the same branch.

For larger changes (new features, architectural shifts, new dependencies),
please open an issue first to discuss the approach before investing
significant work.

## Code of Conduct

Participation in this project is governed by the
[Code of Conduct](CODE_OF_CONDUCT.md). By participating, you agree to abide
by its terms.
