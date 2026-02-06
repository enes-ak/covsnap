#!/bin/bash
set -euo pipefail

# Install the Python package
${PYTHON} -m pip install . --no-deps --no-build-isolation -vvv

# Ensure packaged data files are in the right place
# The gene/exon index files are included in the sdist/wheel via
# package_data in pyproject.toml / setup.cfg.
# Verify they exist after install:
DATA_DIR="${PREFIX}/lib/python${PY_VER}/site-packages/covsnap/data"
if [ ! -f "${DATA_DIR}/hg38_genes.tsv.gz" ]; then
    echo "ERROR: hg38_genes.tsv.gz not found in installed package data" >&2
    exit 1
fi
if [ ! -f "${DATA_DIR}/hg38_genes.tsv.gz.tbi" ]; then
    echo "ERROR: hg38_genes.tsv.gz.tbi not found in installed package data" >&2
    exit 1
fi

echo "covsnap installed successfully."
