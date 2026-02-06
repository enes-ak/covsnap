"""Shared pytest fixtures for covsnap tests."""

import os

import pysam
import pytest


@pytest.fixture(scope="session")
def synthetic_bam(tmp_path_factory):
    """Create a minimal synthetic BAM with reads over BRCA1 locus.

    Generates 200 reads (150bp each, spaced 400bp apart) starting at
    chr17:43044294, covering part of the BRCA1 gene region. Sorted and
    indexed. Session-scoped so it's created only once.

    Coverage pattern:
    - Positions 43044294 to ~43124294: sparse coverage (~0-1x per position)
    - This gives us a realistic "low-coverage" scenario for testing
      classification logic.
    """
    tmp_dir = tmp_path_factory.mktemp("data")
    unsorted = str(tmp_dir / "unsorted.bam")
    sorted_bam = str(tmp_dir / "test.bam")

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr17", "LN": 83257441}],
        "RG": [{"ID": "rg1", "SM": "SYNTH_SAMPLE", "PL": "ILLUMINA"}],
    }

    with pysam.AlignmentFile(unsorted, "wb", header=header) as outf:
        for i in range(200):
            a = pysam.AlignedSegment()
            a.query_name = f"read_{i:04d}"
            a.query_sequence = "A" * 150
            a.flag = 0
            a.reference_id = 0
            a.reference_start = 43044294 + (i * 400)
            a.mapping_quality = 60
            a.cigar = [(0, 150)]  # 150M
            a.query_qualities = pysam.qualitystring_to_array("I" * 150)
            a.set_tag("RG", "rg1")
            outf.write(a)

    pysam.sort("-o", sorted_bam, unsorted)
    pysam.index(sorted_bam)
    os.remove(unsorted)

    return sorted_bam


@pytest.fixture(scope="session")
def synthetic_bam_nochr(tmp_path_factory):
    """Create a synthetic BAM using non-chr contig style (e.g., '17')."""
    tmp_dir = tmp_path_factory.mktemp("data_nochr")
    unsorted = str(tmp_dir / "unsorted.bam")
    sorted_bam = str(tmp_dir / "test_nochr.bam")

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "17", "LN": 83257441}],
        "RG": [{"ID": "rg1", "SM": "NOCHR_SAMPLE", "PL": "ILLUMINA"}],
    }

    with pysam.AlignmentFile(unsorted, "wb", header=header) as outf:
        for i in range(50):
            a = pysam.AlignedSegment()
            a.query_name = f"read_{i:04d}"
            a.query_sequence = "A" * 150
            a.flag = 0
            a.reference_id = 0
            a.reference_start = 43044294 + (i * 400)
            a.mapping_quality = 60
            a.cigar = [(0, 150)]
            a.query_qualities = pysam.qualitystring_to_array("I" * 150)
            a.set_tag("RG", "rg1")
            outf.write(a)

    pysam.sort("-o", sorted_bam, unsorted)
    pysam.index(sorted_bam)
    os.remove(unsorted)

    return sorted_bam


@pytest.fixture
def tmp_output_dir(tmp_path):
    """Provide a temporary directory for output files."""
    return tmp_path


@pytest.fixture
def sample_bed(tmp_path):
    """Create a minimal BED file with two BRCA1 sub-regions."""
    bed = tmp_path / "targets.bed"
    bed.write_text(
        "chr17\t43044294\t43050000\tregion_1\n"
        "chr17\t43090000\t43095000\tregion_2\n"
    )
    return str(bed)


@pytest.fixture
def large_bed(tmp_path):
    """Create a BED file that exceeds default --max-targets (2000)."""
    bed = tmp_path / "large_targets.bed"
    lines = []
    for i in range(3000):
        start = 1000000 + i * 1000
        end = start + 500
        lines.append(f"chr1\t{start}\t{end}\ttarget_{i:04d}")
    bed.write_text("\n".join(lines) + "\n")
    return str(bed)
