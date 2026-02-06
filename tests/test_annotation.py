"""Tests for the gene annotation lookup module."""

import pytest

from covsnap.annotation import (
    _has_full_index,
    detect_contig_style,
    is_region_string,
    lookup_exons,
    lookup_gene,
    parse_region,
    suggest_genes,
    translate_contig,
)


class TestGeneLookup:
    """Test gene name → region resolution from built-in index."""

    def test_known_gene_brca1(self):
        result = lookup_gene("BRCA1")
        assert result is not None
        assert result["contig"] == "chr17"
        assert result["start"] == 43044294  # 0-based
        # end varies: built-in dict says 43125482; full GENCODE v44 index says 43170245
        assert result["end"] >= 43125482
        assert result["strand"] == "-"

    def test_known_gene_brca2(self):
        result = lookup_gene("BRCA2")
        assert result is not None
        assert result["contig"] == "chr13"
        assert result["gene_type"] == "protein_coding"

    def test_known_gene_tp53(self):
        result = lookup_gene("TP53")
        assert result is not None
        assert result["contig"] == "chr17"

    def test_case_insensitive(self):
        for name in ["brca1", "Brca1", "BRCA1", "bRcA1"]:
            result = lookup_gene(name)
            assert result is not None, f"Failed for: {name}"
            assert result["gene_name"] == "BRCA1"

    def test_alias_resolution(self):
        """Common gene aliases resolve correctly."""
        # HER2 → ERBB2
        result = lookup_gene("HER2")
        assert result is not None
        assert result["gene_name"] == "ERBB2"

        # P53 → TP53
        result = lookup_gene("P53")
        assert result is not None
        assert result["gene_name"] == "TP53"

    def test_unknown_gene_returns_none(self):
        result = lookup_gene("NOT_A_REAL_GENE_12345")
        assert result is None

    def test_fuzzy_suggestions(self):
        suggestions = suggest_genes("BRAC1")
        assert "BRCA1" in suggestions

    def test_fuzzy_suggestions_brca(self):
        suggestions = suggest_genes("BRCA")
        assert any("BRCA" in s for s in suggestions)


class TestContigDetection:
    """Test BAM contig style auto-detection."""

    def test_chr_style_detected(self, synthetic_bam):
        style = detect_contig_style(synthetic_bam)
        assert style == "chr"

    def test_nochr_style_detected(self, synthetic_bam_nochr):
        style = detect_contig_style(synthetic_bam_nochr)
        assert style == "nochr"


class TestContigTranslation:
    """Test contig name translation between chr/nochr styles."""

    def test_chr_to_nochr(self):
        assert translate_contig("chr17", "nochr") == "17"
        assert translate_contig("chrX", "nochr") == "X"
        assert translate_contig("chrM", "nochr") == "MT"

    def test_nochr_to_chr(self):
        assert translate_contig("17", "chr") == "chr17"
        assert translate_contig("X", "chr") == "chrX"
        assert translate_contig("MT", "chr") == "chrM"

    def test_already_correct_style(self):
        assert translate_contig("chr17", "chr") == "chr17"
        assert translate_contig("17", "nochr") == "17"


class TestRegionParsing:
    """Test region string parsing and validation."""

    def test_valid_region(self):
        contig, start, end = parse_region("chr17:43044295-43125482")
        assert contig == "chr17"
        assert start == 43044294  # converted to 0-based
        assert end == 43125482  # half-open stays same

    def test_region_with_commas(self):
        contig, start, end = parse_region("chr17:43,044,295-43,125,482")
        assert contig == "chr17"
        assert start == 43044294
        assert end == 43125482

    def test_invalid_region_raises(self):
        with pytest.raises(ValueError):
            parse_region("chr17:abc-def")

    def test_start_gt_end_raises(self):
        with pytest.raises(ValueError, match="start.*less than.*end"):
            parse_region("chr17:50000000-40000000")

    def test_missing_colon(self):
        with pytest.raises(ValueError):
            parse_region("chr17_43044295_43125482")


class TestIsRegionString:
    """Test region string detection."""

    def test_region_detected(self):
        assert is_region_string("chr17:43044295-43125482") is True
        assert is_region_string("chrX:100-200") is True

    def test_gene_not_detected(self):
        assert is_region_string("BRCA1") is False
        assert is_region_string("TP53") is False

    def test_malformed_not_detected(self):
        assert is_region_string("chr17:abc-def") is False


# ---------------------------------------------------------------------------
# Full-index tests (require the built tabix gene + exon index)
# ---------------------------------------------------------------------------

_skip_no_index = pytest.mark.skipif(
    not _has_full_index(),
    reason="Full hg38 gene index not built; run scripts/build_gene_index.py first",
)


@_skip_no_index
class TestFullIndexLookup:
    """Tests that exercise the tabix-indexed gene file."""

    def test_full_index_gene_count(self):
        """Tabix index has tens of thousands of genes."""
        from covsnap.annotation import _load_tabix_index

        tbx = _load_tabix_index()
        assert tbx is not None
        assert len(tbx.contigs) >= 24  # at least chr1-22 + X + Y

    def test_obscure_gene_found(self):
        """A gene NOT in the built-in dict but present in the full index."""
        result = lookup_gene("PCSK9")
        assert result is not None
        assert result["contig"] == "chr1"
        assert result["gene_type"] == "protein_coding"

    def test_lncrna_gene(self):
        """Non-coding gene lookups work."""
        result = lookup_gene("XIST")
        assert result is not None
        assert result["contig"] == "chrX"

    def test_gene_on_chrY(self):
        result = lookup_gene("SRY")
        assert result is not None
        assert result["contig"] == "chrY"
        assert result["strand"] == "-"

    def test_mitochondrial_gene(self):
        result = lookup_gene("MT-CO1")
        assert result is not None
        assert result["contig"] == "chrM"


@_skip_no_index
class TestExonLookup:
    """Test exon-level annotation from the packaged exon index."""

    def test_brca1_has_exons(self):
        exons = lookup_exons("BRCA1")
        assert len(exons) > 0
        # BRCA1 MANE Select transcript has 23 exons
        assert len(exons) == 23

    def test_brca1_exon_fields(self):
        exons = lookup_exons("BRCA1")
        for e in exons:
            assert "contig" in e
            assert "start" in e
            assert "end" in e
            assert "exon_id" in e
            assert "gene_name" in e
            assert e["gene_name"] == "BRCA1"
            assert e["contig"] == "chr17"
            assert e["start"] < e["end"]
            assert e["transcript_id"].startswith("ENST")

    def test_brca1_exons_sorted_by_position(self):
        exons = lookup_exons("BRCA1")
        starts = [e["start"] for e in exons]
        assert starts == sorted(starts)

    def test_brca1_exons_within_gene_bounds(self):
        gene = lookup_gene("BRCA1")
        exons = lookup_exons("BRCA1")
        for e in exons:
            assert e["start"] >= gene["start"]
            assert e["end"] <= gene["end"]

    def test_brca2_exons(self):
        exons = lookup_exons("BRCA2")
        assert len(exons) > 0
        assert all(e["contig"] == "chr13" for e in exons)
        # BRCA2 MANE Select has 27 exons
        assert len(exons) == 27

    def test_tp53_exons(self):
        exons = lookup_exons("TP53")
        assert len(exons) > 0
        assert all(e["contig"] == "chr17" for e in exons)

    def test_egfr_exons(self):
        exons = lookup_exons("EGFR")
        assert len(exons) > 0
        assert all(e["contig"] == "chr7" for e in exons)

    def test_unknown_gene_returns_empty(self):
        exons = lookup_exons("NOT_A_REAL_GENE_XYZ")
        assert exons == []

    def test_exon_ids_are_unique(self):
        exons = lookup_exons("BRCA1")
        ids = [e["exon_id"] for e in exons]
        assert len(ids) == len(set(ids))

    def test_exon_numbers_sequential(self):
        """Exon numbers should be present and reasonable."""
        exons = lookup_exons("BRCA1")
        numbers = [e["exon_number"] for e in exons]
        assert all(isinstance(n, int) and n > 0 for n in numbers)
        assert len(set(numbers)) == len(numbers)  # unique
