"""Gene annotation lookup, region parsing, and contig style detection.

Provides gene-name → genomic-region resolution using either:
1. A tabix-indexed hg38 gene TSV shipped with the package (preferred).
2. A built-in fallback dictionary of ~60 common clinical genes.

No internet access or user-provided GTF is ever required.
"""

from __future__ import annotations

import gzip
import json
import logging
import os
import re
from difflib import get_close_matches
from pathlib import Path
from typing import Any, Optional

import pysam

from covsnap import ANNOTATION_VERSION

logger = logging.getLogger(__name__)

_DATA_DIR = Path(__file__).parent / "data"
_GENE_INDEX = _DATA_DIR / "hg38_genes.tsv.gz"
_GENE_ALIASES = _DATA_DIR / "hg38_gene_aliases.json.gz"
_EXON_INDEX = _DATA_DIR / "hg38_exons.bed.gz"

# ---------------------------------------------------------------------------
# Region string parsing
# ---------------------------------------------------------------------------

_REGION_RE = re.compile(r"^([^:]+):([0-9,]+)-([0-9,]+)$")


def parse_region(region_str: str) -> tuple[str, int, int]:
    """Parse a samtools-style region string (1-based inclusive).

    Returns (contig, start_0based, end_halfopen).
    Raises ValueError on malformed input or start >= end.
    """
    region_str = region_str.strip()
    m = _REGION_RE.match(region_str)
    if not m:
        raise ValueError(
            f"Invalid region format: '{region_str}'. "
            "Expected: contig:start-end (e.g. chr17:43044295-43125482)"
        )
    contig = m.group(1)
    start_1 = int(m.group(2).replace(",", ""))
    end_1 = int(m.group(3).replace(",", ""))
    if start_1 >= end_1:
        raise ValueError(
            f"Invalid region: start ({start_1}) must be less than end ({end_1})."
        )
    # 1-based inclusive → 0-based half-open
    return contig, start_1 - 1, end_1


# ---------------------------------------------------------------------------
# Contig style detection
# ---------------------------------------------------------------------------

_NOCHR_CONTIG_RE = re.compile(r"^(\d{1,2}|X|Y|MT?)$")


def detect_contig_style(bam_path: str) -> str:
    """Detect contig naming convention from alignment header.

    Returns ``"chr"`` or ``"nochr"``.
    Raises ``ValueError`` if the style is ambiguous.
    """
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as af:
        sq = af.header.to_dict().get("SQ", [])
        contigs = [entry["SN"] for entry in sq]

    if not contigs:
        raise ValueError("BAM/CRAM header contains no @SQ lines.")

    chr_count = sum(1 for c in contigs if c.startswith("chr"))
    nochr_count = sum(1 for c in contigs if _NOCHR_CONTIG_RE.match(c))

    if chr_count > 0 and nochr_count == 0:
        return "chr"
    if nochr_count > 0 and chr_count == 0:
        return "nochr"

    raise ValueError(
        "Could not determine contig naming style from BAM header.\n"
        f"Found contigs: {contigs[:10]}{'...' if len(contigs) > 10 else ''}\n"
        "Expected either 'chr1'-style or '1'-style contigs for hg38."
    )


def translate_contig(contig: str, target_style: str) -> str:
    """Translate a contig name between chr/nochr styles.

    The packaged index always uses chr-style.  If the BAM uses nochr-style
    we strip the ``chr`` prefix when querying the BAM, and add it when
    querying the index.
    """
    if target_style == "chr":
        if contig.startswith("chr"):
            return contig
        # nochr → chr
        if contig == "MT":
            return "chrM"
        return f"chr{contig}"
    else:
        # chr → nochr
        if not contig.startswith("chr"):
            return contig
        bare = contig[3:]
        if bare == "M":
            return "MT"
        return bare


# ---------------------------------------------------------------------------
# Built-in gene dictionary (fallback when tabix index not available)
# Coordinates: 0-based half-open, hg38, GENCODE v44
# ---------------------------------------------------------------------------

_BUILTIN_GENES: dict[str, dict[str, Any]] = {
    # --- Hereditary breast/ovarian ---
    "BRCA1": {"contig": "chr17", "start": 43044294, "end": 43125482, "strand": "-",
              "gene_id": "ENSG00000012048", "gene_type": "protein_coding"},
    "BRCA2": {"contig": "chr13", "start": 32315085, "end": 32400268, "strand": "+",
              "gene_id": "ENSG00000139618", "gene_type": "protein_coding"},
    "PALB2": {"contig": "chr16", "start": 23603159, "end": 23641310, "strand": "-",
              "gene_id": "ENSG00000083093", "gene_type": "protein_coding"},
    "BRIP1": {"contig": "chr17", "start": 61679185, "end": 61863558, "strand": "-",
              "gene_id": "ENSG00000136492", "gene_type": "protein_coding"},
    "RAD51C": {"contig": "chr17", "start": 58692572, "end": 58735610, "strand": "-",
               "gene_id": "ENSG00000108384", "gene_type": "protein_coding"},
    "RAD51D": {"contig": "chr17", "start": 35100868, "end": 35120783, "strand": "+",
               "gene_id": "ENSG00000185379", "gene_type": "protein_coding"},
    "ATM": {"contig": "chr11", "start": 108222484, "end": 108369102, "strand": "+",
            "gene_id": "ENSG00000149311", "gene_type": "protein_coding"},
    "CHEK2": {"contig": "chr22", "start": 28687743, "end": 28742422, "strand": "-",
              "gene_id": "ENSG00000183765", "gene_type": "protein_coding"},
    "CDH1": {"contig": "chr16", "start": 68737225, "end": 68835541, "strand": "+",
             "gene_id": "ENSG00000039068", "gene_type": "protein_coding"},
    # --- Tumor suppressors ---
    "TP53": {"contig": "chr17", "start": 7668401, "end": 7687550, "strand": "-",
             "gene_id": "ENSG00000141510", "gene_type": "protein_coding"},
    "PTEN": {"contig": "chr10", "start": 87862624, "end": 87971930, "strand": "+",
             "gene_id": "ENSG00000171862", "gene_type": "protein_coding"},
    "RB1": {"contig": "chr13", "start": 48303747, "end": 48481890, "strand": "+",
            "gene_id": "ENSG00000139687", "gene_type": "protein_coding"},
    "APC": {"contig": "chr5", "start": 112707497, "end": 112846239, "strand": "+",
            "gene_id": "ENSG00000134982", "gene_type": "protein_coding"},
    "VHL": {"contig": "chr3", "start": 10141727, "end": 10153668, "strand": "+",
            "gene_id": "ENSG00000134086", "gene_type": "protein_coding"},
    "NF1": {"contig": "chr17", "start": 31094927, "end": 31377677, "strand": "+",
            "gene_id": "ENSG00000196712", "gene_type": "protein_coding"},
    "NF2": {"contig": "chr22", "start": 29603556, "end": 29698600, "strand": "+",
            "gene_id": "ENSG00000186575", "gene_type": "protein_coding"},
    "WT1": {"contig": "chr11", "start": 32389063, "end": 32435260, "strand": "-",
            "gene_id": "ENSG00000184937", "gene_type": "protein_coding"},
    "CDKN2A": {"contig": "chr9", "start": 21967750, "end": 21995324, "strand": "-",
               "gene_id": "ENSG00000147889", "gene_type": "protein_coding"},
    "STK11": {"contig": "chr19", "start": 1205797, "end": 1228431, "strand": "+",
              "gene_id": "ENSG00000118046", "gene_type": "protein_coding"},
    "SMAD4": {"contig": "chr18", "start": 51028393, "end": 51085042, "strand": "+",
              "gene_id": "ENSG00000141646", "gene_type": "protein_coding"},
    "BAP1": {"contig": "chr3", "start": 52401101, "end": 52411753, "strand": "-",
             "gene_id": "ENSG00000163930", "gene_type": "protein_coding"},
    "MUTYH": {"contig": "chr1", "start": 45329162, "end": 45340399, "strand": "-",
              "gene_id": "ENSG00000132781", "gene_type": "protein_coding"},
    # --- Oncogenes ---
    "EGFR": {"contig": "chr7", "start": 55019016, "end": 55211628, "strand": "+",
             "gene_id": "ENSG00000146648", "gene_type": "protein_coding"},
    "KRAS": {"contig": "chr12", "start": 25204788, "end": 25250936, "strand": "-",
             "gene_id": "ENSG00000133703", "gene_type": "protein_coding"},
    "NRAS": {"contig": "chr1", "start": 114704463, "end": 114716894, "strand": "-",
             "gene_id": "ENSG00000213281", "gene_type": "protein_coding"},
    "BRAF": {"contig": "chr7", "start": 140719326, "end": 140924929, "strand": "-",
             "gene_id": "ENSG00000157764", "gene_type": "protein_coding"},
    "PIK3CA": {"contig": "chr3", "start": 179148113, "end": 179240093, "strand": "+",
               "gene_id": "ENSG00000121879", "gene_type": "protein_coding"},
    "ERBB2": {"contig": "chr17", "start": 39687913, "end": 39730426, "strand": "+",
              "gene_id": "ENSG00000141736", "gene_type": "protein_coding"},
    "MET": {"contig": "chr7", "start": 116672195, "end": 116798386, "strand": "+",
            "gene_id": "ENSG00000105976", "gene_type": "protein_coding"},
    "ALK": {"contig": "chr2", "start": 29192773, "end": 29921566, "strand": "-",
            "gene_id": "ENSG00000171094", "gene_type": "protein_coding"},
    "ROS1": {"contig": "chr6", "start": 117288198, "end": 117425855, "strand": "-",
             "gene_id": "ENSG00000047936", "gene_type": "protein_coding"},
    "RET": {"contig": "chr10", "start": 43077027, "end": 43130351, "strand": "+",
            "gene_id": "ENSG00000165731", "gene_type": "protein_coding"},
    "KIT": {"contig": "chr4", "start": 54657927, "end": 54740715, "strand": "+",
            "gene_id": "ENSG00000157404", "gene_type": "protein_coding"},
    "PDGFRA": {"contig": "chr4", "start": 54229096, "end": 54298245, "strand": "+",
               "gene_id": "ENSG00000134853", "gene_type": "protein_coding"},
    "AKT1": {"contig": "chr14", "start": 104769349, "end": 104795743, "strand": "-",
             "gene_id": "ENSG00000142208", "gene_type": "protein_coding"},
    "CDK4": {"contig": "chr12", "start": 57747727, "end": 57752330, "strand": "-",
             "gene_id": "ENSG00000135446", "gene_type": "protein_coding"},
    # --- Lynch syndrome / MMR ---
    "MLH1": {"contig": "chr3", "start": 36993349, "end": 37050846, "strand": "+",
             "gene_id": "ENSG00000076242", "gene_type": "protein_coding"},
    "MSH2": {"contig": "chr2", "start": 47403067, "end": 47709830, "strand": "+",
             "gene_id": "ENSG00000095002", "gene_type": "protein_coding"},
    "MSH6": {"contig": "chr2", "start": 47695530, "end": 47721031, "strand": "+",
             "gene_id": "ENSG00000116062", "gene_type": "protein_coding"},
    "PMS2": {"contig": "chr7", "start": 5970924, "end": 6009019, "strand": "-",
             "gene_id": "ENSG00000122512", "gene_type": "protein_coding"},
    "EPCAM": {"contig": "chr2", "start": 47345277, "end": 47387687, "strand": "+",
              "gene_id": "ENSG00000119888", "gene_type": "protein_coding"},
    # --- Cardiac ---
    "MYBPC3": {"contig": "chr11", "start": 47331268, "end": 47352957, "strand": "-",
               "gene_id": "ENSG00000134571", "gene_type": "protein_coding"},
    "MYH7": {"contig": "chr14", "start": 23412739, "end": 23435444, "strand": "-",
             "gene_id": "ENSG00000092054", "gene_type": "protein_coding"},
    "TNNT2": {"contig": "chr1", "start": 201328135, "end": 201346803, "strand": "-",
              "gene_id": "ENSG00000118194", "gene_type": "protein_coding"},
    "LMNA": {"contig": "chr1", "start": 156082573, "end": 156109880, "strand": "+",
             "gene_id": "ENSG00000160789", "gene_type": "protein_coding"},
    "SCN5A": {"contig": "chr3", "start": 38548061, "end": 38649670, "strand": "+",
              "gene_id": "ENSG00000183873", "gene_type": "protein_coding"},
    "KCNQ1": {"contig": "chr11", "start": 2444707, "end": 2849109, "strand": "-",
              "gene_id": "ENSG00000053918", "gene_type": "protein_coding"},
    "KCNH2": {"contig": "chr7", "start": 150642042, "end": 150675403, "strand": "-",
              "gene_id": "ENSG00000055118", "gene_type": "protein_coding"},
    # --- Other clinical ---
    "CFTR": {"contig": "chr7", "start": 117479963, "end": 117668665, "strand": "+",
             "gene_id": "ENSG00000001626", "gene_type": "protein_coding"},
    "HBB": {"contig": "chr11", "start": 5225464, "end": 5227071, "strand": "-",
            "gene_id": "ENSG00000244734", "gene_type": "protein_coding"},
    "HEXA": {"contig": "chr15", "start": 72345371, "end": 72380305, "strand": "-",
             "gene_id": "ENSG00000213614", "gene_type": "protein_coding"},
    "GBA1": {"contig": "chr1", "start": 155234452, "end": 155244699, "strand": "-",
             "gene_id": "ENSG00000177628", "gene_type": "protein_coding"},
    "FMR1": {"contig": "chrX", "start": 147911919, "end": 147951125, "strand": "-",
             "gene_id": "ENSG00000102081", "gene_type": "protein_coding"},
    "SMN1": {"contig": "chr5", "start": 70924941, "end": 70953015, "strand": "+",
             "gene_id": "ENSG00000172062", "gene_type": "protein_coding"},
    # --- Pharmacogenomics ---
    "CYP2D6": {"contig": "chr22", "start": 42126498, "end": 42130881, "strand": "-",
               "gene_id": "ENSG00000100197", "gene_type": "protein_coding"},
    "CYP2C19": {"contig": "chr10", "start": 94762680, "end": 94855547, "strand": "+",
                "gene_id": "ENSG00000165841", "gene_type": "protein_coding"},
    "DPYD": {"contig": "chr1", "start": 97543299, "end": 98386615, "strand": "-",
             "gene_id": "ENSG00000188641", "gene_type": "protein_coding"},
    "TPMT": {"contig": "chr6", "start": 18128542, "end": 18155374, "strand": "+",
             "gene_id": "ENSG00000137364", "gene_type": "protein_coding"},
    "UGT1A1": {"contig": "chr2", "start": 233757808, "end": 233773299, "strand": "-",
               "gene_id": "ENSG00000241635", "gene_type": "protein_coding"},
    # --- Common research genes ---
    "MYC": {"contig": "chr8", "start": 127735433, "end": 127742951, "strand": "+",
            "gene_id": "ENSG00000136997", "gene_type": "protein_coding"},
    "FGFR1": {"contig": "chr8", "start": 38411138, "end": 38468834, "strand": "+",
              "gene_id": "ENSG00000077782", "gene_type": "protein_coding"},
    "FGFR2": {"contig": "chr10", "start": 121478332, "end": 121598458, "strand": "-",
              "gene_id": "ENSG00000066468", "gene_type": "protein_coding"},
    "FGFR3": {"contig": "chr4", "start": 1793293, "end": 1808872, "strand": "+",
              "gene_id": "ENSG00000068078", "gene_type": "protein_coding"},
    "IDH1": {"contig": "chr2", "start": 208236227, "end": 208255071, "strand": "-",
             "gene_id": "ENSG00000138413", "gene_type": "protein_coding"},
    "IDH2": {"contig": "chr15", "start": 90083040, "end": 90102468, "strand": "-",
             "gene_id": "ENSG00000182054", "gene_type": "protein_coding"},
    "JAK2": {"contig": "chr9", "start": 4984390, "end": 5129948, "strand": "+",
             "gene_id": "ENSG00000096968", "gene_type": "protein_coding"},
}

# Build case-insensitive lookup and common aliases
_BUILTIN_LOOKUP: dict[str, str] = {}
for _name in _BUILTIN_GENES:
    _BUILTIN_LOOKUP[_name.upper()] = _name

# Common aliases → canonical name
_ALIASES: dict[str, str] = {
    "HER2": "ERBB2",
    "HER-2": "ERBB2",
    "NEU": "ERBB2",
    "P53": "TP53",
    "GBA": "GBA1",
}
for _alias, _canon in _ALIASES.items():
    _BUILTIN_LOOKUP[_alias.upper()] = _canon


# ---------------------------------------------------------------------------
# Tabix-based gene lookup (full index)
# ---------------------------------------------------------------------------

_tabix_genes: Optional[pysam.TabixFile] = None
_alias_map: Optional[dict[str, str]] = None


def _load_tabix_index() -> Optional[pysam.TabixFile]:
    """Load the tabix-indexed gene file if it exists."""
    global _tabix_genes
    if _tabix_genes is not None:
        return _tabix_genes
    if _GENE_INDEX.exists():
        try:
            _tabix_genes = pysam.TabixFile(str(_GENE_INDEX))
            return _tabix_genes
        except Exception as exc:
            logger.warning("Could not load gene index %s: %s", _GENE_INDEX, exc)
    return None


def _load_alias_map() -> dict[str, str]:
    """Load the case-folded alias map if it exists."""
    global _alias_map
    if _alias_map is not None:
        return _alias_map
    if _GENE_ALIASES.exists():
        try:
            with gzip.open(_GENE_ALIASES, "rt") as f:
                _alias_map = json.load(f)
                return _alias_map
        except Exception as exc:
            logger.warning("Could not load alias map: %s", exc)
    _alias_map = {}
    return _alias_map


def _has_full_index() -> bool:
    """Check if the full tabix gene index is available."""
    return _GENE_INDEX.exists() and Path(str(_GENE_INDEX) + ".tbi").exists()


def _search_tabix_by_name(gene_name: str) -> Optional[dict[str, Any]]:
    """Search the full tabix index for a gene name.

    Since the tabix file is indexed by genomic position (not name), we
    need a linear scan of the alias map or a name→position lookup.
    The alias map maps lower(name) → canonical name.  We then need a
    second file or a header-to-coordinates mapping.

    In practice, the alias map gives us the canonical name, and then
    we look up coordinates from a secondary dict loaded from the same
    file on first access.
    """
    tbx = _load_tabix_index()
    if tbx is None:
        return None

    # We need to scan the file for the gene name since tabix
    # indexes by genomic position.  For efficiency, we keep a name→record
    # cache built on first access.
    if not hasattr(_search_tabix_by_name, "_name_cache"):
        cache: dict[str, dict[str, Any]] = {}
        try:
            for contig in tbx.contigs:
                for row in tbx.fetch(contig):
                    fields = row.split("\t")
                    if len(fields) >= 7:
                        name = fields[3]
                        cache[name.upper()] = {
                            "contig": fields[0],
                            "start": int(fields[1]),
                            "end": int(fields[2]),
                            "gene_name": fields[3],
                            "gene_id": fields[4],
                            "strand": fields[5],
                            "gene_type": fields[6],
                        }
        except Exception as exc:
            logger.warning("Error building gene name cache: %s", exc)
        _search_tabix_by_name._name_cache = cache  # type: ignore[attr-defined]

    return _search_tabix_by_name._name_cache.get(gene_name.upper())  # type: ignore[attr-defined]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def lookup_gene(gene_name: str) -> Optional[dict[str, Any]]:
    """Resolve a gene name to its hg38 genomic coordinates.

    Tries the full tabix index first, then falls back to the built-in dict.
    Returns ``None`` if the gene is not found.

    Returned dict keys: contig, start, end, strand, gene_name, gene_id, gene_type.
    """
    # Normalise input
    query = gene_name.strip().upper()

    # Try alias resolution first
    aliases = _load_alias_map()
    canonical = aliases.get(query.lower(), _BUILTIN_LOOKUP.get(query))
    if canonical:
        query = canonical.upper()

    # Try full tabix index
    if _has_full_index():
        result = _search_tabix_by_name(gene_name)
        if result is not None:
            return result

    # Fallback to built-in
    canon_name = _BUILTIN_LOOKUP.get(query)
    if canon_name and canon_name in _BUILTIN_GENES:
        record = dict(_BUILTIN_GENES[canon_name])
        record["gene_name"] = canon_name
        if _has_full_index():
            logger.debug("Gene found in built-in dict (not in full index): %s", gene_name)
        else:
            logger.debug("Using built-in gene dict (full index not available): %s", gene_name)
        return record

    return None


def suggest_genes(gene_name: str, n: int = 5) -> list[str]:
    """Return up to *n* close matches for a possibly-misspelled gene name."""
    query = gene_name.strip().upper()

    # Collect all known names
    known: list[str] = list(_BUILTIN_GENES.keys())
    if _has_full_index():
        # Also check the tabix name cache if built
        if hasattr(_search_tabix_by_name, "_name_cache"):
            known.extend(_search_tabix_by_name._name_cache.keys())  # type: ignore[attr-defined]

    # Add aliases
    known.extend(k.upper() for k in _ALIASES)

    # Deduplicate preserving order
    seen: set[str] = set()
    unique: list[str] = []
    for k in known:
        ku = k.upper()
        if ku not in seen:
            seen.add(ku)
            unique.append(k)

    matches = get_close_matches(query, unique, n=n, cutoff=0.6)
    # Return canonical names, not aliases
    result: list[str] = []
    for m in matches:
        canon = _BUILTIN_LOOKUP.get(m.upper(), m)
        if canon not in result:
            result.append(canon)
    return result


def lookup_exons(gene_name: str) -> list[dict[str, Any]]:
    """Return exon records for a gene from the packaged exon index.

    Returns an empty list if the exon index is not available or the gene
    has no exon records.

    Each dict: contig, start, end, exon_id, gene_name, exon_number, transcript_id.
    """
    if not _EXON_INDEX.exists():
        logger.warning(
            "Exon index not found at %s. Run build_gene_index.py --build-exons first.",
            _EXON_INDEX,
        )
        return []

    gene_info = lookup_gene(gene_name)
    if gene_info is None:
        return []

    exons: list[dict[str, Any]] = []
    try:
        tbx = pysam.TabixFile(str(_EXON_INDEX))
        for row in tbx.fetch(gene_info["contig"], gene_info["start"], gene_info["end"]):
            fields = row.split("\t")
            if len(fields) >= 7 and fields[4].upper() == gene_info["gene_name"].upper():
                exons.append({
                    "contig": fields[0],
                    "start": int(fields[1]),
                    "end": int(fields[2]),
                    "exon_id": fields[3],
                    "gene_name": fields[4],
                    "exon_number": int(fields[5]) if fields[5].isdigit() else 0,
                    "transcript_id": fields[6],
                })
        tbx.close()
    except Exception as exc:
        logger.warning("Error reading exon index: %s", exc)

    exons.sort(key=lambda e: e["start"])
    return exons


def is_region_string(target: str) -> bool:
    """Return True if *target* looks like a region string (contig:start-end)."""
    return bool(_REGION_RE.match(target.strip()))


def get_sample_name(bam_path: str) -> str:
    """Extract sample name from BAM/CRAM @RG SM: tag. Returns 'NA' if absent."""
    try:
        with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as af:
            rg = af.header.to_dict().get("RG", [])
            for entry in rg:
                sm = entry.get("SM")
                if sm:
                    return sm
    except Exception:
        pass
    return "NA"
