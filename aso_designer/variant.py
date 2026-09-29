"""Variant/SNP awareness for ASO candidate evaluation.

Fetches variant data from dbSNP/ClinVar APIs or local VCF files,
annotates ASO candidates with overlapping variants, and computes
variant-aware penalties. Inspired by GuideX's conservation module
but redesigned for ASO-specific variant analysis (seed region
sensitivity, pathogenic flagging, MAF-based penalties).

Scientific note: Variant data is fetched from public databases
(dbSNP, ClinVar) or provided as VCF files. MAF values and clinical
significance are used for candidate ranking, not clinical diagnosis.
"""

from __future__ import annotations

from dataclasses import dataclass

from aso_designer.utils.exceptions import VariantError
from aso_designer.utils.logger import get_logger

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Variant data structures
# ---------------------------------------------------------------------------


@dataclass
class Variant:
    """A single genetic variant (SNP, indel, etc.).

    Attributes
    ----------
    chrom : str
        Chromosome (e.g., "1", "chr1", "X").
    pos : int
        Genomic position (1-based).
    ref : str
        Reference allele.
    alt : str
        Alternate allele(s), comma-separated for multi-allelic.
    rs_id : str
        dbSNP rs ID (e.g., "rs123456") or empty if unknown.
    maf : float
        Minor allele frequency (0.0–1.0). -1.0 if unknown.
    variant_type : str
        Type: "snv", "indel", "mnv", or "unknown".
    clinical_significance : str
        ClinVar significance: "pathogenic", "likely_pathogenic",
        "benign", "uncertain", or empty.
    gene : str
        Gene symbol or empty.
    """

    chrom: str
    pos: int
    ref: str
    alt: str
    rs_id: str = ""
    maf: float = -1.0
    variant_type: str = "unknown"
    clinical_significance: str = ""
    gene: str = ""

    @property
    def is_pathogenic(self) -> bool:
        """True if ClinVar significance is pathogenic or likely pathogenic."""
        return self.clinical_significance.lower() in (
            "pathogenic",
            "likely_pathogenic",
        )

    @property
    def is_snv(self) -> bool:
        """True if single nucleotide variant."""
        return self.variant_type == "snv"

    @property
    def has_maf(self) -> bool:
        """True if MAF is known (≥ 0)."""
        return self.maf >= 0.0

    def to_dict(self) -> dict:
        """Serialize to dictionary."""
        return {
            "chrom": self.chrom,
            "pos": self.pos,
            "ref": self.ref,
            "alt": self.alt,
            "rs_id": self.rs_id,
            "maf": round(self.maf, 4) if self.has_maf else -1.0,
            "variant_type": self.variant_type,
            "clinical_significance": self.clinical_significance,
            "gene": self.gene,
            "is_pathogenic": self.is_pathogenic,
        }


class VariantSet:
    """A collection of variants with position-based indexing.

    Stores variants and provides fast lookup by genomic position
    for overlap detection with ASO binding sites.

    Attributes
    ----------
    variants : list[Variant]
        All variants in the set.
    """

    def __init__(self):
        self._variants: list[Variant] = []
        # Index: (chrom, pos) → list of variant indices
        self._pos_index: dict[tuple[str, int], list[int]] = {}
        # Index: chrom → sorted list of positions
        self._chrom_pos: dict[str, list[int]] = {}

    @property
    def variants(self) -> list[Variant]:
        return self._variants

    @property
    def num_variants(self) -> int:
        return len(self._variants)

    def add_variant(self, variant: Variant) -> None:
        """Add a single variant to the set."""
        idx = len(self._variants)
        self._variants.append(variant)
        key = (variant.chrom, variant.pos)
        self._pos_index.setdefault(key, []).append(idx)
        self._chrom_pos.setdefault(variant.chrom, [])
        self._chrom_pos[variant.chrom].append(variant.pos)

    def add_variants(self, variants: list[Variant]) -> None:
        """Add multiple variants."""
        for v in variants:
            self.add_variant(v)

    def get_at_position(self, chrom: str, pos: int) -> list[Variant]:
        """Get all variants at a specific genomic position.

        Parameters
        ----------
        chrom : str
            Chromosome.
        pos : int
            Position (1-based).

        Returns
        -------
        list[Variant]
            Variants at that position.
        """
        key = (chrom, pos)
        indices = self._pos_index.get(key, [])
        return [self._variants[i] for i in indices]

    def get_in_range(self, chrom: str, start: int, end: int) -> list[Variant]:
        """Get all variants in a genomic range [start, end] (1-based, inclusive).

        Parameters
        ----------
        chrom : str
            Chromosome.
        start : int
            Start position (1-based, inclusive).
        end : int
            End position (1-based, inclusive).

        Returns
        -------
        list[Variant]
            Variants overlapping the range.
        """
        if chrom not in self._chrom_pos:
            return []

        positions = self._chrom_pos[chrom]
        # Binary search for start and end
        import bisect

        lo = bisect.bisect_left(positions, start)
        hi = bisect.bisect_right(positions, end)

        result: list[Variant] = []
        for i in range(lo, hi):
            pos = positions[i]
            for v in self.get_at_position(chrom, pos):
                result.append(v)
        return result

    def filter_by_maf(self, min_maf: float = 0.01) -> list[Variant]:
        """Return variants with MAF ≥ threshold.

        Parameters
        ----------
        min_maf : float
            Minimum MAF (default 0.01).

        Returns
        -------
        list[Variant]
            Variants passing the MAF filter.
        """
        return [v for v in self._variants if v.has_maf and v.maf >= min_maf]

    def filter_pathogenic(self) -> list[Variant]:
        """Return only pathogenic variants."""
        return [v for v in self._variants if v.is_pathogenic]

    def summary(self) -> dict:
        """Return a summary of the variant set."""
        return {
            "total_variants": self.num_variants,
            "snvs": sum(1 for v in self._variants if v.is_snv),
            "pathogenic": sum(1 for v in self._variants if v.is_pathogenic),
            "with_maf": sum(1 for v in self._variants if v.has_maf),
            "chromosomes": len(self._chrom_pos),
        }


# ---------------------------------------------------------------------------
# dbSNP API client
# ---------------------------------------------------------------------------

# NCBI E-utilities base URL
_EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


class DbSNPClient:
    """Client for fetching variant data from NCBI dbSNP via E-utilities.

    Uses NCBI ESearch and ESummary to find variants by gene name
    or genomic region. Inspired by GuideX's EntrezFetcher but
    redesigned for variant-specific queries.

    Parameters
    ----------
    email : str
        Email for NCBI API identification (required by NCBI).
    api_key : str, optional
        NCBI API key for higher rate limits.
    rate_limit : float
        Minimum seconds between requests (default 0.34 with API key,
        0.55 without).
    """

    def __init__(
        self,
        email: str = "",
        api_key: str | None = None,
        rate_limit: float | None = None,
    ):
        self.email = email
        self.api_key = api_key
        self.rate_limit = rate_limit or (0.1 if api_key else 0.34)
        self._last_request_time: float = 0.0

    def fetch_variants_by_gene(
        self,
        gene_name: str,
        organism: str = "Homo sapiens",
        max_results: int = 500,
    ) -> list[Variant]:
        """Fetch variants associated with a gene from dbSNP.

        Parameters
        ----------
        gene_name : str
            Gene symbol (e.g., "DMD").
        organism : str
            Organism name (default "Homo sapiens").
        max_results : int
            Maximum number of variants to fetch (default 500).

        Returns
        -------
        list[Variant]
            Variants from dbSNP.

        Raises
        ------
        VariantError
            If the API request fails.
        """
        import time
        import urllib.parse
        import urllib.request

        # Rate limiting
        elapsed = time.time() - self._last_request_time
        if elapsed < self.rate_limit:
            time.sleep(self.rate_limit - elapsed)

        # ESearch: find dbSNP records for the gene
        search_term = f"{gene_name}[Gene] AND {organism}[Organism] AND snp[SNP]"
        params = {
            "db": "snp",
            "term": search_term,
            "retmax": str(max_results),
            "retmode": "json",
        }
        if self.email:
            params["email"] = self.email
        if self.api_key:
            params["api_key"] = self.api_key

        url = f"{_EUTILS_URL}/esearch.fcgi?" + urllib.parse.urlencode(params)

        try:
            self._last_request_time = time.time()
            req = urllib.request.Request(url)
            with urllib.request.urlopen(req, timeout=30) as resp:
                import json

                data = json.loads(resp.read().decode())
        except Exception as e:
            raise VariantError(f"dbSNP ESearch failed for gene {gene_name}: {e}") from e

        ids = data.get("esearchresult", {}).get("idlist", [])
        if not ids:
            logger.info("No dbSNP variants found for gene %s", gene_name)
            return []

        # ESummary: fetch variant details
        variants = self._fetch_summaries(ids)
        logger.info("Fetched %d dbSNP variants for gene %s", len(variants), gene_name)
        return variants

    def _fetch_summaries(self, ids: list[str]) -> list[Variant]:
        """Fetch variant summaries from ESummary."""
        import json
        import time
        import urllib.parse
        import urllib.request

        # Rate limiting
        elapsed = time.time() - self._last_request_time
        if elapsed < self.rate_limit:
            time.sleep(self.rate_limit - elapsed)

        params = {
            "db": "snp",
            "id": ",".join(ids),
            "retmode": "json",
        }
        if self.email:
            params["email"] = self.email
        if self.api_key:
            params["api_key"] = self.api_key

        url = f"{_EUTILS_URL}/esummary.fcgi?" + urllib.parse.urlencode(params)

        try:
            self._last_request_time = time.time()
            req = urllib.request.Request(url)
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = json.loads(resp.read().decode())
        except Exception as e:
            raise VariantError(f"dbSNP ESummary failed: {e}") from e

        result_list = data.get("result", [])
        variants: list[Variant] = []

        for item in result_list:
            if not isinstance(item, dict):
                continue
            uid = item.get("uid", "")
            if not uid:
                continue

            # Parse variant fields from ESummary
            chrom = item.get("chrname", "")
            pos_str = item.get("chrpos", "")
            pos = int(pos_str) if pos_str else 0
            rs_id = f"rs{uid}"
            ref_alt = item.get("spdi", "")
            maf = _parse_maf_from_summary(item)
            variant_type = _classify_variant_type(item.get("snp_class", ""))
            clinical = item.get("clinical_significance", "")

            if chrom and pos:
                ref, alt = _parse_spdi(ref_alt)
                variants.append(
                    Variant(
                        chrom=chrom,
                        pos=pos,
                        ref=ref,
                        alt=alt,
                        rs_id=rs_id,
                        maf=maf,
                        variant_type=variant_type,
                        clinical_significance=clinical,
                        gene=item.get("gene_name", ""),
                    )
                )

        return variants


def _parse_maf_from_summary(item: dict) -> float:
    """Extract MAF from dbSNP ESummary item."""
    global_maf = item.get("global_maf", "")
    if global_maf:
        try:
            # MAF is often in format "A/G:0.05"
            parts = global_maf.split(":")
            if len(parts) >= 2:
                return float(parts[-1])
        except (ValueError, IndexError):
            pass
    return -1.0


def _classify_variant_type(snp_class: str) -> str:
    """Classify variant type from dbSNP snp_class string."""
    snp_class_lower = snp_class.lower()
    if "snv" in snp_class_lower or "snp" in snp_class_lower:
        return "snv"
    if (
        "indel" in snp_class_lower
        or "deletion" in snp_class_lower
        or "insertion" in snp_class_lower
    ):
        return "indel"
    if "mnv" in snp_class_lower or "mnp" in snp_class_lower:
        return "mnv"
    return "unknown"


def _parse_spdi(spdi: str) -> tuple[str, str]:
    """Parse SPDI notation into ref and alt alleles.

    SPDI format: NC_000001.1:100:A:G
    Returns (ref, alt) = ("A", "G").
    """
    if not spdi:
        return "", ""
    parts = spdi.split(":")
    if len(parts) >= 4:
        return parts[-2], parts[-1]
    return "", ""
