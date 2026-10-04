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
        """List of all variants in the set."""
        return self._variants

    @property
    def num_variants(self) -> int:
        """Number of variants in the set."""
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


# ---------------------------------------------------------------------------
# ClinVar API client
# ---------------------------------------------------------------------------


class ClinVarClient:
    """Client for fetching pathogenic variants from ClinVar via E-utilities.

    Queries NCBI ClinVar for variants with clinical significance
    (pathogenic, likely pathogenic) associated with a gene.

    Parameters
    ----------
    email : str
        Email for NCBI API identification.
    api_key : str, optional
        NCBI API key for higher rate limits.
    rate_limit : float
        Minimum seconds between requests.
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

    def fetch_pathogenic_variants(
        self,
        gene_name: str,
        organism: str = "Homo sapiens",
        max_results: int = 200,
    ) -> list[Variant]:
        """Fetch pathogenic variants for a gene from ClinVar.

        Parameters
        ----------
        gene_name : str
            Gene symbol.
        organism : str
            Organism name (default "Homo sapiens").
        max_results : int
            Maximum results (default 200).

        Returns
        -------
        list[Variant]
            Pathogenic variants from ClinVar.

        Raises
        ------
        VariantError
            If the API request fails.
        """
        import json
        import time
        import urllib.parse
        import urllib.request

        # Rate limiting
        elapsed = time.time() - self._last_request_time
        if elapsed < self.rate_limit:
            time.sleep(self.rate_limit - elapsed)

        # ESearch: find ClinVar records for pathogenic variants in gene
        search_term = (
            f"{gene_name}[Gene] AND {organism}[Organism] "
            f'AND ("pathogenic"[Clinical_Significance] '
            f'OR "likely pathogenic"[Clinical_Significance])'
        )
        params = {
            "db": "clinvar",
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
                data = json.loads(resp.read().decode())
        except Exception as e:
            raise VariantError(f"ClinVar ESearch failed for gene {gene_name}: {e}") from e

        ids = data.get("esearchresult", {}).get("idlist", [])
        if not ids:
            logger.info("No ClinVar pathogenic variants found for %s", gene_name)
            return []

        # ESummary: fetch variant details
        variants = self._fetch_clinvar_summaries(ids)
        logger.info(
            "Fetched %d ClinVar pathogenic variants for gene %s",
            len(variants),
            gene_name,
        )
        return variants

    def _fetch_clinvar_summaries(self, ids: list[str]) -> list[Variant]:
        """Fetch ClinVar variant summaries from ESummary."""
        import json
        import time
        import urllib.parse
        import urllib.request

        elapsed = time.time() - self._last_request_time
        if elapsed < self.rate_limit:
            time.sleep(self.rate_limit - elapsed)

        params = {
            "db": "clinvar",
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
            raise VariantError(f"ClinVar ESummary failed: {e}") from e

        result_list = data.get("result", [])
        variants: list[Variant] = []

        for item in result_list:
            if not isinstance(item, dict):
                continue
            uid = item.get("uid", "")
            if not uid:
                continue

            # Parse ClinVar fields
            chrom = item.get("chromosome", "")
            pos_str = item.get("location", "")
            pos = _parse_clinvar_position(pos_str)
            clinical_sig = item.get("clinical_significance", "")
            # Normalize clinical significance
            clinical_sig = _normalize_clinical_significance(clinical_sig)
            variant_type = _classify_variant_type(item.get("variation_type", ""))
            ref, alt = _parse_clinvar_alleles(item.get("alleles", ""))

            if chrom and pos:
                variants.append(
                    Variant(
                        chrom=chrom,
                        pos=pos,
                        ref=ref,
                        alt=alt,
                        rs_id=item.get("rs_id", ""),
                        maf=-1.0,  # ClinVar doesn't provide MAF
                        variant_type=variant_type,
                        clinical_significance=clinical_sig,
                        gene=item.get("gene_sort", ""),
                    )
                )

        return variants


def _parse_clinvar_position(pos_str: str) -> int:
    """Parse a ClinVar position string into an integer."""
    if not pos_str:
        return 0
    # Positions can be "12345" or "12345-12350" (range)
    pos_part = pos_str.split("-")[0].split(",")[0]
    try:
        return int(pos_part)
    except ValueError:
        return 0


def _normalize_clinical_significance(sig: str) -> str:
    """Normalize ClinVar significance strings."""
    sig_lower = sig.lower().strip()
    if "pathogenic" in sig_lower and "likely" in sig_lower:
        return "likely_pathogenic"
    if "pathogenic" in sig_lower:
        return "pathogenic"
    if "benign" in sig_lower and "likely" in sig_lower:
        return "likely_benign"
    if "benign" in sig_lower:
        return "benign"
    if "uncertain" in sig_lower or "vus" in sig_lower:
        return "uncertain"
    return sig_lower


def _parse_clinvar_alleles(alleles: str) -> tuple[str, str]:
    """Parse ClinVar alleles string into ref and alt."""
    if not alleles:
        return "", ""
    # Format can be "A/G" or "AT/-" etc.
    if "/" in alleles:
        parts = alleles.split("/")
        ref = parts[0] if parts else ""
        alt = parts[1] if len(parts) > 1 else ""
        return ref, alt
    return "", ""


# ---------------------------------------------------------------------------
# Variant fetch manager (API → VCF fallback)
# ---------------------------------------------------------------------------


class VariantFetchManager:
    """Orchestrates variant fetching from multiple sources.

    Tries dbSNP and ClinVar APIs first, then falls back to a
    local VCF file if provided. Merges results into a single
    VariantSet, deduplicating by (chrom, pos, ref, alt).

    Parameters
    ----------
    email : str
        Email for NCBI API identification.
    api_key : str, optional
        NCBI API key.
    vcf_path : str or Path, optional
        Path to a local VCF file for fallback.
    """

    def __init__(
        self,
        email: str = "",
        api_key: str | None = None,
        vcf_path: str | None = None,
    ):
        self.dbsnp = DbSNPClient(email=email, api_key=api_key)
        self.clinvar = ClinVarClient(email=email, api_key=api_key)
        self.vcf_path = vcf_path

    def fetch_variants(
        self,
        gene_name: str,
        organism: str = "Homo sapiens",
        max_results: int = 500,
        use_clinvar: bool = True,
    ) -> VariantSet:
        """Fetch variants for a gene from all available sources.

        Parameters
        ----------
        gene_name : str
            Gene symbol.
        organism : str
            Organism name.
        max_results : int
            Max variants per source.
        use_clinvar : bool
            Whether to query ClinVar (default True).

        Returns
        -------
        VariantSet
            Merged variant set from all sources.
        """
        variant_set = VariantSet()

        # 1. dbSNP API
        try:
            dbsnp_variants = self.dbsnp.fetch_variants_by_gene(
                gene_name, organism=organism, max_results=max_results
            )
            variant_set.add_variants(dbsnp_variants)
            logger.info("dbSNP: %d variants for %s", len(dbsnp_variants), gene_name)
        except VariantError as e:
            logger.warning("dbSNP fetch failed: %s", e)

        # 2. ClinVar API
        if use_clinvar:
            try:
                clinvar_variants = self.clinvar.fetch_pathogenic_variants(
                    gene_name, organism=organism, max_results=max_results
                )
                variant_set.add_variants(clinvar_variants)
                logger.info(
                    "ClinVar: %d variants for %s",
                    len(clinvar_variants),
                    gene_name,
                )
            except VariantError as e:
                logger.warning("ClinVar fetch failed: %s", e)

        # 3. VCF fallback (if API returned nothing or VCF is provided)
        if self.vcf_path and variant_set.num_variants == 0:
            try:
                vcf_variants = _load_vcf_fallback(self.vcf_path)
                variant_set.add_variants(vcf_variants)
                logger.info(
                    "VCF fallback: %d variants from %s",
                    len(vcf_variants),
                    self.vcf_path,
                )
            except VariantError as e:
                logger.warning("VCF fallback failed: %s", e)

        # Deduplicate
        variant_set = _deduplicate_variant_set(variant_set)

        logger.info(
            "Total unique variants for %s: %d",
            gene_name,
            variant_set.num_variants,
        )
        return variant_set


def _load_vcf_fallback(vcf_path: str) -> list[Variant]:
    """Load variants from a VCF file (fallback).

    Uses the VCF parser to load and parse variants from a local
    VCF file when API sources are unavailable or return no results.
    """
    return parse_vcf(vcf_path)


def _deduplicate_variant_set(
    variant_set: VariantSet,
) -> VariantSet:
    """Deduplicate variants by (chrom, pos, ref, alt).

    Keeps the first occurrence of each unique variant, preferring
    variants with more complete information (MAF, clinical significance).
    """
    seen: dict[tuple[str, int, str, str], Variant] = {}

    for v in variant_set.variants:
        key = (v.chrom, v.pos, v.ref, v.alt)
        if key not in seen:
            seen[key] = v
        else:
            # Merge: keep the one with more complete data
            existing = seen[key]
            if v.has_maf and not existing.has_maf:
                existing = v
            if v.clinical_significance and not existing.clinical_significance:
                existing.clinical_significance = v.clinical_significance
            if v.rs_id and not existing.rs_id:
                existing.rs_id = v.rs_id
            seen[key] = existing

    result = VariantSet()
    for v in seen.values():
        result.add_variant(v)
    return result


# ---------------------------------------------------------------------------
# VCF parser
# ---------------------------------------------------------------------------


def parse_vcf(file_path: str) -> list[Variant]:
    """Parse a VCF file and return a list of Variant objects.

    Handles standard VCF 4.x format:
    - Header lines (##) with metadata
    - Column header line (#CHROM POS ID REF ALT QUAL FILTER INFO)
    - Data lines with tab-separated fields
    - INFO field parsing for AF (allele frequency), CLNSIG (clinical
      significance), and other relevant annotations

    Parameters
    ----------
    file_path : str
        Path to the VCF file.

    Returns
    -------
    list[Variant]
        Parsed variants.

    Raises
    ------
    VariantError
        If the file cannot be read or is not valid VCF.
    """
    from pathlib import Path

    path = Path(file_path)
    if not path.exists():
        raise VariantError(f"VCF file not found: {path}")

    variants: list[Variant] = []

    try:
        with open(path, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("##"):
                    continue  # Metadata line
                if line.startswith("#CHROM"):
                    continue  # Column header

                # Parse data line
                fields = line.split("\t")
                if len(fields) < 8:
                    logger.warning("Skipping malformed VCF line: %s", line[:80])
                    continue

                chrom = fields[0]
                # Normalize chromosome (strip "chr" prefix for consistency)
                chrom = chrom.replace("chr", "")
                pos = int(fields[1])
                rs_id = fields[2] if fields[2] != "." else ""
                ref = fields[3]
                alt = fields[4]
                info_str = fields[7]

                # Parse INFO field
                info = _parse_info_field(info_str)
                maf = _extract_maf_from_info(info)
                clinical_sig = _extract_clinical_from_info(info)
                variant_type = _classify_vcf_variant_type(ref, alt)

                variants.append(
                    Variant(
                        chrom=chrom,
                        pos=pos,
                        ref=ref,
                        alt=alt,
                        rs_id=rs_id,
                        maf=maf,
                        variant_type=variant_type,
                        clinical_significance=clinical_sig,
                    )
                )

    except OSError as e:
        raise VariantError(f"Failed to read VCF file {path}: {e}") from e

    logger.info("Parsed %d variants from VCF: %s", len(variants), path)
    return variants


def _parse_info_field(info_str: str) -> dict[str, str]:
    """Parse a VCF INFO field into a dictionary.

    INFO format: key1=value1;key2=value2;key3
    """
    result: dict[str, str] = {}
    if not info_str or info_str == ".":
        return result
    for item in info_str.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
        else:
            result[item] = "true"
    return result


def _extract_maf_from_info(info: dict[str, str]) -> float:
    """Extract MAF from VCF INFO field.

    Looks for AF (allele frequency) or MAF keys.
    """
    for key in ("AF", "MAF", "GMAF"):
        if key in info:
            try:
                return float(info[key])
            except ValueError:
                continue
    return -1.0


def _extract_clinical_from_info(info: dict[str, str]) -> str:
    """Extract clinical significance from VCF INFO field.

    Looks for CLNSIG (ClinVar) or CLNSIGCONF.
    """
    for key in ("CLNSIG", "CLNSIGCONF", "ClinicalSignificance"):
        if key in info:
            return _normalize_clinical_significance(info[key])
    return ""


def _classify_vcf_variant_type(ref: str, alt: str) -> str:
    """Classify variant type from REF and ALT alleles.

    - SNV: single nucleotide change (both ref and alt are single bases)
    - Indel: insertion or deletion
    - MNV: multi-nucleotide change
    """
    # Multi-allelic (comma-separated alt)
    if "," in alt:
        return "mnv"
    if len(ref) == 1 and len(alt) == 1:
        return "snv"
    if len(ref) != len(alt):
        return "indel"
    if len(ref) > 1 and len(alt) > 1:
        return "mnv"
    return "unknown"


# ---------------------------------------------------------------------------
# VCF loading and filtering convenience functions
# ---------------------------------------------------------------------------


def load_vcf(file_path: str) -> VariantSet:
    """Parse a VCF file and return a VariantSet.

    Convenience function that parses a VCF file and builds a
    VariantSet with position-based indexing.

    Parameters
    ----------
    file_path : str
        Path to the VCF file.

    Returns
    -------
    VariantSet
        Variant set with all parsed variants.
    """
    variants = parse_vcf(file_path)
    vs = VariantSet()
    vs.add_variants(variants)
    return vs


def load_vcf_filtered(
    file_path: str,
    min_maf: float = 0.0,
    variant_types: set[str] | None = None,
    pathogenic_only: bool = False,
) -> VariantSet:
    """Parse a VCF file and return a filtered VariantSet.

    Parameters
    ----------
    file_path : str
        Path to the VCF file.
    min_maf : float
        Minimum MAF (default 0.0 = no MAF filter).
    variant_types : set[str] or None
        Allowed variant types (e.g., {"snv"}). None = all types.
    pathogenic_only : bool
        If True, keep only pathogenic variants.

    Returns
    -------
    VariantSet
        Filtered variant set.
    """
    vs = load_vcf(file_path)

    filtered: list[Variant] = []
    for v in vs.variants:
        # MAF filter
        if min_maf > 0 and (not v.has_maf or v.maf < min_maf):
            continue
        # Type filter
        if variant_types and v.variant_type not in variant_types:
            continue
        # Pathogenic filter
        if pathogenic_only and not v.is_pathogenic:
            continue
        filtered.append(v)

    result = VariantSet()
    result.add_variants(filtered)
    logger.info(
        "VCF filtered: %d → %d variants (min_maf=%s, types=%s, pathogenic=%s)",
        vs.num_variants,
        result.num_variants,
        min_maf,
        variant_types,
        pathogenic_only,
    )
    return result


# ---------------------------------------------------------------------------
# Variant-aware annotation
# ---------------------------------------------------------------------------


@dataclass
class VariantAnnotation:
    """Result of variant annotation for an ASO candidate.

    Attributes
    ----------
    aso_start : int
        ASO start position (0-based, transcript coordinates).
    aso_end : int
        ASO end position (0-based, exclusive, transcript coordinates).
    total_variants : int
        Total variants overlapping the ASO binding site.
    seed_variants : int
        Variants in the seed region.
    non_seed_variants : int
        Variants outside the seed region.
    pathogenic_variants : int
        Pathogenic variants in the binding site.
    has_pathogenic : bool
        Whether any pathogenic variant overlaps.
    penalty : float
        Computed variant penalty (0.0 = no penalty).
    variants : list[Variant]
        All overlapping variant objects.
    """

    aso_start: int
    aso_end: int
    total_variants: int
    seed_variants: int
    non_seed_variants: int
    pathogenic_variants: int
    has_pathogenic: bool
    penalty: float
    variants: list[Variant]


class VariantAnnotator:
    """Annotates ASO candidates with overlapping variant information.

    Checks for variants in the ASO binding site, classifies them by
    seed/non-seed region, and computes a variant-aware penalty.

    Parameters
    ----------
    variant_set : VariantSet
        The set of variants to check against.
    seed_start : int
        Seed region start position in the ASO (0-based, default 2).
    seed_end : int
        Seed region end position in the ASO (exclusive, default 8).
    maf_threshold : float
        Minimum MAF to consider a variant (default 0.01).
    seed_region_multiplier : float
        Penalty multiplier for seed-region variants (default 3.0).
    flag_pathogenic : bool
        Whether to flag pathogenic variants (default True).
    """

    def __init__(
        self,
        variant_set: VariantSet,
        seed_start: int = 2,
        seed_end: int = 8,
        maf_threshold: float = 0.01,
        seed_region_multiplier: float = 3.0,
        flag_pathogenic: bool = True,
    ):
        self.variant_set = variant_set
        self.seed_start = seed_start
        self.seed_end = seed_end
        self.maf_threshold = maf_threshold
        self.seed_region_multiplier = seed_region_multiplier
        self.flag_pathogenic = flag_pathogenic

    def annotate(
        self,
        aso_start: int,
        aso_end: int,
        chrom: str = "",
    ) -> VariantAnnotation:
        """Annotate an ASO binding site with overlapping variants.

        Parameters
        ----------
        aso_start : int
            ASO start position (0-based, transcript or genomic coordinates).
        aso_end : int
            ASO end position (0-based, exclusive).
        chrom : str
            Chromosome (if using genomic coordinates). Empty for
            transcript coordinates (searches all chromosomes).

        Returns
        -------
        VariantAnnotation
            Annotation result with variant counts and penalty.
        """
        # Convert to 1-based inclusive range for variant lookup
        # ASO positions are 0-based [start, end), variants are 1-based
        genomic_start = aso_start + 1
        genomic_end = aso_end

        # Get overlapping variants
        if chrom:
            overlapping = self.variant_set.get_in_range(chrom, genomic_start, genomic_end)
        else:
            # Search all chromosomes
            overlapping = [
                v for v in self.variant_set.variants if genomic_start <= v.pos <= genomic_end
            ]

        # Filter by MAF
        significant = [v for v in overlapping if not v.has_maf or v.maf >= self.maf_threshold]

        # Classify by seed/non-seed region
        # Seed in ASO: positions [seed_start, seed_end) (0-based)
        # In genomic coordinates: [aso_start + seed_start + 1, aso_start + seed_end]
        seed_genomic_start = aso_start + self.seed_start + 1
        seed_genomic_end = aso_start + self.seed_end

        seed_count = 0
        non_seed_count = 0
        pathogenic_count = 0

        for v in significant:
            if seed_genomic_start <= v.pos <= seed_genomic_end:
                seed_count += 1
            else:
                non_seed_count += 1
            if v.is_pathogenic:
                pathogenic_count += 1

        # Compute penalty
        penalty = self._compute_penalty(significant, seed_count, non_seed_count)

        return VariantAnnotation(
            aso_start=aso_start,
            aso_end=aso_end,
            total_variants=len(significant),
            seed_variants=seed_count,
            non_seed_variants=non_seed_count,
            pathogenic_variants=pathogenic_count,
            has_pathogenic=pathogenic_count > 0,
            penalty=penalty,
            variants=significant,
        )

    def _compute_penalty(
        self,
        variants: list[Variant],
        seed_count: int,
        non_seed_count: int,
    ) -> float:
        """Compute a variant penalty score.

        Penalty = (seed_variants * seed_multiplier + non_seed_variants) * maf_factor
        where maf_factor is the average MAF of overlapping variants.

        Returns
        -------
        float
            Penalty in [0.0, 2.0]. Higher = worse.
        """
        if not variants:
            return 0.0

        # Base penalty: count-based
        base = seed_count * self.seed_region_multiplier + non_seed_count

        # MAF factor: average MAF of variants with known MAF
        maf_values = [v.maf for v in variants if v.has_maf]
        if maf_values:
            avg_maf = sum(maf_values) / len(maf_values)
            maf_factor = 1.0 + avg_maf  # Range [1.0, 2.0]
        else:
            maf_factor = 1.0

        # Pathogenic bonus
        pathogenic_bonus = 0.0
        if self.flag_pathogenic:
            pathogenic_count = sum(1 for v in variants if v.is_pathogenic)
            pathogenic_bonus = pathogenic_count * 0.5

        penalty = (base * maf_factor + pathogenic_bonus) / 10.0
        return min(penalty, 2.0)  # Cap at 2.0

    def annotate_batch(
        self,
        candidates: list[tuple[int, int]],
        chrom: str = "",
    ) -> list[VariantAnnotation]:
        """Annotate multiple ASO candidates.

        Parameters
        ----------
        candidates : list[tuple[int, int]]
            List of (start, end) positions.
        chrom : str
            Chromosome (if using genomic coordinates).

        Returns
        -------
        list[VariantAnnotation]
            Annotations for each candidate.
        """
        return [self.annotate(start, end, chrom=chrom) for start, end in candidates]


# ---------------------------------------------------------------------------
# Variant scoring and risk classification
# ---------------------------------------------------------------------------


def score_variant_penalty(annotation: VariantAnnotation) -> float:
    """Convert a VariantAnnotation penalty to a [0.0, 1.0] score.

    Score = 1.0 - penalty (clamped). A score of 1.0 means no variants
    (best), 0.0 means many variants (worst).

    Parameters
    ----------
    annotation : VariantAnnotation
        The annotation result.

    Returns
    -------
    float
        Variant score in [0.0, 1.0].
    """
    return max(0.0, 1.0 - annotation.penalty)


def classify_variant_risk(
    annotation: VariantAnnotation,
    high_risk_threshold: float = 0.5,
    moderate_risk_threshold: float = 0.2,
) -> str:
    """Classify variant risk for an ASO candidate.

    Parameters
    ----------
    annotation : VariantAnnotation
        The annotation result.
    high_risk_threshold : float
        Penalty at or above this is "high_risk" (default 0.5).
    moderate_risk_threshold : float
        Penalty at or above this (but below high) is "moderate_risk"
        (default 0.2). Below this is "low_risk".

    Returns
    -------
    str
        One of "low_risk", "moderate_risk", or "high_risk".
    """
    if annotation.has_pathogenic:
        return "high_risk"
    if annotation.penalty >= high_risk_threshold:
        return "high_risk"
    elif annotation.penalty >= moderate_risk_threshold:
        return "moderate_risk"
    else:
        return "low_risk"


def filter_candidates_by_variants(
    candidates: list[tuple[int, int]],
    annotator: VariantAnnotator,
    max_penalty: float = 1.0,
    exclude_pathogenic: bool = False,
    chrom: str = "",
) -> list[tuple[int, int, VariantAnnotation]]:
    """Filter ASO candidates by variant penalty.

    Parameters
    ----------
    candidates : list[tuple[int, int]]
        List of (start, end) positions.
    annotator : VariantAnnotator
        The annotator to use.
    max_penalty : float
        Maximum allowed penalty (default 1.0).
    exclude_pathogenic : bool
        If True, exclude candidates overlapping pathogenic variants.
    chrom : str
        Chromosome (if using genomic coordinates).

    Returns
    -------
    list[tuple[int, int, VariantAnnotation]]
        Filtered candidates with their annotations.
    """
    results: list[tuple[int, int, VariantAnnotation]] = []
    for start, end in candidates:
        ann = annotator.annotate(start, end, chrom=chrom)
        if ann.penalty > max_penalty:
            continue
        if exclude_pathogenic and ann.has_pathogenic:
            continue
        results.append((start, end, ann))
    return results
