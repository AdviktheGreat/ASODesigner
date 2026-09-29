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
