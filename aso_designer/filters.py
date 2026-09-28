"""Pre-filters for ASO candidate screening.

Implements configurable filters that reject candidates before the
expensive feature-computation stage (accessibility, thermodynamics,
off-target). Inspired by GuideX's _validate_basic method but
redesigned for ASO-specific needs: GC content with configurable
bounds, homopolymer detection with chemistry-specific thresholds,
5'/3' end stability, and ambiguous base handling.

All filters operate on ASOCandidate objects and set the
``passes_filters`` flag. Failed candidates are retained (not
discarded) so the user can inspect why they were rejected.
"""

from aso_designer.candidate import ASOCandidate
from aso_designer.utils.logger import get_logger
from aso_designer.utils.sequence import gc_content, has_homopolymer

logger = get_logger(__name__)


def filter_gc_content(
    candidates: list[ASOCandidate],
    gc_min: float = 0.30,
    gc_max: float = 0.60,
) -> int:
    """Filter candidates by GC content.

    Candidates outside the [gc_min, gc_max] range have their
    ``passes_filters`` flag set to False and a filter reason recorded
    in the features dict.

    Parameters
    ----------
    candidates : list[ASOCandidate]
        Candidates to filter (modified in place).
    gc_min : float
        Minimum acceptable GC content fraction.
    gc_max : float
        Maximum acceptable GC content fraction.

    Returns
    -------
    int
        Number of candidates that passed this filter.
    """
    passed = 0
    for c in candidates:
        if not c.passes_filters:
            continue
        if c.gc_content < gc_min:
            c.passes_filters = False
            c.features["filter_reason"] = f"GC content {c.gc_content:.2f} < {gc_min}"
        elif c.gc_content > gc_max:
            c.passes_filters = False
            c.features["filter_reason"] = f"GC content {c.gc_content:.2f} > {gc_max}"
        else:
            passed += 1
    return passed


def filter_homopolymer(
    candidates: list[ASOCandidate],
    max_run: int = 4,
) -> int:
    """Filter candidates by homopolymer run length.

    Candidates with consecutive identical bases longer than ``max_run``
    are rejected.

    Parameters
    ----------
    candidates : list[ASOCandidate]
        Candidates to filter (modified in place).
    max_run : int
        Maximum allowed consecutive identical bases.

    Returns
    -------
    int
        Number of candidates that passed this filter.
    """
    passed = 0
    for c in candidates:
        if not c.passes_filters:
            continue
        if has_homopolymer(c.sequence, max_run=max_run):
            c.passes_filters = False
            c.features["filter_reason"] = f"Homopolymer run > {max_run}"
        else:
            passed += 1
    return passed


def filter_end_stability(
    candidates: list[ASOCandidate],
    end_length: int = 3,
    max_end_gc: float = 1.0,
    min_end_gc: float = 0.0,
) -> int:
    """Filter candidates by 5' and 3' end GC content.

    Avoids extreme GC clamps at the termini, which can cause
    non-specific binding or synthesis issues.

    Parameters
    ----------
    candidates : list[ASOCandidate]
        Candidates to filter (modified in place).
    end_length : int
        Number of nucleotides at each end to check (default 3).
    max_end_gc : float
        Maximum GC fraction allowed at either end (1.0 = no limit).
    min_end_gc : float
        Minimum GC fraction required at either end (0.0 = no limit).

    Returns
    -------
    int
        Number of candidates that passed this filter.
    """
    passed = 0
    for c in candidates:
        if not c.passes_filters:
            continue
        five_end = c.sequence[:end_length]
        three_end = c.sequence[-end_length:]

        five_gc = gc_content(five_end)
        three_gc = gc_content(three_end)

        if five_gc > max_end_gc:
            c.passes_filters = False
            c.features["filter_reason"] = f"5' end GC {five_gc:.2f} > {max_end_gc}"
        elif three_gc > max_end_gc:
            c.passes_filters = False
            c.features["filter_reason"] = f"3' end GC {three_gc:.2f} > {max_end_gc}"
        elif five_gc < min_end_gc:
            c.passes_filters = False
            c.features["filter_reason"] = f"5' end GC {five_gc:.2f} < {min_end_gc}"
        elif three_gc < min_end_gc:
            c.passes_filters = False
            c.features["filter_reason"] = f"3' end GC {three_gc:.2f} < {min_end_gc}"
        else:
            passed += 1
    return passed


def filter_ambiguous_bases(
    candidates: list[ASOCandidate],
    max_ambiguous: int = 0,
) -> int:
    """Filter candidates with ambiguous (IUPAC) bases.

    Parameters
    ----------
    candidates : list[ASOCandidate]
        Candidates to filter (modified in place).
    max_ambiguous : int
        Maximum allowed ambiguous bases (default 0 — reject any).

    Returns
    -------
    int
        Number of candidates that passed this filter.
    """
    ambiguous_chars = set("NRYSWKMBDHV")
    passed = 0
    for c in candidates:
        if not c.passes_filters:
            continue
        count = sum(1 for ch in c.sequence if ch in ambiguous_chars)
        if count > max_ambiguous:
            c.passes_filters = False
            c.features["filter_reason"] = f"Ambiguous bases: {count} > {max_ambiguous}"
        else:
            passed += 1
    return passed


def run_all_filters(
    candidates: list[ASOCandidate],
    gc_min: float = 0.30,
    gc_max: float = 0.60,
    homopolymer_max: int = 4,
    end_length: int = 3,
    max_end_gc: float = 1.0,
    min_end_gc: float = 0.0,
    max_ambiguous: int = 0,
) -> dict[str, int]:
    """Apply all pre-filters in sequence and return pass counts.

    Filters are applied in order: GC content, homopolymer, end
    stability, ambiguous bases. Each filter only evaluates candidates
    that passed all previous filters.

    Parameters
    ----------
    candidates : list[ASOCandidate]
        Candidates to filter (modified in place).
    gc_min, gc_max : float
        GC content bounds.
    homopolymer_max : int
        Max consecutive identical bases.
    end_length : int
        Nucleotides at each end to check for GC clamp.
    max_end_gc, min_end_gc : float
        End GC content bounds.
    max_ambiguous : int
        Max allowed ambiguous bases.

    Returns
    -------
    dict[str, int]
        Pass counts per filter stage: {"gc": N, "homopolymer": N,
        "end_stability": N, "ambiguous": N, "total_passed": N}.
    """
    gc_passed = filter_gc_content(candidates, gc_min, gc_max)
    homo_passed = filter_homopolymer(candidates, homopolymer_max)
    end_passed = filter_end_stability(candidates, end_length, max_end_gc, min_end_gc)
    ambig_passed = filter_ambiguous_bases(candidates, max_ambiguous)

    total = sum(1 for c in candidates if c.passes_filters)
    total_failed = len(candidates) - total

    logger.info(
        "Pre-filters: %d/%d passed (GC: %d, homopolymer: %d, end: %d, ambiguous: %d), %d rejected",
        total,
        len(candidates),
        gc_passed,
        homo_passed,
        end_passed,
        ambig_passed,
        total_failed,
    )

    return {
        "gc": gc_passed,
        "homopolymer": homo_passed,
        "end_stability": end_passed,
        "ambiguous": ambig_passed,
        "total_passed": total,
        "total_failed": total_failed,
    }
