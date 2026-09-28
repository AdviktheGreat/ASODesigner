"""Tests for ASO candidate pre-filters."""

from aso_designer.candidate import ASOCandidate
from aso_designer.filters import (
    filter_ambiguous_bases,
    filter_end_stability,
    filter_gc_content,
    filter_homopolymer,
    run_all_filters,
)


def _make(seq: str, gc: float | None = None) -> ASOCandidate:
    """Helper to create a candidate with computed GC."""
    from aso_designer.utils.sequence import gc_content

    return ASOCandidate(
        seq, 0, len(seq), len(seq), gc_content=gc if gc is not None else gc_content(seq)
    )


# ---------------------------------------------------------------------------
# GC content filter
# ---------------------------------------------------------------------------


def test_gc_filter_passes_valid() -> None:
    """Candidate with GC in range passes."""
    cands = [_make("ATGCATGCATGCATGCATGC")]  # GC = 0.5
    passed = filter_gc_content(cands, gc_min=0.3, gc_max=0.6)
    assert passed == 1
    assert cands[0].passes_filters


def test_gc_filter_rejects_low() -> None:
    """Candidate with GC below minimum is rejected."""
    cands = [_make("ATATATATATATATATATAT")]  # GC = 0.0
    passed = filter_gc_content(cands, gc_min=0.3, gc_max=0.6)
    assert passed == 0
    assert not cands[0].passes_filters
    assert "GC content" in cands[0].features["filter_reason"]


def test_gc_filter_rejects_high() -> None:
    """Candidate with GC above maximum is rejected."""
    cands = [_make("GCGCGCGCGCGCGCGCGCGC")]  # GC = 1.0
    passed = filter_gc_content(cands, gc_min=0.3, gc_max=0.6)
    assert passed == 0
    assert not cands[0].passes_filters


# ---------------------------------------------------------------------------
# Homopolymer filter
# ---------------------------------------------------------------------------


def test_homopolymer_passes_clean() -> None:
    """Candidate without long homopolymers passes."""
    cands = [_make("ATGCATGCATGCATGCATGC")]
    passed = filter_homopolymer(cands, max_run=4)
    assert passed == 1


def test_homopolymer_rejects_run() -> None:
    """Candidate with long homopolymer run is rejected."""
    cands = [_make("AAAAATGCATGCATGCATGC")]  # 5x A
    passed = filter_homopolymer(cands, max_run=4)
    assert passed == 0
    assert not cands[0].passes_filters
    assert "Homopolymer" in cands[0].features["filter_reason"]


def test_homopolymer_boundary() -> None:
    """Homopolymer of exactly max_run passes."""
    cands = [_make("AAAATGCATGCATGCATGC")]  # 4x A
    passed = filter_homopolymer(cands, max_run=4)
    assert passed == 1


# ---------------------------------------------------------------------------
# End stability filter
# ---------------------------------------------------------------------------


def test_end_stability_passes_balanced() -> None:
    """Candidate with balanced ends passes."""
    cands = [_make("ATGCATGCATGCATGCATGC")]  # ATG at both ends
    passed = filter_end_stability(cands, end_length=3, max_end_gc=0.67)
    assert passed == 1


def test_end_stability_rejects_gc_clamp() -> None:
    """Candidate with all-GC 5' end is rejected."""
    cands = [_make("GCGCATATATATATATATAT")]  # GCG at 5'
    passed = filter_end_stability(cands, end_length=3, max_end_gc=0.67)
    assert passed == 0
    assert "5' end GC" in cands[0].features["filter_reason"]


def test_end_stability_rejects_low_gc() -> None:
    """Candidate with all-AT 5' end is rejected when min_end_gc > 0."""
    cands = [_make("ATATGCGCGCGCGCGCGCGC")]  # ATA at 5'
    passed = filter_end_stability(cands, end_length=3, min_end_gc=0.34)
    assert passed == 0
    assert "5' end GC" in cands[0].features["filter_reason"]


# ---------------------------------------------------------------------------
# Ambiguous base filter
# ---------------------------------------------------------------------------


def test_ambiguous_passes_clean() -> None:
    """Candidate without ambiguous bases passes."""
    cands = [_make("ATGCATGCATGCATGCATGC")]
    passed = filter_ambiguous_bases(cands, max_ambiguous=0)
    assert passed == 1


def test_ambiguous_rejects_n() -> None:
    """Candidate with N is rejected when max_ambiguous=0."""
    cands = [_make("ATGNATGCATGCATGCATGC")]
    passed = filter_ambiguous_bases(cands, max_ambiguous=0)
    assert passed == 0
    assert "Ambiguous" in cands[0].features["filter_reason"]


def test_ambiguous_allows_within_limit() -> None:
    """Candidate with 1 ambiguous base passes when max_ambiguous=1."""
    cands = [_make("ATGNATGCATGCATGCATGC")]
    passed = filter_ambiguous_bases(cands, max_ambiguous=1)
    assert passed == 1


# ---------------------------------------------------------------------------
# run_all_filters
# ---------------------------------------------------------------------------


def test_run_all_filters_mixed() -> None:
    """run_all_filters correctly filters a mixed set of candidates."""
    cands = [
        _make("ATGCATGCATGCATGCATGC"),  # good
        _make("ATATATATATATATATATAT"),  # low GC
        _make("GCGCGCGCGCGCGCGCGCGC"),  # high GC
        _make("AAAAATGCATGCATGCATGC"),  # homopolymer
        _make("ATGNATGCATGCATGCATGC"),  # ambiguous
    ]
    result = run_all_filters(cands, gc_min=0.3, gc_max=0.6, homopolymer_max=4, max_ambiguous=0)
    assert result["total_passed"] == 1
    assert result["total_failed"] == 4


def test_run_all_filters_all_pass() -> None:
    """run_all_filters passes all valid candidates."""
    cands = [
        _make("ATGCATGCATGCATGCATGC"),
        _make("GCTAGCTAGCTAGCTAGCTA"),
    ]
    result = run_all_filters(cands, gc_min=0.3, gc_max=0.6, homopolymer_max=4, max_ambiguous=0)
    assert result["total_passed"] == 2
