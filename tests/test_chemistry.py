"""Tests for ASO chemistry modeling, design classes, and annotation."""

import pytest

from aso_designer.candidate import ASOCandidate
from aso_designer.candidate_generator import generate_candidates
from aso_designer.chemistry import (
    FullChemistryDesign,
    GapmerDesign,
    MixmerDesign,
    apply_chemistry,
    chemistry_summary,
    create_chemistry_design,
    get_dna_gap_region,
    has_continuous_dna_gap,
    mixmer_from_scheme,
    parse_mixmer_scheme,
    validate_pattern,
)
from aso_designer.config_schema import CandidateConfig, GapmerConfig
from aso_designer.utils.exceptions import ConfigurationError

# ---------------------------------------------------------------------------
# GapmerDesign tests
# ---------------------------------------------------------------------------


def test_gapmer_build_pattern() -> None:
    """Gapmer produces [wing][DNA gap][wing] layout."""
    design = GapmerDesign(gap_size=10, wing_length=5, wing_chemistry="moe")
    pattern = design.build_pattern(20)
    assert len(pattern) == 20
    assert pattern[:5] == ["MOE"] * 5  # 5' wing
    assert pattern[5:15] == ["DNA"] * 10  # DNA gap
    assert pattern[15:] == ["MOE"] * 5  # 3' wing


def test_gapmer_2ome_wings() -> None:
    """Gapmer with 2'OMe wings uses 2OMe label."""
    design = GapmerDesign(gap_size=10, wing_length=5, wing_chemistry="2ome")
    pattern = design.build_pattern(20)
    assert pattern[0] == "2OMe"


def test_gapmer_auto_adjust_gap() -> None:
    """Gapmer auto-adjusts gap size when length doesn't match."""
    design = GapmerDesign(gap_size=10, wing_length=5, wing_chemistry="moe")
    # Length 18 instead of 20 — gap should shrink to 8
    pattern = design.build_pattern(18)
    assert len(pattern) == 18
    assert pattern[5:13] == ["DNA"] * 8


def test_gapmer_too_short() -> None:
    """Gapmer raises error when too short for wings."""
    design = GapmerDesign(gap_size=10, wing_length=5, wing_chemistry="moe")
    with pytest.raises(ConfigurationError):
        design.build_pattern(10)  # 10 - 2*5 = 0 < 4


def test_gapmer_rnase_h() -> None:
    """Gapmers are RNase H eligible."""
    design = GapmerDesign()
    assert design.is_rnase_h_eligible()


# ---------------------------------------------------------------------------
# MixmerDesign tests
# ---------------------------------------------------------------------------


def test_mixmer_default_pattern() -> None:
    """Default mixmer alternates modified and DNA."""
    design = MixmerDesign(modification="2ome")
    pattern = design.build_pattern(10)
    assert pattern[0] == "2OMe"
    assert pattern[1] == "DNA"
    assert pattern[2] == "2OMe"


def test_mixmer_custom_pattern() -> None:
    """Custom mixmer pattern is respected."""
    design = MixmerDesign(pattern=[True, False, True, False], modification="moe")
    pattern = design.build_pattern(4)
    assert pattern == ["MOE", "DNA", "MOE", "DNA"]


def test_mixmer_pattern_length_mismatch() -> None:
    """Mixmer raises error when pattern length != ASO length."""
    design = MixmerDesign(pattern=[True, False], modification="2ome")
    with pytest.raises(ConfigurationError):
        design.build_pattern(10)


def test_mixmer_not_rnase_h() -> None:
    """Default mixmer is not RNase H eligible."""
    design = MixmerDesign()
    assert not design.is_rnase_h_eligible()


# ---------------------------------------------------------------------------
# FullChemistryDesign tests
# ---------------------------------------------------------------------------


def test_full_chemistry_2ome() -> None:
    """Full 2'OMe chemistry labels all positions as 2OMe."""
    design = FullChemistryDesign(mod_type="2ome")
    pattern = design.build_pattern(15)
    assert all(p == "2OMe" for p in pattern)


def test_full_chemistry_moe() -> None:
    """Full MOE chemistry labels all positions as MOE."""
    design = FullChemistryDesign(mod_type="moe")
    pattern = design.build_pattern(15)
    assert all(p == "MOE" for p in pattern)


def test_full_chemistry_morpholino() -> None:
    """Full morpholino chemistry labels all positions as MORPHOLINO."""
    design = FullChemistryDesign(mod_type="morpholino")
    pattern = design.build_pattern(15)
    assert all(p == "MORPHOLINO" for p in pattern)


def test_full_chemistry_not_rnase_h() -> None:
    """Full-modification chemistries are not RNase H eligible."""
    for mod in ("2ome", "moe", "morpholino"):
        design = FullChemistryDesign(mod_type=mod)
        assert not design.is_rnase_h_eligible()


# ---------------------------------------------------------------------------
# Factory and application tests
# ---------------------------------------------------------------------------


def test_create_chemistry_gapmer() -> None:
    """Factory creates GapmerDesign from config."""
    cfg = CandidateConfig(chemistry="gapmer")
    gapmer_cfg = GapmerConfig(gap_size=10, wing_length=5)
    design = create_chemistry_design(cfg, gapmer_config=gapmer_cfg)
    assert isinstance(design, GapmerDesign)


def test_create_chemistry_gapmer_default_config() -> None:
    """Factory falls back to default GapmerConfig when none is provided."""
    design = create_chemistry_design(CandidateConfig(chemistry="gapmer"))
    assert isinstance(design, GapmerDesign)


def test_create_chemistry_mixmer() -> None:
    """Factory creates MixmerDesign from config."""
    cfg = CandidateConfig(chemistry="mixmer")
    design = create_chemistry_design(cfg)
    assert isinstance(design, MixmerDesign)


def test_create_chemistry_full() -> None:
    """Factory creates FullChemistryDesign for 2ome/moe/morpholino."""
    for chem in ("2ome", "moe", "morpholino"):
        cfg = CandidateConfig(chemistry=chem)
        design = create_chemistry_design(cfg)
        assert isinstance(design, FullChemistryDesign)


def test_apply_chemistry_sets_pattern() -> None:
    """apply_chemistry sets modification_pattern on all candidates."""
    seq = "ATGCATGCATGCATGCATGCATGCATGCATGC" * 2
    cfg = CandidateConfig(length=20, step=10, chemistry="gapmer")
    gapmer_cfg = GapmerConfig(gap_size=10, wing_length=5, wing_chemistry="moe")
    design = create_chemistry_design(cfg, gapmer_config=gapmer_cfg)
    cands = generate_candidates(seq, cfg)
    apply_chemistry(cands, design)
    for c in cands:
        assert len(c.modification_pattern) == 20
        assert c.features["rnase_h_eligible"] is True


# ---------------------------------------------------------------------------
# Mixmer scheme utilities
# ---------------------------------------------------------------------------


def test_parse_mixmer_scheme() -> None:
    """Scheme string parses to boolean pattern."""
    assert parse_mixmer_scheme("MMDM") == [True, True, False, True]
    assert parse_mixmer_scheme("d.m") == [False, False, True]


def test_parse_mixmer_scheme_invalid() -> None:
    """Invalid scheme characters raise ConfigurationError."""
    with pytest.raises(ConfigurationError):
        parse_mixmer_scheme("MXDM")


def test_has_continuous_dna_gap() -> None:
    """DNA gap detection works correctly."""
    assert has_continuous_dna_gap([True, False, False, False, False, True], min_gap=4)
    assert not has_continuous_dna_gap([True, False, True, False, True], min_gap=4)
    assert has_continuous_dna_gap([True, False, False, True], min_gap=2)


def test_mixmer_from_scheme() -> None:
    """mixmer_from_scheme creates MixmerDesign from scheme string."""
    design = mixmer_from_scheme("MMDMMDDMMDM", modification="moe")
    pattern = design.build_pattern(11)
    assert pattern[0] == "MOE"
    assert pattern[2] == "DNA"


# ---------------------------------------------------------------------------
# Annotation utilities
# ---------------------------------------------------------------------------


def test_validate_pattern() -> None:
    """validate_pattern checks pattern length matches sequence length."""
    c = ASOCandidate("ATGCATGCATGCATGCATGC", 0, 20, 20)
    c.modification_pattern = ["DNA"] * 20
    assert validate_pattern(c)
    c.modification_pattern = ["DNA"] * 10
    assert not validate_pattern(c)


def test_get_dna_gap_region() -> None:
    """get_dna_gap_region returns the contiguous DNA stretch."""
    c = ASOCandidate("ATGCATGCATGCATGCATGC", 0, 20, 20)
    c.modification_pattern = ["MOE"] * 5 + ["DNA"] * 10 + ["MOE"] * 5
    gap = get_dna_gap_region(c)
    assert gap == (5, 15)


def test_get_dna_gap_region_none() -> None:
    """get_dna_gap_region returns None for fully modified candidates."""
    c = ASOCandidate("ATGCATGCATGCATGCATGC", 0, 20, 20)
    c.modification_pattern = ["2OMe"] * 20
    assert get_dna_gap_region(c) is None


def test_chemistry_summary() -> None:
    """chemistry_summary returns a readable string."""
    c = ASOCandidate("ATGCATGCATGCATGCATGC", 0, 20, 20, chemistry="gapmer")
    c.modification_pattern = ["MOE"] * 5 + ["DNA"] * 10 + ["MOE"] * 5
    c.features["rnase_h_eligible"] = True
    summary = chemistry_summary(c)
    assert "gapmer" in summary
    assert "5-10-5" in summary
    assert "RNase H: yes" in summary


def test_validate_pattern_no_pattern() -> None:
    """An unmodified candidate fails pattern validation."""
    c = ASOCandidate("ATGCATGCATGCATGCATGC", 0, 20, 20)
    assert validate_pattern(c) is False


def test_get_dna_gap_region_no_pattern() -> None:
    """An unmodified candidate has no DNA gap region."""
    c = ASOCandidate("ATGCATGCATGCATGCATGC", 0, 20, 20)
    assert get_dna_gap_region(c) is None


def test_get_dna_gap_region_best_at_end() -> None:
    """The final DNA stretch counts even without a trailing modification."""
    c = ASOCandidate("ATGCATGCATGCATGCATGC", 0, 20, 20)
    c.modification_pattern = ["2OMe"] * 5 + ["DNA"] * 3
    assert get_dna_gap_region(c) == (5, 8)


def test_chemistry_summary_no_pattern() -> None:
    """An unmodified candidate reports that no pattern is set."""
    c = ASOCandidate("ATGCATGCATGCATGCATGC", 0, 20, 20, chemistry="2ome")
    assert "no pattern" in chemistry_summary(c)


def test_create_chemistry_unknown_raises() -> None:
    """Chemistries without factory support raise ConfigurationError.

    ``full`` is rejected by the schema and the CLI; model_construct
    bypasses schema validation to reach the factory's own guard.
    """
    cfg = CandidateConfig.model_construct(chemistry="full")
    with pytest.raises(ConfigurationError, match="Unknown chemistry"):
        create_chemistry_design(cfg)
