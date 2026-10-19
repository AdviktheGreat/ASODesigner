"""Tests for ASO chemistry modeling, design classes, and annotation."""

import pytest

from aso_designer.candidate import ASOCandidate
from aso_designer.candidate_generator import generate_candidates
from aso_designer.chemistry import (
    CleavageProfile,
    FullChemistryDesign,
    GapmerDesign,
    LnaDesign,
    MixmerDesign,
    apply_chemistry,
    chemistry_summary,
    cleavage_likelihood,
    compute_cleavage_profile,
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


def test_gapmer_rejects_mismatched_length() -> None:
    """Gapmer build_pattern rejects lengths that don't match the layout."""
    design = GapmerDesign(gap_size=10, wing_length=5, wing_chemistry="moe")
    with pytest.raises(ConfigurationError, match="Invalid gapmer geometry"):
        design.build_pattern(18)  # expected 20, no silent gap mutation


def test_gapmer_design_is_immutable() -> None:
    """GapmerDesign instances cannot be mutated after creation."""
    design = GapmerDesign(gap_size=10, wing_length=5, wing_chemistry="moe")
    from dataclasses import FrozenInstanceError

    with pytest.raises(FrozenInstanceError):
        design.gap_size = 8  # type: ignore[misc]


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


def test_create_chemistry_gapmer_geometry_mismatch() -> None:
    """Factory rejects gapmer configs whose length breaks the geometry."""
    cfg = CandidateConfig(chemistry="gapmer", length=22)
    with pytest.raises(ConfigurationError, match="Invalid gapmer geometry"):
        create_chemistry_design(cfg, gapmer_config=GapmerConfig(gap_size=10, wing_length=5))


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


# ---------------------------------------------------------------------------
# RNase H positional cleavage-likelihood model
# ---------------------------------------------------------------------------


class TestCleavageLikelihood:
    def test_short_gaps_are_not_substrates(self):
        assert cleavage_likelihood(0, 1) == 0.0
        assert cleavage_likelihood(0, 2) == 0.0
        # 3-nt gaps are marginal substrates at best
        assert 0.0 < cleavage_likelihood(1, 3) < 0.15

    def test_optimal_gap_size(self):
        # 6-9 nt gaps reach the maximum size factor
        for gap in (6, 7, 8, 9):
            peak = round(0.65 * (gap - 1))  # model peak position
            mid = cleavage_likelihood(peak, gap)
            assert 0.85 <= mid <= 1.0

    def test_position_profile_peaks_center_to_3p(self):
        gap = 10
        scores = [cleavage_likelihood(i, gap) for i in range(gap)]
        assert max(scores) == scores.index(max(scores)) * 0 + max(scores)
        peak_index = scores.index(max(scores))
        assert 0.4 <= peak_index / (gap - 1) <= 0.8
        # terminal gap positions are disfavored
        assert scores[0] < scores[peak_index]
        assert scores[-1] < scores[peak_index]

    def test_5p_terminal_gap_penalized(self):
        # A gap running into the 5' end: first positions nearly zero
        with_wing = [cleavage_likelihood(i, 8) for i in range(3)]
        at_end = [cleavage_likelihood(i, 8, gap_at_5p_end=True) for i in range(3)]
        assert at_end[0] < with_wing[0]
        assert at_end[1] < with_wing[1]
        assert at_end[2] == pytest.approx(with_wing[2])

    def test_3p_terminal_gap_penalized(self):
        with_wing = [cleavage_likelihood(i, 8) for i in range(8)]
        at_end = [cleavage_likelihood(i, 8, gap_at_3p_end=True) for i in range(8)]
        assert at_end[-1] < with_wing[-1]
        assert at_end[-2] < with_wing[-2]
        assert at_end[0] == pytest.approx(with_wing[0])

    def test_out_of_range(self):
        assert cleavage_likelihood(-1, 8) == 0.0
        assert cleavage_likelihood(8, 8) == 0.0
        assert cleavage_likelihood(0, 0) == 0.0


class TestComputeCleavageProfile:
    def test_gapmer_pattern(self):
        pattern = ["MOE"] * 4 + ["DNA"] * 8 + ["MOE"] * 4
        profile = compute_cleavage_profile(pattern)
        assert isinstance(profile, CleavageProfile)
        assert profile.gap_length == 8
        assert len(profile.likelihoods) == 16
        assert all(x == 0.0 for x in profile.likelihoods[:4])
        assert all(x == 0.0 for x in profile.likelihoods[12:])
        assert profile.max_score > 0.8
        assert profile.expected_site is not None
        assert 4 <= profile.expected_site < 12
        assert 0.5 <= profile.mean_gap_score <= 1.0

    def test_full_modification_no_gap(self):
        profile = compute_cleavage_profile(["2OMe"] * 20)
        assert profile.gap_length == 0
        assert profile.max_score == 0.0
        assert profile.expected_site is None
        assert profile.mean_gap_score == 0.0
        assert all(x == 0.0 for x in profile.likelihoods)

    def test_largest_gap_wins(self):
        pattern = ["DNA", "DNA", "2OMe"] + ["DNA"] * 6 + ["2OMe"]
        profile = compute_cleavage_profile(pattern)
        assert profile.gap_length == 6
        assert profile.expected_site is not None
        assert 3 <= profile.expected_site < 9

    def test_gap_at_terminus_scores_lower(self):
        interior = compute_cleavage_profile(["MOE", "MOE"] + ["DNA"] * 6 + ["MOE", "MOE"])
        at_5p = compute_cleavage_profile(["DNA"] * 6 + ["MOE", "MOE"])
        # the 5'-terminal gap positions lose their junction context
        assert at_5p.likelihoods[0] < interior.likelihoods[2]
        assert at_5p.mean_gap_score < interior.mean_gap_score


# ---------------------------------------------------------------------------
# LNA chemistry
# ---------------------------------------------------------------------------


class TestLnaDesign:
    def test_default_pattern_every_third(self):
        design = LnaDesign()
        pattern = design.build_pattern(12)
        assert pattern[1] == "LNA"
        assert pattern[4] == "LNA"
        assert pattern[0] == "DNA"
        assert pattern.count("LNA") == 4

    def test_custom_pattern_length_mismatch(self):
        design = LnaDesign(pattern=[True, False])
        with pytest.raises(ConfigurationError, match="LNA pattern length"):
            design.build_pattern(10)

    def test_custom_pattern_labels(self):
        design = LnaDesign(pattern=[True, True, False, False])
        assert design.build_pattern(4) == ["LNA", "LNA", "DNA", "DNA"]

    def test_rnase_h_eligibility(self):
        # Default sprinkled pattern has no >=4 DNA run
        assert not LnaDesign().is_rnase_h_eligible()
        # Custom pattern with a contiguous 6-nt DNA gap is eligible
        design = LnaDesign(pattern=[True, True] + [False] * 6 + [True, True])
        assert design.is_rnase_h_eligible()

    def test_factory_and_valid_chemistries(self):
        from aso_designer.chemistry import VALID_CHEMISTRIES

        assert "lna" in VALID_CHEMISTRIES
        design = create_chemistry_design(CandidateConfig(chemistry="lna"))
        assert isinstance(design, LnaDesign)


class TestLnaThermodynamics:
    SEQ = "ACGTACGTAC"
    PATTERN = ["DNA", "LNA", "DNA", "LNA", "DNA", "LNA", "DNA", "LNA", "DNA", "LNA"]

    def test_tm_per_lna_site(self):
        from aso_designer.thermodynamics import calculate_tm_with_chemistry

        base = calculate_tm_with_chemistry(self.SEQ, chemistry="lna")
        adjusted = calculate_tm_with_chemistry(
            self.SEQ, chemistry="lna", modification_pattern=self.PATTERN
        )
        assert adjusted == pytest.approx(base + 3.5 * 5)

    def test_dg_per_lna_site(self):
        from aso_designer.thermodynamics import calculate_duplex_dg

        base = calculate_duplex_dg(self.SEQ, chemistry="lna").dg_duplex
        adjusted = calculate_duplex_dg(
            self.SEQ, chemistry="lna", modification_pattern=self.PATTERN
        ).dg_duplex
        assert adjusted == pytest.approx(base - 1.4 * 5)

    def test_lna_uses_hybrid_nn_params(self):
        from aso_designer.thermodynamics import calculate_tm_with_chemistry

        # LNA and gapmer share the RNA/DNA hybrid NN base parameters
        lna_tm = calculate_tm_with_chemistry(self.SEQ, chemistry="lna")
        gapmer_tm = calculate_tm_with_chemistry(self.SEQ, chemistry="gapmer")
        assert lna_tm == pytest.approx(gapmer_tm)

    def test_lna_preset_loads(self):
        from aso_designer.thermodynamics import get_chemistry_tm_adjustment

        assert get_chemistry_tm_adjustment("lna") == 0.0

    def test_summary_threads_pattern(self):
        from aso_designer.thermodynamics import compute_thermodynamic_summary

        plain = compute_thermodynamic_summary(self.SEQ, chemistry="lna")
        adjusted = compute_thermodynamic_summary(
            self.SEQ, chemistry="lna", modification_pattern=self.PATTERN
        )
        assert adjusted.tm == pytest.approx(plain.tm + 3.5 * 5)
        assert adjusted.dg_duplex == pytest.approx(plain.dg_duplex - 1.4 * 5)


# ---------------------------------------------------------------------------
# PS/PO backbone pattern
# ---------------------------------------------------------------------------


class TestBackbonePattern:
    def test_gapmer_full_ps_default(self):
        design = GapmerDesign(gap_size=10, wing_length=5)
        backbone = design.build_backbone(20)
        assert backbone == ["PS"] * 19

    def test_gapmer_po_option(self):
        design = GapmerDesign(gap_size=10, wing_length=5, backbone="po")
        assert design.build_backbone(20) == ["PO"] * 19

    def test_gapmer_invalid_backbone(self):
        design = GapmerDesign(gap_size=10, wing_length=5, backbone="xyz")
        with pytest.raises(ConfigurationError, match="Invalid gapmer backbone"):
            design.build_backbone(20)

    def test_other_chemistries_po(self):
        assert MixmerDesign().build_backbone(10) == ["PO"] * 9
        assert LnaDesign().build_backbone(10) == ["PO"] * 9
        assert FullChemistryDesign().build_backbone(10) == ["PO"] * 9

    def test_apply_chemistry_sets_backbone(self):
        c = ASOCandidate(sequence="A" * 20, start=0, end=20, length=20)
        apply_chemistry([c], GapmerDesign())
        assert c.backbone_pattern == ["PS"] * 19

    def test_backbone_label(self):
        from aso_designer.chemistry import backbone_label

        assert backbone_label([]) == "none"
        assert backbone_label(["PS"] * 5) == "full-PS"
        assert backbone_label(["PO"] * 5) == "full-PO"
        assert backbone_label(["PS", "PO", "PS"]) == "mixed"

    def test_candidate_backbone_roundtrip(self):
        c = ASOCandidate(sequence="A" * 20, start=0, end=20, length=20)
        c.backbone_pattern = ["PS"] * 19
        data = c.to_dict()
        restored = ASOCandidate.from_dict(data)
        assert restored.backbone_pattern == ["PS"] * 19

    def test_gapmer_backbone_config(self):
        from aso_designer.config_schema import GapmerConfig

        assert GapmerConfig().backbone == "ps"
        assert GapmerConfig(backbone="po").backbone == "po"
