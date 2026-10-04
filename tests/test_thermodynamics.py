"""Tests for thermodynamic calculations (Tm, ΔG, self-structure, net binding)."""

from unittest.mock import patch

import pytest

from aso_designer.thermodynamics import (
    _RNA_NN,
    DuplexResult,
    SelfStructureResult,
    TargetStructureResult,
    ThermodynamicSummary,
    _is_self_complementary,
    calculate_dg,
    calculate_duplex_dg,
    calculate_duplex_dg_nn,
    calculate_duplex_dg_rnaduplex,
    calculate_net_binding_dg,
    calculate_nn_thermodynamics,
    calculate_self_structure_dg,
    calculate_self_structure_dg_nn,
    calculate_self_structure_penalty,
    calculate_target_structure_dg,
    calculate_target_structure_dg_nn,
    calculate_tm,
    calculate_tm_with_chemistry,
    classify_binding,
    compute_thermodynamic_summary,
    get_chemistry_self_structure_weight,
    get_chemistry_tm_adjustment,
    score_tm,
)
from aso_designer.utils.exceptions import ThermodynamicsError
from aso_designer.viennarna import RNAduplexResult, RNAfoldResult

# ---------------------------------------------------------------------------
# Nearest-neighbor thermodynamics
# ---------------------------------------------------------------------------


class TestCalculateNNThermodynamics:
    def test_dna_sequence(self):
        dh, ds = calculate_nn_thermodynamics("GCGCGCGC", rna=False)
        assert dh < 0  # Exothermic
        assert ds < 0  # Entropy decreases

    def test_rna_sequence(self):
        dh, ds = calculate_nn_thermodynamics("GCGCGCGC", rna=True)
        assert dh < 0
        assert ds < 0

    def test_turner2004_rna_table_values(self):
        """The 10 unique Turner 2004 (Xia et al. 1998) stacking parameters
        must match the NNDB-published ΔH/ΔS values exactly.

        Source: https://rna.urmc.rochester.edu/NNDB/turner04
        (ΔS derived from published ΔG°37 and ΔH° at 310.15 K).
        """
        expected_unique = {
            "AA": (-6.82, -19.0),
            "AU": (-9.38, -26.7),
            "UA": (-7.69, -20.5),
            "CU": (-10.48, -27.1),
            "CA": (-10.44, -26.9),
            "GU": (-11.40, -29.5),
            "GA": (-12.44, -32.5),
            "CG": (-10.64, -26.7),
            "GG": (-13.39, -32.7),
            "GC": (-14.88, -36.9),
        }
        for pair, (dh, ds) in expected_unique.items():
            actual = _RNA_NN[pair]
            assert actual[0] == pytest.approx(dh, abs=1e-9), f"{pair} ΔH"
            assert actual[1] == pytest.approx(ds, abs=1e-9), f"{pair} ΔS"

    def test_turner2004_reverse_complement_symmetry(self):
        """Symmetric expansions must equal their unique parameter."""
        symmetry_pairs = {
            "AA": "UU",
            "AC": "GU",
            "AG": "CU",
            "CA": "UG",
            "CC": "GG",
            "GA": "UC",
        }
        for a, b in symmetry_pairs.items():
            assert _RNA_NN[a] == _RNA_NN[b], f"{a} and {b} must share parameters"

    def test_turner2004_gc_terminal_initiation(self):
        """GC-terminal RNA duplex: stacks + intermolecular initiation +
        symmetry (self-complementary), no AU end penalty.

        GCGCGCGC: 4×GC + 3×CG stacks, both terminals GC, self-complementary.
        ΔH = 4(-14.88) + 3(-10.64) + 3.61 = -87.83
        ΔS = 4(-36.9) + 3(-26.7) + (-1.5) + (-1.4) = -230.6
        """
        dh, ds = calculate_nn_thermodynamics("GCGCGCGC", rna=True)
        assert dh == pytest.approx(-87.83, abs=1e-9)
        assert ds == pytest.approx(-230.6, abs=1e-9)

    def test_turner2004_au_terminal_initiation(self):
        """AU-terminal RNA duplex: stacks + intermolecular initiation +
        one AU end penalty, non-self-complementary.

        AUGCAUGC: AU×2, UG×2, GC×2, CA×1 stacks; terminals A (AU) and C.
        ΔH = -18.76 - 20.88 - 29.76 - 10.44 + 3.61 + 3.72 = -72.51
        ΔS = -53.4 - 53.8 - 73.8 - 26.9 - 1.5 + 10.5 = -198.9
        """
        dh, ds = calculate_nn_thermodynamics("AUGCAUGC", rna=True)
        assert dh == pytest.approx(-72.51, abs=1e-9)
        assert ds == pytest.approx(-198.9, abs=1e-9)

    def test_dna_vs_rna_different(self):
        dh_dna, _ = calculate_nn_thermodynamics("GCGCGCGC", rna=False)
        dh_rna, _ = calculate_nn_thermodynamics("GCGCGCGC", rna=True)
        assert dh_dna != dh_rna

    def test_all_dna_dinucleotides(self):
        """Every DNA dinucleotide should be in the parameter table."""
        for a in "ACGT":
            for b in "ACGT":
                seq = a + b + "ATGCATGC"  # Ensure length >= 2
                dh, ds = calculate_nn_thermodynamics(seq, rna=False)
                assert isinstance(dh, float)
                assert isinstance(ds, float)

    def test_all_rna_dinucleotides(self):
        """Every RNA dinucleotide should be in the parameter table."""
        for a in "ACGU":
            for b in "ACGU":
                seq = a + b + "AUGCAUGC"
                dh, ds = calculate_nn_thermodynamics(seq, rna=True)
                assert isinstance(dh, float)

    def test_too_short(self):
        with pytest.raises(ThermodynamicsError, match="too short"):
            calculate_nn_thermodynamics("A", rna=False)

    def test_self_complementary_dna(self):
        # A palindromic DNA sequence
        assert _is_self_complementary("GAATTC", rna=False)
        assert not _is_self_complementary("GAATT", rna=False)

    def test_self_complementary_rna(self):
        assert _is_self_complementary("GAAUUC", rna=True)
        assert not _is_self_complementary("GAAUU", rna=True)

    def test_rna_input_with_t(self):
        """DNA input with T should work for RNA calculation (T→U)."""
        dh, ds = calculate_nn_thermodynamics("ATGCATGC", rna=True)
        # Should convert T to U internally
        assert dh < 0


# ---------------------------------------------------------------------------
# Tm calculation
# ---------------------------------------------------------------------------


class TestCalculateTm:
    def test_basic_dna(self):
        tm = calculate_tm("GCGCGCGC", na_conc=0.1, oligo_conc=1e-6, rna=False)
        assert 30 < tm < 100

    def test_basic_rna(self):
        tm = calculate_tm("GCGCGCGC", na_conc=0.1, oligo_conc=1e-6, rna=True)
        assert 30 < tm < 100

    def test_rna_higher_than_dna(self):
        """RNA duplexes typically have higher Tm than DNA."""
        tm_dna = calculate_tm("GCGCGCGC", rna=False)
        tm_rna = calculate_tm("GCGCGCGC", rna=True)
        assert tm_rna > tm_dna

    def test_longer_sequence_higher_tm(self):
        tm_short = calculate_tm("GCGCGCGC", rna=False)
        tm_long = calculate_tm("GCGCGCGCGCGCGCGC", rna=False)
        assert tm_long > tm_short

    def test_salt_correction_santalucia(self):
        tm_low_salt = calculate_tm("GCGCGCGC", na_conc=0.01, rna=False)
        tm_high_salt = calculate_tm("GCGCGCGC", na_conc=1.0, rna=False)
        assert tm_high_salt > tm_low_salt

    def test_salt_correction_log10(self):
        tm = calculate_tm("GCGCGCGC", na_conc=0.1, rna=False, salt_correction="log10")
        assert isinstance(tm, float)

    def test_invalid_salt_correction(self):
        with pytest.raises(ThermodynamicsError, match="Unknown salt correction"):
            calculate_tm("GCGCGCGC", salt_correction="invalid")

    def test_self_complementary_symmetry(self):
        """Self-complementary sequences use x=1, not x=4."""
        tm_sc = calculate_tm("GAATTC", rna=False)
        tm_ns = calculate_tm("GAATTA", rna=False)
        # Both should produce valid Tm values
        assert isinstance(tm_sc, float)
        assert isinstance(tm_ns, float)


# ---------------------------------------------------------------------------
# ΔG calculation
# ---------------------------------------------------------------------------


class TestCalculateDg:
    def test_basic_dna(self):
        dg = calculate_dg("GCGCGCGC", temperature=37.0, rna=False)
        assert dg < 0  # Favorable

    def test_basic_rna(self):
        dg = calculate_dg("GCGCGCGC", temperature=37.0, rna=True)
        assert dg < 0

    def test_temperature_effect(self):
        dg_low = calculate_dg("GCGCGCGC", temperature=10.0, rna=False)
        dg_high = calculate_dg("GCGCGCGC", temperature=80.0, rna=False)
        # At higher temperature, ΔG should be less favorable (less negative)
        assert dg_high > dg_low


# ---------------------------------------------------------------------------
# Chemistry-aware Tm
# ---------------------------------------------------------------------------


class TestChemistryTm:
    def test_gapmer_adjustment(self):
        assert get_chemistry_tm_adjustment("gapmer") == 0.0

    def test_2ome_adjustment(self):
        assert get_chemistry_tm_adjustment("2ome") == 5.0

    def test_moe_adjustment(self):
        assert get_chemistry_tm_adjustment("moe") == 6.0

    def test_morpholino_adjustment(self):
        assert get_chemistry_tm_adjustment("morpholino") == -3.0

    def test_unknown_chemistry_adjustment(self):
        assert get_chemistry_tm_adjustment("unknown") == 0.0

    def test_2ome_higher_than_gapmer(self):
        tm_gapmer = calculate_tm_with_chemistry("ATGCATGCATGCATGC", chemistry="gapmer")
        tm_2ome = calculate_tm_with_chemistry("ATGCATGCATGCATGC", chemistry="2ome")
        assert tm_2ome > tm_gapmer

    def test_morpholino_lower_than_gapmer(self):
        tm_gapmer = calculate_tm_with_chemistry("ATGCATGCATGCATGC", chemistry="gapmer")
        tm_morph = calculate_tm_with_chemistry("ATGCATGCATGCATGC", chemistry="morpholino")
        assert tm_morph < tm_gapmer


# ---------------------------------------------------------------------------
# Tm scoring
# ---------------------------------------------------------------------------


class TestScoreTm:
    def test_optimal(self):
        assert score_tm(52.0, tm_optimal=52.0, tm_min=40.0, tm_max=65.0) == 1.0

    def test_below_min(self):
        assert score_tm(35.0, tm_optimal=52.0, tm_min=40.0, tm_max=65.0) == 0.0

    def test_above_max(self):
        assert score_tm(70.0, tm_optimal=52.0, tm_min=40.0, tm_max=65.0) == 0.0

    def test_at_min(self):
        assert score_tm(40.0, tm_optimal=52.0, tm_min=40.0, tm_max=65.0) == 0.5

    def test_at_max(self):
        assert score_tm(65.0, tm_optimal=52.0, tm_min=40.0, tm_max=65.0) == 0.5

    def test_between_optimal_and_max(self):
        score = score_tm(58.0, tm_optimal=52.0, tm_min=40.0, tm_max=65.0)
        assert 0.5 < score < 1.0

    def test_between_min_and_optimal(self):
        score = score_tm(46.0, tm_optimal=52.0, tm_min=40.0, tm_max=65.0)
        assert 0.5 < score < 1.0


# ---------------------------------------------------------------------------
# Duplex ΔG
# ---------------------------------------------------------------------------


class TestDuplexDg:
    def test_nn_gapmer(self):
        r = calculate_duplex_dg_nn("ATGCATGCATGCATGC", chemistry="gapmer")
        assert r.dg_duplex < 0
        assert r.method == "nearest_neighbor"

    def test_nn_2ome(self):
        r = calculate_duplex_dg_nn("ATGCATGCATGCATGC", chemistry="2ome")
        assert r.dg_duplex < 0
        # 2'OMe uses RNA params → more negative
        r_gapmer = calculate_duplex_dg_nn("ATGCATGCATGCATGC", chemistry="gapmer")
        assert r.dg_duplex < r_gapmer.dg_duplex  # RNA more stable

    def test_convenience_nn(self):
        r = calculate_duplex_dg("ATGCATGCATGCATGC", method="nearest_neighbor")
        assert isinstance(r, DuplexResult)
        assert r.dg_duplex < 0

    @patch("aso_designer.viennarna.run_rnaduplex")
    def test_rnaduplex_mock(self, mock_rnaduplex):
        mock_rnaduplex.return_value = RNAduplexResult(
            sequence1="AUGCAUGCAUGCAUGC",
            sequence2="AUGCAUGCAUGCAUGC",
            structure=".((((...))))...&...((((...)))).",
            mfe=-12.5,
            start1=1,
            end1=15,
            start2=1,
            end2=15,
        )
        r = calculate_duplex_dg_rnaduplex("ATGCATGCATGCATGC", "ATGCATGCATGCATGC")
        assert r.dg_duplex == -12.5
        assert r.method == "rnaduplex"
        assert r.start_target == 1
        assert r.end_target == 15

    @patch("aso_designer.viennarna.run_rnaduplex")
    def test_rnaduplex_failure(self, mock_rnaduplex):
        mock_rnaduplex.side_effect = Exception("tool not found")
        with pytest.raises(ThermodynamicsError, match="RNAduplex failed"):
            calculate_duplex_dg_rnaduplex("ATGC", "ATGC")

    def test_rnaduplex_without_target(self):
        with pytest.raises(ThermodynamicsError, match="target_sequence is required"):
            calculate_duplex_dg("ATGC", method="rnaduplex")

    def test_unknown_method(self):
        with pytest.raises(ThermodynamicsError, match="Unknown duplex method"):
            calculate_duplex_dg("ATGC", method="invalid")


# ---------------------------------------------------------------------------
# Self-structure ΔG
# ---------------------------------------------------------------------------


class TestSelfStructure:
    def test_nn_basic(self):
        r = calculate_self_structure_dg_nn("ATGCATGCATGCATGC", chemistry="gapmer")
        assert r.dg_self < 0
        assert r.method == "nearest_neighbor"
        assert r.structure == ""

    def test_nn_fraction_default_03(self):
        """Default NN fraction scales the perfect-duplex ΔG by 0.3."""
        dg_duplex = calculate_dg("ATGCATGCATGCATGC", temperature=37.0, rna=False)
        r = calculate_self_structure_dg_nn("ATGCATGCATGCATGC", chemistry="gapmer")
        assert r.dg_self == pytest.approx(dg_duplex * 0.3)

    def test_nn_fraction_custom(self):
        """A custom nn_fraction scales linearly."""
        seq = "ATGCATGCATGCATGC"
        r_default = calculate_self_structure_dg_nn(seq, chemistry="gapmer")
        r_half = calculate_self_structure_dg_nn(seq, chemistry="gapmer", nn_fraction=0.6)
        assert r_half.dg_self == pytest.approx(r_default.dg_self * 2.0)

    def test_nn_fraction_threaded_through_summary(self):
        """compute_thermodynamic_summary forwards nn_fraction."""
        seq = "ATGCATGCATGCATGCATGC"
        s_default = compute_thermodynamic_summary(seq, nn_fraction=0.3)
        s_high = compute_thermodynamic_summary(seq, nn_fraction=0.9)
        assert s_high.dg_self == pytest.approx(s_default.dg_self * 3.0)

    def test_convenience_nn(self):
        r = calculate_self_structure_dg("ATGCATGCATGCATGC", chemistry="gapmer")
        assert isinstance(r, SelfStructureResult)
        assert r.penalty >= 0

    @patch("aso_designer.viennarna.run_rnafold")
    def test_rnafold_mock(self, mock_rnafold):
        mock_rnafold.return_value = RNAfoldResult(
            sequence="AUGCAUGCAUGCAUGC",
            structure="((((....))))....",
            mfe=-3.5,
        )
        r = calculate_self_structure_dg("ATGCATGCATGCATGC", chemistry="gapmer", method="rnafold")
        assert r.dg_self == -3.5
        assert r.method == "rnafold"
        assert r.structure == "((((....))))...."
        assert r.penalty > 0

    @patch("aso_designer.viennarna.run_rnafold")
    def test_rnafold_failure(self, mock_rnafold):
        mock_rnafold.side_effect = Exception("tool not found")
        with pytest.raises(ThermodynamicsError, match="RNAfold failed"):
            calculate_self_structure_dg("ATGCATGCATGCATGC", method="rnafold")

    def test_unknown_method(self):
        with pytest.raises(ThermodynamicsError, match="Unknown self-structure method"):
            calculate_self_structure_dg("ATGC", method="invalid")


class TestSelfStructurePenalty:
    def test_no_self_structure(self):
        assert calculate_self_structure_penalty(0.0) == 0.0
        assert calculate_self_structure_penalty(1.0) == 0.0

    def test_at_threshold(self):
        assert calculate_self_structure_penalty(-3.0) == 1.0

    def test_half_threshold(self):
        assert calculate_self_structure_penalty(-1.5) == 0.5

    def test_beyond_threshold(self):
        assert calculate_self_structure_penalty(-6.0) == 2.0

    def test_custom_threshold(self):
        assert calculate_self_structure_penalty(-5.0, self_structure_max=-5.0) == 1.0


class TestChemistrySelfStructureWeight:
    def test_gapmer(self):
        assert get_chemistry_self_structure_weight("gapmer") == 0.8

    def test_2ome(self):
        assert get_chemistry_self_structure_weight("2ome") == 1.2

    def test_mixmer(self):
        assert get_chemistry_self_structure_weight("mixmer") == 1.0

    def test_unknown(self):
        assert get_chemistry_self_structure_weight("unknown") == 1.0

    def test_gapmer_reduces_penalty(self):
        """Gapmer chemistry should reduce the self-structure penalty."""
        r_gapmer = calculate_self_structure_dg("GCGCGCGCGCGCGCGC", chemistry="gapmer")
        r_mixmer = calculate_self_structure_dg("GCGCGCGCGCGCGCGC", chemistry="mixmer")
        assert r_gapmer.penalty < r_mixmer.penalty


# ---------------------------------------------------------------------------
# Target structure ΔG
# ---------------------------------------------------------------------------


class TestTargetStructure:
    def test_nn_basic(self):
        r = calculate_target_structure_dg_nn("ATGCATGCATGCATGC")
        assert r.dg_target < 0
        assert r.method == "nearest_neighbor"
        assert r.structure == ""

    def test_convenience_nn(self):
        r = calculate_target_structure_dg("ATGCATGCATGCATGC")
        assert isinstance(r, TargetStructureResult)

    @patch("aso_designer.viennarna.run_rnafold")
    def test_rnafold_mock(self, mock_rnafold):
        mock_rnafold.return_value = RNAfoldResult(
            sequence="AUGCAUGCAUGCAUGC",
            structure="((((....))))....",
            mfe=-4.2,
        )
        r = calculate_target_structure_dg("ATGCATGCATGCATGC", method="rnafold")
        assert r.dg_target == -4.2
        assert r.method == "rnafold"

    def test_unknown_method(self):
        with pytest.raises(ThermodynamicsError, match="Unknown target structure method"):
            calculate_target_structure_dg("ATGC", method="invalid")


# ---------------------------------------------------------------------------
# Net binding energy
# ---------------------------------------------------------------------------


class TestNetBindingDg:
    def test_basic(self):
        dg = calculate_net_binding_dg(-15.0, -3.0, -2.0)
        assert dg == -10.0

    def test_no_self_structure(self):
        dg = calculate_net_binding_dg(-10.0, 0.0, 0.0)
        assert dg == -10.0

    def test_all_unfavorable(self):
        dg = calculate_net_binding_dg(-5.0, -5.0, -5.0)
        assert dg == 5.0  # Positive = unfavorable


class TestClassifyBinding:
    def test_strong(self):
        assert classify_binding(-12.0) == "strong"

    def test_at_strong_threshold(self):
        assert classify_binding(-10.0) == "strong"

    def test_moderate(self):
        assert classify_binding(-7.0) == "moderate"

    def test_at_moderate_threshold(self):
        assert classify_binding(-5.0) == "moderate"

    def test_weak(self):
        assert classify_binding(-3.0) == "weak"

    def test_positive(self):
        assert classify_binding(5.0) == "weak"

    def test_custom_thresholds(self):
        assert classify_binding(-8.0, strong_threshold=-15.0, moderate_threshold=-7.0) == "moderate"


# ---------------------------------------------------------------------------
# Thermodynamic summary
# ---------------------------------------------------------------------------


class TestThermodynamicSummary:
    def test_nn_summary_with_target(self):
        s = compute_thermodynamic_summary(
            "ATGCATGCATGCATGC",
            target_sequence="GCATGCATGCATGCAT",
            chemistry="gapmer",
            method="nearest_neighbor",
        )
        assert isinstance(s, ThermodynamicSummary)
        assert 30 < s.tm < 100
        assert 0 <= s.tm_score <= 1
        assert s.dg_duplex < 0
        assert s.dg_self < 0
        assert s.dg_target < 0
        assert isinstance(s.binding_class, str)

    def test_nn_summary_without_target(self):
        s = compute_thermodynamic_summary(
            "ATGCATGCATGCATGC",
            chemistry="gapmer",
            method="nearest_neighbor",
        )
        assert s.dg_target == 0.0
        # Without target structure cost, net should be more favorable
        s_with = compute_thermodynamic_summary(
            "ATGCATGCATGCATGC",
            target_sequence="ATGCATGCATGCATGC",
            chemistry="gapmer",
            method="nearest_neighbor",
        )
        assert s.dg_net < s_with.dg_net  # No target cost = more favorable

    @patch("aso_designer.viennarna.run_rnafold")
    def test_rnafold_summary(self, mock_rnafold):
        mock_rnafold.return_value = RNAfoldResult(
            sequence="AUGCAUGCAUGCAUGC",
            structure="((((....))))....",
            mfe=-3.5,
        )
        s = compute_thermodynamic_summary(
            "ATGCATGCATGCATGC",
            target_sequence="ATGCATGCATGCATGC",
            chemistry="gapmer",
            method="rnafold",
        )
        assert s.dg_self == -3.5
        assert s.dg_target == -3.5  # Same mock for both calls
        assert s.binding_class in ("strong", "moderate", "weak")

    @patch("aso_designer.viennarna.run_rnaduplex")
    @patch("aso_designer.viennarna.run_rnafold")
    def test_rnaduplex_summary(self, mock_rnafold, mock_rnaduplex):
        mock_rnaduplex.return_value = RNAduplexResult(
            sequence1="AUGCAUGCAUGCAUGC",
            sequence2="AUGCAUGCAUGCAUGC",
            structure=".((((...))))...&...((((...)))).",
            mfe=-15.0,
            start1=1,
            end1=15,
            start2=1,
            end2=15,
        )
        mock_rnafold.return_value = RNAfoldResult(
            sequence="AUGCAUGCAUGCAUGC",
            structure="((((....))))....",
            mfe=-2.0,
        )
        s = compute_thermodynamic_summary(
            "ATGCATGCATGCATGC",
            target_sequence="ATGCATGCATGCATGC",
            chemistry="gapmer",
            method="rnaduplex",
        )
        assert s.dg_duplex == -15.0
        assert s.dg_self == -2.0  # RNAfold for self-structure
        assert s.dg_target == -2.0  # RNAfold for target

    def test_2ome_uses_rna_params(self):
        """2'OMe chemistry should use RNA NN parameters for Tm."""
        s_gapmer = compute_thermodynamic_summary(
            "ATGCATGCATGCATGC",
            chemistry="gapmer",
            method="nearest_neighbor",
        )
        s_2ome = compute_thermodynamic_summary(
            "ATGCATGCATGCATGC",
            chemistry="2ome",
            method="nearest_neighbor",
        )
        # 2'OMe should have higher Tm (RNA params + adjustment)
        assert s_2ome.tm > s_gapmer.tm

    def test_binding_class_values(self):
        s = compute_thermodynamic_summary(
            "ATGCATGCATGCATGC",
            chemistry="gapmer",
            method="nearest_neighbor",
        )
        assert s.binding_class in ("strong", "moderate", "weak")
