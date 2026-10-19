"""Tests for the ASO self-dimer / self-complementarity screen."""

from unittest.mock import MagicMock, patch

from aso_designer.candidate import ASOCandidate
from aso_designer.dimer import (
    DimerScreenResult,
    find_3p_complementarity,
    run_self_duplex,
    screen_self_dimer,
)
from aso_designer.filters import filter_self_dimer

# Clean sequence: its 3' tail ends in AA, so any paired internal window
# would have to start with TT, which never occurs in the sequence.
CLEAN_SEQ = "CGCGCGCGCAAA"
# CATG is a reverse palindrome: two copies pair via their 3' tails.
DIMMER_SEQ = "AAACATG"


def _candidate(seq: str) -> ASOCandidate:
    return ASOCandidate(
        sequence=seq,
        start=0,
        end=len(seq),
        length=len(seq),
        region="CDS",
    )


# ---------------------------------------------------------------------------
# find_3p_complementarity
# ---------------------------------------------------------------------------


class TestFind3pComplementarity:
    def test_clean_sequence(self):
        assert find_3p_complementarity(CLEAN_SEQ, min_stretch=4) is None

    def test_internal_match(self):
        # 3' tail "TTTT" pairs with the leading "AAAA" window
        assert find_3p_complementarity("AAAAACCCCCAAATTTT", min_stretch=4) == (4, 0)

    def test_palindromic_tail_pairs_tail_to_tail(self):
        assert find_3p_complementarity(DIMMER_SEQ, min_stretch=4) == (4, 3)

    def test_full_palindrome(self):
        assert find_3p_complementarity("ACGT", min_stretch=4) == (4, 0)

    def test_wobble_pairs(self):
        # 3' tail GGGG pairs the internal TTTT window only via G·T wobble
        assert find_3p_complementarity("TTTTAAAAAAGGGG", min_stretch=4) == (4, 0)
        assert find_3p_complementarity("TTTTAAAAAAGGGG", min_stretch=4, allow_wobble=False) is None

    def test_min_stretch_threshold(self):
        # "ATTCAAAACCAT" has a longest 3'-anchored stretch of 2 nt;
        # "GGTCTCTCCGCC" has one of 3 nt.
        assert find_3p_complementarity("ATTCAAAACCAT", min_stretch=4) is None
        assert find_3p_complementarity("ATTCAAAACCAT", min_stretch=2) is not None
        assert find_3p_complementarity("GGTCTCTCCGCC", min_stretch=4) is None
        assert find_3p_complementarity("GGTCTCTCCGCC", min_stretch=3) is not None

    def test_rna_input_normalized(self):
        # U treated as T
        assert find_3p_complementarity("AAAAUUUU", min_stretch=4) is not None

    def test_short_sequence(self):
        assert find_3p_complementarity("A", min_stretch=4) is None
        assert find_3p_complementarity("", min_stretch=4) is None


# ---------------------------------------------------------------------------
# screen_self_dimer
# ---------------------------------------------------------------------------


class TestScreenSelfDimer:
    def test_flags_complementary(self):
        result = screen_self_dimer(DIMMER_SEQ, min_stretch=4)
        assert isinstance(result, DimerScreenResult)
        assert result.flagged
        assert result.three_prime_stretch == 4
        assert result.three_prime_position == 3
        assert "self-complementary" in result.detail

    def test_passes_clean_sequence(self):
        result = screen_self_dimer(CLEAN_SEQ, min_stretch=4)
        assert not result.flagged
        assert result.three_prime_stretch == 0
        assert result.three_prime_position is None
        assert result.homodimer_dg is None

    def test_zero_min_stretch_disables(self):
        result = screen_self_dimer(DIMMER_SEQ, min_stretch=0)
        assert not result.flagged

    @patch("aso_designer.viennarna.run_rnaduplex")
    def test_dg_threshold_flags(self, mock_duplex):
        mock_duplex.return_value = MagicMock(mfe=-7.5)
        result = screen_self_dimer(CLEAN_SEQ, min_stretch=4, compute_dg=True)
        # No sequence-level stretch, but the homodimer MFE trips the threshold
        assert result.flagged
        assert result.homodimer_dg == -7.5
        assert "homodimer MFE" in result.detail

    @patch("aso_designer.viennarna.run_rnaduplex")
    def test_dg_above_threshold_not_flagged(self, mock_duplex):
        mock_duplex.return_value = MagicMock(mfe=-2.0)
        result = screen_self_dimer(CLEAN_SEQ, min_stretch=4, compute_dg=True)
        assert not result.flagged

    @patch("aso_designer.viennarna.run_rnaduplex")
    def test_run_self_duplex_returns_mfe(self, mock_duplex):
        mock_duplex.return_value = MagicMock(mfe=-5.25)
        assert run_self_duplex("ACGTACGT") == -5.25
        # both arguments are the same sequence (homodimer)
        assert mock_duplex.call_args.args == ("ACGTACGT", "ACGTACGT")


# ---------------------------------------------------------------------------
# Filter integration
# ---------------------------------------------------------------------------


class TestFilterSelfDimer:
    def test_flags_dimer_risk(self):
        candidates = [_candidate(DIMMER_SEQ), _candidate(CLEAN_SEQ)]
        passed = filter_self_dimer(candidates, min_stretch=4)
        assert passed == 1
        assert not candidates[0].passes_filters
        assert candidates[0].features["self_dimer_stretch"] == 4
        assert "self-complementary" in candidates[0].features["filter_reason"]
        assert candidates[1].passes_filters
        assert candidates[1].features["self_dimer_stretch"] == 0

    def test_zero_disables(self):
        candidates = [_candidate(DIMMER_SEQ)]
        passed = filter_self_dimer(candidates, min_stretch=0)
        assert passed == 1
        assert candidates[0].passes_filters
        assert candidates[0].features["self_dimer_stretch"] == 4

    def test_report_only_by_default(self):
        """Default min_stretch=0 records the stretch without rejecting."""
        candidates = [_candidate(DIMMER_SEQ), _candidate(CLEAN_SEQ)]
        passed = filter_self_dimer(candidates)
        assert passed == 2
        for c in candidates:
            assert c.passes_filters
        assert candidates[0].features["self_dimer_stretch"] == 4
        assert candidates[1].features["self_dimer_stretch"] == 0

    def test_skips_already_failed(self):
        c = _candidate(DIMMER_SEQ)
        c.passes_filters = False
        assert filter_self_dimer([c], min_stretch=4) == 0
        assert "self_dimer_stretch" not in c.features
