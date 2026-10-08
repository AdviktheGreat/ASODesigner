"""Tests for immunostimulatory-motif and G-quadruplex screening."""

import pytest

from aso_designer.motifs import (
    TLR9_MOTIFS,
    MotifReport,
    assess_motifs,
    cpg_count,
    cpg_fraction,
    find_g_quadruplex_motifs,
    longest_g_run,
)


class TestCpGCount:
    def test_counts_cg_dinucleotides(self):
        # Three CG dinucleotides.
        assert cpg_count("CGATCGATCG") == 3

    def test_no_cpg(self):
        assert cpg_count("ATATATAT") == 0

    def test_case_insensitive(self):
        assert cpg_count("cgATCgatCG") == 3

    def test_empty_and_short(self):
        assert cpg_count("") == 0
        assert cpg_count("C") == 0

    def test_overlapping_not_double_counted(self):
        # "CCG" contains one CG.
        assert cpg_count("CCG") == 1


class TestCpGFraction:
    def test_fraction_of_dinucleotide_positions(self):
        # 3 CpGs over 9 dinucleotide positions.
        assert cpg_fraction("CGATCGATCG") == pytest.approx(3 / 9)

    def test_no_cpg_is_zero(self):
        assert cpg_fraction("ATATAT") == 0.0

    def test_short_sequences(self):
        assert cpg_fraction("") == 0.0
        assert cpg_fraction("A") == 0.0

    def test_two_base_sequence(self):
        assert cpg_fraction("CG") == 1.0


class TestGQuadruplexMotifs:
    def test_canonical_motif_detected(self):
        # Four GGG tracts separated by 3-nt loops.
        seq = "GGGATTTGGGATTTGGGATTTGGG"
        assert find_g_quadruplex_motifs(seq) == [seq]

    def test_minimal_loops_detected(self):
        # Four GGG tracts separated by single bases.
        seq = "GGGAGGGAGGGAGGG"
        assert find_g_quadruplex_motifs(seq) == [seq]

    def test_two_motifs(self):
        seq = "GGGAGGGAGGGAGGG" + "TTTTTTTT" + "GGGAGGGAGGGAGGG"
        assert len(find_g_quadruplex_motifs(seq)) == 2

    def test_three_tracts_insufficient(self):
        # Only three G-tracts — no G-quadruplex motif.
        assert find_g_quadruplex_motifs("GGGAGGGAGGGA") == []

    def test_long_loop_breaks_motif(self):
        # Loop longer than 7 nt between tracts.
        assert (
            find_g_quadruplex_motifs("GGG" + "A" * 8 + "GGG" + "A" * 8 + "GGG" + "A" * 8 + "GGG")
            == []
        )

    def test_two_g_tracts_insufficient(self):
        assert find_g_quadruplex_motifs("GGGAGGGAGGGA") == []

    def test_no_guanines(self):
        assert find_g_quadruplex_motifs("ATATATAT") == []

    def test_case_insensitive(self):
        assert len(find_g_quadruplex_motifs("gggagggagggaggg")) == 1


class TestLongestGRun:
    def test_single_run(self):
        assert longest_g_run("ATGGGCA") == 3

    def test_longest_wins(self):
        assert longest_g_run("GGAGGGGAGG") == 4

    def test_no_guanines(self):
        assert longest_g_run("ATAT") == 0

    def test_empty(self):
        assert longest_g_run("") == 0

    def test_case_insensitive(self):
        assert longest_g_run("gggg") == 4


class TestAssessMotifs:
    def test_report_fields(self):
        report = assess_motifs("CGATCGATCG")
        assert report.sequence == "CGATCGATCG"
        assert report.cpg_count == 3
        assert report.cpg_fraction == pytest.approx(3 / 9)
        assert report.g_quadruplex_motifs == 0
        assert report.longest_g_run == 1
        assert report.tlr9_stimulatory is False

    def test_tlr9_human_motif(self):
        report = assess_motifs("AATCGTCGAA")
        assert report.tlr9_stimulatory is True
        assert report.risk_label() == "elevated"

    def test_tlr9_mouse_motif(self):
        report = assess_motifs("AAGTCGTTAA")
        assert report.tlr9_stimulatory is True

    def test_all_canonical_motifs_listed(self):
        assert TLR9_MOTIFS == ("TCGTCG", "GTCGTT")

    def test_g_quadruplex_is_elevated(self):
        report = assess_motifs("GGGAGGGAGGGAGGG")
        assert report.g_quadruplex_motifs == 1
        assert report.risk_label() == "elevated"

    def test_cpg_rich_without_flagged_motif_is_moderate(self):
        # High CpG density but no TLR9 hexamer and no G4 motif.
        report = assess_motifs("CGCGCGCGCG")
        assert report.cpg_fraction > 0.10
        assert report.tlr9_stimulatory is False
        assert report.risk_label() == "moderate"

    def test_clean_sequence_is_low(self):
        report = assess_motifs("ATATATATAT")
        assert report.risk_label() == "low"

    def test_elevated_outranks_moderate(self):
        # TLR9 motif present and CpG-rich: elevated wins.
        assert assess_motifs("TCGTCGCG").risk_label() == "elevated"

    def test_short_sequence(self):
        report = assess_motifs("CG")
        assert report.cpg_count == 1
        assert report.cpg_fraction == pytest.approx(1.0)
        assert report.risk_label() == "moderate"

    def test_returns_motif_report(self):
        assert isinstance(assess_motifs("ATGC"), MotifReport)


class TestRiskLabel:
    def test_labels(self):
        def make(cpg_fraction=0.0, tlr9=False, g4=0):
            return MotifReport(
                sequence="",
                cpg_count=0,
                cpg_fraction=cpg_fraction,
                g_quadruplex_motifs=g4,
                longest_g_run=0,
                tlr9_stimulatory=tlr9,
            )

        assert make().risk_label() == "low"
        assert make(cpg_fraction=0.11).risk_label() == "moderate"
        assert make(cpg_fraction=0.10).risk_label() == "low"  # at threshold, not above
        assert make(tlr9=True).risk_label() == "elevated"
        assert make(g4=1).risk_label() == "elevated"
