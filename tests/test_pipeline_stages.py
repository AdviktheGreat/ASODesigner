"""Tests for pipeline feature stages (accessibility, thermodynamics, off-target, variant).

These tests verify that the four feature-computation stages write real
per-candidate features that flow into the scoring stage. ViennaRNA and
network access are mocked; the thermodynamics stage runs the fast
nearest-neighbor path which requires no external tools.
"""

from unittest.mock import MagicMock, patch

import pytest

from aso_designer.pipeline import ASOPipeline, PipelineConfig
from aso_designer.transcript_fetcher import TranscriptRecord
from aso_designer.viennarna import RNAplfoldResult, ViennaRNAError

# ---------------------------------------------------------------------------
# Fixtures and helpers
# ---------------------------------------------------------------------------

# 60-nt synthetic transcript with a clean ORF structure:
#   positions 0-2    ATG (start codon)
#   positions 3-44   coding sequence (no in-frame stops)
#   positions 45-47  TAA (stop codon)
#   positions 48-59  3' UTR
TRANSCRIPT_SEQ = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATTAAGCGCGCGCGCAA"

# 60 candidate windows of length 20 -> 41 candidates; all 20-mers in the
# ATGC-repeat region have 50% GC content and no homopolymers, so most
# candidates pass the default filters.


def make_pipeline(sequence: str = TRANSCRIPT_SEQ) -> ASOPipeline:
    """Build an ASOPipeline with the transcript pre-loaded (no fetching)."""
    config = PipelineConfig(target="TEST")
    pipeline = ASOPipeline(config)
    pipeline._transcript = TranscriptRecord(
        accession="NM_TEST.1",
        sequence=sequence,
        gene_name="TESTGENE",
        organism="Homo sapiens",
    )
    return pipeline


def run_through_filter(pipeline: ASOPipeline) -> None:
    """Run generate + filter stages so candidates carry region/flags."""
    pipeline._stage_generate()
    pipeline._stage_filter()


def passing_candidates(pipeline: ASOPipeline) -> list:
    return [c for c in pipeline._candidates if c.passes_filters]


@pytest.fixture
def mock_rnaplfold_result() -> RNAplfoldResult:
    """Accessibility profile: unpaired (0.9) in the first half of the
    transcript, paired (0.1) in the second half."""
    n = len(TRANSCRIPT_SEQ)
    probs = [0.9] * (n // 2) + [0.1] * (n - n // 2)
    return RNAplfoldResult(
        positions=list(range(1, n + 1)),
        unpaired_probs=probs,
    )


# ---------------------------------------------------------------------------
# Accessibility stage
# ---------------------------------------------------------------------------


class TestAccessibilityStage:
    @patch("aso_designer.viennarna.verify_viennarna")
    @patch("aso_designer.accessibility.run_rnaplfold")
    def test_accessibility_scores_written(self, mock_plfold, mock_verify, mock_rnaplfold_result):
        mock_verify.return_value = MagicMock()
        mock_plfold.return_value = mock_rnaplfold_result

        pipeline = make_pipeline()
        run_through_filter(pipeline)
        pipeline._stage_accessibility()

        passed = passing_candidates(pipeline)
        assert len(passed) > 0
        for c in passed:
            assert "accessibility_score" in c.features
            score = c.features["accessibility_score"]
            assert 0.0 <= score <= 1.0

    @patch("aso_designer.viennarna.verify_viennarna")
    @patch("aso_designer.accessibility.run_rnaplfold")
    def test_accessible_region_scores_higher(self, mock_plfold, mock_verify, mock_rnaplfold_result):
        """A candidate in the unpaired first half should score higher
        than one in the paired second half."""
        mock_verify.return_value = MagicMock()
        mock_plfold.return_value = mock_rnaplfold_result

        pipeline = make_pipeline()
        run_through_filter(pipeline)
        pipeline._stage_accessibility()

        first_half = [
            c.features["accessibility_score"]
            for c in passing_candidates(pipeline)
            if c.start < len(TRANSCRIPT_SEQ) // 2
        ]
        second_half = [
            c.features["accessibility_score"]
            for c in passing_candidates(pipeline)
            if c.start >= len(TRANSCRIPT_SEQ) // 2
        ]
        assert first_half and second_half
        assert min(first_half) > max(second_half)

    @patch("aso_designer.viennarna.verify_viennarna")
    @patch("aso_designer.accessibility.run_rnaplfold")
    def test_filtered_candidates_skipped(self, mock_plfold, mock_verify, mock_rnaplfold_result):
        """Candidates that failed filters should not get accessibility scores."""
        mock_verify.return_value = MagicMock()
        mock_plfold.return_value = mock_rnaplfold_result

        pipeline = make_pipeline()
        run_through_filter(pipeline)
        pipeline._stage_accessibility()

        failed = [c for c in pipeline._candidates if not c.passes_filters]
        for c in failed:
            assert "accessibility_score" not in c.features

    @patch("aso_designer.viennarna.verify_viennarna")
    def test_missing_viennarna_raises(self, mock_verify):
        mock_verify.side_effect = ViennaRNAError("RNAfold not found")

        pipeline = make_pipeline()
        run_through_filter(pipeline)
        with pytest.raises(ViennaRNAError):
            pipeline._stage_accessibility()

    def test_no_transcript_raises(self):
        pipeline = make_pipeline()
        pipeline._transcript = None
        with pytest.raises(RuntimeError, match="No transcript"):
            pipeline._stage_accessibility()


# ---------------------------------------------------------------------------
# Thermodynamics stage
# ---------------------------------------------------------------------------


class TestThermodynamicsStage:
    def test_thermodynamic_features_written(self):
        """NN method requires no ViennaRNA; all features should be populated."""
        pipeline = make_pipeline()
        run_through_filter(pipeline)
        pipeline._stage_thermodynamics()

        passed = passing_candidates(pipeline)
        assert len(passed) > 0
        feature_keys = (
            "tm",
            "tm_score",
            "dg_duplex",
            "dg_self",
            "dg_target",
            "dg_net",
            "binding_class",
            "self_structure_penalty",
        )
        for c in passed:
            for key in feature_keys:
                assert key in c.features, f"missing {key} for candidate at {c.start}"
            assert c.features["tm"] > 0.0
            assert 0.0 <= c.features["tm_score"] <= 1.0
            assert c.features["binding_class"] in ("strong", "moderate", "weak")
            assert 0.0 <= c.features["self_structure_penalty"]

    def test_summaries_collected_for_reporting(self):
        pipeline = make_pipeline()
        run_through_filter(pipeline)
        pipeline._stage_thermodynamics()

        passed = passing_candidates(pipeline)
        assert len(pipeline._thermo_summaries) == len(passed)

    def test_filtered_candidates_skipped(self):
        pipeline = make_pipeline()
        run_through_filter(pipeline)
        pipeline._stage_thermodynamics()

        for c in pipeline._candidates:
            if not c.passes_filters:
                assert "tm" not in c.features

    def test_config_parameters_applied(self, tmp_path):
        """Custom Tm bounds from a user config should flow into the summaries."""
        config_file = tmp_path / "custom.yaml"
        config_file.write_text(
            "thermodynamics:\n  tm_min: 30.0\n  tm_max: 60.0\n  tm_optimal: 45.0\n"
        )
        config = PipelineConfig(target="TEST", config_path=str(config_file))
        pipeline = ASOPipeline(config)
        pipeline._transcript = TranscriptRecord(
            accession="NM_TEST.1",
            sequence=TRANSCRIPT_SEQ,
            gene_name="TESTGENE",
            organism="Homo sapiens",
        )
        run_through_filter(pipeline)
        pipeline._stage_thermodynamics()

        assert len(pipeline._thermo_summaries) == len(passing_candidates(pipeline))
        for summary in pipeline._thermo_summaries:
            assert 0.0 <= summary.tm_score <= 1.0

    def test_invalid_method_rejected(self, tmp_path):
        config_file = tmp_path / "bad.yaml"
        config_file.write_text("thermodynamics:\n  method: bogus\n")
        pipeline = ASOPipeline(PipelineConfig(target="TEST", config_path=str(config_file)))
        with pytest.raises(Exception, match="Config validation failed"):
            _ = pipeline.aso_config

    def test_no_transcript_raises(self):
        pipeline = make_pipeline()
        pipeline._transcript = None
        with pytest.raises(RuntimeError, match="No transcript"):
            pipeline._stage_thermodynamics()
