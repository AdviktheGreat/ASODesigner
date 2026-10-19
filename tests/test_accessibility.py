"""Tests for target RNA structure prediction and accessibility scoring."""

import json
from unittest.mock import patch

import pytest
from pydantic import ValidationError

from aso_designer.accessibility import (
    AccessibilityProfile,
    TargetStructure,
    TargetStructurePredictor,
    classify_accessibility,
    export_profile,
    score_position_accessibility,
    score_window_accessibility,
)
from aso_designer.viennarna import RNAfoldResult, RNAplfoldResult

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def mock_profile() -> AccessibilityProfile:
    """AccessibilityProfile with known values for deterministic testing."""
    return AccessibilityProfile(
        sequence="AUGCAUGCAUGCAUGC",
        positions=list(range(1, 17)),
        unpaired_probs=[
            0.9,
            0.1,
            0.8,
            0.2,
            0.7,
            0.3,
            0.6,
            0.4,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
        ],
        window_size=150,
    )


@pytest.fixture
def mock_rnafold_result() -> RNAfoldResult:
    return RNAfoldResult(
        sequence="AUGCAUGCAUGCAUGC",
        structure="((((....))))....",
        mfe=-3.40,
    )


@pytest.fixture
def mock_rnaplfold_result() -> RNAplfoldResult:
    return RNAplfoldResult(
        positions=list(range(1, 17)),
        unpaired_probs=[
            0.9,
            0.1,
            0.8,
            0.2,
            0.7,
            0.3,
            0.6,
            0.4,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
            0.5,
        ],
    )


# ---------------------------------------------------------------------------
# TargetStructure dataclass
# ---------------------------------------------------------------------------


class TestTargetStructure:
    def test_from_rnafold_result(self, mock_rnafold_result):
        ts = TargetStructure.from_rnafold_result(mock_rnafold_result)
        assert ts.sequence == "AUGCAUGCAUGCAUGC"
        assert ts.structure == "((((....))))...."
        assert ts.mfe == -3.40


# ---------------------------------------------------------------------------
# AccessibilityProfile dataclass
# ---------------------------------------------------------------------------


class TestAccessibilityProfile:
    def test_from_rnaplfold_result(self, mock_rnaplfold_result):
        profile = AccessibilityProfile.from_rnaplfold_result(
            mock_rnaplfold_result, "AUGCAUGCAUGCAUGC", 150
        )
        assert profile.sequence == "AUGCAUGCAUGCAUGC"
        assert profile.window_size == 150
        assert len(profile.unpaired_probs) == 16

    def test_get_unpaired_prob_valid(self, mock_profile):
        assert mock_profile.get_unpaired_prob(1) == 0.9
        assert mock_profile.get_unpaired_prob(2) == 0.1

    def test_get_unpaired_prob_out_of_range(self, mock_profile):
        assert mock_profile.get_unpaired_prob(0) == 0.0
        assert mock_profile.get_unpaired_prob(100) == 0.0

    def test_from_rnaplfold_result_multi_length(self):
        result = RNAplfoldResult(
            positions=[1, 2, 3],
            unpaired_probs=[0.9, 0.8, 0.7],
            unpaired_lengths=[1, 2, 3],
            unpaired_probs_by_length={1: [0.9, 0.8, 0.7], 2: [0.6, 0.5, 0.4], 3: [0.3, 0.2, 0.1]},
        )
        profile = AccessibilityProfile.from_rnaplfold_result(result, "AUGCAUGC", 150)
        assert profile.unpaired_probs_by_length == {
            1: [0.9, 0.8, 0.7],
            2: [0.6, 0.5, 0.4],
            3: [0.3, 0.2, 0.1],
        }
        assert profile.get_unpaired_prob(2, length=1) == 0.8
        assert profile.get_unpaired_prob(2, length=2) == 0.5
        assert profile.get_unpaired_prob(2, length=3) == 0.2

    def test_get_unpaired_prob_unknown_length_raises(self, mock_profile):
        mock_profile.unpaired_probs_by_length = {1: list(mock_profile.unpaired_probs)}
        with pytest.raises(KeyError, match="not in profile"):
            mock_profile.get_unpaired_prob(1, length=5)

    def test_get_unpaired_prob_length_fallback(self, mock_profile):
        """Profiles without multi-length data fall back to the base column."""
        assert mock_profile.get_unpaired_prob(1, length=2) == 0.9


# ---------------------------------------------------------------------------
# Per-position accessibility scoring
# ---------------------------------------------------------------------------


class TestScorePositionAccessibility:
    def test_normal_values(self):
        assert score_position_accessibility(0.0) == 0.0
        assert score_position_accessibility(0.5) == 0.5
        assert score_position_accessibility(1.0) == 1.0

    def test_clamping_high(self):
        assert score_position_accessibility(1.5) == 1.0

    def test_clamping_low(self):
        assert score_position_accessibility(-0.5) == 0.0


# ---------------------------------------------------------------------------
# Windowed accessibility scoring
# ---------------------------------------------------------------------------


class TestScoreWindowAccessibility:
    def test_mean_aggregation(self, mock_profile):
        score = score_window_accessibility(mock_profile, 1, 4, "mean")
        expected = (0.9 + 0.1 + 0.8 + 0.2) / 4
        assert score == pytest.approx(expected)

    def test_min_aggregation(self, mock_profile):
        score = score_window_accessibility(mock_profile, 1, 4, "min")
        assert score == 0.1

    def test_max_aggregation(self, mock_profile):
        score = score_window_accessibility(mock_profile, 1, 4, "max")
        assert score == 0.9

    def test_median_odd(self, mock_profile):
        score = score_window_accessibility(mock_profile, 1, 5, "median")
        probs = [0.9, 0.1, 0.8, 0.2, 0.7]
        expected = sorted(probs)[2]
        assert score == expected

    def test_median_even(self, mock_profile):
        score = score_window_accessibility(mock_profile, 1, 4, "median")
        probs = [0.9, 0.1, 0.8, 0.2]
        sorted_probs = sorted(probs)
        expected = (sorted_probs[1] + sorted_probs[2]) / 2
        assert score == expected

    def test_single_position_window(self, mock_profile):
        score = score_window_accessibility(mock_profile, 5, 5, "mean")
        assert score == 0.7

    def test_full_window(self, mock_profile):
        score = score_window_accessibility(mock_profile, 1, 16, "mean")
        expected = sum(mock_profile.unpaired_probs) / 16
        assert score == pytest.approx(expected)

    def test_invalid_start(self, mock_profile):
        with pytest.raises(ValueError, match="Invalid window"):
            score_window_accessibility(mock_profile, 0, 4)

    def test_invalid_end(self, mock_profile):
        with pytest.raises(ValueError, match="Invalid window"):
            score_window_accessibility(mock_profile, 1, 20)

    def test_start_greater_than_end(self, mock_profile):
        with pytest.raises(ValueError, match="Invalid window"):
            score_window_accessibility(mock_profile, 4, 1)

    def test_unknown_aggregation(self, mock_profile):
        with pytest.raises(ValueError, match="Unknown aggregation"):
            score_window_accessibility(mock_profile, 1, 4, "sum")

    def test_default_aggregation_is_mean(self, mock_profile):
        score = score_window_accessibility(mock_profile, 1, 4)
        expected = (0.9 + 0.1 + 0.8 + 0.2) / 4
        assert score == pytest.approx(expected)


# ---------------------------------------------------------------------------
# Accessibility classification
# ---------------------------------------------------------------------------


class TestClassifyAccessibility:
    def test_accessible(self):
        assert classify_accessibility(0.6) == "accessible"

    def test_accessible_at_threshold(self):
        assert classify_accessibility(0.5) == "accessible"

    def test_moderate(self):
        assert classify_accessibility(0.4) == "moderate"

    def test_moderate_at_threshold(self):
        assert classify_accessibility(0.3) == "moderate"

    def test_inaccessible(self):
        assert classify_accessibility(0.2) == "inaccessible"

    def test_zero(self):
        assert classify_accessibility(0.0) == "inaccessible"

    def test_one(self):
        assert classify_accessibility(1.0) == "accessible"

    def test_custom_thresholds(self):
        assert classify_accessibility(0.4, accessible_threshold=0.7) == "moderate"
        assert (
            classify_accessibility(0.2, accessible_threshold=0.7, moderate_threshold=0.1)
            == "moderate"
        )
        assert (
            classify_accessibility(0.05, accessible_threshold=0.7, moderate_threshold=0.1)
            == "inaccessible"
        )


# ---------------------------------------------------------------------------
# Profile export
# ---------------------------------------------------------------------------


class TestExportProfile:
    def test_export_json(self, mock_profile, tmp_path):
        out = tmp_path / "profile.json"
        result = export_profile(mock_profile, out, "json")
        assert result == out
        assert out.exists()
        with open(out) as f:
            data = json.load(f)
        assert data["sequence"] == mock_profile.sequence
        assert data["window_size"] == 150
        assert data["positions"] == list(range(1, 17))
        assert len(data["accessibility_scores"]) == 16
        assert data["accessibility_scores"][0] == 0.9

    def test_export_csv(self, mock_profile, tmp_path):
        out = tmp_path / "profile.csv"
        result = export_profile(mock_profile, out, "csv")
        assert result == out
        assert out.exists()
        content = out.read_text()
        assert "position,unpaired_prob,accessibility_score" in content
        lines = content.strip().split("\n")
        assert len(lines) == 17  # header + 16 data rows

    def test_export_creates_parent_dirs(self, mock_profile, tmp_path):
        out = tmp_path / "subdir" / "deeper" / "profile.json"
        export_profile(mock_profile, out, "json")
        assert out.exists()

    def test_export_unsupported_format(self, mock_profile, tmp_path):
        with pytest.raises(ValueError, match="Unsupported format"):
            export_profile(mock_profile, tmp_path / "out.txt", "txt")


# ---------------------------------------------------------------------------
# TargetStructurePredictor (mocked ViennaRNA)
# ---------------------------------------------------------------------------


class TestTargetStructurePredictor:
    @patch("aso_designer.accessibility.run_rnafold")
    def test_predict_structure(self, mock_rnafold, mock_rnafold_result):
        mock_rnafold.return_value = mock_rnafold_result
        predictor = TargetStructurePredictor()
        ts = predictor.predict_structure("ATGCATGCATGCATGC")
        assert ts.structure == "((((....))))...."
        assert ts.mfe == -3.40
        mock_rnafold.assert_called_once()

    @patch("aso_designer.accessibility.run_rnaplfold")
    def test_predict_accessibility(self, mock_rnaplfold, mock_rnaplfold_result):
        mock_rnaplfold.return_value = mock_rnaplfold_result
        predictor = TargetStructurePredictor(window_size=80, max_pairing_distance=50)
        profile = predictor.predict_accessibility("ATGCATGCATGCATGC")
        assert len(profile.unpaired_probs) == 16
        assert profile.window_size == 80
        mock_rnaplfold.assert_called_once()

    @patch("aso_designer.accessibility.run_rnafold")
    @patch("aso_designer.accessibility.run_rnaplfold")
    def test_predict_all(
        self,
        mock_rnaplfold,
        mock_rnafold,
        mock_rnafold_result,
        mock_rnaplfold_result,
    ):
        mock_rnafold.return_value = mock_rnafold_result
        mock_rnaplfold.return_value = mock_rnaplfold_result
        predictor = TargetStructurePredictor()
        ts, profile = predictor.predict_all("ATGCATGCATGCATGC")
        assert ts.mfe == -3.40
        assert len(profile.unpaired_probs) == 16


# ---------------------------------------------------------------------------
# TargetStructurePredictor caching
# ---------------------------------------------------------------------------


class TestStructureCaching:
    @patch("aso_designer.accessibility.run_rnafold")
    def test_structure_cached(self, mock_rnafold, mock_rnafold_result, tmp_path):
        mock_rnafold.return_value = mock_rnafold_result
        predictor = TargetStructurePredictor(cache_dir=tmp_path)

        # First call: computes and caches
        ts1 = predictor.predict_structure("ATGCATGCATGCATGC")
        assert mock_rnafold.call_count == 1

        # Second call: loads from cache, does not call RNAfold
        ts2 = predictor.predict_structure("ATGCATGCATGCATGC")
        assert mock_rnafold.call_count == 1
        assert ts2.mfe == ts1.mfe
        assert ts2.structure == ts1.structure

    @patch("aso_designer.accessibility.run_rnaplfold")
    def test_accessibility_cached(self, mock_rnaplfold, mock_rnaplfold_result, tmp_path):
        mock_rnaplfold.return_value = mock_rnaplfold_result
        predictor = TargetStructurePredictor(cache_dir=tmp_path)

        # First call
        p1 = predictor.predict_accessibility("ATGCATGCATGCATGC")
        assert mock_rnaplfold.call_count == 1

        # Second call: cached
        p2 = predictor.predict_accessibility("ATGCATGCATGCATGC")
        assert mock_rnaplfold.call_count == 1
        assert p2.unpaired_probs == p1.unpaired_probs

    @patch("aso_designer.accessibility.run_rnafold")
    def test_different_sequences_not_cached(self, mock_rnafold, mock_rnafold_result, tmp_path):
        mock_rnafold.return_value = mock_rnafold_result
        predictor = TargetStructurePredictor(cache_dir=tmp_path)

        predictor.predict_structure("ATGCATGCATGCATGC")
        predictor.predict_structure("TTTTAAAATTTTAAAA")
        assert mock_rnafold.call_count == 2

    @patch("aso_designer.accessibility.run_rnafold")
    def test_no_caching_without_dir(self, mock_rnafold, mock_rnafold_result):
        mock_rnafold.return_value = mock_rnafold_result
        predictor = TargetStructurePredictor(cache_dir=None)

        predictor.predict_structure("ATGCATGCATGCATGC")
        predictor.predict_structure("ATGCATGCATGCATGC")
        assert mock_rnafold.call_count == 2

    @patch("aso_designer.accessibility.run_rnaplfold")
    def test_different_params_different_cache(
        self,
        mock_rnaplfold,
        mock_rnaplfold_result,
        tmp_path,
    ):
        mock_rnaplfold.return_value = mock_rnaplfold_result
        predictor_w80 = TargetStructurePredictor(cache_dir=tmp_path, window_size=80)
        predictor_w150 = TargetStructurePredictor(cache_dir=tmp_path, window_size=150)

        # Different window sizes should use different cache keys
        predictor_w80.predict_accessibility("ATGCATGC")
        predictor_w150.predict_accessibility("ATGCATGC")
        assert mock_rnaplfold.call_count == 2


# ---------------------------------------------------------------------------
# RNAup accessibility method (mocked ViennaRNA)
# ---------------------------------------------------------------------------


@pytest.fixture
def mock_rnaup_result():
    from aso_designer.viennarna import RNAupResult

    return RNAupResult(
        sequence="AUGCAUGCAUGCAUGC",
        positions=list(range(1, 17)),
        unpaired_probs=[0.9, 0.1] * 8,
        region_probs=[0.5, 0.4] * 8,
        opening_dg=[0.11, 0.92] * 8,
        region_length=4,
    )


class TestAccessibilityProfileFromRnaup:
    def test_from_rnaup_result(self, mock_rnaup_result):
        profile = AccessibilityProfile.from_rnaup_result(mock_rnaup_result)
        assert profile.method == "rnaup"
        assert profile.window_size == 0
        assert profile.region_length == 4
        assert profile.unpaired_probs == mock_rnaup_result.unpaired_probs
        assert profile.region_probs == mock_rnaup_result.region_probs
        assert profile.opening_dg == mock_rnaup_result.opening_dg

    def test_get_opening_dg(self, mock_rnaup_result):
        profile = AccessibilityProfile.from_rnaup_result(mock_rnaup_result)
        assert profile.get_opening_dg(1) == 0.11
        assert profile.get_opening_dg(2) == 0.92
        assert profile.get_opening_dg(17) is None
        assert profile.get_opening_dg(0) is None

    def test_get_opening_dg_none_for_rnaplfold(self, mock_profile):
        assert mock_profile.get_opening_dg(1) is None


class TestPredictorRnaupMethod:
    @patch("aso_designer.accessibility.run_rnaup")
    def test_predict_accessibility_rnaup(self, mock_rnaup, mock_rnaup_result, tmp_path):
        mock_rnaup.return_value = mock_rnaup_result
        predictor = TargetStructurePredictor(cache_dir=tmp_path, method="rnaup", unpaired_region=4)
        profile = predictor.predict_accessibility("ATGCATGCATGCATGC")
        assert profile.method == "rnaup"
        assert profile.region_length == 4
        assert profile.opening_dg == mock_rnaup_result.opening_dg
        mock_rnaup.assert_called_once()
        # unpaired_region is passed through to RNAup
        assert mock_rnaup.call_args.kwargs["unpaired_region"] == 4

    @patch("aso_designer.accessibility.run_rnaplfold")
    def test_invalid_method_rejected(self, mock_rnaplfold):
        with pytest.raises(ValueError, match="Unknown accessibility method"):
            TargetStructurePredictor(method="mfold")

    @patch("aso_designer.accessibility.run_rnaup")
    def test_rnaup_result_cached(self, mock_rnaup, mock_rnaup_result, tmp_path):
        mock_rnaup.return_value = mock_rnaup_result
        predictor = TargetStructurePredictor(cache_dir=tmp_path, method="rnaup", unpaired_region=4)
        predictor.predict_accessibility("ATGCATGCATGCATGC")
        predictor.predict_accessibility("ATGCATGCATGCATGC")
        assert mock_rnaup.call_count == 1

    @patch("aso_designer.accessibility.run_rnaplfold")
    @patch("aso_designer.accessibility.run_rnaup")
    def test_methods_use_distinct_cache_keys(
        self, mock_rnaup, mock_rnaplfold, mock_rnaplfold_result, mock_rnaup_result, tmp_path
    ):
        mock_rnaplfold.return_value = mock_rnaplfold_result
        mock_rnaup.return_value = mock_rnaup_result
        p1 = TargetStructurePredictor(cache_dir=tmp_path, method="rnaplfold")
        p2 = TargetStructurePredictor(cache_dir=tmp_path, method="rnaup")
        p1.predict_accessibility("ATGCATGC")
        p2.predict_accessibility("ATGCATGC")
        assert mock_rnaplfold.call_count == 1
        assert mock_rnaup.call_count == 1
        # same-method different-region predictions also miss the cache
        p3 = TargetStructurePredictor(cache_dir=tmp_path, method="rnaup", unpaired_region=6)
        p3.predict_accessibility("ATGCATGC")
        assert mock_rnaup.call_count == 2


# ---------------------------------------------------------------------------
# TargetStructurePredictor multi-length RNAplfold
# ---------------------------------------------------------------------------


class TestPredictorMultiLength:
    @patch("aso_designer.accessibility.run_rnaplfold")
    def test_lengths_passed_to_rnaplfold(self, mock_rnaplfold, mock_rnaplfold_result):
        mock_rnaplfold.return_value = mock_rnaplfold_result
        predictor = TargetStructurePredictor(unpaired_lengths=[1, 2, 3, 4])
        profile = predictor.predict_accessibility("ATGCATGCATGCATGC")
        assert mock_rnaplfold.call_args.kwargs["unpaired_region"] == [1, 2, 3, 4]
        assert profile.unpaired_probs_by_length == {1: list(mock_rnaplfold_result.unpaired_probs)}

    @patch("aso_designer.accessibility.run_rnaplfold")
    def test_length_sets_use_distinct_cache(self, mock_rnaplfold, mock_rnaplfold_result, tmp_path):
        mock_rnaplfold.return_value = mock_rnaplfold_result
        p1 = TargetStructurePredictor(cache_dir=tmp_path, unpaired_lengths=[1])
        p2 = TargetStructurePredictor(cache_dir=tmp_path, unpaired_lengths=[1, 2, 3, 4])
        p1.predict_accessibility("ATGCATGC")
        p2.predict_accessibility("ATGCATGC")
        assert mock_rnaplfold.call_count == 2


# ---------------------------------------------------------------------------
# Accessibility config plumbing
# ---------------------------------------------------------------------------


class TestAccessibilityMethodConfig:
    def test_schema_defaults(self):
        from aso_designer.config_schema import ASOConfig

        cfg = ASOConfig()
        assert cfg.accessibility.method == "rnaplfold"
        assert cfg.accessibility.unpaired_region == 4
        assert cfg.accessibility.window_size == 150

    def test_method_accepted(self):
        from aso_designer.config_schema import ASOConfig

        assert ASOConfig(accessibility={"method": "rnaup"}).accessibility.method == "rnaup"

    def test_method_rejected(self):
        from aso_designer.config_schema import ASOConfig

        with pytest.raises(ValidationError):
            ASOConfig(accessibility={"method": "mfold"})

    def test_unpaired_region_bounds(self):
        from aso_designer.config_schema import ASOConfig

        assert ASOConfig(accessibility={"unpaired_region": 1}).accessibility.unpaired_region == 1
        assert ASOConfig(accessibility={"unpaired_region": 15}).accessibility.unpaired_region == 15
        with pytest.raises(ValidationError):
            ASOConfig(accessibility={"unpaired_region": 0})
        with pytest.raises(ValidationError):
            ASOConfig(accessibility={"unpaired_region": 16})

    def test_unpaired_lengths_default(self):
        from aso_designer.config_schema import ASOConfig

        assert ASOConfig().accessibility.unpaired_lengths == [1]

    def test_unpaired_lengths_accepted_and_sorted(self):
        from aso_designer.config_schema import ASOConfig

        assert ASOConfig(
            accessibility={"unpaired_lengths": [4, 2, 1]}
        ).accessibility.unpaired_lengths == [1, 2, 4]

    def test_unpaired_lengths_rejected(self):
        from aso_designer.config_schema import ASOConfig

        with pytest.raises(ValidationError):
            ASOConfig(accessibility={"unpaired_lengths": []})
        with pytest.raises(ValidationError):
            ASOConfig(accessibility={"unpaired_lengths": [2, 2]})
        with pytest.raises(ValidationError):
            ASOConfig(accessibility={"unpaired_lengths": [0, 1]})
        with pytest.raises(ValidationError):
            ASOConfig(accessibility={"unpaired_lengths": [1, 16]})

    def test_yaml_default_loads(self):
        """The shipped YAML parses and carries the new keys."""
        from aso_designer.config_loader import load_config

        cfg = load_config()
        assert cfg.accessibility.method == "rnaplfold"
        assert cfg.accessibility.unpaired_region == 4
        assert cfg.accessibility.unpaired_lengths == [1]
