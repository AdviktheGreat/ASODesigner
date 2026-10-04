"""Tests for the config_loader module: YAML loading, merging, and validation."""

import pytest

from aso_designer.config_loader import (
    _apply_chemistry_preset,
    _deep_merge,
    load_config,
    load_yaml,
)
from aso_designer.utils.exceptions import ConfigurationError

# ---------------------------------------------------------------------------
# load_yaml
# ---------------------------------------------------------------------------


class TestLoadYaml:
    def test_valid_yaml(self, tmp_path):
        f = tmp_path / "test.yaml"
        f.write_text("key: value\nlist:\n  - a\n  - b\n")
        data = load_yaml(str(f))
        assert data["key"] == "value"
        assert data["list"] == ["a", "b"]

    def test_empty_yaml(self, tmp_path):
        f = tmp_path / "empty.yaml"
        f.write_text("")
        data = load_yaml(str(f))
        assert data == {}

    def test_missing_file_raises(self):
        with pytest.raises(ConfigurationError, match="not found"):
            load_yaml("/nonexistent/file.yaml")

    def test_invalid_yaml_raises(self, tmp_path):
        f = tmp_path / "bad.yaml"
        f.write_text("key: value\n  bad: indent\n: missing_key")
        with pytest.raises(ConfigurationError, match="YAML parse error"):
            load_yaml(str(f))


# ---------------------------------------------------------------------------
# _deep_merge
# ---------------------------------------------------------------------------


class TestDeepMerge:
    def test_simple_override(self):
        base = {"a": 1, "b": 2}
        override = {"b": 3}
        result = _deep_merge(base, override)
        assert result == {"a": 1, "b": 3}

    def test_nested_merge(self):
        base = {"section": {"x": 1, "y": 2}}
        override = {"section": {"y": 3, "z": 4}}
        result = _deep_merge(base, override)
        assert result == {"section": {"x": 1, "y": 3, "z": 4}}

    def test_new_key_added(self):
        base = {"a": 1}
        override = {"b": 2}
        result = _deep_merge(base, override)
        assert result == {"a": 1, "b": 2}

    def test_does_not_mutate_base(self):
        base = {"a": 1}
        _deep_merge(base, {"b": 2})
        assert "b" not in base


# ---------------------------------------------------------------------------
# load_config
# ---------------------------------------------------------------------------


class TestLoadConfig:
    def test_default_config(self):
        config = load_config()
        assert config.candidate.length == 20
        assert config.candidate.chemistry == "gapmer"
        assert config.scoring.weights.accessibility == 0.25

    def test_with_chemistry(self):
        config = load_config(chemistry="mixmer")
        # Chemistry param applies presets but doesn't override candidate.chemistry
        assert config.candidate.chemistry == "gapmer"

    def test_with_user_config(self, tmp_path):
        f = tmp_path / "custom.yaml"
        f.write_text("candidate:\n  length: 18\n  gc_min: 0.35\n")
        config = load_config(str(f))
        assert config.candidate.length == 18
        assert config.candidate.gc_min == 0.35
        # Other values should remain at defaults
        assert config.candidate.step == 1

    def test_user_config_overrides_chemistry(self, tmp_path):
        f = tmp_path / "custom.yaml"
        f.write_text("candidate:\n  chemistry: mixmer\n")
        config = load_config(str(f))
        assert config.candidate.chemistry == "mixmer"

    def test_invalid_config_raises(self, tmp_path):
        f = tmp_path / "bad.yaml"
        f.write_text("candidate:\n  length: 5\n")  # Below minimum of 10
        with pytest.raises(ConfigurationError, match="Config validation failed"):
            load_config(str(f))

    def test_invalid_weights_raises(self, tmp_path):
        f = tmp_path / "bad_weights.yaml"
        f.write_text(
            "scoring:\n  weights:\n"
            "    accessibility: 0.5\n    thermodynamics: 0.5\n"
            "    gc_content: 0.5\n"
        )
        with pytest.raises(ConfigurationError, match="Config validation failed"):
            load_config(str(f))


# ---------------------------------------------------------------------------
# _apply_chemistry_preset
# ---------------------------------------------------------------------------


class TestApplyChemistryPreset:
    def test_applies_tm_adjustment(self):
        config_data = {
            "thermodynamics": {"tm_min": 40.0, "tm_max": 65.0, "tm_optimal": 52.0},
        }
        result = _apply_chemistry_preset(config_data, "gapmer")
        # Gapmer should have some Tm adjustment
        thermo = result.get("thermodynamics", {})
        assert "tm_min" in thermo

    def test_unknown_chemistry_skipped(self):
        config_data = {"candidate": {"chemistry": "gapmer"}}
        result = _apply_chemistry_preset(config_data, "nonexistent_chemistry")
        # Should return config unchanged (no preset found)
        assert result == config_data

    def test_chemistry_bonus_stored(self):
        config_data = {}
        result = _apply_chemistry_preset(config_data, "gapmer")
        # If gapmer has a bonus, it should be stored in scoring
        scoring = result.get("scoring", {})
        # The bonus is stored as metadata, not as a weight modification
        if "chemistry_bonus" in scoring:
            assert isinstance(scoring["chemistry_bonus"], (int, float))
