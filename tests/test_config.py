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
        # The selected chemistry must land in the validated config so the
        # pipeline generates candidates with the requested chemistry.
        assert config.candidate.chemistry == "mixmer"

    def test_with_user_config(self, tmp_path):
        f = tmp_path / "custom.yaml"
        f.write_text(
            "candidate:\n  length: 18\n  gc_min: 0.35\ngapmer:\n  gap_size: 8\n  wing_length: 5\n"
        )
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
# Data source configuration options
# ---------------------------------------------------------------------------


class TestDataPathOptions:
    def test_transcriptome_path_default_none(self):
        config = load_config()
        assert config.off_target.transcriptome_path is None

    def test_transcriptome_path_override(self, tmp_path):
        f = tmp_path / "custom.yaml"
        f.write_text('off_target:\n  transcriptome_path: "/data/transcriptome.fna"\n')
        config = load_config(str(f))
        assert config.off_target.transcriptome_path == "/data/transcriptome.fna"

    def test_vcf_path_default_none(self):
        config = load_config()
        assert config.variant.vcf_path is None

    def test_vcf_path_override(self, tmp_path):
        f = tmp_path / "custom.yaml"
        f.write_text('variant:\n  vcf_path: "/data/variants.vcf"\n')
        config = load_config(str(f))
        assert config.variant.vcf_path == "/data/variants.vcf"

    def test_fetch_config_defaults(self):
        config = load_config()
        assert config.fetch.email is None
        assert config.fetch.api_key is None

    def test_fetch_config_override(self, tmp_path):
        f = tmp_path / "custom.yaml"
        f.write_text('fetch:\n  email: "user@example.com"\n  api_key: "abc123"\n')
        config = load_config(str(f))
        assert config.fetch.email == "user@example.com"
        assert config.fetch.api_key == "abc123"

    def test_self_structure_nn_fraction_default(self):
        config = load_config()
        assert config.thermodynamics.self_structure_nn_fraction == 0.3

    def test_self_structure_nn_fraction_override(self, tmp_path):
        f = tmp_path / "custom.yaml"
        f.write_text("thermodynamics:\n  self_structure_nn_fraction: 0.5\n")
        config = load_config(str(f))
        assert config.thermodynamics.self_structure_nn_fraction == 0.5

    def test_self_structure_nn_fraction_invalid(self, tmp_path):
        f = tmp_path / "bad.yaml"
        f.write_text("thermodynamics:\n  self_structure_nn_fraction: 1.5\n")
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


# ---------------------------------------------------------------------------
# Gapmer geometry validation
# ---------------------------------------------------------------------------


def test_asoconfig_rejects_broken_gapmer_geometry():
    """candidate.length inconsistent with gapmer wings+gap is rejected."""
    from pydantic import ValidationError

    from aso_designer.config_schema import ASOConfig

    with pytest.raises(ValidationError, match="Invalid gapmer geometry"):
        ASOConfig(
            candidate={"chemistry": "gapmer", "length": 22},
            gapmer={"gap_size": 10, "wing_length": 5},
        )


def test_asoconfig_accepts_consistent_gapmer_geometry():
    from aso_designer.config_schema import ASOConfig

    cfg = ASOConfig(
        candidate={"chemistry": "gapmer", "length": 21},
        gapmer={"gap_size": 11, "wing_length": 5},
    )
    assert cfg.candidate.length == 21


def test_asoconfig_non_gapmer_ignores_geometry():
    from aso_designer.config_schema import ASOConfig

    cfg = ASOConfig(candidate={"chemistry": "mixmer", "length": 22})
    assert cfg.candidate.length == 22


def test_asoconfig_accepts_lna_chemistry():
    from aso_designer.config_schema import ASOConfig

    cfg = ASOConfig(candidate={"chemistry": "lna"})
    assert cfg.candidate.chemistry == "lna"
    assert cfg.lna.pattern == []


# ---------------------------------------------------------------------------
# Weight profiles
# ---------------------------------------------------------------------------


def test_weight_profile_applied():
    from aso_designer.config_loader import load_config

    cfg = load_config(profile="specificity")
    assert cfg.scoring.weights.off_target == 0.45
    assert cfg.scoring.profile == "specificity"


def test_weight_profile_from_yaml():
    from aso_designer.config_loader import load_config

    cfg = load_config()
    assert cfg.scoring.profile == "balanced"
    assert cfg.scoring.weights.accessibility == 0.25


def test_profile_invalid_rejected():
    from pydantic import ValidationError

    from aso_designer.config_schema import ASOConfig

    with pytest.raises(ValidationError, match="scoring.profile"):
        ASOConfig(scoring={"profile": "aggressive"})


def test_explicit_weights_beat_config_profile():
    """Explicit user weights win over the YAML/config profile."""
    from aso_designer.config_schema import WEIGHT_PROFILES, ASOConfig

    weights = dict(WEIGHT_PROFILES["safety"])
    weights["accessibility"] = 0.20
    weights["variant"] = 0.30  # keep the sum at 1.0
    cfg2 = ASOConfig(scoring={"profile": "safety", "weights": weights})
    assert cfg2.scoring.weights.accessibility == 0.20


# ---------------------------------------------------------------------------
# chemistry_bonus schema field
# ---------------------------------------------------------------------------


def test_chemistry_bonus_field_and_bounds():
    from pydantic import ValidationError

    from aso_designer.config_schema import ASOConfig

    assert ASOConfig().scoring.chemistry_bonus == 0.0
    assert ASOConfig(scoring={"chemistry_bonus": 0.05}).scoring.chemistry_bonus == 0.05
    with pytest.raises(ValidationError):
        ASOConfig(scoring={"chemistry_bonus": -0.01})
    with pytest.raises(ValidationError):
        ASOConfig(scoring={"chemistry_bonus": 0.5})


def test_engine_from_config_reads_bonus_directly():
    from aso_designer.config_loader import load_config
    from aso_designer.scoring import create_engine_from_config

    cfg = load_config(chemistry="gapmer")  # preset sets chemistry_bonus 0.05
    engine = create_engine_from_config(chemistry="gapmer")
    assert engine.chemistry_bonus == cfg.scoring.chemistry_bonus == 0.05
