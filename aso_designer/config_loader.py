"""Configuration loader for ASODesigner.

Loads YAML configuration files, merges them with defaults, validates
against the Pydantic schema, and optionally applies chemistry-specific
presets. Designed from scratch with layered override support
(user config → chemistry preset → defaults).
"""

from pathlib import Path

import yaml
from pydantic import ValidationError

from aso_designer.config_schema import ASOConfig
from aso_designer.utils.exceptions import ConfigurationError
from aso_designer.utils.logger import get_logger

logger = get_logger(__name__)

# Path to the bundled default config and chemistry presets.
_DEFAULT_CONFIG_PATH = Path(__file__).parent / "config" / "aso_config.yaml"
CHEMISTRY_PRESETS_PATH = Path(__file__).parent / "config" / "chemistry_presets.yaml"


def load_yaml(path: str | Path) -> dict:
    """Load a YAML file into a dictionary.

    Parameters
    ----------
    path : str or Path
        Path to the YAML file.

    Returns
    -------
    dict
        Parsed YAML content.

    Raises
    ------
    ConfigurationError
        If the file cannot be read or parsed.
    """
    path = Path(path)
    if not path.exists():
        raise ConfigurationError(f"Config file not found: {path}")
    try:
        with open(path, encoding="utf-8") as f:
            data = yaml.safe_load(f)
    except yaml.YAMLError as e:
        raise ConfigurationError(f"YAML parse error in {path}: {e}") from e
    return data or {}


def _deep_merge(base: dict, override: dict) -> dict:
    """Recursively merge override into base (override wins)."""
    result = base.copy()
    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = value
    return result


def load_config(
    user_config_path: str | Path | None = None,
    chemistry: str | None = None,
    profile: str | None = None,
) -> ASOConfig:
    """Load and validate ASODesigner configuration.

    Merge order (later overrides earlier):
      1. Bundled default config
      2. Chemistry preset (if specified)
      3. User-provided config file (if specified)

    Parameters
    ----------
    user_config_path : str or Path, optional
        Path to a user YAML config that overrides defaults.
    chemistry : str, optional
        Chemistry name to apply preset adjustments (e.g. "gapmer").
        Overrides the chemistry in the config if provided.
    profile : str, optional
        Named scoring weight profile (balanced, efficacy, specificity,
        safety). Overrides the config's profile selection when provided.

    Returns
    -------
    ASOConfig
        Validated configuration object.

    Raises
    ------
    ConfigurationError
        If loading or validation fails.
    """
    # Start with defaults
    config_data = load_yaml(_DEFAULT_CONFIG_PATH)
    logger.debug("Loaded default config from %s", _DEFAULT_CONFIG_PATH)

    # Apply chemistry preset adjustments
    effective_chemistry = chemistry or config_data.get("candidate", {}).get("chemistry")
    if effective_chemistry:
        config_data = _apply_chemistry_preset(config_data, effective_chemistry)
        # The preset only adjusts parameters; make sure the selected
        # chemistry itself is what the validated config carries (user
        # overrides below win over this value via the deep merge).
        config_data.setdefault("candidate", {})["chemistry"] = effective_chemistry

    # Apply user overrides
    if user_config_path is not None:
        user_data = load_yaml(user_config_path)
        config_data = _deep_merge(config_data, user_data)
        logger.debug("Merged user config from %s", user_config_path)

        # Re-apply chemistry preset if user changed the chemistry
        user_chemistry = user_data.get("candidate", {}).get("chemistry")
        if user_chemistry and user_chemistry != effective_chemistry:
            config_data = _apply_chemistry_preset(config_data, user_chemistry)

    # Resolve the scoring weight profile. An explicit profile argument
    # (e.g. from the CLI) wins; otherwise the config's own profile is
    # used, except when the user supplied explicit weights (those win).
    scoring_data = config_data.setdefault("scoring", {})
    effective_profile = profile or scoring_data.get("profile", "balanced")
    user_set_weights = "weights" in scoring_data and profile is None
    if effective_profile and effective_profile != "balanced" and not user_set_weights:
        from aso_designer.config_schema import WEIGHT_PROFILES

        weights = WEIGHT_PROFILES[effective_profile]
        scoring_data["weights"] = dict(weights)
        scoring_data["profile"] = effective_profile
        logger.info("Applied scoring weight profile '%s'", effective_profile)

    # Validate
    try:
        config = ASOConfig(**config_data)
    except ValidationError as e:
        raise ConfigurationError(f"Config validation failed: {e}") from e

    logger.info("Configuration loaded and validated successfully")
    return config


def load_chemistry_presets() -> dict:
    """Load the bundled chemistry presets YAML.

    Returns
    -------
    dict
        Chemistry presets keyed by chemistry name; empty if the file
        is missing or unreadable.

    Raises
    ------
    ConfigurationError
        If the file cannot be read or parsed (callers may catch this
        to degrade gracefully).
    """
    return load_yaml(CHEMISTRY_PRESETS_PATH)


def _apply_chemistry_preset(config_data: dict, chemistry: str) -> dict:
    """Apply chemistry-specific preset adjustments to config data."""
    try:
        presets = load_chemistry_presets()
    except ConfigurationError:
        logger.warning("Chemistry presets file not found; skipping preset.")
        return config_data

    preset = presets.get(chemistry)
    if preset is None:
        logger.warning("No preset for chemistry '%s'; skipping.", chemistry)
        return config_data

    # Apply Tm adjustment to thermodynamics config
    tm_adj = preset.get("tm_adjustment", 0.0)
    if tm_adj:
        thermo = config_data.setdefault("thermodynamics", {})
        thermo["tm_min"] = thermo.get("tm_min", 40.0) + tm_adj
        thermo["tm_max"] = thermo.get("tm_max", 65.0) + tm_adj
        thermo["tm_optimal"] = thermo.get("tm_optimal", 52.0) + tm_adj
        logger.debug("Applied Tm adjustment %+f for %s", tm_adj, chemistry)

    # Store chemistry bonus as metadata for the scoring module to use
    # (do NOT modify scoring weights — they must always sum to 1.0)
    bonus = preset.get("chemistry_bonus", 0.0)
    if bonus:
        config_data.setdefault("scoring", {})["chemistry_bonus"] = bonus
        logger.debug("Stored chemistry bonus %f for %s", bonus, chemistry)

    # Scale the self-structure weight by the preset factor, then
    # renormalize the weight vector so it still sums to 1.0. Later
    # overrides (user weights block or a weight profile) win wholesale.
    adjust = preset.get("self_structure_weight_adjust")
    if adjust and adjust != 1.0:
        scoring = config_data.setdefault("scoring", {})
        weights = dict(
            scoring.get("weights")
            or {
                "accessibility": 0.25,
                "thermodynamics": 0.20,
                "gc_content": 0.10,
                "off_target": 0.20,
                "variant": 0.10,
                "self_structure": 0.10,
                "positional": 0.05,
            }
        )
        weights["self_structure"] = weights.get("self_structure", 0.10) * adjust
        total = sum(weights.values())
        scoring["weights"] = {k: round(v / total, 6) for k, v in weights.items()}
        logger.debug(
            "Applied self_structure weight adjust %s for %s (renormalized)",
            adjust,
            chemistry,
        )

    # Surface the preset's RNase H eligibility flag (e.g. morpholino:
    # false) as a candidate-level override consumed by the pipeline.
    rnase_h = preset.get("rnase_h_eligible")
    if rnase_h is not None:
        config_data.setdefault("candidate", {})["rnase_h_eligible"] = bool(rnase_h)
        logger.debug("Stored rnase_h_eligible=%s for %s", rnase_h, chemistry)

    return config_data


def get_default_config() -> ASOConfig:
    """Return the default configuration without any overrides."""
    return load_config()
