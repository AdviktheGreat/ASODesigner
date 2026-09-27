"""Configuration loader for ASODesigner.

Loads YAML configuration files, merges them with defaults, validates
against the Pydantic schema, and optionally applies chemistry-specific
presets. Designed from scratch with layered override support
(user config → chemistry preset → defaults).
"""

from pathlib import Path

import yaml

from aso_designer.config_schema import ASOConfig
from aso_designer.utils.exceptions import ConfigurationError
from aso_designer.utils.logger import get_logger

logger = get_logger(__name__)

# Path to the bundled default config and chemistry presets.
_DEFAULT_CONFIG_PATH = Path(__file__).parent / "config" / "aso_config.yaml"
_CHEMISTRY_PRESETS_PATH = Path(__file__).parent / "config" / "chemistry_presets.yaml"


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

    # Apply user overrides
    if user_config_path is not None:
        user_data = load_yaml(user_config_path)
        config_data = _deep_merge(config_data, user_data)
        logger.debug("Merged user config from %s", user_config_path)

        # Re-apply chemistry preset if user changed the chemistry
        user_chemistry = user_data.get("candidate", {}).get("chemistry")
        if user_chemistry and user_chemistry != effective_chemistry:
            config_data = _apply_chemistry_preset(config_data, user_chemistry)

    # Validate
    try:
        config = ASOConfig(**config_data)
    except Exception as e:
        raise ConfigurationError(f"Config validation failed: {e}") from e

    logger.info("Configuration loaded and validated successfully")
    return config


def _apply_chemistry_preset(config_data: dict, chemistry: str) -> dict:
    """Apply chemistry-specific preset adjustments to config data."""
    try:
        presets = load_yaml(_CHEMISTRY_PRESETS_PATH)
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

    return config_data


def get_default_config() -> ASOConfig:
    """Return the default configuration without any overrides."""
    return load_config()
