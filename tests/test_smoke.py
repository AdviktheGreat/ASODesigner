"""Smoke tests verifying that core modules import and basic functionality works."""

from aso_designer.config_loader import load_config
from aso_designer.config_schema import ASOConfig
from aso_designer.utils.sequence import gc_content, reverse_complement, validate_sequence


def test_package_imports() -> None:
    """Verify the package can be imported."""
    import aso_designer

    assert aso_designer.__version__ == "1.0.0"


def test_public_api_exposed() -> None:
    """Every name in __all__ resolves to a real attribute."""
    import aso_designer

    for name in aso_designer.__all__:
        assert hasattr(aso_designer, name), f"missing public API name: {name}"


def test_config_loads() -> None:
    """Verify the default configuration loads and validates."""
    config = load_config()
    assert isinstance(config, ASOConfig)
    assert config.candidate.chemistry == "gapmer"
    assert config.candidate.length == 20


def test_scoring_weights_sum_to_one() -> None:
    """Verify scoring weights sum to 1.0."""
    config = load_config()
    weights = config.scoring.weights
    total = (
        weights.accessibility
        + weights.thermodynamics
        + weights.gc_content
        + weights.off_target
        + weights.variant
        + weights.self_structure
        + weights.positional
    )
    assert abs(total - 1.0) < 0.001


def test_sequence_utils() -> None:
    """Verify basic sequence utilities work."""
    seq = "ATGCATGC"
    assert validate_sequence(seq) == "ATGCATGC"
    assert gc_content(seq) == 0.5
    assert reverse_complement(seq) == "GCATGCAT"


def test_exceptions_importable() -> None:
    """Verify exception classes are importable."""
    from aso_designer.utils.exceptions import ASODesignError

    assert issubclass(ASODesignError, Exception)
