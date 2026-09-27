"""Exception hierarchy for ASODesigner.

All ASODesigner-specific exceptions derive from ASODesignError so callers
can catch the full family with a single except clause.
"""


class ASODesignError(Exception):
    """Base exception for all ASODesigner errors."""


class TranscriptFetchError(ASODesignError):
    """Raised when transcript sequence fetching fails (network, cache, or file)."""


class InvalidSequenceError(ASODesignError):
    """Raised when an input sequence fails validation (invalid characters, length)."""


class ConfigurationError(ASODesignError):
    """Raised when configuration loading or validation fails."""


class ViennaRNAError(ASODesignError):
    """Raised when ViennaRNA tools are missing or produce unexpected output."""


class ThermodynamicsError(ASODesignError):
    """Raised when thermodynamic feature computation fails."""


class OffTargetError(ASODesignError):
    """Raised when off-target analysis fails (index, search, or output)."""


class VariantError(ASODesignError):
    """Raised when variant/SNP data retrieval or processing fails."""


class ScoringError(ASODesignError):
    """Raised when candidate scoring or ranking fails."""


class ReportError(ASODesignError):
    """Raised when report generation fails."""
