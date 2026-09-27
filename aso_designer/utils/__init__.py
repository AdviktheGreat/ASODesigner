"""Utility modules for ASODesigner (logging, exceptions, sequence helpers)."""

from aso_designer.utils.exceptions import (
    ASODesignError,
    ConfigurationError,
    InvalidSequenceError,
    OffTargetError,
    ReportError,
    ScoringError,
    ThermodynamicsError,
    TranscriptFetchError,
    VariantError,
    ViennaRNAError,
)
from aso_designer.utils.logger import get_logger, setup_logger
from aso_designer.utils.sequence import (
    complement,
    dna_to_rna,
    gc_content,
    has_homopolymer,
    reverse_complement,
    rna_to_dna,
    validate_sequence,
)

__all__ = [
    # Logger
    "setup_logger",
    "get_logger",
    # Exceptions
    "ASODesignError",
    "TranscriptFetchError",
    "InvalidSequenceError",
    "ConfigurationError",
    "ViennaRNAError",
    "ThermodynamicsError",
    "OffTargetError",
    "VariantError",
    "ScoringError",
    "ReportError",
    # Sequence
    "validate_sequence",
    "complement",
    "reverse_complement",
    "gc_content",
    "has_homopolymer",
    "dna_to_rna",
    "rna_to_dna",
]
