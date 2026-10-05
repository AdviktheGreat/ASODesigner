"""Shared sequence utilities for ASODesigner.

Provides validation, complement, reverse-complement, and GC content
functions used across modules. Designed from scratch for ASO-specific
needs (RNA/DNA handling, IUPAC ambiguity codes, chemistry-aware operations).
"""

from aso_designer.utils.exceptions import InvalidSequenceError

# IUPAC nucleotide codes
DNA_VALID: set[str] = set("ACGTN")
RNA_VALID: set[str] = set("ACGUN")
IUPAC_DNA: set[str] = set("ACGTNRYSWKMBDHV")
IUPAC_RNA: set[str] = set("ACGUNRYSWKMBDHV")

_COMPLEMENT_DNA = str.maketrans("ACGTNRYKMBWS", "TGCANRMVKVSW")
_COMPLEMENT_RNA = str.maketrans("ACGUNRYKMBWS", "UGCANRMVKVSW")


def validate_sequence(
    sequence: str,
    allow_rna: bool = True,
    allow_iupac: bool = False,
) -> str:
    """Validate and normalize a nucleotide sequence.

    Parameters
    ----------
    sequence : str
        Raw nucleotide sequence (case-insensitive).
    allow_rna : bool
        If True, accept U in addition to T.
    allow_iupac : bool
        If True, accept IUPAC ambiguity codes.

    Returns
    -------
    str
        Uppercase, validated sequence.

    Raises
    ------
    InvalidSequenceError
        If the sequence contains invalid characters or is empty.
    """
    if not sequence:
        raise InvalidSequenceError("Sequence is empty")

    seq = sequence.upper().strip()

    if allow_iupac:
        valid = IUPAC_DNA | IUPAC_RNA if allow_rna else IUPAC_DNA
    else:
        valid = DNA_VALID | RNA_VALID if allow_rna else DNA_VALID

    invalid_chars = set(seq) - valid
    if invalid_chars:
        raise InvalidSequenceError(f"Invalid characters in sequence: {sorted(invalid_chars)}")

    return seq


def complement(sequence: str, rna: bool = False) -> str:
    """Return the complement of a nucleotide sequence.

    Parameters
    ----------
    sequence : str
        DNA or RNA sequence (case-insensitive).
    rna : bool
        If True, treat as RNA (U complement).

    Returns
    -------
    str
        Complemented sequence (uppercase).
    """
    seq = sequence.upper()
    table = _COMPLEMENT_RNA if rna else _COMPLEMENT_DNA
    return seq.translate(table)


def reverse_complement(sequence: str, rna: bool = False) -> str:
    """Return the reverse complement of a nucleotide sequence.

    Parameters
    ----------
    sequence : str
        DNA or RNA sequence (case-insensitive).
    rna : bool
        If True, treat as RNA.

    Returns
    -------
    str
        Reverse-complemented sequence (uppercase).
    """
    return complement(sequence, rna=rna)[::-1]


def gc_content(sequence: str) -> float:
    """Calculate GC content as a fraction (0.0–1.0).

    Excludes ambiguous bases (N, R, S, etc.) from the denominator.

    Parameters
    ----------
    sequence : str
        Nucleotide sequence (case-insensitive).

    Returns
    -------
    float
        Fraction of G+C among non-ambiguous bases.
    """
    seq = sequence.upper()
    gc = seq.count("G") + seq.count("C")
    total = len(seq) - (
        seq.count("N")
        + seq.count("R")
        + seq.count("Y")
        + seq.count("S")
        + seq.count("W")
        + seq.count("K")
        + seq.count("M")
        + seq.count("B")
        + seq.count("D")
        + seq.count("H")
        + seq.count("V")
    )
    if total == 0:
        return 0.0
    return gc / total


def has_homopolymer(sequence: str, max_run: int = 4) -> bool:
    """Check if a sequence contains a homopolymer run longer than max_run.

    Parameters
    ----------
    sequence : str
        Nucleotide sequence (case-insensitive).
    max_run : int
        Maximum allowed consecutive identical bases.

    Returns
    -------
    bool
        True if a run longer than max_run is found.
    """
    seq = sequence.upper()
    run = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            run += 1
            if run > max_run:
                return True
        else:
            run = 1
    return False


def dna_to_rna(sequence: str) -> str:
    """Convert a DNA sequence to RNA (T → U)."""
    return sequence.upper().replace("T", "U")


def rna_to_dna(sequence: str) -> str:
    """Convert an RNA sequence to DNA (U → T)."""
    return sequence.upper().replace("U", "T")


def is_accession(target: str) -> bool:
    """Check if a target string looks like an NCBI accession.

    Accessions match patterns like NM_001123, NR_123456, XM_456789.1.

    Parameters
    ----------
    target : str
        Target string.

    Returns
    -------
    bool
        True if the target looks like an accession.
    """
    import re

    return bool(re.match(r"^[A-Z]{1,4}_\d+(\.\d+)?$", target.strip().upper()))
