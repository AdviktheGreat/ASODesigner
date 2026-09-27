"""Transcript quality validation for ASODesigner.

Validates fetched transcript sequences for ambiguous bases, minimum
length, and CDS (coding sequence) detection. Designed from scratch —
GuideX has basic N-count checks; we add ORF detection, CDS annotation,
and structured quality reporting.
"""

from dataclasses import dataclass, field

from aso_designer.transcript_fetcher import TranscriptRecord
from aso_designer.utils.logger import get_logger
from aso_designer.utils.sequence import validate_sequence

logger = get_logger(__name__)

# Standard genetic code (codon table 1)
_START_CODONS = {"ATG"}
_STOP_CODONS = {"TAA", "TAG", "TGA"}


@dataclass
class QualityReport:
    """Structured quality assessment for a transcript."""

    accession: str
    is_valid: bool
    length: int
    ambiguous_count: int
    ambiguous_fraction: float
    has_cds: bool
    cds_start: int | None = None
    cds_end: int | None = None
    cds_length: int = 0
    orf_count: int = 0
    warnings: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "accession": self.accession,
            "is_valid": self.is_valid,
            "length": self.length,
            "ambiguous_count": self.ambiguous_count,
            "ambiguous_fraction": round(self.ambiguous_fraction, 4),
            "has_cds": self.has_cds,
            "cds_start": self.cds_start,
            "cds_end": self.cds_end,
            "cds_length": self.cds_length,
            "orf_count": self.orf_count,
            "warnings": self.warnings,
            "errors": self.errors,
        }


def find_orfs(sequence: str, min_length: int = 30) -> list[tuple[int, int]]:
    """Find open reading frames in a DNA sequence (all three forward frames).

    Parameters
    ----------
    sequence : str
        DNA sequence (uppercase).
    min_length : int
        Minimum ORF length in codons (default 30 = 90 nt).

    Returns
    -------
    List[tuple[int, int]]
        List of (start, end) positions for each ORF found.
    """
    orfs: list[tuple[int, int]] = []
    seq_len = len(sequence)

    for frame in range(3):
        i = frame
        while i < seq_len - 2:
            codon = sequence[i : i + 3]
            if codon in _START_CODONS:
                # Scan for stop codon
                j = i + 3
                while j < seq_len - 2:
                    codon_j = sequence[j : j + 3]
                    if codon_j in _STOP_CODONS:
                        orf_len = (j + 3 - i) // 3
                        if orf_len >= min_length:
                            orfs.append((i, j + 3))
                        i = j + 3
                        break
                    j += 3
                else:
                    # No stop found before end of sequence
                    break
            i += 3

    return orfs


def detect_cds(sequence: str) -> tuple[int, int] | None:
    """Detect the most likely CDS in a transcript.

    Finds the longest ORF as a proxy for the CDS when annotation
    is not available from the source.

    Parameters
    ----------
    sequence : str
        DNA sequence (uppercase).

    Returns
    -------
    tuple[int, int] or None
        (start, end) of the longest ORF, or None if no ORF found.
    """
    orfs = find_orfs(sequence, min_length=30)
    if not orfs:
        return None
    # Return the longest ORF
    return max(orfs, key=lambda x: x[1] - x[0])


def validate_transcript(
    record: TranscriptRecord,
    min_length: int = 200,
    max_ambiguous_fraction: float = 0.05,
    require_cds: bool = False,
) -> QualityReport:
    """Validate a transcript record and generate a quality report.

    Parameters
    ----------
    record : TranscriptRecord
        The transcript to validate.
    min_length : int
        Minimum acceptable sequence length (default 200 nt).
    max_ambiguous_fraction : float
        Maximum fraction of ambiguous bases (default 5%).
    require_cds : bool
        If True, a CDS must be detected for the transcript to be valid.

    Returns
    -------
    QualityReport
        Structured quality assessment.
    """
    seq = record.sequence.upper()
    length = len(seq)
    warnings: list[str] = []
    errors: list[str] = []

    # Ambiguous base count
    ambiguous_chars = set("NRYSWKMBDHV")
    ambiguous_count = sum(1 for c in seq if c in ambiguous_chars)
    ambiguous_fraction = ambiguous_count / length if length > 0 else 1.0

    # Length check
    if length < min_length:
        errors.append(f"Sequence too short: {length} nt < {min_length} nt minimum")

    # Ambiguous base check
    if ambiguous_fraction > max_ambiguous_fraction:
        errors.append(
            f"Too many ambiguous bases: {ambiguous_count} ({ambiguous_fraction:.1%}) "
            f"> {max_ambiguous_fraction:.1%} threshold"
        )
    elif ambiguous_count > 0:
        warnings.append(f"{ambiguous_count} ambiguous bases ({ambiguous_fraction:.1%})")

    # CDS detection
    orfs = find_orfs(seq, min_length=30)
    cds = detect_cds(seq)
    has_cds = cds is not None

    if require_cds and not has_cds:
        errors.append("No CDS detected (required by configuration)")

    if has_cds:
        cds_start, cds_end = cds
        cds_length = cds_end - cds_start
        if cds_length < 100:
            warnings.append(f"Short CDS: {cds_length} nt")
    else:
        cds_start, cds_end, cds_length = None, None, 0
        warnings.append("No CDS detected — transcript may be non-coding")

    # Sequence character validation
    try:
        validate_sequence(seq, allow_rna=False, allow_iupac=True)
    except Exception as e:
        errors.append(f"Invalid sequence characters: {e}")

    is_valid = len(errors) == 0

    report = QualityReport(
        accession=record.accession,
        is_valid=is_valid,
        length=length,
        ambiguous_count=ambiguous_count,
        ambiguous_fraction=ambiguous_fraction,
        has_cds=has_cds,
        cds_start=cds_start,
        cds_end=cds_end,
        cds_length=cds_length,
        orf_count=len(orfs),
        warnings=warnings,
        errors=errors,
    )

    if not is_valid:
        logger.warning("Transcript %s failed validation: %s", record.accession, errors)
    elif warnings:
        logger.info("Transcript %s passed with warnings: %s", record.accession, warnings)

    return report
