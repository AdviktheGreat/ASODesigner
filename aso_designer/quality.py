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


# ---------------------------------------------------------------------------
# Transcript metadata extraction
# ---------------------------------------------------------------------------


@dataclass
class TranscriptMetadata:
    """Extracted metadata and region annotation for a transcript.

    When GenBank annotation is not available, regions are inferred
    from ORF detection: the longest ORF is treated as the CDS, and
    flanking sequences are classified as 5' and 3' UTR.
    """

    accession: str
    gene_name: str
    organism: str
    length: int
    chromosome: str | None = None
    strand: str | None = None
    five_utr_start: int | None = None
    five_utr_end: int | None = None
    five_utr_length: int = 0
    cds_start: int | None = None
    cds_end: int | None = None
    cds_length: int = 0
    three_utr_start: int | None = None
    three_utr_end: int | None = None
    three_utr_length: int = 0
    region_map: list[tuple[int, int, str]] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "accession": self.accession,
            "gene_name": self.gene_name,
            "organism": self.organism,
            "length": self.length,
            "chromosome": self.chromosome,
            "strand": self.strand,
            "five_utr_start": self.five_utr_start,
            "five_utr_end": self.five_utr_end,
            "five_utr_length": self.five_utr_length,
            "cds_start": self.cds_start,
            "cds_end": self.cds_end,
            "cds_length": self.cds_length,
            "three_utr_start": self.three_utr_start,
            "three_utr_end": self.three_utr_end,
            "three_utr_length": self.three_utr_length,
            "region_map": self.region_map,
        }


def extract_metadata(record: TranscriptRecord) -> TranscriptMetadata:
    """Extract structured metadata and region annotations from a transcript.

    Infers 5' UTR, CDS, and 3' UTR boundaries from ORF detection when
    GenBank annotation is not available. The longest ORF is treated as
    the CDS; flanking sequences are classified as UTRs.

    Parameters
    ----------
    record : TranscriptRecord
        The transcript to annotate.

    Returns
    -------
    TranscriptMetadata
        Structured metadata with region boundaries.
    """
    seq = record.sequence.upper()
    length = len(seq)

    cds = detect_cds(seq)

    meta = TranscriptMetadata(
        accession=record.accession,
        gene_name=record.gene_name,
        organism=record.organism,
        length=length,
        chromosome=record.chromosome,
        strand=record.strand,
    )

    if cds is not None:
        cds_start, cds_end = cds
        meta.cds_start = cds_start
        meta.cds_end = cds_end
        meta.cds_length = cds_end - cds_start

        # 5' UTR: everything before the CDS
        if cds_start > 0:
            meta.five_utr_start = 0
            meta.five_utr_end = cds_start
            meta.five_utr_length = cds_start

        # 3' UTR: everything after the CDS
        if cds_end < length:
            meta.three_utr_start = cds_end
            meta.three_utr_end = length
            meta.three_utr_length = length - cds_end

        # Build region map: (start, end, label) tuples
        if meta.five_utr_length > 0:
            meta.region_map.append((0, cds_start, "5UTR"))
        meta.region_map.append((cds_start, cds_end, "CDS"))
        if meta.three_utr_length > 0:
            meta.region_map.append((cds_end, length, "3UTR"))
    else:
        # No CDS found — entire transcript is non-coding
        meta.region_map.append((0, length, "noncoding"))

    logger.debug(
        "Extracted metadata for %s: CDS=%s, 5'UTR=%d nt, 3'UTR=%d nt",
        record.accession,
        f"{meta.cds_start}-{meta.cds_end}" if meta.cds_start is not None else "none",
        meta.five_utr_length,
        meta.three_utr_length,
    )

    return meta


def classify_position(metadata: TranscriptMetadata, position: int) -> str:
    """Classify a transcript position into a region (5UTR, CDS, 3UTR, noncoding).

    Parameters
    ----------
    metadata : TranscriptMetadata
        Transcript metadata with region annotations.
    position : int
        Zero-based position in the transcript.

    Returns
    -------
    str
        Region label: "5UTR", "CDS", "3UTR", or "noncoding".
    """
    for start, end, label in metadata.region_map:
        if start <= position < end:
            return label
    return "unknown"
