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

# Default minimum ORF length in codons (30 codons = 90 nt) used by
# find_orfs and CDS detection.
_DEFAULT_ORF_MIN_CODONS = 30


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
        """Serialize quality report to a dictionary."""
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


def find_orfs(sequence: str, min_length: int = _DEFAULT_ORF_MIN_CODONS) -> list[tuple[int, int]]:
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
                # Scan for the first downstream stop codon
                j = i + 3
                stop_found = False
                while j < seq_len - 2:
                    if sequence[j : j + 3] in _STOP_CODONS:
                        stop_found = True
                        break
                    j += 3
                if stop_found:
                    orf_len = (j + 3 - i) // 3
                    if orf_len >= min_length:
                        orfs.append((i, j + 3))
                    # Resume scanning immediately after the stop codon
                    # (the outer loop must not advance i a second time,
                    # which previously skipped one codon per ORF)
                    i = j + 3
                    continue
                # No stop codon downstream: no further ORFs in this frame
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
    orfs = find_orfs(sequence, min_length=_DEFAULT_ORF_MIN_CODONS)
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
    orfs = find_orfs(seq, min_length=_DEFAULT_ORF_MIN_CODONS)
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
        """Serialize transcript metadata to a dictionary."""
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


# ---------------------------------------------------------------------------
# Multi-isoform support
# ---------------------------------------------------------------------------


@dataclass
class TranscriptSet:
    """A collection of transcript isoforms for a single gene.

    Stores multiple transcripts (e.g. RefSeq NM_ accessions for the
    same gene), provides comparison utilities, and allows selection
    of a representative transcript for ASO design.

    Attributes
    ----------
    gene_name : str
        Gene symbol shared by all transcripts.
    organism : str
        Organism name.
    transcripts : list[TranscriptRecord]
        All isoform transcripts.
    metadata : list[TranscriptMetadata]
        Extracted metadata for each transcript (populated on first access).
    """

    gene_name: str
    organism: str
    transcripts: list[TranscriptRecord] = field(default_factory=list)
    _metadata: list[TranscriptMetadata] | None = None

    def add_transcript(self, record: TranscriptRecord) -> None:
        """Add a transcript isoform to the set."""
        if record.gene_name and record.gene_name != self.gene_name:
            logger.warning(
                "Adding transcript with gene '%s' to set for gene '%s'",
                record.gene_name,
                self.gene_name,
            )
        self.transcripts.append(record)
        self._metadata = None  # Invalidate cached metadata

    def get_metadata(self) -> list[TranscriptMetadata]:
        """Get metadata for all transcripts (computed once, then cached)."""
        if self._metadata is None:
            self._metadata = [extract_metadata(r) for r in self.transcripts]
        return self._metadata

    def longest_transcript(self) -> TranscriptRecord | None:
        """Return the longest transcript isoform.

        This is the default selection for ASO design since longer
        transcripts provide more candidate sites and typically include
        all exon junctions.

        Returns
        -------
        TranscriptRecord or None
            The longest transcript, or None if the set is empty.
        """
        if not self.transcripts:
            return None
        return max(self.transcripts, key=lambda r: r.length)

    def longest_cds_transcript(self) -> TranscriptRecord | None:
        """Return the transcript with the longest CDS.

        Useful when ASO design should target coding regions.

        Returns
        -------
        TranscriptRecord or None
            The transcript with the longest detected CDS, or None if empty.
        """
        if not self.transcripts:
            return None
        metas = self.get_metadata()
        best_idx = max(range(len(metas)), key=lambda i: metas[i].cds_length)
        return self.transcripts[best_idx]

    def shared_regions(
        self, min_overlap: int = 20, reference: "TranscriptRecord | None" = None
    ) -> list[tuple[int, int]]:
        """Find regions shared across all isoforms (by sequence alignment).

        Uses a simple prefix/suffix matching heuristic to identify
        regions present in all transcripts. This is a simplified
        approach — for rigorous isoform comparison, use an external
        alignment tool.

        Parameters
        ----------
        min_overlap : int
            Minimum shared region length to report (default 20 nt).
        reference : TranscriptRecord or None
            Transcript providing the coordinate frame. Defaults to the
            longest transcript; pass the design transcript so returned
            coordinates line up with candidate positions on it.

        Returns
        -------
        list[tuple[int, int]]
            List of (start, end) shared regions in the reference
            transcript's coordinate system.
        """
        if len(self.transcripts) < 2:
            base = reference or self.longest_transcript()
            if base:
                return [(0, base.length)]
            return []

        base = reference or self.longest_transcript()
        if base is None:
            return []

        ref_seq = base.sequence
        window = min_overlap
        base_idx = next(i for i, t in enumerate(self.transcripts) if t is base)
        other_seqs = [t.sequence for i, t in enumerate(self.transcripts) if i != base_idx]

        # Precompute the set of all min_overlap-mers for each other
        # transcript: per-window membership checks become O(1) hash
        # lookups instead of O(L·m) substring scans, dropping overall
        # cost from O(L²·n·m) to O(L·n) after an O(L·n) preprocessing
        # pass.
        kmer_sets = [
            {seq[j : j + window] for j in range(len(seq) - window + 1)} for seq in other_seqs
        ]

        shared: list[tuple[int, int]] = []
        for i in range(0, len(ref_seq) - window + 1):
            kmer = ref_seq[i : i + window]
            if all(kmer in ks for ks in kmer_sets):
                # Extend the shared region
                if shared and shared[-1][1] == i:
                    shared[-1] = (shared[-1][0], i + window)
                else:
                    shared.append((i, i + window))

        return shared

    def summary(self) -> dict:
        """Return a summary of the isoform set.

        Returns
        -------
        dict
            Summary with isoform count, lengths, and CDS info.
        """
        metas = self.get_metadata()
        return {
            "gene_name": self.gene_name,
            "organism": self.organism,
            "isoform_count": len(self.transcripts),
            "accessions": [r.accession for r in self.transcripts],
            "lengths": [r.length for r in self.transcripts],
            "cds_lengths": [m.cds_length for m in metas],
            "longest_accession": self.longest_transcript().accession if self.transcripts else None,
            "longest_length": self.longest_transcript().length if self.transcripts else 0,
        }


def isoform_coverage(shared_regions: list[tuple[int, int]], start: int, end: int) -> float:
    """Fraction of [start, end) covered by isoform-shared regions.

    Args:
        shared_regions: (start, end) tuples in the design transcript's
            coordinate system (see TranscriptSet.shared_regions).
        start: candidate start (0-based).
        end: candidate end (exclusive).

    Returns:
        Coverage fraction in [0.0, 1.0]; 0.0 when the window is empty.
    """
    if end <= start:
        return 0.0
    length = end - start
    covered = 0
    for r_start, r_end in shared_regions:
        lo = max(r_start, start)
        hi = min(r_end, end)
        if hi > lo:
            covered += hi - lo
    return min(1.0, covered / length)


# ---------------------------------------------------------------------------
# Batch quality report generation
# ---------------------------------------------------------------------------


@dataclass
class BatchQualityReport:
    """Aggregate quality assessment across multiple transcripts.

    Provides summary statistics, per-transcript reports, and flags
    transcripts that may be problematic for ASO design.
    """

    total_count: int
    valid_count: int
    invalid_count: int
    warning_count: int
    min_length: int
    max_length: int
    mean_length: float
    total_ambiguous: int
    transcripts_with_cds: int
    transcripts_without_cds: int
    per_transcript: list[QualityReport] = field(default_factory=list)
    flagged: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        """Serialize batch quality report to a dictionary."""
        return {
            "total_count": self.total_count,
            "valid_count": self.valid_count,
            "invalid_count": self.invalid_count,
            "warning_count": self.warning_count,
            "min_length": self.min_length,
            "max_length": self.max_length,
            "mean_length": round(self.mean_length, 1),
            "total_ambiguous": self.total_ambiguous,
            "transcripts_with_cds": self.transcripts_with_cds,
            "transcripts_without_cds": self.transcripts_without_cds,
            "per_transcript": [r.to_dict() for r in self.per_transcript],
            "flagged": self.flagged,
        }


def generate_quality_report(
    records: list[TranscriptRecord],
    min_length: int = 200,
    max_ambiguous_fraction: float = 0.05,
    require_cds: bool = False,
) -> BatchQualityReport:
    """Generate a batch quality report across multiple transcripts.

    Validates each transcript, aggregates statistics, and flags
    transcripts that may be problematic for ASO design.

    Parameters
    ----------
    records : list[TranscriptRecord]
        Transcript records to assess.
    min_length : int
        Minimum acceptable sequence length (default 200 nt).
    max_ambiguous_fraction : float
        Maximum fraction of ambiguous bases (default 5%).
    require_cds : bool
        If True, a CDS must be detected for validity.

    Returns
    -------
    BatchQualityReport
        Aggregate quality assessment with per-transcript details.
    """
    reports: list[QualityReport] = []
    flagged: list[str] = []

    for record in records:
        report = validate_transcript(
            record,
            min_length=min_length,
            max_ambiguous_fraction=max_ambiguous_fraction,
            require_cds=require_cds,
        )
        reports.append(report)

        # Flag transcripts with specific issues
        if not report.is_valid:
            flagged.append(f"{record.accession}: {'; '.join(report.errors)}")
        elif report.warnings:
            for warning in report.warnings:
                if "ambiguous" in warning.lower():
                    flagged.append(f"{record.accession}: {warning}")
                elif "No CDS" in warning:
                    flagged.append(f"{record.accession}: {warning}")

    lengths = [r.length for r in records] if records else [0]
    valid_count = sum(1 for r in reports if r.is_valid)
    warning_count = sum(1 for r in reports if r.is_valid and r.warnings)
    total_ambiguous = sum(r.ambiguous_count for r in reports)
    with_cds = sum(1 for r in reports if r.has_cds)

    batch = BatchQualityReport(
        total_count=len(records),
        valid_count=valid_count,
        invalid_count=len(records) - valid_count,
        warning_count=warning_count,
        min_length=min(lengths) if lengths else 0,
        max_length=max(lengths) if lengths else 0,
        mean_length=sum(lengths) / len(lengths) if lengths else 0.0,
        total_ambiguous=total_ambiguous,
        transcripts_with_cds=with_cds,
        transcripts_without_cds=len(records) - with_cds,
        per_transcript=reports,
        flagged=flagged,
    )

    logger.info(
        "Batch quality report: %d/%d valid, %d with warnings, %d flagged",
        valid_count,
        len(records),
        warning_count,
        len(flagged),
    )

    return batch
