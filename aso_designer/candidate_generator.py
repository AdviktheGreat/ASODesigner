"""ASO candidate generation via sliding window.

Generates ASO candidates across a target transcript using a configurable
sliding window. Inspired by GuideX's _generate_candidates method but
redesigned for ASO-specific needs: configurable length (15–25 nt),
step size, region targeting, and memory-managed batching for long
transcripts.
"""

from collections.abc import Generator, Sequence
from dataclasses import dataclass

from aso_designer.candidate import ASOCandidate
from aso_designer.config_schema import CandidateConfig
from aso_designer.quality import TranscriptMetadata, classify_position
from aso_designer.utils.exceptions import InvalidSequenceError
from aso_designer.utils.logger import get_logger
from aso_designer.utils.sequence import gc_content, validate_sequence

logger = get_logger(__name__)


def generate_candidates(
    sequence: str,
    config: CandidateConfig,
    metadata: TranscriptMetadata | None = None,
    regions: Sequence[tuple[int, int]] | None = None,
) -> list[ASOCandidate]:
    """Generate ASO candidates from a transcript sequence.

    Slides a window of configurable length across the transcript (or
    specified regions), creating an ASOCandidate for each position.
    Candidates with invalid characters are skipped.

    Parameters
    ----------
    sequence : str
        Target transcript sequence (DNA, uppercase or lowercase).
    config : CandidateConfig
        Candidate generation parameters (length, step, chemistry).
    metadata : TranscriptMetadata, optional
        Transcript metadata for region annotation. If provided, each
        candidate's ``region`` field is set based on its position.
    regions : Sequence of (int, int), optional
        Specific regions to generate candidates within. If None,
        candidates are generated across the entire transcript. Each
        tuple is (start, end) in 0-based transcript coordinates.

    Returns
    -------
    list[ASOCandidate]
        All generated candidates (before filtering).
    """
    seq = sequence.upper()
    seq_len = len(seq)
    aso_len = config.length
    step = config.step

    if aso_len > seq_len:
        logger.warning(
            "ASO length (%d) exceeds transcript length (%d) — no candidates",
            aso_len,
            seq_len,
        )
        return []

    # Determine the regions to scan
    if regions is not None:
        scan_regions = list(regions)
    else:
        scan_regions = [(0, seq_len)]

    candidates: list[ASOCandidate] = []

    for region_start, region_end in scan_regions:
        # Validate region boundaries
        if region_start < 0 or region_end > seq_len or region_start >= region_end:
            logger.warning("Skipping invalid region: (%d, %d)", region_start, region_end)
            continue

        region_len = region_end - region_start
        max_start = region_len - aso_len
        if max_start < 0:
            logger.debug(
                "Region (%d, %d) too short for ASO length %d",
                region_start,
                region_end,
                aso_len,
            )
            continue

        for i in range(0, max_start + 1, step):
            abs_start = region_start + i
            abs_end = abs_start + aso_len
            subseq = seq[abs_start:abs_end]

            # Skip sequences with invalid characters
            try:
                subseq = validate_sequence(subseq, allow_rna=False, allow_iupac=False)
            except InvalidSequenceError:
                logger.debug("Skipping candidate at %d: invalid characters", abs_start)
                continue

            # Compute GC content
            gc = gc_content(subseq)

            # Determine region label
            region_label = "unknown"
            if metadata is not None:
                region_label = classify_position(metadata, abs_start)

            candidate = ASOCandidate(
                sequence=subseq,
                start=abs_start,
                end=abs_end,
                length=aso_len,
                chemistry=config.chemistry,
                gc_content=gc,
                region=region_label,
                passes_filters=True,  # Filtering happens later
            )
            candidates.append(candidate)

    logger.info(
        "Generated %d ASO candidates (length=%d, step=%d, regions=%d)",
        len(candidates),
        aso_len,
        step,
        len(scan_regions),
    )
    return candidates


def generate_candidates_batched(
    sequence: str,
    config: CandidateConfig,
    metadata: TranscriptMetadata | None = None,
    batch_size: int = 10000,
) -> Generator[list[ASOCandidate], None, None]:
    """Generate ASO candidates in memory-managed batches.

    For very long transcripts (e.g. DMD at ~14 kb), generating all
    candidates at once can use significant memory. This generator
    yields candidates in batches.

    Parameters
    ----------
    sequence : str
        Target transcript sequence.
    config : CandidateConfig
        Candidate generation parameters.
    metadata : TranscriptMetadata, optional
        Transcript metadata for region annotation.
    batch_size : int
        Number of candidates per batch (default 10,000).

    Yields
    ------
    list[ASOCandidate]
        Batches of candidates.
    """
    seq = sequence.upper()
    seq_len = len(seq)
    aso_len = config.length
    step = config.step

    if aso_len > seq_len:
        return

    batch: list[ASOCandidate] = []
    max_start = seq_len - aso_len

    for i in range(0, max_start + 1, step):
        abs_start = i
        abs_end = abs_start + aso_len
        subseq = seq[abs_start:abs_end]

        try:
            subseq = validate_sequence(subseq, allow_rna=False, allow_iupac=False)
        except InvalidSequenceError:
            continue

        gc = gc_content(subseq)
        region_label = "unknown"
        if metadata is not None:
            region_label = classify_position(metadata, abs_start)

        candidate = ASOCandidate(
            sequence=subseq,
            start=abs_start,
            end=abs_end,
            length=aso_len,
            chemistry=config.chemistry,
            gc_content=gc,
            region=region_label,
            passes_filters=True,
        )
        batch.append(candidate)

        if len(batch) >= batch_size:
            yield batch
            batch = []

    if batch:
        yield batch


# ---------------------------------------------------------------------------
# Memory management utilities
# ---------------------------------------------------------------------------

# Estimated memory per ASOCandidate (bytes) — empirically measured
# from the dataclass with typical feature/score dicts populated.
_ESTIMATED_BYTES_PER_CANDIDATE = 512


def count_candidates(
    sequence: str,
    config: CandidateConfig,
    regions: Sequence[tuple[int, int]] | None = None,
) -> int:
    """Count the number of candidates that would be generated.

    Useful for estimating memory usage before generation. Does not
    actually create candidate objects.

    Parameters
    ----------
    sequence : str
        Target transcript sequence.
    config : CandidateConfig
        Candidate generation parameters.
    regions : Sequence of (int, int), optional
        Specific regions to count within. If None, counts across the
        entire transcript.

    Returns
    -------
    int
        Number of candidates that would be generated.
    """
    seq_len = len(sequence)
    aso_len = config.length
    step = config.step

    if aso_len > seq_len:
        return 0

    if regions is not None:
        scan_regions = list(regions)
    else:
        scan_regions = [(0, seq_len)]

    total = 0
    for region_start, region_end in scan_regions:
        if region_start < 0 or region_end > seq_len or region_start >= region_end:
            continue
        region_len = region_end - region_start
        max_start = region_len - aso_len
        if max_start < 0:
            continue
        total += (max_start // step) + 1

    return total


def estimate_memory_mb(
    sequence: str,
    config: CandidateConfig,
    regions: Sequence[tuple[int, int]] | None = None,
) -> float:
    """Estimate memory usage for candidate generation in megabytes.

    Parameters
    ----------
    sequence : str
        Target transcript sequence.
    config : CandidateConfig
        Candidate generation parameters.
    regions : Sequence of (int, int), optional
        Specific regions.

    Returns
    -------
    float
        Estimated memory in MB.
    """
    n = count_candidates(sequence, config, regions)
    return (n * _ESTIMATED_BYTES_PER_CANDIDATE) / (1024 * 1024)


def should_use_batching(
    sequence: str,
    config: CandidateConfig,
    memory_threshold_mb: float = 500.0,
    regions: Sequence[tuple[int, int]] | None = None,
) -> bool:
    """Determine whether batched generation is recommended.

    Parameters
    ----------
    sequence : str
        Target transcript sequence.
    config : CandidateConfig
        Candidate generation parameters.
    memory_threshold_mb : float
        Memory threshold above which batching is recommended (default 500 MB).
    regions : Sequence of (int, int), optional
        Specific regions.

    Returns
    -------
    bool
        True if estimated memory exceeds the threshold.
    """
    return estimate_memory_mb(sequence, config, regions) > memory_threshold_mb


# ---------------------------------------------------------------------------
# Exon-aware coordinate mapping
# ---------------------------------------------------------------------------


@dataclass
class ExonMap:
    """Maps transcript coordinates to exon coordinates.

    Stores exon boundaries in transcript coordinates (0-based) and
    provides utilities to determine which exon a position falls in,
    and whether a candidate spans an exon-exon junction.

    Attributes
    ----------
    exons : list[tuple[int, int]]
        List of (start, end) exon boundaries in transcript coordinates.
        Each exon is [start, end) — start inclusive, end exclusive.
    """

    exons: list[tuple[int, int]]

    def __post_init__(self) -> None:
        # Sort exons by start position and validate no overlaps
        self.exons.sort(key=lambda e: e[0])
        for i in range(1, len(self.exons)):
            if self.exons[i][0] < self.exons[i - 1][1]:
                raise ValueError(f"Overlapping exons: {self.exons[i - 1]} and {self.exons[i]}")

    def exon_index(self, position: int) -> int | None:
        """Return the exon index (0-based) containing a position, or None.

        Parameters
        ----------
        position : int
            Transcript coordinate (0-based).

        Returns
        -------
        int or None
            Exon index, or None if the position is in an intron.
        """
        for i, (start, end) in enumerate(self.exons):
            if start <= position < end:
                return i
        return None

    def spans_junction(self, candidate: ASOCandidate) -> bool:
        """Check if a candidate spans an exon-exon junction.

        A candidate spans a junction if its start and end fall in
        different exons (or if it crosses an intron).

        Parameters
        ----------
        candidate : ASOCandidate
            The ASO candidate to check.

        Returns
        -------
        bool
            True if the candidate spans at least one exon junction.
        """
        start_exon = self.exon_index(candidate.start)
        end_exon = self.exon_index(candidate.end - 1)
        return start_exon is not None and end_exon is not None and start_exon != end_exon

    def junction_positions(self, candidate: ASOCandidate) -> list[int]:
        """Return all exon junction positions within a candidate's span.

        Parameters
        ----------
        candidate : ASOCandidate
            The ASO candidate to check.

        Returns
        -------
        list[int]
            Transcript positions where exon-exon junctions occur within
            the candidate's [start, end) range.
        """
        junctions: list[int] = []
        for i in range(1, len(self.exons)):
            junction_pos = self.exons[i][0]  # Start of exon i = end of exon i-1
            if candidate.start < junction_pos < candidate.end:
                junctions.append(junction_pos)
        return junctions

    def annotate_candidate(self, candidate: ASOCandidate) -> None:
        """Annotate a candidate with exon and junction information.

        Adds the following to the candidate's features dict:
        - "exon_index": exon index of the candidate start (or None)
        - "spans_junction": whether the candidate spans a junction
        - "junction_positions": list of junction positions within the candidate

        Parameters
        ----------
        candidate : ASOCandidate
            The candidate to annotate (modified in place).
        """
        candidate.features["exon_index"] = self.exon_index(candidate.start)
        candidate.features["spans_junction"] = self.spans_junction(candidate)
        candidate.features["junction_positions"] = self.junction_positions(candidate)


def annotate_with_exons(
    candidates: list[ASOCandidate],
    exon_map: ExonMap,
) -> None:
    """Annotate all candidates with exon and junction information.

    Parameters
    ----------
    candidates : list[ASOCandidate]
        Candidates to annotate (modified in place).
    exon_map : ExonMap
        Exon boundary map.
    """
    junction_count = 0
    for candidate in candidates:
        exon_map.annotate_candidate(candidate)
        if candidate.features.get("spans_junction"):
            junction_count += 1

    logger.info(
        "Annotated %d candidates with exon map (%d span junctions)",
        len(candidates),
        junction_count,
    )
