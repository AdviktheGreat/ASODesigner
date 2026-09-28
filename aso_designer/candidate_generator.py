"""ASO candidate generation via sliding window.

Generates ASO candidates across a target transcript using a configurable
sliding window. Inspired by GuideX's _generate_candidates method but
redesigned for ASO-specific needs: configurable length (15–25 nt),
step size, region targeting, and memory-managed batching for long
transcripts.
"""

from collections.abc import Generator, Sequence

from aso_designer.candidate import ASOCandidate
from aso_designer.config_schema import CandidateConfig
from aso_designer.quality import TranscriptMetadata, classify_position
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
            except ValueError:
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
        except ValueError:
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
