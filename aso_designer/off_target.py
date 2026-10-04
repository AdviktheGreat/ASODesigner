"""Off-target analysis for ASO candidates.

Searches the transcriptome for potential off-target binding sites
using a seed-and-extend strategy with a k-mer index. Scores candidates
based on off-target hit count, mismatch profile, and seed region
matches.

Scientific note: Off-target predictions are computational estimates
for candidate ranking. The seed-and-extend approach finds sites with
exact seed matches; sites with seed mismatches may be missed. This is
a deliberate trade-off for search speed at transcriptome scale.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from aso_designer.utils.exceptions import OffTargetError
from aso_designer.utils.logger import get_logger
from aso_designer.utils.sequence import reverse_complement, validate_sequence

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Transcriptome reference
# ---------------------------------------------------------------------------


@dataclass
class TranscriptEntry:
    """A single transcript in the reference transcriptome.

    Attributes
    ----------
    transcript_id : str
        Transcript accession (e.g., NM_001123.1).
    gene_name : str
        Gene symbol.
    sequence : str
        Transcript sequence (DNA alphabet, uppercase).
    length : int
        Sequence length in nucleotides.
    """

    transcript_id: str
    gene_name: str
    sequence: str
    length: int

    @classmethod
    def from_record(cls, record) -> TranscriptEntry:
        """Build from a TranscriptRecord."""
        return cls(
            transcript_id=record.accession,
            gene_name=record.gene_name,
            sequence=record.sequence,
            length=len(record.sequence),
        )


class TranscriptomeReference:
    """A searchable transcriptome reference for off-target analysis.

    Stores transcript sequences and builds a k-mer index for fast
    seed-based searching. Inspired by GuideX's transcriptome handling
    but redesigned with a k-mer index for efficient ASO-length queries.

    Parameters
    ----------
    kmer_size : int
        K-mer size for the index (default 12). Larger k = fewer false
        positives but more memory; smaller k = more candidates to verify.
    """

    def __init__(self, kmer_size: int = 12):
        self.kmer_size = kmer_size
        self._transcripts: list[TranscriptEntry] = []
        self._index: dict[str, list[tuple[int, int]]] = defaultdict(list)
        # _index maps k-mer → list of (transcript_idx, position)
        self._indexed = False
        self._total_length: int = 0

    @property
    def num_transcripts(self) -> int:
        """Number of transcripts in the reference."""
        return len(self._transcripts)

    @property
    def total_length(self) -> int:
        """Total length of all transcript sequences combined."""
        return self._total_length

    def add_transcript(self, entry: TranscriptEntry) -> None:
        """Add a single transcript to the reference.

        Parameters
        ----------
        entry : TranscriptEntry
            Transcript to add.

        Raises
        ------
        OffTargetError
            If the sequence is invalid.
        """
        try:
            seq = validate_sequence(entry.sequence, allow_rna=False, allow_iupac=True)
        except Exception as e:
            raise OffTargetError(f"Invalid sequence for {entry.transcript_id}: {e}") from e

        self._transcripts.append(
            TranscriptEntry(
                transcript_id=entry.transcript_id,
                gene_name=entry.gene_name,
                sequence=seq,
                length=len(seq),
            )
        )
        self._total_length += len(seq)
        self._indexed = False  # Invalidate index

    def add_transcripts(self, entries: list[TranscriptEntry]) -> None:
        """Add multiple transcripts."""
        for entry in entries:
            self.add_transcript(entry)

    def add_from_records(self, records: list) -> None:
        """Add transcripts from TranscriptRecord objects.

        Parameters
        ----------
        records : list
            List of TranscriptRecord objects (from transcript_fetcher
            or local_loader).
        """
        for record in records:
            self.add_transcript(TranscriptEntry.from_record(record))

    def load_from_fasta(self, file_path: str | Path) -> None:
        """Load transcripts from a FASTA file.

        Uses the local_loader to parse the file and adds all valid
        transcripts to the reference.

        Parameters
        ----------
        file_path : str or Path
            Path to the FASTA file.
        """
        from aso_designer.local_loader import load_fasta

        records = load_fasta(file_path, min_length=50)
        self.add_from_records(records)
        logger.info(
            "Loaded %d transcripts from %s (total %d nt)",
            len(records),
            file_path,
            self._total_length,
        )

    def build_index(self) -> None:
        """Build the k-mer index from all transcript sequences.

        Must be called before searching. Automatically called by
        the search engine if the index is stale.
        """
        self._index.clear()
        k = self.kmer_size

        for tidx, transcript in enumerate(self._transcripts):
            seq = transcript.sequence
            for i in range(len(seq) - k + 1):
                kmer = seq[i : i + k]
                self._index[kmer].append((tidx, i))

        self._indexed = True
        logger.info(
            "Built k-mer index: %d transcripts, %d unique k-mers (k=%d)",
            len(self._transcripts),
            len(self._index),
            k,
        )

    def ensure_index(self) -> None:
        """Build the index if it hasn't been built yet."""
        if not self._indexed:
            self.build_index()

    def get_transcript(self, index: int) -> TranscriptEntry:
        """Get a transcript by its index in the reference."""
        if index < 0 or index >= len(self._transcripts):
            raise OffTargetError(f"Transcript index {index} out of range")
        return self._transcripts[index]

    def get_sequence(self, transcript_idx: int, start: int, end: int) -> str:
        """Extract a subsequence from a transcript.

        Parameters
        ----------
        transcript_idx : int
            Index of the transcript in the reference.
        start : int
            Start position (0-based, inclusive).
        end : int
            End position (0-based, exclusive).

        Returns
        -------
        str
            Subsequence.
        """
        transcript = self._transcripts[transcript_idx]
        return transcript.sequence[start:end]

    def summary(self) -> dict:
        """Return a summary of the reference."""
        return {
            "num_transcripts": self.num_transcripts,
            "total_length": self._total_length,
            "kmer_size": self.kmer_size,
            "indexed": self._indexed,
            "unique_kmers": len(self._index) if self._indexed else 0,
        }


# ---------------------------------------------------------------------------
# Off-target hit
# ---------------------------------------------------------------------------


@dataclass
class OffTargetHit:
    """A single off-target binding site.

    Attributes
    ----------
    transcript_id : str
        Transcript where the hit was found.
    gene_name : str
        Gene name of the transcript.
    position : int
        Start position of the hit in the transcript (0-based).
    mismatches : int
        Total number of mismatches in the full-length alignment.
    seed_mismatches : int
        Number of mismatches in the seed region.
    aligned_sequence : str
        The transcript subsequence at the hit position.
    is_on_target : bool
        True if this is the intended target site (not an off-target).
    """

    transcript_id: str
    gene_name: str
    position: int
    mismatches: int
    seed_mismatches: int
    aligned_sequence: str
    is_on_target: bool = False


# ---------------------------------------------------------------------------
# Search engine
# ---------------------------------------------------------------------------


def _hamming_distance(seq1: str, seq2: str) -> int:
    """Compute Hamming distance between two equal-length sequences."""
    return sum(a != b for a, b in zip(seq1, seq2, strict=True))


class OffTargetSearchEngine:
    """Searches a transcriptome reference for ASO off-target binding sites.

    Uses a seed-and-extend strategy:
    1. Compute the reverse complement of the ASO (target binding site)
    2. Extract k-mers from the target site and look up in the index
    3. For each candidate position, verify the full-length match
    4. Count total and seed-region mismatches

    Parameters
    ----------
    reference : TranscriptomeReference
        The indexed transcriptome to search against.
    max_mismatches : int
        Maximum allowed mismatches for a hit (default 3).
    seed_start : int
        Seed region start position in the ASO (0-based, default 2).
    seed_end : int
        Seed region end position in the ASO (exclusive, default 8).
    """

    def __init__(
        self,
        reference: TranscriptomeReference,
        max_mismatches: int = 3,
        seed_start: int = 2,
        seed_end: int = 8,
    ):
        self.reference = reference
        self.max_mismatches = max_mismatches
        self.seed_start = seed_start
        self.seed_end = seed_end

    def search(
        self,
        aso_sequence: str,
        on_target_transcript: str | None = None,
        on_target_position: int | None = None,
    ) -> list[OffTargetHit]:
        """Search for off-target binding sites for an ASO.

        The ASO binds to its target by reverse complementarity.
        This function searches the transcriptome for sequences
        similar to the reverse complement of the ASO.

        Parameters
        ----------
        aso_sequence : str
            ASO sequence (DNA alphabet, 5'→3').
        on_target_transcript : str or None
            Transcript ID of the intended target. Hits in this
            transcript at the on-target position are flagged as
            on-target (not off-target).
        on_target_position : int or None
            Position of the on-target site in the target transcript.

        Returns
        -------
        list[OffTargetHit]
            All hits (on-target and off-target), sorted by mismatches.
        """
        self.reference.ensure_index()

        aso_seq = aso_sequence.upper()
        target_site = reverse_complement(aso_seq)
        L = len(target_site)
        k = self.reference.kmer_size

        if L < k:
            raise OffTargetError(f"ASO length {L} is shorter than k-mer size {k}")

        # Extract k-mers from the target site and look up in the index
        candidate_positions: set[tuple[int, int]] = set()
        for i in range(L - k + 1):
            kmer = target_site[i : i + k]
            if kmer in self.reference._index:
                for tidx, pos in self.reference._index[kmer]:
                    # Adjust position to the start of the full-length match
                    start = pos - i
                    if start >= 0 and start + L <= self.reference.get_transcript(tidx).length:
                        candidate_positions.add((tidx, start))

        # Verify each candidate with full-length Hamming distance
        hits: list[OffTargetHit] = []
        seen: set[tuple[int, int]] = set()

        for tidx, start in candidate_positions:
            if (tidx, start) in seen:
                continue
            seen.add((tidx, start))

            transcript = self.reference.get_transcript(tidx)
            aligned = transcript.sequence[start : start + L]
            if len(aligned) < L:
                continue

            mismatches = _hamming_distance(target_site, aligned)
            if mismatches > self.max_mismatches:
                continue

            # Count seed mismatches
            # Seed in target site: positions [L - seed_end, L - seed_start)
            seed_start_target = L - self.seed_end
            seed_end_target = L - self.seed_start
            seed_target = target_site[seed_start_target:seed_end_target]
            seed_aligned = aligned[seed_start_target:seed_end_target]
            seed_mismatches = _hamming_distance(seed_target, seed_aligned)

            # Determine if this is the on-target hit
            is_on_target = (
                on_target_transcript is not None
                and transcript.transcript_id == on_target_transcript
                and on_target_position is not None
                and start == on_target_position
            )

            hits.append(
                OffTargetHit(
                    transcript_id=transcript.transcript_id,
                    gene_name=transcript.gene_name,
                    position=start,
                    mismatches=mismatches,
                    seed_mismatches=seed_mismatches,
                    aligned_sequence=aligned,
                    is_on_target=is_on_target,
                )
            )

        hits.sort(key=lambda h: (h.mismatches, h.seed_mismatches))
        return hits

    def search_batch(
        self,
        aso_sequences: list[str],
        on_target_transcripts: list[str | None] | None = None,
        on_target_positions: list[int | None] | None = None,
    ) -> dict[int, list[OffTargetHit]]:
        """Search for off-targets for multiple ASO candidates.

        Parameters
        ----------
        aso_sequences : list[str]
            List of ASO sequences.
        on_target_transcripts : list or None
            List of on-target transcript IDs (same length as aso_sequences).
        on_target_positions : list or None
            List of on-target positions (same length as aso_sequences).

        Returns
        -------
        dict[int, list[OffTargetHit]]
            Maps ASO index to its list of hits.
        """
        if on_target_transcripts is None:
            on_target_transcripts = [None] * len(aso_sequences)
        if on_target_positions is None:
            on_target_positions = [None] * len(aso_sequences)

        results: dict[int, list[OffTargetHit]] = {}
        for i, seq in enumerate(aso_sequences):
            results[i] = self.search(
                seq,
                on_target_transcript=on_target_transcripts[i],
                on_target_position=on_target_positions[i],
            )
        return results


# ---------------------------------------------------------------------------
# Result processing
# ---------------------------------------------------------------------------


def classify_hit(hit: OffTargetHit) -> str:
    """Classify an off-target hit by mismatch profile.

    Categories:
    - "perfect": 0 mismatches (exact match)
    - "seed_match": 0 seed mismatches, >0 total mismatches
    - "seed_mismatch": seed mismatches > 0
    - "on_target": the intended target site

    Parameters
    ----------
    hit : OffTargetHit
        The hit to classify.

    Returns
    -------
    str
        Classification label.
    """
    if hit.is_on_target:
        return "on_target"
    if hit.mismatches == 0:
        return "perfect"
    if hit.seed_mismatches == 0:
        return "seed_match"
    return "seed_mismatch"


def filter_hits(
    hits: list[OffTargetHit],
    max_mismatches: int = 3,
    max_seed_mismatches: int | None = None,
    exclude_on_target: bool = True,
) -> list[OffTargetHit]:
    """Filter off-target hits by mismatch criteria.

    Parameters
    ----------
    hits : list[OffTargetHit]
        Raw hits from the search engine.
    max_mismatches : int
        Maximum total mismatches (default 3).
    max_seed_mismatches : int or None
        Maximum seed mismatches (None = no seed filter).
    exclude_on_target : bool
        If True, exclude the on-target hit from results.

    Returns
    -------
    list[OffTargetHit]
        Filtered hits.
    """
    filtered = []
    for hit in hits:
        if exclude_on_target and hit.is_on_target:
            continue
        if hit.mismatches > max_mismatches:
            continue
        if max_seed_mismatches is not None and hit.seed_mismatches > max_seed_mismatches:
            continue
        filtered.append(hit)
    return filtered


@dataclass
class OffTargetResult:
    """Summary of off-target analysis for a single ASO candidate.

    Attributes
    ----------
    aso_sequence : str
        The ASO sequence.
    total_hits : int
        Total number of hits (including on-target).
    off_target_hits : int
        Number of off-target hits (excluding on-target).
    perfect_hits : int
        Off-target hits with 0 mismatches.
    seed_match_hits : int
        Off-target hits with 0 seed mismatches but >0 total mismatches.
    seed_mismatch_hits : int
        Off-target hits with seed mismatches.
    min_mismatches : int
        Minimum mismatches among off-target hits.
    hits : list[OffTargetHit]
        All off-target hits (excluding on-target).
    has_on_target : bool
        Whether the on-target site was found.
    """

    aso_sequence: str
    total_hits: int
    off_target_hits: int
    perfect_hits: int
    seed_match_hits: int
    seed_mismatch_hits: int
    min_mismatches: int
    hits: list[OffTargetHit]
    has_on_target: bool


def summarize_hits(
    aso_sequence: str,
    hits: list[OffTargetHit],
    exclude_on_target: bool = True,
) -> OffTargetResult:
    """Summarize off-target hits for an ASO candidate.

    Parameters
    ----------
    aso_sequence : str
        The ASO sequence.
    hits : list[OffTargetHit]
        All hits from the search engine.
    exclude_on_target : bool
        If True, exclude the on-target hit from the summary.

    Returns
    -------
    OffTargetResult
        Summary of the off-target profile.
    """
    off_hits = filter_hits(hits, exclude_on_target=exclude_on_target)
    on_target_hits = [h for h in hits if h.is_on_target]

    perfect = sum(1 for h in off_hits if h.mismatches == 0)
    seed_match = sum(1 for h in off_hits if h.mismatches > 0 and h.seed_mismatches == 0)
    seed_mismatch = sum(1 for h in off_hits if h.seed_mismatches > 0)

    min_mm = min((h.mismatches for h in off_hits), default=-1)

    return OffTargetResult(
        aso_sequence=aso_sequence,
        total_hits=len(hits),
        off_target_hits=len(off_hits),
        perfect_hits=perfect,
        seed_match_hits=seed_match,
        seed_mismatch_hits=seed_mismatch,
        min_mismatches=min_mm,
        hits=off_hits,
        has_on_target=len(on_target_hits) > 0,
    )


# ---------------------------------------------------------------------------
# Off-target scoring
# ---------------------------------------------------------------------------

# Penalty weights per hit category (relative to max_hits)
_PERFECT_PENALTY = 1.0
_SEED_MATCH_PENALTY = 0.5
_SEED_MISMATCH_PENALTY = 0.25


def calculate_off_target_score(
    result: OffTargetResult,
    seed_mismatch_weight: float = 2.0,
    max_hits: int = 50,
) -> float:
    """Compute an off-target specificity score in [0.0, 1.0].

    Score = 1.0 - (weighted_hits / max_hits), clamped to [0, 1].

    Weighted hits:
    - Perfect matches (0 mismatches): weight 1.0
    - Seed matches (0 seed mm, >0 total mm): weight 0.5 × seed_mismatch_weight
    - Seed mismatches: weight 0.25

    A score of 1.0 means no off-targets (highly specific).
    A score of 0.0 means many off-targets (non-specific).

    Parameters
    ----------
    result : OffTargetResult
        Summary of off-target hits for a candidate.
    seed_mismatch_weight : float
        Weight for seed-region matches (from config, default 2.0).
    max_hits : int
        Maximum hits threshold for normalization (from config, default 50).

    Returns
    -------
    float
        Specificity score in [0.0, 1.0].
    """
    weighted = (
        result.perfect_hits * _PERFECT_PENALTY
        + result.seed_match_hits * _SEED_MATCH_PENALTY * seed_mismatch_weight
        + result.seed_mismatch_hits * _SEED_MISMATCH_PENALTY
    )

    if max_hits <= 0:
        return 0.0

    score = 1.0 - (weighted / max_hits)
    return max(0.0, min(1.0, score))


def classify_off_target(
    score: float,
    specific_threshold: float = 0.8,
    moderate_threshold: float = 0.3,
) -> str:
    """Classify a candidate by off-target specificity.

    Parameters
    ----------
    score : float
        Off-target specificity score (0.0–1.0).
    specific_threshold : float
        Scores at or above this are "specific" (default 0.8).
    moderate_threshold : float
        Scores at or above this (but below specific_threshold) are
        "moderate" (default 0.3). Below this is "non-specific".

    Returns
    -------
    str
        One of "specific", "moderate", or "non-specific".
    """
    if score >= specific_threshold:
        return "specific"
    elif score >= moderate_threshold:
        return "moderate"
    else:
        return "non-specific"


def score_off_target(
    aso_sequence: str,
    hits: list[OffTargetHit],
    seed_mismatch_weight: float = 2.0,
    max_hits: int = 50,
    exclude_on_target: bool = True,
) -> tuple[OffTargetResult, float, str]:
    """Complete off-target analysis: summarize, score, and classify.

    Convenience function that chains summarize_hits → calculate score
    → classify.

    Parameters
    ----------
    aso_sequence : str
        The ASO sequence.
    hits : list[OffTargetHit]
        Raw hits from the search engine.
    seed_mismatch_weight : float
        Weight for seed-region matches (default 2.0).
    max_hits : int
        Max hits threshold (default 50).
    exclude_on_target : bool
        If True, exclude on-target hit from scoring.

    Returns
    -------
    tuple[OffTargetResult, float, str]
        (Summary, score, classification)
    """
    result = summarize_hits(aso_sequence, hits, exclude_on_target=exclude_on_target)
    score = calculate_off_target_score(
        result, seed_mismatch_weight=seed_mismatch_weight, max_hits=max_hits
    )
    classification = classify_off_target(score)
    return result, score, classification


def score_off_target_batch(
    aso_sequences: list[str],
    all_hits: dict[int, list[OffTargetHit]],
    seed_mismatch_weight: float = 2.0,
    max_hits: int = 50,
    exclude_on_target: bool = True,
) -> dict[int, tuple[OffTargetResult, float, str]]:
    """Score off-targets for multiple ASO candidates.

    Parameters
    ----------
    aso_sequences : list[str]
        List of ASO sequences.
    all_hits : dict[int, list[OffTargetHit]]
        Maps ASO index to its hits (from search_batch).
    seed_mismatch_weight : float
        Weight for seed-region matches.
    max_hits : int
        Max hits threshold.
    exclude_on_target : bool
        If True, exclude on-target hits.

    Returns
    -------
    dict[int, tuple[OffTargetResult, float, str]]
        Maps ASO index to (summary, score, classification).
    """
    results: dict[int, tuple[OffTargetResult, float, str]] = {}
    for i, seq in enumerate(aso_sequences):
        hits = all_hits.get(i, [])
        results[i] = score_off_target(
            seq,
            hits,
            seed_mismatch_weight=seed_mismatch_weight,
            max_hits=max_hits,
            exclude_on_target=exclude_on_target,
        )
    return results
