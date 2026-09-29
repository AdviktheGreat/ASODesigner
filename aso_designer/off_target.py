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
from aso_designer.utils.sequence import validate_sequence

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
        return len(self._transcripts)

    @property
    def total_length(self) -> int:
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
