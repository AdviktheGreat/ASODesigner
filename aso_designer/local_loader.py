"""Local FASTA file loader for ASODesigner.

Loads transcript sequences from local FASTA files with validation.
Designed from scratch for ASO-specific needs — supports multi-record
FASTA, header parsing for gene/organism metadata, and integration
with the TranscriptRecord dataclass.
"""

from pathlib import Path
from typing import List, Optional, Union

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from aso_designer.transcript_fetcher import TranscriptRecord
from aso_designer.utils.exceptions import InvalidSequenceError, TranscriptFetchError
from aso_designer.utils.logger import get_logger
from aso_designer.utils.sequence import validate_sequence

logger = get_logger(__name__)


def load_fasta(
    file_path: Union[str, Path],
    min_length: int = 100,
    allow_iupac: bool = True,
) -> List[TranscriptRecord]:
    """Load transcript sequences from a local FASTA file.

    Parameters
    ----------
    file_path : str or Path
        Path to the FASTA file.
    min_length : int
        Minimum sequence length to accept (default 100 nt).
    allow_iupac : bool
        Whether to accept IUPAC ambiguity codes.

    Returns
    -------
    List[TranscriptRecord]
        Parsed transcript records with metadata extracted from headers.

    Raises
    ------
    TranscriptFetchError
        If the file cannot be read or contains no valid sequences.
    """
    path = Path(file_path)
    if not path.exists():
        raise TranscriptFetchError(f"FASTA file not found: {path}")

    records: List[TranscriptRecord] = []
    try:
        for seq_record in SeqIO.parse(str(path), "fasta"):
            sequence = str(seq_record.seq).upper()

            # Skip sequences below minimum length
            if len(sequence) < min_length:
                logger.debug("Skipping %s: too short (%d nt < %d)", seq_record.id, len(sequence), min_length)
                continue

            # Validate sequence content
            try:
                sequence = validate_sequence(sequence, allow_rna=False, allow_iupac=allow_iupac)
            except InvalidSequenceError as e:
                logger.warning("Skipping %s: invalid sequence — %s", seq_record.id, e)
                continue

            # Parse metadata from FASTA header
            accession = seq_record.id
            description = seq_record.description
            gene_name = _parse_gene_from_header(description, accession)
            organism = _parse_organism_from_header(description)

            record = TranscriptRecord(
                accession=accession,
                sequence=sequence,
                gene_name=gene_name,
                organism=organism,
                description=description,
            )
            records.append(record)

    except Exception as e:
        raise TranscriptFetchError(f"Failed to parse FASTA file {path}: {e}") from e

    if not records:
        raise TranscriptFetchError(f"No valid sequences found in {path}")

    logger.info("Loaded %d transcript(s) from %s", len(records), path)
    return records


def load_single_fasta(
    file_path: Union[str, Path],
    min_length: int = 100,
) -> TranscriptRecord:
    """Load a single transcript from a FASTA file.

    Convenience wrapper that returns the first valid record.

    Parameters
    ----------
    file_path : str or Path
        Path to the FASTA file.
    min_length : int
        Minimum sequence length.

    Returns
    -------
    TranscriptRecord
        First valid transcript record.

    Raises
    ------
    TranscriptFetchError
        If the file is empty or contains no valid sequences.
    """
    records = load_fasta(file_path, min_length=min_length)
    return records[0]


def _parse_gene_from_header(description: str, accession: str) -> str:
    """Extract gene name from a FASTA header (heuristic)."""
    import re

    # Try to find a gene symbol in parentheses: "dystrophin (DMD)"
    match = re.search(r"\(([A-Z][A-Z0-9]{1,10})\)", description)
    if match:
        return match.group(1)
    # Fall back to accession prefix
    return accession.split(".")[0] if "." in accession else accession


def _parse_organism_from_header(description: str) -> str:
    """Extract organism name from a FASTA header (heuristic)."""
    _common_organisms = (
        "Homo sapiens", "Mus musculus", "Rattus norvegicus",
        "Danio rerio", "Drosophila melanogaster", "Caenorhabditis elegans",
        "Saccharomyces cerevisiae", "Arabidopsis thaliana",
    )
    for org in _common_organisms:
        if org in description:
            return org
    return "unknown"
