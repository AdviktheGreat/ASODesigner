"""Unified transcript fetch manager with fallback chain.

Ties together the NCBI EntrezFetcher, TranscriptCacheManager, and local
FASTA loader into a single interface that tries sources in order:
  1. NCBI API (if network available)
  2. Local cache (if not expired)
  3. Local FASTA file (if provided)

Each step is logged so the user can trace which source was used.
"""

from pathlib import Path

from aso_designer.cache_manager import TranscriptCacheManager
from aso_designer.local_loader import load_fasta
from aso_designer.transcript_fetcher import EntrezFetcher, TranscriptRecord
from aso_designer.utils.exceptions import TranscriptFetchError
from aso_designer.utils.logger import get_logger

logger = get_logger(__name__)


class FetchManager:
    """Unified transcript fetcher with NCBI → cache → local file fallback.

    Parameters
    ----------
    api_key : str, optional
        NCBI API key for higher rate limits.
    email : str, optional
        Email for NCBI identification.
    cache_dir : Path, optional
        Base directory for the transcript cache.
    cache_ttl : int
        Cache TTL in hours (default 168 = 7 days).
    local_fasta : Path, optional
        Path to a local FASTA file to use as a fallback source.
    """

    def __init__(
        self,
        api_key: str | None = None,
        email: str | None = None,
        cache_dir: str | Path | None = None,
        cache_ttl: int = 168,
        local_fasta: str | Path | None = None,
    ) -> None:
        self.cache = TranscriptCacheManager(cache_base=cache_dir, default_ttl=cache_ttl)
        self.local_fasta = Path(local_fasta) if local_fasta else None
        self._entrez: EntrezFetcher | None = None
        self._api_key = api_key
        self._email = email

    @property
    def entrez(self) -> EntrezFetcher:
        """Lazily initialize the EntrezFetcher (only when network is needed)."""
        if self._entrez is None:
            self._entrez = EntrezFetcher(
                api_key=self._api_key,
                email=self._email,
                cache_dir=self.cache.cache_base,
            )
        return self._entrez

    def fetch_by_gene(
        self,
        gene_symbol: str,
        organism: str = "Homo sapiens",
    ) -> list[TranscriptRecord]:
        """Fetch transcripts by gene symbol with fallback chain.

        Tries in order:
        1. Local cache (if valid)
        2. NCBI API (if network available)
        3. Local FASTA file (if provided)

        Parameters
        ----------
        gene_symbol : str
            Gene symbol (e.g. "DMD").
        organism : str
            Organism name.

        Returns
        -------
        List[TranscriptRecord]
            Fetched transcript records.

        Raises
        ------
        TranscriptFetchError
            If all sources fail.
        """
        cache_id = f"gene:{gene_symbol}:{organism}"

        # 1. Try cache
        cached = self.cache.retrieve(cache_id, organism)
        if cached is not None:
            logger.info(
                "Cache hit for gene '%s' (%s) — %d transcript(s)",
                gene_symbol,
                organism,
                len(cached),
            )
            return cached

        # 2. Try NCBI API
        try:
            records = self.entrez.fetch_gene(gene_symbol, organism)
            if records:
                self.cache.store(cache_id, records, organism)
                logger.info(
                    "NCBI fetch successful for gene '%s' — %d transcript(s)",
                    gene_symbol,
                    len(records),
                )
                return records
        except TranscriptFetchError as e:
            logger.warning("NCBI fetch failed for gene '%s': %s", gene_symbol, e)

        # 3. Try local FASTA
        if self.local_fasta and self.local_fasta.exists():
            logger.info("Falling back to local FASTA: %s", self.local_fasta)
            records = load_fasta(self.local_fasta)
            self.cache.store(cache_id, records, organism)
            return records

        raise TranscriptFetchError(
            f"All sources failed for gene '{gene_symbol}' in {organism}. "
            "Provide a local FASTA file or check network connectivity."
        )

    def fetch_by_accession(
        self,
        accession: str,
        organism: str = "unknown",
    ) -> TranscriptRecord:
        """Fetch a single transcript by accession with fallback chain.

        Parameters
        ----------
        accession : str
            RefSeq accession (e.g. "NM_000109.3").
        organism : str
            Organism name for cache organization.

        Returns
        -------
        TranscriptRecord
            Fetched transcript record.

        Raises
        ------
        TranscriptFetchError
            If all sources fail.
        """
        cache_id = f"acc:{accession}"

        # 1. Try cache
        cached = self.cache.retrieve(cache_id, organism)
        if cached is not None and len(cached) > 0:
            logger.info("Cache hit for accession '%s'", accession)
            return cached[0]

        # 2. Try NCBI API
        try:
            record = self.entrez.fetch_by_accession(accession)
            self.cache.store(cache_id, [record], organism)
            logger.info("NCBI fetch successful for accession '%s'", accession)
            return record
        except TranscriptFetchError as e:
            logger.warning("NCBI fetch failed for accession '%s': %s", accession, e)

        # 3. Try local FASTA
        if self.local_fasta and self.local_fasta.exists():
            logger.info("Falling back to local FASTA: %s", self.local_fasta)
            records = load_fasta(self.local_fasta)
            for r in records:
                if r.accession == accession or r.accession.startswith(accession.split(".")[0]):
                    self.cache.store(cache_id, [r], organism)
                    return r

        raise TranscriptFetchError(
            f"All sources failed for accession '{accession}'. "
            "Provide a local FASTA file or check network connectivity."
        )

    def fetch_from_local(self, file_path: str | Path) -> list[TranscriptRecord]:
        """Fetch transcripts directly from a local FASTA file.

        Parameters
        ----------
        file_path : str or Path
            Path to the FASTA file.

        Returns
        -------
        List[TranscriptRecord]
            Parsed transcript records.
        """
        return load_fasta(file_path)
