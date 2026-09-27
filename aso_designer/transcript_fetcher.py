"""NCBI Entrez transcript fetcher for ASODesigner.

Redesigned from GuideX's genome_fetcher.py — uses NCBI E-utilities
directly via the Entrez REST API instead of the datasets CLI. Provides
gene symbol resolution, accession-based fetch, ortholog retrieval, and
local file fallback with caching.

Requires NCBI_API_KEY environment variable for higher rate limits
(10 req/s with key, 3 req/s without).
"""

import hashlib
import os
import time
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from pathlib import Path
from typing import List, Optional, Union
from urllib.parse import urlencode

import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from aso_designer.utils.exceptions import TranscriptFetchError
from aso_designer.utils.logger import get_logger
from aso_designer.utils.sequence import validate_sequence

logger = get_logger(__name__)

# NCBI E-utilities base URL
_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Default cache location
_DEFAULT_CACHE_DIR = Path.home() / ".asodesigner" / "transcript_cache"

# Rate limit: 10 req/s with API key, 3 req/s without
_RATE_WITH_KEY = 10.0
_RATE_WITHOUT_KEY = 3.0


@dataclass
class TranscriptRecord:
    """A fetched transcript with metadata."""

    accession: str
    sequence: str
    gene_name: str
    organism: str
    chromosome: Optional[str] = None
    strand: Optional[str] = None
    description: str = ""
    length: int = 0

    def __post_init__(self) -> None:
        self.length = len(self.sequence)

    def to_seq_record(self) -> SeqRecord:
        """Convert to a BioPython SeqRecord."""
        return SeqRecord(
            Seq(self.sequence),
            id=self.accession,
            description=f"{self.gene_name} | {self.organism} | {self.description}",
        )

    def to_dict(self) -> dict:
        """Serialize to dictionary for JSON storage."""
        return {
            "accession": self.accession,
            "sequence": self.sequence,
            "gene_name": self.gene_name,
            "organism": self.organism,
            "chromosome": self.chromosome,
            "strand": self.strand,
            "description": self.description,
            "length": self.length,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "TranscriptRecord":
        """Deserialize from dictionary."""
        return cls(**data)


class EntrezFetcher:
    """NCBI E-utilities client for transcript sequence fetching.

    Features:
    - Rate limiting with automatic API key detection
    - Retry with exponential backoff
    - Local caching with SHA-based keys and TTL
    - Gene symbol → RefSeq transcript resolution
    - Accession-based direct fetch
    - Ortholog retrieval for conservation analysis

    Parameters
    ----------
    api_key : str, optional
        NCBI API key. Falls back to NCBI_API_KEY env var.
    email : str, optional
        Email for NCBI identification (recommended by NCBI).
    cache_dir : Path, optional
        Directory for caching fetched transcripts.
    cache_ttl : int
        Cache time-to-live in hours (default 168 = 7 days).
    """

    def __init__(
        self,
        api_key: Optional[str] = None,
        email: Optional[str] = None,
        cache_dir: Optional[Union[str, Path]] = None,
        cache_ttl: int = 168,
    ) -> None:
        self.api_key = api_key or os.getenv("NCBI_API_KEY")
        self.email = email or os.getenv("NCBI_EMAIL", "asodesigner@local")
        self.cache_dir = Path(cache_dir) if cache_dir else _DEFAULT_CACHE_DIR
        self.cache_ttl = timedelta(hours=cache_ttl)
        self.rate_limit = _RATE_WITH_KEY if self.api_key else _RATE_WITHOUT_KEY
        self._min_interval = 1.0 / self.rate_limit
        self._last_request_time: float = 0.0

        self.session = requests.Session()
        self.session.headers.update({"User-Agent": "ASODesigner/0.1"})
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        logger.info(
            "EntrezFetcher initialized (rate: %.1f req/s, cache: %s)",
            self.rate_limit,
            self.cache_dir,
        )

    def _rate_limit_wait(self) -> None:
        """Enforce minimum interval between requests."""
        elapsed = time.time() - self._last_request_time
        if elapsed < self._min_interval:
            time.sleep(self._min_interval - elapsed)
        self._last_request_time = time.time()

    def _build_params(self, **kwargs) -> dict:
        """Build E-utilities query params with API key and email."""
        params = {"tool": "asodesigner", "email": self.email}
        if self.api_key:
            params["api_key"] = self.api_key
        params.update(kwargs)
        return params

    def _get(self, endpoint: str, **kwargs) -> dict:
        """Make a rate-limited GET request to E-utilities with retry.

        Parameters
        ----------
        endpoint : str
            E-utilities endpoint (e.g. "esearch.fcgi", "efetch.fcgi").
        **kwargs
            Query parameters.

        Returns
        -------
        dict
            Parsed JSON response (for JSON mode) or raw text content.

        Raises
        ------
        TranscriptFetchError
            If the request fails after all retries.
        """
        url = f"{_EUTILS_BASE}/{endpoint}"
        params = self._build_params(**kwargs)

        max_retries = 3
        for attempt in range(max_retries):
            self._rate_limit_wait()
            try:
                response = self.session.get(url, params=params, timeout=30)
                response.raise_for_status()

                # Try JSON parse, fall back to text
                content_type = response.headers.get("Content-Type", "")
                if "json" in content_type or "retmode=json" in str(kwargs):
                    return response.json()
                return {"text": response.text}

            except requests.RequestException as e:
                wait = 2 ** attempt
                logger.warning(
                    "Request failed (attempt %d/%d): %s — retrying in %ds",
                    attempt + 1,
                    max_retries,
                    e,
                    wait,
                )
                if attempt < max_retries - 1:
                    time.sleep(wait)
                else:
                    raise TranscriptFetchError(
                        f"NCBI request failed after {max_retries} retries: {e}"
                    ) from e

        # Should not reach here
        raise TranscriptFetchError("Unexpected: exhausted retries without raising")

    def _cache_key(self, identifier: str) -> Path:
        """Generate a cache file path for a given identifier."""
        key_hash = hashlib.sha256(identifier.encode()).hexdigest()[:16]
        return self.cache_dir / f"{key_hash}.json"

    def _is_cache_valid(self, cache_path: Path) -> bool:
        """Check if a cache file exists and is within TTL."""
        if not cache_path.exists():
            return False
        mtime = datetime.fromtimestamp(cache_path.stat().st_mtime)
        return datetime.now() - mtime < self.cache_ttl

    def _save_cache(self, identifier: str, records: List[TranscriptRecord]) -> None:
        """Save transcript records to cache."""
        import json

        cache_path = self._cache_key(identifier)
        data = [r.to_dict() for r in records]
        cache_path.write_text(json.dumps(data, indent=2), encoding="utf-8")
        logger.debug("Cached %d records for '%s' at %s", len(records), identifier, cache_path)

    def _load_cache(self, identifier: str) -> Optional[List[TranscriptRecord]]:
        """Load transcript records from cache if valid."""
        import json

        cache_path = self._cache_key(identifier)
        if not self._is_cache_valid(cache_path):
            return None
        try:
            data = json.loads(cache_path.read_text(encoding="utf-8"))
            records = [TranscriptRecord.from_dict(d) for d in data]
            logger.debug("Loaded %d cached records for '%s'", len(records), identifier)
            return records
        except (json.JSONDecodeError, KeyError) as e:
            logger.warning("Cache corruption for '%s': %s", identifier, e)
            return None

    # ------------------------------------------------------------------
    # Gene symbol → RefSeq transcript resolution
    # ------------------------------------------------------------------

    def resolve_gene_to_transcripts(
        self,
        gene_symbol: str,
        organism: str = "Homo sapiens",
    ) -> List[str]:
        """Resolve a gene symbol to RefSeq transcript accession IDs.

        Uses ESearch (nucleotide db) with a gene+organism query to find
        RefSeq mRNA transcripts (NM_ accessions).

        Parameters
        ----------
        gene_symbol : str
            Gene symbol (e.g. "DMD", "SMN1").
        organism : str
            Organism name (default "Homo sapiens").

        Returns
        -------
        List[str]
            List of RefSeq mRNA accession IDs (NM_ or NR_ prefix).

        Raises
        ------
        TranscriptFetchError
            If the search fails or returns no results.
        """
        cache_id = f"gene:{gene_symbol}:{organism}"
        cached = self._load_cache(cache_id)
        if cached is not None:
            return [r.accession for r in cached]

        # ESearch for RefSeq transcripts
        query = f"{gene_symbol}[Gene] AND {organism}[Organism] AND refseq[filter] AND biomol_mRNA[Properties]"
        result = self._get(
            "esearch.fcgi",
            db="nucleotide",
            term=query,
            retmax=50,
            retmode="json",
            sort="length_desc",  # Longest transcripts first
        )

        id_list = result.get("esearchresult", {}).get("idlist", [])
        if not id_list:
            logger.warning("No RefSeq transcripts found for gene '%s' in %s", gene_symbol, organism)
            raise TranscriptFetchError(
                f"No RefSeq transcripts found for gene '{gene_symbol}' in {organism}"
            )

        # ESummary to get accession IDs from GI numbers
        summaries = self._get(
            "esummary.fcgi",
            db="nucleotide",
            id=",".join(id_list),
            retmode="json",
        )

        accessions: List[str] = []
        for gi in id_list:
            entry = summaries.get("result", {}).get(gi, {})
            acc = entry.get("accessionversion", "")
            title = entry.get("title", "")
            # Only keep mRNA transcripts (NM_ prefix)
            if acc.startswith("NM_"):
                accessions.append(acc)
            elif "mRNA" in title and acc.startswith(("NM_", "XM_")):
                accessions.append(acc)

        if not accessions:
            raise TranscriptFetchError(
                f"No mRNA accessions resolved for gene '{gene_symbol}'"
            )

        logger.info("Resolved gene '%s' to %d transcript(s): %s", gene_symbol, len(accessions), accessions[:5])
        return accessions
