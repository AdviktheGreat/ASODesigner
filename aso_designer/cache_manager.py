"""Transcript cache manager for ASODesigner.

Extends the basic caching in EntrezFetcher with organism-specific
directory organization, inspection utilities, and cache size management.
Designed from scratch — GuideX uses a flat cache directory; we organize
by organism for faster bulk operations and easier manual inspection.
"""

import hashlib
import json
from datetime import datetime, timedelta
from pathlib import Path

from aso_designer.transcript_fetcher import TranscriptRecord
from aso_designer.utils.logger import get_logger

logger = get_logger(__name__)

_DEFAULT_CACHE_BASE = Path.home() / ".asodesigner" / "transcript_cache"


class TranscriptCacheManager:
    """Organism-organized transcript cache with TTL and inspection utilities.

    Parameters
    ----------
    cache_base : Path, optional
        Base directory for the cache. Organism subdirectories are created
        beneath this.
    default_ttl : int
        Default cache time-to-live in hours (default 168 = 7 days).
    """

    def __init__(
        self,
        cache_base: str | Path | None = None,
        default_ttl: int = 168,
    ) -> None:
        self.cache_base = Path(cache_base) if cache_base else _DEFAULT_CACHE_BASE
        self.default_ttl = timedelta(hours=default_ttl)
        self.cache_base.mkdir(parents=True, exist_ok=True)

    def _organism_dir(self, organism: str) -> Path:
        """Get or create the cache directory for an organism."""
        safe_name = organism.replace(" ", "_").replace("/", "_")
        org_dir = self.cache_base / safe_name
        org_dir.mkdir(parents=True, exist_ok=True)
        return org_dir

    def _cache_key_path(self, identifier: str, organism: str = "unknown") -> Path:
        """Generate a cache file path organized by organism."""
        org_dir = self._organism_dir(organism)
        key_hash = hashlib.sha256(identifier.encode()).hexdigest()[:16]
        return org_dir / f"{key_hash}.json"

    def store(
        self,
        identifier: str,
        records: list[TranscriptRecord],
        organism: str = "unknown",
    ) -> None:
        """Store transcript records in the cache.

        Parameters
        ----------
        identifier : str
            Cache key (accession, gene symbol, etc.).
        records : List[TranscriptRecord]
            Records to cache.
        organism : str
            Organism name for directory organization.
        """
        cache_path = self._cache_key_path(identifier, organism)
        data = [r.to_dict() for r in records]
        cache_path.write_text(json.dumps(data, indent=2), encoding="utf-8")
        logger.debug("Cached %d records for '%s' (%s)", len(records), identifier, organism)

    def retrieve(
        self,
        identifier: str,
        organism: str = "unknown",
        ttl: timedelta | None = None,
    ) -> list[TranscriptRecord] | None:
        """Retrieve cached records if they exist and are within TTL.

        Parameters
        ----------
        identifier : str
            Cache key.
        organism : str
            Organism name for directory lookup.
        ttl : timedelta, optional
            Override default TTL. If None, uses default_ttl.

        Returns
        -------
        List[TranscriptRecord] or None
            Cached records if valid, None otherwise.
        """
        cache_path = self._cache_key_path(identifier, organism)
        if not cache_path.exists():
            return None

        # Check TTL
        effective_ttl = ttl or self.default_ttl
        mtime = datetime.fromtimestamp(cache_path.stat().st_mtime)
        if datetime.now() - mtime > effective_ttl:
            logger.debug("Cache expired for '%s' (age: %s)", identifier, datetime.now() - mtime)
            return None

        try:
            data = json.loads(cache_path.read_text(encoding="utf-8"))
            return [TranscriptRecord.from_dict(d) for d in data]
        except (json.JSONDecodeError, KeyError) as e:
            logger.warning("Cache corruption for '%s': %s", identifier, e)
            return None

    def invalidate(self, identifier: str, organism: str = "unknown") -> bool:
        """Delete a specific cache entry.

        Returns
        -------
        bool
            True if the entry was deleted, False if it didn't exist.
        """
        cache_path = self._cache_key_path(identifier, organism)
        if cache_path.exists():
            cache_path.unlink()
            logger.info("Invalidated cache for '%s'", identifier)
            return True
        return False

    # ------------------------------------------------------------------
    # Inspection and management utilities
    # ------------------------------------------------------------------

    def inspect(self) -> dict[str, dict]:
        """Inspect the cache and return a summary per organism.

        Returns
        -------
        dict
            Mapping of organism → {files, total_size_bytes, oldest_entry, newest_entry}.
        """
        report: dict[str, dict] = {}
        for org_dir in sorted(self.cache_base.iterdir()):
            if not org_dir.is_dir():
                continue
            organism = org_dir.name.replace("_", " ")
            files = list(org_dir.glob("*.json"))
            if not files:
                continue
            total_size = sum(f.stat().st_size for f in files)
            mtimes = [datetime.fromtimestamp(f.stat().st_mtime) for f in files]
            report[organism] = {
                "files": len(files),
                "total_size_bytes": total_size,
                "total_size_mb": round(total_size / (1024 * 1024), 2),
                "oldest_entry": min(mtimes).isoformat(),
                "newest_entry": max(mtimes).isoformat(),
            }
        return report

    def clear_organism(self, organism: str) -> int:
        """Clear all cached entries for a specific organism.

        Parameters
        ----------
        organism : str
            Organism name (will be normalized to directory name).

        Returns
        -------
        int
            Number of files deleted.
        """
        org_dir = self._organism_dir(organism)
        count = 0
        for f in org_dir.glob("*.json"):
            f.unlink()
            count += 1
        logger.info("Cleared %d cache entries for %s", count, organism)
        return count

    def clear_all(self) -> int:
        """Clear the entire cache.

        Returns
        -------
        int
            Number of files deleted.
        """
        count = 0
        for org_dir in self.cache_base.iterdir():
            if not org_dir.is_dir():
                continue
            for f in org_dir.glob("*.json"):
                f.unlink()
                count += 1
        logger.info("Cleared entire cache (%d files)", count)
        return count

    def total_size_bytes(self) -> int:
        """Return total cache size in bytes."""
        total = 0
        for org_dir in self.cache_base.iterdir():
            if not org_dir.is_dir():
                continue
            total += sum(f.stat().st_size for f in org_dir.glob("*.json"))
        return total

    def total_entries(self) -> int:
        """Return total number of cached entries."""
        count = 0
        for org_dir in self.cache_base.iterdir():
            if not org_dir.is_dir():
                continue
            count += len(list(org_dir.glob("*.json")))
        return count
