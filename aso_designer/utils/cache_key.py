"""Deterministic cache-key generation shared across caching modules.

Modules that persist computed results on disk (transcript cache,
ViennaRNA results, accessibility profiles) hash an identifier into a
fixed-length filename-safe key. A single helper keeps the digest
length and joining convention consistent.
"""

import hashlib


def cache_key(*parts: object, length: int = 16) -> str:
    """Build a short, deterministic, filename-safe cache key.

    Parameters
    ----------
    *parts : object
        Identifier components, joined with ``|`` before hashing.
    length : int
        Number of hex characters to keep (default 16 — 64 bits, enough
        to avoid collisions for repository-scale caches).

    Returns
    -------
    str
        Hex digest prefix of the SHA-256 hash of the joined parts.
    """
    joined = "|".join(str(part) for part in parts)
    return hashlib.sha256(joined.encode()).hexdigest()[:length]
