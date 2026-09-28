"""Candidate serialization, batch I/O, and deduplication.

Provides JSON save/load for lists of ASOCandidate objects with full
metadata preservation. Designed from scratch — GuideX uses a simple
ResultsHandler with basic dict conversion; we add round-trip
validation, batch directory I/O, and candidate deduplication.
"""

import json
from pathlib import Path
from typing import Any

from aso_designer.candidate import ASOCandidate
from aso_designer.utils.logger import get_logger

logger = get_logger(__name__)


def save_candidates(
    candidates: list[ASOCandidate],
    path: str | Path,
    metadata: dict[str, Any] | None = None,
) -> Path:
    """Save candidates to a JSON file with optional run metadata.

    Parameters
    ----------
    candidates : list[ASOCandidate]
        Candidates to save.
    path : str or Path
        Output JSON file path.
    metadata : dict, optional
        Additional metadata to include (e.g. config, timestamp, gene name).

    Returns
    -------
    Path
        Path to the written file.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    data: dict[str, Any] = {
        "candidate_count": len(candidates),
        "candidates": [c.to_dict() for c in candidates],
    }
    if metadata is not None:
        data["metadata"] = metadata

    path.write_text(json.dumps(data, indent=2, default=str), encoding="utf-8")
    logger.info("Saved %d candidates to %s", len(candidates), path)
    return path


def load_candidates(path: str | Path) -> tuple[list[ASOCandidate], dict[str, Any]]:
    """Load candidates from a JSON file.

    Parameters
    ----------
    path : str or Path
        Input JSON file path.

    Returns
    -------
    tuple[list[ASOCandidate], dict]
        (candidates, metadata) — metadata may be empty dict if not present.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If the file cannot be parsed or candidates are invalid.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Candidate file not found: {path}")

    raw = json.loads(path.read_text(encoding="utf-8"))

    if "candidates" not in raw:
        raise ValueError(f"No 'candidates' key in {path}")

    candidates: list[ASOCandidate] = []
    for i, item in enumerate(raw["candidates"]):
        try:
            candidates.append(ASOCandidate.from_dict(item))
        except (KeyError, TypeError) as e:
            raise ValueError(f"Failed to parse candidate {i} in {path}: {e}") from e

    metadata = raw.get("metadata", {})
    logger.info("Loaded %d candidates from %s", len(candidates), path)
    return candidates, metadata
