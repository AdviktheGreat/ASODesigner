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


# ---------------------------------------------------------------------------
# Batch directory I/O
# ---------------------------------------------------------------------------


def save_candidates_batch(
    candidates: list[ASOCandidate],
    output_dir: str | Path,
    batch_size: int = 5000,
    metadata: dict[str, Any] | None = None,
) -> list[Path]:
    """Save candidates to a directory in batched JSON files.

    For large candidate sets, splitting into multiple files avoids
    creating a single oversized JSON file.

    Parameters
    ----------
    candidates : list[ASOCandidate]
        Candidates to save.
    output_dir : str or Path
        Output directory (created if it doesn't exist).
    batch_size : int
        Number of candidates per file (default 5000).
    metadata : dict, optional
        Metadata to include in each file.

    Returns
    -------
    list[Path]
        Paths to all written files.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    paths: list[Path] = []
    total_batches = (len(candidates) + batch_size - 1) // batch_size

    for i in range(0, len(candidates), batch_size):
        batch_num = i // batch_size
        batch = candidates[i : i + batch_size]
        filename = f"candidates_{batch_num:04d}_of_{total_batches:04d}.json"
        filepath = output_dir / filename

        batch_meta = {"batch_index": batch_num, "total_batches": total_batches}
        if metadata is not None:
            batch_meta.update(metadata)

        save_candidates(batch, filepath, metadata=batch_meta)
        paths.append(filepath)

    logger.info(
        "Saved %d candidates in %d batch files to %s", len(candidates), len(paths), output_dir
    )
    return paths


def load_candidates_batch(
    input_dir: str | Path,
) -> tuple[list[ASOCandidate], dict[str, Any]]:
    """Load all candidate JSON files from a directory.

    Parameters
    ----------
    input_dir : str or Path
        Directory containing candidate JSON files.

    Returns
    -------
    tuple[list[ASOCandidate], dict]
        (all_candidates, merged_metadata) — metadata from the first file.
    """
    input_dir = Path(input_dir)
    if not input_dir.is_dir():
        raise FileNotFoundError(f"Directory not found: {input_dir}")

    json_files = sorted(input_dir.glob("*.json"))
    if not json_files:
        raise FileNotFoundError(f"No JSON files found in {input_dir}")

    all_candidates: list[ASOCandidate] = []
    merged_metadata: dict[str, Any] = {}

    for fpath in json_files:
        candidates, meta = load_candidates(fpath)
        all_candidates.extend(candidates)
        if not merged_metadata and meta:
            merged_metadata = meta

    logger.info(
        "Loaded %d candidates from %d files in %s", len(all_candidates), len(json_files), input_dir
    )
    return all_candidates, merged_metadata


# ---------------------------------------------------------------------------
# Candidate comparison and deduplication
# ---------------------------------------------------------------------------


def candidates_equal(a: ASOCandidate, b: ASOCandidate) -> bool:
    """Check if two candidates are equal by sequence and position.

    Two candidates are considered equal if they have the same sequence
    and the same start position (and therefore the same end position).

    Parameters
    ----------
    a, b : ASOCandidate
        Candidates to compare.

    Returns
    -------
    bool
        True if the candidates are equal.
    """
    return a.sequence == b.sequence and a.start == b.start and a.end == b.end


def deduplicate_candidates(
    candidates: list[ASOCandidate],
    keep: str = "first",
) -> list[ASOCandidate]:
    """Remove duplicate candidates by sequence.

    Candidates with identical sequences are considered duplicates.
    The first or last occurrence is kept based on the ``keep`` parameter.

    Parameters
    ----------
    candidates : list[ASOCandidate]
        Candidates to deduplicate.
    keep : str
        Which duplicate to keep: "first" or "last" (default "first").

    Returns
    -------
    list[ASOCandidate]
        Deduplicated candidates, preserving order.
    """
    seen: dict[str, ASOCandidate] = {}
    if keep == "first":
        for c in candidates:
            if c.sequence not in seen:
                seen[c.sequence] = c
    elif keep == "last":
        for c in candidates:
            seen[c.sequence] = c
    else:
        raise ValueError(f"keep must be 'first' or 'last', got '{keep}'")

    result = list(seen.values())
    removed = len(candidates) - len(result)
    if removed > 0:
        logger.info("Deduplicated: removed %d duplicate candidates", removed)
    return result


def validate_roundtrip(candidates: list[ASOCandidate], tmp_path: str | Path) -> bool:
    """Verify that save/load preserves all candidate data.

    Parameters
    ----------
    candidates : list[ASOCandidate]
        Candidates to test.
    tmp_path : str or Path
        Temporary file path for the test.

    Returns
    -------
    bool
        True if all candidates survive the round-trip unchanged.
    """
    save_candidates(candidates, tmp_path)
    loaded, _ = load_candidates(tmp_path)

    if len(loaded) != len(candidates):
        return False

    for orig, load in zip(candidates, loaded, strict=False):
        if not candidates_equal(orig, load):
            return False
        if orig.chemistry != load.chemistry:
            return False
        if orig.gc_content != load.gc_content:
            return False
        if orig.region != load.region:
            return False

    return True
