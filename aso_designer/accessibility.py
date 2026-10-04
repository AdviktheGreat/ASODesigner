"""Target RNA structure prediction and accessibility analysis.

Uses ViennaRNA wrappers (RNAfold, RNAplfold) to compute:
- Full transcript MFE secondary structure
- Local accessibility profile (per-position unpaired probabilities)
- Normalized accessibility scores for ASO binding sites

Results are cached to disk to avoid expensive re-computation.
"""

from __future__ import annotations

import hashlib
import json
import pickle
from dataclasses import dataclass
from pathlib import Path

from aso_designer.utils.logger import get_logger
from aso_designer.utils.sequence import dna_to_rna
from aso_designer.viennarna import RNAfoldResult, RNAplfoldResult, run_rnafold, run_rnaplfold

logger = get_logger(__name__)


@dataclass
class TargetStructure:
    """Predicted secondary structure of a target transcript.

    Attributes
    ----------
    sequence : str
        RNA sequence (U-based).
    structure : str
        Dot-bracket MFE structure from RNAfold.
    mfe : float
        Minimum free energy in kcal/mol.
    """

    sequence: str
    structure: str
    mfe: float

    @classmethod
    def from_rnafold_result(cls, result: RNAfoldResult) -> TargetStructure:
        """Build a TargetStructure from an RNAfoldResult."""
        return cls(
            sequence=result.sequence,
            structure=result.structure,
            mfe=result.mfe,
        )


@dataclass
class AccessibilityProfile:
    """Per-position accessibility profile from RNAplfold.

    Attributes
    ----------
    sequence : str
        RNA sequence (U-based).
    positions : list[int]
        1-based positions in the sequence.
    unpaired_probs : list[float]
        Probability of being unpaired at each position (0.0–1.0).
    window_size : int
        RNAplfold sliding window size used.
    """

    sequence: str
    positions: list[int]
    unpaired_probs: list[float]
    window_size: int

    @classmethod
    def from_rnaplfold_result(
        cls, result: RNAplfoldResult, rna_seq: str, window_size: int
    ) -> AccessibilityProfile:
        """Build an AccessibilityProfile from an RNAplfoldResult."""
        return cls(
            sequence=rna_seq,
            positions=result.positions,
            unpaired_probs=result.unpaired_probs,
            window_size=window_size,
        )

    def get_unpaired_prob(self, position: int) -> float:
        """Get the unpaired probability at a 1-based position.

        Parameters
        ----------
        position : int
            1-based position in the sequence.

        Returns
        -------
        float
            Unpaired probability (0.0–1.0), or 0.0 if out of range.
        """
        if position < 1 or position > len(self.unpaired_probs):
            return 0.0
        return self.unpaired_probs[position - 1]


class TargetStructurePredictor:
    """Predicts target RNA structure and accessibility.

    Wraps RNAfold and RNAplfold with caching to avoid
    re-computation on repeated runs.

    Parameters
    ----------
    cache_dir : Path or str, optional
        Directory for caching structure predictions.
        If None, no caching is performed.
    window_size : int
        RNAplfold sliding window size (default 150).
    max_pairing_distance : int
        Maximum base pair distance for RNAplfold (default 100).
    """

    def __init__(
        self,
        cache_dir: Path | str | None = None,
        window_size: int = 150,
        max_pairing_distance: int = 100,
    ):
        self.cache_dir = Path(cache_dir) if cache_dir else None
        self.window_size = window_size
        self.max_pairing_distance = max_pairing_distance

        if self.cache_dir:
            self.cache_dir.mkdir(parents=True, exist_ok=True)

    def predict_structure(self, sequence: str) -> TargetStructure:
        """Predict the full MFE secondary structure of a sequence.

        Uses RNAfold and caches the result.

        Parameters
        ----------
        sequence : str
            DNA or RNA sequence (case-insensitive).

        Returns
        -------
        TargetStructure
            MFE structure and energy.
        """
        rna_seq = dna_to_rna(sequence)

        cache_key = self._cache_key(rna_seq, "structure")
        cached = self._load_cache(cache_key)
        if cached is not None:
            logger.debug("Using cached structure for sequence of length %d", len(rna_seq))
            return cached

        logger.info("Running RNAfold on sequence of length %d", len(rna_seq))
        result = run_rnafold(rna_seq)
        structure = TargetStructure.from_rnafold_result(result)

        self._save_cache(cache_key, structure)
        return structure

    def predict_accessibility(self, sequence: str) -> AccessibilityProfile:
        """Predict the local accessibility profile of a sequence.

        Uses RNAplfold and caches the result.

        Parameters
        ----------
        sequence : str
            DNA or RNA sequence (case-insensitive).

        Returns
        -------
        AccessibilityProfile
            Per-position unpaired probabilities.
        """
        rna_seq = dna_to_rna(sequence)

        cache_key = self._cache_key(
            rna_seq, f"accessibility_W{self.window_size}_L{self.max_pairing_distance}"
        )
        cached = self._load_cache(cache_key)
        if cached is not None:
            logger.debug("Using cached accessibility for sequence of length %d", len(rna_seq))
            return cached

        logger.info("Running RNAplfold on sequence of length %d", len(rna_seq))
        result = run_rnaplfold(
            rna_seq,
            window_size=self.window_size,
            max_span=self.max_pairing_distance,
        )
        profile = AccessibilityProfile.from_rnaplfold_result(result, rna_seq, self.window_size)

        self._save_cache(cache_key, profile)
        return profile

    def predict_all(self, sequence: str) -> tuple[TargetStructure, AccessibilityProfile]:
        """Predict both MFE structure and accessibility profile.

        Parameters
        ----------
        sequence : str
            DNA or RNA sequence.

        Returns
        -------
        tuple[TargetStructure, AccessibilityProfile]
            (MFE structure, accessibility profile)
        """
        structure = self.predict_structure(sequence)
        accessibility = self.predict_accessibility(sequence)
        return structure, accessibility

    # ------------------------------------------------------------------
    # Caching
    # ------------------------------------------------------------------

    def _cache_key(self, rna_seq: str, suffix: str) -> str:
        """Generate a cache key from sequence and suffix."""
        seq_hash = hashlib.md5(rna_seq.encode()).hexdigest()
        return f"{seq_hash}_{suffix}"

    def _cache_path(self, key: str) -> Path:
        return self.cache_dir / f"{key}.pkl"  # type: ignore[union-attr]

    def _load_cache(self, key: str) -> TargetStructure | AccessibilityProfile | None:
        """Load a cached result if available."""
        if self.cache_dir is None:
            return None
        path = self._cache_path(key)
        if not path.exists():
            return None
        try:
            with open(path, "rb") as f:
                return pickle.load(f)  # noqa: S301
        except Exception as e:
            logger.warning("Failed to load cache %s: %s", path, e)
            return None

    def _save_cache(self, key: str, result: TargetStructure | AccessibilityProfile) -> None:
        """Save a result to cache."""
        if self.cache_dir is None:
            return
        path = self._cache_path(key)
        try:
            with open(path, "wb") as f:
                pickle.dump(result, f)
        except Exception as e:
            logger.warning("Failed to save cache %s: %s", path, e)


# ---------------------------------------------------------------------------
# Accessibility scoring
# ---------------------------------------------------------------------------


def score_position_accessibility(unpaired_prob: float) -> float:
    """Convert an unpaired probability to an accessibility score.

    The unpaired probability from RNAplfold is already in [0, 1],
    so this is essentially a pass-through with clamping.

    Parameters
    ----------
    unpaired_prob : float
        Probability of being unpaired (0.0–1.0).

    Returns
    -------
    float
        Accessibility score in [0.0, 1.0].
    """
    return max(0.0, min(1.0, unpaired_prob))


def score_window_accessibility(
    profile: AccessibilityProfile,
    start: int,
    end: int,
    aggregation: str = "mean",
) -> float:
    """Aggregate accessibility across a binding site window.

    Computes a single accessibility score for the region [start, end]
    (1-based, inclusive) by aggregating per-position scores.

    Parameters
    ----------
    profile : AccessibilityProfile
        Pre-computed accessibility profile.
    start : int
        Start position (1-based, inclusive).
    end : int
        End position (1-based, inclusive).
    aggregation : str
        Aggregation method: "mean", "median", "min", or "max".
        "mean" is recommended as it balances overall accessibility.
        "min" is conservative (worst position determines the score).

    Returns
    -------
    float
        Aggregated accessibility score in [0.0, 1.0].

    Raises
    ------
    ValueError
        If start/end are out of range or aggregation is unknown.
    """
    if start < 1 or end > len(profile.unpaired_probs) or start > end:
        raise ValueError(
            f"Invalid window [{start}, {end}] for sequence of length {len(profile.unpaired_probs)}"
        )

    window_probs = profile.unpaired_probs[start - 1 : end]
    scores = [score_position_accessibility(p) for p in window_probs]

    if aggregation == "mean":
        return sum(scores) / len(scores)
    elif aggregation == "median":
        sorted_scores = sorted(scores)
        n = len(sorted_scores)
        if n % 2 == 0:
            return (sorted_scores[n // 2 - 1] + sorted_scores[n // 2]) / 2
        return sorted_scores[n // 2]
    elif aggregation == "min":
        return min(scores)
    elif aggregation == "max":
        return max(scores)
    else:
        raise ValueError(
            f"Unknown aggregation method: {aggregation}. Use 'mean', 'median', 'min', or 'max'."
        )


def classify_accessibility(
    score: float,
    accessible_threshold: float = 0.5,
    moderate_threshold: float = 0.3,
) -> str:
    """Classify a binding site by accessibility score.

    Parameters
    ----------
    score : float
        Aggregated accessibility score (0.0–1.0).
    accessible_threshold : float
        Scores at or above this are "accessible" (default 0.5).
    moderate_threshold : float
        Scores at or above this (but below accessible_threshold) are
        "moderate" (default 0.3). Below this is "inaccessible".

    Returns
    -------
    str
        One of "accessible", "moderate", or "inaccessible".
    """
    if score >= accessible_threshold:
        return "accessible"
    elif score >= moderate_threshold:
        return "moderate"
    else:
        return "inaccessible"


def export_profile(
    profile: AccessibilityProfile,
    output_path: Path | str,
    fmt: str = "json",
) -> Path:
    """Export an accessibility profile for visualization.

    Parameters
    ----------
    profile : AccessibilityProfile
        Pre-computed accessibility profile.
    output_path : Path or str
        Output file path.
    fmt : str
        Output format: "json" or "csv".

    Returns
    -------
    Path
        Path to the exported file.

    Raises
    ------
    ValueError
        If the format is not supported.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if fmt == "json":
        data = {
            "sequence": profile.sequence,
            "window_size": profile.window_size,
            "positions": profile.positions,
            "unpaired_probs": profile.unpaired_probs,
            "accessibility_scores": [
                score_position_accessibility(p) for p in profile.unpaired_probs
            ],
        }
        with open(output_path, "w") as f:
            json.dump(data, f, indent=2)

    elif fmt == "csv":
        import csv

        with open(output_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["position", "unpaired_prob", "accessibility_score"])
            for pos, prob in zip(profile.positions, profile.unpaired_probs, strict=True):
                writer.writerow([pos, prob, score_position_accessibility(prob)])

    else:
        raise ValueError(f"Unsupported format: {fmt}. Use 'json' or 'csv'.")

    logger.info("Exported accessibility profile to %s", output_path)
    return output_path
