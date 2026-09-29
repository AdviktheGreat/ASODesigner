"""Target RNA structure prediction and accessibility analysis.

Uses ViennaRNA wrappers (RNAfold, RNAplfold) to compute:
- Full transcript MFE secondary structure
- Local accessibility profile (per-position unpaired probabilities)
- Normalized accessibility scores for ASO binding sites

Results are cached to disk to avoid expensive re-computation.
"""

from __future__ import annotations

import hashlib
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
