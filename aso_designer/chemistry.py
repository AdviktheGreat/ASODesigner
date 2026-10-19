"""ASO chemistry modeling and modification annotation.

Defines the five supported ASO chemistry types and their modification
patterns. Designed from scratch for ASO-specific needs — GuideX has no
chemistry concept since Cas13 gRNAs are unmodified RNA.

Supported chemistries:
- gapmer: DNA center (RNase H gap) flanked by modified wings
- mixmer: interspersed modified and unmodified nucleotides
- 2ome: full 2'-O-methyl modified RNA
- moe: full 2'-O-methoxyethyl modified RNA
- morpholino: phosphorodiamidate morpholino oligomers
"""

from dataclasses import dataclass
from typing import Protocol

from aso_designer.candidate import ASOCandidate
from aso_designer.config_schema import CandidateConfig, GapmerConfig, MixmerConfig
from aso_designer.utils.exceptions import ConfigurationError
from aso_designer.utils.logger import get_logger

logger = get_logger(__name__)

# Modification type labels
MOD_DNA = "DNA"
MOD_2OME = "2OMe"
MOD_MOE = "MOE"
MOD_MORPHOLINO = "MORPHOLINO"

# All valid chemistry names
VALID_CHEMISTRIES = {"gapmer", "mixmer", "2ome", "moe", "morpholino"}


class ChemistryDesign(Protocol):
    """Protocol for chemistry design strategies."""

    def build_pattern(self, length: int) -> list[str]:
        """Return per-nucleotide modification labels for the given length."""
        ...

    def is_rnase_h_eligible(self) -> bool:
        """Whether this chemistry can recruit RNase H."""
        ...


@dataclass(frozen=True)
class GapmerDesign:
    """Gapmer chemistry: DNA center + modified wings.

    The DNA center (gap) recruits RNase H for target cleavage.
    The modified wings increase binding affinity and nuclease
    resistance.

    Parameters
    ----------
    gap_size : int
        Number of DNA nucleotides in the center.
    wing_length : int
        Number of modified nucleotides on each side.
    wing_chemistry : str
        Wing modification type: "2ome" or "moe".
    """

    gap_size: int = 10
    wing_length: int = 5
    wing_chemistry: str = "moe"

    def build_pattern(self, length: int) -> list[str]:
        """Build the modification pattern for a gapmer.

        Layout: [wing] [DNA gap] [wing]
        Total length = 2 * wing_length + gap_size

        Parameters
        ----------
        length : int
            Total ASO length (must match 2*wing_length + gap_size).

        Returns
        -------
        list[str]
            Per-nucleotide modification labels.

        Raises
        ------
        ConfigurationError
            If the length doesn't match the expected gapmer layout.
        """
        expected = 2 * self.wing_length + self.gap_size
        if length != expected:
            raise ConfigurationError(
                f"Invalid gapmer geometry: length {length} != "
                f"2 * wing_length {self.wing_length} + gap_size {self.gap_size} "
                f"(expected {expected}). Adjust candidate.length, "
                "gapmer.gap_size, or gapmer.wing_length so they are consistent."
            )

        wing_mod = MOD_MOE if self.wing_chemistry == "moe" else MOD_2OME
        pattern = (
            [wing_mod] * self.wing_length
            + [MOD_DNA] * self.gap_size
            + [wing_mod] * self.wing_length
        )
        return pattern

    def is_rnase_h_eligible(self) -> bool:
        """Whether this gapmer chemistry can recruit RNase H (True for gapmers)."""
        return True


@dataclass(frozen=True)
class MixmerDesign:
    """Mixmer chemistry: interspersed modified and unmodified nucleotides.

    Mixmers do not have a contiguous DNA gap, so they typically do not
    recruit RNase H. They are used for steric blocking and splice
    modulation.

    Parameters
    ----------
    pattern : list[bool]
        Modification pattern (True = modified, False = DNA).
        Length must match the ASO length.
    modification : str
        Modification type for marked positions: "2ome" or "moe".
    """

    pattern: list[bool] | None = None
    modification: str = "2ome"

    def build_pattern(self, length: int) -> list[str]:
        """Build the modification pattern for a mixmer.

        If a custom pattern is provided, it must match the ASO length.
        Otherwise, a default alternating pattern is generated.

        Parameters
        ----------
        length : int
            Total ASO length.

        Returns
        -------
        list[str]
            Per-nucleotide modification labels.
        """
        mod_label = MOD_MOE if self.modification == "moe" else MOD_2OME

        if self.pattern is not None:
            if len(self.pattern) != length:
                raise ConfigurationError(
                    f"Mixmer pattern length ({len(self.pattern)}) does not "
                    f"match ASO length ({length})"
                )
            return [mod_label if p else MOD_DNA for p in self.pattern]

        # Default: alternating every other nucleotide (starting with modified)
        default_pattern = []
        for i in range(length):
            default_pattern.append(mod_label if i % 2 == 0 else MOD_DNA)
        return default_pattern

    def is_rnase_h_eligible(self) -> bool:
        """Whether this mixmer can recruit RNase H (False — no contiguous DNA gap)."""
        return False


@dataclass(frozen=True)
class FullChemistryDesign:
    """Full-modification chemistry (2'-OMe, MOE, or morpholino).

    Every nucleotide is modified with the same chemistry. These do
    not recruit RNase H (steric blocking mechanism).

    Parameters
    ----------
    mod_type : str
        Modification type: "2ome", "moe", or "morpholino".
    """

    mod_type: str = "2ome"

    def build_pattern(self, length: int) -> list[str]:
        """Build the modification pattern (all positions modified).

        Parameters
        ----------
        length : int
            Total ASO length.

        Returns
        -------
        list[str]
            Per-nucleotide modification labels.
        """
        label = {
            "2ome": MOD_2OME,
            "moe": MOD_MOE,
            "morpholino": MOD_MORPHOLINO,
        }.get(self.mod_type, MOD_2OME)
        return [label] * length

    def is_rnase_h_eligible(self) -> bool:
        """Whether full-modification chemistry can recruit RNase H (False — steric blocking)."""
        return False


# ---------------------------------------------------------------------------
# Chemistry factory and application
# ---------------------------------------------------------------------------


def create_chemistry_design(
    config: CandidateConfig,
    gapmer_config: GapmerConfig | None = None,
    mixmer_config: MixmerConfig | None = None,
) -> ChemistryDesign:
    """Create the appropriate chemistry design from configuration.

    Parameters
    ----------
    config : CandidateConfig
        Candidate configuration (contains chemistry type and length).
    gapmer_config : GapmerConfig, optional
        Gapmer-specific parameters (used if chemistry is "gapmer").
    mixmer_config : MixmerConfig, optional
        Mixmer-specific parameters (used if chemistry is "mixmer").

    Returns
    -------
    ChemistryDesign
        The chemistry design strategy for the configured chemistry type.

    Raises
    ------
    ConfigurationError
        If the chemistry type is invalid or required config is missing.
    """
    chemistry = config.chemistry

    if chemistry == "gapmer":
        if gapmer_config is None:
            gapmer_config = GapmerConfig()
        expected = 2 * gapmer_config.wing_length + gapmer_config.gap_size
        if config.length != expected:
            raise ConfigurationError(
                f"Invalid gapmer geometry: candidate.length {config.length} != "
                f"2 * wing_length {gapmer_config.wing_length} + gap_size "
                f"{gapmer_config.gap_size} (expected {expected}). Adjust "
                "candidate.length, gapmer.gap_size, or gapmer.wing_length."
            )
        return GapmerDesign(
            gap_size=gapmer_config.gap_size,
            wing_length=gapmer_config.wing_length,
            wing_chemistry=gapmer_config.wing_chemistry,
        )

    elif chemistry == "mixmer":
        if mixmer_config is None:
            mixmer_config = MixmerConfig()
        return MixmerDesign(
            pattern=mixmer_config.pattern if mixmer_config.pattern else None,
            modification=mixmer_config.modification,
        )

    elif chemistry in ("2ome", "moe", "morpholino"):
        return FullChemistryDesign(mod_type=chemistry)

    else:
        raise ConfigurationError(f"Unknown chemistry type: '{chemistry}'")


def apply_chemistry(
    candidates: list[ASOCandidate],
    design: ChemistryDesign,
) -> None:
    """Apply a chemistry design to all candidates (in place).

    Sets the ``modification_pattern`` on each candidate and stores
    RNase H eligibility in the features dict.

    Parameters
    ----------
    candidates : list[ASOCandidate]
        Candidates to annotate (modified in place).
    design : ChemistryDesign
        The chemistry design strategy to apply.
    """
    rnase_h = design.is_rnase_h_eligible()
    for candidate in candidates:
        candidate.modification_pattern = design.build_pattern(candidate.length)
        candidate.features["rnase_h_eligible"] = rnase_h

    logger.info(
        "Applied %s chemistry to %d candidates (RNase H: %s)",
        candidate.chemistry,
        len(candidates),
        rnase_h,
    )


# ---------------------------------------------------------------------------
# Mixmer custom pattern utilities
# ---------------------------------------------------------------------------

# Scheme characters for custom mixmer patterns
_SCHEME_MODIFIED = {"M", "m"}  # Modified nucleotide
_SCHEME_DNA = {"D", "d", "."}  # DNA (unmodified)


def parse_mixmer_scheme(scheme: str) -> list[bool]:
    """Parse a mixmer scheme string into a boolean modification pattern.

    Scheme characters:
    - 'M' or 'm': modified nucleotide (True)
    - 'D', 'd', or '.': DNA / unmodified (False)

    Example: "MMDMMDDMMDM" → [True, True, False, True, True, False, False, True, True, False, True]

    Parameters
    ----------
    scheme : str
        Scheme string (case-insensitive).

    Returns
    -------
    list[bool]
        Boolean pattern where True = modified.

    Raises
    ------
    ConfigurationError
        If the scheme contains invalid characters.
    """
    pattern: list[bool] = []
    for i, char in enumerate(scheme):
        if char in _SCHEME_MODIFIED:
            pattern.append(True)
        elif char in _SCHEME_DNA:
            pattern.append(False)
        else:
            raise ConfigurationError(
                f"Invalid character '{char}' at position {i} in mixmer scheme. "
                f"Use M (modified), D or . (DNA)."
            )
    return pattern


def has_continuous_dna_gap(pattern: list[bool], min_gap: int = 4) -> bool:
    """Check if a mixmer pattern has a contiguous DNA gap sufficient for RNase H.

    RNase H requires a minimum contiguous DNA stretch (typically ≥4 nt)
    to recruit and cleave the RNA target.

    Parameters
    ----------
    pattern : list[bool]
        Boolean modification pattern (False = DNA).
    min_gap : int
        Minimum contiguous DNA stretch for RNase H (default 4).

    Returns
    -------
    bool
        True if the pattern has a contiguous DNA gap of at least min_gap.
    """
    current_gap = 0
    for is_modified in pattern:
        if not is_modified:
            current_gap += 1
            if current_gap >= min_gap:
                return True
        else:
            current_gap = 0
    return False


def mixmer_from_scheme(
    scheme: str,
    modification: str = "2ome",
) -> MixmerDesign:
    """Create a MixmerDesign from a scheme string.

    Parameters
    ----------
    scheme : str
        Scheme string (e.g. "MMDMMDDMMDM").
    modification : str
        Modification type for modified positions: "2ome" or "moe".

    Returns
    -------
    MixmerDesign
        Configured mixmer design with the parsed pattern.
    """
    pattern = parse_mixmer_scheme(scheme)
    return MixmerDesign(pattern=pattern, modification=modification)


# ---------------------------------------------------------------------------
# Chemistry annotation utilities
# ---------------------------------------------------------------------------


def validate_pattern(candidate: ASOCandidate) -> bool:
    """Validate that a candidate's modification pattern matches its length.

    Parameters
    ----------
    candidate : ASOCandidate
        Candidate with a modification pattern to validate.

    Returns
    -------
    bool
        True if the pattern length matches the sequence length.
    """
    if not candidate.modification_pattern:
        return False
    return len(candidate.modification_pattern) == candidate.length


def get_dna_gap_region(candidate: ASOCandidate) -> tuple[int, int] | None:
    """Get the (start, end) of the contiguous DNA gap in a candidate.

    For gapmers, this is the DNA center region that recruits RNase H.
    For mixmers, this is the longest contiguous DNA stretch (if any).
    For full-modification chemistries, returns None.

    Parameters
    ----------
    candidate : ASOCandidate
        Candidate with a modification pattern.

    Returns
    -------
    tuple[int, int] or None
        (start, end) indices within the ASO sequence, or None if no
        contiguous DNA gap exists.
    """
    if not candidate.modification_pattern:
        return None

    best_start = None
    best_len = 0
    current_start = None
    current_len = 0

    for i, mod in enumerate(candidate.modification_pattern):
        if mod == MOD_DNA:
            if current_start is None:
                current_start = i
            current_len += 1
        else:
            if current_len > best_len:
                best_start = current_start
                best_len = current_len
            current_start = None
            current_len = 0

    # Check the last stretch
    if current_len > best_len:
        best_start = current_start
        best_len = current_len

    if best_start is None or best_len == 0:
        return None

    return (best_start, best_start + best_len)


def chemistry_summary(candidate: ASOCandidate) -> str:
    """Return a human-readable summary of a candidate's chemistry.

    Parameters
    ----------
    candidate : ASOCandidate
        Candidate with a modification pattern.

    Returns
    -------
    str
        Summary string (e.g. "gapmer 5-10-5 (MOE/DNA/MOE), RNase H: yes").
    """
    if not candidate.modification_pattern:
        return f"{candidate.chemistry} (no pattern)"

    # Count modification types
    mod_counts: dict[str, int] = {}
    for mod in candidate.modification_pattern:
        mod_counts[mod] = mod_counts.get(mod, 0) + 1

    parts: list[str] = [candidate.chemistry]

    # For gapmers, show wing-gap-wing layout
    if candidate.chemistry == "gapmer":
        gap = get_dna_gap_region(candidate)
        if gap:
            wing = candidate.length - (gap[1] - gap[0])
            wing_each = wing // 2
            parts.append(f"{wing_each}-{gap[1] - gap[0]}-{wing_each}")

    # Show modification types
    mod_str = "/".join(f"{k}({v})" for k, v in sorted(mod_counts.items()))
    parts.append(mod_str)

    # RNase H status
    rnase_h = candidate.features.get("rnase_h_eligible", False)
    parts.append(f"RNase H: {'yes' if rnase_h else 'no'}")

    return ", ".join(parts)


# ---------------------------------------------------------------------------
# RNase H1 positional cleavage-likelihood model
# ---------------------------------------------------------------------------

# Gap-size factor for RNase H1 cleavage efficiency (heuristic, in [0, 1]).
# Human RNase H1 requires a stretch of ~4-6 contiguous DNA nucleotides for
# productive cleavage; gaps of ~6-10 nt are near-optimal, and very long
# gaps add little while costing affinity/specificity.
_GAP_SIZE_FACTOR: dict[int, float] = {
    1: 0.0,
    2: 0.0,
    3: 0.1,
    4: 0.4,
    5: 0.7,
    6: 0.95,
    7: 1.0,
    8: 1.0,
    9: 1.0,
    10: 0.9,
    11: 0.85,
    12: 0.8,
}
_LONG_GAP_FACTOR = 0.75

# Relative cleavage propensity as a function of position within the DNA
# gap (fraction of the gap from its 5' end). Cleavage sites cluster in
# the center-to-3' portion of the gap, and the terminal gap nucleotides
# are rarely cut because RNase H needs duplex substrate on both sides of
# the scissile phosphate.
_EDGE_POSITION_FACTOR = 0.15
_PEAK_FRACTION = 0.65


def _position_factor(fraction: float) -> float:
    """Within-gap positional cleavage propensity in [0, 1].

    A smooth profile that rises from the 5' gap edge, peaks at
    ``_PEAK_FRACTION`` of the gap, and falls toward the 3' edge, with
    reduced likelihood at both terminal gap positions.
    """
    if fraction <= 0.0 or fraction >= 1.0:
        return _EDGE_POSITION_FACTOR
    distance = abs(fraction - _PEAK_FRACTION)
    half_width = max(_PEAK_FRACTION, 1.0 - _PEAK_FRACTION)
    return max(0.2, 1.0 - (distance / half_width) * 0.8)


def cleavage_likelihood(
    position_in_gap: int,
    gap_length: int,
    gap_at_5p_end: bool = False,
    gap_at_3p_end: bool = False,
) -> float:
    """Positional RNase H1 cleavage likelihood for one gap nucleotide.

    Combines a gap-size factor (cleavage needs a contiguous DNA stretch;
    very short gaps are not substrates) with a within-gap positional
    profile (cleavage sites cluster center-to-3' within the gap, and the
    terminal gap positions are disfavored) and a junction-context
    penalty when the gap runs into the end of the ASO (no flanking
    duplex for RNase H).

    Parameters
    ----------
    position_in_gap : int
        0-based position within the DNA gap (from the gap's 5' end).
    gap_length : int
        Total length of the DNA gap.
    gap_at_5p_end : bool
        Whether the gap touches the 5' terminus of the ASO.
    gap_at_3p_end : bool
        Whether the gap touches the 3' terminus of the ASO.

    Returns
    -------
    float
        Cleavage likelihood in [0.0, 1.0] for the target RNA base
        opposite this gap position.
    """
    if gap_length < 1 or position_in_gap < 0 or position_in_gap >= gap_length:
        return 0.0

    size_factor = _GAP_SIZE_FACTOR.get(gap_length, _LONG_GAP_FACTOR)

    fraction = position_in_gap / (gap_length - 1) if gap_length > 1 else 0.0
    position_score = _position_factor(fraction)

    # Without a flanking wing on one side, cleavage near that junction has
    # no intact heteroduplex substrate beyond the gap.
    if gap_at_5p_end and position_in_gap <= 1:
        position_score *= 0.1
    if gap_at_3p_end and position_in_gap >= gap_length - 2:
        position_score *= 0.1

    return max(0.0, min(1.0, size_factor * position_score))


@dataclass
class CleavageProfile:
    """Positional RNase H1 cleavage profile across an ASO.

    Attributes
    ----------
    likelihoods : list[float]
        Per-ASO-position cleavage likelihood (0.0 at non-DNA positions).
    gap_length : int
        Length of the largest contiguous DNA gap (0 if none).
    max_score : float
        Highest per-position likelihood.
    expected_site : int | None
        0-based ASO position with the highest likelihood, if any.
    """

    likelihoods: list[float]
    gap_length: int
    max_score: float
    expected_site: int | None

    @property
    def mean_gap_score(self) -> float:
        """Mean likelihood over gap positions (0.0 without a gap)."""
        scores = [x for x in self.likelihoods if x > 0.0]
        return sum(scores) / len(scores) if scores else 0.0


def compute_cleavage_profile(modification_pattern: list[str]) -> CleavageProfile:
    """Compute the positional RNase H cleavage profile for a pattern.

    Finds the largest contiguous DNA gap in the modification pattern and
    assigns every gap position a cleavage likelihood via
    ``cleavage_likelihood``; all other positions score 0.0.

    Parameters
    ----------
    modification_pattern : list[str]
        Per-nucleotide modification labels (``MOD_DNA`` = unmodified DNA).

    Returns
    -------
    CleavageProfile
        Per-position likelihoods and gap summary.
    """
    n = len(modification_pattern)
    likelihoods = [0.0] * n

    # Locate all contiguous DNA runs, keep the largest.
    best_start, best_len = None, 0
    run_start, run_len = None, 0
    for i, mod in enumerate(modification_pattern):
        if mod == MOD_DNA:
            if run_start is None:
                run_start = i
            run_len += 1
            if run_len > best_len:
                best_start, best_len = run_start, run_len
        else:
            run_start, run_len = None, 0
    if best_start is None or best_len == 0:
        return CleavageProfile(
            likelihoods=likelihoods, gap_length=0, max_score=0.0, expected_site=None
        )

    gap_at_5p = best_start == 0
    gap_at_3p = best_start + best_len == n
    for offset in range(best_len):
        likelihoods[best_start + offset] = cleavage_likelihood(
            offset, best_len, gap_at_5p_end=gap_at_5p, gap_at_3p_end=gap_at_3p
        )

    max_score = max(likelihoods)
    expected_site = likelihoods.index(max_score) if max_score > 0.0 else None
    return CleavageProfile(
        likelihoods=likelihoods,
        gap_length=best_len,
        max_score=max_score,
        expected_site=expected_site,
    )
