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


@dataclass
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
            # Try to adjust gap size to fit
            remaining = length - 2 * self.wing_length
            if remaining < 4:
                raise ConfigurationError(
                    f"Gapmer length {length} too short for wing_length "
                    f"{self.wing_length} (need >= {2 * self.wing_length + 4})"
                )
            self.gap_size = remaining
            logger.debug("Adjusted gap size to %d for length %d", self.gap_size, length)

        wing_mod = MOD_MOE if self.wing_chemistry == "moe" else MOD_2OME
        pattern = (
            [wing_mod] * self.wing_length
            + [MOD_DNA] * self.gap_size
            + [wing_mod] * self.wing_length
        )
        return pattern

    def is_rnase_h_eligible(self) -> bool:
        return True


@dataclass
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
        # Mixmers may or may not recruit RNase H depending on gap continuity.
        # By default, assume not (no contiguous DNA gap).
        return False


@dataclass
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
        return False


# ---------------------------------------------------------------------------
# Chemistry factory and application
# ---------------------------------------------------------------------------

from aso_designer.candidate import ASOCandidate  # noqa: E402
from aso_designer.config_schema import CandidateConfig, GapmerConfig, MixmerConfig  # noqa: E402


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
