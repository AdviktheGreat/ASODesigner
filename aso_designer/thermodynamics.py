"""Thermodynamic calculations for ASO candidate evaluation.

Implements nearest-neighbor (NN) model for melting temperature (Tm)
prediction, duplex ΔG, self-structure ΔG, target structure ΔG, and
net binding energy. Uses SantaLucia 1998 unified NN parameters for
DNA/DNA and Turner 2004 parameters for RNA/RNA duplexes. ViennaRNA
tools (RNAfold, RNAduplex, RNAcofold) are used for structure-based
ΔG calculations where the NN model is insufficient.

Scientific note: Tm and ΔG values are computational predictions for
candidate ranking, not experimental measurements. The NN model
assumes 1 M NaCl monovalent salt and non-self-complementary sequences
unless otherwise specified.
"""

from dataclasses import dataclass
from pathlib import Path

from aso_designer.utils.exceptions import ThermodynamicsError
from aso_designer.utils.logger import get_logger

logger = get_logger(__name__)

# Gas constant in cal/(mol·K)
R = 1.987

# Reference temperature for ΔG calculations (37 °C in Kelvin)
T_REF = 310.15  # 37 °C

# ---------------------------------------------------------------------------
# SantaLucia 1998 unified NN parameters for DNA/DNA
# Reference: SantaLucia (1998) PNAS 95:1460-1465
# Units: ΔH in kcal/mol, ΔS in cal/(mol·K)
# Key format: "XY" = 5'-X Y-3' / 3'-Y' X'-5' (top strand 5'→3')
# ---------------------------------------------------------------------------

_DNA_NN: dict[str, tuple[float, float]] = {
    # NN pair (5'-3'/3'-5'): (ΔH, ΔS)
    # 10 unique parameters from SantaLucia 1998, expanded to all 16
    # dinucleotides using reverse-complement symmetry.
    "AA": (-7.9, -22.2),
    "AC": (-8.4, -22.4),  # = GT
    "AG": (-7.8, -21.0),  # = CT
    "AT": (-7.2, -20.4),
    "CA": (-8.5, -22.7),
    "CC": (-8.0, -19.9),  # = GG
    "CG": (-10.6, -27.2),
    "CT": (-7.8, -21.0),
    "GA": (-8.2, -22.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9),
    "GT": (-8.4, -22.4),
    "TA": (-7.2, -21.3),
    "TC": (-8.2, -22.2),  # = GA
    "TG": (-8.5, -22.7),  # = CA
    "TT": (-7.9, -22.2),  # = AA
}

# Initiation parameters (SantaLucia 1998)
# Terminal G·C: ΔH = 0.1, ΔS = -2.8
# Terminal A·T: ΔH = 2.3, ΔS = 4.1
# Symmetry correction (self-complementary): ΔH = 0, ΔS = -1.4
_DNA_INIT_GC = (0.1, -2.8)
_DNA_INIT_AT = (2.3, 4.1)
_DNA_SYMMETRY = (0.0, -1.4)

# ---------------------------------------------------------------------------
# Turner 2004 NN parameters for RNA/RNA
# Reference: Turner & Mathews (2010) Nucleic Acids Res 38:D280-D282
# (2004 parameters as curated in the NNDB)
# Units: ΔH in kcal/mol, ΔS in cal/(mol·K)
# ---------------------------------------------------------------------------

_RNA_NN: dict[str, tuple[float, float]] = {
    # NN pair (5'-3'/3'-5'): (ΔH, ΔS)
    "AA": (-6.82, -19.0),
    "AC": (-8.40, -22.4),
    "AG": (-7.60, -19.9),
    "AU": (-7.20, -20.4),
    "CA": (-8.40, -22.4),
    "CC": (-8.00, -19.9),
    "CG": (-10.60, -27.2),
    "CU": (-7.60, -19.9),
    "GA": (-7.60, -20.5),
    "GC": (-13.20, -32.9),
    "GG": (-8.00, -19.9),
    "GU": (-8.40, -22.4),
    "UA": (-7.20, -21.3),
    "UC": (-7.60, -19.9),
    "UG": (-8.40, -22.4),
    "UU": (-6.82, -19.0),
}

# RNA initiation parameters (Turner 2004)
# Terminal G·C: ΔH = 0.0, ΔS = -4.0 (approximate)
# Terminal A·U: ΔH = 0.0, ΔS = -1.0 (approximate)
_RNA_INIT_GC = (0.0, -4.0)
_RNA_INIT_AT = (0.0, -1.0)
# Symmetry correction for self-complementary RNA
_RNA_SYMMETRY = (0.0, -1.4)


def _get_nn_params(rna: bool = False) -> dict[str, tuple[float, float]]:
    """Return the appropriate NN parameter table."""
    return _RNA_NN if rna else _DNA_NN


def _get_init_params(
    rna: bool = False,
) -> tuple[tuple[float, float], tuple[float, float], tuple[float, float]]:
    """Return (init_GC, init_AT, symmetry) parameters."""
    if rna:
        return _RNA_INIT_GC, _RNA_INIT_AT, _RNA_SYMMETRY
    return _DNA_INIT_GC, _DNA_INIT_AT, _DNA_SYMMETRY


def _is_self_complementary(sequence: str, rna: bool = False) -> bool:
    """Check if a sequence is self-complementary (palindromic)."""
    seq = sequence.upper()
    if rna:
        seq = seq.replace("T", "U")
    comp_map = str.maketrans("ACGU", "UGCA") if rna else str.maketrans("ACGT", "TGCA")
    complement = seq.translate(comp_map)
    return seq == complement[::-1]


def calculate_nn_thermodynamics(
    sequence: str,
    rna: bool = False,
) -> tuple[float, float]:
    """Calculate total ΔH and ΔS using the nearest-neighbor model.

    Sums nearest-neighbor stacking enthalpies and entropies, plus
    initiation and symmetry corrections.

    Parameters
    ----------
    sequence : str
        Oligonucleotide sequence (DNA or RNA alphabet).
    rna : bool
        If True, use RNA/RNA NN parameters (Turner 2004).
        If False, use DNA/DNA NN parameters (SantaLucia 1998).

    Returns
    -------
    tuple[float, float]
        (Total ΔH in kcal/mol, Total ΔS in cal/(mol·K))

    Raises
    ------
    ThermodynamicsError
        If the sequence is too short or contains invalid NN pairs.
    """
    seq = sequence.upper()
    if rna:
        seq = seq.replace("T", "U")

    if len(seq) < 2:
        raise ThermodynamicsError(f"Sequence too short for NN calculation: {len(seq)} nt")

    nn_params = _get_nn_params(rna)
    init_gc, init_at, symmetry = _get_init_params(rna)

    total_h = 0.0
    total_s = 0.0

    # Sum NN stacking parameters
    for i in range(len(seq) - 1):
        nn = seq[i : i + 2]
        if nn not in nn_params:
            raise ThermodynamicsError(f"Invalid nearest-neighbor pair: {nn} in sequence {seq}")
        h, s = nn_params[nn]
        total_h += h
        total_s += s

    # Initiation correction: terminal base pair
    # 5' terminal
    first_base = seq[0]
    if first_base in ("G", "C"):
        total_h += init_gc[0]
        total_s += init_gc[1]
    else:
        total_h += init_at[0]
        total_s += init_at[1]

    # 3' terminal
    last_base = seq[-1]
    if last_base in ("G", "C"):
        total_h += init_gc[0]
        total_s += init_gc[1]
    else:
        total_h += init_at[0]
        total_s += init_at[1]

    # Symmetry correction
    if _is_self_complementary(seq, rna=rna):
        total_h += symmetry[0]
        total_s += symmetry[1]

    return total_h, total_s


# ---------------------------------------------------------------------------
# Melting temperature (Tm) calculation
# ---------------------------------------------------------------------------


def calculate_tm(
    sequence: str,
    na_conc: float = 0.1,
    oligo_conc: float = 1.0e-6,
    rna: bool = False,
    salt_correction: str = "santalucia",
) -> float:
    """Calculate the melting temperature (Tm) using the NN model.

    Tm = (ΔH * 1000) / (ΔS + R * ln(C_T / x)) - 273.15

    where:
    - ΔH = total enthalpy (kcal/mol)
    - ΔS = total entropy (cal/(mol·K))
    - R = 1.987 cal/(mol·K)
    - C_T = total strand concentration (M)
    - x = 4 for non-self-complementary, 1 for self-complementary

    Parameters
    ----------
    sequence : str
        Oligonucleotide sequence (DNA or RNA alphabet).
    na_conc : float
        Sodium concentration in M (default 0.1 M = 100 mM).
    oligo_conc : float
        Total oligo strand concentration in M (default 1 µM).
    rna : bool
        If True, use RNA/RNA NN parameters.
    salt_correction : str
        Salt correction method: "santalucia" (1998) or "log10" (simple).

    Returns
    -------
    float
        Melting temperature in °C.

    Raises
    ------
    ThermodynamicsError
        If NN calculation fails or salt correction is unknown.
    """
    import math

    dh, ds = calculate_nn_thermodynamics(sequence, rna=rna)

    # Symmetry factor: 1 for self-complementary, 4 otherwise
    seq = sequence.upper()
    if rna:
        seq = seq.replace("T", "U")
    is_self_comp = _is_self_complementary(seq, rna=rna)
    x = 1 if is_self_comp else 4

    # Tm at 1 M NaCl
    # ΔH in kcal/mol → * 1000 to get cal/mol
    # ΔS in cal/(mol·K)
    # R * ln(C_T / x) in cal/(mol·K)
    denominator = ds + R * math.log(oligo_conc / x)
    if denominator == 0:
        raise ThermodynamicsError("Denominator is zero in Tm calculation")

    tm_1m = (dh * 1000) / denominator - 273.15

    # Apply salt correction
    if salt_correction == "santalucia":
        # SantaLucia 1998: Tm(salt) = Tm(1M) + (1/ΔS) * 0.368 * (N-1) * ln([Na+])
        # N = number of phosphates = length of sequence
        n_phosphates = len(seq)
        if ds != 0:
            salt_term = (1.0 / ds) * 0.368 * (n_phosphates - 1) * math.log(na_conc)
            tm = tm_1m + salt_term
        else:
            tm = tm_1m
    elif salt_correction == "log10":
        # Simple: Tm(salt) = Tm(1M) + 16.6 * log10([Na+])
        tm = tm_1m + 16.6 * math.log10(na_conc)
    else:
        raise ThermodynamicsError(
            f"Unknown salt correction method: {salt_correction}. Use 'santalucia' or 'log10'."
        )

    return tm


def calculate_dg(
    sequence: str,
    temperature: float = 37.0,
    rna: bool = False,
) -> float:
    """Calculate ΔG at a given temperature using the NN model.

    ΔG = ΔH - T * ΔS / 1000

    where T is in Kelvin and ΔS is converted from cal to kcal.

    Parameters
    ----------
    sequence : str
        Oligonucleotide sequence.
    temperature : float
        Temperature in °C (default 37 °C).
    rna : bool
        If True, use RNA/RNA NN parameters.

    Returns
    -------
    float
        ΔG in kcal/mol (negative = favorable).
    """
    dh, ds = calculate_nn_thermodynamics(sequence, rna=rna)
    t_kelvin = temperature + 273.15
    # ΔS is in cal/(mol·K), ΔH is in kcal/mol → convert ΔS to kcal
    dg = dh - (t_kelvin * ds / 1000.0)
    return dg


# ---------------------------------------------------------------------------
# Chemistry-aware Tm calculation and scoring
# ---------------------------------------------------------------------------

# Chemistry types that use RNA NN parameters (target is RNA, ASO is RNA-like)
_RNA_CHEMISTRIES = {"2ome", "moe"}


def get_chemistry_tm_adjustment(chemistry: str) -> float:
    """Get the Tm adjustment for a given chemistry from presets.

    Parameters
    ----------
    chemistry : str
        Chemistry name: gapmer, mixmer, 2ome, moe, morpholino.

    Returns
    -------
    float
        Tm adjustment in °C (added to computed Tm).
    """
    import yaml

    preset_path = Path(__file__).parent / "config" / "chemistry_presets.yaml"
    try:
        with open(preset_path) as f:
            presets = yaml.safe_load(f) or {}
    except Exception:
        logger.warning("Could not load chemistry presets; using 0.0 adjustment")
        return 0.0

    preset = presets.get(chemistry, {})
    return preset.get("tm_adjustment", 0.0)


def calculate_tm_with_chemistry(
    sequence: str,
    chemistry: str = "gapmer",
    na_conc: float = 0.1,
    oligo_conc: float = 1.0e-6,
) -> float:
    """Calculate Tm with chemistry-specific NN parameters and adjustment.

    Selects the appropriate NN parameter set based on chemistry:
    - 2'OMe, MOE: RNA NN parameters (ASO is RNA-like, binds RNA target)
    - Gapmer, mixmer, morpholino: DNA NN parameters (DNA or DNA-like backbone)

    Then applies the chemistry-specific Tm adjustment from presets.

    Parameters
    ----------
    sequence : str
        ASO sequence in DNA alphabet (ACGT, will be converted to RNA if needed).
    chemistry : str
        Chemistry name (default "gapmer").
    na_conc : float
        Sodium concentration in M (default 0.1).
    oligo_conc : float
        Oligo concentration in M (default 1 µM).

    Returns
    -------
    float
        Adjusted Tm in °C.
    """
    use_rna = chemistry in _RNA_CHEMISTRIES
    base_tm = calculate_tm(
        sequence,
        na_conc=na_conc,
        oligo_conc=oligo_conc,
        rna=use_rna,
    )
    adjustment = get_chemistry_tm_adjustment(chemistry)
    adjusted_tm = base_tm + adjustment
    logger.debug(
        "Tm for %s (%s): base=%.2f, adj=%+.2f, final=%.2f",
        sequence[:10],
        chemistry,
        base_tm,
        adjustment,
        adjusted_tm,
    )
    return adjusted_tm


def score_tm(
    tm: float,
    tm_optimal: float = 52.0,
    tm_min: float = 40.0,
    tm_max: float = 65.0,
) -> float:
    """Score a Tm value based on distance from optimal.

    Returns a score in [0.0, 1.0]:
    - 1.0 at tm_optimal
    - Linear decrease to 0.5 at tm_min or tm_max
    - 0.0 outside [tm_min, tm_max]

    Parameters
    ----------
    tm : float
        Computed Tm in °C.
    tm_optimal : float
        Optimal Tm in °C (default 52.0).
    tm_min : float
        Minimum acceptable Tm (default 40.0).
    tm_max : float
        Maximum acceptable Tm (default 65.0).

    Returns
    -------
    float
        Tm score in [0.0, 1.0].
    """
    if tm < tm_min or tm > tm_max:
        return 0.0

    if tm == tm_optimal:
        return 1.0

    # Linear interpolation: 1.0 at optimal, 0.5 at boundaries
    if tm < tm_optimal:
        # Below optimal: interpolate between tm_min (0.5) and tm_optimal (1.0)
        if tm_min == tm_optimal:
            return 1.0
        return 0.5 + 0.5 * (tm - tm_min) / (tm_optimal - tm_min)
    else:
        # Above optimal: interpolate between tm_optimal (1.0) and tm_max (0.5)
        if tm_max == tm_optimal:
            return 1.0
        return 1.0 - 0.5 * (tm - tm_optimal) / (tm_max - tm_optimal)


# ---------------------------------------------------------------------------
# Duplex ΔG computation
# ---------------------------------------------------------------------------


@dataclass
class DuplexResult:
    """Result of duplex ΔG computation.

    Attributes
    ----------
    aso_sequence : str
        ASO sequence (DNA alphabet).
    target_sequence : str
        Target subsequence (DNA alphabet).
    method : str
        Method used: "rnaduplex" or "nearest_neighbor".
    dg_duplex : float
        Duplex formation ΔG in kcal/mol (negative = favorable).
    start_target : int
        Start position of hybrid in target (1-based, 0 if NN method).
    end_target : int
        End position of hybrid in target (1-based, 0 if NN method).
    """

    aso_sequence: str
    target_sequence: str
    method: str
    dg_duplex: float
    start_target: int = 0
    end_target: int = 0


def calculate_duplex_dg_rnaduplex(
    aso_sequence: str,
    target_sequence: str,
) -> DuplexResult:
    """Compute duplex ΔG using ViennaRNA RNAduplex.

    RNAduplex computes the minimum free energy of inter-molecular
    hybridization, ignoring intramolecular structure. This is the
    primary method for ASO:target binding energy.

    Parameters
    ----------
    aso_sequence : str
        ASO sequence (DNA or RNA, T → U conversion applied).
    target_sequence : str
        Target subsequence (DNA or RNA, T → U conversion applied).

    Returns
    -------
    DuplexResult
        Duplex ΔG and hybridization positions.

    Raises
    ------
    ThermodynamicsError
        If RNAduplex fails.
    """
    from aso_designer.viennarna import run_rnaduplex

    try:
        result = run_rnaduplex(aso_sequence, target_sequence)
    except Exception as e:
        raise ThermodynamicsError(f"RNAduplex failed for ASO:target duplex: {e}") from e

    return DuplexResult(
        aso_sequence=aso_sequence.upper(),
        target_sequence=target_sequence.upper(),
        method="rnaduplex",
        dg_duplex=result.mfe,
        start_target=result.start2,
        end_target=result.end2,
    )


def calculate_duplex_dg_nn(
    aso_sequence: str,
    chemistry: str = "gapmer",
    temperature: float = 37.0,
) -> DuplexResult:
    """Compute duplex ΔG using the nearest-neighbor model.

    The NN model computes ΔG for the ASO binding to its perfect
    complement. This is faster than RNAduplex but does not account
    for imperfect matches or bulges.

    Parameters
    ----------
    aso_sequence : str
        ASO sequence (DNA alphabet).
    chemistry : str
        Chemistry name (determines RNA vs DNA NN parameters).
    temperature : float
        Temperature in °C (default 37 °C).

    Returns
    -------
    DuplexResult
        Duplex ΔG (NN-based).
    """
    use_rna = chemistry in _RNA_CHEMISTRIES
    dg = calculate_dg(aso_sequence, temperature=temperature, rna=use_rna)
    return DuplexResult(
        aso_sequence=aso_sequence.upper(),
        target_sequence="",
        method="nearest_neighbor",
        dg_duplex=dg,
    )


def calculate_duplex_dg(
    aso_sequence: str,
    target_sequence: str | None = None,
    chemistry: str = "gapmer",
    method: str = "nearest_neighbor",
    temperature: float = 37.0,
) -> DuplexResult:
    """Compute duplex ΔG using the specified method.

    Parameters
    ----------
    aso_sequence : str
        ASO sequence (DNA alphabet).
    target_sequence : str or None
        Target subsequence. Required for "rnaduplex" method.
    chemistry : str
        Chemistry name (for NN method parameter selection).
    method : str
        "nearest_neighbor" (fast, perfect match only) or
        "rnaduplex" (ViennaRNA, handles mismatches).
    temperature : float
        Temperature for NN calculation (default 37 °C).

    Returns
    -------
    DuplexResult
        Duplex ΔG result.
    """
    if method == "rnaduplex":
        if target_sequence is None:
            raise ThermodynamicsError("target_sequence is required for rnaduplex method")
        return calculate_duplex_dg_rnaduplex(aso_sequence, target_sequence)
    elif method == "nearest_neighbor":
        return calculate_duplex_dg_nn(aso_sequence, chemistry=chemistry, temperature=temperature)
    else:
        raise ThermodynamicsError(
            f"Unknown duplex method: {method}. Use 'nearest_neighbor' or 'rnaduplex'."
        )


# ---------------------------------------------------------------------------
# Self-structure ΔG computation
# ---------------------------------------------------------------------------


@dataclass
class SelfStructureResult:
    """Result of self-structure ΔG computation.

    Attributes
    ----------
    sequence : str
        ASO sequence.
    dg_self : float
        Self-folding ΔG in kcal/mol (more negative = more self-structure).
    structure : str
        Dot-bracket secondary structure (empty if NN method).
    method : str
        Method used: "rnafold" or "nearest_neighbor".
    penalty : float
        Self-structure penalty (0.0 if no significant self-structure).
    """

    sequence: str
    dg_self: float
    structure: str
    method: str
    penalty: float = 0.0


def calculate_self_structure_dg_rnafold(
    aso_sequence: str,
) -> SelfStructureResult:
    """Compute self-structure ΔG using ViennaRNA RNAfold.

    RNAfold computes the MFE secondary structure of the ASO alone.
    More negative ΔG indicates stronger self-folding, which reduces
    the effective concentration available for target binding.

    Parameters
    ----------
    aso_sequence : str
        ASO sequence (DNA or RNA, T → U conversion applied).

    Returns
    -------
    SelfStructureResult
        Self-structure ΔG and dot-bracket structure.

    Raises
    ------
    ThermodynamicsError
        If RNAfold fails.
    """
    from aso_designer.viennarna import run_rnafold

    try:
        result = run_rnafold(aso_sequence)
    except Exception as e:
        raise ThermodynamicsError(f"RNAfold failed for self-structure: {e}") from e

    return SelfStructureResult(
        sequence=aso_sequence.upper(),
        dg_self=result.mfe,
        structure=result.structure,
        method="rnafold",
    )


def calculate_self_structure_dg_nn(
    aso_sequence: str,
    chemistry: str = "gapmer",
    temperature: float = 37.0,
) -> SelfStructureResult:
    """Compute self-structure ΔG using the nearest-neighbor model.

    This computes the ΔG of the ASO folding back on itself (hairpin).
    For short ASOs (15-25 nt), self-structure is typically minimal, so
    the NN model provides a quick approximation. For more accurate
    results, use the RNAfold method.

    Parameters
    ----------
    aso_sequence : str
        ASO sequence (DNA alphabet).
    chemistry : str
        Chemistry name (determines NN parameter set).
    temperature : float
        Temperature in °C (default 37 °C).

    Returns
    -------
    SelfStructureResult
        Self-structure ΔG (NN-based, no structure).
    """
    use_rna = chemistry in _RNA_CHEMISTRIES
    dg = calculate_dg(aso_sequence, temperature=temperature, rna=use_rna)
    # For NN self-structure, we use a simple heuristic:
    # self-structure ΔG is a fraction of the total ΔG, since not all
    # of the ASO can fold back on itself. We use 0.3 as the fraction.
    # This is a rough approximation; RNAfold is preferred.
    dg_self = dg * 0.3
    return SelfStructureResult(
        sequence=aso_sequence.upper(),
        dg_self=dg_self,
        structure="",
        method="nearest_neighbor",
    )


def calculate_self_structure_penalty(
    dg_self: float,
    self_structure_max: float = -3.0,
) -> float:
    """Compute a penalty for self-structure.

    The penalty is 0.0 when ΔG ≥ 0 (no self-structure) and increases
    linearly as ΔG becomes more negative. At the threshold
    (self_structure_max), the penalty is 1.0.

    Parameters
    ----------
    dg_self : float
        Self-structure ΔG in kcal/mol.
    self_structure_max : float
        Threshold ΔG for maximum penalty (default -3.0 kcal/mol).

    Returns
    -------
    float
        Penalty in [0.0, 1.0+]. Values >1.0 indicate ΔG worse than threshold.
    """
    if dg_self >= 0:
        return 0.0
    if self_structure_max >= 0:
        return 0.0
    return abs(dg_self) / abs(self_structure_max)


def get_chemistry_self_structure_weight(chemistry: str) -> float:
    """Get the self-structure weight adjustment for a chemistry.

    Some chemistries (e.g., gapmers with DNA center) are less prone
    to self-structure issues, while others (e.g., 2'OMe) are more.

    Parameters
    ----------
    chemistry : str
        Chemistry name.

    Returns
    -------
    float
        Weight multiplier (1.0 = default, <1.0 = less sensitive,
        >1.0 = more sensitive).
    """
    import yaml

    preset_path = Path(__file__).parent / "config" / "chemistry_presets.yaml"
    try:
        with open(preset_path) as f:
            presets = yaml.safe_load(f) or {}
    except Exception:
        return 1.0

    preset = presets.get(chemistry, {})
    return preset.get("self_structure_weight_adjust", 1.0)


def calculate_self_structure_dg(
    aso_sequence: str,
    chemistry: str = "gapmer",
    method: str = "nearest_neighbor",
    self_structure_max: float = -3.0,
    temperature: float = 37.0,
) -> SelfStructureResult:
    """Compute self-structure ΔG and penalty with chemistry awareness.

    Parameters
    ----------
    aso_sequence : str
        ASO sequence (DNA alphabet).
    chemistry : str
        Chemistry name (for weight adjustment and NN parameter selection).
    method : str
        "nearest_neighbor" (fast) or "rnafold" (ViennaRNA, accurate).
    self_structure_max : float
        Threshold for maximum penalty (default -3.0 kcal/mol).
    temperature : float
        Temperature for NN calculation (default 37 °C).

    Returns
    -------
    SelfStructureResult
        Self-structure ΔG, structure, and chemistry-adjusted penalty.
    """
    if method == "rnafold":
        result = calculate_self_structure_dg_rnafold(aso_sequence)
    elif method == "nearest_neighbor":
        result = calculate_self_structure_dg_nn(
            aso_sequence, chemistry=chemistry, temperature=temperature
        )
    else:
        raise ThermodynamicsError(
            f"Unknown self-structure method: {method}. Use 'nearest_neighbor' or 'rnafold'."
        )

    # Compute base penalty
    base_penalty = calculate_self_structure_penalty(result.dg_self, self_structure_max)

    # Apply chemistry-specific weight adjustment
    weight = get_chemistry_self_structure_weight(chemistry)
    result.penalty = base_penalty * weight

    return result


# ---------------------------------------------------------------------------
# Target structure ΔG and net binding energy
# ---------------------------------------------------------------------------


@dataclass
class TargetStructureResult:
    """Result of target structure ΔG computation at a binding site.

    Attributes
    ----------
    target_sequence : str
        Target subsequence at the binding site.
    dg_target : float
        Target folding ΔG in kcal/mol (more negative = more structured).
    structure : str
        Dot-bracket secondary structure (empty if NN method).
    method : str
        Method used: "rnafold" or "nearest_neighbor".
    """

    target_sequence: str
    dg_target: float
    structure: str
    method: str


def calculate_target_structure_dg_rnafold(
    target_sequence: str,
) -> TargetStructureResult:
    """Compute target structure ΔG at a binding site using RNAfold.

    The target subsequence at the ASO binding site may have local
    secondary structure that must be opened for ASO binding. RNAfold
    computes the MFE of this subsequence.

    Parameters
    ----------
    target_sequence : str
        Target subsequence at the binding site (DNA or RNA).

    Returns
    -------
    TargetStructureResult
        Target structure ΔG and dot-bracket structure.
    """
    from aso_designer.viennarna import run_rnafold

    try:
        result = run_rnafold(target_sequence)
    except Exception as e:
        raise ThermodynamicsError(f"RNAfold failed for target structure: {e}") from e

    return TargetStructureResult(
        target_sequence=target_sequence.upper(),
        dg_target=result.mfe,
        structure=result.structure,
        method="rnafold",
    )


def calculate_target_structure_dg_nn(
    target_sequence: str,
    temperature: float = 37.0,
) -> TargetStructureResult:
    """Compute target structure ΔG using the NN model.

    Uses RNA NN parameters since the target is mRNA.

    Parameters
    ----------
    target_sequence : str
        Target subsequence (DNA alphabet, converted to RNA internally).
    temperature : float
        Temperature in °C (default 37 °C).

    Returns
    -------
    TargetStructureResult
        Target structure ΔG (NN-based).
    """
    dg = calculate_dg(target_sequence, temperature=temperature, rna=True)
    return TargetStructureResult(
        target_sequence=target_sequence.upper(),
        dg_target=dg,
        structure="",
        method="nearest_neighbor",
    )


def calculate_target_structure_dg(
    target_sequence: str,
    method: str = "nearest_neighbor",
    temperature: float = 37.0,
) -> TargetStructureResult:
    """Compute target structure ΔG using the specified method.

    Parameters
    ----------
    target_sequence : str
        Target subsequence at the binding site.
    method : str
        "nearest_neighbor" (fast) or "rnafold" (ViennaRNA, accurate).
    temperature : float
        Temperature for NN calculation (default 37 °C).

    Returns
    -------
    TargetStructureResult
        Target structure ΔG result.
    """
    if method == "rnafold":
        return calculate_target_structure_dg_rnafold(target_sequence)
    elif method == "nearest_neighbor":
        return calculate_target_structure_dg_nn(target_sequence, temperature=temperature)
    else:
        raise ThermodynamicsError(
            f"Unknown target structure method: {method}. Use 'nearest_neighbor' or 'rnafold'."
        )


# ---------------------------------------------------------------------------
# Net binding energy
# ---------------------------------------------------------------------------


@dataclass
class ThermodynamicSummary:
    """Complete thermodynamic profile of an ASO candidate.

    Combines Tm, duplex ΔG, self-structure ΔG, target structure ΔG,
    and net binding energy into a single summary object.

    Attributes
    ----------
    tm : float
        Melting temperature in °C.
    tm_score : float
        Tm score (0–1, 1 = optimal).
    dg_duplex : float
        Duplex formation ΔG in kcal/mol.
    dg_self : float
        Self-structure ΔG in kcal/mol.
    self_structure_penalty : float
        Chemistry-adjusted self-structure penalty.
    dg_target : float
        Target structure ΔG at binding site in kcal/mol.
    dg_net : float
        Net binding ΔG = ΔG_duplex - ΔG_self - ΔG_target.
    binding_class : str
        Binding strength classification: "strong", "moderate", or "weak".
    """

    tm: float
    tm_score: float
    dg_duplex: float
    dg_self: float
    self_structure_penalty: float
    dg_target: float
    dg_net: float
    binding_class: str


def calculate_net_binding_dg(
    dg_duplex: float,
    dg_self: float,
    dg_target: float,
) -> float:
    """Compute net binding ΔG from component energies.

    The thermodynamic cycle for ASO binding:
    1. ASO unfolds (cost: -ΔG_self, since ΔG_self is negative)
    2. Target opens at binding site (cost: -ΔG_target)
    3. ASO:target duplex forms (gain: ΔG_duplex)

    ΔG_net = ΔG_duplex - ΔG_self - ΔG_target

    More negative ΔG_net = more favorable binding.

    Parameters
    ----------
    dg_duplex : float
        Duplex formation ΔG (negative = favorable).
    dg_self : float
        Self-structure ΔG (negative = more self-structure).
    dg_target : float
        Target structure ΔG (negative = more structured target).

    Returns
    -------
    float
        Net binding ΔG in kcal/mol.
    """
    return dg_duplex - dg_self - dg_target


def classify_binding(
    dg_net: float,
    strong_threshold: float = -10.0,
    moderate_threshold: float = -5.0,
) -> str:
    """Classify binding strength by net binding ΔG.

    Parameters
    ----------
    dg_net : float
        Net binding ΔG in kcal/mol.
    strong_threshold : float
        ΔG at or below this is "strong" (default -10.0 kcal/mol).
    moderate_threshold : float
        ΔG at or below this (but above strong_threshold) is "moderate"
        (default -5.0 kcal/mol). Above this is "weak".

    Returns
    -------
    str
        One of "strong", "moderate", or "weak".
    """
    if dg_net <= strong_threshold:
        return "strong"
    elif dg_net <= moderate_threshold:
        return "moderate"
    else:
        return "weak"


def compute_thermodynamic_summary(
    aso_sequence: str,
    target_sequence: str | None = None,
    chemistry: str = "gapmer",
    na_conc: float = 0.1,
    oligo_conc: float = 1.0e-6,
    tm_optimal: float = 52.0,
    tm_min: float = 40.0,
    tm_max: float = 65.0,
    self_structure_max: float = -3.0,
    method: str = "nearest_neighbor",
    temperature: float = 37.0,
) -> ThermodynamicSummary:
    """Compute a complete thermodynamic summary for an ASO candidate.

    Combines Tm, duplex ΔG, self-structure ΔG, target structure ΔG,
    and net binding energy into a single summary.

    Parameters
    ----------
    aso_sequence : str
        ASO sequence (DNA alphabet).
    target_sequence : str or None
        Target subsequence at binding site. Required for "rnaduplex"
        duplex method and "rnafold" target method.
    chemistry : str
        Chemistry name (default "gapmer").
    na_conc : float
        Sodium concentration in M (default 0.1).
    oligo_conc : float
        Oligo concentration in M (default 1 µM).
    tm_optimal : float
        Optimal Tm in °C (default 52.0).
    tm_min : float
        Minimum acceptable Tm (default 40.0).
    tm_max : float
        Maximum acceptable Tm (default 65.0).
    self_structure_max : float
        Self-structure threshold (default -3.0 kcal/mol).
    method : str
        Computation method: "nearest_neighbor" (fast) or
        "rnafold"/"rnaduplex" (ViennaRNA, accurate).
    temperature : float
        Temperature for ΔG calculations (default 37 °C).

    Returns
    -------
    ThermodynamicSummary
        Complete thermodynamic profile.
    """
    # Tm
    tm = calculate_tm_with_chemistry(
        aso_sequence, chemistry=chemistry, na_conc=na_conc, oligo_conc=oligo_conc
    )
    tm_score = score_tm(tm, tm_optimal=tm_optimal, tm_min=tm_min, tm_max=tm_max)

    # Duplex ΔG
    if method == "rnaduplex" and target_sequence is not None:
        duplex = calculate_duplex_dg(
            aso_sequence, target_sequence, chemistry=chemistry, method="rnaduplex"
        )
    else:
        duplex = calculate_duplex_dg(
            aso_sequence,
            chemistry=chemistry,
            method="nearest_neighbor",
            temperature=temperature,
        )

    # Self-structure ΔG
    self_struct = calculate_self_structure_dg(
        aso_sequence,
        chemistry=chemistry,
        method=method if method != "rnaduplex" else "rnafold",
        self_structure_max=self_structure_max,
        temperature=temperature,
    )

    # Target structure ΔG
    if target_sequence is not None:
        target_struct = calculate_target_structure_dg(
            target_sequence,
            method=method if method != "rnaduplex" else "rnafold",
            temperature=temperature,
        )
        dg_target = target_struct.dg_target
    else:
        # Without target sequence, assume no target structure cost
        dg_target = 0.0

    # Net binding energy
    dg_net = calculate_net_binding_dg(
        dg_duplex=duplex.dg_duplex,
        dg_self=self_struct.dg_self,
        dg_target=dg_target,
    )

    # Classify binding
    binding_class = classify_binding(dg_net)

    return ThermodynamicSummary(
        tm=tm,
        tm_score=tm_score,
        dg_duplex=duplex.dg_duplex,
        dg_self=self_struct.dg_self,
        self_structure_penalty=self_struct.penalty,
        dg_target=dg_target,
        dg_net=dg_net,
        binding_class=binding_class,
    )
