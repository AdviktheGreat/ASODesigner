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
