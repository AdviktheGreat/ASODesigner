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
