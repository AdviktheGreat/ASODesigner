"""Self-complementarity and homodimer screening for ASO candidates.

Unintended ASO:ASO hybridization is a classic ASO failure mode, and the
dominant kernel is 3'-terminal complementarity: the 3' end of one ASO
copy pairing with (any part of) another copy, which extended the duplex
and blocks RNase H recruitment. This module detects that motif with a
pure-Python sequence scan and can optionally quantify the homodimer
with RNAduplex.
"""

from dataclasses import dataclass

from aso_designer.utils.logger import get_logger

logger = get_logger(__name__)

# Watson-Crick pairs on the DNA alphabet (U is normalized to T first).
_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

# G·T wobble pairs form in DNA/DNA duplexes and are treated as pairing.
_WOBBLE = {("G", "T"), ("T", "G")}


def _pairs(a: str, b: str, allow_wobble: bool) -> bool:
    """Whether two DNA bases pair in an antiparallel duplex."""
    if _COMPLEMENT.get(a) == b:
        return True
    return allow_wobble and (a, b) in _WOBBLE


def _normalize(sequence: str) -> str:
    return sequence.upper().replace("U", "T").replace(" ", "").strip()


def find_3p_complementarity(
    sequence: str,
    min_stretch: int = 4,
    allow_wobble: bool = True,
) -> tuple[int, int] | None:
    """Find a 3'-anchored self-complementary stretch in a sequence.

    A homodimer forms when the 3'-terminal ``k`` nucleotides of one ASO
    copy pair with ``k`` contiguous nucleotides of another copy. This
    function tests every stretch length ``k`` from the full sequence
    down to ``min_stretch``: the reverse complement of the 3'-terminal
    ``k``-mer must pair with some contiguous window of the sequence
    (including the 3' tail itself, i.e. a palindromic tail, which lets
    two copies pair tail-to-tail).

    Parameters
    ----------
    sequence : str
        DNA or RNA sequence (case-insensitive, U normalized to T).
    min_stretch : int
        Minimum stretch length to report (default 4).
    allow_wobble : bool
        Whether G·T wobble pairs count as pairing (default True).

    Returns
    -------
    tuple[int, int] or None
        ``(stretch_length, window_start)`` where ``window_start`` is the
        0-based start of the paired internal window, or None when no
        stretch of at least ``min_stretch`` exists.
    """
    seq = _normalize(sequence)
    n = len(seq)
    if n < 2 or min_stretch < 1:
        return None
    for k in range(n, max(min_stretch, 1) - 1, -1):
        tail_start = n - k
        for window_start in range(0, n - k + 1):
            # Tail base n-k+t pairs with window base a+k-1-t (antiparallel).
            if all(
                _pairs(seq[tail_start + t], seq[window_start + k - 1 - t], allow_wobble)
                for t in range(k)
            ):
                return (k, window_start)
    return None


@dataclass
class DimerScreenResult:
    """Result of the ASO self-dimer screen.

    Attributes
    ----------
    sequence : str
        Normalized DNA sequence screened.
    three_prime_stretch : int
        Length of the longest 3'-anchored self-complementary stretch
        (0 when none was found).
    three_prime_position : int or None
        0-based start of the paired internal window, or None.
    homodimer_dg : float or None
        RNAduplex MFE of the ASO against itself (kcal/mol), when
        computed.
    flagged : bool
        True when the candidate fails the screen.
    detail : str
        Human-readable explanation of the decision.
    """

    sequence: str
    three_prime_stretch: int
    three_prime_position: int | None
    homodimer_dg: float | None
    flagged: bool
    detail: str


def run_self_duplex(sequence: str) -> float:
    """Compute the homodimer MFE of a sequence with RNAduplex.

    Runs RNAduplex with the sequence against itself (intermolecular
    only, intramolecular folding ignored).

    Parameters
    ----------
    sequence : str
        DNA or RNA sequence (case-insensitive).

    Returns
    -------
    float
        Homodimer MFE in kcal/mol (negative = stable dimer).

    Raises
    ------
    ViennaRNAError
        If RNAduplex fails or output cannot be parsed.
    """
    from aso_designer.viennarna import run_rnaduplex

    result = run_rnaduplex(sequence, sequence)
    logger.debug("Homodimer MFE for %s: %.2f kcal/mol", sequence, result.mfe)
    return result.mfe


def screen_self_dimer(
    sequence: str,
    min_stretch: int = 4,
    dg_threshold: float = -6.0,
    compute_dg: bool = False,
    allow_wobble: bool = True,
) -> DimerScreenResult:
    """Screen a sequence for self-dimer risk.

    The sequence is flagged when its 3' terminus is self-complementary
    over at least ``min_stretch`` nucleotides, or (when ``compute_dg``
    is set) when the RNAduplex homodimer MFE is at or below
    ``dg_threshold``.

    Parameters
    ----------
    sequence : str
        DNA or RNA sequence (case-insensitive).
    min_stretch : int
        Minimum 3'-complementary stretch length that flags the
        candidate (default 4).
    dg_threshold : float
        Homodimer MFE at or below which the candidate is flagged,
        when the energy is computed (default -6.0 kcal/mol).
    compute_dg : bool
        Run RNAduplex in addition to the sequence scan (default False;
        requires ViennaRNA when True).
    allow_wobble : bool
        Whether G·T wobble pairs count as pairing (default True).

    Returns
    -------
    DimerScreenResult
        Screen outcome with details.
    """
    seq = _normalize(sequence)
    hit = find_3p_complementarity(seq, min_stretch=1, allow_wobble=allow_wobble)
    stretch, position = hit if hit else (0, None)

    homodimer_dg: float | None = None
    if compute_dg:
        homodimer_dg = run_self_duplex(seq)

    flagged = min_stretch > 0 and stretch >= min_stretch
    reason: list[str] = []
    if flagged:
        reason.append(f"3' self-complementary stretch of {stretch} nt")
    if homodimer_dg is not None and homodimer_dg <= dg_threshold:
        flagged = True
        reason.append(f"homodimer MFE {homodimer_dg:.2f} kcal/mol <= {dg_threshold}")

    detail = "No 3'-anchored self-complementarity detected" if not flagged else "; ".join(reason)
    return DimerScreenResult(
        sequence=seq,
        three_prime_stretch=stretch,
        three_prime_position=position,
        homodimer_dg=homodimer_dg,
        flagged=flagged,
        detail=detail,
    )
