"""Weight calibration from measured ASO activity data.

Given a CSV of previously tested ASOs with their component scores and a
measured activity column, fit a ridge regression of activity on the
components and convert the coefficients into a recommended weight vector
(non-negative, summing to 1.0).

CSV format (header row required)::

    accessibility,thermodynamics,gc_content,off_target,variant,self_structure,positional,activity
    0.72,0.85,0.60,0.90,1.00,0.40,0.50,68.0
    ...

The ``activity`` column holds the measured effect (for example percent
knockdown); higher is better. Extra columns are ignored.
"""

from __future__ import annotations

import csv
import logging
import math
from dataclasses import dataclass
from pathlib import Path

from aso_designer.utils.exceptions import ConfigurationError

logger = logging.getLogger(__name__)

COMPONENT_KEYS: tuple[str, ...] = (
    "accessibility",
    "thermodynamics",
    "gc_content",
    "off_target",
    "variant",
    "self_structure",
    "positional",
)

_ACTIVITY_KEY = "activity"


@dataclass(frozen=True)
class CalibrationResult:
    """Outcome of a weight calibration run."""

    weights: dict[str, float]
    coefficients: dict[str, float]
    r_squared: float
    n_samples: int
    ridge: float


def load_calibration_csv(path: str | Path) -> list[dict[str, float]]:
    """Read a calibration CSV into a list of component/activity rows.

    Raises
    ------
    ConfigurationError
        If the file is missing, malformed, or missing required columns.
    """
    csv_path = Path(path)
    if not csv_path.exists():
        raise ConfigurationError(f"Calibration CSV not found: {csv_path}")

    with csv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ConfigurationError(f"Calibration CSV is empty: {csv_path}")
        fieldnames = [name.strip() for name in reader.fieldnames]
        missing = [k for k in (*COMPONENT_KEYS, _ACTIVITY_KEY) if k not in fieldnames]
        if missing:
            raise ConfigurationError(f"Calibration CSV missing columns: {', '.join(missing)}")
        rows: list[dict[str, float]] = []
        for line_no, raw in enumerate(reader, start=2):
            try:
                row = {
                    key: float(raw[key])
                    for key in (*COMPONENT_KEYS, _ACTIVITY_KEY)
                    if raw.get(key) not in (None, "")
                }
            except (TypeError, ValueError) as exc:
                raise ConfigurationError(
                    f"Calibration CSV line {line_no}: non-numeric value ({exc})"
                ) from exc
            missing_vals = [k for k in (*COMPONENT_KEYS, _ACTIVITY_KEY) if k not in row]
            if missing_vals:
                raise ConfigurationError(
                    f"Calibration CSV line {line_no}: missing {', '.join(missing_vals)}"
                )
            rows.append(row)
    if len(rows) < len(COMPONENT_KEYS) + 2:
        raise ConfigurationError(
            "Calibration needs at least "
            f"{len(COMPONENT_KEYS) + 2} rows for {len(COMPONENT_KEYS)} components; "
            f"got {len(rows)}."
        )
    return rows


def _solve_ridge(x: list[list[float]], y: list[float], ridge: float) -> list[float]:
    """Solve (XᵀX + λI) w = Xᵀy by Gaussian elimination with pivoting."""
    n_features = len(x[0])
    n = len(x)
    ata = [[0.0] * n_features for _ in range(n_features)]
    aty = [0.0] * n_features
    for i in range(n):
        for a in range(n_features):
            aty[a] += x[i][a] * y[i]
            for b in range(n_features):
                ata[a][b] += x[i][a] * x[i][b]
    for a in range(n_features):
        ata[a][a] += ridge

    m = [row[:] + [aty[j]] for j, row in enumerate(ata)]
    for col in range(n_features):
        pivot = max(range(col, n_features), key=lambda r: abs(m[r][col]))
        if abs(m[pivot][col]) < 1e-12:
            raise ConfigurationError("Calibration matrix is singular; raise --ridge.")
        m[col], m[pivot] = m[pivot], m[col]
        for r in range(n_features):
            if r == col:
                continue
            factor = m[r][col] / m[col][col]
            if factor == 0.0:
                continue
            for c in range(col, n_features + 1):
                m[r][c] -= factor * m[col][c]
    return [m[i][n_features] / m[i][i] for i in range(n_features)]


def calibrate_weights(
    rows: list[dict[str, float]],
    ridge: float = 1.0,
) -> CalibrationResult:
    """Fit ridge regression and derive non-negative weights summing to 1.0.

    Features are standardized before fitting; coefficients are mapped to
    weights by magnitude (absolute value, renormalized) so that any
    direction of influence contributes proportionally.
    """
    if ridge < 0.0:
        raise ConfigurationError(f"Ridge penalty must be >= 0, got {ridge}")
    n = len(rows)
    means = {k: sum(r[k] for r in rows) / n for k in COMPONENT_KEYS}
    stds = {
        k: math.sqrt(sum((r[k] - means[k]) ** 2 for r in rows) / max(n - 1, 1))
        for k in COMPONENT_KEYS
    }
    x = [
        [(r[k] - means[k]) / stds[k] if stds[k] > 1e-12 else 0.0 for k in COMPONENT_KEYS]
        for r in rows
    ]
    mean_y = sum(r[_ACTIVITY_KEY] for r in rows) / n
    # Center the target so the standardized no-intercept model is valid.
    y = [r[_ACTIVITY_KEY] - mean_y for r in rows]

    coefs = _solve_ridge(x, y, ridge)
    magnitudes = [abs(c) for c in coefs]
    total = sum(magnitudes)
    if total < 1e-12:
        raise ConfigurationError(
            "All component coefficients vanished; activity is uncorrelated "
            "with the provided scores."
        )
    weights = {k: round(m / total, 6) for k, m in zip(COMPONENT_KEYS, magnitudes, strict=True)}

    raw_y = [r[_ACTIVITY_KEY] for r in rows]
    ss_tot = sum((v - mean_y) ** 2 for v in raw_y)
    ss_res = sum(
        (yi - sum(xi * ci for xi, ci in zip(xi_row, coefs, strict=True))) ** 2
        for xi_row, yi in zip(x, y, strict=True)
    )
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 1e-12 else 0.0

    return CalibrationResult(
        weights=weights,
        coefficients=dict(zip(COMPONENT_KEYS, coefs, strict=True)),
        r_squared=r_squared,
        n_samples=n,
        ridge=ridge,
    )


def recommend_weights(
    csv_path: str | Path,
    ridge: float = 1.0,
) -> CalibrationResult:
    """Load a calibration CSV and return the recommended weight vector."""
    rows = load_calibration_csv(csv_path)
    result = calibrate_weights(rows, ridge=ridge)
    logger.info("Calibrated weights from %d samples (R^2=%.3f)", result.n_samples, result.r_squared)
    return result
