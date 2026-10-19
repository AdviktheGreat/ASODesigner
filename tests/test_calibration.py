"""Tests for the weight-calibration utility."""

from pathlib import Path
from types import SimpleNamespace

import pytest

from aso_designer.calibration import (
    COMPONENT_KEYS,
    calibrate_weights,
    load_calibration_csv,
    recommend_weights,
)
from aso_designer.utils.exceptions import ConfigurationError


def _write_csv(path: Path, n: int = 30, accessibility_signal: float = 40.0) -> Path:
    """Synthetic data: activity driven almost entirely by accessibility."""
    import csv
    import random

    rng = random.Random(42)
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([*COMPONENT_KEYS, "activity", "notes"])
        for _ in range(n):
            row = {k: round(rng.uniform(0.2, 0.9), 3) for k in COMPONENT_KEYS}
            activity = accessibility_signal * row["accessibility"] + rng.gauss(0, 1)
            writer.writerow([*row.values(), round(activity, 2), "x"])
    return path


def test_load_calibration_csv_reads_rows(tmp_path):
    path = _write_csv(tmp_path / "data.csv")
    rows = load_calibration_csv(path)
    assert len(rows) == 30
    assert set(COMPONENT_KEYS) <= set(rows[0])


def test_load_calibration_csv_missing_columns(tmp_path):
    path = tmp_path / "bad.csv"
    path.write_text("accessibility,activity\n0.5,10\n")
    with pytest.raises(ConfigurationError, match="missing columns"):
        load_calibration_csv(path)


def test_load_calibration_csv_missing_file(tmp_path):
    with pytest.raises(ConfigurationError, match="not found"):
        load_calibration_csv(tmp_path / "nope.csv")


def test_calibrate_weights_favors_informative_component(tmp_path):
    path = _write_csv(tmp_path / "data.csv")
    result = calibrate_weights(load_calibration_csv(path))
    assert abs(sum(result.weights.values()) - 1.0) < 1e-3
    assert all(v >= 0.0 for v in result.weights.values())
    assert result.weights["accessibility"] == max(result.weights.values())
    assert result.r_squared > 0.5
    assert result.n_samples == 30


def test_calibrate_weights_rejects_negative_ridge(tmp_path):
    rows = load_calibration_csv(_write_csv(tmp_path / "data.csv"))
    with pytest.raises(ConfigurationError, match="Ridge penalty"):
        calibrate_weights(rows, ridge=-1.0)


def test_recommend_weights_round_trips_through_schema(tmp_path):
    path = _write_csv(tmp_path / "data.csv")
    result = recommend_weights(path)
    from aso_designer.config_schema import ScoringWeightsConfig

    cfg = ScoringWeightsConfig(**result.weights)
    assert abs(sum(cfg.model_dump().values()) - 1.0) < 0.01


def test_cli_calibrate_writes_yaml(tmp_path, capsys):
    path = _write_csv(tmp_path / "data.csv")
    out = tmp_path / "weights.yaml"
    from aso_designer.cli import _handle_calibrate

    args = SimpleNamespace(csv=str(path), ridge=1.0, out=str(out))

    rc = _handle_calibrate(args)
    assert rc == 0
    text = out.read_text()
    assert text.startswith("scoring:")
    assert "accessibility:" in text
    captured = capsys.readouterr()
    assert "Recommended weights" in captured.out


def test_cli_calibrate_reports_bad_csv(tmp_path, capsys):
    from aso_designer.cli import _handle_calibrate

    args = SimpleNamespace(csv=str(tmp_path / "missing.csv"), ridge=1.0, out=None)

    rc = _handle_calibrate(args)
    assert rc == 1
    assert "not found" in capsys.readouterr().err
