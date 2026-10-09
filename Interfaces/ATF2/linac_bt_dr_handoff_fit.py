"""Fit a measured BT-IPZT to a declared DR-entry transverse handoff for the twin.

The current SAD/RF-Track handoff is a design-optics baseline.  A digital twin
needs to replace it with a map fitted from a dither data set: each row contains
the BT endpoint coordinates inferred from upstream diagnostics and the DR
injection coordinates measured after the injection pulse.  This module fits
the affine 4D map without opening a control-system connection.

    It deliberately does not infer longitudinal timing, RF phase, or physical
apertures.  Those require independent data and remain separate calibration
inputs.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Sequence

import numpy as np

from Interfaces.ATF2.simulate_linac_bt_dr_handoff import TransverseHandoff


@dataclass(frozen=True)
class HandoffFitResult:
    """Auditable least-squares result for a 4D affine injection handoff."""

    handoff: TransverseHandoff
    samples: int
    design_matrix_rank: int
    condition_number: float
    rms_residual_mm_mrad: float
    max_abs_residual_mm_or_mrad: float

    def as_dict(self) -> dict[str, object]:
        return {
            "samples": self.samples,
            "design_matrix_rank": self.design_matrix_rank,
            "condition_number": self.condition_number,
            "rms_residual_mm_mrad": self.rms_residual_mm_mrad,
            "max_abs_residual_mm_or_mrad": self.max_abs_residual_mm_or_mrad,
            "handoff": self.handoff.as_dict(),
        }


@dataclass(frozen=True)
class HandoffEvaluation:
    """Out-of-sample transverse residuals for an exported handoff map.

    Coordinates retain the handoff convention ``[mm, mrad, mm, mrad]``.
    The aggregate norms are therefore diagnostics, not a unitless merit;
    downstream users should set their own BPM/monitor uncertainty weights.
    """

    samples: int
    rms_residual_mm_mrad: float
    max_abs_residual_mm_or_mrad: float
    mean_residual_mm_mrad: np.ndarray
    rms_residual_by_coordinate_mm_mrad: np.ndarray

    def as_dict(self) -> dict[str, object]:
        return {
            "samples": self.samples,
            "rms_residual_mm_mrad": self.rms_residual_mm_mrad,
            "max_abs_residual_mm_or_mrad": self.max_abs_residual_mm_or_mrad,
            "mean_residual_mm_mrad": self.mean_residual_mm_mrad.tolist(),
            "rms_residual_by_coordinate_mm_mrad": (
                self.rms_residual_by_coordinate_mm_mrad.tolist()
            ),
            "coordinate_order": "x_mm,xp_mrad,y_mm,yp_mrad",
        }


def evaluate_transverse_handoff(
    handoff: TransverseHandoff,
    source_coordinates_mm_mrad: Sequence[Sequence[float]],
    observed_target_coordinates_mm_mrad: Sequence[Sequence[float]],
) -> HandoffEvaluation:
    """Score a fitted map on independent BT-endpoint/DR-injection samples.

    This intentionally does no refit.  Supply held-out dither shots (or a
    later machine period) to measure the transverse part of twin-model drift.
    It cannot diagnose longitudinal timing, RF capture or an aperture model.
    """
    source = np.asarray(source_coordinates_mm_mrad, dtype=float)
    target = np.asarray(observed_target_coordinates_mm_mrad, dtype=float)
    if source.ndim != 2 or source.shape[1:] != (4,) or source.shape[0] < 1:
        raise ValueError("evaluation source coordinates must have shape (samples>=1, 4)")
    if target.shape != source.shape:
        raise ValueError("evaluation target coordinates must have the same shape as source")
    if not np.all(np.isfinite(source)) or not np.all(np.isfinite(target)):
        raise ValueError("evaluation source and target coordinates must be finite")
    predicted = source @ np.asarray(handoff.matrix, dtype=float).T
    predicted += np.asarray(handoff.offset_mm_mrad, dtype=float)
    residual = target - predicted
    return HandoffEvaluation(
        samples=int(source.shape[0]),
        rms_residual_mm_mrad=float(np.sqrt(np.mean(np.square(residual))),),
        max_abs_residual_mm_or_mrad=float(np.max(np.abs(residual))),
        mean_residual_mm_mrad=np.mean(residual, axis=0),
        rms_residual_by_coordinate_mm_mrad=np.sqrt(np.mean(np.square(residual), axis=0)),
    )


def fit_transverse_handoff(
    source_coordinates_mm_mrad: Sequence[Sequence[float]],
    observed_target_coordinates_mm_mrad: Sequence[Sequence[float]],
    *,
    baseline_dispersion_handoff: TransverseHandoff | None = None,
    rcond: float | None = None,
    target_location: str = "KII.1",
    provenance: str | None = None,
) -> HandoffFitResult:
    """Fit ``target = matrix @ source + offset`` from a dither data set.

    At least five linearly independent source settings are required: four
    transverse coordinates plus a constant offset.  If the design handoff has
    a stated first-order dispersion treatment, pass it as
    ``baseline_dispersion_handoff`` to retain that *separately calibrated*
    energy mapping; transverse dither data alone cannot identify dispersion.
    ``target_location`` must name the same RF-Track entry boundary used later
    by the tracking pipeline (normally ``KII.1``).
    """
    source = np.asarray(source_coordinates_mm_mrad, dtype=float)
    target = np.asarray(observed_target_coordinates_mm_mrad, dtype=float)
    target_location = str(target_location).upper()
    if source.ndim != 2 or source.shape[1:] != (4,):
        raise ValueError("source coordinates must have shape (samples, 4)")
    if target.shape != source.shape:
        raise ValueError("target coordinates must have the same shape as source")
    if source.shape[0] < 5:
        raise ValueError("at least five dither samples are required for a 4D affine map")
    if not np.all(np.isfinite(source)) or not np.all(np.isfinite(target)):
        raise ValueError("source and target coordinates must be finite")

    design = np.column_stack((source, np.ones(source.shape[0])))
    rank = int(np.linalg.matrix_rank(design))
    if rank < 5:
        raise ValueError(
            "dither source coordinates are rank deficient; vary all four transverse "
            "coordinates independently before fitting the handoff"
        )
    coefficients, _, _, singular_values = np.linalg.lstsq(design, target, rcond=rcond)
    matrix = coefficients[:4, :].T
    offset = coefficients[4, :]
    predicted = source @ matrix.T + offset
    residual = target - predicted
    condition = float(singular_values[0] / singular_values[-1])
    kwargs: dict[str, object] = {}
    if baseline_dispersion_handoff is not None:
        if baseline_dispersion_handoff.target_location.upper() != target_location:
            raise ValueError(
                "baseline dispersion handoff target_location does not match the fitted target: "
                f"{baseline_dispersion_handoff.target_location!r} != {target_location!r}"
            )
        if (
            baseline_dispersion_handoff.source_dispersion_mm_mrad is None
            or baseline_dispersion_handoff.target_dispersion_mm_mrad is None
            or baseline_dispersion_handoff.reference_momentum_mev_c is None
        ):
            raise ValueError("baseline dispersion handoff must contain both dispersions and momentum")
        kwargs = {
            "source_dispersion_mm_mrad": baseline_dispersion_handoff.source_dispersion_mm_mrad,
            "target_dispersion_mm_mrad": baseline_dispersion_handoff.target_dispersion_mm_mrad,
            "reference_momentum_mev_c": baseline_dispersion_handoff.reference_momentum_mev_c,
        }
    resolved_provenance = (
        provenance
        or f"least-squares fitted BT IPZT-to-DR {target_location} transverse handoff"
    )
    if baseline_dispersion_handoff is not None:
        resolved_provenance += "; retained first-order dispersion from supplied baseline"
    handoff = TransverseHandoff(
        matrix=matrix,
        offset_mm_mrad=offset,
        provenance=resolved_provenance,
        target_location=target_location,
        **kwargs,
    )
    return HandoffFitResult(
        handoff=handoff,
        samples=int(source.shape[0]),
        design_matrix_rank=rank,
        condition_number=condition,
        rms_residual_mm_mrad=float(np.sqrt(np.mean(np.square(residual)))),
        max_abs_residual_mm_or_mrad=float(np.max(np.abs(residual))),
    )


__all__ = [
    "HandoffEvaluation", "HandoffFitResult", "evaluate_transverse_handoff",
    "fit_transverse_handoff",
]


def main() -> None:
    """Fit a user-owned dither JSON without modifying controls or files.

    The input schema is ``{"source_coordinates_mm_mrad": [[...], ...],
    "observed_target_coordinates_mm_mrad": [[...], ...],
    "validation_source_coordinates_mm_mrad": [[...], ...]?,
    "validation_observed_target_coordinates_mm_mrad": [[...], ...]?,
    "baseline_dispersion_handoff": {...}?}``.  The optional baseline is a
    ``TransverseHandoff.as_dict()`` payload used only to retain separately
    calibrated first-order dispersion.
    """
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument("input_json", help="Measured BT-endpoint/DR-injection dither data")
    args = parser.parse_args()
    payload = json.loads(Path(args.input_json).read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("handoff-fit input must be a JSON object")
    allowed = {
        "source_coordinates_mm_mrad", "observed_target_coordinates_mm_mrad",
        "validation_source_coordinates_mm_mrad",
        "validation_observed_target_coordinates_mm_mrad",
        "baseline_dispersion_handoff", "target_location",
    }
    unknown = set(payload) - allowed
    required = {"source_coordinates_mm_mrad", "observed_target_coordinates_mm_mrad"}
    if unknown or required - set(payload):
        raise ValueError(
            f"handoff-fit input keys invalid: unknown={sorted(unknown)}, "
            f"missing={sorted(required - set(payload))}"
        )
    baseline_payload = payload.get("baseline_dispersion_handoff")
    baseline = (
        None if baseline_payload is None else TransverseHandoff.from_dict(baseline_payload)
    )
    target_location = str(payload.get("target_location", "KII.1"))
    result = fit_transverse_handoff(
        payload["source_coordinates_mm_mrad"],
        payload["observed_target_coordinates_mm_mrad"],
        baseline_dispersion_handoff=baseline,
        target_location=target_location,
        provenance=(
            f"least-squares fitted BT IPZT-to-DR {target_location} transverse handoff "
            f"from {Path(args.input_json).name}"
        ),
    )
    validation_source = payload.get("validation_source_coordinates_mm_mrad")
    validation_target = payload.get("validation_observed_target_coordinates_mm_mrad")
    if (validation_source is None) != (validation_target is None):
        raise ValueError(
            "validation_source_coordinates_mm_mrad and "
            "validation_observed_target_coordinates_mm_mrad must be supplied together"
        )
    output = result.as_dict()
    if validation_source is not None:
        output["held_out_validation"] = evaluate_transverse_handoff(
            result.handoff, validation_source, validation_target
        ).as_dict()
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
