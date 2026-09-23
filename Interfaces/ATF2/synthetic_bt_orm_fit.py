"""Fit selected BT digital-twin parameters from a synthetic measured ORM.

This offline identifiability exercise uses the calibrated-BPM pseudo-machine
ORM from :mod:`synthetic_bt_calibration_audit`.  It jointly fits three chosen
quadrupole scale factors and all corrector calibration factors by linearizing
the nominal ORM.  The hidden values are printed only to validate the synthetic
procedure; in a real fit, uncertainty and independent observables (dispersion,
phase, magnet measurements) are required to resolve degeneracies.
"""

from __future__ import annotations

import json

import numpy as np

from Interfaces.ATF2.InterfaceATF2_BT_RFTrack import InterfaceATF2_BT_RFTrack
from Interfaces.ATF2.synthetic_bt_calibration_audit import BTPseudoMachine, _measure_response, _physical_orbit


FIT_QUADRUPOLES = ("QD10T", "QF11T", "QF20T")
QUAD_RELATIVE_STEP = 1.0e-3
SVD_RCOND = 1.0e-5


def _model_response(machine: InterfaceATF2_BT_RFTrack) -> np.ndarray:
    return _measure_response(
        lambda command: machine.set_correctors(machine.corrs, command),
        lambda: _physical_orbit(machine),
        len(machine.corrs),
    )


def run_orm_fit(seed: int = 20260902) -> dict[str, object]:
    """Fit ORM-sensitive magnet and actuator parameters in a pseudo experiment."""
    pseudo = BTPseudoMachine(seed)
    measured_response = _measure_response(
        pseudo.set_commands,
        lambda: pseudo.observe(calibrated=True),
        len(pseudo.corrector_names),
    )
    model = InterfaceATF2_BT_RFTrack(reference_particle=True)
    nominal_response = _model_response(model)
    nominal_quad_values = {
        name: float(model.get_quadrupoles(name)["bdes"][0]) for name in FIT_QUADRUPOLES
    }

    derivative_columns = []
    for name in FIT_QUADRUPOLES:
        model.set_quadrupoles(name, nominal_quad_values[name] * (1.0 + QUAD_RELATIVE_STEP))
        perturbed = _model_response(model)
        derivative_columns.append(
            ((perturbed - nominal_response) / QUAD_RELATIVE_STEP).ravel()
        )
        model.set_quadrupoles(name, nominal_quad_values[name])

    # A corrector calibration factor rescales the corresponding ORM column.
    for column in range(nominal_response.shape[1]):
        derivative = np.zeros_like(nominal_response)
        derivative[:, column] = nominal_response[:, column]
        derivative_columns.append(derivative.ravel())
    sensitivity = np.column_stack(derivative_columns)
    fitted_delta = np.linalg.pinv(sensitivity, rcond=SVD_RCOND) @ (
        measured_response - nominal_response
    ).ravel()
    fitted_quad_scales = {
        name: float(1.0 + fitted_delta[index])
        for index, name in enumerate(FIT_QUADRUPOLES)
    }
    fitted_corrector_gain = 1.0 + fitted_delta[len(FIT_QUADRUPOLES):]

    for name, scale in fitted_quad_scales.items():
        model.set_quadrupoles(name, nominal_quad_values[name] * scale)
    fitted_response = _model_response(model) * fitted_corrector_gain[np.newaxis, :]
    initial_mismatch = float(
        np.linalg.norm(measured_response - nominal_response) / np.linalg.norm(measured_response)
    )
    final_mismatch = float(
        np.linalg.norm(measured_response - fitted_response) / np.linalg.norm(measured_response)
    )
    hidden_quad_scales = pseudo.errors.quadrupole_scales
    hidden_gain = pseudo.errors.corrector_gain
    quad_sensitivity_norm = {
        name: float(np.linalg.norm(derivative_columns[index]))
        for index, name in enumerate(FIT_QUADRUPOLES)
    }
    return {
        "simulation_only": True,
        "observables": "calibrated BPM orbit response matrix only",
        "fit_parameters": {
            "quadrupole_scales": list(FIT_QUADRUPOLES),
            "corrector_gain_count": len(fitted_corrector_gain),
            "linearization_step_relative": QUAD_RELATIVE_STEP,
        },
        "relative_orm_mismatch": {
            "before_fit": initial_mismatch,
            "after_fit": final_mismatch,
        },
        "quadrupole_scale": {
            name: {
                "fitted": fitted_quad_scales[name],
                "hidden_synthetic_truth": hidden_quad_scales[name],
                "error": fitted_quad_scales[name] - hidden_quad_scales[name],
                "orm_sensitivity_norm": quad_sensitivity_norm[name],
                "observable_from_this_orm": quad_sensitivity_norm[name] > 1.0e-10,
            }
            for name in FIT_QUADRUPOLES
        },
        "corrector_gain": {
            "rms_error": float(np.sqrt(np.mean((fitted_corrector_gain - hidden_gain) ** 2))),
            "max_abs_error": float(np.max(np.abs(fitted_corrector_gain - hidden_gain))),
        },
        "important_limit": (
            "A real ORM-only fit can have degeneracies.  Add dispersion, phase advance, "
            "independent BPM calibration, and magnet-current priors before interpreting "
            "these parameters as physical machine errors."
        ),
    }


if __name__ == "__main__":
    print(json.dumps(run_orm_fit(), indent=2))
