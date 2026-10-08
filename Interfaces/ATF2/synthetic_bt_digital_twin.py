"""Offline, synthetic model-mismatch exercise for the ATF BT digital twin.

It creates a nominal ``digital twin`` and a separate, hidden ``pseudo
machine``.  The pseudo machine has static corrector offsets and quadrupole
scale errors.  An orbit-response matrix measured on the nominal model is used
to predict and apply one safe correction to the pseudo machine.

This is not evidence of real-machine performance: replace the pseudo-machine
orbit and response measurement with archived/live data only after independent
BPM, corrector, and energy calibration is available.  It is intentionally
reference-particle based; finite-bunch transmission and aperture loss must be
evaluated in a separate experiment.
"""

from __future__ import annotations

import json

import numpy as np

from Interfaces.ATF2.InterfaceATF2_BT_RFTrack import InterfaceATF2_BT_RFTrack


RESPONSE_STEP_T_MM = 0.002
SVD_RCOND = 1e-3
MAX_CORRECTION_T_MM = 0.05
PSEUDO_MACHINE_QUAD_SCALES = {
    "QD10T": 1.03,
    "QF11T": 0.98,
    "QF20T": 1.02,
}
PSEUDO_MACHINE_CORRECTOR_OFFSETS_T_MM = {
    "ZX10T": 0.015,
    "ZV11T": -0.012,
    "ZH30T": 0.010,
    "ZV50T": -0.008,
}


def _orbit(interface: InterfaceATF2_BT_RFTrack) -> np.ndarray:
    bpms = interface.get_bpms()
    return np.r_[bpms["x"][0], bpms["y"][0]]


def _rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(values))))


def measure_response_matrix(
    interface: InterfaceATF2_BT_RFTrack,
    step: float = RESPONSE_STEP_T_MM,
) -> np.ndarray:
    """Central-difference orbit response in mm/(T mm), restoring all settings."""
    correctors = interface.get_correctors()
    names = correctors["names"]
    origin = np.asarray(correctors["bdes"], dtype=float)
    columns = []
    for index in range(len(names)):
        plus = origin.copy()
        plus[index] += step
        interface.set_correctors(names, plus)
        orbit_plus = _orbit(interface)
        minus = origin.copy()
        minus[index] -= step
        interface.set_correctors(names, minus)
        orbit_minus = _orbit(interface)
        columns.append((orbit_plus - orbit_minus) / (2.0 * step))
    interface.set_correctors(names, origin)
    return np.column_stack(columns)


def run_synthetic_experiment() -> dict[str, object]:
    """Run a reproducible model-vs-pseudo-machine correction experiment."""
    model = InterfaceATF2_BT_RFTrack(reference_particle=True)
    pseudo_machine = InterfaceATF2_BT_RFTrack(reference_particle=True)

    for name, scale in PSEUDO_MACHINE_QUAD_SCALES.items():
        original = float(pseudo_machine.get_quadrupoles(name)["bdes"][0])
        pseudo_machine.set_quadrupoles(name, original * scale)

    correctors = pseudo_machine.get_correctors()
    names = correctors["names"]
    plant_settings = np.asarray(correctors["bdes"], dtype=float)
    for index, name in enumerate(names):
        plant_settings[index] += PSEUDO_MACHINE_CORRECTOR_OFFSETS_T_MM.get(name, 0.0)
    pseudo_machine.set_correctors(names, plant_settings)

    model_response = measure_response_matrix(model)
    pseudo_response = measure_response_matrix(pseudo_machine)
    measured_orbit = _orbit(pseudo_machine)
    correction = np.linalg.pinv(model_response, rcond=SVD_RCOND) @ (-measured_orbit)
    correction = np.clip(correction, -MAX_CORRECTION_T_MM, MAX_CORRECTION_T_MM)
    prediction = measured_orbit + model_response @ correction

    pseudo_machine.set_correctors(
        names,
        pseudo_machine.get_correctors()["bdes"] + correction,
    )
    corrected_orbit = _orbit(pseudo_machine)
    response_mismatch = float(
        np.linalg.norm(pseudo_response - model_response) / np.linalg.norm(pseudo_response)
    )
    return {
        "simulation_only": True,
        "observables": "BPM x/y orbit only; reference-particle response",
        "pseudo_machine_errors": {
            "quadrupole_scales": PSEUDO_MACHINE_QUAD_SCALES,
            "hidden_corrector_offsets_t_mm": PSEUDO_MACHINE_CORRECTOR_OFFSETS_T_MM,
        },
        "response_matrix": {
            "shape": list(model_response.shape),
            "rank": int(np.linalg.matrix_rank(model_response)),
            "relative_model_mismatch": response_mismatch,
        },
        "correction": {
            "svd_rcond": SVD_RCOND,
            "max_abs_t_mm": float(np.max(np.abs(correction))),
            "predicted_rms_mm": _rms(prediction),
        },
        "orbit_rms_mm": {
            "before": _rms(measured_orbit),
            "after": _rms(corrected_orbit),
        },
        "finite_bunch_transmission": "not evaluated in this reference-particle ORM test",
        "next_calibration_observables": [
            "BPM offset/gain/roll calibration",
            "corrector calibration",
            "quadrupole strength and alignment",
            "dispersion and phase advance",
            "aperture/loss-monitor transmission",
        ],
    }


if __name__ == "__main__":
    print(json.dumps(run_synthetic_experiment(), indent=2))
