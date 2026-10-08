"""Offline calibration audit for an ATF BT digital-twin shadow.

The hidden pseudo machine contains magnetic-model errors, corrector gain and
offset errors, and BPM gain/offset/roll errors.  Three response-matrix cases
are evaluated on precisely the same pseudo machine:

* nominal model ORM with raw BPM values;
* calibrated BPM values but nominal model ORM;
* calibrated BPM values with an ORM measured on the pseudo machine.

The final case is the offline target for a digital shadow.  In operations its
``measured ORM`` must come from archived or dedicated low-amplitude data, not
from access to hidden parameters as done here.  No live machine access is
performed.
"""

from __future__ import annotations

import json
from dataclasses import dataclass

import numpy as np

from Interfaces.ATF2.InterfaceATF2_BT_RFTrack import InterfaceATF2_BT_RFTrack


RESPONSE_STEP_T_MM = 0.002
SVD_RCOND = 1e-3
MAX_CORRECTION_T_MM = 0.05


def _physical_orbit(machine: InterfaceATF2_BT_RFTrack) -> np.ndarray:
    bpms = machine.get_bpms()
    return np.r_[bpms["x"][0], bpms["y"][0]]


def _rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(values))))


@dataclass(frozen=True)
class PseudoMachineErrors:
    quadrupole_scales: dict[str, float]
    corrector_gain: np.ndarray
    corrector_offset_t_mm: np.ndarray
    bpm_gain_x: np.ndarray
    bpm_gain_y: np.ndarray
    bpm_offset_x_mm: np.ndarray
    bpm_offset_y_mm: np.ndarray
    bpm_roll_rad: np.ndarray


class BTPseudoMachine:
    """Hidden physical errors plus an explicit BPM measurement model."""

    def __init__(self, seed: int = 20260902):
        self.machine = InterfaceATF2_BT_RFTrack(reference_particle=True)
        self.corrector_names = self.machine.get_correctors()["names"]
        self.bpm_names = self.machine.get_bpms()["names"]
        rng = np.random.default_rng(seed)
        self.errors = PseudoMachineErrors(
            quadrupole_scales={"QD10T": 1.03, "QF11T": 0.98, "QF20T": 1.02},
            corrector_gain=1.0 + rng.normal(0.0, 0.02, len(self.corrector_names)),
            corrector_offset_t_mm=rng.normal(0.0, 0.002, len(self.corrector_names)),
            bpm_gain_x=1.0 + rng.normal(0.0, 0.015, len(self.bpm_names)),
            bpm_gain_y=1.0 + rng.normal(0.0, 0.015, len(self.bpm_names)),
            bpm_offset_x_mm=rng.normal(0.0, 0.030, len(self.bpm_names)),
            bpm_offset_y_mm=rng.normal(0.0, 0.030, len(self.bpm_names)),
            bpm_roll_rad=rng.normal(0.0, 1.0e-3, len(self.bpm_names)),
        )
        for name, scale in self.errors.quadrupole_scales.items():
            value = float(self.machine.get_quadrupoles(name)["bdes"][0])
            self.machine.set_quadrupoles(name, value * scale)
        self.commands = np.zeros(len(self.corrector_names), dtype=float)
        self.set_commands(self.commands)

    def set_commands(self, commands: np.ndarray) -> None:
        self.commands = np.asarray(commands, dtype=float).copy()
        actual = (
            self.errors.corrector_gain * self.commands
            + self.errors.corrector_offset_t_mm
        )
        self.machine.set_correctors(self.corrector_names, actual)

    def observe(self, *, calibrated: bool) -> np.ndarray:
        physical = _physical_orbit(self.machine)
        count = len(self.bpm_names)
        x, y = physical[:count], physical[count:]
        cosine, sine = np.cos(self.errors.bpm_roll_rad), np.sin(self.errors.bpm_roll_rad)
        measured_x = self.errors.bpm_gain_x * (cosine * x - sine * y) + self.errors.bpm_offset_x_mm
        measured_y = self.errors.bpm_gain_y * (sine * x + cosine * y) + self.errors.bpm_offset_y_mm
        if not calibrated:
            return np.r_[measured_x, measured_y]
        # In operational use these coefficients come from the calibration fit;
        # this synthetic audit uses the hidden values only as ground truth.
        rotated_x = (measured_x - self.errors.bpm_offset_x_mm) / self.errors.bpm_gain_x
        rotated_y = (measured_y - self.errors.bpm_offset_y_mm) / self.errors.bpm_gain_y
        return np.r_[cosine * rotated_x + sine * rotated_y, -sine * rotated_x + cosine * rotated_y]


def _measure_response(set_commands, observe, n_correctors: int) -> np.ndarray:
    origin = np.zeros(n_correctors, dtype=float)
    columns = []
    for index in range(n_correctors):
        plus = origin.copy()
        plus[index] += RESPONSE_STEP_T_MM
        set_commands(plus)
        orbit_plus = observe()
        minus = origin.copy()
        minus[index] -= RESPONSE_STEP_T_MM
        set_commands(minus)
        orbit_minus = observe()
        columns.append((orbit_plus - orbit_minus) / (2.0 * RESPONSE_STEP_T_MM))
    set_commands(origin)
    return np.column_stack(columns)


def _solve_and_apply(pseudo: BTPseudoMachine, response: np.ndarray, *, calibrated: bool) -> dict[str, float]:
    pseudo.set_commands(np.zeros(len(pseudo.corrector_names)))
    measured_before = pseudo.observe(calibrated=calibrated)
    physical_before = _physical_orbit(pseudo.machine)
    command = np.linalg.pinv(response, rcond=SVD_RCOND) @ (-measured_before)
    command = np.clip(command, -MAX_CORRECTION_T_MM, MAX_CORRECTION_T_MM)
    predicted = measured_before + response @ command
    pseudo.set_commands(command)
    return {
        "physical_before_rms_mm": _rms(physical_before),
        "physical_after_rms_mm": _rms(_physical_orbit(pseudo.machine)),
        "measurement_after_rms_mm": _rms(pseudo.observe(calibrated=calibrated)),
        "prediction_rms_mm": _rms(predicted),
        "max_abs_command_t_mm": float(np.max(np.abs(command))),
    }


def run_calibration_audit(seed: int = 20260902) -> dict[str, object]:
    """Compare nominal and calibrated response-matrix correction fairly."""
    nominal = InterfaceATF2_BT_RFTrack(reference_particle=True)
    pseudo = BTPseudoMachine(seed)
    nominal_response = _measure_response(
        lambda command: nominal.set_correctors(nominal.corrs, command),
        lambda: _physical_orbit(nominal),
        len(nominal.corrs),
    )
    pseudo_raw_response = _measure_response(
        pseudo.set_commands,
        lambda: pseudo.observe(calibrated=False),
        len(pseudo.corrector_names),
    )
    pseudo_calibrated_response = _measure_response(
        pseudo.set_commands,
        lambda: pseudo.observe(calibrated=True),
        len(pseudo.corrector_names),
    )
    raw_nominal = _solve_and_apply(pseudo, nominal_response, calibrated=False)
    calibrated_nominal = _solve_and_apply(pseudo, nominal_response, calibrated=True)
    calibrated_measured = _solve_and_apply(pseudo, pseudo_calibrated_response, calibrated=True)
    return {
        "simulation_only": True,
        "response_shape": list(nominal_response.shape),
        "response_rank": int(np.linalg.matrix_rank(nominal_response)),
        "relative_response_mismatch": {
            "raw_pseudo_vs_nominal": float(
                np.linalg.norm(pseudo_raw_response - nominal_response) / np.linalg.norm(pseudo_raw_response)
            ),
            "calibrated_pseudo_vs_nominal": float(
                np.linalg.norm(pseudo_calibrated_response - nominal_response)
                / np.linalg.norm(pseudo_calibrated_response)
            ),
        },
        "physical_orbit_rms_mm": {
            "nominal_orm_raw_bpm": raw_nominal,
            "nominal_orm_calibrated_bpm": calibrated_nominal,
            "measured_orm_calibrated_bpm": calibrated_measured,
        },
        "hidden_error_summary": {
            "quadrupole_scales": pseudo.errors.quadrupole_scales,
            "corrector_gain_rms_fraction": float(np.std(pseudo.errors.corrector_gain - 1.0)),
            "corrector_offset_rms_t_mm": float(np.std(pseudo.errors.corrector_offset_t_mm)),
            "bpm_offset_rms_mm": float(
                np.sqrt(
                    np.mean(
                        np.r_[pseudo.errors.bpm_offset_x_mm, pseudo.errors.bpm_offset_y_mm] ** 2
                    )
                )
            ),
            "bpm_roll_rms_mrad": float(np.std(pseudo.errors.bpm_roll_rad) * 1e3),
        },
        "interpretation": (
            "The measured-ORM case is a synthetic target for an offline digital shadow, "
            "not an estimate of current ATF performance."
        ),
    }


if __name__ == "__main__":
    print(json.dumps(run_calibration_audit(), indent=2))
