"""SAD-style model-assisted COD correction on an RF-Track virtual ATF DR.

The virtual machine exposes only BPM readbacks and tunes.  The RF-Track model
then fits the two QF families to those tunes, calculates its own response
matrix, chooses six correctors by forward selection and solves a bounded SVD
correction.  It does not use a response matrix measured on the virtual machine.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from digital_twin_random_errors import (
    SEED,
    _machine,
    _read,
    _rms,
)
from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_correction import ATFDRRingCorrection
from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_lattice import build_atf_dr_lattice
from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_twiss import (
    _linear_map,
    _periodic_twiss,
)


QF_FAMILIES = ("QF1R", "QF2R")
TUNE_FIT_STEP = 1e-4
TUNE_FIT_LIMIT = 0.01
MAX_CORRECTORS = 6
CORRECTION_GAIN = 0.5
MAX_COMMAND = 0.01


def tunes(lattice, initial_coordinates=None):
    correction = ATFDRRingCorrection(lattice)
    orbit = correction.find_closed_orbit(
        initial_coordinates=initial_coordinates, max_iterations=30
    )
    matrix = _linear_map(lattice, orbit.initial_coordinates, TUNE_FIT_STEP)
    return np.array(
        (
            _periodic_twiss(matrix[:2, :2])[-1],
            _periodic_twiss(matrix[2:, 2:])[-1],
        )
    )


def change_qf_family(lattice, family, relative_change):
    for quadrupole in lattice.get_quadrupoles():
        if quadrupole.get_name().startswith(family):
            quadrupole.set_strength(
                quadrupole.get_strength() * (1.0 + relative_change)
            )


def fit_model_to_measured_tunes(lattice, measured_tunes):
    """RF-Track equivalent of SAD ``FIT NX ... NY ...`` with QF families."""
    for _ in range(8):
        model_tunes = tunes(lattice)
        if np.max(np.abs(model_tunes - measured_tunes)) < 1e-6:
            break
        response = np.empty((2, len(QF_FAMILIES)))
        for column, family in enumerate(QF_FAMILIES):
            change_qf_family(lattice, family, TUNE_FIT_STEP)
            plus = tunes(lattice)
            change_qf_family(lattice, family, -2.0 * TUNE_FIT_STEP)
            minus = tunes(lattice)
            change_qf_family(lattice, family, TUNE_FIT_STEP)
            response[:, column] = (plus - minus) / (2.0 * TUNE_FIT_STEP)
        update = np.linalg.solve(response, measured_tunes - model_tunes)
        update = np.clip(update, -TUNE_FIT_LIMIT, TUNE_FIT_LIMIT)
        for family, value in zip(QF_FAMILIES, update):
            change_qf_family(lattice, family, value)
    return tunes(lattice)


def propose_sad_style_correction(response, readback):
    selected = ATFDRRingCorrection._greedy_columns(
        response.matrix, -readback, MAX_CORRECTORS, 1e-4
    )
    values, singular_values, rank = ATFDRRingCorrection._svd_solve(
        response.matrix[:, selected], -readback, 1e-4
    )
    commands = np.zeros(response.matrix.shape[1])
    commands[selected] = np.clip(
        CORRECTION_GAIN * values, -MAX_COMMAND, MAX_COMMAND
    )
    return commands, selected, rank, singular_values


def apply_to_virtual_machine(lattice, names, plane, commands, gains):
    for name, command in zip(names, commands):
        item = lattice[name]
        item = item[0] if isinstance(item, list) else item
        strength = np.asarray(item.get_strength(), dtype=float)
        strength[plane] += gains[name] * command
        item.set_strength(*strength)


def run_case(fit_tunes):
    machine_lattice, machine, gains, _, initial_guess = _machine()
    rng = np.random.default_rng(SEED + 100)
    guess, before, readback = _read(machine, rng, initial_guess)
    measured_tunes = tunes(machine_lattice, guess)

    model_lattice = build_atf_dr_lattice()
    model_tunes_before = tunes(model_lattice)
    model_tunes_after = (
        fit_model_to_measured_tunes(model_lattice, measured_tunes)
        if fit_tunes
        else model_tunes_before
    )
    model = ATFDRRingCorrection(model_lattice)
    detail = {}
    for name, index in (("x", 0), ("y", 1)):
        correctors = model.get_corrector_names(name)[::3]
        response = model.compute_orbit_response(
            name, corrector_names=correctors, perturbation=1e-5
        )
        commands, selected, rank, _ = propose_sad_style_correction(
            response, readback[:, index]
        )
        detail[name] = {
            "predicted_rms_mm": _rms(readback[:, index] + response.matrix @ commands),
            "selected_correctors": [correctors[item] for item in selected],
            "rank": rank,
            "max_command": float(np.max(np.abs(commands))),
        }
        apply_to_virtual_machine(
            machine_lattice, correctors, index, commands, gains
        )
    _, after, _ = _read(machine, rng, guess)
    return {
        "measured_tunes": measured_tunes.tolist(),
        "model_tunes_before_fit": model_tunes_before.tolist(),
        "model_tunes_after_fit": model_tunes_after.tolist(),
        "before_cod_rms_mm": {"x": _rms(before[:, 0]), "y": _rms(before[:, 1])},
        "after_cod_rms_mm": {"x": _rms(after[:, 0]), "y": _rms(after[:, 1])},
        "planes": detail,
    }


def main():
    result = {
        "seed": SEED,
        "method": "SAD-style: measured tune + model response + bounded SVD",
        "without_tune_fit": run_case(False),
        "with_tune_fit": run_case(True),
    }
    analysis = Path(__file__).resolve().parents[5] / "analysis" / "DR-RFTrack"
    output = analysis / "sad_style_tune_fit_cod_result.json"
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")

    figure, axes = plt.subplots(1, 2, figsize=(8.8, 3.7), constrained_layout=True)
    for plane, axis in (("x", axes[0]), ("y", axes[1])):
        values = np.array((
            result["with_tune_fit"]["before_cod_rms_mm"][plane],
            result["without_tune_fit"]["after_cod_rms_mm"][plane],
            result["with_tune_fit"]["after_cod_rms_mm"][plane],
        )) * 1e3
        bars = axis.bar(("before", "model only", "tune-fitted model"), values, color=("0.55", "tab:blue", "tab:green"))
        axis.bar_label(bars, labels=[f"{value:.0f}" for value in values], padding=3)
        axis.set(title=f"{plane}-plane COD", ylabel="COD RMS [um]")
        axis.tick_params(axis="x", rotation=15, labelsize=8)
        axis.grid(True, axis="y", alpha=0.25)
    figure.suptitle("RF-Track SAD-style COD correction", fontsize=14)
    figure.savefig(analysis / "sad_style_tune_fit_cod_result.png", dpi=180)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
