"""Nominal-twin versus perturbed-virtual-machine COD validation.

Run from the flight-simulator root:

    MPLBACKEND=Agg PYTHONPATH=. /home/motokisato/rftrack-env/bin/python \
      "Helper scripts and tests/Tests/ATF_DR_RFTrack_Validation/exploratory/digital_twin_random_errors.py"

The nominal lattice supplies the response matrix.  A separately built lattice
is the virtual machine and receives hidden, independently sampled errors.
Only its 1 um-rms noisy BPM readbacks are given to the nominal correction twin.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import RF_Track as rft

from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_correction import ATFDRRingCorrection
from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_lattice import build_atf_dr_lattice


SEED = int(os.environ.get("ATF_DR_ERROR_SEED", "20260913"))
MAGNET_ERROR_SCALE = float(os.environ.get("ATF_DR_MAGNET_ERROR_SCALE", "1.0"))
if MAGNET_ERROR_SCALE <= 0:
    raise ValueError("ATF_DR_MAGNET_ERROR_SCALE must be positive")
SIGMA_ALIGNMENT_M = 100e-6 * MAGNET_ERROR_SCALE
SIGMA_STRENGTH = 2e-3 * MAGNET_ERROR_SCALE
SIGMA_ROLL_RAD = 100e-6 * MAGNET_ERROR_SCALE
SIGMA_BPM_MM = 1e-3


def _one(lattice, name):
    item = lattice[name]
    return item[0] if isinstance(item, list) else item


def _rms(values):
    return float(np.sqrt(np.mean(np.asarray(values, dtype=float) ** 2)))


def _apply_errors(lattice, rng, scale=1.0):
    """Apply specified rms errors. BPMs receive random readout noise only."""
    corrector_gain = {}
    counts = {"quadrupoles": 0, "sextupoles": 0, "bends": 0, "correctors": 0, "bpms": 0}
    for item in lattice["*"]:
        dx, dy = rng.normal(0.0, SIGMA_ALIGNMENT_M * scale, 2)
        roll = rng.normal(0.0, SIGMA_ROLL_RAD * scale)
        dk = rng.normal(0.0, SIGMA_STRENGTH * scale)
        if isinstance(item, rft.Quadrupole):
            item.set_offsets(dx, dy, 0.0, roll, 0.0, 0.0, "center")
            item.set_strength(item.get_strength() * (1.0 + dk))
            counts["quadrupoles"] += 1
        elif isinstance(item, rft.Sextupole):
            item.set_offsets(dx, dy, 0.0, roll, 0.0, 0.0, "center")
            item.set_strength(item.get_strength() * (1.0 + dk))
            counts["sextupoles"] += 1
        elif isinstance(item, rft.SBend):
            item.set_offsets(dx, dy, 0.0, roll, 0.0, 0.0, "center")
            item.set_angle(item.get_angle() * (1.0 + dk))
            item.set_K1(item.get_K1() * (1.0 + dk))
            counts["bends"] += 1
        elif isinstance(item, rft.Corrector):
            item.set_offsets(dx, dy, 0.0, roll, 0.0, 0.0, "center")
            corrector_gain[item.get_name()] = 1.0 + dk
            counts["correctors"] += 1
        elif isinstance(item, rft.Bpm):
            counts["bpms"] += 1
    return corrector_gain, counts


def _read(machine, rng, guess=None):
    orbit = machine.find_closed_orbit(initial_coordinates=guess, max_iterations=30)
    truth = orbit.bpm_positions.copy()
    return orbit.initial_coordinates, truth, truth + rng.normal(0.0, SIGMA_BPM_MM, truth.shape)


def _solve(response, readback, *, rcond, gain, maximum):
    u, singular, vt = np.linalg.svd(response, full_matrices=False)
    keep = singular > rcond * singular[0]
    delta = vt[keep].T @ ((-u[:, keep].T @ readback) / singular[keep])
    return np.clip(gain * delta, -maximum, maximum), int(np.count_nonzero(keep))


def _apply(lattice, names, plane, delta, corrector_gain):
    for name, value in zip(names, delta):
        item = _one(lattice, name)
        strength = np.asarray(item.get_strength(), dtype=float)
        strength[plane] += corrector_gain[name] * value
        item.set_strength(*strength)


def _machine_response(machine, lattice, names, plane, gains, rng, guess, probe=1e-3):
    """Independent micro-kick response, included only as a reference case."""
    base = machine.find_closed_orbit(initial_coordinates=guess, max_iterations=30)
    matrix = np.empty((len(base.bpm_names), len(names)))
    for column, name in enumerate(names):
        item = _one(lattice, name)
        initial = np.asarray(item.get_strength(), dtype=float)
        plus, minus = initial.copy(), initial.copy()
        plus[plane] += gains[name] * probe
        minus[plane] -= gains[name] * probe
        item.set_strength(*plus)
        positive = machine.find_closed_orbit(initial_coordinates=base.initial_coordinates, max_iterations=30)
        item.set_strength(*minus)
        negative = machine.find_closed_orbit(initial_coordinates=base.initial_coordinates, max_iterations=30)
        item.set_strength(*initial)
        noise = rng.normal(0.0, SIGMA_BPM_MM, 2 * len(base.bpm_names))
        matrix[:, column] = (positive.bpm_positions[:, plane] + noise[:len(base.bpm_names)] - negative.bpm_positions[:, plane] - noise[len(base.bpm_names):]) / (2.0 * probe)
    return matrix


def _machine():
    """Build the final virtual machine and obtain its orbit by continuation."""
    guess = np.zeros(4)
    for scale in np.linspace(0.0, 1.0, 11):
        lattice = build_atf_dr_lattice()
        gains, counts = _apply_errors(lattice, np.random.default_rng(SEED), scale)
        machine = ATFDRRingCorrection(lattice)
        guess = machine.find_closed_orbit(
            initial_coordinates=guess, max_iterations=30
        ).initial_coordinates
    return lattice, machine, gains, counts, guess


def main():
    twin = ATFDRRingCorrection(build_atf_dr_lattice())
    names = {plane: twin.get_corrector_names(plane)[::3] for plane in ("x", "y")}
    response = {plane: twin.compute_orbit_response(plane, corrector_names=names[plane]) for plane in ("x", "y")}

    lattice, machine, gains, counts, initial_guess = _machine()
    rng = np.random.default_rng(SEED + 1)
    guess, before, readback = _read(machine, rng, initial_guess)
    model = {}
    for plane_name, plane_index in (("x", 0), ("y", 1)):
        delta, rank = _solve(response[plane_name].matrix, readback[:, plane_index], rcond=2e-2, gain=0.20, maximum=0.01)
        model[plane_name] = {"predicted_rms_mm": _rms(readback[:, plane_index] + response[plane_name].matrix @ delta), "rank": rank, "max_command": float(np.max(np.abs(delta)))}
        _apply(lattice, names[plane_name], plane_index, delta, gains)
    _, model_after, _ = _read(machine, rng, guess)

    lattice, machine, gains, _, initial_guess = _machine()
    rng = np.random.default_rng(SEED + 2)
    guess, reference_before, reference_readback = _read(machine, rng, initial_guess)
    reference = {}
    for plane_name, plane_index in (("x", 0), ("y", 1)):
        measured = _machine_response(machine, lattice, names[plane_name], plane_index, gains, rng, guess)
        delta, rank = _solve(measured, reference_readback[:, plane_index], rcond=1e-2, gain=0.20, maximum=0.02)
        reference[plane_name] = {"rank": rank, "response_cosine": float(np.vdot(response[plane_name].matrix, measured) / (np.linalg.norm(response[plane_name].matrix) * np.linalg.norm(measured)))}
        _apply(lattice, names[plane_name], plane_index, delta, gains)
    _, reference_after, _ = _read(machine, rng, guess)

    result = {"seed": SEED, "magnet_error_scale": MAGNET_ERROR_SCALE, "error_rms": {"dx_dy_um": SIGMA_ALIGNMENT_M * 1e6, "dk_over_k_percent": SIGMA_STRENGTH * 100.0, "roll_urad": SIGMA_ROLL_RAD * 1e6, "bpm_noise_um": 1.0}, "component_counts": counts, "correctors_used": {plane: len(value) for plane, value in names.items()}, "twin_only_cod_rms_mm": {"before": {"x": _rms(before[:, 0]), "y": _rms(before[:, 1])}, "after": {"x": _rms(model_after[:, 0]), "y": _rms(model_after[:, 1])}, "details": model}, "measured_response_reference_cod_rms_mm": {"before": {"x": _rms(reference_before[:, 0]), "y": _rms(reference_before[:, 1])}, "after": {"x": _rms(reference_after[:, 0]), "y": _rms(reference_after[:, 1])}, "details": reference}}
    analysis = Path(__file__).resolve().parents[5] / "analysis" / "DR-RFTrack"
    suffix = "" if MAGNET_ERROR_SCALE == 1.0 else f"_magnet_error_x{MAGNET_ERROR_SCALE:g}".replace(".", "p")
    output = analysis / f"digital_twin_random_error_cod_result{suffix}.json"
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    labels = ("before", "twin prediction", "twin actual", "measured-R actual")
    figure, axes = plt.subplots(1, 2, figsize=(10.5, 4.4), constrained_layout=True)
    for plane_name, plane_index, axis, color in (("x", 0, axes[0], "tab:blue"), ("y", 1, axes[1], "tab:orange")):
        values = np.array((_rms(before[:, plane_index]), model[plane_name]["predicted_rms_mm"], _rms(model_after[:, plane_index]), _rms(reference_after[:, plane_index]))) * 1e3
        bars = axis.bar(labels, values, color=("0.55", "0.72", color, "tab:green"))
        axis.bar_label(bars, labels=[f"{value:.0f}" for value in values], padding=3, fontsize=9)
        axis.set_yscale("log")
        axis.set(title=f"{plane_name}-plane COD", ylabel="noise-free COD RMS [um]")
        axis.tick_params(axis="x", rotation=18, labelsize=8)
        axis.grid(True, axis="y", which="both", alpha=0.25)
    figure.suptitle(f"ATF DR: magnet-error scale {MAGNET_ERROR_SCALE:g}, BPM noise 1 um", fontsize=14)
    figure.savefig(analysis / f"digital_twin_random_error_cod_result{suffix}.png", dpi=180)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
