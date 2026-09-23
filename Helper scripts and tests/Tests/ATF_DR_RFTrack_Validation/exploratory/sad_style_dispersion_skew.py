"""SAD-style RF-Track virtual-machine tests for dispersion and skew correction.

This follows the measurement topology of ``coddispersion.n`` and ``skewcor.n``:
the virtual machine supplies off-momentum BPM data or two ZH-probe data sets;
the independent RF-Track model is first tune-fitted and then supplies the
response matrix and bounded-SVD actuator proposal.
"""

from __future__ import annotations

import json
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from digital_twin_random_errors import SEED, SIGMA_BPM_MM, _machine, _rms
from sad_style_cod import fit_model_to_measured_tunes, tunes
from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_correction import ATFDRRingCorrection
from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_lattice import build_atf_dr_lattice


# SAD coddispersion.n uses an energy offset of order 3.5e-3 for alpha=0.002.
DELTA = 3.5e-3
PROBE = 1e-4
MAX_ACTUATORS = 5
CORRECTION_GAIN = 0.35
MAX_STEERER_COMMAND = 0.01
MAX_SKEW_COMMAND = 2.0e-3


def _orbit(machine, guess=None, momentum=None):
    return machine.find_closed_orbit(
        initial_coordinates=guess, momentum_mev_c=momentum, max_iterations=30
    )


def _measure_dispersion(machine, rng, guess):
    """Virtual ``EtaCorrect`` measurement: (BPM(+d)-BPM(-d))/(2d)."""
    reference = _orbit(machine, guess)
    plus = _orbit(machine, reference.initial_coordinates, machine.momentum_mev_c * (1.0 + DELTA))
    minus = _orbit(machine, reference.initial_coordinates, machine.momentum_mev_c * (1.0 - DELTA))
    truth = (plus.bpm_positions - minus.bpm_positions) / (2.0 * DELTA)
    noise = rng.normal(0.0, SIGMA_BPM_MM, plus.bpm_positions.shape)
    noise -= rng.normal(0.0, SIGMA_BPM_MM, minus.bpm_positions.shape)
    return reference.initial_coordinates, reference.bpm_positions.copy(), truth, truth + noise / (2.0 * DELTA)


def _apply_steerers(lattice, names, plane, commands, gains):
    for name, command in zip(names, commands):
        item = lattice[name]
        item = item[0] if isinstance(item, list) else item
        strength = np.asarray(item.get_strength(), dtype=float)
        strength[plane] += gains[name] * command
        item.set_strength(*strength)


def _apply_skews(machine, names, commands):
    for name, command in zip(names, commands):
        machine.set_skew_strength(name, machine.get_skew_strength(name) + command)


def _choose_and_solve(matrix, residual, maximum, rcond=1e-4):
    selected = ATFDRRingCorrection._greedy_columns(
        matrix, -residual, MAX_ACTUATORS, rcond
    )
    values, _, rank = ATFDRRingCorrection._svd_solve(
        matrix[:, selected], -residual, rcond
    )
    commands = np.zeros(matrix.shape[1])
    commands[selected] = np.clip(CORRECTION_GAIN * values, -maximum, maximum)
    return commands, selected, rank


def _distributed(names, count):
    return tuple(names[index] for index in np.linspace(0, len(names) - 1, count, dtype=int))


def dispersion_case():
    machine_lattice, machine, gains, _, guess = _machine()
    rng = np.random.default_rng(SEED + 200)
    guess, orbit_before, truth_before, measured_before = _measure_dispersion(machine, rng, guess)
    measured_tunes = tunes(machine_lattice, guess)

    model_lattice = build_atf_dr_lattice()
    fit_model_to_measured_tunes(model_lattice, measured_tunes)
    model = ATFDRRingCorrection(model_lattice)
    target = model.measure_dispersion(relative_momentum_step=DELTA).values
    result = {"tunes": measured_tunes.tolist(), "planes": {}}

    # SAD uses all available steerers; five distributed candidates keep this
    # finite-difference validation short while retaining the same solve logic.
    for plane_name, plane in (("x", 0), ("y", 1)):
        names = _distributed(model.get_corrector_names(plane_name), 5)
        response = model.compute_dispersion_response(
            plane_name,
            corrector_names=names,
            relative_momentum_step=DELTA,
            perturbation=1e-5,
        )
        matrix = np.vstack((0.05 * response.dispersion_matrix, response.orbit_matrix))
        residual = np.concatenate((
            0.05 * (measured_before[:, plane] - target[:, plane]),
            orbit_before[:, plane],
        ))
        commands, selected, rank = _choose_and_solve(matrix, residual, MAX_STEERER_COMMAND)
        _apply_steerers(machine_lattice, names, plane, commands, gains)
        result["planes"][plane_name] = {
            "candidates": list(names),
            "selected": [names[index] for index in selected],
            "rank": rank,
            "max_command": float(np.max(np.abs(commands))),
            "dispersion_rms_before_mm": _rms(truth_before[:, plane] - target[:, plane]),
            "orbit_rms_before_mm": _rms(orbit_before[:, plane]),
        }

    _, orbit_after, truth_after, _ = _measure_dispersion(machine, rng, guess)
    for plane_name, plane in (("x", 0), ("y", 1)):
        result["planes"][plane_name].update({
            "dispersion_rms_after_mm": _rms(truth_after[:, plane] - target[:, plane]),
            "orbit_rms_after_mm": _rms(orbit_after[:, plane]),
        })
    return result


def _measure_coupling(machine, rng, guess, probes):
    """Virtual ``CorrectSkewCoupling`` readback from two horizontal probes."""
    baseline = _orbit(machine, guess)
    signal = []
    for name in probes:
        item = machine._single_element(name)
        strength = np.asarray(item.get_strength(), dtype=float)
        plus, minus = strength.copy(), strength.copy()
        plus[0] += machine.actuator_scale * PROBE
        minus[0] -= machine.actuator_scale * PROBE
        item.set_strength(*plus)
        positive = _orbit(machine, baseline.initial_coordinates)
        item.set_strength(*minus)
        negative = _orbit(machine, baseline.initial_coordinates)
        item.set_strength(*strength)
        truth = (positive.y - negative.y) / (2.0 * PROBE)
        noise = (rng.normal(0.0, SIGMA_BPM_MM, truth.shape) - rng.normal(0.0, SIGMA_BPM_MM, truth.shape)) / (2.0 * PROBE)
        signal.append((truth, truth + noise))
    return baseline.initial_coordinates, np.asarray([item[0] for item in signal]), np.asarray([item[1] for item in signal])


def skew_case():
    machine_lattice, machine, _, _, guess = _machine()
    rng = np.random.default_rng(SEED + 300)
    guess, truth_before, measured_before = _measure_coupling(machine, rng, guess, ("ZH1R", "ZH2R"))
    measured_tunes = tunes(machine_lattice, guess)
    model_lattice = build_atf_dr_lattice()
    fit_model_to_measured_tunes(model_lattice, measured_tunes)
    model = ATFDRRingCorrection(model_lattice)
    skews = _distributed(model.get_skew_corrector_names(), 5)
    response = model.compute_coupling_response(
        probe_corrector_names=("ZH1R", "ZH2R"),
        skew_corrector_names=skews,
        probe_perturbation=PROBE,
        skew_perturbation=1e-5,
    )
    commands, selected, rank = _choose_and_solve(
        response.matrix, measured_before.reshape(-1), MAX_SKEW_COMMAND, rcond=0.05
    )
    _apply_skews(machine, skews, commands)
    _, truth_after, _ = _measure_coupling(machine, rng, guess, ("ZH1R", "ZH2R"))
    return {
        "tunes": measured_tunes.tolist(),
        "probes": ["ZH1R", "ZH2R"],
        "candidates": list(skews),
        "selected": [skews[index] for index in selected],
        "rank": rank,
        "max_command": float(np.max(np.abs(commands))),
        "vertical_probe_response_rms_before_mm": _rms(truth_before),
        "vertical_probe_response_rms_after_mm": _rms(truth_after),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--case", choices=("all", "dispersion", "skew"), default="all",
        help="Run one observable family when a short standalone validation is needed.",
    )
    args = parser.parse_args()
    result = {
        "seed": SEED,
        "magnet_error_scale": 0.1,
        "method": "SAD-style measured observables + tune-fitted RF-Track response + bounded SVD",
        "test_scope": (
            "Five distributed actuator candidates per plane are used for a "
            "fast finite-difference validation; the RF-Track correction API "
            "also accepts the full SAD actuator lists."
        ),
    }
    if args.case in ("all", "dispersion"):
        result["dispersion"] = dispersion_case()
    if args.case in ("all", "skew"):
        result["skew"] = skew_case()
    analysis = Path(__file__).resolve().parents[5] / "analysis" / "DR-RFTrack"
    suffix = "" if args.case == "all" else f"_{args.case}"
    output = analysis / f"sad_style_dispersion_skew_result{suffix}.json"
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")

    count = 2 if args.case == "dispersion" else 1 if args.case == "skew" else 3
    figure, axes = plt.subplots(1, count, figsize=(3.7 * count, 3.5), constrained_layout=True)
    axes = np.atleast_1d(axes)
    if "dispersion" in result:
        for axis, plane in zip(axes[:2], ("x", "y")):
            item = result["dispersion"]["planes"][plane]
            values = np.array((item["dispersion_rms_before_mm"], item["dispersion_rms_after_mm"])) * 1e3
            bars = axis.bar(("before", "after"), values, color=("0.55", "tab:green"))
            axis.bar_label(bars, labels=[f"{value:.0f}" for value in values], padding=3)
            axis.set(title=f"{plane}: dispersion residual", ylabel="RMS [um / delta]")
            axis.grid(True, axis="y", alpha=0.25)
    if "skew" in result:
        axis = axes[-1]
        values = np.array((result["skew"]["vertical_probe_response_rms_before_mm"], result["skew"]["vertical_probe_response_rms_after_mm"])) * 1e3
        bars = axis.bar(("before", "after"), values, color=("0.55", "tab:green"))
        axis.bar_label(bars, labels=[f"{value:.0f}" for value in values], padding=3)
        axis.set(title="skew: vertical probe response", ylabel="RMS [um / control unit]")
        axis.grid(True, axis="y", alpha=0.25)
    title = {
        "all": "RF-Track virtual-machine SAD-style correction",
        "dispersion": "SAD-style dispersion correction",
        "skew": "SAD-style skew correction",
    }[args.case]
    figure.suptitle(title, fontsize=12)
    figure.savefig(analysis / f"sad_style_dispersion_skew_result{suffix}.png", dpi=180)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
