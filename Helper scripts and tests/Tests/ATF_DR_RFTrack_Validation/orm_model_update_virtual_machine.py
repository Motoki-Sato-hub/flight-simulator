"""Offline ORM model-update validation on an ATF DR RFTrack virtual machine.

This is the first, intentionally small digital-twin fit.  The virtual machine
has five hidden, physically interpretable deviations:

* QF1R and QF2R family-strength scales;
* a common BPM roll;
* horizontal and vertical corrector calibration gains.

Only a training set of four correctors is used to fit the model.  Four separate
correctors are held out and provide the predictive validation.  The script does
not use BPM offsets, alignment errors, individual magnet errors, or PV access.

Run from ``flight-simulator``:

    PYTHONPATH=. /home/motokisato/rftrack-env/bin/python \
      "Helper scripts and tests/Tests/ATF_DR_RFTrack_Validation/orm_model_update_virtual_machine.py"

The default ``ATF_DR_ORM_MODE=fit`` writes the ORM-fit and hold-out result.
Run ``cod``, ``dispersion``, and ``coupling`` separately to append the three
correction comparisons, then use ``plot`` to create their combined figure.
``all`` is convenient on an unrestricted workstation but intentionally takes
longer than the individual modes.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_correction import ATFDRRingCorrection
from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_lattice import build_atf_dr_lattice
from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_orm import (
    apply_quadrupole_individual_scale,
    bpm_roll_readout,
    capture_quadrupole_k1l,
    fit_linearised_orm,
    flatten_bpm_planes,
    measure_full_orbit_response,
)


PARAMETER_NAMES = (
    "QF1R relative K1L",
    "QF2R relative K1L",
    "common BPM roll [rad]",
    "horizontal corrector relative gain",
    "vertical corrector relative gain",
)
PARAMETER_SIGMA = np.asarray((0.01, 0.01, 2.0e-3, 0.05, 0.05))
FINITE_DIFFERENCE_STEP = np.asarray((2.0e-4, 2.0e-4, 1.0e-4, 1.0e-3, 1.0e-3))
# The two family errors are small enough for a one-step ORM linearisation;
# BPM roll and gain terms remain deliberately visible in the cross-plane ORM.
TRUE_PARAMETERS = np.asarray((4.0e-5, -3.0e-5, 7.0e-4, 2.5e-2, -1.8e-2))
ORM_PROBE_RAD = 1.0e-5
DISPERSION_DELTA = 2.0e-3
COUPLING_PROBE_RAD = 4.05e-5  # SAD: 270 urad x steerchange=0.15.


def _analysis_dir() -> Path:
    return Path(__file__).resolve().parents[4] / "analysis" / "DR-RFTrack"


def _corrector_sets(correction: ATFDRRingCorrection):
    """Use phase-separated training and hold-out x/y correctors."""
    x = correction.get_corrector_names("x")
    y = correction.get_corrector_names("y")
    train_names = (x[4], x[22], y[4], y[23])
    train_planes = ("x", "x", "y", "y")
    holdout_names = (x[12], x[36], y[12], y[35])
    holdout_planes = ("x", "x", "y", "y")
    return train_names, train_planes, holdout_names, holdout_planes


def _build_model(
    parameters: np.ndarray,
    *,
    individual_quadrupole_scale: dict[str, float] | None = None,
    bpm_roll_offsets_rad: np.ndarray | None = None,
    bpm_offsets_mm: np.ndarray | None = None,
    corrector_gain_offsets: dict[str, float] | None = None,
):
    """Build a nominal/trial model or a virtual machine with extra errors.

    The optional terms are deliberately not part of the five-parameter fit.
    They let seed scans test prediction under model-mismatch components.
    BPM offsets are returned for COD readout users but cancel in a noiseless
    central-difference ORM.
    """
    qf1, qf2, bpm_roll, gain_x, gain_y = np.asarray(parameters, dtype=float)
    lattice = build_atf_dr_lattice()
    reference = capture_quadrupole_k1l(lattice)
    component_scale = dict(individual_quadrupole_scale or {})
    for name in reference:
        family = name.split(".", maxsplit=1)[0]
        family_scale = {"QF1R": qf1, "QF2R": qf2}.get(family, 0.0)
        component_scale[name] = (
            (1.0 + family_scale) * (1.0 + component_scale.get(name, 0.0)) - 1.0
        )
    apply_quadrupole_individual_scale(lattice, reference, component_scale)
    correction = ATFDRRingCorrection(lattice)
    rolls = np.full(len(correction.bpm_names), bpm_roll)
    if bpm_roll_offsets_rad is not None:
        rolls += np.asarray(bpm_roll_offsets_rad, dtype=float)
    offsets = (
        np.zeros((len(correction.bpm_names), 2), dtype=float)
        if bpm_offsets_mm is None else np.asarray(bpm_offsets_mm, dtype=float)
    )
    gains = {
        name: (1.0 + (gain_x if plane == "x" else gain_y))
        * (1.0 + float((corrector_gain_offsets or {}).get(name, 0.0)))
        for name, plane in zip(
            correction.get_corrector_names("x") + correction.get_corrector_names("y"),
            ("x",) * len(correction.get_corrector_names("x"))
            + ("y",) * len(correction.get_corrector_names("y")),
        )
    }
    return correction, rolls, gains, offsets


def _measure(parameters, names, planes):
    correction, rolls, gains, _ = _build_model(parameters)
    return measure_full_orbit_response(
        correction,
        names,
        planes,
        bpm_rolls_rad=rolls,
        corrector_gains=gains,
        perturbation_rad=ORM_PROBE_RAD,
    )


def _sensitivity(train_names, train_planes, nominal_matrix):
    columns = []
    for index, step in enumerate(FINITE_DIFFERENCE_STEP):
        plus = np.zeros(len(PARAMETER_NAMES))
        minus = np.zeros(len(PARAMETER_NAMES))
        plus[index] = step
        minus[index] = -step
        derivative = (
            _measure(plus, train_names, train_planes).matrix
            - _measure(minus, train_names, train_planes).matrix
        ) / (2.0 * step)
        if derivative.shape != nominal_matrix.shape:
            raise RuntimeError("Finite-difference ORM shape changed")
        columns.append(derivative.reshape(-1))
    return np.column_stack(columns)


def _relative_error(measured: np.ndarray, predicted: np.ndarray) -> float:
    denominator = float(np.linalg.norm(measured))
    return float(np.linalg.norm(measured - predicted) / denominator) if denominator else 0.0


def _rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.asarray(values, dtype=float) ** 2)))


def _readout_orbit(machine: ATFDRRingCorrection, rolls: np.ndarray):
    closed = machine.find_closed_orbit()
    return closed.initial_coordinates, bpm_roll_readout(closed.bpm_positions, rolls)


def _add_corrector_kick(
    machine: ATFDRRingCorrection,
    name: str,
    plane: str,
    requested_kick_rad: float,
    physical_gains: dict[str, float],
) -> None:
    """Apply a requested kick through the hidden virtual-machine calibration."""
    element = machine._single_element(name)
    kick = machine._get_corrector_kick(element)
    kick[0 if plane == "x" else 1] += physical_gains[name] * requested_kick_rad
    machine._set_corrector_kick(element, kick)


def _least_squares_command(matrix: np.ndarray, residual: np.ndarray, maximum: float):
    command = np.linalg.lstsq(matrix, -residual, rcond=1.0e-3)[0]
    return np.clip(command, -maximum, maximum)


def _measure_vertical_dispersion(machine: ATFDRRingCorrection, rolls: np.ndarray):
    """Return vertical BPM orbit and vertical dispersion in the rolled readout."""
    reference = machine.find_closed_orbit()
    plus = machine.find_closed_orbit(
        momentum_mev_c=machine.momentum_mev_c * (1.0 + DISPERSION_DELTA),
        initial_coordinates=reference.initial_coordinates,
    )
    minus = machine.find_closed_orbit(
        momentum_mev_c=machine.momentum_mev_c * (1.0 - DISPERSION_DELTA),
        initial_coordinates=reference.initial_coordinates,
    )
    readout = bpm_roll_readout(reference.bpm_positions, rolls)
    dispersion = (
        bpm_roll_readout(plus.bpm_positions, rolls)
        - bpm_roll_readout(minus.bpm_positions, rolls)
    ) / (2.0 * DISPERSION_DELTA)
    return reference.initial_coordinates, readout[:, 1], dispersion[:, 1]


def _measure_vertical_dispersion_response(
    model: ATFDRRingCorrection,
    names: tuple[str, ...],
    rolls: np.ndarray,
    gains: dict[str, float],
):
    """Model response of [0.05 Dy, y COD] to vertical requested kicks."""
    _, base_y, base_dy = _measure_vertical_dispersion(model, rolls)
    matrix = np.empty((2 * len(base_y), len(names)))
    for column, name in enumerate(names):
        element = model._single_element(name)
        original = model._get_corrector_kick(element)
        physical_probe = ORM_PROBE_RAD * gains[name]
        plus, minus = original.copy(), original.copy()
        plus[1] += physical_probe
        minus[1] -= physical_probe
        try:
            model._set_corrector_kick(element, plus)
            _, plus_y, plus_dy = _measure_vertical_dispersion(model, rolls)
            model._set_corrector_kick(element, minus)
            _, minus_y, minus_dy = _measure_vertical_dispersion(model, rolls)
        finally:
            model._set_corrector_kick(element, original)
        matrix[:, column] = np.concatenate((
            0.05 * (plus_dy - minus_dy),
            plus_y - minus_y,
        )) / (2.0 * ORM_PROBE_RAD)
    baseline = np.concatenate((0.05 * base_dy, base_y))
    return matrix, baseline


def _measure_coupling_signal(
    machine: ATFDRRingCorrection,
    rolls: np.ndarray,
    gains: dict[str, float],
    probes: tuple[str, str],
) -> np.ndarray:
    """SAD-like vertical response to two horizontal steering probes."""
    baseline = machine.find_closed_orbit()
    signals = []
    for name in probes:
        element = machine._single_element(name)
        original = machine._get_corrector_kick(element)
        physical_probe = COUPLING_PROBE_RAD * gains[name]
        plus, minus = original.copy(), original.copy()
        plus[0] += physical_probe
        minus[0] -= physical_probe
        try:
            machine._set_corrector_kick(element, plus)
            positive = machine.find_closed_orbit(initial_coordinates=baseline.initial_coordinates)
            machine._set_corrector_kick(element, minus)
            negative = machine.find_closed_orbit(initial_coordinates=baseline.initial_coordinates)
        finally:
            machine._set_corrector_kick(element, original)
        response = (
            bpm_roll_readout(positive.bpm_positions, rolls)
            - bpm_roll_readout(negative.bpm_positions, rolls)
        ) / (2.0 * COUPLING_PROBE_RAD)
        signals.append(response[:, 1])
    return np.concatenate(signals)


def _measure_coupling_response(
    model: ATFDRRingCorrection,
    skew_names: tuple[str, ...],
    rolls: np.ndarray,
    gains: dict[str, float],
    probes: tuple[str, str],
) -> tuple[np.ndarray, np.ndarray]:
    """Model response of the two-probe coupling signal to skew K1L settings."""
    baseline = _measure_coupling_signal(model, rolls, gains, probes)
    matrix = np.empty((len(baseline), len(skew_names)))
    skew_probe = 1.0e-5
    for column, name in enumerate(skew_names):
        original = model.get_skew_strength(name)
        try:
            model.set_skew_strength(name, original + skew_probe)
            positive = _measure_coupling_signal(model, rolls, gains, probes)
            model.set_skew_strength(name, original - skew_probe)
            negative = _measure_coupling_signal(model, rolls, gains, probes)
        finally:
            model.set_skew_strength(name, original)
        matrix[:, column] = (positive - negative) / (2.0 * skew_probe)
    return matrix, baseline


def _correction_models(fitted: np.ndarray):
    """Return separate nominal and updated model states for a fair comparison."""
    return {
        "nominal": _build_model(np.zeros(len(PARAMETER_NAMES))),
        "updated": _build_model(fitted),
    }


def _virtual_machine():
    return _build_model(TRUE_PARAMETERS)


def _cod_correction_comparison(fitted: np.ndarray):
    base_model = _build_model(np.zeros(len(PARAMETER_NAMES)))[0]
    x = base_model.get_corrector_names("x")
    y = base_model.get_corrector_names("y")
    candidates = (x[4], x[12], x[22], x[36], y[4], y[12], y[23], y[35])
    planes = ("x", "x", "x", "x", "y", "y", "y", "y")
    results = {}
    for label, (model, rolls, gains, _) in _correction_models(fitted).items():
        machine, vm_rolls, vm_gains, _ = _virtual_machine()
        _add_corrector_kick(machine, x[45], "x", 8.0e-5, vm_gains)
        _add_corrector_kick(machine, y[45], "y", -7.0e-5, vm_gains)
        _, before = _readout_orbit(machine, vm_rolls)
        response = measure_full_orbit_response(
            model, candidates, planes, bpm_rolls_rad=rolls,
            corrector_gains=gains, perturbation_rad=ORM_PROBE_RAD,
        )
        command = _least_squares_command(response.matrix, flatten_bpm_planes(before), 1.0e-3)
        for name, plane, value in zip(candidates, planes, command):
            _add_corrector_kick(machine, name, plane, float(value), vm_gains)
        _, after = _readout_orbit(machine, vm_rolls)
        results[label] = {
            "before_rms_mm": _rms(before),
            "after_rms_mm": _rms(after),
            "maximum_requested_kick_urad": float(np.max(np.abs(command)) * 1.0e6),
        }
    return results


def _dispersion_correction_comparison(fitted: np.ndarray):
    base_model = _build_model(np.zeros(len(PARAMETER_NAMES)))[0]
    y = base_model.get_corrector_names("y")
    candidates = (y[4], y[12], y[23], y[35])
    results = {}
    for label, (model, rolls, gains, _) in _correction_models(fitted).items():
        machine, vm_rolls, vm_gains, _ = _virtual_machine()
        _add_corrector_kick(machine, y[45], "y", -8.0e-5, vm_gains)
        _, before_y, before_dy = _measure_vertical_dispersion(machine, vm_rolls)
        response, _ = _measure_vertical_dispersion_response(model, candidates, rolls, gains)
        command = _least_squares_command(
            response, np.concatenate((0.05 * before_dy, before_y)), 1.0e-3
        )
        for name, value in zip(candidates, command):
            _add_corrector_kick(machine, name, "y", float(value), vm_gains)
        _, after_y, after_dy = _measure_vertical_dispersion(machine, vm_rolls)
        results[label] = {
            "before_vertical_orbit_rms_mm": _rms(before_y),
            "after_vertical_orbit_rms_mm": _rms(after_y),
            "before_vertical_dispersion_rms_mm_per_delta": _rms(before_dy),
            "after_vertical_dispersion_rms_mm_per_delta": _rms(after_dy),
            "maximum_requested_kick_urad": float(np.max(np.abs(command)) * 1.0e6),
        }
    return results


def _coupling_correction_comparison(fitted: np.ndarray):
    base_model = _build_model(np.zeros(len(PARAMETER_NAMES)))[0]
    probes = ("ZH42R", "ZH44R")
    all_skews = base_model.get_skew_corrector_names()
    sf = tuple(name for name in all_skews if name.startswith("SF1R."))
    sd_fault = next(name for name in all_skews if name.startswith("SD1R."))
    # Four phase-distributed SF1R knobs keep this isolated ORM-update validation
    # quick; the SAD-operation mapping itself retains all 34 SF1R channels.
    candidates = tuple(sf[index] for index in np.linspace(0, len(sf) - 1, 4, dtype=int))
    results = {}
    for label, (model, rolls, gains, _) in _correction_models(fitted).items():
        machine, vm_rolls, vm_gains, _ = _virtual_machine()
        machine.set_skew_strength(sd_fault, 1.0e-3)
        before = _measure_coupling_signal(machine, vm_rolls, vm_gains, probes)
        response, _ = _measure_coupling_response(model, candidates, rolls, gains, probes)
        command = _least_squares_command(response, before, 0.034)
        for name, value in zip(candidates, command):
            machine.set_skew_strength(name, machine.get_skew_strength(name) + float(value))
        after = _measure_coupling_signal(machine, vm_rolls, vm_gains, probes)
        results[label] = {
            "before_coupling_rms_mm_per_rad": _rms(before),
            "after_coupling_rms_mm_per_rad": _rms(after),
            "maximum_requested_skew_k1l_m_inv": float(np.max(np.abs(command))),
        }
    return results


def _write_correction_figure(results: dict, analysis: Path) -> Path:
    """Render the three independently evaluated correction comparisons."""
    comparison = results["correction_comparison"]
    required = {"cod", "vertical_cod_dispersion", "coupling"}
    if set(comparison) != required:
        raise ValueError("All three correction cases are required for the comparison figure")
    figure, axes = plt.subplots(1, 3, figsize=(11.5, 3.4), constrained_layout=True)
    panels = (
        ("COD", "cod", "before_rms_mm", "after_rms_mm", "BPM orbit RMS [mm]"),
        ("vertical COD + dispersion", "vertical_cod_dispersion",
         "before_vertical_dispersion_rms_mm_per_delta",
         "after_vertical_dispersion_rms_mm_per_delta", "vertical dispersion RMS [mm/delta]"),
        ("coupling", "coupling", "before_coupling_rms_mm_per_rad",
         "after_coupling_rms_mm_per_rad", "coupling signal RMS [mm/rad]"),
    )
    for axis, (title, key, before_key, after_key, ylabel) in zip(axes, panels):
        case = comparison[key]
        values = [case["nominal"][before_key], case["nominal"][after_key],
                  case["updated"][after_key]]
        bars = axis.bar(("before", "nominal\nmodel", "updated\nmodel"), values,
                        color=("0.55", "tab:orange", "tab:blue"))
        axis.bar_label(bars, labels=[f"{value:.2e}" for value in values], padding=3, fontsize=8)
        axis.set(title=title, ylabel=ylabel)
        axis.tick_params(axis="x", labelsize=8)
        axis.grid(axis="y", alpha=0.3)
    path = analysis / "atf_dr_rftrack_orm_updated_model_correction_comparison.png"
    figure.savefig(path, dpi=180)
    return path


def main():
    mode = os.environ.get("ATF_DR_ORM_MODE", "fit").lower()
    modes = {"fit", "cod", "dispersion", "coupling", "plot", "all"}
    if mode not in modes:
        raise ValueError(f"ATF_DR_ORM_MODE must be one of {sorted(modes)}")
    analysis = _analysis_dir()
    result_path = analysis / "atf_dr_rftrack_orm_model_update_virtual_machine.json"
    case_functions = {
        "cod": _cod_correction_comparison,
        "dispersion": _dispersion_correction_comparison,
        "coupling": _coupling_correction_comparison,
    }
    if mode in case_functions or mode == "plot":
        results = json.loads(result_path.read_text(encoding="utf-8"))
        if mode == "plot":
            print("saved:", _write_correction_figure(results, analysis))
        else:
            key = "vertical_cod_dispersion" if mode == "dispersion" else mode
            results.setdefault("correction_comparison", {})[key] = case_functions[mode](
                np.asarray(results["fitted_parameters"], dtype=float)
            )
            result_path.write_text(json.dumps(results, indent=2), encoding="utf-8")
            print(json.dumps(results["correction_comparison"], indent=2))
            print("saved:", result_path)
        return

    nominal = _build_model(np.zeros(len(PARAMETER_NAMES)))[0]
    train_names, train_planes, holdout_names, holdout_planes = _corrector_sets(nominal)

    model_train = _measure(np.zeros(len(PARAMETER_NAMES)), train_names, train_planes)
    measured_train = _measure(TRUE_PARAMETERS, train_names, train_planes)
    sensitivity = _sensitivity(train_names, train_planes, model_train.matrix)
    fit = fit_linearised_orm(
        measured_train.matrix,
        model_train.matrix,
        sensitivity,
        PARAMETER_SIGMA,
        PARAMETER_NAMES,
    )

    model_holdout = _measure(np.zeros(len(PARAMETER_NAMES)), holdout_names, holdout_planes)
    measured_holdout = _measure(TRUE_PARAMETERS, holdout_names, holdout_planes)
    updated_holdout = _measure(fit.parameter_delta, holdout_names, holdout_planes)

    holdout_before = _relative_error(measured_holdout.matrix, model_holdout.matrix)
    holdout_after = _relative_error(measured_holdout.matrix, updated_holdout.matrix)
    correction_comparison = (
        {
            "cod": _cod_correction_comparison(fit.parameter_delta),
            "vertical_cod_dispersion": _dispersion_correction_comparison(fit.parameter_delta),
            "coupling": _coupling_correction_comparison(fit.parameter_delta),
        }
        if mode == "all" else {}
    )
    results = {
        "scope": "Offline RFTrack virtual-machine ORM model update; no PV access.",
        "training_correctors": list(train_names),
        "holdout_correctors": list(holdout_names),
        "parameter_names": list(PARAMETER_NAMES),
        "true_parameters": TRUE_PARAMETERS.tolist(),
        "fitted_parameters": fit.parameter_delta.tolist(),
        "parameter_prior_sigma": PARAMETER_SIGMA.tolist(),
        "fit": {
            "rank": fit.rank,
            "singular_values": fit.singular_values.tolist(),
            "training_residual_norm_before": fit.residual_norm_before,
            "training_residual_norm_after": fit.residual_norm_after,
            "training_residual_reduction": fit.residual_reduction,
        },
        "heldout_orm_relative_error": {
            "nominal_model": holdout_before,
            "updated_model": holdout_after,
            "reduction": 1.0 - holdout_after / holdout_before,
        },
        "orm_shape": list(measured_train.shape),
        "orm_probe_rad": ORM_PROBE_RAD,
        "correction_comparison": correction_comparison,
    }

    result_path.write_text(json.dumps(results, indent=2), encoding="utf-8")

    figure, axes = plt.subplots(1, 2, figsize=(10.2, 3.6), constrained_layout=True)
    indices = np.arange(len(PARAMETER_NAMES))
    # Display each mixed-unit parameter in a readable physical unit.
    display_scale = np.asarray((1.0e4, 1.0e4, 1.0e3, 1.0e2, 1.0e2))
    axes[0].bar(indices - 0.18, TRUE_PARAMETERS * display_scale, width=0.36, label="hidden truth")
    axes[0].bar(indices + 0.18, fit.parameter_delta * display_scale, width=0.36, label="ORM fit")
    axes[0].set(xticks=indices,
                xticklabels=("QF1R\n[1e-4]", "QF2R\n[1e-4]", "BPM roll\n[mrad]",
                            "x gain\n[%]", "y gain\n[%]"),
                ylabel="displayed parameter unit", title="hidden errors recovered from training ORM")
    axes[0].tick_params(axis="x", rotation=20)
    axes[0].grid(axis="y", alpha=0.3); axes[0].legend(fontsize=8)
    bars = axes[1].bar(("nominal", "updated"), (holdout_before, holdout_after),
                       color=("tab:orange", "tab:blue"))
    axes[1].bar_label(bars, labels=(f"{holdout_before:.3%}", f"{holdout_after:.3%}"), padding=3)
    axes[1].set(ylabel="relative ORM error", title="held-out corrector prediction")
    axes[1].grid(axis="y", alpha=0.3)
    figure_path = analysis / "atf_dr_rftrack_orm_model_update_virtual_machine.png"
    figure.savefig(figure_path, dpi=180)

    print(json.dumps(results, indent=2))
    print("saved:", result_path)
    print("saved:", figure_path)
    if mode == "all":
        print("saved:", _write_correction_figure(results, analysis))


if __name__ == "__main__":
    main()
