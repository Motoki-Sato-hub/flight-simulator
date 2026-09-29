"""Seed and error-scale scan for the ATF DR offline ORM model update.

The scan draws the same five fit parameters used by
``orm_model_update_virtual_machine.py`` from zero-mean Gaussian distributions:
QF1R/QF2R family K1L scales, a common BPM roll, and x/y corrector gains.
For every draw it fits one horizontal and one vertical training corrector and
scores the nominal and updated model on a separate horizontal/vertical pair.
The shorter scan layout makes many independent seeds practical; the companion
single-case study retains four training and four held-out correctors.

``matched`` draws only the five fitted parameters.  ``extended-noiseless`` and
``extended`` additionally draw independent quadrupole-strength errors, BPM
offsets and per-BPM rolls, and per-corrector gain errors.  ``extended`` also
adds independent 1-um BPM readout noise to every central-difference ORM
sample.  This deliberately evaluates
ORM prediction robustness, not the truth of an individual fitted parameter;
there is no live PV access.

Examples, from ``flight-simulator``:

    PYTHONPATH=. /home/motokisato/rftrack-env/bin/python \
      "Helper scripts and tests/Tests/ATF_DR_RFTrack_Validation/orm_model_update_seed_scan.py" \
      --seeds 2003,2004,2005 --scales 0.5

Run individual scales separately when operating under a short job time limit,
then use ``--summary-only`` to generate the combined plot.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_orm import (
    fit_linearised_orm,
    measure_full_orbit_response,
)


SCRIPT_PATH = Path(__file__).with_name("orm_model_update_virtual_machine.py")
ANALYSIS = Path(__file__).resolve().parents[4] / "analysis" / "DR-RFTrack"
RESULT_PATH = ANALYSIS / "atf_dr_rftrack_orm_model_update_fast_seed_scan.json"
FIGURE_PATH = ANALYSIS / "atf_dr_rftrack_orm_model_update_fast_seed_scan.png"
EXTENDED_RESULT_PATH = ANALYSIS / "atf_dr_rftrack_orm_model_update_extended_seed_scan.json"
EXTENDED_FIGURE_PATH = ANALYSIS / "atf_dr_rftrack_orm_model_update_extended_seed_scan.png"
EXTENDED_NOISELESS_RESULT_PATH = ANALYSIS / "atf_dr_rftrack_orm_model_update_extended_noiseless_seed_scan.json"
EXTENDED_NOISELESS_FIGURE_PATH = ANALYSIS / "atf_dr_rftrack_orm_model_update_extended_noiseless_seed_scan.png"

# Scale-dependent unmodelled virtual-machine errors.  The individual
# quadrupole term is relative K1L; the BPM offset is in RFTrack's mm readout
# coordinates.  These terms are deliberately absent from the five-parameter
# fitted model.
INDIVIDUAL_QUAD_SIGMA = 1.0e-4
INDIVIDUAL_BPM_ROLL_SIGMA_RAD = 3.0e-4
INDIVIDUAL_CORRECTOR_GAIN_SIGMA = 1.0e-2
BPM_OFFSET_SIGMA_MM = 2.0e-2


def _load_case_module():
    specification = importlib.util.spec_from_file_location("atf_dr_orm_case", SCRIPT_PATH)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"Could not import {SCRIPT_PATH}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def _parse_csv(value: str, converter):
    return tuple(converter(item) for item in value.split(",") if item.strip())


def _relative_error(measured: np.ndarray, model: np.ndarray) -> float:
    norm = float(np.linalg.norm(measured))
    return float(np.linalg.norm(measured - model) / norm) if norm else 0.0


def _std_or_zero(values) -> float:
    values = list(values)
    return float(np.std(values)) if values else 0.0


def _prepare(module):
    zero = np.zeros(len(module.PARAMETER_NAMES))
    nominal = module._build_model(zero)[0]
    all_train_names, all_train_planes, all_holdout_names, all_holdout_planes = module._corrector_sets(nominal)
    # Fast scan: retain one phase-separated horizontal and vertical channel in
    # each split.  The model parameter count remains five, with 392 ORM values
    # in the training data and 392 independent hold-out values.
    train_names = (all_train_names[0], all_train_names[2])
    train_planes = (all_train_planes[0], all_train_planes[2])
    holdout_names = (all_holdout_names[0], all_holdout_names[2])
    holdout_planes = (all_holdout_planes[0], all_holdout_planes[2])
    model_train = module._measure(zero, train_names, train_planes)
    sensitivity = module._sensitivity(train_names, train_planes, model_train.matrix)
    model_holdout = module._measure(zero, holdout_names, holdout_planes)
    return train_names, train_planes, holdout_names, holdout_planes, model_train, sensitivity, model_holdout


def _noise_orm(matrix: np.ndarray, rng: np.random.Generator, bpm_noise_um: float, probe_rad: float):
    """Add raw-BPM noise propagated through a central ORM difference."""
    if bpm_noise_um <= 0.0:
        return matrix
    sigma_mm_per_rad = bpm_noise_um * 1.0e-3 / (np.sqrt(2.0) * probe_rad)
    return matrix + rng.normal(0.0, sigma_mm_per_rad, size=matrix.shape)


def _extended_errors(module, rng: np.random.Generator, scale: float):
    """Draw error sources not represented by the five fitted parameters."""
    correction, _, _, _ = module._build_model(np.zeros(len(module.PARAMETER_NAMES)))
    quadrupoles = tuple(element.get_name() for element in correction.lattice.get_quadrupoles())
    correctors = correction.get_corrector_names("x") + correction.get_corrector_names("y")
    return {
        "individual_quadrupole_scale": {
            name: float(rng.normal(0.0, INDIVIDUAL_QUAD_SIGMA * scale))
            for name in quadrupoles
        },
        "bpm_roll_offsets_rad": rng.normal(
            0.0, INDIVIDUAL_BPM_ROLL_SIGMA_RAD * scale, len(correction.bpm_names)
        ),
        "bpm_offsets_mm": rng.normal(
            0.0, BPM_OFFSET_SIGMA_MM * scale, (len(correction.bpm_names), 2)
        ),
        "corrector_gain_offsets": {
            name: float(rng.normal(0.0, INDIVIDUAL_CORRECTOR_GAIN_SIGMA * scale))
            for name in correctors
        },
    }


def _measure_virtual_machine(module, parameters, names, planes, extra):
    correction, rolls, gains, offsets = module._build_model(parameters, **extra)
    return measure_full_orbit_response(
        correction,
        names,
        planes,
        bpm_rolls_rad=rolls,
        bpm_offsets_mm=offsets,
        corrector_gains=gains,
        perturbation_rad=module.ORM_PROBE_RAD,
    )


def _case(module, prepared, seed: int, scale: float, scenario: str, bpm_noise_um: float):
    (
        train_names, train_planes, holdout_names, holdout_planes,
        model_train, sensitivity, model_holdout,
    ) = prepared
    # These are the 1-sigma values used for the single-case virtual machine.
    error_sigma = np.asarray((4.0e-5, 3.0e-5, 7.0e-4, 2.5e-2, 1.8e-2))
    rng = np.random.default_rng(seed)
    true = rng.normal(0.0, error_sigma * scale)
    extra = _extended_errors(module, rng, scale) if scenario != "matched" else {}
    measured_train = _measure_virtual_machine(module, true, train_names, train_planes, extra)
    train_matrix = _noise_orm(
        measured_train.matrix, rng, bpm_noise_um if scenario == "extended" else 0.0,
        module.ORM_PROBE_RAD,
    )
    fit = fit_linearised_orm(
        train_matrix,
        model_train.matrix,
        sensitivity,
        module.PARAMETER_SIGMA,
        module.PARAMETER_NAMES,
    )
    measured_holdout = _measure_virtual_machine(module, true, holdout_names, holdout_planes, extra)
    holdout_matrix = _noise_orm(
        measured_holdout.matrix, rng, bpm_noise_um if scenario == "extended" else 0.0,
        module.ORM_PROBE_RAD,
    )
    updated_holdout = module._measure(fit.parameter_delta, holdout_names, holdout_planes)
    nominal_error = _relative_error(holdout_matrix, model_holdout.matrix)
    updated_error = _relative_error(holdout_matrix, updated_holdout.matrix)
    return {
        "seed": seed,
        "scenario": scenario,
        "error_scale": scale,
        "true_parameters": true.tolist(),
        "fitted_parameters": fit.parameter_delta.tolist(),
        "fit_rank": fit.rank,
        "training_residual_reduction": fit.residual_reduction,
        "heldout_orm_relative_error_nominal": nominal_error,
        "heldout_orm_relative_error_updated": updated_error,
        "heldout_orm_improvement": 1.0 - updated_error / nominal_error if nominal_error else 0.0,
        "unmodelled_error_rms": {
            "individual_quadrupole_relative_k1l": _std_or_zero(extra.get("individual_quadrupole_scale", {}).values()),
            "individual_bpm_roll_rad": _std_or_zero(extra.get("bpm_roll_offsets_rad", ())),
            "individual_corrector_gain": _std_or_zero(extra.get("corrector_gain_offsets", {}).values()),
        },
    }


def _paths(scenario: str):
    if scenario == "extended":
        return EXTENDED_RESULT_PATH, EXTENDED_FIGURE_PATH
    if scenario == "extended-noiseless":
        return EXTENDED_NOISELESS_RESULT_PATH, EXTENDED_NOISELESS_FIGURE_PATH
    return RESULT_PATH, FIGURE_PATH


def _load_results(scenario: str, bpm_noise_um: float):
    result_path, _ = _paths(scenario)
    if result_path.exists():
        return json.loads(result_path.read_text(encoding="utf-8"))
    return {
        "scenario": scenario,
        "scope": (
            "Extended scan: five fitted parameters plus independent quadrupole, BPM "
            "and corrector errors and additive BPM readout noise; no PV access."
            if scenario == "extended" else
            "Extended noiseless scan: five fitted parameters plus independent quadrupole, "
            "BPM and corrector errors; no PV access."
            if scenario == "extended-noiseless" else
            "Matched scan of the five-parameter RFTrack virtual-machine ORM update; "
            "no BPM offsets, alignment errors, additive BPM noise, or PV access."
        ),
        "parameter_sigma": [4.0e-5, 3.0e-5, 7.0e-4, 2.5e-2, 1.8e-2],
        "unmodelled_error_sigma": {
            "individual_quadrupole_relative_k1l": INDIVIDUAL_QUAD_SIGMA,
            "individual_bpm_roll_rad": INDIVIDUAL_BPM_ROLL_SIGMA_RAD,
            "individual_corrector_gain": INDIVIDUAL_CORRECTOR_GAIN_SIGMA,
            "bpm_offset_mm": BPM_OFFSET_SIGMA_MM,
            "bpm_readout_noise_um": bpm_noise_um if scenario == "extended" else 0.0,
        },
        "cases": [],
    }


def _write_figure(results: dict, figure_path: Path):
    cases = results["cases"]
    if not cases:
        raise ValueError("No scan cases are available")
    scales = sorted({float(case["error_scale"]) for case in cases})
    figure, axes = plt.subplots(1, 2, figsize=(10.2, 3.5), constrained_layout=True)
    positions = np.arange(len(scales))
    nominal, updated, improvement = [], [], []
    for scale in scales:
        selected = [case for case in cases if float(case["error_scale"]) == scale]
        nominal.append(np.median([case["heldout_orm_relative_error_nominal"] for case in selected]))
        updated.append(np.median([case["heldout_orm_relative_error_updated"] for case in selected]))
        improvement.append(np.median([case["heldout_orm_improvement"] for case in selected]))
    axes[0].plot(positions, nominal, "o-", label="nominal model", color="tab:orange")
    axes[0].plot(positions, updated, "o-", label="updated model", color="tab:blue")
    axes[0].set(xticks=positions, xticklabels=[f"{scale:g}" for scale in scales],
                xlabel="error scale [sigma]", ylabel="median held-out ORM error",
                title="prediction of unused correctors")
    axes[0].set_yscale("log")
    axes[0].grid(alpha=0.3); axes[0].legend(fontsize=8)
    axes[1].plot(positions, np.asarray(improvement) * 100.0, "o-", color="tab:green")
    axes[1].set(xticks=positions, xticklabels=[f"{scale:g}" for scale in scales],
                xlabel="error scale [sigma]", ylabel="median improvement [%]",
                title="held-out ORM improvement")
    axes[1].grid(alpha=0.3)
    figure.savefig(figure_path, dpi=180)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seeds", default="2003,2004,2005")
    parser.add_argument("--scales", default="0.5,1.0,2.0")
    parser.add_argument(
        "--scenario", choices=("matched", "extended-noiseless", "extended"), default="extended"
    )
    parser.add_argument("--bpm-noise-um", type=float, default=1.0)
    parser.add_argument("--refresh", action="store_true", help="replace existing cases for this scenario")
    parser.add_argument("--summary-only", action="store_true")
    args = parser.parse_args()
    results = _load_results(args.scenario, args.bpm_noise_um)
    if args.refresh:
        results["cases"] = []
    result_path, figure_path = _paths(args.scenario)
    if not args.summary_only:
        module = _load_case_module()
        prepared = _prepare(module)
        for scale in _parse_csv(args.scales, float):
            if scale <= 0.0:
                raise ValueError("scales must be positive")
            for seed in _parse_csv(args.seeds, int):
                key = (seed, scale)
                if any((case["seed"], case["error_scale"]) == key for case in results["cases"]):
                    continue
                row = _case(module, prepared, seed, scale, args.scenario, args.bpm_noise_um)
                results["cases"].append(row)
                print(
                    f"seed={seed}, scale={scale:g}: "
                    f"hold-out {row['heldout_orm_relative_error_nominal']:.3%} -> "
                    f"{row['heldout_orm_relative_error_updated']:.3%}",
                    flush=True,
                )
        results["cases"].sort(key=lambda case: (case["error_scale"], case["seed"]))
        result_path.write_text(json.dumps(results, indent=2), encoding="utf-8")
    _write_figure(results, figure_path)
    print("saved:", result_path)
    print("saved:", figure_path)


if __name__ == "__main__":
    main()
