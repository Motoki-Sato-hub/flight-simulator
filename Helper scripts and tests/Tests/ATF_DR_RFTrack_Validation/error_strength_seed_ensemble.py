"""Test error-strength correction completion across independent random seeds.

This deliberately uses the same nominal response matrix and correction knobs
for every case.  It reports a solver failure separately from an unsuccessful
radiation-envelope evaluation.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


FLIGHT_SIMULATOR = Path(__file__).resolve().parents[3]
ANALYSIS = FLIGHT_SIMULATOR.parent / "analysis" / "DR-RFTrack"
ENVELOPE = Path(__file__).with_name("six_d_equilibrium_envelope.py")


def _values(text: str, kind):
    return tuple(kind(value.strip()) for value in text.split(",") if value.strip())


def _label(value: float) -> str:
    return f"{value:g}".replace(".", "p").replace("-", "m")


def _result_path(scale, seed, particles, first_gain, weight, ramp_initial, ramp_max, ramp_steps):
    stem = (
        "kubo2003_rftrack_emittance_envelope_20111111b_SF1R_full_6d_"
        f"seed_{seed}_scale_{_label(scale)}_solver_svd_codrcond_0p01_"
        f"couplingrcond_0p4_qparticles_{particles}"
    )
    if not np.isclose(first_gain, 0.7):
        stem += f"_firstgain_{_label(first_gain)}"
    if not np.isclose(weight, 0.05):
        stem += f"_dyweight_{_label(weight)}"
    stem += f"_ramp_{_label(ramp_initial)}_{_label(ramp_max)}_steps_{ramp_steps}"
    return ANALYSIS / f"{stem}.json"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seeds", default="2003,2004,2005,2006")
    parser.add_argument("--scales", default="0.8,0.9")
    parser.add_argument("--quantum-particles", type=int, default=10)
    parser.add_argument("--first-gain", type=float, default=1.0)
    parser.add_argument("--dispersion-weight", type=float, default=0.1)
    parser.add_argument("--ramp-initial", type=float, default=0.25)
    parser.add_argument("--ramp-max", type=float, default=0.25)
    parser.add_argument("--ramp-feedback-steps", type=int, default=20)
    args = parser.parse_args()
    seeds = _values(args.seeds, int)
    scales = _values(args.scales, float)
    if not seeds or any(not 0.0 < scale <= 1.0 for scale in scales):
        raise ValueError("provide seeds and error scales in (0, 1]")

    environment = os.environ.copy()
    environment.update({
        "PYTHONPATH": str(FLIGHT_SIMULATOR), "MPLBACKEND": "Agg",
        "ATF_DR_KUBO_COD_DISPERSION_SOLVER": "svd",
        "ATF_DR_KUBO_EMITTANCE_MODEL": "full_6d",
        "ATF_DR_KUBO_QUANTUM_PARTICLES": str(args.quantum_particles),
        "ATF_DR_KUBO_FIRST_GAIN": f"{args.first_gain:g}",
        "ATF_DR_KUBO_DISPERSION_WEIGHT": f"{args.dispersion_weight:g}",
        "ATF_DR_KUBO_RAMP_INITIAL_INCREMENT": f"{args.ramp_initial:g}",
        "ATF_DR_KUBO_RAMP_MAX_INCREMENT": f"{args.ramp_max:g}",
        "ATF_DR_KUBO_RAMP_MAX_FEEDBACK_STEPS": str(args.ramp_feedback_steps),
    })
    cases = []
    for scale in scales:
        for seed in seeds:
            environment["ATF_DR_KUBO_SEED"] = str(seed)
            environment["ATF_DR_KUBO_MAGNET_ERROR_SCALE"] = f"{scale:g}"
            result_path = _result_path(
                scale, seed, args.quantum_particles, args.first_gain,
                args.dispersion_weight, args.ramp_initial, args.ramp_max,
                args.ramp_feedback_steps,
            )
            log_path = ANALYSIS / (
                f"error_strength_seed_ensemble_scale_{_label(scale)}_seed_{seed}.log"
            )
            completed = subprocess.run(
                (sys.executable, str(ENVELOPE)), cwd=FLIGHT_SIMULATOR,
                env=environment, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, text=True,
            )
            log_path.write_text(completed.stdout, encoding="utf-8")
            row = {"seed": seed, "error_scale_of_kubo_table_i": scale,
                   "returncode": completed.returncode, "log_path": str(log_path),
                   "result_path": str(result_path)}
            if completed.returncode == 0 and result_path.exists():
                payload = json.loads(result_path.read_text(encoding="utf-8"))
                row["correction_sequence_completed"] = all(
                    stage in payload["correction_observables"]
                    for stage in ("after_cod", "after_cod_dispersion", "after_coupling")
                )
                final = payload["equilibrium_emittance"]["after_coupling"]
                row["final_envelope_status"] = final["status"]
                if final["status"] == "ok":
                    row["final_vertical_like_pm_rad"] = final["vertical_like_eigen_pm_rad"]
                row["status"] = "completed"
            else:
                row["status"] = "process_failed"
            cases.append(row)
            print(row, flush=True)

    payload = {"scope": "Independent-seed scan of the current RF-Track Kubo workflow.",
               "settings": vars(args), "cases": cases}
    output = ANALYSIS / "atf_dr_rftrack_error_strength_seed_ensemble.json"
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    figure, axis = plt.subplots(figsize=(6.5, 3.8), constrained_layout=True)
    for scale in scales:
        rows = [row for row in cases if row["error_scale_of_kubo_table_i"] == scale]
        x = [row["seed"] for row in rows if "final_vertical_like_pm_rad" in row]
        y = [row["final_vertical_like_pm_rad"] for row in rows if "final_vertical_like_pm_rad" in row]
        if x:
            axis.plot(x, y, "o-", label=f"{scale:g} scale")
    axis.set(xlabel="random-error seed", ylabel="final vertical-like emittance [pm rad]")
    axis.grid(alpha=0.3)
    axis.legend(fontsize=8)
    figure.savefig(ANALYSIS / "atf_dr_rftrack_error_strength_seed_ensemble.png", dpi=180)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
