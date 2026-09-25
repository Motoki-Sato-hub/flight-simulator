"""Classify the usable Table-I error range of the RF-Track Kubo workflow.

Each point is classified separately for periodic-orbit/correction completion
and for the three 6D equilibrium-emittance evaluations.  Hence a CPU timeout,
a lost periodic orbit, and an unstable radiation envelope are never presented
as the same kind of failure.
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


def _values(text):
    return tuple(float(value.strip()) for value in text.split(",") if value.strip())


def _label(value):
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
    parser.add_argument("--scales", default="0.1,0.2,0.3,0.5,0.7,1.0")
    parser.add_argument("--seed", type=int, default=2003)
    parser.add_argument("--quantum-particles", type=int, default=32)
    parser.add_argument("--first-gain", type=float, default=1.0)
    parser.add_argument("--dispersion-weight", type=float, default=0.1)
    parser.add_argument("--ramp-initial", type=float, default=0.25)
    parser.add_argument("--ramp-max", type=float, default=0.25)
    parser.add_argument("--ramp-feedback-steps", type=int, default=20)
    parser.add_argument("--skip-run", action="store_true")
    args = parser.parse_args()
    scales = _values(args.scales)
    if any(not 0.0 < value <= 1.0 for value in scales):
        raise ValueError("all error scales must be in (0, 1]")
    if args.quantum_particles < 10:
        raise ValueError("quantum-particles must be at least 10")

    environment = os.environ.copy()
    environment.update({
        "PYTHONPATH": str(FLIGHT_SIMULATOR), "MPLBACKEND": "Agg",
        "ATF_DR_KUBO_SEED": str(args.seed),
        "ATF_DR_KUBO_COD_DISPERSION_SOLVER": "svd",
        "ATF_DR_KUBO_EMITTANCE_MODEL": "full_6d",
        "ATF_DR_KUBO_QUANTUM_PARTICLES": str(args.quantum_particles),
        "ATF_DR_KUBO_FIRST_GAIN": f"{args.first_gain:g}",
        "ATF_DR_KUBO_DISPERSION_WEIGHT": f"{args.dispersion_weight:g}",
        "ATF_DR_KUBO_RAMP_INITIAL_INCREMENT": f"{args.ramp_initial:g}",
        "ATF_DR_KUBO_RAMP_MAX_INCREMENT": f"{args.ramp_max:g}",
        "ATF_DR_KUBO_RAMP_MAX_FEEDBACK_STEPS": str(args.ramp_feedback_steps),
    })
    rows = []
    for scale in scales:
        environment["ATF_DR_KUBO_MAGNET_ERROR_SCALE"] = f"{scale:g}"
        result_path = _result_path(
            scale, args.seed, args.quantum_particles, args.first_gain,
            args.dispersion_weight, args.ramp_initial, args.ramp_max,
            args.ramp_feedback_steps,
        )
        log_path = ANALYSIS / f"error_strength_scan_scale_{_label(scale)}.log"
        if not args.skip_run:
            completed = subprocess.run(
                (sys.executable, str(ENVELOPE)), cwd=FLIGHT_SIMULATOR,
                env=environment, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, text=True,
            )
            log_path.write_text(completed.stdout, encoding="utf-8")
            returncode = completed.returncode
        else:
            returncode = 0
        row = {
            "error_scale_of_kubo_table_i": scale,
            "result_path": str(result_path), "log_path": str(log_path),
            "returncode": returncode,
        }
        if result_path.exists():
            payload = json.loads(result_path.read_text(encoding="utf-8"))
            row["correction_sequence_completed"] = all(
                stage in payload["correction_observables"]
                for stage in ("after_cod", "after_cod_dispersion", "after_coupling")
            )
            for stage, item in payload["equilibrium_emittance"].items():
                row[f"{stage}_envelope_status"] = item["status"]
                if item["status"] == "ok":
                    row[f"{stage}_vertical_like_pm_rad"] = item["vertical_like_eigen_pm_rad"]
            row["status"] = "completed"
        elif args.skip_run:
            row["status"] = "not_available"
        else:
            row["status"] = "process_failed"
        rows.append(row)
        print(row, flush=True)

    payload = {
        "scope": (
            "Current 20111111b daihon error-strength scan. A completed result "
            "does not by itself establish multi-seed robustness."
        ),
        "settings": vars(args), "cases": rows,
    }
    output = ANALYSIS / "atf_dr_rftrack_error_strength_scan.json"
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    figure, axis = plt.subplots(figsize=(7.2, 4.2), constrained_layout=True)
    completed = [row for row in rows if row["status"] == "completed"]
    for stage, marker in zip(("after_cod", "after_cod_dispersion", "after_coupling"), ("o", "s", "^")):
        x = [row["error_scale_of_kubo_table_i"] for row in completed if f"{stage}_vertical_like_pm_rad" in row]
        y = [row[f"{stage}_vertical_like_pm_rad"] for row in completed if f"{stage}_vertical_like_pm_rad" in row]
        if x:
            axis.plot(x, y, marker + "-", label=stage.replace("after_", ""))
    axis.set(xlabel="Table-I magnet-error scale", ylabel="vertical-like emittance [pm rad]")
    axis.grid(alpha=0.3)
    axis.legend(fontsize=8)
    axis.set_title("RF-Track correction validity versus error strength")
    figure.savefig(ANALYSIS / "atf_dr_rftrack_error_strength_scan.png", dpi=180)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
