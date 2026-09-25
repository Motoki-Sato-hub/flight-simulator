"""Grid-scan Kubo correction parameters against RF-Track 6D emittance.

Kubo's operational method deliberately uses BPM observables rather than a
beam-size monitor.  In a simulation, however, RF-Track's radiation-equilibrium
calculation supplies the unobservable vertical-like eigen-emittance.  This
script uses it *offline* to compare otherwise identical Kubo sequences.

It is a model-study tool, not an online correction algorithm and not a claim
that emittance is directly available from BPM data in the real ATF DR.
"""

from __future__ import annotations

import argparse
import itertools
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


def _numbers(text):
    return tuple(float(value.strip()) for value in text.split(",") if value.strip())


def _seeds(text):
    return tuple(int(value.strip()) for value in text.split(",") if value.strip())


def _label(value):
    return f"{value:g}".replace(".", "p").replace("-", "m")


def _result_path(
    seed, first_gain, weight, quantum_particles, ramp_initial, ramp_max, ramp_steps
):
    stem = (
        "kubo2003_rftrack_emittance_envelope_20111111b_SF1R_full_6d_"
        f"seed_{seed}_scale_0p1_solver_svd_codrcond_0p01_couplingrcond_0p4_"
        f"qparticles_{quantum_particles}"
    )
    if not np.isclose(first_gain, 0.7):
        stem += f"_firstgain_{_label(first_gain)}"
    if not np.isclose(weight, 0.05):
        stem += f"_dyweight_{_label(weight)}"
    stem += f"_ramp_{_label(ramp_initial)}_{_label(ramp_max)}_steps_{ramp_steps}"
    return ANALYSIS / f"{stem}.json"


def _vertical_emittance(payload, stage):
    item = payload["equilibrium_emittance"].get(stage, {})
    return item.get("vertical_like_eigen_pm_rad") if item.get("status") == "ok" else None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seeds", default="2003,2004")
    parser.add_argument("--dispersion-weights", default="0.02,0.05,0.1")
    parser.add_argument("--first-gains", default="0.7")
    parser.add_argument("--quantum-particles", type=int, default=32)
    parser.add_argument("--ramp-initial", type=float, default=0.25)
    parser.add_argument("--ramp-max", type=float, default=0.25)
    parser.add_argument("--ramp-feedback-steps", type=int, default=20)
    parser.add_argument("--skip-run", action="store_true")
    args = parser.parse_args()
    seeds = _seeds(args.seeds)
    weights = _numbers(args.dispersion_weights)
    gains = _numbers(args.first_gains)
    if not seeds or not weights or not gains:
        raise ValueError("at least one seed, weight, and first gain are required")
    if any(value <= 0.0 for value in weights) or any(not 0.0 < value <= 1.0 for value in gains):
        raise ValueError("weights must be positive and first gains must be in (0, 1]")
    if args.quantum_particles < 10:
        raise ValueError("quantum-particles must be at least 10")
    if not 0.0 < args.ramp_initial <= args.ramp_max <= 1.0:
        raise ValueError("ramp values must satisfy 0 < initial <= maximum <= 1")
    if args.ramp_feedback_steps < 0:
        raise ValueError("ramp-feedback-steps must be non-negative")

    environment = os.environ.copy()
    environment.update({
        "PYTHONPATH": str(FLIGHT_SIMULATOR),
        "MPLBACKEND": "Agg",
        "ATF_DR_KUBO_MAGNET_ERROR_SCALE": "0.1",
        "ATF_DR_KUBO_COD_DISPERSION_SOLVER": "svd",
        "ATF_DR_KUBO_COD_RCOND": "1e-2",
        "ATF_DR_KUBO_COUPLING_RCOND": "0.4",
        "ATF_DR_KUBO_RESPONSE_KICK_RAD": "1e-5",
        "ATF_DR_KUBO_SKEW_FAMILY": "SF1R",
        "ATF_DR_KUBO_EMITTANCE_MODEL": "full_6d",
        "ATF_DR_KUBO_QUANTUM_PARTICLES": str(args.quantum_particles),
        # Numerical continuation only: four equal 10%-error substeps.
        "ATF_DR_KUBO_RAMP_INITIAL_INCREMENT": f"{args.ramp_initial:g}",
        "ATF_DR_KUBO_RAMP_MAX_INCREMENT": f"{args.ramp_max:g}",
        "ATF_DR_KUBO_RAMP_MAX_FEEDBACK_STEPS": str(args.ramp_feedback_steps),
    })
    rows = []
    for seed, first_gain, weight in itertools.product(seeds, gains, weights):
        environment.update({
            "ATF_DR_KUBO_SEED": str(seed),
            "ATF_DR_KUBO_FIRST_GAIN": f"{first_gain:g}",
            "ATF_DR_KUBO_DISPERSION_WEIGHT": f"{weight:g}",
        })
        result_path = _result_path(
            seed, first_gain, weight, args.quantum_particles,
            args.ramp_initial, args.ramp_max, args.ramp_feedback_steps,
        )
        log_path = ANALYSIS / (
            f"emittance_scan_seed_{seed}_gain_{_label(first_gain)}_"
            f"weight_{_label(weight)}.log"
        )
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
            "seed": seed,
            "first_gain": first_gain,
            "dispersion_weight_r": weight,
            "result_path": str(result_path),
            "log_path": str(log_path),
            "returncode": returncode,
        }
        if returncode == 0 and result_path.exists():
            payload = json.loads(result_path.read_text(encoding="utf-8"))
            for stage in ("after_cod", "after_cod_dispersion", "after_coupling"):
                row[f"vertical_like_{stage}_pm_rad"] = _vertical_emittance(payload, stage)
                row[f"dy_{stage}_mm_per_delta"] = payload["correction_observables"][stage]["true_dy_all_bpm_mm"]
                row[f"cxy_{stage}"] = payload["correction_observables"][stage]["measured_cxy"]
            row["status"] = "ok"
        elif args.skip_run:
            row["status"] = "not_available"
        else:
            row["status"] = "failed"
        rows.append(row)
        print(row, flush=True)

    payload = {
        "scope": (
            "Offline RF-Track 6D-emittance grid scan on the current 20111111b "
            "daihon with 10%-Table-I random magnet errors. The optimization "
            "target is simulation truth, not a real-time BPM observable."
        ),
        "settings": {
            "seeds": list(seeds), "dispersion_weights_r": list(weights),
            "first_gains": list(gains), "quantum_particles": args.quantum_particles,
            "ramp_initial": args.ramp_initial, "ramp_max": args.ramp_max,
            "ramp_feedback_steps": args.ramp_feedback_steps,
            "magnet_error_scale_of_table_i": 0.1,
        },
        "cases": rows,
    }
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    summary = ANALYSIS / "atf_dr_rftrack_emittance_guided_parameter_scan.json"
    summary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    figure, axis = plt.subplots(figsize=(7.2, 4.2), constrained_layout=True)
    for gain in gains:
        subset = [row for row in rows if row["status"] == "ok" and np.isclose(row["first_gain"], gain)]
        for weight in weights:
            values = [
                row["vertical_like_after_coupling_pm_rad"]
                for row in subset if np.isclose(row["dispersion_weight_r"], weight)
                and row["vertical_like_after_coupling_pm_rad"] is not None
            ]
            if values:
                axis.scatter(
                    [weight], [np.median(values)], s=55,
                    label=f"gain={gain:g}" if weight == weights[0] else None,
                )
    axis.set(xscale="log", xlabel="COD/dispersion weight r", ylabel="median vertical-like emittance [pm rad]")
    axis.grid(alpha=0.3, which="both")
    if len(gains) > 1:
        axis.legend(fontsize=8)
    axis.set_title("Offline emittance-guided scan: current daihon, 10% Table-I errors")
    figure.savefig(ANALYSIS / "atf_dr_rftrack_emittance_guided_parameter_scan.png", dpi=180)
    print(f"Wrote {summary}")


if __name__ == "__main__":
    main()
