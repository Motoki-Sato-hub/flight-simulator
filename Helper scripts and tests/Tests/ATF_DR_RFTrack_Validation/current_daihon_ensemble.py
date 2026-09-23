"""Run a small, reproducible current-daihon correction ensemble.

This is deliberately an observable-level ensemble.  The full 6D radiation
envelope is evaluated separately for representative successful seeds because
it is substantially more expensive.  Each subprocess writes its own result
artifact, keyed by random seed.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
ANALYSIS = ROOT.parent / "analysis" / "DR-RFTrack"
PROCEDURE = Path(__file__).with_name("kubo_style_procedure.py")
SEEDS = tuple(int(value) for value in os.environ.get(
    "ATF_DR_CURRENT_ENSEMBLE_SEEDS", "2003,2004,2005,2006,2007"
).split(","))
SKIP_RUN = os.environ.get("ATF_DR_CURRENT_ENSEMBLE_SKIP_RUN", "0") == "1"


def _result_path(seed: int) -> Path:
    return ANALYSIS / (
        "kubo2003_rftrack_procedure_result_20111111b_SF1R_"
        f"probe_0p01mrad_solver_svd_rcond_1em02_seed_{seed}_scale_0p1.json"
    )


def main():
    environment = os.environ.copy()
    environment.update({
        "PYTHONPATH": str(ROOT),
        "MPLBACKEND": "Agg",
        "ATF_DR_KUBO_MAGNET_ERROR_SCALE": "0.1",
        "ATF_DR_KUBO_COD_RCOND": "1e-2",
        "ATF_DR_KUBO_COUPLING_RCOND": "0.4",
        "ATF_DR_KUBO_RESPONSE_KICK_RAD": "1e-5",
        "ATF_DR_KUBO_SKEW_FAMILY": "SF1R",
        "ATF_DR_KUBO_COD_DISPERSION_SOLVER": "svd",
    })
    rows = []
    failures = []
    for seed in SEEDS:
        environment["ATF_DR_KUBO_SEED"] = str(seed)
        log = ANALYSIS / f"current_daihon_ensemble_seed_{seed}.log"
        if SKIP_RUN:
            returncode = 0
        else:
            completed = subprocess.run(
                (sys.executable, str(PROCEDURE)), cwd=ROOT, env=environment,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
            )
            returncode = completed.returncode
            log.write_text(completed.stdout, encoding="utf-8")
        result_path = _result_path(seed)
        if returncode or not result_path.exists():
            failures.append({
                "seed": seed, "returncode": returncode,
                "log": str(log),
            })
            continue
        result = json.loads(result_path.read_text(encoding="utf-8"))
        values = {"seed": seed}
        for stage, metrics in result["stages"].items():
            values[f"{stage}_dy_mm_per_delta"] = metrics["true_dy_all_bpm_mm"]
            values[f"{stage}_cxy"] = metrics["measured_cxy"]
        rows.append(values)
        print(f"seed {seed}: complete", flush=True)

    payload = {
        "scope": (
            "Current 20111111b daihon, 10% Table-I random magnet errors, "
            "one nominal 0.01-mrad response model, observable-level ensemble."
        ),
        "settings": {
            "seeds_requested": list(SEEDS), "magnet_error_scale": 0.1,
            "response_probe_mrad": 0.01, "cod_dispersion_solver": "svd",
        },
        "successful_cases": rows,
        "failures": failures,
    }
    result_path = ANALYSIS / "current_daihon_correction_ensemble.json"
    result_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    stages = ("after_cod", "after_cod_dispersion", "after_coupling")
    labels = ("COD", "COD + Dy", "coupling")
    figure, axes = plt.subplots(1, 2, figsize=(9.8, 3.6), constrained_layout=True)
    for axis, metric, ylabel in (
        (axes[0], "dy_mm_per_delta", "true RMS Dy [mm / delta]"),
        (axes[1], "cxy", "measured coupling observable Cxy"),
    ):
        for row in rows:
            values = [row[f"{stage}_{metric}"] for stage in stages]
            axis.plot(labels, values, "o-", alpha=0.65)
        if rows:
            median = np.median(
                [[row[f"{stage}_{metric}"] for stage in stages] for row in rows], axis=0
            )
            axis.plot(labels, median, "ko--", lw=2.2, label="median")
            axis.legend(fontsize=8)
        axis.set(ylabel=ylabel)
        axis.grid(alpha=0.3)
    figure.suptitle(
        f"Current daihon correction ensemble: {len(rows)}/{len(SEEDS)} successful seeds"
    )
    figure.savefig(ANALYSIS / "current_daihon_correction_ensemble.png", dpi=180)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
