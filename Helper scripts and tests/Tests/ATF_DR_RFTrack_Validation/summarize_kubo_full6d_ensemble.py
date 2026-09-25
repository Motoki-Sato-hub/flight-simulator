"""Summarise compatible RF-Track Kubo full-6D seed results.

This is a small-seed procedural check, not a replacement for Kubo's published
500-seed Table-II statistic.  It only groups artifacts with the explicitly
selected current-daihon settings encoded in ``PATTERN``.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ANALYSIS_DIR = Path(__file__).resolve().parents[4] / "analysis" / "DR-RFTrack"
ERROR_SCALE = float(os.environ.get("ATF_DR_KUBO_SUMMARY_SCALE", "0.1"))
if not 0.0 < ERROR_SCALE <= 1.0:
    raise ValueError("ATF_DR_KUBO_SUMMARY_SCALE must be in (0, 1]")
SCALE_LABEL = f"{ERROR_SCALE:g}".replace(".", "p")
PATTERN = (
    "kubo2003_rftrack_emittance_envelope_v2_20111111b_SF1R_"
    f"full_6d_seed_*_scale_{SCALE_LABEL}_solver_svd_codrcond_0p01_"
    "couplingrcond_0p4_qparticles_128.json"
)
STAGES = ("before_cod", "after_cod", "after_cod_dispersion", "after_coupling")
LABELS = ("before COD", "COD", "COD + Dy", "skew")
PAPER_TABLE_II = np.array((np.nan, 22.8, 16.7, 5.8))


def _stage_matrix(payloads: list[dict], section: str, field: str) -> np.ndarray:
    """Return stage data, retaining failed or absent stages as NaN."""
    return np.asarray([
        [item.get(section, {}).get(stage, {}).get(field, np.nan) for stage in STAGES]
        for item in payloads
    ], dtype=float)


def main() -> None:
    paths = sorted(ANALYSIS_DIR.glob(PATTERN))
    if not paths:
        raise FileNotFoundError(f"No compatible artifacts: {ANALYSIS_DIR / PATTERN}")
    payloads = [json.loads(path.read_text(encoding="utf-8")) for path in paths]
    seeds = [item["model_settings"]["seed"] for item in payloads]
    vertical = _stage_matrix(payloads, "equilibrium_emittance", "vertical_like_eigen_pm_rad")
    dy = _stage_matrix(payloads, "correction_observables", "true_dy_all_bpm_mm")
    cxy = _stage_matrix(payloads, "correction_observables", "measured_cxy")
    expected_seeds = {
        int(seed) for seed in os.environ.get("ATF_DR_KUBO_EXPECTED_SEEDS", "").split(",")
        if seed.strip()
    }

    result = {
        "scope": (
            f"Pilot ensemble on the 2011 daihon: {100 * ERROR_SCALE:g}% Kubo Table-I magnet "
            "errors, fixed BPM errors, full-6D one-turn radiation envelope. "
            "This is not Kubo Table-II's 500-seed historical-lattice statistic."
        ),
        "artifact_pattern": PATTERN,
        "seeds": seeds,
        "sample_count": len(seeds),
        "missing_expected_artifacts": sorted(expected_seeds - set(seeds)),
        "stages": list(STAGES),
        "successful_emittance_samples_per_stage": np.isfinite(vertical).sum(axis=0).tolist(),
        "vertical_like_eigen_pm_rad": {
            "median": np.nanmedian(vertical, axis=0).tolist(),
            "minimum": np.nanmin(vertical, axis=0).tolist(),
            "maximum": np.nanmax(vertical, axis=0).tolist(),
        },
        "true_dy_rms_mm_per_delta": {
            "median": np.nanmedian(dy, axis=0).tolist(),
        },
        "measured_cxy": {"median": np.nanmedian(cxy, axis=0).tolist()},
        "kubo_table_ii_vertical_emittance_pm_rad": [None, 22.8, 16.7, 5.8],
    }
    output = ANALYSIS_DIR / (
        f"kubo2003_rftrack_full6d_v2_scale_{SCALE_LABEL}_ensemble.json"
    )
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")

    figure, axes = plt.subplots(1, 3, figsize=(11.5, 3.3), constrained_layout=True)
    for row, seed in zip(vertical, seeds):
        axes[0].plot(LABELS, row, "o-", alpha=0.5, label=f"seed {seed}")
    axes[0].plot(LABELS, np.nanmedian(vertical, axis=0), "ko-", lw=2, label="median")
    axes[0].plot(LABELS, PAPER_TABLE_II, "k^:", label="Kubo Table II")
    axes[0].set(ylabel="vertical-like emittance [pm rad]", title="equilibrium emittance")
    axes[0].legend(fontsize=7, ncol=2)
    axes[1].plot(LABELS, np.nanmedian(dy, axis=0), "o-", color="tab:blue")
    axes[1].set(ylabel="median RMS Dy [mm/delta]", title="vertical dispersion")
    axes[2].plot(LABELS, np.nanmedian(cxy, axis=0), "o-", color="tab:orange")
    axes[2].set(ylabel="median Cxy", title="coupling")
    for axis in axes:
        axis.grid(alpha=0.3)
        axis.tick_params(axis="x", rotation=15, labelsize=8)
    figure.suptitle(f"RF-Track Kubo full-6D pilot ensemble ({len(seeds)} seeds)")
    figure.savefig(
        ANALYSIS_DIR / (
            f"kubo2003_rftrack_full6d_v2_scale_{SCALE_LABEL}_ensemble.png"
        ),
        dpi=180,
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
