"""Short multi-turn validation of the RF-Track ATF DR radiation/RF model.

This is a smoke test of the equilibrium-emittance infrastructure, not a claim
that three turns have reached damping equilibrium.  Increase ``--turns`` and
``--particles`` together only after a radiation-step convergence study.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

import matplotlib.pyplot as plt

from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_emittance import (
    build_equilibrium_lattice,
    find_synchronous_orbit,
    gaussian_bunch,
    track_emittance,
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--turns", type=int, default=3)
    parser.add_argument("--particles", type=int, default=32)
    parser.add_argument("--radiation-steps", type=int, default=10)
    parser.add_argument("--sample-every", type=int, default=1)
    arguments = parser.parse_args()

    deterministic = build_equilibrium_lattice(
        quantum=False, radiation_steps=arguments.radiation_steps
    )
    synchronous = find_synchronous_orbit(deterministic, max_iterations=50)
    lattice = build_equilibrium_lattice(
        quantum=True, radiation_steps=arguments.radiation_steps
    )
    bunch = gaussian_bunch(
        synchronous,
        particles=arguments.particles,
        emit_x_m_rad=3e-9,
        emit_y_m_rad=1e-11,
        beta_x_m=5.0,
        beta_y_m=2.0,
        sigma_delta=1e-4,
        sigma_ct_mm=1.0,
        seed=3,
    )
    _, samples = track_emittance(
        lattice,
        bunch,
        turns=arguments.turns,
        sample_every=arguments.sample_every,
    )

    analysis = Path(__file__).resolve().parents[5] / "analysis" / "DR-RFTrack"
    result = {
        "scope": (
            "Radiation/RF multi-turn emittance infrastructure smoke test. "
            "The selected turn count is not assumed to be damping equilibrium."
        ),
        "configuration": vars(arguments),
        "synchronous_orbit": {
            "coordinates": synchronous.coordinates.tolist(),
            "residual": synchronous.residual.tolist(),
            "iterations": synchronous.iterations,
        },
        "samples": [
            {
                **asdict(sample),
                "covariance_m_rad": sample.covariance_m_rad.tolist(),
            }
            for sample in samples
        ],
    }
    result_path = analysis / "atf_dr_equilibrium_emittance_smoke_result.json"
    result_path.write_text(json.dumps(result, indent=2), encoding="utf-8")

    turns = [sample.turns for sample in samples]
    projected_x = [sample.projected_x_m_rad * 1e9 for sample in samples]
    projected_y = [sample.projected_y_m_rad * 1e12 for sample in samples]
    eigen_low = [sample.eigen_emittances_m_rad[0] * 1e12 for sample in samples]
    figure, axes = plt.subplots(1, 2, figsize=(9.2, 3.5), constrained_layout=True)
    axes[0].plot(turns, projected_x, "o-", label="projected x")
    axes[0].set(xlabel="turn", ylabel="emittance [nm rad]", title="Horizontal")
    axes[0].grid(alpha=0.3)
    axes[1].plot(turns, projected_y, "o-", label="projected y")
    axes[1].plot(turns, eigen_low, "s--", label="lower eigen-emittance")
    axes[1].set(xlabel="turn", ylabel="emittance [pm rad]", title="Vertical")
    axes[1].grid(alpha=0.3)
    axes[1].legend()
    figure.suptitle("ATF DR RF-Track radiation/RF emittance smoke test")
    figure.savefig(
        analysis / "atf_dr_equilibrium_emittance_smoke_result.png", dpi=180
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
