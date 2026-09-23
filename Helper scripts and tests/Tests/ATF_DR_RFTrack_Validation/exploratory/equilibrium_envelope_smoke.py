"""Evaluate ATF DR radiation equilibrium with the one-turn envelope method."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt

from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_emittance import equilibrium_emittance


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--quantum-particles", type=int, default=32)
    arguments = parser.parse_args()
    result = equilibrium_emittance(quantum_particles=arguments.quantum_particles)

    analysis = Path(__file__).resolve().parents[5] / "analysis" / "DR-RFTrack"
    payload = {
        "scope": (
            "Nominal, error-free ATF DR radiation equilibrium from a numerical "
            "one-turn damping map and quantum-diffusion covariance. Vertical "
            "emittance is expected to vanish without coupling/error sources."
        ),
        "quantum_particles": result.quantum_particles,
        "spectral_radius": result.spectral_radius,
        "damping_turns_one_over_e": result.damping_turns_one_over_e,
        "synchronous_orbit": {
            "coordinates": result.synchronous_orbit.coordinates.tolist(),
            "residual": result.synchronous_orbit.residual.tolist(),
            "iterations": result.synchronous_orbit.iterations,
        },
        "projected_x_m_rad": result.projected_x_m_rad,
        "projected_y_m_rad": result.projected_y_m_rad,
        "eigen_emittances_m_rad": list(result.eigen_emittances_m_rad),
        "one_turn_map": result.one_turn_map.tolist(),
        "diffusion_native": result.diffusion_native.tolist(),
        "covariance_native": result.covariance_native.tolist(),
    }
    (analysis / "atf_dr_equilibrium_envelope_result.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )

    labels = ("projected x", "upper eigen", "projected y", "lower eigen")
    values_pm = (
        result.projected_x_m_rad * 1e12,
        result.eigen_emittances_m_rad[1] * 1e12,
        result.projected_y_m_rad * 1e12,
        result.eigen_emittances_m_rad[0] * 1e12,
    )
    figure, axis = plt.subplots(figsize=(7.6, 3.7), constrained_layout=True)
    bars = axis.bar(labels, values_pm, color=("tab:blue", "tab:blue", "tab:orange", "tab:orange"))
    axis.bar_label(bars, labels=[f"{value:.3g}" for value in values_pm], padding=3)
    axis.set(
        title="ATF DR RF-Track nominal radiation equilibrium",
        ylabel="geometric emittance [pm rad]",
    )
    axis.grid(axis="y", alpha=0.3)
    axis.text(
        0.98,
        0.96,
        f"spectral radius = {result.spectral_radius:.9f}\n"
        f"1/e damping = {result.damping_turns_one_over_e:.0f} turns\n"
        f"quantum sample = {result.quantum_particles}",
        transform=axis.transAxes,
        ha="right",
        va="top",
        bbox={"facecolor": "white", "edgecolor": "0.7"},
    )
    figure.savefig(analysis / "atf_dr_equilibrium_envelope_result.png", dpi=180)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
