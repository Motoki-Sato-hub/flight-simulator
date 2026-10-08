"""End-to-end skew correction check using radiation-equilibrium emittance.

One SF1R skew fault is deliberately outside the SD1R correction family.  A
nominal RF-Track coupling response drives SD1R correctors, then the same fault
and correction are evaluated with radiation/RF equilibrium emittance.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_correction import ATFDRRingCorrection
from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_emittance import (
    build_equilibrium_lattice,
    equilibrium_emittance_for_lattices,
    find_synchronous_orbit,
)
from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_lattice import build_atf_dr_lattice


FAULT_NAME = "SF1R.1$SKEW"
FAULT_K1L = 0.05
PROBES = ("ZH1R", "ZH2R")
# The 2011 reference lattice remains periodic for this 0.01-mrad linear
# probe; 0.1 mrad is the Kubo-paper value but exceeds its horizontal branch.
PROBE_KICK_RAD = 1e-5
# Retain the response modes required by this controlled SF1R-to-SD1R test.
# The Kubo full-ensemble cutoff will be calibrated separately.
RCOND = 0.01
QUANTUM_PARTICLES = 32


def _coupling_signal(correction):
    _, signal = correction._measure_vertical_response_to_horizontal_probes(
        PROBES,
        bpm_names=correction.bpm_names,
        probe_perturbation=PROBE_KICK_RAD,
    )
    return signal.reshape(-1)


def _cxy(correction, signal):
    horizontal = []
    for name in PROBES:
        element = correction._single_element(name)
        original = correction._get_corrector_kick(element)
        plus, minus = original.copy(), original.copy()
        plus[0] += PROBE_KICK_RAD
        minus[0] -= PROBE_KICK_RAD
        correction._set_corrector_kick(element, plus)
        positive = correction.find_closed_orbit().bpm_positions[:, 0]
        correction._set_corrector_kick(element, minus)
        negative = correction.find_closed_orbit().bpm_positions[:, 0]
        correction._set_corrector_kick(element, original)
        horizontal.append((positive - negative) / (2.0 * PROBE_KICK_RAD))
    horizontal = np.asarray(horizontal).reshape(-1)
    return float(np.sqrt(np.mean(signal**2) / np.mean(horizontal**2)))


def _set_skews(correction, names, values):
    for name, value in zip(names, values):
        correction.set_skew_strength(name, correction.get_skew_strength(name) + value)


def _equilibrium(fault, names=(), values=(), initial_coordinates=None):
    deterministic = build_equilibrium_lattice(quantum=False)
    quantum = build_equilibrium_lattice(quantum=True)
    for lattice in (deterministic, quantum):
        correction = ATFDRRingCorrection(lattice)
        correction.set_skew_strength(FAULT_NAME, fault)
        _set_skews(correction, names, values)
    return equilibrium_emittance_for_lattices(
        deterministic,
        quantum,
        quantum_particles=QUANTUM_PARTICLES,
        initial_coordinates=initial_coordinates,
    )


def main():
    analysis = Path(__file__).resolve().parents[5] / "analysis" / "DR-RFTrack"
    cache = np.load(analysis / "kubo2003_nominal_kick_response_cache.npz", allow_pickle=True)
    names = tuple(cache["skew_names"].tolist())
    matrix = cache["coupling"]

    nominal_sync = find_synchronous_orbit(
        build_equilibrium_lattice(quantum=False), max_iterations=50
    ).coordinates
    machine = ATFDRRingCorrection(build_atf_dr_lattice())
    machine.set_skew_strength(FAULT_NAME, FAULT_K1L)
    before_signal = _coupling_signal(machine)
    before_cxy = _cxy(machine, before_signal)
    correction, _, _ = ATFDRRingCorrection._svd_solve(matrix, -before_signal, RCOND)
    _set_skews(machine, names, correction)
    after_signal = _coupling_signal(machine)
    after_cxy = _cxy(machine, after_signal)

    before = _equilibrium(FAULT_K1L, initial_coordinates=nominal_sync)
    after = _equilibrium(
        FAULT_K1L, names, correction, initial_coordinates=nominal_sync
    )
    payload = {
        "scope": (
            "Controlled end-to-end skew test: one SF1R fault outside the SD1R "
            "correction family; nominal coupling response; radiation/RF envelope "
            "emittance. This is not yet the full Kubo Table-I ensemble."
        ),
        "fault": {"name": FAULT_NAME, "skew_k1l": FAULT_K1L},
        "correction": {
            "family": "SD1R",
            "actuators": len(names),
            "rcond": RCOND,
            "rms_skew_k1l": float(np.sqrt(np.mean(correction**2))),
            "max_skew_k1l": float(np.max(np.abs(correction))),
        },
        "before": {
            "cxy": before_cxy,
            "projected_y_pm_rad": before.projected_y_m_rad * 1e12,
            "lower_eigen_pm_rad": before.eigen_emittances_m_rad[0] * 1e12,
        },
        "after": {
            "cxy": after_cxy,
            "projected_y_pm_rad": after.projected_y_m_rad * 1e12,
            "lower_eigen_pm_rad": after.eigen_emittances_m_rad[0] * 1e12,
        },
    }
    (analysis / "atf_dr_skew_emittance_correction_result.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    figure, axes = plt.subplots(1, 2, figsize=(8.6, 3.5), constrained_layout=True)
    labels = ("before", "after")
    bars = axes[0].bar(
        labels,
        (payload["before"]["projected_y_pm_rad"], payload["after"]["projected_y_pm_rad"]),
        color=("tab:red", "tab:green"),
    )
    axes[0].bar_label(bars, fmt="%.3g", padding=3)
    axes[0].set(title="radiation-equilibrium projected y", ylabel="emittance [pm rad]")
    axes[0].grid(axis="y", alpha=0.3)
    bars = axes[1].bar(
        labels, (before_cxy, after_cxy), color=("tab:red", "tab:green")
    )
    axes[1].bar_label(bars, fmt="%.3g", padding=3)
    axes[1].set(title="transverse coupling observable", ylabel="Cxy")
    axes[1].grid(axis="y", alpha=0.3)
    figure.suptitle("ATF DR RF-Track: skew correction to equilibrium emittance")
    figure.savefig(analysis / "atf_dr_skew_emittance_correction_result.png", dpi=180)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
