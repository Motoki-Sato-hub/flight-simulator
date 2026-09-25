"""Audit RF-Track/SAD unit boundaries used by the ATF DR correction study.

The correction public API uses metres for alignment input and radians for
steerer kicks, while RF-Track reports transverse particle coordinates in mm
and exposes Corrector kicks in mrad.  This script tests those boundaries on
the actual 2011 DR lattice and records an auditable JSON result.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_correction import ATFDRRingCorrection
from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_lattice import build_atf_dr_lattice


FLIGHT_SIMULATOR = Path(__file__).resolve().parents[3]
ANALYSIS = FLIGHT_SIMULATOR.parent / "analysis" / "DR-RFTrack"
ALIGNMENT_TEST_M = 30e-6
# Kubo's nominal-response perturbation: 0.1 mrad = 1e-4 rad.
KICK_TEST_RAD = 1e-4
PROBES_RAD = (1e-6, 3e-6, 1e-5)


def _one(lattice, name):
    element = lattice[name]
    if isinstance(element, list):
        if len(element) != 1:
            raise RuntimeError(f"Expected exactly one {name}")
        return element[0]
    return element


def _paper_x_correctors(machine):
    return tuple(
        name for name in machine.get_corrector_names("x")
        if name not in {"ZH100R", "ZH101R"}
    )


def _paper_bpms(machine):
    return tuple(name for name in machine.bpm_names if name.startswith("MB"))


def _relative_difference(reference, candidate):
    return float(np.linalg.norm(candidate - reference) / np.linalg.norm(reference))


def main():
    lattice = build_atf_dr_lattice()
    machine = ATFDRRingCorrection(lattice)

    # 1. RF-Track placement API: input m must become its internal mm frame.
    quadrupole = _one(lattice, "QM1R.1")
    quadrupole.set_offsets(
        ALIGNMENT_TEST_M, -ALIGNMENT_TEST_M, 0.0, 0.0, 0.0, 0.0, "center"
    )
    origin_mm = np.asarray(quadrupole.get_offsets().get_origin(), dtype=float).reshape(-1)
    quadrupole.set_offsets(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, "center")

    # 2. Correction API: public physical rad must become RF-Track mrad only
    # at Corrector.set_kick, and recover unchanged at get_kick.
    corrector_name = "ZH28R"
    corrector = _one(lattice, corrector_name)
    original_kick = machine._get_corrector_kick(corrector)
    requested_kick = original_kick.copy()
    requested_kick[0] += KICK_TEST_RAD
    machine._set_corrector_kick(corrector, requested_kick)
    readback_kick = machine._get_corrector_kick(corrector)
    raw_rftrack_kick_mrad = np.asarray(
        corrector.get_kick(machine._p_over_q), dtype=float
    ).reshape(-1)
    machine._set_corrector_kick(corrector, original_kick)

    # 3. A true unit error produces a probe-dependent factor (typically 10^3).
    # The ORM must instead remain stable across a decade of small physical-rad
    # finite-difference probes.  Sample five well-spaced Kubo COD correctors;
    # this keeps the audit fast enough to be run routinely.
    all_names = _paper_x_correctors(machine)
    names = tuple(all_names[index] for index in (0, 12, 24, 36, 47))
    bpm_names = _paper_bpms(machine)
    matrices = {}
    for probe in PROBES_RAD:
        matrices[str(probe)] = machine.compute_orbit_response(
            "x", corrector_names=names, bpm_names=bpm_names, perturbation=probe
        ).matrix
    reference = matrices[str(PROBES_RAD[-1])]
    response_stability = {
        str(probe): _relative_difference(reference, matrices[str(probe)])
        for probe in PROBES_RAD[:-1]
    }

    result = {
        "purpose": "ATF DR RF-Track correction unit-convention audit",
        "coordinate_conventions": {
            "alignment_input": "m",
            "rftrack_transverse_coordinates": "mm",
            "public_corrector_kick": "rad",
            "rftrack_corrector_kick": "mrad",
            "skew_corrector_strength": "normalized integrated K1L",
        },
        "alignment_round_trip": {
            "requested_dx_dy_m": [ALIGNMENT_TEST_M, -ALIGNMENT_TEST_M],
            "rftrack_frame_origin_mm": origin_mm.tolist(),
            "expected_mm": [0.03, -0.03, 0.0],
            "max_abs_error_mm": float(np.max(np.abs(origin_mm - [0.03, -0.03, 0.0])),),
        },
        "corrector_round_trip": {
            "corrector": corrector_name,
            "requested_increment_rad": KICK_TEST_RAD,
            "public_readback_increment_rad": float(readback_kick[0] - original_kick[0]),
            "raw_rftrack_increment_mrad": float(raw_rftrack_kick_mrad[0] - original_kick[0] * 1e3),
            "public_round_trip_error_rad": float(readback_kick[0] - requested_kick[0]),
            "raw_mrad_error": float(raw_rftrack_kick_mrad[0] - requested_kick[0] * 1e3),
        },
        "horizontal_orm_probe_stability": {
            "shape": list(reference.shape),
            "reference_probe_rad": PROBES_RAD[-1],
            "relative_frobenius_difference_to_reference": response_stability,
            "response_rms_mm_per_rad": float(np.sqrt(np.mean(reference**2))),
        },
        "interpretation": (
            "All three checks must be near floating-point precision or a small "
            "nonlinear finite-difference deviation. A factor-of-1000 error in "
            "either unit boundary would fail the round-trip and probe-stability checks."
        ),
    }
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    output = ANALYSIS / "atf_dr_rftrack_unit_convention_audit.json"
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
