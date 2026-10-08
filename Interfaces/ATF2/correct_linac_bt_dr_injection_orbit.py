"""Offline response-matrix correction baseline for Linac--BT--DR injection.

The script creates a deliberately synthetic corrector error, reconstructs the
linear response from four independent correctors to the DR injection orbit,
then applies an SVD least-squares correction and checks a multi-turn reference
particle.  It never connects to controls.  This is the response-matrix
baseline for a future RL comparison, not a statement about ATF machine errors
or actuator limits.
"""

from __future__ import annotations

import argparse
import json

import numpy as np

from Interfaces.ATF2.ATF2_LinacBTDR_RFTrack import ATF2LinacBTDRRFTrack


ACTUATORS = ("ZH1L", "ZH2L", "ZV1L", "ZV2L")
SYNTHETIC_ERROR = {"ZH5L": 3.0e-4, "ZV5L": -2.0e-4}
MODEL_MATCHED_CAVITY_VOLTAGE_MV = 75.85997836492652


def injection_orbit(machine: ATF2LinacBTDRRFTrack) -> np.ndarray:
    result = machine.track(machine.make_reference_bunch(), dr_turns=0)
    return np.array((
        result.dr_injection.mean_x_mm,
        result.dr_injection.mean_xp_mrad,
        result.dr_injection.mean_y_mm,
        result.dr_injection.mean_yp_mrad,
    ))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--turns", type=int, default=100)
    parser.add_argument("--response-step", type=float, default=1.0e-4)
    parser.add_argument("--svd-rcond", type=float, default=1.0e-8)
    args = parser.parse_args()
    if args.turns < 1:
        raise ValueError("turns must be positive")

    machine = ATF2LinacBTDRRFTrack(
        cavity_voltage_mv=MODEL_MATCHED_CAVITY_VOLTAGE_MV,
        handoff_mode="sad_optics_matched",
        dr_rf_mode="equilibrium",
    )
    initial_actuators = machine.get_linac_bt_correctors(ACTUATORS)
    response, names = machine.linac_bt_dr_orbit_response_matrix(
        ACTUATORS, step=args.response_step
    )
    machine.set_linac_bt_correctors(SYNTHETIC_ERROR)
    before = injection_orbit(machine) - machine.dr_closed_orbit
    delta, _, rank, singular_values = np.linalg.lstsq(
        response, -before, rcond=args.svd_rcond
    )
    machine.set_linac_bt_correctors({
        name: initial_actuators[name] + value
        for name, value in zip(names, delta)
    })
    tracked = machine.track(
        machine.make_reference_bunch(),
        dr_turns=args.turns,
        record_turn_history=True,
    )
    after = np.array((
        tracked.dr_injection.mean_x_mm,
        tracked.dr_injection.mean_xp_mrad,
        tracked.dr_injection.mean_y_mm,
        tracked.dr_injection.mean_yp_mrad,
    )) - machine.dr_closed_orbit
    print(json.dumps({
        "simulation_only": True,
        "purpose": "synthetic response-matrix correction baseline for RL comparison",
        "cavity_voltage_mv": MODEL_MATCHED_CAVITY_VOLTAGE_MV,
        "synthetic_error_rftrack_strength": SYNTHETIC_ERROR,
        "actuators": list(names),
        "response_step_rftrack_strength": args.response_step,
        "response_matrix_mm_mrad_per_strength": response.tolist(),
        "singular_values": singular_values.tolist(),
        "numerical_rank": int(rank),
        "orbit_error_before_mm_mrad": before.tolist(),
        "correction_delta_rftrack_strength": delta.tolist(),
        "orbit_error_after_mm_mrad": after.tolist(),
        "turns_requested": args.turns,
        "turns_completed": tracked.completed_dr_turns,
        "ring_survival_after_turns": tracked.dr_after_turns.survival_fraction_from_input,
        "limitations": [
            "Synthetic errors and uncalibrated RFTrack corrector strengths are used.",
            "No surveyed IPZT-to-RING0 coordinate/pulsed-kicker calibration is loaded.",
            "No complete measured aperture/loss-monitor table is loaded.",
        ],
    }, indent=2))


if __name__ == "__main__":
    main()
