"""Diagnose model energy matching for SAD Linac--BT injection into the DR.

This is an offline RF-Track diagnostic.  It does not recommend an ATF RF
setting: the DR cavity phase, BT time-of-flight, and injection energy must be
calibrated before an operational conclusion is possible.  Its purpose is to
separate a longitudinal model mismatch from transverse handoff effects by
finding the Linac cavity voltage whose *model* reference exit momentum equals
the DR model synchronous momentum, then tracking that reference particle for
a requested multi-turn horizon.
"""

from __future__ import annotations

import argparse
import json

from Interfaces.ATF2.ATF2_LinacBTDR_RFTrack import ATF2LinacBTDRRFTrack


def reference_exit_momentum(voltage_mv: float) -> tuple[ATF2LinacBTDRRFTrack, float]:
    machine = ATF2LinacBTDRRFTrack(
        cavity_voltage_mv=voltage_mv,
        handoff_mode="sad_optics_matched",
        dr_rf_mode="equilibrium",
    )
    result = machine.track(machine.make_reference_bunch(), dr_turns=0)
    return machine, result.linac_bt_exit.mean_p_mev_c


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first-voltage-mv", type=float, default=76.14)
    parser.add_argument("--second-voltage-mv", type=float, default=75.86)
    parser.add_argument("--turns", type=int, default=100)
    args = parser.parse_args()
    if args.turns < 1:
        raise ValueError("turns must be positive")
    if args.first_voltage_mv == args.second_voltage_mv:
        raise ValueError("the two voltage probes must differ")

    first_machine, first_p = reference_exit_momentum(args.first_voltage_mv)
    second_machine, second_p = reference_exit_momentum(args.second_voltage_mv)
    target_p = float(second_machine.dr_synchronous_orbit[5])
    slope = (second_p - first_p) / (args.second_voltage_mv - args.first_voltage_mv)
    if slope == 0.0:
        raise RuntimeError("Linac exit momentum did not change between RF voltage probes")
    matched_voltage = args.second_voltage_mv + (target_p - second_p) / slope

    matched_machine, matched_p = reference_exit_momentum(matched_voltage)
    tracked = matched_machine.track(
        matched_machine.make_reference_bunch(),
        dr_turns=args.turns,
        record_turn_history=True,
    )
    survival_history = [
        item.survival_fraction_from_input for item in tracked.dr_turn_history
    ]
    first_loss_turn = next(
        (index + 1 for index, value in enumerate(survival_history) if value < 1.0),
        None,
    )
    print(json.dumps({
        "simulation_only": True,
        "purpose": "model RF-energy/capture diagnostic; not an operational RF setting",
        "probes": [
            {"cavity_voltage_mv": args.first_voltage_mv, "linac_bt_exit_p_mev_c": first_p},
            {"cavity_voltage_mv": args.second_voltage_mv, "linac_bt_exit_p_mev_c": second_p},
        ],
        "dr_synchronous_p_mev_c": target_p,
        "secant_slope_mev_c_per_mv": slope,
        "model_matched_cavity_voltage_mv": matched_voltage,
        "model_matched_exit_p_mev_c": matched_p,
        "model_matched_delta_p_mev_c": matched_p - target_p,
        "turns_requested": args.turns,
        "turns_completed": tracked.completed_dr_turns,
        "ring_survival_after_turns": tracked.dr_after_turns.survival_fraction_from_input,
        "first_loss_turn": first_loss_turn,
        "limitations": [
            "No surveyed IPZT-to-RING0 coordinate/pulsed-kicker calibration is loaded.",
            "No complete measured aperture/loss-monitor table is loaded.",
            "DR RF phase and BT time-of-flight are model values, not calibrated timing data.",
        ],
    }, indent=2))


if __name__ == "__main__":
    main()
