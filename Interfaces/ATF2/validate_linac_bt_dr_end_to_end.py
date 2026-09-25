"""Fail-fast integration smoke test for the SAD Linac -> BT -> DR pipeline.

This is deliberately a *software/model-consistency* gate.  It verifies the
configured 1.3-GeV RF energy match, reproducible finite-bunch generation,
SAD-seeded/RF-Track-propagated optics handoff, short finite-bunch survival,
and sparse long-turn reference tracking.  It is not a prediction of measured
ATF final transmission, for which survey, timing, and complete aperture data
are still required.
"""

from __future__ import annotations

import argparse
import json

import numpy as np

from Interfaces.ATF2.ATF2_LinacBTDR_RFTrack import (
    ATF2LinacBTDRRFTrack,
    EntranceBunchTwiss,
)
from Interfaces.ATF2.benchmark_linac_bt_dr_capture_proxy import (
    MODEL_ENERGY_MATCHED_CAVITY_VOLTAGE_MV,
)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--long-turns", type=int, default=100)
    parser.add_argument("--long-sample-every", type=int, default=20)
    args = parser.parse_args()
    if args.long_turns < 1 or args.long_sample_every < 1:
        raise ValueError("long-turns and long-sample-every must be positive")

    machine = ATF2LinacBTDRRFTrack(
        cavity_voltage_mv=MODEL_ENERGY_MATCHED_CAVITY_VOLTAGE_MV,
        dr_rf_mode="equilibrium",
        handoff_mode="sad_optics_matched",
    )
    energy_error = (
        machine.reference_exit_momentum_mev_c - float(machine.dr_synchronous_orbit[5])
    )
    _require(abs(energy_error) < 1.0e-4, f"DR energy mismatch too large: {energy_error} MeV/c")

    twiss = EntranceBunchTwiss(
        emitt_x_norm_mm_mrad=1.0,
        emitt_y_norm_mm_mrad=1.0,
        beta_x_m=1.93,
        beta_y_m=1.93,
    )
    generated_a = np.asarray(
        machine.make_entrance_bunch(twiss, particles=16).get_phase_space(), dtype=float
    )
    generated_b = np.asarray(
        machine.make_entrance_bunch(twiss, particles=16).get_phase_space(), dtype=float
    )
    _require(
        np.array_equal(generated_a, generated_b),
        "finite entrance bunch generation is not reproducible in this RF-Track environment",
    )
    finite = machine.track(
        machine.make_entrance_bunch(twiss, particles=16),
        dr_turns=10,
        record_turn_history=True,
    )
    _require(finite.linac_bt_exit.survival_fraction_from_input == 1.0, "Linac+BT lost design test bunch")
    _require(finite.dr_after_turns.survival_fraction_from_input == 1.0, "DR lost design test bunch")
    _require(
        finite.dr_injection_optics.x.mismatch_to_design < 1.2
        and finite.dr_injection_optics.y.mismatch_to_design < 1.2,
        "configured Linac+BT handoff no longer matches the DR design optics",
    )
    reference = machine.track(
        machine.make_reference_bunch(),
        dr_turns=args.long_turns,
        record_turn_history=True,
        turn_history_sample_every=args.long_sample_every,
    )
    _require(reference.completed_dr_turns == args.long_turns, "reference did not complete long horizon")
    _require(reference.dr_after_turns.survival_fraction_from_input == 1.0, "reference lost in long horizon")
    _require(
        reference.dr_turn_history_turns[-1] == args.long_turns,
        "long-horizon history did not retain its final turn",
    )
    print(json.dumps({
        "status": "passed",
        "simulation_only": True,
        "energy_error_mev_c": energy_error,
        "finite_bunch": {
            "particles": finite.input.particles,
            "linac_bt_survival": finite.linac_bt_exit.survival_fraction_from_input,
            "dr_survival_after_10_turns": finite.dr_after_turns.survival_fraction_from_input,
            "mismatch_x": finite.dr_injection_optics.x.mismatch_to_design,
            "mismatch_y": finite.dr_injection_optics.y.mismatch_to_design,
            "reproducible_entrance_generation": True,
        },
        "reference_long_horizon": {
            "turns": reference.completed_dr_turns,
            "survival": reference.dr_after_turns.survival_fraction_from_input,
            "recorded_turns": list(reference.dr_turn_history_turns),
        },
        "limitations": [
            "No surveyed IPZT-to-RING0 handoff or pulsed-kicker calibration.",
            "No complete physical aperture/loss map.",
            "Finite-bunch Twiss/emittance is an explicit study input, not measured ATF data.",
        ],
    }, indent=2))


if __name__ == "__main__":
    main()
