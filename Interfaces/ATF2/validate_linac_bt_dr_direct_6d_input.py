"""Validate the strict direct-6D entrance-bunch path offline.

The bundled input is deliberately a tiny synthetic wiring fixture.  This
test establishes that a future reconstructed IPP1L distribution is preserved
at the API boundary and can run through the same Linac -> BT -> DR path; it
does not validate an ATF measurement or predict transmission.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from Interfaces.ATF2.ATF2_LinacBTDR_RFTrack import (
    ATF2LinacBTDRRFTrack,
    load_entrance_bunch_json,
)
from Interfaces.ATF2.benchmark_linac_bt_dr_capture_proxy import (
    MODEL_ENERGY_MATCHED_CAVITY_VOLTAGE_MV,
)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    fixture = Path(__file__).with_name("linac_bt_dr_entrance_bunch_synthetic_example.json")
    entrance = load_entrance_bunch_json(fixture)
    machine = ATF2LinacBTDRRFTrack(
        cavity_voltage_mv=MODEL_ENERGY_MATCHED_CAVITY_VOLTAGE_MV,
        dr_rf_mode="equilibrium",
        handoff_mode="sad_optics_matched",
    )
    reconstructed = machine.make_entrance_bunch_from_phase_space(entrance)
    _require(
        np.array_equal(
            reconstructed.get_phase_space(), entrance.coordinates_mm_mrad_mm_c_mev_c
        ),
        "direct 6D entrance phase space was altered before Linac tracking",
    )
    result = machine.track(
        machine.make_entrance_bunch_from_phase_space(entrance),
        dr_turns=10,
        record_turn_history=True,
    )
    _require(result.input.particles == entrance.particles, "direct input particle count changed")
    _require(result.completed_dr_turns == 10, "direct 6D fixture did not complete 10 DR turns")
    _require(result.dr_after_turns.survival_fraction_from_input == 1.0, "direct 6D fixture lost in DR")
    print(json.dumps({
        "status": "passed",
        "simulation_only": True,
        "fixture": fixture.name,
        "particles": entrance.particles,
        "input_preserved_exactly": True,
        "dr_turns": result.completed_dr_turns,
        "dr_survival": result.dr_after_turns.survival_fraction_from_input,
        "limitation": "fixture is synthetic; this validates boundary/unit wiring, not measured transmission.",
    }, indent=2))


if __name__ == "__main__":
    main()
