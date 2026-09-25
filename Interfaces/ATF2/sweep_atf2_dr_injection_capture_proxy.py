"""Compare DR-injection endpoint orbit with prompt multi-turn loss.

The sweep applies a synthetic offset to the *BT-to-DR handoff*, not to a real
corrector or control system.  It therefore probes the sensitivity that remains
while the IPZT-to-RING0 survey/kicker map is uncalibrated.  The default uses
only the historically located KIX.1/.2 5-mm half-aperture screen;
``--aperture-screen none`` provides the corresponding dynamic-capture
comparison.  Neither is a complete ATF aperture model or a final-transmission
prediction.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import replace

import numpy as np

from Interfaces.ATF2.ATF2_LinacBTDR_RFTrack import ATF2LinacBTDRRFTrack
from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_lattice import (
    HISTORICAL_ATF_DR_APERTURE_SOURCE,
    get_historical_extraction_kicker_apertures,
)
from Interfaces.ATF2.benchmark_linac_bt_dr_capture_proxy import (
    MODEL_ENERGY_MATCHED_CAVITY_VOLTAGE_MV,
)


def _offsets(value: str) -> list[float]:
    values = [float(item) for item in value.split(",")]
    if not values or not all(np.isfinite(values)):
        raise ValueError("offsets-mm must be a comma-separated finite number list")
    return values


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--offsets-mm",
        default="0,0.1,0.25,0.5,0.75,1,2,3",
        help="Synthetic horizontal handoff offsets in mm.",
    )
    parser.add_argument("--turns", type=int, default=10)
    parser.add_argument(
        "--aperture-screen",
        choices=("historical-kix", "none"),
        default="historical-kix",
        help=(
            "historical-kix: only KIX.1/.2 5-mm half apertures; none: no "
            "explicit aperture, for separating dynamic prompt loss from this "
            "partial aperture screen."
        ),
    )
    parser.add_argument(
        "--cavity-voltage-mv",
        type=float,
        default=MODEL_ENERGY_MATCHED_CAVITY_VOLTAGE_MV,
    )
    args = parser.parse_args()
    if args.turns < 1:
        raise ValueError("turns must be positive")
    offsets = _offsets(args.offsets_mm)
    apertures = (
        get_historical_extraction_kicker_apertures()
        if args.aperture_screen == "historical-kix" else {}
    )
    machine = ATF2LinacBTDRRFTrack(
        cavity_voltage_mv=args.cavity_voltage_mv,
        dr_rf_mode="equilibrium",
        handoff_mode="sad_optics_matched",
        dr_apertures_mm=apertures,
    )
    baseline_handoff = machine.handoff
    results = []
    try:
        for offset in offsets:
            machine.handoff = replace(
                baseline_handoff,
                offset_mm_mrad=(
                    baseline_handoff.offset_mm_mrad
                    + np.array((offset, 0.0, 0.0, 0.0), dtype=float)
                ),
            )
            endpoint = machine.track(machine.make_reference_bunch(), dr_turns=0)
            tracked = machine.track(
                machine.make_reference_bunch(),
                dr_turns=args.turns,
                record_turn_history=True,
            )
            injection_orbit = np.array((
                endpoint.dr_injection.mean_x_mm,
                endpoint.dr_injection.mean_xp_mrad,
                endpoint.dr_injection.mean_y_mm,
                endpoint.dr_injection.mean_yp_mrad,
            )) - machine.dr_closed_orbit
            survival_history = [
                item.survival_fraction_from_input for item in tracked.dr_turn_history
            ]
            results.append({
                "synthetic_handoff_dx_mm": offset,
                "endpoint_injection_orbit_mm_mrad": injection_orbit.tolist(),
                "endpoint_orbit_norm_mixed_units": float(np.linalg.norm(injection_orbit)),
                "completed_dr_turns": tracked.completed_dr_turns,
                "survival_after_turns": tracked.dr_after_turns.survival_fraction_from_input,
                "turn_history_survival": survival_history,
                "first_recorded_loss_turn": next(
                    (index + 1 for index, value in enumerate(survival_history) if value < 1.0),
                    None,
                ),
            })
    finally:
        machine.handoff = baseline_handoff
    print(json.dumps({
        "simulation_only": True,
        "purpose": "endpoint-orbit versus prompt-loss sensitivity; no real controls",
        "synthetic_error": (
            "horizontal offset added to the design-optics IPZT-to-RING0 handoff; "
            "it represents uncalibrated injection-map/kicker error, not a corrector setting"
        ),
        "cavity_voltage_mv": args.cavity_voltage_mv,
        "turns_per_case": args.turns,
        "partial_aperture_screen": {
            "mode": args.aperture_screen,
            "configured": apertures,
            "source": HISTORICAL_ATF_DR_APERTURE_SOURCE if apertures else None,
            "limitation": (
                "Only KIX.1/.2 are mapped.  This does not include arc, wiggler-mask, "
                "or south-straight constraints; circular KIX geometry is an explicit "
                "RF-Track study assumption."
                if apertures else
                "No explicit physical aperture is configured; this is a dynamic-capture comparison."
            ),
        },
        "cases": results,
        "interpretation": {
            "endpoint_proxy": (
                "A large injection centroid error can reject candidates before costly "
                "multi-turn tracking, but cannot replace a complete aperture map."
            ),
            "loss_attribution": (
                "Compare aperture-screen=historical-kix and none before attributing a "
                "loss to the partial KIX screen.  A loss in both modes is a model "
                "dynamic-capture result, not an aperture measurement."
            ),
            "safety": (
                "This script deliberately perturbs the handoff rather than applying "
                "unbounded native corrector strengths to RF-Track."
            ),
        },
    }, indent=2))


if __name__ == "__main__":
    main()
