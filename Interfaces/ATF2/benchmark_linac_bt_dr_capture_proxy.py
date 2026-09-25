"""Benchmark the staged Linac -> BT -> DR capture evaluation.

This executable uses the SAD ``IPP1L`` design Twiss as an entrance *optics*
seed, propagates its endpoint optics numerically through the configured
RF-Track Linac+BT lattice, and then uses the resulting design-optics handoff
to the DR.  The supplied normalised emittance and particle count are a stated
simulation study choice; they are not a substitute for a measured ATF bunch.

The report deliberately separates cheap endpoint diagnostics from a short
multi-turn survival check.  It does not claim final physical transmission:
that needs a calibrated IPZT-to-ring injection map and aperture/loss data.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from Interfaces.ATF2.ATF2_LinacBTDR_RFTrack import (
    ATF2LinacBTDRRFTrack,
    EntranceBunchTwiss,
    load_entrance_bunch_json,
)
from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_lattice import (
    HISTORICAL_ATF_DR_APERTURE_SOURCE,
    get_historical_extraction_kicker_apertures,
)
from Interfaces.ATF2.simulate_linac_bt_dr_handoff import TransverseHandoff


# Offline RF voltage which brings the current SAD-derived reference exit to
# the deterministic RF-Track DR synchronous momentum.  It is model tuning,
# not an operational ATF voltage recommendation; RF phase and BT time of
# flight still require calibration before any real-machine interpretation.
MODEL_ENERGY_MATCHED_CAVITY_VOLTAGE_MV = 75.85997836493
SAD_IPP1L_BETA_M = 1.93


def _summary(result):
    return {
        "linac_bt_survival": result.linac_bt_exit.survival_fraction_from_input,
        "dr_survival": result.dr_after_turns.survival_fraction_from_input,
        "completed_dr_turns": result.completed_dr_turns,
        "dr_injection_mismatch": {
            "x": result.dr_injection_optics.x.mismatch_to_design,
            "y": result.dr_injection_optics.y.mismatch_to_design,
        },
        "dr_injection_rms_mm": {
            "x": result.dr_injection.rms_x_mm,
            "y": result.dr_injection.rms_y_mm,
        },
    }


def _load_handoff(path: str | None) -> TransverseHandoff | None:
    """Load either a handoff payload or a handoff-fit result JSON."""
    if path is None:
        return None
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("handoff JSON must be an object")
    handoff_payload = payload.get("handoff", payload)
    if not isinstance(handoff_payload, dict):
        raise ValueError("handoff JSON 'handoff' field must be an object")
    return TransverseHandoff.from_dict(handoff_payload)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--particles", type=int, default=16)
    parser.add_argument("--short-turns", type=int, default=10)
    parser.add_argument(
        "--long-turns",
        type=int,
        default=0,
        help=(
            "Optional direct long-horizon validation from the same entrance "
            "bunch.  Zero skips it; it is intentionally not run for every "
            "endpoint candidate."
        ),
    )
    parser.add_argument(
        "--turn-history-sample-every",
        type=int,
        default=1,
        help=(
            "Keep every Nth DR-turn summary during the tracking stage.  This "
            "bounds history size for long validation but the reported minimum "
            "is then a sampled, not exact, minimum."
        ),
    )
    parser.add_argument(
        "--long-turn-history-sample-every",
        type=int,
        default=100,
        help="Retain every Nth long-validation turn summary (also retains final/loss turn).",
    )
    parser.add_argument(
        "--normalised-emittance-mm-mrad",
        type=float,
        default=1.0,
        help="Explicit study value for both planes; not a measured-bunch default.",
    )
    parser.add_argument(
        "--entrance-bunch-json",
        help=(
            "Validated direct 6D bunch at SAD IPP1L.  When supplied, it replaces "
            "--particles and --normalised-emittance-mm-mrad, preserving measured "
            "correlations/tails for each endpoint and tracking comparison."
        ),
    )
    parser.add_argument(
        "--cavity-voltage-mv",
        type=float,
        default=MODEL_ENERGY_MATCHED_CAVITY_VOLTAGE_MV,
    )
    parser.add_argument(
        "--historical-kix-screen",
        action="store_true",
        help=(
            "Apply the historically reported 5-mm half-aperture at the two "
            "unambiguously located KIX extraction-kicker elements.  This is "
            "a partial loss screen, not a complete DR aperture model."
        ),
    )
    parser.add_argument(
        "--handoff-json",
        help=(
            "A TransverseHandoff.as_dict() payload or linac_bt_dr_handoff_fit "
            "result JSON.  It replaces the design handoff without any controls access."
        ),
    )
    args = parser.parse_args()
    if args.particles < 2 and args.entrance_bunch_json is None:
        raise ValueError("particles must be at least two for projected optics")
    if args.short_turns < 0:
        raise ValueError("short-turns must be non-negative")
    if args.long_turns < 0:
        raise ValueError("long-turns must be non-negative")
    if args.turn_history_sample_every < 1:
        raise ValueError("turn-history-sample-every must be positive")
    if args.long_turn_history_sample_every < 1:
        raise ValueError("long-turn-history-sample-every must be positive")
    if args.normalised_emittance_mm_mrad <= 0.0:
        raise ValueError("normalised emittance must be positive")

    started = time.perf_counter()
    historical_kix_apertures = (
        get_historical_extraction_kicker_apertures()
        if args.historical_kix_screen else {}
    )
    fitted_handoff = _load_handoff(args.handoff_json)
    machine = ATF2LinacBTDRRFTrack(
        cavity_voltage_mv=args.cavity_voltage_mv,
        dr_rf_mode="equilibrium",
        handoff_mode="sad_optics_matched",
        handoff=fitted_handoff,
        dr_apertures_mm=historical_kix_apertures,
    )
    construction_seconds = time.perf_counter() - started
    direct_entrance = (
        None if args.entrance_bunch_json is None
        else load_entrance_bunch_json(args.entrance_bunch_json)
    )
    if direct_entrance is not None and direct_entrance.particles < 2:
        raise ValueError("direct entrance bunch must contain at least two particles for projected optics")
    if direct_entrance is None:
        twiss = EntranceBunchTwiss(
            emitt_x_norm_mm_mrad=args.normalised_emittance_mm_mrad,
            emitt_y_norm_mm_mrad=args.normalised_emittance_mm_mrad,
            beta_x_m=SAD_IPP1L_BETA_M,
            beta_y_m=SAD_IPP1L_BETA_M,
        )

        def make_entrance_bunch():
            return machine.make_entrance_bunch(twiss, particles=args.particles)

        entrance_distribution = {
            "kind": "synthetic_twiss_study",
            "optics_seed": "SAD combined-lattice IPP1L: beta_x=beta_y=1.93 m, alpha=0",
            "normalised_emittance_mm_mrad": args.normalised_emittance_mm_mrad,
            "particles": args.particles,
            "warning": "emittance and distribution are explicit study inputs, not measured ATF bunch data",
        }
    else:
        def make_entrance_bunch():
            return machine.make_entrance_bunch_from_phase_space(direct_entrance)

        entrance_distribution = {
            "kind": "direct_6d_ippl_input",
            "source_json": str(Path(args.entrance_bunch_json)),
            "location": direct_entrance.location,
            "coordinate_order": "x_mm,xp_mrad,y_mm,yp_mrad,t_mm_c,p_mev_c",
            "particles": direct_entrance.particles,
            "charge_e": direct_entrance.charge_e,
            "provenance": direct_entrance.provenance,
            "warning": "input is validated for units/boundary only; diagnostic reconstruction remains an external calibration.",
        }

    started = time.perf_counter()
    endpoint = machine.track(make_entrance_bunch(), dr_turns=0)
    endpoint_seconds = time.perf_counter() - started
    started = time.perf_counter()
    short = machine.track(
        make_entrance_bunch(),
        dr_turns=args.short_turns,
        record_turn_history=True,
        turn_history_sample_every=args.turn_history_sample_every,
    )
    short_seconds = time.perf_counter() - started
    long = None
    long_seconds = 0.0
    if args.long_turns:
        started = time.perf_counter()
        long = machine.track(
            make_entrance_bunch(),
            dr_turns=args.long_turns,
            record_turn_history=True,
            turn_history_sample_every=args.long_turn_history_sample_every,
        )
        long_seconds = time.perf_counter() - started
    min_turn_survival = min(
        (item.survival_fraction_from_input for item in short.dr_turn_history),
        default=short.dr_after_turns.survival_fraction_from_input,
    )
    print(json.dumps({
        "simulation_only": True,
        "entrance_distribution": entrance_distribution,
        "configured_model": {
            "linac_input_momentum_mev_c": machine.input_momentum_mev_c,
            "cavity_voltage_mv": args.cavity_voltage_mv,
            "linac_bt_reference_exit_momentum_mev_c": machine.reference_exit_momentum_mev_c,
            "dr_synchronous_momentum_mev_c": float(machine.dr_synchronous_orbit[5]),
            "reference_energy_difference_mev_c": (
                machine.reference_exit_momentum_mev_c - float(machine.dr_synchronous_orbit[5])
            ),
            "handoff_provenance": machine.handoff.provenance,
            "handoff_json": args.handoff_json,
            "linac_bt_energy_profile_magnet_adjustments": (
                machine.linac_bt_metadata.get("energy_profile_magnet_adjustments")
            ),
            "linac_bt_scaling_assumption": (
                "SAD normalised quadrupole/bend strengths are evaluated at the "
                "configured local reference momentum (ideal 1.3-GeV field scaling), "
                "not loaded from measured magnet-current settings."
            ),
        },
        "aperture_screen": {
            "enabled": bool(historical_kix_apertures),
            "configured_dr_half_apertures_mm": machine.dr_apertures,
            "source": (
                HISTORICAL_ATF_DR_APERTURE_SOURCE
                if historical_kix_apertures else None
            ),
            "scope": (
                "Only KIX.1/.2 are mapped; arc, wiggler-mask, and south-straight "
                "restrictions are deliberately not inferred.  The KIX circular "
                "cross-section is an explicit RF-Track study assumption, not sourced "
                "engineering geometry."
                if historical_kix_apertures else
                "No aperture/loss screen selected."
            ),
        },
        "endpoint_proxy": _summary(endpoint),
        "short_tracking": {
            **_summary(short),
            "requested_dr_turns": args.short_turns,
            "turn_history_sample_every": short.turn_history_sample_every,
            "minimum_recorded_turn_survival": min_turn_survival,
            "minimum_recorded_turn_survival_is_exact": (
                args.turn_history_sample_every == 1
            ),
        },
        "long_tracking": (
            None if long is None else {
                **_summary(long),
                "requested_dr_turns": args.long_turns,
                "turn_history_sample_every": long.turn_history_sample_every,
                "recorded_turns": list(long.dr_turn_history_turns),
                "recorded_survival": [
                    item.survival_fraction_from_input for item in long.dr_turn_history
                ],
                "first_recorded_loss_turn": next(
                    (
                        turn for turn, item in zip(
                            long.dr_turn_history_turns, long.dr_turn_history
                        ) if item.survival_fraction_from_input < 1.0
                    ),
                    None,
                ),
                "first_recorded_loss_turn_is_exact": (
                    args.long_turn_history_sample_every == 1
                ),
            }
        ),
        "timing_seconds": {
            "construction_including_configured_linearisation": construction_seconds,
            "endpoint_only": endpoint_seconds,
            "short_tracking_total": short_seconds,
            "incremental_short_turn": (
                (short_seconds - endpoint_seconds) / args.short_turns
                if args.short_turns else 0.0
            ),
            "long_tracking_total": long_seconds if long is not None else None,
            "incremental_long_turn": (
                (long_seconds - endpoint_seconds) / args.long_turns
                if long is not None and args.long_turns else None
            ),
        },
        "interpretation": {
            "fast_filter": "Rank corrector/optics candidates by energy error and endpoint mismatch before ring tracking.",
            "validation": "Use short tracking to reject prompt loss, then reserve the optional direct long-turn stage for selected candidates.",
            "limit": "Neither quantity is final physical transmission until the IPZT-to-RING0 map, injection pulse, and complete aperture/loss model are calibrated.",
        },
    }, indent=2))


if __name__ == "__main__":
    main()
