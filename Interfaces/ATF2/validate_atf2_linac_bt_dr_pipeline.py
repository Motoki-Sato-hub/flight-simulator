"""Run a small deterministic Linac->BT->DR pipeline validation.

This is a software-integration and timing check, not a physical transmission
prediction.  The default four-particle numerical probe has deliberately tiny
transverse offsets and no time/energy spread.  Supply a validated injected
6D distribution through :class:`ATF2LinacBTDRRFTrack` for physics studies.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import RF_Track as rft

from Interfaces.ATF2.ATF2_LinacBTDR_RFTrack import ATF2LinacBTDRRFTrack


def numerical_probe(machine: ATF2LinacBTDRRFTrack, particles: int):
    if particles < 4:
        raise ValueError("particles must be at least four for non-singular projected optics")
    # Symmetric, non-collinear offsets avoid changing the reference centroid
    # and make each projected 2D covariance non-singular for >=4 particles.
    # RF phase and momentum are held fixed so this remains an integration test.
    angle = 2.0 * np.pi * np.arange(particles, dtype=float) / particles
    phase_space = np.zeros((particles, 6), dtype=float)
    phase_space[:, 0] = 0.002 * np.cos(angle)
    phase_space[:, 1] = 0.001 * np.sin(angle)
    phase_space[:, 2] = 0.0015 * np.cos(3.0 * angle)
    phase_space[:, 3] = 0.0005 * np.sin(3.0 * angle)
    phase_space[:, 5] = machine.input_momentum_mev_c
    return rft.Bunch6d(rft.electronmass, 1.0e9, -1.0, phase_space)


def dispersion_handoff_check(machine: ATF2LinacBTDRRFTrack):
    """Verify the declared first-order dispersion conversion, if selected."""
    handoff = machine.handoff
    if handoff.source_dispersion_mm_mrad is None:
        return None
    delta = 0.01
    phase_space = np.zeros((1, 6), dtype=float)
    phase_space[0, :4] = (
        machine.reference_exit
        + handoff.source_dispersion_mm_mrad * delta
    )
    phase_space[0, 5] = handoff.reference_momentum_mev_c * (1.0 + delta)
    mapped = handoff.apply_phase_space(phase_space)[0, :4]
    expected = machine.dr_closed_orbit + handoff.target_dispersion_mm_mrad * delta
    return {
        "delta": delta,
        "max_abs_error_mm_or_mrad": float(np.max(np.abs(mapped - expected))),
        "expected_target_coordinates_mm_mrad": expected.tolist(),
    }


def load_apertures(path: str | None):
    """Load an explicit, user-owned aperture table for a transmission study.

    The JSON schema is ``{"linac_bt": {element: [x_mm, y_mm, shape?]},
    "dr": {...}}``.  Omitting a section is allowed; unknown sections catch
    spelling mistakes rather than being silently ignored.
    """
    if path is None:
        return {}, {}
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or set(payload) - {"linac_bt", "dr"}:
        raise ValueError("aperture JSON must contain only optional 'linac_bt' and 'dr' objects")
    linac_bt = payload.get("linac_bt", {})
    dr = payload.get("dr", {})
    if not isinstance(linac_bt, dict) or not isinstance(dr, dict):
        raise ValueError("each aperture JSON section must be an object keyed by element name")
    return linac_bt, dr


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--particles", type=int, default=4)
    parser.add_argument("--turns", type=int, default=20)
    parser.add_argument(
        "--dr-rf-mode",
        choices=("disabled", "equilibrium"),
        default="disabled",
        help=(
            "disabled: transverse ring only; equilibrium: deterministic RF/radiation "
            "ring with its model synchronous phase (not calibrated ATF injection)."
        ),
    )
    parser.add_argument(
        "--handoff-mode",
        choices=("reference_anchored", "sad_optics_matched"),
        default="reference_anchored",
        help=(
            "reference_anchored: identity plus reference offset; "
            "sad_optics_matched: zero-phase Twiss/dispersion design baseline."
        ),
    )
    parser.add_argument(
        "--aperture-json",
        help=(
            "JSON table {'linac_bt': {'IPZT': [x_mm, y_mm]}, "
            "'dr': {'QF2R.1': [x_mm, y_mm, 'rectangular']}}; omitted means no aperture model."
        ),
    )
    args = parser.parse_args()

    linac_bt_apertures, dr_apertures = load_apertures(args.aperture_json)
    machine = ATF2LinacBTDRRFTrack(
        dr_rf_mode=args.dr_rf_mode,
        handoff_mode=args.handoff_mode,
        linac_bt_apertures_mm=linac_bt_apertures,
        dr_apertures_mm=dr_apertures,
    )
    probe = numerical_probe(machine, args.particles)
    start = time.perf_counter()
    endpoint_only = machine.track(probe, dr_turns=0)
    endpoint_seconds = time.perf_counter() - start
    start = time.perf_counter()
    tracked = machine.track(probe, dr_turns=args.turns, record_turn_history=True)
    full_seconds = time.perf_counter() - start
    injection = tracked.dr_injection
    centroid_mismatch = np.array(
        [injection.mean_x_mm, injection.mean_xp_mrad,
         injection.mean_y_mm, injection.mean_yp_mrad]
    ) - machine.dr_closed_orbit
    payload = {
        "simulation_only": True,
        "probe": "symmetric numerical 6D probe; not an ATF measured bunch",
        "particles": args.particles,
        "requested_dr_turns": args.turns,
        "dr_rf_mode": args.dr_rf_mode,
        "handoff_mode": args.handoff_mode,
        "handoff": machine.handoff.as_dict(),
        "aperture_json": args.aperture_json,
        "configured_apertures_mm": {
            "linac_bt": machine.linac_bt_apertures,
            "dr": machine.dr_apertures,
        },
        "dispersion_handoff_check": dispersion_handoff_check(machine),
        "dr_synchronous_orbit": (
            None if machine.dr_synchronous_orbit is None
            else machine.dr_synchronous_orbit.tolist()
        ),
        "linac_bt_exit_minus_dr_synchronous_p_mev_c": (
            None if machine.dr_synchronous_orbit is None
            else tracked.linac_bt_exit.mean_p_mev_c - machine.dr_synchronous_orbit[5]
        ),
        "timing_seconds": {
            "linac_bt_and_handoff": endpoint_seconds,
            "with_dr_turns": full_seconds,
            "incremental_dr_per_turn": (
                (full_seconds - endpoint_seconds) / args.turns if args.turns else 0.0
            ),
        },
        "linac_bt_survival": tracked.linac_bt_exit.survival_fraction_from_input,
        "ring_survival_after_turns": tracked.dr_after_turns.survival_fraction_from_input,
        "dr_injection_optics": tracked.as_dict()["dr_injection_optics"],
        "injection_centroid_minus_dr_closed_orbit_mm_mrad": centroid_mismatch.tolist(),
        "turn_history_survival": [
            summary.survival_fraction_from_input for summary in tracked.dr_turn_history
        ],
        "turn_history_turns": list(tracked.dr_turn_history_turns),
        "interpretation": {
            "endpoint_proxy": (
                "DR injection centroid, dispersion-subtracted projected Twiss, and "
                "mismatch are cheap diagnostics, but cannot prove capture or final transmission."
            ),
            "short_tracking": (
                "Use turn history to expose immediate geometric/dynamic loss once an "
                "injection map and apertures are loaded."
            ),
            "final_transmission": (
                "Not established by this run: the current DR model lacks a surveyed "
                "injection map and aperture/loss monitors"
                + (
                    "; the transverse-only run also lacks longitudinal capture."
                    if args.dr_rf_mode == "disabled" else
                    "; RF phase is only aligned to the model synchronous orbit, not "
                    "to calibrated BT time-of-flight and injection energy."
                )
            ),
        },
    }
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
