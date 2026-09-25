"""Find an offline 80-MeV Linac RF phase setting compatible with DR injection.

The model has sixteen powered SAD structures and eight online-style phase
setpoints.  The explicitly provisional mapping is L1 -> CA1L/CA2L, ...,
L8 -> CA15L/CA16L.  Phase offsets are relative to each structure's SAD phase
(-1.67 degrees from accelerating crest), not absolute hardware phase values.

The default common powered-structure voltage was matched once by direct
RF-Track reference tracking to the DR synchronous momentum for an 80-MeV
entrance.  The resulting Linac+BT magnet fields and the IPZT-to-RING0 handoff
are then held fixed while the eight RF phases are scanned.  This is
intentional: re-matching the optics handoff or re-tuning magnets for every
trial would mask the effect of RF phase errors.

Everything is simulation-only.  It is a calibration-design exercise, not an
ATF RF setpoint recommendation: the klystron-to-structure mapping, amplitude
calibration, timing, and IPZT-to-RING0 injection map remain to be measured.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from typing import Mapping

import numpy as np

from Interfaces.ATF2.ATF2_LinacBTDR_RFTrack import (
    ATF2LinacBTDRRFTrack,
    EntranceBunchTwiss,
    KLYSTRON_CAVITY_GROUPS,
)


INPUT_MOMENTUM_MEV_C = 80.0
# Direct RF-Track model match at 80 MeV/c with all klystron offsets zero.
# It gives a reference exit mismatch of 1.35e-5 MeV/c in the current model.
# This is not an operational amplitude setpoint.
MODEL_MATCHED_80MEV_CAVITY_VOLTAGE_MV = 76.007891
STUDY_TWISS = EntranceBunchTwiss(1.0, 1.0, 1.93, 1.93)
SYNTHETIC_STARTING_PHASE_ERROR_DEG = {
    "L1": 2.0, "L2": -1.5, "L3": 1.0, "L4": -2.0,
    "L5": 1.5, "L6": -1.0, "L7": 2.0, "L8": -1.5,
}


@dataclass(frozen=True)
class Candidate:
    phases_deg: dict[str, float]
    energy_error_mev_c: float
    injection_orbit_change_mm_mrad: np.ndarray
    mismatch_x: float
    mismatch_y: float
    score: float

    def as_dict(self) -> dict[str, object]:
        return {
            "klystron_phase_offsets_from_sad_deg": self.phases_deg,
            "energy_error_mev_c": self.energy_error_mev_c,
            "injection_orbit_change_from_nominal_mm_mrad": self.injection_orbit_change_mm_mrad.tolist(),
            "mismatch_x": self.mismatch_x,
            "mismatch_y": self.mismatch_y,
            "score": self.score,
        }


def _reference_energy_error(machine: ATF2LinacBTDRRFTrack) -> tuple[float, np.ndarray]:
    tracked = machine.track(machine.make_reference_bunch(), dr_turns=0)
    orbit = np.array((
        tracked.dr_injection.mean_x_mm, tracked.dr_injection.mean_xp_mrad,
        tracked.dr_injection.mean_y_mm, tracked.dr_injection.mean_yp_mrad,
    ))
    return (
        tracked.linac_bt_exit.mean_p_mev_c - float(machine.dr_synchronous_orbit[5]),
        orbit,
    )


def _build_machine(
    voltage_mv: float, *, handoff=None
) -> ATF2LinacBTDRRFTrack:
    return ATF2LinacBTDRRFTrack(
        input_momentum_mev_c=INPUT_MOMENTUM_MEV_C,
        cavity_voltage_mv=voltage_mv,
        # A zero reference phase fixes magnets at the nominal 80-MeV operating
        # point.  Later phase trials use set_klystron_phase_offsets_deg(),
        # which does not re-scale those magnets.
        magnet_reference_klystron_phase_offsets_deg={},
        handoff_mode="sad_optics_matched",
        handoff=handoff,
        dr_rf_mode="equilibrium",
    )


def _evaluate(
    machine: ATF2LinacBTDRRFTrack,
    phases_deg: Mapping[str, float],
    *, nominal_orbit: np.ndarray, nominal_mismatch: tuple[float, float],
) -> Candidate:
    machine.set_klystron_phase_offsets_deg(phases_deg)
    energy_error, orbit = _reference_energy_error(machine)
    finite = machine.track(
        machine.make_entrance_bunch(STUDY_TWISS, particles=16), dr_turns=0
    )
    mismatch_x = finite.dr_injection_optics.x.mismatch_to_design
    mismatch_y = finite.dr_injection_optics.y.mismatch_to_design
    orbit_change = orbit - nominal_orbit
    # Scales are deliberately stated instead of concealing a multi-objective
    # choice: 50 keV energy, 0.05 mm/mrad orbit movement, and 0.05 Bmag.
    score = float(np.sqrt(
        (energy_error / 0.05)**2
        + np.sum((orbit_change / 0.05)**2)
        + ((mismatch_x - nominal_mismatch[0]) / 0.05)**2
        + ((mismatch_y - nominal_mismatch[1]) / 0.05)**2
    ))
    return Candidate(dict(phases_deg), energy_error, orbit_change, mismatch_x, mismatch_y, score)


def _phase_payload(text: str | None) -> dict[str, float]:
    if text is None:
        return {}
    payload = json.loads(text)
    if not isinstance(payload, dict) or set(payload) - set(KLYSTRON_CAVITY_GROUPS):
        raise ValueError("initial-phase-offsets-json must be an object containing only L1 through L8")
    result = {name: float(value) for name, value in payload.items()}
    if not all(np.isfinite(value) for value in result.values()):
        raise ValueError("initial phases must be finite")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cavity-voltage-mv", type=float,
        default=MODEL_MATCHED_80MEV_CAVITY_VOLTAGE_MV,
        help="Common structure voltage; default is the direct-RFTrack 80-MeV model match.",
    )
    parser.add_argument("--phase-limit-deg", type=float, default=8.0)
    parser.add_argument("--coordinate-step-deg", type=float, default=1.0)
    parser.add_argument("--coordinate-passes", type=int, default=2)
    parser.add_argument("--initial-phase-offsets-json", help="JSON object of phase offsets relative to SAD, e.g. '{\"L1\":2.0}'.")
    parser.add_argument("--synthetic-starting-phase-error", action="store_true", help="Demonstrate recovery from a labelled synthetic eight-klystron phase error.")
    parser.add_argument("--turns", type=int, default=10)
    args = parser.parse_args()
    if args.phase_limit_deg <= 0.0 or args.coordinate_step_deg <= 0.0 or args.coordinate_passes < 1 or args.turns < 1:
        raise ValueError("phase limit, coordinate step, passes, and turns must be positive")

    voltage_mv = float(args.cavity_voltage_mv)
    baseline = _build_machine(voltage_mv)
    nominal_energy_error, nominal_orbit = _reference_energy_error(baseline)
    nominal_finite = baseline.track(
        baseline.make_entrance_bunch(STUDY_TWISS, particles=16), dr_turns=0
    )
    nominal_mismatch = (
        nominal_finite.dr_injection_optics.x.mismatch_to_design,
        nominal_finite.dr_injection_optics.y.mismatch_to_design,
    )
    # Keep this already-built machine's nominal physical handoff and magnet
    # fields fixed for all phase trials.  Rebuilding RFTrack lattices for each
    # sample is both slow and has shown an intermittent parser-level failure.
    machine = baseline
    phases = {name: 0.0 for name in KLYSTRON_CAVITY_GROUPS}
    phases.update(_phase_payload(args.initial_phase_offsets_json))
    if args.synthetic_starting_phase_error:
        phases.update(SYNTHETIC_STARTING_PHASE_ERROR_DEG)
    if any(abs(value) > args.phase_limit_deg for value in phases.values()):
        raise ValueError("initial phase is outside --phase-limit-deg")
    initial = _evaluate(
        machine, phases, nominal_orbit=nominal_orbit, nominal_mismatch=nominal_mismatch
    )
    best = initial
    scan_history: list[dict[str, object]] = []
    offsets = (-args.coordinate_step_deg, 0.0, args.coordinate_step_deg)
    for scan_pass in range(args.coordinate_passes):
        for name in KLYSTRON_CAVITY_GROUPS:
            candidates: list[Candidate] = []
            for offset in offsets:
                trial_phases = dict(best.phases_deg)
                trial_phases[name] = float(np.clip(
                    trial_phases[name] + offset, -args.phase_limit_deg, args.phase_limit_deg
                ))
                candidates.append(_evaluate(
                    machine, trial_phases, nominal_orbit=nominal_orbit,
                    nominal_mismatch=nominal_mismatch,
                ))
            selected = min(candidates, key=lambda candidate: candidate.score)
            best = selected
            scan_history.append({
                "pass": scan_pass + 1,
                "klystron": name,
                "selected": selected.as_dict(),
            })
    machine.set_klystron_phase_offsets_deg(best.phases_deg)
    final = machine.track(
        machine.make_entrance_bunch(STUDY_TWISS, particles=16),
        dr_turns=args.turns, record_turn_history=True,
    )
    print(json.dumps({
        "simulation_only": True,
        "input_momentum_mev_c": INPUT_MOMENTUM_MEV_C,
        "klystron_cavity_groups_assumption": {
            name: list(cavities) for name, cavities in KLYSTRON_CAVITY_GROUPS.items()
        },
        "phase_convention": "offset in degrees added to the SAD -1.67-degree-from-crest structure phase",
        "common_structure_voltage_mv": voltage_mv,
        "voltage_match": {
            "method": "direct RFTrack reference-track calibration performed once",
            "reference_exit_energy_error_mev_c": nominal_energy_error,
        },
        "magnet_and_handoff_rule": "Fixed at the nominal zero-phase 80-MeV solution during phase trials; only RF phases changed.",
        "nominal_zero_phase": {
            "energy_error_mev_c": nominal_energy_error,
            "mismatch_x": nominal_mismatch[0], "mismatch_y": nominal_mismatch[1],
        },
        "initial_candidate": initial.as_dict(),
        "selected_candidate": best.as_dict(),
        "coordinate_scan_history": scan_history,
        "short_turn_validation": {
            "turns_requested": args.turns,
            "turns_completed": final.completed_dr_turns,
            "ring_survival": final.dr_after_turns.survival_fraction_from_input,
        },
        "limitations": [
            "L1..L8 to CA pairs is a provisional model assumption, not a verified RF distribution table.",
            "Common structure voltage is an offline model match, not a klystron amplitude calibration or operating setpoint.",
            "No measured RF phase/amplitude, beam loading, jitter, timing, surveyed injection map, or complete aperture model is loaded.",
            "The finite bunch is a 1-mm-mrad synthetic study distribution, not a measured ATF bunch.",
        ],
    }, indent=2))


if __name__ == "__main__":
    main()
