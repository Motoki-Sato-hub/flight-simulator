"""Plot model ring survival while scanning 80-MeV Linac RF phases.

The plot is deliberately labelled ``model ring survival`` rather than ATF
transmission.  It uses a fixed nominal Linac magnet configuration and fixed
design IPZT-to-KII.1 handoff, then varies only RF phase offsets.  Thus it is a
useful sensitivity/capture screen for the full Linac -> BT -> DR path, but
cannot predict physical final transmission before aperture, timing, RF, and
injection-map calibration.

The default common-phase scan moves all eight provisional klystron groups
together.  The individual scan moves one group at a time while the other
seven remain at their nominal SAD phase.  Every point starts from a fresh
synthetic finite bunch, so no radiation damping or particle loss leaks from a
previous point into the next.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import numpy as np

from Interfaces.ATF2.ATF2_LinacBTDR_RFTrack import KLYSTRON_CAVITY_GROUPS
from Interfaces.ATF2.scan_linac_klystron_phase_injection import (
    MODEL_MATCHED_80MEV_CAVITY_VOLTAGE_MV,
    STUDY_TWISS,
    _build_machine,
)

# Import the RFTrack/ATF interface first.  Some legacy interface modules
# select their Qt logging backend during import.  Override that selection only
# for this noninteractive plotting executable after those imports complete.
import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt


def _parse_values(text: str) -> tuple[float, ...]:
    values = tuple(float(value) for value in text.split(",") if value.strip())
    if not values or not all(np.isfinite(value) for value in values):
        raise ValueError("phase value list must contain finite comma-separated numbers")
    return values


def _phase_point(machine, phases: dict[str, float], *, particles: int, turns: int) -> dict[str, object]:
    machine.set_klystron_phase_offsets_deg(phases)
    tracked = machine.track(
        machine.make_entrance_bunch(STUDY_TWISS, particles=particles),
        dr_turns=turns,
        record_turn_history=True,
    )
    return {
        "phase_offsets_from_sad_deg": dict(phases),
        "linac_bt_transmission": tracked.linac_bt_exit.survival_fraction_from_input,
        "model_ring_survival": tracked.dr_after_turns.survival_fraction_from_input,
        "completed_dr_turns": tracked.completed_dr_turns,
        "energy_error_mev_c": (
            tracked.linac_bt_exit.mean_p_mev_c - float(machine.dr_synchronous_orbit[5])
        ),
        "mismatch_x": tracked.dr_injection_optics.x.mismatch_to_design,
        "mismatch_y": tracked.dr_injection_optics.y.mismatch_to_design,
    }


def _fresh_phase_scan_machine(voltage_mv: float):
    """Build one nominal machine; its fields/handoff remain fixed in the scan."""
    return _build_machine(voltage_mv)


def common_phase_scan(
    phase_values: Iterable[float], *, voltage_mv: float, particles: int, turns: int
) -> list[dict[str, object]]:
    machine = _fresh_phase_scan_machine(voltage_mv)
    records: list[dict[str, object]] = []
    for phase in phase_values:
        record = _phase_point(
            machine,
            {name: float(phase) for name in KLYSTRON_CAVITY_GROUPS},
            particles=particles, turns=turns,
        )
        record["phase_deg"] = float(phase)
        records.append(record)
    return records


def individual_phase_scan(
    phase_values: Iterable[float], *, voltage_mv: float, particles: int, turns: int
) -> dict[str, list[dict[str, object]]]:
    machine = _fresh_phase_scan_machine(voltage_mv)
    records: dict[str, list[dict[str, object]]] = {}
    for klystron in KLYSTRON_CAVITY_GROUPS:
        group_records: list[dict[str, object]] = []
        for phase in phase_values:
            phases = {name: 0.0 for name in KLYSTRON_CAVITY_GROUPS}
            phases[klystron] = float(phase)
            record = _phase_point(machine, phases, particles=particles, turns=turns)
            record["phase_deg"] = float(phase)
            group_records.append(record)
        records[klystron] = group_records
    return records


def _xy(records: list[dict[str, object]], field: str) -> tuple[np.ndarray, np.ndarray]:
    return (
        np.asarray([record["phase_deg"] for record in records], dtype=float),
        np.asarray([record[field] for record in records], dtype=float),
    )


def make_plot(
    common: list[dict[str, object]], individual: dict[str, list[dict[str, object]]],
    *, turns: int, output_png: Path,
) -> None:
    """Render compact common and per-klystron sensitivity plots."""
    figure, axes = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True)
    phase, energy = _xy(common, "energy_error_mev_c")
    _, survival = _xy(common, "model_ring_survival")
    _, mismatch_x = _xy(common, "mismatch_x")
    _, mismatch_y = _xy(common, "mismatch_y")
    axes[0, 0].plot(
        phase, energy, "o-",
        label="Finite-bunch mean Linac+BT exit − DR synchronous p",
    )
    axes[0, 0].axhline(0.0, color="black", lw=0.8)
    axes[0, 0].set(xlabel="Common L1–L8 phase offset from SAD [deg]", ylabel="Energy error [MeV/c]")
    axes[0, 0].legend(fontsize=8)
    axes[1, 0].plot(phase, survival, "o-", label=f"{turns}-turn model ring survival")
    axes[1, 0].set(xlabel="Common L1–L8 phase offset from SAD [deg]", ylabel="Survival from Linac entrance", ylim=(-0.05, 1.05))
    axes[1, 0].legend(fontsize=8)
    right_energy, right_survival = axes[0, 1], axes[1, 1]
    for name, records in individual.items():
        p, e = _xy(records, "energy_error_mev_c")
        _, s = _xy(records, "model_ring_survival")
        right_energy.plot(p, e, "o-", label=name)
        right_survival.plot(p, s, "o-", label=name)
    right_energy.axhline(0.0, color="black", lw=0.8)
    right_energy.set(xlabel="One-klystron phase offset from SAD [deg]", ylabel="Energy error [MeV/c]")
    right_energy.legend(ncol=2, fontsize=8)
    right_survival.set(xlabel="One-klystron phase offset from SAD [deg]", ylabel=f"{turns}-turn model ring survival", ylim=(-0.05, 1.05))
    right_survival.legend(ncol=2, fontsize=8)
    figure.suptitle(
        "80-MeV Linac → BT → DR RF phase sensitivity\n"
        "Fixed magnets/handoff; synthetic finite bunch; model survival, not physical transmission"
    )
    output_png.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_png, dpi=160)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cavity-voltage-mv", type=float, default=MODEL_MATCHED_80MEV_CAVITY_VOLTAGE_MV)
    parser.add_argument("--particles", type=int, default=4, help="Synthetic finite-bunch macroparticles per scan point.")
    parser.add_argument("--turns", type=int, default=10, help="DR turns per scan point.")
    parser.add_argument("--common-phase-values-deg", default="-12,-8,-4,0,4,8,12")
    parser.add_argument(
        "--individual-phase-values-deg", default="-16,-8,0,8,16",
        help="One-klystron offsets; the default spans the current model's loss onset.",
    )
    parser.add_argument("--output-png", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    args = parser.parse_args()
    if args.particles < 2 or args.turns < 1:
        raise ValueError("particles must be at least two and turns must be positive")
    common_values = _parse_values(args.common_phase_values_deg)
    individual_values = _parse_values(args.individual_phase_values_deg)
    common = common_phase_scan(
        common_values, voltage_mv=args.cavity_voltage_mv,
        particles=args.particles, turns=args.turns,
    )
    individual = individual_phase_scan(
        individual_values, voltage_mv=args.cavity_voltage_mv,
        particles=args.particles, turns=args.turns,
    )
    make_plot(common, individual, turns=args.turns, output_png=args.output_png)
    payload = {
        "simulation_only": True,
        "input_momentum_mev_c": 80.0,
        "common_structure_voltage_mv": args.cavity_voltage_mv,
        "phase_convention": "offset from the SAD -1.67-degree-from-crest structure phase",
        "klystron_cavity_groups_assumption": {
            name: list(cavities) for name, cavities in KLYSTRON_CAVITY_GROUPS.items()
        },
        "finite_bunch": {"particles": args.particles, "turns_per_point": args.turns},
        "common_phase_scan": common,
        "individual_phase_scan": individual,
        "individual_curve_note": (
            "The eight ideal-model traces may overlap exactly: each provisional "
            "L1..L8 group contains two identically calibrated SAD structures. "
            "Measured per-klystron amplitude/phase calibration will break this "
            "symmetry when supplied."
        ),
        "plot": str(args.output_png),
        "interpretation": (
            "model_ring_survival includes the entire tracked Linac->BT->DR path "
            "but is not physical final transmission without a calibrated injection "
            "map, RF/timing, and complete aperture/loss model"
        ),
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({
        "status": "passed", "plot": str(args.output_png), "data": str(args.output_json),
        "common_points": len(common),
        "individual_points": sum(len(records) for records in individual.values()),
    }, indent=2))


if __name__ == "__main__":
    main()
