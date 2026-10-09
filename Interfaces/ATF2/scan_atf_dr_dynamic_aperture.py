"""Offline finite-turn DR acceptance scans, not a fixed physical aperture.

The default includes deterministic RF/radiation, so this is a *damped capture*
map rather than a conservative Hamiltonian dynamic-aperture calculation.
Every probe and injected particle is tracked independently to retain its loss
turn without relying on RF-Track particle ordering after losses. No controls
are accessed and no literature radius is imposed on the tracking.
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
import time

import numpy as np
import RF_Track as rft

from Interfaces.ATF2.ATF2_LinacBTDR_RFTrack import (
    ATF2LinacBTDRRFTrack,
    EntranceBunchTwiss,
    load_entrance_bunch_json,
)
from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_lattice import (
    HISTORICAL_ATF_DR_APERTURE_SOURCE,
    get_historical_extraction_kicker_apertures,
)
from Interfaces.ATF2.linac_bt_dr_model_config import MODEL_ENERGY_MATCHED_CAVITY_VOLTAGE_MV


def track_probe(machine, coordinates, *, turns, guard_position_mm=100.0,
                guard_angle_mrad=100.0):
    """Track one KII.1 particle; guards are numerical screens, NOT walls.

    Guards are checked at the entry and each completed turn, not everywhere
    along the ring. A guard rejection is reported separately from RF-Track
    loss. They do not establish long-term stability or a full physical aperture.
    """
    if isinstance(turns, bool) or not isinstance(turns, (int, np.integer)) or turns < 1:
        raise ValueError("turns must be a positive integer")
    guards = np.asarray((guard_position_mm, guard_angle_mrad), dtype=float)
    if not np.all(np.isfinite(guards)) or np.any(guards <= 0):
        raise ValueError("numerical guards must be finite and positive")
    z = np.asarray(coordinates, dtype=float)
    if z.shape != (6,) or not np.all(np.isfinite(z)) or z[5] <= 0:
        raise ValueError("probe must be a finite six-vector with positive momentum")
    if machine.dr_first_turn_kicks_rad:
        raise ValueError("periodic acceptance scan requires zero pulsed kickers")
    maxima = np.zeros(4)
    for turn in range(turns + 1):
        if not np.all(np.isfinite(z)) or z[5] <= 0:
            status = "nonfinite_or_nonpositive_momentum"
            break
        displacement = np.abs(z[:4] - machine.dr_closed_orbit)
        maxima = np.maximum(maxima, displacement)
        if np.any(displacement[[0, 2]] > guard_position_mm) or np.any(
            displacement[[1, 3]] > guard_angle_mrad
        ):
            status = "numerical_guard"
            break
        if turn == turns:
            status = "survived"
            break
        bunch = rft.Bunch6d(rft.electronmass, 1.0, -1.0, z.copy())
        tracked = machine.dr_lattice.track(bunch)
        if tracked.size() != 1:
            status = "rftrack_loss"
            turn += 1
            break
        z = np.asarray(tracked.get_phase_space()[0], dtype=float)
    return {
        "initial_mm_mrad_mm_c_mev_c": np.asarray(coordinates).tolist(),
        "status": status,
        "survived": status == "survived",
        "first_rejected_turn": None if status == "survived" else int(turn),
        "completed_turns": int(turn),
        "max_entry_displacement_mm_mrad": maxima.tolist(),
    }


def make_scan_coordinates(machine, *, radius_mm, ray_deg, phase_x_deg,
                          phase_y_deg, delta_p, time_offset_mm_c=0.0):
    """Sample transverse Courant--Snyder phases around the dispersive orbit.

    Radius is a local mm amplitude parameter, NOT invariant acceptance.
    At phase zero x'= -alpha*x/beta; mm / m numerically gives mrad.
    Both planes' betatron phases are sampled independently.
    """
    z = np.array(machine.dr_synchronous_orbit if machine.dr_synchronous_orbit is not None
                 else [*machine.dr_closed_orbit, 0.0, machine.dr_momentum_mev_c], dtype=float)
    theta = np.deg2rad(ray_deg)
    amplitudes = (radius_mm * np.cos(theta), radius_mm * np.sin(theta))
    bx, ax, by, ay = machine.dr_start_twiss
    for start, amplitude, beta, alpha, phase in (
        (0, amplitudes[0], bx, ax, phase_x_deg),
        (2, amplitudes[1], by, ay, phase_y_deg),
    ):
        psi = np.deg2rad(phase)
        z[start] += amplitude * np.cos(psi)
        z[start + 1] += amplitude * (-alpha * np.cos(psi) - np.sin(psi)) / beta
    z[:4] += machine.dr_start_dispersion_mm_mrad * delta_p
    z[4] += time_offset_mm_c
    z[5] *= 1 + delta_p
    return z


def scan_acceptance(machine, *, radii_mm, rays_deg, phases_deg, delta_ps, turns,
                    time_offset_mm_c=0.0):
    radii = np.asarray(radii_mm, dtype=float)
    if radii.ndim != 1 or len(radii) == 0 or radii[0] != 0 or np.any(np.diff(radii) <= 0):
        raise ValueError("radii must be increasing, starting at zero")
    if not np.all(np.isfinite(radii)):
        raise ValueError("radii must be finite")
    for values in (rays_deg, phases_deg, delta_ps):
        if len(values) == 0 or not np.all(np.isfinite(values)):
            raise ValueError("ray, phase and momentum lists must be nonempty and finite")
    if not np.isfinite(time_offset_mm_c) or any(dp <= -1 for dp in delta_ps):
        raise ValueError("finite time offset and delta-p > -1 required")
    cases, boundaries = [], []
    for dp, ray, px, py in itertools.product(delta_ps, rays_deg, phases_deg, phases_deg):
        ray_cases = []
        for radius in radii:
            coordinates = make_scan_coordinates(
                machine, radius_mm=float(radius), ray_deg=ray, phase_x_deg=px,
                phase_y_deg=py, delta_p=dp, time_offset_mm_c=time_offset_mm_c,
            )
            outcome = track_probe(machine, coordinates, turns=turns)
            case = dict(radius_mm=float(radius), ray_deg=float(ray),
                        phase_x_deg=float(px), phase_y_deg=float(py), delta_p=float(dp),
                        **outcome)
            cases.append(case)
            ray_cases.append(case)
        first_failure = next((i for i, c in enumerate(ray_cases) if not c["survived"]), None)
        boundaries.append({
            "delta_p": float(dp), "ray_deg": float(ray),
            "phase_x_deg": float(px), "phase_y_deg": float(py),
            "last_contiguous_surviving_radius_mm": (
                None if first_failure == 0 else
                float(radii[-1] if first_failure is None else radii[first_failure - 1])
            ),
            "first_failing_radius_mm": None if first_failure is None else float(radii[first_failure]),
            "boundary_beyond_scan": first_failure is None,
            "surviving_samples_after_first_failure": (
                0 if first_failure is None else sum(c["survived"] for c in ray_cases[first_failure + 1:])
            ),
        })
    return {"cases": cases, "sampled_boundaries": boundaries}


def evaluate_injected_phase_space(machine, phase_space, *, turns):
    """Directly test each injected particle; do not interpolate a coarse DA map.

    Input is already in the DR entry frame; no further rephasing/rematching is
    performed. The reported fraction assumes equal-weight macroparticles.
    """
    phase_space = np.asarray(phase_space, dtype=float)
    if phase_space.ndim != 2 or phase_space.shape[1] != 6:
        raise ValueError("injected phase space must have shape (N, 6)")
    outcomes = [track_probe(machine, z, turns=turns) for z in phase_space]
    survivors = sum(o["survived"] for o in outcomes)
    return {"particles": len(outcomes), "survivors": survivors,
            "equal_weight_survival_fraction": survivors / len(outcomes) if outcomes else None,
            "particles_results": outcomes}


def _numbers(value):
    values = [float(item) for item in value.split(",")]
    if not values or not np.all(np.isfinite(values)):
        raise argparse.ArgumentTypeError("expected comma-separated finite numbers")
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--turns", type=int, default=10)
    parser.add_argument("--radii-mm", type=_numbers, default=[0., 1., 3., 6., 8.])
    parser.add_argument("--rays-deg", type=_numbers, default=[0., 90., 180., 270.])
    parser.add_argument("--phases-deg", type=_numbers, default=[0., 90.])
    parser.add_argument("--delta-p", type=_numbers, default=[0.])
    parser.add_argument("--time-offset-mm-c", type=float, default=0.)
    parser.add_argument("--aperture-screen", choices=("none", "historical-kix"), default="none")
    parser.add_argument("--rf-mode", choices=("equilibrium", "disabled"), default="equilibrium")
    parser.add_argument("--entrance-bunch-json", type=Path)
    parser.add_argument("--pipeline-particles", type=int, default=0,
                        help="Optionally test an equal-weight synthetic Linac entrance bunch.")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.turns < 1 or args.pipeline_particles < 0:
        parser.error("turns must be positive; pipeline-particles must be nonnegative")
    if args.entrance_bunch_json and args.pipeline_particles:
        parser.error("choose direct input or synthetic pipeline particles, not both")
    if not np.isfinite(args.time_offset_mm_c) or any(dp <= -1 for dp in args.delta_p):
        parser.error("time offset must be finite and delta-p must exceed -1")
    started = time.perf_counter()
    apertures = get_historical_extraction_kicker_apertures() if args.aperture_screen == "historical-kix" else {}
    machine = ATF2LinacBTDRRFTrack(
        cavity_voltage_mv=MODEL_ENERGY_MATCHED_CAVITY_VOLTAGE_MV,
        handoff_mode="sad_optics_matched", dr_rf_mode=args.rf_mode,
        dr_apertures_mm=apertures,
    )
    report = {
        "simulation_only": True, "entry": machine.dr_entry, "turns": args.turns,
        "kind": "damped_RF_capture" if args.rf_mode == "equilibrium" else "RF_off_transverse_survival",
        "rf_mode": args.rf_mode, "quantum_radiation": False,
        "synchronous_orbit": None if machine.dr_synchronous_orbit is None else machine.dr_synchronous_orbit.tolist(),
        "twiss_beta_m_alpha": list(machine.dr_start_twiss),
        "aperture_screen": args.aperture_screen, "physical_apertures": apertures,
        "numerical_guards": {"position_mm": 100., "angle_mrad": 100., "checked_at": "entry and turn boundaries"},
        "time_offset_mm_c": args.time_offset_mm_c,
        "literature_reference": {"source": HISTORICAL_ATF_DR_APERTURE_SOURCE,
                                 "historical_dynamic_aperture_mm": 6.,
                                 "used_as_tracking_cut": False},
        "limitations": [
            "Finite-turn sampled model acceptance, not measured ATF dynamic aperture or final transmission.",
            "Default RF/radiation capture is not conservative Hamiltonian dynamic aperture; RF-off mode is not longitudinal capture.",
            "Historical +/-6 mm has an unspecified observation point here and is not a KII.1 acceptance limit.",
            "Physical apertures are absent or a partial assumed-shape historical KIX screen, never a complete chamber model.",
            "No calibrated field-error ensemble, full nonlinear wiggler field, or injection pulse is included.",
            "The optional pipeline uses the provisional BT-to-DR Twiss/dispersion handoff and ideal mean timing alignment.",
            "Each probe is independent, with no charge-dependent collective loss mechanisms; injected macroparticles are assumed equally weighted.",
        ],
    }
    report.update(scan_acceptance(
        machine, radii_mm=args.radii_mm, rays_deg=args.rays_deg,
        phases_deg=args.phases_deg, delta_ps=args.delta_p, turns=args.turns,
        time_offset_mm_c=args.time_offset_mm_c,
    ))
    if args.entrance_bunch_json or args.pipeline_particles:
        bunch = (machine.make_entrance_bunch_from_phase_space(load_entrance_bunch_json(args.entrance_bunch_json))
                 if args.entrance_bunch_json else machine.make_entrance_bunch(
                     EntranceBunchTwiss(emitt_x_norm_mm_mrad=1., emitt_y_norm_mm_mrad=1.,
                                        beta_x_m=1.93, beta_y_m=1.93), particles=args.pipeline_particles))
        initial_charge = abs(float(bunch.get_total_charge()))
        bt_exit = machine.linac_bt_lattice.track(bunch)
        dr_input = machine.handoff_bunch(bt_exit)
        coordinates = np.asarray(dr_input.get_phase_space(), dtype=float).reshape(-1, 6).copy()
        coordinates[:, 4] += args.time_offset_mm_c
        evaluation = evaluate_injected_phase_space(machine, coordinates, turns=args.turns)
        transport = abs(float(dr_input.get_total_charge())) / initial_charge if initial_charge else None
        evaluation["linac_bt_charge_transmission"] = transport
        fraction = evaluation["equal_weight_survival_fraction"]
        evaluation["model_end_to_end_survival"] = (transport * fraction if transport is not None and fraction is not None else
                                                   0. if initial_charge and not len(coordinates) else None)
        evaluation["entrance_source"] = str(args.entrance_bunch_json) if args.entrance_bunch_json else "synthetic Twiss bunch"
        report["pipeline_bunch"] = evaluation
    report["elapsed_seconds"] = time.perf_counter() - started
    args.output_dir.mkdir(parents=True, exist_ok=True)
    path = args.output_dir / "dr_acceptance.json"
    path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    _plot(report, args.output_dir / "dr_acceptance.png")
    print(json.dumps({"report": str(path), "plot": str(args.output_dir / "dr_acceptance.png"),
                      "probes": len(report["cases"]), "survived": sum(c["survived"] for c in report["cases"]),
                      "elapsed_seconds": report["elapsed_seconds"],
                      "pipeline_bunch_survival": report.get("pipeline_bunch", {}).get("model_end_to_end_survival")}, indent=2))


def _plot(report, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    slices = sorted({(c["delta_p"], c["phase_x_deg"], c["phase_y_deg"]) for c in report["cases"]})
    fig, axes = plt.subplots(len(slices), 1, figsize=(6, 4 * len(slices)), squeeze=False)
    for ax, key in zip(axes.flat, slices):
        for survived, marker, color in ((True, "o", "tab:blue"), (False, "x", "tab:red")):
            selected = [c for c in report["cases"] if (c["delta_p"], c["phase_x_deg"], c["phase_y_deg"]) == key and c["survived"] == survived]
            ax.scatter([c["radius_mm"] * np.cos(np.deg2rad(c["ray_deg"])) for c in selected],
                       [c["radius_mm"] * np.sin(np.deg2rad(c["ray_deg"])) for c in selected],
                       marker=marker, color=color, label="survived" if survived else "lost / guard")
        ax.set(xlabel="signed x CS amplitude at entry [mm]", ylabel="signed y CS amplitude at entry [mm]",
               title=f"{report['kind']}: {report['turns']} turns; dp/p={key[0]:g}; phases={key[1]:g},{key[2]:g} deg")
        ax.set_aspect("equal", adjustable="datalim")
        ax.grid(True)
        ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    main()
