"""Connect the Kubo-style correction sequence to RF-Track equilibrium emittance.

This connects the Kubo-style correction sequence to an RF-Track radiation
envelope.  The selectable transverse-envelope mode is explicitly diagnostic:
it remains meaningful while a historical RF thin-pillbox model has a
longitudinally unstable 6D map, but cannot by itself reproduce Table II.
"""

from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_correction import ATFDRRingCorrection
from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_lattice import build_atf_dr_lattice
from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_emittance import (
    build_equilibrium_lattice,
    equilibrium_emittance_for_lattices,
    find_synchronous_orbit,
    transverse_envelope_for_lattices,
)


QUANTUM_PARTICLES = int(os.environ.get("ATF_DR_KUBO_QUANTUM_PARTICLES", "32"))
ENVELOPE_MODEL = os.environ.get("ATF_DR_KUBO_EMITTANCE_MODEL", "transverse")
if ENVELOPE_MODEL not in {"transverse", "full_6d"}:
    raise ValueError("ATF_DR_KUBO_EMITTANCE_MODEL must be transverse or full_6d")


def _load_procedure_module():
    path = Path(__file__).with_name("kubo_style_procedure.py")
    spec = importlib.util.spec_from_file_location("atf_dr_kubo_procedure", path)
    if spec is None or spec.loader is None:  # pragma: no cover
        raise RuntimeError(f"Cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _apply_snapshot(lattice, response, snapshot):
    correction = ATFDRRingCorrection(lattice)
    for name, kick in zip(response["x_names"], snapshot["x_kick_rad"]):
        element = correction._single_element(name)
        current = correction._get_corrector_kick(element)
        current[0] = kick
        correction._set_corrector_kick(element, current)
    for name, kick in zip(response["y_names"], snapshot["y_kick_rad"]):
        element = correction._single_element(name)
        current = correction._get_corrector_kick(element)
        current[1] = kick
        correction._set_corrector_kick(element, current)
    for name, strength in zip(response["skew_names"], snapshot["skew_k1l"]):
        correction.set_skew_strength(name, strength)


def _equilibrium_for_stage(
    procedure, response, snapshot, nominal_synchronous, closed_orbit
):
    path = procedure.LATTICE_DATA_PATH or None
    deterministic = build_equilibrium_lattice(quantum=False, lattice_data_path=path)
    quantum = build_equilibrium_lattice(quantum=True, lattice_data_path=path)
    for lattice in (deterministic, quantum):
        procedure._apply_magnet_errors(lattice, procedure.MAGNET_ERROR_SCALE)
        _apply_snapshot(lattice, response, snapshot)
    solver = (
        transverse_envelope_for_lattices
        if ENVELOPE_MODEL == "transverse" else equilibrium_emittance_for_lattices
    )
    if ENVELOPE_MODEL == "transverse":
        # The correction operates at the nominal 1300-MeV/c transverse
        # reference, whereas the radiation/RF fixed point has a slightly
        # shifted synchronous momentum.  Re-close the same corrected lattice
        # at that momentum before taking a transverse one-turn derivative.
        transverse = build_atf_dr_lattice(
            lattice_data_path=path,
        )
        procedure._apply_magnet_errors(transverse, procedure.MAGNET_ERROR_SCALE)
        _apply_snapshot(transverse, response, snapshot)
        transverse_machine = ATFDRRingCorrection(transverse)
        closed = transverse_machine.find_closed_orbit(
            momentum_mev_c=float(nominal_synchronous[5]),
            initial_coordinates=closed_orbit,
            max_iterations=30,
        )
        reference = np.asarray(nominal_synchronous, dtype=float).copy()
        reference[:4] = closed.initial_coordinates
        return solver(
            deterministic, quantum, quantum_particles=QUANTUM_PARTICLES,
            reference_coordinates=reference,
        )
    return solver(
        deterministic, quantum, quantum_particles=QUANTUM_PARTICLES,
        initial_coordinates=nominal_synchronous,
    )


def main():
    procedure = _load_procedure_module()
    analysis = Path(__file__).resolve().parents[4] / "analysis" / "DR-RFTrack"
    probe_label = f"{procedure.STEERER_RESPONSE_KICK_RAD * 1e3:g}".replace(".", "p")
    response = procedure._load_or_build_responses(analysis / (
        f"kubo2003_nominal_kick_response_cache_{procedure.LATTICE_LABEL}_"
        f"{procedure.SKEW_FAMILY}_probe_{probe_label}mrad.npz"
    ))
    correction = procedure.run_case(response)
    nominal_synchronous = find_synchronous_orbit(
        build_equilibrium_lattice(
            quantum=False, lattice_data_path=procedure.LATTICE_DATA_PATH or None
        ), max_iterations=50
    ).coordinates

    emittance = {}
    stages = (
        "before_cod", "after_cod", "after_cod_dispersion", "after_coupling"
    )
    for stage in stages:
        try:
            result = _equilibrium_for_stage(
                procedure,
                response,
                correction["actuator_settings"][stage],
                nominal_synchronous,
                correction["stages"][stage]["closed_orbit_initial_coordinates_mm_mrad"],
            )
        except Exception as error:  # Preserve the failed stage as evidence.
            emittance[stage] = {"status": "failed", "error": str(error)}
            print(f"{stage}: envelope failed: {error}", flush=True)
            continue
        emittance[stage] = {
            "status": "ok",
            "vertical_like_eigen_pm_rad": result.eigen_emittances_m_rad[0] * 1e12,
            "projected_y_pm_rad": result.projected_y_m_rad * 1e12,
            "horizontal_like_eigen_pm_rad": result.eigen_emittances_m_rad[1] * 1e12,
            "damping_turns_one_over_e": result.damping_turns_one_over_e,
            "spectral_radius": result.spectral_radius,
        }
        print(
            f"{stage}: vertical-like eigen emittance "
            f"{emittance[stage]['vertical_like_eigen_pm_rad']:.3g} pm rad",
            flush=True,
        )

    payload = {
        "scope": (
            "Kubo-2003 sequential correction connected to an RF-Track "
            f"{ENVELOPE_MODEL} radiation envelope for one seed. "
            + (
                "The transverse mode excludes CT/delta covariance and therefore "
                "is not a Table-II emittance reproduction."
                if ENVELOPE_MODEL == "transverse" else
                "The full_6d mode includes CT/delta covariance, but uses the "
                "selected current daihon rather than the historical 2003 optics."
            )
        ),
        "paper_reference": "K. Kubo, Phys. Rev. ST Accel. Beams 6, 092801 (2003), Table II",
        "paper_table_ii_vertical_emittance_pm_rad": {
            "after_cod": 22.8,
            "after_cod_dispersion": 16.7,
            "after_coupling": 5.8,
        },
        "correction_observables": correction["stages"],
        "equilibrium_emittance": emittance,
        "model_settings": {
            "seed": procedure.SEED,
            "magnet_error_scale_of_table_i": procedure.MAGNET_ERROR_SCALE,
            "quantum_particles": QUANTUM_PARTICLES,
            "nominal_response_probe_mrad": procedure.STEERER_RESPONSE_KICK_RAD * 1e3,
            "paper_nominal_response_probe_mrad": 0.1,
            "cod_svd_rcond": procedure.RCOND_COD,
            "coupling_svd_rcond": procedure.RCOND_COUPLING,
            "first_correction_gain": procedure.KUBO_GAIN_FIRST,
            "dispersion_weight_r": procedure.DISPERSION_WEIGHT,
            "cod_dispersion_solver": procedure.COD_DISPERSION_SOLVER,
            "envelope_model": ENVELOPE_MODEL,
            "lattice_data_path": procedure.LATTICE_DATA_PATH or "default",
            "skew_corrector_family": procedure.SKEW_FAMILY,
            "ramp_initial_increment": procedure.PRELIMINARY_INITIAL_INCREMENT,
            "ramp_max_increment": procedure.PRELIMINARY_MAX_INCREMENT,
            "ramp_max_feedback_steps": procedure.PRELIMINARY_MAX_FEEDBACK_STEPS,
            "bpm_error_scale": procedure.BPM_ERROR_SCALE,
            "bpm_offset_um": procedure.SIGMA_BPM_OFFSET_MM
            * procedure.BPM_ERROR_SCALE * 1e3,
            "bpm_roll_mrad": procedure.SIGMA_BPM_ROLL_RAD
            * procedure.BPM_ERROR_SCALE * 1e3,
        },
    }
    # These SVD cutoffs materially change the correction command.  Include
    # them in the artifact name so a later diagnostic cannot silently replace
    # a different response-inversion study.
    cod_label = f"{procedure.RCOND_COD:g}".replace(".", "p").replace("-", "m")
    coupling_label = (
        f"{procedure.RCOND_COUPLING:g}".replace(".", "p").replace("-", "m")
    )
    suffix = (
        f"v2_{procedure.LATTICE_LABEL}_{procedure.SKEW_FAMILY}_"
        f"{ENVELOPE_MODEL}_seed_{procedure.SEED}_scale_{procedure.MAGNET_ERROR_SCALE:g}_"
        f"solver_{procedure.COD_DISPERSION_SOLVER}_codrcond_{cod_label}_"
        f"couplingrcond_{coupling_label}_qparticles_{QUANTUM_PARTICLES}"
    ).replace(".", "p")
    if not np.isclose(procedure.KUBO_GAIN_FIRST, 0.7):
        suffix += f"_firstgain_{procedure.KUBO_GAIN_FIRST:g}".replace(".", "p")
    if not np.isclose(procedure.DISPERSION_WEIGHT, 0.05):
        suffix += f"_dyweight_{procedure.DISPERSION_WEIGHT:g}".replace(".", "p")
    if not (
        np.isclose(procedure.PRELIMINARY_INITIAL_INCREMENT, 0.05)
        and np.isclose(procedure.PRELIMINARY_MAX_INCREMENT, 0.10)
        and procedure.PRELIMINARY_MAX_FEEDBACK_STEPS == 20
    ):
        suffix += (
            f"_ramp_{procedure.PRELIMINARY_INITIAL_INCREMENT:g}_"
            f"{procedure.PRELIMINARY_MAX_INCREMENT:g}_"
            f"steps_{procedure.PRELIMINARY_MAX_FEEDBACK_STEPS}"
        ).replace(".", "p")
    result_path = analysis / f"kubo2003_rftrack_emittance_envelope_{suffix}.json"
    result_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    labels = ("before COD", "COD", "COD + Dy", "coupling")
    vertical = [emittance[stage].get("vertical_like_eigen_pm_rad", np.nan) for stage in stages]
    projected = [emittance[stage].get("projected_y_pm_rad", np.nan) for stage in stages]
    dy = [correction["stages"][stage]["true_dy_all_bpm_mm"] for stage in stages]
    cxy = [correction["stages"][stage]["measured_cxy"] for stage in stages]
    figure, axes = plt.subplots(1, 3, figsize=(11.5, 3.4), constrained_layout=True)
    axes[0].plot(labels, vertical, "o-", label="vertical-like eigen")
    axes[0].plot(labels, projected, "s--", label="projected y")
    paper_labels = labels[1:]
    paper_values = [
        payload["paper_table_ii_vertical_emittance_pm_rad"][stage]
        for stage in stages[1:]
    ]
    axes[0].plot(
        paper_labels, paper_values, "k^:", label="Kubo Table II (reference)"
    )
    axes[0].set(title="radiation equilibrium", ylabel="emittance [pm rad]")
    axes[0].grid(alpha=0.3)
    axes[0].legend(fontsize=8)
    axes[1].plot(labels, dy, "o-", color="tab:blue")
    axes[1].set(title="true vertical dispersion", ylabel="RMS Dy [mm / delta]")
    axes[1].grid(alpha=0.3)
    axes[2].plot(labels, cxy, "o-", color="tab:orange")
    axes[2].set(title="measured coupling", ylabel="Cxy")
    axes[2].grid(alpha=0.3)
    figure.suptitle(
        f"RF-Track Kubo sequence + {ENVELOPE_MODEL.replace('_', ' ')} envelope "
        f"({procedure.MAGNET_ERROR_SCALE:g} Table-I magnet errors)"
    )
    figure_path = analysis / f"kubo2003_rftrack_emittance_envelope_{suffix}.png"
    figure.savefig(figure_path, dpi=180)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
