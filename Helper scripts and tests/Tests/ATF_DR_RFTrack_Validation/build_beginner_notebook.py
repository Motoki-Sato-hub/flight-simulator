"""Generate the beginner-oriented ATF DR RFTrack correction notebook."""

from __future__ import annotations

import json
from pathlib import Path
from textwrap import dedent


REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
NOTEBOOK = REPOSITORY_ROOT / "Jupyter Notebooks/ATF2/DR_RFTrack_Correction_Validation.ipynb"


def markdown(source: str) -> dict:
    return {
        "cell_type": "markdown", "metadata": {},
        "source": dedent(source).lstrip().splitlines(True),
    }


def code(source: str) -> dict:
    return {
        "cell_type": "code", "execution_count": None, "metadata": {},
        "outputs": [], "source": dedent(source).lstrip().splitlines(True),
    }


def main():
    cells = [
        markdown(r"""
        # ATF DR RFTrack correction validation

        Input: `ATF_DR_20111111b_RFTrack_twiss.tfs`, an all-element RFTrack lattice
        generated from the 20111111b SAD daihon. Unlike the legacy target-optics
        `ATF_DR_twiss_file.tws`, it is read directly as an RFTrack lattice.

        This notebook starts with a one-turn periodic orbit, builds response matrices,
        applies COD / dispersion / skew corrections, and then evaluates the radiation
        equilibrium represented by repeated turns.
        """),
        code("""
        from pathlib import Path
        import sys

        import matplotlib.pyplot as plt
        import numpy as np
        import RF_Track as rft

        # Locate flight-simulator from the notebook working directory.
        REPOSITORY_ROOT = Path.cwd().resolve()
        while not (REPOSITORY_ROOT / "Interfaces").is_dir():
            if REPOSITORY_ROOT.parent == REPOSITORY_ROOT:
                raise RuntimeError("Run this notebook inside flight-simulator.")
            REPOSITORY_ROOT = REPOSITORY_ROOT.parent
        if str(REPOSITORY_ROOT) not in sys.path:
            sys.path.insert(0, str(REPOSITORY_ROOT))

        # Keep the SAD provenance and the directly read correction lattice explicit.
        SAD_DAIHON = REPOSITORY_ROOT.parent / "sad/operation/daihon/atfdr-design-20111111b.sad"
        TWISS_TFS = REPOSITORY_ROOT / "Interfaces/ATF2/DR_ATF2/ATF_DR_20111111b_RFTrack_twiss.tfs"
        LATTICE_JSON = REPOSITORY_ROOT / "Interfaces/ATF2/DR_ATF2/ATF_DR_RFTrack_lattice.json"

        from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_correction import ATFDRRingCorrection
        from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_emittance import (
            build_equilibrium_lattice,
            equilibrium_emittance_for_lattices,
            gaussian_bunch,
            track_emittance,
        )

        print("SAD daihon:", SAD_DAIHON)
        print("RFTrack Twiss/TFS lattice:", TWISS_TFS)
        print("RFTrack JSON used only for the radiation model:", LATTICE_JSON)
        assert SAD_DAIHON.exists() and TWISS_TFS.exists() and LATTICE_JSON.exists()


        def load_correction_ring():
            # Load the nominal one-turn correction lattice directly from TFS.
            return rft.Lattice(str(TWISS_TFS))
        """),
        markdown(r"""
        ## 1. One turn: build the ring and find its periodic orbit

        `find_closed_orbit()` repeatedly tracks **one turn** and solves
        \(f(z)-z=0\). It does not perform long-term multi-turn tracking.
        """),
        code("""
        one_turn_lattice = load_correction_ring()
        one_turn = ATFDRRingCorrection(one_turn_lattice)
        orbit = one_turn.find_closed_orbit()
        dispersion = one_turn.measure_dispersion(relative_momentum_step=3.5e-3)

        print(f"circumference = {one_turn_lattice.get_length():.6f} m")
        print(f"BPMs = {len(orbit.bpm_names)}")
        print("closed-orbit initial coordinates [mm, mrad] =", orbit.initial_coordinates)
        print("RMS Dx, Dy [mm/delta] =", np.sqrt(np.mean(dispersion.x**2)), np.sqrt(np.mean(dispersion.y**2)))
        """),
        markdown(r"""
        ## 2. COD correction from a one-turn response matrix

        The nominal lattice is the correction model. A separate lattice is the virtual
        machine and contains an unregistered kick. The response is \(R_{ij}=\partial x_i/\partial\theta_j\).
        """),
        code("""
        def rms(values):
            values = np.asarray(values, dtype=float)
            return float(np.sqrt(np.mean(values**2)))


        def add_corrector_kick(machine, name, plane, kick_rad):
            element = machine._single_element(name)
            kick = np.asarray(element.get_strength(), dtype=float)
            kick[plane] += kick_rad
            element.set_strength(*kick)


        def bounded_svd_response(matrix, residual, maximum_kick, count=6, rcond=1e-4):
            # Pick phase-distributed useful columns, solve only those columns,
            # then impose a physical corrector limit.
            selected = ATFDRRingCorrection._greedy_columns(matrix, -residual, count, rcond)
            selected_commands, _, rank = ATFDRRingCorrection._svd_solve(
                matrix[:, selected], -residual, rcond
            )
            commands = np.zeros(matrix.shape[1])
            commands[selected] = np.clip(selected_commands, -maximum_kick, maximum_kick)
            return commands, selected, rank


        model = ATFDRRingCorrection(load_correction_ring())
        virtual_machine = ATFDRRingCorrection(load_correction_ring())

        # Hidden fault: it exists only in the virtual machine.
        add_corrector_kick(virtual_machine, "ZH1R", plane=0, kick_rad=1.0e-5)
        before = virtual_machine.find_closed_orbit().x

        correctors = model.get_corrector_names("x")[1:13]
        response = model.compute_orbit_response("x", corrector_names=correctors)
        commands, selected, rank = bounded_svd_response(
            response.matrix, before, maximum_kick=5.0e-5
        )
        for name, command in zip(correctors, commands):
            add_corrector_kick(virtual_machine, name, plane=0, kick_rad=float(command))
        after = virtual_machine.find_closed_orbit().x

        print(f"SVD rank = {rank}")
        print(f"COD RMS: {rms(before):.4e} -> {rms(after):.4e} mm")
        plt.plot(before, label="before")
        plt.plot(after, label="after")
        plt.xlabel("BPM index"); plt.ylabel("x [mm]"); plt.grid(alpha=0.3); plt.legend();
        """),
        markdown(r"""
        ## 3. Dispersion correction

        A vertical-steerer fault changes both vertical COD and vertical dispersion.
        The stacked response gives dispersion a relative weight of 0.05, as in the SAD-style procedure.
        """),
        code("""
        model = ATFDRRingCorrection(load_correction_ring())
        virtual_machine = ATFDRRingCorrection(load_correction_ring())
        target = model.measure_dispersion(relative_momentum_step=3.5e-3)

        # Hidden fault in the virtual machine only.
        add_corrector_kick(virtual_machine, "ZV1R", plane=1, kick_rad=1.0e-5)

        orbit_before = virtual_machine.find_closed_orbit()
        dispersion_before = virtual_machine.measure_dispersion(relative_momentum_step=3.5e-3)
        correctors = model.get_corrector_names("y")[1:9]
        response = model.compute_dispersion_response(
            "y", corrector_names=correctors, relative_momentum_step=3.5e-3
        )
        matrix = np.vstack((0.05 * response.dispersion_matrix, response.orbit_matrix))
        residual = np.concatenate((0.05 * (dispersion_before.y - target.y), orbit_before.y))
        commands, selected, rank = bounded_svd_response(
            matrix, residual, maximum_kick=5.0e-5
        )
        for name, command in zip(correctors, commands):
            add_corrector_kick(virtual_machine, name, plane=1, kick_rad=float(command))
        dispersion_after = virtual_machine.measure_dispersion(relative_momentum_step=3.5e-3)

        print(f"SVD rank = {rank}")
        print(f"RMS Dy error: {rms(dispersion_before.y-target.y):.4e} -> {rms(dispersion_after.y-target.y):.4e} mm/delta")
        plt.plot(dispersion_before.y-target.y, label="before")
        plt.plot(dispersion_after.y-target.y, label="after")
        plt.xlabel("BPM index"); plt.ylabel("Dy - target [mm/delta]"); plt.grid(alpha=0.3); plt.legend();
        """),
        markdown(r"""
        ## 4. Skew correction

        Two horizontal correctors are probed. Their vertical BPM responses form the coupling observable.
        A skew-quadrupole response matrix proposes compensation without resetting the hidden skew fault.
        """),
        code("""
        model = ATFDRRingCorrection(load_correction_ring())
        virtual_machine = ATFDRRingCorrection(load_correction_ring())
        probes = ("ZH1R", "ZH2R")
        available_skews = model.get_skew_corrector_names()
        skew_names = tuple(
            available_skews[index]
            for index in np.linspace(0, len(available_skews) - 2, 8, dtype=int)
        )
        bpm_names = model.bpm_names
        response = model.compute_coupling_response(
            probe_corrector_names=probes, skew_corrector_names=skew_names,
            bpm_names=bpm_names, probe_perturbation=1e-5, skew_perturbation=1e-5,
        )

        hidden_skew = available_skews[-1]
        virtual_machine.set_skew_strength(hidden_skew, 1.0e-3)
        _, signal_before_2d = virtual_machine._measure_vertical_response_to_horizontal_probes(
            probes, bpm_names=bpm_names, probe_perturbation=1e-5
        )
        signal_before = signal_before_2d.reshape(-1)
        commands, selected, rank = bounded_svd_response(
            response.matrix, signal_before, maximum_kick=2.0e-3, rcond=0.4
        )
        for name, command in zip(skew_names, commands):
            virtual_machine.set_skew_strength(name, virtual_machine.get_skew_strength(name) + float(command))
        _, signal_after_2d = virtual_machine._measure_vertical_response_to_horizontal_probes(
            probes, bpm_names=bpm_names, probe_perturbation=1e-5
        )
        signal_after = signal_after_2d.reshape(-1)

        print(f"SVD rank = {rank}")
        print(f"coupling signal RMS: {rms(signal_before):.4e} -> {rms(signal_after):.4e} mm/rad")
        """),
        markdown(r"""
        ## 5. From one turn to multi-turn equilibrium

        The correction lattice above is loaded directly from the checked-in Twiss/TFS.
        Radiation equilibrium additionally needs RF and radiation, so it is built from the same
        SAD-derived JSON with those explicit options enabled. It uses two numerical **one-turn** maps: deterministic damping \(M\) and
        quantum diffusion \(Q\). The stationary covariance is calculated from
        \(\Sigma=M\Sigma M^T+Q\). This is the long-time, multi-turn limit without tracking a particle
        one turn at a time for hundreds of thousands of turns.
        """),
        code("""
        deterministic_lattice = build_equilibrium_lattice(quantum=False)
        quantum_lattice = build_equilibrium_lattice(quantum=True)
        equilibrium = equilibrium_emittance_for_lattices(
            deterministic_lattice, quantum_lattice, quantum_particles=64
        )

        print(f"spectral radius = {equilibrium.spectral_radius:.9f}  (< 1 means damping is stable)")
        print(f"1/e damping time = {equilibrium.damping_turns_one_over_e:.0f} turns")
        print("vertical-like, horizontal-like eigen emittance [pm rad] =", np.asarray(equilibrium.eigen_emittances_m_rad) * 1e12)
        print("projected x, y emittance [pm rad] =", equilibrium.projected_x_m_rad * 1e12, equilibrium.projected_y_m_rad * 1e12)
        """),
        markdown(r"""
        ## 6. Explicit multi-turn tracking

        This cell actually calls `lattice.track()` once per turn for a quantum-radiation bunch.
        Twenty turns verify the implementation and show the sampled emittance evolution. They are
        deliberately far fewer than the printed damping time, so the stationary value above—not this
        short demonstration—is the equilibrium estimate. Increase `TRACKED_TURNS` only for a dedicated run.
        """),
        code("""
        TRACKED_TURNS = 20
        TRACKED_PARTICLES = 32

        # Start above the equilibrium emittance to make damping visible in a short run.
        initial_bunch = gaussian_bunch(
            equilibrium.synchronous_orbit,
            particles=TRACKED_PARTICLES,
            emit_x_m_rad=1.0e-9,
            emit_y_m_rad=1.0e-11,
            beta_x_m=5.0,
            beta_y_m=5.0,
            sigma_delta=1.0e-3,
            sigma_ct_mm=1.0,
            seed=20260916,
        )
        tracked_bunch, tracking_samples = track_emittance(
            quantum_lattice, initial_bunch,
            turns=TRACKED_TURNS, sample_every=5,
        )

        sample_turns = [sample.turns for sample in tracking_samples]
        vertical_like_pm = [sample.eigen_emittances_m_rad[0] * 1e12 for sample in tracking_samples]
        horizontal_like_pm = [sample.eigen_emittances_m_rad[1] * 1e12 for sample in tracking_samples]
        print(f"tracked turns = {TRACKED_TURNS}; survived particles = {tracked_bunch.size()}")
        print(f"1/e damping time = {equilibrium.damping_turns_one_over_e:.0f} turns")
        print(f"stationary vertical-like emittance = {equilibrium.eigen_emittances_m_rad[0] * 1e12:.3g} pm rad")

        plt.plot(sample_turns, vertical_like_pm, "o-", label="vertical-like")
        plt.plot(sample_turns, horizontal_like_pm, "s-", label="horizontal-like")
        plt.xlabel("turn"); plt.ylabel("eigen emittance [pm rad]")
        plt.grid(alpha=0.3); plt.legend();
        """),
        markdown("""
        ## Next step

        For a current-daihon study, use `kubo_style_procedure.py` for the full
        COD → dispersion → coupling sequence, `six_d_equilibrium_envelope.py`
        for the correction-stage multi-turn equilibrium emittance, and
        `current_daihon_ensemble.py` for seed statistics.
        """),
    ]
    notebook = {
        "cells": cells,
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python", "version": "3"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    for number, cell in enumerate(notebook["cells"]):
        cell["id"] = f"atf-dr-{number:02d}"
    NOTEBOOK.write_text(json.dumps(notebook, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    print(NOTEBOOK)


if __name__ == "__main__":
    main()
