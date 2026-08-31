"""Smoke tests for the SAD-derived ATF DR RF-Track lattice.

Run from the flight-simulator root with an RF-Track-enabled Python:

    PYTHONPATH=. python "Helper scripts and tests/Tests/TestATFDRRFTrackLattice.py"
"""

import unittest

import numpy as np
import RF_Track as rft

from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_lattice import (
    NOMINAL_MOMENTUM_MEV_C,
    REFERENCE_SAD_DAIHON,
    build_atf_dr_lattice,
    get_lattice_metadata,
)
from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_correction import (
    ATFDRRingCorrection,
)


class TestATFDRRFTrackLattice(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.metadata = get_lattice_metadata()
        cls.lattice = build_atf_dr_lattice()

    def test_reference_and_geometry(self):
        self.assertEqual(
            self.metadata["reference_sad_daihon"], REFERENCE_SAD_DAIHON
        )
        self.assertAlmostEqual(
            self.lattice.get_length(), self.metadata["circumference_m"], places=10
        )

    def test_machine_devices(self):
        correctors = self.lattice.get_correctors()
        self.assertEqual(len(self.lattice.get_bpms()), 98)
        self.assertEqual(
            len([c for c in correctors if c.get_name().startswith("ZH")]), 50
        )
        self.assertEqual(
            len([c for c in correctors if c.get_name().startswith("ZV")]), 51
        )
        self.assertEqual(len(self.lattice.get_quadrupoles()), 100)
        self.assertEqual(len(self.lattice.get_sbends()), 36)
        self.assertEqual(len(self.lattice.get_rf_elements()), 0)
        correction = ATFDRRingCorrection(self.lattice)
        self.assertEqual(len(correction.get_skew_corrector_names()), 68)

    def test_nominal_particle_completes_one_turn_without_energy_change(self):
        phase_space = np.array(
            [0.0, 0.0, 0.0, 0.0, 0.0, NOMINAL_MOMENTUM_MEV_C]
        )
        bunch = rft.Bunch6d(rft.electronmass, 0.0, -1.0, phase_space)
        tracked = self.lattice.track(bunch)
        self.assertEqual(tracked.size(), 1)
        self.assertAlmostEqual(
            tracked.get_phase_space()[0, 5], NOMINAL_MOMENTUM_MEV_C, places=8
        )

    def test_unvalidated_rf_mode_is_rejected(self):
        with self.assertRaises(NotImplementedError):
            build_atf_dr_lattice(rf_mode="longitudinal")

    def test_periodic_closed_orbit_and_dispersion(self):
        correction = ATFDRRingCorrection(build_atf_dr_lattice())
        orbit = correction.find_closed_orbit()
        dispersion = correction.measure_dispersion(relative_momentum_step=1e-3)
        self.assertLess(orbit.residual_norm, 1e-8)
        self.assertEqual(orbit.bpm_positions.shape, (98, 2))
        self.assertEqual(dispersion.values.shape, (98, 2))
        self.assertTrue(np.all(np.isfinite(dispersion.values)))

    def test_sad_style_limited_steerer_correction_in_both_planes(self):
        cases = (
            ("x", "ZH1R", 0, 0.10),
            ("y", "ZV1R", 1, 0.20),
        )
        for plane, error_corrector, strength_index, maximum_ratio in cases:
            with self.subTest(plane=plane):
                lattice = build_atf_dr_lattice()
                correction = ATFDRRingCorrection(lattice)
                strength = np.zeros(2)
                strength[strength_index] = 1e-3
                lattice[error_corrector].set_strength(*strength)

                available = correction.get_corrector_names(plane)
                response = correction.compute_orbit_response(
                    plane,
                    corrector_names=available[1:13],
                    perturbation=1e-5,
                )
                suggestion = correction.propose_orbit_correction(
                    response,
                    max_correctors=6,
                    rcond=1e-4,
                )
                corrected = correction.apply_orbit_correction(suggestion)
                corrected_plane = corrected.x if plane == "x" else corrected.y
                rms_after = float(np.sqrt(np.mean(corrected_plane**2)))

                self.assertLessEqual(len(suggestion.selected_correctors), 6)
                self.assertLess(
                    suggestion.rms_predicted,
                    maximum_ratio * suggestion.rms_before,
                )
                self.assertAlmostEqual(
                    rms_after,
                    suggestion.rms_predicted,
                    delta=max(1e-8, suggestion.rms_predicted * 0.01),
                )

    def test_joint_dispersion_and_orbit_correction(self):
        lattice = build_atf_dr_lattice()
        correction = ATFDRRingCorrection(lattice)
        target = correction.measure_dispersion(relative_momentum_step=1e-3)

        quadrupole = lattice["QM10R.1"]
        p_over_q = -NOMINAL_MOMENTUM_MEV_C
        nominal_k1 = quadrupole.get_K1(p_over_q)
        quadrupole.set_K1(p_over_q, nominal_k1 * 1.01)

        response = correction.compute_dispersion_response(
            "x",
            corrector_names=correction.get_corrector_names("x")[:4],
            relative_momentum_step=1e-3,
            perturbation=1e-5,
        )
        suggestion = correction.propose_dispersion_correction(
            response,
            target_dispersion=target.x,
            dispersion_weight=0.05,
            orbit_weight=1.0,
            max_correctors=4,
            rcond=1e-4,
        )
        orbit, dispersion = correction.apply_dispersion_correction(suggestion)
        actual_dispersion_rms = float(
            np.sqrt(np.mean((dispersion.x - target.x) ** 2))
        )
        actual_orbit_rms = float(np.sqrt(np.mean(orbit.x**2)))

        self.assertLess(
            suggestion.dispersion_rms_predicted,
            0.01 * suggestion.dispersion_rms_before,
        )
        self.assertAlmostEqual(
            actual_dispersion_rms,
            suggestion.dispersion_rms_predicted,
            delta=max(1e-5, suggestion.dispersion_rms_predicted * 0.30),
        )
        self.assertLess(
            actual_dispersion_rms,
            0.01 * suggestion.dispersion_rms_before,
        )
        self.assertLess(actual_orbit_rms, 1e-3)

    def test_sad_style_skew_coupling_correction(self):
        lattice = build_atf_dr_lattice()
        correction = ATFDRRingCorrection(lattice)
        skew_names = correction.get_skew_corrector_names()[:8]
        response = correction.compute_coupling_response(
            probe_corrector_names=("ZH1R", "ZH2R"),
            skew_corrector_names=skew_names,
            bpm_names=correction.bpm_names[:16],
            probe_perturbation=1e-4,
            skew_perturbation=1e-5,
        )

        # A known thin skew-K1L error creates the vertical response measured
        # by SAD skewcor.n after its two horizontal-steerer probe changes.
        correction.set_skew_strength(skew_names[0], 1e-3)
        suggestion = correction.propose_coupling_correction(
            response,
            max_skew_correctors=6,
        )
        corrected_signal = correction.apply_coupling_correction(suggestion)
        actual_rms = float(np.sqrt(np.mean(corrected_signal**2)))

        self.assertLessEqual(len(suggestion.selected_skew_correctors), 6)
        self.assertLess(
            suggestion.rms_predicted,
            0.01 * suggestion.rms_before,
        )
        self.assertLess(actual_rms, 0.001 * suggestion.rms_before)


if __name__ == "__main__":
    unittest.main()
