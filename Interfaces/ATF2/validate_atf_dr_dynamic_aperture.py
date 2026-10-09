"""Fast, offline regression checks for acceptance-scan bookkeeping."""

from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

from Interfaces.ATF2 import scan_atf_dr_dynamic_aperture as scan


class FakeBunch:
    def __init__(self, mass, charge, species, coordinates):
        self.coordinates = np.asarray(coordinates).reshape(-1, 6).copy()

    def size(self):
        return len(self.coordinates)

    def get_phase_space(self):
        return self.coordinates


class FakeLattice:
    def track(self, bunch):
        z = bunch.coordinates.copy()
        z[:, 0] += 1
        if z[0, 0] > 2:
            z = np.empty((0, 6))
        return FakeBunch(None, None, None, z)


def machine():
    return SimpleNamespace(
        dr_first_turn_kicks_rad={}, dr_closed_orbit=np.zeros(4),
        dr_synchronous_orbit=np.array([0., 0., 0., 0., 7., 1300.]),
        dr_momentum_mev_c=1300., dr_start_twiss=(2., 1., 4., -2.),
        dr_start_dispersion_mm_mrad=np.array([10., 1., 2., 3.]),
        dr_lattice=FakeLattice(),
    )


class AcceptanceTests(unittest.TestCase):
    def test_phase_and_dispersion_units(self):
        z = scan.make_scan_coordinates(machine(), radius_mm=2., ray_deg=0.,
                                       phase_x_deg=0., phase_y_deg=90., delta_p=.01,
                                       time_offset_mm_c=3.)
        np.testing.assert_allclose(z, [2.1, -.99, .02, .03, 10., 1313.])

    @patch.object(scan.rft, "Bunch6d", FakeBunch)
    def test_loss_turn_and_independent_particles(self):
        m = machine()
        z = np.array([0., 0., 0., 0., 7., 1300.])
        original = z.copy()
        result = scan.track_probe(m, z, turns=4)
        self.assertEqual(result["first_rejected_turn"], 3)
        self.assertEqual(result["status"], "rftrack_loss")
        np.testing.assert_array_equal(z, original)
        injected = scan.evaluate_injected_phase_space(m, [z, z + [-2., 0., 0., 0., 0., 0.]], turns=4)
        self.assertEqual(injected["equal_weight_survival_fraction"], .5)
        self.assertEqual(scan.evaluate_injected_phase_space(m, np.empty((0, 6)), turns=4)["equal_weight_survival_fraction"], None)

    def test_guard_is_not_physical_loss(self):
        result = scan.track_probe(machine(), [101., 0., 0., 0., 0., 1300.], turns=1)
        self.assertEqual(result["status"], "numerical_guard")
        self.assertEqual(result["first_rejected_turn"], 0)

    def test_no_false_boundary_and_preserve_islands(self):
        def mock_track(m, z, **kwargs):
            return {"survived": not np.isclose(z[0], 1.)}
        with patch.object(scan, "track_probe", mock_track):
            result = scan.scan_acceptance(machine(), radii_mm=[0., 1., 2.],
                                          rays_deg=[0., 180.], phases_deg=[0.],
                                          delta_ps=[0.], turns=2)
        restricted, open_ray = result["sampled_boundaries"]
        self.assertEqual(restricted["last_contiguous_surviving_radius_mm"], 0.)
        self.assertEqual(restricted["surviving_samples_after_first_failure"], 1)
        self.assertTrue(open_ray["boundary_beyond_scan"])
        self.assertIsNone(open_ray["first_failing_radius_mm"])

    def test_invalid_inputs_and_pulses(self):
        for turns in (0, 1.5, True):
            with self.assertRaises(ValueError):
                scan.track_probe(machine(), [0., 0., 0., 0., 0., 1300.], turns=turns)
        m = machine()
        m.dr_first_turn_kicks_rad = {"KII.1": (1e-6, 0)}
        with self.assertRaises(ValueError):
            scan.track_probe(m, [0., 0., 0., 0., 0., 1300.], turns=1)


if __name__ == "__main__":
    unittest.main()
