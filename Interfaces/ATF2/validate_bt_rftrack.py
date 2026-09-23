"""Smoke validation of SAD-derived BT rigidity scaling.

Run from the flight-simulator directory with::

    MPLCONFIGDIR=/tmp/mpl-flight-simulator PYTHONPATH=. \
      /home/motokisato/rftrack-env/bin/python \
      Interfaces/ATF2/validate_bt_rftrack.py

This is a geometric/reference-particle check.  It does not establish that
the historical SAD optical functions, apertures, or live machine calibration
are yet represented accurately enough for digital-twin claims.
"""

from __future__ import annotations

import numpy as np
import RF_Track as rft

from Interfaces.ATF2.InterfaceATF2_BT_RFTrack import (
    InterfaceATF2_BT_RFTrack,
    _default_tfs_path,
)


def _track_reference(lattice: rft.Lattice, momentum_mev_c: float) -> np.ndarray:
    bunch = rft.Bunch6d(
        rft.electronmass,
        1,
        -1,
        np.array([[0.0, 0.0, 0.0, 0.0, 0.0, momentum_mev_c]]),
    )
    return np.asarray(lattice.track(bunch).get_phase_space("%x %xp %y %yp %Pc")).ravel()


def main() -> None:
    # Deliberately show why changing the particle momentum alone is invalid.
    unscaled = _track_reference(rft.Lattice(str(_default_tfs_path())), 1300.0)
    print("unscaled_1300_exit", unscaled.tolist())

    for momentum_mev_c in (1300.0, 1542.282):
        interface = InterfaceATF2_BT_RFTrack(
            nparticles=100,
            momentum_mev_c=momentum_mev_c,
        )
        exit_phase_space = _track_reference(interface.lattice, momentum_mev_c)
        if not np.allclose(exit_phase_space[:4], 0.0, atol=1e-10):
            raise RuntimeError(
                f"scaled {momentum_mev_c:g} MeV/c reference orbit did not close: "
                f"{exit_phase_space.tolist()}"
            )
        if not np.isclose(interface.get_transmission(), 1.0):
            raise RuntimeError("nominal BT bunch did not have unit transmission")
        print(
            "scaled",
            momentum_mev_c,
            "exit",
            exit_phase_space.tolist(),
            "transmission",
            interface.get_transmission(),
        )


if __name__ == "__main__":
    main()
