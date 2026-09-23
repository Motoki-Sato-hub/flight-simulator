"""Regression checks for the historical ATF DR SAD-to-RF-Track exports."""

from __future__ import annotations

from pathlib import Path

from Interfaces.ATF2.DR_ATF2.ATF_DR_RFTrack_lattice import (
    get_magnet_centres,
    load_lattice_data,
)


ROOT = Path(__file__).resolve().parents[4]
DATA_2008 = ROOT / "Interfaces/ATF2/DR_ATF2/ATF_DR_20080526b_RFTrack_lattice.json"
DATA_2011 = ROOT / "Interfaces/ATF2/DR_ATF2/ATF_DR_RFTrack_lattice.json"


def main():
    historical = load_lattice_data(DATA_2008)
    reference = load_lattice_data(DATA_2011)
    # The 2008 source has ``LMZH12R1=(L=2.2 - 1.06933)``.  This assertion
    # protects the arithmetic-expression parser from silently treating it as
    # 2.2 m, which would move the RF harmonic by more than two buckets.
    assert abs(historical["definitions"]["LMZH12R1"]["attributes"]["L"] - 1.13067) < 1e-12
    assert abs(historical["metadata"]["circumference_m"] - 138.55981) < 1e-8
    assert abs(reference["metadata"]["circumference_m"] - 138.55953941176458) < 1e-8
    assert len(get_magnet_centres(DATA_2008)) == 204
    assert len(get_magnet_centres(DATA_2011)) == 204
    print("Historical SAD export regression checks passed")


if __name__ == "__main__":
    main()
