"""Synthetic validation of the BT-IPZT to DR-KII.1 handoff fitter.

This is a software-only recovery test.  The hidden map is deliberately
different from the SAD design baseline and the samples are synthetic; replace
them with BPM/profile/injection observations for a real digital-twin fit.
"""

from __future__ import annotations

import json
from dataclasses import replace

import numpy as np

from Interfaces.ATF2.ATF2_LinacBTDR_RFTrack import ATF2LinacBTDRRFTrack
from Interfaces.ATF2.linac_bt_dr_handoff_fit import (
    evaluate_transverse_handoff,
    fit_transverse_handoff,
)
from Interfaces.ATF2.linac_bt_dr_model_config import (
    MODEL_ENERGY_MATCHED_CAVITY_VOLTAGE_MV,
)
from Interfaces.ATF2.simulate_linac_bt_dr_handoff import TransverseHandoff


def main() -> None:
    baseline = TransverseHandoff.twiss_dispersion_matched(
        np.zeros(4), np.zeros(4),
        source_twiss=(2.0, 0.1, 3.0, -0.2),
        target_twiss=(4.0, -0.3, 5.0, 0.4),
        source_dispersion_mm_mrad=np.array((2.0, 0.5, 0.0, 0.0)),
        target_dispersion_mm_mrad=np.array((6.0, 1.0, 0.0, 0.0)),
        reference_momentum_mev_c=1295.524178,
    )
    source = np.array([
        (-0.8, -0.2, -0.5, -0.1), (0.8, -0.2, -0.5, -0.1),
        (-0.8, 0.2, -0.5, -0.1), (-0.8, -0.2, 0.5, -0.1),
        (-0.8, -0.2, -0.5, 0.1), (0.6, 0.15, 0.45, 0.08),
        (-0.4, 0.12, -0.35, 0.06),
    ], dtype=float)
    hidden_matrix = baseline.matrix + np.array([
        (0.015, -0.020, 0.0, 0.0), (0.010, 0.005, 0.0, 0.0),
        (0.0, 0.0, -0.010, 0.012), (0.0, 0.0, -0.006, 0.008),
    ])
    hidden_offset = baseline.offset_mm_mrad + np.array((0.12, -0.03, -0.08, 0.02))
    observed = source @ hidden_matrix.T + hidden_offset
    result = fit_transverse_handoff(
        source, observed, baseline_dispersion_handoff=baseline,
        provenance="synthetic hidden-map recovery; not machine data",
    )
    wrong_boundary_rejected = False
    try:
        ATF2LinacBTDRRFTrack(
            cavity_voltage_mv=MODEL_ENERGY_MATCHED_CAVITY_VOLTAGE_MV,
            handoff=replace(result.handoff, target_location="RING0$START"),
            handoff_mode="sad_optics_matched",
            dr_rf_mode="equilibrium",
        )
    except ValueError as error:
        wrong_boundary_rejected = "target_location" in str(error)
    if not wrong_boundary_rejected:
        raise AssertionError("a handoff referenced to RING0$START was accepted at KII.1")
    held_out_source = np.array([
        (0.35, -0.07, -0.15, 0.03),
        (-0.55, 0.19, 0.25, -0.05),
    ], dtype=float)
    held_out_target = held_out_source @ hidden_matrix.T + hidden_offset
    held_out = evaluate_transverse_handoff(
        result.handoff, held_out_source, held_out_target
    )
    if held_out.max_abs_residual_mm_or_mrad > 1.0e-12:
        raise AssertionError("recovered handoff failed independent synthetic validation")
    print(json.dumps({
        "simulation_only": True,
        "fit": result.as_dict(),
        "max_abs_matrix_error": float(np.max(np.abs(result.handoff.matrix - hidden_matrix))),
        "max_abs_offset_error_mm_mrad": float(np.max(np.abs(result.handoff.offset_mm_mrad - hidden_offset))),
        "dispersion_retained": result.handoff.source_dispersion_mm_mrad is not None,
        "wrong_boundary_handoff_rejected": wrong_boundary_rejected,
        "held_out_validation": held_out.as_dict(),
        "real_data_required": [
            "BT endpoint coordinate reconstruction and its uncertainty",
            "DR first-turn/injection orbit observations for independent dithers",
            "injection-kicker pulse timing/amplitude calibration",
            "momentum/dispersion and RF time-of-flight calibration",
        ],
    }, indent=2))


if __name__ == "__main__":
    main()
