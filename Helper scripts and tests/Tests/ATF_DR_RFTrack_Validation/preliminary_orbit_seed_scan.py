"""Scan seed dependence of Kubo's rough-COD / periodic-orbit preparation.

This is deliberately a preparation-stage diagnostic, not a correction or
emittance result.  It records whether RF-Track can retain a periodic orbit
while Table-I-strength errors are introduced and the all-steerer rough-COD
search is performed.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import kubo_style_procedure as kubo


FLIGHT_SIMULATOR = Path(__file__).resolve().parents[3]
ANALYSIS = FLIGHT_SIMULATOR.parent / "analysis" / "DR-RFTrack"


def _parse_seeds(text):
    if ":" in text:
        first, last = (int(value) for value in text.split(":", 1))
        if last < first:
            raise ValueError("seed range must be ascending")
        return tuple(range(first, last + 1))
    return tuple(int(value.strip()) for value in text.split(",") if value.strip())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seeds", default="2000:2009")
    parser.add_argument("--continuation-x-mm", type=float, default=0.5)
    parser.add_argument("--continuation-y-mm", type=float, default=0.5)
    args = parser.parse_args()
    if args.continuation_x_mm <= 0 or args.continuation_y_mm <= 0:
        raise ValueError("continuation limits must be positive")

    # The nominal error-free response remains the Kubo convention and is
    # intentionally shared by all seeds, just as in the paper.
    probe_label = f"{kubo.STEERER_RESPONSE_KICK_RAD * 1e3:g}".replace(".", "p")
    cache = ANALYSIS / (
        f"kubo2003_nominal_kick_response_cache_{kubo.LATTICE_LABEL}_"
        f"{kubo.SKEW_FAMILY}_probe_{probe_label}mrad.npz"
    )
    response = kubo._load_or_build_responses(cache)
    kubo.PRELIMINARY_CONTINUATION_X_MM = args.continuation_x_mm
    kubo.PRELIMINARY_CONTINUATION_Y_MM = args.continuation_y_mm
    records = []

    for seed in _parse_seeds(args.seeds):
        kubo.SEED = seed
        try:
            _, _, _, _, _, history, x_commands, y_commands = kubo._make_machine(response)
            final = history[-1]
            records.append({
                "seed": seed,
                "status": "periodic orbit retained",
                "final_error_fraction": final["error_fraction"],
                "final_max_measured_x_mm": final["max_measured_x_mm"],
                "final_max_measured_y_mm": final["max_measured_y_mm"],
                "n_ramp_records": len(history),
                "x_command_rms_mrad": float((x_commands @ x_commands / len(x_commands)) ** 0.5 * 1e3),
                "y_command_rms_mrad": float((y_commands @ y_commands / len(y_commands)) ** 0.5 * 1e3),
            })
        except Exception as error:  # failure itself is the observable here
            records.append({
                "seed": seed,
                "status": "failed",
                "exception": f"{type(error).__name__}: {error}",
            })
        print(records[-1], flush=True)

    success = [item for item in records if item["status"] == "periodic orbit retained"]
    result = {
        "purpose": "Seed dependence of RF-Track Kubo rough-COD preparation",
        "lattice": kubo.LATTICE_LABEL,
        "magnet_error_scale_of_kubo_table_i": kubo.MAGNET_ERROR_SCALE,
        "published_rough_cod_target_mm": {"x": 2.0, "y": 1.0},
        "continuation_trigger_mm": {
            "x": args.continuation_x_mm,
            "y": args.continuation_y_mm,
        },
        "n_seeds": len(records),
        "n_periodic_orbit_retained": len(success),
        "records": records,
        "interpretation": (
            "A failure is a periodic-orbit/preparation failure for this RF-Track "
            "implementation. It is not an emittance or final-correction result."
        ),
    }
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    suffix = f"x{args.continuation_x_mm:g}_y{args.continuation_y_mm:g}".replace(".", "p")
    output = ANALYSIS / f"atf_dr_rftrack_preliminary_seed_scan_{suffix}.json"
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
