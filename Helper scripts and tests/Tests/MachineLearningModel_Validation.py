import numpy as np
from pathlib import Path
import sys, os
project_root_path = Path.cwd().resolve()
while not (project_root_path / "Interfaces").exists() and project_root_path.parent != project_root_path:
    project_root_path = project_root_path.parent
sys.path.insert(0, str(project_root_path))
os.chdir(project_root_path)

from Interfaces.ATF2.InterfaceATF2_Ext_RFTrack import InterfaceATF2_Ext_RFTrack
from MachineLearning.ML_train import MLInterface


QUAD_NAME = "QF17X"
SCREENS = ["OTR0X", "OTR1X", "OTR2X", "OTR3X"]
REFERENCE_SCREEN = "OTR0X"

# QF17X
K1_VALUES = np.array([0.9088758826255798, 0.9593689441680908, 1.0098620653152466, 1.060355305671692, 1.1108481884002686], dtype=float)

# QD18X
# K1_VALUES = np.array( [-0.5077778100967407, -0.5924074053764343, -0.6770371198654175, -0.7616666555404663, -0.8462962508201599], dtype=float)

# # QF17X real machine solution
# RFT_SOLUTION = {
#     "emit_x": 7.8677,
#     "beta_x0": 5.78,
#     "alpha_x0": -2.64 ,
#     "emit_y": 0.388,
#     "beta_y0": 0.639,
#     "alpha_y0": 0.02945,
# }

# QF17X simulation solution
RFT_SOLUTION = {
    "emit_x": 5.2,
    "beta_x0": 4.868043658,
    "alpha_x0": 0.001450331482 ,
    "emit_y": 0.03,
    "beta_y0":  0.6796808288,
    "alpha_y0": 0.1019908616,
}

# # QD18X
# RFT_SOLUTION = {
#     "emit_x": 5.2,
#     "beta_x0": 1.316683546,
#     "alpha_x0": -1.36594821 ,
#     "emit_y": 0.03,
#     "beta_y0": 10.70167088,
#     "alpha_y0": -0.4245954529,
# }


def relative_error(prediction, reference):
    denominator = np.maximum(np.abs(reference), 1e-30)
    return np.abs(prediction - reference) / denominator

def print_plane_results(name, ml_sigma, rft_sigma):
    ml_sigma2 = np.asarray(ml_sigma, dtype=float) ** 2
    rft_sigma2 = np.asarray(rft_sigma, dtype=float) ** 2

    error = relative_error(ml_sigma2, rft_sigma2)

    print(f"\n========== {name} ==========")

    for screen_index, screen in enumerate(SCREENS):
        print(f"\n{screen}")

        for k1_index, k1 in enumerate(K1_VALUES):
            print(
                f"K1={k1: .8f} | "
                f"RFTrack sigma={rft_sigma[k1_index, screen_index]:.8g} | "
                f"ML sigma={ml_sigma[k1_index, screen_index]:.8g} | "
                f"error sigma2={100.0 * error[k1_index, screen_index]:.3f}%"
            )

        screen_error = error[:, screen_index]

        print(
            f"SUMMARY {screen}: "
            f"mean={100.0 * np.mean(screen_error):.3f}%, "
            f"median={100.0 * np.median(screen_error):.3f}%, "
            f"max={100.0 * np.max(screen_error):.3f}%"
        )

    print(
        f"\nTOTAL {name}: "
        f"mean={100.0 * np.mean(error):.3f}%, "
        f"median={100.0 * np.median(error):.3f}%, "
        f"p95={100.0 * np.percentile(error, 95):.3f}%, "
        f"max={100.0 * np.max(error):.3f}%"
    )

    return error


def main():
    if K1_VALUES.size == 0:
        raise RuntimeError(
            "K1L_VALUES is empty. Paste the exact K1L values from the session."
        )

    rft_interface = InterfaceATF2_Ext_RFTrack()

    ml_interface = MLInterface(
        interface=rft_interface,
        quad_name=QUAD_NAME,
        screens=SCREENS,
        machine_name="ATF2",
    )

    common_arguments = {
        "quad_name": QUAD_NAME,
        "screens": SCREENS,
        "K1L_values": K1_VALUES,
        "emit_x": RFT_SOLUTION["emit_x"],
        "beta_x0": RFT_SOLUTION["beta_x0"],
        "alpha_x0": RFT_SOLUTION["alpha_x0"],
        "emit_y": RFT_SOLUTION["emit_y"],
        "beta_y0": RFT_SOLUTION["beta_y0"],
        "alpha_y0": RFT_SOLUTION["alpha_y0"],
        "reference_screen": REFERENCE_SCREEN,
    }

    print("Calculating RF-Track reference...")
    rft_sigx, rft_sigy = (
        rft_interface.predict_emittance_scan_response(**common_arguments)
    )

    print("Calculating ML prediction...")
    ml_sigx, ml_sigy = (
        ml_interface.predict_emittance_scan_response(**common_arguments)
    )

    rft_sigx = np.asarray(rft_sigx, dtype=float)
    rft_sigy = np.asarray(rft_sigy, dtype=float)
    ml_sigx = np.asarray(ml_sigx, dtype=float)
    ml_sigy = np.asarray(ml_sigy, dtype=float)

    expected_shape = (len(K1_VALUES), len(SCREENS))

    for name, array in [
        ("rft_sigx", rft_sigx),
        ("rft_sigy", rft_sigy),
        ("ml_sigx", ml_sigx),
        ("ml_sigy", ml_sigy),
    ]:
        if array.shape != expected_shape:
            raise RuntimeError(
                f"{name} has shape {array.shape}, "
                f"expected {expected_shape}"
            )

        if not np.all(np.isfinite(array)):
            raise RuntimeError(f"{name} contains NaN or inf values.")

    error_x = print_plane_results("SIGMA X", ml_sigx, rft_sigx)
    error_y = print_plane_results("SIGMA Y", ml_sigy, rft_sigy)

    worst_error = max(float(np.max(error_x)), float(np.max(error_y)))
    p95_error = float(
        np.percentile(
            np.concatenate([error_x.ravel(), error_y.ravel()]),
            95,
        )
    )

    print("\n========== FINAL VERDICT ==========")

    if worst_error <= 0.03:
        print(
            "PASS: ML differs from RF-Track by at most 3% "
            "in this fitted region."
        )
    elif p95_error <= 0.10:
        print(
            "WARNING: most points are acceptable, but some local "
            "predictions exceed 3% error."
        )
    else:
        print(
            "FAIL: ML significantly deforms the RF-Track response "
            "near the correct fit."
        )

# import numpy as np
#
# d = np.load("/Users/wiktoriamalek/Desktop/flight-simulator/MachineLearning/ATF2/QF17X/all_screens/EM_dataset.npz")
#
# Y = d["Y"]
#
# print(Y.min(axis=0))
# print(Y.max(axis=0))

if __name__ == "__main__":
    main()
