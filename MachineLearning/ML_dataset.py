'''
Learns what is currently done by interface.predict_emittance_scan_response()

Input parameters per sample:
emit_x_norm, beta_x0, alpha_x0,
emit_y_norm, beta_y0, alpha_y0,
K1L

Output targets per sample:
sigx on selected screens + sigy on selected screens
'''

import sys, time
from pathlib import Path
from Interfaces.interface_setup import INTERFACE_SETUP
import numpy as np
import RF_Track as rft

THIS_FILE = Path(__file__).resolve()
PROJECT_ROOT = THIS_FILE.parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

OUTPUT_FILE = PROJECT_ROOT / "MachineLearning" / "EM_dataset.npz"
N_SAMPLES = 2000
RANDOM_SEED = 2137
CHECKPOINT_ROWS = 1000
PARAMETER_NAMES = [
    "emit_x_norm",
    "beta_x0",
    "alpha_x0",
    "emit_y_norm",
    "beta_y0",
    "alpha_y0",
    "K1L",
]

N_K1_PER_TWISS = 7
TARGET_OBSERVABLE = "sigma"

def _save_dataset_checkpoint(output_file, X, Y, twiss_group_ids, screens, quad_name, k1l_relative_change, bounds, interface, K1L_nominal, n_samples, n_twiss_combinations, n_k1l_per_twiss_set):
    output_file = Path(output_file)
    tmp_file = Path(output_file).with_name("tmp_dataset.npz")
    np.savez(
        tmp_file,
        X=np.asarray(X, dtype=float),
        Y=np.asarray(Y, dtype=float),
        twiss_group_ids=np.asarray(twiss_group_ids, dtype=int),
        param_names=np.array(PARAMETER_NAMES),
        sigma_names=np.array([f"sigx_{screen}" for screen in screens] + [f"sigy_{screen}" for screen in screens]),
        target_observable=np.array(TARGET_OBSERVABLE),
        screens=np.array(screens),
        quad_name=np.array(quad_name),
        reference_screen=np.array(screens[0]),
        K1L_relative_change=np.array(k1l_relative_change, dtype=float),
        bounds=np.array([bounds[name] for name in PARAMETER_NAMES[:-1]], dtype=float),
        bounds_names=np.array(PARAMETER_NAMES[:-1]),
        interface_class_name=np.array(interface.__class__.__name__),
        interface_module=np.array(interface.__class__.__module__),
        K1L_nominal=np.array(K1L_nominal, dtype=float),
        n_requested_samples=np.array(n_samples, dtype=int),
        n_twiss_combinations=np.array(n_twiss_combinations, dtype=int),
        n_k1l_per_twiss_set=np.array(n_k1l_per_twiss_set, dtype=int),
    )
    tmp_file.replace(output_file)


def _get_interface_initial_settings(interface):
    interface_class_name = interface.__class__.__name__
    interface_module_name = interface.__class__.__module__

    for machine_name, machine_interfaces in INTERFACE_SETUP.items():
        for interface_defaults in machine_interfaces:
            if (interface_defaults.get("class_name") == interface_class_name) and (
                    interface_defaults.get("module") == interface_module_name):
                result = dict(interface_defaults)
                result["machine"] = machine_name
                return result
    return None

def _get_interface_bounds(interface):
    interface_defaults =_get_interface_initial_settings(interface)
    if interface_defaults is None:
        return {}
    return dict(interface_defaults.get("bounds", {}))

def get_nominal_K1L(interface, quad_name):
    quads = interface.get_quadrupoles()
    names = list(quads["names"])
    strengths = np.asarray(quads["bdes"], dtype=float)
    if quad_name not in names:
        raise RuntimeError(f"Quad {quad_name} not present in interface. Available examples: {names[:10]}")
    return float(strengths[names.index(quad_name)])

def build_twiss_samples_quasirandom(bounds, n_samples):
    twiss_names = PARAMETER_NAMES[:6]
    lower = np.array([bounds[name][0] for name in twiss_names], dtype=float)
    upper = np.array([bounds[name][1] for name in twiss_names], dtype=float)
    if np.any(upper <= lower):
        raise ValueError("Every Twiss sampling bound must satisfy lower < upper.")
    unit_points = np.asarray(rft.qrand(int(n_samples), 6), dtype=float)
    return lower + unit_points * (upper - lower)

def build_k1l_set(rng, K1L_nominal, relative_k_change, number_of_k1l_per_twiss_set, jitter_fraction=0.05):
    delta_min = float(relative_k_change[0])
    delta_max = float(relative_k_change[1])
    number_of_k1l_per_twiss_set = int(number_of_k1l_per_twiss_set)

    if number_of_k1l_per_twiss_set <= 1:
        relative_values = np.array([float(rng.uniform(delta_min, delta_max))], dtype=float)
    else:
        relative_values = np.linspace(delta_min, delta_max, number_of_k1l_per_twiss_set, dtype=float)
        if jitter_fraction > 0:
            step = (delta_max - delta_min) / max(number_of_k1l_per_twiss_set - 1, 1)
            relative_values += rng.uniform(-jitter_fraction * step, jitter_fraction * step, size = relative_values.shape)
            relative_values = np.clip(relative_values, delta_min, delta_max)
            relative_values[0] = delta_min
            relative_values[-1] = delta_max
    return K1L_nominal * (1.0 + relative_values)


def append_dataset(quad_name, screens, interface, k1l_relative_change, n_samples, existing_file, new_bounds):
    existing = np.load(existing_file, allow_pickle=True)
    tmp_file = Path(existing_file).with_name("tmp_dataset.npz")
    generate_dataset(quad_name=quad_name, screens=screens, interface=interface, k1l_relative_change=k1l_relative_change, n_samples=n_samples, output_file=tmp_file, log_callback=print, progress_callback=None, stop_checker=None, new_bounds=new_bounds)
    new = np.load(tmp_file, allow_pickle=True)

    existing_target = str(existing["target_observable"].item()) if "target_observable" in existing.files else "sigma2"
    if existing_target not in {"sigma", "sigma2"}:
        raise RuntimeError(f"Unknown target observable in {existing_file}: {existing_target}")

    existing_y = np.asarray(existing["Y"], dtype=float)
    if existing_target == "sigma2":
        existing_y = np.sqrt(np.maximum(existing_y, 0.0))

    X = np.concatenate([existing["X"], new["X"]])
    Y = np.concatenate([existing_y, new["Y"]])
    existing_screens = [str(screen) for screen in existing["screens"]]

    np.savez(
        existing_file,
        X=X,
        Y=Y,
        param_names=existing["param_names"],
        sigma_names=np.array([f"sigx_{screen}" for screen in existing_screens] + [f"sigy_{screen}" for screen in existing_screens]),
        target_observable=np.array(TARGET_OBSERVABLE),
        screens=existing["screens"],
        quad_name=existing["quad_name"],
        reference_screen=existing["reference_screen"],
        K1L_relative_change=existing["K1L_relative_change"],
        bounds=existing["bounds"],
        bounds_names=existing["bounds_names"],
        interface_class_name=existing["interface_class_name"],
        interface_module=existing["interface_module"],
        K1L_nominal=existing["K1L_nominal"],
        n_requested_samples=len(X),
        n_twiss_combinations=len(X) // int(existing["n_k1l_per_twiss_set"]),
        n_k1l_per_twiss_set=existing["n_k1l_per_twiss_set"],
    )

    tmp_file.unlink()

    return {
        "old_rows": int(len(existing["X"])),
        "new_rows": int(len(new["X"])),
        "total_rows": int(len(X)),
        "output_file": str(existing_file),
    }

def generate_dataset(quad_name, screens, interface, k1l_relative_change, n_samples, output_file, log_callback, progress_callback, stop_checker, new_bounds = None):

    rng = np.random.default_rng(RANDOM_SEED)
    log = log_callback if callable(log_callback) else print
    print("RF-Track quasi-random diagnostic: rft.qrand(100, 6)")
    print(rft.qrand(100, 6))
    n_samples = int(n_samples)
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    K1L_nominal = get_nominal_K1L(interface, quad_name)

    if new_bounds is not None:
        bounds = dict(new_bounds)
        log = log_callback if callable(log_callback) else print
    else:
        bounds = _get_interface_bounds(interface)

    n_k1l_per_twiss_set = min(N_K1_PER_TWISS, max(1, n_samples))
    n_twiss_combinations = max(1, int(np.ceil(n_samples / n_k1l_per_twiss_set)))
    n_total = n_k1l_per_twiss_set * n_twiss_combinations
    twiss_samples = build_twiss_samples_quasirandom(bounds, n_twiss_combinations)

    log(f"{quad_name} nominal K1L: {K1L_nominal}")
    log(f"Generating {n_samples} rows with {n_twiss_combinations} twiss combinations and {n_k1l_per_twiss_set} K1Ls per sample.")
    log("Twiss sampling: 6D RF-Track quasi-random sequence.")
    log(f"Screens: {screens}")
    log(f"K1L relative change: {k1l_relative_change}")
    log(f"Output file: {output_file}")

    twiss_group_ids = []

    X = []
    Y = []
    failed = 0
    processed_rows = 0
    last_checkpoint_rows = 0
    t0 = time.perf_counter()

    for twiss_combination in range(n_twiss_combinations):
        if callable(stop_checker) and stop_checker():
            log("Dataset generation stopped by user.")
            break
        twiss_parameters = twiss_samples[twiss_combination]
        emit_x_norm, beta_x0, alpha_x0, emit_y_norm, beta_y0, alpha_y0 = twiss_parameters
        K1L_array = build_k1l_set(rng, K1L_nominal, k1l_relative_change, n_k1l_per_twiss_set)
        try:
            pred_sigx, pred_sigy = interface.predict_emittance_scan_response(
                quad_name=quad_name,
                screens=screens,
                K1L_values=np.asarray(K1L_array, dtype=float),
                emit_x=emit_x_norm,
                emit_y=emit_y_norm,
                beta_x0=beta_x0,
                alpha_x0=alpha_x0,
                beta_y0=beta_y0,
                alpha_y0=alpha_y0,
                reference_screen=screens[0],
                stop_checker=None,
            )
            required_shape = (len(K1L_array), len(screens))
            if pred_sigx.shape != required_shape or pred_sigy.shape != required_shape:
                raise RuntimeError(f"Shapes are mismatched.")

            for k_idx, K1L in enumerate(K1L_array):
                parameters = np.array([emit_x_norm, beta_x0, alpha_x0, emit_y_norm, beta_y0, alpha_y0, float(K1L)], dtype=float)
                sigma = np.concatenate([pred_sigx[k_idx,:], pred_sigy[k_idx,:]]).astype(float)

                if not (np.all(np.isfinite(parameters)) and np.all(np.isfinite(sigma))):
                    continue

                X.append(parameters)
                Y.append(sigma)
                twiss_group_ids.append(int(twiss_combination))
                processed_rows += 1

                if processed_rows - last_checkpoint_rows >= CHECKPOINT_ROWS:
                    _save_dataset_checkpoint(
                        output_file=output_file,
                        X=X,
                        Y=Y,
                        twiss_group_ids=twiss_group_ids,
                        screens=screens,
                        quad_name=quad_name,
                        k1l_relative_change=k1l_relative_change,
                        bounds=bounds,
                        interface=interface,
                        K1L_nominal=K1L_nominal,
                        n_samples=n_samples,
                        n_twiss_combinations=n_twiss_combinations,
                        n_k1l_per_twiss_set=n_k1l_per_twiss_set,
                    )

                    last_checkpoint_rows = processed_rows

                    log(
                        f"Checkpoint saved: "
                        f"{processed_rows} valid rows going to {output_file}"
                    )

        except Exception as e:
            failed += 1
            log(f"Twiss combination {twiss_combination+1} failed due to {e}")

        if callable(progress_callback):
            progress_callback(min(processed_rows, n_samples), n_samples)
        if (twiss_combination + 1) % 50 == 0 or twiss_combination == 0:
            elapsed = time.perf_counter() - t0
            log(f"Elapsed time: {elapsed:.2f} seconds = {elapsed / 60:.2f} minutes = {elapsed / 3600} hours")
            log(f"{twiss_combination + 1}/{n_twiss_combinations} Twiss combinations processed")
            log(f"{len(X)} samples are valid")
            log(f"{failed} samples failed")

    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    twiss_group_ids = np.asarray(twiss_group_ids, dtype=int)

    if X.size == 0 or Y.size == 0:
        raise RuntimeError("X and Y are empty")

    _save_dataset_checkpoint(
        output_file=output_file,
        X=X,
        Y=Y,
        twiss_group_ids=twiss_group_ids,
        screens=screens,
        quad_name=quad_name,
        k1l_relative_change=k1l_relative_change,
        bounds=bounds,
        interface=interface,
        K1L_nominal=K1L_nominal,
        n_samples=n_samples,
        n_twiss_combinations=n_twiss_combinations,
        n_k1l_per_twiss_set=n_k1l_per_twiss_set,
    )

    log("Done.")
    log(f"Saved dataset to: {output_file}")
    log(f"X shape: {X.shape}")
    log(f"Y shape: {Y.shape}")
    log(f"Failed Twiss combinations: {failed}")

    return {
        "output_file": str(output_file),
        "X_shape": X.shape,
        "Y_shape": Y.shape,
        "failed": int(failed),
        "valid": int(X.shape[0]),
        "n_twiss_combinations": int(n_twiss_combinations),
        "n_k1l_per_twiss_set": int(n_k1l_per_twiss_set),
    }
