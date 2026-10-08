import sys, time
from pathlib import Path
THIS_FILE = Path(__file__).resolve()
PROJECT_ROOT = THIS_FILE.parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
from Interfaces.interface_setup import INTERFACE_SETUP
import numpy as np

RANDOM_SEED = 2137

PARAMETER_NAMES = [
    "emit_x_norm",
    "beta_x0",
    "alpha_x0",
    "emit_y_norm",
    "beta_y0",
    "alpha_y0",
]

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

def samples_uniform(rng, bounds, name):
    low, high = bounds[name]
    return float(rng.uniform(low, high))

def get_nominal_K1L(interface, quad_name):
    quads = interface.get_quadrupoles()
    names = list(quads["names"])
    strengths = np.asarray(quads["bdes"], dtype=float)
    if quad_name not in names:
        raise RuntimeError(f"Quad {quad_name} not present in interface. Available examples: {names[:10]}")
    return float(strengths[names.index(quad_name)])

def build_sample(rng, bounds):
    emit_x_norm = samples_uniform(rng, bounds, "emit_x_norm")
    beta_x0 = samples_uniform(rng, bounds, "beta_x0")
    alpha_x0 = samples_uniform(rng, bounds, "alpha_x0")
    emit_y_norm = samples_uniform(rng, bounds, "emit_y_norm")
    beta_y0 = samples_uniform(rng, bounds, "beta_y0")
    alpha_y0 = samples_uniform(rng, bounds, "alpha_y0")

    return np.array([
        emit_x_norm, beta_x0, alpha_x0,
        emit_y_norm, beta_y0, alpha_y0,
    ], dtype=float)

def generate_dataset(quad_name, screens, interface, k1l_relative_change, n_samples, output_file, log_callback, progress_callback, stop_checker):
    rng = np.random.default_rng(RANDOM_SEED)
    log = log_callback if callable(log_callback) else print
    n_samples = int(n_samples)
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    K1L_nominal = get_nominal_K1L(interface, quad_name)
    bounds = _get_interface_bounds(interface)
    gamma_rel, beta_rel, beta_gamma = interface.get_beam_factors()

    if not np.isfinite(beta_gamma) or beta_gamma <= 0:
        raise RuntimeError(f"Invalid beta_gamma from interface: {beta_gamma}")
    fixed_k1l = float(K1L_nominal)

    log(f"{quad_name} nominal/fixed K1L: {fixed_k1l}")
    log(f"Generating {n_samples} rows with fixed lattice and fixed K1L.")
    log(f"Only incoming beam parameters are varied.")
    log(f"Screens: {screens}")
    log(f"Output file: {output_file}")

    twiss_group_ids = []

    X = []
    Y = []
    failed = 0
    processed_rows = 0
    t0 = time.perf_counter()

    for sample_idx in range(n_samples):
        if callable(stop_checker) and stop_checker():
            log("Dataset generation stopped by user.")
            break
        twiss_parameters = build_sample(rng, bounds)
        emit_x_norm, beta_x0, alpha_x0, emit_y_norm, beta_y0, alpha_y0 = twiss_parameters
        K1L_array = np.array([fixed_k1l], dtype=float)
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

            sigmas_interleaved = []
            for screen_idx in range(len(screens)):
                sigmas_interleaved.append(float(pred_sigx[0, screen_idx]))
                sigmas_interleaved.append(float(pred_sigy[0, screen_idx]))
            sigmas_interleaved = np.asarray(sigmas_interleaved, dtype=float)

            if not (np.all(np.isfinite(twiss_parameters)) and np.all(np.isfinite(sigmas_interleaved))):
                continue

            X.append(twiss_parameters)
            Y.append(sigmas_interleaved)
            twiss_group_ids.append(int(sample_idx))
            processed_rows += 1

        except Exception as e:
            failed += 1
            log(f"Sample {sample_idx + 1} failed due to {e}")

        if callable(progress_callback):
            progress_callback(min(processed_rows, n_samples), n_samples)
        if (sample_idx + 1) % 50 == 0 or sample_idx == 0:
            elapsed = time.perf_counter() - t0
            log(f"Elapsed time: {elapsed:.2f} seconds = {elapsed / 60:.2f} minutes = {elapsed / 3600} hours")
            log(f"{sample_idx + 1}/{n_samples} samples processed")
            log(f"{len(X)} samples are valid")
            log(f"{failed} samples failed")

    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    twiss_group_ids = np.asarray(twiss_group_ids, dtype=int)

    if X.size == 0 or Y.size == 0:
        raise RuntimeError("X and Y are empty")

    np.savez(
        output_file,
        X=X,
        Y=Y,
        param_names=np.array(PARAMETER_NAMES),
        sigma_names=np.array([name for screen in screens for name in (f"sigx_{screen}", f"sigy_{screen}")]),
        screens=np.array(screens),
        quad_name=np.array(quad_name),
        reference_screen=np.array(screens[0]),
        K1L_fixed=np.array(fixed_k1l, dtype=float),
        K1L_relative_change=np.array(k1l_relative_change, dtype=float),
        bounds=np.array([bounds[name] for name in PARAMETER_NAMES], dtype=float),
        bounds_names=np.array(PARAMETER_NAMES),
        interface_class_name=np.array(interface.__class__.__name__),
        interface_module=np.array(interface.__class__.__module__),
        K1L_nominal=np.array(K1L_nominal, dtype=float),
        n_requested_samples = np.array(n_samples, dtype=int),
        beta_gamma=np.array(beta_gamma, dtype=float),
        gamma_rel=np.array(gamma_rel, dtype=float),
        beta_rel=np.array(beta_rel, dtype=float),
    )

    log("Done.")
    log(f"Saved dataset to: {output_file}")
    log(f"X shape: {X.shape}")
    log(f"Y shape: {Y.shape}")
    log(f"Failed samples: {failed}")

    return {
        "output_file": str(output_file),
        "X_shape": X.shape,
        "Y_shape": Y.shape,
        "failed": int(failed),
        "valid": int(X.shape[0]),
        "K1L_fixed": float(fixed_k1l),
    }
