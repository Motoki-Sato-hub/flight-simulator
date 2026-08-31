import numpy as np
from scipy.optimize import least_squares
import pandas as pd
from Backend.EM_helpers.CheckLinearOptics import estimate_twiss_use_linear_optics_start, CheckLinearOpticsUnavailable

def _finite_diff_jacobian(residual_func, x, low, high, param_scale=None, rel_step=1e-4):
    x = np.asarray(x, dtype=float)
    if param_scale is None:
        scale = np.maximum(np.abs(x), 1.0)
    else:
        scale = np.asarray(param_scale, dtype=float)
        fallback_scale = np.maximum(np.abs(x), 1.0)
        scale = np.where(np.isfinite(scale) & (scale > 0), scale, fallback_scale)
    r0 = np.asarray(residual_func(x), dtype=float)
    J = np.empty((r0.size, x.size), dtype=float)
    for i in range(x.size):
        h = rel_step * scale[i]
        x_step = x.copy()
        if x[i] + h <= high[i]:
            x_step[i] += h
            step = h
        else:
            x_step[i] -= h
            step = -h
        r1 = np.asarray(residual_func(x_step), dtype=float)
        J[:, i] = (r1 - r0) / step
    return J

class OptimizationStopped(Exception):
    def __init__(self, message = "Optimization stopped", solution = None):
        super().__init__(message)
        self.solution = solution

class OptimizationPaused(Exception):
    def __init__(self, message="Optimization paused", solution=None):
        super().__init__(message)
        self.solution = solution

class Optimization:
    def __init__(self, interface, fit_quadrupole_strength=False, fit_quad_offset=False, fit_quad_roll=False, progress_callback=None):
        self.progress_callback = progress_callback
        self.interface = interface
        self._stop_requested = False
        self._pause_requested = False
        self.best_out_so_far = None
        self._last_completed_output = None
        self.fit_quadrupole_strength = bool(fit_quadrupole_strength)
        self.fit_quad_offset = bool(fit_quad_offset)
        self.fit_quad_roll = bool(fit_quad_roll)

    def _emit_progress(self, phase, current, total):
        if self.progress_callback is None:
            return
        self.progress_callback(str(phase), int(current), int(total))

    def request_pause(self):
        self._pause_requested = True

    def clear_pause(self):
        self._pause_requested = False

    def request_stop(self):
        self._stop_requested = True

    def clear_stop(self):
        self._stop_requested = False
        self.best_out_so_far = None
        self._last_completed_output = None
        self._pause_requested = False

    def _calculate_optimalization_errors(self, J, n_params=None, param_scale=None):
        n_default = int(n_params) if n_params is not None else 6

        try:
            J = np.asarray(J, dtype=float) # that's the information about how narrow the minimum region is, it's an array with derivatives
        except Exception:
            return {
                "param_errors": np.full(n_default, np.nan),
                "cov": None,
            }

        npar = J.shape[1] if J.ndim == 2 else n_default

        if param_scale is None:
            scale = np.ones(npar, dtype=float)
        else:
            scale = np.asarray(param_scale, dtype=float)
            scale = np.where(np.isfinite(scale) & (scale > 0), scale, 1.0)

        try:
            J_scaled = J * scale[np.newaxis, :]
            cov_scaled = np.linalg.pinv(J_scaled.T @ J_scaled) # default rcond; narrow region -> small cov -> small errors, wide region -> big cov -> big errors
            cov = cov_scaled * np.outer(scale, scale)
            param_errors = np.sqrt(np.maximum(np.diag(cov), 0.0))
            print("Covariance matrix:")
            print(cov)
        except Exception:
            cov = None
            param_errors = np.full(npar, np.nan)

        return {
            "param_errors": param_errors,
            "cov": cov,
        }

    def fit_from_session(self, session, bounds):
        was_pause_requested = bool(self._pause_requested)
        self.clear_stop()
        self._pause_requested = was_pause_requested
        print("Starting to fit Twiss parameters and emittance...")
        screens = list(session.get("screens", []))
        quad_name = session.get("quad_name")
        K1L_values = np.asarray(session.get("K1L_values", []), dtype=float)
        quad_k1l_0_readback = np.asarray(session.get("K1L_0", K1L_values[len(K1L_values)//2] if K1L_values.size else np.nan))
        sigx = np.asarray(session.get("sigx_mean", []), dtype=float)
        sigy = np.asarray(session.get("sigy_mean", []), dtype=float)
        sigx_shots = np.asarray(session.get("sigx_shots", []), dtype=float)
        sigy_shots = np.asarray(session.get("sigy_shots", []), dtype=float)
        x_shots = np.asarray(session.get("x_shots", []), dtype=float)
        y_shots = np.asarray(session.get("y_shots", []), dtype=float)
        sigxy_shots = np.asarray(session.get("sigxy_shots", []), dtype=float)
        sigma_template_x = np.asarray(sigx if sigx.size else np.empty((0, len(screens))), dtype=float)
        sigma_template_y = np.asarray(sigy if sigy.size else np.empty((0, len(screens))), dtype=float)

        def _plane_no_solution(plane_name, sigma_template):
            return {
                "emit": np.nan,
                "beta0": np.nan,
                "alpha0": np.nan,
                "pred": np.full_like(sigma_template, np.nan, dtype=float),
                "cost": np.nan,
            }
        joint_fit = None
        fit_x = None
        fit_y = None

        try:
            joint_fit = self._fit_6d(screens=screens, quad_name=quad_name, K1L_values=K1L_values, sigx_shots=sigx_shots, sigy_shots=sigy_shots,
                                      x_shots=x_shots, y_shots=y_shots, sigxy_shots=sigxy_shots, bounds = bounds)

            fit_x = {
                "emit": joint_fit["emit_x_geom"],
                "beta0": joint_fit["beta_x0"],
                "alpha0": joint_fit["alpha_x0"],
                "pred": joint_fit["pred_x"],
                "cost": joint_fit["cost"],
            }
            fit_y = {
                "emit": joint_fit["emit_y_geom"],
                "beta0": joint_fit["beta_y0"],
                "alpha0": joint_fit["alpha_y0"],
                "pred": joint_fit["pred_y"],
                "cost": joint_fit["cost"],
            }

            print(
                f"Joint fit done with cost={joint_fit['cost']:.6g}, "
                f"emit_x_geom={joint_fit['emit_x_geom']:.6g}, beta_x0={joint_fit['beta_x0']:.6g}, alpha_x0={joint_fit['alpha_x0']:.6g}, "
                f"emit_y_geom={joint_fit['emit_y_geom']:.6g}, beta_y0={joint_fit['beta_y0']:.6g}, alpha_y0={joint_fit['alpha_y0']:.6g}")

        except (OptimizationStopped, OptimizationPaused) as e:
            if isinstance(getattr(e, "solution", None), dict):
                joint_fit = e.solution
                fit_x = {
                        "emit": joint_fit["emit_x_geom"],
                        "beta0": joint_fit["beta_x0"],
                        "alpha0": joint_fit["alpha_x0"],
                        "pred": joint_fit["pred_x"],
                        "cost": joint_fit["cost"],
                    }
                fit_y = {
                        "emit": joint_fit["emit_y_geom"],
                        "beta0": joint_fit["beta_y0"],
                        "alpha0": joint_fit["alpha_y0"],
                        "pred": joint_fit["pred_y"],
                        "cost": joint_fit["cost"],
                    }
            else:
                fit_x = None
                fit_y = None

        if fit_x is None:
            fit_x = _plane_no_solution("x", sigma_template_x)
        if fit_y is None:
            fit_y = _plane_no_solution("y", sigma_template_y)

        gamma_rel, beta_rel, beta_gamma = self.interface.get_beam_factors()
        emit_x_norm = (beta_gamma * fit_x["emit"] if np.isfinite(fit_x["emit"]) else np.nan)
        emit_y_norm = (beta_gamma * fit_y["emit"] if np.isfinite(fit_y["emit"]) else np.nan)
        emit_x_geom = fit_x["emit"]
        emit_y_geom = fit_y["emit"]

        err_dict = {}
        if isinstance(joint_fit, dict):
            err_dict = dict(joint_fit.get("param_errors") or {})

        emit_x_norm_err = float(err_dict.get("emit_x_norm", np.nan))
        emit_y_norm_err = float(err_dict.get("emit_y_norm", np.nan))
        beta_x0_err = float(err_dict.get("beta_x0", np.nan))
        alpha_x0_err = float(err_dict.get("alpha_x0", np.nan))
        beta_y0_err = float(err_dict.get("beta_y0", np.nan))
        alpha_y0_err = float(err_dict.get("alpha_y0", np.nan))
        quad_k1l_0_err = float(err_dict.get("quad_k1l_0", np.nan))
        quad_dx0_err = float(err_dict.get("quad_dx0", np.nan)) * 1e3 # mm
        quad_dy0_err = float(err_dict.get("quad_dy0", np.nan)) * 1e3 # mm
        quad_roll_err = float(err_dict.get("quad_roll", np.nan)) * 1e3 # mrad

        if np.isfinite(beta_gamma) and beta_gamma > 0:
            emit_x_geom_err = emit_x_norm_err / beta_gamma * 1e3 if np.isfinite(emit_x_norm_err) else np.nan
            emit_y_geom_err = emit_y_norm_err / beta_gamma * 1e3 if np.isfinite(emit_y_norm_err) else np.nan
        else:
            emit_x_geom_err = np.nan
            emit_y_geom_err = np.nan

        result = {
            "screen0": screens[0],
            "quad_name": quad_name,
            "emit_x_norm": emit_x_norm,
            "emit_y_norm": emit_y_norm,
            "emit_x_geom": emit_x_geom * 1e3 , # nm*rad
            "emit_y_geom": emit_y_geom * 1e3, # nm*rad
            "beta_x0": fit_x["beta0"],
            "alpha_x0": fit_x["alpha0"],
            "beta_y0": fit_y["beta0"],
            "alpha_y0": fit_y["alpha0"],
            "fit_x_cost": fit_x["cost"],
            "fit_y_cost": fit_y["cost"],
            "fit_quadrupole_strength": bool(self.fit_quadrupole_strength),
            "fit_quad_offset": bool(self.fit_quad_offset),
            "fit_quad_roll": bool(self.fit_quad_roll),
            "quad_k1l_0": (
                float(joint_fit.get("quad_k1l_0", np.nan))
                if bool(self.fit_quadrupole_strength) and isinstance(joint_fit, dict)
                else quad_k1l_0_readback
            ),
            "quad_dx0": (float(joint_fit.get("quad_dx0", np.nan)) * 1e3 if bool(self.fit_quad_offset) and isinstance(joint_fit, dict) else np.nan),  # mm
            "quad_dy0": (float(joint_fit.get("quad_dy0", np.nan)) * 1e3 if bool(self.fit_quad_offset) and isinstance(joint_fit, dict) else np.nan),  # mm
            "quad_dx0_err": quad_dx0_err,
            "quad_dy0_err": quad_dy0_err,
            "quad_roll": (float(joint_fit.get("quad_roll", np.nan)) * 1e3 if bool(self.fit_quad_roll) and isinstance(joint_fit, dict) else np.nan),  # mrad
            "quad_roll_err": quad_roll_err,
            "emit_x_norm_err": emit_x_norm_err,
            "emit_y_norm_err": emit_y_norm_err,
            "emit_x_geom_err": emit_x_geom_err,
            "emit_y_geom_err": emit_y_geom_err,
            "beta_x0_err": beta_x0_err,
            "alpha_x0_err": alpha_x0_err,
            "beta_y0_err": beta_y0_err,
            "alpha_y0_err": alpha_y0_err,
            "quad_k1l_0_err": quad_k1l_0_err,
        }

        output = {
            "result": result,
            "pred_x": fit_x["pred"],
            "pred_y": fit_y["pred"],
            "screens": list(session.get("screens", [])),
            "K1L_values": K1L_values.tolist(),
        }

        self.best_out_so_far = output
        self._last_completed_output = output
        self._pause_requested = False

        print(
            f"Final result: "
            f"emit_x_norm={result['emit_x_norm']:.6g}, "
            f"emit_y_norm={result['emit_y_norm']:.6g}, "
            f"fit_x_cost={result['fit_x_cost']:.6g}, "
            f"fit_y_cost={result['fit_y_cost']:.6g}"
        )

        return output

    def _build_joint_partial_output(self, screens, sigma_x, sigma_y, pred_x, pred_y, best_row, best_cost):
        if best_row is None or pred_x is None or pred_y is None:
            return None
        gamma_rel, beta_rel, beta_gamma = self.interface.get_beam_factors()
        emit_x_geom = max(float(best_row["emit_x_norm"]) / beta_gamma, 1e-12)
        emit_y_geom = max(float(best_row["emit_y_norm"]) / beta_gamma, 1e-12)

        return {
            "emit_x_geom": emit_x_geom,
            "beta_x0": float(best_row["beta_x0"]),
            "alpha_x0": float(best_row["alpha_x0"]),
            "emit_y_geom": emit_y_geom,
            "beta_y0": float(best_row["beta_y0"]),
            "alpha_y0": float(best_row["alpha_y0"]),
            "pred_x": pred_x,
            "pred_y": pred_y,
            "emit_x_norm": float(best_row["emit_x_norm"]),
            "emit_y_norm": float(best_row["emit_y_norm"]),
            "quad_k1l_0": (float(best_row["quad_k1l_0"]) if "quad_k1l_0" in best_row else np.nan),
            "quad_dx0": (float(best_row["quad_dx0"]) if "quad_dx0" in best_row else np.nan),
            "quad_dy0": (float(best_row["quad_dy0"]) if "quad_dy0" in best_row else np.nan),
            "quad_roll": (float(best_row["quad_roll"]) if "quad_roll" in best_row else np.nan),
            "cost": float(best_cost) if np.isfinite(best_cost) else np.nan,
        }

    def _fit_6d(self, screens, quad_name, K1L_values, sigx_shots, sigy_shots, bounds, x_shots=None, y_shots=None, sigxy_shots=None):

        # Beam size for every individual shot.
        sigma_x_shots = np.asarray(sigx_shots, dtype=float)
        sigma_y_shots = np.asarray(sigy_shots, dtype=float)

        # number of valid shots at each scan point and screen
        n_x_sum = np.sum(np.isfinite(sigma_x_shots), axis=2)
        n_y_sum = np.sum(np.isfinite(sigma_y_shots), axis=2)

        # Measured mean beam size.
        sig_x = np.nanmean(sigma_x_shots, axis=2)
        sig_y = np.nanmean(sigma_y_shots, axis=2)

        # Sample standard deviation of sigma, then standard error of its mean.
        s_sigx = np.nanstd(sigma_x_shots, axis=2, ddof=1)
        s_sigy = np.nanstd(sigma_y_shots, axis=2, ddof=1)

        fallback_u_x = np.maximum(0.08 * np.abs(sig_x), 1e-12)
        fallback_u_y = np.maximum(0.08 * np.abs(sig_y), 1e-12)
        u_x = np.where(n_x_sum >= 2, s_sigx / np.sqrt(n_x_sum), fallback_u_x)
        u_y = np.where(n_y_sum >= 2, s_sigy / np.sqrt(n_y_sum), fallback_u_y)
        # u_x > 0 alone would accept float64 round-off noise (~1e-16) from nanstd on
        # bit-identical shots (a deterministic model) as if it were a real, tiny
        # measurement uncertainty, instead of falling back to the 8% heuristic.
        u_x = np.where(np.isfinite(u_x) & (u_x > 1e-9), u_x, fallback_u_x)
        u_y = np.where(np.isfinite(u_y) & (u_y > 1e-9), u_y, fallback_u_y)

        # Points usable in the weighted least-squares residual vector.
        valid_x = np.isfinite(sig_x) & np.isfinite(u_x) & (u_x > 0) & (n_x_sum >= 1)
        valid_y = np.isfinite(sig_y) & np.isfinite(u_y) & (u_y > 0) & (n_y_sum >= 1)

        if not np.any(valid_x) and not np.any(valid_y):
            raise RuntimeError("No valid sigma_x or sigma_y measurements were available for the fit.")
        fallback_u_dx = 1e-3  # mm
        fallback_u_dy = 1e-3  # mm
        fallback_u_sigxy = 1e-3  # mm^2

        x_shots_arr = np.asarray(x_shots, dtype=float) if x_shots is not None else np.empty((0, 0, 0))
        y_shots_arr = np.asarray(y_shots, dtype=float) if y_shots is not None else np.empty((0, 0, 0))
        sigxy_shots_arr = np.asarray(sigxy_shots, dtype=float) if sigxy_shots is not None else np.empty((0, 0, 0))

        have_dx_data = x_shots_arr.ndim == 3 and x_shots_arr.shape[:2] == sig_x.shape
        have_dy_data = y_shots_arr.ndim == 3 and y_shots_arr.shape[:2] == sig_x.shape
        have_sigxy_data = sigxy_shots_arr.ndim == 3 and sigxy_shots_arr.shape[:2] == sig_x.shape

        if have_dx_data:
            n_dx_sum = np.sum(np.isfinite(x_shots_arr), axis=2)
            dx_meas = np.nanmean(x_shots_arr, axis=2)
            s_dx = np.nanstd(x_shots_arr, axis=2, ddof=1)
            u_dx = np.where(n_dx_sum >= 2, s_dx / np.sqrt(n_dx_sum), fallback_u_dx)
            u_dx = np.where(np.isfinite(u_dx) & (u_dx > 1e-9), u_dx, fallback_u_dx)
            valid_dx = np.isfinite(dx_meas) & np.isfinite(u_dx) & (u_dx > 0) & (n_dx_sum >= 1)
        else:
            dx_meas = np.full(sig_x.shape, np.nan, dtype=float)
            u_dx = np.full(sig_x.shape, fallback_u_dx, dtype=float)
            valid_dx = np.zeros(sig_x.shape, dtype=bool)

        if have_dy_data:
            n_dy_sum = np.sum(np.isfinite(y_shots_arr), axis=2)
            dy_meas = np.nanmean(y_shots_arr, axis=2)
            s_dy = np.nanstd(y_shots_arr, axis=2, ddof=1)
            u_dy = np.where(n_dy_sum >= 2, s_dy / np.sqrt(n_dy_sum), fallback_u_dy)
            u_dy = np.where(np.isfinite(u_dy) & (u_dy > 1e-9), u_dy, fallback_u_dy)
            valid_dy = np.isfinite(dy_meas) & np.isfinite(u_dy) & (u_dy > 0) & (n_dy_sum >= 1)
        else:
            dy_meas = np.full(sig_x.shape, np.nan, dtype=float)
            u_dy = np.full(sig_x.shape, fallback_u_dy, dtype=float)
            valid_dy = np.zeros(sig_x.shape, dtype=bool)

        if have_sigxy_data:
            n_sigxy_sum = np.sum(np.isfinite(sigxy_shots_arr), axis=2)
            sigxy_meas = np.nanmean(sigxy_shots_arr, axis=2)
            s_sigxy = np.nanstd(sigxy_shots_arr, axis=2, ddof=1)
            u_sigxy = np.where(n_sigxy_sum >= 2, s_sigxy / np.sqrt(n_sigxy_sum), fallback_u_sigxy)
            u_sigxy = np.where(np.isfinite(u_sigxy) & (u_sigxy > 1e-9), u_sigxy, fallback_u_sigxy)
            valid_sigxy = np.isfinite(sigxy_meas) & np.isfinite(u_sigxy) & (u_sigxy > 0) & (n_sigxy_sum >= 1)
        else:
            sigxy_meas = np.full(sig_x.shape, np.nan, dtype=float)
            u_sigxy = np.full(sig_x.shape, fallback_u_sigxy, dtype=float)
            valid_sigxy = np.zeros(sig_x.shape, dtype=bool)

        if self.fit_quad_offset and not (np.any(valid_dx) and np.any(valid_dy)):
            print("Fit quad offset is not reliable.")
        if self.fit_quad_roll and not np.any(valid_sigxy):
            print("Fit quad roll is not reliable.")

        gamma_rel, beta_rel, beta_gamma = self.interface.get_beam_factors()
        bounds = dict(bounds or {})
        K1L_values = np.asarray(K1L_values, dtype=float)
        K1L_0_readback = float(K1L_values[len(K1L_values)//2])
        deltas_for_fit = K1L_values / K1L_0_readback - 1.0

        params_order = ["emit_x_norm", "beta_x0", "alpha_x0", "emit_y_norm", "beta_y0", "alpha_y0"]
        n_core_params = len(params_order)
        if self.fit_quadrupole_strength:
            params_order.append("quad_k1l_0")
        if self.fit_quad_offset:
            params_order.extend(["quad_dx0", "quad_dy0"])
        if self.fit_quad_roll:
            params_order.append("quad_roll")

        low_bounds = np.array([bounds[p][0] for p in params_order], dtype=float)
        high_bounds = np.array([bounds[p][1] for p in params_order], dtype=float)

        def predict_from_params(params, allow_stop=True):
            emit_x_norm = float(params["emit_x_norm"])
            beta_x0 = float(params["beta_x0"])
            alpha_x0 = float(params["alpha_x0"])
            emit_y_norm = float(params["emit_y_norm"])
            beta_y0 = float(params["beta_y0"])
            alpha_y0 = float(params["alpha_y0"])

            if self.fit_quadrupole_strength:
                K1L_values_used = float(params["quad_k1l_0"]) * (1.0 + deltas_for_fit)
            else:
                K1L_values_used = K1L_values

            stop_checker = (lambda: self._stop_requested or self._pause_requested) if allow_stop else None
            try:
                if self.fit_quad_offset or self.fit_quad_roll:
                    full = self.interface.predict_emittance_scan_response_full(
                        quad_name=quad_name, screens=screens, K1L_values=K1L_values_used,
                        emit_x=emit_x_norm, emit_y=emit_y_norm, beta_x0=beta_x0, beta_y0=beta_y0,
                        alpha_x0=alpha_x0, alpha_y0=alpha_y0,
                        quad_dx0=(float(params["quad_dx0"]) if self.fit_quad_offset else None),
                        quad_dy0=(float(params["quad_dy0"]) if self.fit_quad_offset else None),
                        quad_roll=(float(params["quad_roll"]) if self.fit_quad_roll else None),
                        reference_screen=screens[0], stop_checker=stop_checker)
                    pred = {k: np.asarray(v, dtype=float) for k, v in full.items()}
                else:
                    pred_sigx, pred_sigy = self.interface.predict_emittance_scan_response(
                        quad_name=quad_name, screens=screens, K1L_values=K1L_values_used,
                        emit_x=emit_x_norm, emit_y=emit_y_norm, beta_x0=beta_x0, beta_y0=beta_y0,
                        alpha_x0=alpha_x0, alpha_y0=alpha_y0,
                        reference_screen=screens[0], stop_checker=stop_checker)
                    pred = {"sigma_x": np.asarray(pred_sigx, dtype=float), "sigma_y": np.asarray(pred_sigy, dtype=float)}
            except RuntimeError as e:
                if str(e) == "__OPTIMIZATION_STOP__":
                    if self._pause_requested:
                        raise OptimizationPaused("Optimization paused.")
                    raise OptimizationStopped("Optimization stopped.")
                raise
            return pred

        def _residual_blocks(pred):
            blocks = [
                (pred["sigma_x"] - sig_x)[valid_x] / u_x[valid_x],
                (pred["sigma_y"] - sig_y)[valid_y] / u_y[valid_y],
            ]
            if self.fit_quad_offset:
                blocks.append((pred["x_mean"] - dx_meas)[valid_dx] / u_dx[valid_dx])
                blocks.append((pred["y_mean"] - dy_meas)[valid_dy] / u_dy[valid_dy])
            if self.fit_quad_roll:
                blocks.append((pred["sigma_xy"] - sigxy_meas)[valid_sigxy] / u_sigxy[valid_sigxy])
            return blocks

        def compute_cost(params, allow_stop=True):
            if allow_stop and (self._stop_requested or self._pause_requested):
                if self._pause_requested:
                    raise OptimizationPaused("Optimization paused.")
                raise OptimizationStopped("Optimization stopped.")
            pred = predict_from_params(params, allow_stop=allow_stop)

            if np.any(valid_x) and not np.all(np.isfinite(pred["sigma_x"][valid_x])): return 1e12, pred
            if np.any(valid_y) and not np.all(np.isfinite(pred["sigma_y"][valid_y])): return 1e12, pred

            residuals = np.concatenate([block.ravel() for block in _residual_blocks(pred)])
            cost = float(np.sum(residuals ** 2))
            return cost, pred

        original_low_bounds = low_bounds.copy()
        original_high_bounds = high_bounds.copy()
        best_row = None
        best_cost = np.inf
        stopped_during_fit = False

        try:
            linear_optics = estimate_twiss_use_linear_optics_start(self.interface, quad_name, screens, K1L_values, sig_x, sig_y, u_x, u_y, valid_x, valid_y, beta_gamma)
        except CheckLinearOpticsUnavailable as e:
            raise RuntimeError(f"Fit will not converge: {e} Try a different K1L scan range/number of steps, or scan a different quadrupole.") from e
        except Exception as e:
            raise RuntimeError(f"Fit will not converge: could not compute a linear-optics starting estimate. {e} Try a different K1L scan range/number of steps, or scan a different quadrupole.") from e

        x0_linear_optics_values = [linear_optics[p] for p in ("emit_x_norm", "beta_x0", "alpha_x0", "emit_y_norm", "beta_y0", "alpha_y0")]
        if self.fit_quadrupole_strength: x0_linear_optics_values.append(K1L_0_readback)
        if self.fit_quad_offset: x0_linear_optics_values.extend([0.0, 0.0])  # start from "no offset", the fit pulls away from it if the data supports it
        if self.fit_quad_roll: x0_linear_optics_values.append(0.0)  # start from "no roll"
        x0_linear_optics = np.clip(np.array(x0_linear_optics_values, dtype=float), original_low_bounds, original_high_bounds)
        cost_linear_optics, _ = compute_cost(dict(zip(params_order, x0_linear_optics)), allow_stop=False)

        if not np.isfinite(cost_linear_optics):
            raise RuntimeError("Fit will not converge: the linear-optics starting estimate does not reproduce the measured beam sizes. Try a different K1L scan range/number of steps, or scan a different quadrupole.")

        row_values = dict(zip(params_order, x0_linear_optics))
        row_values["f"] = float(cost_linear_optics)
        best_row = pd.Series(row_values)
        best_cost = float(cost_linear_optics)
        rel_window = 0.15
        half_width = rel_window * (original_high_bounds - original_low_bounds)
        half_width[n_core_params:] = (original_high_bounds - original_low_bounds)[n_core_params:]
        low_bounds = np.maximum(original_low_bounds, x0_linear_optics - half_width)
        high_bounds = np.minimum(original_high_bounds, x0_linear_optics + half_width)
        low_bounds = np.minimum(low_bounds, x0_linear_optics - 1e-9)
        high_bounds = np.maximum(high_bounds, x0_linear_optics + 1e-9)
        print(f"Linear optics cost: cost={cost_linear_optics:.6g}, x0={row_values}")

        if self._stop_requested or self._pause_requested:
            pred_partial = predict_from_params(best_row[params_order].to_dict(), allow_stop=False)
            solution = self._build_joint_partial_output(screens=screens, sigma_x=sig_x, sigma_y=sig_y, pred_x=pred_partial["sigma_x"], pred_y=pred_partial["sigma_y"], best_row=best_row, best_cost=best_cost)
            if self._pause_requested:
                raise OptimizationPaused("Optimization paused.", solution=solution)
            return solution
        pred = predict_from_params(best_row[params_order].to_dict(), allow_stop=True)
        print(f"Starting local optimization from f={best_cost:.4g}...")
        x0 = np.array([float(best_row[p]) for p in params_order], dtype=float)
        x0 = np.clip(x0, low_bounds, high_bounds)

        def _move_away_from_bounds_edges(point):
            point = np.asarray(point, dtype=float)
            margin = 0.03 * (high_bounds - low_bounds)
            return np.clip(point, low_bounds + margin, high_bounds - margin)

        x0_try = _move_away_from_bounds_edges(x0)
        ls_best_cost = [float(best_cost)]
        ls_best_params = [x0.copy()]
        ls_stopped = [False]
        ls_eval = [0]

        n_x = max(int(np.count_nonzero(valid_x)), 1)
        n_y = max(int(np.count_nonzero(valid_y)), 1)
        n_dx = max(int(np.count_nonzero(valid_dx)), 1) if self.fit_quad_offset else 0
        n_dy = max(int(np.count_nonzero(valid_dy)), 1) if self.fit_quad_offset else 0
        n_sigxy = max(int(np.count_nonzero(valid_sigxy)), 1) if self.fit_quad_roll else 0
        n_residuals_total = n_x + n_y + n_dx + n_dy + n_sigxy

        def _raw_residuals(p_c):
            params = dict(zip(params_order, np.asarray(p_c, dtype=float)))
            pred = predict_from_params(params, allow_stop=False)
            return np.concatenate([block.ravel() for block in _residual_blocks(pred)])

        def _ls_residuals(z):
            if self._stop_requested or self._pause_requested:
                ls_stopped[0] = True
                raise StopIteration("Local least-squares stop requested")

            p_c = np.asarray(z, dtype=float)
            try:
                residuals = _raw_residuals(p_c)
            except Exception:
                return np.full(n_residuals_total, 1e3, dtype=float)

            f = float(np.sum(residuals ** 2))
            if np.isfinite(f) and f < ls_best_cost[0]:
                ls_best_cost[0] = f
                ls_best_params[0] = p_c.copy()

            ls_eval[0] += 1
            self._emit_progress("Least squares", min(ls_eval[0], 200), 200)
            params = dict(zip(params_order, p_c))
            print(f"Least squares {ls_eval[0]}: best_f={ls_best_cost[0]:.4g}, " + ", ".join(f"{k}={v:.6g}" for k, v in params.items()))

            return residuals

        stagnation_patience = 25
        min_improvement_of_cost = 1e-3
        good_fit_final_cost = 5e-11
        best_cost_in_this_fit = [np.inf]
        steps_without_improvement = [0]
        reason_to_stop = [None]

        def exit_ls_if_no_improvement_or_reached_goal(intermediate_result):
            current_cost = 2.0 * float(intermediate_result.cost)
            if not np.isfinite(current_cost):
                return
            if current_cost <= good_fit_final_cost:
                reason_to_stop[0] = "target cost reached"
                raise StopIteration
            if not np.isfinite(best_cost_in_this_fit[0]):
                best_cost_in_this_fit[0] = current_cost
                return
            if current_cost < best_cost_in_this_fit[0]:
                relative_improvement = ((best_cost_in_this_fit[0] - current_cost) / max(abs(best_cost_in_this_fit[0]), 1e-12))
                best_cost_in_this_fit[0] = current_cost
                if relative_improvement >= min_improvement_of_cost:
                    steps_without_improvement[0] = 0
                else:
                    steps_without_improvement[0] += 1
            else:
                steps_without_improvement[0] += 1

            if steps_without_improvement[0] >= stagnation_patience:
                reason_to_stop[0] = "no meaningful improvement"
                raise StopIteration

        try:
            res_try = least_squares(_ls_residuals, x0_try, bounds=(low_bounds, high_bounds), method="trf", loss="linear", f_scale=1.0, max_nfev=200, x_scale=np.maximum(high_bounds - low_bounds, 1e-12), ftol=1e-8, xtol=1e-8, gtol=1e-8, callback = exit_ls_if_no_improvement_or_reached_goal)
            p_try = np.asarray(res_try.x, dtype=float)
            f_try, _ = compute_cost(dict(zip(params_order, p_try)), allow_stop=False)
            if np.isfinite(f_try) and f_try < ls_best_cost[0]:
                ls_best_cost[0] = float(f_try)
                ls_best_params[0] = p_try.copy()
            print(
                f"Least squares fit: cost={float(f_try):.4g}")
            if reason_to_stop[0] is not None:
                print(f"Stopping LS: {reason_to_stop[0]}.")

        except StopIteration:
            ls_stopped[0] = True
            print("  LS interrupted.")
        except Exception as e:
            print(f"  LS failed ({e}).")

        p_final = ls_best_params[0]
        best_cost_final = ls_best_cost[0]

        best_row = best_row.copy()
        for i, name in enumerate(params_order):
            best_row[name] = float(p_final[i])

        stopped_during_fit = stopped_during_fit or ls_stopped[0]

        if stopped_during_fit or self._stop_requested or self._pause_requested:
            pred_p = predict_from_params(best_row[params_order].to_dict(), allow_stop=False)
            solution = self._build_joint_partial_output(screens=screens, sigma_x=sig_x, sigma_y=sig_y, pred_x=pred_p["sigma_x"], pred_y=pred_p["sigma_y"], best_row=best_row, best_cost=best_cost_final)
            if self._pause_requested:
                raise OptimizationPaused("Optimization paused.", solution=solution)
            return solution

        pred_final = predict_from_params(best_row[params_order].to_dict(), allow_stop=True)

        solution = self._build_joint_partial_output(screens=screens, sigma_x=sig_x, sigma_y=sig_y, pred_x=pred_final["sigma_x"], pred_y=pred_final["sigma_y"], best_row=best_row, best_cost=best_cost_final)
        try:
            param_scale = np.maximum(high_bounds - low_bounds, 1e-12)  # same scale least_squares used as x_scale
            J_at_p_final = _finite_diff_jacobian(_raw_residuals, p_final, low_bounds, high_bounds, param_scale=param_scale)
            fit_error = self._calculate_optimalization_errors(J_at_p_final, n_params=len(params_order), param_scale=param_scale)
        except Exception as e:
            print(f"Couldn't calculate optimalization error: {e}.")
            fit_error = {"param_errors": np.full(len(params_order), np.nan), "cov": None}
        param_errors = fit_error["param_errors"]
        if param_errors is None or len(param_errors) != len(params_order):
            err_dict = {p: np.nan for p in params_order}
        else:
            err_dict = {p: float(e) for p, e in zip(params_order, param_errors)}

        solution["param_errors"] = err_dict
        solution["param_cov"] = fit_error["cov"]

        print(
            f"Fit parameter errors: "
            + ", ".join(f"{k}={v:.4g}" for k, v in err_dict.items())
        )

        return solution
