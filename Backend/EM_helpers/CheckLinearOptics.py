import numpy as np
from scipy.optimize import least_squares

class CheckLinearOpticsUnavailable(Exception):
    pass

def _solve_plane(R1, R2, sigma, u, valid):
    R1 = np.asarray(R1, dtype=float)[valid]
    R2 = np.asarray(R2, dtype=float)[valid]
    sigma_v = np.asarray(sigma, dtype=float)[valid] # beam size
    sigma2 = sigma_v ** 2
    u_v = np.asarray(u, dtype=float)[valid] # error of measurement

    if R1.size < 3: raise CheckLinearOpticsUnavailable("Not enough valid scan points for a first estimate (need >= 3).")

    u_sigma2 = np.sqrt((2.0 * np.abs(sigma_v) * u_v) ** 2)
    u_sigma2 = np.where(np.isfinite(u_sigma2) & (u_sigma2 > 0), u_sigma2, np.maximum(np.abs(sigma2), 1e-12))

    coefficients = np.column_stack([R1 ** 2,
                                    2.0 * R1 * R2,
                                    R2 ** 2])
    twiss_parameters_joint, _, rank, _ = np.linalg.lstsq(coefficients / u_sigma2[:, None], sigma2 / u_sigma2, rcond=None)
    if rank < 3: raise CheckLinearOpticsUnavailable("K1L scan is too small or too few distinct points. Please do a wider scan, or with more steps.")

    eps_beta, eps_alpha, eps_gamma = twiss_parameters_joint
    # beta*gamma - alpha^2 = 1 ---- * eps^2
    # eps^2 * beta * gamma - alpha^2*eps^2 = eps^2 (geom)
    eps_squared = eps_beta * eps_gamma - eps_alpha ** 2
    unconstrained_is_physical = bool(np.isfinite(eps_squared) and eps_squared > 0 and eps_beta > 0 and eps_gamma > 0)

    # sigma^2 = (a*R1)^2 + 2*a*b*R1*R2 + (b^2+c^2)*R2^2 stays a physical beam matrix for any a > 0, c > 0
    def residuals(parameters):
        a, b, c = parameters
        model = (a * R1) ** 2 + 2.0 * a * b * R1 * R2 + (b ** 2 + c ** 2) * R2 ** 2
        return (model - sigma2) / u_sigma2

    if unconstrained_is_physical:
        a_start = np.sqrt(eps_beta)
        start = [a_start, eps_alpha / a_start, np.sqrt(max(eps_gamma - eps_alpha ** 2 / eps_beta, 1e-12))]
    else:
        sigma_typical = float(np.median(np.abs(sigma_v)))
        start = [sigma_typical, 0.0, sigma_typical / max(float(np.median(np.abs(R2))), 1e-9)]

    solution = least_squares(residuals, start, bounds=([1e-12, -np.inf, 1e-12], [np.inf, np.inf, np.inf]))
    a, b, c = solution.x
    emit_geom = float(abs(a * c))
    beta0 = float(a / c)
    alpha0 = float(-b / c)
    reduced_chi2 = float(np.sum(solution.fun ** 2) / max(R1.size - 3, 1))

    if not np.isfinite(emit_geom) or emit_geom <= 0 or not np.isfinite(beta0) or beta0 <= 0 or not np.isfinite(alpha0):
        raise CheckLinearOpticsUnavailable(
            f"Estimate is unphysical (eps_beta={eps_beta:.4g}, eps_gamma={eps_gamma:.4g}, det={eps_squared:.4g}):")
    return emit_geom, beta0, alpha0, reduced_chi2, unconstrained_is_physical

def estimate_twiss_use_linear_optics_start(interface, quad_name, screens, K1L_values, sig_x, sig_y, u_x, u_y, valid_x, valid_y, beta_gamma):

    transport = interface.get_R_matrix_scan(quad_name=quad_name, screens=screens, K1L_values=K1L_values)
    R11, R12 = np.asarray(transport["R11"], dtype=float), np.asarray(transport["R12"], dtype=float)
    R33, R34 = np.asarray(transport["R33"], dtype=float), np.asarray(transport["R34"], dtype=float)

    finite_R_x = np.isfinite(R11) & np.isfinite(R12)
    finite_R_y = np.isfinite(R33) & np.isfinite(R34)

    emit_x_geom, beta_x0, alpha_x0, reduced_chi2_x, physical_x = _solve_plane(R11, R12, sig_x, u_x, valid_x & finite_R_x)
    emit_y_geom, beta_y0, alpha_y0, reduced_chi2_y, physical_y = _solve_plane(R33, R34, sig_y, u_y, valid_y & finite_R_y)

    return {
        "emit_x_norm": emit_x_geom * beta_gamma,
        "beta_x0": beta_x0,
        "alpha_x0": alpha_x0,
        "emit_y_norm": emit_y_geom * beta_gamma,
        "beta_y0": beta_y0,
        "alpha_y0": alpha_y0,
        "reduced_chi2_x": reduced_chi2_x,
        "reduced_chi2_y": reduced_chi2_y,
        "unconstrained_is_physical_x": physical_x,
        "unconstrained_is_physical_y": physical_y,
    }
