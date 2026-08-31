import numpy as np

class CheckLinearOpticsUnavailable(Exception):
    pass

def _solve_plane(R1, R2, sigma, u, valid):
    R1 = np.asarray(R1, dtype=float)[valid]
    R2 = np.asarray(R2, dtype=float)[valid]
    sigma_v = np.asarray(sigma, dtype=float)[valid] # beam size
    sigma2 = sigma_v ** 2
    u_v = np.asarray(u, dtype=float)[valid] # error of measurement

    if R1.size < 3: raise CheckLinearOpticsUnavailable("Not enough valid scan points for a linear estimate (need >= 3).")

    coefficients = np.column_stack([R1 ** 2,
                                    2.0 * R1 * R2,
                                    R2 ** 2])
    relative_uncertainty = u_v / np.maximum(np.abs(sigma_v), 1e-9)
    row_weight = 1.0 / np.maximum(relative_uncertainty, 1e-3) # it checks if the point is more precise
    coefficients_w = coefficients * row_weight[:, None]
    sigma_w = sigma2 * row_weight

    twiss_parameters_joint, residuals, rank, _ = np.linalg.lstsq(coefficients_w, sigma_w, rcond=None)
    if rank < 3: raise CheckLinearOpticsUnavailable("K1L scan is too small or too few distinct points. Please do a wider scan, or with more steps.")

    predicted_w_sigmas2 = coefficients_w @ twiss_parameters_joint
    degrees_of_freedom = max(len(sigma_w) - 3, 1) # because we fit eps*beta, -eps*alpha, eps*gamma
    reduced_chi2 = float(np.sum((sigma_w - predicted_w_sigmas2) ** 2) / degrees_of_freedom) # sigma measured - sigma modeled
    if not np.isfinite(reduced_chi2) or reduced_chi2 > 5.0: # threshold for how well a linear (no sextupole nonlinearity) model must fit the scan
        raise CheckLinearOpticsUnavailable(f"Linear model doesn't fit the scan well (chi_2 = {reduced_chi2:.3g}). There is probably dispersion or nonlinearities (e.g. sextupoles) between the elements: try a narrower K1L range, or scan a different quadrupole.")

    eps_beta, eps_alpha, eps_gamma = twiss_parameters_joint
    # beta*gamma - alpha^2 = 1 ---- * eps^2
    # eps^2 * beta * gamma - alpha^2*eps^2 = eps^2 (geom)
    eps_squared = eps_beta * eps_gamma - eps_alpha ** 2
    if not np.isfinite(eps_squared) or eps_squared <= 0 or eps_beta <= 0 or eps_gamma <= 0: raise CheckLinearOpticsUnavailable(f"Linear estimate is unphysical (eps_beta={eps_beta:.4g}, eps_gamma={eps_gamma:.4g}, det={eps_squared:.4g}): try a different K1L scan range, more steps, or a different quadrupole.")

    emit_geom = float(np.sqrt(eps_squared))
    beta0 = float(eps_beta / emit_geom)
    alpha0 = float(-eps_alpha / emit_geom)
    return emit_geom, beta0, alpha0

def estimate_twiss_use_linear_optics_start(interface, quad_name, screens, K1L_values, sig_x, sig_y, u_x, u_y, valid_x, valid_y, beta_gamma):

    transport = interface.get_R_matrix_scan(quad_name=quad_name, screens=screens, K1L_values=K1L_values)
    R11, R12 = np.asarray(transport["R11"], dtype=float), np.asarray(transport["R12"], dtype=float)
    R33, R34 = np.asarray(transport["R33"], dtype=float), np.asarray(transport["R34"], dtype=float)

    finite_R_x = np.isfinite(R11) & np.isfinite(R12)
    finite_R_y = np.isfinite(R33) & np.isfinite(R34)

    emit_x_geom, beta_x0, alpha_x0 = _solve_plane(R11, R12, sig_x, u_x, valid_x & finite_R_x)
    emit_y_geom, beta_y0, alpha_y0 = _solve_plane(R33, R34, sig_y, u_y, valid_y & finite_R_y)

    return {
        "emit_x_norm": emit_x_geom * beta_gamma,
        "beta_x0": beta_x0,
        "alpha_x0": alpha_x0,
        "emit_y_norm": emit_y_geom * beta_gamma,
        "beta_y0": beta_y0,
        "alpha_y0": alpha_y0,
    }
