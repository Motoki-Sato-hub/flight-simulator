import os
import numpy as np

'''
LinearResponse - pure mathematics behind linear response approach
linear_response_engine.py - using this mathematics on session from GUI
sigma² = R11² * <x²> + 2 R11 R12 * <x x'> + R12² * <x'²>
<x²>   = emit * beta
<x x'> = -emit * alpha
<x'²>  = emit * gamma
gamma = (1 + alpha²) / beta
'''
class LinearResponse:
    def __init__(self, coefficients_path, dataset_path, beta_gamma):
        self.coefficients_path = coefficients_path
        self.dataset_path = dataset_path
        self.beta_gamma = beta_gamma
        if os.path.isfile(self.coefficients_path):
            self._load_coefficients(self.coefficients_path)
        elif os.path.isfile(self.dataset_path):
            self.fit_coefficients_from_R_dataset(self.dataset_path)
            self.save_coefficients(self.coefficients_path)
        else:
            raise FileNotFoundError(f"Linear response data not found. Expected either {self.coefficients_path} or {self.dataset_path}.")

    def _load_coefficients(self, path):
        data = np.load(path, allow_pickle=True)
        self.screens = [str(s) for s in data["screens"]]
        self.Rx = np.asarray(data["Rx_fit"], dtype=float)
        self.Ry = np.asarray(data["Ry_fit"], dtype=float)
        if "beta_gamma" in data.files:
            self.beta_gamma = float(data["beta_gamma"])

    def save_coefficients(self, path=None):
        if path is None:
            path = self.coefficients_path
        os.makedirs(os.path.dirname(path), exist_ok=True)
        np.savez(path, screens=np.asarray(self.screens), Rx_fit=self.Rx, Ry_fit=self.Ry, beta_gamma=np.array(self.beta_gamma, dtype=float))

    def fit_coefficients_from_R_dataset(self, dataset_path=None):
        if dataset_path is None:
            dataset_path = self.dataset_path
        dataset = np.load(dataset_path, allow_pickle=True)
        if "beta_gamma" in dataset.files:
            self.beta_gamma = float(dataset["beta_gamma"])
        elif self.beta_gamma is not None:
            self.beta_gamma = float(self.beta_gamma)
            print(f"Warning: dataset {dataset_path} does not contain beta_gamma, using fallback from interface.")
        else:
            raise RuntimeError(
                f"Dataset {dataset_path} does not contain beta_gamma "
                "and no fallback beta_gamma was provided."
            )
        if not np.isfinite(self.beta_gamma) or self.beta_gamma <= 0:
            raise RuntimeError(f"Invalid beta_gamma: {self.beta_gamma}")

        I = np.asarray(dataset["X"], dtype=float) # twiss parameters at the entrance
        # I[i] = [
        #     emit_x_norm,
        #     beta_x,
        #     alpha_x,
        #     emit_y_norm,
        #     beta_y,
        #     alpha_y
        # ]

        # O[i] = [
        #     [
        #         sigx_screen0,
        #         sigy_screen0,
        #         sigx_screen1,
        #         sigy_screen1,
        #         sigx_screen2,
        #         sigy_screen2,
        #         sigx_screen3,
        #         sigy_screen3,
        #     ]
        # ]

        O = np.asarray(dataset["Y"], dtype=float) # response on screens, sigma**2
        self.screens = [str(s) for s in dataset["screens"]]

        emit_x_geom = I[:, 0] / self.beta_gamma
        beta_x = I[:, 1]
        alpha_x = I[:, 2]
        gamma_x = (1.0 + alpha_x**2) / beta_x

        emit_y_geom = I[:, 3] / self.beta_gamma
        beta_y = I[:, 4]
        alpha_y = I[:, 5]
        gamma_y = (1.0 + alpha_y**2) / beta_y

        ones = np.ones(len(I))
        Cx = np.column_stack((beta_x * emit_x_geom, -alpha_x * emit_x_geom, gamma_x * emit_x_geom, ones))
        Cy = np.column_stack((beta_y * emit_y_geom, -alpha_y * emit_y_geom, gamma_y * emit_y_geom, ones))

        Bx = O[:, 0::2]**2
        By = O[:, 1::2]**2

        self.Rx = np.linalg.lstsq(Cx, Bx, rcond=None)[0].T
        self.Ry = np.linalg.lstsq(Cy, By, rcond=None)[0].T

        return self.Rx, self.Ry

    def predict_sigma2_from_twiss_set(self, screens, emit_x_norm, beta_x0, alpha_x0, emit_y_norm, beta_y0, alpha_y0):

        emit_x_geom = float(emit_x_norm) / float(self.beta_gamma)
        emit_y_geom = float(emit_y_norm) / float(self.beta_gamma)

        beta_x0 = float(beta_x0)
        alpha_x0 = float(alpha_x0)
        beta_y0 = float(beta_y0)
        alpha_y0 = float(alpha_y0)

        gamma_x = (1.0 + alpha_x0**2) / beta_x0
        gamma_y = (1.0 + alpha_y0**2) / beta_y0

        Cx = np.array([beta_x0 * emit_x_geom, -alpha_x0 * emit_x_geom, gamma_x * emit_x_geom, 1.0])
        Cy = np.array([beta_y0 * emit_y_geom, -alpha_y0 * emit_y_geom, gamma_y * emit_y_geom, 1.0])

        idx = [self.screens.index(str(screen)) for screen in screens]

        pred_x = Cx @ self.Rx[idx].T
        pred_y = Cy @ self.Ry[idx].T

        return pred_x.reshape(1, -1), pred_y.reshape(1, -1)

    def solve_twiss_from_measured_sigma2(self, screens, sigma2_x, sigma2_y):
        idx = [self.screens.index(str(screen)) for screen in screens]

        if len(idx) < 3:
            raise RuntimeError("At least 3 screens are required for direct linear R-response fit.")

        Rx = self.Rx[idx]
        Ry = self.Ry[idx]

        Mx = Rx[:, :3] # R11^2 2R111R12 R12^2
        My = Ry[:, :3]

        yx = np.asarray(sigma2_x, dtype=float).reshape(-1) - Rx[:, 3]
        yy = np.asarray(sigma2_y, dtype=float).reshape(-1) - Ry[:, 3]

        px = np.linalg.lstsq(Mx, yx, rcond=None)[0]
        py = np.linalg.lstsq(My, yy, rcond=None)[0]

        emit_x_geom = np.sqrt(px[0] * px[2] - px[1] ** 2)
        beta_x0 = px[0] / emit_x_geom
        alpha_x0 = -px[1] / emit_x_geom

        emit_y_geom = np.sqrt(py[0] * py[2] - py[1] ** 2)
        beta_y0 = py[0] / emit_y_geom
        alpha_y0 = -py[1] / emit_y_geom

        return {
            "emit_x_geom": emit_x_geom,
            "emit_y_geom": emit_y_geom,
            "emit_x_norm": emit_x_geom * self.beta_gamma,
            "emit_y_norm": emit_y_geom * self.beta_gamma,
            "beta_x0": beta_x0,
            "alpha_x0": alpha_x0,
            "beta_y0": beta_y0,
            "alpha_y0": alpha_y0,
            "pred_x": (Mx @ px + Rx[:, 3]).reshape(1, -1),
            "pred_y": (My @ py + Ry[:, 3]).reshape(1, -1),
        }