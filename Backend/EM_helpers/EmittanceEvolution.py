try:
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import QVBoxLayout, QDialog, QLabel, QSizePolicy
except ImportError:
    from PyQt5.QtCore import Qt
    from PyQt5.QtWidgets import QVBoxLayout, QDialog, QLabel, QSizePolicy
import numpy as np
import matplotlib
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

FIT_PARAMS = ("emit_x_norm", "emit_y_norm", "beta_x0", "beta_y0", "alpha_x0", "alpha_y0")


class EmittanceEvolution(QDialog):
    def __init__(self, interface, parent=None, screens=None):
        super().__init__(parent)
        self.setWindowTitle("Emittance evolution")
        self.setWindowFlag(Qt.WindowType.WindowContextHelpButtonHint, False)
        self.setMinimumSize(800, 500)
        self.resize(1000, 650)
        self.setSizeGripEnabled(True)
        self.interface = interface
        self.screens = list(screens or [])
        self.figure = Figure(figsize=(10, 6), constrained_layout=True)
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.header = QLabel("", self)
        self.header.setWordWrap(True)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(6)
        layout.addWidget(self.header, 0)
        layout.addWidget(self.canvas, 1)

    def _model_interface(self):
        return getattr(self.interface, "tracking_interface", self.interface)

    def _quadrupole_value(self, session, result):
        for value in (result.get("quad_k1l_0"), session.get("K1L_0")):
            value = np.asarray(value, dtype=float).ravel()
            if value.size and np.isfinite(value[0]):
                return float(value[0])
        scan_values = np.asarray(session.get("K1L_values", []), dtype=float)
        if scan_values.size:
            return float(scan_values[scan_values.size // 2])
        return np.nan

    def _track_fitted_beam(self, names, session, result):
        quad_name = session.get("quad_name")
        quad_value = self._quadrupole_value(session, result)
        energy_pref = float(result.get("energy_pref", np.nan))
        quad_dx0 = float(result.get("quad_dx0", np.nan))
        quad_dy0 = float(result.get("quad_dy0", np.nan))
        quad_roll = float(result.get("quad_roll", np.nan))
        offsets_fitted = bool(result.get("fit_quad_offset")) and np.isfinite(quad_dx0) and np.isfinite(quad_dy0)

        full = self._model_interface().predict_emittance_scan_response_full(
            quad_name=quad_name, screens=names, K1L_values=np.array([quad_value], dtype=float),
            emit_x=float(result["emit_x_norm"]), emit_y=float(result["emit_y_norm"]),
            beta_x0=float(result["beta_x0"]), beta_y0=float(result["beta_y0"]),
            alpha_x0=float(result["alpha_x0"]), alpha_y0=float(result["alpha_y0"]),
            quad_dx0=(quad_dx0 * 1e-3 if offsets_fitted else None),
            quad_dy0=(quad_dy0 * 1e-3 if offsets_fitted else None),
            quad_roll=(quad_roll * 1e-3 if bool(result.get("fit_quad_roll")) and np.isfinite(quad_roll) else None),
            energy_pref=(energy_pref if bool(result.get("fit_energy_pref")) and np.isfinite(energy_pref) else None),
            with_twiss=True, reference_screen=names[0])
        return np.asarray(full["emitt_x"][0], dtype=float), np.asarray(full["emitt_y"][0], dtype=float), quad_value

    def _plot_plane(self, axes, s, fitted_curve, default_curve, fitted, fitted_error, names, label, colour):
        if fitted_curve is not None:
            axes.plot(s, fitted_curve, "o-", color=colour, label=f"fitted beam tracked through the model, {label}")
            for position, value, name in zip(s, fitted_curve, names):
                if np.isfinite(value):
                    axes.annotate(f"{name}\n{value:.4g}", (position, value), textcoords="offset points",
                                  xytext=(0, 8), ha="center", fontsize=7, color=colour)
        axes.plot(s, default_curve, "s--", color="0.6", markersize=4, linewidth=1.0, label=f"model default beam, {label}")
        if fitted_curve is None:
            for position, value, name in zip(s, default_curve, names):
                if np.isfinite(value):
                    axes.annotate(f"{name}\n{value:.4g}", (position, value), textcoords="offset points",
                                  xytext=(0, 8), ha="center", fontsize=7, color="0.4")
        if np.isfinite(fitted):
            axes.axhline(fitted, color="k", linestyle=":", linewidth=1.0, label = None)
            if np.isfinite(fitted_error) and fitted_error > 0:
                axes.axhspan(fitted - fitted_error, fitted + fitted_error, color="k", alpha=0.08)
        axes.set_ylabel(f"{label} [mm·mrad]")
        axes.margins(x=0.08, y=0.22)
        axes.grid(True, alpha=0.3)
        axes.legend(fontsize=7, loc="best")

    def _display_emittance_evolution(self, screens=None, session=None, result=None):
        if screens is not None:
            self.screens = list(screens)
        self.figure.clear()
        axes_x = self.figure.add_subplot(2, 1, 1)
        axes_y = self.figure.add_subplot(2, 1, 2, sharex=axes_x)
        session = session if isinstance(session, dict) else {}
        result = result if isinstance(result, dict) else {}

        try:
            default = self._model_interface().get_emittance_at_screens(names=self.screens or None)
        except Exception as e:
            self.canvas.draw_idle()
            return

        names = list(default.get("names", []))
        if not names:
            self.canvas.draw_idle()
            return

        s = np.asarray(default["S"], dtype=float)
        default_x = np.asarray(default["emitt_x"], dtype=float)
        default_y = np.asarray(default["emitt_y"], dtype=float)

        fit_available = all(np.isfinite(float(result.get(parameter, np.nan))) for parameter in FIT_PARAMS)
        fitted_x = fitted_y = None
        quad_value = np.nan
        note = ""
        if fit_available:
            fitted_x, fitted_y, quad_value = self._track_fitted_beam(names, session, result)

        self._plot_plane(axes_x, s, fitted_x, default_x, float(result.get("emit_x_norm", np.nan)), float(result.get("emit_x_norm_err", np.nan)), names, "εₓ (norm.)", "tab:blue")
        self._plot_plane(axes_y, s, fitted_y, default_y, float(result.get("emit_y_norm", np.nan)), float(result.get("emit_y_norm_err", np.nan)), names, "εᵧ (norm.)", "tab:red")
        axes_y.set_xlabel("S [m]")
        axes_x.set_title("Normalised emittance along the selected screens")

        curve = fitted_x if fitted_x is not None else default_x
        growth = curve[-1] / curve[0] - 1.0 if curve.size > 1 and curve[0] > 0 else np.nan
        curve_y = fitted_y if fitted_y is not None else default_y
        growth_y = curve_y[-1] / curve_y[0] - 1.0 if curve_y.size > 1 and curve_y[0] > 0 else np.nan
        quadrupole = f"{session.get('quad_name')} at {quad_value:.4g}" if np.isfinite(quad_value) else "-"
        self.canvas.draw_idle()
