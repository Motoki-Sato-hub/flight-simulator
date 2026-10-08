try:
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import QVBoxLayout, QHBoxLayout, QDialog, QLabel, QWidget, QComboBox
except ImportError:
    from PyQt5.QtCore import Qt
    from PyQt5.QtWidgets import QVBoxLayout, QHBoxLayout, QDialog, QLabel, QWidget, QComboBox
import numpy as np
import matplotlib
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.lines import Line2D


## if beam size at for example OTR2X is sigma,
## phase space on OTR0X has to fit between the constrain lines

class PhaseSpaces(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Phase Space Graphs")
        self.setWindowFlag(Qt.WindowType.WindowContextHelpButtonHint, False)
        self.setMinimumSize(800, 500)
        self.resize(1000, 650)
        self.setSizeGripEnabled(True)
        self._plot_data = None
        self.figure = Figure(figsize=(10, 6), constrained_layout=True)
        self.canvas = FigureCanvas(self.figure)
        plot_widget = QWidget(self)
        plot_layout = QVBoxLayout(plot_widget)
        plot_layout.setContentsMargins(0, 0, 0, 0)
        plot_layout.addWidget(self.canvas)
        self.location_combo = QComboBox(self)
        self.location_combo.currentIndexChanged.connect(self._redraw_current_location)
        header = QWidget(self)
        header_layout = QHBoxLayout(header)
        header_layout.setContentsMargins(0, 0, 0, 0)
        header_layout.addStretch(1)
        header_layout.addWidget(QLabel("Location:", self))
        header_layout.addWidget(self.location_combo)
        layout = QVBoxLayout(self)
        layout.addWidget(header)
        layout.addWidget(plot_widget)

    @staticmethod
    def _ellipse_points(beta0, alpha0):
        '''
        gamma * x^2 + 2*alpha*x*x'+ beta x'^2 = emittance
        gamma = (1+alpha^2)/beta
        In normalised coordinates:
        u = x / sqrt(emittance)
        u' = x' / sqrt(emittance)
        And then:
        gamma * u^2 + 2 * alpha * u * u' + beta * u'^2 = 1
        u = sqrt(beta) * cos(theta)
        u' = -(alpha * cos(theta) + sin(theta)) / sqrt(beta)

        If you substitute u and u' above gamma * u^2 + 2 * alpha * u * u' + beta * u'^2 = 1,
        one gets cos(theta)^2 + sin(theta)^2 = 1, so an ellipse.
        '''

        beta0 = float(beta0)
        alpha0 = float(alpha0)

        theta = np.linspace(0.0, 2.0 * np.pi, 500)
        u = np.sqrt(beta0) * np.cos(theta)
        up = -(alpha0 * np.cos(theta) + np.sin(theta)) / np.sqrt(beta0)

        return u, up

    @staticmethod
    def _is_valid(beta0, alpha0):
        try:
            beta0 = float(beta0)
            alpha0 = float(alpha0)
        except Exception:
            return False

        return np.isfinite(beta0) and np.isfinite(alpha0) and beta0 > 0.0

    @staticmethod
    def propagate_twiss(beta0, alpha0, R11, R12, R21, R22):
        try:
            beta0 = float(beta0)
            alpha0 = float(alpha0)
            R11 = float(R11)
            R12 = float(R12)
            R21 = float(R21)
            R22 = float(R22)
        except Exception:
            return np.nan, np.nan
        gamma0 = (1.0 + alpha0 ** 2) / beta0
        beta_i = R11 ** 2 * beta0 - 2.0 * R11 * R12 * alpha0 + R12 ** 2 * gamma0
        alpha_i = -R11 * R21 * beta0 + (R11 * R22 + R12 * R21) * alpha0 - R12 * R22 * gamma0
        return float(beta_i), float(alpha_i)

    def plot_from_result(self, result, reference_name=None, screen_parameters = None):
        self.figure.clear()

        if not isinstance(result, dict):
            ax = self.figure.add_subplot(111)
            ax.text(0.5, 0.5, "No optimization result provided.", ha="center", va="center", transform=ax.transAxes)
            ax.axis("off")
            self.canvas.draw_idle()
            return

        screens_twiss = []
        reference_label = reference_name or "fitted scan reference"
        screens_twiss.append((
            reference_label,
            {
                "beta_x": result.get("beta_x0", np.nan),
                "beta_y": result.get("beta_y0", np.nan),
                "alpha_x": result.get("alpha_x0", np.nan),
                "alpha_y": result.get("alpha_y0", np.nan),
                "emit_x_norm": result.get("emit_x_norm", np.nan),
                "emit_y_norm": result.get("emit_y_norm", np.nan),
            }))

        if isinstance(screen_parameters, dict):
            screens_dict = screen_parameters.get("screens", screen_parameters)
            if isinstance(screens_dict, dict):
                for screen_name, values in screens_dict.items():
                    if isinstance(values, dict):
                        screens_twiss.append((str(screen_name), values))

        self._plot_data = {
            "result": result,
            "screens_twiss": screens_twiss,
        }

        self.location_combo.blockSignals(True)
        self.location_combo.clear()
        for name, _ in screens_twiss:
            self.location_combo.addItem(str(name))
        self.location_combo.blockSignals(False)
        self._redraw_current_location()

    def _redraw_current_location(self):
        self.figure.clear()
        if not isinstance(self._plot_data, dict):
            self.canvas.draw_idle()
            return
        locations = self._plot_data.get("screens_twiss", [])
        if not locations:
            ax = self.figure.add_subplot(111)
            ax.text(0.5, 0.5, "No data to build phase-space ellipse available.", ha="center", va="center", transform=ax.transAxes)
            ax.axis("off")
            self.canvas.draw_idle()
            return

        idx = int(self.location_combo.currentIndex())
        if idx < 0 or idx >= len(locations):
            idx = 0
        location_name, values = locations[idx]

        beta_x = values.get("beta_x", values.get("beta_x0", np.nan))
        alpha_x = values.get("alpha_x", values.get("alpha_x0", np.nan))
        beta_y = values.get("beta_y", values.get("beta_y0", np.nan))
        alpha_y = values.get("alpha_y", values.get("alpha_y0", np.nan))
        emit_x_norm = values.get("emit_x_norm", np.nan)
        emit_y_norm = values.get("emit_y_norm", np.nan)
        ax_x = self.figure.add_subplot(121)
        ax_y = self.figure.add_subplot(122)

        def draw_one(ax, beta0, alpha0, emit_norm, plane_label):
            ax.axhline(0.0, color="0.5", linestyle=":", linewidth=1.0)
            ax.axvline(0.0, color="0.5", linestyle=":", linewidth=1.0)
            ax.grid(True, alpha=0.3)
            ax.set_aspect("equal", adjustable="box")
            ax.set_xlabel("Position")
            ax.set_ylabel("Angle")
            ax.set_title(plane_label)

            try:
                ax.set_box_aspect(1.0)
            except Exception:
                pass
            ax.set_aspect("equal", adjustable="datalim")
            if not self._is_valid(beta0, alpha0):
                ax.text(0.5, 0.5, "Invalid fitted Twiss", ha="center", va="center", transform=ax.transAxes)
                return

            beta0 = float(beta0)
            alpha0 = float(alpha0)
            u, up = self._ellipse_points(beta0, alpha0)

            theta_ref = np.linspace(0.0, 2.0 * np.pi, 500)
            ax.plot(np.cos(theta_ref), np.sin(theta_ref), color="black", linestyle="--", linewidth=1.2, alpha=0.7)
            ax.plot(u, up, color="black", linewidth=2.0)

            reference_name = str(self._plot_data.get("projection_reference", ""))
            show_projection_lines = str(location_name) == reference_name
            projection_lines = self._plot_data.get("projection_lines") if show_projection_lines else None
            plane_key = "x" if "Horizontal" in str(plane_label) else "y"
            legend_handles = []

            lines_for_plane = []
            if isinstance(projection_lines, dict):
                lines_for_plane = list(projection_lines.get(plane_key, []))

            x_values_for_limits = [np.asarray(u, dtype=float), np.cos(theta_ref)]
            y_values_for_limits = [np.asarray(up, dtype=float), np.sin(theta_ref)]

            for line in lines_for_plane:
                try:
                    R11_l = float(line.get("R11", np.nan))
                    R12_l = float(line.get("R12", np.nan))
                    c_l = float(line.get("sigma_over_sqrt_emit", np.nan))
                except Exception:
                    continue
                if not (np.isfinite(R11_l) and np.isfinite(R12_l) and np.isfinite(c_l)):
                    continue
                if abs(R12_l) < 1e-14:
                    if abs(R11_l) < 1e-14:
                        continue
                    x_const = c_l / R11_l
                    x_values_for_limits.append(np.asarray([+x_const, -x_const], dtype=float))
                    y_values_for_limits.append(np.asarray([-1.5, 1.5], dtype=float))
                else:
                    x_probe = np.linspace(np.nanmin(u), np.nanmax(u), 200)
                    y_plus = (+c_l - R11_l * x_probe) / R12_l
                    y_minus = (-c_l - R11_l * x_probe) / R12_l
                    x_values_for_limits.append(x_probe)
                    x_values_for_limits.append(x_probe)
                    y_values_for_limits.append(y_plus)
                    y_values_for_limits.append(y_minus)

            all_x = np.concatenate([arr[np.isfinite(arr)] for arr in x_values_for_limits if np.asarray(arr).size > 0])
            all_y = np.concatenate([arr[np.isfinite(arr)] for arr in y_values_for_limits if np.asarray(arr).size > 0])
            if all_x.size > 0 and all_y.size > 0:
                x_abs = max(float(np.nanmax(np.abs(all_x))), 1.0)
                y_abs = max(float(np.nanmax(np.abs(all_y))), 1.0)
                x_abs *= 1.10
                y_abs *= 1.10
                ax.set_xlim(-x_abs, x_abs)
                ax.set_ylim(-y_abs, y_abs)
            else:
                pad = 1.25
                ax.set_xlim(np.nanmin(u) * pad, np.nanmax(u) * pad)
                ax.set_ylim(np.nanmin(up) * pad, np.nanmax(up) * pad)

            if isinstance(projection_lines, dict):
                color_cycle = ["blue", "lime", "red", "magenta", "orange", "cyan", "purple", "brown"]
                for line_index, line in enumerate(lines_for_plane):
                    screen_label = str(line.get("screen", ""))
                    screen_colors = self._plot_data.get("screen_colors", {})
                    color = screen_colors.get(screen_label, color_cycle[line_index % len(color_cycle)])
                    self._draw_projection_lines(
                        ax,
                        R11=line.get("R11", np.nan),
                        R12=line.get("R12", np.nan),
                        sigma_over_sqrt_emit=line.get("sigma_over_sqrt_emit", np.nan),
                        color=color,
                        label=None,
                    )
                    legend_handles.append(Line2D([0], [0], color=color, linewidth=1.5, label=screen_label))
            if legend_handles:
                ax.legend(handles=legend_handles, fontsize=8)

            try:
                emit_norm_float = float(emit_norm)
            except Exception:
                emit_norm_float = np.nan

        draw_one(ax_x, beta_x, alpha_x, emit_x_norm, "Normalised Horizontal Phase Space")
        draw_one(ax_y, beta_y, alpha_y, emit_y_norm, "Normalised Vertical Phase Space")
        self.canvas.draw_idle()

    @staticmethod
    def _draw_projection_lines(ax, R11, R12, sigma_over_sqrt_emit, color=None, label=None):
        '''
        Draw screen projection constraints in normalised phase space:

            R11 * U + R12 * U' = ± sigma_screen / sqrt(epsilon)
        '''
        try:
            R11 = float(R11)
            R12 = float(R12)
            c = float(sigma_over_sqrt_emit)
        except Exception:
            return

        if not (np.isfinite(R11) and np.isfinite(R12) and np.isfinite(c)):
            return

        if abs(R12) < 1e-14:
            if abs(R11) < 1e-14:
                return
            x_const = c / R11
            ax.axvline(+x_const, color=color, linewidth=1.5, label=label)
            ax.axvline(-x_const, color=color, linewidth=1.5)
            return

        x_min, x_max = ax.get_xlim()
        if not (np.isfinite(x_min) and np.isfinite(x_max)) or x_min == x_max:
            x_min, x_max = -2.0, 2.0
        U = np.linspace(x_min, x_max, 200)
        Up_plus = (+c - R11 * U) / R12
        Up_minus = (-c - R11 * U) / R12
        ax.plot(U, Up_plus, linewidth=1.5, color=color, label=label)
        ax.plot(U, Up_minus, linewidth=1.5, color=color)

    def plot_projection_constraints(self, result, session, interface=None):
        self.figure.clear()
        screens = list(session.get("screens", []))
        reference_screen = session.get("reference_screen", screens[0])
        transport = session.get("measured_transport")
        if not isinstance(transport, dict) and interface is not None:
            get_transport = getattr(interface, "get_phase_space_transport_to_screens", None)
            if callable(get_transport):
                transport = get_transport(reference_screen=reference_screen, screens=screens)

        if not isinstance(transport, dict):
            self.plot_from_result(result, reference_name=session.get("reference_screen", session.get("quad_name", "fitted scan reference")))
            return

        sigx = np.asarray(session.get("sigx_mean", []), dtype=float)
        sigy = np.asarray(session.get("sigy_mean", []), dtype=float)

        if sigx.ndim == 2:
            sigx = sigx[0, :]
        if sigy.ndim == 2:
            sigy = sigy[0, :]

        emit_x_norm = float(result.get("emit_x_norm", np.nan))
        emit_y_norm = float(result.get("emit_y_norm", np.nan))

        gamma_rel, beta_rel, beta_gamma = interface.get_beam_factors()
        emit_x_geom = emit_x_norm / beta_gamma
        emit_y_geom = emit_y_norm / beta_gamma

        tx = transport["x"]
        ty = transport["y"]

        projection_lines = {"x": [], "y": []}
        color_cycle = ["blue", "lime", "red", "magenta", "orange", "cyan", "purple", "brown"]
        screen_colors = {str(screen): color_cycle[i % len(color_cycle)] for i, screen in enumerate(screens)}

        for i, screen in enumerate(screens):
            projection_lines["x"].append({
                "screen": screen,
                "R11": tx["R11"][i],
                "R12": tx["R12"][i],
                "sigma_over_sqrt_emit": sigx[i] / np.sqrt(emit_x_geom),
            })

            projection_lines["y"].append({
                "screen": screen,
                "R11": ty["R33"][i],
                "R12": ty["R34"][i],
                "sigma_over_sqrt_emit": sigy[i] / np.sqrt(emit_y_geom),
            })

        propagate_twiss = self.propagate_twiss

        beta_x0 = result.get("beta_x0", np.nan)
        alpha_x0 = result.get("alpha_x0", np.nan)
        beta_y0 = result.get("beta_y0", np.nan)
        alpha_y0 = result.get("alpha_y0", np.nan)

        # beta_x0/alpha_x0 are fitted at the quadrupole entrance (Optimization._fit_6d builds the
        # synthetic bunch there), not at reference_screen, so this propagation must start there too -
        # it is a separate transport call from the one used for the projection_lines above, which are
        # deliberately expressed in the reference_screen frame.
        quad_name = session.get("quad_name")
        transport_from_quad = transport
        if quad_name and quad_name != reference_screen:
            get_transport = getattr(interface, "get_phase_space_transport_to_screens", None)
            if callable(get_transport):
                transport_from_quad = get_transport(reference_screen=quad_name, screens=screens)
        txq, tyq = transport_from_quad["x"], transport_from_quad["y"]

        R21 = np.asarray(txq.get("R21", []), dtype=float)
        R22 = np.asarray(txq.get("R22", []), dtype=float)
        R43 = np.asarray(tyq.get("R43", []), dtype=float)
        R44 = np.asarray(tyq.get("R44", []), dtype=float)

        screens_twiss = []
        for i, screen in enumerate(screens):
            beta_x_i, alpha_x_i = np.nan, np.nan
            beta_y_i, alpha_y_i = np.nan, np.nan
            if i < len(txq["R11"]) and i < len(txq["R12"]) and i < R21.size and i < R22.size:
                beta_x_i, alpha_x_i = propagate_twiss(beta_x0, alpha_x0, txq["R11"][i], txq["R12"][i], R21[i], R22[i])
            if i < len(tyq["R33"]) and i < len(tyq["R34"]) and i < R43.size and i < R44.size:
                beta_y_i, alpha_y_i = propagate_twiss(beta_y0, alpha_y0, tyq["R33"][i], tyq["R34"][i], R43[i], R44[i])
            screens_twiss.append((
                str(screen),
                {
                    "beta_x": beta_x_i,
                    "beta_y": beta_y_i,
                    "alpha_x": alpha_x_i,
                    "alpha_y": alpha_y_i,
                    "emit_x_norm": result.get("emit_x_norm", np.nan),
                    "emit_y_norm": result.get("emit_y_norm", np.nan),
                }
            ))

        if len(screens_twiss) == 0:
            screens_twiss = [(
                str(quad_name or reference_screen),
                {
                    "beta_x": beta_x0,
                    "beta_y": beta_y0,
                    "alpha_x": alpha_x0,
                    "alpha_y": alpha_y0,
                    "emit_x_norm": result.get("emit_x_norm", np.nan),
                    "emit_y_norm": result.get("emit_y_norm", np.nan),
                }
            )]

        self._plot_data = {
            "result": result,
            "projection_lines": projection_lines,
            "projection_reference": str(reference_screen),
            "screen_colors": screen_colors,
            "screens_twiss": screens_twiss,
        }

        self.location_combo.blockSignals(True)
        self.location_combo.clear()
        for screen_name, _ in screens_twiss:
            self.location_combo.addItem(str(screen_name))
        if str(reference_screen) in [str(screen_name) for screen_name, _ in screens_twiss]:
            self.location_combo.setCurrentText(str(reference_screen))
        elif self.location_combo.count() > 0:
            self.location_combo.setCurrentIndex(0)
        self.location_combo.blockSignals(False)

        self._redraw_current_location()