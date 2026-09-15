import RF_Track as rft  # do not touch this import!
import os, sys, time
import numpy as np
from datetime import datetime, timezone, timezone
from enum import Enum

try:
    pyqt_version = 6
    from PyQt6 import uic
    from PyQt6.QtWidgets import QApplication, QMainWindow, QMessageBox, QVBoxLayout, QHBoxLayout, QListWidgetItem, \
        QStyledItemDelegate, QScrollArea, QFrame, QDialog, QDialogButtonBox, QFileDialog
    from PyQt6.QtCore import Qt, QTimer, QRect, QObject, QThread, pyqtSignal
    from PyQt6.QtGui import QPainter, QPixmap, QFont
except ImportError:
    pyqt_version = 5
    from PyQt5 import uic
    from PyQt5.QtWidgets import QApplication, QMainWindow, QMessageBox, QVBoxLayout, QHBoxLayout, QListWidgetItem, \
        QStyledItemDelegate, QScrollArea, QFrame, QDialog, QDialogButtonBox, QFormLayout, QFileDialog
    from PyQt5.QtCore import Qt, QTimer, QRect, QObject, QThread, pyqtSignal
    from PyQt5.QtGui import QPainter, QPixmap, QFont
import matplotlib

matplotlib.use("QtAgg")
from Interfaces.interface_setup import INTERFACE_SETUP
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from Backend.SaveOrLoad import SaveOrLoad
from Backend.State import State
from Backend.EmittanceComputingEngines.select_engine import EmittanceComputingEngineSelector
from Backend.EM_helpers.QuadrupoleScan import QuadrupoleScan
from Backend.LogConsole import LogConsole
# from Backend.EM_helpers.PhaseSpaceGraphs import PhaseSpaces
from Backend.EM_helpers.ShowBeamline import ShowBeamline
from Backend.EM_helpers.DisplayScreenImages import DisplayScreenImages
from Backend.EM_helpers.EmittanceEvolution import EmittanceEvolution
from Backend.EM_helpers.FitBounds import BoundsForParameter
from Backend.EM_helpers.ScanCurrentRanges import ScanCurrentRanges
from Backend.EM_helpers.ScanPointSelection import ScanPointSelection
from Backend.EM_helpers.QuadrupoleModelUpdate import load_model_quadrupoles, quadrupole_tracking_model, quadrupole_status_metadata, save_quadrupole_status

class ComputationMode(Enum):
    LRM = "Linear R-response model"
    ML = "Machine learning model"
    RFT = "RF-Track tracking"

class SPositionDelegate(QStyledItemDelegate):
    S_ROLE = int(Qt.ItemDataRole.UserRole) + 1

    def paint(self, painter: QPainter, option, index):
        painter.save()
        try:
            opt = option
            self.initStyleOption(opt, index)
            style = opt.widget.style() if opt.widget is not None else None
            if style is not None:
                opt_no_text = opt
                opt_no_text.text = ""
                style.drawControl(style.ControlElement.CE_ItemViewItem, opt_no_text, painter, opt.widget)

            device_name = str(index.data(Qt.ItemDataRole.UserRole) or index.data(Qt.ItemDataRole.DisplayRole) or "")
            s_text = str(index.data(self.S_ROLE) or "")
            r = opt.rect
            margin = 8
            painter.setFont(opt.font)
            painter.setPen(opt.palette.color(opt.palette.ColorRole.Text))

            fm = painter.fontMetrics()
            s_column_width = max(fm.horizontalAdvance("S = 000.000 m"), 90)

            left_rect = QRect(r.left() + margin, r.top(), max(10, r.width() - s_column_width - 3 * margin), r.height())
            right_rect = QRect(r.left() + r.width() - s_column_width - margin, r.top(), s_column_width, r.height())
            elided_name = fm.elidedText(device_name, Qt.TextElideMode.ElideRight, max(10, left_rect.width()))
            painter.drawText(left_rect, int(Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft), elided_name)
            if s_text:
                painter.drawText(right_rect, int(Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft), s_text)
        finally:
            painter.restore()


class MatplotlibWidget(FigureCanvas):
    def __init__(self, parent=None):
        self.figure = Figure(figsize=(6, 4), tight_layout=True)
        super().__init__(self.figure)
        self.setParent(parent)
 
class OptimizationWorker(QObject):
    finished = pyqtSignal(object)
    error = pyqtSignal(str)
    optimizer_ready = pyqtSignal(object)
    done = pyqtSignal()
    progress = pyqtSignal(str, int, int)
    info = pyqtSignal(str)
    ml_not_found_rft_fallback = pyqtSignal()

    def __init__(self, interface, session, selected_screens=None, bounds=None, fit_quadrupole_strength=False,
                 fit_quad_offset=False, fit_quad_roll=False, fit_energy_pref = None, computing_method="Linear R-response model"):
        super().__init__()
        self.interface = interface
        self.session = session
        self.selected_screens = list(selected_screens or [])
        self.bounds = bounds
        self.fit_quadrupole_strength = bool(fit_quadrupole_strength)
        self.fit_quad_offset = bool(fit_quad_offset)
        self.fit_quad_roll = bool(fit_quad_roll)
        self.fit_energy_pref = bool(fit_energy_pref)
        self.computing_method = computing_method

    def _emit_progress(self, phase, current, total):
        self.progress.emit(str(phase), int(current), int(total))

    def _get_interface_initial_settings(self):
        interface = getattr(self.interface, "interface", self.interface)
        interface_class_name = interface.__class__.__name__
        interface_module_name = interface.__class__.__module__

        for machine_name, machine_interfaces in INTERFACE_SETUP.items():
            for interface_defaults in machine_interfaces:
                if (interface_defaults.get("class_name") == interface_class_name) and (
                        interface_defaults.get("module") == interface_module_name):
                    return dict(interface_defaults, machine_name=str(machine_name))
        return None

    def _get_interface_bounds(self):
        if self.bounds is not None:
            return dict(self.bounds)
        interface_defaults = self._get_interface_initial_settings()
        if interface_defaults is None:
            return {}
        return dict(interface_defaults.get("bounds", {}))

    def _cut_session_to_detected_devices(self):
        if self.session is None:
            return None
        selected_screens = self.selected_screens
        if not selected_screens:
            raise ValueError("Select at least one screen")
        session_screens = list(self.session.get("screens", []))
        matched_screens = [self.interface.match_screen_name(screen, session_screens) for screen in selected_screens]
        selected_indices = [session_screens.index(screen) for screen in matched_screens if screen is not None]
        if not selected_indices:
            raise ValueError(
                f"None of the selected screens {list(selected_screens)} are present in the loaded session data {session_screens}.")
        cut_session = dict(self.session)
        cut_session["screens"] = [session_screens[i] for i in selected_indices]

        for key in ("sigx_mean", "sigy_mean", "sigxy_mean", "sigx_std", "sigy_std", "sigxy_std", "sigx_shots",
                    "sigy_shots"):
            values = np.asarray(self.session[key], dtype=float)
            cut_session[key] = values[:, selected_indices, ...].tolist()
        for key in ("sigxy_shots", "x_mean", "y_mean", "x_std", "y_std", "x_shots", "y_shots"):
            if key not in self.session:
                continue
            values = np.asarray(self.session[key], dtype=float)
            cut_session[key] = values[:, selected_indices, ...].tolist()
        cut_session["images"] = [[self.session["images"][step_index][screen_index] for screen_index in selected_indices]
                                 for step_index in range(len(self.session["images"]))]
        cut_session["nscreens"] = len(selected_indices)
        reference_screen = cut_session.get("reference_screen")
        if reference_screen not in cut_session["screens"]:
            cut_session["reference_screen"] = cut_session["screens"][0]
        return cut_session

    def run(self):
        try:
            interface_defaults = self._get_interface_initial_settings() or {}
            machine_name = str(interface_defaults.get("machine_name", ""))
            bounds = self._get_interface_bounds()
            session_for_opt = self._cut_session_to_detected_devices()
            tool = EmittanceComputingEngineSelector.create(method=self.computing_method, interface=self.interface,
                                                           session=session_for_opt, machine_name=machine_name,
                                                           info_callback=self.info.emit,
                                                           fit_quadrupole_strength=self.fit_quadrupole_strength,
                                                           fit_quad_offset=self.fit_quad_offset,
                                                           fit_quad_roll=self.fit_quad_roll,
                                                           fit_energy_pref = self.fit_energy_pref,
                                                           progress_callback=self._emit_progress,
                                                           fallback_callback=self.ml_not_found_rft_fallback.emit)
            self.optimizer_ready.emit(tool)
            output = tool.fit_from_session(session_for_opt, bounds=bounds)
            self.finished.emit(output)
        except Exception as e:
            self.error.emit(str(e))
        finally:
            self.done.emit()


class MainWindow(QMainWindow, QuadrupoleScan):
    def __init__(self, interface, dir_name, is_simulation, bg_shots=10):
        super().__init__()
        self.interface = interface
        self.dir_name = dir_name
        self._general_session_dir = dir_name
        self.is_simulation = is_simulation
        self.session = None
        ui_path = os.path.join(os.path.dirname(__file__), "UI files/Emittance_Measurement_GUI.ui")
        uic.loadUi(ui_path, self)
        self._make_settings_panel_scrollable()
        QTimer.singleShot(0, self.showMaximized)
        self._load_logo()
        self.load_session_button.clicked.connect(self.load_scan_and_optimization_settings)
        self.session_directory.setText(dir_name)
        self.start_optimization_button.clicked.connect(self._run_optimization)
        self.stop_optimization_button.clicked.connect(self._stop_optimization)
        self.setWindowTitle("Emittance Measurement GUI")
        self.fitResultsVBox.setStretch(0, 0)
        self.fitResultsVBox.setStretch(1, 1)
        self.progressBar.setValue(0)
        self.quadrupoles_list.setItemDelegate(SPositionDelegate(self.quadrupoles_list))
        self.screens_list.setItemDelegate(SPositionDelegate(self.screens_list))
        self._optimization_t0 = None
        self._scan_stop_requested = False
        self._is_scanning = False
        self._is_optimizing = False
        self._current_optimizer = None
        self._optimization_thread = None
        self._optimization_worker = None
        self.canvas = MatplotlibWidget(self.plotPlaceholder)
        layout = self.plotPlaceholder.layout()
        if layout is None:
            layout = QVBoxLayout(self.plotPlaceholder)
        layout.addWidget(self.canvas)
        quadrupoles = list(getattr(self.interface, "quadrupoles", []))
        screens = list(getattr(self.interface, "screens", []))
        screen_order, screen_order_type = self._get_element_order_values(screens)
        screen_pairs = sorted(zip(screens, screen_order),
                              key=lambda x: x[1] if np.isfinite(x[1]) else np.inf)  # assigns S position to each screen
        screens_sorted = [name for name, _ in screen_pairs]  # only names
        self._show_s_values_and_device_lists(self.quadrupoles_list, quadrupoles)
        self._show_s_values_and_device_lists(self.screens_list, screens_sorted)
        self.show_scan_on_all_screens.toggled.connect(self._on_show_all_screens_toggled)
        self.screen_on_plot.setEnabled(not self.show_scan_on_all_screens.isChecked())
        self.quad_on_plot.clear()
        self.quad_on_plot.addItems(quadrupoles)
        self.screen_on_plot.clear()
        self.screen_on_plot.addItems(screens_sorted)
        self.start_button_scan.clicked.connect(lambda _checked=False: self._run_scan())
        self.stop_button_scan.clicked.connect(self._stop_scan)
        self.quad_on_plot.currentIndexChanged.connect(lambda _=None: self._draw_live_scan(self.session))
        self.screen_on_plot.currentIndexChanged.connect(lambda _=None: self._draw_live_scan(self.session))
        self._set_progress(0)
        self._clear_fit_panel()
        self._reset_canvas()
        self.screens_list.itemSelectionChanged.connect(self._screen_selection_changed)
        self._last_selected_quadrupoles = []
        self._filter_quadrupoles_in_gui()
        self.clear_plots_button.clicked.connect(self._clear_plots)
        self.log_console = None
        self.phase_spaces = None
        self.screen_images = None
        self.beamline_view = None
        self.beta_evolution_window = None
        self.emittance_evolution_window = None
        self.log_console_button.clicked.connect(self._show_console_log)
        # self.phase_spaces_button.clicked.connect(self._show_phase_spaces)
        self.display_screen_images_button.clicked.connect(self._show_screen_images)
        # self.beta_function_button.clicked.connect(self._show_beta_function_evolution)
        self.emittance_evolution_button.clicked.connect(self._show_emittance_evolution)
        self.pause_button.clicked.connect(self._pause_task)
        self.resume_button.clicked.connect(self._resume_task)
        self._scan_pause_requested = False
        self._scan_is_paused = False
        self._optimization_paused = False
        self._interrupted_scan = None
        self._last_scan_status = None
        self.fit_core_params = ("emit_x_norm", "beta_x0", "alpha_x0", "emit_y_norm", "beta_y0", "alpha_y0")
        self.additional_params = ("quad_k1l_0", "quad_dx0", "quad_dy0", "quad_roll", "energy_pref")
        self.additional_params_scales = {"quad_k1l_0": 1.0, "quad_dx0": 1e-3, "quad_dy0": 1e-3, "quad_roll": 1e-3, "energy_pref": 1}
        self.additional_params_defaults = {"quad_k1l_0": (0.0, 0.0), "quad_dx0": (-2.0, 2.0), "quad_dy0": (-2.0, 2.0), "quad_roll": (-50.0, 50.0), "energy_pref": (150, 250)}
        self._setup_bounds_buttons()
        self._populate_default_bounds()
        self.fit_quadrupole_strength_checkbox.toggled.connect(self._update_additional_fit_controls)
        self.fit_quad_offset_checkbox.toggled.connect(self._update_additional_fit_controls)
        self.fit_quad_roll_checkbox.toggled.connect(self._update_additional_fit_controls)
        self.fit_energy_checkbox.toggled.connect(self._update_additional_fit_controls)
        self.computation_mode = ComputationMode(self.computing_method_combo.currentText())
        self.computing_method_combo.currentTextChanged.connect(self._on_computation_mode_changed)
        self.steps_settings.valueChanged.connect(self._on_nsteps_scan_changed)
        self._on_computation_mode_changed(self.computing_method_combo.currentText())
        self._on_nsteps_scan_changed(self.steps_settings.value())
        self.load_screens_data_button.clicked.connect(self._load_screens_data)
        self.background_shots.setValue(bg_shots)
        self.interface.bg_shots = int(self.background_shots.value())
        self.background_shots.valueChanged.connect(self._on_bg_shots_changed)
        self.show_beamline_button.clicked.connect(self._show_beamline)
        if self.is_simulation == True:
            self.download_quads_button.setEnabled(False)
        else:
            self.download_quads_button.setEnabled(True)
            self.download_quads_button.clicked.connect(self._download_all_quads_status)
        self.quadrupoles_status = None
        self.model_quadrupoles_status = None
        self._quadrupoles_status_payload = None
        self._model_quadrupoles_status_payload = None
        self.load_machine_status_button.clicked.connect(self._load_model_quadrupoles_status)
        try:
            quadrupole_tracking_model(self.interface)
        except ValueError as e:
            self.load_machine_status_button.setEnabled(False)
            self.load_machine_status_button.setToolTip(str(e))
        self._screen_current_ranges = {}
        self._excluded_points = set()
        self.delete_point_button.clicked.connect(self._edit_excluded_points)
        self.per_screen_ranges_button.clicked.connect(self._edit_per_screen_current_ranges)
        self.quadrupoles_list.itemSelectionChanged.connect(self._update_quad_readback_label)
        self.minimum_current.valueChanged.connect(lambda _=None: self._update_per_screen_ranges_button())
        self.maximum_current.valueChanged.connect(lambda _=None: self._update_per_screen_ranges_button())
        self._update_per_screen_ranges_button()
        self._update_quad_readback_label()

    def _make_settings_panel_scrollable(self):
        main_layout = self.centralwidget.layout()
        settings_panel = self.leftGroup
        settings_scroll = QScrollArea(self.centralwidget)
        settings_scroll.setWidgetResizable(True)
        settings_scroll.setFrameShape(QFrame.Shape.NoFrame if pyqt_version == 6 else QFrame.NoFrame)
        settings_scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAsNeeded if pyqt_version == 6 else Qt.ScrollBarAsNeeded)
        settings_scroll.setVerticalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAsNeeded if pyqt_version == 6 else Qt.ScrollBarAsNeeded)
        main_layout.replaceWidget(settings_panel, settings_scroll)
        settings_scroll.setWidget(settings_panel)
        main_layout.setStretch(main_layout.indexOf(settings_scroll), 0)
        main_layout.setStretch(main_layout.indexOf(self.tabs), 1)
        self.settings_scroll = settings_scroll
        self._settings_sections = []
        while self.leftVBox.count():
            item = self.leftVBox.takeAt(0)
            if item.widget() is not None:
                self._settings_sections.append(item.widget())
        for section in self._settings_sections:
            self.leftVBox.addWidget(section)
        self._update_settings_panel_width()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._update_settings_panel_width()

    def _update_settings_panel_width(self):
        if not hasattr(self, "settings_scroll"):
            return
        if self.width() >= 1900:
            width = self.leftGroup.minimumSizeHint().width()
            self.settings_scroll.setMinimumWidth(width)
            self.settings_scroll.setMaximumWidth(width)
        else:
            self.settings_scroll.setMinimumWidth(0)
            self.settings_scroll.setMaximumWidth(16777215)

    def _download_all_quads_status(self):
        try:
            target_dir = self.session_directory.text().strip() or self.dir_name
            output_file_name = os.path.join(target_dir, "quadrupoles_status.npz")
            quadrupoles = list(getattr(self.interface, "quadrupoles", []))
            self.log(f"Reading {len(quadrupoles)} quadrupole currents from the machine.")
            payload = dict(self.interface.get_quadrupoles(names=quadrupoles))
            payload["captured_at_utc"] = datetime.now(timezone.utc).isoformat()
            payload["source_interface"] = self.interface.get_name()
            momentum = getattr(getattr(self.interface, "tracking_interface", self.interface), "Pref", None)
            if momentum is not None and np.isfinite(float(momentum)):
                payload["reference_momentum_MeV_c"] = float(momentum)
            output_file_name = save_quadrupole_status(payload, output_file_name)
            self._quadrupoles_status_payload = payload
            self.quadrupoles_status = quadrupole_status_metadata(payload, output_file_name)
            self.quadrupoles_status["saved_path"] = output_file_name
            if self.session is not None:
                self.session["quadrupoles_status"] = dict(self.quadrupoles_status)
            self.log(f"Saved quadrupole readbacks to {output_file_name} at {payload['captured_at_utc']}.")

        except Exception as e:
            QMessageBox.information(self, "Save quadrupoles status",
                                    f"An error occured while trying to save quadrupoles status. {e}")
            return

    def _load_model_quadrupoles_status(self):
        print("Quadrupoles currents before RFTrack model update:")
        print(self.interface.get_quadrupoles()["bact"])
        if self._is_scanning or self._is_optimizing:
            QMessageBox.information(self, "Load quadrupoles status",
                                    "Stop the scan or fit before changing the model quadrupoles.")
            return
        default_dir = self.load_screens_data_database.text().strip() or self.session_directory.text().strip() or self.dir_name
        source_path, _ = QFileDialog.getOpenFileName(self, "Load quadrupoles status", default_dir,
                                                     "Quadrupole snapshots (*.npz)")
        if not source_path:
            return
        try:
            payload, metadata = load_model_quadrupoles(self.interface, source_path)
        except Exception as e:
            QMessageBox.warning(self, "Load quadrupoles status", str(e))
            return
        self._model_quadrupoles_status_payload = payload
        self.model_quadrupoles_status = metadata
        if self.session is not None:
            self.session["model_quadrupoles_status"] = dict(metadata)
        timestamp = metadata["captured_at_utc"] or "None"
        print(f"Loaded {metadata['quadrupole_count']} quadrupoles into the RF-Track model")

        print("Quadrupoles currents after RFTrack model update:")
        print(self.interface.get_quadrupoles()["bact"])

    def _preserve_quadrupole_status_files(self, target_dir):
        for attribute, filename in (("quadrupoles_status", "quadrupoles_status.npz"),
                                    ("model_quadrupoles_status", "model_quadrupoles_status.npz")):
            payload = getattr(self, f"_{attribute}_payload", None)
            metadata = getattr(self, attribute, None)
            if payload is None or metadata is None:
                continue
            target_path = os.path.join(target_dir, filename)
            previous_path = metadata.get("saved_path")
            if previous_path and os.path.dirname(previous_path) == os.path.abspath(target_dir):
                target_path = previous_path
            saved_path = save_quadrupole_status(payload, target_path)
            metadata["saved_path"] = saved_path
            if self.session is not None:
                self.session[attribute] = dict(metadata)

    def _show_beamline(self):
        selected_quadrupole, screens = self._get_selection()
        if self.beamline_view is None:
            self.beamline_view = ShowBeamline(interface=self.interface, parent=self, quad_selected=selected_quadrupole,
                                              screens=screens)
        self.beamline_view._display_beamline_view()
        self.beamline_view.show()
        self.beamline_view.raise_()
        self.beamline_view.activateWindow()

    def _show_emittance_evolution(self):
        _, screens = self._get_selection()
        if self.emittance_evolution_window is None:
            self.emittance_evolution_window = EmittanceEvolution(interface=self.interface, parent=self, screens=screens)
        session = self.session if isinstance(self.session, dict) else {}
        self.emittance_evolution_window._display_emittance_evolution(
            screens=screens, session=session, result=session.get("optimization_result"))
        self.emittance_evolution_window.show()
        self.emittance_evolution_window.raise_()
        self.emittance_evolution_window.activateWindow()

    def _on_bg_shots_changed(self, value):
        self.interface.bg_shots = max(0, int(value))

    def _load_screens_data(self):
        loaded_states = self.load_screens_data()
        if loaded_states is None:
            return
        self.session = self._get_session_data_from_database()
        if self.session is None:
            QMessageBox.information(self, "Emittance Measurement Session Error", "Session not found.")
            return
        self._set_default_quad_strength_bounds_from_session(self.session)
        self._refresh_plot_comboboxes_from_session(self.session)
        self._draw_live_scan(self.session)
        self.steps_settings.setValue(int(self.session["steps"]))
        self.meas_per_step.setValue(int(self.session["nshots"]))
        self.steps_settings.setEnabled(False)
        self.meas_per_step.setEnabled(False)
        self.quadrupoles_list.setEnabled(False)

    def _on_nsteps_scan_changed(self, nsteps_settings):
        n_scan_steps = nsteps_settings
        is_steps_zero = bool(n_scan_steps == 0)
        self.quadrupoles_list.setEnabled(not is_steps_zero)
        if is_steps_zero:
            self.start_button_scan.setText("ACQUIRE SCREEN DATA")
        else:
            self.start_button_scan.setText("START SCAN")

    def _get_interface_initial_settings(self):
        interface_class_name = self.interface.__class__.__name__
        interface_module_name = self.interface.__class__.__module__

        for machine_interfaces in INTERFACE_SETUP.values():
            for interface_defaults in machine_interfaces:
                if (interface_defaults.get("class_name") == interface_class_name) and (
                        interface_defaults.get("module") == interface_module_name):
                    return interface_defaults
        return None

    def _get_interface_units(self):
        interface_defaults = self._get_interface_initial_settings()
        if interface_defaults is None:
            return {}, 0.01, "mm", ""
        units_settings = interface_defaults.get("units", {})
        em_sigma_unit = units_settings.get("em_sigma_unit", "mm")
        return em_sigma_unit

    def _get_interface_default_bounds(self):
        interface_defaults = self._get_interface_initial_settings()
        if interface_defaults is None:
            return {}
        return dict(interface_defaults.get("bounds", {}))

    def _setup_bounds_buttons(self):
        self.bounds_buttons = {
            "emit_x_norm": self.eps_h_bounds,
            "beta_x0": self.beta_h_bounds,
            "alpha_x0": self.alpha_h_bounds,
            "emit_y_norm": self.eps_v_bounds,
            "beta_y0": self.beta_v_bounds,
            "alpha_y0": self.alpha_v_bounds,
            "quad_k1l_0": self.quad_strength_bounds,
            "quad_dx0": self.quad_dx_bounds,
            "quad_dy0": self.quad_dy_bounds,
            "quad_roll": self.quad_roll_bounds,
            "energy_pref": self.energy_pref_bounds,
        }
        self.bounds_units = {
            "emit_x_norm": "mm·mrad", "beta_x0": "m", "alpha_x0": "",
            "emit_y_norm": "mm·mrad", "beta_y0": "m", "alpha_y0": "",
            "quad_k1l_0": self._quadrupole_value_unit(), "quad_dx0": "mm", "quad_dy0": "mm", "quad_roll": "mrad",
            "energy_pref": "MeV/c"

        }
        self.bounds_display_names = {
            "emit_x_norm": "Horizontal ε (norm.)",
            "beta_x0": "Horizontal β",
            "alpha_x0": "Horizontal α",
            "emit_y_norm": "Vertical ε (norm.)",
            "beta_y0": "Vertical β",
            "alpha_y0": "Vertical α",
            "quad_k1l_0": "Quadrupole strength",
            "quad_dx0": "Quadrupole Δx",
            "quad_dy0": "Quadrupole Δy",
            "quad_roll": "Quadrupole roll",
            "energy_pref": "Beam energy"
        }
        self._bounds_button_base_text = {param: button.text() for param, button in self.bounds_buttons.items()}
        quad_value_unit = self.bounds_units["quad_k1l_0"]
        self._bounds_button_base_text[
            "quad_k1l_0"] = "I₀ [A]" if quad_value_unit == "A" else f"K1L₀ [{quad_value_unit}]"
        self._bounds_values = {}
        for param, button in self.bounds_buttons.items():
            button.clicked.connect(lambda _checked=False, parameter=param: self._edit_bounds(parameter))

    def _set_bounds_value(self, param, low, high):
        low, high = float(low), float(high)
        self._bounds_values[param] = [min(low, high), max(low, high)]
        self._refresh_bounds_button(param)

    def _refresh_bounds_button(self, param):
        button = self.bounds_buttons.get(param)
        if button is None:
            return
        low, high = self._bounds_values.get(param, (0.0, 0.0))
        base_text = self._bounds_button_base_text.get(param, param)
        button.setText(f"{base_text} ({low:.4g}, {high:.4g})")
        unit = self.bounds_units.get(param, "")
        button.setToolTip(
            f"{self.bounds_display_names.get(param, param)}: ({low:g}, {high:g}) {unit}".strip() + "\nClick to change the fit search bounds.")

    def _bounds_hint(self, param):
        if param == "quad_k1l_0":
            return self._quad_readback_text()
        return ""

    def _edit_bounds(self, param):
        if param == "quad_k1l_0":
            self.bounds_units[param] = self._quadrupole_value_unit()
        low, high = self._bounds_values.get(param, (0.0, 0.0))
        dialog = BoundsForParameter(self.bounds_display_names.get(param, param), parent=self,
                                    unit=self.bounds_units.get(param, ""), hint=self._bounds_hint(param))
        dialog.set_values(float(low), float(high))
        if not dialog.exec():
            return
        low, high = dialog.get_values()
        if np.isclose(low, high):
            QMessageBox.information(self, "Fit bounds", "Lower and upper bound must be different.")
            return
        self._set_bounds_value(param, low, high)

    def _populate_default_bounds(self):
        defaults = self._get_interface_default_bounds()
        for param in self.fit_core_params:
            low, high = defaults.get(param, (0.0, 0.0))
            self._set_bounds_value(param, low, high)
        for param in self.additional_params:
            low, high = defaults.get(param, self.additional_params_defaults[param])
            self._set_bounds_value(param, low, high)

    def _set_default_quad_strength_bounds_from_session(self, session):
        if session is None:
            return
        unit = self._quadrupole_value_unit(session)
        self.bounds_units["quad_k1l_0"] = unit
        self._bounds_button_base_text["quad_k1l_0"] = "I₀ [A]" if unit == "A" else f"K1L₀ [{unit}]"
        k1l_0 = float(session.get("K1L_0", np.nan))
        if not np.isfinite(k1l_0) or np.isclose(k1l_0, 0.0):
            scanned_values = np.asarray(session.get("K1L_values", []), dtype=float)
            k1l_0 = float(np.nanmedian(scanned_values)) if scanned_values.size else np.nan
        if not np.isfinite(k1l_0):
            self._refresh_bounds_button("quad_k1l_0")
            return
        low, high = sorted((0.7 * k1l_0, 1.3 * k1l_0))
        self._set_bounds_value("quad_k1l_0", low, high)

    def _set_bound_row_enabled(self, param, enabled):
        button = self.bounds_buttons.get(param)
        if button is not None:
            button.setEnabled(bool(enabled))

    def _update_additional_fit_controls(self, _checked=None):
        is_linear_mode = self.computation_mode == ComputationMode.LRM
        for checkbox in (self.fit_quadrupole_strength_checkbox, self.fit_quad_offset_checkbox, self.fit_quad_roll_checkbox, self.fit_energy_checkbox):
            checkbox.setEnabled(not is_linear_mode)
        self._set_bound_row_enabled("quad_k1l_0", not is_linear_mode and self.fit_quadrupole_strength_checkbox.isChecked())
        offset_bounds_enabled = not is_linear_mode and self.fit_quad_offset_checkbox.isChecked()
        self._set_bound_row_enabled("quad_dx0", offset_bounds_enabled)
        self._set_bound_row_enabled("quad_dy0", offset_bounds_enabled)
        self._set_bound_row_enabled("quad_roll", not is_linear_mode and self.fit_quad_roll_checkbox.isChecked())
        self._set_bound_row_enabled("energy_pref", not is_linear_mode and self.fit_energy_checkbox.isChecked())

    def _set_bounds_from_saved_settings(self, saved_bounds):
        for param in self.fit_core_params:
            if param in saved_bounds:
                low, high = saved_bounds[param]
                self._set_bounds_value(param, low, high)

        for param in self.additional_params:
            if param in saved_bounds:
                low, high = saved_bounds[param]
                scale = self.additional_params_scales[param]
                self._set_bounds_value(param, float(low) / scale, float(high) / scale)

    def _get_bounds_from_gui(self):
        bounds = {param: list(self._bounds_values.get(param, (0.0, 0.0))) for param in self.fit_core_params}
        checkbox_by_param = {
            "quad_k1l_0": self.fit_quadrupole_strength_checkbox,
            "quad_dx0": self.fit_quad_offset_checkbox,
            "quad_dy0": self.fit_quad_offset_checkbox,
            "quad_roll": self.fit_quad_roll_checkbox,
            "energy_pref": self.fit_energy_checkbox,
        }
        for param in self.additional_params:
            if checkbox_by_param[param].isChecked():
                scale = self.additional_params_scales[param]
                low, high = self._bounds_values.get(param, (0.0, 0.0))
                bounds[param] = [float(low) * scale, float(high) * scale]
        return bounds

    def _scan_points_of_session(self, session):
        session = self._get_session_for_selected_quad(session)
        if not isinstance(session, dict):
            return []
        quad_values = np.asarray(session.get("K1L_values", []), dtype=float)
        currents = np.asarray(session.get("current_values", []), dtype=float)
        screens = list(session.get("screens", []))
        sigx = np.asarray(session.get("sigx_mean", []), dtype=float)
        sigy = np.asarray(session.get("sigy_mean", []), dtype=float)
        unit = self._quadrupole_value_unit(session)
        points = []
        for step_index in range(len(quad_values)):
            for screen_index, screen in enumerate(screens):
                if sigx.ndim != 2 or step_index >= sigx.shape[0] or screen_index >= sigx.shape[1]:
                    continue
                if not np.isfinite(sigx[step_index, screen_index]) and not np.isfinite(sigy[step_index, screen_index]):
                    if (step_index, screen_index) not in self._excluded_points:
                        continue
                quad_value = currents[step_index] if unit != "A" and step_index < currents.size else quad_values[
                    step_index]
                points.append({"step": step_index, "screen_index": screen_index, "screen": screen,
                               "quad_value": quad_values[step_index], "quad_unit": unit,
                               "sigx": sigx[step_index, screen_index], "sigy": sigy[step_index, screen_index]})
        return points

    def _edit_excluded_points(self):
        if self.session is None:
            QMessageBox.information(self, "Delete scan point", "Run or load a scan first.")
            return
        points = self._scan_points_of_session(self.session)
        if not points:
            QMessageBox.information(self, "Delete scan point", "This session has no scan points to choose from.")
            return
        dialog = ScanPointSelection(points, self._excluded_points, parent=self)
        if not dialog.exec():
            return
        self._excluded_points = dialog.get_excluded()
        self._update_delete_point_button()
        self._draw_live_scan(self.session)

    def _update_delete_point_button(self):
        excluded = getattr(self, "_excluded_points", set())
        if not excluded:
            self.delete_point_button.setText("Delete chosen point")
            return
        self.delete_point_button.setText(f"Delete chosen point ({len(excluded)} removed)")

    def _session_without_excluded_points(self, session):
        excluded = getattr(self, "_excluded_points", set())
        if not isinstance(session, dict) or not excluded:
            return session
        masked = dict(session)
        target = masked
        if masked.get("mode") == "multi_quad_scan":
            per_quad = list(masked.get("per_quad_sessions", []))
            index = int(self.quad_on_plot.currentIndex())
            if not (0 <= index < len(per_quad)):
                return masked
            target = dict(per_quad[index])
            per_quad[index] = target
            masked["per_quad_sessions"] = per_quad
        for key in ("sigx_mean", "sigy_mean", "sigxy_mean", "x_mean", "y_mean",
                    "sigx_std", "sigy_std", "sigxy_std", "x_std", "y_std",
                    "sigx_shots", "sigy_shots", "sigxy_shots", "x_shots", "y_shots"):
            values = target.get(key)
            if values is None:
                continue
            array = np.asarray(values, dtype=float)
            if array.ndim < 2:
                continue
            for step_index, screen_index in excluded:
                if step_index < array.shape[0] and screen_index < array.shape[1]:
                    array[step_index, screen_index] = np.nan
            target[key] = array.tolist()
        return masked

    def _quad_readback_text(self, quad_names=None):
        if quad_names is None:
            quad_names, _ = self._get_selection()
        quad_names = list(quad_names or [])
        if not quad_names:
            return "Quadrupole readback: no quadrupole selected."
        quad_name = quad_names[0]
        try:
            quadrupoles = self.interface.get_quadrupoles([quad_name])
        except Exception as e:
            return f"{quad_name} readback: unavailable ({e})"
        unit = str(quadrupoles.get("value_unit", "1/m"))
        bdes = np.asarray(quadrupoles.get("bdes", []), dtype=float)
        if bdes.size == 0 or not np.isfinite(bdes[0]):
            return f"{quad_name} readback: not available."
        value = float(bdes[0])
        if unit == "A":
            return f"{quad_name} readback: {value:.3f} A"
        iact = np.asarray(quadrupoles.get("iact", []), dtype=float)
        current = float(iact[0]) if iact.size and np.isfinite(iact[0]) else self.quad_setpoint_to_current(quad_name,
                                                                                                          value, unit)
        current_text = f" = {current:.3f} A" if np.isfinite(current) else " (no A calibration)"
        return f"{quad_name} readback: K1L = {value:.5g} {unit}{current_text}"

    def _update_quad_readback_label(self):
        self.quad_readback_label.setText(self._quad_readback_text())

    def _edit_per_screen_current_ranges(self):
        _, screens = self._get_selection()
        if not screens:
            QMessageBox.information(self, "Scan current range", "No screens available.")
            return
        default_range = (float(self.minimum_current.value()), float(self.maximum_current.value()))
        dialog = ScanCurrentRanges(screens=screens, default_range=default_range,
                                   ranges=getattr(self, "_screen_current_ranges", {}),
                                   steps=int(self.steps_settings.value()), parent=self)
        if not dialog.exec():
            return
        self._screen_current_ranges = dialog.get_values()
        self._update_per_screen_ranges_button()

    def _update_per_screen_ranges_button(self):
        _, screens = self._get_selection()
        default_range = (float(self.minimum_current.value()), float(self.maximum_current.value()))
        ranges = {screen: values for screen, values in getattr(self, "_screen_current_ranges", {}).items()
                  if screen in screens and not np.allclose(values, default_range)}
        if not ranges:
            self.per_screen_ranges_button.setText("Current range per screen...")
            self.per_screen_ranges_button.setToolTip("Every screen is scanned over the current range set above.")
            return
        self.per_screen_ranges_button.setText(f"Current range per screen ({len(ranges)} changed)...")
        self.per_screen_ranges_button.setToolTip(
            "\n".join(f"{screen}: ({values[0]:g}, {values[1]:g}) A" for screen, values in ranges.items()))

    def _on_computation_mode_changed(self, text):
        self.computation_mode = ComputationMode(text)
        is_linear_mode = self.computation_mode == ComputationMode.LRM
        if self.steps_settings.isEnabled():
            if is_linear_mode:
                if self.steps_settings.value() > 0:
                    self._scan_steps_before_linear_mode = self.steps_settings.value()
                self.steps_settings.setValue(0)
                self.quadrupoles_list.setEnabled(True)
            elif self.steps_settings.value() == 0:
                self.steps_settings.setValue(getattr(self, "_scan_steps_before_linear_mode", 5))

        widgets_to_disable = [self.boundsSettingsGroup, self.bounds_quad_group, self.localOptimizationSettingsGroup]
        for widget in widgets_to_disable:
            widget.setEnabled(not is_linear_mode)
        self._update_additional_fit_controls()

    def _on_show_all_screens_toggled(self, checked):
        self.screen_on_plot.setEnabled(not bool(checked))
        self._draw_live_scan(self.session)

    def _pause_task(self):
        if self._is_scanning:
            self.log("Pausing scan...")
            self._scan_pause_requested = True
            return

        if self._is_optimizing and self._current_optimizer is not None:
            self.log("Pausing optimization...")
            self._current_optimizer.request_pause()
            return

    def _resume_task(self):
        if self._is_scanning and (self._scan_pause_requested or self._scan_is_paused):
            self.log("Resuming scan...")
            self._scan_pause_requested = False
            self._scan_is_paused = False
            return
        if not self._is_scanning and getattr(self, "_interrupted_scan", None):
            self._resume_interrupted_scan()
            return
        if self._optimization_paused and not self._is_optimizing and self.session is not None:
            self.log("Resuming optimization...")
            self._optimization_paused = False
            self._run_optimization()
            return

    def _stop_scan(self):
        self.log("Stopping scan...")
        if self._is_scanning:
            self._scan_stop_requested = True
            self._scan_pause_requested = False
            self._scan_is_paused = False

    def _clear_plots(self):
        self.session = None
        self._interrupted_scan = None
        self._excluded_points = set()
        self._update_delete_point_button()
        self._scan_stop_requested = False
        self._scan_pause_requested = False
        self._optimization_paused = False
        self._scan_is_paused = False
        self.quad_on_plot.clear()
        self.screen_on_plot.clear()
        quadrupoles = list(getattr(self.interface, "quadrupoles", []))
        # screens_data = self.interface.get_screens()
        screens = list(getattr(self.interface, "screens", []))
        screen_order, screen_order_type = self._get_element_order_values(screens)
        screen_pairs = sorted(zip(screens, screen_order), key=lambda x: x[1] if np.isfinite(x[1]) else np.inf)
        screens_sorted = [name for name, _ in screen_pairs]
        self.quad_on_plot.addItems(quadrupoles)
        self.screen_on_plot.addItems(screens_sorted)
        self._clear_fit_panel()
        self._reset_canvas()
        self._set_progress(0)

    def _stop_optimization(self):
        self.log("Stopping optimization...")
        self._optimization_paused = False
        if self._is_optimizing and self._current_optimizer is not None:
            self._current_optimizer.request_stop()

    def _get_element_order_values(self, names):
        names = list(names)
        s_positions = self._get_twiss_s_positions(names)
        s_values = []
        for value in s_positions:
            try:
                value = float(value)
            except (ValueError, TypeError):
                value = np.nan
            s_values.append(value)

        if any(np.isfinite(s_values)):
            return s_values, "S"
        sequence_indices = None
        try:
            sequence_indices = self.interface.get_elements_indices(names)
        except Exception:
            sequence_indices = None
        if sequence_indices is not None:
            try:
                sequence_indices = list(sequence_indices)
            except TypeError:
                sequence_indices = None
        if sequence_indices is not None and len(sequence_indices) == len(names):
            index_values = []
            for value in sequence_indices:
                try:
                    value = float(value)
                except (ValueError, TypeError):
                    value = np.nan
                index_values.append(value)
            if any(np.isfinite(index_values)):
                return index_values, "index"
        return [np.nan] * len(names), ""

    def _show_s_values_and_device_lists(self, list_widget, names):
        names = list(names)
        order_values, order_kind = self._get_element_order_values(names)
        list_widget.clear()
        for name, order_value in zip(names, order_values):
            item = QListWidgetItem(str(name))
            item.setData(Qt.ItemDataRole.UserRole, str(name))
            list_widget.addItem(item)
            try:
                order_value = float(order_value)
            except (ValueError, TypeError):
                order_value = np.nan

            if np.isfinite(order_value):
                if order_kind == "index":
                    order_text = f"index = {int(order_value)}"
                else:
                    order_text = f"S = {order_value:.3f} m"
            else:
                order_text = ""
            item.setData(SPositionDelegate.S_ROLE, order_text)

    def _add_missing_screens_to_list(self, screens):
        existing = []
        for i in range(self.screens_list.count()):
            item = self.screens_list.item(i)
            existing.append(str(item.data(Qt.ItemDataRole.UserRole) or item.text()))
        missing = [str(screen) for screen in screens if self.interface.match_screen_name(screen, existing) is None]
        if not missing:
            return
        combined = existing + missing
        order_values, _ = self._get_element_order_values(combined)
        pairs = sorted(zip(combined, order_values), key=lambda pair: pair[1] if np.isfinite(pair[1]) else np.inf)
        self.screens_list.blockSignals(True)
        self._show_s_values_and_device_lists(self.screens_list, [name for name, _ in pairs])
        self.screens_list.blockSignals(False)

    def _load_logo(self):
        self.logo_label.setText("")
        self.logo_label.setScaledContents(False)

        transform_mode = (Qt.TransformationMode.SmoothTransformation
                          if pyqt_version == 6
                          else Qt.SmoothTransformation
                          )
        base_dir = os.path.dirname(os.path.abspath(__file__))
        logo_path = os.path.join(base_dir, "UI files", "Assets", "CERN_logo.png")
        if not os.path.isfile(logo_path):
            return
        pixmap = QPixmap(logo_path)
        if pixmap.isNull():
            return
        scaled = pixmap.scaledToHeight(80, transform_mode)
        self.logo_label.setPixmap(scaled)
        self.logo_label.setToolTip(logo_path)

    def _set_progress(self, value):
        self.progressBar.setRange(0, 100)
        self.progressBar.setValue(int(max(0, min(100, value))))
        QApplication.processEvents()

    def _clear_fit_panel(self):
        self.result_quad.setText("-")
        self.result_quad_strength.setText("-")
        self.result_emit_x_norm.setText("-")
        self.result_emit_y_norm.setText("-")
        self.result_emit_x_geom.setText("-")
        self.result_emit_y_geom.setText("-")
        self.result_beta_x0.setText("-")
        self.result_alpha_x0.setText("-")
        self.result_beta_y0.setText("-")
        self.result_alpha_y0.setText("-")
        self.result_quad_dx0.setText("-")
        self.result_quad_dy0.setText("-")
        self.result_quad_roll.setText("-")
        self.result_energy_pref.setText("-")
        self.result_reference_screen.setText("-")

    def _interface_quad_value_unit(self):
        cached_unit = getattr(self, "_cached_quad_value_unit", None)
        if cached_unit is not None:
            return cached_unit
        units = (self._get_interface_initial_settings() or {}).get("units", {})
        unit = str(units.get("quadrupole_strength", ""))
        if not unit:
            try:
                quadrupoles = list(getattr(self.interface, "quadrupoles", []))
                unit = str(self.interface.get_quadrupoles(quadrupoles[:1]).get("value_unit", "1/m"))
            except Exception:
                unit = "1/m"
        self._cached_quad_value_unit = unit
        return unit

    def _quadrupole_value_unit(self, session=None):
        session = self.session if session is None else session
        if isinstance(session, dict) and session.get("mode") == "multi_quad_scan":
            session = self._get_session_for_selected_quad(session) or session
        if not isinstance(session, dict):
            return self._interface_quad_value_unit()
        return str(session.get("quad_value_unit", "1/m"))

    def _quadrupole_scan_axis_label(self, session=None):
        unit = self._quadrupole_value_unit(session)
        if unit == "A":
            return "Quadrupole current [A]"
        return f"K1L [{unit}]"

    def _update_fit_panel(self, result):
        self.result_quad.setText(str(result["quad_name"]))
        self.result_reference_screen.setText(result["screen0"])

        def fmt_value(value, suffix=""):  # formats numbers to text
            try:
                value = float(value)
            except Exception:
                return "-"
            if not np.isfinite(value):
                return "-"
            return f"{value:.3f}{suffix}"

        def formatted_result(value, error, unit=""):
            value_text = fmt_value(value)
            error_text = fmt_value(error)
            if value_text == "-":
                return "-"
            if error_text == "-":
                return f"{value_text} {unit}".rstrip()
            return f"{value_text} ± {error_text} {unit}".rstrip()

        quad_strength_text = fmt_value(
            result.get("quad_k1l_0"), f" {self._quadrupole_value_unit()}"
        )
        if result.get("fit_quadrupole_strength", False) and quad_strength_text != "-":
            quad_strength_text += " (fit)"
        elif quad_strength_text != "-":
            quad_strength_text += " (nominal)"
        self.result_quad_strength.setText(quad_strength_text)
        self.result_emit_x_norm.setText(
            formatted_result(result.get("emit_x_norm"), result.get("emit_x_norm_err"), "mm·mrad"))
        self.result_emit_y_norm.setText(
            formatted_result(result.get("emit_y_norm"), result.get("emit_y_norm_err"), "mm·mrad"))
        self.result_emit_x_geom.setText(
            formatted_result(result.get("emit_x_geom"), result.get("emit_x_geom_err"), "nm·rad"))
        self.result_emit_y_geom.setText(
            formatted_result(result.get("emit_y_geom"), result.get("emit_y_geom_err"), "nm·rad"))
        self.result_beta_x0.setText(formatted_result(result.get("beta_x0"), result.get("beta_x0_err"), "m"))
        self.result_alpha_x0.setText(formatted_result(result.get("alpha_x0"), result.get("alpha_x0_err")))
        self.result_beta_y0.setText(formatted_result(result.get("beta_y0"), result.get("beta_y0_err"), "m"))
        self.result_alpha_y0.setText(formatted_result(result.get("alpha_y0"), result.get("alpha_y0_err")))
        self.result_quad_dx0.setText(formatted_result(result.get("quad_dx0"), result.get("quad_dx0_err"), "mm") if result.get("fit_quad_offset") else "-")
        self.result_quad_dy0.setText(formatted_result(result.get("quad_dy0"), result.get("quad_dy0_err"), "mm") if result.get("fit_quad_offset") else "-")
        self.result_quad_roll.setText(formatted_result(result.get("quad_roll"), result.get("quad_roll_err"), "mrad") if result.get("fit_quad_roll") else "-")
        self.result_energy_pref.setText(formatted_result(result.get("energy_pref"), result.get("energy_pref_err"), "MeV/c") if result.get("fit_energy_pref") else "-")
        self.result_reference_screen.setText(result["screen0"])

        print("Errors of the fit:")
        print("Error of emit_x_norm: ", result.get("emit_x_norm_err"))
        print("Error of emit_y_norm: ", result.get("emit_y_norm_err"))
        print("Error of emit_x_geom: ", result.get("emit_x_geom_err"))
        print("Error of emit_y_geom: ", result.get("emit_y_geom_err"))
        print("Error of beta_x0: ", result.get("beta_x0_err"))
        print("Error of alpha_x0: ", result.get("alpha_x0_err"))
        print("Error of beta_y0: ", result.get("beta_y0_err"))
        print("Error of alpha_y0: ", result.get("alpha_y0_err"))

    def _reset_canvas(self):
        fig = self.canvas.figure
        fig.clear()
        ax = fig.add_subplot(111)
        ax.set_title("Quadrupole scan")
        ax.set_xlabel(self._quadrupole_scan_axis_label())
        ax.set_ylabel("Beam size")
        ax.grid(True, alpha=0.3)
        self.canvas.draw()

    def _get_selection(self):
        quadrupoles_all = []
        for i in range(self.quadrupoles_list.count()):
            it = self.quadrupoles_list.item(i)
            quadrupoles_all.append(it.data(Qt.ItemDataRole.UserRole) or it.text())

        screens_all = []
        for i in range(self.screens_list.count()):
            it = self.screens_list.item(i)
            screens_all.append(it.data(Qt.ItemDataRole.UserRole) or it.text())

        selected_quadrupoles = []
        for i in range(self.quadrupoles_list.count()):
            it = self.quadrupoles_list.item(i)
            if it.isSelected():
                selected_quadrupoles.append(it.data(Qt.ItemDataRole.UserRole) or it.text())

        selected_screens = []
        for i in range(self.screens_list.count()):
            it = self.screens_list.item(i)
            if it.isSelected():
                selected_screens.append(it.data(Qt.ItemDataRole.UserRole) or it.text())

        quadrupoles = selected_quadrupoles or quadrupoles_all
        screens = selected_screens or screens_all

        return quadrupoles, screens

    def _draw_live_scan(self, session):
        if session is None:
            return
        self._refresh_plot_comboboxes_from_session(session)
        session_to_plot = self._get_session_for_selected_quad(self._session_without_excluded_points(session))
        if session_to_plot is None:
            return
        K1L_values = np.asarray(session_to_plot["K1L_values"], dtype=float)
        sigx = np.asarray(session_to_plot["sigx_mean"], dtype=float)
        sigy = np.asarray(session_to_plot["sigy_mean"], dtype=float)
        screens = list(session_to_plot["screens"])
        quad_name = session_to_plot.get("quad_name", "-")
        em_sigma_unit = session_to_plot.get("sigma_unit", self._get_interface_units())
        quad_axis_label = self._quadrupole_scan_axis_label(session_to_plot)
        fig = self.canvas.figure
        fig.clear()

        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212, sharex=ax1)

        color_cycle = [
            "#E69F00",  # orange
            "#0072B2",  # blue
            "#009E73",  # green
            "#D55E00",  # red
            "#CC79A7",  # reddish purple
            "#56B4E9",  # light blue
            "#F0E442",  # yellow
            "#000000",  # black
        ]

        if self.show_scan_on_all_screens.isChecked():
            screen_indices = list(range(len(screens)))
        else:
            selected_screen = self.screen_on_plot.currentText().strip()
            if selected_screen in screens:
                screen_indices = [screens.index(selected_screen)]
            else:
                screen_indices = list(range(len(screens)))

        for i in screen_indices:
            screen = screens[i]
            mask_x = np.isfinite(sigx[:, i])
            mask_y = np.isfinite(sigy[:, i])
            color = color_cycle[i % len(color_cycle)]
            ax1.plot(K1L_values[mask_x], sigx[mask_x, i], 'o--', color=color, label=screen)
            ax2.plot(K1L_values[mask_y], sigy[mask_y, i], 'o--', color=color, label=screen)
        ax1.set_title(f"Quadrupole scan: {quad_name}")
        ax1.set_ylabel(f"sigx [{em_sigma_unit}]")
        ax2.set_ylabel(f"sigy [{em_sigma_unit}]")
        ax2.set_xlabel(quad_axis_label)

        ax1.grid(True, alpha=0.3)
        ax2.grid(True, alpha=0.3)
        ax1.legend(fontsize=8, ncol=2)
        ax2.legend(fontsize=8, ncol=2)

        fig.tight_layout()
        self.canvas.draw()

    def _plot_fit_overlay(self, pred_x, pred_y, result=None, screens=None, fit_k1l_values=None):
        if self.session is None:
            return
        plotted_session = self._session_without_excluded_points(self.session)
        session_screens = list(self.session.get("screens", []))
        if screens is None:
            screens = session_screens
        else:
            screens = [screen for screen in screens if screen in session_screens]
        pred_x = np.asarray(pred_x, dtype=float)
        pred_y = np.asarray(pred_y, dtype=float)
        prediction_observable = str((result or {}).get("prediction_observable", "sigma"))
        if prediction_observable not in {"sigma", "sigma2"}:
            raise ValueError(f"Unknown prediction observable: {prediction_observable}")
        n_screens = min(len(screens), pred_x.shape[1], pred_y.shape[1])
        screens = screens[:n_screens]
        K1L_values = np.asarray(plotted_session["K1L_values"], dtype=float)
        sigx = np.asarray(plotted_session["sigx_mean"], dtype=float)
        sigy = np.asarray(plotted_session["sigy_mean"], dtype=float)
        if fit_k1l_values is not None and len(fit_k1l_values) == pred_x.shape[0]:
            fit_K1L_values = np.asarray(fit_k1l_values, dtype=float)
        else:
            fit_K1L_values = K1L_values
        fig = self.canvas.figure
        fig.clear()

        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212, sharex=ax1)

        color_cycle = [
            "#E69F00",  # orange
            "#0072B2",  # blue
            "#009E73",  # green
            "#D55E00",  # red
            "#CC79A7",  # reddish purple
            "#56B4E9",  # light blue
            "#F0E442",  # yellow
            "#000000",  # black
        ]

        for prediction_i, screen in enumerate(screens):
            session_i = session_screens.index(screen)
            color = color_cycle[session_i % len(color_cycle)]
            ax1.plot(K1L_values, sigx[:, session_i], "o--", color=color, linewidth=1.0, label=f"{screen} data")
            fit_x = pred_x[:, prediction_i] if prediction_observable == "sigma" else np.sqrt(
                np.maximum(pred_x[:, prediction_i], 0.0))
            ax1.plot(fit_K1L_values, fit_x, "-", color=color, linewidth=2.0, label=f"{screen} fit")
            ax2.plot(K1L_values, sigy[:, session_i], "o--", color=color, linewidth=1.0, label=f"{screen} data")
            fit_y = pred_y[:, prediction_i] if prediction_observable == "sigma" else np.sqrt(
                np.maximum(pred_y[:, prediction_i], 0.0))
            ax2.plot(fit_K1L_values, fit_y, "-", color=color, linewidth=2.0, label=f"{screen} fit")

        unit = self.session.get("sigma_unit", self._get_interface_units())
        quad_axis_label = self._quadrupole_scan_axis_label(self.session)

        ax1.set_title(f"Quadrupole scan: {self.session.get('quad_name', '-')}")
        ax1.set_ylabel(f"sigx [{unit}]")
        ax2.set_ylabel(f"sigy [{unit}]")
        ax2.set_xlabel(quad_axis_label)

        ax1.grid(True, alpha=0.3)
        ax2.grid(True, alpha=0.3)
        ax1.legend(fontsize=7, ncol=2)
        ax2.legend(fontsize=7, ncol=2)

        fig.tight_layout()
        self.canvas.draw()

    def _get_session_data_from_database(self):
        self._excluded_points = set()
        self._update_delete_point_button()
        states = list(getattr(self, "loaded_states_from_scan", []))
        files = list(getattr(self, "loaded_state_files", []))
        if not states:
            return
        folder = self.load_screens_data_database.text().strip()

        is_quad_scan = bool(self.emittance_settings.get("is_quad_scan", True))
        steps_requested = int(self.emittance_settings["scan_steps"])
        quad_name = self.emittance_settings.get("quad_name")
        if not quad_name:
            raise ValueError("Choose a quadrupole before doing a fixed-K1L session.")
        quad_value_unit = "1/m"

        if is_quad_scan:
            current_A_min = float(self.emittance_settings.get("current_A_min", 0.0))
            current_A_max = float(self.emittance_settings.get("current_A_max", 0.0))
            nsteps_scan = max(int(os.path.basename(path).replace(".pkl", "").split("_")[3]) for path in files) + 1
            K1L_values = np.full(nsteps_scan, np.nan)
            for path, state in zip(files, states):
                filename = os.path.basename(path)
                step_i = int(filename.split("_")[3])  # screen_0000_step_0003_shot_0000.pkl -> 0003
                quad = state.get_quadrupoles()
                quad_value_unit = str(quad.get("value_unit", quad_value_unit))
                names = list(quad.get("names", []))
                bdes = np.ravel(np.asarray(quad.get("bdes", []), dtype=float))
                index = names.index(quad_name) if quad_name in names else 0
                if index < bdes.size:
                    K1L_values[step_i] = float(bdes[index])
            K1L_0 = float(self.emittance_settings.get("K1L_0", np.nan))
            if not np.isfinite(K1L_0):
                K1L_0 = float(np.nanmedian(K1L_values))
            if np.isfinite(K1L_0) and not np.isclose(K1L_0, 0.0):
                deltas = K1L_values / K1L_0 - 1.0
            else:
                deltas = np.full(nsteps_scan, np.nan)

        else:
            strengths = []
            for state in states:
                quadrupoles = state.get_quadrupoles()
                quad_value_unit = str(quadrupoles.get("value_unit", quad_value_unit))
                names = list(quadrupoles.get("names", []))
                bdes = np.asarray(quadrupoles.get("bdes", []), dtype=float)
                if quad_name in names:
                    strengths.append(float(bdes[names.index(quad_name)]))
            if not strengths or not np.any(np.isfinite(strengths)):
                raise ValueError("The loaded fixed-K1L session does not contain the selected quadrupole strength.")
            K1L_0 = float(np.nanmean(strengths))
            current_A_min, current_A_max, nsteps_scan = 0.0, 0.0, 1
            deltas = np.array([0.0])
            K1L_values = np.array([K1L_0])

        nscreens = int(self.emittance_settings["nscreens"])
        screens = list(self.emittance_settings.get("screens", []))
        if not screens:
            _, screens = self._get_selection()
        screens = screens[:nscreens]

        nshots = int(self.emittance_settings["nshots"])
        sigx_samples = np.full((nsteps_scan, nscreens, nshots), np.nan)
        sigy_samples = np.full((nsteps_scan, nscreens, nshots), np.nan)
        sigxy_samples = np.full((nsteps_scan, nscreens, nshots), np.nan)
        x_samples = np.full((nsteps_scan, nscreens, nshots), np.nan)
        y_samples = np.full((nsteps_scan, nscreens, nshots), np.nan)
        images = [[[None for _ in range(nshots)] for _ in range(nscreens)] for _ in range(nsteps_scan)]
        hedges = [[[None for _ in range(nshots)] for _ in range(nscreens)] for _ in range(nsteps_scan)]
        vedges = [[[None for _ in range(nshots)] for _ in range(nscreens)] for _ in range(nsteps_scan)]
        print(f"GUI Nshots: {nshots}, GUI Scan steps: {nsteps_scan}")

        for path, state in zip(files, states):
            filename = os.path.basename(path)
            parts = filename.replace(".pkl", "").split("_")
            screen_i = int(parts[1])
            step_i = int(parts[3])
            shot_i = int(parts[5])
            screen_data = state.get_screens()
            sigx_samples[step_i, screen_i, shot_i] = float(np.ravel(screen_data["sigx"])[0])  # / 1000.0
            sigy_samples[step_i, screen_i, shot_i] = float(np.ravel(screen_data["sigy"])[0])  # / 1000.0
            sigxy_samples[step_i, screen_i, shot_i] = float(np.ravel(screen_data.get("sigxy", [np.nan]))[0])  # / 1000.0
            x_samples[step_i, screen_i, shot_i] = float(np.ravel(screen_data.get("x", [np.nan]))[0])
            y_samples[step_i, screen_i, shot_i] = float(np.ravel(screen_data.get("y", [np.nan]))[0])
            screen_images = state.get_screens().get("images", [])
            if len(screen_images) > 0:
                images[step_i][screen_i][shot_i] = np.asarray(screen_images[0])
            screen_hedges = screen_data.get("hedges", [])
            if len(screen_hedges) > 0:
                hedges[step_i][screen_i][shot_i] = np.asarray(screen_hedges[0], dtype=float)
            screen_vedges = screen_data.get("vedges", [])
            if len(screen_vedges) > 0:
                vedges[step_i][screen_i][shot_i] = np.asarray(screen_vedges[0], dtype=float)

        sigx_mean = np.nanmedian(sigx_samples, axis=2)
        sigy_mean = np.nanmedian(sigy_samples, axis=2)
        sigxy_mean = np.nanmean(sigxy_samples, axis=2)
        sigx_std = np.nanstd(sigx_samples, axis=2)
        sigy_std = np.nanstd(sigy_samples, axis=2)
        sigxy_std = np.nanstd(sigxy_samples, axis=2)
        x_mean = np.nanmean(x_samples, axis=2)
        y_mean = np.nanmean(y_samples, axis=2)
        x_std = np.nanstd(x_samples, axis=2)
        y_std = np.nanstd(y_samples, axis=2)

        scan_steps = []
        for i in range(nsteps_scan):
            state_files = [path for path in files if int(os.path.basename(path).split("_")[3]) == i]

            scan_steps.append({
                "step_index": int(i),
                "delta": float(deltas[i]),
                "K1L": float(K1L_values[i]),
                "quad_value": float(K1L_values[i]),
                "state_files": state_files,
            })

        session = {
            "current_A_min": current_A_min,
            "current_A_max": current_A_max,
            "is_quad_scan": is_quad_scan,
            "steps": steps_requested,
            "nshots": int(self.emittance_settings["nshots"]),
            "sigma_unit": "mm",
            "quad_name": quad_name,
            "quadrupoles": [quad_name],
            "screens": screens,
            "reference_screen": screens[0] if screens else "",
            "quad_value_unit": quad_value_unit,
            "K1L_0": float(K1L_0),
            "sigx_mean": sigx_mean.tolist(),
            "sigy_mean": sigy_mean.tolist(),
            "sigxy_mean": sigxy_mean.tolist(),
            "sigx_std": sigx_std.tolist(),
            "sigy_std": sigy_std.tolist(),
            "sigx_shots": sigx_samples.tolist(),
            "sigy_shots": sigy_samples.tolist(),
            "sigxy_shots": sigxy_samples.tolist(),
            "sigxy_std": sigxy_std.tolist(),
            "x_mean": x_mean.tolist(),
            "y_mean": y_mean.tolist(),
            "x_std": x_std.tolist(),
            "y_std": y_std.tolist(),
            "x_shots": x_samples.tolist(),
            "y_shots": y_samples.tolist(),
            "deltas": deltas.tolist(),
            "K1L_values": K1L_values.tolist(),
            "current_values": [self.quad_setpoint_to_current(quad_name, value, quad_value_unit) for value in
                               K1L_values],
            "scan_steps": scan_steps,
            "states_dir": folder,
            "cancelled": False,
            "nsteps_scan": int(nsteps_scan),
            "images": images,
            "hedges": hedges,
            "vedges": vedges,
        }

        print(f"Quadrupole values [{quad_value_unit}]:", session["K1L_values"])
        print("sigx:", session["sigx_mean"])
        print("sigy:", session["sigy_mean"])
        print("unit:", session.get("sigma_unit"))

        return session

    def _run_optimization(self):
        self.log("Fitting emittance and twiss parameters at scanned quadrupole started...")
        session_was_missing = self.session is None
        if self.session is None:
            data_folder = self.load_screens_data_database.text().strip()
            if data_folder and os.path.isdir(data_folder):
                self.session = self._get_session_data_from_database()
                self._refresh_plot_comboboxes_from_session(self.session)
                self._draw_live_scan(self.session)
            if self.session is None:
                QMessageBox.information(self, "Optimization", "No session.")
                return
        if session_was_missing:
            self._set_default_quad_strength_bounds_from_session(self.session)
        bounds = self._get_bounds_from_gui()
        if self._is_optimizing:
            return
        self._is_optimizing = True
        self._current_optimizer = None
        self._set_progress(0)
        self._optimization_t0 = time.perf_counter()
        thread = QThread(self)
        computing_method = self.computing_method_combo.currentText().strip()
        _, selected_screens = self._get_selection()
        worker = OptimizationWorker(self.interface, self._session_without_excluded_points(self.session),
                                    selected_screens=selected_screens, bounds=bounds,
                                    fit_quadrupole_strength=bool(self.fit_quadrupole_strength_checkbox.isChecked()),
                                    fit_quad_offset=bool(self.fit_quad_offset_checkbox.isChecked()),
                                    fit_quad_roll=bool(self.fit_quad_roll_checkbox.isChecked()),
                                    fit_energy_pref=bool(self.fit_energy_checkbox.isChecked()),
                                    computing_method=computing_method)
        worker.info.connect(self.log)

        worker.moveToThread(thread)
        worker.optimizer_ready.connect(self._store_current_optimizer)
        worker.finished.connect(self._on_optimization_output)
        worker.error.connect(self._on_optimization_error)
        worker.progress.connect(self._on_optimization_progress)
        worker.done.connect(thread.quit)
        worker.done.connect(worker.deleteLater)

        thread.finished.connect(self._on_optimization_finished)
        thread.finished.connect(thread.deleteLater)
        thread.started.connect(worker.run)
        self._optimization_thread = thread
        self._optimization_worker = worker
        worker.ml_not_found_rft_fallback.connect(
            lambda: self.computing_method_combo.setCurrentText(ComputationMode.RFT.value))
        self._set_progress(30)
        thread.start()

    def _on_optimization_progress(self, phase, current, total):
        total = max(int(total), 1)
        current = max(0, min(int(current), total))
        value = 30 + 65 * current / total
        self._set_progress(value)
        self.progressBar.setFormat(f"{phase}: {current}/{total}")

    def _store_current_optimizer(self, optimizer):
        self._current_optimizer = optimizer

    def _on_optimization_output(self, output):
        self._set_progress(85)
        result = output["result"]
        pred_x = np.asarray(output["pred_x"], dtype=float)
        pred_y = np.asarray(output["pred_y"], dtype=float)
        optimization_screens = list(output.get("screens", self.session.get("screens", [])))
        fit_k1l_values = output.get("K1L_values")
        self.session["optimization_result"] = result
        self.session["optimization_pred_x"] = pred_x.tolist()
        self.session["optimization_pred_y"] = pred_y.tolist()
        self._update_fit_panel(result)
        self._plot_fit_overlay(pred_x, pred_y, result, screens=optimization_screens, fit_k1l_values=fit_k1l_values)
        self.save_emittance_measurement_session(session=self.session, is_fit_quad_strength_checked=bool(
            self.fit_quadrupole_strength_checkbox.isChecked()), bounds=self._get_bounds_from_gui(),
                                                target_dir=getattr(self, "_general_session_dir", None))
        self._set_progress(100)

        elapsed = time.perf_counter() - self._optimization_t0

        joint_found = bool(
            np.isfinite(result.get("emit_x_norm", np.nan)) and np.isfinite(result.get("emit_y_norm", np.nan)))
        paused = bool(result.get("paused", False))

        if paused:
            if joint_found:
                message = (
                    "Best joint solution found so far.\n\n"
                    f"εₓ = {result['emit_x_norm']:.4f} mm·mrad\n"
                    f"εᵧ = {result['emit_y_norm']:.4f} mm·mrad\n"
                    f"βₓ0 = {result['beta_x0']:.4f} m, αₓ0 = {result['alpha_x0']:.4f}\n"
                    f"βᵧ0 = {result['beta_y0']:.4f} m, αᵧ0 = {result['alpha_y0']:.4f}"
                )
            else:
                message = "Optimization was paused before any joint solution was found."

            self._optimization_paused = True
            QMessageBox.information(self, "Optimization paused", message)

        elif result.get("stopped", False):
            self._optimization_paused = False

            if joint_found:
                message = (
                    "Best joint solution found so far.\n\n"
                    f"εₓ = {result['emit_x_norm']:.4f} mm·mrad\n"
                    f"εᵧ = {result['emit_y_norm']:.4f} mm·mrad\n"
                    f"βₓ0 = {result['beta_x0']:.4f} m, αₓ0 = {result['alpha_x0']:.4f}\n"
                    f"βᵧ0 = {result['beta_y0']:.4f} m, αᵧ0 = {result['alpha_y0']:.4f}"
                )
            else:
                message = "Optimization was stopped before any joint solution was found."

            QMessageBox.information(self, "Optimization stopped", message)

        else:
            self._optimization_paused = False
            QMessageBox.information(
                self,
                "Optimization complete",
                f"εₓ = {result['emit_x_norm']:.4f} mm·mrad\n"
                f"εᵧ = {result['emit_y_norm']:.4f} mm·mrad\n"
                f"βₓ0 = {result['beta_x0']:.4f} m, αₓ0 = {result['alpha_x0']:.4f}\n"
                f"βᵧ0 = {result['beta_y0']:.4f} m, αᵧ0 = {result['alpha_y0']:.4f}"
            )

        print(f"Elapsed time: {elapsed}s = {elapsed / 60}min")

    def _on_optimization_error(self, message):
        self._set_progress(0)
        self._optimization_paused = False
        if message == "Optimization stopped.":
            QMessageBox.information(self, "Optimization stopped",
                                    "Optimization was stopped before any solution was found.")
            self.log("Optimization was stopped before any solution was found.")
        else:
            QMessageBox.information(self, "Optimization", message)

    def _on_optimization_finished(self):
        self._is_optimizing = False
        self._current_optimizer = None
        self._optimization_worker = None
        self._optimization_thread = None

    def _scan_progress_callback(self, session_partial, current_step, total_steps):  # refreshes plot in the gui
        if self._scan_stop_requested:
            raise KeyboardInterrupt("Scan stopped by user.")
        self.session = session_partial
        self._draw_live_scan(session_partial)
        if total_steps:
            self._set_progress(100.0 * float(current_step) / float(total_steps))
        QApplication.processEvents()
        if self._scan_stop_requested:
            raise KeyboardInterrupt("Scan stopped by user.")

    def _run_scan(self, resume=None):
        resume_states_dir = resume["states_dir"] if resume else None
        quadrupoles, _ = self._get_selection()
        if resume:
            quadrupoles = list(resume["quadrupoles"])
        if len(quadrupoles) == 0:
            QMessageBox.information(self, "Scan error", "No quadrupole selected.")
            return
        quad_label = quadrupoles[0] if len(quadrupoles) == 1 else f"multi-quad scan ({len(quadrupoles)} quadrupoles)"
        steps_preview = int(self.steps_settings.value())
        if steps_preview == 0:
            self.log(f"Gathering screen data...")
        else:
            self.log(f"Running quadrupole scan for {quad_label}...")
        self.quad_on_plot.blockSignals(True)
        self.quad_on_plot.clear()
        self.quad_on_plot.addItems(quadrupoles)
        if quadrupoles:
            self.quad_on_plot.setCurrentIndex(0)
        self.quad_on_plot.blockSignals(False)
        selected_items = self.screens_list.selectedItems()
        if not selected_items:
            self.screens_list.blockSignals(True)
            self.screens_list.selectAll()
            self.screens_list.blockSignals(False)
        _, screens = self._get_selection()
        if resume:
            screens = list(resume["screens"])
        self.screen_on_plot.blockSignals(True)
        self.screen_on_plot.clear()
        self.screen_on_plot.addItems(screens)
        if screens:
            self.screen_on_plot.setCurrentIndex(0)
        self.screen_on_plot.blockSignals(False)
        if not screens:
            QMessageBox.information(self, "Scan error", "No screens available.")
            return

        current_min = float(self.minimum_current.value())
        current_max = float(self.maximum_current.value())
        steps = int(self.steps_settings.value())
        nshots = int(self.meas_per_step.value())
        screen_current_ranges = {screen: values for screen, values in
                                 getattr(self, "_screen_current_ranges", {}).items() if screen in screens}
        if resume:
            current_min = float(resume["current_min"])
            current_max = float(resume["current_max"])
            steps = int(resume["steps"])
            nshots = int(resume["nshots"])
            screen_current_ranges = dict(resume["screen_current_ranges"])
        if steps > 0 and current_max <= current_min:
            QMessageBox.information(self, "Scan error",
                                    "Maximum scan current must be larger than the minimum scan current.")
            return

        self._last_scan_status = None
        self._scan_stop_requested = False
        self._is_scanning = True
        if resume is None:
            self._excluded_points = set()
            self._update_delete_point_button()
        scan_settings = {"quadrupoles": quadrupoles, "screens": screens, "current_min": current_min,
                         "current_max": current_max, "steps": steps, "nshots": nshots,
                         "screen_current_ranges": screen_current_ranges}

        self._clear_fit_panel()
        self._set_progress(0)
        try:
            self.session = self.run_scan(quad_name=quadrupoles, current_min=current_min, current_max=current_max,
                                         steps=steps, nshots=nshots,
                                         screens=screens, screen_current_ranges=screen_current_ranges,
                                         reference_screen=screens[0],
                                         progress_callback=self._scan_progress_callback,
                                         resume_states_dir=resume_states_dir)
            self._interrupted_scan = None
            self._set_default_quad_strength_bounds_from_session(self.session)
            self._update_quad_readback_label()
            if steps == 0:
                self.log("Finished gathering data from the screens.")
            else:
                self.log("Quadrupole scan finished.")
        except KeyboardInterrupt as e:
            self._set_progress(0)
            self._remember_interrupted_scan(scan_settings, str(e))
            QMessageBox.information(self, "Scan", str(e))
            return
        except TypeError as e:
            self._set_progress(0)
            self._remember_interrupted_scan(scan_settings, f"Type Error: {e}")
            QMessageBox.information(self, "Scan error", f"Type Error: {e}")
            return
        except Exception as e:
            self._set_progress(0)
            self._remember_interrupted_scan(scan_settings, str(e))
            QMessageBox.information(self, "Scan error",
                                    f"{e}\n\nPress Resume to continue this scan from where it stopped."
                                    if getattr(self, "_interrupted_scan", None) else str(e))
            return
        finally:
            self._is_scanning = False
        QMessageBox.information(self, "Scan", f"Scan completed.")
        self._set_progress(100)

    def _remember_interrupted_scan(self, scan_settings, reason):
        states_dir = getattr(self, "_last_scan_states_dir", None)
        if not states_dir or not os.path.isdir(states_dir):
            self._interrupted_scan = None
            return
        if len(list(scan_settings.get("quadrupoles", []))) != 1 or int(scan_settings.get("steps", 0)) <= 0:
            self._interrupted_scan = None
            return
        self._interrupted_scan = dict(scan_settings, states_dir=states_dir, reason=str(reason))
        self.log(
            f"Scan interrupted ({reason}). Measured points are kept in {states_dir}. You can click resume to continue.")

    def _resume_interrupted_scan(self):
        interrupted = getattr(self, "_interrupted_scan", None)
        if not interrupted:
            return False
        self.log(f"Resuming the scan of {interrupted['quadrupoles'][0]} from {interrupted['states_dir']}...")
        self._scan_stop_requested = False
        self._scan_pause_requested = False
        self._scan_is_paused = False
        self._run_scan(resume=interrupted)
        return True

    def _get_twiss_s_positions(self, names):
        names = list(names)
        positions = [np.nan] * len(names)

        if not hasattr(self.interface, "_get_elements_positions"):
            return positions

        try:
            pos = self.interface._get_elements_positions()
            pos_names = [str(name) for name in pos.get("names", [])]
            s_values = np.asarray(pos.get("S", []), dtype=float)
            lookup = {name: float(s_values[i]) for i, name in enumerate(pos_names) if
                      i < s_values.size and np.isfinite(s_values[i])}
            positions = []
            for requested_name in names:
                requested_name = str(requested_name)
                matched_name = self.interface.match_screen_name(requested_name, lookup)
                if matched_name is not None:
                    positions.append(lookup[matched_name])
                    continue
                quad_part_positions = []
                for lattice_name in (f"{requested_name}_1", f"{requested_name}_2"):
                    if lattice_name in lookup:
                        quad_part_positions.append(lookup[lattice_name])
                if quad_part_positions:
                    positions.append(min(quad_part_positions))
                else:
                    positions.append(np.nan)
        except Exception:
            positions = [np.nan] * len(names)

        return positions

    def _get_session_for_selected_quad(self, session):
        if not isinstance(session, dict):
            return session
        if session.get("mode") != "multi_quad_scan":
            return session
        per_quad = list(session.get("per_quad_sessions", []))
        if len(per_quad) == 0:
            return None
        combo_index = int(self.quad_on_plot.currentIndex())
        if 0 <= combo_index < len(per_quad):
            return per_quad[combo_index]
        selected_quad = self.quad_on_plot.currentText().strip()
        if selected_quad:
            for quad_session in per_quad:
                if str(quad_session.get("quad_name", "")).strip() == selected_quad:
                    return quad_session
        return per_quad[0]

    def _refresh_plot_comboboxes_from_session(self, session):
        if not isinstance(session, dict):
            return
        if session.get("mode") == "multi_quad_scan":
            per_quad = list(session.get("per_quad_sessions", []))
            quad_names = [str(qs.get("quad_name", "")).strip() for qs in per_quad if isinstance(qs, dict)]
            if quad_names:
                current_quad = self.quad_on_plot.currentText().strip()
                self.quad_on_plot.blockSignals(True)
                self.quad_on_plot.clear()
                self.quad_on_plot.addItems(quad_names)
                if current_quad in quad_names:
                    self.quad_on_plot.setCurrentText(current_quad)
                else:
                    self.quad_on_plot.setCurrentIndex(0)
                self.quad_on_plot.blockSignals(False)

            session_to_plot = self._get_session_for_selected_quad(session)
            if isinstance(session_to_plot, dict):
                screens = list(session_to_plot.get("screens", []))
                current_screen = self.screen_on_plot.currentText().strip()
                self.screen_on_plot.blockSignals(True)
                self.screen_on_plot.clear()
                self.screen_on_plot.addItems(screens)
                if current_screen in screens:
                    self.screen_on_plot.setCurrentText(current_screen)
                elif screens:
                    self.screen_on_plot.setCurrentIndex(0)
                self.screen_on_plot.blockSignals(False)

        else:
            quad_name = str(session.get("quad_name", "")).strip()
            screens = list(session.get("screens", []))

            if quad_name:
                self.quad_on_plot.blockSignals(True)
                self.quad_on_plot.clear()
                self.quad_on_plot.addItem(quad_name)
                self.quad_on_plot.setCurrentIndex(0)
                self.quad_on_plot.blockSignals(False)

            current_screen = self.screen_on_plot.currentText().strip()
            self.screen_on_plot.blockSignals(True)
            self.screen_on_plot.clear()
            self.screen_on_plot.addItems(screens)
            if current_screen in screens:
                self.screen_on_plot.setCurrentText(current_screen)
            elif screens:
                self.screen_on_plot.setCurrentIndex(0)
            self.screen_on_plot.blockSignals(False)

        self.screen_on_plot.setEnabled(not self.show_scan_on_all_screens.isChecked())

    def _filter_quadrupoles_in_gui(self):
        previously_selected = []
        for i in range(self.quadrupoles_list.count()):
            item = self.quadrupoles_list.item(i)
            if item.isSelected():
                previously_selected.append(item.data(Qt.ItemDataRole.UserRole) or item.text())

        if previously_selected:
            self._last_selected_quadrupoles = list(previously_selected)
        else:
            self._last_selected_quadrupoles = list(getattr(self, "_last_selected_quadrupoles", []))

        _, selected_screens = self._get_selection()
        if not selected_screens:
            return

        screen_position, screen_order_kind = self._get_element_order_values(selected_screens)
        finite_screen_positions = [float(s) for s in screen_position if np.isfinite(s)]
        if not finite_screen_positions:
            return

        first_screen_position = min(finite_screen_positions)
        all_quadrupoles = list(getattr(self.interface, "quadrupoles", []))
        quad_order, quad_order_kind = self._get_element_order_values(all_quadrupoles)

        if quad_order_kind != screen_order_kind:  # S [m] or index
            return

        quad_pos = {name: float(s) for name, s in zip(all_quadrupoles, quad_order) if np.isfinite(s)}
        valid_quadrupoles = [name for name in all_quadrupoles if
                             name in quad_pos and quad_pos[name] < first_screen_position]
        valid_previous = [q for q in getattr(self, "_last_selected_quadrupoles", []) if q in valid_quadrupoles]

        self.quadrupoles_list.blockSignals(True)
        self._show_s_values_and_device_lists(self.quadrupoles_list, valid_quadrupoles)
        if valid_previous:
            for i in range(self.quadrupoles_list.count()):
                item = self.quadrupoles_list.item(i)
                item_name = item.data(Qt.ItemDataRole.UserRole) or item.text()
                if item_name in valid_previous:
                    item.setSelected(True)
            self._last_selected_quadrupoles = list(valid_previous)
        elif valid_quadrupoles:
            closest_quad = max(valid_quadrupoles, key=lambda name: quad_pos[name])
            for i in range(self.quadrupoles_list.count()):
                item = self.quadrupoles_list.item(i)
                item_name = item.data(Qt.ItemDataRole.UserRole) or item.text()
                if item_name == closest_quad:
                    item.setSelected(True)
                    self._last_selected_quadrupoles = [closest_quad]
                    break
        elif self.quadrupoles_list.count() > 0:
            item = self.quadrupoles_list.item(0)
            if item is not None:
                item.setSelected(True)
                item_name = item.data(Qt.ItemDataRole.UserRole) or item.text()
                self._last_selected_quadrupoles = [item_name]
        else:
            self._last_selected_quadrupoles = []

        self.quadrupoles_list.blockSignals(False)

    def _screen_selection_changed(self):
        self._filter_quadrupoles_in_gui()
        self._update_per_screen_ranges_button()
        self._update_quad_readback_label()

    def _show_console_log(self):
        if self.log_console is None:
            self.log_console = LogConsole(self)
        self.log_console.show()
        self.log_console.raise_()
        self.log_console.activateWindow()

    def _show_phase_spaces(self):
        result = None
        reference_name = None
        if isinstance(self.session, dict):
            result = self.session.get("optimization_result")
            reference_name = self.session.get("quad_name") or self.session.get("current_quadrupole")
        if self.phase_spaces is None:
            pass
            # self.phase_spaces = PhaseSpaces(self)
        screens = []
        session_to_plot = None
        if isinstance(self.session, dict):
            session_to_plot = self._get_session_for_selected_quad(self.session)
            if isinstance(session_to_plot, dict):
                screens = list(session_to_plot.get("screens", []))
        if not isinstance(result, dict):
            QMessageBox.information(self, "Phase Space", "Run the emittance/Twiss optimization first.")
            return
        session_to_plot = self._get_session_for_selected_quad(self.session) if isinstance(self.session, dict) else None
        if isinstance(session_to_plot, dict) and screens:
            self.phase_spaces.plot_projection_constraints(result, session_to_plot, interface=self.interface)
        else:
            self.phase_spaces.plot_from_result(result, reference_name=reference_name)
        self.phase_spaces.show()
        self.phase_spaces.raise_()
        self.phase_spaces.activateWindow()

    def _show_screen_images(self):
        if self.screen_images is None:
            self.screen_images = DisplayScreenImages(self)
        if self.session is None:
            QMessageBox.information(self, "No screen images", "No data to display as screen image.")
            return
        else:
            self.screen_images._plot_screen_image(session=self.session)
        self.screen_images.show()
        self.screen_images.raise_()
        self.screen_images.activateWindow()

    def log(self, text):
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        line = f"[{timestamp}] {text}"
        if self.log_console is None:
            self.log_console = LogConsole(self)
        self.log_console.log(line)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    from Backend import SelectInterface

    dialog = SelectInterface.choose_acc_and_interface()
    if dialog is None:
        print("Selection cancelled.")
        sys.exit(1)

    I = dialog

    project_name = I.get_name()
    is_simulation = bool(getattr(I, "is_simulation"))
    bg_shots = 10

    if is_simulation:
        bg_shots = 0
    print(f"Selected interface: {project_name}")
    time_str = datetime.now().strftime("%Y%m%d_%H%M%S")
    dir_name = f"~/CERN-Flight_Simulator-Data/EM_{I.get_name()}{time_str}_session_settings"
    dir_name = os.path.expanduser(os.path.expandvars(dir_name))

    w = MainWindow(I, dir_name=dir_name, is_simulation=is_simulation, bg_shots=bg_shots)
    w.show()
    sys.exit(app.exec())
