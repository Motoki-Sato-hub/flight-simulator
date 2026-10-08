from datetime import datetime
import argparse, glob, os, re, sys, threading
import matplotlib
import numpy as np
from Backend.Response import Response
from zoneinfo import ZoneInfo
from Interfaces.interface_setup import INTERFACE_SETUP

try:
    from PyQt6 import uic
    from PyQt6.QtCore import Qt, QMetaObject, QObject, QProcess, QProcessEnvironment, QThread, QTimer, pyqtSignal, pyqtSlot
    from PyQt6.QtGui import QPixmap
    from PyQt6.QtWidgets import (QApplication, QCheckBox, QDialog, QDialogButtonBox, QDoubleSpinBox, QFileDialog,
                                 QFormLayout, QHeaderView, QLabel, QMainWindow, QMessageBox, QPushButton, QSizePolicy,
                                 QVBoxLayout, QWidget, QTableWidgetItem)

    pyqt_version = 6
except ImportError:
    from PyQt5 import uic
    from PyQt5.QtCore import Qt, QMetaObject, QObject, QProcess, QProcessEnvironment, QThread, QTimer, pyqtSignal, pyqtSlot
    from PyQt5.QtGui import QPixmap
    from PyQt5.QtWidgets import (QApplication, QCheckBox, QDialog, QDialogButtonBox, QDoubleSpinBox, QFileDialog,
                                 QFormLayout, QHeaderView, QLabel, QMainWindow, QMessageBox, QPushButton, QSizePolicy,
                                 QVBoxLayout, QWidget)

    pyqt_version = 5

matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from Backend.SaveOrLoad import SaveOrLoad
from Backend.State import State
from Backend.ResponseMatrix_DFS_WFS import ResponseMatrix_DFS_WFS


class OrbitPoller(QObject):
    bpms_ready = pyqtSignal(object)

    def __init__(self, interface, lock, interval_ms):
        super().__init__()
        self._interface = interface
        self._lock = lock
        self._interval_ms = interval_ms
        self._timer = None

    @pyqtSlot()
    def start(self):
        self._timer = QTimer(self)
        self._timer.setInterval(self._interval_ms)
        self._timer.timeout.connect(self._poll)
        self._timer.start()

    @pyqtSlot()
    def stop(self):
        if self._timer is not None:
            self._timer.stop()

    def _poll(self):
        if not self._lock.acquire(blocking=False):
            return
        try:
            raw = self._interface.get_bpms()
        except Exception as exc:  # noqa: BLE001
            print(f"Live BPM read skipped: {exc}")
            return
        finally:
            self._lock.release()
        self.bpms_ready.emit(raw)


class TwinPlot(FigureCanvas):
    def __init__(self, parent=None):
        figure = Figure(tight_layout={"rect": (0.0, 0.07, 1.0, 1.0)})
        super().__init__(figure)
        self.setParent(parent)
        self.ax_x = figure.add_subplot(211)
        self.ax_y = figure.add_subplot(212, sharex=self.ax_x)

    def clear(self):
        for axis in (self.ax_x, self.ax_y):
            axis.clear()
            axis.grid(True, alpha=0.3)

    def set_time_caption(self, text):
        if getattr(self, "_time_caption", None) is None:
            self._time_caption = self.figure.text(0.5, 0.015, text, ha="center", va="bottom", fontsize=12, fontweight="bold", color="black")
        else:
            self._time_caption.set_text(text)
            self._time_caption.set_visible(bool(text))


class PlotPopup(QMainWindow):
    def __init__(self, title, parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(1000, 700)
        central = QWidget(self)
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)
        layout.setContentsMargins(6, 6, 6, 6)
        self.plot = TwinPlot(central)
        self.plot.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        layout.addWidget(self.plot)

def install_canvas(host):
    layout = host.layout() or QVBoxLayout(host)
    layout.setContentsMargins(0, 0, 0, 0)
    canvas = TwinPlot(host)
    canvas.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
    layout.addWidget(canvas)
    return canvas

def dialog_accepted(result):
    accepted = QDialog.DialogCode.Accepted if pyqt_version == 6 else QDialog.Accepted
    return result == accepted

class BPMTargetDialog(QDialog):
    def __init__(self, name, reference, current_target, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Desired orbit - {name}")
        self.setModal(True)
        form = QFormLayout(self)
        ref_x, ref_y = reference
        if np.isfinite(ref_x) and np.isfinite(ref_y):
            ref_text = f"Reference: x = {ref_x:+.4f} mm, y = {ref_y:+.4f} mm"
        else:
            ref_text = "Reference has not been read yet; displayed values use 0 mm."
        reference_label = QLabel(ref_text)
        reference_label.setWordWrap(True)
        form.addRow(reference_label)
        target_x, target_y = current_target
        self.x_enabled = QCheckBox("Set desired horizontal x")
        self.x_enabled.setChecked(target_x is not None)
        self.x_value = self._spin(target_x if target_x is not None else ref_x)
        form.addRow(self.x_enabled, self.x_value)

        self.y_enabled = QCheckBox("Set desired vertical y")
        self.y_enabled.setChecked(target_y is not None)
        self.y_value = self._spin(target_y if target_y is not None else ref_y)
        form.addRow(self.y_enabled, self.y_value)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Cancel | QDialogButtonBox.StandardButton.Ok)
        clear_button = QPushButton("Clear")
        buttons.addButton(clear_button, QDialogButtonBox.ButtonRole.ResetRole)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        clear_button.clicked.connect(self._set_both_to_zero)
        form.addRow(buttons)

    @staticmethod
    def _spin(value):
        spin = QDoubleSpinBox()
        spin.setRange(-1_000_000.0, 1_000_000.0)
        spin.setDecimals(6)
        spin.setSingleStep(0.05)
        spin.setSuffix(" mm")
        spin.setValue(float(value) if np.isfinite(value) else 0.0)
        return spin

    def _set_both_to_zero(self):
        self.x_value.setValue(0.0)
        self.y_value.setValue(0.0)

    def target(self):
        return (
            self.x_value.value() if self.x_enabled.isChecked() else None,
            self.y_value.value() if self.y_enabled.isChecked() else None,
        )

class MainWindow(QMainWindow, SaveOrLoad, ResponseMatrix_DFS_WFS):
    def __init__(self, interface=None, dir_name=None, start_state=None):
        super().__init__()
        ui_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "UI files", "BumpsApplication_GUI.ui")
        uic.loadUi(ui_path, self)
        self.interface = interface
        self.cwd = os.getcwd()
        self.dir_name = dir_name
        self.initial_state = start_state
        self.response_matrix_data = None
        self.reference = None
        self.targets = {}
        self._session_dir = None
        self._session_kick_count = 0
        self._drag = None
        self._desired_lines = {}
        self._desired_bpms = []
        self._desired_popup = None
        self._result_popup = None
        self._result_data = None
        self._load_logo()
        self.desired_plot = install_canvas(self.desired_plot_widget)
        self.result_plot = install_canvas(self.result_plot_widget)
        self._configure_readonly_table()
        self._set_bpm_samples()
        self._setup_corrector_controls()
        self._wire_signals()
        self._populate_devices_from_interface()
        self.loadsave_session_edit.setText(dir_name or "")
        self._draw_result_placeholder()
        self._refresh_desired_plot()
        self.setWindowTitle(f"Bumps Application — {self.interface.get_name()}")
        self.compute_button.clicked.connect(self._compute_orbit_bump)
        self.pinv_value = float(self.pinv_edit.text())
        self.beta_value = float(self.beta_edit.text())
        self.start_state = start_state if start_state is not None else interface.get_state()
        self.restore_state = None
        self._hist_desired_orbit, self._hist_predicted_orbit, self._hist_measured_orbit = [], [], []
        self._hist_desired_orbit_error, self._hist_predicted_orbit_error, self._hist_measured_orbit_error = [], [], []
        self.apply_button.clicked.connect(self._apply_orbit_bump)
        self.stop_bump_button.clicked.connect(lambda: self._stop_bump_loop("stopped by user"))
        self.restore_button.clicked.connect(self._restore_reference)
        self.restore_button.setEnabled(False)
        self.apply_button.setEnabled(False)
        self.stop_bump_button.setEnabled(False)
        self.current_delta = None
        self.corrector_names = None
        self.R = None
        self._R_matrix = None
        self._bump_running = False
        self._bump_iter = 0
        self._bump_prev_rms = None
        self._bump_stall = 0
        self._bump_tol_mm = 0.02
        self._bump_max_iter = 20
        self._bump_step_delay_ms = 1500
        self._bump_min_improvement_mm = 0.002
        self._bump_timer = QTimer(self)
        self._bump_timer.setSingleShot(True)
        self._bump_timer.timeout.connect(self._bump_step)
        self.show_response_matrix.clicked.connect(self._display_response_matrix)
        self._clock_zone_name = self._get_clock_zone()
        self._setup_machine_clock()
        self.clear_plots_button.clicked.connect(self.clear_graphs)
        self._live_orbit = None
        self._live_busy = False
        self._compute_time = None
        self._measured_time = None
        self._live_interval_ms = 2000
        self._machine_lock = threading.Lock()
        self._live_thread = QThread(self)
        self._poller = OrbitPoller(self.interface, self._machine_lock, self._live_interval_ms)
        self._poller.moveToThread(self._live_thread)
        self._live_thread.started.connect(self._poller.start)
        self._poller.bpms_ready.connect(self._on_live_bpms)
        self._live_thread.start()

    def clear_graphs(self):
        self._hist_desired_orbit.clear(), self._hist_predicted_orbit.clear(), self._hist_measured_orbit.clear()
        self._hist_desired_orbit_error.clear(), self._hist_predicted_orbit_error.clear(), self._hist_measured_orbit_error.clear()
        self._drag = None
        self._desired_lines = {}
        self._desired_bpms = []
        self._live_orbit = None
        self._compute_time = None
        self._measured_time = None
        canvases = [self.desired_plot, self.result_plot]
        if self._desired_popup is not None and self._desired_popup.isVisible():
            canvases.append(self._desired_popup.plot)
        if self._result_popup is not None and self._result_popup.isVisible():
            canvases.append(self._result_popup.plot)
        for canvas in canvases:
            canvas.clear()
            canvas.set_time_caption("")
            canvas.draw_idle()

    def _display_response_matrix(self):
        orbit_dir = self.response_dir_edit.text()
        self.handling('ComputeResponseMatrix_GUI.py', cwd=self.cwd, args=["--dir1", orbit_dir, "--compute"])

    def _clock_now(self):
        return datetime.now(ZoneInfo(self._clock_zone_name))

    def _update_machine_clock(self):
        now = self._clock_now()
        self.machine_clock_label.setText(f"{now:%Y-%m-%d %H:%M:%S} {now.tzname()}")

    def _setup_machine_clock(self):
        self._clock_timer = QTimer(self)
        self._clock_timer.setInterval(250)
        self._clock_timer.timeout.connect(self._update_machine_clock)
        self._update_machine_clock()
        self._clock_timer.start()

    def handling(self, app_name, cwd=None, args=None):
        try:
            path = os.path.join(os.path.dirname(__file__), app_name)
            workdir = os.path.expanduser(os.path.expandvars(cwd))
            proc = QProcess(self)
            proc.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
            proc.setWorkingDirectory(workdir)
            env = QProcessEnvironment.systemEnvironment()
            proc.setProcessEnvironment(env)
            argv = [path] + list(args or [])
            proc.start(sys.executable, argv)
            proc.readyReadStandardOutput.connect(
                lambda p=proc: print(bytes(p.readAllStandardOutput()).decode(errors="ignore")))
            self._procs.append(proc)
            print(workdir)

        except Exception as e:
            print(f"Error in opening Compute Response Matrix: {e}")

    def _apply_orbit_bump(self):
        if self._bump_running:
            return
        if self.R is None or self.current_delta is None or self.corrector_names is None or self._result_data is None:
            QMessageBox.information(self, "Apply", "Compute the bump first.")
            return
        if not self._set_bpm_samples():
            QMessageBox.warning(self, "BPM samples", "Enter a positive integer number of BPM samples.")
            return
        self._start_bump_loop()

    def _apply_current_delta(self):
        with self._machine_lock:
            current_corrector_settings = np.asarray(self.interface.get_correctors(self.corrector_names)["bdes"], dtype=float)
            final_current = current_corrector_settings + self.current_delta
            try:
                if self.interface.set_correctors(self.corrector_names, final_current) is False:
                    mismatch = True
                else:
                    mismatch = False
            except Exception as e:
                print(f"Error in setting correctors: {e}")
                return False
            if not mismatch:
                final_corrector_values = self.interface.get_correctors(self.corrector_names)["bdes"]
                final_state = self.interface.get_state()
        if mismatch:
            if self._bump_running:
                print("There's a mismatch in corrector values.")
            else:
                QMessageBox.warning(self, "Apply", "Not every corrector was set at the requested current. Check them on the machine.")
            return False
        self._update_corrector_table(self.corrector_names, self.current_delta, current_corrector_settings, final_corrector_values)
        orbit_measured = final_state.get_orbit(self.R.bpms)
        self._measured_time = self._clock_now()
        measured_x = np.asarray(orbit_measured["x"], dtype=float).reshape(-1)
        measured_y = np.asarray(orbit_measured["y"], dtype=float).reshape(-1)
        n_measured = max(int(orbit_measured.get("nshots", 1) or 1), 1)
        measured_x_err = np.asarray(orbit_measured["stdx"], dtype=float).reshape(-1) / np.sqrt(n_measured)
        measured_y_err = np.asarray(orbit_measured["stdy"], dtype=float).reshape(-1) / np.sqrt(n_measured)
        self._draw_result_plot(self._result_data["bpms"], self._result_data["desired_x"], self._result_data["desired_y"], self._result_data["predicted_x"], self._result_data["predicted_y"], measured_x, measured_y, measured_x_err, measured_y_err)
        self._session_kick_count += 1
        self._save_bumps_session()
        self._append_bumps_kicks(current_corrector_settings, final_corrector_values)
        final_state.timestamp = self._clock_now()
        final_state.save(filename=os.path.join(self._session_dir, "machine_status_after_bump.pkl"))
        return True

    def _bump_residual_rms(self):
        data = self._result_data
        if not data or data.get("measured_x") is None:
            return float("nan")
        dx = np.asarray(data["measured_x"], dtype=float) - np.asarray(data["desired_x"], dtype=float)
        dy = np.asarray(data["measured_y"], dtype=float) - np.asarray(data["desired_y"], dtype=float)
        diff = np.concatenate([dx, dy])
        diff = diff[np.isfinite(diff)]
        return float(np.sqrt(np.mean(diff ** 2))) if diff.size else float("nan")

    def _set_bump_controls_running(self, running):
        self.compute_button.setEnabled(not running)
        self.apply_button.setEnabled(not running and self.current_delta is not None)
        self.stop_bump_button.setEnabled(running)
        self.clear_plots_button.setEnabled(not running)
        self.read_reference_button.setEnabled(not running)
        self.restore_button.setEnabled(not running and self.restore_state is not None)

    def _start_bump_loop(self):
        self._bump_running = True
        self._bump_iter = 0
        self._bump_prev_rms = None
        self._bump_stall = 0
        self._set_bump_controls_running(True)
        self._bump_step()

    def _bump_step(self):
        if not self._bump_running:
            return
        try:
            self._bump_iter += 1
            if self.current_delta is None or not self._apply_current_delta():
                self._stop_bump_loop("Bump was not excited.")
                return
            rms = self._bump_residual_rms()
            if np.isfinite(rms) and rms <= self._bump_tol_mm:
                self._stop_bump_loop(f"Bump was set.")
                return
            if self._bump_iter >= self._bump_max_iter:
                self._stop_bump_loop(f"Bump was not set.")
                return
            if self._bump_prev_rms is not None and rms >= self._bump_prev_rms - self._bump_min_improvement_mm:
                self._bump_stall += 1
                if self._bump_stall >= 2:
                    self._stop_bump_loop(f"Bump was not set.")
                    return
            else:
                self._bump_stall = 0
            self._bump_prev_rms = rms
            if self._compute_orbit_bump(rebuild=False) is not True:
                self._stop_bump_loop("Compute failed")
                return
            self._bump_timer.start(self._bump_step_delay_ms)
        except Exception as exc:
            self._stop_bump_loop(f"error: {exc}")

    def _stop_bump_loop(self, reason="stopped"):
        was_running = self._bump_running
        self._bump_running = False
        self._bump_timer.stop()
        self._set_bump_controls_running(False)
        if was_running:
            print(f"{reason}")

    def _wire_signals(self):
        self.browse_response_button.clicked.connect(self._pick_response_directory)
        self.browse_session_button.clicked.connect(self._pick_and_load_bumps_session)
        self.loadsave_session_edit.returnPressed.connect(
            lambda: self._load_bumps_session(self.loadsave_session_edit.text()))
        self.corr_status_load.clicked.connect(self._pick_and_load_bumps_machine_status_file)
        self.corr_status_edit.returnPressed.connect(
            lambda: self._load_bumps_machine_status_file(self.corr_status_edit.text()))
        self.bpm_samples_edit.textChanged.connect(self._set_bpm_samples)
        self.response_dir_edit.returnPressed.connect(
            lambda: self._load_response_directory(self.response_dir_edit.text()))
        self.save_correctors_button.clicked.connect(
            lambda: self._saving_func(self.correctors_list, "bump_correctors.txt", "Save Correctors"))
        self.load_correctors_button.clicked.connect(self._load_correctors)
        self.clear_correctors_button.clicked.connect(self.correctors_list.clearSelection)
        self.save_bpms_button.clicked.connect(lambda: self._saving_func(self.bpms_list, "bump_bpms.txt", "Save BPMs"))
        self.load_bpms_button.clicked.connect(self._load_bpms)
        self.clear_bpms_button.clicked.connect(self.bpms_list.clearSelection)
        self.bpms_list.itemSelectionChanged.connect(self._refresh_desired_plot)
        self.bpms_list.itemDoubleClicked.connect(self._edit_bpm_target)
        self.read_reference_button.clicked.connect(self._read_reference_orbit)
        self.clear_targets_button.clicked.connect(self._clear_targets)
        self._connect_desired_plot_events(self.desired_plot)
        self.desired_plot.mpl_connect("button_press_event",
                                      lambda event: self._handle_plot_double_click(event, "desired"))
        self.result_plot.mpl_connect("button_press_event",
                                     lambda event: self._handle_plot_double_click(event, "result"))

    def _set_bpm_samples(self, _=None):
        try:
            samples = int(self.bpm_samples_edit.text())
        except ValueError:
            return False
        if samples < 1:
            return False
        self.interface.nsamples = samples
        return True

    def _setup_corrector_controls(self):
        correctors = self.interface.get_correctors()
        correctors_list = correctors['names']
        self.hcorrector_names = set(map(str, self.interface.get_hcorrectors_names() or []))
        self.vcorrector_names = set(map(str, self.interface.get_vcorrectors_names() or []))
        units_settings, sysid_kick, bpm_unit, corrs_unit = self._get_interface_units()
        self.sysid_kick = sysid_kick
        self.bpm_unit = bpm_unit
        self.corrs_unit = corrs_unit
        max_curr_h = 0.0
        max_curr_v = 0.0
        if correctors_list is not None:
            hcorrs = self.interface.get_hcorrectors_names()
            vcorrs = self.interface.get_vcorrectors_names()
            hcorr_indexes = np.array([index for index, string in enumerate(correctors_list) if string in hcorrs])
            vcorr_indexes = np.array([index for index, string in enumerate(correctors_list) if string in vcorrs])

            def clean_array(a):
                a = np.array([0 if x is None else x for x in a], dtype=float)
                a[np.isnan(a)] = 0
                return a

            max_curr_h = 1.15 * np.max(np.abs(clean_array(np.array(correctors['bdes'])[hcorr_indexes])))
            max_curr_v = 1.15 * np.max(np.abs(clean_array(np.array(correctors['bdes'])[vcorr_indexes])))
            if "bba_max_h_strength" in units_settings: max_curr_h = units_settings["bba_max_h_strength"]
            if "bba_max_v_strength" in units_settings: max_curr_v = units_settings["bba_max_v_strength"]
        self.max_horizontal_current_spinbox.setValue(max_curr_h)
        self.max_horizontal_current_spinbox.setSingleStep(0.01)
        self.max_vertical_current_spinbox.setValue(max_curr_v)
        self.max_vertical_current_spinbox.setSingleStep(0.01)

    def _get_interface_units(self):
        interface_defaults = self._get_interface_initial_settings()
        if interface_defaults is None:
            return {}, 0.01, "mm", ""
        units_settings = interface_defaults.get("units", {})
        sysid_kick = units_settings.get("sysid_corrector_kick", 0.01)
        bpm_unit = units_settings.get("bpm_position", "mm")
        corrs_unit = units_settings.get("corrector_strength", "T*mm")

        return units_settings, sysid_kick, bpm_unit, corrs_unit

    def _connect_desired_plot_events(self, canvas):
        canvas.mpl_connect("button_press_event", lambda event: self._desired_press(event, canvas))
        canvas.mpl_connect("motion_notify_event", lambda event: self._desired_motion(event, canvas))
        canvas.mpl_connect("button_release_event", lambda event: self._desired_release(event, canvas))

    def _handle_plot_double_click(self, event, plot_name):
        if getattr(event, "dblclick", False) and getattr(event, "button", None) == 1:
            self._open_plot_popup(plot_name)

    def _open_plot_popup(self, plot_name):
        popup_attr = f"_{plot_name}_popup"
        popup = getattr(self, popup_attr)
        if popup is None:
            title = "Desired bump" if plot_name == "desired" else "Resulting orbit"
            popup = PlotPopup(title, parent=self)
            setattr(self, popup_attr, popup)
            if plot_name == "desired":
                self._connect_desired_plot_events(popup.plot)
                popup.plot.mpl_connect("button_press_event",
                                       lambda event: self._handle_plot_double_click(event, "desired"))
        popup.show()
        if plot_name == "desired":
            self._refresh_desired_plot()
        else:
            self._refresh_result_plot()
        popup.raise_()
        popup.activateWindow()

    def _configure_readonly_table(self):
        mode = self.corr_delta_table.EditTrigger.NoEditTriggers if pyqt_version == 6 else self.corr_delta_table.NoEditTriggers
        self.corr_delta_table.setEditTriggers(mode)
        header = self.corr_delta_table.horizontalHeader()
        stretch = QHeaderView.ResizeMode.Stretch if pyqt_version == 6 else QHeaderView.Stretch
        header.setSectionResizeMode(stretch)

    def _load_logo(self):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "UI files", "Assets", "CERN_logo.png")
        if not os.path.isfile(path):
            return
        pixmap = QPixmap(path)
        if pixmap.isNull():
            return
        mode = Qt.TransformationMode.SmoothTransformation if pyqt_version == 6 else Qt.SmoothTransformation
        self.logo_label.setPixmap(pixmap.scaledToHeight(80, mode))

    def _populate_devices_from_interface(self):
        if self.interface is None:
            return
        state = self.initial_state or self.interface.get_state()
        self.initial_state = state
        correctors = [str(name) for name in state.get_correctors()["names"]]
        bpms = [str(name) for name in state.get_bpms()["names"]]

        self._set_list_items(self.correctors_list, correctors)
        self._set_list_items(self.bpms_list, bpms)

    @staticmethod
    def _set_list_items(widget, names):
        widget.blockSignals(True)
        widget.clear()
        widget.addItems([str(name) for name in names])
        widget.blockSignals(False)

    @staticmethod
    def _selected_names(widget):
        return [widget.item(index).text() for index in range(widget.count()) if widget.item(index).isSelected()]

    @staticmethod
    def _set_selected_names(widget, names):
        names = {str(name) for name in names}
        widget.blockSignals(True)
        widget.clearSelection()
        for index in range(widget.count()):
            item = widget.item(index)
            item.setSelected(item.text() in names)
        widget.blockSignals(False)

    def _pick_response_directory(self):
        start = self.response_dir_edit.text().strip() or self.dir_name
        directory = QFileDialog.getExistingDirectory(self, "Select response matrix data directory", start)
        if directory:
            self._load_response_directory(directory)

    @staticmethod
    def _find_response_pairs(directory):
        pairs = []
        pattern = re.compile(r"DATA_(.+)_(p|m)(\d+)\.pkl$")
        for plus_path in sorted(glob.glob(os.path.join(directory, "DATA*.pkl"))):
            match = pattern.search(os.path.basename(plus_path))
            if match is None or match.group(2) != "p":
                continue
            device, index = match.group(1), match.group(3)
            minus_path = os.path.join(directory, f"DATA_{device}_m{index}.pkl")
            if os.path.isfile(minus_path):
                pairs.append((plus_path, minus_path, device))
        return pairs

    def _load_response_directory(self, directory):
        directory = os.path.abspath(os.path.expanduser((directory or "").strip()))
        if not os.path.isdir(directory):
            QMessageBox.warning(self, "Response matrix data", "Choose an existing response matrix data directory.")
            return
        pairs = self._find_response_pairs(directory)
        if not pairs:
            QMessageBox.warning(self, "Response matrix data", "Failed to find any valid data.")
            return
        try:
            first_state = State(filename=pairs[0][0])
            interface_id = f"{type(self.interface).__module__}.{type(self.interface).__name__}"
            if first_state.get_interface_id() != interface_id:
                QMessageBox.warning(self, "Response matrix data", "This response-matrix data was created on a different interface/machine. Choose appropriate one or change interface.")
                return
            bpms = [str(name) for name in first_state.get_bpms()["names"]]
            correctors = list(dict.fromkeys(str(device) for _, _, device in pairs))
            horizontal = {str(name) for name in first_state.hcorrectors_names}
            vertical = {str(name) for name in first_state.vcorrectors_names}
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(self, "Response matrix data", f"Could not read the Response matrix data files:\n{exc}")
            return

        self.response_matrix_data = {
            "directory": directory, "pairs": pairs, "bpms": bpms, "correctors": correctors,
            "hcorrs": [name for name in correctors if name in horizontal],
            "vcorrs": [name for name in correctors if name in vertical],
        }
        self.reference = None
        self.restore_state = None
        self.restore_button.setEnabled(False)
        self.targets.clear()
        self._result_data = None
        self.response_dir_edit.setText(directory)
        self._select_devices_present_in_response_matrix()
        self._refresh_desired_plot()
        self._draw_result_placeholder()

    def _select_devices_present_in_response_matrix(self):
        response_matrix_bpms = set(self.response_matrix_data["bpms"])
        response_matrix_correctors = set(self.response_matrix_data["correctors"])
        if self.interface is None:
            self._set_list_items(self.bpms_list, self.response_matrix_data["bpms"])
            self._set_list_items(self.correctors_list, self.response_matrix_data["correctors"])

        self._set_selected_names(self.bpms_list, response_matrix_bpms)
        self._set_selected_names(self.correctors_list, response_matrix_correctors)
        selected_bpms = set(self._selected_names(self.bpms_list))
        selected_correctors = set(self._selected_names(self.correctors_list))
        missing_bpms = response_matrix_bpms - selected_bpms
        missing_correctors = response_matrix_correctors - selected_correctors
        if missing_bpms or missing_correctors:
            print(f" Devices unavailable in this interface: {len(missing_bpms)} BPMs, {len(missing_correctors)} correctors.")

    def _ordered_selected_bpms(self):
        selected = set(self._selected_names(self.bpms_list))
        if self.response_matrix_data is None:
            return []
        return [name for name in self.response_matrix_data["bpms"] if name in selected]

    def _read_reference_orbit(self):
        bpms = self._ordered_selected_bpms()
        if not bpms:
            QMessageBox.warning(self, "Reference orbit", "Load Response Matrix data and select at least one BPM first.")
            return
        if not self._set_bpm_samples():
            QMessageBox.warning(self, "BPM samples", "Enter a positive integer number of BPM samples.")
            return
        try:
            with self._machine_lock:
                state = self.interface.get_state()
            orbit = state.get_orbit(bpms)
        except Exception as exc:
            QMessageBox.critical(self, "Reference orbit", f"Could not read the BPMs:\n{exc}")
            return
        self.reference = {
            "names": [str(name) for name in orbit["names"]],
            "x": np.asarray(orbit["x"], dtype=float),
            "y": np.asarray(orbit["y"], dtype=float),
            "stdx": np.asarray(orbit["stdx"], dtype=float),
            "stdy": np.asarray(orbit["stdy"], dtype=float),
            "nshots": max(int(orbit.get("nshots", 1) or 1), 1),
        }
        self.restore_state = state
        self.restore_button.setEnabled(True)
        self._save_bumps_machine_status(state)
        state.save(filename=os.path.join(self._session_dir, "reference_machine_status.pkl"))
        self._save_bumps_session(state)
        self._refresh_desired_plot()

    def _reference_values(self, bpms):
        if self.reference is None:
            return np.full(len(bpms), np.nan), np.full(len(bpms), np.nan)
        lookup = dict(zip(self.reference["names"], zip(self.reference["x"], self.reference["y"])))
        x = np.asarray([lookup.get(name, (np.nan, np.nan))[0] for name in bpms], dtype=float)
        y = np.asarray([lookup.get(name, (np.nan, np.nan))[1] for name in bpms], dtype=float)
        return x, y

    def _target_values(self, bpms):
        x = np.asarray([self.targets.get(name, {}).get("x", np.nan) for name in bpms], dtype=float)
        y = np.asarray([self.targets.get(name, {}).get("y", np.nan) for name in bpms], dtype=float)
        return x, y

    def _reference_errors(self, bpms):
        if self.reference is None or "stdx" not in self.reference:
            nan = np.full(len(bpms), np.nan)
            return nan, nan
        n = max(int(self.reference.get("nshots", 1) or 1), 1)
        lookup = dict(zip(self.reference["names"], zip(self.reference["stdx"], self.reference["stdy"])))
        sx = np.asarray([lookup.get(name, (np.nan, np.nan))[0] for name in bpms], dtype=float)
        sy = np.asarray([lookup.get(name, (np.nan, np.nan))[1] for name in bpms], dtype=float)
        return sx / np.sqrt(n), sy / np.sqrt(n)

    def _live_values(self, bpms):
        nan = np.full(len(bpms), np.nan)
        if not self._live_orbit:
            return nan, nan, nan, nan
        n = max(int(self._live_orbit.get("nshots", 1) or 1), 1)
        lookup = dict(zip(self._live_orbit["names"],
                          zip(self._live_orbit["x"], self._live_orbit["y"],
                              self._live_orbit["stdx"], self._live_orbit["stdy"])))
        cols = [lookup.get(name, (np.nan, np.nan, np.nan, np.nan)) for name in bpms]
        x = np.asarray([c[0] for c in cols], dtype=float)
        y = np.asarray([c[1] for c in cols], dtype=float)
        sx = np.asarray([c[2] for c in cols], dtype=float) / np.sqrt(n)
        sy = np.asarray([c[3] for c in cols], dtype=float) / np.sqrt(n)
        return x, y, sx, sy

    def _edit_bpm_target(self, item):
        item.setSelected(True)
        name = item.text()
        reference_x, reference_y = self._reference_values([name])
        current = self.targets.get(name, {})
        dialog = BPMTargetDialog(name, (reference_x[0], reference_y[0]), (current.get("x"), current.get("y")), self)
        if not dialog_accepted(dialog.exec()):
            return
        target_x, target_y = dialog.target()
        if target_x is None and target_y is None:
            self.targets.pop(name, None)
        else:
            self.targets[name] = {"x": target_x, "y": target_y}
        self._refresh_desired_plot()

    def _clear_targets(self):
        self.targets.clear()
        self._refresh_desired_plot()

    def _refresh_desired_plot(self):
        if self._drag is not None:
            return
        bpms = self._ordered_selected_bpms()
        self._desired_bpms = bpms
        canvases = [self.desired_plot]
        if self._desired_popup is not None and self._desired_popup.isVisible():
            canvases.append(self._desired_popup.plot)
        self._desired_lines = {}
        for canvas in canvases:
            self._draw_desired_plot(canvas, bpms)

    def _draw_desired_plot(self, canvas, bpms):
        canvas.clear()
        if not bpms:
            canvas.ax_x.set_title("Load Response matrix data to see the selected BPMs.")
            canvas.draw_idle()
            return
        reference_x, reference_y = self._reference_values(bpms)
        target_x, target_y = self._target_values(bpms)
        live_x, live_y, live_sx, live_sy = self._live_values(bpms)
        ref_sx, ref_sy = self._reference_errors(bpms)
        index = np.arange(len(bpms))
        lines = {}
        for axis, reference, target, live, live_sem, ref_sem, plane in (
                (canvas.ax_x, reference_x, target_x, live_x, live_sx, ref_sx, "x"),
                (canvas.ax_y, reference_y, target_y, live_y, live_sy, ref_sy, "y")):
            base = np.where(np.isfinite(reference), reference, 0.0)
            desired_rel = np.where(np.isfinite(target), target - base, 0.0)
            measured_rel = live - base
            diff_err = np.sqrt(np.nan_to_num(live_sem) ** 2 + np.nan_to_num(ref_sem) ** 2)
            axis.axhline(0.0, color="0.55", linestyle="--", linewidth=1, label="reference", zorder=1)
            if self.reference is not None and np.any(np.isfinite(measured_rel)):
                axis.errorbar(index, measured_rel, yerr=diff_err, color="tab:green", marker="o",
                              markersize=4, capsize=3, elinewidth=1, label="measured − reference orbit", zorder=3)
            line, = axis.plot(index, desired_rel, color="tab:blue", marker="o", markersize=7,
                              label="desired bump", zorder=4)
            lines[plane] = line
            axis.set_ylabel(f"Δ{plane} [mm]")
            axis.legend(fontsize=7, loc="best")
            span = np.concatenate([desired_rel, measured_rel[np.isfinite(measured_rel)]])
            scale = max(0.5, float(np.max(np.abs(span))) if span.size else 0.5)
            axis.set_ylim(-1.25 * scale, 1.25 * scale)
        canvas.ax_y.set_xlabel("BPM")
        self._set_bpm_ticks(canvas, bpms)
        self._desired_lines[canvas] = lines
        stamp = self._live_orbit.get("time") if self._live_orbit else None
        canvas.set_time_caption(f"{stamp:%Y-%m-%d %H:%M:%S %Z}" if stamp is not None else "waiting for first BPM read…")
        canvas.draw_idle()

    @staticmethod
    def _set_bpm_ticks(canvas, bpms):
        step = max(1, len(bpms) // 12)
        ticks = np.arange(0, len(bpms), step)
        labels = [bpms[index] for index in ticks]
        for axis in (canvas.ax_x, canvas.ax_y):
            axis.set_xticks(ticks)
            axis.set_xticklabels(labels, rotation=55, ha="right", fontsize=6)
            axis.tick_params(axis="x", which="both", labelbottom=True)

    def _nearest_marker(self, event, canvas):
        if event.inaxes is canvas.ax_x:
            plane = "x"
        elif event.inaxes is canvas.ax_y:
            plane = "y"
        else:
            return None
        line = self._desired_lines.get(canvas, {}).get(plane)
        if line is None or not self._desired_bpms:
            return None
        coordinates = event.inaxes.transData.transform(np.column_stack([line.get_xdata(), line.get_ydata()]))
        distance = np.hypot(coordinates[:, 0] - event.x, coordinates[:, 1] - event.y)
        row = int(np.argmin(distance))
        return (plane, row) if distance[row] <= 12.0 else None

    def _desired_press(self, event, canvas):
        if event.inaxes is None or event.x is None or event.y is None:
            return
        if getattr(event, "dblclick", False):
            return
        hit = self._nearest_marker(event, canvas)
        if hit is None:
            return
        plane, row = hit
        name = self._desired_bpms[row]
        if event.button == 3:
            self.targets.setdefault(name, {})[plane] = None
            if all(value is None for value in self.targets[name].values()):
                self.targets.pop(name, None)
            self._refresh_desired_plot()
        elif event.button == 1:
            self._drag = (canvas, plane, row)

    def _desired_motion(self, event, canvas):
        if self._drag is None or event.x is None or event.y is None:
            return
        drag_canvas, plane, row = self._drag
        if drag_canvas is not canvas:
            return
        line = self._desired_lines.get(canvas, {}).get(plane)
        if line is None:
            return
        axis = line.axes
        value = float(axis.transData.inverted().transform((event.x, event.y))[1])
        y_values = np.asarray(line.get_ydata(), dtype=float)
        y_values[row] = value
        line.set_ydata(y_values)
        lower, upper = axis.get_ylim()
        if value < lower + 0.1 * (upper - lower) or value > upper - 0.1 * (upper - lower):
            scale = 1.25 * max(abs(value), abs(lower), abs(upper), 0.5)
            axis.set_ylim(-scale, scale)
        axis.set_title(f"{self._desired_bpms[row]}: Δ = {value:+.3f} mm", fontsize=8)
        canvas.draw_idle()

    def _desired_release(self, _event, canvas):
        if self._drag is None:
            return
        drag_canvas, plane, row = self._drag
        if drag_canvas is not canvas:
            return
        self._drag = None
        line = self._desired_lines.get(canvas, {}).get(plane)
        if line is None:
            return
        value = float(np.asarray(line.get_ydata(), dtype=float)[row])
        name = self._desired_bpms[row]
        reference = self._reference_values([name])[0 if plane == "x" else 1][0]
        base = float(reference) if np.isfinite(reference) else 0.0
        self.targets.setdefault(name, {})[plane] = base + value
        line.axes.set_title("")
        self._refresh_desired_plot()

    def _draw_result_placeholder(self):
        canvases = [self.result_plot]
        if self._result_popup is not None and self._result_popup.isVisible():
            canvases.append(self._result_popup.plot)
        for canvas in canvases:
            canvas.clear()
            for axis, plane in ((canvas.ax_x, "x"), (canvas.ax_y, "y")):
                axis.set_ylabel(f"{plane} [mm]")
            canvas.ax_y.set_xlabel("BPM")
            canvas.draw_idle()

    def _refresh_result_plot(self):
        if self._result_data is None:
            self._draw_result_placeholder()
            return
        self._draw_result_plot(**self._result_data)

    def _on_live_bpms(self, raw):
        if self._live_busy or self._drag is not None:
            return
        self._live_busy = True
        try:
            bpms = self._ordered_selected_bpms()
            if self.R is not None:
                bpms = list(dict.fromkeys(list(self.R.bpms) + bpms))
            if not bpms:
                return
            orbit = State(bpms=raw).get_orbit(bpms)
            self._measured_time = self._clock_now()
            self._live_orbit = {
                "names": [str(name) for name in orbit["names"]],
                "x": np.asarray(orbit["x"], dtype=float),
                "y": np.asarray(orbit["y"], dtype=float),
                "stdx": np.asarray(orbit["stdx"], dtype=float),
                "stdy": np.asarray(orbit["stdy"], dtype=float),
                "nshots": max(int(orbit.get("nshots", 1) or 1), 1),
                "time": self._measured_time,
            }
            self._refresh_desired_plot()
            if self._result_data is not None:
                mx, my, ex, ey = self._live_values(self._result_data["bpms"])
                data = dict(self._result_data)
                data.update(measured_x=mx, measured_y=my, measured_x_err=ex, measured_y_err=ey)
                self._draw_result_plot(**data)
        except Exception as exc:  # noqa: BLE001
            print(f"Live orbit update skipped: {exc}")
        finally:
            self._live_busy = False

    def _expand_path(self, path):
        expanded_path = (path or "").strip()
        expanded_path = os.path.expandvars(os.path.expanduser(expanded_path))
        return os.path.abspath(os.path.normpath(expanded_path))

    def _get_response_matrix(self, directory):
        directory = self._expand_path(directory)
        datafiles = sorted(glob.glob(os.path.join(directory, "DATA*.pkl")))
        if not datafiles:
            QMessageBox.warning(self, "Error", "No data files found")
            return
        S = State(filename=datafiles[0])
        interface_id = f"{type(self.interface).__module__}.{type(self.interface).__name__}"
        if S.get_interface_id() != interface_id:
            QMessageBox.warning(self, "Response matrix data", "This response-matrix data was created on a different interface/machine. Choose appropriate one or change interface.")
            return
        self.sequence = S.get_sequence()
        correctors = [self.correctors_list.item(i).text() for i in range(self.correctors_list.count()) if
                      self.correctors_list.item(i).isSelected()]
        bpms = [self.bpms_list.item(i).text() for i in range(self.bpms_list.count()) if
                self.bpms_list.item(i).isSelected()]

        if not correctors:
            for i in range(self.correctors_list.count()):
                self.correctors_list.item(i).setSelected(True)
            correctors = self.correctors

        if not bpms:
            for i in range(self.bpms_list.count()):
                self.bpms_list.item(i).setSelected(True)
            bpms = self.bpms

        Rxx, Ryy, Rxy, Ryx, Bx, By, hcorrs, vcorrs, bpms = self._compute_response_matrix_from_directory(
            directory=directory, correctors=correctors, bpms=bpms,
            triangular=bool(self.triangular_checkbox.isChecked()))
        R = Response()
        R.bpms = bpms
        R.hcorrs = hcorrs
        R.vcorrs = vcorrs
        R.Rxx = Rxx
        R.Rxy = Rxy
        R.Ryx = Ryx
        R.Ryy = Ryy
        print(R.Ryy)
        R.Bx = Bx
        R.By = By
        return R

    def _compute_orbit_bump(self, *_args, rebuild=True):
        if not self._set_bpm_samples():
            QMessageBox.warning(self, "BPM samples", "Enter valid number of BPM samples.")
            return False
        try:
            self.pinv_value = float(self.pinv_edit.text())
        except ValueError:
            QMessageBox.warning(self, "PINV tolerance", "Enter a valid PINV tolerance.")
            return False
        try:
            self.beta_value = float(self.beta_edit.text())
        except ValueError:
            QMessageBox.warning(self, "Beta parameter", "Enter a valid beta parameter.")
            return False

        if rebuild or self.R is None or self._R_matrix is None:
            self._save_bumps_machine_status()
            self.R = self._get_response_matrix(self.response_dir_edit.text())
            if self.R is None:
                QMessageBox.warning(self, "Error", "No data files found")
                return False
            self._R_matrix = np.block([
                [self.R.Rxx, self.R.Rxy],
                [self.R.Ryx, self.R.Ryy],
            ])
            self.corrector_names = self.R.hcorrs + self.R.vcorrs
        R_matrix = self._R_matrix
        max_curr_h = self.max_horizontal_current_spinbox.value()
        max_curr_v = self.max_vertical_current_spinbox.value()

        def clamp(val, max_val):
            val = np.asarray(val, dtype=float)
            max_val = np.asarray(max_val, dtype=float)
            result = val.copy()
            finite = np.isfinite(max_val) & (max_val > 0.0)
            result[finite] = np.clip(result[finite], -max_val[finite], max_val[finite])
            return result

        with self._machine_lock:
            raw_bpms = self.interface.get_bpms()
        current_orbit = State(bpms=raw_bpms).get_orbit(self.R.bpms)
        current_x = np.asarray(current_orbit["x"], dtype=float).reshape(-1)
        current_y = np.asarray(current_orbit["y"], dtype=float).reshape(-1)

        target_x, target_y = self._target_values(self.R.bpms)
        reference_x, reference_y = self._reference_values(self.R.bpms)
        baseline_x = np.where(np.isfinite(reference_x), reference_x, 0.0)
        baseline_y = np.where(np.isfinite(reference_y), reference_y, 0.0)

        desired_x = np.where(np.isfinite(target_x), target_x, baseline_x)
        desired_y = np.where(np.isfinite(target_y), target_y, baseline_y)

        delta_orbit = np.concatenate([
            desired_x - current_x,
            desired_y - current_y,
        ])

        if self.beta_value > 0:
            delta = np.linalg.solve(
                R_matrix.T @ R_matrix + self.beta_value * np.eye(R_matrix.shape[1]),
                R_matrix.T @ delta_orbit,
            )
        else:
            delta = np.linalg.pinv(R_matrix, rcond=self.pinv_value) @ delta_orbit

        nh = len(self.R.hcorrs)
        delta_x = np.asarray(delta[:nh], dtype=float)
        delta_y = np.asarray(delta[nh:], dtype=float)
        max_vals_x = np.full(delta_x.shape, max_curr_h, dtype=float)
        max_vals_y = np.full(delta_y.shape, max_curr_v, dtype=float)
        max_vals = np.concatenate([max_vals_x, max_vals_y])

        with self._machine_lock:
            current_corrector_settings = np.asarray(self.interface.get_correctors(self.corrector_names)["bdes"], dtype=float)
        new_bdes = current_corrector_settings + delta
        new_bdes = clamp(new_bdes, max_vals)
        delta = new_bdes - current_corrector_settings
        self.current_delta = delta
        self._update_corrector_table(self.corrector_names, delta, current_corrector_settings, current_corrector_settings + delta)
        # predicted
        orbit_predicted = R_matrix @ delta
        predicted_x = current_x + orbit_predicted[:len(self.R.bpms)]
        predicted_y = current_y + orbit_predicted[len(self.R.bpms):]

        prev = self._result_data if self._result_data else {}
        mx, my = prev.get("measured_x"), prev.get("measured_y")
        mex, mey = prev.get("measured_x_err"), prev.get("measured_y_err")
        if mx is not None and len(mx) != len(self.R.bpms):
            mx = my = mex = mey = None
        self._compute_time = self._clock_now()
        self._draw_result_plot(self.R.bpms, desired_x, desired_y, predicted_x, predicted_y, mx, my, mex, mey)
        self._save_bumps_session()
        self.apply_button.setEnabled(not self._bump_running)
        return True

    def _update_corrector_table(self, corrector_names, delta, reference, bdes):
        # corrector | kick | reference | new setpoint
        self.corr_delta_table.setRowCount(len(corrector_names))
        for row, (corrector, kick, reference, bdes) in enumerate(zip(corrector_names, delta, reference, bdes)):
            self.corr_delta_table.setItem(row, 0, QTableWidgetItem(str(corrector)))
            self.corr_delta_table.setItem(row, 1, QTableWidgetItem(f"{kick:.6g}"))
            self.corr_delta_table.setItem(row, 2, QTableWidgetItem(f"{float(reference):.6g}"))
            self.corr_delta_table.setItem(row, 3, QTableWidgetItem(f"{float(bdes):.6g}"))

    def _restore_reference(self):
        if self.restore_state is None:
            QMessageBox.information(self, "Restore reference", "Read the reference orbit first.")
            return
        print("Restoring reference corrector settings...")
        with self._machine_lock:
            restored = self.interface.restore_correctors_state(self.restore_state)
        if restored is False:
            print("Warning: not every corrector was confirmed back at its reference current.")
            QMessageBox.warning(self, "Restore initial settings",
                                "Some correctors were not confirmed back at their reference current within the readback tolerance. Check them on the machine before the next correction.")
        self.reset_ref_orb = True
        print("Reference corrector settings restored.")

    def _draw_result_plot(self, bpms, desired_x, desired_y, predicted_x, predicted_y, measured_x=None, measured_y=None,
                          measured_x_err=None, measured_y_err=None):
        self._result_data = {
            "bpms": list(bpms),
            "desired_x": np.asarray(desired_x, dtype=float).copy(),
            "desired_y": np.asarray(desired_y, dtype=float).copy(),
            "predicted_x": np.asarray(predicted_x, dtype=float).copy(),
            "predicted_y": np.asarray(predicted_y, dtype=float).copy(),
            "measured_x": None if measured_x is None else np.asarray(measured_x, dtype=float).copy(),
            "measured_y": None if measured_y is None else np.asarray(measured_y, dtype=float).copy(),
            "measured_x_err": None if measured_x_err is None else np.asarray(measured_x_err, dtype=float).copy(),
            "measured_y_err": None if measured_y_err is None else np.asarray(measured_y_err, dtype=float).copy(),
        }
        canvases = [self.result_plot]
        if self._result_popup is not None and self._result_popup.isVisible():
            canvases.append(self._result_popup.plot)

        index = np.arange(len(bpms))
        for canvas in canvases:
            canvas.clear()
            for axis, desired, predicted, measured, merr, plane in (
                    (canvas.ax_x, desired_x, predicted_x, measured_x, measured_x_err, "x"),
                    (canvas.ax_y, desired_y, predicted_y, measured_y, measured_y_err, "y"),
            ):
                axis.plot(index, desired, "--o", label="desired")
                axis.plot(index, predicted, ":o", label="predicted")
                if measured is not None:
                    if merr is not None:
                        axis.errorbar(index, measured, yerr=merr, fmt="-o", capsize=3, elinewidth=1, label="measured")
                    else:
                        axis.plot(index, measured, "-o", label="measured")
                axis.set_ylabel(f"{plane} [mm]")
                axis.legend(fontsize=7, loc="best")
            canvas.ax_y.set_xlabel("BPM")
            self._set_bpm_ticks(canvas, bpms)
            parts = []
            if (measured_x is not None or measured_y is not None) and self._measured_time is not None:
                parts.append(f"{self._measured_time:%Y-%m-%d %H:%M:%S %Z}")
            canvas.set_time_caption("      ".join(parts))
            canvas.draw_idle()

    def _clear_graphs(self):
        self._hist_desired_orbit.clear(), self._hist_predicted_orbit.clear(), self._hist_measured_orbit.clear()
        self._hist_desired_orbit_error.clear(), self._hist_predicted_orbit_error.clear(), self._hist_measured_orbit_error.clear()
        self._plot_series(self._hist_desired_orbit_ax, self._hist_desired_orbit_canvas, values_x=[], values_y=[], title=None)
        self._plot_series(self._hist_predicted_orbit_ax, self._hist_predicted_orbit_canvas, values_x=[], values_y=[], title=None)
        self._plot_series(self._hist_measured_orbit_ax, self._hist_measured_orbit_canvas, values_x=[], values_y=[], title=None)
        self._refresh_all_plot_popups()

    def closeEvent(self, event):
        self._stop_bump_loop("window closed")
        if self._live_thread.isRunning():
            blocking = Qt.ConnectionType.BlockingQueuedConnection if pyqt_version == 6 else Qt.BlockingQueuedConnection
            QMetaObject.invokeMethod(self._poller, "stop", blocking)
            self._live_thread.quit()
            self._live_thread.wait(5000)
        super().closeEvent(event)

    def _get_interface_initial_settings(self):
        interface_class_name = self.interface.__class__.__name__
        interface_module_name = self.interface.__class__.__module__

        for machine_interfaces in INTERFACE_SETUP.values():
            for interface_defaults in machine_interfaces:
                if (interface_defaults.get("class_name") == interface_class_name) and (
                        interface_defaults.get("module") == interface_module_name):
                    return interface_defaults
        return None

    def _get_clock_zone(self):
        interface_defaults = self._get_interface_initial_settings() or {}
        timezone = interface_defaults.get("clock_timezone", "Europe/Zurich")
        return timezone

if __name__ == "__main__":
    app = QApplication(sys.argv)
    from Backend import SelectInterface

    dialog = SelectInterface.choose_acc_and_interface()
    if dialog is None:
        print("Selection cancelled.")
        sys.exit(1)

    I = dialog
    project_name = I.get_name()
    nominal_state = None
    print(f"Selected interface: {project_name}")
    time_str = datetime.now().strftime("%Y%m%d_%H%M%S")
    dir_name = os.path.expanduser("~/CERN-Flight_Simulator-Data")
    w = MainWindow(interface=I, dir_name=dir_name)

    w.show()
    sys.exit(app.exec())