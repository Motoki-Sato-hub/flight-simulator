"""SysID application dedicated to quadrupole-mover measurements."""

import os
import sys
import time
from datetime import datetime

import numpy as np

try:
    from PyQt6.QtWidgets import QApplication, QMessageBox
except ImportError:
    from PyQt5.QtWidgets import QApplication, QMessageBox

from SysID_GUI import (
    MainWindow as CorrectorSysIDWindow,
    Worker as CorrectorWorker,
    finite_abs_max,
    update_amplitude,
)
from Backend.ActuatorMode import ActuatorMode


class QMWorker(CorrectorWorker):
    """Measure response by moving each quadrupole in X and Y."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.actuator_mode = ActuatorMode.QM

    def run(self):
        self.running = True
        self.paused = False
        self.progress_value = 0
        interface = self.interface
        hkicks = self.hkicks
        vkicks = self.vkicks

        pending_steps = sum(
            1
            for iteration in range(self.Niter)
            for magnet in self.correctors
            for axis in ("x", "y")
            if iteration > 0
            or not (
                os.path.isfile(self._data_filename(magnet, axis, "+", iteration))
                and os.path.isfile(self._data_filename(magnet, axis, "-", iteration))
            )
        )
        total_steps = max(pending_steps, 1)

        for iteration in range(self.Niter):
            if not self.running:
                break
            if self.paused:
                self._await_user()

            for index, magnet in enumerate(self.correctors):
                if not self.running:
                    break
                if self.paused:
                    self._await_user()

                state0 = self._quadrupole_state(interface, magnet)
                if state0 is None:
                    continue
                x0, y0, roll0 = state0

                for axis in ("x", "y"):
                    if not self.running:
                        break
                    if self.paused:
                        self._await_user()

                    amplitude = float(hkicks[index] if axis == "x" else vkicks[index])
                    target = float(self.max_osc_h if axis == "x" else self.max_osc_v)
                    max_range = float(self.max_curr_h if axis == "x" else self.max_curr_v)
                    if max_range > 0:
                        amplitude = min(amplitude, max_range)

                    plus_file = self._data_filename(magnet, axis, "+", iteration)
                    minus_file = self._data_filename(magnet, axis, "-", iteration)
                    measured = False
                    try:
                        plus, plus_measured = self._measure_state(
                            interface, magnet, axis, x0, y0, roll0, amplitude, plus_file, iteration, "+"
                        )
                        minus, minus_measured = self._measure_state(
                            interface, magnet, axis, x0, y0, roll0, amplitude, minus_file, iteration, "-"
                        )
                        measured = plus_measured or minus_measured
                    finally:
                        try:
                            interface.apply_qmag_xyroll(magnet, x0, y0, roll0)
                        except Exception as exc:
                            print(f"WARNING: failed to restore {magnet} mover state ({exc})")

                    orbit_plus = plus.get_orbit(self.bpms)
                    orbit_minus = minus.get_orbit(self.bpms)
                    diff_x = (orbit_plus["x"] - orbit_minus["x"]) / 2.0
                    diff_y = (orbit_plus["y"] - orbit_minus["y"]) / 2.0
                    nsamples = max(1, np.asarray(orbit_plus["stdx"]).size)
                    err_x = np.sqrt(np.square(orbit_plus["stdx"]) + np.square(orbit_minus["stdx"])) / np.sqrt(nsamples)
                    err_y = np.sqrt(np.square(orbit_plus["stdy"]) + np.square(orbit_minus["stdy"])) / np.sqrt(nsamples)

                    if measured:
                        self.plot_data.emit(orbit_plus, diff_x, err_x, diff_y, err_y, self.bpms, f"{magnet}:{axis}")
                        self.progress_value += 1
                        self.progress.emit(int(self.progress_value / total_steps * 100))

                    observed = max(finite_abs_max(diff_x), finite_abs_max(diff_y))
                    amplitude = 0.8 * update_amplitude(amplitude, observed, target, max_range) + 0.2 * amplitude
                    if axis == "x":
                        hkicks[index] = amplitude
                    else:
                        vkicks[index] = amplitude
                    self._save_kicks(hkicks, vkicks)

                    if measured:
                        time.sleep(0.2)

        self.running = False
        self.finished.emit()

    def _data_filename(self, magnet, axis, sign, iteration):
        suffix = "p" if sign == "+" else "m"
        return os.path.join(self.output_dir, f"DATA_{magnet}_{axis}_{suffix}{iteration:04d}.pkl")

    @staticmethod
    def _quadrupole_state(interface, magnet):
        try:
            quadrupole = interface.get_quadrupoles(magnet)
        except TypeError:
            quadrupole = interface.get_quadrupoles([magnet])
        except Exception as exc:
            print(f"Skipping {magnet}: failed to read quadrupole state ({exc})")
            return None
        if not len(quadrupole.get("names", [])) or not all(
            key in quadrupole for key in ("xdes", "ydes", "rolldes")
        ):
            print(f"Skipping {magnet}: mover readback is incomplete.")
            return None
        return tuple(float(np.asarray(quadrupole[key])[0]) for key in ("xdes", "ydes", "rolldes"))

    def _measure_state(self, interface, magnet, axis, x0, y0, roll0, amplitude, filename, iteration, sign):
        measure = iteration > 0 or not os.path.isfile(filename)
        if not measure:
            return self.state_class(filename=filename), False
        displacement = amplitude if sign == "+" else -amplitude
        if axis == "x":
            interface.apply_qmag_xyroll(magnet, x0 + displacement, y0, roll0)
        else:
            interface.apply_qmag_xyroll(magnet, x0, y0 + displacement, roll0)
        state = interface.get_state(include_screens=False)
        state.save(filename=filename)
        return state, True

    def _save_kicks(self, hkicks, vkicks):
        with open(os.path.join(self.output_dir, "kicks.txt"), "w") as file:
            for index, magnet in enumerate(self.correctors):
                file.write(f"{magnet} {hkicks[index]} {vkicks[index]}\n")


class MainWindow(CorrectorSysIDWindow):
    """SysID window with quadrupole movers as its only actuator type."""

    ACTUATOR_MODE = ActuatorMode.QM
    WORKER_CLASS = QMWorker

    def __init__(self, interface, dir_name):
        self.qm_corrs = [str(name) for name in interface.get_quadrupole_movers_names()]
        super().__init__(interface, dir_name)
        self.setWindowTitle(f"SYSID QM - {self.interface.__class__.__name__} [Idle]")
        self.save_correctors_button.setText("Save quadrupoles")
        self.load_correctors_button.setText("Load quadrupoles")
        self.clear_correctors_button.setText("Clear quadrupoles")
        self.pattern_corrs_input.setPlaceholderText("e.g. QF*, QD*")

    def _window_title_prefix(self):
        return "SYSID QM"

    def _available_actuators(self):
        return self._sort_elements(self.qm_corrs, which="sequence")

    def _apply_corrector_checkbox_selection(self):
        pass

    def _refresh_actuator_labels(self):
        self.correctorsGroup.setTitle("Quadrupole movers")
        self.initial_hkick_label.setText("Initial X [um]")
        self.initial_vkick_label.setText("Initial Y [um]")
        self.current_label.setText("Max mover range [um]")
        self.horizontal_current_label.setText("X:")
        self.vertical_current_label.setText("Y:")
        self.max_horizontal_current_spinbox.setMaximum(1e6)
        self.max_vertical_current_spinbox.setMaximum(1e6)
        self.max_horizontal_current_spinbox.setSingleStep(10.0)
        self.max_vertical_current_spinbox.setSingleStep(10.0)
        self.max_horizontal_current_spinbox.setValue(1000.0)
        self.max_vertical_current_spinbox.setValue(1000.0)
        self.initial_hkick_settings.setText("100")
        self.initial_vkick_settings.setText("100")
        self.excursion_label.setText(f"Target orbit excursion ({self.bpm_unit})")
        self.choose_mode.setCurrentText("Orbit Correction")
        self.choose_mode.setEnabled(False)
        self.select_h_corrs_checkbox.setEnabled(False)
        self.select_v_corrs_checkbox.setEnabled(False)

    def _validate_start(self):
        if not hasattr(self.interface, "apply_qmag_xyroll") or not hasattr(self.interface, "get_quadrupoles"):
            QMessageBox.critical(
                self,
                "QM mode not available",
                "This interface does not expose quadrupole-mover controls.",
            )
            return False
        if not self.qm_corrs:
            QMessageBox.critical(self, "QM mode not available", "No quadrupole movers are available.")
            return False
        return True

    def _sort_actuators(self, names):
        return self._sort_elements(names, which="sequence")

    def _restore_actuators_state(self, machine_state):
        # QMWorker restores each mover immediately after its +/- excitation.
        return None

    def _actuator_selection_filename(self):
        return "quadrupole_movers.txt"

    def _actuator_label(self):
        return "Quadrupole mover"


def main():
    app = QApplication(sys.argv)
    from Backend import SelectInterface

    interface = SelectInterface.choose_acc_and_interface()
    if interface is None:
        return 1

    time_str = datetime.now().strftime("%Y%m%d_%H%M%S")
    dir_name = os.path.expanduser(
        f"~/CERN-Flight_Simulator-Data/{interface.get_name()}_{time_str}_QM"
    )
    window = MainWindow(interface, dir_name)
    window.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
