import os, sys, matplotlib, pickle, time
from datetime import datetime
from Backend.ResponseMatrix_DFS_WFS import ResponseMatrix_DFS_WFS
from Backend.State import State
import numpy as np
matplotlib.use("QtAgg")
try:
    pyqt_version = 6
    from PyQt6.QtWidgets import (
        QDialog, QVBoxLayout, QDialogButtonBox,
        QRadioButton, QLabel,QFileDialog, QApplication
        )
    from PyQt6.QtCore import QEvent, Qt
except ImportError:
    pyqt_version = 5
    from PyQt5.QtWidgets import (
        QDialog, QVBoxLayout, QDialogButtonBox,
        QRadioButton, QLabel, QFileDialog, QApplication
    )
    from PyQt5.QtCore import QEvent, Qt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from scipy.optimize import least_squares
from Backend.SaveOrLoad import SaveOrLoad

class QuadrupoleScan(SaveOrLoad):
    def _new_scan_session_dir(self, quad_names, is_quad_scan):
        time_str = datetime.now().strftime("%Y%m%d_%H%M%S")
        root_dir = os.path.expanduser(os.path.expandvars("~/CERN-Flight_Simulator-Data/"))
        interface_name = self.interface.get_name()
        if is_quad_scan and quad_names:
            name = f"Quadrupole_scan_{interface_name}_{'-'.join(quad_names)}_{time_str}"
        else:
            name = f"EM_{interface_name}{time_str}_session_settings"
        session_dir = os.path.join(root_dir, name)
        os.makedirs(session_dir, exist_ok=True)
        return session_dir

    def _quad_value_unit(self, quad_name):
        quadrupoles = self.interface.get_quadrupoles([quad_name])
        return str(quadrupoles.get("value_unit", "1/m"))

    def current_to_quad_setpoint(self, quad_name, current_A, quad_value_unit):
        if quad_value_unit == "A":
            return float(current_A)
        current_to_k1l = getattr(self.interface, "current_to_k1l", None)
        if not callable(current_to_k1l):
            raise ValueError(
                f"{self.interface.get_name()} sets quadrupoles in {quad_value_unit}, but it has no "
                f"current_to_k1l() calibration, so a scan range given in amperes cannot be applied.")
        return float(current_to_k1l(quad_name, float(current_A)))

    def quad_setpoint_to_current(self, quad_name, value, quad_value_unit):
        if quad_value_unit == "A":
            return float(value)
        k1l_to_current = getattr(self.interface, "k1l_to_current", None)
        if not callable(k1l_to_current):
            return np.nan
        try:
            return float(k1l_to_current(quad_name, float(value)))
        except Exception:
            return np.nan

    @staticmethod
    def _unique_scan_values(values, relative_tolerance=1e-4):
        values = np.asarray(sorted(float(value) for value in np.ravel(values)), dtype=float)
        if values.size == 0:
            return values
        tolerance = max(abs(float(values[-1] - values[0])) * relative_tolerance, 1e-12)
        grouped = [float(values[0])]
        for value in values[1:]:
            if float(value) - grouped[-1] > tolerance:
                grouped.append(float(value))
        return np.asarray(grouped, dtype=float)

    def _scan_grid_per_screen(self, quad_name, screens, current_min, current_max, steps, quad_value_unit, screen_current_ranges=None):
        screen_current_ranges = dict(screen_current_ranges or {})
        currents_per_screen, values_per_screen = [], []
        for screen_name in screens:
            low, high = screen_current_ranges.get(screen_name, (current_min, current_max))
            low, high = float(low), float(high)
            if high <= low:
                raise ValueError(f"Maximum current must be larger than minimum current for screen {screen_name}")
            currents = np.linspace(low, high, int(steps))
            currents_per_screen.append(currents)
            values_per_screen.append(np.array([self.current_to_quad_setpoint(quad_name, current, quad_value_unit) for current in currents], dtype=float))

        quad_values = self._unique_scan_values(np.concatenate(values_per_screen))
        rows_per_screen = [np.array([int(np.argmin(np.abs(quad_values - value))) for value in values]) for values in values_per_screen]
        return quad_values, values_per_screen, currents_per_screen, rows_per_screen

    def run_scan(self, quad_name, screens, current_min, current_max, steps, nshots, screen_current_ranges=None, reference_screen=None, progress_callback=None, resume_states_dir=None):
        steps = int(steps)
        if isinstance(quad_name, str):
            quad_names = [quad_name]
        else:
            quad_names = list(quad_name or [])

        if len(quad_names) == 0:
            raise ValueError("At least one quadrupole must be provided")

        if resume_states_dir:
            self.dir_name = os.path.dirname(os.path.normpath(resume_states_dir))
        else:
            self.dir_name = self._new_scan_session_dir(quad_names=quad_names, is_quad_scan=(steps > 0))
        self.session_directory.setText(self.dir_name)
        preserve_status = getattr(self, "_preserve_quadrupole_status_files", None)
        if callable(preserve_status):
            preserve_status(self.dir_name)

        if steps == 0:
            if len(quad_names) != 1:
                raise ValueError("For Linear Response, choose exactly one quadrupole.")
            return self._run_single_scan(quad_name=quad_names[0], screens=screens,
                current_min=0.0, current_max=0.0, steps=0, nshots=nshots,
                reference_screen=reference_screen, progress_callback=progress_callback)

        if len(quad_names) == 1:
            print(f"{quad_names[0]} is going to be scanned")
            return self._run_single_scan(quad_name=quad_names[0], screens=screens, current_min=current_min, current_max=current_max,
                steps=steps, nshots=nshots, screen_current_ranges=screen_current_ranges,
                reference_screen=reference_screen, progress_callback=progress_callback, resume_states_dir=resume_states_dir)

        per_quad_sessions = []
        cancelled = False
        skipped_quadrupoles = []
        completed_quadrupoles = []

        total_quads = len(quad_names)
        for quad_idx, quad_name in enumerate(quad_names):
            if getattr(self, "_scan_stop_requested", False) or getattr(self, "_cancel", False):
                cancelled = True
                break

            def _wrapped_progress(session_partial, i, nsteps, _quad_idx=quad_idx, _quad_name=quad_name):
                merged_partial = {
                    "mode": "multi_quad_scan",
                    "quadrupoles": quad_names,
                    "current_quadrupole": _quad_name,
                    "current_quadrupole_index": int(_quad_idx),
                    "total_quadrupoles": int(total_quads),
                    "completed_quadrupoles": list(completed_quadrupoles),
                    "skipped_quadrupoles": list(skipped_quadrupoles),
                    "per_quad_sessions": per_quad_sessions + [session_partial],
                    "cancelled": bool(session_partial.get("cancelled", False)),
                }
                if progress_callback is not None:
                    progress_callback(merged_partial, i, nsteps)

            try:
                single_session = self._run_single_scan(quad_name=quad_name, screens=screens, current_min=current_min, current_max=current_max,
                    steps=steps, nshots=nshots, screen_current_ranges=screen_current_ranges,
                    reference_screen=reference_screen, progress_callback=_wrapped_progress)

            except ValueError as e:
                msg = str(e)
                if "zero quadrupole value" in msg:
                    skipped_quadrupoles.append({"quad_name": quad_name, "reason": msg})
                    continue
                raise

            per_quad_sessions.append(single_session)
            completed_quadrupoles.append(quad_name)

            if bool(single_session.get("cancelled", False)):
                cancelled = True
                break

        return {
            "mode": "multi_quad_scan",
            "quadrupoles": quad_names,
            "screens": list(screens),
            "reference_screen": reference_screen if reference_screen is not None else (
                list(screens)[0] if len(list(screens)) > 0 else None),
            "current_A_min": float(current_min),
            "current_A_max": float(current_max),
            "screen_current_ranges": {str(k): [float(v[0]), float(v[1])] for k, v in dict(screen_current_ranges or {}).items()},
            "steps": int(steps),
            "is_quad_scan": bool(int(steps) > 0),
            "nsteps_scan": 1 if int(steps) == 0 else int(steps),
            "nshots": int(nshots),
            "per_quad_sessions": per_quad_sessions,
            "completed_quadrupoles": list(completed_quadrupoles),
            "skipped_quadrupoles": list(skipped_quadrupoles),
            "cancelled": bool(cancelled),
        }

    def _run_single_scan(self, quad_name, screens, current_min, current_max, steps, nshots, screen_current_ranges=None, reference_screen=None, progress_callback=None, resume_states_dir=None):
        screens = list(screens)
        if reference_screen is None:
            reference_screen = screens[0]
        steps_requested = int(steps)
        screens = [reference_screen] + [s for s in screens if s != reference_screen] # so that reference screen is first on the list

        quadrupoles = self.interface.get_quadrupoles([quad_name])
        quad_value_unit = str(quadrupoles.get("value_unit", "1/m"))
        quad_value_label = "current" if quad_value_unit == "A" else "K1L"
        bdes = np.asarray(quadrupoles.get("bdes", []), dtype=float)
        K1L_0 = float(bdes[0])
        current_0 = self.quad_setpoint_to_current(quad_name, K1L_0, quad_value_unit)

        if steps_requested == 0:
            quad_values = np.array([K1L_0], dtype=float)
            values_per_screen = [quad_values.copy() for _ in screens]
            currents_per_screen = [np.array([current_0], dtype=float) for _ in screens]
            rows_per_screen = [np.array([0]) for _ in screens]
        else:
            quad_values, values_per_screen, currents_per_screen, rows_per_screen = self._scan_grid_per_screen(
                quad_name=quad_name, screens=screens, current_min=current_min, current_max=current_max,
                steps=steps_requested, quad_value_unit=quad_value_unit, screen_current_ranges=screen_current_ranges)

        K1L_values = quad_values
        currents_scanned = np.array([self.quad_setpoint_to_current(quad_name, value, quad_value_unit) for value in quad_values], dtype=float)
        nsteps_scan = len(K1L_values)
        nscreens = len(screens)
        deltas = K1L_values / K1L_0 - 1.0 if np.isfinite(K1L_0) and not np.isclose(K1L_0, 0.0) else np.full(nsteps_scan, np.nan)

        sigx_mean = np.full((nsteps_scan, nscreens), np.nan, dtype=float)
        sigy_mean = np.full((nsteps_scan, nscreens), np.nan, dtype=float)
        sigxy_mean = np.full((nsteps_scan, nscreens), np.nan, dtype=float)
        sigx_std = np.full((nsteps_scan, nscreens), np.nan, dtype=float)
        sigy_std = np.full((nsteps_scan, nscreens), np.nan, dtype=float)
        sigxy_std = np.full((nsteps_scan, nscreens), np.nan, dtype=float)
        sigx_shots = np.full((nsteps_scan, nscreens, nshots), np.nan, dtype=float)
        sigy_shots = np.full((nsteps_scan, nscreens, nshots), np.nan, dtype=float)
        sigxy_shots = np.full((nsteps_scan, nscreens, nshots), np.nan, dtype=float)
        x_mean = np.full((nsteps_scan, nscreens), np.nan, dtype=float)
        y_mean = np.full((nsteps_scan, nscreens), np.nan, dtype=float)
        x_std = np.full((nsteps_scan, nscreens), np.nan, dtype=float)
        y_std = np.full((nsteps_scan, nscreens), np.nan, dtype=float)
        x_shots = np.full((nsteps_scan, nscreens, nshots), np.nan, dtype=float)
        y_shots = np.full((nsteps_scan, nscreens, nshots), np.nan, dtype=float)

        scan_steps = []
        images = [[[None for _ in range(nshots)] for _ in range(nscreens)] for _ in range(nsteps_scan)]
        hedges = [[[None for _ in range(nshots)] for _ in range(nscreens)] for _ in range(nsteps_scan)]
        vedges = [[[None for _ in range(nshots)] for _ in range(nscreens)] for _ in range(nsteps_scan)]
        output_dir = resume_states_dir if resume_states_dir and os.path.isdir(resume_states_dir) else self._get_scan_dir(quad_name, steps_requested)
        cancel_requested = False
        self._last_scan_states_dir = output_dir
        self.load_screens_data_database.setText(output_dir)
        total_measurements = int(sum(len(values) for values in values_per_screen))
        completed_measurements = 0

        try:
            for k, screen_name in enumerate(screens):
                while getattr(self, "_scan_pause_requested", False) and not getattr(self, "_scan_stop_requested", False):
                    setattr(self, "_scan_is_paused", True)
                    QApplication.processEvents()
                    time.sleep(0.05)

                setattr(self, "_scan_is_paused", False)
                if getattr(self, "_scan_stop_requested", False):
                    raise KeyboardInterrupt("Scan stopped by user.")
                if getattr(self, "_cancel", False):
                    cancel_requested = True
                    break
                try:
                    for i_local, K1L in enumerate(values_per_screen[k]):
                        i = int(rows_per_screen[k][i_local])
                        while getattr(self, "_scan_pause_requested", False) and not getattr(self, "_scan_stop_requested", False):
                            setattr(self, "_scan_is_paused", True)
                            QApplication.processEvents()
                            time.sleep(0.05)
                        setattr(self, "_scan_is_paused", False)
                        if getattr(self, "_scan_stop_requested", False):
                            raise KeyboardInterrupt("Scan stopped by user.")
                        if getattr(self, "_cancel", False):
                            cancel_requested = True
                            break
                        sx_shots = np.full(nshots, np.nan, dtype=float)
                        sy_shots = np.full(nshots, np.nan, dtype=float)
                        sxy_shots = np.full(nshots, np.nan, dtype=float)
                        dx_shots = np.full(nshots, np.nan, dtype=float)
                        dy_shots = np.full(nshots, np.nan, dtype=float)
                        state_files = []
                        already_measured = self._saved_measurement_files(output_dir, k, i, nshots) if resume_states_dir else None
                        if already_measured is not None:
                            print(f"Resuming: screen {screen_name} step {i} is already on disk, reusing it")
                            self._fill_shots_from_saved_files(already_measured, screen_name, sx_shots, sy_shots, sxy_shots,
                                dx_shots, dy_shots, images, hedges, vedges, i, k)
                            state_files = list(already_measured)
                        else:
                            if steps_requested > 0:
                                print("Before set_quadrupoles")
                                reached = self.interface.set_quadrupoles([quad_name], [float(K1L)])
                                if reached is False:
                                    reached = self.interface.set_quadrupoles([quad_name], [float(K1L)])
                                if reached is False:
                                    raise RuntimeError(
                                        f"{quad_name} did not reach requested {quad_value_label}="
                                        f"{K1L:.6g} {quad_value_unit} (step {i}) after retry")
                                print("After set_quadrupoles")
                            print("before calling get_quadrupoles")
                            quad_data = self.interface.get_quadrupoles([quad_name])
                            print("after calling get_quadrupoles")
                            for j in range(nshots):
                                while getattr(self, "_scan_pause_requested", False) and not getattr(self, "_scan_stop_requested", False):
                                    setattr(self, "_scan_is_paused", True)
                                    QApplication.processEvents()
                                    time.sleep(0.05)
                                setattr(self, "_scan_is_paused", False)
                                if getattr(self, "_scan_stop_requested", False):
                                    raise KeyboardInterrupt("Scan stopped by user.")
                                if getattr(self, "_cancel", False):
                                    cancel_requested = True
                                    break
                                print("before calling get_screens")
                                screens_data = self.interface.get_screens([screen_name])
                                sigma_scale = 1.0
                                sigxy_scale = 1.0
                                if self._get_interface_units() == "um":
                                    sigma_scale = 1.0 / 1000.0
                                    sigxy_scale = 1.0 / 1000000.0
                                print("after calling get_screens")
                                idx_map = {name: idx for idx, name in enumerate(screens_data["names"])}
                                idx = idx_map.get(screen_name)
                                if idx is not None:
                                    sx_shots[j] = float(screens_data["sigx"][idx]) * sigma_scale
                                    sy_shots[j] = float(screens_data["sigy"][idx]) * sigma_scale
                                    images[i][k][j]=np.asarray(screens_data["images"][idx])
                                    if idx < len(screens_data.get("hedges", [])):
                                        hedges[i][k][j] = np.asarray(screens_data["hedges"][idx], dtype=float)
                                    if idx < len(screens_data.get("vedges", [])):
                                        vedges[i][k][j] = np.asarray(screens_data["vedges"][idx], dtype=float)
                                    if "sigxy" in screens_data:
                                        sxy_shots[j] = float(screens_data["sigxy"][idx]) * sigxy_scale
                                    if "x" in screens_data:
                                        dx_shots[j] = float(screens_data["x"][idx]) * sigma_scale
                                    if "y" in screens_data:
                                        dy_shots[j] = float(screens_data["y"][idx]) * sigma_scale
                                state_for_scan = State(sextupoles=None, correctors=None, bpms=None,
                                    icts=None, sequence=self.interface.get_sequence(), hcorrectors_names=None,
                                    vcorrectors_names=None, screens=screens_data, quadrupoles=quad_data)

                                state_filename = os.path.join(output_dir, f"screen_{k:04d}_step_{i:04d}_shot_{j:04d}.pkl")
                                state_for_scan.save(filename=state_filename)
                                state_files.append(state_filename)

                        if state_files:
                            sigx_mean[i, k] = np.nanmedian(sx_shots)
                            sigy_mean[i, k] = np.nanmedian(sy_shots)
                            sigxy_mean[i, k] = np.nanmean(sxy_shots)
                            sigx_std[i, k] = np.nanstd(sx_shots)
                            sigy_std[i, k] = np.nanstd(sy_shots)
                            sigxy_std[i, k] = np.nanstd(sxy_shots)
                            sigx_shots[i, k, :] = sx_shots
                            sigy_shots[i, k, :] = sy_shots
                            sigxy_shots[i, k, :] = sxy_shots
                            x_mean[i, k] = np.nanmean(dx_shots)
                            y_mean[i, k] = np.nanmean(dy_shots)
                            x_std[i, k] = np.nanstd(dx_shots)
                            y_std[i, k] = np.nanstd(dy_shots)
                            x_shots[i, k, :] = dx_shots
                            y_shots[i, k, :] = dy_shots

                        existing_step = next((step for step in scan_steps if int(step.get("step_index", -1)) == int(i)), None)
                        if existing_step is None:
                            existing_step = {
                                "step_index": int(i),
                                "delta": float(deltas[i]),
                                "K1L": float(K1L),
                                "quad_value": float(K1L),
                                "current_A": float(currents_per_screen[k][i_local]),
                                "state_files": [],
                            }
                            scan_steps.append(existing_step)
                        existing_step["state_files"].extend(state_files)
                        completed_measurements += 1

                        session_partial = {
                            "current_A_min": float(current_min),
                            "current_A_max": float(current_max),
                            "screen_current_ranges": {str(name): [float(currents[0]), float(currents[-1])] for name, currents in zip(screens, currents_per_screen)},
                            "steps": int(steps_requested),
                            "is_quad_scan": bool(steps_requested>0),
                            "nshots": int(nshots),
                            "quad_name": quad_name,
                            "quadrupoles": [quad_name],
                            "screens": screens,
                            "reference_screen": reference_screen,
                            "quad_value_unit": quad_value_unit,
                            "K1L_0": float(K1L_0),
                            "current_0": float(current_0),
                            "sigx_mean": sigx_mean.tolist(),
                            "sigy_mean": sigy_mean.tolist(),
                            "sigxy_mean": sigxy_mean.tolist(),
                            "sigx_std": sigx_std.tolist(),
                            "sigy_std": sigy_std.tolist(),
                            "sigxy_std": sigxy_std.tolist(),
                            "sigx_shots": sigx_shots.tolist(),
                            "sigy_shots": sigy_shots.tolist(),
                            "sigxy_shots": sigxy_shots.tolist(),
                            "x_mean": x_mean.tolist(),
                            "y_mean": y_mean.tolist(),
                            "x_std": x_std.tolist(),
                            "y_std": y_std.tolist(),
                            "x_shots": x_shots.tolist(),
                            "y_shots": y_shots.tolist(),
                            "deltas": deltas.tolist(),
                            "K1L_values": K1L_values.tolist(),
                            "current_values": currents_scanned.tolist(),
                            "scan_steps": scan_steps,
                            "states_dir": output_dir,
                            "cancelled": bool(cancel_requested),
                            "current_screen": screen_name,
                            "current_screen_index": int(k),
                            "nsteps_scan": int(nsteps_scan),
                            "images": images,
                            "hedges": hedges,
                            "vedges": vedges,
                            "sigma_unit": "mm",
                        }

                        if progress_callback is not None:
                            progress_callback(session_partial, completed_measurements, total_measurements)
                        if cancel_requested:
                            break
                    if cancel_requested:
                        break

                finally:
                    extract_screen = getattr(self.interface, "extract_screen", None)
                    if callable(extract_screen):
                        extract_screen(screen_name)
        finally:
            if steps_requested > 0 and np.isfinite(K1L_0):
                try:
                    restored = self.interface.set_quadrupoles([quad_name], [float(K1L_0)])
                    if restored is False:
                        print(
                            f"{quad_name} did not confirm returning to its original "
                            f"{quad_value_label}={K1L_0:.6g} {quad_value_unit} after the scan."
                        )
                except Exception as exc:
                    print(
                        f"Failed to restore {quad_name} to its original "
                        f"{quad_value_label}={K1L_0:.6g} {quad_value_unit}: {exc}"
                    )

        session = {
            "current_A_min": float(current_min),
            "current_A_max": float(current_max),
            "screen_current_ranges": {str(name): [float(currents[0]), float(currents[-1])] for name, currents in zip(screens, currents_per_screen)},
            "steps": int(steps_requested),
            "is_quad_scan": bool(steps_requested > 0),
            "nshots": int(nshots),
            "quad_name": quad_name,
            "quadrupoles": [quad_name],
            "screens": screens,
            "reference_screen": reference_screen,
            "quad_value_unit": quad_value_unit,
            "K1L_0": float(K1L_0),
            "current_0": float(current_0),
            "sigx_mean": sigx_mean.tolist(),
            "sigy_mean": sigy_mean.tolist(),
            "sigxy_mean": sigxy_mean.tolist(),
            "sigx_std": sigx_std.tolist(),
            "sigy_std": sigy_std.tolist(),
            "sigxy_std": sigxy_std.tolist(),
            "sigx_shots": sigx_shots.tolist(),
            "sigy_shots": sigy_shots.tolist(),
            "sigxy_shots": sigxy_shots.tolist(),
            "x_mean": x_mean.tolist(),
            "y_mean": y_mean.tolist(),
            "x_std": x_std.tolist(),
            "y_std": y_std.tolist(),
            "x_shots": x_shots.tolist(),
            "y_shots": y_shots.tolist(),
            "deltas": deltas.tolist(),
            "K1L_values": K1L_values.tolist(),
            "current_values": currents_scanned.tolist(),
            "scan_steps": scan_steps,
            "states_dir": output_dir,
            "cancelled": bool(cancel_requested),
            "nsteps_scan": int(nsteps_scan), # number of measurements at screens (even if steps=0, nsteps_scan = 1)
            "images": images,
            "hedges": hedges,
            "vedges": vedges,
            "sigma_unit": "mm",
            "nscreens": len(screens),
        }

        set_default_quad_strength_bounds = getattr(self, "_set_default_quad_strength_bounds_from_session", None)
        if callable(set_default_quad_strength_bounds):
            set_default_quad_strength_bounds(session)
        self.save_emittance_measurement_session(session)
        return session

    @staticmethod
    def _saved_measurement_files(output_dir, screen_index, step_index, nshots):
        if not output_dir or not os.path.isdir(output_dir):
            return None
        paths = [os.path.join(output_dir, f"screen_{screen_index:04d}_step_{step_index:04d}_shot_{shot:04d}.pkl")
                 for shot in range(int(nshots))]
        return paths if all(os.path.isfile(path) for path in paths) else None

    def _fill_shots_from_saved_files(self, paths, screen_name, sx_shots, sy_shots, sxy_shots, dx_shots, dy_shots, images, hedges, vedges, step_index, screen_index):
        sigma_scale = 1.0
        sigxy_scale = 1.0
        if self._get_interface_units() == "um":
            sigma_scale = 1.0 / 1000.0
            sigxy_scale = 1.0 / 1000000.0
        for shot, path in enumerate(paths):
            screens_data = State(filename=path).get_screens()
            names = list(screens_data.get("names", []))
            idx = names.index(screen_name) if screen_name in names else (0 if names else None)
            if idx is None:
                continue
            sx_shots[shot] = float(np.ravel(screens_data["sigx"])[idx]) * sigma_scale
            sy_shots[shot] = float(np.ravel(screens_data["sigy"])[idx]) * sigma_scale
            saved_images = screens_data.get("images", [])
            if idx < len(saved_images):
                images[step_index][screen_index][shot] = np.asarray(saved_images[idx])
            if idx < len(screens_data.get("hedges", [])):
                hedges[step_index][screen_index][shot] = np.asarray(screens_data["hedges"][idx], dtype=float)
            if idx < len(screens_data.get("vedges", [])):
                vedges[step_index][screen_index][shot] = np.asarray(screens_data["vedges"][idx], dtype=float)
            if "sigxy" in screens_data:
                sxy_shots[shot] = float(np.ravel(screens_data["sigxy"])[idx]) * sigxy_scale
            if "x" in screens_data:
                dx_shots[shot] = float(np.ravel(screens_data["x"])[idx]) * sigma_scale
            if "y" in screens_data:
                dy_shots[shot] = float(np.ravel(screens_data["y"])[idx]) * sigma_scale

    def _get_scan_dir(self,quad_name, steps_requested): # saves state files for each quadrupole
        time_str = datetime.now().strftime("%Y%m%d_%H%M%S")
        base_dir=getattr(self,"dir_name",None)
        if not base_dir:
            base_dir=os.path.join(os.getcwd(),"emittance_measurement_session")
        base_dir = os.path.expanduser(os.path.expandvars(base_dir))

        if int(steps_requested) == 0:
            scan_dir = os.path.join(base_dir, f"screens_data_{time_str}")
        else:
            scan_dir=os.path.join(base_dir,f"states_{quad_name}_{time_str}")
        os.makedirs(scan_dir,exist_ok=True) # if exists, no error while trying to create a folder
        return scan_dir
