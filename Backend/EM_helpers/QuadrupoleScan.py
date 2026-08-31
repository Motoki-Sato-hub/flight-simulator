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

    def run_scan(self, quad_name, screens, delta_min, delta_max, steps, nshots, reference_screen=None, progress_callback=None):
        steps = int(steps)
        if isinstance(quad_name, str):
            quad_names = [quad_name]
        else:
            quad_names = list(quad_name or [])

        if len(quad_names) == 0:
            raise ValueError("At least one quadrupole must be provided")

        self.dir_name = self._new_scan_session_dir(quad_names=quad_names, is_quad_scan=(steps > 0))
        self.session_directory.setText(self.dir_name)

        if steps == 0:
            if len(quad_names) != 1:
                raise ValueError("For Linear Response, choose exactly one quadrupole.")
            return self._run_single_scan(quad_name=quad_names[0], screens=screens,
                delta_min=0.0, delta_max=0.0, steps=0, nshots=nshots,
                reference_screen=reference_screen, progress_callback=progress_callback)

        if len(quad_names) == 1:
            print(f"{quad_names[0]} is going to be scanned")
            return self._run_single_scan(quad_name=quad_names[0], screens=screens, delta_min=delta_min, delta_max=delta_max, steps=steps, nshots=nshots, reference_screen=reference_screen, progress_callback=progress_callback)

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
                single_session = self._run_single_scan(quad_name=quad_name, screens=screens, delta_min=delta_min, delta_max=delta_max,
                    steps=steps, nshots=nshots, reference_screen=reference_screen, progress_callback=_wrapped_progress)

            except ValueError as e:
                msg = str(e)
                if "zero K1L_0" in msg:
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
            "delta_min": float(delta_min),
            "delta_max": float(delta_max),
            "steps": int(steps),
            "is_quad_scan": bool(int(steps) > 0),
            "nsteps_scan": 1 if int(steps) == 0 else int(steps),
            "nshots": int(nshots),
            "per_quad_sessions": per_quad_sessions,
            "completed_quadrupoles": list(completed_quadrupoles),
            "skipped_quadrupoles": list(skipped_quadrupoles),
            "cancelled": bool(cancelled),
        }

    def _run_single_scan(self, quad_name, screens, delta_min, delta_max, steps, nshots, reference_screen=None, progress_callback=None):
        screens = list(screens)
        if reference_screen is None:
            reference_screen = screens[0]
        steps_requested = int(steps)
        if steps_requested > 0 and delta_max <= delta_min:
            raise ValueError("delta_max must be larger than delta_min")
        screens = [reference_screen] + [s for s in screens if s != reference_screen] # so that reference screen is first on the list

        quad_names = list(getattr(self.interface, "quadrupoles", []))
        quadrupoles = self.interface.get_quadrupoles([quad_name])
        bdes = np.asarray(quadrupoles.get("bdes", []), dtype=float)
        K1L_0 = float(bdes[0])

        if steps_requested > 0 and np.isclose(K1L_0, 0.0):
            raise ValueError("The quadrupole has zero K1L_0. You should choose another one.")

        if steps_requested == 0:
            deltas = np.array([0.0], dtype=float)
            K1L_values = np.array([K1L_0], dtype=float)
        else:
            deltas = np.linspace(float(delta_min), float(delta_max), steps_requested)
            K1L_values = K1L_0 * (1 + deltas)
        nsteps_scan = len(K1L_values)
        nscreens = len(screens)

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
        output_dir = self._get_scan_dir(quad_name, steps_requested)
        cancel_requested = False
        self.load_screens_data_database.setText(output_dir)

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
                    for i, K1L in enumerate(K1L_values):
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
                        if steps_requested > 0:
                            print("Before set_quadrupoles")
                            self.interface.set_quadrupoles([quad_name], [float(K1L)])
                            #time.sleep(1)
                            print("After set_quadrupoles")
                        sx_shots = np.full(nshots, np.nan, dtype=float)
                        sy_shots = np.full(nshots, np.nan, dtype=float)
                        sxy_shots = np.full(nshots, np.nan, dtype=float)
                        dx_shots = np.full(nshots, np.nan, dtype=float)
                        dy_shots = np.full(nshots, np.nan, dtype=float)
                        state_files = []
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
                            sigx_mean[i, k] = np.nanmean(sx_shots)
                            sigy_mean[i, k] = np.nanmean(sy_shots)
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
                                "state_files": [],
                            }
                            scan_steps.append(existing_step)
                        existing_step["state_files"].extend(state_files)

                        session_partial = {
                            "delta_min": float(delta_min),
                            "delta_max": float(delta_max),
                            "steps": int(steps_requested),
                            "is_quad_scan": bool(steps_requested>0),
                            "nshots": int(nshots),
                            "quad_name": quad_name,
                            "quadrupoles": [quad_name],
                            "screens": screens,
                            "reference_screen": reference_screen,
                            "K1L_0": float(K1L_0),
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
                            completed = k * len(K1L_values) + i + 1
                            total = len(screens) * len(K1L_values)
                            progress_callback(session_partial, completed, total)
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
                self.interface.set_quadrupoles([quad_name], [float(K1L_0)])

        session = {
            "delta_min": float(delta_min),
            "delta_max": float(delta_max),
            "steps": int(steps_requested),
            "is_quad_scan": bool(steps_requested > 0),
            "nshots": int(nshots),
            "quad_name": quad_name,
            "quadrupoles": [quad_name],
            "screens": screens,
            "reference_screen": reference_screen,
            "K1L_0": float(K1L_0),
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
