import os, pickle, re, matplotlib, glob
import numpy as np
matplotlib.use("QtAgg")
from Backend.State import State
from Knobs.jitter_subtraction import apply_jitter_subtraction

'''
Mutual response matrix calculation algorithm.
'''

class AdaptiveResponseMatrix:
    def __init__(self, R0, forgetting=0.995, rank=3, p0_scale=100.0):
        self.R0 = np.array(R0, dtype=float)
        self.delta = np.zeros_like(self.R0)
        self.P = None
        self.p0_scale = p0_scale
        self.forgetting = forgetting
        self.rank = rank

    def update(self, u, dx):
        u = np.asarray(u, dtype=float).ravel() # vector with corrector kicks
        dx = np.asarray(dx, dtype=float).ravel() # change of orbits on BPMs [on x, on y]
        if self.P is None:
            self.P = self.p0_scale / max(u @ u, 1e-12) * np.eye(len(u)) # how much RLS knows the matrix, if its big, it moves the R faster
            # self.R0 + self.delta is a currently used response matrix
        residual = dx - (self.R0 + self.delta) @ u # what bpms saw - what new R predicted, the residual is used to correct delta
        Pu = self.P @ u # which fragments are not yet fully explored
        denom = self.forgetting + u @ Pu # normalizes the size of update
        self.P = (self.P - np.outer(Pu, Pu) / denom) / self.forgetting
        self.delta = self.delta + np.outer(residual, Pu) / denom # change in R calculated by RLS
        U, S, Vt = np.linalg.svd(self.delta, full_matrices=False) # strongest known directions
        r = min(self.rank, S.size) # 3 or biggest S
        self.delta = (U[:, :r] * S[:r]) @ Vt[:r] # update calculated by RLS

    @property
    def R(self):
        return self.R0 + self.delta

class ResponseMatrix_DFS_WFS():

    def _compute_response_matrix_from_directory(self, directory, correctors, bpms, triangular=False, actuator_mode="correctors", rcond=1e-3):
        info=self._find_useful_files(directory)
        if not info["ok"]:
            raise RuntimeError(f"Could not find any valid DATA pairs in {directory}")

        return self._compute_response_matrix(pairs=info["pairs"], correctors=correctors, bpms=bpms, triangular=triangular, actuator_mode=actuator_mode, rcond=rcond)

    def _find_useful_files(self, directory):
        datafiles=sorted(glob.glob(os.path.join(directory, 'DATA*.pkl')))
        pairs=[]
        pair_re = re.compile(r"DATA_(.+)_(p|m)(\d+)\.pkl$")
        for fp in datafiles:
            basename=os.path.basename(fp)
            match = pair_re.search(basename)
            if not match or match.group(2) != "p":
                continue
            tag = match.group(1)
            shot_index = match.group(3)
            fm = os.path.join(os.path.dirname(fp), f"DATA_{tag}_m{shot_index}.pkl")
            if os.path.exists(fm):
                pairs.append((fp, fm, tag))

        return {"ok":bool(pairs), "dir":directory, "pairs":pairs}

    def _compute_response_matrix(self, pairs, correctors, bpms, triangular=False, actuator_mode="correctors", rcond=1e-3):
        if not hasattr(self, 'sequence'):
            file = pairs[0][0]
            S = State(filename=file)
            self.sequence = S.get_sequence()

        if actuator_mode == "quadrupole_movers":
            return self._compute_qm_response_matrix(pairs=pairs, qcorrs=correctors, bpms=bpms, triangular=triangular)
        else:
            # hcorrs = [string for string in correctors if self._is_h_corrector(string)]
            # vcorrs = [string for string in correctors if self._is_v_corrector(string)]
            file = pairs[0][0]
            S = State(filename=file)
            hcorrs = [c for c in correctors if c in S.hcorrectors_names]
            vcorrs = [c for c in correctors if c in S.vcorrectors_names]

        # Pick all correctors preceding the last bpm
        hcorrs = [corr for corr in hcorrs if self.sequence.index(corr) < self.sequence.index(bpms[-1])]
        vcorrs = [corr for corr in vcorrs if self.sequence.index(corr) < self.sequence.index(bpms[-1])]

        # Pick all bpms following the first corrector
        if hcorrs:
            bpms = [bpm for bpm in bpms if self.sequence.index(bpm) > self.sequence.index(hcorrs[0])]
        if vcorrs:
            bpms = [bpm for bpm in bpms if self.sequence.index(bpm) > self.sequence.index(vcorrs[0])]

        # Read all orbits
        Bx = np.empty((0, len(bpms)))
        By = np.empty((0, len(bpms)))
        Cx = np.empty((0, len(hcorrs)))
        Cy = np.empty((0, len(vcorrs)))
        excited_hcorrs = set()
        excited_vcorrs = set()

        for pair in pairs:
            if len(pair) == 3:
                fp, fm, tag = pair
            else:
                fp, fm = pair
                tag = ""
            Sp = State(filename=fp)
            Sm = State(filename=fm)

            # jitter_model = getattr(self, "jitter_model", None)
            # if jitter_model is not None:
            #     Sp.bpms = apply_jitter_subtraction(Sp.get_bpms(), jitter_model)
            #     Sm.bpms = apply_jitter_subtraction(Sm.get_bpms(), jitter_model)

            Op = Sp.get_orbit(bpms)
            Om = Sm.get_orbit(bpms)

            all_not_finite = not np.any(np.isfinite(Op['x']))
            all_not_finite |= not np.any(np.isfinite(Om['x']))
            all_not_finite |= not np.any(np.isfinite(Op['y']))
            all_not_finite |= not np.any(np.isfinite(Om['y']))

            if all_not_finite:
                print(f"Skipping all-NaN files: {os.path.basename(fp)} / {os.path.basename(fm)}")
                continue

            if actuator_mode == "quadrupole_movers":
                if str(tag).endswith("_x"):
                    qname = str(tag)[:-2]
                    Cx_p = np.zeros(len(hcorrs), dtype=float)
                    Cx_m = np.zeros(len(hcorrs), dtype=float)
                    Cy_p = np.zeros(len(vcorrs), dtype=float)
                    Cy_m = np.zeros(len(vcorrs), dtype=float)
                    if qname in hcorrs:
                        qidx = hcorrs.index(qname)
                        qxp = Sp.get_quadrupoles([qname])
                        qxm = Sm.get_quadrupoles([qname])
                        Cx_p[qidx] = np.asarray(qxp.get('xact', qxp.get('xdes', [0.0])), dtype=float)[0]
                        Cx_m[qidx] = np.asarray(qxm.get('xact', qxm.get('xdes', [0.0])), dtype=float)[0]
                    else:
                        continue
                elif str(tag).endswith("_y"):
                    qname = str(tag)[:-2]
                    Cx_p = np.zeros(len(hcorrs), dtype=float)
                    Cx_m = np.zeros(len(hcorrs), dtype=float)
                    Cy_p = np.zeros(len(vcorrs), dtype=float)
                    Cy_m = np.zeros(len(vcorrs), dtype=float)
                    if qname in vcorrs:
                        qidx = vcorrs.index(qname)
                        qyp = Sp.get_quadrupoles([qname])
                        qym = Sm.get_quadrupoles([qname])
                        Cy_p[qidx] = np.asarray(qyp.get('yact', qyp.get('ydes', [0.0])), dtype=float)[0]
                        Cy_m[qidx] = np.asarray(qym.get('yact', qym.get('ydes', [0.0])), dtype=float)[0]
                    else:
                        continue
                else:
                    continue
            else:
                Cx_p = Sp.get_correctors(hcorrs)['bact']
                Cy_p = Sp.get_correctors(vcorrs)['bact']
                Cx_m = Sm.get_correctors(hcorrs)['bact']
                Cy_m = Sm.get_correctors(vcorrs)['bact']

                if tag in hcorrs:
                    index = hcorrs.index(tag)
                    requested = Sp.get_correctors([tag])['bdes'][0] - Sm.get_correctors([tag])['bdes'][0]
                    measured = Cx_p[index] - Cx_m[index]
                elif tag in vcorrs:
                    index = vcorrs.index(tag)
                    requested = Sp.get_correctors([tag])['bdes'][0] - Sm.get_correctors([tag])['bdes'][0]
                    measured = Cy_p[index] - Cy_m[index]
                else:
                    requested = measured = np.nan

                if np.isfinite(requested) and abs(requested) > 1e-12 and abs(measured) < 0.5 * abs(requested):
                    print(f"Skipping unexecuted excitation {tag}: Δbdes={requested:.6g}, Δbact={measured:.6g}")
                    continue
                if tag in hcorrs:
                    excited_hcorrs.add(tag)
                elif tag in vcorrs:
                    excited_vcorrs.add(tag)


            Bx = np.vstack((Bx, Op['x']))
            Bx = np.vstack((Bx, Om['x']))
            By = np.vstack((By, Op['y']))
            By = np.vstack((By, Om['y']))
            Cx = np.vstack((Cx, Cx_p))
            Cx = np.vstack((Cx, Cx_m))
            Cy = np.vstack((Cy, Cy_p))
            Cy = np.vstack((Cy, Cy_m))

        Bx, By = Bx.astype(float), By.astype(float)
        Cx, Cy = Cx.astype(float), Cy.astype(float)

        # Drop unavailable correctors and incomplete shots before fitting.
        # A corrector needs at least one successfully executed excitation;
        # otherwise its constant readback is collinear with the offset term.
        hcol_mask = np.array([corr in excited_hcorrs for corr in hcorrs]) & np.any(np.isfinite(Cx), axis=0)
        vcol_mask = np.array([corr in excited_vcorrs for corr in vcorrs]) & np.any(np.isfinite(Cy), axis=0)
        hcorrs = [corr for corr, keep in zip(hcorrs, hcol_mask) if keep]
        vcorrs = [corr for corr, keep in zip(vcorrs, vcol_mask) if keep]
        Cx, Cy = Cx[:, hcol_mask], Cy[:, vcol_mask]

        row_mask = np.all(np.isfinite(Cx), axis=1) & np.all(np.isfinite(Cy), axis=1)
        if not np.any(row_mask):
            raise RuntimeError("No complete corrector readbacks available for response-matrix fitting")
        Bx, By, Cx, Cy = Bx[row_mask], By[row_mask], Cx[row_mask], Cy[row_mask]

        B_mask = np.all(np.isfinite(Bx), axis=0) & np.all(np.isfinite(By), axis=0)
        if not np.any(B_mask):
            raise RuntimeError("No BPM has complete orbit data for response-matrix fitting")

        Cx = np.hstack((Cx, np.ones((Cx.shape[0], 1))))
        Cy = np.hstack((Cy, np.ones((Cy.shape[0], 1))))

        def lstsq(C, B):
            return np.transpose(np.linalg.lstsq(C, B[:, B_mask], rcond=None)[0])

        Rxx_ = lstsq(Cx, Bx)
        Rxy_ = lstsq(Cy, Bx)
        Ryx_ = lstsq(Cx, By)
        Ryy_ = lstsq(Cy, By)

        # Response matrices: remove offset column
        Rxx_ = Rxx_[:, :-1]
        Rxy_ = Rxy_[:, :-1]
        Ryx_ = Ryx_[:, :-1]
        Ryy_ = Ryy_[:, :-1]

        # Restore nan columns
        k = B_mask.size

        Rxx = np.full((k,Rxx_.shape[1]), np.nan)
        Rxy = np.full((k,Rxy_.shape[1]), np.nan)
        Ryx = np.full((k,Ryx_.shape[1]), np.nan)
        Ryy = np.full((k,Ryy_.shape[1]), np.nan)

        Rxx[B_mask,:] = Rxx_
        Rxy[B_mask,:] = Rxy_
        Ryx[B_mask,:] = Ryx_
        Ryy[B_mask,:] = Ryy_

        # Reference trajectory
        Bx = np.mean(Bx, axis=0).reshape(-1, 1)
        By = np.mean(By, axis=0).reshape(-1, 1)

        # Zero the response of all bpms preceeding the correctors
        if triangular:
            for corr in hcorrs:
                bpm_indexes = [bpms.index(bpm) for bpm in bpms if self.sequence.index(bpm) < self.sequence.index(corr)]
                Rxx[bpm_indexes, hcorrs.index(corr)] = 0
                Ryx[bpm_indexes, hcorrs.index(corr)] = 0

                if actuator_mode == "quadrupole_movers":
                    forbidden = f"M{corr}"
                    bpm_indexes = [bpms.index(bpm) for bpm in bpms if forbidden in bpm]
                    Rxx[bpm_indexes, hcorrs.index(corr)] = 0
                    Ryx[bpm_indexes, hcorrs.index(corr)] = 0

            for corr in vcorrs:
                bpm_indexes = [bpms.index(bpm) for bpm in bpms if self.sequence.index(bpm) < self.sequence.index(corr)]
                Rxy[bpm_indexes, vcorrs.index(corr)] = 0
                Ryy[bpm_indexes, vcorrs.index(corr)] = 0
                if actuator_mode == "quadrupole_movers":
                    forbidden = f"M{corr}"
                    bpm_indexes = [bpms.index(bpm) for bpm in bpms if forbidden in bpm]
                    Rxy[bpm_indexes, vcorrs.index(corr)] = 0
                    Ryy[bpm_indexes, vcorrs.index(corr)] = 0

        return Rxx, Ryy, Rxy, Ryx, Bx, By, hcorrs, vcorrs, bpms

    def _get_data_from_loaded_directories(self, selected_bpms, selected_corrs, _force_triangular=False):

        info_traj = self._data_dirs.get("traj")
        info_dfs = self._data_dirs.get("dfs")
        info_wfs = self._data_dirs.get("wfs")

        w1, w2, w3, rcond, iters, gain, beta, transmission_threshold = self._read_params()

        if not (info_traj and info_traj["ok"]):
            raise RuntimeError("Please select a trajectory data directory")
        if w2 > 0 and not (info_dfs and info_dfs["ok"]):
            raise RuntimeError("Please select a dispersion data directory")
        if w3 > 0 and not (info_wfs and info_wfs["ok"]):
            raise RuntimeError("Please select a wakefield data directory")

        triangular = bool(self._force_triangular() or _force_triangular)

        R0xx, R0yy, R0xy, R0yx, B0x, B0y, hcorrs0, vcorrs0, bpms0 = self._compute_response_matrix(
            pairs=info_traj["pairs"], correctors=selected_corrs, bpms=selected_bpms, triangular=triangular
        )

        if w2 > 0:
            R1xx, R1yy, R1xy, R1yx, B1x, B1y, hcorrs1, vcorrs1, bpms1 = self._compute_response_matrix(
                pairs=info_dfs["pairs"], correctors=selected_corrs, bpms=selected_bpms, triangular=triangular
            )
        else:
            R1xx, R1yy, R1xy, R1yx = R0xx.copy(), R0yy.copy(), R0xy.copy(), R0yx.copy()
            B1x, B1y = B0x.copy(), B0y.copy()
            hcorrs1, vcorrs1, bpms1 = list(hcorrs0), list(vcorrs0), list(bpms0)

        if w3 > 0:
            R2xx, R2yy, R2xy, R2yx, B2x, B2y, hcorrs2, vcorrs2, bpms2 = self._compute_response_matrix(
                pairs=info_wfs["pairs"], correctors=selected_corrs, bpms=selected_bpms, triangular=triangular
            )
        else:
            R2xx, R2yy, R2xy, R2yx = R0xx.copy(), R0yy.copy(), R0xy.copy(), R0yx.copy()
            B2x, B2y = B0x.copy(), B0y.copy()
            hcorrs2, vcorrs2, bpms2 = list(hcorrs0), list(vcorrs0), list(bpms0)

        return (
            R0xx, R0yy, R0xy, R0yx, B0x, B0y,
            R1xx, R1yy, R1xy, R1yx, B1x, B1y,
            R2xx, R2yy, R2xy, R2yx, B2x, B2y,
            hcorrs0, vcorrs0, hcorrs1, vcorrs1, hcorrs2, vcorrs2,
            bpms0, bpms1, bpms2,
        )
    def _creating_response_matrices(self, selected_corrs = None, selected_bpms = None):

        if selected_corrs is None or selected_bpms is None:
            corrs, bpms = self._get_selection()
        else:
            corrs = list(selected_corrs)
            bpms = list(selected_bpms)


        w1, w2, w3, rcond, iters, gain, beta, transmission_threshold = self._read_params()
        wgt_orb, wgt_dfs, wgt_wfs = w1, w2, w3

        #corrs, bpms = self._get_selection()

        (
            R0xx, R0yy, R0xy, R0yx, B0x, B0y,
            R1xx, R1yy, R1xy, R1yx, B1x, B1y,
            R2xx, R2yy, R2xy, R2yx, B2x, B2y,
            hcorrs0, vcorrs0, hcorrs1, vcorrs1, hcorrs2, vcorrs2,
            bpms0, bpms1, bpms2,
        ) = self._get_data_from_loaded_directories(
            selected_corrs=corrs, selected_bpms=bpms, _force_triangular=self._force_triangular()
        )
        def _handle_missing_corrector(R,corrs,corrs_ref):
            if corrs== corrs_ref:
                return R
            else:
                index={c: i for i, c in enumerate(corrs)}
                final_matrix=np.zeros((R.shape[0],len(corrs_ref)), dtype=float)
                for j,c in enumerate(corrs_ref):
                    i=index.get(c)
                    if i is not None:
                        final_matrix[:,j] = R[:,i]
                return final_matrix

        def _handle_missing_bpm(Rxx,Ryy,Rxy,Ryx,Bx,By,bpms,bpms_ref):
            if bpms== bpms_ref:
                return Rxx,Ryy,Rxy,Ryx,Bx,By
            index = {b: i for i, b in enumerate(bpms)}
            bpm_indexes=[index[b] for b in bpms_ref]

            return (
                Rxx[bpm_indexes, :],
                Ryy[bpm_indexes, :],
                Rxy[bpm_indexes, :],
                Ryx[bpm_indexes, :],
                Bx[bpm_indexes, :],
                By[bpm_indexes, :],
            )

        #intersection
        hcorrs=[c for c in hcorrs0 if (c in hcorrs1 and c in hcorrs2)]
        vcorrs=[c for c in vcorrs0 if (c in vcorrs1 and c in vcorrs2)]
        bpms=[b for b in bpms0 if (b in bpms1 and b in bpms2)]

        R0xx=_handle_missing_corrector(R0xx, hcorrs0,hcorrs)
        R0yy=_handle_missing_corrector(R0yy, vcorrs0,vcorrs)
        R0xy=_handle_missing_corrector(R0xy, vcorrs0,vcorrs)
        R0yx=_handle_missing_corrector(R0yx, hcorrs0,hcorrs)
        R1xx=_handle_missing_corrector(R1xx, hcorrs1,hcorrs)
        R1yy=_handle_missing_corrector(R1yy, vcorrs1,vcorrs)
        R1xy=_handle_missing_corrector(R1xy, vcorrs1,vcorrs)
        R1yx=_handle_missing_corrector(R1yx, hcorrs1,hcorrs)
        R2xx=_handle_missing_corrector(R2xx, hcorrs2,hcorrs)
        R2yy=_handle_missing_corrector(R2yy, vcorrs2,vcorrs)
        R2xy=_handle_missing_corrector(R2xy, vcorrs2,vcorrs)
        R2yx=_handle_missing_corrector(R2yx, hcorrs2,hcorrs)

        R0xx, R0yy, R0xy, R0yx, B0x, B0y = _handle_missing_bpm(R0xx, R0yy, R0xy, R0yx, B0x, B0y, bpms0, bpms)
        R1xx, R1yy, R1xy, R1yx, B1x, B1y = _handle_missing_bpm(R1xx, R1yy, R1xy, R1yx, B1x, B1y, bpms1, bpms)
        R2xx, R2yy, R2xy, R2yx, B2x, B2y = _handle_missing_bpm(R2xx, R2yy, R2xy, R2yx, B2x, B2y, bpms2, bpms)

        Axx=[]
        Ayy=[]
        Axy=[]
        Ayx=[]

        if wgt_orb > 0:
            Axx.append(wgt_orb * R0xx)
            Ayy.append(wgt_orb * R0yy)
            Axy.append(wgt_orb * R0xy)
            Ayx.append(wgt_orb * R0yx)

        if wgt_dfs > 0:
            Axx.append(wgt_dfs * (R1xx - R0xx))
            Ayy.append(wgt_dfs * (R1yy - R0yy))
            Axy.append(wgt_dfs * (R1xy - R0xy))
            Ayx.append(wgt_dfs * (R1yx - R0yx))

        if wgt_wfs > 0:
            Axx.append(wgt_wfs * (R2xx - R0xx))
            Ayy.append(wgt_wfs * (R2yy - R0yy))
            Axy.append(wgt_wfs * (R2xy - R0xy))
            Ayx.append(wgt_wfs * (R2yx - R0yx))

        Axx = np.vstack(Axx)
        Ayy = np.vstack(Ayy)
        Axy = np.vstack(Axy)
        Ayx = np.vstack(Ayx)

        return Axx, Ayy, Axy, Ayx, B0x, B0y, hcorrs, vcorrs, bpms

    def _creating_qm_response_matrices(self, selected_corrs, selected_bpms, triangular=True):
        info_traj = self._data_dirs.get("traj")
        if not (info_traj and info_traj.get("ok") and info_traj.get("pairs")):
            return None

        Rxx, Ryy, Rxy, Ryx, _B0x, _B0y, hcorrs, vcorrs, bpms_common = self._compute_response_matrix(
            pairs=info_traj["pairs"],
            correctors=selected_corrs,
            bpms=selected_bpms,
            triangular=triangular,
            actuator_mode="quadrupole_movers",
        )

        common_qcorrs = [str(q) for q in selected_corrs if str(q) in hcorrs and str(q) in vcorrs]
        if not common_qcorrs:
            raise RuntimeError("No quadrupole movers with both x and y measured response data")

        h_index = {str(name): i for i, name in enumerate(hcorrs)}
        v_index = {str(name): i for i, name in enumerate(vcorrs)}

        h_cols = [h_index[q] for q in common_qcorrs]
        v_cols = [v_index[q] for q in common_qcorrs]

        R_xx = np.asarray(Rxx[:, h_cols], dtype=float)
        R_xy = np.asarray(Rxy[:, v_cols], dtype=float)
        R_yx = np.asarray(Ryx[:, h_cols], dtype=float)
        R_yy = np.asarray(Ryy[:, v_cols], dtype=float)

        T_xx = np.zeros_like(R_xx)
        T_yy = np.zeros_like(R_yy)

        for j, q in enumerate(common_qcorrs):
            attached = f"M{q}"
            for i, bpm in enumerate(bpms_common):
                if attached in str(bpm):
                    T_xx[i, j] = -1.0
                    T_yy[i, j] = -1.0

        return {
            "qcorrs": common_qcorrs,
            "bpms": list(bpms_common),
            "R_xx": R_xx,
            "R_xy": R_xy,
            "R_yx": R_yx,
            "R_yy": R_yy,
            "T_xx": T_xx,
            "T_yy": T_yy,
        }

    def _compute_qm_response_matrix(self, pairs, qcorrs, bpms, triangular = False):
        if not pairs:
            raise RuntimeError("No DATA pairs available for quadrupole mover response matrix")

        if not hasattr(self, "sequence"):
            self.sequence = State(filename=pairs[0][0]).get_sequence()

        qcorrs = [str(q) for q in qcorrs]
        bpms = [str(bpm) for bpm in bpms]

        available_x = {str(tag)[:-2] for _fp, _fm, tag in pairs if str(tag).endswith("_x")}
        available_y = {str(tag)[:-2] for _fp, _fm, tag in pairs if str(tag).endswith("_y")}

        qcorrs_x = [q for q in qcorrs if q in available_x]
        qcorrs_y = [q for q in qcorrs if q in available_y]

        nb = len(bpms)
        Rxx_samples = {q: [] for q in qcorrs_x}
        Ryx_samples = {q: [] for q in qcorrs_x}
        Rxy_samples = {q: [] for q in qcorrs_y}
        Ryy_samples = {q: [] for q in qcorrs_y}

        Bx_rows, By_rows = [], []

        def qval(state, qname, axis):
            qdata = state.get_quadrupoles([qname])
            keys = ("xact", "xdes") if axis == "x" else ("yact", "ydes")
            for key in keys:
                if key in qdata:
                    arr = np.asarray(qdata[key], dtype=float).ravel()
                    if arr.size:
                        return float(arr[0])
            raise RuntimeError(f"No {axis} mover value for {qname}")
        for pair in pairs:
            fp, fm = pair[0], pair[1]
            tag = str(pair[2]) if len(pair) > 2 else ""

            if not (tag.endswith("_x") or tag.endswith("_y")):
                continue

            axis = tag[-1]
            qname = tag[:-2]
            Sp = State(filename=fp)
            Sm = State(filename=fm)
            # jitter_model = getattr(self, "jitter_model", None)
            # if jitter_model is not None:
            #     Sp.bpms = apply_jitter_subtraction(Sp.get_bpms(), jitter_model)
            #     Sm.bpms = apply_jitter_subtraction(Sm.get_bpms(), jitter_model)
            Op = Sp.get_orbit(bpms)
            Om = Sm.get_orbit(bpms)
            den = qval(Sp, qname, axis) - qval(Sm, qname, axis)
            if not np.isfinite(den) or den == 0:
                continue
            dx = (np.asarray(Op["x"], dtype=float) - np.asarray(Om["x"], dtype=float)) / den
            dy = (np.asarray(Op["y"], dtype=float) - np.asarray(Om["y"], dtype=float)) / den
            if axis == "x" and qname in Rxx_samples:
                Rxx_samples[qname].append(dx)
                Ryx_samples[qname].append(dy)
            if axis == "y" and qname in Ryy_samples:
                Rxy_samples[qname].append(dx)
                Ryy_samples[qname].append(dy)

            Bx_rows += [np.asarray(Op["x"], dtype=float), np.asarray(Om["x"], dtype=float)]
            By_rows += [np.asarray(Op["y"], dtype=float), np.asarray(Om["y"], dtype=float)]

        def mean_or_nan(samples):
            if not samples:
                return np.full(nb, np.nan)
            return np.nanmean(np.vstack(samples), axis=0)

        Rxx = np.column_stack([mean_or_nan(Rxx_samples[q]) for q in qcorrs_x])
        Ryx = np.column_stack([mean_or_nan(Ryx_samples[q]) for q in qcorrs_x])

        Rxy = np.column_stack([mean_or_nan(Rxy_samples[q]) for q in qcorrs_y])
        Ryy = np.column_stack([mean_or_nan(Ryy_samples[q]) for q in qcorrs_y])

        if not Bx_rows or not By_rows:
            raise RuntimeError("No valid quadrupole mover DATA pairs found")

        Bx = np.nanmean(np.vstack(Bx_rows), axis=0).reshape(-1, 1)
        By = np.nanmean(np.vstack(By_rows), axis=0).reshape(-1, 1)

        if triangular:
            sequence_index = {str(name): i for i, name in enumerate(self.sequence)}

            def sequence_pos(name):
                return sequence_index.get(str(name), np.inf)

            for j, q in enumerate(qcorrs_x):
                qpos = sequence_pos(q)
                for i, bpm in enumerate(bpms):
                    if sequence_pos(bpm) < qpos or str(bpm) == f"M{q}" or f"M{q}" in str(bpm):
                        Rxx[i, j] = 0.0
                        Ryx[i, j] = 0.0

            for j, q in enumerate(qcorrs_y):
                qpos = sequence_pos(q)
                for i, bpm in enumerate(bpms):
                    if sequence_pos(bpm) < qpos or str(bpm) == f"M{q}" or f"M{q}" in str(bpm):
                        Rxy[i, j] = 0.0
                        Ryy[i, j] = 0.0

        return Rxx, Ryy, Rxy, Ryx, Bx, By, qcorrs_x, qcorrs_y, bpms
