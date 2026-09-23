import RF_Track as rft
import sys
import numpy as np
import time, math
from Interfaces.AbstractMachineInterface import AbstractMachineInterface
# must run on SLAC controls network
sys.path.append('/usr/local/facet/tools/python/F2_live_model/')
sys.path.append('/usr/local/facet/tools/python/F2_pytools/')
# sys.path.append('/home/fphysics/zack/workspace/F2_pytools/')
try:
    from bmad import BmadLiveModel
    from epics import get_pv
    from F2_pytools.controls_jurisdiction import is_SLC
    from F2_pytools.f2bsaBuffer import make_bpm_buffer, get_bpmdata2
    from F2_pytools.mags import set_magnets
except ImportError:
    bmad, epics, F2_pytools  = None, None, None
from traceback import print_exception
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "FACET2"))
import FACET2
'''
TO IMPLEMENT!!!

def get_screens(self):
    pass
def set_quadrupoles(self):
    pass

'''


# missing/deffered elements in the model to ignore
BPM_BLACKLIST = [
    'BPM19851','BPM19871'
    ]
CORRECTOR_BLACKLIST = [
    'YC57145','YC57146',
    'XCB2LE','XCB3RE','XCB3LE','XCB2RE','XC1EX',
    ]

def devname_swap_micro_primary(device):
    ds = device.split(':')
    return f'{ds[1]}:{ds[0]}:{ds[2]}'

class InterfaceFACET2_Linac(AbstractMachineInterface):
    def get_name(self):
        return 'FACET2_Linac'

    def __init__(self, nsamples=10, livemodel=False, tao_initfile=None):
        self.log = print
        self.lattice = FACET2.load_FACET()
        self.nsamples = nsamples
        if livemodel:
            self.f2m = BmadLiveModel(instanced=True, init_filename=tao_initfile)
            # Brho = self.f2m.live.Brho
        else:
            self.f2m = BmadLiveModel(design_only=True)
            # Brho = self.f2m.design.Brho
        _ix = self.f2m.ix
        self.sequence = list(self.f2m.elements[np.sort(np.append(_ix['BPMS'], _ix['COR']))])
        self.bpms = list(self.f2m.elements[_ix['BPMS']])
        self.xcorrs = list(self.f2m.elements[_ix['XCOR']])
        self.ycorrs = list(self.f2m.elements[_ix['YCOR']])
        self.corrs = list(self.f2m.elements[_ix['COR']])
        for bpmname in BPM_BLACKLIST:
            self.sequence.remove(bpmname)
            self.bpms.remove(bpmname)
        for corname in CORRECTOR_BLACKLIST:
            self.sequence.remove(corname)
            self.corrs.remove(corname)
            if corname.startswith('X'): self.xcorrs.remove(corname)
            if corname.startswith('Y'): self.ycorrs.remove(corname)
        self.bpmdevs =  [self.f2m.device_names[_ix[bpm]] for bpm in self.bpms]
        self.bpms_s =   [self.f2m.S[_ix[bpm]]  for bpm  in self.bpms]
        self.xcorrs_s = [self.f2m.S[_ix[xcor]] for xcor in self.xcorrs]
        self.ycorrs_s = [self.f2m.S[_ix[ycor]] for ycor in self.ycorrs]
        self.corrs_s =  [self.f2m.S[_ix[cor]]  for cor  in self.corrs]
        # self.Brho =     [Brho[ix[elem]] for elem in self.sequence]
        self.PVs = {
            'Q_setpoint':     get_pv('SIOC:SYS1:ML03:AO518'),
            'Q_readback':     get_pv('TORO:IN10:431:TMIT_PC'),
            'UVWP_angle':     get_pv('WPLT:LT10:150:WP_ANGLE'),
            'dl10_en':        get_pv('BEND:IN10:751:BDES'),
            'bc11_en':        get_pv('BEND:LI11:314:BDES'),
            'bc14_en':        get_pv('BEND:LI14:720:BDES'),
            'bc20_en':        get_pv('LI20:LGPS:1990:BDES'),
            'dl10e_setpoint': get_pv('PHYS:SYS1:1:F2LFB_DL10E_VERN'),
            'bc11e_setpoint': get_pv('PHYS:SYS1:1:F2LFB_BC11E_VERN'),
            'bc14e_setpoint': get_pv('PHYS:SYS1:1:F2LFB_BC14E_VERN'),
            'bc20e_setpoint': get_pv('PHYS:SYS1:1:F2LFB_BC20E_VERN'),
            }
        # initial bunch charge setpoint for reset_intensity
        self.init_charge_setpoint = self.PVs['Q_setpoint'].get()
        self.UVWP_init = self.PVs['UVWP_angle'].get()
        self.bba_bc11_energy_offset_mev = -3.0
        self.bba_downstream_energy_offset_mev = -40.0
        self.bba_uvwp_offset_deg = -2.5
        # initialize bpm data buffer
        self.bpm_buffer = make_bpm_buffer(self.f2m, self.bpms, Npts=self.nsamples)
        print('InterfaceFACET2_Linac is ready')
        self.machine_name = "FACET"

    def log_messages(self,console):
        self.log=console or print

    @staticmethod
    def _safe_float(value, default=np.nan):
        try:
            if value is None:
                return float(default)
            arr = np.asarray(value)
            if arr.size == 0:
                return float(default)
            return float(arr.flat[0])
        except Exception:
            return float(default)

    def _give_elements_to_show_beamline(self, quad_selected):
        start_quad_element_name = quad_selected
        return start_quad_element_name

    def _get_elements_positions_show_beamline(self, names=None):
        if isinstance(names, str):
            names = [names]
        all_names = []
        all_s = []
        all_l = []
        s_pos = 0.0
        for element in self.lattice['*']:
            element_name = element.get_name()
            try:
                element_length = float(element.get_length())
            except Exception:
                element_length = 0.0
            if names is None or element_name in names:
                all_names.append(element_name)
                all_s.append(s_pos)
                all_l.append(element_length)
            s_pos += element_length
        return {
            "names": all_names,
            "S": np.array(all_s, dtype=float),
        }

    def _wait_for_magnet_readback(self, devname, target, tolerance=1e-4, timeout=1.0, poll_interval=0.05):
        bact_pv = get_pv(f'{devname}:BACT')
        t0 = time.perf_counter()
        last_value = np.nan

        while time.perf_counter() - t0 < timeout:
            try:
                last_value = self._safe_float(bact_pv.get(), default=np.nan)
            except Exception:
                last_value = np.nan

            if np.isfinite(last_value) and abs(last_value - float(target)) <= tolerance:
                return True

            time.sleep(poll_interval)

        self.log(
            f'Warning: {devname}:BACT did not reach target {float(target):.6g} '
            f'within {timeout:.2f}s. Last readback = {last_value:.6g}'
        )
        return False

    def _wait_for_pv_readback(self, pv, target, *, description, tolerance=1e-4, timeout=10.0):
        return self._wait_for_readback(lambda: self._safe_float(pv.get(), default=np.nan), target, description=description, tolerance=tolerance, timeout=timeout)

    def _meascharge(self):
        qraw = []
        for j in range(self.nsamples):
            qraw.append(self.PVs['Q_readback'].get())
            time.sleep(0.1)
        return np.nanmean(qraw)

    def get_beam_factors(self):
        # TO BE REPLACED WITH A PV OF REAL BEAM ENERGY
        Pref = np.sqrt(125.0**2 - 0.51099895**2)
        gamma_rel = np.sqrt((Pref / 0.51099895) ** 2 + 1.0)
        beta_rel = np.sqrt(1.0 - 1.0 / gamma_rel ** 2)
        return gamma_rel, beta_rel

    def _pv_readback(self, pv):
        return lambda: self._safe_float(pv.get(), default=np.nan)

    def change_energy(self):
        if np.isclose(float(self.bba_bc11_energy_offset_mev), 0.0) or np.isclose(float(self.bba_downstream_energy_offset_mev), 0.0):
            raise RuntimeError(
                "bba_bc11_energy_offset_mev / bba_downstream_energy_offset_mev is 0; reset_energy() "
                "sets these PVs to 0.0 too, so a DFS would be meaningless."
            )
        print('Lowering beam energy starting from BC11')
        targets = {
            'bc11e_setpoint': self.bba_bc11_energy_offset_mev,
            'bc14e_setpoint': self.bba_downstream_energy_offset_mev,
            'bc20e_setpoint': self.bba_downstream_energy_offset_mev,
        }
        for name, target in targets.items():
            self._set_and_verify(lambda name=name, target=target: self.PVs[name].put(target), self._pv_readback(self.PVs[name]), target, description=name)
        # get_pv(f'PHYS:SYS1:1:F2LFB_BC11BL_TARGET').put(4400)
        # get_pv(f'PHYS:SYS1:1:F2LFB_BC14BL_TARGET').put(5000)
        return float(self.bba_bc11_energy_offset_mev) / 335.0

    def reset_energy(self):
        print('Restoring beam energy at BC11, re-enabling feedbacks')
        for name in ('bc11e_setpoint', 'bc14e_setpoint', 'bc20e_setpoint'):
            self._set_and_verify(lambda name=name: self.PVs[name].put(0.0), self._pv_readback(self.PVs[name]), 0.0, description=name)
        # get_pv(f'PHYS:SYS1:1:F2LFB_BC11BL_TARGET').put(5000)
        # get_pv(f'PHYS:SYS1:1:F2LFB_BC14BL_TARGET').put(8000)

    def change_intensity(self):
        if np.isclose(float(self.bba_uvwp_offset_deg), 0.0):
            raise RuntimeError("bba_uvwp_offset_deg is 0; a WFS would be meaningless.")
        self.UVWP_init = self.PVs['UVWP_angle'].get()
        self.Q_init = self._meascharge()
        uvwp_target = self.UVWP_init + self.bba_uvwp_offset_deg
        self._set_and_verify(lambda: self.PVs['UVWP_angle'].put(uvwp_target), self._pv_readback(self.PVs['UVWP_angle']), uvwp_target, description='UVWP_angle')
        self.Q_new = self._meascharge()
        self._set_and_verify(lambda: self.PVs['Q_setpoint'].put(self.Q_new), self._pv_readback(self.PVs['Q_setpoint']), self.Q_new, description='Q_setpoint')
        print(f'Charge changed: {self.Q_init:.1f}  {self.Q_new:.1f} pC')
        return self

    def reset_intensity(self):
        print(f'restoring charge setpoint to {self.init_charge_setpoint} pC')
        self._set_and_verify(lambda: self.PVs['UVWP_angle'].put(self.UVWP_init), self._pv_readback(self.PVs['UVWP_angle']), self.UVWP_init, description='UVWP_angle')
        self._set_and_verify(lambda: self.PVs['Q_setpoint'].put(self.init_charge_setpoint), self._pv_readback(self.PVs['Q_setpoint']), self.init_charge_setpoint, description='Q_setpoint')
        return self

    def get_beam_settings(self):
        energy_pvs = {
            "dl10_feedback_vernier": "dl10e_setpoint",
            "bc11_feedback_vernier": "bc11e_setpoint",
            "bc14_feedback_vernier": "bc14e_setpoint",
            "bc20_feedback_vernier": "bc20e_setpoint",
        }
        intensity_pvs = {
            "uvwp_angle": "UVWP_angle",
            "charge_setpoint": "Q_setpoint",
        }
        return {
            "energy": {name: self._safe_float(self.PVs[pv].get()) for name, pv in energy_pvs.items()},
            "intensity": {name: self._safe_float(self.PVs[pv].get()) for name, pv in intensity_pvs.items()},
        }

    def restore_beam_settings(self, settings):
        settings = settings or {}
        energy_pvs = {
            "dl10_feedback_vernier": "dl10e_setpoint",
            "bc11_feedback_vernier": "bc11e_setpoint",
            "bc14_feedback_vernier": "bc14e_setpoint",
            "bc20_feedback_vernier": "bc20e_setpoint",
        }
        intensity_pvs = {
            "uvwp_angle": "UVWP_angle",
            "charge_setpoint": "Q_setpoint",
        }
        for values, pv_names in ((settings.get("energy", {}), energy_pvs),
                                 (settings.get("intensity", {}), intensity_pvs)):
            for name, pv_name in pv_names.items():
                value = self._safe_float(values.get(name))
                if np.isfinite(value):
                    self.PVs[pv_name].put(value)
                    self._wait_for_pv_readback(self.PVs[pv_name], value, description=pv_name)
        return True

    def get_icts(self, names=None):
        return {"names": [], "charge": np.array([])}

    def get_sequence(self):
        return self.sequence

    def get_hcorrectors_names(self): return self.xcorrs

    def get_target_dispersion(self, names=None):
        if names is None:
            names = self.get_bpms()["names"]
        return [0. for _ in names], [0. for _ in names]

    def get_vcorrectors_names(self): return self.ycorrs

    def get_elements_position(self, names): return [self.f2m.S[self.f2m.ix[name]] for name in names]

    def get_correctors(self, names=None):
        print("Reading correctors' strengths...")
        corr_names, bdes, bact = [], [], []
        for corname in self.corrs:
            devname = self.f2m.device_names[self.f2m.ix[corname]]
            if devname == '':
                print(f'skipping {corname}')
                continue
            corr_names.append(corname)
            bdes.append(get_pv(f'{devname}:BDES').get())
            bact.append(get_pv(f'{devname}:BACT').get())

        correctors = {
            "names": corr_names,
            "bdes": np.array(bdes),
            "bact": np.array(bact),
        }

        if isinstance(names, str):
            names = [names]
        if names is not None:
            idx = [i for i, s in enumerate(correctors["names"]) if s in names]
            correctors = {
                "names": [correctors["names"][i] for i in idx],
                "bdes": np.asarray(correctors["bdes"])[idx],
                "bact": np.asarray(correctors["bact"])[idx],
            }

        return correctors

    def get_bpms(self, names=None):
        itry = 0
        while itry < 10:
            try:
                bpmdata = get_bpmdata2(self.bpm_buffer)
                if isinstance(names, str):
                    names = [names]
                if names is not None:
                    idx = [i for i, s in enumerate(bpmdata["names"]) if s in names]
                    bpmdata = {
                        "names": [bpmdata["names"][i] for i in idx],
                        "x": np.asarray(bpmdata["x"])[:, idx],
                        "y": np.asarray(bpmdata["y"])[:, idx],
                        "tmit": np.asarray(bpmdata["tmit"])[:, idx],
                    }
                return bpmdata
            except Exception as e:
                print('BPM data acquisition failed! trying again ...')
                print_exception(e)
                itry += 1
        else:
            print(f'unable to get BPMs after {itry} tries')
            bpmdata = self._dummy_bpmdata()
            if isinstance(names, str):
                names = [names]
            if names is not None:
                idx = [i for i, s in enumerate(bpmdata["names"]) if s in names]
                bpmdata = {
                    "names": [bpmdata["names"][i] for i in idx],
                    "x": np.asarray(bpmdata["x"])[:, idx],
                    "y": np.asarray(bpmdata["y"])[:, idx],
                    "tmit": np.asarray(bpmdata["tmit"])[:, idx],
                }
            return bpmdata

    def _dummy_bpmdata(self):
        """ mock BPM data, all nans """
        N, M = len(self.bpm_buffer.devicelist), self.bpm_buffer.EPICS_Npts
        xraw = np.empty((N, M))
        yraw = np.empty((N, M))
        qraw = np.empty((N, M))
        xraw.fill(np.nan)
        yraw.fill(np.nan)
        qraw.fill(np.nan)
        return {
            'names': self.bpm_buffer.elemlist,
            'x': xraw.transpose(),
            'y': yraw.transpose(),
            'tmit': qraw.transpose(),
        }

    def set_correctors(self, names, corr_vals):
        """ write BDES to correctors """
        if isinstance(names, str):
            names = [names]
        if not isinstance(corr_vals, (list, tuple, np.ndarray)):
            corr_vals = [corr_vals]
        if len(names) != len(corr_vals):
            raise ValueError('len(names) != len(corr_vals) in set_correctors(names, corr_vals)')
        devnames = [self.f2m.device_names[self.f2m.ix[name]] for name in names]
        set_magnets(devnames, corr_vals, perturb=True)
        for devname, target in zip(devnames, corr_vals):
            self._wait_for_magnet_readback(devname, target)

    def vary_correctors(self, names, corr_vals):
        """ write deltas to correctors """
        if isinstance(names, str):
            names = [names]
        if not isinstance(corr_vals, (list, tuple, np.ndarray)):
            corr_vals = [corr_vals]
        if len(names) != len(corr_vals):
            raise ValueError('len(names) != len(corr_vals) in vary_correctors(names, corr_vals)')
        devnames = [self.f2m.device_names[self.f2m.ix[name]] for name in names]
        bdes_init = [get_pv(f'{devname}:BDES').get() for devname in devnames]
        updated_corr_vals = [bdes+delta for bdes,delta in zip(bdes_init, corr_vals)]
        set_magnets(devnames, updated_corr_vals, perturb=True)
        for devname, target in zip(devnames, updated_corr_vals):
            self._wait_for_magnet_readback(devname, target)

    def get_elements_indices(self, names):
        if isinstance(names, str):
            names = [names]
        name_to_index = {string: index for index, string in enumerate(self.sequence)}
        return [name_to_index.get(name, np.nan) for name in names]
