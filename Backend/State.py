from datetime import datetime
from zoneinfo import ZoneInfo
import numpy as np
import pickle
import gzip

def reject_large_outliers(values, factor=10.0):
    arr = np.asarray(values, dtype=float).copy()
    if arr.ndim != 2 or arr.size == 0:
        return arr
    med_abs = np.nanmedian(np.abs(arr), axis=0)
    threshold = factor * med_abs
    for j, thr in enumerate(threshold):
        if not np.isfinite(thr) or thr <= 0:
            continue
        mask = np.isfinite(arr[:, j]) & (np.abs(arr[:, j]) > thr)
        arr[mask, j] = np.nan

    return arr

class State:
    def __init__(self, sextupoles = None, correctors=None,bpms=None, icts=None,sequence=None,hcorrectors_names=None,vcorrectors_names=None,screens=None,quadrupoles=None,beam_settings=None,timestamp=None,filename=None,interface_id=None):
        if filename is not None:
            self.load(filename)
            return
        self.correctors = correctors if correctors is not None else {"names": [], "bdes": np.array([]), "bact": np.array([])}
        self.bpms = bpms if bpms is not None else {"names": [], "x": np.empty((0, 0)), "y": np.empty((0, 0)), "tmit": np.empty((0, 0))}
        self.icts = icts if icts is not None else {"names": [], "charge": np.array([])}
        self.sequence = sequence if sequence is not None else []
        self.hcorrectors_names = hcorrectors_names if hcorrectors_names is not None else []
        self.vcorrectors_names = vcorrectors_names if vcorrectors_names is not None else []
        self.screens = screens if screens is not None else {"names": [], "hpixel": np.array([]), "vpixel": np.array([]), "x":np.array([]),"y":np.array([]), "sigx":np.array([]), "sigy":np.array([]),"sum":np.array([]),"hedges":[],"vedges":[],"images":[],"S":np.array([])}
        self.quadrupoles = quadrupoles if quadrupoles is not None else {"names": [], "bdes": np.array([]), "bact": np.array([]), "xdes": np.array([]), "ydes": np.array([]), "rolldes": np.array([])}
        self.sextupoles = sextupoles if sextupoles is not None else {"names": [], "bdes": np.array([]), "bact": np.array([])}
        self.beam_settings = beam_settings if beam_settings is not None else {}
        self.interface_id = interface_id
        self.timestamp = timestamp if timestamp is not None else datetime.now()

    def get_sequence(self):
        return self.sequence #from rf track

    def get_correctors(self, names=None):
        if isinstance(names, str):
            names = [names]
        if names is not None:
            corr_indexes = np.array([index for index, string in enumerate(self.correctors['names']) if string in names], dtype=int)
            correctors = {
                "names": np.array(self.correctors['names'])[corr_indexes],
                "bdes": np.array(self.correctors['bdes'])[corr_indexes],
                "bact": np.array(self.correctors['bact'])[corr_indexes]
            }
        else:
            correctors = self.correctors
        return correctors

    def get_bpms(self, names=None):
        if isinstance(names, str):
            names = [names]
        if names is not None:
            bpm_indexes = np.array([index for index, string in enumerate(self.bpms['names']) if string in names], dtype=int)

            bpms = {
                "names": np.array(self.bpms['names'])[bpm_indexes],
                "x": np.array(self.bpms['x'])[:,bpm_indexes],
                "y": np.array(self.bpms['y'])[:,bpm_indexes],
                "tmit": np.array(self.bpms['tmit'])[:,bpm_indexes],
            }
        else:
            bpms = self.bpms
        return bpms

    def get_icts(self, names=None):
        if isinstance(names, str):
            names = [names]
        icts = self.icts
        if names is not None:
            ict_indexes = np.array([index for index, string in enumerate(icts['names']) if string in names], dtype=int)
            icts = {
                "names": np.array(self.icts['names'])[ict_indexes],
                "charge": np.array(self.icts['charge'])[ict_indexes]
            }
        return icts

    def get_quadrupoles(self, names=None):
        if isinstance(names, str):
            names = [names]
        quadrupoles = self.quadrupoles
        if names is not None:
            quadrupole_indexes = np.array([index for index, string in enumerate(quadrupoles.get('names', [])) if string in names], dtype=int)
            subset = {"names": np.array(quadrupoles.get('names', []))[quadrupole_indexes]}
            for key, value in quadrupoles.items():
                if key == "names":
                    continue
                try:
                    subset[key] = np.array(value)[quadrupole_indexes]
                except Exception:
                    subset[key] = np.array([])
            quadrupoles = subset
        return quadrupoles

    def get_sextupoles(self, names=None):
        if isinstance(names, str):
            names = [names]
        sextupoles=self.sextupoles
        if names is not None:
            sextupole_indexes = np.array([index for index, string in enumerate(sextupoles['names']) if string in names], dtype=int)
            sextupoles = {
                "names": np.array(self.sextupoles['names'])[sextupole_indexes],
                "bdes": np.array(self.sextupoles['bdes'])[sextupole_indexes],
                "bact": np.array(self.sextupoles['bact'])[sextupole_indexes],
            }
        return sextupoles

    def get_orbit(self, names=None):
        bpms = self.get_bpms(names)
        # x = np.mean(bpms['x'],axis=0) # mm
        # y = np.mean(bpms['y'],axis=0) # mm
        # stdx = np.std(bpms['x'],axis=0) # mm #standard deviation
        # stdy = np.std(bpms['y'],axis=0) # mm
        # tmit = np.mean(bpms['tmit'],axis=0)
        # nshots=int(np.shape(bpms['x'])[0])
        # faulty = (x == 0.0) & (y == 0.0)
        x_all = reject_large_outliers(bpms['x'])
        y_all = reject_large_outliers(bpms['y'])
        tmit_all = np.asarray(bpms['tmit'], dtype=float)
        x = np.nanmean(x_all, axis=0)  # mm
        y = np.nanmean(y_all, axis=0)  # mm
        stdx = np.nanstd(x_all, axis=0)  # mm
        stdy = np.nanstd(y_all, axis=0)  # mm
        tmit = np.nanmean(tmit_all, axis=0)
        nshots = int(np.shape(x_all)[0])
        faulty = np.isnan(x) | np.isnan(y) | ((x == 0.0) & (y == 0.0))
        x[faulty] = np.nan
        y[faulty] = np.nan
        orbit = {
            "names": bpms["names"],
            "x": x,
            "y": y,
            "stdx": stdx,
            "stdy": stdy,
            "tmit": tmit,
            "faulty": faulty,
            "nbpms": len(bpms["names"]),
            "nshots": nshots,
        }
        return orbit

    def get_beam_settings(self):
        return self.beam_settings

    def get_interface_id(self):
        return self.interface_id

    def get_timestamp(self):
        return self.timestamp

    def get_screens(self,names=None):
        if isinstance(names, str):
            names = [names]
        if names is not None:
            screen_indexes = np.array([index for index, string in enumerate(self.screens['names']) if string in names], dtype=int)
            screens = {"names": np.array(self.screens['names'])[screen_indexes],
                       "hpixel": np.array(self.screens['hpixel'])[screen_indexes],
                       "vpixel": np.array(self.screens['vpixel'])[screen_indexes],
                       "x": np.array(self.screens['x'])[screen_indexes],
                       "y": np.array(self.screens['y'])[screen_indexes],
                       "sigx": np.array(self.screens['sigx'])[screen_indexes],
                       "sigy": np.array(self.screens['sigy'])[screen_indexes],
                       "sum": np.array(self.screens['sum'])[screen_indexes],
                        "S": np.array(self.screens['S'])[screen_indexes],
                       "hedges": [self.screens['hedges'][i] for i in screen_indexes],
                       "vedges": [self.screens['vedges'][i] for i in screen_indexes],
                       "images": [self.screens['images'][i] for i in screen_indexes],
                       }
        else:
            screens = self.screens
        return screens

    def load(self, filename):
        from glob import glob
        f = glob(f'{filename}*')
        if len(f)==0:
            raise FileNotFoundError(f"Couldn't find state file matching {filename}")
        try:
            try:
                with gzip.open(f[0], "rb") as pickle_file:
                    data = pickle.load(pickle_file)
            except OSError:
                with open(f[0], "rb") as pickle_file:
                    data = pickle.load(pickle_file)
            self.sequence = data['sequence']
            self.correctors = data['correctors']
            self.bpms = data['bpms']
            self.icts = data['icts']
            self.screens = data.get('screens',
                                    {"names": [], "hpixel": np.array([]), "vpixel": np.array([]), "x": np.array([]),
                                     "y": np.array([]), "sigx": np.array([]), "sigy": np.array([]), "sum": np.array([]),
                                     "hedges": [], "vedges": [], "images": [], "S": np.array([])})
            self.quadrupoles = data.get('quadrupoles', {"names": [], "bdes": np.array([]), "bact": np.array([]), "xdes": np.array([]), "ydes": np.array([]), "rolldes": np.array([])})
            self.sextupoles = data.get('sextupoles', {"names": [], "bdes": np.array([]), "bact": np.array([])})
            self.beam_settings = data.get('beam_settings', {})
            self.interface_id = data.get('interface_id')
            self.hcorrectors_names = data['hcorrectors_names']
            self.vcorrectors_names = data['vcorrectors_names']
            timestamp_iso = data.get("timestamp_iso")
            if timestamp_iso:
                self.timestamp = datetime.fromisoformat(timestamp_iso)
                timestamp_timezone = data.get("timestamp_timezone")
                if timestamp_timezone:
                    self.timestamp = self.timestamp.astimezone(ZoneInfo(timestamp_timezone))
            else:
                self.timestamp = datetime.strptime(data['timestamp'], "%Y/%m/%d, %H:%M:%S")
        except Exception:
            raise Exception(f"Could not load {filename}")

    def save(self, basename=None, filename=None):
        if basename is not None:
            time_str = self.timestamp.strftime("%Y%m%d_%H%M%S")
            filename = f"{basename}_{time_str}.pkl"
        correctors = {
            'names': self.correctors['names'],
            'bdes': self.correctors['bdes'], #setpoint for a corrector
            'bact': self.correctors['bact'] #readback for a corr
        }
        bpms = {
            'names': self.bpms['names'],
            'x': self.bpms['x'],
            'y': self.bpms['y'],
            'tmit': self.bpms['tmit'] #it's a local intensity at each bpm
        }
        icts = { #tells us about intensity
            'names': self.icts['names'],
            'charge': self.icts['charge'] #intensity
        }

        quadrupoles = dict(self.quadrupoles)

        sextupoles = {
            'names': self.sextupoles['names'],
            "bact": self.sextupoles['bact'],
            "bdes": self.sextupoles['bdes'],
        }

        screens={
            'names': self.screens['names'],
            'hpixel': self.screens['hpixel'],
            'vpixel': self.screens['vpixel'],
            'x': self.screens['x'],
            'y': self.screens['y'],
            'sigx': self.screens['sigx'],
            'sigy': self.screens['sigy'],
            'sum': self.screens['sum'],
            'hedges': self.screens['hedges'],
            'vedges': self.screens['vedges'],
            'images': self.screens['images'],
            'S': self.screens['S'],
        }

        state = {
            "sequence": self.sequence,
            "correctors": correctors,
            "bpms": bpms,
            "icts": icts,
            "screens": screens,
            "quadrupoles": quadrupoles,
            "sextupoles": sextupoles,
            "beam_settings": self.beam_settings,
            "interface_id": self.interface_id,
            "hcorrectors_names": self.hcorrectors_names,
            "vcorrectors_names": self.vcorrectors_names,
            "timestamp": self.timestamp.strftime("%Y/%m/%d, %H:%M:%S"),
            "timestamp_iso": self.timestamp.isoformat(timespec="seconds"),
            "timestamp_timezone": getattr(self.timestamp.tzinfo, "key", None),
        }
        if filename is None and basename is None:
            raise ValueError("Either filename or basename is required")
        with gzip.open(filename, "wb", compresslevel=6) as file:
            pickle.dump(state, file)
        return filename
