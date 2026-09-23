import os, ctypes

class MagKiWrapper:
    MODE_K_TO_I = 1
    MODE_I_TO_K = 2

    def __init__(self, library_path):
        self.library_path = os.path.abspath(str(library_path))
        self.lib = ctypes.CDLL(self.library_path)
        self._func = self.lib.mag_ki_main
        self._func.argtypes = [
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_float,
            ctypes.POINTER(ctypes.c_float),
            ctypes.POINTER(ctypes.c_float),
            ctypes.POINTER(ctypes.c_float),
            ctypes.POINTER(ctypes.c_float),
        ]
        self._func.restype = ctypes.c_int

    def _call(self, mode, name, energy_GeV, k_main=0.0, current_main=0.0):
        kvalue = (ctypes.c_float * 2)(float(k_main), 0.0)
        current = (ctypes.c_float * 2)(float(current_main), 0.0)
        field = (ctypes.c_float * 2)(0.0, 0.0)
        efflen = ctypes.c_float(0.0)
        status = int(
            self._func(int(mode), str(name).encode("ascii"), ctypes.c_float(float(energy_GeV)), kvalue, current, ctypes.byref(efflen), field))
        if status != 1:
            raise RuntimeError(f"mag_ki_main failed for {name}, mode={mode}, status={status}")
        return {
            "k": float(kvalue[0]),
            "current": float(current[0]),
            "efflen": float(efflen.value),
            "field": float(field[0]),
        }

    def current_to_k1l(self, name, current_A, energy_GeV):
        return self._call(self.MODE_I_TO_K, name, energy_GeV, current_main=current_A)["k"]

    def k1l_to_current(self, name, k1, energy_GeV):
        return self._call(self.MODE_K_TO_I, name, energy_GeV, k_main=k1)["current"]

def load_mag_ki():
    candidates = []
    env_path = os.environ.get("ATF2_MAG_KI_LIB", "")
    if env_path:
        candidates.append(env_path)
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    libmagnet_dir = os.path.join(repo_root, "Machine specifics, user implementations", "ATF2", "libmagnet")
    candidates.extend([
        os.path.join(libmagnet_dir, "libmagnet.dylib"),
        os.path.join(libmagnet_dir, "libmagnet.so"),
    ])

    for library_path in candidates:
        if not library_path or not os.path.exists(library_path):
            continue
        try:
            mag_ki = MagKiWrapper(library_path)
            print(f"Loaded ATF2 mag_ki library: {library_path}")
            return mag_ki
        except Exception as exc:
            print(f"ATF2 mag_ki library '{library_path}': {exc} not loaded")

    print("ATF2 mag_ki library not loaded.")
    return None
