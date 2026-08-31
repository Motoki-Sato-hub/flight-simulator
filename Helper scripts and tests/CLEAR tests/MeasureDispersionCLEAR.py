import os, sys, time
import pyda, pyda_japc, pyjapc
import os, sys
from pathlib import Path
import numpy as np
project_root_path = Path.cwd().resolve()
while (not (project_root_path / "Interfaces").exists() and project_root_path.parent != project_root_path):
    project_root_path = project_root_path.parent
sys.path.insert(0, str(project_root_path))
os.chdir(project_root_path)
from Interfaces.CLEAR.InterfaceCLEAR import CLEAR_real_machine

#from Interfaces.CLEAR.InterfaceCLEAR_RFTrack import InterfaceCLEAR_RFTrack

import matplotlib.pyplot as plt

#I = InterfaceCLEAR_RFTrack()
I = CLEAR_real_machine(nsamples=20)
bpms = I.bpms
bpm0 = I.get_bpms(bpms)
output_dir = Path.home() / "CERN-Flight_Simulator-Data" / "CLEAR_dispersion"
output_dir.mkdir(parents=True, exist_ok=True)

for i in range(5):
    bpm0 = I.get_bpms(bpms)
    I.change_energy()
    bpm1 = I.get_bpms(bpms)
    I.reset_energy()

    x0_samples = np.asarray(bpm0["x"], dtype=float)
    y0_samples = np.asarray(bpm0["y"], dtype=float)
    x1_samples = np.asarray(bpm1["x"], dtype=float)
    y1_samples = np.asarray(bpm1["y"], dtype=float)

    x0 = np.nanmean(x0_samples, axis=0)
    y0 = np.nanmean(y0_samples, axis=0)
    x1 = np.nanmean(x1_samples, axis=0)
    y1 = np.nanmean(y1_samples, axis=0)

    dx = x1 - x0
    dy = y1 - y0

    output_file = output_dir / f"disp_measurement_{i:03d}.txt"

    with open(output_file, "w") as file:
        file.write("# BPM\tx0\ty0\tx1\ty1\tdx\tdy\n")
        for name, a, b, c, d, e, f in zip(bpms, x0, y0, x1, y1, dx, dy):
            file.write(
                f"{name}\t{a:.8g}\t{b:.8g}\t{c:.8g}\t{d:.8g}\t"
                f"{e:.8g}\t{f:.8g}\n"
            )