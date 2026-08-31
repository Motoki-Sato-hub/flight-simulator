import os, sys, pyjapc, time
import pyda, pyda_japc
from pathlib import Path
import numpy as np
project_root_path = Path.cwd().resolve()
while (not (project_root_path / "Interfaces").exists() and project_root_path.parent != project_root_path):
    project_root_path = project_root_path.parent
sys.path.insert(0, str(project_root_path))
os.chdir(project_root_path)
from Interfaces.CLEAR.InterfaceCLEAR import CLEAR_real_machine
import matplotlib.pyplot as plt

I = CLEAR_real_machine()

client = pyda.SimpleClient(provider=pyda_japc.JapcProvider())
context_acquisition = "SCT.USER.SETUP"
context_empty = ""

"""Test of one BPM and method get_bpms"""

print("Testing BPM CA.BPM0530H-SA, reading its values directly from japc/pyda... ")
bpm = "CA.BPM0530"
hsamples = client.get(f"{bpm}H-SA/SamplesFromTrigger", context = context_acquisition).data
print(f"{bpm}H-SA/SamplesFromTrigger = {hsamples}")

vsamples = client.get(f"{bpm}V-SA/SamplesFromTrigger", context = context_acquisition).data
print(f"{bpm}V-SA/SamplesFromTrigger = {vsamples}")

Ssamples = client.get(f"{bpm}S-SA/SamplesFromTrigger", context = context_acquisition).data
print(f"{bpm}S-SA/SamplesFromTrigger = {Ssamples}")

h = np.asarray(hsamples["samples"], dtype=float).ravel()
v = np.asarray(vsamples["samples"], dtype=float).ravel()
s = np.asarray(Ssamples["samples"], dtype=float).ravel()
Hpos = np.sum(h)/np.sum(s)
Vpos = np.sum(v)/np.sum(s)

print("Hpos:", Hpos)
print("Vpos:", Vpos)
# plt.figure()
# plt.plot(h, label="H")
# plt.plot(v, label="V")
# plt.plot(s, label="S")
# plt.legend()
# plt.show()

print("================================================================================")
print("Testing BPM CA.BPM0530H-SA, reading its values directly from interface, using get_bpms method... ")
result = I.get_bpms(names=["CA.BPM0530", "CA.BPM0595"])
print(f"names from get_bpms: {result['names']}")
print(f"x from get_bpms: {result['x']}")
print(f"y from get_bpms: {result['y']}")
print(f"tmit from get_bpms: {result['tmit']}")
print("================================================================================")

