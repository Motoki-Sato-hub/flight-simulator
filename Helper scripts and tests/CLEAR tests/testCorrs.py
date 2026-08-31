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

"""Test of one corrector and method get_correctors"""

print("Testing corrector CA.DVG0130, reading its values directly from japc/pyda... ")

corrector = "CA.DVG0130"

setting = client.get(f"{corrector}/SettingPPM", context = context_empty).data["current"]
print(f"{corrector}/SettingPPM#current = {setting}")
client.set(f"{corrector}/SettingPPM", data={"current": 0.0})
time.sleep(5)
print(f"{corrector}/SettingPPM#current = {setting}")


acquisition = client.get(f"{corrector}/Acquisition", context = context_acquisition).data["currentAverage"]
print(f"{corrector}/Acquisition#currentAverage = {acquisition}")

status = client.get(f"{corrector}/Status", context = context_empty).data
print(f"{corrector}/Status = {status}")
print("================================================================================")
print("Testing corrector CA.DHG0130, reading its values directly from interface, using get_correctors method... ")

result = I.get_correctors(names=["CA.DHG0130","CA.DVG0130" ])
print(f"names from method get_correctors: {result['names']}")
print(f"bdes  from method get_correctors: {result['bdes']}")
print(f"bact  from method get_correctors: {result['bact']}")
result2 = I.set_correctors(names=["CA.DHG0130", "CA.DVG0130"], corr_vals = [0.46, -0.81])
# result2 = I.set_correctors(names=["CA.DHG0130", "CA.DVG0130"], corr_vals = [0.0, -0.1])
result2 = I.get_correctors(names=["CA.DHG0130","CA.DVG0130" ])
print(f"names from method get_correctors: {result2['names']}")
print(f"bdes  from method get_correctors: {result2['bdes']}")
print(f"bact  from method get_correctors: {result2['bact']}")
print("================================================================================")
