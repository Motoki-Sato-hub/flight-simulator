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

"""Test of one quadrupole and method get_quadrupoles"""

print("Testing quadrupole CA.QFD0350, reading its values directly from japc/pyda... ")
quadrupole = "CA.QFD0350"
setting = client.get(f"{quadrupole}/SettingPPM", context=context_empty).data["current"]
print(f"{quadrupole}/SettingPPM#current = {setting}")
#client.set(f"{quadrupole}/SettingPPM", data={"current": 5.0})
setting2 = client.get(f"{quadrupole}/SettingPPM", context=context_empty).data["current"]
print(f"{quadrupole}/SettingPPM#current = {setting2}")
#client.set(f"{quadrupole}/SettingPPM", data={"current": 10.699})
setting2 = client.get(f"{quadrupole}/SettingPPM", context=context_empty).data["current"]
print(f"{quadrupole}/SettingPPM#current = {setting2}")


acquisition = client.get(f"{quadrupole}/Acquisition", context = context_acquisition).data["currentAverage"]
print(f"{quadrupole}/Acquisition#currentAverage = {acquisition}")

status = client.get(f"{quadrupole}/Status", context = "").data
print(f"{quadrupole}/Status = {status}")


print("================================================================================")
print("Testing quadrupole CA.QFD0350, reading its values directly from interface, using get_quadrupoles method... ")

result = I.get_quadrupoles(names=["CA.QFD0350"])
print(f"names from method get_quadrupoles: {result['names']}")
print(f"bdes  from method get_quadrupoles: {result['bdes']}")
print(f"bact  from method get_quadrupoles: {result['bact']}")

result2= I.set_quadrupoles(names=["CA.QFD0350"], currents_A=10.699)
result = I.get_quadrupoles(names=["CA.QFD0350"])
print(f"names from method set_quadrupoles: {result['names']}")
print(f"bdes  from method set_quadrupoles: {result['bdes']}")
print(f"bact  from method set_quadrupoles: {result['bact']}")

print("================================================================================")


