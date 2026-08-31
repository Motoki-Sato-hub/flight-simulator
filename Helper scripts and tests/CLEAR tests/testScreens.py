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
#
# screens = I.get_screens(names=["CA.BTV0390"])
# print(screens)

"""Test of one screen and method get_screens"""

# print("Testing screen CA.BTV0390L, reading its values directly from japc/pyda... ")
# camera = client.get("CA.BTV0390.DigiCam/LastImage", context=context_empty).data
#
# image = camera["image2D"]
# print(f"image from PyDa = {image}")
#
# hpixel = camera["pixelCalSet1"]
# print(f"hpixel from PyDa = {hpixel}")
#
# vpixel = camera["pixelCalSet2"]
# print(f"vpixel from PyDa = {vpixel}")

# proj_x = camera["projDataSet1"]
# print(f"proj_x.shape = {proj_x.shape}")
#
# proj_y = camera["projDataSet2"]
# print(f"proj_y.shape = {proj_y.shape}")

# x_positions = camera["imagePositionSet1"]
# print(f"x_positions.shape = {x_positions.shape}")
#
# y_positions = camera["imagePositionSet2"]
# print(f"y_positions.shape = {y_positions.shape}")
#
# inout = client.get("CA.BTV0390_CAS.BTV0420/OPSettingSystem1", context = context_empty).data["positionChannel1"]
# print(f"inout from PyDa = {inout}")

print("================================================================================")
print("Testing screen CA.BTV0390L, reading its values directly from interface, using get_screens method... ")
# result = I.get_screens(names=["CA.BTV0390L"])
# print(f"names from get_screens: {result['names']}")
# print(f"hpixel from get_screens: {result['hpixel']}")
# print(f"vpixel from get_screens: {result['vpixel']}")
# print(f"x from get_screens: {result['x']}")
# print(f"y from get_screens: {result['y']}")
# print(f"sigx from get_screens:  {result['sigx']}")
# print(f"sigy from get_screens:  {result['sigy']}")
# print(f"sum from get_screens: {result['sum']}")
# print(f"images from get_screens: {result['images']}")
# print(f"inout from get_screens: {result['inout']}")
hpixel, vpixel = I._get_screen_pixel_calibration("CA.BTV0390L")
print(f"hpixel = {hpixel} mm/px")
print(f"vpixel = {vpixel} mm/px")
# r =client.get('CS.BTV0120_CS.BTV0305/OPSettingSystem2').data['positionChannel5']
# print(r.value) # prints just a value
#                 # screen out -> value=0
#                 # screen in -> value - 1
# print("=================================")
# r2 =client.get('CS.BTV0120_CS.BTV0305/Description').data['dcm3DriverNames']
# print(r2) # prints : ['Screen OUT', 'Screen IN']
#
# print("================================================================================")

# print("Testing inserting the screen using the interface method...")
# I.extract_screen("BTV0215")
# time.sleep(10)
#
#
# client.set(f"CS.BTV0305/OPSettingSystem2", data={f"info['get_set_field']": 1})
#
#
#
# client.set('CS.BTV0120_CS.BTV0305/OPSettingSystem2', data={"positionChannel5":1})
# r3 =client.get('CS.BTV0120_CS.BTV0305/Description').data['dcm3DriverNames']
# print(r3)





