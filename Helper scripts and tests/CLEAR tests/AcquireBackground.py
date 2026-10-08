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
from pprint import pprint

I = CLEAR_real_machine()

client = pyda.SimpleClient(provider=pyda_japc.JapcProvider())
context_acquisition = "SCT.USER.SETUP"
context_empty = ""

print("Getting background image...")
background_image = I.acquire_screen_background("CA.BTV0390", frames = 10)
print("Got background image!")
pprint(background_image)

print("Getting beam image, subtracted image and background image...")
subtracted_img, background_img, beam_img = I.acquire_screen_image("CA.BTV0390")
print("Got beam image, subtracted image and background image!")
pprint("Subtracted image:", subtracted_img)
pprint("Background image:", background_img)
pprint("Beam image:", beam_img)

print("Acquiring one raw beam image...")
I.insert_screen(screen_name)
camera_data = I._acquire_screen_data(screen_name)
beam_img = np.asarray(camera_data["image2D"], dtype=float)
print("Beam image data:", beam_img)


