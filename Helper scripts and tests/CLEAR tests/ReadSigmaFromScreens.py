import numpy as np
from Backend.State import State
from Interfaces.CLEAR.InterfaceCLEAR import CLEAR_real_machine

state = State(filename="/Users/wiktoriamalek/CERN-Flight_Simulator-Data/CLEAR_EM_260821/Quadrupole_scan_CLEAR_CA.QFD0350_20260821_190813/states_CA.QFD0350_20260821_190813/screen_0000_step_0001_shot_0009.pkl")
screen = state.get_screens()

image = np.asarray(screen["images"][0], dtype=float)
pixel_x_mm = float(screen["hpixel"][0])
pixel_y_mm = float(screen["vpixel"][0])

image = np.flipud(image)
ny, nx = image.shape
x_mm = pixel_x_mm * np.linspace(-nx / 2, nx / 2, nx)
y_mm = pixel_y_mm * np.linspace(-ny / 2, ny / 2, ny)

_, sigx = CLEAR_real_machine._fit_projection(x_mm, image.mean(axis=0))
_, sigy = CLEAR_real_machine._fit_projection(y_mm, image.mean(axis=1))

print(f"Gaussian σx = {sigx:.3f} mm")
print(f"Gaussian σy = {sigy:.3f} mm")