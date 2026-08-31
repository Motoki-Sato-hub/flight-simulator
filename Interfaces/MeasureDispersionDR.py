from epics import PV, ca, caget
from Interfaces.ATF2.InterfaceATF2_DR import InterfaceATF2_DR
import numpy as np
import sys
from matplotlib import pyplot as plt
I = InterfaceATF2_DR(nsamples= 3)

# Reading initial BPM readings
state = I.get_state()
bpms0 = state.get_orbit()

plt.clf()
plt.ion()
plt.plot(bpms0['x'], label = "x")
plt.plot(bpms0['y'], label = "y")
sys.exit(0)

# Changing the energy
dP_P = I.change_energy()

# Reading BPM readings after changing the energy
bmps1 = I.get_bpms()

# Resetting the energy
I.reset_energy()

x0 = np.mean(bmps0["x"], axis=0)
y0 = np.mean(bmps0["y"], axis=0)

x1 = np.mean(bmps1["x"], axis=0)
y1 = np.mean(bmps1["y"], axis=0)

Dx = (x1 - x0) / dP_P
Dy = (y1 - y0) / dP_P

for bpm, dx, dy in zip(bpms0["names"], Dx, Dy):
    print(f"BPM: {bpm}, Dx: {dx}, Dy: {dy}")
print("Arc dispersion =", I.get_arc_dispersion())
