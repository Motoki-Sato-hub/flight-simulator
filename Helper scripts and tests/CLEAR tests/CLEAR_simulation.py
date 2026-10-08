import csv, os, sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import RF_Track as rft
PROJECT_ROOT = Path.cwd().resolve()
while not (PROJECT_ROOT / "Interfaces").exists() and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
sys.path.insert(0, str(PROJECT_ROOT))
os.chdir(PROJECT_ROOT)
from Interfaces.CLEAR.InterfaceCLEAR_RFTrack import InterfaceCLEAR_RFTrack

interface = InterfaceCLEAR_RFTrack()
lattice = interface.lattice
lattice.set_tt_nsteps(1000)
lattice.track(interface.B0)

T = lattice.get_transport_table("%S %sigma_x %sigma_y %emitt_x %emitt_y %beta_x %beta_y %alpha_x %alpha_y")
s, sigma_x, sigma_y, emit_x, emit_y, beta_x, beta_y, alpha_x, alpha_y = T.T
gamma_x = (1.0 + alpha_x**2) / beta_x
gamma_y = (1.0 + alpha_y**2) / beta_y

def print_twiss_at(element_name: str) -> None:
    s_target = lattice[element_name].get_S("entrance")
    index = np.argmin(np.abs(s - s_target))
    print(f"\nTwiss parameters at {element_name} (s = {s_target:.4f} m)")
    print(f"  epsilon_x = {emit_x[index]:.4f} mm mrad")
    print(f"  beta_x    = {beta_x[index]:.4f} m")
    print(f"  alpha_x   = {alpha_x[index]:.4f}")
    print(f"  epsilon_y = {emit_y[index]:.4f} mm mrad")
    print(f"  beta_y    = {beta_y[index]:.4f} m")
    print(f"  alpha_y   = {alpha_y[index]:.4f}")


# The notebook reference Twiss values should be recovered at QFD350.
print_twiss_at("CA.QFD0350")
print_twiss_at("CA.QDD0515")
quad_s = lattice["CA.QDD0515"].get_S("entrance")

fig, axes = plt.subplots(4, 1, figsize=(11, 11), sharex=True)
ax_sigma, ax_beta, ax_alpha, ax_emit = axes

ax_sigma.plot(s, sigma_x, label=r"$\sigma_x$")
ax_sigma.plot(s, sigma_y, label=r"$\sigma_y$")
ax_sigma.set_ylabel("beam size [mm]")

ax_beta.plot(s, beta_x, label=r"$\beta_x$")
ax_beta.plot(s, beta_y, label=r"$\beta_y$")
ax_beta.set_ylabel(r"$\beta$ [m]")

ax_alpha.plot(s, alpha_x, label=r"$\alpha_x$")
ax_alpha.plot(s, alpha_y, label=r"$\alpha_y$")
ax_alpha.axhline(0.0, color="black", linewidth=0.8, alpha=0.5)
ax_alpha.set_ylabel(r"$\alpha$")

ax_emit.plot(s, emit_x, label=r"$\varepsilon_x$")
ax_emit.plot(s, emit_y, label=r"$\varepsilon_y$")
ax_emit.set_ylabel(r"$\varepsilon_n$ [mm mrad]")
ax_emit.set_xlabel("s [m]")
fig.suptitle("CLEAR beamline reconstruction - sigma_t = 0, sigma_pt = 0")

# for axis in axes:
#     axis.grid(alpha=0.3)
#     axis.legend(loc="best")
#     axis.axvline(quad_s, color="crimson", linestyle="--", linewidth=1.5)
# fig.suptitle(f"CLEAR RF-Track simulation")
# fig.tight_layout()

plt.show()
