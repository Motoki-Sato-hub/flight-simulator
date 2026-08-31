import os, sys, h5py, matplotlib
import numpy as np
from pprint import pprint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pathlib import Path
matplotlib.use('QtAgg')
project_root_path = Path.cwd().resolve()
while not (project_root_path / "Interfaces").exists() and project_root_path.parent != project_root_path:
    project_root_path = project_root_path.parent
sys.path.insert(0, str(project_root_path))
os.chdir(project_root_path)
from Interfaces.CLEAR.Setup_files.CLEAR_BPM_getHV import change_inverted_bpm_polarity, process_bpm_signal

path_to_datafiles = "/Users/wiktoriamalek/CERN-Flight_Simulator-Data/20260729_BPMtests"
bpms = ["BPM0530", "BPM0595", "BPM0690", "BPM0820", "BPM0890"]

def analyze_data(filename, bpm, mode):
    file_path = os.path.join(path_to_datafiles, filename)
    with h5py.File(file_path, "r") as f:

        H_SA = change_inverted_bpm_polarity(f["CLEAREventData"][f"CA.{bpm}H-SA"]["SamplesFromTrigger"]["samples"][0], bpm)
        V_SA = change_inverted_bpm_polarity(f["CLEAREventData"][f"CA.{bpm}V-SA"]["SamplesFromTrigger"]["samples"][0], bpm)
        S_SA = change_inverted_bpm_polarity(f["CLEAREventData"][f"CA.{bpm}S-SA"]["SamplesFromTrigger"]["samples"][0], bpm)

        # H = process_bpm_signal(H_SA, mode)
        # V = process_bpm_signal(V_SA, mode)
        # S = process_bpm_signal(S_SA, mode)
        H = np.sum(H_SA[320:330])
        V = np.sum(V_SA[320:330])
        S = np.sum(S_SA[320:330])

        print(f"{bpm}: H = {H}, V = {V}")

    return {
        "file": filename,
        "bpm": bpm,
        "H": H,
        "V": V,
        "S": S
    }

list_to_datafiles = "/Users/wiktoriamalek/CERN-Flight_Simulator-Data/20260729_BPMtests"
list_of_files = sorted(filename for filename in os.listdir(list_to_datafiles) if filename.endswith(".h5"))
all_results = []

'''Here you change the method from Sara's script'''
# done: baseline_peak, integral, integral_window, integral_threshold, peak
mode = "peak"

for orbit_number, filename in enumerate(list_of_files, start=1):
    print(f"=====================ORBIT {orbit_number}====================")
    for bpm in bpms:
        result = analyze_data(filename, bpm, mode=mode)
        all_results.append(result)

plt.figure()
for bpm in bpms:
    bpm_results = [result for result in all_results if result["bpm"] == bpm]
    H_values = [result["H"] for result in bpm_results]
    event_numbers = range(1, len(bpm_results) + 1)
    plt.plot(event_numbers, H_values, label=f"{bpm} H")
    plt.title("H values [mV]")
    plt.xlabel("Event number")
    plt.ylabel("Signal [mV]")
    plt.legend()
plt.show()


plt.figure()
for bpm in bpms:
    bpm_results = [result for result in all_results if result["bpm"] == bpm]
    V_values = [result["V"] for result in bpm_results]
    event_numbers = range(1, len(bpm_results) + 1)
    plt.plot(event_numbers, V_values, label=f"{bpm} V")
    plt.title("V values [mV]")
    plt.xlabel("Event number")
    plt.ylabel("Signal [mV]")
    plt.legend()
plt.show()

plt.figure()
for bpm in bpms:
    bpm_results = [result for result in all_results if result["bpm"] == bpm]
    S_values = [result["S"] for result in bpm_results]
    event_numbers = range(1, len(bpm_results) + 1)
    plt.plot(event_numbers, S_values, label=f"{bpm} S")
    plt.title("S values")
    plt.xlabel("Event number")
    plt.ylabel("Signal [mV]")
    plt.legend()
plt.show()

M_h = np.array([
    [next(result["H"] for result in all_results if result["bpm"] == bpm and result["file"]==filename)
                 for bpm in bpms]
                for filename in list_of_files])

M_v = np.array([
    [next(result["V"] for result in all_results if result["bpm"] == bpm and result["file"]==filename)
                 for bpm in bpms]
                for filename in list_of_files])

M_h = M_h - np.mean(M_h, axis=0)
U_h, singular_values_h, Vt_h = np.linalg.svd(M_h, full_matrices=False)
print(f"Singular values for H: {singular_values_h}")
variance_percent_h = singular_values_h**2 / np.sum(singular_values_h**2) * 100
print("Singular values for H in %: ", variance_percent_h)
print("First mode:", Vt_h[0])
print("Second mode:", Vt_h[1])

M_v = M_v - np.mean(M_v, axis=0)
U_v, singular_values_v, Vt_v = np.linalg.svd(M_v, full_matrices=False)
print(f"Singular values for V: {singular_values_v}")
variance_percent_v = singular_values_v**2 / np.sum(singular_values_v**2) * 100
print("Singular values for V in %: ", variance_percent_v)
print("First mode:", Vt_v[0])
print("Second mode:", Vt_v[1])

fig1 = plt.figure()
nh_events, nh_bpms = M_h.shape
ax1 = fig1.add_subplot(1,2,1, projection='3d')
ax2 = fig1.add_subplot(1,2,2, projection='3d')
X_h, Y_h = np.meshgrid(np.arange(nh_bpms), np.arange(nh_events))

ax1.plot_surface(X_h, Y_h, M_h, cmap='viridis')
ax1.set_title('$M_h - beam matrix$')
ax1.set_xlabel('BPM')
ax1.set_ylabel('Orbit [#]')
ax1.set_zlabel('H')

nv_events, nv_bpms = M_v.shape
X_v, Y_v = np.meshgrid(np.arange(nv_bpms), np.arange(nv_events))
ax2.plot_surface(X_v, Y_v, M_v, cmap='viridis')
ax2.set_title('$M_v - beam matrix$')
ax2.set_xlabel('BPM')
ax2.set_ylabel('Orbit [#]')
ax2.set_zlabel('V')
plt.show()

std_h = np.std(M_h, axis=0)
std_v = np.std(M_v, axis=0)
plt.figure()
x = np.arange(len(bpms))
width = 0.35
plt.plot(x, std_h, label='H')
plt.plot(x, std_v, label='V')
plt.xticks(x, bpms)
plt.xlabel("BPM")
plt.ylabel("Std [mV * sample]")
plt.title("Std of BPM signals")
plt.legend()
plt.show()

'''========U==========='''
fig = plt.figure()
ax1 = fig.add_subplot(1, 2, 1, projection="3d")
ax2 = fig.add_subplot(1, 2, 2, projection="3d")
n_events_h, n_modes_h = U_h.shape
X_h, Y_h = np.meshgrid(np.arange(1, n_modes_h + 1), np.arange(1, n_events_h + 1))
ax1.plot_surface(X_h, Y_h, U_h, cmap="viridis")
ax1.set_title(r"$U_h$ matrix from SVD")
ax1.set_xlabel("SVD mode")
ax1.set_ylabel("Event number")
ax1.set_zlabel(r"$U_h$ coefficient value")
n_events_v, n_modes_v = U_v.shape
X_v, Y_v = np.meshgrid(np.arange(1, n_modes_v + 1), np.arange(1, n_events_v + 1))
ax2.plot_surface(X_v, Y_v, U_v, cmap="viridis")
ax2.set_title(r"$U_v$ matrix from SVD")
ax2.set_xlabel("SVD mode")
ax2.set_ylabel("Event number")
ax2.set_zlabel(r"$U_v$ coefficient value")
plt.tight_layout()
plt.show()

fig, (ax1, ax2) = plt.subplots(1, 2)
event_numbers_h = np.arange(1, U_h.shape[0] + 1)
event_numbers_v = np.arange(1, U_v.shape[0] + 1)
ax1.plot(event_numbers_h, U_h[:, 0])
ax1.set_title(r"$U_h$ first column")
ax1.set_xlabel("Event number")
ax1.set_ylabel(r"$U_h[:, 0]$")
ax1.grid(True)
ax2.plot(event_numbers_v, U_v[:, 0])
ax2.set_title(r"$U_v$ first column")
ax2.set_xlabel("Event number")
ax2.set_ylabel(r"$U_v[:, 0]$")
ax2.grid(True)
plt.tight_layout()
plt.show()

fig = plt.figure()
ax1 = fig.add_subplot(1, 2, 1, projection="3d")
ax2 = fig.add_subplot(1, 2, 2, projection="3d")
n_events_h, n_modes_h = Vt_h.shape
X_h, Y_h = np.meshgrid(np.arange(1, n_modes_h + 1), np.arange(1, n_events_h + 1))
ax1.plot_surface(X_h, Y_h, Vt_h, cmap="viridis")
ax1.set_title(r"$Vt_h$ matrix from SVD")
ax1.set_xlabel("SVD mode")
ax1.set_ylabel("BPM")
ax1.set_zlabel(r"$Vt_h$ coefficient value")

n_events_v, n_modes_v = Vt_v.shape
X_v, Y_v = np.meshgrid(np.arange(1, n_modes_v + 1), np.arange(1, n_events_v + 1))
ax2.plot_surface(X_v, Y_v, Vt_v, cmap="viridis")
ax2.set_title(r"$Vt_v$ matrix from SVD")
ax2.set_xlabel("SVD mode")
ax2.set_ylabel("BPM")
ax2.set_zlabel(r"$Vt_v$ coefficient value")
plt.tight_layout()
plt.show()

fig, (ax1, ax2) = plt.subplots(1, 2)
x = np.arange(len(bpms))
ax1.plot(x, Vt_h[0], marker="o")
ax1.set_title(r"$V_h^T$ first row")
ax1.set_xlabel("BPM")
ax1.set_ylabel("Mode coefficient value")
ax1.set_xticks(x)
ax1.set_xticklabels(bpms)
ax1.axhline(0, linewidth=1)

ax2.plot(x, Vt_v[0], marker="o")
ax2.set_title(r"$V_v^T$ first row")
ax2.set_xlabel("BPM")
ax2.set_ylabel("Mode coefficient value")
ax2.set_xticks(x)
ax2.set_xticklabels(bpms)
ax2.axhline(0, linewidth=1)
plt.tight_layout()
plt.show()

plt.figure()
modes = np.arange(1, len(singular_values_v) + 1)
plt.plot(modes, singular_values_h, label='singular values H')
plt.plot(modes, singular_values_v, label='singular values V')
plt.xlabel("SVD mode")
plt.ylabel("Singular value")
plt.title("Singular values")
plt.xticks(modes)
plt.grid(True)
plt.legend()
plt.show()

mean_s = []
std_s = []
plt.figure()
x = np.arange(len(bpms))
plt.xticks(x, bpms)
for bpm in bpms:
    bpm_results = [result for result in all_results if result["bpm"] == bpm]
    s_values = np.asarray([result["S"] for result in bpm_results], dtype=float)
    mean_s.append(np.mean(s_values))
    std_s.append(np.std(s_values))
plt.errorbar(bpms, mean_s, yerr=std_s, fmt='o-', capsize = 5, linewidth=1.5)
plt.xlabel("BPM")
plt.ylabel("Mean S integral from 97 orbits [mV * sample]")
plt.title(r"Mean BPM total signal $\pm 1\sigma$")
plt.grid(True)
plt.show()
