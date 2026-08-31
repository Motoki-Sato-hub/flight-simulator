import os, sys, h5py
from pathlib import Path
import matplotlib
import numpy as np
matplotlib.use("QtAgg")
import matplotlib.pyplot as plt
project_root_path = Path.cwd().resolve()
while not (project_root_path / "Interfaces").exists() and project_root_path.parent != project_root_path:
    project_root_path = project_root_path.parent
sys.path.insert(0, str(project_root_path))
os.chdir(project_root_path)
from Interfaces.CLEAR.Setup_files.CLEAR_BPM_getHV import change_inverted_bpm_polarity, process_bpm_signal
path_to_datafiles = "/Users/wiktoriamalek/CERN-Flight_Simulator-Data/20260729_BPMtests"
bpms = ["BPM0530", "BPM0595", "BPM0690", "BPM0820", "BPM0890"]

'''Here you change the method from Sara's script'''
mode = "peak"

def analyze_data(filename, bpm, mode, show=False):
    file_path = os.path.join(path_to_datafiles, filename)
    with h5py.File(file_path, "r") as f:
        event_data = f["CLEAREventData"]
        H_SA = change_inverted_bpm_polarity(event_data[f"CA.{bpm}H-SA"]["SamplesFromTrigger"]["samples"][0], bpm)
        V_SA = change_inverted_bpm_polarity(event_data[f"CA.{bpm}V-SA"]["SamplesFromTrigger"]["samples"][0], bpm)
        S_SA = change_inverted_bpm_polarity(event_data[f"CA.{bpm}S-SA"]["SamplesFromTrigger"]["samples"][0], bpm)

    H = process_bpm_signal(H_SA, mode)
    V = process_bpm_signal(V_SA, mode)
    S = process_bpm_signal(S_SA, mode)
    print(f"{bpm}: H = {H}, V = {V}")

    if show:
        plt.figure()
        plt.plot(H_SA, label="H")
        plt.plot(V_SA, label="V")
        plt.plot(S_SA, label="S")
        plt.xlabel("Sample number")
        plt.ylabel("Signal amplitude [mV]")
        plt.title(f"{bpm}: H = {H}, V = {V}")
        plt.legend()
        plt.show()

    return {"file": filename, "bpm": bpm, "H": H, "V": V, "S": S}

list_of_files = sorted(filename for filename in os.listdir(path_to_datafiles) if filename.endswith(".h5"))
all_results = []

for orbit_number, filename in enumerate(list_of_files, start=1):
    print(f"===================== ORBIT {orbit_number} =====================")
    for bpm in bpms:
        all_results.append(analyze_data(filename, bpm, mode=mode))

M_h = np.array([
    [next(result["H"] for result in all_results if result["bpm"] == bpm and result["file"]==filename)
                 for bpm in bpms]
                for filename in list_of_files])

M_v = np.array([
    [next(result["V"] for result in all_results if result["bpm"] == bpm and result["file"]==filename)
                 for bpm in bpms]
                for filename in list_of_files])

M_s = np.array([
    [next(result["S"] for result in all_results if result["bpm"] == bpm and result["file"]==filename)
                 for bpm in bpms]
                for filename in list_of_files])

M_raw = np.hstack((M_h, M_v, M_s))
M = M_raw - M_raw.mean(axis=0, keepdims=True)
channel_labels = [f"{bpm} H" for bpm in bpms] + [f"{bpm} V" for bpm in bpms] + [f"{bpm} S" for bpm in bpms]

U, singular_values, Vt = np.linalg.svd(M, full_matrices=False)
variance_percent = singular_values**2 / np.sum(singular_values**2) * 100

print("Joint H+V+S singular values:", singular_values)
print("Joint H+V+S explained variance [%]:", variance_percent)
print("First joint mode:", Vt[0])
print("Second joint mode:", Vt[1])

for coordinate, matrix in (("H", M_h), ("V", M_v), ("S", M_s)):
    plt.figure()
    for bpm_index, bpm in enumerate(bpms):
        plt.plot(np.arange(1, len(list_of_files) + 1), matrix[:, bpm_index], label=bpm)
    plt.title(f"{coordinate} values [mV]")
    plt.xlabel("Event number")
    plt.ylabel("Signal amplitude [mV]")
    plt.legend()
    plt.tight_layout()
    plt.show()

fig = plt.figure()
ax = fig.add_subplot(projection="3d")
n_events, n_channels = M.shape
X, Y = np.meshgrid(np.arange(n_channels), np.arange(n_events))
ax.plot_surface(X, Y, M, cmap="viridis")
ax.set_title(r"Joint matrix with data $M = [M_H\;M_V\;M_S]$")
ax.set_xlabel("BPM coordinate")
ax.set_ylabel("Orbit number")
ax.set_zlabel("Centred signal [mV]")
ax.set_xticks(np.arange(n_channels))
ax.set_xticklabels(channel_labels, rotation=60, ha="right")
plt.tight_layout()
plt.show()

std_h = np.std(M[:, :5], axis=0)
std_v = np.std(M[:, 5:10], axis=0)
std_s = np.std(M[:, 10:15], axis=0)
plt.figure()
x = np.arange(len(bpms))
plt.plot(x, std_h, "o-", label="H")
plt.plot(x, std_v, "o-", label="V")
plt.plot(x, std_s, "o-", label="S")

plt.xticks(x, bpms)
plt.xlabel("BPM")
plt.ylabel("Std [mV]")
plt.title("Standard deviation of joint-matrix M (H+V+S)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

fig = plt.figure()
ax = fig.add_subplot(projection="3d")
n_events, n_modes = U.shape
X, Y = np.meshgrid(np.arange(1, n_modes + 1), np.arange(1, n_events + 1))
ax.plot_surface(X, Y, U, cmap="viridis")
ax.set_title(r"$U$ from joint H+V+S SVD")
ax.set_xlabel("SVD mode")
ax.set_ylabel("Event number")
ax.set_zlabel("Coefficient")
plt.tight_layout()
plt.show()

plt.figure()
event_numbers = np.arange(1, U.shape[0] + 1)
plt.plot(event_numbers, U[:, 0], label="mode 1")
plt.plot(event_numbers, U[:, 1], label="mode 2")
plt.title(r"First two temporal modes of joint SVD")
plt.xlabel("Event number")
plt.ylabel(r"$U[:, k]$")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

fig = plt.figure()
ax = fig.add_subplot(projection="3d")
n_modes, n_channels = Vt.shape
X, Y = np.meshgrid(np.arange(n_channels), np.arange(1, n_modes + 1))
ax.plot_surface(X, Y, Vt, cmap="viridis")
ax.set_title(r"$V^T$ from joint H+V+S SVD")
ax.set_xlabel("BPM coordinate (H, then V)")
ax.set_ylabel("SVD mode")
ax.set_zlabel("Mode coefficient")
ax.set_xticks(np.arange(n_channels))
ax.set_xticklabels(channel_labels, rotation=60, ha="right")
plt.tight_layout()
plt.show()

plt.figure()
x = np.arange(len(channel_labels))
plt.plot(x, Vt[0], "o-", label="mode 1")
plt.plot(x, Vt[1], "o-", label="mode 2")
plt.xticks(x, channel_labels, rotation=45, ha="right")
plt.axhline(0, color="black", linewidth=1)
plt.xlabel("BPM coordinate")
plt.ylabel("Mode coefficient")
plt.title(r"First 2 modes $V^T$ of joint SVD")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

plt.figure()
x = np.arange(len(channel_labels))
plt.plot(x, Vt[2], "o-", label="mode 3")
plt.plot(x, Vt[3], "o-", label="mode 4")
plt.xticks(x, channel_labels, rotation=45, ha="right")
plt.axhline(0, color="black", linewidth=1)
plt.xlabel("BPM coordinate")
plt.ylabel("Mode coefficient")
plt.title(r"3rd and 4th mode $V^T$ of joint SVD")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

plt.figure()
modes = np.arange(1, len(singular_values) + 1)
plt.plot(modes, singular_values, "o-")
plt.xlabel("SVD mode")
plt.ylabel("Singular value")
plt.title("Singular values of joint H+V+S matrix")
plt.xticks(modes)
plt.grid(True)
plt.tight_layout()
plt.show()

plt.figure()
plt.bar(modes, variance_percent)
plt.xlabel("SVD mode")
plt.ylabel("Variance [%]")
plt.title("Variance of joint H+V+S SVD")
plt.xticks(modes)
plt.grid(axis="y")
plt.tight_layout()
plt.show()

plt.figure()
for bpm in bpms:
    bpm_results = [result for result in all_results if result["bpm"] == bpm]
    s_values = np.asarray([result["S"] for result in bpm_results], dtype=float)
    v_values = np.asarray([result["V"] for result in bpm_results], dtype=float)
    plt.scatter(s_values, v_values, label=bpm)
plt.xlabel("S [mV]")
plt.ylabel("V [mV]")
plt.title("Correlation between S and V")
plt.legend()
plt.tight_layout()
plt.show()

plt.figure()
for bpm in bpms:
    bpm_results = [result for result in all_results if result["bpm"] == bpm]
    s_values = np.asarray([result["S"] for result in bpm_results], dtype=float)
    h_values = np.asarray([result["H"] for result in bpm_results], dtype=float)
    plt.scatter(s_values, h_values, label=bpm)

plt.xlabel("S [mV]")
plt.ylabel("H [mV]")
plt.title("Correlation between S and H")
plt.legend()
plt.tight_layout()
plt.show()
