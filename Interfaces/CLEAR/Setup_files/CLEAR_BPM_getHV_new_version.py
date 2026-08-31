# Initailly made for Wiktoria and Andrea to be able to test their flight simulator at CLEAR
# Sara S. July 2026

import numpy as np
import matplotlib.pyplot as plt
from pyjapc import PyJapc
from scipy.integrate import trapezoid


# Subtract baseline from signal
def baseline_correct(samples, n_baseline=100):
    baseline = np.mean(samples[:n_baseline])
    corrected_samples = samples - baseline
    return corrected_samples

# Peak seek
def find_peak(samples):
    idx = np.argmax(np.abs(samples))
    return samples[idx], idx

# 5% of peak integration threshold region
def threshold_integral(samples, threshold_fraction=0.05):
    # Find the peak (largest magnitude)
    peak_idx = np.argmax(np.abs(samples))
    peak = samples[peak_idx]
    pulse_sign = np.sign(peak)

    # Threshold relative to the peak
    threshold = threshold_fraction * abs(peak)

    # Search backwards for the start
    start = peak_idx
    while start > 0 and pulse_sign * samples[start] > threshold:
        start -= 1

    # Search forwards for the end
    end = peak_idx
    while end < len(samples) - 1 and pulse_sign * samples[end] > threshold:
        end += 1

    # Integrate the pulse
    integral = trapezoid(samples[start:end + 1])

    return integral, start, end, peak_idx

def plot_peak(signals, indices, peaks, labels, BPM):

    fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

    for axis, signal, idx, peak, label in zip(ax, signals, indices, peaks, labels):

        axis.plot(signal, label=f"{label} signal")
        axis.plot(idx, peak, "ro", label=f"Peak = {peak:.2f} mV")

        axis.set_title(f"{BPM} - {label}")
        axis.set_ylabel("Amplitude [mV]")
        axis.grid(True)
        axis.legend()

    ax[-1].set_xlabel("Sample")
    plt.tight_layout()
    plt.show()

def plot_integral(signals, integrals, labels, BPM, starts=None,ends=None,):
    fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

    for i, (axis, signal, integral, label) in enumerate(zip(ax, signals, integrals, labels)):

        x = np.arange(len(signal))
        axis.plot(x, signal, label=f"{label} signal")

        if starts is None:
            # Shade the whole signal
            axis.fill_between(x, signal, 0, alpha=0.3, label=f"Integral = {integral:.2f}",)

        else:
            start = starts[i]
            end = ends[i]
            axis.fill_between(x[start:end + 1], signal[start:end + 1], 0, alpha=0.3, label=f"Integral = {integral:.2f}",)
            axis.axvline(start, color="green", linestyle="--", label="Start")
            axis.axvline(end, color="red", linestyle="--", label="End")

        axis.set_title(f"{BPM} - {label}")
        axis.set_ylabel("Amplitude [mV]")
        axis.grid(True)
        axis.legend()

    ax[-1].set_xlabel("Sample")
    plt.tight_layout()
    plt.show()


scaling_factors = {
    "BPM0530": {"H": -0.403, "V": -0.396},
    "BPM0595": {"H": -0.351, "V": -0.374},
    "BPM0690": {"H": -0.400, "V": -0.410},
    "BPM0820": {"H": -0.346, "V": -0.392},
    "BPM0890": {"H": -0.391, "V": -0.417},
}


japc = PyJapc(incaAcceleratorName='CTF', selector="SCT.USER.SETUP")
DEFAULT_WINDOW = (260, 360)

def get_bpm_hv(BPM, mode, plot=False, window=DEFAULT_WINDOW,):
    
    """
    Acquire and process the horizontal (H) and vertical (V) BPM signals.

    Parameters
    ----------
    BPM : str
        BPM name: "BPM0530", "BPM0595", "BPM0690", "BPM0820", "BPM0890".
    mode : str
        - "peak"
            Return the peak value of the raw signal.
        - "baseline_peak"
            Return the peak value after baseline correction.
        - "integral"
            Return the integral of the baseline-corrected signal.
        - "integral_window"
            Return the integral of the baseline-corrected signal over the
            specified integration window, [min=0, max=730]. If not specified, 
            running with default window.
        - "integral_threshold"
            Return the integral of the baseline-corrected signal between the
            5% peak threshold crossings.
    plot : bool, optional
        Display processing plots.
        Default is False.

    Returns
    -------
    H, V
        Processed H and V values.
    """
    
    print(f"Reading {BPM}")
    H_data = japc.getParam(f"CA.{BPM}H-SA/SamplesFromTrigger")
    V_data = japc.getParam(f"CA.{BPM}V-SA/SamplesFromTrigger")
    S_data = japc.getParam(f"CA.{BPM}S-SA/SamplesFromTrigger")
    H_samples = H_data["samples"]
    V_samples = V_data["samples"]
    S_samples = S_data["samples"]
    H_b_samples = baseline_correct(H_samples)
    V_b_samples = baseline_correct(V_samples)
    S_b_samples = baseline_correct(S_samples)


    # Find peak (largest magnitude, keeping the sign)
    if mode == "peak":
        H, H_idx = find_peak(H_samples)
        V, V_idx = find_peak(V_samples)

        if plot: # Plot the signals
                plot_peak([H_samples, V_samples], [H_idx, V_idx], [H, V], ["H", "V"], BPM,)


    # Find peak on baseline corrected signal (largest magnitude, keeping the sign)
    elif mode == "baseline_peak":
        H, H_idx = find_peak(H_b_samples)
        V, V_idx = find_peak(V_b_samples)

        if plot: # Plot the signals
            plot_peak([H_b_samples, V_b_samples], [H_idx, V_idx], [H, V], ["H", "V"], BPM,)


    # Integration of full BPM signal with baseline correction
    elif mode == "integral":
        H = trapezoid(H_b_samples)
        V = trapezoid(V_b_samples)

        if plot:
            plot_integral(signals=[H_b_samples, V_b_samples], integrals=[H, V], labels=["H", "V"], BPM=BPM,)


    # Integration of fixed window BPM signal with baseline correction
    elif mode == "integral_window":
        window_start, window_end = window
        H = trapezoid(H_b_samples[window_start:window_end])
        V = trapezoid(V_b_samples[window_start:window_end])

        if plot:
            plot_integral(signals=[H_b_samples, V_b_samples], integrals=[H, V], labels=["H", "V"], BPM=BPM, starts=[window_start, window_start], ends=[window_end, window_end],)

    # Integration between 5% of the peak threshold region with baseline correction
    elif mode == "integral_threshold":
        H, H_start, H_end, H_peak_idx = threshold_integral(H_b_samples)
        V, V_start, V_end, V_peak_idx = threshold_integral(V_b_samples)
        S, S_start, S_end, S_peak_idx = threshold_integral(S_b_samples)

        H = (H/S) / scaling_factors[BPM]["H"]
        V = (V/S) / scaling_factors[BPM]["V"]

        if plot:
            plot_integral(signals=[H_b_samples, V_b_samples], integrals=[H, V], labels=["H", "V"], BPM=BPM, starts=[H_start, V_start], ends=[H_end, V_end],)

    else:
        raise ValueError(f"Unknown mode '{mode}'. " "Choose from: 'peak', 'baseline_peak', " "'integral', 'integral_window', 'integral_threshold'.")

    return H, V #S