import numpy as np
import obspy
from scipy.signal import correlate, welch
from scipy.stats import entropy

def spectral_entropy(trace, nperseg=256):
    """
    Computes the spectral entropy of a seismic trace.
    
    Parameters:
    trace (obspy.Trace): The seismic trace to analyze.
    nperseg (int): Number of samples per segment for the Welch method.
    
    Returns:
    float: Spectral entropy value (lower values indicate more periodic signals).
    """
    fs = trace.stats.sampling_rate
    f, Pxx = welch(trace.data, fs=fs, nperseg=nperseg)
    
    # Normalize power spectral density
    Pxx_norm = Pxx / np.sum(Pxx)
    
    # Compute spectral entropy
    spec_entropy = entropy(Pxx_norm)
    return spec_entropy

def compute_easy_snr(trace):
    """
    Computes the 'easy SNR' by comparing the loudest and quietest second.
    
    Parameters:
    trace (obspy.Trace): The seismic trace to analyze.
    
    Returns:
    float: The easy SNR value (higher means stronger event presence).
    """
    data = np.abs(trace.data)
    sample_rate = int(trace.stats.sampling_rate)
    num_samples = len(data)
    
    # Compute per-second median amplitudes
    second_medians = [np.median(data[i:i + sample_rate]) for i in range(0, num_samples, sample_rate)]
    
    if len(second_medians) < 2:
        return 0  # Not enough data to compare

    # Easy SNR is the ratio of the median of the loudest second to the median of the quietest second
    easy_snr = np.max(second_medians) / (np.min(second_medians) + 1e-10)  # Avoid division by zero
    return easy_snr

def identify_low_variability_trace(trace, std_threshold=0.02, autocorr_threshold=0.8, entropy_threshold=2.0):
    """
    Identifies a trace with little amplitude variation and quasi-repetitive signals.

    Parameters:
    trace (obspy.Trace): The seismic trace to analyze.
    std_threshold (float): Normalized standard deviation threshold for low variability.
    autocorr_threshold (float): Autocorrelation similarity threshold for repetition.
    entropy_threshold (float): Spectral entropy threshold for artificial signals.

    Returns:
    bool: True if the trace is likely an unresponsive station, False otherwise.
    """
    print(trace)
    data = trace.data.astype(float)
    
    # Compute standard deviation and mean absolute amplitude
    std_dev = np.std(data)
    mean_amp = np.mean(np.abs(data)) + 1e-10  # Avoid division by zero
    
    # Normalized variability metric
    variability_ratio = std_dev / mean_amp
    print(variability_ratio)

    # Compute autocorrelation (normalized)
    autocorr = correlate(data, data, mode='full')  # Full autocorrelation
    autocorr /= np.max(np.abs(autocorr))  # Normalize
    mid_idx = len(autocorr) // 2
    max_autocorr = np.max(autocorr[mid_idx + 1 :])  # Peak outside the central peak
    print(max_autocorr)

    # Compute spectral entropy
    spec_entropy = spectral_entropy(trace)
    print(spec_entropy)

    print((variability_ratio * spec_entropy)*(1.0-max_autocorr))

    # Compute easy SNR
    easy_snr = compute_easy_snr(trace)
    print(easy_snr)
    print((easy_snr * variability_ratio * spec_entropy)*(1.0-max_autocorr))
    print('\n')

    # Decision logic
    if variability_ratio < std_threshold and max_autocorr > autocorr_threshold and spec_entropy < entropy_threshold:
        return True  # This is likely a non-responsive station
    return False  # This trace has expected seismic behavior

# Example Usage:
st = obspy.read("/shares/hal9000/share/data/ASL_DB/20010301T141727/stream.mseed")
#st = obspy.read("/shares/hal9000/share/data/ASL_DB/20010301T150002/stream.mseed")

for trace in st:
    if identify_low_variability_trace(trace, std_threshold=0.5, autocorr_threshold=0.2, entropy_threshold=2.0):
        print(f"Station {trace.id} is likely malfunctioning.")
