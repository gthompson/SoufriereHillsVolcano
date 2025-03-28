#!/usr/bin/env python
#import sys
#sys.path.append('/Users/thompsong/src/kitchensinkGT/LIB')
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import describe
from scipy.interpolate import interp1d
from obspy import Stream, Trace
from libseisGT import add_to_trace_history
import pandas as pd
from math import pi
from scipy.signal import savgol_filter
from obspy.core.event import Catalog, Event, Origin, Magnitude, Amplitude, StationMagnitude, CreationInfo, Comment, EventDescription
from obspy.core.event.resourceid import ResourceIdentifier
from obspy.core.event.event import EventType
from obspy.core.event.waveform import WaveformStreamID
from obspy import UTCDateTime, read
from collections import defaultdict
from obspy.geodetics.base import gps2dist_azimuth
"""
Functions for computing data quality metrics and statistical metrics (such as amplitude, energy and frequency) 
on Stream/Trace objects.

In terms of order of application:

1. Read the raw data.
2. Fix Trace IDs.
3. Compute QC metrics (and potentially remove bad Traces).
4. Correct the data (and save corrected data as MSEED/StationXML).
5. Compute statistical metrics.

"""
def estimate_snr(trace_or_stream, method='std', window_length=1.0, split_time=None,
                 verbose=False, spectral_kwargs=None, freq_band=None):
    """
    Estimate the signal-to-noise ratio (SNR) from a Trace or Stream.

    Parameters
    ----------
    trace_or_stream : obspy.Trace or obspy.Stream
        Input signal(s).
    method : str
        SNR method: 'max', 'std', 'rms', or 'spectral'.
    window_length : float
        Window length in seconds for window-based methods.
    split_time : UTCDateTime or tuple/list of UTCDateTime
        Time or (start, end) tuple to define signal/noise intervals.
    verbose : bool
        Print detailed output.
    spectral_kwargs : dict, optional
        Arguments passed to compute_amplitude_ratios().
    freq_band : tuple(float, float), optional
        Frequency range (Hz) to compute average SNR in 'spectral' mode.

    Returns
    -------
    snr, signal_val, noise_val : float
        SNR and underlying values.
    """
    if isinstance(trace_or_stream, Stream):
        return zip(*[
            estimate_snr(tr, method, window_length, split_time, verbose, spectral_kwargs, freq_band)
            for tr in trace_or_stream
        ])

    trace = trace_or_stream
    fs = trace.stats.sampling_rate

    if split_time:
        try:
            if isinstance(split_time, (list, tuple)) and len(split_time) == 2:
                signal_trace = trace.copy().trim(starttime=split_time[0], endtime=split_time[1], pad=True, fill_value=0)
                noise_trace = trace.copy().trim(endtime=split_time[0], pad=True, fill_value=0)
            else:
                noise_trace = trace.copy().trim(endtime=split_time, pad=True, fill_value=0)
                signal_trace = trace.copy().trim(starttime=split_time, pad=True, fill_value=0)

            signal_data = signal_trace.data
            noise_data = noise_trace.data

            if method == 'max':
                signal_val = np.nanmax(np.abs(signal_data))
                noise_val = np.nanmax(np.abs(noise_data))

            elif method in ('std', 'rms'):
                signal_val = np.nanstd(signal_data)
                noise_val = np.nanstd(noise_data)

            elif method == 'spectral':
                spectral_kwargs = spectral_kwargs or {}
                freqs, avg_ratio, *_ = compute_amplitude_ratios(
                    signal_trace, noise_trace, **spectral_kwargs, verbose=verbose
                )
                if avg_ratio is None:
                    snr = signal_val = noise_val = np.nan
                else:
                    if freq_band:
                        fmin, fmax = freq_band
                        band_mask = (freqs >= fmin) & (freqs <= fmax)
                        band_vals = avg_ratio[band_mask]
                        signal_val = np.nanmax(band_vals)
                        noise_val = np.nanmin(band_vals)
                        snr = np.nanmean(band_vals)
                    else:
                        signal_val = np.nanmax(avg_ratio)
                        noise_val = np.nanmin(avg_ratio)
                        snr = np.nanmean(avg_ratio)
            else:
                raise ValueError(f"Unknown SNR method: {method}")

            snr = snr if noise_val != 0 else np.inf

        except Exception as e:
            if verbose:
                print(f"Error during split-time SNR estimation: {e}")
            return np.nan, np.nan, np.nan

    else:
        data = trace.data
        npts = trace.stats.npts
        samples_per_window = int(fs * window_length)
        num_windows = int(npts // samples_per_window)

        if num_windows < 2:
            if verbose:
                print("Trace too short for SNR estimation.")
            return np.nan, np.nan, np.nan

        reshaped = data[:samples_per_window * num_windows].reshape((num_windows, samples_per_window))

        if method == 'max':
            values = np.nanmax(np.abs(reshaped), axis=1)
        elif method in ('std', 'rms'):
            values = np.nanstd(reshaped, axis=1)
        else:
            raise ValueError(f"Unknown SNR method: {method}")

        signal_val = np.max(values)
        noise_val = np.min(values)
        snr = signal_val / noise_val if noise_val != 0 else np.inf

    if 'metrics' not in trace.stats:
        trace.stats.metrics = {}
    trace.stats.metrics[f'snr_{method}'] = snr
    trace.stats.metrics[f'signal_level_{method}'] = signal_val
    trace.stats.metrics[f'noise_level_{method}'] = noise_val

    if verbose:
        print(f"[{method}] SNR = {snr:.2f} (signal={signal_val:.2f}, noise={noise_val:.2f})")

    return snr, signal_val, noise_val


def _check_spectral_qc(freqs, ratios, threshold=2.0, min_fraction_pass=0.5):
    """
    Check if a spectral ratio passes quality control.

    Parameters
    ----------
    freqs : np.ndarray
        Frequency bins.
    ratios : np.ndarray
        Spectral amplitude ratios.
    threshold : float
        Minimum acceptable amplitude ratio.
    min_fraction_pass : float
        Fraction of frequency bins that must exceed threshold.

    Returns
    -------
    passed : bool
        Whether the QC check passed.
    failed_fraction : float
        Fraction of frequencies below threshold.
    """    
    failed = np.sum(ratios < threshold)
    frac_failed = failed / len(ratios)
    return frac_failed <= (1 - min_fraction_pass), frac_failed

def compute_amplitude_ratios(signal_stream, noise_stream, smooth_window=None,
                              verbose=False, average='geometric',
                              qc_threshold=None, qc_fraction=0.5):
    """
    Compute spectral amplitude ratios between signal and noise traces.

    Parameters
    ----------
    signal_stream : obspy.Stream or obspy.Trace
        Signal trace(s). If a single Trace, it will be wrapped in a Stream.
    noise_stream : obspy.Stream or obspy.Trace
        Noise trace(s). If a single Trace, it will be wrapped in a Stream.
    smooth_window : int, optional
        Length of moving average window for smoothing the amplitude ratio.
    verbose : bool, optional
        If True, prints processing information.
    average : str, optional
        Method for computing the average spectral ratio. One of:
        'geometric', 'median', or 'mean'.
    qc_threshold : float, optional
        Minimum acceptable amplitude ratio for QC check.
    qc_fraction : float, optional
        Minimum fraction of frequency bins that must exceed `qc_threshold`
        for the QC to pass.

    Returns
    -------
    avg_freqs : np.ndarray or None
        Frequencies (Hz) corresponding to the average spectral ratio.
    avg_spectral_ratio : np.ndarray or None
        Averaged amplitude ratio spectrum (signal / noise).
    individual_ratios : list of np.ndarray
        List of spectral amplitude ratios for each matching trace pair.
    freqs_list : list of np.ndarray
        List of frequency arrays corresponding to each ratio.
    trace_ids : list of str
        Trace IDs for which amplitude ratios were computed.
    qc_results : dict
        Dictionary of QC results per trace. Each entry contains:
        {'passed': bool, 'failed_fraction': float}
    """
    if isinstance(signal_stream, Trace):
        signal_stream = Stream([signal_stream])
    if isinstance(noise_stream, Trace):
        noise_stream = Stream([noise_stream])

    noise_dict = {tr.id: tr for tr in noise_stream}
    individual_ratios = []
    freqs_list = []
    trace_ids = []
    qc_results = {}

    for sig_tr in signal_stream:
        trace_id = sig_tr.id
        if trace_id not in noise_dict:
            if verbose:
                print(f"Skipping {trace_id}: No matching noise trace.")
            continue

        noise_tr = noise_dict[trace_id]

        max_len = max(len(sig_tr.data), len(noise_tr.data))
        sig_data = np.pad(sig_tr.data, (0, max_len - len(sig_tr.data)), mode='constant')
        noise_data = np.pad(noise_tr.data, (0, max_len - len(noise_tr.data)), mode='constant')

        dt = sig_tr.stats.delta
        N = max_len

        fft_signal = np.fft.fft(sig_data)
        fft_noise = np.fft.fft(noise_data)
        freqs = np.fft.fftfreq(N, d=dt)

        amp_signal = np.abs(fft_signal)
        amp_noise = np.abs(fft_noise)
        amp_noise[amp_noise == 0] = 1e-10  # Avoid division by zero

        amplitude_ratio = amp_signal / amp_noise

        if smooth_window:
            kernel = np.ones(smooth_window) / smooth_window
            amplitude_ratio = np.convolve(amplitude_ratio, kernel, mode="same")

        ratio_half = amplitude_ratio[:N//2]
        freq_half = freqs[:N//2]

        individual_ratios.append(ratio_half)
        freqs_list.append(freq_half)
        trace_ids.append(trace_id)

        if qc_threshold is not None:
            passed, frac_failed = _check_spectral_qc(freq_half, ratio_half,
                                                     threshold=qc_threshold,
                                                     min_fraction_pass=qc_fraction)
            qc_results[trace_id] = {'passed': passed, 'failed_fraction': frac_failed}

    if not individual_ratios:
        return None, None, [], [], [], {}

    if average == 'geometric':
        avg_spectral_ratio = np.exp(np.mean(np.log(np.vstack(individual_ratios) + 1e-10), axis=0))
    elif average == 'median':
        avg_spectral_ratio = np.median(np.array(individual_ratios), axis=0)
    else:
        avg_spectral_ratio = np.mean(np.array(individual_ratios), axis=0)

    avg_freqs = freqs_list[0]

    if verbose:
        print(f"Computed amplitude ratios for {len(individual_ratios)} traces.")

    return avg_freqs, avg_spectral_ratio, individual_ratios, freqs_list, trace_ids, qc_results

def plot_amplitude_ratios(avg_freqs, avg_spectral_ratio,
                          individual_ratios=None, freqs_list=None, trace_ids=None,
                          log_scale=False, outfile=None, max_freq=50, threshold=None):
    """
    Plot spectral amplitude ratios (signal / noise), including optional individual traces.

    Parameters
    ----------
    avg_freqs : np.ndarray
        Frequency bins for the averaged spectral ratio.
    avg_spectral_ratio : np.ndarray
        Averaged amplitude ratio (signal / noise) spectrum.
    individual_ratios : list of np.ndarray, optional
        List of individual trace amplitude ratio spectra to overlay (default: None).
    freqs_list : list of np.ndarray, optional
        List of frequency bins corresponding to each entry in `individual_ratios`.
        Must match in length and shape (default: None).
    trace_ids : list of str, optional
        Labels for individual traces to include in the legend (default: index numbers).
    log_scale : bool, optional
        If True, plot the log10 of (amplitude ratio + 1) to better show weak signals (default: False).
    outfile : str, optional
        Path to save the figure instead of displaying it interactively (default: None).
    max_freq : float, optional
        Upper x-axis limit for frequency (Hz) (default: 50).
    threshold : float, optional
        If provided, plot a horizontal line at this amplitude ratio to indicate a QC threshold.

    Returns
    -------
    None
        Displays or saves the plot. Does not return anything.
    """


    plt.figure(figsize=(10, 6))

    if individual_ratios and freqs_list:
        for i, ratio in enumerate(individual_ratios):
            freqs = freqs_list[i]
            label = trace_ids[i] if trace_ids else f"Trace {i}"
            y = np.log10(ratio + 1) if log_scale else ratio
            plt.plot(freqs, y, label=label, alpha=0.5, linewidth=1)

    if avg_spectral_ratio is not None:
        y_avg = np.log10(avg_spectral_ratio + 1) if log_scale else avg_spectral_ratio
        plt.plot(avg_freqs, y_avg, color='black', linewidth=2.5, label='Average')

    if threshold:
        y_thresh = np.log10(threshold + 1) if log_scale else threshold
        plt.axhline(y_thresh, color='red', linestyle='--', label='SNR Threshold')

    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Log Amplitude Ratio" if log_scale else "Amplitude Ratio")
    plt.title("Amplitude Spectrum Ratio (Signal/Noise)")
    plt.xlim(0, max_freq)
    plt.ylim(bottom=0)
    plt.grid(True)
    plt.legend()

    if outfile:
        plt.savefig(outfile, bbox_inches='tight')
    else:
        plt.show()


def compute_amplitude_spectra(stream):
    """
    Compute amplitude spectra for all traces in a stream and store results in trace.stats.spectral.

    Parameters
    ----------
    stream : obspy.Stream
        Stream of Trace objects.

    Returns
    -------
    stream : obspy.Stream
        Same stream with each trace tagged with:
        - tr.stats.spectral['freqs']
        - tr.stats.spectral['amplitudes']
    """
    for tr in stream:
        dt = tr.stats.delta
        N = len(tr.data)
        fft_vals = np.fft.fft(tr.data)
        freqs = np.fft.fftfreq(N, d=dt)
        amp = np.abs(fft_vals)
        pos_freqs = freqs[:N//2]
        pos_amps = amp[:N//2]

        if not hasattr(tr.stats, 'spectral'):
            tr.stats.spectral = {}
        tr.stats.spectral['freqs'] = pos_freqs
        tr.stats.spectral['amplitudes'] = pos_amps

    return stream


def plot_amplitude_spectra(stream, max_freq=50, outfile=None):
    """
    Plot amplitude spectra for all traces in a stream, assuming they have spectral data.

    Parameters
    ----------
    stream : obspy.Stream
        Stream of Trace objects. Each trace must have .stats.spectral['freqs'] and ['amplitudes'].
    max_freq : float, optional
        Upper limit of x-axis (Hz). Default is 50 Hz.
    outfile : str, optional
        If provided, save the plot to this file instead of showing it.

    Returns
    -------
    None
    """
    import matplotlib.pyplot as plt

    plt.figure(figsize=(10, 6))

    for tr in stream:
        try:
            f = tr.stats.spectral['freqs']
            a = tr.stats.spectral['amplitudes']
            plt.plot(f, a, label=tr.id)
        except (AttributeError, KeyError):
            continue

    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Amplitude")
    plt.title("Amplitude Spectrum of Seismic Signals")
    plt.legend()
    plt.grid()
    plt.xlim(0, max_freq)

    if outfile:
        plt.savefig(outfile, bbox_inches='tight')
    else:
        plt.show()


def get_bandwidth(frequencies, amplitudes, threshold=0.707,
                  window_length=9, polyorder=2, trace=None):
    """
    Estimate peak frequency, bandwidth, and cutoff frequencies from
    a smoothed amplitude ratio spectrum. Optionally store results in
    trace.stats.metrics if a Trace is provided.

    Parameters
    ----------
    frequencies : np.ndarray
        Frequency values (Hz).
    amplitudes : np.ndarray
        Amplitude (or amplitude ratio) spectrum.
    threshold : float
        Fraction of peak amplitude to define bandwidth cutoff (e.g. 0.707 for -3 dB).
    window_length : int
        Smoothing window length for Savitzky-Golay filter.
    polyorder : int
        Polynomial order for Savitzky-Golay filter.
    trace : obspy.Trace, optional
        If given, store metrics in trace.stats.metrics.

    Returns
    -------
    dict
        Dictionary with keys: 'f_peak', 'A_peak', 'f_low', 'f_high', 'bandwidth'
    """
    smoothed = savgol_filter(amplitudes, window_length=window_length, polyorder=polyorder)

    peak_index = np.argmax(smoothed)
    f_peak = frequencies[peak_index]
    A_peak = smoothed[peak_index]
    cutoff_level = A_peak * threshold

    lower = np.where(smoothed[:peak_index] < cutoff_level)[0]
    f_low = frequencies[lower[-1]] if lower.size else frequencies[0]

    upper = np.where(smoothed[peak_index:] < cutoff_level)[0]
    f_high = frequencies[peak_index + upper[0]] if upper.size else frequencies[-1]

    bandwidth = f_high - f_low

    metrics = {
        "f_peak": f_peak,
        "A_peak": A_peak,
        "f_low": f_low,
        "f_high": f_high,
        "bandwidth": bandwidth
    }

    if trace is not None:
        if not hasattr(trace.stats, 'metrics') or not isinstance(trace.stats.metrics, dict):
            trace.stats.metrics = {}
        for key, val in metrics.items():
            trace.stats.metrics[key] = val

    return metrics

def ampengfft(stream, threshold=0.707, window_length=9, polyorder=2, differentiate=True):
    """
    Compute amplitude, energy, and frequency metrics for all traces in a stream.
    Note that we generally will want stream to be a displacement seismogram, and 
    then compute amplitude on this, and differentiate to compute energy. This is because
    we usually want to estimate a magnitude based on a displacement amplitude, but we 
    compute energy magnitude based on a velocity seismogram.
    However, this function by default will also then use displacement seismograms for
    all the frequency metrics, which effectively divides by frequency. If this is not
    the behaviour wanted, pass a velocity seismogram and pass differentiate=False.

    For each trace:
    - Computes FFT and amplitude spectrum.
    - Stores frequency and amplitude in trace.stats.spectral.
    - Computes bandwidth metrics using get_bandwidth().
    - Computes SSAM band-averaged amplitudes.
    - Computes dominant and mean frequency.
    - Computes band ratio metrics.
    - Computes amplitude, energy, and scipy.stats metrics.
    - Stores all results in trace.stats.metrics

    Parameters
    ----------
    stream : obspy.Stream
        Stream of Trace objects to analyze.
    threshold : float, optional
        Fraction of peak amplitude to define bandwidth cutoff (default is 0.707 for -3 dB).
    window_length : int, optional
        Length of Savitzky-Golay smoothing window (must be odd).
    polyorder : int, optional
        Polynomial order for smoothing filter.
    differentiate: bool, optional
        If stream contains displacement traces, then we want to differentiate to a velocity traces to compute energy

    Returns
    -------
    stream : obspy.Stream
        Modified stream with frequency metrics stored per trace.
    """
    if not isinstance(stream, Stream):
        if isinstance(stream, Trace):
            stream = Stream([stream])
        else:
            return
        
    if len(stream)==0:
        return
    
    if not 'detrended' in stream[0].stats.history: # change this to use better detrend from libseis?
        for tr in stream:
            tr.detrend(type='linear')
            add_to_trace_history(tr, 'detrended')    

    if not hasattr(stream[0].stats, 'spectral'):
        compute_amplitude_spectra(stream)

    for tr in stream:
        dt = tr.stats.delta
        N = len(tr.data)
        fft_vals = np.fft.fft(tr.data)
        freqs = np.fft.fftfreq(N, d=dt)
        amps = np.abs(fft_vals)

        pos_freqs = freqs[:N//2]
        pos_amps = amps[:N//2]

        if not hasattr(tr.stats, 'spectral'):
            tr.stats.spectral = {}
        tr.stats.spectral['freqs'] = pos_freqs
        tr.stats.spectral['amplitudes'] = pos_amps

        if not hasattr(tr.stats, 'metrics'):
            tr.stats.metrics = {}

        # Amplitude and energy metrics
        y = np.abs(tr.data)
        tr.stats.metrics['peakamp'] = np.max(y)
        tr.stats.metrics['peaktime'] = np.argmax(y) / tr.stats.sampling_rate + tr.stats.starttime
        if differentiate:
            y = np.diff(y)
        tr.stats.metrics['energy'] = np.sum(y**2) / tr.stats.sampling_rate

        # Magnitude - should be a separate after function
        # we need to convert tr.stats.metrics['energy'] to a source energy, Eseismic
        # based on source and station location, and an attenuation/decay model
        #tr.stats.metrics['Me'] = Eseismic2magnitude(Eseismic, correction=3.7)

        # Bandwidth
        try:
            get_bandwidth(pos_freqs, pos_amps,
                          threshold=threshold,
                          window_length=window_length,
                          polyorder=polyorder,
                          trace=tr)
        except Exception as e:
            print(f"[{tr.id}] Skipping bandwidth metrics: {e}")

        # SSAM
        try:
            _ssam(tr)
        except Exception as e:
            print(f"[{tr.id}] SSAM computation failed: {e}")

        # Band ratios
        try:
            _band_ratio(tr, freqlims=[1.0, 6.0, 11.0])
            _band_ratio(tr, freqlims=[0.5, 3.0, 18.0])
        except Exception as e:
            print(f"[{tr.id}] Band ratio computation failed: {e}")

        # Dominant and mean frequency
        try:
            f = tr.stats.spectral['freqs']
            A = tr.stats.spectral['amplitudes']
            tr.stats.metrics['peakf'] = f[np.argmax(A)]
            tr.stats.metrics['meanf'] = np.sum(f * A) / np.sum(A) if np.sum(A) > 0 else np.nan
        except Exception as e:
            print(f"[{tr.id}] Dominant/mean frequency computation failed: {e}")

        # Scipy descriptive stats
        try:
            from scipy.stats import describe
            stats = describe(tr.data, nan_policy='omit')._asdict()
            tr.stats.metrics['skewness'] = stats['skewness']
            tr.stats.metrics['kurtosis'] = stats['kurtosis']
        except Exception as e:
            print(f"[{tr.id}] scipy.stats failed: {e}")

        add_to_trace_history(tr, 'ampengfft')    

    return stream


def _ssam(tr, freq_bins=None):
    """
    Compute single-station amplitude measurements (SSAM) by binning
    the amplitude spectrum already stored in tr.stats.spectral.

    Parameters
    ----------
    tr : obspy.Trace
        Trace with .stats.spectral['freqs'] and ['amplitudes'].
    freq_bins : array-like, optional
        Frequency bin edges (e.g., np.arange(0, 16, 1.0)).
        Default is 0–16 Hz in 1-Hz bins.

    Returns
    -------
    None
        Stores {'f': bin_centers, 'A': band_averaged_amplitudes} in tr.stats.metrics.ssam
    """
    if freq_bins is None:
        freq_bins = np.arange(0.0, 16.0, 1.0)

    if not hasattr(tr.stats, 'spectral'):
        print(f"[{tr.id}] Missing tr.stats.spectral — run compute_amplitude_spectra() first.")
        return

    f = tr.stats.spectral.get('freqs')
    A = tr.stats.spectral.get('amplitudes')

    if f is None or A is None:
        print(f"[{tr.id}] Missing spectral frequencies or amplitudes.")
        return

    f = np.asarray(f)
    A = np.asarray(A)

    bin_centers = []
    ssam_values = []

    for i in range(len(freq_bins) - 1):
        fmin = freq_bins[i]
        fmax = freq_bins[i+1]
        idx = np.where((f >= fmin) & (f < fmax))[0]

        bin_centers.append((fmin + fmax) / 2.0)
        ssam_values.append(np.nanmean(A[idx]) if idx.size else np.nan)

    tr.stats.metrics['ssam'] = {
        'f': np.array(bin_centers),
        'A': np.array(ssam_values)
    }


    
def _band_ratio(tr, freqlims=[1, 6, 11]):
    """
    Compute band ratio as log2(amplitude above split frequency / amplitude below),
    using frequency limits defined by freqlims.

    Parameters
    ----------
    tr : obspy.Trace
        Trace object with spectral data.
    freqlims : list of float
        Frequency limits in Hz: [low, split, high]
        Some values used before are:
            [1, 6, 11]
            [0.8, 4, 18]
        Better values might be:
            [0.5, 4.0, 32.0]
            [0.5, 3.0, 18.0]
        Or choose two sets of values to better differentiate LPs from hybrids from VTs?
            [0.5, 2.0, 8.0]  LF vs HF
            [2.0, 6.0, 18.0] hybrid vs VT?
        Need to play around with data to determine     

    Stores
    -------
    tr.stats.metrics.bandratio : list of dict
        Appends a dictionary with 'freqlims', 'RSAM_low', 'RSAM_high', 'RSAM_ratio'.
    """
    A = None
    f = None

    # Preferred: new spectral storage
    if hasattr(tr.stats, 'spectral'):
        f = tr.stats.spectral.get('freqs')
        A = tr.stats.spectral.get('amplitudes')

    # Legacy support
    elif hasattr(tr.stats, 'spectrum'):
        f = tr.stats.spectrum.get('F')
        A = tr.stats.spectrum.get('A')

    elif hasattr(tr.stats, 'ssam'):
        f = tr.stats.ssam.get('f')
        A = np.array(tr.stats.ssam.get('A'))

    # Proceed if valid spectral data found
    if A is not None and f is not None and len(A) > 0:
        f = np.array(f)
        A = np.array(A)

        idx_low = np.where((f > freqlims[0]) & (f < freqlims[1]))[0]
        idx_high = np.where((f > freqlims[1]) & (f < freqlims[2]))[0]

        A_low = A[idx_low]
        A_high = A[idx_high]

        sum_low = np.sum(A_low)
        sum_high = np.sum(A_high)

        ratio = np.log2(sum_high / sum_low) if sum_low > 0 else np.nan

        br = {
            'freqlims': freqlims,
            'RSAM_low': sum_low,
            'RSAM_high': sum_high,
            'RSAM_ratio': ratio
        }

        if not hasattr(tr.stats.metrics, 'bandratio'):
            tr.stats.metrics.bandratio = []

        tr.stats.metrics.bandratio.append(br)




def export_trace_metrics_to_dataframe(stream, include_id=True, include_starttime=True):
    """
    Export all trace.stats.metrics and associated metadata into a pandas DataFrame.

    Parameters
    ----------
    stream : obspy.Stream
        Stream with populated trace.stats.metrics dictionaries.
    include_id : bool
        Whether to include trace.id as a column.
    include_starttime : bool
        Whether to include trace.stats.starttime as a column.

    Returns
    -------
    df : pandas.DataFrame
        DataFrame of metrics, one row per trace.
    """
    rows = []

    for tr in stream:
        s = tr.stats
        row = {}

        if include_id:
            row['id'] = tr.id
        if include_starttime:
            row['starttime'] = s.starttime

        row['Fs'] = s.sampling_rate
        row['calib'] = getattr(s, 'calib', None)
        row['units'] = getattr(s, 'units', None)
        row['quality'] = getattr(s, 'quality_factor', None)

        if hasattr(s, 'spectrum'):
            for item in ['medianF', 'peakF', 'peakA', 'bw_min', 'bw_max']:
                row[item] = s.spectrum.get(item, None)

        if hasattr(s, 'metrics'):
            m = s.metrics
            for item in [
                'snr', 'signal_level', 'noise_level', 'twin',
                'peakamp', 'peaktime', 'energy', 'RSAM_high', 'RSAM_low',
                'sample_min', 'sample_max', 'sample_mean', 'sample_median',
                'sample_lower_quartile', 'sample_upper_quartile', 'sample_rms',
                'sample_stdev', 'percent_availability', 'num_gaps', 'skewness', 'kurtosis'
            ]:
                row[item] = m.get(item, None)
            for key, value in m.items():
                if isinstance(value, dict):
                    for subkey, subval in value.items():
                        row[f"{key}_{subkey}"] = subval

        if 'bandratio' in m:
            for dictitem in m['bandratio']:
                label = 'bandratio_' + "".join(str(dictitem['freqlims'])).replace(', ', '_')
                row[label] = dictitem['RSAM_ratio']

        if 'lon' in s:
            row['lon'] = s['lon']
            row['lat'] = s['lat']
            row['elev'] = s['elev']

        rows.append(row)

    df = pd.DataFrame(rows)
    df = df.round({
        'Fs': 2, 'secs': 2, 'quality': 2, 'medianF': 1, 'peakF': 1,
        'bw_max': 1, 'bw_min': 1, 'peaktime': 2, 'twin': 2,
        'skewness': 2, 'kurtosis': 2
    })
    return df

        
def load_trace_metrics_from_dataframe(df, stream, match_on='id'):
    """
    Load metrics from a DataFrame into each trace.stats.metrics.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame with one row per trace and metrics as columns.
    stream : obspy.Stream
        Stream of Trace objects to update.
    match_on : str
        Column to match traces by ('id' or 'starttime').

    Returns
    -------
    stream : obspy.Stream
        Stream with updated .stats.metrics per trace.
    """
    for tr in stream:
        key = tr.id if match_on == 'id' else tr.stats.starttime
        row = df[df[match_on] == key]
        if row.empty:
            continue

        tr.stats.metrics = {}

        for col in row.columns:
            if col in ['id', 'starttime']:
                continue
            val = row.iloc[0][col]

            # Optional: split compound keys like 'ssam_A' back into nested dicts
            if '_' in col:
                main_key, sub_key = col.split('_', 1)
                if main_key not in tr.stats.metrics:
                    tr.stats.metrics[main_key] = {}
                tr.stats.metrics[main_key][sub_key] = val
            else:
                tr.stats.metrics[col] = val

    return stream

def save_enhanced_stream(st, enhanced_wavpath, save_pickle=False):
    """
    Save a stream to MiniSEED along with a CSV of enhanced metrics.

    Parameters
    ----------
    st : obspy.Stream
        Stream with .stats.metrics, .stats.spectrum, etc.
    enhanced_wavpath : str
        Path without extension (.mseed, .csv will be added).
    save_pickle : bool
        Whether to also save as ObsPy .pickle format (optional).
    """
    if enhanced_wavpath.endswith('.mseed'):
        enhanced_wavpath = enhanced_wavpath[:-6]

    os.makedirs(os.path.dirname(enhanced_wavpath), exist_ok=True)

    # Save MiniSEED
    st.write(enhanced_wavpath + '.mseed', format='MSEED')

    # Save metrics CSV
    df = export_trace_metrics_to_dataframe(st)
    df.to_csv(enhanced_wavpath + '.csv', index=False)

    # Optional: also save full stream with attributes
    if save_pickle:
        st.write(enhanced_wavpath + '.pickle', format='PICKLE')


def read_enhanced_stream(enhanced_wavpath):
    """
    Read a stream from MiniSEED + CSV metrics and restore enhanced trace.stats.

    Parameters
    ----------
    enhanced_wavpath : str
        Path without extension (.mseed, .csv will be added).

    Returns
    -------
    st : obspy.Stream
        Stream with restored .stats.metrics and other attributes.
    """
    if enhanced_wavpath.endswith('.mseed'):
        enhanced_wavpath = enhanced_wavpath[:-6]

    st = read(enhanced_wavpath + '.mseed')
    df = pd.read_csv(enhanced_wavpath + '.csv')

    # First pass: restore metrics
    st = load_trace_metrics_from_dataframe(df, st)

    # Second pass: restore extra fields not in .metrics
    for tr in st:
        row = df[df['id'] == tr.id]
        if row.empty:
            continue
        row = row.iloc[0]

        s = tr.stats
        s.units = row.get('units', None)
        s.calib = row.get('calib', None)
        s.quality_factor = row.get('quality', None)

        # Restore spectrum if present
        spectrum_keys = ['medianF', 'peakF', 'peakA', 'bw_min', 'bw_max']
        s.spectrum = {k: row[k] for k in spectrum_keys if k in row and not pd.isna(row[k])}

        # Restore bandratio (if any)
        bandratios = []
        for col in row.index:
            if col.startswith('bandratio_') and not pd.isna(row[col]):
                try:
                    freqlims_str = col.replace('bandratio_', '').replace('[', '').replace(']', '')
                    freqlims = [float(x) for x in freqlims_str.split('_')]
                    bandratios.append({'freqlims': freqlims, 'RSAM_ratio': row[col]})
                except Exception as e:
                    print(f"Warning: couldn't parse {col}: {e}")
        if bandratios:
            s.bandratio = bandratios

        for coord in ['lon', 'lat', 'elev']:
            if coord in row:
                s[coord] = row[coord]

    return st




def estimate_distance(trace, source_coords):
    """
    Compute hypocentral distance R (in meters) from trace coordinates to source.

    Parameters
    ----------
    trace : obspy.Trace
        Must have .stats.coordinates = {'latitude': ..., 'longitude': ..., 'elevation': ...}
    source_coords : dict
        {'latitude': ..., 'longitude': ..., 'depth': ...} in meters

    Returns
    -------
    R : float
        Hypocentral distance in meters
    """
    sta = trace.stats.coordinates
    lat1, lon1, elev = sta['latitude'], sta['longitude'], sta.get('elevation', 0)
    lat2, lon2, depth = source_coords['latitude'], source_coords['longitude'], source_coords.get('depth', 0)

    epic_dist, _, _ = gps2dist_azimuth(lat1, lon1, lat2, lon2)
    dz = (depth + elev)  # add elevation since depth is below surface
    return np.sqrt(epic_dist**2 + dz**2)



def estimate_source_energy(trace, R, model='body', Q=50, c_earth=2500):
    """
    Estimate source energy by correcting station energy for geometric spreading and attenuation.

    Parameters
    ----------
    trace : obspy.Trace
        Must have .stats.metrics['energy'] and .stats.spectral['freqs'] / 'peakF'
    R : float
        Distance in meters.
    model : str
        'body' or 'surface'
    Q : float
        Quality factor for attenuation.
    c_earth : float
        Seismic wave speed (m/s)

    Returns
    -------
    Eseismic : float
        Estimated source energy in Joules
    """
    E_obs = trace.stats.metrics.get('energy')
    if E_obs is None:
        return None

    # Use peakF if available
    f_peak = None
    if 'peakF' in trace.stats.spectral:
        f_peak = trace.stats.spectral['peakF']
    else:
        freqs = trace.stats.spectral.get('freqs')
        amps = trace.stats.spectral.get('amplitudes')
        if freqs is not None and amps is not None:
            f_peak = freqs[np.argmax(amps)]

    if f_peak is None:
        return None

    # Attenuation factor
    A_att = np.exp(-np.pi * f_peak * R / (Q * c_earth))

    # Geometric spreading correction
    if model == 'body':
        geom = R**2
    elif model == 'surface':
        wavelength = c_earth / f_peak
        geom = R * wavelength
    else:
        raise ValueError("Model must be 'body' or 'surface'.")

    Eseismic = E_obs * geom / A_att  # undo attenuation and spreading
    return Eseismic


def compute_station_magnitudes(stream, inventory, source_coords,
                                 model='body', Q=50, c_earth=2500, correction=3.7,
                                 a=1.6, b=-0.15, g=0,
                                 use_boatwright=True,
                                 rho_earth=2000, S=1.0, A=1.0,
                                 rho_atmos=1.2, c_atmos=340, z=100000,
                                 attach_coords=True, compute_distances=True):
    """
    Attach coordinates and estimate energy-based and local magnitudes for all traces.

    Parameters
    ----------
    stream : obspy.Stream
    inventory : obspy.Inventory
        Station metadata
    source_coords : dict
        {'latitude': ..., 'longitude': ..., 'depth': ...}
    model : str
        'body' or 'surface' wave geometric spreading
    Q : float
        Attenuation quality factor
    c_earth : float
        Seismic wave speed (m/s)
    correction : float
        Correction factor in Hanks & Kanamori formula
    a, b, g : float
        ML Richter coefficients and station correction
    use_boatwright : bool
        If True, use Boatwright formulation for energy estimation
    rho_earth, S, A : float
        Parameters for seismic Boatwright energy model
    rho_atmos, c_atmos, z : float
        Parameters for infrasound Boatwright energy model
    attach_coords : bool
        Whether to attach station coordinates from inventory
    compute_distances : bool
        Whether to compute and store distance (in meters)

    Returns
    -------
    stream : obspy.Stream
        Updated stream with .stats.metrics['energy_magnitude'] and ['local_magnitude']
    """
    if attach_coords:
        attach_station_coordinates_from_inventory(inventory, stream)

    for tr in stream:
        try:
            R = estimate_distance(tr, source_coords)

            if compute_distances:
                tr.stats['distance'] = R

            if not hasattr(tr.stats, 'metrics'):
                tr.stats.metrics = {}

            # Choose energy model based on data type and user preference
            if use_boatwright:
                if tr.stats.channel[1].upper() == 'D':  # Infrasound
                    E0 = Eacoustic_Boatwright(tr, R, rho_atmos=rho_atmos, c_atmos=c_atmos, z=z)
                else:  # Seismic
                    E0 = Eseismic_Boatwright(tr, R, rho_earth=rho_earth, c_earth=c_earth, S=S, A=A)

                    # Apply Q correction using spectral peak frequency if available
                    if 'spectral' in tr.stats and 'freqs' in tr.stats.spectral and 'amplitudes' in tr.stats.spectral:
                        freqs = tr.stats.spectral['freqs']
                        amps = tr.stats.spectral['amplitudes']
                        if np.any(amps > 0):
                            f_peak = freqs[np.argmax(amps)]
                            A_att = np.exp(-np.pi * f_peak * R / (Q * c_earth))
                            E0 /= A_att
            else:
                E0 = estimate_source_energy(tr, R, model=model, Q=Q, c_earth=c_earth)

            ME = Eseismic2magnitude(E0, correction=correction)

            tr.stats.metrics['source_energy'] = E0
            tr.stats.metrics['energy_magnitude'] = ME

            if tr.stats.channel[1].upper() in ('H', 'L'):
                R_km = R / 1000
                ML = estimate_local_magnitude(tr, R_km, a=a, b=b, g=g)
                tr.stats.metrics['local_magnitude'] = ML

        except Exception as e:
            print(f"[{tr.id}] Magnitude estimation failed: {e}")

    return stream


# could just import this from InventoryTools
def attach_station_coordinates_from_inventory(inventory, st):
    """ attach_station_coordinates_from_inventory """
    from obspy.core.util import AttribDict
    for tr in st:
        for netw in inventory.networks:
            for sta in netw.stations:
                if tr.stats.station == sta.code and netw.code == tr.stats.network:
                    for cha in sta.channels:
                        if tr.stats.location == cha.location_code:
                            tr.stats.coordinates = AttribDict({
                                'latitude':cha.latitude,
                                'longitude':cha.longitude,
                                'elevation':cha.elevation})
                            #tr.stats.latitude = cha.latitude
                            #tr.stats.longitude = cha.longitude  



def summarize_magnitudes(stream, include_network=True):
    """
    Summarize trace-level and network-averaged magnitude values.

    Parameters
    ----------
    stream : obspy.Stream
    include_network : bool
        If True, append a summary row with network-average magnitudes.

    Returns
    -------
    df : pandas.DataFrame
        Columns: id, starttime, distance, local_magnitude (ML), energy_magnitude (ME),
                 plus network_mean_* and network_std_* if include_network is True.
    """
    rows = []

    for tr in stream:
        row = {
            'id': tr.id,
            'starttime': tr.stats.starttime
        }

        row['distance_m'] = tr.stats.get('distance', np.nan)
        metrics = tr.stats.get('metrics', {})

        row['ML'] = metrics.get('local_magnitude', np.nan)
        row['ME'] = metrics.get('energy_magnitude', np.nan)
        row['source_energy'] = metrics.get('source_energy', np.nan)
        row['peakamp'] = metrics.get('peakamp', np.nan)
        row['energy'] = metrics.get('energy', np.nan)

        rows.append(row)

    df = pd.DataFrame(rows)

    if include_network:
        network_stats = {}
        for col in ['ML', 'ME']:
            valid = df[col].dropna()
            if not valid.empty:
                network_stats[f'network_mean_{col}'] = valid.mean()
                network_stats[f'network_mean_{col}'] = valid.median()
                network_stats[f'network_std_{col}'] = valid.std()
                network_stats[f'network_n_{col}'] = valid.count()

        # Add a summary row at the bottom
        summary_row = {col: np.nan for col in df.columns}
        summary_row.update(network_stats)
        df = pd.concat([df, pd.DataFrame([summary_row])], ignore_index=True)

    return df


def estimate_local_magnitude(trace, R_km, a=1.6, b=-0.15, g=0):
    """
    Estimate and store ML from peak amplitude and distance.

    Parameters
    ----------
    trace : obspy.Trace
    R_km : float
        Distance from source in km
    a, b, g : float
        Richter scaling parameters

    Returns
    -------
    ml : float
    """
    peakamp = trace.stats.metrics.get('peakamp')
    if peakamp is None or R_km <= 0:
        return None

    ml = Mlrichter(peakamp, R_km, a=a, b=b, g=g)
    trace.stats.metrics['local_magnitude'] = ml
    return ml

def estimate_network_magnitude(stream, key='local_magnitude'):
    """
    Compute network-average magnitude with uncertainty.

    Parameters
    ----------
    stream : obspy.Stream
    key : str
        Metric to average (e.g., 'local_magnitude', 'energy_magnitude')

    Returns
    -------
    tuple of (mean_mag, std_dev, n)
    """
    mags = []
    for tr in stream:
        mag = tr.stats.metrics.get(key)
        if mag is not None:
            mags.append(mag)

    if len(mags) == 0:
        return None, None, 0

    return np.mean(mags), np.std(mags), len(mags)


def Eseismic_Boatwright(val, R, rho_earth=2000, c_earth=2500, S=1.0, A=1.0):
    # val can be a Stream, Trace, a stationEnergy or a list of stationEnergy
    # R in m
    # Following values assumed by Johnson and Aster, 2005:
    # rho_earth 2000 kg/m^3
    # c_earth 2500 m/s
    # A is attenuation = 1
    # S is site response = 1
    #
    # These equations seem to be valid for body waves only, that spread like hemispherical waves in a flat earth.
    # But if surface waves dominate, they would spread like ripples on a pond, so energy density of wavefront like 2*pi*R
    if isinstance(val,Stream): # Stream
        Eseismic = []
        for tr in val:
            Eseismic.append(Eseismic_Boatwright(tr, R, rho_earth, c_earth, S, A))
    elif isinstance(val,Trace): # Trace
        stationEnergy = compute_stationEnergy(val) 
        Eseismic = Eseismic_Boatwright(stationEnergy, R, rho_earth, c_earth, S, A)
    elif isinstance(val, list): # list of stationEnergy
        Eseismic = []
        for thisval in val:
            Eseismic.append(Eseismic_Boatwright(thisval, R, rho_earth, c_earth, S, A))
    else: # stationEnergy
        Eseismic = 2 * pi * (R ** 2) * rho_earth * c_earth * (S ** 2) * val / A 
    return Eseismic

def Eacoustic_Boatwright(val, R, rho_atmos=1.2, c_atmos=340, z=100000):
    # val can be a Stream, Trace, a stationEnergy or a list of stationEnergy
    # R in m
    # Following values assumed by Johnson and Aster, 2005:
    # rho_atmos 1.2 kg/m^3
    # c_atmos 340 m/s  
    # z is just an estimate of the atmospheric vertical scale length - the height of ripples of infrasound energy spreading globally
    if isinstance(val,Stream): # Stream
        Eacoustic = []
        for tr in val:
            Eacoustic.append(Eacoustic_Boatwright(tr, R, rho_atmos, c_atmos))
    elif isinstance(val, Trace): # Trace
        stationEnergy = compute_stationEnergy(val) 
        Eacoustic = Eacoustic_Boatwright(stationEnergy, R, rho_atmos, c_atmos)
    elif isinstance(val, list): # list of stationEnergy
        for thisval in val:
            Eacoustic.append(Eacoustic_Boatwright(thisval, R, rho_atmos, c_atmos, S, A))
    else:
        if R > 100000: # beyond distance z (e.g. 100 km), assume spreading like 2*pi*R
            E_if_station_were_at_z = 2 * pi * (z ** 2) / (rho_atmos * c_atmos) * val
            Eacoustic = E_if_station_were_at_z* R/1e5
        else:
            Eacoustic = 2 * pi * R ** 2 / (rho_atmos * c_atmos) * val
    return Eacoustic
    
def VASR(Eacoustic, Eseismic):
    # From Johnson and Aster 2005
    eta = Eacoustic / Eseismic
    return eta

def attenuation(tr, R, Q=50, c_earth=2500):
    s = tr.stats
    if 'spectrum' in s: 
        peakF = s['spectrum']['peakF']
        exponent = - ((pi) * peakF * R) / (c_earth * Q)
        A = np.exp(exponent)
        return A
    else:
        return None

def Eseismic2magnitude(Eseismic, correction=3.7):
    # after equation 7 in Hanks and Kanamori 1979, where moment is substitute with energy
    # energy in Joules rather than ergs, so correction is 3.7 rather than 10.7
    if isinstance(Eseismic, list): # list of stationEnergy
        mag = [] 
        for thisE in Eseismic:
            mag.append(Eseismic2magnitude(thisE, correction=correction))
    else:
        mag = np.log10(Eseismic)/1.5 - correction
    return mag

def magnitude2Eseismic(mag, correction=3.7):
    # after equation 7 in Hanks and Kanamori 1979, where moment is substitute with energy
    # energy in Joules rather than ergs, so correction is 3.7 rather than 10.7   
    if isinstance(mag, list): # list of stationEnergy
        Eseismic = [] 
        for thismag in mag:
            Eseismic.append(magnitude2Eseismic(thismag, correction=correction))
    else:
        Eseismic = np.power(10, 1.5 * mag + correction)
    return Eseismic

def Mlrichter(peakA, R, a=1.6, b=-0.15, g=0):
    """
    Compute Richter local magnitude (ML) from peak amplitude and distance.

    Parameters
    ----------
    peakA : float
        Peak amplitude (in mm or nm depending on calibration).
    R : float
        Epicentral distance in km.
    a : float
        Log-distance scaling (default 1.6).
    b : float
        Offset (default -0.15).
    g : float
        Station correction (default 0).

    Returns
    -------
    ml : float
        Local magnitude.
    """
    return np.log10(peakA) + a * np.log10(R) + b + g


##########################################
# Main wrapper - a level above ampengfft #
##########################################
def ampengfftmag(stream, inventory, source_coords,
                 model='body', Q=50, c_earth=2500, correction=3.7,
                 a=1.6, b=-0.15, g=0,
                 threshold=0.707, window_length=9, polyorder=2,
                 differentiate=True, verbose=True, snr_method='std',
                 snr_split_time=None, snr_window_length=1.0,
                 snr_min=None):
    """
    Wrapper to compute amplitude, spectral metrics, SNR, and magnitudes for a Stream.

    This function:
    - Computes all amplitude and spectral metrics using `ampengfft()`.
    - Estimates signal-to-noise ratios using `estimate_snr()`.
    - Attaches coordinates using the provided inventory.
    - Computes distances, energy-based magnitude (Me) and local magnitude (Ml).
    - Returns updated stream and summary DataFrame.

    Parameters
    ----------
    stream : obspy.Stream
        Input traces (usually displacement). Use `differentiate=True` to get velocity for energy.
    inventory : obspy.Inventory
        Station metadata.
    source_coords : dict
        {'latitude': ..., 'longitude': ..., 'depth': ...}
    model : str
        'body' or 'surface' wave spreading.
    Q : float
        Attenuation quality factor.
    c_earth : float
        Wave speed (m/s).
    correction : float
        Hanks & Kanamori energy magnitude offset (3.7 for SI).
    a, b, g : float
        Local magnitude coefficients.
    threshold : float
        Bandwidth -3 dB threshold.
    window_length : int
        Savitzky-Golay smoothing window length (must be odd).
    polyorder : int
        Polynomial order for smoothing filter.
    differentiate : bool
        If True, differentiate displacement trace for energy estimate.
    verbose : bool
        Print summary progress.
    snr_method : str
        SNR method ('std', 'rms', 'max', 'spectral').
    snr_split_time : UTCDateTime or tuple
        Optional start/end for signal window in SNR.
    snr_window_length : float
        Window size in seconds for SNR estimation.
    snr_min : float or None
        Minimum SNR to keep trace. If set, removes low-SNR traces.

    Returns
    -------
    stream : obspy.Stream
        Updated traces with metrics.
    df : pandas.DataFrame
        Magnitude summary table (one row per trace + optional network row).
    """

    if verbose:
        print("[2] Estimating SNR values...")
    for tr in stream:
        if not hasattr(tr.stats, 'metrics'):
            tr.stats['metrics'] = {}
        if not hasattr(tr.stats.metrics, 'snr'):
            try:
                snr, signal_val, noise_val = estimate_snr(tr,
                                                        method=snr_method,
                                                        window_length=snr_window_length,
                                                        split_time=snr_split_time,
                                                        verbose=False)
                tr.stats.metrics['snr'] = snr
                tr.stats.metrics['signal'] = signal_val
                tr.stats.metrics['noise'] = noise_val
            except Exception as e:
                if verbose:
                    print(f"[{tr.id}] SNR estimation failed: {e}")

    # filter stream based on minimum snr
    if snr_min is not None:
        if verbose:
            print(f"[3] Filtering traces with SNR < {snr_min:.1f}...")
        stream = Stream([tr for tr in stream if tr.stats.metrics.get('snr', 0) >= snr_min])

    # ampengfft
    if verbose:
        print("[1] Computing amplitude and spectral metrics...")
    ampengfft(stream, threshold=threshold, window_length=window_length,
              polyorder=polyorder, differentiate=differentiate)


    # station magnitudes
    if verbose:
        print("[4] Estimating station magnitudes...")
    compute_station_magnitudes(stream, inventory, source_coords,
                               model=model, Q=Q, c_earth=c_earth, correction=correction,
                               a=a, b=b, g=g,
                               attach_coords=True, compute_distances=True)

    # network magnitudes
    if verbose:
        print("[5] Summarizing network magnitude statistics...")
    df = summarize_magnitudes(stream, include_network=True)

    if verbose and not df.empty:
        net_ml = df.iloc[-1].get('network_mean_ML', np.nan)
        net_me = df.iloc[-1].get('network_mean_ME', np.nan)
        print(f"    → Network ML: {net_ml:.2f} | Network ME: {net_me:.2f}")

    return stream, df

from obspy.core.event import Event, Origin, Magnitude, StationMagnitude, StationMagnitudeContribution, Comment
from obspy.core.event.resourceid import ResourceIdentifier
from obspy.core.event.base import WaveformStreamID
from obspy.core import UTCDateTime
import json
import pandas as pd

def stream_to_event(stream, source_coords, df=None,
                    event_id=None, creation_time=None,
                    event_type=None, mainclass=None, subclass=None):
    """
    Convert a Stream and source metadata into an ObsPy Event object.

    Parameters
    ----------
    stream : obspy.Stream
        Stream of traces with metrics in .stats.metrics
    source_coords : dict
        {'latitude': ..., 'longitude': ..., 'depth': ...} in meters
    df : pandas.DataFrame, optional
        Output from summarize_magnitudes()
    event_id : str, optional
        Resource ID for the Event
    creation_time : str or UTCDateTime, optional
        Event creation time
    event_type : str, optional
        QuakeML-compatible event type (e.g. "earthquake", "volcanic eruption")
    mainclass : str, optional
        High-level classification (e.g., "volcano-seismic", "tectonic")
    subclass : str, optional
        Detailed classification (e.g., "hybrid", "long-period")

    Returns
    -------
    event : obspy.core.event.Event
    """
    event = Event()

    # Resource ID and creation time
    if event_id:
        event.resource_id = ResourceIdentifier(id=event_id)
    if creation_time:
        event.creation_info = {'creation_time': UTCDateTime(creation_time)}

    # Event type
    if event_type:
        event.event_type = event_type

    # Add origin
    origin = Origin(
        latitude=source_coords['latitude'],
        longitude=source_coords['longitude'],
        depth=source_coords.get('depth', 0)
    )
    event.origins.append(origin)

    # Add network-averaged magnitudes from DataFrame
    if df is not None and hasattr(df, 'iloc'):
        net_ml = df.iloc[-1].get('network_mean_ML')
        net_me = df.iloc[-1].get('network_mean_ME')

        if pd.notnull(net_ml):
            mag_ml = Magnitude(
                mag=net_ml,
                magnitude_type="ML",
                origin_id=origin.resource_id
            )
            event.magnitudes.append(mag_ml)

        if pd.notnull(net_me):
            mag_me = Magnitude(
                mag=net_me,
                magnitude_type="Me",
                origin_id=origin.resource_id
            )
            event.magnitudes.append(mag_me)

    # Add station magnitudes
    for tr in stream:
        ml = tr.stats.metrics.get('local_magnitude')
        if ml is None:
            continue

        wid = WaveformStreamID(
            network_code=tr.stats.network,
            station_code=tr.stats.station,
            location_code=tr.stats.location,
            channel_code=tr.stats.channel
        )

        smag = StationMagnitude(
            mag=ml,
            magnitude_type="ML",
            station_magnitude_type="ML",
            waveform_id=wid
        )
        event.station_magnitudes.append(smag)

        if event.magnitudes:
            contrib = StationMagnitudeContribution(
                station_magnitude_id=smag.resource_id,
                weight=1.0
            )
            event.magnitudes[0].station_magnitude_contributions.append(contrib)

    # Add trace-level metrics as JSON comment
    try:
        metrics_dict = {tr.id: tr.stats.metrics for tr in stream if hasattr(tr.stats, 'metrics')}
        comment = Comment(
            text=json.dumps(metrics_dict, indent=2),
            force_resource_id=False
        )
        event.comments.append(comment)
    except Exception:
        pass

    # Add mainclass and subclass if given
    if mainclass or subclass:
        class_comment = {
            "mainclass": mainclass,
            "subclass": subclass
        }
        comment = Comment(
            text=json.dumps(class_comment, indent=2),
            force_resource_id=False
        )
        event.comments.append(comment)

    return event




################################################################
'''
def max_3c(st):
    """ max of a 3-component seismogram """
    N = len(st)/3
    m = []

    if N.is_integer():
        st.detrend()
        for c in range(int(N)):
            y1 = st[c*3+0].data
            y2 = st[c*3+1].data
            y3 = st[c*3+2].data
            y = np.sqrt(np.square(y1) + np.square(y2) + np.square(y3))
            m.append(max(y))
    return m 
'''
def peak_amplitudes(st):   
    """ Peak Ground Motion. Should rename it peakGroundMotion """
    
    seismic1d_list = []
    seismic3d_list = []
    infrasound_list = []
    
    #ls.clean_stream(st, taper_fraction=0.05, freqmin=0.05, causal=True)
               
    # velocity, displacement, acceleration seismograms
    stV = st.select(channel='[ESBH]H?') 
    stD = stV.copy().integrate()
    for tr in stD:
        add_to_trace_history(tr, 'integrated')    
    stA = stV.copy().differentiate()
    for tr in stA:
        add_to_trace_history(tr, 'differentiated') 
     
    # Seismic vector data  
    stZ = stV.select(channel="[ESBH]HZ")
    for tr in stZ:
        thisID = tr.id[:-1]
        st3cV = stV.select(id = '%s[ENZ12RT]' % thisID)
        if len(st3cV)==3:
            st3cD = stD.select(id = '%s[ENZ12RT]' % thisID)
            st3cA = stA.select(id = '%s[ENZ12RT]' % thisID)
            md = ls.max_3c(st3cD)
            mv = ls.max_3c(st3cV)
            ma = ls.max_3c(st3cA)
            d = {'traceID':thisID, 'PGD':md[0], 'PGV':mv[0], 'PGA':ma[0], 'calib':tr.stats.calib, 'units':tr.stats.units}
            seismic3d_list.append(d)              
    seismic3d = pd.DataFrame(seismic3d_list)
    
    # Seismic 1-c data
    peakseismo1cfile = os.path.join(eventdir, 'summary_seismic_1c.csv')
    for c in range(len(stV)):
        md = max(abs(stD[c].data))
        mv = max(abs(stV[c].data))
        ma = max(abs(stA[c].data))  
        d = {'traceID':stV[c].id, 'PGD':md[0], 'PGV':mv[0], 'PGA':ma[0], 'calib':stV[c].stats.calib, 'units':stV[c].stats.units}
        seismic1d_list.append(d)    
    seismic1d = pd.DataFrame(seismic1d_list)        
            
    # Infrasound data
    peakinfrafile = os.path.join(eventdir, 'summary_infrasound.csv')
    stP = st.select(channel="[ESBH]D?")
    stPF = stP.copy().filter('bandpass', freqmin=1.0, freqmax=20.0, corners=2, zerophase=True)    
    for c in range(len(stP)):
        mb = max(abs(stP[c].data))
        mn = max(abs(stPF[c].data)) 
        d = {'traceID':stP[c].id, 'PP':mb[0], 'PPF':mn[0], 'calib':stP[c].stats.calib, 'units':stP[c].stats.units}
        infrasound_list.append(d)  
    infrasound = pd.DataFrame(infrasound_list)    
    
    return (seismic3d, seismic1d, infrasound)

from collections import defaultdict

def max_3c(st):
    """Compute maximum 3-component vector amplitude per station."""
    grouped = defaultdict(list)
    for tr in st:
        key = tr.id[:-1]  # Strip component letter
        grouped[key].append(tr)

    max_vals = []
    for traces in grouped.values():
        if len(traces) == 3:
            traces.sort(key=lambda t: t.stats.channel)  # Optional: sort by N/E/Z
            y1, y2, y3 = traces[0].data, traces[1].data, traces[2].data
            y = np.sqrt(y1**2 + y2**2 + y3**2)
            max_vals.append(np.max(y))

    return max_vals

###########################################################
#####  Helper functions for machine learning workflow #####
def choose_best_traces(st, MAX_TRACES=8, include_seismic=True, include_infrasound=False, include_uncorrected=False):

    priority = np.array([float(tr.stats.quality_factor) for tr in st])      
    for i, tr in enumerate(st):           
        if tr.stats.channel[1]=='H':
            if include_seismic:
                if tr.stats.channel[2] == 'Z':
                    priority[i] *= 2
            else:
                priority[i] = 0
        if tr.stats.channel[1]=='D':
            if include_infrasound:
                priority[i] *= 2 
            else:
                priority[i] = 0
        if not include_uncorrected:
            if 'units' in tr.stats:
                if tr.stats.units == 'Counts':
                    priority[i] = 0
            else:
                priority[i] = 0

    n = np.count_nonzero(priority > 0.0)
    n = min([n, MAX_TRACES])
    j = np.argsort(priority)
    chosen = j[-n:]  
    return chosen        
        
def select_by_index_list(st, chosen):
    st2 = Stream()
    for i, tr in enumerate(st):
        if i in chosen:
            st2.append(tr)
    return st2 

'''
def compute_stationEnergy(val):
    # seismic: eng is sum of velocity trace (in m/s) squared, divided by samples per second
    # infrasound: eng is sum of pressure trace (in Pa) squared, divided by samples per second
    if isinstance(val, Stream):
        stationEnergy =[]
        for tr in val:
            stationEnergy.append(compute_stationEnergy(tr))
    if isinstance(val, Trace):
        tr2 = val.copy()
        tr2.detrend()
        y = tr2.data
        stationEnergy = np.sum(y ** 2)*tr2.stats.delta
    return stationEnergy
'''

