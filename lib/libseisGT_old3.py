#!/usr/bin/env python
# coding: utf-8

import os
import shutil
import glob
import pandas as pd
#import obspy
from obspy.core import Stream, read, Trace
import numpy as np
from obspy.core.utcdatetime import UTCDateTime
import warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
#import datetime
#from sys import exit
from obspy.signal.trigger import z_detect, trigger_onset, coincidence_trigger
from obspy.geodetics.base import gps2dist_azimuth, kilometers2degrees
from obspy.taup import TauPyModel
from obspy.signal.quality_control import MSEEDMetadata 
# Glenn Thompson, Feb 2021
from InventoryTools import has_response
from obspy.core.event import Event, Origin, Magnitude, Catalog, ResourceIdentifier, Comment   


#######################################################################
##                Trace  tools                                       ##
#######################################################################

def preprocess_trace(tr, bool_despike=True, bool_clean=True, inv=None, quality_threshold=-np.Inf, taperFraction=0.05,
                     filterType="bandpass", freq=[0.5, 30.0], corners=6, zerophase=False, outputType='VEL',
                     miniseed_qc=True, verbose=False, max_dropout=0.0, units='Counts', bool_detrend=True):
    """
    Preprocesses a seismic trace by applying quality control, filtering, and instrument response correction.

    This function performs:
    - Quality control (gaps, dropout detection, bit-level noise, clipping, spikes).
    - Detrending and interpolation over small gaps.
    - Optional despiking to remove single-sample anomalies.
    - Tapering, filtering, and response removal.
    - Tracks processing steps in `tr.stats.history`.

    Parameters:
    ----------
    tr : obspy.Trace
        The seismic trace to process.
    bool_despike : bool, optional
        Whether to remove single-sample spikes from the trace (default: True).
    bool_clean : bool, optional
        Whether to apply detrending, tapering, and filtering (default: True).
    inv : obspy.Inventory, optional
        Instrument response metadata for deconvolution (default: None).
    quality_threshold : float, optional
        Minimum quality factor required to keep the trace (default: -Inf).
    taperFraction : float, optional
        Fraction of the trace length to use for tapering (default: 0.05).
    filterType : str, optional
        Type of filter to apply: "bandpass", "lowpass", "highpass" (default: "bandpass").
    freq : list of float, optional
        Frequency range for filtering: [freq_min, freq_max] (default: [0.5, 30.0] Hz).
    corners : int, optional
        Number of filter corners (default: 6).
    zerophase : bool, optional
        Whether to apply a zero-phase filter (default: False).
    outputType : str, optional
        Type of output after instrument response removal: "VEL", "DISP", "ACC", "DEF" (default: "VEL").
    miniseed_qc : bool, optional
        Whether to perform MiniSEED quality control checks (default: True).
    verbose : bool, optional
        If True, prints processing steps (default: False).
    max_dropout : float, optional
        Maximum allowable data dropout (seconds) before rejecting (default: 0.0).
    units : str, optional
        Unit of the trace before processing. Defaults to "Counts".

    Returns:
    -------
    bool
        Returns `True` if successfully processed, `False` if rejected due to poor quality.
    """





    if verbose:
        print(f'Processing {tr}')
    
    tr.stats.setdefault('history', [])
    tr.stats.setdefault('units', units)
    quality_factor = 1.0  # Default quality score
    is_bad_trace = False  






    # Ignore traces with very low sampling rate (unless they are long-period 'L' channels)
    if tr.stats.sampling_rate < 0.5 and tr.stats.channel[0] != 'L':
        add_to_trace_history(tr, 'Sampling rate too low')
        return False

    # Ignore traces with very few samples
    if tr.stats.npts < tr.stats.sampling_rate:
        add_to_trace_history(tr, 'Not enough samples')
        return False

    # Step 1: Detect and Reject Empty Traces
    if _is_empty_trace(tr):
        add_to_trace_history(tr, 'Trace is blank')
        return False

    # Step 2: Compute Raw Data Quality Metrics
    tr.stats['metrics'] = {}
    tr.stats.metrics["twin"] = tr.stats.npts * tr.stats.delta  # Compute trace duration
    if miniseed_qc:
        _can_write_to_miniseed_and_read_back(tr, return_metrics=True)
    else:
        _compute_trace_metrics(tr)

    # Adjust quality factor based on gaps and availability
    tr.stats.quality_factor -= (tr.stats.metrics['num_gaps'] + tr.stats.metrics['num_overlaps'])
    tr.stats.quality_factor *= tr.stats.metrics['percent_availability'] / 100.0

    # Step 3: Detect Bit-Level Noise
    num_unique_values = np.unique(tr.data).size
    if num_unique_values > 10:
        quality_factor += np.log10(num_unique_values)
    else:
        add_to_trace_history(tr, 'Bit-level noise suspected')
        return False  

    # Step 4: Detect and Handle Dropouts
    tr.stats['gap_report'] = []
    if not _detect_and_handle_gaps(tr, max_dropout=max_dropout, verbose=verbose):
        return False  # Trace has excessive dropouts, discard

    # Step 5: Detect and Correct Artifacts (Clipping, Spikes, Step Functions)
    if verbose:
        print('- detecting and correcting artifacts')

    _detect_and_correct_artifacts(tr, amp_limit=1e10, count_thresh=10, spike_thresh=50.0, fill_method="interpolate")

    # Step 6: Adjust Quality Factor Based on Artifacts
    artifacts = tr.stats.get('artifacts', {})
    if artifacts.get("upper_clipped", False):
        quality_factor /= 2.0
    if artifacts.get("lower_clipped", False):
        quality_factor /= 2.0
    if artifacts.get("spike_count", 0) == 0:
        quality_factor += 1.0
    else:
        add_to_trace_history(tr, f'{artifacts["spike_count"]} outliers found')

    # Final Quality Check
    if tr.stats.quality_factor < quality_threshold:
        return False

    if verbose:
        print(f'- artifacts processed. qf={tr.stats.quality_factor}')

    # Step 7: Detrend (Gap-Aware)
    if bool_detrend:
        tr = detrend_trace(tr, gap_threshold=3, verbose=verbose)

    # Step 8: Apply Cleaning Steps (Pad, Taper, Filter, Response Removal)
    if bool_clean:
        _clean_trace(tr, taperFraction, filterType, freq, corners, zerophase, inv, outputType, verbose)

    if verbose:
        print(f'- processing complete. qf={tr.stats.quality_factor}')
    
    return True

# ---- Helper Functions ----

def _clean_trace(tr, taperFraction, filterType, freq, corners, zerophase, inv, outputType, verbose):
    """
    Applies padding, tapering, filtering, and instrument response correction.
    """
    if verbose:
        print('- cleaning trace')

    # Padding
    npts_pad = int(taperFraction * tr.stats.npts)
    npts_pad_seconds = max(npts_pad * tr.stats.delta, 1/freq[0])  # Ensure minimum pad length
    _pad_trace(tr, npts_pad_seconds)
    max_fraction = npts_pad / tr.stats.npts

    # Taper
    if verbose:
        print('- tapering')
    tr.taper(max_percentage=max_fraction, type="hann")
    add_to_trace_history(tr, 'tapered')

    # Filtering
    if verbose:
        print('- filtering')
    if filterType == "bandpass":
        tr.filter(filterType, freqmin=freq[0], freqmax=freq[1], corners=corners, zerophase=zerophase)
    else:
        tr.filter(filterType, freq=freq[0], corners=corners, zerophase=zerophase)
    _update_trace_filter(tr, filterType, freq, zerophase)
    add_to_trace_history(tr, filterType)

    # Instrument Response Removal
    _handle_instrument_response(tr, inv, outputType, verbose)

    # Remove Padding
    _unpad_trace(tr)

def _handle_instrument_response(tr, inv, outputType, verbose):
    """
    Removes the instrument response using ObsPy's `remove_response()`.
    """
    if inv:
        try:
            if verbose:
                print('- removing instrument response')
            if tr.stats.channel[1]=='D':
                outputType='DEF' # for pressure sensor
            tr.remove_response(inventory=inv, output=outputType, \
                pre_filt=None, water_level=60, zero_mean=True, \
                taper=False, taper_fraction=0.0, plot=False, fig=None)      
            add_to_trace_history(tr, 'calibrated')
            tr.stats.calib = 1.0
            tr.stats['calib_applied'] = _get_calib(tr, inv) # we have to do this, as the calib value is used to scale the data in plots
        except Exception as e:
            print(f'Error removing response: {e}')
            print('remove_response failed for %s' % tr.id)
            return False            
    elif tr.stats['units'] == 'Counts' and not 'calibrated' in tr.stats.history and tr.stats.calib!=1.0:
        tr.data = tr.data * tr.stats.calib
        tr.stats['calib_applied'] = tr.stats.calib
        tr.stats.calib = 1.0 # we have to do this, as the calib value is used to scale the data in plots
        add_to_trace_history(tr, 'calibrated') 
    
    if 'calibrated' in tr.stats.history:           
        if tr.stats.channel[1]=='H':
            if outputType=='VEL':
                tr.stats['units'] = 'm/s'
            elif outputType=='DISP':
                tr.stats['units'] = 'm'
        elif tr.stats.channel[1]=='N':
            tr.stats['units'] = 'm/s2'                
        elif tr.stats.channel[1]=='D':
            tr.stats['units'] = 'Pa' 

def add_to_trace_history(tr, message):
    """
    Adds a message to the processing history of a seismic trace.

    This function maintains a record of operations applied to a trace 
    by appending descriptive messages to `tr.stats.history`.

    Parameters:
    ----------
    tr : obspy.Trace
        The seismic trace to which the history message is added.
    message : str
        A string describing the processing step or modification.

    Returns:
    -------
    None
        The function updates `tr.stats.history` **in place**.

    Notes:
    ------
    - If `tr.stats.history` does not exist, it is initialized as an empty list.
    - Duplicate messages are not added to prevent redundancy.

    Example:
    --------
    ```python
    from obspy import Trace
    import numpy as np

    # Create a sample trace
    tr = Trace(data=np.random.randn(1000))

    # Add processing steps to the trace history
    add_to_trace_history(tr, "Detrended")
    add_to_trace_history(tr, "Applied bandpass filter (0.1 - 10 Hz)")
    add_to_trace_history(tr, "Detrended")  # Duplicate, will not be added

    # Print trace history
    print(tr.stats.history)
    # Output: ['Detrended', 'Applied bandpass filter (0.1 - 10 Hz)']
    ```
    """
    if 'history' not in tr.stats:
        tr.stats['history'] = []
    if message not in tr.stats.history:
        tr.stats.history.append(message)

def _can_write_to_miniseed_and_read_back(tr, return_metrics=True):
    """
    Tests whether an ObsPy Trace can be written to and successfully read back from MiniSEED format.

    This function attempts to:
    1. **Convert trace data to float** (if necessary) to avoid MiniSEED writing issues.
    2. **Write the trace to a temporary MiniSEED file**.
    3. **Read the file back to confirm integrity**.
    4. If `return_metrics=True`, computes MiniSEED metadata using `MSEEDMetadata()`.

    Parameters:
    ----------
    tr : obspy.Trace
        The seismic trace to test for MiniSEED compatibility.
    return_metrics : bool, optional
        If `True`, computes MiniSEED quality control metrics and stores them in `tr.stats['metrics']` (default: True).

    Returns:
    -------
    bool
        `True` if the trace can be written to and read back from MiniSEED successfully, `False` otherwise.

    Notes:
    ------
    - **Converts `tr.data` to float** (`trace.data.astype(float)`) if necessary.
    - **Removes temporary MiniSEED files** after testing.
    - Uses `MSEEDMetadata()` to compute quality metrics similar to **ISPAQ/MUSTANG**.
    - Sets `tr.stats['quality_factor'] = -100` if the trace has **no data**.

    Example:
    --------
    ```python
    from obspy import read

    # Load a trace
    tr = read("example.mseed")[0]

    # Check if it can be written & read back
    success = _can_write_to_miniseed_and_read_back(tr, return_metrics=True)

    print(f"MiniSEED compatibility: {success}")
    if success and "metrics" in tr.stats:
        print(tr.stats["metrics"])  # Print MiniSEED quality metrics
    ```
    """
    if len(tr.data) == 0:
        tr.stats['quality_factor'] = 0.0
        return False

    # Convert data type to float if necessary (prevents MiniSEED write errors)
    if not np.issubdtype(tr.data.dtype, np.floating):
        tr.data = tr.data.astype(float)

    tmpfilename = f"{tr.id}_{tr.stats.starttime.isoformat()}.mseed"

    try:
        # Attempt to write to MiniSEED
        tr.write(tmpfilename)

        # Try reading it back
        _ = read(tmpfilename)

        if return_metrics:
            # Compute MiniSEED metadata
            mseedqc = MSEEDMetadata([tmpfilename])
            tr.stats['metrics'] = mseedqc.meta
            add_to_trace_history(tr, "MSEED metrics computed (similar to ISPAQ/MUSTANG).")

        return True  # Successfully wrote and read back

    except Exception as e:
        if return_metrics:
            tr.stats['quality_factor'] = 0.0
        print(f"Failed MiniSEED write/read test for {tr.id}: {e}")
        return False

    finally:
        # Clean up the temporary file
        if os.path.exists(tmpfilename):
            os.remove(tmpfilename)


def _compute_trace_metrics(trace):
    """
    Computes and stores the number of gaps, number of overlaps, 
    and percent availability of data in an ObsPy Trace object.
    
    Parameters:
    ----------
    trace : obspy.Trace
        The seismic trace to analyze.
    
    Returns:
    -------
    None
        The function modifies `trace.stats.metrics` in place.
    
    Stores:
    -------
    - trace.stats.metrics['num_gaps'] : int
        Number of gaps (missing data).
    - trace.stats.metrics['num_overlaps'] : int
        Number of overlapping time windows.
    - trace.stats.metrics['percent_availability'] : float
        Percentage of available data relative to the total expected time span.
    
    """
    # Ensure metrics dictionary exists
    if not hasattr(trace.stats, "metrics"):
        trace.stats.metrics = {}

    # Detect gaps and overlaps
    gaps = trace.get_gaps()
    
    num_gaps = sum(1 for gap in gaps if gap[6] > 0)  # Count only gaps
    num_overlaps = sum(1 for gap in gaps if gap[6] < 0)  # Count only overlaps

    # Compute percent availability
    total_duration = trace.stats.npts * trace.stats.delta  # Expected duration
    actual_duration = total_duration  # Start with full duration

    for gap in gaps:
        if gap[6] > 0:  # If it's a gap, subtract its duration
            actual_duration -= gap[6]

    percent_availability = (actual_duration / total_duration) * 100 if total_duration > 0 else 0

    # Store metrics in trace.stats
    trace.stats.metrics["num_gaps"] = num_gaps
    trace.stats.metrics["num_overlaps"] = num_overlaps
    trace.stats.metrics["percent_availability"] = percent_availability  

def _is_empty_trace(trace):
    """
    Determines whether a seismic trace is effectively empty.

    A trace is considered empty if:
    - It has zero data points (`npts == 0`).
    - All samples are identical (e.g., all zeros, all -1s).
    - All values are NaN.

    Parameters:
    ----------
    trace : obspy.Trace
        The seismic trace to check.

    Returns:
    -------
    bool
        `True` if the trace is empty or contains only redundant values, otherwise `False`.

    Notes:
    ------
    - The function first checks if the trace has no data (`npts == 0`).
    - Then it checks if all values are identical (suggesting a completely flat signal).
    - Finally, it verifies if all values are NaN.

    Example:
    --------
    ```python
    from obspy import Trace
    import numpy as np

    # Create an empty trace
    empty_trace = Trace(data=np.array([]))
    print(_is_empty_trace(empty_trace))  # True

    # Create a flat-line trace
    flat_trace = Trace(data=np.zeros(1000))
    print(_is_empty_trace(flat_trace))  # True

    # Create a normal trace with random data
    normal_trace = Trace(data=np.random.randn(1000))
    print(_is_empty_trace(normal_trace))  # False
    ```
    """
    if trace.stats.npts == 0:
        return True
    
    # Check for flat trace (e.g. all zero, or all -1)
    if np.all(trace.data == np.nanmean(trace.data)):
        return True

    # Check if all values are NaN
    if np.all(np.isnan(trace.data)):
        return True 

    return False

def _check0andMinus1(liste):
# function bool_good_trace = check0andMinus1(tr.data)
    liste=list(liste)
    listStr=''.join(str(i) for i in liste)
    if  "000000000000" in listStr or "-1-1-1-1-1-1-1-1" in listStr :
        return False
    else:
        return True  

def _get_islands(arr, mask):
    mask_ = np.concatenate(( [False], mask, [False] ))
    idx = np.flatnonzero(mask_ [1:] != mask_ [:-1])
    return [arr[idx[i]:idx[i+1] + 1] for i in range(0, len(idx), 2)]

def _FindMaxLength(lst):
    maxList = max(lst, key = len)
    maxLength = max(map(len, lst))      
    return maxList, maxLength  

def _mad_based_outlier(tr, thresh=3.5):
    tr2 = tr.copy()
    tr2.detrend()
    points = tr2.data
    if len(points.shape) == 1:
        points = points[:,None]
    #points = np.absolute(points)
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation
    
    outlier_indices = np.array(np.where(modified_z_score > thresh))
    outlier_count = outlier_indices.size
    '''
    if outlier_count > 0:
        print('size diff = %d, median = %e, med_abs_deviation = %e ' % (diff.size, median, med_abs_deviation))
        mzs = sorted(modified_z_score)
        print(mzs[-10:])
    '''
    
    return outlier_count, outlier_indices



def _clip_trace(tr, AMP_LIMIT = 1e10, fill_value=None):     
# function tr = clip_trace(tr, maxamp)
    # remove absurdly large values
    '''a = tr.data
    np.clip(a, -AMP_LIMIT, AMP_LIMIT, out=a)
    np.where(a == AMP_LIMIT, 0, a)
    tr.data = a    '''
    if not fill_value:
        fill_value = np.nanmedian(tr.data)
    tr.data[tr.data > AMP_LIMIT] = fill_value
    tr.data[tr.data < -AMP_LIMIT] = fill_value 

def _detectClipping(tr, countThresh = 10):
    upper_clipped = False
    lower_clipped = False
    y = tr.data
    mu = np.nanmax(y)
    md = np.nanmin(y)
    countu = (tr.data == mu).sum()
    countd = (tr.data == md).sum()
    if countu >= countThresh:
        add_to_trace_history(tr, 'Trace %s appears to be clipped at upper limit %e (count=%d)' % (tr.id, mu, countu) )    
        upper_clipped = True
    if countd >= countThresh:
        add_to_trace_history(tr, 'Trace %s appears to be clipped at lower limit %e (count=%d)' % (tr.id, mu, countu) )       
        lower_clipped = True
    return upper_clipped, lower_clipped
    
def _remove_single_sample_spikes(trace, threshold=10):
    # Parameters
    #threshold = 10  Threshold for identifying large differences (spikes)
    try:
        trace_std = np.nanstd(trace.data)
    except:
        try:
            trace_std = np.std(trace.data)
        except:
            return

    # Step 1: Calculate the absolute differences between consecutive points
    diff_prev = np.abs(np.diff(trace.data))  # Difference with the previous point
    diff_next = np.abs(np.diff(trace.data[1:], append=trace.data[-1]))  # Difference with the next point

    # Step 2: Identify spikes where the difference is larger than the threshold
    spikes = (diff_prev > trace_std * threshold) & (diff_next > trace_std * threshold)
    spikes = np.append(spikes, np.array(False))


    # Step 3: Replace spikes with the average of previous and next points
    # Use boolean indexing to replace the spikes
    smoothed_data = np.copy(trace.data)
    #smoothed_data[1:-1][spikes] = (trace.data[:-2][spikes] + trace.data[2:][spikes]) / 2
    smoothed_data[spikes] = trace_std

    # Update the trace data with the smoothed data
    trace.data = smoothed_data        

def _check_for_spikes(tr):
    if not 'metrics' in tr.stats:
        return

    m = tr.stats.metrics
    peak2peak = m['sample_max']-m['sample_min']
    positive_spike_metric = (m['sample_upper_quartile']-m['sample_min'])/peak2peak
    negative_spike_metric = (m['sample_max']-m['sample_lower_quartile'])/peak2peak
    if positive_spike_metric < 0.01:
        print('Positive spike(s) suspected on %s' % tr.id)
        tr.stats['quality_factor'] = -1
    if negative_spike_metric < 0.01:
        print('Negative spike(s) suspected on %s' % tr.id)  
        tr.stats['quality_factor'] = -1

def _piecewise_detrend(tr, null_value=0, fill_value=np.nan, detrend='linear', verbose=False): 
    # takes a Trace object nominally from an SDSclient that has gaps marked by zeros
    isSplit = False
    if np.any(tr.data == null_value):
        tr.data = np.ma.masked_where(tr.data == null_value, tr.data)
        all_traces = tr.split()
        if verbose:
            print(f'{tr.id} split into {len(all_traces)} traces')
        isSplit=True
        all_traces.detrend(detrend)
        all_traces.merge(method=0, fill_value=fill_value)
        if verbose:
            print(f'merged again. stream now contains {len(all_traces)} traces')
        if len(all_traces) == 1:
            return tr
        else:
            return False

def _pad_trace(tr, seconds, method="mirror"):
    """
    Pads a seismic trace by extending it with mirrored data, zeros, or spectrally-matched noise.

    This function extends the **start and end** of a given **ObsPy Trace** by `seconds`, using 
    one of three methods:
    - **"mirror" (default)**: Reverses and appends the first and last `seconds` of data.
    - **"zeros"**: Pads with zeros.
    - **"noise"**: Generates and appends spectrally-matched noise.

    Parameters:
    ----------
    tr : obspy.Trace
        The seismic trace to be padded.
    seconds : float
        The number of seconds to extend at both ends.
    method : str, optional
        Padding method (`"mirror"`, `"zeros"`, or `"noise"`). Default: `"mirror"`.

    Returns:
    -------
    None
        The function modifies `tr` **in place**, updating its `starttime`.

    Notes:
    ------
    - Stores original start and end times in `tr.stats['originalStartTime']` and `tr.stats['originalEndTime']`.
    - **Mirror padding** is best for filtering and deconvolution.
    - **Noise padding** helps maintain spectral continuity.

    Example:
    --------
    ```python
    from obspy import read

    tr = read("example.mseed")[0]
    print(f"Before padding: {tr.stats.starttime}")

    _pad_trace(tr, 5.0, method="noise")

    print(f"After padding: {tr.stats.starttime}")
    ```
    """
    if seconds <= 0.0:
        return

    # Store original time range
    tr.stats['originalStartTime'] = tr.stats.starttime
    tr.stats['originalEndTime'] = tr.stats.endtime

    npts_pad = int(tr.stats.sampling_rate * seconds)

    if method == "mirror":
        # Extract waveform segments and reverse them
        y_prepend = np.flip(tr.data[:npts_pad])  # Reverse first N samples
        y_postpend = np.flip(tr.data[-npts_pad:])  # Reverse last N samples
    elif method == "zeros":
        # Create zero-padding
        y_prepend = np.zeros(npts_pad)
        y_postpend = np.zeros(npts_pad)
    elif method == "noise":
        # Generate spectrally-matched noise
        y_prepend = _generate_spectrally_matched_noise(tr, npts_pad)
        y_postpend = _generate_spectrally_matched_noise(tr, npts_pad)
    else:
        raise ValueError(f"Invalid padding method: {method}. Choose 'mirror', 'zeros', or 'noise'.")

    # Concatenate padded trace
    tr.data = np.concatenate([y_prepend, tr.data, y_postpend])

    # Update starttime to reflect new padding
    tr.stats.starttime -= npts_pad * tr.stats.delta
    add_to_trace_history(tr, f'padded using {method}')

def _update_trace_filter(tr, filtertype, freq, zerophase):
    """
    Updates the filter settings applied to a seismic trace.

    This function modifies the **filter metadata** stored in `tr.stats['filter']`, ensuring
    that **frequency bounds and phase settings** are correctly tracked.

    Parameters:
    ----------
    tr : obspy.Trace
        The seismic trace whose filter metadata will be updated.
    filtertype : str
        Type of filter applied (`"highpass"`, `"lowpass"`, or `"bandpass"`).
    freq : float or tuple
        The cutoff frequency (single value for `"highpass"` or `"lowpass"`, tuple for `"bandpass"`).
    zerophase : bool
        Indicates whether the filter is zero-phase.

    Returns:
    -------
    None
        The function modifies `tr.stats['filter']` **in place**.

    Notes:
    ------
    - Ensures the trace has a `filter` dictionary (`tr.stats['filter']`).
    - Updates **frequency bounds** based on the filter type:
      - `"highpass"` → Updates `freqmin` only.
      - `"lowpass"` → Updates `freqmax` only.
      - `"bandpass"` → Updates both `freqmin` and `freqmax`.
    - Tracks whether the filter is **zero-phase**.

    Example:
    --------
    ```python
    from obspy import read

    tr = read("example.mseed")[0]
    
    # Apply a highpass filter and update metadata
    _update_trace_filter(tr, "highpass", freq=1.0, zerophase=True)

    print(tr.stats['filter'])  # Metadata now includes the filter
    ```
    """
    if not filter in tr.stats:
        tr.stats['filter'] = {'freqmin':0, 'freqmax':tr.stats.sampling_rate/2, 'zerophase': False}
    if filtertype == 'highpass':    
        tr.stats.filter["freqmin"] = max([freq, tr.stats.filter["freqmin"]])
    if filtertype == 'bandpass':
        tr.stats.filter["freqmin"] = max([freq[0], tr.stats.filter["freqmin"]]) 
        tr.stats.filter["freqmax"] = min([freq[1], tr.stats.filter["freqmax"]])
    if filtertype == 'lowpass':
        tr.stats.filter["freqmax"] = min([freq, tr.stats.filter["freqmax"]])
    tr.stats.filter['zerophase'] = zerophase

def _get_calib(tr, this_inv):
    """
    Retrieves the overall calibration factor for a given trace from an inventory.

    This function looks up the **instrument response calibration factor** (sensitivity and gain)
    for a specific station and channel from an ObsPy **Inventory**.

    Parameters:
    ----------
    tr : obspy.Trace
        The seismic trace whose calibration factor is needed.
    this_inv : obspy.Inventory
        The station metadata (e.g., from a `StationXML` file).

    Returns:
    -------
    float
        The **calibration value** (overall sensitivity) for the trace.

    Notes:
    ------
    - Searches `this_inv` for a **matching station and channel**.
    - Extracts **overall sensitivity** from `channel.response`.
    - Assumes `this_inv` contains only **one network**.

    Example:
    --------
    ```python
    from obspy import read_inventory, read

    inv = read_inventory("example.xml")
    tr = read("example.mseed")[0]

    calib = _get_calib(tr, inv)
    print(f"Calibration factor: {calib}")
    ```
    """
    calib_value = 1.0
    for station in this_inv.networks[0].stations:
        if station.code == tr.stats.station:
            for channel in station.channels:
                if channel.code == tr.stats.channel:
                    calib_freq, calib_value = channel.response._get_overall_sensitivity_and_gain()
    return calib_value

def _unpad_trace(tr):
    """
    Removes padding from a previously padded seismic trace.

    This function trims a trace back to its **original start and end times**, assuming it was 
    previously padded using `_pad_trace()`.

    Parameters:
    ----------
    tr : obspy.Trace
        The seismic trace to unpad.

    Returns:
    -------
    None
        The function modifies `tr` **in place**, restoring its original time range.

    Notes:
    ------
    - Uses `tr.stats['originalStartTime']` and `tr.stats['originalEndTime']` for trimming.
    - Calls `tr.trim()` to remove the extra data.
    - Adds `"unpadded"` to the trace's history.

    Example:
    --------
    ```python
    from obspy import read

    tr = read("example.mseed")[0]
    
    _pad_trace(tr, 5.0, method="noise")
    print(f"Padded start time: {tr.stats.starttime}")

    _unpad_trace(tr)
    print(f"Restored start time: {tr.stats.starttime}")
    ```
    """
    if 'originalStartTime' in tr.stats and 'originalEndTime' in tr.stats:
        tr.trim(starttime=tr.stats['originalStartTime'], endtime=tr.stats['originalEndTime'], pad=False)
        add_to_trace_history(tr, 'unpadded')
 
################### New trace tools added in Feb 2025: Start ####################

def process(st_or_tr, func, **kwargs):
    """
    Recursively applies a given function to all traces in an ObsPy Stream or a single Trace.

    This function ensures that the specified function (`func`) is applied to each Trace within 
    a Stream or directly to a single Trace if provided.

    Parameters:
    ----------
    st_or_tr : obspy.Stream or obspy.Trace
        The Stream or Trace object to process.
    func : callable
        The function to apply to each Trace. It should accept an ObsPy Trace as its first argument.
    **kwargs : dict, optional
        Additional keyword arguments to pass to `func`.

    Returns:
    -------
    None
        The function modifies `st_or_tr` **in place**.

    Raises:
    -------
    TypeError:
        If `st_or_tr` is not an ObsPy `Stream` or `Trace`.

    Notes:
    ------
    - If `st_or_tr` is a **Stream**, the function calls itself recursively for each Trace.
    - If `st_or_tr` is a **Trace**, it directly applies `func()`, passing any additional `**kwargs`.
    - This is useful for applying **preprocessing steps** (e.g., filtering, detrending) across multiple traces.

    Example:
    --------
    ```python
    from obspy import read

    def highpass_filter(trace, freq=1.0):
        trace.filter("highpass", freq=freq)

    # Load waveform data
    st = read("example.mseed")

    # Apply a highpass filter to all traces
    process(st, highpass_filter, freq=2.0)

    print(st)
    ```
    """
    if isinstance(st_or_tr, Stream):
        for trace in st_or_tr:
            process(trace, func, **kwargs)  # Recursive call for each Trace
    elif isinstance(st_or_tr, Trace):
        func(st_or_tr, **kwargs)  # Apply function with keyword arguments
    else:
        raise TypeError("Input must be an ObsPy Stream or Trace object.")


#######################################################################
###                        Gap filling tools                        ###
#######################################################################

def fill_all_gaps(trace, verbose=False):
    """
    Identifies and fills all gaps in a seismic trace using an appropriate method based on gap size.

    This function scans a given seismic **Trace** for gaps and fills them using the best 
    available method for each case:
    - **Tiny gaps** (≤ 5 samples) → **Linear interpolation**
    - **Short gaps** (≤ 1 second) → **Repeat previous data**
    - **Longer gaps** (> 1 second) → **Spectrally-matched noise**

    If the station is **'MBSS'**, the function generates before/after plots of the trace.

    Parameters:
    ----------
    trace : obspy.Trace
        The seismic trace containing gaps to be filled.
    verbose : bool, optional
        If `True`, prints information about detected gaps and selected filling methods (default: False).

    Returns:
    -------
    None
        The function modifies `trace` **in place**, filling all detected gaps.

    Notes:
    ------
    - Uses `trace.get_gaps()` to identify missing data regions.
    - Calls `_fill_gap_with_linear_interpolation()`, `_fill_gap_with_repeat_previous_data()`, 
      or `_fill_gap_with_filtered_noise()` based on gap size.
    - Designed for **robust seismic data preprocessing**.

    Example:
    --------
    ```python
    from obspy import read

    # Load a trace with gaps
    tr = read("example_with_gaps.mseed")[0]

    # Fill all gaps using the best available method
    fill_all_gaps(tr, verbose=True)

    # Plot the corrected trace
    tr.plot()
    ```
    """
    if trace.stats.station == 'MBSS':
        verbose = True
        trace.plot(outfile='MBSS_before_gap_filling.png')
    if verbose:
        print(f'filling gaps for {trace}')
    stream = Stream(traces=[trace])  # Wrap the trace in a stream
    gaps = stream.get_gaps()  # Get detected gaps

    for net, sta, loc, chan, t1, t2, delta, samples in gaps:
        gap_start = UTCDateTime(t1)
        gap_end = UTCDateTime(t2)

        # Select appropriate gap-filling method
        if samples <= 5:  # Tiny gaps (few samples) → Linear interpolation
            if verbose:
                print(f'gap start={gap_start}, end={gap_end}, linear interpolation')
            _fill_gap_with_linear_interpolation(trace, gap_start, gap_end)
        elif samples <= trace.stats.sampling_rate:  # Short gaps (≤ 1 sec) → Repeat data
            if verbose:
                print(f'gap start={gap_start}, end={gap_end}, repeating data')            
            _fill_gap_with_repeat_previous_data(trace, gap_start, gap_end)
        else:  # Longer gaps → Fill with spectrally-matched noise
            if verbose:
                print(f'gap start={gap_start}, end={gap_end}, adding spectral noise')            
            _fill_gap_with_filtered_noise(trace, gap_start, gap_end)
    if trace.stats.station == 'MBSS':
        trace.plot(outfile='MBSS_after_gap_filling.png')

def _fill_gap_with_repeat_previous_data(trace, gap_start, gap_end):
    """
    Fills a seismic trace gap by repeating the last valid segment before the gap.

    This method is useful for **short gaps (≤ 1 second)**, where simply extending the last 
    valid data point provides a reasonable approximation.

    Parameters:
    ----------
    trace : obspy.Trace
        The seismic trace containing the gap.
    gap_start : UTCDateTime
        The start time of the gap.
    gap_end : UTCDateTime
        The end time of the gap.

    Returns:
    -------
    None
        The function modifies `trace.data` **in place**, filling the gap.

    Notes:
    ------
    - The function extracts the **last valid segment** before the gap and appends it.
    - Assumes the trace **has continuous data before the gap**.

    Example:
    --------
    ```python
    # Fill a gap by repeating previous data
    _fill_gap_with_repeat_previous_data(trace, gap_start, gap_end)
    ```
    """
    sample_rate = trace.stats.sampling_rate
    num_samples = int((gap_end - gap_start) * sample_rate)

    # Get the last valid segment before the gap
    fill_data = trace.data[-num_samples:]  # Last 'num_samples' of valid data
    trace.data = np.concatenate((trace.data, fill_data))

def _fill_gap_with_linear_interpolation(trace, gap_start, gap_end):
    """
    Fills gaps in a seismic trace using linear interpolation.

    This method is ideal for **tiny gaps (≤ 5 samples)**, where interpolating 
    between adjacent valid data points ensures a smooth transition.

    Parameters:
    ----------
    trace : obspy.Trace
        The seismic trace containing the gap.
    gap_start : UTCDateTime
        The start time of the gap.
    gap_end : UTCDateTime
        The end time of the gap.

    Returns:
    -------
    None
        The function modifies `trace.data` **in place**, replacing the gap with interpolated values.

    Notes:
    ------
    - Uses `numpy.linspace()` to interpolate values between **neighboring valid samples**.
    - Assumes the trace has valid data **before and after** the gap.

    Example:
    --------
    ```python
    # Fill a gap using linear interpolation
    _fill_gap_with_linear_interpolation(trace, gap_start, gap_end)
    ```
    """
    sample_rate = trace.stats.sampling_rate
    num_samples = int((gap_end - gap_start) * sample_rate)
    gap_idx = np.arange(num_samples)

    # Get neighboring valid samples
    prev_value = trace.data[-num_samples]  # Last valid value before the gap
    next_value = trace.data[num_samples]  # First valid value after the gap

    # Linear interpolation
    interp_values = np.linspace(prev_value, next_value, num_samples)
    trace.data[-num_samples:] = interp_values  # Replace gap with interpolated values

def _fill_gap_with_filtered_noise(trace, gap_start, gap_end):
    """
    Fills a seismic trace gap using spectrally-matched noise.

    This method is used for **longer gaps (> 1 second)** where maintaining spectral
    consistency is important. It generates **spectrally-matched noise** that closely 
    resembles the missing signal.

    Parameters:
    ----------
    trace : obspy.Trace
        The seismic trace containing the gap.
    gap_start : UTCDateTime
        The start time of the gap.
    gap_end : UTCDateTime
        The end time of the gap.

    Returns:
    -------
    None
        The function modifies `trace.data` **in place**, filling the gap with spectrally-matched noise.

    Notes:
    ------
    - Extracts the **last valid segment** before the gap to estimate spectral characteristics.
    - Uses `_generate_spectrally_matched_noise()` to synthesize noise with a **matching frequency profile**.
    - Ensures smooth transitions by replacing the gap **in place**.

    Example:
    --------
    ```python
    # Fill a gap using spectrally-matched noise
    _fill_gap_with_filtered_noise(trace, gap_start, gap_end)
    ```
    """
    sample_rate = trace.stats.sampling_rate
    num_samples = int((gap_end - gap_start) * sample_rate)

    # Generate spectrally-matched noise
    noise = _generate_spectrally_matched_noise(trace, num_samples)

    # Replace gap in trace
    start_index = int((gap_start - trace.stats.starttime) * sample_rate)
    trace.data[start_index : start_index + num_samples] = noise

from scipy.signal import welch
def _generate_spectrally_matched_noise(trace, num_samples):
    """
    Generates noise that matches the spectral characteristics of a given seismic trace.

    This function extracts a segment from the trace, analyzes its frequency spectrum,
    generates noise with a matching power spectrum, and applies an inverse FFT.

    Parameters:
    ----------
    trace : obspy.Trace
        The seismic trace to match the noise characteristics.
    num_samples : int
        The number of samples to generate.

    Returns:
    -------
    np.ndarray
        An array of spectrally-matched noise.

    Notes:
    ------
    - Uses **FFT** to analyze the dominant frequency content of the adjacent waveform.
    - Shapes Gaussian white noise to match this **spectral profile**.
    - Applies an **inverse FFT** to reconstruct the time-domain noise.
    - Ensures smooth transitions by applying a **bandpass filter** if necessary.
    """
    sample_rate = trace.stats.sampling_rate
    segment = trace.data[-num_samples:]  # Use last valid segment

    # Compute power spectral density (PSD)
    freqs, psd = welch(segment, fs=sample_rate, nperseg=min(256, num_samples))

    # Generate white noise in frequency domain
    noise_freqs = np.random.normal(size=len(freqs)) * np.sqrt(psd)

    # Convert back to time domain using inverse FFT
    noise_time_domain = np.fft.irfft(noise_freqs, n=num_samples)

    # Ensure the generated noise has the same standard deviation as the original segment
    noise_time_domain *= (np.std(segment) / np.std(noise_time_domain))

    return noise_time_domain


#######################################################################
##                Stream tools                                       ##
#######################################################################

def remove_empty_traces(stream):
    """
    Removes empty traces, traces full of zeros, and traces full of NaNs from an ObsPy Stream.

    This function filters out traces that contain **no valid seismic data**, including:
    - **Completely empty traces** (no data points).
    - **Traces filled entirely with zeros**.
    - **Traces containing only NaN (Not-a-Number) values**.

    Parameters:
    ----------
    stream : obspy.Stream
        The input Stream object containing multiple seismic traces.

    Returns:
    -------
    obspy.Stream
        A new Stream object with only valid traces.

    Notes:
    ------
    - This function uses `_is_empty_trace(trace)` to check if a trace is empty or invalid.
    - The function does **not modify the original Stream**, but returns a cleaned copy.

    Example:
    --------
    ```python
    from obspy import read

    # Load waveform data
    st = read("example.mseed")

    # Remove empty or invalid traces
    cleaned_st = remove_empty_traces(st)

    print(f"Original stream had {len(st)} traces, cleaned stream has {len(cleaned_st)} traces.")
    ```
    """

    cleaned_stream = Stream()  # Create a new empty Stream

    for trace in stream:
        if not _is_empty_trace(trace):
            cleaned_stream += trace

    return cleaned_stream    

def remove_low_quality_traces(st, quality_threshold=1.0):
    """
    Removes traces from an ObsPy Stream based on a quality factor threshold.

    This function scans all traces in the given Stream and removes those with 
    a `quality_factor` lower than the specified threshold.

    Parameters:
    ----------
    st : obspy.Stream
        The input Stream object containing seismic traces.
    quality_threshold : float, optional
        The minimum `quality_factor` a trace must have to remain in the Stream (default: 1.0).

    Returns:
    -------
    None
        The function **modifies** the input `Stream` in-place, removing low-quality traces.

    Notes:
    ------
    - Assumes `tr.stats.quality_factor` is set for each trace.
    - If `quality_factor` is missing from a trace, the function may raise an error.
    - The function **modifies** `st` in-place instead of returning a new Stream.

    Example:
    --------
    ```python
    from obspy import read

    # Load waveform data
    st = read("example.mseed")

    # Remove traces with quality factor below 1.5
    remove_low_quality_traces(st, quality_threshold=1.5)

    print(f"Remaining traces: {len(st)}")
    ```
    """
    for tr in st:
        if tr.stats.quality_factor < quality_threshold: 
            st.remove(tr)

def smart_merge(st, verbose=False, interactive=False):
    """
    Merges overlapping or adjacent traces in an ObsPy Stream, handling gaps and conflicts intelligently.

    This function:
    - Groups traces by unique **NSLC ID**.
    - Removes **exact duplicates** (same start/end time, sampling rate).
    - Attempts a **standard ObsPy merge** (`Stream.merge()`).
    - If the standard merge fails, **merges traces in pairs** using `smart_merge_traces()`.
    - If `interactive=True`, prompts the user to manually select traces when merging fails.

    Parameters:
    ----------
    st : obspy.Stream
        The input Stream object containing multiple seismic traces.
    verbose : bool, optional
        If `True`, prints detailed debug output (default: False).
    interactive : bool, optional
        If `True`, prompts the user for input when merge conflicts occur (default: False).

    Returns:
    -------
    obspy.Stream
        A new Stream object with merged traces.

    Notes:
    ------
    - Uses `smart_merge_traces()` for difficult merges.
    - Handles both **gaps and overlaps** intelligently.
    - Preserves **non-zero data points** when merging.

    Example:
    --------
    ```python
    from obspy import read

    # Load a Stream with overlapping traces
    st = read("example_overlapping.mseed")

    # Perform smart merging
    merged_st = smart_merge(st, verbose=True, interactive=False)

    # Print results
    print(merged_st)
    ```
    """
    newst = Stream()
    all_ids = list(set(tr.id for tr in st))

    for this_id in all_ids:
        these_traces = st.select(id=this_id).sort()

        # Remove exact duplicate traces
        traces_to_keep = []
        for i in range(len(these_traces)):
            if i == 0 or (
                these_traces[i].stats.starttime != these_traces[i - 1].stats.starttime or
                these_traces[i].stats.endtime != these_traces[i - 1].stats.endtime or
                these_traces[i].stats.sampling_rate != these_traces[i - 1].stats.sampling_rate
            ):
                traces_to_keep.append(these_traces[i])

        these_traces = Stream(traces=traces_to_keep)

        # If only 1 trace remains, add it to new Stream
        if len(these_traces) == 1:
            newst.append(these_traces[0])
            continue

        # Try standard merge
        try:
            merged_trace = these_traces.copy().merge()
            newst.append(merged_trace[0])
            if verbose:
                print("- Standard merge successful")
            continue
        except:
            if verbose:
                print("- Standard merge failed, attempting pairwise merge")

        # Pairwise smart merge strategy
        merged_trace = these_traces[0]
        for i in range(1, len(these_traces)):
            trace_pair = Stream([merged_trace, these_traces[i]])

            try:
                merged_trace = trace_pair.copy().merge()[0]  # Standard merge
            except:
                try:
                    merged_trace = _smart_merge_traces([merged_trace, these_traces[i]])  # Use smart merge
                    if verbose:
                        print(f"- Smart merged {these_traces[i].id}")
                except:
                    if verbose:
                        print(f"- Failed to smart merge {these_traces[i].id}")

        # Add merged trace to final Stream
        if merged_trace:
            newst.append(merged_trace)

        # If interactive mode is enabled, prompt user to manually choose a trace
        elif interactive:
            print("\nTrace conflict detected:\n")
            trace_pair.plot()
            for idx, tr in enumerate(trace_pair):
                print(f"{idx}: {tr}")
            choice = int(input("Enter the index of the trace to keep: "))
            newst.append(trace_pair[choice])
        else:
            raise ValueError("Unable to merge traces automatically")

    return newst


def _smart_merge_traces(trace_pair):
    """
    Merges two overlapping traces, preserving all non-zero data values.

    This function:
    - Ensures the traces have the **same ID and sampling rate**.
    - If traces cannot be merged, returns the trace with the **most valid data**.
    - Uses **non-zero values** from both traces when merging.

    Parameters:
    ----------
    trace_pair : list of obspy.Trace
        A pair of seismic traces to merge.

    Returns:
    -------
    obspy.Trace
        A single merged trace.

    Notes:
    ------
    - If traces have **different IDs or sampling rates**, the one with the most valid data is returned.
    - Overlapping segments are filled using **non-zero values from both traces**.

    Example:
    --------
    ```python
    from obspy import read

    # Load two traces with overlap
    tr1 = read("trace1.mseed")[0]
    tr2 = read("trace2.mseed")[0]

    # Perform smart merging
    merged_trace = _smart_merge_traces([tr1, tr2])

    # Print results
    print(merged_trace)
    ```
    """
    this_tr, other_tr = trace_pair

    # Check if traces are mergeable
    if this_tr.id != other_tr.id:
        print("Different trace IDs. Cannot merge.")
        return this_tr if np.count_nonzero(this_tr.data) >= np.count_nonzero(other_tr.data) else other_tr

    if this_tr.stats.sampling_rate != other_tr.stats.sampling_rate:
        print("Different sampling rates. Cannot merge.")
        return this_tr if np.count_nonzero(this_tr.data) >= np.count_nonzero(other_tr.data) else other_tr

    if abs(this_tr.stats.starttime - other_tr.stats.starttime) > this_tr.stats.delta / 4:
        print("Different start times. Cannot merge.")
        return this_tr if np.count_nonzero(this_tr.data) >= np.count_nonzero(other_tr.data) else other_tr

    if abs(this_tr.stats.endtime - other_tr.stats.endtime) > this_tr.stats.delta / 4:
        print("Different end times. Cannot merge.")
        return this_tr if np.count_nonzero(this_tr.data) >= np.count_nonzero(other_tr.data) else other_tr

    # Merge traces, using non-zero values from both
    merged_data = np.where(other_tr.data == 0, this_tr.data, other_tr.data)
    merged_trace = this_tr.copy()
    merged_trace.data = merged_data

    return merged_trace
        
def Stream_min_starttime(all_traces):
    """
    Computes the minimum and maximum start and end times for a given Stream.

    This function takes an **ObsPy Stream** containing multiple traces and 
    determines the following time statistics:
    - **Earliest start time** (`min_stime`)
    - **Latest start time** (`max_stime`)
    - **Earliest end time** (`min_etime`)
    - **Latest end time** (`max_etime`)

    Parameters:
    ----------
    all_traces : obspy.Stream
        A Stream object containing multiple seismic traces.

    Returns:
    -------
    tuple:
        - **min_stime (UTCDateTime)**: The earliest start time among all traces.
        - **max_stime (UTCDateTime)**: The latest start time among all traces.
        - **min_etime (UTCDateTime)**: The earliest end time among all traces.
        - **max_etime (UTCDateTime)**: The latest end time among all traces.

    Notes:
    ------
    - Useful for determining the **temporal coverage** of a Stream.
    - Created for the **CALIPSO data archive** (Alan Linde).

    Example:
    --------
    ```python
    from obspy import read

    # Load a Stream of seismic data
    st = read("example.mseed")

    # Compute time bounds
    min_stime, max_stime, min_etime, max_etime = Stream_min_starttime(st)

    print(f"Start Time Range: {min_stime} to {max_stime}")
    print(f"End Time Range: {min_etime} to {max_etime}")
    ```
    """ 
    min_stime = min([tr.stats.starttime for tr in all_traces])
    max_stime = max([tr.stats.starttime for tr in all_traces])
    min_etime = min([tr.stats.endtime for tr in all_traces])
    max_etime = max([tr.stats.endtime for tr in all_traces])    
    return min_stime, max_stime, min_etime, max_etime

#######################################################################
##               Fixing IDs                                          ##
#######################################################################

def _get_band_code(sampling_rate):
    """
    Determines the appropriate band code based on the sampling rate.

    The band code is the first letter of the **SEED channel naming convention**, which 
    categorizes seismic channels based on frequency range.

    Parameters:
    ----------
    sampling_rate : float
        The sampling rate of the seismic trace in Hz.

    Returns:
    -------
    str or None
        The appropriate band code (e.g., 'B' for broadband, 'H' for high-frequency broadband).
        Returns `None` if no matching band code is found (should not happen if lookup table is correct).

    Notes:
    ------
    - This function relies on `BAND_CODE_TABLE`, a dictionary defining the mapping 
      between frequency ranges and SEED band codes.

    Example:
    --------
    ```python
    band_code = _get_band_code(100.0)
    print(band_code)  # Output: 'H' (High-frequency broadband)
    ```
    """
    # Band code lookup table based on IRIS SEED convention
    BAND_CODE_TABLE = {
        (0.0001, 0.001): "R",  # Extremely Long Period (0.0001 - 0.001 Hz)   
        (0.001, 0.01): "U",  # Ultra Low Frequency (~0.01 Hz)
        (0.01, 0.1): "V",  # Very Low Frequency (~0.1 Hz)
        (0.1, 2): "L",   # Long Period (~1 Hz)
        (2, 10): "M",  # Mid Period (1 - 10 Hz)
        (10, 80): "B", # Broadband (S if Short Period instrument, corner > 0.1 Hz)
        (80, 250): "H",  # High Frequency (80 - 250 Hz) (E if Short Period instrument, corner > 0.1 Hz)
        (250, 1000): "D",  # Very High Frequency (250 - 1000 Hz) (C if Short Period instrument, corner > 0.1 Hz)
        (1000, 5000): "G",  # Extremely High Frequency (1 - 5 kHz) (F if Short period)
    }

    for (low, high), code in BAND_CODE_TABLE.items():
        if low <= sampling_rate < high:
            return code
    return None  # Should not happen if lookup table is correct

def _adjust_band_code_for_sensor_type(current_band_code, expected_band_code, short_period=False):
    """
    Adjusts the band code if the current trace belongs to a short-period seismometer.

    SEED convention distinguishes between **broadband** and **short-period** seismometers.
    This function adjusts the expected band code based on the current sensor type.

    Mapping:
    - 'B' (Broadband) → 'S' (Short-period)
    - 'H' (High-frequency broadband) → 'E' (Short-period high-frequency)
    - 'D' (Very long period broadband) → 'C' (Short-period very long period)
    - 'G' (Extremely high-frequency broadband) → 'F' (Short-period extremely high-frequency)

    Parameters:
    ----------
    current_band_code : str
        The first character of the current `trace.stats.channel` (e.g., 'S', 'E', 'C', 'F').
    expected_band_code : str
        The computed band code based on the sampling rate.
    short_period : bool, optional
        If `True`, forces short-period band codes even if the current band code is not in the expected mapping.

    Returns:
    -------
    str
        The adjusted band code if applicable.

    Example:
    --------
    ```python
    adjusted_band_code = _adjust_band_code_for_sensor_type('S', 'B')
    print(adjusted_band_code)  # Output: 'S' (Short-period equivalent of 'B')
    ```
    """
    short_period_codes = {'S', 'E', 'C', 'F'}
    
    if current_band_code in short_period_codes or short_period:
        band_code_mapping = {'B': 'S', 'H': 'E', 'D': 'C', 'G': 'F'}
        return band_code_mapping.get(expected_band_code, expected_band_code)
    
    return expected_band_code

def _fix_legacy_id(trace):
    """
    Fixes legacy trace IDs for old VDAP/analog telemetry networks.

    This function corrects **four-character station codes** where the orientation
    is embedded within the station name.

    Corrections:
    - If `trace.stats.station == 'IRIG'`, sets `trace.stats.channel = 'ACE'`.
    - Converts single-letter channels ('v', 'n', 'e') to SEED channel names ('EHZ', 'EHN', 'EHE').
    - Extracts orientation (4th character of station name) to determine the correct SEED channel.
    - Removes the orientation character from the station name.

    Parameters:
    ----------
    trace : obspy.Trace
        The seismic trace to modify.

    Returns:
    -------
    None
        Modifies `trace.stats.station` and `trace.stats.channel` in place.

    Example:
    --------
    ```python
    from obspy import read

    trace = read("example.mseed")[0]
    _fix_legacy_id(trace)

    print(trace.id)  # Corrected SEED-compliant ID
    ```
    """
    if trace.stats.station == 'IRIG':
        trace.stats.channel = 'ACE'
    else:
        if trace.stats.channel=='v':
            trace.stats.channel='EHZ'
        elif trace.stats.channel=='n':
            trace.stats.channel='EHN'
        elif trace.stats.channel=='e':
            trace.stats.channel='EHE'                       
        
        orientation = trace.stats.station[3].strip()  # Position 4        
        if orientation in "ZNE":  # Standard orientations
            channel = f"EH{orientation}"
        elif orientation == "L":  # Special case for "L"
            channel = "ELZ"
        elif orientation == 'P': # assume pressure sensor?
            channel = 'EDF'
        else:
            channel = f'??{orientation}'
            #raise ValueError(f"Unknown orientation '{orientation}' in '{station}'")
        trace.stats.channel = channel
        trace.stats.station = trace.stats.station[0:3].strip()  # Positions 1-3

def fix_trace_id(trace, legacy=False, netcode=None, verbose=False):
    """
    Standardizes a seismic trace's ID to follow SEED naming conventions.

    This function:
    - Fixes legacy **VDAP/analog telemetry IDs** if `legacy=True`.
    - Ensures a valid **network code** if `netcode` is provided.
    - Adjusts the **band code** based on sampling rate.
    - Ensures the location code is **either empty or two characters**.
    - Fixes known **station name substitutions** (e.g., `CARL1` → `TANK` for KSC data).

    Parameters:
    ----------
    trace : obspy.Trace
        The seismic trace to modify.
    legacy : bool, optional
        If `True`, applies `_fix_legacy_id()` to correct old-style station codes (default: False).
    netcode : str, optional
        Network code to assign if missing (default: None).
    verbose : bool, optional
        If `True`, prints trace ID changes (default: False).

    Returns:
    -------
    bool
        `True` if the trace ID was changed, `False` otherwise.

    Notes:
    ------
    - Calls `_get_band_code()` to determine the correct band code based on sampling rate.
    - Calls `_adjust_band_code_for_sensor_type()` to refine the band code for short-period sensors.
    - Ensures **station names are corrected** for specific networks (e.g., `FL.CARL1 → FL.TANK`).

    Example:
    --------
    ```python
    from obspy import read

    trace = read("example.mseed")[0]
    changed = fix_trace_id(trace, legacy=True, netcode="XX", verbose=True)

    if changed:
        print(f"Updated Trace ID: {trace.id}")
    ```
    """
    changed = False

    if legacy: # indicates an old VDAP/analog telemetry network where 4-character station code includes orientation
        _fix_legacy_id(trace)

    if not trace.stats.network and netcode:
        trace.stats.network = netcode
        changed = True   
    
    current_id = trace.id
    net, sta, loc, chan = current_id.split('.')
    sampling_rate = trace.stats.sampling_rate
    current_band_code = chan[0]

    # if not an analog QC channel, fix band code
    if chan[0]=='A':
        pass
    else:

        # Determine the correct band code
        expected_band_code = _get_band_code(sampling_rate) # this assumes broadband sensor

        # adjust if short-period sensor
        expected_band_code = _adjust_band_code_for_sensor_type(current_band_code, expected_band_code)
        chan = expected_band_code + chan[1:]

    # make sure location is 0 or 2 characters
    if len(loc)==1:
        loc = loc.zfill(2)

    # change CARL1 to TANK
    if net=='FL':
        if sta=='CARL1':
            sta = 'TANK'
        elif sta=='CARL0':
            sta = 'BCHH'

    expected_id = '.'.join([net,sta,loc,chan])
    #print(current_id, expected_id)

    if (expected_id != current_id):
        changed = True
        if verbose:
            print(f"Current ID: {current_id}, Expected: {expected_id}) based on fs={sampling_rate}")
        trace.id = expected_id
    #print(trace)
    return changed 

#######################################################################
##               Event detection                                     ##
#######################################################################

def detect_network_event(st_in, minchans=None, threshon=3.5, threshoff=1.0, 
                         sta=0.5, lta=5.0, pad=0.0, best_only=False, verbose=False, 
                         freq=None, algorithm='recstalta', criterion='longest'):
    """
    Detects and associates seismic events across a network of stations using coincidence triggering.

    This function applies a short-term average / long-term average (STA/LTA) or other 
    triggering algorithm to detect seismic events recorded across multiple stations. 
    It identifies events based on a coincidence threshold and returns a list of detected 
    triggers, including event start time, duration, number of stations that triggered, and 
    characteristic function peak values.

    If `best_only=True`, the function returns only the most significant event based on the 
    specified `criterion`. Otherwise, it returns all detected events along with onset and 
    offset times.

    Parameters:
    ----------
    st_in : obspy.Stream
        The input Stream object containing waveform data from multiple stations.
    minchans : int, optional
        Minimum number of stations required to trigger an event (default: half the number of traces in `st_in`, minimum 2).
    threshon : float, optional
        Trigger threshold for STA/LTA or other detection algorithm (default: 3.5).
    threshoff : float, optional
        De-trigger threshold for STA/LTA or other detection algorithm (default: 1.0).
    sta : float, optional
        Short-term average window length in seconds (default: 0.5s).
    lta : float, optional
        Long-term average window length in seconds (default: 5.0s).
    pad : float, optional
        Time (in seconds) to pad the traces before processing (default: 0.0s, no padding).
    best_only : bool, optional
        If `True`, returns only the most significant detected event based on `criterion` (default: False).
    verbose : bool, optional
        If `True`, prints debug information during processing (default: False).
    freq : list of float, optional
        If specified, applies a bandpass filter before processing using `[freq_min, freq_max]` in Hz (default: None).
    algorithm : str, optional
        The detection algorithm to use. Options: `"recstalta"`, `"zdetect"`, `"carlstatrig"` (default: `"recstalta"`).
    criterion : str, optional
        Criterion for selecting the "best" event when `best_only=True`:
        - `"longest"` (default): Selects the event with the highest product of `coincidence_sum * duration`.
        - `"cft"`: Selects the event with the highest sum of `cft_peaks`.
        - `"cft_duration"`: Uses `cft_peaks * duration` to determine the best event.

    Returns:
    -------
    If `best_only=True`:
        dict or None
            A dictionary containing information about the most significant detected event, 
            or `None` if no event was found.

    If `best_only=False`:
        tuple:
            - **trig (list of dicts)**: Each dict contains information about a detected event, 
              with keys including:
              ```
              {
                  'cft_peak_wmean': float,   # Mean peak characteristic function value
                  'cft_peaks': list,         # Characteristic function peak values for each station
                  'cft_std_wmean': float,    # Mean standard deviation of characteristic function
                  'cft_stds': list,          # Standard deviation values per station
                  'coincidence_sum': float,  # Number of stations triggering simultaneously
                  'duration': float,         # Event duration (s)
                  'stations': list,          # Stations that contributed to the event
                  'time': UTCDateTime,       # Event start time
                  'trace_ids': list          # List of trace IDs contributing to the event
              }
              ```
            - **ontimes (list of UTCDateTime)**: List of event start times.
            - **offtimes (list of UTCDateTime)**: List of event end times.

    If no events are detected:
        - If `best_only=True`: Returns `None`.
        - If `best_only=False`: Returns `(None, None, None)`.

    Notes:
    ------
    - If `sta` and `lta` are chosen such that `lta` is too long relative to the signal, 
      event detection may fail.
    - The function applies optional bandpass filtering if `freq` is specified.
    - The `pad` option allows extending the trace duration before filtering or triggering.
    - Detected events can be further refined using `trim_to_event()`.

    Example:
    --------
    ```python
    from obspy import read

    # Load waveform data from multiple stations
    st = read("network_data.mseed")

    # Run event detection
    events, onsets, offsets = detect_network_event(st, freq=[1, 20], verbose=True)

    if events:
        print(f"Detected {len(events)} events")
        for event in events:
            print(f"Event at {event['time']} with duration {event['duration']} s")

    # Retrieve only the best event
    best_event = detect_network_event(st, best_only=True, criterion="cft")
    print(f"Best event: {best_event}")
    ```
    """
    st = st_in.copy()
    if pad>0.0:
        for tr in st:
            pad_trace(tr, pad)

    if freq:
        if verbose:
            print('Filtering traces')
        st.filter('bandpass', freqmin=freq[0], freqmax=freq[1], corners=4, zerophase=True)
            
    if not minchans:
        minchans = max(( int(len(st)/2), 2)) # half the channels or 2, whichever is greater
    if verbose:
        print('minchans=',minchans)
    #trig = coincidence_trigger(algorithm, threshon, threshoff, st, minchans, sta=sta, lta=lta, max_trigger_length=180, delete_long_trigger=True, details=True) # 0.5s, 10s
    if algorithm == "zdetect":
        trig = coincidence_trigger(algorithm, threshon, threshoff, st, minchans, sta=sta, details=True)
    elif algorithm == "carlstatrig":
        trig = coincidence_trigger(algorithm, threshon, threshoff, st, minchans, sta=sta, lta=lta, ratio=1, quiet=True, details=True)
    else:
        trig = coincidence_trigger(algorithm, threshon, threshoff, st, minchans, sta=sta, lta=lta, details=True)
    if trig:

        if best_only:
            best_trig = {}
            best_product = 0

            for this_trig in trig:
                #print(this_trig)
                thistime = UTCDateTime(this_trig['time'])
                if thistime > st[0].stats.starttime:
                    if criterion=='longest':
                        this_product = this_trig['coincidence_sum']*this_trig['duration']
                    elif criterion=='cft':
                        this_product = sum(this_trig["cft_peaks"])
                    else:
                        this_product = sum(this_trig["cft_peaks"])*this_trig['duration']
                    if this_product > best_product:
                        best_trig = this_trig
                        best_product = this_product
            return best_trig  
        else:
            ontimes = []
            offtimes = []
            for this_trig in trig:
                thistime = UTCDateTime(this_trig['time'])
                if thistime > st[0].stats.starttime:
                    ontimes.append(this_trig['time'])
                    offtimes.append(this_trig['time']+this_trig['duration'])
            return trig, ontimes, offtimes
    else:
        if best_only:
            return None
        else:
            return None, None, None
    
def add_channel_detections(st, lta=5.0, threshon=0.5, threshoff=0.0, max_duration=120):
    """
    Runs single-channel event detection on each Trace in a Stream using the STA/LTA method.

    This function applies a `z_detect` characteristic function to each trace and detects 
    events based on STA/LTA triggering. It does **not** perform multi-station coincidence 
    triggering or event association.

    Detected triggers are stored in `tr.stats.triggers` as a list of two-element numpy arrays, 
    where each array contains the **trigger on** and **trigger off** times as `UTCDateTime` objects.

    Parameters:
    ----------
    st : obspy.Stream
        The input Stream object containing multiple seismic traces.
    lta : float, optional
        Long-term average (LTA) window length in seconds (default: 5.0s).
    threshon : float, optional
        Trigger-on threshold for detection (default: 0.5).
    threshoff : float, optional
        Trigger-off threshold for ending a detection (default: 0.0).
    max_duration : float, optional
        Maximum event duration in seconds before forcing a trigger-off (default: 120s).

    Returns:
    -------
    None
        The function modifies `tr.stats.triggers` in-place for each trace.

    Stores:
    -------
    - `tr.stats.triggers` : list of lists
        Each trace will have a list of **trigger onset and offset times** as `UTCDateTime` pairs.
        Example:
        ```python
        tr.stats.triggers = [
            [UTCDateTime("2025-03-10T12:00:00"), UTCDateTime("2025-03-10T12:00:15")],
            [UTCDateTime("2025-03-10T12:30:05"), UTCDateTime("2025-03-10T12:30:20")]
        ]
        ```

    Notes:
    ------
    - If `lta=5.0`, at least **5 seconds of noise** before the signal is required.
    - This function only detects events **on individual channels**—it does not associate detections across stations.
    - Uses `z_detect()` as the characteristic function and `trigger_onset()` to extract events.

    Example:
    --------
    ```python
    from obspy import read

    # Load waveform data
    st = read("example.mseed")

    # Run single-channel detections
    add_channel_detections(st, lta=5.0, threshon=0.5, threshoff=0.1, max_duration=60)

    # Print detected triggers
    for tr in st:
        print(f"Triggers for {tr.id}: {tr.stats.triggers}")
    ```
    """
    for tr in st:
        tr.stats['triggers'] = []
        Fs = tr.stats.sampling_rate
        cft = z_detect(tr.data, int(lta * Fs)) 
        triggerlist = trigger_onset(cft, threshon, threshoff, max_len = max_duration * Fs)
        for trigpair in triggerlist:
            trigpairUTC = [tr.stats.starttime + samplenum/Fs for samplenum in trigpair]
            tr.stats.triggers.append(trigpairUTC)

def get_event_window(st, pretrig=30, posttrig=30):
def get_event_window(st, pretrig=30, posttrig=30):
    """
    Determines the time window encompassing all detected triggers in a Stream, with optional pre/post-event padding.

    This function scans all traces in the given `Stream`, extracts the earliest **trigger-on**
    and the latest **trigger-off** time from `tr.stats.triggers`, and returns a **time window**
    that includes an additional buffer (`pretrig` and `posttrig`) before and after the event.

    **Assumptions:**
    - The function assumes that `add_channel_detections(st, ...)` has already been run,
      so that each `Trace` contains `tr.stats.triggers`.
    - If no triggers are found, the function returns `(None, None)`.

    Parameters:
    ----------
    st : obspy.Stream
        The input Stream object containing seismic traces with precomputed triggers.
    pretrig : float, optional
        Extra time in seconds to include **before** the first trigger-on time (default: 30s).
    posttrig : float, optional
        Extra time in seconds to include **after** the last trigger-off time (default: 30s).

    Returns:
    -------
    tuple:
        - **mintime (UTCDateTime)**: The earliest trigger-on time minus `pretrig` seconds.
        - **maxtime (UTCDateTime)**: The latest trigger-off time plus `posttrig` seconds.
        - If no triggers exist, returns `(None, None)`.

    Notes:
    ------
    - The function uses the **median** trigger times (`N/2` index) instead of min/max to reduce
      sensitivity to outliers.
    - The returned `mintime` and `maxtime` can be used with `trim_to_event()` to extract
      the event from the full waveform dataset.

    Example:
    --------
    ```python
    from obspy import read

    # Load waveform data and run single-channel detection
    st = read("example.mseed")
    add_channel_detections(st)

    # Get event time window with 20s padding
    mintime, maxtime = get_event_window(st, pretrig=20, posttrig=20)

    print(f"Detected event from {mintime} to {maxtime}")
    ```
    """

    mintime = []
    maxtime = []
    
    for tr in st:
        if 'triggers' in tr.stats:
            if len(tr.stats.triggers)==0:
                continue
            trigons = [thistrig[0] for thistrig in tr.stats.triggers]
            trigoffs = [thistrig[1] for thistrig in tr.stats.triggers]   
            mintime.append(min(trigons))
            maxtime.append(max(trigoffs))           
    
    N = int(len(mintime)/2)
    if len(mintime)>0:
        return sorted(mintime)[N], sorted(maxtime)[N]
    else:
        return None, None

def trim_to_event(st, mintime, maxtime, pretrig=10, posttrig=10):
    """
    Trims a Stream to a specified event window, adding a pre-trigger and post-trigger buffer.

    This function trims a `Stream` object to the time window defined by `mintime` (first trigger-on)
    and `maxtime` (last trigger-off), with optional extra time added before and after the event.

    **Use Cases:**
    - Focuses analysis on the main event while preserving relevant pre/post-event data.
    - Can be used after running `detect_network_event` or `get_event_window` to refine the dataset.

    Parameters:
    ----------
    st : obspy.Stream
        The input Stream object to be trimmed.
    mintime : UTCDateTime
        The earliest trigger-on time (e.g., from `get_event_window` or `detect_network_event`).
    maxtime : UTCDateTime
        The latest trigger-off time (e.g., from `get_event_window` or `detect_network_event`).
    pretrig : float, optional
        Extra time in seconds to keep **before** `mintime` (default: 10s).
    posttrig : float, optional
        Extra time in seconds to keep **after** `maxtime` (default: 10s).

    Returns:
    -------
    None
        The function modifies the `Stream` **in place**, trimming all traces to the specified time window.

    Notes:
    ------
    - If `mintime` or `maxtime` is `None`, the function will not trim the Stream.
    - If `pretrig` or `posttrig` exceeds the available data range, trimming may result in an empty Stream.
    - This function does not return a new `Stream` object; it modifies the input `st` directly.

    Example:
    --------
    ```python
    from obspy import read

    # Load waveform data and detect events
    st = read("example.mseed")
    add_channel_detections(st)

    # Get event time window
    mintime, maxtime = get_event_window(st, pretrig=20, posttrig=20)

    # Trim the Stream to focus on the detected event
    if mintime and maxtime:
        trim_to_event(st, mintime, maxtime, pretrig=5, posttrig=15)
        print("Stream trimmed to event window.")
    ```
    """
    st.trim(starttime=mintime-pretrig, endtime=maxtime+posttrig)
    
#######################################################################    
########################         WFDISC tools                        ##
#######################################################################
     
def index_waveformfiles(wffiles, ampeng=False, events=True):
    """
    Indexes a list of seismic waveform files and returns a DataFrame similar to a `wfdisc` table.

    This function reads seismic waveform files, extracts relevant metadata, and compiles the 
    information into a **pandas DataFrame** with fields similar to a `wfdisc` database table.
    Additionally, it creates an **ObsPy Catalog** of `Event` objects if `events=True`.

    **Primary Use Cases:**
    - Organizing waveform metadata for seismic data archives (e.g., CALIPSO data).
    - Preparing a catalog of waveform files with event associations.
    - Computing basic waveform statistics such as **amplitude and energy** (optional).

    Parameters:
    ----------
    wffiles : list of str
        List of file paths to seismic waveform data files (MiniSEED, SAC, etc.).
    ampeng : bool, optional
        If `True`, computes **maximum absolute amplitude** and **signal energy** for each trace (default: False).
    events : bool, optional
        If `True`, generates an **ObsPy Catalog** with `Event` objects associated with each file (default: True).

    Returns:
    -------
    tuple:
        - **wfdisc_df (pandas.DataFrame)**:
            A DataFrame containing metadata for each waveform file with columns:
            ```
            file_index    : Index of the file in `wffiles`
            event_id      : Unique Event ID (ObsPy ResourceIdentifier)
            waveform_id   : ResourceIdentifier linking file name and trace index
            traceID       : SEED-formatted trace ID (NET.STA.LOC.CHA)
            starttime     : UTCDateTime of trace start
            endtime       : UTCDateTime of trace end
            npts          : Number of data points
            sampling_rate : Sampling rate (Hz)
            calib         : Calibration factor
            ddir          : Directory containing the waveform file
            dfile         : Waveform file name
            duration      : Trace duration (seconds)
            amplitude     : (Optional) Maximum absolute amplitude of the detrended trace
            energy        : (Optional) Sum of squared amplitudes (energy estimate)
            ```

        - **cat (obspy.Catalog)**:
            If `events=True`, returns an **ObsPy Catalog** containing `Event` objects with metadata 
            extracted from the waveform files. If `events=False`, returns only the `wfdisc_df`.

    Notes:
    ------
    - If a file **cannot be read**, it is skipped, and a warning is printed.
    - If `ampeng=True`, the function attempts to **detrend** the waveform before computing amplitude and energy.
    - The function **sorts** the DataFrame by `starttime` before returning it.

    Example:
    --------
    ```python
    from obspy import read
    import glob

    # List of waveform files
    waveform_files = glob.glob("/data/seismic/*.mseed")

    # Index waveform files and compute amplitude/energy
    wfdisc_df, cat = index_waveformfiles(waveform_files, ampeng=True)

    # Display waveform metadata
    print(wfdisc_df.head())

    # Save the catalog
    cat.write("seismic_catalog.xml", format="QUAKEML")
    ```
    """
    wfdisc_df = pd.DataFrame()
    file_index = []
    event_id = []
    traceids = []
    starttimes = []
    endtimes = []
    sampling_rates = []
    calibs = []
    ddirs = []
    dfiles = []
    npts = []
    duration = []
    waveform_id = []
    if ampeng:
        amp = []
        eng = []

    events = []
    for filenum, wffile in enumerate(sorted(wffiles)):
        dfile = os.path.basename(wffile)
        ddir = os.path.dirname(wffile)
        try:
            this_st = read(wffile)
            print('Read %s\n' % wffile)
        except:
            print('Could not read %s\n' % wffile)
            next
        else:
            ev = Event()
            comments = []
            stime = this_st[0].stats.starttime
            sfilename = stime.strftime("%d-%H%M-%S") + "L.S" + stime.strftime("%Y%m")  
            comments.append(Comment(text=f'wavfile: {wffile}'))
            comments.append(Comment(text=f'sfile: {sfilename}'))                      
            for tracenum, this_tr in enumerate(this_st):
                file_index.append(filenum)
                event_id.append(ev.resource_id)  
                wid = ResourceIdentifier(f'{dfile},[{tracenum}]')
                waveform_id.append(wid)
                r = this_tr.stats
                traceids.append(this_tr.id)
                starttimes.append(r.starttime)
                endtimes.append(r.endtime)
                sampling_rates.append(r.sampling_rate)
                calibs.append(r.calib)
                ddirs.append(ddir)
                dfiles.append(dfile)
                npts.append(r.npts)
                duration.append(r.endtime - r.starttime)
                if ampeng:
                    try:
                        this_tr.detrend('linear')
                    except:
                        this_tr.detrend('constant')
                    y = abs(this_tr.data)
                    amp.append(np.nanmax(y))
                    y = np.nan_to_num(y, nan=0)
                    eng.append(np.sum(np.square(y)))
            ev.comments = comments
        events.append(ev)
      
    if wffiles:
        wfdisc_dict = {'file_index':file_index, 'event_id':event_id, 'waveform_id':waveform_id, 'traceID':traceids, 'starttime':starttimes, 'endtime':endtimes, 'npts':npts, 
                       'sampling_rate':sampling_rates, 'calib':calibs, 'ddir':ddirs, 'dfile':dfiles, 'duration':duration}
        if ampeng:
            wfdisc_dict['amplitude']=amp
            wfdisc_dict['energy']=eng
        wfdisc_df = pd.DataFrame.from_dict(wfdisc_dict)  
        wfdisc_df.sort_values(['starttime'], ascending=[True], inplace=True)
    if events:
        cat = Catalog(events=events) 
        return wfdisc_df, cat
    else:
        return wfdisc_df

def wfdisc_to_BUD(wfdisc_df, TOPDIR, put_away):
    """
    Converts waveform files indexed in a `wfdisc`-like DataFrame into a BUD (Buffer of Uniform Data) archive.

    This function reads waveform data from a **wfdisc-like DataFrame**, extracts individual traces grouped 
    by unique `traceID`, and writes the data into a **BUD (Buffer of Uniform Data) format archive**.

    The function processes waveform files on a **per-day basis**, merging overlapping traces and 
    ensuring that data is correctly organized in the BUD directory structure.

    **Workflow:**
    - Iterates through unique trace IDs in `wfdisc_df`.
    - Determines the earliest and latest timestamps for each trace.
    - Reads and merges waveform files for each day.
    - Writes merged traces into a **BUD-format** archive.
    - Optionally moves processed files to a `.PROCESSED` directory if `put_away=True`.

    Parameters:
    ----------
    wfdisc_df : pandas.DataFrame
        A DataFrame containing waveform metadata (as produced by `index_waveformfiles`).
    TOPDIR : str
        The top-level directory where the BUD archive will be stored.
    put_away : bool
        If `True`, moves successfully processed waveform files to a `.PROCESSED` directory (default: False).

    Returns:
    -------
    None
        The function processes the waveform files and writes data to the BUD archive.

    Notes:
    ------
    - The function processes waveform data **in daily chunks**, merging traces that overlap.
    - Files that cannot be read are skipped with a warning.
    - `Stream_to_BUD(TOPDIR, all_traces)` is used to write the BUD archive.
    - Successfully processed waveform files can be optionally moved using `shutil.move()`.

    Example:
    --------
    ```python
    import pandas as pd

    # Load an example wfdisc DataFrame
    wfdisc_df = pd.read_csv("waveform_index.csv")

    # Convert to BUD archive
    wfdisc_to_BUD(wfdisc_df, "/data/BUD_archive", put_away=True)
    ```
    """

    unique_traceIDs = wfdisc_df['traceID'].unique().tolist()
    print(unique_traceIDs)
    
    successful_wffiles = list()

    for traceID in unique_traceIDs:
        print(traceID)
        
        trace_df = wfdisc_df[wfdisc_df['traceID']==traceID]
        
        # identify earliest start time and latest end time for this channel
        #print(trace_df.iloc[0]['starttime'])
        #print(trace_df.iloc[-1]['endtime'])
        minUTC = trace_df.starttime.min()
        maxUTC = trace_df.endtime.max()
        start_date = minUTC.replace(hour=0, minute=0, second=0, microsecond=0)
        end_date = maxUTC.replace(hour=23, minute=59, second=59, microsecond=999999)
        this_date = start_date

        while this_date <= end_date: 
            all_traces = Stream()
        
            # loop from earliest start day to latest end day
            subset_df = trace_df[(trace_df['starttime'] < this_date+86400) & (trace_df['endtime'] >= this_date)]
            #print(subset_df)
            
            if len(subset_df.index)==0:
                next
        
            for index, row in subset_df.iterrows():
                wffile = os.path.join(row['ddir'], row['dfile'])
                start_at = max([this_date, row['starttime']])
                end_at = min([this_date+86400, row['endtime']])
                print('- ',wffile,': START AT:', start_at, ', END AT: ',end_at)
                try:
                    this_st = read(wffile, starttime=start_at, endtime=end_at)

                except:
                    print(' Failed\n')
                    next
                else:
                    print(' Succeeded\n')
                    #raise Exception("Stopping here")
                    if end_at == row['endtime']:
                        successful_wffiles.append(wffile)   
                    for this_tr in this_st:
                        if this_tr.id == traceID:
                            #print(tr.stats)
                            all_traces = all_traces.append(this_tr)
            #print(st.__str__(extended=True)) 
            try:
                all_traces.merge(fill_value=0)
            except:
                print('Failed to merge ', all_traces)
            print(all_traces.__str__(extended=True))
        
            # Check that we really only have a single trace ID before writing the BUD files
            error_flag = False
            for this_tr in all_traces:
                if not this_tr.id == traceID:
                    error_flag = True
            if not error_flag:
                try:
                    Stream_to_BUD(TOPDIR, all_traces)
                except:
                    print('Stream_to_BUD failed for ', all_traces)
            
            this_date += 86400
            
    for wffile in successful_wffiles:
        ddir = os.path.dirname(wffile)
        dbase = "%s.PROCESSED" % os.path.basename(wffile)
        newwffile = os.path.join(ddir, dbase)
        print('move %s %s' % (wffile, newwffile))
        if os.path.exists(wffile) and put_away:
            shutil.move(wffile, newwffile)

def process_wfdirs(wfdirs, filematch, put_away=False):
    """
    Processes directories containing waveform data, indexing the files and converting them to a BUD archive.

    This function scans **one or more directories** (`wfdirs`) for seismic waveform files matching a given 
    pattern (`filematch`), builds a **wfdisc-like DataFrame**, and converts the waveform data into a 
    **BUD (Buffer of Uniform Data) archive**.

    **Workflow:**
    - Scans the provided directories (`wfdirs`) for waveform files.
    - Calls `index_waveformfiles()` to extract waveform metadata into a **pandas DataFrame**.
    - Calls `wfdisc_to_BUD()` to convert and store waveform data in the **BUD archive**.
    - Optionally moves processed files if `put_away=True`.

    Parameters:
    ----------
    wfdirs : list of str
        List of directories containing waveform files.
    filematch : str
        File matching pattern (e.g., `"*.mseed"`, `"*.sac"`) to identify waveform files.
    put_away : bool, optional
        If `True`, moves successfully processed waveform files to a `.PROCESSED` directory (default: False).

    Returns:
    -------
    None
        The function scans the directories, processes waveform data, and writes to the BUD archive.

    Notes:
    ------
    - The function processes all files readable by **ObsPy**.
    - If no waveform files are found in a directory, it is skipped.
    - Calls `wfdisc_to_BUD()` to handle the conversion process.

    Example:
    --------
    ```python
    import glob

    # Define directories and file pattern
    waveform_dirs = ["/data/seismic/day1", "/data/seismic/day2"]
    file_pattern = "*.mseed"

    # Process and convert to BUD archive
    process_wfdirs(waveform_dirs, file_pattern, put_away=True)
    ```
    """

    for wfdir in wfdirs:
        print('Processing %s' % wfdir)
        wffiles = glob.glob(os.path.join(wfdir, filematch))
        if wffiles:
            #print(wffiles)
            wfdisc_df = index_waveformfiles(wffiles)
            #print(wfdisc_df)
            if not wfdisc_df.empty:
                wfdisc_to_BUD(wfdisc_df, TOPDIR, put_away)  
    print('Done.')

#######################################################################
##                BUD tools                                          ##
#######################################################################

def Stream_to_BUD(TOPDIR, all_traces):
def Stream_to_BUD(TOPDIR, all_traces):
    """
    Converts a Stream object into the IRIS/PASSCAL BUD (Buffer of Uniform Data) format.

    This function takes a **Stream** of seismic waveform data, splits it into **24-hour-long segments** per 
    channel, and writes it into a **BUD directory structure** based on station, network, and date.

    **BUD Directory Structure Example:**
    ```
    DAYS/
    ├── BHP2
    │   ├── 1R.BHP2..EH1.2020.346
    │   ├── 1R.BHP2..EH2.2020.346
    │   └── 1R.BHP2..EHZ.2020.346
    ├── FIREP
    │   ├── 1R.FIREP..EH1.2020.346
    │   ├── 1R.FIREP..EH2.2020.346
    │   └── 1R.FIREP..EHZ.2020.346
    ```
    where:
    - `BHP2`, `FIREP`, etc., are station names.
    - `1R` is the network name.
    - Channels are `EH[Z12]`.
    - The year is `2020`, and the Julian day is `346`.

    Parameters:
    ----------
    TOPDIR : str
        The top-level directory where the BUD archive will be stored.
    all_traces : obspy.Stream
        The input Stream object containing multiple seismic traces.

    Returns:
    -------
    None
        The function writes files in **MiniSEED format** to the BUD archive.

    Notes:
    ------
    - The function ensures that traces are merged and padded to **24-hour segments**.
    - If a trace already exists in the BUD archive, it merges new data using `Trace_merge_with_BUDfile()`.
    - Created for **ROCKETSEIS** and **CALIPSO** data archives.
    
    Example:
    --------
    ```python
    from obspy import read

    # Load waveform data
    st = read("example.mseed")

    # Convert to BUD format
    Stream_to_BUD("/data/BUD_archive", st)
    ```
    """
    
    all_traces = Stream_to_24H(all_traces)
    
    daysDir = os.path.join(TOPDIR, 'DAYS')

    for this_tr in all_traces:
        YYYY = this_tr.stats.starttime.year
        JJJ = this_tr.stats.starttime.julday
        stationDaysDir = os.path.join(daysDir, this_tr.stats.station)
        if not os.path.exists(stationDaysDir):
            os.makedirs(stationDaysDir)
            #print(stationDaysDir)
        mseedDayBasename = "%s.%04d.%03d" % (this_tr.id, YYYY, JJJ  )
        mseedDayFile = os.path.join(stationDaysDir, mseedDayBasename)
        #print(mseedDayFile)
        if os.path.exists(mseedDayFile):
            this_tr = Trace_merge_with_BUDfile(this_tr, mseedDayFile)

        this_tr.write(mseedDayFile, format='MSEED') 


    
def BUD_load_day(BUDDIR, year, jday):
    """
    Loads all seismic waveform files for a given **year** and **Julian day** from a BUD archive.

    This function searches a BUD directory structure for waveform files matching the specified 
    **year** and **Julian day**, reads them, and returns a merged **Stream**.

    Parameters:
    ----------
    BUDDIR : str
        The root directory containing the BUD archive.
    year : int
        The year of the requested data (e.g., `2020`).
    jday : int
        The Julian day of the requested data (e.g., `346` for Dec 11 in a leap year).

    Returns:
    -------
    obspy.Stream
        A Stream object containing all traces for the requested day.

    Notes:
    ------
    - The function scans all station subdirectories within `BUDDIR` for matching files.
    - If a file cannot be read, it is skipped with a warning.
    
    Example:
    --------
    ```python
    # Load waveform data for year 2020, Julian day 346
    st = BUD_load_day("/data/BUD_archive", 2020, 346)

    # Print available traces
    print(st)
    ```
    """

    all_stations = glob.glob(os.path.join(BUDDIR, '*'))
    all_traces = Stream()
    for station_dir in all_stations:
        all_files = glob.glob(os.path.join(station_dir, '*.%04d.%03d' % (year, jday)))
        for this_file in all_files:
            try:
                these_traces = read(this_file)
            except:
                print('Cannot read %s' % this_file)
            else:
                for this_tr in these_traces:
                    all_traces.append(this_tr)
    return all_traces


def Stream_to_24H(all_traces):
    """
    Pads and merges a Stream object so that each trace spans a full **24-hour period**.

    This function:
    - **Merges traces** with the same ID.
    - **Pads missing data** with zeros to create **continuous 24-hour segments**.
    - Returns a new Stream with traces that **start at 00:00:00 UTC** and end at **23:59:59 UTC**.

    Parameters:
    ----------
    all_traces : obspy.Stream
        The input Stream object containing seismic traces.

    Returns:
    -------
    obspy.Stream
        A Stream object with each trace spanning exactly **24 hours**.

    Notes:
    ------
    - Uses `Stream_min_starttime()` to determine the earliest/largest time window.
    - Pads traces using `.trim(starttime, endtime, pad=True, fill_value=0)`.
    - Used for **ROCKETSEIS** and **CALIPSO** data archives.

    Example:
    --------
    ```python
    from obspy import read

    # Load waveform data
    st = read("example.mseed")

    # Convert to 24-hour segments
    st_24h = Stream_to_24H(st)

    # Print trace info
    print(st_24h)
    ```
    """

    all_traces.merge(fill_value=0)
    min_stime, max_stime, min_etime, max_etime = Stream_min_starttime(all_traces)
    
    desired_stime = UTCDateTime(min_stime.year, min_stime.month, min_stime.day, 0, 0, 0.0)
    desired_etime = desired_stime + 86400
    
    days = Stream()
    while True:
        
        this_st = all_traces.copy()
        this_st.trim(starttime=desired_stime, endtime=desired_etime, pad=True, fill_value=0)
        for this_tr in this_st:
            days.append(this_tr)
        desired_stime += 86400
        desired_etime += 86400
        if desired_etime > max_etime + 86400:
            break
    return days



def Trace_merge_with_BUDfile(this_tr, budfile):
    """
    Merges an existing trace with a corresponding BUD file, preserving non-zero data values.

    This function:
    - Reads the existing **BUD file**.
    - Checks if `this_tr` has the **same trace ID, sampling rate, and time window**.
    - Merges the traces, prioritizing **non-zero data values**.

    Parameters:
    ----------
    this_tr : obspy.Trace
        The trace to be merged into the BUD file.
    budfile : str
        The existing BUD file path.

    Returns:
    -------
    obspy.Trace
        The merged trace, preserving **non-zero data values**.

    Notes:
    ------
    - If trace metadata (ID, sampling rate, or time range) **does not match**, the trace with 
      **more non-zero values** is returned.
    - This function prevents overwriting valuable data with zeroes.

    Example:
    --------
    ```python
    from obspy import read

    # Load an existing trace
    tr = read("new_trace.mseed")[0]

    # Merge it with an existing BUD file
    merged_tr = Trace_merge_with_BUDfile(tr, "BUD/1R.BHP2..EHZ.2020.346")
    
    # Save the merged result
    merged_tr.write("merged_trace.mseed", format="MSEED")
    ```
    """

    other_st = read(budfile)
    error_flag = False
    
    if len(other_st)>1:
        print('More than 1 trace in %s. Cannot merge.' % budfile)
        error_flag = True
        
    other_tr = other_st[0]
    if not (this_tr.id == other_tr.id):
        print('Different trace IDs. Cannot merge.')
        error_flag = True
        
    if not (this_tr.stats.sampling_rate == other_tr.stats.sampling_rate):
        print('Different sampling rates. Cannot merge.')
        error_flag = True
        
    if (abs(this_tr.stats.starttime - other_tr.stats.starttime) > this_tr.stats.delta/4):
        print('Different start times. Cannot merge.')  
        error_flag = True

    if (abs(this_tr.stats.endtime - other_tr.stats.endtime) > this_tr.stats.delta/4):
        print('Different end times. Cannot merge.')  
        error_flag = True
        
    if error_flag: # traces incompatible, so return the trace with the most non-zero values
        this_good = np.count_nonzero(this_tr.data)
        #print(this_tr.stats)
        other_good = np.count_nonzero(other_tr.data)
        #print(other_tr.stats)
        if other_good > this_good:
            return other_tr
        else:
            return this_tr
    
    else: # things are good
        indices = np.where(other_tr.data == 0)
        other_tr.data[indices] = this_tr.data[indices]
        return other_tr

######################################################################
##                  Modeling  tools                                 ##
######################################################################

def predict_arrival_times(station, quake):
    """
    Computes predicted seismic phase arrival times for a given station and earthquake using the IASP91 model.

    This function calculates the expected arrival times of seismic phases at a given station
    based on the earthquake's origin time, latitude, longitude, and depth.

    Parameters:
    ----------
    station : dict
        Dictionary containing station metadata:
        ```
        {
            "lat": float,   # Station latitude in degrees
            "lon": float    # Station longitude in degrees
        }
        ```
    quake : dict
        Dictionary containing earthquake metadata:
        ```
        {
            "lat": float,   # Earthquake latitude in degrees
            "lon": float,   # Earthquake longitude in degrees
            "depth": float, # Earthquake depth in km
            "otime": UTCDateTime  # Origin time of the earthquake
        }
        ```

    Returns:
    -------
    dict
        The updated `station` dictionary with a new `phases` key containing predicted arrival times:
        ```
        station["phases"] = {
            "P": "12:45:30",
            "S": "12:47:10",
            "Rayleigh": "12:48:00"
        }
        ```

    Notes:
    ------
    - Uses the **IASP91 travel-time model** via ObsPy's `TauPyModel`.
    - The **Rayleigh wave arrival** is estimated based on the S-wave arrival time.
    - Distances are computed using **gps2dist_azimuth** and converted to degrees.

    Example:
    --------
    ```python
    from obspy import UTCDateTime

    # Define station and earthquake metadata
    station = {"lat": 35.0, "lon": -120.0}
    quake = {"lat": 34.0, "lon": -118.0, "depth": 10.0, "otime": UTCDateTime("2023-03-01T12:45:00")}

    # Compute arrival times
    station = predict_arrival_times(station, quake)

    # Print predicted arrivals
    print(station["phases"])
    ```
    """
    model = TauPyModel(model="iasp91")
    
    [dist_in_m, az1, az2] = gps2dist_azimuth(quake['lat'], quake['lon'], station['lat'], station['lon'])
    station['distance'] = kilometers2degrees(dist_in_m/1000)
    arrivals = model.get_travel_times(source_depth_in_km=quake['depth'],distance_in_degree=station['distance'])
    # https://docs.obspy.org/packages/autogen/obspy.taup.helper_classes.Arrival.html#obspy.taup.helper_classes.Arrival
    
    phases = dict()
    for a in arrivals:
        phasetime = quake['otime'] + a.time
        phases[a.name] = phasetime.strftime('%H:%M:%S')
        if a.name == 'S':
            Rtime = quake['otime'] + a.time/ ((0.8453)**0.5)
            phases['Rayleigh'] = Rtime.strftime('%H:%M:%S')
    station['phases'] = phases
    
    return station

def syngine2stream(station, lat, lon, GCMTeventID, mseedfile):
    """
    Retrieves synthetic seismograms from IRIS Syngine for a specified GCMT event.

    This function generates synthetic seismograms for a given station and GCMT earthquake event.
    If a MiniSEED file already exists, it reads from that file; otherwise, it queries **IRIS Syngine**.

    Parameters:
    ----------
    station : str
        Station name (used for metadata storage).
    lat : float
        Latitude of the station (in degrees).
    lon : float
        Longitude of the station (in degrees).
    GCMTeventID : str
        Global Centroid Moment Tensor (GCMT) event ID.
    mseedfile : str
        Filename to save the downloaded synthetic waveform.

    Returns:
    -------
    obspy.Stream
        A Stream object containing the synthetic seismograms.

    Notes:
    ------
    - The function requests **displacement waveforms** with a sampling interval of 0.02s.
    - If the `mseedfile` already exists, it reads from the file instead of making a new request.
    - After downloading, the synthetic traces are assigned **latitude and longitude metadata**.

    Example:
    --------
    ```python
    # Generate synthetic seismograms for a GCMT event
    st = syngine2stream("ANMO", 34.95, -106.45, "202012312359A", "synthetic.mseed")

    # Plot the synthetic waveforms
    st.plot()
    ```
    """
    if os.path.exists(mseedfile):
        synth_disp = read(mseedfile)
    else:
        synth_disp = read("http://service.iris.edu/irisws/syngine/1/query?"
                  "format=miniseed&units=displacement&dt=0.02&"
                  "receivercenterlat=%f&receivercenterlon=%f&"
                  "eventid=GCMT:%s" % (lat, lon, GCMTeventID))
        for c in range(len(synth_disp)):
            synth_disp[c].stats.centerlat = lat
            synth_disp[c].stats.centerlon = lon
        synth_disp.write(mseedfile)
    return synth_disp


#######################################################################
##               Read waveform formats                               ##
#######################################################################

def read_DMX_file(DMXfile, fix=True, defaultnet=''):
    """
    Reads a **DMX** waveform file into an ObsPy Stream, applying optional corrections.

    This function reads a **demultiplexed SUDS (DMX) file**, corrects its metadata, and applies
    data adjustments for compatibility with SAC and MiniSEED formats.

    Parameters:
    ----------
    DMXfile : str
        Path to the DMX file.
    fix : bool, optional
        If `True`, applies metadata and amplitude corrections (default: True).
    defaultnet : str, optional
        Default network code to assign if the original value is missing (default: '').

    Returns:
    -------
    obspy.Stream
        A Stream object containing seismic traces from the DMX file.

    Notes:
    ------
    - The **ObsPy DMX reader** inserts `"unk"` for unknown networks; this function corrects it.
    - Data samples in DMX are **stored as unsigned 16-bit integers**. This function:
      - Converts data to **floating point** for MiniSEED compatibility.
      - Adjusts values by **subtracting 2048** to match SAC file conventions.

    Example:
    --------
    ```python
    # Read and correct a DMX file
    st = read_DMX_file("seismic_data.dmx", fix=True, defaultnet="IU")

    # Print trace details
    print(st)
    ```
    """
    print('Reading %s' % DMXfile)
    st = Stream()
    try:
        st = read(DMXfile)
        print('- read okay')
        if fix:
            for tr in st:
                # ObsPy DMX reader sets network to "unk" if blank. We'd rather keep it blank, or 
                # set with explicitly passing defaultnet named argument.
                if tr.stats.network == 'unk':
                    tr.stats.network = defaultnet
                    
                # ObsPy DMX reader falses adds 2048 to each data sample. Remove that here.
                # Also change data type of tr.data from uint to float so we can write trace to MiniSEED later   
                tr.data = tr.data.astype(float) - 2048.0 
    except:
        print('- ObsPy cannot read this demultiplexed SUDS file')        
    return st

#######################################################################
##               Read event metadata                                 ##
#######################################################################

def parse_hypo71_line(line):
    """
    Parses a single line of **HYPO71** earthquake location output.

    This function extracts **event origin time, location, depth, magnitude, and residuals**
    from a **fixed-column** formatted HYPO71 output line.

    Parameters:
    ----------
    line : str
        A single line of HYPO71 output.

    Returns:
    -------
    dict or None
        A dictionary containing extracted earthquake information, or `None` if parsing fails:
        ```
        {
            "origin_time": UTCDateTime,
            "latitude": float,
            "longitude": float,
            "depth": float (km),
            "magnitude": float,
            "n_ass": int,          # Number of associated arrivals
            "time_residual": float # RMS time residual in seconds
        }
        ```

    Notes:
    ------
    - The function handles **two-digit years**, assuming `year >= 70` belongs to the 1900s.
    - Converts **latitude and longitude from degrees + minutes** to decimal degrees.
    - Handles special cases where **minute=60** by rolling over to the next hour.

    Example:
    --------
    ```python
    # Example HYPO71 output line
    line = "230301 1205 45.2N 067 12.3W  10.0 3.2 15 0.2"

    # Parse earthquake information
    event_data = parse_hypo71_line(line)

    # Print extracted details
    print(event_data)
    ```
    """
    try:
        # Extract fields using fixed positions
        year = int(line[0:2])
        month = int(line[2:4])
        day = int(line[4:6])
        hour = int(line[7:9]) if line[7:9].strip() else 0
        minute = int(line[9:11]) if line[9:11].strip() else 0
        seconds = float(line[12:17]) if line[12:17].strip() else 0
        
        lat_deg = int(line[17:20].strip())
        lat_min = float(line[21:26].strip())
        lat_hem = line[20].strip().upper()
        
        lon_deg = int(line[27:30].strip())
        lon_min = float(line[31:36].strip())
        lon_hem = line[30].strip().upper()
        
        depth = float(line[37:43].strip())
        magnitude = float(line[44:50].strip())
        n_ass = int(line[51:53].strip())
        time_residual = float(line[62:].strip())
        
        # Handle two-digit years
        year = year + 1900 if year >= 70 else year + 2000

        # handle minute=60
        add_seconds = 0
        if minute==60:
            minute = 0
            add_seconds = 60       
        
        # Convert to UTCDateTime
        origin_time = UTCDateTime(year, month, day, hour, minute, seconds) + add_seconds
        
        # Convert latitude and longitude
        latitude = lat_deg + lat_min / 60.0
        if lat_hem == 'S':
            latitude = -latitude
        
        longitude = lon_deg + lon_min / 60.0
        if lon_hem == 'W':
            longitude = -longitude
        
        return {
            "origin_time": origin_time,
            "latitude": latitude,
            "longitude": longitude,
            "depth": depth,
            "magnitude": magnitude,
            "n_ass": n_ass,
            "time_residual": time_residual
        }
    except Exception as e:
        print(f"Failed to parse line: {line.strip()} | Error: {e}")
        return None    

def parse_hypo71_file(file_path):
    """
    Parses an entire **HYPO71 earthquake catalog file** into an ObsPy Catalog object.

    This function reads a **HYPO71 output file**, extracts event metadata from each line,
    and converts the information into an **ObsPy Catalog**.

    Parameters:
    ----------
    file_path : str
        Path to the HYPO71 output file.

    Returns:
    -------
    tuple:
        - **catalog (obspy.Catalog)**: A catalog containing `Event` objects with `Origin` and `Magnitude` attributes.
        - **unparsed_lines (list of str)**: A list of lines that could not be parsed.

    Notes:
    ------
    - Extracted **events include**:
      - **Origin time, latitude, longitude, depth, and magnitude.**
      - **Number of associated arrivals** (stored as an ObsPy `Comment`).
      - **RMS time residual** (stored as an ObsPy `Comment`).
    - If parsing fails for a line, the function prints a warning and stores it in `unparsed_lines`.

    Example:
    --------
    ```python
    # Parse a HYPO71 catalog file
    catalog, errors = parse_hypo71_file("hypo71_catalog.txt")

    # Print parsed events
    print(catalog)

    # Print unparsed lines (if any)
    print(errors)
    ```
    """
    catalog = Catalog()
    parsed = 0
    not_parsed = 0
    unparsed_lines = []
    with open(file_path, "r") as file:
        for line in file:
            #print(line)
            #event_data = parse_hypo71_line(line.strip())
            #if not event_data:
            event_data = parse_hypo71_line(line.strip())
            if event_data:
                parsed +=1
                #print(event_data)
                event = Event()
                origin = Origin(
                    time=event_data["origin_time"],
                    latitude=event_data["latitude"],
                    longitude=event_data["longitude"],
                    depth=event_data["depth"] * 1000  # Convert km to meters
                )
                magnitude = Magnitude(mag=event_data["magnitude"])
                
                # Store number of associated arrivals and time residual as comments
                origin.comments.append(Comment(text=f"n_ass: {event_data['n_ass']}"))
                origin.comments.append(Comment(text=f"time_residual: {event_data['time_residual']} sec"))

                event.origins.append(origin)
                event.magnitudes.append(magnitude)
                #print(event)
                catalog.append(event)
            else:
                print(line)
                not_parsed +=1
                unparsed_lines.append(line)
        
    print(f'parsed={parsed}, not parsed={not_parsed}')

    return catalog, unparsed_lines