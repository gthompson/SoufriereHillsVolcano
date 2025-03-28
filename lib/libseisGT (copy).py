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
##               Event detection                                     ##
#######################################################################



    








   

#not sure where this goes from ChapGpt










if __name__ == "__main__":

    # Create a sample trace with gaps
    data = np.concatenate([np.random.randn(100), np.full(5, np.nan), np.random.randn(100)])
    tr = Trace(data=data)

    # Run the detrending function
    tr = detrend_trace(tr, gap_threshold=3, verbose=True)


    # Load a seismic trace with gaps
    st = read("example_with_gaps.mseed")
    tr = st[0]  # Extract first trace

    # Process gaps
    tr = detect_and_handle_gaps(tr, gap_threshold=5, verbose=True)

    # Plot the result
    tr.plot()         
