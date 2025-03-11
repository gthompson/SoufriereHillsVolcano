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

def add_to_trace_history(tr, str):
    if not 'history' in tr.stats:
        tr.stats['history'] = list()
    if not str in tr.stats.history:
        tr.stats.history.append(str)

def update_trace_filter(tr, filtertype, freq, zerophase):
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

def clip_trace(tr, AMP_LIMIT = 1e10, fill_value=None):     
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

    
def smart_merge_traces(trace_pair):
    """
    Clever way to merge overlapping traces. Uses all non-zero data values from both.
    """
    this_tr = trace_pair[0] 
    other_tr = trace_pair[1]

    error_flag = False


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
        
def pad_trace(tr, seconds):
    if seconds>0.0:
        y = tr.data
        tr.stats['originalStartTime'] = tr.stats.starttime
        tr.stats['originalEndTime'] = tr.stats.endtime
        npts_pad = int(tr.stats.sampling_rate * seconds)
        y_prepend = np.flip(y[0:npts_pad])
        y_postpend = np.flip(y[-npts_pad:])
        y = np.concatenate( [y_prepend, y, y_postpend ] )
        padStartTime = tr.stats.starttime - npts_pad * tr.stats.delta
        tr.data = y
        add_to_trace_history(tr, 'padded')
        tr.stats.starttime = padStartTime

def unpad_trace(tr):
    s = tr.stats
    if 'originalStartTime' in s:
        tr.trim(starttime=s.originalStartTime, endtime=s.originalEndTime, pad=False)
        add_to_trace_history(tr, 'unpadded') 

def remove_single_sample_spikes(trace, threshold=10):
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

'''
def clean_trace(tr, taperFraction=0.05, filterType="bandpass", freq=[0.1, 20.0], \
                corners=2, zerophase=True, inv=None, outputType='VEL'):
    """
    Clean Trace object in place.
    clean_trace(tr, taperFraction=0.05, filterType="bandpass", freq=[0.1, 20.0], corners=2, zerophase=True, inv=None)
    """

    print('clean_trace is Deprecated. using simple_clean instead')

    if not 'history' in tr.stats:
        tr.stats['history'] = list()    
    
    # remove absurd values
    print('- clipping')
    clip_trace(tr) # could add function here to correct for clipping - algorithms exist

    # remove single sample spikes
    print('- removing single sample spikes')
    remove_single_sample_spikes(tr)
    
    # save the start and end times for later 
    startTime = tr.stats.starttime
    endTime = tr.stats.endtime
    #print(tr.id, startTime, endTime)
    
    # get trace max and min
    #amp_before = (max(tr.data)-min(tr.data))/2
    #stdev_before = np.std(tr.data)
        
    # pad the Trace
    #y = tr.data
    npts = tr.stats.npts
    npts_pad = int(taperFraction * npts)
    npts_pad_seconds = npts_pad * tr.stats.delta
    if npts_pad_seconds < 1/freq[0]: # impose a minimum pad length of 10-seconds
        npts_pad_seconds = 1/freq[0]
    
    print('- padding')
    pad_trace(tr, npts_pad_seconds)
    max_fraction = npts_pad / tr.stats.npts
    
    # clean
    if not 'detrended' in tr.stats.history:
        try:
            print('- detrending')
            tr.detrend('linear')
            
        except:
            try:
                tr.data = tr.data - np.nanmedian(tr.data)
                if np.ma.is_masked(tr.data):
                    # Replace masked values with 0
                    tr.data = tr.data.filled(0)   
            except:
                tr.detrend('simple')
        add_to_trace_history(tr, 'detrended')
    
        
    if not 'tapered' in tr.stats.history:
        print('- tapering')
        tr.taper(max_percentage=max_fraction, type="hann") 
        add_to_trace_history(tr, 'tapered')        
    
    print('- filtering')
    if filterType == 'bandpass':
        tr.filter(filterType, freqmin=freq[0], freqmax=freq[1], corners=corners, zerophase=zerophase)
    else:    
        tr.filter(filterType, freq=freq, corners=corners, zerophase=zerophase)
    update_trace_filter(tr, filterType, freq, zerophase)
    add_to_trace_history(tr, filterType)    
        

    '
    if not 'deconvolved' in tr.stats.history:
        if not 'units' in tr.stats:
            tr.stats['units'] = 'Counts'   
    if inv and tr.stats['units'] == 'Counts':
        tr.remove_response(inventory=inv, output="VEL") 
        tr.stats['units'] = 'm/s'
        add_to_trace_history(tr, 'deconvolved')
    if not inv and tr.stats.calib and not 'calibrated' in tr.stats.history:
        if not tr.stats.calib==1.0:
            tr.data = tr.data * tr.stats.calib
            tr.stats['units'] = 'm/s'
            add_to_trace_history(tr, 'calibrated')
          
    '
    
    if not 'units' in tr.stats:
        tr.stats['units'] = 'Counts'   
        
    if tr.stats['units'] == 'Counts' and not 'calibrated' in tr.stats.history:
        if inv:
            try:
                print('- calling removeInstrumentResponse')
                removeInstrumentResponse(tr, None, outputType = outputType, inventory = inv, taperFraction=0)
                #tr.remove_response(inventory=inv)
            except:
                print('No matching response info found for %s' % tr.id)
            else:
                add_to_trace_history(tr, 'calibrated')
                # update the calib value
                tr.stats.calib = _get_calib(tr, inv)
                
        elif not tr.stats.calib==1.0:
            tr.data = tr.data * tr.stats.calib
            add_to_trace_history(tr, 'calibrated') 
        if 'calibrated' in tr.stats.history:           
            if tr.stats.channel[1]=='H':
                tr.stats['units'] = 'm/s'
            if tr.stats.channel[1]=='N':
                tr.stats['units'] = 'm/s2'                
            if tr.stats.channel[1]=='D':
                tr.stats['units'] = 'Pa'  
              
    # remove the pad
    #tr.trim(starttime=startTime, endtime=endTime, pad=False)
    #add_to_trace_history(tr, 'unpadded')
    print('- unpadding')
    unpad_trace(tr)
    
    #amp_after = (max(tr.data)-min(tr.data))/2
    #stdev_after = np.std(tr.data)
    #tr.stats.calib = amp_after / amp_before
    #tr.stats.calib = stdev_before / stdev_after
    #if 'calibrated' in tr.stats.history:
    #    tr.stats.calib = amp_after / amp_before
'''

def simple_clean(tr, taperFraction=0.05, filterType="bandpass", freq=[0.1, 20.0], \
                corners=2, zerophase=True, inv=None, outputType='VEL', verbose=False):


    if not 'history' in tr.stats:
        tr.stats['history'] = list()   

    if not 'units' in tr.stats:
        tr.stats['units'] = 'Counts'           
    
    # remove absurd values
    if verbose:
        print('- clipping')
    clip_trace(tr) # could add function here to correct for clipping - algorithms exist

    # remove single sample spikes
    if verbose:
        print('- removing single sample spikes')
    remove_single_sample_spikes(tr)
    
    # save the start and end times for later 
    startTime = tr.stats.starttime
    endTime = tr.stats.endtime
        
    # pad the Trace
    npts = tr.stats.npts
    npts_pad = int(taperFraction * npts)
    npts_pad_seconds = npts_pad * tr.stats.delta
    if npts_pad_seconds < 1/freq[0]: # impose a minimum pad length of 10-seconds
        npts_pad_seconds = 1/freq[0]
    
    if verbose:
        print('- padding')
    pad_trace(tr, npts_pad_seconds)
    max_fraction = npts_pad / tr.stats.npts

    # filter or fully correct
    if inv:
        # fully correct
        if verbose:
            print('- removing instrument response')
        #removeInstrumentResponse(tr, None, outputType = outputType, inventory = inv, taperFraction=0)
        tr.remove_response(inventory=inv, output=outputType, \
                        pre_filt=(freq[0]/1.5, freq[0], freq[1], freq[1]*1.5), \
                        water_level=60, zero_mean=True, \
                        taper=False, taper_fraction=taperFraction, plot=False, fig=None)
        add_to_trace_history(tr, 'calibrated')
        #tr.stats.calib = _get_calib(tr, inv)        
    else:
        # detrend, taper, filter
        if not 'detrended' in tr.stats.history:
            try:
                if verbose:
                    print('- detrending')
                tr.detrend('linear')
                
            except:
                try:
                    tr.data = tr.data - np.nanmedian(tr.data)
                    if np.ma.is_masked(tr.data):
                        # Replace masked values with 0
                        tr.data = tr.data.filled(0)   
                except:
                    tr.detrend('simple')
            add_to_trace_history(tr, 'detrended')
     
        if not 'tapered' in tr.stats.history:
            if verbose:
                print('- tapering')
            tr.taper(max_percentage=max_fraction, type="hann") 
            add_to_trace_history(tr, 'tapered')        
    
        if verbose:
            print('- filtering')
        if filterType == 'bandpass':
            tr.filter(filterType, freqmin=freq[0], freqmax=freq[1], corners=corners, zerophase=zerophase)
        else:    
            tr.filter(filterType, freq=freq, corners=corners, zerophase=zerophase)
        update_trace_filter(tr, filterType, freq, zerophase)
        add_to_trace_history(tr, filterType)  

    # remove the pad
    if verbose:
        print('- unpadding')
    unpad_trace(tr)          

    '''
    if tr.stats['units'] == 'Counts' and not 'calibrated' in tr.stats.history and tr.stats.calib!=1.0:
        tr.data = tr.data * tr.stats.calib
        add_to_trace_history(tr, 'calibrated') 
    '''

    if 'calibrated' in tr.stats.history:           
        if tr.stats.channel[1]=='H':
            if outputType=='VEL':
                tr.stats['units'] = 'm/s'
            elif outputType=='DISP':
                tr.stats['units'] = 'm'
        '''
        if tr.stats.channel[1]=='N':
            tr.stats['units'] = 'm/s2'                
        if tr.stats.channel[1]=='D':
            tr.stats['units'] = 'Pa'  
        '''


                        
def _get_calib(tr, this_inv):
    calib = 1.0
    for station in this_inv.networks[0].stations:
        if station.code == tr.stats.station:
            for channel in station.channels:
                if channel.code == tr.stats.channel:
                    calib_freq, calib_value = channel.response._get_overall_sensitivity_and_gain()
    return calib_value
        
  ################### New trace tools added in Feb 2025: Start ####################
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

def _get_band_code(sampling_rate):
    """Determine the appropriate band code based on sampling rate."""
    for (low, high), code in BAND_CODE_TABLE.items():
        if low <= sampling_rate < high:
            return code
    return None  # Should not happen if lookup table is correct

def _adjust_band_code_for_sensor_type(current_band_code, expected_band_code, short_period=False):
    """
    Adjusts the expected band code if the current band code indicates a short-period seismometer.

    Short-period seismometers use different band codes compared to broadband seismometers.
    If the current_band_code is one of 'S', 'E', 'C', or 'F', then:
      - 'B' (Broadband) -> 'S' (Short-period)
      - 'H' (High-frequency broadband) -> 'E' (Short-period high-frequency)
      - 'D' (Very long period broadband) -> 'C' (Short-period very long period)
      - 'G' (Extremely high frequency broadband) -> 'F' (Short-period extremely high frequency)

    :param current_band_code: The first letter of the current trace.stats.channel (e.g., 'S', 'E', 'C', 'F')
    :param expected_band_code: The computed band code based on sampling rate
    :return: Adjusted band code if necessary
    """
    short_period_codes = {'S', 'E', 'C', 'F'}
    
    if current_band_code in short_period_codes or short_period:
        band_code_mapping = {'B': 'S', 'H': 'E', 'D': 'C', 'G': 'F'}
        return band_code_mapping.get(expected_band_code, expected_band_code)
    
    return expected_band_code

def _fix_legacy_id(trace):
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
    # note: for MVO data call libMVO.fix_trace_id_mvo
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

def check_write_read(tr):
    mfile = '/tmp/tmpminiseedfile'
    try:
        tr.write(mfile, format='MSEED')
        tr2 = read(mfile, format='MSEED')
        #tr2.plot();
        #plt.show()
        os.remove(mfile)
        return True
    except Exception as e:
        print(e)
        return False

def process(st_or_tr, func):
    """Recursively calls a given function on all traces in a Stream, or a single Trace if passed."""
    if isinstance(st_or_tr, Stream):
        for trace in st_or_tr:
            process(trace, func)  # Recursive call for each Trace
    elif isinstance(st_or_tr, Trace):
        func(st_or_tr)  # Base case: process a single Trace
    else:
        raise TypeError("Input must be an ObsPy Stream or Trace object.")

################### New trace tools added in Feb 2025: End ####################

# QC tools moved from metrics library:


"""
Functions for computing data quality metrics and statistical metrics (such as amplitude, energy and frequency) 
on Stream/Trace objects.

In terms of order of application:

1. Read the raw data.
2. Fix Trace IDs.
3. process_trace


"""

def process_trace(tr, bool_despike=True, bool_clean=True, inv=None, quality_threshold=0.0, taperFraction=0.05, \
                  filterType="bandpass", freq=[0.5, 30.0], corners=6, zerophase=False, outputType='VEL', \
                    miniseed_qc=True, verbose=False, max_dropout=0.0):
    if verbose:
        print(f'Processing {tr}')
    if not 'history' in tr.stats:
        tr.stats['history'] = list()    
        
    """ RAW DATA QC METRICS """
    if miniseed_qc:
        try:
            qcTrace(tr)
        except:
            print('qcTrace failed on %s for raw trace' % tr.id)
            tr.stats['quality_factor'] = -1
        else:        
            tr.stats["quality_factor"] = trace_quality_factor(tr, max_dropout=max_dropout) #0 = blank trace, 1 = has some 0s and -1s, 3 = all looks good
            tr.stats.quality_factor -= tr.stats.metrics['num_gaps']
            tr.stats.quality_factor -= tr.stats.metrics['num_overlaps']
            tr.stats.quality_factor *= tr.stats.metrics['percent_availability']/100.0
            tr.stats.metrics["twin"] = tr.stats.npts /  tr.stats.sampling_rate # before or after detrending  
    else:
        tr.stats['quality_factor'] = 3
    
    if tr.stats.quality_factor > quality_threshold: # only clean traces better than the threshold

        # Check for spikes - been seeing this in 1998 data
        if bool_despike:
            check_for_spikes(tr)
    elif verbose:
        print(f'- not despiking. qf={tr.stats.quality_factor}')

    if tr.stats.quality_factor > quality_threshold:

        """ CLEAN (DETREND, BANDPASS, CORRECT) TRACE """
        if bool_clean:
            try:
                #clean_trace(tr, taperFraction=taperFraction, filterType=filterType, freq=freq, corners=corners, zerophase=zerophase, inv=inv, outputType=outputType)
                simple_clean(tr, taperFraction=taperFraction, filterType=filterType, freq=freq, corners=corners, zerophase=zerophase, inv=inv, outputType=outputType)
                #tr.remove_response(inventory=inv, output=outputType, pre_filt=(freq[0]/2, freq[0], freq[1], freq[1]*2), water_level=60, zero_mean=True, taper=True, taper_fraction=taperFraction, plot=False, fig=None)
            except Exception as e:
                print(e)
                return False
        ''' 
        The logic here is probably that clean trace could change the zeros in the trace
        But is that what we want?
        # Update other stats
        if miniseed_qc:
            try:
                qcTrace(tr) # swap this out for miniseed write/read
                check_write_read(tr)
            except:
                print('qcTrace failed on %s for cleaned trace' % tr.id)
        if not check_write_read(tr):
            # fi cannot write to miniseed and read back again, lower the quality_factor?
        '''
    elif verbose:
        print(f'- not cleaning. qf={tr.stats.quality_factor}')
    return True

def check_for_spikes(tr):
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
    
def qcTrace(tr):
    """ qcTrace(tr) DATA QUALITY CHECKS """
    
    """ Useful MSEED metrics
    {'start_gap': None, 'end_gap': None, 'num_gaps': 0, 
     'sum_gaps': 0, 'max_gap': None, 'num_overlaps': 0, 
     'sum_overlaps': 0, 'max_overlap': None, 'quality': 'D', 
     'sample_min': -22404, 'sample_max': 9261, 
     'sample_mean': -3854.7406382978725, 'sample_median': -3836.0, 
     'sample_lower_quartile': -4526.0, 'sample_upper_quartile': -3105.0, 
     'sample_rms': 4426.1431329789848, 
     'sample_stdev': 2175.2511682727431, 
     'percent_availability': 100.0}           
    """
    if len(tr.data)>0:
        tmpfilename = '%s%s.mseed' % (tr.id, tr.stats.starttime.isoformat())
        tr.write(tmpfilename)
        mseedqc = MSEEDMetadata([tmpfilename]) 
        tr.stats['metrics'] = mseedqc.meta
        os.remove(tmpfilename)
        add_to_trace_history(tr, 'MSEED metrics computed (similar to ISPAQ/MUSTANG).')
    else:
        tr.stats['quality_factor'] = -100

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

    
def _get_islands(arr, mask):
    mask_ = np.concatenate(( [False], mask, [False] ))
    idx = np.flatnonzero(mask_ [1:] != mask_ [:-1])
    return [arr[idx[i]:idx[i+1] + 1] for i in range(0, len(idx), 2)]

def _FindMaxLength(lst):
    maxList = max(lst, key = len)
    maxLength = max(map(len, lst))      
    return maxList, maxLength  

def trace_quality_factor(tr, min_sampling_rate=19.99, max_dropout=0.0):
    # trace_quality_factor(tr)
    # a good trace has quality factor 3, one with 0s and -1s has 1, bad trace has 0
    quality_factor = 1.0
    is_bad_trace = False
    
    # ignore traces with few samples
    if tr.stats.npts < 100:
        add_to_trace_history(tr, 'Not enough samples')
        is_bad_trace = True
    
    # ignore traces with weirdly low sampling rates
    if tr.stats.sampling_rate < min_sampling_rate:
        add_to_trace_history(tr, 'Sampling rate too low')
        is_bad_trace = True

    # ignore blank trace
    anyData = np.count_nonzero(tr.data)
    if anyData==0:
        add_to_trace_history(tr, 'Trace is blank')
        is_bad_trace = True
    
    # check for bit level noise
    u = np.unique(tr.data)
    num_unique_values = u.size
    if num_unique_values > 10:
        quality_factor += np.log10(num_unique_values)
    else:
        add_to_trace_history(tr, 'bit level noise suspected')
        is_bad_trace = True

    # check for sequences of 0 or 1
    if max_dropout == 0.0:
        trace_good_flag = _check0andMinus1(tr.data)
        if not trace_good_flag:
            add_to_trace_history(tr, 'sequences of 0 or -1 found')
            is_bad_trace = True
    
    # replacement for check0andMinus1
    seq = tr.data
    islands = _get_islands(seq, np.r_[np.diff(seq) == 0, False]) 
    try:
        maxList, maxLength = _FindMaxLength(islands)
        add_to_trace_history(tr, 'longest flat sequence found: %d samples' % maxLength)
        if maxLength >= tr.stats.sampling_rate * max_dropout:
            is_bad_trace = True
        else:
            fill_all_gaps(tr)
    except:
        is_bad_trace = True
        
    # time to exit?
    if is_bad_trace:
        return 0.0
        
    # check if trace clipped - but so far I don't see clipped trace as terminal
    upperClipped, lowerClipped = _detectClipping(tr) # can add another function to interpolate clipped values
    if upperClipped:
        quality_factor /= 2.0
    if lowerClipped:
        quality_factor /= 2.0
         
    # check for outliers - but I haven't tuned this yet, so not making it a decision between quality 0.0 and continue
    outlier_count, outlier_indices = _mad_based_outlier(tr, thresh=50.0)
    #print('Outliers: %d' % outlier_count)
    if outlier_count == 0:
        quality_factor += 1.0
    else:
        add_to_trace_history(tr, '%d outliers found' % outlier_count)
        tr.stats['outlier_indices'] = outlier_indices    
       
    return quality_factor    

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


def _check0andMinus1(liste):
# function bool_good_trace = check0andMinus1(tr.data)
    liste=list(liste)
    listStr=''.join(str(i) for i in liste)
    if  "000000000000" in listStr or "-1-1-1-1-1-1-1-1" in listStr :
        return False
    else:
        return True  

#######################################################################
###                        Gap filling tools                        ###
#######################################################################


def fill_all_gaps(trace, verbose=False):
    """
    Finds and fills all gaps in a seismic trace using the best method for each case.
    
    Parameters:
    trace (obspy.Trace): The seismic trace to process.
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
            linear_interpolation(trace, gap_start, gap_end)
        elif samples <= trace.stats.sampling_rate:  # Short gaps (≤ 1 sec) → Repeat data
            if verbose:
                print(f'gap start={gap_start}, end={gap_end}, repeating data')            
            repeat_previous_data(trace, gap_start, gap_end)
        else:  # Longer gaps → Fill with spectrally-matched noise
            if verbose:
                print(f'gap start={gap_start}, end={gap_end}, adding spectral noise')            
            fill_gap_with_filtered_noise(trace, gap_start, gap_end)
    if trace.stats.station == 'MBSS':
        trace.plot(outfile='MBSS_after_gap_filling.png')
def repeat_previous_data(trace, gap_start, gap_end):
    """
    Fills gaps by repeating the last valid segment before the gap.
    """
    sample_rate = trace.stats.sampling_rate
    num_samples = int((gap_end - gap_start) * sample_rate)

    # Get the last valid segment before the gap
    fill_data = trace.data[-num_samples:]  # Last 'num_samples' of valid data
    trace.data = np.concatenate((trace.data, fill_data))

def linear_interpolation(trace, gap_start, gap_end):
    """
    Fills gaps in the seismic trace using linear interpolation.
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

def fill_gap_with_filtered_noise(trace, gap_start, gap_end):
    """
    Fills a gap with noise that matches the spectral characteristics of the original signal
    using ObsPy's built-in filtering.
    """
    sample_rate = trace.stats.sampling_rate
    num_samples = int((gap_end - gap_start) * sample_rate)

    # Get last valid segment before the gap
    prev_segment = trace.data[-num_samples:]

    # Generate white noise
    noise = np.random.normal(scale=np.std(prev_segment), size=num_samples)

    # Convert noise to an ObsPy trace and apply the same filtering as the original trace
    noise_trace = Trace(data=noise, header=trace.stats)
    noise_trace.filter("bandpass", freqmin=0.01, freqmax=0.5, corners=4, zerophase=True)

    # Append the filtered noise to the trace
    trace.data = np.concatenate((trace.data, noise_trace.data))

    
#######################################################################
##                Stream tools                                       ##
#######################################################################


def remove_empty_traces(stream):
    """Removes empty traces, traces full of zeros, and traces full of NaNs from an ObsPy Stream."""
    cleaned_stream = Stream()  # Create a new empty Stream

    for trace in stream:
        # Check if trace is empty (npts == 0)
        if trace.stats.npts == 0:
            continue
        
        # Check for flat trace (e.g. all zero, or all -1)
        if np.all(trace.data == np.nanmean(trace.data)):
            continue

        # Check if all values are NaN
        if np.all(np.isnan(trace.data)):
            continue

        # If trace passes all checks, add it to the cleaned stream
        cleaned_stream += trace

    return cleaned_stream    

def remove_low_quality_traces(st, quality_threshold=1.0):
    for tr in st:
        if tr.stats.quality_factor < quality_threshold: 
            st.remove(tr)

def piecewise_detrend(st, null_value=0, fill_value=np.nan, detrend='linear', verbose=False): 
    # takes a Stream object nominally from an SDSclient that has gaps marked by zeros, and applies

    if not detrend and not highpass:
        return

    # split into contiguous segments
    isSplit = False
    for tr in st:
        if np.any(tr.data == null_value):
            tr.data = np.ma.masked_where(tr.data == null_value, tr.data)
            all_traces = tr.split()
            if verbose:
                print(f'{tr.id} split into {len(all_traces)} traces')
            isSplit=True
            all_traces.detrend(detrend)
            all_traces.merge(method=0, fill_value=fill_value)
            print(f'merged again. stream now contains {len(all_traces)} traces')
            if len(all_traces) == 1:
                tr = all_traces.copy()
            else:
                st.remove(tr)
                for newtr in all_traces:
                    st.append(newtr)

    # detrend
    #if detrend:
    #    st.detrend(detrend)        
    
    # recombine
    #if isSplit:
    #    st.merge(method=0, fill_value=fill_value)
    #    if verbose:
    #        print('merged again. stream now contains {len(st)} traces')


def get_seed_band_code(sr, shortperiod = False):
    # adjusted Feb 2025 to exploit new code above
    
    # Determine the correct band code
    bc = _get_band_code(sr) # this assumes broadband sensor

    # adjust if short-period sensor
    bc = _adjust_band_code_for_sensor_type(None, bc, shortperiod)  
    '''
    bc = '_'
    if sr >= 1 and sr < 10:
        bc = 'M'
    if sr >= 10 and sr < 80:
        if shortperiod:
            bc = 'S'
        else:
            bc = 'B'
    if sr >= 80:
        if shortperiod:
            bc = 'E'
        else:
            bc = 'H'
    '''
    return bc
        
def fix_seed_band_code(st, shortperiod = False):  
    print('libseisGT.fix_seed_band_code is deprecated. please use libseisGT.fix_trace_id')
    for tr in st:
        sr = tr.stats.sampling_rate
        bc = get_seed_band_code(sr, shortperiod=shortperiod)
        if not bc == tr.stats.channel[0]:
            tr.stats.channel = bc + tr.stats.channel[1:]
            add_to_trace_history(tr, 'bandcode_fixed') 

                  
def smart_merge(st, verbose=False, interactive=False):
    print('also see smart_merge_traces(trace_pair)')
    # need to loop over st and find traces with same ids
    ##### GOT HERE
    newst = Stream()
    all_ids = []
    for tr in st:
        if not tr.id in all_ids:
            all_ids.append(tr.id)

    for this_id in all_ids: # loop over all nsl combinations
        these_traces = st.copy().select(id=this_id).sort() # find all traces for this nsl combination
        
        # remove duplicates
        traces_to_remove = []
        for c in range(len(these_traces)-1):
            s0 = these_traces[c].stats
            s1 = these_traces[c].stats
            if s0.starttime == s1.starttime and s0.endtime == s1.endtime and s0.sampling_rate == s1.sampling_rate:
                traces_to_remove.append(c)
        if verbose:
            print(these_traces)
            print(traces_to_remove)
        if traces_to_remove:
            for c in traces_to_remove:
                these_traces.remove(these_traces[c])
                
        if len(these_traces)==1: # if only 1 trace, append it, and go to next trace id
            newst.append(these_traces[0]) 
            continue
        
        # must have more than 1 trace
        try: # try regular merge now duplicates removed
            merged_trace = these_traces.copy().merge()
            if verbose:
                print('- regular merge of these traces success')
        except:
            if verbose:
                print('- regular merge of these traces failed')   
            # need to try merging traces in pairs instead
            N = len(these_traces)
            these_traces.sort() # sort the traces
            for c in range(N-1): # loop over traces in pairs
                appended = False
                
                # choose pair
                if c==0:
                    trace_pair = these_traces[0:2]
                else:
                    trace_pair = Stream(traces=[merged_trace, these_traces[c+1] ] )
                                        
                # merge these two traces together    
                try: # standard merge
                    merged_trace = trace_pair.copy().merge()
                    if verbose:
                        print('- regular merge of trace pair success')
                except: # smart merge
                    if verbose:
                        print('- regular merge of trace pair failed')
                    try:
                        min_stime, max_stime, min_etime, max_etime = ls.Stream_min_starttime(trace_pair)
                        trace_pair.trim(starttime=min_stime, endtime=max_etime, pad=True, fill_value=0)
                        merged_trace = Stream.append(smart_merge_traces(trace_pair)) # this is a trace, not a Stream
                    except:
                        print('- smart_merge of trace pair failed')
                        
            # we have looped over all pairs and merged_trace should now contain everything
            # we should only have 1 trace in merged_trace
            if verbose:
                print(merged_trace)
            if len(merged_trace)==1:
                try:
                    newst.append(merged_trace[0])
                    appended = True
                except:
                    pass
                
            if not appended:
                if interactive:
                    print('\n\nTrace conflict\n')
                    trace_pair.plot()
                    for c in range(len(trace_pair)):
                        print(c, trace_pair[c])
                    choice = int(input('Keep which trace ? '))
                    newst.append(trace_pair[choice])  
                    appended = True 
                else:
                    raise('not able to merge')

                             
    return newst 
        
        
def Stream_min_starttime(all_traces):
    """
    Take a Stream object, and return the minimum starttime

    Created for CALIPSO data archive from Alan Linde.
    """ 

    min_stime = UTCDateTime(2099, 12, 31, 0, 0, 0.0)
    max_stime = UTCDateTime(1900, 1, 1, 0, 0, 0.0)
    min_etime = UTCDateTime(2099, 12, 31, 0, 0, 0.0)
    max_etime = UTCDateTime(1900, 1, 1, 0, 0, 0.0)    
    for this_tr in all_traces:
        if this_tr.stats.starttime < min_stime:
            min_stime = this_tr.stats.starttime
        if this_tr.stats.starttime > max_stime:
            max_stime = this_tr.stats.starttime  
        if this_tr.stats.endtime < min_etime:
            min_etime = this_tr.stats.endtime
        if this_tr.stats.endtime > max_etime:
            max_etime = this_tr.stats.endtime              
    return min_stime, max_stime, min_etime, max_etime

'''
def removeInstrumentResponse(st_or_tr, preFilter = (1, 1.5, 30.0, 45.0), outputType = "VEL", inventory = None, taperFraction=0.05):  
    """
    Remove instrument response - note inventories may have been added to Stream object
    Written for Miami Lakes
    
    This function may be obsolete. Could improve clean_trace instead.

    Note that taper_fraction=0.05, water_level=60 are other default parameters
    """
    print('removeInstrumentResponse is Deprecated. For some reason this does not seem to work well. Use process_trace instead')
    if isinstance(st_or_tr, Stream):
        for trace in st_or_tr:
            try:
                removeInstrumentResponse(trace, preFilter = preFilter, outputType = outputType, inventory = inventory, taperFraction=taperFraction)
            except:
                st_or_tr.remove(trace)
            
    elif isinstance(st_or_tr, Trace):  
        tr = st_or_tr 
        bool_taper = False
        if taperFraction:
            bool_taper = True
        if has_response(inventory, tr):
            #print(f"Response exists for {tr.id}")
            print(f'Correcting {tr.id}')        
            tr.remove_response(output=outputType, pre_filt=preFilter, inventory = None, water_level=60, zero_mean=True, taper=bool_taper, taper_fraction=taperFraction, plot=False, fig=None)
        else:
            print(f"No response found for {tr.id}")         
    #OLDER STUFF
    try:
        st.remove_response(output=outputType, pre_filt=preFilter, inventory = None, water_level=60, zero_mean=True, taper=bool_taper, taper_fraction=taperFraction, plot=False, fig=None)
    except:
        for tr in st:
            try:
                tr.remove_response(output=outputType, pre_filt=preFilter, inventory=inventory)
            except:
                print("- Not able to correct data for %s " %  tr.id)
                st.remove(tr)
    return
    '''

def detect_network_event(st_in, minchans=None, threshon=3.5, threshoff=1.0, \
                         sta=0.5, lta=5.0, pad=0.0, best_only=False, verbose=False, freq=None, algorithm='recstalta', criterion='longest'):
    """
    Run a full network event detector/associator 
    
    Note that if you run a 5-s LTA, you need at least 5-s of noise before the signal.
    
    Output is a list of dicts like:
    
    {'cft_peak_wmean': 19.561900329259956,
 'cft_peaks': [19.535644192544272,
               19.872432918501264,
               19.622171410201297,
               19.217352795792998],
 'cft_std_wmean': 5.4565629691954713,
 'cft_stds': [5.292458320417178,
              5.6565387957966404,
              5.7582248973698507,
              5.1190298631982163],
 'coincidence_sum': 4.0,
 'duration': 4.5299999713897705,
 'stations': ['UH3', 'UH2', 'UH1', 'UH4'],
 'time': UTCDateTime(2010, 5, 27, 16, 24, 33, 190000),
 'trace_ids': ['BW.UH3..SHZ', 'BW.UH2..SHZ', 'BW.UH1..SHZ', 'BW.UH4..SHZ']}
 
 
    Any trimming of the Stream object can then by done with trim_to_event.
 
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
    Runs a single channel detection on each Trace. No coincidence trigger/association into event.
        
    Take a Stream object and run an STA/LTA on each channel, adding a triggers list to Trace.stats.
    This should be a list of 2-element numpy arrays with trigger times on and off as UTCDateTime
    
    Note that if you run a 5-s LTA, you need at least 5-s of noise before the signal.

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
    """ 
    Take a Stream object and run an STA/LTA on each channel. Return first trigger on and last trigger off time.
    Assumes that tr.stats has triggers lists added already with add_channel_detections
    Any trimming of the Stream taking into account pretrig and posttrig seconds can by done with trim_to_event
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
    """ Trims a Stream based on mintime and maxtime which could come from detect_network_event or get_event_window """
    st.trim(starttime=mintime-pretrig, endtime=maxtime+posttrig)
    
#######################################################################    
########################         WFDISC tools                        ##
#######################################################################


     
def index_waveformfiles(wffiles, ampeng=False, events=True):
    """ 
    Take a list of seismic waveform data files and return a dataframe similar to a wfdisc table

    Created for CALIPSO data archive from Alan Linde.
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
    Read a wfdisc-like dataframe, read the associated files, and write out as BUD format


    Created for CALIPSO data archive from Alan Linde.
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
    Process a directory containing waveform data files in any format readable by ObsPy.
    Build a wfdisc-like dataframe indexing those waveform files.
    Convert them to a BUD archive.

    Created for CALIPSO data archive from Alan Linde.
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
    """ 
    Take a Stream object and write it out in IRIS/PASSCAL BUD format. 
    
    Example:
    
        Stream_to_BUD('RAW', all_traces)
    Where all_traces is a Stream object with traces from 2020346

    Creates a BUD directory structure that looks like:

        DAYS
        ├── BHP2
        │   ├── 1R.BHP2..EH1.2020.346
        │   ├── 1R.BHP2..EH2.2020.346
        │   └── 1R.BHP2..EHZ.2020.346
        ├── BHP4
        │   ├── 1R.BHP4..EH1.2020.346
        │   ├── 1R.BHP4..EH2.2020.346
        │   └── 1R.BHP4..EHZ.2020.346
        ├── FIREP
        │   ├── 1R.FIREP..EH1.2020.346
        │   ├── 1R.FIREP..EH2.2020.346
        │   └── 1R.FIREP..EHZ.2020.346
        └── TANKP
            ├── 1R.TANKP..EH1.2020.346
            ├── 1R.TANKP..EH2.2020.346
            └── 1R.TANKP..EHZ.2020.346
        
    where BHP2, BHP4, FIREP and TANKP are station names, 1R is network name, 
    location is blank, channels are EH[Z12], year is 2020 and day of year is 346.

    Created for ROCKETSEIS data conversion and modified for CALIPSO data archive from Alan Linde.   
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
    Load all files corresponding to this year and day from a BUD archive


    Created for CALIPSO data archive from Alan Linde.
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



def Stream_to_dayplot(TOPDIR, all_traces):
    """ 
    Take a Stream object, pad it to 24-hours, plot it, and save to a PNG file. 
    
    Example: 
        Stream_to_dayplot('RAW', all_traces)

    Creates: 
    
        DAYS
        ├── 1R.2020.346.png


    Make sure that all_traces[0] contains full trace-id metadata. 

    Created for ROCKETSEIS project.

    """    

    daysDir = os.path.join(TOPDIR, 'DAYPLOTS')
    os.makedirs(daysDir)
    NETWORK = all_traces[0].stats.network
    stime = all_traces[0].stats.starttime
    YYYY = stime.year
    JJJ = stime.yearday
    pngfile = os.path.join(daysDir, '%s.%s.%s.png' % (NETWORK, YYYYJJJ[0:4], YYYYJJJ[4:])  )   
    all_traces.plot(equal_scale=False, outfile=pngfile);
    return


def Stream_to_24H(all_traces):
    """
    Take a Stream object, merge all traces with common ids and pad out to 24-hour-long traces

    Created for ROCKETSEIS data conversion and modified for CALIPSO data archive from Alan Linde. 
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
    Clever way to merge overlapping traces into a BUD file. Uses all non-zero data values from both.

    Created for CALIPSO data archive from Alan Linde, when needed to upgrade Stream_to_BUD.
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
    """ calculate predicted travel times based on IASP91 model  - see https://docs.obspy.org/packages/obspy.taup.html
        Input: station and quake both are dicts with lat and lon keys
        Output: a phases dict is added to statihttps://www.facebook.com/on, with phase name keys and predicted arrival times """
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
    """ Generate synthetics for a GCMT event, save into an mseedfile, return as Stream object """
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

def read_DMX_file(DMXfile, fix=True, defaultnet=''):
    # DMX read support now (2023) included in ObsPy. Was not available for the Montserrat ASN conversion in 2019.
    # This produces same result as converting DMX to SAC with sud2sac.exe in Win-SUDS, and then reading into ObsPy
    # Has also been tested against sud2gse.exe.
    # sud2msed.exe is messier, because that program completely loses all tr.id info when converting, so all tr.id set to ...
    # Tested on data from Montserrat 1995-6 and Pinatubo 1991
    #
    # ObsPy DMX reader inserts "unk" in place of an unknown network. We do not want this.
    #
    # ObsPy DMX reader reads DMXfile as uint16 and so is all +ve. 
    # sud2sac.exe converts to numbers either side of 0. 
    # Subtracting 2048 from each sample of tr.data corrects data read in using ObsPy DMX reader to match that from SAC
    #
    # Obspy Miniseed writer needs float, not int, so recast as float.
    #
    # Passing fix=False will just run ObsPy DMX reader without applying any corrections.

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



def parse_hypo71_line(line):
    """
    Parses a single line of HYPO71 output format using fixed column positions.
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
    Parses a HYPO71 earthquake catalog file into an ObsPy Catalog object.
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