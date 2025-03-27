#!/usr/bin/env python
import sys
#import numpy as np
import os
import pandas as pd
from glob import glob
from obspy import read_inventory, read, Stream, Trace

LIBpath = os.path.join( os.getenv('HOME'),'src','kitchensinkGT', 'LIB')
sys.path.append(LIBpath)
from lib.libseisGT_old3 import get_seed_band_code, fix_trace_id, remove_empty_traces
from metrics import process_trace, ampengfft
sys.path.append(os.path.join( os.getenv('HOME'),'src', 'icewebPy') )
import IceWeb


def change_last_sample(tr):
    # For some SEISAN files from Montserrat - possibly from SEISLOG conversion,
    # the last value in the time series was always some absurdly large value
    # So remove the last sample
    tr.data = tr.data[0:-2]

def swap32(i):
    # Change the endianess
    return struct.unpack("<i", struct.pack(">i", i))[0]

'''
def inventory_fix_id_mvo(inv):
    inv[0].code='MV'
    net = inv[0].code
    for station in inv[0].stations:
        sta = station.code
        for channel in station.channels:
            chan = channel.code
            if chan[0] in 'ES':
                shortperiod = True
            if chan[0] in 'BH':
                shortperiod = False
            Fs = channel.sample_rate
            nslc = net + '.' + sta + '..' + chan
            nslc = correct_nslc(nslc, Fs, shortperiod=shortperiod)
            net, sta, loc, chan = nslc.split('.')
            channel.code = chan
        station.code = sta
    inv[0].code = net
    return inv

'''


def load_mvo_inventory(tr, CALDIR):
    this_inv = None
    matchcode = None
    if len(tr.stats.channel)<3:
        return this_inv
    if tr.stats.channel[0] in 'ES':
        matchcode = '[ES]'
    elif tr.stats.channel[0] in 'BH':
        matchcode = '[BH]'
    if not matchcode:
        print("Cannot match trace ID %s ",tr.id)
        return this_inv
    #xmlfilepattern = os.path.join(CALDIR, "station.MV.%s..%s*%s.xml" % (tr.stats.station, matchcode, tr.stats.channel[2]) )
    xmlfiles = glob(os.path.join(CALDIR, "station.MV.%s..%s*%s.xml" % (tr.stats.station, matchcode, tr.stats.channel[2]) ))
    N = len(xmlfiles)
    if N==1:
        xmlfile = xmlfiles[0]
        print('Correcting %s with %s' % (tr.id, xmlfile))
        this_inv = read_inventory(xmlfile)  
    return this_inv

 
# new functions from enhanced stream object workflow 2021/11/23
def metrics2df(st):
    tracedf = pd.DataFrame()
    list_of_tracerows = []
    for tr in st:
        s = tr.stats
        tracerow = {'id':tr.id, 'starttime':s.starttime, 
               'Fs':s.sampling_rate, 
               'calib':s.calib, 'units':s.units, 
               'quality':s.quality_factor}
        if 'spectrum' in s: 
            for item in ['medianF', 'peakF', 'peakA', 'bw_min', 'bw_max']:
                try:
                    tracerow[item] = s.spectrum[item]
                except:
                    pass
        if 'metrics' in s:
            m = s.metrics
            for item in ['snr', 'signal_level', 'noise_level', 'twin',
                         'peakamp', 'peaktime', 'energy', 'RSAM_high', 'RSAM_low',
                         'sample_min', 'sample_max', 'sample_mean', 'sample_median', 
                         'sample_lower_quartile', 'sample_upper_quartile', 'sample_rms', 
                         'sample_stdev', 'percent_availability', 'num_gaps', 'skewness', 'kurtosis']:
                         #'start_gap', 'num_gaps', 'end_gap', 'sum_gaps', 'max_gap', 
                         #'num_overlaps', 'sum_overlaps', 'num_records', 'record_length', 
                try:
                    tracerow[item] = m[item]
                except:
                    pass 
        if 'bandratio' in s:
            for dictitem in s['bandratio']:
                label = 'bandratio_' +  "".join(str(dictitem['freqlims'])).replace(', ','_')
                tracerow[label] = dictitem['RSAM_ratio']
        if 'lon' in s:
            tracerow['lon'] = s['lon']
            tracerow['lat'] = s['lat']
            tracerow['elev'] = s['elev']

        list_of_tracerows.append(tracerow)
    tracedf = pd.DataFrame(list_of_tracerows)
    tracedf = tracedf.round({'Fs': 2, 'secs': 2, 'quality':2, 'medianF':1, 'peakF':1, 'bw_max':1, 'bw_min':1, 'peaktime':2, 'twin':2, 'skewness':2, 'kurtosis':2})
    #tracedf.set_index('id')
    return tracedf



def enhance_stream(stream_in, CALDIR=None, master_inv=False, quality_threshold=0.0):
    # if CALDIR provided, will try to correct traces
    st = stream_in.copy()
    for tr in st:
        this_inv = None
        if CALDIR: # this means user wants to correct traces
            if master_inv:
                this_inv = master_inv
            else:
                this_inv = load_mvo_inventory(tr, CALDIR)
            if tr.stats.channel[0:2]=='SN': # accelerometer channels
                tr.stats.calib = 27500000
                tr.data = tr.data / tr.stats.calib # approx calib from comparing with co-located differentiated LA100 waveforms
                tr.stats.units = 'm/s2'            
        process_trace(tr, inv=this_inv, quality_threshold=quality_threshold)
            
    # remove bad traces
    for tr in st:    
        if tr.stats.quality_factor <= quality_threshold:
            st.remove(tr)
            
    if CALDIR: # remove traces not corrected
        physical_units = ['m', 'm/s', 'm/s2', 'Pa']
        for tr in st:    
            if not tr.stats.units in physical_units:
                st.remove(tr)        
        
    # add AEF metrics        
    iwsobj = IceWeb.icewebSpectrogram(stream=st)
    iwsobj = iwsobj.precompute() # spectrograms data added
    iwsobj.compute_amplitude_spectrum(compute_bandwidth=True) # adds tr.stats.spectrum
    for tr in iwsobj.stream:
        ampengfft(tr) # add peaktime, peakamp, energy, frequency metrics 
        #tr.stats.pop('spectrogramdata', None) # Remove spectrogramdata as it is large
    enhanced_stream = iwsobj.stream 
    
    return enhanced_stream

#

def save_enhanced_stream(st, eventdf, enhanced_wavpath, save_pickle=False):

    if enhanced_wavpath[-6:] == '.mseed':
        enhanced_wavpath = enhanced_wavpath[:-6]

    parentdir = os.path.dirname(enhanced_wavpath)
    if not os.path.exists(parentdir):
        os.makedirs(parentdir)

    # save Miniseed
    st.write(enhanced_wavpath + '.mseed', 'MSEED')
    
    # metrics dataframe
    eventdf = metrics2df(st)
    eventdf.to_csv(enhanced_wavpath + '.csv',index=False)
        
    # write to pickle file?
    if save_pickle:
        st.write(enhanced_wavpath + '.pickle', 'PICKLE') # Rewrite pickle file with extra attributes     

#        
        
def read_enhanced_stream(enhanced_wavpath):
    enhanced_wavpath = enhanced_wavpath.replace('.mseed','')    

    # read Miniseed
    st = read(enhanced_wavpath + '.mseed', 'MSEED')
    
    # read metrics dataframe
    eventdf = pd.read_csv(enhanced_wavpath + '.csv', index_col=False)

    for tr in st:
        s = tr.stats
        row = eventdf[eventdf.id==tr.id].iloc[0].to_dict()
        s.units = row['units']
        s.calib = row['calib']
        s['quality_factor'] = row['quality']
        s['spectrum'] = {}
        for item in ['medianF', 'peakF', 'peakA', 'bw_min', 'bw_max']:
            try:
                s.spectrum[item] = row[item]
            except:
                pass
        s['metrics']={}
        for item in ['snr', 'signal_level', 'noise_level', 'twin',
                     'peakamp', 'peaktime', 'energy', 'RSAM_high', 'RSAM_low',
                     'sample_min', 'sample_max', 'sample_mean', 'sample_median', 
                     'sample_lower_quartile', 'sample_upper_quartile', 'sample_rms', 
                     'sample_stdev', 'percent_availability', 'num_gaps', 'skewness', 'kurtosis']:
            try:
                s.metrics[item] = row[item]
            except:
                pass 
        s['bandratio']=[{'freqlims': [1.0, 6.0, 11.0], 'RSAM_ratio': row['bandratio_[1.0_6.0_11.0]']}, 
                        {'freqlims': [0.8, 4.0, 16.0], 'RSAM_ratio': row['bandratio_[0.8_4.0_16.0]']}]
        for item in ['lon', 'lat', 'elev']:
            try:
                s[item] = row[item]
            except:
                pass
    
    return st        

#

def respfiles2masterstationxml(SEISAN_DATA, xmlfile):# Merge Montserrat RESP files
    # A function we are only likely to use once, but an important one to keep
    from obspy.io.xseed.core import _read_resp
    from lib.libseisGT_old3 import create_dummy_inventory, merge_inventories
    station0hypfile = os.path.join(SEISAN_DATA, 'DAT', 'STATION0_MVO.HYP')
    station_locationsDF = parse_STATION0HYP(station0hypfile) 
    respfiles = glob.glob(os.path.join(SEISAN_DATA, 'CAL', 'RESP*'))
    master_inv = create_dummy_inventory()
    for respfile in respfiles:
        this_inv = None
        print('RESP file = ',respfile)
        this_inv = _read_resp(respfile)
        this_inv = inventory_fix_id_mvo(this_inv)
        sta_code = this_inv.networks[0].stations[0].code
        location = station_locationsDF[station_locationsDF['name']==sta_code].iloc[0].to_dict()
        for station in this_inv.networks[0].stations:
            station.latitude = location['lat']
            station.longitude = location['lon']
            station.elevation = location['elev']
            for channel in station.channels:
                channel.latitude = location['lat']
                channel.longitude = location['lon']
                channel.elevation = location['elev']
                channel.depth = 0.0
                if channel.sample_rate==75.19:
                    channel.sample_rate=75.0
        merge_inventories(master_inv, this_inv)  
    master_inv.write(xmlfile,format="STATIONXML", validate=True)


if __name__ == '__main__':
    pass
