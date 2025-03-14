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
def fix_trace_id(st, shortperiod=False):
    # convenience method to wrap correct_nslc

    for tr in st:
        if not 'history' in tr.stats:
            tr.stats['history'] = list() 
        nslc = correct_nslc(tr.id, tr.stats.sampling_rate, shortperiod=shortperiod)
        tr.id = nslc
        if 'deconvolved' in tr.stats.history or 'calibrated' in tr.stats.history:
            if tr.stats.channel[1] == 'D':
                tr.stats['units'] = 'Pa'
            if tr.stats.channel[1] == 'H':
                tr.stats['units'] = 'm/s'    
        else:
            tr.stats['units'] = 'Counts'
'''

def fix_trace_id_mvo(trace, legacy=False, netcode='MV'):
    fix_times(trace)
    fix_sample_rate(trace)
    trace.id = correct_nslc(trace.id, trace.stats.sampling_rate)
    fix_trace_id(trace, legacy=legacy, netcode=netcode)

def correct_nslc(traceID, Fs, shortperiod=None):
    # Montserrat trace IDs are often bad. return correct trace ID
    # also see fix_nslc_montserrat in /home/thompsong/Developer/SoufriereHillsVolcano/AnalogSeismicNetworkPaper/LIB/fix_mvoe_traceid.ipynb
    # special case - based on waveform analysis, this trace is either noise or a copy of MV.MBLG..SHZ
    if traceID == '.MBLG.M.DUM':
        traceID= 'MV.MBLG.10.SHZ'
    traceID = traceID.replace("?", "x")

    oldnet, oldsta, oldloc, oldcha = traceID.split('.')

    net = 'MV'    
    sta = oldsta.strip()
    loc = oldloc.strip()
    chan = oldcha.strip()

    if loc == '--' or loc == 'J' or loc=='I':
        loc = ''    

    if not shortperiod:
        if 'SB' in chan or chan[0] in 'BH':
            shortperiod=False
        else:
            shortperiod = True
        bandcode = get_seed_band_code(Fs, shortperiod=shortperiod)

    bandcode = get_seed_band_code(Fs, shortperiod=shortperiod) # not sure here if BB or SP sensor
    instrumentcode = 'H'
    orientationcode = 'x'
    if 'Z' in loc or 'Z' in chan:
        orientationcode = 'Z'
    elif 'N' in loc or 'N' in chan:
        orientationcode = 'N'
    elif 'E' in loc or 'E' in chan:
        orientationcode = 'E'    


    if (sta=='MBLY' and chan[0]=='P') or 'AP' in chan or 'PR' in chan or 'PH' in chan or chan=='S A':
        instrumentcode = 'D'
        orientationcode = 'F'
        if chan[-1].isnumeric():
            loc = chan[-1].zfill(2)
        elif loc.isnumeric():
            loc = loc.zfill(2)
        else:
            loc = ''


    elif len(chan)>1:
        if chan[1].strip():
            instrumentcode = chan[1] 
        if len(chan)>2:
            orientationcode = chan[2]
            if orientationcode=='H':
                orientationcode='x'
        
    #print(f'chan here {chan}: {bandcode},{instrumentcode},{orientationcode}')  
    
    # channel code is bandcode + instrumentcode + orientationcode
    '''
    if len(chan)>0:
        if chan[0]=='E' or chan[0]=='S':
            shortperiod=True
        if chan[0]=='B' or chan[0]=='H':  
            shortperiod=False    
    

 
    if not instrumentcode  in 'HLD':
        instrumentcode = 'H' # 'H' seismic velocity sensor is default
    '''
    
    # Montserrat BB network 1996-2004 had weirdness like
    # BB stations having channels 'SB[Z,N,E]' and
    # SP stations having channels 'S [Z,N,E]'
    # location code was usually 'J' for seismic, 'E' for pressure
    # channel was 'PRS' for pressure
    # there were also 'A N' channels co-located with single-component Integra LA100s, so perhaps those were some other
    # type of seismometer, oriented North?
    # let's handle these directly here
    if len(chan)==2:

        # could be a 2006 era waveform trace ID where given as .STAT.[ZNE].[BS]H
        if len(loc)==1:
            #chan=chan+loc # now length 3
            if loc in 'ZNE':
                orientationcode = loc
                loc = ''
            #if not loc.isnumeric():
            #    loc='' 
        elif len(loc)==0:
            # could be arrival row from an Sfile, where the "H" is omitted
            # or an AEF line where trace dead and orientation missing
            #instrumentcode = 'H'
            if chan[1] in 'ZNE':
                orientationcode = chan[1]
            #else:
            #    orientationcode = '' # sometimes get two-character chans from AEF lines which omit component when trace is dead, e.g. 01-0954-24L.S200601, 


    
    elif len(chan)==3:
        
        if chan[0:2]=='SB':
             # just because we know it is BB sensor
            instrumentcode = 'H'
        """    
        elif chan[0:2]=='S ':
            bandcode = get_seed_band_code(Fs, shortperiod=True) # just because we know it is SP sensor
            instrumentcode = 'H'
        
        if chan[0:3]=='PRS':
            instrumentcode = 'D'
            orientationcode = 'F'
        elif chan[0:3]=='A N':
            instrumentcode = 'N'
            orientationcode = 'Z'
            #loc = '10'
        elif chan[0]=='P' or chan[0:2]=='AP' or chan=='PRS': # just because we know it is SP sensor
            instrumentcode ='D' # infrasound/acoustic
            orientationcode = 'F'
            if chan[1].isnumeric(): # e.g. channel like P5
                loc = chan[1].zfill(2)
            if chan[2].isnumeric(): # e.g. channel like P5
                loc = chan[2].zfill(2)                   
        
        elif chan[1] in 'ZNE': # seismic component in wrong position
            bandcode = get_seed_band_code(Fs, shortperiod=shortperiod)
            instrumentcode = 'H'  
            orientationcode = chan[1]
        
        elif len(chan)==2 and instrumentcode == 'H': # e.g. .MBLY.5.SH
            orientationcode = loc
            loc = ''  
        """


    '''
    if len(loc)>0:
        if loc[0]=='J' or loc[0]=='E' or loc == "--":
            if len(loc)>1:
                loc=loc[1:]
            else:
                loc=''   
    '''

    chan = bandcode + instrumentcode + orientationcode

    newID = net + "." + sta + "." + loc + "." + chan
    #print(traceID,'->',newID)
    return newID
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

def inventory_fix_ids(inv, netcode='MV'):
    for network in inv.networks:
        if not network.code:
            network.code = netcode
        for station in network.stations:
            for channel in station.channels:
                if channel.code[0] in 'FCES':
                    shortperiod = True
                if channel.code[0] in 'GDBH':
                    shortperiod = False
                Fs = channel.sample_rate
                nslc = network.code + '.' + station.code + '..' + channel.code
                nslc = correct_nslc(nslc, Fs, shortperiod=shortperiod)
                net, sta, loc, chan = nslc.split('.')
                channel.code = chan
            station.code = sta
#    return inv


def load_mvo_master_inventory(XMLDIR):
    master_station_xml = os.path.join(XMLDIR, 'MontserratDigitalSeismicNetwork.xml')
    if os.path.exists(master_station_xml):
        print('Loading ',master_station_xml)
        return read_inventory(master_station_xml)
    else:
        print('Could not find ',master_station_xml)        
        return None

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


def fix_sample_rate(st, Fs=75.0):
    if isinstance(st, Stream):
        for tr in st:
            fix_sample_rate(tr, Fs=Fs)  # Recursive call for each Trace
    elif isinstance(st, Trace):
        tr = st
        if tr.stats.sampling_rate > Fs * 0.99 and tr.stats.sampling_rate < Fs * 1.01:
            tr.stats.sampling_rate = Fs 
    else:
        raise TypeError("Input must be an ObsPy Stream or Trace object.")    


#from obspy.core.utcdatetime import UTCDateTime
def fix_times(st):
    if isinstance(st, Stream):
        for tr in st:
            fix_sample_rate(tr, Fs=75.0)  # Recursive call for each Trace
    elif isinstance(st, Trace):
        tr = st
        yyyy = tr.stats.starttime.year
        if yyyy == 1991 or yyyy == 1992 or yyyy == 1993: # digitize
            # OS9/Seislog for a while subtracted 8 years to avoid Y2K problem with OS9
            # should be already fixed but you never know
            tr.stats.starttime._set_year(yyyy+8)
        if yyyy < 1908:
            # getting some files exactly 100 years off
            tr.stats.starttime._set_year(yyyy+100) 
    else:
        raise TypeError("Input must be an ObsPy Stream or Trace object.")  

                      
            
# bool_ASN: set True only if data are from MVO analog seismic network
def read_monty_wavfile_and_correct_traceIDs(wavpath, bool_ASN=False, verbose=False, \
                                            seismic_only=False, vertical_only=False):
    if os.path.exists(wavpath):
        st = read(wavpath)
    else:
        print('ERROR. %s not found.' % wavpath)
        return Stream()
    
    if vertical_only:
        st = st.select(component='Z')
    elif seismic_only:
        for tr in st:
            if not tr.stats.channel[1] in 'HL':
                st.remove(tr)
    remove_empty_traces(st)
    for tr in st:
        fix_trace_id_mvo(tr)
    '''
    if verbose:
        print('- checking sampling rate')
    fix_sample_rate(st) # set any sample rates around 74.3 - 75.7 Hz to 75.0 Hz
    if verbose:
        print('- checking trace ID')    
    fix_trace_id(st, shortperiod=bool_ASN) 
    if verbose:
        print('- checking year not before 1995')    
    fix_times(st)
    '''

    return st

#

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
