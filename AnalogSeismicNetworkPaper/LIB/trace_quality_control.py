import numpy as np
def check0andMinus1(liste):
# function bool_good_trace = check0andMinus1(tr.data)
    liste=list(liste)
    listStr=''.join(str(i) for i in liste)
    if  "000000000000" in listStr or "-1-1-1-1-1-1-1-1" in listStr :
        return False
    else:
        return True

def signaltonoise(tr):
# function snr, highval, lowval = signaltonoise(tr)
    # Here we just make an estimate of the signal-to-noise ratio
    #
    # Normally the trace should be pre-processed before passing to this routine, e.g.
    # * remove ridiculously large values
    # * remove any sequences of 0 from start or end
    # * detrend
    # * bandpass filter
    #
    # Processing:
    #    1. ensure we still have at least 10 seconds
    #    2. take absolute values
    #    3. compute the maximum of each 1-s of data, call this time series M
    #    4. compute 95th and 5th percentile of M, call these M95 and M5
    #    5. estimate signal-to-noise ratio as M95/M5
    
    highval = -1
    lowval = -1
    snr = -1 
    a = tr.data

    fsamp = int(tr.stats.sampling_rate)
    npts = tr.stats.npts
    numseconds = int(npts/fsamp)
    if numseconds > 10:
        a = a[0:int(fsamp * numseconds - 1)]             # remove any fractional second from end of trace
        abs = np.absolute(a)                             # take absolute value of a        
        abs_resize = np.resize(abs, (fsamp, numseconds)) # resize so that each second is one row
        M = np.max(abs_resize,axis=0)                    # find max of each second / row
        highval = np.nanpercentile(M,95)                    # set highval to 95th percentile, to represent signal amplitude
        if highval < 1:
            highval = -1
            return (snr, highval, lowval)
        lowval = np.nanpercentile(M,5)                      # set lowval to 5th percentile, to represent noise amplitude
        snr = highval / lowval
        print(abs_resize.shape)
        print(M.shape)
    return (snr, highval, lowval,)

def compute_metrics(tr):
# function tr, quality_factor, snr = compute_metrics(tr)
# This function wraps all others
    # Here we compute simple metrics on each trace and write them to NET.STA.CHAN files. 
    # These metrics are:
    #     1. duration of signal
    #     2. signal amplitude
    #     3. noise amplitude
    #     4. signal-to-noise ratio

    duration = tr.stats.npts /  tr.stats.sampling_rate
    quality_factor = trace_quality_factor(tr)
    snr = (-1, -1, -1)
    if quality_factor > 0:

        tr = clean_trace(tr)
        tr = fix_nslc_montserrat(tr)

        # estimate signal-to-noise ratio
        snr = signaltonoise(tr)
        if snr[0] <= 1:
            quality_factor = quality_factor * 0.5

    return (tr, quality_factor, snr)

 
def clip_trace(tr):     
# function tr = clip_trace(tr)
    # remove absurdly large values
    AMP_LIMIT = 10000000
    a = tr.data
    np.clip(a, -AMP_LIMIT, AMP_LIMIT, out=a)
    np.where(a == AMP_LIMIT, 0, a)
    tr.data = a
    return tr

def change_last_sample(tr):
# function tr = change_last_sample(tr)
    # For some SEISAN files from Montserrat - possibly from SEISLOG conversion,
    # the last value in the time series was always some absurdly large value
    # This sections was to change last value to be the mean of the rest
    #a = tr.data
    #a = np.delete(a,[np.size(a) - 1])
    #m = np.mean(a)
    #tr.data[-1] = m
    return tr

def trace_quality_factor(tr):
# function quality_factor = trace_quality_factor(tr)
    quality_factor = 1

    # ignore blank trace
    anyData = np.count_nonzero(tr.data)
    if anyData==0:
        quality_factor = 0 
        return quality_factor

    # check for sequences of 0 or 1
    trace_good_flag = check0andMinus1(tr.data)

    if not trace_good_flag:
        quality_factor = 0
        return quality_factor

    return quality_factor

def clean_trace(tr):
# function tr = clean_trace(tr)

    # remove absurd values
    tr = clip_trace(tr)

    # Now high-pass filter, plot and compute_metrics on each trace
    try:
        tr.filter('highpass', freq=0.5, corners=2, zerophase=True)
    except:
        print('Filtering failed')
        print(tr.stats)

    # remove linear trend
    try:
        tr.detrend(type='linear')
    except:
        print('Detrending failed')
        print(tr.stats)

    return tr

def fix_nslc_montserrat(tr):
# function tr = fix_nslc_montserrat(tr)
    # FIX NET, STA, LOC, CHAN CODES ###
    # fix the network, channel and location
    network = 'MV'
    tr.stats['network']=network
    sta = tr.stats['station'].strip()
    chan = tr.stats['channel'].strip()
    if chan=='PRS' or chan[0:2]=='AP':
        chan='BDF'
    else:
        if chan[0]=='A':
           if tr.stats['location'] == 'J':
               bandcode = 'S'
           else:
               bandcode = chan[0]
        else:
            try:
                if chan[1]=='B':
                    bandcode = 'B'
                else:
                    bandcode = chan[0]
            except:
                bandcode = chan[0]
        instrumentcode = 'H'
        if len(chan)<3:
            orientationcode = tr.stats['location']
        else:
            orientationcode = chan[2]
                
        chan = bandcode + instrumentcode + orientationcode

    if chan[0]=='A':
        print(tr.stats)
        print(chan)
        sys.exit('bugger!')
    tr.stats['channel'] = chan
    tr.stats['location']='--'
    return tr

def swap32(i):
# function y = swap(tr.data)
    return struct.unpack("<i", struct.pack(">i", i))[0]
