#!/usr/bin/env python
#import obspy.core
def fix_trace_ids(tr):
    # fix the network, channel and location
    network = 'MV'
    tr.stats['network']=network
    sta = tr.stats['station'].strip()
    chan = tr.stats['channel'].strip()
    if chan=='PRS' or chan=='APS':
        chan='BDF'
    else:
        if chan[0]=='A':
            if tr.stats['location'] == 'J':
                bandcode = 'S'
            #else:
                #bandcode = 'B'
        else:
            if chan[1]=='B':
                bandcode = 'B'
            else:
                bandcode = chan[0]
            instrumentcode = 'H'
            if len(chan)==2:
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
    
