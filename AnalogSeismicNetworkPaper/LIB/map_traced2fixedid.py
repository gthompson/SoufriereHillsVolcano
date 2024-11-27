#!/usr/bin/env python
import pandas as pd
import csv
import obspy
import matplotlib.pyplot as plt
import numpy as np
import trace_quality_control as qc
import os
import glob
import datetime as dt

def fix_nslc_montserrat(traceid):
# function tr = fix_nslc_montserrat(tr)
    # FIX NET, STA, LOC, CHAN CODES ###
    # fix the network, channel and location
    oldnet, oldsta, oldloc, oldcha = traceid.split('.')
    print(oldnet, oldsta, oldloc, oldcha)
    return traceid
    network = 'MV'
    chan = oldcha
    sta = oldsta
    loc = oldloc
    if chan=='PRS' or chan[0:2]=='AP':
        chan='BDF'
    else:
        if chan[0]=='A':
           if loc == 'J':
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
            orientationcode = oldloc
        else:
            orientationcode = chan[2]

        chan = bandcode + instrumentcode + orientationcode

    if chan[0]=='A':
        print(network, sta, chan, loc)
        sys.exit('bugger!')
    return network + "." + sta + "." + loc + ".' + chan

days_df = pd.DataFrame()
list_of_csv_files = sorted(glob.glob('./MVOE_wavfiles*.csv'))
print(list_of_csv_files)
all_trace_ids = list()
# get a list of all traceids
for csvfile in list_of_csv_files:
    print(csvfile)
    df = pd.read_csv(csvfile)
    #print(df)
    df.columns = df.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')', '')
    # loop over every day in this month & make a summary row for each channel
    old_trace_ids = df.traceid.unique()
    all_trace_ids = list(set().union(all_trace_ids, old_trace_ids))   
print(all_trace_ids)

for oldid in all_trace_ids:
    newid = fix_nslc_montserrat(oldid)
    print(oldid, '->', newid)

