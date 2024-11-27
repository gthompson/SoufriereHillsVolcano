#!/usr/bin/env python
import pandas as pd
import csv
import obspy
import matplotlib.pyplot as plt
import numpy as np
import trace_quality_control as qc 
import os
import glob


# Load the catalogue
SEISAN_DATA = os.environ['SEISAN_DATA']
CSVpath = os.path.dirname(SEISAN_DATA)
os.chdir(CSVpath)
list_of_csv_files = sorted(glob.glob('ASNE_catalog??????.csv'))
#frames = []
for csvfile in list_of_csv_files:
    df = pd.read_csv(csvfile)
    #frames.append(df)

    # fix column names
    df.columns = df.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')', '')

    nrows, ncolumns = df.shape
    print(nrows)
    quality_index = np.zeros(nrows)
    df['quality_index'] = quality_index
    fixed_ids = df.traceid.tolist() # a list of all trace ids in the csvfile

    for i in range(nrows):
        if i % 100 == 0:
            print('Done %d of %d' % ( i, nrows, ))

        # code from fix_montserrat_traceids.ipynb
        traceid = fixed_ids[i]
        wavfilepath = mvocat.iloc[i].wavfilepath
        #print(wavfilepath)
        if os.path.exists(wavfilepath):
            st = obspy.read(wavfilepath)
            tr = st[df.iloc[i].tracenum]
            print(tr)
            tr = qc.compute_metrics(tr)
            tr2 = fix.fix_trace_ids(tr)
            fixed_ids[i]=tr2.id # fix the list
        else:
            print('does not exist')

        # code from qc_traces_from_CSV.ipynb
        #traceid = mvocat.iloc[i].traceid
        #wavfile = mvocat.iloc[i].wavfilepath
        #if wavfile.empty:
        #    continue
        #try:
        #    os.path.exists(wavfile)
        #except:
        #    print(wavfile[40:])
        #try:
        #    if not os.path.exists(wavfile):
        #        continue
        #except:
        #    continue
        #st = obspy.read(wavfile).select(id = traceid)
        if not len(st)==1:
            continue
        data = st[0].data
        if check0andMinus1(data):
            quality_index[i] = 1
        else:
            pass
            #print("There is consecutive 0 or -1. in this recording")
    print('Done!')
    df['quality_index'] = quality_index
    qc_csvfile = csvfile[0:12] + "_qc_.csv"
    df.to_csv(qc_csvfile)
