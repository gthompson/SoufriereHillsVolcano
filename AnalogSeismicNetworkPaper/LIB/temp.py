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
    all_trace_ids = list(set().union(all_trace_ids, new_trace_ids))   
print(all_trace_ids)

for oldid in all_trace_ids:
    newid = qc.fix_trace_id(oldid)
    print(oldid, '->', newid)

