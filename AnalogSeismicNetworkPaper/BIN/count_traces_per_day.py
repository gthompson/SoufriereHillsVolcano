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




# Load each monthly catalog
#SEISAN_DATA = os.environ['SEISAN_DATA']
#CSVpath = os.path.dirname('/media/sdd1/seismo')
#os.chdir(CSVpath)
#list_of_csv_files = sorted(glob.glob('/media/sdd1/seismo/MVOE_catalog*.csv'))
#list_of_csv_files = sorted(glob.glob('/media/sdd1/seismo/ASNE_wavfiles*.csv'))
list_of_csv_files = sorted(glob.glob('CSVDB/ASNE_wavfiles*.csv'))
print(list_of_csv_files)
all_trace_ids = list()
# get a list of all traceids
for csvfile in list_of_csv_files:
    print(csvfile)
    df = pd.read_csv(csvfile)
    #print(df)
    df.columns = df.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')', '')
    # loop over every day in this month & make a summary row for each channel
    new_trace_ids = df.traceid.unique()
    all_trace_ids = list(set().union(all_trace_ids, new_trace_ids))   
print(all_trace_ids)
dailytraceid_df = pd.DataFrame(columns=['yyyymmdd', 'traceid', 'count']) 
for csvfile in list_of_csv_files:
    print(csvfile)
    df = pd.read_csv(csvfile)
    df.columns = df.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')', '')
    yyyy = int(csvfile[-10:-6])
    mm = int(csvfile[-6:-4])
    start_date = dt.date(yyyy, mm, 1)
    emm = mm+1
    eyyy=yyyy
    if emm==13:
       emm=1
       eyyy = yyyy+1
    end_date = dt.date(eyyy, emm, 1)
    delta = dt.timedelta(days=1)
    while start_date < end_date: # loop over all days in this month
        yyyymmdd = start_date.strftime("%Y%m%d")
        print(yyyymmdd)
        sfiledatetime = pd.to_datetime(df.datetime)
        daycat = df.loc[sfiledatetime.dt.date == start_date] # get rows just for this day
        #daytraceid_hash = daycat.traceid.value_counts() 
        print(daycat)
        print(daycat.shape)
        for thistraceid in all_trace_ids:
            if not daycat.empty: 
                print(thistraceid)
                try:
                    thisdaytrace_df = daycat.loc[daycat['traceid'] == thistraceid]
                    thiscount = thisdaytrace_df.shape[0] 
                    #print('*** it works ***')
                except:
                    print('thistraceid = %s' % thistraceid)
                    print(daycat)
                    barf 
            else:
                thiscount = 0
            #thisrow = pd.DataFrame(yyyymmdd, thistraceid, thiscount)
            thisrow = {'yyyymmdd':yyyymmdd, 'traceid':thistraceid, 'count':thiscount}
            dailytraceid_df = dailytraceid_df.append(thisrow, ignore_index=True)
            #print(thisrow)
                
        start_date += delta # add a day  
    #print('Done with this csvfile')
dailytraceid_df.to_csv('../DATA/ASNE_dailytraceid_wavfiles_df.csv')  
