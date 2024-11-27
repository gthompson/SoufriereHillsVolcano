# -*- coding: utf-8 -*-
"""
Created on Mon Mar 02 13:02:50 2015

@author: thompsong
"""

from obspy import read
import os.path
import numpy as np
import datetime
import sys
import matplotlib.pyplot as plt
import struct
import os
import glob
def qnxsei2obspy(file, antelopetop, db ):
    st = read(file)
    for tr in st:
        data = tr.data
        for i in range(0, len(data) -1):
            data[i] = swap32(data[i])
        tr.data = data
        tr.detrend(type='constant')
        #tr.plot()
        # remove blank traces check for all 0
        a = tr.data
        a = np.delete(a,[np.size(a) - 1])
        m = np.mean(a)
        # change last value to be the mean of the rest
        tr.data[-1] = m
        anyData = np.count_nonzero(tr.data)
        if anyData==0:
            #st.remove(tr)
            continue        
        #raw_input("Press Enter to continue...")
        #st.plot()

        # fix the network, channel and location
        #print " "
        #print tr.stats
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
            print tr.stats
            print chan
            sys.exit('bugger!')
        tr.stats['channel'] = chan
        tr.stats['location']='--'
        #continue

        # write out to SDS structure - but this only works for full day files
        #filepath = "/raid/data/SDS/" + year + "/" + network + "/" + tr.stats['station'] + "/" + tr.stats['channel'] + "D"
        #filename = tr.id + ".D." + year + "." + julday
        #tr.write(filename + ".MSEED", format="MSEED")

        # get sta, chan
        sta = tr.stats['station']
        chan = tr.stats['channel']

        # get time strings for this trace
        starttime = tr.stats['starttime']
        jjj = str(starttime.julday).zfill(3)
        yyyy = str(starttime.year).zfill(4)
        mm = str(starttime.month).zfill(2)
        dd = str(starttime.day).zfill(2)
        hh = str(starttime.hour).zfill(2)
        ii = str(starttime.minute).zfill(2)
        ss = str(starttime.second).zfill(2)

        # Write Miniseed files to Antelope structure YYYY/JJJ/sta.chan.yyyy:jjj:hh:mm:ss
        antelopedir = os.path.join(antelopetop, db, yyyy, jjj)
        if not os.path.exists(antelopedir):
            os.makedirs(antelopedir)
        #antelopefile = sta + "." + chan + "." + yyyy + ":" + jjj + ":" + hh + ":" + ii + ":" + ss
        antelopefile = sta + "." + chan + "." + yyyy + "_" + jjj + "_" + hh + "_" + ii + "_" + ss
        antelopefullpath = os.path.join(antelopedir, antelopefile)
        print antelopefullpath
        tr.write(antelopefullpath, format="MSEED")
        #input("Press Enter to continue...")
    #dirbase, ext = os.path.splitext(file)
    #outfile = dirbase + '.mseed'
    #st.write(outfile, 'MSEED')
    # write SAC or MSEED files to a new Seisan database
    #(filename,ext) = os.path.splitext(seisanbasename)
    #st.write(filename + ".MSEED", format="MSEED")



def swap32(i):
    return struct.unpack("<i", struct.pack(">i", i))[0]

def main():
    # read in the wav file name from the command line
    antelopetop = '/raid/data/antelope'
    #antelopetop = 'D:\\roxio\\seisan_from_qnx'
    db = 'dbDSNC_SEISLOG'
    os.chdir('D:\\roxio\\seislog\\2000\\03')
    print os.getcwd()
    years = glob.glob('[0-9]' * 4)
    for year in years:
        os.chdir(year)
        print os.getcwd()
        months = glob.glob('[0-9]' * 2)
        for month in months:
            os.chdir(month)
            print os.getcwd()
            #seislogfiles = glob.glob('*T.MVO*')
            seislogfiles = glob.glob('*S.MVO*')
            print("Found %d Seislog files" % (len(seislogfiles)))
            for seislogfile in seislogfiles:
                qnxsei2obspy(seislogfile, antelopetop, db)
            os.chdir('..')
        os.chdir('..')

if __name__ == "__main__":
	main()