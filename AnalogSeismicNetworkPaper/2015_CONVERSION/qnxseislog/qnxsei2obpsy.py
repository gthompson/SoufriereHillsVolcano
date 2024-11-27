from obspy.core import read
import numpy as np
import matplotlib.pyplot as plt
import struct
import os
import glob
def qnxsei2obspy(file):
    st = read(file)
    for tr in st:
        data = tr.data
        for i in range(0, len(data) -1):
            data[i] = swap32(data[i])
        data[-1] = 0
        #data = np.where(data < 360000, data, 0)
        tr.data = data
        tr.detrend(type='constant')
        #tr.plot()
    raw_input("Press Enter to continue...")
    st.plot()
    dirbase, ext = os.path.splitext(file)
    outfile = dirbase + '.mseed'
    #st.write(outfile, 'MSEED')
         
         

def swap32(i):
    return struct.unpack("<i", struct.pack(">i", i))[0]

def main():
        # read in the wav file name from the command line
    seisantop = 'D:\\roxio\\seisan_from_qnx'
    os.chdir('D:\\roxio\\seislog\\2000\\03')
    print os.getcwd()
    #years = glob.glob('[0-9]' * 4)
    #for year in years:
        #os.chdir(year)
        #print os.getcwd()
        #months = glob.glob('[0-9]' * 2)
        #for month in months:
        #		os.chdir(month)
    print os.getcwd()
	#		seislogfiles = glob.glob('*T.MVO*')
    seislogfiles = glob.glob('*S.MVO*')
    print("Found %d Seislog files" % (len(seislogfiles)))
    for seislogfile in seislogfiles:
        qnxsei2obspy(seislogfile)
        #os.chdir('..')
        #os.chdir('..')


if __name__ == "__main__":
	main()
