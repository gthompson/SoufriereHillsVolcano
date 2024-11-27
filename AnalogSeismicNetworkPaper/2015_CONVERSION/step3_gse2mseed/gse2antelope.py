# Walk through the directory structure GSE_TOP_DIR/YYYY/M
# Convert all the GSE files found into Miniseed format in ANTELOPE_TOP_DIR/YYYY/JJJ directories suitable for Antelope
# You must have ObsPy installed
#
# The environment variables GSE_TOP_DIR and ANTELOPE_TOP_DIR need to be set

# Glenn Thompson 2015 University of South Florida
from obspy import read
import os.path
import numpy as np
import datetime
import sys
from os import walk
import os
import glob
#sys.path.insert(0, '/home/t/thompsong/Desktop')

def gse2antelope(thisfile, antelopedirectoryroot):
	network = 'MV'
	print thisfile
	#continue

	#seisandir = '/raid/data/seisan/WAV/MVOE_'
	#year = '2002'
	#month = '01'
	#julday = '04'
	#seisanbasename = '2002-01-04-1424-06S.MVO___015'
	#st = read(seisandir + '/' + year + '/' + month + '/' + seisanbasename)
	#ntraces = np.size(st)

	try:
		st = read(thisfile)
	except:
		print "Cannot read " + thisfile
		return
	
	for tr in st:
	
		# remove the final sample because reading an MVO Seisan file into ObsPy seems to cause the last sample to be random and large
		#length_data = np.size(tr.data)
		#tr.data = np.delete(tr.data,[length_data - 1])
	
		# check for all 0
		a = tr.copy()
		a = a.detrend(type='demean')
		anyData = np.count_nonzero(a.data)
		if anyData==0:
			st.remove(tr)
			continue
		
		# fix the network, channel and location
		#print " "
		#print tr.stats
		tr.stats['network']=network
		sta = tr.stats['station'].strip()
		chan = 'SH' + tr.stats['channel'].strip()
		'''
		if chan=='PRS' or chan=='APS':
			chan='BDF'
		else:
			if chan[0]=='A':
				if tr.stats['location'] == 'J':
					bandcode = 'S'
				#else:
				#	bandcode = 'B'
				
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
		'''
		tr.stats['location']='--'
		#continue
	
		# write out to SDS structure - but this only works for full day files
		#filepath = "/raid/data/SDS/" + year + "/" + network + "/" + tr.stats['station'] + "/" + tr.stats['channel'] + "D"
		#filename = tr.id + ".D." + year + "." + julday
		#tr.write(filename + ".MSEED", format="MSEED")
	
		# get sta, chan
		#sta = tr.stats['station']
		#chan = tr.stats['channel']
	
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
		antelopedir = os.path.join(antelopedirectoryroot, yyyy, jjj)
		if not os.path.exists(antelopedir):
	    		os.makedirs(antelopedir)
		#antelopefile = sta + "." + chan + "." + yyyy + ":" + jjj + ":" + hh + ":" + ii + ":" + ss
		antelopefile = sta + "." + chan + "." + yyyy + "_" + jjj + "_" + hh + "_" + ii + "_" + ss
		antelopefullpath = os.path.join(antelopedir, antelopefile)
		print antelopefullpath
		tr.write(antelopefullpath, format="MSEED")
		#input("Press Enter to continue...")			
	# write SAC or MSEED files to a new Seisan database
	#(filename,ext) = os.path.splitext(seisanbasename)
	#st.write(filename + ".MSEED", format="MSEED")
	
	



def main():
        # read in the wav file name from the command line
        #fileroot = sys.argv[1]
	#fileroot = os.path.join('D:', 'waveforms')
	#print fileroot
        #print os.getcwd()
	#os.chdir(fileroot)
	GSE_TOP = os.environ("GSE_TOP")
	ANTELOPE_TOP = os.environ("ANTELOPE_TOP")
	print os.getcwd()
        os.chdir(GSE_TOP)
	years = glob.glob('[0-9]' * 4)
	for year in years:
		os.chdir(year)
		print os.getcwd()
		months = glob.glob('[0-9]' * 2)
		for month in months:
			os.chdir(month)
			print os.getcwd()
			gsefiles = glob.glob('*.gse')
			print("Found %d GSE files" % (len(gsefiles)))
			for gsefile in gsefiles:
				gse2antelope(gsefile, ANTELOPE_TOP)
			os.chdir('..')
		os.chdir('..')

if __name__ == "__main__":
	main()
