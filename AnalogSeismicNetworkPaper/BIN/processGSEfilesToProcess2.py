from obspy.core import read
import os
infile = 'c:/Users/thompsong/SUDS/GSEfilesToProcess2.txt'
f = open(infile, 'r')
linelist = f.readlines()
for line in linelist:
    gsefile = 'C:/Seismo/WAV/ASNE_/' + line[:-1] + '.gse'
    mseedfile = 'C:/Seismo/WAV/ASNE_/' + line[:-1] + '.mseed'
    if os.path.exists(gsefile):
        print('%s yes' % gsefile)
        if not os.path.exists(mseedfile):
            st = read(gsefile)
            print('Saving ' + mseedfile)
            st.write(mseedfile)
    else:
        print('%s no' % gsefile)
f.close()
