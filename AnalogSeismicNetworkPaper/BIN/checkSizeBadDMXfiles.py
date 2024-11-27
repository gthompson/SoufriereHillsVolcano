import os, sys
f = open('badDMXfileList.txt','r')
lines = f.readlines()
nonzero = 0
numzero = 0
numnonexisting = 0
numgse = 0
nummseed = 0
numboth = 0
for line in lines:
    line = line[:-1]
    if not os.path.exists(line):
        line = 'c:/seismo/wav/asne_/' + line
    if os.path.exists(line):
        statinfo = os.stat(line)
        if statinfo.st_size > 0:
            nonzero += 1
            gsefile = line[:-4] + '.gse'
            mseedfile = line[:-4] + '.mseed'
            if os.path.exists(gsefile):
                if os.path.exists(mseedfile):
                    numboth += 1
                    statinfo = os.stat(mseedfile)
                    if statinfo.st_size == 0:
                        print('bad MSEED file: ' + mseedfile)
                else:
                    numgse += 1
            elif os.path.exists(mseedfile):
                nummseed += 1
                statinfo = os.stat(mseedfile)
                if statinfo.st_size == 0:
                    print('bad MSEED file and no GSE file: ' + mseedfile)
            else:
                print('No GSE or MSEED file: ' + line)
        else:
            numzero += 1
            print('Zero size: ' + line)
    else:
        numnonexisting += 1

f.close()
print('zero-sized files: %d\nnon-zero-sized files: %d\nnot existing: %d\n' % (numzero, nonzero, numnonexisting))
print('of non-zero-sized files: %d have GSE files, %d have MSEED files, %d have both\n' % (numgse, nummseed, numboth))

print('FIRST PART DONE')

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
        statinfo = os.stat(mseedfile)
        if statinfo.st_size == 0:
            print('bad MSEED file from GSE files to process: ' + mseedfile)
    else:
        print('%s no' % gsefile)
f.close()

