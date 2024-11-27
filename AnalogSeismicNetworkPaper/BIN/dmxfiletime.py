import os
import sys
if len(sys.argv) > 1:
    dmxfile = sys.argv[1]
else:
    dmxfile =  'c:/seismo/WAV/ASNE_/1995/12/95123109.dmx'

prog = 'c:/users/thompsong/src/SUDS2SEISAN/winsuds/bin/sud2asc.exe'
ascfile = dmxfile[-12:-4] + '.asc'
#print(prog)
#print(dmxfile)
#print(ascfile)
os.system('%s %s > %s' % (prog, dmxfile, ascfile))
if os.path.exists(ascfile):
    f = open(ascfile, 'r')
    all_lines = f.readlines()
    for line in all_lines:
        if line.find('effective') > -1:
            #print(line);
            MM = line[0:2]
            DD = line[3:5]
            YY = line[6:8]
            hh = line[9:11]
            mi = line[12:14]
            ss = line[15:21]
            print(YY + MM + DD + hh + mi + ss)
            break;
    f.close()
