import os
import pandas as pd
import datetime as dt
import glob
startdate = dt.date(1995,7,27)
enddate = dt.date(2001,3,1)
delta = dt.timedelta(days=1)
date0=startdate
wavfilesperday = dict()
counterfilesperday = dict()
seisandb = "ASNE_"
seisanpath = os.path.join(os.environ['SEISAN_TOP'], 'WAV', seisandb)
while date0 <= enddate:
    yyyymmdir = os.path.join(seisanpath, date0.strftime("%Y"), date0.strftime("%m") )
    fileroot = date0.strftime("%y%m%d")
    keyname = date0.strftime("%Y-%m-%d")
    globstring = os.path.join(yyyymmdir, fileroot) + '*.mseed'
    #print(globstring)
    allmseedfiles = sorted(glob.glob(globstring))
    wavfilesperday[keyname] = len(allmseedfiles)
    counterfilesperday[keyname] = len(allmseedfiles)
    if len(allmseedfiles) > 0:
        lastfile = allmseedfiles[-1]
        fname = os.path.basename(lastfile)
        CC = fname[6:8]
        if CC[0] > '9':
            counter = (ord(CC[0]) - 65 + 10) * 36
        else:
            counter = (ord(CC[0]) - 48) * 36
        if CC[1] > '9':
            counter += ord(CC[1]) - 65 + 10
        else:
            counter += ord(CC[1]) - 48
        counterfilesperday[keyname] = counter + 1 # since 00 is file 1
    date0 += delta
    print(keyname, wavfilesperday[keyname], counterfilesperday[keyname], counterfilesperday[keyname] - wavfilesperday[keyname])
df = pd.DataFrame.from_dict(wavfilesperday, orient='index')
df.to_csv('wavfilesperday.csv')
df2 = pd.DataFrame.from_dict(counterfilesperday, orient='index')
df2.to_csv('counterfilesperday.csv')
