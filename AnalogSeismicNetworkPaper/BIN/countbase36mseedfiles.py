import os
import pandas as pd
rootDir = 'C:/seismo/WAV/ASNE_'
df = pd.DataFrame(columns=['date', 'maxcounter'])
for dirName, subdirList, fileList in os.walk(rootDir):
    print('Found directory: %s' % dirName)
    for fname in fileList:
        if fname[-5:] == 'mseed':

            YY = fname[0:2]
            if YY[0] == '9':
                YYYY = '19' + YY
            else:
                YYYY = '20' + YY
            MM = fname[2:4]
            DD = fname[4:6]
            CC = fname[6:8]
            filedate = YYYY + MM + DD
            if CC[0] > '9':
                counter = (ord(CC[0]) - 65 + 10) * 36
            else:
                counter = (ord(CC[0]) - 48) * 36
            if CC[1] > '9':
                counter += ord(CC[1]) - 65 + 10
            else:
                counter += ord(CC[1]) - 48

            print('\t%s %d' % (fname, counter))
            if filedate in df:
                if counter > df.loc[filedate,'maxcounter']:
                    df.ix[filedate] = counter
            else:
                df.ix[filedate] = counter

df.to_csv('mseedfilesperday2.csv')
