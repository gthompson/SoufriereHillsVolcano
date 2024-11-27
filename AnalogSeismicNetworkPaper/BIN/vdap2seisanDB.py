#!/usr/bin/env python
# coding: utf-8

# # SUDS to MiniSEED converter
# This Python3 notebook demonstrates how to convert a PC-SUDS seismic waveform data file into Miniseed
# Requirements:
# 1. Win-SUDS Utilities (for programs demux.exe, irig.exe and sud2sac.exe).
# 2. ObsPy (to read SAC files and combine them into a Miniseed file).
# 3. Seisan (for programs makerea, dirf, autoreg & hypnor)
#
# (c) Glenn Thompson, November 2019 thompsong@mail.usf.edu

# Need 2 command line arguments for AGENCY and SEISANDBNAME
import sys
if len(sys.argv) != 3:
    print("Usage:\n\t %s 3-letter-agency-code 5-letter-seisan-db-name\n" % sys.argv[0])
    exit(0)

agency = sys.argv[1]
seisanDBname = sys.argv[2]
if len(agency) < 3:
    agency = agency.ljust(3)
if len(agency) > 3:
    agency = agency[0:3]
if len(seisanDBname) < 5:
    seisanDBname = seisanDBname + "_" * (5-len(seisanDBname))
if len(seisanDBname) > 5:
    seisanDBname = seisanDBname[0:5]

# Do all library imports first
import obspy.core as op
import glob
import matplotlib.pyplot as plt
import os

# Set up all paths that will be used
cwd = os.getcwd()
WINSUDSPATH = os.path.join(cwd,'WINAPPS','winsuds','bin')
demux = os.path.join(WINSUDSPATH, 'demux.exe')
irig = os.path.join(WINSUDSPATH, 'irig.exe')
sud2sac = os.path.join(WINSUDSPATH, 'sud2sac.exe')
sud2msed = os.path.join(WINSUDSPATH, 'sud2msed.exe')
SEISAN_TOP = os.environ['SEISAN_TOP']
print('Paths setup okay')

# Define all functions 

def findLatestFile(dirpath):
    print("Looking for last file in " + dirpath)
    lists = os.listdir(dirpath)
    if len(lists)>0:
        lists.sort(key=lambda fn:os.path.getmtime(os.path.join(dirpath, fn)))
        return os.path.join(dirpath,lists[-1])
    else:
        return ""

def displayFile(dirpath):
    if os.path.exists(dirpath):
        fptr = open(dirpath, 'r')
        str = fptr.read()
        print(str)
    else:
        print(dirpath + ' not found')

def appendPicks(hypnoroutfile, sfilepath):
    if os.path.exists(hypnoroutfile) and os.path.exists(sfilepath):
        pass
    else:
        print('Cannot combine files, one of the inputs does not exist')
    print("Appending " + hypnoroutfile + " to " + sfilepath)
    newlines = list()
    with open(sfilepath, "r") as fs:
        for line in fs:
            #print(line[-2])
            if line[-2] == '7':
                print("- Inserting " + hypnoroutfile)
                #fh = open(hypnoroutfile, 'r')
                with open(hypnoroutfile, 'r') as fh:
                    for hnline in fh:
                        if hnline[-2] != '1':
                            newlines.append(hnline)
                #str = fh.read()
                #fh.close()
                #newlines.append(str)
            else:
                newlines.append(line)
    fs2 = open(sfilepath, "w")
    for newline in newlines:
        fs2.write(newline)
    fs2.close()
    return


def HSUMNOR(PUNfile, agency):
    # Based on hsumnor.for in Seisan (which did not work, but this does!)
    # I should add something for the high accuracy lines too, as this only records origin to nearest tenth of second
    outstr = ' ' * 79 + '1'
    if os.path.exists(PUNfile):
        #print('0123456789'*8)
        #displayFile(PUNfile)
        with open(PUNfile, 'r') as f1:
            for line in f1:
                pass
            
        str = line
        yy = float(str[0:2].strip())
        mm = float(str[2:4].strip())
        dd = float(str[4:6].strip())
        hr = float(str[7:9].strip())
        mi = float(str[9:11].strip())
        sec = float(str[12:17].strip())
        lat = float(str[17:20].strip())
        mlat = float(str[21:26].strip())
        lon = float(str[27:30].strip())
        mlon = float(str[31:36].strip())
        depth = float(str[37:43].strip())
        magstr = str[45:49].strip()
        if magstr == "":
            magstr = " " * 3
        else:
            magstr = "%3.1f" % float(magstr)
        nphase = float(str[50:53].strip())
        rms = float(str[62:66].strip())
        dlat = float(lat) + float(mlat)/60.0
        dlon = float(lon) + float(mlon)/60.0
        if yy < 80:
            cen = 20
        else:
            cen = 19
        newstr = ' %2d%02d %02d%02d %02d%02d %4.1f L ' % (cen,yy,mm,dd,hr,mi,sec)    
        newstr = newstr + '%7.3f%8.3f%5.1f  %s%3d%4.1f %s' % (dlat, dlon, depth, agency, nphase, rms, magstr)
        outstr = newstr + outstr[len(newstr):]
    return outstr


def insertHYPO71Summary(PUNfile, sfilepath, agency, seisanDBname):
    if os.path.exists(PUNfile) and os.path.exists(sfilepath):
        pass
    else:
        print('Cannot combine files, one of the inputs does not exist')
    print("Converting " + PUNfile + " and inserting into " + sfilepath)
    newlines = list()
    try:
        hypo71sumstr = HSUMNOR(PUNfile, agency)
    except:
        print('HSUMNOR failed on %s. Adding to bad PUN list.' % PUNfile)
        fbad = open("badPUNfileList.txt","a+")
        fbad.write(PUNfile + "\n")
        fbad.close()
        return
    if hypo71sumstr != "":
        newlines.append(hypo71sumstr + "\n")
    with open(sfilepath, "r") as fs:
        for line in fs:
            newlines.append(line)
    fs2 = open(sfilepath, "w")
    for newline in newlines:
        fs2.write(newline)
    fs2.close()
    os.system("move %s %s" % (PUNfile, os.path.join(os.environ['SEISAN_TOP'], "WAV", seisanDBname, yyyy, mm, os.path.basename(PUNfile) )) )
    return

def processPHAfile(PHAfile, SFILE, seisanDBname, yyyy):
    if os.path.exists(PHAfile):
        yy = yyyy[2:4]
        century = '19'
        if yy < '80':
            century = '20'
        fptr = open('hypnor_wrapper.txt','w')
        fptr.write(PHAfile + '\n')
        fptr.write(century + '\n')
        fptr.close()
        #displayFile('hypnor_wrapper.txt')
        os.system('hypnor < hypnor_wrapper.txt')
        if os.path.exists('hypnor.out'):
            print('Running HYPNOR on ' + PHAfile)
            appendPicks('hypnor.out', SFILE)
            os.system("move %s %s" % (PHAfile, os.path.join(os.environ['SEISAN_TOP'], "WAV", seisanDBname, yyyy, mm, os.path.basename(PHAfile)) ) )
            print(' - Success')
        else:
            print(' - Failed')
            print('appendPicks failed on %s. Adding to bad PHA list.' % PHAfile)
            fbad = open("badPHAfileList.txt","a+")
            fbad.write(PHAfile + "\n")
            fbad.close()

def findMatchingSfile(file, seisanDBname, yyyy, mm):
    # find SFILE related to this PHAfile
    MSEEDfile = os.path.basename(file)[0:-4] + '.mseed'
    yy = yyyy[2:4]
    dd = MSEEDfile[4:6]
    SFILE = ''
    possible_sfiles = glob.glob(os.path.join(os.environ['SEISAN_TOP'], 'REA', seisanDBname, yyyy, mm, dd + '*') ) 
    for possible_sfile in possible_sfiles:
        fs = open(possible_sfile, 'r')
        sfile_contents = fs.read()
        fs.close()
        if sfile_contents.find(MSEEDfile) != -1:
            SFILE = possible_sfile
    return SFILE

#### MAIN PROGRAM STARTS HERE #############
# Loop over year directories
yyyydirs = glob.glob('[1-2][0-9][0-9][0-9]')
for yyyydir in yyyydirs:
    yyyy = os.path.basename(yyyydir)
    mmdirs = glob.glob(os.path.join(yyyydir,'[0-1][0-9]'))
    for mmdir in mmdirs:
        mm = os.path.basename(mmdir)

        # make directories if needed
        if not os.path.exists(os.path.join(os.environ['SEISAN_TOP'], "WAV", seisanDBname, yyyy, mm)):
            fptr = open('makerea_wrapper.txt','w')
            fptr.write(seisanDBname + '\n')
            fptr.write(yyyy + mm + '\n')
            fptr.write('\n')
            fptr.write('BOTH\n')
            fptr.close()
            os.system('makerea < makerea_wrapper.txt')

        WVMfiles = glob.glob(os.path.join(mmdir, '????????.WVM'))
        print('Found %d WVM files in %s' % (len(WVMfiles), mmdir))
        for WVMfile in WVMfiles:
            SUDSdirbase = WVMfile[0:-4]
            SUDSbasename = os.path.basename(SUDSdirbase)
            MSEEDfile = SUDSdirbase + '.mseed' # produced by recombining SAC file
            MSEEDdbfile = os.path.join(os.environ['SEISAN_TOP'], "WAV", seisanDBname, yyyy, mm, os.path.basename(MSEEDfile) )
            #print(MSEEDdbfile + "\n")
            if os.path.exists(MSEEDdbfile):
                #print(MSEEDdbfile + ' already exists. Not processing ' + WVMfile)
                continue
            else:
                print(MSEEDdbfile + ' not found. Processing ' + WVMfile)
                #continue
            DMXfile = SUDSdirbase + '.dmx'
            if os.path.exists(DMXfile):
                pass
            else:
                print('Demultiplexing ' + WVMfile)
                os.system(demux + ' ' + WVMfile)
                print('Done')
                if (os.path.exists(DMXfile)):
                    os.system("move %s %s" % (WVMfile, os.path.join(os.environ['SEISAN_TOP'], "WAV", seisanDBname, yyyy, mm, os.path.basename(WVMfile)) ) )
                else:
                    print('Demultiplexing seems to have failed as ' + DMXfile + ' not created. Adding to list of bad WVM files.')
                    fbad = open('badWVMfileList.txt','a+')
                    fbad.write(WVMfile + '\n')
                    fbad.close()

        DMXfiles = glob.glob(os.path.join(mmdir, '????????.dmx'))
        print('Found %d dmx files in %s' % (len(DMXfiles), mmdir))
        for DMXfile in DMXfiles:
            SUDSdirbase = DMXfile[0:-4]
            SUDSbasename = os.path.basename(SUDSdirbase)
            SACdirbase = SUDSdirbase + '.sac' # produced  by sud2sac.exe
            MSEEDfile = SUDSdirbase + '.mseed' # produced by recombining SAC file
            PHAfile = SUDSdirbase + '.PHA' # this might exist if HYPO71 was run to locate the event
            PUNfile = SUDSdirbase + '.PUN' # this might exist if HYPO71 was run and generated a hypocenter
            #yyyy = SUDSdirbase[0:4]
            #mm = SUDSdirbase[5:7]
            MSEEDdbfile = os.path.join(os.environ['SEISAN_TOP'], "WAV", seisanDBname, yyyy, mm, os.path.basename(MSEEDfile) )
            if os.path.exists(MSEEDdbfile):
                #print(MSEEDdbfile + ' already exists. Not processing ' + DMXfile)
                continue
            else:
                print(MSEEDdbfile + ' not found. Processing ' + DMXfile)
                #continue

            print('Time correcting ' + DMXfile)
            os.system(irig + ' ' + DMXfile)
    
            # This produces one file per trace
            print('Converting ' + DMXfile + ' to SAC files')
            try:
                os.system(sud2sac + ' ' + DMXfile)
            except:
                print('Crashed on converting')

            sacfilelist = glob.glob(SACdirbase + '-???')
            if len(sacfilelist) == 0:
                # DMX file is likely a bad file. Check for WVM file and demultiplex it
                if os.path.exists(WVMfile):
                    os.system(demux + ' ' + WVMfile)
                    os.system(sud2sac + ' ' + DMXfile)
                    os.system(irig + ' ' + DMXfile)
            print('Done')

            # Now merge the SAC files into a single valid Miniseed file    
            st = op.Stream()
            if len(sacfilelist) > 0:
                for sacfile in sacfilelist:
                    print('Combining ' + sacfile + ' into Miniseed file ' + MSEEDfile)
                    try:
                        tr = op.read(sacfile)
                        st = st + tr
                    except:
                        pass
                if len(st) == 0:
                    print('No good SAC files')
                    continue
                st.write(MSEEDfile)
                print('Done')
                #os.remove(SACdirbase + '-???')
                os.system("del " + SACdirbase + '-???')
            else:
                print('No SAC files found matching ' + SACdirbase + '*. Must have been bad DMX file ' + DMXfile)
                fbad = open('badDMXfileList.txt','a+')
                fbad.write(DMXfile + '\n')
                fbad.close()

            if not os.path.exists(MSEEDfile):
                print('Failed to create a Miniseed file for %s' % MSEEDfile)
                continue

            # autoregister the MSEEDfile if it does not already appear in an SFILE
            MSEEDbase = os.path.basename(MSEEDfile)
            os.system('move ' + MSEEDfile + ' ' + MSEEDbase)
            SFILE = findMatchingSfile(PHAfile, seisanDBname, yyyy, mm)
            if not SFILE:
                os.system('dirf ' + MSEEDbase)
                fptr = open('autoreg_wrapper.txt','w')
                fptr.write('L\n')
                fptr.write('m\n')
                fptr.write(seisanDBname + '\n')
                fptr.write('gt\n')
                fptr.close()
                os.system('autoreg < autoreg_wrapper.txt')


            ## try to guess name of the Sfile that was just created
            #SFILE = findLatestFile(os.path.join(os.environ['SEISAN_TOP'],'REA',seisanDBname,yyyy,mm))
            #if len(SFILE)==0:
            #    print('no file found')
            #    exit(0)
            #print(SFILE)
            ##displayFile(SFILE)

            # If BUDSPICK generated a phase picks file (file extension *.PHA), add that into the S-file using the Seisan program HYPNOR
            #if os.path.exists(PHAfile):
            #    processPHAfile(PHAfile, SFILE, seisanDBname, yyyy)
            #else:
            #    print('No PHA file for ' + SUDSdirbase)

            #if os.path.exists(PUNfile):
            #    insertHYPO71Summary(PUNfile, SFILE, agency, seisanDBname)

            # display the final S-file    
            displayFile(SFILE)
        
        # Process any other PHA or PUN files, in case they had no SFILE
        PHAfiles = glob.glob(os.path.join(mmdir, '????????.PHA'))
        print('Found %d PHA files in %s' % (len(PHAfiles), mmdir))
        for PHAfile in PHAfiles:
            SFILE = findMatchingSfile(PHAfile, seisanDBname, yyyy, mm)
            if SFILE:
                # check if PHA file already in WAV directory. If not, run processPHAfile
                PHAdbfile = os.path.join(os.environ['SEISAN_TOP'], 'WAV', seisanDBname, yyyy, mm, os.path.basename(PHAfile) )
                if os.path.exists(PHAdbfile):
                    pass
                else:
                    processPHAfile(PHAfile, SFILE, seisanDBname, yyyy)
                    displayFile(SFILE)
            else:
                print('No SFILE exists yet matching %s. Ignoring for now.' % PHAfile)
            
        PUNfiles = glob.glob(os.path.join(mmdir, '????????.PUN'))
        print('Found %d PUN files in %s' % (len(PUNfiles), mmdir))
        for PUNfile in PUNfiles:
            SFILE = findMatchingSfile(PUNfile, seisanDBname, yyyy, mm)
            if SFILE:
                # check if PUN file already in WAV directory. If not, run insertHYPO71Summary
                PUNdbfile = os.path.join(os.environ['SEISAN_TOP'], 'WAV', seisanDBname, yyyy, mm, os.path.basename(PUNfile) )
                if os.path.exists(PUNdbfile):
                    pass
                else:
                    insertHYPO71Summary(PUNfile, SFILE, agency, seisanDBname)
                    displayFile(SFILE)
            else:
                print('No SFILE exists yet matching %s. Ignoring for now.' % PUNfile)
