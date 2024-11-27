#!/shares/gly4930.018f19.gt/miniconda3/bin/python
## Allow student to check out 1 day of Seisan WAV files from the MVO ASNE db
# Requirements:
# 1. Each student can only check-out 1 day at a time
# 2. When a day is checked-out, it is locked
# 3. When a day is complete, it is check back in & the student is allowed to check out a new day
# 4. The status for each day needs to be tracked. Available. Checked-out (and assigned to who). And finished (and who checked-in by).
# 5. The status for each student needs to be tracked. Should be no day, or day (and which day). Also days done.
# 6. A student should also be able to release a day, to make it available again.
#
# Data types:
# daystatus['yyyymmdd'] = 'Available'|'student-name'|'Finished'
# methods: check-out, check-in, release
# These processes should all send an email to me and CC to everyone else
# studentstatus['student-name'] = 'None'|'yyyymmdd'
# methods: 
#
# Also need methods to check how many days and how many S-files each student has processed.
# This needs to be based on author field in the S-file.
# Do this as a weekly update.
# 
# Several levels: tin (1 k), bronze (2.5k), silver (5k), gold (10k), platinum (20k)
# Needs to work with the second way of doing the processing too
import sqlite3
from sqlite3 import Error
import pandas as pd 
import datetime as dt
import os
import sys
import glob
import re
#SEISAN_TOP = '/shares/gly4930.018f19.gt/seisan'

def count_sfiles(yyyymmdd):
    #globstr = seisan_top + '/REA/ASNE_/' + yyyymmdd[0:4] + '/' + yyyymmdd[4:6] + '/' + yyyymmdd[6:8]
    globstr = '../../REA/ASNE_/' + yyyymmdd[0:4] + '/' + yyyymmdd[4:6] + '/' + yyyymmdd[6:8]
    #print(globstr)
    sfiles = sorted(glob.glob(globstr + '*[LRD].S*'))
    num_sfiles = len(sfiles)
    #print(num_sfiles, ' S-files found for ',yyyymmdd)
    return num_sfiles

def count_wavfiles(yyyymmdd):
    #globstr = seisan_top + '/WAV/ASNE_/' + yyyymmdd[0:4] + '/' + yyyymmdd[4:6] + '/' + yyyymmdd[0:4] + '-' + yyyymmdd[4:6] + '-' + yyyymmdd[6:8]
    globstr = '../../WAV/ASNE_/' + yyyymmdd[0:4] + '/' + yyyymmdd[4:6] + '/' + yyyymmdd[0:4] + '-' + yyyymmdd[4:6] + '-' + yyyymmdd[6:8]
    #print(globstr)
    wavfiles = sorted(glob.glob(globstr + '*0??'))
    num_wavfiles = len(wavfiles)
    #print(num_wavfiles, ' WAV-files found for ',yyyymmdd)
    return num_wavfiles

def sfilename(wavfile):
    spath = os.path.dirname(wavfile)
    spath = spath.replace('WAV','REA',1)
    wavbase = os.path.basename(wavfile)
    #print(wavbase)
    yyyy = wavbase[0:4]
    mm = wavbase[5:7]
    dd = wavbase[8:10]
    HH = wavbase[11:13]
    MI = wavbase[13:15]
    SS = wavbase[16:18]
    sbase = dd + "-" + HH + MI + "-" + SS + "L.S" + yyyy + mm
    sbaseminute = dd + "-" + HH + MI + "-00L.S" + yyyy + mm
    sfile = spath + "/" + sbase
    sfileminute = spath + "/" + sbaseminute
    return sfile, sfileminute

def dirf(yyyymmdd, netid1):
    #globstr = seisan_top + '/WAV/ASNE_/' + yyyymmdd[0:4] + '/' + yyyymmdd[4:6] + '/' + yyyymmdd[0:4] + '-' + yyyymmdd[4:6] + '-' + yyyymmdd[6:8]
    globstr = '../../WAV/ASNE_/' + yyyymmdd[0:4] + '/' + yyyymmdd[4:6] + '/' + yyyymmdd[0:4] + '-' + yyyymmdd[4:6] + '-' + yyyymmdd[6:8]
    #print(globstr)
    wavfiles = sorted(glob.glob(globstr + '*0??'))
    num_wavfiles = len(wavfiles)
    dirflist = list()
    #print(num_wavfiles, ' WAV-files found for ',yyyymmdd)
    for wavfile in wavfiles:
        sfile, sfileminute = sfilename(wavfile) # assumes L local event. 
        if os.path.exists(sfile) or os.path.exists(sfileminute):
            # event already processed
            #print("WAV %s already processed. %s exists" % (wavfile, sfile))
            pass
        else:
            # Could still be a regional or teleseism
            sfileR = sfile.replace('L.S','R.S')
            sfileminuteR = sfileminute.replace('L.S','R.S')
            if os.path.exists(sfileR) or os.path.exists(sfileminuteR):
                pass
            else:
                sfileD = sfile.replace('L.S','D.S')
                sfileminuteD = sfileminute.replace('L.S','D.S')
                if os.path.exists(sfileD) or os.path.exists(sfileminuteD):
                    pass
                else:
                    dirflist.append(wavfile)
    listnum = 1
    f = open(seisan_top + "/" + "WOR" + "/" + netid1 + "/filenr.lis", "w+")
    for dirfelement in dirflist:
        print(" #%3d  %s" % (listnum, dirfelement))
        f.write(" #%3d  %s\n" % (listnum, dirfelement))
        listnum += 1
    f.close()
    return len(dirflist)

def initialize_db(db_file, netid1):
    if os.path.exists(db_file):
        if netid1 == 'thompsong':
            print(db_file,' exists. Overwrite (y/n) ? ')
            yn = input('')
            if yn == "y":
                os.remove(db_file)
        else:
            print('Not authorized to delete existing database')
            return

    """ create a database connection to a SQLite database """
    dbptr = None
    try:
        dbptr = sqlite3.connect(db_file)
        cursor = dbptr.cursor()
        print('Creating USER table')
        cursor.execute('''
            CREATE TABLE user(netid TEXT PRIMARY KEY, 
                current_day INTEGER, days_finished INTEGER)
            ''')
        print('Creating DAY_STATUS table')
        cursor.execute('''
            CREATE TABLE day_status(ymd INTEGER PRIMARY KEY, status TEXT,
                netid TEXT)
            ''')


        # initialize the day_status table 
        start_date = dt.datetime(1995,7,18)
        end_date = dt.datetime(1996,12,31)
        days_finished_tshea = 0
        for thisdate in pd.date_range(start_date, end_date):
            yyyymmdd = thisdate.strftime("%Y%m%d")
            netid1 = ''
            num_sfiles = count_sfiles(yyyymmdd)
            num_wavfiles = count_wavfiles(yyyymmdd)
            print('Adding ',yyyymmdd,' to DAY_STATUS')
            if num_sfiles == 0:
                if num_wavfiles > 0:
                    status1 = 'Available'
                else:
                    status1 = 'GlennCopyWAVS'
            elif num_sfiles >= num_wavfiles:
                status1 = 'Done'
                netid1 = 'tshea'
                days_finished_tshea += 1
            else:
                status1 = 'Incomplete'
                netid1 = 'tshea'
            cursor.execute('''INSERT INTO day_status(ymd, status, netid)
                VALUES(:ymd,:status,:netid)''',
                {'ymd':yyyymmdd, 'status':status1, 'netid':netid1})

        # initialize the user table 
        userlist = ['thompsong', 'mahsaafra', 'cesil', 'dinovav', 'jferreira', 'dbgraybeal', 'mshastings1', 'bcmack', 'lmartins1', 'elhammoslemi', 'andries', 'aprevatt', 'taha4', 'tshea', 'candicesmith']
        for thisuser in userlist:    
            print('Adding ',thisuser,' to USERS')
            days_done = 0
            if thisuser == 'tshea':
                days_done = days_finished_tshea
            cursor.execute('''INSERT INTO user(netid, current_day, days_finished) 
                VALUES(:netid,:current_day,:days_finished)''',
                {'netid':thisuser, 'current_day':0, 'days_finished':days_done})

        print('Saving new database')
        dbptr.commit() 

    except Error as e:
        print(e)
 
    finally:
        if dbptr:
            dbptr.close()


def list_available_days(db_file):
    dbptr = None
    try:
        dbptr = sqlite3.connect(db_file)
        cursor = dbptr.cursor()
        print('Date   #WAVs  #S-files  Status    NetID')    
        start_date = dt.datetime(1995,7,18)
        end_date = dt.datetime(1996,12,31)
        for thisdate in pd.date_range(start_date, end_date):
            yyyymmdd = thisdate.strftime("%Y%m%d")
            num_sfiles = count_sfiles(yyyymmdd)
            num_wavfiles = count_wavfiles(yyyymmdd)

            cursor.execute('''SELECT status,netid FROM day_status WHERE ymd=?''',(yyyymmdd,)) 
            status1_tup = cursor.fetchone() #retrieve the first row
            status1 = status1_tup[0]
            netid1 = status1_tup[1]
            if num_sfiles >= num_wavfiles and status1 != "Done":
                status1 = "Done"
                #fix_daystatus_table(db_file, "thompsong", yyyymmdd, status1, netid1 )
                if netid1 is not None:
                    cursor.execute('''SELECT days_finished FROM user WHERE netid=?''',(netid1,))
                    daysfinished_tup = cursor.fetchone()
                    if daysfinished_tup is not None:
                        #print(netid1,' ',daysfinished_tup)
                        daysfinished = daysfinished_tup[0] 
                        daysfinished += 1
                        #print(netid1,' days finished = ',daysfinished)
                        #fix_user_table(db_file, "thompsong", netid1, yyyymmdd, daysfinished )

            if status1 == 'Available' or status1 == 'Incomplete':
                print('%s, %4d, %4d, %s, %s' % (yyyymmdd, num_wavfiles, num_sfiles, status1, netid1))

        print('\nListing complete')
    except Error as e:
        print(e)
 
    finally:
        if dbptr:
            dbptr.close()

# CHECKOUT DAY
# If day_status:status = 'Available' for ymd=yyyymmdd and user:current_day = 0, 
# then set status = netid and current_day = yyyymmdd
def checkout_day(db_file, netid1):
    try:
        yyyy = int(input("Enter year as 4 digits, e.g. 1996 :?  "))
    except ValueError:
        print("not an integer")
        return
    try:
        mm = int(input("Enter month as 2 digits, e.g. 09 :?  "))
    except ValueError:
        print("not an integer")
        return
    try:
        dd = int(input("Enter day as 2 digits, e.g. 28 :?  "))
    except ValueError:
        print("not an integer")
        return

    correctDate = None
    try:
        newDate = dt.datetime(yyyy,mm,dd)
        correctDate = True
        yyyymmdd = newDate.strftime("%Y%m%d")
    except ValueError:
        correctDate = False
        print("The date you entered is invalid")
        return
    if yyyymmdd < "19950718" or yyyymmdd > "19961231":
        print('Date needs to be between 19950718 and 19961231')
        return

    dbptr = None
    try:
        dbptr = sqlite3.connect(db_file)
        cursor = dbptr.cursor()
    
        # Check if this day is available
        cursor.execute('''SELECT status FROM day_status WHERE ymd=?''',(yyyymmdd,)) 
        try:
            status1_tup = cursor.fetchone() #retrieve the first row
            status1 = status1_tup[0]
        except Error as e:
            print('No matching row found')
            print(e)
            return
        for row in cursor:
            print(row)
        if status1 == 'Available' or status1 == 'Incomplete': 
            # Check if this user has a day checked-out already
            cursor.execute('''SELECT current_day FROM user WHERE netid=?''',(netid1,)) 
            checkedoutday_tup = cursor.fetchone() #retrieve the first row
            checkedoutday = checkedoutday_tup[0]
            if checkedoutday == 0:
                # Assign to this user
                #cursor.execute('''UPDATE day_status SET status = ? WHERE ymd = ? ''', (netid1, yyyymmdd,))
		#cursor.execute("UPDATE day_status SET status=%s, netid=%s WHERE ymd=%s" % (netid1, netid1, yyyymmdd))
                cursor.execute("UPDATE day_status SET status=?, netid=? WHERE ymd=? ", (netid1, netid1, yyyymmdd))
                #cursor.execute('''UPDATE day_status SET netid = ? WHERE ymd = ? ''', (netid1, yyyymmdd,))
                cursor.execute('''UPDATE user SET current_day = ? WHERE netid = ? ''', (yyyymmdd, netid1,))
                cursor = dbptr.cursor()
                dbptr.commit() 
                dirf(yyyymmdd, netid1)
                print('\nCheckout complete')
            else:
                print('You need to complete day ',checkedoutday, ' before you can check out another day')

        elif status1 == 'Done':
            print('Day ',yyyymmdd,' has already been Done')
        #elif status1 == 'Incomplete':
        #    print('Day ',yyyymmdd,' is Incomplete')
        else:
            print('Day ',yyyymmdd,' is not available. It is being processed by', status1)

    except Error as e:
        print(e)

    finally:
        if dbptr:
            dbptr.close()


# CHECK A DAY BACK IN, ONCE COMPLETE
# user:current_day -> NULL
# day_status: status -> "Done"
def checkin_day(db_file, netid1, intent):
    dbptr = None
    try:
        dbptr = sqlite3.connect(db_file)
        cursor = dbptr.cursor()
    
        # Check which day the user has checked out
        cursor.execute('''SELECT current_day FROM user 
            WHERE netid=?''',(netid1,)) 
        current_day_tup = cursor.fetchone() #retrieve the first row
        current_day = str(current_day_tup[0])
        if current_day == '0':
            print('You currently have no day checked out. So nothing to check in or release')
            return

        num_sfiles = count_sfiles(current_day)
        num_wavfiles = count_wavfiles(current_day)
        #print('For ',current_day,' there are ',num_wavfiles, ' WAV files and ',num_sfiles,' S-files')
        status1 = 'Available'
        num_wavs_not_processed = dirf(current_day, netid1)
        if num_sfiles >= num_wavfiles:
            status1 = 'Done'
        elif num_sfiles > 0:
            status1 = 'Incomplete'
        print("Day=%s, new status=%s, #S-files=%d, %d of %d WAV files remaining" % (current_day, status1, num_sfiles, num_wavs_not_processed, num_wavfiles))

        # Set user current_day to none and day_status status to finished (or Available)
        if (intent == "checkin" and status1 != "Done"):
            print("You cannot check-in an incomplete day. You must finish it, or release it (surrender)")
            return
        cursor.execute('''UPDATE day_status SET status = ?  
            WHERE ymd = ? ''', (status1, current_day, ))
        cursor.execute('''UPDATE user SET current_day = 0 
            WHERE netid = ? ''', (netid1,))

        if status1 == 'Done':
            print("")
            print("Check the following table against the event logbook")
            generate_eventlogbook(db_file, netid1)
            print("")
            print("Are all the events in the event logbook already in the table above?")
            choice = input("(y/n)")
            if choice[0]!='y':
                print("Make a note of any missing events between 20:00 the previous day, and 19:59 on %s" % (current_day))
                print("and then use the fastclassify program to classify them")
                return

            
            cursor.execute('''SELECT days_finished FROM user 
                WHERE netid=?''',(netid1,)) 
            days_finished_tup = cursor.fetchone() #retrieve the first row
            days_finished = days_finished_tup[0]
            days_finished += 1
            cursor.execute('''UPDATE user SET days_finished = ?  
                WHERE netid = ? ''', (days_finished, netid1,))

            print('\nCheckin complete')
        else:
            print('\nRelease complete')

        dbptr.commit() 

    except Error as e:
        print(e)

    finally:
        if dbptr:
            dbptr.close()

def show_tables(db_file):
    dbptr = None
    try:
        dbptr = sqlite3.connect(db_file)
        cursor = dbptr.cursor()

        cursor.execute('''SELECT * FROM day_status''')
        print("Table DAY_STATUS")
        rows = cursor.fetchall()
        for row in rows:
            print(row)

        cursor.execute('''SELECT * FROM user''')
        print("Table USER")
        rows = cursor.fetchall()
        for row in rows:
            print(row)

    except Error as e:
        print(e)
 
    finally:
        if dbptr:
            dbptr.close()

def generate_eventlogbook(db_file, netid1):
    dbptr = None
    try:
        dbptr = sqlite3.connect(db_file)
        cursor = dbptr.cursor()
    
        # Check which day the user has checked out
        cursor.execute('''SELECT current_day FROM user 
            WHERE netid=?''',(netid1,)) 
        current_day_tup = cursor.fetchone() #retrieve the first row
        current_day = str(current_day_tup[0])
        if current_day == '0':
            print('You currently have no day checked out.')
            return

        yyyymmdd = current_day
        globstr = '../../REA/ASNE_/' + yyyymmdd[0:4] + '/' + yyyymmdd[4:6] + '/' + yyyymmdd[6:8]
        #print(globstr)
        last_subclass = "dummy"
        subclass = "dummy"
        sfiles = sorted(glob.glob(globstr + '*L.S*'))
        days_per_month = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
        for sfile in sfiles:
            #print(sfile)
            sbase = os.path.basename(sfile)
            dd = int(sbase[0:2])
            hh = int(sbase[3:5])
            mi = int(sbase[5:7])
            localhh = hh - 4
            localdd = dd
            localmm = int(yyyymmdd[4:6])
            localyyyy = int(yyyymmdd[0:4])
            if localhh < 0:
                localdd -= 1
                localhh += 24
                if localdd < 1:
                    localmm -= 1
                    if localmm < 1:
                        localyyyy -= 1
                        localmm = 12
                    localdd = days_per_month[localmm-1]
            events_per_line = 0 
            with open(sfile) as fp:
                for cnt, line in enumerate(fp):
                    before_after = line.split("VOLC MAIN")
                    num_matches = len(before_after) - 1
                    #print(num_matches)
                    if num_matches > 0:
                        subclass = before_after[num_matches][0:9].strip() 
                        #print(subclass)
                        
                        if last_subclass == subclass:
                            events_per_line += 1
                            print(" %02d:%02d" % (localhh, mi), end = "")
                        else:
                            if last_subclass!="dummy":
                                while events_per_line < 6:
                                    #print("      ", end = "")
                                    events_per_line += 1
                                print("| %s" % (last_subclass))
                            events_per_line = 1
                            last_subclass = subclass
                            print("\n%02d %02d %02d | %02d:%02d " % (localdd, localmm, localyyyy%100, localhh, mi), end = "")

        print("| %s" % (subclass))
    except Error as e:
        print(e)

    finally:
        if dbptr:
            dbptr.close()


def fix_user_table(db_file, netid1, netid2, currentday, daysfinished ):
    if netid1 == 'thompsong':
        print('Set USER ',netid2,' currentday = ',currentday, ' daysfinished = ',daysfinished,' (y/n) ?')
        yn = input();
        if yn!='y':
            return
        print('*%s*' % netid2)
        if netid2 is None or not netid2:
            print('Enter correct netid for this day')
            netid2 = input()
        """ create a database connection to a SQLite database """
        dbptr = None
        return
        try:
            dbptr = sqlite3.connect(db_file)
            cursor = dbptr.cursor()
            cursor.execute('''UPDATE user SET current_day = ? WHERE netid = ? ''', (currentday, netid2,))
            cursor.execute('''UPDATE user SET days_finished = ? WHERE netid = ? ''', (daysfinished, netid2,))
            dbptr.commit() 
        except Error as e:
            print(e)
 
        finally:
            if dbptr:
                dbptr.close()

def fix_daystatus_table(db_file, netid1, yyyymmdd, status1, netid2 ):
    if netid1 == 'thompsong':

        print('Set DAY_STATUS ',yyyymmdd, ' status = ',status1, ' netid = ',netid2,' (y/n) ?')
        yn = input()
        if yn!='y':
            return
        print('*%s*' % netid2)
        if netid2 is None or not netid2:
            print('Enter correct netid for this day')
            netid2 = input()
        return
        """ create a database connection to a SQLite database """
        dbptr = None
        try:
            dbptr = sqlite3.connect(db_file)
            cursor = dbptr.cursor()
            cursor.execute('''UPDATE day_status SET status = ? WHERE ymd = ? ''', (status1, yyyymmdd,))
            cursor.execute('''UPDATE day_status SET netid = ? WHERE ymd = ? ''', (netid2, yyyymmdd,))
            dbptr.commit() 

        except Error as e:
            print(e)
 
        finally:
            if dbptr:
                dbptr.close()

if __name__ == '__main__':
    print('\n*** Some diagnostic information to send to Glenn if you have problems ***')
    # set netid1 to be value of $USER
    netid1 = os.environ['USER']
    print('USER = ',netid1) # the username
    os.system('pwd') # current working directory
    seisan_top = os.environ['SEISAN_TOP']
    print('SEISAN_TOP = ',seisan_top) # Seisan directory
    dbpath = "/shares/gly4930.018f19.gt/DB/project1_track_days.db"
    print('Using database at ',dbpath)
    os.system('ls -l %s' % dbpath) # database permissions
    os.system('groups') # groups this user belongs to
    currentDT = dt.datetime.now() # computer (local) time now
    print('Time is',str(currentDT))
    print('*** ***')
    #fix_user_table(dbpath, netid1, "lmartins1", "19960717", 0 )
    #fix_daystatus_table(dbpath, netid1, "19960711", "Incomplete", "lmartins1")
    #fix_daystatus_table(dbpath, netid1, "19960712", "Incomplete", "lmartins1")
    #fix_daystatus_table(dbpath, netid1, "19960713", "Incomplete", "lmartins1")
    #fix_daystatus_table(dbpath, netid1, "19960714", "Incomplete", "lmartins1")
    #fix_daystatus_table(dbpath, netid1, "19960715", "Incomplete", "lmartins1")
    #fix_daystatus_table(dbpath, netid1, "19960716", "Incomplete", "lmartins1")
    #fix_daystatus_table(dbpath, netid1, "19960717", "lmartins1", "lmartins1")
    print('\nMENU\n')
    loop_condition = True
    while loop_condition == True:
        print("\nPlease enter a number for what you want to do.\n")
        print("1. See available days.")
        print("2. Check-out a new day.")
        print("3. Check-in a day that is fully processed.")
        print("4. Release a day that is not fully processed.")
        print("5. Check against event logbook.")
        print("6. Quit.\n")

        try:
            main_input= int(input('Choice ? '))
        except ValueError:
            print('You need to type an integer in the range 1-6')
            main_input = 0
        if main_input == 6:
            loop_condition = False
            break
        else:
            if main_input == 1:
                list_available_days(dbpath)
            elif main_input == 2:
                checkout_day(dbpath, netid1)
            elif main_input == 3:
                checkin_day(dbpath, netid1, "checkin")
            elif main_input == 4:
                checkin_day(dbpath, netid1, "release")
            elif main_input == 5:
                generate_eventlogbook(dbpath, netid1)



