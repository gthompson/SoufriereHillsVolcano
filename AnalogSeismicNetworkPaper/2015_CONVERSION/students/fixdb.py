#!/shares/gly4930.018f19.gt/miniconda3/bin/python
#import sqlite3
#from sqlite3 import Error
#import pandas as pd 
#import datetime as dt
#import os
#import sys
#import glob
#import re
import project1
#SEISAN_TOP = '/shares/gly4930.018f19.gt/seisan'

if __name__ == '__main__':
    dbpath = "/shares/gly4930.018f19.gt/DB/project1_track_days.db"
    print('Using database at ',dbpath)
    print('\nMENU\n')
    loop_condition = True
    while loop_condition == True:
        print("\nPlease enter a number for what you want to do.\n")
        print("1. Create new DB. (Glenn only).")
        print("2. Show tables.")
        print("3. Fix record in USER table.")
        print("4. Fix record in DAY_STATUS table.")
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
                project1.initialize_db(dbpath, netid1)
            elif main_input == 2:
                project1.show_tables(dbpath)
            elif main_input == 3:
                print("Enter NETID of person whose record you want to change")
                netid2 = input()
                print("Enter YYYYMMDD of date that %s is working on. Or hit ENTER if they are not working on any days." % (netid2,))
                yyyymmdd = input()
                print("Enter number of days that %s has finished processing" % (netid2,))
                daysfinished = input()
                print("Okay, do you want to change record for %s to current working on day %s and days done = %d" % (netid2, yyyymmdd, int(daysfinished),) )
                yn = input('y or n?')
                if yn[0]=='y':
                    project1.fix_user_table(dbpath, "thompsong", netid2, yyyymmdd, daysfinished)
            elif main_input == 4:
                print("Enter YYYYMMDD of date you wish to change")
                yyyymmdd = input()
                print("Enter the new status for day %s. Should be Available, Incomplete, Done or Checkedout" % (yyyymmdd,) )
                status = input()
                print("Okay, do you want to change record for %s to status %s ?" % (yyyymmdd, status,) )
                yn = input('y or n?')
                if yn[0]=='y':
                    project1.fix_daystatus_table(dbpath, "thompsong", yyyymmdd, status, "")
            elif main_input == 5:
                project1.generate_eventlogbook(dbpath, netid1)



