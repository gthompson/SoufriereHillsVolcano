@ECHO OFF

REM PROCESS.BAT, Thu 04-Jan-1996 13:35, RB

REM This is the process batch file that is called by MERGE every time it
REM adds a record to the merge buffer, merges a record into the merge
REM buffer, or purges the merge buffer.
                                                 
REM The file specification is passed as the first three arguments to this
REM batch file.  These are the pathname, filename, and extension 
REM respectively.  

REM The fourth argument is a status indication passed as a string as follows:

REM  "%4" = "MERGED:n" = This file made up of n merged records.  Note 
REM                      that for a complete record, the file is 
REM                      processed at the final merge and then again 
REM                      when purged as a complete file.
REM  "%4" = "COMPLETE" = This file is being purge because it is 
REM                      complete. It contains records .GE. the
REM                      MergeCompleteCount entry in MERGE.INI.
REM  "%4" = "STALE"    = This file is being purged because it has been
REM                      in the merge buffer longer than allowed by the
REM                      StaleFile entry in MERGE.INI.
              
REM Assert that a complete and valid filespec was passed in.
IF "%1" == "" GOTO ASSERT
IF "%2" == "" GOTO ASSERT
IF "%3" == "" GOTO ASSERT                  

REM Assert that the file exists
IF EXIST %1\%2.%3 GOTO START
ECHO ERROR: File does not exist: %1\%2.%3
PAUSE > NUL
GOTO END

:ASSERT
ECHO ERROR: Incomplete or invalid filespec passed was to PROCESS.BAT.
ECHO Argument 1: %1
ECHO Argument 2: %2
ECHO Argument 3: %3
PAUSE > NUL
GOTO END

:START

ECHO.
ECHO %4 - %1\%2.%3
ECHO.               

REM - If this file was purged as a complete record...
IF "%4" == "COMPLETE" GOTO COMPLETE        
REM - If this file was purged as an incomplete record...
IF "%4" == "STALE" GOTO STALEFILE

REM -------------------------------------------------------------------
REM - Perform routine processing for each merged record here...
ECHO Process merged record: %1\%2.%3


GOTO END

REM -------------------------------------------------------------------
REM - Process complete records here...
:COMPLETE
ECHO Process complete record: %1\%2.%3        


GOTO END

REM -------------------------------------------------------------------
REM - Handle stale files here...
:STALEFILE
ECHO Handle stale file: %1\%2.%3        


GOTO END      
                                                   
:END
