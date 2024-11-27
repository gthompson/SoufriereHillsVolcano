@ECHO OFF

REM SOURCE1.BAT, Thu 04-Jan-1996 13:35, RB

REM This is a preprocess batch file that is called by MERGE every time it
REM moves a file from a data source.  The file should be processed into
REM a time-corrected, demuxed PC-SUDS data file and then moved
REM to the working directory (see MERGE.INI).
                                                 
REM The file specification is passed as the first three arguments to this
REM batch file.  These are the; pathname, filename, and extension 
REM respectively.
              
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
ECHO ERROR: Incomplete or invalid filespec passed was to PREPROC1.BAT.
ECHO Argument 1: %1
ECHO Argument 2: %2
ECHO Argument 3: %3
PAUSE > NUL
GOTO END

:START
REM -------------------------------------------------------------------
REM All user processing starts here.

REM Convert data to PC-SUDS format

:END
