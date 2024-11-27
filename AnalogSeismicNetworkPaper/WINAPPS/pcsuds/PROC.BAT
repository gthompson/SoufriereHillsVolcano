@ECHO OFF
REM This batch file is spawned by SUDSMAN every time it moves a file
REM from it's source directory.

REM %1 = The pathname of the destination file.
REM %2 = The filename of the destination file without path or extension.
REM %3 = The filename extension of the destination file.

REM Assert that a complete and valid filespec was passed in.
IF "%1" == "" GOTO ASSERT
IF "%2" == "" GOTO ASSERT
IF "%3" == "" GOTO ASSERT
IF EXIST %1\%2.%3 GOTO START

:ASSERT
ECHO Assertion Failed: Incomplete or invalid filespec passed to PROC.BAT.
ECHO                   %1\%2.%3
PAUSE > NUL
GOTO END

:START
REM -------------------------------------------------------------------
REM All user processing starts here.

ECHO PROC.BAT was passed: %1\%2.%3

REM If source and destination are the same, you must move or delete the file.
REM MOVE %1\%2.%3 C:\SUDSUTIL\SUDSMAN\DEST > NUL



REM All user processing ends here.
REM -------------------------------------------------------------------
:END
