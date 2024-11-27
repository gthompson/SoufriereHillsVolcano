@if "%1"=="" goto help
@if "%1"=="?" goto help
sg2eq %1.bin temp.eq
eq2seis temp.eq /N%2
del temp.eq
goto end
:help
REM
REM     Calling syntax:   SG2SEIS InputFile n
REM              where:           InputFile has extension .BIN
REM                               n is a current number between 1-999
REM
:end
