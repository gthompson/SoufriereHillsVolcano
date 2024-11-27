@ECHO OFF
rem CSS2SUD1.BAT v1.1 1-component converter

if "%1"=="" goto help
if "%2"=="" goto help
if "%3"=="" goto help

css2sga %1 temp.inp /C%3
sg2suds temp.inp temp.out /fA > nul
sudsrw temp.out %2 > nul
if exist temp.* del temp.* > nul
goto end

:help
echo
echo Usage: CSS2SUD1 InpFile OutFile Component (1-comp. converter)
echo        (where Component can be Z,V,N,E)
echo

:end
