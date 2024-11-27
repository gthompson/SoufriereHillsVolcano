@ECHO OFF
rem CSS2SUD3.BAT v1.1 3-component converter

if "%1"=="" goto help
if "%2"=="" goto help
if "%3"=="" goto help
if "%4"=="" goto help

css2sga %1 temp.z /CZ
css2sga %2 temp.n /CN
css2sga %3 temp.e /CE
sga2sga temp.z temp.n temp.e temp.inp > nul
sg2suds temp.inp temp.out /fA > nul
sudsrw temp.out %4 > nul
if exist temp.* del temp.* > nul
goto end

:help
echo
echo Usage: CSS2SUD3 InpFileZ InpFileNS InpFileEW OutFile (3-comp. converter)
echo

:end
