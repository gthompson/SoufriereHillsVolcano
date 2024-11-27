rem

@if "%1"=="?" goto help

rem  SDP2EQ v4.0 -  SDP to PCEQ conversion in progress...

@echo off
if not exist arriv.dat goto help
if not exist event.dat goto help
if not exist waveform.dat goto help
sdp1 arriv.dat
if errorlevel 1 goto exit1
sdp1 event.dat
if errorlevel 1 goto exit2
sdp1 waveform.dat
if errorlevel 1 goto exit3
sdp2
del arriv.* >null
del event.* >null
del waveform.* >null
if exist nul del null
@goto exit

:exit1
@echo on
rem  Error while converting file ARRIV.DAT
goto exit
:exit2
@echo on
rem  Error while converting file EVENT.DAT
goto exit
:exit3
@echo on
rem  Error while converting file WAVEFORM.DAT
goto exit

:help
@echo on
rem  SDP2EQ v4.0  -  SDP to PCEQ converter
rem
rem  INTERACTIVE MODE:
rem  Call:    SDP2EQ
rem  The input files in SDP Teledyne Geotech format are:
rem  EVENT.DAT, ARRIV.DAT  & WAVEFORM.DAT
:exit

