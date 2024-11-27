rem

@rem  Make sure we have a filename as first argument
@if "%1"=="" goto help
@if "%1"=="?" goto help

@rem  Make sure we have at least one selected channel
@rem @if "%2"=="" goto help

rem   WV2ASC  -  PC-Quake to ASCII conversion in progress...

rem   Run conversion PC-Quake -> PCEQ
@wv2eq %1.WVM %1.EQ
@rem del pceq

rem   Run conversion PCEQ     -> ASCII
@eq2asc %1.EQ %1.ASC %2 %3 %4 %5 %6 %7 %8 %9
@del %1.EQ

rem   Done!
@rem
@goto end

:help
rem  WV2ASC v4.0  -  PC-Quake to ASCII (2-step) converter
rem  COMMAND MODE:
rem  Call : WV2ASC InFile /Options
rem  Where:
rem  InFile : input file name without extension;
rem           extension of the input file should be .WVM
rem  Options: /CChn1 ... /CChnN
rem  CChnI  : name of selected channel is ChnI (I=1,8) [none]
rem  The output file name is InFile.ASC .
:end
@rem del ascii
