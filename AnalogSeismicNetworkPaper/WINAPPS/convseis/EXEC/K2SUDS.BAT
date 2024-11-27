@ECHO OFF
rem K2SUDS.BAT
rem Requires KW2ASC program from KINEMETRICS, INC.

rem K2&MountWhitney to PC-SUDS conversion in progress ...

kw2asc %1 -H
k22suds %1 %2 %3 %4 %5 %6 %7 %8 %9
del *.0??
del *.shd
del *.hdr
