# NMAKE description file for the SUDS library

# Specify library filename and destination below
libname = SUDS145
libdir  = C:\MSVC\USERLIBS

.SUFFIXES: .for .c .asm .obj .exe .lib
.c.obj:
   @CL $(coptions) $<

# C compiler options
coptions = /c /FPi /AL /Ox

$(libname).lib: sudsio.obj sudsinit.obj sudstime.obj
   if exist $(libname).lib del $(libname).lib
   lib @<<
$(libname).lib
y
+sudsio +sudstime +sudsinit
NUL
<<
   copy $(libname).lib $(libdir)

sudsio.obj:   sudsio.c suds.h

sudsinit.obj: sudsinit.c suds.h sudsform.h

sudstime.obj: sudstime.c suds.h
