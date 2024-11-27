ORIGIN = PWB
ORIGIN_VER = 2.1.49
PROJ = testf
PROJFILE = testf.mak
DEBUG = 0

FOR  = fl
FORFLAGS_R  = /Ot /Ol
FORFLAGS_D  = /Od /Zi /4Yb
CC  = cl
CFLAGS_G  = /W2 /BATCH
CFLAGS_D  = /f /Zi /Od
CFLAGS_R  = /f- /Ot /Oi /Ol /Oe /Og /Gs
CXX  = cl
CXXFLAGS_G  = /W2 /BATCH
CXXFLAGS_D  = /f /Zi /Od
CXXFLAGS_R  = /f- /Ot /Oi /Ol /Oe /Og /Gs
ASM  = ml
H2INC  = h2inc
AFLAGS_G  = /Cx /W2 /WX
AFLAGS_D  = /Zi
AFLAGS_R  = /nologo
MAPFILE_D  = NUL
MAPFILE_R  = NUL
LFLAGS_G  = /NOI /NOD /BATCH /ONERROR:NOEXE
LFLAGS_D  = /CO /FAR /PACKC
LFLAGS_R  = /EXE /FAR /PACKC
LINKER  = link
ILINK  = ilink
LRF  = echo > NUL
ILFLAGS  = /a /e
LLIBS_G  = suds144 oldnames llibce llibfore

FILES  = TESTF.FOR READTEST.FOR WRITETST.FOR
OBJS  = TESTF.obj READTEST.obj WRITETST.obj

all: $(PROJ).exe

.SUFFIXES:
.SUFFIXES:
.SUFFIXES: .obj .for

TESTF.obj : TESTF.FOR suds.fi suds.inc
!IF $(DEBUG)
        $(FOR) /c $(FORFLAGS_G) $(FORFLAGS_D) /FoTESTF.obj TESTF.FOR
!ELSE
        $(FOR) /c $(FORFLAGS_G) $(FORFLAGS_R) /FoTESTF.obj TESTF.FOR
!ENDIF

READTEST.obj : READTEST.FOR suds.fi suds.inc
!IF $(DEBUG)
        $(FOR) /c $(FORFLAGS_G) $(FORFLAGS_D) /FoREADTEST.obj READTEST.FOR
!ELSE
        $(FOR) /c $(FORFLAGS_G) $(FORFLAGS_R) /FoREADTEST.obj READTEST.FOR
!ENDIF

WRITETST.obj : WRITETST.FOR suds.fi suds.inc
!IF $(DEBUG)
        $(FOR) /c $(FORFLAGS_G) $(FORFLAGS_D) /FoWRITETST.obj WRITETST.FOR
!ELSE
        $(FOR) /c $(FORFLAGS_G) $(FORFLAGS_R) /FoWRITETST.obj WRITETST.FOR
!ENDIF


$(PROJ).exe : $(OBJS)
!IF $(DEBUG)
        $(LRF) @<<$(PROJ).lrf
$(RT_OBJS: = +^
) $(OBJS: = +^
)
$@
$(MAPFILE_D)
$(LIBS: = +^
) +
$(LLIBS_G: = +^
) +
$(LLIBS_D: = +^
)
$(DEF_FILE) $(LFLAGS_G) $(LFLAGS_D);
<<
!ELSE
        $(LRF) @<<$(PROJ).lrf
$(RT_OBJS: = +^
) $(OBJS: = +^
)
$@
$(MAPFILE_R)
$(LIBS: = +^
) +
$(LLIBS_G: = +^
) +
$(LLIBS_R: = +^
)
$(DEF_FILE) $(LFLAGS_G) $(LFLAGS_R);
<<
!ENDIF
        $(LINKER) @$(PROJ).lrf


.for.obj :
!IF $(DEBUG)
        $(FOR) /c $(FORFLAGS_G) $(FORFLAGS_D) /Fo$@ $<
!ELSE
        $(FOR) /c $(FORFLAGS_G) $(FORFLAGS_R) /Fo$@ $<
!ENDIF


run: $(PROJ).exe
        $(PROJ).exe $(RUNFLAGS)

debug: $(PROJ).exe
        CV $(CVFLAGS) $(PROJ).exe $(RUNFLAGS)
