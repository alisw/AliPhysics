noinst_LTLIBRARIES = libLHAPDFCore.la
EXTRA_DIST = parmsetup.inc binreloc.h \
             common.inc commonlhacontrol.inc commonlhaglsta.inc  \
             commonlhapdfc.inc commonlhapdf.inc commonlhasets.inc

#MODULECODE   = modules.f
MODULECODE   = commoninit.f
LHAPDFCODE   = version.cc LHpdflib.F description.f evolution.f \
               parameter.F QCDparams.f alphas.f inputPDF.F
CUSTOMCODE   = wrapevolve.F Sqcdnum.f Szeus.f uncertainties.f 
PHOTONCODE   = wrapsasg.f wrapgrvg.f wrapdog.f wrapdgg.f wraplacg.f wrapgsg.f \
               wrapgsg96.f wrapacfgpg.f wrapwhitg.f
PIONCODE     = wrapowpi.f wrapsmrspi.f wrapgrvpi.f wrapabfkwpi.f 
MRSTCODE     = wrapmrst.f  
MRST06CODE   = wrapmrst2006.f  
MRST98CODE   = wrapmrst98.f  
MRSTQEDCODE  = wrapmrstqed.f  
MRSTSCODE    = Smrst.f  
CTEQCODE     = wrapcteq5.f wrapcteq65.f wrapcteq6.f wrapcteq6lg.f wrapct12.f
NNPDFCODE    = wrapNNPDF.f wrapNNPDFgrid.f wrapNNPDF20grid.f wrapNNPDF20qedgrid.f
MSTWCODE     = wrapmstw.f
GJRCODE      = wrapgjr.f
GRVCODE      = wrapgrv.f
H1CODE       = wraph1.f
ZEUSCODE     = wrapzeus.f 
HERACODE     = wraphera.f wrapheragrid.f 
ALEKHINCODE  = wrapa02m.f wrapabkm09.f wrapabm11.f
HKNCODE      = wraphkn.f
USERCODE     = wrapUSER.f wrapusergrid.f
MRSTLITE     = wrapmrst-lite.f
MRST06LITE   = wrapmrst2006-lite.f
MRSTQEDLITE  = wrapmrstqed-lite.f
MRSTSLITE    = Smrst-lite.f
CTEQLITE     = wrapcteq5.f wrapcteq65-lite.f wrapcteq6-lite.f wrapcteq6lg.f wrapct12.f
NNPDFLITE    = wrapNNPDF.f wrapNNPDFgrid-lite.f wrapNNPDF20grid-lite.f wrapNNPDF20qedgrid-lite.f
MSTWLITE     = wrapmstw-lite.f
GJRLITE      = wrapgjr-lite.f
H1LITE       = wraph1-lite.f
HERALITE     = wraphera.f wrapheragrid-lite.f 
ALEKHINLITE  = wrapa02m-lite.f wrapabkm09-lite.f wrapabm11-lite.f
HKNLITE      = wraphkn-lite.f
QCDNUMCODE   = QCDNUM.F wrapQCDNUM.F wrapQCDNUM3.f wrapQCDNUM4.f
EVLCTEQCODE  = wrapEVLCTEQ.f EVLCTEQ.f
EVLNNPDFCODE = EVLNNPDF.f wrapXNN.f
LHAGLUECODE  = lhaglue.f
EKSCODE      = eks98.f eksarp.f eps08.f eps09.f
PATHFINDCODE = getdatapath.cc binreloc.c

libLHAPDFCore_la_SOURCES = $(MODULECODE) $(LHAPDFCODE) \
                            $(EVLCTEQCODE) $(CUSTOMCODE) \
			   $(EKSCODE) $(PATHFINDCODE)
if ENABLE_LHAGLUE
libLHAPDFCore_la_SOURCES += $(LHAGLUECODE)
endif

if WANT_ALL
AM_CPPFLAGS += -DALL
endif

if WANT_QCDNUM
libLHAPDFCore_la_SOURCES += $(QCDNUMCODE)
AM_CPPFLAGS += -DQCDNUM
endif

if WANT_PIONS
libLHAPDFCore_la_SOURCES += $(PIONCODE)
AM_CPPFLAGS += -DPIONS
endif

if WANT_PHOTONS
libLHAPDFCore_la_SOURCES += $(PHOTONCODE)
AM_CPPFLAGS += -DPHOTONS
endif

if WANT_USER
libLHAPDFCore_la_SOURCES += $(USERCODE)
AM_CPPFLAGS += -DUSER
endif

if ENABLE_LOWMEM
libLHAPDFCore_la_SOURCES += $(MRSTSLITE)
AM_CPPFLAGS += -DLOW_MEMORY
if WANT_MRST
libLHAPDFCore_la_SOURCES += $(MRSTLITE)
AM_CPPFLAGS += -DMRST
endif
if WANT_MRST06
libLHAPDFCore_la_SOURCES += $(MRST06LITE)
AM_CPPFLAGS += -DMRST06
endif
if WANT_MRST98
libLHAPDFCore_la_SOURCES += $(MRST98CODE)
AM_CPPFLAGS += -DMRST98
endif
if WANT_MRSTQED
libLHAPDFCore_la_SOURCES += $(MRSTQEDLITE)
AM_CPPFLAGS += -DMRSTQED
endif
if WANT_CTEQ
libLHAPDFCore_la_SOURCES += $(CTEQLITE)
AM_CPPFLAGS += -DCTEQ
endif
if WANT_GRV
libLHAPDFCore_la_SOURCES += $(GRVCODE)
AM_CPPFLAGS += -DGRV
endif
if WANT_NNPDF
libLHAPDFCore_la_SOURCES += $(NNPDFLITE) $(EVLNNPDFCODE)
AM_CPPFLAGS += -DNNPDF
endif
if WANT_MSTW
libLHAPDFCore_la_SOURCES += $(MSTWLITE)
AM_CPPFLAGS += -DMSTW
endif
if WANT_GJR
libLHAPDFCore_la_SOURCES += $(GJRLITE)
AM_CPPFLAGS += -DGJR
endif
if WANT_H1
libLHAPDFCore_la_SOURCES += $(H1LITE)
AM_CPPFLAGS += -DH1
endif
if WANT_ZEUS
libLHAPDFCore_la_SOURCES += $(ZEUSCODE)
AM_CPPFLAGS += -DZEUS
endif
if WANT_HERA
libLHAPDFCore_la_SOURCES += $(HERALITE)
AM_CPPFLAGS += -DHERA
endif
if WANT_ALEKHIN
libLHAPDFCore_la_SOURCES += $(ALEKHINLITE)
AM_CPPFLAGS += -DALEKHIN
endif
if WANT_HKN
libLHAPDFCore_la_SOURCES += $(HKNLITE)
AM_CPPFLAGS += -DHKN
endif
else
libLHAPDFCore_la_SOURCES += $(MRSTSCODE)
AM_CPPFLAGS += -DFULL_MEMORY
if WANT_MRST
libLHAPDFCore_la_SOURCES += $(MRSTCODE)
AM_CPPFLAGS += -DMRST
endif
if WANT_MRST06
libLHAPDFCore_la_SOURCES += $(MRST06CODE)
AM_CPPFLAGS += -DMRST06
endif
if WANT_MRST98
libLHAPDFCore_la_SOURCES += $(MRST98CODE)
AM_CPPFLAGS += -DMRST98
endif
if WANT_MRSTQED
libLHAPDFCore_la_SOURCES += $(MRSTQEDCODE)
AM_CPPFLAGS += -DMRSTQED
endif
if WANT_CTEQ
libLHAPDFCore_la_SOURCES += $(CTEQCODE)
AM_CPPFLAGS += -DCTEQ
endif
if WANT_GRV
libLHAPDFCore_la_SOURCES += $(GRVCODE)
AM_CPPFLAGS += -DGRV
endif
if WANT_NNPDF
libLHAPDFCore_la_SOURCES += $(NNPDFCODE) $(EVLNNPDFCODE)
AM_CPPFLAGS += -DNNPDF
endif
if WANT_MSTW
libLHAPDFCore_la_SOURCES += $(MSTWCODE)
AM_CPPFLAGS += -DMSTW
endif
if WANT_GJR
libLHAPDFCore_la_SOURCES += $(GJRCODE)
AM_CPPFLAGS += -DGJR
endif
if WANT_H1
libLHAPDFCore_la_SOURCES += $(H1CODE)
AM_CPPFLAGS += -DH1
endif
if WANT_ZEUS
libLHAPDFCore_la_SOURCES += $(ZEUSCODE)
AM_CPPFLAGS += -DZEUS
endif
if WANT_HERA
libLHAPDFCore_la_SOURCES += $(HERACODE)
AM_CPPFLAGS += -DHERA
endif
if WANT_ALEKHIN
libLHAPDFCore_la_SOURCES += $(ALEKHINCODE)
AM_CPPFLAGS += -DALEKHIN
endif
if WANT_BOTJE
AM_CPPFLAGS += -DBOTJE
endif
if WANT_FERMI
AM_CPPFLAGS += -DFERMI
endif
if WANT_HKN
libLHAPDFCore_la_SOURCES += $(HKNCODE)
AM_CPPFLAGS += -DHKN
endif
endif

if USING_WIN32
AM_CPPFLAGS += -DWIN32
else
AM_CPPFLAGS += $(BINRELOC_CFLAGS) -DENABLE_BINRELOC
endif
AM_FFLAGS += -I. -I$(srcdir)
