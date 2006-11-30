# Makefile for ALICE MDC program for Linux

include $(ROOTSYS)/test/Makefile.arch

LDFLAGS      += -Wl,-rpath,/ROOT/Linux/CurrentRelease/root/lib \
                -Wl,-rpath,/date/smi/linux \
                -Wl,-rpath,/date/dim/linux \
                -Wl,-rpath,/adcRoot/alimdc5

ifneq ($(findstring rdm,$(MAKECMDGOALS)),)
CXXFLAGS += -DUSE_RDM
endif

ifneq ($(findstring debug,$(MAKECMDGOALS)),)
CXXFLAGS += -DUSE_DEBUG -g
endif

ifneq ($(findstring smi,$(MAKECMDGOALS)),)
#SMIDIR    = ./smi
#DIMDIR    = ./dim
SMIDIR    = /date/smi
DIMDIR    = /date/dim
CXXFLAGS += -DUSE_SMI -I$(SMIDIR)/smixx -I$(DIMDIR)/dim
LIBS     += -L$(SMIDIR)/linux -lsmiui -lsmi -L$(DIMDIR)/linux -ldim
endif

ifneq ($(findstring evb,$(MAKECMDGOALS)),)
ifneq ($(findstring rdm,$(MAKECMDGOALS)),)
EBDIR     = ./eb
BMDIR     = ./eb
DBMDIR    = ./eb
DBDIR     = ./eb
FIFODIR   = ./eb
INFODIR   = ./eb
else
EBDIR     = /date/eventBuilder
BMDIR     = /date/banksManager
DBMDIR    = /date/bufferManager
DBDIR     = /date/db
FIFODIR   = /date/simpleFifo
INFODIR   = /date/infoLogger
endif
CXXFLAGS += -DUSE_EB -I$(EBDIR)
LIBS     += -L$(EBDIR)/Linux -L$(BMDIR)/Linux -L$(DBMDIR)/Linux \
            -L$(DBDIR)/Linux -L$(FIFODIR)/Linux -L$(INFODIR)/Linux \
            -lDateEb -lBanksManager -lDateBufferManager \
            -lDb -lFifo -lInfo
endif

ifneq ($(findstring hlt,$(MAKECMDGOALS)),)
include ../HLT/hlt.conf
CXXFLAGS += -DUSE_HLT ${HLTDEFS} -I${ALICE_ROOT}/include -I../HLT/src -I../HLT/hough -I../HLT/comp
LIBS     += -lEG -lVMC -lGeom -lThread -L../lib/tgt_linuxDEBUG -lSTEER -lRAWData -lTPCbase -lTPCsim -lTPCrec -lCONTAINERS -lAliHLTSrc -lAliHLTMisc -lAliHLTHough -lAliHLTComp
endif

##### DATE #####

ifdef DATE_ROOT
DATEFLAGS = -DALI_DATE -D${DATE_SYS} -DDATE_SYS=${DATE_SYS} -Dlong32=${DATE_LONG32} -Dlong64='${DATE_LONG64}' -DdatePointer=${DATE_POINTER} -I${DATE_COMMON_DEFS}
CXXFLAGS += $(DATEFLAGS)
endif

##### MACROS and TARGETS #####

include Make-macros


##### DEPENDENCIES #####

#include Make-depend
