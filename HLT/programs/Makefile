#$Id$
###############################
# Makefile for HLT programs.  #
#                             #
# Author: Constantin Loizides #                    
###############################


ifndef ALIHLT_TOPDIR
ALIHLT_TOPDIR = $(shell pwd)/..
endif

include $(ALIHLT_TOPDIR)/Makefile.conf

DEFSTR = -Dno_root
INCLUDES = -I$(ALIHLT_TOPDIR)/src -I$(ALIHLT_TOPDIR)/hough -I$(ALIHLT_TOPDIR)/comp -I$(ALIHLT_TOPDIR)/misc -I$(ALIHLT_TOPDIR)/trigger -I$(ALIHLT_TOPDIR)/kalman

LIBS = -L$(ALIHLT_LIBDIR) -lAliHLTSrc -lAliHLTHough -lAliHLTComp -lAliHLTMisc 
STATICLIBS = $(ALIHLT_LIBDIR)/libAliHLTSrc.a  $(ALIHLT_LIBDIR)/libAliHLTHough.a  $(ALIHLT_LIBDIR)/libAliHLTComp.a  $(ALIHLT_LIBDIR)/libAliHLTMisc.a 

ifneq ($(ALIHLT_USEPACKAGE),STANDALONE)
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)
ROOTGLIBS  := $(shell root-config --glibs)
endif

ifeq ($(ALIHLT_USEPACKAGE),ROOT) 
DEFSTR = -Duse_root
INCLUDES += -I$(ROOTSYS)/include 
RLIBS = $(ROOTLIBS) -lMinuit -lEG 
LIBS += $(RLIBS)
STATICLIBS += $(RLIBS)
endif

ifeq ($(ALIHLT_USEPACKAGE),ALIROOT)
DEFSTR = -Duse_aliroot -Duse_root
ifeq ($(USENEWIO),1)
DEFSTR += -Duse_newio
ALICELIBS = -L$(ALICE_ROOT)/lib/tgt_Linux -lTPC -lITS -lCONTAINERS -lSTEER -lRAWDatabase -lRAWDatarec -lRAWDatasim -lVMC -lGeom
else
ALICELIBS = -L$(ALICE_ROOT)/lib/tgt_Linux -lTPC -lITS -lCONTAINERS -lSTEER -lMC
endif
INCLUDES += -I$(ROOTSYS)/include -I$(ALICE_ROOT)/include/ -I$(ALICE_ROOT)/TPC -I$(ALICE_ROOT)/CONTAINERS -I$(ALICE_ROOT)/STEER
RLIBS = $(ROOTLIBS) -lMinuit -lEG
LIBS += $(ALICELIBS) $(RLIBS)
STATICLIBS += $(ALICELIBS) $(RLIBS)
endif

ifeq ($(ALIHLT_USEPACKAGE),STANDALONE)
DEFSTR = -Dno_root
RLIBS =  -lpthread -lm -ldl -rdynamic
LIBS += $(RLIBS)
STATICLIBS += $(RLIBS)
endif

ifeq ($(DOMC),1)
DEFSTR += -Ddo_mc
endif

#Use logging classes
ifndef NOLOGGING
DEFSTR += -Duse_logging
ifdef ALIHLT_MLUCDIR
INCLUDES += -I$(ALIHLT_MLUCDIR)/include
MLUCLIBS = -L$(ALIHLT_MLUCDIR)/lib/ -lMLUC 
else
INCLUDES += -I/prog/alice/level3/kip/MLUC/include
MLUCLIBS += -L/prog/alice/level3/kip/MLUC/lib/linux-i386 -lMLUC
endif
LIBS += $(MLUCLIBS)
STATICLIBS += $(MLUCLIBS)
endif

DEFSTR += -D$(ARCH) $(EXTRADEF)

CXX          = g++
LD           = $(CXX)
#CXXGCC3FLAGS += -pedantic
#CXXGCC3FLAGS += -Wno-deprecated
#CXXGCC3FLAGS += -Woverloaded-virtual
CXXFLAGS       = -O2 -Wall -ggdb $(CXXGCC3FLAGS) $(EXTRACXXFLAGS) 
LDFLAGS        = -O2 $(EXTRALDFLAGS) $(LIBS)
STATICLDFLAGS  = -O2 $(PROFILEFLAGS) $(EXTRALDFLAGS) $(STATICLIBS)


PRGS	=  runit read runvhdlcf runhough runvhdlhough runtracker tpcbeamtesttracker tpcbeamtestdisplay
STPRGS	= 
UTILS   = ccfile gettransform
STUTILS = gettransform_st runit_st
UTILS_not_working  = convbin speedtest 

all: $(PRGS) $(UTILS)

programs: $(PRGS)

utils: $(UTILS)

stprograms: $(STPRGS)

stutils: $(STUTILS)

#run hough transform on l3 data.
runhough: runhough.o
	$(LD) $< $(LDFLAGS) -o $@

#run hough transform on l3 data.
runvhdlhough: runvhdlhough.o
	$(LD) $< $(LDFLAGS) -o $@

#run conventional clusterfinder on altro data
runvhdlcf: runvhdlcf.o
	$(LD) $< $(LDFLAGS) -o $@

#run conventional clusterfinder
runit: runit.o
	$(LD) $< $(LDFLAGS) -o $@

runit_st: runit.sto
	$(LD) $< $(LDFLAGS) -o $@

#run conventional tracker
runtracker: runtracker.o
	$(LD) $< $(STATICLDFLAGS) -o $@

#TPC Test cluster finder and tracker
tpcbeamtesttracker: tpcbeamtesttracker.o
	$(LD) $< $(LDFLAGS) -o $@

#TPC Test display
tpcbeamtestdisplay: tpcbeamtestdisplay.o
	$(LD) $< $(LDFLAGS) -o $@

#read l3 raw data, print it and convert to altro data
read: read.o
	$(LD) $< $(LDFLAGS) -o $@

#convert big <-> little endian cosmics data file
ccfile: convcosmicsfile.o
	$(LD) $< $(LDFLAGS) -o $@

#convert to binary
convbin: ali2raw.o
	$(LD) $< $(LDFLAGS) -o $@

#get transform values
gettransform: gettransform.o
	$(LD) $< $(LDFLAGS) -o $@

gettransform_st: gettransform.sto
	$(LD) $< $(STATICLDFLAGS) -o $@


speedtest: speedtest.cxx
	g++ -O -Wall -fPIC -ffast-math $(CXXFLAGS) speedtest.cxx -o speedtest
#       g++ -O -Wall -fPIC -mcpu=ev6   $(CXXFLAGS) speedtest.cxx -o speedtest
#       cxx -O -Wall -fPIC -lm         $(CXXFLAGS) speedtest.cxx -o speedtest 

%.o:%.cxx 
	$(CXX) $(CXXFLAGS) $(DEFSTR) -c $(INCLUDES) -o $@ $<

%.sto:%.cxx 
	$(CXX) $(PROFILEFLAGS) $(CXXFLAGS) $(DEFSTR) -c $(INCLUDES) -o $@ $<

clean :
	rm -f *.o *.sto gmon.out $(PRGS) $(UTILS) $(STPRGS) $(STUTILS) $(UTILS_not_working)

