#$Id$
###############################
# Makefile for HLT programs.  #
#                             #
# Author: Constantin Loizides #                    
###############################


ifndef TOPDIR
TOPDIR = `pwd`/..
endif

include $(TOPDIR)/Makefile.conf

DEFSTR = -Dno_root
INCLUDES = -I$(TOPDIR)/hough -I$(TOPDIR)/src -I$(TOPDIR)/comp -I$(TOPDIR)/misc
LIBS = -L$(LIBDIR) -lAliL3 -lAliL3Hough -lAliL3Comp -lAliL3Misc 

ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)
ROOTGLIBS  := $(shell root-config --glibs)

ifeq ($(USEPACKAGE),ROOT) 
DEFSTR = -Duse_root
INCLUDES += -I$(ROOTSYS)/include 
RLIBS = $(ROOTLIBS) -lMinuit -lEG 
LIBS += $(RLIBS)
endif

ifeq ($(USEPACKAGE),ALIROOT)
DEFSTR = -Duse_aliroot -Duse_root
ALICELIBS = -L$(ALICE_ROOT)/lib/tgt_Linux -lCONTAINERS -lTPC -lSTEER
#ALICELIBS := $(shell aliroot-config --libs)
INCLUDES += -I$(ROOTSYS)/include -I$(ALICE_ROOT)/include/ -I$(ALICE_ROOT)/TPC -I$(ALICE_ROOT)/CONTAINERS -I$(ALICE_ROOT)/STEER
RLIBS = $(ROOTLIBS) -lMinuit -lEG
LIBS += $(ALICELIBS)
endif

ifdef DO_MC
DEFSTR += -Ddo_mc
endif

#Use logging classes
ifndef NOLOGGING
DEFSTR += -Duse_logging

ifdef MLUCDIR
INCLUDES += -I$(MLUCDIR)/include
MLUCLIBS = -L$(MLUCDIR)/lib/linux-i386 -lMLUC 
else
INCLUDES += -I/prog/alice/level3/kip/MLUC/include
MLUCLIBS += -L/prog/alice/level3/kip/MLUC/lib/linux-i386 -lMLUC
endif
LIBS += $(MLUCLIBS)
endif

DEFSTR += $(EXTRADEF)

CXXFLAGS = -O2 -Wall -ggdb $(EXTRACXXFLAGS) $(INCLUDES)
LDFLAGS  = -O2 $(EXTRALDFLAGS) $(LIBS)

PRGS	= runit read  runvhdlcf runhough
UTILS   = convbin speedtest gettransform


all : $(PRGS)

#run hough transform on l3 data.
runhough : runhough.o
	$(LD) $< $(LDFLAGS) -o $@

#run conventional clusterfinder on altro data
runvhdlcf : runvhdlcf.o
	$(LD) $< $(LDFLAGS) -o $@

#run conventional clusterfinder
runit : runit.o
	$(LD) $< $(LDFLAGS) -o $@

convbin : ali2raw.o
	$(LD) $< $(LDFLAGS) -o $@

#read l3 raw data, print it and convert to altro data
read : read.o
	$(LD) $< $(LDFLAGS) -o $@

gettransform : gettransform.o
	$(LD) $< $(LDFLAGS) -o $@

speedtest : speedtest.cxx
	g++ -O -Wall -fPIC -ffast-math $(CXXFLAGS) speedtest.cxx -o speedtest
	# g++ -O -Wall -fPIC -mcpu=ev6 $(CXXFLAGS) speedtest.cxx -o speedtest
	# cxx -O -Wall -fPIC  -lm  $(CXXFLAGS) speedtest.cxx -o speedtest 

%.o : %.cxx 
	$(CXX) $(CXXFLAGS) $(DEFSTR) -c $(INCLUDES) -o $@ $<

clean :
	rm -f *.o $(PRGS) $(UTILS)


