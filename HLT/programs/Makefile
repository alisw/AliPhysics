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

ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)
ROOTGLIBS  := $(shell root-config --glibs)

ifeq ($(USEPACKAGE),ROOT) 
DEFSTR = -Duse_root
INCLUDES = -I$(ROOTSYS)/include 
RLIBS = $(ROOTLIBS) -lMinuit -lEG 
endif

ifeq ($(USEPACKAGE),ALIROOT)
DEFSTR = -Duse_aliroot
ALICELIBS = -L$(ALICE_ROOT)/lib/tgt_Linux -lCONTAINERS -lTPC -lSTEER
#ALICELIBS := $(shell aliroot-config --libs)
INCLUDES = -I$(ROOTSYS)/include -I$(ALICE_ROOT)/include/ -I$(ALICE_ROOT)/TPC -I$(ALICE_ROOT)/CONTAINERS -I$(ALICE_ROOT)/STEER
RLIBS = $(ROOTLIBS) -lMinuit -lEG
endif

#Save the particle id's
#DEFSTR += -Ddo_mc

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
endif

INCLUDES += -I$(TOPDIR)/hough -I$(TOPDIR)/src -I$(TOPDIR)/comp -I$(TOPDIR)/misc
LIBS = -L$(LIBDIR) -lAliL3 -lAliL3Hough -lAliL3Comp -lAliL3Misc $(MLUCLIBS) $(ALICELIBS) $(RLIBS)

CXXFLAGS = -O2 -Wall -ggdb $(INCLUDES)
LDFLAGS  = -O2 $(LIBS)


PRGS	= runit read speedtest convbin gettransform runvhdlcf


all : $(PRGS)

runit : runit.o
	$(LD) $< $(LDFLAGS) -o $@

convbin : ali2raw.o
	$(LD) $< $(LDFLAGS) -o $@

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
	rm -f *.o $(PRGS)































