#$Id$
###############################
# Makefile for HLT programs.  #
#                             #
# Author: Constantin Loizides #                    
###############################

ifndef USEPACKAGE
#USEPACKAGE = ROOT
USEPACKAGE = ALIROOT
endif

#used if USEPACKAGE=anything else
DEFSTR = -Dno_root

CXX      = g++
LD       = g++
CXXFLAGS = -O2 -Wall
LDFLAGS  = -O2

ifndef TOPDIR
TOPDIR = $(HOME)/level3code
endif
ifndef LIBDIR
LIBDIR = $(LEVEL3)/lib_$(USERNAME)
endif

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
ALICELIBS = -L$(ALICE_ROOT)/lib -lCONTAINERS -lTPC -lSTEER
#ALICELIBS := $(shell aliroot-config --libs)
INCLUDES = -I$(ROOTSYS)/include -I$(ALICE_ROOT)/include/ -I$(ALICE_ROOT)/TPC -I$(ALICE_ROOT)/CONTAINERS -I$(ALICE_ROOT)/STEER
RLIBS = $(ROOTLIBS) -lMinuit -lEG
endif

#Save the particle id's
#DEFSTR += -Ddo_mc

#Use logging classes
#DEFSTR += -Duse_logging

ifdef MLUCDIR
INCLUDES += -I$(MLUCDIR)/include
MLUCLIBS = -L$(MLUCDIR)/lib/linux-i386 -lMLUC 
else
INCLUDES += -I/prog/alice/level3/kip/MLUC/include
MLUCLIBS += -L/prog/alice/level3/kip/MLUC/lib/linux-i386 -lMLUC
endif

INCLUDES += -I$(TOPDIR)/hough -I$(TOPDIR)/src -I$(TOPDIR)/comp -I$(TOPDIR)/misc
LIBS = -L$(LIBDIR) -lAliL3 -lAliL3Hough -lAliL3Comp -lAliL3Misc $(MLUCLIBS) $(ALICELIBS) $(RLIBS)

CXXFLAGS = -O2 -Wall -ggdb $(INCLUDES)
LDFLAGS  = -O2 $(LIBS)


PRGS	= runit



all : $(PRGS)

runit : runit.o $(OBJS)
	$(LD) runit.o $(LDFLAGS) $(OBJS) -o runit


all : $(OBJDIR) $(PRGS)

%.o : %.cxx 
	$(CXX) $(CXXFLAGS) $(DEFSTR) -c $(INCLUDES) -o $@ $<


clean :
	rm -f *.o $(PRGS)































