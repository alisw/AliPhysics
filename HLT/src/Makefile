#$Id$
###########################
# Makefile for HLT code.  #
#                         #
# Author: Anders Vestbo   #                    
###########################

ifndef USEPACKAGE
#USEPACKAGE = ROOT
USEPACKAGE = ALIROOT
endif

#used if USEPACKAGE=anything else
DEFSTR = -Dno_root

CXX         = g++
LD          = g++
CXXFLAGS    = -O2 -Wall -fPIC
LDFLAGS     = -O2
SOFLAGS     = -shared

ifndef TOPDIR
TOPDIR = $(HOME)/level3code
endif
ifndef LIBDIR
LIBDIR = $(LEVEL3)/lib_$(USERNAME)
endif
OBJDIR = lib

ifeq ($(USEPACKAGE),ROOT) 
INCLUDES = -I$(ROOTSYS)/include -I$(TOPDIR)/hough -I$(TOPDIR)/src -I$(TOPDIR)/comp
DEFSTR = -Duse_root
OBJDIR = lib_ROOT
endif

ifeq ($(USEPACKAGE),ALIROOT)
INCLUDES = -I$(ROOTSYS)/include -I$(TOPDIR)/hough -I$(TOPDIR)/src -I$(TOPDIR)/comp -I$(ALICE_ROOT)/include/ -I$(ALICE_ROOT)/TPC -I$(ALICE_ROOT)/CONTAINERS -I$(ALICE_ROOT)/STEER
DEFSTR = -Duse_aliroot
endif

#Save the particle id's
#DEFSTR += -Ddo_mc

#Use logging classes
ifndef NOLOGGING
DEFSTR += -Duse_logging
endif

ifdef MLUCDIR
INCLUDES += -I$(MLUCDIR)/include
else
INCLUDES += -I/prog/alice/level3/kip/MLUC/include
endif

SRCS	= AliL3ConfMapper.cxx AliL3ConfMapPoint.cxx AliLevel3.cxx\
          AliL3ConfMapTrack.cxx AliL3ConfMapFit.cxx AliL3Track.cxx\
          AliL3Transform.cxx AliL3TrackMerger.cxx AliL3VertexArray.cxx \
          AliL3Vertex.cxx AliL3VertexFinder.cxx \
          AliL3Merger.cxx AliL3GlobalMerger.cxx AliL3ClustFinderNew.cxx \
          AliL3TrackArray.cxx AliL3InterMerger.cxx AliL3Logger.cxx \
          AliL3MemHandler.cxx  AliL3Benchmark.cxx AliL3Display.cxx 

DICT = AliL3Cint.cxx
DICTH = AliL3Cint.h
DICTO = $(OBJDIR)/AliL3Cint.o

ifeq ($(USEPACKAGE),ROOT)
HDRS = $(SRCS:.cxx=.h) AliL3LinkDef.h
OBJS = $(patsubst %.cxx,$(OBJDIR)/%.o,$(SRCS)) $(DICTO)
endif

ifeq ($(USEPACKAGE),ALIROOT)
HDRS = $(SRCS:.cxx=.h) AliL3FileHandler.h AliL3Evaluate.h AliL3LinkDef.h
OBJS = $(patsubst %.cxx,$(OBJDIR)/%.o,$(SRCS)) $(OBJDIR)/AliL3Evaluate.o $(OBJDIR)/AliL3FileHandler.o $(DICTO)
endif

LEVELSO = $(LIBDIR)/libAliL3.so

all : $(OBJDIR) $(LIBDIR) $(LEVELSO) 

$(LEVELSO) : $(OBJS)
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ -o $@

$(DICT) : $(HDRS)
	@echo "Generating dictionary..."
	@rootcint -f $(DICT) -c $(INCLUDES) $(HDRS)

$(OBJDIR)/%.o : %.cxx 
	$(CXX) $(CXXFLAGS) $(DEFSTR) -c $(INCLUDES) -o $@ $<

$(OBJDIR) : 
	test ! -e $(OBJDIR) && mkdir -p $(OBJDIR)

$(LIBDIR) :
	test ! -e $(LIBDIR) && mkdir -p $(LIBDIR)

clean :
	rm -f $(OBJDIR)/*.o
	rm -f $(LIBDIR)/libAliL3.so
	rm -f $(DICT) $(DICTH) 

so :
	rm -f $(LIBDIR)/libAliL3.so


























