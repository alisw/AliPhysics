#$Id$
###########################
# Makefile for HLT code.  #
#                         #
# Author: Anders Vestbo   #                    
###########################

USEPACKAGE = ROOT
#USEPACKAGE = ALIROOT

CXX         = g++
LD          = g++
CXXFLAGS    = -O2 -Wall -fPIC
LDFLAGS     = -O2
SOFLAGS     = -shared

TOPDIR = $(HOME)/level3code
OBJDIR = lib
LIBDIR = $(LEVEL3)/lib_$(USERNAME)

ifeq ($(USEPACKAGE),ROOT) 
INCLUDES = -I$(ROOTSYS)/include -I$(TOPDIR)/hough -I$(TOPDIR)/src -I$(TOPDIR)/comp
endif

ifeq ($(USEPACKAGE),ALIROOT)
INCLUDES = -I$(ROOTSYS)/include -I$(TOPDIR)/hough -I$(TOPDIR)/src -I$(TOPDIR)/comp -I$(ALICE_ROOT)/include/ -I$(ALICE_ROOT)/TPC -I$(ALICE_ROOT)/CONTAINERS -I$(ALICE_ROOT)/STEER
endif

SRCS	= AliL3ConfMapper.cxx AliL3ConfMapPoint.cxx AliLevel3.cxx\
          AliL3ConfMapTrack.cxx AliL3ConfMapFit.cxx AliL3Track.cxx\
          AliL3Transform.cxx AliL3TrackMerger.cxx AliL3VertexArray.cxx \
          AliL3Vertex.cxx AliL3VertexFinder.cxx \
          AliL3Merger.cxx AliL3GlobalMerger.cxx AliL3ClustFinderNew.cxx \
          AliL3TrackArray.cxx AliL3InterMerger.cxx AliL3Logger.cxx \
          AliL3MemHandler.cxx  AliL3Benchmark.cxx AliL3Display.cxx\

DICT = AliL3Cint.cxx
DICTH = AliL3Cint.h
DICTO = $(OBJDIR)/AliL3Cint.o

ifeq ($(USEPACKAGE),ROOT)
HDRS = $(SRCS:.cxx=.h) AliL3LinkDef.h
OBJS = $(patsubst %.cxx,$(OBJDIR)/%.o,$(SRCS)) $(DICTO)
endif

ifeq ($(USEPACKAGE),ALIROOT)
HDRS = $(SRCS:.cxx=.h) AliL3FileHandler.h AliL3LinkDef.h
OBJS = $(patsubst %.cxx,$(OBJDIR)/%.o,$(SRCS)) $(OBJDIR)/AliL3FileHandler.o $(DICTO)
endif

LEVELSO = $(LIBDIR)/libAliL3.so

all : $(OBJDIR) $(LEVELSO) 

$(LEVELSO) : $(OBJS)
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ -o $@

$(DICT) : $(HDRS)
	@echo "Generating dictionary..."
	@rootcint -f $(DICT) -c $(INCLUDES) $(HDRS)

$(OBJDIR)/%.o : %.cxx 
	$(CXX) $(CXXFLAGS) -c $(INCLUDES) -o $@ $<

$(OBJDIR) : 
	test ! -e $(OBJDIR) && mkdir -p $(OBJDIR)

clean :
	rm -f $(OBJDIR)/*.o
	rm -f $(LIBDIR)/libAliL3.so
	rm $(DICT) $(DICTH) 

so :
	rm -f $(LIBDIR)/libAliL3.so