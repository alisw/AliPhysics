#
# Makefile for ROOT version of HLT code.
#
# Author: Anders Vestbo                       
#

CXX         = g++
LD          = g++
CXXFLAGS    = -O2 -Wall -fPIC
LDFLAGS     = -O2
SOFLAGS     = -shared

TOPDIR = $(HOME)/level3code
OBJDIR = tgt_Linux

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

HDRS = $(SRCS:.cxx=.h) AliL3LinkDef.h

OBJS = $(patsubst %.cxx,$(OBJDIR)/%.o,$(SRCS)) $(DICTO)

LEVELSO = $(LEVEL3)/lib_$(USERNAME)/libAliL3.so

all : $(OBJDIR) $(LEVELSO)

$(LEVELSO) : $(OBJS)
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ -o $@

$(DICT) : $(HDRS)
	@echo "Generating dictionary..."
	@rootcint -f $(DICT) -c $(HDRS)

$(OBJDIR)/%.o : %.cxx 
	$(CXX) $(CXXFLAGS) -c -I$(ROOTSYS)/include -I$(TOPDIR)/hough -I$(TOPDIR)/src -o $@ $<

$(OBJDIR) : 
	test ! -e $(OBJDIR) && mkdir -p $(OBJDIR)

clean :
	rm -f $(OBJDIR)/*.o
	rm -f $(LEVEL3)/lib_$(USERNAME)/libAliL3.so
	rm $(DICT) $(DICTH) 