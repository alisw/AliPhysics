#$Id$
###########################
# Makefile for HLT code.  #
#                         #
# Author: Anders Vestbo   #                    
###########################

CXX         = g++
LD          = g++
CXXFLAGS    = -O2 -Wall -fPIC
LDFLAGS     = -O2
SOFLAGS     = -shared

TOPDIR = $(HOME)/level3code
OBJDIR = lib
LIBDIR = $(LEVEL3)/lib_$(USERNAME)

SRCS	= AliL3HoughTransformer.cxx AliL3Hough.cxx AliL3HoughTrack.cxx\
	  AliL3HoughMaxFinder.cxx AliL3HoughEval.cxx AliL3HoughMerger.cxx \
	  AliL3Histogram.cxx AliL3Histogram1D.cxx \
	  AliL3HoughIntMerger.cxx AliL3HoughGlobalMerger.cxx 

DICT = AliL3HoughCint.cxx
DICTH = AliL3HoughCint.h
DICTO = $(OBJDIR)/AliL3HoughCint.o

HDRS = $(SRCS:.cxx=.h) AliL3HoughLinkDef.h

OBJS = $(patsubst %.cxx,$(OBJDIR)/%.o,$(SRCS)) $(DICTO)

LEVELSO = $(LIBDIR)/libAliL3Hough.so

all : $(OBJDIR) $(LIBDIR) $(LEVELSO)

$(LEVELSO) : $(OBJS)
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ -o $@

$(DICT) : $(HDRS)
	@echo "Generating dictionary..."
	@rootcint -f $(DICT) -c -I$(ROOTSYS) -I$(TOPDIR)/src $(HDRS)


$(OBJDIR)/%.o : %.cxx 
	$(CXX) $(CXXFLAGS) -c -I$(ROOTSYS)/include -I$(TOPDIR)/src -I$(TOPDIR)/hough -o $@ $<

$(OBJDIR) : 
	test ! -e $(OBJDIR) && mkdir -p $(OBJDIR)
$(LIBDIR) :
	test ! -e $(LIBDIR) && mkdir -p $(LIBDIR)

clean :
	rm -f $(OBJDIR)/*.o
	rm -f $(LIBDIR)/libAliL3Hough.so
	rm $(DICT) $(DICTH) 

so: 
	rm -f $(LIBDIR)/libAliL3Hough.so