#$Id$
###########################
# Makefile for HLT code.  #
#                         #
# Author: Anders Vestbo   #                    
###########################

#USEPACKAGE = ROOT
USEPACKAGE = ALIROOT

DEFSTR = -Dno_root

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
DEFSTR = -Duse_root
endif

ifeq ($(USEPACKAGE),ALIROOT)
INCLUDES = -I$(ROOTSYS)/include -I$(TOPDIR)/hough -I$(TOPDIR)/src -I$(TOPDIR)/comp -I$(ALICE_ROOT)/include/ -I$(ALICE_ROOT)/TPC -I$(ALICE_ROOT)/CONTAINERS -I$(ALICE_ROOT)/STEER
DEFSTR = -Duse_aliroot -Duse_root
endif

#Save the particle id's
#DEFSTR += -Ddo_mc

SRCS	= AliL3HoughTransformer.cxx AliL3Hough.cxx AliL3HoughTrack.cxx\
	  AliL3HoughMaxFinder.cxx AliL3HoughEval.cxx AliL3HoughMerger.cxx \
	  AliL3Histogram.cxx AliL3Histogram1D.cxx \
	  AliL3HoughIntMerger.cxx AliL3HoughGlobalMerger.cxx GetGoodParticles.cxx

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
	@rootcint -f $(DICT) -c $(INCLUDES) $(HDRS)


$(OBJDIR)/%.o : %.cxx 
	$(CXX) $(CXXFLAGS) $(DEFSTR) -c $(INCLUDES) -o $@ $<

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