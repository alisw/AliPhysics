#$Id$
###########################
# Makefile for HLT code.  #
#                         #
# Author: Anders Vestbo   #                    
###########################

ifndef TOPDIR
TOPDIR = `pwd`/..
endif

include $(TOPDIR)/Makefile.conf

DEFSTR = -Dno_root

ifeq ($(USEPACKAGE),ROOT) 
INCLUDES = -I$(ROOTSYS)/include -I$(TOPDIR)/hough -I$(TOPDIR)/src -I$(TOPDIR)/comp -I$(TOPDIR)/misc
DEFSTR = -Duse_root
OBJDIR = lib_ROOT
endif

ifeq ($(USEPACKAGE),ALIROOT)
INCLUDES = -I$(ROOTSYS)/include -I$(TOPDIR)/hough -I$(TOPDIR)/src -I$(TOPDIR)/comp -I$(TOPDIR)/misc -I$(ALICE_ROOT)/include/ -I$(ALICE_ROOT)/TPC -I$(ALICE_ROOT)/CONTAINERS -I$(ALICE_ROOT)/STEER
DEFSTR = -Duse_aliroot -Duse_root
endif

#Save the particle id's
ifdef DO_MC
DEFSTR += -Ddo_mc
endif
DEFSTR += $(EXTRADEF)

#Use logging classes
ifndef NOLOGGING
DEFSTR += -Duse_logging
endif

ifdef MLUCDIR
INCLUDES += -I$(MLUCDIR)/include
else
INCLUDES += -I/prog/alice/level3/kip/MLUC/include
endif

SRCS	= AliL3HoughTransformer.cxx AliL3Hough.cxx AliL3HoughTrack.cxx\
	  AliL3HoughMaxFinder.cxx AliL3HoughEval.cxx AliL3HoughMerger.cxx \
	  AliL3Histogram.cxx AliL3Histogram1D.cxx AliL3HoughBaseTransformer.cxx \
	  AliL3HoughIntMerger.cxx AliL3HoughGlobalMerger.cxx \
	  AliL3HoughTransformerVhdl.cxx AliL3HoughDisplay.cxx

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
	rm -f $(DICT) $(DICTH) 

so: 
	rm -f $(LIBDIR)/libAliL3Hough.so
