#$Id$
###########################
# Makefile for HLT code.  #
#                         #
# Author: Anders Vestbo   #                    
# changes by C. Loizides  #
###########################

ifndef TOPDIR
TOPDIR = `pwd`/..
endif

include $(TOPDIR)/Makefile.conf

DEFSTR = -Dno_root
OBJDIR = lib_alone
INCLUDES = -I$(TOPDIR)/hough -I$(TOPDIR)/src -I$(TOPDIR)/comp -I$(TOPDIR)/misc

ifeq ($(USEPACKAGE),ROOT) 
INCLUDES = -I$(ROOTSYS)/include -I$(TOPDIR)/hough -I$(TOPDIR)/src -I$(TOPDIR)/comp
DEFSTR = -Duse_root
OBJDIR = lib_ROOT
endif

ifeq ($(USEPACKAGE),ALIROOT)
INCLUDES = -I$(ROOTSYS)/include -I$(TOPDIR)/hough -I$(TOPDIR)/src -I$(TOPDIR)/comp -I$(ALICE_ROOT)/include/ -I$(ALICE_ROOT)/TPC -I$(ALICE_ROOT)/CONTAINERS -I$(ALICE_ROOT)/STEER
DEFSTR = -Duse_aliroot
OBJDIR = lib
endif

ifdef DO_MC
DEFSTR += -Ddo_mc
endif

ifndef NOLOGGING
DEFSTR += -Duse_logging
endif

DEFSTR += $(EXTRADEF)

ifdef MLUCDIR
INCLUDES += -I$(MLUCDIR)/include
else
INCLUDES += -I/prog/alice/level3/kip/MLUC/include
endif

SRCS	= AliTransBit.cxx AliL3AltroMemHandler.cxx AliL3DataHandler.cxx AliL3VHDLClusterFinder.cxx \
          AliL3FFloat.cxx

DICT  = AliL3MiscCint-$(USEPACKAGE).cxx
DICTH = AliL3MiscCint-$(USEPACKAGE).h
DICTO = $(OBJDIR)/AliL3MiscCint-$(USEPACKAGE).o

HDRS = $(SRCS:.cxx=.h) AliL3MiscLinkDef.h
OBJS = $(patsubst %.cxx,$(OBJDIR)/%.o,$(SRCS)) $(DICTO)

ifeq ($(OBJDIR),lib_alone)
HDRS = $(SRCS:.cxx=.h) 
OBJS = $(patsubst %.cxx,$(OBJDIR)/%.o,$(SRCS))
endif


LEVELSO = $(LIBDIR)/libAliL3Misc.so

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
	rm -f $(LIBDIR)/libAliL3Misc.so
	rm -f $(DICT) $(DICTH) 

so :
	rm -f $(LIBDIR)/libAliL3Misc.so
