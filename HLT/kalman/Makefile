#$Id$
###########################
# Makefile for HLT code.  #
#                         #
# Author: Anders Vestbo   #                    
# changes by C. Loizides  #
###########################


ifndef TOPDIR
TOPDIR = $(shell pwd)/..
endif

include $(TOPDIR)/Makefile.conf

DEFSTR = -Dno_root
OBJDIR = lib_alone
INCLUDES = -I$(TOPDIR)/hough -I$(TOPDIR)/src -I$(TOPDIR)/comp -I$(TOPDIR)/misc -I$(TOPDIR)/kalman

ifeq ($(USEPACKAGE),ROOT) 
INCLUDES = -I$(ROOTSYS)/include -I$(TOPDIR)/hough -I$(TOPDIR)/src -I$(TOPDIR)/comp -I$(TOPDIR)/misc -I$(TOPDIR)/kalman
DEFSTR = -Duse_root
OBJDIR = lib_ROOT
endif

ifeq ($(USEPACKAGE),ALIROOT)
INCLUDES = -I$(ROOTSYS)/include -I$(TOPDIR)/hough -I$(TOPDIR)/src -I$(TOPDIR)/comp -I$(TOPDIR)/misc -I$(TOPDIR)/kalman -I$(ALICE_ROOT)/include/ -I$(ALICE_ROOT)/TPC -I$(ALICE_ROOT)/CONTAINERS -I$(ALICE_ROOT)/STEER 
DEFSTR = -Duse_aliroot -Duse_root
OBJDIR = lib
endif

#Save the particle id's
ifdef DO_MC
DEFSTR += -Ddo_mc
endif

#Use logging classes
ifndef NOLOGGING
DEFSTR += -Duse_logging
endif

DEFSTR += $(EXTRADEF)

ifdef MLUCDIR
INCLUDES += -I$(MLUCDIR)/include
else
INCLUDES += -I/prog/alice/level3/kip/MLUC/include
endif

SRCS = AliL3Kalman.cxx AliL3KalmanTrack.cxx AliL3KalmanDisplay.cxx
DICT  = AliL3KalmanCint-$(USEPACKAGE).cxx
DICTH = AliL3KalmanCint-$(USEPACKAGE).h
DICTO = $(OBJDIR)/AliL3KalmanCint-$(USEPACKAGE).o

ifeq ($(OBJDIR),lib_alone)
HDRS = $(SRCS:.cxx=.h) 
OBJS = $(patsubst %.cxx,$(OBJDIR)/%.o,$(SRCS))
endif

ifeq ($(USEPACKAGE),ROOT)
HDRS = $(SRCS:.cxx=.h) AliL3KalmanLinkDef.h
OBJS = $(patsubst %.cxx,$(OBJDIR)/%.o,$(SRCS)) $(DICTO)
endif

ifeq ($(USEPACKAGE),ALIROOT)
HDRS = $(SRCS:.cxx=.h) AliL3KalmanLinkDef.h
OBJS = $(patsubst %.cxx,$(OBJDIR)/%.o,$(SRCS)) $(DICTO)
endif


LEVELSO = $(LIBDIR)/libAliL3Kalman.so

all : $(OBJDIR) $(LIBDIR) $(LEVELSO)

$(LEVELSO) : $(OBJS)
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ -o $@

$(DICT) : $(HDRS)
	@echo "Generating dictionary..."
	@rootcint -f $(DICT) -c $(INCLUDES) -DGCCVERSION=$(GCCVERSION) $(DEFS) -include AliL3StandardIncludes.h $(HDRS)

$(OBJDIR)/%.o : %.cxx 
	$(CXX) $(CXXFLAGS) $(DEFSTR) -c $(INCLUDES) -o $@ $<

$(OBJDIR) : 
	test ! -e $(OBJDIR) && mkdir -p $(OBJDIR)

$(LIBDIR) :
	test ! -e $(LIBDIR) && mkdir -p $(LIBDIR)

clean :
	rm -f $(OBJDIR)/*.o
	rm -f $(LIBDIR)/libAliL3Kalman.so
	rm -f $(DICT) $(DICTH) 

so: 
	rm -f $(LIBDIR)/libAliL3Kalman.so
