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
DEFSTR = -Duse_aliroot -Duse_root
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

SRCS =  AliL3Modeller.cxx AliL3ModelTrack.cxx 

DICT  = AliL3CompCint-$(USEPACKAGE).cxx
DICTH = AliL3CompCint-$(USEPACKAGE).h
DICTO = $(OBJDIR)/AliL3CompCint-$(USEPACKAGE).o

ifeq ($(USEPACKAGE),ROOT)
SRCS += AliL3Compress.cxx
HDRS = $(SRCS:.cxx=.h) bitio.h errhand.h AliL3CompLinkDef.h
OBJS = $(patsubst %.cxx,$(OBJDIR)/%.o,$(SRCS)) $(OBJDIR)/bitio.o $(OBJDIR)/errhand.o $(DICTO)
endif

ifeq ($(USEPACKAGE),ALIROOT)
SRCS += AliL3Compress.cxx
HDRS = $(SRCS:.cxx=.h) bitio.h errhand.h AliL3CompLinkDef.h
OBJS = $(patsubst %.cxx,$(OBJDIR)/%.o,$(SRCS)) $(OBJDIR)/bitio.o $(OBJDIR)/errhand.o $(DICTO)
endif

ifeq ($(OBJDIR),lib_alone)
HDRS = $(SRCS:.cxx=.h) bitio.h errhand.h
OBJS = $(patsubst %.cxx,$(OBJDIR)/%.o,$(SRCS)) $(OBJDIR)/bitio.o $(OBJDIR)/errhand.o
endif


LEVELSO = $(LIBDIR)/libAliL3Comp.so

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

$(OBJDIR)/bitio.o : bitio.c
	$(CXX) -c -o $(OBJDIR)/bitio.o bitio.c

$(OBJDIR)/errhand.o : errhand.c
	$(CXX) -c -o $(OBJDIR)/errhand.o errhand.c

clean :
	rm -f $(OBJDIR)/*.o
	rm -f $(LIBDIR)/libAliL3Comp.so
	rm -f $(DICT) $(DICTH)
so :
	rm -f $(LIBDIR)/libAliL3Comp.so
