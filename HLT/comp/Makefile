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
INCLUDES = -I$(ROOTSYS)/include -I$(TOPDIR)/hough -I$(TOPDIR)/src -I$(TOPDIR)/comp
DEFSTR = -Duse_root
OBJDIR = lib_ROOT
endif

ifeq ($(USEPACKAGE),ALIROOT)
INCLUDES = -I$(ROOTSYS)/include -I$(TOPDIR)/hough -I$(TOPDIR)/src -I$(TOPDIR)/comp -I$(ALICE_ROOT)/include/ -I$(ALICE_ROOT)/TPC -I$(ALICE_ROOT)/CONTAINERS -I$(ALICE_ROOT)/STEER
DEFSTR = -Duse_aliroot -Duse_root
endif

ifdef DO_MC
DEFSTR += -Ddo_mc
endif

ifndef NOLOGGING
DEFSTR += -Duse_logging
endif

ifdef MLUCDIR
INCLUDES += -I$(MLUCDIR)/include
else
INCLUDES += -I/prog/alice/level3/kip/MLUC/include
endif

SRCS	= AliL3Modeller.cxx AliL3ModelTrack.cxx AliL3Compress.cxx

DICT = AliL3CompCint.cxx
DICTH = AliL3CompCint.h
DICTO = $(OBJDIR)/AliL3CompCint.o

HDRS = $(SRCS:.cxx=.h) bitio.h errhand.h AliL3CompLinkDef.h

OBJS = $(patsubst %.cxx,$(OBJDIR)/%.o,$(SRCS)) $(OBJDIR)/bitio.o $(OBJDIR)/errhand.o $(DICTO)

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