#$Id$
################################
# Makefile for HLT Comp code.  #
#                              #
# Author: Anders Vestbo        #                    
#         Constantin Loizides  #
################################

MODNAME = AliL3Comp

ifndef ALIHLT_TOPDIR
ALIHLT_TOPDIR = $(shell pwd)/..
endif

include $(ALIHLT_TOPDIR)/Makefile.conf

SRCS =  AliL3Modeller.cxx AliL3ClusterFitter.cxx \
        AliL3Compress.cxx  AliL3DataCompressor.cxx AliL3CompressAC.cxx

ifneq ($(ARCH),Darwin)
SRCS += AliL3DataCompressorHelper.cxx AliL3ModelTrack.cxx
endif

ADD_C_OBJS = $(OBJDIR)/bitio.o $(OBJDIR)/errhand.o $(OBJDIR)/AliL3FitUtilities.o

ifeq ($(ALIHLT_USEPACKAGE),ROOT)
HDRS = $(SRCS:.cxx=.h) $(MODNAME)LinkDef.h
OBJS = $(patsubst %.cxx,$(OBJDIR)/%.o,$(SRCS)) $(ADD_C_OBJS) $(DICTO)
endif

ifeq ($(ALIHLT_USEPACKAGE),ALIROOT)
SRCS += AliL3OfflineDataCompressor.cxx
HDRS = $(SRCS:.cxx=.h) $(MODNAME)LinkDef.h
OBJS = $(patsubst %.cxx,$(OBJDIR)/%.o,$(SRCS)) $(ADD_C_OBJS) $(DICTO)
endif

ifeq ($(ALIHLT_USEPACKAGE),STANDALONE)
OBJS = $(patsubst %.cxx,$(OBJDIR)/%.o,$(SRCS)) $(ADD_C_OBJS)
endif

include $(ALIHLT_TOPDIR)/Makefile.rules


$(OBJDIR)/bitio.o: bitio.c
	$(CXX) $(CXXFLAGS) -c -o $(OBJDIR)/bitio.o bitio.c

$(OBJDIR)/errhand.o: errhand.c
	$(CXX) $(CXXFLAGS) -c -o $(OBJDIR)/errhand.o errhand.c

$(OBJDIR)/AliL3FitUtilities.o: AliL3FitUtilities.c
	$(CXX) $(CXXFLAGS) -c -o $(OBJDIR)/AliL3FitUtilities.o AliL3FitUtilities.c

$(OBJDIR)/bitio.sto: bitio.c
	$(CXX) $(PROFILEFLAGS) $(CXXFLAGS) -c -o $(OBJDIR)/bitio.sto bitio.c

$(OBJDIR)/errhand.sto: errhand.c
	$(CXX) $(PROFILEFLAGS) $(CXXFLAGS) -c -o $(OBJDIR)/errhand.sto errhand.c

$(OBJDIR)/AliL3FitUtilities.sto: AliL3FitUtilities.c
	$(CXX) $(PROFILEFLAGS) $(CXXFLAGS) -c -o $(OBJDIR)/AliL3FitUtilities.sto AliL3FitUtilities.c
