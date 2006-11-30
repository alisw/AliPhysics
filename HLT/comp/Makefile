#$Id$
################################
# Makefile for HLT Comp code.  #
#                              #
# Author: Anders Vestbo        #                    
#         Constantin Loizides  #
################################

MODNAME = AliHLTComp

ifndef ALIHLT_TOPDIR
ALIHLT_TOPDIR = $(shell pwd)/..
endif

include $(ALIHLT_TOPDIR)/Makefile.conf

SRCS =  AliHLTModeller.cxx AliHLTClusterFitter.cxx \
        AliHLTCompress.cxx  AliHLTDataCompressor.cxx AliHLTCompressAC.cxx

ifneq ($(ARCH),macosx)
## for macosx put into src as symbolic link
SRCS += AliHLTDataCompressorHelper.cxx AliHLTModelTrack.cxx
endif

ADD_C_OBJS = $(OBJDIR)/bitio.o $(OBJDIR)/errhand.o $(OBJDIR)/AliHLTFitUtilities.o

ifeq ($(ALIHLT_USEPACKAGE),ROOT)
HDRS = $(SRCS:.cxx=.h) $(MODNAME)LinkDef.h
OBJS = $(patsubst %.cxx,$(OBJDIR)/%.o,$(SRCS)) $(ADD_C_OBJS) $(DICTO)
endif

ifeq ($(ALIHLT_USEPACKAGE),ALIROOT)
SRCS += AliHLTOfflineDataCompressor.cxx
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

$(OBJDIR)/AliHLTFitUtilities.o: AliHLTFitUtilities.c
	$(CXX) $(CXXFLAGS) -c -o $(OBJDIR)/AliHLTFitUtilities.o AliHLTFitUtilities.c

$(OBJDIR)/bitio.sto: bitio.c
	$(CXX) $(PROFILEFLAGS) $(CXXFLAGS) -c -o $(OBJDIR)/bitio.sto bitio.c

$(OBJDIR)/errhand.sto: errhand.c
	$(CXX) $(PROFILEFLAGS) $(CXXFLAGS) -c -o $(OBJDIR)/errhand.sto errhand.c

$(OBJDIR)/AliHLTFitUtilities.sto: AliHLTFitUtilities.c
	$(CXX) $(PROFILEFLAGS) $(CXXFLAGS) -c -o $(OBJDIR)/AliHLTFitUtilities.sto AliHLTFitUtilities.c
