#$Id$
################################
# Makefile for HLT Hough code. #
#                              #
# Author: Anders Vestbo,       #                    
#         Constantin Loizides  #
################################

MODNAME = AliL3Hough

ifndef ALIHLT_TOPDIR
ALIHLT_TOPDIR = $(shell pwd)/..
endif

include $(ALIHLT_TOPDIR)/Makefile.conf

SRCS = AliL3HoughTransformer.cxx AliL3HoughClusterTransformer.cxx \
       AliL3HoughEval.cxx AliL3HoughMerger.cxx AliL3HoughBaseTransformer.cxx \
       AliL3HoughIntMerger.cxx AliL3HoughGlobalMerger.cxx AliL3HoughTransformerVhdl.cxx \
       AliL3Histogram.cxx AliL3Histogram1D.cxx AliL3HoughMaxFinder.cxx AliL3Hough.cxx \
       AliL3HoughTransformerLUT.cxx AliL3HoughTransformerRow.cxx AliL3HistogramAdaptive.cxx 

ifeq ($(ARCH),macosx)
## AliL3HoughTrack put into src as symbolic link
EXTRA_SRCS = AliL3HoughDisplay.cxx AliL3HoughTransformerNew.cxx
else
EXTRA_SRCS = AliL3HoughDisplay.cxx AliL3HoughTransformerNew.cxx  AliL3HoughTrack.cxx AliL3HoughTest.cxx
endif

ifeq ($(ALIHLT_USEPACKAGE),ROOT)
SRCS += $(EXTRA_SRCS)
HDRS = $(SRCS:.cxx=.h) $(MODNAME)LinkDef.h
OBJS = $(patsubst %.cxx,$(OBJDIR)/%.o,$(SRCS)) $(DICTO)
endif

ifeq ($(ALIHLT_USEPACKAGE),ALIROOT)
SRCS += $(EXTRA_SRCS) AliL3HoughTransformerGlobal.cxx
HDRS = $(SRCS:.cxx=.h) $(MODNAME)LinkDef.h
OBJS = $(patsubst %.cxx,$(OBJDIR)/%.o,$(SRCS)) $(DICTO)
endif

include $(ALIHLT_TOPDIR)/Makefile.rules
