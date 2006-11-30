#$Id$
################## ##############
# Makefile for HLT Hough code. #
#                              #
# Author: Anders Vestbo,       #                    
#         Constantin Loizides  #
################################

MODNAME = AliHLTHough

ifndef ALIHLT_TOPDIR
ALIHLT_TOPDIR = $(shell pwd)/..
endif

include $(ALIHLT_TOPDIR)/Makefile.conf

SRCS = AliHLTHoughTransformer.cxx AliHLTHoughClusterTransformer.cxx \
       AliHLTHoughEval.cxx AliHLTHoughMerger.cxx AliHLTHoughBaseTransformer.cxx \
       AliHLTHoughIntMerger.cxx AliHLTHoughGlobalMerger.cxx AliHLTHoughTransformerVhdl.cxx \
       AliHLTHistogram.cxx AliHLTHistogram1D.cxx AliHLTHoughMaxFinder.cxx AliHLTHough.cxx \
       AliHLTHoughTransformerLUT.cxx AliHLTHoughTransformerRow.cxx AliHLTHistogramAdaptive.cxx 

ifneq ($(ARCH),macosx)
## for macosx put into src as symbolic link
SRCS +=  AliHLTHoughTrack.cxx 
endif

ifeq ($(ALIHLT_USEPACKAGE),ROOT)
SRCS += AliHLTHoughDisplay.cxx AliHLTHoughTest.cxx AliHLTHoughTransformerNew.cxx AliHLTHoughTransformerGlobal.cxx
HDRS = $(SRCS:.cxx=.h) $(MODNAME)LinkDef.h
OBJS = $(patsubst %.cxx,$(OBJDIR)/%.o,$(SRCS)) $(DICTO)
endif

ifeq ($(ALIHLT_USEPACKAGE),ALIROOT)
SRCS += AliHLTHoughDisplay.cxx AliHLTHoughTest.cxx AliHLTHoughTransformerNew.cxx AliHLTHoughTransformerGlobal.cxx
HDRS = $(SRCS:.cxx=.h) $(MODNAME)LinkDef.h
OBJS = $(patsubst %.cxx,$(OBJDIR)/%.o,$(SRCS)) $(DICTO)
endif

include $(ALIHLT_TOPDIR)/Makefile.rules
