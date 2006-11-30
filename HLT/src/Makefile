#$Id$
#################################
# Makefile for HLT Src code.    #
#                               #
# Author: Anders Vestbo,        #                     
#         Constantin Loizides   #
#################################

MODNAME = AliHLTSrc

ifndef ALIHLT_TOPDIR
ALIHLT_TOPDIR = $(shell pwd)/..
endif

include $(ALIHLT_TOPDIR)/Makefile.conf

SRCS	= AliHLTConfMapper.cxx AliHLTConfMapPoint.cxx \
          AliHLTConfMapTrack.cxx AliHLTConfMapFit.cxx AliHLTTrack.cxx\
          AliHLTTransform.cxx AliHLTTrackMerger.cxx AliHLTVertexArray.cxx \
          AliHLTVertex.cxx AliHLTVertexFinder.cxx \
          AliHLTMerger.cxx AliHLTGlobalMerger.cxx AliHLTClustFinderNew.cxx \
          AliHLTTrackArray.cxx AliHLTInterMerger.cxx AliHLTLogger.cxx \
          AliHLTMemHandler.cxx AliLevel3.cxx AliHLTBenchmark.cxx  \
          AliHLTFitter.cxx AliHLTRawDataFileHandler.cxx AliHLTLog.cxx \
          AliHLTTPCBeamTestMemHandler.cxx

ifeq ($(ARCH),macosx)
EXTRA_SRCS = AliHLTHoughTrack.cxx AliHLTModelTrack.cxx \
             AliHLTDataCompressorHelper.cxx AliHLTDataHandler.cxx \
             AliHLTTransBit.cxx  AliHLTDDLDataFileHandler.cxx
endif

SRCS += $(EXTRA_SRCS)

ifeq ($(ALIHLT_USEPACKAGE),ROOT)
SRCS += AliHLTDisplay.cxx
HDRS = $(SRCS:.cxx=.h) AliHLTSpacePointData.h $(MODNAME)LinkDef.h
OBJS = $(patsubst %.cxx,$(OBJDIR)/%.o,$(SRCS)) $(DICTO)
endif

ifeq ($(ALIHLT_USEPACKAGE),ALIROOT)
SRCS += AliHLTDisplay.cxx AliHLTFileHandler.cxx AliHLTEvaluate.cxx AliHLTReconstructor.cxx
HDRS = $(SRCS:.cxx=.h) AliHLTSpacePointData.h $(MODNAME)LinkDef.h
OBJS = $(patsubst %.cxx,$(OBJDIR)/%.o,$(SRCS)) $(DICTO)
endif

include $(ALIHLT_TOPDIR)/Makefile.rules
