#$Id$
#################################
# Makefile for HLT Src code.    #
#                               #
# Author: Anders Vestbo,        #                     
#         Constantin Loizides   #
#################################

MODNAME = AliL3Src

ifndef ALIHLT_TOPDIR
ALIHLT_TOPDIR = $(shell pwd)/..
endif

include $(ALIHLT_TOPDIR)/Makefile.conf

SRCS	= AliL3ConfMapper.cxx AliL3ConfMapPoint.cxx \
          AliL3ConfMapTrack.cxx AliL3ConfMapFit.cxx AliL3Track.cxx\
          AliL3Transform.cxx AliL3TrackMerger.cxx AliL3VertexArray.cxx \
          AliL3Vertex.cxx AliL3VertexFinder.cxx \
          AliL3Merger.cxx AliL3GlobalMerger.cxx AliL3ClustFinderNew.cxx \
          AliL3TrackArray.cxx AliL3InterMerger.cxx AliL3Logger.cxx \
          AliL3MemHandler.cxx AliLevel3.cxx AliL3Benchmark.cxx  \
          AliL3Fitter.cxx AliL3RawDataFileHandler.cxx AliL3Log.cxx \
          AliL3TPCBeamTestMemHandler.cxx

ifeq ($(ARCH),macosx)
EXTRA_SRCS = AliL3HoughTrack.cxx AliL3ModelTrack.cxx \
             AliL3DataCompressorHelper.cxx AliL3DataHandler.cxx \
             AliL3TransBit.cxx  AliL3DDLDataFileHandler.cxx
endif

SRCS += $(EXTRA_SRCS)

ifeq ($(ALIHLT_USEPACKAGE),ROOT)
SRCS += AliL3Display.cxx
HDRS = $(SRCS:.cxx=.h) AliL3SpacePointData.h $(MODNAME)LinkDef.h
OBJS = $(patsubst %.cxx,$(OBJDIR)/%.o,$(SRCS)) $(DICTO)
endif

ifeq ($(ALIHLT_USEPACKAGE),ALIROOT)
SRCS += AliL3Display.cxx AliL3FileHandler.cxx AliL3Evaluate.cxx AliHLTReconstructor.cxx
HDRS = $(SRCS:.cxx=.h) AliL3SpacePointData.h $(MODNAME)LinkDef.h
OBJS = $(patsubst %.cxx,$(OBJDIR)/%.o,$(SRCS)) $(DICTO)
endif

include $(ALIHLT_TOPDIR)/Makefile.rules
