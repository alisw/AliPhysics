#$Id$
##########################################
# Makefile for AliRoot HLT Sample        #
#                                        #
# Author: Matthias Richter,              #                     
##########################################

MODNAME = AliHLTTPCRef

ifndef ALIHLT_TOPDIR
ALIHLT_TOPDIR = $(shell pwd)/../..
endif

USENEWIO=1
#EXTRACXXFLAGS=-ggdb3 -O0
#EXTRALDFLAGS=-ggdb3 -O0

include $(ALIHLT_TOPDIR)/Makefile.conf

SRCS	= bitio.cxx \
          errhand.cxx \
          AliHLTTPCLog.cxx \
          AliHLTTPCTransform.cxx \
          AliHLTTPCMemHandler.cxx \
          AliHLTTPCCompress.cxx \
          AliHLTTPCDataCompressor.cxx \
          AliHLTTPCDataCompressorHelper.cxx \
          AliHLTTPCFitter.cxx \
          AliHLTTPCFileHandler.cxx \
          AliHLTTPCBenchmark.cxx \
          AliHLTTPCDDLDataFileHandler.cxx \
          AliHLTTPCClustFinderNew.cxx \
          AliHLTTPCVertex.cxx \
          AliHLTTPCVertexArray.cxx \
          AliHLTTPCVertexFinder.cxx \
          AliHLTTPCTrackArray.cxx \
          AliHLTTPCTrack.cxx \
          AliHLTTPCConfMapFit.cxx \
          AliHLTTPCConfMapTrack.cxx \
          AliHLTTPCConfMapPoint.cxx \
          AliHLTTPCConfMapper.cxx \
          AliHLTTPCModelTrack.cxx \
          AliHLTTPCHoughTrack.cxx \
          AliHLTTPCMerger.cxx \
          AliHLTTPCTrackMerger.cxx \
          AliHLTTPCGlobalMerger.cxx \
          AliHLTTPCInterMerger.cxx \
          AliHLTTPC.cxx \
          AliHLTTPCDisplay.cxx \
          AliHLTTPCDefinitions.cxx \
          AliHLTTPCRawDataUnpackerComponent.cxx \
          AliHLTTPCClusterFinderComponent.cxx \
          AliHLTTPCVertexFinderComponent.cxx \
          AliHLTTPCSliceTrackerComponent.cxx \
          AliHLTTPCGlobalMergerComponent.cxx \
 



          #AliRawReaderMemory.cxx

          #AliHLTTPCHistogram.cxx \
          #AliHLTTPCHistogramAdaptive.cxx \
          #AliHLTTPCHoughBaseTransformer.cxx \
          #AliHLTTPCHoughTransformerRow.cxx \

EXTRA_INCLUDES = $(MODNAME)StandardIncludes.h AliHLTTPCSpacePointData.h

SRCS += $(EXTRA_SRCS)


include $(ALIHLT_TOPDIR)/Makefile.complib.rules
