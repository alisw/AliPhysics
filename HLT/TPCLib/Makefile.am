# $Id$
# Makefile template for the Alice HLT TPC library

MODULE 				= AliHLTTPC

AM_CPPFLAGS			= -I$(top_srcdir)/BASE       \
				  -I$(ALICE_ROOT)/include    \
				  -I$(ALICE_ROOT)/TPC        \
				  -I$(ALICE_ROOT)/CONTAINERS \
				  -I$(ALICE_ROOT)/STEER      \
				  -I$(ALICE_ROOT)/RAW


# library definition
lib_LTLIBRARIES			=  libAliHLTTPC.la

# library sources
libAliHLTTPC_la_SOURCES		= AliHLTTPCLog.cxx \
				  AliHLTTPCTransform.cxx \
				  AliHLTTPCMemHandler.cxx \
				  AliHLTTPCDataCompressorHelper.cxx \
				  AliHLTTPCFitter.cxx \
				  AliHLTTPCFileHandler.cxx \
				  AliHLTTPCBenchmark.cxx \
				  AliHLTTPCDDLDataFileHandler.cxx \
				  AliHLTTPCClustFinderNew.cxx \
				  AliHLTTPCClusterFinder.cxx \
				  AliHLTTPCDigitReader.cxx \
				  AliHLTTPCDigitReaderPacked.cxx \
				  AliHLTTPCDigitReaderUnpacked.cxx \
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
				  AliHLTTPCGlobalMergerComponent.cxx


#				  AliRawReaderMemory.cxx


#				  AliHLTTPCHoughTrack.cxx \
#				  AliHLTTPCHistogram.cxx \
#				  AliHLTTPCHistogramAdaptive.cxx \
#				  AliHLTTPCHoughBaseTransformer.cxx \
#				  AliHLTTPCHoughTransformerRow.cxx


# class header files, the link definition for the root dictionary
# will be created from the names of the header files
CLASS_HDRS			= AliHLTTPCLog.h \
				  AliHLTTPCTransform.h \
				  AliHLTTPCMemHandler.h \
				  AliHLTTPCDataCompressorHelper.h \
				  AliHLTTPCFitter.h \
				  AliHLTTPCFileHandler.h \
				  AliHLTTPCBenchmark.h \
				  AliHLTTPCDDLDataFileHandler.h \
				  AliHLTTPCClustFinderNew.h \
				  AliHLTTPCClusterFinder.h \
				  AliHLTTPCDigitReader.h \
				  AliHLTTPCDigitReaderPacked.h \
				  AliHLTTPCDigitReaderUnpacked.h \
				  AliHLTTPCVertex.h \
				  AliHLTTPCVertexArray.h \
				  AliHLTTPCVertexFinder.h \
				  AliHLTTPCTrackArray.h \
				  AliHLTTPCTrack.h \
				  AliHLTTPCConfMapFit.h \
				  AliHLTTPCConfMapTrack.h \
				  AliHLTTPCConfMapPoint.h \
				  AliHLTTPCConfMapper.h \
				  AliHLTTPCModelTrack.h \
				  AliHLTTPCMerger.h \
				  AliHLTTPCTrackMerger.h \
				  AliHLTTPCGlobalMerger.h \
				  AliHLTTPCInterMerger.h \
				  AliHLTTPC.h \
				  AliHLTTPCSpacePointData.h \
				  AliHLTTPCDisplay.h \
				  AliHLTTPCDefinitions.h \
				  AliHLTTPCRawDataUnpackerComponent.h \
				  AliHLTTPCClusterFinderComponent.h \
				  AliHLTTPCVertexFinderComponent.h \
				  AliHLTTPCSliceTrackerComponent.h \
				  AliHLTTPCGlobalMergerComponent.h


pkginclude_HEADERS		= $(CLASS_HDRS) \
				  AliHLTTPCLogging.h \
				  AliHLTTPCStandardIncludes.h \
				  AliHLTTPCRootTypes.h \
				  AliHLTTPCDigitData.h \
				  AliHLTTPCTrackSegmentData.h \
				  AliHLTTPCVertexData.h \
				  AliHLTTPCTrackletDataFormat.h \
				  AliHLTTPCRawDataFormat.h \
				  AliHLTTPCClusterDataFormat.h \
				  AliHLTTPCModels.h


# version info for the library
libAliHLTTPC_la_LDFLAGS	=  -version-info 2:0:0

# set the file name for the generated root dictionary
DICTCPP				=  AliHLTTPC-DICT.cxx
nodist_libAliHLTTPC_la_SOURCES  =  $(DICTCPP)

CLEANFILES			=

include $(top_srcdir)/make.dict

if EN_HLT_TPCREF
REFDIR				= Ref
endif
SUBDIRS 			= . $(REFDIR)
