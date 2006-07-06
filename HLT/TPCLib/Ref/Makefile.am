# $Id$
# Makefile template for the Alice HLT TPC library

MODULE 				= AliHLTTPCRef

AM_CPPFLAGS			= -I$(top_srcdir)/src        \
				  -I$(top_srcdir)/BASE       \
				  -I$(ALICE_ROOT)/include    \
				  -I$(ALICE_ROOT)/TPC        \
				  -I$(ALICE_ROOT)/CONTAINERS \
				  -I$(ALICE_ROOT)/STEER      \
				  -I$(ALICE_ROOT)/RAW


# library definition
lib_LTLIBRARIES			=  libAliHLTTPCRef.la

# library sources
libAliHLTTPCRef_la_SOURCES	= bitio.cxx \
				  errhand.cxx \
				  AliHLTTPCLog.cxx \
				  AliHLTTPCTransform.cxx \
				  AliHLTTPCMemHandler.cxx \
			          AliHLTTPCCompress.cxx \
			          AliHLTTPCDataCompressor.cxx \
				  AliHLTTPCClusterFitter.cxx \
				  AliHLTTPCModeller.cxx \
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
			          AliHLTTPCCompress.h \
			          AliHLTTPCDataCompressor.h \
				  AliHLTTPCClusterFitter.h \
				  AliHLTTPCModeller.h \
				  AliHLTTPCDataCompressorHelper.h \
				  AliHLTTPCFitter.h \
				  AliHLTTPCFileHandler.h \
				  AliHLTTPCBenchmark.h \
				  AliHLTTPCDDLDataFileHandler.h \
				  AliHLTTPCClustFinderNew.h \
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
			          AliHLTTPCHoughTrack.h \
				  AliHLTTPCMerger.h \
				  AliHLTTPCTrackMerger.h \
				  AliHLTTPCGlobalMerger.h \
				  AliHLTTPCInterMerger.h \
				  AliHLTTPC.h \
				  AliHLTTPCDisplay.h \
				  AliHLTTPCDefinitions.h \
				  AliHLTTPCRawDataUnpackerComponent.h \
				  AliHLTTPCClusterFinderComponent.h \
				  AliHLTTPCVertexFinderComponent.h \
				  AliHLTTPCSliceTrackerComponent.h \
				  AliHLTTPCGlobalMergerComponent.h


pkginclude_HEADERS		= $(CLASS_HDRS) \
				  bitio.h \
				  errhand.h \
				  AliHLTTPCLogging.h \
				  AliHLTTPCStandardIncludes.h \
				  AliHLTTPCRootTypes.h \
				  AliHLTTPCDigitData.h \
				  AliHLTTPCTrackSegmentData.h \
				  AliHLTTPCSpacePointData.h \
				  AliHLTTPCVertexData.h \
				  AliHLTTPCTrackletDataFormat.h \
				  AliHLTTPCRawDataFormat.h \
				  AliHLTTPCClusterDataFormat.h \
				  AliHLTTPCModels.h \
				  AliHLTTPCFitUtilities.h


# version info for the library
libAliHLTTPCRef_la_LDFLAGS	=  -version-info 1:0:0

# set the file name for the generated root dictionary
DICTCPP				=  AliHLTTPCRef-DICT.cxx
nodist_libAliHLTTPCRef_la_SOURCES  =  $(DICTCPP)

CLEANFILES			=

include $(top_srcdir)/make.dict
