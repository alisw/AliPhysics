# $Id$
# Makefile template for the Alice HLT TPC library

MODULE 				= AliHLTTPC

AM_CPPFLAGS			= @ALIROOT_CPPFLAGS@         \
				  -I$(top_srcdir)/BASE

# generation of the TPC Pad mapping tables for the 
# AliHLTTPCDigitReaderRaw reader
if USE_TPC_MAPPING
TPC_MAPPING_NOINST_SRC		= mapping2array.cxx
TPC_MAPPING_INC			= mapping_array_out.inc
TPC_MAPPING_NOINST_PGM		= mapping2array
TPC_MAPPING_COMPONENTS		= AliHLTTPCDigitReaderRaw.cxx
endif

noinst_PROGRAMS			= $(TPC_MAPPING_NOINST_PGM)
if USE_TPC_MAPPING
mapping2array_SOURCES		= $(TPC_MAPPING_NOINST_SRC)
endif

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
				  AliHLTTPCClusterFinder.cxx \
				  AliHLTTPCDigitReader.cxx \
				  AliHLTTPCDigitReaderPacked.cxx \
				  AliHLTTPCDigitReaderUnpacked.cxx \
				  AliHLTTPCDigitReaderRaw.cxx \
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
				  AliHLTTPCPad.cxx \
				  AliHLTTPCDefinitions.cxx \
				  AliHLTTPCRawDataUnpackerComponent.cxx \
				  AliHLTTPCClusterFinderComponent.cxx \
				  AliHLTTPCVertexFinderComponent.cxx \
				  AliHLTTPCSliceTrackerComponent.cxx \
				  AliHLTTPCGlobalMergerComponent.cxx


#				  AliHLTTPCClustFinderNew.cxx \
#				  AliHLTTPC.cxx \
# 				  AliHLTTPCDisplay.cxx \
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
				  AliHLTTPCClusterFinder.h \
				  AliHLTTPCDigitReader.h \
				  AliHLTTPCDigitReaderPacked.h \
				  AliHLTTPCDigitReaderUnpacked.h \
				  AliHLTTPCDigitReaderRaw.h \
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
				  AliHLTTPCPad.h \
				  AliHLTTPCSpacePointData.h \
				  AliHLTTPCDefinitions.h \
				  AliHLTTPCRawDataUnpackerComponent.h \
				  AliHLTTPCClusterFinderComponent.h \
				  AliHLTTPCVertexFinderComponent.h \
				  AliHLTTPCSliceTrackerComponent.h \
				  AliHLTTPCGlobalMergerComponent.h

#				  AliHLTTPCClustFinderNew.h \
#				  AliHLTTPC.h


pkginclude_HEADERS		= $(CLASS_HDRS) \
				  AliHLTTPCLogging.h \
				  AliHLTTPCRootTypes.h \
				  AliHLTTPCDigitData.h \
				  AliHLTTPCTrackSegmentData.h \
				  AliHLTTPCVertexData.h \
				  AliHLTTPCTrackletDataFormat.h \
				  AliHLTTPCRawDataFormat.h \
				  AliHLTTPCClusterDataFormat.h \
				  AliHLTTPCModels.h


# version info for the library
libAliHLTTPC_la_LDFLAGS		= -L@ROOTLIBDIR@ \
				  @ROOTLIBS@ \
				  @ALIROOT_LDFLAGS@ \
				  @ALIROOT_LIBS@ \
				  -version-info 3:0:0

# set the file name for the generated root dictionary
DICTCPP				=  AliHLTTPC-DICT.cxx
nodist_libAliHLTTPC_la_SOURCES  =  $(DICTCPP)

BUILT_SOURCES			= $(TPC_MAPPING_INC)
CLEANFILES			= $(BUILT_SOURCES)

include $(top_srcdir)/make.dict

SUBDIRS 			= . OnlineDisplay

$(TPC_MAPPING_INC): $(TPC_MAPPING_NOINST_PGM)
	./$<
