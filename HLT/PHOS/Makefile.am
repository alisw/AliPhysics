# $Id$
# Makefile template for the Alice HLT PHOS library

MODULE 				= AliHLTPHOS

EXTRA_DIST			= HLTPHOSLinkDef.h

AM_CPPFLAGS			= @ALIROOT_CPPFLAGS@ \
				  -I@ALICE_ROOT@/PHOS \
				  -I$(top_srcdir)/BASE \
				  -I$(top_srcdir)/BASE/util


# library definition
lib_LTLIBRARIES			=  libAliHLTPHOS.la

# library sources
libAliHLTPHOS_la_SOURCES	= AliHLTPHOSPulseGenerator.cxx \
			          AliHLTPHOSRawAnalyzerPeakFinder.cxx \
				  AliHLTPHOSRawAnalyzerCrude.cxx \
				  AliHLTPHOSRawAnalyzerKLevel.cxx \
		                  AliHLTPHOSRawAnalyzerLMS.cxx \
				  AliHLTPHOSRawAnalyzerChiSquareFit.cxx \
				  AliHLTPHOSRawAnalyzer.cxx \
				  AliHLTPHOSEmcCalibData.cxx \
				  AliHLTPHOSRawAnalyzerComponent.cxx \
				  AliHLTPHOSRawAnalyzerCrudeComponent.cxx \
				  AliHLTPHOSRawAnalyzerPeakFinderComponent.cxx \
			          AliHLTPHOSDefinitions.cxx \
				  AliHLTPHOSFileWriterComponent.cxx \
				  AliHLTPHOSFileWriter.cxx \
				  AliHLTPHOSDDLPackedFileWriter.cxx \
				  AliHLTPHOSCellEnergiesFileWriter.cxx \
				  AliHLTPHOSRcuHistogramProducer.cxx \
				  AliHLTPHOSRcuHistogramProducerComponent.cxx \
				  AliHLTPHOSDataCorruptor.cxx \
				  AliHLTPHOSDataQualityMonitor.cxx \
				  AliHLTPHOSRcuDataQualityMonitorComponent.cxx \
				  AliHLTPHOSProcessor.cxx \
				  AliHLTPHOSClusterizer.cxx \
				  AliHLTPHOSClusterizerComponent.cxx \
				  AliHLTPHOSPhysicsAnalyzerSpectrum.cxx \
				  AliHLTPHOSPhysicsAnalyzerSpectrumComponent.cxx \
				  AliHLTPHOSPhysicsAnalyzer.cxx \
				  AliHLTPHOSPhysicsAnalyzerPeakFitter.cxx \
			          AliHLTDDLDecoder.cxx \
				  AliHLTAltroData.cxx \
				  AliHLTAltroBunch.cxx \
			          AliHLTPHOSMapper.cxx \
				  AliHLTPHOSRcuProcessor.cxx \
				  AliHLTPHOSDDLDecoderComponent.cxx \
				  AliHLTPHOSBase.cxx



# class header files, the link definition for the root dictionary
# will be created from the names of the header files
CLASS_HDRS			= AliHLTPHOSPulseGenerator.h \
				  AliHLTPHOSRawAnalyzerPeakFinder.h \
				  AliHLTPHOSRawAnalyzerCrude.h \
				  AliHLTPHOSRawAnalyzerKLevel.h \
				  AliHLTPHOSRawAnalyzerLMS.h \
				  AliHLTPHOSRawAnalyzerChiSquareFit.h \
				  AliHLTPHOSEmcCalibData.h \
				  AliHLTPHOSRawAnalyzerComponent.h \
				  AliHLTPHOSRawAnalyzerCrudeComponent.h \
				  AliHLTPHOSRawAnalyzerPeakFinderComponent.h \
                                  AliHLTPHOSFileWriterComponent.h \
				  AliHLTPHOSFileWriter.h \
				  AliHLTPHOSDDLPackedFileWriter.h \
				  AliHLTPHOSCellEnergiesFileWriter.h \
				  AliHLTPHOSRcuHistogramProducer.h \
				  AliHLTPHOSRcuHistogramProducerComponent.h \
				  AliHLTPHOSDataCorruptor.h \
				  AliHLTPHOSDataQualityMonitor.h \
				  AliHLTPHOSRcuDataQualityMonitorComponent.h \
				  AliHLTPHOSProcessor.h \
				  AliHLTPHOSClusterizer.h \
				  AliHLTPHOSClusterizerComponent.h \
				  AliHLTPHOSPhysicsAnalyzerSpectrum.h \
				  AliHLTPHOSPhysicsAnalyzerSpectrumComponent.h \
				  AliHLTPHOSPhysicsAnalyzer.h \
				  AliHLTPHOSPhysicsAnalyzerPeakFitter.h \
				  AliHLTDDLDecoder.h \
				  AliHLTAltroData.h \
				  AliHLTPHOSMapper.h \
				  AliHLTPHOSRcuProcessor.h \
				  AliHLTPHOSBase.h \
				  AliHLTPHOSDDLDecoderComponent.h \
				  AliHLTPHOSMapper.h

pkginclude_HEADERS		= $(CLASS_HDRS) \
				  AliHLTPHOSCommonDefs.h \
				  AliHLTPHOSConstants.h \
				  AliHLTAltroBunch.h \
				  AliHLTPHOSRecPointListDataStruct.h \
				  AliHLTPHOSRecPointDataStruct.h \
				  AliHLTPHOSClusterDataStruct.h \
	 			  AliHLTPHOSRawAnalyzer.h \
			          AliHLTPHOSDefinitions.h \
				  AliHLTPHOSModuleCellEnergyDataStruct.h\
				  AliHLTPHOSRcuCellEnergyDataStruct.h \
				  AliHLTPHOSValidCellDataStruct.h\
				  AliHLTPHOSDataHeaderStruct.h \
				  AliHLTPHOSModuleCellAccumulatedEnergyDataStruct.h \
				  AliHLTPHOSModuleCellAverageEnergyDataStruct.h \
				  AliHLTPHOSRcuCellAccumulatedEnergyDataStruct.h \
			 	  AliHLTPHOSRcuCellAverageEnergyDataStruct.h \
				  AliHLTPHOSRcuChannelDataStruct.h \
				  AliHLTPHOSValidChannelDataStruct.h
				  AliHLTPHOSClusterDataStruct.h \
				  AliHLTPHOSRecPointDataStruct.h \
				  AliHLTPHOSPhysicsDefinitions.h \
				  AliHLTPHOSRecPointListDataStruct.h

# version info and linking flags for the library
libAliHLTPHOS_la_LIBADD		= $(top_builddir)/BASE/util/libAliHLTUtil.la
libAliHLTPHOS_la_LDFLAGS	= -L@ROOTLIBDIR@ \
				  @ROOTLIBS@ \
				  @ALIROOT_LDFLAGS@ \
				  @ALIROOT_LIBS@ \
				  -lRAWDatasim \
				  -lPHOSbase \
				  -lPHOSsim \
				  -lPHOSrec \
				  -version-info 0:0:0

# automatic generation of data and time of library build
COMPILE_INFO			=  AliHLTPHOSCompileInfo.cxx

# set the file name for the generated root dictionary
DICTCPP				=  AliHLTPHOS-DICT.cxx
# add additional include files which are necessary for the compilation of the 
# dictionary files and which are not part of the CLASS_HDRS
DICTINCLUDE			= 

nodist_libAliHLTPHOS_la_SOURCES = $(COMPILE_INFO) \
				  $(DICTCPP)

CLEANFILES			= $(COMPILE_INFO)

include $(top_srcdir)/make.dict

$(COMPILE_INFO): $(libAliHLTPHOS_la_SOURCES) $(pkginclude_HEADERS) $(noinst_HEADERS) Makefile.am
	@echo '//automatically generated compilation info' > $@
	@echo '//!!! DO NOT EDIT THIS FILE !!!' >> $@
	@echo '//add changes in Makefile.am' >> $@
	@echo '#include <iostream>' >> $@
	@echo 'extern "C" void CompileInfo( char*& date, char*& time)' >> $@
	@echo '{date=__DATE__; time=__TIME__; return;}' >> $@
