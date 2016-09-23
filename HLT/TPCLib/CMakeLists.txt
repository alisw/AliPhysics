# **************************************************************************
# * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
# *                                                                        *
# * Author: The ALICE Off-line Project.                                    *
# * Contributors are mentioned in the code where appropriate.              *
# *                                                                        *
# * Permission to use, copy, modify and distribute this software and its   *
# * documentation strictly for non-commercial purposes is hereby granted   *
# * without fee, provided that the above copyright notice appears in all   *
# * copies and that both the copyright notice and this permission notice   *
# * appear in the supporting documentation. The authors make no claims     *
# * about the suitability of this software for any purpose. It is          *
# * provided "as is" without express or implied warranty.                  *
# **************************************************************************

# Other TPC libraries
# Submodules have to be included before the main library otherwise
# they will inherit parent folder properties

if(ROOT_HASOPENGL)
    add_subdirectory(EVE)
endif()

if(HLT_CUDA)
    add_subdirectory(tracking-ca/cagpu/cuda)
endif(HLT_CUDA)

if(HLT_OPENCL)
    add_subdirectory(tracking-ca/cagpu/opencl)
endif()

# Module
set(MODULE AliHLTTPC)

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/HLT/TPCLib
                    ${AliRoot_SOURCE_DIR}/HLT/TPCLib/HWCFemulator/
                    ${AliRoot_SOURCE_DIR}/HLT/TPCLib/comp
                    ${AliRoot_SOURCE_DIR}/HLT/TPCLib/merger-ca
                    ${AliRoot_SOURCE_DIR}/HLT/TPCLib/offline
                    ${AliRoot_SOURCE_DIR}/HLT/TPCLib/tracking-ca
                    ${AliRoot_SOURCE_DIR}/HLT/TPCLib/transform
)

# Additional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIR}
                    ${AliRoot_SOURCE_DIR}/ANALYSIS/ANALYSISalice
                    ${AliRoot_SOURCE_DIR}/HLT/BASE
                    ${AliRoot_SOURCE_DIR}/HLT/BASE/util
                    ${AliRoot_SOURCE_DIR}/HLT/RCU
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                    ${AliRoot_SOURCE_DIR}/TPC/TPCbase
                    ${AliRoot_SOURCE_DIR}/TPC/TPCcalib
                    ${AliRoot_SOURCE_DIR}/TPC/TPCrec
                    ${AliRoot_SOURCE_DIR}/TPC/TPCsim
                   )

# Enable Vc
ALICE_UseVc()

# Sources in alphabetical order
set(SRCS
    AliHLTTPCAgent.cxx
    AliHLTTPCCFComparisonComponent.cxx
    AliHLTTPCClusterAccessHLTOUT.cxx    
    AliHLTTPCClusterDumpComponent.cxx
    AliHLTTPCClusterFinderComponent.cxx
    AliHLTTPCClusterFinder.cxx
    AliHLTTPCClusterHistoComponent.cxx
    AliHLTTPCClusters.cxx
    AliHLTTPCClusterTransformationComponent.cxx
    AliHLTTPCClusterTransformationPrepareComponent.cxx
    AliHLTTPCOfflinePreprocessorWrapperComponent.cxx
    AliHLTTPCClusterTransformation.cxx
    AliHLTTPCDataCheckerComponent.cxx
    AliHLTTPCDataPublisherComponent.cxx
    AliHLTTPCdEdxComponent.cxx
    AliHLTTPCdEdxMonitoringComponent.cxx
    AliHLTTPCDefinitions.cxx
    AliHLTTPCDigitDumpComponent.cxx
    AliHLTTPCDigitReader32Bit.cxx
    AliHLTTPCDigitReader.cxx
    AliHLTTPCDigitReaderUnpacked.cxx
    AliHLTTPCFileHandler.cxx
    AliHLTTPCHistogramHandlerComponent.cxx
    AliHLTTPCHWCFData.cxx
    AliHLTTPCHWCFDataReverterComponent.cxx
    AliHLTTPCHWCFSpacePointContainer.cxx
    AliHLTTPCHWClusterDecoderComponent.cxx
    AliHLTTPCLog.cxx
    AliHLTTPCMapping.cxx
    AliHLTTPCMemHandler.cxx
    AliHLTTPCPad.cxx
    AliHLTTPCRawSpacePointContainer.cxx
    AliHLTTPCSpacePointContainer.cxx
    AliHLTTPCTrackGeometry.cxx
    AliHLTTPCTrackHistoComponent.cxx
    AliHLTTPCTrackMCMarkerComponent.cxx
    AliHLTTPCGeometry.cxx
    comp/AliHLTTPCClusterStatComponent.cxx
    comp/AliHLTTPCDataCompressionComponent.cxx
    comp/AliHLTTPCDataCompressionDecoder.cxx
    comp/AliHLTTPCDataCompressionFilterComponent.cxx
    comp/AliHLTTPCDataCompressionIDMap.cxx
    comp/AliHLTTPCDataCompressionMonitorComponent.cxx
    comp/AliHLTTPCDataCompressionUnpackerComponent.cxx
    HWCFemulator/AliHLTTPCHWCFConsistencyControlComponent.cxx
    HWCFemulator/AliHLTTPCHWCFDivisionUnit.cxx
    HWCFemulator/AliHLTTPCHWCFEmulatorComponent.cxx
    HWCFemulator/AliHLTTPCHWCFEmulator.cxx
    HWCFemulator/AliHLTTPCHWCFExtractorUnit.cxx
    HWCFemulator/AliHLTTPCHWCFMergerUnit.cxx
    HWCFemulator/AliHLTTPCHWCFPeakFinderUnit.cxx
    HWCFemulator/AliHLTTPCHWCFProcessorUnit.cxx
    HWCFemulator/AliHLTTPCHWCFSupport.cxx
    HWCFemulator/AliHLTTPCHWClusterMerger.cxx
    HWCFemulator/AliHLTTPCHWClusterMergerV1.cxx
    merger-ca/AliHLTTPCGMMerger.cxx
    merger-ca/AliHLTTPCGMSliceTrack.cxx
    merger-ca/AliHLTTPCGMTrackParam.cxx
    offline/AliHLTTPCDigitPublisherComponent.cxx
    offline/AliHLTTPCOfflineCluster.cxx
    tracking-ca/AliHLTTPCCABaseTrackParam.cxx
    tracking-ca/AliHLTTPCCAClusterData.cxx
    tracking-ca/AliHLTTPCCADisplay.cxx
    tracking-ca/AliHLTTPCCAGlobalMergerComponent.cxx
    tracking-ca/AliHLTTPCCAGPUTracker.cxx
    tracking-ca/AliHLTTPCCAGrid.cxx
    tracking-ca/AliHLTTPCCAHitArea.cxx
    tracking-ca/AliHLTTPCCAInputDataCompressorComponent.cxx
    tracking-ca/AliHLTTPCCAMCPoint.cxx
    tracking-ca/AliHLTTPCCAMCTrack.cxx
    tracking-ca/AliHLTTPCCAMerger.cxx
    tracking-ca/AliHLTTPCCANeighboursCleaner.cxx
    tracking-ca/AliHLTTPCCANeighboursFinder.cxx
    tracking-ca/AliHLTTPCCAParam.cxx
    tracking-ca/AliHLTTPCCAPerformance.cxx
    tracking-ca/AliHLTTPCCARow.cxx
    tracking-ca/AliHLTTPCCASliceData.cxx
    tracking-ca/AliHLTTPCCASliceOutput.cxx
    tracking-ca/AliHLTTPCCAStandaloneFramework.cxx
    tracking-ca/AliHLTTPCCAStartHitsFinder.cxx
    tracking-ca/AliHLTTPCCATrackConvertor.cxx
    tracking-ca/AliHLTTPCCATrackerComponent.cxx
    tracking-ca/AliHLTTPCCATrackerFramework.cxx
    tracking-ca/AliHLTTPCCATracker.cxx
    tracking-ca/AliHLTTPCCATrackerOutputConverter.cxx
    tracking-ca/AliHLTTPCCATrack.cxx
    tracking-ca/AliHLTTPCCATrackletConstructor.cxx
    tracking-ca/AliHLTTPCCATracklet.cxx
    tracking-ca/AliHLTTPCCATrackletSelector.cxx
    tracking-ca/AliHLTTPCCATrackParam.cxx
    tracking-ca/AliTPCtrackerCA.cxx
    transform/AliHLTTPCFastTransform.cxx
    transform/AliHLTTPCSpline2D3D.cxx
    transform/AliHLTTPCSpline2D3DObject.cxx
    transform/AliHLTTPCFastTransformObject.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")
set(HDRS ${HDRS}
    AliHLTTPCLog.h
    AliHLTTPCLogging.h
    AliHLTTPCSpacePointData.h
    AliHLTTPCRawCluster.h
    AliHLTTPCRawClustersDescriptor.h
    AliHLTTPCDigitData.h
    AliHLTTPCClusterDataFormat.h
    AliHLTTPCClusterMCData.h
    AliHLTTPCClusterFlagsData.h
    HWCFemulator/AliHLTTPCHWCFDataTypes.h
    tracking-ca/AliHLTTPCCADef.h
    tracking-ca/AliHLTTPCCAMath.h
)

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Generate the ROOT map
# Dependecies
set(LIBDEPS STEERBase STEER ESD CDB RAWDatabase RAWDatarec ANALYSISalice HLTbase AliHLTUtil AliHLTRCU TPCbase TPCrec TPCsim TPCcalib Eve Graf Gpad Recorder)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# NOTE: Vc should not stay in the rootmap!
set(LIBDEPS ${LIBDEPS} Vc)

# Add a library to the project using the specified source files
add_library_tested(${MODULE} SHARED ${SRCS} G__${MODULE}.cxx)
target_link_libraries(${MODULE} ${LIBDEPS})

# Additional compilation flags
set_target_properties(${MODULE} PROPERTIES COMPILE_FLAGS "")

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Installation
install(TARGETS ${MODULE}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

install(FILES ${HDRS} DESTINATION include)
