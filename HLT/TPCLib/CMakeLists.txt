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

set(LIBDEPS)
if(ALIROOT_ONLINE_MODE)
    message(WARNING "Online build for non-HLT: Disabling GPU")
else()
    set(ALIGPU_BUILD_TYPE "ALIROOT")
    add_subdirectory(${AliRoot_SOURCE_DIR}/GPU/GPUTracking AliGPUTracking)
    add_subdirectory(${AliRoot_SOURCE_DIR}/GPU/TPCFastTransformation TPCFastTransformation)
    include_directories(${AliRoot_SOURCE_DIR}/GPU/Common
                        ${AliRoot_SOURCE_DIR}/GPU/Utils
                        ${AliRoot_SOURCE_DIR}/GPU/GPUTracking/SliceTracker
                        ${AliRoot_SOURCE_DIR}/GPU/GPUTracking/Merger
                        ${AliRoot_SOURCE_DIR}/GPU/GPUTracking/Global
                        ${AliRoot_SOURCE_DIR}/GPU/GPUTracking/Base
                        ${AliRoot_SOURCE_DIR}/GPU/GPUTracking/DataCompression
                        ${AliRoot_SOURCE_DIR}/GPU/TPCFastTransformation)
    add_definitions(-DHAVE_ALIGPU)
    add_definitions(-DGPUCA_ALIROOT_LIB)
    set(LIBDEPS ${LIBDEPS} AliGPUTracking)
    set(LIBDEPS ${LIBDEPS} AliTPCFastTransformation)
endif()

# Module
set(MODULE AliHLTTPC)

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/HLT/TPCLib
                    ${AliRoot_SOURCE_DIR}/HLT/TPCLib/HWCFemulator/
                    ${AliRoot_SOURCE_DIR}/HLT/TPCLib/comp
                    ${AliRoot_SOURCE_DIR}/HLT/TPCLib/offline
                    ${AliRoot_SOURCE_DIR}/HLT/TPCLib/transform
                    ${AliRoot_SOURCE_DIR}/HLT/TPCLib/transform-new
)

# Additional include folders in alphabetical order except ROOT
include_directories(SYSTEM ${ROOT_INCLUDE_DIR})
include_directories(${AliRoot_SOURCE_DIR}/ANALYSIS/ANALYSISalice
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
    AliHLTTPCClusterAccessHLTOUT.cxx
    AliHLTTPCClusterDumpComponent.cxx
    AliHLTTPCClusterHistoComponent.cxx
    AliHLTTPCClusterTransformationComponent.cxx
    AliHLTTPCClusterTransformationPrepareComponent.cxx
    AliHLTTPCOfflinePreprocessorWrapperComponent.cxx
    AliHLTTPCClusterTransformation.cxx
    AliHLTTPCDataCheckerComponent.cxx
    AliHLTTPCDataPublisherComponent.cxx
    AliHLTTPCdEdxComponent.cxx
    AliHLTTPCFastdEdxComponent.cxx
    AliHLTTPCdEdxMonitoringComponent.cxx
    AliHLTTPCDefinitions.cxx
    AliHLTTPCDigitDumpComponent.cxx
    AliHLTTPCDigitReader32Bit.cxx
    AliHLTTPCDigitReader.cxx
    AliHLTTPCDigitReaderUnpacked.cxx
    AliHLTTPCFileHandler.cxx
    AliHLTTPCHistogramHandlerComponent.cxx
    AliHLTTPCHWCFData.cxx
    AliHLTTPCHWCFSpacePointContainer.cxx
    AliHLTTPCHWClusterDecoderComponent.cxx
    AliHLTTPCLog.cxx
    AliHLTTPCMapping.cxx
    AliHLTTPCMemHandler.cxx
    AliHLTTPCRawSpacePointContainer.cxx
    AliHLTTPCSpacePointContainer.cxx
    AliHLTTPCTrackGeometry.cxx
    AliHLTTPCTrackMCMarkerComponent.cxx
    AliHLTTPCGeometry.cxx
    AliHLTTPCPolynomFit.cxx
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
    offline/AliHLTTPCDigitPublisherComponent.cxx
    offline/AliHLTTPCOfflineCluster.cxx
    transform/AliHLTTPCFastTransform.cxx
    transform/AliHLTTPCSpline2D3D.cxx
    transform/AliHLTTPCSpline2D3DObject.cxx
    transform/AliHLTTPCFastTransformObject.cxx
    transform-new/AliFlexibleSpline1D.cxx
    transform-new/AliFlexibleSpline2D3D.cxx
    test/AliHLTTPCRawClusterDumpComponent.cxx
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
    transform/AliHLTTPCReverseTransformInfoV1.h
)

#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O0 -ggdb")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary_flat("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Generate the ROOT map
# Dependecies
set(LIBDEPS ${LIBDEPS} STEERBase STEER ESD CDB RAWDatabase RAWDatarec ANALYSISalice HLTbase AliHLTUtil AliHLTRCU TPCbase TPCrec TPCsim TPCcalib Graf Gpad Recorder)
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
