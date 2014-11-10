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

# HLTbase + subfolders (end of the cmake file)
# Subfolders hve to be includes at the begining otherwise
# they will inherit properties from the parent
add_subdirectory(HOMER)
add_subdirectory(interface)
add_subdirectory(util)

# Module
set(MODULE HLTbase)

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/HLT/BASE
                    ${AliRoot_SOURCE_DIR}/HLT/BASE/HOMER
                   )

# Additional include folders in alphabetical order except ROOT
include_directories(
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )

# Sources in alphabetical order
set(SRCS
    AliHLTBlockDataCollection.cxx
    AliHLTCalibrationProcessor.cxx
    AliHLT_C_Component_WrapperInterface.cxx
    AliHLT_C_Component_WrapperInterface.cxx
    AliHLTComponentBenchmark.cxx
    AliHLTComponentConfiguration.cxx
    AliHLTComponent.cxx
    AliHLTComponentHandler.cxx
    AliHLTConfiguration.cxx
    AliHLTConfigurationHandler.cxx
    AliHLTConsumerDescriptor.cxx
    AliHLTControlTask.cxx
    AliHLTCTPData.cxx
    AliHLTDAQ.cxx
    AliHLTDataBuffer.cxx
    AliHLTDataDeflater.cxx
    AliHLTDataDeflaterHuffman.cxx
    AliHLTDataDeflaterSimple.cxx
    AliHLTDataInflater.cxx
    AliHLTDataInflaterHuffman.cxx
    AliHLTDataInflaterSimple.cxx
    AliHLTDataSink.cxx
    AliHLTDataSource.cxx
    AliHLTDataTypes.cxx
    AliHLTDataTypes.cxx
    AliHLTDimServer.cxx
    AliHLTDomainEntry.cxx
    AliHLTDumpTask.cxx
    AliHLTErrorGuard.cxx
    AliHLTEsdManager.cxx
    AliHLTEventStatistics.cxx
    AliHLTGlobalTriggerDecision.cxx
    AliHLTHOMERBlockDesc.cxx
    AliHLTHOMERLibManager.cxx
    AliHLTHOMERManager.cxx
    AliHLTHOMERProxyHandler.cxx
    AliHLTHOMERSourceDesc.cxx
    AliHLTHuffman.cxx
    AliHLTLogging.cxx
    AliHLTMemoryFile.cxx
    AliHLTMessage.cxx
    AliHLTMisc.cxx
    AliHLTModuleAgent.cxx
    AliHLTModulePreprocessor.cxx
    AliHLTOfflineDataSource.cxx
    AliHLTOfflineInterface.cxx
    AliHLTOnlineConfiguration.cxx
    AliHLTOUT.cxx
    AliHLTOUTComponent.cxx
    AliHLTOUTHandlerChain.cxx
    AliHLTOUTHandlerDetectorDDL.cxx
    AliHLTOUTHandlerEquId.cxx
    AliHLTOUTHandlerEsdBranch.cxx
    AliHLTOUTHandler.cxx
    AliHLTOUTHandlerIgnore.cxx
    AliHLTOUTHomerBuffer.cxx
    AliHLTOUTTask.cxx
    AliHLTPluginBase.cxx
    AliHLTProcessor.cxx
    AliHLTReadoutList.cxx
    AliHLTReconstructorBase.cxx
    AliHLTRunStatistics.cxx
    AliHLTScalars.cxx
    AliHLTShuttleInterface.cxx
    AliHLTSpacePointContainer.cxx
    AliHLTSystem.cxx
    AliHLTTask.cxx
    AliHLTTrackGeometry.cxx
    AliHLTTriggerDecision.cxx
    AliHLTTriggerDomain.cxx
    AliHLTTriggerMenu.cxx
    AliHLTTriggerMenuItem.cxx
    AliHLTTriggerMenuSymbol.cxx
    AliHLTTTreeProcessor.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")
set(HDRS ${HDRS}
    AliHLTCommonCDBEntries.h
    AliHLTDefinitions.h
    AliHLTExternalTrackParam.h
    AliHLTLoggingVariadicFree.h
    AliHLTStdIncludes.h
   )

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Generate the ROOT map
# Dependecies
set(LIBDEPS STEERBase RAWDatabase AliHLTHOMER)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Add a library to the project using the specified source files
add_library(${MODULE} SHARED ${SRCS} G__${MODULE}.cxx)

# Additional compilation flags
set_target_properties(${MODULE} PROPERTIES COMPILE_FLAGS "-O -g")

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Installation
install(TARGETS ${MODULE}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

install(FILES ${HDRS} DESTINATION include)