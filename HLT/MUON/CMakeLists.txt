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

# Module
set(MODULE AliHLTMUON)

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/HLT/MUON
                    ${AliRoot_SOURCE_DIR}/HLT/MUON/OnlineAnalysis
                    ${AliRoot_SOURCE_DIR}/HLT/MUON/utils
                   )

# Additional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIR}
                    ${AliRoot_SOURCE_DIR}/HLT/BASE
                    ${AliRoot_SOURCE_DIR}/MUON/MUONbase
                    ${AliRoot_SOURCE_DIR}/MUON/MUONcalib
                    ${AliRoot_SOURCE_DIR}/MUON/MUONcore
                    ${AliRoot_SOURCE_DIR}/MUON//MUONevaluation
                    ${AliRoot_SOURCE_DIR}/MUON//MUONgeometry
                    ${AliRoot_SOURCE_DIR}/MUON/MUONmapping
                    ${AliRoot_SOURCE_DIR}/MUON/MUONraw
                    ${AliRoot_SOURCE_DIR}/MUON/MUONrec
                    ${AliRoot_SOURCE_DIR}/MUON/MUONsim
                    ${AliRoot_SOURCE_DIR}/MUON/MUONtrigger
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatasim
                    ${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )

# Sources in alphabetical order
set(SRCS
    AliHLTMUONChannelsBlockStruct.cxx
    AliHLTMUONClustersBlockStruct.cxx
    AliHLTMUONConstants.cxx
    AliHLTMUONDecision.cxx
    AliHLTMUONEvent.cxx
    AliHLTMUONMansoCandidatesBlockStruct.cxx
    AliHLTMUONMansoTrack.cxx
    AliHLTMUONMansoTracksBlockStruct.cxx
    AliHLTMUONPairsDecisionBlockStruct.cxx
    AliHLTMUONProcessor.cxx
    AliHLTMUONRecHit.cxx
    AliHLTMUONRecHitsBlockStruct.cxx
    AliHLTMUONSinglesDecisionBlockStruct.cxx
    AliHLTMUONTrack.cxx
    AliHLTMUONTracksBlockStruct.cxx
    AliHLTMUONTriggerRecord.cxx
    AliHLTMUONTriggerRecordsBlockStruct.cxx
    AliHLTMUONTrigRecsDebugBlockStruct.cxx
    AliHLTMUONUtils.cxx
    OfflineInterface/AliHLTMUONAgent.cxx
    OfflineInterface/AliHLTMUONDigitPublisherComponent.cxx
    OfflineInterface/AliHLTMUONESDMaker.cxx
    OfflineInterface/AliHLTMUONRecHitsSource.cxx
    OfflineInterface/AliHLTMUONRootifierComponent.cxx
    OfflineInterface/AliHLTMUONTriggerRecordsSource.cxx
    OnlineAnalysis/AliHLTMUONCalculations.cxx
    OnlineAnalysis/AliHLTMUONClusterFinderComponent.cxx
    OnlineAnalysis/AliHLTMUONDecisionComponent.cxx
    OnlineAnalysis/AliHLTMUONFullTrackerComponent.cxx
    OnlineAnalysis/AliHLTMUONFullTracker.cxx
    OnlineAnalysis/AliHLTMUONHitReconstructorComponent.cxx
    OnlineAnalysis/AliHLTMUONHitReconstructor.cxx
    OnlineAnalysis/AliHLTMUONMansoTrackerFSMComponent.cxx
    OnlineAnalysis/AliHLTMUONMansoTrackerFSM.cxx
    OnlineAnalysis/AliHLTMUONTriggerReconstructorComponent.cxx
    OnlineAnalysis/AliHLTMUONTriggerReconstructor.cxx
    utils/AliHLTMUONClusterHistoComponent.cxx
    utils/AliHLTMUONDataCheckerComponent.cxx
    utils/AliHLTMUONEmptyEventFilterComponent.cxx
    utils/AliHLTMUONRawDataHistoComponent.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")
set(HDRS ${HDRS}
    AliHLTMUONDataBlockReader.h
    AliHLTMUONDataBlockWriter.h
    AliHLTMUONDataTypes.h
    AliHLTMUONList.h
    AliHLTMUONCountedList.h
    OnlineAnalysis/AliHLTMUONMansoTrackerFSMCallback.h
  )

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Generate the ROOT map
# Dependecies
set(LIBDEPS STEERBase STEER CDB ESD RAWDatabase RAWDatarec RAWDatasim HLTbase MUONraw MUONcore MUONcalib MUONmapping MUONsim MUONevaluation MUONtrigger MUONgeometry MUONrec)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

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

# Compilation of some executables
add_subdirectory(utils)
      
