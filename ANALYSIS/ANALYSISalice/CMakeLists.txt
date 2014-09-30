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
set(MODULE ANALYSISalice)

# Module include folders in alphabetical order except ROOT
include_directories(${CMAKE_SOURCE_DIR}/ANALYSIS/${MODULE})

# Additional include directories
include_directories(${ROOT_INCLUDE_DIRS}
                    ${CMAKE_SOURCE_DIR}/ANALYSIS/ANALYSIS
                    ${CMAKE_SOURCE_DIR}/OADB
                    ${CMAKE_SOURCE_DIR}/STEER/AOD
                    ${CMAKE_SOURCE_DIR}/STEER/ESD
                    ${CMAKE_SOURCE_DIR}/STEER/STEERBase
                   )

# Sources in alphabetical order 
set(SRCS 
    AliAnalysisCuts.cxx
    AliAnalysisFilter.cxx
    AliAnalysisTaskBadChunkID.cxx
    AliAnalysisTaskBaseLine.cxx
    AliAnalysisTaskMCParticleFilter.cxx
    AliAnalysisTaskME.cxx
    AliAnalysisTaskPIDCombined.cxx
    AliAnalysisTaskPIDqa.cxx
    AliAnalysisTaskPIDResponse.cxx
    AliAnalysisTaskSE.cxx
    AliAnalysisTaskStat.cxx
    AliAnalysisTaskTagCreator.cxx
    AliAnalysisUtils.cxx
    AliAODv0KineCuts.cxx
    AliBackgroundSelection.cxx
    AliCentralitySelectionTask.cxx
    AliCollisionNormalization.cxx
    AliCollisionNormalizationTask.cxx
    AliEPSelectionTask.cxx
    AliESDpidCuts.cxx
    AliESDtrackCuts.cxx
    AliESDv0Cuts.cxx
    AliESDv0KineCuts.cxx
    AliEventPoolLoop.cxx
    AliEventPoolManager.cxx
    AliEventPoolOTF.cxx
    AliEventPoolSparse.cxx
    AliHEPDataParser.cxx
    AliKineTrackCuts.cxx
    AliMultiEventInputHandler.cxx
    AliMultiInputEventHandler.cxx
    AliPhysicsSelection.cxx
    AliPhysicsSelectionTask.cxx
    AliPIDResponseInputHandler.cxx
    AliPPVsMultUtils.cxx
    AliTagAnalysis.cxx
    AliTriggerAnalysis.cxx
    AliUnfolding.cxx
    AliVZEROEPSelectionTask.cxx
    AliXMLCollection.cxx
  )

# If ROOT was build with no AliEn support we disable the AliROOT AliEn support
if(ROOT_HASALIEN STREQUAL "no")
    message(STATUS "ROOT was build with no alien support. Disabling AliROOT AliEn support")
else()
    set(SRCS ${SRCS} AliAnalysisAlien.cxx)
endif(ROOT_HASALIEN STREQUAL "no")

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

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
