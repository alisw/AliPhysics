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
# **************************************************************************/


#except for true debug mode, disable AliDebug (and the AliCodeTimers)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS} -DLOG_NO_DEBUG ")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -DLOG_NO_DEBUG ")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} ")

#Module
set(MODULE PWGmuon)
add_definitions(-D_MODULE_="${MODULE}")

add_subdirectory(aaf)
message(STATUS "PWG muon aaf enabled")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWG/muon)

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliPhysics_SOURCE_DIR}/CORRFW
                    ${AliPhysics_SOURCE_DIR}/EVENTMIX
                    ${AliPhysics_SOURCE_DIR}/OADB
                    ${AliPhysics_SOURCE_DIR}/OADB/COMMON/MULTIPLICITY
  )

# Sources - alphabetical order
set(SRCS
  AliAnalysisCountTriggers.cxx
  AliAnalysisMuMuBase.cxx
  AliAnalysisMuMuBinning.cxx
  AliAnalysisMuMuCutCombination.cxx
  AliAnalysisMuMuCutElement.cxx
  AliAnalysisMuMuCutRegistry.cxx
  AliAnalysisMuMuEventCutter.cxx
  AliAnalysisMuMuGlobal.cxx
  AliAnalysisMuMuMCGene.cxx
  AliAnalysisMuMuMinv.cxx
  AliAnalysisMuMuNch.cxx
  AliAnalysisMuMuSingle.cxx
  AliAnalysisMuMuTriggerResponse.cxx
  AliAnalysisMuonUtility.cxx
  AliAnalysisNonMuonTrackCuts.cxx
  AliAnalysisNonPrimaryVertices.cxx
  AliAnalysisTaskAOD2MuonAOD.cxx
  AliAnalysisTaskAOD2MuonAOD.cxx
  AliAnalysisTaskCreateMixedDimuons.cxx
  AliAnalysisTaskDimuonCFContainerBuilder.cxx
  AliAnalysisTaskESDMuonFilter.cxx
  AliAnalysisTaskFlowSingleMu.cxx
  AliAnalysisTaskFromStandardToMuonAOD.cxx
  AliAnalysisTaskGenTuner.cxx
  AliAnalysisTaskLinkToMC.cxx
  AliAnalysisTaskLUT.cxx
  AliAnalysisTaskMTRResponse.cxx
  AliAnalysisTaskMuMu.cxx
  AliAnalysisTaskMuMu.cxx
  AliAnalysisTaskMuonAODCreation.cxx
  AliAnalysisTaskMuonAODfromGeneral.cxx
  AliAnalysisTaskMuonCollisionMultiplicity.cxx
  AliAnalysisTaskMuonDistributions.cxx
  AliAnalysisTaskMuonTreeBuilder.cxx
  AliAnalysisTaskNorm.cxx
  AliAnalysisTaskOnTheFlyAliAODDimuon.cxx
  AliAnalysisTaskSEMuonsHF.cxx
  AliAnalysisTaskSingleMu.cxx
  AliAnalysisTaskSingleMuESD.cxx
  AliAnalysisTaskWeightMTRResponse.cxx
  AliAODEventInfo.cxx
  AliAODMuonReplicator.cxx
  AliCFMuonResTask1.cxx
  AliCFMuonResUpsilon.cxx
  AliCFMuonSingleTask1.cxx
  AliDimuInfoStoreMC.cxx
  AliDimuInfoStoreRD.cxx
  AliESDMuonTrackCuts.cxx
  AliEventPoolMuon.cxx
  AliHistogramCollection.cxx
  AliMTRParameterizedResponse.cxx
  AliMuonAnalysisOutput.cxx
  AliMuonEventCuts.cxx
  AliMuonInfoStoreMC.cxx
  AliMuonInfoStoreRD.cxx
  AliMuonPairCuts.cxx
  AliMuonsHFHeader.cxx
  AliMuonTrackCuts.cxx
  AliMuonTrackSmearing.cxx
  AliMuonTriggerCombo.cxx
  AliOADBMuonTrackCutsParam.cxx
  AliTaskMuonTrackSmearing.cxx
  AliUtilityDimuonSource.cxx
  AliUtilityMuonAncestor.cxx
  AliVAnalysisMuon.cxx
  AliAnalysisTaskLMREventFilter.cxx
  AliLMREvent.cxx
  AliLMRMuon.cxx
)

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Add a shared library
add_library_tested(${MODULE} SHARED  ${SRCS} G__${MODULE}.cxx)

# Generate the ROOT map
# Dependecies
set(LIBDEPS ANALYSISalice AOD CORRFW OADB STEERBase EventMixing)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Generate a PARfile target for this library
add_target_parfile(${MODULE} "${SRCS}" "${HDRS}" "${MODULE}LinkDef.h" "${LIBDEPS}" "../OADB/COMMON/MULTIPLICITY")

# Linking the library
target_link_libraries(${MODULE} ${LIBDEPS})

# Public include folders that will be propagated to the dependecies
target_include_directories(${MODULE} PUBLIC ${incdirs})

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Installation
install(TARGETS ${MODULE}
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
install(FILES ${HDRS} DESTINATION include)

# Installing the macros
install(DIRECTORY . DESTINATION PWG/muon FILES_MATCHING PATTERN "*.C")
