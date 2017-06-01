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

#Module
set(MODULE PWGmuondep)
add_definitions(-D_MODULE_="${MODULE}")

add_subdirectory(aaf)
message(STATUS "PWG muondep aaf enabled")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWG/muondep)

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliPhysics_SOURCE_DIR}/PWG/muon
  )

# Sources - alphabetical order
set(SRCS
  AliAnalysisMuMu.cxx
  AliAnalysisMuMuConfig.cxx
  AliAnalysisMuMuFnorm.cxx
  AliAnalysisMuMuGraphUtil.cxx
  AliAnalysisMuMuJpsiResult.cxx
  AliAnalysisMuMuResult.cxx
  AliAnalysisMuMuSpectra.cxx
  AliAnalysisMuMuSpectraProcessor.cxx
  AliAnalysisMuMuSpectraProcessorPP.cxx
  AliAnalysisMuMuSpectraProcessorPbP.cxx
  AliAnalysisMuMuSpectraProcessorPbPb.cxx
  AliAnalysisRunList.cxx
  AliAnalysisTaskESDMCLabelAddition.cxx 
  AliAnalysisTaskMuonCDBConnect.cxx
  AliAnalysisTaskMuonRefit.cxx
  AliAnalysisTaskMuonRefitVtx.cxx
  AliAnalysisTaskPileup.cxx
  AliAnalysisTriggerScalers.cxx
  AliMuonAccEffSubmitter.cxx
  AliMuonCompactCluster.cxx
  AliMuonCompactEvent.cxx
  AliMuonCompactManuStatus.cxx
  AliMuonCompactMapping.cxx
  AliMuonCompactQuickAccEff.cxx
  AliMuonCompactQuickAccEffChecker.cxx
  AliMuonCompactTrack.cxx
  AliMuonCompactTreeMaker.cxx
  AliMuonGridSubmitter.cxx
  AliMuonQAMergeSubmitter.cxx
  AliMuonOCDBSnapshotGenerator.cxx
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
set(LIBDEPS ANALYSISalice CDB MUONmapping MUONevaluation MUONrec PWGmuon STEERBase STEER)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Generate a PARfile target for this library
add_target_parfile(${MODULE} "${SRCS}" "${HDRS}" "${MODULE}LinkDef.h" "${LIBDEPS}")

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
install(DIRECTORY . DESTINATION PWG/muondep FILES_MATCHING PATTERN "*.C")
install(DIRECTORY . DESTINATION PWG/muondep FILES_MATCHING PATTERN "*.sh")
install(DIRECTORY . DESTINATION PWG/muondep FILES_MATCHING PATTERN "*.input")
