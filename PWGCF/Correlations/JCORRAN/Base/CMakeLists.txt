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
set(MODULE PWGCFCorrelationsJCORRAN)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWGCF/Correlations/JCORRAN/Base)

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliPhysics_SOURCE_DIR}/CORRFW
                    ${AliPhysics_SOURCE_DIR}/PWG/Tools
                    ${AliPhysics_SOURCE_DIR}/PWG/EMCAL/EMCALbase
                    ${AliPhysics_SOURCE_DIR}/OADB
                    ${AliPhysics_SOURCE_DIR}/OADB/COMMON/MULTIPLICITY
                    ${AliPhysics_SOURCE_DIR}/PWGUD/base
  )

# Sources - alphabetical order
set(SRCS
  AliJAcceptanceCorrection.cxx
  AliJAcceptanceFunctions.cxx
  AliJCORRAN.cxx
  AliJTrackCounter.cxx
  AliJCORRANTask.cxx
  AliJEfficiencyTask.cxx
  AliJDiHadronCorrTask.cxx
  AliJFilter.cxx
  AliJEfficiencyScanner.cxx
  AliJBaseEventHeader.cxx
  AliJBaseTrack.cxx
  AliJEventHeader.cxx
  AliJJet.cxx
  AliJMCTrack.cxx
  AliJPhoton.cxx
  AliJPiZero.cxx
  AliJRunHeader.cxx
  AliJTrack.cxx
  AliJHistManager.cxx
  AliJCorrelations.cxx
  AliJDataManager.cxx
  AliJRunTable.cxx
  AliJHistos.cxx
  AliJEventPool.cxx
  AliJEfficiency.cxx
  AliJTrackCut.cxx
  AliJBaseCard.cxx
  AliJCard.cxx
  AliJFFlucTask.cxx
  AliJFFlucAnalysis.cxx
  AliAnalysisAnaTwoMultiCorrelations.cxx
  AliJHSCTask.cxx
  AliJXtTask.cxx
  AliJXtAnalysis.cxx
  AliJHistogramInterface.cxx
  AliJCorrelationInterface.cxx
  AliJESE.cxx
  jtAnalysis/AliJDiHadronJtTask.cxx
  jtAnalysis/AliJJtAnalysis.cxx
  jtAnalysis/AliJJtHistograms.cxx
  jtAnalysis/AliJJtCorrelations.cxx
  iaaAnalysis/AliJDiHadronIaaTask.cxx
  iaaAnalysis/AliJIaaAnalysis.cxx
  iaaAnalysis/AliJIaaCorrelations.cxx
  iaaAnalysis/AliJIaaHistograms.cxx
  AliJPartLifetime.cxx
  AliJCatalystTask.cxx
  AliJCorrectionMapTask.cxx
  jtAnalysis/AliJJtTask.cxx
  jtAnalysis/AliJJtAna.cxx
  iaaAnalysis/AliJCIaaTask.cxx
  iaaAnalysis/AliJIaaAna.cxx
  )

# Headers from sources
  string(REPLACE ".cxx" ".h" HDRS "${SRCS}")
  set(HDRS ${HDRS} AliJConst.h) #add the remaining headers

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
	generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Add a shared library
add_library_tested(${MODULE} SHARED  ${SRCS} G__${MODULE}.cxx)

# Generate the ROOT map
# Dependecies
set(LIBDEPS ANALYSISalice CORRFW EMCALUtils OADB PHOSUtils)
	generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

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
