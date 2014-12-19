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
set (MODULE PWGEMCAL)

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/PWG/EMCAL)

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliRoot_SOURCE_DIR}/ANALYSIS/ANALYSIS
                    ${AliRoot_SOURCE_DIR}/ANALYSIS/ANALYSISalice
                    ${AliRoot_SOURCE_DIR}/ANALYSIS/Tender
                    ${AliRoot_SOURCE_DIR}/ANALYSIS/TenderSupplies
                    ${AliRoot_SOURCE_DIR}/EMCAL/EMCALbase
                    ${AliRoot_SOURCE_DIR}/EMCAL/EMCALrec
                    ${AliRoot_SOURCE_DIR}/EMCAL/EMCALUtils
                    ${AliRoot_SOURCE_DIR}/PWG/Tools
                    ${AliRoot_SOURCE_DIR}/STEER/AOD
                    ${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
  )

# Sources - alphabetical order
set(SRCS
  AliAnalysisTaskEMCALClusterizeFast.cxx
  AliAnalysisTaskEmcal.cxx
  AliAnalysisTaskEmcalSample.cxx
  AliClusterContainer.cxx
  AliEMCALClusterParams.cxx
  AliEmcalAodTrackFilterTask.cxx
  AliEmcalClusTrackMatcherTask.cxx
  AliEmcalClusterMaker.cxx
  AliEmcalCompatTask.cxx
  AliEmcalContainer.cxx
  AliEmcalDebugTask.cxx
  AliEmcalEsdTrackFilterTask.cxx
  AliEmcalMCTrackSelector.cxx
  AliEmcalParticle.cxx
  AliEmcalParticleMaker.cxx
  AliEmcalPatchFromCellMaker.cxx
  AliEmcalPhysicsSelection.cxx
  AliEmcalPhysicsSelectionTask.cxx
  AliEmcalPicoTrackMaker.cxx
  AliEmcalSetupTask.cxx
  AliEmcalTenderTask.cxx
  AliEmcalTrackingQATask.cxx
  AliEmcalTrackPropagatorTask.cxx
  AliEmcalTriggerBitConfig.cxx
  AliEmcalTriggerMaker.cxx
  AliEmcalTriggerPatchInfo.cxx
  AliEmcalTriggerSetupInfo.cxx
  AliEsdSkimTask.cxx
  AliEsdTrackExt.cxx
  AliHadCorrTask.cxx
  AliParticleContainer.cxx
  AliPicoTrack.cxx
  AliEMCALConfiguration.cxx
  AliEMCALConfigHandler.cxx
  AliEMCALConfigurationMatcher.cxx
  AliEmcalTriggerDecision.cxx
  AliEmcalTriggerDecisionContainer.cxx
  AliEmcalTriggerSelectionCuts.cxx
  AliEmcalTriggerSelection.cxx
  AliAnalysisTaskEmcalTriggerSelection.cxx
  )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Add a shared library
add_library(${MODULE} SHARED  ${SRCS} G__${MODULE}.cxx)

# Generate the ROOT map
# Dependecies
set(LIBDEPS ANALYSIS ANALYSISalice AOD CDB EMCALrec EMCALUtils ESD PWGTools STEER STEERBase Tender)
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

# Installing the macros
install (DIRECTORY macros DESTINATION PWG/EMCAL)
install (DIRECTORY data DESTINATION PWG/EMCAL)
