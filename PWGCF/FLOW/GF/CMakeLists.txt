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
set(MODULE PWGCFFLOWGF)
#add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(
    ${AliPhysics_SOURCE_DIR}/PWGCF/FLOW/GF
)

# Additional includes - alphabetical order except ROOT
include_directories(
    ${ROOT_INCLUDE_DIRS}
    ${AliPhysics_SOURCE_DIR}/CORRFW
    ${AliPhysics_SOURCE_DIR}/OADB
	  ${AliPhysics_SOURCE_DIR}/OADB/COMMON/MULTIPLICITY
    ${AliPhysics_SOURCE_DIR}/PWG/EMCAL/EMCALbase
    ${AliPhysics_SOURCE_DIR}/PWG/Tools
    ${AliPhysics_SOURCE_DIR}/PWGHF/vertexingHF
    ${AliPhysics_SOURCE_DIR}/PWGCF/FLOW/GF
    ${AliPhysics_SOURCE_DIR}/PWGLF/FORWARD/analysis2
		${AliPhysics_SOURCE_DIR}/EVENTMIX
)

# Sources - alphabetical order
set(SRCS
  AliAnalysisTaskChargedFlow.cxx
  AliAnalysisTaskCorrForFlow.cxx
  AliAnalysisTaskCorrForFlowFMD.cxx
  AliAnalysisTaskCorrForFlowMaster.cxx
  AliAnalysisTaskCorrForFlowEta.cxx
  AliAnalysisTaskDeform.cxx
  AliAnalysisTaskFlowExample.cxx
  AliAnalysisTaskGFWPIDFlow.cxx
  AliAnalysisTaskPtCorr.cxx
  AliAnalysisTaskUniFlow.cxx
  AliAnalysisTaskUniFlowWithSphericity.cxx
  AliAnalysisTaskUniFlowMultiStrange.cxx
  AliDecorrFlowCorrTask.cxx
  AliGFW.cxx
  AliGFWCumulant.cxx
  AliGFWCuts.cxx
  AliGFWFlowContainer.cxx
  AliGFWWeights.cxx
  AliLWUtils.cxx
  AliPartSimpleForCorr.cxx
  AliProfileSubset.cxx
  AliProfileBS.cxx
  AliPtContainer.cxx
  AliCkContainer.cxx
  AliUniFlowCorrTask.cxx
  AliAnalysisDecorrTask.cxx
  AliAnalysisTaskESEFlow.cxx
  AliAnalysisTaskFlowSquareBracket.cxx
  AliAnalysisTaskCMWESETrkSyst.cxx
  AliAnalysisTaskPtFlowCorrelation.cxx
  AliAnalysisTaskMomentumFlowCorrelation.cxx
  AliPtSubEventContainer.cxx
  AliAnalysisTaskDensity.cxx
  )

if(${ROOT_VERSION} GREATER_EQUAL 6.0)
  list(APPEND SRCS AliAnalysisTaskCreateNUA.cxx)
  list(APPEND SRCS AliAnalysisTaskCreatePeriodNUA.cxx)
  list(APPEND SRCS AliAnalysisTaskCreateNUE.cxx)
  list(APPEND SRCS AliAnalysisTaskEffFDExample.cxx)
  list(APPEND SRCS AliAnalysisTaskFlowPPTask.cxx)
  list(APPEND SRCS AliAnalysisTaskFlowOnTheFly.cxx)
  list(APPEND SRCS AliAnalysisTaskGFWFlow.cxx)
  list(APPEND SRCS AliAnalysisTaskJetQ.cxx)
  list(APPEND SRCS AliAnalysisTaskLWTree.cxx)
  list(APPEND SRCS AliAnalysisTaskMeanPtV2Corr.cxx)
  list(APPEND SRCS AliAnalysisTaskNonlinearFlow.cxx)
  list(APPEND SRCS AliAnalysisTaskXDeptFlow.cxx)
  list(APPEND SRCS AliAnalysisTaskCorrForNonlinearFlow.cxx)
  list(APPEND SRCS AliAnalysisTaskAcceptanceCorr.cxx)
  list(APPEND SRCS AliEffFDContainer.cxx)
  list(APPEND SRCS AliGFWFilter.cxx)
  list(APPEND SRCS AliGFWFilterTask.cxx)
  list(APPEND SRCS CorrelationCalculator.cxx)
  list(APPEND SRCS AliGFWMCuts.cxx)
  list(APPEND SRCS AliGFWNFCuts.cxx)
  list(APPEND SRCS AliGFWXXCuts.cxx)
  list(APPEND SRCS AliGFWPbpass23Cuts.cxx)
endif()

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")
# Also add the header with GFW flags
list(APPEND HDRS GFWFlags.h)

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Add a shared library
add_library_tested(${MODULE} SHARED  ${SRCS} G__${MODULE}.cxx)

# Generate the ROOT map
# Dependecies
set(LIBDEPS ANALYSISalice OADB EMCALbase PWGTools)
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
