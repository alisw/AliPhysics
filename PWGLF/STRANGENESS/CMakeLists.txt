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
set(MODULE PWGLFSTRANGENESS)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWGLF/STRANGENESS)

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliPhysics_SOURCE_DIR}/EVENTMIX
                    ${AliPhysics_SOURCE_DIR}/CORRFW
                    ${AliPhysics_SOURCE_DIR}/OADB
                    ${AliPhysics_SOURCE_DIR}/OADB/COMMON/MULTIPLICITY
                    ${AliPhysics_SOURCE_DIR}/PWG/Tools
                    ${AliPhysics_SOURCE_DIR}/PWGCF/Correlations/Base
		                ${AliPhysics_SOURCE_DIR}/PWGLF/FORWARD/analysis2
                    ${AliPhysics_SOURCE_DIR}/PWGLF/STRANGENESS/Cascades/lightvertexers
                    ${AliPhysics_SOURCE_DIR}/PWGLF/STRANGENESS/Cascades/Run2
                    ${AliPhysics_SOURCE_DIR}/PWGUD/base
  )

# Sources - alphabetical order
set(SRCS
  Cascades/AliAnalysisTaskCheckCascade.cxx
  Cascades/AliAnalysisTaskCheckCascadePbPb.cxx
  Cascades/AliAnalysisTaskCheckPerformanceCascade.cxx
  Cascades/AliAnalysisTaskCheckPerformanceCascadePbPb.cxx
  LambdaK0PbPb/AliAnalysisCentralitySelector.cxx
  LambdaK0PbPb/AliAnalysisTaskCTauPbPb.cxx
  LambdaK0PbPb/AliAnalysisTaskCTauPbPbaod.cxx
  LambdaK0PbPb/AliAnalysisTaskLukeAOD.cxx
  LambdaK0PbPb/AliAnalysisTaskLukeV0.cxx
  LambdaK0PbPb/AliAnalysisTaskPerformanceStrange.cxx
  LambdaK0PbPb/AliAnalysisTaskV0ForRAA.cxx
  LambdaK0/AliAnalysisTaskExtractV0.cxx
  LambdaK0/AliAnalysisTaskExtractV0pPb.cxx
  LambdaK0/AliAnalysisTaskExtractV0AOD.cxx
  LambdaK0/AliAnalysisTaskExtractV0AODRun2.cxx
  LambdaK0/AliAnalysisTaskExtractPerformanceV0.cxx
  LambdaK0/AliAnalysisTaskExtractPerformanceV0pPb.cxx
  LambdaK0PbPb/AliAnalysisTaskLK0Spectra.cxx
  Lifetimes/AliAnalysisTaskStrangenessLifetimes.cxx
  Cascades/AliAnalysisTaskExtractCascade.cxx
  Cascades/AliAnalysisTaskExtractCascadePbPbRun2.cxx
  Cascades/AliAnalysisTaskExtractPerformanceCascade.cxx
  Cascades/AliAnalysisTaskStrangenessVsMultiplicity.cxx
  Cascades/AliAnalysisTaskStrangenessVsMultiplicityMC.cxx
  Cascades/Run2/AliAnalysisTaskStrangenessVsMultiplicityRun2.cxx
  Cascades/Run2/AliAnalysisTaskStrangenessVsMultiplicityEERun2.cxx
  Cascades/Run2/AliAnalysisTaskStrangenessVsMultiplicityAODRun2.cxx
  Cascades/Run2/AliAnalysisTaskStrangenessVsMultiplicityRsnLikeBgSub.cxx
  Cascades/Run2/AliAnalysisTaskStrangenessVsMultiplicityMCRun2.cxx
  Cascades/Run2/AliAnalysisTaskStrangenessVsMultiplicityEEMCRun2.cxx
  Cascades/Run2/AliAnalysisTaskStrangenessVsMultiplicityRun2pPb.cxx
  Cascades/Run2/AliAnalysisTaskStrangenessVsMultiplicityMCRun2pPb.cxx
  Cascades/AliAnalysisTaskPPVsMultCrossCheckMC.cxx
  Cascades/AliAnalysisTaskpANormalizationCheckMC.cxx
  Correlations/AliAnalysisTaskLambdaOverK0sJets.cxx
  Correlations/AliAnalysisTaskV0ChCorrelations.cxx
  Correlations/AliAnalysisTaskhStrCorr.cxx
  Cascades/AliAnalysisTaskCheckCascadepp.cxx
  Cascades/AliAnalysisTaskCheckPerformanceCascadepp.cxx
  Cascades/lightvertexers/AliLightV0vertexer.cxx
  Cascades/lightvertexers/AliLightCascadeVertexer.cxx
  Cascades/lightvertexers/AliCascadeVertexerUncheckedCharges.cxx
  Cascades/lightvertexers/AliV0vertexerUncheckedCharges.cxx
  Cascades/Run2/AliVWeakResult.cxx
  Cascades/Run2/AliV0Result.cxx
  Cascades/Run2/AliCascadeResult.cxx
  Cascades/Run2/AliStrangenessModule.cxx
  Cascades/Run2/AliAnalysisTaskWeakDecayVertexer.cxx
  Cascades/Run2/AliAnalysisTaskStrEffStudy.cxx
  Cascades/Run2/AliAnalysisTaskMCPredictions.cxx
  Cascades/Run2/AliAnalysisTaskStrangeCascadesDiscrete.cxx
  Correlationpp/AliAnalysisCorrelationEventCollection.cxx
  Correlationpp/AliAnalysisTaskCorrelationhhK0s.cxx
  Correlationpp/AliAnalysisTaskCorrelationhCasc.cxx
  HeavyIonsRun1Run2/AliAnalysisTaskStrVsMult.cxx
  )

if(ROOT_VERSION_MAJOR EQUAL 6)
	set(SRCS
    ${SRCS}
    Ratios/AliAnalysisTaskStrangenessRatios.cxx
  )
endif()

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Add header only files
set(HDRS
  ${HDRS}
  Lifetimes/MiniEvent.h
  Lifetimes/HyperTriton2Body.h
  Lifetimes/MCparticle.h
  Lifetimes/MiniV0.h
  Lifetimes/Utils.h
)

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

set(ROOT_DEPENDENCIES)
set(ALIROOT_DEPENDENCIES ANALYSISalice CORRFW OADB PWGUDbase)

# Generate the ROOT map
# Dependecies
set(LIBDEPS ${ALIROOT_DEPENDENCIES} ${ROOT_DEPENDENCIES})
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Generate a PARfile target for this library. Note the extra includes ("Tracks UserTasks")
add_target_parfile(${MODULE} "${SRCS}" "${HDRS}" "${MODULE}LinkDef.h" "${LIBDEPS}" "Tracks UserTasks")

# Create an object to be reused in case of static libraries
# Otherwise the sources will be compiled twice
add_library(${MODULE}-object OBJECT ${SRCS} G__${MODULE}.cxx)

# Add a library to the project using the object
add_library_tested(${MODULE} SHARED $<TARGET_OBJECTS:${MODULE}-object>)
target_link_libraries(${MODULE} ${ALIROOT_DEPENDENCIES} ${ROOT_DEPENDENCIES})

# Setting the correct headers for the object as gathered from the dependencies
target_include_directories(${MODULE}-object PUBLIC $<TARGET_PROPERTY:${MODULE},INCLUDE_DIRECTORIES>)
set_target_properties(${MODULE}-object PROPERTIES COMPILE_DEFINITIONS $<TARGET_PROPERTY:${MODULE},COMPILE_DEFINITIONS>)

# Public include folders that will be propagated to the dependecies
target_include_directories(${MODULE} PUBLIC ${incdirs})

set(MODULE_COMPILE_FLAGS)
set(MODULE_LINK_FLAGS)

# Additional compilation and linking flags
set(MODULE_COMPILE_FLAGS " ${MODULE_COMPILE_FLAGS}")

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set(MODULE_LINK_FLAGS "-undefined dynamic_lookup ${MODULE_LINK_FLAGS}")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Setting compilation flags for the object
set_target_properties(${MODULE}-object PROPERTIES COMPILE_FLAGS "${MODULE_COMPILE_FLAGS}")
# Setting the linking flags for the library
set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

# Installation
install(TARGETS ${MODULE}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

install(FILES ${HDRS} DESTINATION include)

# install macros
install(DIRECTORY Cascades/macros DESTINATION PWGLF/STRANGENESS/Cascades)
install(DIRECTORY Cascades/Run2/macros DESTINATION PWGLF/STRANGENESS/Cascades/Run2)
install(DIRECTORY Correlations/macros DESTINATION PWGLF/STRANGENESS/Correlations)
install(DIRECTORY Lifetimes/macros DESTINATION PWGLF/STRANGENESS/Lifetimes)
install(DIRECTORY LambdaK0/macros DESTINATION PWGLF/STRANGENESS/LambdaK0)
install(DIRECTORY LambdaK0PbPb DESTINATION  PWGLF/STRANGENESS)
install(DIRECTORY Correlationpp DESTINATION  PWGLF/STRANGENESS)
install(DIRECTORY HeavyIonsRun1Run2/macros DESTINATION PWGLF/STRANGENESS/HeavyIonsRun1Run2)
install(DIRECTORY Ratios DESTINATION PWGLF/STRANGENESS/ FILES_MATCHING PATTERN "*.C")

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/Cascades/corrections)
  install(DIRECTORY Cascades/corrections DESTINATION PWGLF/STRANGENESS/Cascades)
endif()
