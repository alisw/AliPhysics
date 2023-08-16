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
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                    ${AliPhysics_SOURCE_DIR}/EVENTMIX
                    ${AliPhysics_SOURCE_DIR}/CORRFW
                    ${AliPhysics_SOURCE_DIR}/OADB
                    ${AliPhysics_SOURCE_DIR}/OADB/COMMON/MULTIPLICITY
                    ${AliPhysics_SOURCE_DIR}/PWG/Tools
                    ${AliPhysics_SOURCE_DIR}/PWGCF/Correlations/Base
                    ${AliPhysics_SOURCE_DIR}/PWGHF/vertexingHF
                    ${AliPhysics_SOURCE_DIR}/PWGLF/FORWARD/analysis2
                    ${AliPhysics_SOURCE_DIR}/PWGLF/RESONANCES
                    ${AliPhysics_SOURCE_DIR}/PWGLF/STRANGENESS/Cascades/lightvertexers
                    ${AliPhysics_SOURCE_DIR}/PWGLF/STRANGENESS/Cascades/Run2
                    ${AliPhysics_SOURCE_DIR}/PWGUD/base
  )

# pulled in for HF vertexing
if(FASTJET_FOUND)
    include_directories(SYSTEM ${FASTJET_INCLUDE_DIR})
    link_directories(${FASTJET_LIBS_DIR})
    add_definitions(${FASTJET_DEFINITIONS})
    list(APPEND HDRS FJ_includes.h)
endif(FASTJET_FOUND)

if(ROOT_VERSION_MAJOR EQUAL 6)
# Additional includes for ROOT 6 only dependencies - alphabetical order except ROOT
  include_directories(${AliPhysics_SOURCE_DIR}/ML
                    ${TREELITE_ROOT}/include
                    ${TREELITE_ROOT}/runtime/native/include
                    ${YAML_CPP_SOURCE_DIR}/include
    )
endif()

# Sources - alphabetical order
set(SRCS
  Cascades/AliAnalysisTaskCheckCascade.cxx
  Cascades/AliAnalysisTaskCheckCascadePbPb.cxx
  Cascades/AliAnalysisTaskCheckPerformanceCascade.cxx
  Cascades/AliAnalysisTaskCheckPerformanceCascadePbPb.cxx
  Cascades/AliAnalysisTaskXiNucleusInteraction.cxx
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
  Cascades/Run2/AliTaskLeadingMC.cxx
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
  Cascades/Run2/AliAnalysisTaskMCPredictions2pc.cxx
  Cascades/Run2/AliAnalysisTaskMCPredictionsEE.cxx
  Cascades/Run2/AliAnalysisTaskMCPredictionsStrgVsMultVsZDC.cxx
  Cascades/Run2/AliAnalysisTaskStrangeCascadesDiscrete.cxx
  Correlationpp/AliAnalysisCorrelationEventCollection.cxx
  Correlationpp/AliAnalysisTaskCorrelationhhK0s.cxx
  Correlationpp/AliAnalysisTaskCorrelationhCasc.cxx
  Correlationpp/AliAnalysisTaskCorrelationhCascDATA.cxx
  Correlationpp/AliAnalysisTaskCorrelationhK0sXi_PureMCOnly.cxx
  HeavyIonsRun1Run2/AliAnalysisTaskStrVsMult.cxx
  HeavyIonsRun1Run2/AliAnalysisTaskStrVsMult_BumpStudies.cxx
  HeavyIonsRun1Run2/AliAnalysisTaskStrVsMult_PileUpStudies.cxx
  HeavyIonsRun1Run2/AliAnalysisTaskStrVsMult_SelectionStudies.cxx
  HeavyIonsRun1Run2/AliAnalysisTaskStrVsMult_RStudies.cxx
  Cascades/Run2/AliEventContainer.cxx
  Cascades/Run2/AliCascadeContainer.cxx
  Cascades/Run2/AliV0Container.cxx
  Cascades/Run2/AliTrackInfoContainer.cxx
  Cascades/Run2/AliAnalysisTaskStrangeCascadesTriggerAODRun2.cxx
  Cascades/Run2/AliMCEventContainer.cxx
  Cascades/Run2/AliMCCascadeContainer.cxx
  Cascades/Run2/AliMCV0Container.cxx
  Cascades/Run2/AliMCTrackInfoContainer.cxx
  Cascades/Run2/AliMCGenParticleContainer.cxx
  Cascades/Run2/AliAnalysisTaskStrangeCascadesTriggerAODRun2MC.cxx
  StrPDF/AliAnalysisTaskV0multspec.cxx
  Cascades/Run2/AliAnalysisTaskStrangenessVsMultVsEffEnergyAODRun2.cxx
  Cascades/Run2/AliAnalysisTaskStrangenessVsMultVsEffEnergyAODMCRun2.cxx
  Cascades/Run2/AliMachineLearning.cxx
  Cascades/Run2/AliNeuralNetwork.cxx
  Cascades/Run2/AliBDT.cxx
  Cascades/Run2/AliAnalysisTaskStrangenessML.cxx
  StrPDF/AliAnalysisTaskMultspec.cxx
  StrPDF/AliAnalysisTaskMultspec_MCpred.cxx
  StrPDF/AliAnalysisTaskMultspec_evbased.cxx
  )

if(ROOT_VERSION_MAJOR EQUAL 6)
	set(SRCS
    ${SRCS}
    Ratios/AliAnalysisTaskStrangenessRatios.cxx
    KaonXiCorrelation/AliAnalysisTaskKaonXiCorrelation.cxx
  )
endif()

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

set(ROOT_DEPENDENCIES)
set(ALIROOT_DEPENDENCIES ANALYSISalice CORRFW OADB PWGUDbase PWGLFresonances PWGHFvertexingHF)

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
install(DIRECTORY LambdaK0/macros DESTINATION PWGLF/STRANGENESS/LambdaK0)
install(DIRECTORY LambdaK0PbPb DESTINATION  PWGLF/STRANGENESS)
install(DIRECTORY Correlationpp DESTINATION  PWGLF/STRANGENESS)
install(DIRECTORY HeavyIonsRun1Run2/macros DESTINATION PWGLF/STRANGENESS/HeavyIonsRun1Run2)
install(DIRECTORY StrPDF/macros DESTINATION PWGLF/STRANGENESS/StrPDF)
install(DIRECTORY Ratios DESTINATION PWGLF/STRANGENESS/ FILES_MATCHING PATTERN "*.C")
install(DIRECTORY KaonXiCorrelation DESTINATION PWGLF/STRANGENESS/ FILES_MATCHING PATTERN "*.C")

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/Cascades/corrections)
  install(DIRECTORY Cascades/corrections DESTINATION PWGLF/STRANGENESS/Cascades)
endif()
