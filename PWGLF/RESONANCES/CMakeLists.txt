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
set(MODULE PWGLFresonances)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWGLF/RESONANCES)

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliPhysics_SOURCE_DIR}/EVENTMIX
                    ${AliPhysics_SOURCE_DIR}/CORRFW
                    ${AliPhysics_SOURCE_DIR}/OADB
                    ${AliPhysics_SOURCE_DIR}/PWG/Tools
                    ${AliPhysics_SOURCE_DIR}/PWG/DevNanoAOD
                    ${AliPhysics_SOURCE_DIR}/OADB/COMMON/MULTIPLICITY
                    ${AliPhysics_SOURCE_DIR}/PWGPP/EVCHAR/FlowVectorCorrections/QnCorrections
                    ${AliPhysics_SOURCE_DIR}/PWGPP/EVCHAR/FlowVectorCorrections/QnCorrectionsInterface
  )

# Sources - alphabetical order
set(SRCS
  AliRsnDaughter.cxx
  AliRsnMother.cxx
  AliRsnEvent.cxx
  AliRsnDaughterDef.cxx
  AliRsnPairDef.cxx
  AliRsnTarget.cxx
  AliRsnValue.cxx
  AliRsnAction.cxx
  AliRsnCut.cxx
  AliRsnValueDaughter.cxx
  AliRsnValuePair.cxx
  AliRsnValueEvent.cxx
  AliRsnCutEventUtils.cxx
  AliRsnEventCuts.cxx
  AliRsnCutPrimaryVertex.cxx
  AliRsnCutTrackQuality.cxx
  AliRsnCutTrue.cxx
  AliRsnCutValue.cxx
  AliRsnCutPID.cxx
  AliRsnCutPIDITS.cxx
  AliRsnCutPIDTPC.cxx
  AliRsnCutPIDTOF.cxx
  AliRsnPIDRange.cxx
  AliRsnCutPIDNSigma.cxx
  AliRsnCutMomentumComparison.cxx
  AliRsnCutTOFMatch.cxx
  AliRsnCutPhi.cxx
  AliRsnCutDelta.cxx
  AliRsnCutKaonForPhi2010.cxx
  AliRsnCutKaonForPhi2010PP.cxx
  AliRsnCutPion2010PP.cxx
  AliRsnCutProton2010PP.cxx
  AliRsnCutDaughterKStar2010PP.cxx 
  AliRsnCutDaughterLStar2010.cxx
  AliRsnCutDaughterSigmaStar2010PP.cxx 
  AliRsnCutSetDaughterParticle.cxx 
  AliRsnCutPhiRAA.cxx
  AliRsnCutDaughterD0.cxx
  AliRsnCutV0.cxx 
  AliRsnCutCascade.cxx 
  AliRsnCutSet.cxx
  AliRsnExpression.cxx
  AliRsnVariableExpression.cxx
  AliRsnCutManager.cxx
  AliRsnListOutput.cxx
  AliRsnLoop.cxx
  AliRsnLoopDaughter.cxx
  AliRsnLoopPair.cxx
  AliRsnLoopEvent.cxx
  AliRsnLoopEff.cxx
  AliRsnLoopEffPair.cxx
  AliRsnDaughterSelector.cxx
  AliRsnInputHandler.cxx
  AliRsnAnalysisTask.cxx
  AliRsnMiniParticle.cxx
  AliRsnMiniPair.cxx
  AliRsnCutMiniPair.cxx
  AliRsnMiniEvent.cxx
  AliRsnMiniAxis.cxx
  AliRsnMiniOutput.cxx
  AliRsnMiniValue.cxx
  AliRsnMiniMonitor.cxx
  AliRsnMiniAnalysisTask.cxx
  AliRsnMiniMonitorTask.cxx
  AliRsnMiniResonanceFinder.cxx
  AliRsnTrainManager.cxx
  AliRsnMiniTaskPhiVn.cxx
  AliRsnCutPIDkink.cxx
  PostProcessing/Xi1530/AliAnalysisTaskXi1530.cxx
  PostProcessing/Xi1530/NanoCheck/AliAnalysisTaskNanoCheck.cxx
  PostProcessing/Sigma1385/AliAnalysisTaskSigma1385PM.cxx
  PostProcessing/Sigma1385/AliAnalysisTaskTrackMixer.cxx
  PostProcessing/f0980/AliAnalysisTaskInclusivef0f2.cxx
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
if(ROOT_VERSION_MAJOR EQUAL 6)
    set(ROOT_DEPENDENCIES PyROOT)
else()
    set(ROOT_DEPENDENCIES)
endif()
set(ALIROOT_DEPENDENCIES ANALYSISalice CORRFW EventMixing PWGPPevcharQnInterface PWGTools PWGDevNanoAOD)
set(LIBDEPS ${ALIROOT_DEPENDENCIES} ${ROOT_DEPENDENCIES})
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

# Install macros
install(DIRECTORY macros DESTINATION PWGLF/RESONANCES)
install(DIRECTORY PostProcessing DESTINATION PWGLF/RESONANCES)
