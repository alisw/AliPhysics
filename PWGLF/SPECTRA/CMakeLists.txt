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
set(MODULE PWGLFspectra)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWGLF/SPECTRA)

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliPhysics_SOURCE_DIR}/CORRFW
                    ${AliPhysics_SOURCE_DIR}/OADB
                    ${AliPhysics_SOURCE_DIR}/OADB/COMMON/MULTIPLICITY
                    ${AliPhysics_SOURCE_DIR}/PWG/Tools
                    ${AliPhysics_SOURCE_DIR}/PWGLF/SPECTRA/XtAnalysis
                    ${AliPhysics_SOURCE_DIR}/PWGCF/Correlations/JCORRAN/Base
                    ${AliPhysics_SOURCE_DIR}/PWGUD/base
                    ${AliPhysics_SOURCE_DIR}/EVENTMIX
                    ${AliPhysics_SOURCE_DIR}/PWG/TRD
  )

# Sources - alphabetical order
set(SRCS
  ChargedHadrons/dNdPt/AlidNdPtHelper.cxx
  ChargedHadrons/dNdPt/AliAnalysisTaskCutTest.cxx
  ChargedHadrons/dNdPt/AlidNdPtUnifiedAnalysisTask.cxx
  ChargedHadrons/dNdPt/MeanPtRaw/AliAnalysisTaskMeanPtRaw.cxx
  ChargedHadrons/dNdPt/AlidNdPtAcceptanceCuts.cxx
  ChargedHadrons/dNdPt/AlidNdPtAnalysis.cxx
  ChargedHadrons/dNdPt/AlidNdPtAnalysisPbPb.cxx
  ChargedHadrons/dNdPt/AlidNdPtAnalysisPbPb2011.cxx
  ChargedHadrons/dNdPt/AlidNdPtAnalysisPbPbAOD.cxx
  ChargedHadrons/dNdPt/AlidNdPtAnalysispPb.cxx
  ChargedHadrons/dNdPt/AlidNdPtBackgroundCuts.cxx
  ChargedHadrons/dNdPt/AlidNdPtCorrection.cxx
  ChargedHadrons/dNdPt/AlidNdPtCutAnalysis.cxx
  ChargedHadrons/dNdPt/AlidNdPtCutAnalysisPbPb.cxx
  ChargedHadrons/dNdPt/AlidNdPt.cxx
  ChargedHadrons/dNdPt/AlidNdPtEfficiency.cxx
  ChargedHadrons/dNdPt/AlidNdPtEventCuts.cxx
  ChargedHadrons/dNdPt/AlidNdPtHelper.cxx
  ChargedHadrons/dNdPt/AlidNdPtTask.cxx
  ChargedHadrons/AliMCSpectraWeights/AliMCSpectraWeightsAnalysisTask.cxx
  ChargedHadrons/dNdPt/Framework/AlidNdPtTools.cxx
  ChargedHadrons/dNdPt/Framework/AliAnalysisTaskMKBase.cxx
  ChargedHadrons/dNdPt/Framework/AliAnalysisTaskMKTest.cxx
  ChargedHadrons/dNdPt/Framework/AliAnalysisTaskTPCMatchEff.cxx
  ChargedHadrons/dNdPt/Framework/AliAnalysisTaskDCArStudy.cxx
  ChargedHadrons/dNdPt/Framework/AliAnalysisTaskPtResStudy.cxx
  ChargedHadrons/dNdPt/Framework/AliAnalysisTaskEffContStudy.cxx
  ChargedHadrons/dNdPt/Framework/AliAnalysisTaskUEStudy.cxx
  ChargedHadrons/dNdPt/Framework/AliAnalysisTaskSpectra.cxx
  ChargedHadrons/dNdPt/Framework/AliAnalysisTaskSpectraV0M.cxx
  ChargedHadrons/dNdPt/Framework/AliAnalysisTaskFilterEventTPCdEdx.cxx
  ChargedHadrons/dNdPt/AlidNdPtTrackDumpTask.cxx
  ChargedHadrons/dNdPt/AliPtResolAnalysis.cxx
  ChargedHadrons/dNdPt/AliPtResolAnalysisPbPb.cxx
  ChargedHadrons/dNdPt/AlimPtAnalysis.cxx
  ChargedHadrons/dNdPt/AlimPtMatrix.cxx
  ChargedHadrons/dNdPt/AlimPtMatrixZDC.cxx
  ChargedHadrons/dNdPt/AliAnalysisTaskSpectraDPhi.cxx
  ChargedHadrons/dNdPtVsMultpp/AliAnalysisTaskPPvsMultINEL0.cxx
  PiKaPr/ITSsa/AliAnalysisTaskSEITSsaSpectra.cxx
  PiKaPr/ITSsa/AliAnalysisTaskSEITSsaSpectraMultiplicity.cxx
  PiKaPr/TPC/Run3/AliAnalysisTaskSpectraTPCRun3.cxx
  PiKaPr/TPC/rTPC/AliAnalysisTaskQAHighPtDeDxTest.cxx
  PiKaPr/TPC/rTPC/AliAnalysisTaskPPvsMult.cxx
  PiKaPr/TPC/rTPC/AliAnalysisTaskSpherocity.cxx
  PiKaPr/TPC/rTPC/AliAnalysisTaskPPvsRT.cxx
  PiKaPr/TPC/rTPC/AliAnalysisTaskPPvsRT_TPCTOF.cxx
  PiKaPr/TPC/rTPC/AliAnalysisTaskSpectraRT.cxx
  PiKaPr/TPC/rTPC/AliAnalysisTaskSpectraMC.cxx
  PiKaPr/TPCTOF/AliAnalysisCombinedHadronSpectra.cxx
  PiKaPr/TPCTOFpA/AliAnalysisTPCTOFpA.cxx
  PiKaPr/TPCTOFfits/AliAnalysisPIDEvent.cxx
  PiKaPr/TPCTOFfits/AliAnalysisPIDParticle.cxx
  PiKaPr/TPCTOFfits/AliAnalysisPIDTrack.cxx
  PiKaPr/TPCTOFfits/AliAnalysisPIDV0.cxx
  PiKaPr/TPCTOFfits/AliAnalysisTaskTPCTOFPID.cxx
  PiKaPr/TPCTOFCascades/AliAnalysisPIDCascadeEvent.cxx
  PiKaPr/TPCTOFCascades/AliAnalysisPIDCascadeParticle.cxx
  PiKaPr/TPCTOFCascades/AliAnalysisPIDCascadeTrack.cxx
  PiKaPr/TPCTOFCascades/AliAnalysisPIDCascadeV0.cxx
  PiKaPr/TPCTOFCascades/AliAnalysisPIDCascade.cxx
  PiKaPr/TPCTOFCascades/AliAnalysisTaskTPCTOFCascade.cxx
  PiKaPr/TOF/pp13/AliAnalysisTaskTOFppSpectra.cxx
  PiKaPr/TOF/pp13/AliAnalysisTaskTOFMC.cxx
  PiKaPr/TOF/pp7/TOFSpectrappAnalysis.cxx
  PiKaPr/TOF/pp7/AliAnalysisCombinedHadronSpectra2MC.cxx
  PiKaPr/TOF/PbPb502/task/AliUtilTOFParams.cxx
  PiKaPr/TOF/PbPb502/task/AliAnTOFtrack.cxx
  PiKaPr/TOF/PbPb502/task/AliAnTOFevent.cxx
  PiKaPr/TOF/PbPb502/task/AliAnalysisTaskTOFSpectra.cxx
  ParticleEfficiency/AliAnalysisTaskParticleEfficiencyLF.cxx
  PiKaPr/TestAOD/AliSpectraAODTrackCuts.cxx
  PiKaPr/TestAOD/AliSpectraAODEventCuts.cxx
  PiKaPr/TestAOD/AliSpectraAODHistoManager.cxx
  PiKaPr/TestAOD/AliSpectraAODPID.cxx
  PiKaPr/TestAOD/AliAnalysisTaskSpectraAOD.cxx
  PiKaPr/TestAOD/AliAnalysisTaskSpectraAllChAOD.cxx
  PiKaPr/TestAOD/AliAnalysisTaskV2AllChAOD.cxx
  Kinks/AliAnalysisKinkESDat.cxx
  Kinks/AliAnalysisKinkESDMC.cxx
  Kinks/AliESDkinkCuts.cxx
  Kinks/AliAnalysisTaskKinksFilimon.cxx
  Kinks/AliAnalysisKinkESDat13.cxx
  Kinks/AliAnalysisKinkESDMC13.cxx
  Kinks/AliAnalysisKinkTaskMult13pp.cxx
  Kinks/AliAnalysisKinkTaskMult13ppMC.cxx
  Kinks/AliAnalysisTaskKinkpp5TeV.cxx
  Kinks/AliAnalysisTaskKinkpp5TeVMC.cxx
  Kinks/AliAnalysisTaskKinkPbPb.cxx
  Kinks/AliAnalysisTaskKinkPbPbMC.cxx
  Kinks/AliAnalysisKinkTaskPbPb2018.cxx
  PiKaPr/TestAOD/AliSpectraBothTrackCuts.cxx
  PiKaPr/TestAOD/AliSpectraBothEventCuts.cxx
  PiKaPr/TestAOD/AliSpectraBothHistoManager.cxx
  PiKaPr/TestAOD/AliSpectraBothPID.cxx
  PiKaPr/TestAOD/AliAnalysisTaskSpectraBoth.cxx
  PiKaPr/TestAOD/AliMultiplictyLoaderTask.cxx
  XtAnalysis/AliXtAnalysis.cxx
  XtAnalysis/AliJXtHistos.cxx
  IdentifiedHighPt/train/AliAnalysisTaskHighPtDeDx.cxx
  IdentifiedHighPt/train/AliAnaTaskV0EffDecomposition.cxx
  IdentifiedHighPt/train/DebugClassesMultESA2013.cxx
  MultEvShape/AliMESeventInfo.cxx
  MultEvShape/AliMEStrackInfo.cxx
  MultEvShape/AliMESbaseTask.cxx
  MultEvShape/AliMEStender.cxx
  MultEvShape/AliMESpidTask.cxx
  MultEvShape/AliMESchgTask.cxx
  MultEvShape/AliMESppColTask.cxx
  AntiprotonToProton/AliProtonAnalysis.cxx
  AntiprotonToProton/AliProtonAnalysisBase.cxx
  AntiprotonToProton/train/AliAnalysisTaskProton.cxx
  AntiprotonToProton/train/AliAnalysisTaskProtonQA.cxx
  AntiprotonToProton/train/AliAnalysisTaskProtonAbsorbtion.cxx
  Spherocity/AliTransverseEventShape.cxx
  Spherocity/AliSpherocityEstimator.cxx
  Spherocity/AliSpherocityUtils.cxx
  Spherocity/macros/AliAnalysisSphericityTask.cxx
)
if(${ROOT_VERSION} GREATER_EQUAL 6.0)
  list(APPEND SRCS ChargedHadrons/MultDepSpec/AliMultDepSpecAnalysisTask.cxx)
  list(APPEND SRCS ChargedHadrons/MultDepSpec/AliMultDepSpecAnalysisTaskUE.cxx)
  list(APPEND SRCS ChargedHadrons/dNdPt/Framework/AliAnalysisTaskCutStudies.cxx)
  list(APPEND SRCS ChargedHadrons/dNdPt/Framework/AliAnalysisTaskSpectraEtaPhi.cxx)
  list(APPEND SRCS ChargedHadrons/dNdPt/Framework/AliAnalysisTaskSpectraINEL0.cxx)
  list(APPEND SRCS ChargedHadrons/dNdPt/Framework/AliAnalysisTaskBaseWeights.cxx)
  list(APPEND SRCS ChargedHadrons/AliMCSpectraWeights/AliMCWeightsTask.cxx)
endif()

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Extra files to include in the PARfile
set(PAR_EXTRAFILES_CXX
  PiKaPr/TestAOD/HistogramNamesBoth.cxx
  PiKaPr/TestAOD/HistogramNames.cxx
)

set(PAR_EXTRAFILES_H
  PiKaPr/TestAOD/HistogramNamesBoth.h
  PiKaPr/TestAOD/HistogramNames.h
  PiKaPr/TestAOD/HistogramsBoth.h
  PiKaPr/TestAOD/Histograms.h
)

set(PAR_EXTRAFILES
   ${PAR_EXTRAFILES_H}
   ${PAR_EXTRAFILES_CXX}
)

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs} ${PAR_EXTRAFILES_H}")

set(ROOT_DEPENDENCIES)
set(ALIROOT_DEPENDENCIES ANALYSISalice PWGTools PWGUDbase EventMixing PWGCFCorrelationsJCORRAN PWGLFnuclex)

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
install(FILES PiKaPr/TestAOD/HistogramsBoth.h PiKaPr/TestAOD/Histograms.h DESTINATION include)

if(${ROOT_VERSION} GREATER_EQUAL 6.0)
  install(FILES ChargedHadrons/MultDepSpec/AliAnalysisHelpersHist.h DESTINATION include)
endif()

# install macros
install(DIRECTORY ChargedHadrons/MultDepSpec/macros DESTINATION  PWGLF/SPECTRA/ChargedHadrons/MultDepSpec)
install(DIRECTORY ChargedHadrons/dNdPt/macros DESTINATION  PWGLF/SPECTRA/ChargedHadrons/dNdPt)
install(DIRECTORY ChargedHadrons/dNdPt/MeanPtRaw/macros DESTINATION
  PWGLF/SPECTRA/ChargedHadrons/dNdPt/MeanPtRaw)
install(DIRECTORY ChargedHadrons/dNdPtVsMultpp/macros DESTINATION  PWGLF/SPECTRA/ChargedHadrons/dNdPtVsMultpp)
install(FILES ChargedHadrons/AliMCSpectraWeights/AddTask_phuhn_MCSpectraWeights.C DESTINATION PWGLF/SPECTRA/ChargedHadrons/AliMCSpectraWeights)
install(FILES ChargedHadrons/AliMCSpectraWeights/AddTask_MCWeights.C DESTINATION PWGLF/SPECTRA/ChargedHadrons/AliMCSpectraWeights)
install(DIRECTORY IdentifiedHighPt/train/macros DESTINATION  PWGLF/SPECTRA/IdentifiedHighPt/train)
install(DIRECTORY IdentifiedHighPt/train/macros1 DESTINATION  PWGLF/SPECTRA/IdentifiedHighPt/train)
install(DIRECTORY Kinks/macros DESTINATION  PWGLF/SPECTRA/Kinks)
install(DIRECTORY PiKaPr/TestAOD DESTINATION PWGLF/SPECTRA/PiKaPr)
install(DIRECTORY PiKaPr/ITSsa DESTINATION PWGLF/SPECTRA/PiKaPr)
install(DIRECTORY PiKaPr/TPC/rTPC DESTINATION PWGLF/SPECTRA/PiKaPr/TPC)
install(DIRECTORY PiKaPr/TPC/Run3 DESTINATION PWGLF/SPECTRA/PiKaPr/TPC)
install(DIRECTORY PiKaPr/TPCTOFfits DESTINATION PWGLF/SPECTRA/PiKaPr)
install(DIRECTORY PiKaPr/TPCTOFCascades DESTINATION PWGLF/SPECTRA/PiKaPr)

install(FILES PiKaPr/TOF/PbPb502/task/AddTaskTOFSpectra.C DESTINATION  PWGLF/SPECTRA/PiKaPr/TOF)
install(FILES PiKaPr/TOF/PbPb502/task/AliUtilTOFParams.h DESTINATION  PWGLF/SPECTRA/PiKaPr/TOF/utilities)
install(FILES PiKaPr/TOF/PbPb502/task/AliUtilTOFParams.cxx DESTINATION  PWGLF/SPECTRA/PiKaPr/TOF/utilities)
install(FILES PiKaPr/TOF/pp13/AddTaskTOFppSpectra.C DESTINATION  PWGLF/SPECTRA/PiKaPr/TOF)
install(FILES PiKaPr/TOF/pp13/AddTaskTOFMC.C DESTINATION  PWGLF/SPECTRA/PiKaPr/TOF)
install(FILES PiKaPr/TPCTOF/AddTaskCombinedHadronSpectra.C DESTINATION  PWGLF/SPECTRA/PiKaPr/TPCTOF)
install(FILES PiKaPr/TPCTOFpA/AddTaskTPCTOFpA.C DESTINATION  PWGLF/SPECTRA/PiKaPr/TPCTOFpA)
install(DIRECTORY XtAnalysis/macros DESTINATION PWGLF/SPECTRA/XtAnalysis)
install(DIRECTORY MultEvShape DESTINATION PWGLF/SPECTRA FILES_MATCHING PATTERN "*.C")
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/MultEvShape/priorsDist_data_LHC10d_newAliroot.root)
  install(FILES MultEvShape/priorsDist_data_LHC10d_newAliroot.root DESTINATION PWGLF/SPECTRA/MultEvShape)
endif()
install(DIRECTORY Spherocity DESTINATION  PWGLF/SPECTRA)
install(DIRECTORY AntiprotonToProton/train DESTINATION PWGLF/SPECTRA/AntiprotonToProton)
