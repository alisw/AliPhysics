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
set(MODULE PWGCFebye)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWGCF/EBYE)

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliPhysics_SOURCE_DIR}/EVENTMIX
                    ${AliPhysics_SOURCE_DIR}/CORRFW
                    ${AliPhysics_SOURCE_DIR}/PWG/Tools
                    ${AliPhysics_SOURCE_DIR}/PWGCF/EBYE/LRC
                    ${AliPhysics_SOURCE_DIR}/PWGPP
		    ${AliPhysics_SOURCE_DIR}/PWGPP/EVCHAR/FlowVectorCorrections/QnCorrections
                    ${AliPhysics_SOURCE_DIR}/PWGPP/EVCHAR/FlowVectorCorrections/QnCorrectionsInterface
                    ${AliPhysics_SOURCE_DIR}/OADB
                    ${AliPhysics_SOURCE_DIR}/OADB/COMMON/MULTIPLICITY
  )

# Sources - alphabetical order
set(SRCS
  BalanceFunctions/AliAnalysisTaskToyModel.cxx
  BalanceFunctions/AliAnalysisTaskBF.cxx
  BalanceFunctions/AliAnalysisTaskBFPsi.cxx
  BalanceFunctions/AliAnalysisTaskSignedBF.cxx
  BalanceFunctions/AliAnalysisTaskSignedBFMC.cxx
  BalanceFunctions/AliAnalysisTaskEventMixingBF.cxx
  BalanceFunctions/AliAnalysisTaskTriggeredBF.cxx
  BalanceFunctions/AliAnalysisTaskEfficiencyBF.cxx
  BalanceFunctions/AliAnalysisTaskEfficiencyBFPsi.cxx
  BalanceFunctions/AliAnalysisTaskEffContBF.cxx
  BalanceFunctions/AliAnalysisTaskEffContPIDBF.cxx
  BalanceFunctions/AliAnalysisTaskRecEff.cxx
  BalanceFunctions/AliBalance.cxx
  BalanceFunctions/AliBalancePsi.cxx
  BalanceFunctions/AliBalanceEbyE.cxx
  BalanceFunctions/AliBalanceEventMixing.cxx
  BalanceFunctions/AliBalanceTriggered.cxx
  BalanceFunctions/AliAnalysisTaskAODFilterBitQA.cxx
  BalanceFunctions/AliAnalysisTaskParticleStudies.cxx
  BalanceFunctions/AliAnalysisTaskPIDBF.cxx
  BalanceFunctions/AliPidBFBase.cxx
  BalanceFunctions/Corrections/AliAnalysisTaskAccCont.cxx
  BalanceFunctions/Corrections/AliAnalysisTaskPIDMCEffCont.cxx
  BalanceFunctions/Corrections/AliAnalysisTaskPIDPerformCombIDPtDep.cxx
  BalanceFunctions/HFtest/AliAnalysisTaskLcToD0MC.cxx
  KaonIsospinFluctuations/task/AliAnalysisTaskFluctMCTEPOS.cxx
  LRC/AliLRCBase.cxx
  LRC/AliAnalysisTaskLRC.cxx
  LRC/AliLRCAnalysis.cxx
  LRC/AliLRCFit.cxx
  LRC/AliLRCNN.cxx
  LRC/AliLRCPtN.cxx
  LRC/AliLRCPtPt.cxx
  LRC/AliLRCProcess.cxx
  LRC/AliForwardBackwardAnalysis.cxx
  Fluctuations/AliEbyEFluctuationAnalysisTask.cxx
  Fluctuations/AliAnalysisTaskChargeFluctuations.cxx
  Fluctuations/AliEbyEHigherMomentsTask.cxx
  Fluctuations/AliEbyEHigherMomentsEffContTask.cxx
  Fluctuations/AliHigherMomentsToyModel.cxx
  Fluctuations/AliEbyEMultFluctuationTask.cxx
  Fluctuations/AliEbyEHigherMomentsTaskPID.cxx
  Fluctuations/AliEbyEPidEfficiencyContamination.cxx
  MeanPtFluctuations/AliAnalysisTaskPtFluc.cxx
  MeanPtFluctuations/AliAnalysisTaskPtFlucPbPb.cxx
  MeanPtFluctuations/pp13TeV_MB/AliAnalysisMeanPt.cxx
  MeanPtFluctuations/pp13TeV_MB/AliAnalysisMeanPtFluc.cxx
  MeanPtFluctuations/pp13TeV_MB/AliAnalysisHMMeanPt.cxx
  MeanPtFluctuations/pp13TeV_MB/AliAnalysisMCMeanPt.cxx
  MeanPtFluctuations_PbPb5020GeV2018/AliAnalysisMeanPtdata.cxx
  NetChargeFluctuations/tasks/AliAnalysisTaskEbyeCharge.cxx
  NetChargeFluctuations/tasks/AliAnalysisTaskEbyeNetChargeFluctuations.cxx
  NetChargeFluctuations/tasks/AliAnalysisTaskEbyeNetChargeMCPbPbESD.cxx
  NetChargeFluctuations/tasks/AliAnalysisTaskEbyeChargeFlucPbPb.cxx
  NetChargeFluctuations/tasks/AliAnalysisTaskEbyeChargeFlucpPb.cxx
  NetChargeFluctuations/tasks/AliAnalysisTaskEbyeChargeFlucpp.cxx
  NetChargeFluctuations/tasks/AliAnalysisTaskEbyeChargeFlucPbPbQA.cxx
  NetChargeFluctuations/tasks/AliAnalysisTaskEbyeChargeFlucLHC10h.cxx
  NetLambda/AliAnalysisTaskNetLambdaIdent.cxx
  NetLambda/AliAnalysisTaskNetLambdaMCTrad.cxx
  NetLambda/AliAnalysisTaskNetLambdaTrad.cxx
  PIDFluctuation/task/AliAnalysisTaskPIDFluctuation.cxx
  PIDFluctuation/task/AliEbyENetChargeFluctuationTask.cxx
  PIDFluctuation/task/AliEbyEParticleRatioFluctuationTask.cxx
  PIDFluctuation/task/AliEbyEPidTaskFastGen.cxx
  PIDFluctuation/task/AliEbyEPidTTask.cxx
  PIDFluctuation/task/AliEbyEPidTTaskExPid.cxx
  PIDFluctuation/task/AliEbyEPidTTaskMC.cxx
  PIDFluctuation/task/AliEbyEPidRatioBase.cxx
  PIDFluctuation/task/AliEbyEPidRatioDCA.cxx
  PIDFluctuation/task/AliEbyEPidRatioPhy.cxx
  PIDFluctuation/task/AliEbyEPidRatioEffCont.cxx
  PIDFluctuation/task/AliEbyEPidRatioEffContExtra.cxx
  PIDFluctuation/task/AliEbyEPidRatioHelper.cxx
  PIDFluctuation/task/AliEbyEPidRatioQA.cxx
  PIDFluctuation/task/AliEbyEPidRatioTaskOnFly.cxx
  PIDFluctuation/task/AliEbyEPidQATaskOnFlyKine.cxx
  PIDFluctuation/task/AliEbyEPidRatioTask.cxx
  IdentityMethodEbyeFluctuations/tasks/AliAnalysisTaskEbyeIterPID.cxx
  NetParticle/AliAnalysisNetParticleBase.cxx
  NetParticle/AliAnalysisNetParticleDCA.cxx
  NetParticle/AliAnalysisNetParticleDistribution.cxx
  NetParticle/AliAnalysisNetParticleEffCont.cxx
  NetParticle/AliAnalysisNetParticleHelper.cxx
  NetParticle/AliAnalysisTaskNetParticle.cxx
  NetParticle/AliAnalysisNetParticleQA.cxx
  NetParticle/AliEbyEPhiDistNew.cxx
  TempFluctuations/AliAnalysisTempFluc.cxx
  LongAsymmetry/AliAnalysisFBMultFluct.cxx
  MultPt/AliAnalysisMultPt.cxx
  LongAsymmetry/AliAnalysisTaskLongFluctuations.cxx
  LongAsymmetry/AliAnalysisTaskLegendreCoef.cxx
  LongAsymmetry/AliAnalysisTaskLegendreCoef_local.cxx
  LongAsymmetry/AliAnalysisTaskEtaDist.cxx
  StronglyIntensiveCorrelations/AliAnalysisTaskStronglyIntensiveCorrTree.cxx
  FBCorrelationsWithPID/AliAnalysisFBCorrelationsWithPID.cxx
  ParticleYieldRatioCorrelations/AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency.cxx
  ParticleYieldRatioCorrelations/AliAnalysisTaskParticleYieldRatioCorrelations.cxx
  ################################
  HMPTF/AliAnalysisTaskResonanceVsMultiplicity.cxx
  HMPTF/AliAnalysisTaskMeanptFluctuationPP.cxx
  HMPTF/AliAnalysisTaskResonanceVsMultiplicityMC.cxx
  HMPTF/AliAnalysisTaskResonanceVsMultiplicity2.cxx
  HMPTF/AliAnalysisTaskResonanceVsMultiplicityROOT6.cxx
  HMPTF/AliAnalysisTaskResonanceVsMultiplicity_MCnoESD.cxx
  HMPTF/AliAnalysisTaskResonanceVsMultiplicity_MCnoESD_v2.cxx
  HMPTF/AliAnalysisTaskResonanceVsMultiplicity_MCnoESD_trackinfo.cxx
  HMPTF/AliAnalysisTaskCorrPbPb.cxx
  HMPTF/AliAnalysisTaskCorrPbPb_withoutTree.cxx
  HMPTF/AliAnalysisTaskCorrPbPb_treemerge.cxx
  HMPTF/AliAnalysisTaskCorrPbPb_treemerge_v2.cxx
  HMPTF/AliAnalysisTaskCorrPbPbMC.cxx
  HMPTF/AliAnalysisTaskResonanceVsMultiplicityROOT6_new.cxx
  HMPTF/AliAnalysisTaskNetPiKProtFluc_MCnoESD_gen.cxx
  HMPTF/AliAnalysisTaskMeanptFluctMC.cxx
  HMPTF/AliAnalysisTaskMeanptFluctMC_effcheck.cxx
  HMPTF/AliAnalysisTaskNetProtonCumulants_pp.cxx
  HMPTF/AliAnalysisTaskEffCalc_MC.cxx
  DiffPTF/AliAnalysisTaskDiffPtFluc_PIDhadrons.cxx
  DiffPTF/AliAnalysisTaskDiffPtFluc_PIDhadrons_v2.cxx
  DiffPTF/AliAnalysisTaskDiffPtFluc_PiKaProtHad_pTmax5.cxx
  DiffPTF/AliAnalysisTask_v02pT.cxx
  DiffPTF/AliAnalysisTaskDiffPtFluc_PiKaProtHad_pTmax5_HalfPhi.cxx
  DiffPTF/AliAnalysisTaskDiffPtFluc_PiKaProtHad_pTmax5_v2.cxx
  DiffPTF/AliAnalysisTaskDiffPtFluc_PIDhadrons_ParticleLossEffect.cxx
  DiffPTF/AliAnalysisTaskDiffPtFluc_MCnoESD_detect.cxx
  DiffPTF/AliAnalysisTaskDiffPtFluc_MCnoESD_gen.cxx
  DiffPTF/AliAnalysisTaskDiffPtFluc_MCnoESD_correction.cxx
  DiffPTF/AliAnalysisTaskDiffPtFluc_MCnoESD_detect_v2.cxx
  DiffPTF/AliAnalysisTaskDiffPtFluc_MCnoESD_gen_v2.cxx
  DiffPTF/AliAnalysisTaskDiffPtFluc_MCnoESD_correct_v2.cxx
  DiffPTF/AliAnalysisTaskHptcorrId_tree.cxx
  DiffPTF/AliAnalysisTaskDiffPtFluc_V0particles_pTmax5.cxx
  DiffPTF/AliAnalysisTaskDiffPtFluc_V0particles_pTmax5_v2.cxx
  IntermittencyAnalysis/AliAnalysisTaskNFactorialMoments.cxx
  IntermittencyAnalysis/AliAnalysisTaskNFactorialMomentsPID.cxx
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
set(LIBDEPS ANALYSISalice OADB PWGTools PWGPPevcharQnInterface EventMixing)
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

# Install the macros
install(DIRECTORY macros DESTINATION PWGCF/EBYE)
install(DIRECTORY PIDFluctuation/macros DESTINATION PWGCF/EBYE/PIDFluctuation)
install(DIRECTORY NetParticle/macros DESTINATION PWGCF/EBYE/NetParticle)
install(DIRECTORY NetLambda/macros DESTINATION PWGCF/EBYE/NetLambda)
install(DIRECTORY TempFluctuations/macros DESTINATION PWGCF/EBYE/TempFluctuations)
install(DIRECTORY LongAsymmetry/macros DESTINATION PWGCF/EBYE/LongAsymmetry)
install(DIRECTORY MultPt/macros DESTINATION PWGCF/EBYE/MultPt)
install(DIRECTORY IdentityMethodEbyeFluctuations/macros DESTINATION PWGCF/EBYE/IdentityMethodEbyeFluctuations)
install(DIRECTORY NetChargeFluctuations/macros DESTINATION PWGCF/EBYE/NetChargeFluctuations)
install(DIRECTORY StronglyIntensiveCorrelations DESTINATION PWGCF/EBYE/StronglyIntensiveCorrelations)
install(DIRECTORY FBCorrelationsWithPID/macros DESTINATION PWGCF/EBYE/FBCorrelationsWithPID)
install(DIRECTORY KaonIsospinFluctuations/macros DESTINATION PWGCF/EBYE/KaonIsospinFluctuations)
install(DIRECTORY ParticleYieldRatioCorrelations/macros DESTINATION PWGCF/EBYE/ParticleYieldRatioCorrelations)
install(DIRECTORY HMPTF/macros DESTINATION PWGCF/EBYE/HMPTF)
install(DIRECTORY DiffPTF/macros DESTINATION PWGCF/EBYE/DiffPTF)
