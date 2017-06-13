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
set(MODULE PWGHFvertexingHF)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWGHF/vertexingHF)

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliPhysics_SOURCE_DIR}/CORRFW
		    ${AliPhysics_SOURCE_DIR}/OADB
		    ${AliPhysics_SOURCE_DIR}/OADB/COMMON/MULTIPLICITY
                    ${AliPhysics_SOURCE_DIR}/PWGPP/EVCHAR/FlowVectorCorrections/QnCorrections
		    ${AliPhysics_SOURCE_DIR}/PWGPP/EVCHAR/FlowVectorCorrections/QnCorrectionsInterface/
                    ${AliPhysics_SOURCE_DIR}/PWG/FLOW/Base
                    ${AliPhysics_SOURCE_DIR}/PWG/FLOW/Tasks
                    ${AliPhysics_SOURCE_DIR}/PWG/muon
                    ${AliPhysics_SOURCE_DIR}/PWG/TRD
  )

# Sources - alphabetical order
set(SRCS
  AliAnalysisTaskSEB0toDStarPi.cxx
  AliRDHFCutsB0toDStarPi.cxx
  AliAODHFUtil.cxx
  AliAODRecoDecayHF.cxx
  AliAODRecoDecayHF2Prong.cxx
  AliAODRecoDecayHF3Prong.cxx
  AliAODRecoDecayHF4Prong.cxx
  AliAODRecoCascadeHF.cxx
  AliAODRecoCascadeHF3Prong.cxx
  AliAODPidHF.cxx
  AliRDHFCuts.cxx
  AliVertexingHFUtils.cxx
  AliHFSystErr.cxx
  AliRDHFCutsD0toKpi.cxx
  AliRDHFCutsJpsitoee.cxx
  AliRDHFCutsDplustoKpipi.cxx
  AliRDHFCutsDplustoK0spi.cxx
  AliRDHFCutsDstoKKpi.cxx
  AliRDHFCutsDstoK0sK.cxx
  AliRDHFCutsLctopKpi.cxx
  AliRDHFCutsLctoV0.cxx
  AliRDHFCutsD0toKpipipi.cxx
  AliRDHFCutsDStartoKpipi.cxx
  AliRDHFCutsLctoeleLambdafromAODtracks.cxx
  AliRDHFCutsLctopK0sfromAODtracks.cxx
  AliRDHFCutsOmegactoeleOmegafromAODtracks.cxx
  AliRDHFCutsXicPlustoXiPiPifromAODtracks.cxx
  AliRDHFCutsXictoeleXifromAODtracks.cxx
  AliAnalysisVertexingHF.cxx
  AliAnalysisTaskSEVertexingHF.cxx
  AliAnalysisTaskMEVertexingHF.cxx
  AliAnalysisTaskSESelectHF.cxx
  AliAnalysisTaskSECleanupVertexingHF.cxx
  AliAnalysisTaskSECompareHF.cxx
  AliAnalysisTaskSEDplus.cxx
  AliAnalysisTaskSEDs.cxx
  AliAnalysisTaskSELambdac.cxx
  AliAnalysisTaskSED0Mass.cxx
  AliAnalysisTaskSELambdacTMVA.cxx
  AliAnalysisTaskSED0MassCutVarFDSub.cxx
  AliAnalysisTaskSECharmFraction.cxx
  AliAnalysisTaskSED0MassNonPromptFraction.cxx
  AliAnalysisTaskSELc2V0bachelor.cxx
  AliAnalysisTaskSEDvsMultiplicity.cxx
  AliAnalysisTaskCheckHFMCProd.cxx
  AliAnalysisTaskDmesonMCPerform.cxx
  AliAnalysisTaskSELc2pK0sfromAODtracks.cxx
  AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks.cxx
  AliCFVertexingHF.cxx
  AliCFVertexingHF2Prong.cxx
  AliCFVertexingHF3Prong.cxx
  AliCFVertexingHFCascade.cxx
  AliCFVertexingHFLctoV0bachelor.cxx
  AliCFTaskVertexingHF.cxx
  AliCFTaskVertexingHFCutVarFDSub.cxx
  AliCFTaskForDStarAnalysis.cxx
  AliAnalysisTaskSEDStarJets.cxx
  AliAnalysisTaskSEDStarSpectra.cxx
  AliAnalysisTaskSEDStarCharmFraction.cxx
  AliAnalysisTaskSESignificance.cxx
  AliAnalysisTaskSEHFQA.cxx
  AliAnalysisTaskTrackingSysPropagation.cxx
  AliMultiDimVector.cxx
  AliSignificanceCalculator.cxx
  AliHFMassFitter.cxx
  AliHFMassFitterVAR.cxx
  AliHFInvMassFitter.cxx
  AliHFMultiTrials.cxx
  AliHFInvMassMultiTrialFit.cxx
  AliHFPtSpectrum.cxx
  AliHFsubtractBFDcuts.cxx
  AliNormalizationCounter.cxx
  AliAnalysisTaskSEMonitNorm.cxx
  AliAnalysisTaskSEBkgLikeSignD0.cxx
  AliAnalysisTaskSEImproveITS.cxx
  AliAnalysisTaskCombinHF.cxx
  charmFlow/AliEventPlaneResolutionHandler.cxx
  charmFlow/AliAnalysisTaskSEHFv2.cxx
  charmFlow/AliAnalysisTaskSEHFvn.cxx
  charmFlow/AliAnalysisTaskHFv1.cxx
  charmFlow/AliHFAfterBurner.cxx
  charmFlow/AliAnalysisTaskFlowD2H.cxx
  cutVarFDsub/AliHFCutVarFDsubAnalysisManager.cxx
  cutVarFDsub/AliHFCutVarFDsubAnalysisManagerD0.cxx
  cutVarFDsub/AliHFCutVarFDsubAnalysisManagerDplus.cxx
  cutVarFDsub/AliHFCutVarFDsubAxis.cxx
  cutVarFDsub/AliHFCutVarFDsubCut.cxx
  cutVarFDsub/AliHFCutVarFDsubCutSet.cxx
  cutVarFDsub/AliHFCutVarFDsubEfficiency.cxx
  cutVarFDsub/AliHFCutVarFDsubMassFitter.cxx
  cutVarFDsub/AliHFCutVarFDsubMinimiser.cxx
  upgrade/AliAnalysisTaskSELambdacUp.cxx
  upgrade/AliAnalysisTaskCountLcEta.cxx
  AliAnalysisTaskSELc2eleLambdafromAODtracks.cxx
  AliAnalysisTaskSELc2V0bachelorTMVA.cxx
  AliAnalysisTaskSEOmegac2eleOmegafromAODtracks.cxx
  AliAnalysisTaskSEXic2eleXifromAODtracks.cxx
  AliAnalysisTaskCharmBaryonsMC.cxx
  AliAnalysisTaskSEDvsEventShapes.cxx
  AliAnalysisTaskSEB0toDminuspi.cxx
  AliAnalysisTaskSEDstoK0sK.cxx
  AliHFVnVsMassFitter.cxx
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
set(LIBDEPS ANALYSISalice PWGflowTasks PWGTRD PWGPPevcharQn PWGPPevcharQnInterface)
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

install(FILES
            AliAODRecoDecayHF.h
            AliAODRecoDecayHF2Prong.h
            AddD2HTrain.C
	    AddGoodRuns.C
            ConfigVertexingHF.C
            ConfigVertexingHF_pPbRHF.C
            ConfigVertexingHF_highmult.C
            ConfigVertexingHF_pp_LambdaC.C
            ConfigVertexingHF_TrackRel.C
            ConfigVertexingHF_Pb_Cent2080.C
            ConfigVertexingHF_Pb_AllCent.C
            ConfigVertexingHF_Pb_AllCent_NoLS_PIDLc.C
	    ConfigVertexingHF_Pb_AllCent_NoLS_PIDLc_PtDepSel_LooseIP.C
	    ConfigVertexingHF_Pb_AllCent_NoLS_PIDLc_PtDepSel_LooseIP_Central.C
	    ConfigVertexingHF_Pb_AllCent_NoLS_PIDLc_PtDepSel_LooseIP_SemiCentral.C
            ConfigVertexingHF_Pb_AllCent_NoLS_NoDstar_16var.C
            ConfigVertexingHF_Pb_AllCent_NoLS_NoDstar.C
            ConfigVertexingHF_Pb_AllCent_NoLS_16var.C
            ConfigVertexingHF_Pb_AllCent_NoLS.C
            ConfigVertexingHF_DplustoK0s_DstoK0s.C
            ConfigVertexingHF_pPbRHF_DplustoK0s_DstoK0s.C
        DESTINATION
            PWGHF/vertexingHF/)

install(DIRECTORY
            macros
        DESTINATION PWGHF/vertexingHF/)

install(DIRECTORY charmFlow DESTINATION PWGHF/vertexingHF)
