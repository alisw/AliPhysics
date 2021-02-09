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
# **************************************************************************

# Module
set(MODULE PWGJEEMCALJetTasks)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWGJE/EMCALJetTasks
                    ${AliPhysics_SOURCE_DIR}/PWGJE/EMCALJetTasks/UserTasks
                   )

# Additional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    $ENV{KFPARTICLE_ROOT}/include
                    ${AliPhysics_SOURCE_DIR}/CORRFW
                    ${AliPhysics_SOURCE_DIR}/OADB
                    ${AliPhysics_SOURCE_DIR}/PWG/DevNanoAOD
                    ${AliPhysics_SOURCE_DIR}/PWG/EMCAL/EMCALbase
                    ${AliPhysics_SOURCE_DIR}/PWG/EMCAL/EMCALtrigger
                    ${AliPhysics_SOURCE_DIR}/PWG/JETFW
                    ${AliPhysics_SOURCE_DIR}/PWG/Tools
                    ${AliPhysics_SOURCE_DIR}/PWG/Tools/yaml-cpp/include
                    ${AliPhysics_SOURCE_DIR}/PWG/CaloTrackCorrBase
                    ${AliPhysics_SOURCE_DIR}/PWGHF/jetsHF
                    ${AliPhysics_SOURCE_DIR}/PWGHF/vertexingHF
                    ${AliPhysics_SOURCE_DIR}/PWGPP/vertexingHF
                    ${AliPhysics_SOURCE_DIR}/PWGPP/EVCHAR/FlowVectorCorrections/QnCorrections
                    ${AliPhysics_SOURCE_DIR}/PWGPP/EVCHAR/FlowVectorCorrections/QnCorrectionsInterface
                             )

# Sources in alphabetical order
set(SRCS
    AliAnalysisTaskDeltaPt.cxx
    AliAnalysisTaskEmcalJetSample.cxx
    AliAnalysisTaskEmcalJetQA.cxx
    AliAnalysisTaskEmcalJetSpectraQA.cxx
    AliAnalysisTaskEmcalJetTagger.cxx
    AliAnalysisTaskPWGJEQA.cxx
    AliAnalysisTaskLocalRho.cxx
    AliAnalysisTaskRhoAverage.cxx
    AliAnalysisTaskRhoBase.cxx
    AliAnalysisTaskRho.cxx
    AliAnalysisTaskRhoFlow.cxx
    AliAnalysisTaskRhoMassBase.cxx
    AliAnalysisTaskRhoMass.cxx
    AliAnalysisTaskRhoMassSparse.cxx
    AliAnalysisTaskRhoSparse.cxx
    AliAnalysisTaskJetUE.cxx
    AliAnalysisTaskRhoBaseDev.cxx
    AliAnalysisTaskRhoDev.cxx
    AliAnalysisTaskRhoTransDev.cxx
    AliAnalysisTaskScale.cxx
    AliEmcalJetByJetCorrection.cxx
    AliEmcalJetTaggerTaskFast.cxx
    AliEmcalPicoTrackInGridMaker.cxx
    AliJetConstituentTagCopier.cxx
    AliJetEmbeddingFromGenTask.cxx
    AliJetEmbeddingTask.cxx
    AliJetFastSimulation.cxx
    AliJetModelBaseTask.cxx
    AliJetModelCopyTracks.cxx
    AliJetModelMergeBranches.cxx
    AliJetRandomizerTask.cxx
    AliJetResponseMaker.cxx
    AliJetTriggerSelectionTask.cxx
    AliNanoAODArrayMaker.cxx
    AliPWGJETrainHelpers.cxx
    Tracks/AliAnalysisTaskEmcalTriggerBase.cxx
    Tracks/AliAnalysisTaskEmcalTriggerPosition.cxx
    Tracks/AliAnalysisTaskPtEMCalTrigger.cxx
    Tracks/AliAnalysisTaskPtEMCalTriggerV1.cxx
    Tracks/AliCutValueRange.cxx
    Tracks/AliEmcalTriggerOfflineSelection.cxx
    Tracks/AliEMCalTriggerExtraCuts.cxx
    Tracks/AliEMCalTriggerAnaTriggerDecision.cxx
    Tracks/AliEMCalTriggerAnaTriggerDecisionConfig.cxx
    Tracks/AliEMCalTriggerBinningComponent.cxx
    Tracks/AliEMCalTriggerBinningFactory.cxx
    Tracks/AliEMCalTriggerClusterAnalysisComponent.cxx
    Tracks/AliEMCalTriggerEventCounterAnalysisComponent.cxx
    Tracks/AliEMCalTriggerEventData.cxx
    Tracks/AliEMCalTriggerEventSelection.cxx
    Tracks/AliEMCalTriggerKineCuts.cxx
    Tracks/AliEMCalTriggerMCJetAnalysisComponent.cxx
    Tracks/AliEMCalTriggerMCParticleAnalysisComponent.cxx
    Tracks/AliEMCalTriggerPatchAnalysisComponent.cxx
    Tracks/AliEMCalTriggerRecJetAnalysisComponent.cxx
    Tracks/AliEMCalTriggerRecTrackAnalysisComponent.cxx
    Tracks/AliEMCalTriggerTaskGroup.cxx
    Tracks/AliEMCalTriggerTracksAnalysisComponent.cxx
    Tracks/AliEMCalTriggerWeightHandler.cxx
    Tracks/AliEMCalTriggerAnaTriggerClass.cxx
    Tracks/AliEMCalTriggerAnaClassManager.cxx
    Tracks/AliParticleMap.cxx
    Tracks/AliReducedHighPtEventCreator.cxx
    Tracks/AliReducedEmcalCluster.cxx
    Tracks/AliReducedPatchInfo.cxx
    Tracks/AliReducedPatchContainer.cxx
    Tracks/AliReducedMCHeader.cxx
    Tracks/AliReducedGeneratedParticle.cxx
    Tracks/AliReducedReconstructedTrack.cxx
    Tracks/AliReducedHighPtEvent.cxx
    Tracks/AliReducedJetEvent.cxx
    Tracks/AliReducedJetInfo.cxx
    Tracks/AliReducedJetParticle.cxx
    Tracks/AliReducedJetConstituent.cxx
    Tracks/AliReducedMatchedTrack.cxx
    Tracks/AliAnalysisTaskTrackDCA.cxx
    Tracks/AliEmcalAnalysisFactory.cxx
    Tracks/AliAnalysisTaskChargedParticlesRef.cxx
    Tracks/AliAnalysisTaskChargedParticlesRefMC.cxx
    Tracks/AliAnalysisTaskChargedParticlesMCTriggerMimic.cxx
    Tracks/AliAnalysisTaskEmcalTriggerMultiplicity.cxx
    Tracks/AliAnalysisTaskEmcalClustersRef.cxx
    Tracks/AliAnalysisTaskEmcalClusterMatched.cxx
    Tracks/AliAnalysisTaskEmcalPatchesRef.cxx
    Tracks/AliAnalysisTaskEmcalOnlinePatchesRef.cxx
    Tracks/AliAnalysisTaskEmcalOfflinePatchesRef.cxx
    Tracks/AliAnalysisTaskEmcalRecalcPatchesRef.cxx
    Tracks/AliAnalysisTaskEmcalTriggerCorrelation.cxx
    Tracks/AliAnalysisTaskEventSelectionRef.cxx
    Tracks/AliAnalysisTaskEventFilter.cxx
    Tracks/AliAnalysisTaskEtaPhiEfficiency.cxx
    Tracks/AliAnalysisTaskParticleInJet.cxx
    Tracks/AliAnalysisTaskTrackDensity.cxx
    Tracks/AliAnalysisTaskTrackDensityData.cxx
    Tracks/AliEmcalTriggerEmulation.cxx
    Tracks/AliAnalysisTaskEmcalTriggerEmulation.cxx
    Tracks/AliAnalysisTaskEGAMonitor.cxx
    Tracks/AliAnalysisTaskEmcalPatchMasking.cxx
    Tracks/AliAnalysisTaskCountITStracks.cxx
    Tracks/AliAnalysisTaskEmcalTriggerJets.cxx
    Tracks/AliAnalysisTaskEmcalTriggerJetsIDcorr.cxx
    Tracks/AliAnalysisTaskEmcalJetConstituentQA.cxx
    Tracks/AliAnalysisEmcalTriggerSelectionHelper.cxx
    Tracks/AliAnalysisTaskEmcalJetEnergyScale.cxx
    Tracks/AliAnalysisTaskEmcalJetEnergySpectrum.cxx
    Tracks/AliAnalysisTaskEmcalJetSpectrumSDPart.cxx
    Tracks/AliAnalysisTaskEmcalResponseOutliers.cxx
    Tracks/AliAnalysisTaskEmcalOutliersGen.cxx
    Tracks/AliAnalysisTaskEmcalPatchMultCorr.cxx
    Tracks/AliAnalysisTaskEmcalClustersInJets.cxx
    Tracks/AliAnalysisTaskEmcalTriggerCorrelationMC.cxx
    Tracks/AliAnalysisTaskEmcalHighEClusterTree.cxx
    Tracks/AliAnalysisTaskEmcalEG1Correlation.cxx
    Tracks/AliAnalysisTaskEmcalTriggerBackground.cxx
    Tracks/AliAnalysisTaskEmcalQoverPtShift.cxx
    Tracks/AliAnalysisTaskEmcalJetIterativeDeclustering.cxx
    Tracks/AliAnalysisTaskEmcalFastorMultiplicity.cxx
    Tracks/AliAnalysisTaskEmcalTriggerSelectionTest.cxx
    Tracks/AliEmcalTRDTrackingTask.cxx
    UserTasks/AliAnalysisTaskJetUEStudies.cxx
    UserTasks/AliAnalysisTaskBackFlucRandomCone.cxx
    UserTasks/AliAnalysisTaskChargedJetsPA.cxx
    UserTasks/AliAnalysisTaskChargedJetsHadronCF.cxx
    UserTasks/AliAnalysisTaskChargedJetsHadronToy.cxx
    UserTasks/Cumulants.cxx
    UserTasks/AliAnalysisTaskCLQA.cxx
    UserTasks/AliAnalysisTaskDcalDijetPerf.cxx
    UserTasks/AliAnalysisTaskDeltaPtJEmb.cxx
    UserTasks/AliAnalysisTaskDijetHadron.cxx
    UserTasks/AliAnalysisTaskParticleRandomizer.cxx
    UserTasks/AliAnalysisTaskEmcalBadCells.cxx
    UserTasks/AliAnalysisTaskEmcalDiJetAna.cxx
    UserTasks/AliAnalysisTaskEmcalDiJetBase.cxx
    UserTasks/AliAnalysisTaskEmcalDiJetResponse.cxx
    UserTasks/AliAnalysisTaskEmcalDijetImbalance.cxx
    UserTasks/AliAnalysisTaskEmcalVsPhos.cxx
    UserTasks/AliAnalysisTaskEmcalJetPerformance.cxx
    UserTasks/AliAnalysisTaskGenerateThermalBackgroundMC.cxx
    UserTasks/AliAnalysisTaskEmcalHighMultTrigger.cxx
    UserTasks/AliAnalysisTaskEmcalHJetMass.cxx
    UserTasks/AliAnalysisTaskIDFragmentationFunction.cxx
    UserTasks/AliAnalysisTaskEmcalJetCorrection.cxx
    UserTasks/AliAnalysisTaskEmcalJetHadCorQA.cxx
    UserTasks/AliAnalysisTaskEmcalJetHadEPpid.cxx
    UserTasks/AliAnalysisTaskEmcalJetHUtils.cxx
    UserTasks/AliAnalysisTaskEmcalJetHPerformance.cxx
    UserTasks/AliAnalysisTaskEmcalJetHCorrelations.cxx
    UserTasks/AliAnalysisTaskEmcalJetMassBkg.cxx
    UserTasks/AliAnalysisTaskEmcalJetMass.cxx
    UserTasks/AliAnalysisTaskEmcalJetMassResponse.cxx
    UserTasks/AliAnalysisTaskEmcalJetMassStructure.cxx
    UserTasks/AliAnalysisTaskEmcalJetPatchTriggerQA.cxx
    UserTasks/AliAnalysisTaskEmcalJetSparseMaker.cxx
    UserTasks/AliAnalysisTaskEmcalJetSpectra.cxx
    UserTasks/AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA.cxx
    UserTasks/AliAnalysisTaskEmcalJetSpectraMECpA.cxx
    UserTasks/AliAnalysisTaskEmcalJetTriggerMatcher.cxx
    UserTasks/AliAnalysisTaskEmcalJetTriggerQA.cxx
    UserTasks/AliAnalysisTaskEmcalJetv2QA.cxx
    UserTasks/AliAnalysisTaskEmcalRun2QA.cxx
    UserTasks/AliAnalysisTaskEmcalHfeTagging.cxx
    UserTasks/AliAnalysisTaskFakeJets.cxx
    UserTasks/AliAnalysisTaskEmcalMissingEnergy.cxx
    UserTasks/AliAnalysisTaskEmcalNeutralJets.cxx
    UserTasks/AliAnalysisTaskEmcalTriggerPatchJetMatch.cxx
    UserTasks/AliAnalysisTaskEmcalTmpSparseMaker.cxx
    UserTasks/AliAnalysisTaskEmcalTriggerInfoQA.cxx
    UserTasks/AliAnalysisTaskEmcalTriggerTreeWriter.cxx
    UserTasks/AliAnalysisTaskEventExtractor.cxx
    UserTasks/AliAnalysisTaskHJetEmbed.cxx
    UserTasks/AliAnalysisTaskJetJTJT.cxx
    UserTasks/AliAnalysisTaskJetMassResponseDet.cxx
    UserTasks/AliAnalysisTaskJetMatching.cxx
    UserTasks/AliAnalysisTaskJetsEvshape.cxx
    UserTasks/AliAnalysisTaskJetShapeBase.cxx
    UserTasks/AliAnalysisTaskJetShapeConst.cxx
    UserTasks/AliAnalysisTaskJetShapeDeriv.cxx
    UserTasks/AliAnalysisTaskJetShapeGR.cxx
    UserTasks/AliAnalysisTaskJetV2.cxx
    UserTasks/AliAnalysisTaskJetV3.cxx
    UserTasks/AliAnalysisTaskJetVn.cxx
    UserTasks/AliAnalysisTaskPID.cxx
    UserTasks/AliAnalysisTaskPIDV0base.cxx
    UserTasks/AliAnalysisTaskPSHFE.cxx
    UserTasks/AliAnalysisTaskPrepareInputForEmbedding.cxx
    UserTasks/AliAnalysisTaskRhoMassScale.cxx
    UserTasks/AliAnalysisTaskSOH.cxx

    UserTasks/AliAnalysisTaskTriggerRejection.cxx
    UserTasks/AliAnalysisTaskV0sInJetsEmcal.cxx
    UserTasks/AliEmcalPicoTrackFromJetMaker.cxx
    UserTasks/AliAnalysisTaskPi0Hadron.cxx
    UserTasks/AliAnalysisTaskEmcalJetTree.cxx
    UserTasks/AliJetEmbeddingSelRhoTask.cxx
    UserTasks/AliTrackContainerV0.cxx
    UserTasks/AliTrackContainerToyModel.cxx
    UserTasks/AliMCParticleContainerToyModel.cxx)

# Tasks that need direct access to RooUnfold
if(RooUnfold_FOUND)
    get_target_property(ROOUNFOLD_INCLUDE_DIR RooUnfold::RooUnfold INTERFACE_INCLUDE_DIRECTORIES)
    include_directories(${ROOUNFOLD_INCLUDE_DIR})
    add_definitions("-DWITH_ROOUNFOLD")
endif(RooUnfold_FOUND)

# Tasks that need direct access to Fastjet
if(FASTJET_FOUND)
    set(SRCS ${SRCS}
        AliFJWrapper.cxx
        AliEmcalJetUtility.cxx
        AliEmcalJetUtilityGenSubtractor.cxx
        AliEmcalJetUtilityConstSubtractor.cxx
	    AliEmcalJetUtilityEventSubtractor.cxx
        AliEmcalJetUtilitySoftDrop.cxx
        AliEmcalJetTask.cxx
        AliEmcalJetFinder.cxx
        AliJetEmbeddingFromAODTask.cxx
	    AliJetEmbeddingFromPYTHIATask.cxx
        AliJetShape.cxx
        AliLundPlaneHelper.cxx
      AliAnalysisTaskJetCharge.cxx
      AliAnalysisTaskJetChargeFlavourTemplates.cxx
	    UserTasks/AliAnalysisTaskEmcalQGTagging.cxx
	    UserTasks/AliAnalysisTaskNewJetSubstructure.cxx
            UserTasks/AliAnalysisTaskJetPull.cxx
	    UserTasks/AliAnalysisTaskCheckDeadCone.cxx
	    UserTasks/AliAnalysisTaskJetPull.cxx
	    UserTasks/AliAnalysisTaskLundPlane.cxx
	    UserTasks/AliAnalysisTaskJetDynamicalGrooming.cxx
	    UserTasks/AliAnalysisTaskJetHardestKt.cxx
	    UserTasks/AliAnalysisTaskHardestBranch.cxx
        UserTasks/AliAnalysisTaskEmcalJetShapesMC.cxx
        UserTasks/AliAnalysisTaskEmcalJetShapeExtra.cxx
	    UserTasks/AliAnalysisTaskFullpAJets.cxx
        UserTasks/AliAnalysisTaskFullppJet.cxx
        UserTasks/AliAnalysisTaskHJetDphi.cxx
        UserTasks/AliAnalysisTaskHJetSpectra.cxx
        UserTasks/AliAnalysisTaskJetExtractor.cxx
        UserTasks/AliAnalysisTaskJetPP.cxx
        UserTasks/AliAnalysisTaskSVtaskMCFilter.cxx
        UserTasks/AliAnalysisTaskEA.cxx
        UserTasks/AliAnalysisEmcalJetHelperEA.cxx
        UserTasks/AliAnalysisTaskPythiaBranchEA.cxx
        UserTasks/AliAnalysisTaskSubJetFraction.cxx
	    UserTasks/AliAnalysisTaskRecoilJetYield.cxx
	    UserTasks/AliAnalysisTaskRecursiveSoftDrop.cxx
        UserTasks/AliAnalysisTaskJetSubstructure.cxx
        UserTasks/AliAnalysisTaskSoftDrop.cxx
        UserTasks/AliAnalysisTaskSoftDropResponse.cxx
        UserTasks/AliAnalysisTaskSDKL.cxx
        UserTasks/AliAnalysisTaskSDKLResponse.cxx
        UserTasks/AliAnalysisTaskEmcalSubjet.cxx
	    UserTasks/AliAnalysisTaskJetCoreEmcal.cxx
        UserTasks/AliAnalysisTaskEmcalJetCDF.cxx
        Tracks/AliAnalysisTaskPtEfficiencyJets.cxx
        Tracks/AliHighPtReconstructionEfficiency.cxx
        Tracks/AliAnalysisTaskTracksInJet.cxx
	    Tracks/AliAnalysisTaskEmcalJetSubstructureTree.cxx
        Tracks/AliAnalysisEmcalSoftdropHelper.cxx
        Tracks/AliAnalysisTaskEmcalSoftDropData.cxx
        )
    # Tasks that need both Fastjet and RooUnfold
    if(RooUnfold_FOUND)
        set(SRCS ${SRCS}
        Tracks/AliAnalysisTaskEmcalSoftDropResponse.cxx)
    endif(RooUnfold_FOUND)
    include_directories(${AliPhysics_SOURCE_DIR}/JETAN/JETAN)
    include_directories(SYSTEM ${FASTJET_INCLUDE_DIR})
    link_directories(${FASTJET_LIBS_DIR})
    add_definitions(${FASTJET_DEFINITIONS})
    list(APPEND HDRS FJ_includes.h)
endif(FASTJET_FOUND)

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")
set(HDRS ${HDRS} Tracks/AliEMCalTriggerAnaHelper.h)

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}" "${FASTJET_ROOTDICT_OPTS}")

if(ROOT_VERSION_MAJOR EQUAL 6)
    set(ROOT_DEPENDENCIES PyROOT)
else()
    set(ROOT_DEPENDENCIES)
endif()

set(ALIROOT_DEPENDENCIES ANALYSIS CORRFW PHOSbase PWGEMCALbase PWGEMCALtrigger PWGJETFW PWGCaloTrackCorrBase PWGDevNanoAOD PWGPPevcharQn PWGPPevcharQnInterface KFParticle)
if(FASTJET_FOUND)
    set(ALIROOT_DEPENDENCIES JETAN ${ALIROOT_DEPENDENCIES})
endif(FASTJET_FOUND)

# Generate the ROOT map
# Dependecies
set(LIBDEPS ${ALIROOT_DEPENDENCIES} ${ROOT_DEPENDENCIES})
if(FASTJET_FOUND)
  set(LIBDEPS
    ${LIBDEPS}
    fastjetcontribfragile
    fastjetplugins
    siscone_spherical
    siscone
    fastjettools
    fastjet
    )
endif(FASTJET_FOUND)
set(ROOUNFOLD_DEPENDENCIES)
if(RooUnfold_FOUND)
    get_target_property(ROOUNFOLD_LIBRARY RooUnfold::RooUnfold IMPORTED_LOCATION)
    set(LIBDEPS ${LIBDEPS} ${ROOUNFOLD_LIBRARY})
    set(ROOUNFOLD_DEPENDENCIES ${ROOUNFOLD_LIBRARY})
endif(RooUnfold_FOUND)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Generate a PARfile target for this library. Note the extra includes ("Tracks UserTasks")
add_target_parfile(${MODULE} "${SRCS}" "${HDRS}" "${MODULE}LinkDef.h" "${LIBDEPS}" "Tracks UserTasks")

# Create an object to be reused in case of static libraries
# Otherwise the sources will be compiled twice
add_library(${MODULE}-object OBJECT ${SRCS} G__${MODULE}.cxx)
# Add a library to the project using the object
add_library_tested(${MODULE} SHARED $<TARGET_OBJECTS:${MODULE}-object>)
target_link_libraries(${MODULE} ${LIBDEPS})
target_link_libraries(${MODULE} \-L$ENV{KFPARTICLE_ROOT}/lib ${LIBDEPS})

# Setting the correct headers for the object as gathered from the dependencies
target_include_directories(${MODULE}-object PUBLIC $<TARGET_PROPERTY:${MODULE},INCLUDE_DIRECTORIES>)
set_target_properties(${MODULE}-object PROPERTIES COMPILE_DEFINITIONS $<TARGET_PROPERTY:${MODULE},COMPILE_DEFINITIONS>)

# Public include folders that will be propagated to the dependecies
target_include_directories(${MODULE} PUBLIC ${incdirs})

set(MODULE_COMPILE_FLAGS)
set(MODULE_LINK_FLAGS)

if(DATE_FOUND)
    set(MODULE_COMPILE_FLAGS "${DATE_CFLAGS}")
    set(MODULE_LINK_FLAGS "${DATE_LDFLAGS} ${DATE_LIBS}")
endif(DATE_FOUND)

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

# Installing the macros
install (DIRECTORY macros DESTINATION PWGJE/EMCALJetTasks)
