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
set(MODULE PWGHFhfe)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWGHF/hfe)

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    $ENV{KFPARTICLE_ROOT}/include
                    ${AliPhysics_SOURCE_DIR}/CORRFW
                    ${AliPhysics_SOURCE_DIR}/OADB
                    ${AliPhysics_SOURCE_DIR}/PWG/CaloTrackCorrBase
                    ${AliPhysics_SOURCE_DIR}/PWG/EMCAL/Trigger
                    ${AliPhysics_SOURCE_DIR}/PWG/FLOW/Base
                    ${AliPhysics_SOURCE_DIR}/PWG/FLOW/Tasks
                    ${AliPhysics_SOURCE_DIR}/PWG/muon
                    ${AliPhysics_SOURCE_DIR}/PWG/TRD
                    ${AliPhysics_SOURCE_DIR}/PWGPP/EVCHAR/FlowVectorCorrections/QnCorrections
                    ${AliPhysics_SOURCE_DIR}/PWGPP/EVCHAR/FlowVectorCorrections/QnCorrectionsInterface
		    ${AliPhysics_SOURCE_DIR}/PWGHF/vertexingHF
		    ${AliPhysics_SOURCE_DIR}/PWGHF/vertexingHF/charmFlow)

# Sources - alphabetical order
set(SRCS
  AliHFEtools.cxx 
  AliHFEparamBag.cxx 
  AliHFEcollection.cxx 
  AliHFEvarManager.cxx 
  AliHFEcontainer.cxx 
  AliHFEmcQA.cxx 
  AliHFEpairs.cxx 
  AliHFEsecVtxs.cxx 
  AliHFEsecVtx.cxx 
  AliHFEpriVtx.cxx 
  AliHFEspectrum.cxx 
  AliHFEV0info.cxx 
  AliHFEV0pid.cxx 
  AliHFEV0taginfo.cxx
  AliHFEpidQA.cxx 
  AliHFEtrdPIDqa.cxx 
  AliHFEitsPIDqa.cxx
  AliHFEdetPIDqa.cxx 
  AliHFEtpcPIDqa.cxx 
  AliHFEtofPIDqa.cxx 
  AliHFEemcalPIDqa.cxx 
  AliHFEtrdPIDqaV1.cxx 
  AliHFEpidQAmanager.cxx 
  AliHFEpid.cxx 
  AliHFEpidBase.cxx 
  AliHFEpidITS.cxx 
  AliHFEpidTPC.cxx 
  AliHFEpidTRD.cxx 
  AliHFEpidTOF.cxx 
  AliHFEpidEMCAL.cxx 
  AliHFEpidMC.cxx 
  AliHFEextraCuts.cxx 
  AliHFEextraEventCuts.cxx 
  AliHFEcutStep.cxx 
  AliHFEtrackFilter.cxx 
  AliHFEcuts.cxx 
  AliHFEsignalCuts.cxx 
  AliHFEelecbackground.cxx 
  AliHFEtaggedTrackAnalysis.cxx 
  AliHFEpostAnalysis.cxx 
  AliAnalysisTaskHFEpidQA.cxx 
  AliHFEefficiency.cxx 
  AliAnalysisTaskHFE.cxx 
  AliHFEV0pidMC.cxx 
  AliHFEV0cuts.cxx
  AliHFEOADBThresholdsTRD.cxx
  AliHFEbayesPIDqa.cxx
  AliHFEpidBayes.cxx
  AliHFEpidObject.cxx
  AliAnalysisTaskEHCorrel.cxx
  AliAnalysisTaskFlowTPCEMCalEP.cxx
  AliAnalysisTaskFlowTPCEMCalRun2.cxx
  AliAnalysisTaskHFECal.cxx 
  AliAnalysisTaskEMCalHFEpA.cxx 
  AliAnalysisTaskHFEpACorrelation.cxx
  AliHFEdebugTreeTask.cxx
  AliHFEVZEROEventPlane.cxx
  AliAnalysisTaskFlowTPCTOFEPSP.cxx
  AliSelectNonHFE.cxx
  AliHFENonPhotonicElectron.cxx
  AliHFEdebugTreeTaskAOD.cxx
  AliHFECorrectSpectrumBase.cxx
  AliHFEInclusiveSpectrum.cxx
  AliHFEInclusiveSpectrumQA.cxx 
  AliHFEBeautySpectrum.cxx
  AliHFEBeautySpectrumQA.cxx 
  AliHFEsmearDCA.cxx
  AliAnalysisTaskFlowTPCEMCalQCSP.cxx
  AliHFEreducedEventCreatorAOD.cxx
  AliHFEreducedEventCreatorESD.cxx
  AliHFEreducedEvent.cxx
  AliHFEreducedTrack.cxx
  AliHFEreducedMCParticle.cxx
  AliHFEminiTrack.cxx
  AliHFEminiEvent.cxx
  AliHFEminiEventCreator.cxx
  AliAnalysisTaskHFEQA.cxx
  AliAnalysisTaskHFEemcQA.cxx
  AliAnalysisTaskBeautyCal.cxx
  AliAnalysisTaskTPCCalBeauty.cxx
  AliAnalysisTaskFlowITSTPCTOFQCSP.cxx
  AliAnalysisTaskHFEIPDistribution.cxx
  AliAnalysisTaskHFEBeautyMCTemplates.cxx
  AliAnalysisTaskHFEBeautyMCTemplatesRun2.cxx
  AliAnalysisTaskHFEIPCorrection.cxx
  AliAnalysisTaskHFEMulti.cxx
  AliAnalysisTaskHFEMultiplicity.cxx
  AliAnalysisTaskHFEtemplate.cxx
  AliHFEhadronicbackground.cxx
  AliHFEUnfolding.cxx	
  AliAnalysisTaskHFEEfficiency.cxx
  AliAnalysisTaskHaHFECorrel.cxx
  AliAnalysisHFETPCTOFNew.cxx
  AliAnalysisHFETPCTOFBeauty.cxx
  AliAnalysisHFEppTPCTOFBeauty.cxx
  AliAnalysisHFEppTPCTOFBeauty5TeVNew.cxx	
  AliAnalysisTaskCaloHFEpPbRun2.cxx
  AliAnalysisTaskHFEmultTPCTOF.cxx
  AliAnalysisTaskCaloHFEpp.cxx
  AliAnalysisTaskHFEBESpectraEMC.cxx
  AliAnalysisTaskHFETPCTOFMultiplicity.cxx
  AliAnalysisTaskQAHFE.cxx
  AliAnalysisTask_QA_EMCALElectrons.cxx
  AliAnalysisTaskHFEmultTPCEMCAL.cxx
  AliAnalysisTaskHFEBeautyMultiplicity.cxx 
  AliAnalysisHFEppEMCalBeauty.cxx
  AliAnalysisTaskIPResolBeautyppCal.cxx
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
set(LIBDEPS OADB ANALYSISalice CORRFW PWGflowTasks PWGTRD MLP PWGPPevcharQn PWGPPevcharQnInterface PWGHFvertexingHF KFParticle)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Generate a PARfile target for this library
add_target_parfile(${MODULE} "${SRCS}" "${HDRS}" "${MODULE}LinkDef.h" "${LIBDEPS}")

# Linking the library
target_link_libraries(${MODULE} ${LIBDEPS})
target_link_libraries(${MODULE} \-L$ENV{KFPARTICLE_ROOT}/lib ${LIBDEPS})

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
install(DIRECTORY macros DESTINATION PWGHF/hfe)
