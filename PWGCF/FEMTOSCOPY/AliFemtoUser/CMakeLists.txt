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
set(MODULE PWGCFfemtoscopyUser)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliPhysics_SOURCE_DIR}/PWGCF/FEMTOSCOPY/AliFemtoUser)

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${AliPhysics_SOURCE_DIR}/OADB
                    ${AliPhysics_SOURCE_DIR}/PWGCF/FEMTOSCOPY/AliFemto
)

# Sources - alphabetical order
set(SRCS
  AliFemtoShareQualityPairCut.cxx
  AliFemtoShareQualityKTPairCut.cxx
  AliFemtoShareQualityTPCEntranceSepPairCut.cxx
  AliFemtoESDTrackCut.cxx
  AliFemtoEPSelectionTaskThird.cxx
  AliFemtoEventCutCentrality.cxx
  AliFemtoMJTrackCut.cxx
  AliFemtoKKTrackCut.cxx
  AliFemtoKKTrackCutFull.cxx
  AliFemtoKpmTrackCut.cxx
  AliFemtoKpm45TrackCut.cxx
  AliFemtoKtBinnedCorrFunc.cxx
  AliFemtoKTPairCutThird.cxx
  AliFemtoTPCInnerCorrFctn.cxx
  AliFemtoChi2CorrFctn.cxx
  AliFemtoModelAllHiddenInfo.cxx
  AliFemtoModelBPLCMSCorrFctn.cxx
  AliFemtoModelBPLCMSCorrFctnKK.cxx
  AliFemtoModelGausRinvFreezeOutGenerator.cxx
  AliFemtoModelCorrFctn3DSpherical.cxx
#  Gleb:
   AliFemtoModelCorrFctn3DKKGR.cxx
   AliFemtoModelCorrFctn3DKKGRQinv.cxx
   AliFemtoModelGausLCMSFreezeOutGeneratorGR.cxx
#  Ludmila
   AliFemtoCutMonitorDphiDeta1.cxx
   AliFemtoESDTrackCutPlusJets25.cxx
#
  AliFemtoModelCorrFctn3DLCMSSpherical.cxx
  AliFemtoModelCorrFctnDEtaDPhi.cxx
  AliFemtoModelCorrFctnDEtaDPhiStar.cxx
  AliFemtoCorrFctnDEtaDPhiWithWeights.cxx
  AliFemtoQinvCorrFctnWithWeights.cxx
  AliFemtoModelCorrFctnDEtaDPhiWithWeights.cxx
  AliFemtoModelCorrFctnDirectYlm.cxx
  AliFemtoModelCorrFctnSource.cxx
  AliFemtoModelCorrFctnKStar.cxx
  AliFemtoModelCorrFctnNonIdDR.cxx
  AliFemtoModelCorrFctnQinv.cxx
  AliFemtoModelCorrFctnQinvPtSum.cxx
  AliFemtoModelCorrFctnTrueQ.cxx
  AliFemtoModelCorrFctnTrueQ3D.cxx
  AliFemtoModelCorrFctnTrueQ3DByParent.cxx
  AliFemtoModelCorrFctnTrueQ6D.cxx
  AliFemtoModelCorrFctnKK.cxx
  AliFemtoModelCorrFctnKpKp.cxx
  AliFemtoNSigmaFilter.cxx
  AliFemtoCutAttrEvent.cxx
  AliFemtoCutAttrPairTrack.cxx
  AliFemtoCutAttrTrack.cxx
  AliFemtoAnalysisPionPionCuts.cxx
  AliFemtoAnalysisPionPionObjectConstructor.cxx
  AliFemtoCutMonitorParticleEtCorr.cxx
  AliFemtoCutMonitorParticleMomRes.cxx
  AliFemtoCutMonitorParticlePtPDG.cxx
  AliFemtoCutMonitorParticlePtPDGV0.cxx
  AliFemtoCutMonitorParticleNumber.cxx
  AliFemtoCutMonitorEventNumber.cxx
  AliFemtoPionLambdaCutMonitor.cxx
  AliFemtoCutMonitorPionPion.cxx
  AliFemtoCutMonitorTrackTPCchiNdof.cxx
  AliFemtoCutMonitorTrackTPCncls.cxx
  AliFemtoQATrackCut.cxx
  AliFemtoQAEventCut.cxx
  AliFemtoShareQualityTPCEntranceSepQAPairCut.cxx
  AliFemtoShareQualityQAPairCut.cxx
  AliFemtoShareQualityCorrFctn.cxx
  AliFemtoYlm.cxx
  AliFemtoCorrFctn3DLCMSPosQuad.cxx
  AliFemtoCorrFctnDirectYlm.cxx
  AliFemtoCorrFctnTPCNcls.cxx
  AliFemtoCorrFctnDEtaDPhi.cxx
  AliFemtoCorrFctnDEtaDPhiStar.cxx
  AliFemtoCorrFctnDEtaDPhiTHn.cxx
  AliFemtoCorrFctnDEtaDPhiCorrections.cxx
  AliFemtoCorrFctnDYDPhi.cxx
  AliFemtoCorrFctnKStar.cxx
  AliFemtoCorrFctnGammaMonitor.cxx
  AliFemtoCorrFctnQLCMS.cxx
  AliFemtoCorrFctn3DSphericalEMCIC.cxx
  AliFemtoCorrFctnPtKstar.cxx
  AliFemtoQinvCorrFctnEMCIC.cxx
  AliFemtoBPLCMS3DCorrFctnEMCIC.cxx
  AliFemtoPairCutAntiGamma.cxx
  AliFemtoPairCutRadialDistance.cxx
  AliFemtoPairCutPt.cxx
  AliFemtoAvgSepCorrFctn.cxx
  AliFemtoAnalysisAzimuthalPbPb.cxx
  AliFemtoAnalysisAzimuthalPbPbThird.cxx
  AliFemtoAnalysisAzimuthalPbPb2Order.cxx
  AliFemtoAnalysisPionLambda.cxx
  AliFemtoAnalysisPionPion.cxx
  AliFemtoAnalysisPionPionObjectConstructor.cxx
  AliFemtoPairCutRadialDistanceLM.cxx
  AliFemtoPicoEventRP.cxx
  AliFemtoPairCutRadialDistanceKK.cxx
  AliFemtoPairCutResonances.cxx
  AliFemtoCorrFctnMinvMonitor.cxx
  AliFemtoCorrFctnPairFractions.cxx
  AliAnalysisTaskParticleEfficiency.cxx
  AliAnalysisTaskParticleFractions.cxx
  AliFemtoPairOriginMonitor.cxx
  AliFemtoPairCutAntiGammaAlpha.cxx
  AliFemtoCorrFctnGammaMonitorAlpha.cxx
  AliFemtoCorrFctn3DPRF.cxx
  AliFemtoCorrFctn3DPRF_qosl_q.cxx
  AliFemtoPairCutRadialDistanceAsymmetric.cxx
  AliFemtoPairCutGeneralizedRadialDistance.cxx
  AliFemtoPairCutMergedFraction.cxx
  AliFemtoCorrFctnDEtaDPhiSimple.cxx
  AliFemtoCorrFctnDEtaDPhiSimpleWithCorrections.cxx
  AliFemtoCorrFctnDYDPhiSimple.cxx
  AliFemtoV0TrackCutNSigmaFilter.cxx
  AliFemtoXiTrackCutNSigmaFilter.cxx
  AliFemtoCutMonitorEventPartCollSize.cxx
  AliFemtoModelCorrFctnWithWeights.cxx
  AliFemtoCutMonitorParticleYPtWithWeights.cxx
  Faddeeva.cxx
  AliFemtoESDTrackCutNSigmaFilter.cxx
  AliFemtoModelCorrFctnKStarFull.cxx
  AliFemtoModelWeightGeneratorBasicLednicky.cxx
  AliFemtoAnalysisLambdaKaon.cxx
  AliFemtoModelCorrFctnDEtaDPhiAK.cxx
  AliFemtoModelCorrFctnDEtaDPhiRM.cxx
  AliFemtoCutMonitorPairKT.cxx
  AliFemtoPairCutMInv.cxx
  AliFemtoV0PurityBgdEstimator.cxx
  AliFemtoPairCutRadialDistanceKKdist.cxx
  AliFemtoSphericityEventCutKK.cxx #  Sphericity for KK analysis in pp13TeV
  AliFemtoSphericityEventCutKK1.cxx #  3may2019!! Sphericity for KK analysis in pp13TeV
  AliFemtoCutMonitorEventSphericity.cxx # Sphericity monitor
  AliFemtoKKTrackCutTest.cxx # new class
  AliFemtoCorrFctnInvMass.cxx
  AliFemtoCutMonitorV0CosPointingAngle.cxx
  AliFemtoCutMonitorPairMomRes.cxx
  AliFemtoCorrFctnPairsForCorrFit.cxx
  AliFemtoWRzTrackCut.cxx
  )



# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

# Additional headers
set(HDRS ${HDRS}
  AliFemtoPairCutDetaDphi.h
  AliFemtoCorrFctnQ3D.h
  )

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Add a shared library
add_library_tested(${MODULE} SHARED  ${SRCS} G__${MODULE}.cxx)

# Generate the ROOT map
# Dependecies
set(LIBDEPS ANALYSISalice PWGCFfemtoscopy)
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
