/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Friederike Bock 									              *
 * Version 1                                                              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//

#include <memory>
#include <vector>
#include "TParticle.h"
#include "TPDGCode.h"
#include "TMCProcess.h"
#include "TDatabasePDG.h"
#include "TList.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1F.h"
#include "THnSparse.h"
#include "TH2F.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliPID.h"
#include "AliLog.h"
#include "AliESDtrackCuts.h"
#include "AliESDpidCuts.h"
#include "AliMCEvent.h"
#include "AliESDv0.h"
#include "AliESDEvent.h"
#include "AliESDpid.h"
#include "AliKFParticle.h"
#include "AliMCEventHandler.h"
#include "AliKFVertex.h"
#include "AliTriggerAnalysis.h"
#include "AliCentrality.h"
#include "AliMultiplicity.h"
#include "AliAODEvent.h"
#include "AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson.h"
#include "AliCaloTrackMatcher.h"
#include <vector>

ClassImp( AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson )

//-----------------------------------------------------------------------------------------------
AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson():
  fV0Reader(nullptr),
  fV0ReaderName("V0ReaderV1"),
  fPionSelector(NULL),
  fPionSelectorName("PionSelector"),
  fBGHandlerPiPl(NULL),
  fBGHandlerPiMi(NULL),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fCutFolder(NULL),
  fESDList(NULL),
  fTrueList(NULL),
  fTrueTreeList(NULL),
  fMCList(NULL),
  fOutputContainer(0),
  fReaderGammas(nullptr),
  fSelectorNegPionIndex(0),
  fSelectorPosPionIndex(0),
  fGoodConvGammas(nullptr),
  fClusterCandidates(nullptr),
  fNeutralDecayParticleCandidates(nullptr),
  fNeutralDecayParticleSidebandCandidates(nullptr),
  fNeutralDecayParticleSwappCandidates(nullptr),
  fPosPionCandidates(nullptr),
  fNegPionCandidates(nullptr),
  fEventCutArray(nullptr),
  fGammaCutArray(nullptr),
  fClusterCutArray(nullptr),
  fPionCutArray(nullptr),
  fNeutralDecayMesonCutArray(nullptr),
  fMesonCutArray(nullptr),
  fEventCuts(nullptr),
  fConversionCuts(nullptr),
  fClusterCuts(nullptr),
  fOutlierJetReader(nullptr),
  fTreePiPiSameMother(nullptr),
  fTreePiPiPiSameMother(nullptr),
  fTreeEventInfoHNM(nullptr),
  fCasePiPi(-1),
  fSamePiPiMotherID(-1),
  fSamePiPiMotherInvMass(-1),
  fSamePiPiMotherPt(-1),
  fSamePiPiPiMotherID(-1),
  fSamePiPiPiMotherInvMass(-1),
  fSamePiPiPiMotherPt(-1),
  fV0MultiplicityHNMEvent(-1),
  fTrackMultiplicityHNMEvent(-1),
  fZVertexHNMEvent(-1),
  fPtHNM(-1),
  fPDGMassNDM(-1),
  fNDMMinPtPossible(0.),
  fPDGMassChargedPion(-1),
  fPDGCodeNDM(-1),
  fPDGCodeAnalyzedMeson(-1),
  enableDalitzAllPt(kFALSE),
  enableDalitzLowPt(kFALSE),
  enableDalitzMidPt(kFALSE),
  enableDalitzHighPt(kFALSE),
  HistoDalitzPtRangeMin_LowPt(0.),
  HistoDalitzPtRangeMax_LowPt(0.),
  HistoDalitzPtRangeMin_MidPt(0.),
  HistoDalitzPtRangeMax_MidPt(0.),
  HistoDalitzPtRangeMin_HighPt(0.),
  HistoDalitzPtRangeMax_HighPt(0.),
  fHistoConvGammaPt(nullptr),
  fHistoConvGammaEta(nullptr),
  fHistoClusterGammaPt(nullptr),
  fHistoClusterGammaEta(nullptr),
  fHistoClusterGammaE(nullptr),
  fHistoNumberClusterGamma(nullptr),
  fHistoNegPionPt(nullptr),
  fHistoPosPionPt(nullptr),
  fHistoNegPionPhi(nullptr),
  fHistoPosPionPhi(nullptr),
  fHistoNegPionEta(nullptr),
  fHistoPosPionEta(nullptr),
  fHistoNegPionClsTPC(nullptr),
  fHistoPosPionClsTPC(nullptr),
  fHistoPionDCAxy(nullptr),
  fHistoPionDCAz(nullptr),
  fHistoPionTPCdEdxNSigma(nullptr),
  fHistoPionTPCdEdx(nullptr),
  fHistoPionPionInvMassPt(nullptr),
  fHistoGammaGammaInvMassPt(nullptr),
  fHistoGammaGammaInvMassPtBeforeCuts(nullptr),
  fHistoSwappingGammaGammaInvMassPt(nullptr),
  fHistoMotherInvMassPt(nullptr),
  fHistoMotherInvMassPtRejectedKinematic(nullptr),
  fHistoDalitzPlotPosFixedPzNDM(nullptr),
  fHistoDalitzPlotNegFixedPzNDM(nullptr),
  fHistoDalitzPlotPosSubNDM(nullptr),
  fHistoDalitzPlotNegSubNDM(nullptr),
  fHistoDalitzPlotPosFixedPzNDM_LowPt(nullptr),
  fHistoDalitzPlotNegFixedPzNDM_LowPt(nullptr),
  fHistoDalitzPlotPosSubNDM_LowPt(nullptr),
  fHistoDalitzPlotNegSubNDM_LowPt(nullptr),
  fHistoDalitzPlotPosFixedPzNDM_MidPt(nullptr),
  fHistoDalitzPlotNegFixedPzNDM_MidPt(nullptr),
  fHistoDalitzPlotPosSubNDM_MidPt(nullptr),
  fHistoDalitzPlotNegSubNDM_MidPt(nullptr),
  fHistoDalitzPlotPosFixedPzNDM_HighPt(nullptr),
  fHistoDalitzPlotNegFixedPzNDM_HighPt(nullptr),
  fHistoDalitzPlotPosSubNDM_HighPt(nullptr),
  fHistoDalitzPlotNegSubNDM_HighPt(nullptr),
  fHistoBackInvMassPt(nullptr),
  fHistoMotherLikeSignBackInvMassPt(nullptr),
  fHistoAngleHNMesonPiPlPiMi(nullptr),
  fHistoAngleHNMesonNDM(nullptr),
  fHistoAngleHNMesonPiPl(nullptr),
  fHistoAngleHNMesonPiMi(nullptr),
  fHistoAnglePiPlPiMi(nullptr),
  fHistoAngleNDMPiMi(nullptr),
  fHistoAnglePiPlNDM(nullptr),
  fHistoAngleSum(nullptr),
  fHistoTrueAngleSum(nullptr),
  fHistoTrueHNMesonPtvsNDMPt(nullptr),
  fHistoMotherInvMassSubNDM(nullptr),
  fHistoBackInvMassPtSubNDM(nullptr),
  fHistoMotherLikeSignBackInvMassSubNDMPt(nullptr),
  fHistoMotherInvMassFixedPzNDM(nullptr),
  fHistoBackInvMassPtFixedPzNDM(nullptr),
  fHistoMotherLikeSignBackInvMassFixedPzNDMPt(nullptr),
  fHistoMCAllGammaPt(nullptr),
  fHistoMCConvGammaPt(nullptr),
  fHistoMCAllMesonPt(nullptr),
  fHistoMCAllMesonEta(nullptr),
  fHistoMCAllMesonPhi(nullptr),
  fHistoMCMesonFromNeutralMesonPt(nullptr),
  fHistoMCMesonFromNeutralMesonEta(nullptr),
  fHistoMCMesonFromNeutralMesonPhi(nullptr),
  fHistoMCAllPosPionsPt(nullptr),
  fHistoMCAllPosPionsEta(nullptr),
  fHistoMCAllPosPionsPhi(nullptr),
  fHistoMCAllNegPionsPt(nullptr),
  fHistoMCAllNegPionsEta(nullptr),
  fHistoMCAllNegPionsPhi(nullptr),
  fHistoMCGammaFromNeutralMesonPt(nullptr),
  fHistoMCPosPionsFromNeutralMesonPt(nullptr),
  fHistoMCPosPionsFromNeutralMesonEta(nullptr),
  fHistoMCPosPionsFromNeutralMesonPhi(nullptr),
  fHistoMCNegPionsFromNeutralMesonPt(nullptr),
  fHistoMCNegPionsFromNeutralMesonEta(nullptr),
  fHistoMCNegPionsFromNeutralMesonPhi(nullptr),
  fHistoMCHNMPiPlPiMiNDMPt(nullptr),
  fHistoMCHNMPiPlPiMiNDMEta(nullptr),
  fHistoMCHNMPiPlPiMiNDMPhi(nullptr),
  fHistoMCHNMPiPlPiMiNDMEtavsPt(nullptr),
  fHistoMCHNMPiPlPiMiNDMInAccPt(nullptr),
  fHistoMCHNMInAccVsNDMPt(nullptr),
  fHistoMCHeavyAllPt(nullptr),
  fHistoMCHeavyAllEta(nullptr),
  fHistoMCHeavyAllPhi(nullptr),
  fHistoMCHeavyChannelPt(nullptr),
  fHistoMCHeavyChannelEta(nullptr),
  fHistoMCHeavyChannelPhi(nullptr),
  fHistMCChannelNDMFromHeavyPt(nullptr),
  fHistMCChannelNDMFromHeavyEta(nullptr),
  fHistMCChannelNDMFromHeavyPhi(nullptr),
  fHistMCChannelPiPlusFromHeavyPt(nullptr),
  fHistMCChannelPiPlusFromHeavyEta(nullptr),
  fHistMCChannelPiPlusFromHeavyPhi(nullptr),
  fHistMCChannelPiMinusFromHeavyPt(nullptr),
  fHistMCChannelPiMinusFromHeavyEta(nullptr),
  fHistMCChannelPiPMinusFromHeavyPhi(nullptr),
  fHistMCChannelNDMPtHeavyPt(nullptr),
  fHistMCChannelPiPlusPtHeavyPt(nullptr),
  fHistMCChannelPiMinusPtHeavyPt(nullptr),
  fHistoMCHeavyReconstructiblePt(nullptr),
  fHistoMCHeavyReconstructibleEta(nullptr),
  fHistoMCHeavyReconstructiblePhi(nullptr),
  fHistMCReconstructibleNDMFromHeavyPt(nullptr),
  fHistMCReconstructibleNDMFromHeavyEta(nullptr),
  fHistMCReconstructibleNDMFromHeavyPhi(nullptr),
  fHistMCReconstructiblePiPlusFromHeavyPt(nullptr),
  fHistMCReconstructiblePiPlusFromHeavyEta(nullptr),
  fHistMCReconstructiblePiPlusFromHeavyPhi(nullptr),
  fHistMCReconstructiblePiMinusFromHeavyPt(nullptr),
  fHistMCReconstructiblePiMinusFromHeavyEta(nullptr),
  fHistMCReconstructiblePiPMinusFromHeavyPhi(nullptr),
  fHistMCReconstructibleNDMPtHeavyPt(nullptr),
  fHistMCReconstructiblePiPlusPtHeavyPt(nullptr),
  fHistMCReconstructiblePiMinusPtHeavyPt(nullptr),
  fHistoTrueMesonFlags(nullptr),
  fHistoTrueMotherPiPlPiMiNDMInvMassPt(nullptr),
  fHistoTrueMotherPiPlPiMiNDMInvMassPtSubNDM(nullptr),
  fHistoTrueMotherPiPlPiMiNDMInvMassPtFixedPzNDM(nullptr),
  fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromDifferent(nullptr),
  fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromEtaOmega(nullptr),
  fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromRho(nullptr),
  fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromK0s(nullptr),
  fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromK0l(nullptr),
  fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromEtaPrime(nullptr),
  fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromOther(nullptr),
  fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM(nullptr),
  fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM(nullptr),
  fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM(nullptr),
  fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM(nullptr),
  fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM_LowPt(nullptr),
  fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM_LowPt(nullptr),
  fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM_LowPt(nullptr),
  fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM_LowPt(nullptr),
  fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM_MidPt(nullptr),
  fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM_MidPt(nullptr),
  fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM_MidPt(nullptr),
  fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM_MidPt(nullptr),
  fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM_HighPt(nullptr),
  fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM_HighPt(nullptr),
  fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM_HighPt(nullptr),
  fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM_HighPt(nullptr),
  fHistoTrueMotherGammaGammaInvMassPt(nullptr),
  fHistoTrueMotherGammaGammaFromHNMInvMassPt(nullptr),
  fHistoTrueConvGammaPt(nullptr),
  fHistoTrueConvGammaFromNeutralMesonPt(nullptr),
  fHistoTrueClusterGammaPt(nullptr),
  fHistoTrueClusterGammaFromNeutralMesonPt(nullptr),
  fHistoTruePosPionPt(nullptr),
  fHistoTruePosPionFromNeutralMesonPt(nullptr),
  fHistoTrueNegPionPt(nullptr),
  fHistoTrueNegPionFromNeutralMesonPt(nullptr),
  fHistoTruePionPionInvMassPt(nullptr),
  fHistoTruePionPionFromSameMotherInvMassPt(nullptr),
  fHistoTruePionPionFromHNMInvMassPt(nullptr),
  fHistoTruePiPlPiMiSameMotherFromEtaInvMassPt(nullptr),
  fHistoTruePiPlPiMiSameMotherFromOmegaInvMassPt(nullptr),
  fHistoTruePiPlPiMiSameMotherFromRhoInvMassPt(nullptr),
  fHistoTruePiPlPiMiSameMotherFromEtaPrimeInvMassPt(nullptr),
  fHistoTruePiPlPiMiSameMotherFromK0sInvMassPt(nullptr),
  fHistoTruePiPlPiMiSameMotherFromK0lInvMassPt(nullptr),
  fHistoTruePiPlPiMiSameMotherFromOtherlInvMassPt(nullptr),
  fHistoTruePiMiPiZeroSameMotherFromEtaInvMassPt(nullptr),
  fHistoTruePiMiPiZeroSameMotherFromOmegaInvMassPt(nullptr),
  fHistoTruePiMiPiZeroSameMotherFromRhoInvMassPt(nullptr),
  fHistoTruePiMiPiZeroSameMotherFromK0lInvMassPt(nullptr),
  fHistoTruePiMiPiZeroSameMotherFromOtherlInvMassPt(nullptr),
  fHistoTruePiPlPiZeroSameMotherFromEtaInvMassPt(nullptr),
  fHistoTruePiPlPiZeroSameMotherFromOmegaInvMassPt(nullptr),
  fHistoTruePiPlPiZeroSameMotherFromRhoInvMassPt(nullptr),
  fHistoTruePiPlPiZeroSameMotherFromK0lInvMassPt(nullptr),
  fHistoTruePiPlPiZeroSameMotherFromOtherInvMassPt(nullptr),
  fHistoTruePiPlPiMiNDMPureCombinatoricalInvMassPt(nullptr),
  fHistoTruePiPlPiMiNDMCombinatoricalInvMassPt(nullptr),
  fHistoTruePiPlPiMiNDMContaminationInvMassPt(nullptr),
  fHistoTruePiPlPiMiNDMContamination_Pi0_InvMassPt(nullptr),
  fHistoTruePiPlPiMiNDMContamination_PiPl_InvMassPt(nullptr),
  fHistoTruePiPlPiMiNDMContamination_PiMi_InvMassPt(nullptr),
  fHistoTruePiPlPiMiNDMContamination_Crosscheck_InvMassPt(nullptr),
  fHistoTruePiPlPiMiNDMContamination_multipel_InvMassPt(nullptr),
  fHistoDoubleCountTruePi0InvMassPt(nullptr),
  fHistoDoubleCountTrueHNMInvMassPt(nullptr),
  fHistoDoubleCountTrueConvGammaRPt(nullptr),
  fVectorDoubleCountTruePi0s(0),
  fVectorDoubleCountTrueHNMs(0),
  fVectorDoubleCountTrueConvGammas(0),
  fHistoNEvents(nullptr),
  fHistoNEventsWOWeight(nullptr),
  fProfileJetJetXSection(nullptr),
  fHistoJetJetNTrials(nullptr),
  fHistoNGoodESDTracks(nullptr),
  fProfileEtaShift(nullptr),
  fHistoSPDClusterTrackletBackground(nullptr),
  fHistovParticleChi2PerNDF(nullptr),
  fHistovParticleChi2PerNDFBothConstrained(nullptr),
  fHistovParticleChi2PerNDFOneConstrained(nullptr),
  fHistovParticledS(nullptr),
  fHistovParticledSBothConstrained(nullptr),
  fHistovParticledSOneConstrained(nullptr),
  fHistoTruevParticleChi2PerNDF(nullptr),
  fHistoTruevParticleFromSameMotherChi2PerNDF(nullptr),
  fHistoTruevParticleFromHNMChi2PerNDF(nullptr),
  fHistoTruevParticledS(nullptr),
  fHistoTruevParticleFromSameMotherdS(nullptr),
  fHistoTruevParticleFromHNMdS(nullptr),
  fRandom(0),
  fnCuts(0),
  fiCut(0),
  fNumberOfESDTracks(0),
  fMoveParticleAccordingToVertex(kFALSE),
  fIsHeavyIon(0),
  fDoMesonAnalysis(kTRUE),
  fDoMesonQA(0),
  fIsFromMBHeader(kTRUE),
  fIsFromDesiredHeader(kTRUE),
  fIsOverlappingWithOtherHeader(kFALSE),
  fAllowOverlapHeaders(kTRUE),
  fIsMC(kFALSE),
  fSelectedHeavyNeutralMeson(kFALSE),
  fDoLightOutput(kFALSE),
  fNDMRecoMode(0),
  fTolerance(-1),
  fWeightJetJetMC(1.),
  fTrackMatcherRunningMode(0),
  fMCEventPos(),
  fMCEventNeg(),
  fESDArrayPos(),
  fESDArrayNeg()
{

}

//-----------------------------------------------------------------------------------------------
AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson( const char* name ):
  AliAnalysisTaskSE(name),
  fV0Reader(nullptr),
  fV0ReaderName("V0ReaderV1"),
  fPionSelector(NULL),
  fPionSelectorName("PionSelector"),
  fBGHandlerPiPl(NULL),
  fBGHandlerPiMi(NULL),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fCutFolder(NULL),
  fESDList(NULL),
  fTrueList(NULL),
  fTrueTreeList(NULL),
  fMCList(NULL),
  fOutputContainer(0),
  fReaderGammas(nullptr),
  fSelectorNegPionIndex(0),
  fSelectorPosPionIndex(0),
  fGoodConvGammas(nullptr),
  fClusterCandidates(nullptr),
  fNeutralDecayParticleCandidates(nullptr),
  fNeutralDecayParticleSidebandCandidates(nullptr),
  fNeutralDecayParticleSwappCandidates(nullptr),
  fPosPionCandidates(nullptr),
  fNegPionCandidates(nullptr),
  fEventCutArray(nullptr),
  fGammaCutArray(nullptr),
  fClusterCutArray(nullptr),
  fPionCutArray(nullptr),
  fNeutralDecayMesonCutArray(nullptr),
  fMesonCutArray(nullptr),
  fEventCuts(nullptr),
  fConversionCuts(nullptr),
  fClusterCuts(nullptr),
  fOutlierJetReader(nullptr),
  fTreePiPiSameMother(nullptr),
  fTreePiPiPiSameMother(nullptr),
  fTreeEventInfoHNM(nullptr),
  fCasePiPi(-1),
  fSamePiPiMotherID(-1),
  fSamePiPiMotherInvMass(-1),
  fSamePiPiMotherPt(-1),
  fSamePiPiPiMotherID(-1),
  fSamePiPiPiMotherInvMass(-1),
  fSamePiPiPiMotherPt(-1),
  fV0MultiplicityHNMEvent(-1),
  fTrackMultiplicityHNMEvent(-1),
  fZVertexHNMEvent(-1),
  fPtHNM(-1),
  fPDGMassNDM(-1),
  fNDMMinPtPossible(0.),
  fPDGMassChargedPion(-1),
  fPDGCodeNDM(-1),
  fPDGCodeAnalyzedMeson(-1),
  enableDalitzAllPt(kFALSE),
  enableDalitzLowPt(kFALSE),
  enableDalitzMidPt(kFALSE),
  enableDalitzHighPt(kFALSE),
  HistoDalitzPtRangeMin_LowPt(0.),
  HistoDalitzPtRangeMax_LowPt(0.),
  HistoDalitzPtRangeMin_MidPt(0.),
  HistoDalitzPtRangeMax_MidPt(0.),
  HistoDalitzPtRangeMin_HighPt(0.),
  HistoDalitzPtRangeMax_HighPt(0.),
  fHistoConvGammaPt(nullptr),
  fHistoConvGammaEta(nullptr),
  fHistoClusterGammaPt(nullptr),
  fHistoClusterGammaEta(nullptr),
  fHistoClusterGammaE(nullptr),
  fHistoNumberClusterGamma(nullptr),
  fHistoNegPionPt(nullptr),
  fHistoPosPionPt(nullptr),
  fHistoNegPionPhi(nullptr),
  fHistoPosPionPhi(nullptr),
  fHistoNegPionEta(nullptr),
  fHistoPosPionEta(nullptr),
  fHistoNegPionClsTPC(nullptr),
  fHistoPosPionClsTPC(nullptr),
  fHistoPionDCAxy(nullptr),
  fHistoPionDCAz(nullptr),
  fHistoPionTPCdEdxNSigma(nullptr),
  fHistoPionTPCdEdx(nullptr),
  fHistoPionPionInvMassPt(nullptr),
  fHistoGammaGammaInvMassPt(nullptr),
  fHistoGammaGammaInvMassPtBeforeCuts(nullptr),
  fHistoSwappingGammaGammaInvMassPt(nullptr),
  fHistoMotherInvMassPt(nullptr),
  fHistoMotherInvMassPtRejectedKinematic(nullptr),
  fHistoDalitzPlotPosFixedPzNDM(nullptr),
  fHistoDalitzPlotNegFixedPzNDM(nullptr),
  fHistoDalitzPlotPosSubNDM(nullptr),
  fHistoDalitzPlotNegSubNDM(nullptr),
  fHistoDalitzPlotPosFixedPzNDM_LowPt(nullptr),
  fHistoDalitzPlotNegFixedPzNDM_LowPt(nullptr),
  fHistoDalitzPlotPosSubNDM_LowPt(nullptr),
  fHistoDalitzPlotNegSubNDM_LowPt(nullptr),
  fHistoDalitzPlotPosFixedPzNDM_MidPt(nullptr),
  fHistoDalitzPlotNegFixedPzNDM_MidPt(nullptr),
  fHistoDalitzPlotPosSubNDM_MidPt(nullptr),
  fHistoDalitzPlotNegSubNDM_MidPt(nullptr),
  fHistoDalitzPlotPosFixedPzNDM_HighPt(nullptr),
  fHistoDalitzPlotNegFixedPzNDM_HighPt(nullptr),
  fHistoDalitzPlotPosSubNDM_HighPt(nullptr),
  fHistoDalitzPlotNegSubNDM_HighPt(nullptr),
  fHistoBackInvMassPt(nullptr),
  fHistoMotherLikeSignBackInvMassPt(nullptr),
  fHistoAngleHNMesonPiPlPiMi(nullptr),
  fHistoAngleHNMesonNDM(nullptr),
  fHistoAngleHNMesonPiPl(nullptr),
  fHistoAngleHNMesonPiMi(nullptr),
  fHistoAnglePiPlPiMi(nullptr),
  fHistoAngleNDMPiMi(nullptr),
  fHistoAnglePiPlNDM(nullptr),
  fHistoAngleSum(nullptr),
  fHistoTrueAngleSum(nullptr),
  fHistoTrueHNMesonPtvsNDMPt(nullptr),
  fHistoMotherInvMassSubNDM(nullptr),
  fHistoBackInvMassPtSubNDM(nullptr),
  fHistoMotherLikeSignBackInvMassSubNDMPt(nullptr),
  fHistoMotherInvMassFixedPzNDM(nullptr),
  fHistoBackInvMassPtFixedPzNDM(nullptr),
  fHistoMotherLikeSignBackInvMassFixedPzNDMPt(nullptr),
  fHistoMCAllGammaPt(nullptr),
  fHistoMCConvGammaPt(nullptr),
  fHistoMCAllMesonPt(nullptr),
  fHistoMCAllMesonEta(nullptr),
  fHistoMCAllMesonPhi(nullptr),
  fHistoMCMesonFromNeutralMesonPt(nullptr),
  fHistoMCMesonFromNeutralMesonEta(nullptr),
  fHistoMCMesonFromNeutralMesonPhi(nullptr),
  fHistoMCAllPosPionsPt(nullptr),
  fHistoMCAllPosPionsEta(nullptr),
  fHistoMCAllPosPionsPhi(nullptr),
  fHistoMCAllNegPionsPt(nullptr),
  fHistoMCAllNegPionsEta(nullptr),
  fHistoMCAllNegPionsPhi(nullptr),
  fHistoMCGammaFromNeutralMesonPt(nullptr),
  fHistoMCPosPionsFromNeutralMesonPt(nullptr),
  fHistoMCPosPionsFromNeutralMesonEta(nullptr),
  fHistoMCPosPionsFromNeutralMesonPhi(nullptr),
  fHistoMCNegPionsFromNeutralMesonPt(nullptr),
  fHistoMCNegPionsFromNeutralMesonEta(nullptr),
  fHistoMCNegPionsFromNeutralMesonPhi(nullptr),
  fHistoMCHNMPiPlPiMiNDMPt(nullptr),
  fHistoMCHNMPiPlPiMiNDMEta(nullptr),
  fHistoMCHNMPiPlPiMiNDMPhi(nullptr),
  fHistoMCHNMPiPlPiMiNDMEtavsPt(nullptr),
  fHistoMCHNMPiPlPiMiNDMInAccPt(nullptr),
  fHistoMCHNMInAccVsNDMPt(nullptr),
  fHistoMCHeavyAllPt(nullptr),
  fHistoMCHeavyAllEta(nullptr),
  fHistoMCHeavyAllPhi(nullptr),
  fHistoMCHeavyChannelPt(nullptr),
  fHistoMCHeavyChannelEta(nullptr),
  fHistoMCHeavyChannelPhi(nullptr),
  fHistMCChannelNDMFromHeavyPt(nullptr),
  fHistMCChannelNDMFromHeavyEta(nullptr),
  fHistMCChannelNDMFromHeavyPhi(nullptr),
  fHistMCChannelPiPlusFromHeavyPt(nullptr),
  fHistMCChannelPiPlusFromHeavyEta(nullptr),
  fHistMCChannelPiPlusFromHeavyPhi(nullptr),
  fHistMCChannelPiMinusFromHeavyPt(nullptr),
  fHistMCChannelPiMinusFromHeavyEta(nullptr),
  fHistMCChannelPiPMinusFromHeavyPhi(nullptr),
  fHistMCChannelNDMPtHeavyPt(nullptr),
  fHistMCChannelPiPlusPtHeavyPt(nullptr),
  fHistMCChannelPiMinusPtHeavyPt(nullptr),
  fHistoMCHeavyReconstructiblePt(nullptr),
  fHistoMCHeavyReconstructibleEta(nullptr),
  fHistoMCHeavyReconstructiblePhi(nullptr),
  fHistMCReconstructibleNDMFromHeavyPt(nullptr),
  fHistMCReconstructibleNDMFromHeavyEta(nullptr),
  fHistMCReconstructibleNDMFromHeavyPhi(nullptr),
  fHistMCReconstructiblePiPlusFromHeavyPt(nullptr),
  fHistMCReconstructiblePiPlusFromHeavyEta(nullptr),
  fHistMCReconstructiblePiPlusFromHeavyPhi(nullptr),
  fHistMCReconstructiblePiMinusFromHeavyPt(nullptr),
  fHistMCReconstructiblePiMinusFromHeavyEta(nullptr),
  fHistMCReconstructiblePiPMinusFromHeavyPhi(nullptr),
  fHistMCReconstructibleNDMPtHeavyPt(nullptr),
  fHistMCReconstructiblePiPlusPtHeavyPt(nullptr),
  fHistMCReconstructiblePiMinusPtHeavyPt(nullptr),
  fHistoTrueMesonFlags(nullptr),
  fHistoTrueMotherPiPlPiMiNDMInvMassPt(nullptr),
  fHistoTrueMotherPiPlPiMiNDMInvMassPtSubNDM(nullptr),
  fHistoTrueMotherPiPlPiMiNDMInvMassPtFixedPzNDM(nullptr),
  fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromDifferent(nullptr),
  fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromEtaOmega(nullptr),
  fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromRho(nullptr),
  fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromK0s(nullptr),
  fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromK0l(nullptr),
  fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromEtaPrime(nullptr),
  fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromOther(nullptr),
  fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM(nullptr),
  fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM(nullptr),
  fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM(nullptr),
  fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM(nullptr),
  fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM_LowPt(nullptr),
  fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM_LowPt(nullptr),
  fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM_LowPt(nullptr),
  fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM_LowPt(nullptr),
  fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM_MidPt(nullptr),
  fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM_MidPt(nullptr),
  fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM_MidPt(nullptr),
  fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM_MidPt(nullptr),
  fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM_HighPt(nullptr),
  fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM_HighPt(nullptr),
  fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM_HighPt(nullptr),
  fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM_HighPt(nullptr),
  fHistoTrueMotherGammaGammaInvMassPt(nullptr),
  fHistoTrueMotherGammaGammaFromHNMInvMassPt(nullptr),
  fHistoTrueConvGammaPt(nullptr),
  fHistoTrueConvGammaFromNeutralMesonPt(nullptr),
  fHistoTrueClusterGammaPt(nullptr),
  fHistoTrueClusterGammaFromNeutralMesonPt(nullptr),
  fHistoTruePosPionPt(nullptr),
  fHistoTruePosPionFromNeutralMesonPt(nullptr),
  fHistoTrueNegPionPt(nullptr),
  fHistoTrueNegPionFromNeutralMesonPt(nullptr),
  fHistoTruePionPionInvMassPt(nullptr),
  fHistoTruePionPionFromSameMotherInvMassPt(nullptr),
  fHistoTruePionPionFromHNMInvMassPt(nullptr),
  fHistoTruePiPlPiMiSameMotherFromEtaInvMassPt(nullptr),
  fHistoTruePiPlPiMiSameMotherFromOmegaInvMassPt(nullptr),
  fHistoTruePiPlPiMiSameMotherFromRhoInvMassPt(nullptr),
  fHistoTruePiPlPiMiSameMotherFromEtaPrimeInvMassPt(nullptr),
  fHistoTruePiPlPiMiSameMotherFromK0sInvMassPt(nullptr),
  fHistoTruePiPlPiMiSameMotherFromK0lInvMassPt(nullptr),
  fHistoTruePiPlPiMiSameMotherFromOtherlInvMassPt(nullptr),
  fHistoTruePiMiPiZeroSameMotherFromEtaInvMassPt(nullptr),
  fHistoTruePiMiPiZeroSameMotherFromOmegaInvMassPt(nullptr),
  fHistoTruePiMiPiZeroSameMotherFromRhoInvMassPt(nullptr),
  fHistoTruePiMiPiZeroSameMotherFromK0lInvMassPt(nullptr),
  fHistoTruePiMiPiZeroSameMotherFromOtherlInvMassPt(nullptr),
  fHistoTruePiPlPiZeroSameMotherFromEtaInvMassPt(nullptr),
  fHistoTruePiPlPiZeroSameMotherFromOmegaInvMassPt(nullptr),
  fHistoTruePiPlPiZeroSameMotherFromRhoInvMassPt(nullptr),
  fHistoTruePiPlPiZeroSameMotherFromK0lInvMassPt(nullptr),
  fHistoTruePiPlPiZeroSameMotherFromOtherInvMassPt(nullptr),
  fHistoTruePiPlPiMiNDMPureCombinatoricalInvMassPt(nullptr),
  fHistoTruePiPlPiMiNDMCombinatoricalInvMassPt(nullptr),
  fHistoTruePiPlPiMiNDMContaminationInvMassPt(nullptr),
  fHistoTruePiPlPiMiNDMContamination_Pi0_InvMassPt(nullptr),
  fHistoTruePiPlPiMiNDMContamination_PiPl_InvMassPt(nullptr),
  fHistoTruePiPlPiMiNDMContamination_PiMi_InvMassPt(nullptr),
  fHistoTruePiPlPiMiNDMContamination_Crosscheck_InvMassPt(nullptr),
  fHistoTruePiPlPiMiNDMContamination_multipel_InvMassPt(nullptr),
  fHistoDoubleCountTruePi0InvMassPt(nullptr),
  fHistoDoubleCountTrueHNMInvMassPt(nullptr),
  fHistoDoubleCountTrueConvGammaRPt(nullptr),
  fVectorDoubleCountTruePi0s(0),
  fVectorDoubleCountTrueHNMs(0),
  fVectorDoubleCountTrueConvGammas(0),
  fHistoNEvents(nullptr),
  fHistoNEventsWOWeight(nullptr),
  fProfileJetJetXSection(nullptr),
  fHistoJetJetNTrials(nullptr),
  fHistoNGoodESDTracks(nullptr),
  fProfileEtaShift(nullptr),
  fHistoSPDClusterTrackletBackground(nullptr),
  fHistovParticleChi2PerNDF(nullptr),
  fHistovParticleChi2PerNDFBothConstrained(nullptr),
  fHistovParticleChi2PerNDFOneConstrained(nullptr),
  fHistovParticledS(nullptr),
  fHistovParticledSBothConstrained(nullptr),
  fHistovParticledSOneConstrained(nullptr),
  fHistoTruevParticleChi2PerNDF(nullptr),
  fHistoTruevParticleFromSameMotherChi2PerNDF(nullptr),
  fHistoTruevParticleFromHNMChi2PerNDF(nullptr),
  fHistoTruevParticledS(nullptr),
  fHistoTruevParticleFromSameMotherdS(nullptr),
  fHistoTruevParticleFromHNMdS(nullptr),
  fRandom(0),
  fnCuts(0),
  fiCut(0),
  fNumberOfESDTracks(0),
  fMoveParticleAccordingToVertex(kFALSE),
  fIsHeavyIon(0),
  fDoMesonAnalysis(kTRUE),
  fDoMesonQA(0),
  fIsFromMBHeader(kTRUE),
  fIsFromDesiredHeader(kTRUE),
  fIsOverlappingWithOtherHeader(kFALSE),
  fAllowOverlapHeaders(kTRUE),
  fIsMC(kFALSE),
  fSelectedHeavyNeutralMeson(kFALSE),
  fDoLightOutput(kFALSE),
  fNDMRecoMode(0),
  fTolerance(-1),
  fWeightJetJetMC(1.),
  fTrackMatcherRunningMode(0),
  fMCEventPos(),
  fMCEventNeg(),
  fESDArrayPos(),
  fESDArrayNeg()
{
  DefineOutput(1, TList::Class());
}

//-----------------------------------------------------------------------------------------------
AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::~AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson()
{
  //
  // virtual destructor
  //
  cout<<"Destructor"<<endl;
  if(fGoodConvGammas){
    delete fGoodConvGammas;
  }
  if(fClusterCandidates){
    delete fClusterCandidates;
  }

  if(fNeutralDecayParticleCandidates){
    delete fNeutralDecayParticleCandidates;
  }

  if(fNeutralDecayParticleSidebandCandidates){
    delete fNeutralDecayParticleSidebandCandidates;
  }

  if(fNeutralDecayParticleSwappCandidates){
    delete fNeutralDecayParticleSwappCandidates;
  }

  if(fPosPionCandidates){
    delete fPosPionCandidates;
  }

  if(fNegPionCandidates){
    delete fNegPionCandidates;
  }

  if(fBGHandlerPiPl){
    for(int icut = 0; icut < fnCuts; icut++)
      delete fBGHandlerPiPl[icut];
    delete[] fBGHandlerPiPl;
  }

  if(fBGHandlerPiMi){
    for(int icut = 0; icut < fnCuts; icut++)
      delete fBGHandlerPiMi[icut];
    delete[] fBGHandlerPiMi;
  }
}
//___________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::InitBack(){

  fBGHandlerPiPl = new AliGammaConversionAODBGHandler*[fnCuts];
  fBGHandlerPiMi = new AliGammaConversionAODBGHandler*[fnCuts];

  for(Int_t iCut = 0; iCut<fnCuts;iCut++){

    TString cutstringEvent		= ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
    TString cutstringPion		= ((AliPrimaryPionCuts*)fPionCutArray->At(iCut))->GetCutNumber();
    TString cutstringConvGamma = "";
    if (fNDMRecoMode < 2)  cutstringConvGamma = ((AliConversionPhotonCuts*)fGammaCutArray->At(iCut))->GetCutNumber();
    TString cutstringCaloGamma = "";
    if (fNDMRecoMode > 0)  cutstringCaloGamma = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
    TString cutstringNeutralPion= ((AliConversionMesonCuts*)fNeutralDecayMesonCutArray->At(iCut))->GetCutNumber();
    TString cutstringMeson		= ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();

    TString fullCutString = "";
    if (fNDMRecoMode == 0) fullCutString = Form("%i_%s_%s_%s_%s_%s",fNDMRecoMode,cutstringEvent.Data(),cutstringConvGamma.Data(),cutstringNeutralPion.Data(), cutstringPion.Data(),cutstringMeson.Data());
    else if (fNDMRecoMode == 1) fullCutString = Form("%i_%s_%s_%s_%s_%s_%s",fNDMRecoMode,cutstringEvent.Data(),cutstringConvGamma.Data(),cutstringCaloGamma.Data(),cutstringNeutralPion.Data(), cutstringPion.Data(),cutstringMeson.Data());
    else if (fNDMRecoMode == 2) fullCutString = Form("%i_%s_%s_%s_%s_%s",fNDMRecoMode,cutstringEvent.Data(),cutstringCaloGamma.Data(),cutstringNeutralPion.Data(), cutstringPion.Data(),cutstringMeson.Data());

    Int_t collisionSystem = atoi((TString)(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber())(0,1));
    Int_t centMin = atoi((TString)(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber())(1,1));
    Int_t centMax = atoi((TString)(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber())(2,1));

    if(collisionSystem == 1 || collisionSystem == 2 ||
      collisionSystem == 5 || collisionSystem == 8 ||
      collisionSystem == 9){
      centMin = centMin*10;
      centMax = centMax*10;
    }
    else if(collisionSystem == 3 || collisionSystem == 6){
      centMin = centMin*5;
      centMax = centMax*5;
    }
    else if(collisionSystem == 4 || collisionSystem == 7){
      centMin = ((centMin*5)+45);
      centMax = ((centMax*5)+45);
    }

    fBGHandlerPiPl[iCut] = new AliGammaConversionAODBGHandler(	collisionSystem,centMin,centMax,
                                ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetNumberOfBGEvents(),
                                ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity(),
                                4,8,5);

    fBGHandlerPiMi[iCut] = new AliGammaConversionAODBGHandler(	collisionSystem,centMin,centMax,
                                ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetNumberOfBGEvents(),
                                ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity(),
                                4,8,5);
  }
}

//______________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::UserCreateOutputObjects()
{
  //
  // Create ouput objects
  //

  if(((AliConvEventCuts*)fEventCutArray->At(0))->GetUseJetFinderForOutliers()){
    fOutlierJetReader=(AliAnalysisTaskJetOutlierRemoval*)AliAnalysisManager::GetAnalysisManager()->GetTask("AliAnalysisTaskJetOutlierRemoval");
    if(!fOutlierJetReader){AliFatal("Error: No AliAnalysisTaskJetOutlierRemoval");} // GetV0Reader
    else{printf("Found AliAnalysisTaskJetOutlierRemoval used for outlier removal!\n");}
  }

  // Set pT and mass ranges
  Double_t HistoNMassBins                             = 600;
  Double_t HistoMassRange[2]                          = {0.4,1.0};
  Double_t HistoNMassBinsSub                          = 600;
  Double_t HistoMassRangeSub[2]                       = {0.4,1.0};
  Double_t HistoNPtBins                               = 800;
  Double_t HistoPtRange[2]                            = {0.,80.};
  Double_t HistoNMassBinsDecayMeson                   = 450;
  Double_t HistoMassRangeNDM[2]                       = {0.0,0.45};
  Double_t HistoNMassBinsPiPlusPiMinus                = 250;
  Double_t HistoMassRangePiPlusPiMinus[2]             = {0.0,2.0};
  TString NameNeutralMesonAnalyzed                    = "not set";
  TString NameNeutralMesonAnalyzedLatex               = "not set";
  TString NameNDM                                     = "not set";
  TString NameNDMLatex                                = "not set";
  Double_t HistoMassRangeDalitzMin                    = 0.0;
  Double_t HistoMassRangeDalitz                       = 3.0;
  //fNDMRecoMode: 0=PCM-PCM, 1=PCM-Calo, 2=Calo-Calo
  //fClusterCutArray->At(iCut))->GetClusterType(): 1=EMCAL, 2=PHOS, 3=DCAL, 4=EMCAL+DCAL, 0=All
  if (fNDMRecoMode == 0){ //PCM-PCM
    HistoDalitzPtRangeMin_LowPt                         = 1.6;
    HistoDalitzPtRangeMax_LowPt                         = 2.0;
    HistoDalitzPtRangeMin_MidPt                         = 5.0;
    HistoDalitzPtRangeMax_MidPt                         = 6.0;
    HistoDalitzPtRangeMin_HighPt                        = 10.0;
    HistoDalitzPtRangeMax_HighPt                        = 12.0;
  } else if (fNDMRecoMode == 1){ //PCM-Calo
    if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType()==2){ //PCM-PHOS
      HistoDalitzPtRangeMin_LowPt                         = 2.5;
      HistoDalitzPtRangeMax_LowPt                         = 3.0;
      HistoDalitzPtRangeMin_MidPt                         = 5.0;
      HistoDalitzPtRangeMax_MidPt                         = 6.0;
      HistoDalitzPtRangeMin_HighPt                        = 10.0;
      HistoDalitzPtRangeMax_HighPt                        = 12.0;
    } else { //PCM-EMC
      HistoDalitzPtRangeMin_LowPt                         = 2.5;
      HistoDalitzPtRangeMax_LowPt                         = 3.0;
      HistoDalitzPtRangeMin_MidPt                         = 5.0;
      HistoDalitzPtRangeMax_MidPt                         = 6.0;
      HistoDalitzPtRangeMin_HighPt                        = 10.0;
      HistoDalitzPtRangeMax_HighPt                        = 12.0;
    }
  } else { //Calo-Calo
    if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType()==2){ //PHOS-PHOS
      HistoDalitzPtRangeMin_LowPt                         = 3.5;
      HistoDalitzPtRangeMax_LowPt                         = 4.0;
      HistoDalitzPtRangeMin_MidPt                         = 5.0;
      HistoDalitzPtRangeMax_MidPt                         = 6.0;
      HistoDalitzPtRangeMin_HighPt                        = 8.0;
      HistoDalitzPtRangeMax_HighPt                        = 10.0;
    } else { //EMC-EMC
      HistoDalitzPtRangeMin_LowPt                         = 3.0;
      HistoDalitzPtRangeMax_LowPt                         = 4.0;
      HistoDalitzPtRangeMin_MidPt                         = 8.0;
      HistoDalitzPtRangeMax_MidPt                         = 10.0;
      HistoDalitzPtRangeMin_HighPt                        = 16.0;
      HistoDalitzPtRangeMax_HighPt                        = 20.0;
    }
  }

  //Enable Histograms
  if( fDoMesonQA>0 ) {
    //fNDMRecoMode: 0=PCM-PCM, 1=PCM-Calo, 2=Calo-Calo
    //fClusterCutArray->At(iCut))->GetClusterType(): 1=EMCAL, 2=PHOS, 3=DCAL, 4=EMCAL+DCAL, 0=All
    if (fNDMRecoMode == 0){ //PCM-PCM
        enableDalitzLowPt=kTRUE;
        enableDalitzMidPt=kTRUE;
        enableDalitzHighPt=kTRUE;
    } else if (fNDMRecoMode == 1){ //PCM-Calo
        enableDalitzLowPt=kTRUE;
        enableDalitzMidPt=kTRUE;
        enableDalitzHighPt=kTRUE;
    } else { //Calo-Calo
        enableDalitzLowPt=kTRUE;
        enableDalitzMidPt=kTRUE;
        enableDalitzHighPt=kTRUE;
    }
  }

  switch( fSelectedHeavyNeutralMeson ) {
  case 0: // ETA MESON
    HistoNMassBins                                    = 400;
    HistoMassRange[0]                                 = 0.3;
    HistoMassRange[1]                                 = 0.7;
    HistoNMassBinsSub                                 = 400;
    HistoMassRangeSub[0]                              = 0.3;
    HistoMassRangeSub[1]                              = 0.7;
    HistoNMassBinsDecayMeson                          = 450;
    HistoMassRangeNDM[0]                              = 0.0;
    HistoMassRangeNDM[1]                              = 0.45;
    NameNeutralMesonAnalyzed                          = "Eta";
    NameNeutralMesonAnalyzedLatex                     = "#eta";
    NameNDM                                           = "NeutralPion";
    NameNDMLatex                                      = "#pi^{0}";
    fPDGMassNDM                                       = 0.1349766; // hard coded PDG value to keep results reproducable later
    fPDGCodeNDM                                       = 111; // PDG pi0
    fPDGMassChargedPion                               = 0.1395706; // hard coded PDG 2018 value to keep results reproducable later
    fPDGCodeAnalyzedMeson                             = 221; // PDG omega
    break;
  case 1: // OMEGA MESON
    HistoNMassBins                                    = 500;
    HistoMassRange[0]                                 = 0.5;
    HistoMassRange[1]                                 = 1.0;
    HistoNMassBinsSub                                 = 500;
    HistoMassRangeSub[0]                              = 0.5;
    HistoMassRangeSub[1]                              = 1.0;
    HistoNMassBinsDecayMeson                          = 450;
    HistoMassRangeNDM[0]                              = 0.0;
    HistoMassRangeNDM[1]                              = 0.45;
    NameNeutralMesonAnalyzed                          = "Omega";
    NameNeutralMesonAnalyzedLatex                     = "#omega";
    NameNDM                                           = "NeutralPion";
    NameNDMLatex                                      = "#pi^{0}";
    fPDGMassNDM                                       = 0.1349766; // hard coded PDG value to keep results reproducable later
    fPDGCodeNDM                                       = 111; // PDG pi0
    fPDGMassChargedPion                               = 0.1395706; // hard coded PDG 2018 value to keep results reproducable later
    fPDGCodeAnalyzedMeson                             = 223; // PDG omega
    break;
  case 2: // ETA PRIME MESON
    HistoNMassBins                                    = 600;
    HistoMassRange[0]                                 = 0.6;
    HistoMassRange[1]                                 = 1.2;
    HistoNMassBinsSub                                 = 600;
    HistoMassRangeSub[0]                              = 0.6;
    HistoMassRangeSub[1]                              = 1.2;
    HistoNMassBinsDecayMeson                          = 450;
    HistoMassRangeNDM[0]                              = 0.4;
    HistoMassRangeNDM[1]                              = 0.85;
    NameNeutralMesonAnalyzed                          = "EtaPrime";
    NameNeutralMesonAnalyzedLatex                     = "#eta'";
    NameNDM                                           = "EtaMeson";
    NameNDMLatex                                      = "#eta";
    fPDGMassNDM                                       = 0.547862; // hard coded PDG value to keep results reproducable later
    fPDGCodeNDM                                       = 221; // PDG value eta
    fPDGMassChargedPion                               = 0.1395706; // hard coded PDG 2018 value to keep results reproducable later
    fPDGCodeAnalyzedMeson                             = 331; // PDG value eta prime
    HistoMassRangeDalitz                              = 3.0;
    break;
  case 3:         // D0 MESON
    HistoNMassBins                                    = 600;
    HistoMassRange[0]                                 = 1.4;
    HistoMassRange[1]                                 = 2.0;
    HistoNMassBinsSub                                 = 400;      // MF: Background subtraction?
    HistoMassRangeSub[0]                              = 1.6;
    HistoMassRangeSub[1]                              = 2.0;
    HistoNMassBinsDecayMeson                          = 450;
    HistoMassRangeNDM[0]                              = 0.0;
    HistoMassRangeNDM[1]                              = 0.45;
    HistoNPtBins                                      = 500;
    HistoPtRange[1]                                   = 50.;
    NameNeutralMesonAnalyzed                          = "D0";
    NameNeutralMesonAnalyzedLatex                     = "D^{0}";
    NameNDM                                           = "NeutralPion";
    NameNDMLatex                                      = "#pi^{0}";
    fPDGMassNDM                                       = 0.1349766; // hard coded PDG value to keep results reproducable later
    fPDGCodeNDM                                       = 111; // PDG value pi0
    fPDGMassChargedPion                               = 0.1395706; // hard coded PDG 2018 value to keep results reproducable later
    fPDGCodeAnalyzedMeson                             = 421; // PDG value D0
    break;
  default:
    AliError(Form("Heavy neutral meson not correctly selected (only 0,1,2,3 valid)... selected: %d",fSelectedHeavyNeutralMeson));
  }

  // Create the output container
  if(fOutputContainer != nullptr){
    delete fOutputContainer;
    fOutputContainer            = nullptr;
  }
  if(fOutputContainer == nullptr){
    fOutputContainer            = new TList();
    fOutputContainer->SetOwner(kTRUE);
  }

  fGoodConvGammas               = new TList();
  fClusterCandidates            = new TList();
  fClusterCandidates->SetOwner(kTRUE);

  fNeutralDecayParticleCandidates        = new TList();
  fNeutralDecayParticleCandidates->SetOwner(kTRUE);

  fNeutralDecayParticleSidebandCandidates        = new TList();
  fNeutralDecayParticleSidebandCandidates->SetOwner(kTRUE);

  fNeutralDecayParticleSwappCandidates        = new TList();
  fNeutralDecayParticleSwappCandidates->SetOwner(kTRUE);


  fPosPionCandidates            = new TList();
  fPosPionCandidates->SetOwner(kTRUE);
  fNegPionCandidates            = new TList();
  fNegPionCandidates->SetOwner(kTRUE);

  fCutFolder                    = new TList*[fnCuts];
  fESDList                      = new TList*[fnCuts];
  fHistoNEvents                 = new TH1F*[fnCuts];
  if(fIsMC > 1){
    fHistoNEventsWOWeight       = new TH1F*[fnCuts];
  }
  fHistoNGoodESDTracks          = new TH1I*[fnCuts];
  if(fIsMC == 2) {
    fProfileJetJetXSection       = new TProfile*[fnCuts];
    fHistoJetJetNTrials          = new TH1F*[fnCuts];
  }
  if (fNDMRecoMode > 0){
    fHistoClusterGammaPt        = new TH1F*[fnCuts];
    fHistoClusterGammaEta       = new TH1F*[fnCuts];
    fHistoClusterGammaE         = new TH1F*[fnCuts];
    if( fDoMesonQA>1 ) {
      fHistoNumberClusterGamma         = new TH1I*[fnCuts];
    }
  }
  if(!fDoLightOutput){
    fProfileEtaShift              = new TProfile*[fnCuts];
    fHistoSPDClusterTrackletBackground = new TH2F*[fnCuts];

    fHistovParticleChi2PerNDF = new TH1F*[fnCuts];
    fHistovParticleChi2PerNDFBothConstrained  = new TH1F*[fnCuts];
    fHistovParticleChi2PerNDFOneConstrained = new TH1F*[fnCuts];
    fHistovParticledS  = new TH1F*[fnCuts];
    fHistovParticledSBothConstrained  = new TH1F*[fnCuts];
    fHistovParticledSOneConstrained  = new TH1F*[fnCuts];
    if (fNDMRecoMode < 2){
      fHistoConvGammaPt           = new TH1F*[fnCuts];
      fHistoConvGammaEta          = new TH1F*[fnCuts];
    }
    fHistoNegPionPt               = new TH1F*[fnCuts];
    fHistoPosPionPt               = new TH1F*[fnCuts];
    fHistoNegPionPhi              = new TH1F*[fnCuts];
    fHistoPosPionPhi              = new TH1F*[fnCuts];
    fHistoPionPionInvMassPt       = new TH2F*[fnCuts];

    if( fDoMesonQA>0 ) {
      fHistoNegPionEta            = new TH1F*[fnCuts];
      fHistoPosPionEta            = new TH1F*[fnCuts];
      fHistoNegPionClsTPC         = new TH2F*[fnCuts];
      fHistoPosPionClsTPC         = new TH2F*[fnCuts];
      fHistoPionDCAxy             = new TH2F*[fnCuts];
      fHistoPionDCAz              = new TH2F*[fnCuts];
      fHistoPionTPCdEdxNSigma     = new TH2F*[fnCuts];
      fHistoPionTPCdEdx           = new TH2F*[fnCuts];
      if (enableDalitzAllPt){
        fHistoDalitzPlotPosFixedPzNDM          = new TH2F*[fnCuts];
        fHistoDalitzPlotNegFixedPzNDM          = new TH2F*[fnCuts];
        fHistoDalitzPlotPosSubNDM              = new TH2F*[fnCuts];
        fHistoDalitzPlotNegSubNDM              = new TH2F*[fnCuts];
      }
      if (enableDalitzLowPt){
        fHistoDalitzPlotPosFixedPzNDM_LowPt       = new TH2F*[fnCuts];
        fHistoDalitzPlotNegFixedPzNDM_LowPt       = new TH2F*[fnCuts];
        fHistoDalitzPlotPosSubNDM_LowPt           = new TH2F*[fnCuts];
        fHistoDalitzPlotNegSubNDM_LowPt           = new TH2F*[fnCuts];
      }
      if (enableDalitzMidPt){
        fHistoDalitzPlotPosFixedPzNDM_MidPt       = new TH2F*[fnCuts];
        fHistoDalitzPlotNegFixedPzNDM_MidPt       = new TH2F*[fnCuts];
        fHistoDalitzPlotPosSubNDM_MidPt           = new TH2F*[fnCuts];
        fHistoDalitzPlotNegSubNDM_MidPt           = new TH2F*[fnCuts];
      }
      if (enableDalitzHighPt){
        fHistoDalitzPlotPosFixedPzNDM_HighPt      = new TH2F*[fnCuts];
        fHistoDalitzPlotNegFixedPzNDM_HighPt      = new TH2F*[fnCuts];
        fHistoDalitzPlotPosSubNDM_HighPt          = new TH2F*[fnCuts];
        fHistoDalitzPlotNegSubNDM_HighPt          = new TH2F*[fnCuts];
      }
    }


    fHistoAngleHNMesonPiPlPiMi    = new TH2F*[fnCuts];
    fHistoAngleHNMesonNDM      = new TH2F*[fnCuts];
    fHistoAngleHNMesonPiPl        = new TH2F*[fnCuts];
    fHistoAngleHNMesonPiMi        = new TH2F*[fnCuts];
    fHistoAngleNDMPiMi       = new TH2F*[fnCuts];
    fHistoAnglePiPlPiMi         = new TH2F*[fnCuts];
    fHistoAnglePiPlNDM       = new TH2F*[fnCuts];
    fHistoAngleSum              = new TH2F*[fnCuts];
  }

  fHistoGammaGammaInvMassPt               = new TH2F*[fnCuts];
  fHistoGammaGammaInvMassPtBeforeCuts     = new TH2F*[fnCuts];
  if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetBackgroundMode() == 7){
    fHistoSwappingGammaGammaInvMassPt       = new TH2F*[fnCuts];
  }
  fHistoMotherInvMassPt                   = new TH2F*[fnCuts];
  fHistoMotherInvMassPtRejectedKinematic  = new TH2F*[fnCuts];
  fHistoBackInvMassPt                     = new TH2F*[fnCuts];

  fHistoMotherLikeSignBackInvMassPt       = new TH2F*[fnCuts];

  fHistoMotherInvMassSubNDM                               = new TH2F*[fnCuts];
  fHistoBackInvMassPtSubNDM = new TH2F*[fnCuts];

  fHistoMotherLikeSignBackInvMassSubNDMPt       = new TH2F*[fnCuts];

  fHistoMotherInvMassFixedPzNDM                     = new TH2F*[fnCuts];
  fHistoBackInvMassPtFixedPzNDM = new TH2F*[fnCuts];
  fHistoMotherLikeSignBackInvMassFixedPzNDMPt       = new TH2F*[fnCuts];

  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    TString cutstringEvent        = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
    TString cutstringPion         = ((AliPrimaryPionCuts*)fPionCutArray->At(iCut))->GetCutNumber();
    TString cutstringConvGamma    = "";
    if (fNDMRecoMode < 2)
      cutstringConvGamma          = ((AliConversionPhotonCuts*)fGammaCutArray->At(iCut))->GetCutNumber();
    TString cutstringCaloGamma    = "";
    if (fNDMRecoMode > 0)
      cutstringCaloGamma          = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
    TString cutstringNeutralPion  = ((AliConversionMesonCuts*)fNeutralDecayMesonCutArray->At(iCut))->GetCutNumber();
    TString cutstringMeson        = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();

    TString fullCutString         = "";
    if (fNDMRecoMode == 0)
      fullCutString               = Form("%i_%s_%s_%s_%s_%s",fNDMRecoMode,cutstringEvent.Data(),cutstringConvGamma.Data(),cutstringNeutralPion.Data(), cutstringPion.Data(),cutstringMeson.Data());
    else if (fNDMRecoMode == 1)
      fullCutString               = Form("%i_%s_%s_%s_%s_%s_%s",fNDMRecoMode,cutstringEvent.Data(),cutstringConvGamma.Data(),cutstringCaloGamma.Data(),cutstringNeutralPion.Data(), cutstringPion.Data(),cutstringMeson.Data());
    else if (fNDMRecoMode == 2)
      fullCutString               = Form("%i_%s_%s_%s_%s_%s",fNDMRecoMode,cutstringEvent.Data(),cutstringCaloGamma.Data(),cutstringNeutralPion.Data(), cutstringPion.Data(),cutstringMeson.Data());
    TString nameCutFolder         = Form("Cut Number %s", fullCutString.Data());
    TString nameESDList           = Form("%s ESD histograms", fullCutString.Data());

    // Set min pt of pi0 that each method is able to reconstruct -> will be used for MC studies
    if(fSelectedHeavyNeutralMeson == 1){ // omega
      if(fNDMRecoMode==0){ // PCM
        fNDMMinPtPossible                              = 0.3;
      } else if (fNDMRecoMode==1){ // mixed
        if(cutstringCaloGamma(0,1).String().EqualTo("1") || cutstringCaloGamma(0,1).String().EqualTo("4")){
          fNDMMinPtPossible                              = 0.8; //PCM-EMC
        } else{
          fNDMMinPtPossible                              = 0.4; //PCM-PHOS
        }
      } else if (fNDMRecoMode==2){ // pure
        if(cutstringCaloGamma(0,1).String().EqualTo("1") || cutstringCaloGamma(0,1).String().EqualTo("4")){
          fNDMMinPtPossible                              = 1.4; //EMC
        } else{
          fNDMMinPtPossible                              = 1.6; //PHOS
        }
      } else{
        fNDMMinPtPossible = 0.;
      }
    }

    fCutFolder[iCut]              = new TList();
    fCutFolder[iCut]->SetName(nameCutFolder.Data());
    fCutFolder[iCut]->SetOwner(kTRUE);
    fOutputContainer->Add(fCutFolder[iCut]);

    fESDList[iCut]                = new TList();
    fESDList[iCut]->SetName(nameESDList.Data());
    fESDList[iCut]->SetOwner(kTRUE);

    fHistoNEvents[iCut]           = new TH1F("NEvents","NEvents",14,-0.5,13.5);
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(3,"Miss. MC or inc. ev.");
    if (((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger() > 1 ){
      TString TriggerNames  = "Not Trigger: ";
      TriggerNames          = TriggerNames+ ( (AliConvEventCuts*)fEventCutArray->At(iCut))->GetSpecialTriggerName();
      fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(4,TriggerNames.Data());
    } else {
      fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(4,"Trigger");
    }
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(5,"Vertex Z");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(7,"Pile-Up");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(8,"no SDD");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(9,"no V0AND");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(10,"EMCAL/TPC problem");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(14,"Pileup V0M-TPCout Tracks");
    fHistoNEvents[iCut]->GetYaxis()->SetTitle("N_{events}");
    fESDList[iCut]->Add(fHistoNEvents[iCut]);

    if(fIsMC > 1){
      fHistoNEventsWOWeight[iCut]           = new TH1F("NEventsWOWeight","NEventsWOWeight",14,-0.5,13.5);
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(3,"Miss. MC or inc. ev.");
      if (((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger() > 1 ){
        TString TriggerNames  = "Not Trigger: ";
        TriggerNames          = TriggerNames+ ( (AliConvEventCuts*)fEventCutArray->At(iCut))->GetSpecialTriggerName();
        fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(4,TriggerNames.Data());
      } else {
        fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(4,"Trigger");
      }
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(5,"Vertex Z");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(7,"Pile-Up");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(8,"no SDD");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(9,"no V0AND");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(10,"EMCAL problem");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(14,"Pileup V0M-TPCout Tracks");
      fHistoNEventsWOWeight[iCut]->GetYaxis()->SetTitle("N_{events}");
      fESDList[iCut]->Add(fHistoNEventsWOWeight[iCut]);
    }

    if (fIsMC == 2){
      fProfileJetJetXSection[iCut]  = new TProfile("XSection", "XSection", 1, -0.5, 0.5);
      fESDList[iCut]->Add(fProfileJetJetXSection[iCut]);
      fHistoJetJetNTrials[iCut]     = new TH1F("NTrials", "#sum{NTrials}", 1, 0, 1);
      fHistoJetJetNTrials[iCut]->GetXaxis()->SetBinLabel(1,"#sum{NTrials}");
      fESDList[iCut]->Add(fHistoJetJetNTrials[iCut]);
    }

    if(fIsHeavyIon>0)
      fHistoNGoodESDTracks[iCut]  = new TH1I("GoodESDTracks","GoodESDTracks",3000,0,3000);
    else
      fHistoNGoodESDTracks[iCut]  = new TH1I("GoodESDTracks","GoodESDTracks",200,0,200);
    fHistoNGoodESDTracks[iCut]->GetXaxis()->SetTitle("N_{good ESD tracks}");
    fHistoNGoodESDTracks[iCut]->GetYaxis()->SetTitle("N_{events}");
    fHistoNGoodESDTracks[iCut]->Sumw2();
    fESDList[iCut]->Add(fHistoNGoodESDTracks[iCut]);

    if (fNDMRecoMode > 0){
      fHistoClusterGammaPt[iCut]  = new TH1F("ESD_ClusterGamma_Pt","ESD_ClusterGamma_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
      fHistoClusterGammaPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      fHistoClusterGammaPt[iCut]->GetYaxis()->SetTitle("N_{#gamma,cluster}");
      fHistoClusterGammaPt[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoClusterGammaPt[iCut]);
      fHistoClusterGammaEta[iCut] = new TH1F("ESD_ClusterGamma_Eta","ESD_ClusterGamma_Eta",600,-1.5,1.5);
      fHistoClusterGammaEta[iCut]->GetXaxis()->SetTitle("#eta");
      fHistoClusterGammaEta[iCut]->GetYaxis()->SetTitle("N_{#gamma,cluster}");
      fHistoClusterGammaEta[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoClusterGammaEta[iCut]);
      fHistoClusterGammaE[iCut]  = new TH1F("ESD_ClusterGamma_E","ESD_ClusterGamma_E",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
      fHistoClusterGammaE[iCut]->GetXaxis()->SetTitle("E (GeV)");
      fHistoClusterGammaE[iCut]->GetYaxis()->SetTitle("N_{#gamma,cluster}");
      fHistoClusterGammaE[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoClusterGammaE[iCut]);
      if ( fDoMesonQA>1 ) {
        fHistoNumberClusterGamma[iCut]  = new TH1I("ESD_ClusterGamma","ESD_ClusterGamma",10,0,10);
        fHistoNumberClusterGamma[iCut]->GetXaxis()->SetTitle("N_{#gamma}");
        fHistoNumberClusterGamma[iCut]->GetYaxis()->SetTitle("Events");
        fHistoNumberClusterGamma[iCut]->Sumw2();
        fESDList[iCut]->Add(fHistoNumberClusterGamma[iCut]);
      }
    }
    if(!fDoLightOutput){
      fProfileEtaShift[iCut]        = new TProfile("Eta Shift","Eta Shift",1, -0.5,0.5);
      fESDList[iCut]->Add(fProfileEtaShift[iCut]);
      fHistoSPDClusterTrackletBackground[iCut] = new TH2F("SPD tracklets vs SPD clusters","SPD tracklets vs SPD clusters",100,0,200,250,0,1000);
      fHistoSPDClusterTrackletBackground[iCut]->GetXaxis()->SetTitle("N_{SPD tracklets}");
      fHistoSPDClusterTrackletBackground[iCut]->GetYaxis()->SetTitle("N_{SPD clusters}");
      fHistoSPDClusterTrackletBackground[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoSPDClusterTrackletBackground[iCut]);

      fHistovParticleChi2PerNDF[iCut] = new TH1F("fHistovParticleChi2PerNDF","fHistovParticleChi2PerNDF",300,0,300);
      fHistovParticleChi2PerNDF[iCut]->GetXaxis()->SetTitle("#chi^{2}/ndf");
      fHistovParticleChi2PerNDF[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistovParticleChi2PerNDF[iCut]);

      fHistovParticleChi2PerNDFBothConstrained[iCut] = new TH1F("fHistovParticleChi2PerNDFBothConstrained","fHistovParticleChi2PerNDFBothConstrained",300,0,300);
      fHistovParticleChi2PerNDFBothConstrained[iCut]->GetXaxis()->SetTitle("#chi^{2}/ndf");
      fHistovParticleChi2PerNDFBothConstrained[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistovParticleChi2PerNDFBothConstrained[iCut]);

      fHistovParticleChi2PerNDFOneConstrained[iCut] = new TH1F("fHistovParticleChi2PerNDFOneConstrained","fHistovParticleChi2PerNDFOneConstrained",300,0,300);
      fHistovParticleChi2PerNDFOneConstrained[iCut]->GetXaxis()->SetTitle("#chi^{2}/ndf");
      fHistovParticleChi2PerNDFOneConstrained[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistovParticleChi2PerNDFOneConstrained[iCut]);

      fHistovParticledS[iCut] = new TH1F("fHistovParticledS","fHistovParticledS",400,-4,4);
      fHistovParticledS[iCut]->GetXaxis()->SetTitle("dS");
      fHistovParticledS[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistovParticledS[iCut]);

      fHistovParticledSBothConstrained[iCut] = new TH1F("fHistovParticledSBothConstrained","fHistovParticledSBothConstrained",400,-4,4);
      fHistovParticledSBothConstrained[iCut]->GetXaxis()->SetTitle("dS");
      fHistovParticledSBothConstrained[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistovParticledSBothConstrained[iCut]);

      fHistovParticledSOneConstrained[iCut] = new TH1F("fHistovParticledSOneConstrained","fHistovParticledSOneConstrained",400,-4,4);
      fHistovParticledSOneConstrained[iCut]->GetXaxis()->SetTitle("dS");
      fHistovParticledSOneConstrained[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistovParticledSOneConstrained[iCut]);

      if (fNDMRecoMode < 2){
        fHistoConvGammaPt[iCut]     = new TH1F("ESD_ConvGamma_Pt","ESD_ConvGamma_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
        fHistoConvGammaPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoConvGammaPt[iCut]->GetYaxis()->SetTitle("N_{#gamma,conv}");
        fHistoConvGammaPt[iCut]->Sumw2();
        fESDList[iCut]->Add(fHistoConvGammaPt[iCut]);
        fHistoConvGammaEta[iCut]    = new TH1F("ESD_ConvGamma_Eta","ESD_ConvGamma_Eta",600,-1.5,1.5);
        fHistoConvGammaEta[iCut]->GetXaxis()->SetTitle("#eta");
        fHistoConvGammaEta[iCut]->GetYaxis()->SetTitle("N_{#gamma,conv}");
        fHistoConvGammaEta[iCut]->Sumw2();
        fESDList[iCut]->Add(fHistoConvGammaEta[iCut]);
      }

      fHistoNegPionPt[iCut]         = new TH1F("ESD_PrimaryNegPions_Pt","ESD_PrimaryNegPions_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
      fHistoNegPionPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      fHistoNegPionPt[iCut]->GetYaxis()->SetTitle("N_{#pi^{-}}");
      fHistoNegPionPt[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoNegPionPt[iCut]);
      fHistoPosPionPt[iCut]         = new TH1F("ESD_PrimaryPosPions_Pt","ESD_PrimaryPosPions_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
      fHistoPosPionPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      fHistoPosPionPt[iCut]->GetYaxis()->SetTitle("N_{#pi^{+}}");
      fHistoPosPionPt[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoPosPionPt[iCut]);
      fHistoNegPionPhi[iCut]        = new TH1F("ESD_PrimaryNegPions_Phi","ESD_PrimaryNegPions_Phi",360,0,2*TMath::Pi());
      fHistoNegPionPhi[iCut]->GetXaxis()->SetTitle("#phi");
      fHistoNegPionPhi[iCut]->GetYaxis()->SetTitle("N_{#pi^{-}}");
      fHistoNegPionPhi[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoNegPionPhi[iCut]);
      fHistoPosPionPhi[iCut]        = new TH1F("ESD_PrimaryPosPions_Phi","ESD_PrimaryPosPions_Phi",360,0,2*TMath::Pi());
      fHistoPosPionPhi[iCut]->GetXaxis()->SetTitle("#phi");
      fHistoPosPionPhi[iCut]->GetYaxis()->SetTitle("N_{#pi^{+}}");
      fHistoPosPionPhi[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoPosPionPhi[iCut]);
      fHistoPionPionInvMassPt[iCut] = new TH2F("ESD_PiPlusPiNeg_InvMassPt","ESD_PiPlusPiNeg_InvMassPt",HistoNMassBinsPiPlusPiMinus,HistoMassRangePiPlusPiMinus[0],HistoMassRangePiPlusPiMinus[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
      fHistoPionPionInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} (GeV/c^{2})");
      fHistoPionPionInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
      fHistoPionPionInvMassPt[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoPionPionInvMassPt[iCut]);

      if ( fDoMesonQA>0 ) {
        fHistoNegPionEta[iCut]        = new TH1F("ESD_PrimaryNegPions_Eta","ESD_PrimaryNegPions_Eta",600,-1.5,1.5);
        fHistoNegPionEta[iCut]->GetXaxis()->SetTitle("#eta");
        fHistoNegPionEta[iCut]->GetYaxis()->SetTitle("N_{#pi^{-}}");
        fHistoNegPionEta[iCut]->Sumw2();
        fESDList[iCut]->Add(fHistoNegPionEta[iCut]);
        fHistoPosPionEta[iCut]        = new TH1F("ESD_PrimaryPosPions_Eta","ESD_PrimaryPosPions_Eta",600,-1.5,1.5);
        fHistoPosPionEta[iCut]->GetXaxis()->SetTitle("#eta");
        fHistoPosPionEta[iCut]->GetYaxis()->SetTitle("N_{#pi^{+}}");
        fHistoPosPionEta[iCut]->Sumw2();
        fESDList[iCut]->Add(fHistoPosPionEta[iCut]);
        fHistoNegPionClsTPC[iCut]     = new TH2F("ESD_PrimaryNegPions_ClsTPC","ESD_PrimaryNegPions_ClsTPC",100,0,1,400,0.,10.);
        fHistoNegPionClsTPC[iCut]->GetXaxis()->SetTitle("N_{findable cls. TPC #pi^{-}}");
        fHistoNegPionClsTPC[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoNegPionClsTPC[iCut]->Sumw2();
        fESDList[iCut]->Add(fHistoNegPionClsTPC[iCut]);
        fHistoPosPionClsTPC[iCut]     = new TH2F("ESD_PrimaryPosPions_ClsTPC","ESD_PrimaryPosPions_ClsTPC",100,0,1,400,0.,10.);
        fHistoPosPionClsTPC[iCut]->GetXaxis()->SetTitle("N_{findable cls. TPC #pi^{+}}");
        fHistoPosPionClsTPC[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoPosPionClsTPC[iCut]->Sumw2();
        fESDList[iCut]->Add(fHistoPosPionClsTPC[iCut]);
        fHistoPionDCAxy[iCut]         = new TH2F("ESD_PrimaryPions_DCAxy","ESD_PrimaryPions_DCAxy",800,-4.0,4.0,400,0.,10.);
        fHistoPionDCAxy[iCut]->GetXaxis()->SetTitle("DCA_{xy}");
        fHistoPionDCAxy[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoPionDCAxy[iCut]->Sumw2();
        fESDList[iCut]->Add(fHistoPionDCAxy[iCut]);
        fHistoPionDCAz[iCut]          = new TH2F("ESD_PrimaryPions_DCAz","ESD_PrimaryPions_DCAz",800,-4.0,4.0,400,0.,10.);
        fHistoPionDCAz[iCut]->GetXaxis()->SetTitle("DCA_{z}");
        fHistoPionDCAz[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoPionDCAz[iCut]->Sumw2();
        fESDList[iCut]->Add(fHistoPionDCAz[iCut]);
        fHistoPionTPCdEdxNSigma[iCut] = new TH2F("ESD_PrimaryPions_TPCdEdx","ESD_PrimaryPions_TPCdEdx",150,0.05,20,400,-10,10);
        fHistoPionTPCdEdxNSigma[iCut]->GetXaxis()->SetTitle("p (GeV/c)");
        fHistoPionTPCdEdxNSigma[iCut]->GetYaxis()->SetTitle("#sigma_{PID,TPC}");
        fHistoPionTPCdEdxNSigma[iCut]->Sumw2();
        fESDList[iCut]->Add(fHistoPionTPCdEdxNSigma[iCut]);
        fHistoPionTPCdEdx[iCut]       = new TH2F("ESD_PrimaryPions_TPCdEdxSignal","ESD_PrimaryPions_TPCdEdxSignal" ,150,0.05,20.0,800,0.0,200);
        fHistoPionTPCdEdx[iCut]->GetXaxis()->SetTitle("p (GeV/c)");
        fHistoPionTPCdEdx[iCut]->GetYaxis()->SetTitle("dE/dx signal (au)");
        fHistoPionTPCdEdx[iCut]->Sumw2();
        fESDList[iCut]->Add(fHistoPionTPCdEdx[iCut]);
        //---------------------------------------
        //Dalitz All Pt
        if (enableDalitzAllPt){
          fHistoDalitzPlotPosFixedPzNDM[iCut]          = new TH2F("ESD_DalitzPlotPos_FixedPz","ESD_DalitzPlotPos_FixedPz", HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz, HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz);
          fHistoDalitzPlotPosFixedPzNDM[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-}} (GeV/c^{2})");
          fHistoDalitzPlotPosFixedPzNDM[iCut]->GetYaxis()->SetTitle(Form("M_{#pi^{+} %s} (GeV/c^{2}",NameNDMLatex.Data()));
          fHistoDalitzPlotPosFixedPzNDM[iCut]->Sumw2();
          fESDList[iCut]->Add(fHistoDalitzPlotPosFixedPzNDM[iCut]);
          fHistoDalitzPlotNegFixedPzNDM[iCut]          = new TH2F("ESD_DalitzPlotNeg_FixedPz","ESD_DalitzPlotNeg_FixedPz", HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz, HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz);
          fHistoDalitzPlotNegFixedPzNDM[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-}} (GeV/c^{2})");
          fHistoDalitzPlotNegFixedPzNDM[iCut]->GetYaxis()->SetTitle(Form("M_{#pi^{-} %s} (GeV/c^{2}",NameNDMLatex.Data()));
          fHistoDalitzPlotNegFixedPzNDM[iCut]->Sumw2();
          fESDList[iCut]->Add(fHistoDalitzPlotNegFixedPzNDM[iCut]);
          fHistoDalitzPlotPosSubNDM[iCut]              = new TH2F("ESD_DalitzPlotPos_Sub","ESD_DalitzPlotPos_Sub", HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz, HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz);
          fHistoDalitzPlotPosSubNDM[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-}} (GeV/c^{2})");
          fHistoDalitzPlotPosSubNDM[iCut]->GetYaxis()->SetTitle(Form("M_{#pi^{+} %s} (GeV/c^{2}",NameNDMLatex.Data()));
          fHistoDalitzPlotPosSubNDM[iCut]->Sumw2();
          fESDList[iCut]->Add(fHistoDalitzPlotPosSubNDM[iCut]);
          fHistoDalitzPlotNegSubNDM[iCut]              = new TH2F("ESD_DalitzPlotNeg_Sub","ESD_DalitzPlotNeg_Sub", HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz, HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz);
          fHistoDalitzPlotNegSubNDM[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-}} (GeV/c^{2})");
          fHistoDalitzPlotNegSubNDM[iCut]->GetYaxis()->SetTitle(Form("M_{#pi^{-} %s} (GeV/c^{2}",NameNDMLatex.Data()));
          fHistoDalitzPlotNegSubNDM[iCut]->Sumw2();
          fESDList[iCut]->Add(fHistoDalitzPlotNegSubNDM[iCut]);
        }
        //---------------------------------------
        //Dalitz Low Pt
        if (enableDalitzLowPt){
          fHistoDalitzPlotPosFixedPzNDM_LowPt[iCut]          = new TH2F("ESD_DalitzPlotPos_FixedPz_LowPt", Form("ESD_DalitzPlotPos_FixedPz_LowPt (%3.1f to %3.1f)Gev/c", HistoDalitzPtRangeMin_LowPt, HistoDalitzPtRangeMax_LowPt), HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz, HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz);
          fHistoDalitzPlotPosFixedPzNDM_LowPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-}} (GeV/c^{2})");
          fHistoDalitzPlotPosFixedPzNDM_LowPt[iCut]->GetYaxis()->SetTitle(Form("M_{#pi^{+} %s} (GeV/c^{2}",NameNDMLatex.Data()));
          fHistoDalitzPlotPosFixedPzNDM_LowPt[iCut]->Sumw2();
          fESDList[iCut]->Add(fHistoDalitzPlotPosFixedPzNDM_LowPt[iCut]);
          fHistoDalitzPlotNegFixedPzNDM_LowPt[iCut]          = new TH2F("ESD_DalitzPlotNeg_FixedPz_LowPt", Form("ESD_DalitzPlotNeg_FixedPz_LowPt (%3.1f to %3.1f)Gev/c", HistoDalitzPtRangeMin_LowPt, HistoDalitzPtRangeMax_LowPt), HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz, HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz);
          fHistoDalitzPlotNegFixedPzNDM_LowPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-}} (GeV/c^{2})");
          fHistoDalitzPlotNegFixedPzNDM_LowPt[iCut]->GetYaxis()->SetTitle(Form("M_{#pi^{-} %s} (GeV/c^{2}",NameNDMLatex.Data()));
          fHistoDalitzPlotNegFixedPzNDM_LowPt[iCut]->Sumw2();
          fESDList[iCut]->Add(fHistoDalitzPlotNegFixedPzNDM_LowPt[iCut]);
          fHistoDalitzPlotPosSubNDM_LowPt[iCut]              = new TH2F("ESD_DalitzPlotPos_Sub_LowPt", Form("ESD_DalitzPlotPos_Sub_LowPt (%3.1f to %3.1f)Gev/c", HistoDalitzPtRangeMin_LowPt, HistoDalitzPtRangeMax_LowPt), HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz, HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz);
          fHistoDalitzPlotPosSubNDM_LowPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-}} (GeV/c^{2})");
          fHistoDalitzPlotPosSubNDM_LowPt[iCut]->GetYaxis()->SetTitle(Form("M_{#pi^{+} %s} (GeV/c^{2}",NameNDMLatex.Data()));
          fHistoDalitzPlotPosSubNDM_LowPt[iCut]->Sumw2();
          fESDList[iCut]->Add(fHistoDalitzPlotPosSubNDM_LowPt[iCut]);
          fHistoDalitzPlotNegSubNDM_LowPt[iCut]              = new TH2F("ESD_DalitzPlotNeg_Sub_LowPt", Form("ESD_DalitzPlotNeg_Sub_LowPt (%3.1f to %3.1f)Gev/c", HistoDalitzPtRangeMin_LowPt, HistoDalitzPtRangeMax_LowPt), HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz, HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz);
          fHistoDalitzPlotNegSubNDM_LowPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-}} (GeV/c^{2})");
          fHistoDalitzPlotNegSubNDM_LowPt[iCut]->GetYaxis()->SetTitle(Form("M_{#pi^{-} %s} (GeV/c^{2}",NameNDMLatex.Data()));
          fHistoDalitzPlotNegSubNDM_LowPt[iCut]->Sumw2();
          fESDList[iCut]->Add(fHistoDalitzPlotNegSubNDM_LowPt[iCut]);
        }
        //---------------------------------------
        //Dalitz Mid Pt
        if (enableDalitzMidPt){
          fHistoDalitzPlotPosFixedPzNDM_MidPt[iCut]          = new TH2F("ESD_DalitzPlotPos_FixedPz_MidPt", Form("ESD_DalitzPlotPos_FixedPz_MidPt (%3.1f to %3.1f)Gev/c", HistoDalitzPtRangeMin_MidPt, HistoDalitzPtRangeMax_MidPt), HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz, HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz);
          fHistoDalitzPlotPosFixedPzNDM_MidPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-}} (GeV/c^{2})");
          fHistoDalitzPlotPosFixedPzNDM_MidPt[iCut]->GetYaxis()->SetTitle(Form("M_{#pi^{+} %s} (GeV/c^{2}",NameNDMLatex.Data()));
          fHistoDalitzPlotPosFixedPzNDM_MidPt[iCut]->Sumw2();
          fESDList[iCut]->Add(fHistoDalitzPlotPosFixedPzNDM_MidPt[iCut]);
          fHistoDalitzPlotNegFixedPzNDM_MidPt[iCut]          = new TH2F("ESD_DalitzPlotNeg_FixedPz_MidPt", Form("ESD_DalitzPlotNeg_FixedPz_MidPt (%3.1f to %3.1f)Gev/c", HistoDalitzPtRangeMin_MidPt, HistoDalitzPtRangeMax_MidPt), HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz, HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz);
          fHistoDalitzPlotNegFixedPzNDM_MidPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-}} (GeV/c^{2})");
          fHistoDalitzPlotNegFixedPzNDM_MidPt[iCut]->GetYaxis()->SetTitle(Form("M_{#pi^{-} %s} (GeV/c^{2}",NameNDMLatex.Data()));
          fHistoDalitzPlotNegFixedPzNDM_MidPt[iCut]->Sumw2();
          fESDList[iCut]->Add(fHistoDalitzPlotNegFixedPzNDM_MidPt[iCut]);
          fHistoDalitzPlotPosSubNDM_MidPt[iCut]              = new TH2F("ESD_DalitzPlotPos_Sub_MidPt", Form("ESD_DalitzPlotPos_Sub_MidPt (%3.1f to %3.1f)Gev/c", HistoDalitzPtRangeMin_MidPt, HistoDalitzPtRangeMax_MidPt), HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz, HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz);
          fHistoDalitzPlotPosSubNDM_MidPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-}} (GeV/c^{2})");
          fHistoDalitzPlotPosSubNDM_MidPt[iCut]->GetYaxis()->SetTitle(Form("M_{#pi^{+} %s} (GeV/c^{2}",NameNDMLatex.Data()));
          fHistoDalitzPlotPosSubNDM_MidPt[iCut]->Sumw2();
          fESDList[iCut]->Add(fHistoDalitzPlotPosSubNDM_MidPt[iCut]);
          fHistoDalitzPlotNegSubNDM_MidPt[iCut]              = new TH2F("ESD_DalitzPlotNeg_Sub_MidPt", Form("ESD_DalitzPlotNeg_Sub_MidPt (%3.1f to %3.1f)Gev/c", HistoDalitzPtRangeMin_MidPt, HistoDalitzPtRangeMax_MidPt), HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz, HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz);
          fHistoDalitzPlotNegSubNDM_MidPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-}} (GeV/c^{2})");
          fHistoDalitzPlotNegSubNDM_MidPt[iCut]->GetYaxis()->SetTitle(Form("M_{#pi^{-} %s} (GeV/c^{2}",NameNDMLatex.Data()));
          fHistoDalitzPlotNegSubNDM_MidPt[iCut]->Sumw2();
          fESDList[iCut]->Add(fHistoDalitzPlotNegSubNDM_MidPt[iCut]);
        }
        //---------------------------------------
        //Dalitz High Pt
        if (enableDalitzHighPt){
          fHistoDalitzPlotPosFixedPzNDM_HighPt[iCut]          = new TH2F("ESD_DalitzPlotPos_FixedPz_HighPt", Form("ESD_DalitzPlotPos_FixedPz_HighPt (%3.1f to %3.1f)Gev/c", HistoDalitzPtRangeMin_HighPt, HistoDalitzPtRangeMax_HighPt), HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz, HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz);
          fHistoDalitzPlotPosFixedPzNDM_HighPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-}} (GeV/c^{2})");
          fHistoDalitzPlotPosFixedPzNDM_HighPt[iCut]->GetYaxis()->SetTitle(Form("M_{#pi^{+} %s} (GeV/c^{2}",NameNDMLatex.Data()));
          fHistoDalitzPlotPosFixedPzNDM_HighPt[iCut]->Sumw2();
          fESDList[iCut]->Add(fHistoDalitzPlotPosFixedPzNDM_HighPt[iCut]);
          fHistoDalitzPlotNegFixedPzNDM_HighPt[iCut]          = new TH2F("ESD_DalitzPlotNeg_FixedPz_HighPt", Form("ESD_DalitzPlotNeg_FixedPz_HighPt (%3.1f to %3.1f)Gev/c", HistoDalitzPtRangeMin_HighPt, HistoDalitzPtRangeMax_HighPt), HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz, HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz);
          fHistoDalitzPlotNegFixedPzNDM_HighPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-}} (GeV/c^{2})");
          fHistoDalitzPlotNegFixedPzNDM_HighPt[iCut]->GetYaxis()->SetTitle(Form("M_{#pi^{-} %s} (GeV/c^{2}",NameNDMLatex.Data()));
          fHistoDalitzPlotNegFixedPzNDM_HighPt[iCut]->Sumw2();
          fESDList[iCut]->Add(fHistoDalitzPlotNegFixedPzNDM_HighPt[iCut]);
          fHistoDalitzPlotPosSubNDM_HighPt[iCut]              = new TH2F("ESD_DalitzPlotPos_Sub_HighPt", Form("ESD_DalitzPlotPos_Sub_HighPt (%3.1f to %3.1f)Gev/c", HistoDalitzPtRangeMin_HighPt, HistoDalitzPtRangeMax_HighPt), HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz, HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz);
          fHistoDalitzPlotPosSubNDM_HighPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-}} (GeV/c^{2})");
          fHistoDalitzPlotPosSubNDM_HighPt[iCut]->GetYaxis()->SetTitle(Form("M_{#pi^{+} %s} (GeV/c^{2}",NameNDMLatex.Data()));
          fHistoDalitzPlotPosSubNDM_HighPt[iCut]->Sumw2();
          fESDList[iCut]->Add(fHistoDalitzPlotPosSubNDM_HighPt[iCut]);
          fHistoDalitzPlotNegSubNDM_HighPt[iCut]              = new TH2F("ESD_DalitzPlotNeg_Sub_HighPt", Form("ESD_DalitzPlotNeg_Sub_HighPt (%3.1f to %3.1f)Gev/c", HistoDalitzPtRangeMin_HighPt, HistoDalitzPtRangeMax_HighPt), HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz, HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz);
          fHistoDalitzPlotNegSubNDM_HighPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-}} (GeV/c^{2})");
          fHistoDalitzPlotNegSubNDM_HighPt[iCut]->GetYaxis()->SetTitle(Form("M_{#pi^{-} %s} (GeV/c^{2}",NameNDMLatex.Data()));
          fHistoDalitzPlotNegSubNDM_HighPt[iCut]->Sumw2();
          fESDList[iCut]->Add(fHistoDalitzPlotNegSubNDM_HighPt[iCut]);
        }
        //---------------------------------------
        //End of Dalitz
      }
      fHistoGammaGammaInvMassPt[iCut]               = new TH2F("ESD_GammaGamma_InvMass_Pt","ESD_GammaGamma_InvMass_Pt",HistoNMassBinsDecayMeson,HistoMassRangeNDM[0],HistoMassRangeNDM[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
      fHistoGammaGammaInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#gamma #gamma} (GeV/c^{2})");
      fHistoGammaGammaInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
      fHistoGammaGammaInvMassPt[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoGammaGammaInvMassPt[iCut]);

      fHistoGammaGammaInvMassPtBeforeCuts[iCut]               = new TH2F("ESD_GammaGamma_InvMass_Pt_Before_Cuts","ESD_GammaGamma_InvMass_Pt_Before_Cuts",HistoNMassBinsDecayMeson,HistoMassRangeNDM[0],HistoMassRangeNDM[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
      fHistoGammaGammaInvMassPtBeforeCuts[iCut]->GetXaxis()->SetTitle("M_{#gamma #gamma} (GeV/c^{2})");
      fHistoGammaGammaInvMassPtBeforeCuts[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
      fHistoGammaGammaInvMassPtBeforeCuts[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoGammaGammaInvMassPtBeforeCuts[iCut]);
    }
    fHistoMotherInvMassPt[iCut]                   = new TH2F("ESD_Mother_InvMass_Pt","ESD_Mother_InvMass_Pt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
    fHistoMotherInvMassPt[iCut]->GetXaxis()->SetTitle(Form("M_{#pi^{+} #pi^{-} %s} (GeV/c^{2})",NameNDMLatex.Data()));
    fHistoMotherInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    fHistoMotherInvMassPt[iCut]->Sumw2();
    fESDList[iCut]->Add(fHistoMotherInvMassPt[iCut]);
    fHistoMotherInvMassPtRejectedKinematic[iCut]  = new TH2F("ESD_Mother_InvMass_Pt_KinematicRejected","ESD_Mother_InvMass_Pt_KinematicRejected",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
    fHistoMotherInvMassPtRejectedKinematic[iCut]->GetXaxis()->SetTitle(Form("M_{#pi^{+} #pi^{-} %s} (GeV/c^{2})",NameNDMLatex.Data()));
    fHistoMotherInvMassPtRejectedKinematic[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    fHistoMotherInvMassPtRejectedKinematic[iCut]->Sumw2();
    fESDList[iCut]->Add(fHistoMotherInvMassPtRejectedKinematic[iCut]);

    fHistoBackInvMassPt[iCut] = new TH2F("ESD_Background_InvMass_Pt","ESD_Background_InvMass_Pt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
    fHistoBackInvMassPt[iCut]->Sumw2();
    fHistoBackInvMassPt[iCut]->GetXaxis()->SetTitle(Form("M_{#pi^{+} #pi^{-} %s} (GeV/c^{2})",NameNDMLatex.Data()));
    fHistoBackInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    if(!(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseLikeSignMixing())&&(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->DoBGCalculation())){
       fESDList[iCut]->Add(fHistoBackInvMassPt[iCut]);
    }

    fHistoMotherLikeSignBackInvMassPt[iCut]  = new TH2F("ESD_Background_LikeSign_InvMass_Pt","ESD_Background_LikeSign_InvMass_Pt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
    fHistoMotherLikeSignBackInvMassPt[iCut]->GetXaxis()->SetTitle(Form("M_{#pi^{#pm} #pi^{#pm} %s} (GeV/c^{2})",NameNDMLatex.Data()));
    fHistoMotherLikeSignBackInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    fHistoMotherLikeSignBackInvMassPt[iCut]->Sumw2();
    if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseLikeSignMixing()){
      fESDList[iCut]->Add(fHistoMotherLikeSignBackInvMassPt[iCut]);
    }

    if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetBackgroundMode() == 7){
      fHistoSwappingGammaGammaInvMassPt[iCut]  = new TH2F("ESD_Swapping_GammaGamma_InvMass_Pt","ESD_Swapping_GammaGamma_InvMass_Pt",HistoNMassBinsDecayMeson,HistoMassRangeNDM[0],HistoMassRangeNDM[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
      fHistoSwappingGammaGammaInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#gamma #gamma} (GeV/c^{2})");
      fHistoSwappingGammaGammaInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
      fHistoSwappingGammaGammaInvMassPt[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoSwappingGammaGammaInvMassPt[iCut]);
    }

    fHistoMotherInvMassSubNDM[iCut]                       = new TH2F("ESD_InvMass_Mother_Sub_InvMass_Neutral_Pt","ESD_InvMass_Mother_Sub_InvMass_Neutral_Pt",HistoNMassBinsSub,HistoMassRangeSub[0],HistoMassRangeSub[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
    fHistoMotherInvMassSubNDM[iCut]->GetXaxis()->SetTitle(Form("M_{#pi^{+} #pi^{-} %s} - (M_{%s}-M_{%s},PDG}) (GeV/c^{2})",NameNDMLatex.Data(),NameNDMLatex.Data(),NameNDMLatex.Data()));
    fHistoMotherInvMassSubNDM[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    fHistoMotherInvMassSubNDM[iCut]->Sumw2();
    fESDList[iCut]->Add(fHistoMotherInvMassSubNDM[iCut]);

    fHistoBackInvMassPtSubNDM[iCut]   = new TH2F("ESD_Background_InvMass_Sub_InvMass_Neutral_Pt","ESD_Background_InvMass_Sub_InvMass_Neutral_Pt",
                                                                     HistoNMassBinsSub,HistoMassRangeSub[0],HistoMassRangeSub[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
    fHistoBackInvMassPtSubNDM[iCut]->GetXaxis()->SetTitle(Form("M_{#pi^{+} #pi^{-} %s} - (M_{%s}-M_{%s},PDG}) (GeV/c^{2})",NameNDMLatex.Data(),NameNDMLatex.Data(),NameNDMLatex.Data()));
    fHistoBackInvMassPtSubNDM[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    fHistoBackInvMassPtSubNDM[iCut]->Sumw2();
    if(!(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseLikeSignMixing())&&(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->DoBGCalculation())){
      fESDList[iCut]->Add(fHistoBackInvMassPtSubNDM[iCut]);
    }
    fHistoMotherLikeSignBackInvMassSubNDMPt[iCut]    = new TH2F("ESD_Background_LikeSign_InvMass_Sub_InvMass_Neutral_Pt","ESD_Background_LikeSign_InvMass_Sub_InvMass_Neutral_Pt",
                                                                HistoNMassBinsSub,HistoMassRangeSub[0],HistoMassRangeSub[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
    fHistoMotherLikeSignBackInvMassSubNDMPt[iCut]->GetXaxis()->SetTitle(Form("M_{#pi^{#pm} #pi^{#pm} %s} - M_{%s} (GeV/c^{2})",NameNDMLatex.Data(),NameNDMLatex.Data()));
    fHistoMotherLikeSignBackInvMassSubNDMPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    fHistoMotherLikeSignBackInvMassSubNDMPt[iCut]->Sumw2();
    if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseLikeSignMixing()){
      fESDList[iCut]->Add(fHistoMotherLikeSignBackInvMassSubNDMPt[iCut]);
    }
    fHistoMotherInvMassFixedPzNDM[iCut]                     = new TH2F("ESD_InvMass_Mother_FixedPz_Neutral_Pt","ESD_Mother_InvMass_FixedPz_Neutral_Pt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
    fHistoMotherInvMassFixedPzNDM[iCut]->GetXaxis()->SetTitle(Form("M_{#pi^{+} #pi^{-} %s} (GeV/c^{2})",NameNDMLatex.Data()));
    fHistoMotherInvMassFixedPzNDM[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    fHistoMotherInvMassFixedPzNDM[iCut]->Sumw2();
    fESDList[iCut]->Add(fHistoMotherInvMassFixedPzNDM[iCut]);

    fHistoBackInvMassPtFixedPzNDM[iCut] = new TH2F("ESD_Background_InvMass_FixedPz_Neutral_Pt","ESD_Background_InvMass_FixedPz_Neutral_Pt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
    fHistoBackInvMassPtFixedPzNDM[iCut]->GetXaxis()->SetTitle(Form("M_{#pi^{+} #pi^{-} %s} (GeV/c^{2})",NameNDMLatex.Data()));
    fHistoBackInvMassPtFixedPzNDM[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    fHistoBackInvMassPtFixedPzNDM[iCut]->Sumw2();
    if(!(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseLikeSignMixing()) && (((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->DoBGCalculation())){
      fESDList[iCut]->Add(fHistoBackInvMassPtFixedPzNDM[iCut]);
    }
    fHistoMotherLikeSignBackInvMassFixedPzNDMPt[iCut]  = new TH2F("ESD_Background_LikeSign_InvMass_FixedPz_Neutral_Pt","ESD_Background_LikeSign_InvMass_FixedPz_Neutral_Pt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
    fHistoMotherLikeSignBackInvMassFixedPzNDMPt[iCut]->GetXaxis()->SetTitle(Form("M_{#pi^{#pm} #pi^{#pm} %s} (GeV/c^{2})",NameNDMLatex.Data()));
    fHistoMotherLikeSignBackInvMassFixedPzNDMPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    fHistoMotherLikeSignBackInvMassFixedPzNDMPt[iCut]->Sumw2();
    if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseLikeSignMixing()){
      fESDList[iCut]->Add(fHistoMotherLikeSignBackInvMassFixedPzNDMPt[iCut]);
    }
    if(!fDoLightOutput){
      fHistoAngleHNMesonPiPlPiMi[iCut]      = new TH2F(Form("ESD_Mother_Angle%sNegPionsPosPions_Pt",NameNeutralMesonAnalyzed.Data()),Form("ESD_Mother_Angle%sNegPionsPosPions_Pt",NameNeutralMesonAnalyzed.Data()),HistoNPtBins,HistoPtRange[0],HistoPtRange[1],360,0,TMath::Pi());
      fHistoAngleHNMesonPiPlPiMi[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      fHistoAngleHNMesonPiPlPiMi[iCut]->GetYaxis()->SetTitle("#angle (meson,#pi^{+}#pi^{-})");
      fHistoAngleHNMesonPiPlPiMi[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoAngleHNMesonPiPlPiMi[iCut]);
      fHistoAngleHNMesonPiMi[iCut]          = new TH2F(Form("ESD_Mother_Angle%sNegPions_Pt",NameNeutralMesonAnalyzed.Data()),Form("ESD_Mother_Angle%sNegPions_Pt",NameNeutralMesonAnalyzed.Data()),HistoNPtBins,HistoPtRange[0],HistoPtRange[1],360,0,TMath::Pi());
      fHistoAngleHNMesonPiMi[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      fHistoAngleHNMesonPiMi[iCut]->GetYaxis()->SetTitle("#angle (meson,#pi^{-})");
      fHistoAngleHNMesonPiMi[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoAngleHNMesonPiMi[iCut]);
      fHistoAngleHNMesonPiPl[iCut]          = new TH2F(Form("ESD_Mother_Angle%sPosPions_Pt",NameNeutralMesonAnalyzed.Data()),Form("ESD_Mother_Angle%sPosPions_Pt",NameNeutralMesonAnalyzed.Data()),HistoNPtBins,HistoPtRange[0],HistoPtRange[1],360,0,TMath::Pi());
      fHistoAngleHNMesonPiPl[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      fHistoAngleHNMesonPiPl[iCut]->GetYaxis()->SetTitle("#angle (meson,#pi^{+})");
      fHistoAngleHNMesonPiPl[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoAngleHNMesonPiPl[iCut]);
      fHistoAngleHNMesonNDM[iCut]        = new TH2F(Form("ESD_Mother_Angle%s%s_Pt",NameNeutralMesonAnalyzed.Data(),NameNDM.Data()),Form("ESD_Mother_Angle%s%s_Pt",NameNeutralMesonAnalyzed.Data(),NameNDM.Data()),HistoNPtBins,HistoPtRange[0],HistoPtRange[1],360,0,TMath::Pi());
      fHistoAngleHNMesonNDM[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      fHistoAngleHNMesonNDM[iCut]->GetYaxis()->SetTitle(Form("#angle (meson,%s)",NameNDMLatex.Data()));
      fHistoAngleHNMesonNDM[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoAngleHNMesonNDM[iCut]);
      fHistoAnglePiPlNDM[iCut]         = new TH2F(Form("ESD_Mother_AnglePosPions%s_Pt",NameNDM.Data()),Form("ESD_Mother_AnglePosPions%s_Pt",NameNDM.Data()),HistoNPtBins,HistoPtRange[0],HistoPtRange[1],360,0,TMath::Pi());
      fHistoAnglePiPlNDM[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      fHistoAnglePiPlNDM[iCut]->GetYaxis()->SetTitle(Form("#angle (#pi^{+},%s)",NameNDMLatex.Data()));
      fHistoAnglePiPlNDM[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoAnglePiPlNDM[iCut]);
      fHistoAnglePiPlPiMi[iCut]           = new TH2F("ESD_Mother_AnglePosPionsNegPions_Pt","ESD_Mother_AnglePosPionsNegPions_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1],360,0,TMath::Pi());
      fHistoAnglePiPlPiMi[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      fHistoAnglePiPlPiMi[iCut]->GetYaxis()->SetTitle("#angle (#pi^{+},#pi^{-})");
      fHistoAnglePiPlPiMi[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoAnglePiPlPiMi[iCut]);
      fHistoAngleNDMPiMi[iCut]         = new TH2F(Form("ESD_Mother_Angle%sNegPions_Pt",NameNDM.Data()),Form("ESD_Mother_Angle%sNegPions_Pt",NameNDM.Data()),HistoNPtBins,HistoPtRange[0],HistoPtRange[1],360,0,TMath::Pi());
      fHistoAngleNDMPiMi[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      fHistoAngleNDMPiMi[iCut]->GetYaxis()->SetTitle(Form("#angle (%s,#pi^{-})",NameNDMLatex.Data()));
      fHistoAngleNDMPiMi[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoAngleNDMPiMi[iCut]);
      fHistoAngleSum[iCut]                = new TH2F("ESD_Mother_AngleSum_Pt","ESD_Mother_AngleSum_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1],720,0,2*TMath::Pi());
      fHistoAngleSum[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      fHistoAngleSum[iCut]->GetYaxis()->SetTitle("#sum #angle");
      fHistoAngleSum[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoAngleSum[iCut]);
    }
    if ( fDoMesonQA>0 && (!fDoLightOutput) ) {
      TAxis *AxisAfter        = fHistoPionTPCdEdxNSigma[iCut]->GetXaxis();
      Int_t bins              = AxisAfter->GetNbins();
      Double_t from           = AxisAfter->GetXmin();
      Double_t to             = AxisAfter->GetXmax();
      Double_t *newBins       = new Double_t[bins+1];
      newBins[0]              = from;
      Double_t factor         = TMath::Power(to/from, 1./bins);
      for(Int_t i=1; i<=bins; ++i) newBins[i] = factor * newBins[i-1];

      AxisAfter->Set(bins, newBins);
      AxisAfter               = fHistoPionTPCdEdx[iCut]->GetXaxis();
      AxisAfter->Set(bins, newBins);
      delete [] newBins;
    }

    fCutFolder[iCut]->Add(fESDList[iCut]);

  }

  if( fIsMC ){
    // MC Histogramms
    fMCList                                 = new TList*[fnCuts];
    // True Histogramms
    fTrueList                               = new TList*[fnCuts];
    if(!fDoLightOutput){
      if (fNDMRecoMode < 2){
        fHistoTrueConvGammaPt                 = new TH1F*[fnCuts];
        fHistoDoubleCountTrueConvGammaRPt     = new TH2F*[fnCuts];
        fHistoTrueConvGammaFromNeutralMesonPt = new TH1F*[fnCuts];
      }
      if (fNDMRecoMode > 0){
        fHistoTrueClusterGammaPt                  = new TH1F*[fnCuts];
        fHistoTrueClusterGammaFromNeutralMesonPt  = new TH1F*[fnCuts];
      }
      fHistoTruePosPionPt                     = new TH1F*[fnCuts];
      fHistoTrueNegPionPt                     = new TH1F*[fnCuts];
      fHistoTruePosPionFromNeutralMesonPt     = new TH1F*[fnCuts];
      fHistoTrueNegPionFromNeutralMesonPt     = new TH1F*[fnCuts];


      fHistoMCAllGammaPt                      = new TH1F*[fnCuts];
      if (fNDMRecoMode < 2){
        fHistoMCConvGammaPt                   = new TH1F*[fnCuts];
      }
      fHistoMCAllPosPionsPt                   = new TH1F*[fnCuts];
      fHistoMCAllNegPionsPt                   = new TH1F*[fnCuts];
      fHistoMCGammaFromNeutralMesonPt         = new TH1F*[fnCuts];
      fHistoMCPosPionsFromNeutralMesonPt      = new TH1F*[fnCuts];
      fHistoMCNegPionsFromNeutralMesonPt      = new TH1F*[fnCuts];

      if(fDoMesonQA > 0){
        fHistoMCAllMesonPt                      = new TH1F*[fnCuts];
        fHistoMCAllMesonEta                     = new TH1F*[fnCuts];
        fHistoMCAllMesonPhi                     = new TH1F*[fnCuts];
        fHistoMCMesonFromNeutralMesonPt         = new TH1F*[fnCuts];
        fHistoMCMesonFromNeutralMesonEta        = new TH1F*[fnCuts];
        fHistoMCMesonFromNeutralMesonPhi        = new TH1F*[fnCuts];
        fHistoMCAllPosPionsEta                        = new TH1F*[fnCuts];
        fHistoMCAllPosPionsPhi                        = new TH1F*[fnCuts];
        fHistoMCAllNegPionsEta                        = new TH1F*[fnCuts];
        fHistoMCAllNegPionsPhi                        = new TH1F*[fnCuts];
        fHistoMCPosPionsFromNeutralMesonEta           = new TH1F*[fnCuts];
        fHistoMCPosPionsFromNeutralMesonPhi           = new TH1F*[fnCuts];
        fHistoMCNegPionsFromNeutralMesonEta           = new TH1F*[fnCuts];
        fHistoMCNegPionsFromNeutralMesonPhi           = new TH1F*[fnCuts];
        fHistoMCHeavyAllPt                            = new TH1F*[fnCuts];
        fHistoMCHeavyAllEta                           = new TH1F*[fnCuts];
        fHistoMCHeavyAllPhi                           = new TH1F*[fnCuts];
        fHistoMCHeavyChannelPt                        = new TH1F*[fnCuts];
        fHistoMCHeavyChannelEta                       = new TH1F*[fnCuts];
        fHistoMCHeavyChannelPhi                       = new TH1F*[fnCuts];
        fHistMCChannelNDMFromHeavyPt                  = new TH1F*[fnCuts];
        fHistMCChannelNDMFromHeavyEta                 = new TH1F*[fnCuts];
        fHistMCChannelNDMFromHeavyPhi                 = new TH1F*[fnCuts];
        fHistMCChannelPiPlusFromHeavyPt               = new TH1F*[fnCuts];
        fHistMCChannelPiPlusFromHeavyEta              = new TH1F*[fnCuts];
        fHistMCChannelPiPlusFromHeavyPhi              = new TH1F*[fnCuts];
        fHistMCChannelPiMinusFromHeavyPt              = new TH1F*[fnCuts];
        fHistMCChannelPiMinusFromHeavyEta             = new TH1F*[fnCuts];
        fHistMCChannelPiPMinusFromHeavyPhi            = new TH1F*[fnCuts];
        fHistMCChannelNDMPtHeavyPt                    = new TH2F*[fnCuts];
        fHistMCChannelPiPlusPtHeavyPt                 = new TH2F*[fnCuts];
        fHistMCChannelPiMinusPtHeavyPt                = new TH2F*[fnCuts];
        fHistoMCHeavyReconstructiblePt                = new TH1F*[fnCuts];
        fHistoMCHeavyReconstructibleEta               = new TH1F*[fnCuts];
        fHistoMCHeavyReconstructiblePhi               = new TH1F*[fnCuts];
        fHistMCReconstructibleNDMFromHeavyPt          = new TH1F*[fnCuts];
        fHistMCReconstructibleNDMFromHeavyEta         = new TH1F*[fnCuts];
        fHistMCReconstructibleNDMFromHeavyPhi         = new TH1F*[fnCuts];
        fHistMCReconstructiblePiPlusFromHeavyPt       = new TH1F*[fnCuts];
        fHistMCReconstructiblePiPlusFromHeavyEta      = new TH1F*[fnCuts];
        fHistMCReconstructiblePiPlusFromHeavyPhi      = new TH1F*[fnCuts];
        fHistMCReconstructiblePiMinusFromHeavyPt      = new TH1F*[fnCuts];
        fHistMCReconstructiblePiMinusFromHeavyEta     = new TH1F*[fnCuts];
        fHistMCReconstructiblePiPMinusFromHeavyPhi    = new TH1F*[fnCuts];
        fHistMCReconstructibleNDMPtHeavyPt            = new TH2F*[fnCuts];
        fHistMCReconstructiblePiPlusPtHeavyPt         = new TH2F*[fnCuts];
        fHistMCReconstructiblePiMinusPtHeavyPt        = new TH2F*[fnCuts];
      }
    }
    fHistoMCHNMPiPlPiMiNDMPt                     = new TH1F*[fnCuts];
    fHistoMCHNMPiPlPiMiNDMInAccPt                = new TH1F*[fnCuts];
    if(!fDoLightOutput){
      fHistoMCHNMInAccVsNDMPt                      = new TH2F*[fnCuts];
      fHistoMCHNMPiPlPiMiNDMEta                    = new TH1F*[fnCuts];
      fHistoMCHNMPiPlPiMiNDMPhi                    = new TH1F*[fnCuts];
      fHistoDoubleCountTruePi0InvMassPt               = new TH2F*[fnCuts];
      fHistoDoubleCountTrueHNMInvMassPt               = new TH2F*[fnCuts];
      if(fDoMesonQA > 1){
        fHistoMCHNMPiPlPiMiNDMEtavsPt                 = new TH2F*[fnCuts];
      }
    }
    fHistoTrueMotherPiPlPiMiNDMInvMassPt              = new TH2F*[fnCuts];
    fHistoTrueMotherPiPlPiMiNDMInvMassPtSubNDM        = new TH2F*[fnCuts];
    fHistoTrueMotherPiPlPiMiNDMInvMassPtFixedPzNDM    = new TH2F*[fnCuts];
    fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromDifferent  = new TH2F*[fnCuts];
    fHistoTruePiPlPiMiNDMCombinatoricalInvMassPt        = new TH2F*[fnCuts];
    fHistoTruePiPlPiMiNDMContaminationInvMassPt         = new TH2F*[fnCuts];
    if(!fDoLightOutput){
      fHistoTrueMotherGammaGammaInvMassPt             = new TH2F*[fnCuts];
      fHistoTrueMotherGammaGammaFromHNMInvMassPt      = new TH2F*[fnCuts];
      fHistoTrueAngleSum                              = new TH2F*[fnCuts];
      fHistoTrueHNMesonPtvsNDMPt                      = new TH2F*[fnCuts];
    }
    if(!fDoLightOutput){
      if (fDoMesonQA>0){
        fHistoTruePionPionInvMassPt                               = new TH2F*[fnCuts];
        fHistoTruePionPionFromSameMotherInvMassPt                 = new TH2F*[fnCuts];
        fHistoTruePionPionFromHNMInvMassPt                        = new TH2F*[fnCuts];

        fHistoTrueMesonFlags                                      = new TH1F*[fnCuts];

        fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromEtaOmega        = new TH2F*[fnCuts];
        fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromRho             = new TH2F*[fnCuts];
        fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromK0s             = new TH2F*[fnCuts];
        fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromK0l             = new TH2F*[fnCuts];
        fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromEtaPrime        = new TH2F*[fnCuts];
        fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromOther           = new TH2F*[fnCuts];

        fHistoTruePiPlPiMiSameMotherFromEtaInvMassPt              = new TH2F*[fnCuts];
        fHistoTruePiPlPiMiSameMotherFromOmegaInvMassPt            = new TH2F*[fnCuts];
        fHistoTruePiPlPiMiSameMotherFromRhoInvMassPt              = new TH2F*[fnCuts];
        fHistoTruePiPlPiMiSameMotherFromEtaPrimeInvMassPt         = new TH2F*[fnCuts];
        fHistoTruePiPlPiMiSameMotherFromK0sInvMassPt              = new TH2F*[fnCuts];
        fHistoTruePiPlPiMiSameMotherFromK0lInvMassPt              = new TH2F*[fnCuts];
        fHistoTruePiPlPiMiSameMotherFromOtherlInvMassPt           = new TH2F*[fnCuts];
        fHistoTruePiMiPiZeroSameMotherFromEtaInvMassPt            = new TH2F*[fnCuts];
        fHistoTruePiMiPiZeroSameMotherFromOmegaInvMassPt          = new TH2F*[fnCuts];
        fHistoTruePiMiPiZeroSameMotherFromRhoInvMassPt            = new TH2F*[fnCuts];
        fHistoTruePiMiPiZeroSameMotherFromK0lInvMassPt            = new TH2F*[fnCuts];
        fHistoTruePiMiPiZeroSameMotherFromOtherlInvMassPt         = new TH2F*[fnCuts];
        fHistoTruePiPlPiZeroSameMotherFromEtaInvMassPt            = new TH2F*[fnCuts];
        fHistoTruePiPlPiZeroSameMotherFromOmegaInvMassPt          = new TH2F*[fnCuts];
        fHistoTruePiPlPiZeroSameMotherFromRhoInvMassPt            = new TH2F*[fnCuts];
        fHistoTruePiPlPiZeroSameMotherFromK0lInvMassPt            = new TH2F*[fnCuts];
        fHistoTruePiPlPiZeroSameMotherFromOtherInvMassPt          = new TH2F*[fnCuts];
        fHistoTruePiPlPiMiNDMPureCombinatoricalInvMassPt          = new TH2F*[fnCuts];

        fHistoTruePiPlPiMiNDMContamination_Pi0_InvMassPt          = new TH2F*[fnCuts];
        fHistoTruePiPlPiMiNDMContamination_PiPl_InvMassPt         = new TH2F*[fnCuts];
        fHistoTruePiPlPiMiNDMContamination_PiMi_InvMassPt         = new TH2F*[fnCuts];
        fHistoTruePiPlPiMiNDMContamination_Crosscheck_InvMassPt   = new TH2F*[fnCuts];
        fHistoTruePiPlPiMiNDMContamination_multipel_InvMassPt     = new TH2F*[fnCuts];
        //Dalitz All Pt
        if (enableDalitzAllPt){
          fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM  = new TH2F*[fnCuts];
          fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM  = new TH2F*[fnCuts];
          fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM    = new TH2F*[fnCuts];
          fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM    = new TH2F*[fnCuts];
        }
        //Dalitz Low Pt
        if (enableDalitzLowPt){
          fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM_LowPt    = new TH2F*[fnCuts];
          fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM_LowPt    = new TH2F*[fnCuts];
          fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM_LowPt        = new TH2F*[fnCuts];
          fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM_LowPt        = new TH2F*[fnCuts];
        }
        //Dalitz Mid Pt
        if (enableDalitzMidPt){
          fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM_MidPt    = new TH2F*[fnCuts];
          fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM_MidPt    = new TH2F*[fnCuts];
          fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM_MidPt        = new TH2F*[fnCuts];
          fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM_MidPt        = new TH2F*[fnCuts];
        }
        if (enableDalitzHighPt){
        //Dalitz High Pt
          fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM_HighPt   = new TH2F*[fnCuts];
          fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM_HighPt   = new TH2F*[fnCuts];
          fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM_HighPt       = new TH2F*[fnCuts];
          fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM_HighPt       = new TH2F*[fnCuts];
        }

        fHistoTruevParticleChi2PerNDF = new TH1F*[fnCuts];
        fHistoTruevParticleFromSameMotherChi2PerNDF = new TH1F*[fnCuts];
        fHistoTruevParticleFromHNMChi2PerNDF = new TH1F*[fnCuts];
        fHistoTruevParticledS = new TH1F*[fnCuts];
        fHistoTruevParticleFromSameMotherdS = new TH1F*[fnCuts];
        fHistoTruevParticleFromHNMdS = new TH1F*[fnCuts];


        if (fDoMesonQA>1){
          fTrueTreeList                                           = new TList*[fnCuts];
          fTreePiPiSameMother                                     = new TTree*[fnCuts];
          fTreePiPiPiSameMother                                   = new TTree*[fnCuts];
          fTreeEventInfoHNM                                       = new TTree*[fnCuts];
        }
      }
    }

    for(Int_t iCut = 0; iCut<fnCuts;iCut++){
      TString cutstringEvent            = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
      TString cutstringPion             = ((AliPrimaryPionCuts*)fPionCutArray->At(iCut))->GetCutNumber();
      TString cutstringConvGamma        = "";
      if (fNDMRecoMode < 2)
        cutstringConvGamma              = ((AliConversionPhotonCuts*)fGammaCutArray->At(iCut))->GetCutNumber();
      TString cutstringCaloGamma        = "";
      if (fNDMRecoMode > 0)
        cutstringCaloGamma              = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
      TString cutstringNeutralPion      = ((AliConversionMesonCuts*)fNeutralDecayMesonCutArray->At(iCut))->GetCutNumber();
      TString cutstringMeson            = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();

      TString fullCutString             = "";
      if (fNDMRecoMode == 0)
        fullCutString                   = Form("%i_%s_%s_%s_%s_%s",fNDMRecoMode,cutstringEvent.Data(),cutstringConvGamma.Data(),cutstringNeutralPion.Data(), cutstringPion.Data(),
                                               cutstringMeson.Data());
      else if (fNDMRecoMode == 1)
        fullCutString                   = Form("%i_%s_%s_%s_%s_%s_%s",fNDMRecoMode,cutstringEvent.Data(),cutstringConvGamma.Data(),cutstringCaloGamma.Data(), cutstringNeutralPion.Data(),
                                               cutstringPion.Data(), cutstringMeson.Data());
      else if (fNDMRecoMode == 2)
        fullCutString                   = Form("%i_%s_%s_%s_%s_%s",fNDMRecoMode,cutstringEvent.Data(), cutstringCaloGamma.Data(), cutstringNeutralPion.Data(), cutstringPion.Data(),
                                               cutstringMeson.Data());
      TString nameMCList                = Form("%s MC histograms", fullCutString.Data());
      TString nameTrueRecList           = Form("%s True histograms", fullCutString.Data());
      TString nameTrueRecTTreeList      = Form("%s True TTrees", fullCutString.Data());

      fMCList[iCut]                     = new TList();
      fMCList[iCut]->SetName(nameMCList.Data());
      fMCList[iCut]->SetOwner(kTRUE);
      fCutFolder[iCut]->Add(fMCList[iCut]);

      if(!fDoLightOutput){
        fHistoMCAllGammaPt[iCut]          = new TH1F("MC_AllGamma_Pt","MC_AllGamma_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
        fHistoMCAllGammaPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoMCAllGammaPt[iCut]->GetYaxis()->SetTitle("N_{#gamma}");
        fHistoMCAllGammaPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCAllGammaPt[iCut]);
        if (fNDMRecoMode < 2){
          fHistoMCConvGammaPt[iCut]       = new TH1F("MC_ConvGamma_Pt","MC_ConvGamma_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoMCConvGammaPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
          fHistoMCConvGammaPt[iCut]->GetYaxis()->SetTitle("N_{#gamma,conv}");
          fHistoMCConvGammaPt[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCConvGammaPt[iCut]);
        }
        fHistoMCGammaFromNeutralMesonPt[iCut]     = new TH1F("MC_GammaFromNeutralMeson_Pt","MC_GammaFromNeutralMeson_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
        fHistoMCGammaFromNeutralMesonPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoMCGammaFromNeutralMesonPt[iCut]->GetYaxis()->SetTitle("N_{#gamma}");
        fHistoMCGammaFromNeutralMesonPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCGammaFromNeutralMesonPt[iCut]);
        fHistoMCAllPosPionsPt[iCut]               = new TH1F("MC_AllPosPions_Pt","MC_AllPosPions_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
        fHistoMCAllPosPionsPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoMCAllPosPionsPt[iCut]->GetYaxis()->SetTitle("N_{#pi^{+}}");
        fHistoMCAllPosPionsPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCAllPosPionsPt[iCut]);
        fHistoMCAllNegPionsPt[iCut]               = new TH1F("MC_AllNegPions_Pt","MC_AllNegPions_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
        fHistoMCAllNegPionsPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoMCAllNegPionsPt[iCut]->GetYaxis()->SetTitle("N_{#pi^{-}}");
        fHistoMCAllNegPionsPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCAllNegPionsPt[iCut]);
        fHistoMCNegPionsFromNeutralMesonPt[iCut]  = new TH1F("MC_NegPionsFromNeutralMeson_Pt","MC_NegPionsFromNeutralMeson_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
        fHistoMCNegPionsFromNeutralMesonPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoMCNegPionsFromNeutralMesonPt[iCut]->GetYaxis()->SetTitle("N_{#pi^{-}}");
        fHistoMCNegPionsFromNeutralMesonPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCNegPionsFromNeutralMesonPt[iCut]);
        fHistoMCPosPionsFromNeutralMesonPt[iCut]  = new TH1F("MC_PosPionsFromNeutralMeson_Pt","MC_PosPionsFromNeutralMeson_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
        fHistoMCPosPionsFromNeutralMesonPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoMCPosPionsFromNeutralMesonPt[iCut]->GetYaxis()->SetTitle("N_{#pi^{+}}");
        fHistoMCPosPionsFromNeutralMesonPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPosPionsFromNeutralMesonPt[iCut]);
        if (fDoMesonQA>0){
          fHistoMCAllMesonPt[iCut]                  = new TH1F("MC_AllNDM_Pt", "MC_AllNDM_Pt", HistoNPtBins, HistoPtRange[0], HistoPtRange[1]);
          fHistoMCAllMesonPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
          fHistoMCAllMesonPt[iCut]->GetYaxis()->SetTitle(Form("N_{%s}}", NameNDMLatex.Data()));
          fHistoMCAllMesonPt[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCAllMesonPt[iCut]);
          fHistoMCAllMesonEta[iCut]                 = new TH1F("MC_AllNDM_Eta", "MC_AllNDM_Eta", 200, -2, 2);
          fHistoMCAllMesonEta[iCut]->GetXaxis()->SetTitle("#eta");
          fHistoMCAllMesonEta[iCut]->GetYaxis()->SetTitle(Form("N_{%s}}", NameNDMLatex.Data()));
          fHistoMCAllMesonEta[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCAllMesonEta[iCut]);
          fHistoMCAllMesonPhi[iCut]                 = new TH1F("MC_AllNDM_Phi", "MC_AllNDM_Phi", 200, 0, TMath::TwoPi());
          fHistoMCAllMesonPhi[iCut]->GetXaxis()->SetTitle("#varphi");
          fHistoMCAllMesonPhi[iCut]->GetYaxis()->SetTitle(Form("N_{%s}}", NameNDMLatex.Data()));
          fHistoMCAllMesonPhi[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCAllMesonPhi[iCut]);
          fHistoMCMesonFromNeutralMesonPt[iCut]     = new TH1F("MC_NDMFromNeutralMeson_Pt", "MC_NDMFormNeutralMeson_Pt", HistoNPtBins, HistoPtRange[0], HistoPtRange[1]);
          fHistoMCMesonFromNeutralMesonPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
          fHistoMCMesonFromNeutralMesonPt[iCut]->GetYaxis()->SetTitle(Form("N_{%s}}", NameNDMLatex.Data()));
          fHistoMCMesonFromNeutralMesonPt[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCMesonFromNeutralMesonPt[iCut]);
          fHistoMCMesonFromNeutralMesonEta[iCut]    = new TH1F("MC_NDMFromNeutralMeson_Eta", "MC_NDMFromNeutralMeson_Eta", 200, -2, 2);
          fHistoMCMesonFromNeutralMesonEta[iCut]->GetXaxis()->SetTitle("#eta");
          fHistoMCMesonFromNeutralMesonEta[iCut]->GetYaxis()->SetTitle(Form("N_{%s}}", NameNDMLatex.Data()));
          fHistoMCMesonFromNeutralMesonEta[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCMesonFromNeutralMesonEta[iCut]);
          fHistoMCMesonFromNeutralMesonPhi[iCut]    = new TH1F("MC_NDMFromNeutralMeson_Phi", "MC_NDMFromNeutralMeson_Phi", 200, 0, TMath::TwoPi());
          fHistoMCMesonFromNeutralMesonPhi[iCut]->GetXaxis()->SetTitle("#varphi");
          fHistoMCMesonFromNeutralMesonPhi[iCut]->GetYaxis()->SetTitle(Form("N_{%s}}", NameNDMLatex.Data()));
          fHistoMCMesonFromNeutralMesonPhi[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCMesonFromNeutralMesonPhi[iCut]);
          fHistoMCAllPosPionsEta[iCut]               = new TH1F("MC_AllPosPions_Eta","MC_AllPosPions_Eta", 200, -2., 2.);
          fHistoMCAllPosPionsEta[iCut]->GetXaxis()->SetTitle("#eta");
          fHistoMCAllPosPionsEta[iCut]->GetYaxis()->SetTitle("N_{#pi^{+}}");
          fHistoMCAllPosPionsEta[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCAllPosPionsEta[iCut]);
          fHistoMCAllPosPionsPhi[iCut]               = new TH1F("MC_AllPosPions_Phi","MC_AllPosPions_Phi", 200, 0., TMath::TwoPi());
          fHistoMCAllPosPionsPhi[iCut]->GetXaxis()->SetTitle("#phi");
          fHistoMCAllPosPionsPhi[iCut]->GetYaxis()->SetTitle("N_{#pi^{+}}");
          fHistoMCAllPosPionsPhi[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCAllPosPionsPhi[iCut]);
          fHistoMCAllNegPionsEta[iCut]               = new TH1F("MC_AllNegPions_Eta","MC_AllNegPions_Eta", 200, -2., 2.);
          fHistoMCAllNegPionsEta[iCut]->GetXaxis()->SetTitle("#eta");
          fHistoMCAllNegPionsEta[iCut]->GetYaxis()->SetTitle("N_{#pi^{-}}");
          fHistoMCAllNegPionsEta[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCAllNegPionsEta[iCut]);
          fHistoMCAllNegPionsPhi[iCut]               = new TH1F("MC_AllNegPions_Phi","MC_AllNegPions_Phi", 200, 0., TMath::TwoPi());
          fHistoMCAllNegPionsPhi[iCut]->GetXaxis()->SetTitle("#phi");
          fHistoMCAllNegPionsPhi[iCut]->GetYaxis()->SetTitle("N_{#pi^{-}}");
          fHistoMCAllNegPionsPhi[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCAllNegPionsPhi[iCut]);
          fHistoMCPosPionsFromNeutralMesonEta[iCut]  = new TH1F("MC_PosPionsFromNeutralMeson_Eta","MC_PosPionsFromNeutralMeson_Eta", 200, -2., 2.);
          fHistoMCPosPionsFromNeutralMesonEta[iCut]->GetXaxis()->SetTitle("#eta");
          fHistoMCPosPionsFromNeutralMesonEta[iCut]->GetYaxis()->SetTitle("N_{#pi^{+}}");
          fHistoMCPosPionsFromNeutralMesonEta[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCPosPionsFromNeutralMesonEta[iCut]);
          fHistoMCPosPionsFromNeutralMesonPhi[iCut]  = new TH1F("MC_PosPionsFromNeutralMeson_Phi","MC_PosPionsFromNeutralMeson_Phi", 200, 0., TMath::TwoPi());
          fHistoMCPosPionsFromNeutralMesonPhi[iCut]->GetXaxis()->SetTitle("#phi");
          fHistoMCPosPionsFromNeutralMesonPhi[iCut]->GetYaxis()->SetTitle("N_{#pi^{+}}");
          fHistoMCPosPionsFromNeutralMesonPhi[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCPosPionsFromNeutralMesonPhi[iCut]);
          fHistoMCNegPionsFromNeutralMesonEta[iCut]  = new TH1F("MC_NegPionsFromNeutralMeson_Eta","MC_NegPionsFromNeutralMeson_Eta", 200, -2., 2.);
          fHistoMCNegPionsFromNeutralMesonEta[iCut]->GetXaxis()->SetTitle("#eta");
          fHistoMCNegPionsFromNeutralMesonEta[iCut]->GetYaxis()->SetTitle("N_{#pi^{-}}");
          fHistoMCNegPionsFromNeutralMesonEta[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCNegPionsFromNeutralMesonEta[iCut]);
          fHistoMCNegPionsFromNeutralMesonPhi[iCut]  = new TH1F("MC_NegPionsFromNeutralMeson_Phi","MC_NegPionsFromNeutralMeson_Phi", 200, 0., TMath::TwoPi());
          fHistoMCNegPionsFromNeutralMesonPhi[iCut]->GetXaxis()->SetTitle("#phi");
          fHistoMCNegPionsFromNeutralMesonPhi[iCut]->GetYaxis()->SetTitle("N_{#pi^{-}}");
          fHistoMCNegPionsFromNeutralMesonPhi[iCut]->Sumw2();

          fMCList[iCut]->Add(fHistoMCNegPionsFromNeutralMesonPhi[iCut]);
          fHistoMCHeavyAllPt[iCut]                            = new TH1F("MC_HeavyAll_Pt", "MC_HeavyAll_Pt", HistoNPtBins, HistoPtRange[0], HistoPtRange[1]);
          fHistoMCHeavyAllPt[iCut]->SetXTitle("p_{t} (GeV/c)");
          fHistoMCHeavyAllPt[iCut]->SetYTitle(Form("N_{%s}", NameNeutralMesonAnalyzedLatex.Data()));
          fHistoMCHeavyAllPt[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCHeavyAllPt[iCut]);
          fHistoMCHeavyAllEta[iCut]                           = new TH1F("MC_HeavyAll_Eta", "MC_HeavyAll_Eta", 200, -2., 2.);
          fHistoMCHeavyAllEta[iCut]->SetXTitle("#eta");
          fHistoMCHeavyAllEta[iCut]->SetYTitle(Form("N_{%s}", NameNeutralMesonAnalyzedLatex.Data()));
          fHistoMCHeavyAllEta[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCHeavyAllEta[iCut]);
          fHistoMCHeavyAllPhi[iCut]                           = new TH1F("MC_HeavyAll_Phi", "MC_HeavyAll_Phi", 200, 0., TMath::TwoPi());
          fHistoMCHeavyAllPhi[iCut]->SetXTitle("#phi{t} (GeV/c");
          fHistoMCHeavyAllPhi[iCut]->SetYTitle(Form("N_{%s}", NameNeutralMesonAnalyzedLatex.Data()));
          fHistoMCHeavyAllPhi[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCHeavyAllPhi[iCut]);
          fHistoMCHeavyChannelPt[iCut]                        = new TH1F("MC_HeavyChannel_Pt", "MC_HeavyChannel_Pt", HistoNPtBins, HistoPtRange[0], HistoPtRange[1]);
          fHistoMCHeavyChannelPt[iCut]->SetXTitle("p_{t} (GeV/c)");
          fHistoMCHeavyChannelPt[iCut]->SetYTitle(Form("N_{%s}", NameNeutralMesonAnalyzedLatex.Data()));
          fHistoMCHeavyChannelPt[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCHeavyChannelPt[iCut]);
          fHistoMCHeavyChannelEta[iCut]                       = new TH1F("MC_HeavyChannel_Eta", "MC_HeavyChannel_Eta", 200, -2., 2.);
          fHistoMCHeavyAllEta[iCut]->SetXTitle("#eta");
          fHistoMCHeavyAllEta[iCut]->SetYTitle(Form("N_{%s}", NameNeutralMesonAnalyzedLatex.Data()));
          fHistoMCHeavyAllEta[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCHeavyChannelEta[iCut]);
          fHistoMCHeavyChannelPhi[iCut]                       = new TH1F("MC_HeavyChannel_Phi", "MC_HeavyChannel_Phi", 200, 0., TMath::TwoPi());
          fHistoMCHeavyChannelPhi[iCut]->SetXTitle("#phi{t} (GeV/c");
          fHistoMCHeavyChannelPhi[iCut]->SetYTitle(Form("N_{%s}", NameNeutralMesonAnalyzedLatex.Data()));
          fHistoMCHeavyChannelPhi[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCHeavyChannelPhi[iCut]);
          fHistMCChannelNDMFromHeavyPt[iCut]                  = new TH1F("MC_NDMFromHeavyChannel_Pt", "MC_NDMFromHeavyChannel_Pt", HistoNPtBins, HistoPtRange[0], HistoPtRange[1]);
          fHistMCChannelNDMFromHeavyPt[iCut]->SetXTitle("p_{t} (GeV/c)");
          fHistMCChannelNDMFromHeavyPt[iCut]->SetYTitle(Form("N_{%s}", NameNDMLatex.Data()));
          fHistMCChannelNDMFromHeavyPt[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistMCChannelNDMFromHeavyPt[iCut]);
          fHistMCChannelNDMFromHeavyEta[iCut]                 = new TH1F("MC_NDMFromHeavyChannel_Eta", "MC_NDMFromHeavyChannel_Eta", 200, -2., 2.);
          fHistMCChannelNDMFromHeavyEta[iCut]->SetXTitle("#eta");
          fHistMCChannelNDMFromHeavyEta[iCut]->SetYTitle(Form("N_{%s}", NameNDMLatex.Data()));
          fHistMCChannelNDMFromHeavyEta[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistMCChannelNDMFromHeavyEta[iCut]);
          fHistMCChannelNDMFromHeavyPhi[iCut]                 = new TH1F("MC_NDMFromHeavyChannel_Phi", "MC_NDMFromHeavyChannel_Phi", 200, 0., TMath::TwoPi());
          fHistMCChannelNDMFromHeavyPhi[iCut]->SetXTitle("#phi");
          fHistMCChannelNDMFromHeavyPhi[iCut]->SetYTitle(Form("N_{%s}", NameNDMLatex.Data()));
          fHistMCChannelNDMFromHeavyPhi[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistMCChannelNDMFromHeavyPhi[iCut]);
          fHistMCChannelPiPlusFromHeavyPt[iCut]               = new TH1F("MC_PiPlusFromHeavyChannel_Pt", "MC_PiPlusFromHeavyChannel_Pt", HistoNPtBins, HistoPtRange[0], HistoPtRange[1]);
          fHistMCChannelPiPlusFromHeavyPt[iCut]->SetXTitle("p_{t} (GeV/c)");
          fHistMCChannelPiPlusFromHeavyPt[iCut]->SetYTitle("N_{#pi^{+}}");
          fHistMCChannelPiPlusFromHeavyPt[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistMCChannelPiPlusFromHeavyPt[iCut]);
          fHistMCChannelPiPlusFromHeavyEta[iCut]              = new TH1F("MC_PiPlusFromHeavyChannel_Eta", "MC_PiPlusFromHeavyChannel_Eta", 200, -2., 2.);
          fHistMCChannelPiPlusFromHeavyEta[iCut]->SetXTitle("#eta");
          fHistMCChannelPiPlusFromHeavyEta[iCut]->SetYTitle("N_{#pi^{+}}");
          fHistMCChannelPiPlusFromHeavyEta[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistMCChannelPiPlusFromHeavyEta[iCut]);
          fHistMCChannelPiPlusFromHeavyPhi[iCut]              = new TH1F("MC_PiPlusFromHeavyChannel_Phi", "MC_PiPlusFromHeavyChannel_Phi", 200, 0., TMath::TwoPi());
          fHistMCChannelPiPlusFromHeavyPhi[iCut]->SetXTitle("#phi");
          fHistMCChannelPiPlusFromHeavyPhi[iCut]->SetYTitle("N_{#pi^{+}}");
          fHistMCChannelPiPlusFromHeavyPhi[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistMCChannelPiPlusFromHeavyPhi[iCut]);
          fHistMCChannelPiMinusFromHeavyPt[iCut]              = new TH1F("MC_PiMinusFromHeavyChannel_Pt", "MC_PiMinusFromHeavyChannel_Pt", HistoNPtBins, HistoPtRange[0], HistoPtRange[1]);
          fHistMCChannelPiMinusFromHeavyPt[iCut]->SetXTitle("p_{t} (GeV/c)");
          fHistMCChannelPiMinusFromHeavyPt[iCut]->SetYTitle("N_{#pi^{-}}");
          fHistMCChannelPiMinusFromHeavyPt[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistMCChannelPiMinusFromHeavyPt[iCut]);
          fHistMCChannelPiMinusFromHeavyEta[iCut]             = new TH1F("MC_PiMinusFromHeavyChannel_Eta", "MC_PiMinusFromHeavyChannel_Eta", 200, -2., 2.);
          fHistMCChannelPiMinusFromHeavyEta[iCut]->SetXTitle("#eta");
          fHistMCChannelPiMinusFromHeavyEta[iCut]->SetYTitle("N_{#pi^{-}}");
          fHistMCChannelPiMinusFromHeavyEta[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistMCChannelPiMinusFromHeavyEta[iCut]);
          fHistMCChannelPiPMinusFromHeavyPhi[iCut]            = new TH1F("MC_PiMinusFromHeavyChannel_Phi", "MC_PiMinusFromHeavyChannel_Phi", 200, 0., TMath::TwoPi());
          fHistMCChannelPiPMinusFromHeavyPhi[iCut]->SetXTitle("#phi");
          fHistMCChannelPiPMinusFromHeavyPhi[iCut]->SetYTitle("N_{#pi^{-}}");
          fHistMCChannelPiPMinusFromHeavyPhi[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistMCChannelPiPMinusFromHeavyPhi[iCut]);

          fHistMCChannelNDMPtHeavyPt[iCut]                    = new TH2F("MC_CorrPtNDMHeavyChannel", "MC_CorrPtNDMHeavyChannel", HistoNPtBins, HistoPtRange[0], HistoPtRange[1], HistoNPtBins, HistoPtRange[0], HistoPtRange[1]);
          fHistMCChannelNDMPtHeavyPt[iCut]->SetXTitle(Form("p_{t, %s} (GeV/c)", NameNeutralMesonAnalyzedLatex.Data()));
          fHistMCChannelNDMPtHeavyPt[iCut]->SetYTitle(Form("p_{t, %s} (GeV/c)", NameNDMLatex.Data()));
          fHistMCChannelNDMPtHeavyPt[iCut]->SetZTitle(Form("N_{%s}", NameNeutralMesonAnalyzedLatex.Data()));
          fHistMCChannelNDMPtHeavyPt[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistMCChannelNDMPtHeavyPt[iCut]);
          fHistMCChannelPiPlusPtHeavyPt[iCut]                 = new TH2F("MC_CorrPtPiPlusHeavyChannel", "MC_CorrPtPiPlusHeavyChannel", HistoNPtBins, HistoPtRange[0], HistoPtRange[1], HistoNPtBins, HistoPtRange[0], HistoPtRange[1]);
          fHistMCChannelPiPlusPtHeavyPt[iCut]->SetXTitle(Form("p_{t, %s} (GeV/c)", NameNeutralMesonAnalyzedLatex.Data()));
          fHistMCChannelPiPlusPtHeavyPt[iCut]->SetYTitle("p_{t, #pi^{-}%s} (GeV/c)");
          fHistMCChannelPiPlusPtHeavyPt[iCut]->SetZTitle(Form("N_{%s}", NameNeutralMesonAnalyzedLatex.Data()));
          fHistMCChannelPiPlusPtHeavyPt[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistMCChannelPiPlusPtHeavyPt[iCut]);
          fHistMCChannelPiMinusPtHeavyPt[iCut]                = new TH2F("MC_CorrPtPiMinusHeavyChannel", "MC_CorrPtPiMinusHeavyChannel", HistoNPtBins, HistoPtRange[0], HistoPtRange[1], HistoNPtBins, HistoPtRange[0], HistoPtRange[1]);
          fHistMCChannelPiMinusPtHeavyPt[iCut]->SetXTitle(Form("p_{t, %s} (GeV/c)", NameNeutralMesonAnalyzedLatex.Data()));
          fHistMCChannelPiMinusPtHeavyPt[iCut]->SetYTitle("p_{t, #pi^{-}%s} (GeV/c)");
          fHistMCChannelPiMinusPtHeavyPt[iCut]->SetZTitle(Form("N_{%s}", NameNeutralMesonAnalyzedLatex.Data()));
          fHistMCChannelPiMinusPtHeavyPt[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistMCChannelPiMinusPtHeavyPt[iCut]);

          fHistoMCHeavyReconstructiblePt[iCut]                        = new TH1F("MC_HeavyReconstructible_Pt", "MC_HeavyReconstructible_Pt", HistoNPtBins, HistoPtRange[0], HistoPtRange[1]);
          fHistoMCHeavyReconstructiblePt[iCut]->SetXTitle("p_{t} (GeV/c)");
          fHistoMCHeavyReconstructiblePt[iCut]->SetYTitle(Form("N_{%s}", NameNeutralMesonAnalyzedLatex.Data()));
          fHistoMCHeavyReconstructiblePt[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCHeavyReconstructiblePt[iCut]);
          fHistoMCHeavyReconstructibleEta[iCut]                       = new TH1F("MC_HeavyReconstructible_Eta", "MC_HeavyReconstructible_Eta", 200, -2., 2.);
          fHistoMCHeavyAllEta[iCut]->SetXTitle("#eta");
          fHistoMCHeavyAllEta[iCut]->SetYTitle(Form("N_{%s}", NameNeutralMesonAnalyzedLatex.Data()));
          fHistoMCHeavyAllEta[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCHeavyReconstructibleEta[iCut]);
          fHistoMCHeavyReconstructiblePhi[iCut]                       = new TH1F("MC_HeavyReconstructible_Phi", "MC_HeavyReconstructible_Phi", 200, 0., TMath::TwoPi());
          fHistoMCHeavyReconstructiblePhi[iCut]->SetXTitle("#phi{t} (GeV/c");
          fHistoMCHeavyReconstructiblePhi[iCut]->SetYTitle(Form("N_{%s}", NameNeutralMesonAnalyzedLatex.Data()));
          fHistoMCHeavyReconstructiblePhi[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCHeavyReconstructiblePhi[iCut]);
          fHistMCReconstructibleNDMFromHeavyPt[iCut]                  = new TH1F("MC_NDMFromHeavyReconstructible_Pt", "MC_NDMFromHeavyReconstructible_Pt", HistoNPtBins, HistoPtRange[0], HistoPtRange[1]);
          fHistMCReconstructibleNDMFromHeavyPt[iCut]->SetXTitle("p_{t} (GeV/c)");
          fHistMCReconstructibleNDMFromHeavyPt[iCut]->SetYTitle(Form("N_{%s}", NameNDMLatex.Data()));
          fHistMCReconstructibleNDMFromHeavyPt[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistMCReconstructibleNDMFromHeavyPt[iCut]);
          fHistMCReconstructibleNDMFromHeavyEta[iCut]                 = new TH1F("MC_NDMFromHeavyReconstructible_Eta", "MC_NDMFromHeavyReconstructible_Eta", 200, -2., 2.);
          fHistMCReconstructibleNDMFromHeavyEta[iCut]->SetXTitle("#eta");
          fHistMCReconstructibleNDMFromHeavyEta[iCut]->SetYTitle(Form("N_{%s}", NameNDMLatex.Data()));
          fHistMCReconstructibleNDMFromHeavyEta[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistMCReconstructibleNDMFromHeavyEta[iCut]);
          fHistMCReconstructibleNDMFromHeavyPhi[iCut]                 = new TH1F("MC_NDMFromHeavyReconstructible_Phi", "MC_NDMFromHeavyReconstructible_Phi", 200, 0., TMath::TwoPi());
          fHistMCReconstructibleNDMFromHeavyPhi[iCut]->SetXTitle("#phi");
          fHistMCReconstructibleNDMFromHeavyPhi[iCut]->SetYTitle(Form("N_{%s}", NameNDMLatex.Data()));
          fHistMCReconstructibleNDMFromHeavyPhi[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistMCReconstructibleNDMFromHeavyPhi[iCut]);
          fHistMCReconstructiblePiPlusFromHeavyPt[iCut]               = new TH1F("MC_PiPlusFromHeavyReconstructible_Pt", "MC_PiPlusFromHeavyReconstructible_Pt", HistoNPtBins, HistoPtRange[0], HistoPtRange[1]);
          fHistMCReconstructiblePiPlusFromHeavyPt[iCut]->SetXTitle("p_{t} (GeV/c)");
          fHistMCReconstructiblePiPlusFromHeavyPt[iCut]->SetYTitle("N_{#pi^{+}}");
          fHistMCReconstructiblePiPlusFromHeavyPt[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistMCReconstructiblePiPlusFromHeavyPt[iCut]);
          fHistMCReconstructiblePiPlusFromHeavyEta[iCut]              = new TH1F("MC_PiPlusFromHeavyReconstructible_Eta", "MC_PiPlusFromHeavyReconstructible_Eta", 200, -2., 2.);
          fHistMCReconstructiblePiPlusFromHeavyEta[iCut]->SetXTitle("#eta");
          fHistMCReconstructiblePiPlusFromHeavyEta[iCut]->SetYTitle("N_{#pi^{+}}");
          fHistMCReconstructiblePiPlusFromHeavyEta[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistMCReconstructiblePiPlusFromHeavyEta[iCut]);
          fHistMCReconstructiblePiPlusFromHeavyPhi[iCut]              = new TH1F("MC_PiPlusFromHeavyReconstructible_Phi", "MC_PiPlusFromHeavyReconstructible_Phi", 200, 0., TMath::TwoPi());
          fHistMCReconstructiblePiPlusFromHeavyPhi[iCut]->SetXTitle("#phi");
          fHistMCReconstructiblePiPlusFromHeavyPhi[iCut]->SetYTitle("N_{#pi^{+}}");
          fHistMCReconstructiblePiPlusFromHeavyPhi[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistMCReconstructiblePiPlusFromHeavyPhi[iCut]);
          fHistMCReconstructiblePiMinusFromHeavyPt[iCut]              = new TH1F("MC_PiMinusFromHeavyReconstructible_Pt", "MC_PiMinusFromHeavyReconstructible_Pt", HistoNPtBins, HistoPtRange[0], HistoPtRange[1]);
          fHistMCReconstructiblePiMinusFromHeavyPt[iCut]->SetXTitle("p_{t} (GeV/c)");
          fHistMCReconstructiblePiMinusFromHeavyPt[iCut]->SetYTitle("N_{#pi^{-}}");
          fHistMCReconstructiblePiMinusFromHeavyPt[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistMCReconstructiblePiMinusFromHeavyPt[iCut]);
          fHistMCReconstructiblePiMinusFromHeavyEta[iCut]             = new TH1F("MC_PiMinusFromHeavyReconstructible_Eta", "MC_PiMinusFromHeavyReconstructible_Eta", 200, -2., 2.);
          fHistMCReconstructiblePiMinusFromHeavyEta[iCut]->SetXTitle("#eta");
          fHistMCReconstructiblePiMinusFromHeavyEta[iCut]->SetYTitle("N_{#pi^{-}}");
          fHistMCReconstructiblePiMinusFromHeavyEta[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistMCReconstructiblePiMinusFromHeavyEta[iCut]);
          fHistMCReconstructiblePiPMinusFromHeavyPhi[iCut]            = new TH1F("MC_PiMinusFromHeavyReconstructible_Phi", "MC_PiMinusFromHeavyReconstructible_Phi", 200, 0., TMath::TwoPi());
          fHistMCReconstructiblePiPMinusFromHeavyPhi[iCut]->SetXTitle("#phi");
          fHistMCReconstructiblePiPMinusFromHeavyPhi[iCut]->SetYTitle("N_{#pi^{-}}");
          fHistMCReconstructiblePiPMinusFromHeavyPhi[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistMCReconstructiblePiPMinusFromHeavyPhi[iCut]);

          fHistMCReconstructibleNDMPtHeavyPt[iCut]                    = new TH2F("MC_CorrPtNDMHeavyReconstructible", "MC_CorrPtNDMHeavyReconstructible", HistoNPtBins, HistoPtRange[0], HistoPtRange[1], HistoNPtBins, HistoPtRange[0], HistoPtRange[1]);
          fHistMCReconstructibleNDMPtHeavyPt[iCut]->SetXTitle(Form("p_{t, %s} (GeV/c)", NameNeutralMesonAnalyzedLatex.Data()));
          fHistMCReconstructibleNDMPtHeavyPt[iCut]->SetYTitle(Form("p_{t, %s} (GeV/c)", NameNDMLatex.Data()));
          fHistMCReconstructibleNDMPtHeavyPt[iCut]->SetZTitle(Form("N_{%s}", NameNeutralMesonAnalyzedLatex.Data()));
          fHistMCReconstructibleNDMPtHeavyPt[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistMCReconstructibleNDMPtHeavyPt[iCut]);
          fHistMCReconstructiblePiPlusPtHeavyPt[iCut]                 = new TH2F("MC_CorrPtPiPlusHeavyReconstructible", "MC_CorrPtPiPlusHeavyReconstructible", HistoNPtBins, HistoPtRange[0], HistoPtRange[1], HistoNPtBins, HistoPtRange[0], HistoPtRange[1]);
          fHistMCReconstructiblePiPlusPtHeavyPt[iCut]->SetXTitle(Form("p_{t, %s} (GeV/c)", NameNeutralMesonAnalyzedLatex.Data()));
          fHistMCReconstructiblePiPlusPtHeavyPt[iCut]->SetYTitle("p_{t, #pi^{-}%s} (GeV/c)");
          fHistMCReconstructiblePiPlusPtHeavyPt[iCut]->SetZTitle(Form("N_{%s}", NameNeutralMesonAnalyzedLatex.Data()));
          fHistMCReconstructiblePiPlusPtHeavyPt[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistMCReconstructiblePiPlusPtHeavyPt[iCut]);
          fHistMCReconstructiblePiMinusPtHeavyPt[iCut]                = new TH2F("MC_CorrPtPiMinusHeavyReconstructible", "MC_CorrPtPiMinusHeavyReconstructible", HistoNPtBins, HistoPtRange[0], HistoPtRange[1], HistoNPtBins, HistoPtRange[0], HistoPtRange[1]);
          fHistMCReconstructiblePiMinusPtHeavyPt[iCut]->SetXTitle(Form("p_{t, %s} (GeV/c)", NameNeutralMesonAnalyzedLatex.Data()));
          fHistMCReconstructiblePiMinusPtHeavyPt[iCut]->SetYTitle("p_{t, #pi^{-}%s} (GeV/c)");
          fHistMCReconstructiblePiMinusPtHeavyPt[iCut]->SetZTitle(Form("N_{%s}", NameNeutralMesonAnalyzedLatex.Data()));
          fHistMCReconstructiblePiMinusPtHeavyPt[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistMCReconstructiblePiMinusPtHeavyPt[iCut]);
        }
      }
      fHistoMCHNMPiPlPiMiNDMPt[iCut]         = new TH1F("MC_HNM_Pt","MC_HNM_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
      fHistoMCHNMPiPlPiMiNDMPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      fHistoMCHNMPiPlPiMiNDMPt[iCut]->GetYaxis()->SetTitle("N_{HNM}");
      fHistoMCHNMPiPlPiMiNDMPt[iCut]->Sumw2();
      fMCList[iCut]->Add(fHistoMCHNMPiPlPiMiNDMPt[iCut]);

      fHistoMCHNMPiPlPiMiNDMInAccPt[iCut]    = new TH1F("MC_HNMInAcc_Pt","MC_HNMInAcc_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
      fHistoMCHNMPiPlPiMiNDMInAccPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      fHistoMCHNMPiPlPiMiNDMInAccPt[iCut]->GetYaxis()->SetTitle("A #times N_{HNM}");
      fHistoMCHNMPiPlPiMiNDMInAccPt[iCut]->Sumw2();
      fMCList[iCut]->Add(fHistoMCHNMPiPlPiMiNDMInAccPt[iCut]);

      if(!fDoLightOutput){

        fHistoMCHNMInAccVsNDMPt[iCut]    = new TH2F("MC_HNMInAccVsNDMPt","MC_HNMInAccVsNDMPt",200,HistoPtRange[0],HistoPtRange[1],200,HistoPtRange[0],HistoPtRange[1]);
        fHistoMCHNMInAccVsNDMPt[iCut]->GetXaxis()->SetTitle("p_{T} of HNM (GeV/c)");
        fHistoMCHNMInAccVsNDMPt[iCut]->GetYaxis()->SetTitle("p_{T} of NDM (GeV/c)");
        fHistoMCHNMInAccVsNDMPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCHNMInAccVsNDMPt[iCut]);
        fHistoMCHNMPiPlPiMiNDMEta[iCut]         = new TH1F("MC_HNM_Eta","MC_HNM_Eta",200,-2,2);
        fHistoMCHNMPiPlPiMiNDMEta[iCut]->GetXaxis()->SetTitle("#eta");
        fHistoMCHNMPiPlPiMiNDMEta[iCut]->GetYaxis()->SetTitle("N_{HNM}");
        fHistoMCHNMPiPlPiMiNDMEta[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCHNMPiPlPiMiNDMEta[iCut]);

        fHistoMCHNMPiPlPiMiNDMPhi[iCut]         = new TH1F("MC_HNM_Phi","MC_HNM_Phi",200,0,2 * TMath::Pi());
        fHistoMCHNMPiPlPiMiNDMPhi[iCut]->GetXaxis()->SetTitle("#phi");
        fHistoMCHNMPiPlPiMiNDMPhi[iCut]->GetYaxis()->SetTitle("N_{HNM}");
        fHistoMCHNMPiPlPiMiNDMPhi[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCHNMPiPlPiMiNDMPhi[iCut]);

        if(fDoMesonQA > 1){
          fHistoMCHNMPiPlPiMiNDMEtavsPt[iCut]         = new TH2F("MC_HNM_EtavsPt","MC_HNM_EtavsPt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1],100,-1,1);
          fHistoMCHNMPiPlPiMiNDMEtavsPt[iCut]->GetXaxis()->SetTitle("#it{p}_{T}");
          fHistoMCHNMPiPlPiMiNDMEtavsPt[iCut]->GetYaxis()->SetTitle("#eta");
          fHistoMCHNMPiPlPiMiNDMEtavsPt[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCHNMPiPlPiMiNDMEtavsPt[iCut]);
        }
      }

      fTrueList[iCut]                           = new TList();
      fTrueList[iCut]->SetName(nameTrueRecList.Data());
      fTrueList[iCut]->SetOwner(kTRUE);
      fCutFolder[iCut]->Add(fTrueList[iCut]);

      if(!fDoLightOutput){
        if (fNDMRecoMode < 2){
          fHistoTrueConvGammaPt[iCut]                 = new TH1F("ESD_TrueConvGamma_Pt","ESD_TrueConvGamma_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTrueConvGammaPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
          fHistoTrueConvGammaPt[iCut]->GetYaxis()->SetTitle("N_{#gamma,conv}");
          fHistoTrueConvGammaPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTrueConvGammaPt[iCut]);
          fHistoDoubleCountTrueConvGammaRPt[iCut]     = new TH2F("ESD_TrueDoubleCountConvGamma_R_Pt","ESD_TrueDoubleCountConvGamma_R_Pt",800,0,200,HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoDoubleCountTrueConvGammaRPt[iCut]->GetXaxis()->SetTitle("R_{conv} (cm)");
          fHistoDoubleCountTrueConvGammaRPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          fHistoDoubleCountTrueConvGammaRPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoDoubleCountTrueConvGammaRPt[iCut]);
          fHistoTrueConvGammaFromNeutralMesonPt[iCut] = new TH1F("ESD_TrueConvGammaFromNeutralMeson_Pt","ESD_TrueConvGammaFromNeutralMeson_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTrueConvGammaFromNeutralMesonPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
          fHistoTrueConvGammaFromNeutralMesonPt[iCut]->GetYaxis()->SetTitle("N_{#gamma,conv}");
          fHistoTrueConvGammaFromNeutralMesonPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTrueConvGammaFromNeutralMesonPt[iCut]);
        }
        if (fNDMRecoMode > 0){
          fHistoTrueClusterGammaPt[iCut]                  = new TH1F("ESD_TrueClusterGamma_Pt","ESD_TrueClusterGamma_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTrueClusterGammaPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
          fHistoTrueClusterGammaPt[iCut]->GetYaxis()->SetTitle("N_{#gamma,cluster}");
          fHistoTrueClusterGammaPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTrueClusterGammaPt[iCut]);
          fHistoTrueClusterGammaFromNeutralMesonPt[iCut]  = new TH1F("ESD_TrueClusterGammaFromNeutralMeson_Pt","ESD_TrueClusterGammaFromNeutralMeson_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTrueClusterGammaFromNeutralMesonPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
          fHistoTrueClusterGammaFromNeutralMesonPt[iCut]->GetYaxis()->SetTitle("N_{#gamma,cluster}");
          fHistoTrueClusterGammaFromNeutralMesonPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTrueClusterGammaFromNeutralMesonPt[iCut]);
        }
        fHistoTruePosPionPt[iCut]                       = new TH1F("ESD_TruePosPion_Pt","ESD_TruePosPion_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
        fHistoTruePosPionPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoTruePosPionPt[iCut]->GetYaxis()->SetTitle("N_{#pi^{+}}");
        fHistoTruePosPionPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTruePosPionPt[iCut]);
        fHistoTrueNegPionPt[iCut]                       = new TH1F("ESD_TrueNegPion_Pt","ESD_TrueNegPion_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
        fHistoTrueNegPionPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoTrueNegPionPt[iCut]->GetYaxis()->SetTitle("N_{#pi^{-}}");
        fHistoTrueNegPionPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTrueNegPionPt[iCut]);

        fHistoTrueNegPionFromNeutralMesonPt[iCut]       = new TH1F("ESD_TrueNegPionFromNeutralMeson_Pt","ESD_TrueNegPionFromNeutralMeson_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
        fHistoTrueNegPionFromNeutralMesonPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoTrueNegPionFromNeutralMesonPt[iCut]->GetYaxis()->SetTitle("N_{#pi^{-}}");
        fHistoTrueNegPionFromNeutralMesonPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTrueNegPionFromNeutralMesonPt[iCut]);
        fHistoTruePosPionFromNeutralMesonPt[iCut]       = new TH1F("ESD_TruePosPionFromNeutralMeson_Pt","ESD_TruePosPionFromNeutralMeson_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
        fHistoTruePosPionFromNeutralMesonPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoTruePosPionFromNeutralMesonPt[iCut]->GetYaxis()->SetTitle("N_{#pi^{+}}");
        fHistoTruePosPionFromNeutralMesonPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTruePosPionFromNeutralMesonPt[iCut]);

        fHistoDoubleCountTruePi0InvMassPt[iCut]         = new TH2F("ESD_TrueDoubleCountPi0_InvMass_Pt","ESD_TrueDoubleCountPi0_InvMass_Pt",800,0,0.8,HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
        fHistoDoubleCountTruePi0InvMassPt[iCut]->GetXaxis()->SetTitle("M_{#gamma #gamma} (GeV/c^{2})");
        fHistoDoubleCountTruePi0InvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoDoubleCountTruePi0InvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoDoubleCountTruePi0InvMassPt[iCut]);
        fHistoDoubleCountTrueHNMInvMassPt[iCut]         = new TH2F("ESD_TrueDoubleCountHNM_InvMass_Pt","ESD_TrueDoubleCountHNM_InvMass_Pt",800,0,0.8,HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
        fHistoDoubleCountTrueHNMInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#eta} (GeV/c^{2})");
        fHistoDoubleCountTrueHNMInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoDoubleCountTrueHNMInvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoDoubleCountTrueHNMInvMassPt[iCut]);
      }
      fHistoTrueMotherPiPlPiMiNDMInvMassPt[iCut]       = new TH2F("ESD_TrueMotherPiPlPiMiNDM_InvMass_Pt","ESD_TrueMotherPiPlPiMiNDM_InvMass_Pt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
      fHistoTrueMotherPiPlPiMiNDMInvMassPt[iCut]->Sumw2();
      fHistoTrueMotherPiPlPiMiNDMInvMassPt[iCut]->GetXaxis()->SetTitle(Form("M_{#pi^{+} #pi^{-} %s} (GeV/c^{2})",NameNDMLatex.Data()));
      fHistoTrueMotherPiPlPiMiNDMInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
      fTrueList[iCut]->Add(fHistoTrueMotherPiPlPiMiNDMInvMassPt[iCut]);

      fHistoTrueMotherPiPlPiMiNDMInvMassPtSubNDM[iCut]       = new TH2F("ESD_TrueMotherPiPlPiMiNDM_InvMass_SubNDM_Pt","ESD_TrueMotherPiPlPiMiNDM_InvMass_SubNDM_Pt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
      fHistoTrueMotherPiPlPiMiNDMInvMassPtSubNDM[iCut]->Sumw2();
      fHistoTrueMotherPiPlPiMiNDMInvMassPtSubNDM[iCut]->GetXaxis()->SetTitle(Form("M_{#pi^{+} #pi^{-} %s} (GeV/c^{2})",NameNDMLatex.Data()));
      fHistoTrueMotherPiPlPiMiNDMInvMassPtSubNDM[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
      fTrueList[iCut]->Add(fHistoTrueMotherPiPlPiMiNDMInvMassPtSubNDM[iCut]);

      fHistoTrueMotherPiPlPiMiNDMInvMassPtFixedPzNDM[iCut]       = new TH2F("ESD_TrueMotherPiPlPiMiNDM_InvMass_FixedPzNDM_Pt","ESD_TrueMotherPiPlPiMiNDM_InvMass_FixedPzNDM_Pt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
      fHistoTrueMotherPiPlPiMiNDMInvMassPtFixedPzNDM[iCut]->Sumw2();
      fHistoTrueMotherPiPlPiMiNDMInvMassPtFixedPzNDM[iCut]->GetXaxis()->SetTitle(Form("M_{#pi^{+} #pi^{-} %s} (GeV/c^{2})",NameNDMLatex.Data()));
      fHistoTrueMotherPiPlPiMiNDMInvMassPtFixedPzNDM[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
      fTrueList[iCut]->Add(fHistoTrueMotherPiPlPiMiNDMInvMassPtFixedPzNDM[iCut]);

      fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromDifferent[iCut]       = new TH2F("ESD_TrueMotherPiPlPiMiNDMInvMassPt_FromDifferent","ESD_TrueMotherPiPlPiMiNDMInvMassPt_FromDifferent",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
      fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromDifferent[iCut]->Sumw2();
      fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromDifferent[iCut]->GetXaxis()->SetTitle(Form("M_{#pi^{+} #pi^{-} %s} (GeV/c^{2})",NameNDMLatex.Data()));
      fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromDifferent[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
      fTrueList[iCut]->Add(fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromDifferent[iCut]);

      fHistoTruePiPlPiMiNDMCombinatoricalInvMassPt[iCut]  = new TH2F("ESD_TruePiPlPiMiNDMCombinatorical_InvMassPt","ESD_TruePiPlPiMiNDMCombinatorical_InvMassPt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
      fHistoTruePiPlPiMiNDMCombinatoricalInvMassPt[iCut]->GetXaxis()->SetTitle(Form("M_{#pi^{+}#pi^{-}%s} (GeV/c^{2})",NameNDMLatex.Data()));
      fHistoTruePiPlPiMiNDMCombinatoricalInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
      if (fIsMC>1) fHistoTruePiPlPiMiNDMCombinatoricalInvMassPt[iCut]->Sumw2();
      fTrueList[iCut]->Add(fHistoTruePiPlPiMiNDMCombinatoricalInvMassPt[iCut]);

      fHistoTruePiPlPiMiNDMContaminationInvMassPt[iCut]  = new TH2F("ESD_TruePiPlPiMiNDMContamination_InvMassPt","ESD_TruePiPlPiMiNDMContamination_InvMassPt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
      fHistoTruePiPlPiMiNDMContaminationInvMassPt[iCut]->GetXaxis()->SetTitle(Form("M_{#pi^{+}#pi^{-}%s} (GeV/c^{2})",NameNDMLatex.Data()));
      fHistoTruePiPlPiMiNDMContaminationInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
      if (fIsMC>1) fHistoTruePiPlPiMiNDMContaminationInvMassPt[iCut]->Sumw2();
      fTrueList[iCut]->Add(fHistoTruePiPlPiMiNDMContaminationInvMassPt[iCut]);

      if(!fDoLightOutput){
        fHistoTrueMotherGammaGammaInvMassPt[iCut]           = new TH2F("ESD_TrueMotherGG_InvMass_Pt","ESD_TrueMotherGG_InvMass_Pt",HistoNMassBinsDecayMeson,HistoMassRangeNDM[0],HistoMassRangeNDM[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
        fHistoTrueMotherGammaGammaInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#gamma #gamma} (GeV/c^{2})");
        fHistoTrueMotherGammaGammaInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoTrueMotherGammaGammaInvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTrueMotherGammaGammaInvMassPt[iCut]);
        fHistoTrueMotherGammaGammaFromHNMInvMassPt[iCut]    = new TH2F("ESD_TrueMotherGGFromHNM_InvMass_Pt","ESD_TrueMotherGGFromHNM_InvMass_Pt",HistoNMassBinsDecayMeson,HistoMassRangeNDM[0],HistoMassRangeNDM[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
        fHistoTrueMotherGammaGammaFromHNMInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#gamma #gamma} (GeV/c^{2})");
        fHistoTrueMotherGammaGammaFromHNMInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoTrueMotherGammaGammaFromHNMInvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTrueMotherGammaGammaFromHNMInvMassPt[iCut]);
        fHistoTrueAngleSum[iCut]                            = new TH2F("ESD_TrueMother_AngleSum_Pt","ESD_TrueMother_AngleSum_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1],720,0,2*TMath::Pi());
        fHistoTrueAngleSum[iCut]->GetXaxis()->SetTitle("#sum #angle");
        fHistoTrueAngleSum[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoTrueAngleSum[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTrueAngleSum[iCut]);

        fHistoTrueHNMesonPtvsNDMPt[iCut]                            = new TH2F("ESD_TrueMother_HNMesonPtvsNDMPt","ESD_TrueMother_HNMesonPtvsNDMPt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
      fHistoTrueHNMesonPtvsNDMPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c) of HNM");
        fHistoTrueHNMesonPtvsNDMPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c) of NDM");
        fHistoTrueHNMesonPtvsNDMPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTrueHNMesonPtvsNDMPt[iCut]);

        if (fDoMesonQA>0){
          fHistoTruePionPionInvMassPt[iCut]                 = new TH2F("ESD_TruePiPlusPiNeg_InvMassPt","ESD_TruePiPlusPiNeg_InvMassPt",HistoNMassBinsPiPlusPiMinus,HistoMassRangePiPlusPiMinus[0],HistoMassRangePiPlusPiMinus[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTruePionPionInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} (GeV/c^{2})");
          fHistoTruePionPionInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          fHistoTruePionPionInvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePionPionInvMassPt[iCut]);
          fHistoTruePionPionFromSameMotherInvMassPt[iCut]   = new TH2F("ESD_TruePiPlusPiNegFromSameMother_InvMassPt","ESD_TruePiPlusPiNegFromSameMother_InvMassPt",HistoNMassBinsPiPlusPiMinus,HistoMassRangePiPlusPiMinus[0],HistoMassRangePiPlusPiMinus[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTruePionPionFromSameMotherInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} (GeV/c^{2})");
          fHistoTruePionPionFromSameMotherInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          fHistoTruePionPionFromSameMotherInvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePionPionFromSameMotherInvMassPt[iCut]);
          fHistoTruePionPionFromHNMInvMassPt[iCut]          = new TH2F("ESD_TruePiPlusPiNegFromHNM_InvMassPt","ESD_TruePiPlusPiNegFromHNM_InvMassPt",HistoNMassBinsPiPlusPiMinus,HistoMassRangePiPlusPiMinus[0],HistoMassRangePiPlusPiMinus[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTruePionPionFromHNMInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} (GeV/c^{2})");
          fHistoTruePionPionFromHNMInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          fHistoTruePionPionFromHNMInvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePionPionFromHNMInvMassPt[iCut]);

          fHistoTruevParticleChi2PerNDF[iCut] = new TH1F("fHistoTruevParticleChi2PerNDF","fHistoTruevParticleChi2PerNDF",300,0,300);
          fHistoTruevParticleChi2PerNDF[iCut]->GetXaxis()->SetTitle("#chi^{2}/ndf");
          fHistoTruevParticleChi2PerNDF[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruevParticleChi2PerNDF[iCut]);

          fHistoTruevParticleFromSameMotherChi2PerNDF[iCut] = new TH1F("fHistoTruevParticleFromSameMotherChi2PerNDF","fHistoTruevParticleFromSameMotherChi2PerNDF",300,0,300);
          fHistoTruevParticleFromSameMotherChi2PerNDF[iCut]->GetXaxis()->SetTitle("#chi^{2}/ndf");
          fHistoTruevParticleFromSameMotherChi2PerNDF[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruevParticleFromSameMotherChi2PerNDF[iCut]);

          fHistoTruevParticleFromHNMChi2PerNDF[iCut] = new TH1F("fHistoTruevParticleFromHNMChi2PerNDF","fHistoTruevParticleFromHNMChi2PerNDF",300,0,300);
          fHistoTruevParticleFromHNMChi2PerNDF[iCut]->GetXaxis()->SetTitle("#chi^{2}/ndf");
          fHistoTruevParticleFromHNMChi2PerNDF[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruevParticleFromHNMChi2PerNDF[iCut]);

          fHistoTruevParticledS[iCut] = new TH1F("fHistoTruevParticledS","fHistoTruevParticledS",500,-4,4);
          fHistoTruevParticledS[iCut]->GetXaxis()->SetTitle("dS");
          fHistoTruevParticledS[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruevParticledS[iCut]);

          fHistoTruevParticleFromSameMotherdS[iCut] = new TH1F("fHistoTruevParticleFromSameMotherdS","fHistoTruevParticleFromSameMotherdS",400,-4,4);
          fHistoTruevParticleFromSameMotherdS[iCut]->GetXaxis()->SetTitle("dS");
          fHistoTruevParticleFromSameMotherdS[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruevParticleFromSameMotherdS[iCut]);

          fHistoTruevParticleFromHNMdS[iCut] = new TH1F("fHistoTruevParticleFromHNMdS","fHistoTruevParticleFromHNMdS",400,-4,4);
          fHistoTruevParticleFromHNMdS[iCut]->GetXaxis()->SetTitle("dS");
          fHistoTruevParticleFromHNMdS[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruevParticleFromHNMdS[iCut]);

          fHistoTrueMesonFlags[iCut]           = new TH1F("TrueMesonFlags","TrueMesonFlags",11,0.5,11.5);
          fHistoTrueMesonFlags[iCut]->GetXaxis()->SetBinLabel(1,  "All candidates");
          fHistoTrueMesonFlags[iCut]->GetXaxis()->SetBinLabel(2,  "Same mother");
          fHistoTrueMesonFlags[iCut]->GetXaxis()->SetBinLabel(3,  "True");
          fHistoTrueMesonFlags[iCut]->GetXaxis()->SetBinLabel(4,  "Not same mother");
          fHistoTrueMesonFlags[iCut]->GetXaxis()->SetBinLabel(5,  "Wrongly identified pions");
          fHistoTrueMesonFlags[iCut]->GetXaxis()->SetBinLabel(6,  "Wrongly identified pi0");
          fHistoTrueMesonFlags[iCut]->GetXaxis()->SetBinLabel(7,  "Wrongly identified pi+");
          fHistoTrueMesonFlags[iCut]->GetXaxis()->SetBinLabel(8,  "Wrongly identified pi-");
          fHistoTrueMesonFlags[iCut]->GetXaxis()->SetBinLabel(9,  "Wrongly identified multiple");
          fHistoTrueMesonFlags[iCut]->GetXaxis()->SetBinLabel(10, "Problem with pi0 flag");
          fHistoTrueMesonFlags[iCut]->GetXaxis()->SetBinLabel(11, "Problem with meson declaration flag");
          fTrueList[iCut]->Add(fHistoTrueMesonFlags[iCut]);

          TString TStr_MesonType_TrueMotherPiPlPiMiNDMInvMassPt_FromEtaOmega="Eta";
          if (fSelectedHeavyNeutralMeson==0){
              TStr_MesonType_TrueMotherPiPlPiMiNDMInvMassPt_FromEtaOmega="Omega";
          }
          fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromEtaOmega[iCut]    = new TH2F("ESD_TrueMotherPiPlPiMiNDMInvMassPt_FromEtaOmega", Form("ESD_TrueMotherPiPlPiMiNDMInvMassPt_From%s", TStr_MesonType_TrueMotherPiPlPiMiNDMInvMassPt_FromEtaOmega.Data()), HistoNMassBins, HistoMassRange[0], HistoMassRange[1], HistoNPtBins, HistoPtRange[0], HistoPtRange[1]);
          fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromEtaOmega[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} (GeV/c^{2})");
          fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromEtaOmega[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          if (fIsMC>1) fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromEtaOmega[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromEtaOmega[iCut]);

          fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromRho[iCut]    = new TH2F("ESD_TrueMotherPiPlPiMiNDMInvMassPt_FromRho","ESD_TrueMotherPiPlPiMiNDMInvMassPt_FromRho", HistoNMassBins, HistoMassRange[0], HistoMassRange[1], HistoNPtBins, HistoPtRange[0], HistoPtRange[1]);
          fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromRho[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} (GeV/c^{2})");
          fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromRho[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          if (fIsMC>1) fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromRho[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromRho[iCut]);

          fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromK0s[iCut]    = new TH2F("ESD_TrueMotherPiPlPiMiNDMInvMassPt_FromK0s","ESD_TrueMotherPiPlPiMiNDMInvMassPt_FromK0s", HistoNMassBins, HistoMassRange[0], HistoMassRange[1], HistoNPtBins, HistoPtRange[0], HistoPtRange[1]);
          fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromK0s[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} (GeV/c^{2})");
          fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromK0s[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          if (fIsMC>1) fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromK0s[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromK0s[iCut]);

          fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromK0l[iCut]    = new TH2F("ESD_TrueMotherPiPlPiMiNDMInvMassPt_FromK0l","ESD_TrueMotherPiPlPiMiNDMInvMassPt_FromK0l", HistoNMassBins, HistoMassRange[0], HistoMassRange[1], HistoNPtBins, HistoPtRange[0], HistoPtRange[1]);
          fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromK0l[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} (GeV/c^{2})");
          fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromK0l[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          if (fIsMC>1) fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromK0l[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromK0l[iCut]);

          fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromEtaPrime[iCut]    = new TH2F("ESD_TrueMotherPiPlPiMiNDMInvMassPt_FromEtaPrime","ESD_TrueMotherPiPlPiMiNDMInvMassPt_FromEtaPrime", HistoNMassBins, HistoMassRange[0], HistoMassRange[1], HistoNPtBins, HistoPtRange[0], HistoPtRange[1]);
          fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromEtaPrime[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} (GeV/c^{2})");
          fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromEtaPrime[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          if (fIsMC>1) fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromEtaPrime[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromEtaPrime[iCut]);

          fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromOther[iCut]    = new TH2F("ESD_TrueMotherPiPlPiMiNDMInvMassPt_FromOther","ESD_TrueMotherPiPlPiMiNDMInvMassPt_FromOther", HistoNMassBins, HistoMassRange[0], HistoMassRange[1], HistoNPtBins, HistoPtRange[0], HistoPtRange[1]);
          fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromOther[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} (GeV/c^{2})");
          fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromOther[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          if (fIsMC>1) fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromOther[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromOther[iCut]);

          fHistoTruePiPlPiMiSameMotherFromEtaInvMassPt[iCut]    = new TH2F("ESD_TruePiPlPiMiSameMotherFromEta_InvMassPt","ESD_TruePiPlPiMiSameMotherFromEta_InvMassPt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTruePiPlPiMiSameMotherFromEtaInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} (GeV/c^{2})");
          fHistoTruePiPlPiMiSameMotherFromEtaInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          if (fIsMC>1) fHistoTruePiPlPiMiSameMotherFromEtaInvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePiPlPiMiSameMotherFromEtaInvMassPt[iCut]);
          fHistoTruePiPlPiMiSameMotherFromOmegaInvMassPt[iCut]  = new TH2F("ESD_TruePiPlPiMiSameMotherFromOmega_InvMassPt","ESD_TruePiPlPiMiSameMotherFromOmega_InvMassPt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTruePiPlPiMiSameMotherFromOmegaInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} (GeV/c^{2})");
          fHistoTruePiPlPiMiSameMotherFromOmegaInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          if (fIsMC>1) fHistoTruePiPlPiMiSameMotherFromOmegaInvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePiPlPiMiSameMotherFromOmegaInvMassPt[iCut]);

          fHistoTruePiPlPiMiSameMotherFromRhoInvMassPt[iCut]    = new TH2F("ESD_TruePiPlPiMiSameMotherFromRho_InvMassPt","ESD_TruePiPlPiMiSameMotherFromRho_InvMassPt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTruePiPlPiMiSameMotherFromRhoInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} (GeV/c^{2})");
          fHistoTruePiPlPiMiSameMotherFromRhoInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          if (fIsMC>1) fHistoTruePiPlPiMiSameMotherFromRhoInvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePiPlPiMiSameMotherFromRhoInvMassPt[iCut]);

          fHistoTruePiPlPiMiSameMotherFromEtaPrimeInvMassPt[iCut] = new TH2F("ESD_TruePiPlPiMiSameMotherFromEtaPrime_InvMassPt","ESD_TruePiPlPiMiSameMotherFromEtaPrime_InvMassPt",
                                                                                 HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTruePiPlPiMiSameMotherFromEtaPrimeInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} (GeV/c^{2})");
          fHistoTruePiPlPiMiSameMotherFromEtaPrimeInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          if (fIsMC>1) fHistoTruePiPlPiMiSameMotherFromEtaPrimeInvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePiPlPiMiSameMotherFromEtaPrimeInvMassPt[iCut]);

          fHistoTruePiPlPiMiSameMotherFromK0sInvMassPt[iCut]    = new TH2F("ESD_TruePiPlPiMiSameMotherFromK0s_InvMassPt","ESD_TruePiPlPiMiSameMotherFromK0s_InvMassPt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTruePiPlPiMiSameMotherFromK0sInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} (GeV/c^{2})");
          fHistoTruePiPlPiMiSameMotherFromK0sInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          if (fIsMC>1) fHistoTruePiPlPiMiSameMotherFromK0sInvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePiPlPiMiSameMotherFromK0sInvMassPt[iCut]);

          fHistoTruePiPlPiMiSameMotherFromK0lInvMassPt[iCut]    = new TH2F("ESD_TruePiPlPiMiSameMotherFromK0l_InvMassPt","ESD_TruePiPlPiMiSameMotherFromK0l_InvMassPt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTruePiPlPiMiSameMotherFromK0lInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} (GeV/c^{2})");
          fHistoTruePiPlPiMiSameMotherFromK0lInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          if (fIsMC>1) fHistoTruePiPlPiMiSameMotherFromK0lInvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePiPlPiMiSameMotherFromK0lInvMassPt[iCut]);

          fHistoTruePiPlPiMiSameMotherFromOtherlInvMassPt[iCut]    = new TH2F("ESD_TruePiPlPiMiSameMotherFromOther_InvMassPt","ESD_TruePiPlPiMiSameMotherFromOther_InvMassPt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTruePiPlPiMiSameMotherFromOtherlInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} (GeV/c^{2})");
          fHistoTruePiPlPiMiSameMotherFromOtherlInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          if (fIsMC>1) fHistoTruePiPlPiMiSameMotherFromOtherlInvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePiPlPiMiSameMotherFromOtherlInvMassPt[iCut]);

          fHistoTruePiMiPiZeroSameMotherFromEtaInvMassPt[iCut]  = new TH2F("ESD_TruePiMiPiZeroSameMotherFromEta_InvMassPt","ESD_TruePiMiPiZeroSameMotherFromEta_InvMassPt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTruePiMiPiZeroSameMotherFromEtaInvMassPt[iCut]->GetXaxis()->SetTitle(Form("M_{#pi^{-}%s} (GeV/c^{2})",NameNDMLatex.Data()));
          fHistoTruePiMiPiZeroSameMotherFromEtaInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          if (fIsMC>1) fHistoTruePiMiPiZeroSameMotherFromEtaInvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePiMiPiZeroSameMotherFromEtaInvMassPt[iCut]);
          fHistoTruePiMiPiZeroSameMotherFromOmegaInvMassPt[iCut] = new TH2F("ESD_TruePiMiPiZeroSameMotherFromOmega_InvMassPt","ESD_TruePiMiPiZeroSameMotherFromOmega_InvMassPt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTruePiMiPiZeroSameMotherFromOmegaInvMassPt[iCut]->GetXaxis()->SetTitle(Form("M_{#pi^{-}%s} (GeV/c^{2})",NameNDMLatex.Data()));
          fHistoTruePiMiPiZeroSameMotherFromOmegaInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          if (fIsMC>1) fHistoTruePiMiPiZeroSameMotherFromOmegaInvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePiMiPiZeroSameMotherFromOmegaInvMassPt[iCut]);
          fHistoTruePiMiPiZeroSameMotherFromRhoInvMassPt[iCut]  = new TH2F("ESD_TruePiMiPiZeroSameMotherFromRho_InvMassPt","ESD_TruePiMiPiZeroSameMotherFromRho_InvMassPt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTruePiMiPiZeroSameMotherFromRhoInvMassPt[iCut]->GetXaxis()->SetTitle(Form("M_{#pi^{-}%s} (GeV/c^{2})",NameNDMLatex.Data()));
          fHistoTruePiMiPiZeroSameMotherFromRhoInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          if (fIsMC>1) fHistoTruePiMiPiZeroSameMotherFromRhoInvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePiMiPiZeroSameMotherFromRhoInvMassPt[iCut]);

          fHistoTruePiMiPiZeroSameMotherFromK0lInvMassPt[iCut]  = new TH2F("ESD_TruePiMiPiZeroSameMotherFromK0l_InvMassPt","ESD_TruePiMiPiZeroSameMotherFromK0l_InvMassPt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTruePiMiPiZeroSameMotherFromK0lInvMassPt[iCut]->GetXaxis()->SetTitle(Form("M_{#pi^{-}%s} (GeV/c^{2})",NameNDMLatex.Data()));
          fHistoTruePiMiPiZeroSameMotherFromK0lInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          if (fIsMC>1) fHistoTruePiMiPiZeroSameMotherFromK0lInvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePiMiPiZeroSameMotherFromK0lInvMassPt[iCut]);

          fHistoTruePiMiPiZeroSameMotherFromOtherlInvMassPt[iCut]  = new TH2F("ESD_TruePiMiPiZeroSameMotherFromOther_InvMassPt","ESD_TruePiMiPiZeroSameMotherFromOther_InvMassPt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTruePiMiPiZeroSameMotherFromOtherlInvMassPt[iCut]->GetXaxis()->SetTitle(Form("M_{#pi^{-}%s} (GeV/c^{2})",NameNDMLatex.Data()));
          fHistoTruePiMiPiZeroSameMotherFromOtherlInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          if (fIsMC>1) fHistoTruePiMiPiZeroSameMotherFromOtherlInvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePiMiPiZeroSameMotherFromOtherlInvMassPt[iCut]);

          fHistoTruePiPlPiZeroSameMotherFromEtaInvMassPt[iCut]  = new TH2F("ESD_TruePiPlPiZeroSameMotherFromEta_InvMassPt","ESD_TruePiPlPiZeroSameMotherFromEta_InvMassPt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTruePiPlPiZeroSameMotherFromEtaInvMassPt[iCut]->GetXaxis()->SetTitle(Form("M_{#pi^{+}%s} (GeV/c^{2})",NameNDMLatex.Data()));
          fHistoTruePiPlPiZeroSameMotherFromEtaInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          if (fIsMC>1) fHistoTruePiPlPiZeroSameMotherFromEtaInvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePiPlPiZeroSameMotherFromEtaInvMassPt[iCut]);
          fHistoTruePiPlPiZeroSameMotherFromOmegaInvMassPt[iCut] = new TH2F("ESD_TruePiPlPiZeroSameMotherFromOmega_InvMassPt","ESD_TruePiPlPiZeroSameMotherFromOmega_InvMassPt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTruePiPlPiZeroSameMotherFromOmegaInvMassPt[iCut]->GetXaxis()->SetTitle(Form("M_{#pi^{+}%s} (GeV/c^{2})",NameNDMLatex.Data()));
          fHistoTruePiPlPiZeroSameMotherFromOmegaInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          if (fIsMC>1) fHistoTruePiPlPiZeroSameMotherFromOmegaInvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePiPlPiZeroSameMotherFromOmegaInvMassPt[iCut]);
          fHistoTruePiPlPiZeroSameMotherFromRhoInvMassPt[iCut]  = new TH2F("ESD_TruePiPlPiZeroSameMotherFromRho_InvMassPt","ESD_TruePiPlPiZeroSameMotherFromRho_InvMassPt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTruePiPlPiZeroSameMotherFromRhoInvMassPt[iCut]->GetXaxis()->SetTitle(Form("M_{#pi^{+}%s} (GeV/c^{2})",NameNDMLatex.Data()));
          fHistoTruePiMiPiZeroSameMotherFromRhoInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          if (fIsMC>1) fHistoTruePiMiPiZeroSameMotherFromRhoInvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePiPlPiZeroSameMotherFromRhoInvMassPt[iCut]);

          fHistoTruePiPlPiZeroSameMotherFromK0lInvMassPt[iCut]  = new TH2F("ESD_TruePiPlPiZeroSameMotherFromK0l_InvMassPt","ESD_TruePiPlPiZeroSameMotherFromK0l_InvMassPt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTruePiPlPiZeroSameMotherFromK0lInvMassPt[iCut]->GetXaxis()->SetTitle(Form("M_{#pi^{+}%s} (GeV/c^{2})",NameNDMLatex.Data()));
          fHistoTruePiPlPiZeroSameMotherFromK0lInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          if (fIsMC>1) fHistoTruePiPlPiZeroSameMotherFromK0lInvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePiPlPiZeroSameMotherFromK0lInvMassPt[iCut]);

          fHistoTruePiPlPiZeroSameMotherFromOtherInvMassPt[iCut]  = new TH2F("ESD_TruePiPlPiZeroSameMotherFromOther_InvMassPt","ESD_TruePiPlPiZeroSameMotherFromOther_InvMassPt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTruePiPlPiZeroSameMotherFromOtherInvMassPt[iCut]->GetXaxis()->SetTitle(Form("M_{#pi^{+}%s} (GeV/c^{2})",NameNDMLatex.Data()));
          fHistoTruePiPlPiZeroSameMotherFromOtherInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          if (fIsMC>1) fHistoTruePiPlPiZeroSameMotherFromOtherInvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePiPlPiZeroSameMotherFromOtherInvMassPt[iCut]);

          fHistoTruePiPlPiMiNDMPureCombinatoricalInvMassPt[iCut]  = new TH2F("ESD_TruePiPlPiMiNDMPureCombinatorical_InvMassPt","ESD_TruePiPlPiMiNDMPureCombinatorical_InvMassPt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTruePiPlPiMiNDMPureCombinatoricalInvMassPt[iCut]->GetXaxis()->SetTitle(Form("M_{#pi^{+}#pi^{-}%s} (GeV/c^{2})",NameNDMLatex.Data()));
          fHistoTruePiPlPiMiNDMPureCombinatoricalInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          if (fIsMC>1) fHistoTruePiPlPiMiNDMPureCombinatoricalInvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePiPlPiMiNDMPureCombinatoricalInvMassPt[iCut]);

          fHistoTruePiPlPiMiNDMContamination_Pi0_InvMassPt[iCut]  = new TH2F("ESD_TruePiPlPiMiNDMContamination_Pi0_InvMassPt","ESD_TruePiPlPiMiNDMContamination_Pi0_InvMassPt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTruePiPlPiMiNDMContamination_Pi0_InvMassPt[iCut]->GetXaxis()->SetTitle(Form("M_{#pi^{+}#pi^{-}%s} (GeV/c^{2})",NameNDMLatex.Data()));
          fHistoTruePiPlPiMiNDMContamination_Pi0_InvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          if (fIsMC>1) fHistoTruePiPlPiMiNDMContamination_Pi0_InvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePiPlPiMiNDMContamination_Pi0_InvMassPt[iCut]);
          fHistoTruePiPlPiMiNDMContamination_PiPl_InvMassPt[iCut]  = new TH2F("ESD_TruePiPlPiMiNDMContamination_PiPl_InvMassPt","ESD_TruePiPlPiMiNDMContamination_PiPl_InvMassPt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTruePiPlPiMiNDMContamination_PiPl_InvMassPt[iCut]->GetXaxis()->SetTitle(Form("M_{#pi^{+}#pi^{-}%s} (GeV/c^{2})",NameNDMLatex.Data()));
          fHistoTruePiPlPiMiNDMContamination_PiPl_InvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          if (fIsMC>1) fHistoTruePiPlPiMiNDMContamination_PiPl_InvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePiPlPiMiNDMContamination_PiPl_InvMassPt[iCut]);
          fHistoTruePiPlPiMiNDMContamination_PiMi_InvMassPt[iCut]  = new TH2F("ESD_TruePiPlPiMiNDMContamination_PiMi_InvMassPt","ESD_TruePiPlPiMiNDMContamination_PiMi_InvMassPt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTruePiPlPiMiNDMContamination_PiMi_InvMassPt[iCut]->GetXaxis()->SetTitle(Form("M_{#pi^{+}#pi^{-}%s} (GeV/c^{2})",NameNDMLatex.Data()));
          fHistoTruePiPlPiMiNDMContamination_PiMi_InvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          if (fIsMC>1) fHistoTruePiPlPiMiNDMContamination_PiMi_InvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePiPlPiMiNDMContamination_PiMi_InvMassPt[iCut]);
          fHistoTruePiPlPiMiNDMContamination_Crosscheck_InvMassPt[iCut]  = new TH2F("ESD_TruePiPlPiMiNDMContamination_Crosscheck_InvMassPt","ESD_TruePiPlPiMiNDMContamination_Crosscheck_InvMassPt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTruePiPlPiMiNDMContamination_Crosscheck_InvMassPt[iCut]->GetXaxis()->SetTitle(Form("M_{#pi^{+}#pi^{-}%s} (GeV/c^{2})",NameNDMLatex.Data()));
          fHistoTruePiPlPiMiNDMContamination_Crosscheck_InvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          if (fIsMC>1) fHistoTruePiPlPiMiNDMContamination_Crosscheck_InvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePiPlPiMiNDMContamination_Crosscheck_InvMassPt[iCut]);
          fHistoTruePiPlPiMiNDMContamination_multipel_InvMassPt[iCut]  = new TH2F("ESD_TruePiPlPiMiNDMContamination_multipel_InvMassPt","ESD_TruePiPlPiMiNDMContamination_multipel_InvMassPt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTruePiPlPiMiNDMContamination_multipel_InvMassPt[iCut]->GetXaxis()->SetTitle(Form("M_{#pi^{+}#pi^{-}%s} (GeV/c^{2})",NameNDMLatex.Data()));
          fHistoTruePiPlPiMiNDMContamination_multipel_InvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          if (fIsMC>1) fHistoTruePiPlPiMiNDMContamination_multipel_InvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePiPlPiMiNDMContamination_multipel_InvMassPt[iCut]);
          //---------------------------------------
          //Dalitz All Pt
          if (enableDalitzAllPt){
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM[iCut]    = new TH2F("ESD_TrueMotherPiPlPiMiNDM_DalitzPlot_FixedPzNDM_Pos","ESD_TrueMotherPiPlPiMiNDM_DalitzPlot_FixedPzNDM_Pos", HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz, HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz);
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM[iCut]->Sumw2();
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-}} (GeV/c^{2})");
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM[iCut]->GetYaxis()->SetTitle(Form("M_{#pi^{+} %s} (GeV/c^{2})", NameNDMLatex.Data()));
            fTrueList[iCut]->Add(fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM[iCut]);

            fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM[iCut]    = new TH2F("ESD_TrueMotherPiPlPiMiNDM_DalitzPlot_FixedPzNDM_Neg","ESD_TrueMotherPiPlPiMiNDM_DalitzPlot_FixedPzNDM_Neg", HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz, HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz);
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM[iCut]->Sumw2();
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-}} (GeV/c^{2})");
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM[iCut]->GetYaxis()->SetTitle(Form("M_{#pi^{-} %s} (GeV/c^{2})", NameNDMLatex.Data()));
            fTrueList[iCut]->Add(fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM[iCut]);

            fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM[iCut]    = new TH2F("ESD_TrueMotherPiPlPiMiNDM_DalitzPlot_SubNDM_Pos","ESD_TrueMotherPiPlPiMiNDM_DalitzPlot_SubNDM_Pos", HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz, HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz);
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM[iCut]->Sumw2();
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-}} (GeV/c^{2})");
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM[iCut]->GetYaxis()->SetTitle(Form("M_{#pi^{+} %s} (GeV/c^{2})", NameNDMLatex.Data()));
            fTrueList[iCut]->Add(fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM[iCut]);

            fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM[iCut]    = new TH2F("ESD_TrueMotherPiPlPiMiNDM_DalitzPlot_SubNDM_Neg","ESD_TrueMotherPiPlPiMiNDM_DalitzPlot_SubNDM_Neg", HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz, HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz);
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM[iCut]->Sumw2();
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-}} (GeV/c^{2})");
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM[iCut]->GetYaxis()->SetTitle(Form("M_{#pi^{-} %s} (GeV/c^{2})", NameNDMLatex.Data()));
            fTrueList[iCut]->Add(fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM[iCut]);
          }
          //---------------------------------------
          //Dalitz Low Pt
          if (enableDalitzLowPt){
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM_LowPt[iCut]    = new TH2F("ESD_TrueMotherPiPlPiMiNDM_DalitzPlot_FixedPzNDM_Pos_LowPt", Form("ESD_TrueMotherPiPlPiMiNDM_DalitzPlot_FixedPzNDM_Pos_LowPt (%3.1f to %3.1f)Gev/c", HistoDalitzPtRangeMin_LowPt, HistoDalitzPtRangeMax_LowPt), HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz, HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz);
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM_LowPt[iCut]->Sumw2();
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM_LowPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-}} (GeV/c^{2})");
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM_LowPt[iCut]->GetYaxis()->SetTitle(Form("M_{#pi^{+} %s} (GeV/c^{2})", NameNDMLatex.Data()));
            fTrueList[iCut]->Add(fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM_LowPt[iCut]);

            fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM_LowPt[iCut]    = new TH2F("ESD_TrueMotherPiPlPiMiNDM_DalitzPlot_FixedPzNDM_Neg_LowPt", Form("ESD_TrueMotherPiPlPiMiNDM_DalitzPlot_FixedPzNDM_Neg_LowPt (%3.1f to %3.1f)Gev/c", HistoDalitzPtRangeMin_LowPt, HistoDalitzPtRangeMax_LowPt), HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz, HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz);
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM_LowPt[iCut]->Sumw2();
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM_LowPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-}} (GeV/c^{2})");
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM_LowPt[iCut]->GetYaxis()->SetTitle(Form("M_{#pi^{-} %s} (GeV/c^{2})", NameNDMLatex.Data()));
            fTrueList[iCut]->Add(fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM_LowPt[iCut]);

            fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM_LowPt[iCut]    = new TH2F("ESD_TrueMotherPiPlPiMiNDM_DalitzPlot_SubNDM_Pos_LowPt", Form("ESD_TrueMotherPiPlPiMiNDM_DalitzPlot_SubNDM_Pos_LowPt (%3.1f to %3.1f)Gev/c", HistoDalitzPtRangeMin_LowPt, HistoDalitzPtRangeMax_LowPt), HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz, HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz);
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM_LowPt[iCut]->Sumw2();
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM_LowPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-}} (GeV/c^{2})");
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM_LowPt[iCut]->GetYaxis()->SetTitle(Form("M_{#pi^{+} %s} (GeV/c^{2})", NameNDMLatex.Data()));
            fTrueList[iCut]->Add(fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM_LowPt[iCut]);

            fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM_LowPt[iCut]    = new TH2F("ESD_TrueMotherPiPlPiMiNDM_DalitzPlot_SubNDM_Neg_LowPt", Form("ESD_TrueMotherPiPlPiMiNDM_DalitzPlot_SubNDM_Neg_LowPt (%3.1f to %3.1f)Gev/c", HistoDalitzPtRangeMin_LowPt, HistoDalitzPtRangeMax_LowPt), HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz, HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz);
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM_LowPt[iCut]->Sumw2();
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM_LowPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-}} (GeV/c^{2})");
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM_LowPt[iCut]->GetYaxis()->SetTitle(Form("M_{#pi^{-} %s} (GeV/c^{2})", NameNDMLatex.Data()));
            fTrueList[iCut]->Add(fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM_LowPt[iCut]);
          }
          //---------------------------------------
          //Dalitz Mid Pt
          if (enableDalitzMidPt){
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM_MidPt[iCut]    = new TH2F("ESD_TrueMotherPiPlPiMiNDM_DalitzPlot_FixedPzNDM_Pos_MidPt", Form("ESD_TrueMotherPiPlPiMiNDM_DalitzPlot_FixedPzNDM_Pos_MidPt (%3.1f to %3.1f)Gev/c", HistoDalitzPtRangeMin_MidPt, HistoDalitzPtRangeMax_MidPt), HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz, HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz);
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM_MidPt[iCut]->Sumw2();
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM_MidPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-}} (GeV/c^{2})");
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM_MidPt[iCut]->GetYaxis()->SetTitle(Form("M_{#pi^{+} %s} (GeV/c^{2})", NameNDMLatex.Data()));
            fTrueList[iCut]->Add(fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM_MidPt[iCut]);

            fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM_MidPt[iCut]    = new TH2F("ESD_TrueMotherPiPlPiMiNDM_DalitzPlot_FixedPzNDM_Neg_MidPt", Form("ESD_TrueMotherPiPlPiMiNDM_DalitzPlot_FixedPzNDM_Neg_MidPt (%3.1f to %3.1f)Gev/c", HistoDalitzPtRangeMin_MidPt, HistoDalitzPtRangeMax_MidPt), HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz, HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz);
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM_MidPt[iCut]->Sumw2();
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM_MidPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-}} (GeV/c^{2})");
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM_MidPt[iCut]->GetYaxis()->SetTitle(Form("M_{#pi^{-} %s} (GeV/c^{2})", NameNDMLatex.Data()));
            fTrueList[iCut]->Add(fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM_MidPt[iCut]);

            fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM_MidPt[iCut]    = new TH2F("ESD_TrueMotherPiPlPiMiNDM_DalitzPlot_SubNDM_Pos_MidPt", Form("ESD_TrueMotherPiPlPiMiNDM_DalitzPlot_SubNDM_Pos_MidPt (%3.1f to %3.1f)Gev/c", HistoDalitzPtRangeMin_MidPt, HistoDalitzPtRangeMax_MidPt), HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz, HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz);
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM_MidPt[iCut]->Sumw2();
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM_MidPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-}} (GeV/c^{2})");
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM_MidPt[iCut]->GetYaxis()->SetTitle(Form("M_{#pi^{+} %s} (GeV/c^{2})", NameNDMLatex.Data()));
            fTrueList[iCut]->Add(fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM_MidPt[iCut]);

            fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM_MidPt[iCut]    = new TH2F("ESD_TrueMotherPiPlPiMiNDM_DalitzPlot_SubNDM_Neg_MidPt", Form("ESD_TrueMotherPiPlPiMiNDM_DalitzPlot_SubNDM_Neg_MidPt (%3.1f to %3.1f)Gev/c", HistoDalitzPtRangeMin_MidPt, HistoDalitzPtRangeMax_MidPt), HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz, HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz);
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM_MidPt[iCut]->Sumw2();
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM_MidPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-}} (GeV/c^{2})");
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM_MidPt[iCut]->GetYaxis()->SetTitle(Form("M_{#pi^{-} %s} (GeV/c^{2})", NameNDMLatex.Data()));
            fTrueList[iCut]->Add(fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM_MidPt[iCut]);
          }
          //---------------------------------------
          //Dalitz High Pt
          if (enableDalitzHighPt){
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM_HighPt[iCut]    = new TH2F("ESD_TrueMotherPiPlPiMiNDM_DalitzPlot_FixedPzNDM_Pos_HighPt", Form("ESD_TrueMotherPiPlPiMiNDM_DalitzPlot_FixedPzNDM_Pos_HighPt (%3.1f to %3.1f)Gev/c", HistoDalitzPtRangeMin_HighPt, HistoDalitzPtRangeMax_HighPt), HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz, HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz);
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM_HighPt[iCut]->Sumw2();
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM_HighPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-}} (GeV/c^{2})");
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM_HighPt[iCut]->GetYaxis()->SetTitle(Form("M_{#pi^{+} %s} (GeV/c^{2})", NameNDMLatex.Data()));
            fTrueList[iCut]->Add(fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM_HighPt[iCut]);

            fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM_HighPt[iCut]    = new TH2F("ESD_TrueMotherPiPlPiMiNDM_DalitzPlot_FixedPzNDM_Neg_HighPt", Form("ESD_TrueMotherPiPlPiMiNDM_DalitzPlot_FixedPzNDM_Neg_HighPt (%3.1f to %3.1f)Gev/c", HistoDalitzPtRangeMin_HighPt, HistoDalitzPtRangeMax_HighPt), HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz, HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz);
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM_HighPt[iCut]->Sumw2();
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM_HighPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-}} (GeV/c^{2})");
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM_HighPt[iCut]->GetYaxis()->SetTitle(Form("M_{#pi^{-} %s} (GeV/c^{2})", NameNDMLatex.Data()));
            fTrueList[iCut]->Add(fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM_HighPt[iCut]);

            fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM_HighPt[iCut]    = new TH2F("ESD_TrueMotherPiPlPiMiNDM_DalitzPlot_SubNDM_Pos_HighPt", Form("ESD_TrueMotherPiPlPiMiNDM_DalitzPlot_SubNDM_Pos_HighPt (%3.1f to %3.1f)Gev/c", HistoDalitzPtRangeMin_HighPt, HistoDalitzPtRangeMax_HighPt), HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz, HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz);
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM_HighPt[iCut]->Sumw2();
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM_HighPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-}} (GeV/c^{2})");
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM_HighPt[iCut]->GetYaxis()->SetTitle(Form("M_{#pi^{+} %s} (GeV/c^{2})", NameNDMLatex.Data()));
            fTrueList[iCut]->Add(fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM_HighPt[iCut]);

            fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM_HighPt[iCut]    = new TH2F("ESD_TrueMotherPiPlPiMiNDM_DalitzPlot_SubNDM_Neg_HighPt", Form("ESD_TrueMotherPiPlPiMiNDM_DalitzPlot_SubNDM_Neg_HighPt (%3.1f to %3.1f)Gev/c", HistoDalitzPtRangeMin_HighPt, HistoDalitzPtRangeMax_HighPt), HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz, HistoNMassBins, HistoMassRangeDalitzMin, HistoMassRangeDalitz);
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM_HighPt[iCut]->Sumw2();
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM_HighPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-}} (GeV/c^{2})");
            fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM_HighPt[iCut]->GetYaxis()->SetTitle(Form("M_{#pi^{-} %s} (GeV/c^{2})", NameNDMLatex.Data()));
            fTrueList[iCut]->Add(fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM_HighPt[iCut]);
          }

          if(fDoMesonQA>1){
            fTrueTreeList[iCut]                               = new TList();
            fTrueTreeList[iCut]->SetName(nameTrueRecTTreeList.Data());
            fTrueTreeList[iCut]->SetOwner(kTRUE);
            fCutFolder[iCut]->Add(fTrueTreeList[iCut]);

            fTreePiPiSameMother[iCut]                         = new TTree("TreePiPiSameMother","TreePiPiSameMother");
            fTreePiPiSameMother[iCut]->Branch("fCasePiPi", &fCasePiPi, "fCasePiPi/S");
            fTreePiPiSameMother[iCut]->Branch("fSamePiPiMotherID", &fSamePiPiMotherID, "fSamePiPiMotherID/F");
            fTreePiPiSameMother[iCut]->Branch("fSamePiPiMotherInvMass", &fSamePiPiMotherInvMass, "fSamePiPiMotherInvMass/F");
            fTreePiPiSameMother[iCut]->Branch("fSamePiPiMotherPt", &fSamePiPiMotherPt, "fSamePiPiMotherPt/F");
            fTrueTreeList[iCut]->Add(fTreePiPiSameMother[iCut]);

            fTreePiPiPiSameMother[iCut]                         = new TTree("TreePiPiPiSameMother","TreePiPiPiSameMother");
            fTreePiPiPiSameMother[iCut]->Branch("fSamePiPiPiMotherID", &fSamePiPiPiMotherID, "fSamePiPiPiMotherID/F");
            fTreePiPiPiSameMother[iCut]->Branch("fSamePiPiPiMotherInvMass", &fSamePiPiPiMotherInvMass, "fSamePiPiPiMotherInvMass/F");
            fTreePiPiPiSameMother[iCut]->Branch("fSamePiPiPiMotherPt", &fSamePiPiPiMotherPt, "fSamePiPiPiMotherPt/F");
            fTrueTreeList[iCut]->Add(fTreePiPiPiSameMother[iCut]);

            fTreeEventInfoHNM[iCut]                         = new TTree("TreeEventInfoHNM","TreeEventInfoHNM");
            fTreeEventInfoHNM[iCut]->Branch("fV0MultiplicityHNMEvent", &fV0MultiplicityHNMEvent, "fV0MultiplicityHNMEvent/F");
            fTreeEventInfoHNM[iCut]->Branch("fTrackMultiplicityHNMEvent", &fTrackMultiplicityHNMEvent, "fTrackMultiplicityHNMEvent/F");
            fTreeEventInfoHNM[iCut]->Branch("fZVertexHNMEvent", &fZVertexHNMEvent, "fZVertexHNMEvent/F");
            fTreeEventInfoHNM[iCut]->Branch("fPtHNM", &fPtHNM, "fPtHNM/F");
            fTrueTreeList[iCut]->Add(fTreeEventInfoHNM[iCut]);
          }
        }
      }
    }
  }

  fVectorDoubleCountTruePi0s.clear();
  fVectorDoubleCountTrueHNMs.clear();
  fVectorDoubleCountTrueConvGammas.clear();

  InitBack(); // Init Background Handler

  fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
  if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader

  if(fV0Reader){
    if((AliConvEventCuts*)fV0Reader->GetEventCuts()){
      if(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms()){
        fOutputContainer->Add(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms());
      }
    }

    if((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts()){
      if(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms()){
        fOutputContainer->Add(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms());
      }
    }

  }

  for(Int_t iMatcherTask = 0; iMatcherTask < 3; iMatcherTask++){
    AliCaloTrackMatcher* temp = (AliCaloTrackMatcher*) (AliAnalysisManager::GetAnalysisManager()->GetTask(Form("CaloTrackMatcher_%i_%i",iMatcherTask,fTrackMatcherRunningMode)));
    if(temp && (!fDoLightOutput)) fOutputContainer->Add(temp->GetCaloTrackMatcherHistograms());
  }

  fPionSelector=(AliPrimaryPionSelector*)AliAnalysisManager::GetAnalysisManager()->GetTask(fPionSelectorName.Data());
  if(!fPionSelector){printf("Error: No PionSelector");return;} // GetV0Reader

  if( fPionSelector && (!fDoLightOutput)){
    if ( ((AliPrimaryPionCuts*)fPionSelector->GetPrimaryPionCuts())->GetCutHistograms() ){
      fOutputContainer->Add( ((AliPrimaryPionCuts*)fPionSelector->GetPrimaryPionCuts())->GetCutHistograms() );
    }
  }

  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if( fEventCutArray) {
      if( ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms() ) {
        fCutFolder[iCut]->Add( ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms());
      }
    }

    if( fPionCutArray){
      if( ((AliPrimaryPionCuts*)fPionCutArray->At(iCut))->GetCutHistograms() ) {
        fCutFolder[iCut]->Add( ((AliPrimaryPionCuts*)fPionCutArray->At(iCut))->GetCutHistograms() );
      }
    }
    if (fNDMRecoMode < 2){
      if( fGammaCutArray ) {
        if( ((AliConversionPhotonCuts*)fGammaCutArray->At(iCut))->GetCutHistograms() ) {
          fCutFolder[iCut]->Add( ((AliConversionPhotonCuts*)fGammaCutArray->At(iCut))->GetCutHistograms()  );
        }
      }
    }
    if (fNDMRecoMode > 0){
      if( fClusterCutArray ) {
        if( ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutHistograms() ) {
          fCutFolder[iCut]->Add( ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutHistograms()  );
        }
      }
    }
    if( fNeutralDecayMesonCutArray ) {
      if( ((AliConversionMesonCuts*)fNeutralDecayMesonCutArray->At(iCut))->GetCutHistograms() ) {
        fCutFolder[iCut]->Add( ((AliConversionMesonCuts*)fNeutralDecayMesonCutArray->At(iCut))->GetCutHistograms());
      }
    }
    if( fMesonCutArray ) {
      if( ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms() ) {
        fCutFolder[iCut]->Add( ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms());
      }
    }
  }

  PostData(1, fOutputContainer);

}

//______________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::UserExec(Option_t *){

  //
  // Execute analysis for current event
  //
  fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
  if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader

  Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();
  if(InputEvent()->IsIncompleteDAQ()==kTRUE) eventQuality = 2;  // incomplete event
  if(eventQuality == 2 || eventQuality == 3){// Event Not Accepted due to MC event missing or wrong trigger for V0ReaderV1 or because it is incomplete
    for(Int_t iCut = 0; iCut<fnCuts; iCut++){
      fHistoNEvents[iCut]->Fill(eventQuality);
      if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventQuality);
    }
    return;
  }

  fPionSelector=(AliPrimaryPionSelector*)AliAnalysisManager::GetAnalysisManager()->GetTask(fPionSelectorName.Data());
  if(!fPionSelector){printf("Error: No PionSelector");return;} // GetV0Reader

  if(fIsMC) fMCEvent     =  MCEvent();
  fInputEvent        = InputEvent();
  fReaderGammas    = fV0Reader->GetReconstructedGammas(); // Gammas from default Cut
  fSelectorNegPionIndex = fPionSelector->GetReconstructedNegPionIndex(); // Electrons from default Cut
  fSelectorPosPionIndex = fPionSelector->GetReconstructedPosPionIndex(); // Positrons from default Cut

  fNumberOfESDTracks = fV0Reader->GetNumberOfPrimaryTracks();
  //AddTaskContainers(); //Add conatiner
  if(fIsMC && fInputEvent->IsA()==AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled())){
    RelabelAODPhotonCandidates(kTRUE);    // In case of AODMC relabeling MC
    fV0Reader->RelabelAODs(kTRUE);
  }
  for(Int_t iCut = 0; iCut<fnCuts; iCut++){
    fiCut = iCut;

    Bool_t isRunningEMCALrelAna = kFALSE;
    if (fNDMRecoMode > 0){
      if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 1) isRunningEMCALrelAna = kTRUE;
    }

    Int_t eventNotAccepted = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon, isRunningEMCALrelAna);

    if(fIsMC==2){
      Float_t xsection      = -1.;
      Float_t ntrials       = -1.;
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetXSectionAndNTrials(fMCEvent,xsection,ntrials, fInputEvent );
      if((xsection==-1.) || (ntrials==-1.)) AliFatal(Form("ERROR: GetXSectionAndNTrials returned invalid xsection/ntrials, periodName from V0Reader: '%s'",fV0Reader->GetPeriodName().Data()));
      fProfileJetJetXSection[iCut]->Fill(0.,xsection);
      fHistoJetJetNTrials[iCut]->Fill("#sum{NTrials}", ntrials);
    }

    if (fIsMC > 1){
      fWeightJetJetMC       = 1;
      Float_t maxjetpt      = -1.;
      Float_t pthard = -1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetUseJetFinderForOutliers()) maxjetpt = fOutlierJetReader->GetMaxJetPt();
      Bool_t isMCJet        = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsJetJetMCEventAccepted( fMCEvent, fWeightJetJetMC ,pthard, fInputEvent, maxjetpt);
      if (fIsMC == 3){
        Double_t weightMult   = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetWeightForMultiplicity(fV0Reader->GetNumberOfPrimaryTracks());
        fWeightJetJetMC       = fWeightJetJetMC*weightMult;
      }

      if (!isMCJet){
        fHistoNEvents[iCut]->Fill(10,fWeightJetJetMC);
        if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(10);
        continue;
      }
    }

    Bool_t triggered = kTRUE;
    if(eventNotAccepted!=0){
      // 			cout << "event rejected due to wrong trigger: " <<eventNotAccepted << endl;
      fHistoNEvents[iCut]->Fill(eventNotAccepted,fWeightJetJetMC); // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1
      if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventNotAccepted);
      if (eventNotAccepted==3 && fIsMC > 0){
        triggered = kFALSE;
      }else {
        continue;
      }
    }

    if(eventQuality != 0 && triggered== kTRUE){// Event Not Accepted
      //cout << "event rejected due to: " <<eventQuality << endl;
      fHistoNEvents[iCut]->Fill(eventQuality, fWeightJetJetMC);
      if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventQuality); // Should be 0 here
      continue;
    }

    if (triggered == kTRUE) {
      fHistoNEvents[iCut]->Fill(eventQuality,fWeightJetJetMC);
      if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventQuality); // Should be 0 here
      fHistoNGoodESDTracks[iCut]->Fill(fNumberOfESDTracks,fWeightJetJetMC);
      if(!fDoLightOutput){
        fHistoSPDClusterTrackletBackground[iCut]->Fill(fInputEvent->GetMultiplicity()->GetNumberOfTracklets(),(fInputEvent->GetNumberOfITSClusters(0)+fInputEvent->GetNumberOfITSClusters(1)),fWeightJetJetMC);
      }
    }
    if(fIsMC> 0){ // Process MC Particle
      if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetSignalRejection() != 0){

        if(fInputEvent->IsA()==AliESDEvent::Class()){
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetNotRejectedParticles(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetSignalRejection(),
                                          ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetAcceptedHeader(),
                                          fMCEvent);
        } else if(fInputEvent->IsA()==AliAODEvent::Class()){
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetNotRejectedParticles(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetSignalRejection(),
                                          ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetAcceptedHeader(),
                                          fInputEvent);
        }

      }
      if(fInputEvent->IsA()==AliESDEvent::Class()) ProcessMCParticles();
      if(fInputEvent->IsA()==AliAODEvent::Class()) ProcessAODMCParticles();
    }
    // continue is called after processing of MC particle to get correct acceptance
    // and include trigger efficiency in efficiency (relevant when mimicking)
    if (triggered==kFALSE) continue;

    if (fNDMRecoMode < 2){
      ProcessConversionPhotonCandidates(); // Process this cuts conversion gammas
    }
    if (fNDMRecoMode > 0){
      ProcessCaloPhotonCandidates(); // Process this cuts calo gammas
    }

    if (fNDMRecoMode == 0 ){
      ProcessNeutralDecayMesonCandidatesPureConversions(); // Process neutral pion candidates purely from conversions
    }
    if (fNDMRecoMode == 1){
      ProcessNeutralPionCandidatesMixedConvCalo(); // Process neutral pion candidates mixed conv and calo
    }
    if (fNDMRecoMode == 2){
      ProcessNeutralPionCandidatesPureCalo(); // Process neutral pion candidates purely from calo
    }

    if(fInputEvent->IsA()==AliESDEvent::Class()) ProcessPionCandidates(); // Process this cuts gammas
    if(fInputEvent->IsA()==AliAODEvent::Class()) ProcessPionCandidatesAOD();

    //CalculateMesonCandidates();

    if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->DoBGCalculation()){
      if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseLikeSignMixing()){
        CalculateBackground(5);
      } else{
        CalculateBackground(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetBackgroundMode());
      }
    }
    UpdateEventByEventData();

    fVectorDoubleCountTruePi0s.clear();
    fVectorDoubleCountTrueHNMs.clear();
    fVectorDoubleCountTrueConvGammas.clear();

    fGoodConvGammas->Clear();
    fClusterCandidates->Clear();
    fNeutralDecayParticleCandidates->Clear();
    if(fNeutralDecayParticleSidebandCandidates) fNeutralDecayParticleSidebandCandidates->Clear();
    if(fNeutralDecayParticleSwappCandidates) fNeutralDecayParticleSwappCandidates->Clear();
    fPosPionCandidates->Clear();
    fNegPionCandidates->Clear();
  }

  fSelectorNegPionIndex.clear();
  fSelectorPosPionIndex.clear();

  if( fIsMC && fInputEvent->IsA()==AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled())){
    RelabelAODPhotonCandidates(kFALSE); // Back to ESDMC Label
    fV0Reader->RelabelAODs(kFALSE);
  }

  PostData( 1, fOutputContainer );
}
//________________________________________________________________________
Bool_t AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::Notify(){
  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetPeriodEnum() == AliConvEventCuts::kNoPeriod && ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum() != AliConvEventCuts::kNoPeriod){
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->SetPeriodEnumExplicit(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum());
    } else if (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetPeriodEnum() == AliConvEventCuts::kNoPeriod ){
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->SetPeriodEnum(fV0Reader->GetPeriodName());
    }

    if( !((AliConvEventCuts*)fEventCutArray->At(iCut))->GetDoEtaShift() ){
      if(!fDoLightOutput){
        fProfileEtaShift[iCut]->Fill(0.,0.);
      }
      continue; // No Eta Shift requested, continue
    }
    if( ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift() == 0.0){ // Eta Shift requested but not set, get shift automatically
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCorrectEtaShiftFromPeriod();
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once
      ((AliPrimaryPionCuts*)fPionCutArray->At(iCut))->SetEtaShift( ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift() );
      if(!fDoLightOutput){
        fProfileEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
      }
      continue;
    } else {
      printf(" Eta t PiPlusPiMinus Gamma Task %s :: Eta Shift Manually Set to %f \n\n",
      (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber()).Data(),((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift());
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once
      ((AliPrimaryPionCuts*)fPionCutArray->At(iCut))->SetEtaShift( ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift() );
      if(!fDoLightOutput){
        fProfileEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
      }
    }
  }
  return kTRUE;
}


void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::Terminate(const Option_t *){
///Grid
}


//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::ProcessCaloPhotonCandidates()
{

  Int_t nclus = 0;
  nclus = fInputEvent->GetNumberOfCaloClusters();

  // 	cout << nclus << endl;

  if(nclus == 0)	return;
    // plotting histograms on cell/tower level, only if extendedMatchAndQA > 1
  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FillHistogramsExtendedQA(fInputEvent,fIsMC);

  // match tracks to clusters
  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchTracksToClusters(fInputEvent,fWeightJetJetMC,kTRUE, fMCEvent);


  // vertex
  Double_t vertex[3] = {0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  // Loop over EMCal clusters
  for(Long_t i = 0; i < nclus; i++){
    Double_t tempClusterWeight        = fWeightJetJetMC;
    Double_t tempPhotonWeight         = fWeightJetJetMC;

    std::unique_ptr<AliVCluster> clus;
    if(fInputEvent->IsA()==AliESDEvent::Class()) clus = std::unique_ptr<AliVCluster>(new AliESDCaloCluster(*(AliESDCaloCluster*)fInputEvent->GetCaloCluster(i)));
    else if(fInputEvent->IsA()==AliAODEvent::Class()) clus = std::unique_ptr<AliVCluster>(new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(i)));

    if (!clus) continue;
    // Set the jetjet weight to 1 in case the cluster orignated from the minimum bias header
    if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4){
      if( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(clus->GetLabelAt(0), fMCEvent, fInputEvent) == 2)
        tempClusterWeight = 1;
    }

    if(!((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelected(clus.get(),fInputEvent,fMCEvent,fIsMC, tempClusterWeight,i)){
      continue;
    }
    // TLorentzvector with cluster
    TLorentzVector clusterVector;
    clus->GetMomentum(clusterVector,vertex);

    TLorentzVector tmpvec;
    tmpvec.SetPxPyPzE(clusterVector.Px(),clusterVector.Py(),clusterVector.Pz(),clusterVector.E());

    // convert to AODConversionPhoton
    AliAODConversionPhoton *PhotonCandidate=new AliAODConversionPhoton(&tmpvec);
    if(!PhotonCandidate){ continue;}

    // Flag Photon as CaloPhoton
    PhotonCandidate->SetIsCaloPhoton(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType());
    PhotonCandidate->SetCaloClusterRef(i);
    // get MC label
    if(fIsMC>0){
      Int_t* mclabelsCluster = clus->GetLabels();
      PhotonCandidate->SetNCaloPhotonMCLabels(clus->GetNLabels());
      // cout << clus->GetNLabels() << endl;
      if (clus->GetNLabels()>0){
        for (Int_t k =0; k< (Int_t)clus->GetNLabels(); k++){
          PhotonCandidate->SetCaloPhotonMCLabel(k,mclabelsCluster[k]);
          // Int_t pdgCode = fMCEvent->Particle(mclabelsCluster[k])->GetPdgCode();
          // cout << "label " << k << "\t" << mclabelsCluster[k] << " pdg code: " << pdgCode << endl;
        }
      }
    }

    fIsFromDesiredHeader            = kTRUE;
    fIsOverlappingWithOtherHeader   = kFALSE;
    //TString periodName         = fV0Reader->GetPeriodName();
    // test whether largest contribution to cluster orginates in added signals
    if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() > 0){
      if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 0){
        fIsFromDesiredHeader = kFALSE;
      }
      if (clus->GetNLabels()>1){
        Int_t* mclabelsCluster = clus->GetLabels();
        for (Int_t l = 1; l < (Int_t)clus->GetNLabels(); l++ ){
          if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(mclabelsCluster[l], fMCEvent, fInputEvent) == 0){
            fIsOverlappingWithOtherHeader = kTRUE;
          }
        }
      }
    }

    if ( (fIsFromDesiredHeader && !fIsOverlappingWithOtherHeader && !fAllowOverlapHeaders) || (fIsFromDesiredHeader && fAllowOverlapHeaders)){
      fHistoClusterGammaPt[fiCut]->Fill(PhotonCandidate->Pt(),tempPhotonWeight);
      fHistoClusterGammaEta[fiCut]->Fill(PhotonCandidate->Eta(),tempPhotonWeight);
      fHistoClusterGammaE[fiCut]->Fill(PhotonCandidate->E(),tempPhotonWeight);
      if(fIsMC>0){
        if(fInputEvent->IsA()==AliESDEvent::Class()){
          ProcessTrueCaloPhotonCandidates(PhotonCandidate);
        }else {
          ProcessTrueCaloPhotonCandidatesAOD(PhotonCandidate);
        }
      }
      fClusterCandidates->Add(PhotonCandidate); // if no second loop is required add to events good gammas
    } else{
      delete PhotonCandidate;
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::ProcessTrueCaloPhotonCandidates(AliAODConversionPhoton *TruePhotonCandidate)
{
  TParticle *Photon = nullptr;
  if (TruePhotonCandidate->GetIsCaloPhoton() == 0) AliFatal("CaloPhotonFlag has not been set task will abort");
        // fHistoTrueNLabelsInClus[fiCut]->Fill(TruePhotonCandidate->GetNCaloPhotonMCLabels());

  const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX 	= primVtxMC->GetX();
  Double_t mcProdVtxY 	= primVtxMC->GetY();
  Double_t mcProdVtxZ 	= primVtxMC->GetZ();

  if (TruePhotonCandidate->GetNCaloPhotonMCLabels()>0)Photon = fMCEvent->Particle(TruePhotonCandidate->GetCaloPhotonMCLabel(0));
    else return;

  if(Photon == nullptr){
  //    cout << "no photon" << endl;
    return;
  }

        // Int_t pdgCodeParticle = Photon->GetPdgCode();
  TruePhotonCandidate->SetCaloPhotonMCFlags(fMCEvent, kFALSE);

  // True Photon
  if(!fDoLightOutput){
    Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, TruePhotonCandidate->GetCaloPhotonMCLabel(0), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
    if(isPrimary){
      if (TruePhotonCandidate->IsLargestComponentPhoton()){
        fHistoTrueClusterGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
        if (GammaIsNeutralMesonPiPlPiMiNDMDaughter(TruePhotonCandidate->GetCaloPhotonMCLabel(0))){
          fHistoTrueClusterGammaFromNeutralMesonPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
        }
      }
      if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
          fHistoTrueClusterGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
          if (GammaIsNeutralMesonPiPlPiMiNDMDaughter(TruePhotonCandidate->GetCaloPhotonMCLabel(0))){
          fHistoTrueClusterGammaFromNeutralMesonPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
        }
      }
    }
  }
  return;
}

//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::ProcessTrueCaloPhotonCandidatesAOD(AliAODConversionPhoton *TruePhotonCandidate)
{
  AliAODMCParticle *Photon = NULL;
  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));

  if (AODMCTrackArray){
    if (TruePhotonCandidate->GetIsCaloPhoton() == 0) AliFatal("CaloPhotonFlag has not been set task will abort");
  }else{
    AliInfo("AODMCTrackArray could not be loaded");
    return;
  }
    // fHistoTrueNLabelsInClus[fiCut]->Fill(TruePhotonCandidate->GetNCaloPhotonMCLabels());

  const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX 	= primVtxMC->GetX();
  Double_t mcProdVtxY 	= primVtxMC->GetY();
  Double_t mcProdVtxZ 	= primVtxMC->GetZ();

  if (TruePhotonCandidate->GetNCaloPhotonMCLabels()>0)Photon = (AliAODMCParticle*) AODMCTrackArray->At(TruePhotonCandidate->GetCaloPhotonMCLabel(0));
  else return;

  if(Photon == NULL){
    return;
  }

    // Int_t pdgCodeParticle = Photon->GetPdgCode();
  TruePhotonCandidate->SetCaloPhotonMCFlagsAOD(AODMCTrackArray, kFALSE);

  // True Photon
  if(!fDoLightOutput){
    Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD( fInputEvent, Photon, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
    if(isPrimary){
      if (TruePhotonCandidate->IsLargestComponentPhoton()){
        fHistoTrueClusterGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
        if (GammaIsNeutralMesonPiPlPiMiNDMDaughterAOD(AODMCTrackArray,TruePhotonCandidate->GetCaloPhotonMCLabel(0))){
          fHistoTrueClusterGammaFromNeutralMesonPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
        }
      }
      if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
          fHistoTrueClusterGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
          if (GammaIsNeutralMesonPiPlPiMiNDMDaughterAOD(AODMCTrackArray,TruePhotonCandidate->GetCaloPhotonMCLabel(0))){
          fHistoTrueClusterGammaFromNeutralMesonPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
        }
      }
    }
  }
  return;
}


//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::ProcessConversionPhotonCandidates(){
  Int_t nV0 = 0;
  TList GoodGammasStepOne;
  TList GoodGammasStepTwo;
  // Loop over Photon Candidates allocated by ReaderV1

  for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){
    AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(i);
    if(!PhotonCandidate) continue;

    fIsFromMBHeader = kTRUE;

    if( fMCEvent && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0 ){
      Int_t isPosFromMBHeader
        = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCEvent, fInputEvent);
      if(isPosFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
      Int_t isNegFromMBHeader
        = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCEvent,fInputEvent);
      if(isNegFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
      if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
    }

    if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelected(PhotonCandidate,fInputEvent)) continue;

    if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseElecSharingCut() &&
      !((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseToCloseV0sCut()){ // if no post reader loop is required add to events good gammas

      fGoodConvGammas->Add(PhotonCandidate);

      if(fIsFromMBHeader && (!fDoLightOutput)){
        fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt(), fWeightJetJetMC);
        fHistoConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta(), fWeightJetJetMC);
      }

      if(fMCEvent){
        if(fInputEvent->IsA()==AliESDEvent::Class())
          ProcessTrueConversionPhotonCandidates(PhotonCandidate);
        if(fInputEvent->IsA()==AliAODEvent::Class())
          ProcessTrueConversionPhotonCandidatesAOD(PhotonCandidate);
      }
    } else if(((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseElecSharingCut()){ // if Shared Electron cut is enabled, Fill array, add to step one
      ((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->FillElectonLabelArray(PhotonCandidate,nV0);
      nV0++;
      GoodGammasStepOne.Add(PhotonCandidate);
    } else if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseElecSharingCut() &&
        ((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseToCloseV0sCut()){ // shared electron is disabled, step one not needed -> step two
      GoodGammasStepTwo.Add(PhotonCandidate);
    }
  }


  if(((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseElecSharingCut()){
    for(Int_t i = 0;i<GoodGammasStepOne.GetEntries();i++){
      AliAODConversionPhoton *PhotonCandidate= (AliAODConversionPhoton*) GoodGammasStepOne.At(i);
      if(!PhotonCandidate) continue;
      fIsFromMBHeader = kTRUE;
      if(fMCEvent && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        Int_t isPosFromMBHeader
        = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCEvent,fInputEvent);
        Int_t isNegFromMBHeader
        = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCEvent,fInputEvent);
        if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
      }
      if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->RejectSharedElectronV0s(PhotonCandidate,i,GoodGammasStepOne.GetEntries())) continue;
      if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseToCloseV0sCut()){ // To Colse v0s cut diabled, step two not needed
        fGoodConvGammas->Add(PhotonCandidate);
        if(fIsFromMBHeader && (!fDoLightOutput)){
          fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt(), fWeightJetJetMC);
          fHistoConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta(), fWeightJetJetMC);
        }
        if(fMCEvent){
          if(fInputEvent->IsA()==AliESDEvent::Class())
            ProcessTrueConversionPhotonCandidates(PhotonCandidate);
          if(fInputEvent->IsA()==AliAODEvent::Class())
            ProcessTrueConversionPhotonCandidatesAOD(PhotonCandidate);
        }
      }
      else GoodGammasStepTwo.Add(PhotonCandidate); // Close v0s cut enabled -> add to list two
    }
  }
  if(((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseToCloseV0sCut()){
    for(Int_t i = 0;i<GoodGammasStepTwo.GetEntries();i++){
      AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) GoodGammasStepTwo.At(i);
      if(!PhotonCandidate) continue;

      if(fMCEvent && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        Int_t isPosFromMBHeader
        = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCEvent,fInputEvent);
        Int_t isNegFromMBHeader
        = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCEvent,fInputEvent);
        if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
      }

      if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->RejectToCloseV0s(PhotonCandidate,&GoodGammasStepTwo,i)) continue;
      fGoodConvGammas->Add(PhotonCandidate); // Add gamma to current cut TList

      if(fIsFromMBHeader && (!fDoLightOutput)){
        fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt(), fWeightJetJetMC); // Differences to old V0Reader in p_t due to conversion KF->TLorentzVector
        fHistoConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta(), fWeightJetJetMC);
      }

      if(fMCEvent){
        if(fInputEvent->IsA()==AliESDEvent::Class())
          ProcessTrueConversionPhotonCandidates(PhotonCandidate);
        if(fInputEvent->IsA()==AliAODEvent::Class())
          ProcessTrueConversionPhotonCandidatesAOD(PhotonCandidate);
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::ProcessTrueConversionPhotonCandidates(AliAODConversionPhoton *TruePhotonCandidate)
{
  // Process True Photons
  TParticle *posDaughter = TruePhotonCandidate->GetPositiveMCDaughter(fMCEvent);
  TParticle *negDaughter = TruePhotonCandidate->GetNegativeMCDaughter(fMCEvent);

  const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX 	= primVtxMC->GetX();
  Double_t mcProdVtxY 	= primVtxMC->GetY();
  Double_t mcProdVtxZ 	= primVtxMC->GetZ();


  if(posDaughter == NULL || negDaughter == NULL){
    return; // One particle does not exist
  }
  if(posDaughter->GetMother(0) != negDaughter->GetMother(0)){  // Not Same Mother == Combinatorial Bck
    return;
  }

  else if (posDaughter->GetMother(0) == -1){
    return;
  }

  if(TMath::Abs(posDaughter->GetPdgCode())!=11 || TMath::Abs(negDaughter->GetPdgCode())!=11){
    return; //One Particle is not electron
  }
  if(posDaughter->GetPdgCode()==negDaughter->GetPdgCode()){
    return; // Same Charge
  }
  if(posDaughter->GetUniqueID() != 5 || negDaughter->GetUniqueID() !=5){
    return;// check if the daughters come from a conversion
  }
  TParticle *Photon = TruePhotonCandidate->GetMCParticle(fMCEvent);
  if(Photon->GetPdgCode() != 22){
    return; // Mother is no Photon
  }
  // True Photon

  if (CheckVectorForDoubleCount(fVectorDoubleCountTrueConvGammas,posDaughter->GetMother(0)) && (!fDoLightOutput)) fHistoDoubleCountTrueConvGammaRPt[fiCut]->Fill(TruePhotonCandidate->GetConversionRadius(),TruePhotonCandidate->Pt());

  Int_t labelGamma = TruePhotonCandidate->GetMCParticleLabel(fMCEvent);
  Bool_t gammaIsPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, labelGamma, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
  if( gammaIsPrimary ){
    if( fIsFromMBHeader && (!fDoLightOutput) ){
      fHistoTrueConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC);
      if (GammaIsNeutralMesonPiPlPiMiNDMDaughter(labelGamma)){
        fHistoTrueConvGammaFromNeutralMesonPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::ProcessTrueConversionPhotonCandidatesAOD(AliAODConversionPhoton *TruePhotonCandidate)
{
  // Process True Photons
  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (AODMCTrackArray == NULL || TruePhotonCandidate == NULL){
    return;
  }
  AliAODMCParticle *posDaughter = (AliAODMCParticle*) AODMCTrackArray->At(TruePhotonCandidate->GetMCLabelPositive());
  AliAODMCParticle *negDaughter = (AliAODMCParticle*) AODMCTrackArray->At(TruePhotonCandidate->GetMCLabelNegative());

  const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX 	= primVtxMC->GetX();
  Double_t mcProdVtxY 	= primVtxMC->GetY();
  Double_t mcProdVtxZ 	= primVtxMC->GetZ();


  if(posDaughter == NULL || negDaughter == NULL) {
    return; // One particle does not exist
  }
  if(posDaughter->GetMother() != negDaughter->GetMother()){  // Not Same Mother == Combinatorial Bck
    return;
  }
  else if (posDaughter->GetMother() == -1){
    return;
  }

  if(TMath::Abs(posDaughter->GetPdgCode())!=11 || TMath::Abs(negDaughter->GetPdgCode())!=11){
    return; //One Particle is not electron
  }
  if(posDaughter->GetPdgCode()==negDaughter->GetPdgCode()){
    return; // Same Charge
  }
  if(posDaughter->GetMCProcessCode() != 5 || negDaughter->GetMCProcessCode() !=5){
    return;// check if the daughters come from a conversion
  }
  AliAODMCParticle *Photon = (AliAODMCParticle*) AODMCTrackArray->At(posDaughter->GetMother());
  if(Photon->GetPdgCode() != 22){
    return; // Mother is no Photon
}
  // True Photon

  if (CheckVectorForDoubleCount(fVectorDoubleCountTrueConvGammas,posDaughter->GetMother()) && (!fDoLightOutput)) fHistoDoubleCountTrueConvGammaRPt[fiCut]->Fill(TruePhotonCandidate->GetConversionRadius(),TruePhotonCandidate->Pt());

  // for AOD you have to ask electron for mother to get label
  Int_t labelGamma = posDaughter->GetMother();
  Bool_t gammaIsPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD( fInputEvent, Photon, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
  if( gammaIsPrimary ){
    if( fIsFromMBHeader && (!fDoLightOutput) ){
      fHistoTrueConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
      if (GammaIsNeutralMesonPiPlPiMiNDMDaughterAOD(AODMCTrackArray, labelGamma)){
        fHistoTrueConvGammaFromNeutralMesonPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::ProcessNeutralDecayMesonCandidatesPureConversions(){
  // Conversion Gammas
  if(fGoodConvGammas->GetEntries()>1){
    for(Int_t firstGammaIndex=0;firstGammaIndex<fGoodConvGammas->GetEntries()-1;firstGammaIndex++){
      AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGoodConvGammas->At(firstGammaIndex));
      if (gamma0==nullptr) continue;
      for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fGoodConvGammas->GetEntries();secondGammaIndex++){
        AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fGoodConvGammas->At(secondGammaIndex));
        //Check for same Electron ID
        if (gamma1==nullptr) continue;
        if(gamma0->GetTrackLabelPositive() == gamma1->GetTrackLabelPositive() ||
        gamma0->GetTrackLabelNegative() == gamma1->GetTrackLabelNegative() ||
        gamma0->GetTrackLabelNegative() == gamma1->GetTrackLabelPositive() ||
        gamma0->GetTrackLabelPositive() == gamma1->GetTrackLabelNegative() ) continue;

        AliAODConversionMother *NDMcand = new AliAODConversionMother(gamma0,gamma1);
        NDMcand->SetLabels(firstGammaIndex,secondGammaIndex);

        NDMcand->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());

        if(!fDoLightOutput){
            fHistoGammaGammaInvMassPtBeforeCuts[fiCut]->Fill(NDMcand->M(),NDMcand->Pt(), fWeightJetJetMC);
        }

        if( ((AliConversionMesonCuts*)fNeutralDecayMesonCutArray->At(fiCut))->DoGammaMinEnergyCut() ){
          Int_t minDaughters        = ((AliConversionMesonCuts*)fNeutralDecayMesonCutArray->At(fiCut))->GetNDaughterEnergyCut();
          Float_t minDaughterEnergy = ((AliConversionMesonCuts*)fNeutralDecayMesonCutArray->At(fiCut))->GetSingleDaughterMinE();
          if(minDaughters==1){ // at least one over threshold
             if( (gamma0->E() < minDaughterEnergy)  && (gamma1->E() < minDaughterEnergy)) continue;
          } else if (minDaughters==2){ // both over threshold
             if( (gamma0->E() < minDaughterEnergy)  || (gamma1->E() < minDaughterEnergy)) continue;
          }
        }

        if((((AliConversionMesonCuts*)fNeutralDecayMesonCutArray->At(fiCut))->MesonIsSelected(NDMcand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
          if(fIsMC){
            if(fInputEvent->IsA()==AliESDEvent::Class())
              ProcessTrueNeutralPionCandidatesPureConversions(NDMcand,gamma0,gamma1);
            if(fInputEvent->IsA()==AliAODEvent::Class())
              ProcessTrueNeutralPionCandidatesPureConversionsAOD(NDMcand,gamma0,gamma1);
          }
          if (((AliConversionMesonCuts*)fNeutralDecayMesonCutArray->At(fiCut))->MesonIsSelectedByMassCut(NDMcand, 0)){
            fNeutralDecayParticleCandidates->Add(NDMcand);

            if(!fDoLightOutput){
                fHistoGammaGammaInvMassPt[fiCut]->Fill(NDMcand->M(),NDMcand->Pt(), fWeightJetJetMC);
            }
          } else if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseSidebandMixing()) &&
                    (((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedByMassCut(NDMcand, 1))){
            fNeutralDecayParticleSidebandCandidates->Add(NDMcand);

            if(!fDoLightOutput){
                fHistoGammaGammaInvMassPt[fiCut]->Fill(NDMcand->M(),NDMcand->Pt(), fWeightJetJetMC);
            }
          } else if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseSidebandMixingBothSides()) &&
                    ((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedByMassCut(NDMcand, 2)) ||
                     ((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedByMassCut(NDMcand, 3))))){
            fNeutralDecayParticleSidebandCandidates->Add(NDMcand);

            if(!fDoLightOutput){
                fHistoGammaGammaInvMassPt[fiCut]->Fill(NDMcand->M(),NDMcand->Pt(), fWeightJetJetMC);
            }
          } else{
            delete NDMcand;
          }
        } else{
          delete NDMcand;
        }
      }
    }
  }
}


//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::ProcessNeutralPionCandidatesPureCalo(){

  if(fDoMesonQA>1 ){
    fHistoNumberClusterGamma[fiCut]->Fill(fClusterCandidates->GetEntries());
  }
  // Calo Gammas
  if(fClusterCandidates->GetEntries()>0){

    // vertex
    Double_t vertex[3] = {0};
    InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

    // Create two variables for min E cut (used later)
    Int_t minDaughters        = ((AliConversionMesonCuts*)fNeutralDecayMesonCutArray->At(fiCut))->GetNDaughterEnergyCut();
    Float_t minDaughterEnergy = ((AliConversionMesonCuts*)fNeutralDecayMesonCutArray->At(fiCut))->GetSingleDaughterMinE();

    for(Int_t firstGammaIndex=0;firstGammaIndex<fClusterCandidates->GetEntries();firstGammaIndex++){
      AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(firstGammaIndex));
      if (gamma0==nullptr) continue;

      for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fClusterCandidates->GetEntries();secondGammaIndex++){
        AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(secondGammaIndex));
        if (gamma1==nullptr) continue;

        AliAODConversionMother *NDMcand = new AliAODConversionMother(gamma0,gamma1);
        NDMcand->SetLabels(firstGammaIndex,secondGammaIndex);

        if(!fDoLightOutput){
            fHistoGammaGammaInvMassPtBeforeCuts[fiCut]->Fill(NDMcand->M(),NDMcand->Pt(), fWeightJetJetMC);
        }

        if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetBackgroundMode() == 7) && (fClusterCandidates->GetEntries()>2) ){
          TVector3 vRotationPi0;
          TLorentzVector lvRotationgamma0; // First rotated gamma
          TLorentzVector lvRotationgamma1; // Second rotated gamma

          for(int iSwapp = 0; iSwapp < ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetNumberOfSwappsForBg(); ++iSwapp){
            Bool_t acceptfirstgamma = true; // Initiate with true, when any cut "fails" it's set to false
            Bool_t acceptsecondgamma = true;  // Reset to true for every rotation

            vRotationPi0.SetX(NDMcand->Px());
            vRotationPi0.SetY(NDMcand->Py());
            vRotationPi0.SetZ(NDMcand->Pz());

            lvRotationgamma0.SetX(gamma0->Px());
            lvRotationgamma0.SetY(gamma0->Py());
            lvRotationgamma0.SetZ(gamma0->Pz());
            lvRotationgamma0.SetE(gamma0->E());

            lvRotationgamma1.SetX(gamma1->Px());
            lvRotationgamma1.SetY(gamma1->Py());
            lvRotationgamma1.SetZ(gamma1->Pz());
            lvRotationgamma1.SetE(gamma1->E());

            lvRotationgamma0.Rotate(TMath::Pi()/2.0, vRotationPi0);
            lvRotationgamma1.Rotate(TMath::Pi()/2.0, vRotationPi0);

            std::unique_ptr<AliAODConversionPhoton> gamma0swapped = std::unique_ptr<AliAODConversionPhoton>(new AliAODConversionPhoton(&lvRotationgamma0));
            std::unique_ptr<AliAODConversionPhoton> gamma1swapped = std::unique_ptr<AliAODConversionPhoton>(new AliAODConversionPhoton(&lvRotationgamma1));

            //__________________________________________________________________
            // Cuts on rotated gammas
            //  Cuts for first swapped pi0
            Int_t cellIDRotatedgamma0 = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetCaloCellIdFromEtaPhi(lvRotationgamma0.Eta(), static_cast<double>((lvRotationgamma0.Phi()<0) ? lvRotationgamma0.Phi() + TMath::Pi()*2. : lvRotationgamma0.Phi()));
            if(!(!(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedgamma0, lvRotationgamma0.Phi(), fInputEvent)) && lvRotationgamma0.E() > ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetMinClusterEnergy())) acceptfirstgamma = false;

            //  Cuts for second swapped pi0
            Int_t cellIDRotatedgamma1 = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetCaloCellIdFromEtaPhi(lvRotationgamma1.Eta(), static_cast<double>((lvRotationgamma1.Phi()<0) ? lvRotationgamma1.Phi() + TMath::Pi()*2. : lvRotationgamma1.Phi()));
            if(!(!(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedgamma1, lvRotationgamma1.Phi(), fInputEvent)) && lvRotationgamma1.E() > ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetMinClusterEnergy())) acceptsecondgamma = false;

            //__________________________________________________________________

            if(!acceptfirstgamma && !acceptsecondgamma) continue; // When no rotated gamma is accepted, dont even enter third gamma loop. Look for other pi0 rotations again instead

            for(Int_t thirdGammaIndex=0;thirdGammaIndex<fClusterCandidates->GetEntries();thirdGammaIndex++){
              AliAODConversionPhoton *gamma2=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(thirdGammaIndex));
              if (gamma2==nullptr || thirdGammaIndex==firstGammaIndex || thirdGammaIndex==secondGammaIndex) continue;

              if(acceptfirstgamma){ // Pair third gamma with first gamma only when first gamma is accepted
                std::unique_ptr<AliAODConversionMother> NDMcandswapp1(new AliAODConversionMother(gamma0swapped.get(), ((AliAODConversionPhoton*) gamma2)));
                NDMcandswapp1->SetLabels(thirdGammaIndex,firstGammaIndex);

                if( (! (((AliConversionMesonCuts*)fNeutralDecayMesonCutArray->At(fiCut))->DoGammaMinEnergyCut())) //Not using min E cut
                  || ( (minDaughters==1) && ( (gamma0swapped->E() > minDaughterEnergy)  || (gamma2->E() > minDaughterEnergy)) ) // Or require one daughter above min E
                  || ( (minDaughters==2) && ( (gamma0swapped->E() > minDaughterEnergy)  && (gamma2->E() > minDaughterEnergy)) ) ){// require both above min E
                    if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetBackgroundMode() == 7){
                      fHistoSwappingGammaGammaInvMassPt[fiCut]->Fill(NDMcandswapp1.get()->M(),NDMcandswapp1.get()->Pt(), fWeightJetJetMC);
                    }
                    if (((AliConversionMesonCuts*)fNeutralDecayMesonCutArray->At(fiCut))->MesonIsSelectedByMassCut(NDMcandswapp1.get(), 0)){
                      fNeutralDecayParticleSwappCandidates->Add(NDMcandswapp1.get());
                  }
                }
              } // End of first gamma with third gamma if


              if(acceptsecondgamma){ // Pair third gamma with second gamma only when second gamma is accepted
                std::unique_ptr<AliAODConversionMother> NDMcandswapp2(new AliAODConversionMother(gamma1swapped.get(), ((AliAODConversionPhoton*) gamma2)));
                NDMcandswapp2->SetLabels(thirdGammaIndex,secondGammaIndex);

                if( (! (((AliConversionMesonCuts*)fNeutralDecayMesonCutArray->At(fiCut))->DoGammaMinEnergyCut())) //Not using min E cut
                  || ( (minDaughters==1) && ( (gamma1swapped->E() > minDaughterEnergy)  || (gamma2->E() > minDaughterEnergy)) ) // Or require one daughter above min E
                  || ( (minDaughters==2) && ( (gamma1swapped->E() > minDaughterEnergy)  && (gamma2->E() > minDaughterEnergy)) ) ){// require both above min E
                    if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetBackgroundMode() == 7){
                        fHistoSwappingGammaGammaInvMassPt[fiCut]->Fill(NDMcandswapp2.get()->M(),NDMcandswapp2.get()->Pt(), fWeightJetJetMC);
                    }
                    if (((AliConversionMesonCuts*)fNeutralDecayMesonCutArray->At(fiCut))->MesonIsSelectedByMassCut(NDMcandswapp2.get(), 0)){
                      fNeutralDecayParticleSwappCandidates->Add(NDMcandswapp2.get());
                  }
                }
              } // End of second gamma with third gamma if
            }
          }
        }

        if( ((AliConversionMesonCuts*)fNeutralDecayMesonCutArray->At(fiCut))->DoGammaMinEnergyCut() ){
          Int_t minDaughters        = ((AliConversionMesonCuts*)fNeutralDecayMesonCutArray->At(fiCut))->GetNDaughterEnergyCut();
          Float_t minDaughterEnergy = ((AliConversionMesonCuts*)fNeutralDecayMesonCutArray->At(fiCut))->GetSingleDaughterMinE();
          if(minDaughters==1){ // at least one over threshold
             if( (gamma0->E() < minDaughterEnergy)  && (gamma1->E() < minDaughterEnergy)) continue;
          } else if (minDaughters==2){ // both over threshold
             if( (gamma0->E() < minDaughterEnergy)  || (gamma1->E() < minDaughterEnergy)) continue;
          }
        }

        if((((AliConversionMesonCuts*)fNeutralDecayMesonCutArray->At(fiCut))->MesonIsSelected(NDMcand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
          if(fIsMC){
            if(fInputEvent->IsA()==AliESDEvent::Class())
              ProcessTrueNeutralPionCandidatesPureCalo(NDMcand,gamma0,gamma1);
            if(fInputEvent->IsA()==AliAODEvent::Class())
              ProcessTrueNeutralPionCandidatesPureCaloAOD(NDMcand,gamma0,gamma1);
          }

          if (((AliConversionMesonCuts*)fNeutralDecayMesonCutArray->At(fiCut))->MesonIsSelectedByMassCut(NDMcand, 0)){
            fNeutralDecayParticleCandidates->Add(NDMcand);

            if(!fDoLightOutput){
                fHistoGammaGammaInvMassPt[fiCut]->Fill(NDMcand->M(),NDMcand->Pt(), fWeightJetJetMC);
            }
          } else if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseSidebandMixing()) &&
                    (((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedByMassCut(NDMcand, 1))){
            fNeutralDecayParticleSidebandCandidates->Add(NDMcand);

            if(!fDoLightOutput){
                fHistoGammaGammaInvMassPt[fiCut]->Fill(NDMcand->M(),NDMcand->Pt(), fWeightJetJetMC);
            }
          } else if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseSidebandMixingBothSides()) &&
                    ((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedByMassCut(NDMcand, 2)) ||
                      ((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedByMassCut(NDMcand, 3))))){
            fNeutralDecayParticleSidebandCandidates->Add(NDMcand);

            if(!fDoLightOutput){
                fHistoGammaGammaInvMassPt[fiCut]->Fill(NDMcand->M(),NDMcand->Pt(), fWeightJetJetMC);
            }
          }else {
            delete NDMcand;
          }
        } else{
          delete NDMcand;
        }
      }
    }
  }
}

//______________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::ProcessTrueNeutralPionCandidatesPureCalo( AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1)
{
  // Process True Mesons

  Bool_t isTrueNDM = kFALSE;
  Int_t gamma0MCLabel = TrueGammaCandidate0->GetCaloPhotonMCLabel(0); 	// get most probable MC label
  Int_t gamma0MotherLabel = -1;

  if(gamma0MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
    TParticle * gammaMC0 = (TParticle*)fMCEvent->Particle(gamma0MCLabel);
    if (TrueGammaCandidate0->IsLargestComponentPhoton() || TrueGammaCandidate0->IsLargestComponentElectron()){		// largest component is electro magnetic
      // get mother of interest (pi0 or eta)
      if (TrueGammaCandidate0->IsLargestComponentPhoton()){														// for photons its the direct mother
        gamma0MotherLabel=gammaMC0->GetMother(0);
      } else if (TrueGammaCandidate0->IsLargestComponentElectron()){ 												// for electrons its either the direct mother or for conversions the grandmother
                if (TrueGammaCandidate0->IsConversion() && gammaMC0->GetMother(0)>-1){
          gamma0MotherLabel=fMCEvent->Particle(gammaMC0->GetMother(0))->GetMother(0);
        } else {
          gamma0MotherLabel=gammaMC0->GetMother(0);
        }
      }
    }
  }

  if (TrueGammaCandidate1->GetIsCaloPhoton() == 0) AliFatal("CaloPhotonFlag has not been set. Aborting");

  Int_t gamma1MCLabel = TrueGammaCandidate1->GetCaloPhotonMCLabel(0); 	// get most probable MC label
  Int_t gamma1MotherLabel = -1;
  // check if
  if(gamma1MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
    // Daughters Gamma 1
    TParticle * gammaMC1 = (TParticle*)fMCEvent->Particle(gamma1MCLabel);
    if (TrueGammaCandidate1->IsLargestComponentPhoton() || TrueGammaCandidate1->IsLargestComponentElectron()){		// largest component is electro magnetic
      // get mother of interest (pi0 or eta)
      if (TrueGammaCandidate1->IsLargestComponentPhoton()){														// for photons its the direct mother
        gamma1MotherLabel=gammaMC1->GetMother(0);
      } else if (TrueGammaCandidate1->IsLargestComponentElectron()){ 												// for electrons its either the direct mother or for conversions the grandmother
                if (TrueGammaCandidate1->IsConversion() && gammaMC1->GetMother(0)>-1) gamma1MotherLabel=fMCEvent->Particle(gammaMC1->GetMother(0))->GetMother(0);
        else gamma1MotherLabel=gammaMC1->GetMother(0);
      }
    }
  }

  if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
    if(((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->GetPdgCode() == fPDGCodeNDM){
      isTrueNDM=kTRUE;
      if (CheckVectorForDoubleCount(fVectorDoubleCountTruePi0s,gamma0MotherLabel) && (!fDoLightOutput)) fHistoDoubleCountTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
    }
  }

  if(isTrueNDM){// True Pion
    Pi0Candidate->SetTrueMesonValue(1);
    Pi0Candidate->SetMCLabel(gamma0MotherLabel);
    if(!fDoLightOutput){
      fHistoTrueMotherGammaGammaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
      switch( fSelectedHeavyNeutralMeson ) {
      case 0: // ETA MESON
        if( IsEtaPiPlPiMiPiZeroDaughter(gamma0MotherLabel) )
            fHistoTrueMotherGammaGammaFromHNMInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
        break;
      case 1: // OMEGA MESON
        if( IsOmegaPiPlPiMiPiZeroDaughter(gamma0MotherLabel) )
            fHistoTrueMotherGammaGammaFromHNMInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
        break;
      case 2: // ETA PRIME MESON
        if( IsEtaPrimePiPlPiMiEtaDaughter(gamma0MotherLabel) )
            fHistoTrueMotherGammaGammaFromHNMInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
        break;
      case 3: // D0 MESON
        if( IsD0PiPlPiMiPiZeroDaughter(gamma0MotherLabel) )
            fHistoTrueMotherGammaGammaFromHNMInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
        break;
      default:
        AliError(Form("Heavy neutral meson not correctly selected (only 0,1,2,3 valid)... selected: %d",fSelectedHeavyNeutralMeson));
      }
    }
  }
}

//______________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::ProcessTrueNeutralPionCandidatesPureCaloAOD( AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1)
{
  // Process True Mesons
  // Process True Mesons
  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (AODMCTrackArray == NULL) return;

  Bool_t isTrueNDM = kFALSE;
  Int_t gamma0MCLabel = TrueGammaCandidate0->GetCaloPhotonMCLabel(0); 	// get most probable MC label
  Int_t gamma0MotherLabel = -1;

  if(gamma0MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
    AliAODMCParticle * gammaMC0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MCLabel));
    if (TrueGammaCandidate0->IsLargestComponentPhoton() || TrueGammaCandidate0->IsLargestComponentElectron()){		// largest component is electro magnetic
      // get mother of interest (pi0 or eta)
      if (TrueGammaCandidate0->IsLargestComponentPhoton()){														// for photons its the direct mother
        gamma0MotherLabel=gammaMC0->GetMother();
      } else if (TrueGammaCandidate0->IsLargestComponentElectron()){ 												// for electrons its either the direct mother or for conversions the grandmother
        if (TrueGammaCandidate0->IsConversion() && gammaMC0->GetMother()>-1){
          AliAODMCParticle * gammaGrandMotherMC0 =  static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gammaMC0->GetMother()));
          gamma0MotherLabel=gammaGrandMotherMC0->GetMother();
        } else {
          gamma0MotherLabel=gammaMC0->GetMother();
        }
      }
    }
  }

  if (TrueGammaCandidate1->GetIsCaloPhoton() == 0) AliFatal("CaloPhotonFlag has not been set. Aborting");

  Int_t gamma1MCLabel = TrueGammaCandidate1->GetCaloPhotonMCLabel(0); 	// get most probable MC label
  Int_t gamma1MotherLabel = -1;
  // check if
  if(gamma1MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
    // Daughters Gamma 1
    AliAODMCParticle *  gammaMC1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MCLabel));
    if (TrueGammaCandidate1->IsLargestComponentPhoton() || TrueGammaCandidate1->IsLargestComponentElectron()){		// largest component is electro magnetic
      // get mother of interest (pi0 or eta)
      if (TrueGammaCandidate1->IsLargestComponentPhoton()){														// for photons its the direct mother
        gamma1MotherLabel=gammaMC1->GetMother();
      } else if (TrueGammaCandidate1->IsLargestComponentElectron()){ 												// for electrons its either the direct mother or for conversions the grandmother
        if (TrueGammaCandidate1->IsConversion() && gammaMC1->GetMother()>-1){
          AliAODMCParticle * gammaGrandMotherMC1 =  static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gammaMC1->GetMother()));
          gamma1MotherLabel=gammaGrandMotherMC1->GetMother();
        }else gamma1MotherLabel=gammaMC1->GetMother();
      }
    }
  }

  if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
    if(((AliAODMCParticle*)AODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == fPDGCodeNDM){
      isTrueNDM=kTRUE;
      if (CheckVectorForDoubleCount(fVectorDoubleCountTruePi0s,gamma0MotherLabel) && (!fDoLightOutput)) fHistoDoubleCountTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
    }
  }

  if(isTrueNDM){// True Pion
    Pi0Candidate->SetTrueMesonValue(1);
    Pi0Candidate->SetMCLabel(gamma0MotherLabel);
    if(!fDoLightOutput){
      fHistoTrueMotherGammaGammaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
      switch( fSelectedHeavyNeutralMeson ) {
      case 0: // ETA MESON
        if( IsEtaPiPlPiMiPiZeroDaughterAOD(AODMCTrackArray,gamma0MotherLabel) )
            fHistoTrueMotherGammaGammaFromHNMInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
        break;
      case 1: // OMEGA MESON
        if( IsOmegaPiPlPiMiPiZeroDaughterAOD(AODMCTrackArray,gamma0MotherLabel) )
            fHistoTrueMotherGammaGammaFromHNMInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
        break;
      case 2: // ETA PRIME MESON
        if( IsEtaPrimePiPlPiMiEtaDaughterAOD(AODMCTrackArray,gamma0MotherLabel) )
            fHistoTrueMotherGammaGammaFromHNMInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
        break;
      case 3: // D0 MESON
        if( IsD0PiPlPiMiPiZeroDaughterAOD(AODMCTrackArray,gamma0MotherLabel) )
            fHistoTrueMotherGammaGammaFromHNMInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
        break;
      default:
        AliError(Form("Heavy neutral meson not correctly selected (only 0,1,2,3 valid)... selected: %d",fSelectedHeavyNeutralMeson));
      }
    }
  }
}



//______________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::ProcessTrueNeutralPionCandidatesPureConversions(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1)
{
  // Process True Mesons
  if(TrueGammaCandidate0->GetV0Index()<fInputEvent->GetNumberOfV0s()){
    Bool_t isTrueNDM = kFALSE;
    Bool_t isTruePi0Dalitz = kFALSE;
    Bool_t gamma0DalitzCand = kFALSE;
    Bool_t gamma1DalitzCand = kFALSE;
    Int_t gamma0MCLabel = TrueGammaCandidate0->GetMCParticleLabel(fMCEvent);
    Int_t gamma0MotherLabel = -1;
    Int_t motherRealLabel = -1;
    if(gamma0MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
      // Daughters Gamma 0
      TParticle * negativeMC = (TParticle*)TrueGammaCandidate0->GetNegativeMCDaughter(fMCEvent);
      TParticle * positiveMC = (TParticle*)TrueGammaCandidate0->GetPositiveMCDaughter(fMCEvent);
      TParticle * gammaMC0 = (TParticle*)fMCEvent->Particle(gamma0MCLabel);
      if(TMath::Abs(negativeMC->GetPdgCode())==11 && TMath::Abs(positiveMC->GetPdgCode())==11){  // Electrons ...
        if(negativeMC->GetUniqueID() == 5 && positiveMC->GetUniqueID() ==5){ // ... From Conversion ...
          if(gammaMC0->GetPdgCode() == 22){ // ... with Gamma Mother
            gamma0MotherLabel=gammaMC0->GetFirstMother();
            motherRealLabel=gammaMC0->GetFirstMother();
          }
        }
        if(gammaMC0->GetPdgCode() ==111){ // Dalitz candidate
          gamma0DalitzCand = kTRUE;
          gamma0MotherLabel=-111;
          motherRealLabel=gamma0MCLabel;
        }
      }
    }
    if(TrueGammaCandidate1->GetV0Index()<fInputEvent->GetNumberOfV0s()){
      Int_t gamma1MCLabel = TrueGammaCandidate1->GetMCParticleLabel(fMCEvent);
      Int_t gamma1MotherLabel = -1;
      if(gamma1MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
        // Daughters Gamma 1
        TParticle * negativeMC = (TParticle*)TrueGammaCandidate1->GetNegativeMCDaughter(fMCEvent);
        TParticle * positiveMC = (TParticle*)TrueGammaCandidate1->GetPositiveMCDaughter(fMCEvent);
        TParticle * gammaMC1 = (TParticle*)fMCEvent->Particle(gamma1MCLabel);
        if(TMath::Abs(negativeMC->GetPdgCode())==11 && TMath::Abs(positiveMC->GetPdgCode())==11){  // Electrons ...
          if(negativeMC->GetUniqueID() == 5 && positiveMC->GetUniqueID() ==5){ // ... From Conversion ...
            if(gammaMC1->GetPdgCode() == 22){ // ... with Gamma Mother
              gamma1MotherLabel=gammaMC1->GetFirstMother();
            }
          }
          if(gammaMC1->GetPdgCode() ==111 ){ // Dalitz candidate
            gamma1DalitzCand = kTRUE;
            gamma1MotherLabel=-111;
          }
        }
      }
      if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
        if(((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->GetPdgCode() == fPDGCodeNDM){
          isTrueNDM=kTRUE;
          if (CheckVectorForDoubleCount(fVectorDoubleCountTruePi0s,gamma0MotherLabel) && (!fDoLightOutput)) fHistoDoubleCountTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
        }
      }

      //Identify Dalitz candidate
      if (gamma1DalitzCand || gamma0DalitzCand){
        if (gamma0DalitzCand && gamma0MCLabel >=0 && gamma0MCLabel==gamma1MotherLabel){
          if (gamma0MotherLabel == -111) isTruePi0Dalitz = kTRUE;
        }
        if (gamma1DalitzCand && gamma1MCLabel >=0 && gamma1MCLabel==gamma0MotherLabel){
          if (gamma1MotherLabel == -111) isTruePi0Dalitz = kTRUE;
        }
      }


      if(isTrueNDM || isTruePi0Dalitz){// True Pion
        Pi0Candidate->SetTrueMesonValue(1);
        Pi0Candidate->SetMCLabel(motherRealLabel);
        if(!fDoLightOutput){
          fHistoTrueMotherGammaGammaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
          switch( fSelectedHeavyNeutralMeson ) {
          case 0: // ETA MESON
            if( IsEtaPiPlPiMiPiZeroDaughter(motherRealLabel) )
                fHistoTrueMotherGammaGammaFromHNMInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
            break;
          case 1: // OMEGA MESON
            if( IsOmegaPiPlPiMiPiZeroDaughter(motherRealLabel) )
                fHistoTrueMotherGammaGammaFromHNMInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
            break;
          case 2: // ETA PRIME MESON
            if( IsEtaPrimePiPlPiMiEtaDaughter(motherRealLabel) )
                fHistoTrueMotherGammaGammaFromHNMInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
            break;
          case 3: // D0 MESON
            if( IsD0PiPlPiMiPiZeroDaughter(motherRealLabel) )
                fHistoTrueMotherGammaGammaFromHNMInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
            break;
          default:
            AliError(Form("Heavy neutral meson not correctly selected (only 0,1,2,3 valid)... selected: %d",fSelectedHeavyNeutralMeson));
          }
        }
      }
    }
  }
}

//______________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::ProcessTrueNeutralPionCandidatesPureConversionsAOD(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1)
{

  // Process True Mesons
  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  Bool_t isTruePi0 = kFALSE;
  Bool_t isTruePi0Dalitz = kFALSE;
  Bool_t gamma0DalitzCand = kFALSE;
  Bool_t gamma1DalitzCand = kFALSE;
  Int_t motherRealLabel = -1;

  if (AODMCTrackArray!=nullptr && TrueGammaCandidate0 != nullptr){
    AliAODMCParticle *positiveMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelPositive()));
    AliAODMCParticle *negativeMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelNegative()));

    Int_t gamma0MCLabel = -1;
    Int_t gamma0MotherLabel = -1;
    if(!positiveMC||!negativeMC)
      return;

    if(positiveMC->GetMother()>-1&&(negativeMC->GetMother() == positiveMC->GetMother())){
      gamma0MCLabel = positiveMC->GetMother();
    }

    if(gamma0MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
      // Daughters Gamma 0
      AliAODMCParticle * gammaMC0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MCLabel));
      if(TMath::Abs(negativeMC->GetPdgCode())==11 && TMath::Abs(positiveMC->GetPdgCode())==11){  // Electrons ...
        if(((positiveMC->GetMCProcessCode())) == 5 && ((negativeMC->GetMCProcessCode())) == 5){ // ... From Conversion ...
          if(gammaMC0->GetPdgCode() == 22){ // ... with Gamma Mother
            gamma0MotherLabel=gammaMC0->GetMother();
            motherRealLabel=gammaMC0->GetMother();
          }
        }
        if(gammaMC0->GetPdgCode() ==111){ // Dalitz candidate
          gamma0DalitzCand = kTRUE;
          gamma0MotherLabel=-111;
          motherRealLabel=gamma0MCLabel;
        }
      }
    }
    positiveMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate1->GetMCLabelPositive()));
    negativeMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate1->GetMCLabelNegative()));

    Int_t gamma1MCLabel = -1;
    Int_t gamma1MotherLabel = -1;
    if(!positiveMC||!negativeMC)
      return;

    if(positiveMC->GetMother()>-1&&(negativeMC->GetMother() == positiveMC->GetMother())){
      gamma1MCLabel = positiveMC->GetMother();
    }
    if(gamma1MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
      // Daughters Gamma 1
      AliAODMCParticle * gammaMC1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MCLabel));
      if(TMath::Abs(negativeMC->GetPdgCode())==11 && TMath::Abs(positiveMC->GetPdgCode())==11){  // Electrons ...
        if(((positiveMC->GetMCProcessCode())) == 5 && ((negativeMC->GetMCProcessCode())) == 5){ // ... From Conversion ...
          if(gammaMC1->GetPdgCode() == 22){ // ... with Gamma Mother
          gamma1MotherLabel=gammaMC1->GetMother();
          }
        }
        if(gammaMC1->GetPdgCode() ==111 ){ // Dalitz candidate
            gamma1DalitzCand = kTRUE;
            gamma1MotherLabel=-111;
        }
      }
    }
    if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
      if(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == fPDGCodeNDM){
        isTruePi0=kTRUE;
        if (CheckVectorForDoubleCount(fVectorDoubleCountTruePi0s,gamma0MotherLabel) &&(!fDoLightOutput)) fHistoDoubleCountTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
      }
    }

    //Identify Dalitz candidate
    if (gamma1DalitzCand || gamma0DalitzCand){
      if (gamma0DalitzCand && gamma0MCLabel >=0 && gamma0MCLabel==gamma1MotherLabel){
        if (gamma0MotherLabel == -111) isTruePi0Dalitz = kTRUE;
      }
      if (gamma1DalitzCand && gamma1MCLabel >=0 && gamma1MCLabel==gamma0MotherLabel){
        if (gamma1MotherLabel == -111) isTruePi0Dalitz = kTRUE;
      }
    }

    if(isTruePi0 || isTruePi0Dalitz){// True Pion
      Pi0Candidate->SetTrueMesonValue(1);
      Pi0Candidate->SetMCLabel(motherRealLabel);
      if(!fDoLightOutput){
        fHistoTrueMotherGammaGammaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
        switch( fSelectedHeavyNeutralMeson ) {
        case 0: // ETA MESON
          if( IsEtaPiPlPiMiPiZeroDaughterAOD(AODMCTrackArray,motherRealLabel) )
              fHistoTrueMotherGammaGammaFromHNMInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
          break;
        case 1: // OMEGA MESON
          if( IsOmegaPiPlPiMiPiZeroDaughterAOD(AODMCTrackArray,motherRealLabel) )
              fHistoTrueMotherGammaGammaFromHNMInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
          break;
        case 2: // ETA PRIME MESON
          if( IsEtaPrimePiPlPiMiEtaDaughterAOD(AODMCTrackArray, motherRealLabel) )
              fHistoTrueMotherGammaGammaFromHNMInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
          break;
        case 3: // D0 MESON
          if( IsD0PiPlPiMiPiZeroDaughterAOD(AODMCTrackArray, motherRealLabel) )
              fHistoTrueMotherGammaGammaFromHNMInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
          break;
        default:
          AliError(Form("Heavy neutral meson not correctly selected (only 0,1,2,3 valid)... selected: %d",fSelectedHeavyNeutralMeson));
        }
      }
    }
  }
  return;
}


//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::ProcessNeutralPionCandidatesMixedConvCalo(){

  // Conversion Gammas
  if(fGoodConvGammas->GetEntries()>0){
    // vertex
    Double_t vertex[3] = {0};
    InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

    for(Int_t firstGammaIndex=0;firstGammaIndex<fGoodConvGammas->GetEntries();firstGammaIndex++){
      AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGoodConvGammas->At(firstGammaIndex));
      if (gamma0==nullptr) continue;

      for(Int_t secondGammaIndex=0;secondGammaIndex<fClusterCandidates->GetEntries();secondGammaIndex++){
        Bool_t matched = kFALSE;
        AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(secondGammaIndex));
        if (gamma1==nullptr) continue;

        if (gamma1->GetIsCaloPhoton() > 0){
          AliVCluster* cluster = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef());
          matched = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gamma0,cluster, fInputEvent );
        }

        AliAODConversionMother *NDMcand = new AliAODConversionMother(gamma0,gamma1);
        NDMcand->SetLabels(firstGammaIndex,secondGammaIndex);

        if(!fDoLightOutput){
          fHistoGammaGammaInvMassPtBeforeCuts[fiCut]->Fill(NDMcand->M(),NDMcand->Pt(), fWeightJetJetMC);
        }

        if( ((AliConversionMesonCuts*)fNeutralDecayMesonCutArray->At(fiCut))->DoGammaMinEnergyCut() ){
          Int_t minDaughters        = ((AliConversionMesonCuts*)fNeutralDecayMesonCutArray->At(fiCut))->GetNDaughterEnergyCut();
          Float_t minDaughterEnergy = ((AliConversionMesonCuts*)fNeutralDecayMesonCutArray->At(fiCut))->GetSingleDaughterMinE();
          if(minDaughters==1){ // at least one over threshold
             if( (gamma0->E() < minDaughterEnergy)  && (gamma1->E() < minDaughterEnergy)) continue;
          } else if (minDaughters==2){ // both over threshold
             if( (gamma0->E() < minDaughterEnergy)  || (gamma1->E() < minDaughterEnergy)) continue;
          }
        }

        if((((AliConversionMesonCuts*)fNeutralDecayMesonCutArray->At(fiCut))->MesonIsSelected(NDMcand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
          if (!matched){
            if(fIsMC){
              if(fInputEvent->IsA()==AliESDEvent::Class())
                ProcessTrueNeutralPionCandidatesMixedConvCalo(NDMcand,gamma0,gamma1);
              if(fInputEvent->IsA()==AliAODEvent::Class())
                ProcessTrueNeutralPionCandidatesMixedConvCaloAOD(NDMcand,gamma0,gamma1);
            }
            if (((AliConversionMesonCuts*)fNeutralDecayMesonCutArray->At(fiCut))->MesonIsSelectedByMassCut(NDMcand, 0)){
              fNeutralDecayParticleCandidates->Add(NDMcand);

              if(!fDoLightOutput){
                fHistoGammaGammaInvMassPt[fiCut]->Fill(NDMcand->M(),NDMcand->Pt(), fWeightJetJetMC);
              }
            } else if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseSidebandMixing()) &&
                      (((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedByMassCut(NDMcand, 1))){
              fNeutralDecayParticleSidebandCandidates->Add(NDMcand);

              if(!fDoLightOutput){
                fHistoGammaGammaInvMassPt[fiCut]->Fill(NDMcand->M(),NDMcand->Pt(), fWeightJetJetMC);
              }
            } else if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseSidebandMixingBothSides()) &&
                      ((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedByMassCut(NDMcand, 2)) ||
                      ((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedByMassCut(NDMcand, 3))))){
              fNeutralDecayParticleSidebandCandidates->Add(NDMcand);

              if(!fDoLightOutput){
                fHistoGammaGammaInvMassPt[fiCut]->Fill(NDMcand->M(),NDMcand->Pt(), fWeightJetJetMC);
              }
            } else{
              delete NDMcand;
              NDMcand=0x0;
            }
          }else{
            delete NDMcand;
            NDMcand=0x0;
          }
        }else{
          delete NDMcand;
          NDMcand=0x0;
        }
      }
    }
  }
}

//______________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::ProcessTrueNeutralPionCandidatesMixedConvCalo( AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1)
{
  // Process True Mesons
  if(TrueGammaCandidate0->GetV0Index()<fInputEvent->GetNumberOfV0s()){
    Bool_t isTruePi0 = kFALSE;
    Bool_t isTruePi0Dalitz = kFALSE;
    Bool_t gamma0DalitzCand = kFALSE;

    Int_t gamma0MCLabel = TrueGammaCandidate0->GetMCParticleLabel(fMCEvent);
    Int_t gamma0MotherLabel = -1;
    Int_t motherRealLabel = -1;
    if(gamma0MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
      // Daughters Gamma 0
      TParticle * negativeMC = (TParticle*)TrueGammaCandidate0->GetNegativeMCDaughter(fMCEvent);
      TParticle * positiveMC = (TParticle*)TrueGammaCandidate0->GetPositiveMCDaughter(fMCEvent);
      TParticle * gammaMC0 = (TParticle*)fMCEvent->Particle(gamma0MCLabel);
      if(TMath::Abs(negativeMC->GetPdgCode())==11 && TMath::Abs(positiveMC->GetPdgCode())==11){  // Electrons ...
        if(negativeMC->GetUniqueID() == 5 && positiveMC->GetUniqueID() ==5){ // ... From Conversion ...
          if(gammaMC0->GetPdgCode() == 22){ // ... with Gamma Mother
            gamma0MotherLabel=gammaMC0->GetFirstMother();
            motherRealLabel=gammaMC0->GetFirstMother();
          }
        }
        if(gammaMC0->GetPdgCode() ==111){ // Dalitz candidate
          gamma0DalitzCand = kTRUE;
          gamma0MotherLabel=-111;
          motherRealLabel=gamma0MCLabel;
        }

      }
    }

    if (TrueGammaCandidate1->GetIsCaloPhoton() == 0) AliFatal("CaloPhotonFlag has not been set. Aborting");

    Int_t gamma1MCLabel = TrueGammaCandidate1->GetCaloPhotonMCLabel(0); 	// get most probable MC label
    Int_t gamma1MotherLabel = -1;
    // check if

    if(gamma1MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
      // Daughters Gamma 1
      TParticle * gammaMC1 = (TParticle*)fMCEvent->Particle(gamma1MCLabel);
      if (TrueGammaCandidate1->IsLargestComponentPhoton() || TrueGammaCandidate1->IsLargestComponentElectron()){		// largest component is electro magnetic
        // get mother of interest (pi0 or eta)
        if (TrueGammaCandidate1->IsLargestComponentPhoton()){														// for photons its the direct mother
          gamma1MotherLabel=gammaMC1->GetMother(0);
        } else if (TrueGammaCandidate1->IsLargestComponentElectron()){ 												// for electrons its either the direct mother or for conversions the grandmother
                    if (TrueGammaCandidate1->IsConversion() && gammaMC1->GetMother(0)>-1) gamma1MotherLabel=fMCEvent->Particle(gammaMC1->GetMother(0))->GetMother(0);
          else gamma1MotherLabel=gammaMC1->GetMother(0);
        }
      }
    }

    if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
      if(((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->GetPdgCode() == fPDGCodeNDM){
        isTruePi0=kTRUE;
        if (CheckVectorForDoubleCount(fVectorDoubleCountTruePi0s,gamma0MotherLabel) && (!fDoLightOutput)) fHistoDoubleCountTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
      }
    }

    if (gamma0DalitzCand ){
      if (gamma0DalitzCand && gamma0MCLabel >=0 && gamma0MCLabel==gamma1MotherLabel){
        if (gamma0MotherLabel == -111) isTruePi0Dalitz = kTRUE;
      }
    }

    if(isTruePi0 || isTruePi0Dalitz ){
      Pi0Candidate->SetTrueMesonValue(1);
      Pi0Candidate->SetMCLabel(motherRealLabel);
      if(!fDoLightOutput){
        fHistoTrueMotherGammaGammaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
        switch( fSelectedHeavyNeutralMeson ) {
        case 0: // ETA MESON
          if( IsEtaPiPlPiMiPiZeroDaughter(motherRealLabel) )
              fHistoTrueMotherGammaGammaFromHNMInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
          break;
        case 1: // OMEGA MESON
          if( IsOmegaPiPlPiMiPiZeroDaughter(motherRealLabel) )
              fHistoTrueMotherGammaGammaFromHNMInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
          break;
        case 2: // ETA PRIME MESON
          if( IsEtaPrimePiPlPiMiEtaDaughter(motherRealLabel) )
              fHistoTrueMotherGammaGammaFromHNMInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
          break;
        case 3: // D0 MESON
          if( IsD0PiPlPiMiPiZeroDaughter(motherRealLabel) )
              fHistoTrueMotherGammaGammaFromHNMInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
          break;
        default:
          AliError(Form("Heavy neutral meson not correctly selected (only 0,1,2,3 valid)... selected: %d",fSelectedHeavyNeutralMeson));
        }
      }
    }
  }
}

//______________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::ProcessTrueNeutralPionCandidatesMixedConvCaloAOD( AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1)
{
  // Process True Mesons
  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (AODMCTrackArray == NULL) return;

  AliAODMCParticle * negativeMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelNegative()));
  AliAODMCParticle * positiveMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelPositive()));
  if(TrueGammaCandidate0->GetV0Index()<fInputEvent->GetNumberOfV0s()){
    Bool_t isTruePi0 = kFALSE;
    Bool_t isTruePi0Dalitz = kFALSE;
    Bool_t gamma0DalitzCand = kFALSE;

    Int_t gamma0MCLabel = positiveMC->GetMother();// check that this always works
    Int_t gamma0MotherLabel = -1;

    Int_t motherRealLabel = -1;
    if(gamma0MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
      // Daughters Gamma 0
      AliAODMCParticle * gammaMC0   = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MCLabel));
      if(TMath::Abs(negativeMC->GetPdgCode())==11 && TMath::Abs(positiveMC->GetPdgCode())==11){  // Electrons ...
        if(negativeMC->GetMCProcessCode() == 5 && positiveMC->GetMCProcessCode() ==5){ // ... From Conversion ...
          if(gammaMC0->GetPdgCode() == 22){ // ... with Gamma Mother
            gamma0MotherLabel=gammaMC0->GetMother();
            motherRealLabel=gammaMC0->GetMother();
          }
        }
        if(gammaMC0->GetPdgCode() ==111){ // Dalitz candidate
          gamma0DalitzCand = kTRUE;
          gamma0MotherLabel=-111;
          motherRealLabel=gamma0MCLabel;
        }

      }
    }

    if (TrueGammaCandidate1->GetIsCaloPhoton() == 0) AliFatal("CaloPhotonFlag has not been set. Aborting");

    Int_t gamma1MCLabel = TrueGammaCandidate1->GetCaloPhotonMCLabel(0); 	// get most probable MC label
    Int_t gamma1MotherLabel = -1;
    // check if

    if(gamma1MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
      // Daughters Gamma 1
      AliAODMCParticle *gammaMC1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MCLabel));
      if (TrueGammaCandidate1->IsLargestComponentPhoton() || TrueGammaCandidate1->IsLargestComponentElectron()){		// largest component is electro magnetic
        // get mother of interest (pi0 or eta)
        if (TrueGammaCandidate1->IsLargestComponentPhoton()){														// for photons its the direct mother
          gamma1MotherLabel=gammaMC1->GetMother();
        } else if (TrueGammaCandidate1->IsLargestComponentElectron()){ 												// for electrons its either the direct mother or for conversions the grandmother
          if (TrueGammaCandidate1->IsConversion() && gammaMC1->GetMother()>-1) gamma1MotherLabel= (static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gammaMC1->GetMother())))->GetMother();
          else gamma1MotherLabel=gammaMC1->GetMother();
        }
      }
    }

    if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
      if((static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel)))->GetPdgCode() == fPDGCodeNDM){
        isTruePi0=kTRUE;
        if (CheckVectorForDoubleCount(fVectorDoubleCountTruePi0s,gamma0MotherLabel) && (!fDoLightOutput)) fHistoDoubleCountTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
      }
    }

    if (gamma0DalitzCand ){
      if (gamma0DalitzCand && gamma0MCLabel >=0 && gamma0MCLabel==gamma1MotherLabel){
        if (gamma0MotherLabel == -111) isTruePi0Dalitz = kTRUE;
      }
    }

    if(isTruePi0 || isTruePi0Dalitz ){
      Pi0Candidate->SetTrueMesonValue(1);
      Pi0Candidate->SetMCLabel(motherRealLabel);
      if(!fDoLightOutput){
        fHistoTrueMotherGammaGammaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
        switch( fSelectedHeavyNeutralMeson ) {
        case 0: // ETA MESON
          if( IsEtaPiPlPiMiPiZeroDaughterAOD(AODMCTrackArray,motherRealLabel) )
              fHistoTrueMotherGammaGammaFromHNMInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
          break;
        case 1: // OMEGA MESON
          if( IsOmegaPiPlPiMiPiZeroDaughterAOD(AODMCTrackArray, motherRealLabel) )
              fHistoTrueMotherGammaGammaFromHNMInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
          break;
        case 2: // ETA PRIME MESON
          if( IsEtaPrimePiPlPiMiEtaDaughterAOD(AODMCTrackArray, motherRealLabel) )
              fHistoTrueMotherGammaGammaFromHNMInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
          break;
        case 3: // D0 MESON
          if( IsD0PiPlPiMiPiZeroDaughterAOD(AODMCTrackArray, motherRealLabel) )
              fHistoTrueMotherGammaGammaFromHNMInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
          break;
        default:
          AliError(Form("Heavy neutral meson not correctly selected (only 0,1,2,3 valid)... selected: %d",fSelectedHeavyNeutralMeson));
        }
      }
    }
  }
}


//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::ProcessPionCandidates(){

  Bool_t use4vecformass = ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->Use4VecForMass();

  Double_t magField = fInputEvent->GetMagneticField();
  if( magField  < 0.0 ){
    magField =  1.0;
  } else {
    magField =  -1.0;
  }

  vector<Int_t> lGoodNegPionIndexPrev(0);
  vector<Int_t> lGoodPosPionIndexPrev(0);
  for(UInt_t i = 0; i < fSelectorNegPionIndex.size(); i++){

    AliESDtrack* negPionCandidate =dynamic_cast<AliESDtrack*>(fInputEvent->GetTrack(fSelectorNegPionIndex[i]));
    if(! ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelected(negPionCandidate) ) continue;
    lGoodNegPionIndexPrev.push_back(   fSelectorNegPionIndex[i] );

    TLorentzVector* negPionforHandler = new TLorentzVector();
    negPionforHandler->SetPxPyPzE(negPionCandidate->Px(), negPionCandidate->Py(), negPionCandidate->Pz(), negPionCandidate->E());
    FixPzVecToMatchPDGInvMass(negPionforHandler);
    AliAODConversionPhoton *negPionHandler = new AliAODConversionPhoton(negPionforHandler);
    delete negPionforHandler;

    fNegPionCandidates->Add(negPionHandler);
    if(!fDoLightOutput){
        fHistoNegPionPt[fiCut]->Fill(negPionCandidate->Pt(), fWeightJetJetMC);
        fHistoNegPionPhi[fiCut]->Fill(negPionCandidate->Phi(), fWeightJetJetMC);
    }

    if( fMCEvent ) {
      const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
      Double_t mcProdVtxX 	= primVtxMC->GetX();
      Double_t mcProdVtxY 	= primVtxMC->GetY();
      Double_t mcProdVtxZ 	= primVtxMC->GetZ();

      Int_t labelNegPion = TMath::Abs( negPionCandidate->GetLabel() );
      Bool_t negPionIsPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, labelNegPion, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
            if( labelNegPion>-1 && labelNegPion < fMCEvent->GetNumberOfTracks() ){
        TParticle* negPion = fMCEvent->Particle(labelNegPion);
        if( negPion->GetPdgCode() ==  -211 ){
          if(!fDoLightOutput){
            if( negPionIsPrimary ){
                fHistoTrueNegPionPt[fiCut]->Fill(negPionCandidate->Pt(), fWeightJetJetMC);    //primary negPion
            }
            switch( fSelectedHeavyNeutralMeson ) {
            case 0: // ETA MESON
              if( IsEtaPiPlPiMiPiZeroDaughter(labelNegPion) && negPionIsPrimary )
                fHistoTrueNegPionFromNeutralMesonPt[fiCut]->Fill(negPionCandidate->Pt(), fWeightJetJetMC);
              break;
            case 1: // OMEGA MESON
              if( IsOmegaPiPlPiMiPiZeroDaughter(labelNegPion) && negPionIsPrimary)
                fHistoTrueNegPionFromNeutralMesonPt[fiCut]->Fill(negPionCandidate->Pt(), fWeightJetJetMC);
              break;
            case 2: // ETA PRIME MESON
              if( IsEtaPrimePiPlPiMiEtaDaughter(labelNegPion) &&  negPionIsPrimary)
                fHistoTrueNegPionFromNeutralMesonPt[fiCut]->Fill(negPionCandidate->Pt(), fWeightJetJetMC);
              break;
            case 3: // D0 MESON
              if( IsD0PiPlPiMiPiZeroDaughter(labelNegPion)& negPionIsPrimary)
                fHistoTrueNegPionFromNeutralMesonPt[fiCut]->Fill(negPionCandidate->Pt(), fWeightJetJetMC);
              break;
            default:
              AliError(Form("Heavy neutral meson not correctly selected (only 0,1,2,3 valid)... selected: %d",fSelectedHeavyNeutralMeson));
            }
          }
        }
      }
    }
  }
  for(UInt_t i = 0; i < fSelectorPosPionIndex.size(); i++){
    AliESDtrack* posPionCandidate = dynamic_cast<AliESDtrack*>(fInputEvent->GetTrack(fSelectorPosPionIndex[i]));
    if(! ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelected(posPionCandidate) ) continue;
    lGoodPosPionIndexPrev.push_back(   fSelectorPosPionIndex[i]  );

    TLorentzVector* posPionforHandler = new TLorentzVector();
    posPionforHandler->SetPxPyPzE(posPionCandidate->Px(), posPionCandidate->Py(), posPionCandidate->Pz(), posPionCandidate->E());
    FixPzVecToMatchPDGInvMass(posPionforHandler);
    AliAODConversionPhoton *posPionHandler = new AliAODConversionPhoton(posPionforHandler);
    delete posPionforHandler;

    fPosPionCandidates->Add(posPionHandler);
    if(!fDoLightOutput){
        fHistoPosPionPt[fiCut]->Fill(posPionCandidate->Pt(), fWeightJetJetMC);
        fHistoPosPionPhi[fiCut]->Fill(posPionCandidate->Phi(), fWeightJetJetMC);
    }
    if( fMCEvent ) {
      const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
      Double_t mcProdVtxX 	= primVtxMC->GetX();
      Double_t mcProdVtxY 	= primVtxMC->GetY();
      Double_t mcProdVtxZ 	= primVtxMC->GetZ();

      Int_t labelPosPion = TMath::Abs( posPionCandidate->GetLabel() );
      Bool_t posPionIsPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, labelPosPion, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
            if( labelPosPion>-1 && labelPosPion < fMCEvent->GetNumberOfTracks() ) {
        TParticle* posPion = fMCEvent->Particle(labelPosPion);
        if( posPion->GetPdgCode() ==  211 ){
          if(!fDoLightOutput){
            if( posPionIsPrimary ){
              fHistoTruePosPionPt[fiCut]->Fill(posPionCandidate->Pt(), fWeightJetJetMC);
            }
            switch( fSelectedHeavyNeutralMeson ) {
            case 0: // ETA MESON
              if( IsEtaPiPlPiMiPiZeroDaughter(labelPosPion) && posPionIsPrimary )
                fHistoTruePosPionFromNeutralMesonPt[fiCut]->Fill(posPionCandidate->Pt(), fWeightJetJetMC);
              break;
            case 1: // OMEGA MESON
              if( IsOmegaPiPlPiMiPiZeroDaughter(labelPosPion) && posPionIsPrimary)
                fHistoTruePosPionFromNeutralMesonPt[fiCut]->Fill(posPionCandidate->Pt(), fWeightJetJetMC);
              break;
            case 2: // ETA PRIME MESON
              if( IsEtaPrimePiPlPiMiEtaDaughter(labelPosPion) &&  posPionIsPrimary)
                fHistoTruePosPionFromNeutralMesonPt[fiCut]->Fill(posPionCandidate->Pt(), fWeightJetJetMC);
              break;
            case 3: // D0 MESON
              if( IsD0PiPlPiMiPiZeroDaughter(labelPosPion) && posPionIsPrimary )
                fHistoTruePosPionFromNeutralMesonPt[fiCut]->Fill(posPionCandidate->Pt(), fWeightJetJetMC);
              break;
            default:
              AliError(Form("Heavy neutral meson not correctly selected (only 0,1,2,3 valid)... selected: %d",fSelectedHeavyNeutralMeson));
            }
          }
        }
      }
    }
  }
 // AliInfo(Form("Number of good neg pions: %i \t pos pions = %i", lGoodNegPionIndexPrev.size(), lGoodPosPionIndexPrev.size()));

  for(UInt_t i = 0; i < lGoodNegPionIndexPrev.size(); i++){
    AliVTrack* negPionCandidate = dynamic_cast<AliVTrack*>(fInputEvent->GetTrack(lGoodNegPionIndexPrev[i]));

    //AliKFParticle negPionCandidateKF( *negPionCandidate->GetConstrainedParam(), 211 );
    AliKFParticle negPionCandidateKF( *negPionCandidate, 211 );
    for(UInt_t j = 0; j < lGoodPosPionIndexPrev.size(); j++){
      AliVTrack *posPionCandidate = dynamic_cast<AliVTrack*>(fInputEvent->GetTrack(lGoodPosPionIndexPrev[j]));
      //AliKFParticle posPionCandidateKF( *posPionCandidate->GetConstrainedParam(), 211 );
      AliKFParticle posPionCandidateKF( *posPionCandidate, 211 );
      // AliInfo(Form("virtualPhoton distance before pi+ pi- = %f", negPionCandidateKF.GetDistanceFromParticle(posPionCandidateKF)));

      AliKFConversionPhoton virtualPhoton(negPionCandidateKF,posPionCandidateKF);

     // AliInfo(Form(" GetPrimaryVertex() x=%f y=%f z=%f",primx,primy,primz));
     // virtualPhoton->SetProductionVertex(primaryVertex);
      virtualPhoton.SetTrackLabels( lGoodPosPionIndexPrev[j], lGoodNegPionIndexPrev[i]);

      TLorentzVector posPionVec4;
      TLorentzVector negPionVec4;
      TLorentzVector virtPionVec4;
      TLorentzVector posKFPionVec4;
      TLorentzVector negKFPionVec4;
      TLorentzVector virtKFPionVec4;

      if(use4vecformass){
        posKFPionVec4.SetPxPyPzE(posPionCandidateKF.Px(),posPionCandidateKF.Py(),posPionCandidateKF.Pz(),posPionCandidateKF.E());
        negKFPionVec4.SetPxPyPzE(negPionCandidateKF.Px(),negPionCandidateKF.Py(),negPionCandidateKF.Pz(),negPionCandidateKF.E());
        virtKFPionVec4 = posKFPionVec4 + negKFPionVec4;
      }

      Int_t labeln=0;
      Int_t labelp=0;
      Int_t motherlabelp = 0;
      Int_t motherlabeln = 0;
      TParticle *fNegativeMCParticle =nullptr;
      TParticle *fPositiveMCParticle =nullptr;
      if( fMCEvent ) {
        labeln=TMath::Abs(negPionCandidate->GetLabel());
        labelp=TMath::Abs(posPionCandidate->GetLabel());
        if(labeln>-1) fNegativeMCParticle = fMCEvent->Particle(labeln);
        if(labelp>-1) fPositiveMCParticle = fMCEvent->Particle(labelp);
        // check whether MC particles exist, else abort
        if (fNegativeMCParticle == nullptr || fPositiveMCParticle == nullptr) return;

        motherlabeln = fNegativeMCParticle->GetMother(0);
        motherlabelp = fPositiveMCParticle->GetMother(0);
        virtualPhoton.SetMCLabelPositive(labelp);
        virtualPhoton.SetMCLabelNegative(labeln);

      }
      AliAODConversionPhoton *vParticle = new AliAODConversionPhoton(&virtualPhoton); //To apply mass 2 pion mass cut
      if(use4vecformass){
        vParticle->SetPxPyPzE(virtKFPionVec4.Px(),virtKFPionVec4.Py(),virtKFPionVec4.Pz(),virtKFPionVec4.E());
        // set mass to one calculated from four vector, not from parameters
        vParticle->SetMass(vParticle->M());
      }
      if(!fDoLightOutput){
        Double_t ds,dsp;
        posPionCandidateKF.GetDStoParticle(negPionCandidateKF,ds,dsp);
        Float_t chi2 = virtualPhoton.GetChi2perNDF();
        if(chi2>299) chi2 = 299; // to illustrate overflow bin
        fHistovParticleChi2PerNDF[fiCut]->Fill(chi2,fWeightJetJetMC);
        fHistovParticledS[fiCut]->Fill(ds,fWeightJetJetMC);
        if (fMCEvent &&(fDoMesonQA>0)){
          if (fPositiveMCParticle && fNegativeMCParticle ) {
            if (((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->DoMassCut()){
              if (vParticle->GetMass() < ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetMassCut()){
                if(TMath::Abs(fNegativeMCParticle->GetPdgCode())==211 && TMath::Abs(fPositiveMCParticle->GetPdgCode())==211){  // Pions ...
                  fHistoTruePionPionInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt(), fWeightJetJetMC);
                  fHistoTruevParticleChi2PerNDF[fiCut]->Fill(chi2,fWeightJetJetMC);
                  fHistoTruevParticledS[fiCut]->Fill(ds,fWeightJetJetMC);
                  if (motherlabeln == motherlabelp){
                    fHistoTruePionPionFromSameMotherInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt(), fWeightJetJetMC);
                    fHistoTruevParticleFromSameMotherChi2PerNDF[fiCut]->Fill(chi2,fWeightJetJetMC);
                    fHistoTruevParticleFromSameMotherdS[fiCut]->Fill(ds,fWeightJetJetMC);
                    switch( fSelectedHeavyNeutralMeson ) {
                    case 0: // ETA MESON
                      if( IsEtaPiPlPiMiPiZeroDaughter(labeln) ){
                        fHistoTruePionPionFromHNMInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt(), fWeightJetJetMC);
                        fHistoTruevParticleFromHNMChi2PerNDF[fiCut]->Fill(chi2,fWeightJetJetMC);
                        fHistoTruevParticleFromHNMdS[fiCut]->Fill(ds,fWeightJetJetMC);
                      }
                      break;
                    case 1: // OMEGA MESON
                      if( IsOmegaPiPlPiMiPiZeroDaughter(labeln) ){
                        fHistoTruePionPionFromHNMInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt(), fWeightJetJetMC);
                        fHistoTruevParticleFromHNMChi2PerNDF[fiCut]->Fill(chi2,fWeightJetJetMC);
                        fHistoTruevParticleFromHNMdS[fiCut]->Fill(ds,fWeightJetJetMC);
                      }
                      break;
                    case 2: // ETA PRIME MESON
                      if( IsEtaPrimePiPlPiMiEtaDaughter(labeln) ){
                        fHistoTruePionPionFromHNMInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt(), fWeightJetJetMC);
                        fHistoTruevParticleFromHNMChi2PerNDF[fiCut]->Fill(chi2,fWeightJetJetMC);
                        fHistoTruevParticleFromHNMdS[fiCut]->Fill(ds,fWeightJetJetMC);
                      }
                      break;
                    case 3: // D0 MESON
                      if( IsD0PiPlPiMiPiZeroDaughter(labeln) ){
                        fHistoTruePionPionFromHNMInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt(), fWeightJetJetMC);
                        fHistoTruevParticleFromHNMChi2PerNDF[fiCut]->Fill(chi2,fWeightJetJetMC);
                        fHistoTruevParticleFromHNMdS[fiCut]->Fill(ds,fWeightJetJetMC);
                      }
                      break;
                    default:
                      AliError(Form("Heavy neutral meson not correctly selected (only 0,1,2,3 valid)... selected: %d",fSelectedHeavyNeutralMeson));
                    }
                  }
                }
              }
            } else {
              if(TMath::Abs(fNegativeMCParticle->GetPdgCode())==211 && TMath::Abs(fPositiveMCParticle->GetPdgCode())==211){  // Pions ...
                fHistoTruePionPionInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt(), fWeightJetJetMC);
                fHistoTruevParticleChi2PerNDF[fiCut]->Fill(chi2,fWeightJetJetMC);
                fHistoTruevParticledS[fiCut]->Fill(ds,fWeightJetJetMC);
                if (motherlabeln == motherlabelp){
                  fHistoTruePionPionFromSameMotherInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt(), fWeightJetJetMC);
                  fHistoTruevParticleFromSameMotherChi2PerNDF[fiCut]->Fill(chi2,fWeightJetJetMC);
                  fHistoTruevParticleFromSameMotherdS[fiCut]->Fill(ds,fWeightJetJetMC);
                  switch( fSelectedHeavyNeutralMeson ) {
                  case 0: // ETA MESON
                    if( IsEtaPiPlPiMiPiZeroDaughter(labeln) ){
                      fHistoTruePionPionFromHNMInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt(), fWeightJetJetMC);
                      fHistoTruevParticleFromHNMChi2PerNDF[fiCut]->Fill(chi2,fWeightJetJetMC);
                      fHistoTruevParticleFromHNMdS[fiCut]->Fill(ds,fWeightJetJetMC);
                    }
                    break;
                  case 1: // OMEGA MESON
                    if( IsOmegaPiPlPiMiPiZeroDaughter(labeln) ){
                      fHistoTruePionPionFromHNMInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt(), fWeightJetJetMC);
                      fHistoTruevParticleFromHNMChi2PerNDF[fiCut]->Fill(chi2,fWeightJetJetMC);
                      fHistoTruevParticleFromHNMdS[fiCut]->Fill(ds,fWeightJetJetMC);
                    }
                    break;
                  case 2: // ETA PRIME MESON
                    if( IsEtaPrimePiPlPiMiEtaDaughter(labeln) ){
                      fHistoTruePionPionFromHNMInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt(), fWeightJetJetMC);
                      fHistoTruevParticleFromHNMChi2PerNDF[fiCut]->Fill(chi2,fWeightJetJetMC);
                      fHistoTruevParticleFromHNMdS[fiCut]->Fill(ds,fWeightJetJetMC);
                    }
                    break;
                  case 3: // D0 MESON
                    if( IsD0PiPlPiMiPiZeroDaughter(labeln) ){
                      fHistoTruePionPionFromHNMInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt(), fWeightJetJetMC);
                      fHistoTruevParticleFromHNMChi2PerNDF[fiCut]->Fill(chi2,fWeightJetJetMC);
                      fHistoTruevParticleFromHNMdS[fiCut]->Fill(ds,fWeightJetJetMC);
                    }
                    break;
                  default:
                    AliError(Form("Heavy neutral meson not correctly selected (only 0,1,2,3 valid)... selected: %d",fSelectedHeavyNeutralMeson));
                  }
                }
              }
            }
          }
        }
      }

      Bool_t survivesMassCut = kFALSE;

      if (((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->DoMassCut()){
        if (vParticle->GetMass() < ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetMassCut()){
          survivesMassCut = kTRUE;
        }
      } else{
        survivesMassCut = kTRUE;
      }


      if(survivesMassCut){
          //fGoodVirtualParticles->Add( vParticle );
          if(!fDoLightOutput){
            fHistoPionPionInvMassPt[fiCut]->Fill( vParticle->GetMass(),vParticle->Pt(), fWeightJetJetMC);
          }
          CalculateMesonCandidates(vParticle);
          delete vParticle;
          vParticle=0x0;
      }else{
          delete vParticle;
          vParticle=0x0;
      }
    }
  }

  Double_t clsToFPos = -1.0;
  Double_t clsToFNeg = -1.0;

  Float_t dcaToVertexXYPos = -1.0;
  Float_t dcaToVertexZPos  = -1.0;
  Float_t dcaToVertexXYNeg = -1.0;
  Float_t dcaToVertexZNeg  = -1.0;

  if ( fDoMesonQA>0 ) {
    for(UInt_t i = 0; i < lGoodNegPionIndexPrev.size(); i++){
      AliVTrack* negPionCandidate = dynamic_cast<AliVTrack*>(fInputEvent->GetTrack(lGoodNegPionIndexPrev[i]));

      clsToFNeg = ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetNFindableClustersTPC(negPionCandidate);

      Float_t bNeg[2];
      Float_t bCovNeg[3];
      negPionCandidate->GetImpactParameters(bNeg,bCovNeg);
      if (bCovNeg[0]<=0 || bCovNeg[2]<=0) {
        AliDebug(1, "Estimated b resolution lower or equal zero!");
        bCovNeg[0]=0; bCovNeg[2]=0;
      }

      dcaToVertexXYNeg = bNeg[0];
      dcaToVertexZNeg  = bNeg[1];

      if(!fDoLightOutput){
        fHistoNegPionEta[fiCut]->Fill(negPionCandidate->Eta(), fWeightJetJetMC);
        fHistoNegPionClsTPC[fiCut]->Fill(clsToFNeg,negPionCandidate->Pt(), fWeightJetJetMC);

        fHistoPionDCAxy[fiCut]->Fill(  dcaToVertexXYNeg, negPionCandidate->Pt(), fWeightJetJetMC );
        fHistoPionDCAz[fiCut]->Fill(   dcaToVertexZNeg,  negPionCandidate->Pt(), fWeightJetJetMC );

        fHistoPionTPCdEdxNSigma[fiCut]->Fill( negPionCandidate->P(),((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetPIDResponse()->NumberOfSigmasTPC(negPionCandidate, AliPID::kPion), fWeightJetJetMC );

        fHistoPionTPCdEdx[fiCut]->Fill(negPionCandidate->P(), TMath::Abs(negPionCandidate->GetTPCsignal()), fWeightJetJetMC);
      }
    }

    for(UInt_t i = 0; i < lGoodPosPionIndexPrev.size(); i++){
      AliVTrack* posPionCandidate = dynamic_cast<AliVTrack*>(fInputEvent->GetTrack(lGoodPosPionIndexPrev[i]));

      clsToFPos = ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetNFindableClustersTPC(posPionCandidate);

      Float_t bPos[2];
      Float_t bCovPos[3];
      posPionCandidate->GetImpactParameters(bPos,bCovPos);
      if (bCovPos[0]<=0 || bCovPos[2]<=0) {
        AliDebug(1, "Estimated b resolution lower or equal zero!");
        bCovPos[0]=0; bCovPos[2]=0;
      }

      dcaToVertexXYPos = bPos[0];
      dcaToVertexZPos  = bPos[1];

      if(!fDoLightOutput){
        fHistoPosPionEta[fiCut]->Fill(posPionCandidate->Eta(), fWeightJetJetMC);
        fHistoPosPionClsTPC[fiCut]->Fill(clsToFPos,posPionCandidate->Pt(), fWeightJetJetMC);

        fHistoPionDCAxy[fiCut]->Fill(  dcaToVertexXYPos, posPionCandidate->Pt(), fWeightJetJetMC );
        fHistoPionDCAz[fiCut]->Fill(   dcaToVertexZPos,  posPionCandidate->Pt(), fWeightJetJetMC );

        fHistoPionTPCdEdxNSigma[fiCut]->Fill( posPionCandidate->P(),((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetPIDResponse()->NumberOfSigmasTPC(posPionCandidate, AliPID::kPion), fWeightJetJetMC );

        fHistoPionTPCdEdx[fiCut]->Fill(posPionCandidate->P(), TMath::Abs(posPionCandidate->GetTPCsignal()), fWeightJetJetMC);
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::ProcessPionCandidatesAOD(){

  Bool_t use4vecformass = ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->Use4VecForMass();

  Double_t magField = fInputEvent->GetMagneticField();
  if( magField  < 0.0 ){
    magField =  1.0;
  } else {
    magField =  -1.0;
  }

  vector<Int_t> lGoodNegPionIndexPrev(0);
  vector<Int_t> lGoodPosPionIndexPrev(0);

  TClonesArray *AODMCTrackArray = NULL;
  if(fMCEvent){
     AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  }

  for(UInt_t i = 0; i < fSelectorNegPionIndex.size(); i++){

    AliAODTrack* negPionCandidate =dynamic_cast<AliAODTrack*>(fInputEvent->GetTrack(fSelectorNegPionIndex[i]));
    if(! ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedAOD(negPionCandidate) ) continue;
    lGoodNegPionIndexPrev.push_back(   fSelectorNegPionIndex[i] );

    TLorentzVector* negPionforHandler = new TLorentzVector();
    negPionforHandler->SetPxPyPzE(negPionCandidate->Px(), negPionCandidate->Py(), negPionCandidate->Pz(), negPionCandidate->E());
    FixPzVecToMatchPDGInvMass(negPionforHandler);
    AliAODConversionPhoton *negPionHandler = new AliAODConversionPhoton(negPionforHandler);
    delete negPionforHandler;
    fNegPionCandidates->Add(negPionHandler);
    if(!fDoLightOutput){
        fHistoNegPionPt[fiCut]->Fill(negPionCandidate->Pt(), fWeightJetJetMC);
        fHistoNegPionPhi[fiCut]->Fill(negPionCandidate->Phi(), fWeightJetJetMC);
    }

    if( fMCEvent ) {
      const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
      Double_t mcProdVtxX 	= primVtxMC->GetX();
      Double_t mcProdVtxY 	= primVtxMC->GetY();
      Double_t mcProdVtxZ 	= primVtxMC->GetZ();

      Int_t labelNegPion = TMath::Abs( negPionCandidate->GetLabel() );
        if( labelNegPion>-1 && labelNegPion < AODMCTrackArray->GetEntriesFast() ){
        AliAODMCParticle* negPion =  static_cast<AliAODMCParticle*>(AODMCTrackArray->At(labelNegPion));
        Bool_t negPionIsPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD( fInputEvent, negPion, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
        if( negPion->GetPdgCode() ==  -211 ){
          if(!fDoLightOutput){
            if( negPionIsPrimary ){
                fHistoTrueNegPionPt[fiCut]->Fill(negPionCandidate->Pt(), fWeightJetJetMC);    //primary negPion
            }
            switch( fSelectedHeavyNeutralMeson ) {
            case 0: // ETA MESON
              if( IsEtaPiPlPiMiPiZeroDaughterAOD(AODMCTrackArray, labelNegPion) && negPionIsPrimary )
                fHistoTrueNegPionFromNeutralMesonPt[fiCut]->Fill(negPionCandidate->Pt(), fWeightJetJetMC);
              break;
            case 1: // OMEGA MESON
              if( IsOmegaPiPlPiMiPiZeroDaughterAOD(AODMCTrackArray, labelNegPion) && negPionIsPrimary)
                fHistoTrueNegPionFromNeutralMesonPt[fiCut]->Fill(negPionCandidate->Pt(), fWeightJetJetMC);
              break;
            case 2: // ETA PRIME MESON
              if( IsEtaPrimePiPlPiMiEtaDaughterAOD(AODMCTrackArray,labelNegPion) &&  negPionIsPrimary)
                fHistoTrueNegPionFromNeutralMesonPt[fiCut]->Fill(negPionCandidate->Pt(), fWeightJetJetMC);
              break;
            case 3: // D0 MESON
              if( IsD0PiPlPiMiPiZeroDaughterAOD(AODMCTrackArray,labelNegPion)&& negPionIsPrimary)
                fHistoTrueNegPionFromNeutralMesonPt[fiCut]->Fill(negPionCandidate->Pt(), fWeightJetJetMC);
              break;
            default:
              AliError(Form("Heavy neutral meson not correctly selected (only 0,1,2,3 valid)... selected: %d",fSelectedHeavyNeutralMeson));
            }
          }
        }
      }
    }
  }

  for(UInt_t i = 0; i < fSelectorPosPionIndex.size(); i++){
    AliAODTrack* posPionCandidate = dynamic_cast<AliAODTrack*>(fInputEvent->GetTrack(fSelectorPosPionIndex[i]));
    if(! ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedAOD(posPionCandidate) ) continue;
    lGoodPosPionIndexPrev.push_back(   fSelectorPosPionIndex[i]  );

    TLorentzVector* posPionforHandler = new TLorentzVector();
    posPionforHandler->SetPxPyPzE(posPionCandidate->Px(), posPionCandidate->Py(), posPionCandidate->Pz(), posPionCandidate->E());
    FixPzVecToMatchPDGInvMass(posPionforHandler);
    AliAODConversionPhoton *posPionHandler = new AliAODConversionPhoton(posPionforHandler);
    delete posPionforHandler;

    fPosPionCandidates->Add(posPionHandler);
    if(!fDoLightOutput){
        fHistoPosPionPt[fiCut]->Fill( posPionCandidate->Pt(), fWeightJetJetMC );
        fHistoPosPionPhi[fiCut]->Fill( posPionCandidate->Phi(), fWeightJetJetMC );
    }
    if( fMCEvent ) {
      const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
      Double_t mcProdVtxX 	= primVtxMC->GetX();
      Double_t mcProdVtxY 	= primVtxMC->GetY();
      Double_t mcProdVtxZ 	= primVtxMC->GetZ();

      Int_t labelPosPion = TMath::Abs( posPionCandidate->GetLabel() );
        if( labelPosPion>-1 && labelPosPion < fMCEvent->GetNumberOfTracks() ) {
        AliAODMCParticle* posPion = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(labelPosPion));
        Bool_t posPionIsPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD( fInputEvent, posPion, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
        if( posPion->GetPdgCode() ==  211 ){
          if(!fDoLightOutput){
            if( posPionIsPrimary ){
              fHistoTruePosPionPt[fiCut]->Fill(posPionCandidate->Pt(), fWeightJetJetMC);
            }
            switch( fSelectedHeavyNeutralMeson ) {
            case 0: // ETA MESON
              if( IsEtaPiPlPiMiPiZeroDaughterAOD(AODMCTrackArray,labelPosPion) && posPionIsPrimary )
                fHistoTruePosPionFromNeutralMesonPt[fiCut]->Fill(posPionCandidate->Pt(), fWeightJetJetMC);
              break;
            case 1: // OMEGA MESON
              if( IsOmegaPiPlPiMiPiZeroDaughterAOD(AODMCTrackArray, labelPosPion) && posPionIsPrimary)
                fHistoTruePosPionFromNeutralMesonPt[fiCut]->Fill(posPionCandidate->Pt(), fWeightJetJetMC);
              break;
            case 2: // ETA PRIME MESON
              if( IsEtaPrimePiPlPiMiEtaDaughterAOD(AODMCTrackArray, labelPosPion) &&  posPionIsPrimary)
                fHistoTruePosPionFromNeutralMesonPt[fiCut]->Fill(posPionCandidate->Pt(), fWeightJetJetMC);
              break;
            case 3: // D0 MESON
              if( IsD0PiPlPiMiPiZeroDaughterAOD(AODMCTrackArray,labelPosPion) && posPionIsPrimary )
                fHistoTruePosPionFromNeutralMesonPt[fiCut]->Fill(posPionCandidate->Pt(), fWeightJetJetMC);
              break;
            default:
              AliError(Form("Heavy neutral meson not correctly selected (only 0,1,2,3 valid)... selected: %d",fSelectedHeavyNeutralMeson));
            }
          }
        }
      }
    }
  }

//AliInfo(Form("Number of good neg pions: %i \t pos pions = %i", lGoodNegPionIndexPrev.size(), lGoodPosPionIndexPrev.size()));

  for(UInt_t i = 0; i < lGoodNegPionIndexPrev.size(); i++){
    AliVTrack* negPionCandidate = dynamic_cast<AliVTrack*>(fInputEvent->GetTrack(lGoodNegPionIndexPrev[i]));
    AliAODTrack* negPionCandidateAOD = dynamic_cast<AliAODTrack*>(fInputEvent->GetTrack(lGoodNegPionIndexPrev[i]));
    AliKFParticle negPionCandidateKF( *negPionCandidate, 211 );

    for(UInt_t j = 0; j < lGoodPosPionIndexPrev.size(); j++){
      AliVTrack *posPionCandidate = dynamic_cast<AliVTrack*>(fInputEvent->GetTrack(lGoodPosPionIndexPrev[j]));
      AliAODTrack* posPionCandidateAOD = dynamic_cast<AliAODTrack*>(fInputEvent->GetTrack(lGoodPosPionIndexPrev[j]));
      AliKFParticle posPionCandidateKF( *posPionCandidate, 211 );
      AliKFConversionPhoton* virtualPhoton = NULL;
      virtualPhoton = new AliKFConversionPhoton(negPionCandidateKF,posPionCandidateKF);

      //AliKFVertex primaryVertex(*fInputEvent->GetPrimaryVertex());
      // primaryVertexImproved+=*virtualPhoton;
     // virtualPhoton->SetProductionVertex(primaryVertex);

      virtualPhoton->SetTrackLabels( lGoodPosPionIndexPrev[j], lGoodNegPionIndexPrev[i]);

      TLorentzVector posPionVec4;
      TLorentzVector negPionVec4;
      TLorentzVector virtPionVec4;
      TLorentzVector posKFPionVec4;
      TLorentzVector negKFPionVec4;
      TLorentzVector virtKFPionVec4;

      if(use4vecformass){
        posKFPionVec4.SetPxPyPzE(posPionCandidateKF.Px(),posPionCandidateKF.Py(),posPionCandidateKF.Pz(),posPionCandidateKF.E());
        negKFPionVec4.SetPxPyPzE(negPionCandidateKF.Px(),negPionCandidateKF.Py(),negPionCandidateKF.Pz(),negPionCandidateKF.E());
        virtKFPionVec4 = posKFPionVec4 + negKFPionVec4;
      }

      Int_t labeln=0;
      Int_t labelp=0;
      Int_t motherlabelp = 0;
      Int_t motherlabeln = 0;
      AliAODMCParticle *fNegativeMCParticle =NULL;
      AliAODMCParticle *fPositiveMCParticle =NULL;
      if( fMCEvent ) {
        labeln=TMath::Abs(negPionCandidate->GetLabel());
        labelp=TMath::Abs(posPionCandidate->GetLabel());
                if(labeln>-1) fNegativeMCParticle = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(labeln));
                if(labelp>-1) fPositiveMCParticle = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(labelp));
        // check whether MC particles exist, else abort
        if (fNegativeMCParticle == NULL || fPositiveMCParticle == NULL) return;

        motherlabeln = fNegativeMCParticle->GetMother();
        motherlabelp = fPositiveMCParticle->GetMother();
        virtualPhoton->SetMCLabelPositive(labelp);
        virtualPhoton->SetMCLabelNegative(labeln);

      }
      // AliInfo(Form("virtualPhoton chi2 = %f", virtualPhoton->GetChi2perNDF()));
      //AliInfo(Form("virtualPhoton distance pi+ pi- = %f", negPionCandidateKF.GetDeviationFromParticle(posPionCandidateKF)));

      AliAODConversionPhoton *vParticle = new AliAODConversionPhoton(virtualPhoton); //To apply mass 2 pion mass cut

      if(use4vecformass){
        vParticle->SetPxPyPzE(virtKFPionVec4.Px(),virtKFPionVec4.Py(),virtKFPionVec4.Pz(),virtKFPionVec4.E());
        vParticle->SetMass(vParticle->M());
      }
      if(!fDoLightOutput){
        Bool_t isPiMiGlobalC = negPionCandidateAOD->IsGlobalConstrained();
        Bool_t isPiPlGlobalC = posPionCandidateAOD->IsGlobalConstrained();

        Double_t ds,dsp;
        posPionCandidateKF.GetDStoParticle(negPionCandidateKF,ds,dsp);
        //AliInfo(Form("ds = %f dsp = %f", ds,dsp));
        //AliInfo(Form("Is pi+ constrained = %i Is pi- constrained = %i",isPiPlGlobalC,isPiMiGlobalC));
        Float_t chi2 = virtualPhoton->GetChi2perNDF();

        if(chi2>299) chi2 = 299; // to illustrate overflow bin
        fHistovParticleChi2PerNDF[fiCut]->Fill(chi2);
        fHistovParticledS[fiCut]->Fill(ds);

        if(isPiMiGlobalC == isPiPlGlobalC){
          fHistovParticleChi2PerNDFBothConstrained[fiCut]->Fill(chi2);
          fHistovParticledSBothConstrained[fiCut]->Fill(ds);
        } else{
          fHistovParticleChi2PerNDFOneConstrained[fiCut]->Fill(chi2);
          fHistovParticledSOneConstrained[fiCut]->Fill(ds);
        }
        if (fMCEvent &&(fDoMesonQA>0)){
          if (fPositiveMCParticle && fNegativeMCParticle ) {
            if (((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->DoMassCut()){
              if (vParticle->GetMass() < ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetMassCut()){
                if(TMath::Abs(fNegativeMCParticle->GetPdgCode())==211 && TMath::Abs(fPositiveMCParticle->GetPdgCode())==211){  // Pions ...
                  fHistoTruePionPionInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt(), fWeightJetJetMC);
                  fHistoTruevParticleChi2PerNDF[fiCut]->Fill(chi2);
                  fHistoTruevParticledS[fiCut]->Fill(ds);
                  if (motherlabeln == motherlabelp){
                    fHistoTruePionPionFromSameMotherInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt(), fWeightJetJetMC);
                    fHistoTruevParticleFromSameMotherChi2PerNDF[fiCut]->Fill(chi2);
                    fHistoTruevParticleFromSameMotherdS[fiCut]->Fill(ds);
                    switch( fSelectedHeavyNeutralMeson ) {
                    case 0: // ETA MESON
                      if( IsEtaPiPlPiMiPiZeroDaughterAOD(AODMCTrackArray,labeln) ){
                        fHistoTruePionPionFromHNMInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt(), fWeightJetJetMC);
                        fHistoTruevParticleFromHNMChi2PerNDF[fiCut]->Fill(chi2);
                        fHistoTruevParticleFromHNMdS[fiCut]->Fill(ds);
                      }
                      break;
                    case 1: // OMEGA MESON
                      if( IsOmegaPiPlPiMiPiZeroDaughterAOD(AODMCTrackArray,labeln) ){
                        fHistoTruePionPionFromHNMInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt(), fWeightJetJetMC);
                        fHistoTruevParticleFromHNMChi2PerNDF[fiCut]->Fill(chi2);
                        fHistoTruevParticleFromHNMdS[fiCut]->Fill(ds);
                      }
                      break;
                    case 2: // ETA PRIME MESON
                      if( IsEtaPrimePiPlPiMiEtaDaughterAOD(AODMCTrackArray,labeln) ){
                        fHistoTruePionPionFromHNMInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt(), fWeightJetJetMC);
                        fHistoTruevParticleFromHNMChi2PerNDF[fiCut]->Fill(chi2);
                        fHistoTruevParticleFromHNMdS[fiCut]->Fill(ds);
                      }
                      break;
                    case 3: // D0 MESON
                      if( IsD0PiPlPiMiPiZeroDaughterAOD(AODMCTrackArray,labeln) ){
                        fHistoTruePionPionFromHNMInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt(), fWeightJetJetMC);
                        fHistoTruevParticleFromHNMChi2PerNDF[fiCut]->Fill(chi2);
                        fHistoTruevParticleFromHNMdS[fiCut]->Fill(ds);
                      }
                      break;
                    default:
                      AliError(Form("Heavy neutral meson not correctly selected (only 0,1,2,3 valid)... selected: %d",fSelectedHeavyNeutralMeson));
                    }
                  }
                }
              }
            } else {
              if(TMath::Abs(fNegativeMCParticle->GetPdgCode())==211 && TMath::Abs(fPositiveMCParticle->GetPdgCode())==211){  // Pions ...
                fHistoTruePionPionInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt(), fWeightJetJetMC);
                fHistoTruevParticleChi2PerNDF[fiCut]->Fill(chi2);
                fHistoTruevParticledS[fiCut]->Fill(ds);
                if (motherlabeln == motherlabelp){
                  fHistoTruePionPionFromSameMotherInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt(), fWeightJetJetMC);
                  fHistoTruevParticleFromSameMotherChi2PerNDF[fiCut]->Fill(chi2);
                  fHistoTruevParticleFromSameMotherdS[fiCut]->Fill(ds);
                  switch( fSelectedHeavyNeutralMeson ) {
                  case 0: // ETA MESON
                    if( IsEtaPiPlPiMiPiZeroDaughterAOD(AODMCTrackArray,labeln) ){
                      fHistoTruePionPionFromHNMInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt(), fWeightJetJetMC);
                      fHistoTruevParticleFromHNMChi2PerNDF[fiCut]->Fill(chi2);
                      fHistoTruevParticleFromHNMdS[fiCut]->Fill(ds);
                    }
                    break;
                  case 1: // OMEGA MESON
                    if( IsOmegaPiPlPiMiPiZeroDaughterAOD(AODMCTrackArray,labeln) ){
                      fHistoTruePionPionFromHNMInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt(), fWeightJetJetMC);
                      fHistoTruevParticleFromHNMChi2PerNDF[fiCut]->Fill(chi2);
                      fHistoTruevParticleFromHNMdS[fiCut]->Fill(ds);
                    }
                    break;
                  case 2: // ETA PRIME MESON
                    if( IsEtaPrimePiPlPiMiEtaDaughterAOD(AODMCTrackArray,labeln) ){
                      fHistoTruePionPionFromHNMInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt(), fWeightJetJetMC);
                      fHistoTruevParticleFromHNMChi2PerNDF[fiCut]->Fill(chi2);
                      fHistoTruevParticleFromHNMdS[fiCut]->Fill(ds);
                    }
                    break;
                  case 3: // D0 MESON
                    if( IsD0PiPlPiMiPiZeroDaughterAOD(AODMCTrackArray,labeln) ){
                      fHistoTruePionPionFromHNMInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt(), fWeightJetJetMC);
                      fHistoTruevParticleFromHNMChi2PerNDF[fiCut]->Fill(chi2);
                      fHistoTruevParticleFromHNMdS[fiCut]->Fill(ds);
                    }
                    break;
                  default:
                    AliError(Form("Heavy neutral meson not correctly selected (only 0,1,2,3 valid)... selected: %d",fSelectedHeavyNeutralMeson));
                  }
                }
              }
            }
          }
        }
      }
      Bool_t survivesMassCut = kFALSE;

      if (((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->DoMassCut()){
        if (vParticle->GetMass() < ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetMassCut()){
          survivesMassCut = kTRUE;
        }
      } else{
        survivesMassCut = kTRUE;
      }

      if(survivesMassCut){
          if(!fDoLightOutput){
            fHistoPionPionInvMassPt[fiCut]->Fill( vParticle->GetMass(),vParticle->Pt(), fWeightJetJetMC);
          }
          CalculateMesonCandidates(vParticle);
          delete vParticle;
          vParticle=0x0;
      }else{
          delete vParticle;
          vParticle=0x0;
      }
      delete virtualPhoton;
      virtualPhoton=0x0;
    }
  }

  Double_t clsToFPos = -1.0;
  Double_t clsToFNeg = -1.0;

  Float_t dcaToVertexXYPos = -1.0;
  Float_t dcaToVertexZPos  = -1.0;
  Float_t dcaToVertexXYNeg = -1.0;
  Float_t dcaToVertexZNeg  = -1.0;

  if ( fDoMesonQA>0 ) {
    for(UInt_t i = 0; i < lGoodNegPionIndexPrev.size(); i++){
      AliVTrack* negPionCandidate = dynamic_cast<AliVTrack*>(fInputEvent->GetTrack(lGoodNegPionIndexPrev[i]));

      clsToFNeg = ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetNFindableClustersTPC(negPionCandidate);

      Float_t bNeg[2];
      Float_t bCovNeg[3];
      negPionCandidate->GetImpactParameters(bNeg,bCovNeg);
      if (bCovNeg[0]<=0 || bCovNeg[2]<=0) {
        AliDebug(1, "Estimated b resolution lower or equal zero!");
        bCovNeg[0]=0; bCovNeg[2]=0;
      }

      dcaToVertexXYNeg = bNeg[0];
      dcaToVertexZNeg  = bNeg[1];

      if(!fDoLightOutput){
        fHistoNegPionEta[fiCut]->Fill(negPionCandidate->Eta(), fWeightJetJetMC);
        fHistoNegPionClsTPC[fiCut]->Fill(clsToFNeg,negPionCandidate->Pt(), fWeightJetJetMC);

        fHistoPionDCAxy[fiCut]->Fill(  dcaToVertexXYNeg, negPionCandidate->Pt(), fWeightJetJetMC );
        fHistoPionDCAz[fiCut]->Fill(   dcaToVertexZNeg,  negPionCandidate->Pt(), fWeightJetJetMC );

        fHistoPionTPCdEdxNSigma[fiCut]->Fill( negPionCandidate->P(),((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetPIDResponse()->NumberOfSigmasTPC(negPionCandidate, AliPID::kPion), fWeightJetJetMC );

        fHistoPionTPCdEdx[fiCut]->Fill(negPionCandidate->P(), TMath::Abs(negPionCandidate->GetTPCsignal()), fWeightJetJetMC);
      }
    }

    for(UInt_t i = 0; i < lGoodPosPionIndexPrev.size(); i++){
      AliVTrack* posPionCandidate = dynamic_cast<AliVTrack*>(fInputEvent->GetTrack(lGoodPosPionIndexPrev[i]));

      clsToFPos = ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetNFindableClustersTPC(posPionCandidate);

      Float_t bPos[2];
      Float_t bCovPos[3];
      posPionCandidate->GetImpactParameters(bPos,bCovPos);
      if (bCovPos[0]<=0 || bCovPos[2]<=0) {
        AliDebug(1, "Estimated b resolution lower or equal zero!");
        bCovPos[0]=0; bCovPos[2]=0;
      }

      dcaToVertexXYPos = bPos[0];
      dcaToVertexZPos  = bPos[1];

      if(!fDoLightOutput){
        fHistoPosPionEta[fiCut]->Fill(posPionCandidate->Eta(), fWeightJetJetMC);
        fHistoPosPionClsTPC[fiCut]->Fill(clsToFPos,posPionCandidate->Pt(), fWeightJetJetMC);

        fHistoPionDCAxy[fiCut]->Fill(  dcaToVertexXYPos, posPionCandidate->Pt(), fWeightJetJetMC );
        fHistoPionDCAz[fiCut]->Fill(   dcaToVertexZPos,  posPionCandidate->Pt(), fWeightJetJetMC );

        fHistoPionTPCdEdxNSigma[fiCut]->Fill( posPionCandidate->P(),((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetPIDResponse()->NumberOfSigmasTPC(posPionCandidate, AliPID::kPion), fWeightJetJetMC );

        fHistoPionTPCdEdx[fiCut]->Fill(posPionCandidate->P(), TMath::Abs(posPionCandidate->GetTPCsignal()), fWeightJetJetMC);
      }
    }
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::ProcessMCParticles(){

  // Loop over all primary MC particle
  const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX 	= primVtxMC->GetX();
  Double_t mcProdVtxY 	= primVtxMC->GetY();
  Double_t mcProdVtxZ 	= primVtxMC->GetZ();

  for(Int_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {
    if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, i, mcProdVtxX, mcProdVtxY, mcProdVtxZ)){
      Double_t tempParticleWeight       = fWeightJetJetMC;

      TParticle* particle = (TParticle *)fMCEvent->Particle(i);
      if (!particle) continue;

      Int_t isMCFromMBHeader = -1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        isMCFromMBHeader
          = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
        if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
        // Set the jetjet weight to 1 in case the particle orignated from the minimum bias header
        if(isMCFromMBHeader == 2 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4) tempParticleWeight = 1;
      }

      if(!fDoLightOutput){
        if(fDoMesonQA > 0){
          // Fill kinematics for heavy particle
          // This also contains not-reconstructed particles
          if(TMath::Abs(particle->GetPdgCode()) == fPDGCodeAnalyzedMeson){
            fHistoMCHeavyAllPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
            fHistoMCHeavyAllEta[fiCut]->Fill(particle->Eta(), tempParticleWeight);
            fHistoMCHeavyAllPhi[fiCut]->Fill(TVector2::Phi_0_2pi(particle->Phi()), tempParticleWeight);
            // check for Decay kinematics
            if(particle->GetNDaughters() == 3){
              AliVParticle *neutralMeson(nullptr), *piplus(nullptr), *piminus(nullptr);
              Int_t indexpiplus(-1), indexpiminus(-1);
              for(int idaug = particle->GetFirstDaughter(); idaug <= particle->GetLastDaughter(); idaug++) {
                AliVParticle *daughter = fMCEvent->GetTrack(idaug);
                if(daughter->PdgCode() == kPiPlus){
                  piplus = daughter;
                  indexpiplus = idaug;
                }  else if(daughter->PdgCode() == kPiMinus) {
                  piminus = daughter;
                  indexpiminus = idaug;
                } else if(TMath::Abs(daughter->PdgCode()) == fPDGCodeNDM){
                  neutralMeson = daughter;
                }
              }
              if(neutralMeson && piplus && piminus) {
                fHistoMCHeavyChannelPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
                fHistoMCHeavyChannelEta[fiCut]->Fill(particle->Eta(), tempParticleWeight);
                fHistoMCHeavyChannelPhi[fiCut]->Fill(TVector2::Phi_0_2pi(particle->Phi()), tempParticleWeight);
                // Fill kinematics for daughter particles
                fHistMCChannelNDMFromHeavyPt[fiCut]->Fill(neutralMeson->Pt(), tempParticleWeight);
                fHistMCChannelNDMFromHeavyEta[fiCut]->Fill(neutralMeson->Eta(), tempParticleWeight);
                fHistMCChannelNDMFromHeavyPhi[fiCut]->Fill(TVector2::Phi_0_2pi(neutralMeson->Phi()), tempParticleWeight);
                fHistMCChannelPiPlusFromHeavyPt[fiCut]->Fill(piplus->Pt(), tempParticleWeight);
                fHistMCChannelPiPlusFromHeavyEta[fiCut]->Fill(piplus->Eta(), tempParticleWeight);
                fHistMCChannelPiPlusFromHeavyPhi[fiCut]->Fill(TVector2::Phi_0_2pi(piplus->Phi()), tempParticleWeight);
                fHistMCChannelPiMinusFromHeavyPt[fiCut]->Fill(piminus->Pt(), tempParticleWeight);
                fHistMCChannelPiMinusFromHeavyEta[fiCut]->Fill(piminus->Eta(), tempParticleWeight);
                fHistMCChannelPiPMinusFromHeavyPhi[fiCut]->Fill(TVector2::Phi_0_2pi(piminus->Phi()), tempParticleWeight);
                fHistMCChannelNDMPtHeavyPt[fiCut]->Fill(particle->Pt(), neutralMeson->Pt(), tempParticleWeight);
                fHistMCChannelPiPlusPtHeavyPt[fiCut]->Fill(particle->Pt(), piplus->Pt(), tempParticleWeight);
                fHistMCChannelPiMinusPtHeavyPt[fiCut]->Fill(particle->Pt(), piminus->Pt(), tempParticleWeight);

                // check if particle is reconstructible
                bool reconstructible(true);
                if(!((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedMC(indexpiminus,fMCEvent)) reconstructible = false;
                if(!((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedMC(indexpiplus,fMCEvent)) reconstructible = false;
                if(neutralMeson->GetNDaughters() == 3) {
                  // exclude Dalitz-decays
                  reconstructible = false;
                } else {
                  AliVParticle *photon1 = fMCEvent->GetTrack(neutralMeson->GetDaughterFirst()), *photon2 = fMCEvent->GetTrack(neutralMeson->GetDaughterLast());
                  if(!(photon1 && photon2)) reconstructible = false;
                  else {
                    switch(fNDMRecoMode) {
                      case 0 : {
                        if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedMC(photon1->Particle(),fMCEvent,kFALSE)) reconstructible = false;
                        if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedMC(photon2->Particle(),fMCEvent,kFALSE)) reconstructible = false;
                        break;
                      }
                      case 1: {
                        if(!(((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedMC(photon1->Particle(),fMCEvent,kFALSE)  &&
                             ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(photon1->Particle(),fMCEvent)) ||
                           !(((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedMC(photon2->Particle(),fMCEvent,kFALSE)  &&
                             ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(photon1->Particle(),fMCEvent))
                        ) reconstructible = false;
                        break;
                      }
                      case 2: {
                        if(!((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(photon1->Particle(),fMCEvent)) reconstructible = false;
                        if(!((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(photon1->Particle(),fMCEvent)) reconstructible = false;
                        break;
                      }
                    };
                  }
                }
                if(reconstructible) {
                  fHistoMCHeavyReconstructiblePt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
                  fHistoMCHeavyReconstructibleEta[fiCut]->Fill(particle->Eta(), tempParticleWeight);
                  fHistoMCHeavyReconstructiblePhi[fiCut]->Fill(TVector2::Phi_0_2pi(particle->Phi()), tempParticleWeight);
                  // Fill kinematics for daughter particles
                  fHistMCReconstructibleNDMFromHeavyPt[fiCut]->Fill(neutralMeson->Pt(), tempParticleWeight);
                  fHistMCReconstructibleNDMFromHeavyEta[fiCut]->Fill(neutralMeson->Eta(), tempParticleWeight);
                  fHistMCReconstructibleNDMFromHeavyPhi[fiCut]->Fill(TVector2::Phi_0_2pi(neutralMeson->Phi()), tempParticleWeight);
                  fHistMCReconstructiblePiPlusFromHeavyPt[fiCut]->Fill(piplus->Pt(), tempParticleWeight);
                  fHistMCReconstructiblePiPlusFromHeavyEta[fiCut]->Fill(piplus->Eta(), tempParticleWeight);
                  fHistMCReconstructiblePiPlusFromHeavyPhi[fiCut]->Fill(TVector2::Phi_0_2pi(piplus->Phi()), tempParticleWeight);
                  fHistMCReconstructiblePiMinusFromHeavyPt[fiCut]->Fill(piminus->Pt(), tempParticleWeight);
                  fHistMCReconstructiblePiMinusFromHeavyEta[fiCut]->Fill(piminus->Eta(), tempParticleWeight);
                  fHistMCReconstructiblePiPMinusFromHeavyPhi[fiCut]->Fill(TVector2::Phi_0_2pi(piminus->Phi()), tempParticleWeight);
                  fHistMCReconstructibleNDMPtHeavyPt[fiCut]->Fill(particle->Pt(), neutralMeson->Pt(), tempParticleWeight);
                  fHistMCReconstructiblePiPlusPtHeavyPt[fiCut]->Fill(particle->Pt(), piplus->Pt(), tempParticleWeight);
                  fHistMCReconstructiblePiMinusPtHeavyPt[fiCut]->Fill(particle->Pt(), piminus->Pt(), tempParticleWeight);
                }
              }
            }
          }
        }

        if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent,fInputEvent)){
          // find MC photons
          if (fNDMRecoMode < 2 ){
            if(((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCEvent,kFALSE)){
              fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt(), tempParticleWeight); // All MC Gamma
              if(particle->GetMother(0) >-1){
                if (fMCEvent->Particle(particle->GetMother(0))->GetPdgCode() ==fPDGCodeNDM){
                  TParticle *particleNDM = fMCEvent->Particle(particle->GetMother(0));
                  if(fDoMesonQA){
                    fHistoMCAllMesonPt[fiCut]->Fill(particleNDM->Pt(), tempParticleWeight);
                    fHistoMCAllMesonEta[fiCut]->Fill(particleNDM->Eta(), tempParticleWeight);
                    fHistoMCAllMesonPhi[fiCut]->Fill(TVector2::Phi_0_2pi(particleNDM->Phi()), tempParticleWeight);
                  }
                  if (fMCEvent->Particle(particle->GetMother(0))->GetMother(0) > -1){
                    if ( fMCEvent->Particle((fMCEvent->Particle(particle->GetMother(0)))->GetMother(0))->GetPdgCode() == fPDGCodeAnalyzedMeson ){
                      if ( fMCEvent->Particle((fMCEvent->Particle(particle->GetMother(0)))->GetMother(0))->GetNDaughters()==3 ){
                        fHistoMCGammaFromNeutralMesonPt[fiCut]->Fill(particle->Pt(), tempParticleWeight); // All photons from eta or omega via pi0
                        if(fDoMesonQA > 0) {
                          fHistoMCMesonFromNeutralMesonPt[fiCut]->Fill(particleNDM->Pt(), tempParticleWeight);   // ALl meson mothers (pi0/eta) from analyzed heavy meson
                          fHistoMCMesonFromNeutralMesonEta[fiCut]->Fill(particleNDM->Eta(), tempParticleWeight);   // All meson mothers (pi0/eta) from analyzed heavy meson
                          fHistoMCMesonFromNeutralMesonPhi[fiCut]->Fill(TVector2::Phi_0_2pi(particleNDM->Phi()), tempParticleWeight);   // All meson mothers (pi0/eta) from analyzed heavy meson
                        }
                      }
                    }
                  }
                }
              }
            }
          } else if (fNDMRecoMode == 2){
            if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(particle,fMCEvent)){
              fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt(), tempParticleWeight); // All MC Gamma
              if(particle->GetMother(0) >-1){
                if (fMCEvent->Particle(particle->GetMother(0))->GetPdgCode() ==fPDGCodeNDM){
                  TParticle *particleNDM = fMCEvent->Particle(particle->GetMother(0));
                  if(fDoMesonQA > 0){
                    fHistoMCAllMesonPt[fiCut]->Fill(particleNDM->Pt(), tempParticleWeight);
                    fHistoMCAllMesonEta[fiCut]->Fill(particleNDM->Eta(), tempParticleWeight);
                    fHistoMCAllMesonPhi[fiCut]->Fill(TVector2::Phi_0_2pi(particleNDM->Phi()), tempParticleWeight);
                  }
                  if (fMCEvent->Particle(particle->GetMother(0))->GetMother(0) > -1){
                    if ( fMCEvent->Particle((fMCEvent->Particle(particle->GetMother(0)))->GetMother(0))->GetPdgCode() == fPDGCodeAnalyzedMeson ){
                      if ( fMCEvent->Particle((fMCEvent->Particle(particle->GetMother(0)))->GetMother(0))->GetNDaughters()==3 ){
                        fHistoMCGammaFromNeutralMesonPt[fiCut]->Fill(particle->Pt(), tempParticleWeight); // All photons from analyzed meson via pi0 or eta from decay
                        if(fDoMesonQA > 0){
                          fHistoMCMesonFromNeutralMesonPt[fiCut]->Fill(particleNDM->Pt(), tempParticleWeight);   // All meson mothers (pi0/eta) from analyzed heavy meson
                          fHistoMCMesonFromNeutralMesonEta[fiCut]->Fill(particleNDM->Eta(), tempParticleWeight);   // All meson mothers (pi0/eta) from analyzed heavy meson
                          fHistoMCMesonFromNeutralMesonPhi[fiCut]->Fill(TVector2::Phi_0_2pi(particleNDM->Phi()), tempParticleWeight);   // All meson mothers (pi0/eta) from analyzed heavy meson
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          if (fNDMRecoMode < 2){
            if (((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCEvent,kTRUE)){
              fHistoMCConvGammaPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
            } // Converted MC Gamma
          }
          if(((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedMC(i,fMCEvent)){
            if( particle->GetPdgCode() == 211){
              fHistoMCAllPosPionsPt[fiCut]->Fill(particle->Pt(), tempParticleWeight); // All pos pions
              if(fDoMesonQA > 0){
                fHistoMCAllPosPionsEta[fiCut]->Fill(particle->Eta(), tempParticleWeight); // All pos pions
                fHistoMCAllPosPionsPhi[fiCut]->Fill(TVector2::Phi_0_2pi(particle->Phi()), tempParticleWeight); // All pos pions
              }
              if(particle->GetMother(0) >-1){
                if (fMCEvent->Particle(particle->GetMother(0))->GetPdgCode() ==fPDGCodeAnalyzedMeson){
                  fHistoMCPosPionsFromNeutralMesonPt[fiCut]->Fill(particle->Pt(), tempParticleWeight); // All pos from neutral heavy meson (omega, eta OR eta prime)
                  if(fDoMesonQA > 0){
                    fHistoMCPosPionsFromNeutralMesonEta[fiCut]->Fill(particle->Eta(), tempParticleWeight); // All pos from neutral heavy meson (omega, eta OR eta prime)
                    fHistoMCPosPionsFromNeutralMesonPhi[fiCut]->Fill(TVector2::Phi_0_2pi(particle->Phi()), tempParticleWeight); // All pos from neutral heavy meson (omega, eta OR eta prime)
                  }
                }
              }
            }
            if( particle->GetPdgCode() == -211){
              fHistoMCAllNegPionsPt[fiCut]->Fill(particle->Pt(), tempParticleWeight); // All neg pions
              if(fDoMesonQA > 0) {
                fHistoMCAllNegPionsEta[fiCut]->Fill(particle->Eta(), tempParticleWeight); // All neg pions
                fHistoMCAllNegPionsPhi[fiCut]->Fill(TVector2::Phi_0_2pi(particle->Phi()), tempParticleWeight); // All neg pions
              }
              if(particle->GetMother(0) >-1){
                if (fMCEvent->Particle(particle->GetMother(0))->GetPdgCode() ==fPDGCodeAnalyzedMeson){
                  fHistoMCNegPionsFromNeutralMesonPt[fiCut]->Fill(particle->Pt(), tempParticleWeight); // All pos from neutral heavy meson (omega, eta OR eta prime)
                  if(fDoMesonQA > 0){
                    fHistoMCNegPionsFromNeutralMesonEta[fiCut]->Fill(particle->Eta(), tempParticleWeight); // All pos from neutral heavy meson (omega, eta OR eta prime)
                    fHistoMCNegPionsFromNeutralMesonPhi[fiCut]->Fill(TVector2::Phi_0_2pi(particle->Phi()), tempParticleWeight); // All pos from neutral heavy meson (omega, eta OR eta prime)
                  }
                }
              }
            }
          }
        }
      }

      // \eta -> pi+ pi- \gamma
      Int_t labelNDM = -1;
      Int_t labelNegPion = -1;
      Int_t labelPosPion = -1;

      if( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedMCPiPlPiMiPiZero(particle,fMCEvent,labelNegPion,labelPosPion,labelNDM,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()) ||   ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedMCPiPlPiMiEta(particle,fMCEvent,labelNegPion,labelPosPion,labelNDM,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
        Float_t weighted= 1.;
        if( ((AliPrimaryPionCuts*) fPionCutArray->At(fiCut))->DoWeights() ) {
          if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent,fInputEvent)){
            if (particle->Pt()>0.005){
              weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(i, fMCEvent,fInputEvent);
            }
          }
        }
        if(particle->GetPdgCode() == fPDGCodeAnalyzedMeson){
          fHistoMCHNMPiPlPiMiNDMPt[fiCut]->Fill(particle->Pt(), weighted* tempParticleWeight); 	// All MC eta, omega OR eta prime in respective decay channel
          if(!fDoLightOutput){
            fHistoMCHNMPiPlPiMiNDMEta[fiCut]->Fill(particle->Eta(),weighted* tempParticleWeight);
            fHistoMCHNMPiPlPiMiNDMPhi[fiCut]->Fill(particle->Phi(),weighted* tempParticleWeight);
            if(fDoMesonQA > 1){
              fHistoMCHNMPiPlPiMiNDMEtavsPt[fiCut]->Fill(particle->Pt(),particle->Eta(),weighted* tempParticleWeight);
            }
          }
        }
        if(labelNDM>-1){
          TParticle *particleNDM    = fMCEvent->Particle(labelNDM);
          if(particleNDM->GetDaughter(0)>-1 && particleNDM->GetDaughter(1)>-1){
            TParticle *gamma1 = fMCEvent->Particle(particleNDM->GetDaughter(0));
            TParticle *gamma2 = fMCEvent->Particle(particleNDM->GetDaughter(1));
            Bool_t kDaughter0IsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, particleNDM->GetDaughter(0), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
            Bool_t kDaughter1IsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, particleNDM->GetDaughter(1), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
            Bool_t kNegPionIsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, labelNegPion, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
            Bool_t kPosPionIsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, labelPosPion, mcProdVtxX, mcProdVtxY, mcProdVtxZ);

            if (fNDMRecoMode == 0){
              if( kDaughter0IsPrim && kDaughter1IsPrim && kNegPionIsPrim && kPosPionIsPrim &&
                  ((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedMC(gamma1,fMCEvent,kFALSE) &&					// test first daugther of pi0
                  ((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedMC(gamma2,fMCEvent,kFALSE) &&					// test second daughter of pi0
                  ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedMC(labelNegPion,fMCEvent) &&								// test negative pion
                  ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedMC(labelPosPion,fMCEvent) 								// test positive pion
              ) {
                if(particle->GetPdgCode() == fPDGCodeAnalyzedMeson){
                  fHistoMCHNMPiPlPiMiNDMInAccPt[fiCut]->Fill(particle->Pt(), weighted* tempParticleWeight ); 		// MC Eta, omega or eta prime with pi+ pi- pi0 with gamma's and e+e- in acc

                  // check relation between HNM pt and NDM, while respecting pT cutoff for NDM
                  if(!fDoLightOutput){
                    if(particleNDM->Pt() >= fNDMMinPtPossible) fHistoMCHNMInAccVsNDMPt[fiCut]->Fill(particle->Pt(),particleNDM->Pt(),weighted* tempParticleWeight);
                  }
                }

              }
            } else if (fNDMRecoMode == 1){ // mixed mode
              // check acceptamce of pions firs
              if( kDaughter0IsPrim && kDaughter1IsPrim && kNegPionIsPrim && kPosPionIsPrim &&
                      ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedMC(labelNegPion,fMCEvent) &&								// test negative pion
                      ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedMC(labelPosPion,fMCEvent) 								// test positive pion
                ) {
                // check acceptance of clusters and PCM photons
                if((((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedMC(gamma1,fMCEvent,kFALSE)	&&
                   ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(gamma2,fMCEvent)) ||
                   (((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedMC(gamma2,fMCEvent,kFALSE)	&&
                   ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(gamma1,fMCEvent))
                   ){
                      if(particle->GetPdgCode() == fPDGCodeAnalyzedMeson){
                        fHistoMCHNMPiPlPiMiNDMInAccPt[fiCut]->Fill(particle->Pt(), weighted* tempParticleWeight  ); 		// MC Eta, omega or eta prime with pi+ pi- pi0 with gamma's and e+e- in acc

                        // check relation between HNM pt and NDM, while respecting pT cutoff for NDM
                        if(!fDoLightOutput){
                          if(particleNDM->Pt() >= fNDMMinPtPossible) fHistoMCHNMInAccVsNDMPt[fiCut]->Fill(particle->Pt(),particleNDM->Pt(),weighted* tempParticleWeight );
                        }
                      }
                   }
              }
            } else if (fNDMRecoMode == 2){
              if( kDaughter0IsPrim && kDaughter1IsPrim && kNegPionIsPrim && kPosPionIsPrim &&
                  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(gamma1,fMCEvent) &&					// test first daugther of pi0
                  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(gamma2,fMCEvent) &&					// test second daughter of pi0
                  ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedMC(labelNegPion,fMCEvent) &&								// test negative pion
                  ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedMC(labelPosPion,fMCEvent) 								// test positive pion
              ) {
                if(particle->GetPdgCode() == fPDGCodeAnalyzedMeson){
                  fHistoMCHNMPiPlPiMiNDMInAccPt[fiCut]->Fill(particle->Pt(), weighted* tempParticleWeight  ); 		// MC Eta pi+ pi- pi0 with gamma's and e+e- in acc

                  // check relation between HNM pt and NDM, while respecting pT cutoff for NDM
                  if(!fDoLightOutput){
                    if(particleNDM->Pt() >= fNDMMinPtPossible) fHistoMCHNMInAccVsNDMPt[fiCut]->Fill(particle->Pt(),particleNDM->Pt(),weighted* tempParticleWeight );
                  }
                }
              }
            }
          }
        }
      }
    }

  }
}

//_____________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::ProcessAODMCParticles(){
  // Loop over all primary MC particle
  const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX 	= primVtxMC->GetX();
  Double_t mcProdVtxY 	= primVtxMC->GetY();
  Double_t mcProdVtxZ 	= primVtxMC->GetZ();

  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));

  if (AODMCTrackArray){
    for(Int_t i = 0; i < AODMCTrackArray->GetEntriesFast(); i++) {
      Double_t tempParticleWeight       = fWeightJetJetMC;
      AliAODMCParticle* particle =  static_cast<AliAODMCParticle*>(AODMCTrackArray->At(i));
      if (!particle) continue;

      if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD( fInputEvent, particle, mcProdVtxX, mcProdVtxY, mcProdVtxZ)){
        Int_t isMCFromMBHeader = -1;
        if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
          isMCFromMBHeader
              = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
          if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
          // Set the jetjet weight to 1 in case the particle orignated from the minimum bias header
          if(isMCFromMBHeader == 2 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4) tempParticleWeight = 1;
        }
        if(!fDoLightOutput){
          // Fill kinematics for heavy particle
          // This also contains not-reconstructed particles
          if(fDoMesonQA > 0){
            if(TMath::Abs(particle->GetPdgCode()) == fPDGCodeAnalyzedMeson){
              fHistoMCHeavyAllPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
              fHistoMCHeavyAllEta[fiCut]->Fill(particle->Eta(), tempParticleWeight);
              fHistoMCHeavyAllPhi[fiCut]->Fill(TVector2::Phi_0_2pi(particle->Phi()), fWeightJetJetMC);
              // check for Decay kinematics
              if(particle->GetNDaughters() == 3){
                AliAODMCParticle *neutralMeson(nullptr), *piplus(nullptr), *piminus(nullptr);
                Int_t indexpiplus(-1), indexpiminus(-1);
                for(int idaug = particle->GetDaughterFirst(); idaug <= particle->GetDaughterLast(); idaug++) {
                  AliAODMCParticle *daughter = static_cast<AliAODMCParticle *>(AODMCTrackArray->At(idaug));
                  if(daughter->PdgCode() == kPiPlus){
                    piplus = daughter;
                    indexpiplus = idaug;
                  } else if(daughter->PdgCode() == kPiMinus) {
                    piminus = daughter;
                    indexpiminus = idaug;
                  } else if(TMath::Abs(daughter->PdgCode()) == fPDGCodeNDM){
                    neutralMeson = daughter;
                  }
                }
                if(neutralMeson && piplus && piminus) {
                  // Meson is in the expected channel
                  fHistoMCHeavyChannelPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
                  fHistoMCHeavyChannelEta[fiCut]->Fill(particle->Eta(), fWeightJetJetMC);
                  fHistoMCHeavyChannelPhi[fiCut]->Fill(TVector2::Phi_0_2pi(particle->Phi()), fWeightJetJetMC);
                  // Fill kinematics for daughter particles
                  fHistMCChannelNDMFromHeavyPt[fiCut]->Fill(neutralMeson->Pt(), fWeightJetJetMC);
                  fHistMCChannelNDMFromHeavyEta[fiCut]->Fill(neutralMeson->Eta(), fWeightJetJetMC);
                  fHistMCChannelNDMFromHeavyPhi[fiCut]->Fill(TVector2::Phi_0_2pi(neutralMeson->Phi()), fWeightJetJetMC);
                  fHistMCChannelPiPlusFromHeavyPt[fiCut]->Fill(piplus->Pt(), fWeightJetJetMC);
                  fHistMCChannelPiPlusFromHeavyEta[fiCut]->Fill(piplus->Eta(), fWeightJetJetMC);
                  fHistMCChannelPiPlusFromHeavyPhi[fiCut]->Fill(TVector2::Phi_0_2pi(piplus->Phi()), fWeightJetJetMC);
                  fHistMCChannelPiMinusFromHeavyPt[fiCut]->Fill(piminus->Pt(), fWeightJetJetMC);
                  fHistMCChannelPiMinusFromHeavyEta[fiCut]->Fill(piminus->Eta(), fWeightJetJetMC);
                  fHistMCChannelPiPMinusFromHeavyPhi[fiCut]->Fill(TVector2::Phi_0_2pi(piminus->Phi()), fWeightJetJetMC);
                  fHistMCChannelNDMPtHeavyPt[fiCut]->Fill(particle->Pt(), neutralMeson->Pt(), fWeightJetJetMC);
                  fHistMCChannelPiPlusPtHeavyPt[fiCut]->Fill(particle->Pt(), piplus->Pt(), fWeightJetJetMC);
                  fHistMCChannelPiMinusPtHeavyPt[fiCut]->Fill(particle->Pt(), piminus->Pt(), fWeightJetJetMC);

                  // check if particle is reconstructible
                  bool reconstructible(true);
                  if(!((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedAODMC(indexpiminus,AODMCTrackArray)) reconstructible = false;
                  if(!((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedAODMC(indexpiplus,AODMCTrackArray)) reconstructible = false;
                  if(neutralMeson->GetNDaughters() == 3) {
                    // exclude Dalitz-decays
                    reconstructible = false;
                  } else {
                    AliAODMCParticle *photon1 = static_cast<AliAODMCParticle *>(AODMCTrackArray->At(neutralMeson->GetDaughterFirst())),
                                     *photon2 = static_cast<AliAODMCParticle *>(AODMCTrackArray->At(neutralMeson->GetDaughterLast()));
                    if(!(photon1 && photon2)) reconstructible = false;
                    else {
                      switch(fNDMRecoMode) {
                        case 0 : {
                          if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedAODMC(photon1, AODMCTrackArray, kFALSE)) reconstructible = false;
                          if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedAODMC(photon2, AODMCTrackArray, kFALSE)) reconstructible = false;
                          break;
                        }
                        case 1: {
                          if(!(((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedAODMC(photon1, AODMCTrackArray, kFALSE)  &&
                              ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(photon1, AODMCTrackArray)) ||
                             !(((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedAODMC(photon2, AODMCTrackArray, kFALSE)  &&
                              ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(photon1, AODMCTrackArray))
                          ) reconstructible = false;
                          break;
                        }
                        case 2: {
                          if(!((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(photon1, AODMCTrackArray)) reconstructible = false;
                          if(!((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(photon1, AODMCTrackArray)) reconstructible = false;
                          break;
                        }
                      };
                    }
                  }

                  if(reconstructible) {
                    fHistoMCHeavyReconstructiblePt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
                    fHistoMCHeavyReconstructibleEta[fiCut]->Fill(particle->Eta(), fWeightJetJetMC);
                    fHistoMCHeavyReconstructiblePhi[fiCut]->Fill(TVector2::Phi_0_2pi(particle->Phi()), fWeightJetJetMC);
                    // Fill kinematics for daughter particles
                    fHistMCReconstructibleNDMFromHeavyPt[fiCut]->Fill(neutralMeson->Pt(), fWeightJetJetMC);
                    fHistMCReconstructibleNDMFromHeavyEta[fiCut]->Fill(neutralMeson->Eta(), fWeightJetJetMC);
                    fHistMCReconstructibleNDMFromHeavyPhi[fiCut]->Fill(TVector2::Phi_0_2pi(neutralMeson->Phi()), fWeightJetJetMC);
                    fHistMCReconstructiblePiPlusFromHeavyPt[fiCut]->Fill(piplus->Pt(), fWeightJetJetMC);
                    fHistMCReconstructiblePiPlusFromHeavyEta[fiCut]->Fill(piplus->Eta(), fWeightJetJetMC);
                    fHistMCReconstructiblePiPlusFromHeavyPhi[fiCut]->Fill(TVector2::Phi_0_2pi(piplus->Phi()), fWeightJetJetMC);
                    fHistMCReconstructiblePiMinusFromHeavyPt[fiCut]->Fill(piminus->Pt(), fWeightJetJetMC);
                    fHistMCReconstructiblePiMinusFromHeavyEta[fiCut]->Fill(piminus->Eta(), fWeightJetJetMC);
                    fHistMCReconstructiblePiPMinusFromHeavyPhi[fiCut]->Fill(TVector2::Phi_0_2pi(piminus->Phi()), fWeightJetJetMC);
                    fHistMCReconstructibleNDMPtHeavyPt[fiCut]->Fill(particle->Pt(), neutralMeson->Pt(), fWeightJetJetMC);
                    fHistMCReconstructiblePiPlusPtHeavyPt[fiCut]->Fill(particle->Pt(), piplus->Pt(), fWeightJetJetMC);
                    fHistMCReconstructiblePiMinusPtHeavyPt[fiCut]->Fill(particle->Pt(), piminus->Pt(), fWeightJetJetMC);
                  }
                }
              }
            }
          }
          // find MC photons
          if (fNDMRecoMode < 2 ){
            if(((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedAODMC(particle,AODMCTrackArray,kFALSE)){
              fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC); // All MC Gamma
              if(particle->GetMother() >-1){
                if ((static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetMother())))->GetPdgCode() ==fPDGCodeNDM){
                  AliAODMCParticle *particleNDM = static_cast<AliAODMCParticle *>(AODMCTrackArray->At(particle->GetMother()));
                  if (fDoMesonQA>0){
                    fHistoMCAllMesonPt[fiCut]->Fill(particleNDM->Pt(), fWeightJetJetMC);
                    fHistoMCAllMesonEta[fiCut]->Fill(particleNDM->Eta(), fWeightJetJetMC);
                    fHistoMCAllMesonPhi[fiCut]->Fill(TVector2::Phi_0_2pi(particleNDM->Phi()), fWeightJetJetMC);
                  }
                  if ((static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetMother())))->GetMother() > -1){
                    if ( (static_cast<AliAODMCParticle*>(AODMCTrackArray->At( (static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetMother())))->GetMother() )))->GetPdgCode() == fPDGCodeAnalyzedMeson ){
                      if ( (static_cast<AliAODMCParticle*>(AODMCTrackArray->At( (static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetMother())))->GetMother() )))->GetNDaughters() ==3 ) {
                        fHistoMCGammaFromNeutralMesonPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC); // All photons from eta or omega via pi0
                        if(fDoMesonQA > 0) {
                          fHistoMCMesonFromNeutralMesonPt[fiCut]->Fill(particleNDM->Pt(), fWeightJetJetMC);   // ALl meson mothers (pi0/eta) from analyzed heavy meson
                          fHistoMCMesonFromNeutralMesonEta[fiCut]->Fill(particleNDM->Eta(), fWeightJetJetMC);   // All meson mothers (pi0/eta) from analyzed heavy meson
                          fHistoMCMesonFromNeutralMesonPhi[fiCut]->Fill(TVector2::Phi_0_2pi(particleNDM->Phi()), fWeightJetJetMC);   // All meson mothers (pi0/eta) from analyzed heavy meson
                        }
                      }
                    }
                  }
                }
              }
            }
          } else if (fNDMRecoMode == 2){
            if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(particle,AODMCTrackArray)){
              fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC); // All MC Gamma
              if(particle->GetMother() >-1){
                if ( (static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetMother())))->GetPdgCode() == fPDGCodeNDM){
                  AliAODMCParticle *particleNDM = static_cast<AliAODMCParticle *>(AODMCTrackArray->At(particle->GetMother()));
                  if(fDoMesonQA > 0){
                    fHistoMCAllMesonPt[fiCut]->Fill(particleNDM->Pt(), fWeightJetJetMC);
                    fHistoMCAllMesonEta[fiCut]->Fill(particleNDM->Eta(), fWeightJetJetMC);
                    fHistoMCAllMesonPhi[fiCut]->Fill(TVector2::Phi_0_2pi(particleNDM->Phi()), fWeightJetJetMC);
                  }
                  if ((static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetMother())))->GetMother() > -1){
                    if ( (static_cast<AliAODMCParticle*>(AODMCTrackArray->At( (static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetMother())))->GetMother() )))->GetPdgCode() == fPDGCodeAnalyzedMeson ){
                      if ( (static_cast<AliAODMCParticle*>(AODMCTrackArray->At( (static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetMother())))->GetMother() )))->GetNDaughters() ==3 ) {
                        fHistoMCGammaFromNeutralMesonPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC); // All photons from analyzed meson via pi0 or eta from decay
                        if(fDoMesonQA > 0){
                          fHistoMCMesonFromNeutralMesonPt[fiCut]->Fill(particleNDM->Pt(), fWeightJetJetMC);   // All meson mothers (pi0/eta) from analyzed heavy meson
                          fHistoMCMesonFromNeutralMesonEta[fiCut]->Fill(particleNDM->Eta(), fWeightJetJetMC);   // All meson mothers (pi0/eta) from analyzed heavy meson
                          fHistoMCMesonFromNeutralMesonPhi[fiCut]->Fill(TVector2::Phi_0_2pi(particleNDM->Phi()), fWeightJetJetMC);   // All meson mothers (pi0/eta) from analyzed heavy meson
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          if (fNDMRecoMode < 2){
            if (((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedAODMC(particle,AODMCTrackArray,kTRUE)){
              fHistoMCConvGammaPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
            } // Converted MC Gamma
          }
          if(((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedAODMC(i,AODMCTrackArray)){
            if( particle->GetPdgCode() == 211){
              fHistoMCAllPosPionsPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC); // All pos pions
              if(fDoMesonQA > 0){
                fHistoMCAllPosPionsEta[fiCut]->Fill(particle->Eta(), fWeightJetJetMC); // All pos pions
                fHistoMCAllPosPionsPhi[fiCut]->Fill(TVector2::Phi_0_2pi(particle->Phi()), fWeightJetJetMC); // All pos pions
              }
              if(particle->GetMother() >-1){
                if ( (static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetMother())))->GetPdgCode() ==fPDGCodeAnalyzedMeson){
                  fHistoMCPosPionsFromNeutralMesonPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC); // All pos from neutral heavy meson (omega, eta OR eta prime)
                  if(fDoMesonQA > 0){
                    fHistoMCPosPionsFromNeutralMesonEta[fiCut]->Fill(particle->Eta(), fWeightJetJetMC); // All pos from neutral heavy meson (omega, eta OR eta prime)
                    fHistoMCPosPionsFromNeutralMesonPhi[fiCut]->Fill(TVector2::Phi_0_2pi(particle->Phi()), fWeightJetJetMC); // All pos from neutral heavy meson (omega, eta OR eta prime)
                  }
                }
              }
            }
            if( particle->GetPdgCode() == -211){
              fHistoMCAllNegPionsPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC); // All neg pions
              if(fDoMesonQA){
                fHistoMCAllNegPionsEta[fiCut]->Fill(particle->Eta(), fWeightJetJetMC); // All neg pions
                fHistoMCAllNegPionsPhi[fiCut]->Fill(TVector2::Phi_0_2pi(particle->Phi()), fWeightJetJetMC); // All neg pions
              }
              if(particle->GetMother() >-1){
                if ((static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetMother())))->GetPdgCode() ==fPDGCodeAnalyzedMeson) {
                  fHistoMCNegPionsFromNeutralMesonPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC); // All pos from neutral heavy meson (omega, eta OR eta prime)
                  if(fDoMesonQA > 0) {
                    fHistoMCNegPionsFromNeutralMesonEta[fiCut]->Fill(particle->Eta(), fWeightJetJetMC); // All pos from neutral heavy meson (omega, eta OR eta prime)
                    fHistoMCNegPionsFromNeutralMesonPhi[fiCut]->Fill(TVector2::Phi_0_2pi(particle->Phi()), fWeightJetJetMC); // All pos from neutral heavy meson (omega, eta OR eta prime)
                  }
                }
              }
            }
          }
        }

        // \eta -> pi+ pi- \gamma
        Int_t labelNDM = -1;
        Int_t labelNegPion = -1;
        Int_t labelPosPion = -1;

        if( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedAODMCPiPlPiMiPiZero(particle,AODMCTrackArray,labelNegPion,labelPosPion,labelNDM,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()) ||   ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedAODMCPiPlPiMiEta(particle,AODMCTrackArray,labelNegPion,labelPosPion,labelNDM,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
          Float_t weighted= 1.;
          if( ((AliPrimaryPionCuts*) fPionCutArray->At(fiCut))->DoWeights() ) {
            if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent,fInputEvent)){
              if (particle->Pt()>0.005){
                weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(i, 0x0,fInputEvent);
              }
            }
          }
          if(particle->GetPdgCode() == fPDGCodeAnalyzedMeson)fHistoMCHNMPiPlPiMiNDMPt[fiCut]->Fill(particle->Pt(), tempParticleWeight); 	// All MC eta, omega OR eta prime in respective decay channel
          if(!fDoLightOutput){
            fHistoMCHNMPiPlPiMiNDMEta[fiCut]->Fill(particle->Eta(),tempParticleWeight);
            fHistoMCHNMPiPlPiMiNDMPhi[fiCut]->Fill(particle->Phi(),tempParticleWeight);
            if(fDoMesonQA > 1){
              fHistoMCHNMPiPlPiMiNDMEtavsPt[fiCut]->Fill(particle->Pt(),particle->Eta(),weighted* tempParticleWeight);
            }
          }
          if(labelNDM>-1){
            AliAODMCParticle* particleNDM    = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(labelNDM));
            if(particleNDM->GetDaughterLabel(0)>-1 && particleNDM->GetDaughterLabel(1)>-1){
              AliAODMCParticle *gamma1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particleNDM->GetDaughterLabel(0)));
              AliAODMCParticle *gamma2 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particleNDM->GetDaughterLabel(1)));
              AliAODMCParticle *negpion = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(labelNegPion));
              AliAODMCParticle *pospion = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(labelPosPion));
              Bool_t kDaughter0IsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD( fInputEvent, gamma1, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
              Bool_t kDaughter1IsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD( fInputEvent, gamma2, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
              Bool_t kNegPionIsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD( fInputEvent, negpion, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
              Bool_t kPosPionIsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD( fInputEvent, pospion, mcProdVtxX, mcProdVtxY, mcProdVtxZ);

              if (fNDMRecoMode == 0){
                if( kDaughter0IsPrim && kDaughter1IsPrim && kNegPionIsPrim && kPosPionIsPrim &&
                    ((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedAODMC(gamma1,AODMCTrackArray,kFALSE) &&					// test first daugther of pi0
                    ((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedAODMC(gamma2,AODMCTrackArray,kFALSE) &&					// test second daughter of pi0
                    ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedAODMC(labelNegPion,AODMCTrackArray) &&								// test negative pion
                    ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedAODMC(labelPosPion,AODMCTrackArray) 								// test positive pion
                    ) {
                  if(particle->GetPdgCode() == fPDGCodeAnalyzedMeson) {
                    fHistoMCHNMPiPlPiMiNDMInAccPt[fiCut]->Fill(particle->Pt(), tempParticleWeight ); 		// MC Eta, omega or eta prime with pi+ pi- pi0 with gamma's and e+e- in acc

                    // check relation between HNM pt and NDM, while respecting pT cutoff for NDM
                    if(!fDoLightOutput){
                      if(particleNDM->Pt() >= fNDMMinPtPossible) fHistoMCHNMInAccVsNDMPt[fiCut]->Fill(particle->Pt(),particleNDM->Pt(),tempParticleWeight);
                    }
                  }
                }
              } else if (fNDMRecoMode == 1){ // mixed mode
                // check if within PCM acceptance first
                if( kDaughter0IsPrim && kDaughter1IsPrim && kNegPionIsPrim && kPosPionIsPrim &&
                    ((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedAODMC(gamma1,AODMCTrackArray,kFALSE) &&					// test first daugther of pi0
                    ((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedAODMC(gamma2,AODMCTrackArray,kFALSE) &&					// test second daughter of pi0
                    ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedAODMC(labelNegPion,AODMCTrackArray) &&								// test negative pion
                    ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedAODMC(labelPosPion,AODMCTrackArray) 								// test positive pion
                    ) {
                  // check acceptance of clusters as well, true if one of them points into the Calo acceptance
                  if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(gamma1,AODMCTrackArray) ||
                      ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(gamma2,AODMCTrackArray) ){
                    if(particle->GetPdgCode() == fPDGCodeAnalyzedMeson){
                      fHistoMCHNMPiPlPiMiNDMInAccPt[fiCut]->Fill(particle->Pt(), weighted* tempParticleWeight ); 		// MC Eta, omega or eta prime with pi+ pi- pi0 with gamma's and e+e- in acc

                      // check relation between HNM pt and NDM, while respecting pT cutoff for NDM
                      if(!fDoLightOutput){
                        if(particleNDM->Pt() >= fNDMMinPtPossible) fHistoMCHNMInAccVsNDMPt[fiCut]->Fill(particle->Pt(),particleNDM->Pt(),weighted* tempParticleWeight);
                      }
                    }
                  }
                }
              } else if (fNDMRecoMode == 2){
                if( kDaughter0IsPrim && kDaughter1IsPrim && kNegPionIsPrim && kPosPionIsPrim &&
                    ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(gamma1,AODMCTrackArray) &&					// test first daugther of pi0
                    ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(gamma2,AODMCTrackArray) &&					// test second daughter of pi0
                    ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedAODMC(labelNegPion,AODMCTrackArray) &&								// test negative pion
                    ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedAODMC(labelPosPion,AODMCTrackArray) 								// test positive pion
                    ) {
                  if(particle->GetPdgCode() == fPDGCodeAnalyzedMeson){
                    fHistoMCHNMPiPlPiMiNDMInAccPt[fiCut]->Fill(particle->Pt(), weighted* tempParticleWeight ); 		// MC Eta pi+ pi- pi0 with gamma's and e+e- in acc

                    // check relation between HNM pt and NDM, while respecting pT cutoff for NDM
                    if(!fDoLightOutput){
                      if(particleNDM->Pt() >= fNDMMinPtPossible) fHistoMCHNMInAccVsNDMPt[fiCut]->Fill(particle->Pt(),particleNDM->Pt(),weighted* tempParticleWeight);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::CalculateMesonCandidates(AliAODConversionPhoton *vParticle){

  // Conversion Gammas

  if( fNeutralDecayParticleCandidates->GetEntries() > 0){
    for(Int_t mesonIndex=0; mesonIndex<fNeutralDecayParticleCandidates->GetEntries(); mesonIndex++){
      AliAODConversionMother *neutralDecayMeson= (AliAODConversionMother*) fNeutralDecayParticleCandidates->At(mesonIndex);
      if (neutralDecayMeson==nullptr) continue;
        if (vParticle==nullptr) continue;
        //Check for same Electron ID
        AliAODConversionMother* mesoncand = new AliAODConversionMother(neutralDecayMeson,vParticle);
        //mesoncand->SetLabels(mesonIndex,virtualParticleIndex);
        if( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(mesoncand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
          AliVTrack *negPionCandidatetmp = dynamic_cast<AliVTrack*>(fInputEvent->GetTrack(vParticle->GetTrackLabel(1)));
          AliVTrack *posPionCandidatetmp = dynamic_cast<AliVTrack*>(fInputEvent->GetTrack(vParticle->GetTrackLabel(0)));
          if(!(negPionCandidatetmp || posPionCandidatetmp)) continue;
          AliAODConversionMother NegPiontmp, PosPiontmp;
          NegPiontmp.SetPxPyPzE(negPionCandidatetmp->Px(), negPionCandidatetmp->Py(), negPionCandidatetmp->Pz(), negPionCandidatetmp->E());
          PosPiontmp.SetPxPyPzE(posPionCandidatetmp->Px(), posPionCandidatetmp->Py(), posPionCandidatetmp->Pz(), posPionCandidatetmp->E());


          // Fix Pz of pi0 candidate to match pi0 PDG mass
          AliAODConversionMother NDMtmp;
          NDMtmp.SetPxPyPzE(neutralDecayMeson->Px(), neutralDecayMeson->Py(), neutralDecayMeson->Pz(), neutralDecayMeson->Energy());
          FixPzToMatchPDGInvMassNDM(&NDMtmp);

          //Variables for Dalitz plot and Pi0 Pi+- Mass Cut
          AliKFParticle PosPionKFtmp( *posPionCandidatetmp, 211 );
          AliKFParticle NegPionKFtmp( *negPionCandidatetmp, 211 );

          TLorentzVector PosPionTLVtmp;
          TLorentzVector NegPionTLVtmp;
          TLorentzVector PosNegPionTLVtmp;

          PosPionTLVtmp.SetPxPyPzE (PosPionKFtmp.Px(), PosPionKFtmp.Py(), PosPionKFtmp.Pz(), PosPionKFtmp.E() );
          NegPionTLVtmp.SetPxPyPzE (NegPionKFtmp.Px(), NegPionKFtmp.Py(), NegPionKFtmp.Pz(), NegPionKFtmp.E() );
          PosNegPionTLVtmp = PosPionTLVtmp + NegPionTLVtmp;

          TLorentzVector NDMTLVtmp;
          TLorentzVector NDMSubTLVtmp;
          TLorentzVector PosPionNDMTLVtmp;
          TLorentzVector NegPionNDMTLVtmp;
          TLorentzVector PosPionNDMSubTLVtmp;
          TLorentzVector NegPionNDMSubTLVtmp;

          NDMTLVtmp.SetPxPyPzE( NDMtmp.Px(), NDMtmp.Py(), NDMtmp.Pz(), NDMtmp.E() ); //Fixed Pz
          NDMSubTLVtmp.SetPxPyPzE (neutralDecayMeson->Px(), neutralDecayMeson->Py(), neutralDecayMeson->Pz(), neutralDecayMeson->Energy());
          PosPionNDMTLVtmp = PosPionTLVtmp + NDMTLVtmp;
          NegPionNDMTLVtmp = NegPionTLVtmp + NDMTLVtmp;
          PosPionNDMSubTLVtmp = PosPionTLVtmp + NDMSubTLVtmp;
          NegPionNDMSubTLVtmp = NegPionTLVtmp + NDMSubTLVtmp;

          //Double_t  Mass_PiPl_NDM_FixPz         = PosPionNDMTLVtmp.M(); //Outcommented as currently not needed
          //Double_t  Mass_PiMi_NDM_FixPz         = NegPionNDMTLVtmp.M(); //Outcommented as currently not needed
          Double_t  Mass_PiPl_NDM_Sub           = PosPionNDMSubTLVtmp.M() - (NDMSubTLVtmp.M() - fPDGMassNDM);
          Double_t  Mass_PiMi_NDM_Sub           = NegPionNDMSubTLVtmp.M() - (NDMSubTLVtmp.M() - fPDGMassNDM);

          Double_t  MassCutValue_PiPlMi_NDM     = ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetMassCut_WithNDM();
          Bool_t    useMassCut_WithNDM          = ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->DoMassCut_WithNDM();
          if ((!useMassCut_WithNDM)||
                  ((useMassCut_WithNDM)
                   &&(Mass_PiPl_NDM_Sub<MassCutValue_PiPlMi_NDM)
                   &&(Mass_PiMi_NDM_Sub<MassCutValue_PiPlMi_NDM)
                   )){
            if(KinematicCut(&NegPiontmp, &PosPiontmp, neutralDecayMeson, mesoncand)){
              if(!fDoLightOutput){
                fHistoAngleHNMesonNDM[fiCut]->Fill(mesoncand->Pt(),neutralDecayMeson->Angle(mesoncand->Vect()), fWeightJetJetMC);
                fHistoAngleHNMesonPiPl[fiCut]->Fill(mesoncand->Pt(),PosPiontmp.Angle(mesoncand->Vect()), fWeightJetJetMC);
                fHistoAngleHNMesonPiMi[fiCut]->Fill(mesoncand->Pt(),NegPiontmp.Angle(mesoncand->Vect()), fWeightJetJetMC);
                fHistoAngleNDMPiMi[fiCut]->Fill(mesoncand->Pt(),NegPiontmp.Angle(neutralDecayMeson->Vect()), fWeightJetJetMC);
                fHistoAnglePiPlPiMi[fiCut]->Fill(mesoncand->Pt(),NegPiontmp.Angle(PosPiontmp.Vect()), fWeightJetJetMC);
                fHistoAnglePiPlNDM[fiCut]->Fill(mesoncand->Pt(),PosPiontmp.Angle(neutralDecayMeson->Vect()), fWeightJetJetMC);
                fHistoAngleHNMesonPiPlPiMi[fiCut]->Fill(mesoncand->Pt(),vParticle->Angle(mesoncand->Vect()), fWeightJetJetMC);
                fHistoAngleSum[fiCut]->Fill(mesoncand->Pt(),((PosPiontmp.Angle(mesoncand->Vect()))+(NegPiontmp.Angle(PosPiontmp.Vect()))+(PosPiontmp.Angle(neutralDecayMeson->Vect()))), fWeightJetJetMC);
              }

              // Subtract mass of used pi0 candidate and then add PDG mass to get to right range again
              fHistoMotherInvMassSubNDM[fiCut]->Fill(mesoncand->M()-(neutralDecayMeson->M()-fPDGMassNDM),mesoncand->Pt(), fWeightJetJetMC);

              AliAODConversionMother mesontmp(&NDMtmp,vParticle);
              fHistoMotherInvMassFixedPzNDM[fiCut]->Fill(mesontmp.M(),mesontmp.Pt(), fWeightJetJetMC);
              fHistoMotherInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);

              if(fDoMesonQA>0){
                //Dalitz All Pt
                if (enableDalitzAllPt){
                  fHistoDalitzPlotPosFixedPzNDM[fiCut]->Fill(PosNegPionTLVtmp.M(), PosPionNDMTLVtmp.M() );
                  fHistoDalitzPlotNegFixedPzNDM[fiCut]->Fill(PosNegPionTLVtmp.M(), NegPionNDMTLVtmp.M() );
                  fHistoDalitzPlotPosSubNDM[fiCut]->Fill( PosNegPionTLVtmp.M(), PosPionNDMSubTLVtmp.M() - (NDMSubTLVtmp.M() - fPDGMassNDM));
                  fHistoDalitzPlotNegSubNDM[fiCut]->Fill( PosNegPionTLVtmp.M(), NegPionNDMSubTLVtmp.M() - (NDMSubTLVtmp.M() - fPDGMassNDM));
                }
                //Dalitz Low Pt
                if (enableDalitzLowPt){
                  if ((mesoncand->Pt()>HistoDalitzPtRangeMin_LowPt)&&(mesoncand->Pt()<HistoDalitzPtRangeMax_LowPt)){
                    fHistoDalitzPlotPosFixedPzNDM_LowPt[fiCut]->Fill(PosNegPionTLVtmp.M(), PosPionNDMTLVtmp.M() );
                    fHistoDalitzPlotNegFixedPzNDM_LowPt[fiCut]->Fill(PosNegPionTLVtmp.M(), NegPionNDMTLVtmp.M() );
                    fHistoDalitzPlotPosSubNDM_LowPt[fiCut]->Fill( PosNegPionTLVtmp.M(), PosPionNDMSubTLVtmp.M() - (NDMSubTLVtmp.M() - fPDGMassNDM));
                    fHistoDalitzPlotNegSubNDM_LowPt[fiCut]->Fill( PosNegPionTLVtmp.M(), NegPionNDMSubTLVtmp.M() - (NDMSubTLVtmp.M() - fPDGMassNDM));
                  }
                }
                //Dalitz Mid Pt
                if (enableDalitzMidPt){
                  if ((mesoncand->Pt()>HistoDalitzPtRangeMin_MidPt)&&(mesoncand->Pt()<HistoDalitzPtRangeMax_MidPt)){
                    fHistoDalitzPlotPosFixedPzNDM_MidPt[fiCut]->Fill(PosNegPionTLVtmp.M(), PosPionNDMTLVtmp.M() );
                    fHistoDalitzPlotNegFixedPzNDM_MidPt[fiCut]->Fill(PosNegPionTLVtmp.M(), NegPionNDMTLVtmp.M() );
                    fHistoDalitzPlotPosSubNDM_MidPt[fiCut]->Fill( PosNegPionTLVtmp.M(), PosPionNDMSubTLVtmp.M() - (NDMSubTLVtmp.M() - fPDGMassNDM));
                    fHistoDalitzPlotNegSubNDM_MidPt[fiCut]->Fill( PosNegPionTLVtmp.M(), NegPionNDMSubTLVtmp.M() - (NDMSubTLVtmp.M() - fPDGMassNDM));
                  }
                }
                //Dalitz High Pt
                if (enableDalitzHighPt){
                  if ((mesoncand->Pt()>HistoDalitzPtRangeMin_HighPt)&&(mesoncand->Pt()<HistoDalitzPtRangeMax_HighPt)){
                    fHistoDalitzPlotPosFixedPzNDM_HighPt[fiCut]->Fill(PosNegPionTLVtmp.M(), PosPionNDMTLVtmp.M() );
                    fHistoDalitzPlotNegFixedPzNDM_HighPt[fiCut]->Fill(PosNegPionTLVtmp.M(), NegPionNDMTLVtmp.M() );
                    fHistoDalitzPlotPosSubNDM_HighPt[fiCut]->Fill( PosNegPionTLVtmp.M(), PosPionNDMSubTLVtmp.M() - (NDMSubTLVtmp.M() - fPDGMassNDM));
                    fHistoDalitzPlotNegSubNDM_HighPt[fiCut]->Fill( PosNegPionTLVtmp.M(), NegPionNDMSubTLVtmp.M() - (NDMSubTLVtmp.M() - fPDGMassNDM));
                  }
                }
              } //end fDoMesonQA>0
              if(fMCEvent){
                if(fInputEvent->IsA()==AliESDEvent::Class())
                  ProcessTrueMesonCandidates(mesoncand,neutralDecayMeson,vParticle);
                if(fInputEvent->IsA()==AliAODEvent::Class())
                  ProcessTrueMesonCandidatesAOD(mesoncand,neutralDecayMeson,vParticle);
              }
            }else{ //else KinematicCut
              if(!fDoLightOutput){
                fHistoMotherInvMassPtRejectedKinematic[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
              }
            } //end KinematicCut
          }
        } //end MesonIsSelected
       delete mesoncand;
       mesoncand=0x0;
    }
  }

}

//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::CalculateBackground(Int_t mode = 4){

  /*
  * Mode 1 => pi+ and pi- from same event
  * Mode 2 => pi+ and pi0 from same event
  * Mode 3 => pi- and pi0 from same event
  * Mode 4 => no pions from same event (default)
  * Mode 5 => Ligesign mixing
  * Mode 6 => Sideband
  * Mode 7 => Rotation around pi0
  */
  // Get multiplicity and zbin from fBGHandler
  Int_t zbin = fBGHandlerPiMi[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
  Int_t mbin = 0;

  // Multiplicity can be determined either by number of cluster candidates or track mulitiplicity
  if (((AliConversionMesonCuts *)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()) {
    mbin = fBGHandlerPiMi[fiCut]->GetMultiplicityBinIndex(fNumberOfESDTracks);
  } else {
    if (fNDMRecoMode < 2)
      mbin = fBGHandlerPiMi[fiCut]->GetMultiplicityBinIndex(fGoodConvGammas->GetEntries());
    else
      mbin = fBGHandlerPiMi[fiCut]->GetMultiplicityBinIndex(fClusterCandidates->GetEntries());
  }

  AliGammaConversionAODBGHandler::GammaConversionVertex *bgEventVertexPl = nullptr;
  AliGammaConversionAODBGHandler::GammaConversionVertex *bgEventVertexMi = nullptr;

  // Get N of Pi0 according to chosen mix mode
  Int_t NNDMCandidates = 0;
  if ((((AliConversionMesonCuts *)fMesonCutArray->At(fiCut))->UseSidebandMixing()) || (((AliConversionMesonCuts *)fMesonCutArray->At(fiCut))->UseSidebandMixingBothSides())) {
    NNDMCandidates = fNeutralDecayParticleSidebandCandidates->GetEntries();
  } else if(mode == 7){ // Swapp (Rotation) method entries
    NNDMCandidates = fNeutralDecayParticleSwappCandidates->GetEntries();
  } else {
    NNDMCandidates = fNeutralDecayParticleCandidates->GetEntries();
  }
  //
  // ─── LOOP OVER ALL NDM FROM CURRENT EVENT ───────────────────────────────────────
  //
  for (Int_t iCurrentNDM = 0; iCurrentNDM < NNDMCandidates; iCurrentNDM++) {
    AliAODConversionMother *EventNDMGoodMeson;
    if ((((AliConversionMesonCuts *)fMesonCutArray->At(fiCut))->UseSidebandMixing()) || (((AliConversionMesonCuts *)fMesonCutArray->At(fiCut))->UseSidebandMixingBothSides())) {
      EventNDMGoodMeson = (AliAODConversionMother *)(fNeutralDecayParticleSidebandCandidates->At(iCurrentNDM));
    } else if(mode == 7){  // Swapp (Rotation) method entries
      EventNDMGoodMeson = (AliAODConversionMother *)(fNeutralDecayParticleSwappCandidates->At(iCurrentNDM));
    } else {
      EventNDMGoodMeson = (AliAODConversionMother *)(fNeutralDecayParticleCandidates->At(iCurrentNDM));
    }
    //
    // ─── Pi+ PI- from same event───
    //
    if(mode==1){
      // Begin loop over BG events for Pi+
      for (Int_t nEventsInBGPl = 0; nEventsInBGPl < fBGHandlerPiPl[fiCut]->GetNBGEvents(); nEventsInBGPl++) {

        // Store all Pi+ of current event in right binning in vector
        AliGammaConversionMotherAODVector *EventPiPlMeson = fBGHandlerPiPl[fiCut]->GetBGGoodMesons(zbin, mbin, nEventsInBGPl);
        if(!EventPiPlMeson) continue;

        // Begin loop over BG events for Pi-
        for (Int_t nEventsInBGMi = 0; nEventsInBGMi < fBGHandlerPiMi[fiCut]->GetNBGEvents(); nEventsInBGMi++) {
          AliGammaConversionMotherAODVector *EventPiMiMeson = fBGHandlerPiMi[fiCut]->GetBGGoodMesons(zbin, mbin, nEventsInBGMi);

          // If one of the events isn't found skip to next one
          if(!EventPiMiMeson) continue;

          // If events are unequal, skip:
          if (nEventsInBGMi != nEventsInBGPl) continue;


          // Determine Background event vertex
          if (fMoveParticleAccordingToVertex == kTRUE) {
            bgEventVertexPl = fBGHandlerPiPl[fiCut]->GetBGEventVertex(zbin, mbin, nEventsInBGPl);
            bgEventVertexMi = fBGHandlerPiMi[fiCut]->GetBGEventVertex(zbin, mbin, nEventsInBGMi);
          }
          // Loop over all Pi+
          for (UInt_t iCurrentPiPl = 0; iCurrentPiPl < EventPiPlMeson->size(); iCurrentPiPl++) {
            AliAODConversionMother EventPiPlGoodMeson = (AliAODConversionMother)(*(EventPiPlMeson->at(iCurrentPiPl)));

            // Move Vertex
            if (fMoveParticleAccordingToVertex == kTRUE) {
              MoveParticleAccordingToVertex(&EventPiPlGoodMeson, bgEventVertexPl);
            }

            for (UInt_t iCurrentPiMi = 0; iCurrentPiMi < EventPiMiMeson->size(); iCurrentPiMi++) {
              AliAODConversionMother EventPiMiGoodMeson = (AliAODConversionMother)(*(EventPiMiMeson->at(iCurrentPiMi)));

              // Move Vertex
              if (fMoveParticleAccordingToVertex == kTRUE) {
                MoveParticleAccordingToVertex(&EventPiMiGoodMeson, bgEventVertexMi);
              }

              // create momentum vector for all three particles
              TLorentzVector vec4PiPlus, vec4PiMinus, vec4NDM;
              vec4PiPlus.SetPxPyPzE(EventPiPlGoodMeson.Px(),EventPiPlGoodMeson.Py(),EventPiPlGoodMeson.Pz(),EventPiPlGoodMeson.Energy());
              vec4PiMinus.SetPxPyPzE(EventPiMiGoodMeson.Px(),EventPiMiGoodMeson.Py(),EventPiMiGoodMeson.Pz(),EventPiMiGoodMeson.Energy());
              vec4NDM.SetPxPyPzE(EventNDMGoodMeson->Px(),EventNDMGoodMeson->Py(),EventNDMGoodMeson->Pz(),EventNDMGoodMeson->Energy());

              // Mass cut (pi+pi-)
              if (((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->DoMassCut()) {
                TLorentzVector vec4PiPlusPiMinus = vec4PiPlus + vec4PiMinus;
                if (vec4PiPlusPiMinus.M() >= ((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->GetMassCut()) {
                  continue;
                }
              }
              // Mass cut (pi0pi+-)
              if (((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->DoMassCut_WithNDM()) {
                TLorentzVector vec4PiPlusPiZero = vec4PiPlus + vec4NDM;
                TLorentzVector vec4PiMinusPiZero = vec4PiMinus + vec4NDM;
                Double_t  Mass_PiPlus_PiZero_Sub           = vec4PiPlusPiZero.M() - (vec4NDM.M() - fPDGMassNDM);
                Double_t  Mass_PiMinus_PiZero_Sub           = vec4PiMinusPiZero.M() - (vec4NDM.M() - fPDGMassNDM);
                if ((Mass_PiPlus_PiZero_Sub >= ((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->GetMassCut_WithNDM())||
                     (Mass_PiMinus_PiZero_Sub >= ((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->GetMassCut_WithNDM())
                     ) {
                  continue;
                }
              }

              // Create Pi+Pi- pair (only at this stage after cuts were applied to save time, before only vectors)
              AliAODConversionMother backPiPlPiMiCandidate(&EventPiPlGoodMeson, &EventPiMiGoodMeson);
              AliAODConversionMother PiPlPiMiNDMBackgroundCandidate(&backPiPlPiMiCandidate, EventNDMGoodMeson);

              // Check if candidate survives meson cut
              if (((AliConversionMesonCuts *)fMesonCutArray->At(fiCut))->MesonIsSelected(&PiPlPiMiNDMBackgroundCandidate, kFALSE, ((AliConvEventCuts *)fEventCutArray->At(fiCut))->GetEtaShift())) {

                // Check if candidate survives kinematic cut
                if (KinematicCut(&EventPiMiGoodMeson, &EventPiPlGoodMeson, EventNDMGoodMeson, &PiPlPiMiNDMBackgroundCandidate)) {
                  // Create temporary mesons to be able to fix pz
                  AliAODConversionMother NDMtmp;
                  NDMtmp.SetPxPyPzE(EventNDMGoodMeson->Px(), EventNDMGoodMeson->Py(), EventNDMGoodMeson->Pz(), EventNDMGoodMeson->Energy());
                  FixPzToMatchPDGInvMassNDM(&NDMtmp);
                  AliAODConversionMother PiMiNDMtmp(&EventPiMiGoodMeson, &NDMtmp);
                  AliAODConversionMother PiPlPiMiNDMtmp(&EventPiPlGoodMeson, &PiMiNDMtmp); // Must be two separate lines since second instance depends on first and execution order is not guaranteed

                  // Pi+ and Pi- don't come from the same event (but different than pi0 event)
                  // Fill histograms
                  fHistoBackInvMassPt[fiCut]->Fill(PiPlPiMiNDMBackgroundCandidate.M(), PiPlPiMiNDMBackgroundCandidate.Pt(), fWeightJetJetMC);
                  fHistoBackInvMassPtSubNDM[fiCut]->Fill(PiPlPiMiNDMBackgroundCandidate.M() - (EventNDMGoodMeson->M() - fPDGMassNDM), PiPlPiMiNDMBackgroundCandidate.Pt(), fWeightJetJetMC);
                  fHistoBackInvMassPtFixedPzNDM[fiCut]->Fill(PiPlPiMiNDMtmp.M(), PiPlPiMiNDMtmp.Pt(), fWeightJetJetMC);
                }
              }
            } // end pi- loop
          }   // end pi+ loop
        }     // end loop over all pi- event
      }       // end loop over pi+ events
    //
    // ─── NO PIONS FROM SAME EVENT ───
    //
    } else if(mode==4){
      // Begin loop over BG events for Pi+
      for (Int_t nEventsInBGPl = 0; nEventsInBGPl < fBGHandlerPiPl[fiCut]->GetNBGEvents(); nEventsInBGPl++) {

        // Store all Pi+ of current event in right binning in vector
        AliGammaConversionMotherAODVector *EventPiPlMeson = fBGHandlerPiPl[fiCut]->GetBGGoodMesons(zbin, mbin, nEventsInBGPl);
        if(!EventPiPlMeson) continue;

        // Begin loop over BG events for Pi-
        for (Int_t nEventsInBGMi = 0; nEventsInBGMi < fBGHandlerPiMi[fiCut]->GetNBGEvents(); nEventsInBGMi++) {
          AliGammaConversionMotherAODVector *EventPiMiMeson = fBGHandlerPiMi[fiCut]->GetBGGoodMesons(zbin, mbin, nEventsInBGMi);

          // If one of the events isn't found skip to next one
          if(!EventPiMiMeson) continue;

          // If events are equal, skip:
          if (nEventsInBGMi == nEventsInBGPl) continue;

          // Determine Background event vertex
          if (fMoveParticleAccordingToVertex == kTRUE) {
            bgEventVertexPl = fBGHandlerPiPl[fiCut]->GetBGEventVertex(zbin, mbin, nEventsInBGPl);
            bgEventVertexMi = fBGHandlerPiMi[fiCut]->GetBGEventVertex(zbin, mbin, nEventsInBGMi);
          }
          // Loop over all Pi+
          for (UInt_t iCurrentPiPl = 0; iCurrentPiPl < EventPiPlMeson->size(); iCurrentPiPl++) {
            AliAODConversionMother EventPiPlGoodMeson = (AliAODConversionMother)(*(EventPiPlMeson->at(iCurrentPiPl)));

            // Move Vertex
            if (fMoveParticleAccordingToVertex == kTRUE) {
              MoveParticleAccordingToVertex(&EventPiPlGoodMeson, bgEventVertexPl);
            }


            for (UInt_t iCurrentPiMi = 0; iCurrentPiMi < EventPiMiMeson->size(); iCurrentPiMi++) {
              AliAODConversionMother EventPiMiGoodMeson = (AliAODConversionMother)(*(EventPiMiMeson->at(iCurrentPiMi)));

              // Move Vertex
              if (fMoveParticleAccordingToVertex == kTRUE) {
                MoveParticleAccordingToVertex(&EventPiMiGoodMeson, bgEventVertexMi);
              }

              // create momentum vector for all three particles
              TLorentzVector vec4PiPlus, vec4PiMinus, vec4NDM;
              vec4PiPlus.SetPxPyPzE(EventPiPlGoodMeson.Px(),EventPiPlGoodMeson.Py(),EventPiPlGoodMeson.Pz(),EventPiPlGoodMeson.Energy());
              vec4PiMinus.SetPxPyPzE(EventPiMiGoodMeson.Px(),EventPiMiGoodMeson.Py(),EventPiMiGoodMeson.Pz(),EventPiMiGoodMeson.Energy());
              vec4NDM.SetPxPyPzE(EventNDMGoodMeson->Px(),EventNDMGoodMeson->Py(),EventNDMGoodMeson->Pz(),EventNDMGoodMeson->Energy());

              // Mass cut (pi+pi-)
              if (((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->DoMassCut()) {
                TLorentzVector vec4PiPlusPiMinus = vec4PiPlus + vec4PiMinus;
                if (vec4PiPlusPiMinus.M() >= ((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->GetMassCut()) {
                  continue;
                }
              }
              // Mass cut (pi0pi+-)
              if (((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->DoMassCut_WithNDM()) {
                TLorentzVector vec4PiPlusPiZero = vec4PiPlus + vec4NDM;
                TLorentzVector vec4PiMinusPiZero = vec4PiMinus + vec4NDM;
                Double_t  Mass_PiPlus_PiZero_Sub           = vec4PiPlusPiZero.M() - (vec4NDM.M() - fPDGMassNDM);
                Double_t  Mass_PiMinus_PiZero_Sub           = vec4PiMinusPiZero.M() - (vec4NDM.M() - fPDGMassNDM);
                if ((Mass_PiPlus_PiZero_Sub >= ((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->GetMassCut_WithNDM())||
                     (Mass_PiMinus_PiZero_Sub >= ((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->GetMassCut_WithNDM())
                     ) {
                  continue;
                }
              }

              // Create Pi+Pi- pair (only at this stage after cuts were applied to save time, before only vectors)
              AliAODConversionMother backPiPlPiMiCandidate(&EventPiPlGoodMeson, &EventPiMiGoodMeson);
              AliAODConversionMother PiPlPiMiNDMBackgroundCandidate(&backPiPlPiMiCandidate, EventNDMGoodMeson);

              // Check if candidate survives meson cut
              if (((AliConversionMesonCuts *)fMesonCutArray->At(fiCut))->MesonIsSelected(&PiPlPiMiNDMBackgroundCandidate, kFALSE, ((AliConvEventCuts *)fEventCutArray->At(fiCut))->GetEtaShift())) {

                // Check if candidate survives kinematic cut
                if (KinematicCut(&EventPiMiGoodMeson, &EventPiPlGoodMeson, EventNDMGoodMeson, &PiPlPiMiNDMBackgroundCandidate)) {
                  // Create temporary mesons to be able to fix pz
                  AliAODConversionMother NDMtmp;
                  NDMtmp.SetPxPyPzE(EventNDMGoodMeson->Px(), EventNDMGoodMeson->Py(), EventNDMGoodMeson->Pz(), EventNDMGoodMeson->Energy());
                  FixPzToMatchPDGInvMassNDM(&NDMtmp);
                  AliAODConversionMother PiMiNDMtmp(&EventPiMiGoodMeson, &NDMtmp);
                  AliAODConversionMother PiPlPiMiNDMtmp(&EventPiPlGoodMeson, &PiMiNDMtmp); // Must be two separate lines since second instance depends on first and execution order is not guaranteed

                  // Pi+ and Pi- don't come from the same event (but different than pi0 event)
                  // Fill histograms
                  fHistoBackInvMassPt[fiCut]->Fill(PiPlPiMiNDMBackgroundCandidate.M(), PiPlPiMiNDMBackgroundCandidate.Pt(), fWeightJetJetMC);
                  fHistoBackInvMassPtSubNDM[fiCut]->Fill(PiPlPiMiNDMBackgroundCandidate.M() - (EventNDMGoodMeson->M() - fPDGMassNDM), PiPlPiMiNDMBackgroundCandidate.Pt(), fWeightJetJetMC);
                  fHistoBackInvMassPtFixedPzNDM[fiCut]->Fill(PiPlPiMiNDMtmp.M(), PiPlPiMiNDMtmp.Pt(), fWeightJetJetMC);
                }
              }
            } // end pi- loop
          }   // end pi+ loop
        }     // end loop over all pi- event
      }       // end loop over pi+ events

    //
    // ─── PIPL AND PIZERO FROM SAME EVENT ─────────────────────────────
    //
    } else if(mode==2){
      // Loop over PiPl from current event
      for (Int_t iCurrentPiPl = 0; iCurrentPiPl < fPosPionCandidates->GetEntries(); iCurrentPiPl++) {
         AliAODConversionMother EventPiPlGoodMeson = *(AliAODConversionMother *)(fPosPionCandidates->At(iCurrentPiPl));

        // Begin loop over BG events for Pi-
        for (Int_t nEventsInBGMi = 0; nEventsInBGMi < fBGHandlerPiMi[fiCut]->GetNBGEvents(); nEventsInBGMi++) {
          AliGammaConversionMotherAODVector *EventPiMiMeson = fBGHandlerPiMi[fiCut]->GetBGGoodMesons(zbin, mbin, nEventsInBGMi);

          // If one of the events isn't found skip to next one
          if(!EventPiMiMeson) continue;

          // Determine Background event vertex
          if (fMoveParticleAccordingToVertex == kTRUE) {
            bgEventVertexMi = fBGHandlerPiMi[fiCut]->GetBGEventVertex(zbin, mbin, nEventsInBGMi);
          }

          for (UInt_t iCurrentPiMi = 0; iCurrentPiMi < EventPiMiMeson->size(); iCurrentPiMi++) {
            AliAODConversionMother EventPiMiGoodMeson = (AliAODConversionMother)(*(EventPiMiMeson->at(iCurrentPiMi)));

            // Move Vertex
            if (fMoveParticleAccordingToVertex == kTRUE) {
              MoveParticleAccordingToVertex(&EventPiMiGoodMeson, bgEventVertexMi);
            }

            // create momentum vector for all three particles
            TLorentzVector vec4PiPlus, vec4PiMinus, vec4NDM;
            vec4PiPlus.SetPxPyPzE(EventPiPlGoodMeson.Px(),EventPiPlGoodMeson.Py(),EventPiPlGoodMeson.Pz(),EventPiPlGoodMeson.Energy());
            vec4PiMinus.SetPxPyPzE(EventPiMiGoodMeson.Px(),EventPiMiGoodMeson.Py(),EventPiMiGoodMeson.Pz(),EventPiMiGoodMeson.Energy());
            vec4NDM.SetPxPyPzE(EventNDMGoodMeson->Px(),EventNDMGoodMeson->Py(),EventNDMGoodMeson->Pz(),EventNDMGoodMeson->Energy());

            // Mass cut (pi+pi-)
            if (((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->DoMassCut()) {
              TLorentzVector vec4PiPlusPiMinus = vec4PiPlus + vec4PiMinus;
              if (vec4PiPlusPiMinus.M() >= ((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->GetMassCut()) {
                continue;
              }
            }
            // Mass cut (pi0pi+-)
            if (((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->DoMassCut_WithNDM()) {
              TLorentzVector vec4PiPlusPiZero = vec4PiPlus + vec4NDM;
              TLorentzVector vec4PiMinusPiZero = vec4PiMinus + vec4NDM;
              Double_t  Mass_PiPlus_PiZero_Sub           = vec4PiPlusPiZero.M() - (vec4NDM.M() - fPDGMassNDM);
              Double_t  Mass_PiMinus_PiZero_Sub           = vec4PiMinusPiZero.M() - (vec4NDM.M() - fPDGMassNDM);
              if ((Mass_PiPlus_PiZero_Sub >= ((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->GetMassCut_WithNDM())||
                   (Mass_PiMinus_PiZero_Sub >= ((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->GetMassCut_WithNDM())
                   ) {
                continue;
              }
            }

            // Create Pi+Pi- pair (only at this stage after cuts were applied to save time, before only vectors)
            AliAODConversionMother backPiPlPiMiCandidate(&EventPiPlGoodMeson, &EventPiMiGoodMeson);
            AliAODConversionMother PiPlPiMiNDMBackgroundCandidate(&backPiPlPiMiCandidate, EventNDMGoodMeson);

            // Check if candidate survives meson cut
            if (((AliConversionMesonCuts *)fMesonCutArray->At(fiCut))->MesonIsSelected(&PiPlPiMiNDMBackgroundCandidate, kFALSE, ((AliConvEventCuts *)fEventCutArray->At(fiCut))->GetEtaShift())) {

              // Check if candidate survives kinematic cut
              if (KinematicCut(&EventPiMiGoodMeson, &EventPiPlGoodMeson, EventNDMGoodMeson, &PiPlPiMiNDMBackgroundCandidate)) {
                // Create temporary mesons to be able to fix pz
                AliAODConversionMother NDMtmp;
                NDMtmp.SetPxPyPzE(EventNDMGoodMeson->Px(), EventNDMGoodMeson->Py(), EventNDMGoodMeson->Pz(), EventNDMGoodMeson->Energy());
                FixPzToMatchPDGInvMassNDM(&NDMtmp);
                AliAODConversionMother PiMiNDMtmp(&EventPiMiGoodMeson, &NDMtmp);
                AliAODConversionMother PiPlPiMiNDMtmp(&EventPiPlGoodMeson, &PiMiNDMtmp); // Must be two separate lines since second instance depends on first and execution order is not guaranteed

                // Pi+ and Pi- don't come from the same event (but different than pi0 event)
                // Fill histograms
                fHistoBackInvMassPt[fiCut]->Fill(PiPlPiMiNDMBackgroundCandidate.M(), PiPlPiMiNDMBackgroundCandidate.Pt(), fWeightJetJetMC);
                fHistoBackInvMassPtSubNDM[fiCut]->Fill(PiPlPiMiNDMBackgroundCandidate.M() - (EventNDMGoodMeson->M() - fPDGMassNDM), PiPlPiMiNDMBackgroundCandidate.Pt(), fWeightJetJetMC);
                fHistoBackInvMassPtFixedPzNDM[fiCut]->Fill(PiPlPiMiNDMtmp.M(), PiPlPiMiNDMtmp.Pt(), fWeightJetJetMC);
              }
            }
          } // end pi- loop

        }     // end loop over all pi- event

      }
    } else if(mode==3){
      // Loop over PiMi from current event
      for (Int_t iCurrentPiMi = 0; iCurrentPiMi < fNegPionCandidates->GetEntries(); iCurrentPiMi++) {
         AliAODConversionMother EventPiMiGoodMeson = *(AliAODConversionMother *)(fNegPionCandidates->At(iCurrentPiMi));

        // Begin loop over BG events for Pi+
        for (Int_t nEventsInBGPl = 0; nEventsInBGPl < fBGHandlerPiPl[fiCut]->GetNBGEvents(); nEventsInBGPl++) {
          AliGammaConversionMotherAODVector *EventPiPlMeson = fBGHandlerPiPl[fiCut]->GetBGGoodMesons(zbin, mbin, nEventsInBGPl);

          // If one of the events isn't found skip to next one
          if(!EventPiPlMeson) continue;

          // Determine Background event vertex
          if (fMoveParticleAccordingToVertex == kTRUE) {
            bgEventVertexPl = fBGHandlerPiPl[fiCut]->GetBGEventVertex(zbin, mbin, nEventsInBGPl);
          }

          for (UInt_t iCurrentPiPl = 0; iCurrentPiPl < EventPiPlMeson->size(); iCurrentPiPl++) {
            AliAODConversionMother EventPiPlGoodMeson = (AliAODConversionMother)(*(EventPiPlMeson->at(iCurrentPiPl)));

            // Move Vertex
            if (fMoveParticleAccordingToVertex == kTRUE) {
              MoveParticleAccordingToVertex(&EventPiPlGoodMeson, bgEventVertexPl);
            }

            // create momentum vector for all three particles
            TLorentzVector vec4PiPlus, vec4PiMinus, vec4NDM;
            vec4PiPlus.SetPxPyPzE(EventPiPlGoodMeson.Px(),EventPiPlGoodMeson.Py(),EventPiPlGoodMeson.Pz(),EventPiPlGoodMeson.Energy());
            vec4PiMinus.SetPxPyPzE(EventPiMiGoodMeson.Px(),EventPiMiGoodMeson.Py(),EventPiMiGoodMeson.Pz(),EventPiMiGoodMeson.Energy());
            vec4NDM.SetPxPyPzE(EventNDMGoodMeson->Px(),EventNDMGoodMeson->Py(),EventNDMGoodMeson->Pz(),EventNDMGoodMeson->Energy());

            // Mass cut (pi+pi-)
            if (((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->DoMassCut()) {
              TLorentzVector vec4PiPlusPiMinus = vec4PiPlus + vec4PiMinus;
              if (vec4PiPlusPiMinus.M() >= ((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->GetMassCut()) {
                continue;
              }
            }
            // Mass cut (pi0pi+-)
            if (((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->DoMassCut_WithNDM()) {
              TLorentzVector vec4PiPlusPiZero = vec4PiPlus + vec4NDM;
              TLorentzVector vec4PiMinusPiZero = vec4PiMinus + vec4NDM;
              Double_t  Mass_PiPlus_PiZero_Sub           = vec4PiPlusPiZero.M() - (vec4NDM.M() - fPDGMassNDM);
              Double_t  Mass_PiMinus_PiZero_Sub           = vec4PiMinusPiZero.M() - (vec4NDM.M() - fPDGMassNDM);
              if ((Mass_PiPlus_PiZero_Sub >= ((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->GetMassCut_WithNDM())||
                   (Mass_PiMinus_PiZero_Sub >= ((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->GetMassCut_WithNDM())
                   ) {
                continue;
              }
            }

            // Create Pi+Pi- pair (only at this stage after cuts were applied to save time, before only vectors)
            AliAODConversionMother backPiPlPiMiCandidate(&EventPiPlGoodMeson, &EventPiMiGoodMeson);
            AliAODConversionMother PiPlPiMiNDMBackgroundCandidate(&backPiPlPiMiCandidate, EventNDMGoodMeson);

            // Check if candidate survives meson cut
            if (((AliConversionMesonCuts *)fMesonCutArray->At(fiCut))->MesonIsSelected(&PiPlPiMiNDMBackgroundCandidate, kFALSE, ((AliConvEventCuts *)fEventCutArray->At(fiCut))->GetEtaShift())) {

              // Check if candidate survives kinematic cut
              if (KinematicCut(&EventPiMiGoodMeson, &EventPiPlGoodMeson, EventNDMGoodMeson, &PiPlPiMiNDMBackgroundCandidate)) {
                // Create temporary mesons to be able to fix pz
                AliAODConversionMother NDMtmp;
                NDMtmp.SetPxPyPzE(EventNDMGoodMeson->Px(), EventNDMGoodMeson->Py(), EventNDMGoodMeson->Pz(), EventNDMGoodMeson->Energy());
                FixPzToMatchPDGInvMassNDM(&NDMtmp);
                AliAODConversionMother PiMiNDMtmp(&EventPiMiGoodMeson, &NDMtmp);
                AliAODConversionMother PiPlPiMiNDMtmp(&EventPiPlGoodMeson, &PiMiNDMtmp); // Must be two separate lines since second instance depends on first and execution order is not guaranteed

                // Pi+ and Pi- don't come from the same event (but different than pi0 event)
                // Fill histograms
                fHistoBackInvMassPt[fiCut]->Fill(PiPlPiMiNDMBackgroundCandidate.M(), PiPlPiMiNDMBackgroundCandidate.Pt(), fWeightJetJetMC);
                fHistoBackInvMassPtSubNDM[fiCut]->Fill(PiPlPiMiNDMBackgroundCandidate.M() - (EventNDMGoodMeson->M() - fPDGMassNDM), PiPlPiMiNDMBackgroundCandidate.Pt(), fWeightJetJetMC);
                fHistoBackInvMassPtFixedPzNDM[fiCut]->Fill(PiPlPiMiNDMtmp.M(), PiPlPiMiNDMtmp.Pt(), fWeightJetJetMC);
              }
            }
          } // end pi- loop

        }     // end loop over all pi- event

      }

    //
    // ─── LIKESIGN MIXING ─────────────────────────────────────────────
    //
    } else if (mode == 5){
      // Loops for Pi0Pi+Pi+ LikeSign mixing
      for (Int_t iCurrentPiPl = 0; iCurrentPiPl < fPosPionCandidates->GetEntries(); iCurrentPiPl++) {

        AliAODConversionMother EventPiPlGoodMeson = *(AliAODConversionMother *)(fPosPionCandidates->At(iCurrentPiPl));

        for (Int_t iCurrentPiPl2 = iCurrentPiPl; iCurrentPiPl2 < fPosPionCandidates->GetEntries(); iCurrentPiPl2++) {

          if (iCurrentPiPl == iCurrentPiPl2) continue;
            AliAODConversionMother EventPiPlGoodMeson2 = *(AliAODConversionMother *)(fPosPionCandidates->At(iCurrentPiPl2));

            // Mass cut on pi+pi+
            if (((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->DoMassCut()) {
              AliAODConversionMother backPiPlPiPlCandidate(&EventPiPlGoodMeson, &EventPiPlGoodMeson2);
              if (backPiPlPiPlCandidate.M() >= ((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->GetMassCut()) {
                continue;
              }
            }


            // Combine Pi+ and Pi0
            AliAODConversionMother PiPlNDMBackgroundCandidate(&EventPiPlGoodMeson, EventNDMGoodMeson);

            // Mass cut (pi0pi+-)
            if (((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->DoMassCut_WithNDM()) {
              AliAODConversionMother PiPlNDMBackgroundCandidate2(&EventPiPlGoodMeson2, EventNDMGoodMeson);
              Double_t  Mass_PiPlus_PiZero_Sub           = PiPlNDMBackgroundCandidate.M() - (EventNDMGoodMeson->M() - fPDGMassNDM);
              Double_t  Mass_PiPlus2_PiZero_Sub           = PiPlNDMBackgroundCandidate2.M() - (EventNDMGoodMeson->M() - fPDGMassNDM);
              if ((Mass_PiPlus_PiZero_Sub >= ((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->GetMassCut_WithNDM())||
                   (Mass_PiPlus2_PiZero_Sub >= ((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->GetMassCut_WithNDM())
                   ) {
                continue;
              }
            }

            // Create (final) Candidate
            AliAODConversionMother PiPlPiPlNDMBackgroundCandidate(&PiPlNDMBackgroundCandidate, &EventPiPlGoodMeson2);

            // Check if candidate survives meson cut
            if (((AliConversionMesonCuts *)fMesonCutArray->At(fiCut))->MesonIsSelected(&PiPlNDMBackgroundCandidate, kFALSE, ((AliConvEventCuts *)fEventCutArray->At(fiCut))->GetEtaShift())) {

              // Check if candidate survives kinematic cut
              if (KinematicCut(&EventPiPlGoodMeson, &EventPiPlGoodMeson2, EventNDMGoodMeson, &PiPlPiPlNDMBackgroundCandidate)) {

                // Create temporary mesons to be able to fix pz
                AliAODConversionMother NDMtmp;
                NDMtmp.SetPxPyPzE(EventNDMGoodMeson->Px(), EventNDMGoodMeson->Py(), EventNDMGoodMeson->Pz(), EventNDMGoodMeson->Energy());
                FixPzToMatchPDGInvMassNDM(&NDMtmp);
                AliAODConversionMother PiPlNDMtmp(&EventPiPlGoodMeson, &NDMtmp);
                AliAODConversionMother PiPlPiPlNDMtmp(&EventPiPlGoodMeson2, &PiPlNDMtmp); // Must be two separate lines since second instance depends on first and execution order is not guaranteed

                // Fill histograms (likesign)
                fHistoMotherLikeSignBackInvMassPt[fiCut]->Fill(PiPlPiPlNDMBackgroundCandidate.M(), PiPlPiPlNDMBackgroundCandidate.Pt(), fWeightJetJetMC);
                fHistoMotherLikeSignBackInvMassSubNDMPt[fiCut]->Fill(PiPlPiPlNDMBackgroundCandidate.M() - (EventNDMGoodMeson->M() - fPDGMassNDM), PiPlPiPlNDMBackgroundCandidate.Pt(), fWeightJetJetMC);
                fHistoMotherLikeSignBackInvMassFixedPzNDMPt[fiCut]->Fill(PiPlPiPlNDMtmp.M(), PiPlPiPlNDMtmp.Pt(), fWeightJetJetMC);
              }
            }

        } // end of iCurrentPiPl2
      }   // end of iCurrenPiPl

      // Loops for Pi0Pi-Pi- LikeSign mixing
      for (Int_t iCurrentPiMi = 0; iCurrentPiMi < fNegPionCandidates->GetEntries(); iCurrentPiMi++) {

        AliAODConversionMother EventPiMiGoodMeson = *(AliAODConversionMother *)(fNegPionCandidates->At(iCurrentPiMi));

        for (Int_t iCurrentPiMi2 = iCurrentPiMi; iCurrentPiMi2 < fNegPionCandidates->GetEntries(); iCurrentPiMi2++){

          if (iCurrentPiMi == iCurrentPiMi2) continue;
            AliAODConversionMother EventPiMiGoodMeson2 = *(AliAODConversionMother *)(fNegPionCandidates->At(iCurrentPiMi2));

            // Mass cut on pi-pi-
            if (((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->DoMassCut()) {
              AliAODConversionMother backPiMiPiMiCandidate(&EventPiMiGoodMeson, &EventPiMiGoodMeson2);
              if (backPiMiPiMiCandidate.M() >= ((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->GetMassCut()) {
                continue;
              }
            }

            // Combine Pi- and Pi0
            AliAODConversionMother PiMiNDMBackgroundCandidate(&EventPiMiGoodMeson, EventNDMGoodMeson);

            // Mass cut (pi0pi+-)
            if (((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->DoMassCut_WithNDM()) {
              AliAODConversionMother PiMiNDMBackgroundCandidate2(&EventPiMiGoodMeson2, EventNDMGoodMeson);
              Double_t  Mass_PiMinus_PiZero_Sub           = PiMiNDMBackgroundCandidate.M() - (EventNDMGoodMeson->M() - fPDGMassNDM);
              Double_t  Mass_PiMinus2_PiZero_Sub           = PiMiNDMBackgroundCandidate2.M() - (EventNDMGoodMeson->M() - fPDGMassNDM);
              if ((Mass_PiMinus_PiZero_Sub >= ((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->GetMassCut_WithNDM())||
                   (Mass_PiMinus2_PiZero_Sub >= ((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->GetMassCut_WithNDM())
                   ) {
                continue;
              }
            }

            // Create (final) Candidate
            AliAODConversionMother PiMiPiMiNDMBackgroundCandidate(&PiMiNDMBackgroundCandidate, &EventPiMiGoodMeson2);

            // Check if candidate survives meson cut
            if (((AliConversionMesonCuts *)fMesonCutArray->At(fiCut))->MesonIsSelected(&PiMiPiMiNDMBackgroundCandidate, kFALSE, ((AliConvEventCuts *)fEventCutArray->At(fiCut))->GetEtaShift())) {

              // Check if candidate survives kinematic cut
              if (KinematicCut(&EventPiMiGoodMeson, &EventPiMiGoodMeson2, EventNDMGoodMeson, &PiMiPiMiNDMBackgroundCandidate)) {

                // Create temporary mesons to be able to fix pz
                AliAODConversionMother NDMtmp;
                NDMtmp.SetPxPyPzE(EventNDMGoodMeson->Px(), EventNDMGoodMeson->Py(), EventNDMGoodMeson->Pz(), EventNDMGoodMeson->Energy());
                FixPzToMatchPDGInvMassNDM(&NDMtmp);
                AliAODConversionMother PiMiNDMtmp(&EventPiMiGoodMeson, &NDMtmp);
                AliAODConversionMother PiMiPiMiNDMtmp(&EventPiMiGoodMeson2, &PiMiNDMtmp);  // Must be two separate lines since second instance depends on first and execution order is not guaranteed

                // Fill histograms (likesign)
                fHistoMotherLikeSignBackInvMassPt[fiCut]->Fill(PiMiPiMiNDMBackgroundCandidate.M(), PiMiPiMiNDMBackgroundCandidate.Pt(), fWeightJetJetMC);
                fHistoMotherLikeSignBackInvMassSubNDMPt[fiCut]->Fill(PiMiPiMiNDMBackgroundCandidate.M() - (EventNDMGoodMeson->M() - fPDGMassNDM), PiMiPiMiNDMBackgroundCandidate.Pt(), fWeightJetJetMC);
                fHistoMotherLikeSignBackInvMassFixedPzNDMPt[fiCut]->Fill(PiMiPiMiNDMtmp.M(), PiMiPiMiNDMtmp.Pt(), fWeightJetJetMC);
              }
            }

        } // end of iCurrentPiMi2
      }   // end of iCurrenPiMi
    } else if (mode == 6 || mode == 7){
        // Loops for Pi0Pi+Pi- Sideband mixing and swapping method (Uses same loop with different pi0 candidates)
        for (Int_t iCurrentPiPl = 0; iCurrentPiPl < fPosPionCandidates->GetEntries(); iCurrentPiPl++) {

          AliAODConversionMother EventPiPlGoodMeson = *(AliAODConversionMother *)(fPosPionCandidates->At(iCurrentPiPl));

          for (Int_t iCurrentPiMi = 0; iCurrentPiMi < fNegPionCandidates->GetEntries(); iCurrentPiMi++) {

              AliAODConversionMother EventPiMiGoodMeson = *(AliAODConversionMother *)(fNegPionCandidates->At(iCurrentPiMi));

              // Mass cut on pi+pi-
              if (((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->DoMassCut()) {
                AliAODConversionMother backPiPlPiMiCandidate(&EventPiPlGoodMeson, &EventPiMiGoodMeson);
                if (backPiPlPiMiCandidate.M() >= ((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->GetMassCut()) {
                  continue;
                }
              }

              // Combine Pi+ and Pi0
              AliAODConversionMother PiPlNDMBackgroundCandidate(&EventPiPlGoodMeson, EventNDMGoodMeson);

              // Mass cut (pi0pi+-)
              if (((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->DoMassCut_WithNDM()) {
                AliAODConversionMother PiPlNDMBackgroundCandidate2(&EventPiMiGoodMeson, EventNDMGoodMeson);
                Double_t  Mass_PiPlus_PiZero_Sub           = PiPlNDMBackgroundCandidate.M() - (EventNDMGoodMeson->M() - fPDGMassNDM);
                Double_t  Mass_PiPlus2_PiZero_Sub           = PiPlNDMBackgroundCandidate2.M() - (EventNDMGoodMeson->M() - fPDGMassNDM);
                if ((Mass_PiPlus_PiZero_Sub >= ((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->GetMassCut_WithNDM())||
                     (Mass_PiPlus2_PiZero_Sub >= ((AliPrimaryPionCuts *)fPionCutArray->At(fiCut))->GetMassCut_WithNDM())
                     ) {
                  continue;
                }
              }

              // Create (final) Candidate
              AliAODConversionMother PiPlPiMiNDMBackgroundCandidate(&PiPlNDMBackgroundCandidate, &EventPiMiGoodMeson);

              // Check if candidate survives meson cut
              if (((AliConversionMesonCuts *)fMesonCutArray->At(fiCut))->MesonIsSelected(&PiPlPiMiNDMBackgroundCandidate, kFALSE, ((AliConvEventCuts *)fEventCutArray->At(fiCut))->GetEtaShift())) {

                // Check if candidate survives kinematic cut
                if (KinematicCut(&EventPiPlGoodMeson, &EventPiMiGoodMeson, EventNDMGoodMeson, &PiPlPiMiNDMBackgroundCandidate)) {

                  // Create temporary mesons to be able to fix pz
                  AliAODConversionMother NDMtmp;
                  NDMtmp.SetPxPyPzE(EventNDMGoodMeson->Px(), EventNDMGoodMeson->Py(), EventNDMGoodMeson->Pz(), EventNDMGoodMeson->Energy());
                  FixPzToMatchPDGInvMassNDM(&NDMtmp);
                  AliAODConversionMother PiPlNDMtmp(&EventPiPlGoodMeson, &NDMtmp);
                  AliAODConversionMother PiPlPiMiNDMtmp(&EventPiMiGoodMeson, &PiPlNDMtmp); // Must be two separate lines since second instance depends on first and execution order is not guaranteed

                  // Fill histograms
                  fHistoBackInvMassPt[fiCut]->Fill(PiPlPiMiNDMBackgroundCandidate.M(), PiPlPiMiNDMBackgroundCandidate.Pt(), fWeightJetJetMC);
                  fHistoBackInvMassPtSubNDM[fiCut]->Fill(PiPlPiMiNDMBackgroundCandidate.M() - (EventNDMGoodMeson->M() - fPDGMassNDM), PiPlPiMiNDMBackgroundCandidate.Pt(), fWeightJetJetMC);
                  fHistoBackInvMassPtFixedPzNDM[fiCut]->Fill(PiPlPiMiNDMtmp.M(), PiPlPiMiNDMtmp.Pt(), fWeightJetJetMC);
                }
              }

          } // end of iCurrentPiMi
        }   // end of iCurrenPiPl

    } // end of mode if
  } // end of NDM from current event loop
}

//______________________________________________________________________
Bool_t AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::KinematicCut(AliAODConversionMother *negpion, AliAODConversionMother *pospion, AliAODConversionMother *neutpion, AliAODConversionMother *omega){

  if(fTolerance == -1) return kTRUE;
  if((omega->Pt())<=5.){
    if( (omega->Angle(pospion->Vect()))    < ((2.78715*(TMath::Exp(-0.589934*(omega->Pt()))+0.0519574))*fTolerance) &&
        (omega->Angle(negpion->Vect()))    < ((5.94216*(TMath::Exp(-0.444428*(omega->Pt()))-0.0574076))*fTolerance) &&
        (omega->Angle(neutpion->Vect()))   < ((2.79529*(TMath::Exp(-0.565999*(omega->Pt()))+0.0413576))*fTolerance) &&
        (pospion->Angle(negpion->Vect()))  < ((3.14446*(TMath::Exp(-0.666433*(omega->Pt()))+0.0964309))*fTolerance) &&
        (pospion->Angle(neutpion->Vect())) < ((3.08241*(TMath::Exp(-0.650657*(omega->Pt()))+0.0997539))*fTolerance) &&
        (negpion->Angle(neutpion->Vect())) < ((3.18536*(TMath::Exp(-0.752847*(omega->Pt()))+0.1262780))*fTolerance)
      ){
        return kTRUE;
    }
  }else{
    if( (omega->Angle(pospion->Vect()))    < ((0.459270*(TMath::Exp(-0.126007*(omega->Pt()))+0.100475))*fTolerance) &&
        (omega->Angle(negpion->Vect()))    < ((0.521250*(TMath::Exp(-0.152532*(omega->Pt()))+0.114617))*fTolerance) &&
        (omega->Angle(neutpion->Vect()))   < ((0.409766*(TMath::Exp(-0.108566*(omega->Pt()))+0.103594))*fTolerance) &&
        (pospion->Angle(negpion->Vect()))  < ((0.709206*(TMath::Exp(-0.149072*(omega->Pt()))+0.111345))*fTolerance) &&
        (pospion->Angle(neutpion->Vect())) < ((0.662184*(TMath::Exp(-0.123397*(omega->Pt()))+0.104675))*fTolerance) &&
        (negpion->Angle(neutpion->Vect())) < ((0.730228*(TMath::Exp(-0.120859*(omega->Pt()))+0.105522))*fTolerance)
      ){
        return kTRUE;
    }
  }
  return kFALSE;
}


//______________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::ProcessTrueMesonCandidates(AliAODConversionMother *mesoncand, AliAODConversionMother *TrueNeutralDecayMesonCandidate, AliAODConversionPhoton *TrueVirtualParticleCandidate)
{

  // Process True Mesons

  Bool_t isSameMotherPiPlPiMiNDM   = kFALSE;   // pi+ pi- and pi0 have the same mother
  Bool_t isSameMotherPiPlPiMi      = kFALSE;   // pi+ and pi- have the same mother
  Bool_t isSameMotherPiPlNDM       = kFALSE;   // pi+ and pi0 have the same mother
  Bool_t isSameMotherPiMiNDM       = kFALSE;   // pi- and pi0 have the same mother
  Bool_t isNoSameMother            = kFALSE;   // none of the pions have the same mother

  Bool_t areAllPionsCorrectlyIdentified     = kFALSE;   // All Pion Identifications correct
  Bool_t isPiPlWronglyIdentified            = kFALSE;   // Pi+ Identification not correct
  Bool_t isPiMiWronglyIdentified            = kFALSE;   // Pi- Identification not correct
  Bool_t isPiZeroWronglyIdentified          = kFALSE;   // Pi0 Identification not correct
  Bool_t isMultipleWronglyIdentified        = kFALSE;   // more than one Pion Identification not correct

  Bool_t isTrueMeson                        = kFALSE;   //True analyzed meson
  Bool_t isDifferentMesonContribution       = kFALSE;   //True meson, but NOT analyzed meson
  Bool_t isCombinatoricsMeson               = kFALSE;   //Combinatorics candidate
  Bool_t isContaminationMeson               = kFALSE;   //Contamination candidate

  Bool_t NDMMC_PDGCheck                     = kFALSE;

  Int_t virtualParticleMCLabel = -1;
  virtualParticleMCLabel = TrueVirtualParticleCandidate->GetMCParticleLabel(fMCEvent);
  Int_t virtualParticleMotherLabel = -1;
  //Is set when:
  //if(((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->GetPdgCode() == fPDGCodeNDM) isTrueNDM=kTRUE;
  //if(isTrueNDM){// True Pion
    //Pi0Candidate->SetTrueMesonValue(1);
  Int_t trueMesonFlag  = TrueNeutralDecayMesonCandidate->GetTrueMesonValue();
  Int_t NDMMCLabel     = TrueNeutralDecayMesonCandidate->GetMCLabel();

  Float_t weighted= fWeightJetJetMC;

  if (fDoMesonQA>0){fHistoTrueMesonFlags[fiCut]->Fill(1);} //All candidates

  if ( !(trueMesonFlag == 1 && NDMMCLabel != -1)){ //more understandable: (trueMesonFlag != 1 || NDMMCLabel == -1)
    if((fDoMesonQA>0 ) && (!fDoLightOutput)){
     fHistoTruePiPlPiMiNDMContamination_Crosscheck_InvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
    }
    isPiZeroWronglyIdentified   = kTRUE;
    isContaminationMeson        = kTRUE;
  }
  Int_t NDMMotherLabel =  0;
  TParticle * negativeMC = (TParticle*)TrueVirtualParticleCandidate->GetNegativeMCDaughter(fMCEvent);
  TParticle * positiveMC = (TParticle*)TrueVirtualParticleCandidate->GetPositiveMCDaughter(fMCEvent);

  if (NDMMCLabel == -1){
    NDMMC_PDGCheck = kFALSE;
  } else {
    NDMMotherLabel =  fMCEvent->Particle(NDMMCLabel)->GetMother(0);
    NDMMC_PDGCheck=fMCEvent->Particle(NDMMCLabel)->GetPdgCode()==fPDGCodeNDM;
  }

  Int_t posMotherLabelMC = positiveMC->GetMother(0);
  Int_t negMotherLabelMC = negativeMC->GetMother(0);

  if ( (isPiZeroWronglyIdentified)&&(NDMMC_PDGCheck) ){
    if (fDoMesonQA>0){fHistoTrueMesonFlags[fiCut]->Fill(10);} //Problem with pi0 flag
  }
  // Check case present
  if((TMath::Abs(negativeMC->GetPdgCode())==211) && (TMath::Abs(positiveMC->GetPdgCode())==211) && (NDMMC_PDGCheck)){
    // three pion decay
    areAllPionsCorrectlyIdentified = kTRUE;
    if(virtualParticleMCLabel!=-1){
      // pi+ pi- have same mother
      isSameMotherPiPlPiMi = kTRUE;
      virtualParticleMotherLabel  = virtualParticleMCLabel;
      if(virtualParticleMotherLabel==NDMMotherLabel){
        // all pions from same mother
        isSameMotherPiPlPiMiNDM  = kTRUE;
      } else{
        // only pi+ pi- from same mother
        isCombinatoricsMeson = kTRUE;
      }
    } else{
      //pi+ and pi- do not have same mother -> Combinatorics
      isCombinatoricsMeson = kTRUE;
      if(NDMMotherLabel==negMotherLabelMC && negMotherLabelMC != -1){
        // pi0 and pi- same mother
        isSameMotherPiMiNDM      = kTRUE;
      } else if(NDMMotherLabel==posMotherLabelMC && posMotherLabelMC != -1){
        // pi0 and pi+ same mother
        isSameMotherPiPlNDM      = kTRUE;
      } else{
        // all pions different mother
        isNoSameMother              = kTRUE;
      }
    }
  } else{
    // not a three pion decay, Contamination
    isContaminationMeson        = kTRUE;
    if (!(TMath::Abs(negativeMC->GetPdgCode())==211)){
        isPiMiWronglyIdentified     = kTRUE;
    }
    if (!(TMath::Abs(positiveMC->GetPdgCode())==211)){
        isPiPlWronglyIdentified     = kTRUE;
        if (isPiMiWronglyIdentified){
            isMultipleWronglyIdentified = kTRUE;
        }
    }
    if (!(NDMMC_PDGCheck)){
        isPiZeroWronglyIdentified     = kTRUE;
        if ((isPiMiWronglyIdentified)||(isPiPlWronglyIdentified)){
            isMultipleWronglyIdentified = kTRUE;
        }
    }
  }

  if(areAllPionsCorrectlyIdentified&&isSameMotherPiPlPiMiNDM){
    if(fMCEvent->Particle(NDMMotherLabel)->GetPdgCode()                        == fPDGCodeAnalyzedMeson){
      isTrueMeson                   = kTRUE;
    } else {
      isDifferentMesonContribution  = kTRUE;
    }
  }

  Int_t iNumberOfDeclarationFlags=0;
  if (isTrueMeson){iNumberOfDeclarationFlags++;}
  if (isDifferentMesonContribution){iNumberOfDeclarationFlags++;}
  if (isCombinatoricsMeson){iNumberOfDeclarationFlags++;}
  if (isContaminationMeson){iNumberOfDeclarationFlags++;}
  if (iNumberOfDeclarationFlags!=1){
    if (fDoMesonQA>0){fHistoTrueMesonFlags[fiCut]->Fill(11);} //Problem with meson declaration flag
  }

  if(areAllPionsCorrectlyIdentified&&isSameMotherPiPlPiMiNDM){
    if (fDoMesonQA>0){fHistoTrueMesonFlags[fiCut]->Fill(2);} //Same mother
    if (isTrueMeson){
        if (fDoMesonQA>0){fHistoTrueMesonFlags[fiCut]->Fill(3);} //True
    }
  } else if (areAllPionsCorrectlyIdentified){
    if (fDoMesonQA>0){fHistoTrueMesonFlags[fiCut]->Fill(4);} //Not same mother
  } else if (isContaminationMeson){
    if (fDoMesonQA>0){fHistoTrueMesonFlags[fiCut]->Fill(5);} //Wrongly identified pions
    if (!isMultipleWronglyIdentified){
      if (isPiZeroWronglyIdentified){
        if (fDoMesonQA>0){fHistoTrueMesonFlags[fiCut]->Fill(6);} //Wrongly identified pi0
      } else if (isPiPlWronglyIdentified){
          if (fDoMesonQA>0){fHistoTrueMesonFlags[fiCut]->Fill(7);} //Wrongly identified pi+
      } else if (isPiMiWronglyIdentified){
        if (fDoMesonQA>0){fHistoTrueMesonFlags[fiCut]->Fill(8);} //Wrongly identified pi-
      }
    } else {
      if (fDoMesonQA>0){fHistoTrueMesonFlags[fiCut]->Fill(9);} //Wrongly identified multiple
    }
  }

  // Do things for each case
  if(isTrueMeson){
    // neutral meson was found

    fHistoTrueMotherPiPlPiMiNDMInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);

    // Subtract mass of used NDM candidate and then add PDG mass
    fHistoTrueMotherPiPlPiMiNDMInvMassPtSubNDM[fiCut]->Fill(mesoncand->M()-(TrueNeutralDecayMesonCandidate->M()-fPDGMassNDM),mesoncand->Pt(),weighted);

    // Fix Pz of pi0 candidate to match pi0 PDG mass
    AliAODConversionMother NDMtmp;
    NDMtmp.SetPxPyPzE(TrueNeutralDecayMesonCandidate->Px(), TrueNeutralDecayMesonCandidate->Py(), TrueNeutralDecayMesonCandidate->Pz(), TrueNeutralDecayMesonCandidate->Energy());
    FixPzToMatchPDGInvMassNDM(&NDMtmp);
    AliAODConversionMother mesontmp(&NDMtmp,TrueVirtualParticleCandidate);

    fHistoTrueMotherPiPlPiMiNDMInvMassPtFixedPzNDM[fiCut]->Fill(mesontmp.M(),mesontmp.Pt(),weighted);


    AliAODConversionMother PosPiontmp, NegPiontmp;
    PosPiontmp.SetPxPyPzE(positiveMC->Px(), positiveMC->Py(), positiveMC->Pz(), positiveMC->Energy());
    NegPiontmp.SetPxPyPzE(negativeMC->Px(), negativeMC->Py(), negativeMC->Pz(), negativeMC->Energy());
    if(!fDoLightOutput){
      fHistoTrueAngleSum[fiCut]->Fill(mesoncand->Pt(),((PosPiontmp.Angle(mesoncand->Vect()))+(NegPiontmp.Angle(PosPiontmp.Vect()))+(PosPiontmp.Angle(TrueNeutralDecayMesonCandidate->Vect()))));
      fHistoTrueHNMesonPtvsNDMPt[fiCut]->Fill(mesoncand->Pt(),TrueNeutralDecayMesonCandidate->Pt(),weighted);
    }


    // Fill tree to get info about event that the eta was found in
    if(fDoMesonQA>1 && (!fDoLightOutput)){
      fV0MultiplicityHNMEvent = fMCEvent->GetNumberOfV0s();
      fTrackMultiplicityHNMEvent = fMCEvent->GetNumberOfTracks();
      fZVertexHNMEvent = fMCEvent->GetPrimaryVertex()->GetZ();
      fPtHNM = mesoncand->Pt();

      fTreeEventInfoHNM[fiCut]->Fill();
    }
    if (CheckVectorForDoubleCount(fVectorDoubleCountTrueHNMs,NDMMotherLabel) && (!fDoLightOutput)) fHistoDoubleCountTrueHNMInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt());
  } else if (isDifferentMesonContribution) {
    //True Meson, but the analyzed one
    fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromDifferent[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
    if((fDoMesonQA>0 ) && (!fDoLightOutput)){
      if((fMCEvent->Particle(posMotherLabelMC)->GetPdgCode()                     == 223)||
         (fMCEvent->Particle(posMotherLabelMC)->GetPdgCode()                     == 221)){
        // pi+pi- come from eta
        fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromEtaOmega[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      } else if(fMCEvent->Particle(posMotherLabelMC)->GetPdgCode()              == 113){
        // pi+pi- come from rho0
        fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromRho[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      } else if(fMCEvent->Particle(posMotherLabelMC)->GetPdgCode()              == 331){
        // pi+pi- come from eta prime
        fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromEtaPrime[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      } else if(fMCEvent->Particle(posMotherLabelMC)->GetPdgCode()              == 310){
        // pi+pi- come from K0 short
        fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromK0s[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      } else if(fMCEvent->Particle(posMotherLabelMC)->GetPdgCode()              == 130){
        // pi+pi- come from K0 long
        fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromK0l[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      } else{
        // pi+pi- come from something else
        fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromOther[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      }
    }
  } else if (isCombinatoricsMeson) {
    fHistoTruePiPlPiMiNDMCombinatoricalInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
    if(isSameMotherPiPlPiMi &&  (fDoMesonQA>0 ) && (!fDoLightOutput)){
      if(fMCEvent->Particle(posMotherLabelMC)->GetPdgCode()                     == 221){
        // pi+pi- come from eta
        fHistoTruePiPlPiMiSameMotherFromEtaInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      } else if(fMCEvent->Particle(posMotherLabelMC)->GetPdgCode()              == 223){
        // pi+pi- come from omega
        fHistoTruePiPlPiMiSameMotherFromOmegaInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      } else if(fMCEvent->Particle(posMotherLabelMC)->GetPdgCode()              == 113){
       // pi+pi- come from rho0
        fHistoTruePiPlPiMiSameMotherFromRhoInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      } else if(fMCEvent->Particle(posMotherLabelMC)->GetPdgCode()              == 331){
        // pi+pi- come from eta prime
        fHistoTruePiPlPiMiSameMotherFromEtaPrimeInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      } else if(fMCEvent->Particle(posMotherLabelMC)->GetPdgCode()              == 310){
        // pi+pi- come from K0 short
        fHistoTruePiPlPiMiSameMotherFromK0sInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      } else if(fMCEvent->Particle(posMotherLabelMC)->GetPdgCode()              == 130){
        // pi+pi- come from K0 short
        fHistoTruePiPlPiMiSameMotherFromK0lInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      } else {
        // pi+pi- come from something else
        fHistoTruePiPlPiMiSameMotherFromOtherlInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      }
    } else if(isSameMotherPiMiNDM  &&  (fDoMesonQA>0 ) && (!fDoLightOutput)){
      if(fMCEvent->Particle(NDMMotherLabel)->GetPdgCode()                       == 221){
        // pi0pi- come from eta
        fHistoTruePiMiPiZeroSameMotherFromEtaInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      } else if(fMCEvent->Particle(NDMMotherLabel)->GetPdgCode()                == 223){
        // pi0pi- come from omega
        fHistoTruePiMiPiZeroSameMotherFromOmegaInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      } else if(fMCEvent->Particle(NDMMotherLabel)->GetPdgCode()                ==-213){
        // pi0pi- come from rho-
        fHistoTruePiMiPiZeroSameMotherFromRhoInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      } else if(fMCEvent->Particle(NDMMotherLabel)->GetPdgCode()                == 130){
        // pi0pi- come from rho-
        fHistoTruePiMiPiZeroSameMotherFromK0lInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      } else {
        // pi0pi- come from something else
        fHistoTruePiMiPiZeroSameMotherFromOtherlInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      }
    } else if(isSameMotherPiPlNDM  &&  (fDoMesonQA>0 ) && (!fDoLightOutput)){
      if(fMCEvent->Particle(posMotherLabelMC)->GetPdgCode()                     == 221){
        // pi+pi0 come from eta
        fHistoTruePiPlPiZeroSameMotherFromEtaInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      } else if(fMCEvent->Particle(posMotherLabelMC)->GetPdgCode()              == 223){
        // pi+pi0 come from omega
        fHistoTruePiPlPiZeroSameMotherFromOmegaInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      } else if(fMCEvent->Particle(posMotherLabelMC)->GetPdgCode()              == 213) {
        // pi+pi0 come from rho+
        fHistoTruePiPlPiZeroSameMotherFromRhoInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      } else if(fMCEvent->Particle(posMotherLabelMC)->GetPdgCode()              == 130) {
        // pi+pi0 come from rho+
        fHistoTruePiPlPiZeroSameMotherFromK0lInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      } else {
        // pi+pi0 come from something else
        fHistoTruePiPlPiZeroSameMotherFromOtherInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      }
    } else if(isNoSameMother  &&  (fDoMesonQA>0 ) && (!fDoLightOutput)){
      // no same mother purecombinatorical
      fHistoTruePiPlPiMiNDMPureCombinatoricalInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
    }
  } else if (isContaminationMeson) {
    // no pi pi pi decay contamination
    fHistoTruePiPlPiMiNDMContaminationInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
    if((fDoMesonQA>0 ) && (!fDoLightOutput)){
      if (!isMultipleWronglyIdentified){
        if (isPiPlWronglyIdentified){
          fHistoTruePiPlPiMiNDMContamination_PiPl_InvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
        }
        if (isPiMiWronglyIdentified){
          fHistoTruePiPlPiMiNDMContamination_PiMi_InvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
        }
        if (isPiZeroWronglyIdentified){
          fHistoTruePiPlPiMiNDMContamination_Pi0_InvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
        }
      } else {
        fHistoTruePiPlPiMiNDMContamination_multipel_InvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      }
    }
  }
}

//______________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::ProcessTrueMesonCandidatesAOD(AliAODConversionMother *mesoncand, AliAODConversionMother *TrueNeutralDecayMesonCandidate, AliAODConversionPhoton *TrueVirtualParticleCandidate)
{

  // Process True Mesons

  Bool_t isSameMotherPiPlPiMiNDM   = kFALSE;   // pi+ pi- and pi0 have the same mother
  Bool_t isSameMotherPiPlPiMi      = kFALSE;   // pi+ and pi- have the same mother
  Bool_t isSameMotherPiPlNDM       = kFALSE;   // pi+ and pi0 have the same mother
  Bool_t isSameMotherPiMiNDM       = kFALSE;   // pi- and pi0 have the same mother
  Bool_t isNoSameMother            = kFALSE;   // none of the pions have the same mother

  Bool_t areAllPionsCorrectlyIdentified     = kFALSE;   // All Pion Identifications correct
  Bool_t isPiPlWronglyIdentified            = kFALSE;   // Pi+ Identification not correct
  Bool_t isPiMiWronglyIdentified            = kFALSE;   // Pi- Identification not correct
  Bool_t isPiZeroWronglyIdentified          = kFALSE;   // Pi0 Identification not correct
  Bool_t isMultipleWronglyIdentified        = kFALSE;   // more than one Pion Identification not correct

  Bool_t isTrueMeson                        = kFALSE;   //True analyzed meson
  Bool_t isDifferentMesonContribution       = kFALSE;   //True meson, but NOT analyzed meson
  Bool_t isCombinatoricsMeson               = kFALSE;   //Combinatorics candidate
  Bool_t isContaminationMeson               = kFALSE;   //Contamination candidate

  Bool_t NDMMC_PDGCheck                     = kFALSE;

  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));

  Int_t virtualParticleMCLabel = -1;
  Int_t virtualParticleMotherLabel = -1;
  //Is set when:
  //if(((AliAODMCParticle*)AODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == fPDGCodeNDM) isTrueNDM=kTRUE;
  //if(isTrueNDM){// True Pion
    //Pi0Candidate->SetTrueMesonValue(1);
  Int_t trueMesonFlag  = TrueNeutralDecayMesonCandidate->GetTrueMesonValue();
  Int_t NDMMCLabel     = TrueNeutralDecayMesonCandidate->GetMCLabel();

  Float_t weighted= fWeightJetJetMC;

  if (fDoMesonQA>0){fHistoTrueMesonFlags[fiCut]->Fill(1);} //All candidates

  if ( !(trueMesonFlag == 1 && NDMMCLabel != -1)){ //more understandable: (trueMesonFlag != 1 || NDMMCLabel == -1)
    if((fDoMesonQA>0 ) && (!fDoLightOutput)){
      fHistoTruePiPlPiMiNDMContamination_Crosscheck_InvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
    }
    isPiZeroWronglyIdentified   = kTRUE;
    isContaminationMeson        = kTRUE;
    //return;
  }
  Int_t NDMMotherLabel = 0;
  AliAODMCParticle *negativeMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueVirtualParticleCandidate->GetMCLabelNegative())); // pi-
  AliAODMCParticle *positiveMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueVirtualParticleCandidate->GetMCLabelPositive())); // pi+
  AliAODMCParticle *NDMMC      = NULL;
  if (NDMMCLabel == -1){
    NDMMC_PDGCheck = kFALSE;
  } else {
    NDMMotherLabel = (static_cast<AliAODMCParticle*>(AODMCTrackArray->At(NDMMCLabel)))->GetMother();
    NDMMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(NDMMCLabel)); // pi0
    NDMMC_PDGCheck=NDMMC->GetPdgCode()==fPDGCodeNDM;
  }


  if(positiveMC->GetMother()>-1&&(negativeMC->GetMother() == positiveMC->GetMother())){
    virtualParticleMCLabel = positiveMC->GetMother();
  }


  Int_t posMotherLabelMC = positiveMC->GetMother();
  Int_t negMotherLabelMC = negativeMC->GetMother();

  if ( (isPiZeroWronglyIdentified)&&((NDMMC_PDGCheck)) ){
    if (fDoMesonQA>0){fHistoTrueMesonFlags[fiCut]->Fill(10);} //Problem with pi0 flag
  }
  // Check case present
  if((TMath::Abs(negativeMC->GetPdgCode())==211) && (TMath::Abs(positiveMC->GetPdgCode())==211) && (NDMMC_PDGCheck)){
    // three pion decay, Combinatorics and trues
    areAllPionsCorrectlyIdentified = kTRUE;
    if(virtualParticleMCLabel!=-1){
      // pi+ pi- have same mother
      isSameMotherPiPlPiMi = kTRUE;
      virtualParticleMotherLabel  = virtualParticleMCLabel;
      if(virtualParticleMotherLabel==NDMMotherLabel){
        // all pions from same mother
        isSameMotherPiPlPiMiNDM  = kTRUE;
      } else{
        // only pi+ pi- from same mother -> Combinatorics
        isCombinatoricsMeson = kTRUE;
      }
    } else{
      //pi+ and pi- do not have same mother -> Combinatorics
      isCombinatoricsMeson = kTRUE;
      if(NDMMotherLabel==negMotherLabelMC && negMotherLabelMC != -1){
        // pi0 and pi- same mother
        isSameMotherPiMiNDM      = kTRUE;
      } else if(NDMMotherLabel==posMotherLabelMC && posMotherLabelMC != -1){
        // pi0 and pi+ same mother
        isSameMotherPiPlNDM      = kTRUE;
      } else{
        // all pions different mother
        isNoSameMother              = kTRUE;
      }
    }
  } else{
    // not a three pion decay, Contamination
    isContaminationMeson = kTRUE;
    if (!(TMath::Abs(negativeMC->GetPdgCode())==211)){
        isPiMiWronglyIdentified     = kTRUE;
    }
    if (!(TMath::Abs(positiveMC->GetPdgCode())==211)){
        isPiPlWronglyIdentified     = kTRUE;
        if (isPiMiWronglyIdentified){
            isMultipleWronglyIdentified = kTRUE;
        }
    }
    if (!(NDMMC_PDGCheck)){
        isPiZeroWronglyIdentified     = kTRUE;
        if ((isPiMiWronglyIdentified)||(isPiPlWronglyIdentified)){
            isMultipleWronglyIdentified = kTRUE;
        }
    }
  }

  if(areAllPionsCorrectlyIdentified&&isSameMotherPiPlPiMiNDM){
    if((static_cast<AliAODMCParticle*>(AODMCTrackArray->At(NDMMotherLabel)))->GetPdgCode()                        == fPDGCodeAnalyzedMeson){
      isTrueMeson                   = kTRUE;
    } else {
      isDifferentMesonContribution  = kTRUE;
    }
  }

  Int_t iNumberOfDeclarationFlags=0;
  if (isTrueMeson){iNumberOfDeclarationFlags++;}
  if (isDifferentMesonContribution){iNumberOfDeclarationFlags++;}
  if (isCombinatoricsMeson){iNumberOfDeclarationFlags++;}
  if (isContaminationMeson){iNumberOfDeclarationFlags++;}
  if (iNumberOfDeclarationFlags!=1){
      if (fDoMesonQA>0){fHistoTrueMesonFlags[fiCut]->Fill(11);} //Problem with meson declaration flag
  }

  if(areAllPionsCorrectlyIdentified&&isSameMotherPiPlPiMiNDM){
    if (fDoMesonQA>0){fHistoTrueMesonFlags[fiCut]->Fill(2);} //Same mother
    if (isTrueMeson){
        if (fDoMesonQA>0){fHistoTrueMesonFlags[fiCut]->Fill(3);} //True
    }
  } else if (areAllPionsCorrectlyIdentified){
    if (fDoMesonQA>0){fHistoTrueMesonFlags[fiCut]->Fill(4);} //Not same mother
  } else if (isContaminationMeson){
    if (fDoMesonQA>0){fHistoTrueMesonFlags[fiCut]->Fill(5);} //Wrongly identified pions
    if (!isMultipleWronglyIdentified){
      if (isPiZeroWronglyIdentified){
        if (fDoMesonQA>0){fHistoTrueMesonFlags[fiCut]->Fill(6);} //Wrongly identified pi0
      } else if (isPiPlWronglyIdentified){
          if (fDoMesonQA>0){fHistoTrueMesonFlags[fiCut]->Fill(7);} //Wrongly identified pi+
      } else if (isPiMiWronglyIdentified){
        if (fDoMesonQA>0){fHistoTrueMesonFlags[fiCut]->Fill(8);} //Wrongly identified pi-
      }
    } else {
      if (fDoMesonQA>0){fHistoTrueMesonFlags[fiCut]->Fill(9);} //Wrongly identified multiple
    }
  }

  // Do things for each case
  if(isTrueMeson){
    // neutral meson was found
    fHistoTrueMotherPiPlPiMiNDMInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);

    // Subtract mass of used NDM candidate and then add PDG mass
    fHistoTrueMotherPiPlPiMiNDMInvMassPtSubNDM[fiCut]->Fill(mesoncand->M()-(TrueNeutralDecayMesonCandidate->M()-fPDGMassNDM),mesoncand->Pt(),weighted);

    // Fix Pz of pi0 candidate to match pi0 PDG mass
    AliAODConversionMother NDMtmp;
    NDMtmp.SetPxPyPzE(TrueNeutralDecayMesonCandidate->Px(), TrueNeutralDecayMesonCandidate->Py(), TrueNeutralDecayMesonCandidate->Pz(), TrueNeutralDecayMesonCandidate->Energy());
    FixPzToMatchPDGInvMassNDM(&NDMtmp);
    AliAODConversionMother mesontmp(&NDMtmp,TrueVirtualParticleCandidate);

    fHistoTrueMotherPiPlPiMiNDMInvMassPtFixedPzNDM[fiCut]->Fill(mesontmp.M(),mesontmp.Pt(),weighted);

    if(fDoMesonQA>0){
      //Dalitz plot
      TLorentzVector PosPionTLVtmp;
      TLorentzVector NegPionTLVtmp;
      TLorentzVector PosNegPionTLVtmp;

      PosPionTLVtmp.SetPxPyPzE (positiveMC->Px(), positiveMC->Py(), positiveMC->Pz(), positiveMC->E() );
      NegPionTLVtmp.SetPxPyPzE (negativeMC->Px(), negativeMC->Py(), negativeMC->Pz(), negativeMC->E() );
      PosNegPionTLVtmp = PosPionTLVtmp + NegPionTLVtmp;

      TLorentzVector NDMTLVtmp;
      TLorentzVector NDMSubTLVtmp;
      TLorentzVector PosPionNDMTLVtmp;
      TLorentzVector NegPionNDMTLVtmp;
      TLorentzVector PosPionNDMSubTLVtmp;
      TLorentzVector NegPionNDMSubTLVtmp;

      NDMTLVtmp.SetPxPyPzE( NDMtmp.Px(), NDMtmp.Py(), NDMtmp.Pz(), NDMtmp.E() );
      NDMSubTLVtmp.SetPxPyPzE (TrueNeutralDecayMesonCandidate->Px(), TrueNeutralDecayMesonCandidate->Py(), TrueNeutralDecayMesonCandidate->Pz(), TrueNeutralDecayMesonCandidate->Energy());
      PosPionNDMTLVtmp = PosPionTLVtmp + NDMTLVtmp;
      NegPionNDMTLVtmp = NegPionTLVtmp + NDMTLVtmp;
      PosPionNDMSubTLVtmp = PosPionTLVtmp + NDMSubTLVtmp;
      NegPionNDMSubTLVtmp = NegPionTLVtmp + NDMSubTLVtmp;

      //Dalitz All Pt
      if (enableDalitzAllPt){
        fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM[fiCut]->Fill(PosNegPionTLVtmp.M(), PosPionNDMTLVtmp.M() );
        fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM[fiCut]->Fill(PosNegPionTLVtmp.M(), PosPionNDMSubTLVtmp.M() - (NDMSubTLVtmp.M() - fPDGMassNDM) );
        fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM[fiCut]->Fill(PosNegPionTLVtmp.M(), NegPionNDMTLVtmp.M() );
        fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM[fiCut]->Fill(PosNegPionTLVtmp.M(), NegPionNDMSubTLVtmp.M() - (NDMSubTLVtmp.M() - fPDGMassNDM) );
      }
      //Dalitz Low Pt
      if (enableDalitzLowPt){
        if ((mesoncand->Pt()>HistoDalitzPtRangeMin_LowPt)&&(mesoncand->Pt()<HistoDalitzPtRangeMax_LowPt)){
          fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM_LowPt[fiCut]->Fill(PosNegPionTLVtmp.M(), PosPionNDMTLVtmp.M() );
          fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM_LowPt[fiCut]->Fill(PosNegPionTLVtmp.M(), PosPionNDMSubTLVtmp.M() - (NDMSubTLVtmp.M() - fPDGMassNDM) );
          fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM_LowPt[fiCut]->Fill(PosNegPionTLVtmp.M(), NegPionNDMTLVtmp.M() );
          fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM_LowPt[fiCut]->Fill(PosNegPionTLVtmp.M(), NegPionNDMSubTLVtmp.M() - (NDMSubTLVtmp.M() - fPDGMassNDM) );
        }
      }
      //Dalitz Mid Pt
      if (enableDalitzMidPt){
        if ((mesoncand->Pt()>HistoDalitzPtRangeMin_MidPt)&&(mesoncand->Pt()<HistoDalitzPtRangeMax_MidPt)){
          fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM_MidPt[fiCut]->Fill(PosNegPionTLVtmp.M(), PosPionNDMTLVtmp.M() );
          fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM_MidPt[fiCut]->Fill(PosNegPionTLVtmp.M(), PosPionNDMSubTLVtmp.M() - (NDMSubTLVtmp.M() - fPDGMassNDM) );
          fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM_MidPt[fiCut]->Fill(PosNegPionTLVtmp.M(), NegPionNDMTLVtmp.M() );
          fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM_MidPt[fiCut]->Fill(PosNegPionTLVtmp.M(), NegPionNDMSubTLVtmp.M() - (NDMSubTLVtmp.M() - fPDGMassNDM) );
        }
      }
      //Dalitz High Pt
      if (enableDalitzHighPt){
        if ((mesoncand->Pt()>HistoDalitzPtRangeMin_HighPt)&&(mesoncand->Pt()<HistoDalitzPtRangeMax_HighPt)){
          fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosFixedPzNDM_HighPt[fiCut]->Fill(PosNegPionTLVtmp.M(), PosPionNDMTLVtmp.M() );
          fHistoTrueMotherPiPlPiMiNDMDalitzPlotPosSubNDM_HighPt[fiCut]->Fill(PosNegPionTLVtmp.M(), PosPionNDMSubTLVtmp.M() - (NDMSubTLVtmp.M() - fPDGMassNDM) );
          fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegFixedPzNDM_HighPt[fiCut]->Fill(PosNegPionTLVtmp.M(), NegPionNDMTLVtmp.M() );
          fHistoTrueMotherPiPlPiMiNDMDalitzPlotNegSubNDM_HighPt[fiCut]->Fill(PosNegPionTLVtmp.M(), NegPionNDMSubTLVtmp.M() - (NDMSubTLVtmp.M() - fPDGMassNDM) );
        }
      }
    }

    AliAODConversionMother PosPiontmp, NegPiontmp;
    PosPiontmp.SetPxPyPzE(positiveMC->Px(), positiveMC->Py(), positiveMC->Pz(), positiveMC->E());
    NegPiontmp.SetPxPyPzE(negativeMC->Px(), negativeMC->Py(), negativeMC->Pz(), negativeMC->E());
    if(!fDoLightOutput){
      fHistoTrueAngleSum[fiCut]->Fill(mesoncand->Pt(),((PosPiontmp.Angle(mesoncand->Vect()))+(NegPiontmp.Angle(PosPiontmp.Vect()))+(PosPiontmp.Angle(TrueNeutralDecayMesonCandidate->Vect()))));
      fHistoTrueHNMesonPtvsNDMPt[fiCut]->Fill(mesoncand->Pt(),TrueNeutralDecayMesonCandidate->Pt(),weighted);
    }
    // Fill tree to get info about event that the eta was found in
    if(fDoMesonQA>1 && (!fDoLightOutput)){
      fV0MultiplicityHNMEvent = fMCEvent->GetNumberOfV0s();
      fTrackMultiplicityHNMEvent = fMCEvent->GetNumberOfTracks();
      fZVertexHNMEvent = fMCEvent->GetPrimaryVertex()->GetZ();
      fPtHNM = mesoncand->Pt();

      fTreeEventInfoHNM[fiCut]->Fill();
    }
    if (CheckVectorForDoubleCount(fVectorDoubleCountTrueHNMs,NDMMotherLabel) && (!fDoLightOutput)) fHistoDoubleCountTrueHNMInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt());
  } else if (isDifferentMesonContribution){
    //True Meson, but the analyzed one
    fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromDifferent[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
    if((fDoMesonQA>0 ) && (!fDoLightOutput)){
      if(((static_cast<AliAODMCParticle*>(AODMCTrackArray->At(posMotherLabelMC)))->GetPdgCode()                     == 223)||
              ((static_cast<AliAODMCParticle*>(AODMCTrackArray->At(posMotherLabelMC)))->GetPdgCode()                     == 221)){
        // pi+pi- come from eta
        fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromEtaOmega[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      } else if((static_cast<AliAODMCParticle*>(AODMCTrackArray->At(posMotherLabelMC)))->GetPdgCode()              == 113){
        // pi+pi- come from rho0
        fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromRho[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      } else if((static_cast<AliAODMCParticle*>(AODMCTrackArray->At(posMotherLabelMC)))->GetPdgCode()              == 331){
        // pi+pi- come from eta prime
        fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromEtaPrime[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      } else if((static_cast<AliAODMCParticle*>(AODMCTrackArray->At(posMotherLabelMC)))->GetPdgCode()              == 310){
        // pi+pi- come from K0 short
        fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromK0s[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      } else if((static_cast<AliAODMCParticle*>(AODMCTrackArray->At(posMotherLabelMC)))->GetPdgCode()              == 130){
        // pi+pi- come from K0 long
        fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromK0l[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      } else{
        // pi+pi- come from something else
        fHistoTrueMotherPiPlPiMiNDMInvMassPt_FromOther[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      }
    }
  } else if (isCombinatoricsMeson) {
      fHistoTruePiPlPiMiNDMCombinatoricalInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      if(isSameMotherPiPlPiMi &&  (fDoMesonQA>0 ) && (!fDoLightOutput)){
        if((static_cast<AliAODMCParticle*>(AODMCTrackArray->At(posMotherLabelMC)))->GetPdgCode()                     == 221){
          // pi+pi- come from eta
          fHistoTruePiPlPiMiSameMotherFromEtaInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
        } else if((static_cast<AliAODMCParticle*>(AODMCTrackArray->At(posMotherLabelMC)))->GetPdgCode()              == 223){
          // pi+pi- come from omega
          fHistoTruePiPlPiMiSameMotherFromOmegaInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
        } else if((static_cast<AliAODMCParticle*>(AODMCTrackArray->At(posMotherLabelMC)))->GetPdgCode()              == 113){
          // pi+pi- come from rho0
          fHistoTruePiPlPiMiSameMotherFromRhoInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
        } else if((static_cast<AliAODMCParticle*>(AODMCTrackArray->At(posMotherLabelMC)))->GetPdgCode()              == 331){
          // pi+pi- come from eta prime
          fHistoTruePiPlPiMiSameMotherFromEtaPrimeInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
        } else if((static_cast<AliAODMCParticle*>(AODMCTrackArray->At(posMotherLabelMC)))->GetPdgCode()              == 310){
          // pi+pi- come from K0 short
          fHistoTruePiPlPiMiSameMotherFromK0sInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
        } else if((static_cast<AliAODMCParticle*>(AODMCTrackArray->At(posMotherLabelMC)))->GetPdgCode()              == 130){
          // pi+pi- come from K0 long
          fHistoTruePiPlPiMiSameMotherFromK0lInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
        } else{
          // pi+pi- come from something else
          fHistoTruePiPlPiMiSameMotherFromOtherlInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
        }
      } else if(isSameMotherPiMiNDM  &&  (fDoMesonQA>0 ) && (!fDoLightOutput)){
        if((static_cast<AliAODMCParticle*>(AODMCTrackArray->At(NDMMotherLabel)))->GetPdgCode()                       == 221){
          // pi0pi- come from eta
          fHistoTruePiMiPiZeroSameMotherFromEtaInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
        } else if((static_cast<AliAODMCParticle*>(AODMCTrackArray->At(NDMMotherLabel)))->GetPdgCode()                == 223){
          // pi0pi- come from omega
          fHistoTruePiMiPiZeroSameMotherFromOmegaInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
        } else if((static_cast<AliAODMCParticle*>(AODMCTrackArray->At(NDMMotherLabel)))->GetPdgCode()                ==-213){
          // pi0pi- come from rho-
          fHistoTruePiMiPiZeroSameMotherFromRhoInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
        } else if((static_cast<AliAODMCParticle*>(AODMCTrackArray->At(NDMMotherLabel)))->GetPdgCode()                == 130){
          // pi0pi- come from K0l
          fHistoTruePiMiPiZeroSameMotherFromK0lInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
        } else{
          // pi0pi- come from something else
          fHistoTruePiMiPiZeroSameMotherFromOtherlInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
        }
      } else if(isSameMotherPiPlNDM  &&  (fDoMesonQA>0 ) && (!fDoLightOutput)){
        if((static_cast<AliAODMCParticle*>(AODMCTrackArray->At(posMotherLabelMC)))->GetPdgCode()                     == 221){
          // pi+pi0 come from eta
          fHistoTruePiPlPiZeroSameMotherFromEtaInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
        } else if((static_cast<AliAODMCParticle*>(AODMCTrackArray->At(posMotherLabelMC)))->GetPdgCode()              == 223){
          // pi+pi0 come from omega
          fHistoTruePiPlPiZeroSameMotherFromOmegaInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
        } else if((static_cast<AliAODMCParticle*>(AODMCTrackArray->At(posMotherLabelMC)))->GetPdgCode()              == 213) {
          // pi+pi0 come from rho+
          fHistoTruePiPlPiZeroSameMotherFromRhoInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
        } else if((static_cast<AliAODMCParticle*>(AODMCTrackArray->At(posMotherLabelMC)))->GetPdgCode()              == 130) {
          // pi+pi0 come from K0l
          fHistoTruePiPlPiZeroSameMotherFromK0lInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
        } else{
          // pi+pi0 come from something else
          fHistoTruePiPlPiZeroSameMotherFromOtherInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
        }
      } else if(isNoSameMother  &&  (fDoMesonQA>0 ) && (!fDoLightOutput)){
        // no same mother purecombinatorical
        fHistoTruePiPlPiMiNDMPureCombinatoricalInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      }
  } else if (isContaminationMeson) {
      // no pi pi pi decay contamination
      fHistoTruePiPlPiMiNDMContaminationInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
    if((fDoMesonQA>0 ) && (!fDoLightOutput)){
      if (!isMultipleWronglyIdentified){
        if (isPiPlWronglyIdentified){
          fHistoTruePiPlPiMiNDMContamination_PiPl_InvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
        }
        if (isPiMiWronglyIdentified){
          fHistoTruePiPlPiMiNDMContamination_PiMi_InvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
        }
        if (isPiZeroWronglyIdentified){
          fHistoTruePiPlPiMiNDMContamination_Pi0_InvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
        }
      } else {
        fHistoTruePiPlPiMiNDMContamination_multipel_InvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      }
    }
  }
}


//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::UpdateEventByEventData(){
  //see header file for documentation

  Int_t method = 1;
  if( method == 1 ) {
    if(fPosPionCandidates->GetEntries() >0 && fNegPionCandidates->GetEntries() >0){
      if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
        fBGHandlerPiPl[fiCut]->AddMesonEvent(fPosPionCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fV0Reader->GetNumberOfPrimaryTracks(),0);
        fBGHandlerPiMi[fiCut]->AddMesonEvent(fNegPionCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fV0Reader->GetNumberOfPrimaryTracks(),0);
      } else { // means we use #V0s for multiplicity
        if (fNDMRecoMode < 2){
          fBGHandlerPiPl[fiCut]->AddMesonEvent(fPosPionCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fGoodConvGammas->GetEntries(),0);
          fBGHandlerPiMi[fiCut]->AddMesonEvent(fNegPionCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fGoodConvGammas->GetEntries(),0);
        }else {
          fBGHandlerPiPl[fiCut]->AddMesonEvent(fPosPionCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fClusterCandidates->GetEntries(),0);
          fBGHandlerPiMi[fiCut]->AddMesonEvent(fNegPionCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fClusterCandidates->GetEntries(),0);
        }
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::MoveParticleAccordingToVertex(AliAODConversionMother* particle,const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex){
  //see header file for documentation
  Double_t dx = vertex->fX - fInputEvent->GetPrimaryVertex()->GetX();
  Double_t dy = vertex->fY - fInputEvent->GetPrimaryVertex()->GetY();
  Double_t dz = vertex->fZ - fInputEvent->GetPrimaryVertex()->GetZ();

  Double_t movedPlace[3] = {particle->GetProductionX() - dx,particle->GetProductionY() - dy,particle->GetProductionZ() - dz};
  particle->SetProductionPoint(movedPlace);
}

//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::FixPzToMatchPDGInvMassNDM(AliAODConversionMother* particle) {

  Double_t px = particle->Px();
  Double_t py = particle->Py();
  Int_t signPz = particle->Pz()<0?-1:1;
  Double_t energy = particle->Energy();
  Double_t pz = signPz*TMath::Sqrt(TMath::Abs(pow(fPDGMassNDM,2)-pow(energy,2)+pow(px,2)+pow(py,2)));
  particle->SetPxPyPzE(px,py,pz,energy);

  return;
}
//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::FixPzVecToMatchPDGInvMass(TLorentzVector* track) {

  Double_t px = track->Px();
  Double_t py = track->Py();
  Int_t signPz = track->Pz()<0?-1:1;
  Double_t energy = track->E();
  Double_t pz = signPz*TMath::Sqrt(TMath::Abs(pow(fPDGMassChargedPion,2)-pow(energy,2)+pow(px,2)+pow(py,2)));
  track->SetPxPyPzE(px,py,pz,energy);
  return;
}

//_____________________________________________________________________________________
Bool_t AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::IsEtaPrimePiPlPiMiEtaDaughter( Int_t label ) const {
  //
  // Returns true if the particle comes from eta -> pi+ pi- gamma
  //
  if(label<0) return kFALSE;
  Int_t motherLabel = fMCEvent->Particle( label )->GetMother(0);
  if( motherLabel < 0 || motherLabel >= fMCEvent->GetNumberOfTracks() ) return kFALSE;

  TParticle* mother = fMCEvent->Particle( motherLabel );
  if( mother->GetPdgCode() != 331 ) return kFALSE;
  if( IsPiPlPiMiEtaDecay( mother ) ) return kTRUE;
  return kFALSE;
}

//_____________________________________________________________________________________
Bool_t AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::IsEtaPrimePiPlPiMiEtaDaughterAOD(TClonesArray* trackArray, Int_t label ) const {
  //
  // Returns true if the particle comes from eta -> pi+ pi- gamma
  //
  if(label<0) return kFALSE;
  Int_t motherLabel = (static_cast<AliAODMCParticle*>(trackArray->At(label)))->GetMother();
  if( motherLabel < 0 || motherLabel >= trackArray->GetEntriesFast() ) return kFALSE;

  AliAODMCParticle* mother = static_cast<AliAODMCParticle*>(trackArray->At(motherLabel));
  if( mother->GetPdgCode() != 331 ) return kFALSE;
  if( IsPiPlPiMiEtaDecayAOD(trackArray,mother) ) return kTRUE;
  return kFALSE;
}

//_____________________________________________________________________________________
Bool_t AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::IsEtaPiPlPiMiPiZeroDaughter( Int_t label ) const {
  //
  // Returns true if the particle comes from eta -> pi+ pi- gamma
  //
    if(label<0) return kFALSE;
  Int_t motherLabel = fMCEvent->Particle( label )->GetMother(0);
  if( motherLabel < 0 || motherLabel >= fMCEvent->GetNumberOfTracks() ) return kFALSE;

  TParticle* mother = fMCEvent->Particle( motherLabel );
  if( mother->GetPdgCode() != 221 ) return kFALSE;
  if( IsPiPlPiMiPiZeroDecay( mother ) ) return kTRUE;
  return kFALSE;
}
//_____________________________________________________________________________________
Bool_t AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::IsEtaPiPlPiMiPiZeroDaughterAOD(TClonesArray* trackArray, Int_t label ) const {
  //
  // Returns true if the particle comes from eta -> pi+ pi- gamma
  //
    if(label<0) return kFALSE;
    Int_t motherLabel = (static_cast<AliAODMCParticle*>(trackArray->At(label)))->GetMother();
  if( motherLabel < 0 || motherLabel >= trackArray->GetEntriesFast()  ) return kFALSE;

  AliAODMCParticle* mother = static_cast<AliAODMCParticle*>(trackArray->At(motherLabel));
  if( mother->GetPdgCode() != 221 ) return kFALSE;
  if( IsPiPlPiMiPiZeroDecayAOD(trackArray,mother)) return kTRUE;
  return kFALSE;
}

//_____________________________________________________________________________________
Bool_t AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::IsOmegaPiPlPiMiPiZeroDaughter( Int_t label ) const {
  //
  // Returns true if the particle comes from eta -> pi+ pi- gamma
  //
  if(label<0) return kFALSE;
  Int_t motherLabel = fMCEvent->Particle( label )->GetMother(0);
  if( motherLabel < 0 || motherLabel >= fMCEvent->GetNumberOfTracks() ) return kFALSE;

  TParticle* mother = fMCEvent->Particle( motherLabel );
  if( mother->GetPdgCode() != 223 ) return kFALSE;
  if( IsPiPlPiMiPiZeroDecay( mother ) ) return kTRUE;
  return kFALSE;
}

//_____________________________________________________________________________________
Bool_t AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::IsOmegaPiPlPiMiPiZeroDaughterAOD(TClonesArray* trackArray, Int_t label ) const {
  //
  // Returns true if the particle comes from eta -> pi+ pi- gamma
  //
  if(label<0) return kFALSE;
  Int_t motherLabel = (static_cast<AliAODMCParticle*>(trackArray->At(label)))->GetMother();
  if( motherLabel < 0 || motherLabel >= trackArray->GetEntriesFast() ) return kFALSE;

  AliAODMCParticle* mother = static_cast<AliAODMCParticle*>(trackArray->At(motherLabel));
  if( mother->GetPdgCode() != 223 ) return kFALSE;
  if( IsPiPlPiMiPiZeroDecayAOD(trackArray, mother ) ) return kTRUE;
  return kFALSE;
}

//_____________________________________________________________________________________
Bool_t AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::IsD0PiPlPiMiPiZeroDaughter( Int_t label ) const {
  //
  // Returns true if the particle comes from eta -> pi+ pi- gamma
  //
  if(label<0) return kFALSE;
  Int_t motherLabel = fMCEvent->Particle( label )->GetMother(0);
  if( motherLabel < 0 || motherLabel >= fMCEvent->GetNumberOfTracks() ) return kFALSE;

  TParticle* mother = fMCEvent->Particle( motherLabel );
  if( mother->GetPdgCode() != 421 ) return kFALSE;
  if( IsPiPlPiMiPiZeroDecay( mother ) ) return kTRUE;
  return kFALSE;
}

//_____________________________________________________________________________________
Bool_t AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::IsD0PiPlPiMiPiZeroDaughterAOD(TClonesArray* trackArray, Int_t label ) const {
  //
  // Returns true if the particle comes from eta -> pi+ pi- gamma
  //
  if(label<0) return kFALSE;
  Int_t motherLabel = (static_cast<AliAODMCParticle*>(trackArray->At(label)))->GetMother();
  if( motherLabel < 0 || motherLabel >=  trackArray->GetEntriesFast()  ) return kFALSE;

  AliAODMCParticle* mother = static_cast<AliAODMCParticle*>(trackArray->At(motherLabel));
  if( mother->GetPdgCode() != 421 ) return kFALSE;
  if( IsPiPlPiMiPiZeroDecayAOD(trackArray, mother ) ) return kTRUE;
  return kFALSE;
}


//_____________________________________________________________________________
Bool_t AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::IsPiPlPiMiPiZeroDecay(TParticle *fMCMother) const
{
  if( fMCMother->GetNDaughters() != 3 ) return kFALSE;
  if( !(fMCMother->GetPdgCode() == 221 || fMCMother->GetPdgCode() == 223 || fMCMother->GetPdgCode() == 421)  ) return kFALSE;

  TParticle *posPion = 0x0;
  TParticle *negPion = 0x0;
  TParticle *neutPion    = 0x0;

  for(Int_t index= fMCMother->GetFirstDaughter();index<= fMCMother->GetLastDaughter();index++){
    if(index<0) continue;
    TParticle* temp = (TParticle*)fMCEvent->Particle( index );

    switch( temp->GetPdgCode() ) {
    case 211:
      posPion =  temp;
      break;
    case -211:
      negPion =  temp;
      break;
    case 111:
      neutPion = temp;
      break;
    }
  }
  if( posPion && negPion && neutPion) return kTRUE;

  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::IsPiPlPiMiPiZeroDecayAOD(TClonesArray* trackArray, AliAODMCParticle *fMCMother) const
{
  if( fMCMother->GetNDaughters() != 3 ) return kFALSE;
  if( !(fMCMother->GetPdgCode() == 221 || fMCMother->GetPdgCode() == 223 || fMCMother->GetPdgCode() == 421)  ) return kFALSE;

  AliAODMCParticle *posPion = 0x0;
  AliAODMCParticle *negPion = 0x0;
  AliAODMCParticle *neutPion    = 0x0;

  for(Int_t index= fMCMother->GetDaughterFirst();index<= fMCMother->GetDaughterLast();index++){
    if(index<0) continue;
    AliAODMCParticle* temp =  static_cast<AliAODMCParticle*>(trackArray->At( index ));

    switch( temp->GetPdgCode() ) {
    case 211:
      posPion =  temp;
      break;
    case -211:
      negPion =  temp;
      break;
    case 111:
      neutPion = temp;
      break;
    }
  }
  if( posPion && negPion && neutPion) return kTRUE;

  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::IsPiPlPiMiEtaDecay(TParticle *fMCMother) const
{
  if( fMCMother->GetNDaughters() != 3 ) return kFALSE;
  if( !(fMCMother->GetPdgCode() == 331)  ) return kFALSE;

  TParticle *posPion = 0x0;
  TParticle *negPion = 0x0;
  TParticle *etaMeson    = 0x0;

  for(Int_t index= fMCMother->GetFirstDaughter();index<= fMCMother->GetLastDaughter();index++){
    if(index<0) continue;
    TParticle* temp = (TParticle*)fMCEvent->Particle( index );

    switch( temp->GetPdgCode() ) {
    case 211:
      posPion =  temp;
      break;
    case -211:
      negPion =  temp;
      break;
    case 221:
      etaMeson = temp;
      break;
    }
  }
  if( posPion && negPion && etaMeson) return kTRUE;

  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::IsPiPlPiMiEtaDecayAOD(TClonesArray* trackArray, AliAODMCParticle *fMCMother) const
{
  if( fMCMother->GetNDaughters() != 3 ) return kFALSE;
  if( !(fMCMother->GetPdgCode() == 331)  ) return kFALSE;

  AliAODMCParticle *posPion  = 0x0;
  AliAODMCParticle *negPion  = 0x0;
  AliAODMCParticle *etaMeson = 0x0;

  for(Int_t index= fMCMother->GetDaughterFirst();index<= fMCMother->GetDaughterLast();index++){
    if(index<0) continue;
    AliAODMCParticle * temp =static_cast<AliAODMCParticle*>(trackArray->At( index ));

    switch( temp->GetPdgCode() ) {
    case 211:
      posPion =  temp;
      break;
    case -211:
      negPion =  temp;
      break;
    case 221:
      etaMeson = temp;
      break;
    }
  }
  if( posPion && negPion && etaMeson) return kTRUE;

  return kFALSE;
}

//_____________________________________________________________________________________
Bool_t AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::GammaIsNeutralMesonPiPlPiMiNDMDaughter( Int_t label ) const {
  //
  // Returns true if the particle comes from eta -> pi+ pi- gamma
  //
    if(label<0) return kFALSE;
  Int_t motherLabel = fMCEvent->Particle( label )->GetMother(0);
  if( motherLabel < 0 || motherLabel >= fMCEvent->GetNumberOfTracks() ) return kFALSE;

  TParticle* mother = fMCEvent->Particle( motherLabel );
  if( mother->GetPdgCode() != fPDGCodeNDM ) return kFALSE;
  Int_t grandMotherLabel = mother->GetMother(0);
  if( grandMotherLabel < 0 || grandMotherLabel >= fMCEvent->GetNumberOfTracks() ) return kFALSE;
  TParticle* grandmother = fMCEvent->Particle( grandMotherLabel );

  switch( fSelectedHeavyNeutralMeson ) {
  case 0: // ETA MESON
  case 1: // OMEGA MESON
    if( IsPiPlPiMiPiZeroDecay( grandmother ) ) return kTRUE;
    break;
  case 2: // ETA PRIME MESON
    if( IsPiPlPiMiEtaDecay( grandmother ) ) return kTRUE;
    break;
  case 3: // D0 MESON
    if( IsPiPlPiMiPiZeroDecay(grandmother) )  return kTRUE;
    break;
  default:
    AliError(Form("Heavy neutral meson not correctly selected (only 0,1,2,3 valid)... selected: %d",fSelectedHeavyNeutralMeson));
  }

  return kFALSE;
}

//_____________________________________________________________________________________
Bool_t AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::GammaIsNeutralMesonPiPlPiMiNDMDaughterAOD(TClonesArray* trackArray, Int_t label ) const {
  //
  // Returns true if the particle comes from eta -> pi+ pi- gamma
  //
  if(label<0) return kFALSE;

  Int_t motherLabel =  (static_cast<AliAODMCParticle*>(trackArray->At(label)))->GetMother();
  if( motherLabel < 0 || motherLabel >= trackArray->GetEntriesFast()) return kFALSE;
  AliAODMCParticle* mother = static_cast<AliAODMCParticle*>(trackArray->At(motherLabel));
  if( mother->GetPdgCode() != fPDGCodeNDM ) return kFALSE;
  Int_t grandMotherLabel = mother->GetMother();
  if( grandMotherLabel < 0 || grandMotherLabel >= trackArray->GetEntriesFast()) return kFALSE;
  AliAODMCParticle* grandmother = static_cast<AliAODMCParticle*>(trackArray->At(grandMotherLabel));

  switch( fSelectedHeavyNeutralMeson ) {
  case 0: // ETA MESON
  case 1: // OMEGA MESON
    if( IsPiPlPiMiPiZeroDecayAOD(trackArray, grandmother ) ) return kTRUE;
    break;
  case 2: // ETA PRIME MESON
    if( IsPiPlPiMiEtaDecayAOD(trackArray, grandmother ) ) return kTRUE;
    break;
  case 3: // D0 MESON
    if( IsPiPlPiMiPiZeroDecayAOD(trackArray, grandmother) )  return kTRUE;
    break;
  default:
    AliError(Form("Heavy neutral meson not correctly selected (only 0,1,2,3 valid)... selected: %d",fSelectedHeavyNeutralMeson));
  }

  return kFALSE;
}

//_________________________________________________________________________________
Bool_t AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::CheckVectorForDoubleCount(vector<Int_t> &vec, Int_t tobechecked)
{
  if(tobechecked > -1)
  {
    vector<Int_t>::iterator it;
    it = find (vec.begin(), vec.end(), tobechecked);
    if (it != vec.end()) return true;
    else{
      vec.push_back(tobechecked);
      return false;
    }
  }
  return false;
}
//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::RelabelAODPhotonCandidates(Bool_t mode){

  // Relabeling For AOD Event
  // ESDiD -> AODiD
  // MCLabel -> AODMCLabel

  if(mode){
    fMCEventPos.Set(fReaderGammas->GetEntries());
    fMCEventNeg.Set(fReaderGammas->GetEntries());
    fESDArrayPos.Set(fReaderGammas->GetEntries());
    fESDArrayNeg.Set(fReaderGammas->GetEntries());
  }

  for(Int_t iGamma = 0;iGamma<fReaderGammas->GetEntries();iGamma++){
    AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(iGamma);
    if(!PhotonCandidate) continue;
    if(!mode){// Back to ESD Labels
      PhotonCandidate->SetMCLabelPositive(fMCEventPos[iGamma]);
      PhotonCandidate->SetMCLabelNegative(fMCEventNeg[iGamma]);
      PhotonCandidate->SetLabelPositive(fESDArrayPos[iGamma]);
      PhotonCandidate->SetLabelNegative(fESDArrayNeg[iGamma]);
      continue;
    }
    fMCEventPos[iGamma] =  PhotonCandidate->GetMCLabelPositive();
    fMCEventNeg[iGamma] =  PhotonCandidate->GetMCLabelNegative();
    fESDArrayPos[iGamma] = PhotonCandidate->GetTrackLabelPositive();
    fESDArrayNeg[iGamma] = PhotonCandidate->GetTrackLabelNegative();

    Bool_t AODLabelPos = kFALSE;
    Bool_t AODLabelNeg = kFALSE;

    for(Int_t i = 0; i<fInputEvent->GetNumberOfTracks();i++){
      AliAODTrack *tempDaughter = static_cast<AliAODTrack*>(fInputEvent->GetTrack(i));
      if(!AODLabelPos){
        if( tempDaughter->GetID() == PhotonCandidate->GetTrackLabelPositive() ){
        PhotonCandidate->SetMCLabelPositive(TMath::Abs(tempDaughter->GetLabel()));
        PhotonCandidate->SetLabelPositive(i);
        AODLabelPos = kTRUE;
        }
      }
      if(!AODLabelNeg){
        if( tempDaughter->GetID() == PhotonCandidate->GetTrackLabelNegative()){
        PhotonCandidate->SetMCLabelNegative(TMath::Abs(tempDaughter->GetLabel()));
        PhotonCandidate->SetLabelNegative(i);
        AODLabelNeg = kTRUE;
        }
      }
      if(AODLabelNeg && AODLabelPos){
        break;
      }
    }
    if(!AODLabelPos || !AODLabelNeg){
      cout<<"WARNING!!! AOD TRACKS NOT FOUND FOR"<<endl;
      if(!AODLabelNeg){
        PhotonCandidate->SetMCLabelNegative(-999999);
        PhotonCandidate->SetLabelNegative(-999999);
      }
      if(!AODLabelPos){
        PhotonCandidate->SetMCLabelPositive(-999999);
        PhotonCandidate->SetLabelPositive(-999999);
      }
    }
  }
}
//________________________________________________________________________
AliExternalTrackParam* AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::GetConstrainedParameterAOD(const AliAODTrack* aodTr, const AliAODVertex* vtx, double bz)
{
  double chi2;
  AliExternalTrackParam* par = new AliExternalTrackParam();
  par->CopyFromVTrack(aodTr);
  double dz[2];
  if (!par->PropagateToDCA(vtx,bz,999.,dz,0)) {
    delete par;
    chi2 = 1e9;
    return 0;
  }
  Double_t covar[6]; vtx->GetCovarianceMatrix(covar);
  Double_t p[2]= { par->GetParameter()[0]-dz[0], par->GetParameter()[1]-dz[1]};
  Double_t c[3]= { covar[2],0.,covar[5] };
  chi2 = par->GetPredictedChi2(p,c);
  if (chi2>1e9 || !par->Update(p,c)) {
    AliFatal(Form("Propagation failed with chi2 = %f",chi2));
    delete par;
    return 0;
  }
  return par;
}

//_____________________________________________________________________________
Double32_t AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::CalculateP2(Double_t xyz[3],Double_t pxpypz[3])
{
  //
  // create external track parameters from the global parameters
  // x,y,z,px,py,pz and their 6x6 covariance matrix
  // A.Dainese 10.10.08

  // Calculate alpha: the rotation angle of the corresponding local system.
  //
  // For global radial position inside the beam pipe, alpha is the
  // azimuthal angle of the momentum projected on (x,y).
  //
  // For global radial position outside the ITS, alpha is the
  // azimuthal angle of the centre of the TPC sector in which the point
  // xyz lies
  //
  const double kSafe = 1e-5;
  Double32_t p2;
  Double32_t fAlpha;
  Double_t radPos2 = xyz[0]*xyz[0]+xyz[1]*xyz[1];
  Double_t radMax  = 45.; // approximately ITS outer radius
  if (radPos2 < radMax*radMax) { // inside the ITS
     fAlpha = TMath::ATan2(pxpypz[1],pxpypz[0]);
  } else { // outside the ITS
     Float_t phiPos = TMath::Pi()+TMath::ATan2(-xyz[1], -xyz[0]);
     fAlpha =
     TMath::DegToRad()*(20*((((Int_t)(phiPos*TMath::RadToDeg()))/20))+10);
  }
  //
  Double_t cs=TMath::Cos(fAlpha), sn=TMath::Sin(fAlpha);
  // protection:  avoid alpha being too close to 0 or +-pi/2
  if (TMath::Abs(sn)<2*kSafe) {
    if (fAlpha>0) fAlpha += fAlpha< TMath::Pi()/2. ?  2*kSafe : -2*kSafe;
    else          fAlpha += fAlpha>-TMath::Pi()/2. ? -2*kSafe :  2*kSafe;
    cs=TMath::Cos(fAlpha);
    sn=TMath::Sin(fAlpha);
  }
  else if (TMath::Abs(cs)<2*kSafe) {
    if (fAlpha>0) fAlpha += fAlpha> TMath::Pi()/2. ? 2*kSafe : -2*kSafe;
    else          fAlpha += fAlpha>-TMath::Pi()/2. ? 2*kSafe : -2*kSafe;
    cs=TMath::Cos(fAlpha);
    sn=TMath::Sin(fAlpha);
  }
  // Get the vertex of origin and the momentum
  TVector3 ver(xyz[0],xyz[1],xyz[2]);
  TVector3 mom(pxpypz[0],pxpypz[1],pxpypz[2]);
  //
  // Rotate to the local coordinate system
  ver.RotateZ(-fAlpha);
  mom.RotateZ(-fAlpha);

  //
  // x of the reference plane

   p2= TMath::Sin(mom.Phi());
  //
  if      (TMath::Abs( 1-p2) < kSafe) p2= 1.- kSafe; //Protection
  else if (TMath::Abs(-1-p2) < kSafe) p2=-1.+ kSafe; //Protection
  return p2;
}
