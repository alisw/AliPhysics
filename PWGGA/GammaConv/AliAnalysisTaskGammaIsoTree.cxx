/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Florian Jonas                                                 *
 * Version 1.0                                                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliAnalysisTaskGammaIsoTree.h"
#include "TChain.h"
#include "TRandom.h"
#include "AliAnalysisManager.h"
#include "TParticle.h"
#include "TVectorF.h"
#include "AliPIDResponse.h"
#include "TFile.h"
#include "AliESDtrackCuts.h"
#include "AliAODMCParticle.h"
#include "AliAODConversionPhoton.h"
#include "AliAODMCHeader.h"
#include "AliAODEvent.h"
#include "AliMultSelection.h"
#include "AliAODCaloCluster.h"
#include "AliAODTrack.h"
#include "AliVTrack.h"
#include "AliEMCALRecoUtilsBase.h"
#include "AliAODConversionMother.h"
#include "TObjectTable.h"
#include <chrono>

ClassImp(AliAnalysisTaskGammaIsoTree)
//________________________________________________________________________
AliAnalysisTaskGammaIsoTree::AliAnalysisTaskGammaIsoTree() : AliAnalysisTaskSE(),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fWeightJetJetMC(1),
  fOutputList(NULL),
  fConvFolderRec(NULL),
  fConvFolderTrue(NULL),
  fCaloFolderRec(NULL),
  fCaloFolderTrue(NULL),
  fGeneralFolder(NULL),
  fQAFolder(NULL),
  fGeneratorFolder(NULL),
  fAnalysisTree(NULL),
  fIsMC(0),
  fIsHeavyIon(0),
  fV0Reader(NULL),
  fV0ReaderName(""),
  fReaderGammas(NULL),
  fConversionCandidates(0),
  fClusterEMCalCandidates(0),
  fClusterEMCalCandidatesIsolation(0),
  fClusterEMCalCandidatesTagging(0),
  fClusterPHOSCandidates(0),
  fTracks(0),
  fMCParticles(0),
  fAODMCTrackArray(NULL),
  fGeomEMCAL(NULL),
  fCorrTaskSetting(""),
  fEventCuts(NULL),
  fClusterCutsEMC(NULL),
  fClusterCutsEMCTrackMatching(NULL),
  fClusterCutsIsolationEMC(NULL),
  fClusterCutsTaggingEMC(NULL),
  fClusterCutsPHOS(NULL),
  fConvCuts(NULL),
  fCaloUtils(NULL),
  fMCAnalysisUtils(NULL),
  fMinClsTPC(0),
  fChi2PerClsTPC(9999),
  fMinClsITS(0),
  fEtaCut(9999),
  fPtCut(0),
  fYMCCut(9999),
  fMatchingParamsPhi(),
  fMatchingParamsEta(),
  fMatchingEOverP(),
  fDoBackgroundTrackMatching(kFALSE),
  fDoOwnTrackMatching(kFALSE),
  fDoTrackIsolation(kFALSE),
  fTrackIsolationR(),
  fTrackIsolationE(),
  fDoNeutralIsolation(kFALSE),
  fNeutralIsolationR(),
  fNeutralIsolationE(),
  fDoCellIsolation(kTRUE),
  fDoTagging(kFALSE),
  fPi0TaggingWindow(),
  fEtaTaggingWindow(),
  fSaveConversions(kTRUE),
  fSaveEMCClusters(kTRUE),
  fSavePHOSClusters(kTRUE),
  fSaveTracks(kTRUE),
  fUseHistograms(kFALSE),
  fUseTree(0),
  fRecPtCut(0),
  fGenPtCut(0),
  fHistoNEvents(NULL),
  fHistoNEventsWOWeight(NULL),
  fHistoChargedIso(NULL),
  fHistoTaggingPCMPCM(NULL),
  fHistoTaggingPCMEMC(NULL),
  fHistoTaggingEMCPCM(NULL),
  fHistoTaggingEMCEMC(NULL),
  fConvPt(NULL),
  fConvPtBeforeAcc(NULL),
  fConvPtTaggedCalo(NULL),
  fConvPtTaggedAsDecayCalo(NULL),
  fConvIsoCharged(),
  fConvIsoNeutral(),
  fConvIsoFull(),
  fConvIsoCell(),
  fConvIsoCorr(),
  fConvTruePt(NULL),
  fConvTruePtPrimary(NULL),
  fConvTruePtDecay(NULL),
  fConvTruePtDecayFoundOtherInCluster(NULL),
  fConvTruePtDecayOtherInAcc(NULL),
  fConvTruePtDecayOtherInAccAboveMinEnergy(NULL),
  fConvTruePtTaggedCalo(NULL),
  fConvTruePtTaggedAsDecayCalo(NULL),
  fConvTrueRecPt(NULL),
  fConvTrueRecPtPrimary(NULL),
  fConvTrueRecPtDecay(NULL),
  fConvTrueRecPtDecayFoundOtherInCluster(NULL),
  fConvTrueRecPtDecayOtherInAcc(NULL),
  fConvTrueRecPtDecayOtherInAccAboveMinEnergy(NULL),
  fConvTrueRecPtTaggedCalo(NULL),
  fConvTrueRecPtTaggedAsDecayCalo(NULL),
  fConvTrueIsoCharged(),
  fConvTrueIsoNeutral(),
  fConvTrueIsoFull(),
  fConvTrueIsoCorr(),
  fConvTrueIsoCell(),
  fConvTrueIsoCharged_FromDecay(),
  fConvTrueIsoNeutral_FromDecay(),
  fConvTrueIsoFull_FromDecay(),
  fConvTrueIsoCell_FromDecay(),
  fConvTrueIsoCharged_FromDirect(),
  fConvTrueIsoNeutral_FromDirect(),
  fConvTrueIsoFull_FromDirect(),
  fConvTrueIsoCell_FromDirect(),
  fConvPtIsoCharged(),
  fConvPtIsoNeutral(),
  fConvPtIsoFull(),
  fConvPtIsoCell(),
  fConvPtTaggedCaloIsoCharged(),
  fConvPtTaggedCaloIsoNeutral(),
  fConvPtTaggedCaloIsoFull(),
  fConvPtTaggedCaloIsoCell(),
  fConvTruePtIsoCharged(),
  fConvTruePtIsoNeutral(),
  fConvTruePtIsoFull(),
  fConvTruePtIsoCell(),
  fConvTruePtTaggedCaloIsoCharged(),
  fConvTruePtTaggedCaloIsoNeutral(),
  fConvTruePtTaggedCaloIsoFull(),
  fConvTruePtTaggedCaloIsoCell(),
   fConvTruePtIsoChargedFromDirect(),
  fConvTruePtIsoNeutralFromDirect(),
  fConvTruePtIsoFullFromDirect(),
  fConvTruePtIsoCellFromDirect(),
  fConvTruePtTaggedCaloIsoChargedFromDirect(),
  fConvTruePtTaggedCaloIsoNeutralFromDirect(),
  fConvTruePtTaggedCaloIsoFullFromDirect(),
  fConvTruePtTaggedCaloIsoCellFromDirect(),
  fConvTrueRecPtIsoCharged(),
  fConvTrueRecPtIsoNeutral(),
  fConvTrueRecPtIsoFull(),
  fConvTrueRecPtIsoCell(),
  fConvTrueRecPtIsoChargedFromDirect(),
  fConvTrueRecPtIsoNeutralFromDirect(),
  fConvTrueRecPtIsoFullFromDirect(),
  fConvTrueRecPtIsoCellFromDirect(),
  fConvTrueRecPtTaggedCaloIsoCharged(),
  fConvTrueRecPtTaggedCaloIsoNeutral(),
  fConvTrueRecPtTaggedCaloIsoFull(),
  fConvTrueRecPtTaggedCaloIsoCell(),
  fConvTruePtMCIsoCharged(),
  fConvTruePtMCIsoNeutral(),
  fConvTruePtMCIsoFull(),
  fConvTruePtMCIsoCell(),
  fConvTruePtTaggedCaloMCIsoCharged(),
  fConvTruePtTaggedCaloMCIsoNeutral(),
  fConvTruePtTaggedCaloMCIsoFull(),
  fConvTruePtTaggedCaloMCIsoCell(),
  fConvTrueRecPtMCIsoCharged(),
  fConvTrueRecPtMCIsoNeutral(),
  fConvTrueRecPtMCIsoFull(),
  fConvTrueRecPtMCIsoCell(),
  fConvTrueRecPtTaggedCaloMCIsoCharged(),
  fConvTrueRecPtTaggedCaloMCIsoNeutral(),
  fConvTrueRecPtTaggedCaloMCIsoFull(),
  fConvTrueRecPtTaggedCaloMCIsoCell(),
  fConvInvMass(NULL),
  fConvInvMassChargedIsolated(),
  fConvInvMassAntiChargedIsolated(),
  fConvInvMassNeutralIsolated(),
  fConvInvMassAntiNeutralIsolated(),
  fConvInvMassCellIsolated(),
  fConvInvMassAntiCellIsolated(),
  fConvInvMassFullIsolated(),
  fConvInvMassAntiFullIsolated(),
  fConvTrueInvMass(NULL),
  fConvTrueInvMass_FromDecay(NULL),
  fConvTrueInvMass_FromDirect(NULL),
  fCaloPt(NULL),
  fCaloPtBeforeAcc(NULL),
  fCaloE(NULL),
  fCaloPtTaggedCalo(NULL),
  fCaloPtTaggedAsDecayCalo(NULL),
  fCaloIsoCharged(),
  fCaloIsoNeutral(),
  fCaloIsoFull(),
  fCaloIsoCell(),
  fCaloIsoCorr(),
  fCaloRho(NULL),
  fCaloRhoTimesArea(),
  fCaloRhoTimesAreaLeft(),
  fCaloRhoTimesAreaRight(),
  fCaloRhoTimesAreaBack(),
  fCaloTruePt(NULL),
  fCaloTruePtNotProper(NULL),
  fCaloTruePtNoIsPrimary(NULL),
  fCaloTrueRecPtNoIsPrimary(NULL),
  fCaloTrueWithoutConvPt(NULL),
  fCaloTruePtPrimary(NULL),
  fCaloTruePtDecay(NULL),
  fCaloTruePtDecayFoundOtherInCluster(NULL),
  fCaloTruePtDecayOtherInAcc(NULL),
  fCaloTruePtDecayOtherInAccAboveMinEnergy(NULL),
  fCaloTruePtTaggedCalo(NULL),
  fCaloTruePtTaggedAsDecayCalo(NULL),
  fCaloTrueRecPt(NULL),
  fCaloTrueRecPtvsTruePt(NULL),
  fCaloTrueRecPtvsTruePtNotProper(NULL),
  fCaloTrueWithoutConvRecPt(NULL),
  fCaloTrueRecPtPrimary(NULL),
  fCaloTrueRecPtDecay(NULL),
  fCaloTrueRecPtDecayFoundOtherInCluster(NULL),
  fCaloTrueRecPtDecayOtherInAcc(NULL),
  fCaloTrueRecPtDecayOtherInAccAboveMinEnergy(NULL),
  fCaloTrueRecPtTaggedCalo(NULL),
  fCaloTrueRecPtTaggedAsDecayCalo(NULL),
  fCaloTrueIsoCharged(),
  fCaloTrueIsoMCCharged(),
  fCaloTrueIsoNeutral(),
  fCaloTrueIsoFull(),
  fCaloTrueIsoCorr(),
  fCaloTrueIsoCell(),
  fCaloTrueIsoCharged_FromDecay(),
  fCaloTrueIsoNeutral_FromDecay(),
  fCaloTrueIsoFull_FromDecay(),
  fCaloTrueIsoCell_FromDecay(),
  fCaloTrueIsoCharged_FromDirect(),
  fCaloTrueIsoNeutral_FromDirect(),
  fCaloTrueIsoFull_FromDirect(),
  fCaloTrueIsoCell_FromDirect(),
  fCaloTrueSignalPtClusterCuts(NULL),
  fCaloTrueSignalPtClusterCutsOnlyPhoton(NULL),
  fCaloTrueSignalPtClusterCutsOnlyConv(NULL),
  fCaloTrueSignalPtClusterCutsM02Cuts(NULL),
  fCaloTrueSignalPtClusterCutsM02CutsIsoCuts(),
  fCaloTrueSignalPtClusterCutsM02CutsIsoCutsOnlyPhoton(),
  fCaloTrueSignalPtClusterCutsM02CutsIsoCutsOnlyConv(),
  fCaloTrueSignalRecPtClusterCuts(NULL),
  fCaloTrueSignalRecPtClusterCutsOnlyPhoton(NULL),
  fCaloTrueSignalRecPtClusterCutsOnlyConv(NULL),
  fCaloTrueSignalRecPtClusterCutsM02Cuts(NULL),
  fCaloTrueSignalRecPtClusterCutsM02CutsIsoCuts(),
  fCaloTrueSignalRecPtClusterCutsM02CutsIsoCutsOnlyPhoton(),
  fCaloTrueSignalRecPtClusterCutsM02CutsIsoCutsOnlyConv(),
  fCaloPtIsoCharged(),
  fCaloPtIsoNeutral(),
  fCaloPtIsoFull(),
  fCaloPtIsoCell(),
  fCaloPtTaggedCaloIsoCharged(),
  fCaloPtTaggedCaloIsoNeutral(),
  fCaloPtTaggedCaloIsoFull(),
  fCaloPtTaggedCaloIsoCell(),
  fCaloTruePtIsoCharged(),
  fCaloTruePtIsoNeutral(),
  fCaloTruePtIsoFull(),
  fCaloTruePtIsoCell(),
  fCaloTruePtTaggedCaloIsoCharged(),
  fCaloTruePtTaggedCaloIsoNeutral(),
  fCaloTruePtTaggedCaloIsoFull(),
  fCaloTruePtTaggedCaloIsoCell(),
  fCaloTruePtIsoChargedFromDirect(),
  fCaloTruePtIsoChargedAndMCIsoChargedFromDirect(),
  fCaloTruePtIsoNeutralFromDirect(),
  fCaloTruePtIsoFullFromDirect(),
  fCaloTruePtIsoCellFromDirect(),
  fCaloTruePtTaggedCaloIsoChargedFromDirect(),
  fCaloTruePtTaggedCaloIsoNeutralFromDirect(),
  fCaloTruePtTaggedCaloIsoFullFromDirect(),
  fCaloTruePtTaggedCaloIsoCellFromDirect(),
  fCaloTrueRecPtIsoCharged(),
  fCaloTrueRecPtIsoNeutral(),
  fCaloTrueRecPtIsoFull(),
  fCaloTrueRecPtIsoCell(),
  fCaloTrueRecPtIsoChargedFromDirect(),
  fCaloTrueRecPtIsoChargedAndMCIsoChargedFromDirect(),
  fCaloTrueRecPtIsoNeutralFromDirect(),
  fCaloTrueRecPtIsoFullFromDirect(),
  fCaloTrueRecPtIsoCellFromDirect(),
  fCaloTrueRecPtTaggedCaloIsoCharged(),
  fCaloTrueRecPtTaggedCaloIsoNeutral(),
  fCaloTrueRecPtTaggedCaloIsoFull(),
  fCaloTrueRecPtTaggedCaloIsoCell(),
  fCaloTruePtMCIsoCharged(),
  fCaloTruePtMCIsoNeutral(),
  fCaloTruePtMCIsoFull(),
  fCaloTruePtMCIsoCell(),
  fCaloTruePtTaggedCaloMCIsoCharged(),
  fCaloTruePtTaggedCaloMCIsoNeutral(),
  fCaloTruePtTaggedCaloMCIsoFull(),
  fCaloTruePtTaggedCaloMCIsoCell(),
  fCaloTrueRecPtMCIsoCharged(),
  fCaloTrueRecPtMCIsoNeutral(),
  fCaloTrueRecPtMCIsoFull(),
  fCaloTrueRecPtMCIsoCell(),
  fCaloTrueRecPtTaggedCaloMCIsoCharged(),
  fCaloTrueRecPtTaggedCaloMCIsoNeutral(),
  fCaloTrueRecPtTaggedCaloMCIsoFull(),
  fCaloTrueRecPtTaggedCaloMCIsoCell(),
  fCaloInvMass(),
  fCaloInvMassChargedIsolated(),
  fCaloInvMassAntiChargedIsolated(),
  fCaloInvMassNeutralIsolated(),
  fCaloInvMassAntiNeutralIsolated(),
  fCaloInvMassCellIsolated(),
  fCaloInvMassAntiCellIsolated(),
  fCaloInvMassFullIsolated(),
  fCaloInvMassAntiFullIsolated(),
  fCaloTrueInvMass(),
  fCaloTrueInvMass_FromDecay(),
  fCaloTrueInvMass_FromDirect(),
  fCaloM02(),
  fCaloM02ChargedIsolated(),
  fCaloM02ChargedIsolated_FromDirect(),
  fCaloM02ChargedIsolated_FromDirectConvOnly(),
  fCaloM02ChargedIsolated_FromDirectNormOnly(),
  fCaloM02AntiChargedIsolated(),
  fCaloM02NeutralIsolated(),
  fCaloM02AntiNeutralIsolated(),
  fCaloM02CellIsolated(),
  fCaloM02AntiCellIsolated(),
  fCaloM02FullIsolated(),
  fCaloM02AntiFullIsolated(),
  fCaloM02vsIsoVsPt(),
  fCaloTrueM02(),
  fCaloTrueM02_FromDecay(),
  fCaloTrueM02_FromDirect(),
  fHistoMCHeaders(NULL),
  fGenPhotonPt(NULL),
  fGenPhotonPt_FromDecay(NULL),
  fGenPhotonPt_FromDirect(NULL),
  fGenPhotonPtInEMCalAcc(NULL),
  fGenPhotonPtInEMCalAcc_FromDecay(NULL),
  fGenPhotonPtInEMCalAcc_FromDirect(NULL),
  fGenPhotonPtInEMCalAcc_FromDirect_OnlyPhoton(NULL),
  fGenPhotonPtInEMCalAcc_FromDirect_OnlyConv(NULL),
  fGenPhotonPtInEMCalAccChargedMCIso_FromDirect(),
  fGenPhotonPtInEMCalAccChargedMCIso_FromDirect_OnlyPhoton(),
  fGenPhotonPtInEMCalAccChargedMCIso_FromDirect_OnlyConv(),
  fGenPhotonPtFoundNormCluster(NULL),
  fGenPhotonPtFoundTaggingCluster(NULL),
  fGenPhotonPtFoundIsoCluster(NULL),
  fGenPhotonEFoundNoClusterVsCellE(NULL),
  fGenPi0Pt(NULL),
  fGenPi0PtInEMCalAcc(NULL),
  fGenPi0PtInEMCalAcc_BothGammaInEMCal(NULL),
  fGenPi0PtInEMCalAcc_BothGammaInClusters(NULL),
  fRhoOutName("Rho"),
  fRhoOutNameMC("RhoMC"),
  fTreeBuffSize(60*1024*1024),
  fMemCountAOD(0),
  fTrackMatcherRunningMode(0),
  fAntiIsolationE(),
  fMinM02(0),
  fMaxM02(9999),
  fChargedRho(0),
  fChargedRhoMC(0),
  fChargedRhoTimesArea(),
  fExclusionRadius(0.4),
  fDebug(0),
  fOutlierJetReader(0),
  fIsFromDesiredHeader(kTRUE),
  fIsOverlappingWithOtherHeader(kFALSE),
  fAllowOverlapHeaders(kTRUE),
  fBuffer_EventRho(0),
  fBuffer_EventRhoMC(0),
  fBuffer_EventWeight(0),
  fBuffer_EventXsection(0),
  fBuffer_EventNtrials(0),
  fBuffer_EventIsTriggered(0),
  fBuffer_ClusterE(0), 
  fBuffer_ClusterPx(0), 
  fBuffer_ClusterPy(0), 
  fBuffer_ClusterPz(0), 
  fBuffer_ClusterM02(0), 
  fBuffer_ClusterM02Recalc(0), 
  fBuffer_ClusterM20(0), 
  fBuffer_ClusterV1SplitMass(0), 
  fBuffer_ClusterNLM(0), 
  fBuffer_ClusterSM(0), 
  fBuffer_ClusterEFrac(0), 
  fBuffer_ClusterIsoCharged1(0), 
  fBuffer_ClusterIsoCharged2(0), 
  fBuffer_ClusterIsoCharged3(0), 
  fBuffer_ClusterIsoBckLeft(0), 
  fBuffer_ClusterMatchTrackdEta(0), 
  fBuffer_ClusterMatchTrackdPhi(0), 
  fBuffer_ClusterMatchTrackP(0), 
  fBuffer_ClusterMatchTrackPt(0), 
  fBuffer_ClusterMatchTrackIsConv(0), 
  fBuffer_TrueClusterE(0), 
  fBuffer_TrueClusterPx(0), 
  fBuffer_TrueClusterPy(0), 
  fBuffer_TrueClusterPz(0), 
  fBuffer_TrueClusterLeadingEFrac(0), 
  fBuffer_TrueClusterMCIsoCharged1(0), 
  fBuffer_TrueClusterMCIsoCharged2(0), 
  fBuffer_TrueClusterMCIsoCharged3(0), 
  fBuffer_TrueClusterMCIsoBckLeft(0), 
  fBuffer_TrueClusterMCTag(0), 
  fBuffer_TrueClusterIsConv(0),
  fBuffer_GenPhotonE(0),
  fBuffer_GenPhotonPx(0),
  fBuffer_GenPhotonPy(0),
  fBuffer_GenPhotonPz(0),
  fBuffer_GenPhotonMCIsoCharged1(0),
  fBuffer_GenPhotonMCIsoCharged2(0),
  fBuffer_GenPhotonMCIsoCharged3(0),
  fBuffer_GenPhotonMCIsoBckLeft(0),
  fBuffer_GenPhotonIsConv(0),
  fBuffer_GenPhotonMCTag(0)
{

  SetEtaMatching(0.010,4.07,-2.5);
  SetPhiMatching(0.015,3.65,3.65);
  SetEOverP(1.75);

  SetPi0TaggingWindow(0.120,0.145);
  SetEtaTaggingWindow(0.5,0.6);

  SetAntiIsolationE(5.,10);

  fChargedRhoTimesArea[0] = 0.;
  fChargedRhoTimesArea[1] = 0.;
  fChargedRhoTimesArea[2] = 0.;
  fChargedRhoTimesArea[3] = 0.;
  fChargedRhoTimesArea[4] = 0.;

}

AliAnalysisTaskGammaIsoTree::AliAnalysisTaskGammaIsoTree(const char *name) : AliAnalysisTaskSE(name),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fWeightJetJetMC(1),
  fOutputList(NULL),
  fConvFolderRec(NULL),
  fConvFolderTrue(NULL),
  fCaloFolderRec(NULL),
  fCaloFolderTrue(NULL),
  fGeneralFolder(NULL),
  fQAFolder(NULL),
  fGeneratorFolder(NULL),
  fAnalysisTree(NULL),
  fIsMC(0),
  fIsHeavyIon(0),
  fV0Reader(NULL),
  fV0ReaderName(""),
  fReaderGammas(NULL),
  fConversionCandidates(0),
  fClusterEMCalCandidates(0),
  fClusterEMCalCandidatesIsolation(0),
  fClusterEMCalCandidatesTagging(0),
  fClusterPHOSCandidates(0),
  fTracks(0),
  fMCParticles(0),
  fAODMCTrackArray(NULL),
  fGeomEMCAL(NULL),
  fCorrTaskSetting(""),
  fEventCuts(NULL),
  fClusterCutsEMC(NULL),
  fClusterCutsEMCTrackMatching(NULL),
  fClusterCutsIsolationEMC(NULL),
  fClusterCutsTaggingEMC(NULL),
  fClusterCutsPHOS(NULL),
  fConvCuts(NULL),
  fCaloUtils(NULL),
  fMCAnalysisUtils(NULL),
  fMinClsTPC(0),
  fChi2PerClsTPC(9999),
  fMinClsITS(0),
  fEtaCut(9999),
  fPtCut(0),
  fYMCCut(99999),
  fMatchingParamsPhi(),
  fMatchingParamsEta(),
  fMatchingEOverP(),
  fDoBackgroundTrackMatching(kFALSE),
  fDoOwnTrackMatching(kFALSE),
  fDoTrackIsolation(kFALSE),
  fTrackIsolationR(),
  fTrackIsolationE(),
  fDoNeutralIsolation(kFALSE),
  fNeutralIsolationR(),
  fNeutralIsolationE(),
  fDoCellIsolation(kTRUE),
  fDoTagging(kFALSE),
  fPi0TaggingWindow(),
  fEtaTaggingWindow(),
  fSaveConversions(kTRUE),
  fSaveEMCClusters(kTRUE),
  fSavePHOSClusters(kTRUE),
  fSaveTracks(kTRUE),
  fUseHistograms(kFALSE),
  fUseTree(0),
  fRecPtCut(0),
  fGenPtCut(0),
  fHistoNEvents(NULL),
  fHistoNEventsWOWeight(NULL),
  fHistoChargedIso(NULL),
  fHistoTaggingPCMPCM(NULL),
  fHistoTaggingPCMEMC(NULL),
  fHistoTaggingEMCPCM(NULL),
  fHistoTaggingEMCEMC(NULL),
  fConvPt(NULL),
  fConvPtBeforeAcc(NULL),
  fConvPtTaggedCalo(NULL),
  fConvPtTaggedAsDecayCalo(NULL),
  fConvIsoCharged(),
  fConvIsoNeutral(),
  fConvIsoFull(),
  fConvIsoCell(),
  fConvIsoCorr(),
  fConvTruePt(NULL),
  fConvTruePtPrimary(NULL),
  fConvTruePtDecay(NULL),
  fConvTruePtDecayFoundOtherInCluster(NULL),
  fConvTruePtDecayOtherInAcc(NULL),
  fConvTruePtDecayOtherInAccAboveMinEnergy(NULL),
  fConvTruePtTaggedCalo(NULL),
  fConvTruePtTaggedAsDecayCalo(NULL),
  fConvTrueRecPt(NULL),
  fConvTrueRecPtPrimary(NULL),
  fConvTrueRecPtDecay(NULL),
  fConvTrueRecPtDecayFoundOtherInCluster(NULL),
  fConvTrueRecPtDecayOtherInAcc(NULL),
  fConvTrueRecPtDecayOtherInAccAboveMinEnergy(NULL),
  fConvTrueRecPtTaggedCalo(NULL),
  fConvTrueRecPtTaggedAsDecayCalo(NULL),
  fConvTrueIsoCharged(),
  fConvTrueIsoNeutral(),
  fConvTrueIsoFull(),
  fConvTrueIsoCorr(),
  fConvTrueIsoCell(),
  fConvTrueIsoCharged_FromDecay(),
  fConvTrueIsoNeutral_FromDecay(),
  fConvTrueIsoFull_FromDecay(),
  fConvTrueIsoCell_FromDecay(),
  fConvTrueIsoCharged_FromDirect(),
  fConvTrueIsoNeutral_FromDirect(),
  fConvTrueIsoFull_FromDirect(),
  fConvTrueIsoCell_FromDirect(),
  fConvPtIsoCharged(),
  fConvPtIsoNeutral(),
  fConvPtIsoFull(),
  fConvPtIsoCell(),
  fConvPtTaggedCaloIsoCharged(),
  fConvPtTaggedCaloIsoNeutral(),
  fConvPtTaggedCaloIsoFull(),
  fConvPtTaggedCaloIsoCell(),
  fConvTruePtIsoCharged(),
  fConvTruePtIsoNeutral(),
  fConvTruePtIsoFull(),
  fConvTruePtIsoCell(),
  fConvTruePtTaggedCaloIsoCharged(),
  fConvTruePtTaggedCaloIsoNeutral(),
  fConvTruePtTaggedCaloIsoFull(),
  fConvTruePtTaggedCaloIsoCell(),
  fConvTruePtIsoChargedFromDirect(),
  fConvTruePtIsoNeutralFromDirect(),
  fConvTruePtIsoFullFromDirect(),
  fConvTruePtIsoCellFromDirect(),
  fConvTruePtTaggedCaloIsoChargedFromDirect(),
  fConvTruePtTaggedCaloIsoNeutralFromDirect(),
  fConvTruePtTaggedCaloIsoFullFromDirect(),
  fConvTruePtTaggedCaloIsoCellFromDirect(),
  fConvTrueRecPtIsoCharged(),
  fConvTrueRecPtIsoNeutral(),
  fConvTrueRecPtIsoFull(),
  fConvTrueRecPtIsoCell(),
  fConvTrueRecPtIsoChargedFromDirect(),
  fConvTrueRecPtIsoNeutralFromDirect(),
  fConvTrueRecPtIsoFullFromDirect(),
  fConvTrueRecPtIsoCellFromDirect(),
  fConvTrueRecPtTaggedCaloIsoCharged(),
  fConvTrueRecPtTaggedCaloIsoNeutral(),
  fConvTrueRecPtTaggedCaloIsoFull(),
  fConvTrueRecPtTaggedCaloIsoCell(),
  fConvTruePtMCIsoCharged(),
  fConvTruePtMCIsoNeutral(),
  fConvTruePtMCIsoFull(),
  fConvTruePtMCIsoCell(),
  fConvTruePtTaggedCaloMCIsoCharged(),
  fConvTruePtTaggedCaloMCIsoNeutral(),
  fConvTruePtTaggedCaloMCIsoFull(),
  fConvTruePtTaggedCaloMCIsoCell(),
  fConvTrueRecPtMCIsoCharged(),
  fConvTrueRecPtMCIsoNeutral(),
  fConvTrueRecPtMCIsoFull(),
  fConvTrueRecPtMCIsoCell(),
  fConvTrueRecPtTaggedCaloMCIsoCharged(),
  fConvTrueRecPtTaggedCaloMCIsoNeutral(),
  fConvTrueRecPtTaggedCaloMCIsoFull(),
  fConvTrueRecPtTaggedCaloMCIsoCell(),
  fConvInvMass(NULL),
  fConvInvMassChargedIsolated(),
  fConvInvMassAntiChargedIsolated(),
  fConvInvMassNeutralIsolated(),
  fConvInvMassAntiNeutralIsolated(),
  fConvInvMassCellIsolated(),
  fConvInvMassAntiCellIsolated(),
  fConvInvMassFullIsolated(),
  fConvInvMassAntiFullIsolated(),
  fConvTrueInvMass(NULL),
  fConvTrueInvMass_FromDecay(NULL),
  fConvTrueInvMass_FromDirect(NULL),
  fCaloPt(NULL),
  fCaloPtBeforeAcc(NULL),
  fCaloE(NULL),
  fCaloPtTaggedCalo(NULL),
  fCaloPtTaggedAsDecayCalo(NULL),
  fCaloIsoCharged(),
  fCaloIsoRawCharged(),
  fCaloIsoNeutral(),
  fCaloIsoFull(),
  fCaloIsoCell(),
  fCaloIsoCorr(),
  fCaloRho(NULL),
  fCaloRhoTimesArea(),
  fCaloRhoTimesAreaLeft(),
  fCaloRhoTimesAreaRight(),
  fCaloRhoTimesAreaBack(),
  fCaloTruePt(NULL),
  fCaloTruePtNotProper(NULL),
  fCaloTruePtNoIsPrimary(NULL),
  fCaloTrueRecPtNoIsPrimary(NULL),
  fCaloTrueWithoutConvPt(NULL),
  fCaloTruePtPrimary(NULL),
  fCaloTruePtDecay(NULL),
  fCaloTruePtDecayFoundOtherInCluster(NULL),
  fCaloTruePtDecayOtherInAcc(NULL),
  fCaloTruePtDecayOtherInAccAboveMinEnergy(NULL),
  fCaloTruePtTaggedCalo(NULL),
  fCaloTruePtTaggedAsDecayCalo(NULL),
  fCaloTrueRecPt(NULL),
  fCaloTrueRecPtvsTruePt(NULL),
  fCaloTrueRecPtvsTruePtNotProper(NULL),
  fCaloTrueWithoutConvRecPt(NULL),
  fCaloTrueRecPtPrimary(NULL),
  fCaloTrueRecPtDecay(NULL),
  fCaloTrueRecPtDecayFoundOtherInCluster(NULL),
  fCaloTrueRecPtDecayOtherInAcc(NULL),
  fCaloTrueRecPtDecayOtherInAccAboveMinEnergy(NULL),
  fCaloTrueRecPtTaggedCalo(NULL),
  fCaloTrueRecPtTaggedAsDecayCalo(NULL),
  fCaloTrueIsoCharged(),
  fCaloTrueIsoMCCharged(),
  fCaloTrueIsoNeutral(),
  fCaloTrueIsoFull(),
  fCaloTrueIsoCorr(),
  fCaloTrueIsoCell(),
  fCaloTrueIsoCharged_FromDecay(),
  fCaloTrueIsoNeutral_FromDecay(),
  fCaloTrueIsoFull_FromDecay(),
  fCaloTrueIsoCell_FromDecay(),
  fCaloTrueIsoCharged_FromDirect(),
  fCaloTrueIsoNeutral_FromDirect(),
  fCaloTrueIsoFull_FromDirect(),
  fCaloTrueIsoCell_FromDirect(),
  fCaloTrueSignalPtClusterCuts(NULL),
  fCaloTrueSignalPtClusterCutsOnlyPhoton(NULL),
  fCaloTrueSignalPtClusterCutsOnlyConv(NULL),
  fCaloTrueSignalPtClusterCutsM02Cuts(NULL),
  fCaloTrueSignalPtClusterCutsM02CutsIsoCuts(),
  fCaloTrueSignalPtClusterCutsM02CutsIsoCutsOnlyPhoton(),
  fCaloTrueSignalPtClusterCutsM02CutsIsoCutsOnlyConv(),
  fCaloTrueSignalRecPtClusterCuts(NULL),
  fCaloTrueSignalRecPtClusterCutsOnlyPhoton(NULL),
  fCaloTrueSignalRecPtClusterCutsOnlyConv(NULL),
  fCaloTrueSignalRecPtClusterCutsM02Cuts(NULL),
  fCaloTrueSignalRecPtClusterCutsM02CutsIsoCuts(),
  fCaloTrueSignalRecPtClusterCutsM02CutsIsoCutsOnlyPhoton(),
  fCaloTrueSignalRecPtClusterCutsM02CutsIsoCutsOnlyConv(),
  fCaloPtIsoCharged(),
  fCaloPtIsoNeutral(),
  fCaloPtIsoFull(),
  fCaloPtIsoCell(),
  fCaloPtTaggedCaloIsoCharged(),
  fCaloPtTaggedCaloIsoNeutral(),
  fCaloPtTaggedCaloIsoFull(),
  fCaloPtTaggedCaloIsoCell(),
  fCaloTruePtIsoCharged(),
  fCaloTruePtIsoNeutral(),
  fCaloTruePtIsoFull(),
  fCaloTruePtIsoCell(),
  fCaloTruePtTaggedCaloIsoCharged(),
  fCaloTruePtTaggedCaloIsoNeutral(),
  fCaloTruePtTaggedCaloIsoFull(),
  fCaloTruePtTaggedCaloIsoCell(),
  fCaloTruePtIsoChargedFromDirect(),
  fCaloTruePtIsoChargedAndMCIsoChargedFromDirect(),
  fCaloTruePtIsoNeutralFromDirect(),
  fCaloTruePtIsoFullFromDirect(),
  fCaloTruePtIsoCellFromDirect(),
  fCaloTruePtTaggedCaloIsoChargedFromDirect(),
  fCaloTruePtTaggedCaloIsoNeutralFromDirect(),
  fCaloTruePtTaggedCaloIsoFullFromDirect(),
  fCaloTruePtTaggedCaloIsoCellFromDirect(),
  fCaloTrueRecPtIsoCharged(),
  fCaloTrueRecPtIsoNeutral(),
  fCaloTrueRecPtIsoFull(),
  fCaloTrueRecPtIsoCell(),
  fCaloTrueRecPtIsoChargedFromDirect(),
  fCaloTrueRecPtIsoChargedAndMCIsoChargedFromDirect(),
  fCaloTrueRecPtIsoNeutralFromDirect(),
  fCaloTrueRecPtIsoFullFromDirect(),
  fCaloTrueRecPtIsoCellFromDirect(),
  fCaloTrueRecPtTaggedCaloIsoCharged(),
  fCaloTrueRecPtTaggedCaloIsoNeutral(),
  fCaloTrueRecPtTaggedCaloIsoFull(),
  fCaloTrueRecPtTaggedCaloIsoCell(),
  fCaloTruePtMCIsoCharged(),
  fCaloTruePtMCIsoNeutral(),
  fCaloTruePtMCIsoFull(),
  fCaloTruePtMCIsoCell(),
  fCaloTruePtTaggedCaloMCIsoCharged(),
  fCaloTruePtTaggedCaloMCIsoNeutral(),
  fCaloTruePtTaggedCaloMCIsoFull(),
  fCaloTruePtTaggedCaloMCIsoCell(),
  fCaloTrueRecPtMCIsoCharged(),
  fCaloTrueRecPtMCIsoNeutral(),
  fCaloTrueRecPtMCIsoFull(),
  fCaloTrueRecPtMCIsoCell(),
  fCaloTrueRecPtTaggedCaloMCIsoCharged(),
  fCaloTrueRecPtTaggedCaloMCIsoNeutral(),
  fCaloTrueRecPtTaggedCaloMCIsoFull(),
  fCaloTrueRecPtTaggedCaloMCIsoCell(),
  fCaloInvMass(NULL),
  fCaloInvMassChargedIsolated(),
  fCaloInvMassAntiChargedIsolated(),
  fCaloInvMassNeutralIsolated(),
  fCaloInvMassAntiNeutralIsolated(),
  fCaloInvMassCellIsolated(),
  fCaloInvMassAntiCellIsolated(),
  fCaloInvMassFullIsolated(),
  fCaloInvMassAntiFullIsolated(),
  fCaloTrueInvMass(NULL),
  fCaloTrueInvMass_FromDecay(NULL),
  fCaloTrueInvMass_FromDirect(NULL),
  fCaloM02(),
  fCaloM02ChargedIsolated(),
  fCaloM02ChargedIsolated_FromDirect(),
  fCaloM02ChargedIsolated_FromDirectConvOnly(),
  fCaloM02ChargedIsolated_FromDirectNormOnly(),
  fCaloM02AntiChargedIsolated(),
  fCaloM02NeutralIsolated(),
  fCaloM02AntiNeutralIsolated(),
  fCaloM02CellIsolated(),
  fCaloM02AntiCellIsolated(),
  fCaloM02FullIsolated(),
  fCaloM02AntiFullIsolated(),
  fCaloM02vsIsoVsPt(),
  fCaloTrueM02(),
  fCaloTrueM02_FromDecay(),
  fCaloTrueM02_FromDirect(),
  fHistoMCHeaders(NULL),
  fGenPhotonPt(NULL),
  fGenPhotonPt_FromDecay(NULL),
  fGenPhotonPt_FromDirect(NULL),
  fGenPhotonPtInEMCalAcc(NULL),
  fGenPhotonPtInEMCalAcc_FromDecay(NULL),
  fGenPhotonPtInEMCalAcc_FromDirect(NULL),
  fGenPhotonPtInEMCalAcc_FromDirect_OnlyPhoton(NULL),
  fGenPhotonPtInEMCalAcc_FromDirect_OnlyConv(NULL),
  fGenPhotonPtInEMCalAccChargedMCIso_FromDirect(),
  fGenPhotonPtInEMCalAccChargedMCIso_FromDirect_OnlyPhoton(),
  fGenPhotonPtInEMCalAccChargedMCIso_FromDirect_OnlyConv(),
  fGenPhotonChargedMCIsoInEMCalAcc_FromDirect(),
  fGenPhotonPtFoundNormCluster(NULL),
  fGenPhotonPtFoundTaggingCluster(NULL),
  fGenPhotonPtFoundIsoCluster(NULL),
  fGenPhotonEFoundNoClusterVsCellE(NULL),
  fGenPi0Pt(NULL),
  fGenPi0PtInEMCalAcc(NULL),
  fGenPi0PtInEMCalAcc_BothGammaInEMCal(NULL),
  fGenPi0PtInEMCalAcc_BothGammaInClusters(NULL),
  fRhoOutName("Rho"),
  fRhoOutNameMC("RhoMC"),
  fTreeBuffSize(60*1024*1024),
  fMemCountAOD(0),
  fTrackMatcherRunningMode(0),
  fAntiIsolationE(),
  fMinM02(0),
  fMaxM02(9999),
  fChargedRho(0),
  fChargedRhoMC(0),
  fChargedRhoTimesArea(),
  fExclusionRadius(0.4),
  fDebug(0),
  fOutlierJetReader(0),
  fIsFromDesiredHeader(kTRUE),
  fIsOverlappingWithOtherHeader(kFALSE),
  fAllowOverlapHeaders(kTRUE),
  fBuffer_EventRho(0),
  fBuffer_EventRhoMC(0),
  fBuffer_EventWeight(0),
  fBuffer_EventXsection(0),
  fBuffer_EventNtrials(0),
  fBuffer_EventIsTriggered(0),
  fBuffer_ClusterE(0), 
  fBuffer_ClusterPx(0), 
  fBuffer_ClusterPy(0), 
  fBuffer_ClusterPz(0), 
  fBuffer_ClusterM02(0), 
  fBuffer_ClusterM02Recalc(0), 
  fBuffer_ClusterM20(0), 
  fBuffer_ClusterV1SplitMass(0), 
  fBuffer_ClusterNLM(0), 
  fBuffer_ClusterSM(0), 
  fBuffer_ClusterEFrac(0), 
  fBuffer_ClusterIsoCharged1(0), 
  fBuffer_ClusterIsoCharged2(0), 
  fBuffer_ClusterIsoCharged3(0), 
  fBuffer_ClusterIsoBckLeft(0), 
  fBuffer_ClusterMatchTrackdEta(0), 
  fBuffer_ClusterMatchTrackdPhi(0), 
  fBuffer_ClusterMatchTrackP(0), 
  fBuffer_ClusterMatchTrackPt(0), 
  fBuffer_ClusterMatchTrackIsConv(0), 
  fBuffer_TrueClusterE(0), 
  fBuffer_TrueClusterPx(0), 
  fBuffer_TrueClusterPy(0), 
  fBuffer_TrueClusterPz(0), 
  fBuffer_TrueClusterLeadingEFrac(0),
  fBuffer_TrueClusterMCIsoCharged1(0), 
  fBuffer_TrueClusterMCIsoCharged2(0), 
  fBuffer_TrueClusterMCIsoCharged3(0), 
  fBuffer_TrueClusterMCIsoBckLeft(0), 
  fBuffer_TrueClusterMCTag(0), 
  fBuffer_TrueClusterIsConv(0),
  fBuffer_GenPhotonE(0),
  fBuffer_GenPhotonPx(0),
  fBuffer_GenPhotonPy(0),
  fBuffer_GenPhotonPz(0),
  fBuffer_GenPhotonMCIsoCharged1(0),
  fBuffer_GenPhotonMCIsoCharged2(0),
  fBuffer_GenPhotonMCIsoCharged3(0),
  fBuffer_GenPhotonMCIsoBckLeft(0),
  fBuffer_GenPhotonIsConv(0),
  fBuffer_GenPhotonMCTag(0)
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());

  SetEtaMatching(0.010,4.07,-2.5);
  SetPhiMatching(0.015,3.65,3.65);
  SetEOverP(1.75);

  SetPi0TaggingWindow(0.120,0.145);
  SetEtaTaggingWindow(0.5,0.6);

  SetAntiIsolationE(5.,10);
  fChargedRhoTimesArea[0] = 0.;
  fChargedRhoTimesArea[1] = 0.;
  fChargedRhoTimesArea[2] = 0.;
  fChargedRhoTimesArea[3] = 0.;
  fChargedRhoTimesArea[4] = 0.;
}

//________________________________________________________________________
AliAnalysisTaskGammaIsoTree::~AliAnalysisTaskGammaIsoTree()
{
  // default deconstructor
  // if(fConversionCandidates){
  //   delete fConversionCandidates;
  //   fConversionCandidates = 0x0;
  // }
  if(fClusterEMCalCandidates){
    delete fClusterEMCalCandidates;
    fClusterEMCalCandidates = 0x0;
  }
  if(fClusterEMCalCandidatesIsolation){
    delete fClusterEMCalCandidatesIsolation;
    fClusterEMCalCandidatesIsolation = 0x0;
  }
  if(fClusterEMCalCandidatesTagging){
    delete fClusterEMCalCandidatesTagging;
    fClusterEMCalCandidatesTagging = 0x0;
  }
  if(fClusterPHOSCandidates){
    delete fClusterPHOSCandidates;
    fClusterPHOSCandidates = 0x0;
  }
  // if(fTracks){
  //   delete fTracks;
  //   fTracks = 0x0;
  // }
  // if(fMCParticles){
  //   delete fMCParticles;
  //   fMCParticles = 0x0;
  // }

}
//________________________________________________________________________
void AliAnalysisTaskGammaIsoTree::UserCreateOutputObjects()
{
  fOutputList                         = new TList();
  fOutputList->SetOwner(kTRUE);

  if(((AliConvEventCuts*)fEventCuts)->GetCutHistograms()){
    fOutputList->Add(((AliConvEventCuts*)fEventCuts)->GetCutHistograms());
  }

  if(((AliCaloPhotonCuts*)fClusterCutsEMC)->GetCutHistograms()){
    fOutputList->Add(((AliCaloPhotonCuts*)fClusterCutsEMC)->GetCutHistograms());
  }

  if(((AliCaloPhotonCuts*)fClusterCutsEMCTrackMatching)->GetCutHistograms()){
    fOutputList->Add(((AliCaloPhotonCuts*)fClusterCutsEMCTrackMatching)->GetCutHistograms());
  }

  if(((AliCaloPhotonCuts*)fClusterCutsIsolationEMC)->GetCutHistograms()){
    fOutputList->Add(((AliCaloPhotonCuts*)fClusterCutsIsolationEMC)->GetCutHistograms());
  }

  if(((AliCaloPhotonCuts*)fClusterCutsPHOS)->GetCutHistograms()){
    fOutputList->Add(((AliCaloPhotonCuts*)fClusterCutsPHOS)->GetCutHistograms());
  }

  if(((AliConversionPhotonCuts*)fConvCuts)->GetCutHistograms()){
    fOutputList->Add(((AliConversionPhotonCuts*)fConvCuts)->GetCutHistograms());
  }

  if(!fDoOwnTrackMatching){
    TObjArray* tasks = AliAnalysisManager::GetAnalysisManager()->GetTasks();
    for(Int_t i =0; i < tasks->GetEntriesFast();i++){
      auto task = tasks->At(i);
      TString name = task->GetName();
      cout << "NAME OF TASK = " << name.Data() << endl;
      if(name.Contains("CaloTrackMatcher_Signal")){
          fOutputList->Add(((AliCaloTrackMatcher*) AliAnalysisManager::GetAnalysisManager()->GetTask(name.Data()))->GetCaloTrackMatcherHistograms());
      }
    }

  }

  // If requested, load Jet finder for outlier removal
    if(fEventCuts->GetUseJetFinderForOutliers()){
    fOutlierJetReader=(AliAnalysisTaskJetOutlierRemoval*)AliAnalysisManager::GetAnalysisManager()->GetTask("AliAnalysisTaskJetOutlierRemoval");
    if(!fOutlierJetReader){AliFatal("Error: No AliAnalysisTaskJetOutlierRemoval");} // GetV0Reader
    else{printf("Found AliAnalysisTaskJetOutlierRemoval used for outlier removal!\n");}
  }


  fGeneralFolder          = new TList();
  fGeneralFolder->SetName("general");
  fGeneralFolder->SetOwner(kTRUE);
  fOutputList->Add(fGeneralFolder);

  fHistoNEvents           = new TH1F("NEvents","NEvents",14,-0.5,13.5);
  fHistoNEvents->GetXaxis()->SetBinLabel(1,"Accepted");
  fHistoNEvents->GetXaxis()->SetBinLabel(2,"Centrality");
  fHistoNEvents->GetXaxis()->SetBinLabel(3,"Miss. MC or inc. ev.");
  if (((AliConvEventCuts*)fEventCuts)->IsSpecialTrigger() > 1 ){
    TString TriggerNames  = "Not Trigger: ";
    TriggerNames          = TriggerNames+ ( (AliConvEventCuts*)fEventCuts)->GetSpecialTriggerName();
    fHistoNEvents->GetXaxis()->SetBinLabel(4,TriggerNames.Data());
  } else {
    fHistoNEvents->GetXaxis()->SetBinLabel(4,"Trigger");
  }
  fHistoNEvents->GetXaxis()->SetBinLabel(5,"Vertex Z");
  fHistoNEvents->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
  fHistoNEvents->GetXaxis()->SetBinLabel(7,"Pile-Up");
  fHistoNEvents->GetXaxis()->SetBinLabel(8,"no SDD");
  fHistoNEvents->GetXaxis()->SetBinLabel(9,"no V0AND");
  fHistoNEvents->GetXaxis()->SetBinLabel(10,"EMCAL/TPC problem");
  fHistoNEvents->GetXaxis()->SetBinLabel(11,"rejectedForJetJetMC");
  fHistoNEvents->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
  fHistoNEvents->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
  fHistoNEvents->GetXaxis()->SetBinLabel(14,"Pileup V0M-TPCout Tracks");
  fHistoNEvents->GetYaxis()->SetTitle("N_{events}");
  fGeneralFolder->Add(fHistoNEvents);
  fHistoNEvents->Sumw2();

  fHistoChargedIso           = new TH1F("fHistoChargedIso","fHistoChargedIso",500,-0.5,50);
  fOutputList->Add(fHistoChargedIso);
  fHistoChargedIso->Sumw2();

  fHistoTaggingPCMPCM           = new TH2F("fHistoTaggingPCMPCM","fHistoTaggingPCMPCM;M (GeV/c^2); photon p_{T} (GeV/c)",500,0.,1.,100,0,50.);
  fGeneralFolder->Add(fHistoTaggingPCMPCM);
  fHistoTaggingPCMEMC           = new TH2F("fHistoTaggingPCMEMC","fHistoTaggingPCMEMC;M (GeV/c^2); photon p_{T} (GeV/c)",500,0.,1.,100,0,50.);
  fGeneralFolder->Add(fHistoTaggingPCMEMC);
  fHistoTaggingEMCPCM           = new TH2F("fHistoTaggingEMCPCM","fHistoTaggingEMCPCM;M (GeV/c^2); photon p_{T} (GeV/c)",500,0.,1.,100,0,50.);
  fGeneralFolder->Add(fHistoTaggingEMCPCM);
  fHistoTaggingEMCEMC           = new TH2F("fHistoTaggingEMCEMC","fHistoTaggingEMCEMC;M (GeV/c^2); photon p_{T} (GeV/c)",500,0.,1.,100,0,50.);
  fGeneralFolder->Add(fHistoTaggingEMCEMC);

  fHistoTaggingPCMPCM->Sumw2();
  fHistoTaggingPCMEMC->Sumw2();
  fHistoTaggingEMCPCM->Sumw2();
  fHistoTaggingEMCEMC->Sumw2();

  if(fIsMC > 1){
    fHistoNEventsWOWeight           = new TH1F("NEventsWOWeight","NEventsWOWeight",14,-0.5,13.5);
    fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(1,"Accepted");
    fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(2,"Centrality");
    fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(3,"Miss. MC or inc. ev.");
    if (((AliConvEventCuts*)fEventCuts)->IsSpecialTrigger() > 1 ){
      TString TriggerNames  = "Not Trigger: ";
      TriggerNames          = TriggerNames+ ( (AliConvEventCuts*)fEventCuts)->GetSpecialTriggerName();
      fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(4,TriggerNames.Data());
    } else {
      fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(4,"Trigger");
    }
    fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(5,"Vertex Z");
    fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
    fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(7,"Pile-Up");
    fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(8,"no SDD");
    fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(9,"no V0AND");
    fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(11,"rejectedForJetJetMC");
    fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(10,"EMCAL problem");
    fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
    fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
    fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(14,"Pileup V0M-TPCout Tracks");
    fHistoNEventsWOWeight->GetYaxis()->SetTitle("N_{events}");
    fHistoNEventsWOWeight->Sumw2();
    fGeneralFolder->Add(fHistoNEventsWOWeight);
  }

  //
  // ─── CONVERSION HISOGRAMS ───────────────────────────────────────────────────────
  //
  Double_t minPt = 0;
  Double_t maxPt = 150;
  Int_t  nPtBins = 300;

  Double_t minMass = 0.;
  Double_t maxMass = 2.;
  Int_t nMassBins = 200;

  fConvFolderRec          = new TList();
  fConvFolderRec->SetName("convPhotonsRec");
  fConvFolderRec->SetOwner(kTRUE);
  fOutputList->Add(fConvFolderRec);

  if(fUseHistograms){

    fConvPt = new TH1F("fConvPt","conversion photons in EMC acc;p_{T} (GeV/c); counts",nPtBins,minPt,maxPt);
    fConvPtBeforeAcc = new TH1F("fConvPtBeforeAcc", "conversion photons all acc;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);

    fConvPtTaggedCalo = new TH1F("fConvPtTaggedCalo", "conversion photons that survived tagging;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
    fConvPtTaggedAsDecayCalo = new TH1F("fConvPtTaggedAsDecayCalo", "conversion photons that survived tagging;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);



    fConvPt->Sumw2();
    fConvPtBeforeAcc->Sumw2();
    fConvPtTaggedCalo->Sumw2();
    fConvPtTaggedAsDecayCalo->Sumw2();

    fConvFolderRec->Add(fConvPt);
    fConvFolderRec->Add(fConvPtBeforeAcc);
    fConvFolderRec->Add(fConvPtTaggedCalo);
    fConvFolderRec->Add(fConvPtTaggedAsDecayCalo);

    for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
    {
      TH2F *convIsoCharged = new TH2F(Form("convIsoCharged%i",r), Form("charged track ISO in R < %1.1f;#sum p_{T} (GeV/c); conv. p_{T} (GeV/c)",fTrackIsolationR.at(r)), nPtBins,-10.,maxPt,nPtBins,minPt,maxPt);
      fConvIsoCharged[r] = (TH2F*) convIsoCharged->Clone(Form("fConvIsoCharged_R%1.1f",fTrackIsolationR.at(r)));
      fConvFolderRec->Add(fConvIsoCharged[r]);
      fConvIsoCharged[r]->Sumw2();
    }
    for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
    {
      TH2F *convIsoNeutral = new TH2F(Form("fConvIsoNeutral_%i",r), Form("Neutral track ISO in R < %1.1f;#sum p_{T} (GeV/c); conv. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,-10.,maxPt,nPtBins,minPt,maxPt);
      fConvIsoNeutral[r] = (TH2F*) convIsoNeutral->Clone(Form("fConvIsoNeutral_R%1.1f",fNeutralIsolationR.at(r)));
      fConvIsoNeutral[r]->Sumw2();
      fConvFolderRec->Add(fConvIsoNeutral[r]);

      TH2F *convIsoFull = new TH2F(Form("fConvIsoFull_%i",r), Form("Full track ISO in R < %1.1f;#sum p_{T} (GeV/c); conv. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,-10.,maxPt,nPtBins,minPt,maxPt);
      fConvIsoFull[r] = (TH2F*) convIsoFull->Clone(Form("fConvIsoFull_R%1.1f",fNeutralIsolationR.at(r)));
      fConvIsoFull[r]->Sumw2();
      fConvFolderRec->Add(fConvIsoFull[r]);

      TH2F *convIsoCell = new TH2F(Form("fConvIsoCell_%i",r), Form("Full track ISO in R < %1.1f;#sum p_{T} (GeV/c); conv. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,-10.,maxPt,nPtBins,minPt,maxPt);
      fConvIsoCell[r] = (TH2F*) convIsoCell->Clone(Form("fConvIsoCell_R%1.1f",fNeutralIsolationR.at(r)));
      fConvIsoCell[r]->Sumw2();
      fConvFolderRec->Add(fConvIsoCell[r]);
    }

    fConvFolderTrue          = new TList();
    fConvFolderTrue->SetName("convPhotonsTrue");
    fConvFolderTrue->SetOwner(kTRUE);
    if(fIsMC > 0){
      fOutputList->Add(fConvFolderTrue);
      fConvTruePt = new TH1F("fConvTruePt", "validated conversion photons in EMC acceptance;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fConvTruePtPrimary = new TH1F("fConvTruePtPrimary", "conversion photon that has not a pi0 etc. as mother;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fConvTruePtDecay = new TH1F("fConvTruePtDecay", "conversion photon from decay;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fConvTruePtDecayFoundOtherInCluster = new TH1F("fConvTruePtDecayFoundOtherInCluster", "conversion photon from decay, where the other decay particle was found in EMC;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fConvTruePtDecayOtherInAcc = new TH1F("fConvTruePtDecayOtherInAcc", "conversion photon from decay, where the other decay particle was found;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fConvTruePtDecayOtherInAccAboveMinEnergy = new TH1F("fConvTruePtDecayOtherInAccAboveMinEnergy", "conversion photon from decay, where the other decay particle is in acc. and above 0.7 GeV;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);

      // iso and tagging studies
      fConvTruePtTaggedCalo = new TH1F("fConvTruePtTaggedCalo", "conversion photons that survived tagging;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fConvTruePtTaggedAsDecayCalo = new TH1F("fConvTruePtTaggedAsDecayCalo", "conversion photons that survived tagging;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);

      // True with rec pT
      fConvTrueRecPt = new TH1F("fConvTrueRecPt", "validated conversion photons in EMC acceptance;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fConvTrueRecPtPrimary = new TH1F("fConvTrueRecPtPrimary", "conversion photon that has not a pi0 etc. as mother;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fConvTrueRecPtDecay = new TH1F("fConvTrueRecPtDecay", "conversion photon from decay;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fConvTrueRecPtDecayFoundOtherInCluster = new TH1F("fConvTrueRecPtDecayFoundOtherInCluster", "conversion photon from decay, where the other decay particle was found in EMC;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fConvTrueRecPtDecayOtherInAcc = new TH1F("fConvTrueRecPtDecayOtherInAcc", "conversion photon from decay, where the other decay particle was found;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fConvTrueRecPtDecayOtherInAccAboveMinEnergy = new TH1F("fConvTrueRecPtDecayOtherInAccAboveMinEnergy", "conversion photon from decay, where the other decay particle is in acc. and above 0.7 GeV;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);

      // iso and tagging studies
      fConvTrueRecPtTaggedCalo = new TH1F("fConvTrueRecPtTaggedCalo", "conversion photons that survived tagging;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fConvTrueRecPtTaggedAsDecayCalo = new TH1F("fConvTrueRecPtTaggedAsDecayCalo", "conversion photons that survived tagging;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);

      fConvTruePt->Sumw2();
      fConvTruePtPrimary->Sumw2();
      fConvTruePtDecay->Sumw2();
      fConvTruePtDecayFoundOtherInCluster->Sumw2();
      fConvTruePtDecayOtherInAcc->Sumw2();
      fConvTruePtDecayOtherInAccAboveMinEnergy->Sumw2();
      fConvTruePtTaggedCalo->Sumw2();
      fConvTruePtTaggedAsDecayCalo->Sumw2();
      fConvTrueRecPt->Sumw2();
      fConvTrueRecPtPrimary->Sumw2();
      fConvTrueRecPtDecay->Sumw2();
      fConvTrueRecPtDecayFoundOtherInCluster->Sumw2();
      fConvTrueRecPtDecayOtherInAcc->Sumw2();
      fConvTrueRecPtDecayOtherInAccAboveMinEnergy->Sumw2();
      fConvTrueRecPtTaggedCalo->Sumw2();
      fConvTrueRecPtTaggedAsDecayCalo->Sumw2();

      // add to folders
      fConvFolderTrue->Add(fConvTruePt);
      fConvFolderTrue->Add(fConvTruePtPrimary);
      fConvFolderTrue->Add(fConvTruePtDecay);
      fConvFolderTrue->Add(fConvTruePtDecayFoundOtherInCluster);
      fConvFolderTrue->Add(fConvTruePtDecayOtherInAcc);
      fConvFolderTrue->Add(fConvTruePtDecayOtherInAccAboveMinEnergy);

      // iso and tagging studies
      fConvFolderTrue->Add(fConvTruePtTaggedCalo);
      fConvFolderTrue->Add(fConvTruePtTaggedAsDecayCalo);

      // True with rec pT
      fConvFolderTrue->Add(fConvTrueRecPt);
      fConvFolderTrue->Add(fConvTrueRecPtPrimary);
      fConvFolderTrue->Add(fConvTrueRecPtDecay);
      fConvFolderTrue->Add(fConvTrueRecPtDecayFoundOtherInCluster);
      fConvFolderTrue->Add(fConvTrueRecPtDecayOtherInAcc);
      fConvFolderTrue->Add(fConvTrueRecPtDecayOtherInAccAboveMinEnergy);

      // iso and tagging studies
      fConvFolderTrue->Add(fConvTrueRecPtTaggedCalo);
      fConvFolderTrue->Add(fConvTrueRecPtTaggedAsDecayCalo);



      for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
      {
        TH2F *convTrueIsoCharged = new TH2F(Form("fConvTrueIsoCharged_%i",r), Form("Charged track ISO in R < %1.1f;#sum p_{T} (GeV/c); counts",fTrackIsolationR.at(r)), nPtBins,minPt,maxPt,nPtBins,minPt,maxPt);
        fConvTrueIsoCharged[r] = (TH2F*) convTrueIsoCharged->Clone(Form("fConvTrueIsoCharged_R%1.1f",fTrackIsolationR.at(r)));
        fConvFolderTrue->Add(fConvTrueIsoCharged[r]);
        fConvTrueIsoCharged[r]->Sumw2();

        TH2F *convTrueIsoCharged_FromDecay = new TH2F(Form("fConvTrueIsoCharged_FromDecay_%i",r), Form("Charged track ISO in R < %1.1f;#sum p_{T} (GeV/c); counts",fTrackIsolationR.at(r)), nPtBins,minPt,maxPt,nPtBins,minPt,maxPt);
        fConvTrueIsoCharged_FromDecay[r] = (TH2F*) convTrueIsoCharged_FromDecay->Clone(Form("fConvTrueIsoCharged_FromDecay_R%1.1f",fTrackIsolationR.at(r)));
        fConvFolderTrue->Add(fConvTrueIsoCharged_FromDecay[r]);
        fConvTrueIsoCharged_FromDecay[r]->Sumw2();

        TH2F *convTrueIsoCharged_FromDirect = new TH2F(Form("fConvTrueIsoCharged_FromDirect_%i",r), Form("Charged track ISO in R < %1.1f;#sum p_{T} (GeV/c); counts",fTrackIsolationR.at(r)), nPtBins,minPt,maxPt,nPtBins,minPt,maxPt);
        fConvTrueIsoCharged_FromDirect[r] = (TH2F*) convTrueIsoCharged_FromDirect->Clone(Form("fConvTrueIsoCharged_FromDirect_R%1.1f",fTrackIsolationR.at(r)));
        fConvTrueIsoCharged_FromDirect[r]->Sumw2();
        fConvFolderTrue->Add(fConvTrueIsoCharged_FromDirect[r]);
      }

      for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
      {
        TH2F *convTrueIsoNeutral = new TH2F(Form("fConvTrueIsoNeutral_%i",r), Form("Neutral track ISO in R < %1.1f;#sum p_{T} (GeV/c); conv. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,nPtBins,minPt,maxPt);
        fConvTrueIsoNeutral[r] = (TH2F*) convTrueIsoNeutral->Clone(Form("fConvTrueIsoNeutral_R%1.1f",fNeutralIsolationR.at(r)));
        fConvTrueIsoNeutral[r]->Sumw2();
        fConvFolderTrue->Add(fConvTrueIsoNeutral[r]);

        TH2F *convTrueIsoFull = new TH2F(Form("fConvTrueIsoFull_%i",r), Form("Full track ISO in R < %1.1f;#sum p_{T} (GeV/c); conv. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,nPtBins,minPt,maxPt);
        fConvTrueIsoFull[r] = (TH2F*) convTrueIsoFull->Clone(Form("fConvTrueIsoFull_R%1.1f",fNeutralIsolationR.at(r)));
        fConvTrueIsoFull[r]->Sumw2();
        fConvFolderTrue->Add(fConvTrueIsoFull[r]);

        TH2F *convTrueIsoCell = new TH2F(Form("fConvTrueIsoCell_%i",r), Form("Full track ISO in R < %1.1f;#sum p_{T} (GeV/c); conv. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,nPtBins,minPt,maxPt);
        fConvTrueIsoCell[r] = (TH2F*) convTrueIsoCell->Clone(Form("fConvTrueIsoCell_R%1.1f",fNeutralIsolationR.at(r)));
        fConvTrueIsoCell[r]->Sumw2();
        fConvFolderTrue->Add(fConvTrueIsoCell[r]);

        TH2F *convTrueIsoNeutral_FromDecay = new TH2F(Form("fConvTrueIsoNeutral_FromDecay_%i",r), Form("Neutral track ISO in R < %1.1f;#sum p_{T} (GeV/c); conv. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,nPtBins,minPt,maxPt);
        fConvTrueIsoNeutral_FromDecay[r] = (TH2F*) convTrueIsoNeutral_FromDecay->Clone(Form("fConvTrueIsoNeutral_FromDecay_R%1.1f",fNeutralIsolationR.at(r)));
        fConvTrueIsoNeutral_FromDecay[r]->Sumw2();
        fConvFolderTrue->Add(fConvTrueIsoNeutral_FromDecay[r]);

        TH2F *convTrueIsoFull_FromDecay = new TH2F(Form("fConvTrueIsoFull_FromDecay_%i",r), Form("Full track ISO in R < %1.1f;#sum p_{T} (GeV/c); conv. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,nPtBins,minPt,maxPt);
        fConvTrueIsoFull_FromDecay[r] = (TH2F*) convTrueIsoFull_FromDecay->Clone(Form("fConvTrueIsoFull_FromDecay_R%1.1f",fNeutralIsolationR.at(r)));
        fConvTrueIsoFull_FromDecay[r]->Sumw2();
        fConvFolderTrue->Add(fConvTrueIsoFull_FromDecay[r]);

        TH2F *convTrueIsoCell_FromDecay = new TH2F(Form("fConvTrueIsoCell_FromDecay_%i",r), Form("Full track ISO in R < %1.1f;#sum p_{T} (GeV/c); conv. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,nPtBins,minPt,maxPt);
        fConvTrueIsoCell_FromDecay[r] = (TH2F*) convTrueIsoCell_FromDecay->Clone(Form("fConvTrueIsoCell_FromDecay_R%1.1f",fNeutralIsolationR.at(r)));
        fConvTrueIsoCell_FromDecay[r]->Sumw2();
        fConvFolderTrue->Add(fConvTrueIsoCell_FromDecay[r]);

        TH2F *convTrueIsoNeutral_FromDirect = new TH2F(Form("fConvTrueIsoNeutral_FromDirect_%i",r), Form("Neutral track ISO in R < %1.1f;#sum p_{T} (GeV/c); conv. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,nPtBins,minPt,maxPt);
        fConvTrueIsoNeutral_FromDirect[r] = (TH2F*) convTrueIsoNeutral_FromDirect->Clone(Form("fConvTrueIsoNeutral_FromDirect_R%1.1f",fNeutralIsolationR.at(r)));
        fConvTrueIsoNeutral_FromDirect[r]->Sumw2();
        fConvFolderTrue->Add(fConvTrueIsoNeutral_FromDirect[r]);

        TH2F *convTrueIsoFull_FromDirect = new TH2F(Form("fConvTrueIsoFull_FromDirect_%i",r), Form("Full track ISO in R < %1.1f;#sum p_{T} (GeV/c); conv. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,nPtBins,minPt,maxPt);
        fConvTrueIsoFull_FromDirect[r] = (TH2F*) convTrueIsoFull_FromDirect->Clone(Form("fConvTrueIsoFull_FromDirect_R%1.1f",fNeutralIsolationR.at(r)));
        fConvTrueIsoFull_FromDirect[r]->Sumw2();
        fConvFolderTrue->Add(fConvTrueIsoFull_FromDirect[r]);

        TH2F *convTrueIsoCell_FromDirect = new TH2F(Form("fConvTrueIsoCell_FromDirect_%i",r), Form("Full track ISO in R < %1.1f;#sum p_{T} (GeV/c); conv. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,nPtBins,minPt,maxPt);
        fConvTrueIsoCell_FromDirect[r] = (TH2F*) convTrueIsoCell_FromDirect->Clone(Form("fConvTrueIsoCell_FromDirect_R%1.1f",fNeutralIsolationR.at(r)));
        fConvTrueIsoCell_FromDirect[r]->Sumw2();
        fConvFolderTrue->Add(fConvTrueIsoCell_FromDirect[r]);
      }
    }
    for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
    {
      for (UInt_t e = 0; e < fTrackIsolationE.size(); e++)
      {
        TH1F *convPtIsoCharged = new TH1F(Form("fConvPtIsoCharged_%i_%i",r,e), Form("conversion photons with charged track ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
        fConvPtIsoCharged[r][e] = (TH1F*) convPtIsoCharged->Clone(Form("fConvPtIsoCharged_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
        fConvPtIsoCharged[r][e]->Sumw2();
        fConvFolderRec->Add(fConvPtIsoCharged[r][e]);
        TH1F *convPtTaggedCaloIsoCharged = new TH1F(Form("fConvPtTaggedCaloIsoCharged_%i_%i",r,e), Form("conversion photons with charged track ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
        fConvPtTaggedCaloIsoCharged[r][e] = (TH1F*) convPtTaggedCaloIsoCharged->Clone(Form("fConvPtTaggedCaloIsoCharged_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
        fConvPtTaggedCaloIsoCharged[r][e]->Sumw2();
        fConvFolderRec->Add(fConvPtTaggedCaloIsoCharged[r][e]);

        if(fIsMC>0){
          // true pT
          TH1F *convTruePtIsoCharged = new TH1F(Form("fConvTruePtIsoCharged_%i_%i",r,e), Form("conversion photons with charged track ISO < %1.1f GeV in R < %1.1f;gen. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtIsoCharged[r][e] = (TH1F*) convTruePtIsoCharged->Clone(Form("fConvTruePtIsoCharged_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fConvTruePtIsoCharged[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTruePtIsoCharged[r][e]);
          TH1F *convTruePtTaggedCaloIsoCharged = new TH1F(Form("fConvTruePtTaggedCaloIsoCharged_%i_%i",r,e), Form("conversion photons with charged track ISO < %1.1f GeV in R < %1.1f + tagging;gen. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtTaggedCaloIsoCharged[r][e] = (TH1F*) convTruePtTaggedCaloIsoCharged->Clone(Form("fConvTruePtTaggedCaloIsoCharged_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fConvTruePtTaggedCaloIsoCharged[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTruePtTaggedCaloIsoCharged[r][e]);

          TH1F *convTruePtIsoChargedFromDirect = new TH1F(Form("fConvTruePtIsoChargedFromDirect_%i_%i",r,e), Form("conversion photons with chargedFromDirect track ISO < %1.1f GeV in R < %1.1f;gen. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtIsoChargedFromDirect[r][e] = (TH1F*) convTruePtIsoChargedFromDirect->Clone(Form("fConvTruePtIsoChargedFromDirect_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fConvTruePtIsoChargedFromDirect[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTruePtIsoChargedFromDirect[r][e]);
          TH1F *convTruePtTaggedCaloIsoChargedFromDirect = new TH1F(Form("fConvTruePtTaggedCaloIsoChargedFromDirect_%i_%i",r,e), Form("conversion photons with chargedFromDirect track ISO < %1.1f GeV in R < %1.1f + tagging;gen. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtTaggedCaloIsoChargedFromDirect[r][e] = (TH1F*) convTruePtTaggedCaloIsoChargedFromDirect->Clone(Form("fConvTruePtTaggedCaloIsoChargedFromDirect_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fConvTruePtTaggedCaloIsoChargedFromDirect[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTruePtTaggedCaloIsoChargedFromDirect[r][e]);

          //mc iso
          TH1F *convTruePtMCIsoCharged = new TH1F(Form("fConvTruePtMCIsoCharged_%i_%i",r,e), Form("conversion photons with charged track MCIso < %1.1f GeV in R < %1.1f;gen. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtMCIsoCharged[r][e] = (TH1F*) convTruePtMCIsoCharged->Clone(Form("fConvTruePtMCIsoCharged_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fConvTruePtMCIsoCharged[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTruePtMCIsoCharged[r][e]);
          TH1F *convTruePtTaggedCaloMCIsoCharged = new TH1F(Form("fConvTruePtTaggedCaloMCIsoCharged_%i_%i",r,e), Form("conversion photons with charged track MCIso < %1.1f GeV in R < %1.1f + tagging;gen. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtTaggedCaloMCIsoCharged[r][e] = (TH1F*) convTruePtTaggedCaloMCIsoCharged->Clone(Form("fConvTruePtTaggedCaloMCIsoCharged_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fConvTruePtTaggedCaloMCIsoCharged[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTruePtTaggedCaloMCIsoCharged[r][e]);

          // rec Pt
          TH1F *convTrueRecPtIsoCharged = new TH1F(Form("fConvTrueRecPtIsoCharged_%i_%i",r,e), Form("conversion photons with charged track ISO < %1.1f GeV in R < %1.1f;rec. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtIsoCharged[r][e] = (TH1F*) convTrueRecPtIsoCharged->Clone(Form("fConvTrueRecPtIsoCharged_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fConvTrueRecPtIsoCharged[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTrueRecPtIsoCharged[r][e]);
          TH1F *convTrueRecPtIsoChargedFromDirect = new TH1F(Form("fConvTrueRecPtIsoChargedFromDirect_%i_%i",r,e), Form("conversion photons with charged track ISO < %1.1f GeV in R < %1.1f;rec. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtIsoChargedFromDirect[r][e] = (TH1F*) convTrueRecPtIsoChargedFromDirect->Clone(Form("fConvTrueRecPtIsoChargedFromDirect_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fConvTrueRecPtIsoChargedFromDirect[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTrueRecPtIsoChargedFromDirect[r][e]);
          TH1F *convTrueRecPtTaggedCaloIsoCharged = new TH1F(Form("fConvTrueRecPtTaggedCaloIsoCharged_%i_%i",r,e), Form("conversion photons with charged track ISO < %1.1f GeV in R < %1.1f + tagging;rec. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtTaggedCaloIsoCharged[r][e] = (TH1F*) convTrueRecPtTaggedCaloIsoCharged->Clone(Form("fConvTrueRecPtTaggedCaloIsoCharged_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fConvTrueRecPtTaggedCaloIsoCharged[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTrueRecPtTaggedCaloIsoCharged[r][e]);

          TH1F *convTrueRecPtMCIsoCharged = new TH1F(Form("fConvTrueRecPtMCIsoCharged_%i_%i",r,e), Form("conversion photons with charged track MCIso < %1.1f GeV in R < %1.1f;rec. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtMCIsoCharged[r][e] = (TH1F*) convTrueRecPtMCIsoCharged->Clone(Form("fConvTrueRecPtMCIsoCharged_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fConvTrueRecPtMCIsoCharged[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTrueRecPtMCIsoCharged[r][e]);
          TH1F *convTrueRecPtTaggedCaloMCIsoCharged = new TH1F(Form("fConvTrueRecPtTaggedCaloMCIsoCharged_%i_%i",r,e), Form("conversion photons with charged track MCIso < %1.1f GeV in R < %1.1f + tagging;rec. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtTaggedCaloMCIsoCharged[r][e] = (TH1F*) convTrueRecPtTaggedCaloMCIsoCharged->Clone(Form("fConvTrueRecPtTaggedCaloMCIsoCharged_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fConvTrueRecPtTaggedCaloMCIsoCharged[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTrueRecPtTaggedCaloMCIsoCharged[r][e]);
        }

      }

    }
    for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
    {
      for (UInt_t e = 0; e < fNeutralIsolationE.size(); e++)
      {
        TH1F *convPtIsoNeutral = new TH1F(Form("fConvPtIsoNeutral_%i_%i",r,e), Form("conversion photons with Neutral ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
        fConvPtIsoNeutral[r][e] = (TH1F*) convPtIsoNeutral->Clone(Form("fConvPtIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
        fConvPtIsoNeutral[r][e]->Sumw2();
        fConvFolderRec->Add(fConvPtIsoNeutral[r][e]);
        TH1F *convPtIsoFull = new TH1F(Form("fConvPtIsoFull_%i_%i",r,e), Form("conversion photons with Full ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
        fConvPtIsoFull[r][e] = (TH1F*) convPtIsoFull->Clone(Form("fConvPtIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
        fConvPtIsoFull[r][e]->Sumw2();
        fConvFolderRec->Add(fConvPtIsoFull[r][e]);
        TH1F *convPtIsoCell = new TH1F(Form("fConvPtIsoCell_%i_%i",r,e), Form("conversion photons with Cell ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
        fConvPtIsoCell[r][e] = (TH1F*) convPtIsoCell->Clone(Form("fConvPtIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
        fConvPtIsoCell[r][e]->Sumw2();
        fConvFolderRec->Add(fConvPtIsoCell[r][e]);

        TH1F *convPtTaggedCaloIsoNeutral = new TH1F(Form("fConvPtTaggedCaloIsoNeutral_%i_%i",r,e), Form("conversion photons with Neutral ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
        fConvPtTaggedCaloIsoNeutral[r][e] = (TH1F*) convPtTaggedCaloIsoNeutral->Clone(Form("fConvPtTaggedCaloIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
        fConvPtTaggedCaloIsoNeutral[r][e]->Sumw2();
        fConvFolderRec->Add(fConvPtTaggedCaloIsoNeutral[r][e]);
        TH1F *convPtTaggedCaloIsoFull = new TH1F(Form("fConvTaggedCaloIsoFull_%i_%i",r,e), Form("conversion photons with Full ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
        fConvPtTaggedCaloIsoFull[r][e] = (TH1F*) convPtTaggedCaloIsoFull->Clone(Form("fConvTaggedCaloIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
        fConvPtTaggedCaloIsoFull[r][e]->Sumw2();
        fConvFolderRec->Add(fConvPtTaggedCaloIsoFull[r][e]);
        TH1F *convPtTaggedCaloIsoCell = new TH1F(Form("fConvTaggedCaloIsoCell_%i_%i",r,e), Form("conversion photons with Cell ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
        fConvPtTaggedCaloIsoCell[r][e] = (TH1F*) convPtTaggedCaloIsoCell->Clone(Form("fConvTaggedCaloIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
        fConvPtTaggedCaloIsoCell[r][e]->Sumw2();
        fConvFolderRec->Add(fConvPtTaggedCaloIsoCell[r][e]);

        // True Pt
        if(fIsMC>0){
          TH1F *convTruePtIsoNeutral = new TH1F(Form("fConvTruePtIsoNeutral_%i_%i",r,e), Form("conversion photons with Neutral ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtIsoNeutral[r][e] = (TH1F*) convTruePtIsoNeutral->Clone(Form("fConvTruePtIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTruePtIsoNeutral[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTruePtIsoNeutral[r][e]);
          TH1F *convTruePtIsoFull = new TH1F(Form("fConvTruePtIsoFull_%i_%i",r,e), Form("conversion photons with Full ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtIsoFull[r][e] = (TH1F*) convTruePtIsoFull->Clone(Form("fConvTruePtIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTruePtIsoFull[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTruePtIsoFull[r][e]);
          TH1F *convTruePtIsoCell = new TH1F(Form("fConvTruePtIsoCell_%i_%i",r,e), Form("conversion photons with Cell ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtIsoCell[r][e] = (TH1F*) convTruePtIsoCell->Clone(Form("fConvTruePtIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTruePtIsoCell[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTruePtIsoCell[r][e]);

          TH1F *convTruePtTaggedCaloIsoNeutral = new TH1F(Form("fConvTruePtTaggedCaloIsoNeutral_%i_%i",r,e), Form("conversion photons with Neutral ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtTaggedCaloIsoNeutral[r][e] = (TH1F*) convTruePtTaggedCaloIsoNeutral->Clone(Form("fConvTruePtTaggedCaloIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTruePtTaggedCaloIsoNeutral[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTruePtTaggedCaloIsoNeutral[r][e]);
          TH1F *convTruePtTaggedCaloIsoFull = new TH1F(Form("fConvTruePtTaggedCaloIsoFull_%i_%i",r,e), Form("conversion photons with Full ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtTaggedCaloIsoFull[r][e] = (TH1F*) convTruePtTaggedCaloIsoFull->Clone(Form("fConvTruePtTaggedCaloIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTruePtTaggedCaloIsoFull[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTruePtTaggedCaloIsoFull[r][e]);
          TH1F *convTruePtTaggedCaloIsoCell = new TH1F(Form("fConvTruePtTaggedCaloIsoCell_%i_%i",r,e), Form("conversion photons with Cell ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtTaggedCaloIsoCell[r][e] = (TH1F*) convTruePtTaggedCaloIsoCell->Clone(Form("fConvTruePtTaggedCaloIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTruePtTaggedCaloIsoCell[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTruePtTaggedCaloIsoCell[r][e]);

          // from direct
          TH1F *convTruePtIsoNeutralFromDirect = new TH1F(Form("fConvTruePtIsoNeutralFromDirect_%i_%i",r,e), Form("conversion photons with Neutral ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtIsoNeutralFromDirect[r][e] = (TH1F*) convTruePtIsoNeutralFromDirect->Clone(Form("fConvTruePtIsoNeutralFromDirect_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTruePtIsoNeutralFromDirect[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTruePtIsoNeutralFromDirect[r][e]);
          TH1F *convTruePtIsoFullFromDirect = new TH1F(Form("fConvTruePtIsoFullFromDirect_%i_%i",r,e), Form("conversion photons with Full ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtIsoFullFromDirect[r][e] = (TH1F*) convTruePtIsoFullFromDirect->Clone(Form("fConvTruePtIsoFullFromDirect_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTruePtIsoFullFromDirect[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTruePtIsoFullFromDirect[r][e]);
          TH1F *convTruePtIsoCellFromDirect = new TH1F(Form("fConvTruePtIsoCellFromDirect_%i_%i",r,e), Form("conversion photons with Cell ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtIsoCellFromDirect[r][e] = (TH1F*) convTruePtIsoCellFromDirect->Clone(Form("fConvTruePtIsoCellFromDirect_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTruePtIsoCellFromDirect[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTruePtIsoCellFromDirect[r][e]);

          TH1F *convTruePtTaggedCaloIsoNeutralFromDirect = new TH1F(Form("fConvTruePtTaggedCaloIsoNeutralFromDirect_%i_%i",r,e), Form("conversion photons with Neutral ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtTaggedCaloIsoNeutralFromDirect[r][e] = (TH1F*) convTruePtTaggedCaloIsoNeutralFromDirect->Clone(Form("fConvTruePtTaggedCaloIsoNeutralFromDirect_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTruePtTaggedCaloIsoNeutralFromDirect[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTruePtTaggedCaloIsoNeutralFromDirect[r][e]);
          TH1F *convTruePtTaggedCaloIsoFullFromDirect = new TH1F(Form("fConvTruePtTaggedCaloIsoFullFromDirect_%i_%i",r,e), Form("conversion photons with Full ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtTaggedCaloIsoFullFromDirect[r][e] = (TH1F*) convTruePtTaggedCaloIsoFullFromDirect->Clone(Form("fConvTruePtTaggedCaloIsoFullFromDirect_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTruePtTaggedCaloIsoFullFromDirect[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTruePtTaggedCaloIsoFullFromDirect[r][e]);
          TH1F *convTruePtTaggedCaloIsoCellFromDirect = new TH1F(Form("fConvTruePtTaggedCaloIsoCellFromDirect_%i_%i",r,e), Form("conversion photons with Cell ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtTaggedCaloIsoCellFromDirect[r][e] = (TH1F*) convTruePtTaggedCaloIsoCellFromDirect->Clone(Form("fConvTruePtTaggedCaloIsoCellFromDirect_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTruePtTaggedCaloIsoCellFromDirect[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTruePtTaggedCaloIsoCellFromDirect[r][e]);

          //mc iso

          TH1F *convTruePtMCIsoNeutral = new TH1F(Form("fConvTruePtMCIsoNeutral_%i_%i",r,e), Form("conversion photons with Neutral MCIso < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtMCIsoNeutral[r][e] = (TH1F*) convTruePtMCIsoNeutral->Clone(Form("fConvTruePtMCIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTruePtMCIsoNeutral[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTruePtMCIsoNeutral[r][e]);
          TH1F *convTruePtMCIsoFull = new TH1F(Form("fConvTruePtMCIsoFull_%i_%i",r,e), Form("conversion photons with Full MCIso < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtMCIsoFull[r][e] = (TH1F*) convTruePtMCIsoFull->Clone(Form("fConvTruePtMCIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTruePtMCIsoFull[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTruePtMCIsoFull[r][e]);
          TH1F *convTruePtMCIsoCell = new TH1F(Form("fConvTruePtMCIsoCell_%i_%i",r,e), Form("conversion photons with Cell MCIso < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtMCIsoCell[r][e] = (TH1F*) convTruePtMCIsoCell->Clone(Form("fConvTruePtMCIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTruePtMCIsoCell[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTruePtMCIsoCell[r][e]);

          TH1F *convTruePtTaggedCaloMCIsoNeutral = new TH1F(Form("fConvTruePtTaggedCaloMCIsoNeutral_%i_%i",r,e), Form("conversion photons with Neutral MCIso < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtTaggedCaloMCIsoNeutral[r][e] = (TH1F*) convTruePtTaggedCaloMCIsoNeutral->Clone(Form("fConvTruePtTaggedCaloMCIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTruePtTaggedCaloMCIsoNeutral[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTruePtTaggedCaloMCIsoNeutral[r][e]);
          TH1F *convTruePtTaggedCaloMCIsoFull = new TH1F(Form("fConvTruePtTaggedCaloMCIsoFull_%i_%i",r,e), Form("conversion photons with Full MCIso < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtTaggedCaloMCIsoFull[r][e] = (TH1F*) convTruePtTaggedCaloMCIsoFull->Clone(Form("fConvTruePtTaggedCaloMCIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTruePtTaggedCaloMCIsoFull[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTruePtTaggedCaloMCIsoFull[r][e]);
          TH1F *convTruePtTaggedCaloMCIsoCell = new TH1F(Form("fConvTruePtTaggedCaloMCIsoCell_%i_%i",r,e), Form("conversion photons with Cell MCIso < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtTaggedCaloMCIsoCell[r][e] = (TH1F*) convTruePtTaggedCaloMCIsoCell->Clone(Form("fConvTruePtTaggedCaloMCIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTruePtTaggedCaloMCIsoCell[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTruePtTaggedCaloMCIsoCell[r][e]);

          // rec Pt

          TH1F *convTrueRecPtIsoNeutral = new TH1F(Form("fConvTrueRecPtIsoNeutral_%i_%i",r,e), Form("conversion photons with Neutral ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtIsoNeutral[r][e] = (TH1F*) convTrueRecPtIsoNeutral->Clone(Form("fConvTrueRecPtIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTrueRecPtIsoNeutral[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTrueRecPtIsoNeutral[r][e]);
          TH1F *convTrueRecPtIsoFull = new TH1F(Form("fConvTrueRecPtIsoFull_%i_%i",r,e), Form("conversion photons with Full ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtIsoFull[r][e] = (TH1F*) convTrueRecPtIsoFull->Clone(Form("fConvTrueRecPtIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTrueRecPtIsoFull[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTrueRecPtIsoFull[r][e]);
          TH1F *convTrueRecPtIsoCell = new TH1F(Form("fConvTrueRecPtIsoCell_%i_%i",r,e), Form("conversion photons with Cell ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtIsoCell[r][e] = (TH1F*) convTrueRecPtIsoCell->Clone(Form("fConvTrueRecPtIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTrueRecPtIsoCell[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTrueRecPtIsoCell[r][e]);

          TH1F *convTrueRecPtIsoNeutralFromDirect = new TH1F(Form("fConvTrueRecPtIsoNeutralFromDirect_%i_%i",r,e), Form("conversion photons with Neutral ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtIsoNeutralFromDirect[r][e] = (TH1F*) convTrueRecPtIsoNeutralFromDirect->Clone(Form("fConvTrueRecPtIsoNeutralFromDirect_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTrueRecPtIsoNeutralFromDirect[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTrueRecPtIsoNeutralFromDirect[r][e]);
          TH1F *convTrueRecPtIsoFullFromDirect = new TH1F(Form("fConvTrueRecPtIsoFullFromDirect_%i_%i",r,e), Form("conversion photons with Full ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtIsoFullFromDirect[r][e] = (TH1F*) convTrueRecPtIsoFullFromDirect->Clone(Form("fConvTrueRecPtIsoFullFromDirect_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTrueRecPtIsoFullFromDirect[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTrueRecPtIsoFullFromDirect[r][e]);
          TH1F *convTrueRecPtIsoCellFromDirect = new TH1F(Form("fConvTrueRecPtIsoCellFromDirect_%i_%i",r,e), Form("conversion photons with Cell ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtIsoCellFromDirect[r][e] = (TH1F*) convTrueRecPtIsoCellFromDirect->Clone(Form("fConvTrueRecPtIsoCellFromDirect_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTrueRecPtIsoCellFromDirect[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTrueRecPtIsoCellFromDirect[r][e]);


          TH1F *convTrueRecPtTaggedCaloIsoNeutral = new TH1F(Form("fConvTrueRecPtTaggedCaloIsoNeutral_%i_%i",r,e), Form("conversion photons with Neutral ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtTaggedCaloIsoNeutral[r][e] = (TH1F*) convTrueRecPtTaggedCaloIsoNeutral->Clone(Form("fConvTrueRecPtTaggedCaloIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTrueRecPtTaggedCaloIsoNeutral[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTrueRecPtTaggedCaloIsoNeutral[r][e]);
          TH1F *convTrueRecPtTaggedCaloIsoFull = new TH1F(Form("fConvTrueRecPtTaggedCaloIsoFull_%i_%i",r,e), Form("conversion photons with Full ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtTaggedCaloIsoFull[r][e] = (TH1F*) convTrueRecPtTaggedCaloIsoFull->Clone(Form("fConvTrueRecPtTaggedCaloIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTrueRecPtTaggedCaloIsoFull[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTrueRecPtTaggedCaloIsoFull[r][e]);
          TH1F *convTrueRecPtTaggedCaloIsoCell = new TH1F(Form("fConvTrueRecPtTaggedCaloIsoCell_%i_%i",r,e), Form("conversion photons with Cell ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtTaggedCaloIsoCell[r][e] = (TH1F*) convTrueRecPtTaggedCaloIsoCell->Clone(Form("fConvTrueRecPtTaggedCaloIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTrueRecPtTaggedCaloIsoCell[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTrueRecPtTaggedCaloIsoCell[r][e]);

          TH1F *convTrueRecPtMCIsoNeutral = new TH1F(Form("fConvTrueRecPtMCIsoNeutral_%i_%i",r,e), Form("conversion photons with Neutral ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtMCIsoNeutral[r][e] = (TH1F*) convTrueRecPtMCIsoNeutral->Clone(Form("fConvTrueRecPtMCIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTrueRecPtMCIsoNeutral[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTrueRecPtMCIsoNeutral[r][e]);
          TH1F *convTrueRecPtMCIsoFull = new TH1F(Form("fConvTrueRecPtMCIsoFull_%i_%i",r,e), Form("conversion photons with Full ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtMCIsoFull[r][e] = (TH1F*) convTrueRecPtMCIsoFull->Clone(Form("fConvTrueRecPtMCIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTrueRecPtMCIsoFull[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTrueRecPtMCIsoFull[r][e]);
          TH1F *convTrueRecPtMCIsoCell = new TH1F(Form("fConvTrueRecPtMCIsoCell_%i_%i",r,e), Form("conversion photons with Cell ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtMCIsoCell[r][e] = (TH1F*) convTrueRecPtMCIsoCell->Clone(Form("fConvTrueRecPtMCIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTrueRecPtMCIsoCell[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTrueRecPtMCIsoCell[r][e]);

          TH1F *convTrueRecPtTaggedCaloMCIsoNeutral = new TH1F(Form("fConvTrueRecPtTaggedCaloMCIsoNeutral_%i_%i",r,e), Form("conversion photons with Neutral ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtTaggedCaloMCIsoNeutral[r][e] = (TH1F*) convTrueRecPtTaggedCaloMCIsoNeutral->Clone(Form("fConvTrueRecPtTaggedCaloMCIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTrueRecPtTaggedCaloMCIsoNeutral[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTrueRecPtTaggedCaloMCIsoNeutral[r][e]);
          TH1F *convTrueRecPtTaggedCaloMCIsoFull = new TH1F(Form("fConvTrueRecPtTaggedCaloMCIsoFull_%i_%i",r,e), Form("conversion photons with Full ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtTaggedCaloMCIsoFull[r][e] = (TH1F*) convTrueRecPtTaggedCaloMCIsoFull->Clone(Form("fConvTrueRecPtTaggedCaloMCIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTrueRecPtTaggedCaloMCIsoFull[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTrueRecPtTaggedCaloMCIsoFull[r][e]);
          TH1F *convTrueRecPtTaggedCaloMCIsoCell = new TH1F(Form("fConvTrueRecPtTaggedCaloMCIsoCell_%i_%i",r,e), Form("conversion photons with Cell ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtTaggedCaloMCIsoCell[r][e] = (TH1F*) convTrueRecPtTaggedCaloMCIsoCell->Clone(Form("fConvTrueRecPtTaggedCaloMCIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvTrueRecPtTaggedCaloMCIsoCell[r][e]->Sumw2();
          fConvFolderTrue->Add(fConvTrueRecPtTaggedCaloMCIsoCell[r][e]);
        }

      }
    }

    // Inv mass histos
    fConvInvMass = new TH2F("fConvInvMass","fConvInvMass;m_{#gamma #gamma} (GeV/c^{2});conv p_{T} (GeV/c)",nMassBins,minMass,maxMass,nPtBins,minPt,maxPt);
    fConvFolderRec->Add(fConvInvMass);
    if(fIsMC>0){
       fConvTrueInvMass = new TH2F("fConvTrueInvMass","fConvTrueInvMass;m_{#gamma #gamma} (GeV/c^{2});conv p_{T} (GeV/c)",nMassBins,minMass,maxMass,nPtBins,minPt,maxPt);
       fConvTrueInvMass_FromDecay = new TH2F("fConvTrueInvMass_FromDecay","fConvTrueInvMass_FromDecay;m_{#gamma #gamma} (GeV/c^{2});conv p_{T} (GeV/c)",nMassBins,minMass,maxMass,nPtBins,minPt,maxPt);
       fConvTrueInvMass_FromDirect = new TH2F("fConvTrueInvMass_FromDirect","fConvTrueInvMass_FromDirect;m_{#gamma #gamma} (GeV/c^{2});conv p_{T} (GeV/c)",nMassBins,minMass,maxMass,nPtBins,minPt,maxPt);
      fConvTrueInvMass->Sumw2();
       fConvTrueInvMass_FromDecay->Sumw2();
       fConvTrueInvMass_FromDirect->Sumw2();
       fConvFolderTrue->Add(fConvTrueInvMass);
       fConvFolderTrue->Add(fConvTrueInvMass_FromDecay);
       fConvFolderTrue->Add(fConvTrueInvMass_FromDirect);
    }
    for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
    {
        TH2F *convInvMassAntiChargedIsolated = new TH2F(Form("convInvMassAntiChargedIsolated_R%1.1f",fTrackIsolationR.at(r)),Form("fConvInvMassAntiChargedIsolated_R%1.1f;m_{#gamma #gamma} (GeV/c^{2});conv p_{T} (GeV/c)",fTrackIsolationR.at(r)),nMassBins,minMass,maxMass,nPtBins,minPt,maxPt);
        fConvInvMassAntiChargedIsolated[r] = (TH2F*) convInvMassAntiChargedIsolated->Clone(Form("fConvInvMassAntiChargedIsolated_R%1.1f",fTrackIsolationR.at(r)));
        fConvInvMassAntiChargedIsolated[r]->Sumw2();
        fConvFolderRec->Add(fConvInvMassAntiChargedIsolated[r]);
        for (UInt_t e = 0; e < fTrackIsolationE.size(); e++)
        {
            TH2F *convInvMassChargedIsolated = new TH2F(Form("convInvMassChargedIsolated_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)),Form("fConvInvMassChargedIsolated_R%1.1f_E%1.1f;m_{#gamma #gamma} (GeV/c^{2});conv p_{T} (GeV/c)",fTrackIsolationR.at(r),fTrackIsolationE.at(e)),nMassBins,minMass,maxMass,nPtBins,minPt,maxPt);
            fConvInvMassChargedIsolated[r][e] = (TH2F*) convInvMassChargedIsolated->Clone(Form("fConvInvMassChargedIsolated_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
            fConvInvMassChargedIsolated[r][e]->Sumw2();
            fConvFolderRec->Add(fConvInvMassChargedIsolated[r][e]);
        }
    }
    for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
    {
        TH2F *convInvMassAntiNeutralIsolated = new TH2F(Form("convInvMassAntiNeutralIsolated_R%1.1f",fNeutralIsolationR.at(r)),Form("fConvInvMassAntiNeutralIsolated_R%1.1f;m_{#gamma #gamma} (GeV/c^{2});conv p_{T} (GeV/c)",fNeutralIsolationR.at(r)),nMassBins,minMass,maxMass,nPtBins,minPt,maxPt);
        fConvInvMassAntiNeutralIsolated[r] = (TH2F*) convInvMassAntiNeutralIsolated->Clone(Form("fConvInvMassAntiNeutralIsolated_R%1.1f",fNeutralIsolationR.at(r)));
        fConvInvMassAntiNeutralIsolated[r]->Sumw2();
        fConvFolderRec->Add(fConvInvMassAntiNeutralIsolated[r]);

        TH2F *convInvMassAntiCellIsolated = new TH2F(Form("convInvMassAntiCellIsolated_R%1.1f",fNeutralIsolationR.at(r)),Form("fConvInvMassAntiCellIsolated_R%1.1f;m_{#gamma #gamma} (GeV/c^{2});conv p_{T} (GeV/c)",fNeutralIsolationR.at(r)),nMassBins,minMass,maxMass,nPtBins,minPt,maxPt);
        fConvInvMassAntiCellIsolated[r] = (TH2F*) convInvMassAntiCellIsolated->Clone(Form("fConvInvMassAntiCellIsolated_R%1.1f",fNeutralIsolationR.at(r)));
        fConvInvMassAntiCellIsolated[r]->Sumw2();
        fConvFolderRec->Add(fConvInvMassAntiCellIsolated[r]);

        TH2F *convInvMassAntiFullIsolated = new TH2F(Form("convInvMassAntiFullIsolated_R%1.1f",fNeutralIsolationR.at(r)),Form("fConvInvMassAntiFullIsolated_R%1.1f;m_{#gamma #gamma} (GeV/c^{2});conv p_{T} (GeV/c)",fNeutralIsolationR.at(r)),nMassBins,minMass,maxMass,nPtBins,minPt,maxPt);
        fConvInvMassAntiFullIsolated[r] = (TH2F*) convInvMassAntiFullIsolated->Clone(Form("fConvInvMassAntiFullIsolated_R%1.1f",fNeutralIsolationR.at(r)));
        fConvInvMassAntiFullIsolated[r]->Sumw2();
        fConvFolderRec->Add(fConvInvMassAntiFullIsolated[r]);
        for (UInt_t e = 0; e < fNeutralIsolationE.size(); e++)
        {
            TH2F *convInvMassNeutralIsolated = new TH2F(Form("convInvMassNeutralIsolated_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)),Form("fConvInvMassNeutralIsolated_R%1.1f_E%1.1f;m_{#gamma #gamma} (GeV/c^{2});conv p_{T} (GeV/c)",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)),nMassBins,minMass,maxMass,nPtBins,minPt,maxPt);
            fConvInvMassNeutralIsolated[r][e] = (TH2F*) convInvMassNeutralIsolated->Clone(Form("fConvInvMassNeutralIsolated_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
            fConvInvMassNeutralIsolated[r][e]->Sumw2();
            fConvFolderRec->Add(fConvInvMassNeutralIsolated[r][e]);

            TH2F *convInvMassCellIsolated = new TH2F(Form("convInvMassCellIsolated_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)),Form("fConvInvMassCellIsolated_R%1.1f_E%1.1f;m_{#gamma #gamma} (GeV/c^{2});conv p_{T} (GeV/c)",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)),nMassBins,minMass,maxMass,nPtBins,minPt,maxPt);
            fConvInvMassCellIsolated[r][e] = (TH2F*) convInvMassCellIsolated->Clone(Form("fConvInvMassCellIsolated_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
            fConvInvMassCellIsolated[r][e]->Sumw2();
            fConvFolderRec->Add(fConvInvMassCellIsolated[r][e]);

            TH2F *convInvMassFullIsolated = new TH2F(Form("convInvMassFullIsolated_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)),Form("fConvInvMassFullIsolated_R%1.1f_E%1.1f;m_{#gamma #gamma} (GeV/c^{2});conv p_{T} (GeV/c)",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)),nMassBins,minMass,maxMass,nPtBins,minPt,maxPt);
            fConvInvMassFullIsolated[r][e] = (TH2F*) convInvMassFullIsolated->Clone(Form("fConvInvMassFullIsolated_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
            fConvInvMassFullIsolated[r][e]->Sumw2();
            fConvFolderRec->Add(fConvInvMassFullIsolated[r][e]);
        }
    }

  }


  //
  // ─── CALO HISTOGRAMS ────────────────────────────────────────────────────────────
  //

  fCaloFolderRec          = new TList();
  fCaloFolderRec->SetName("caloPhotonsRec");
  fCaloFolderRec->SetOwner(kTRUE);
  fOutputList->Add(fCaloFolderRec);
  // always fill this for trigger rejection factor
  fCaloPt = new TH1F("fCaloPt","calo photons in EMC acc;p_{T} (GeV/c); counts",nPtBins,minPt,maxPt);
  fCaloPt->Sumw2();
  fCaloFolderRec->Add(fCaloPt);
  if(fUseHistograms){

    fCaloPtBeforeAcc = new TH1F("fCaloPtBeforeAcc", "calo photons all acc;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);

    fCaloE = new TH1F("fCaloE", "calo photons in EMC;E_{clus} (GeV); counts", nPtBins,minPt,maxPt);

    fCaloPtTaggedCalo = new TH1F("fCaloPtTaggedCalo", "calo photons that survived tagging;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
    fCaloPtTaggedAsDecayCalo = new TH1F("fCaloPtTaggedAsDecayCalo", "calo photons that survived tagging;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);


    fCaloRho = new TH2F("fCaloRho", "charged event density;#rho; counts", nPtBins,minPt,maxPt, nPtBins,minPt,maxPt);
    for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
    {
      TH2F *caloRhoTimesArea = new TH2F(Form("fCaloRhoTimesArea_%i",r), Form("Rho times area for R =  %1.1f;#rho; counts",fTrackIsolationR.at(r)), nPtBins,minPt,maxPt, nPtBins,minPt,maxPt);
      fCaloRhoTimesArea[r] = (TH2F*) caloRhoTimesArea->Clone(Form("fCaloRhoTimesArea_R%1.1f",fTrackIsolationR.at(r)));
      fCaloRhoTimesArea[r]->Sumw2();
      fCaloFolderRec->Add(fCaloRhoTimesArea[r]);

      TH2F *caloRhoTimesAreaLeft = new TH2F(Form("fCaloRhoTimesAreaLeft_%i",r), Form("Rho times area for R =  %1.1f;#rho; counts",fTrackIsolationR.at(r)), nPtBins,minPt,maxPt, nPtBins,minPt,maxPt);
      fCaloRhoTimesAreaLeft[r] = (TH2F*) caloRhoTimesAreaLeft->Clone(Form("fCaloRhoTimesAreaLeft_R%1.1f",fTrackIsolationR.at(r)));
      fCaloRhoTimesAreaLeft[r]->Sumw2();
      fCaloFolderRec->Add(fCaloRhoTimesAreaLeft[r]);

      TH2F *caloRhoTimesAreaRight = new TH2F(Form("fCaloRhoTimesAreaRight_%i",r), Form("Rho times area for R =  %1.1f;#rho; counts",fTrackIsolationR.at(r)), nPtBins,minPt,maxPt, nPtBins,minPt,maxPt);
      fCaloRhoTimesAreaRight[r] = (TH2F*) caloRhoTimesAreaRight->Clone(Form("fCaloRhoTimesAreaRight_R%1.1f",fTrackIsolationR.at(r)));
      fCaloRhoTimesAreaRight[r]->Sumw2();
      fCaloFolderRec->Add(fCaloRhoTimesAreaRight[r]);

      TH2F *caloRhoTimesAreaBack = new TH2F(Form("fCaloRhoTimesAreaBack_%i",r), Form("Rho times area for R =  %1.1f;#rho; counts",fTrackIsolationR.at(r)), nPtBins,minPt,maxPt, nPtBins,minPt,maxPt);
      fCaloRhoTimesAreaBack[r] = (TH2F*) caloRhoTimesAreaBack->Clone(Form("fCaloRhoTimesAreaBack_R%1.1f",fTrackIsolationR.at(r)));
      fCaloRhoTimesAreaBack[r]->Sumw2();
      fCaloFolderRec->Add(fCaloRhoTimesAreaBack[r]);
    }

    fCaloPtBeforeAcc->Sumw2();
    fCaloE->Sumw2();
    fCaloPtTaggedCalo->Sumw2();
    fCaloPtTaggedAsDecayCalo->Sumw2();
    fCaloRho->Sumw2();

    fCaloFolderRec->Add(fCaloPtBeforeAcc);
    fCaloFolderRec->Add(fCaloE);
    fCaloFolderRec->Add(fCaloPtTaggedCalo);
    fCaloFolderRec->Add(fCaloPtTaggedAsDecayCalo);
    fCaloFolderRec->Add(fCaloRho);

    for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
    {
      TH2F *caloIsoCharged = new TH2F(Form("caloIsoCharged_%i",r), Form("charged track ISO in R < %1.1f;#sum p_{T} (GeV/c); calo. p_{T} (GeV/c)",fTrackIsolationR.at(r)), 250,-10,50,nPtBins,minPt,maxPt);
      fCaloIsoCharged[r] = (TH2F*) caloIsoCharged->Clone(Form("fCaloIsoCharged_R%1.1f",fTrackIsolationR.at(r)));
      fCaloIsoCharged[r]->Sumw2();
      fCaloFolderRec->Add(fCaloIsoCharged[r]);

      TH2F *caloIsoRawCharged = new TH2F(Form("caloIsoRawCharged_%i",r), Form("charged track ISO in R < %1.1f;#sum p_{T} (GeV/c); calo. p_{T} (GeV/c)",fTrackIsolationR.at(r)), 250,-10,50,nPtBins,minPt,maxPt);
      fCaloIsoRawCharged[r] = (TH2F*) caloIsoRawCharged->Clone(Form("fCaloIsoRawCharged_R%1.1f",fTrackIsolationR.at(r)));
      fCaloIsoRawCharged[r]->Sumw2();
      fCaloFolderRec->Add(fCaloIsoRawCharged[r]);
    }
    for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
    {
      TH2F *caloIsoNeutral = new TH2F(Form("caloIsoNeutral_%i",r), Form("Neutral track ISO in R < %1.1f;#sum p_{T} (GeV/c); calo. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), 250,-10,50,nPtBins,minPt,maxPt);
      fCaloIsoNeutral[r] = (TH2F*) caloIsoNeutral->Clone(Form("fCaloIsoNeutral_R%1.1f",fNeutralIsolationR.at(r)));
      fCaloIsoNeutral[r]->Sumw2();
      fCaloFolderRec->Add(fCaloIsoNeutral[r]);

      TH2F *caloIsoFull = new TH2F(Form("caloIsoFull_%i",r), Form("Full track ISO in R < %1.1f;#sum p_{T} (GeV/c); calo. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), 250,-10.,50,nPtBins,minPt,maxPt);
      fCaloIsoFull[r] = (TH2F*) caloIsoFull->Clone(Form("fCaloIsoFull_R%1.1f",fNeutralIsolationR.at(r)));
      fCaloIsoFull[r]->Sumw2();
      fCaloFolderRec->Add(fCaloIsoFull[r]);

      TH2F *caloIsoCell = new TH2F(Form("caloIsoCell_%i",r), Form("Full track ISO in R < %1.1f;#sum p_{T} (GeV/c); calo. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), 250,-10.,50,nPtBins,minPt,maxPt);
      fCaloIsoCell[r] = (TH2F*) caloIsoCell->Clone(Form("fCaloIsoCell_R%1.1f",fNeutralIsolationR.at(r)));
      fCaloIsoCell[r]->Sumw2();
      fCaloFolderRec->Add(fCaloIsoCell[r]);
    }

    fCaloFolderTrue          = new TList();
    fCaloFolderTrue->SetName("caloPhotonsTrue");
    fCaloFolderTrue->SetOwner(kTRUE);
    if(fIsMC > 0){
      fOutputList->Add(fCaloFolderTrue);
      fCaloTruePt = new TH1F("fCaloTruePt", "validated calo photons in EMC acceptance;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fCaloTruePtNotProper = new TH1F("fCaloTruePtNotProper", "validated calo photons in EMC acceptance;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fCaloTruePtNoIsPrimary = new TH1F("fCaloTruePtNoIsPrimary", "validated calo photons in EMC acceptance;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fCaloTrueWithoutConvPt = new TH1F("fCaloTrueWithoutConvPt", "validated calo photons in EMC acceptance without conv a truth;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fCaloTruePtPrimary = new TH1F("fCaloTruePtPrimary", "calo photon that has not a pi0 etc. as mother;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fCaloTruePtDecay = new TH1F("fCaloTruePtDecay", "calo photon from decay;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fCaloTruePtDecayFoundOtherInCluster = new TH1F("fCaloTruePtDecayFoundOtherInCluster", "calo photon from decay, where the other decay particle was found in EMC;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fCaloTruePtDecayOtherInAcc = new TH1F("fCaloTruePtDecayOtherInAcc", "calo photon from decay, where the other decay particle was found;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fCaloTruePtDecayOtherInAccAboveMinEnergy = new TH1F("fCaloTruePtDecayOtherInAccAboveMinEnergy", "calo photon from decay, where the other decay particle is in acc. and above 0.7 GeV;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);

      // effifciency studies
      fCaloTrueSignalPtClusterCuts = new TH1F("fCaloTrueSignalPtClusterCuts", "validated photons fulfilling cluster cuts;gen p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fCaloTrueSignalPtClusterCutsOnlyPhoton = new TH1F("fCaloTrueSignalPtClusterCutsOnlyPhoton", "validated photons fulfilling cluster cuts;gen p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fCaloTrueSignalPtClusterCutsOnlyConv = new TH1F("fCaloTrueSignalPtClusterCutsOnlyConv", "validated photons fulfilling cluster cuts;gen p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fCaloTrueSignalPtClusterCutsM02Cuts = new TH1F("fCaloTrueSignalPtClusterCutsM02Cuts", "validated photons fulfilling cluster cuts;gen p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      
      fCaloTrueSignalRecPtClusterCuts = new TH1F("fCaloTrueSignalRecPtClusterCuts", "validated photons fulfilling cluster cuts;rec p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fCaloTrueSignalRecPtClusterCutsOnlyPhoton = new TH1F("fCaloTrueSignalRecPtClusterCutsOnlyPhoton", "validated photons fulfilling cluster cuts;rec p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fCaloTrueSignalRecPtClusterCutsOnlyConv = new TH1F("fCaloTrueSignalRecPtClusterCutsOnlyConv", "validated photons fulfilling cluster cuts;rec p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fCaloTrueSignalRecPtClusterCutsM02Cuts = new TH1F("fCaloTrueSignalRecPtClusterCutsM02Cuts", "validated photons fulfilling cluster cuts;rec p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
        
      fCaloTrueSignalPtClusterCuts->Sumw2();
      fCaloTrueSignalPtClusterCutsOnlyPhoton->Sumw2();
      fCaloTrueSignalPtClusterCutsOnlyConv->Sumw2();
      fCaloTrueSignalPtClusterCutsM02Cuts->Sumw2();
      
      fCaloTrueSignalRecPtClusterCuts->Sumw2();
      fCaloTrueSignalRecPtClusterCutsOnlyPhoton->Sumw2();
      fCaloTrueSignalRecPtClusterCutsOnlyConv->Sumw2();
      fCaloTrueSignalRecPtClusterCutsM02Cuts->Sumw2();

      fCaloFolderTrue->Add(fCaloTrueSignalPtClusterCuts);
      fCaloFolderTrue->Add(fCaloTrueSignalPtClusterCutsOnlyPhoton);
      fCaloFolderTrue->Add(fCaloTrueSignalPtClusterCutsOnlyConv);
      fCaloFolderTrue->Add(fCaloTrueSignalPtClusterCutsM02Cuts);
      
      fCaloFolderTrue->Add(fCaloTrueSignalRecPtClusterCuts);
      fCaloFolderTrue->Add(fCaloTrueSignalRecPtClusterCutsOnlyPhoton);
      fCaloFolderTrue->Add(fCaloTrueSignalRecPtClusterCutsOnlyConv);
      fCaloFolderTrue->Add(fCaloTrueSignalRecPtClusterCutsM02Cuts);
  
      for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
      {
        for (UInt_t e = 0; e < fTrackIsolationE.size(); e++)
        {
          TH1F *caloTrueSignalPtClusterCutsM02CutsIsoCuts = new TH1F(Form("fCaloTrueSignalPtClusterCutsM02CutsIsoCuts_%i_%i",r,e), Form("calo photons with charged track ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueSignalPtClusterCutsM02CutsIsoCuts[r][e] = (TH1F*) caloTrueSignalPtClusterCutsM02CutsIsoCuts->Clone(Form("fCaloTrueSignalPtClusterCutsM02CutsIsoCuts_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fCaloTrueSignalPtClusterCutsM02CutsIsoCuts[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTrueSignalPtClusterCutsM02CutsIsoCuts[r][e]);

          TH1F *caloTrueSignalPtClusterCutsM02CutsIsoCutsOnlyPhoton = new TH1F(Form("fCaloTrueSignalPtClusterCutsM02CutsIsoCutsOnlyPhoton_%i_%i",r,e), Form("calo photons with charged track ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueSignalPtClusterCutsM02CutsIsoCutsOnlyPhoton[r][e] = (TH1F*) caloTrueSignalPtClusterCutsM02CutsIsoCutsOnlyPhoton->Clone(Form("fCaloTrueSignalPtClusterCutsM02CutsIsoCutsOnlyPhoton_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fCaloTrueSignalPtClusterCutsM02CutsIsoCutsOnlyPhoton[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTrueSignalPtClusterCutsM02CutsIsoCutsOnlyPhoton[r][e]);

          TH1F *caloTrueSignalPtClusterCutsM02CutsIsoCutsOnlyConv = new TH1F(Form("fCaloTrueSignalPtClusterCutsM02CutsIsoCutsOnlyConv_%i_%i",r,e), Form("calo photons with charged track ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueSignalPtClusterCutsM02CutsIsoCutsOnlyConv[r][e] = (TH1F*) caloTrueSignalPtClusterCutsM02CutsIsoCutsOnlyConv->Clone(Form("fCaloTrueSignalPtClusterCutsM02CutsIsoCutsOnlyConv_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fCaloTrueSignalPtClusterCutsM02CutsIsoCutsOnlyConv[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTrueSignalPtClusterCutsM02CutsIsoCutsOnlyConv[r][e]);

          TH1F *caloTrueSignalRecPtClusterCutsM02CutsIsoCuts = new TH1F(Form("fCaloTrueSignalRecPtClusterCutsM02CutsIsoCuts_%i_%i",r,e), Form("calo photons with charged track ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueSignalRecPtClusterCutsM02CutsIsoCuts[r][e] = (TH1F*) caloTrueSignalRecPtClusterCutsM02CutsIsoCuts->Clone(Form("fCaloTrueSignalRecPtClusterCutsM02CutsIsoCuts_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fCaloTrueSignalRecPtClusterCutsM02CutsIsoCuts[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTrueSignalRecPtClusterCutsM02CutsIsoCuts[r][e]);

          TH1F *caloTrueSignalRecPtClusterCutsM02CutsIsoCutsOnlyPhoton = new TH1F(Form("fCaloTrueSignalRecPtClusterCutsM02CutsIsoCutsOnlyPhoton_%i_%i",r,e), Form("calo photons with charged track ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueSignalRecPtClusterCutsM02CutsIsoCutsOnlyPhoton[r][e] = (TH1F*) caloTrueSignalRecPtClusterCutsM02CutsIsoCutsOnlyPhoton->Clone(Form("fCaloTrueSignalRecPtClusterCutsM02CutsIsoCutsOnlyPhoton_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fCaloTrueSignalRecPtClusterCutsM02CutsIsoCutsOnlyPhoton[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTrueSignalRecPtClusterCutsM02CutsIsoCutsOnlyPhoton[r][e]);

          TH1F *caloTrueSignalRecPtClusterCutsM02CutsIsoCutsOnlyConv = new TH1F(Form("fCaloTrueSignalRecPtClusterCutsM02CutsIsoCutsOnlyConv_%i_%i",r,e), Form("calo photons with charged track ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueSignalRecPtClusterCutsM02CutsIsoCutsOnlyConv[r][e] = (TH1F*) caloTrueSignalRecPtClusterCutsM02CutsIsoCutsOnlyConv->Clone(Form("fCaloTrueSignalRecPtClusterCutsM02CutsIsoCutsOnlyConv_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fCaloTrueSignalRecPtClusterCutsM02CutsIsoCutsOnlyConv[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTrueSignalRecPtClusterCutsM02CutsIsoCutsOnlyConv[r][e]);
        }
      }

      // iso and tagging studies
      fCaloTruePtTaggedCalo = new TH1F("fCaloTruePtTaggedCalo", "calo photons that survived tagging;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fCaloTruePtTaggedAsDecayCalo = new TH1F("fCaloTruePtTaggedAsDecayCalo", "calo photons that survived tagging;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);

      // True with rec pT
      fCaloTrueRecPt = new TH1F("fCaloTrueRecPt", "validated calo photons in EMC acceptance;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fCaloTrueRecPtNoIsPrimary = new TH1F("fCaloTrueRecPtNoIsPrimary", "validated calo photons in EMC acceptance;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fCaloTrueRecPtvsTruePt = new TH2F("fCaloTrueRecPtvsTruePt", "rec. vs true;true p_{T} (GeV/c); rec. p_{T} (GeV/c);", nPtBins,minPt,maxPt, nPtBins,minPt,maxPt);
      fCaloTrueRecPtvsTruePtNotProper = new TH2F("fCaloTrueRecPtvsTruePtNotProper", "rec. vs true;true p_{T} (GeV/c); rec. p_{T} (GeV/c);", nPtBins,minPt,maxPt, nPtBins,minPt,maxPt);
      fCaloTrueWithoutConvRecPt = new TH1F("fCaloTrueWithoutConvRecPt", "validated calo photons in EMC acceptance without conv as truth;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fCaloTrueRecPtPrimary = new TH1F("fCaloTrueRecPtPrimary", "calo photon that has not a pi0 etc. as mother;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fCaloTrueRecPtDecay = new TH1F("fCaloTrueRecPtDecay", "calo photon from decay;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fCaloTrueRecPtDecayFoundOtherInCluster = new TH1F("fCaloTrueRecPtDecayFoundOtherInCluster", "calo photon from decay, where the other decay particle was found in EMC;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fCaloTrueRecPtDecayOtherInAcc = new TH1F("fCaloTrueRecPtDecayOtherInAcc", "calo photon from decay, where the other decay particle was found;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fCaloTrueRecPtDecayOtherInAccAboveMinEnergy = new TH1F("fCaloTrueRecPtDecayOtherInAccAboveMinEnergy", "calo photon from decay, where the other decay particle is in acc. and above 0.7 GeV;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);

      // iso and tagging studies
      fCaloTrueRecPtTaggedCalo = new TH1F("fCaloTrueRecPtTaggedCalo", "calo photons that survived tagging;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);
      fCaloTrueRecPtTaggedAsDecayCalo = new TH1F("fCaloTrueRecPtTaggedAsDecayCalo", "calo photons that survived tagging;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);

      fCaloTruePt->Sumw2();
      fCaloTruePtNotProper->Sumw2();
      fCaloTruePtNoIsPrimary->Sumw2();
      fCaloTrueWithoutConvPt->Sumw2();
      fCaloTruePtPrimary->Sumw2();
      fCaloTruePtDecay->Sumw2();
      fCaloTruePtDecayFoundOtherInCluster->Sumw2();
      fCaloTruePtDecayOtherInAcc->Sumw2();
      fCaloTruePtDecayOtherInAccAboveMinEnergy->Sumw2();
      fCaloTruePtTaggedCalo->Sumw2();
      fCaloTruePtTaggedAsDecayCalo->Sumw2();
      fCaloTrueRecPt->Sumw2();
      fCaloTrueRecPtNoIsPrimary->Sumw2();
      fCaloTrueRecPtvsTruePt->Sumw2();
      fCaloTrueRecPtvsTruePtNotProper->Sumw2();
      fCaloTrueWithoutConvRecPt->Sumw2();
      fCaloTrueRecPtPrimary->Sumw2();
      fCaloTrueRecPtDecay->Sumw2();
      fCaloTrueRecPtDecayFoundOtherInCluster->Sumw2();
      fCaloTrueRecPtDecayOtherInAcc->Sumw2();
      fCaloTrueRecPtDecayOtherInAccAboveMinEnergy->Sumw2();
      fCaloTrueRecPtTaggedCalo->Sumw2();
      fCaloTrueRecPtTaggedAsDecayCalo->Sumw2();


      // add to folders
      fCaloFolderTrue->Add(fCaloTruePt);
      fCaloFolderTrue->Add(fCaloTruePtNotProper);
      fCaloFolderTrue->Add(fCaloTrueWithoutConvPt);
      fCaloFolderTrue->Add(fCaloTruePtNoIsPrimary);
      fCaloFolderTrue->Add(fCaloTruePtPrimary);
      fCaloFolderTrue->Add(fCaloTruePtDecay);
      fCaloFolderTrue->Add(fCaloTruePtDecayFoundOtherInCluster);
      fCaloFolderTrue->Add(fCaloTruePtDecayOtherInAcc);
      fCaloFolderTrue->Add(fCaloTruePtDecayOtherInAccAboveMinEnergy);

      // iso and tagging studies
      fCaloFolderTrue->Add(fCaloTruePtTaggedCalo);
      fCaloFolderTrue->Add(fCaloTruePtTaggedAsDecayCalo);

      // True with rec pT
      fCaloFolderTrue->Add(fCaloTrueRecPt);
      fCaloFolderTrue->Add(fCaloTrueRecPtNoIsPrimary);
      fCaloFolderTrue->Add(fCaloTrueRecPtvsTruePt);
      fCaloFolderTrue->Add(fCaloTrueRecPtvsTruePtNotProper);
      fCaloFolderTrue->Add(fCaloTrueWithoutConvRecPt);
      fCaloFolderTrue->Add(fCaloTrueRecPtPrimary);
      fCaloFolderTrue->Add(fCaloTrueRecPtDecay);
      fCaloFolderTrue->Add(fCaloTrueRecPtDecayFoundOtherInCluster);
      fCaloFolderTrue->Add(fCaloTrueRecPtDecayOtherInAcc);
      fCaloFolderTrue->Add(fCaloTrueRecPtDecayOtherInAccAboveMinEnergy);

      // iso and tagging studies
      fCaloFolderTrue->Add(fCaloTrueRecPtTaggedCalo);
      fCaloFolderTrue->Add(fCaloTrueRecPtTaggedAsDecayCalo);

      for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
      {
        TH2F *caloTrueIsoCharged = new TH2F(Form("fCaloTrueIsoCharged_%i",r), Form("Charged track ISO in R < %1.1f;#sum p_{T} (GeV/c); counts",fTrackIsolationR.at(r)), 100,0,50,nPtBins,minPt,maxPt);
        fCaloTrueIsoCharged[r] = (TH2F*) caloTrueIsoCharged->Clone(Form("fCaloTrueIsoCharged_R%1.1f",fTrackIsolationR.at(r)));
        fCaloTrueIsoCharged[r]->Sumw2();
        fCaloFolderTrue->Add(fCaloTrueIsoCharged[r]);

        TH2F *caloTrueIsoCharged_FromDecay = new TH2F(Form("fCaloTrueIsoCharged_FromDecay_%i",r), Form("Charged track ISO in R < %1.1f;#sum p_{T} (GeV/c); counts",fTrackIsolationR.at(r)), 100,0,50,nPtBins,minPt,maxPt);
        fCaloTrueIsoCharged_FromDecay[r] = (TH2F*) caloTrueIsoCharged_FromDecay->Clone(Form("fCaloTrueIsoCharged_FromDecay_R%1.1f",fTrackIsolationR.at(r)));
        fCaloTrueIsoCharged_FromDecay[r]->Sumw2();
        fCaloFolderTrue->Add(fCaloTrueIsoCharged_FromDecay[r]);

        TH2F *caloTrueIsoCharged_FromDirect = new TH2F(Form("fCaloTrueIsoCharged_FromDirect_%i",r), Form("Charged track ISO in R < %1.1f;#sum p_{T} (GeV/c); counts",fTrackIsolationR.at(r)), 100,0,50,nPtBins,minPt,maxPt);
        fCaloTrueIsoCharged_FromDirect[r] = (TH2F*) caloTrueIsoCharged_FromDirect->Clone(Form("fCaloTrueIsoCharged_FromDirect_R%1.1f",fTrackIsolationR.at(r)));
        fCaloTrueIsoCharged_FromDirect[r]->Sumw2();
        fCaloFolderTrue->Add(fCaloTrueIsoCharged_FromDirect[r]);

        TH2F *caloTrueIsoMCCharged = new TH2F(Form("fCaloTrueIsoMCCharged_%i",r), Form("gen level Charged track ISO in R < %1.1f;#sum p_{T} (GeV/c); counts",fTrackIsolationR.at(r)), 100,0,50,nPtBins,minPt,maxPt);
        fCaloTrueIsoMCCharged[r] = (TH2F*) caloTrueIsoMCCharged->Clone(Form("fCaloTrueIsoMCCharged_R%1.1f",fTrackIsolationR.at(r)));
        fCaloTrueIsoMCCharged[r]->Sumw2();
        fCaloFolderTrue->Add(fCaloTrueIsoMCCharged[r]);
      }

      for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
      {
        TH2F *caloTrueIsoNeutral = new TH2F(Form("fCaloTrueIsoNeutral_%i",r), Form("Neutral track ISO in R < %1.1f;#sum p_{T} (GeV/c); calo. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), 100,0,50,nPtBins,minPt,maxPt);
        fCaloTrueIsoNeutral[r] = (TH2F*) caloTrueIsoNeutral->Clone(Form("fCaloTrueIsoNeutral_R%1.1f",fNeutralIsolationR.at(r)));
        fCaloTrueIsoNeutral[r]->Sumw2();
        fCaloFolderTrue->Add(fCaloTrueIsoNeutral[r]);

        TH2F *caloTrueIsoFull = new TH2F(Form("fCaloTrueIsoFull_%i",r), Form("Full track ISO in R < %1.1f;#sum p_{T} (GeV/c); calo. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), 100,0,50,nPtBins,minPt,maxPt);
        fCaloTrueIsoFull[r] = (TH2F*) caloTrueIsoFull->Clone(Form("fCaloTrueIsoFull_R%1.1f",fNeutralIsolationR.at(r)));
        fCaloTrueIsoFull[r]->Sumw2();
        fCaloFolderTrue->Add(fCaloTrueIsoFull[r]);

        TH2F *caloTrueIsoCell = new TH2F(Form("fCaloTrueIsoCell_%i",r), Form("Full track ISO in R < %1.1f;#sum p_{T} (GeV/c); calo. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), 100,0,50,nPtBins,minPt,maxPt);
        fCaloTrueIsoCell[r] = (TH2F*) caloTrueIsoCell->Clone(Form("fCaloTrueIsoCell_R%1.1f",fNeutralIsolationR.at(r)));
        fCaloTrueIsoCell[r]->Sumw2();
        fCaloFolderTrue->Add(fCaloTrueIsoCell[r]);

        TH2F *caloTrueIsoNeutral_FromDecay = new TH2F(Form("fCaloTrueIsoNeutral_FromDecay_%i",r), Form("Neutral track ISO in R < %1.1f;#sum p_{T} (GeV/c); calo. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), 100,0,50,nPtBins,minPt,maxPt);
        fCaloTrueIsoNeutral_FromDecay[r] = (TH2F*) caloTrueIsoNeutral_FromDecay->Clone(Form("fCaloTrueIsoNeutral_FromDecay_R%1.1f",fNeutralIsolationR.at(r)));
        fCaloTrueIsoNeutral_FromDecay[r]->Sumw2();
        fCaloFolderTrue->Add(fCaloTrueIsoNeutral_FromDecay[r]);

        TH2F *caloTrueIsoFull_FromDecay = new TH2F(Form("fCaloTrueIsoFull_FromDecay_%i",r), Form("Full track ISO in R < %1.1f;#sum p_{T} (GeV/c); calo. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), 100,0,50,nPtBins,minPt,maxPt);
        fCaloTrueIsoFull_FromDecay[r] = (TH2F*) caloTrueIsoFull_FromDecay->Clone(Form("fCaloTrueIsoFull_FromDecay_R%1.1f",fNeutralIsolationR.at(r)));
        fCaloTrueIsoFull_FromDecay[r]->Sumw2();
        fCaloFolderTrue->Add(fCaloTrueIsoFull_FromDecay[r]);

        TH2F *caloTrueIsoCell_FromDecay = new TH2F(Form("fCaloTrueIsoCell_FromDecay_%i",r), Form("Full track ISO in R < %1.1f;#sum p_{T} (GeV/c); calo. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), 100,0,50,nPtBins,minPt,maxPt);
        fCaloTrueIsoCell_FromDecay[r] = (TH2F*) caloTrueIsoCell_FromDecay->Clone(Form("fCaloTrueIsoCel_FromDecay_R%1.1f",fNeutralIsolationR.at(r)));
        fCaloTrueIsoCell_FromDecay[r]->Sumw2();
        fCaloFolderTrue->Add(fCaloTrueIsoCell_FromDecay[r]);

        TH2F *caloTrueIsoNeutral_FromDirect = new TH2F(Form("fCaloTrueIsoNeutral_FromDirect_%i",r), Form("Neutral track ISO in R < %1.1f;#sum p_{T} (GeV/c); calo. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), 100,0,50,nPtBins,minPt,maxPt);
        fCaloTrueIsoNeutral_FromDirect[r] = (TH2F*) caloTrueIsoNeutral_FromDirect->Clone(Form("fCaloTrueIsoNeutral_FromDirect_R%1.1f",fNeutralIsolationR.at(r)));
        fCaloTrueIsoNeutral_FromDirect[r]->Sumw2();
        fCaloFolderTrue->Add(fCaloTrueIsoNeutral_FromDirect[r]);

        TH2F *caloTrueIsoFull_FromDirect = new TH2F(Form("fCaloTrueIsoFull_FromDirect_%i",r), Form("Full track ISO in R < %1.1f;#sum p_{T} (GeV/c); calo. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), 100,0,50,nPtBins,minPt,maxPt);
        fCaloTrueIsoFull_FromDirect[r] = (TH2F*) caloTrueIsoFull_FromDirect->Clone(Form("fCaloTrueIsoFull_FromDirect_R%1.1f",fNeutralIsolationR.at(r)));
        fCaloTrueIsoFull_FromDirect[r]->Sumw2();
        fCaloFolderTrue->Add(fCaloTrueIsoFull_FromDirect[r]);

        TH2F *caloTrueIsoCell_FromDirect = new TH2F(Form("fCaloTrueIsoCell_FromDirect_%i",r), Form("Full track ISO in R < %1.1f;#sum p_{T} (GeV/c); calo. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), 100,0,50,nPtBins,minPt,maxPt);
        fCaloTrueIsoCell_FromDirect[r] = (TH2F*) caloTrueIsoCell_FromDirect->Clone(Form("fCaloTrueIsoCel_FromDirect_R%1.1f",fNeutralIsolationR.at(r)));
        fCaloTrueIsoCell_FromDirect[r]->Sumw2();
        fCaloFolderTrue->Add(fCaloTrueIsoCell_FromDirect[r]);
      }
    }
    for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
    {
      for (UInt_t e = 0; e < fTrackIsolationE.size(); e++)
      {
        TH1F *caloPtIsoCharged = new TH1F(Form("fCaloPtIsoCharged_%i_%i",r,e), Form("calo photons with charged track ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
        fCaloPtIsoCharged[r][e] = (TH1F*) caloPtIsoCharged->Clone(Form("fCaloPtIsoCharged_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
        fCaloPtIsoCharged[r][e]->Sumw2();
        fCaloFolderRec->Add(fCaloPtIsoCharged[r][e]);
        TH1F *caloPtTaggedCaloIsoCharged = new TH1F(Form("fCaloPtTaggedCaloIsoCharged_%i_%i",r,e), Form("calo photons with charged track ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
        fCaloPtTaggedCaloIsoCharged[r][e] = (TH1F*) caloPtTaggedCaloIsoCharged->Clone(Form("fCaloPtTaggedCaloIsoCharged_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
        fCaloPtTaggedCaloIsoCharged[r][e]->Sumw2();
        fCaloFolderRec->Add(fCaloPtTaggedCaloIsoCharged[r][e]);

        if(fIsMC>0){
          // true pT

          TH1F *caloTruePtIsoCharged = new TH1F(Form("fCaloTruePtIsoCharged_%i_%i",r,e), Form("calo photons with charged track ISO < %1.1f GeV in R < %1.1f;gen. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtIsoCharged[r][e] = (TH1F*) caloTruePtIsoCharged->Clone(Form("fCaloTruePtIsoCharged_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fCaloTruePtIsoCharged[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTruePtIsoCharged[r][e]);
          TH1F *caloTruePtTaggedCaloIsoCharged = new TH1F(Form("fCaloTruePtTaggedCaloIsoCharged_%i_%i",r,e), Form("calo photons with charged track ISO < %1.1f GeV in R < %1.1f + tagging;gen. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtTaggedCaloIsoCharged[r][e] = (TH1F*) caloTruePtTaggedCaloIsoCharged->Clone(Form("fCaloTruePtTaggedCaloIsoCharged_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fCaloTruePtTaggedCaloIsoCharged[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTruePtTaggedCaloIsoCharged[r][e]);

          TH1F *caloTruePtIsoChargedFromDirect = new TH1F(Form("fCaloTruePtIsoChargedFromDirect_%i_%i",r,e), Form("calo photons with chargedFromDirect track ISO < %1.1f GeV in R < %1.1f;gen. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtIsoChargedFromDirect[r][e] = (TH1F*) caloTruePtIsoChargedFromDirect->Clone(Form("fCaloTruePtIsoChargedFromDirect_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fCaloTruePtIsoChargedFromDirect[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTruePtIsoChargedFromDirect[r][e]);

          TH1F *caloTruePtIsoChargedAndMCIsoChargedFromDirect = new TH1F(Form("fCaloTruePtIsoChargedAndMCIsoChargedFromDirect_%i_%i",r,e), Form("calo photons with chargedFromDirect track ISO < %1.1f GeV in R < %1.1f;gen. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtIsoChargedAndMCIsoChargedFromDirect[r][e] = (TH1F*) caloTruePtIsoChargedAndMCIsoChargedFromDirect->Clone(Form("fCaloTruePtIsoChargedAndMCIsoChargedFromDirect_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fCaloTruePtIsoChargedAndMCIsoChargedFromDirect[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTruePtIsoChargedAndMCIsoChargedFromDirect[r][e]);

          TH1F *caloTruePtTaggedCaloIsoChargedFromDirect = new TH1F(Form("fCaloTruePtTaggedCaloIsoChargedFromDirect_%i_%i",r,e), Form("calo photons with chargedFromDirect track ISO < %1.1f GeV in R < %1.1f + tagging;gen. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtTaggedCaloIsoChargedFromDirect[r][e] = (TH1F*) caloTruePtTaggedCaloIsoChargedFromDirect->Clone(Form("fCaloTruePtTaggedCaloIsoChargedFromDirect_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fCaloTruePtTaggedCaloIsoChargedFromDirect[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTruePtTaggedCaloIsoChargedFromDirect[r][e]);


          TH1F *caloTruePtMCIsoCharged = new TH1F(Form("fCaloTruePtMCIsoCharged_%i_%i",r,e), Form("calo photons with charged track ISO < %1.1f GeV in R < %1.1f;gen. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtMCIsoCharged[r][e] = (TH1F*) caloTruePtMCIsoCharged->Clone(Form("fCaloTruePtMCIsoCharged_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fCaloTruePtMCIsoCharged[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTruePtMCIsoCharged[r][e]);

          TH1F *caloTruePtTaggedCaloMCIsoCharged = new TH1F(Form("fCaloTruePtTaggedCaloMCIsoCharged_%i_%i",r,e), Form("calo photons with charged track ISO < %1.1f GeV in R < %1.1f + tagging;gen. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtTaggedCaloMCIsoCharged[r][e] = (TH1F*) caloTruePtTaggedCaloMCIsoCharged->Clone(Form("fCaloTruePtTaggedCaloMCIsoCharged_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fCaloTruePtTaggedCaloMCIsoCharged[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTruePtTaggedCaloMCIsoCharged[r][e]);

          // rec Pt
          TH1F *caloTrueRecPtIsoCharged = new TH1F(Form("fCaloTrueRecPtIsoCharged_%i_%i",r,e), Form("calo photons with charged track ISO < %1.1f GeV in R < %1.1f;rec. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtIsoCharged[r][e] = (TH1F*) caloTrueRecPtIsoCharged->Clone(Form("fCaloTrueRecPtIsoCharged_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fCaloTrueRecPtIsoCharged[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTrueRecPtIsoCharged[r][e]);

          TH1F *caloTrueRecPtIsoChargedFromDirect = new TH1F(Form("fCaloTrueRecPtIsoChargedFromDirect_%i_%i",r,e), Form("calo photons with charged track ISO < %1.1f GeV in R < %1.1f;rec. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtIsoChargedFromDirect[r][e] = (TH1F*) caloTrueRecPtIsoChargedFromDirect->Clone(Form("fCaloTrueRecPtIsoChargedFromDirect_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fCaloTrueRecPtIsoChargedFromDirect[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTrueRecPtIsoChargedFromDirect[r][e]);

          TH1F *caloTrueRecPtIsoChargedAndMCIsoChargedFromDirect = new TH1F(Form("fCaloTrueRecPtIsoChargedAndMCIsoChargedFromDirect_%i_%i",r,e), Form("calo photons with charged track ISO < %1.1f GeV in R < %1.1f;rec. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtIsoChargedAndMCIsoChargedFromDirect[r][e] = (TH1F*) caloTrueRecPtIsoChargedAndMCIsoChargedFromDirect->Clone(Form("fCaloTrueRecPtIsoChargedAndMCIsoChargedFromDirect_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fCaloTrueRecPtIsoChargedAndMCIsoChargedFromDirect[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTrueRecPtIsoChargedAndMCIsoChargedFromDirect[r][e]);

          TH1F *caloTrueRecPtTaggedCaloIsoCharged = new TH1F(Form("fCaloTrueRecPtTaggedCaloIsoCharged_%i_%i",r,e), Form("calo photons with charged track ISO < %1.1f GeV in R < %1.1f + tagging;rec. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtTaggedCaloIsoCharged[r][e] = (TH1F*) caloTrueRecPtTaggedCaloIsoCharged->Clone(Form("fCaloTrueRecPtTaggedCaloIsoCharged_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fCaloTrueRecPtTaggedCaloIsoCharged[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTrueRecPtTaggedCaloIsoCharged[r][e]);

          TH1F *caloTrueRecPtMCIsoCharged = new TH1F(Form("fCaloTrueRecPtMCIsoCharged_%i_%i",r,e), Form("calo photons with charged track ISO < %1.1f GeV in R < %1.1f;rec. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtMCIsoCharged[r][e] = (TH1F*) caloTrueRecPtMCIsoCharged->Clone(Form("fCaloTrueRecPtMCIsoCharged_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fCaloTrueRecPtMCIsoCharged[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTrueRecPtMCIsoCharged[r][e]);
          TH1F *caloTrueRecPtTaggedCaloMCIsoCharged = new TH1F(Form("fCaloTrueRecPtTaggedCaloMCIsoCharged_%i_%i",r,e), Form("calo photons with charged track ISO < %1.1f GeV in R < %1.1f + tagging;rec. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtTaggedCaloMCIsoCharged[r][e] = (TH1F*) caloTrueRecPtTaggedCaloMCIsoCharged->Clone(Form("fCaloTrueRecPtTaggedCaloMCIsoCharged_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fCaloTrueRecPtTaggedCaloMCIsoCharged[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTrueRecPtTaggedCaloMCIsoCharged[r][e]);
        }

      }

    }
    for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
    {
      for (UInt_t e = 0; e < fNeutralIsolationE.size(); e++)
      {
        TH1F *caloPtIsoNeutral = new TH1F(Form("fCaloPtIsoNeutral_%i_%i",r,e), Form("calo photons with Neutral ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
        fCaloPtIsoNeutral[r][e] = (TH1F*) caloPtIsoNeutral->Clone(Form("fCaloPtIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
        fCaloPtIsoNeutral[r][e]->Sumw2();
        fCaloFolderRec->Add(fCaloPtIsoNeutral[r][e]);
        TH1F *caloPtIsoFull = new TH1F(Form("fCaloPtIsoFull_%i_%i",r,e), Form("calo photons with Full ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
        fCaloPtIsoFull[r][e] = (TH1F*) caloPtIsoFull->Clone(Form("fCaloPtIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
        fCaloPtIsoFull[r][e]->Sumw2();
        fCaloFolderRec->Add(fCaloPtIsoFull[r][e]);
        TH1F *caloPtIsoCell = new TH1F(Form("fCaloPtIsoCell_%i_%i",r,e), Form("calo photons with Cell ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
        fCaloPtIsoCell[r][e] = (TH1F*) caloPtIsoCell->Clone(Form("fCaloPtIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
        fCaloPtIsoCell[r][e]->Sumw2();
        fCaloFolderRec->Add(fCaloPtIsoCell[r][e]);

        TH1F *caloPtTaggedCaloIsoNeutral = new TH1F(Form("fCaloPtTaggedCaloIsoNeutral_%i_%i",r,e), Form("calo photons with Neutral ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
        fCaloPtTaggedCaloIsoNeutral[r][e] = (TH1F*) caloPtTaggedCaloIsoNeutral->Clone(Form("fCaloPtTaggedCaloIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
        fCaloPtTaggedCaloIsoNeutral[r][e]->Sumw2();
        fCaloFolderRec->Add(fCaloPtTaggedCaloIsoNeutral[r][e]);
        TH1F *caloPtTaggedCaloIsoFull = new TH1F(Form("fCaloPtTaggedCaloIsoFull_%i_%i",r,e), Form("calo photons with Full ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
        fCaloPtTaggedCaloIsoFull[r][e] = (TH1F*) caloPtTaggedCaloIsoFull->Clone(Form("fCaloPtTaggedCaloIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
        fCaloPtTaggedCaloIsoFull[r][e]->Sumw2();
        fCaloFolderRec->Add(fCaloPtTaggedCaloIsoFull[r][e]);
        TH1F *caloPtTaggedCaloIsoCell = new TH1F(Form("fCaloPtTaggedCaloIsoCell_%i_%i",r,e), Form("calo photons with Cell ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
        fCaloPtTaggedCaloIsoCell[r][e] = (TH1F*) caloPtTaggedCaloIsoCell->Clone(Form("fCaloPtTaggedCaloIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
        fCaloPtTaggedCaloIsoCell[r][e]->Sumw2();
        fCaloFolderRec->Add(fCaloPtTaggedCaloIsoCell[r][e]);

        // True Pt
        if(fIsMC>0){
          TH1F *caloTruePtIsoNeutral = new TH1F(Form("fCaloTruePtIsoNeutral_%i_%i",r,e), Form("calo photons with Neutral ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtIsoNeutral[r][e] = (TH1F*) caloTruePtIsoNeutral->Clone(Form("fCaloTruePtIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTruePtIsoNeutral[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTruePtIsoNeutral[r][e]);
          TH1F *caloTruePtIsoFull = new TH1F(Form("fCaloTruePtIsoFull_%i_%i",r,e), Form("calo photons with Full ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtIsoFull[r][e] = (TH1F*) caloTruePtIsoFull->Clone(Form("fCaloTruePtIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTruePtIsoFull[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTruePtIsoFull[r][e]);
          TH1F *caloTruePtIsoCell = new TH1F(Form("fCaloTruePtIsoCell_%i_%i",r,e), Form("calo photons with Cell ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtIsoCell[r][e] = (TH1F*) caloTruePtIsoCell->Clone(Form("fCaloTruePtIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTruePtIsoCell[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTruePtIsoCell[r][e]);

          TH1F *caloTruePtTaggedCaloIsoNeutral = new TH1F(Form("fCaloTruePtTaggedCaloIsoNeutral_%i_%i",r,e), Form("calo photons with Neutral ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtTaggedCaloIsoNeutral[r][e] = (TH1F*) caloTruePtTaggedCaloIsoNeutral->Clone(Form("fCaloTruePtTaggedCaloIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTruePtTaggedCaloIsoNeutral[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTruePtTaggedCaloIsoNeutral[r][e]);
          TH1F *caloTruePtTaggedCaloIsoFull = new TH1F(Form("fCaloTruePtTaggedCaloIsoFull_%i_%i",r,e), Form("calo photons with Full ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtTaggedCaloIsoFull[r][e] = (TH1F*) caloTruePtTaggedCaloIsoFull->Clone(Form("fCaloTruePtTaggedCaloIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTruePtTaggedCaloIsoFull[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTruePtTaggedCaloIsoFull[r][e]);
          TH1F *caloTruePtTaggedCaloIsoCell = new TH1F(Form("fCaloTruePtTaggedCaloIsoCell_%i_%i",r,e), Form("calo photons with Cell ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtTaggedCaloIsoCell[r][e] = (TH1F*) caloTruePtTaggedCaloIsoCell->Clone(Form("fCaloTruePtTaggedCaloIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTruePtTaggedCaloIsoCell[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTruePtTaggedCaloIsoCell[r][e]);

          // from direct
          TH1F *caloTruePtIsoNeutralFromDirect = new TH1F(Form("fCaloTruePtIsoNeutralFromDirect_%i_%i",r,e), Form("calo photons with Neutral ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtIsoNeutralFromDirect[r][e] = (TH1F*) caloTruePtIsoNeutralFromDirect->Clone(Form("fCaloTruePtIsoNeutralFromDirect_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTruePtIsoNeutralFromDirect[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTruePtIsoNeutralFromDirect[r][e]);
          TH1F *caloTruePtIsoFullFromDirect = new TH1F(Form("fCaloTruePtIsoFullFromDirect_%i_%i",r,e), Form("calo photons with Full ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtIsoFullFromDirect[r][e] = (TH1F*) caloTruePtIsoFullFromDirect->Clone(Form("fCaloTruePtIsoFullFromDirect_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTruePtIsoFullFromDirect[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTruePtIsoFullFromDirect[r][e]);
          TH1F *caloTruePtIsoCellFromDirect = new TH1F(Form("fCaloTruePtIsoCellFromDirect_%i_%i",r,e), Form("calo photons with Cell ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtIsoCellFromDirect[r][e] = (TH1F*) caloTruePtIsoCellFromDirect->Clone(Form("fCaloTruePtIsoCellFromDirect_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTruePtIsoCellFromDirect[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTruePtIsoCellFromDirect[r][e]);

          TH1F *caloTruePtTaggedCaloIsoNeutralFromDirect = new TH1F(Form("fCaloTruePtTaggedCaloIsoNeutralFromDirect_%i_%i",r,e), Form("calo photons with Neutral ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtTaggedCaloIsoNeutralFromDirect[r][e] = (TH1F*) caloTruePtTaggedCaloIsoNeutralFromDirect->Clone(Form("fCaloTruePtTaggedCaloIsoNeutralFromDirect_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTruePtTaggedCaloIsoNeutralFromDirect[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTruePtTaggedCaloIsoNeutralFromDirect[r][e]);
          TH1F *caloTruePtTaggedCaloIsoFullFromDirect = new TH1F(Form("fCaloTruePtTaggedCaloIsoFullFromDirect_%i_%i",r,e), Form("calo photons with Full ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtTaggedCaloIsoFullFromDirect[r][e] = (TH1F*) caloTruePtTaggedCaloIsoFullFromDirect->Clone(Form("fCaloTruePtTaggedCaloIsoFullFromDirect_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTruePtTaggedCaloIsoFullFromDirect[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTruePtTaggedCaloIsoFullFromDirect[r][e]);
          TH1F *caloTruePtTaggedCaloIsoCellFromDirect = new TH1F(Form("fCaloTruePtTaggedCaloIsoCellFromDirect_%i_%i",r,e), Form("calo photons with Cell ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtTaggedCaloIsoCellFromDirect[r][e] = (TH1F*) caloTruePtTaggedCaloIsoCellFromDirect->Clone(Form("fCaloTruePtTaggedCaloIsoCellFromDirect_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTruePtTaggedCaloIsoCellFromDirect[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTruePtTaggedCaloIsoCellFromDirect[r][e]);

          // mc iso
          TH1F *caloTruePtMCIsoNeutral = new TH1F(Form("fCaloTruePtMCIsoNeutral_%i_%i",r,e), Form("calo photons with Neutral ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtMCIsoNeutral[r][e] = (TH1F*) caloTruePtMCIsoNeutral->Clone(Form("fCaloTruePtMCIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTruePtMCIsoNeutral[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTruePtMCIsoNeutral[r][e]);
          TH1F *caloTruePtMCIsoFull = new TH1F(Form("fCaloTruePtMCIsoFull_%i_%i",r,e), Form("calo photons with Full ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtMCIsoFull[r][e] = (TH1F*) caloTruePtMCIsoFull->Clone(Form("fCaloTruePtMCIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTruePtMCIsoFull[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTruePtMCIsoFull[r][e]);
          TH1F *caloTruePtMCIsoCell = new TH1F(Form("fCaloTruePtMCIsoCell_%i_%i",r,e), Form("calo photons with Cell ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtMCIsoCell[r][e] = (TH1F*) caloTruePtMCIsoCell->Clone(Form("fCaloTruePtMCIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTruePtMCIsoCell[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTruePtMCIsoCell[r][e]);

          TH1F *caloTruePtTaggedCaloMCIsoNeutral = new TH1F(Form("fCaloTruePtTaggedCaloMCIsoNeutral_%i_%i",r,e), Form("calo photons with Neutral ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtTaggedCaloMCIsoNeutral[r][e] = (TH1F*) caloTruePtTaggedCaloMCIsoNeutral->Clone(Form("fCaloTruePtTaggedCaloMCIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTruePtTaggedCaloMCIsoNeutral[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTruePtTaggedCaloMCIsoNeutral[r][e]);
          TH1F *caloTruePtTaggedCaloMCIsoFull = new TH1F(Form("fCaloTruePtTaggedCaloMCIsoFull_%i_%i",r,e), Form("calo photons with Full ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtTaggedCaloMCIsoFull[r][e] = (TH1F*) caloTruePtTaggedCaloMCIsoFull->Clone(Form("fCaloTruePtTaggedCaloMCIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTruePtTaggedCaloMCIsoFull[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTruePtTaggedCaloMCIsoFull[r][e]);
          TH1F *caloTruePtTaggedCaloMCIsoCell = new TH1F(Form("fCaloTruePtTaggedCaloMCIsoCell_%i_%i",r,e), Form("calo photons with Cell ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtTaggedCaloMCIsoCell[r][e] = (TH1F*) caloTruePtTaggedCaloMCIsoCell->Clone(Form("fCaloTruePtTaggedCaloMCIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTruePtTaggedCaloMCIsoCell[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTruePtTaggedCaloMCIsoCell[r][e]);

          // rec pT

          TH1F *caloTrueRecPtIsoNeutral = new TH1F(Form("fCaloTrueRecPtIsoNeutral_%i_%i",r,e), Form("calo photons with Neutral ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtIsoNeutral[r][e] = (TH1F*) caloTrueRecPtIsoNeutral->Clone(Form("fCaloTrueRecPtIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTrueRecPtIsoNeutral[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTrueRecPtIsoNeutral[r][e]);
          TH1F *caloTrueRecPtIsoFull = new TH1F(Form("fCaloTrueRecPtIsoFull_%i_%i",r,e), Form("calo photons with Full ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtIsoFull[r][e] = (TH1F*) caloTrueRecPtIsoFull->Clone(Form("fCaloTrueRecPtIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTrueRecPtIsoFull[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTrueRecPtIsoFull[r][e]);
          TH1F *caloTrueRecPtIsoCell = new TH1F(Form("fCaloTrueRecPtIsoCell_%i_%i",r,e), Form("calo photons with Cell ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtIsoCell[r][e] = (TH1F*) caloTrueRecPtIsoCell->Clone(Form("fCaloTrueRecPtIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTrueRecPtIsoCell[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTrueRecPtIsoCell[r][e]);

          TH1F *caloTrueRecPtIsoNeutralFromDirect = new TH1F(Form("fCaloTrueRecPtIsoNeutralFromDirect_%i_%i",r,e), Form("calo photons with Neutral ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtIsoNeutralFromDirect[r][e] = (TH1F*) caloTrueRecPtIsoNeutralFromDirect->Clone(Form("fCaloTrueRecPtIsoNeutralFromDirect_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTrueRecPtIsoNeutralFromDirect[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTrueRecPtIsoNeutralFromDirect[r][e]);
          TH1F *caloTrueRecPtIsoFullFromDirect = new TH1F(Form("fCaloTrueRecPtIsoFullFromDirect_%i_%i",r,e), Form("calo photons with Full ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtIsoFullFromDirect[r][e] = (TH1F*) caloTrueRecPtIsoFullFromDirect->Clone(Form("fCaloTrueRecPtIsoFullFromDirect_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTrueRecPtIsoFullFromDirect[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTrueRecPtIsoFullFromDirect[r][e]);
          TH1F *caloTrueRecPtIsoCellFromDirect = new TH1F(Form("fCaloTrueRecPtIsoCellFromDirect_%i_%i",r,e), Form("calo photons with Cell ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtIsoCellFromDirect[r][e] = (TH1F*) caloTrueRecPtIsoCellFromDirect->Clone(Form("fCaloTrueRecPtIsoCellFromDirect_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTrueRecPtIsoCellFromDirect[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTrueRecPtIsoCellFromDirect[r][e]);


          TH1F *caloTrueRecPtTaggedCaloIsoNeutral = new TH1F(Form("fCaloTrueRecPtTaggedCaloIsoNeutral_%i_%i",r,e), Form("calo photons with Neutral ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtTaggedCaloIsoNeutral[r][e] = (TH1F*) caloTrueRecPtTaggedCaloIsoNeutral->Clone(Form("fCaloTrueRecPtTaggedCaloIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTrueRecPtTaggedCaloIsoNeutral[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTrueRecPtTaggedCaloIsoNeutral[r][e]);
          TH1F *caloTrueRecPtTaggedCaloIsoFull = new TH1F(Form("fCaloTrueRecPtTaggedCaloIsoFull_%i_%i",r,e), Form("calo photons with Full ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtTaggedCaloIsoFull[r][e] = (TH1F*) caloTrueRecPtTaggedCaloIsoFull->Clone(Form("fCaloTrueRecPtTaggedCaloIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTrueRecPtTaggedCaloIsoFull[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTrueRecPtTaggedCaloIsoFull[r][e]);
          TH1F *caloTrueRecPtTaggedCaloIsoCell = new TH1F(Form("fCaloTrueRecPtTaggedCaloIsoCell_%i_%i",r,e), Form("calo photons with Cell ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtTaggedCaloIsoCell[r][e] = (TH1F*) caloTrueRecPtTaggedCaloIsoCell->Clone(Form("fCaloTrueRecPtTaggedCaloIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTrueRecPtTaggedCaloIsoCell[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTrueRecPtTaggedCaloIsoCell[r][e]);

          TH1F *caloTrueRecPtMCIsoNeutral = new TH1F(Form("fCaloTrueRecPtMCIsoNeutral_%i_%i",r,e), Form("calo photons with Neutral ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtMCIsoNeutral[r][e] = (TH1F*) caloTrueRecPtMCIsoNeutral->Clone(Form("fCaloTrueRecPtMCIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTrueRecPtMCIsoNeutral[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTrueRecPtMCIsoNeutral[r][e]);
          TH1F *caloTrueRecPtMCIsoFull = new TH1F(Form("fCaloTrueRecPtMCIsoFull_%i_%i",r,e), Form("calo photons with Full ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtMCIsoFull[r][e] = (TH1F*) caloTrueRecPtMCIsoFull->Clone(Form("fCaloTrueRecPtMCIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTrueRecPtMCIsoFull[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTrueRecPtMCIsoFull[r][e]);
          TH1F *caloTrueRecPtMCIsoCell = new TH1F(Form("fCaloTrueRecPtMCIsoCell_%i_%i",r,e), Form("calo photons with Cell ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtMCIsoCell[r][e] = (TH1F*) caloTrueRecPtMCIsoCell->Clone(Form("fCaloTrueRecPtMCIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTrueRecPtMCIsoCell[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTrueRecPtMCIsoCell[r][e]);

          TH1F *caloTrueRecPtTaggedCaloMCIsoNeutral = new TH1F(Form("fCaloTrueRecPtTaggedCaloMCIsoNeutral_%i_%i",r,e), Form("calo photons with Neutral ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtTaggedCaloMCIsoNeutral[r][e] = (TH1F*) caloTrueRecPtTaggedCaloMCIsoNeutral->Clone(Form("fCaloTrueRecPtTaggedCaloMCIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTrueRecPtTaggedCaloMCIsoNeutral[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTrueRecPtTaggedCaloMCIsoNeutral[r][e]);
          TH1F *caloTrueRecPtTaggedCaloMCIsoFull = new TH1F(Form("fCaloTrueRecPtTaggedCaloMCIsoFull_%i_%i",r,e), Form("calo photons with Full ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtTaggedCaloMCIsoFull[r][e] = (TH1F*) caloTrueRecPtTaggedCaloMCIsoFull->Clone(Form("fCaloTrueRecPtTaggedCaloMCIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTrueRecPtTaggedCaloMCIsoFull[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTrueRecPtTaggedCaloMCIsoFull[r][e]);
          TH1F *caloTrueRecPtTaggedCaloMCIsoCell = new TH1F(Form("fCaloTrueRecPtTaggedCaloMCIsoCell_%i_%i",r,e), Form("calo photons with Cell ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtTaggedCaloMCIsoCell[r][e] = (TH1F*) caloTrueRecPtTaggedCaloMCIsoCell->Clone(Form("fCaloTrueRecPtTaggedCaloMCIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloTrueRecPtTaggedCaloMCIsoCell[r][e]->Sumw2();
          fCaloFolderTrue->Add(fCaloTrueRecPtTaggedCaloMCIsoCell[r][e]);
        }

      }
    }

    // Inv mass histos
    fCaloInvMass = new TH2F("fCaloInvMass","fCaloInvMass;m_{#gamma #gamma} (GeV/c^{2});conv p_{T} (GeV/c)",nMassBins,minMass,maxMass,nPtBins,minPt,maxPt);
    fCaloFolderRec->Add(fCaloInvMass);
    if(fIsMC>0){
       fCaloTrueInvMass = new TH2F("fCaloTrueInvMass","fCaloTrueInvMass;m_{#gamma #gamma} (GeV/c^{2});conv p_{T} (GeV/c)",nMassBins,minMass,maxMass,nPtBins,minPt,maxPt);
       fCaloTrueInvMass_FromDecay = new TH2F("fCaloTrueInvMass_FromDecay","fCaloTrueInvMass_FromDecay;m_{#gamma #gamma} (GeV/c^{2});conv p_{T} (GeV/c)",nMassBins,minMass,maxMass,nPtBins,minPt,maxPt);
       fCaloTrueInvMass_FromDirect = new TH2F("fCaloTrueInvMass_FromDirect","fCaloTrueInvMass_FromDirect;m_{#gamma #gamma} (GeV/c^{2});conv p_{T} (GeV/c)",nMassBins,minMass,maxMass,nPtBins,minPt,maxPt);
       fCaloTrueInvMass->Sumw2();
       fCaloTrueInvMass_FromDecay->Sumw2();
       fCaloTrueInvMass_FromDirect->Sumw2();
       fCaloFolderTrue->Add(fCaloTrueInvMass);
       fCaloFolderTrue->Add(fCaloTrueInvMass_FromDecay);
       fCaloFolderTrue->Add(fCaloTrueInvMass_FromDirect);
    }
    for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
    {
        TH2F *caloInvMassAntiChargedIsolated = new TH2F(Form("caloInvMassAntiChargedIsolated_R%1.1f",fTrackIsolationR.at(r)),Form("fCaloInvMassAntiChargedIsolated_R%1.1f;m_{#gamma #gamma} (GeV/c^{2});conv p_{T} (GeV/c)",fTrackIsolationR.at(r)),nMassBins,minMass,maxMass,nPtBins,minPt,maxPt);
        fCaloInvMassAntiChargedIsolated[r] = (TH2F*) caloInvMassAntiChargedIsolated->Clone(Form("fCaloInvMassAntiChargedIsolated_R%1.1f",fTrackIsolationR.at(r)));
        fCaloInvMassAntiChargedIsolated[r]->Sumw2();
        fCaloFolderRec->Add(fCaloInvMassAntiChargedIsolated[r]);

        for (UInt_t e = 0; e < fTrackIsolationE.size(); e++)
        {
            TH2F *caloInvMassChargedIsolated = new TH2F(Form("caloInvMassChargedIsolated_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)),Form("fCaloInvMassChargedIsolated_R%1.1f_E%1.1f;m_{#gamma #gamma} (GeV/c^{2});conv p_{T} (GeV/c)",fTrackIsolationR.at(r),fTrackIsolationE.at(e)),nMassBins,minMass,maxMass,nPtBins,minPt,maxPt);
            fCaloInvMassChargedIsolated[r][e] = (TH2F*) caloInvMassChargedIsolated->Clone(Form("fCaloInvMassChargedIsolated_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
            fCaloInvMassChargedIsolated[r][e]->Sumw2();
            fCaloFolderRec->Add(fCaloInvMassChargedIsolated[r][e]);
        }
    }
    for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
    {
        TH2F *caloInvMassAntiNeutralIsolated = new TH2F(Form("caloInvMassAntiNeutralIsolated_R%1.1f",fNeutralIsolationR.at(r)),Form("fCaloInvMassAntiNeutralIsolated_R%1.1f;m_{#gamma #gamma} (GeV/c^{2});conv p_{T} (GeV/c)",fNeutralIsolationR.at(r)),nMassBins,minMass,maxMass,nPtBins,minPt,maxPt);
        fCaloInvMassAntiNeutralIsolated[r] = (TH2F*) caloInvMassAntiNeutralIsolated->Clone(Form("fCaloInvMassAntiNeutralIsolated_R%1.1f",fNeutralIsolationR.at(r)));
        fCaloInvMassAntiNeutralIsolated[r]->Sumw2();
        fCaloFolderRec->Add(fCaloInvMassAntiNeutralIsolated[r]);

        TH2F *caloInvMassAntiCellIsolated = new TH2F(Form("caloInvMassAntiCellIsolated_R%1.1f",fNeutralIsolationR.at(r)),Form("fCaloInvMassAntiCellIsolated_R%1.1f;m_{#gamma #gamma} (GeV/c^{2});conv p_{T} (GeV/c)",fNeutralIsolationR.at(r)),nMassBins,minMass,maxMass,nPtBins,minPt,maxPt);
        fCaloInvMassAntiCellIsolated[r] = (TH2F*) caloInvMassAntiCellIsolated->Clone(Form("fCaloInvMassAntiCellIsolated_R%1.1f",fNeutralIsolationR.at(r)));
        fCaloInvMassAntiCellIsolated[r]->Sumw2();
        fCaloFolderRec->Add(fCaloInvMassAntiCellIsolated[r]);

        TH2F *caloInvMassAntiFullIsolated = new TH2F(Form("caloInvMassAntiFullIsolated_R%1.1f",fNeutralIsolationR.at(r)),Form("fCaloInvMassAntiFullIsolated_R%1.1f;m_{#gamma #gamma} (GeV/c^{2});conv p_{T} (GeV/c)",fNeutralIsolationR.at(r)),nMassBins,minMass,maxMass,nPtBins,minPt,maxPt);
        fCaloInvMassAntiFullIsolated[r] = (TH2F*) caloInvMassAntiFullIsolated->Clone(Form("fCaloInvMassAntiFullIsolated_R%1.1f",fNeutralIsolationR.at(r)));
        fCaloInvMassAntiFullIsolated[r]->Sumw2();
        fCaloFolderRec->Add(fCaloInvMassAntiFullIsolated[r]);

        for (UInt_t e = 0; e < fNeutralIsolationE.size(); e++)
        {
            TH2F *caloInvMassNeutralIsolated = new TH2F(Form("caloInvMassNeutralIsolated_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)),Form("fCaloInvMassNeutralIsolated_R%1.1f_E%1.1f;m_{#gamma #gamma} (GeV/c^{2});conv p_{T} (GeV/c)",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)),nMassBins,minMass,maxMass,nPtBins,minPt,maxPt);
            fCaloInvMassNeutralIsolated[r][e] = (TH2F*) caloInvMassNeutralIsolated->Clone(Form("fCaloInvMassNeutralIsolated_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
            fCaloInvMassNeutralIsolated[r][e]->Sumw2();
            fCaloFolderRec->Add(fCaloInvMassNeutralIsolated[r][e]);

            TH2F *caloInvMassCellIsolated = new TH2F(Form("caloInvMassCellIsolated_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)),Form("fCaloInvMassCellIsolated_R%1.1f_E%1.1f;m_{#gamma #gamma} (GeV/c^{2});conv p_{T} (GeV/c)",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)),nMassBins,minMass,maxMass,nPtBins,minPt,maxPt);
            fCaloInvMassCellIsolated[r][e] = (TH2F*) caloInvMassCellIsolated->Clone(Form("fCaloInvMassCellIsolated_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
            fCaloInvMassCellIsolated[r][e]->Sumw2();
            fCaloFolderRec->Add(fCaloInvMassCellIsolated[r][e]);

            TH2F *caloInvMassFullIsolated = new TH2F(Form("caloInvMassFullIsolated_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)),Form("fCaloInvMassFullIsolated_R%1.1f_E%1.1f;m_{#gamma #gamma} (GeV/c^{2});conv p_{T} (GeV/c)",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)),nMassBins,minMass,maxMass,nPtBins,minPt,maxPt);
            fCaloInvMassFullIsolated[r][e] = (TH2F*) caloInvMassFullIsolated->Clone(Form("fCaloInvMassFullIsolated_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
            fCaloInvMassFullIsolated[r][e]->Sumw2();
            fCaloFolderRec->Add(fCaloInvMassFullIsolated[r][e]);
        }
    }


    // M02
    Int_t bins[3] = {nMassBins,nPtBins,nMassBins};
    Double_t xmin[3] = {minMass,minPt,minMass};
    Double_t xmax[3] = {maxMass,maxPt,maxMass};
    fCaloM02 = new THnSparseF("fCaloM02","fCaloM02;M02;calo p_{T} (GeV/c); sub clus mass (GeV/c^{2})",3,bins,xmin,xmax);
    fCaloFolderRec->Add(fCaloM02);
    if(fIsMC>0){
       fCaloTrueM02 = new THnSparseF("fCaloTrueM02","fCaloTrueM02;M02;calo p_{T} (GeV/c); sub clus mass (GeV/c^{2})",3,bins,xmin,xmax);
       fCaloTrueM02_FromDecay = new THnSparseF("fCaloTrueM02_FromDecay","fCaloTrueM02_FromDecay;M02;calo p_{T} (GeV/c); sub clus mass (GeV/c^{2})",3,bins,xmin,xmax);
       fCaloTrueM02_FromDirect = new THnSparseF("fCaloTrueM02_FromDirect","fCaloTrueM02_FromDirect;M02;calo p_{T} (GeV/c); sub clus mass (GeV/c^{2})",3,bins,xmin,xmax);
       fCaloTrueM02->Sumw2();
       fCaloTrueM02_FromDecay->Sumw2();
       fCaloTrueM02_FromDirect->Sumw2();
       fCaloFolderTrue->Add(fCaloTrueM02);
       fCaloFolderTrue->Add(fCaloTrueM02_FromDecay);
       fCaloFolderTrue->Add(fCaloTrueM02_FromDirect);
    }
    for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
    {
        THnSparseF *caloM02AntiChargedIsolated = new THnSparseF(Form("caloM02AntiChargedIsolated_R%1.1f",fTrackIsolationR.at(r)),Form("fCaloM02AntiChargedIsolated_R%1.1f;M02;calo p_{T} (GeV/c); sub clus mass (GeV/c^{2})",fTrackIsolationR.at(r)),3,bins,xmin,xmax);
        fCaloM02AntiChargedIsolated[r] = (THnSparseF*) caloM02AntiChargedIsolated->Clone(Form("fCaloM02AntiChargedIsolated_R%1.1f",fTrackIsolationR.at(r)));
        fCaloM02AntiChargedIsolated[r]->Sumw2();
        fCaloFolderRec->Add(fCaloM02AntiChargedIsolated[r]);

        // used for ABCD method
        THnSparseF *caloM02vsIsoVsPt = new THnSparseF(Form("caloM02vsIsoVsPt_R%1.1f",fTrackIsolationR.at(r)),Form("fCaloM02vsIsoVsPt_R%1.1f;M02;charged ISO ;calo p_{T} (GeV/c)",fTrackIsolationR.at(r)),3,bins,xmin,xmax);
        fCaloM02vsIsoVsPt[r] = (THnSparseF*) caloM02vsIsoVsPt->Clone(Form("fCaloM02vsIsoVsPt_R%1.1f",fTrackIsolationR.at(r)));
        fCaloM02vsIsoVsPt[r]->Sumw2();
        fCaloFolderRec->Add(fCaloM02vsIsoVsPt[r]);

        for (UInt_t e = 0; e < fTrackIsolationE.size(); e++)
        {
            THnSparseF *caloM02ChargedIsolated = new THnSparseF(Form("caloM02ChargedIsolated_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)),Form("fCaloM02ChargedIsolated_R%1.1f_E%1.1f;M02;calo p_{T} (GeV/c); sub clus mass (GeV/c^{2})",fTrackIsolationR.at(r),fTrackIsolationE.at(e)),3,bins,xmin,xmax);
            fCaloM02ChargedIsolated[r][e] = (THnSparseF*) caloM02ChargedIsolated->Clone(Form("fCaloM02ChargedIsolated_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
            fCaloM02ChargedIsolated[r][e]->Sumw2();
            fCaloFolderRec->Add(fCaloM02ChargedIsolated[r][e]);

            THnSparseF *caloM02ChargedIsolated_FromDirect = new THnSparseF(Form("caloM02ChargedIsolated_FromDirect_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)),Form("fCaloM02ChargedIsolated_FromDirect_R%1.1f_E%1.1f;M02;calo p_{T} (GeV/c); sub clus mass (GeV/c^{2})",fTrackIsolationR.at(r),fTrackIsolationE.at(e)),3,bins,xmin,xmax);
            fCaloM02ChargedIsolated_FromDirect[r][e] = (THnSparseF*) caloM02ChargedIsolated_FromDirect->Clone(Form("fCaloM02ChargedIsolated_FromDirect_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
            fCaloM02ChargedIsolated_FromDirect[r][e]->Sumw2();
            fCaloFolderRec->Add(fCaloM02ChargedIsolated_FromDirect[r][e]);

            THnSparseF *caloM02ChargedIsolated_FromDirectConvOnly = new THnSparseF(Form("caloM02ChargedIsolated_FromDirectConvOnly_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)),Form("fCaloM02ChargedIsolated_FromDirectConvOnly_R%1.1f_E%1.1f;M02;calo p_{T} (GeV/c); sub clus mass (GeV/c^{2})",fTrackIsolationR.at(r),fTrackIsolationE.at(e)),3,bins,xmin,xmax);
            fCaloM02ChargedIsolated_FromDirectConvOnly[r][e] = (THnSparseF*) caloM02ChargedIsolated_FromDirectConvOnly->Clone(Form("fCaloM02ChargedIsolated_FromDirectConvOnly_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
            fCaloM02ChargedIsolated_FromDirectConvOnly[r][e]->Sumw2();
            fCaloFolderRec->Add(fCaloM02ChargedIsolated_FromDirectConvOnly[r][e]);

            THnSparseF *caloM02ChargedIsolated_FromDirectNormOnly = new THnSparseF(Form("caloM02ChargedIsolated_FromDirectNormOnly_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)),Form("fCaloM02ChargedIsolated_FromDirectNormOnly_R%1.1f_E%1.1f;M02;calo p_{T} (GeV/c); sub clus mass (GeV/c^{2})",fTrackIsolationR.at(r),fTrackIsolationE.at(e)),3,bins,xmin,xmax);
            fCaloM02ChargedIsolated_FromDirectNormOnly[r][e] = (THnSparseF*) caloM02ChargedIsolated_FromDirectNormOnly->Clone(Form("fCaloM02ChargedIsolated_FromDirectNormOnly_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
            fCaloM02ChargedIsolated_FromDirectNormOnly[r][e]->Sumw2();
            fCaloFolderRec->Add(fCaloM02ChargedIsolated_FromDirectNormOnly[r][e]);
        }
    }
    for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
    {
        THnSparseF *caloM02AntiNeutralIsolated = new THnSparseF(Form("caloM02AntiNeutralIsolated_R%1.1f",fNeutralIsolationR.at(r)),Form("fCaloM02AntiNeutralIsolated_R%1.1f;M02;calo p_{T} (GeV/c); sub clus mass (GeV/c^{2})",fNeutralIsolationR.at(r)),3,bins,xmin,xmax);
        fCaloM02AntiNeutralIsolated[r] = (THnSparseF*) caloM02AntiNeutralIsolated->Clone(Form("fCaloM02AntiNeutralIsolated_R%1.1f",fNeutralIsolationR.at(r)));
        fCaloM02AntiNeutralIsolated[r]->Sumw2();
        fCaloFolderRec->Add(fCaloM02AntiNeutralIsolated[r]);

        THnSparseF *caloM02AntiCellIsolated = new THnSparseF(Form("caloM02AntiCellIsolated_R%1.1f",fNeutralIsolationR.at(r)),Form("fCaloM02AntiCellIsolated_R%1.1f;M02;calo p_{T} (GeV/c); sub clus mass (GeV/c^{2})",fNeutralIsolationR.at(r)),3,bins,xmin,xmax);
        fCaloM02AntiCellIsolated[r] = (THnSparseF*) caloM02AntiCellIsolated->Clone(Form("fCaloM02AntiCellIsolated_R%1.1f",fNeutralIsolationR.at(r)));
        fCaloM02AntiCellIsolated[r]->Sumw2();
        fCaloFolderRec->Add(fCaloM02AntiCellIsolated[r]);

        THnSparseF *caloM02AntiFullIsolated = new THnSparseF(Form("caloM02AntiFullIsolated_R%1.1f",fNeutralIsolationR.at(r)),Form("fCaloM02AntiFullIsolated_R%1.1f;M02;calo p_{T} (GeV/c); sub clus mass (GeV/c^{2})",fNeutralIsolationR.at(r)),3,bins,xmin,xmax);
        fCaloM02AntiFullIsolated[r] = (THnSparseF*) caloM02AntiFullIsolated->Clone(Form("fCaloM02AntiFullIsolated_R%1.1f",fNeutralIsolationR.at(r)));
        fCaloM02AntiFullIsolated[r]->Sumw2();
        fCaloFolderRec->Add(fCaloM02AntiFullIsolated[r]);


        for (UInt_t e = 0; e < fNeutralIsolationE.size(); e++)
        {
            THnSparseF *caloM02NeutralIsolated = new THnSparseF(Form("caloM02NeutralIsolated_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)),Form("fCaloM02NeutralIsolated_R%1.1f_E%1.1f;M02;calo p_{T} (GeV/c); sub clus mass (GeV/c^{2})",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)),3,bins,xmin,xmax);
            fCaloM02NeutralIsolated[r][e] = (THnSparseF*) caloM02NeutralIsolated->Clone(Form("fCaloM02NeutralIsolated_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
            fCaloM02NeutralIsolated[r][e]->Sumw2();
            fCaloFolderRec->Add(fCaloM02NeutralIsolated[r][e]);

            THnSparseF *caloM02CellIsolated = new THnSparseF(Form("caloM02CellIsolated_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)),Form("fCaloM02CellIsolated_R%1.1f_E%1.1f;M02;calo p_{T} (GeV/c); sub clus mass (GeV/c^{2})",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)),3,bins,xmin,xmax);
            fCaloM02CellIsolated[r][e] = (THnSparseF*) caloM02CellIsolated->Clone(Form("fCaloM02CellIsolated_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
            fCaloM02CellIsolated[r][e]->Sumw2();
            fCaloFolderRec->Add(fCaloM02CellIsolated[r][e]);

            THnSparseF *caloM02FullIsolated = new THnSparseF(Form("caloM02FullIsolated_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)),Form("fCaloM02FullIsolated_R%1.1f_E%1.1f;M02;calo p_{T} (GeV/c); sub clus mass (GeV/c^{2})",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)),3,bins,xmin,xmax);
            fCaloM02FullIsolated[r][e] = (THnSparseF*) caloM02FullIsolated->Clone(Form("fCaloM02FullIsolated_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
            fCaloM02FullIsolated[r][e]->Sumw2();
            fCaloFolderRec->Add(fCaloM02FullIsolated[r][e]);
        }
    }

    //
    // ─── GENERATOR LEVEL ─────────────────────────────────────────────
    //
    if(fIsMC > 0){
      fGeneratorFolder          = new TList();
      fGeneratorFolder->SetName("genLevel");
      fGeneratorFolder->SetOwner(kTRUE);
      fOutputList->Add(fGeneratorFolder);

      fHistoMCHeaders = new TH1I("MC_Headers", "MC_Headers", 20, 0, 20);
      fGenPhotonPt = new TH1F("fGenPhotonPt","fGenPhotonPt;gen. p_{T} (GeV/c); counts",nPtBins,minPt,maxPt);
      fGenPhotonPt_FromDecay  = new TH1F("fGenPhotonPt_FromDecay","fGenPhotonPt_FromDecay;gen. p_{T} (GeV/c); counts",nPtBins,minPt,maxPt);
      fGenPhotonPt_FromDirect  = new TH1F("fGenPhotonPt_FromDirect","fGenPhotonPt_FromDirect;gen. p_{T} (GeV/c); counts",nPtBins,minPt,maxPt);
      fGenPhotonPtInEMCalAcc  = new TH1F("fGenPhotonPtInEMCalAcc","fGenPhotonPtInEMCalAcc;gen. p_{T} (GeV/c); counts",nPtBins,minPt,maxPt);
      fGenPhotonPtInEMCalAcc_FromDecay  = new TH1F("fGenPhotonPtInEMCalAcc_FromDecay","fGenPhotonPtInEMCalAcc_FromDecay;gen. p_{T} (GeV/c); counts",nPtBins,minPt,maxPt);
      fGenPhotonPtInEMCalAcc_FromDirect  = new TH1F("fGenPhotonPtInEMCalAcc_FromDirect","fGenPhotonPtInEMCalAcc_FromDirect;gen. p_{T} (GeV/c); counts",nPtBins,minPt,maxPt);
      for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
      {
        for (UInt_t e = 0; e < fTrackIsolationE.size(); e++)
        {
          TH1F *genPhotonPtInEMCalAccChargedMCIso_FromDirect = new TH1F(Form("genPhotonPtInEMCalAccChargedMCIso_FromDirect_%i_%i",r,e), Form("direct gen photon with MC ISO < %1.1f GeV in R < %1.1f;gen. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fGenPhotonPtInEMCalAccChargedMCIso_FromDirect[r][e] = (TH1F*) genPhotonPtInEMCalAccChargedMCIso_FromDirect->Clone(Form("fGenPhotonPtInEMCalAccChargedMCIso_FromDirect_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fGenPhotonPtInEMCalAccChargedMCIso_FromDirect[r][e]->Sumw2();
          fGeneratorFolder->Add(fGenPhotonPtInEMCalAccChargedMCIso_FromDirect[r][e]);
        }   
        TH2F *genPhotonChargedMCIsoInEMCalAcc_FromDirect = new TH2F(Form("genPhotonChargedMCIsoInEMCalAcc_FromDirect_%i",r), Form("direct gen photon with MC ISO in R < %1.1f;ISO (GeV); p_{T} (GeV/c)",fTrackIsolationR.at(r)), nPtBins,minPt,maxPt,nPtBins,minPt,maxPt);
        fGenPhotonChargedMCIsoInEMCalAcc_FromDirect[r] = (TH2F*) genPhotonChargedMCIsoInEMCalAcc_FromDirect->Clone(Form("fGenPhotonChargedMCIsoInEMCalAcc_FromDirect_R%1.1f",fTrackIsolationR.at(r)));
        fGenPhotonChargedMCIsoInEMCalAcc_FromDirect[r]->Sumw2();
        fGeneratorFolder->Add(fGenPhotonChargedMCIsoInEMCalAcc_FromDirect[r]);
      }
      fGenPhotonPtFoundNormCluster  = new TH1F("fGenPhotonPtFoundNormCluster","fGenPhotonPtFoundNormCluster;gen. p_{T} (GeV/c); counts",nPtBins,minPt,maxPt);
      fGenPhotonPtFoundTaggingCluster  = new TH1F("fGenPhotonPtFoundTaggingCluster","fGenPhotonPtFoundTaggingCluster;gen. p_{T} (GeV/c); counts",nPtBins,minPt,maxPt);
      fGenPhotonPtFoundIsoCluster  = new TH1F("fGenPhotonPtFoundIsoCluster","fGenPhotonPtFoundIsoCluster;gen. p_{T} (GeV/c); counts",nPtBins,minPt,maxPt);
      fGenPhotonEFoundNoClusterVsCellE  = new TH2F("fGenPhotonEFoundNoClusterVsCellE","fGenPhotonEFoundNoClusterVsCellE;gen. E (GeV/c), cell E; counts",nPtBins,minPt,maxPt,nPtBins,minPt,maxPt);


      fGenPi0Pt  = new TH1F("fGenPi0Pt","fGenPi0Pt;gen. p_{T} (GeV/c); counts",nPtBins,minPt,maxPt);
      fGenPi0PtInEMCalAcc  = new TH1F("fGenPi0PtInEMCalAcc","fGenPi0PtInEMCalAcc;gen. p_{T} (GeV/c); counts",nPtBins,minPt,maxPt);
      fGenPi0PtInEMCalAcc_BothGammaInEMCal  = new TH1F("fGenPi0PtInEMCalAcc_BothGammaInEMCal","fGenPi0PtInEMCalAcc_BothGammaInEMCal;gen. p_{T} (GeV/c); counts",nPtBins,minPt,maxPt);
      fGenPi0PtInEMCalAcc_BothGammaInClusters  = new TH1F("fGenPi0PtInEMCalAcc_BothGammaInClusters","fGenPi0PtInEMCalAcc_BothGammaInClusters;gen. p_{T} (GeV/c); counts",nPtBins,minPt,maxPt);

      fHistoMCHeaders->Sumw2();
      fGenPhotonPt->Sumw2();
      fGenPhotonPt_FromDecay->Sumw2();
      fGenPhotonPt_FromDirect->Sumw2();
      fGenPhotonPtInEMCalAcc->Sumw2();
      fGenPhotonPtInEMCalAcc_FromDecay->Sumw2();
      fGenPhotonPtInEMCalAcc_FromDirect->Sumw2();
      fGenPhotonPtFoundNormCluster->Sumw2();
      fGenPhotonPtFoundTaggingCluster->Sumw2();
      fGenPhotonPtFoundIsoCluster->Sumw2();
      fGenPhotonEFoundNoClusterVsCellE->Sumw2();
      fGenPi0Pt->Sumw2();
      fGenPi0PtInEMCalAcc->Sumw2();
      fGenPi0PtInEMCalAcc_BothGammaInEMCal->Sumw2();
      fGenPi0PtInEMCalAcc_BothGammaInClusters->Sumw2();

      fGeneratorFolder->Add(fHistoMCHeaders);
      fGeneratorFolder->Add(fGenPhotonPt);
      fGeneratorFolder->Add(fGenPhotonPt_FromDecay);
      fGeneratorFolder->Add(fGenPhotonPt_FromDirect);
      fGeneratorFolder->Add(fGenPhotonPtInEMCalAcc);
      fGeneratorFolder->Add(fGenPhotonPtInEMCalAcc_FromDecay);
      fGeneratorFolder->Add(fGenPhotonPtInEMCalAcc_FromDirect);
      fGeneratorFolder->Add(fGenPhotonPtFoundNormCluster);
      fGeneratorFolder->Add(fGenPhotonPtFoundTaggingCluster);
      fGeneratorFolder->Add(fGenPhotonPtFoundIsoCluster);
      fGeneratorFolder->Add(fGenPhotonEFoundNoClusterVsCellE);
      fGeneratorFolder->Add(fGenPi0Pt);
      fGeneratorFolder->Add(fGenPi0PtInEMCalAcc);
      fGeneratorFolder->Add(fGenPi0PtInEMCalAcc_BothGammaInEMCal);
      fGeneratorFolder->Add(fGenPi0PtInEMCalAcc_BothGammaInClusters);
    }

  }

  PostData(1, fOutputList);

  OpenFile(2);
  // fConversionCandidates = new TClonesArray("AliAODConversionPhoton",50);
  fClusterEMCalCandidates = new TList();
  fClusterEMCalCandidatesIsolation = new TList();
  fClusterEMCalCandidatesTagging = new TList();
  fClusterPHOSCandidates = new TList();
  fMCParticles = new TList();

  TString treename = "CaloTree";
  if(fCorrTaskSetting.CompareTo("")){
      treename = Form("CaloTree_%s",fCorrTaskSetting.Data());
  }
  fAnalysisTree = new TTree(treename,treename);
  if(!fUseHistograms){  // full tree

    fAnalysisTree->Branch("Event_Rho", &fBuffer_EventRho,"Event_Rho/F");
    fAnalysisTree->Branch("Event_RhoMC", &fBuffer_EventRhoMC,"Event_RhoMC/F");
    if(fIsMC>0) fAnalysisTree->Branch("Event_Weight", &fBuffer_EventWeight,"Event_Weight/D");
    if(fIsMC>0) fAnalysisTree->Branch("Event_Xsection", &fBuffer_EventXsection,"Event_Xsection/F");
    if(fIsMC>0) fAnalysisTree->Branch("Event_Ntrials", &fBuffer_EventNtrials,"Event_Ntrials/s");
    fAnalysisTree->Branch("Event_IsTriggered", &fBuffer_EventIsTriggered,"Event_IsTriggered/O");
    fAnalysisTree->Branch("Cluster_E","std::vector<Float_t>",&fBuffer_ClusterE);
    fAnalysisTree->Branch("Cluster_Px","std::vector<Float_t>",&fBuffer_ClusterPx);
    fAnalysisTree->Branch("Cluster_Py","std::vector<Float_t>",&fBuffer_ClusterPy);
    fAnalysisTree->Branch("Cluster_Pz","std::vector<Float_t>",&fBuffer_ClusterPz);
    fAnalysisTree->Branch("Cluster_M02","std::vector<Float_t>",&fBuffer_ClusterM02);
    fAnalysisTree->Branch("Cluster_M02Recalc","std::vector<Float_t>",&fBuffer_ClusterM02Recalc);
    fAnalysisTree->Branch("Cluster_M20","std::vector<Float_t>",&fBuffer_ClusterM20);
    fAnalysisTree->Branch("Cluster_V1SplitMass","std::vector<Float_t>",&fBuffer_ClusterV1SplitMass);
    fAnalysisTree->Branch("Cluster_NLM","std::vector<UShort_t>",&fBuffer_ClusterNLM);
    fAnalysisTree->Branch("Cluster_SM","std::vector<UShort_t>",&fBuffer_ClusterSM);
    fAnalysisTree->Branch("Cluster_EFrac","std::vector<Float_t>",&fBuffer_ClusterEFrac);
    fAnalysisTree->Branch("Cluster_IsoCharged1","std::vector<Float_t>",&fBuffer_ClusterIsoCharged1);
    fAnalysisTree->Branch("Cluster_IsoCharged2","std::vector<Float_t>",&fBuffer_ClusterIsoCharged2);
    fAnalysisTree->Branch("Cluster_IsoCharged3","std::vector<Float_t>",&fBuffer_ClusterIsoCharged3);
    fAnalysisTree->Branch("Cluster_IsoBckLeft","std::vector<Float_t>",&fBuffer_ClusterIsoBckLeft);
    fAnalysisTree->Branch("Cluster_MatchTrackdEta","std::vector<Float_t>",&fBuffer_ClusterMatchTrackdEta);
    fAnalysisTree->Branch("Cluster_MatchTrackdPhi","std::vector<Float_t>",&fBuffer_ClusterMatchTrackdPhi);
    fAnalysisTree->Branch("Cluster_MatchTrackP","std::vector<Float_t>",&fBuffer_ClusterMatchTrackP);
    fAnalysisTree->Branch("Cluster_MatchTrackPt","std::vector<Float_t>",&fBuffer_ClusterMatchTrackPt);
    fAnalysisTree->Branch("Cluster_MatchTrackIsConv","std::vector<Bool_t>",&fBuffer_ClusterMatchTrackIsConv);
    if(fIsMC>0){
      fAnalysisTree->Branch("TrueCluster_E","std::vector<Float_t>",&fBuffer_TrueClusterE);
      fAnalysisTree->Branch("TrueCluster_Px","std::vector<Float_t>",&fBuffer_TrueClusterPx);
      fAnalysisTree->Branch("TrueCluster_Py","std::vector<Float_t>",&fBuffer_TrueClusterPy);
      fAnalysisTree->Branch("TrueCluster_Pz","std::vector<Float_t>",&fBuffer_TrueClusterPz);
      fAnalysisTree->Branch("TrueCluster_LeadingEFrac","std::vector<Float_t>",&fBuffer_TrueClusterLeadingEFrac);
      fAnalysisTree->Branch("TrueCluster_MCIsoCharged1","std::vector<Float_t>",&fBuffer_TrueClusterMCIsoCharged1);
      fAnalysisTree->Branch("TrueCluster_MCIsoCharged2","std::vector<Float_t>",&fBuffer_TrueClusterMCIsoCharged2);
      fAnalysisTree->Branch("TrueCluster_MCIsoCharged3","std::vector<Float_t>",&fBuffer_TrueClusterMCIsoCharged3);
      fAnalysisTree->Branch("TrueCluster_MCIsoBckLeft","std::vector<Float_t>",&fBuffer_TrueClusterMCIsoBckLeft);
      fAnalysisTree->Branch("TrueCluster_MCTag","std::vector<Int_t>",&fBuffer_TrueClusterMCTag);
      fAnalysisTree->Branch("TrueCluster_IsConv","std::vector<Bool_t>",&fBuffer_TrueClusterIsConv);
    
      fAnalysisTree->Branch("GenPhoton_E","std::vector<Float_t>",&fBuffer_GenPhotonE);
      fAnalysisTree->Branch("GenPhoton_Px","std::vector<Float_t>",&fBuffer_GenPhotonPx);
      fAnalysisTree->Branch("GenPhoton_Py","std::vector<Float_t>",&fBuffer_GenPhotonPy);
      fAnalysisTree->Branch("GenPhoton_Pz","std::vector<Float_t>",&fBuffer_GenPhotonPz);
      fAnalysisTree->Branch("GenPhoton_MCIsoCharged1","std::vector<Float_t>",&fBuffer_GenPhotonMCIsoCharged1);
      fAnalysisTree->Branch("GenPhoton_MCIsoCharged2","std::vector<Float_t>",&fBuffer_GenPhotonMCIsoCharged2);
      fAnalysisTree->Branch("GenPhoton_MCIsoCharged3","std::vector<Float_t>",&fBuffer_GenPhotonMCIsoCharged3);
      fAnalysisTree->Branch("GenPhoton_MCIsoBckLeft","std::vector<Float_t>",&fBuffer_GenPhotonMCIsoBckLeft);
      fAnalysisTree->Branch("GenPhoton_IsConv","std::vector<Bool_t>",&fBuffer_GenPhotonIsConv);  
      fAnalysisTree->Branch("GenPhoton_MCTag","std::vector<Int_t>",&fBuffer_GenPhotonMCTag);  
    }
  }

  PostData(2, fAnalysisTree);

}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskGammaIsoTree::Notify()
{
    return kTRUE;
}
//________________________________________________________________________
void AliAnalysisTaskGammaIsoTree::UserExec(Option_t *){
  // auto start = std::chrono::high_resolution_clock::now();
  fInputEvent                         = InputEvent();
  ((AliCaloPhotonCuts*)fClusterCutsEMC)->InitializeEMCAL(fInputEvent);
  if(fIsMC>0) fMCEvent                  = MCEvent();
  if((fIsMC) > 0 && (!fMCEvent)) {printf("Error: No MC Event");return;}
  // Get V0 reader
  fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
  if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader

  if(fIsMC > 0 && fInputEvent->IsA()==AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled())){
     RelabelAODPhotonCandidates(kTRUE);    // In case of AODMC relabeling MC
     fV0Reader->RelabelAODs(kTRUE);
  }

  Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();
  if(InputEvent()->IsIncompleteDAQ()==kTRUE) eventQuality = 2;  // incomplete event
  if(eventQuality == 2 || eventQuality == 3){// Event Not Accepted due to MC event missing or wrong trigger for V0ReaderV1 or because it is incomplete
    fHistoNEvents->Fill(eventQuality);
    if (fIsMC>1) fHistoNEventsWOWeight->Fill(eventQuality);
    return;
  }

  fReaderGammas = fV0Reader->GetReconstructedGammas(); // Gammas from default Cut

  Bool_t isRunningEMCALrelAna = kFALSE;
  if (fClusterCutsEMC->GetClusterType() == 1) isRunningEMCALrelAna = kTRUE;

  Int_t eventNotAccepted              = fEventCuts->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon,isRunningEMCALrelAna);
  // if(eventNotAccepted) return; // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1


  AliRhoParameter* outrho= (AliRhoParameter*) InputEvent()->FindListObject(fRhoOutName.Data());
  if(!outrho) AliInfo("could not find rho container!");

  AliRhoParameter* outrhoMC = NULL;
  if(fIsMC>0){
    outrhoMC= (AliRhoParameter*) InputEvent()->FindListObject(fRhoOutNameMC.Data());
    if(!outrhoMC) AliInfo("could not find rho MC container!");
  }

  Float_t xsection = 0;
  Float_t ntrials = 0;
  if (fIsMC > 0){
      fWeightJetJetMC       = 1;
      Float_t maxjetpt      = -1.;
      Float_t pthard = -1;
      
      if(fEventCuts->GetUseJetFinderForOutliers()) maxjetpt = fOutlierJetReader->GetMaxJetPt();
      Bool_t isMCJet        = ((AliConvEventCuts*)fEventCuts)->IsJetJetMCEventAccepted( fMCEvent, fWeightJetJetMC ,pthard, fInputEvent, maxjetpt);
      fEventCuts->GetXSectionAndNTrials(fMCEvent,xsection,ntrials,fInputEvent);
      
      if (fIsMC == 3){
        Double_t weightMult   = ((AliConvEventCuts*)fEventCuts)->GetWeightForMultiplicity(fV0Reader->GetNumberOfPrimaryTracks());
        fWeightJetJetMC       = fWeightJetJetMC*weightMult;
      }

      if (!isMCJet){
        fHistoNEvents->Fill(10,fWeightJetJetMC);
        if (fIsMC>1) fHistoNEventsWOWeight->Fill(10);
        return;
      }

  }

  Bool_t triggered = kTRUE;
  if(eventNotAccepted!=0){
      fHistoNEvents->Fill(eventNotAccepted,fWeightJetJetMC); // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1
      if (fIsMC>1) fHistoNEventsWOWeight->Fill(eventNotAccepted);
      if (eventNotAccepted==3 && fIsMC > 0){
        triggered = kFALSE;
      }else {
        return;
      }
  }

  if(eventQuality != 0 && triggered== kTRUE){// Event Not Accepted
    fHistoNEvents->Fill(eventQuality, fWeightJetJetMC);
    if (fIsMC>1) fHistoNEventsWOWeight->Fill(eventQuality); // Should be 0 here
    return;
  }

  if (triggered == kTRUE) {
    fHistoNEvents->Fill(eventQuality,fWeightJetJetMC);
    if (fIsMC>1) fHistoNEventsWOWeight->Fill(eventQuality); // Should be 0 here
  }

  fGeomEMCAL                          = AliEMCALGeometry::GetInstance();
  if(!fGeomEMCAL){ AliFatal("EMCal geometry not initialized!");}

  if(fIsMC> 0){
  // Process MC Particle
  if(((AliConvEventCuts*)fEventCuts)->GetSignalRejection() != 0){
    if(fInputEvent->IsA()==AliESDEvent::Class()){
    ((AliConvEventCuts*)fEventCuts)->GetNotRejectedParticles(((AliConvEventCuts*)fEventCuts)->GetSignalRejection(),
                                        ((AliConvEventCuts*)fEventCuts)->GetAcceptedHeader(),
                                        fMCEvent);
    }
    else if(fInputEvent->IsA()==AliAODEvent::Class()){
    ((AliConvEventCuts*)fEventCuts)->GetNotRejectedParticles(((AliConvEventCuts*)fEventCuts)->GetSignalRejection(),
                                      ((AliConvEventCuts*)fEventCuts)->GetAcceptedHeader(),
                                      fInputEvent);
    }
    if(((AliConvEventCuts*)fEventCuts)->GetAcceptedHeader()){
      for(Int_t i = 0;i<(((AliConvEventCuts*)fEventCuts)->GetAcceptedHeader())->GetEntries();i++){
        TString nameBin= fHistoMCHeaders->GetXaxis()->GetBinLabel(i+1);
        if (nameBin.CompareTo("")== 0){
          TString nameHeader = ((TObjString*)((TList*)((AliConvEventCuts*)fEventCuts)
                            ->GetAcceptedHeader())->At(i))->GetString();
          fHistoMCHeaders->GetXaxis()->SetBinLabel(i+1,nameHeader.Data());
        }
      }
    }
  }
}


  //
  // ─── MAIN PROCESSING ────────────────────────────────────────────────────────────
  //
  if(outrho){ // ugly workaround until problem is fixed
    fChargedRho = outrho->GetVal();
  } else{
    fChargedRho = 0.;
  }

  fChargedRhoMC = 0;
  if(outrhoMC){ // ugly workaround until problem is fixed
    fChargedRhoMC = outrhoMC->GetVal();
  } else{
    fChargedRhoMC = 0.;
  }
  for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
  {
    fChargedRhoTimesArea[r] = fChargedRho * TMath::Pi() * pow(fTrackIsolationR.at(r),2);
  }
  // auto startMCPart = std::chrono::high_resolution_clock::now();
  if(fIsMC>0) ProcessMCParticles();
  // auto endMCPart = std::chrono::high_resolution_clock::now();
  if(!fUseHistograms){
    fBuffer_EventWeight = fWeightJetJetMC;
    fBuffer_EventXsection = xsection;
    fBuffer_EventNtrials = ntrials;
    fBuffer_EventIsTriggered = kFALSE;
  }
  if (triggered==kFALSE){
    // fill output from MC
    if(!fUseHistograms){
      fBuffer_EventRho = fChargedRho;
      fBuffer_EventRhoMC = fChargedRhoMC;
      fAnalysisTree->Fill();
      PostData(2, fAnalysisTree);
    }
    ResetBuffer();
    return;
  }

  if(!fUseHistograms) fBuffer_EventIsTriggered = kTRUE;
  // vertex
  Double_t vertex[3] = {0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);


  //ProcessTracks(); // always run ProcessTracks before calo photons! (even if save tracks is false)
  if(fSaveConversions) ProcessConversionPhotons();
  // auto startCalo = std::chrono::high_resolution_clock::now();
  ProcessCaloPhotons(); // track matching is done here as well
  // auto endCalo = std::chrono::high_resolution_clock::now();
  // ReduceTrackInfo(); // track matching is done, we can remove cov matrix etc now

  // do not fill these in tree if user wants it
  // processing was needed anyways because of track matching
  // and isolation

  if( fIsMC > 0 && fInputEvent->IsA()==AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled())){
    RelabelAODPhotonCandidates(kFALSE); // Back to ESDMC Label
    fV0Reader->RelabelAODs(kFALSE);
  }
  // fill output
  if(!fUseHistograms){
    fBuffer_EventRho = fChargedRho;
    fBuffer_EventRhoMC = fChargedRhoMC;
    fAnalysisTree->Fill();
    PostData(2, fAnalysisTree);
  }
  // std::chrono::duration<double> eventLoading = startMCPart - start;
  // std::chrono::duration<double> mcLoop = endMCPart - startMCPart;
  // std::chrono::duration<double> caloLoop = endCalo - startCalo;
  // cout << "Event: " << eventLoading.count() << endl;
  // cout << "MC: " << mcLoop.count() << endl;
  // cout << "Calo: " << caloLoop.count() << endl;
  ResetBuffer();
  //gObjectTable->Print();
  PostData(1,fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskGammaIsoTree::Terminate(Option_t *){

}

//________________________________________________________________________
void AliAnalysisTaskGammaIsoTree::ResetBuffer(){
  // for vectors of pointers, memory needs to be freed manually

  // fConversionCandidates->Delete();
  // fClusterEMCalCandidates->Delete();
  // fClusterEMCalCandidatesIsolation->Delete();
  // fClusterEMCalCandidatesTagging->Delete();
  // fClusterPHOSCandidates->Delete();
  // fTracks->Delete();
  // fMCParticles->Delete();
  // fExtraClusterInfo->Delete();
  // fExtraClusterInfoBackground->Delete();


  // cout << "Entry = " << fAnalysisTree->GetEntriesFast() << endl;
  // cout << "Autosave = " << fAnalysisTree->GetAutoSave() << endl;
  // cout << "Autoflush = " << fAnalysisTree->GetAutoFlush() << endl;

  // fConversionCandidates->Delete();
  fClusterEMCalCandidates->Clear();
  fClusterEMCalCandidatesIsolation->Clear();
  fClusterEMCalCandidatesTagging->Clear();
  fClusterPHOSCandidates->Clear();
  fMCParticles->Clear();
  fBuffer_ClusterE.clear();     //!<! array buffer
  fBuffer_ClusterPx.clear(); 
  fBuffer_ClusterPy.clear(); 
  fBuffer_ClusterPz.clear(); 
  fBuffer_ClusterM02.clear(); 
  fBuffer_ClusterM02Recalc.clear(); 
  fBuffer_ClusterM20.clear(); 
  fBuffer_ClusterV1SplitMass.clear(); 
  fBuffer_ClusterNLM.clear(); 
  fBuffer_ClusterSM.clear(); 
  fBuffer_ClusterEFrac.clear(); 
  fBuffer_ClusterIsoCharged1.clear(); 
  fBuffer_ClusterIsoCharged2.clear(); 
  fBuffer_ClusterIsoCharged3.clear(); 
  fBuffer_ClusterIsoBckLeft.clear(); 
  fBuffer_ClusterMatchTrackdEta.clear(); 
  fBuffer_ClusterMatchTrackdPhi.clear(); 
  fBuffer_ClusterMatchTrackP.clear(); 
  fBuffer_ClusterMatchTrackPt.clear(); 
  fBuffer_ClusterMatchTrackIsConv.clear(); 
  fBuffer_TrueClusterE.clear(); 
  fBuffer_TrueClusterPx.clear(); 
  fBuffer_TrueClusterPy.clear(); 
  fBuffer_TrueClusterPz.clear(); 
  fBuffer_TrueClusterLeadingEFrac.clear();
  fBuffer_TrueClusterMCIsoCharged1.clear(); 
  fBuffer_TrueClusterMCIsoCharged2.clear(); 
  fBuffer_TrueClusterMCIsoCharged3.clear(); 
  fBuffer_TrueClusterMCIsoBckLeft.clear(); 
  fBuffer_TrueClusterMCTag.clear();
  fBuffer_TrueClusterIsConv.clear();
    
  fBuffer_GenPhotonE.clear();
  fBuffer_GenPhotonPx.clear();
  fBuffer_GenPhotonPy.clear();
  fBuffer_GenPhotonPz.clear();
  fBuffer_GenPhotonMCIsoCharged1.clear();
  fBuffer_GenPhotonMCIsoCharged2.clear();
  fBuffer_GenPhotonMCIsoCharged3.clear();
  fBuffer_GenPhotonMCIsoBckLeft.clear();
  fBuffer_GenPhotonIsConv.clear();
  fBuffer_GenPhotonMCTag.clear();


}
//________________________________________________________________________
void AliAnalysisTaskGammaIsoTree::ProcessConversionPhotons(){
   fReaderGammas    = fV0Reader->GetReconstructedGammas();
   Int_t pos = 0;
   for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){
    AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(i);
    if(!PhotonCandidate) continue;
    fIsFromDesiredHeader = kTRUE;

    if(fIsMC>0 && (fEventCuts->GetSignalRejection() != 0)){
      Int_t isPosFromMBHeader = fEventCuts->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCEvent, fInputEvent);
      if(isPosFromMBHeader == 0 && (fEventCuts->GetSignalRejection() != 3)) continue;
      Int_t isNegFromMBHeader = fEventCuts->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCEvent, fInputEvent);
      if(isNegFromMBHeader == 0 && (fEventCuts->GetSignalRejection() != 3)) continue;
      if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromDesiredHeader = kFALSE;
    }


    if(!fConvCuts->PhotonIsSelected(PhotonCandidate,fInputEvent)) continue;
    Bool_t isWithinTPC =IsWithinRadiusTPC(PhotonCandidate->Eta(),PhotonCandidate->Phi(),fExclusionRadius);
    if(!isWithinTPC) continue;
    if(fIsFromDesiredHeader){
      // new((*fConversionCandidates)[pos]) AliAODConversionPhoton(*PhotonCandidate);


      Int_t      tmp_tag= 0;
      Double32_t tmp_isoNeutral[2] = {0,0};
      Double32_t tmp_isoCell[2] = {0,0};
      Double32_t tmp_isoCharged[2] = {0,0};
      vector<Double32_t> isoCharged;
      vector<Double32_t> isoNeutral;
      vector<Double32_t> isoCell;

      if(fDoTrackIsolation) isoCharged = ProcessChargedIsolation(PhotonCandidate);
      if(fDoNeutralIsolation){
        isoNeutral = ProcessNeutralIsolation(PhotonCandidate);
      } else{
          for (UInt_t r = 0; r < isoNeutral.size(); r++)
          {
             isoNeutral.push_back(0.);
          }
      }
      if(fDoCellIsolation) {
        isoCell =  ProcessCellIsolation(PhotonCandidate);
      } else{
          for (UInt_t r = 0; r < isoNeutral.size(); r++)
          {
             isoCell.push_back(0.);
          }
      }
      if(fDoTagging) tmp_tag = ProcessTagging(PhotonCandidate);

      // when writing to tree, only fill first two R
      for (UInt_t r = 0; r < isoCharged.size(); r++)
      {
        if(r<2) tmp_isoCharged[r] = isoCharged.at(r);
      }
      for (UInt_t r = 0; r < isoNeutral.size(); r++)
      {
        if(r<2){
          tmp_isoNeutral[r] = isoCharged.at(r);
          tmp_isoCell[r] = isoCell.at(r);
        }
      }

      if(fUseHistograms) FillConversionHistos(PhotonCandidate,isoCharged,isoNeutral,isoCell,tmp_tag);

      if((fIsMC>0) && fUseHistograms) ProcessMCConversionPhoton(PhotonCandidate,isoCharged,isoNeutral,isoCell,tmp_tag);

      pos++;
    }

   }
}
//________________________________________________________________________
void AliAnalysisTaskGammaIsoTree::ProcessMCConversionPhoton(AliAODConversionPhoton* photon,vector<Double32_t> isoCharged,vector<Double32_t> isoNeutral,vector<Double32_t> isoCell,Int_t tmptag){
  if(!photon) return;
  if (!IsInEMCalAcceptance(photon)) return;
  Bool_t isTrueConv = IsTrueConversionPhoton(photon);
  if(!isTrueConv) return;
  if(!fIsFromDesiredHeader) return;
  Bool_t isDecayPhoton = IsDecayPhoton(photon);
  Bool_t isPromptPhoton = IsPromptPhoton(photon);
  Bool_t isFragPhoton = IsFragPhoton(photon);
  Bool_t isSignal = isPromptPhoton || isFragPhoton;

  if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));

  Int_t pconvlabel = GetConvPhotonMCLabel(photon);
  AliAODMCParticle* convPhotonMC = (AliAODMCParticle *) fAODMCTrackArray->At(pconvlabel);

  isoValues mcIso;
  vector<Double32_t> mcIsoCharged;
  vector<Double32_t> mcIsoNeutral;
  vector<Double32_t> mcIsoFull;
  mcIso = ProcessMCIsolation(pconvlabel);

  for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
  {
    Double_t isoCh = mcIso.isolationCone.at(r) - fChargedRhoTimesArea[r];
    mcIsoCharged.push_back(isoCh);
  }
  for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
  {
    mcIsoNeutral.push_back(mcIso.isolationCone.at(r+fTrackIsolationR.size()));
  }
  for (UInt_t r = 0; r < mcIsoCharged.size(); r++)
  {
    mcIsoFull.push_back(mcIsoCharged.at(r)+ mcIsoNeutral.at(r)); // they should always have same length
  }



  fConvTrueRecPt->Fill(photon->Pt(), fWeightJetJetMC);
  fConvTruePt->Fill(convPhotonMC->Pt(), fWeightJetJetMC);
  if(tmptag<2){
      fConvTrueRecPtTaggedCalo->Fill(photon->Pt(), fWeightJetJetMC);
      fConvTruePtTaggedCalo->Fill(convPhotonMC->Pt(), fWeightJetJetMC);
  } else{
      fConvTrueRecPtTaggedAsDecayCalo->Fill(photon->Pt(), fWeightJetJetMC);
      fConvTruePtTaggedAsDecayCalo->Fill(convPhotonMC->Pt(), fWeightJetJetMC);
  }
  if(isDecayPhoton){
    fConvTruePtDecay->Fill(convPhotonMC->Pt(), fWeightJetJetMC);
    fConvTrueRecPtDecay->Fill(photon->Pt(), fWeightJetJetMC);

    // checkout mother
    Int_t labelMother = convPhotonMC->GetMother();
    AliAODMCParticle *convPhotonMother = (AliAODMCParticle *) fAODMCTrackArray->At(labelMother);
    Int_t nDaughters = convPhotonMother->GetNDaughters();
    Int_t otherDaughterLabel = -1;
    for (Int_t d = 0; d < nDaughters; d++)
    {
        Int_t tmp = convPhotonMother->GetDaughterLabel(d);
        if (tmp == pconvlabel)
            continue;
        otherDaughterLabel = tmp;
    }
    if (otherDaughterLabel != -1)
    {
        if (IsInEMCalAcceptance((AliAODMCParticle*)fAODMCTrackArray->At(otherDaughterLabel)))
        {
            fConvTruePtDecayOtherInAcc->Fill(convPhotonMC->Pt(), fWeightJetJetMC);
            fConvTrueRecPtDecayOtherInAcc->Fill(photon->Pt(), fWeightJetJetMC);
            if (((AliAODMCParticle*)fAODMCTrackArray->At(otherDaughterLabel))->E() >= 0.7){
                fConvTruePtDecayOtherInAccAboveMinEnergy->Fill(convPhotonMC->Pt(), fWeightJetJetMC);
                fConvTrueRecPtDecayOtherInAccAboveMinEnergy->Fill(photon->Pt(), fWeightJetJetMC);
            }
        }
        Int_t clusterLabel = CheckClustersForMCContribution(otherDaughterLabel, fClusterEMCalCandidatesTagging);
        if (clusterLabel != -1)
        {
            fConvTruePtDecayFoundOtherInCluster->Fill(convPhotonMC->Pt(), fWeightJetJetMC);
            fConvTrueRecPtDecayFoundOtherInCluster->Fill(photon->Pt(), fWeightJetJetMC);
        }
    }


  } else if (isSignal){
    fConvTruePtPrimary->Fill(convPhotonMC->Pt(), fWeightJetJetMC);
    fConvTrueRecPtPrimary->Fill(photon->Pt(), fWeightJetJetMC);
  }

  // charged
  for (UInt_t i = 0; i < isoCharged.size(); i++)
  {
    if(i<5){
        fConvTrueIsoCharged[i]->Fill(isoCharged.at(i),photon->Pt(), fWeightJetJetMC);
        if(isDecayPhoton){
          fConvTrueIsoCharged_FromDecay[i]->Fill(isoCharged.at(i),photon->Pt(), fWeightJetJetMC);
        } else if (isSignal){
          fConvTrueIsoCharged_FromDirect[i]->Fill(isoCharged.at(i),photon->Pt(), fWeightJetJetMC);
        }
    }
  }

  // neutral
  for (UInt_t i = 0; i < isoNeutral.size(); i++)
  {
    if(i<5){
        fConvTrueIsoNeutral[i]->Fill(isoNeutral.at(i),photon->Pt(), fWeightJetJetMC);
        if(isDecayPhoton){
          fConvTrueIsoNeutral_FromDecay[i]->Fill(isoNeutral.at(i),photon->Pt(), fWeightJetJetMC);
        } else if (isSignal){
          fConvTrueIsoNeutral_FromDirect[i]->Fill(isoNeutral.at(i),photon->Pt(), fWeightJetJetMC);
        }
    }
  }
  // cell
  for (UInt_t i = 0; i < isoCell.size(); i++)
  {
    if(i<5){
        fConvTrueIsoCell[i]->Fill(isoCell.at(i),photon->Pt(), fWeightJetJetMC);
        if(isDecayPhoton){
          fConvTrueIsoCell_FromDecay[i]->Fill(isoCell.at(i),photon->Pt(), fWeightJetJetMC);
        } else if (isSignal){
          fConvTrueIsoCell_FromDirect[i]->Fill(isoCell.at(i),photon->Pt(), fWeightJetJetMC);
        }
    }
  }
  // full
  for (UInt_t i = 0; i < isoNeutral.size(); i++)
  {
    if(i<5){
      fConvTrueIsoFull[i]->Fill(isoNeutral.at(i) + isoCharged.at(i),photon->Pt(), fWeightJetJetMC);
      if(isDecayPhoton){
        fConvTrueIsoFull_FromDecay[i]->Fill(isoNeutral.at(i) + isoCharged.at(i),photon->Pt(), fWeightJetJetMC);
      } else if (isSignal){
        fConvTrueIsoFull_FromDirect[i]->Fill(isoNeutral.at(i) + isoCharged.at(i),photon->Pt(), fWeightJetJetMC);
      }
    }
  }

// Pt Isolated and tagged

  for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
  {
    for (UInt_t e = 0; e < fTrackIsolationE.size(); e++)
    {
      if(isoCharged.at(r) < fTrackIsolationE.at(e)){
        fConvTruePtIsoCharged[r][e]->Fill(convPhotonMC->Pt(),fWeightJetJetMC);
        if(tmptag<2 )fConvTruePtTaggedCaloIsoCharged[r][e]->Fill(convPhotonMC->Pt(),fWeightJetJetMC); // not tagged by calo
        fConvTrueRecPtIsoCharged[r][e]->Fill(photon->Pt(),fWeightJetJetMC);
        if(tmptag<2 )fConvTrueRecPtTaggedCaloIsoCharged[r][e]->Fill(photon->Pt(),fWeightJetJetMC); // not tagged by calo
      }

      if(isSignal){
         if(isoCharged.at(r) < fTrackIsolationE.at(e)){
          fConvTruePtIsoChargedFromDirect[r][e]->Fill(convPhotonMC->Pt(),fWeightJetJetMC);
          fConvTrueRecPtIsoChargedFromDirect[r][e]->Fill(photon->Pt(),fWeightJetJetMC);
          if(tmptag<2 )fConvTruePtTaggedCaloIsoChargedFromDirect[r][e]->Fill(convPhotonMC->Pt(),fWeightJetJetMC); // not tagged by calo
         }
      }

      // isolation on gen level
      if(mcIsoCharged.at(r) < fTrackIsolationE.at(e)){
        fConvTruePtMCIsoCharged[r][e]->Fill(convPhotonMC->Pt(),fWeightJetJetMC);
        if(tmptag<2 )fConvTruePtTaggedCaloMCIsoCharged[r][e]->Fill(convPhotonMC->Pt(),fWeightJetJetMC); // not tagged by calo
        fConvTrueRecPtMCIsoCharged[r][e]->Fill(photon->Pt(),fWeightJetJetMC);
        if(tmptag<2 )fConvTrueRecPtTaggedCaloMCIsoCharged[r][e]->Fill(photon->Pt(),fWeightJetJetMC); // not tagged by calo
      }
    }

  }

  for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
  {
    for (UInt_t e = 0; e < fNeutralIsolationE.size(); e++)
    {
      if(isoNeutral.at(r) < fNeutralIsolationE.at(e))fConvTruePtIsoNeutral[r][e]->Fill(convPhotonMC->Pt(),fWeightJetJetMC);
      if(isoCell.at(r) < fNeutralIsolationE.at(e)) fConvTruePtIsoCell[r][e]->Fill(convPhotonMC->Pt(),fWeightJetJetMC);
      if((isoNeutral.at(r) + isoCharged.at(r)) < fNeutralIsolationE.at(e)) fConvTruePtIsoFull[r][e]->Fill(convPhotonMC->Pt(),fWeightJetJetMC);

      if(isoNeutral.at(r) < fNeutralIsolationE.at(e))fConvTrueRecPtIsoNeutral[r][e]->Fill(photon->Pt(),fWeightJetJetMC);
      if(isoCell.at(r) < fNeutralIsolationE.at(e)) fConvTrueRecPtIsoCell[r][e]->Fill(photon->Pt(),fWeightJetJetMC);
      if((isoNeutral.at(r) + isoCharged.at(r)) < fNeutralIsolationE.at(e)) fConvTrueRecPtIsoFull[r][e]->Fill(photon->Pt(),fWeightJetJetMC);

      if(tmptag<2){
          if(isoNeutral.at(r) < fNeutralIsolationE.at(e))fConvTruePtTaggedCaloIsoNeutral[r][e]->Fill(convPhotonMC->Pt(),fWeightJetJetMC);
          if(isoCell.at(r) < fNeutralIsolationE.at(e)) fConvTruePtTaggedCaloIsoCell[r][e]->Fill(convPhotonMC->Pt(),fWeightJetJetMC);
          if((isoNeutral.at(r) + isoCharged.at(r)) < fNeutralIsolationE.at(e)) fConvTruePtTaggedCaloIsoFull[r][e]->Fill(convPhotonMC->Pt(),fWeightJetJetMC);

          if(isoNeutral.at(r) < fNeutralIsolationE.at(e))fConvTrueRecPtTaggedCaloIsoNeutral[r][e]->Fill(photon->Pt(),fWeightJetJetMC);
          if(isoCell.at(r) < fNeutralIsolationE.at(e)) fConvTrueRecPtTaggedCaloIsoCell[r][e]->Fill(photon->Pt(),fWeightJetJetMC);
          if((isoNeutral.at(r) + isoCharged.at(r)) < fNeutralIsolationE.at(e)) fConvTrueRecPtTaggedCaloIsoFull[r][e]->Fill(photon->Pt(),fWeightJetJetMC);
      }

      // from direct
      if(isSignal){
        if(isoNeutral.at(r) < fNeutralIsolationE.at(e))fConvTruePtIsoNeutralFromDirect[r][e]->Fill(convPhotonMC->Pt(),fWeightJetJetMC);
        if(isoCell.at(r) < fNeutralIsolationE.at(e)) fConvTruePtIsoCellFromDirect[r][e]->Fill(convPhotonMC->Pt(),fWeightJetJetMC);
        if((isoNeutral.at(r) + isoCharged.at(r)) < fNeutralIsolationE.at(e)) fConvTruePtIsoFullFromDirect[r][e]->Fill(convPhotonMC->Pt(),fWeightJetJetMC);
        if(isoNeutral.at(r) < fNeutralIsolationE.at(e))fConvTrueRecPtIsoNeutralFromDirect[r][e]->Fill(photon->Pt(),fWeightJetJetMC);
        if(isoCell.at(r) < fNeutralIsolationE.at(e)) fConvTrueRecPtIsoCellFromDirect[r][e]->Fill(photon->Pt(),fWeightJetJetMC);
        if((isoNeutral.at(r) + isoCharged.at(r)) < fNeutralIsolationE.at(e)) fConvTrueRecPtIsoFullFromDirect[r][e]->Fill(photon->Pt(),fWeightJetJetMC);

        if(tmptag<2){
          if(isoNeutral.at(r) < fNeutralIsolationE.at(e))fConvTruePtTaggedCaloIsoNeutralFromDirect[r][e]->Fill(convPhotonMC->Pt(),fWeightJetJetMC);
          if(isoCell.at(r) < fNeutralIsolationE.at(e)) fConvTruePtTaggedCaloIsoCellFromDirect[r][e]->Fill(convPhotonMC->Pt(),fWeightJetJetMC);
          if((isoNeutral.at(r) + isoCharged.at(r)) < fNeutralIsolationE.at(e)) fConvTruePtTaggedCaloIsoFullFromDirect[r][e]->Fill(convPhotonMC->Pt(),fWeightJetJetMC);
        }
      }

      // Isolation on gen level
      if(mcIsoNeutral.at(r) < fNeutralIsolationE.at(e))fConvTruePtMCIsoNeutral[r][e]->Fill(convPhotonMC->Pt(),fWeightJetJetMC);
      if(mcIsoFull.at(r) < fNeutralIsolationE.at(e)) fConvTruePtMCIsoFull[r][e]->Fill(convPhotonMC->Pt(),fWeightJetJetMC);

      if(mcIsoNeutral.at(r) < fNeutralIsolationE.at(e))fConvTrueRecPtMCIsoNeutral[r][e]->Fill(photon->Pt(),fWeightJetJetMC);
      if(mcIsoFull.at(r) < fNeutralIsolationE.at(e)) fConvTrueRecPtMCIsoFull[r][e]->Fill(photon->Pt(),fWeightJetJetMC);

      if(tmptag<2){
          if(mcIsoNeutral.at(r) < fNeutralIsolationE.at(e))fConvTruePtTaggedCaloMCIsoNeutral[r][e]->Fill(convPhotonMC->Pt(),fWeightJetJetMC);
          if(mcIsoFull.at(r) < fNeutralIsolationE.at(e)) fConvTruePtTaggedCaloMCIsoFull[r][e]->Fill(convPhotonMC->Pt(),fWeightJetJetMC);

          if(mcIsoNeutral.at(r) < fNeutralIsolationE.at(e))fConvTrueRecPtTaggedCaloMCIsoNeutral[r][e]->Fill(photon->Pt(),fWeightJetJetMC);
          if(mcIsoFull.at(r) < fNeutralIsolationE.at(e)) fConvTrueRecPtTaggedCaloMCIsoFull[r][e]->Fill(photon->Pt(),fWeightJetJetMC);
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskGammaIsoTree::ProcessMCCaloPhoton(AliAODCaloCluster* clus,AliAODConversionPhoton* photon,vector<Double32_t> isoCharged,vector<Double32_t> isoNeutral,vector<Double32_t> isoCell,Int_t tmptag, Double_t weight){
  if(!clus) return;

  Bool_t isTruePhoton = kFALSE;
  Bool_t isTruePhotonWithoutConv = kFALSE; // old definition as crosscheck
  Bool_t isDecay = kFALSE;
  Bool_t isPrompt = kFALSE;
  Bool_t isFrag = kFALSE;
  Bool_t isSignal = kFALSE;
  AliAODMCParticle *MCPhoton = NULL;
  AliAODMCParticle *MCPhotonOld = NULL;

  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  Int_t photonlabel = -1;
  if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (fAODMCTrackArray){
      if (photon->GetIsCaloPhoton() == 0) AliFatal("CaloPhotonFlag has not been set task will abort");
      if (photon->GetNCaloPhotonMCLabels()>0) {
        photonlabel = photon->GetCaloPhotonMCLabel(0);
        MCPhoton = (AliAODMCParticle*) fAODMCTrackArray->At(photonlabel);
        MCPhotonOld = (AliAODMCParticle*) fAODMCTrackArray->At(photonlabel);
      }else{ return;}
  } else {
    AliInfo("AODMCTrackArray could not be loaded");
    return;
  }
  Bool_t isPrimary = fEventCuts->IsConversionPrimaryAOD(fInputEvent, MCPhoton, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
  if (photon->IsLargestComponentPhoton() || (photon->IsLargestComponentElectron() && photon->IsConversion())) {
      isTruePhoton = kTRUE;
      if(isPrimary && photon->IsLargestComponentPhoton()) isTruePhotonWithoutConv = kTRUE;


      // check for proper label, in case of conversion return photon label,otherwise do nothing
      Int_t tmplabel = GetProperLabel(MCPhoton);
      if(tmplabel != -99) photonlabel = tmplabel;

      if((tmplabel != -99) && (!photon->IsConversion())) cout << "photonLabel = " << photonlabel << " motherlabel=" <<  MCPhoton->GetMother() << endl;
      MCPhoton = (AliAODMCParticle*) fAODMCTrackArray->At(photonlabel);

      // if i accept conversion i have to ask about mothers here
      isDecay = IsDecayPhoton(photonlabel);
      isPrompt = IsPromptPhoton(photonlabel);
      isFrag = IsFragPhoton(photonlabel);
      isSignal = isPrompt || isFrag;
  }
  
  if(isTruePhoton){
    fCaloTrueRecPtNoIsPrimary->Fill(photon->Pt(), weight);
    fCaloTruePtNoIsPrimary->Fill(MCPhoton->Pt(), weight);
  }
  
  if(!isTruePhoton || !isPrimary) return;

  // if(MCPhoton->Pt()!=MCPhotonOld->Pt())cout << "PtRec = " << photon->Pt() << "   PtTrueProper =" << MCPhoton->Pt() << "PDG="<<MCPhoton->GetPdgCode() << " PtTrueNoProper=" << MCPhotonOld->Pt() << " PDGCode="<<MCPhotonOld->GetPdgCode() << endl;

  if(isTruePhoton){
     if( ! fClusterCutsEMC->ClusterIsSelectedAODMC(MCPhoton,fAODMCTrackArray)){
         // true photon but would not have been selected on gen levelk
         cout << "true photon but would not have been selected on gen level" << endl;
     }
  }

  isoValues mcIso;
  vector<Double32_t> mcIsoCharged;
  vector<Double32_t> mcIsoNeutral;
  vector<Double32_t> mcIsoFull;
  mcIso = ProcessMCIsolation(photonlabel);

  for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
  {
    Double_t isoCh = mcIso.isolationCone.at(r) - fChargedRhoTimesArea[r];
    mcIsoCharged.push_back(isoCh);
  }
  for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
  {
    mcIsoNeutral.push_back(mcIso.isolationCone.at(r+fTrackIsolationR.size()));
  }
  for (UInt_t r = 0; r < mcIsoCharged.size(); r++)
  {
    mcIsoFull.push_back(mcIsoCharged.at(r)+ mcIsoNeutral.at(r)); // they should always have same length
  }

  fCaloTrueRecPt->Fill(photon->Pt(), weight);
  fCaloTruePt->Fill(MCPhoton->Pt(), weight);
  fCaloTruePtNotProper->Fill(MCPhotonOld->Pt(), weight);

  fCaloTrueRecPtvsTruePt->Fill(MCPhoton->Pt(),photon->Pt(),weight);
  fCaloTrueRecPtvsTruePtNotProper->Fill(MCPhotonOld->Pt(),photon->Pt(),weight);

  if(isTruePhotonWithoutConv) fCaloTrueWithoutConvRecPt->Fill(photon->Pt(), weight);
  if(isTruePhotonWithoutConv) fCaloTrueWithoutConvPt->Fill(MCPhoton->Pt(), weight);

  if(tmptag<2){
      fCaloTrueRecPtTaggedCalo->Fill(photon->Pt(), weight);
      fCaloTruePtTaggedCalo->Fill(MCPhoton->Pt(), weight);
  } else{
      fCaloTrueRecPtTaggedAsDecayCalo->Fill(photon->Pt(), weight);
      fCaloTruePtTaggedAsDecayCalo->Fill(MCPhoton->Pt(), weight);
  }
  if(isDecay){
    fCaloTruePtDecay->Fill(MCPhoton->Pt(), weight);
    fCaloTrueRecPtDecay->Fill(photon->Pt(), weight);

    // checkout mother
    Int_t labelMother = MCPhoton->GetMother();
    AliAODMCParticle *calophotonMother = (AliAODMCParticle *) fAODMCTrackArray->At(labelMother);
    Int_t nDaughters = calophotonMother->GetNDaughters();
    Int_t otherDaughterLabel = -1;
    for (Int_t d = 0; d < nDaughters; d++)
    {
        Int_t tmp = calophotonMother->GetDaughterLabel(d);
        if (tmp == photon->GetCaloPhotonMCLabel(0))
            continue;
        otherDaughterLabel = tmp;
    }
    if (otherDaughterLabel != -1)
    {
        if (IsInEMCalAcceptance((AliAODMCParticle*)fAODMCTrackArray->At(otherDaughterLabel)))
        {
            fCaloTruePtDecayOtherInAcc->Fill(MCPhoton->Pt(), weight);
            fCaloTrueRecPtDecayOtherInAcc->Fill(photon->Pt(), weight);
            if (((AliAODMCParticle*)fAODMCTrackArray->At(otherDaughterLabel))->E() >= 0.7){
                fCaloTruePtDecayOtherInAccAboveMinEnergy->Fill(MCPhoton->Pt(), weight);
                fCaloTrueRecPtDecayOtherInAccAboveMinEnergy->Fill(photon->Pt(), weight);
            }
        }
        Int_t clusterLabel = CheckClustersForMCContribution(otherDaughterLabel, fClusterEMCalCandidatesTagging);
        if (clusterLabel != -1)
        {
            fCaloTruePtDecayFoundOtherInCluster->Fill(MCPhoton->Pt(), weight);
            fCaloTrueRecPtDecayFoundOtherInCluster->Fill(photon->Pt(), weight);
        }
    }


  } else if (isSignal){
    fCaloTruePtPrimary->Fill(MCPhoton->Pt(), weight);
    fCaloTrueRecPtPrimary->Fill(photon->Pt(), weight);
  }

  // charged
  for (UInt_t i = 0; i < isoCharged.size(); i++)
  {
    if(i<5){
        fCaloTrueIsoCharged[i]->Fill(isoCharged.at(i),photon->Pt(), weight);

        fCaloTrueIsoMCCharged[i]->Fill(mcIsoCharged.at(i),photon->Pt(), weight);
        if(isDecay){
          fCaloTrueIsoCharged_FromDecay[i]->Fill(isoCharged.at(i),photon->Pt(), weight);
        } else if (isSignal){
          fCaloTrueIsoCharged_FromDirect[i]->Fill(isoCharged.at(i),photon->Pt(), weight);
        }
    }
  }

  // neutral
  for (UInt_t i = 0; i < isoNeutral.size(); i++)
  {
    if(i<5){
        fCaloTrueIsoNeutral[i]->Fill(isoNeutral.at(i),photon->Pt(), weight);
        if(isDecay){
          fCaloTrueIsoNeutral_FromDecay[i]->Fill(isoNeutral.at(i),photon->Pt(), weight);
        } else if (isSignal){
          fCaloTrueIsoNeutral_FromDirect[i]->Fill(isoNeutral.at(i),photon->Pt(), weight);
        }
    }
  }
  // cell
  for (UInt_t i = 0; i < isoCell.size(); i++)
  {
    if(i<5){
        fCaloTrueIsoCell[i]->Fill(isoCell.at(i),photon->Pt(), weight);
        if(isDecay){
          fCaloTrueIsoCell_FromDecay[i]->Fill(isoCell.at(i),photon->Pt(), weight);
        } else if (isSignal){
          fCaloTrueIsoCell_FromDirect[i]->Fill(isoCell.at(i),photon->Pt(), weight);
        }
    }
  }
  // full
  for (UInt_t i = 0; i < isoNeutral.size(); i++)
  {
    if(i<5){
      fCaloTrueIsoFull[i]->Fill(isoNeutral.at(i) + isoCharged.at(i),photon->Pt(), weight);
      if(isDecay){
        fCaloTrueIsoFull_FromDecay[i]->Fill(isoNeutral.at(i) + isoCharged.at(i),photon->Pt(), weight);
      } else if(isSignal){
        fCaloTrueIsoFull_FromDirect[i]->Fill(isoNeutral.at(i) + isoCharged.at(i),photon->Pt(), weight);
      }
    }
  }


  // Pt Isolated and tagged

  for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
  {
    for (UInt_t e = 0; e < fTrackIsolationE.size(); e++)
    {
      if(isoCharged.at(r) < fTrackIsolationE.at(e)){
        fCaloTruePtIsoCharged[r][e]->Fill(MCPhoton->Pt(),weight);
        if(tmptag<2 )fCaloTruePtTaggedCaloIsoCharged[r][e]->Fill(MCPhoton->Pt(),weight); // not tagged by calo
        fCaloTrueRecPtIsoCharged[r][e]->Fill(photon->Pt(),weight);
        if(tmptag<2 )fCaloTrueRecPtTaggedCaloIsoCharged[r][e]->Fill(photon->Pt(),weight); // not tagged by calo

        if(isSignal){
           fCaloTruePtIsoChargedFromDirect[r][e]->Fill(MCPhoton->Pt(),weight);
           fCaloTrueRecPtIsoChargedFromDirect[r][e]->Fill(photon->Pt(),weight);
           if(tmptag<2 )fCaloTruePtTaggedCaloIsoChargedFromDirect[r][e]->Fill(MCPhoton->Pt(),weight); // not tagged by calo
           if(mcIsoCharged.at(r) < fTrackIsolationE.at(e)){
               fCaloTrueRecPtIsoChargedAndMCIsoChargedFromDirect[r][e]->Fill(photon->Pt(),weight);
               fCaloTruePtIsoChargedAndMCIsoChargedFromDirect[r][e]->Fill(MCPhoton->Pt(),weight);
           }
        }

      }

      // iso on gen level
      if(mcIsoCharged.at(r) < fTrackIsolationE.at(e)){
        fCaloTruePtMCIsoCharged[r][e]->Fill(MCPhoton->Pt(),weight);
        if(tmptag<2 )fCaloTruePtTaggedCaloMCIsoCharged[r][e]->Fill(MCPhoton->Pt(),weight); // not tagged by calo
        fCaloTrueRecPtMCIsoCharged[r][e]->Fill(photon->Pt(),weight);
        if(tmptag<2 )fCaloTrueRecPtTaggedCaloMCIsoCharged[r][e]->Fill(photon->Pt(),weight); // not tagged by calo
      }
    }

  }

  for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
  {
    for (UInt_t e = 0; e < fNeutralIsolationE.size(); e++)
    {
      if(isoNeutral.at(r) < fNeutralIsolationE.at(e))fCaloTruePtIsoNeutral[r][e]->Fill(MCPhoton->Pt(),weight);
      if(isoCell.at(r) < fNeutralIsolationE.at(e)) fCaloTruePtIsoCell[r][e]->Fill(MCPhoton->Pt(),weight);
      if((isoNeutral.at(r) + isoCharged.at(r)) < fNeutralIsolationE.at(e)) fCaloTruePtIsoFull[r][e]->Fill(MCPhoton->Pt(),weight);

      if(isoNeutral.at(r) < fNeutralIsolationE.at(e))fCaloTrueRecPtIsoNeutral[r][e]->Fill(photon->Pt(),weight);
      if(isoCell.at(r) < fNeutralIsolationE.at(e)) fCaloTrueRecPtIsoCell[r][e]->Fill(photon->Pt(),weight);
      if((isoNeutral.at(r) + isoCharged.at(r)) < fNeutralIsolationE.at(e)) fCaloTrueRecPtIsoFull[r][e]->Fill(photon->Pt(),weight);

      if(tmptag<2){
          if(isoNeutral.at(r) < fNeutralIsolationE.at(e))fCaloTruePtTaggedCaloIsoNeutral[r][e]->Fill(MCPhoton->Pt(),weight);
          if(isoCell.at(r) < fNeutralIsolationE.at(e)) fCaloTruePtTaggedCaloIsoCell[r][e]->Fill(MCPhoton->Pt(),weight);
          if((isoNeutral.at(r) + isoCharged.at(r)) < fNeutralIsolationE.at(e)) fCaloTruePtTaggedCaloIsoFull[r][e]->Fill(MCPhoton->Pt(),weight);

          if(isoNeutral.at(r) < fNeutralIsolationE.at(e))fCaloTrueRecPtTaggedCaloIsoNeutral[r][e]->Fill(photon->Pt(),weight);
          if(isoCell.at(r) < fNeutralIsolationE.at(e)) fCaloTrueRecPtTaggedCaloIsoCell[r][e]->Fill(photon->Pt(),weight);
          if((isoNeutral.at(r) + isoCharged.at(r)) < fNeutralIsolationE.at(e)) fCaloTrueRecPtTaggedCaloIsoFull[r][e]->Fill(photon->Pt(),weight);
      }

      if(isSignal){
        if(isoNeutral.at(r) < fNeutralIsolationE.at(e))fCaloTruePtIsoNeutralFromDirect[r][e]->Fill(MCPhoton->Pt(),weight);
        if(isoCell.at(r) < fNeutralIsolationE.at(e)) fCaloTruePtIsoCellFromDirect[r][e]->Fill(MCPhoton->Pt(),weight);
        if((isoNeutral.at(r) + isoCharged.at(r)) < fNeutralIsolationE.at(e)) fCaloTruePtIsoFullFromDirect[r][e]->Fill(MCPhoton->Pt(),weight);
          if(isoNeutral.at(r) < fNeutralIsolationE.at(e))fCaloTrueRecPtIsoNeutralFromDirect[r][e]->Fill(photon->Pt(),weight);
        if(isoCell.at(r) < fNeutralIsolationE.at(e)) fCaloTrueRecPtIsoCellFromDirect[r][e]->Fill(photon->Pt(),weight);
        if((isoNeutral.at(r) + isoCharged.at(r)) < fNeutralIsolationE.at(e)) fCaloTrueRecPtIsoFullFromDirect[r][e]->Fill(photon->Pt(),weight);
        if(tmptag<2){
          if(isoNeutral.at(r) < fNeutralIsolationE.at(e))fCaloTruePtTaggedCaloIsoNeutralFromDirect[r][e]->Fill(MCPhoton->Pt(),weight);
          if(isoCell.at(r) < fNeutralIsolationE.at(e)) fCaloTruePtTaggedCaloIsoCellFromDirect[r][e]->Fill(MCPhoton->Pt(),weight);
          if((isoNeutral.at(r) + isoCharged.at(r)) < fNeutralIsolationE.at(e)) fCaloTruePtTaggedCaloIsoFullFromDirect[r][e]->Fill(MCPhoton->Pt(),weight);
        }
      }

      // iso on gen level
      if(mcIsoNeutral.at(r) < fNeutralIsolationE.at(e))fCaloTruePtMCIsoNeutral[r][e]->Fill(MCPhoton->Pt(),weight);
      if(mcIsoFull.at(r) < fNeutralIsolationE.at(e)) fCaloTruePtMCIsoFull[r][e]->Fill(MCPhoton->Pt(),weight);

      if(mcIsoNeutral.at(r) < fNeutralIsolationE.at(e))fCaloTrueRecPtMCIsoNeutral[r][e]->Fill(photon->Pt(),weight);
      if(mcIsoFull.at(r) < fNeutralIsolationE.at(e)) fCaloTrueRecPtMCIsoFull[r][e]->Fill(photon->Pt(),weight);

      if(tmptag<2){
          if(mcIsoNeutral.at(r) < fNeutralIsolationE.at(e))fCaloTruePtTaggedCaloMCIsoNeutral[r][e]->Fill(MCPhoton->Pt(),weight);
          if(mcIsoFull.at(r) < fNeutralIsolationE.at(e)) fCaloTruePtTaggedCaloMCIsoFull[r][e]->Fill(MCPhoton->Pt(),weight);

          if(mcIsoNeutral.at(r) < fNeutralIsolationE.at(e))fCaloTrueRecPtTaggedCaloMCIsoNeutral[r][e]->Fill(photon->Pt(),weight);
          if(mcIsoFull.at(r) < fNeutralIsolationE.at(e)) fCaloTrueRecPtTaggedCaloMCIsoFull[r][e]->Fill(photon->Pt(),weight);
      }
    }

  }
}
//________________________________________________________________________
void AliAnalysisTaskGammaIsoTree::ProcessCaloPhotons(){
   Int_t nclus                         = 0;

   Int_t posEMC = 0;
   Int_t posEMCIso = 0;
   Int_t posEMCTag = 0;
   Int_t posPHOS = 0;
   Double_t vertex[3] = {0};
   InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

   TClonesArray * arrClustersProcess   = NULL;
   std::vector<Double_t> clusWeights;
   std::vector<Int_t> clusterPos;
   std::vector<Bool_t> isFromDesiredHeader;
   if(!fCorrTaskSetting.CompareTo("")){
     nclus = fInputEvent->GetNumberOfCaloClusters();
   } else {
    arrClustersProcess                = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));
    if(!arrClustersProcess)
      AliFatal(Form("%sClustersBranch was not found in AliAnalysisTaskGammaCalo! Check the correction framework settings!",fCorrTaskSetting.Data()));
    nclus                         = arrClustersProcess->GetEntries();
  }

  if(nclus == 0)  return;
  ((AliCaloPhotonCuts*)fClusterCutsEMC)->FillHistogramsExtendedQA(fInputEvent,fIsMC);
  // ((AliCaloPhotonCuts*)fClusterCutsPHOS)->FillHistogramsExtendedQA(fInputEvent,fIsMC);

  // in case user wants to use default track matching
  AliAODCaloCluster* clus                       = NULL;
 // 1000 in 42 sec
  if(!fDoOwnTrackMatching){
    ((AliCaloPhotonCuts*)fClusterCutsEMCTrackMatching)->MatchTracksToClusters(fInputEvent,fWeightJetJetMC,kTRUE, fMCEvent);
    //((AliCaloPhotonCuts*)fClusterCutsIsolationEMC)->MatchTracksToClusters(fInputEvent,fWeightJetJetMC,kTRUE, fMCEvent);
    //((AliCaloPhotonCuts*)fClusterCutsTaggingEMC)->MatchTracksToClusters(fInputEvent,fWeightJetJetMC,kTRUE, fMCEvent);
  }
  if(arrClustersProcess){
     // EMCal correction framework was used
     // we need to loop over this for EMCal clusters and
     // over the others for PHOS cluster

      // Loop over EMCal clusters
      for(Long_t i = 0; i < nclus; i++){
        clus                                = (AliAODCaloCluster*)arrClustersProcess->At(i);

        if(!clus) continue;

        // Do all needed checks
        Double_t tempClusterWeight        = fWeightJetJetMC;
        Double_t tempPhotonWeight         = fWeightJetJetMC;


        // Set the jetjet weight to 1 in case the cluster orignated from the minimum bias header
        if (fIsMC>0 && (fEventCuts->GetSignalRejection() == 4)){
          if(fEventCuts->IsParticleFromBGEvent(clus->GetLabelAt(0), fMCEvent, fInputEvent) == 2)
            tempClusterWeight = 1;
        }

        // Header check
        fIsFromDesiredHeader          = kTRUE;
        fIsOverlappingWithOtherHeader = kFALSE;
        // test whether largest contribution to cluster orginates in added signals
        if (fIsMC>0 && fEventCuts->GetSignalRejection() > 0){
          Int_t* mclabelsCluster = clus->GetLabels();
          // Set the jetjet weight to 1 in case the photon candidate orignated from the minimum bias header
          if ( fEventCuts->IsParticleFromBGEvent(mclabelsCluster[0], fMCEvent, fInputEvent) == 2 && fEventCuts->GetSignalRejection() == 4) tempPhotonWeight = 1;
          if ( fEventCuts->IsParticleFromBGEvent(mclabelsCluster[0], fMCEvent, fInputEvent) == 0) fIsFromDesiredHeader = kFALSE;
          if (clus->GetNLabels()>1){
            // Int_t* mclabelsCluster = clus->GetLabels();
            // if (fLocalDebugFlag > 1)   cout << "testing if other labels in cluster belong to different header, need to test " << (Int_t)clus->GetNLabels()-1 << " additional labels" << endl;
              for (Int_t l = 1; l < (Int_t)clus->GetNLabels(); l++ ){
                if (fEventCuts->IsParticleFromBGEvent(mclabelsCluster[l], fMCEvent, fInputEvent, 0) == 0) fIsOverlappingWithOtherHeader = kTRUE;
              }
            // if (fLocalDebugFlag > 1 && fIsOverlappingWithOtherHeader) cout << "found overlapping header: " << endl;
          }
        }

        isFromDesiredHeader.push_back(fIsFromDesiredHeader);

        if ( !clus->IsEMCAL()){ // for PHOS: cluster->GetType() == AliVCluster::kPHOSNeutral
          continue;
        }



        if ( (fIsFromDesiredHeader && !fIsOverlappingWithOtherHeader && !fAllowOverlapHeaders) || (fIsFromDesiredHeader && fAllowOverlapHeaders) ){

          // get additional cluster info
          Short_t nLM = ((AliCaloPhotonCuts*)fClusterCutsEMC)->GetNumberOfLocalMaxima(clus, fInputEvent);
          Short_t matchIndex = -1;
          // if(fDoOwnTrackMatching){
          //     matchIndex = ProcessTrackMatching(clus,fTracks);
          // }
          Float_t eFrac = GetExoticEnergyFraction(clus,fInputEvent);
          // if(((AliCaloPhotonCuts*)fClusterCutsIsolationEMC)->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC, tempClusterWeight,i)){
          //   if(!IsMatchedWithConv(clus,fClusterCutsIsolationEMC)){
          //     new((*fExtraClusterInfoBackground)[posEMCIso]) AliExtraClusterInfoHelper(nLM,matchIndex,eFrac);
          //     new((*fClusterEMCalCandidatesIsolation)[posEMCIso]) AliAODCaloCluster(*clus);
          //     posEMCIso++;
          //   }
          // }

          // if(((AliCaloPhotonCuts*)fClusterCutsTaggingEMC)->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC, tempClusterWeight,i)){
          //   if(!IsMatchedWithConv(clus,fClusterCutsTaggingEMC)){
          //     new((*fClusterEMCalCandidatesTagging)[posEMCTag]) AliAODCaloCluster(*clus);
          //     posEMCTag++;
          //   }
          // }

          // check if given EMC cuts are fulfilled
          if(!((AliCaloPhotonCuts*)fClusterCutsEMC)->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC, tempClusterWeight,i)){
            continue;
          }
          // if(IsMatchedWithConv(clus,fClusterCutsEMC)){
          //   delete clus;
          //   continue;
          // }
          fClusterEMCalCandidates->Add(clus);
          clusWeights.push_back(tempPhotonWeight);
          clusterPos.push_back(i);
          isFromDesiredHeader.push_back(fIsFromDesiredHeader);
          posEMC++;
        }
      } // end of initial cluster loop
  }




  // Loop over normal clusters
  for(Long_t i = 0; i < nclus; i++){
    clus                                = (AliAODCaloCluster*)fInputEvent->GetCaloCluster(i);

    if(!clus) continue;

    // Do all needed checks
    Double_t tempClusterWeight        = fWeightJetJetMC;
    Double_t tempPhotonWeight         = fWeightJetJetMC;

    // Set the jetjet weight to 1 in case the cluster orignated from the minimum bias header
    if (fIsMC>0 && (fEventCuts->GetSignalRejection() == 4)){
      if(fEventCuts->IsParticleFromBGEvent(clus->GetLabelAt(0), fMCEvent, fInputEvent) == 2)
        tempClusterWeight = 1;
    }

    // Header check
    fIsFromDesiredHeader          = kTRUE;
    fIsOverlappingWithOtherHeader = kFALSE;
    // test whether largest contribution to cluster orginates in added signals
    if (fIsMC>0 && fEventCuts->GetSignalRejection() > 0){
      Int_t* mclabelsCluster = clus->GetLabels();
      // Set the jetjet weight to 1 in case the photon candidate orignated from the minimum bias header
      if ( fEventCuts->IsParticleFromBGEvent(mclabelsCluster[0], fMCEvent, fInputEvent) == 2 && fEventCuts->GetSignalRejection() == 4) tempPhotonWeight = 1;
      if ( fEventCuts->IsParticleFromBGEvent(mclabelsCluster[0], fMCEvent, fInputEvent) == 0) fIsFromDesiredHeader = kFALSE;
      if (clus->GetNLabels()>1){
        // Int_t* mclabelsCluster = clus->GetLabels();
        // if (fLocalDebugFlag > 1)   cout << "testing if other labels in cluster belong to different header, need to test " << (Int_t)clus->GetNLabels()-1 << " additional labels" << endl;
          for (Int_t l = 1; l < (Int_t)clus->GetNLabels(); l++ ){
            if (fEventCuts->IsParticleFromBGEvent(mclabelsCluster[l], fMCEvent, fInputEvent, 0) == 0) fIsOverlappingWithOtherHeader = kTRUE;
          }
        // if (fLocalDebugFlag > 1 && fIsOverlappingWithOtherHeader) cout << "found overlapping header: " << endl;
      }
    }

    if(!arrClustersProcess && clus->IsEMCAL()){ // if is was not saved already
      if ( (fIsFromDesiredHeader && !fIsOverlappingWithOtherHeader && !fAllowOverlapHeaders) || (fIsFromDesiredHeader && fAllowOverlapHeaders) ){
        // get additional cluster info
        Short_t nLM = ((AliCaloPhotonCuts*)fClusterCutsEMC)->GetNumberOfLocalMaxima(clus, fInputEvent);
        Short_t matchIndex = -1;
        // if(fDoOwnTrackMatching){
        //       matchIndex = ProcessTrackMatching(clus,fTracks);
        // }
        Float_t eFrac = GetExoticEnergyFraction(clus,fInputEvent);

        // if(((AliCaloPhotonCuts*)fClusterCutsIsolationEMC)->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC, tempClusterWeight,i)){
        //   if(!IsMatchedWithConv(clus,fClusterCutsIsolationEMC)){
        //     new((*fExtraClusterInfoBackground)[posEMCIso]) AliExtraClusterInfoHelper(nLM,matchIndex,eFrac);
        //     new((*fClusterEMCalCandidatesIsolation)[posEMCIso]) AliAODCaloCluster(*clus);
        //     posEMCIso++;
        //   }
        // }

        // if(((AliCaloPhotonCuts*)fClusterCutsTaggingEMC)->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC, tempClusterWeight,i)){
        //   if(!IsMatchedWithConv(clus,fClusterCutsTaggingEMC)){
        //     new((*fClusterEMCalCandidatesTagging)[posEMCTag]) AliAODCaloCluster(*clus);
        //     posEMCTag++;
        //   }
        // }
        if(!((AliCaloPhotonCuts*)fClusterCutsEMC)->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC, tempClusterWeight,i)){
          // delete clus;
          continue;
        }
        // if(IsMatchedWithConv(clus,fClusterCutsEMC)){
        //   delete clus;
        //   continue;
        // }
        
        // no track matching at this point, do manually fi tree

        fClusterEMCalCandidates->Add(clus);
        clusWeights.push_back(tempPhotonWeight);
        clusterPos.push_back(i);
        isFromDesiredHeader.push_back(fIsFromDesiredHeader);
        posEMC++;
      }
      continue;
    }
    if(clus->IsPHOS() && fSavePHOSClusters){
    // if(clus->GetType() == AliVCluster::kPHOSNeutral){
      if(!((AliCaloPhotonCuts*)fClusterCutsPHOS)->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC, tempClusterWeight,i)){
        continue;
      }
      // if(!fSavePHOSClusters) new((*fClusterPHOSCandidates)[posPHOS]) AliAODCaloCluster(*clus);
      posPHOS++;
      continue;
    }
  }

  // Loop again over selected cluster candidates to do isolation and tagging
  for (Int_t c = 0; c < fClusterEMCalCandidates->GetEntries(); c++)
  {
     AliAODCaloCluster* clus = (AliAODCaloCluster*) fClusterEMCalCandidates->At(c);
     Double32_t tmp_isoCharged[2] = {0,0};
     Double32_t tmp_isoNeutral[2] = {0,0};
     Double32_t tmp_isoCell[2] = {0,0};
     Int_t      tmp_tag= 0;
     isoValues isoCharged;
     vector<Double32_t> isoNeutral;
     vector<Double32_t> isoCell;
     if(fDoTrackIsolation) isoCharged = ProcessChargedIsolation(clus);
     if(fDoNeutralIsolation){
       isoNeutral =ProcessNeutralIsolation(clus);
     } else{
          for (UInt_t r = 0; r < isoNeutral.size(); r++)
          {
             isoNeutral.push_back(0.);
          }
     }
     if(fDoCellIsolation){
        isoCell = ProcessCellIsolation(clus);
     } else{
          for (UInt_t r = 0; r < isoNeutral.size(); r++)
          {
             isoCell.push_back(0.);
          }
     }
     if(fDoTagging) tmp_tag = ProcessTagging(clus);

     // Subtract correction of charged isolation
     // others not implemented yes


     // when writing to tree, only fill first two R
     if(fUseHistograms){
      for (UInt_t r = 0; r < isoCharged.isolationCone.size(); r++)
      {
        // Fill charged isolation corrected for background
        if(r<2) tmp_isoCharged[r] = isoCharged.isolationCone.at(r) - fChargedRhoTimesArea[r];
        if(r<2) isoCharged.isolationCone.at(r) =  isoCharged.isolationCone.at(r) - fChargedRhoTimesArea[r];
      }
      for (UInt_t r = 0; r < isoNeutral.size(); r++)
      {
        if(r<2){
          tmp_isoNeutral[r] = isoNeutral.at(r);
          tmp_isoCell[r] = isoCell.at(r);
        }
      }
     }

     // convert to AliAODCOnversionPhoton
     // TLorentzvector with cluster
     TLorentzVector clusterVector;
     clus->GetMomentum(clusterVector,vertex);


     TLorentzVector* tmpvec = new TLorentzVector();
     tmpvec->SetPxPyPzE(clusterVector.Px(),clusterVector.Py(),clusterVector.Pz(),clusterVector.E());
     // convert to AODConversionPhoton
     AliAODConversionPhoton *PhotonCandidate=new AliAODConversionPhoton(tmpvec);
     if(!PhotonCandidate){ delete tmpvec; continue;}

     // Flag Photon as CaloPhoton
     PhotonCandidate->SetIsCaloPhoton(fClusterCutsEMC->GetClusterType());
     PhotonCandidate->SetCaloClusterRef(clusterPos.at(c));
     PhotonCandidate->SetLeadingCellID(fClusterCutsEMC->FindLargestCellInCluster(clus,fInputEvent));
     // get MC label
     if(fIsMC> 0){
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
       if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
       PhotonCandidate->SetCaloPhotonMCFlagsAOD(fAODMCTrackArray, kFALSE);

     }

     // Only consider clusters that allow for a cone of 0.4 with tracking
     Bool_t isWithinTPC =IsWithinRadiusTPC(PhotonCandidate->Eta(),PhotonCandidate->Phi(),fExclusionRadius);
    if(!isWithinTPC){
      if(tmpvec)  delete tmpvec;
      if(PhotonCandidate) delete PhotonCandidate;
      continue;
     }
     if(fUseHistograms) FillCaloHistosPurity(clus,PhotonCandidate,isoCharged.isolationCone,isoNeutral,isoCell,tmp_tag,clusWeights.at(c));
     
     if(!fUseHistograms) {
       // only fill tree for clusters fulfilling cut
       if(PhotonCandidate->Pt()>=fRecPtCut) FillCaloTree(clus,PhotonCandidate,isoCharged,isoNeutral,isoCell,tmp_tag,clusWeights.at(c));
     }
     if((clus->GetM02() < fMinM02) || (clus->GetM02() > fMaxM02)){
      if(tmpvec)  delete tmpvec;
      if(PhotonCandidate) delete PhotonCandidate;
      continue;
     }

     if(fUseHistograms){
       FillCaloHistos(clus,PhotonCandidate,isoCharged,isoNeutral,isoCell,tmp_tag,clusWeights.at(c));
     } else{ // fill even in case of tree for trigger rejection
       fCaloPt->Fill(PhotonCandidate->Pt(),clusWeights.at(c));
     }
     if((fIsMC>0) && fUseHistograms) ProcessMCCaloPhoton(clus,PhotonCandidate,isoCharged.isolationCone,isoNeutral,isoCell,tmp_tag,clusWeights.at(c));
     if(PhotonCandidate) delete PhotonCandidate;
     if(tmpvec)  delete tmpvec;
  }
  clusWeights.clear();
  clusterPos.clear();
  isFromDesiredHeader.clear();
}

///________________________________________________________________________
Bool_t AliAnalysisTaskGammaIsoTree::TrackIsSelectedAOD(AliAODTrack* lTrack) {
  // apply filter bits
  if( ! lTrack->IsHybridGlobalConstrainedGlobal()){
    return kFALSE;
  }

	// Absolute TPC Cluster cut
	if(lTrack->GetTPCNcls()<fMinClsTPC) return kFALSE;
	if(lTrack->GetTPCchi2perCluster()>fChi2PerClsTPC) return kFALSE;
  // DCA cut
  //if(!IsDCACutAccepted(lTrack)) return kFALSE;

  // ITS Cluster Cut
	// SetClusterRequirementITS and SetRequireITSRefit can
	// not be set for AODs after filtering
	if(lTrack->GetITSNcls()<fMinClsITS) return kFALSE;

  if(  TMath::Abs(lTrack->Eta()) > fEtaCut ) {
      return kFALSE;
  }

  if( lTrack->Pt() < fPtCut ) {
    return kFALSE;
  }
  return kTRUE;
}

//_____________________________________________________________________________
void AliAnalysisTaskGammaIsoTree::ProcessTracks(){
  Int_t pos = 0;
  for(Int_t t=0;t<fInputEvent->GetNumberOfTracks();t++){
      AliAODTrack *fCurrentTrack = static_cast<AliAODTrack*> (fInputEvent->GetTrack(t));
      //if(!TrackIsSelectedAOD(fCurrentTrack)){
        // save empty track to preserve position
      //  new((*fTracks)[pos]) AliAODTrack();
      //} else{

      // we need to save all tracks in order to identify tracks from conv for isolation
      // new((*fTracks)[pos]) AliAODTrack(*fCurrentTrack);
      // pos++;
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskGammaIsoTree::ProcessMCParticles(){
  // Loop over all primary MC particle
    const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  Int_t pos = 0;
  if (fAODMCTrackArray){
    for(Int_t i = 0; i < fAODMCTrackArray->GetEntriesFast(); i++) {
      Double_t tempParticleWeight       = fWeightJetJetMC;
      AliAODMCParticle* particle =  static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(i));
      if(!particle) continue;

      if(!fUseHistograms && (particle->Pt()<fGenPtCut)) continue;
      // if(TMath::Abs(particle->Y())< fYMCCut){
      //   new((*fMCParticles)[pos]) AliAODMCParticle(*particle);
      //   pos++;

      // } else
      // {
      //   new((*fMCParticles)[pos]) AliAODMCParticle();
      //   pos++;
      // }
      // check if primary
      Bool_t isPrimary = fEventCuts->IsConversionPrimaryAOD(fInputEvent, particle, mcProdVtxX, mcProdVtxY, mcProdVtxZ);


      // check header
      Int_t isMCFromMBHeader = -1;
      if(fEventCuts->GetSignalRejection() != 0){
        isMCFromMBHeader = fEventCuts->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
        if(isMCFromMBHeader == 0 && fEventCuts->GetSignalRejection() != 3) continue;
        // Set the jetjet weight to 1 in case the particle orignated from the minimum bias header
        if(isMCFromMBHeader == 2 && fEventCuts->GetSignalRejection() == 4) {
          tempParticleWeight = 1;
        }
      }

      //
      // ─── CHECK PI0 CASE ──────────────────────────────────────────────
      //
      if(fUseHistograms){
        if(particle->PdgCode() == 111){
          fGenPi0Pt->Fill(particle->Pt(),tempParticleWeight);
          if (IsInEMCalAcceptance(particle))
          {
              fGenPi0PtInEMCalAcc->Fill(particle->Pt(), tempParticleWeight);

              Int_t nDaughters = particle->GetNDaughters();
              Int_t nPhotonsFound = 0;
              Int_t photonLabels[2] = {0};
              if(nDaughters == 2){ // only look at this for simplicity
                  for (Int_t d = 0; d < nDaughters; d++)
                  {
                      Int_t daughterLabel = particle->GetDaughterLabel(d);
                      AliAODMCParticle *photon = static_cast<AliAODMCParticle *>(fAODMCTrackArray->At(daughterLabel));

                      if(photon->GetPdgCode() == 22){
                          nPhotonsFound++;
                          photonLabels[d] = daughterLabel;
                      }
                  }

                  // Decay to two photons
                  if(nPhotonsFound == 2){
                      if (photonLabels[0] < fAODMCTrackArray->GetEntriesFast()){
                          AliAODMCParticle *daughter1 = static_cast<AliAODMCParticle *>(fAODMCTrackArray->At(photonLabels[0]));
                          AliAODMCParticle *daughter2 = static_cast<AliAODMCParticle *>(fAODMCTrackArray->At(photonLabels[1]));
                          if (IsInEMCalAcceptance(daughter1) &&
                              IsInEMCalAcceptance(daughter2)){
                              fGenPi0PtInEMCalAcc_BothGammaInEMCal->Fill(particle->Pt(), tempParticleWeight);
                          }
                          // Int_t clusIndex[2] = {-1,-1};
                          // clusIndex[0] = CheckClustersForMCContribution(photonLabels[0], fClusterEMCalCandidatesTagging);
                          // clusIndex[1] = CheckClustersForMCContribution(photonLabels[1], fClusterEMCalCandidatesTagging);
                          // if ( (clusIndex[0] != -1) && (clusIndex[1] != -1))
                          // {
                          //     fGenPi0PtInEMCalAcc_BothGammaInClusters->Fill(particle->Pt(), tempParticleWeight);
                          // }
                      }
                  }
              }

          }
        }  // end pi0
      }

      //
      // ─── PHOTON CASE ─────────────────────────────────────────────────
      //
      if(particle->PdgCode()!= 22 ) continue;

      if(fUseHistograms) fGenPhotonPt->Fill(particle->Pt(), tempParticleWeight);

      Bool_t isDecay = IsDecayPhoton(i); // gamma as mother would survive this
      Bool_t isPrompt = IsPromptPhoton(i); // gamma as mother would survive this
      Bool_t isFrag = IsFragPhoton(i); // gamma as mother would survive this

      Bool_t isSignal = isPrompt || isFrag;
      if(fUseHistograms){
        if(isDecay){
          fGenPhotonPt_FromDecay->Fill(particle->Pt(), tempParticleWeight);
        } else if(isSignal){
          fGenPhotonPt_FromDirect->Fill(particle->Pt(), tempParticleWeight);
        }
      }
      if(!isPrimary) continue;

      // in EMC acceptance and not gamma as mother to avoid double counting
      // this should only count photon highest up the chain, technically could have wrong pT
      Bool_t isWithinTPC =IsWithinRadiusTPC(particle->Eta(),particle->Phi(),fExclusionRadius);
      if (fClusterCutsEMC->ClusterIsSelectedAODMC(particle,fAODMCTrackArray) && isWithinTPC)
      {
         // if(!particle->IsPhysicalPrimary()) AliInfo("none physical primary selected!");
         isoValues mcIso;
          vector<Double32_t> mcIsoCharged;
          vector<Double32_t> mcIsoNeutral;
          vector<Double32_t> mcIsoFull;
         
          if(fUseHistograms || (!fUseHistograms && (particle->Pt()>=fGenPtCut))){
            mcIso = ProcessMCIsolation(i);

            for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
            {
              Double_t isoCh = mcIso.isolationCone.at(r) - fChargedRhoTimesArea[r];
              mcIsoCharged.push_back(isoCh);
            }
            for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
            {
              mcIsoNeutral.push_back(mcIso.isolationCone.at(r+fTrackIsolationR.size()));
            }
            for (UInt_t r = 0; r < mcIsoCharged.size(); r++)
            {
              mcIsoFull.push_back(mcIsoCharged.at(r)+ mcIsoNeutral.at(r)); // they should always have same length
            }
          }

          // fill for tree (if genPt cut is fulfilled)
          if(!fUseHistograms && (particle->Pt()>=fGenPtCut)){
            TString headerName = fEventCuts->GetParticleHeaderName(i, fMCEvent, fInputEvent);
            if(((AliConvEventCuts*)fEventCuts)->GetSignalRejection() == 0){
              headerName = "";
            }
            Int_t tag = GetMCAnalysisUtils()->CheckOrigin(i, fMCEvent,headerName,1.);
            fBuffer_GenPhotonE.push_back(particle->E());
            fBuffer_GenPhotonPx.push_back(particle->Px());
            fBuffer_GenPhotonPy.push_back(particle->Py());
            fBuffer_GenPhotonPz.push_back(particle->Pz());
            fBuffer_GenPhotonMCIsoCharged1.push_back(mcIso.isolationCone.at(0));
            if(mcIso.isolationCone.size()>1) fBuffer_GenPhotonMCIsoCharged2.push_back(mcIso.isolationCone.at(1));
            if(mcIso.isolationCone.size()>2) fBuffer_GenPhotonMCIsoCharged3.push_back(mcIso.isolationCone.at(2));
            fBuffer_GenPhotonMCIsoBckLeft.push_back(mcIso.backgroundLeft.at(2));
            Bool_t isConv = kFALSE;
            
            Int_t nDaughters = particle->GetNDaughters();
            Int_t nElecsFound = 0;
            // cout << "nDaughter =" << nDaughters << endl;
            if(nDaughters == 2){ // only look at this for simplicity
              for (Int_t d = 0; d < nDaughters; d++)
              {
                  Int_t daughterLabel = particle->GetDaughterLabel(d);
                  AliAODMCParticle *temp = static_cast<AliAODMCParticle *>(fAODMCTrackArray->At(daughterLabel));
                  Double_t Xv = temp->Xv();
                  Double_t Yv = temp->Yv();
                  Double_t Zv = temp->Zv();
                  Double_t R = TMath::Sqrt(Xv*Xv + Yv*Yv + Zv*Zv);
                  if(TMath::Abs(temp->GetPdgCode()) == 11){
                      if(R<430) nElecsFound++; // make sure its not a conversion in EMCal itself that is still on the stack
                  }
              }
            }
            if(nElecsFound==2) isConv=kTRUE;

            fBuffer_GenPhotonIsConv.push_back(isConv);
            fBuffer_GenPhotonMCTag.push_back(tag);
          }


          if(fUseHistograms) fGenPhotonPtInEMCalAcc->Fill(particle->Pt(), tempParticleWeight);
          if(isDecay){
             if(fUseHistograms) fGenPhotonPtInEMCalAcc_FromDecay->Fill(particle->Pt(), tempParticleWeight);
          } else if (isSignal){
              if(fUseHistograms) fGenPhotonPtInEMCalAcc_FromDirect->Fill(particle->Pt(), tempParticleWeight);
              if(fUseHistograms){
                for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
                {
                  fGenPhotonChargedMCIsoInEMCalAcc_FromDirect[r]->Fill(mcIsoCharged.at(r) ,particle->Pt(), tempParticleWeight);
                  for (UInt_t e = 0; e < fTrackIsolationE.size(); e++)
                  { 
                    if(mcIsoCharged.at(r) < fTrackIsolationE.at(e))fGenPhotonPtInEMCalAccChargedMCIso_FromDirect[r][e]->Fill(particle->Pt(), tempParticleWeight);
                  }
                }
              }
          }
      }
    }
  }
}

//_____________________________________________________________________________
Int_t AliAnalysisTaskGammaIsoTree::ProcessTrackMatching(AliAODCaloCluster* clus, TList* tracks){
     Int_t nModules = fGeomEMCAL->GetNumberOfSuperModules();
     Int_t highestMatchIndex = -1;
     AliExternalTrackParam *trackParam = 0;
     for (Int_t t = 0; t < tracks->GetEntries(); t++)
     {
      AliAODTrack *aodt = static_cast<AliAODTrack*>(tracks->At(t));
        if(!aodt) continue;
        if(!aodt->IsHybridGlobalConstrainedGlobal()) continue;
        if(aodt->Pt()<0.5) continue;

        Double_t xyz[3] = {0}, pxpypz[3] = {0}, cv[21] = {0};
        aodt->GetPxPyPz(pxpypz);
        aodt->GetXYZ(xyz);
        aodt->GetCovarianceXYZPxPyPz(cv);

        trackParam = new AliExternalTrackParam(xyz,pxpypz,cv,aodt->Charge());

        AliExternalTrackParam emcParam(*trackParam);
        Float_t eta, phi, pt;
        //propagate tracks to emc surfaces
        if (!AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(&emcParam, 440., 0.139, 20., eta, phi, pt)) {
          delete trackParam;
          continue;
        }
        if( TMath::Abs(eta) > 0.75 ) {
          delete trackParam;
          continue;
        }
        // Save some time and memory in case of no DCal present
        if( nModules < 13 && ( phi < 70*TMath::DegToRad() || phi > 190*TMath::DegToRad())){
          delete trackParam;
          continue;
        }
        // Save some time and memory in case of run2
        if( nModules > 12 ){
          if (( phi < 70*TMath::DegToRad() || phi > 190*TMath::DegToRad()) && ( phi < 250*TMath::DegToRad() || phi > 340*TMath::DegToRad())){
            delete trackParam;
            continue;
          }
        }
        Float_t dEta=-999, dPhi=-999;
        Double_t trkPos[3] = {0.,0.,0.};
        if (!emcParam.GetXYZ(trkPos)){
          delete trackParam;
          continue;
        }

        AliExternalTrackParam trackParamTmp(emcParam);//Retrieve the starting point every time before the extrapolation
        if(!AliEMCALRecoUtils::ExtrapolateTrackToCluster(&trackParamTmp, clus, 0.139, 5., dEta, dPhi)) continue;
        if(TMath::Abs(dEta) > (fMatchingParamsEta[0] + pow(aodt->Pt() + fMatchingParamsEta[1],fMatchingParamsEta[2]))) continue;
        if(TMath::Abs(dPhi) > (fMatchingParamsPhi[0] + pow(aodt->Pt() + fMatchingParamsPhi[1],fMatchingParamsPhi[2]))) continue;
        if((clus->E()/aodt->P()) > fMatchingEOverP) continue;
        if(highestMatchIndex == -1){ // this is the first match
           highestMatchIndex = t;
        } else{
           if(aodt->P()>((AliAODTrack*)tracks->At(highestMatchIndex))->P()){
             highestMatchIndex = t;
           }
        }

        delete trackParam;

     }

     // If cluster was not matched, check also conversion sample to be sure nothing from there is missing
     for (Int_t c = 0; c < fReaderGammas->GetEntriesFast(); c++)
     {
        AliAODConversionPhoton* photon = (AliAODConversionPhoton*) fReaderGammas->At(c);
        if(!photon) continue;
        if(!((AliConversionPhotonCuts*)fConvCuts)->PhotonIsSelected(photon,fInputEvent)) continue;

        for (Int_t iElec = 0;iElec < 2;iElec++){
          Int_t tracklabel = photon->GetLabel(iElec);
          AliAODTrack *convtrack = static_cast<AliAODTrack*>(tracks->At(tracklabel));
          if(!convtrack) continue;
          if(convtrack->IsHybridGlobalConstrainedGlobal()) continue; // that means we already treated it
          if(convtrack->Pt()<0.5) continue;

          Double_t xyz[3] = {0}, pxpypz[3] = {0}, cv[21] = {0};
          convtrack->GetPxPyPz(pxpypz);
          convtrack->GetXYZ(xyz);
          convtrack->GetCovarianceXYZPxPyPz(cv);

          trackParam = new AliExternalTrackParam(xyz,pxpypz,cv,convtrack->Charge());

          AliExternalTrackParam emcParam(*trackParam);
          Float_t eta, phi, pt;
          //propagate tracks to emc surfaces
          if (!AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(&emcParam, 440., 0.139, 20., eta, phi, pt)) {
            delete trackParam;
            continue;
          }
          if( TMath::Abs(eta) > 0.75 ) {
            delete trackParam;
            continue;
          }
          // Save some time and memory in case of no DCal present
          if( nModules < 13 && ( phi < 70*TMath::DegToRad() || phi > 190*TMath::DegToRad())){
            delete trackParam;
            continue;
          }
          // Save some time and memory in case of run2
          if( nModules > 12 ){
            if (( phi < 70*TMath::DegToRad() || phi > 190*TMath::DegToRad()) && ( phi < 250*TMath::DegToRad() || phi > 340*TMath::DegToRad())){
              delete trackParam;
              continue;
            }
          }
          Float_t dEta=-999, dPhi=-999;
          Double_t trkPos[3] = {0.,0.,0.};
          if (!emcParam.GetXYZ(trkPos)){
            delete trackParam;
            continue;
          }

          AliExternalTrackParam trackParamTmp(emcParam);//Retrieve the starting point every time before the extrapolation
          if(!AliEMCALRecoUtils::ExtrapolateTrackToCluster(&trackParamTmp, clus, 0.139, 5., dEta, dPhi)) continue;
          if(TMath::Abs(dEta) > (fMatchingParamsEta[0] + pow(convtrack->Pt() + fMatchingParamsEta[1],fMatchingParamsEta[2]))) continue;
          if(TMath::Abs(dPhi) > (fMatchingParamsPhi[0] + pow(convtrack->Pt() + fMatchingParamsPhi[1],fMatchingParamsPhi[2]))) continue;
          if((clus->E()/convtrack->P()) > fMatchingEOverP) continue;
          if(highestMatchIndex == -1){ // this is the first match
            highestMatchIndex = tracklabel;
          } else{
            if(convtrack->P()>((AliAODTrack*)tracks->At(highestMatchIndex))->P()){
              highestMatchIndex = tracklabel;
            }
          }
          delete trackParam;
        }
     }
     return highestMatchIndex;
}

// Charged isolation for conversion photons
//_____________________________________________________________________________
vector<Double32_t> AliAnalysisTaskGammaIsoTree::ProcessChargedIsolation(AliAODConversionPhoton* photon){
    TLorentzVector* v4photon = new TLorentzVector();

    vector<Double32_t> vecIso; // (radius, energy)
    for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
    {
        vecIso.push_back(0.);
    }


    v4photon->SetPxPyPzE(photon->Px(),photon->Py(),photon->Pz(),photon->E());
    for (Int_t t = 0; t < fInputEvent->GetNumberOfTracks(); t++)
    {
        AliAODTrack *aodt = static_cast<AliAODTrack*>(fInputEvent->GetTrack(t));
        if(!aodt) continue;
        if(!TrackIsSelectedAOD(aodt)) continue;
        TLorentzVector v4track;
        v4track.SetPxPyPzE(aodt->Px(),aodt->Py(),aodt->Pz(),aodt->E());

        Double_t trackEta = v4track.Eta();
        Double_t trackPhi = v4track.Phi();
        if (trackPhi < 0) trackPhi += 2*TMath::Pi();

        Double_t photonEta = v4photon->Eta();
        Double_t photonPhi = v4photon->Phi();
        if (photonPhi < 0) photonPhi += 2*TMath::Pi();

        Double_t dEta = trackEta - photonEta;
        Double_t dPhi = trackPhi - photonPhi;
        Double_t dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);

        //if((dR > fTrackIsolationR[0]) && (dR > fTrackIsolationR[1])) continue;

        // track is in cone
        // check if track comes from pi0
         Bool_t trackIsFromV0 = kFALSE;
        // check that track is not from conversion
        for (Int_t iElec = 0;iElec < 2;iElec++){
          Int_t tracklabel = photon->GetLabel(iElec);
          AliAODTrack *convtrack = static_cast<AliAODTrack*>(fInputEvent->GetTrack(tracklabel));
          if(!convtrack) continue;
          if(IsSameTrack(convtrack->GetID(),aodt->GetID())){
            trackIsFromV0 = kTRUE;
          }
        }
        if(trackIsFromV0) continue;
        for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
        {
           if(dR <= fTrackIsolationR.at(r)) vecIso.at(r) += v4track.Pt();
        }
    }
    delete v4photon;
    return vecIso;
}

// charged isolation for clusters
//_____________________________________________________________________________
isoValues AliAnalysisTaskGammaIsoTree::ProcessChargedIsolation(AliAODCaloCluster* cluster){
    Double_t vertex[3] = {0,0,0};
    InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

    TLorentzVector v4cluster;
    cluster->GetMomentum(v4cluster,vertex);

    Double_t clusterEta = v4cluster.Eta();
    Double_t clusterPhi = v4cluster.Phi();
    if (clusterPhi < 0) clusterPhi += 2*TMath::Pi();
    isoValues isoV;
    vector<Double32_t> vecIso;
    for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
    {
        isoV.isolationCone.push_back(0.);
        isoV.backgroundLeft.push_back(0.);
        isoV.backgroundRight.push_back(0.);
        isoV.backgroundBack.push_back(0.);
    }
    TLorentzVector v4track;
    for (Int_t t = 0; t < fInputEvent->GetNumberOfTracks(); t++)
    {
        AliAODTrack *aodt = static_cast<AliAODTrack*>(fInputEvent->GetTrack(t));
        if(!aodt) continue;
        if(!TrackIsSelectedAOD(aodt)) continue;
        v4track.SetPxPyPzE(aodt->Px(),aodt->Py(),aodt->Pz(),aodt->E());
        Double_t trackEta = v4track.Eta();
        Double_t trackPhi = v4track.Phi();
        if (trackPhi < 0) trackPhi += 2*TMath::Pi();

        Double_t dEta = trackEta - clusterEta;
        Double_t dPhi = trackPhi - clusterPhi;
        Double_t dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);
        for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
        {
           if(dR <= fTrackIsolationR.at(r)) isoV.isolationCone.at(r) += v4track.Pt();
        }

        // do background estimates
        // left (positive) 90 degrees
        Double_t clusterPhiMod = clusterPhi + (TMath::Pi()/2);
        if(clusterPhiMod> 2*TMath::Pi()) clusterPhiMod-=2 * TMath::Pi();
        Double_t dPhiMod = trackPhi - clusterPhiMod;
        Double_t dRMod = TMath::Sqrt(dEta*dEta + dPhiMod*dPhiMod);
        for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
        {
           if(dRMod <= fTrackIsolationR.at(r)) isoV.backgroundLeft.at(r) += v4track.Pt();
        }

        // right (negative) 90 degrees
        clusterPhiMod = clusterPhi - (TMath::Pi()/2);
        if(clusterPhiMod< 0) clusterPhiMod+=2 * TMath::Pi();
        dPhiMod = trackPhi - clusterPhiMod;
        dRMod = TMath::Sqrt(dEta*dEta + dPhiMod*dPhiMod);
        for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
        {
           if(dRMod <= fTrackIsolationR.at(r)) isoV.backgroundRight.at(r) += v4track.Pt();
        }

        // backwards 180 degrees
        clusterPhiMod = clusterPhi - TMath::Pi();
        if(clusterPhiMod< 0) clusterPhiMod+= 2*TMath::Pi();
        dPhiMod = trackPhi - clusterPhiMod;
        dRMod = TMath::Sqrt(dEta*dEta + dPhiMod*dPhiMod);
        for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
        {
           if(dRMod <= fTrackIsolationR.at(r)) isoV.backgroundBack.at(r) += v4track.Pt();
        }

    }
    if(fUseHistograms) fHistoChargedIso->Fill(isoV.isolationCone.at(0)); // debug only
    return isoV;
}

//_____________________________________________________________________________
vector<Double32_t> AliAnalysisTaskGammaIsoTree::ProcessNeutralIsolation(AliAODConversionPhoton* photon){
    Double_t vertex[3] = {0,0,0};
    InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
    vector<Double32_t> vecIso;
    TLorentzVector* v4photon = new TLorentzVector();
    for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
    {
        vecIso.push_back(0.);
    }

    v4photon->SetPxPyPzE(photon->Px(),photon->Py(),photon->Pz(),photon->E());
    for (Int_t c = 0; c < fClusterEMCalCandidatesIsolation->GetEntries(); c++)
    {
        AliAODCaloCluster* clusterE = (AliAODCaloCluster*) fClusterEMCalCandidatesIsolation->At(c);
        if(!clusterE) continue;

        // AliExtraClusterInfoHelper* clusInfo = (AliExtraClusterInfoHelper*) fExtraClusterInfoBackground->At(c);
        // check if cluster is neutral
        // if(fDoBackgroundTrackMatching){
        //    if(clusInfo->isMatched()) continue;
        // }


        TLorentzVector v4cluster;
        clusterE->GetMomentum(v4cluster,vertex);
        Double_t photonEta = v4photon->Eta();
        Double_t photonPhi = v4photon->Phi();
        if (photonPhi < 0) photonPhi += 2*TMath::Pi();

        Double_t clusterEta = v4cluster.Eta();
        Double_t clusterPhi = v4cluster.Phi();
        if (clusterPhi < 0) clusterPhi += 2*TMath::Pi();

        Double_t dEta = photonEta - clusterEta;
        Double_t dPhi = photonPhi - clusterPhi;
        Double_t dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);

        for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
        {
           if(dR <= fNeutralIsolationR.at(r)) vecIso.at(r) += v4cluster.Et();
        }
    }
    delete v4photon;
    return vecIso;
}

//____Experimental isolation using EMC cells
vector<Double32_t> AliAnalysisTaskGammaIsoTree::ProcessCellIsolation(AliAODConversionPhoton* photon){
    Double_t vertex[3] = {0,0,0};
    InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

    TLorentzVector* v4photon = new TLorentzVector();
    v4photon->SetPxPyPzE(photon->Px(),photon->Py(),photon->Pz(),photon->E());
    vector<Double32_t> vecIso;
    for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
    {
        vecIso.push_back(0.);
    }

    AliVCaloCells *cells = InputEvent()->GetEMCALCells();
    const Short_t nCells = cells->GetNumberOfCells();

    // count cells above threshold per sm
    Int_t bunchCrossNo = InputEvent()->GetBunchCrossNumber();
    for(Int_t iCell=0; iCell<nCells; ++iCell) {

              // Define necessary variables
        Short_t cellId                   = 0;
        Double_t cellE = 0,  cellTime = 0, cellEFrac = 0;
        Int_t cellMCLabel = 0;

        // Get Cell
        cells->GetCell(iCell,cellId,cellE,cellTime,cellMCLabel,cellEFrac);


        UShort_t cellMax[]  = {(UShort_t) cellId};
        Bool_t   badCell    = GetCaloUtils()->GetEMCALRecoUtils()->ClusterContainsBadChannel(GetCaloUtils()->GetEMCALGeometry(),cellMax,1);
        if(badCell) continue;
        Bool_t   exoticCell = GetCaloUtils()->GetEMCALRecoUtils()->IsExoticCell(cellId,cells,bunchCrossNo);
        if(exoticCell) continue;
        // Int_t sm       = cellId / (24*48);


        // energy cut
        if(cellE < 0.1) continue;

        Float_t cellEta = 0;
        Float_t cellPhi = 0;

        fGeomEMCAL->EtaPhiFromIndex(cellId,cellEta,cellPhi);

        Double_t photonEta = v4photon->Eta();
        Double_t photonPhi = v4photon->Phi();
        if (photonPhi < 0) photonPhi += 2*TMath::Pi();

        if (cellPhi < 0) cellPhi += 2*TMath::Pi();

        Double_t dEta = photonEta - cellEta;
        Double_t dPhi = photonPhi - cellPhi;
        Double_t dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);

        for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
        {
           if(dR <= fNeutralIsolationR.at(r)) vecIso.at(r) += cellE;
        }
    }
    // // Go through all conversion and subtract E of leg if cluster was found
    // fReaderGammas    = fV0Reader->GetReconstructedGammas();

    // for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){
    //   AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(i);
    //   if(!PhotonCandidate) continue;

    //   Bool_t trackIsFromV0 = kFALSE;
    //   for (Int_t iElec = 0; iElec < 2; iElec++)
    //   {
    //     Int_t tracklabel = photon->GetTrackLabel(iElec);
    //     AliAODTrack *convt = static_cast<AliAODTrack*>(fInputEvent->GetTrack(tracklabel));
    //     //cout << "TrackID = " << aodt->GetID() <<"ConvID =" <<convt->GetID() << endl;

    //     // Loop over all clusters (no cuts)
    //     for (Int_t i = 0; i < count; i++)
    //     {
    //       /* code */
    //     }

    //   }
    // }

    delete v4photon;
    return vecIso;
}

//____Experimental isolation using EMC cells
vector<Double32_t> AliAnalysisTaskGammaIsoTree::ProcessCellIsolation(AliAODCaloCluster* cluster){
    Double_t vertex[3] = {0,0,0};
    InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

    TLorentzVector tmp;
    cluster->GetMomentum(tmp,vertex);

    TLorentzVector* v4thiscluster = new TLorentzVector(tmp);

    const Short_t nClusterCells = cluster->GetNCells();
    UShort_t* idClusterCells = cluster->GetCellsAbsId();

    vector<Double32_t> vecIso;
    for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
    {
        vecIso.push_back(0.);
    }

    AliVCaloCells *cells = InputEvent()->GetEMCALCells();
    const Short_t nCells = cells->GetNumberOfCells();

    // count cells above threshold per sm
    Int_t bunchCrossNo = InputEvent()->GetBunchCrossNumber();


    for(Int_t iCell=0; iCell<nCells; ++iCell) {

        Short_t cellId                   = 0;
        Double_t cellE = 0,  cellTime = 0, cellEFrac = 0;
        Int_t cellMCLabel = 0;

        // Get Cell
        cells->GetCell(iCell,cellId,cellE,cellTime,cellMCLabel,cellEFrac);

        UShort_t cellMax[]  = {(UShort_t) cellId};
        Bool_t   badCell    = GetCaloUtils()->GetEMCALRecoUtils()->ClusterContainsBadChannel(GetCaloUtils()->GetEMCALGeometry(),cellMax,1);
        if(badCell) continue;
        Bool_t   exoticCell = GetCaloUtils()->GetEMCALRecoUtils()->IsExoticCell(cellId,cells,bunchCrossNo);
        if(exoticCell) continue;
        // Int_t sm       = cellId / (24*48);

        // energy cut
        if(cellE < 0.1) continue;

        // check that cell is not contained in cluster
        Bool_t cellInCluster = kFALSE;
        for (Short_t c = 0; c < nClusterCells; c++)
        {
            if(idClusterCells[c] == cellId) cellInCluster = kTRUE;
        }

        if(cellInCluster) continue;

        Float_t cellEta = 0.;
        Float_t cellPhi = 0.;
        fGeomEMCAL->EtaPhiFromIndex(cellId,cellEta,cellPhi);

        Double_t clusterEta = v4thiscluster->Eta();
        Double_t clusterPhi = v4thiscluster->Phi();
        if (clusterPhi < 0) clusterPhi += 2*TMath::Pi();

        if (cellPhi < 0) cellPhi += 2*TMath::Pi();

        Double_t dEta = clusterEta - cellEta;
        Double_t dPhi = clusterPhi - cellPhi;
        Double_t dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);

        for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
        {
           if(dR <= fNeutralIsolationR.at(r)) vecIso.at(r) += cellE;
        }
    }
    delete v4thiscluster;
    return vecIso;
}
//_____________________________________________________________________________
vector<Double32_t> AliAnalysisTaskGammaIsoTree::ProcessNeutralIsolation(AliAODCaloCluster* cluster){
    Double_t vertex[3] = {0,0,0};
    InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
    vector<Double32_t> vecIso;
    for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
    {
        vecIso.push_back(0.);
    }

    TLorentzVector tmp;
    cluster->GetMomentum(tmp,vertex);

    TLorentzVector* v4thiscluster = new TLorentzVector(tmp);
    for (Int_t c = 0; c < fClusterEMCalCandidatesIsolation->GetEntries(); c++)
    {
        AliAODCaloCluster* clusterE = (AliAODCaloCluster*) fClusterEMCalCandidatesIsolation->At(c);
        if(!clusterE) continue;
        if(clusterE->GetID() == cluster->GetID()) continue;
        // AliExtraClusterInfoHelper* clusInfo = (AliExtraClusterInfoHelper*) fExtraClusterInfoBackground->At(c);
        // if(fDoBackgroundTrackMatching){
        //    if(clusInfo->isMatched()) continue;
        // }

        TLorentzVector v4othercluster;
        clusterE->GetMomentum(v4othercluster,vertex);

        Double_t thisclusterEta = v4thiscluster->Eta();
        Double_t thisclusterPhi = v4thiscluster->Phi();
        if (thisclusterPhi < 0) thisclusterPhi += 2*TMath::Pi();

        Double_t otherclusterEta = v4othercluster.Eta();
        Double_t otherclusterPhi = v4othercluster.Phi();
        if (otherclusterPhi < 0) otherclusterPhi += 2*TMath::Pi();

        Double_t dEta = thisclusterEta - otherclusterEta;
        Double_t dPhi = thisclusterPhi - otherclusterPhi;
        Double_t dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);

        for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
        {
           if(dR <= fNeutralIsolationR.at(r)) vecIso.at(r) += v4othercluster.Et();
        }
    }
    delete v4thiscluster;
    return vecIso;
}
//_____________________________________________________________________________
isoValues AliAnalysisTaskGammaIsoTree::ProcessMCIsolation(Int_t mclabel){
  vector<Double32_t> vecIso;
  isoValues isoV;
  if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
  {
      isoV.isolationCone.push_back(0.);
      isoV.backgroundLeft.push_back(0.);
      isoV.backgroundRight.push_back(0.);
      isoV.backgroundBack.push_back(0.);
  }
  for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
  {
      isoV.isolationCone.push_back(0.);
      isoV.backgroundLeft.push_back(0.);
      isoV.backgroundRight.push_back(0.);
      isoV.backgroundBack.push_back(0.);
  }
  AliAODMCParticle *thisParticle = (AliAODMCParticle*) fAODMCTrackArray->At(mclabel);
  Double_t thisEta = thisParticle->Eta();
  Double_t thisPhi = thisParticle->Phi();

  // Get potential conv daughters to exclude them from iso in cone);
  Int_t nDaughters = thisParticle->GetNDaughters();
  Int_t daughterLabels[2];
  Int_t nFound = 0;
  AliAODMCParticle *daughter;
  if(thisParticle->GetPdgCode()==22 && (nDaughters==2)){
    for (Int_t d = 0; d < nDaughters; d++)
    {
      
        Int_t tmp = thisParticle->GetDaughterLabel(d);
        daughter = static_cast<AliAODMCParticle *>(fAODMCTrackArray->At(tmp));
        if(!daughter) continue;
        if(TMath::Abs(daughter->GetPdgCode())==11){
          daughterLabels[nFound] = tmp;
          nFound++;
        }
    }
  }

  AliAODMCParticle *pmc = NULL;
  for (Int_t p = 0; p < fAODMCTrackArray->GetEntriesFast(); p++)
  {
      if (p == mclabel)
          continue;
      
      pmc = static_cast<AliAODMCParticle *>(fAODMCTrackArray->At(p));

      Bool_t isPrimary = fEventCuts->IsConversionPrimaryAOD(fInputEvent, pmc, mcProdVtxX, mcProdVtxY, mcProdVtxZ); 
      if(!isPrimary) continue;
   
      // found conversion, exclude e daughters from iso
      if(nFound == 2){
           if((p==daughterLabels[0]) || (p==daughterLabels[1])) continue;
      }

      Double_t otherEta = pmc->Eta();
      Double_t otherPhi = pmc->Phi();

      if (thisPhi < 0)
        thisPhi += TMath::TwoPi();
      if (otherPhi < 0)
        otherPhi += TMath::TwoPi();
      Double_t dEta = otherEta - thisEta;
      Double_t dPhi = otherPhi - thisPhi;
      Double_t dR = TMath::Sqrt((dEta * dEta) + (dPhi * dPhi));

      if (pmc->Charge() != 0){ //charged
        if(pmc->Pt() > fPtCut){
          for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
          {
            if(dR<fTrackIsolationR.at(r)) isoV.isolationCone.at(r)+= pmc->Pt();
          }
        }
      }
      if ((pmc->Charge() == 0) && (IsInEMCalAcceptance(pmc))){ //neutral

        for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
        {
           if(dR<fNeutralIsolationR.at(r)) isoV.isolationCone.at(r+fTrackIsolationR.size()) += pmc->Pt();
        }
      }
      
      // do background estimates only for charged
      if(pmc->Charge() == 0) continue;
      if(pmc->Pt() < fPtCut) continue;

      // do background estimates
      // left (positive) 90 degrees
      Double_t thisPhiMod = thisPhi + (TMath::Pi()/2);
      if(thisPhiMod> 2*TMath::Pi()) thisPhiMod-=2 * TMath::Pi();
      Double_t dPhiMod = otherPhi - thisPhiMod;
      Double_t dRMod = TMath::Sqrt(dEta*dEta + dPhiMod*dPhiMod);
      for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
      {
          if(dRMod <= fTrackIsolationR.at(r)) isoV.backgroundLeft.at(r) += pmc->Pt();
      }

      // right (negative) 90 degrees
      thisPhiMod = thisPhi - (TMath::Pi()/2);
      if(thisPhiMod< 0) thisPhiMod+=2 * TMath::Pi();
      dPhiMod = otherPhi - thisPhiMod;
      dRMod = TMath::Sqrt(dEta*dEta + dPhiMod*dPhiMod);
      for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
      {
          if(dRMod <= fTrackIsolationR.at(r)) isoV.backgroundRight.at(r) += pmc->Pt();
      }

      // backwards 180 degrees
      thisPhiMod = thisPhi - TMath::Pi();
      if(thisPhiMod< 0) thisPhiMod+= 2*TMath::Pi();
      dPhiMod = otherPhi - thisPhiMod;
      dRMod = TMath::Sqrt(dEta*dEta + dPhiMod*dPhiMod);
      for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
      {
          if(dRMod <= fTrackIsolationR.at(r)) isoV.backgroundBack.at(r) += pmc->Pt();
      }


  }
  return isoV;
}

//_____________________________________________________________________________
Int_t AliAnalysisTaskGammaIsoTree::ProcessTagging(AliAODConversionPhoton* photon){
  Int_t taggedConv = 0;
  Int_t taggedClus = 0;
  AliAODConversionMother *pi0cand = NULL;
  fReaderGammas    = fV0Reader->GetReconstructedGammas();
  for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){
    AliAODConversionPhoton* otherPhoton = (AliAODConversionPhoton*) fReaderGammas->At(i);
    if(!otherPhoton) continue;
    if(!((AliConversionPhotonCuts*)fConvCuts)->PhotonIsSelected(otherPhoton,fInputEvent)) continue;

    pi0cand = new AliAODConversionMother(photon,otherPhoton);

    // check mass window
    Double_t mass = pi0cand->M();
    fHistoTaggingPCMPCM->Fill(mass,photon->Pt(),fWeightJetJetMC);
    if((mass > fPi0TaggingWindow[0]) && (mass < fPi0TaggingWindow[1])){ // pi0
        taggedConv = 1;
        delete pi0cand;
        break;
    } else if((mass > fEtaTaggingWindow[0]) && (mass < fEtaTaggingWindow[1])){ // eta
        taggedConv = 1;
        delete pi0cand;
        break;
    }
    delete pi0cand;
  }

  Double_t vertex[3] = {0,0,0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  for (Int_t c = 0; c < fClusterEMCalCandidatesTagging->GetEntries(); c++)
  {
    // AliExtraClusterInfoHelper* clusInfo = (AliExtraClusterInfoHelper*) fExtraClusterInfoBackground->At(c);
    // if(fDoBackgroundTrackMatching){
      //  if(clusInfo->isMatched()) continue;
    // }
    // TLorentzvector with cluster
    TLorentzVector clusterVector;
    ((AliAODCaloCluster*)fClusterEMCalCandidatesTagging->At(c))->GetMomentum(clusterVector,vertex);

    TLorentzVector* tmpvec = new TLorentzVector();
    tmpvec->SetPxPyPzE(clusterVector.Px(),clusterVector.Py(),clusterVector.Pz(),clusterVector.E());

    // convert to AODConversionPhoton
    AliAODConversionPhoton *PhotonCluster = new AliAODConversionPhoton(tmpvec);
    if(!PhotonCluster){ delete tmpvec; continue;}

    pi0cand = new AliAODConversionMother(photon,PhotonCluster);

    // check mass window
    Double_t mass = pi0cand->M();
    fHistoTaggingPCMEMC->Fill(mass,photon->Pt(),fWeightJetJetMC);
   // cout << "checking if mass " << mass <<" is between " <<  fPi0TaggingWindow[0] << " and " << fPi0TaggingWindow[1] << endl;
    if((mass > fPi0TaggingWindow[0]) && (mass < fPi0TaggingWindow[1])){ // pi0
        taggedClus = 2;
        delete pi0cand;
        delete PhotonCluster;
        delete tmpvec;
        break;
    } else  if((mass > fEtaTaggingWindow[0]) && (mass < fEtaTaggingWindow[1])){ // eta
        delete pi0cand;
        delete PhotonCluster;
        delete tmpvec;
        taggedClus = 2;
        break;
    }
    delete pi0cand;
    delete PhotonCluster;
    delete tmpvec;
  }
  //cout << "tagging = " << taggedConv+taggedClus << endl;
  return (taggedConv+taggedClus);
}

// tagging of calo clusters
Int_t AliAnalysisTaskGammaIsoTree::ProcessTagging(AliAODCaloCluster* cluster){
  Int_t taggedConv = 0;
  Int_t taggedClus = 0;

  Double_t vertex[3] = {0,0,0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  TLorentzVector tmp;
  cluster->GetMomentum(tmp,vertex);

  TLorentzVector* thisclustervec = new TLorentzVector();
  thisclustervec->SetPxPyPzE(tmp.Px(),tmp.Py(),tmp.Pz(),tmp.E());

  // convert to AODConversionPhoton
  AliAODConversionPhoton *ThisPhotonCluster = new AliAODConversionPhoton(thisclustervec);
  if(!ThisPhotonCluster){ delete thisclustervec; return -1;}

  AliAODConversionMother *pi0cand = NULL;
  fReaderGammas    = fV0Reader->GetReconstructedGammas();

  for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){
    AliAODConversionPhoton* otherPhoton = (AliAODConversionPhoton*) fReaderGammas->At(i);
    if(!otherPhoton) continue;
    if(!((AliConversionPhotonCuts*)fConvCuts)->PhotonIsSelected(otherPhoton,fInputEvent)) continue;

    pi0cand = new AliAODConversionMother(ThisPhotonCluster,otherPhoton);

    // check mass window
    Double_t mass = pi0cand->M();
    fHistoTaggingEMCPCM->Fill(mass,ThisPhotonCluster->Pt(),fWeightJetJetMC);
    if((mass > fPi0TaggingWindow[0]) && (mass < fPi0TaggingWindow[1])){ // pi0
        taggedConv = 1;
        delete pi0cand;
        break;
    } else if((mass > fEtaTaggingWindow[0]) && (mass < fEtaTaggingWindow[1])){ // eta
        taggedConv = 1;
        delete pi0cand;
        break;
    }
    delete pi0cand;
  }


  for (Int_t c = 0; c < fClusterEMCalCandidatesTagging->GetEntries(); c++)
  {
    // TLorentzvector with cluster
    TLorentzVector clusterVector;

    if(((AliAODCaloCluster*)fClusterEMCalCandidatesTagging->At(c))->GetID() == cluster->GetID()) continue;

    ((AliAODCaloCluster*)fClusterEMCalCandidatesTagging->At(c))->GetMomentum(clusterVector,vertex);

    TLorentzVector* tmpvec = new TLorentzVector();
    tmpvec->SetPxPyPzE(clusterVector.Px(),clusterVector.Py(),clusterVector.Pz(),clusterVector.E());

    // convert to AODConversionPhoton
    AliAODConversionPhoton *PhotonCluster = new AliAODConversionPhoton(tmpvec);
    if(!PhotonCluster){ delete tmpvec; continue;}

    pi0cand = new AliAODConversionMother(ThisPhotonCluster,PhotonCluster);

    // check mass window
    Double_t mass = pi0cand->M();
    fHistoTaggingEMCEMC->Fill(mass,ThisPhotonCluster->Pt(),fWeightJetJetMC);
    if((mass > fPi0TaggingWindow[0]) && (mass < fPi0TaggingWindow[1])){ // pi0
        taggedClus = 2;
        delete pi0cand;
        delete PhotonCluster;
        delete tmpvec;
        break;
    } else  if((mass > fEtaTaggingWindow[0]) && (mass < fEtaTaggingWindow[1])){ // eta
        delete pi0cand;
        delete PhotonCluster;
        delete tmpvec;
        taggedClus = 2;
        break;
    }
    delete pi0cand;
    delete PhotonCluster;
    delete tmpvec;
  }
  delete ThisPhotonCluster;
  delete thisclustervec;
  return taggedConv+taggedClus;
}

// delete covariance matrix etc. to save space
void AliAnalysisTaskGammaIsoTree::ReduceTrackInfo(){
    for (Int_t i = 0; i < fTracks->GetEntries(); i++)
    {
      // fTracks.at(i)->SetTOFchi2(0);
      // fTracks.at(i)->SetTOFLabel(NULL);
      // fTracks.at(i)->SetTOFsignalDx(0);
      // fTracks.at(i)->SetTOFsignalDz(0);
      ((AliAODTrack*)fTracks->At(i))->Clear();
    }
}

//________________________________________________________________________
void AliAnalysisTaskGammaIsoTree::RelabelAODPhotonCandidates(Bool_t mode){

  // Relabeling For AOD Event
  // ESDiD -> AODiD
  // MCLabel -> AODMCLabel

  Int_t* fMCEventPos = nullptr;
  Int_t* fMCEventNeg = nullptr;
  Int_t* fESDArrayPos = nullptr;
  Int_t* fESDArrayNeg = nullptr;
  if(mode){
    fMCEventPos = new Int_t[fReaderGammas->GetEntries()];
    fMCEventNeg = new Int_t[fReaderGammas->GetEntries()];
    fESDArrayPos = new Int_t[fReaderGammas->GetEntries()];
    fESDArrayNeg = new Int_t[fReaderGammas->GetEntries()];
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


  if(!mode){
    delete[] fMCEventPos;
    delete[] fMCEventNeg;
    delete[] fESDArrayPos;
    delete[] fESDArrayNeg;
  }
}

void AliAnalysisTaskGammaIsoTree::FillConversionHistos(AliAODConversionPhoton* photon,vector<Double32_t> isoCharged,vector<Double32_t> isoNeutral,vector<Double32_t> isoCell,Int_t tmptag){
    fConvPtBeforeAcc->Fill(photon->Pt(),fWeightJetJetMC);
    if (!IsInEMCalAcceptance(photon)) return;

    fConvPt->Fill(photon->Pt(),fWeightJetJetMC);
    if(tmptag<2 ){
       fConvPtTaggedCalo->Fill(photon->Pt(),fWeightJetJetMC);
    } else{
       fConvPtTaggedAsDecayCalo->Fill(photon->Pt(),fWeightJetJetMC);
    }

    for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
    {
      fConvIsoCharged[r]->Fill(isoCharged.at(r),photon->Pt(),fWeightJetJetMC);
      for (UInt_t e = 0; e < fTrackIsolationE.size(); e++)
      {
        if(isoCharged.at(r) < fTrackIsolationE.at(e)){
          fConvPtIsoCharged[r][e]->Fill(photon->Pt(),fWeightJetJetMC);
          if(tmptag<2 )fConvPtTaggedCaloIsoCharged[r][e]->Fill(photon->Pt(),fWeightJetJetMC); // not tagged by calo
        }
      }

    }

    for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
    {
      fConvIsoNeutral[r]->Fill(isoNeutral.at(r),photon->Pt(),fWeightJetJetMC);
      fConvIsoFull[r]->Fill(isoNeutral.at(r) + isoCharged.at(r),photon->Pt(),fWeightJetJetMC);
      fConvIsoCell[r]->Fill(isoCell.at(r),photon->Pt(),fWeightJetJetMC);
      for (UInt_t e = 0; e < fNeutralIsolationE.size(); e++)
      {
        if(isoNeutral.at(r) < fNeutralIsolationE.at(e))fConvPtIsoNeutral[r][e]->Fill(photon->Pt(),fWeightJetJetMC);
        if(isoCell.at(r) < fNeutralIsolationE.at(e)) fConvPtIsoCell[r][e]->Fill(photon->Pt(),fWeightJetJetMC);
        if((isoNeutral.at(r) + isoCharged.at(r)) < fNeutralIsolationE.at(e)) fConvPtIsoFull[r][e]->Fill(photon->Pt(),fWeightJetJetMC);
        if(tmptag<2){
            if(isoNeutral.at(r) < fNeutralIsolationE.at(e))fConvPtTaggedCaloIsoNeutral[r][e]->Fill(photon->Pt(),fWeightJetJetMC);
            if(isoCell.at(r) < fNeutralIsolationE.at(e)) fConvPtTaggedCaloIsoCell[r][e]->Fill(photon->Pt(),fWeightJetJetMC);
            if((isoNeutral.at(r) + isoCharged.at(r)) < fNeutralIsolationE.at(e)) fConvPtTaggedCaloIsoFull[r][e]->Fill(photon->Pt(),fWeightJetJetMC);
        }
      }

    }


    //
    // ─── FILL INV MASS HISTOS ────────────────────────────────────────
    //

    Double_t vertex[3] = {0,0,0};
    InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

    Bool_t isConv = kFALSE;
    Bool_t isDecay = kFALSE;
    Bool_t isPrompt = kFALSE;
    Bool_t isFrag = kFALSE;
    Bool_t isSignal = kFALSE;

    if(fIsMC>0){
      isConv = IsTrueConversionPhoton(photon);
      if(isConv) {
        isDecay = IsDecayPhoton(photon);
        isPrompt = IsPromptPhoton(photon);
        isFrag = IsFragPhoton(photon);
        isSignal = isPrompt || isFrag;
      }
    }

    for (Int_t c = 0; c < fClusterEMCalCandidates->GetEntries(); c++)
    {
      // TLorentzvector with cluster
      TLorentzVector clusterVector;
      AliAODCaloCluster* clus = (AliAODCaloCluster*)fClusterEMCalCandidates->At(c);
      if(!clus) continue;
      clus->GetMomentum(clusterVector,vertex);
      if((clus->GetM02() < fMinM02) || (clus->GetM02() > fMaxM02)) continue;
      TLorentzVector* tmpvec = new TLorentzVector();
      tmpvec->SetPxPyPzE(clusterVector.Px(),clusterVector.Py(),clusterVector.Pz(),clusterVector.E());

      // convert to AODConversionPhoton
      AliAODConversionPhoton *PhotonCluster = new AliAODConversionPhoton(tmpvec);
      if(!PhotonCluster){ delete tmpvec; continue;}

      AliAODConversionMother* pi0cand = new AliAODConversionMother(photon,PhotonCluster);

      // check mass window
      Double_t mass = pi0cand->M();

      fConvInvMass->Fill(mass,photon->Pt(),fWeightJetJetMC);

      if(fIsMC>0){
          if(isConv){
            fConvTrueInvMass->Fill(mass,photon->Pt(),fWeightJetJetMC);
            if(isDecay){
              fConvTrueInvMass_FromDecay->Fill(mass,photon->Pt(),fWeightJetJetMC);
            } else if(isSignal){
              fConvTrueInvMass_FromDirect->Fill(mass,photon->Pt(),fWeightJetJetMC);
            }
          }
      }

      for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
      {
        if((isoCharged.at(r) > fAntiIsolationE[0]) &&(isoCharged.at(r) < fAntiIsolationE[1]) ) fConvInvMassAntiChargedIsolated[r]->Fill(mass,photon->Pt(),fWeightJetJetMC);
        for (UInt_t e = 0; e < fTrackIsolationE.size(); e++)
        {
          if(isoCharged.at(r) < fTrackIsolationE.at(e)){
            fConvInvMassChargedIsolated[r][e]->Fill(mass,photon->Pt(),fWeightJetJetMC);
          }
        }
      }

      for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
      {
        if((isoNeutral.at(r) > fAntiIsolationE[0]) && (isoNeutral.at(r) < fAntiIsolationE[1])) fConvInvMassAntiNeutralIsolated[r]->Fill(mass,photon->Pt(),fWeightJetJetMC);
        if((isoCell.at(r) > fAntiIsolationE[0]) && (isoCell.at(r) < fAntiIsolationE[1])) fConvInvMassAntiCellIsolated[r]->Fill(mass,photon->Pt(),fWeightJetJetMC);
        if(((isoNeutral.at(r) + isoCharged.at(r)) > fAntiIsolationE[0]) && ((isoNeutral.at(r) + isoCharged.at(r)) < fAntiIsolationE[1])) fConvInvMassAntiFullIsolated[r]->Fill(mass,photon->Pt(),fWeightJetJetMC);
        for (UInt_t e = 0; e < fNeutralIsolationE.size(); e++)
        {
          if(isoNeutral.at(r) < fNeutralIsolationE.at(e)) fConvInvMassNeutralIsolated[r][e]->Fill(mass,photon->Pt(),fWeightJetJetMC);
          if(isoCell.at(r) < fNeutralIsolationE.at(e)) fConvInvMassCellIsolated[r][e]->Fill(mass,photon->Pt(),fWeightJetJetMC);
          if((isoNeutral.at(r) + isoCharged.at(r)) < fNeutralIsolationE.at(e)) fConvInvMassFullIsolated[r][e]->Fill(mass,photon->Pt(),fWeightJetJetMC);

        }
      }

      delete tmpvec;
      delete PhotonCluster;
      delete pi0cand;
  }
}
void AliAnalysisTaskGammaIsoTree::FillCaloHistosPurity(AliAODCaloCluster* clus,AliAODConversionPhoton* photon,vector<Double32_t> isoCharged,vector<Double32_t> isoNeutral,vector<Double32_t> isoCell,Int_t tmptag, Double_t weight){
  if(!clus) return;

  //if (!IsInEMCalAcceptance(photon)) return;

  Bool_t isTruePhoton = kFALSE;
  Bool_t isDecay = kFALSE;
  Bool_t isPrompt = kFALSE;
  Bool_t isFrag = kFALSE;
  Bool_t isSignal = kFALSE;
  AliAODMCParticle *MCPhoton = NULL;

  Double_t vertex[3] = {0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  Bool_t truePhotonFromConv = kFALSE;
  Bool_t truePhotonNormal = kFALSE;
  // check MC properties
  isoValues mcIso;
  vector<Double32_t> mcIsoCharged;
  vector<Double32_t> mcIsoNeutral;
  vector<Double32_t> mcIsoFull;
  if(fIsMC>0){
      const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
      Double_t mcProdVtxX   = primVtxMC->GetX();
      Double_t mcProdVtxY   = primVtxMC->GetY();
      Double_t mcProdVtxZ   = primVtxMC->GetZ();
      Int_t photonlabel = -1;
      if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
      if (fAODMCTrackArray){
          if (photon->GetIsCaloPhoton() == 0) AliFatal("CaloPhotonFlag has not been set task will abort");
          if (photon->GetNCaloPhotonMCLabels()>0){
            photonlabel = photon->GetCaloPhotonMCLabel(0);
            MCPhoton = (AliAODMCParticle*) fAODMCTrackArray->At(photon->GetCaloPhotonMCLabel(0));
          }
      } else {
        AliInfo("AODMCTrackArray could not be loaded");
        return;
      }
      if (photon->IsLargestComponentPhoton() || (photon->IsLargestComponentElectron() && photon->IsConversion())) {
         Bool_t isPrimary = fEventCuts->IsConversionPrimaryAOD(fInputEvent, MCPhoton, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
         if(isPrimary){
           isTruePhoton = kTRUE;
           if(photon->IsLargestComponentElectron() && photon->IsConversion()){
             truePhotonFromConv = kTRUE;
           } else if(photon->IsLargestComponentPhoton()){
             truePhotonNormal = kTRUE;
           }
         }

         // check for proper label, in case of conversion return photon label,otherwise do nothing
         Int_t tmplabel = GetProperLabel(MCPhoton);
         if(tmplabel != -99) photonlabel = tmplabel;
         MCPhoton = (AliAODMCParticle*) fAODMCTrackArray->At(photonlabel);

         isDecay = IsDecayPhoton(photonlabel);
         isPrompt = IsPromptPhoton(photonlabel);
         isFrag = IsFragPhoton(photonlabel);
         isSignal = isPrompt || isFrag;
      }
      
      mcIso = ProcessMCIsolation(photonlabel);

      for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
      {
        Double_t isoCh = mcIso.isolationCone.at(r) - fChargedRhoTimesArea[r];
        mcIsoCharged.push_back(isoCh);
      }
      for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
      {
        mcIsoNeutral.push_back(mcIso.isolationCone.at(r+fTrackIsolationR.size()));
      }
      for (UInt_t r = 0; r < mcIsoCharged.size(); r++)
      {
        mcIsoFull.push_back(mcIsoCharged.at(r)+ mcIsoNeutral.at(r)); // they should always have same length
      }
  }

  // Caluclate M02
  Double_t m02 = clus->GetM02();

  const Int_t   nc = clus->GetNCells();
  Int_t   absCellIdList[nc];
  Float_t   maxEList[nc];

  // GetNLM
  Int_t nlm = fClusterCutsEMC->GetNumberOfLocalMaxima(clus,fInputEvent,absCellIdList,maxEList);

  // Split cluster
  AliAODCaloCluster* clusSub1 = new AliAODCaloCluster();
  AliAODCaloCluster* clusSub2 = new AliAODCaloCluster();
  // split clusters according to their shares in the cluster (NLM == 1) needs to be treated differently
  if (nlm == 1){
    Int_t absCellIdFirst    = ((AliCaloPhotonCuts*)fClusterCutsEMC)->FindLargestCellInCluster(clus, fInputEvent);
    Int_t absCellIdSecond   = ((AliCaloPhotonCuts*)fClusterCutsEMC)->FindSecondLargestCellInCluster(clus, fInputEvent);

    ((AliCaloPhotonCuts*)fClusterCutsEMC)->SplitEnergy(absCellIdFirst, absCellIdSecond, clus, fInputEvent, fIsMC, clusSub1, clusSub2);
  } else if (nlm > 1 ){
    ((AliCaloPhotonCuts*)fClusterCutsEMC)->SplitEnergy(absCellIdList[0], absCellIdList[1], clus, fInputEvent, fIsMC, clusSub1, clusSub2);
  }

  TLorentzVector clusterVector1;
  clusSub1->GetMomentum(clusterVector1,vertex);
  TLorentzVector* tmpvec1 = new TLorentzVector();
  tmpvec1->SetPxPyPzE(clusterVector1.Px(),clusterVector1.Py(),clusterVector1.Pz(),clusterVector1.E());
  // convert to AODConversionPhoton
  AliAODConversionPhoton *PhotonCandidate1=new AliAODConversionPhoton(tmpvec1);
  if(!PhotonCandidate1){
    delete clusSub1;
    delete tmpvec1;
  }
  // TLorentzvector with sub cluster 2
  TLorentzVector clusterVector2;
  clusSub2->GetMomentum(clusterVector2,vertex);
  TLorentzVector* tmpvec2 = new TLorentzVector();
  tmpvec2->SetPxPyPzE(clusterVector2.Px(),clusterVector2.Py(),clusterVector2.Pz(),clusterVector2.E());
  // convert to AODConversionPhoton
  AliAODConversionPhoton *PhotonCandidate2=new AliAODConversionPhoton(tmpvec2);
  if(!PhotonCandidate2){
    delete clusSub2;
    delete tmpvec2;
  }

  Float_t mass = -1;
  AliAODConversionMother* pi0cand = NULL;
  if(PhotonCandidate1 && PhotonCandidate2){
      pi0cand = new AliAODConversionMother(PhotonCandidate1,PhotonCandidate2);
      if(pi0cand) mass = pi0cand->M();
  }

  // Take the trash out
  if(pi0cand)          delete pi0cand;
  if(PhotonCandidate1) delete PhotonCandidate1;
  if(PhotonCandidate2) delete PhotonCandidate2;
  if(clusSub1)         delete clusSub1;
  if(clusSub2)         delete clusSub2;
  if(tmpvec1)          delete tmpvec1;
  if(tmpvec2)          delete tmpvec2;



  // Fill histos
  Double_t fillArr[3] = {m02,photon->Pt(),mass};
  fCaloM02->Fill(fillArr,weight);



  if(fIsMC>0){
    if(isTruePhoton){
      fCaloTrueM02->Fill(fillArr,weight);
      if(isDecay){
        fCaloTrueM02_FromDecay->Fill(fillArr,weight);
      } else if(isSignal){
        fCaloTrueM02_FromDirect->Fill(fillArr,weight);
        
        // for split efficiency
        fCaloTrueSignalPtClusterCuts->Fill(MCPhoton->Pt(),weight);
        fCaloTrueSignalRecPtClusterCuts->Fill(MCPhoton->Pt(),weight);
      }
    }
  }

  for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
  {
    Double_t fillArrABCD[3] = {m02,isoCharged.at(r),photon->Pt()};
    fCaloM02vsIsoVsPt[r]->Fill(fillArrABCD,weight);
    if((isoCharged.at(r) > fAntiIsolationE[0]) && (isoCharged.at(r) < fAntiIsolationE[1])) fCaloM02AntiChargedIsolated[r]->Fill(fillArr,weight);
    for (UInt_t e = 0; e < fTrackIsolationE.size(); e++)
    {
      if(isoCharged.at(r) < fTrackIsolationE.at(e)){
        fCaloM02ChargedIsolated[r][e]->Fill(fillArr,weight);

        if(isTruePhoton){
           if(isSignal){
            fCaloM02ChargedIsolated_FromDirect[r][e]->Fill(fillArr,weight);
            if(truePhotonFromConv) fCaloM02ChargedIsolated_FromDirectConvOnly[r][e]->Fill(fillArr,weight);
            if(truePhotonNormal) fCaloM02ChargedIsolated_FromDirectNormOnly[r][e]->Fill(fillArr,weight);
           }
        }
      }
    }
  }

  for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
  {
    if((isoNeutral.at(r) > fAntiIsolationE[0]) && (isoNeutral.at(r) < fAntiIsolationE[1])) fCaloM02AntiNeutralIsolated[r]->Fill(fillArr,weight);
    if((isoCell.at(r) > fAntiIsolationE[0]) && (isoCell.at(r) < fAntiIsolationE[1])) fCaloM02AntiCellIsolated[r]->Fill(fillArr,weight);
    if(((isoNeutral.at(r) + isoCharged.at(r)) > fAntiIsolationE[0]) && ((isoNeutral.at(r) + isoCharged.at(r)) < fAntiIsolationE[1])) fCaloM02AntiFullIsolated[r]->Fill(fillArr,weight);
    for (UInt_t e = 0; e < fNeutralIsolationE.size(); e++)
    {
      if(isoNeutral.at(r) < fNeutralIsolationE.at(e)) fCaloM02NeutralIsolated[r][e]->Fill(fillArr,weight);
      if(isoCell.at(r) < fNeutralIsolationE.at(e)) fCaloM02CellIsolated[r][e]->Fill(fillArr,weight);
      if((isoNeutral.at(r) + isoCharged.at(r)) < fNeutralIsolationE.at(e)) fCaloM02FullIsolated[r][e]->Fill(fillArr,weight);
    }
  }
}

void AliAnalysisTaskGammaIsoTree::FillCaloTree(AliAODCaloCluster* clus,AliAODConversionPhoton* photon,isoValues isoCharged,vector<Double32_t> isoNeutral,vector<Double32_t> isoCell,Int_t tmptag, Double_t weight){
  if(!clus) return;

  //if (!IsInEMCalAcceptance(photon)) return;

  Bool_t isTruePhoton = kFALSE;
  Bool_t isDecay = kFALSE;
  Bool_t isPrompt = kFALSE;
  Bool_t isFrag = kFALSE;
  Bool_t isSignal = kFALSE;
  AliAODMCParticle *MCPhoton = NULL;
  AliAODMCParticle *MCPhotonOld = NULL;

  Double_t vertex[3] = {0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  Bool_t truePhotonFromConv = kFALSE;
  Bool_t truePhotonNormal = kFALSE;
  // check MC properties
  isoValues mcIso;
  vector<Double32_t> mcIsoCharged;
  vector<Double32_t> mcIsoNeutral;
  vector<Double32_t> mcIsoFull;

  Int_t photonlabel = -1;
  // ────────────────────────────────────────────────────────────────────────────────
  // Process track matching info

  Bool_t isConvMatch = kFALSE; // kTrue if from conv sampe, else from primary sample
  Float_t trackdEta = 0;
  Float_t trackdPhi = 0;
  Float_t trackP = -1;
  Float_t trackPt = -1;
  
  AliAODTrack* aodt = NULL;
  if(fClusterCutsEMCTrackMatching->CheckClusterForTrackMatch(clus)){
    Int_t labelTrackClosest = -1;
    if(fClusterCutsEMCTrackMatching->GetClosestMatchedTrackToCluster(fInputEvent,clus,labelTrackClosest)){
        Int_t properLabel = labelTrackClosest;
        // if(labelTrack<0) properLabel = (-1 * labelTrack) - 1; // conversion hybrid none hybrid
        aodt = static_cast<AliAODTrack*> (fInputEvent->GetTrack(properLabel));
        if(TrackIsSelectedAOD(aodt)){
          Float_t dEta = -9;
          Float_t dPhi = -9;
          // Get dEta and dPhi of track on EMCal surface!
          ((AliCaloTrackMatcher*) fClusterCutsEMCTrackMatching->GetCaloTrackMatcherInstance())->GetTrackClusterMatchingResidual(aodt->GetID(),clus->GetID(),dEta,dPhi);
          trackdEta = dEta;
          trackdPhi = dPhi;
          trackP = aodt->P();
          trackPt = aodt->Pt();
          isConvMatch = kFALSE;
        }
    }
  } else {
    // check for conv matching
    if(!fReaderGammas) fReaderGammas    = fV0Reader->GetReconstructedGammas();
    for (Int_t conv = 0; conv < fReaderGammas->GetEntriesFast(); conv++)
    {
      AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(conv);
      for (Int_t i = 0;i < 2;i++){ // loop over daughters
         Int_t tracklabel = PhotonCandidate->GetLabel(i);
         AliVTrack *inTrack = 0x0;
         if(((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data()))->AreAODsRelabeled()){
            inTrack = static_cast<AliVTrack*>(fInputEvent->GetTrack(tracklabel));
         } else {
            for(Int_t ii=0;ii<fInputEvent->GetNumberOfTracks();ii++) {
              inTrack = static_cast<AliVTrack*>(fInputEvent->GetTrack(ii));
              if(inTrack){
                if(inTrack->GetID() == tracklabel) {
                break;
              }
            }
          }
        } 
        Float_t dEtaConv = 0;
        Float_t dPhiConv = 0;
        Bool_t propagated = ((AliCaloTrackMatcher*) fClusterCutsEMCTrackMatching->GetCaloTrackMatcherInstance())->PropagateV0TrackToClusterAndGetMatchingResidual(inTrack,clus,fInputEvent,dEtaConv,dPhiConv);
        if (propagated){
           trackdEta = dEtaConv;
           trackdPhi = dPhiConv;
           trackP = inTrack->P();
           trackPt = inTrack->Pt();
           isConvMatch = kTRUE;
        }
      } // end daughter loop
    } // end conv photon loop 
  } // end conv matching case
  
  // ────────────────────────────────────────────────────────────────────────────────

  if(fIsMC>0){
      const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
      Double_t mcProdVtxX   = primVtxMC->GetX();
      Double_t mcProdVtxY   = primVtxMC->GetY();
      Double_t mcProdVtxZ   = primVtxMC->GetZ();
      photonlabel = -1;
      if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
      if (fAODMCTrackArray){
          if (photon->GetIsCaloPhoton() == 0) AliFatal("CaloPhotonFlag has not been set task will abort");
          if (photon->GetNCaloPhotonMCLabels()>0){
            photonlabel = photon->GetCaloPhotonMCLabel(0);
            MCPhoton = (AliAODMCParticle*) fAODMCTrackArray->At(photon->GetCaloPhotonMCLabel(0));
            MCPhotonOld = (AliAODMCParticle*) fAODMCTrackArray->At(photon->GetCaloPhotonMCLabel(0));
          }
      } else {
        AliInfo("AODMCTrackArray could not be loaded");
        return;
      }
      if (photon->IsLargestComponentPhoton() || (photon->IsLargestComponentElectron() && photon->IsConversion())) {
         Bool_t isPrimary = fEventCuts->IsConversionPrimaryAOD(fInputEvent, MCPhoton, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
         if(isPrimary){
           isTruePhoton = kTRUE;
           if(photon->IsLargestComponentElectron() && photon->IsConversion()){
             truePhotonFromConv = kTRUE;
           } else if(photon->IsLargestComponentPhoton()){
             truePhotonNormal = kTRUE;
           }

         }

         // check for proper label, in case of conversion return photon label,otherwise do nothing
         Int_t tmplabel = GetProperLabel(MCPhoton);
         if(tmplabel != -99) photonlabel = tmplabel;
         MCPhoton = (AliAODMCParticle*) fAODMCTrackArray->At(photonlabel);

         isDecay = IsDecayPhoton(photonlabel);
         isPrompt = IsPromptPhoton(photonlabel);
         isFrag = IsFragPhoton(photonlabel);
         isSignal = isPrompt || isFrag;

      }
      
      mcIso = ProcessMCIsolation(photonlabel);

      for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
      {
        Double_t isoCh = mcIso.isolationCone.at(r) - fChargedRhoTimesArea[r];
        mcIsoCharged.push_back(isoCh);
      }
      for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
      {
        mcIsoNeutral.push_back(mcIso.isolationCone.at(r+fTrackIsolationR.size()));
      }
      for (UInt_t r = 0; r < mcIsoCharged.size(); r++)
      {
        mcIsoFull.push_back(mcIsoCharged.at(r)+ mcIsoNeutral.at(r)); // they should always have same length
      }
  }

  // Caluclate M02
  Double_t m02 = clus->GetM02();
  Double_t m20 = clus->GetM20();

  // Recalculate M02 (should only be used when using V2 clusterizer!)
  Int_t   nMaxima5x5 = 0; // output: number of local maxima on NxN region
  Float_t m02_5x5 = 0; // output: main shower shape eigenvalue
  Float_t m20_5x5 = 0; // output: second shower shape eigenvalue
  Float_t dispp= 0.; // output: dispersion
  Float_t dEta = 0.; // output: dispersion in eta (cols) direction
  Float_t dPhi = 0.; // output: dispersion in phi (cols) direction
  Float_t sEta = 0.; // output: shower shape in eta  (cols) direction
  Float_t sPhi = 0.; // output: shower shape in phi (rows) direction
  Float_t sEtaPhi = 0.; // output: shower shape on phi / eta directions term
  Float_t energy5x5 = 0; // output: sum of energy in NxN region
 
  AliVCaloCells* cells = fInputEvent->GetEMCALCells(); // input

  Bool_t onlyNeighbours = kTRUE; // input: Make sure all cells are adjacent to another cell in the cluster centred in absIdMax.
                                 // very important to have this option on true!

  Int_t cellDiff = 2; // input: NxN window to consider! 1 =3x3 2=5x5 3=7x7 
  Float_t cellEcut =0.1;// minimum cell energy to be considered in the shower shape recalculation
                        // should correspond to min cell cut of clusterizer
                      

  // Do recalculation for V2 clusterizer in 5x5 window
        GetCaloUtils()->GetEMCALRecoUtils()->RecalculateClusterShowerShapeParametersNxNCells
      (fGeomEMCAL, cells, clus,
       onlyNeighbours, cellDiff, cellEcut, 1000000,
       energy5x5, nMaxima5x5, m02_5x5, m20_5x5, dispp, dEta, dPhi, sEta, sPhi, sEtaPhi);

  cellDiff = 3;
  Int_t   nMaxima7x7 = 0; // output: number of local maxima on NxN region
  Float_t m02_7x7 = 0; // output: main shower shape eigenvalue
  Float_t m20_7x7 = 0; // output: second shower shape eigenvalue
  Float_t energy7x7 = 0; // output: sum of energy in NxN region
  
  // Do recalculation for V2 clusterizer in 7x7 window
        GetCaloUtils()->GetEMCALRecoUtils()->RecalculateClusterShowerShapeParametersNxNCells
      (fGeomEMCAL, cells, clus,
       onlyNeighbours, cellDiff, cellEcut, 1000000,
       energy7x7, nMaxima7x7, m02_7x7, m20_7x7, dispp, dEta, dPhi, sEta, sPhi, sEtaPhi);

  // ────────────────────────────────────────────────────────────────────────────────
   // get Max cell ID to get SM number
  Int_t absid    = ((AliCaloPhotonCuts*)fClusterCutsEMC)->FindLargestCellInCluster(clus, fInputEvent);
  UShort_t SMNumber = fGeomEMCAL->GetSuperModuleNumber(absid);

  // ────────────────────────────────────────────────────────────────────────────────

  const Int_t   nc = clus->GetNCells();
  Int_t   absCellIdList[nc];
  Float_t   maxEList[nc];

  // GetNLM
  Int_t nlm = fClusterCutsEMC->GetNumberOfLocalMaxima(clus,fInputEvent,absCellIdList,maxEList);
  
  // Get Efrac
  Float_t eFrac = GetExoticEnergyFraction(clus,fInputEvent);


  // Split cluster
  AliAODCaloCluster* clusSub1 = new AliAODCaloCluster();
  AliAODCaloCluster* clusSub2 = new AliAODCaloCluster();
  // split clusters according to their shares in the cluster (NLM == 1) needs to be treated differently
  if (nlm == 1){
    Int_t absCellIdFirst    = ((AliCaloPhotonCuts*)fClusterCutsEMC)->FindLargestCellInCluster(clus, fInputEvent);
    Int_t absCellIdSecond   = ((AliCaloPhotonCuts*)fClusterCutsEMC)->FindSecondLargestCellInCluster(clus, fInputEvent);

    ((AliCaloPhotonCuts*)fClusterCutsEMC)->SplitEnergy(absCellIdFirst, absCellIdSecond, clus, fInputEvent, fIsMC, clusSub1, clusSub2);
  } else if (nlm > 1 ){
    ((AliCaloPhotonCuts*)fClusterCutsEMC)->SplitEnergy(absCellIdList[0], absCellIdList[1], clus, fInputEvent, fIsMC, clusSub1, clusSub2);
  }

  TLorentzVector clusterVector1;
  clusSub1->GetMomentum(clusterVector1,vertex);
  TLorentzVector tmpvec1;
  tmpvec1.SetPxPyPzE(clusterVector1.Px(),clusterVector1.Py(),clusterVector1.Pz(),clusterVector1.E());
  
  // TLorentzvector with sub cluster 2
  TLorentzVector clusterVector2;
  clusSub2->GetMomentum(clusterVector2,vertex);
  TLorentzVector tmpvec2;
  tmpvec2.SetPxPyPzE(clusterVector2.Px(),clusterVector2.Py(),clusterVector2.Pz(),clusterVector2.E());
  // convert to AODConversionPhoton

  TLorentzVector vSum = tmpvec1 + tmpvec2;

  Float_t mass = -1;
  mass = vSum.M();

  // Take the trash out
  if(clusSub1)         delete clusSub1;
  if(clusSub2)         delete clusSub2;


  Double_t clusPx = photon->Px();
  Double_t clusPy = photon->Py();
  Double_t clusPz = photon->Pz();
  Double_t clusE = photon->E();
  Double_t clusM02 = m02;
  Double_t clusM20 = m20;
  Double_t clusV1SplitMass = mass;
  Double_t clusterIsoCharged = isoCharged.isolationCone.at(0);
  Double_t clusterIsoChargedLeft = isoCharged.backgroundLeft.at(2);

  //cout << "IsoCharged =" << isoCharged.isolationCone.at(0) << endl;
  
  Double_t trueClusterE = 0.;
  Double_t trueClusterPx = 0.;
  Double_t trueClusterPy = 0.;
  Double_t trueClusterPz = 0.;
  Double_t trueClusterMCIsoCharged1 = 0.;
  Double_t trueClusterMCIsoCharged2 = 0.;
  Double_t trueClusterMCIsoCharged3 = 0.;
  Double_t trueClusterMCIsoBckLeft = 0.;
  Int_t trueClusterIsSignal = 0;
  Bool_t trueClusterIsConv = kFALSE;
  Float_t fracLeadingLabel = 0;

  fBuffer_ClusterE.push_back(clusE);
  fBuffer_ClusterPx.push_back(clusPx);
  fBuffer_ClusterPy.push_back(clusPy);
  fBuffer_ClusterPz.push_back(clusPz);
  fBuffer_ClusterM02.push_back(clusM02);
  fBuffer_ClusterM02Recalc.push_back(m02_5x5);
  fBuffer_ClusterM20.push_back(clusM20);
  fBuffer_ClusterV1SplitMass.push_back(clusV1SplitMass);
  fBuffer_ClusterNLM.push_back(nlm);
  fBuffer_ClusterSM.push_back(SMNumber);
  fBuffer_ClusterEFrac.push_back(eFrac);
  fBuffer_ClusterIsoCharged1.push_back(isoCharged.isolationCone.at(0));
  if(isoCharged.isolationCone.size()>1) fBuffer_ClusterIsoCharged2.push_back(isoCharged.isolationCone.at(1));
  if(isoCharged.isolationCone.size()>2) fBuffer_ClusterIsoCharged3.push_back(isoCharged.isolationCone.at(2));
  fBuffer_ClusterIsoBckLeft.push_back(clusterIsoChargedLeft);

  fBuffer_ClusterMatchTrackdEta.push_back(trackdEta); 
  fBuffer_ClusterMatchTrackdPhi.push_back(trackdPhi); 
  fBuffer_ClusterMatchTrackP.push_back(trackP); // if no match this is -1
  fBuffer_ClusterMatchTrackPt.push_back(trackPt); 
  fBuffer_ClusterMatchTrackIsConv.push_back(isConvMatch); 
  
  
  if(fIsMC>0){
    if(isTruePhoton){
      trueClusterE = MCPhoton->E();
      trueClusterPx = MCPhoton->Px();
      trueClusterPy = MCPhoton->Py();
      trueClusterPz = MCPhoton->Pz();
      fracLeadingLabel = clus->GetClusterMCEdepFraction(0);
      trueClusterMCIsoCharged1 = mcIso.isolationCone.at(0);
      if(mcIso.isolationCone.size()>1)trueClusterMCIsoCharged2 = mcIso.isolationCone.at(1);
      if(mcIso.isolationCone.size()>2)trueClusterMCIsoCharged3 = mcIso.isolationCone.at(2);
      trueClusterMCIsoBckLeft = mcIso.backgroundLeft.at(2);

      TString headerName = fEventCuts->GetParticleHeaderName(photonlabel, fMCEvent, fInputEvent);
      if(((AliConvEventCuts*)fEventCuts)->GetSignalRejection() == 0){
        headerName = "";
      }
      Int_t tag = GetMCAnalysisUtils()->CheckOrigin(photonlabel, fMCEvent,headerName,1.);
      trueClusterIsSignal = tag;
      trueClusterIsConv = truePhotonFromConv;
    }
    fBuffer_TrueClusterE.push_back(trueClusterE);
    fBuffer_TrueClusterPx.push_back(trueClusterPx);
    fBuffer_TrueClusterPy.push_back(trueClusterPy);
    fBuffer_TrueClusterPz.push_back(trueClusterPz);
    fBuffer_TrueClusterLeadingEFrac.push_back(fracLeadingLabel);
    fBuffer_TrueClusterMCIsoCharged1.push_back(trueClusterMCIsoCharged1); 
    fBuffer_TrueClusterMCIsoCharged2.push_back(trueClusterMCIsoCharged2); 
    fBuffer_TrueClusterMCIsoCharged3.push_back(trueClusterMCIsoCharged3); 
    fBuffer_TrueClusterMCIsoBckLeft.push_back(trueClusterMCIsoBckLeft); 
    fBuffer_TrueClusterMCTag.push_back(trueClusterIsSignal); 
    fBuffer_TrueClusterIsConv.push_back(trueClusterIsConv);
  }



}
void AliAnalysisTaskGammaIsoTree::FillCaloHistos(AliAODCaloCluster* clus, AliAODConversionPhoton* photon, isoValues isoCharged,vector<Double32_t> isoNeutral,vector<Double32_t> isoCell,Int_t tmptag, Double_t weight){
    if(!clus) return;

    Double_t vertex[3] = {0};
    InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

    fCaloPtBeforeAcc->Fill(photon->Pt(),weight);
    //if (!IsInEMCalAcceptance(photon)) return;

    Bool_t isTruePhoton = kFALSE;
    Bool_t isDecay = kFALSE;
    Bool_t isPrompt = kFALSE;
    Bool_t isFrag = kFALSE;
    Bool_t isSignal = kFALSE;
    AliAODMCParticle *MCPhoton = NULL;

    // check MC properties
    if(fIsMC>0){
        const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
        Double_t mcProdVtxX   = primVtxMC->GetX();
        Double_t mcProdVtxY   = primVtxMC->GetY();
        Double_t mcProdVtxZ   = primVtxMC->GetZ();
        Int_t photonlabel = -1;
        if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
        if (fAODMCTrackArray){
            if (photon->GetIsCaloPhoton() == 0) AliFatal("CaloPhotonFlag has not been set task will abort");
            if (photon->GetNCaloPhotonMCLabels()>0) {
              photonlabel = photon->GetCaloPhotonMCLabel(0);
              MCPhoton = (AliAODMCParticle*) fAODMCTrackArray->At(photonlabel);
            }
        } else {
          AliInfo("AODMCTrackArray could not be loaded");
          return;
        }
        if (photon->IsLargestComponentPhoton() || (photon->IsLargestComponentElectron() && photon->IsConversion())) {
          Bool_t isPrimary = fEventCuts->IsConversionPrimaryAOD(fInputEvent, MCPhoton, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
          if(isPrimary) isTruePhoton = kTRUE;

          // check for proper label, in case of conversion return photon label,otherwise do nothing
          Int_t tmplabel = GetProperLabel(MCPhoton);
          if(tmplabel != -99) photonlabel = tmplabel;
          MCPhoton = (AliAODMCParticle*) fAODMCTrackArray->At(photonlabel);

          isDecay = IsDecayPhoton(photonlabel);
          isPrompt = IsPromptPhoton(photonlabel);
          isFrag = IsFragPhoton(photonlabel);
          isSignal = isPrompt || isFrag;
        }
    }

    fCaloPt->Fill(photon->Pt(),weight);
    fCaloE->Fill(clus->E(),weight);
    if(tmptag<2 ){
       fCaloPtTaggedCalo->Fill(photon->Pt(),weight);
    } else{
       fCaloPtTaggedAsDecayCalo->Fill(photon->Pt(),weight);
    }

    fCaloRho->Fill(fChargedRho,photon->Pt(),fWeightJetJetMC);

    for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
    {
      fCaloRhoTimesArea[r]->Fill(fChargedRhoTimesArea[r], photon->Pt(),fWeightJetJetMC);
      fCaloRhoTimesAreaLeft[r]->Fill(isoCharged.backgroundLeft.at(r), photon->Pt(),fWeightJetJetMC);
      fCaloRhoTimesAreaRight[r]->Fill(isoCharged.backgroundRight.at(r), photon->Pt(),fWeightJetJetMC);
      fCaloRhoTimesAreaBack[r]->Fill(isoCharged.backgroundBack.at(r), photon->Pt(),fWeightJetJetMC);
      fCaloIsoCharged[r]->Fill(isoCharged.isolationCone.at(r),photon->Pt(),weight);
      fCaloIsoRawCharged[r]->Fill(isoCharged.isolationCone.at(r)+fChargedRhoTimesArea[r],photon->Pt(),weight);

      // Do check if it is in similar range to photon
      if(fIsMC && fDebug){ // debug
           if((isoCharged.isolationCone.at(r)< (photon->Pt()*1.2)) && (isoCharged.isolationCone.at(r)> (photon->Pt()*0.8)) && (photon->Pt()>15.)){
          //  if(kFALSE){
              cout << "Photon pt = " << photon->Pt() << "Iso is = " << isoCharged.isolationCone.at(r) << endl;

              // do MC checks
              Int_t photonlabel = -1;
                const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
              Double_t mcProdVtxX   = primVtxMC->GetX();
              Double_t mcProdVtxY   = primVtxMC->GetY();
              Double_t mcProdVtxZ   = primVtxMC->GetZ();
              if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
              if (fAODMCTrackArray){
                  if (photon->GetIsCaloPhoton() == 0) AliFatal("CaloPhotonFlag has not been set task will abort");
                  if (photon->GetNCaloPhotonMCLabels()>0) {
                    photonlabel = photon->GetCaloPhotonMCLabel(0);
                    MCPhoton = (AliAODMCParticle*) fAODMCTrackArray->At(photonlabel);

                    AliAODMCParticle* MCMother = (AliAODMCParticle*) fAODMCTrackArray->At(MCPhoton->GetMother());
                    cout << "Leading PDG is = " << MCPhoton->PdgCode() << endl;
                    if(MCMother){
                        cout << "Mother PDG is = " << MCMother->GetPdgCode() << endl;
                        cout << "Mother Pt is = " << MCMother->Pt() << endl;
                        Double_t vx= MCMother->Xv();
                        Double_t vy= MCMother->Yv();
                        Double_t vz= MCMother->Zv();
                        Double_t r = TMath::Sqrt(vx*vx + vy*vy + vz*vz);
                        cout << "Mother R is = " << r << endl;

                        AliAODMCParticle* Daughter1 = (AliAODMCParticle*) fAODMCTrackArray->At(MCMother->GetDaughterFirst());
                        AliAODMCParticle* Daughter2 = (AliAODMCParticle*) fAODMCTrackArray->At(MCMother->GetDaughterLast());
                        cout << "Daughter 1; pos = " << MCMother->GetDaughterFirst() << " PDG = " << Daughter1->GetPdgCode() << " Pt = " << Daughter1->Pt()<< endl;
                        cout << "Daughter 2; pos = " << MCMother->GetDaughterLast()  << " PDG = " << Daughter2->GetPdgCode()<< " Pt = " << Daughter2->Pt()<<endl;
                        cout << "----------------" << endl;
                       ((AliCaloPhotonCuts*)fClusterCutsEMC)->MatchTracksToClusters(fInputEvent,fWeightJetJetMC,kTRUE, fMCEvent);
                        cout << "-------------------" << endl;
                        // Check distance and check track matching condition
                        Double_t trackEta = Daughter2->Eta();
                        Double_t trackPhi = Daughter2->Phi();
                        Double_t photonEta = photon->Eta();
                        Double_t photonPhi = photon->Phi();

                        Double_t deta = TMath::Abs(trackEta - photonEta);
                        Double_t dphi = TMath::Abs(trackPhi - photonPhi);

                        TF1* fFuncPtDepEta = new TF1("funcEta15", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
                        fFuncPtDepEta->SetParameters(0.04, 0.010, 2.5);
                        TF1* fFuncPtDepPhi = new TF1("funcPhi15", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
                        fFuncPtDepPhi->SetParameters(0.09, 0.015, 2.);

                        Double_t maxEta = fFuncPtDepEta->Eval(Daughter2->Pt());
                        Double_t maxPhi = fFuncPtDepPhi->Eval(Daughter2->Pt());

                        cout << "dEta = " << deta << " maxdEta = " << maxEta << endl;
                        cout << "dPhi = " << dphi << " maxdPhi = " << maxPhi << endl;

                        cout << "Number of conv:" << fReaderGammas->GetEntriesFast() << endl;
                        cout << "CheckClusterForTrackMatch = " << fClusterCutsEMC->CheckClusterForTrackMatch(clus) << endl;
                        Int_t labelTrackHighest = -1;
                        if(fClusterCutsEMC->GetHighestPtMatchedTrackToCluster(fInputEvent,clus,labelTrackHighest)){
                            Double_t pt = ((AliAODTrack*)fInputEvent->GetTrack(labelTrackHighest))->Pt();
                            cout << "Matched track with pt " << pt << endl;
                        }
                        if(fClusterCutsEMC->GetClosestMatchedTrackToCluster(fInputEvent,clus,labelTrackHighest)){
                            Double_t pt = ((AliAODTrack*)fInputEvent->GetTrack(labelTrackHighest))->Pt();
                            cout << "Matched closest track with pt " << pt << endl;
                        }
                        cout << "IsMatchedWithConv = " << IsMatchedWithConv(clus,fClusterCutsEMC) << endl;
                        cout << "IsMatchedWithConv = " << IsMatchedWithConv(clus,fClusterCutsEMC) << endl;
                        AliAODTrack *tempTrack = 0x0;
                        for(Int_t t=0;t<fInputEvent->GetNumberOfTracks();t++) {
                          tempTrack = static_cast<AliAODTrack*>(fInputEvent->GetTrack(t));
                          if(tempTrack){
                            if((tempTrack->Pt()< (1.2* photon->Pt())) && (tempTrack->Pt() > (photon->Pt()*0.8)) && (photon->Pt()>15.)){
                            if(TrackIsSelectedAOD(tempTrack)){
                              cout << "Track Pt = " << tempTrack->Pt() <<"TrackID="<<tempTrack->GetID() << endl;  
                              vector<Int_t> matchedClus(0);
                              Int_t id = tempTrack->GetID();
                              matchedClus = fClusterCutsEMC->GetCaloTrackMatcherInstance()->GetMatchedClusterIDsForTrack(fInputEvent,id,0.4);
                              cout << "matched clus size" << matchedClus.size() << endl;


                              // do own propagation
                              Double_t xyz[3] = {0}, pxpypz[3] = {0}, cv[21] = {0};
                              tempTrack->GetPxPyPz(pxpypz);
                              tempTrack->GetXYZ(xyz);
                              tempTrack->GetCovarianceXYZPxPyPz(cv);

                              AliExternalTrackParam* trackParam = new AliExternalTrackParam(xyz,pxpypz,cv,tempTrack->Charge());

                              AliExternalTrackParam emcParam(*trackParam);
                              Float_t eta, phi, pt;
                              //propagate tracks to emc surfaces
                              if (AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(&emcParam, 440., 0.139, 20., eta, phi, pt)) {
                                cout << "managed to Propagate" << endl;
                                Double_t etaDiff = TMath::Abs(eta - photonEta);
                                Double_t phiDiff = TMath::Abs(phi - photonPhi);
                                cout << "etaDiff=" << etaDiff << endl;
                                cout << "phiDiff=" << phiDiff << endl;

                                cout << "E/p = " << clus->E()/tempTrack->P() << endl;
                              }
                              delete trackParam;

                            }
                            }
                          }
                        }
                        vector<Int_t> labelsMatched(0);
                        labelsMatched = fClusterCutsEMC->GetCaloTrackMatcherInstance()->GetMatchedTrackIDsForCluster(fInputEvent,clus->GetID(),0.4,0,0.4,0);
                       cout <<"SIZE = " << labelsMatched.size() << endl;
                       for (Int_t conv = 0; conv < fReaderGammas->GetEntriesFast(); conv++)
                       {
                          AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(conv);
                            for (Int_t i = 0;i < 2;i++){
                              Int_t tracklabel = PhotonCandidate->GetLabel(i);
                              AliVTrack *inTrack = 0x0;
                              if(((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data()))->AreAODsRelabeled()){
                                inTrack = dynamic_cast<AliVTrack*>(fInputEvent->GetTrack(tracklabel));
                              } else {
                                for(Int_t ii=0;ii<fInputEvent->GetNumberOfTracks();ii++) {
                                  inTrack = dynamic_cast<AliVTrack*>(fInputEvent->GetTrack(ii));
                                  if(inTrack){
                                    if(inTrack->GetID() == tracklabel) {
                                      break;
                                    }
                                  }
                                }
                              }
                              
                          }
                       }
//
                    }

                  }
              }
              if (photon->IsLargestComponentPhoton() || (photon->IsLargestComponentElectron() && photon->IsConversion())) {
                  Bool_t isPrimary = fEventCuts->IsConversionPrimaryAOD(fInputEvent, MCPhoton, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
                  if(isPrimary) isTruePhoton = kTRUE;

                  // check for proper label, in case of conversion return photon label,otherwise do nothing
                  Int_t tmplabel = GetProperLabel(MCPhoton);
                  if(tmplabel != -99) photonlabel = tmplabel;
                  MCPhoton = (AliAODMCParticle*) fAODMCTrackArray->At(photonlabel);

                  isDecay = IsDecayPhoton(photonlabel);
                  isPrompt = IsPromptPhoton(photonlabel);
                  isFrag = IsFragPhoton(photonlabel);
                  isSignal = isPrompt || isFrag;

                  cout <<"IsDecay = " << isDecay << endl;
                  cout <<"IsPrompt = " << isPrompt << endl;
                  cout <<"IsFrag = " << isFrag << endl;
              }
           }

      }

      
      for (UInt_t e = 0; e < fTrackIsolationE.size(); e++)
      {
        if(isoCharged.isolationCone.at(r) < fTrackIsolationE.at(e)){
          fCaloPtIsoCharged[r][e]->Fill(photon->Pt(),weight);
          if(tmptag<2 )fCaloPtTaggedCaloIsoCharged[r][e]->Fill(photon->Pt(),weight); // not tagged by calo
        }
      }
    }

    for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
    {
      fCaloIsoNeutral[r]->Fill(isoNeutral.at(r),photon->Pt(),weight);
      fCaloIsoFull[r]->Fill(isoNeutral.at(r) + isoCharged.isolationCone.at(r),photon->Pt(),weight);
      fCaloIsoCell[r]->Fill(isoCell.at(r),photon->Pt(),weight);
      for (UInt_t e = 0; e < fNeutralIsolationE.size(); e++)
      {
        if(isoNeutral.at(r) < fNeutralIsolationE.at(e))fCaloPtIsoNeutral[r][e]->Fill(photon->Pt(),weight);
        if(isoCell.at(r) < fNeutralIsolationE.at(e)) fCaloPtIsoCell[r][e]->Fill(photon->Pt(),weight);
        if((isoNeutral.at(r) + isoCharged.isolationCone.at(r)) < fNeutralIsolationE.at(e)) fCaloPtIsoFull[r][e]->Fill(photon->Pt(),weight);
        if(tmptag<2){
            if(isoNeutral.at(r) < fNeutralIsolationE.at(e))fCaloPtTaggedCaloIsoNeutral[r][e]->Fill(photon->Pt(),weight);
            if(isoCell.at(r) < fNeutralIsolationE.at(e)) fCaloPtTaggedCaloIsoCell[r][e]->Fill(photon->Pt(),weight);
            if((isoNeutral.at(r) + isoCharged.isolationCone.at(r)) < fNeutralIsolationE.at(e)) fCaloPtTaggedCaloIsoFull[r][e]->Fill(photon->Pt(),weight);
        }
      }
    }

    //
    // ─── FILL INV MASS HISTOS ────────────────────────────────────────
    //
    for (Int_t c = 0; c < fClusterEMCalCandidates->GetEntries(); c++)
    {
      if(((AliAODCaloCluster*)fClusterEMCalCandidates->At(c))->GetID() == clus->GetID()) continue;
      // TLorentzvector with cluster
      TLorentzVector clusterVector;
      ((AliAODCaloCluster*)fClusterEMCalCandidates->At(c))->GetMomentum(clusterVector,vertex);

      TLorentzVector* tmpvec = new TLorentzVector();
      tmpvec->SetPxPyPzE(clusterVector.Px(),clusterVector.Py(),clusterVector.Pz(),clusterVector.E());

      // convert to AODConversionPhoton
      AliAODConversionPhoton *othercluster = new AliAODConversionPhoton(tmpvec);
      if(!othercluster){ delete tmpvec; continue;}

      AliAODConversionMother* pi0cand = new AliAODConversionMother(photon,othercluster);

      // check mass window
      Double_t mass = pi0cand->M();

      fCaloInvMass->Fill(mass,photon->Pt(),weight);

      if(fIsMC>0){
        if(isTruePhoton){
          fCaloTrueInvMass->Fill(mass,photon->Pt(),weight);
          if(isDecay){
            fCaloTrueInvMass_FromDecay->Fill(mass,photon->Pt(),weight);
          } else if (isSignal){
            fCaloTrueInvMass_FromDirect->Fill(mass,photon->Pt(),weight);
          }
        }
      }

      for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
      {
        if((isoCharged.isolationCone.at(r) > fAntiIsolationE[0]) && (isoCharged.isolationCone.at(r) < fAntiIsolationE[1])) fCaloInvMassAntiChargedIsolated[r]->Fill(mass,photon->Pt(),weight);
        for (UInt_t e = 0; e < fTrackIsolationE.size(); e++)
        {
          if(isoCharged.isolationCone.at(r) < fTrackIsolationE.at(e)){
            fCaloInvMassChargedIsolated[r][e]->Fill(mass,photon->Pt(),weight);
          }
        }
      }

      for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
      {
        if((isoNeutral.at(r) > fAntiIsolationE[0]) && (isoNeutral.at(r) < fAntiIsolationE[1])) fCaloInvMassAntiNeutralIsolated[r]->Fill(mass,photon->Pt(),weight);
        if((isoCell.at(r) > fAntiIsolationE[0]) && (isoCell.at(r) < fAntiIsolationE[1])) fCaloInvMassAntiCellIsolated[r]->Fill(mass,photon->Pt(),weight);
        if(((isoNeutral.at(r) + isoCharged.isolationCone.at(r)) > fAntiIsolationE[0]) && ((isoNeutral.at(r) + isoCharged.isolationCone.at(r)) < fAntiIsolationE[1])) fCaloInvMassAntiFullIsolated[r]->Fill(mass,photon->Pt(),weight);
        for (UInt_t e = 0; e < fNeutralIsolationE.size(); e++)
        {
          if(isoNeutral.at(r) < fNeutralIsolationE.at(e)) fCaloInvMassNeutralIsolated[r][e]->Fill(mass,photon->Pt(),weight);
          if(isoCell.at(r) < fNeutralIsolationE.at(e)) fCaloInvMassCellIsolated[r][e]->Fill(mass,photon->Pt(),weight);
          if((isoNeutral.at(r) + isoCharged.isolationCone.at(r)) < fNeutralIsolationE.at(e)) fCaloInvMassFullIsolated[r][e]->Fill(mass,photon->Pt(),weight);

        }
      }

      delete tmpvec;
      delete othercluster;
      delete pi0cand;
  }
}

Float_t AliAnalysisTaskGammaIsoTree::GetExoticEnergyFraction(AliVCluster *cluster, AliVEvent *event){
    Float_t exoticEnergyFrac =-999;
    if (!cluster) {
      AliInfo("Cluster pointer null!");
      return 999;
    }
    AliVCaloCells* cells    = event->GetEMCALCells();

    Int_t largestCellID     = ((AliCaloPhotonCuts*)fClusterCutsEMC)->FindLargestCellInCluster(cluster,event);
    Float_t ecell1          = cells->GetCellAmplitude(largestCellID);
    Float_t eCross          = ((AliCaloPhotonCuts*)fClusterCutsEMC)->GetECross(largestCellID,cells);
    exoticEnergyFrac = 1-eCross/ecell1;

    return exoticEnergyFrac;

}
Bool_t AliAnalysisTaskGammaIsoTree::IsMatchedWithConv(AliAODCaloCluster* clus, AliCaloPhotonCuts* cuts){
  Bool_t matched = kFALSE;
  if(!fReaderGammas) fReaderGammas    = fV0Reader->GetReconstructedGammas();
  for (Int_t conv = 0; conv < fReaderGammas->GetEntriesFast(); conv++)
  {
     AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(conv);
     matched = cuts->MatchConvPhotonToCluster(PhotonCandidate,clus, fInputEvent );
     if(matched) return matched;
  }
  return matched;
}
Bool_t AliAnalysisTaskGammaIsoTree::IsSameTrack(Int_t id1, Int_t id2){
    if((id1 == -999) || (id2 == -999)){
        cout << "ERROR: Track info is missing for one track!" << endl;
        return kFALSE;
    }
    Int_t esdID1 = id1;
    Int_t esdID2 = id2;
    if(id1<0) esdID1 = (-1 * id1) - 1;
    if(id2<0) esdID2 = (-1 * id2) - 1;
    if(esdID1 == esdID2){
        return kTRUE;
    } else{
        return kFALSE;
    }

}

Bool_t AliAnalysisTaskGammaIsoTree::IsInEMCalAcceptance(AliAODConversionPhoton *photon)
{
    Double_t eta = photon->GetPhotonEta();
    Double_t phi = photon->GetPhotonPhi();

    Double_t etaMin = fClusterCutsEMC->GetMinEtaCut();
    Double_t etaMax = fClusterCutsEMC->GetMaxEtaCut();
    Double_t phiMin = fClusterCutsEMC->GetMinPhiCut();
    Double_t phiMax = fClusterCutsEMC->GetMaxPhiCut();
    // cout << phi << endl;
    // cout << eta << endl;
    if (phi < 0)
        phi += 2 * TMath::Pi();
    if ((eta < etaMin) || (eta > etaMax))
        return kFALSE;
    if ((phi < phiMin) || (phi > phiMax))
        return kFALSE;
    return kTRUE;
}
Bool_t AliAnalysisTaskGammaIsoTree::IsInEMCalAcceptance(AliAODMCParticle *part)
{
    Double_t eta = part->Eta();
    Double_t phi = part->Phi();

    Double_t etaMin = fClusterCutsEMC->GetMinEtaCut();
    Double_t etaMax = fClusterCutsEMC->GetMaxEtaCut();
    Double_t phiMin = fClusterCutsEMC->GetMinPhiCut();
    Double_t phiMax = fClusterCutsEMC->GetMaxPhiCut();

    if (phi < 0)
        phi += 2 * TMath::Pi();
    // cout << phi << endl;
    // cout << eta << endl;
    if ((eta < etaMin) || (eta > etaMax))
        return kFALSE;
    if ((phi < phiMin) || (phi > phiMax))
        return kFALSE;
    return kTRUE;
}
Bool_t AliAnalysisTaskGammaIsoTree::IsTrueConversionPhoton(AliAODConversionPhoton* photon)
{
    const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
    Double_t mcProdVtxX   = primVtxMC->GetX();
    Double_t mcProdVtxY   = primVtxMC->GetY();
    Double_t mcProdVtxZ   = primVtxMC->GetZ();

    if ((photon->GetMCLabelPositive() == -1) || (photon->GetMCLabelNegative() == -1))
        return false;
    if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));

    AliAODMCParticle *posDaughter = static_cast<AliAODMCParticle *>(fAODMCTrackArray->At(photon->GetMCLabelPositive()));
    AliAODMCParticle *negDaughter = static_cast<AliAODMCParticle *>(fAODMCTrackArray->At(photon->GetMCLabelNegative()));

    if (posDaughter == NULL || negDaughter == NULL)
        return kFALSE;
    Int_t pdgCode[2] = {TMath::Abs(posDaughter->GetPdgCode()), TMath::Abs(negDaughter->GetPdgCode())};
    if (posDaughter->GetMother() != negDaughter->GetMother())
        return kFALSE;
    if (posDaughter->GetMother() == -1)
        return kFALSE;
    if (pdgCode[0] != 11 || pdgCode[1] != 11)
        return kFALSE;
    if (posDaughter->GetPdgCode() == negDaughter->GetPdgCode())
        return kFALSE;

    // cout << "Get Mother" << posDaughter->GetMother() << endl;
    AliAODMCParticle *mother = static_cast<AliAODMCParticle *>(fAODMCTrackArray->At(posDaughter->GetMother()));
    if (mother->GetPdgCode() != 22)
        return kFALSE;
    if(((posDaughter->GetMCProcessCode())) != 5 || ((negDaughter->GetMCProcessCode())) != 5){
        return kFALSE;// check if the daughters come from a conversion
    }
    Bool_t isPrimary = fEventCuts->IsConversionPrimaryAOD(fInputEvent, mother, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
    if(!isPrimary) return kFALSE;

    return kTRUE;
}
Int_t AliAnalysisTaskGammaIsoTree::GetConvPhotonMCLabel(AliAODConversionPhoton *photon)
{
    if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    AliAODMCParticle *posDaughter = static_cast<AliAODMCParticle *>(fAODMCTrackArray->At(photon->GetMCLabelPositive()));
    AliAODMCParticle *negDaughter = static_cast<AliAODMCParticle *>(fAODMCTrackArray->At(photon->GetMCLabelNegative()));

    if (posDaughter->GetMother() != negDaughter->GetMother())
        return -1;

    return posDaughter->GetMother();
}
Bool_t AliAnalysisTaskGammaIsoTree::IsDecayPhoton(Int_t label){ // i.e. not direct photon
    Bool_t isFromDecay = kFALSE;
    TString headerName = fEventCuts->GetParticleHeaderName(label, fMCEvent, fInputEvent);
    if(((AliConvEventCuts*)fEventCuts)->GetSignalRejection() == 0){
       headerName = "";
    }
    Int_t tag = GetMCAnalysisUtils()->CheckOrigin(label, fMCEvent,headerName,1.);
    Bool_t isPhoton = GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPhoton);
    if(!isPhoton) return kFALSE;

    Bool_t isPi0Decay = GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPi0Decay);
    Bool_t isEtaDecay = GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCEtaDecay);
    Bool_t isOtherDecay = GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCOtherDecay);

    isFromDecay = (isPi0Decay || isEtaDecay || isOtherDecay);

    return isFromDecay;

}
Bool_t AliAnalysisTaskGammaIsoTree::IsDecayPhoton(AliAODConversionPhoton *photon){
    Bool_t isFromDecay = kFALSE;

    if (!IsTrueConversionPhoton(photon))
        return kFALSE;
    Int_t motherlabel = GetConvPhotonMCLabel(photon);
    TString headerName = fEventCuts->GetParticleHeaderName(motherlabel, fMCEvent, fInputEvent);
    if(((AliConvEventCuts*)fEventCuts)->GetSignalRejection() == 0){
       headerName = "";
    }
    Int_t tag = GetMCAnalysisUtils()->CheckOrigin(motherlabel, fMCEvent,headerName,photon->E());
    Bool_t isPhoton = GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPhoton);
    if(!isPhoton) return kFALSE;

    Bool_t isPi0Decay = GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPi0Decay);
    Bool_t isEtaDecay = GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCEtaDecay);
    Bool_t isOtherDecay = GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCOtherDecay);

    isFromDecay = (isPi0Decay || isEtaDecay || isOtherDecay);

    return isFromDecay;
}

Int_t AliAnalysisTaskGammaIsoTree::CheckClustersForMCContribution(Int_t mclabel, TList *vclus)
{
    Int_t clusterLabel = -1; // position of cluster in array where mc label was found as contribution
    for (Int_t p = 0; p < vclus->GetEntries(); p++)
    {
        AliAODCaloCluster *clus = (AliAODCaloCluster *)vclus->At(p);
        if (!clus)
            continue;
        Int_t *mclabelsCluster = clus->GetLabels();
        if (clus->GetNLabels() > 0)
        {
            for (Int_t k = 0; k < (Int_t)clus->GetNLabels(); k++)
            {
                if (mclabelsCluster[k] == mclabel)
                    clusterLabel = p;
            }
        }
    }

    return clusterLabel;
}
Int_t AliAnalysisTaskGammaIsoTree::CheckConvForMCContribution(Int_t mclabel, TList *vconv)
{
  Int_t convLabel = -1; // position of conversion in array where mc label was found as contribution
   for (Int_t p = 0; p < vconv->GetEntries(); p++)
    {
        AliAODConversionPhoton *photon = (AliAODConversionPhoton *)vconv->At(p);
        if (!photon)
            continue;

        Int_t label = GetConvPhotonMCLabel(photon);
        if (mclabel == label)
            convLabel = p;
    }
    return convLabel;
}
Bool_t AliAnalysisTaskGammaIsoTree::IsPromptPhoton(AliAODConversionPhoton *photon){
    Bool_t isPromptPhoton = kFALSE;
    if (!IsTrueConversionPhoton(photon))
        return kFALSE;
    Int_t motherlabel = GetConvPhotonMCLabel(photon);

    // find header of this particle
    TString headerName = fEventCuts->GetParticleHeaderName(motherlabel, fMCEvent, fInputEvent);
    if(((AliConvEventCuts*)fEventCuts)->GetSignalRejection() == 0){
       headerName = "";
    }

    Int_t tag = GetMCAnalysisUtils()->CheckOrigin(motherlabel, fMCEvent,headerName,photon->E());
    if(!GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPhoton)) return kFALSE;
    isPromptPhoton = GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPrompt);
    return isPromptPhoton;
}
Bool_t AliAnalysisTaskGammaIsoTree::IsPromptPhoton(Int_t label){
    Bool_t isPromptPhoton = kFALSE;

    // find header of this particle
    TString headerName = fEventCuts->GetParticleHeaderName(label, fMCEvent, fInputEvent);
    if(((AliConvEventCuts*)fEventCuts)->GetSignalRejection() == 0){
       headerName = "";
    }
    Int_t tag = GetMCAnalysisUtils()->CheckOrigin(label, fMCEvent,headerName,1.);
    if(!GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPhoton)) return kFALSE;
    isPromptPhoton = GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPrompt);


    return isPromptPhoton;
}
Bool_t AliAnalysisTaskGammaIsoTree::IsFragPhoton(AliAODConversionPhoton *photon){
    Bool_t isFragPhoton = kFALSE;
    if (!IsTrueConversionPhoton(photon))
        return kFALSE;
    Int_t motherlabel = GetConvPhotonMCLabel(photon);

    // find header of this particle
    TString headerName = fEventCuts->GetParticleHeaderName(motherlabel, fMCEvent, fInputEvent);
    if(((AliConvEventCuts*)fEventCuts)->GetSignalRejection() == 0){
       headerName = "";
    }

    Int_t tag = GetMCAnalysisUtils()->CheckOrigin(motherlabel, fMCEvent,headerName,photon->E());
    if(!GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPhoton)) return kFALSE;
    isFragPhoton = GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCFragmentation);

    return isFragPhoton;
}
Bool_t AliAnalysisTaskGammaIsoTree::IsFragPhoton(Int_t label){
    Bool_t isFragPhoton = kFALSE;

    // find header of this particle
    TString headerName = fEventCuts->GetParticleHeaderName(label, fMCEvent, fInputEvent);
    if(((AliConvEventCuts*)fEventCuts)->GetSignalRejection() == 0){
       headerName = "";
    }

    Int_t tag = GetMCAnalysisUtils()->CheckOrigin(label, fMCEvent,headerName,1.);
    if(!GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCPhoton)) return kFALSE;
    isFragPhoton = GetMCAnalysisUtils()->CheckTagBit(tag, AliMCAnalysisUtils::kMCFragmentation);

    return isFragPhoton;
}
Bool_t AliAnalysisTaskGammaIsoTree::IsWithinRadiusEMCal(Double_t eta, Double_t phi, Double_t riso){

    // get max used area of emcal and subtract R
    Double_t etamin = fClusterCutsIsolationEMC->GetMinEtaCut()+riso;
    Double_t etamax = fClusterCutsIsolationEMC->GetMaxEtaCut()-riso;
    Double_t phimin = fClusterCutsIsolationEMC->GetMinPhiCut()+riso;
    Double_t phimax = fClusterCutsIsolationEMC->GetMaxPhiCut()-riso;

    // make sure phi is correct format
    if (phi < 0) phi += 2*TMath::Pi();

    if((eta>etamax) || (eta<etamin)) return kFALSE;
    if((phi>phimax) || (phi<phimin)) return kFALSE;

    return kTRUE;
}
Bool_t AliAnalysisTaskGammaIsoTree::IsWithinRadiusTPC(Double_t eta, Double_t phi, Double_t riso){
    Double_t etaRange= fEtaCut-riso;
    if(TMath::Abs(eta)>etaRange) return kFALSE;
    return kTRUE;
}

// function to return proper label to photon in case of conversion electron
Int_t AliAnalysisTaskGammaIsoTree::GetProperLabel(AliAODMCParticle* mcpart){
  if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  AliAODMCParticle* currentParticle = mcpart;
  Int_t label = -99;
  if (TMath::Abs(currentParticle->GetPdgCode()) != 11 ) return label;
  if (TMath::Abs(currentParticle->GetPdgCode()) == 11 ){
    if (currentParticle->GetMother() != -1){
      AliAODMCParticle* particleMother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(currentParticle->GetMother()));
      if (particleMother){
        if (TMath::Abs(particleMother->GetPdgCode()) == 22){
          label = currentParticle->GetMother();
          currentParticle = particleMother;
        }
      }
    }
  }
  return label;
}
