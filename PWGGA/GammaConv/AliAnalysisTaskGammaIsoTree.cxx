/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Nicolas Schmidt                                               *
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
  fExtraClusterInfo(),
  fExtraClusterInfoBackground(),
  fDataEvtHeader(),
  fMCEvtHeader(),
  fConvIsoInfo(),
  fCaloIsoInfo(),
  fGeomEMCAL(NULL),
  fCorrTaskSetting(""),
  fEventCuts(NULL),
  fClusterCutsEMC(NULL),
  fClusterCutsIsolationEMC(NULL),
  fClusterCutsTaggingEMC(NULL),
  fClusterCutsPHOS(NULL),
  fConvCuts(NULL),
  fCaloUtils(NULL),
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
  fConvIsoRawCharged(),
  fConvIsoRawNeutral(),
  fConvIsoRawFull(),
  fConvIsoCell(),
  fConvIsoCorr(),
  fConvRho(NULL),
  fConvRhoTimesArea(NULL),
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
  fConvTrueIsoRawCharged(),
  fConvTrueIsoRawNeutral(),
  fConvTrueIsoRawFull(),
  fConvTrueIsoCorr(),
  fConvTrueIsoCell(),
  fConvTrueIsoRawCharged_FromDecay(),
  fConvTrueIsoRawNeutral_FromDecay(),
  fConvTrueIsoRawFull_FromDecay(),
  fConvTrueIsoCell_FromDecay(),
  fConvTrueIsoRawCharged_FromDirect(),
  fConvTrueIsoRawNeutral_FromDirect(),
  fConvTrueIsoRawFull_FromDirect(),
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
  fConvTrueRecPtIsoCharged(),
  fConvTrueRecPtIsoNeutral(),
  fConvTrueRecPtIsoFull(),
  fConvTrueRecPtIsoCell(),
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
  fCaloPt(NULL),
  fCaloPtBeforeAcc(NULL),
  fCaloPtTaggedCalo(NULL),
  fCaloPtTaggedAsDecayCalo(NULL),
  fCaloIsoRawCharged(),
  fCaloIsoRawNeutral(),
  fCaloIsoRawFull(),
  fCaloIsoCell(),
  fCaloIsoCorr(),
  fCaloRho(NULL),
  fCaloRhoTimesArea(NULL),
  fCaloTruePt(NULL),
  fCaloTruePtPrimary(NULL),
  fCaloTruePtDecay(NULL),
  fCaloTruePtDecayFoundOtherInCluster(NULL),
  fCaloTruePtDecayOtherInAcc(NULL),
  fCaloTruePtDecayOtherInAccAboveMinEnergy(NULL),
  fCaloTruePtTaggedCalo(NULL),
  fCaloTruePtTaggedAsDecayCalo(NULL),
  fCaloTrueRecPt(NULL),
  fCaloTrueRecPtPrimary(NULL),
  fCaloTrueRecPtDecay(NULL),
  fCaloTrueRecPtDecayFoundOtherInCluster(NULL),
  fCaloTrueRecPtDecayOtherInAcc(NULL),
  fCaloTrueRecPtDecayOtherInAccAboveMinEnergy(NULL),
  fCaloTrueRecPtTaggedCalo(NULL),
  fCaloTrueRecPtTaggedAsDecayCalo(NULL),
  fCaloTrueIsoRawCharged(),
  fCaloTrueIsoRawNeutral(),
  fCaloTrueIsoRawFull(),
  fCaloTrueIsoCorr(),
  fCaloTrueIsoCell(),
  fCaloTrueIsoRawCharged_FromDecay(),
  fCaloTrueIsoRawNeutral_FromDecay(),
  fCaloTrueIsoRawFull_FromDecay(),
  fCaloTrueIsoCell_FromDecay(),
  fCaloTrueIsoRawCharged_FromDirect(),
  fCaloTrueIsoRawNeutral_FromDirect(),
  fCaloTrueIsoRawFull_FromDirect(),
  fCaloTrueIsoCell_FromDirect(),
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
  fCaloTrueRecPtIsoCharged(),
  fCaloTrueRecPtIsoNeutral(),
  fCaloTrueRecPtIsoFull(),
  fCaloTrueRecPtIsoCell(),
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
  fRhoOutName("Rho"),
  fTreeBuffSize(60*1024*1024),
  fMemCountAOD(0),
  fTrackMatcherRunningMode(0)
{
  fDataEvtHeader.pVtxX = -9999;
  fDataEvtHeader.pVtxY = -9999;
  fDataEvtHeader.pVtxZ = -9999;
  fDataEvtHeader.runnumber = -1;
  fDataEvtHeader.numberESDtracks = -1;
  
  fMCEvtHeader.pVtxX = -9999;
  fMCEvtHeader.pVtxY = -9999;
  fMCEvtHeader.pVtxZ = -9999;
  fMCEvtHeader.runnumber = -1;
  fMCEvtHeader.numberESDtracks = -1;
  fMCEvtHeader.evtType = 0;

  SetEtaMatching(0.010,4.07,-2.5);
  SetPhiMatching(0.015,3.65,3.65);
  SetEOverP(1.75);

  SetPi0TaggingWindow(0.120,0.145);
  SetEtaTaggingWindow(0.5,0.6);

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
  fExtraClusterInfo(),
  fExtraClusterInfoBackground(),
  fDataEvtHeader(),
  fMCEvtHeader(),
  fConvIsoInfo(),
  fCaloIsoInfo(),
  fGeomEMCAL(NULL),
  fCorrTaskSetting(""),
  fEventCuts(NULL),
  fClusterCutsEMC(NULL),
  fClusterCutsIsolationEMC(NULL),
  fClusterCutsTaggingEMC(NULL),
  fClusterCutsPHOS(NULL),
  fConvCuts(NULL),
  fCaloUtils(NULL),
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
  fUseHistograms(kTRUE),
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
  fConvIsoRawCharged(),
  fConvIsoRawNeutral(),
  fConvIsoRawFull(),
  fConvIsoCell(),
  fConvIsoCorr(),
  fConvRho(NULL),
  fConvRhoTimesArea(NULL),
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
  fConvTrueIsoRawCharged(),
  fConvTrueIsoRawNeutral(),
  fConvTrueIsoRawFull(),
  fConvTrueIsoCorr(),
  fConvTrueIsoCell(),
  fConvTrueIsoRawCharged_FromDecay(),
  fConvTrueIsoRawNeutral_FromDecay(),
  fConvTrueIsoRawFull_FromDecay(),
  fConvTrueIsoCell_FromDecay(),
  fConvTrueIsoRawCharged_FromDirect(),
  fConvTrueIsoRawNeutral_FromDirect(),
  fConvTrueIsoRawFull_FromDirect(),
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
  fConvTrueRecPtIsoCharged(),
  fConvTrueRecPtIsoNeutral(),
  fConvTrueRecPtIsoFull(),
  fConvTrueRecPtIsoCell(),
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
  fCaloPt(NULL),
  fCaloPtBeforeAcc(NULL),
  fCaloPtTaggedCalo(NULL),
  fCaloPtTaggedAsDecayCalo(NULL),
  fCaloIsoRawCharged(),
  fCaloIsoRawNeutral(),
  fCaloIsoRawFull(),
  fCaloIsoCell(),
  fCaloIsoCorr(),
  fCaloRho(NULL),
  fCaloRhoTimesArea(NULL),
  fCaloTruePt(NULL),
  fCaloTruePtPrimary(NULL),
  fCaloTruePtDecay(NULL),
  fCaloTruePtDecayFoundOtherInCluster(NULL),
  fCaloTruePtDecayOtherInAcc(NULL),
  fCaloTruePtDecayOtherInAccAboveMinEnergy(NULL),
  fCaloTruePtTaggedCalo(NULL),
  fCaloTruePtTaggedAsDecayCalo(NULL),
  fCaloTrueRecPt(NULL),
  fCaloTrueRecPtPrimary(NULL),
  fCaloTrueRecPtDecay(NULL),
  fCaloTrueRecPtDecayFoundOtherInCluster(NULL),
  fCaloTrueRecPtDecayOtherInAcc(NULL),
  fCaloTrueRecPtDecayOtherInAccAboveMinEnergy(NULL),
  fCaloTrueRecPtTaggedCalo(NULL),
  fCaloTrueRecPtTaggedAsDecayCalo(NULL),
  fCaloTrueIsoRawCharged(),
  fCaloTrueIsoRawNeutral(),
  fCaloTrueIsoRawFull(),
  fCaloTrueIsoCorr(),
  fCaloTrueIsoCell(),
  fCaloTrueIsoRawCharged_FromDecay(),
  fCaloTrueIsoRawNeutral_FromDecay(),
  fCaloTrueIsoRawFull_FromDecay(),
  fCaloTrueIsoCell_FromDecay(),
  fCaloTrueIsoRawCharged_FromDirect(),
  fCaloTrueIsoRawNeutral_FromDirect(),
  fCaloTrueIsoRawFull_FromDirect(),
  fCaloTrueIsoCell_FromDirect(),
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
  fCaloTrueRecPtIsoCharged(),
  fCaloTrueRecPtIsoNeutral(),
  fCaloTrueRecPtIsoFull(),
  fCaloTrueRecPtIsoCell(),
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
  fRhoOutName("Rho"),
  fTreeBuffSize(60*1024*1024),
  fMemCountAOD(0),
  fTrackMatcherRunningMode(0)
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  
  fDataEvtHeader.pVtxX = -9999;
  fDataEvtHeader.pVtxY = -9999;
  fDataEvtHeader.pVtxZ = -9999;
  fDataEvtHeader.runnumber = -1;
  fDataEvtHeader.numberESDtracks = -1;
  
  fMCEvtHeader.pVtxX = -9999;
  fMCEvtHeader.pVtxY = -9999;
  fMCEvtHeader.pVtxZ = -9999;
  fMCEvtHeader.runnumber = -1;
  fMCEvtHeader.numberESDtracks = -1;
  fMCEvtHeader.evtType = 0;

  SetEtaMatching(0.010,4.07,-2.5);
  SetPhiMatching(0.015,3.65,3.65);
  SetEOverP(1.75);

  SetPi0TaggingWindow(0.120,0.145);
  SetEtaTaggingWindow(0.5,0.6);
}

//________________________________________________________________________
AliAnalysisTaskGammaIsoTree::~AliAnalysisTaskGammaIsoTree()
{
  // default deconstructor
  if(fConversionCandidates){
    delete fConversionCandidates;
    fConversionCandidates = 0x0;
  }
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
  if(fTracks){
    delete fTracks;
    fTracks = 0x0;
  }
  if(fMCParticles){
    delete fMCParticles;
    fMCParticles = 0x0;
  }
  if(fMCParticles){
    delete fMCParticles;
    fMCParticles = 0x0;
  }
  if(fMCParticles){
    delete fMCParticles;
    fMCParticles = 0x0;
  }
  if(fExtraClusterInfo){
    delete fExtraClusterInfo;
    fExtraClusterInfo = 0x0;
  }
  if(fExtraClusterInfoBackground){
    delete fExtraClusterInfoBackground;
    fExtraClusterInfoBackground = 0x0;
  }
  if(fConvIsoInfo){
    delete fConvIsoInfo;
    fConvIsoInfo = 0x0;
  }
  if(fCaloIsoInfo){
    delete fCaloIsoInfo;
    fCaloIsoInfo = 0x0;
  }
}
//________________________________________________________________________
void AliAnalysisTaskGammaIsoTree::UserCreateOutputObjects()
{
  // Create User Output Objects
  fOutputList                         = new TList();
  fOutputList->SetOwner(kTRUE);

  if(((AliConvEventCuts*)fEventCuts)->GetCutHistograms()){
    fOutputList->Add(((AliConvEventCuts*)fEventCuts)->GetCutHistograms());
  }

  if(((AliCaloPhotonCuts*)fClusterCutsEMC)->GetCutHistograms()){
    fOutputList->Add(((AliCaloPhotonCuts*)fClusterCutsEMC)->GetCutHistograms());
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
    for(Int_t iMatcherTask = 0; iMatcherTask < 5; iMatcherTask++){
      AliCaloTrackMatcher* temp = 0x0;
      if(!fCorrTaskSetting.CompareTo("")){
        temp = (AliCaloTrackMatcher*) (AliAnalysisManager::GetAnalysisManager()->GetTask(Form("CaloTrackMatcherSignal_%i_%i",iMatcherTask,fTrackMatcherRunningMode)));
      } else {
        temp = (AliCaloTrackMatcher*) (AliAnalysisManager::GetAnalysisManager()->GetTask(Form("CaloTrackMatcherSignal_%i_%i_%s",iMatcherTask,fTrackMatcherRunningMode,fCorrTaskSetting.Data())));
      }
      if(temp) fOutputList->Add(temp->GetCaloTrackMatcherHistograms());
    }
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
  fHistoNEvents->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
  fHistoNEvents->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
  fHistoNEvents->GetXaxis()->SetBinLabel(14,"Pileup V0M-TPCout Tracks");
  fHistoNEvents->GetYaxis()->SetTitle("N_{events}");
  fGeneralFolder->Add(fHistoNEvents);

  fHistoChargedIso           = new TH1F("fHistoChargedIso","fHistoChargedIso",500,-0.5,50);
  fOutputList->Add(fHistoChargedIso);

  fHistoTaggingPCMPCM           = new TH2F("fHistoTaggingPCMPCM","fHistoTaggingPCMPCM;M (GeV/c^2); photon p_{T} (GeV/c)",500,0.,1.,100,0,50.);
  fGeneralFolder->Add(fHistoTaggingPCMPCM);
  fHistoTaggingPCMEMC           = new TH2F("fHistoTaggingPCMEMC","fHistoTaggingPCMEMC;M (GeV/c^2); photon p_{T} (GeV/c)",500,0.,1.,100,0,50.);
  fGeneralFolder->Add(fHistoTaggingPCMEMC);
  fHistoTaggingEMCPCM           = new TH2F("fHistoTaggingEMCPCM","fHistoTaggingEMCPCM;M (GeV/c^2); photon p_{T} (GeV/c)",500,0.,1.,100,0,50.);
  fGeneralFolder->Add(fHistoTaggingEMCPCM);
  fHistoTaggingEMCEMC           = new TH2F("fHistoTaggingEMCEMC","fHistoTaggingEMCEMC;M (GeV/c^2); photon p_{T} (GeV/c)",500,0.,1.,100,0,50.);
  fGeneralFolder->Add(fHistoTaggingEMCEMC);


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
    fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(10,"EMCAL problem");
    fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
    fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
    fHistoNEventsWOWeight->GetXaxis()->SetBinLabel(14,"Pileup V0M-TPCout Tracks");
    fHistoNEventsWOWeight->GetYaxis()->SetTitle("N_{events}");
    fGeneralFolder->Add(fHistoNEventsWOWeight);
  }

  //
  // ─── CONVERSION HISOGRAMS ───────────────────────────────────────────────────────
  //
  Double_t minPt = 0;
  Double_t maxPt = 30;
  Int_t  nPtBins = 600;

  fConvFolderRec          = new TList();
  fConvFolderRec->SetName("convPhotonsRec");
  fConvFolderRec->SetOwner(kTRUE);
  fOutputList->Add(fConvFolderRec);

  if(fUseHistograms){

    fConvPt = new TH1F("fConvPt","conversion photons in EMC acc;p_{T} (GeV/c); counts",nPtBins,minPt,maxPt);
    fConvPtBeforeAcc = new TH1F("fConvPtBeforeAcc", "conversion photons all acc;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);

    fConvPtTaggedCalo = new TH1F("fConvPtTaggedCalo", "conversion photons that survived tagging;p_{T} (GeV/c); counts", 500, 0, 50.);
    fConvPtTaggedAsDecayCalo = new TH1F("fConvPtTaggedAsDecayCalo", "conversion photons that survived tagging;p_{T} (GeV/c); counts", 500, 0, 50.);


    fConvRho = new TH1F("fConvRho", "charged event density;#rho; counts", 500, 0, 50.);
    fConvRhoTimesArea = new TH1F("fConvRhoTimesArea", "charged event density;#rho #times jet Area; counts", 500, 0, 50.);
    
    fConvFolderRec->Add(fConvPt);
    fConvFolderRec->Add(fConvPtBeforeAcc);
    fConvFolderRec->Add(fConvPtTaggedCalo);
    fConvFolderRec->Add(fConvPtTaggedAsDecayCalo);
    fConvFolderRec->Add(fConvRho);
    fConvFolderRec->Add(fConvRhoTimesArea);

    for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
    {
      TH2F *convIsoRawCharged = new TH2F(Form("convIsoRawCharged%i",r), Form("charged track ISO in R < %1.1f;#sum p_{T} (GeV/c); conv. p_{T} (GeV/c)",fTrackIsolationR.at(r)), nPtBins,minPt,maxPt,100,0,50);
      fConvIsoRawCharged[r] = (TH2F*) convIsoRawCharged->Clone(Form("fConvIsoRawCharged_R%1.1f",fTrackIsolationR.at(r)));
      fConvFolderRec->Add(fConvIsoRawCharged[r]);
    }
    for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
    {
      TH2F *convIsoRawNeutral = new TH2F(Form("fConvIsoRawNeutral_%i",r), Form("Neutral track ISO in R < %1.1f;#sum p_{T} (GeV/c); conv. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,100,0,50);
      fConvIsoRawNeutral[r] = (TH2F*) convIsoRawNeutral->Clone(Form("fConvIsoRawNeutral_R%1.1f",fNeutralIsolationR.at(r)));
      fConvFolderRec->Add(fConvIsoRawNeutral[r]);

      TH2F *convIsoRawFull = new TH2F(Form("fConvIsoRawFull_%i",r), Form("Full track ISO in R < %1.1f;#sum p_{T} (GeV/c); conv. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,100,0,50);
      fConvIsoRawFull[r] = (TH2F*) convIsoRawFull->Clone(Form("fConvIsoRawFull_R%1.1f",fNeutralIsolationR.at(r)));
      fConvFolderRec->Add(fConvIsoRawFull[r]);

      TH2F *convIsoCell = new TH2F(Form("fConvIsoCell_%i",r), Form("Full track ISO in R < %1.1f;#sum p_{T} (GeV/c); conv. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,100,0,50);
      fConvIsoCell[r] = (TH2F*) convIsoCell->Clone(Form("fConvIsoCell_R%1.1f",fNeutralIsolationR.at(r)));
      fConvFolderRec->Add(fConvIsoCell[r]);
    }

    fConvFolderTrue          = new TList();
    fConvFolderTrue->SetName("convPhotonsTrue");
    fConvFolderTrue->SetOwner(kTRUE);
    if(fIsMC > 0){
      fOutputList->Add(fConvFolderTrue);
      fConvTruePt = new TH1F("fConvTruePt", "validated conversion photons in EMC acceptance;p_{T} (GeV/c); counts", 500, 0, 50.);
      fConvTruePtPrimary = new TH1F("fConvTruePtPrimary", "conversion photon that has not a pi0 etc. as mother;p_{T} (GeV/c); counts", 500, 0, 50.);
      fConvTruePtDecay = new TH1F("fConvTruePtDecay", "conversion photon from decay;p_{T} (GeV/c); counts", 500, 0, 50.);
      fConvTruePtDecayFoundOtherInCluster = new TH1F("fConvTruePtDecayFoundOtherInCluster", "conversion photon from decay, where the other decay particle was found in EMC;p_{T} (GeV/c); counts", 500, 0, 50.);
      fConvTruePtDecayOtherInAcc = new TH1F("fConvTruePtDecayOtherInAcc", "conversion photon from decay, where the other decay particle was found;p_{T} (GeV/c); counts", 500, 0, 50.);
      fConvTruePtDecayOtherInAccAboveMinEnergy = new TH1F("fConvTruePtDecayOtherInAccAboveMinEnergy", "conversion photon from decay, where the other decay particle is in acc. and above 0.7 GeV;p_{T} (GeV/c); counts", 500, 0, 50.);
  
      // iso and tagging studies
      fConvTruePtTaggedCalo = new TH1F("fConvTruePtTaggedCalo", "conversion photons that survived tagging;p_{T} (GeV/c); counts", 500, 0, 50.);
      fConvTruePtTaggedAsDecayCalo = new TH1F("fConvTruePtTaggedAsDecayCalo", "conversion photons that survived tagging;p_{T} (GeV/c); counts", 500, 0, 50.);

      // True with rec pT
      fConvTrueRecPt = new TH1F("fConvTrueRecPt", "validated conversion photons in EMC acceptance;p_{T} (GeV/c); counts", 500, 0, 50.);
      fConvTrueRecPtPrimary = new TH1F("fConvTrueRecPtPrimary", "conversion photon that has not a pi0 etc. as mother;p_{T} (GeV/c); counts", 500, 0, 50.);
      fConvTrueRecPtDecay = new TH1F("fConvTrueRecPtDecay", "conversion photon from decay;p_{T} (GeV/c); counts", 500, 0, 50.);
      fConvTrueRecPtDecayFoundOtherInCluster = new TH1F("fConvTrueRecPtDecayFoundOtherInCluster", "conversion photon from decay, where the other decay particle was found in EMC;p_{T} (GeV/c); counts", 500, 0, 50.);
      fConvTrueRecPtDecayOtherInAcc = new TH1F("fConvTrueRecPtDecayOtherInAcc", "conversion photon from decay, where the other decay particle was found;p_{T} (GeV/c); counts", 500, 0, 50.);
      fConvTrueRecPtDecayOtherInAccAboveMinEnergy = new TH1F("fConvTrueRecPtDecayOtherInAccAboveMinEnergy", "conversion photon from decay, where the other decay particle is in acc. and above 0.7 GeV;p_{T} (GeV/c); counts", 500, 0, 50.);

      // iso and tagging studies
      fConvTrueRecPtTaggedCalo = new TH1F("fConvTrueRecPtTaggedCalo", "conversion photons that survived tagging;p_{T} (GeV/c); counts", 500, 0, 50.);
      fConvTrueRecPtTaggedAsDecayCalo = new TH1F("fConvTrueRecPtTaggedAsDecayCalo", "conversion photons that survived tagging;p_{T} (GeV/c); counts", 500, 0, 50.);

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
        TH2F *convTrueIsoRawCharged = new TH2F(Form("fConvTrueIsoRawCharged_%i",r), Form("Charged track ISO in R < %1.1f;#sum p_{T} (GeV/c); counts",fTrackIsolationR.at(r)), nPtBins,minPt,maxPt,100,0,50);
        fConvTrueIsoRawCharged[r] = (TH2F*) convTrueIsoRawCharged->Clone(Form("fConvTrueIsoRawCharged_R%1.1f",fTrackIsolationR.at(r)));
        fConvFolderTrue->Add(fConvTrueIsoRawCharged[r]); 

        TH2F *convTrueIsoRawCharged_FromDecay = new TH2F(Form("fConvTrueIsoRawCharged_FromDecay_%i",r), Form("Charged track ISO in R < %1.1f;#sum p_{T} (GeV/c); counts",fTrackIsolationR.at(r)), nPtBins,minPt,maxPt,100,0,50);
        fConvTrueIsoRawCharged_FromDecay[r] = (TH2F*) convTrueIsoRawCharged_FromDecay->Clone(Form("fConvTrueIsoRawCharged_FromDecay_R%1.1f",fTrackIsolationR.at(r)));
        fConvFolderTrue->Add(fConvTrueIsoRawCharged_FromDecay[r]); 

        TH2F *convTrueIsoRawCharged_FromDirect = new TH2F(Form("fConvTrueIsoRawCharged_FromDirect_%i",r), Form("Charged track ISO in R < %1.1f;#sum p_{T} (GeV/c); counts",fTrackIsolationR.at(r)), nPtBins,minPt,maxPt,100,0,50);
        fConvTrueIsoRawCharged_FromDirect[r] = (TH2F*) convTrueIsoRawCharged_FromDirect->Clone(Form("fConvTrueIsoRawCharged_FromDirect_R%1.1f",fTrackIsolationR.at(r)));
        fConvFolderTrue->Add(fConvTrueIsoRawCharged_FromDirect[r]); 
      }
      
      for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
      {
        TH2F *convTrueIsoRawNeutral = new TH2F(Form("fConvTrueIsoRawNeutral_%i",r), Form("Neutral track ISO in R < %1.1f;#sum p_{T} (GeV/c); conv. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,100,0,50);
        fConvTrueIsoRawNeutral[r] = (TH2F*) convTrueIsoRawNeutral->Clone(Form("fConvTrueIsoRawNeutral_R%1.1f",fNeutralIsolationR.at(r)));
        fConvFolderTrue->Add(fConvTrueIsoRawNeutral[r]);

        TH2F *convTrueIsoRawFull = new TH2F(Form("fConvTrueIsoRawFull_%i",r), Form("Full track ISO in R < %1.1f;#sum p_{T} (GeV/c); conv. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,100,0,50);
        fConvTrueIsoRawFull[r] = (TH2F*) convTrueIsoRawFull->Clone(Form("fConvTrueIsoRawFull_R%1.1f",fNeutralIsolationR.at(r)));
        fConvFolderTrue->Add(fConvTrueIsoRawFull[r]);

        TH2F *convTrueIsoCell = new TH2F(Form("fConvTrueIsoCell_%i",r), Form("Full track ISO in R < %1.1f;#sum p_{T} (GeV/c); conv. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,100,0,50);
        fConvTrueIsoCell[r] = (TH2F*) convTrueIsoCell->Clone(Form("fConvTrueIsoCell_R%1.1f",fNeutralIsolationR.at(r)));
        fConvFolderTrue->Add(fConvTrueIsoCell[r]);

        TH2F *convTrueIsoRawNeutral_FromDecay = new TH2F(Form("fConvTrueIsoRawNeutral_FromDecay_%i",r), Form("Neutral track ISO in R < %1.1f;#sum p_{T} (GeV/c); conv. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,100,0,50);
        fConvTrueIsoRawNeutral_FromDecay[r] = (TH2F*) convTrueIsoRawNeutral_FromDecay->Clone(Form("fConvTrueIsoRawNeutral_FromDecay_R%1.1f",fNeutralIsolationR.at(r)));
        fConvFolderTrue->Add(fConvTrueIsoRawNeutral_FromDecay[r]);

        TH2F *convTrueIsoRawFull_FromDecay = new TH2F(Form("fConvTrueIsoRawFull_FromDecay_%i",r), Form("Full track ISO in R < %1.1f;#sum p_{T} (GeV/c); conv. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,100,0,50);
        fConvTrueIsoRawFull_FromDecay[r] = (TH2F*) convTrueIsoRawFull_FromDecay->Clone(Form("fConvTrueIsoRawFull_FromDecay_R%1.1f",fNeutralIsolationR.at(r)));
        fConvFolderTrue->Add(fConvTrueIsoRawFull_FromDecay[r]);

        TH2F *convTrueIsoCell_FromDecay = new TH2F(Form("fConvTrueIsoCell_FromDecay_%i",r), Form("Full track ISO in R < %1.1f;#sum p_{T} (GeV/c); conv. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,100,0,50);
        fConvTrueIsoCell_FromDecay[r] = (TH2F*) convTrueIsoCell_FromDecay->Clone(Form("fConvTrueIsoCell_FromDecay_R%1.1f",fNeutralIsolationR.at(r)));
        fConvFolderTrue->Add(fConvTrueIsoCell_FromDecay[r]);

        TH2F *convTrueIsoRawNeutral_FromDirect = new TH2F(Form("fConvTrueIsoRawNeutral_FromDirect_%i",r), Form("Neutral track ISO in R < %1.1f;#sum p_{T} (GeV/c); conv. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,100,0,50);
        fConvTrueIsoRawNeutral_FromDirect[r] = (TH2F*) convTrueIsoRawNeutral_FromDirect->Clone(Form("fConvTrueIsoRawNeutral_FromDirect_R%1.1f",fNeutralIsolationR.at(r)));
        fConvFolderTrue->Add(fConvTrueIsoRawNeutral_FromDirect[r]);

        TH2F *convTrueIsoRawFull_FromDirect = new TH2F(Form("fConvTrueIsoRawFull_FromDirect_%i",r), Form("Full track ISO in R < %1.1f;#sum p_{T} (GeV/c); conv. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,100,0,50);
        fConvTrueIsoRawFull_FromDirect[r] = (TH2F*) convTrueIsoRawFull_FromDirect->Clone(Form("fConvTrueIsoRawFull_FromDirect_R%1.1f",fNeutralIsolationR.at(r)));
        fConvFolderTrue->Add(fConvTrueIsoRawFull_FromDirect[r]);

        TH2F *convTrueIsoCell_FromDirect = new TH2F(Form("fConvTrueIsoCell_FromDirect_%i",r), Form("Full track ISO in R < %1.1f;#sum p_{T} (GeV/c); conv. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,100,0,50);
        fConvTrueIsoCell_FromDirect[r] = (TH2F*) convTrueIsoCell_FromDirect->Clone(Form("fConvTrueIsoCell_FromDirect_R%1.1f",fNeutralIsolationR.at(r)));
        fConvFolderTrue->Add(fConvTrueIsoCell_FromDirect[r]);
      }
    }
    for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
    {
      for (UInt_t e = 0; e < fTrackIsolationE.size(); e++)
      {
        TH1F *convPtIsoCharged = new TH1F(Form("fConvPtIsoCharged_%i_%i",r,e), Form("conversion photons with charged track ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
        fConvPtIsoCharged[r][e] = (TH1F*) convPtIsoCharged->Clone(Form("fConvPtIsoCharged_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
        fConvFolderRec->Add(fConvPtIsoCharged[r][e]);
        TH1F *convPtTaggedCaloIsoCharged = new TH1F(Form("fConvPtTaggedCaloIsoCharged_%i_%i",r,e), Form("conversion photons with charged track ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
        fConvPtTaggedCaloIsoCharged[r][e] = (TH1F*) convPtTaggedCaloIsoCharged->Clone(Form("fConvPtTaggedCaloIsoCharged_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
        fConvFolderRec->Add(fConvPtTaggedCaloIsoCharged[r][e]);

        if(fIsMC>0){
          // true pT
          TH1F *convTruePtIsoCharged = new TH1F(Form("fConvTruePtIsoCharged_%i_%i",r,e), Form("conversion photons with charged track ISO < %1.1f GeV in R < %1.1f;gen. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtIsoCharged[r][e] = (TH1F*) convTruePtIsoCharged->Clone(Form("fConvTruePtIsoCharged_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fConvFolderTrue->Add(fConvTruePtIsoCharged[r][e]);
          TH1F *convTruePtTaggedCaloIsoCharged = new TH1F(Form("fConvTruePtTaggedCaloIsoCharged_%i_%i",r,e), Form("conversion photons with charged track ISO < %1.1f GeV in R < %1.1f + tagging;gen. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtTaggedCaloIsoCharged[r][e] = (TH1F*) convTruePtTaggedCaloIsoCharged->Clone(Form("fConvTruePtTaggedCaloIsoCharged_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fConvFolderTrue->Add(fConvTruePtTaggedCaloIsoCharged[r][e]);

          //mc iso
          TH1F *convTruePtMCIsoCharged = new TH1F(Form("fConvTruePtMCIsoCharged_%i_%i",r,e), Form("conversion photons with charged track MCIso < %1.1f GeV in R < %1.1f;gen. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtMCIsoCharged[r][e] = (TH1F*) convTruePtMCIsoCharged->Clone(Form("fConvTruePtMCIsoCharged_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fConvFolderTrue->Add(fConvTruePtMCIsoCharged[r][e]);
          TH1F *convTruePtTaggedCaloMCIsoCharged = new TH1F(Form("fConvTruePtTaggedCaloMCIsoCharged_%i_%i",r,e), Form("conversion photons with charged track MCIso < %1.1f GeV in R < %1.1f + tagging;gen. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtTaggedCaloMCIsoCharged[r][e] = (TH1F*) convTruePtTaggedCaloMCIsoCharged->Clone(Form("fConvTruePtTaggedCaloMCIsoCharged_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fConvFolderTrue->Add(fConvTruePtTaggedCaloMCIsoCharged[r][e]);

          // rec Pt
          TH1F *convTrueRecPtIsoCharged = new TH1F(Form("fConvTrueRecPtIsoCharged_%i_%i",r,e), Form("conversion photons with charged track ISO < %1.1f GeV in R < %1.1f;rec. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtIsoCharged[r][e] = (TH1F*) convTrueRecPtIsoCharged->Clone(Form("fConvTrueRecPtIsoCharged_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fConvFolderTrue->Add(fConvTrueRecPtIsoCharged[r][e]);
          TH1F *convTrueRecPtTaggedCaloIsoCharged = new TH1F(Form("fConvTrueRecPtTaggedCaloIsoCharged_%i_%i",r,e), Form("conversion photons with charged track ISO < %1.1f GeV in R < %1.1f + tagging;rec. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtTaggedCaloIsoCharged[r][e] = (TH1F*) convTrueRecPtTaggedCaloIsoCharged->Clone(Form("fConvTrueRecPtTaggedCaloIsoCharged_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fConvFolderTrue->Add(fConvTrueRecPtTaggedCaloIsoCharged[r][e]);

          TH1F *convTrueRecPtMCIsoCharged = new TH1F(Form("fConvTrueRecPtMCIsoCharged_%i_%i",r,e), Form("conversion photons with charged track MCIso < %1.1f GeV in R < %1.1f;rec. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtMCIsoCharged[r][e] = (TH1F*) convTrueRecPtMCIsoCharged->Clone(Form("fConvTrueRecPtMCIsoCharged_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
          fConvFolderTrue->Add(fConvTrueRecPtMCIsoCharged[r][e]);
          TH1F *convTrueRecPtTaggedCaloMCIsoCharged = new TH1F(Form("fConvTrueRecPtTaggedCaloMCIsoCharged_%i_%i",r,e), Form("conversion photons with charged track MCIso < %1.1f GeV in R < %1.1f + tagging;rec. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtTaggedCaloMCIsoCharged[r][e] = (TH1F*) convTrueRecPtTaggedCaloMCIsoCharged->Clone(Form("fConvTrueRecPtTaggedCaloMCIsoCharged_R%1.1f_E%1.1f",fTrackIsolationR.at(r),fTrackIsolationE.at(e)));
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
        fConvFolderRec->Add(fConvPtIsoNeutral[r][e]);
        TH1F *convPtIsoFull = new TH1F(Form("fConvPtIsoFull_%i_%i",r,e), Form("conversion photons with Full ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
        fConvPtIsoFull[r][e] = (TH1F*) convPtIsoFull->Clone(Form("fConvPtIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
        fConvFolderRec->Add(fConvPtIsoFull[r][e]);
        TH1F *convPtIsoCell = new TH1F(Form("fConvPtIsoCell_%i_%i",r,e), Form("conversion photons with Cell ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
        fConvPtIsoCell[r][e] = (TH1F*) convPtIsoCell->Clone(Form("fConvPtIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
        fConvFolderRec->Add(fConvPtIsoCell[r][e]);

        TH1F *convPtTaggedCaloIsoNeutral = new TH1F(Form("fConvPtTaggedCaloIsoNeutral_%i_%i",r,e), Form("conversion photons with Neutral ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
        fConvPtTaggedCaloIsoNeutral[r][e] = (TH1F*) convPtTaggedCaloIsoNeutral->Clone(Form("fConvPtTaggedCaloIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
        fConvFolderRec->Add(fConvPtTaggedCaloIsoNeutral[r][e]);
        TH1F *convPtTaggedCaloIsoFull = new TH1F(Form("fConvTaggedCaloIsoFull_%i_%i",r,e), Form("conversion photons with Full ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
        fConvPtTaggedCaloIsoFull[r][e] = (TH1F*) convPtTaggedCaloIsoFull->Clone(Form("fConvTaggedCaloIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
        fConvFolderRec->Add(fConvPtTaggedCaloIsoFull[r][e]);
        TH1F *convPtTaggedCaloIsoCell = new TH1F(Form("fConvTaggedCaloIsoCell_%i_%i",r,e), Form("conversion photons with Cell ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
        fConvPtTaggedCaloIsoCell[r][e] = (TH1F*) convPtTaggedCaloIsoCell->Clone(Form("fConvTaggedCaloIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
        fConvFolderRec->Add(fConvPtTaggedCaloIsoCell[r][e]);

        // True Pt
        if(fIsMC>0){
          TH1F *convTruePtIsoNeutral = new TH1F(Form("fConvTruePtIsoNeutral_%i_%i",r,e), Form("conversion photons with Neutral ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtIsoNeutral[r][e] = (TH1F*) convTruePtIsoNeutral->Clone(Form("fConvTruePtIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvFolderTrue->Add(fConvTruePtIsoNeutral[r][e]);
          TH1F *convTruePtIsoFull = new TH1F(Form("fConvTruePtIsoFull_%i_%i",r,e), Form("conversion photons with Full ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtIsoFull[r][e] = (TH1F*) convTruePtIsoFull->Clone(Form("fConvTruePtIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvFolderTrue->Add(fConvTruePtIsoFull[r][e]);
          TH1F *convTruePtIsoCell = new TH1F(Form("fConvTruePtIsoCell_%i_%i",r,e), Form("conversion photons with Cell ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtIsoCell[r][e] = (TH1F*) convTruePtIsoCell->Clone(Form("fConvTruePtIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvFolderTrue->Add(fConvTruePtIsoCell[r][e]);

          TH1F *convTruePtTaggedCaloIsoNeutral = new TH1F(Form("fConvTruePtTaggedCaloIsoNeutral_%i_%i",r,e), Form("conversion photons with Neutral ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtTaggedCaloIsoNeutral[r][e] = (TH1F*) convTruePtTaggedCaloIsoNeutral->Clone(Form("fConvTruePtTaggedCaloIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvFolderTrue->Add(fConvTruePtTaggedCaloIsoNeutral[r][e]);
          TH1F *convTruePtTaggedCaloIsoFull = new TH1F(Form("fConvTruePtTaggedCaloIsoFull_%i_%i",r,e), Form("conversion photons with Full ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtTaggedCaloIsoFull[r][e] = (TH1F*) convTruePtTaggedCaloIsoFull->Clone(Form("fConvTruePtTaggedCaloIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvFolderTrue->Add(fConvTruePtTaggedCaloIsoFull[r][e]);
          TH1F *convTruePtTaggedCaloIsoCell = new TH1F(Form("fConvTruePtTaggedCaloIsoCell_%i_%i",r,e), Form("conversion photons with Cell ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtTaggedCaloIsoCell[r][e] = (TH1F*) convTruePtTaggedCaloIsoCell->Clone(Form("fConvTruePtTaggedCaloIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvFolderTrue->Add(fConvTruePtTaggedCaloIsoCell[r][e]);

          TH1F *convTruePtMCIsoNeutral = new TH1F(Form("fConvTruePtMCIsoNeutral_%i_%i",r,e), Form("conversion photons with Neutral MCIso < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtMCIsoNeutral[r][e] = (TH1F*) convTruePtMCIsoNeutral->Clone(Form("fConvTruePtMCIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvFolderTrue->Add(fConvTruePtMCIsoNeutral[r][e]);
          TH1F *convTruePtMCIsoFull = new TH1F(Form("fConvTruePtMCIsoFull_%i_%i",r,e), Form("conversion photons with Full MCIso < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtMCIsoFull[r][e] = (TH1F*) convTruePtMCIsoFull->Clone(Form("fConvTruePtMCIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvFolderTrue->Add(fConvTruePtMCIsoFull[r][e]);
          TH1F *convTruePtMCIsoCell = new TH1F(Form("fConvTruePtMCIsoCell_%i_%i",r,e), Form("conversion photons with Cell MCIso < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtMCIsoCell[r][e] = (TH1F*) convTruePtMCIsoCell->Clone(Form("fConvTruePtMCIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvFolderTrue->Add(fConvTruePtMCIsoCell[r][e]);

          TH1F *convTruePtTaggedCaloMCIsoNeutral = new TH1F(Form("fConvTruePtTaggedCaloMCIsoNeutral_%i_%i",r,e), Form("conversion photons with Neutral MCIso < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtTaggedCaloMCIsoNeutral[r][e] = (TH1F*) convTruePtTaggedCaloMCIsoNeutral->Clone(Form("fConvTruePtTaggedCaloMCIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvFolderTrue->Add(fConvTruePtTaggedCaloMCIsoNeutral[r][e]);
          TH1F *convTruePtTaggedCaloMCIsoFull = new TH1F(Form("fConvTruePtTaggedCaloMCIsoFull_%i_%i",r,e), Form("conversion photons with Full MCIso < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtTaggedCaloMCIsoFull[r][e] = (TH1F*) convTruePtTaggedCaloMCIsoFull->Clone(Form("fConvTruePtTaggedCaloMCIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvFolderTrue->Add(fConvTruePtTaggedCaloMCIsoFull[r][e]);
          TH1F *convTruePtTaggedCaloMCIsoCell = new TH1F(Form("fConvTruePtTaggedCaloMCIsoCell_%i_%i",r,e), Form("conversion photons with Cell MCIso < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTruePtTaggedCaloMCIsoCell[r][e] = (TH1F*) convTruePtTaggedCaloMCIsoCell->Clone(Form("fConvTruePtTaggedCaloMCIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvFolderTrue->Add(fConvTruePtTaggedCaloMCIsoCell[r][e]);

          // rec Pt


          TH1F *convTrueRecPtIsoNeutral = new TH1F(Form("fConvTrueRecPtIsoNeutral_%i_%i",r,e), Form("conversion photons with Neutral ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtIsoNeutral[r][e] = (TH1F*) convTrueRecPtIsoNeutral->Clone(Form("fConvTrueRecPtIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvFolderTrue->Add(fConvTrueRecPtIsoNeutral[r][e]);
          TH1F *convTrueRecPtIsoFull = new TH1F(Form("fConvTrueRecPtIsoFull_%i_%i",r,e), Form("conversion photons with Full ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtIsoFull[r][e] = (TH1F*) convTrueRecPtIsoFull->Clone(Form("fConvTrueRecPtIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvFolderTrue->Add(fConvTrueRecPtIsoFull[r][e]);
          TH1F *convTrueRecPtIsoCell = new TH1F(Form("fConvTrueRecPtIsoCell_%i_%i",r,e), Form("conversion photons with Cell ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtIsoCell[r][e] = (TH1F*) convTrueRecPtIsoCell->Clone(Form("fConvTrueRecPtIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvFolderTrue->Add(fConvTrueRecPtIsoCell[r][e]);

          TH1F *convTrueRecPtTaggedCaloIsoNeutral = new TH1F(Form("fConvTrueRecPtTaggedCaloIsoNeutral_%i_%i",r,e), Form("conversion photons with Neutral ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtTaggedCaloIsoNeutral[r][e] = (TH1F*) convTrueRecPtTaggedCaloIsoNeutral->Clone(Form("fConvTrueRecPtTaggedCaloIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvFolderTrue->Add(fConvTrueRecPtTaggedCaloIsoNeutral[r][e]);
          TH1F *convTrueRecPtTaggedCaloIsoFull = new TH1F(Form("fConvTrueRecPtTaggedCaloIsoFull_%i_%i",r,e), Form("conversion photons with Full ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtTaggedCaloIsoFull[r][e] = (TH1F*) convTrueRecPtTaggedCaloIsoFull->Clone(Form("fConvTrueRecPtTaggedCaloIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvFolderTrue->Add(fConvTrueRecPtTaggedCaloIsoFull[r][e]);
          TH1F *convTrueRecPtTaggedCaloIsoCell = new TH1F(Form("fConvTrueRecPtTaggedCaloIsoCell_%i_%i",r,e), Form("conversion photons with Cell ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtTaggedCaloIsoCell[r][e] = (TH1F*) convTrueRecPtTaggedCaloIsoCell->Clone(Form("fConvTrueRecPtTaggedCaloIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvFolderTrue->Add(fConvTrueRecPtTaggedCaloIsoCell[r][e]);

          TH1F *convTrueRecPtMCIsoNeutral = new TH1F(Form("fConvTrueRecPtMCIsoNeutral_%i_%i",r,e), Form("conversion photons with Neutral ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtMCIsoNeutral[r][e] = (TH1F*) convTrueRecPtMCIsoNeutral->Clone(Form("fConvTrueRecPtMCIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvFolderTrue->Add(fConvTrueRecPtMCIsoNeutral[r][e]);
          TH1F *convTrueRecPtMCIsoFull = new TH1F(Form("fConvTrueRecPtMCIsoFull_%i_%i",r,e), Form("conversion photons with Full ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtMCIsoFull[r][e] = (TH1F*) convTrueRecPtMCIsoFull->Clone(Form("fConvTrueRecPtMCIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvFolderTrue->Add(fConvTrueRecPtMCIsoFull[r][e]);
          TH1F *convTrueRecPtMCIsoCell = new TH1F(Form("fConvTrueRecPtMCIsoCell_%i_%i",r,e), Form("conversion photons with Cell ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtMCIsoCell[r][e] = (TH1F*) convTrueRecPtMCIsoCell->Clone(Form("fConvTrueRecPtMCIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvFolderTrue->Add(fConvTrueRecPtMCIsoCell[r][e]);

          TH1F *convTrueRecPtTaggedCaloMCIsoNeutral = new TH1F(Form("fConvTrueRecPtTaggedCaloMCIsoNeutral_%i_%i",r,e), Form("conversion photons with Neutral ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtTaggedCaloMCIsoNeutral[r][e] = (TH1F*) convTrueRecPtTaggedCaloMCIsoNeutral->Clone(Form("fConvTrueRecPtTaggedCaloMCIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvFolderTrue->Add(fConvTrueRecPtTaggedCaloMCIsoNeutral[r][e]);
          TH1F *convTrueRecPtTaggedCaloMCIsoFull = new TH1F(Form("fConvTrueRecPtTaggedCaloMCIsoFull_%i_%i",r,e), Form("conversion photons with Full ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtTaggedCaloMCIsoFull[r][e] = (TH1F*) convTrueRecPtTaggedCaloMCIsoFull->Clone(Form("fConvTrueRecPtTaggedCaloMCIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvFolderTrue->Add(fConvTrueRecPtTaggedCaloMCIsoFull[r][e]);
          TH1F *convTrueRecPtTaggedCaloMCIsoCell = new TH1F(Form("fConvTrueRecPtTaggedCaloMCIsoCell_%i_%i",r,e), Form("conversion photons with Cell ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fConvTrueRecPtTaggedCaloMCIsoCell[r][e] = (TH1F*) convTrueRecPtTaggedCaloMCIsoCell->Clone(Form("fConvTrueRecPtTaggedCaloMCIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fConvFolderTrue->Add(fConvTrueRecPtTaggedCaloMCIsoCell[r][e]);
        }

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

  if(fUseHistograms){

    fCaloPt = new TH1F("fCaloPt","calo photons in EMC acc;p_{T} (GeV/c); counts",nPtBins,minPt,maxPt);
    fCaloPtBeforeAcc = new TH1F("fCaloPtBeforeAcc", "calo photons all acc;p_{T} (GeV/c); counts", nPtBins,minPt,maxPt);

    fCaloPtTaggedCalo = new TH1F("fCaloPtTaggedCalo", "calo photons that survived tagging;p_{T} (GeV/c); counts", 500, 0, 50.);
    fCaloPtTaggedAsDecayCalo = new TH1F("fCaloPtTaggedAsDecayCalo", "calo photons that survived tagging;p_{T} (GeV/c); counts", 500, 0, 50.);


    fCaloRho = new TH1F("fCaloRho", "charged event density;#rho; counts", 500, 0, 50.);
    fCaloRhoTimesArea = new TH1F("fCaloRhoTimesArea", "charged event density;#rho #times jet Area; counts", 500, 0, 50.);
    
    fCaloFolderRec->Add(fCaloPt);
    fCaloFolderRec->Add(fCaloPtBeforeAcc);
    fCaloFolderRec->Add(fCaloPtTaggedCalo);
    fCaloFolderRec->Add(fCaloPtTaggedAsDecayCalo);
    fCaloFolderRec->Add(fCaloRho);
    fCaloFolderRec->Add(fCaloRhoTimesArea);

    for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
    {
      TH2F *caloIsoRawCharged = new TH2F(Form("caloIsoRawCharged_%i",r), Form("charged track ISO in R < %1.1f;#sum p_{T} (GeV/c); calo. p_{T} (GeV/c)",fTrackIsolationR.at(r)), nPtBins,minPt,maxPt,100,0,50);
      fCaloIsoRawCharged[r] = (TH2F*) caloIsoRawCharged->Clone(Form("fCaloIsoRawCharged_R%1.1f",fTrackIsolationR.at(r)));
      fCaloFolderRec->Add(fCaloIsoRawCharged[r]);
    }
    for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
    {
      TH2F *caloIsoRawNeutral = new TH2F(Form("caloIsoRawNeutral_%i",r), Form("Neutral track ISO in R < %1.1f;#sum p_{T} (GeV/c); calo. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,100,0,50);
      fCaloIsoRawNeutral[r] = (TH2F*) caloIsoRawNeutral->Clone(Form("fCaloIsoRawNeutral_R%1.1f",fNeutralIsolationR.at(r)));
      fCaloFolderRec->Add(fCaloIsoRawNeutral[r]);

      TH2F *caloIsoRawFull = new TH2F(Form("caloIsoRawFull_%i",r), Form("Full track ISO in R < %1.1f;#sum p_{T} (GeV/c); calo. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,100,0,50);
      fCaloIsoRawFull[r] = (TH2F*) caloIsoRawFull->Clone(Form("fCaloIsoRawFull_R%1.1f",fNeutralIsolationR.at(r)));
      fCaloFolderRec->Add(fCaloIsoRawFull[r]);

      TH2F *caloIsoCell = new TH2F(Form("caloIsoCell_%i",r), Form("Full track ISO in R < %1.1f;#sum p_{T} (GeV/c); calo. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,100,0,50);
      fCaloIsoCell[r] = (TH2F*) caloIsoCell->Clone(Form("fCaloIsoCell_R%1.1f",fNeutralIsolationR.at(r)));
      fCaloFolderRec->Add(fCaloIsoCell[r]);
    }

    fCaloFolderTrue          = new TList();
    fCaloFolderTrue->SetName("caloPhotonsTrue");
    fCaloFolderTrue->SetOwner(kTRUE);
    if(fIsMC > 0){
      fOutputList->Add(fCaloFolderTrue);
      fCaloTruePt = new TH1F("fCaloTruePt", "validated calo photons in EMC acceptance;p_{T} (GeV/c); counts", 500, 0, 50.);
      fCaloTruePtPrimary = new TH1F("fCaloTruePtPrimary", "calo photon that has not a pi0 etc. as mother;p_{T} (GeV/c); counts", 500, 0, 50.);
      fCaloTruePtDecay = new TH1F("fCaloTruePtDecay", "calo photon from decay;p_{T} (GeV/c); counts", 500, 0, 50.);
      fCaloTruePtDecayFoundOtherInCluster = new TH1F("fCaloTruePtDecayFoundOtherInCluster", "calo photon from decay, where the other decay particle was found in EMC;p_{T} (GeV/c); counts", 500, 0, 50.);
      fCaloTruePtDecayOtherInAcc = new TH1F("fCaloTruePtDecayOtherInAcc", "calo photon from decay, where the other decay particle was found;p_{T} (GeV/c); counts", 500, 0, 50.);
      fCaloTruePtDecayOtherInAccAboveMinEnergy = new TH1F("fCaloTruePtDecayOtherInAccAboveMinEnergy", "calo photon from decay, where the other decay particle is in acc. and above 0.7 GeV;p_{T} (GeV/c); counts", 500, 0, 50.);
  
      // iso and tagging studies
      fCaloTruePtTaggedCalo = new TH1F("fCaloTruePtTaggedCalo", "calo photons that survived tagging;p_{T} (GeV/c); counts", 500, 0, 50.);
      fCaloTruePtTaggedAsDecayCalo = new TH1F("fCaloTruePtTaggedAsDecayCalo", "calo photons that survived tagging;p_{T} (GeV/c); counts", 500, 0, 50.);

      // True with rec pT
      fCaloTrueRecPt = new TH1F("fCaloTrueRecPt", "validated calo photons in EMC acceptance;p_{T} (GeV/c); counts", 500, 0, 50.);
      fCaloTrueRecPtPrimary = new TH1F("fCaloTrueRecPtPrimary", "calo photon that has not a pi0 etc. as mother;p_{T} (GeV/c); counts", 500, 0, 50.);
      fCaloTrueRecPtDecay = new TH1F("fCaloTrueRecPtDecay", "calo photon from decay;p_{T} (GeV/c); counts", 500, 0, 50.);
      fCaloTrueRecPtDecayFoundOtherInCluster = new TH1F("fCaloTrueRecPtDecayFoundOtherInCluster", "calo photon from decay, where the other decay particle was found in EMC;p_{T} (GeV/c); counts", 500, 0, 50.);
      fCaloTrueRecPtDecayOtherInAcc = new TH1F("fCaloTrueRecPtDecayOtherInAcc", "calo photon from decay, where the other decay particle was found;p_{T} (GeV/c); counts", 500, 0, 50.);
      fCaloTrueRecPtDecayOtherInAccAboveMinEnergy = new TH1F("fCaloTrueRecPtDecayOtherInAccAboveMinEnergy", "calo photon from decay, where the other decay particle is in acc. and above 0.7 GeV;p_{T} (GeV/c); counts", 500, 0, 50.);

      // iso and tagging studies
      fCaloTrueRecPtTaggedCalo = new TH1F("fCaloTrueRecPtTaggedCalo", "calo photons that survived tagging;p_{T} (GeV/c); counts", 500, 0, 50.);
      fCaloTrueRecPtTaggedAsDecayCalo = new TH1F("fCaloTrueRecPtTaggedAsDecayCalo", "calo photons that survived tagging;p_{T} (GeV/c); counts", 500, 0, 50.);

      // add to folders
      fCaloFolderTrue->Add(fCaloTruePt);
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
        TH2F *caloTrueIsoRawCharged = new TH2F(Form("fCaloTrueIsoRawCharged_%i",r), Form("Charged track ISO in R < %1.1f;#sum p_{T} (GeV/c); counts",fTrackIsolationR.at(r)), nPtBins,minPt,maxPt,100,0,50);
        fCaloTrueIsoRawCharged[r] = (TH2F*) caloTrueIsoRawCharged->Clone(Form("fCaloTrueIsoRawCharged_R%1.1f",fTrackIsolationR.at(r)));
        fCaloFolderTrue->Add(fCaloTrueIsoRawCharged[r]); 

        TH2F *caloTrueIsoRawCharged_FromDecay = new TH2F(Form("fCaloTrueIsoRawCharged_FromDecay_%i",r), Form("Charged track ISO in R < %1.1f;#sum p_{T} (GeV/c); counts",fTrackIsolationR.at(r)), nPtBins,minPt,maxPt,100,0,50);
        fCaloTrueIsoRawCharged_FromDecay[r] = (TH2F*) caloTrueIsoRawCharged_FromDecay->Clone(Form("fCaloTrueIsoRawCharged_FromDecay_R%1.1f",fTrackIsolationR.at(r)));
        fCaloFolderTrue->Add(fCaloTrueIsoRawCharged_FromDecay[r]); 

        TH2F *caloTrueIsoRawCharged_FromDirect = new TH2F(Form("fCaloTrueIsoRawCharged_FromDirect_%i",r), Form("Charged track ISO in R < %1.1f;#sum p_{T} (GeV/c); counts",fTrackIsolationR.at(r)), nPtBins,minPt,maxPt,100,0,50);
        fCaloTrueIsoRawCharged_FromDirect[r] = (TH2F*) caloTrueIsoRawCharged_FromDirect->Clone(Form("fCaloTrueIsoRawCharged_FromDirect_R%1.1f",fTrackIsolationR.at(r)));
        fCaloFolderTrue->Add(fCaloTrueIsoRawCharged_FromDirect[r]); 
      }
      
      for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
      {
        TH2F *caloTrueIsoRawNeutral = new TH2F(Form("fCaloTrueIsoRawNeutral_%i",r), Form("Neutral track ISO in R < %1.1f;#sum p_{T} (GeV/c); calo. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,100,0,50);
        fCaloTrueIsoRawNeutral[r] = (TH2F*) caloTrueIsoRawNeutral->Clone(Form("fCaloTrueIsoRawNeutral_R%1.1f",fNeutralIsolationR.at(r)));
        fCaloFolderTrue->Add(fCaloTrueIsoRawNeutral[r]);

        TH2F *caloTrueIsoRawFull = new TH2F(Form("fCaloTrueIsoRawFull_%i",r), Form("Full track ISO in R < %1.1f;#sum p_{T} (GeV/c); calo. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,100,0,50);
        fCaloTrueIsoRawFull[r] = (TH2F*) caloTrueIsoRawFull->Clone(Form("fCaloTrueIsoRawFull_R%1.1f",fNeutralIsolationR.at(r)));
        fCaloFolderTrue->Add(fCaloTrueIsoRawFull[r]);

        TH2F *caloTrueIsoCell = new TH2F(Form("fCaloTrueIsoCell_%i",r), Form("Full track ISO in R < %1.1f;#sum p_{T} (GeV/c); calo. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,100,0,50);
        fCaloTrueIsoCell[r] = (TH2F*) caloTrueIsoCell->Clone(Form("fCaloTrueIsoCell_R%1.1f",fNeutralIsolationR.at(r)));
        fCaloFolderTrue->Add(fCaloTrueIsoCell[r]);

        TH2F *caloTrueIsoRawNeutral_FromDecay = new TH2F(Form("fCaloTrueIsoRawNeutral_FromDecay_%i",r), Form("Neutral track ISO in R < %1.1f;#sum p_{T} (GeV/c); calo. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,100,0,50);
        fCaloTrueIsoRawNeutral_FromDecay[r] = (TH2F*) caloTrueIsoRawNeutral_FromDecay->Clone(Form("fCaloTrueIsoRawNeutral_FromDecay_R%1.1f",fNeutralIsolationR.at(r)));
        fCaloFolderTrue->Add(fCaloTrueIsoRawNeutral_FromDecay[r]);

        TH2F *caloTrueIsoRawFull_FromDecay = new TH2F(Form("fCaloTrueIsoRawFull_FromDecay_%i",r), Form("Full track ISO in R < %1.1f;#sum p_{T} (GeV/c); calo. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,100,0,50);
        fCaloTrueIsoRawFull_FromDecay[r] = (TH2F*) caloTrueIsoRawFull_FromDecay->Clone(Form("fCaloTrueIsoRawFull_FromDecay_R%1.1f",fNeutralIsolationR.at(r)));
        fCaloFolderTrue->Add(fCaloTrueIsoRawFull_FromDecay[r]);

        TH2F *caloTrueIsoCell_FromDecay = new TH2F(Form("fCaloTrueIsoCell_FromDecay_%i",r), Form("Full track ISO in R < %1.1f;#sum p_{T} (GeV/c); calo. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,100,0,50);
        fCaloTrueIsoCell_FromDecay[r] = (TH2F*) caloTrueIsoCell_FromDecay->Clone(Form("fCaloTrueIsoCel_FromDecay_R%1.1f",fNeutralIsolationR.at(r)));
        fCaloFolderTrue->Add(fCaloTrueIsoCell_FromDecay[r]);

        TH2F *caloTrueIsoRawNeutral_FromDirect = new TH2F(Form("fCaloTrueIsoRawNeutral_FromDirect_%i",r), Form("Neutral track ISO in R < %1.1f;#sum p_{T} (GeV/c); calo. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,100,0,50);
        fCaloTrueIsoRawNeutral_FromDirect[r] = (TH2F*) caloTrueIsoRawNeutral_FromDirect->Clone(Form("fCaloTrueIsoRawNeutral_FromDirect_R%1.1f",fNeutralIsolationR.at(r)));
        fCaloFolderTrue->Add(fCaloTrueIsoRawNeutral_FromDirect[r]);

        TH2F *caloTrueIsoRawFull_FromDirect = new TH2F(Form("fCaloTrueIsoRawFull_FromDirect_%i",r), Form("Full track ISO in R < %1.1f;#sum p_{T} (GeV/c); calo. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,100,0,50);
        fCaloTrueIsoRawFull_FromDirect[r] = (TH2F*) caloTrueIsoRawFull_FromDirect->Clone(Form("fCaloTrueIsoRawFull_FromDirect_R%1.1f",fNeutralIsolationR.at(r)));
        fCaloFolderTrue->Add(fCaloTrueIsoRawFull_FromDirect[r]);

        TH2F *caloTrueIsoCell_FromDirect = new TH2F(Form("fCaloTrueIsoCell_FromDirect_%i",r), Form("Full track ISO in R < %1.1f;#sum p_{T} (GeV/c); calo. p_{T} (GeV/c)",fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt,100,0,50);
        fCaloTrueIsoCell_FromDirect[r] = (TH2F*) caloTrueIsoCell_FromDirect->Clone(Form("fCaloTrueIsoCel_FromDirect_R%1.1f",fNeutralIsolationR.at(r)));
        fCaloFolderTrue->Add(fCaloTrueIsoCell_FromDirect[r]);
      }
    }
    for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
    {
      for (UInt_t e = 0; e < fTrackIsolationE.size(); e++)
      {
        TH1F *caloPtIsoCharged = new TH1F(Form("fCaloPtIsoCharged_%i_%i",r,e), Form("calo photons with charged track ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
        fCaloPtIsoCharged[r][e] = (TH1F*) caloPtIsoCharged->Clone(Form("fCaloPtIsoCharged_R%1.1f_E%1.1f",fTrackIsolationE.at(e),fTrackIsolationR.at(r)));
        fCaloFolderRec->Add(fCaloPtIsoCharged[r][e]);
        TH1F *caloPtTaggedCaloIsoCharged = new TH1F(Form("fCaloPtTaggedCaloIsoCharged_%i_%i",r,e), Form("calo photons with charged track ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
        fCaloPtTaggedCaloIsoCharged[r][e] = (TH1F*) caloPtTaggedCaloIsoCharged->Clone(Form("fCaloPtTaggedCaloIsoCharged_R%1.1f_E%1.1f",fTrackIsolationE.at(e),fTrackIsolationR.at(r)));
        fCaloFolderRec->Add(fCaloPtTaggedCaloIsoCharged[r][e]);

        if(fIsMC>0){
          // true pT
          TH1F *caloTruePtIsoCharged = new TH1F(Form("fCaloTruePtIsoCharged_%i_%i",r,e), Form("calo photons with charged track ISO < %1.1f GeV in R < %1.1f;gen. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtIsoCharged[r][e] = (TH1F*) caloTruePtIsoCharged->Clone(Form("fCaloTruePtIsoCharged_R%1.1f_E%1.1f",fTrackIsolationE.at(e),fTrackIsolationR.at(r)));
          fCaloFolderTrue->Add(fCaloTruePtIsoCharged[r][e]);
          TH1F *caloTruePtTaggedCaloIsoCharged = new TH1F(Form("fCaloTruePtTaggedCaloIsoCharged_%i_%i",r,e), Form("calo photons with charged track ISO < %1.1f GeV in R < %1.1f + tagging;gen. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtTaggedCaloIsoCharged[r][e] = (TH1F*) caloTruePtTaggedCaloIsoCharged->Clone(Form("fCaloTruePtTaggedCaloIsoCharged_R%1.1f_E%1.1f",fTrackIsolationE.at(e),fTrackIsolationR.at(r)));
          fCaloFolderTrue->Add(fCaloTruePtTaggedCaloIsoCharged[r][e]);

          TH1F *caloTruePtMCIsoCharged = new TH1F(Form("fCaloTruePtMCIsoCharged_%i_%i",r,e), Form("calo photons with charged track ISO < %1.1f GeV in R < %1.1f;gen. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtMCIsoCharged[r][e] = (TH1F*) caloTruePtMCIsoCharged->Clone(Form("fCaloTruePtMCIsoCharged_R%1.1f_E%1.1f",fTrackIsolationE.at(e),fTrackIsolationR.at(r)));
          fCaloFolderTrue->Add(fCaloTruePtMCIsoCharged[r][e]);
          TH1F *caloTruePtTaggedCaloMCIsoCharged = new TH1F(Form("fCaloTruePtTaggedCaloMCIsoCharged_%i_%i",r,e), Form("calo photons with charged track ISO < %1.1f GeV in R < %1.1f + tagging;gen. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtTaggedCaloMCIsoCharged[r][e] = (TH1F*) caloTruePtTaggedCaloMCIsoCharged->Clone(Form("fCaloTruePtTaggedCaloMCIsoCharged_R%1.1f_E%1.1f",fTrackIsolationE.at(e),fTrackIsolationR.at(r)));
          fCaloFolderTrue->Add(fCaloTruePtTaggedCaloMCIsoCharged[r][e]);

          // rec Pt
          TH1F *caloTrueRecPtIsoCharged = new TH1F(Form("fCaloTrueRecPtIsoCharged_%i_%i",r,e), Form("calo photons with charged track ISO < %1.1f GeV in R < %1.1f;rec. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtIsoCharged[r][e] = (TH1F*) caloTrueRecPtIsoCharged->Clone(Form("fCaloTrueRecPtIsoCharged_R%1.1f_E%1.1f",fTrackIsolationE.at(e),fTrackIsolationR.at(r)));
          fCaloFolderTrue->Add(fCaloTrueRecPtIsoCharged[r][e]);
          TH1F *caloTrueRecPtTaggedCaloIsoCharged = new TH1F(Form("fCaloTrueRecPtTaggedCaloIsoCharged_%i_%i",r,e), Form("calo photons with charged track ISO < %1.1f GeV in R < %1.1f + tagging;rec. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtTaggedCaloIsoCharged[r][e] = (TH1F*) caloTrueRecPtTaggedCaloIsoCharged->Clone(Form("fCaloTrueRecPtTaggedCaloIsoCharged_R%1.1f_E%1.1f",fTrackIsolationE.at(e),fTrackIsolationR.at(r)));
          fCaloFolderTrue->Add(fCaloTrueRecPtTaggedCaloIsoCharged[r][e]);

          TH1F *caloTrueRecPtMCIsoCharged = new TH1F(Form("fCaloTrueRecPtMCIsoCharged_%i_%i",r,e), Form("calo photons with charged track ISO < %1.1f GeV in R < %1.1f;rec. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtMCIsoCharged[r][e] = (TH1F*) caloTrueRecPtMCIsoCharged->Clone(Form("fCaloTrueRecPtMCIsoCharged_R%1.1f_E%1.1f",fTrackIsolationE.at(e),fTrackIsolationR.at(r)));
          fCaloFolderTrue->Add(fCaloTrueRecPtMCIsoCharged[r][e]);
          TH1F *caloTrueRecPtTaggedCaloMCIsoCharged = new TH1F(Form("fCaloTrueRecPtTaggedCaloMCIsoCharged_%i_%i",r,e), Form("calo photons with charged track ISO < %1.1f GeV in R < %1.1f + tagging;rec. p_{T} (GeV/c); counts",fTrackIsolationE.at(e),fTrackIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtTaggedCaloMCIsoCharged[r][e] = (TH1F*) caloTrueRecPtTaggedCaloMCIsoCharged->Clone(Form("fCaloTrueRecPtTaggedCaloMCIsoCharged_R%1.1f_E%1.1f",fTrackIsolationE.at(e),fTrackIsolationR.at(r)));
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
        fCaloFolderRec->Add(fCaloPtIsoNeutral[r][e]);
        TH1F *caloPtIsoFull = new TH1F(Form("fCaloPtIsoFull_%i_%i",r,e), Form("calo photons with Full ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
        fCaloPtIsoFull[r][e] = (TH1F*) caloPtIsoFull->Clone(Form("fCaloPtIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
        fCaloFolderRec->Add(fCaloPtIsoFull[r][e]);
        TH1F *caloPtIsoCell = new TH1F(Form("fCaloPtIsoCell_%i_%i",r,e), Form("calo photons with Cell ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
        fCaloPtIsoCell[r][e] = (TH1F*) caloPtIsoCell->Clone(Form("fCaloPtIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
        fCaloFolderRec->Add(fCaloPtIsoCell[r][e]);

        TH1F *caloPtTaggedCaloIsoNeutral = new TH1F(Form("fCaloPtTaggedCaloIsoNeutral_%i_%i",r,e), Form("calo photons with Neutral ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
        fCaloPtTaggedCaloIsoNeutral[r][e] = (TH1F*) caloPtTaggedCaloIsoNeutral->Clone(Form("fCaloPtTaggedCaloIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
        fCaloFolderRec->Add(fCaloPtTaggedCaloIsoNeutral[r][e]);
        TH1F *caloPtTaggedCaloIsoFull = new TH1F(Form("fCaloPtTaggedCaloIsoFull_%i_%i",r,e), Form("calo photons with Full ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
        fCaloPtTaggedCaloIsoFull[r][e] = (TH1F*) caloPtTaggedCaloIsoFull->Clone(Form("fCaloPtTaggedCaloIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
        fCaloFolderRec->Add(fCaloPtTaggedCaloIsoFull[r][e]);
        TH1F *caloPtTaggedCaloIsoCell = new TH1F(Form("fCaloPtTaggedCaloIsoCell_%i_%i",r,e), Form("calo photons with Cell ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
        fCaloPtTaggedCaloIsoCell[r][e] = (TH1F*) caloPtTaggedCaloIsoCell->Clone(Form("fCaloPtTaggedCaloIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
        fCaloFolderRec->Add(fCaloPtTaggedCaloIsoCell[r][e]);

        // True Pt
        if(fIsMC>0){
          TH1F *caloTruePtIsoNeutral = new TH1F(Form("fCaloTruePtIsoNeutral_%i_%i",r,e), Form("calo photons with Neutral ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtIsoNeutral[r][e] = (TH1F*) caloTruePtIsoNeutral->Clone(Form("fCaloTruePtIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloFolderTrue->Add(fCaloTruePtIsoNeutral[r][e]);
          TH1F *caloTruePtIsoFull = new TH1F(Form("fCaloTruePtIsoFull_%i_%i",r,e), Form("calo photons with Full ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtIsoFull[r][e] = (TH1F*) caloTruePtIsoFull->Clone(Form("fCaloTruePtIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloFolderTrue->Add(fCaloTruePtIsoFull[r][e]);
          TH1F *caloTruePtIsoCell = new TH1F(Form("fCaloTruePtIsoCell_%i_%i",r,e), Form("calo photons with Cell ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtIsoCell[r][e] = (TH1F*) caloTruePtIsoCell->Clone(Form("fCaloTruePtIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloFolderTrue->Add(fCaloTruePtIsoCell[r][e]);

          TH1F *caloTruePtTaggedCaloIsoNeutral = new TH1F(Form("fCaloTruePtTaggedCaloIsoNeutral_%i_%i",r,e), Form("calo photons with Neutral ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtTaggedCaloIsoNeutral[r][e] = (TH1F*) caloTruePtTaggedCaloIsoNeutral->Clone(Form("fCaloTruePtTaggedCaloIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloFolderTrue->Add(fCaloTruePtTaggedCaloIsoNeutral[r][e]);
          TH1F *caloTruePtTaggedCaloIsoFull = new TH1F(Form("fCaloTruePtTaggedCaloIsoFull_%i_%i",r,e), Form("calo photons with Full ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtTaggedCaloIsoFull[r][e] = (TH1F*) caloTruePtTaggedCaloIsoFull->Clone(Form("fCaloTruePtTaggedCaloIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloFolderTrue->Add(fCaloTruePtTaggedCaloIsoFull[r][e]);
          TH1F *caloTruePtTaggedCaloIsoCell = new TH1F(Form("fCaloTruePtTaggedCaloIsoCell_%i_%i",r,e), Form("calo photons with Cell ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtTaggedCaloIsoCell[r][e] = (TH1F*) caloTruePtTaggedCaloIsoCell->Clone(Form("fCaloTruePtTaggedCaloIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloFolderTrue->Add(fCaloTruePtTaggedCaloIsoCell[r][e]);

          TH1F *caloTruePtMCIsoNeutral = new TH1F(Form("fCaloTruePtMCIsoNeutral_%i_%i",r,e), Form("calo photons with Neutral ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtMCIsoNeutral[r][e] = (TH1F*) caloTruePtMCIsoNeutral->Clone(Form("fCaloTruePtMCIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloFolderTrue->Add(fCaloTruePtMCIsoNeutral[r][e]);
          TH1F *caloTruePtMCIsoFull = new TH1F(Form("fCaloTruePtMCIsoFull_%i_%i",r,e), Form("calo photons with Full ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtMCIsoFull[r][e] = (TH1F*) caloTruePtMCIsoFull->Clone(Form("fCaloTruePtMCIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloFolderTrue->Add(fCaloTruePtMCIsoFull[r][e]);
          TH1F *caloTruePtMCIsoCell = new TH1F(Form("fCaloTruePtMCIsoCell_%i_%i",r,e), Form("calo photons with Cell ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtMCIsoCell[r][e] = (TH1F*) caloTruePtMCIsoCell->Clone(Form("fCaloTruePtMCIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloFolderTrue->Add(fCaloTruePtMCIsoCell[r][e]);

          TH1F *caloTruePtTaggedCaloMCIsoNeutral = new TH1F(Form("fCaloTruePtTaggedCaloMCIsoNeutral_%i_%i",r,e), Form("calo photons with Neutral ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtTaggedCaloMCIsoNeutral[r][e] = (TH1F*) caloTruePtTaggedCaloMCIsoNeutral->Clone(Form("fCaloTruePtTaggedCaloMCIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloFolderTrue->Add(fCaloTruePtTaggedCaloMCIsoNeutral[r][e]);
          TH1F *caloTruePtTaggedCaloMCIsoFull = new TH1F(Form("fCaloTruePtTaggedCaloMCIsoFull_%i_%i",r,e), Form("calo photons with Full ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtTaggedCaloMCIsoFull[r][e] = (TH1F*) caloTruePtTaggedCaloMCIsoFull->Clone(Form("fCaloTruePtTaggedCaloMCIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloFolderTrue->Add(fCaloTruePtTaggedCaloMCIsoFull[r][e]);
          TH1F *caloTruePtTaggedCaloMCIsoCell = new TH1F(Form("fCaloTruePtTaggedCaloMCIsoCell_%i_%i",r,e), Form("calo photons with Cell ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTruePtTaggedCaloMCIsoCell[r][e] = (TH1F*) caloTruePtTaggedCaloMCIsoCell->Clone(Form("fCaloTruePtTaggedCaloMCIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloFolderTrue->Add(fCaloTruePtTaggedCaloMCIsoCell[r][e]);

          // rec pT

          TH1F *caloTrueRecPtIsoNeutral = new TH1F(Form("fCaloTrueRecPtIsoNeutral_%i_%i",r,e), Form("calo photons with Neutral ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtIsoNeutral[r][e] = (TH1F*) caloTrueRecPtIsoNeutral->Clone(Form("fCaloTrueRecPtIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloFolderTrue->Add(fCaloTrueRecPtIsoNeutral[r][e]);
          TH1F *caloTrueRecPtIsoFull = new TH1F(Form("fCaloTrueRecPtIsoFull_%i_%i",r,e), Form("calo photons with Full ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtIsoFull[r][e] = (TH1F*) caloTrueRecPtIsoFull->Clone(Form("fCaloTrueRecPtIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloFolderTrue->Add(fCaloTrueRecPtIsoFull[r][e]);
          TH1F *caloTrueRecPtIsoCell = new TH1F(Form("fCaloTrueRecPtIsoCell_%i_%i",r,e), Form("calo photons with Cell ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtIsoCell[r][e] = (TH1F*) caloTrueRecPtIsoCell->Clone(Form("fCaloTrueRecPtIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloFolderTrue->Add(fCaloTrueRecPtIsoCell[r][e]);

          TH1F *caloTrueRecPtTaggedCaloIsoNeutral = new TH1F(Form("fCaloTrueRecPtTaggedCaloIsoNeutral_%i_%i",r,e), Form("calo photons with Neutral ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtTaggedCaloIsoNeutral[r][e] = (TH1F*) caloTrueRecPtTaggedCaloIsoNeutral->Clone(Form("fCaloTrueRecPtTaggedCaloIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloFolderTrue->Add(fCaloTrueRecPtTaggedCaloIsoNeutral[r][e]);
          TH1F *caloTrueRecPtTaggedCaloIsoFull = new TH1F(Form("fCaloTrueRecPtTaggedCaloIsoFull_%i_%i",r,e), Form("calo photons with Full ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtTaggedCaloIsoFull[r][e] = (TH1F*) caloTrueRecPtTaggedCaloIsoFull->Clone(Form("fCaloTrueRecPtTaggedCaloIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloFolderTrue->Add(fCaloTrueRecPtTaggedCaloIsoFull[r][e]);
          TH1F *caloTrueRecPtTaggedCaloIsoCell = new TH1F(Form("fCaloTrueRecPtTaggedCaloIsoCell_%i_%i",r,e), Form("calo photons with Cell ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtTaggedCaloIsoCell[r][e] = (TH1F*) caloTrueRecPtTaggedCaloIsoCell->Clone(Form("fCaloTrueRecPtTaggedCaloIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloFolderTrue->Add(fCaloTrueRecPtTaggedCaloIsoCell[r][e]);

          TH1F *caloTrueRecPtMCIsoNeutral = new TH1F(Form("fCaloTrueRecPtMCIsoNeutral_%i_%i",r,e), Form("calo photons with Neutral ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtMCIsoNeutral[r][e] = (TH1F*) caloTrueRecPtMCIsoNeutral->Clone(Form("fCaloTrueRecPtMCIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloFolderTrue->Add(fCaloTrueRecPtMCIsoNeutral[r][e]);
          TH1F *caloTrueRecPtMCIsoFull = new TH1F(Form("fCaloTrueRecPtMCIsoFull_%i_%i",r,e), Form("calo photons with Full ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtMCIsoFull[r][e] = (TH1F*) caloTrueRecPtMCIsoFull->Clone(Form("fCaloTrueRecPtMCIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloFolderTrue->Add(fCaloTrueRecPtMCIsoFull[r][e]);
          TH1F *caloTrueRecPtMCIsoCell = new TH1F(Form("fCaloTrueRecPtMCIsoCell_%i_%i",r,e), Form("calo photons with Cell ISO < %1.1f GeV in R < %1.1f;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtMCIsoCell[r][e] = (TH1F*) caloTrueRecPtMCIsoCell->Clone(Form("fCaloTrueRecPtMCIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloFolderTrue->Add(fCaloTrueRecPtMCIsoCell[r][e]);

          TH1F *caloTrueRecPtTaggedCaloMCIsoNeutral = new TH1F(Form("fCaloTrueRecPtTaggedCaloMCIsoNeutral_%i_%i",r,e), Form("calo photons with Neutral ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtTaggedCaloMCIsoNeutral[r][e] = (TH1F*) caloTrueRecPtTaggedCaloMCIsoNeutral->Clone(Form("fCaloTrueRecPtTaggedCaloMCIsoNeutral_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloFolderTrue->Add(fCaloTrueRecPtTaggedCaloMCIsoNeutral[r][e]);
          TH1F *caloTrueRecPtTaggedCaloMCIsoFull = new TH1F(Form("fCaloTrueRecPtTaggedCaloMCIsoFull_%i_%i",r,e), Form("calo photons with Full ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtTaggedCaloMCIsoFull[r][e] = (TH1F*) caloTrueRecPtTaggedCaloMCIsoFull->Clone(Form("fCaloTrueRecPtTaggedCaloMCIsoFull_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloFolderTrue->Add(fCaloTrueRecPtTaggedCaloMCIsoFull[r][e]);
          TH1F *caloTrueRecPtTaggedCaloMCIsoCell = new TH1F(Form("fCaloTrueRecPtTaggedCaloMCIsoCell_%i_%i",r,e), Form("calo photons with Cell ISO < %1.1f GeV in R < %1.1f + tagging;p_{T} (GeV/c); counts",fNeutralIsolationE.at(e),fNeutralIsolationR.at(r)), nPtBins,minPt,maxPt);
          fCaloTrueRecPtTaggedCaloMCIsoCell[r][e] = (TH1F*) caloTrueRecPtTaggedCaloMCIsoCell->Clone(Form("fCaloTrueRecPtTaggedCaloMCIsoCell_R%1.1f_E%1.1f",fNeutralIsolationR.at(r),fNeutralIsolationE.at(e)));
          fCaloFolderTrue->Add(fCaloTrueRecPtTaggedCaloMCIsoCell[r][e]);
        }

      }
    }
  }   
  
  PostData(1, fOutputList);
  
  OpenFile(2);
  fConversionCandidates = new TClonesArray("AliAODConversionPhoton",50);  
  fClusterEMCalCandidates = new TClonesArray("AliAODCaloCluster",50);   
  fClusterEMCalCandidatesIsolation = new TClonesArray("AliAODCaloCluster",50);  
  fClusterEMCalCandidatesTagging = new TClonesArray("AliAODCaloCluster",50);  
  fClusterPHOSCandidates = new TClonesArray("AliAODCaloCluster",50);  
  fTracks = new TClonesArray("AliAODTrack",10000);   
  fMCParticles = new TClonesArray("AliAODMCParticle",50000);  

  fExtraClusterInfo = new TClonesArray("AliExtraClusterInfoHelper",50);
  fExtraClusterInfoBackground = new TClonesArray("AliExtraClusterInfoHelper",50);
  fConvIsoInfo = new TClonesArray("AliIsoInfoHelper",50);
  fCaloIsoInfo = new TClonesArray("AliIsoInfoHelper",50);
  
  fAnalysisTree = new TTree("AnalysisTree","AnalysisTree");
  if(!fUseHistograms){ 
    
    Int_t split = 1;
    fAnalysisTree->Branch("fDataEvtHeader",&fDataEvtHeader,"pVtxX/d:pVtxY/d:pVtxZ/d:runnumber/I:numberESDtracks/I:rho/D",32000);
    fAnalysisTree->Branch("fConversionCandidates",&fConversionCandidates,32000,split);
    fAnalysisTree->Branch("fClusterEMCalCandidates",&fClusterEMCalCandidates,32000,split);
    fAnalysisTree->Branch("fClusterPHOSCandidates",&fClusterPHOSCandidates,32000,split);
    fAnalysisTree->Branch("fTracks",&fTracks,32000,split);
    fAnalysisTree->Branch("fExtraClusterInfo",&fExtraClusterInfo,32000,split);
    
    fAnalysisTree->Branch("fConvIsoInfo",&fConvIsoInfo,32000,split);
    fAnalysisTree->Branch("fCaloIsoInfo",&fCaloIsoInfo,32000,split);
    if(fIsMC>0){
      fAnalysisTree->Branch("fMCParticles",&fMCParticles,32000,split);
      fAnalysisTree->Branch("fMCEvtHeader",&fMCEvtHeader,"pVtxX/d:pVtxY/d:pVtxZ/d:runnumber/I:numberESDtracks/I:weightJJ/F:rho/D:evtType/i",32000);
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
  Int_t eventNotAccepted              = fEventCuts->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon,kFALSE);
  if(eventNotAccepted) return; // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1



  AliRhoParameter* outrho= (AliRhoParameter*) InputEvent()->FindListObject(fRhoOutName.Data());
  if(!outrho) AliFatal("could not find rho container!");
  if (fIsMC > 1){
      fWeightJetJetMC       = 1;
      Float_t maxjetpt      = -1.;
      Float_t pthard = -1;
      Bool_t isMCJet        = ((AliConvEventCuts*)fEventCuts)->IsJetJetMCEventAccepted( fMCEvent, fWeightJetJetMC ,pthard, fInputEvent, maxjetpt);
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


  //
  // ─── MAIN PROCESSING ────────────────────────────────────────────────────────────
  //

  ProcessTracks(); // always run ProcessTracks before calo photons! (even if save tracks is false)
  ProcessCaloPhotons(); // track matching is done here as well
  if(fSaveConversions)
    ProcessConversionPhotons();
  ReduceTrackInfo(); // track matching is done, we can remove cov matrix etc now
  if(fIsMC>0) ProcessMCParticles();
  // vertex
  Double_t vertex[3] = {0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
 
  fDataEvtHeader.pVtxX = vertex[0];
  fDataEvtHeader.pVtxY = vertex[1];
  fDataEvtHeader.pVtxZ = vertex[2];
  fDataEvtHeader.rho = outrho->GetVal();
  fDataEvtHeader.runnumber = InputEvent()->GetRunNumber();
  fDataEvtHeader.numberESDtracks = InputEvent()->GetNumberOfESDTracks();
  if(fIsMC>0){
    if(fMCEvent){
      fMCEvtHeader.pVtxX = fMCEvent->GetPrimaryVertex()->GetX();
      fMCEvtHeader.pVtxY = fMCEvent->GetPrimaryVertex()->GetY();
      fMCEvtHeader.pVtxZ = fMCEvent->GetPrimaryVertex()->GetZ();
    
      fMCEvtHeader.runnumber = fMCEvent->GetRunNumber();
      fMCEvtHeader.numberESDtracks = fMCEvent->GetNumberOfESDTracks();
      fMCEvtHeader.weightJJ = fWeightJetJetMC;
      fMCEvtHeader.rho = outrho->GetVal();
    // cout << "Event type = "  << fMCEvent->GetEventType() << endl;
      fMCEvtHeader.evtType = fMCEvent->GetEventType();
    }else{
      AliInfo("Could not find fMCEvent Header Information! Not saving any header info ...");
    }
  }


  // do not fill these in tree if user wants it
  // processing was needed anyways because of track matching
  // and isolation
  if(!fSaveTracks) fTracks->Delete();
  if(!fSaveEMCClusters) fClusterEMCalCandidates->Delete();
  if(!fSavePHOSClusters) fClusterPHOSCandidates->Delete();

  if( fIsMC > 0 && fInputEvent->IsA()==AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled())){
    RelabelAODPhotonCandidates(kFALSE); // Back to ESDMC Label
    fV0Reader->RelabelAODs(kFALSE);
  }
  // fill output
  if(!fUseHistograms){
    fAnalysisTree->Fill();
    if(fAnalysisTree->GetEntriesFast()%1000==0){
      fAnalysisTree->OptimizeBaskets();
    }
    PostData(2, fAnalysisTree);
  }
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
  // fConvIsoInfo->Delete();
  // fCaloIsoInfo->Delete();
  

  // cout << "Entry = " << fAnalysisTree->GetEntriesFast() << endl;
  // cout << "Autosave = " << fAnalysisTree->GetAutoSave() << endl;
  // cout << "Autoflush = " << fAnalysisTree->GetAutoFlush() << endl;

  fConversionCandidates->Delete();
  fClusterEMCalCandidates->Delete();
  fClusterEMCalCandidatesIsolation->Delete();
  fClusterEMCalCandidatesTagging->Delete();
  fClusterPHOSCandidates->Delete();
  fTracks->Delete();
  fMCParticles->Delete();
  fExtraClusterInfo->Delete();
  fExtraClusterInfoBackground->Delete();
  fConvIsoInfo->Delete();
  fCaloIsoInfo->Delete();

  fDataEvtHeader.pVtxX = -9999;
  fDataEvtHeader.pVtxY = -9999;
  fDataEvtHeader.pVtxZ = -9999;
  fDataEvtHeader.runnumber = -1;
  fDataEvtHeader.numberESDtracks = -1;
  
  fMCEvtHeader.pVtxX = -9999;
  fMCEvtHeader.pVtxY = -9999;
  fMCEvtHeader.pVtxZ = -9999;
  fMCEvtHeader.runnumber = -1;
  fMCEvtHeader.numberESDtracks = -1;
  fMCEvtHeader.evtType = 0;
}
//________________________________________________________________________
void AliAnalysisTaskGammaIsoTree::ProcessConversionPhotons(){
   fReaderGammas    = fV0Reader->GetReconstructedGammas();
   Int_t pos = 0;
   for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){
    AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(i);
    if(!PhotonCandidate) continue;
    if(!((AliConversionPhotonCuts*)fConvCuts)->PhotonIsSelected(PhotonCandidate,fInputEvent)) continue;
    new((*fConversionCandidates)[pos]) AliAODConversionPhoton(*PhotonCandidate);

    
    Int_t      tmp_tag= 0;
    Double32_t tmp_isoNeutral[2] = {0,0};
    Double32_t tmp_isoCell[2] = {0,0};
    Double32_t tmp_isoCharged[2] = {0,0};
    vector<Double32_t> isoCharged;
    vector<Double32_t> isoNeutral;
    vector<Double32_t> isoCell;
    
    if(fDoTrackIsolation) isoCharged = ProcessChargedIsolation(PhotonCandidate);
    if(fDoNeutralIsolation) isoNeutral = ProcessNeutralIsolation(PhotonCandidate);
    if(fDoCellIsolation) isoCell =  ProcessCellIsolation(PhotonCandidate);
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
    
    new((*fConvIsoInfo)[pos])AliIsoInfoHelper(tmp_isoCharged,tmp_isoNeutral,tmp_isoCell,tmp_tag);
    
    if(fUseHistograms) FillConversionHistos(PhotonCandidate,isoCharged,isoNeutral,isoCell,tmp_tag);

    if((fIsMC>0) && fUseHistograms) ProcessMCConversionPhoton(PhotonCandidate,isoCharged,isoNeutral,isoCell,tmp_tag);
   
    pos++;
    
   }
}
//________________________________________________________________________
void AliAnalysisTaskGammaIsoTree::ProcessMCConversionPhoton(AliAODConversionPhoton* photon,vector<Double32_t> isoCharged,vector<Double32_t> isoNeutral,vector<Double32_t> isoCell,Int_t tmptag){
  if(!photon) return;
  Bool_t isTrueConv = IsTrueConversionPhoton(photon);
  if(!isTrueConv) return;
  Bool_t isDecayPhoton = IsDecayPhoton(photon);

  if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  
  Int_t pconvlabel = GetConvPhotonMCLabel(photon);
  AliAODMCParticle* convPhotonMC = (AliAODMCParticle *) fAODMCTrackArray->At(pconvlabel);

  vector<Double32_t> mcIso;
  vector<Double32_t> mcIsoCharged;
  vector<Double32_t> mcIsoNeutral;
  vector<Double32_t> mcIsoFull;
  mcIso = ProcessMCIsolation(pconvlabel);

  for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
  {
    mcIsoCharged.push_back(mcIso.at(r));
  }
  for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
  {
    mcIsoNeutral.push_back(mcIso.at(r+fTrackIsolationR.size()));
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


  } else{
    fConvTruePtPrimary->Fill(convPhotonMC->Pt(), fWeightJetJetMC);
    fConvTrueRecPtPrimary->Fill(photon->Pt(), fWeightJetJetMC); 
  }

  // charged
  for (UInt_t i = 0; i < isoCharged.size(); i++)
  {
    if(i<5){
        fConvTrueIsoRawCharged[i]->Fill(isoCharged.at(i),photon->Pt(), fWeightJetJetMC);
        if(isDecayPhoton){
          fConvTrueIsoRawCharged_FromDecay[i]->Fill(isoCharged.at(i),photon->Pt(), fWeightJetJetMC);
        } else{
          fConvTrueIsoRawCharged_FromDirect[i]->Fill(isoCharged.at(i),photon->Pt(), fWeightJetJetMC);
        }
    }
  }
  
  // neutral
  for (UInt_t i = 0; i < isoNeutral.size(); i++)
  {
    if(i<5){
        fConvTrueIsoRawNeutral[i]->Fill(isoNeutral.at(i),photon->Pt(), fWeightJetJetMC);
        if(isDecayPhoton){
          fConvTrueIsoRawNeutral_FromDecay[i]->Fill(isoNeutral.at(i),photon->Pt(), fWeightJetJetMC);
        } else{
          fConvTrueIsoRawNeutral_FromDirect[i]->Fill(isoNeutral.at(i),photon->Pt(), fWeightJetJetMC);
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
        } else{
          fConvTrueIsoCell_FromDirect[i]->Fill(isoCell.at(i),photon->Pt(), fWeightJetJetMC);
        }
    }
  }
  // full
  for (UInt_t i = 0; i < isoNeutral.size(); i++)
  {
    if(i<5){
      fConvTrueIsoRawFull[i]->Fill(isoNeutral.at(i) + isoCharged.at(i),photon->Pt(), fWeightJetJetMC);
      if(isDecayPhoton){
        fConvTrueIsoRawFull_FromDecay[i]->Fill(isoNeutral.at(i) + isoCharged.at(i),photon->Pt(), fWeightJetJetMC);
      } else{
        fConvTrueIsoRawFull_FromDirect[i]->Fill(isoNeutral.at(i) + isoCharged.at(i),photon->Pt(), fWeightJetJetMC);
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
void AliAnalysisTaskGammaIsoTree::ProcessMCCaloPhoton(AliAODCaloCluster* clus,vector<Double32_t> isoCharged,vector<Double32_t> isoNeutral,vector<Double32_t> isoCell,Int_t tmptag){
  if(!clus) return;
  Double_t vertex[3] = {0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
  TLorentzVector v4cluster;
  clus->GetMomentum(v4cluster,vertex);

  if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));

  // check if photon cluster
  Bool_t isTruePhoton = kFALSE;
  Int_t *mclabelsCluster = clus->GetLabels();
  Int_t clusmclabel = 0;
  if (clus->GetNLabels() > 0)
  {
      // for (Int_t k = 0; k < (Int_t)clusterE->GetNLabels(); k++)
      // {
          if (((AliAODMCParticle* )fAODMCTrackArray->At(mclabelsCluster[0]))->PdgCode() == 22) {
              isTruePhoton = kTRUE;
              clusmclabel = mclabelsCluster[0];
          }
      // }
  }
  if(!isTruePhoton) return;
  
  // check if decay photon
  AliAODMCParticle* mcphoton = (AliAODMCParticle *)fAODMCTrackArray->At(clusmclabel);

  vector<Double32_t> mcIso;
  vector<Double32_t> mcIsoCharged;
  vector<Double32_t> mcIsoNeutral;
  vector<Double32_t> mcIsoFull;
  mcIso = ProcessMCIsolation(clusmclabel);

  for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
  {
    mcIsoCharged.push_back(mcIso.at(r));
  }
  for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
  {
    mcIsoNeutral.push_back(mcIso.at(r+fTrackIsolationR.size()));
  }
  for (UInt_t r = 0; r < mcIsoCharged.size(); r++)
  {
    mcIsoFull.push_back(mcIsoCharged.at(r)+ mcIsoNeutral.at(r)); // they should always have same length
  }
  

  Bool_t isDecayPhoton = IsDecayPhoton(mcphoton);

  fCaloTrueRecPt->Fill(v4cluster.Pt(), fWeightJetJetMC);
  fCaloTruePt->Fill(mcphoton->Pt(), fWeightJetJetMC);

  if(tmptag<2){
      fCaloTrueRecPtTaggedCalo->Fill(v4cluster.Pt(), fWeightJetJetMC);
      fCaloTruePtTaggedCalo->Fill(mcphoton->Pt(), fWeightJetJetMC);
  } else{
      fCaloTrueRecPtTaggedAsDecayCalo->Fill(v4cluster.Pt(), fWeightJetJetMC);
      fCaloTruePtTaggedAsDecayCalo->Fill(mcphoton->Pt(), fWeightJetJetMC);
  }
  if(isDecayPhoton){
    fCaloTruePtDecay->Fill(mcphoton->Pt(), fWeightJetJetMC);
    fCaloTrueRecPtDecay->Fill(v4cluster.Pt(), fWeightJetJetMC);

    // checkout mother
    Int_t labelMother = mcphoton->GetMother();
    AliAODMCParticle *calophotonMother = (AliAODMCParticle *) fAODMCTrackArray->At(labelMother);
    Int_t nDaughters = calophotonMother->GetNDaughters();
    Int_t otherDaughterLabel = -1;
    for (Int_t d = 0; d < nDaughters; d++)
    {
        Int_t tmp = calophotonMother->GetDaughterLabel(d);
        if (tmp == clusmclabel)
            continue;
        otherDaughterLabel = tmp;
    }
    if (otherDaughterLabel != -1)
    {
        if (IsInEMCalAcceptance((AliAODMCParticle*)fAODMCTrackArray->At(otherDaughterLabel)))
        {
            fCaloTruePtDecayOtherInAcc->Fill(mcphoton->Pt(), fWeightJetJetMC);
            fCaloTrueRecPtDecayOtherInAcc->Fill(v4cluster.Pt(), fWeightJetJetMC);
            if (((AliAODMCParticle*)fAODMCTrackArray->At(otherDaughterLabel))->E() >= 0.7){
                fCaloTruePtDecayOtherInAccAboveMinEnergy->Fill(mcphoton->Pt(), fWeightJetJetMC);
                fCaloTrueRecPtDecayOtherInAccAboveMinEnergy->Fill(v4cluster.Pt(), fWeightJetJetMC);
            }
        }
        Int_t clusterLabel = CheckClustersForMCContribution(otherDaughterLabel, fClusterEMCalCandidatesTagging);
        if (clusterLabel != -1)
        {
            fCaloTruePtDecayFoundOtherInCluster->Fill(mcphoton->Pt(), fWeightJetJetMC);
            fCaloTrueRecPtDecayFoundOtherInCluster->Fill(v4cluster.Pt(), fWeightJetJetMC);
        }
    }


  } else{
    fCaloTruePtPrimary->Fill(mcphoton->Pt(), fWeightJetJetMC);
    fCaloTrueRecPtPrimary->Fill(v4cluster.Pt(), fWeightJetJetMC); 
  }

  // charged
  for (UInt_t i = 0; i < isoCharged.size(); i++)
  {
    if(i<5){
        fCaloTrueIsoRawCharged[i]->Fill(isoCharged.at(i),v4cluster.Pt(), fWeightJetJetMC);
        if(isDecayPhoton){
          fCaloTrueIsoRawCharged_FromDecay[i]->Fill(isoCharged.at(i),v4cluster.Pt(), fWeightJetJetMC);
        } else{
          fCaloTrueIsoRawCharged_FromDirect[i]->Fill(isoCharged.at(i),v4cluster.Pt(), fWeightJetJetMC);
        }
    }
  }
  
  // neutral
  for (UInt_t i = 0; i < isoNeutral.size(); i++)
  {
    if(i<5){
        fCaloTrueIsoRawNeutral[i]->Fill(isoNeutral.at(i),v4cluster.Pt(), fWeightJetJetMC);
        if(isDecayPhoton){
          fCaloTrueIsoRawNeutral_FromDecay[i]->Fill(isoNeutral.at(i),v4cluster.Pt(), fWeightJetJetMC);
        } else{
          fCaloTrueIsoRawNeutral_FromDirect[i]->Fill(isoNeutral.at(i),v4cluster.Pt(), fWeightJetJetMC);
        }
    }
  }
  // cell
  for (UInt_t i = 0; i < isoCell.size(); i++)
  {
    if(i<5){
        fCaloTrueIsoCell[i]->Fill(isoCell.at(i),v4cluster.Pt(), fWeightJetJetMC);
        if(isDecayPhoton){
          fCaloTrueIsoCell_FromDecay[i]->Fill(isoCell.at(i),v4cluster.Pt(), fWeightJetJetMC);
        } else{
          fCaloTrueIsoCell_FromDirect[i]->Fill(isoCell.at(i),v4cluster.Pt(), fWeightJetJetMC);
        }
    }
  }
  // full
  for (UInt_t i = 0; i < isoNeutral.size(); i++)
  {
    if(i<5){
      fCaloTrueIsoRawFull[i]->Fill(isoNeutral.at(i) + isoCharged.at(i),v4cluster.Pt(), fWeightJetJetMC);
      if(isDecayPhoton){
        fCaloTrueIsoRawFull_FromDecay[i]->Fill(isoNeutral.at(i) + isoCharged.at(i),v4cluster.Pt(), fWeightJetJetMC);
      } else{
        fCaloTrueIsoRawFull_FromDirect[i]->Fill(isoNeutral.at(i) + isoCharged.at(i),v4cluster.Pt(), fWeightJetJetMC);
      }
    }
  }


  // Pt Isolated and tagged

  for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
  {
    for (UInt_t e = 0; e < fTrackIsolationE.size(); e++)
    {
      if(isoCharged.at(r) < fTrackIsolationE.at(e)){
        fCaloTruePtIsoCharged[r][e]->Fill(mcphoton->Pt(),fWeightJetJetMC);
        if(tmptag<2 )fCaloTruePtTaggedCaloIsoCharged[r][e]->Fill(mcphoton->Pt(),fWeightJetJetMC); // not tagged by calo
        fCaloTrueRecPtIsoCharged[r][e]->Fill(v4cluster.Pt(),fWeightJetJetMC);
        if(tmptag<2 )fCaloTrueRecPtTaggedCaloIsoCharged[r][e]->Fill(v4cluster.Pt(),fWeightJetJetMC); // not tagged by calo
      }

      // iso on gen level
      if(mcIsoCharged.at(r) < fTrackIsolationE.at(e)){
        fCaloTruePtMCIsoCharged[r][e]->Fill(mcphoton->Pt(),fWeightJetJetMC);
        if(tmptag<2 )fCaloTruePtTaggedCaloMCIsoCharged[r][e]->Fill(mcphoton->Pt(),fWeightJetJetMC); // not tagged by calo
        fCaloTrueRecPtMCIsoCharged[r][e]->Fill(v4cluster.Pt(),fWeightJetJetMC);
        if(tmptag<2 )fCaloTrueRecPtTaggedCaloMCIsoCharged[r][e]->Fill(v4cluster.Pt(),fWeightJetJetMC); // not tagged by calo
      }
    }
    
  }
  
  for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
  {
    for (UInt_t e = 0; e < fNeutralIsolationE.size(); e++)
    {
      if(isoNeutral.at(r) < fNeutralIsolationE.at(e))fCaloTruePtIsoNeutral[r][e]->Fill(mcphoton->Pt(),fWeightJetJetMC);
      if(isoCell.at(r) < fNeutralIsolationE.at(e)) fCaloTruePtIsoCell[r][e]->Fill(mcphoton->Pt(),fWeightJetJetMC);
      if((isoNeutral.at(r) + isoCharged.at(r)) < fNeutralIsolationE.at(e)) fCaloTruePtIsoFull[r][e]->Fill(mcphoton->Pt(),fWeightJetJetMC);
      
      if(isoNeutral.at(r) < fNeutralIsolationE.at(e))fCaloTrueRecPtIsoNeutral[r][e]->Fill(v4cluster.Pt(),fWeightJetJetMC);
      if(isoCell.at(r) < fNeutralIsolationE.at(e)) fCaloTrueRecPtIsoCell[r][e]->Fill(v4cluster.Pt(),fWeightJetJetMC);
      if((isoNeutral.at(r) + isoCharged.at(r)) < fNeutralIsolationE.at(e)) fCaloTrueRecPtIsoFull[r][e]->Fill(v4cluster.Pt(),fWeightJetJetMC);
      
      if(tmptag<2){
          if(isoNeutral.at(r) < fNeutralIsolationE.at(e))fCaloTruePtTaggedCaloIsoNeutral[r][e]->Fill(mcphoton->Pt(),fWeightJetJetMC);
          if(isoCell.at(r) < fNeutralIsolationE.at(e)) fCaloTruePtTaggedCaloIsoCell[r][e]->Fill(mcphoton->Pt(),fWeightJetJetMC);
          if((isoNeutral.at(r) + isoCharged.at(r)) < fNeutralIsolationE.at(e)) fCaloTruePtTaggedCaloIsoFull[r][e]->Fill(mcphoton->Pt(),fWeightJetJetMC);

          if(isoNeutral.at(r) < fNeutralIsolationE.at(e))fCaloTrueRecPtTaggedCaloIsoNeutral[r][e]->Fill(v4cluster.Pt(),fWeightJetJetMC);
          if(isoCell.at(r) < fNeutralIsolationE.at(e)) fCaloTrueRecPtTaggedCaloIsoCell[r][e]->Fill(v4cluster.Pt(),fWeightJetJetMC);
          if((isoNeutral.at(r) + isoCharged.at(r)) < fNeutralIsolationE.at(e)) fCaloTrueRecPtTaggedCaloIsoFull[r][e]->Fill(v4cluster.Pt(),fWeightJetJetMC);
      }

      // iso on gen level
      if(mcIsoNeutral.at(r) < fNeutralIsolationE.at(e))fCaloTruePtMCIsoNeutral[r][e]->Fill(mcphoton->Pt(),fWeightJetJetMC);
      if(mcIsoFull.at(r) < fNeutralIsolationE.at(e)) fCaloTruePtMCIsoFull[r][e]->Fill(mcphoton->Pt(),fWeightJetJetMC);
      
      if(mcIsoNeutral.at(r) < fNeutralIsolationE.at(e))fCaloTrueRecPtMCIsoNeutral[r][e]->Fill(v4cluster.Pt(),fWeightJetJetMC);
      if(mcIsoFull.at(r) < fNeutralIsolationE.at(e)) fCaloTrueRecPtMCIsoFull[r][e]->Fill(v4cluster.Pt(),fWeightJetJetMC);
      
      if(tmptag<2){
          if(mcIsoNeutral.at(r) < fNeutralIsolationE.at(e))fCaloTruePtTaggedCaloMCIsoNeutral[r][e]->Fill(mcphoton->Pt(),fWeightJetJetMC);
          if(mcIsoFull.at(r) < fNeutralIsolationE.at(e)) fCaloTruePtTaggedCaloMCIsoFull[r][e]->Fill(mcphoton->Pt(),fWeightJetJetMC);

          if(mcIsoNeutral.at(r) < fNeutralIsolationE.at(e))fCaloTrueRecPtTaggedCaloMCIsoNeutral[r][e]->Fill(v4cluster.Pt(),fWeightJetJetMC);
          if(mcIsoFull.at(r) < fNeutralIsolationE.at(e)) fCaloTrueRecPtTaggedCaloMCIsoFull[r][e]->Fill(v4cluster.Pt(),fWeightJetJetMC);
      }
    }
    
  }
}
//________________________________________________________________________
void AliAnalysisTaskGammaIsoTree::ProcessCaloPhotons(){
   Int_t nclus                         = 0;
   Int_t nclusCorr                     = 0;

   Int_t posEMC = 0;
   Int_t posEMCIso = 0;
   Int_t posEMCTag = 0;
   Int_t posPHOS = 0;
   TClonesArray * arrClustersProcess   = NULL;
   if(!fCorrTaskSetting.CompareTo("")){
     nclus = fInputEvent->GetNumberOfCaloClusters();

     nclusCorr = nclus;
   } else {
    arrClustersProcess                = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));
    if(!arrClustersProcess)
      AliFatal(Form("%sClustersBranch was not found in AliAnalysisTaskGammaCalo! Check the correction framework settings!",fCorrTaskSetting.Data()));
    nclusCorr                            = arrClustersProcess->GetEntries();
    nclus = fInputEvent->GetNumberOfCaloClusters();
  }
  if(nclus == 0)  return;
  // ((AliCaloPhotonCuts*)fClusterCutsEMC)->FillHistogramsExtendedQA(fInputEvent,fIsMC);
  // ((AliCaloPhotonCuts*)fClusterCutsPHOS)->FillHistogramsExtendedQA(fInputEvent,fIsMC);
  
  // in case user wants to use default track matching
  AliAODCaloCluster* clus                       = NULL;   
  if(arrClustersProcess){ 
     // EMCal correction framework was used
     // we need to loop over this for EMCal clusters and 
     // over the others for PHOS cluster

      // Loop over EMCal clusters
      for(Long_t i = 0; i < nclusCorr; i++){
        Double_t tempClusterWeight              =  fWeightJetJetMC;
        clus                                = new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersProcess->At(i));
        
        if(!clus) continue;
        if ( !clus->IsEMCAL()){ // for PHOS: cluster->GetType() == AliVCluster::kPHOSNeutral
          delete clus;
          continue;
        }

        if(!fDoOwnTrackMatching){
           ((AliCaloPhotonCuts*)fClusterCutsEMC)->MatchTracksToClusters(fInputEvent,fWeightJetJetMC,kTRUE, fMCEvent);
           ((AliCaloPhotonCuts*)fClusterCutsIsolationEMC)->MatchTracksToClusters(fInputEvent,fWeightJetJetMC,kTRUE, fMCEvent);
           //((AliCaloPhotonCuts*)fClusterCutsTaggingEMC)->MatchTracksToClusters(fInputEvent,fWeightJetJetMC,kTRUE, fMCEvent);
        }

        // get additional cluster info
        Short_t nLM = ((AliCaloPhotonCuts*)fClusterCutsEMC)->GetNumberOfLocalMaxima(clus, fInputEvent);
        Short_t matchIndex = -1;
        if(fDoOwnTrackMatching){
            matchIndex = ProcessTrackMatching(clus,fTracks);
        }
        Float_t eFrac = GetExoticEnergyFraction(clus,fInputEvent);
        if(((AliCaloPhotonCuts*)fClusterCutsIsolationEMC)->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC, tempClusterWeight,i)){
          if(!IsMatchedWithConv(clus,fClusterCutsIsolationEMC)){
            new((*fExtraClusterInfoBackground)[posEMCIso]) AliExtraClusterInfoHelper(nLM,matchIndex,eFrac);
            new((*fClusterEMCalCandidatesIsolation)[posEMCIso]) AliAODCaloCluster(*clus);
            posEMCIso++;
          }
        }

        if(((AliCaloPhotonCuts*)fClusterCutsTaggingEMC)->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC, tempClusterWeight,i)){
          if(!IsMatchedWithConv(clus,fClusterCutsTaggingEMC)){
            new((*fClusterEMCalCandidatesTagging)[posEMCTag]) AliAODCaloCluster(*clus);
            posEMCTag++;
          }
        }

        // check if given EMC cuts are fulfilled
        if(!((AliCaloPhotonCuts*)fClusterCutsEMC)->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC, tempClusterWeight,i)){
          delete clus;
          continue;
        }
        if(IsMatchedWithConv(clus,fClusterCutsEMC)){
          delete clus;
          continue;
        }
        new((*fExtraClusterInfo)[posEMC]) AliExtraClusterInfoHelper(nLM,matchIndex,eFrac);
        new((*fClusterEMCalCandidates)[posEMC])AliAODCaloCluster(*clus);
        
        posEMC++;
        delete clus;
      } // end of initial cluster loop
  }

  


  // Loop over normal clusters
  for(Long_t i = 0; i < nclus; i++){
    Double_t tempClusterWeight              =  fWeightJetJetMC;                   
    clus                                = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(i));
    
    if(!clus) continue;

    if(!arrClustersProcess && clus->IsEMCAL()){ // if is was not saved already
      // get additional cluster info
      if(!fDoOwnTrackMatching){
        ((AliCaloPhotonCuts*)fClusterCutsEMC)->MatchTracksToClusters(fInputEvent,fWeightJetJetMC,kTRUE, fMCEvent);
        ((AliCaloPhotonCuts*)fClusterCutsIsolationEMC)->MatchTracksToClusters(fInputEvent,fWeightJetJetMC,kTRUE, fMCEvent);
      }
      Short_t nLM = ((AliCaloPhotonCuts*)fClusterCutsEMC)->GetNumberOfLocalMaxima(clus, fInputEvent);
      Short_t matchIndex = -1;
      if(fDoOwnTrackMatching){
            matchIndex = ProcessTrackMatching(clus,fTracks);
      }
      Float_t eFrac = GetExoticEnergyFraction(clus,fInputEvent);

      if(((AliCaloPhotonCuts*)fClusterCutsIsolationEMC)->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC, tempClusterWeight,i)){
        if(!IsMatchedWithConv(clus,fClusterCutsIsolationEMC)){
          new((*fExtraClusterInfoBackground)[posEMCIso]) AliExtraClusterInfoHelper(nLM,matchIndex,eFrac);
          new((*fClusterEMCalCandidatesIsolation)[posEMCIso]) AliAODCaloCluster(*clus);
          posEMCIso++;
        }
      }

      if(((AliCaloPhotonCuts*)fClusterCutsTaggingEMC)->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC, tempClusterWeight,i)){
        if(!IsMatchedWithConv(clus,fClusterCutsTaggingEMC)){
          new((*fClusterEMCalCandidatesTagging)[posEMCTag]) AliAODCaloCluster(*clus);
          posEMCTag++;
        }
      }
      if(!((AliCaloPhotonCuts*)fClusterCutsEMC)->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC, tempClusterWeight,i)){
        delete clus;
        continue;
      }
      if(IsMatchedWithConv(clus,fClusterCutsEMC)){
        delete clus;
        continue;
      }

      new((*fExtraClusterInfo)[posEMC]) AliExtraClusterInfoHelper(nLM,matchIndex,eFrac);
      new((*fClusterEMCalCandidates)[posEMC]) AliAODCaloCluster(*clus);
      posEMC++;
      
      delete clus;
      continue;
    }
    if(clus->IsPHOS() && fSavePHOSClusters){
    // if(clus->GetType() == AliVCluster::kPHOSNeutral){
      if(!((AliCaloPhotonCuts*)fClusterCutsPHOS)->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC, tempClusterWeight,i)){
        delete clus;
        continue;
      }
      if(!fSavePHOSClusters) new((*fClusterPHOSCandidates)[posPHOS]) AliAODCaloCluster(*clus);
      posPHOS++;
      delete clus;
      continue;
    }
    delete clus;
  }

  // Loop again over selected cluster candidates to do isolation and tagging
  for (Int_t c = 0; c < fClusterEMCalCandidates->GetEntriesFast(); c++)
  {
     AliAODCaloCluster* clus = (AliAODCaloCluster*) fClusterEMCalCandidates->At(c);
     Double32_t tmp_isoCharged[2] = {0,0};
     Double32_t tmp_isoNeutral[2] = {0,0};
     Double32_t tmp_isoCell[2] = {0,0};
     Int_t      tmp_tag= 0;
     vector<Double32_t> isoCharged;
     vector<Double32_t> isoNeutral;
     vector<Double32_t> isoCell;
     if(fDoTrackIsolation) isoCharged = ProcessChargedIsolation(clus);
     if(fDoNeutralIsolation) isoNeutral =ProcessNeutralIsolation(clus);
     if(fDoCellIsolation) isoCell = ProcessCellIsolation(clus);
     if(fDoTagging) tmp_tag = ProcessTagging(clus);  

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
     new((*fCaloIsoInfo)[c])AliIsoInfoHelper(tmp_isoCharged,tmp_isoNeutral,tmp_isoCell,tmp_tag);

     if(fUseHistograms) FillCaloHistos(clus,isoCharged,isoNeutral,isoCell,tmp_tag);

     if((fIsMC>0) && fUseHistograms) ProcessMCCaloPhoton(clus,isoCharged,isoNeutral,isoCell,tmp_tag);
  }
  
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
      AliAODTrack *fCurrentTrack = dynamic_cast<AliAODTrack*> (fInputEvent->GetTrack(t));
      //if(!TrackIsSelectedAOD(fCurrentTrack)){
        // save empty track to preserve position
      //  new((*fTracks)[pos]) AliAODTrack();
      //} else{

      // we need to save all tracks in order to identify tracks from conv for isolation 
      new((*fTracks)[pos]) AliAODTrack(*fCurrentTrack);
      pos++;
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskGammaIsoTree::ProcessMCParticles(){
  // Loop over all primary MC particle
  if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  Int_t pos = 0;
  if (fAODMCTrackArray){
    for(Int_t i = 0; i < fAODMCTrackArray->GetEntriesFast(); i++) {
      AliAODMCParticle* particle =  static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(i));
      if(TMath::Abs(particle->Y())< fYMCCut){
        new((*fMCParticles)[pos]) AliAODMCParticle(*particle);
        pos++;

      } else
      {
        new((*fMCParticles)[pos]) AliAODMCParticle();
        pos++;
      }
    }
  }
}

//_____________________________________________________________________________
Int_t AliAnalysisTaskGammaIsoTree::ProcessTrackMatching(AliAODCaloCluster* clus, TClonesArray* tracks){
     Int_t nModules = fGeomEMCAL->GetNumberOfSuperModules();
     Int_t highestMatchIndex = -1;
     AliExternalTrackParam *trackParam = 0;
     for (Int_t t = 0; t < tracks->GetEntriesFast(); t++)
     {
      AliAODTrack *aodt = dynamic_cast<AliAODTrack*>(tracks->At(t));
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
          AliAODTrack *convtrack = dynamic_cast<AliAODTrack*>(tracks->At(tracklabel));
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
        AliAODTrack *aodt = dynamic_cast<AliAODTrack*>(fInputEvent->GetTrack(t));
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
          AliAODTrack *convtrack = dynamic_cast<AliAODTrack*>(fInputEvent->GetTrack(tracklabel));
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
vector<Double32_t> AliAnalysisTaskGammaIsoTree::ProcessChargedIsolation(AliAODCaloCluster* cluster){
    Double_t vertex[3] = {0,0,0};
    InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
    
    TLorentzVector v4cluster;
    cluster->GetMomentum(v4cluster,vertex);
    vector<Double32_t> vecIso; 
    for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
    {
        vecIso.push_back(0.);
    }
    for (Int_t t = 0; t < fInputEvent->GetNumberOfTracks(); t++)
    {
        AliAODTrack *aodt = dynamic_cast<AliAODTrack*>(fInputEvent->GetTrack(t));
        if(!aodt) continue;
        if(!TrackIsSelectedAOD(aodt)) continue;
        TLorentzVector v4track;
        v4track.SetPxPyPzE(aodt->Px(),aodt->Py(),aodt->Pz(),aodt->E());        
        Double_t trackEta = v4track.Eta();
        Double_t trackPhi = v4track.Phi();
        if (trackPhi < 0) trackPhi += 2*TMath::Pi();

        Double_t clusterEta = v4cluster.Eta();
        Double_t clusterPhi = v4cluster.Phi();
        if (clusterPhi < 0) clusterPhi += 2*TMath::Pi();

        Double_t dEta = trackEta - clusterEta;
        Double_t dPhi = trackPhi - clusterPhi;
        Double_t dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);
        for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
        {
           if(dR <= fTrackIsolationR.at(r)) vecIso.at(r) += v4track.Pt();
        }
    }
    fHistoChargedIso->Fill(vecIso.at(0)); // debug only
    return vecIso;
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
    for (Int_t c = 0; c < fClusterEMCalCandidatesIsolation->GetEntriesFast(); c++)
    {
        AliAODCaloCluster* clusterE = (AliAODCaloCluster*) fClusterEMCalCandidatesIsolation->At(c);
        if(!clusterE) continue;

        AliExtraClusterInfoHelper* clusInfo = (AliExtraClusterInfoHelper*) fExtraClusterInfoBackground->At(c);
        // check if cluster is neutral
        if(fDoBackgroundTrackMatching){
           if(clusInfo->isMatched()) continue;
        }
     
    
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
    //     AliAODTrack *convt = dynamic_cast<AliAODTrack*>(fInputEvent->GetTrack(tracklabel));
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
    for (Int_t c = 0; c < fClusterEMCalCandidatesIsolation->GetEntriesFast(); c++)
    {
        AliAODCaloCluster* clusterE = (AliAODCaloCluster*) fClusterEMCalCandidatesIsolation->At(c);
        if(!clusterE) continue;
        if(clusterE->GetID() == cluster->GetID()) continue;
        AliExtraClusterInfoHelper* clusInfo = (AliExtraClusterInfoHelper*) fExtraClusterInfoBackground->At(c);
        if(fDoBackgroundTrackMatching){
           if(clusInfo->isMatched()) continue;
        }
 
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
vector<Double32_t> AliAnalysisTaskGammaIsoTree::ProcessMCIsolation(Int_t mclabel){
  vector<Double32_t> vecIso;
  if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));

  for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
  {
      vecIso.push_back(0.);
  }
  for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
  {
      vecIso.push_back(0.);
  }
  AliAODMCParticle *thisParticle = (AliAODMCParticle*) fAODMCTrackArray->At(mclabel);


  Double_t thisEta = thisParticle->Eta();
  Double_t thisPhi = thisParticle->Phi();

  for (Int_t p = 0; p < fAODMCTrackArray->GetEntriesFast(); p++)
  {
      if (p == mclabel)
          continue;
      AliAODMCParticle *pmc = (AliAODMCParticle *)fAODMCTrackArray->At(p);
      if(pmc->GetPdgCode()==0) continue;
      if (!pmc)
          continue;
      if (pmc->MCStatusCode() != 1)
          continue;
      if(!pmc->IsPhysicalPrimary()) continue;
      if(pmc->IsSecondaryFromMaterial()) continue;

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
        for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
        {
           if(dR<fTrackIsolationR.at(r)) vecIso.at(r)+= pmc->Pt(); 
        }
      }
      if ((pmc->Charge() == 0) && (IsInEMCalAcceptance(pmc))){ //neutral

        for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
        {
           if(dR<fNeutralIsolationR.at(r)) vecIso.at(r+fTrackIsolationR.size()) += pmc->Pt(); 
        }
      }
  }
  return vecIso;
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

  for (Int_t c = 0; c < fClusterEMCalCandidatesTagging->GetEntriesFast(); c++)
  {
    AliExtraClusterInfoHelper* clusInfo = (AliExtraClusterInfoHelper*) fExtraClusterInfoBackground->At(c);
    if(fDoBackgroundTrackMatching){
       if(clusInfo->isMatched()) continue;
    }
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


  for (Int_t c = 0; c < fClusterEMCalCandidatesTagging->GetEntriesFast(); c++)
  {
    // TLorentzvector with cluster
    TLorentzVector clusterVector;
    
    if(((AliAODCaloCluster*)fClusterEMCalCandidatesTagging->At(c))->GetID() == cluster->GetID()) continue;
    
    AliExtraClusterInfoHelper* clusInfo = (AliExtraClusterInfoHelper*) fExtraClusterInfoBackground->At(c);
    if(fDoBackgroundTrackMatching){
       if(clusInfo->isMatched()) continue;
    }
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
    for (Int_t i = 0; i < fTracks->GetEntriesFast(); i++)
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
      fConvIsoRawCharged[r]->Fill(isoCharged.at(r),photon->Pt(),fWeightJetJetMC);
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
      fConvIsoRawNeutral[r]->Fill(isoNeutral.at(r),photon->Pt(),fWeightJetJetMC);
      fConvIsoRawFull[r]->Fill(isoNeutral.at(r) + isoCharged.at(r),photon->Pt(),fWeightJetJetMC);
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
    

}
void AliAnalysisTaskGammaIsoTree::FillCaloHistos(AliAODCaloCluster* clus,vector<Double32_t> isoCharged,vector<Double32_t> isoNeutral,vector<Double32_t> isoCell,Int_t tmptag){
    if(!clus) return;
    Double_t vertex[3] = {0};
    InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
    TLorentzVector v4cluster;
    clus->GetMomentum(v4cluster,vertex);
    
    fCaloPtBeforeAcc->Fill(v4cluster.Pt(),fWeightJetJetMC);
    //if (!IsInEMCalAcceptance(photon)) return;

    fCaloPt->Fill(v4cluster.Pt(),fWeightJetJetMC);
    if(tmptag<2 ){
       fCaloPtTaggedCalo->Fill(v4cluster.Pt(),fWeightJetJetMC);
    } else{
       fCaloPtTaggedAsDecayCalo->Fill(v4cluster.Pt(),fWeightJetJetMC);
    }
    
    for (UInt_t r = 0; r < fTrackIsolationR.size(); r++)
    {
      fCaloIsoRawCharged[r]->Fill(isoCharged.at(r),v4cluster.Pt(),fWeightJetJetMC);
      for (UInt_t e = 0; e < fTrackIsolationE.size(); e++)
      {
        if(isoCharged.at(r) < fTrackIsolationE.at(e)){
          fCaloPtIsoCharged[r][e]->Fill(v4cluster.Pt(),fWeightJetJetMC);
          if(tmptag<2 )fCaloPtTaggedCaloIsoCharged[r][e]->Fill(v4cluster.Pt(),fWeightJetJetMC); // not tagged by calo
        }
      }
      
    }
    
    for (UInt_t r = 0; r < fNeutralIsolationR.size(); r++)
    {
      fCaloIsoRawNeutral[r]->Fill(isoNeutral.at(r),v4cluster.Pt(),fWeightJetJetMC);
      fCaloIsoRawFull[r]->Fill(isoNeutral.at(r) + isoCharged.at(r),v4cluster.Pt(),fWeightJetJetMC);
      fCaloIsoCell[r]->Fill(isoCell.at(r),v4cluster.Pt(),fWeightJetJetMC);
      for (UInt_t e = 0; e < fNeutralIsolationE.size(); e++)
      {
        if(isoNeutral.at(r) < fNeutralIsolationE.at(e))fCaloPtIsoNeutral[r][e]->Fill(v4cluster.Pt(),fWeightJetJetMC);
        if(isoCell.at(r) < fNeutralIsolationE.at(e)) fCaloPtIsoCell[r][e]->Fill(v4cluster.Pt(),fWeightJetJetMC);
        if((isoNeutral.at(r) + isoCharged.at(r)) < fNeutralIsolationE.at(e)) fCaloPtIsoFull[r][e]->Fill(v4cluster.Pt(),fWeightJetJetMC);
        if(tmptag<2){
            if(isoNeutral.at(r) < fNeutralIsolationE.at(e))fCaloPtTaggedCaloIsoNeutral[r][e]->Fill(v4cluster.Pt(),fWeightJetJetMC);
            if(isoCell.at(r) < fNeutralIsolationE.at(e)) fCaloPtTaggedCaloIsoCell[r][e]->Fill(v4cluster.Pt(),fWeightJetJetMC);
            if((isoNeutral.at(r) + isoCharged.at(r)) < fNeutralIsolationE.at(e)) fCaloPtTaggedCaloIsoFull[r][e]->Fill(v4cluster.Pt(),fWeightJetJetMC);
        }
      }
      
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
    // cout << phi << endl;
    // cout << eta << endl;
    if (phi < 0)
        phi += 2 * TMath::Pi();
    if ((eta < -0.6687) || (eta > 0.66465))
        return kFALSE;
    if ((phi < 1.39626) || (phi > 3.15))
        return kFALSE;
    return kTRUE;
}
Bool_t AliAnalysisTaskGammaIsoTree::IsInEMCalAcceptance(AliAODMCParticle *part)
{
    Double_t eta = part->Eta();
    Double_t phi = part->Phi();

    if (phi < 0)
        phi += 2 * TMath::Pi();
    // cout << phi << endl;
    // cout << eta << endl;
    if ((eta < -0.6687) || (eta > 0.66465))
        return kFALSE;
    if ((phi < 1.39626) || (phi > 3.15))
        return kFALSE;
    return kTRUE;
}
Bool_t AliAnalysisTaskGammaIsoTree::IsTrueConversionPhoton(AliAODConversionPhoton* photon)
{
    if ((photon->GetMCLabelPositive() == -1) || (photon->GetMCLabelNegative() == -1))
        return false;
    if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    
    AliAODMCParticle *posDaughter = (AliAODMCParticle *)fAODMCTrackArray->At(photon->GetMCLabelPositive());
    AliAODMCParticle *negDaughter = (AliAODMCParticle *)fAODMCTrackArray->At(photon->GetMCLabelNegative());

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
    AliAODMCParticle *mother = (AliAODMCParticle *)fAODMCTrackArray->At(posDaughter->GetMother());
    if (mother->GetPdgCode() != 22)
        return kFALSE;

    return kTRUE;
}
Int_t AliAnalysisTaskGammaIsoTree::GetConvPhotonMCLabel(AliAODConversionPhoton *photon)
{
    if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    AliAODMCParticle *posDaughter = (AliAODMCParticle *)fAODMCTrackArray->At(photon->GetMCLabelPositive());
    AliAODMCParticle *negDaughter = (AliAODMCParticle *)fAODMCTrackArray->At(photon->GetMCLabelNegative());

    if (posDaughter->GetMother() != negDaughter->GetMother())
        return -1;

    return posDaughter->GetMother();
}
Bool_t AliAnalysisTaskGammaIsoTree::IsDecayPhoton(AliAODMCParticle *mcphoton){
    if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    AliAODMCParticle *mcphotonmother = (AliAODMCParticle *)fAODMCTrackArray->At(mcphoton->GetMother());
    Int_t pdgMom = mcphotonmother->GetPdgCode();
    if (TMath::Abs(pdgMom) > 100 && TMath::Abs(pdgMom) < 1000)
    {
        return kTRUE;
    }
    else
    {
        return kFALSE;
    }

}
Bool_t AliAnalysisTaskGammaIsoTree::IsDecayPhoton(AliAODConversionPhoton *photon){
    Bool_t isFromDecay = kFALSE;
    if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!IsTrueConversionPhoton(photon))
        return kFALSE;
    Int_t motherlabel = GetConvPhotonMCLabel(photon);
    AliAODMCParticle *mcphoton = (AliAODMCParticle *)fAODMCTrackArray->At(motherlabel);
    if (mcphoton->GetPdgCode() != 22)
        return kFALSE;

    isFromDecay = IsDecayPhoton(mcphoton);

    return isFromDecay;
}

Int_t AliAnalysisTaskGammaIsoTree::CheckClustersForMCContribution(Int_t mclabel, TClonesArray *vclus)
{
    Int_t clusterLabel = -1; // position of cluster in array where mc label was found as contribution
    for (Int_t p = 0; p < vclus->GetEntriesFast(); p++)
    {
        AliAODCaloCluster *clus = (AliAODCaloCluster *)vclus->At(p);
        if (!clus)
            continue;
        Int_t *mclabelsCluster = clus->GetLabels();
        if (clus->GetNLabels() > 0)
        {
            for (Int_t k = 0; k < (Int_t)clus->GetNLabels(); k++)
            {
                if (mclabelsCluster[p] == mclabel)
                    clusterLabel = p;
            }
        }
    }

    return clusterLabel;
}
Int_t AliAnalysisTaskGammaIsoTree::CheckConvForMCContribution(Int_t mclabel, TClonesArray *vconv)
{
  Int_t convLabel = -1; // position of conversion in array where mc label was found as contribution
   for (Int_t p = 0; p < vconv->GetEntriesFast(); p++)
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


