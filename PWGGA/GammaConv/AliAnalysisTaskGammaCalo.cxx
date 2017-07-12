/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                       *
 * Author: Baldo Sahlmueller, Friederike Bock                     *
 * Version 1.0                                 *
 *                                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its    *
 * documentation strictly for non-commercial purposes is hereby granted    *
 * without fee, provided that the above copyright notice appears in all    *
 * copies and that both the copyright notice and this permission notice    *
 * appear in the supporting documentation. The authors make no claims    *
 * about the suitability of this software for any purpose. It is      *
 * provided "as is" without express or implied warranty.               *
 **************************************************************************/

//////////////////////////////////////////////////////////////////
//----------------------------------------------------------------
// Class used to do analysis on conversion photons + calo photons
//----------------------------------------------------------------
//////////////////////////////////////////////////////////////////
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliCentrality.h"
#include "AliESDVZERO.h"
#include "AliESDpid.h"
#include "AliAnalysisTaskGammaCalo.h"
#include "AliVParticle.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliKFVertex.h"
#include "AliV0ReaderV1.h"
#include "AliGenCocktailEventHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliEventplane.h"
#include "AliAnalysisTaskEMCALClusterizeFast.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliInputEventHandler.h"
#include "AliCaloTrackMatcher.h"
#include <vector>
#include <map>
#include <fstream>

ClassImp(AliAnalysisTaskGammaCalo)

//________________________________________________________________________
AliAnalysisTaskGammaCalo::AliAnalysisTaskGammaCalo(): AliAnalysisTaskSE(),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fBGHandler(NULL),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fCutFolder(NULL),
  fESDList(NULL),
  fBackList(NULL),
  fMotherList(NULL),
  fTrueList(NULL),
  fMCList(NULL),
  fTreeList(NULL),
  fClusterTreeList(NULL),
  fOutputContainer(NULL),
  fClusterCandidates(NULL),
  fEventCutArray(NULL),
  fEventCuts(NULL),
  fClusterCutArray(NULL),
  fCaloPhotonCuts(NULL),
  fMesonCutArray(NULL),
  fMesonCuts(NULL),
  fHistoMotherInvMassPt(NULL),
  fSparseMotherInvMassPtZM(NULL),
  fHistoMotherBackInvMassPt(NULL),
  fSparseMotherBackInvMassPtZM(NULL),
  fHistoMotherInvMassPtAlpha(NULL),
  fHistoMotherBackInvMassPtAlpha(NULL),
  fHistoMotherPi0PtY(NULL),
  fHistoMotherEtaPtY(NULL),
  fHistoMotherPi0PtAlpha(NULL),
  fHistoMotherEtaPtAlpha(NULL),
  fHistoMotherPi0PtOpenAngle(NULL),
  fHistoMotherEtaPtOpenAngle(NULL),
  fHistoClusGammaPt(NULL),
  fHistoClusGammaE(NULL),
  fHistoClusOverlapHeadersGammaPt(NULL),
  fHistoClusGammaPtM02(NULL),
  fHistoMCHeaders(NULL),
  fHistoMCAllGammaPt(NULL),
  fHistoMCAllSecondaryGammaPt(NULL),
  fHistoMCDecayGammaPi0Pt(NULL),
  fHistoMCDecayGammaRhoPt(NULL),
  fHistoMCDecayGammaEtaPt(NULL),
  fHistoMCDecayGammaOmegaPt(NULL),
  fHistoMCDecayGammaEtapPt(NULL),
  fHistoMCDecayGammaPhiPt(NULL),
  fHistoMCDecayGammaSigmaPt(NULL),
  fHistoMCPi0Pt(NULL),
  fHistoMCPi0WOWeightPt(NULL),
  fHistoMCPi0WOEvtWeightPt(NULL),
  fHistoMCEtaPt(NULL),
  fHistoMCEtaWOWeightPt(NULL),
  fHistoMCEtaWOEvtWeightPt(NULL),
  fHistoMCPi0InAccPt(NULL),
  fHistoMCEtaInAccPt(NULL),
  fHistoMCPi0WOEvtWeightInAccPt(NULL),
  fHistoMCEtaWOEvtWeightInAccPt(NULL),
  fHistoMCPi0PtY(NULL),
  fHistoMCEtaPtY(NULL),
  fHistoMCPi0PtAlpha(NULL),
  fHistoMCEtaPtAlpha(NULL),
  fHistoMCPrimaryPtvsSource(NULL),
  fHistoMCSecPi0PtvsSource(NULL),
  fHistoMCSecPi0Source(NULL),
  fHistoMCSecPi0InAccPtvsSource(NULL),
  fHistoMCSecEtaPt(NULL),
  fHistoMCSecEtaSource(NULL),
  fHistoMCPi0PtJetPt(NULL),
  fHistoMCEtaPtJetPt(NULL),
  fHistoTruePi0InvMassPt(NULL),
  fHistoTrueEtaInvMassPt(NULL),
  fHistoTruePi0CaloPhotonInvMassPt(NULL),
  fHistoTrueEtaCaloPhotonInvMassPt(NULL),
  fHistoTruePi0CaloConvertedPhotonInvMassPt(NULL),
  fHistoTrueEtaCaloConvertedPhotonInvMassPt(NULL),
  fHistoTruePi0CaloMixedPhotonConvPhotonInvMassPt(NULL),
  fHistoTrueEtaCaloMixedPhotonConvPhotonInvMassPt(NULL),
  fHistoTruePi0CaloElectronInvMassPt(NULL),
  fHistoTrueEtaCaloElectronInvMassPt(NULL),
  fHistoTruePi0CaloMergedClusterInvMassPt(NULL),
  fHistoTrueEtaCaloMergedClusterInvMassPt(NULL),
  fHistoTruePi0CaloMergedClusterPartConvInvMassPt(NULL),
  fHistoTrueEtaCaloMergedClusterPartConvInvMassPt(NULL),
  fHistoTruePi0NonMergedElectronPhotonInvMassPt(NULL),
  fHistoTruePi0NonMergedElectronMergedPhotonInvMassPt(NULL),
  fHistoTruePi0Category1(NULL),
  fHistoTrueEtaCategory1(NULL),
  fHistoTruePi0Category2(NULL),
  fHistoTrueEtaCategory2(NULL),
  fHistoTruePi0Category3(NULL),
  fHistoTrueEtaCategory3(NULL),
  fHistoTruePi0Category4_6(NULL),
  fHistoTrueEtaCategory4_6(NULL),
  fHistoTruePi0Category5(NULL),
  fHistoTrueEtaCategory5(NULL),
  fHistoTruePi0Category7(NULL),
  fHistoTrueEtaCategory7(NULL),
  fHistoTruePi0Category8(NULL),
  fHistoTrueEtaCategory8(NULL),
  fHistoTruePrimaryPi0InvMassPt(NULL),
  fHistoTruePrimaryEtaInvMassPt(NULL),
  fHistoTruePrimaryPi0W0WeightingInvMassPt(NULL),
  fHistoTruePrimaryEtaW0WeightingInvMassPt(NULL),
  fProfileTruePrimaryPi0WeightsInvMassPt(NULL),
  fProfileTruePrimaryEtaWeightsInvMassPt(NULL),
  fHistoTruePrimaryPi0MCPtResolPt(NULL),
  fHistoTruePrimaryEtaMCPtResolPt(NULL),
  fHistoTrueSecondaryPi0InvMassPt(NULL),
  fHistoTrueSecondaryPi0FromK0sInvMassPt(NULL),
  fHistoTrueSecondaryPi0FromK0lInvMassPt(NULL),
  fHistoTrueK0sWithPi0DaughterMCPt(NULL),
  fHistoTrueK0lWithPi0DaughterMCPt(NULL),
  fHistoTrueSecondaryPi0FromEtaInvMassPt(NULL),
  fHistoTrueEtaWithPi0DaughterMCPt(NULL),
  fHistoTrueSecondaryPi0FromLambdaInvMassPt(NULL),
  fHistoTrueLambdaWithPi0DaughterMCPt(NULL),
  fHistoTrueBckGGInvMassPt(NULL),
  fHistoTrueBckFullMesonContainedInOneClusterInvMassPt(NULL),
  fHistoTrueBckAsymEClustersInvMassPt(NULL),
  fHistoTrueBckContInvMassPt(NULL),
  fHistoTruePi0PtY(NULL),
  fHistoTrueEtaPtY(NULL),
  fHistoTruePi0PtAlpha(NULL),
  fHistoTrueEtaPtAlpha(NULL),
  fHistoTruePi0PtOpenAngle(NULL),
  fHistoTrueEtaPtOpenAngle(NULL),
  fHistoClusPhotonBGPt(NULL),
  fHistoClusPhotonPlusConvBGPt(NULL),
  fHistoClustPhotonElectronBGPtM02(NULL),
  fHistoClustPhotonPionBGPtM02(NULL),
  fHistoClustPhotonKaonBGPtM02(NULL),
  fHistoClustPhotonK0lBGPtM02(NULL),
  fHistoClustPhotonNeutronBGPtM02(NULL),
  fHistoClustPhotonRestBGPtM02(NULL),
  fHistoClustPhotonPlusConvElectronBGPtM02(NULL),
  fHistoClustPhotonPlusConvPionBGPtM02(NULL),
  fHistoClustPhotonPlusConvKaonBGPtM02(NULL),
  fHistoClustPhotonPlusConvK0lBGPtM02(NULL),
  fHistoClustPhotonPlusConvNeutronBGPtM02(NULL),
  fHistoClustPhotonPlusConvRestBGPtM02(NULL),
  fHistoTrueClusGammaPt(NULL),
  fHistoTrueClusUnConvGammaPt(NULL),
  fHistoTrueClusUnConvGammaMCPt(NULL),
  fHistoTrueClusGammaPtM02(NULL),
  fHistoTrueClusUnConvGammaPtM02(NULL),
  fHistoTrueClusElectronPt(NULL),
  fHistoTrueClusConvGammaPt(NULL),
  fHistoTrueClusConvGammaMCPt(NULL),
  fHistoTrueClusConvGammaFullyPt(NULL),
  fHistoTrueClusMergedGammaPt(NULL),
  fHistoTrueClusMergedPartConvGammaPt(NULL),
  fHistoTrueClusDalitzPt(NULL),
  fHistoTrueClusDalitzMergedPt(NULL),
  fHistoTrueClusPhotonFromElecMotherPt(NULL),
  fHistoTrueClusShowerPt(NULL),
  fHistoTrueClusSubLeadingPt(NULL),
  fHistoTrueClusNParticles(NULL),
  fHistoTrueClusEMNonLeadingPt(NULL),
  fHistoTrueNLabelsInClus(NULL),
  fHistoTruePrimaryClusGammaPt(NULL),
  fHistoTruePrimaryClusGammaESDPtMCPt(NULL),
  fHistoTruePrimaryClusConvGammaPt(NULL),
  fHistoTruePrimaryClusConvGammaESDPtMCPt(NULL),
  fHistoTrueSecondaryClusGammaPt(NULL),
  fHistoTrueSecondaryClusConvGammaPt(NULL),
  fHistoTrueSecondaryClusGammaMCPt(NULL),
  fHistoTrueSecondaryClusConvGammaMCPt(NULL),
  fHistoTrueSecondaryClusGammaFromXFromK0sMCPtESDPt(NULL),
  fHistoTrueSecondaryClusConvGammaFromXFromK0sMCPtESDPt(NULL),
  fHistoTrueSecondaryClusGammaFromXFromK0lMCPtESDPt(NULL),
  fHistoTrueSecondaryClusConvGammaFromXFromK0lMCPtESDPt(NULL),
  fHistoTrueSecondaryClusGammaFromXFromLambdaMCPtESDPt(NULL),
  fHistoTrueSecondaryClusConvGammaFromXFromLambdaMCPtESDPt(NULL),
  fHistoDoubleCountTruePi0InvMassPt(NULL),
  fHistoDoubleCountTrueEtaInvMassPt(NULL),
  fHistoDoubleCountTrueClusterGammaPt(NULL),
  fHistoMultipleCountTrueClusterGamma(NULL),
  fVectorDoubleCountTruePi0s(0),
  fVectorDoubleCountTrueEtas(0),
  fVectorDoubleCountTrueClusterGammas(0),
  fMapMultipleCountTrueClusterGammas(),
  fHistoTruePi0InvMassPtAlpha(NULL),
  fHistoTruePi0PureGammaInvMassPtAlpha(NULL),
  fHistCellIDvsClusterEnergy(NULL),
  fHistCellIDvsClusterEnergyMax(NULL),
  fHistoNEvents(NULL),
  fHistoNEventsWOWeight(NULL),
  fHistoNGoodESDTracks(NULL),
  fHistoVertexZ(NULL),
  fHistoNGammaCandidates(NULL),
  fHistoNGoodESDTracksVsNGammaCandidates(NULL),
  fHistoSPDClusterTrackletBackground(NULL),
  fHistoNV0Tracks(NULL),
  fProfileEtaShift(NULL),
  fProfileJetJetXSection(NULL),
  fHistoJetJetNTrials(NULL),
  tTrueInvMassROpenABPtFlag(NULL),
  tSigInvMassPtAlphaTheta(NULL),
  tBckInvMassPtAlphaTheta(NULL),
  fInvMassTreeInvMass(0),
  fInvMassTreePt(0),
  fInvMassTreeAlpha(0),
  fInvMassTreeTheta(0),
  fInvMassTreeMixPool(0),
  fInvMassTreeZVertex(0),
  fInvMassTreeEta(0),
  fInvMass(-1),
  fRconv(-1),
  fOpenRPrim(-1),
  fInvMassRTOF(-1),
  fPt(-1),
  iFlag(3),
  tClusterEOverP(NULL),
  fClusterE(0),
  fClusterM02(0),
  fClusterM20(0),
  fClusterEP(0),
  fClusterLeadCellID(0),
  fClusterClassification(0),
  fDeltaEta(0),
  fDeltaPhi(0),
  fTrackPt(0),
  fTrackPID_e(0),
  fTrackPID_Pi(0),
  fTrackPID_K(0),
  fTrackPID_P(0),
  fClusterIsoSumClusterEt(0),
  fClusterIsoSumTrackEt(0),
//  fHistoTruePi0NonLinearity(NULL),
//  fHistoTrueEtaNonLinearity(NULL),
  fEventPlaneAngle(-100),
  fRandom(0),
  fnCuts(0),
  fiCut(0),
  fIsHeavyIon(0),
  fDoLightOutput(kFALSE),
  fDoMesonAnalysis(kTRUE),
  fDoMesonQA(0),
  fDoClusterQA(0),
  fIsFromMBHeader(kTRUE),
  fIsOverlappingWithOtherHeader(kFALSE),
  fIsMC(0),
  fDoTHnSparse(kTRUE),
  fSetPlotHistsExtQA(kFALSE),
  fWeightJetJetMC(1),
  fDoInOutTimingCluster(kFALSE),
  fMinTimingCluster(0),
  fMaxTimingCluster(0),
  fEnableSortForClusMC(kFALSE),
  fProduceTreeEOverP(kFALSE),
  fProduceCellIDPlots(kFALSE),
  tBrokenFiles(NULL),
  fFileNameBroken(NULL),
  fCloseHighPtClusters(NULL),
  fLocalDebugFlag(0)
{
  
}

//________________________________________________________________________
AliAnalysisTaskGammaCalo::AliAnalysisTaskGammaCalo(const char *name):
  AliAnalysisTaskSE(name),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fBGHandler(NULL),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fCutFolder(NULL),
  fESDList(NULL),
  fBackList(NULL),
  fMotherList(NULL),
  fTrueList(NULL),
  fClusterTreeList(NULL),
  fMCList(NULL),
  fTreeList(NULL),
  fOutputContainer(0),
  fClusterCandidates(NULL),
  fEventCutArray(NULL),
  fEventCuts(NULL),
  fClusterCutArray(NULL),
  fCaloPhotonCuts(NULL),
  fMesonCutArray(NULL),
  fMesonCuts(NULL),
  fHistoMotherInvMassPt(NULL),
  fSparseMotherInvMassPtZM(NULL),
  fHistoMotherBackInvMassPt(NULL),
  fSparseMotherBackInvMassPtZM(NULL),
  fHistoMotherInvMassPtAlpha(NULL),
  fHistoMotherBackInvMassPtAlpha(NULL),
  fHistoMotherPi0PtY(NULL),
  fHistoMotherEtaPtY(NULL),
  fHistoMotherPi0PtAlpha(NULL),
  fHistoMotherEtaPtAlpha(NULL),
  fHistoMotherPi0PtOpenAngle(NULL),
  fHistoMotherEtaPtOpenAngle(NULL),
  fHistoClusGammaPt(NULL),
  fHistoClusGammaE(NULL),
  fHistoClusOverlapHeadersGammaPt(NULL),
  fHistoClusGammaPtM02(NULL),
  fHistoMCHeaders(NULL),
  fHistoMCAllGammaPt(NULL),
  fHistoMCAllSecondaryGammaPt(NULL),
  fHistoMCDecayGammaPi0Pt(NULL),
  fHistoMCDecayGammaRhoPt(NULL),
  fHistoMCDecayGammaEtaPt(NULL),
  fHistoMCDecayGammaOmegaPt(NULL),
  fHistoMCDecayGammaEtapPt(NULL),
  fHistoMCDecayGammaPhiPt(NULL),
  fHistoMCDecayGammaSigmaPt(NULL),
  fHistoMCPi0Pt(NULL),
  fHistoMCPi0WOWeightPt(NULL),
  fHistoMCPi0WOEvtWeightPt(NULL),
  fHistoMCEtaPt(NULL),
  fHistoMCEtaWOWeightPt(NULL),
  fHistoMCEtaWOEvtWeightPt(NULL),
  fHistoMCPi0InAccPt(NULL),
  fHistoMCEtaInAccPt(NULL),
  fHistoMCPi0WOEvtWeightInAccPt(NULL),
  fHistoMCEtaWOEvtWeightInAccPt(NULL),
  fHistoMCPi0PtY(NULL),
  fHistoMCEtaPtY(NULL),
  fHistoMCPi0PtAlpha(NULL),
  fHistoMCEtaPtAlpha(NULL),
  fHistoMCPrimaryPtvsSource(NULL),
  fHistoMCSecPi0PtvsSource(NULL),
  fHistoMCSecPi0Source(NULL),
  fHistoMCSecPi0InAccPtvsSource(NULL),
  fHistoMCSecEtaPt(NULL),
  fHistoMCSecEtaSource(NULL),
  fHistoMCPi0PtJetPt(NULL),
  fHistoMCEtaPtJetPt(NULL),  
  fHistoTruePi0InvMassPt(NULL),
  fHistoTrueEtaInvMassPt(NULL),
  fHistoTruePi0CaloPhotonInvMassPt(NULL),
  fHistoTrueEtaCaloPhotonInvMassPt(NULL),
  fHistoTruePi0CaloConvertedPhotonInvMassPt(NULL),
  fHistoTrueEtaCaloConvertedPhotonInvMassPt(NULL),
  fHistoTruePi0CaloMixedPhotonConvPhotonInvMassPt(NULL),
  fHistoTrueEtaCaloMixedPhotonConvPhotonInvMassPt(NULL),
  fHistoTruePi0CaloElectronInvMassPt(NULL),
  fHistoTrueEtaCaloElectronInvMassPt(NULL),
  fHistoTruePi0CaloMergedClusterInvMassPt(NULL),
  fHistoTrueEtaCaloMergedClusterInvMassPt(NULL),
  fHistoTruePi0CaloMergedClusterPartConvInvMassPt(NULL),
  fHistoTrueEtaCaloMergedClusterPartConvInvMassPt(NULL),
  fHistoTruePi0NonMergedElectronPhotonInvMassPt(NULL),
  fHistoTruePi0NonMergedElectronMergedPhotonInvMassPt(NULL),
  fHistoTruePi0Category1(NULL),
  fHistoTrueEtaCategory1(NULL),
  fHistoTruePi0Category2(NULL),
  fHistoTrueEtaCategory2(NULL),
  fHistoTruePi0Category3(NULL),
  fHistoTrueEtaCategory3(NULL),
  fHistoTruePi0Category4_6(NULL),
  fHistoTrueEtaCategory4_6(NULL),
  fHistoTruePi0Category5(NULL),
  fHistoTrueEtaCategory5(NULL),
  fHistoTruePi0Category7(NULL),
  fHistoTrueEtaCategory7(NULL),
  fHistoTruePi0Category8(NULL),
  fHistoTrueEtaCategory8(NULL),
  fHistoTruePrimaryPi0InvMassPt(NULL),
  fHistoTruePrimaryEtaInvMassPt(NULL),
  fHistoTruePrimaryPi0W0WeightingInvMassPt(NULL),
  fHistoTruePrimaryEtaW0WeightingInvMassPt(NULL),
  fProfileTruePrimaryPi0WeightsInvMassPt(NULL),
  fProfileTruePrimaryEtaWeightsInvMassPt(NULL),
  fHistoTruePrimaryPi0MCPtResolPt(NULL),
  fHistoTruePrimaryEtaMCPtResolPt(NULL),
  fHistoTrueSecondaryPi0InvMassPt(NULL),
  fHistoTrueSecondaryPi0FromK0sInvMassPt(NULL),
  fHistoTrueSecondaryPi0FromK0lInvMassPt(NULL),
  fHistoTrueK0sWithPi0DaughterMCPt(NULL),
  fHistoTrueK0lWithPi0DaughterMCPt(NULL),
  fHistoTrueSecondaryPi0FromEtaInvMassPt(NULL),
  fHistoTrueEtaWithPi0DaughterMCPt(NULL),
  fHistoTrueSecondaryPi0FromLambdaInvMassPt(NULL),
  fHistoTrueLambdaWithPi0DaughterMCPt(NULL),
  fHistoTrueBckGGInvMassPt(NULL),
  fHistoTrueBckFullMesonContainedInOneClusterInvMassPt(NULL),
  fHistoTrueBckAsymEClustersInvMassPt(NULL),
  fHistoTrueBckContInvMassPt(NULL),
  fHistoTruePi0PtY(NULL),
  fHistoTrueEtaPtY(NULL),
  fHistoTruePi0PtAlpha(NULL),
  fHistoTrueEtaPtAlpha(NULL),
  fHistoTruePi0PtOpenAngle(NULL),
  fHistoTrueEtaPtOpenAngle(NULL),
  fHistoClusPhotonBGPt(NULL),
  fHistoClusPhotonPlusConvBGPt(NULL),
  fHistoClustPhotonElectronBGPtM02(NULL),
  fHistoClustPhotonPionBGPtM02(NULL),
  fHistoClustPhotonKaonBGPtM02(NULL),
  fHistoClustPhotonK0lBGPtM02(NULL),
  fHistoClustPhotonNeutronBGPtM02(NULL),
  fHistoClustPhotonRestBGPtM02(NULL),
  fHistoClustPhotonPlusConvElectronBGPtM02(NULL),
  fHistoClustPhotonPlusConvPionBGPtM02(NULL),
  fHistoClustPhotonPlusConvKaonBGPtM02(NULL),
  fHistoClustPhotonPlusConvK0lBGPtM02(NULL),
  fHistoClustPhotonPlusConvNeutronBGPtM02(NULL),
  fHistoClustPhotonPlusConvRestBGPtM02(NULL),
  fHistoTrueClusGammaPt(NULL),
  fHistoTrueClusUnConvGammaPt(NULL),
  fHistoTrueClusUnConvGammaMCPt(NULL),
  fHistoTrueClusGammaPtM02(NULL),
  fHistoTrueClusUnConvGammaPtM02(NULL),
  fHistoTrueClusElectronPt(NULL),
  fHistoTrueClusConvGammaPt(NULL),
  fHistoTrueClusConvGammaMCPt(NULL),
  fHistoTrueClusConvGammaFullyPt(NULL),
  fHistoTrueClusMergedGammaPt(NULL),
  fHistoTrueClusMergedPartConvGammaPt(NULL),
  fHistoTrueClusDalitzPt(NULL),
  fHistoTrueClusDalitzMergedPt(NULL),
  fHistoTrueClusPhotonFromElecMotherPt(NULL),
  fHistoTrueClusShowerPt(NULL),
  fHistoTrueClusSubLeadingPt(NULL),
  fHistoTrueClusNParticles(NULL),
  fHistoTrueClusEMNonLeadingPt(NULL),
  fHistoTrueNLabelsInClus(NULL),
  fHistoTruePrimaryClusGammaPt(NULL),
  fHistoTruePrimaryClusGammaESDPtMCPt(NULL),
  fHistoTruePrimaryClusConvGammaPt(NULL),
  fHistoTruePrimaryClusConvGammaESDPtMCPt(NULL),
  fHistoTrueSecondaryClusGammaPt(NULL),
  fHistoTrueSecondaryClusConvGammaPt(NULL),
  fHistoTrueSecondaryClusGammaMCPt(NULL),
  fHistoTrueSecondaryClusConvGammaMCPt(NULL),
  fHistoTrueSecondaryClusGammaFromXFromK0sMCPtESDPt(NULL),
  fHistoTrueSecondaryClusConvGammaFromXFromK0sMCPtESDPt(NULL),
  fHistoTrueSecondaryClusGammaFromXFromK0lMCPtESDPt(NULL),
  fHistoTrueSecondaryClusConvGammaFromXFromK0lMCPtESDPt(NULL),
  fHistoTrueSecondaryClusGammaFromXFromLambdaMCPtESDPt(NULL),
  fHistoTrueSecondaryClusConvGammaFromXFromLambdaMCPtESDPt(NULL),
  fHistoDoubleCountTruePi0InvMassPt(NULL),
  fHistoDoubleCountTrueEtaInvMassPt(NULL),
  fHistoDoubleCountTrueClusterGammaPt(NULL),
  fHistoMultipleCountTrueClusterGamma(NULL),
  fVectorDoubleCountTruePi0s(0),
  fVectorDoubleCountTrueEtas(0),
  fVectorDoubleCountTrueClusterGammas(0),
  fMapMultipleCountTrueClusterGammas(),
  fHistoTruePi0InvMassPtAlpha(NULL),
  fHistoTruePi0PureGammaInvMassPtAlpha(NULL),
  fHistCellIDvsClusterEnergy(NULL),
  fHistCellIDvsClusterEnergyMax(NULL),
  fHistoNEvents(NULL),
  fHistoNEventsWOWeight(NULL),
  fHistoNGoodESDTracks(NULL),
  fHistoVertexZ(NULL),
  fHistoNGammaCandidates(NULL),
  fHistoNGoodESDTracksVsNGammaCandidates(NULL),
  fHistoSPDClusterTrackletBackground(NULL),
  fHistoNV0Tracks(NULL),
  fProfileEtaShift(NULL),
  fProfileJetJetXSection(NULL),
  fHistoJetJetNTrials(NULL),
  tTrueInvMassROpenABPtFlag(NULL),
  tSigInvMassPtAlphaTheta(NULL),
  tBckInvMassPtAlphaTheta(NULL),
  fInvMassTreeInvMass(0),
  fInvMassTreePt(0),
  fInvMassTreeAlpha(0),
  fInvMassTreeTheta(0),
  fInvMassTreeMixPool(0),
  fInvMassTreeZVertex(0),
  fInvMassTreeEta(0),
  fInvMass(-1),
  fRconv(-1),
  fOpenRPrim(-1),
  fInvMassRTOF(-1),
  fPt(-1),
  iFlag(3),
  tClusterEOverP(NULL),
  fClusterE(0),
  fClusterM02(0),
  fClusterM20(0),
  fClusterEP(0),
  fClusterLeadCellID(0),
  fClusterClassification(0),
  fDeltaEta(0),
  fDeltaPhi(0),
  fTrackPt(0),
  fTrackPID_e(0),
  fTrackPID_Pi(0),
  fTrackPID_K(0),
  fTrackPID_P(0),
  fClusterIsoSumClusterEt(0),
  fClusterIsoSumTrackEt(0),
//  fHistoTruePi0NonLinearity(NULL),
//  fHistoTrueEtaNonLinearity(NULL),
  fEventPlaneAngle(-100),
  fRandom(0),
  fnCuts(0),
  fiCut(0),
  fIsHeavyIon(0),
  fDoLightOutput(kFALSE),
  fDoMesonAnalysis(kTRUE),
  fDoMesonQA(0),
  fDoClusterQA(0),
  fIsFromMBHeader(kTRUE),
  fIsOverlappingWithOtherHeader(kFALSE),
  fIsMC(0),
  fDoTHnSparse(kTRUE),
  fSetPlotHistsExtQA(kFALSE),
  fWeightJetJetMC(1),
  fDoInOutTimingCluster(kFALSE),
  fMinTimingCluster(0),
  fMaxTimingCluster(0),
  fEnableSortForClusMC(kFALSE),
  fProduceTreeEOverP(kFALSE),
  fProduceCellIDPlots(kFALSE),
  tBrokenFiles(NULL),
  fFileNameBroken(NULL),
  fCloseHighPtClusters(NULL),
  fLocalDebugFlag(0)
{
  // Define output slots here
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskGammaCalo::~AliAnalysisTaskGammaCalo()
{
  if(fClusterCandidates){
    delete fClusterCandidates;
    fClusterCandidates = 0x0;
  }
  if(fBGHandler){
    delete[] fBGHandler;
    fBGHandler = 0x0;
  }
}
//___________________________________________________________
void AliAnalysisTaskGammaCalo::InitBack(){
  
  const Int_t nDim = 4;
  Int_t nBins[nDim] = {800,350,7,6};
  Double_t xMin[nDim] = {0,0, 0,0};
  Double_t xMax[nDim] = {0.8,35,7,6};
  
    if(fDoTHnSparse){
        fSparseMotherInvMassPtZM = new THnSparseF*[fnCuts];
        fSparseMotherBackInvMassPtZM = new THnSparseF*[fnCuts];
    }

  fBGHandler = new AliGammaConversionAODBGHandler*[fnCuts];

  
  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if (((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->DoBGCalculation()){
      TString cutstringEvent  = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
      TString cutstringCalo   = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
      TString cutstringMeson  = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();
      
      Int_t collisionSystem   = atoi((TString)(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber())(0,1));
      Int_t centMin           = atoi((TString)(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber())(1,1));
      Int_t centMax           = atoi((TString)(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber())(2,1));
      
      if(collisionSystem == 1 || collisionSystem == 2 ||
        collisionSystem == 5 || collisionSystem == 8 ||
        collisionSystem == 9){
        centMin   = centMin*10;
        centMax   = centMax*10;
        if(centMax ==0 && centMax!=centMin) centMax=100;
      } else if(collisionSystem == 3 || collisionSystem == 6){
        centMin   = centMin*5;
        centMax   = centMax*5;
      } else if(collisionSystem == 4 || collisionSystem == 7){
        centMin   = ((centMin*5)+45);
        centMax   = ((centMax*5)+45);
      }
      
      if(fDoTHnSparse){
        fBackList[iCut] = new TList();
        fBackList[iCut]->SetName(Form("%s_%s_%s Back histograms",cutstringEvent.Data(),cutstringCalo.Data(), cutstringMeson.Data()));
        fBackList[iCut]->SetOwner(kTRUE);
        fCutFolder[iCut]->Add(fBackList[iCut]);

        fSparseMotherBackInvMassPtZM[iCut] = new THnSparseF("Back_Back_InvMass_Pt_z_m","Back_Back_InvMass_Pt_z_m",nDim,nBins,xMin,xMax);
        fBackList[iCut]->Add(fSparseMotherBackInvMassPtZM[iCut]);

        fMotherList[iCut] = new TList();
        fMotherList[iCut]->SetName(Form("%s_%s_%s Mother histograms",cutstringEvent.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
        fMotherList[iCut]->SetOwner(kTRUE);
        fCutFolder[iCut]->Add(fMotherList[iCut]);

        fSparseMotherInvMassPtZM[iCut] = new THnSparseF("Back_Mother_InvMass_Pt_z_m","Back_Mother_InvMass_Pt_z_m",nDim,nBins,xMin,xMax);
        fMotherList[iCut]->Add(fSparseMotherInvMassPtZM[iCut]);
      }

      if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->BackgroundHandlerType() == 0){
        fBGHandler[iCut] = new AliGammaConversionAODBGHandler(
                                  collisionSystem,centMin,centMax,
                                  ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetNumberOfBGEvents(),
                                  ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseTrackMultiplicity(),
                                  4,8,7);
      }
    }
  }
}
//________________________________________________________________________
void AliAnalysisTaskGammaCalo::UserCreateOutputObjects(){
  
  fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
  if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader

  if (fDoClusterQA == 2) fProduceCellIDPlots = kTRUE;
  if (fIsMC == 2){
    fDoClusterQA      = 0;
    fDoTHnSparse      = kFALSE;
  } else if (fIsMC == 3){
    fDoTHnSparse      = kFALSE;
  }    

  // set common binning in pT for mesons and photons
  Int_t nBinsPt               = 250;
  Float_t minPt               = 0;
  Float_t maxPt               = 25;
  Int_t nBinsQAPt             = 200;
  Float_t minQAPt             = 0.4;
  Float_t maxQAPt             = 25;
  Int_t nBinsClusterPt        = 500;
  Float_t minClusterPt        = 0;
  Float_t maxClusterPt        = 50;
  if (((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::k8TeV ){
    nBinsPt                   = 400;
    minPt                     = 0;
    maxPt                     = 40;
    nBinsQAPt                 = 300;
    minQAPt                   = 0.4;
    maxQAPt                   = 40;
    nBinsClusterPt            = 800;
    minClusterPt              = 0;
    maxClusterPt              = 80;  
  }  

  
  // Create histograms
  if(fOutputContainer != NULL){
    delete fOutputContainer;
    fOutputContainer  = NULL;
  }
  if(fOutputContainer == NULL){
    fOutputContainer  = new TList();
    fOutputContainer->SetOwner(kTRUE);
  }
  
  // Array of current cut's gammas
  fClusterCandidates  = new TList();
  fClusterCandidates->SetOwner(kTRUE);
  
  fCutFolder          = new TList*[fnCuts];
  fESDList            = new TList*[fnCuts];
    if(fDoTHnSparse){
        fBackList     = new TList*[fnCuts];
        fMotherList   = new TList*[fnCuts];
    }
  fHistoNEvents       = new TH1F*[fnCuts];
  if(fIsMC > 1){
    fHistoNEventsWOWeight   = new TH1F*[fnCuts];
  }  
  if(fIsMC == 2){  
    fProfileJetJetXSection  = new TProfile*[fnCuts];
    fHistoJetJetNTrials     = new TH1F*[fnCuts];
  }
  
  fHistoNGoodESDTracks      = new TH1F*[fnCuts];
  fHistoVertexZ             = new TH1F*[fnCuts];
  fHistoNGammaCandidates    = new TH1F*[fnCuts];
  if(!fDoLightOutput){
    fHistoNGoodESDTracksVsNGammaCandidates  = new TH2F*[fnCuts];
    fHistoSPDClusterTrackletBackground      = new TH2F*[fnCuts];
    fHistoNV0Tracks                         = new TH1F*[fnCuts];
  }
  if(fIsHeavyIon==2) fProfileEtaShift          = new TProfile*[fnCuts];
    
  if(fDoMesonAnalysis){
    fHistoMotherInvMassPt           = new TH2F*[fnCuts];
    fHistoMotherBackInvMassPt       = new TH2F*[fnCuts];
    if(!fDoLightOutput){
      fHistoMotherInvMassPtAlpha      = new TH2F*[fnCuts];
      fHistoMotherBackInvMassPtAlpha  = new TH2F*[fnCuts];
    }
    if (fDoMesonQA > 0 && fDoMesonQA < 3){
      fHistoMotherPi0PtY            = new TH2F*[fnCuts];
      fHistoMotherEtaPtY            = new TH2F*[fnCuts];
      fHistoMotherPi0PtAlpha        = new TH2F*[fnCuts];
      fHistoMotherEtaPtAlpha        = new TH2F*[fnCuts];
      fHistoMotherPi0PtOpenAngle    = new TH2F*[fnCuts];
      fHistoMotherEtaPtOpenAngle    = new TH2F*[fnCuts];
    }
  }

  if(fProduceCellIDPlots){
    fHistCellIDvsClusterEnergy      = new TH2F*[fnCuts];
    fHistCellIDvsClusterEnergyMax   = new TH2F*[fnCuts];
  }
    
  fHistoClusGammaPt                 = new TH1F*[fnCuts];
  fHistoClusGammaE                  = new TH1F*[fnCuts];
  fHistoClusOverlapHeadersGammaPt   = new TH1F*[fnCuts];
  if(!fDoLightOutput && fDoClusterQA > 0)
    fHistoClusGammaPtM02            = new TH2F*[fnCuts];

  if (fDoMesonQA == 4 && fIsMC == 0){
    fTreeList                       = new TList*[fnCuts];
    tSigInvMassPtAlphaTheta         = new TTree*[fnCuts];
    tBckInvMassPtAlphaTheta         = new TTree*[fnCuts];
  }

  if (fProduceTreeEOverP){
    fClusterTreeList                = new TList*[fnCuts];
    tClusterEOverP                  = new TTree*[fnCuts];
  }
  
  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    TString cutstringEvent    = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
    TString cutstringCalo     = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
    TString cutstringMeson    = "NoMesonCut";
    if(fDoMesonAnalysis)
      cutstringMeson          = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();
    
    fCutFolder[iCut]        = new TList();
    fCutFolder[iCut]->SetName(Form("Cut Number %s_%s_%s",cutstringEvent.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
    fCutFolder[iCut]->SetOwner(kTRUE);
    fOutputContainer->Add(fCutFolder[iCut]);
    fESDList[iCut]          = new TList();
    fESDList[iCut]->SetName(Form("%s_%s_%s ESD histograms",cutstringEvent.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
    fESDList[iCut]->SetOwner(kTRUE);
    fCutFolder[iCut]->Add(fESDList[iCut]);
    
    fHistoNEvents[iCut]     = new TH1F("NEvents","NEvents",14,-0.5,13.5);
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
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(7,"SPD Pile-Up");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(8,"no SDD");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(9,"no V0AND");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(10,"EMCAL problems");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(11,"rejectedForJetJetMC");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(14,"Pileup V0M-TPCout Tracks");
    fESDList[iCut]->Add(fHistoNEvents[iCut]);
    
    if (fIsMC > 1){
      fHistoNEventsWOWeight[iCut] = new TH1F("NEventsWOWeight","NEventsWOWeight",14,-0.5,13.5);
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
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(11,"rejectedForJetJetMC");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(14,"Pileup V0M-TPCout Tracks");
      fESDList[iCut]->Add(fHistoNEventsWOWeight[iCut]);
    }
    if (fIsMC == 2){  
      fProfileJetJetXSection[iCut]  = new TProfile("XSection","XSection",1,-0.5,0.5);
      fESDList[iCut]->Add(fProfileJetJetXSection[iCut]);
      fHistoJetJetNTrials[iCut]     = new TH1F("NTrials","#sum{NTrials}",1,0,1);
      fHistoJetJetNTrials[iCut]->GetXaxis()->SetBinLabel(1,"#sum{NTrials}");
      fESDList[iCut]->Add(fHistoJetJetNTrials[iCut]);
    }

    if(fIsHeavyIon == 1) 
      fHistoNGoodESDTracks[iCut]    = new TH1F("GoodESDTracks","GoodESDTracks",4000,0,4000);
    else if(fIsHeavyIon == 2) 
      fHistoNGoodESDTracks[iCut]    = new TH1F("GoodESDTracks","GoodESDTracks",400,0,400);
    else 
      fHistoNGoodESDTracks[iCut]    = new TH1F("GoodESDTracks","GoodESDTracks",200,0,200);
    fESDList[iCut]->Add(fHistoNGoodESDTracks[iCut]);

    fHistoVertexZ[iCut]             = new TH1F("VertexZ","VertexZ",1000,-50,50);
    fESDList[iCut]->Add(fHistoVertexZ[iCut]);

    if(fIsHeavyIon == 1) 
      fHistoNGammaCandidates[iCut]  = new TH1F("GammaCandidates","GammaCandidates",600,0,600);
    else if(fIsHeavyIon == 2) 
      fHistoNGammaCandidates[iCut]  = new TH1F("GammaCandidates","GammaCandidates",400,0,400);
    else 
      fHistoNGammaCandidates[iCut]  = new TH1F("GammaCandidates","GammaCandidates",100,0,100);
    fESDList[iCut]->Add(fHistoNGammaCandidates[iCut]);

    if(!fDoLightOutput){
      if(fIsHeavyIon == 1)
        fHistoNGoodESDTracksVsNGammaCandidates[iCut]  = new TH2F("GoodESDTracksVsGammaCandidates","GoodESDTracksVsGammaCandidates",4000,0,4000,600,0,600);
      else if(fIsHeavyIon == 2)
        fHistoNGoodESDTracksVsNGammaCandidates[iCut]  = new TH2F("GoodESDTracksVsGammaCandidates","GoodESDTracksVsGammaCandidates",400,0,400,400,0,400);
      else
        fHistoNGoodESDTracksVsNGammaCandidates[iCut]  = new TH2F("GoodESDTracksVsGammaCandidates","GoodESDTracksVsGammaCandidates",200,0,200,100,0,100);
      fESDList[iCut]->Add(fHistoNGoodESDTracksVsNGammaCandidates[iCut]);

      fHistoSPDClusterTrackletBackground[iCut]        = new TH2F("SPD tracklets vs SPD clusters","SPD tracklets vs SPD clusters",100,0,200,250,0,1000);
      fESDList[iCut]->Add(fHistoSPDClusterTrackletBackground[iCut]);

      if(fIsHeavyIon == 1)
        fHistoNV0Tracks[iCut]       = new TH1F("V0 Multiplicity","V0 Multiplicity",30000,0,30000);
      else if(fIsHeavyIon == 2)
        fHistoNV0Tracks[iCut]       = new TH1F("V0 Multiplicity","V0 Multiplicity",2500,0,2500);
      else
        fHistoNV0Tracks[iCut]       = new TH1F("V0 Multiplicity","V0 Multiplicity",1500,0,1500);
      fESDList[iCut]->Add(fHistoNV0Tracks[iCut]);
    }

    if(fIsHeavyIon==2) {
      fProfileEtaShift[iCut]        = new TProfile("Eta Shift","Eta Shift",1, -0.5,0.5);
      fESDList[iCut]->Add(fProfileEtaShift[iCut]);
    }

    if (fIsMC > 1){
      fHistoNEvents[iCut]->Sumw2();
      fHistoNGoodESDTracks[iCut]->Sumw2();
      fHistoVertexZ[iCut]->Sumw2();
      fHistoNGammaCandidates[iCut]->Sumw2();
      if(!fDoLightOutput){
        fHistoNGoodESDTracksVsNGammaCandidates[iCut]->Sumw2();
        fHistoSPDClusterTrackletBackground[iCut]->Sumw2();
        fHistoNV0Tracks[iCut]->Sumw2();
      }
      if(fIsHeavyIon==2) fProfileEtaShift[iCut]->Sumw2();
    }

    fHistoClusGammaPt[iCut]               = new TH1F("ClusGamma_Pt","ClusGamma_Pt",nBinsClusterPt, minClusterPt, maxClusterPt);
    fESDList[iCut]->Add(fHistoClusGammaPt[iCut]);
    fHistoClusGammaE[iCut]               = new TH1F("ClusGamma_E","ClusGamma_E",nBinsClusterPt, minClusterPt, maxClusterPt);
    fESDList[iCut]->Add(fHistoClusGammaE[iCut]);
    fHistoClusOverlapHeadersGammaPt[iCut] = new TH1F("ClusGammaOverlapHeaders_Pt","ClusGammaOverlapHeaders_Pt",nBinsClusterPt, minClusterPt, maxClusterPt);
    fESDList[iCut]->Add(fHistoClusOverlapHeadersGammaPt[iCut]);
    if(!fDoLightOutput && fDoClusterQA > 0){
      fHistoClusGammaPtM02[iCut]               = new TH2F("ClusGamma_Pt_M02","ClusGamma_Pt_M02",nBinsClusterPt, minClusterPt, maxClusterPt,100,0,1);
      fESDList[iCut]->Add(fHistoClusGammaPtM02[iCut]);
    }
    
    if (fIsMC > 1){
      fHistoClusGammaPt[iCut]->Sumw2();
      fHistoClusGammaE[iCut]->Sumw2();
      fHistoClusOverlapHeadersGammaPt[iCut]->Sumw2();
      if(!fDoLightOutput && fDoClusterQA > 0)fHistoClusGammaPtM02[iCut]->Sumw2();
    }
    
    if(fDoMesonAnalysis){
      fHistoMotherInvMassPt[iCut]           = new TH2F("ESD_Mother_InvMass_Pt","ESD_Mother_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
      fESDList[iCut]->Add(fHistoMotherInvMassPt[iCut]);
      fHistoMotherBackInvMassPt[iCut]       = new TH2F("ESD_Background_InvMass_Pt","ESD_Background_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
      fESDList[iCut]->Add(fHistoMotherBackInvMassPt[iCut]);
      if(!fDoLightOutput){
        fHistoMotherInvMassPtAlpha[iCut]      = new TH2F("ESD_Mother_InvMass_vs_Pt_Alpha","ESD_Mother_InvMass_vs_Pt_Alpha",800,0,0.8,nBinsPt, minPt, maxPt);
        fESDList[iCut]->Add(fHistoMotherInvMassPtAlpha[iCut]);
        fHistoMotherBackInvMassPtAlpha[iCut]  = new TH2F("ESD_Background_InvMass_vs_Pt_Alpha","ESD_Background_InvMass_vs_Pt_Alpha",800,0,0.8,nBinsPt, minPt, maxPt);
        fESDList[iCut]->Add(fHistoMotherBackInvMassPtAlpha[iCut]);
      }

      if (fIsMC > 1){
        fHistoMotherInvMassPt[iCut]->Sumw2();
        fHistoMotherBackInvMassPt[iCut]->Sumw2();
        if(!fDoLightOutput){
          fHistoMotherInvMassPtAlpha[iCut]->Sumw2();
          fHistoMotherBackInvMassPtAlpha[iCut]->Sumw2();
        }
      }
      
      if (fDoMesonQA > 0 && fDoMesonQA < 3 ){
        fHistoMotherPi0PtY[iCut]          = new TH2F("ESD_MotherPi0_Pt_Y","ESD_MotherPi0_Pt_Y",nBinsQAPt, minQAPt, maxQAPt,150,-1.5,1.5);
        SetLogBinningXTH2(fHistoMotherPi0PtY[iCut]);
        fESDList[iCut]->Add(fHistoMotherPi0PtY[iCut]);
        fHistoMotherEtaPtY[iCut]          = new TH2F("ESD_MotherEta_Pt_Y","ESD_MotherEta_Pt_Y",nBinsQAPt, minQAPt, maxQAPt,150,-1.5,1.5);
        SetLogBinningXTH2(fHistoMotherEtaPtY[iCut]);
        fESDList[iCut]->Add(fHistoMotherEtaPtY[iCut]);
        fHistoMotherPi0PtAlpha[iCut]      = new TH2F("ESD_MotherPi0_Pt_Alpha","ESD_MotherPi0_Pt_Alpha",nBinsQAPt, minQAPt, maxQAPt,100,0,1);
        SetLogBinningXTH2(fHistoMotherPi0PtAlpha[iCut]);
        fESDList[iCut]->Add(fHistoMotherPi0PtAlpha[iCut]);
        fHistoMotherEtaPtAlpha[iCut]      = new TH2F("ESD_MotherEta_Pt_Alpha","ESD_MotherEta_Pt_Alpha",nBinsQAPt, minQAPt, maxQAPt,100,0,1);
        SetLogBinningXTH2(fHistoMotherEtaPtAlpha[iCut]);
        fESDList[iCut]->Add(fHistoMotherEtaPtAlpha[iCut]);
        fHistoMotherPi0PtOpenAngle[iCut]  = new TH2F("ESD_MotherPi0_Pt_OpenAngle","ESD_MotherPi0_Pt_OpenAngle",nBinsQAPt, minQAPt, maxQAPt,100,0, 0.5);
        SetLogBinningXTH2(fHistoMotherPi0PtOpenAngle[iCut]);
        fESDList[iCut]->Add(fHistoMotherPi0PtOpenAngle[iCut]);
        fHistoMotherEtaPtOpenAngle[iCut]  = new TH2F("ESD_MotherEta_Pt_OpenAngle","ESD_MotherEta_Pt_OpenAngle",nBinsQAPt, minQAPt, maxQAPt,180,0, 1.8);
        SetLogBinningXTH2(fHistoMotherEtaPtOpenAngle[iCut]);
        fESDList[iCut]->Add(fHistoMotherEtaPtOpenAngle[iCut]);
      }
      
      if (fIsMC > 1 && fDoMesonQA > 0 && fDoMesonQA < 3){
        fHistoMotherPi0PtY[iCut]->Sumw2();
        fHistoMotherEtaPtY[iCut]->Sumw2();
        fHistoMotherPi0PtAlpha[iCut]->Sumw2();
        fHistoMotherEtaPtAlpha[iCut]->Sumw2();
        fHistoMotherPi0PtOpenAngle[iCut]->Sumw2();
        fHistoMotherEtaPtOpenAngle[iCut]->Sumw2();
      }

      if (fProduceCellIDPlots){
        Int_t nMaxCells      = 12*48*24;
        if(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetClusterType() == 2) nMaxCells = 5*56*64;
        fHistCellIDvsClusterEnergy[iCut]         = new TH2F("CellIDvsClusterEnergy","CellIDvsClusterEnergy",100,0.5,100.,nMaxCells,0,nMaxCells);
        SetLogBinningXTH2(fHistCellIDvsClusterEnergy[iCut]);
        fESDList[iCut]->Add(fHistCellIDvsClusterEnergy[iCut]);
        fHistCellIDvsClusterEnergyMax[iCut]      = new TH2F("CellIDvsClusterEnergyMax","CellIDvsClusterEnergyMax",100,0.5,100.,nMaxCells,0,nMaxCells);
        SetLogBinningXTH2(fHistCellIDvsClusterEnergyMax[iCut]);
        fESDList[iCut]->Add(fHistCellIDvsClusterEnergyMax[iCut]);
      }

      if (fDoMesonQA == 4 && fIsMC == 0){
        fTreeList[iCut] = new TList();
        fTreeList[iCut]->SetName(Form("%s_%s_%s InvMass Tree",cutstringEvent.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
        fTreeList[iCut]->SetOwner(kTRUE);
        fCutFolder[iCut]->Add(fTreeList[iCut]);

        tSigInvMassPtAlphaTheta[iCut] = new TTree("Sig_InvMass_Pt_Alpha_Theta_MixPool","Sig_InvMass_Pt_Alpha_Theta_MixPool");
        tSigInvMassPtAlphaTheta[iCut]->Branch("InvMass",&fInvMassTreeInvMass,"fInvMassTreeInvMass/F");
        tSigInvMassPtAlphaTheta[iCut]->Branch("Pt",&fInvMassTreePt,"fInvMassTreePt/F");
        tSigInvMassPtAlphaTheta[iCut]->Branch("Alpha",&fInvMassTreeAlpha,"fInvMassTreeAlpha/F");
        tSigInvMassPtAlphaTheta[iCut]->Branch("Theta",&fInvMassTreeTheta,"fInvMassTreeTheta/F");
        tSigInvMassPtAlphaTheta[iCut]->Branch("MixPool",&fInvMassTreeMixPool,"fInvMassTreeMixPool/I");
        tSigInvMassPtAlphaTheta[iCut]->Branch("zVtx",&fInvMassTreeZVertex,"fInvMassTreeZVertex/F");
        tSigInvMassPtAlphaTheta[iCut]->Branch("Eta",&fInvMassTreeEta,"fInvMassTreeEta/F");
        fTreeList[iCut]->Add(tSigInvMassPtAlphaTheta[iCut]);

        tBckInvMassPtAlphaTheta[iCut] = new TTree("Bck_InvMass_Pt_Alpha_Theta_MixPool","Bck_InvMass_Pt_Alpha_Theta_MixPool");
        tBckInvMassPtAlphaTheta[iCut]->Branch("InvMass",&fInvMassTreeInvMass,"fInvMassTreeInvMass/F");
        tBckInvMassPtAlphaTheta[iCut]->Branch("Pt",&fInvMassTreePt,"fInvMassTreePt/F");
        tBckInvMassPtAlphaTheta[iCut]->Branch("Alpha",&fInvMassTreeAlpha,"fInvMassTreeAlpha/F");
        tBckInvMassPtAlphaTheta[iCut]->Branch("Theta",&fInvMassTreeTheta,"fInvMassTreeTheta/F");
        tBckInvMassPtAlphaTheta[iCut]->Branch("MixPool",&fInvMassTreeMixPool,"fInvMassTreeMixPool/I");
        tBckInvMassPtAlphaTheta[iCut]->Branch("zVtx",&fInvMassTreeZVertex,"fInvMassTreeZVertex/F");
        tBckInvMassPtAlphaTheta[iCut]->Branch("Eta",&fInvMassTreeEta,"fInvMassTreeEta/F");
        fTreeList[iCut]->Add(tBckInvMassPtAlphaTheta[iCut]);
      }

      if (fProduceTreeEOverP ){
        fClusterTreeList[iCut] = new TList();
        fClusterTreeList[iCut]->SetName(Form("%s_%s EoverP Tree",cutstringEvent.Data(),cutstringCalo.Data()));
        fClusterTreeList[iCut]->SetOwner(kTRUE);
        fCutFolder[iCut]->Add(fClusterTreeList[iCut]);

        tClusterEOverP[iCut] = new TTree("EOverP_ClusE_ClusM02_ClusM20_TrackP_TrackPt","EOverP_ClusE_ClusM02_ClusM20_TrackP_TrackPt");
        tClusterEOverP[iCut]->Branch("ClusE",&fClusterE,"fClusterE/F");
        tClusterEOverP[iCut]->Branch("ClusM02",&fClusterM02,"fClusterM02/F");
        tClusterEOverP[iCut]->Branch("ClusM20",&fClusterM20,"fClusterM20/F");
        tClusterEOverP[iCut]->Branch("ClusEP",&fClusterEP,"fClusterEP/F");
        tClusterEOverP[iCut]->Branch("ClusLeadCellID",&fClusterLeadCellID,"fClusterLeadCellID/I");
        if(fIsMC > 0) tClusterEOverP[iCut]->Branch("ClusClassification",&fClusterClassification,"fClusterClassification/I");
        tClusterEOverP[iCut]->Branch("ClusTrackDeltaEta",&fDeltaEta,"fDeltaEta/F");
        tClusterEOverP[iCut]->Branch("ClusTrackDeltaPhi",&fDeltaPhi,"fDeltaPhi/F");
        tClusterEOverP[iCut]->Branch("TrackPt",&fTrackPt,"fTrackPt/F");
        tClusterEOverP[iCut]->Branch("TrackPID_e",&fTrackPID_e,"fTrackPID_e/I");
        tClusterEOverP[iCut]->Branch("TrackPID_Pi",&fTrackPID_Pi,"fTrackPID_Pi/I");
        tClusterEOverP[iCut]->Branch("TrackPID_K",&fTrackPID_K,"fTrackPID_K/I");
        tClusterEOverP[iCut]->Branch("TrackPID_P",&fTrackPID_P,"fTrackPID_P/I");
        tClusterEOverP[iCut]->Branch("ClusIsoSumClusEt",&fClusterIsoSumClusterEt,"fClusterIsoSumClusterEt/F");
        tClusterEOverP[iCut]->Branch("ClusIsoSumTrackEt",&fClusterIsoSumTrackEt,"fClusterIsoSumTrackEt/F");
        fClusterTreeList[iCut]->Add(tClusterEOverP[iCut]);
      }
    }
  }
  if(fDoMesonAnalysis){
    InitBack(); // Init Background Handler
  }
  
  if(fIsMC> 0){
    // MC Histogramms
    fMCList     = new TList*[fnCuts];
    // True Histogramms
    fTrueList   = new TList*[fnCuts];
    // Selected Header List
    if (fDoMesonQA ==3){
      fTreeList                     = new TList*[fnCuts];
      tTrueInvMassROpenABPtFlag     = new TTree*[fnCuts];
    }

    if(!fDoLightOutput){
      fHistoMCHeaders                 = new TH1I*[fnCuts];
      fHistoMCAllGammaPt              = new TH1F*[fnCuts];
      fHistoMCAllSecondaryGammaPt     = new TH2F*[fnCuts];
      fHistoMCDecayGammaPi0Pt         = new TH1F*[fnCuts];
      fHistoMCDecayGammaRhoPt         = new TH1F*[fnCuts];
      fHistoMCDecayGammaEtaPt         = new TH1F*[fnCuts];
      fHistoMCDecayGammaOmegaPt       = new TH1F*[fnCuts];
      fHistoMCDecayGammaEtapPt        = new TH1F*[fnCuts];
      fHistoMCDecayGammaPhiPt         = new TH1F*[fnCuts];
      fHistoMCDecayGammaSigmaPt       = new TH1F*[fnCuts];
      fHistoClusPhotonBGPt            = new TH2F*[fnCuts];
      fHistoClusPhotonPlusConvBGPt    = new TH2F*[fnCuts];
      if (fDoClusterQA > 0) {
        fHistoClustPhotonElectronBGPtM02          = new TH2F*[fnCuts];
        fHistoClustPhotonPionBGPtM02              = new TH2F*[fnCuts];
        fHistoClustPhotonKaonBGPtM02              = new TH2F*[fnCuts];
        fHistoClustPhotonK0lBGPtM02               = new TH2F*[fnCuts];
        fHistoClustPhotonNeutronBGPtM02           = new TH2F*[fnCuts];
        fHistoClustPhotonRestBGPtM02              = new TH2F*[fnCuts];
        fHistoClustPhotonPlusConvElectronBGPtM02  = new TH2F*[fnCuts];
        fHistoClustPhotonPlusConvPionBGPtM02      = new TH2F*[fnCuts];
        fHistoClustPhotonPlusConvKaonBGPtM02      = new TH2F*[fnCuts];
        fHistoClustPhotonPlusConvK0lBGPtM02       = new TH2F*[fnCuts];
        fHistoClustPhotonPlusConvNeutronBGPtM02   = new TH2F*[fnCuts];
        fHistoClustPhotonPlusConvRestBGPtM02      = new TH2F*[fnCuts];
      }
    }
  
    fHistoTrueClusGammaPt                                       = new TH1F*[fnCuts];
    if(!fDoLightOutput){
      if (fDoClusterQA > 0) fHistoTrueClusGammaPtM02            = new TH2F*[fnCuts];
      fHistoTruePrimaryClusGammaPt                              = new TH1F*[fnCuts];
      fHistoTruePrimaryClusGammaESDPtMCPt                       = new TH2F*[fnCuts];
      fHistoTruePrimaryClusConvGammaPt                          = new TH1F*[fnCuts];
      fHistoTruePrimaryClusConvGammaESDPtMCPt                   = new TH2F*[fnCuts];
      fHistoTrueSecondaryClusGammaPt                            = new TH2F*[fnCuts];
      fHistoTrueSecondaryClusConvGammaPt                        = new TH2F*[fnCuts];
      fHistoTrueSecondaryClusGammaMCPt                          = new TH2F*[fnCuts];
      fHistoTrueSecondaryClusConvGammaMCPt                      = new TH2F*[fnCuts];
      fHistoTrueSecondaryClusGammaFromXFromK0sMCPtESDPt         = new TH2F*[fnCuts];
      fHistoTrueSecondaryClusConvGammaFromXFromK0sMCPtESDPt     = new TH2F*[fnCuts];
      fHistoTrueSecondaryClusGammaFromXFromK0lMCPtESDPt         = new TH2F*[fnCuts];
      fHistoTrueSecondaryClusConvGammaFromXFromK0lMCPtESDPt     = new TH2F*[fnCuts];
      fHistoTrueSecondaryClusGammaFromXFromLambdaMCPtESDPt      = new TH2F*[fnCuts];
      fHistoTrueSecondaryClusConvGammaFromXFromLambdaMCPtESDPt  = new TH2F*[fnCuts];
      fHistoTrueNLabelsInClus                                   = new TH1F*[fnCuts];
    }
    fHistoDoubleCountTrueClusterGammaPt                         = new TH2F*[fnCuts];
    fHistoMultipleCountTrueClusterGamma                         = new TH1F*[fnCuts];

//    fHistoTruePi0NonLinearity      = new TH2F*[fnCuts];
//    fHistoTrueEtaNonLinearity      = new TH2F*[fnCuts];

    if (fDoClusterQA > 0){
      fHistoTrueClusUnConvGammaPt           = new TH1F*[fnCuts];
      fHistoTrueClusUnConvGammaMCPt         = new TH1F*[fnCuts];
      if (!fDoLightOutput)
        fHistoTrueClusUnConvGammaPtM02      = new TH2F*[fnCuts];
      fHistoTrueClusElectronPt              = new TH1F*[fnCuts];
      fHistoTrueClusConvGammaPt             = new TH1F*[fnCuts];
      fHistoTrueClusConvGammaMCPt           = new TH1F*[fnCuts];
      fHistoTrueClusConvGammaFullyPt        = new TH1F*[fnCuts];
      fHistoTrueClusMergedGammaPt           = new TH1F*[fnCuts];
      fHistoTrueClusMergedPartConvGammaPt   = new TH1F*[fnCuts];
      fHistoTrueClusDalitzPt                = new TH1F*[fnCuts];
      fHistoTrueClusDalitzMergedPt          = new TH1F*[fnCuts];
      fHistoTrueClusPhotonFromElecMotherPt  = new TH1F*[fnCuts];
      fHistoTrueClusShowerPt                = new TH1F*[fnCuts];
      fHistoTrueClusSubLeadingPt            = new TH1F*[fnCuts];
      fHistoTrueClusNParticles              = new TH1F*[fnCuts];
      fHistoTrueClusEMNonLeadingPt          = new TH1F*[fnCuts];
    }
    
    if(fDoMesonAnalysis){
      fHistoMCPi0Pt             = new TH1F*[fnCuts];
      fHistoMCPi0WOWeightPt     = new TH1F*[fnCuts];
      fHistoMCEtaPt             = new TH1F*[fnCuts];
      fHistoMCEtaWOWeightPt     = new TH1F*[fnCuts];
      fHistoMCPi0InAccPt        = new TH1F*[fnCuts];
      fHistoMCEtaInAccPt        = new TH1F*[fnCuts];
      
      if (fIsMC > 1){
        fHistoMCPi0WOEvtWeightPt       = new TH1F*[fnCuts];
        fHistoMCEtaWOEvtWeightPt       = new TH1F*[fnCuts];
        fHistoMCPi0WOEvtWeightInAccPt  = new TH1F*[fnCuts];
        fHistoMCEtaWOEvtWeightInAccPt  = new TH1F*[fnCuts];
      }
      
      fHistoTruePi0InvMassPt                    = new TH2F*[fnCuts];
      fHistoTrueEtaInvMassPt                    = new TH2F*[fnCuts];
      fHistoDoubleCountTruePi0InvMassPt         = new TH2F*[fnCuts];
      fHistoDoubleCountTrueEtaInvMassPt         = new TH2F*[fnCuts];
      fHistoTruePrimaryPi0InvMassPt             = new TH2F*[fnCuts];
      fHistoTruePrimaryEtaInvMassPt             = new TH2F*[fnCuts];
      fHistoTruePrimaryPi0W0WeightingInvMassPt  = new TH2F*[fnCuts];
      fHistoTruePrimaryEtaW0WeightingInvMassPt  = new TH2F*[fnCuts];
      fProfileTruePrimaryPi0WeightsInvMassPt    = new TProfile2D*[fnCuts];
      fProfileTruePrimaryEtaWeightsInvMassPt    = new TProfile2D*[fnCuts];
      fHistoTrueSecondaryPi0InvMassPt           = new TH2F*[fnCuts];
      fHistoTrueSecondaryPi0FromK0sInvMassPt    = new TH2F*[fnCuts];
      fHistoTrueSecondaryPi0FromK0lInvMassPt    = new TH2F*[fnCuts];
      fHistoTrueSecondaryPi0FromEtaInvMassPt    = new TH2F*[fnCuts];
      fHistoTrueSecondaryPi0FromLambdaInvMassPt = new TH2F*[fnCuts];
      if(!fDoLightOutput){
        fHistoTruePi0InvMassPtAlpha               = new TH2F*[fnCuts];
        fHistoTruePi0PureGammaInvMassPtAlpha      = new TH2F*[fnCuts];
      }
      fHistoMCPrimaryPtvsSource                 = new TH2F*[fnCuts];
      fHistoMCSecPi0PtvsSource                  = new TH2F*[fnCuts];
      fHistoMCSecPi0InAccPtvsSource             = new TH2F*[fnCuts];
      fHistoMCSecPi0Source                      = new TH1F*[fnCuts];
      fHistoMCSecEtaPt                          = new TH1F*[fnCuts];
      fHistoMCSecEtaSource                      = new TH1F*[fnCuts];

      if (fDoMesonQA > 0 && fDoMesonQA < 3 ){
        fHistoMCPi0PtY                    = new TH2F*[fnCuts];
        fHistoMCEtaPtY                    = new TH2F*[fnCuts];
        fHistoMCPi0PtAlpha                = new TH2F*[fnCuts];
        fHistoMCEtaPtAlpha                = new TH2F*[fnCuts];
        if (fIsMC == 2){
          fHistoMCPi0PtJetPt              = new TH2F*[fnCuts];
          fHistoMCEtaPtJetPt              = new TH2F*[fnCuts];
        }

        if (fIsMC < 2){
          fHistoTruePi0CaloPhotonInvMassPt                = new TH2F*[fnCuts];
          fHistoTrueEtaCaloPhotonInvMassPt                = new TH2F*[fnCuts];
          fHistoTruePi0CaloMixedPhotonConvPhotonInvMassPt = new TH2F*[fnCuts];
          fHistoTrueEtaCaloMixedPhotonConvPhotonInvMassPt = new TH2F*[fnCuts];
          fHistoTruePi0CaloConvertedPhotonInvMassPt       = new TH2F*[fnCuts];
          fHistoTrueEtaCaloConvertedPhotonInvMassPt       = new TH2F*[fnCuts];
          fHistoTruePi0CaloElectronInvMassPt              = new TH2F*[fnCuts];
          fHistoTrueEtaCaloElectronInvMassPt              = new TH2F*[fnCuts];
          fHistoTruePi0CaloMergedClusterInvMassPt         = new TH2F*[fnCuts];
          fHistoTrueEtaCaloMergedClusterInvMassPt         = new TH2F*[fnCuts];
          fHistoTruePi0CaloMergedClusterPartConvInvMassPt = new TH2F*[fnCuts];
          fHistoTrueEtaCaloMergedClusterPartConvInvMassPt = new TH2F*[fnCuts];
          fHistoTruePi0NonMergedElectronPhotonInvMassPt   = new TH2F*[fnCuts];
          fHistoTruePi0NonMergedElectronMergedPhotonInvMassPt   = new TH2F*[fnCuts];
          fHistoTruePrimaryPi0MCPtResolPt                 = new TH2F*[fnCuts];
          fHistoTruePrimaryEtaMCPtResolPt                 = new TH2F*[fnCuts];
          fHistoTrueK0sWithPi0DaughterMCPt                = new TH1F*[fnCuts];
          fHistoTrueK0lWithPi0DaughterMCPt                = new TH1F*[fnCuts];
          fHistoTrueEtaWithPi0DaughterMCPt                = new TH1F*[fnCuts];
          fHistoTrueLambdaWithPi0DaughterMCPt             = new TH1F*[fnCuts];
        }
        fHistoTruePi0PtY                  = new TH2F*[fnCuts];
        fHistoTrueEtaPtY                  = new TH2F*[fnCuts];
        fHistoTruePi0PtAlpha              = new TH2F*[fnCuts];
        fHistoTrueEtaPtAlpha              = new TH2F*[fnCuts];
        fHistoTruePi0PtOpenAngle          = new TH2F*[fnCuts];
        fHistoTrueEtaPtOpenAngle          = new TH2F*[fnCuts];
        fHistoTrueBckGGInvMassPt          = new TH2F*[fnCuts];
        fHistoTrueBckFullMesonContainedInOneClusterInvMassPt = new TH2F*[fnCuts];
        fHistoTrueBckAsymEClustersInvMassPt                  = new TH2F*[fnCuts];
        fHistoTrueBckContInvMassPt        = new TH2F*[fnCuts];
      }
      if (fDoMesonQA==2){
        fHistoTruePi0Category1            = new TH2F*[fnCuts];
        fHistoTrueEtaCategory1            = new TH2F*[fnCuts];
        fHistoTruePi0Category2            = new TH2F*[fnCuts];
        fHistoTrueEtaCategory2            = new TH2F*[fnCuts];
        fHistoTruePi0Category3            = new TH2F*[fnCuts];
        fHistoTrueEtaCategory3            = new TH2F*[fnCuts];
        fHistoTruePi0Category4_6          = new TH2F*[fnCuts];
        fHistoTrueEtaCategory4_6          = new TH2F*[fnCuts];
        fHistoTruePi0Category5            = new TH2F*[fnCuts];
        fHistoTrueEtaCategory5            = new TH2F*[fnCuts];
        fHistoTruePi0Category7            = new TH2F*[fnCuts];
        fHistoTrueEtaCategory7            = new TH2F*[fnCuts];
        fHistoTruePi0Category8            = new TH2F*[fnCuts];
        fHistoTrueEtaCategory8            = new TH2F*[fnCuts];
      }
    }
    
    
    
    for(Int_t iCut = 0; iCut<fnCuts;iCut++){
      TString cutstringEvent    = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
      TString cutstringCalo     = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
      TString cutstringMeson    = "NoMesonCut";
      if(fDoMesonAnalysis)
        cutstringMeson          = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();
      
      fMCList[iCut]                   = new TList();
      fMCList[iCut]->SetName(Form("%s_%s_%s MC histograms",cutstringEvent.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
      fMCList[iCut]->SetOwner(kTRUE);
      fCutFolder[iCut]->Add(fMCList[iCut]);

      if(!fDoLightOutput){
        fHistoMCHeaders[iCut]             = new TH1I("MC_Headers","MC_Headers",20,0,20);
        fMCList[iCut]->Add(fHistoMCHeaders[iCut]);
        fHistoMCAllGammaPt[iCut]          = new TH1F("MC_AllGamma_Pt","MC_AllGamma_Pt",nBinsClusterPt, minClusterPt, maxClusterPt);
        fMCList[iCut]->Add(fHistoMCAllGammaPt[iCut]);
        fHistoMCAllSecondaryGammaPt[iCut] = new TH2F("MC_AllSecondaryGamma_Pt","MC_AllSecondaryGamma_Pt",nBinsClusterPt, minClusterPt, maxClusterPt, 5, -0.5, 4.5);
        fHistoMCAllSecondaryGammaPt[iCut]->GetYaxis()->SetBinLabel(1,"K0s");
        fHistoMCAllSecondaryGammaPt[iCut]->GetYaxis()->SetBinLabel(2,"K0l");
        fHistoMCAllSecondaryGammaPt[iCut]->GetYaxis()->SetBinLabel(3,"Lambda");
        fHistoMCAllSecondaryGammaPt[iCut]->GetYaxis()->SetBinLabel(4,"Eta");
        fHistoMCAllSecondaryGammaPt[iCut]->GetYaxis()->SetBinLabel(5,"rest");
        fMCList[iCut]->Add(fHistoMCAllSecondaryGammaPt[iCut]);
        fHistoMCDecayGammaPi0Pt[iCut]     = new TH1F("MC_DecayGammaPi0_Pt","MC_DecayGammaPi0_Pt",nBinsClusterPt, minClusterPt, maxClusterPt);
        fMCList[iCut]->Add(fHistoMCDecayGammaPi0Pt[iCut]);
        fHistoMCDecayGammaRhoPt[iCut]     = new TH1F("MC_DecayGammaRho_Pt","MC_DecayGammaRho_Pt",nBinsClusterPt, minClusterPt, maxClusterPt);
        fMCList[iCut]->Add(fHistoMCDecayGammaRhoPt[iCut]);
        fHistoMCDecayGammaEtaPt[iCut]     = new TH1F("MC_DecayGammaEta_Pt","MC_DecayGammaEta_Pt",nBinsClusterPt, minClusterPt, maxClusterPt);
        fMCList[iCut]->Add(fHistoMCDecayGammaEtaPt[iCut]);
        fHistoMCDecayGammaOmegaPt[iCut]   = new TH1F("MC_DecayGammaOmega_Pt","MC_DecayGammaOmmega_Pt",nBinsClusterPt, minClusterPt, maxClusterPt);
        fMCList[iCut]->Add(fHistoMCDecayGammaOmegaPt[iCut]);
        fHistoMCDecayGammaEtapPt[iCut]    = new TH1F("MC_DecayGammaEtap_Pt","MC_DecayGammaEtap_Pt",nBinsClusterPt, minClusterPt, maxClusterPt);
        fMCList[iCut]->Add(fHistoMCDecayGammaEtapPt[iCut]);
        fHistoMCDecayGammaPhiPt[iCut]     = new TH1F("MC_DecayGammaPhi_Pt","MC_DecayGammaPhi_Pt",nBinsClusterPt, minClusterPt, maxClusterPt);
        fMCList[iCut]->Add(fHistoMCDecayGammaPhiPt[iCut]);
        fHistoMCDecayGammaSigmaPt[iCut]   = new TH1F("MC_DecayGammaSigma_Pt","MC_DecayGammaSigma_Pt",nBinsClusterPt, minClusterPt, maxClusterPt);
        fMCList[iCut]->Add(fHistoMCDecayGammaSigmaPt[iCut]);

        if (fIsMC > 1){
          fHistoMCAllGammaPt[iCut]->Sumw2();
          fHistoMCAllSecondaryGammaPt[iCut]->Sumw2();
          fHistoMCDecayGammaPi0Pt[iCut]->Sumw2();
          fHistoMCDecayGammaRhoPt[iCut]->Sumw2();
          fHistoMCDecayGammaEtaPt[iCut]->Sumw2();
          fHistoMCDecayGammaOmegaPt[iCut]->Sumw2();
          fHistoMCDecayGammaEtapPt[iCut]->Sumw2();
          fHistoMCDecayGammaPhiPt[iCut]->Sumw2();
          fHistoMCDecayGammaSigmaPt[iCut]->Sumw2();
        }
      }

      if(fDoMesonAnalysis){
        fHistoMCPi0Pt[iCut]           = new TH1F("MC_Pi0_Pt","MC_Pi0_Pt",nBinsPt, minPt, maxPt);
        fHistoMCPi0Pt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0Pt[iCut]);
        fHistoMCPi0WOWeightPt[iCut]   = new TH1F("MC_Pi0_WOWeights_Pt","MC_Pi0_WOWeights_Pt",nBinsPt, minPt, maxPt);
        fHistoMCPi0WOWeightPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0WOWeightPt[iCut]);
        
        fHistoMCEtaPt[iCut]           = new TH1F("MC_Eta_Pt","MC_Eta_Pt",nBinsPt, minPt, maxPt);
        fHistoMCEtaPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCEtaPt[iCut]);
        fHistoMCEtaWOWeightPt[iCut]   = new TH1F("MC_Eta_WOWeights_Pt","MC_Eta_WOWeights_Pt",nBinsPt, minPt, maxPt);
        fHistoMCEtaWOWeightPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCEtaWOWeightPt[iCut]);
        fHistoMCPi0InAccPt[iCut]      = new TH1F("MC_Pi0InAcc_Pt","MC_Pi0InAcc_Pt",nBinsPt, minPt, maxPt);
        fHistoMCPi0InAccPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0InAccPt[iCut]);
        fHistoMCEtaInAccPt[iCut]      = new TH1F("MC_EtaInAcc_Pt","MC_EtaInAcc_Pt",nBinsPt, minPt, maxPt);
        fHistoMCEtaInAccPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCEtaInAccPt[iCut]);
        if (fIsMC > 1){
          fHistoMCPi0WOWeightPt[iCut]->Sumw2();
          fHistoMCEtaWOWeightPt[iCut]->Sumw2();
          fHistoMCPi0WOEvtWeightPt[iCut] = new TH1F("MC_Pi0_WOEventWeights_Pt","MC_Pi0_WOEventWeights_Pt",nBinsPt, minPt, maxPt);
          fMCList[iCut]->Add(fHistoMCPi0WOEvtWeightPt[iCut]);
          fHistoMCEtaWOEvtWeightPt[iCut] = new TH1F("MC_Eta_WOEventWeights_Pt","MC_Eta_WOEventWeights_Pt",nBinsPt, minPt, maxPt);
          fMCList[iCut]->Add(fHistoMCEtaWOEvtWeightPt[iCut]);
          fHistoMCPi0WOEvtWeightInAccPt[iCut] = new TH1F("MC_Pi0WOEvtWeightInAcc_Pt","MC_Pi0WOEvtWeightInAcc_Pt",nBinsPt, minPt, maxPt);
          fMCList[iCut]->Add(fHistoMCPi0WOEvtWeightInAccPt[iCut]);
          fHistoMCEtaWOEvtWeightInAccPt[iCut] = new TH1F("MC_EtaWOEvtWeightInAcc_Pt","MC_EtaWOEvtWeightInAcc_Pt",nBinsPt, minPt, maxPt);
          fMCList[iCut]->Add(fHistoMCEtaWOEvtWeightInAccPt[iCut]);
          if (fDoMesonQA > 0  && fDoMesonQA < 3 && fIsMC == 2){
            fHistoMCPi0PtJetPt[iCut]  = new TH2F("MC_Pi0_Pt_JetPt","MC_Pi0_Pt_JetPt",nBinsQAPt, minQAPt, maxQAPt,200,0,200);
            fHistoMCPi0PtJetPt[iCut]->Sumw2();
            SetLogBinningXTH2(fHistoMCPi0PtJetPt[iCut]);
            fMCList[iCut]->Add(fHistoMCPi0PtJetPt[iCut]);
            fHistoMCEtaPtJetPt[iCut]  = new TH2F("MC_Eta_Pt_JetPt","MC_Eta_Pt_JetPt",nBinsQAPt, minQAPt, maxQAPt,200,0,200);
            fHistoMCEtaPtJetPt[iCut]->Sumw2();
            SetLogBinningXTH2(fHistoMCEtaPtJetPt[iCut]);
            fMCList[iCut]->Add(fHistoMCEtaPtJetPt[iCut]);
          }
        }

        fHistoMCPrimaryPtvsSource[iCut]   = new TH2F("MC_Primary_Pt_Source","MC_Primary_Pt_Source",nBinsPt, minPt, maxPt, 7, -0.5, 6.5);
        fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(1,"Pi+");
        fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(2,"Pi-");
        fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(3,"K+");
        fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(4,"K-");
        fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(5,"K0s");
        fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(6,"K0l");
        fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(7,"Lambda");
        fMCList[iCut]->Add(fHistoMCPrimaryPtvsSource[iCut]);

        fHistoMCSecPi0Source[iCut]      = new TH1F("MC_SecPi0_Source","MC_SecPi0_Source",5000,0.,5000);
        fMCList[iCut]->Add(fHistoMCSecPi0Source[iCut]);
        fHistoMCSecEtaSource[iCut]      = new TH1F("MC_SecEta_Source","MC_SecEta_Source",5000,0,5000);
        fMCList[iCut]->Add(fHistoMCSecEtaSource[iCut]);
        fHistoMCSecPi0PtvsSource[iCut]  = new TH2F("MC_SecPi0_Pt_Source","MC_SecPi0_Pt_Source",nBinsPt, minPt, maxPt,16,-0.5,15.5);
        fMCList[iCut]->Add(fHistoMCSecPi0PtvsSource[iCut]);
        fHistoMCSecPi0InAccPtvsSource[iCut]  = new TH2F("MC_SecPi0InAcc_Pt_Source","MC_SecPi0InAcc_Pt_Source",nBinsPt, minPt, maxPt, 16,-0.5,15.5);
        fMCList[iCut]->Add(fHistoMCSecPi0InAccPtvsSource[iCut]);
        fHistoMCSecEtaPt[iCut]          = new TH1F("MC_SecEta_Pt","MC_SecEta_Pt",nBinsPt, minPt, maxPt);
        fMCList[iCut]->Add(fHistoMCSecEtaPt[iCut]);
        if (fIsMC == 2) {
          fHistoMCPrimaryPtvsSource[iCut]->Sumw2();
          fHistoMCSecPi0PtvsSource[iCut]->Sumw2();
          fHistoMCSecPi0InAccPtvsSource[iCut]->Sumw2();
          fHistoMCSecEtaPt[iCut]->Sumw2();
        }

        
        if (fDoMesonQA > 0 && fDoMesonQA < 3){
          fHistoMCPi0PtY[iCut]        = new TH2F("MC_Pi0_Pt_Y","MC_Pi0_Pt_Y",nBinsQAPt, minQAPt, maxQAPt,150,-1.5,1.5);
          fHistoMCPi0PtY[iCut]->Sumw2();
          SetLogBinningXTH2(fHistoMCPi0PtY[iCut]);
          fMCList[iCut]->Add(fHistoMCPi0PtY[iCut]);
          fHistoMCEtaPtY[iCut]        = new TH2F("MC_Eta_Pt_Y","MC_Eta_Pt_Y",nBinsQAPt, minQAPt, maxQAPt,150,-1.5,1.5);
          fHistoMCEtaPtY[iCut]->Sumw2();
          SetLogBinningXTH2(fHistoMCEtaPtY[iCut]);
          fMCList[iCut]->Add(fHistoMCEtaPtY[iCut]);
          fHistoMCPi0PtAlpha[iCut]    = new TH2F("MC_Pi0_Pt_Alpha","MC_Pi0_Pt_Alpha",nBinsQAPt, minQAPt, maxQAPt,100,0,1);
          SetLogBinningXTH2(fHistoMCPi0PtAlpha[iCut]);
          fMCList[iCut]->Add(fHistoMCPi0PtAlpha[iCut]);
          fHistoMCEtaPtAlpha[iCut]    = new TH2F("MC_Eta_Pt_Alpha","MC_Eta_Pt_Alpha",nBinsQAPt, minQAPt, maxQAPt,100,0,1);
          SetLogBinningXTH2(fHistoMCEtaPtAlpha[iCut]);
          fMCList[iCut]->Add(fHistoMCEtaPtAlpha[iCut]);

          if (fIsMC == 2) {
            fHistoMCPi0PtAlpha[iCut]->Sumw2();
            fHistoMCEtaPtAlpha[iCut]->Sumw2();
          }
          
        }
        
      }
      fTrueList[iCut]                 = new TList();
      fTrueList[iCut]->SetName(Form("%s_%s_%s True histograms",cutstringEvent.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
      fTrueList[iCut]->SetOwner(kTRUE);
      fCutFolder[iCut]->Add(fTrueList[iCut]);
                  
      if(!fDoLightOutput){
        fHistoClusPhotonBGPt[iCut]          = new TH2F("ESD_TrueClusPhotonBG_Pt","ESD_TrueClusPhotonBG_Pt",nBinsClusterPt, minClusterPt, maxClusterPt,10,-0.5,9.5);
        fHistoClusPhotonBGPt[iCut]->GetYaxis()->SetBinLabel( 1,"Elec");
        fHistoClusPhotonBGPt[iCut]->GetYaxis()->SetBinLabel( 2,"Pion");
        fHistoClusPhotonBGPt[iCut]->GetYaxis()->SetBinLabel( 3,"Proton");
        fHistoClusPhotonBGPt[iCut]->GetYaxis()->SetBinLabel( 4,"Kaon");
        fHistoClusPhotonBGPt[iCut]->GetYaxis()->SetBinLabel( 5,"Neutron");
        fHistoClusPhotonBGPt[iCut]->GetYaxis()->SetBinLabel( 6,"K0s");
        fHistoClusPhotonBGPt[iCut]->GetYaxis()->SetBinLabel( 7,"Lambda");
        fHistoClusPhotonBGPt[iCut]->GetYaxis()->SetBinLabel( 8,"Muon");
        fHistoClusPhotonBGPt[iCut]->GetYaxis()->SetBinLabel( 9,"K0l");
        fHistoClusPhotonBGPt[iCut]->GetYaxis()->SetBinLabel( 10,"Rest");
        fTrueList[iCut]->Add(fHistoClusPhotonBGPt[iCut]);
        fHistoClusPhotonPlusConvBGPt[iCut]  = new TH2F("ESD_TrueClusPhotonPlusConvBG_Pt","ESD_TrueClusPhotonPlusConvBG_Pt",nBinsClusterPt, minClusterPt, maxClusterPt,10,-0.5,9.5);
        fHistoClusPhotonPlusConvBGPt[iCut]->GetYaxis()->SetBinLabel( 1,"Elec");
        fHistoClusPhotonPlusConvBGPt[iCut]->GetYaxis()->SetBinLabel( 2,"Pion");
        fHistoClusPhotonPlusConvBGPt[iCut]->GetYaxis()->SetBinLabel( 3,"Proton");
        fHistoClusPhotonPlusConvBGPt[iCut]->GetYaxis()->SetBinLabel( 4,"Kaon");
        fHistoClusPhotonPlusConvBGPt[iCut]->GetYaxis()->SetBinLabel( 5,"Neutron");
        fHistoClusPhotonPlusConvBGPt[iCut]->GetYaxis()->SetBinLabel( 6,"K0s");
        fHistoClusPhotonPlusConvBGPt[iCut]->GetYaxis()->SetBinLabel( 7,"Lambda");
        fHistoClusPhotonPlusConvBGPt[iCut]->GetYaxis()->SetBinLabel( 8,"Muon");
        fHistoClusPhotonPlusConvBGPt[iCut]->GetYaxis()->SetBinLabel( 9,"K0l");
        fHistoClusPhotonPlusConvBGPt[iCut]->GetYaxis()->SetBinLabel(10,"Rest");
        fTrueList[iCut]->Add(fHistoClusPhotonPlusConvBGPt[iCut]);
        
        if (fDoClusterQA > 0) {
          fHistoClustPhotonElectronBGPtM02[iCut]        = new TH2F("ESD_TrueClusPhotonElectronBG_Pt_M02","ESD_TrueClusPhotonElectronBG_Pt_M02",nBinsClusterPt, minClusterPt, maxClusterPt,100,0,1);
          fTrueList[iCut]->Add(fHistoClustPhotonElectronBGPtM02[iCut]);
          fHistoClustPhotonPionBGPtM02[iCut]            = new TH2F("ESD_TrueClusPhotonPionBG_Pt_M02","ESD_TrueClusPhotonPionBG_Pt_M02",nBinsClusterPt, minClusterPt, maxClusterPt,100,0,1);
          fTrueList[iCut]->Add(fHistoClustPhotonPionBGPtM02[iCut]);
          fHistoClustPhotonKaonBGPtM02[iCut]            = new TH2F("ESD_TrueClusPhotonKaonBG_Pt_M02","ESD_TrueClusPhotonKaonBG_Pt_M02",nBinsClusterPt, minClusterPt, maxClusterPt,100,0,1);
          fTrueList[iCut]->Add(fHistoClustPhotonKaonBGPtM02[iCut]);
          fHistoClustPhotonK0lBGPtM02[iCut]             = new TH2F("ESD_TrueClusPhotonK0lBG_Pt_M02","ESD_TrueClusPhotonK0lBG_Pt_M02",nBinsClusterPt, minClusterPt, maxClusterPt,100,0,1);
          fTrueList[iCut]->Add(fHistoClustPhotonK0lBGPtM02[iCut]);
          fHistoClustPhotonNeutronBGPtM02[iCut]         = new TH2F("ESD_TrueClusPhotonNeutronBG_Pt_M02","ESD_TrueClusPhotonNeutronBG_Pt_M02",nBinsClusterPt, minClusterPt, maxClusterPt,100,0,1);
          fTrueList[iCut]->Add(fHistoClustPhotonNeutronBGPtM02[iCut]);
          fHistoClustPhotonRestBGPtM02[iCut]            = new TH2F("ESD_TrueClusPhotonRestBG_Pt_M02","ESD_TrueClusPhotonRestBG_Pt_M02",nBinsClusterPt, minClusterPt, maxClusterPt,100,0,1);
          fTrueList[iCut]->Add(fHistoClustPhotonRestBGPtM02[iCut]);
          fHistoClustPhotonPlusConvElectronBGPtM02[iCut]= new TH2F("ESD_TrueClusPhotonPlusConvElectronBG_Pt_M02","ESD_TrueClusPhotonPlusConvElectronBG_Pt_M02",nBinsClusterPt, minClusterPt, maxClusterPt,100,0,1);
          fTrueList[iCut]->Add(fHistoClustPhotonPlusConvElectronBGPtM02[iCut]);
          fHistoClustPhotonPlusConvPionBGPtM02[iCut]    = new TH2F("ESD_TrueClusPhotonPlusConvPionBG_Pt_M02","ESD_TrueClusPhotonPlusConvPionBG_Pt_M02",nBinsClusterPt, minClusterPt, maxClusterPt,100,0,1);
          fTrueList[iCut]->Add(fHistoClustPhotonPlusConvPionBGPtM02[iCut]);
          fHistoClustPhotonPlusConvKaonBGPtM02[iCut]    = new TH2F("ESD_TrueClusPhotonPlusConvKaonBG_Pt_M02","ESD_TrueClusPhotonPlusConvKaonBG_Pt_M02",nBinsClusterPt, minClusterPt, maxClusterPt,100,0,1);
          fTrueList[iCut]->Add(fHistoClustPhotonPlusConvKaonBGPtM02[iCut]);
          fHistoClustPhotonPlusConvK0lBGPtM02[iCut]     = new TH2F("ESD_TrueClusPhotonPlusConvK0lBG_Pt_M02","ESD_TrueClusPhotonPlusConvK0lBG_Pt_M02",nBinsClusterPt, minClusterPt, maxClusterPt,100,0,1);
          fTrueList[iCut]->Add(fHistoClustPhotonPlusConvK0lBGPtM02[iCut]);
          fHistoClustPhotonPlusConvNeutronBGPtM02[iCut] = new TH2F("ESD_TrueClusPhotonPlusConvNeutronBG_Pt_M02","ESD_TrueClusPhotonPlusConvNeutronBG_Pt_M02",nBinsClusterPt, minClusterPt, maxClusterPt,100,0,1);
          fTrueList[iCut]->Add(fHistoClustPhotonPlusConvNeutronBGPtM02[iCut]);
          fHistoClustPhotonPlusConvRestBGPtM02[iCut]    = new TH2F("ESD_TrueClusPhotonPlusConvRestBG_Pt_M02","ESD_TrueClusPhotonPlusConvRestBG_Pt_M02",nBinsClusterPt, minClusterPt, maxClusterPt,100,0,1);
          fTrueList[iCut]->Add(fHistoClustPhotonPlusConvRestBGPtM02[iCut]);
        }
      }
    
      fHistoTrueClusGammaPt[iCut]                   = new TH1F("TrueClusGamma_Pt","ESD_TrueClusGamma_Pt",nBinsClusterPt, minClusterPt, maxClusterPt);
      fTrueList[iCut]->Add(fHistoTrueClusGammaPt[iCut]);
      if(!fDoLightOutput){
        if (fDoClusterQA > 0) {
          fHistoTrueClusGammaPtM02[iCut]              = new TH2F("TrueClusGamma_Pt_M02","TrueClusGamma_Pt_M02",nBinsClusterPt, minClusterPt, maxClusterPt,100,0,1);
          fTrueList[iCut]->Add(fHistoTrueClusGammaPtM02[iCut]);
        }
        fHistoTruePrimaryClusGammaPt[iCut]            = new TH1F("TruePrimaryClusGamma_Pt","ESD_TruePrimaryClusGamma_Pt",nBinsClusterPt, minClusterPt, maxClusterPt);
        fTrueList[iCut]->Add(fHistoTruePrimaryClusGammaPt[iCut]);
        fHistoTruePrimaryClusGammaESDPtMCPt[iCut]     = new TH2F("TruePrimaryClusGamma_Pt_MCPt","ESD_TruePrimaryClusGamma_MCPt",nBinsClusterPt, minClusterPt, maxClusterPt,nBinsClusterPt, minClusterPt, maxClusterPt);
        fTrueList[iCut]->Add(fHistoTruePrimaryClusGammaESDPtMCPt[iCut]);
        fHistoTruePrimaryClusConvGammaPt[iCut]        = new TH1F("TruePrimaryClusConvGamma_Pt","ESD_TruePrimaryClusConvGamma_Pt",nBinsClusterPt, minClusterPt, maxClusterPt);
        fTrueList[iCut]->Add(fHistoTruePrimaryClusConvGammaPt[iCut]);
        fHistoTruePrimaryClusConvGammaESDPtMCPt[iCut] = new TH2F("TruePrimaryClusConvGamma_Pt_MCPt","ESD_TruePrimaryClusConvGamma_MCPt",nBinsClusterPt, minClusterPt, maxClusterPt,nBinsClusterPt, minClusterPt, maxClusterPt);
        fTrueList[iCut]->Add(fHistoTruePrimaryClusConvGammaESDPtMCPt[iCut]);
        fHistoTrueSecondaryClusGammaPt[iCut]          = new TH2F("ESD_TrueSecondaryClusGamma_Pt","ESD_TrueSecondaryClusGamma_Pt",nBinsClusterPt, minClusterPt, maxClusterPt, 5, -0.5, 4.5);
        fHistoTrueSecondaryClusGammaPt[iCut]->GetYaxis()->SetBinLabel( 1,"K0s");
        fHistoTrueSecondaryClusGammaPt[iCut]->GetYaxis()->SetBinLabel( 2,"K0l");
        fHistoTrueSecondaryClusGammaPt[iCut]->GetYaxis()->SetBinLabel( 3,"Lambda");
        fHistoTrueSecondaryClusGammaPt[iCut]->GetYaxis()->SetBinLabel( 4,"Eta");
        fHistoTrueSecondaryClusGammaPt[iCut]->GetYaxis()->SetBinLabel( 5,"rest");
        fTrueList[iCut]->Add(fHistoTrueSecondaryClusGammaPt[iCut]);
        fHistoTrueSecondaryClusConvGammaPt[iCut]      = new TH2F("ESD_TrueSecondaryClusConvGamma_Pt","ESD_TrueSecondaryClusConvGamma_Pt",nBinsClusterPt, minClusterPt, maxClusterPt, 5, -0.5, 4.5);
        fHistoTrueSecondaryClusConvGammaPt[iCut]->GetYaxis()->SetBinLabel( 1,"K0s");
        fHistoTrueSecondaryClusConvGammaPt[iCut]->GetYaxis()->SetBinLabel( 2,"K0l");
        fHistoTrueSecondaryClusConvGammaPt[iCut]->GetYaxis()->SetBinLabel( 3,"Lambda");
        fHistoTrueSecondaryClusConvGammaPt[iCut]->GetYaxis()->SetBinLabel( 4,"Eta");
        fHistoTrueSecondaryClusConvGammaPt[iCut]->GetYaxis()->SetBinLabel( 5,"rest");
        fTrueList[iCut]->Add(fHistoTrueSecondaryClusConvGammaPt[iCut]);
        fHistoTrueSecondaryClusGammaMCPt[iCut]          = new TH2F("ESD_TrueSecondaryClusGamma_MCPt","ESD_TrueSecondaryClusGamma_MCPt",nBinsClusterPt, minClusterPt, maxClusterPt, 5, -0.5, 4.5);
        fHistoTrueSecondaryClusGammaMCPt[iCut]->GetYaxis()->SetBinLabel( 1,"K0s");
        fHistoTrueSecondaryClusGammaMCPt[iCut]->GetYaxis()->SetBinLabel( 2,"K0l");
        fHistoTrueSecondaryClusGammaMCPt[iCut]->GetYaxis()->SetBinLabel( 3,"Lambda");
        fHistoTrueSecondaryClusGammaMCPt[iCut]->GetYaxis()->SetBinLabel( 4,"Eta");
        fHistoTrueSecondaryClusGammaMCPt[iCut]->GetYaxis()->SetBinLabel( 5,"rest");
        fTrueList[iCut]->Add(fHistoTrueSecondaryClusGammaMCPt[iCut]);
        fHistoTrueSecondaryClusConvGammaMCPt[iCut]      = new TH2F("ESD_TrueSecondaryClusConvGamma_MCPt","ESD_TrueSecondaryClusConvGamma_MCPt",nBinsClusterPt, minClusterPt, maxClusterPt, 5, -0.5, 4.5);
        fHistoTrueSecondaryClusConvGammaMCPt[iCut]->GetYaxis()->SetBinLabel( 1,"K0s");
        fHistoTrueSecondaryClusConvGammaMCPt[iCut]->GetYaxis()->SetBinLabel( 2,"K0l");
        fHistoTrueSecondaryClusConvGammaMCPt[iCut]->GetYaxis()->SetBinLabel( 3,"Lambda");
        fHistoTrueSecondaryClusConvGammaMCPt[iCut]->GetYaxis()->SetBinLabel( 4,"Eta");
        fHistoTrueSecondaryClusConvGammaMCPt[iCut]->GetYaxis()->SetBinLabel( 5,"rest");
        fTrueList[iCut]->Add(fHistoTrueSecondaryClusConvGammaMCPt[iCut]);

        fHistoTrueSecondaryClusGammaFromXFromK0sMCPtESDPt[iCut] = new TH2F("ESD_TrueSecondaryClusGammaFromXFromK0s_MCPt_Pt", "ESD_TrueSecondaryClusGammaFromXFromK0s_MCPt_Pt",nBinsClusterPt, minClusterPt, maxClusterPt,nBinsClusterPt, minClusterPt, maxClusterPt);
        fTrueList[iCut]->Add(fHistoTrueSecondaryClusGammaFromXFromK0sMCPtESDPt[iCut]);
        fHistoTrueSecondaryClusConvGammaFromXFromK0sMCPtESDPt[iCut] = new TH2F("ESD_TrueSecondaryClusConvGammaFromXFromK0s_MCPt_Pt", "ESD_TrueSecondaryClusConvGammaFromXFromK0s_MCPt_Pt",nBinsClusterPt, minClusterPt, maxClusterPt,nBinsClusterPt, minClusterPt, maxClusterPt);
        fTrueList[iCut]->Add(fHistoTrueSecondaryClusConvGammaFromXFromK0sMCPtESDPt[iCut]);
        fHistoTrueSecondaryClusGammaFromXFromK0lMCPtESDPt[iCut] = new TH2F("ESD_TrueSecondaryClusGammaFromXFromK0l_MCPt_Pt", "ESD_TrueSecondaryClusGammaFromXFromK0l_MCPt_Pt",nBinsClusterPt, minClusterPt, maxClusterPt,nBinsClusterPt, minClusterPt, maxClusterPt);
        fTrueList[iCut]->Add(fHistoTrueSecondaryClusGammaFromXFromK0lMCPtESDPt[iCut]);
        fHistoTrueSecondaryClusConvGammaFromXFromK0lMCPtESDPt[iCut] = new TH2F("ESD_TrueSecondaryClusConvGammaFromXFromK0l_MCPt_Pt", "ESD_TrueSecondaryClusConvGammaFromXFromK0l_MCPt_Pt",nBinsClusterPt, minClusterPt, maxClusterPt,nBinsClusterPt, minClusterPt, maxClusterPt);
        fTrueList[iCut]->Add(fHistoTrueSecondaryClusConvGammaFromXFromK0lMCPtESDPt[iCut]);
        fHistoTrueSecondaryClusGammaFromXFromLambdaMCPtESDPt[iCut] = new TH2F("ESD_TrueSecondaryClusGammaFromXFromLambda_MCPt_Pt", "ESD_TrueSecondaryClusGammaFromXFromLambda_MCPt_Pt",nBinsClusterPt, minClusterPt, maxClusterPt,nBinsClusterPt, minClusterPt, maxClusterPt);
        fTrueList[iCut]->Add(fHistoTrueSecondaryClusGammaFromXFromLambdaMCPtESDPt[iCut]);
        fHistoTrueSecondaryClusConvGammaFromXFromLambdaMCPtESDPt[iCut] = new TH2F("ESD_TrueSecondaryClusConvGammaFromXFromLambda_MCPt_Pt", "ESD_TrueSecondaryClusConvGammaFromXFromLambda_MCPt_Pt",nBinsClusterPt, minClusterPt, maxClusterPt,nBinsClusterPt, minClusterPt, maxClusterPt);
        fTrueList[iCut]->Add(fHistoTrueSecondaryClusConvGammaFromXFromLambdaMCPtESDPt[iCut]);
        
        fHistoTrueNLabelsInClus[iCut]                 = new TH1F("TrueNLabelsInClus","TrueNLabelsInClus",100,-0.5,99.5);
        fTrueList[iCut]->Add(fHistoTrueNLabelsInClus[iCut]);
      }

      fHistoDoubleCountTrueClusterGammaPt[iCut]     = new TH2F("TrueDoubleCountClusterGamma_Pt","TrueDoubleCountClusterGamma_Pt",nBinsClusterPt, minClusterPt, maxClusterPt,2,0,2);
      fTrueList[iCut]->Add(fHistoDoubleCountTrueClusterGammaPt[iCut]);
      fHistoMultipleCountTrueClusterGamma[iCut]     = new TH1F("TrueMultipleCountClusterGamma","TrueMultipleCountClusterGamma",10,1,11);
      fTrueList[iCut]->Add(fHistoMultipleCountTrueClusterGamma[iCut]);

      if (fIsMC > 1){
        fHistoTrueClusGammaPt[iCut]->Sumw2();
        fHistoDoubleCountTrueClusterGammaPt[iCut]->Sumw2();
        fHistoMultipleCountTrueClusterGamma[iCut]->Sumw2();
        if(!fDoLightOutput){
          fHistoTrueNLabelsInClus[iCut]->Sumw2();
          if (fDoClusterQA > 0) fHistoTrueClusGammaPtM02[iCut]->Sumw2();
          fHistoTruePrimaryClusGammaPt[iCut]->Sumw2();
          fHistoTruePrimaryClusGammaESDPtMCPt[iCut]->Sumw2();
          fHistoTruePrimaryClusConvGammaPt[iCut]->Sumw2();
          fHistoTruePrimaryClusConvGammaESDPtMCPt[iCut]->Sumw2();
          fHistoTrueSecondaryClusGammaPt[iCut]->Sumw2();
          fHistoTrueSecondaryClusConvGammaPt[iCut]->Sumw2();
          fHistoTrueSecondaryClusGammaMCPt[iCut]->Sumw2();
          fHistoTrueSecondaryClusConvGammaMCPt[iCut]->Sumw2();
          fHistoTrueSecondaryClusGammaFromXFromK0sMCPtESDPt[iCut]->Sumw2();
          fHistoTrueSecondaryClusConvGammaFromXFromK0sMCPtESDPt[iCut]->Sumw2();
          fHistoTrueSecondaryClusGammaFromXFromK0lMCPtESDPt[iCut]->Sumw2();
          fHistoTrueSecondaryClusConvGammaFromXFromK0lMCPtESDPt[iCut]->Sumw2();
          fHistoTrueSecondaryClusGammaFromXFromLambdaMCPtESDPt[iCut]->Sumw2();
          fHistoTrueSecondaryClusConvGammaFromXFromLambdaMCPtESDPt[iCut]->Sumw2();
        }
      }

//      fHistoTruePi0NonLinearity[iCut] = new TH2F("TruePi0: E_truth / E_rec Vs E_rec","TruePi0: E_truth / E_rec Vs E_rec",700,0,35,200,0,2);
//      fTrueList[iCut]->Add(fHistoTruePi0NonLinearity[iCut]);
//      fHistoTrueEtaNonLinearity[iCut] = new TH2F("TrueEta: E_truth / E_rec Vs E_rec","TrueEta: E_truth / E_rec Vs E_rec",700,0,35,200,0,2);
//      fTrueList[iCut]->Add(fHistoTrueEtaNonLinearity[iCut]);

      if (fDoClusterQA > 0){
        fHistoTrueClusUnConvGammaPt[iCut]           = new TH1F("TrueClusUnConvGamma_Pt","TrueClusUnConvGamma_Pt",nBinsClusterPt, minClusterPt, maxClusterPt);
        fTrueList[iCut]->Add(fHistoTrueClusUnConvGammaPt[iCut]);
        fHistoTrueClusUnConvGammaMCPt[iCut]         = new TH1F("TrueClusUnConvGamma_MCPt","TrueClusUnConvGamma_MCPt",nBinsClusterPt, minClusterPt, maxClusterPt);
        fTrueList[iCut]->Add(fHistoTrueClusUnConvGammaMCPt[iCut]);
        if (!fDoLightOutput) {
          fHistoTrueClusUnConvGammaPtM02[iCut]      = new TH2F("TrueClusUnConvGamma_Pt_M02","TrueClusUnConvGamma_Pt_M02",nBinsClusterPt, minClusterPt, maxClusterPt,100,0,1);
          fTrueList[iCut]->Add(fHistoTrueClusUnConvGammaPtM02[iCut]);
        }
        fHistoTrueClusElectronPt[iCut]              = new TH1F("TrueClusElectron_Pt","TrueElectronGamma_Pt",nBinsClusterPt, minClusterPt, maxClusterPt);
        fTrueList[iCut]->Add(fHistoTrueClusElectronPt[iCut]);
        fHistoTrueClusConvGammaPt[iCut]             = new TH1F("TrueClusConvGamma_Pt","TrueClusConvGamma_Pt",nBinsClusterPt, minClusterPt, maxClusterPt);
        fTrueList[iCut]->Add(fHistoTrueClusConvGammaPt[iCut]);
        fHistoTrueClusConvGammaMCPt[iCut]           = new TH1F("TrueClusConvGamma_MCPt","TrueClusConvGamma_MCPt",nBinsClusterPt, minClusterPt, maxClusterPt);
        fTrueList[iCut]->Add(fHistoTrueClusConvGammaMCPt[iCut]);
        fHistoTrueClusConvGammaFullyPt[iCut]        = new TH1F("TrueClusConvGammaFullyContained_Pt","TrueClusConvGammaFullyContained_Pt",nBinsClusterPt, minClusterPt, maxClusterPt);
        fTrueList[iCut]->Add(fHistoTrueClusConvGammaFullyPt[iCut]);
        fHistoTrueClusMergedGammaPt[iCut]           = new TH1F("TrueClusMergedGamma_Pt","TrueClusMergedGamma_Pt",nBinsClusterPt, minClusterPt, maxClusterPt);
        fTrueList[iCut]->Add(fHistoTrueClusMergedGammaPt[iCut]);
        fHistoTrueClusMergedPartConvGammaPt[iCut]   = new TH1F("TrueClusMergedPartConvGamma_Pt","TrueClusMergedPartConvGamma_Pt",nBinsClusterPt, minClusterPt, maxClusterPt);
        fTrueList[iCut]->Add(fHistoTrueClusMergedPartConvGammaPt[iCut]);
        fHistoTrueClusDalitzPt[iCut]                = new TH1F("TrueClusDalitz_Pt","TrueClusDalitz_Pt",nBinsClusterPt, minClusterPt, maxClusterPt);
        fTrueList[iCut]->Add(fHistoTrueClusDalitzPt[iCut]);
        fHistoTrueClusDalitzMergedPt[iCut]          = new TH1F("TrueClusDalitzMerged_Pt","TrueClusDalitzMerged_Pt",nBinsClusterPt, minClusterPt, maxClusterPt);
        fTrueList[iCut]->Add(fHistoTrueClusDalitzMergedPt[iCut]);
        fHistoTrueClusPhotonFromElecMotherPt[iCut]  = new TH1F("TrueClusPhotonFromElecMother_Pt","TrueClusPhotonFromElecMother_Pt",nBinsClusterPt, minClusterPt, maxClusterPt);
        fTrueList[iCut]->Add(fHistoTrueClusPhotonFromElecMotherPt[iCut]);
        fHistoTrueClusShowerPt[iCut]                = new TH1F("TrueClusShower_Pt","TrueClusShower_Pt",nBinsClusterPt, minClusterPt, maxClusterPt);
        fTrueList[iCut]->Add(fHistoTrueClusShowerPt[iCut]);
        fHistoTrueClusSubLeadingPt[iCut]            = new TH1F("TrueClusSubleading_Pt","TrueClusSubleading_Pt",nBinsClusterPt, minClusterPt, maxClusterPt);
        fTrueList[iCut]->Add(fHistoTrueClusSubLeadingPt[iCut]);
        fHistoTrueClusNParticles[iCut]              = new TH1F("TrueClusNParticles","TrueClusNParticles",20,0,20);
        fTrueList[iCut]->Add(fHistoTrueClusNParticles[iCut]);
        fHistoTrueClusEMNonLeadingPt[iCut]          = new TH1F("TrueClusEMNonLeading_Pt","TrueClusEMNonLeading_Pt",nBinsClusterPt, minClusterPt, maxClusterPt);
        fTrueList[iCut]->Add(fHistoTrueClusEMNonLeadingPt[iCut]);

        if (fIsMC > 1){
            fHistoTrueClusUnConvGammaPt[iCut]->Sumw2();
            fHistoTrueClusUnConvGammaMCPt[iCut]->Sumw2();
            if (!fDoLightOutput)
              fHistoTrueClusUnConvGammaPtM02[iCut]->Sumw2();
            fHistoTrueClusElectronPt[iCut]->Sumw2();
            fHistoTrueClusConvGammaPt[iCut]->Sumw2();
            fHistoTrueClusConvGammaMCPt[iCut]->Sumw2();
            fHistoTrueClusConvGammaFullyPt[iCut]->Sumw2();
            fHistoTrueClusMergedGammaPt[iCut]->Sumw2();
            fHistoTrueClusMergedPartConvGammaPt[iCut]->Sumw2();
            fHistoTrueClusDalitzPt[iCut]->Sumw2();
            fHistoTrueClusDalitzMergedPt[iCut]->Sumw2();
            fHistoTrueClusPhotonFromElecMotherPt[iCut]->Sumw2();
            fHistoTrueClusShowerPt[iCut]->Sumw2();
            fHistoTrueClusSubLeadingPt[iCut]->Sumw2();
            fHistoTrueClusNParticles[iCut]->Sumw2();
            fHistoTrueClusEMNonLeadingPt[iCut]->Sumw2();
        }
      }

      if(fDoMesonAnalysis){
        fHistoTruePi0InvMassPt[iCut]                    = new TH2F("ESD_TruePi0_InvMass_Pt","ESD_TruePi0_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
        fTrueList[iCut]->Add(fHistoTruePi0InvMassPt[iCut]);
        fHistoTrueEtaInvMassPt[iCut]                    = new TH2F("ESD_TrueEta_InvMass_Pt","ESD_TrueEta_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
        fTrueList[iCut]->Add(fHistoTrueEtaInvMassPt[iCut]);

        fHistoDoubleCountTruePi0InvMassPt[iCut]         = new TH2F("ESD_TrueDoubleCountPi0_InvMass_Pt","ESD_TrueDoubleCountPi0_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
        fTrueList[iCut]->Add(fHistoDoubleCountTruePi0InvMassPt[iCut]);
        fHistoDoubleCountTrueEtaInvMassPt[iCut]         = new TH2F("ESD_TrueDoubleCountEta_InvMass_Pt","ESD_TrueDoubleCountEta_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
        fTrueList[iCut]->Add(fHistoDoubleCountTrueEtaInvMassPt[iCut]);

        fHistoTruePrimaryPi0InvMassPt[iCut]             = new TH2F("ESD_TruePrimaryPi0_InvMass_Pt", "ESD_TruePrimaryPi0_InvMass_Pt", 800,0,0.8,nBinsPt, minPt, maxPt);
        fHistoTruePrimaryPi0InvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTruePrimaryPi0InvMassPt[iCut]);
        fHistoTruePrimaryEtaInvMassPt[iCut]             = new TH2F("ESD_TruePrimaryEta_InvMass_Pt", "ESD_TruePrimaryEta_InvMass_Pt", 800,0,0.8,nBinsPt, minPt, maxPt);
        fHistoTruePrimaryEtaInvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTruePrimaryEtaInvMassPt[iCut]);

        fHistoTruePrimaryPi0W0WeightingInvMassPt[iCut]  = new TH2F("ESD_TruePrimaryPi0W0Weights_InvMass_Pt", "ESD_TruePrimaryPi0W0Weights_InvMass_Pt", 800,0,0.8,nBinsPt, minPt, maxPt);
        fHistoTruePrimaryPi0W0WeightingInvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTruePrimaryPi0W0WeightingInvMassPt[iCut]);
        fHistoTruePrimaryEtaW0WeightingInvMassPt[iCut]  = new TH2F("ESD_TruePrimaryEtaW0Weights_InvMass_Pt", "ESD_TruePrimaryEtaW0Weights_InvMass_Pt", 800,0,0.8,nBinsPt, minPt, maxPt);
        fHistoTruePrimaryEtaW0WeightingInvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTruePrimaryEtaW0WeightingInvMassPt[iCut]);

        fProfileTruePrimaryPi0WeightsInvMassPt[iCut]    = new TProfile2D("ESD_TruePrimaryPi0Weights_InvMass_Pt", "ESD_TruePrimaryPi0Weights_InvMass_Pt", 800,0,0.8,nBinsPt, minPt, maxPt);
        fProfileTruePrimaryPi0WeightsInvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fProfileTruePrimaryPi0WeightsInvMassPt[iCut]);
        fProfileTruePrimaryEtaWeightsInvMassPt[iCut]    = new TProfile2D("ESD_TruePrimaryEtaWeights_InvMass_Pt", "ESD_TruePrimaryEtaWeights_InvMass_Pt", 800,0,0.8,nBinsPt, minPt, maxPt);
        fProfileTruePrimaryEtaWeightsInvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fProfileTruePrimaryEtaWeightsInvMassPt[iCut]);

        fHistoTrueSecondaryPi0InvMassPt[iCut]           = new TH2F("ESD_TrueSecondaryPi0_InvMass_Pt", "ESD_TrueSecondaryPi0_InvMass_Pt", 800,0,0.8,nBinsPt, minPt, maxPt);
        fHistoTrueSecondaryPi0InvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTrueSecondaryPi0InvMassPt[iCut]);

        if(!fDoLightOutput){
          fHistoTruePi0InvMassPtAlpha[iCut]               = new TH2F("ESD_TruePi0_InvMass_vs_Pt_Alpha","ESD_TruePi0_InvMass_vs_Pt_Alpha",800,0,0.8,nBinsPt, minPt, maxPt);
          fHistoTruePi0InvMassPtAlpha[iCut]->Sumw2();
          fESDList[iCut]->Add(fHistoTruePi0InvMassPtAlpha[iCut]);
          fHistoTruePi0PureGammaInvMassPtAlpha[iCut]      = new TH2F("ESD_TruePi0PureGamma_InvMass_vs_Pt_Alpha","ESD_TruePi0PureGamma_InvMass_vs_Pt_Alpha",800,0,0.8,nBinsPt, minPt, maxPt);
          fHistoTruePi0PureGammaInvMassPtAlpha[iCut]->Sumw2();
          fESDList[iCut]->Add(fHistoTruePi0PureGammaInvMassPtAlpha[iCut]);
        }

        fHistoTrueSecondaryPi0FromK0sInvMassPt[iCut]    = new TH2F("ESD_TrueSecondaryPi0FromK0s_InvMass_Pt","ESD_TrueSecondaryPi0FromK0s_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
        fHistoTrueSecondaryPi0FromK0sInvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTrueSecondaryPi0FromK0sInvMassPt[iCut]);
        fHistoTrueSecondaryPi0FromK0lInvMassPt[iCut]    = new TH2F("ESD_TrueSecondaryPi0FromK0l_InvMass_Pt","ESD_TrueSecondaryPi0FromK0l_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
        fHistoTrueSecondaryPi0FromK0lInvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTrueSecondaryPi0FromK0lInvMassPt[iCut]);
        fHistoTrueSecondaryPi0FromEtaInvMassPt[iCut]    = new TH2F("ESD_TrueSecondaryPi0FromEta_InvMass_Pt","ESD_TrueSecondaryPi0FromEta_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
        fTrueList[iCut]->Add(fHistoTrueSecondaryPi0FromEtaInvMassPt[iCut]);
        fHistoTrueSecondaryPi0FromLambdaInvMassPt[iCut] = new TH2F("ESD_TrueSecondaryPi0FromLambda_InvMass_Pt","ESD_TrueSecondaryPi0FromLambda_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
        fTrueList[iCut]->Add(fHistoTrueSecondaryPi0FromLambdaInvMassPt[iCut]);

        if (fIsMC > 1){
          fHistoTruePi0InvMassPt[iCut]->Sumw2();
          fHistoTrueEtaInvMassPt[iCut]->Sumw2();
          fHistoDoubleCountTruePi0InvMassPt[iCut]->Sumw2();
          fHistoDoubleCountTrueEtaInvMassPt[iCut]->Sumw2();
          fHistoTruePrimaryPi0W0WeightingInvMassPt[iCut]->Sumw2();
          fHistoTruePrimaryEtaW0WeightingInvMassPt[iCut]->Sumw2();
          fHistoTrueSecondaryPi0FromEtaInvMassPt[iCut]->Sumw2();
          fHistoTrueSecondaryPi0FromLambdaInvMassPt[iCut]->Sumw2();
        }


        if (fDoMesonQA > 0 && fDoMesonQA < 3){
          if (fIsMC < 2){

            fHistoTruePi0CaloPhotonInvMassPt[iCut]      = new TH2F("ESD_TruePi0CaloPhoton_InvMass_Pt","ESD_TruePi0CaloPhoton_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
            fTrueList[iCut]->Add(fHistoTruePi0CaloPhotonInvMassPt[iCut]);
            fHistoTrueEtaCaloPhotonInvMassPt[iCut]      = new TH2F("ESD_TrueEtaCaloPhoton_InvMass_Pt","ESD_TrueEtaCaloPhoton_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
            fTrueList[iCut]->Add(fHistoTrueEtaCaloPhotonInvMassPt[iCut]);

            fHistoTruePi0CaloMixedPhotonConvPhotonInvMassPt[iCut]   = new TH2F("ESD_TruePi0CaloMixedPhotonConvertedPhoton_InvMass_Pt","ESD_TruePi0CaloMixedPhotonConvertedPhoton_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
            fTrueList[iCut]->Add(fHistoTruePi0CaloMixedPhotonConvPhotonInvMassPt[iCut]);
            fHistoTrueEtaCaloMixedPhotonConvPhotonInvMassPt[iCut]   = new TH2F("ESD_TrueEtaCaloMixedPhotonConvertedPhoton_InvMass_Pt","ESD_TrueEtaCaloMixedPhotonConvertedPhoton_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
            fTrueList[iCut]->Add(fHistoTrueEtaCaloMixedPhotonConvPhotonInvMassPt[iCut]);

            fHistoTruePi0CaloConvertedPhotonInvMassPt[iCut]         = new TH2F("ESD_TruePi0CaloConvertedPhoton_InvMass_Pt","ESD_TruePi0CaloConvertedPhoton_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
            fTrueList[iCut]->Add(fHistoTruePi0CaloConvertedPhotonInvMassPt[iCut]);
            fHistoTrueEtaCaloConvertedPhotonInvMassPt[iCut]         = new TH2F("ESD_TrueEtaCaloConvertedPhoton_InvMass_Pt","ESD_TrueEtaCaloConvertedPhoton_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
            fTrueList[iCut]->Add(fHistoTrueEtaCaloConvertedPhotonInvMassPt[iCut]);

            fHistoTruePi0CaloElectronInvMassPt[iCut]    = new TH2F("ESD_TruePi0CaloElectron_InvMass_Pt","ESD_TruePi0CaloElectron_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
            fTrueList[iCut]->Add(fHistoTruePi0CaloElectronInvMassPt[iCut]);
            fHistoTrueEtaCaloElectronInvMassPt[iCut]    = new TH2F("ESD_TrueEtaCaloElectron_InvMass_Pt","ESD_TrueEtaCaloElectron_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
            fTrueList[iCut]->Add(fHistoTrueEtaCaloElectronInvMassPt[iCut]);

            fHistoTruePi0CaloMergedClusterInvMassPt[iCut]           = new TH2F("ESD_TruePi0CaloMergedCluster_InvMass_Pt","ESD_TruePi0CaloMergedCluster_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
            fTrueList[iCut]->Add(fHistoTruePi0CaloMergedClusterInvMassPt[iCut]);
            fHistoTrueEtaCaloMergedClusterInvMassPt[iCut]           = new TH2F("ESD_TrueEtaCaloMergedCluster_InvMass_Pt","ESD_TrueEtaCaloMergedCluster_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
            fTrueList[iCut]->Add(fHistoTrueEtaCaloMergedClusterInvMassPt[iCut]);

            fHistoTruePi0CaloMergedClusterPartConvInvMassPt[iCut]   = new TH2F("ESD_TruePi0CaloMergedClusterPartConv_InvMass_Pt","ESD_TruePi0CaloMergedClusterPartConv_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
            fTrueList[iCut]->Add(fHistoTruePi0CaloMergedClusterPartConvInvMassPt[iCut]);
            fHistoTrueEtaCaloMergedClusterPartConvInvMassPt[iCut]   = new TH2F("ESD_TrueEtaCaloMergedClusterPartConv_InvMass_Pt","ESD_TrueEtaCaloMergedClusterPartConv_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
            fTrueList[iCut]->Add(fHistoTrueEtaCaloMergedClusterPartConvInvMassPt[iCut]);

            fHistoTruePi0NonMergedElectronPhotonInvMassPt[iCut]     = new TH2F("ESD_TruePi0NonMergedElectronPhoton_InvMass_Pt","ESD_TruePi0NonMergedElectronPhoton_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
            fTrueList[iCut]->Add(fHistoTruePi0NonMergedElectronPhotonInvMassPt[iCut]);
            fHistoTruePi0NonMergedElectronMergedPhotonInvMassPt[iCut] = new TH2F("ESD_TruePi0NonMergedElectronMergedPhoton_InvMass_Pt","ESD_TruePi0NonMergedElectronMergedPhoton_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
            fTrueList[iCut]->Add(fHistoTruePi0NonMergedElectronMergedPhotonInvMassPt[iCut]);
            
            fHistoTruePrimaryPi0MCPtResolPt[iCut]       = new TH2F("ESD_TruePrimaryPi0_MCPt_ResolPt","ESD_TruePrimaryPi0_ResolPt_MCPt",500,0.03,35,1000,-1.,1.);
            fHistoTruePrimaryPi0MCPtResolPt[iCut]->Sumw2();
            SetLogBinningXTH2(fHistoTruePrimaryPi0MCPtResolPt[iCut]);
            fTrueList[iCut]->Add(fHistoTruePrimaryPi0MCPtResolPt[iCut]);
            fHistoTruePrimaryEtaMCPtResolPt[iCut]       = new TH2F("ESD_TruePrimaryEta_MCPt_ResolPt","ESD_TruePrimaryEta_ResolPt_MCPt",500,0.03,35,1000,-1.,1.);
            fHistoTruePrimaryEtaMCPtResolPt[iCut]->Sumw2();
            SetLogBinningXTH2(fHistoTruePrimaryEtaMCPtResolPt[iCut]);
            fTrueList[iCut]->Add(fHistoTruePrimaryEtaMCPtResolPt[iCut]);
            fHistoTrueK0sWithPi0DaughterMCPt[iCut]      = new TH1F("ESD_TrueK0sWithPi0Daughter_MCPt","ESD_TrueK0sWithPi0Daughter_MCPt",nBinsPt, minPt, maxPt);
            fTrueList[iCut]->Add(fHistoTrueK0sWithPi0DaughterMCPt[iCut]);
            fHistoTrueK0lWithPi0DaughterMCPt[iCut]      = new TH1F("ESD_TrueK0lWithPi0Daughter_MCPt","ESD_TrueK0lWithPi0Daughter_MCPt",nBinsPt, minPt, maxPt);
            fTrueList[iCut]->Add(fHistoTrueK0lWithPi0DaughterMCPt[iCut]);
            fHistoTrueEtaWithPi0DaughterMCPt[iCut]      = new TH1F("ESD_TrueEtaWithPi0Daughter_MCPt","ESD_TrueEtaWithPi0Daughter_MCPt",nBinsPt, minPt, maxPt);
            fTrueList[iCut]->Add(fHistoTrueEtaWithPi0DaughterMCPt[iCut]);
            fHistoTrueLambdaWithPi0DaughterMCPt[iCut]   = new TH1F("ESD_TrueLambdaWithPi0Daughter_MCPt","ESD_TrueLambdaWithPi0Daughter_MCPt",nBinsPt, minPt, maxPt);
            fTrueList[iCut]->Add(fHistoTrueLambdaWithPi0DaughterMCPt[iCut]);
          }
          
          fHistoTruePi0PtY[iCut]          = new TH2F("ESD_TruePi0_Pt_Y","ESD_TruePi0_Pt_Y",nBinsQAPt, minQAPt, maxQAPt,150,-1.5,1.5);
          SetLogBinningXTH2(fHistoTruePi0PtY[iCut]);
          fTrueList[iCut]->Add(fHistoTruePi0PtY[iCut]);
          fHistoTrueEtaPtY[iCut]          = new TH2F("ESD_TrueEta_Pt_Y","ESD_TrueEta_Pt_Y",nBinsQAPt, minQAPt, maxQAPt,150,-1.5,1.5);
          SetLogBinningXTH2(fHistoTrueEtaPtY[iCut]);
          fTrueList[iCut]->Add(fHistoTrueEtaPtY[iCut]);
          fHistoTruePi0PtAlpha[iCut]      = new TH2F("ESD_TruePi0_Pt_Alpha","ESD_TruePi0_Pt_Alpha",nBinsQAPt, minQAPt, maxQAPt,100,0,1);
          SetLogBinningXTH2(fHistoTruePi0PtAlpha[iCut]);
          fTrueList[iCut]->Add(fHistoTruePi0PtAlpha[iCut]);
          fHistoTrueEtaPtAlpha[iCut]      = new TH2F("ESD_TrueEta_Pt_Alpha","ESD_TrueEta_Pt_Alpha",nBinsQAPt, minQAPt, maxQAPt,100,0,1);
          SetLogBinningXTH2(fHistoTrueEtaPtAlpha[iCut]);
          fTrueList[iCut]->Add(fHistoTrueEtaPtAlpha[iCut]);
          
          fHistoTruePi0PtOpenAngle[iCut]  = new TH2F("ESD_TruePi0_Pt_OpenAngle","ESD_TruePi0_Pt_OpenAngle",nBinsQAPt, minQAPt, maxQAPt,100,0,0.5);
          SetLogBinningXTH2(fHistoTruePi0PtOpenAngle[iCut]);
          fTrueList[iCut]->Add(fHistoTruePi0PtOpenAngle[iCut]);
          fHistoTrueEtaPtOpenAngle[iCut]  = new TH2F("ESD_TrueEta_Pt_OpenAngle","ESD_TrueEta_Pt_OpenAngle",nBinsQAPt, minQAPt, maxQAPt,180,0,1.8);
          SetLogBinningXTH2(fHistoTrueEtaPtOpenAngle[iCut]);
          fTrueList[iCut]->Add(fHistoTrueEtaPtOpenAngle[iCut]);

          fHistoTrueBckGGInvMassPt[iCut]              = new TH2F("ESD_TrueBckGG_InvMass_Pt","ESD_TrueBckGG_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
          fTrueList[iCut]->Add(fHistoTrueBckGGInvMassPt[iCut]);
          fHistoTrueBckFullMesonContainedInOneClusterInvMassPt[iCut] = new TH2F("ESD_TrueBckFullMesonContained_InvMass_Pt","ESD_TrueBckFullMesonContained_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
          fTrueList[iCut]->Add(fHistoTrueBckFullMesonContainedInOneClusterInvMassPt[iCut]);
          fHistoTrueBckAsymEClustersInvMassPt[iCut]   = new TH2F("ESD_TrueBckAsymEClus_InvMass_Pt","ESD_TrueBckAsymEClus_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
          fTrueList[iCut]->Add(fHistoTrueBckAsymEClustersInvMassPt[iCut]);
          fHistoTrueBckContInvMassPt[iCut]            = new TH2F("ESD_TrueBckCont_InvMass_Pt","ESD_TrueBckCont_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
          fTrueList[iCut]->Add(fHistoTrueBckContInvMassPt[iCut]);
          
          if (fIsMC > 1){
            fHistoTruePi0PtY[iCut]->Sumw2();
            fHistoTrueEtaPtY[iCut]->Sumw2();
            fHistoTruePi0PtAlpha[iCut]->Sumw2();
            fHistoTrueEtaPtAlpha[iCut]->Sumw2();
            fHistoTruePi0PtOpenAngle[iCut]->Sumw2();
            fHistoTrueEtaPtOpenAngle[iCut]->Sumw2();
            fHistoTrueBckGGInvMassPt[iCut]->Sumw2();
            fHistoTrueBckFullMesonContainedInOneClusterInvMassPt[iCut]->Sumw2();
            fHistoTrueBckAsymEClustersInvMassPt[iCut]->Sumw2();
            fHistoTrueBckContInvMassPt[iCut]->Sumw2();
          }
          
        }
        
        if (fDoMesonQA == 2 && fIsMC < 2){
          fHistoTruePi0Category1[iCut]    = new TH2F("ESD_TruePi0Category1_InvMass_Pt","ESD_TruePi0Category1_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
          fTrueList[iCut]->Add(fHistoTruePi0Category1[iCut]);
          fHistoTrueEtaCategory1[iCut]    = new TH2F("ESD_TrueEtaCategory1_InvMass_Pt","ESD_TrueEtaCategory1_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
          fTrueList[iCut]->Add(fHistoTrueEtaCategory1[iCut]);
          fHistoTruePi0Category2[iCut]    = new TH2F("ESD_TruePi0Category2_InvMass_Pt","ESD_TruePi0Category2_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
          fTrueList[iCut]->Add(fHistoTruePi0Category2[iCut]);
          fHistoTrueEtaCategory2[iCut]    = new TH2F("ESD_TrueEtaCategory2_InvMass_Pt","ESD_TrueEtaCategory2_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
          fTrueList[iCut]->Add(fHistoTrueEtaCategory2[iCut]);
          fHistoTruePi0Category3[iCut]    = new TH2F("ESD_TruePi0Category3_InvMass_Pt","ESD_TruePi0Category3_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
          fTrueList[iCut]->Add(fHistoTruePi0Category3[iCut]);
          fHistoTrueEtaCategory3[iCut]    = new TH2F("ESD_TrueEtaCategory3_InvMass_Pt","ESD_TrueEtaCategory3_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
          fTrueList[iCut]->Add(fHistoTrueEtaCategory3[iCut]);
          fHistoTruePi0Category4_6[iCut]  = new TH2F("ESD_TruePi0Category4_6_InvMass_Pt","ESD_TruePi0Category4_6_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
          fTrueList[iCut]->Add(fHistoTruePi0Category4_6[iCut]);
          fHistoTrueEtaCategory4_6[iCut]  = new TH2F("ESD_TrueEtaCategory4_6_InvMass_Pt","ESD_TrueEtaCategory4_6_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
          fTrueList[iCut]->Add(fHistoTrueEtaCategory4_6[iCut]);
          fHistoTruePi0Category5[iCut]    = new TH2F("ESD_TruePi0Category5_InvMass_Pt","ESD_TruePi0Category5_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
          fTrueList[iCut]->Add(fHistoTruePi0Category5[iCut]);
          fHistoTrueEtaCategory5[iCut]    = new TH2F("ESD_TrueEtaCategory5_InvMass_Pt","ESD_TrueEtaCategory5_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
          fTrueList[iCut]->Add(fHistoTrueEtaCategory5[iCut]);
          fHistoTruePi0Category7[iCut]    = new TH2F("ESD_TruePi0Category7_InvMass_Pt","ESD_TruePi0Category7_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
          fTrueList[iCut]->Add(fHistoTruePi0Category7[iCut]);
          fHistoTrueEtaCategory7[iCut]    = new TH2F("ESD_TrueEtaCategory7_InvMass_Pt","ESD_TrueEtaCategory7_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
          fTrueList[iCut]->Add(fHistoTrueEtaCategory7[iCut]);
          fHistoTruePi0Category8[iCut]    = new TH2F("ESD_TruePi0Category8_InvMass_Pt","ESD_TruePi0Category8_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
          fTrueList[iCut]->Add(fHistoTruePi0Category8[iCut]);
          fHistoTrueEtaCategory8[iCut]    = new TH2F("ESD_TrueEtaCategory8_InvMass_Pt","ESD_TrueEtaCategory8_InvMass_Pt",800,0,0.8,nBinsPt, minPt, maxPt);
          fTrueList[iCut]->Add(fHistoTrueEtaCategory8[iCut]);
        }
        
        if (fDoMesonQA == 3){
          fTreeList[iCut] = new TList();
          fTreeList[iCut]->SetName(Form("%s_%s_%s True ClusterComb tree",cutstringEvent.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
          fTreeList[iCut]->SetOwner(kTRUE);
          fCutFolder[iCut]->Add(fTreeList[iCut]);
            
          tTrueInvMassROpenABPtFlag[iCut] = new TTree("True_InvMass_R_OpenA_OpenB_Pt_Flag","True_InvMass_R_OpenA_OpenB_Pt_Flag");   
          tTrueInvMassROpenABPtFlag[iCut]->Branch("InvMass",&fInvMass,"fInvMass/F");
          tTrueInvMassROpenABPtFlag[iCut]->Branch("RConv",&fRconv,"fRconv/F");
          tTrueInvMassROpenABPtFlag[iCut]->Branch("OpenAngleRPrimVtx",&fOpenRPrim,"fOpenRPrim/F");
          tTrueInvMassROpenABPtFlag[iCut]->Branch("InvMassRTOF",&fInvMassRTOF,"fInvMassRTOF/F");
          tTrueInvMassROpenABPtFlag[iCut]->Branch("Pt",&fPt,"fPt/F");
//           tTrueInvMassROpenABPtFlag[iCut]->Branch("Weight",&fWeightJetJetMC,"fWeightJetJetMC/F");
          tTrueInvMassROpenABPtFlag[iCut]->Branch("cat",&iFlag,"iFlag/b");
          fTreeList[iCut]->Add(tTrueInvMassROpenABPtFlag[iCut]);
        }
      }
    }
  }
    
  fVectorDoubleCountTruePi0s.clear();
  fVectorDoubleCountTrueEtas.clear();
  fVectorDoubleCountTrueClusterGammas.clear();

  fMapMultipleCountTrueClusterGammas.clear();
  
  if(fV0Reader)
    if((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())
      if(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms())
        fOutputContainer->Add(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms());
  if(fV0Reader)
    if((AliConvEventCuts*)fV0Reader->GetEventCuts())
      if(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms())
        fOutputContainer->Add(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms());

  for(Int_t iMatcherTask = 0; iMatcherTask < 3; iMatcherTask++){
    AliCaloTrackMatcher* temp = (AliCaloTrackMatcher*) (AliAnalysisManager::GetAnalysisManager()->GetTask(Form("CaloTrackMatcher_%i",iMatcherTask)));
    if(temp) fOutputContainer->Add(temp->GetCaloTrackMatcherHistograms());
  }
      
  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if(!((AliConvEventCuts*)fEventCutArray->At(iCut))) continue;
    if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms()){
      fCutFolder[iCut]->Add(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms());
    }
    if(!((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))) continue;
    if(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutHistograms()){
      fCutFolder[iCut]->Add(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutHistograms());
    }
    if(fSetPlotHistsExtQA){
      if(!((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))) continue;
      if(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetExtQAHistograms()){
        fCutFolder[iCut]->Add(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetExtQAHistograms());
      }
    }
    if(fDoMesonAnalysis){
      if(!((AliConversionMesonCuts*)fMesonCutArray->At(iCut))) continue;
      if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms()){
        fCutFolder[iCut]->Add(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms());
      }
    }
  }
  
  if (fIsMC > 0 || fDoClusterQA > 0){
    tBrokenFiles = new TTree("BrokenFiles","BrokenFiles");   
    tBrokenFiles->Branch("fileName",&fFileNameBroken);
    tBrokenFiles->Branch("closeHighPtClusters",&fCloseHighPtClusters);
    fOutputContainer->Add(tBrokenFiles);
  }
  
  if(fLocalDebugFlag > 0){
    fstream fOutputLocalDebug;
    fOutputLocalDebug.open("debugOutput.txt",ios::out);
    fOutputLocalDebug.close();
  }

  PostData(1, fOutputContainer);
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskGammaCalo::Notify()
{
  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetPeriodEnum() == AliConvEventCuts::kNoPeriod && ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum() != AliConvEventCuts::kNoPeriod){        
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->SetPeriodEnumExplicit(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum());
    } else if (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetPeriodEnum() == AliConvEventCuts::kNoPeriod ){
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->SetPeriodEnum(fV0Reader->GetPeriodName());
    }  
    if(fIsHeavyIon==2) {
      if(!((AliConvEventCuts*)fEventCutArray->At(iCut))->GetDoEtaShift()){
        fProfileEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
        continue; // No Eta Shift requested, continue
      }
      if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift() == 0.0){ // Eta Shift requested but not set, get shift automatically
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCorrectEtaShiftFromPeriod();
        fProfileEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once
        continue;
      }
      else{
        printf(" Gamma Conversion Task %s :: Eta Shift Manually Set to %f \n\n",
            (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber()).Data(),((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift());
        fProfileEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once
      }
    }
  }
  
  return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskGammaCalo::UserExec(Option_t *)
{
  //
  // Called for each event
  //
  fInputEvent = InputEvent();
  
  if(fIsMC> 0) fMCEvent = MCEvent();

  Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();
  if(fInputEvent->IsIncompleteDAQ()==kTRUE) eventQuality = 2;  // incomplete event
  if(eventQuality == 2 || eventQuality == 3){// Event Not Accepted due to MC event missing or wrong trigger for V0ReaderV1 or because it is incomplete
    // write out name of broken file for first event
    if (fIsMC > 0){
      if (fInputEvent->IsA()==AliESDEvent::Class()){
        if (((AliESDEvent*)fInputEvent)->GetEventNumberInFile() == 0){  
          fFileNameBroken = new TObjString(Form("%s",((TString)fV0Reader->GetCurrentFileName()).Data()));
          if (tBrokenFiles) tBrokenFiles->Fill();
          delete fFileNameBroken;
        } 
      }  
    }  

    for(Int_t iCut = 0; iCut<fnCuts; iCut++){
      fHistoNEvents[iCut]->Fill(eventQuality);
      if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventQuality);
    }
    return;
  }  
  
  // ------------------- BeginEvent ----------------------------
  
  AliEventplane *EventPlane = fInputEvent->GetEventplane();
  if(fIsHeavyIon ==1)fEventPlaneAngle = EventPlane->GetEventplane("V0",fInputEvent,2);
  else fEventPlaneAngle=0.0;
  
  for(Int_t iCut = 0; iCut<fnCuts; iCut++){
    
    fiCut = iCut;
    
    Bool_t isRunningEMCALrelAna = kFALSE;
    if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 1) isRunningEMCALrelAna = kTRUE;
    
    Int_t eventNotAccepted = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon, isRunningEMCALrelAna);
    
    if(fIsMC==2){
      Float_t xsection      = -1.; 
      Float_t ntrials       = -1.;
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetXSectionAndNTrials(fMCEvent,xsection,ntrials);
      if((xsection==-1.) || (ntrials==-1.)) AliFatal(Form("ERROR: GetXSectionAndNTrials returned invalid xsection/ntrials, periodName from V0Reader: '%s'",fV0Reader->GetPeriodName().Data()));
      fProfileJetJetXSection[iCut]->Fill(0.,xsection);
      fHistoJetJetNTrials[iCut]->Fill("#sum{NTrials}",ntrials);
    }

    if (fIsMC > 0){
      fWeightJetJetMC       = 1;
      Bool_t isMCJet        = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsJetJetMCEventAccepted( fMCEvent, fWeightJetJetMC );
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
    if(eventNotAccepted){
    // cout << "event rejected due to wrong trigger: " <<eventNotAccepted << endl;
      fHistoNEvents[iCut]->Fill(eventNotAccepted, fWeightJetJetMC); // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1
      if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventNotAccepted);
      if (eventNotAccepted==3 && fIsMC>0){
        triggered = kFALSE;
      } else {  
        continue;
      }
    }

    if(eventQuality != 0){// Event Not Accepted
      //cout << "event rejected due to: " <<eventQuality << endl;
      fHistoNEvents[iCut]->Fill(eventQuality, fWeightJetJetMC);
      if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventQuality); // Should be 0 here
      continue;
    }
    if (triggered == kTRUE) {
      fHistoNEvents[iCut]->Fill(eventQuality, fWeightJetJetMC); // Should be 0 here
      if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventQuality); // Should be 0 here

      fHistoNGoodESDTracks[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks(), fWeightJetJetMC);
      fHistoVertexZ[iCut]->Fill(fInputEvent->GetPrimaryVertex()->GetZ(), fWeightJetJetMC);
      if(!fDoLightOutput){
        fHistoSPDClusterTrackletBackground[iCut]->Fill(fInputEvent->GetMultiplicity()->GetNumberOfTracklets(),(fInputEvent->GetNumberOfITSClusters(0)+fInputEvent->GetNumberOfITSClusters(1)), fWeightJetJetMC);
        if(((AliConvEventCuts*)fEventCutArray->At(iCut))->IsHeavyIon() == 2)  fHistoNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A(), fWeightJetJetMC);
          else fHistoNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A()+fInputEvent->GetVZEROData()->GetMTotV0C(), fWeightJetJetMC);
      }
    }
    if(fIsMC> 0){
      // Process MC Particle
      if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetSignalRejection() != 0){
        if(fInputEvent->IsA()==AliESDEvent::Class()){
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetNotRejectedParticles(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetSignalRejection(),
                                          ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetAcceptedHeader(),
                                          fMCEvent);
        }
        else if(fInputEvent->IsA()==AliAODEvent::Class()){
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetNotRejectedParticles(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetSignalRejection(),
                                          ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetAcceptedHeader(),
                                          fInputEvent);
        }

        if(!fDoLightOutput){
          if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetAcceptedHeader()){
            for(Int_t i = 0;i<(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetAcceptedHeader())->GetEntries();i++){
              TString nameBin= fHistoMCHeaders[iCut]->GetXaxis()->GetBinLabel(i+1);
              if (nameBin.CompareTo("")== 0){
                TString nameHeader = ((TObjString*)((TList*)((AliConvEventCuts*)fEventCutArray->At(iCut))
                                  ->GetAcceptedHeader())->At(i))->GetString();
                fHistoMCHeaders[iCut]->GetXaxis()->SetBinLabel(i+1,nameHeader.Data());
              }
            }
          }
        }
      }
    }
    if(fIsMC> 0){
    if(fInputEvent->IsA()==AliESDEvent::Class())
      ProcessMCParticles();
    if(fInputEvent->IsA()==AliAODEvent::Class())
      ProcessAODMCParticles();
    }
    
    if (triggered==kFALSE) continue;
    
    // it is in the loop to have the same conversion cut string (used also for MC stuff that should be same for V0 and Cluster)
    ProcessClusters();            // process calo clusters

    fHistoNGammaCandidates[iCut]->Fill(fClusterCandidates->GetEntries(), fWeightJetJetMC);
    if(!fDoLightOutput) fHistoNGoodESDTracksVsNGammaCandidates[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks(),fClusterCandidates->GetEntries(), fWeightJetJetMC);
    if(fDoMesonAnalysis){ // Meson Analysis
      
      CalculatePi0Candidates(); // Combine Gammas from conversion and from calo
      if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->DoBGCalculation()){
        if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->BackgroundHandlerType() == 0){
          
          CalculateBackground(); // Combinatorial Background
          UpdateEventByEventData(); // Store Event for mixed Events
        }
        
      }
      fVectorDoubleCountTruePi0s.clear();
      fVectorDoubleCountTrueEtas.clear();
    }

    if(fIsMC> 0){
      fVectorDoubleCountTrueClusterGammas.clear();
      FillMultipleCountHistoAndClear(fMapMultipleCountTrueClusterGammas,fHistoMultipleCountTrueClusterGamma[iCut]);
    }

    fClusterCandidates->Clear(); // delete cluster candidates
  }
  
  PostData(1, fOutputContainer);
}

//________________________________________________________________________
void AliAnalysisTaskGammaCalo::ProcessClusters()
{
  
  Int_t nclus = 0;
  nclus = fInputEvent->GetNumberOfCaloClusters();
  
//   cout << nclus << endl;
  
  if(nclus == 0)  return;
  
  // plotting histograms on cell/tower level, only if extendedMatchAndQA > 1
  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FillHistogramsExtendedQA(fInputEvent,fIsMC);

  // match tracks to clusters
  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchTracksToClusters(fInputEvent,fWeightJetJetMC);

  // vertex
  Double_t vertex[3] = {0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
  
  Double_t maxClusterEnergy = -1;
  Int_t maxClusterID        = -1;
  map<Long_t,Int_t> mapIsClusterAccepted;
  map<Long_t,Int_t> mapIsClusterAcceptedWithoutTrackMatch;
  // Loop over EMCal clusters
  for(Long_t i = 0; i < nclus; i++){
    
    AliVCluster* clus = NULL;
    if(fInputEvent->IsA()==AliESDEvent::Class()) clus = new AliESDCaloCluster(*(AliESDCaloCluster*)fInputEvent->GetCaloCluster(i));
    else if(fInputEvent->IsA()==AliAODEvent::Class()) clus = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(i));

    if(!clus) continue;
    if(!((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC,fWeightJetJetMC,i)){
      if(fProduceTreeEOverP && ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedBeforeTrackMatch() ) mapIsClusterAcceptedWithoutTrackMatch[i] = 1;
      delete clus;
      continue;
    }
    // TLorentzvector with cluster
    TLorentzVector clusterVector;
    clus->GetMomentum(clusterVector,vertex);
    
    TLorentzVector* tmpvec = new TLorentzVector();
    tmpvec->SetPxPyPzE(clusterVector.Px(),clusterVector.Py(),clusterVector.Pz(),clusterVector.E());
    
    // convert to AODConversionPhoton
    AliAODConversionPhoton *PhotonCandidate=new AliAODConversionPhoton(tmpvec);
    if(!PhotonCandidate){ delete clus; delete tmpvec; continue;}
    
    //determine maximum cluster energy in event
    if(fProduceCellIDPlots && (clus->E() > maxClusterEnergy)){
      maxClusterEnergy  = clus->E();
      maxClusterID      = (Int_t) i;
    }  
    if(fProduceTreeEOverP || fProduceCellIDPlots){
      mapIsClusterAccepted[i] = 1;
      mapIsClusterAcceptedWithoutTrackMatch[i] = 1;
    }

    // Flag Photon as CaloPhoton
    PhotonCandidate->SetIsCaloPhoton();
    PhotonCandidate->SetCaloClusterRef(i);
    PhotonCandidate->SetLeadingCellID(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FindLargestCellInCluster(clus,fInputEvent));
    // get MC label
    if(fIsMC> 0){
      Int_t* mclabelsCluster = clus->GetLabels();
      PhotonCandidate->SetNCaloPhotonMCLabels(clus->GetNLabels());
//       cout << clus->GetNLabels() << endl;
      if (clus->GetNLabels()>0){
        for (Int_t k =0; k< (Int_t)clus->GetNLabels(); k++){
          if (k< 50)PhotonCandidate->SetCaloPhotonMCLabel(k,mclabelsCluster[k]);
//           Int_t pdgCode = fMCEvent->Particle(mclabelsCluster[k])->GetPdgCode();
//           cout << "label " << k << "\t" << mclabelsCluster[k] << " pdg code: " << pdgCode << endl;
        }
      }
    }
    
    fIsFromMBHeader = kTRUE; 
    fIsOverlappingWithOtherHeader = kFALSE;
    // test whether largest contribution to cluster orginates in added signals
    if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() > 0){
      if ( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 0) fIsFromMBHeader = kFALSE;
      if (clus->GetNLabels()>1){
        Int_t* mclabelsCluster = clus->GetLabels();
        for (Int_t l = 1; l < (Int_t)clus->GetNLabels(); l++ ){
          if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(mclabelsCluster[l], fMCEvent, fInputEvent) == 0) fIsOverlappingWithOtherHeader = kTRUE;
        }
      } 
    }


    if (fIsOverlappingWithOtherHeader) fHistoClusOverlapHeadersGammaPt[fiCut]->Fill(PhotonCandidate->Pt(), fWeightJetJetMC);
    if (fIsFromMBHeader && !fIsOverlappingWithOtherHeader){
      fHistoClusGammaPt[fiCut]->Fill(PhotonCandidate->Pt(), fWeightJetJetMC);
      fHistoClusGammaE[fiCut]->Fill(PhotonCandidate->E(), fWeightJetJetMC);
      if (!fDoLightOutput && fDoClusterQA > 0) fHistoClusGammaPtM02[fiCut]->Fill(PhotonCandidate->Pt(), clus->GetM02(), fWeightJetJetMC);
      if(fIsMC> 0){
        if(fInputEvent->IsA()==AliESDEvent::Class()){
          ProcessTrueClusterCandidates(PhotonCandidate,clus);
        } else {
          ProcessTrueClusterCandidatesAOD(PhotonCandidate,clus);
        }
      }

      fClusterCandidates->Add(PhotonCandidate); // if no second loop is required add to events good gammas
    } else{
      delete PhotonCandidate;
    }
    
    delete clus;
    delete tmpvec;
  }
  
  if(fProduceCellIDPlots){
    for(Long_t i = 0; i < nclus; i++){
      if( mapIsClusterAccepted[i] != 1 ) continue;

      AliVCluster* clus = NULL;
      if(fInputEvent->IsA()==AliESDEvent::Class()) clus = new AliESDCaloCluster(*(AliESDCaloCluster*)fInputEvent->GetCaloCluster(i));
      else if(fInputEvent->IsA()==AliAODEvent::Class()) clus = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(i));
      if(!clus) continue;

      Int_t cellID = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FindLargestCellInCluster(clus,fInputEvent);
      fHistCellIDvsClusterEnergy[fiCut]->Fill(clus->E(),cellID);
      if (maxClusterID == i && maxClusterID > -1 ) fHistCellIDvsClusterEnergyMax[fiCut]->Fill(maxClusterEnergy,cellID);
      delete clus;
    }
  }

  if(fProduceTreeEOverP){
    AliESDEvent *esdev  = dynamic_cast<AliESDEvent*>(fInputEvent);
    AliAODEvent *aodev  = 0;
    Bool_t isESD        = kTRUE;
    if (!esdev) {
      isESD             = kFALSE;
      aodev             = dynamic_cast<AliAODEvent*>(fInputEvent);
      if (!aodev) {
        AliError("Task needs AOD or ESD event...");
      }
    }

    AliESDtrackCuts *EsdTrackCuts = 0x0;
    if(esdev){
      // Using standard function for setting Cuts
      Int_t runNumber = fInputEvent->GetRunNumber();
      // if LHC11a or earlier or if LHC13g or if LHC12a-i -> use 2010 cuts
      if( (runNumber<=146860) || (runNumber>=197470 && runNumber<=197692) || (runNumber>=172440 && runNumber<=193766) ){
        EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();
      // else if run2 data use 2015 PbPb cuts
      }else if (runNumber>=209122){
        EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2015PbPb();
      // else use 2011 version of track cuts
      }else{
        EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
      }
      EsdTrackCuts->SetMaxDCAToVertexZ(2);
      EsdTrackCuts->SetEtaRange(-0.8, 0.8);
      EsdTrackCuts->SetPtRange(0.15);
    }

    for(Long_t i = 0; i < nclus; i++){
      if( mapIsClusterAcceptedWithoutTrackMatch[i] != 1 ) continue;

      AliVCluster* clus = NULL;
      if(fInputEvent->IsA()==AliESDEvent::Class()) clus = new AliESDCaloCluster(*(AliESDCaloCluster*)fInputEvent->GetCaloCluster(i));
      else if(fInputEvent->IsA()==AliAODEvent::Class()) clus = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(i));
      if(!clus) continue;

      fClusterE = clus->E();
      fClusterM02 = clus->GetM02();
      fClusterM20 = clus->GetM20();
      fClusterLeadCellID = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FindLargestCellInCluster(clus,fInputEvent);

      Int_t labelTrackMatch = -1;
      if(!((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetHighestPtMatchedTrackToCluster(fInputEvent,clus,labelTrackMatch)){
        delete clus;
        continue;
      }

      AliVTrack* currTrack  = dynamic_cast<AliVTrack*>(fInputEvent->GetTrack(labelTrackMatch));
      if(!currTrack){
        delete clus;
        continue;
      }
      if(esdev){
        AliESDtrack *esdt = dynamic_cast<AliESDtrack*>(currTrack);
        if(!EsdTrackCuts->AcceptTrack(esdt)){
          delete clus;
          continue;
        }
        if(esdt->Pt()<1.){
          delete clus;
          continue;
        }
        fClusterEP = fClusterE/esdt->P();
        fTrackPt = esdt->Pt();
      }else if(aodev){
        AliAODTrack *aodt = dynamic_cast<AliAODTrack*>(currTrack);
        if(!aodt->IsHybridGlobalConstrainedGlobal()){
          delete clus;
          continue;
        }
        if(TMath::Abs(aodt->Eta())>0.8){
          delete clus;
          continue;
        }
        if(aodt->Pt()<1.){
          delete clus;
          continue;
        }
        fClusterEP = fClusterE/aodt->P();
        fTrackPt = aodt->Pt();
      }

      AliPIDResponse* pidResponse = ((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetPIDResponse();
      if(!pidResponse){
        delete clus;
        continue;
      }

      Float_t temp = TMath::Abs(pidResponse->NumberOfSigmasTPC(currTrack,AliPID::kElectron));
      if(temp<10.){
        fTrackPID_e = temp*10;
      }else fTrackPID_e = 99;

      temp = TMath::Abs(pidResponse->NumberOfSigmasTPC(currTrack,AliPID::kPion));
      if(temp<10.){
        fTrackPID_Pi = temp*10;
      }else fTrackPID_Pi = 99;

      temp = TMath::Abs(pidResponse->NumberOfSigmasTPC(currTrack,AliPID::kKaon));
      if(temp<10.){
        fTrackPID_K = temp*10;
      }else fTrackPID_K = 99;

      temp = TMath::Abs(pidResponse->NumberOfSigmasTPC(currTrack,AliPID::kProton));
      if(temp<10.){
        fTrackPID_P = temp*10;
      }else fTrackPID_P = 99;

      Float_t tempEta = -99999;
      Float_t tempPhi = -99999;
      ((AliCaloTrackMatcher*)((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetCaloTrackMatcherInstance())->GetTrackClusterMatchingResidual(currTrack->GetID(),clus->GetID(),tempEta,tempPhi);
      fDeltaEta = tempEta;
      fDeltaPhi = tempPhi;

      //determine isolation in cluster Et
      Float_t clsPos[3] = {0.,0.,0.};
      Float_t secondClsPos[3] = {0.,0.,0.};
      TLorentzVector clusterVector;

      clus->GetPosition(clsPos);
      TVector3 clsPosVec(clsPos);

      Float_t sum_Et = 0;
      for(Int_t j=0; j<nclus; j++){
        if( i == j ) continue;
        if( mapIsClusterAcceptedWithoutTrackMatch[j] != 1 ) continue;

        AliVCluster* secondClus = NULL;
        if(fInputEvent->IsA()==AliESDEvent::Class()) secondClus = new AliESDCaloCluster(*(AliESDCaloCluster*)fInputEvent->GetCaloCluster(j));
        else if(fInputEvent->IsA()==AliAODEvent::Class()) secondClus = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(j));
        if(!secondClus) continue;
        secondClus->GetPosition(secondClsPos);
        TVector3 secondClsPosVec(secondClsPos);

        Float_t dPhi = clsPosVec.DeltaPhi(secondClsPosVec);
        Float_t dEta = clsPosVec.Eta()-secondClsPosVec.Eta();
        if(TMath::Sqrt(dEta*dEta + dPhi*dPhi) < 0.2){
          secondClus->GetMomentum(clusterVector,vertex);
          sum_Et += clusterVector.Et();
        }
        delete secondClus;
      }
      fClusterIsoSumClusterEt = sum_Et;

      //determine isolation in track Et
      fClusterIsoSumTrackEt = ((AliCaloTrackMatcher*)((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetCaloTrackMatcherInstance())->SumTrackEtAroundCluster(fInputEvent,clus->GetID(),0.2);
      //remove Et from matched track
      TLorentzVector vecTrack;
      vecTrack.SetPxPyPzE(currTrack->Px(),currTrack->Py(),currTrack->Pz(),currTrack->E());
      fClusterIsoSumTrackEt -= vecTrack.Et();

      //get cluster classification
      if(fIsMC > 0) fClusterClassification = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClassifyClusterForTMEffi(clus,fInputEvent,fMCEvent,isESD);

      tClusterEOverP[fiCut]->Fill();
      delete clus;
    }
    mapIsClusterAcceptedWithoutTrackMatch.clear();
  }

  if(fProduceCellIDPlots || fProduceTreeEOverP) mapIsClusterAccepted.clear();

  if(fLocalDebugFlag == 2) EventDebugMethod();
  return;
}

//________________________________________________________________________
void AliAnalysisTaskGammaCalo::ProcessTrueClusterCandidates(AliAODConversionPhoton *TruePhotonCandidate, AliVCluster* clus)
{
    
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  TParticle *Photon = NULL;
  if (!TruePhotonCandidate->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set task will abort");
  if (!fDoLightOutput && fDoClusterQA > 0) fHistoTrueNLabelsInClus[fiCut]->Fill(TruePhotonCandidate->GetNCaloPhotonMCLabels(), fWeightJetJetMC);
  
  if (TruePhotonCandidate->GetCaloPhotonMCLabel(0) < 0) return;
  if (TruePhotonCandidate->GetNCaloPhotonMCLabels() > 0) Photon = fMCEvent->Particle(TruePhotonCandidate->GetCaloPhotonMCLabel(0));
    else return;
    
  if(Photon == NULL){
  //    cout << "no photon" << endl;
    return;
  }

  Int_t pdgCodeParticle = Photon->GetPdgCode();
  TruePhotonCandidate->SetCaloPhotonMCFlags(fMCEvent, fEnableSortForClusMC);
  
  // True Photon
  if(fIsFromMBHeader && !fIsOverlappingWithOtherHeader){
    if (TruePhotonCandidate->IsLargestComponentPhoton() || (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion())) {
      fHistoTrueClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
      if (!fDoLightOutput && fDoClusterQA > 0) fHistoTrueClusGammaPtM02[fiCut]->Fill(TruePhotonCandidate->Pt(), clus->GetM02(), fWeightJetJetMC);
    } else if (fDoClusterQA > 0) fHistoTrueClusEMNonLeadingPt[fiCut]->Fill(TruePhotonCandidate->Pt());
    if (fDoClusterQA > 0){
      if (TruePhotonCandidate->IsLargestComponentPhoton()){ 
        fHistoTrueClusUnConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
        fHistoTrueClusUnConvGammaMCPt[fiCut]->Fill(Photon->Pt(), fWeightJetJetMC);
        if (!fDoLightOutput) fHistoTrueClusUnConvGammaPtM02[fiCut]->Fill(TruePhotonCandidate->Pt(), clus->GetM02(), fWeightJetJetMC);
      }
      if (TruePhotonCandidate->IsLargestComponentElectron()) 
        fHistoTrueClusElectronPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
      if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
        fHistoTrueClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
        if(Photon->GetMother(0) > -1){
          TParticle* motherPart = (TParticle*)fMCEvent->Particle(Photon->GetMother(0));
          fHistoTrueClusConvGammaMCPt[fiCut]->Fill(motherPart->Pt(), fWeightJetJetMC);
        }
      }
      if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion() && TruePhotonCandidate->IsConversionFullyContained()) 
        fHistoTrueClusConvGammaFullyPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
      if (TruePhotonCandidate->IsMerged() || TruePhotonCandidate->IsMergedPartConv() || TruePhotonCandidate->IsDalitzMerged())
        fHistoTrueClusMergedGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
      if (TruePhotonCandidate->IsMergedPartConv())
        fHistoTrueClusMergedPartConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
      if (TruePhotonCandidate->IsDalitz()) 
        fHistoTrueClusDalitzPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
      if (TruePhotonCandidate->IsDalitzMerged()) 
        fHistoTrueClusDalitzMergedPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
      if (TruePhotonCandidate->IsPhotonWithElecMother()) 
        fHistoTrueClusPhotonFromElecMotherPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
      if (TruePhotonCandidate->IsShower()) 
        fHistoTrueClusShowerPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
      if (TruePhotonCandidate->IsSubLeadingEM())
        fHistoTrueClusSubLeadingPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
      fHistoTrueClusNParticles[fiCut]->Fill(TruePhotonCandidate->GetNCaloPhotonMCLabels(), fWeightJetJetMC);
      if (!TruePhotonCandidate->IsLargestComponentPhoton()) {
        FillPhotonBackgroundHist(TruePhotonCandidate,pdgCodeParticle);
        if (!fDoLightOutput) FillPhotonBackgroundM02Hist(TruePhotonCandidate,clus,pdgCodeParticle);
      }
      if (!(TruePhotonCandidate->IsLargestComponentPhoton() || (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion())) ) {
        FillPhotonPlusConversionBackgroundHist(TruePhotonCandidate,pdgCodeParticle);
        if (!fDoLightOutput) FillPhotonPlusConversionBackgroundM02Hist(TruePhotonCandidate,clus,pdgCodeParticle);
      }
      Int_t motherLab = Photon->GetMother(0);
      if (motherLab > -1){
        if (TruePhotonCandidate->IsLargestComponentPhoton()){
          if (CheckVectorForDoubleCount(fVectorDoubleCountTrueClusterGammas,motherLab)){
            fHistoDoubleCountTrueClusterGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),(Double_t)0,fWeightJetJetMC);
            FillMultipleCountMap(fMapMultipleCountTrueClusterGammas,motherLab);
          }
        }
        Int_t grandMotherLab = fMCEvent->Particle(motherLab)->GetMother(0);
        if (grandMotherLab > -1){
          if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
            if (CheckVectorForDoubleCount(fVectorDoubleCountTrueClusterGammas,grandMotherLab)){
              fHistoDoubleCountTrueClusterGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),(Double_t)1,fWeightJetJetMC);
              FillMultipleCountMap(fMapMultipleCountTrueClusterGammas,grandMotherLab);
            }
          }
        }
      }
    }
    
    if(!fDoLightOutput){
      Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, TruePhotonCandidate->GetCaloPhotonMCLabel(0), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
      if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion())
        isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, Photon->GetMother(0), mcProdVtxX, mcProdVtxY, mcProdVtxZ);

      if(isPrimary){
        // filling primary histograms
        if (TruePhotonCandidate->IsLargestComponentPhoton()){
          fHistoTruePrimaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
          fHistoTruePrimaryClusGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt(), fWeightJetJetMC); // Allways Filled
        }
        if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
          fHistoTruePrimaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
          fHistoTruePrimaryClusConvGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt(), fWeightJetJetMC); // Allways Filled
        }

      } else {
        // filling secondary histograms
        Int_t secondaryClass    = -1;
        if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion())
          secondaryClass        = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->SecondaryClassificationPhoton( Photon, fMCEvent, kTRUE);
        else
          secondaryClass        = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->SecondaryClassificationPhoton( Photon, fMCEvent, kFALSE);

        // all secondaries
        if (TruePhotonCandidate->IsLargestComponentPhoton()){
          if (secondaryClass == 2) {
            fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),0.,fWeightJetJetMC);
            fHistoTrueSecondaryClusGammaMCPt[fiCut]->Fill(Photon->Pt(),0.,fWeightJetJetMC);
            fHistoTrueSecondaryClusGammaFromXFromK0sMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),fWeightJetJetMC);
          } else if (secondaryClass == 5) {
            fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),1.,fWeightJetJetMC);
            fHistoTrueSecondaryClusGammaMCPt[fiCut]->Fill(Photon->Pt(),1.,fWeightJetJetMC);
            fHistoTrueSecondaryClusGammaFromXFromK0lMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),fWeightJetJetMC);
          } else if (secondaryClass == 3) {
            fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),2.,fWeightJetJetMC);
            fHistoTrueSecondaryClusGammaMCPt[fiCut]->Fill(Photon->Pt(),2.,fWeightJetJetMC);
            fHistoTrueSecondaryClusGammaFromXFromLambdaMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),fWeightJetJetMC);
          } else if (secondaryClass == 4) {
            fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3.,fWeightJetJetMC);
            fHistoTrueSecondaryClusGammaMCPt[fiCut]->Fill(Photon->Pt(),3.,fWeightJetJetMC);
          } else {
            fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),4.,fWeightJetJetMC);
            fHistoTrueSecondaryClusGammaMCPt[fiCut]->Fill(Photon->Pt(),4.,fWeightJetJetMC);
          }
        }
        if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
          if (secondaryClass == 2) {
            fHistoTrueSecondaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),0.,fWeightJetJetMC);
            fHistoTrueSecondaryClusConvGammaMCPt[fiCut]->Fill(Photon->Pt(),0.,fWeightJetJetMC);
            fHistoTrueSecondaryClusConvGammaFromXFromK0sMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),fWeightJetJetMC);
          } else if (secondaryClass == 5) {
            fHistoTrueSecondaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),1.,fWeightJetJetMC);
            fHistoTrueSecondaryClusConvGammaMCPt[fiCut]->Fill(Photon->Pt(),1.,fWeightJetJetMC);
            fHistoTrueSecondaryClusConvGammaFromXFromK0lMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),fWeightJetJetMC);
          } else if (secondaryClass == 3) {
            fHistoTrueSecondaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),2.,fWeightJetJetMC);
            fHistoTrueSecondaryClusConvGammaMCPt[fiCut]->Fill(Photon->Pt(),2.,fWeightJetJetMC);
            fHistoTrueSecondaryClusConvGammaFromXFromLambdaMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),fWeightJetJetMC);
          } else if (secondaryClass == 4) {
            fHistoTrueSecondaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3.,fWeightJetJetMC);
            fHistoTrueSecondaryClusConvGammaMCPt[fiCut]->Fill(Photon->Pt(),3.,fWeightJetJetMC);
          } else {
            fHistoTrueSecondaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),4.,fWeightJetJetMC);
            fHistoTrueSecondaryClusConvGammaMCPt[fiCut]->Fill(Photon->Pt(),4.,fWeightJetJetMC);
          }
        }
      }
    }
  }
  return;
}


//________________________________________________________________________
void AliAnalysisTaskGammaCalo::ProcessTrueClusterCandidatesAOD(AliAODConversionPhoton *TruePhotonCandidate, AliVCluster* clus)
{
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  AliAODMCParticle *Photon = NULL;
  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (!fDoLightOutput && fDoClusterQA > 0) fHistoTrueNLabelsInClus[fiCut]->Fill(TruePhotonCandidate->GetNCaloPhotonMCLabels(), fWeightJetJetMC);
  if (AODMCTrackArray){
    if (!TruePhotonCandidate->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set task will abort");
    if (TruePhotonCandidate->GetNCaloPhotonMCLabels()>0) Photon = (AliAODMCParticle*) AODMCTrackArray->At(TruePhotonCandidate->GetCaloPhotonMCLabel(0));
      else return;
  } else {
    AliInfo("AODMCTrackArray could not be loaded");
    return;
  }

  if(Photon == NULL){
  //  cout << "no photon" << endl;
    return;
  }
  Int_t pdgCodeParticle = Photon->GetPdgCode();
  TruePhotonCandidate->SetCaloPhotonMCFlagsAOD(fInputEvent, fEnableSortForClusMC);
  
  // True Photon
  if(fIsFromMBHeader && !fIsOverlappingWithOtherHeader){
    if (TruePhotonCandidate->IsLargestComponentPhoton() || (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion())) {
      fHistoTrueClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
      if (!fDoLightOutput && fDoClusterQA > 0) fHistoTrueClusGammaPtM02[fiCut]->Fill(TruePhotonCandidate->Pt(), clus->GetM02(), fWeightJetJetMC);
    } else if (fDoClusterQA > 0) fHistoTrueClusEMNonLeadingPt[fiCut]->Fill(TruePhotonCandidate->Pt());
    if (fDoClusterQA > 0){
      if (TruePhotonCandidate->IsLargestComponentPhoton()) {
        fHistoTrueClusUnConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
        fHistoTrueClusUnConvGammaMCPt[fiCut]->Fill(Photon->Pt(), fWeightJetJetMC);
        if (!fDoLightOutput) fHistoTrueClusUnConvGammaPtM02[fiCut]->Fill(TruePhotonCandidate->Pt(), clus->GetM02(), fWeightJetJetMC);
      }
      if (TruePhotonCandidate->IsLargestComponentElectron()) 
        fHistoTrueClusElectronPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
      if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()) {
        fHistoTrueClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
        AliAODMCParticle *motherPart = (AliAODMCParticle*) AODMCTrackArray->At(Photon->GetMother());
        fHistoTrueClusConvGammaMCPt[fiCut]->Fill(motherPart->Pt(), fWeightJetJetMC);
      }
      if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion() && TruePhotonCandidate->IsConversionFullyContained()) 
        fHistoTrueClusConvGammaFullyPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
      if (TruePhotonCandidate->IsMerged() || TruePhotonCandidate->IsMergedPartConv() || TruePhotonCandidate->IsDalitzMerged())
        fHistoTrueClusMergedGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
      if (TruePhotonCandidate->IsMergedPartConv())
        fHistoTrueClusMergedPartConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
      if (TruePhotonCandidate->IsDalitz()) 
        fHistoTrueClusDalitzPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
      if (TruePhotonCandidate->IsDalitzMerged()) 
        fHistoTrueClusDalitzMergedPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
      if (TruePhotonCandidate->IsPhotonWithElecMother()) 
        fHistoTrueClusPhotonFromElecMotherPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
      if (TruePhotonCandidate->IsShower()) 
        fHistoTrueClusShowerPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
      if (TruePhotonCandidate->IsSubLeadingEM())
        fHistoTrueClusSubLeadingPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
      fHistoTrueClusNParticles[fiCut]->Fill(TruePhotonCandidate->GetNCaloPhotonMCLabels(), fWeightJetJetMC);
      
      if (!TruePhotonCandidate->IsLargestComponentPhoton()) {
        FillPhotonBackgroundHist(TruePhotonCandidate,pdgCodeParticle);
        if (!fDoLightOutput) FillPhotonBackgroundM02Hist(TruePhotonCandidate,clus,pdgCodeParticle);
      }
      if (!(TruePhotonCandidate->IsLargestComponentPhoton() || (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion())) ) {
        FillPhotonPlusConversionBackgroundHist(TruePhotonCandidate,pdgCodeParticle);
        if (!fDoLightOutput) FillPhotonPlusConversionBackgroundM02Hist(TruePhotonCandidate,clus,pdgCodeParticle);
      }
      Int_t motherLab = Photon->GetMother();
      if (motherLab > -1){
        if (TruePhotonCandidate->IsLargestComponentPhoton()){
          if (CheckVectorForDoubleCount(fVectorDoubleCountTrueClusterGammas,motherLab)){
            fHistoDoubleCountTrueClusterGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),(Double_t)0,fWeightJetJetMC);
            FillMultipleCountMap(fMapMultipleCountTrueClusterGammas,motherLab);
          }
        }
        Int_t grandMotherLab = ((AliAODMCParticle*) AODMCTrackArray->At(motherLab))->GetMother();
        if (grandMotherLab > -1){
          if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
            if (CheckVectorForDoubleCount(fVectorDoubleCountTrueClusterGammas,grandMotherLab)){
              fHistoDoubleCountTrueClusterGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),(Double_t)1,fWeightJetJetMC);
              FillMultipleCountMap(fMapMultipleCountTrueClusterGammas,grandMotherLab);
            }
          }
        }
      }
    }
    
    if(!fDoLightOutput){
      Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, Photon, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
      if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
        if (Photon->GetMother()> -1){
          AliAODMCParticle *Mother  = (AliAODMCParticle*) AODMCTrackArray->At(Photon->GetMother());
          isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, Mother, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
        }
      }
      if(isPrimary){
        if (TruePhotonCandidate->IsLargestComponentPhoton()){
          fHistoTruePrimaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
          fHistoTruePrimaryClusGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt(), fWeightJetJetMC); // Allways Filled
        }
        if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
          fHistoTruePrimaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
          fHistoTruePrimaryClusConvGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt(), fWeightJetJetMC); // Allways Filled
        }

      } else {
        // filling secondary histograms
        Int_t secondaryClass    = -1;
        if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion())
          secondaryClass        = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->SecondaryClassificationPhotonAOD( Photon, AODMCTrackArray, kTRUE);
        else
          secondaryClass        = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->SecondaryClassificationPhotonAOD( Photon, AODMCTrackArray, kFALSE);

        // all secondaries
        if (TruePhotonCandidate->IsLargestComponentPhoton()){
          if (secondaryClass == 2) {
            fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),0.,fWeightJetJetMC);
            fHistoTrueSecondaryClusGammaMCPt[fiCut]->Fill(Photon->Pt(),0.,fWeightJetJetMC);
            fHistoTrueSecondaryClusGammaFromXFromK0sMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),fWeightJetJetMC);
          } else if (secondaryClass == 5) {
            fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),1.,fWeightJetJetMC);
            fHistoTrueSecondaryClusGammaMCPt[fiCut]->Fill(Photon->Pt(),1.,fWeightJetJetMC);
            fHistoTrueSecondaryClusGammaFromXFromK0lMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),fWeightJetJetMC);
          } else if (secondaryClass == 3) {
            fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),2.,fWeightJetJetMC);
            fHistoTrueSecondaryClusGammaMCPt[fiCut]->Fill(Photon->Pt(),2.,fWeightJetJetMC);
            fHistoTrueSecondaryClusGammaFromXFromLambdaMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),fWeightJetJetMC);
          } else if (secondaryClass == 4) {
            fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3.,fWeightJetJetMC);
            fHistoTrueSecondaryClusGammaMCPt[fiCut]->Fill(Photon->Pt(),3.,fWeightJetJetMC);
          } else {
            fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),4.,fWeightJetJetMC);
            fHistoTrueSecondaryClusGammaMCPt[fiCut]->Fill(Photon->Pt(),4.,fWeightJetJetMC);
          }
        }
        if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
          if (secondaryClass == 2) {
            fHistoTrueSecondaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),0.,fWeightJetJetMC);
            fHistoTrueSecondaryClusConvGammaMCPt[fiCut]->Fill(Photon->Pt(),0.,fWeightJetJetMC);
            fHistoTrueSecondaryClusConvGammaFromXFromK0sMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),fWeightJetJetMC);
          } else if (secondaryClass == 5) {
            fHistoTrueSecondaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),1.,fWeightJetJetMC);
            fHistoTrueSecondaryClusConvGammaMCPt[fiCut]->Fill(Photon->Pt(),1.,fWeightJetJetMC);
            fHistoTrueSecondaryClusConvGammaFromXFromK0lMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),fWeightJetJetMC);
          } else if (secondaryClass == 3) {
            fHistoTrueSecondaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),2.,fWeightJetJetMC);
            fHistoTrueSecondaryClusConvGammaMCPt[fiCut]->Fill(Photon->Pt(),2.,fWeightJetJetMC);
            fHistoTrueSecondaryClusConvGammaFromXFromLambdaMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),fWeightJetJetMC);
          } else if (secondaryClass == 4) {
            fHistoTrueSecondaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3.,fWeightJetJetMC);
            fHistoTrueSecondaryClusConvGammaMCPt[fiCut]->Fill(Photon->Pt(),3.,fWeightJetJetMC);
          } else {
            fHistoTrueSecondaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),4.,fWeightJetJetMC);
            fHistoTrueSecondaryClusConvGammaMCPt[fiCut]->Fill(Photon->Pt(),4.,fWeightJetJetMC);
          }
        }
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskGammaCalo::ProcessAODMCParticles()
{
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (AODMCTrackArray == NULL) return;
  
  // Loop over all primary MC particle
  for(Long_t i = 0; i < AODMCTrackArray->GetEntriesFast(); i++) {
    
    AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(i));
    if (!particle) continue;

    Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, particle, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
    if (isPrimary){
      
      Int_t isMCFromMBHeader = -1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        isMCFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
        if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
      }
      
      if(!fDoLightOutput){
        if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(particle,AODMCTrackArray)){
          fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC); // All MC Gamma
          if(particle->GetMother() >-1){ // Meson Decay Gamma
            switch((static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetMother())))->GetPdgCode()){
            case 111: // Pi0
              fHistoMCDecayGammaPi0Pt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
              break;
            case 113: // Rho0
              fHistoMCDecayGammaRhoPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
              break;
            case 221: // Eta
              fHistoMCDecayGammaEtaPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
              break;
            case 223: // Omega
              fHistoMCDecayGammaOmegaPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
              break;
            case 331: // Eta'
              fHistoMCDecayGammaEtapPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
              break;
            case 333: // Phi
              fHistoMCDecayGammaPhiPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
              break;
            case 3212: // Sigma
              fHistoMCDecayGammaSigmaPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
              break;
            }
          }
        }
      }
      // Converted MC Gamma
      if(fDoMesonAnalysis){
        Double_t mesonY = 10.;
        if(particle->E() - particle->Pz() == 0 || particle->E() + particle->Pz() == 0){
          mesonY=10.-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
        } else {
          mesonY = 0.5*(TMath::Log((particle->E()+particle->Pz()) / (particle->E()-particle->Pz())))-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
        }
        
        if (TMath::Abs(mesonY) < ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetRapidityCutValue()){
          if ( particle->GetPdgCode() == 211 ){  // positve pions
            fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 0., fWeightJetJetMC); 
          } else if ( particle->GetPdgCode() == -211 ){  // negative pions
            fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 1., fWeightJetJetMC); 
          } else if ( particle->GetPdgCode() == 321 ){  // positve kaons
            fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 2., fWeightJetJetMC); 
          } else if ( particle->GetPdgCode() == -321 ){  // negative kaons
            fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 3., fWeightJetJetMC); 
          } else if ( TMath::Abs(particle->GetPdgCode()) == 310 ){  // K0s
            fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 4., fWeightJetJetMC); 
          } else if ( TMath::Abs(particle->GetPdgCode()) == 130 ){  // K0l
            fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 5., fWeightJetJetMC); 
          } else if ( TMath::Abs(particle->GetPdgCode()) == 3122 ){  // Lambda/ AntiLambda
            fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 6., fWeightJetJetMC); 
          }
        }
      
        if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
          ->MesonIsSelectedAODMC(particle,AODMCTrackArray,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
          AliAODMCParticle* daughter0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetDaughter(0)));
          AliAODMCParticle* daughter1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetDaughter(1)));
          Float_t weighted= 1;
          if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent)){
            if (particle->Pt()>0.005){
              weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(i, 0x0, fInputEvent);
            }
          }
          
          Double_t mesonY = 10.;
          if(particle->E() - particle->Pz() == 0 || particle->E() + particle->Pz() == 0){
            mesonY=10.-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
          } else{
            mesonY = 0.5*(TMath::Log((particle->E()+particle->Pz()) / (particle->E()-particle->Pz())))-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
          }
          Double_t alpha = -1;
          if (particle->GetPdgCode() == 111 || particle->GetPdgCode() == 221){
            alpha = TMath::Abs((daughter0->E() - daughter1->E()))/(daughter0->E() + daughter1->E());
          }
          
          if(particle->GetPdgCode() == 111){
            fHistoMCPi0Pt[fiCut]->Fill(particle->Pt(),weighted* fWeightJetJetMC); // All MC Pi0
            fHistoMCPi0WOWeightPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
            if (fIsMC > 1) fHistoMCPi0WOEvtWeightPt[fiCut]->Fill(particle->Pt());
            if (fDoMesonQA > 0 && fDoMesonQA < 3){
              fHistoMCPi0PtY[fiCut]->Fill(particle->Pt(),mesonY,weighted* fWeightJetJetMC); 
              fHistoMCPi0PtAlpha[fiCut]->Fill(particle->Pt(),alpha, fWeightJetJetMC);
              if (fIsMC == 2) fHistoMCPi0PtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),fWeightJetJetMC);
            }
          } else if(particle->GetPdgCode() == 221){
            fHistoMCEtaPt[fiCut]->Fill(particle->Pt(),weighted* fWeightJetJetMC); // All MC Eta
            fHistoMCEtaWOWeightPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
            if (fIsMC > 1) fHistoMCEtaWOEvtWeightPt[fiCut]->Fill(particle->Pt());
            if (fDoMesonQA > 0 && fDoMesonQA < 3){
              fHistoMCEtaPtY[fiCut]->Fill(particle->Pt(),mesonY,weighted* fWeightJetJetMC); 
              fHistoMCEtaPtAlpha[fiCut]->Fill(particle->Pt(),alpha, fWeightJetJetMC);
              if (fIsMC == 2) fHistoMCEtaPtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),fWeightJetJetMC);
            }
          }
          
          // Check the acceptance for both gammas
          if( ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(daughter0,AODMCTrackArray) &&
              ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(daughter1,AODMCTrackArray) ){
            
            if(particle->GetPdgCode() == 111){
              fHistoMCPi0InAccPt[fiCut]->Fill(particle->Pt(),weighted* fWeightJetJetMC); // MC Pi0 with gamma in acc
              if(fIsMC > 1) fHistoMCPi0WOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Pi0 with gamma in acc
            } else if(particle->GetPdgCode() == 221){
              fHistoMCEtaInAccPt[fiCut]->Fill(particle->Pt(),weighted* fWeightJetJetMC); // MC Eta with gamma in acc
              if(fIsMC > 1) fHistoMCEtaWOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Eta with gamma in acc
            }
          }
        }
      }
    // fill secondaries
    } else {
      Int_t isMCFromMBHeader = -1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        isMCFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
        if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
      }
      
      if(!fDoLightOutput) {
        if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(particle,AODMCTrackArray)){
          if(particle->GetMother() >-1){
            AliAODMCParticle *tmpMother = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetMother()));
            if(tmpMother->GetMother() >-1){
              AliAODMCParticle *tmpGrandMother = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(tmpMother->GetMother()));
              if(tmpGrandMother->GetPdgCode() == 310) {
                fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),0.,fWeightJetJetMC);
              } else if (tmpGrandMother->GetPdgCode() == 130) {
                fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),1.,fWeightJetJetMC);
              } else if (tmpGrandMother->GetPdgCode() == 3122) {
                fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),2.,fWeightJetJetMC);
              } else if (tmpGrandMother->GetPdgCode() == 221) {
                fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),3.,fWeightJetJetMC);
              } else {
                if( !(TMath::Abs(tmpMother->GetPdgCode()) == 11 && tmpGrandMother->GetPdgCode() == 22) )
                  fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),4.,fWeightJetJetMC);
              }
            } else {
              fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),4.,fWeightJetJetMC);
            }
          } else {
            fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),4.,fWeightJetJetMC);
          }
        }
      }
      
      if(fDoMesonAnalysis){
        if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedAODMC(particle,AODMCTrackArray,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
          AliAODMCParticle* daughter0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetDaughter(0)));
          AliAODMCParticle* daughter1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetDaughter(1)));
          AliAODMCParticle* mother = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetMother()));
          Int_t pdgCode = mother->GetPdgCode();
          if(particle->GetPdgCode() == 111){              
            Int_t source = GetSourceClassification(111,pdgCode);
            fHistoMCSecPi0PtvsSource[fiCut]->Fill(particle->Pt(),source,fWeightJetJetMC); // All MC Pi0
            fHistoMCSecPi0Source[fiCut]->Fill(pdgCode);
          }else if(particle->GetPdgCode() == 221){
            
            fHistoMCSecEtaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // All MC Pi0
            fHistoMCSecEtaSource[fiCut]->Fill(pdgCode);
          }
          
          // check if conversion where within acceptance
          if( ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(daughter0,AODMCTrackArray) &&
              ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(daughter1,AODMCTrackArray)){
            if(particle->GetPdgCode() == 111){              
              Int_t source = GetSourceClassification(111,pdgCode);
              fHistoMCSecPi0InAccPtvsSource[fiCut]->Fill(particle->Pt(),source,fWeightJetJetMC); // All MC Pi0
            }
          }    
        }
      }      
    }  
  }
}
//________________________________________________________________________
void AliAnalysisTaskGammaCalo::ProcessMCParticles()
{
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();
  
  // Loop over all primary MC particle
  for(Long_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {
    if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, i, mcProdVtxX, mcProdVtxY, mcProdVtxZ)){

      TParticle* particle = (TParticle *)fMCEvent->Particle(i);
      if (!particle) continue;
      
      Int_t isMCFromMBHeader = -1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        isMCFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
        if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
      }
      
      if(!fDoLightOutput){
        if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(particle,fMCEvent)){
          fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC); // All MC Gamma
          if(particle->GetMother(0) >-1){ // Meson Decay Gamma
            switch(fMCEvent->Particle(particle->GetMother(0))->GetPdgCode()){
            case 111: // Pi0
              fHistoMCDecayGammaPi0Pt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
              break;
            case 113: // Rho0
              fHistoMCDecayGammaRhoPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
              break;
            case 221: // Eta
              fHistoMCDecayGammaEtaPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
              break;
            case 223: // Omega
              fHistoMCDecayGammaOmegaPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
              break;
            case 331: // Eta'
              fHistoMCDecayGammaEtapPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
              break;
            case 333: // Phi
              fHistoMCDecayGammaPhiPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
              break;
            case 3212: // Sigma
              fHistoMCDecayGammaSigmaPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
              break;
            }
          }
        }
      }
      if(fDoMesonAnalysis){
      
        Double_t mesonY = 10.;
        if(particle->Energy() - particle->Pz() == 0 || particle->Energy() + particle->Pz() == 0){
          mesonY=10.-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
        } else{
          mesonY = 0.5*(TMath::Log((particle->Energy()+particle->Pz()) / (particle->Energy()-particle->Pz())))-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
        }
        
        if (TMath::Abs(mesonY) < ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetRapidityCutValue()){
          if ( particle->GetPdgCode() == 211 ){  // positve pions
            fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 0., fWeightJetJetMC); 
          } else if ( particle->GetPdgCode() == -211 ){  // negative pions
            fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 1., fWeightJetJetMC); 
          } else if ( particle->GetPdgCode() == 321 ){  // positve kaons
            fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 2., fWeightJetJetMC); 
          } else if ( particle->GetPdgCode() == -321 ){  // negative kaons
            fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 3., fWeightJetJetMC); 
          } else if ( TMath::Abs(particle->GetPdgCode()) == 310 ){  // K0s
            fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 4., fWeightJetJetMC); 
          } else if ( TMath::Abs(particle->GetPdgCode()) == 130 ){  // K0l
            fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 5., fWeightJetJetMC); 
          } else if ( TMath::Abs(particle->GetPdgCode()) == 3122 ){  // Lambda/ AntiLambda
            fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 6., fWeightJetJetMC); 
          }
        }
        
        if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
          ->MesonIsSelectedMC(particle,fMCEvent,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
          TParticle* daughter0 = (TParticle*)fMCEvent->Particle(particle->GetFirstDaughter());
          TParticle* daughter1 = (TParticle*)fMCEvent->Particle(particle->GetLastDaughter());
          
          Float_t weighted= 1;
          if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent)){
            if (particle->Pt()>0.005){
              weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(i, fMCEvent, fInputEvent);
            }
          }
          Double_t mesonY = 10.;
          if(particle->Energy() - particle->Pz() == 0 || particle->Energy() + particle->Pz() == 0){
            mesonY=10.-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
          } else{
            mesonY = 0.5*(TMath::Log((particle->Energy()+particle->Pz()) / (particle->Energy()-particle->Pz())))-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
          }
    
          Double_t alpha = -1;
          if (particle->GetPdgCode() == 111 || particle->GetPdgCode() == 221){
            alpha = TMath::Abs((daughter0->Energy() - daughter1->Energy()))/(daughter0->Energy() + daughter1->Energy());
          }
    
          if(particle->GetPdgCode() == 111){
            fHistoMCPi0Pt[fiCut]->Fill(particle->Pt(),weighted* fWeightJetJetMC); // All MC Pi0
            fHistoMCPi0WOWeightPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
            if (fIsMC > 1)fHistoMCPi0WOEvtWeightPt[fiCut]->Fill(particle->Pt());
            if (fDoMesonQA > 0 && fDoMesonQA < 3){
              fHistoMCPi0PtY[fiCut]->Fill(particle->Pt(),mesonY,weighted* fWeightJetJetMC); 
              fHistoMCPi0PtAlpha[fiCut]->Fill(particle->Pt(),alpha, fWeightJetJetMC); 
              if (fIsMC == 2) fHistoMCPi0PtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),fWeightJetJetMC);
            }
          } else if(particle->GetPdgCode() == 221){
            fHistoMCEtaPt[fiCut]->Fill(particle->Pt(),weighted* fWeightJetJetMC); // All MC Eta
            fHistoMCEtaWOWeightPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
            if (fIsMC > 1)fHistoMCEtaWOEvtWeightPt[fiCut]->Fill(particle->Pt());
            if (fDoMesonQA > 0 && fDoMesonQA < 3){
              fHistoMCEtaPtY[fiCut]->Fill(particle->Pt(),mesonY,weighted* fWeightJetJetMC); 
              fHistoMCEtaPtAlpha[fiCut]->Fill(particle->Pt(),alpha, fWeightJetJetMC); 
              if (fIsMC == 2) fHistoMCEtaPtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),fWeightJetJetMC);
            }
          }
        
          // Check the acceptance for both gammas & whether they are counted as primaries as well
          Bool_t kDaughter0IsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, particle->GetFirstDaughter(), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
          Bool_t kDaughter1IsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, particle->GetLastDaughter(), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
          
          if( kDaughter0IsPrim && kDaughter1IsPrim &&
            ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter0,fMCEvent) &&
            ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter1,fMCEvent) ){
            if(particle->GetPdgCode() == 111){
              fHistoMCPi0InAccPt[fiCut]->Fill(particle->Pt(),weighted* fWeightJetJetMC); // MC Pi0 with gamma in acc
              if(fIsMC > 1) fHistoMCPi0WOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Pi0 with gamma in acc
            } else if(particle->GetPdgCode() == 221){
              fHistoMCEtaInAccPt[fiCut]->Fill(particle->Pt(),weighted* fWeightJetJetMC); // MC Eta with gamma in acc
              if(fIsMC > 1) fHistoMCEtaWOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Eta with gamma in acc
            }
          }
        }
      }
    // fill secondary histograms  
    } else {
      TParticle* particle = (TParticle *)fMCEvent->Particle(i);
      if (!particle) continue;
  
      Int_t isMCFromMBHeader = -1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        isMCFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
        if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
      }
      
      if(!fDoLightOutput) {
        if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(particle,fMCEvent)){
          if (particle->GetMother(0) > -1 && fMCEvent->Particle(particle->GetMother(0))->GetMother(0) > -1) {
            if (fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 310){
              fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),0.,fWeightJetJetMC);
            } else if (fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 130) {
              fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),1.,fWeightJetJetMC);
            } else if (fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 3122) {
              fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),2.,fWeightJetJetMC);
            } else if (fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 221) {
              fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),3.,fWeightJetJetMC);
            } else {
              if ( !(TMath::Abs(fMCEvent->Particle(particle->GetMother(0))->GetPdgCode()) == 11 && fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 22) )
                fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),4.,fWeightJetJetMC);
            }
          } else {
            fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),4.,fWeightJetJetMC);
          }
        }
      }
  
      if(fDoMesonAnalysis){
        if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedMC(particle,fMCEvent,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
          TParticle* daughter0  = (TParticle*)fMCEvent->Particle(particle->GetFirstDaughter());
          TParticle* daughter1  = (TParticle*)fMCEvent->Particle(particle->GetLastDaughter());
          Int_t pdgCode = -1;
          if(particle->GetFirstMother() > -1) pdgCode = ((TParticle*)fMCEvent->Particle( particle->GetFirstMother() ))->GetPdgCode();
          
          if(particle->GetPdgCode() == 111){
            Int_t source = GetSourceClassification(111,pdgCode);
            fHistoMCSecPi0PtvsSource[fiCut]->Fill(particle->Pt(),source,fWeightJetJetMC); 
            fHistoMCSecPi0Source[fiCut]->Fill(pdgCode);
          } else if(particle->GetPdgCode() == 221){
            fHistoMCSecEtaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); 
            fHistoMCSecEtaSource[fiCut]->Fill(pdgCode);
          }
          
          // check if photons where within acceptance
          if( ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter0,fMCEvent) &&
                ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter1,fMCEvent)){
            if(particle->GetPdgCode() == 111){
              Int_t source = GetSourceClassification(111,pdgCode);
              fHistoMCSecPi0InAccPtvsSource[fiCut]->Fill(particle->Pt(),source,fWeightJetJetMC); 
            }  
          }    
        }
      }
    
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskGammaCalo::CalculatePi0Candidates(){
  
  // Conversion Gammas
  if(fClusterCandidates->GetEntries()>0){

    for(Int_t firstGammaIndex=0;firstGammaIndex<fClusterCandidates->GetEntries();firstGammaIndex++){
      AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(firstGammaIndex));
      if (gamma0==NULL) continue;
      if ( fDoInOutTimingCluster ){
        Double_t tof = fInputEvent->GetCaloCluster(gamma0->GetCaloClusterRef())->GetTOF();
        if ( tof < fMinTimingCluster || tof > fMaxTimingCluster ) continue;
      }
      
      for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fClusterCandidates->GetEntries();secondGammaIndex++){
        AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(secondGammaIndex));
        if (gamma1==NULL) continue;
        if ( fDoInOutTimingCluster ){
          Double_t tof = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef())->GetTOF();
          if ( tof > fMinTimingCluster && tof < fMaxTimingCluster ) continue;
        }
        
        AliAODConversionMother *pi0cand = new AliAODConversionMother(gamma0,gamma1);
        pi0cand->SetLabels(firstGammaIndex,secondGammaIndex);
                
        if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),gamma0->GetLeadingCellID(),gamma1->GetLeadingCellID()))){
          if(fLocalDebugFlag == 1) DebugMethodPrint1(pi0cand,gamma0,gamma1);
          fHistoMotherInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(), fWeightJetJetMC);
          // fill new histograms
          if(!fDoLightOutput && TMath::Abs(pi0cand->GetAlpha())<0.1)
            fHistoMotherInvMassPtAlpha[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(), fWeightJetJetMC);
          
          if (fDoMesonQA > 0 && fDoMesonQA < 3){
            if ( pi0cand->M() > 0.05 && pi0cand->M() < 0.17){
              fHistoMotherPi0PtY[fiCut]->Fill(pi0cand->Pt(),pi0cand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), fWeightJetJetMC);
              fHistoMotherPi0PtAlpha[fiCut]->Fill(pi0cand->Pt(),TMath::Abs(pi0cand->GetAlpha()), fWeightJetJetMC);
              fHistoMotherPi0PtOpenAngle[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetOpeningAngle(), fWeightJetJetMC);
            }
            if ( pi0cand->M() > 0.45 && pi0cand->M() < 0.65){
              fHistoMotherEtaPtY[fiCut]->Fill(pi0cand->Pt(),pi0cand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), fWeightJetJetMC);
              fHistoMotherEtaPtAlpha[fiCut]->Fill(pi0cand->Pt(),TMath::Abs(pi0cand->GetAlpha()), fWeightJetJetMC);
              fHistoMotherEtaPtOpenAngle[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetOpeningAngle(),  fWeightJetJetMC);
            }
          }
          if(fDoTHnSparse && ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoBGCalculation()){
            Int_t zbin = 0;
            Int_t mbin = 0;
            
            if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->BackgroundHandlerType() == 0){
              zbin = fBGHandler[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
              if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
                mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
              } else {
                mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fClusterCandidates->GetEntries());
              }
            } 
            Double_t sparesFill[4] = {pi0cand->M(),pi0cand->Pt(),(Double_t)zbin,(Double_t)mbin};
            fSparseMotherInvMassPtZM[fiCut]->Fill(sparesFill,1);
          }

          if(fDoMesonQA == 4  && fIsMC == 0 && (pi0cand->Pt() > 13.) ){
            Int_t zbin = 0;
            Int_t mbin = 0;
            if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->BackgroundHandlerType() == 0){
              zbin = fBGHandler[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
              if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
                mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
              } else {
                mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fClusterCandidates->GetEntries());
              }
            }
            fInvMassTreeInvMass = pi0cand->M();
            fInvMassTreePt = pi0cand->Pt();
            fInvMassTreeAlpha = TMath::Abs(pi0cand->GetAlpha());
            fInvMassTreeTheta = pi0cand->GetOpeningAngle();
            fInvMassTreeMixPool = zbin*100 + mbin;
            fInvMassTreeZVertex = fInputEvent->GetPrimaryVertex()->GetZ();
            fInvMassTreeEta = pi0cand->Eta();
            tSigInvMassPtAlphaTheta[fiCut]->Fill();
          }
        
          if(fIsMC> 0){
            if(fInputEvent->IsA()==AliESDEvent::Class())
              ProcessTrueMesonCandidates(pi0cand,gamma0,gamma1);
            if(fInputEvent->IsA()==AliAODEvent::Class())
              ProcessTrueMesonCandidatesAOD(pi0cand,gamma0,gamma1);
          }          

          if((pi0cand->GetOpeningAngle() < 0.017) && (pi0cand->Pt() > 15.) && fDoClusterQA > 0){
            fCloseHighPtClusters = new TObjString(Form("%s",((TString)fV0Reader->GetCurrentFileName()).Data()));
            if (tBrokenFiles) tBrokenFiles->Fill();
            delete fCloseHighPtClusters;
          }
        }        
        delete pi0cand;
        pi0cand=0x0;
      }
    }
  }
}
//______________________________________________________________________
void AliAnalysisTaskGammaCalo::ProcessTrueMesonCandidates(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1)
{
  // Process True Mesons
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();
  fRconv = -1;
  
  Bool_t isTruePi0         = kFALSE;
  Bool_t isTrueEta         = kFALSE;
  Bool_t isTrueGamma        = kFALSE;
  Bool_t isSameConvertedGamma   = kFALSE;
  Int_t convertedPhotonLabel0    = -1;
  Int_t convertedPhotonLabel1    = -1;
  
  Int_t gamma0MCLabel = TrueGammaCandidate0->GetCaloPhotonMCLabel(0);   // get most probable MC label
  Int_t gamma0MotherLabel = -1;

  TParticle * gammaMC0 = 0x0;
  if(gamma0MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
    gammaMC0 = (TParticle*)fMCEvent->Particle(gamma0MCLabel);
    if (TrueGammaCandidate0->IsLargestComponentPhoton() || TrueGammaCandidate0->IsLargestComponentElectron()){    // largest component is electro magnetic
      // get mother of interest (pi0 or eta)
      if (TrueGammaCandidate0->IsLargestComponentPhoton()){                            // for photons its the direct mother 
        gamma0MotherLabel=gammaMC0->GetMother(0);
      } else if (TrueGammaCandidate0->IsLargestComponentElectron()){                         // for electrons its either the direct mother or for conversions the grandmother
        if (TrueGammaCandidate0->IsConversion() && (gammaMC0->GetMother(0) > -1)){
          convertedPhotonLabel0 = gammaMC0->GetMother(0);
          gamma0MotherLabel=fMCEvent->Particle(gammaMC0->GetMother(0))->GetMother(0);
        } else {
          gamma0MotherLabel=gammaMC0->GetMother(0); 
        }
      }
    }   
  }
  if (!TrueGammaCandidate1->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
  
  Int_t gamma1MCLabel = TrueGammaCandidate1->GetCaloPhotonMCLabel(0);   // get most probable MC label
  Int_t gamma1MotherLabel = -1;
  // check if 

  TParticle * gammaMC1 = 0x0;
  if(gamma1MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
    // Daughters Gamma 1
    gammaMC1 = (TParticle*)fMCEvent->Particle(gamma1MCLabel);
    if (TrueGammaCandidate1->IsLargestComponentPhoton() || TrueGammaCandidate1->IsLargestComponentElectron()){    // largest component is electro magnetic
      // get mother of interest (pi0 or eta)
      if (TrueGammaCandidate1->IsLargestComponentPhoton()){                            // for photons its the direct mother 
        gamma1MotherLabel      = gammaMC1->GetMother(0);
        fRconv             = gammaMC1->R();
      } else if (TrueGammaCandidate1->IsLargestComponentElectron()){                         // for electrons its either the direct mother or for conversions the grandmother
        if (TrueGammaCandidate1->IsConversion()){
          convertedPhotonLabel1   = gammaMC1->GetMother(0);
          fRconv           = gammaMC1->R();
          if(convertedPhotonLabel1 > -1) gamma1MotherLabel    = fMCEvent->Particle(convertedPhotonLabel1)->GetMother(0);
        } else {
          gamma1MotherLabel=gammaMC1->GetMother(0);
          fRconv           = gammaMC1->R();
        }
      }
    }
  }
  
  if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
    if(((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->GetPdgCode() == 111){
      isTruePi0=kTRUE;
    }
    if(((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->GetPdgCode() == 221){
      isTrueEta=kTRUE;
    }
  }

  if (fDoMesonQA == 3){
    // calculate meson properties
    fPt       = Pi0Candidate->Pt();
    fOpenRPrim   = Pi0Candidate->GetOpeningAngle();
    fInvMass   = Pi0Candidate->M();
    
    Double_t vertex[3] = {0};
    InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
    
//     cout << vertex[0] << "\t" << vertex[1] << "\t" << vertex[2] << "\t" << Pi0Candidate->Px() << "\t" << Pi0Candidate->Py()  << "\t" << Pi0Candidate->Pz() << "\t" 
//     << Pi0Candidate->Phi() << endl;  
    
    Double_t scaling = (375 - TMath::Sqrt(vertex[0]*vertex[0]+vertex[1]*vertex[1]))/(TMath::Sqrt(Pi0Candidate->Px()*Pi0Candidate->Px()+Pi0Candidate->Py()*Pi0Candidate->Py()));
    vertex[0] = vertex[0] + scaling*Pi0Candidate->Px();
    vertex[1] = vertex[1] + scaling*Pi0Candidate->Py();
    vertex[2] = vertex[2] + scaling*Pi0Candidate->Pz();
    
//     cout << vertex[0] << "\t" << vertex[1] << "\t" << vertex[2] << endl;
    
    AliVCluster* clus1 = NULL;
    clus1 = fInputEvent->GetCaloCluster(TrueGammaCandidate0->GetCaloClusterRef());
    if (!clus1) return;
    TLorentzVector clusterVector1;
    clus1->GetMomentum(clusterVector1,vertex);

    TLorentzVector* tmpvec1 = new TLorentzVector();
    tmpvec1->SetPxPyPzE(clusterVector1.Px(),clusterVector1.Py(),clusterVector1.Pz(),clusterVector1.E());
    // convert to AODConversionPhoton
    AliAODConversionPhoton *PhotonCandidate1=new AliAODConversionPhoton(tmpvec1);
    if(!PhotonCandidate1) return;

    
    AliVCluster* clus2 = NULL;
    clus2 = fInputEvent->GetCaloCluster(TrueGammaCandidate1->GetCaloClusterRef());
    if (!clus2) return;
    TLorentzVector clusterVector2;
    clus2->GetMomentum(clusterVector2,vertex);
    TLorentzVector* tmpvec2 = new TLorentzVector();
    tmpvec2->SetPxPyPzE(clusterVector2.Px(),clusterVector2.Py(),clusterVector2.Pz(),clusterVector2.E());
    // convert to AODConversionPhoton
    AliAODConversionPhoton *PhotonCandidate2=new AliAODConversionPhoton(tmpvec2);
    if(!PhotonCandidate2) return;

    AliAODConversionMother *pi0cand2 = new AliAODConversionMother(PhotonCandidate1,PhotonCandidate2);
    fInvMassRTOF = pi0cand2->M();

    delete tmpvec1;
    delete tmpvec2;
    delete PhotonCandidate1;
    delete PhotonCandidate2;
    delete pi0cand2;
    
//     cout << fOpenRPrim << "\t" << fInvMassRTOF << endl;
    
  }
  
  if (convertedPhotonLabel0 > -1 && convertedPhotonLabel1 > -1){
    if (convertedPhotonLabel0==convertedPhotonLabel1){
      isSameConvertedGamma = kTRUE;
      if (fDoMesonQA == 3 ){
        iFlag     = 0;
        tTrueInvMassROpenABPtFlag[fiCut]->Fill();
      }
    }
  }
  
  if(isTruePi0 || isTrueEta){// True Pion or Eta
    if (isTruePi0){
      fHistoTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),  fWeightJetJetMC);
//      fHistoTruePi0NonLinearity[fiCut]->Fill(TrueGammaCandidate0->E(),gammaMC0->Energy()/TrueGammaCandidate0->E());
//      fHistoTruePi0NonLinearity[fiCut]->Fill(TrueGammaCandidate1->E(),gammaMC1->Energy()/TrueGammaCandidate1->E());
      if (!fDoLightOutput && TMath::Abs(Pi0Candidate->GetAlpha())< 0.1){
        fHistoTruePi0InvMassPtAlpha[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),  fWeightJetJetMC);
        if (TrueGammaCandidate1->IsLargestComponentPhoton() && TrueGammaCandidate0->IsLargestComponentPhoton())
          fHistoTruePi0PureGammaInvMassPtAlpha[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),  fWeightJetJetMC);
      }
    }
    if (isTrueEta){
      fHistoTrueEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),  fWeightJetJetMC);
//      fHistoTrueEtaNonLinearity[fiCut]->Fill(TrueGammaCandidate0->E(),gammaMC0->Energy()/TrueGammaCandidate0->E());
//      fHistoTrueEtaNonLinearity[fiCut]->Fill(TrueGammaCandidate1->E(),gammaMC1->Energy()/TrueGammaCandidate1->E());
    }
    if (fDoMesonQA > 0 && fDoMesonQA < 3 && fIsMC < 2){
//       if(isTruePi0) {
//         fHistoPi0OpenAngleInvMassPt[fiCut]->Fill(Pi0Candidate->GetOpeningAngle(),Pi0Candidate->M(),Pi0Candidate->Pt(),1);
//       }
//       if(isTrueEta) {
//         fHistoEtaOpenAngleInvMassPt[fiCut]->Fill(Pi0Candidate->GetOpeningAngle(),Pi0Candidate->M(),Pi0Candidate->Pt(),1);
//       }
      // both gammas are real gammas
      if (TrueGammaCandidate0->IsLargestComponentPhoton() && TrueGammaCandidate1->IsLargestComponentPhoton()) {
        if (isTruePi0) fHistoTruePi0CaloPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        if (isTrueEta) fHistoTrueEtaCaloPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
      // both particles are electrons
      if (TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate1->IsLargestComponentElectron() ) {
        if (isTruePi0) fHistoTruePi0CaloElectronInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        if (isTrueEta) fHistoTrueEtaCaloElectronInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
      // both particles are converted electrons
      if ((TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate0->IsConversion()) && (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion()) ){
        if (isTruePi0 )fHistoTruePi0CaloConvertedPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        if (isTrueEta )fHistoTrueEtaCaloConvertedPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
      // 1 gamma is converted the other one is normals
      if ( (TrueGammaCandidate0->IsLargestComponentPhoton() && TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion()) ||
        (TrueGammaCandidate1->IsLargestComponentPhoton() && TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate0->IsConversion())
      ) {
        if (isTruePi0) fHistoTruePi0CaloMixedPhotonConvPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        if (isTrueEta) fHistoTrueEtaCaloMixedPhotonConvPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
      
      if ( (TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate0->IsConversion() && TrueGammaCandidate1->IsLargestComponentPhoton() && !TrueGammaCandidate1->IsMerged()) ||
        (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion() && TrueGammaCandidate0->IsLargestComponentPhoton() && !TrueGammaCandidate0->IsMerged()) 
      ) {
        if (isTruePi0) fHistoTruePi0NonMergedElectronPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }

      if ( (TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate0->IsConversion() && TrueGammaCandidate1->IsLargestComponentPhoton() && TrueGammaCandidate1->IsMerged()) ||
        (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion() && TrueGammaCandidate0->IsLargestComponentPhoton() && TrueGammaCandidate0->IsMerged()) 
      ) {
        if (isTruePi0) fHistoTruePi0NonMergedElectronMergedPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
      
      // at least one of the photon is merged
      if (TrueGammaCandidate0->IsMerged() || TrueGammaCandidate0->IsMergedPartConv() || TrueGammaCandidate0->IsDalitzMerged() || TrueGammaCandidate1->IsMerged() || TrueGammaCandidate1->IsMergedPartConv() || TrueGammaCandidate1->IsDalitzMerged() ){
        if (isTruePi0) fHistoTruePi0CaloMergedClusterInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        if (isTrueEta) fHistoTruePi0CaloMergedClusterInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
      // at least one of the photon is merged and part conv
      if (TrueGammaCandidate1->IsMergedPartConv() || TrueGammaCandidate0->IsMergedPartConv()) {  
        if (isTruePi0) fHistoTruePi0CaloMergedClusterPartConvInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        if (isTrueEta) fHistoTrueEtaCaloMergedClusterPartConvInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
    }
  
    if (fDoMesonQA == 2 && fIsMC < 2){
      // category 1: 2 real photons unmerged
      if (TrueGammaCandidate0->IsLargestComponentPhoton() && !TrueGammaCandidate0->IsMerged() && TrueGammaCandidate1->IsLargestComponentPhoton() && !TrueGammaCandidate1->IsMerged()) {
        if (isTruePi0) fHistoTruePi0Category1[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        if (isTrueEta) fHistoTrueEtaCategory1[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
      // category 2, 3: 1 real photons unmerged,  1 electron (category 2 merged, category 3 unmerged )
      // -> photon 0 is unconverted
      if ( (TrueGammaCandidate0->IsLargestComponentPhoton() && !TrueGammaCandidate0->IsMerged()) && (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion())) {
        if (isTruePi0){
          if (TrueGammaCandidate1->IsMergedPartConv())  fHistoTruePi0Category2[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
          if (!TrueGammaCandidate1->IsMergedPartConv()){
            fHistoTruePi0Category3[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
          }
        }
        if (isTrueEta){
          if (TrueGammaCandidate1->IsMergedPartConv())  fHistoTrueEtaCategory2[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
          if (!TrueGammaCandidate1->IsMergedPartConv())  fHistoTrueEtaCategory3[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        }
      }
      // -> photon 1 is unconverted
      if ( ( TrueGammaCandidate1->IsLargestComponentPhoton() && !TrueGammaCandidate1->IsMerged()) && (TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate0->IsConversion())) {
        if (isTruePi0){
          if (TrueGammaCandidate0->IsMergedPartConv())  fHistoTruePi0Category2[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
          if (!TrueGammaCandidate0->IsMergedPartConv()){
            fHistoTruePi0Category3[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
          }
        }
        if (isTrueEta){
          if (TrueGammaCandidate0->IsMergedPartConv())  fHistoTrueEtaCategory2[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
          if (!TrueGammaCandidate0->IsMergedPartConv())  fHistoTrueEtaCategory3[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        }
      }
      
      // category 4 & 6, 5, 7, 8
      if ( (TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate0->IsConversion()) && (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion()) ){
        if (isTruePi0){
          // category 4: both electrons are from same conversion 
          if (isSameConvertedGamma && !TrueGammaCandidate0->IsMergedPartConv() && !TrueGammaCandidate1->IsMergedPartConv() ) fHistoTruePi0Category4_6[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
          if (!isSameConvertedGamma ){
            if (!TrueGammaCandidate0->IsMergedPartConv() && !TrueGammaCandidate1->IsMergedPartConv()){     // category 5: both electrons from different converted photons, electrons not merged
              fHistoTruePi0Category5[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
            } else if (TrueGammaCandidate0->IsMergedPartConv() && TrueGammaCandidate1->IsMergedPartConv()){ // category 8: both electrons from different converted photons, both electrons merged    
              fHistoTruePi0Category8[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());  
            } else {                                     // category 7: both electrons from different converted photons, 1 electrons not merged, 1 electron merged    
              fHistoTruePi0Category7[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
            }
          }
        }
        if (isTrueEta){
          // category 4: both electrons are from same conversion 
          if (isSameConvertedGamma && !TrueGammaCandidate0->IsMergedPartConv() && !TrueGammaCandidate1->IsMergedPartConv()) fHistoTrueEtaCategory4_6[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
          if (!isSameConvertedGamma ){
            if (!TrueGammaCandidate0->IsMergedPartConv() && !TrueGammaCandidate1->IsMergedPartConv()){     // category 5: both electrons from different converted photons, electrons not merged
              fHistoTrueEtaCategory5[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
            } else if (TrueGammaCandidate0->IsMergedPartConv() && TrueGammaCandidate1->IsMergedPartConv()){ // category 8: both electrons from different converted photons, both electrons merged    
              fHistoTrueEtaCategory8[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());  
            } else {                                     // category 7: both electrons from different converted photons, 1 electrons not merged, 1 electron merged    
              fHistoTrueEtaCategory7[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
            }
          }
        }
      }
    }
  
    if (fDoMesonQA > 0 && fDoMesonQA < 3){
      if (isTruePi0){
        if ( Pi0Candidate->M() > 0.05 && Pi0Candidate->M() < 0.17){
          fHistoTruePi0PtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), fWeightJetJetMC);
          fHistoTruePi0PtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),TMath::Abs(Pi0Candidate->GetAlpha()), fWeightJetJetMC);
          fHistoTruePi0PtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle(), fWeightJetJetMC);
        }
      } else if (isTrueEta){
        if ( Pi0Candidate->M() > 0.45 && Pi0Candidate->M() < 0.65){
          fHistoTrueEtaPtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), fWeightJetJetMC);
          fHistoTrueEtaPtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),TMath::Abs(Pi0Candidate->GetAlpha()), fWeightJetJetMC);
          fHistoTrueEtaPtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle(), fWeightJetJetMC);
        }
      }
    }
    Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, gamma0MotherLabel, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
    if(!isPrimary){ // Secondary Meson
      // filling secondary histograms
      Long_t secMotherLabel = -1;
      if(gamma0MotherLabel > -1) secMotherLabel = ((TParticle*)fMCEvent->Particle(gamma0MotherLabel))->GetMother(0);

      Float_t weightedSec= 1;
      if((secMotherLabel > -1) && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(secMotherLabel, fMCEvent, fInputEvent) && fMCEvent->Particle(secMotherLabel)->GetPdgCode()==310){
        weightedSec= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(secMotherLabel, fMCEvent, fInputEvent)/2.; //invariant mass is additive thus the weight for the daughters has to be devide by two for the K0s at a certain pt
        //cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
      }
      if (isTruePi0) fHistoTrueSecondaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* fWeightJetJetMC);
      if (secMotherLabel >-1){
        if(fMCEvent->Particle(secMotherLabel)->GetPdgCode()==310 && isTruePi0){
          fHistoTrueSecondaryPi0FromK0sInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* fWeightJetJetMC);
          if (fDoMesonQA > 0 && fDoMesonQA < 3 && fIsMC < 2)fHistoTrueK0sWithPi0DaughterMCPt[fiCut]->Fill(fMCEvent->Particle(secMotherLabel)->Pt(), fWeightJetJetMC);
          if (fDoMesonQA == 3 ){
            iFlag     = 3;
            if (!isSameConvertedGamma)tTrueInvMassROpenABPtFlag[fiCut]->Fill();
          }
        } else if(fMCEvent->Particle(secMotherLabel)->GetPdgCode()==130 && isTruePi0){
          fHistoTrueSecondaryPi0FromK0lInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* fWeightJetJetMC);
          if (fDoMesonQA > 0 && fDoMesonQA < 3 && fIsMC < 2)fHistoTrueK0lWithPi0DaughterMCPt[fiCut]->Fill(fMCEvent->Particle(secMotherLabel)->Pt(), fWeightJetJetMC);

        } else if(fMCEvent->Particle(secMotherLabel)->GetPdgCode()==221 && isTruePi0){
          fHistoTrueSecondaryPi0FromEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* fWeightJetJetMC);
          if (fDoMesonQA > 0 && fDoMesonQA < 3 && fIsMC < 2)fHistoTrueEtaWithPi0DaughterMCPt[fiCut]->Fill(fMCEvent->Particle(secMotherLabel)->Pt(), fWeightJetJetMC);
          if (fDoMesonQA == 3 ){
            iFlag     = 4;
            if (!isSameConvertedGamma)tTrueInvMassROpenABPtFlag[fiCut]->Fill();
          }
        } else if(fMCEvent->Particle(secMotherLabel)->GetPdgCode()==3122 && isTruePi0){
          fHistoTrueSecondaryPi0FromLambdaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* fWeightJetJetMC);
          if (fDoMesonQA > 0 && fDoMesonQA < 3 && fIsMC < 2)fHistoTrueLambdaWithPi0DaughterMCPt[fiCut]->Fill(fMCEvent->Particle(secMotherLabel)->Pt(), fWeightJetJetMC);
          if (fDoMesonQA == 3 ){
            iFlag     = 5;
            if (!isSameConvertedGamma)tTrueInvMassROpenABPtFlag[fiCut]->Fill();
          }
        } else if (isTruePi0){
          if (fDoMesonQA == 3 ){
            iFlag     = 6;
            if (!isSameConvertedGamma)tTrueInvMassROpenABPtFlag[fiCut]->Fill();
          }
        } else if (isTrueEta){
          if (fDoMesonQA == 3 ){
            iFlag     = 7;
            if (!isSameConvertedGamma)tTrueInvMassROpenABPtFlag[fiCut]->Fill();
          }
        }
      }
    } else { // Only primary pi0 for efficiency calculation
      // filling primary histograms
      Float_t weighted= 1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma1MotherLabel, fMCEvent, fInputEvent)){
        if (((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt()>0.005){
          weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(gamma1MotherLabel, fMCEvent, fInputEvent);
          //                      cout << "rec \t " <<gamma1MotherLabel << "\t" <<  weighted << endl;
        }
      }
      if (isTruePi0){
        if (fDoMesonQA == 3 ){
          iFlag     = 1;
          if (!isSameConvertedGamma)tTrueInvMassROpenABPtFlag[fiCut]->Fill();
        }

        fHistoTruePrimaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* fWeightJetJetMC);
        fHistoTruePrimaryPi0W0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
        fProfileTruePrimaryPi0WeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* fWeightJetJetMC);
        if (CheckVectorForDoubleCount(fVectorDoubleCountTruePi0s,gamma0MotherLabel)) fHistoDoubleCountTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), weighted*fWeightJetJetMC);
      } else if (isTrueEta){
        if (fDoMesonQA == 3 ){
          iFlag     = 2;
          if (!isSameConvertedGamma)tTrueInvMassROpenABPtFlag[fiCut]->Fill();
        }

        fHistoTruePrimaryEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* fWeightJetJetMC);
        fHistoTruePrimaryEtaW0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
        fProfileTruePrimaryEtaWeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* fWeightJetJetMC);
        if (CheckVectorForDoubleCount(fVectorDoubleCountTrueEtas,gamma0MotherLabel)) fHistoDoubleCountTrueEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), weighted*fWeightJetJetMC);
      }
        
      if (fDoMesonQA > 0 && fDoMesonQA < 3 && fIsMC<2){
        if(isTruePi0){ // Only primary pi0 for resolution
          fHistoTruePrimaryPi0MCPtResolPt[fiCut]->Fill(((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt(),(Pi0Candidate->Pt()-((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt())/((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt(),weighted* fWeightJetJetMC);
        }
        if (isTrueEta){ // Only primary eta for resolution
          fHistoTruePrimaryEtaMCPtResolPt[fiCut]->Fill(((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt(),(Pi0Candidate->Pt()-((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt())/((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt(),weighted* fWeightJetJetMC);
        }
      }
    }
  } else if(!isTruePi0 && !isTrueEta){ // Background
    if (fDoMesonQA > 0 && fDoMesonQA < 3){
      if(gamma0MotherLabel>-1 && gamma1MotherLabel>-1){ // Both Tracks are Photons and have a mother but not Pi0 or Eta
        fHistoTrueBckGGInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);

        if( ((((TParticle*)fMCEvent->Particle(gamma0MotherLabel))->GetPdgCode() == 111
            || ((TParticle*)fMCEvent->Particle(gamma0MotherLabel))->GetPdgCode() == 221)
            && (TrueGammaCandidate0->IsMerged() || TrueGammaCandidate0->IsMergedPartConv()))
            ||
            ((((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->GetPdgCode() == 111
            || ((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->GetPdgCode() == 221)
            && (TrueGammaCandidate1->IsMerged() || TrueGammaCandidate1->IsMergedPartConv()))
        ){
          fHistoTrueBckFullMesonContainedInOneClusterInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
        }else if( (TrueGammaCandidate0->E()/Pi0Candidate->E() > 0.7) || (TrueGammaCandidate1->E()/Pi0Candidate->E() > 0.7) ){
          fHistoTrueBckAsymEClustersInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
        }
      } else { // No photon or without mother
        fHistoTrueBckContInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
      }
    }
  }

}
//______________________________________________________________________
void AliAnalysisTaskGammaCalo::ProcessTrueMesonCandidatesAOD(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1)
{
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  // Process True Mesons
  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (AODMCTrackArray == NULL) return;
  
  Bool_t isTruePi0         = kFALSE;
  Bool_t isTrueEta         = kFALSE;
  Bool_t isSameConvertedGamma   = kFALSE;
  Int_t convertedPhotonLabel0    = -1;
  Int_t convertedPhotonLabel1    = -1;
    
  Int_t gamma0MCLabel       = TrueGammaCandidate0->GetCaloPhotonMCLabel(0);   // get most probable MC label
  Int_t gamma0MotherLabel     = -1;
  
  // check if 
  AliAODMCParticle * gammaMC0 = 0x0;
  if(gamma0MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
    // Daughters Gamma 0
    gammaMC0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MCLabel));
    if (TrueGammaCandidate0->IsLargestComponentPhoton() || TrueGammaCandidate0->IsLargestComponentElectron()){    // largest component is electro magnetic
      // get mother of interest (pi0 or eta)
      if (TrueGammaCandidate0->IsLargestComponentPhoton()){                            // for photons its the direct mother 
        gamma0MotherLabel=gammaMC0->GetMother();
      } else if (TrueGammaCandidate0->IsLargestComponentElectron()){                         // for electrons its either the direct mother or for conversions the grandmother
        if (TrueGammaCandidate0->IsConversion()){
          convertedPhotonLabel0 = gammaMC0->GetMother();
          AliAODMCParticle * gammaGrandMotherMC0 =  static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gammaMC0->GetMother()));
          gamma0MotherLabel=gammaGrandMotherMC0->GetMother();
        } else gamma0MotherLabel=gammaMC0->GetMother(); 
      }
    }
  }

  Int_t gamma1MCLabel         = TrueGammaCandidate1->GetCaloPhotonMCLabel(0);   // get most probable MC label
  Int_t gamma1MotherLabel     = -1;
  
  // check if 
  AliAODMCParticle *gammaMC1  = 0x0;
  if(gamma1MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
    // Daughters Gamma 1
    gammaMC1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MCLabel));
    if (TrueGammaCandidate1->IsLargestComponentPhoton() || TrueGammaCandidate1->IsLargestComponentElectron()){    // largest component is electro magnetic
      // get mother of interest (pi0 or eta)
      if (TrueGammaCandidate1->IsLargestComponentPhoton()){                            // for photons its the direct mother 
        gamma1MotherLabel=gammaMC1->GetMother();
      } else if (TrueGammaCandidate1->IsLargestComponentElectron()){                         // for electrons its either the direct mother or for conversions the grandmother
        if (TrueGammaCandidate1->IsConversion()){
          convertedPhotonLabel1 = gammaMC1->GetMother();
          AliAODMCParticle * gammaGrandMotherMC1 =  static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gammaMC1->GetMother()));
          gamma1MotherLabel=gammaGrandMotherMC1->GetMother();
        } else gamma1MotherLabel=gammaMC1->GetMother(); 
      }
    }
  }
      
  if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
    if(((AliAODMCParticle*)AODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 111){
      isTruePi0=kTRUE;
    }
    if(((AliAODMCParticle*)AODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 221){
      isTrueEta=kTRUE;
    }
  }

  if (convertedPhotonLabel0 > -1 && convertedPhotonLabel1 > 1){
    if (convertedPhotonLabel0==convertedPhotonLabel1){
      isSameConvertedGamma = kTRUE;
      if (fDoMesonQA == 3 ){
//         fHistoGammaOpenAngleInvMassPt[fiCut]->Fill(Pi0Candidate->GetOpeningAngle(),Pi0Candidate->M(),Pi0Candidate->Pt(),1);
      }
    }
  }

  
  if(isTruePi0 || isTrueEta){// True Pion or Eta
    if(isTruePi0){
      fHistoTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
//      fHistoTruePi0NonLinearity[fiCut]->Fill(TrueGammaCandidate0->E(),gammaMC0->E()/TrueGammaCandidate0->E());
//      fHistoTruePi0NonLinearity[fiCut]->Fill(TrueGammaCandidate1->E(),gammaMC1->E()/TrueGammaCandidate1->E());
      if (!fDoLightOutput && TMath::Abs(Pi0Candidate->GetAlpha())< 0.1){
        fHistoTruePi0InvMassPtAlpha[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),  fWeightJetJetMC);
        if (TrueGammaCandidate1->IsLargestComponentPhoton() && TrueGammaCandidate0->IsLargestComponentPhoton())
          fHistoTruePi0PureGammaInvMassPtAlpha[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),  fWeightJetJetMC);
      }
    }
    if(isTrueEta){
      fHistoTrueEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
//      fHistoTrueEtaNonLinearity[fiCut]->Fill(TrueGammaCandidate0->E(),gammaMC0->E()/TrueGammaCandidate0->E());
//      fHistoTrueEtaNonLinearity[fiCut]->Fill(TrueGammaCandidate1->E(),gammaMC1->E()/TrueGammaCandidate1->E());
    }
    if (fDoMesonQA > 0 && fDoMesonQA < 3 && fIsMC < 2){
//       if(isTruePi0) {
//         fHistoPi0OpenAngleInvMassPt[fiCut]->Fill(Pi0Candidate->GetOpeningAngle(),Pi0Candidate->M(),Pi0Candidate->Pt(),1);
//       }
//       if(isTrueEta) {
//         fHistoEtaOpenAngleInvMassPt[fiCut]->Fill(Pi0Candidate->GetOpeningAngle(),Pi0Candidate->M(),Pi0Candidate->Pt(),1);
//       }
      // both gammas are real gammas
      if (TrueGammaCandidate0->IsLargestComponentPhoton() && TrueGammaCandidate1->IsLargestComponentPhoton()) {
        if (isTruePi0)fHistoTruePi0CaloPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        if (isTrueEta)fHistoTrueEtaCaloPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
      // both particles are electrons
      if (TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate1->IsLargestComponentElectron() ) {
        if (isTruePi0) fHistoTruePi0CaloElectronInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        if (isTrueEta) fHistoTrueEtaCaloElectronInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
      // both particles are converted electrons
      if ((TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate0->IsConversion()) && (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion()) ){
        if (isTruePi0 )fHistoTruePi0CaloConvertedPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        if (isTrueEta )fHistoTrueEtaCaloConvertedPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
      // 1 gamma is converted the other one is normals
      if ( (TrueGammaCandidate0->IsLargestComponentPhoton() && TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion()) ||
        (TrueGammaCandidate1->IsLargestComponentPhoton() && TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate0->IsConversion())
      ) {
        if (isTruePi0) fHistoTruePi0CaloMixedPhotonConvPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        if (isTrueEta) fHistoTrueEtaCaloMixedPhotonConvPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
      
      if ( (TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate0->IsConversion() && TrueGammaCandidate1->IsLargestComponentPhoton() && !TrueGammaCandidate1->IsMerged()) ||
        (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion() && TrueGammaCandidate0->IsLargestComponentPhoton() && !TrueGammaCandidate0->IsMerged()) 
      ) {
        if (isTruePi0) fHistoTruePi0NonMergedElectronPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }

      if ( (TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate0->IsConversion() && TrueGammaCandidate1->IsLargestComponentPhoton() && TrueGammaCandidate1->IsMerged()) ||
        (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion() && TrueGammaCandidate0->IsLargestComponentPhoton() && TrueGammaCandidate0->IsMerged()) 
      ) {
        if (isTruePi0) fHistoTruePi0NonMergedElectronMergedPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }

      // at least one of the photon is merged
      if (TrueGammaCandidate0->IsMerged() || TrueGammaCandidate0->IsMergedPartConv() || TrueGammaCandidate0->IsDalitzMerged() || TrueGammaCandidate1->IsMerged() || TrueGammaCandidate1->IsMergedPartConv() || TrueGammaCandidate1->IsDalitzMerged() ){
        if (isTruePi0) fHistoTruePi0CaloMergedClusterInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        if (isTrueEta) fHistoTrueEtaCaloMergedClusterInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
      // at least one of the photon is merged and part conv
      if (TrueGammaCandidate1->IsMergedPartConv() || TrueGammaCandidate0->IsMergedPartConv()) {
        if (isTruePi0)fHistoTruePi0CaloMergedClusterPartConvInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        if (isTrueEta)fHistoTrueEtaCaloMergedClusterPartConvInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
    }

    if (fDoMesonQA == 2 && fIsMC < 2){
      // category 1: 2 real photons unmerged
      if (TrueGammaCandidate0->IsLargestComponentPhoton() && !TrueGammaCandidate0->IsMerged() && TrueGammaCandidate1->IsLargestComponentPhoton() && !TrueGammaCandidate1->IsMerged()) {
        if (isTruePi0) fHistoTruePi0Category1[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        if (isTrueEta) fHistoTrueEtaCategory1[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
      // category 2, 3: 1 real photons unmerged,  1 electron (category 2 merged, category 3 unmerged )
      // -> photon 0 is unconverted
      if ( (TrueGammaCandidate0->IsLargestComponentPhoton() && !TrueGammaCandidate0->IsMerged()) && (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion())) {
        if (isTruePi0){
          if (TrueGammaCandidate1->IsMergedPartConv())  fHistoTruePi0Category2[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
          if (!TrueGammaCandidate1->IsMergedPartConv())  fHistoTruePi0Category3[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        }
        if (isTrueEta){
          if (TrueGammaCandidate1->IsMergedPartConv())  fHistoTrueEtaCategory2[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
          if (!TrueGammaCandidate1->IsMergedPartConv())  fHistoTrueEtaCategory3[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        }
      }
      // -> photon 1 is unconverted
      if ( ( TrueGammaCandidate1->IsLargestComponentPhoton() && !TrueGammaCandidate1->IsMerged()) && (TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate0->IsConversion())) {
        if (isTruePi0){
          if (TrueGammaCandidate0->IsMergedPartConv())  fHistoTruePi0Category2[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
          if (!TrueGammaCandidate0->IsMergedPartConv())  fHistoTruePi0Category3[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        }
        if (isTrueEta){
          if (TrueGammaCandidate0->IsMergedPartConv())  fHistoTrueEtaCategory2[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
          if (!TrueGammaCandidate0->IsMergedPartConv())  fHistoTrueEtaCategory3[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        }
      }
      
      // category 4 & 6, 5, 7, 8
      if ( (TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate0->IsConversion()) && (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion()) ){
        if (isTruePi0){
          // category 4: both electrons are from same conversion 
          if (isSameConvertedGamma && !TrueGammaCandidate0->IsMergedPartConv() && !TrueGammaCandidate1->IsMergedPartConv() ) fHistoTruePi0Category4_6[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
          if (!isSameConvertedGamma ){
            if (!TrueGammaCandidate0->IsMergedPartConv() && !TrueGammaCandidate1->IsMerged()){     // category 5: both electrons from different converted photons, electrons not merged
              fHistoTruePi0Category5[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
            } else if (TrueGammaCandidate0->IsMergedPartConv() && TrueGammaCandidate1->IsMerged()){ // category 8: both electrons from different converted photons, both electrons merged    
              fHistoTruePi0Category8[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());  
            } else {                                     // category 7: both electrons from different converted photons, 1 electrons not merged, 1 electron merged    
              fHistoTruePi0Category7[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
            }
          }
        }
        if (isTrueEta){
          // category 4: both electrons are from same conversion 
          if (isSameConvertedGamma && !TrueGammaCandidate0->IsMergedPartConv() && !TrueGammaCandidate1->IsMergedPartConv()) fHistoTrueEtaCategory4_6[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
          if (!isSameConvertedGamma ){
            if (!TrueGammaCandidate0->IsMergedPartConv() && !TrueGammaCandidate1->IsMergedPartConv()){     // category 5: both electrons from different converted photons, electrons not merged
              fHistoTrueEtaCategory5[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
            } else if (TrueGammaCandidate0->IsMergedPartConv() && TrueGammaCandidate1->IsMergedPartConv()){ // category 8: both electrons from different converted photons, both electrons merged    
              fHistoTrueEtaCategory8[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());  
            } else {                                     // category 7: both electrons from different converted photons, 1 electrons not merged, 1 electron merged    
              fHistoTrueEtaCategory7[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
            }
          }
        }
      }
    }
    
    if (fDoMesonQA > 0 && fDoMesonQA < 3){
      if (isTruePi0){
        if ( Pi0Candidate->M() > 0.05 && Pi0Candidate->M() < 0.17){
        fHistoTruePi0PtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), fWeightJetJetMC);
        fHistoTruePi0PtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),TMath::Abs(Pi0Candidate->GetAlpha()), fWeightJetJetMC);
        fHistoTruePi0PtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle(), fWeightJetJetMC);
        }
      } else if (isTrueEta){
        if ( Pi0Candidate->M() > 0.45 && Pi0Candidate->M() < 0.65){
        fHistoTrueEtaPtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), fWeightJetJetMC);
        fHistoTrueEtaPtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),TMath::Abs(Pi0Candidate->GetAlpha()), fWeightJetJetMC);
        fHistoTrueEtaPtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle(), fWeightJetJetMC);
        }
      }
    }

    Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MotherLabel)), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
    if(!isPrimary){ // Secondary Meson
      Long_t secMotherLabel = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->GetMother();
      Float_t weightedSec= 1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(secMotherLabel, 0x0, fInputEvent) && static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->GetPdgCode()==310){
        weightedSec= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(secMotherLabel, 0x0, fInputEvent)/2.; //invariant mass is additive thus the weight for the daughters has to be devide by two for the K0s at a certain pt
        //cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
      }
      if (isTruePi0) fHistoTrueSecondaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* fWeightJetJetMC);
      if (secMotherLabel >-1){
        if(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->GetPdgCode()==310 && isTruePi0 ){
          fHistoTrueSecondaryPi0FromK0sInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* fWeightJetJetMC);
          if (fDoMesonQA > 0 && fDoMesonQA < 3 && fIsMC < 2 )fHistoTrueK0sWithPi0DaughterMCPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->Pt(), fWeightJetJetMC);
        }
        if(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->GetPdgCode()==130 && isTruePi0 ){
          fHistoTrueSecondaryPi0FromK0lInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* fWeightJetJetMC);
          if (fDoMesonQA > 0 && fDoMesonQA < 3 && fIsMC < 2 )fHistoTrueK0lWithPi0DaughterMCPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->Pt(), fWeightJetJetMC);
        }
        if(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->GetPdgCode()==221 && isTruePi0){
          fHistoTrueSecondaryPi0FromEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* fWeightJetJetMC);
          if (fDoMesonQA > 0 && fDoMesonQA < 3 && fIsMC < 2)fHistoTrueEtaWithPi0DaughterMCPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->Pt(), fWeightJetJetMC);
        }
        if(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->GetPdgCode()==3122 && isTruePi0){
          fHistoTrueSecondaryPi0FromLambdaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* fWeightJetJetMC);
          if (fDoMesonQA > 0 && fDoMesonQA < 3 && fIsMC < 2)fHistoTrueLambdaWithPi0DaughterMCPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->Pt(), fWeightJetJetMC);
        }
      }
    } else{ // Only primary pi0 for efficiency calculation
      Float_t weighted= 1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma1MotherLabel, 0x0, fInputEvent)){
        if (static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt()>0.005){
        weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(gamma1MotherLabel, 0x0, fInputEvent);
        //                      cout << "rec \t " <<gamma1MotherLabel << "\t" <<  weighted << endl;
        }
      }
      if (isTruePi0){
        fHistoTruePrimaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* fWeightJetJetMC);
        fHistoTruePrimaryPi0W0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
        fProfileTruePrimaryPi0WeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* fWeightJetJetMC);
        if (CheckVectorForDoubleCount(fVectorDoubleCountTruePi0s,gamma0MotherLabel)) fHistoDoubleCountTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), weighted*fWeightJetJetMC);
      } else if (isTrueEta){
        fHistoTruePrimaryEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* fWeightJetJetMC);
        fHistoTruePrimaryEtaW0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
        fProfileTruePrimaryEtaWeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* fWeightJetJetMC);
        if (CheckVectorForDoubleCount(fVectorDoubleCountTrueEtas,gamma0MotherLabel)) fHistoDoubleCountTrueEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), weighted*fWeightJetJetMC);
      }
      if (fDoMesonQA > 0 && fDoMesonQA < 3 && fIsMC<2){
        if(isTruePi0){ // Only primary pi0 for resolution
          fHistoTruePrimaryPi0MCPtResolPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt(),
                              (Pi0Candidate->Pt()-static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt())/static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt(),weighted* fWeightJetJetMC);
        
        }
        if (isTrueEta){ // Only primary eta for resolution
          fHistoTruePrimaryEtaMCPtResolPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt(),
                              (Pi0Candidate->Pt()-static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt())/static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt(),weighted* fWeightJetJetMC);
        }
      }
    }
  } else if(!isTruePi0 && !isTrueEta) { // Background
    if (fDoMesonQA > 0 && fDoMesonQA < 3){
      if(gamma0MotherLabel>-1 && gamma1MotherLabel>-1){ // Both Tracks are Photons and have a mother but not Pi0 or Eta
        fHistoTrueBckGGInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);

        if( ((((AliAODMCParticle*)AODMCTrackArray->At(gamma0MotherLabel))->GetPdgCode() == 111
            || ((AliAODMCParticle*)AODMCTrackArray->At(gamma0MotherLabel))->GetPdgCode() == 221)
            && (TrueGammaCandidate0->IsMerged() || TrueGammaCandidate0->IsMergedPartConv()))
            ||
            ((((AliAODMCParticle*)AODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 111
            || ((AliAODMCParticle*)AODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 221)
            && (TrueGammaCandidate1->IsMerged() || TrueGammaCandidate1->IsMergedPartConv()))
        ){
          fHistoTrueBckFullMesonContainedInOneClusterInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
        }else if( (TrueGammaCandidate0->E()/Pi0Candidate->E() > 0.7) || (TrueGammaCandidate1->E()/Pi0Candidate->E() > 0.7) ){
          fHistoTrueBckAsymEClustersInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
        }
      } else { // No photon or without mother
        fHistoTrueBckContInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskGammaCalo::CalculateBackground(){
  
  Int_t zbin= fBGHandler[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
  Int_t mbin = 0;
  
  if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
    mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
  } else {
    mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fClusterCandidates->GetEntries());
  }
  
  if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
    for(Int_t nEventsInBG=0;nEventsInBG<fBGHandler[fiCut]->GetNBGEvents();nEventsInBG++){
      AliGammaConversionAODVector *previousEventV0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,nEventsInBG);
      for(Int_t iCurrent=0;iCurrent<fClusterCandidates->GetEntries();iCurrent++){
        AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fClusterCandidates->At(iCurrent));
        for(Int_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
          AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));
          AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
          backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
        
          if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
            ->MesonIsSelected(backgroundCandidate,kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()),currentEventGoodV0.GetLeadingCellID(),previousGoodV0.GetLeadingCellID())){
            fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(), fWeightJetJetMC);
            if(fDoTHnSparse){
              Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
              fSparseMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
            }
            if(!fDoLightOutput && TMath::Abs(backgroundCandidate->GetAlpha())<0.1)
              fHistoMotherBackInvMassPtAlpha[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(), fWeightJetJetMC);

            if(fDoMesonQA == 4 && fIsMC == 0 && (backgroundCandidate->Pt() > 13.) ){
              fInvMassTreeInvMass = backgroundCandidate->M();
              fInvMassTreePt = backgroundCandidate->Pt();
              fInvMassTreeAlpha = TMath::Abs(backgroundCandidate->GetAlpha());
              fInvMassTreeTheta = backgroundCandidate->GetOpeningAngle();
              fInvMassTreeMixPool = zbin*100 + mbin;
              fInvMassTreeZVertex = fInputEvent->GetPrimaryVertex()->GetZ();
              fInvMassTreeEta = backgroundCandidate->Eta();
              tBckInvMassPtAlphaTheta[fiCut]->Fill();
            }
          }
          delete backgroundCandidate;
          backgroundCandidate = 0x0;
        }
      }
    }
  } else if ( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UsePtmaxMethod() ) {

    Double_t currentPtMax = 0; Double_t previousPtMax = 0;
    Double_t currentEta = 0;   Double_t previousEta = 0;
    Double_t currentPhi = 0;   Double_t previousPhi = 0;
    Bool_t acceptedPtMax = kFALSE;

    for(Int_t nEventsInBG=0;nEventsInBG <fBGHandler[fiCut]->GetNBGEvents();nEventsInBG++){
      AliGammaConversionAODVector *previousEventV0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,nEventsInBG);
      if(previousEventV0s){
        acceptedPtMax = kFALSE;
        currentPtMax = 0; previousPtMax = 0; currentEta = 0; previousEta = 0; currentPhi = 0; previousPhi = 0;
        for(Int_t iCurrent=0;iCurrent<fClusterCandidates->GetEntries();iCurrent++){
          AliAODConversionPhoton *currentV0 = (AliAODConversionPhoton*)(fClusterCandidates->At(iCurrent));
          if(currentV0->GetPhotonPt() > currentPtMax){ currentPtMax = currentV0->GetPhotonPt(); currentEta = currentV0->GetPhotonEta(); currentPhi = currentV0->GetPhotonPhi(); }
          for(Int_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
            AliAODConversionPhoton *previousV0 = (AliAODConversionPhoton*)(previousEventV0s->at(iPrevious));
            if(previousV0->GetPhotonPt() > previousPtMax){ previousPtMax = previousV0->GetPhotonPt(); previousEta = previousV0->GetPhotonEta(); previousPhi = previousV0->GetPhotonPhi(); }
          }
        }
        if(currentPtMax > 0 && previousPtMax > 0){
         if(TMath::Sqrt(pow((currentEta-previousEta),2)+pow((currentPhi-previousPhi),2)) < 0.2) acceptedPtMax = kTRUE;
        }
        if(acceptedPtMax){
          for(Int_t iCurrent=0;iCurrent<fClusterCandidates->GetEntries();iCurrent++){
            AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fClusterCandidates->At(iCurrent));
            for(Int_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
              AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));
              AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
              backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());

              if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
                ->MesonIsSelected(backgroundCandidate,kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()),currentEventGoodV0.GetLeadingCellID(),previousGoodV0.GetLeadingCellID())){
                fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(), fWeightJetJetMC);
              }
              delete backgroundCandidate;
              backgroundCandidate = 0x0;
            }
          }
        }
      }
    }
  } else {
    for(Int_t nEventsInBG=0;nEventsInBG <fBGHandler[fiCut]->GetNBGEvents();nEventsInBG++){
      AliGammaConversionAODVector *previousEventV0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,nEventsInBG);
      if(previousEventV0s){
        for(Int_t iCurrent=0;iCurrent<fClusterCandidates->GetEntries();iCurrent++){
          AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fClusterCandidates->At(iCurrent));
          for(Int_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
        
            AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));
            AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
            backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
        
            if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate,kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),currentEventGoodV0.GetLeadingCellID(),previousGoodV0.GetLeadingCellID()))){
              fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(), fWeightJetJetMC);
              if(fDoTHnSparse){
                Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
                fSparseMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
              }
              if(!fDoLightOutput && TMath::Abs(backgroundCandidate->GetAlpha())<0.1)
                fHistoMotherBackInvMassPtAlpha[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(), fWeightJetJetMC);
              
              if(fDoMesonQA == 4 && fIsMC == 0 && (backgroundCandidate->Pt() > 13.) ){
                fInvMassTreeInvMass = backgroundCandidate->M();
                fInvMassTreePt = backgroundCandidate->Pt();
                fInvMassTreeAlpha = TMath::Abs(backgroundCandidate->GetAlpha());
                fInvMassTreeTheta = backgroundCandidate->GetOpeningAngle();
                fInvMassTreeMixPool = zbin*100 + mbin;
                fInvMassTreeZVertex = fInputEvent->GetPrimaryVertex()->GetZ();
                fInvMassTreeEta = backgroundCandidate->Eta();
                tBckInvMassPtAlphaTheta[fiCut]->Fill();
              }
            }
            delete backgroundCandidate;
            backgroundCandidate = 0x0;
          }
        }
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskGammaCalo::RotateParticle(AliAODConversionPhoton *gamma){
  Int_t fNDegreesPMBackground= ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->NDegreesRotation();
  Double_t nRadiansPM = fNDegreesPMBackground*TMath::Pi()/180;
  Double_t rotationValue = fRandom.Rndm()*2*nRadiansPM + TMath::Pi()-nRadiansPM;
  gamma->RotateZ(rotationValue);
}

//________________________________________________________________________
void AliAnalysisTaskGammaCalo::UpdateEventByEventData(){
  //see header file for documentation
  if(fClusterCandidates->GetEntries() >0 ){
    if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
      fBGHandler[fiCut]->AddEvent(fClusterCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fV0Reader->GetNumberOfPrimaryTracks(),fEventPlaneAngle);
    } else { // means we use #V0s for multiplicity
      fBGHandler[fiCut]->AddEvent(fClusterCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fClusterCandidates->GetEntries(),fEventPlaneAngle);
    }
  }
}


//________________________________________________________________________
void AliAnalysisTaskGammaCalo::FillPhotonBackgroundHist(AliAODConversionPhoton *TruePhotonCandidate, Int_t pdgCode)
{
  // Bck = 0 e+-, 1 pi+-, 2 p+-, 3 K+-, 4 n, 5 K0s, 6 Lambda, 7 mu+-, 8 K0l, 9 rest
  if(!fDoLightOutput){
    if(fIsFromMBHeader){
      if(TMath::Abs(pdgCode) == 11)        fHistoClusPhotonBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),0., fWeightJetJetMC);
      else if( TMath::Abs(pdgCode)==211)   fHistoClusPhotonBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),1., fWeightJetJetMC);
      else if( TMath::Abs(pdgCode)==2212)  fHistoClusPhotonBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),2., fWeightJetJetMC);
      else if( TMath::Abs(pdgCode)==321)   fHistoClusPhotonBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3., fWeightJetJetMC);
      else if( TMath::Abs(pdgCode)==2112)  fHistoClusPhotonBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),4., fWeightJetJetMC);
      else if( TMath::Abs(pdgCode)==310)   fHistoClusPhotonBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),5., fWeightJetJetMC);
      else if( TMath::Abs(pdgCode)==3122)  fHistoClusPhotonBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),6., fWeightJetJetMC);
      else if( TMath::Abs(pdgCode)==13)    fHistoClusPhotonBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),7., fWeightJetJetMC);
      else if( TMath::Abs(pdgCode)==130)   fHistoClusPhotonBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),8., fWeightJetJetMC);
      else                          fHistoClusPhotonBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),9., fWeightJetJetMC);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskGammaCalo::FillPhotonPlusConversionBackgroundHist(AliAODConversionPhoton *TruePhotonCandidate, Int_t pdgCode)
{
  // Bck = 0 e+-, 1 pi+-, 2 p+-, 3 K+-, 4 n, 5 K0s, 6 Lambda, 7 mu+-, 8 K0l, 9 rest
  if(!fDoLightOutput){
    if(fIsFromMBHeader){
      if(TMath::Abs(pdgCode) == 11)        fHistoClusPhotonPlusConvBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),0., fWeightJetJetMC);
      else if( TMath::Abs(pdgCode)==211)   fHistoClusPhotonPlusConvBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),1., fWeightJetJetMC);
      else if( TMath::Abs(pdgCode)==2212)  fHistoClusPhotonPlusConvBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),2., fWeightJetJetMC);
      else if( TMath::Abs(pdgCode)==321)   fHistoClusPhotonPlusConvBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3., fWeightJetJetMC);
      else if( TMath::Abs(pdgCode)==2112)  fHistoClusPhotonPlusConvBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),4., fWeightJetJetMC);
      else if( TMath::Abs(pdgCode)==310)   fHistoClusPhotonPlusConvBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),5., fWeightJetJetMC);
      else if( TMath::Abs(pdgCode)==3122)  fHistoClusPhotonPlusConvBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),6., fWeightJetJetMC);
      else if( TMath::Abs(pdgCode)==13)    fHistoClusPhotonPlusConvBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),7., fWeightJetJetMC);
      else if( TMath::Abs(pdgCode)==130)   fHistoClusPhotonPlusConvBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),8., fWeightJetJetMC);
      else                          fHistoClusPhotonPlusConvBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),9., fWeightJetJetMC);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskGammaCalo::FillPhotonBackgroundM02Hist(AliAODConversionPhoton *TruePhotonCandidate, AliVCluster* clus, Int_t pdgCode)
{
  // Bck = e+-, pi+-, n, K+-, K0l, rest
  if(!fDoLightOutput){
    if(fIsFromMBHeader){
      if(TMath::Abs(pdgCode) == 11)        fHistoClustPhotonElectronBGPtM02[fiCut]->Fill(TruePhotonCandidate->Pt(),clus->GetM02(), fWeightJetJetMC);
      else if( TMath::Abs(pdgCode)==211)   fHistoClustPhotonPionBGPtM02[fiCut]->Fill(    TruePhotonCandidate->Pt(),clus->GetM02(), fWeightJetJetMC);
      else if( TMath::Abs(pdgCode)==2212)  fHistoClustPhotonNeutronBGPtM02[fiCut]->Fill( TruePhotonCandidate->Pt(),clus->GetM02(), fWeightJetJetMC);
      else if( TMath::Abs(pdgCode)==321)   fHistoClustPhotonKaonBGPtM02[fiCut]->Fill(    TruePhotonCandidate->Pt(),clus->GetM02(), fWeightJetJetMC);
      else if( TMath::Abs(pdgCode)==130)   fHistoClustPhotonK0lBGPtM02[fiCut]->Fill(     TruePhotonCandidate->Pt(),clus->GetM02(), fWeightJetJetMC);
      else                          fHistoClustPhotonRestBGPtM02[fiCut]->Fill(    TruePhotonCandidate->Pt(),clus->GetM02(), fWeightJetJetMC);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskGammaCalo::FillPhotonPlusConversionBackgroundM02Hist(AliAODConversionPhoton *TruePhotonCandidate, AliVCluster* clus, Int_t pdgCode)
{
  // Bck = e+-, pi+-, n, K+-, K0l, rest
  if(!fDoLightOutput){
    if(fIsFromMBHeader){
      if(TMath::Abs(pdgCode) == 11)        fHistoClustPhotonPlusConvElectronBGPtM02[fiCut]->Fill(TruePhotonCandidate->Pt(),clus->GetM02(), fWeightJetJetMC);
      else if( TMath::Abs(pdgCode)==211)   fHistoClustPhotonPlusConvPionBGPtM02[fiCut]->Fill(    TruePhotonCandidate->Pt(),clus->GetM02(), fWeightJetJetMC);
      else if( TMath::Abs(pdgCode)==2212)  fHistoClustPhotonPlusConvNeutronBGPtM02[fiCut]->Fill( TruePhotonCandidate->Pt(),clus->GetM02(), fWeightJetJetMC);
      else if( TMath::Abs(pdgCode)==321)   fHistoClustPhotonPlusConvKaonBGPtM02[fiCut]->Fill(    TruePhotonCandidate->Pt(),clus->GetM02(), fWeightJetJetMC);
      else if( TMath::Abs(pdgCode)==130)   fHistoClustPhotonPlusConvK0lBGPtM02[fiCut]->Fill(     TruePhotonCandidate->Pt(),clus->GetM02(), fWeightJetJetMC);
      else                          fHistoClustPhotonPlusConvRestBGPtM02[fiCut]->Fill(    TruePhotonCandidate->Pt(),clus->GetM02(), fWeightJetJetMC);
    }
  }
}

//_________________________________________________________________________________
void AliAnalysisTaskGammaCalo::SetLogBinningXTH2(TH2* histoRebin){
  TAxis *axisafter = histoRebin->GetXaxis();
  Int_t bins = axisafter->GetNbins();
  Double_t from = axisafter->GetXmin();
  Double_t to = axisafter->GetXmax();
  Double_t *newbins = new Double_t[bins+1];
  newbins[0] = from;
  Double_t factor = TMath::Power(to/from, 1./bins);
  for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
  axisafter->Set(bins, newbins);
  delete [] newbins;
}

//_________________________________________________________________________________
Bool_t AliAnalysisTaskGammaCalo::CheckVectorForDoubleCount(vector<Int_t> &vec, Int_t tobechecked)
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
void AliAnalysisTaskGammaCalo::Terminate(const Option_t *)
{
  
  //fOutputContainer->Print(); // Will crash on GRID
}

//________________________________________________________________________
Int_t AliAnalysisTaskGammaCalo::GetSourceClassification(Int_t daughter, Int_t pdgCode){
  
  if (daughter == 111) {
    if (TMath::Abs(pdgCode) == 310) return 1; // k0s
    else if (TMath::Abs(pdgCode) == 3122) return 2; // Lambda
    else if (TMath::Abs(pdgCode) == 130) return 3; // K0L
    else if (TMath::Abs(pdgCode) == 2212) return 4; // proton
    else if (TMath::Abs(pdgCode) == 2112) return 5; // neutron
    else if (TMath::Abs(pdgCode) == 211) return 6; // pion
    else if (TMath::Abs(pdgCode) == 321) return 7; // kaon
    else if (TMath::Abs(pdgCode) == 113 || TMath::Abs(pdgCode) == 213 ) return 8; // rho 0,+,-
    else if (TMath::Abs(pdgCode) == 3222 || TMath::Abs(pdgCode) == 3212 || TMath::Abs(pdgCode) == 3112  ) return 9; // Sigma
    else if (TMath::Abs(pdgCode) == 2224 || TMath::Abs(pdgCode) == 2214 || TMath::Abs(pdgCode) == 2114 || TMath::Abs(pdgCode) == 1114  ) return 10; // Delta
    else if (TMath::Abs(pdgCode) == 313 || TMath::Abs(pdgCode) == 323   ) return 11; // K*
    else return 15;
  }
  return 15;
  
}

//_________________________________________________________________________________
void AliAnalysisTaskGammaCalo::FillMultipleCountMap(map<Int_t,Int_t> &ma, Int_t tobechecked){
  if( ma.find(tobechecked) != ma.end() ) ma[tobechecked] += 1;
  else ma[tobechecked] = 2;
  return;
}

//_________________________________________________________________________________
void AliAnalysisTaskGammaCalo::FillMultipleCountHistoAndClear(map<Int_t,Int_t> &ma, TH1F* hist){
  map<Int_t, Int_t>::iterator it;
  for (it = ma.begin(); it != ma.end(); it++){
    hist->Fill(it->second, fWeightJetJetMC);
  }
  ma.clear();
  return;
}

//_________________________________________________________________________________
void AliAnalysisTaskGammaCalo::DebugMethodPrint1(AliAODConversionMother *pi0cand, AliAODConversionPhoton *gamma0, AliAODConversionPhoton *gamma1){
  if(pi0cand->GetOpeningAngle() < 0.02) DebugMethod(pi0cand,gamma0,gamma1);
  return;
}

//_________________________________________________________________________________
void AliAnalysisTaskGammaCalo::DebugMethod(AliAODConversionMother *pi0cand, AliAODConversionPhoton *gamma0, AliAODConversionPhoton *gamma1){
  if(fLocalDebugFlag != 1) return;
  fstream fOutputLocalDebug;
  fOutputLocalDebug.open("debugOutput.txt",ios::out|ios::app);
  if(!fOutputLocalDebug.is_open()) return;

  AliVCaloCells* cells = fInputEvent->GetEMCALCells();
  fOutputLocalDebug << "--pi0cand--" << endl;
  fOutputLocalDebug << "openingAngle " << pi0cand->GetOpeningAngle() << endl;
  fOutputLocalDebug << "pT " << pi0cand->Pt() << endl;
  fOutputLocalDebug << "invMass " << pi0cand->M() << endl;
  fOutputLocalDebug << "--cluster1--" << endl;
  Float_t clusPos[3]={0,0,0};
  fInputEvent->GetCaloCluster(gamma0->GetCaloClusterRef())->GetPosition(clusPos);
  TVector3 clusterVector(clusPos[0],clusPos[1],clusPos[2]);
  Double_t etaCluster = clusterVector.Eta();
  Double_t phiCluster = clusterVector.Phi();

  Int_t nCellCluster = fInputEvent->GetCaloCluster(gamma0->GetCaloClusterRef())->GetNCells();
  for(Int_t iCell=0;iCell<nCellCluster;iCell++){
    Int_t nSupMod, nModule, nIphi, nIeta;
    ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetGeomEMCAL()->GetCellIndex(
          fInputEvent->GetCaloCluster(gamma0->GetCaloClusterRef())->GetCellAbsId(iCell),
          nSupMod,
          nModule,
          nIphi,
          nIeta);
    Int_t nphi, neta;
    ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetGeomEMCAL()->GetCellPhiEtaIndexInSModule(
          nSupMod,
          nModule,
          nIphi,
          nIeta,
          nphi,
          neta);
    fOutputLocalDebug << fInputEvent->GetCaloCluster(gamma0->GetCaloClusterRef())->GetCellAbsId(iCell) << " " << nSupMod << " " << nphi << " " << neta << " " << cells->GetCellAmplitude(fInputEvent->GetCaloCluster(gamma0->GetCaloClusterRef())->GetCellAbsId(iCell)) << endl;
  }
  fOutputLocalDebug << "phi " << phiCluster << endl;
  fOutputLocalDebug << "eta " << etaCluster << endl;

  fOutputLocalDebug << "--cluster2--" << endl;
  fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef())->GetPosition(clusPos);
  TVector3 clusterVector2(clusPos[0],clusPos[1],clusPos[2]);
  Double_t etaCluster2 = clusterVector2.Eta();
  Double_t phiCluster2 = clusterVector2.Phi();
  nCellCluster = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef())->GetNCells();
  for(Int_t iCell=0;iCell<nCellCluster;iCell++){
    Int_t nSupMod, nModule, nIphi, nIeta;
    ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetGeomEMCAL()->GetCellIndex(
          fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef())->GetCellAbsId(iCell),
          nSupMod,
          nModule,
          nIphi,
          nIeta);
    Int_t nphi, neta;
    ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetGeomEMCAL()->GetCellPhiEtaIndexInSModule(
          nSupMod,
          nModule,
          nIphi,
          nIeta,
          nphi,
          neta);
    fOutputLocalDebug << fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef())->GetCellAbsId(iCell) << " " << nSupMod << " " << nphi << " " << neta << " " << cells->GetCellAmplitude(fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef())->GetCellAbsId(iCell)) << endl;
  }
  fOutputLocalDebug << "phi " << phiCluster2 << endl;
  fOutputLocalDebug << "eta " << etaCluster2 << endl;
  fOutputLocalDebug << "---------------------------------------" << endl;
  fOutputLocalDebug.close();

  return;
}

//_________________________________________________________________________________
void AliAnalysisTaskGammaCalo::EventDebugMethod(){
  if(fLocalDebugFlag != 2) return;
  fstream fOutputLocalDebug;
  fOutputLocalDebug.open("debugOutput.txt",ios::out|ios::app);
  if(!fOutputLocalDebug.is_open()) return;

  Long_t nclus = fInputEvent->GetNumberOfCaloClusters();
  AliVCaloCells* cells = fInputEvent->GetEMCALCells();
  fOutputLocalDebug << "--event--" << endl;
  fOutputLocalDebug << "nclusters " << nclus << endl;

  // vertex
  Double_t vertex[3] = {0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  for(Long_t i = 0; i < nclus; i++){
    AliVCluster* clus = NULL;
    if(fInputEvent->IsA()==AliESDEvent::Class()) clus = new AliESDCaloCluster(*(AliESDCaloCluster*)fInputEvent->GetCaloCluster(i));
    else if(fInputEvent->IsA()==AliAODEvent::Class()) clus = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(i));

    if(!clus) continue;
    if(!((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC,fWeightJetJetMC,i)){
      delete clus;
      continue;
    }
    // TLorentzvector with cluster
    TLorentzVector clusterVector;
    clus->GetMomentum(clusterVector,vertex);

    TLorentzVector* tmpvec = new TLorentzVector();
    tmpvec->SetPxPyPzE(clusterVector.Px(),clusterVector.Py(),clusterVector.Pz(),clusterVector.E());

    // convert to AODConversionPhoton
    AliAODConversionPhoton *PhotonCandidate=new AliAODConversionPhoton(tmpvec);
    if(!PhotonCandidate){ delete clus; delete tmpvec; continue;}

    fOutputLocalDebug << "--cluster" << i << "--" << endl;

    Double_t etaCluster = clusterVector.Eta();
    Double_t phiCluster = clusterVector.Phi();
    Int_t nCellCluster = clus->GetNCells();
    for(Int_t iCell=0;iCell<nCellCluster;iCell++){
      Int_t nSupMod, nModule, nIphi, nIeta;
      ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetGeomEMCAL()->GetCellIndex(
            clus->GetCellAbsId(iCell),
            nSupMod,
            nModule,
            nIphi,
            nIeta);
      Int_t nphi, neta;
      ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetGeomEMCAL()->GetCellPhiEtaIndexInSModule(
            nSupMod,
            nModule,
            nIphi,
            nIeta,
            nphi,
            neta);
      fOutputLocalDebug << clus->GetCellAbsId(iCell) << " " << nSupMod << " " << nphi << " " << neta << " " << cells->GetCellAmplitude(clus->GetCellAbsId(iCell)) << endl;
    }
    fOutputLocalDebug << "phi " << phiCluster << endl;
    fOutputLocalDebug << "eta " << etaCluster << endl;

    delete clus;
    delete tmpvec;
  }
  fOutputLocalDebug << "---------------------------------------" << endl;
  fOutputLocalDebug.close();

  return;
}
