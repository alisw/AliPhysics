/**************************************************************************
 * Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Daniel Muehlheim, Friederike Bock                              *
 * Contributors: Adrian Mechler                                           *
 * Version 1.0                                                            *
 *                                                                        *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
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
#include "AliAnalysisTaskGammaConvCalo.h"
#include "AliVParticle.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliKFVertex.h"
#include "AliGenCocktailEventHeader.h"
#include "AliConversionAODBGHandlerRP.h"
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
#include "AliCaloTriggerMimicHelper.h"
#include <vector>
#include <map>

ClassImp(AliAnalysisTaskGammaConvCalo)

//________________________________________________________________________
AliAnalysisTaskGammaConvCalo::AliAnalysisTaskGammaConvCalo(): AliAnalysisTaskSE(),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fCaloTriggerHelperName(""),
  fCorrTaskSetting(""),
  fBGHandler(NULL),
  fBGHandlerRP(NULL),
  fBGClusHandler(NULL),
  fBGClusHandlerRP(NULL),
  fBGHBTTrueGammaHandler(NULL),
  fBGHBTGenGammaHandler(NULL),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fCutFolder(NULL),
  fESDList(NULL),
  fBackList(NULL),
  fMotherList(NULL),
  fGammaERM02(NULL),
  fInvMassShowerShape(NULL),
  fTrueList(NULL),
  fMCList(NULL),
  fClusterOutputList(NULL),
  fOutputContainer(NULL),
  fReaderGammas(NULL),
  fGammaCandidates(NULL),
  fClusterCandidates(NULL),
  fEventCutArray(NULL),
  fCutArray(NULL),
  fClusterCutArray(NULL),
  fMesonCutArray(NULL),
  fConvJetReader(NULL),
  fOutlierJetReader(NULL),
  fDoJetAnalysis(kFALSE),
  fDoJetQA(kFALSE),
  fJetHistograms(NULL),
  fTrueJetHistograms(NULL),
  fMaxPtNearEMCalPlace(0),
  fJetNearEMCal(kFALSE),
  fTrueGammaCandidatesConv(NULL),
  fTrueGammaCandidatesCalo(NULL),
  fMCGammaCandidates(NULL),
  fAODMCTrackArray(NULL),
  fCaloTriggerMimicHelper(NULL),
  fSetEventCutsOutputlist(),
  fHistoConvGammaPt(NULL),
  fTreeConvGammaPtDcazCat(NULL),
  fPtGamma(0),
  fDCAzPhoton(0),
  fCharCatPhoton(0),
  fCharPhotonMCInfo(0),
  tESDGammaERM02(NULL),
  tESDClusE(0),
  tESDGammaConvR(0),
  tESDClusterM02(0),
  tESDClusterM20(0),
  tESDClusterEta(0),
  tESDClusterPhi(0),
  tESDClusterNCells(0),
  tESDClusterMaxECell(0),
  tESDClusterNLM(0),
  tESDInvMassShowerShape(0),
  tESDIMMesonInvMass(0),
  tESDIMMesonPt(0),
  tESDIMClusE(0),
  tESDIMClusterM02(0),
  tESDIMClusterM20(0),
  tESDIMClusterLeadCellID(0),
  tESDIMClusterClassification(0),
  tESDIMClusMatchedTrackPt(0),
  tESDIMClusTrackDeltaEta(0),
  tESDIMClusTrackDeltaPhi(0),
  tESDIMClusterIsoSumClusterEt(0),
  tESDIMClusterIsoSumTrackEt(0),
  tESDmapIsClusterAcceptedWithoutTrackMatch(),
  fHistoMotherInvMassPt(NULL),
  fHistoMotherMatchedInvMassPt(NULL),
  fSparseMotherInvMassPtZM(NULL),
  fHistoMotherBackInvMassPt(NULL),
  fSparseMotherBackInvMassPtZM(NULL),
  fHistoMotherInvMassPtAlpha(NULL),
  fHistoMotherPi0PtY(NULL),
  fHistoMotherEtaPtY(NULL),
  fHistoMotherPi0PtAlpha(NULL),
  fHistoMotherEtaPtAlpha(NULL),
  fHistoMotherPi0PtOpenAngle(NULL),
  fHistoMotherEtaPtOpenAngle(NULL),
  fHistoMotherPi0ConvPhotonEtaPhi(NULL),
  fHistoMotherEtaConvPhotonEtaPhi(NULL),
  fHistoMotherInvMassECalib(NULL),
  fHistoMotherBackInvMassECalib(NULL),
  fHistoMotherInvMassECalibPCM(NULL),
  fHistoMotherBackInvMassECalibPCM(NULL),
  fHistoPhotonPairPtconv(NULL),
  fHistoPhotonPairMixedEventPtconv(NULL),
  fHistoClusGammaPt(NULL),
  fHistoClusGammaE(NULL),
  fHistoClusGammaE_BothBM(NULL),
  fHistoClusGammaE_BothBM_highestE(NULL),
  fHistoClusGammaE_AnaBM_highestE(NULL),
  fHistoClusGammaE_onlyTriggered(NULL),
  fHistoGoodMesonClusters(NULL),
  fHistoClusOverlapHeadersGammaPt(NULL),
  fHistoClusAllHeadersGammaPt(NULL),
  fHistoClusRejectedHeadersGammaPt(NULL),
  fHistoMCHeaders(NULL),
  fHistoMCAllGammaPt(NULL),
  fHistoMCAllGammaEMCALAccPt(NULL),
  fHistoMCAllSecondaryGammaPt(NULL),
  fHistoMCDecayGammaPi0Pt(NULL),
  fHistoMCDecayGammaRhoPt(NULL),
  fHistoMCDecayGammaEtaPt(NULL),
  fHistoMCDecayGammaOmegaPt(NULL),
  fHistoMCDecayGammaEtapPt(NULL),
  fHistoMCDecayGammaPhiPt(NULL),
  fHistoMCDecayGammaSigmaPt(NULL),
  fHistoMCConvGammaPt(NULL),
  fHistoMCSecondaryConvGammaPt(NULL),
  fHistoMCPi0Pt(NULL),
  fHistoMCPi0WOWeightPt(NULL),
  fHistoMCPi0WOEvtWeightPt(NULL),
  fHistoMCEtaPt(NULL),
  fHistoMCEtaWOWeightPt(NULL),
  fHistoMCEtaWOEvtWeightPt(NULL),
  fHistoMCPi0InAccPt(NULL),
  fHistoMCPi0WOWeightInAccPt(NULL),
  fHistoMCPi0WOEvtWeightInAccPt(NULL),
  fHistoMCEtaInAccPt(NULL),
  fHistoMCEtaWOWeightInAccPt(NULL),
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
  fHistoMCPi0PtGammaLeg(NULL),
  fHistoMCPi0WOWeightPtGammaLeg(NULL),
  fHistoMCPi0InAccPtGammaLeg(NULL),
  fHistoMCPi0WOWeightInAccPtGammaLeg(NULL),
  fHistoMCSecPi0PtGamma1vsSource(NULL),
  fHistoMCSecPi0InAccPtGamma1vsSource(NULL),
  fHistoMCSecPi0PtGamma2vsSource(NULL),
  fHistoMCSecPi0InAccPtGamma2vsSource(NULL),
  fHistoTruePi0InvMassPt(NULL),
  fHistoTrueEtaInvMassPt(NULL),
  fHistoTruePi0InvMassPtAdditional(NULL),
  fHistoTrueEtaInvMassPtAdditional(NULL),
  fHistoTruePi0MatchedInvMassPt(NULL),
  fHistoTrueEtaMatchedInvMassPt(NULL),
  fHistoTruePi0CaloPhotonInvMassPt(NULL),
  fHistoTrueEtaCaloPhotonInvMassPt(NULL),
  fHistoTruePi0CaloConvertedPhotonInvMassPt(NULL),
  fHistoTruePi0CaloConvertedPhotonMatchedInvMassPt(NULL),
  fHistoTruePi0CaloConvPhotonConvRPt(NULL),
  fHistoTruePi0CaloConvPhotonConvRAlpha(NULL),
  fHistoTruePi0CaloConvPhotonPtAlpha(NULL),
  fHistoTrueEtaCaloConvertedPhotonInvMassPt(NULL),
  fHistoTrueEtaCaloConvertedPhotonMatchedInvMassPt(NULL),
  fHistoTrueEtaCaloConvPhotonConvRPt(NULL),
  fHistoTrueEtaCaloConvPhotonConvRAlpha(NULL),
  fHistoTrueEtaCaloConvPhotonPtAlpha(NULL),
  fHistoTruePi0CaloElectronInvMassPt(NULL),
  fHistoTrueEtaCaloElectronInvMassPt(NULL),
  fHistoTruePi0CaloMergedClusterInvMassPt(NULL),
  fHistoTrueEtaCaloMergedClusterInvMassPt(NULL),
  fHistoTrueMotherCaloEMNonLeadingInvMassPt(NULL),
  fHistoTruePi0CaloMergedClusterPartConvInvMassPt(NULL),
  fHistoTrueEtaCaloMergedClusterPartConvInvMassPt(NULL),
  fHistoTruePrimaryPi0InvMassPt(NULL),
  fHistoTruePrimaryEtaInvMassPt(NULL),
  fHistoTruePrimaryPi0W0WeightingInvMassPt(NULL),
  fHistoTruePrimaryEtaW0WeightingInvMassPt(NULL),
  fProfileTruePrimaryPi0WeightsInvMassPt(NULL),
  fProfileTruePrimaryEtaWeightsInvMassPt(NULL),
  fHistoTruePrimaryPi0MCPtResolPt(NULL),
  fHistoTruePrimaryEtaMCPtResolPt(NULL),
  fHistoTrueMotherPi0ConvPhotonEtaPhi(NULL),
  fHistoTrueMotherEtaConvPhotonEtaPhi(NULL),
  fHistoTrueSecondaryPi0InvMassPt(NULL),
  fHistoTrueSecondaryPi0FromK0sInvMassPt(NULL),
  fHistoTrueK0sWithPi0DaughterMCPt(NULL),
  fHistoTrueSecondaryPi0FromK0lInvMassPt(NULL),
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
  fHistoTrueConvGammaPt(NULL),
  fHistoCombinatorialPt(NULL),
  fHistoTruePrimaryConvGammaPt(NULL),
  fHistoTruePrimaryConvGammaESDPtMCPt(NULL),
  fHistoTrueSecondaryConvGammaPt(NULL),
  fHistoTrueSecondaryConvGammaMCPt(NULL),
  fHistoTrueSecondaryConvGammaFromXFromK0sMCPtESDPt(NULL),
  fHistoTrueSecondaryConvGammaFromXFromK0lMCPtESDPt(NULL),
  fHistoTrueSecondaryConvGammaFromXFromLambdaMCPtESDPt(NULL),
  fHistoTrueClusGammaPt(NULL),
  fHistoTrueClusElectronPt(NULL),
  fHistoTrueClusConvGammaPt(NULL),
  fHistoTrueClusConvGammaFullyPt(NULL),
  fHistoTrueClusMergedGammaPt(NULL),
  fHistoTrueClusMergedPartConvGammaPt(NULL),
  fHistoTrueClusDalitzPt(NULL),
  fHistoTrueClusDalitzMergedPt(NULL),
  fHistoTrueClusPhotonFromElecMotherPt(NULL),
  fHistoTrueClusShowerPt(NULL),
  fHistoTrueClusSubLeadingPt(NULL),
  fHistoTrueClusNMothers(NULL),
  fHistoTrueClusEMNonLeadingPt(NULL),
  fHistoTrueNLabelsInClusPt(NULL),
  fHistoTruePrimaryClusGammaPt(NULL),
  fHistoTruePrimaryClusGammaESDPtMCPt(NULL),
  fHistoTruePrimaryClusConvGammaPt(NULL),
  fHistoTruePrimaryClusConvGammaESDPtMCPt(NULL),
  fHistoTrueSecondaryClusGammaPt(NULL),
  fHistoTrueSecondaryClusGammaFromK0sPt(NULL),
  fHistoTrueSecondaryClusGammaFromK0lPt(NULL),
  fHistoTrueSecondaryClusGammaFromLambdaPt(NULL),
  fHistoTruePrimaryPi0PhotonPairPtconv(NULL),
  fHistoTruePrimaryPi0W0WeightsPhotonPairPtconv(NULL),
  fHistoTruePrimaryPi0DCPtconv(NULL),
  fHistoTruePrimaryPi0MissingPtconv(NULL),
  fHistoTrueSecondaryPi0PhotonPairPtconv(NULL),
  fHistoTrueSecondaryPi0FromK0sPhotonPairPtconv(NULL),
  fHistoTrueSecondaryPi0FromK0lPhotonPairPtconv(NULL),
  fHistoTrueSecondaryPi0FromLambdaPhotonPairPtconv(NULL),
  fHistoTrueSecondaryPi0DCPtconvSource(NULL),
  fHistoTrueSecondaryPi0MissingPtconvSource(NULL),
  fHistoTrueGammaInvMassPt(NULL),
  fHistoMCGammaHBTInvMassPt(NULL),
  fHistoBckTrueGammaInvMassPt(NULL),
  fHistoBckMCGammaHBTInvMassPt(NULL),
  fHistoHBTOpeningAnglePt(NULL),
  fHistoBckHBTOpeningAnglePt(NULL),
  fHistoHBTDeltaEPt(NULL),
  fHistoBckHBTDeltaEPt(NULL),
  fHistoMCGammaHBTOpeningAnglePt(NULL),
  fHistoBckMCGammaHBTOpeningAnglePt(NULL),
  fHistoMCGammaHBTDeltaEPt(NULL),
  fHistoBckMCGammaHBTDeltaEPt(NULL),
  fHistoPtJet(NULL),
  fHistoJetEta(NULL),
  fHistoJetPhi(NULL),
  fHistoJetArea(NULL),
  fHistoNJets(NULL),
  fHistoEventwJets(NULL),
  fHistoJetPi0PtRatio(NULL),
  fHistoDoubleCounting(NULL),
  fHistoJetMotherInvMassPt(NULL),
  fHistoPi0InJetMotherInvMassPt(NULL),
  fHistoMotherBackJetInvMassPt(NULL),
  fHistoRJetPi0Cand(NULL),
  fHistoEtaPhiJetPi0Cand(NULL),
  fHistoEtaPhiJetWithPi0Cand(NULL),
  fHistoJetFragmFunc(NULL),
  fHistoJetFragmFuncZInvMass(NULL),
  fHistoTruevsRecJetPt(NULL),
  fHistoTruePi0JetMotherInvMassPt(NULL),
  fHistoTruePi0InJetMotherInvMassPt(NULL),
  fHistoTruePrimaryPi0JetInvMassPt(NULL),
  fHistoTruePrimaryPi0inJetInvMassPt(NULL),
  fHistoTruePrimaryPi0InJetInvMassTruePt(NULL),
  fHistoTrueDoubleCountingPi0Jet(NULL),
  fHistoTrueEtaJetMotherInvMassPt(NULL),
  fHistoTrueEtaInJetMotherInvMassPt(NULL),
  fHistoTruePrimaryEtaJetInvMassPt(NULL),
  fHistoTruePrimaryEtainJetInvMassPt(NULL),
  fHistoTrueDoubleCountingEtaJet(NULL),
  fHistoTruePi0JetFragmFunc(NULL),
  fHistoTruePi0JetFragmFuncZInvMass(NULL),
  fHistoTrueEtaJetFragmFunc(NULL),
  fHistoTrueEtaJetFragmFuncZInvMass(NULL),
  fHistoMCPi0JetInAccPt(NULL),
  fHistoMCPi0inJetInAccPt(NULL),
  fHistoMCEtaJetInAccPt(NULL),
  fHistoMCEtainJetInAccPt(NULL),
  fHistoMCPi0JetEventGenerated(NULL),
  fHistoMCPi0inJetGenerated(NULL),
  fHistoMCEtaJetEventGenerated(NULL),
  fHistoMCEtainJetGenerated(NULL),
  fHistoTrueSecondaryPi0FromK0sJetInvMassPt(NULL),
  fHistoTrueSecondaryPi0FromK0sinJetInvMassPt(NULL),
  fHistoTrueSecondaryPi0FromLambdaJetInvMassPt(NULL),
  fHistoTrueSecondaryPi0FromLambdainJetInvMassPt(NULL),
  fHistoTrueSecondaryPi0FromK0lJetInvMassPt(NULL),
  fHistoTrueSecondaryPi0FromK0linJetInvMassPt(NULL),
  fHistoTrueSecondaryPi0InvJetMassPt(NULL),
  fHistoTrueSecondaryPi0InvinJetMassPt(NULL),
  fHistoMotherPi0inJetPtY(NULL),
  fHistoMotherEtainJetPtY(NULL),
  fHistoMotherPi0inJetPtPhi(NULL),
  fHistoMotherEtainJetPtPhi(NULL),
  fHistoUnfoldingAsData(NULL),
  fHistoUnfoldingMissed(NULL),
  fHistoUnfoldingReject(NULL),
  fHistoUnfoldingAsDataInvMassZ(NULL),
  fHistoUnfoldingMissedInvMassZ(NULL),
  fHistoUnfoldingRejectInvMassZ(NULL),
  fVectorJetPt(0),
  fVectorJetPx(0),
  fVectorJetPy(0),
  fVectorJetPz(0),
  fVectorJetEta(0),
  fVectorJetPhi(0),
  fVectorJetArea(0),
  fTrueVectorJetPt(0),
  fTrueVectorJetPx(0),
  fTrueVectorJetPy(0),
  fTrueVectorJetPz(0),
  fTrueVectorJetEta(0),
  fTrueVectorJetPhi(0),
  fVectorRecTruePi0s(0),
  fVectorRecTrueEtas(0),
  fHistoDoubleCountTruePi0InvMassPt(NULL),
  fHistoDoubleCountTrueEtaInvMassPt(NULL),
  fHistoDoubleCountTrueConvGammaRPt(NULL),
  fHistoDoubleCountTrueClusterGammaPt(NULL),
  fVectorDoubleCountTruePi0s(0),
  fVectorDoubleCountTrueEtas(0),
  fVectorDoubleCountTrueConvGammas(0),
  fVectorDoubleCountTrueClusterGammas(0),
  fHistoMultipleCountTruePi0(NULL),
  fHistoMultipleCountTrueEta(NULL),
  fHistoMultipleCountTrueConvGamma(NULL),
  fHistoMultipleCountTrueClusterGamma(NULL),
  fMapMultipleCountTruePi0s(),
  fMapMultipleCountTrueEtas(),
  fMapMultipleCountTrueConvGammas(),
  fMapMultipleCountTrueClusterGammas(),
  fHistoTrueClusGammaEM02(NULL),
  fHistoTrueClusPi0EM02(NULL),
  fHistoTruePi0InvMassECalib(NULL),
  fHistoTruePi0PureGammaInvMassECalib(NULL),
  fHistoTruePi0InvMassECalibPCM(NULL),
  fHistoNEvents(NULL),
  fHistoNEventsWOWeight(NULL),
  fHistoNGoodESDTracks(NULL),
  fHistoVertexZ(NULL),
  fHistoVertexX(NULL),
  fHistoVertexY(NULL),
  fHistoNGammaCandidates(NULL),
  fHistoNGoodESDTracksVsNGammaCandidates(NULL),
  fHistoSPDClusterTrackletBackground(NULL),
  fHistoNV0Tracks(NULL),
  fProfileEtaShift(NULL),
  fProfileJetJetXSection(NULL),
  fHistoJetJetNTrials(NULL),
  fHistoEventSphericity(NULL),
  fEventPlaneAngle(-100),
  tESDMesonsInvMassPtDcazMinDcazMaxFlag(NULL),
  fInvMass(0),
  fPt(0),
  fDCAzGammaMin(0),
  fDCAzGammaMax(0),
  iFlag(0),
  iMesonMCInfo(0),
  fRandom(0),
  fNGammaCandidates(0),
  fUnsmearedPx(NULL),
  fUnsmearedPy(NULL),
  fUnsmearedPz(NULL),
  fUnsmearedE(NULL),
  fMCEventPos(NULL),
  fMCEventNeg(NULL),
  fESDArrayPos(NULL),
  fESDArrayNeg(NULL),
  fnCuts(0),
  fiCut(0),
  fMoveParticleAccordingToVertex(kTRUE),
  fIsHeavyIon(0),
  fDoLightOutput(0),
  fPi0EtaSwitch(0),
  fDoECalibOutput(kFALSE),
  fDoMesonAnalysis(kTRUE),
  fDoMesonQA(0),
  fDoPhotonQA(0),
  fDoClusterQA(0),
  fIsFromDesiredHeader(kTRUE),
  fIsOverlappingWithOtherHeader(kFALSE),
  fIsMC(0),
  fDoTHnSparse(kTRUE),
  fSetPlotHistsExtQA(kFALSE),
  fWeightJetJetMC(1),
  fDoConvGammaShowerShapeTree(kFALSE),
  fEnableSortForClusMC(kFALSE),
  fDoPrimaryTrackMatching(kFALSE),
  fDoInvMassShowerShapeTree(kFALSE),
  fDoMaterialBudgetWeightingOfGammasForTrueMesons(kFALSE),
  tBrokenFiles(NULL),
  fFileNameBroken(NULL),
  fAllowOverlapHeaders(kTRUE),
  fTrackMatcherRunningMode(0),
  fDoHBTHistoOutput(kFALSE)
{

}

//________________________________________________________________________
AliAnalysisTaskGammaConvCalo::AliAnalysisTaskGammaConvCalo(const char *name):
  AliAnalysisTaskSE(name),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fCaloTriggerHelperName(""),
  fCorrTaskSetting(""),
  fBGHandler(NULL),
  fBGHandlerRP(NULL),
  fBGClusHandler(NULL),
  fBGClusHandlerRP(NULL),
  fBGHBTTrueGammaHandler(NULL),
  fBGHBTGenGammaHandler(NULL),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fCutFolder(NULL),
  fESDList(NULL),
  fBackList(NULL),
  fMotherList(NULL),
  fGammaERM02(NULL),
  fInvMassShowerShape(NULL),
  fTrueList(NULL),
  fMCList(NULL),
  fClusterOutputList(NULL),
  fOutputContainer(0),
  fReaderGammas(NULL),
  fGammaCandidates(NULL),
  fClusterCandidates(NULL),
  fEventCutArray(NULL),
  fCutArray(NULL),
  fClusterCutArray(NULL),
  fMesonCutArray(NULL),
  fConvJetReader(NULL),
  fOutlierJetReader(NULL),
  fDoJetAnalysis(kFALSE),
  fDoJetQA(kFALSE),
  fJetHistograms(NULL),
  fTrueJetHistograms(NULL),
  fMaxPtNearEMCalPlace(0),
  fJetNearEMCal(kFALSE),
  fTrueGammaCandidatesConv(NULL),
  fTrueGammaCandidatesCalo(NULL),
  fMCGammaCandidates(NULL),
  fAODMCTrackArray(NULL),
  fCaloTriggerMimicHelper(NULL),
  fSetEventCutsOutputlist(),
  fHistoConvGammaPt(NULL),
  fTreeConvGammaPtDcazCat(NULL),
  fPtGamma(0),
  fDCAzPhoton(0),
  fCharCatPhoton(0),
  fCharPhotonMCInfo(0),
  tESDGammaERM02(NULL),
  tESDClusE(0),
  tESDGammaConvR(0),
  tESDClusterM02(0),
  tESDClusterM20(0),
  tESDClusterEta(0),
  tESDClusterPhi(0),
  tESDClusterNCells(0),
  tESDClusterMaxECell(0),
  tESDClusterNLM(0),
  tESDInvMassShowerShape(0),
  tESDIMMesonInvMass(0),
  tESDIMMesonPt(0),
  tESDIMClusE(0),
  tESDIMClusterM02(0),
  tESDIMClusterM20(0),
  tESDIMClusterLeadCellID(0),
  tESDIMClusterClassification(0),
  tESDIMClusMatchedTrackPt(0),
  tESDIMClusTrackDeltaEta(0),
  tESDIMClusTrackDeltaPhi(0),
  tESDIMClusterIsoSumClusterEt(0),
  tESDIMClusterIsoSumTrackEt(0),
  tESDmapIsClusterAcceptedWithoutTrackMatch(),
  fHistoMotherInvMassPt(NULL),
  fHistoMotherMatchedInvMassPt(NULL),
  fSparseMotherInvMassPtZM(NULL),
  fHistoMotherBackInvMassPt(NULL),
  fSparseMotherBackInvMassPtZM(NULL),
  fHistoMotherInvMassPtAlpha(NULL),
  fHistoMotherPi0PtY(NULL),
  fHistoMotherEtaPtY(NULL),
  fHistoMotherPi0PtAlpha(NULL),
  fHistoMotherEtaPtAlpha(NULL),
  fHistoMotherPi0PtOpenAngle(NULL),
  fHistoMotherEtaPtOpenAngle(NULL),
  fHistoMotherPi0ConvPhotonEtaPhi(NULL),
  fHistoMotherEtaConvPhotonEtaPhi(NULL),
  fHistoMotherInvMassECalib(NULL),
  fHistoMotherBackInvMassECalib(NULL),
  fHistoMotherInvMassECalibPCM(NULL),
  fHistoMotherBackInvMassECalibPCM(NULL),
  fHistoPhotonPairPtconv(NULL),
  fHistoPhotonPairMixedEventPtconv(NULL),
  fHistoClusGammaPt(NULL),
  fHistoClusGammaE(NULL),
  fHistoClusGammaE_BothBM(NULL),
  fHistoClusGammaE_BothBM_highestE(NULL),
  fHistoClusGammaE_AnaBM_highestE(NULL),
  fHistoClusGammaE_onlyTriggered(NULL),
  fHistoGoodMesonClusters(NULL),
  fHistoClusOverlapHeadersGammaPt(NULL),
  fHistoClusAllHeadersGammaPt(NULL),
  fHistoClusRejectedHeadersGammaPt(NULL),
  fHistoMCHeaders(NULL),
  fHistoMCAllGammaPt(NULL),
  fHistoMCAllGammaEMCALAccPt(NULL),
  fHistoMCAllSecondaryGammaPt(NULL),
  fHistoMCDecayGammaPi0Pt(NULL),
  fHistoMCDecayGammaRhoPt(NULL),
  fHistoMCDecayGammaEtaPt(NULL),
  fHistoMCDecayGammaOmegaPt(NULL),
  fHistoMCDecayGammaEtapPt(NULL),
  fHistoMCDecayGammaPhiPt(NULL),
  fHistoMCDecayGammaSigmaPt(NULL),
  fHistoMCConvGammaPt(NULL),
  fHistoMCSecondaryConvGammaPt(NULL),
  fHistoMCPi0Pt(NULL),
  fHistoMCPi0WOWeightPt(NULL),
  fHistoMCPi0WOEvtWeightPt(NULL),
  fHistoMCEtaPt(NULL),
  fHistoMCEtaWOWeightPt(NULL),
  fHistoMCEtaWOEvtWeightPt(NULL),
  fHistoMCPi0InAccPt(NULL),
  fHistoMCPi0WOWeightInAccPt(NULL),
  fHistoMCPi0WOEvtWeightInAccPt(NULL),
  fHistoMCEtaInAccPt(NULL),
  fHistoMCEtaWOWeightInAccPt(NULL),
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
  fHistoMCPi0PtGammaLeg(NULL),
  fHistoMCPi0WOWeightPtGammaLeg(NULL),
  fHistoMCPi0InAccPtGammaLeg(NULL),
  fHistoMCPi0WOWeightInAccPtGammaLeg(NULL),
  fHistoMCSecPi0PtGamma1vsSource(NULL),
  fHistoMCSecPi0InAccPtGamma1vsSource(NULL),
  fHistoMCSecPi0PtGamma2vsSource(NULL),
  fHistoMCSecPi0InAccPtGamma2vsSource(NULL),
  fHistoTruePi0InvMassPt(NULL),
  fHistoTrueEtaInvMassPt(NULL),
  fHistoTruePi0InvMassPtAdditional(NULL),
  fHistoTrueEtaInvMassPtAdditional(NULL),
  fHistoTruePi0MatchedInvMassPt(NULL),
  fHistoTrueEtaMatchedInvMassPt(NULL),
  fHistoTruePi0CaloPhotonInvMassPt(NULL),
  fHistoTrueEtaCaloPhotonInvMassPt(NULL),
  fHistoTruePi0CaloConvertedPhotonInvMassPt(NULL),
  fHistoTruePi0CaloConvertedPhotonMatchedInvMassPt(NULL),
  fHistoTruePi0CaloConvPhotonConvRPt(NULL),
  fHistoTruePi0CaloConvPhotonConvRAlpha(NULL),
  fHistoTruePi0CaloConvPhotonPtAlpha(NULL),
  fHistoTrueEtaCaloConvertedPhotonInvMassPt(NULL),
  fHistoTrueEtaCaloConvertedPhotonMatchedInvMassPt(NULL),
  fHistoTrueEtaCaloConvPhotonConvRPt(NULL),
  fHistoTrueEtaCaloConvPhotonConvRAlpha(NULL),
  fHistoTrueEtaCaloConvPhotonPtAlpha(NULL),
  fHistoTruePi0CaloElectronInvMassPt(NULL),
  fHistoTrueEtaCaloElectronInvMassPt(NULL),
  fHistoTruePi0CaloMergedClusterInvMassPt(NULL),
  fHistoTrueEtaCaloMergedClusterInvMassPt(NULL),
  fHistoTrueMotherCaloEMNonLeadingInvMassPt(NULL),
  fHistoTruePi0CaloMergedClusterPartConvInvMassPt(NULL),
  fHistoTrueEtaCaloMergedClusterPartConvInvMassPt(NULL),
  fHistoTruePrimaryPi0InvMassPt(NULL),
  fHistoTruePrimaryEtaInvMassPt(NULL),
  fHistoTruePrimaryPi0W0WeightingInvMassPt(NULL),
  fHistoTruePrimaryEtaW0WeightingInvMassPt(NULL),
  fProfileTruePrimaryPi0WeightsInvMassPt(NULL),
  fProfileTruePrimaryEtaWeightsInvMassPt(NULL),
  fHistoTruePrimaryPi0MCPtResolPt(NULL),
  fHistoTruePrimaryEtaMCPtResolPt(NULL),
  fHistoTrueMotherPi0ConvPhotonEtaPhi(NULL),
  fHistoTrueMotherEtaConvPhotonEtaPhi(NULL),
  fHistoTrueSecondaryPi0InvMassPt(NULL),
  fHistoTrueSecondaryPi0FromK0sInvMassPt(NULL),
  fHistoTrueK0sWithPi0DaughterMCPt(NULL),
  fHistoTrueSecondaryPi0FromK0lInvMassPt(NULL),
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
  fHistoTrueConvGammaPt(NULL),
  fHistoCombinatorialPt(NULL),
  fHistoTruePrimaryConvGammaPt(NULL),
  fHistoTruePrimaryConvGammaESDPtMCPt(NULL),
  fHistoTrueSecondaryConvGammaPt(NULL),
  fHistoTrueSecondaryConvGammaMCPt(NULL),
  fHistoTrueSecondaryConvGammaFromXFromK0sMCPtESDPt(NULL),
  fHistoTrueSecondaryConvGammaFromXFromK0lMCPtESDPt(NULL),
  fHistoTrueSecondaryConvGammaFromXFromLambdaMCPtESDPt(NULL),
  fHistoTrueClusGammaPt(NULL),
  fHistoTrueClusElectronPt(NULL),
  fHistoTrueClusConvGammaPt(NULL),
  fHistoTrueClusConvGammaFullyPt(NULL),
  fHistoTrueClusMergedGammaPt(NULL),
  fHistoTrueClusMergedPartConvGammaPt(NULL),
  fHistoTrueClusDalitzPt(NULL),
  fHistoTrueClusDalitzMergedPt(NULL),
  fHistoTrueClusPhotonFromElecMotherPt(NULL),
  fHistoTrueClusShowerPt(NULL),
  fHistoTrueClusSubLeadingPt(NULL),
  fHistoTrueClusNMothers(NULL),
  fHistoTrueClusEMNonLeadingPt(NULL),
  fHistoTrueNLabelsInClusPt(NULL),
  fHistoTruePrimaryClusGammaPt(NULL),
  fHistoTruePrimaryClusGammaESDPtMCPt(NULL),
  fHistoTruePrimaryClusConvGammaPt(NULL),
  fHistoTruePrimaryClusConvGammaESDPtMCPt(NULL),
  fHistoTrueSecondaryClusGammaPt(NULL),
  fHistoTrueSecondaryClusGammaFromK0sPt(NULL),
  fHistoTrueSecondaryClusGammaFromK0lPt(NULL),
  fHistoTrueSecondaryClusGammaFromLambdaPt(NULL),
  fHistoTruePrimaryPi0PhotonPairPtconv(NULL),
  fHistoTruePrimaryPi0W0WeightsPhotonPairPtconv(NULL),
  fHistoTruePrimaryPi0DCPtconv(NULL),
  fHistoTruePrimaryPi0MissingPtconv(NULL),
  fHistoTrueSecondaryPi0PhotonPairPtconv(NULL),
  fHistoTrueSecondaryPi0FromK0sPhotonPairPtconv(NULL),
  fHistoTrueSecondaryPi0FromK0lPhotonPairPtconv(NULL),
  fHistoTrueSecondaryPi0FromLambdaPhotonPairPtconv(NULL),
  fHistoTrueSecondaryPi0DCPtconvSource(NULL),
  fHistoTrueSecondaryPi0MissingPtconvSource(NULL),
  fHistoTrueGammaInvMassPt(NULL),
  fHistoMCGammaHBTInvMassPt(NULL),
  fHistoBckTrueGammaInvMassPt(NULL),
  fHistoBckMCGammaHBTInvMassPt(NULL),
  fHistoHBTOpeningAnglePt(NULL),
  fHistoBckHBTOpeningAnglePt(NULL),
  fHistoHBTDeltaEPt(NULL),
  fHistoBckHBTDeltaEPt(NULL),
  fHistoMCGammaHBTOpeningAnglePt(NULL),
  fHistoBckMCGammaHBTOpeningAnglePt(NULL),
  fHistoMCGammaHBTDeltaEPt(NULL),
  fHistoBckMCGammaHBTDeltaEPt(NULL),
  fHistoPtJet(NULL),
  fHistoJetEta(NULL),
  fHistoJetPhi(NULL),
  fHistoJetArea(NULL),
  fHistoNJets(NULL),
  fHistoEventwJets(NULL),
  fHistoJetPi0PtRatio(NULL),
  fHistoDoubleCounting(NULL),
  fHistoJetMotherInvMassPt(NULL),
  fHistoPi0InJetMotherInvMassPt(NULL),
  fHistoMotherBackJetInvMassPt(NULL),
  fHistoRJetPi0Cand(NULL),
  fHistoEtaPhiJetPi0Cand(NULL),
  fHistoEtaPhiJetWithPi0Cand(NULL),
  fHistoJetFragmFunc(NULL),
  fHistoJetFragmFuncZInvMass(NULL),
  fHistoTruevsRecJetPt(NULL),
  fHistoTruePi0JetMotherInvMassPt(NULL),
  fHistoTruePi0InJetMotherInvMassPt(NULL),
  fHistoTruePrimaryPi0JetInvMassPt(NULL),
  fHistoTruePrimaryPi0inJetInvMassPt(NULL),
  fHistoTruePrimaryPi0InJetInvMassTruePt(NULL),
  fHistoTrueDoubleCountingPi0Jet(NULL),
  fHistoTrueEtaJetMotherInvMassPt(NULL),
  fHistoTrueEtaInJetMotherInvMassPt(NULL),
  fHistoTruePrimaryEtaJetInvMassPt(NULL),
  fHistoTruePrimaryEtainJetInvMassPt(NULL),
  fHistoTrueDoubleCountingEtaJet(NULL),
  fHistoTruePi0JetFragmFunc(NULL),
  fHistoTruePi0JetFragmFuncZInvMass(NULL),
  fHistoTrueEtaJetFragmFunc(NULL),
  fHistoTrueEtaJetFragmFuncZInvMass(NULL),
  fHistoMCPi0JetInAccPt(NULL),
  fHistoMCPi0inJetInAccPt(NULL),
  fHistoMCEtaJetInAccPt(NULL),
  fHistoMCEtainJetInAccPt(NULL),
  fHistoMCPi0JetEventGenerated(NULL),
  fHistoMCPi0inJetGenerated(NULL),
  fHistoMCEtaJetEventGenerated(NULL),
  fHistoMCEtainJetGenerated(NULL),
  fHistoTrueSecondaryPi0FromK0sJetInvMassPt(NULL),
  fHistoTrueSecondaryPi0FromK0sinJetInvMassPt(NULL),
  fHistoTrueSecondaryPi0FromLambdaJetInvMassPt(NULL),
  fHistoTrueSecondaryPi0FromLambdainJetInvMassPt(NULL),
  fHistoTrueSecondaryPi0FromK0lJetInvMassPt(NULL),
  fHistoTrueSecondaryPi0FromK0linJetInvMassPt(NULL),
  fHistoTrueSecondaryPi0InvJetMassPt(NULL),
  fHistoTrueSecondaryPi0InvinJetMassPt(NULL),
  fHistoMotherPi0inJetPtY(NULL),
  fHistoMotherEtainJetPtY(NULL),
  fHistoMotherPi0inJetPtPhi(NULL),
  fHistoMotherEtainJetPtPhi(NULL),
  fHistoUnfoldingAsData(NULL),
  fHistoUnfoldingMissed(NULL),
  fHistoUnfoldingReject(NULL),
  fHistoUnfoldingAsDataInvMassZ(NULL),
  fHistoUnfoldingMissedInvMassZ(NULL),
  fHistoUnfoldingRejectInvMassZ(NULL),
  fVectorJetPt(0),
  fVectorJetPx(0),
  fVectorJetPy(0),
  fVectorJetPz(0),
  fVectorJetEta(0),
  fVectorJetPhi(0),
  fVectorJetArea(0),
  fTrueVectorJetPt(0),
  fTrueVectorJetPx(0),
  fTrueVectorJetPy(0),
  fTrueVectorJetPz(0),
  fTrueVectorJetEta(0),
  fTrueVectorJetPhi(0),
  fVectorRecTruePi0s(0),
  fVectorRecTrueEtas(0),
  fHistoDoubleCountTruePi0InvMassPt(NULL),
  fHistoDoubleCountTrueEtaInvMassPt(NULL),
  fHistoDoubleCountTrueConvGammaRPt(NULL),
  fHistoDoubleCountTrueClusterGammaPt(NULL),
  fVectorDoubleCountTruePi0s(0),
  fVectorDoubleCountTrueEtas(0),
  fVectorDoubleCountTrueConvGammas(0),
  fVectorDoubleCountTrueClusterGammas(0),
  fHistoMultipleCountTruePi0(NULL),
  fHistoMultipleCountTrueEta(NULL),
  fHistoMultipleCountTrueConvGamma(NULL),
  fHistoMultipleCountTrueClusterGamma(NULL),
  fMapMultipleCountTruePi0s(),
  fMapMultipleCountTrueEtas(),
  fMapMultipleCountTrueConvGammas(),
  fMapMultipleCountTrueClusterGammas(),
  fHistoTrueClusGammaEM02(NULL),
  fHistoTrueClusPi0EM02(NULL),
  fHistoTruePi0InvMassECalib(NULL),
  fHistoTruePi0PureGammaInvMassECalib(NULL),
  fHistoTruePi0InvMassECalibPCM(NULL),
  fHistoNEvents(NULL),
  fHistoNEventsWOWeight(NULL),
  fHistoNGoodESDTracks(NULL),
  fHistoVertexZ(NULL),
  fHistoVertexX(NULL),
  fHistoVertexY(NULL),
  fHistoNGammaCandidates(NULL),
  fHistoNGoodESDTracksVsNGammaCandidates(NULL),
  fHistoSPDClusterTrackletBackground(NULL),
  fHistoNV0Tracks(NULL),
  fProfileEtaShift(NULL),
  fProfileJetJetXSection(NULL),
  fHistoJetJetNTrials(NULL),
  fHistoEventSphericity(NULL),
  fEventPlaneAngle(-100),
  tESDMesonsInvMassPtDcazMinDcazMaxFlag(NULL),
  fInvMass(0),
  fPt(0),
  fDCAzGammaMin(0),
  fDCAzGammaMax(0),
  iFlag(0),
  iMesonMCInfo(0),
  fRandom(0),
  fNGammaCandidates(0),
  fUnsmearedPx(NULL),
  fUnsmearedPy(NULL),
  fUnsmearedPz(NULL),
  fUnsmearedE(NULL),
  fMCEventPos(NULL),
  fMCEventNeg(NULL),
  fESDArrayPos(NULL),
  fESDArrayNeg(NULL),
  fnCuts(0),
  fiCut(0),
  fMoveParticleAccordingToVertex(kTRUE),
  fIsHeavyIon(0),
  fDoLightOutput(0),
  fPi0EtaSwitch(0),
  fDoECalibOutput(kFALSE),
  fDoMesonAnalysis(kTRUE),
  fDoMesonQA(0),
  fDoPhotonQA(0),
  fDoClusterQA(0),
  fIsFromDesiredHeader(kTRUE),
  fIsOverlappingWithOtherHeader(kFALSE),
  fIsMC(0),
  fDoTHnSparse(kTRUE),
  fSetPlotHistsExtQA(kFALSE),
  fWeightJetJetMC(1),
  fDoConvGammaShowerShapeTree(kFALSE),
  fEnableSortForClusMC(kFALSE),
  fDoPrimaryTrackMatching(kFALSE),
  fDoInvMassShowerShapeTree(kFALSE),
  fDoMaterialBudgetWeightingOfGammasForTrueMesons(kFALSE),
  tBrokenFiles(NULL),
  fFileNameBroken(NULL),
  fAllowOverlapHeaders(kTRUE),
  fTrackMatcherRunningMode(0),
  fDoHBTHistoOutput(kFALSE)
{
  // Define output slots here
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());
  DefineOutput(4, TTree::Class());
  DefineOutput(5, TTree::Class());
  DefineOutput(6, TTree::Class());
  DefineOutput(7, TTree::Class());
  DefineOutput(8, TTree::Class());
  DefineOutput(9, TTree::Class());
}

AliAnalysisTaskGammaConvCalo::~AliAnalysisTaskGammaConvCalo()
{
  if(fGammaCandidates){
    delete fGammaCandidates;
    fGammaCandidates = 0x0;
  }
  if(fClusterCandidates){
    delete fClusterCandidates;
    fClusterCandidates = 0x0;
  }
  if(fBGHandler){
    delete[] fBGHandler;
    fBGHandler = 0x0;
  }
  if(fBGHandlerRP){
    delete[] fBGHandlerRP;
    fBGHandlerRP = 0x0;
  }
  if(fBGClusHandler){
    delete[] fBGClusHandler;
    fBGClusHandler = 0x0;
  }
  if(fBGClusHandlerRP){
    delete[] fBGClusHandlerRP;
    fBGClusHandlerRP = 0x0;
  }
  if(fBGHBTTrueGammaHandler){
    delete[] fBGHBTTrueGammaHandler;
    fBGHBTTrueGammaHandler = 0x0;
  }
  if(fBGHBTGenGammaHandler){
    delete[] fBGHBTGenGammaHandler;
    fBGHBTGenGammaHandler = 0x0;
  }
  if(fTrueGammaCandidatesConv){
    delete[] fTrueGammaCandidatesConv;
    fTrueGammaCandidatesConv = 0x0;
  }
  if(fTrueGammaCandidatesCalo){
    delete[] fTrueGammaCandidatesCalo;
    fTrueGammaCandidatesCalo = 0x0;
  }
  if(fMCGammaCandidates){
    delete[] fMCGammaCandidates;
    fMCGammaCandidates = 0x0;
  }
}
//___________________________________________________________
void AliAnalysisTaskGammaConvCalo::InitBack(){

  const Int_t nDim = 4;
  Int_t nBins[nDim] = {800,300,7,4};
  Double_t xMin[nDim] = {0,0, 0,0};
  Double_t xMax[nDim] = {0.8,30,7,4};

  if(fDoTHnSparse){
      fSparseMotherInvMassPtZM = new THnSparseF*[fnCuts];
      fSparseMotherBackInvMassPtZM = new THnSparseF*[fnCuts];
  }

  fBGHandler = new AliGammaConversionAODBGHandler*[fnCuts];
  fBGHandlerRP = new AliConversionAODBGHandlerRP*[fnCuts];

  fBGClusHandler = new AliGammaConversionAODBGHandler*[fnCuts];
  fBGClusHandlerRP = new AliConversionAODBGHandlerRP*[fnCuts];

  if(fIsMC>0 && fDoHBTHistoOutput){
    fBGHBTTrueGammaHandler  = new AliGammaConversionAODBGHandler*[fnCuts];
    fBGHBTGenGammaHandler   = new AliGammaConversionAODBGHandler*[fnCuts];
  }

  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if (((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->DoBGCalculation()){
      TString cutstringEvent   = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
      TString cutstringPhoton  = ((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutNumber();
      TString cutstringCalo   = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
      TString cutstringMeson   = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();

      Int_t collisionSystem = atoi((TString)(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber())(0,1));
      Int_t centMin = atoi((TString)(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber())(1,1));
      Int_t centMax = atoi((TString)(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber())(2,1));

      if(collisionSystem == 1 || collisionSystem == 2 ||
        collisionSystem == 5 || collisionSystem == 8 ||
        collisionSystem == 9){
        centMin = centMin*10;
        centMax = centMax*10;
        if(centMax ==0 && centMax!=centMin) centMax=100;
      }else if(collisionSystem == 3 || collisionSystem == 6){
        centMin = centMin*5;
        centMax = centMax*5;
      }else if(collisionSystem == 4 || collisionSystem == 7){
        centMin = ((centMin*5)+45);
        centMax = ((centMax*5)+45);
      }

      if(fDoTHnSparse){
        fBackList[iCut] = new TList();
        fBackList[iCut]->SetName(Form("%s_%s_%s_%s Back histograms",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringCalo.Data(), cutstringMeson.Data()));
        fBackList[iCut]->SetOwner(kTRUE);
        fCutFolder[iCut]->Add(fBackList[iCut]);

        fSparseMotherBackInvMassPtZM[iCut] = new THnSparseF("Back_Back_InvMass_Pt_z_m", "Back_Back_InvMass_Pt_z_m",nDim,nBins,xMin,xMax);
        fBackList[iCut]->Add(fSparseMotherBackInvMassPtZM[iCut]);

        fMotherList[iCut] = new TList();
        fMotherList[iCut]->SetName(Form("%s_%s_%s_%s Mother histograms",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
        fMotherList[iCut]->SetOwner(kTRUE);
        fCutFolder[iCut]->Add(fMotherList[iCut]);

        fSparseMotherInvMassPtZM[iCut] = new THnSparseF("Back_Mother_InvMass_Pt_z_m", "Back_Mother_InvMass_Pt_z_m",nDim,nBins,xMin,xMax);
        fMotherList[iCut]->Add(fSparseMotherInvMassPtZM[iCut]);
      }

      if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->BackgroundHandlerType() == 0){
        fBGHandler[iCut] = new AliGammaConversionAODBGHandler(
                                  collisionSystem,centMin,centMax,
                                  ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetNumberOfBGEvents(),
                                  ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseTrackMultiplicity(),
                                  2,8,5);
        fBGClusHandler[iCut] = new AliGammaConversionAODBGHandler(
                                  collisionSystem,centMin,centMax,
                                  ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetNumberOfBGEvents(),
                                  ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseTrackMultiplicity(),
                                  2,8,5);
        fBGHandlerRP[iCut] = NULL;
        if(fIsMC>0 && fDoHBTHistoOutput){
          fBGHBTTrueGammaHandler[iCut] = new AliGammaConversionAODBGHandler(
                                  collisionSystem,centMin,centMax,
                                  ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetNumberOfBGEvents(),
                                  ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseTrackMultiplicity(),
                                  2,8,5);
          fBGHBTGenGammaHandler[iCut]  = new AliGammaConversionAODBGHandler(
                                  collisionSystem,centMin,centMax,
                                  1,
                                  ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseTrackMultiplicity(),
                                  2,8,5);
        }
      }else if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->BackgroundHandlerType() != 2){
        fBGHandlerRP[iCut] = new AliConversionAODBGHandlerRP(
                                  ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsHeavyIon(),
                                  ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity(),
                                  ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetNumberOfBGEvents());
        fBGClusHandlerRP[iCut] = new AliConversionAODBGHandlerRP(
                                  ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsHeavyIon(),
                                  ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity(),
                                  ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetNumberOfBGEvents());
        fBGHandler[iCut] = NULL;
      }
    }
  }
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::UserCreateOutputObjects(){

  fV0Reader = (AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
  if(!fV0Reader){printf("Error: No V0 Reader");return;}// GetV0Reader

  if(fDoMesonAnalysis){ //Same Jet Finder MUST be used within same trainconfig
    if( ((AliConversionMesonCuts*)fMesonCutArray->At(0))->DoJetAnalysis())  fDoJetAnalysis = kTRUE;
    if( ((AliConversionMesonCuts*)fMesonCutArray->At(0))->DoJetQA())        fDoJetQA       = kTRUE;
  }

  if(fDoJetAnalysis){
    fConvJetReader=(AliAnalysisTaskConvJet*)AliAnalysisManager::GetAnalysisManager()->GetTask("AliAnalysisTaskConvJet");
    if(!fConvJetReader){printf("Error: No AliAnalysisTaskConvJet");return;} // GetV0Reader
  }
  if(((AliConvEventCuts*)fEventCutArray->At(0))->GetUseJetFinderForOutliers()){
    fOutlierJetReader=(AliAnalysisTaskJetOutlierRemoval*)AliAnalysisManager::GetAnalysisManager()->GetTask("AliAnalysisTaskJetOutlierRemoval");
    if(!fOutlierJetReader){AliFatal("Error: No AliAnalysisTaskJetOutlierRemoval");} // Get jet outlier task
    else{printf("Found AliAnalysisTaskJetOutlierRemoval used for outlier removal!\n");}
  }
  if (fIsMC == 2){
    fDoClusterQA      = 0;
    fDoTHnSparse      = kFALSE;
  } else if (fIsMC == 3){
    fDoTHnSparse      = kFALSE;
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
  fGammaCandidates    = new TList();
  fClusterCandidates  = new TList();
  fClusterCandidates->SetOwner(kTRUE);

  fCutFolder          = new TList*[fnCuts];
  fESDList            = new TList*[fnCuts];

  if(fDoTHnSparse){
    fBackList         = new TList*[fnCuts];
    fMotherList       = new TList*[fnCuts];
  }

  fHistoNEvents               = new TH1F*[fnCuts];
  if (fIsMC > 1){
    fHistoNEventsWOWeight     = new TH1F*[fnCuts];
  }
  if (fIsMC == 2){
    fProfileJetJetXSection    = new TProfile*[fnCuts];
    fHistoJetJetNTrials       = new TH1F*[fnCuts];
  }

  if(fDoLightOutput!=2){
    fHistoNGoodESDTracks        = new TH1F*[fnCuts];
    fHistoVertexZ               = new TH1F*[fnCuts];
    fHistoNGammaCandidates      = new TH1F*[fnCuts];
  }
  if(!fDoLightOutput){
    fHistoVertexX               = new TH1F*[fnCuts];
    fHistoVertexY               = new TH1F*[fnCuts];
  }
  fHistoEventSphericity       = new TH1F*[fnCuts];

  if(fIsHeavyIon==2) fProfileEtaShift = new TProfile*[fnCuts];
  if(!fDoLightOutput){
    fHistoNGoodESDTracksVsNGammaCandidates  = new TH2F*[fnCuts];
    fHistoSPDClusterTrackletBackground      = new TH2F*[fnCuts];
    fHistoNV0Tracks                         = new TH1F*[fnCuts];
    fHistoConvGammaPt                       = new TH1F*[fnCuts];
  }

  if (fDoPhotonQA == 2){
    fTreeConvGammaPtDcazCat   = new TTree*[fnCuts];
  }

  if(fDoMesonAnalysis){
    fHistoMotherInvMassPt             = new TH2F*[fnCuts];
    fHistoMotherBackInvMassPt         = new TH2F*[fnCuts];
    if(!fDoLightOutput){
      fHistoMotherMatchedInvMassPt      = new TH2F*[fnCuts];
      fHistoMotherInvMassPtAlpha        = new TH2F*[fnCuts];
      fHistoPhotonPairPtconv            = new TH2F*[fnCuts];
      fHistoPhotonPairMixedEventPtconv  = new TH2F*[fnCuts];
    }
    if(!fDoLightOutput || fDoECalibOutput){
        fHistoMotherInvMassECalib         = new TH2F*[fnCuts];
        fHistoMotherBackInvMassECalib     = new TH2F*[fnCuts];
        fHistoMotherInvMassECalibPCM         = new TH2F*[fnCuts];
        fHistoMotherBackInvMassECalibPCM     = new TH2F*[fnCuts];
    }
    if (fDoMesonQA > 0){
      fHistoMotherPi0PtY              = new TH2F*[fnCuts];
      fHistoMotherEtaPtY              = new TH2F*[fnCuts];
      fHistoMotherPi0PtAlpha          = new TH2F*[fnCuts];
      fHistoMotherEtaPtAlpha          = new TH2F*[fnCuts];
      fHistoMotherPi0PtOpenAngle      = new TH2F*[fnCuts];
      fHistoMotherEtaPtOpenAngle      = new TH2F*[fnCuts];
      fHistoMotherPi0ConvPhotonEtaPhi = new TH2F*[fnCuts];
      fHistoMotherEtaConvPhotonEtaPhi = new TH2F*[fnCuts];
    }
    if(fDoMesonQA == 2){
      tESDMesonsInvMassPtDcazMinDcazMaxFlag = new TTree*[fnCuts];
    }
  }
  if(fDoJetAnalysis){
    fJetHistograms            = new TList*[fnCuts];

    fHistoPtJet               = new TH1F*[fnCuts];
    fHistoJetEta              = new TH1F*[fnCuts];
    fHistoJetPhi              = new TH1F*[fnCuts];
    fHistoJetArea             = new TH1F*[fnCuts];
    fHistoNJets               = new TH1F*[fnCuts];
    fHistoEventwJets          = new TH1F*[fnCuts];
    if(!fDoLightOutput){
      fHistoJetPi0PtRatio       = new TH1F*[fnCuts];
      fHistoDoubleCounting      = new TH1F*[fnCuts];
      fHistoJetMotherInvMassPt                  = new TH2F*[fnCuts];
      fHistoPi0InJetMotherInvMassPt             = new TH2F*[fnCuts];
      fHistoMotherBackJetInvMassPt              = new TH2F*[fnCuts];
      fHistoRJetPi0Cand                         = new TH2F*[fnCuts];
      fHistoEtaPhiJetPi0Cand                    = new TH2F*[fnCuts];
      fHistoEtaPhiJetWithPi0Cand                = new TH2F*[fnCuts];
      fHistoJetFragmFunc                        = new TH2F*[fnCuts];
      fHistoJetFragmFuncZInvMass                = new TH2F*[fnCuts];
    }
  }
  if(fDoHBTHistoOutput){
    fHistoHBTOpeningAnglePt             = new TH2F*[fnCuts];
    fHistoBckHBTOpeningAnglePt          = new TH2F*[fnCuts];
    fHistoHBTDeltaEPt                   = new TH2F*[fnCuts];
    fHistoBckHBTDeltaEPt                = new TH2F*[fnCuts];
  }

  if(fIsMC>0 && fDoHBTHistoOutput){
    fTrueGammaCandidatesConv    = new TList();
    fTrueGammaCandidatesCalo    = new TList();
    fMCGammaCandidates      = new TList();
  }

  fClusterOutputList                  = new TList*[fnCuts];
  fHistoClusGammaPt                   = new TH1F*[fnCuts];
  fHistoClusGammaE                    = new TH1F*[fnCuts];

  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetClusterType()==2){
      if (fHistoClusGammaE_BothBM_highestE==NULL){
        fHistoClusGammaE_BothBM_highestE = new TH1F*[fnCuts];
      }
      //--------------------------------------------------
      //Only MB
      //if ( (((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger()==1)||(((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger()==0) ){
        if (fHistoClusGammaE_BothBM==NULL){
          fHistoClusGammaE_BothBM = new TH1F*[fnCuts];
        }
        if (fHistoClusGammaE_AnaBM_highestE==NULL){
          fHistoClusGammaE_AnaBM_highestE = new TH1F*[fnCuts];
        }
      //}
      //--------------------------------------------------
      //Only PHI7
      if ( ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger()==6 ){
        if (fHistoClusGammaE_onlyTriggered==NULL) {
          fHistoClusGammaE_onlyTriggered = new TH1F*[fnCuts];
        }
      }
    }
  }

  fHistoClusOverlapHeadersGammaPt     = new TH1F*[fnCuts];
  fHistoClusAllHeadersGammaPt         = new TH1F*[fnCuts];
  fHistoClusRejectedHeadersGammaPt    = new TH1F*[fnCuts];
  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetClusterType()==2){
      if (fCaloTriggerMimicHelper == NULL){
        fCaloTriggerMimicHelper     = new AliCaloTriggerMimicHelper*[fnCuts];
      }
      if ( ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger()==6 ){
        if (fHistoGoodMesonClusters == NULL){
          fHistoGoodMesonClusters     = new TH1I*[fnCuts];
        }
      }
    }
  }

  if(fDoConvGammaShowerShapeTree){
    fGammaERM02               = new TList*[fnCuts];
    tESDGammaERM02            = new TTree*[fnCuts];
  }

  if(fDoInvMassShowerShapeTree){
    fInvMassShowerShape       = new TList*[fnCuts];
    tESDInvMassShowerShape    = new TTree*[fnCuts];
  }

  // set common binning in pT for mesons and photons
  Int_t nBinsPt               = 200;
  Float_t binWidthPt          = 0.1;
  Float_t minPt               = 0;
  Float_t maxPt               = 20;
  Int_t nBinsQAPt             = 170;
  Float_t maxQAPt             = 20;
  Int_t nBinsClusterPt        = 500;
  Float_t minClusterPt        = 0;
  Float_t maxClusterPt        = 50;
  Int_t nBinsMinv             = 800;
  Float_t minMinv             = 0;
  Float_t maxMinv             = 0.8;
  Double_t *arrPtBinning      = new Double_t[1200];
  Double_t *arrQAPtBinning    = new Double_t[1200];
  Double_t *arrClusPtBinning  = new Double_t[1200];

  // Set special pt binning for pp 8TeV
  if (((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::kpPb5TeV ||
      ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::kpPb5TeVR2 ||
      ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::k5TeV ||
      ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::k7TeV ||
      ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::k8TeV ||
      ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::kpPb8TeV){
    nBinsMinv                 = 400;
    nBinsPt                   = 179;
    minPt                     = 0;
    maxPt                     = 60;
    for(Int_t i=0; i<nBinsPt+1;i++){
      if (i < 100) arrPtBinning[i]            = 0.10*i;
      else if(i<140) arrPtBinning[i]          = 10.+0.25*(i-100);
      else if(i<180) arrPtBinning[i]          = 20.+1.00*(i-140);
      else arrPtBinning[i]                    = maxPt;
    }
    nBinsQAPt                 = 179;
    maxQAPt                   = 60;
    for(Int_t i=0; i<nBinsQAPt+1;i++){
      if (i < 100) arrQAPtBinning[i]            = 0.10*i;
      else if(i<140) arrQAPtBinning[i]          = 10.+0.25*(i-100);
      else if(i<180) arrQAPtBinning[i]          = 20.+1.00*(i-140);
      else arrQAPtBinning[i]                    = maxQAPt;
    }
    nBinsClusterPt            = 301;
    minClusterPt              = 0;
    maxClusterPt              = 100;
    for(Int_t i=0; i<nBinsClusterPt+1;i++){
      if (i < 1) arrClusPtBinning[i]          = 0.3*i;
      else if(i<55) arrClusPtBinning[i]       = 0.3+0.05*(i-1);
      else if(i<125) arrClusPtBinning[i]      = 3.+0.1*(i-55);
      else if(i<155) arrClusPtBinning[i]      = 10.+0.2*(i-125);
      else if(i<211) arrClusPtBinning[i]      = 16.+0.25*(i-155);
      else if(i<251) arrClusPtBinning[i]      = 30.+0.5*(i-211);
      else if(i<301) arrClusPtBinning[i]      = 50.+1.0*(i-251);
      else arrClusPtBinning[i]                = maxClusterPt;
    }
  // Set special pt binning for pp 13TeV, pPb 8TeV
  } else if ( ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::k13TeV ||
              ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::k13TeVLowB ){

    nBinsMinv = 400;
    maxMinv   = 0.8;


    nBinsPt                   = 310;
    minPt                     = 0;
    maxPt                     = 100;
    binWidthPt                = 0.1;
    for(Int_t i=0; i<nBinsPt+1;i++){
      if (i < 1) arrPtBinning[i]              = 0.3*i;
      else if(i<198) arrPtBinning[i]          = 0.3+0.1*(i-1);
      else if(i<238) arrPtBinning[i]          = 20.+0.25*(i-198);
      else if(i<278) arrPtBinning[i]          = 30.+0.5*(i-238);
      else if(i<298) arrPtBinning[i]          = 50.+1.0*(i-278);
      else if(i<310) arrPtBinning[i]          = 70.+2.5*(i-298);
      else  arrPtBinning[i]                   = maxPt;
    }
    nBinsQAPt                 = 240;
    maxQAPt                   = 100;
    for(Int_t i=0; i<nBinsQAPt+1;i++){
      if(i<100) arrQAPtBinning[i]             = 0.1*i;
      else if(i<140) arrQAPtBinning[i]        = 10.+0.25*(i-100);
      else if(i<180) arrQAPtBinning[i]        = 20.+0.5*(i-140);
      else if(i<240) arrQAPtBinning[i]        = 40.+1.0*(i-180);
      else arrQAPtBinning[i]                  = maxQAPt;
    }
    nBinsClusterPt            = 310;
    minClusterPt              = 0;
    maxClusterPt              = 100;
    for(Int_t i=0; i<nBinsClusterPt+1;i++){
      if (i < 1) arrClusPtBinning[i]          = 0.3*i;
      else if(i<198) arrClusPtBinning[i]      = 0.3+0.1*(i-1);
      else if(i<238) arrClusPtBinning[i]      = 20.+0.25*(i-198);
      else if(i<278) arrClusPtBinning[i]      = 30.+0.5*(i-238);
      else if(i<298) arrClusPtBinning[i]      = 50.+1.0*(i-278);
      else if(i<310) arrClusPtBinning[i]      = 70.+2.5*(i-298);
      else arrClusPtBinning[i]                = maxClusterPt;
    }
  } else if (((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::kXeXe5440GeV  ){
    nBinsPt                   = 90;
    minPt                     = 0;
    maxPt                     = 20;
    for(Int_t i=0; i<nBinsPt+1;i++){
      if (i < 1) arrPtBinning[i]              = 0.3*i;
      else if(i<58) arrPtBinning[i]           = 0.3+0.1*(i-1);
      else if(i<82) arrPtBinning[i]           = 6.+0.25*(i-58);
      else if(i<90) arrPtBinning[i]           = 12.+1.0*(i-82);
      else arrPtBinning[i]                    = maxPt;
    }
    nBinsQAPt                 = 92;
    maxQAPt                   = 20;
    for(Int_t i=0; i<nBinsQAPt+1;i++){
      if(i<60) arrQAPtBinning[i]              = 0.1*i;
      else if(i<84) arrQAPtBinning[i]         = 6.+0.25*(i-60);
      else if(i<92) arrQAPtBinning[i]         = 12.+1.0*(i-84);
      else arrQAPtBinning[i]                  = maxQAPt;
    }
    nBinsClusterPt            = 148;
    minClusterPt              = 0;
    maxClusterPt              = 40;
    for(Int_t i=0; i<nBinsClusterPt+1;i++){
      if (i < 1) arrClusPtBinning[i]          = 0.3*i;
      else if(i<98) arrClusPtBinning[i]       = 0.3+0.1*(i-1);
      else if(i<123) arrClusPtBinning[i]      = 10.+0.2*(i-98);
      else if(i<148) arrClusPtBinning[i]      = 15.+1.0*(i-123);
      else arrClusPtBinning[i]                = maxClusterPt;
    }
  } else if (((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::kPbPb5TeV  ){
    nBinsPt                   = 110;
    minPt                     = 0;
    maxPt                     = 40;
    for(Int_t i=0; i<nBinsPt+1;i++){
      if (i < 1) arrPtBinning[i]              = 0.3*i;
      else if(i<58) arrPtBinning[i]           = 0.3+0.1*(i-1);
      else if(i<82) arrPtBinning[i]           = 6.+0.25*(i-58);
      else if(i<110) arrPtBinning[i]          = 12.+1.0*(i-82);
      else arrPtBinning[i]                    = maxPt;
    }
    nBinsQAPt                 = 92;
    maxQAPt                   = 20;
    for(Int_t i=0; i<nBinsQAPt+1;i++){
      if(i<60) arrQAPtBinning[i]              = 0.1*i;
      else if(i<84) arrQAPtBinning[i]         = 6.+0.25*(i-60);
      else if(i<92) arrQAPtBinning[i]         = 12.+1.0*(i-84);
      else arrQAPtBinning[i]                  = maxQAPt;
    }
    nBinsClusterPt            = 148;
    minClusterPt              = 0;
    maxClusterPt              = 40;
    for(Int_t i=0; i<nBinsClusterPt+1;i++){
      if (i < 1) arrClusPtBinning[i]          = 0.3*i;
      else if(i<98) arrClusPtBinning[i]       = 0.3+0.1*(i-1);
      else if(i<123) arrClusPtBinning[i]      = 10.+0.2*(i-98);
      else if(i<148) arrClusPtBinning[i]      = 15.+1.0*(i-123);
      else arrClusPtBinning[i]                = maxClusterPt;
    }
    // default binning
  } else {
    for(Int_t i=0; i<nBinsPt+1;i++){
      arrPtBinning[i]         = ((maxPt-minPt)/nBinsPt)*i;
    }
    for(Int_t i=0; i<nBinsClusterPt+1;i++){
      arrClusPtBinning[i]     = ((maxClusterPt-minClusterPt)/nBinsClusterPt)*i;
    }
    for(Int_t i=0; i<nBinsQAPt+1;i++){
      if(i<60) arrQAPtBinning[i]              = 0.05*i;
      else if(i<130) arrQAPtBinning[i]        = 3.+0.1*(i-60);
      else if(i<170) arrQAPtBinning[i]        = 10.+0.25*(i-130);
      else arrQAPtBinning[i]                  = maxQAPt;
    }
  }
   Double_t* arrLogBinning = new Double_t[51]{1.00e-03, 1.20e-03, 1.45e-03, 1.74e-03, 2.01e-03, 2.51e-03, 3.02e-03, 3.63e-03, 4.37e-03, 5.25e-03, 6.31e-03, 7.60e-03,
                                     9.12e-03, 1.01e-02, 1.32e-02, 1.58e-02, 1.91e-02, 2.29e-02, 2.75e-02, 3.31e-02, 3.98e-02, 4.79e-02, 5.75e-02, 6.91e-02,
                                     8.32e-02, 1.00e-01, 1.20e-01, 1.45e-01, 1.74e-01, 2.09e-01, 2.51e-01, 3.02e-01, 3.63e-01, 4.37e-01, 5.25e-01, 6.31e-01,
                                     7.59e-01, 9.12e-01, 1.10e+00, 1.32e+00, 1.58e+00, 1.91e+00, 2.29e+00, 2.75e+00, 3.31e+00, 3.98e+00, 4.79e+00, 5.75e+00,
                                     6.92e+00, 8.32e+00, 1.00e+01};

  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    TString cutstringEvent    = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
    TString cutstringPhoton   = ((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutNumber();
    TString cutstringCalo     = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
    TString cutstringMeson    = "NoMesonCut";
    if(fDoMesonAnalysis)
      cutstringMeson          = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();

    fCutFolder[iCut]          = new TList();
    fCutFolder[iCut]->SetName(Form("Cut Number %s_%s_%s_%s",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
    fCutFolder[iCut]->SetOwner(kTRUE);
    fOutputContainer->Add(fCutFolder[iCut]);
    fESDList[iCut]            = new TList();
    fESDList[iCut]->SetName(Form("%s_%s_%s_%s ESD histograms",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
    fESDList[iCut]->SetOwner(kTRUE);
    fCutFolder[iCut]->Add(fESDList[iCut]);

    fHistoNEvents[iCut]       = new TH1F("NEvents", "NEvents", 14, -0.5, 13.5);
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(3,"Miss. MC or inc. ev.");
    if (((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger() > 1 ){
      TString TriggerNames    = "Not Trigger: ";
      TriggerNames            = TriggerNames+ ( (AliConvEventCuts*)fEventCutArray->At(iCut))->GetSpecialTriggerName();
      fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(4,TriggerNames.Data());
    }else {
      fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(4,"Trigger");
    }
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(5,"Vertex Z");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(7,"Pile-Up");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(8,"no SDD");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(9,"no V0AND");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(10,"EMCAL/TPC problem");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(11,"rejectedForJetJetMC");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(14,"Pileup V0M-TPCout Tracks");
    fESDList[iCut]->Add(fHistoNEvents[iCut]);

    if (fIsMC > 1){
      fHistoNEventsWOWeight[iCut]   = new TH1F("NEventsWOWeight", "NEventsWOWeight", 14, -0.5, 13.5);
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(3,"Miss. MC or inc. ev.");
      if (((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger() > 1 ){
        TString TriggerNames        = "Not Trigger: ";
        TriggerNames                = TriggerNames+ ( (AliConvEventCuts*)fEventCutArray->At(iCut))->GetSpecialTriggerName();
        fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(4,TriggerNames.Data());
      }else {
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
      fProfileJetJetXSection[iCut]  = new TProfile("XSection", "XSection", 1, -0.5, 0.5);
      fESDList[iCut]->Add(fProfileJetJetXSection[iCut]);
      fHistoJetJetNTrials[iCut]     = new TH1F("NTrials", "#sum{NTrials}", 1, 0, 1);
      fHistoJetJetNTrials[iCut]->GetXaxis()->SetBinLabel(1,"#sum{NTrials}");
      fESDList[iCut]->Add(fHistoJetJetNTrials[iCut]);
    }

    if(fDoLightOutput!=2){
      if(fIsHeavyIon == 1)
      fHistoNGoodESDTracks[iCut]    = new TH1F("GoodESDTracks", "GoodESDTracks", 4000, 0, 4000);
      else if(fIsHeavyIon == 2)
      fHistoNGoodESDTracks[iCut]    = new TH1F("GoodESDTracks", "GoodESDTracks", 400, 0, 400);
      else
      fHistoNGoodESDTracks[iCut]    = new TH1F("GoodESDTracks", "GoodESDTracks", 200, 0, 200);
      fHistoNGoodESDTracks[iCut]->SetXTitle("# TPC tracks");
      fESDList[iCut]->Add(fHistoNGoodESDTracks[iCut]);
    }

    if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetUseSphericity()!=0){
      fHistoEventSphericity[iCut]     = new TH1F("EventSphericity", "EventSphericity", 100, 0, 1);
      fHistoEventSphericity[iCut]->GetXaxis()->SetTitle("S");
      fESDList[iCut]->Add(fHistoEventSphericity[iCut]);
      fV0Reader->SetCalcSphericity(kTRUE);
    }
    if(fDoLightOutput!=2){
        fHistoVertexZ[iCut]             = new TH1F("VertexZ", "VertexZ", 200, -10, 10);
        fESDList[iCut]->Add(fHistoVertexZ[iCut]);
    }
    if(!fDoLightOutput){
      fHistoVertexX[iCut]             = new TH1F("VertexX", "VertexX", 100, -0.5, 0.5);
      fESDList[iCut]->Add(fHistoVertexX[iCut]);
      fHistoVertexY[iCut]             = new TH1F("VertexY", "VertexY", 100, -0.1, 0.9);
      fESDList[iCut]->Add(fHistoVertexY[iCut]);
    }

    if(fDoLightOutput!=2){
        if(fIsHeavyIon == 1)
          fHistoNGammaCandidates[iCut]  = new TH1F("GammaCandidates", "GammaCandidates", 150, 0, 150);
        else if(fIsHeavyIon == 2)
          fHistoNGammaCandidates[iCut]  = new TH1F("GammaCandidates", "GammaCandidates", 50, 0, 50);
        else
          fHistoNGammaCandidates[iCut]  = new TH1F("GammaCandidates", "GammaCandidates", 50, 0, 50);
        fHistoNGammaCandidates[iCut]->SetXTitle("# accepted #gamma_{conv}");
        fESDList[iCut]->Add(fHistoNGammaCandidates[iCut]);
    }

    if(!fDoLightOutput){
      if(fIsHeavyIon == 1)
        fHistoNGoodESDTracksVsNGammaCandidates[iCut]  = new TH2F("GoodESDTracksVsGammaCandidates", "GoodESDTracksVsGammaCandidates", 4000, 0, 4000, 100, 0, 100);
      else if(fIsHeavyIon == 2)
        fHistoNGoodESDTracksVsNGammaCandidates[iCut]  = new TH2F("GoodESDTracksVsGammaCandidates", "GoodESDTracksVsGammaCandidates", 400, 0, 400, 50, 0, 50);
      else
        fHistoNGoodESDTracksVsNGammaCandidates[iCut]  = new TH2F("GoodESDTracksVsGammaCandidates", "GoodESDTracksVsGammaCandidates", 200, 0, 200, 50, 0, 50);
      fHistoNGoodESDTracksVsNGammaCandidates[iCut]->SetXTitle("# TPC tracks");
      fHistoNGoodESDTracksVsNGammaCandidates[iCut]->SetYTitle("# accepted #gamma_{conv}");
      fESDList[iCut]->Add(fHistoNGoodESDTracksVsNGammaCandidates[iCut]);

      fHistoSPDClusterTrackletBackground[iCut]        = new TH2F("SPD tracklets vs SPD clusters", "SPD tracklets vs SPD clusters", 100, 0, 200, 250, 0, 1000);
      fESDList[iCut]->Add(fHistoSPDClusterTrackletBackground[iCut]);

      if(fIsHeavyIon == 1)
        fHistoNV0Tracks[iCut]         = new TH1F("V0 Multiplicity", "V0 Multiplicity", 30000, 0, 30000);
      else if(fIsHeavyIon == 2)
        fHistoNV0Tracks[iCut]         = new TH1F("V0 Multiplicity", "V0 Multiplicity", 2500, 0, 2500);
      else
        fHistoNV0Tracks[iCut]         = new TH1F("V0 Multiplicity", "V0 Multiplicity", 1500, 0, 1500);
      fHistoNV0Tracks[iCut]->SetXTitle("VZERO amp [arb. units]");
      fESDList[iCut]->Add(fHistoNV0Tracks[iCut]);

      fHistoConvGammaPt[iCut]         = new TH1F("ESD_ConvGamma_Pt", "ESD_ConvGamma_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
      fHistoConvGammaPt[iCut]->SetXTitle("p_{T,conv}(GeV/c)");
      fESDList[iCut]->Add(fHistoConvGammaPt[iCut]);
    }

    if(fIsHeavyIon == 2){
      fProfileEtaShift[iCut]          = new TProfile("Eta Shift", "Eta Shift", 1, -0.5, 0.5);
      fESDList[iCut]->Add(fProfileEtaShift[iCut]);
    }

    if (fIsMC > 1){
      if(fHistoNEvents[iCut]) fHistoNEvents[iCut]->Sumw2();
      if(fHistoNGoodESDTracks && fHistoNGoodESDTracks[iCut]) fHistoNGoodESDTracks[iCut]->Sumw2();
      if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetUseSphericity()!=0){
        fHistoEventSphericity[iCut]->Sumw2();
      }
      if(fDoLightOutput!=2) {
        fHistoVertexZ[iCut]->Sumw2();
        fHistoNGammaCandidates[iCut]->Sumw2();
      }
      if(!fDoLightOutput){
        fHistoVertexX[iCut]->Sumw2();
        fHistoVertexY[iCut]->Sumw2();
        fHistoNGoodESDTracksVsNGammaCandidates[iCut]->Sumw2();
        fHistoSPDClusterTrackletBackground[iCut]->Sumw2();
        fHistoNV0Tracks[iCut]->Sumw2();
        fHistoConvGammaPt[iCut]->Sumw2();
      }
    }

    if (  (fIsMC>0 && fDoMaterialBudgetWeightingOfGammasForTrueMesons) ) {
      if(!fDoLightOutput){
	    fHistoConvGammaPt[iCut]->Sumw2();
      }
    }

    if (fDoPhotonQA == 2 ){
      fTreeConvGammaPtDcazCat[iCut] = new TTree(Form("%s_%s_%s_%s Photon DCA tree",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringCalo.Data(),cutstringMeson.Data()), "ESD_ConvGamma_Pt_Dcaz_Cat");
      fTreeConvGammaPtDcazCat[iCut]->Branch("Pt",&fPtGamma,"fPtGamma/s");
      fTreeConvGammaPtDcazCat[iCut]->Branch("DcaZPhoton",&fDCAzPhoton,"fDCAzPhoton/S");
      fTreeConvGammaPtDcazCat[iCut]->Branch("cat",&fCharCatPhoton,"fCharCatPhoton/b");
      if (fIsMC > 1){
        fTreeConvGammaPtDcazCat[iCut]->Branch("weightEvent",&fWeightJetJetMC,"fWeightJetJetMC/f");
      }
    }

    fClusterOutputList[iCut]        = new TList();
    fClusterOutputList[iCut]->SetName(Form("%s_%s_%s_%s Cluster Output",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
    fClusterOutputList[iCut]->SetOwner(1);
    fCutFolder[iCut]->Add(fClusterOutputList[iCut]);

    fHistoClusGammaPt[iCut]         = new TH1F("ClusGamma_Pt", "ClusGamma_Pt", nBinsClusterPt, arrClusPtBinning);
    fHistoClusGammaPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
    fClusterOutputList[iCut]->Add(fHistoClusGammaPt[iCut]);
    fHistoClusGammaE[iCut]          = new TH1F("ClusGamma_E", "ClusGamma_E", nBinsClusterPt, arrClusPtBinning);
    fHistoClusGammaE[iCut]->SetXTitle("E_{clus} (GeV)");
    fClusterOutputList[iCut]->Add(fHistoClusGammaE[iCut]);
    if(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetClusterType()==2){
      if ( ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger()==6 ){
        fHistoGoodMesonClusters[iCut]     = new TH1I( "fHistoGoodMesonClusters", "fHistoGoodMesonClusters", 7, 0.5, 7.5);
        fHistoGoodMesonClusters[iCut]->GetXaxis()->SetBinLabel(1,"All Meson Candidates Candidates");
        fHistoGoodMesonClusters[iCut]->GetXaxis()->SetBinLabel(2,"Triggered Meson Candidates");
        fHistoGoodMesonClusters[iCut]->GetXaxis()->SetBinLabel(3,"Cluster Not Triggered");
        fHistoGoodMesonClusters[iCut]->GetXaxis()->SetBinLabel(4,"Cluster E passed");
        fHistoGoodMesonClusters[iCut]->GetXaxis()->SetBinLabel(5,"Cluster E not passed");
        Bool_t FlagMaybeBadDDLs=((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetReduceTriggeredPhiDueBadDDLs();
        if (FlagMaybeBadDDLs == kTRUE){
            fHistoGoodMesonClusters[iCut]->GetXaxis()->SetBinLabel(6,"DDL passed (may.bad)");
            fHistoGoodMesonClusters[iCut]->GetXaxis()->SetBinLabel(7,"DDL not passed (may.bad)");
        } else {
            fHistoGoodMesonClusters[iCut]->GetXaxis()->SetBinLabel(6,"DDL passed");
            fHistoGoodMesonClusters[iCut]->GetXaxis()->SetBinLabel(7,"DDL not passed");
        }
        fESDList[iCut]->Add(fHistoGoodMesonClusters[iCut]);
      }
    }
    if(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetClusterType()==2){
      fHistoClusGammaE_BothBM_highestE[iCut] = new TH1F("HistoClusGammaE_BothBM_highestE", "HistoClusGammaE_BothBM_highestE", nBinsClusterPt, arrClusPtBinning);
      fHistoClusGammaE_BothBM_highestE[iCut]->SetXTitle("E_{clus} (GeV/c)");
      fESDList[iCut]->Add(fHistoClusGammaE_BothBM_highestE[iCut]);
      //--------------------------------------------------
      //Only MB
      //if ( (((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger()==1)||(((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger()==0) ){
        fHistoClusGammaE_BothBM[iCut] = new TH1F("HistoClusGammaE_BothBM", "HistoClusGammaE_BothBM", nBinsClusterPt, arrClusPtBinning);
        fHistoClusGammaE_BothBM[iCut]->SetXTitle("E_{clus} (GeV/c)");
        fESDList[iCut]->Add(fHistoClusGammaE_BothBM[iCut]);
        //----------
        fHistoClusGammaE_AnaBM_highestE[iCut] = new TH1F("HistoClusGammaE_AnaBM_highestE", "HistoClusGammaE_AnaBM_highestE", nBinsClusterPt, arrClusPtBinning);
        fHistoClusGammaE_AnaBM_highestE[iCut]->SetXTitle("E_{clus} (GeV/c)");
        fESDList[iCut]->Add(fHistoClusGammaE_AnaBM_highestE[iCut]);
      //}
      //--------------------------------------------------
      //Only PHI7
      if ( ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger()==6 ){
        fHistoClusGammaE_onlyTriggered[iCut] = new TH1F("ClusGamma_E_onlyTriggered", "ClusGamma_E_onlyTriggered", nBinsClusterPt, arrClusPtBinning);
        fHistoClusGammaE_onlyTriggered[iCut]->SetXTitle("E_{clus} (GeV/c)");
        fESDList[iCut]->Add(fHistoClusGammaE_onlyTriggered[iCut]);
      }
    }
    fHistoClusOverlapHeadersGammaPt[iCut]   = new TH1F("ClusGammaOverlapHeaders_Pt", "ClusGammaOverlapHeaders_Pt", nBinsClusterPt, arrClusPtBinning);
    fHistoClusOverlapHeadersGammaPt[iCut]->SetXTitle("p_{T,clus} (GeV/c), selected header w/ overlap");
    fClusterOutputList[iCut]->Add(fHistoClusOverlapHeadersGammaPt[iCut]);
    fHistoClusAllHeadersGammaPt[iCut]       = new TH1F("ClusGammaAllHeaders_Pt", "ClusGammaAllHeaders_Pt", nBinsClusterPt, arrClusPtBinning);
    fHistoClusAllHeadersGammaPt[iCut]->SetXTitle("p_{T,clus} (GeV/c), all headers");
    fClusterOutputList[iCut]->Add(fHistoClusAllHeadersGammaPt[iCut]);
    fHistoClusRejectedHeadersGammaPt[iCut]  = new TH1F("ClusGammaRejectedHeaders_Pt", "ClusGammaRejectedHeaders_Pt", nBinsClusterPt, arrClusPtBinning);
    fHistoClusRejectedHeadersGammaPt[iCut]->SetXTitle("p_{T,clus} (GeV/c), rejected headers");
    fClusterOutputList[iCut]->Add(fHistoClusRejectedHeadersGammaPt[iCut]);


    if(fDoConvGammaShowerShapeTree){
      fGammaERM02[iCut]           = new TList();
      fGammaERM02[iCut]->SetName(Form("%s_%s_%s_%s ConvGamma-Cluster Matched",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
      fGammaERM02[iCut]->SetOwner(kTRUE);
      fCutFolder[iCut]->Add(fGammaERM02[iCut]);

      tESDGammaERM02[iCut]        = new TTree("ESD_ConvGamma_E_ConvR_M02_M20", "ESD_ConvGamma_E_ConvR_M02_M20");
      tESDGammaERM02[iCut]->Branch("ClusterE",&tESDClusE,"tESDClusE/F");
      tESDGammaERM02[iCut]->Branch("ConvR",&tESDGammaConvR,"tESDGammaConvR/F");
      tESDGammaERM02[iCut]->Branch("M02",&tESDClusterM02,"tESDClusterM02/F");
      tESDGammaERM02[iCut]->Branch("M20",&tESDClusterM20,"tESDClusterM20/F");
      tESDGammaERM02[iCut]->Branch("Eta",&tESDClusterEta,"tESDClusterEta/F");
      tESDGammaERM02[iCut]->Branch("Phi",&tESDClusterPhi,"tESDClusterPhi/F");
      tESDGammaERM02[iCut]->Branch("NCells",&tESDClusterNCells,"tESDClusterNCells/F");
      tESDGammaERM02[iCut]->Branch("MaxECell",&tESDClusterMaxECell,"tESDClusterMaxECell/F");
      tESDGammaERM02[iCut]->Branch("NLM",&tESDClusterNLM,"tESDClusterNLM/F");
      fGammaERM02[iCut]->Add(tESDGammaERM02[iCut]);
    }

    if(fDoInvMassShowerShapeTree){
      fInvMassShowerShape[iCut]           = new TList();
      fInvMassShowerShape[iCut]->SetName(Form("%s_%s_%s_%s InvMass_ShowerShape",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
      fInvMassShowerShape[iCut]->SetOwner(kTRUE);
      fCutFolder[iCut]->Add(fInvMassShowerShape[iCut]);
      tESDInvMassShowerShape[iCut] = new TTree("ESD_Meson_InvMass_Pt_ClusE_ClusM02_ClusM20", "ESD_Meson_InvMass_Pt_ClusE_ClusM02_ClusM20");
      tESDInvMassShowerShape[iCut]->Branch("MesonInvMass",&tESDIMMesonInvMass,"tESDIMMesonInvMass/F");
      tESDInvMassShowerShape[iCut]->Branch("MesonPt",&tESDIMMesonPt,"tESDIMMesonPt/F");
      tESDInvMassShowerShape[iCut]->Branch("ClusE",&tESDIMClusE,"tESDIMClusE/F");
      tESDInvMassShowerShape[iCut]->Branch("ClusM02",&tESDIMClusterM02,"tESDIMClusterM02/F");
      tESDInvMassShowerShape[iCut]->Branch("ClusM20",&tESDIMClusterM20,"tESDIMClusterM20/F");
      tESDInvMassShowerShape[iCut]->Branch("ClusLeadCellID",&tESDIMClusterLeadCellID,"tESDIMClusterLeadCellID/I");
      if(fIsMC>0) tESDInvMassShowerShape[iCut]->Branch("ClusClassification",&tESDIMClusterClassification,"tESDIMClusterClassification/I");
      tESDInvMassShowerShape[iCut]->Branch("ClusMatchedTrackPt",&tESDIMClusMatchedTrackPt,"tESDIMClusMatchedTrackPt/F");
      tESDInvMassShowerShape[iCut]->Branch("ClusTrackDeltaEta",&tESDIMClusTrackDeltaEta,"tESDIMClusTrackDeltaEta/F");
      tESDInvMassShowerShape[iCut]->Branch("ClusTrackDeltaPhi",&tESDIMClusTrackDeltaPhi,"tESDIMClusTrackDeltaPhi/F");
      tESDInvMassShowerShape[iCut]->Branch("ClusIsoSumClusterEt",&tESDIMClusterIsoSumClusterEt,"tESDIMClusterIsoSumClusterEt/F");
      tESDInvMassShowerShape[iCut]->Branch("ClusIsoSumTrackEt",&tESDIMClusterIsoSumTrackEt,"tESDIMClusterIsoSumTrackEt/F");
      fInvMassShowerShape[iCut]->Add(tESDInvMassShowerShape[iCut]);
    }

    if (fIsMC > 1 ){
      fHistoClusGammaPt[iCut]->Sumw2();
      fHistoClusGammaE[iCut]->Sumw2();
      fHistoClusOverlapHeadersGammaPt[iCut]->Sumw2();
      fHistoClusAllHeadersGammaPt[iCut]->Sumw2();
      fHistoClusRejectedHeadersGammaPt[iCut]->Sumw2();
      if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType()==2){
        if ( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsSpecialTrigger()==6 ){
          if(!fDoLightOutput){fHistoGoodMesonClusters[iCut]->Sumw2();}
        }
      }
      if(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetClusterType()==2){
        fHistoClusGammaE_BothBM_highestE[iCut]->Sumw2();
        //if ( (((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger()==1)||(((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger()==0) ){
          fHistoClusGammaE_BothBM[iCut]->Sumw2();
          fHistoClusGammaE_AnaBM_highestE[iCut]->Sumw2();
        //}
        if ( ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger()==6 ){
          fHistoClusGammaE_onlyTriggered[iCut]->Sumw2();
        }
      }
    }

    if(fDoLightOutput==2 && ((((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 1) || (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 3) || (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 4))) {
        nBinsMinv=200;
        if ( ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger() > 1 &&
            ( ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::k13TeV ||
                    ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::k13TeVLowB )){
            nBinsPt                   = 310;
        } else {
            nBinsPt                   = 238;
        }
        if(fPi0EtaSwitch == 1){ // pi0 only
            nBinsMinv=200*0.3/0.8;
            maxMinv = 0.3;
        } else if(fPi0EtaSwitch == 2){ // eta only
            nBinsMinv=200*0.3/0.8;
            minMinv = 0.3;
        }
    } else {
        if(fPi0EtaSwitch == 1){ // pi0 only
            nBinsMinv=400*0.3/0.8;
            maxMinv = 0.3;
        } else if(fPi0EtaSwitch == 2){ // eta only
            nBinsMinv=400*0.3/0.8;
            minMinv = 0.3;
        }
    }


    if(fDoMesonAnalysis){
      fHistoMotherInvMassPt[iCut]             = new TH2F("ESD_Mother_InvMass_Pt", "ESD_Mother_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
      fHistoMotherInvMassPt[iCut]->SetXTitle("M_{inv} (GeV/c^{2})");
      fHistoMotherInvMassPt[iCut]->SetYTitle("p_{T,pair} (GeV/c)");
      fESDList[iCut]->Add(fHistoMotherInvMassPt[iCut]);

      if(!fDoLightOutput){
        fHistoMotherMatchedInvMassPt[iCut]      = new TH2F("ESD_MotherMatched_InvMass_Pt", "ESD_MotherMatched_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
        fHistoMotherMatchedInvMassPt[iCut]->SetXTitle("M_{inv} (GeV/c^{2}) matched conv e^{+/-}to cluster");
        fHistoMotherMatchedInvMassPt[iCut]->SetYTitle("p_{T,pair} (GeV/c)");
        fESDList[iCut]->Add(fHistoMotherMatchedInvMassPt[iCut]);
      }

      fHistoMotherBackInvMassPt[iCut]         = new TH2F("ESD_Background_InvMass_Pt", "ESD_Background_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
      fHistoMotherBackInvMassPt[iCut]->SetXTitle("M_{inv, mxed}(GeV/c^{2})");
      fHistoMotherBackInvMassPt[iCut]->SetYTitle("p_{T,BG pair} (GeV/c)");
      fESDList[iCut]->Add(fHistoMotherBackInvMassPt[iCut]);

      if(!fDoLightOutput){
        fHistoMotherInvMassPtAlpha[iCut]        = new TH2F("ESD_Mother_InvMass_vs_Pt_Alpha", "ESD_Mother_InvMass_vs_Pt_Alpha", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
        fHistoMotherInvMassPtAlpha[iCut]->SetXTitle("M_{inv} (GeV/c^{2})");
        fHistoMotherInvMassPtAlpha[iCut]->SetYTitle("p_{T,pair} (GeV/c)");
        fESDList[iCut]->Add(fHistoMotherInvMassPtAlpha[iCut]);

        fHistoPhotonPairPtconv[iCut]            = new TH2F("ESD_Mother_InvMass_PtConv", "", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
        fHistoPhotonPairPtconv[iCut]->SetXTitle("M_{inv} (GeV/c^{2})");
        fHistoPhotonPairPtconv[iCut]->SetYTitle("#gamma^{conv} p_{T} (GeV/c)");
        fESDList[iCut]->Add(fHistoPhotonPairPtconv[iCut]);

        fHistoPhotonPairMixedEventPtconv[iCut]  = new TH2F("ESD_Background_InvMass_PtConv", "", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
        fHistoPhotonPairMixedEventPtconv[iCut]->SetXTitle("M_{inv,mixed}(GeV/c^{2})");
        fHistoPhotonPairMixedEventPtconv[iCut]->SetYTitle("#gamma^{conv} p_{T} (GeV/c)");
        fESDList[iCut]->Add(fHistoPhotonPairMixedEventPtconv[iCut]);
      }

      if (fIsMC > 1){
        fHistoMotherInvMassPt[iCut]->Sumw2();
        fHistoMotherBackInvMassPt[iCut]->Sumw2();
        if(!fDoLightOutput){
          fHistoMotherMatchedInvMassPt[iCut]->Sumw2();
          fHistoMotherInvMassPtAlpha[iCut]->Sumw2();
          fHistoPhotonPairPtconv[iCut]->Sumw2();
          fHistoPhotonPairMixedEventPtconv[iCut]->Sumw2();
        }
      }
      if(fDoMesonQA == 2){
        tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut] = new TTree(Form("%s_%s_%s_%s Meson DCA tree",cutstringEvent.Data() ,cutstringPhoton.Data(),cutstringCalo.Data(),cutstringMeson.Data()), "ESD_Mesons_InvMass_Pt_DcazMin_DcazMax_Flag");
        tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut]->Branch("InvMass",&fInvMass,"fInvMass/F");
        tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut]->Branch("Pt",&fPt,"fPt/F");
        tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut]->Branch("DcaZMin",&fDCAzGammaMin,"fDCAzGammaMin/F");
        tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut]->Branch("DcaZMax",&fDCAzGammaMax,"fDCAzGammaMax/F");
        tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut]->Branch("kind",&iFlag,"iFlag/b");
        if(fIsMC>0){
          tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut]->Branch("mesonMCInfo",&iMesonMCInfo,"iMesonMCInfo/b");
        }
      }
      if(!fDoLightOutput || fDoECalibOutput){
        fHistoMotherInvMassECalib[iCut]         = new TH2F("ESD_Mother_InvMass_E_Calib", "ESD_Mother_InvMass_E_Calib", 300, 0, 0.3, nBinsPt, arrPtBinning);
        fHistoMotherInvMassECalib[iCut]->SetXTitle("M_{inv} (GeV/c^{2})");
        fHistoMotherInvMassECalib[iCut]->SetYTitle("E_{cluster}(GeV)");
        fESDList[iCut]->Add(fHistoMotherInvMassECalib[iCut]);

        fHistoMotherBackInvMassECalib[iCut]     = new TH2F("ESD_Background_InvMass_E_Calib", "ESD_Background_InvMass_E_Calib", 300, 0, 0.3, nBinsPt, arrPtBinning);
        fHistoMotherBackInvMassECalib[iCut]->SetXTitle("M_{inv} (GeV/c^{2})");
        fHistoMotherBackInvMassECalib[iCut]->SetYTitle("E_{cluster}(GeV)");
        fESDList[iCut]->Add(fHistoMotherBackInvMassECalib[iCut]);

        fHistoMotherInvMassECalibPCM[iCut]         = new TH2F("ESD_Mother_InvMass_E_Calib_PCM", "ESD_Mother_InvMass_E_Calib_PCM", 300, 0, 0.3, nBinsPt, arrPtBinning);
        fHistoMotherInvMassECalibPCM[iCut]->SetXTitle("M_{inv} (GeV/c^{2})");
        fHistoMotherInvMassECalibPCM[iCut]->SetYTitle("E_{PCM #gamma}(GeV)");
        fESDList[iCut]->Add(fHistoMotherInvMassECalibPCM[iCut]);

        fHistoMotherBackInvMassECalibPCM[iCut]     = new TH2F("ESD_Background_InvMass_E_Calib_PCM", "ESD_Background_InvMass_E_Calib_PCM", 300, 0, 0.3, nBinsPt, arrPtBinning);
        fHistoMotherBackInvMassECalibPCM[iCut]->SetXTitle("M_{inv} (GeV/c^{2})");
        fHistoMotherBackInvMassECalibPCM[iCut]->SetYTitle("E_{PCM #gamma}(GeV)");
        fESDList[iCut]->Add(fHistoMotherBackInvMassECalibPCM[iCut]);

        if (fIsMC > 1){
          fHistoMotherInvMassECalib[iCut]->Sumw2();
          fHistoMotherBackInvMassECalib[iCut]->Sumw2();
          fHistoMotherInvMassECalibPCM[iCut]->Sumw2();
          fHistoMotherBackInvMassECalibPCM[iCut]->Sumw2();
        }
      }

      if (fDoMesonQA > 0 ){
        fHistoMotherPi0PtY[iCut]              = new TH2F("ESD_MotherPi0_Pt_Y", "ESD_MotherPi0_Pt_Y", nBinsQAPt, arrQAPtBinning, 150, -1.5, 1.5);
        fHistoMotherPi0PtY[iCut]->SetXTitle("p_{T, #pi^{0} cand} (GeV/c)");
        fHistoMotherPi0PtY[iCut]->SetYTitle("y_{#pi^{0} cand}");
        fESDList[iCut]->Add(fHistoMotherPi0PtY[iCut]);
        fHistoMotherEtaPtY[iCut]              = new TH2F("ESD_MotherEta_Pt_Y", "ESD_MotherEta_Pt_Y", nBinsQAPt, arrQAPtBinning, 150, -1.5, 1.5);
        fHistoMotherEtaPtY[iCut]->SetXTitle("p_{T, #eta cand} (GeV/c)");
        fHistoMotherEtaPtY[iCut]->SetYTitle("y_{ #eta cand}");
        fESDList[iCut]->Add(fHistoMotherEtaPtY[iCut]);
        fHistoMotherPi0PtAlpha[iCut]          = new TH2F("ESD_MotherPi0_Pt_Alpha", "ESD_MotherPi0_Pt_Alpha", nBinsQAPt, arrQAPtBinning, 200, -1, 1);
        fHistoMotherPi0PtAlpha[iCut]->SetXTitle("p_{T, #pi^{0} cand} (GeV/c)");
        fHistoMotherPi0PtAlpha[iCut]->SetYTitle("#alpha_{#pi^{0} cand}");
        fESDList[iCut]->Add(fHistoMotherPi0PtAlpha[iCut]);
        fHistoMotherEtaPtAlpha[iCut]          = new TH2F("ESD_MotherEta_Pt_Alpha", "ESD_MotherEta_Pt_Alpha", nBinsQAPt, arrQAPtBinning, 200, -1, 1);
        fHistoMotherEtaPtAlpha[iCut]->SetXTitle("p_{T, #eta cand} (GeV/c)");
        fHistoMotherEtaPtAlpha[iCut]->SetYTitle("#alpha_{#eta cand}");
        fESDList[iCut]->Add(fHistoMotherEtaPtAlpha[iCut]);
        fHistoMotherPi0PtOpenAngle[iCut]      = new TH2F("ESD_MotherPi0_Pt_OpenAngle", "ESD_MotherPi0_Pt_OpenAngle", nBinsQAPt, arrQAPtBinning, 100, 0, 1);
        fHistoMotherPi0PtOpenAngle[iCut]->SetXTitle("p_{T, #pi^{0} cand} (GeV/c)");
        fHistoMotherPi0PtOpenAngle[iCut]->SetYTitle("#theta_{#pi^{0} cand}");
        fESDList[iCut]->Add(fHistoMotherPi0PtOpenAngle[iCut]);
        fHistoMotherEtaPtOpenAngle[iCut]      = new TH2F("ESD_MotherEta_Pt_OpenAngle", "ESD_MotherEta_Pt_OpenAngle", nBinsQAPt, arrQAPtBinning, 200, 0, TMath::Pi());
        fHistoMotherEtaPtOpenAngle[iCut]->SetXTitle("p_{T, #eta cand} (GeV/c)");
        fHistoMotherEtaPtOpenAngle[iCut]->SetYTitle("#theta_{#eta cand}");
        fESDList[iCut]->Add(fHistoMotherEtaPtOpenAngle[iCut]);
        fHistoMotherPi0ConvPhotonEtaPhi[iCut] = new TH2F("ESD_MotherPi0ConvPhoton_Eta_Phi", "ConvPhoton under #pi^{0}peak", 600, 0, 2*TMath::Pi(), 200, -1, 1);
        fHistoMotherPi0ConvPhotonEtaPhi[iCut]->SetXTitle("#phi_{#gamma_{conv}}(rad)");
        fHistoMotherPi0ConvPhotonEtaPhi[iCut]->SetYTitle("#eta_{#gamma_{conv}}");
        fESDList[iCut]->Add(fHistoMotherPi0ConvPhotonEtaPhi[iCut]);
        fHistoMotherEtaConvPhotonEtaPhi[iCut] = new TH2F("ESD_MotherEtaConvPhoton_Eta_Phi", "ConvPhoton under #eta peak", 600, 0, 2*TMath::Pi(), 200, -1, 1);
        fHistoMotherEtaConvPhotonEtaPhi[iCut]->SetXTitle("#phi_{#gamma_{conv}}(rad)");
        fHistoMotherEtaConvPhotonEtaPhi[iCut]->SetYTitle("#eta_{#gamma_{conv}}");
        fESDList[iCut]->Add(fHistoMotherEtaConvPhotonEtaPhi[iCut]);
        if (fIsMC > 1){
          fHistoMotherPi0PtY[iCut]->Sumw2();
          fHistoMotherEtaPtY[iCut]->Sumw2();
          fHistoMotherPi0PtAlpha[iCut]->Sumw2();
          fHistoMotherEtaPtAlpha[iCut]->Sumw2();
          fHistoMotherPi0PtOpenAngle[iCut]->Sumw2();
          fHistoMotherEtaPtOpenAngle[iCut]->Sumw2();
          fHistoMotherPi0ConvPhotonEtaPhi[iCut]->Sumw2();
          fHistoMotherEtaConvPhotonEtaPhi[iCut]->Sumw2();
        }
      }
    }
    if(fDoJetAnalysis){

      fJetHistograms[iCut] = new TList();
      fJetHistograms[iCut]->SetOwner(kTRUE);
      fJetHistograms[iCut]->SetName(Form("%s_%s_%s_%s Jet histograms", cutstringEvent.Data(),cutstringPhoton.Data(),cutstringCalo.Data(),cutstringMeson.Data()));

      fHistoPtJet[iCut] = new TH1F("JetPt", "JetPt", 150, 0, 150);
      fJetHistograms[iCut]->Add(fHistoPtJet[iCut]);
      fHistoJetEta[iCut] = new TH1F("JetEta", "JetEta", 100, -1, 1);
      fJetHistograms[iCut]->Add(fHistoJetEta[iCut]);
      fHistoJetPhi[iCut] = new TH1F("JetPhi", "JetPhi", 70, 0, 7);
      fJetHistograms[iCut]->Add(fHistoJetPhi[iCut]);
      fHistoJetArea[iCut] = new TH1F("JetArea", "JetArea", 50, 0, 1);
      fJetHistograms[iCut]->Add(fHistoJetArea[iCut]);
      fHistoNJets[iCut] = new TH1F("NJets", "NJets", 10, 0, 10);
      fJetHistograms[iCut]->Add(fHistoNJets[iCut]);
      fHistoEventwJets[iCut] = new TH1F("NEvents_with_Jets", "NEvents_with_Jets", 5, 0, 5);
      fJetHistograms[iCut]->Add(fHistoEventwJets[iCut]);
      if(!fDoLightOutput){
        fHistoJetPi0PtRatio[iCut] = new TH1F("Ratio_Pt_Pi0_Jet", "Ratio_Pt_Pi0_Jet", 20, 0, 1.5);
        fJetHistograms[iCut]->Add(fHistoJetPi0PtRatio[iCut]);
        fHistoJetMotherInvMassPt[iCut] = new TH2F("ESD_Pi0Jet_Mother_InvMass_Pt", "ESD_Pi0Jet_Mother_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fJetHistograms[iCut]->Add(fHistoJetMotherInvMassPt[iCut]);
        fHistoEtaPhiJetPi0Cand[iCut] = new TH2F("Eta_Phi_Distr_Pi0Jet", "Eta_Phi_Distr_Pi0Jet", 20, 0, M_PI, 20, -1, 1);
        fJetHistograms[iCut]->Add(fHistoEtaPhiJetPi0Cand[iCut]);
        fHistoRJetPi0Cand[iCut] = new TH2F("ESD_RPi0Jet_Pt", "ESD_RPi0Jet_Pt", 35, 0, 3.5, nBinsPt, arrPtBinning);
        fJetHistograms[iCut]->Add(fHistoRJetPi0Cand[iCut]);
        fHistoPi0InJetMotherInvMassPt[iCut] = new TH2F("ESD_Pi0inJet_Mother_InvMass_Pt", "ESD_Pi0inJet_Mother_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fJetHistograms[iCut]->Add(fHistoPi0InJetMotherInvMassPt[iCut]);
        fHistoMotherBackJetInvMassPt[iCut] = new TH2F("ESD_Jet_Background_InvMass_Pt", "ESD_Jet_Background_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fJetHistograms[iCut]->Add(fHistoMotherBackJetInvMassPt[iCut]);
        fHistoEtaPhiJetWithPi0Cand[iCut] = new TH2F("Eta_Phi_Distr_Pi0inJet", "Eta_Phi_Distr_Pi0inJet", 15, 0, 0.4, 15, -0.4, 0.4);
        fJetHistograms[iCut]->Add(fHistoEtaPhiJetWithPi0Cand[iCut]);
        fHistoDoubleCounting[iCut] = new TH1F("Double_Counting_Mesons_Jets", "Double_Counting_Mesons_Jets", 6, 0, 6);
        fJetHistograms[iCut]->Add(fHistoDoubleCounting[iCut]);
        fHistoJetFragmFunc[iCut] = new TH2F("ESD_Pi0inJetPt_FragmentationFunc", "ESD_Pi0inJetPt_FragmentationFunc", 50, arrLogBinning, 150, 0., 150.);
        fJetHistograms[iCut]->Add(fHistoJetFragmFunc[iCut]);
        fHistoJetFragmFuncZInvMass[iCut] = new TH2F("ESD_Pi0inJetPt_Fragm_Z_InvMass", "ESD_Pi0inJetPt_Fragm_Z_InvMass", nBinsMinv, 0, maxMinv, 50, arrLogBinning);
        fJetHistograms[iCut]->Add(fHistoJetFragmFuncZInvMass[iCut]);
      }
    }
    if(fDoHBTHistoOutput){
      fHistoHBTOpeningAnglePt[iCut]                = new TH2F("ESD_OpeningAngle_Pt", "ESD_OpeningAngle_Pt", 700, 0, 3.5, nBinsPt, arrPtBinning);
      fHistoHBTOpeningAnglePt[iCut]->SetXTitle("#theta");
      fHistoHBTOpeningAnglePt[iCut]->SetYTitle("p_{T}(GeV/c)");
      if (fIsMC > 1) fHistoHBTOpeningAnglePt[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoHBTOpeningAnglePt[iCut]);

      fHistoBckHBTOpeningAnglePt[iCut]                = new TH2F("ESD_Bck_OpeningAngle_Pt", "ESD_Bck_OpeningAngle_Pt", 700, 0, 3.5, nBinsPt, arrPtBinning);
      fHistoBckHBTOpeningAnglePt[iCut]->SetXTitle("#theta");
      fHistoBckHBTOpeningAnglePt[iCut]->SetYTitle("p_{T}(GeV/c)");
      if (fIsMC > 1) fHistoBckHBTOpeningAnglePt[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoBckHBTOpeningAnglePt[iCut]);

      fHistoHBTDeltaEPt[iCut]                = new TH2F("ESD_DeltaE_Pt", "ESD_DeltaE_Pt", 800, 0, 4.0, nBinsPt, arrPtBinning);
      fHistoHBTDeltaEPt[iCut]->SetXTitle("#Delta E (GeV)");
      fHistoHBTDeltaEPt[iCut]->SetYTitle("p_{T}(GeV/c)");
      if (fIsMC > 1) fHistoHBTDeltaEPt[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoHBTDeltaEPt[iCut]);

      fHistoBckHBTDeltaEPt[iCut]                = new TH2F("ESD_Bck_DeltaE_Pt", "ESD_Bck_DeltaE_Pt", 800, 0, 4.0, nBinsPt, arrPtBinning);
      fHistoBckHBTDeltaEPt[iCut]->SetXTitle("#Delta E (GeV)");
      fHistoBckHBTDeltaEPt[iCut]->SetYTitle("p_{T}(GeV/c)");
      if (fIsMC > 1) fHistoBckHBTDeltaEPt[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoBckHBTDeltaEPt[iCut]);
    }

    if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 2){
      fCaloTriggerMimicHelper[iCut] = NULL;
      if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetTriggerMimicking() > 0){
        fCaloTriggerMimicHelper[iCut] = (AliCaloTriggerMimicHelper*) (AliAnalysisManager::GetAnalysisManager()->GetTask(fCaloTriggerHelperName.Data()));
        if(fCaloTriggerMimicHelper[iCut]){
          if ( fSetEventCutsOutputlist[cutstringEvent] == kFALSE ) {
            fSetEventCutsOutputlist[cutstringEvent]=kTRUE;
            fOutputContainer->Add(fCaloTriggerMimicHelper[iCut]->GetTriggerMimicHelperHistograms());
          }
        }
      }
    }
  }
  if(fDoMesonAnalysis){
    InitBack(); // Init Background Handler
  }

  if(fIsMC>0){
    // MC Histogramms
    fMCList                                         = new TList*[fnCuts];
    // True Histogramms
    fTrueList                                       = new TList*[fnCuts];

    if(!fDoLightOutput){
      fHistoMCHeaders                                 = new TH1I*[fnCuts];
      fHistoMCAllGammaPt                              = new TH1F*[fnCuts];
      fHistoMCAllGammaEMCALAccPt                      = new TH1F*[fnCuts];
      fHistoMCAllSecondaryGammaPt                     = new TH2F*[fnCuts];
      fHistoMCDecayGammaPi0Pt                         = new TH1F*[fnCuts];
      fHistoMCDecayGammaRhoPt                         = new TH1F*[fnCuts];
      fHistoMCDecayGammaEtaPt                         = new TH1F*[fnCuts];
      fHistoMCDecayGammaOmegaPt                       = new TH1F*[fnCuts];
      fHistoMCDecayGammaEtapPt                        = new TH1F*[fnCuts];
      fHistoMCDecayGammaPhiPt                         = new TH1F*[fnCuts];
      fHistoMCDecayGammaSigmaPt                       = new TH1F*[fnCuts];
      fHistoMCConvGammaPt                             = new TH1F*[fnCuts];
      fHistoMCSecondaryConvGammaPt                    = new TH2F*[fnCuts];
      fHistoTrueConvGammaPt                           = new TH1F*[fnCuts];
      fHistoDoubleCountTrueConvGammaRPt               = new TH2F*[fnCuts];
      fHistoMultipleCountTrueConvGamma                = new TH1F*[fnCuts];

      fHistoCombinatorialPt                           = new TH2F*[fnCuts];
      fHistoTruePrimaryConvGammaPt                    = new TH1F*[fnCuts];
      fHistoTruePrimaryConvGammaESDPtMCPt             = new TH2F*[fnCuts];
      fHistoTrueSecondaryConvGammaFromXFromK0sMCPtESDPt     = new TH2F*[fnCuts];
      fHistoTrueSecondaryConvGammaFromXFromK0lMCPtESDPt     = new TH2F*[fnCuts];
      fHistoTrueSecondaryConvGammaFromXFromLambdaMCPtESDPt  = new TH2F*[fnCuts];
      fHistoTrueSecondaryConvGammaPt                        = new TH2F*[fnCuts];
      fHistoTrueSecondaryConvGammaMCPt                      = new TH2F*[fnCuts];
    }

    fHistoTrueClusGammaPt                           = new TH1F*[fnCuts];
    if(!fDoLightOutput){
      fHistoTrueClusConvGammaPt                       = new TH1F*[fnCuts];
      fHistoTruePrimaryClusGammaPt                    = new TH1F*[fnCuts];
      fHistoTruePrimaryClusGammaESDPtMCPt             = new TH2F*[fnCuts];
      fHistoTruePrimaryClusConvGammaPt                = new TH1F*[fnCuts];
      fHistoTruePrimaryClusConvGammaESDPtMCPt         = new TH2F*[fnCuts];
      fHistoTrueSecondaryClusGammaPt                  = new TH1F*[fnCuts];
      fHistoTrueSecondaryClusGammaFromK0sPt           = new TH1F*[fnCuts];
      fHistoTrueSecondaryClusGammaFromK0lPt           = new TH1F*[fnCuts];
      fHistoTrueSecondaryClusGammaFromLambdaPt        = new TH1F*[fnCuts];
      fHistoTrueClusElectronPt                        = new TH1F*[fnCuts];
    }
    fHistoDoubleCountTrueClusterGammaPt             = new TH2F*[fnCuts];
    fHistoMultipleCountTrueClusterGamma             = new TH1F*[fnCuts];

    if(!fDoLightOutput){
      fHistoTrueNLabelsInClusPt                       = new TH2F*[fnCuts];
      fHistoTrueClusGammaEM02                         = new TH2F*[fnCuts];
      fHistoTrueClusPi0EM02                           = new TH2F*[fnCuts];
      fHistoTrueClusEMNonLeadingPt                    = new TH1F*[fnCuts];
    }

    if(!fDoLightOutput){
      fHistoTruePi0InvMassECalib                      = new TH2F*[fnCuts];
      fHistoTruePi0PureGammaInvMassECalib             = new TH2F*[fnCuts];
      fHistoTruePi0InvMassECalibPCM                   = new TH2F*[fnCuts];
    }

    if (fDoClusterQA > 0 && fDoLightOutput!=2){
      fHistoTrueClusConvGammaFullyPt                = new TH1F*[fnCuts];
      fHistoTrueClusMergedGammaPt                   = new TH1F*[fnCuts];
      fHistoTrueClusMergedPartConvGammaPt           = new TH1F*[fnCuts];
      fHistoTrueClusDalitzPt                        = new TH1F*[fnCuts];
      fHistoTrueClusDalitzMergedPt                  = new TH1F*[fnCuts];
      fHistoTrueClusPhotonFromElecMotherPt          = new TH1F*[fnCuts];
      fHistoTrueClusShowerPt                        = new TH1F*[fnCuts];
      fHistoTrueClusSubLeadingPt                    = new TH1F*[fnCuts];
      fHistoTrueClusNMothers                        = new TH1F*[fnCuts];
    }

    if(fDoMesonAnalysis){
      fHistoMCPi0Pt                                 = new TH1F*[fnCuts];
      fHistoMCPi0WOWeightPt                         = new TH1F*[fnCuts];
      fHistoMCEtaPt                                 = new TH1F*[fnCuts];
      fHistoMCEtaWOWeightPt                         = new TH1F*[fnCuts];
      fHistoMCPi0InAccPt                            = new TH1F*[fnCuts];
      fHistoMCPi0WOWeightInAccPt                    = new TH1F*[fnCuts];
      fHistoMCEtaInAccPt                            = new TH1F*[fnCuts];
      fHistoMCEtaWOWeightInAccPt                    = new TH1F*[fnCuts];
      if (fIsMC > 1){
        fHistoMCPi0WOEvtWeightPt                    = new TH1F*[fnCuts];
        fHistoMCEtaWOEvtWeightPt                    = new TH1F*[fnCuts];
        fHistoMCPi0WOEvtWeightInAccPt               = new TH1F*[fnCuts];
        fHistoMCEtaWOEvtWeightInAccPt               = new TH1F*[fnCuts];
      }
      fHistoMCSecPi0PtvsSource                      = new TH2F*[fnCuts];
      fHistoMCSecPi0InAccPtvsSource                 = new TH2F*[fnCuts];
      if (fDoLightOutput!=2){
        fHistoMCPrimaryPtvsSource                     = new TH2F*[fnCuts];
        fHistoMCSecPi0Source                          = new TH1F*[fnCuts];
        fHistoMCSecEtaPt                              = new TH1F*[fnCuts];
        fHistoMCSecEtaSource                          = new TH1F*[fnCuts];
      }
      if (!fDoLightOutput){
        fHistoMCPi0PtGammaLeg                       = new TH2F*[fnCuts];
        fHistoMCPi0WOWeightPtGammaLeg               = new TH2F*[fnCuts];
        fHistoMCPi0InAccPtGammaLeg                  = new TH2F*[fnCuts];
        fHistoMCPi0WOWeightInAccPtGammaLeg          = new TH2F*[fnCuts];
        fHistoMCSecPi0PtGamma1vsSource              = new TH2F*[fnCuts];
        fHistoMCSecPi0PtGamma2vsSource              = new TH2F*[fnCuts];
        fHistoMCSecPi0InAccPtGamma1vsSource         = new TH2F*[fnCuts];
        fHistoMCSecPi0InAccPtGamma2vsSource         = new TH2F*[fnCuts];
      }

      fHistoTruePi0InvMassPt                        = new TH2F*[fnCuts];
      fHistoTrueEtaInvMassPt                        = new TH2F*[fnCuts];
      if (fDoLightOutput!=2){
          fHistoTruePi0InvMassPtAdditional            = new TH2F*[fnCuts];
          fHistoTrueEtaInvMassPtAdditional            = new TH2F*[fnCuts];
        fHistoTruePi0MatchedInvMassPt                 = new TH2F*[fnCuts];
        fHistoTrueEtaMatchedInvMassPt                 = new TH2F*[fnCuts];
        fHistoDoubleCountTruePi0InvMassPt             = new TH2F*[fnCuts];
        fHistoMultipleCountTruePi0                    = new TH1F*[fnCuts];
        fHistoDoubleCountTrueEtaInvMassPt             = new TH2F*[fnCuts];
        fHistoMultipleCountTrueEta                    = new TH1F*[fnCuts];
      }
      fHistoTruePrimaryPi0InvMassPt                 = new TH2F*[fnCuts];
      fHistoTruePrimaryEtaInvMassPt                 = new TH2F*[fnCuts];
      fHistoTruePrimaryPi0W0WeightingInvMassPt      = new TH2F*[fnCuts];
      fHistoTruePrimaryEtaW0WeightingInvMassPt      = new TH2F*[fnCuts];
      fProfileTruePrimaryPi0WeightsInvMassPt        = new TProfile2D*[fnCuts];
      fProfileTruePrimaryEtaWeightsInvMassPt        = new TProfile2D*[fnCuts];
      fHistoTrueSecondaryPi0InvMassPt               = new TH2F*[fnCuts];
      fHistoTrueSecondaryPi0FromK0sInvMassPt        = new TH2F*[fnCuts];
      fHistoTrueSecondaryPi0FromK0lInvMassPt        = new TH2F*[fnCuts];
      fHistoTrueSecondaryPi0FromEtaInvMassPt        = new TH2F*[fnCuts];
      fHistoTrueSecondaryPi0FromLambdaInvMassPt     = new TH2F*[fnCuts];
      if(!fDoLightOutput) {
        fHistoTruePrimaryPi0PhotonPairPtconv          = new TH2F*[fnCuts];
        fHistoTruePrimaryPi0W0WeightsPhotonPairPtconv = new TH2F*[fnCuts];
        fHistoTrueSecondaryPi0PhotonPairPtconv        = new TH2F*[fnCuts];
        fHistoTrueSecondaryPi0FromK0sPhotonPairPtconv = new TH2F*[fnCuts];
        fHistoTrueSecondaryPi0FromK0lPhotonPairPtconv = new TH2F*[fnCuts];
        fHistoTrueSecondaryPi0FromLambdaPhotonPairPtconv= new TH2F*[fnCuts];
        fHistoTruePrimaryPi0DCPtconv                  = new TH1F*[fnCuts];
        fHistoTrueSecondaryPi0DCPtconvSource          = new TH2F*[fnCuts];
        fHistoTruePrimaryPi0MissingPtconv             = new TH1F*[fnCuts];
        fHistoTrueSecondaryPi0MissingPtconvSource     = new TH2F*[fnCuts];
      }
      if(fDoHBTHistoOutput){
        fHistoTrueGammaInvMassPt                      = new TH2F*[fnCuts];
        fHistoMCGammaHBTInvMassPt                     = new TH2F*[fnCuts];
        fHistoBckTrueGammaInvMassPt                   = new TH2F*[fnCuts];
        fHistoBckMCGammaHBTInvMassPt                  = new TH2F*[fnCuts];
        fHistoMCGammaHBTOpeningAnglePt                = new TH2F*[fnCuts];
        fHistoBckMCGammaHBTOpeningAnglePt             = new TH2F*[fnCuts];
        fHistoMCGammaHBTDeltaEPt                      = new TH2F*[fnCuts];
        fHistoBckMCGammaHBTDeltaEPt                   = new TH2F*[fnCuts];
      }


      if (fDoMesonQA > 0){
        fHistoMCPi0PtY                              = new TH2F*[fnCuts];
        fHistoMCEtaPtY                              = new TH2F*[fnCuts];
        fHistoMCPi0PtAlpha                          = new TH2F*[fnCuts];
        fHistoMCEtaPtAlpha                          = new TH2F*[fnCuts];
        if (fIsMC == 2){
          fHistoMCPi0PtJetPt                        = new TH2F*[fnCuts];
          fHistoMCEtaPtJetPt                        = new TH2F*[fnCuts];
        }

        if (fIsMC < 2){
          fHistoTruePi0CaloPhotonInvMassPt                  = new TH2F*[fnCuts];
          fHistoTrueEtaCaloPhotonInvMassPt                  = new TH2F*[fnCuts];
          fHistoTruePi0CaloConvertedPhotonInvMassPt         = new TH2F*[fnCuts];
          fHistoTruePi0CaloConvertedPhotonMatchedInvMassPt  = new TH2F*[fnCuts];
          fHistoTruePi0CaloConvPhotonConvRPt                = new TH2F*[fnCuts];
          fHistoTruePi0CaloConvPhotonConvRAlpha             = new TH2F*[fnCuts];
          fHistoTruePi0CaloConvPhotonPtAlpha                = new TH2F*[fnCuts];
          fHistoTrueEtaCaloConvertedPhotonInvMassPt         = new TH2F*[fnCuts];
          fHistoTrueEtaCaloConvertedPhotonMatchedInvMassPt  = new TH2F*[fnCuts];
          fHistoTrueEtaCaloConvPhotonConvRPt                = new TH2F*[fnCuts];
          fHistoTrueEtaCaloConvPhotonConvRAlpha             = new TH2F*[fnCuts];
          fHistoTrueEtaCaloConvPhotonPtAlpha                = new TH2F*[fnCuts];
          fHistoTruePi0CaloElectronInvMassPt                = new TH2F*[fnCuts];
          fHistoTrueEtaCaloElectronInvMassPt                = new TH2F*[fnCuts];
          fHistoTruePi0CaloMergedClusterInvMassPt           = new TH2F*[fnCuts];
          fHistoTrueEtaCaloMergedClusterInvMassPt           = new TH2F*[fnCuts];
          fHistoTruePi0CaloMergedClusterPartConvInvMassPt   = new TH2F*[fnCuts];
          fHistoTrueEtaCaloMergedClusterPartConvInvMassPt   = new TH2F*[fnCuts];
          fHistoTrueMotherCaloEMNonLeadingInvMassPt         = new TH2F*[fnCuts];
          fHistoTruePrimaryPi0MCPtResolPt                   = new TH2F*[fnCuts];
          fHistoTruePrimaryEtaMCPtResolPt                   = new TH2F*[fnCuts];
          fHistoTrueK0sWithPi0DaughterMCPt                  = new TH1F*[fnCuts];
          fHistoTrueK0lWithPi0DaughterMCPt                  = new TH1F*[fnCuts];
          fHistoTrueEtaWithPi0DaughterMCPt                  = new TH1F*[fnCuts];
          fHistoTrueLambdaWithPi0DaughterMCPt               = new TH1F*[fnCuts];
        }
        if(fDoMesonQA > 1){
          fHistoTrueBckGGInvMassPt                          = new TH2F*[fnCuts];
          fHistoTrueBckFullMesonContainedInOneClusterInvMassPt = new TH2F*[fnCuts];
          fHistoTrueBckAsymEClustersInvMassPt                  = new TH2F*[fnCuts];
          fHistoTrueBckContInvMassPt                        = new TH2F*[fnCuts];
        }
        fHistoTruePi0PtY                            = new TH2F*[fnCuts];
        fHistoTrueEtaPtY                            = new TH2F*[fnCuts];
        fHistoTruePi0PtAlpha                        = new TH2F*[fnCuts];
        fHistoTrueEtaPtAlpha                        = new TH2F*[fnCuts];
        fHistoTruePi0PtOpenAngle                    = new TH2F*[fnCuts];
        fHistoTrueEtaPtOpenAngle                    = new TH2F*[fnCuts];
        fHistoTrueMotherPi0ConvPhotonEtaPhi         = new TH2F*[fnCuts];
        fHistoTrueMotherEtaConvPhotonEtaPhi         = new TH2F*[fnCuts];
      }
    }
    if(fDoJetAnalysis && !fDoLightOutput) {
      fTrueJetHistograms                      = new TList*[fnCuts];
      fHistoTruevsRecJetPt                    = new TH2F*[fnCuts];
      fHistoTruePi0JetMotherInvMassPt         = new TH2F*[fnCuts];
      fHistoTruePi0InJetMotherInvMassPt       = new TH2F*[fnCuts];
      fHistoTruePrimaryPi0JetInvMassPt        = new TH2F*[fnCuts];
      fHistoTruePrimaryPi0inJetInvMassPt      = new TH2F*[fnCuts];
      fHistoTruePrimaryPi0InJetInvMassTruePt  = new TH2F*[fnCuts];
      fHistoTrueDoubleCountingPi0Jet          = new TH1F*[fnCuts];
      fHistoTruePi0JetFragmFunc               = new TH2F*[fnCuts];
      fHistoTruePi0JetFragmFuncZInvMass       = new TH2F*[fnCuts];
      fHistoMCPi0JetInAccPt                   = new TH1F*[fnCuts];
      fHistoMCPi0inJetInAccPt                 = new TH1F*[fnCuts];
      fHistoMCPi0JetEventGenerated            = new TH1F*[fnCuts];
      fHistoMCPi0inJetGenerated               = new TH1F*[fnCuts];
      fHistoTrueSecondaryPi0FromK0sJetInvMassPt      = new TH2F*[fnCuts];
      fHistoTrueSecondaryPi0FromK0sinJetInvMassPt    = new TH2F*[fnCuts];
      fHistoTrueSecondaryPi0FromLambdaJetInvMassPt   = new TH2F*[fnCuts];
      fHistoTrueSecondaryPi0FromLambdainJetInvMassPt = new TH2F*[fnCuts];
      fHistoTrueSecondaryPi0FromK0lJetInvMassPt      = new TH2F*[fnCuts];
      fHistoTrueSecondaryPi0FromK0linJetInvMassPt    = new TH2F*[fnCuts];
      fHistoTrueSecondaryPi0InvJetMassPt             = new TH2F*[fnCuts];
      fHistoTrueSecondaryPi0InvinJetMassPt           = new TH2F*[fnCuts];
      fHistoMotherPi0inJetPtY                        = new TH2F*[fnCuts];
      fHistoMotherPi0inJetPtPhi                      = new TH2F*[fnCuts];
      fHistoTrueDoubleCountingEtaJet          = new TH1F*[fnCuts];
      fHistoTrueEtaJetMotherInvMassPt         = new TH2F*[fnCuts];
      fHistoTrueEtaInJetMotherInvMassPt       = new TH2F*[fnCuts];
      fHistoTruePrimaryEtaJetInvMassPt        = new TH2F*[fnCuts];
      fHistoTruePrimaryEtainJetInvMassPt      = new TH2F*[fnCuts];
      fHistoTrueEtaJetFragmFunc               = new TH2F*[fnCuts];
      fHistoTrueEtaJetFragmFuncZInvMass       = new TH2F*[fnCuts];
      fHistoMCEtaJetInAccPt                   = new TH1F*[fnCuts];
      fHistoMCEtainJetInAccPt                 = new TH1F*[fnCuts];
      fHistoMCEtaJetEventGenerated            = new TH1F*[fnCuts];
      fHistoMCEtainJetGenerated               = new TH1F*[fnCuts];
      fHistoMotherEtainJetPtY                 = new TH2F*[fnCuts];
      fHistoMotherEtainJetPtPhi               = new TH2F*[fnCuts];
    }
    if(fDoJetQA){
        if(fDoLightOutput){
          fTrueJetHistograms                           = new TList*[fnCuts];
        }
        fHistoUnfoldingAsData                          = new TH2F*[fnCuts];
        fHistoUnfoldingMissed                          = new TH2F*[fnCuts];
        fHistoUnfoldingReject                          = new TH2F*[fnCuts];
        fHistoUnfoldingAsDataInvMassZ                  = new TH2F*[fnCuts];
        fHistoUnfoldingMissedInvMassZ                  = new TH2F*[fnCuts];
        fHistoUnfoldingRejectInvMassZ                  = new TH2F*[fnCuts];
      }

    for(Int_t iCut = 0; iCut<fnCuts;iCut++){
      TString cutstringEvent    = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
      TString cutstringPhoton   = ((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutNumber();
      TString cutstringCalo     = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
      TString cutstringMeson    = "NoMesonCut";
      if(fDoMesonAnalysis)
        cutstringMeson          = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();

      fMCList[iCut]                     = new TList();
      fMCList[iCut]->SetName(Form("%s_%s_%s_%s MC histograms",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
      fMCList[iCut]->SetOwner(kTRUE);
      fCutFolder[iCut]->Add(fMCList[iCut]);
      if(!fDoLightOutput){
        fHistoMCHeaders[iCut]             = new TH1I("MC_Headers", "MC_Headers", 20, 0, 20);
        fMCList[iCut]->Add(fHistoMCHeaders[iCut]);
        fHistoMCAllGammaPt[iCut]          = new TH1F("MC_AllGamma_Pt", "MC_AllGamma_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fMCList[iCut]->Add(fHistoMCAllGammaPt[iCut]);
        fHistoMCAllGammaEMCALAccPt[iCut]  = new TH1F("MC_AllGammaEMCALAcc_Pt", "MC_AllGammaEMCALAcc_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fMCList[iCut]->Add(fHistoMCAllGammaEMCALAccPt[iCut]);
        fHistoMCAllSecondaryGammaPt[iCut]    = new TH2F("MC_AllSecondaryGamma_Pt", "MC_AllSecondaryGamma_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt, 4, -0.5, 3.5);
        fHistoMCAllSecondaryGammaPt[iCut]->GetYaxis()->SetBinLabel(1,"K0s");
        fHistoMCAllSecondaryGammaPt[iCut]->GetYaxis()->SetBinLabel(2,"K0l");
        fHistoMCAllSecondaryGammaPt[iCut]->GetYaxis()->SetBinLabel(3,"Lambda");
        fHistoMCAllSecondaryGammaPt[iCut]->GetYaxis()->SetBinLabel(4,"rest");
        fMCList[iCut]->Add(fHistoMCAllSecondaryGammaPt[iCut]);

        fHistoMCDecayGammaPi0Pt[iCut]     = new TH1F("MC_DecayGammaPi0_Pt", "MC_DecayGammaPi0_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fMCList[iCut]->Add(fHistoMCDecayGammaPi0Pt[iCut]);
        fHistoMCDecayGammaRhoPt[iCut]     = new TH1F("MC_DecayGammaRho_Pt", "MC_DecayGammaRho_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fMCList[iCut]->Add(fHistoMCDecayGammaRhoPt[iCut]);
        fHistoMCDecayGammaEtaPt[iCut]     = new TH1F("MC_DecayGammaEta_Pt", "MC_DecayGammaEta_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fMCList[iCut]->Add(fHistoMCDecayGammaEtaPt[iCut]);
        fHistoMCDecayGammaOmegaPt[iCut]   = new TH1F("MC_DecayGammaOmega_Pt", "MC_DecayGammaOmmega_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fMCList[iCut]->Add(fHistoMCDecayGammaOmegaPt[iCut]);
        fHistoMCDecayGammaEtapPt[iCut]    = new TH1F("MC_DecayGammaEtap_Pt", "MC_DecayGammaEtap_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fMCList[iCut]->Add(fHistoMCDecayGammaEtapPt[iCut]);
        fHistoMCDecayGammaPhiPt[iCut]     = new TH1F("MC_DecayGammaPhi_Pt", "MC_DecayGammaPhi_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fMCList[iCut]->Add(fHistoMCDecayGammaPhiPt[iCut]);
        fHistoMCDecayGammaSigmaPt[iCut]   = new TH1F("MC_DecayGammaSigma_Pt", "MC_DecayGammaSigma_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fMCList[iCut]->Add(fHistoMCDecayGammaSigmaPt[iCut]);
        fHistoMCConvGammaPt[iCut]         = new TH1F("MC_ConvGamma_Pt", "MC_ConvGamma_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fMCList[iCut]->Add(fHistoMCConvGammaPt[iCut]);
        fHistoMCSecondaryConvGammaPt[iCut]  = new TH2F("MC_SecondaryConvGamma_Pt", "MC_SecondaryConvGamma_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt, 4, -0.5, 3.5);
        fHistoMCSecondaryConvGammaPt[iCut]->GetYaxis()->SetBinLabel(1,"K0s");
        fHistoMCSecondaryConvGammaPt[iCut]->GetYaxis()->SetBinLabel(2,"K0l");
        fHistoMCSecondaryConvGammaPt[iCut]->GetYaxis()->SetBinLabel(3,"Lambda");
        fHistoMCSecondaryConvGammaPt[iCut]->GetYaxis()->SetBinLabel(4,"rest");
        fMCList[iCut]->Add(fHistoMCSecondaryConvGammaPt[iCut]);

        if (fIsMC > 1){
          fHistoMCAllGammaPt[iCut]->Sumw2();
          fHistoMCAllGammaEMCALAccPt[iCut]->Sumw2();
          fHistoMCAllSecondaryGammaPt[iCut]->Sumw2();
          fHistoMCDecayGammaPi0Pt[iCut]->Sumw2();
          fHistoMCDecayGammaRhoPt[iCut]->Sumw2();
          fHistoMCDecayGammaEtaPt[iCut]->Sumw2();
          fHistoMCDecayGammaOmegaPt[iCut]->Sumw2();
          fHistoMCDecayGammaEtapPt[iCut]->Sumw2();
          fHistoMCDecayGammaPhiPt[iCut]->Sumw2();
          fHistoMCDecayGammaSigmaPt[iCut]->Sumw2();
          fHistoMCConvGammaPt[iCut]->Sumw2();
          fHistoMCSecondaryConvGammaPt[iCut]->Sumw2();
        }
      }

      if(fDoMesonAnalysis){
        fHistoMCPi0Pt[iCut]               = new TH1F("MC_Pi0_Pt", "MC_Pi0_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fHistoMCPi0Pt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0Pt[iCut]);
        fHistoMCPi0WOWeightPt[iCut]       = new TH1F("MC_Pi0_WOWeights_Pt", "MC_Pi0_WOWeights_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fMCList[iCut]->Add(fHistoMCPi0WOWeightPt[iCut]);

        fHistoMCEtaPt[iCut]               = new TH1F("MC_Eta_Pt", "MC_Eta_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fHistoMCEtaPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCEtaPt[iCut]);
        fHistoMCEtaWOWeightPt[iCut]       = new TH1F("MC_Eta_WOWeights_Pt", "MC_Eta_WOWeights_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fMCList[iCut]->Add(fHistoMCEtaWOWeightPt[iCut]);

        fHistoMCPi0InAccPt[iCut]          = new TH1F("MC_Pi0InAcc_Pt", "MC_Pi0InAcc_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fHistoMCPi0InAccPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0InAccPt[iCut]);
        fHistoMCPi0WOWeightInAccPt[iCut]  = new TH1F("MC_Pi0WOWeightInAcc_Pt", "MC_Pi0WOWeightInAcc_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fMCList[iCut]->Add(fHistoMCPi0WOWeightInAccPt[iCut]);
        fHistoMCEtaInAccPt[iCut]          = new TH1F("MC_EtaInAcc_Pt", "MC_EtaInAcc_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fHistoMCEtaInAccPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCEtaInAccPt[iCut]);
        fHistoMCEtaWOWeightInAccPt[iCut]  = new TH1F("MC_EtaWOWeightInAcc_Pt", "MC_EtaWOWeightInAcc_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fMCList[iCut]->Add(fHistoMCEtaWOWeightInAccPt[iCut]);

        if (fIsMC > 1){
          fHistoMCPi0WOWeightPt[iCut]->Sumw2();
          fHistoMCEtaWOWeightPt[iCut]->Sumw2();
          fHistoMCPi0WOWeightInAccPt[iCut]->Sumw2();
          fHistoMCEtaWOWeightInAccPt[iCut]->Sumw2();
          fHistoMCPi0WOEvtWeightPt[iCut]  = new TH1F("MC_Pi0_WOEventWeights_Pt", "MC_Pi0_WOEventWeights_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
          fMCList[iCut]->Add(fHistoMCPi0WOEvtWeightPt[iCut]);
          fHistoMCEtaWOEvtWeightPt[iCut]  = new TH1F("MC_Eta_WOEventWeights_Pt", "MC_Eta_WOEventWeights_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
          fMCList[iCut]->Add(fHistoMCEtaWOEvtWeightPt[iCut]);
          fHistoMCPi0WOEvtWeightInAccPt[iCut]  = new TH1F("MC_Pi0_WOEventWeightsInAcc_Pt", "MC_Pi0_WOEventWeightsInAcc_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
          fMCList[iCut]->Add(fHistoMCPi0WOEvtWeightInAccPt[iCut]);
          fHistoMCEtaWOEvtWeightInAccPt[iCut]  = new TH1F("MC_Eta_WOEventWeightsInAcc_Pt", "MC_Eta_WOEventWeightsInAcc_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
          fMCList[iCut]->Add(fHistoMCEtaWOEvtWeightInAccPt[iCut]);

          if (fDoMesonQA > 0 && fIsMC == 2){
            fHistoMCPi0PtJetPt[iCut]      = new TH2F("MC_Pi0_Pt_JetPt", "MC_Pi0_Pt_JetPt", nBinsQAPt, arrQAPtBinning, 200, 0, 200);
            fHistoMCPi0PtJetPt[iCut]->Sumw2();
            fMCList[iCut]->Add(fHistoMCPi0PtJetPt[iCut]);
            fHistoMCEtaPtJetPt[iCut]      = new TH2F("MC_Eta_Pt_JetPt", "MC_Eta_Pt_JetPt", nBinsQAPt, arrQAPtBinning, 200, 0, 200);
            fHistoMCEtaPtJetPt[iCut]->Sumw2();
            fMCList[iCut]->Add(fHistoMCEtaPtJetPt[iCut]);
          }
        }
        if(fDoLightOutput!=2){
            fHistoMCPrimaryPtvsSource[iCut]   = new TH2F("MC_Primary_Pt_Source", "MC_Primary_Pt_Source", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt, 7, -0.5, 6.5);
            fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(1,"Pi+");
            fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(2,"Pi-");
            fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(3,"K+");
            fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(4,"K-");
            fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(5,"K0s");
            fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(6,"K0l");
            fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(7,"Lambda");
            fMCList[iCut]->Add(fHistoMCPrimaryPtvsSource[iCut]);

            fHistoMCSecPi0Source[iCut]    = new TH1F("MC_SecPi0_Source", "MC_SecPi0_Source", 5000, 0., 5000);
            fMCList[iCut]->Add(fHistoMCSecPi0Source[iCut]);
            fHistoMCSecEtaSource[iCut]    = new TH1F("MC_SecEta_Source", "MC_SecEta_Source", 5000, 0, 5000);
            fMCList[iCut]->Add(fHistoMCSecEtaSource[iCut]);
            fHistoMCSecEtaPt[iCut]          = new TH1F("MC_SecEta_Pt", "MC_SecEta_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
            fMCList[iCut]->Add(fHistoMCSecEtaPt[iCut]);
        }
        fHistoMCSecPi0PtvsSource[iCut]  = new TH2F("MC_SecPi0_Pt_Source", "MC_SecPi0_Pt_Source", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt, 16, -0.5, 15.5);
        fMCList[iCut]->Add(fHistoMCSecPi0PtvsSource[iCut]);
        fHistoMCSecPi0InAccPtvsSource[iCut]  = new TH2F("MC_SecPi0InAcc_Pt_Source", "MC_SecPi0InAcc_Pt_Source", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt, 16, -0.5, 15.5);
        fMCList[iCut]->Add(fHistoMCSecPi0InAccPtvsSource[iCut]);

        if (fIsMC == 2){
          if(fDoLightOutput!=2){
              fHistoMCPrimaryPtvsSource[iCut]->Sumw2();
              fHistoMCSecEtaPt[iCut]->Sumw2();
          }
          fHistoMCSecPi0PtvsSource[iCut]->Sumw2();
          fHistoMCSecPi0InAccPtvsSource[iCut]->Sumw2();
        }

        // book histograms for pure MC handling of PCM-Calo dir gamma reco
        if (!fDoLightOutput){
          fHistoMCPi0PtGammaLeg[iCut]               = new TH2F("MC_Pi0_PtGamma_Leg", "MC_Pi0_PtGamma_Leg", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt, 3, -0.5, 2.5);
          fHistoMCPi0PtGammaLeg[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCPi0PtGammaLeg[iCut]);
          fHistoMCPi0WOWeightPtGammaLeg[iCut]       = new TH2F("MC_Pi0_WOWeights_PtGamma_Leg", "MC_Pi0_WOWeights_PtGamma_Leg", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt, 3, -0.5, 2.5);
          fMCList[iCut]->Add(fHistoMCPi0WOWeightPtGammaLeg[iCut]);

          fHistoMCPi0InAccPtGammaLeg[iCut]          = new TH2F("MC_Pi0InAcc_PtGamma_Leg", "MC_Pi0InAcc_PtGamma_Leg", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt, 3, -0.5, 2.5);
          fHistoMCPi0InAccPtGammaLeg[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCPi0InAccPtGammaLeg[iCut]);
          fHistoMCPi0WOWeightInAccPtGammaLeg[iCut]  = new TH2F("MC_Pi0WOWeightInAcc_PtGamma_Leg", "MC_Pi0WOWeightInAcc_PtGamma_Leg", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt, 3, -0.5, 2.5);
          fMCList[iCut]->Add(fHistoMCPi0WOWeightInAccPtGammaLeg[iCut]);

          fHistoMCSecPi0PtGamma1vsSource[iCut]      = new TH2F("MC_SecPi0_PtGamma1_Source", "MC_SecPi0_PtGamma1_Source", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt, 16, -0.5, 15.5);
          fMCList[iCut]->Add(fHistoMCSecPi0PtGamma1vsSource[iCut]);
          fHistoMCSecPi0PtGamma2vsSource[iCut]      = new TH2F("MC_SecPi0_PtGamma2_Source", "MC_SecPi0_PtGamma2_Source", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt, 16, -0.5, 15.5);
          fMCList[iCut]->Add(fHistoMCSecPi0PtGamma2vsSource[iCut]);
          fHistoMCSecPi0InAccPtGamma1vsSource[iCut] = new TH2F("MC_SecPi0InAcc_PtGamma1_Source", "MC_SecPi0InAcc_PtGamma1_Source", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt, 16, -0.5, 15.5);
          fMCList[iCut]->Add(fHistoMCSecPi0InAccPtGamma1vsSource[iCut]);
          fHistoMCSecPi0InAccPtGamma2vsSource[iCut] = new TH2F("MC_SecPi0InAcc_PtGamma2_Source", "MC_SecPi0InAcc_PtGamma2_Source", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt, 16, -0.5, 15.5);
          fMCList[iCut]->Add(fHistoMCSecPi0InAccPtGamma2vsSource[iCut]);

          if (fIsMC > 1){
            fHistoMCPi0WOWeightPtGammaLeg[iCut]->Sumw2();
            fHistoMCPi0WOWeightInAccPtGammaLeg[iCut]->Sumw2();
            fHistoMCSecPi0PtGamma1vsSource[iCut]->Sumw2();
            fHistoMCSecPi0PtGamma2vsSource[iCut]->Sumw2();
            fHistoMCSecPi0InAccPtGamma1vsSource[iCut]->Sumw2();
            fHistoMCSecPi0InAccPtGamma2vsSource[iCut]->Sumw2();
          }
        }

        // book additional MC QA histograms
        if (fDoMesonQA > 0){
          fHistoMCPi0PtY[iCut]            = new TH2F("MC_Pi0_Pt_Y", "MC_Pi0_Pt_Y", nBinsQAPt, arrQAPtBinning, 150, -1.5, 1.5);
          fHistoMCPi0PtY[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCPi0PtY[iCut]);
          fHistoMCEtaPtY[iCut]            = new TH2F("MC_Eta_Pt_Y", "MC_Eta_Pt_Y", nBinsQAPt, arrQAPtBinning, 150, -1.5, 1.5);
          fHistoMCEtaPtY[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCEtaPtY[iCut]);
          fHistoMCPi0PtAlpha[iCut]        = new TH2F("MC_Pi0_Pt_Alpha", "MC_Pi0_Pt_Alpha", nBinsQAPt, arrQAPtBinning, 200, -1, 1);
          fMCList[iCut]->Add(fHistoMCPi0PtAlpha[iCut]);
          fHistoMCEtaPtAlpha[iCut]        = new TH2F("MC_Eta_Pt_Alpha", "MC_Eta_Pt_Alpha", nBinsQAPt, arrQAPtBinning, 200, -1, 1);
          fMCList[iCut]->Add(fHistoMCEtaPtAlpha[iCut]);

          if (fIsMC == 2){
            fHistoMCPi0PtAlpha[iCut]->Sumw2();
            fHistoMCEtaPtAlpha[iCut]->Sumw2();
          }
        }
      }

      fTrueList[iCut]                           = new TList();
      fTrueList[iCut]->SetName(Form("%s_%s_%s_%s True histograms",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
      fTrueList[iCut]->SetOwner(kTRUE);
      fCutFolder[iCut]->Add(fTrueList[iCut]);

      if(fDoJetAnalysis && !fDoLightOutput){
        fTrueJetHistograms[iCut] = new TList();
        fTrueJetHistograms[iCut]->SetName(Form("%s_%s_%s_%s True Jet histograms", cutstringEvent.Data(),cutstringPhoton.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
        fTrueJetHistograms[iCut]->SetOwner(kTRUE);
        fCutFolder[iCut]->Add(fTrueJetHistograms[iCut]);
        fHistoTruevsRecJetPt[iCut] = new TH2F("True_JetPt_vs_Rec_JetPt", "True_JetPt_vs_Rec_JetPt", 150, 0, 150, 150, 0, 150);
        fTrueJetHistograms[iCut]->Add(fHistoTruevsRecJetPt[iCut]);
        fHistoTruePi0JetMotherInvMassPt[iCut] = new TH2F("ESD_TruePi0_Jet_InvMass_Pt", "ESD_TruePi0_Jet_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoTruePi0JetMotherInvMassPt[iCut]);
        fHistoTruePrimaryPi0JetInvMassPt[iCut] = new TH2F("ESD_TruePrimaryPi0Jet_InvMass_Pt", "ESD_TruePrimaryPi0Jet_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoTruePrimaryPi0JetInvMassPt[iCut]);
        fHistoTruePi0InJetMotherInvMassPt[iCut] = new TH2F("ESD_TruePi0_Pi0inJet_InvMass_Pt", "ESD_TruePi0_Pi0inJet_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoTruePi0InJetMotherInvMassPt[iCut]);
        fHistoTruePrimaryPi0inJetInvMassPt[iCut] = new TH2F("ESD_TruePrimaryPi0inJet_InvMass_Pt", "ESD_TruePrimaryPi0inJet_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoTruePrimaryPi0inJetInvMassPt[iCut]);
        fHistoTruePrimaryPi0InJetInvMassTruePt[iCut] = new TH2F("ESD_TruePrimaryPi0inJet_InvMass_TruePt", "ESD_TruePrimaryPi0inJet_InvMass_TruePt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoTruePrimaryPi0InJetInvMassTruePt[iCut]);
        fHistoTrueDoubleCountingPi0Jet[iCut] = new TH1F("Double_Counting_True_Pi0inJet", "Double_Counting_True_Pi0inJet", 6, 0, 6);
        fTrueJetHistograms[iCut]->Add(fHistoTrueDoubleCountingPi0Jet[iCut]);
        fHistoTruePi0JetFragmFunc[iCut] = new TH2F("ESD_TruePi0inJetPt_FragmentationFunc", "ESD_TruePi0inJetPt_FragmentationFunc", 50, arrLogBinning, 150, 0., 150.);
        fTrueJetHistograms[iCut]->Add(fHistoTruePi0JetFragmFunc[iCut]);
        fHistoTruePi0JetFragmFuncZInvMass[iCut] = new TH2F("ESD_TruePi0inJetPt_Fragm_Z_InvMass", "ESD_TruePi0inJetPt_Fragm_Z_InvMass", nBinsMinv, 0, maxMinv, 50, arrLogBinning);
        fTrueJetHistograms[iCut]->Add(fHistoTruePi0JetFragmFuncZInvMass[iCut]);
        fHistoMCPi0JetInAccPt[iCut]      = new TH1F("MC_Pi0JetInAcc_Pt", "MC_Pi0JetInAcc_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fTrueJetHistograms[iCut]->Add(fHistoMCPi0JetInAccPt[iCut]);
        fHistoMCPi0inJetInAccPt[iCut]      = new TH1F("MC_Pi0inJetInAcc_Pt", "MC_Pi0inJetInAcc_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fTrueJetHistograms[iCut]->Add(fHistoMCPi0inJetInAccPt[iCut]);
        fHistoMCPi0JetEventGenerated[iCut]    = new TH1F("MC_Pi0_JetEvent_Generated", "MC_Pi0_JetEvent_Generated", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fTrueJetHistograms[iCut]->Add(fHistoMCPi0JetEventGenerated[iCut]);
        fHistoMCPi0inJetGenerated[iCut]    = new TH1F("MC_Pi0_inJet_Generated", "MC_Pi0_inJet_Generated", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fTrueJetHistograms[iCut]->Add(fHistoMCPi0inJetGenerated[iCut]);
        fHistoTrueSecondaryPi0FromK0sJetInvMassPt[iCut] = new TH2F("ESD_TrueSecondaryPi0FromK0sJet_InvMass_Pt", "ESD_TrueSecondaryPi0FromK0sJet_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoTrueSecondaryPi0FromK0sJetInvMassPt[iCut]);
        fHistoTrueSecondaryPi0FromK0sinJetInvMassPt[iCut] = new TH2F("ESD_TrueSecondaryPi0FromK0s_inJet_InvMass_Pt", "ESD_TrueSecondaryPi0FromK0s_inJet_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoTrueSecondaryPi0FromK0sinJetInvMassPt[iCut]);
        fHistoTrueSecondaryPi0FromLambdaJetInvMassPt[iCut] = new TH2F("ESD_TrueSecondaryPi0FromLambdaJet_InvMass_Pt", "ESD_TrueSecondaryPi0FromLambdaJet_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoTrueSecondaryPi0FromLambdaJetInvMassPt[iCut]);
        fHistoTrueSecondaryPi0FromLambdainJetInvMassPt[iCut] = new TH2F("ESD_TrueSecondaryPi0FromLambda_inJet_InvMass_Pt", "ESD_TrueSecondaryPi0FromLambda_inJet_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoTrueSecondaryPi0FromLambdainJetInvMassPt[iCut]);
        fHistoTrueSecondaryPi0FromK0lJetInvMassPt[iCut] = new TH2F("ESD_TrueSecondaryPi0FromK0lJet_InvMass_Pt", "ESD_TrueSecondaryPi0FromK0lJet_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoTrueSecondaryPi0FromK0lJetInvMassPt[iCut]);
        fHistoTrueSecondaryPi0FromK0linJetInvMassPt[iCut] = new TH2F("ESD_TrueSecondaryPi0FromK0l_inJet_InvMass_Pt", "ESD_TrueSecondaryPi0FromK0l_inJet_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoTrueSecondaryPi0FromK0linJetInvMassPt[iCut]);
        fHistoTrueSecondaryPi0InvJetMassPt[iCut] = new TH2F("ESD_TrueSecondaryPi0Jet_InvMass_Pt", "ESD_TrueSecondaryPi0Jet_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoTrueSecondaryPi0InvJetMassPt[iCut]);
        fHistoTrueSecondaryPi0InvinJetMassPt[iCut] = new TH2F("ESD_TrueSecondaryPi0_inJet_InvMass_Pt", "ESD_TrueSecondaryPi0_inJet_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoTrueSecondaryPi0InvinJetMassPt[iCut]);
        fHistoMotherPi0inJetPtY[iCut] = new TH2F("ESD_MotherPi0inJet_Pt_Y", "ESD_MotherPi0inJet_Pt_Y", nBinsQAPt, arrQAPtBinning, 150, -1.5, 1.5);
        fTrueJetHistograms[iCut]->Add(fHistoMotherPi0inJetPtY[iCut]);
        fHistoMotherPi0inJetPtPhi[iCut] = new TH2F("ESD_MotherPi0inJet_Pt_Phi", "ESD_MotherPi0inJet_Pt_Phi", nBinsQAPt, arrQAPtBinning, 150, 0, 6.5);
        fTrueJetHistograms[iCut]->Add(fHistoMotherPi0inJetPtPhi[iCut]);
        fHistoTrueEtaJetMotherInvMassPt[iCut] = new TH2F("ESD_TrueEta_Jet_InvMass_Pt", "ESD_TruePi0_Jet_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoTrueEtaJetMotherInvMassPt[iCut]);
        fHistoTruePrimaryEtaJetInvMassPt[iCut] = new TH2F("ESD_TruePrimaryEtaJet_InvMass_Pt", "ESD_TruePrimaryEtaJet_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoTruePrimaryEtaJetInvMassPt[iCut]);
        fHistoTrueEtaInJetMotherInvMassPt[iCut] = new TH2F("ESD_TrueEta_EtainJet_InvMass_Pt", "ESD_TruePi0_Pi0inJet_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoTrueEtaInJetMotherInvMassPt[iCut]);
        fHistoTruePrimaryEtainJetInvMassPt[iCut] = new TH2F("ESD_TruePrimaryEtainJet_InvMass_Pt", "ESD_TruePrimaryEtainJet_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoTruePrimaryEtainJetInvMassPt[iCut]);
        fHistoTrueDoubleCountingEtaJet[iCut] = new TH1F("Double_Counting_True_EtainJet", "Double_Counting_True_EtainJet", 6, 0, 6);
        fTrueJetHistograms[iCut]->Add(fHistoTrueDoubleCountingEtaJet[iCut]);
        fHistoTrueEtaJetFragmFunc[iCut] = new TH2F("ESD_TrueEtainJet_FragmentationFunc", "ESD_TrueEtainJet_FragmentationFunc", 50, arrLogBinning, 150, 0., 150.);
        fTrueJetHistograms[iCut]->Add(fHistoTrueEtaJetFragmFunc[iCut]);
        fHistoTrueEtaJetFragmFuncZInvMass[iCut] = new TH2F("ESD_TrueEtainJetPt_Fragm_Z_InvMass", "ESD_TrueEtainJetPt_Fragm_Z_InvMass", nBinsMinv, 0, maxMinv, 50, arrLogBinning);
        fTrueJetHistograms[iCut]->Add(fHistoTrueEtaJetFragmFuncZInvMass[iCut]);
        fHistoMCEtaJetInAccPt[iCut]      = new TH1F("MC_EtaJetInAcc_Pt", "MC_EtaJetInAcc_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fTrueJetHistograms[iCut]->Add(fHistoMCEtaJetInAccPt[iCut]);
        fHistoMCEtainJetInAccPt[iCut]      = new TH1F("MC_EtainJetInAcc_Pt", "MC_EtainJetInAcc_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fTrueJetHistograms[iCut]->Add(fHistoMCEtainJetInAccPt[iCut]);
        fHistoMCEtaJetEventGenerated[iCut]    = new TH1F("MC_Eta_JetEvent_Generated", "MC_Eta_JetEvent_Generated", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fTrueJetHistograms[iCut]->Add(fHistoMCEtaJetEventGenerated[iCut]);
        fHistoMCEtainJetGenerated[iCut]    = new TH1F("MC_Eta_inJet_Generated", "MC_Eta_inJet_Generated", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fTrueJetHistograms[iCut]->Add(fHistoMCEtainJetGenerated[iCut]);
        fHistoMotherEtainJetPtY[iCut] = new TH2F("ESD_MotherEtainJet_Pt_Y", "ESD_MotherEtainJet_Pt_Y", nBinsQAPt, arrQAPtBinning, 150, -1.5, 1.5);
        fTrueJetHistograms[iCut]->Add(fHistoMotherEtainJetPtY[iCut]);
        fHistoMotherEtainJetPtPhi[iCut] = new TH2F("ESD_MotherEtainJet_Pt_Phi", "ESD_MotherEtainJet_Pt_Phi", nBinsQAPt, arrQAPtBinning, 150, 0, 6.5);
        fTrueJetHistograms[iCut]->Add(fHistoMotherEtainJetPtPhi[iCut]);
      }
      if(fDoJetQA){
        if(fDoLightOutput){
          fTrueJetHistograms[iCut] = new TList();
          fTrueJetHistograms[iCut]->SetName(Form("%s_%s_%s_%s True Jet histograms", cutstringEvent.Data(),cutstringPhoton.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
          fTrueJetHistograms[iCut]->SetOwner(kTRUE);
          fCutFolder[iCut]->Add(fTrueJetHistograms[iCut]);
        }
        fHistoUnfoldingAsData[iCut]      = new TH2F("Unfolding_AsData", "Unfolding_AsData", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoUnfoldingAsData[iCut]);
        fHistoUnfoldingMissed[iCut]      = new TH2F("Unfolding_Missed", "Unfolding_Missed", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoUnfoldingMissed[iCut]);
        fHistoUnfoldingReject[iCut]      = new TH2F("Unfolding_Reject", "Unfolding_Reject", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoUnfoldingReject[iCut]);
        fHistoUnfoldingAsDataInvMassZ[iCut]      = new TH2F("Unfolding_AsData_InvMass_Z", "Unfolding_AsData_InvMass_Z", nBinsMinv, 0, maxMinv, 50, arrLogBinning);
        fTrueJetHistograms[iCut]->Add(fHistoUnfoldingAsDataInvMassZ[iCut]);
        fHistoUnfoldingMissedInvMassZ[iCut]      = new TH2F("Unfolding_Missed_InvMass_Z", "Unfolding_Missed_InvMass_Z", nBinsMinv, 0, maxMinv, 50, arrLogBinning);
        fTrueJetHistograms[iCut]->Add(fHistoUnfoldingMissedInvMassZ[iCut]);
        fHistoUnfoldingRejectInvMassZ[iCut]      = new TH2F("Unfolding_Reject_InvMass_Z", "Unfolding_Reject_InvMass_Z", nBinsMinv, 0, maxMinv, 50, arrLogBinning);
        fTrueJetHistograms[iCut]->Add(fHistoUnfoldingRejectInvMassZ[iCut]);
      }

      if(!fDoLightOutput){
        fHistoTrueConvGammaPt[iCut]               = new TH1F("ESD_TrueConvGamma_Pt", "ESD_TrueConvGamma_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fHistoTrueConvGammaPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
        fTrueList[iCut]->Add(fHistoTrueConvGammaPt[iCut]);

        fHistoDoubleCountTrueConvGammaRPt[iCut]   = new TH2F("ESD_TrueDoubleCountConvGamma_R_Pt", "ESD_TrueDoubleCountConvGamma_R_Pt", 800, 0, 200, (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fHistoDoubleCountTrueConvGammaRPt[iCut]->SetXTitle("R (cm)");
        fHistoDoubleCountTrueConvGammaRPt[iCut]->SetYTitle("p_{T,clus} (GeV/c)");
        fTrueList[iCut]->Add(fHistoDoubleCountTrueConvGammaRPt[iCut]);
        fHistoMultipleCountTrueConvGamma[iCut]    = new TH1F("ESD_TrueMultipleCountConvGamma", "ESD_TrueMultipleCountConvGamma", 10, 1, 11);
        fTrueList[iCut]->Add(fHistoMultipleCountTrueConvGamma[iCut]);

        fHistoCombinatorialPt[iCut]               = new TH2F("ESD_TrueCombinatorial_Pt", "ESD_TrueCombinatorial_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt, 16, -0.5, 15.5);
        fHistoCombinatorialPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
        fHistoCombinatorialPt[iCut]->SetYTitle("origin");
        fHistoCombinatorialPt[iCut]->GetYaxis()->SetBinLabel( 1,"Elec+Elec");
        fHistoCombinatorialPt[iCut]->GetYaxis()->SetBinLabel( 2,"Elec+Pion");
        fHistoCombinatorialPt[iCut]->GetYaxis()->SetBinLabel( 3,"Elec+Kaon");
        fHistoCombinatorialPt[iCut]->GetYaxis()->SetBinLabel( 4,"Elec+Proton");
        fHistoCombinatorialPt[iCut]->GetYaxis()->SetBinLabel( 5,"Elec+Muon");
        fHistoCombinatorialPt[iCut]->GetYaxis()->SetBinLabel( 6,"Pion+Pion");
        fHistoCombinatorialPt[iCut]->GetYaxis()->SetBinLabel( 7,"Pion+Kaon");
        fHistoCombinatorialPt[iCut]->GetYaxis()->SetBinLabel( 8,"Pion+Proton");
        fHistoCombinatorialPt[iCut]->GetYaxis()->SetBinLabel( 9,"Pion+Muon");
        fHistoCombinatorialPt[iCut]->GetYaxis()->SetBinLabel(10,"Kaon+Kaon");
        fHistoCombinatorialPt[iCut]->GetYaxis()->SetBinLabel(11,"Kaon+Proton");
        fHistoCombinatorialPt[iCut]->GetYaxis()->SetBinLabel(12,"Kaon+Muon");
        fHistoCombinatorialPt[iCut]->GetYaxis()->SetBinLabel(13,"Proton+Proton");
        fHistoCombinatorialPt[iCut]->GetYaxis()->SetBinLabel(14,"Proton+Muon");
        fHistoCombinatorialPt[iCut]->GetYaxis()->SetBinLabel(15,"Muon+Muon");
        fHistoCombinatorialPt[iCut]->GetYaxis()->SetBinLabel(16,"Rest");
        fTrueList[iCut]->Add(fHistoCombinatorialPt[iCut]);

        fHistoTruePrimaryConvGammaPt[iCut]        = new TH1F("ESD_TruePrimaryConvGamma_Pt", "ESD_TruePrimaryConvGamma_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fHistoTruePrimaryConvGammaPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
        fTrueList[iCut]->Add(fHistoTruePrimaryConvGammaPt[iCut]);
        fHistoTrueSecondaryConvGammaPt[iCut]      = new TH2F("ESD_TrueSecondaryConvGamma_Pt", "ESD_TrueSecondaryConvGamma_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt, 4, -0.5, 3.5);
        fHistoTrueSecondaryConvGammaPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
        fHistoTrueSecondaryConvGammaPt[iCut]->SetYTitle("source");
        fHistoTrueSecondaryConvGammaPt[iCut]->GetYaxis()->SetBinLabel(1,"K0s");
        fHistoTrueSecondaryConvGammaPt[iCut]->GetYaxis()->SetBinLabel(2,"K0l");
        fHistoTrueSecondaryConvGammaPt[iCut]->GetYaxis()->SetBinLabel(3,"Lambda");
        fHistoTrueSecondaryConvGammaPt[iCut]->GetYaxis()->SetBinLabel(4,"rest");
        fTrueList[iCut]->Add(fHistoTrueSecondaryConvGammaPt[iCut]);
        fHistoTrueSecondaryConvGammaMCPt[iCut]    = new TH2F("ESD_TrueSecondaryConvGamma_MCPt", "ESD_TrueSecondaryConvGamma_MCPt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt, 4, -0.5, 3.5);
        fHistoTrueSecondaryConvGammaMCPt[iCut]->SetXTitle("p_{T, MC, clus} (GeV/c)");
        fHistoTrueSecondaryConvGammaMCPt[iCut]->SetYTitle("source");
        fHistoTrueSecondaryConvGammaMCPt[iCut]->GetYaxis()->SetBinLabel(1,"K0s");
        fHistoTrueSecondaryConvGammaMCPt[iCut]->GetYaxis()->SetBinLabel(2,"K0l");
        fHistoTrueSecondaryConvGammaMCPt[iCut]->GetYaxis()->SetBinLabel(3,"Lambda");
        fHistoTrueSecondaryConvGammaMCPt[iCut]->GetYaxis()->SetBinLabel(4,"rest");
        fTrueList[iCut]->Add(fHistoTrueSecondaryConvGammaMCPt[iCut]);

        fHistoTruePrimaryConvGammaESDPtMCPt[iCut] = new TH2F("ESD_TruePrimaryConvGammaESD_PtMCPt", "ESD_TruePrimaryConvGammaESD_PtMCPt",
                                                             (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt, (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fHistoTruePrimaryConvGammaESDPtMCPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
        fHistoTruePrimaryConvGammaESDPtMCPt[iCut]->SetYTitle("p_{T, MC, clus} (GeV/c)");
        fTrueList[iCut]->Add(fHistoTruePrimaryConvGammaESDPtMCPt[iCut]);

        fHistoTrueSecondaryConvGammaFromXFromK0sMCPtESDPt[iCut]     = new TH2F("ESD_TrueSecondaryConvGammaFromXFromK0sESD_MCPtPt", "ESD_TrueSecondaryConvGammaFromXFromK0sESD_MCPtPt",
                                                                               (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt, (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fHistoTrueSecondaryConvGammaFromXFromK0sMCPtESDPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
        fHistoTrueSecondaryConvGammaFromXFromK0sMCPtESDPt[iCut]->SetYTitle("p_{T, MC, clus} (GeV/c)");
        fTrueList[iCut]->Add(fHistoTrueSecondaryConvGammaFromXFromK0sMCPtESDPt[iCut]);
        fHistoTrueSecondaryConvGammaFromXFromK0lMCPtESDPt[iCut]     = new TH2F("ESD_TrueSecondaryConvGammaFromXFromK0lESD_MCPtPt", "ESD_TrueSecondaryConvGammaFromXFromK0lESD_MCPtPt",
                                                                               (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt, (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fHistoTrueSecondaryConvGammaFromXFromK0lMCPtESDPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
        fHistoTrueSecondaryConvGammaFromXFromK0lMCPtESDPt[iCut]->SetYTitle("p_{T, MC, clus} (GeV/c)");
        fTrueList[iCut]->Add(fHistoTrueSecondaryConvGammaFromXFromK0lMCPtESDPt[iCut]);
        fHistoTrueSecondaryConvGammaFromXFromLambdaMCPtESDPt[iCut]  = new TH2F("ESD_TrueSecondaryConvGammaFromXFromLambdaESD_MCPtPt", "ESD_TrueSecondaryConvGammaFromXFromLambdaESD_MCPtPt",
                                                                               (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt, (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fHistoTrueSecondaryConvGammaFromXFromLambdaMCPtESDPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
        fHistoTrueSecondaryConvGammaFromXFromLambdaMCPtESDPt[iCut]->SetYTitle("p_{T, MC, clus} (GeV/c)");
        fTrueList[iCut]->Add(fHistoTrueSecondaryConvGammaFromXFromLambdaMCPtESDPt[iCut]);
      }

      fHistoTrueClusGammaPt[iCut]               = new TH1F("TrueClusGamma_Pt", "ESD_TrueClusGamma_Pt", nBinsClusterPt, arrClusPtBinning);
      fHistoTrueClusGammaPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
      fClusterOutputList[iCut]->Add(fHistoTrueClusGammaPt[iCut]);
      if(!fDoLightOutput){
        fHistoTruePrimaryClusGammaPt[iCut]        = new TH1F("TruePrimaryClusGamma_Pt", "ESD_TruePrimaryClusGamma_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTruePrimaryClusGammaPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
        fClusterOutputList[iCut]->Add(fHistoTruePrimaryClusGammaPt[iCut]);
        fHistoTruePrimaryClusGammaESDPtMCPt[iCut] = new TH2F("TruePrimaryClusGamma_Pt_MCPt", "ESD_TruePrimaryClusGamma_Pt_MCPt", nBinsClusterPt, arrClusPtBinning, nBinsClusterPt, arrClusPtBinning);
        fHistoTruePrimaryClusGammaESDPtMCPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
        fHistoTruePrimaryClusGammaESDPtMCPt[iCut]->SetYTitle("p_{T, MC, clus} (GeV/c)");
        fClusterOutputList[iCut]->Add(fHistoTruePrimaryClusGammaESDPtMCPt[iCut]);
        fHistoTrueClusElectronPt[iCut]            = new TH1F("TrueClusElectron_Pt", "TrueElectronGamma_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTrueClusElectronPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
        fClusterOutputList[iCut]->Add(fHistoTrueClusElectronPt[iCut]);
        fHistoTrueClusConvGammaPt[iCut]           = new TH1F("TrueClusConvGamma_Pt", "TrueClusConvGamma_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTrueClusConvGammaPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
        fClusterOutputList[iCut]->Add(fHistoTrueClusConvGammaPt[iCut]);
        fHistoTruePrimaryClusConvGammaPt[iCut]    = new TH1F("TruePrimaryClusConvGamma_Pt", "ESD_TruePrimaryClusConvGamma_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTruePrimaryClusConvGammaPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
        fClusterOutputList[iCut]->Add(fHistoTruePrimaryClusConvGammaPt[iCut]);
        fHistoTruePrimaryClusConvGammaESDPtMCPt[iCut]       = new TH2F("TruePrimaryClusConvGamma_Pt_MCPt", "ESD_TruePrimaryClusConvGamma_Pt_MCPt", nBinsClusterPt, arrClusPtBinning, nBinsClusterPt, arrClusPtBinning);
        fHistoTruePrimaryClusConvGammaESDPtMCPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
        fHistoTruePrimaryClusConvGammaESDPtMCPt[iCut]->SetYTitle("p_{T, MC, clus} (GeV/c)");
        fClusterOutputList[iCut]->Add(fHistoTruePrimaryClusConvGammaESDPtMCPt[iCut]);

        fHistoTrueSecondaryClusGammaPt[iCut]        = new TH1F("TrueSecondaryClusGamma_Pt", "ESD_TrueSecondaryClusGamma_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTrueSecondaryClusGammaPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
        fClusterOutputList[iCut]->Add(fHistoTrueSecondaryClusGammaPt[iCut]);
        fHistoTrueSecondaryClusGammaFromK0sPt[iCut]        = new TH1F("TrueSecondaryClusGammaFromK0s_Pt", "ESD_TrueSecondaryClusGammaFromK0s_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTrueSecondaryClusGammaFromK0sPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
        fClusterOutputList[iCut]->Add(fHistoTrueSecondaryClusGammaFromK0sPt[iCut]);
        fHistoTrueSecondaryClusGammaFromK0lPt[iCut]        = new TH1F("TrueSecondaryClusGammaFromK0l_Pt", "ESD_TrueSecondaryClusGammaFromK0l_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTrueSecondaryClusGammaFromK0lPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
        fClusterOutputList[iCut]->Add(fHistoTrueSecondaryClusGammaFromK0lPt[iCut]);
        fHistoTrueSecondaryClusGammaFromLambdaPt[iCut]        = new TH1F("TrueSecondaryClusGammaFromLambda_Pt", "ESD_TrueSecondaryClusGammaFromLambda_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTrueSecondaryClusGammaFromLambdaPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
        fClusterOutputList[iCut]->Add(fHistoTrueSecondaryClusGammaFromLambdaPt[iCut]);
      }

      fHistoDoubleCountTrueClusterGammaPt[iCut] = new TH2F("TrueDoubleCountClusterGamma_Pt", "TrueDoubleCountClusterGamma_Pt", nBinsClusterPt, arrClusPtBinning, 2, 0, 2);
      fHistoDoubleCountTrueClusterGammaPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
      fClusterOutputList[iCut]->Add(fHistoDoubleCountTrueClusterGammaPt[iCut]);
      fHistoMultipleCountTrueClusterGamma[iCut] = new TH1F("TrueMultipleCountClusterGamma", "TrueMultipleCountClusterGamma", 10, 1, 11);
      fClusterOutputList[iCut]->Add(fHistoMultipleCountTrueClusterGamma[iCut]);
      if(!fDoLightOutput) {
        fHistoTrueNLabelsInClusPt[iCut]           = new TH2F("TrueNLabelsInClus_Pt", "TrueNLabelsInClus_Pt", 100, -0.5, 99.5, nBinsClusterPt, arrClusPtBinning);
        fHistoTrueNLabelsInClusPt[iCut]->SetYTitle("p_{T,clus} (GeV/c)");
        fClusterOutputList[iCut]->Add(fHistoTrueNLabelsInClusPt[iCut]);
        fHistoTrueClusGammaEM02[iCut]             = new TH2F("TrueClusGammaEM02", "TrueClusGammaEM02", nBinsClusterPt, arrClusPtBinning, 400, 0, 5);
        fClusterOutputList[iCut]->Add(fHistoTrueClusGammaEM02[iCut]);
        fHistoTrueClusPi0EM02[iCut]               = new TH2F("TrueClusPi0EM02", "TrueClusPi0EM02", nBinsClusterPt, arrClusPtBinning, 400, 0, 5);
        fClusterOutputList[iCut]->Add(fHistoTrueClusPi0EM02[iCut]);
        fHistoTrueClusEMNonLeadingPt[iCut]        = new TH1F("TrueClusEMNonLeading_Pt", "TrueClusEMNonLeading_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTrueClusEMNonLeadingPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
        fClusterOutputList[iCut]->Add(fHistoTrueClusEMNonLeadingPt[iCut]);
      }

      if ( (fDoMaterialBudgetWeightingOfGammasForTrueMesons && fIsMC > 0) ) {
        if(!fDoLightOutput){
          fHistoTrueConvGammaPt[iCut]->Sumw2();
          fHistoDoubleCountTrueConvGammaRPt[iCut]->Sumw2();
          fHistoMultipleCountTrueConvGamma[iCut]->Sumw2();
          fHistoCombinatorialPt[iCut]->Sumw2();
          fHistoTruePrimaryConvGammaPt[iCut]->Sumw2();
          fHistoTrueSecondaryConvGammaPt[iCut]->Sumw2();
          fHistoTrueSecondaryConvGammaMCPt[iCut]->Sumw2();
          fHistoTruePrimaryConvGammaESDPtMCPt[iCut]->Sumw2();
          fHistoTrueSecondaryConvGammaFromXFromK0sMCPtESDPt[iCut]->Sumw2();
          fHistoTrueSecondaryConvGammaFromXFromK0lMCPtESDPt[iCut]->Sumw2();
          fHistoTrueSecondaryConvGammaFromXFromLambdaMCPtESDPt[iCut]->Sumw2();
        }
      }

      if (fIsMC > 1){
        if(!fDoLightOutput){
          fHistoTrueConvGammaPt[iCut]->Sumw2();
          fHistoDoubleCountTrueConvGammaRPt[iCut]->Sumw2();
          fHistoMultipleCountTrueConvGamma[iCut]->Sumw2();
          fHistoCombinatorialPt[iCut]->Sumw2();
          fHistoTruePrimaryConvGammaPt[iCut]->Sumw2();
          fHistoTrueSecondaryConvGammaPt[iCut]->Sumw2();
          fHistoTrueSecondaryConvGammaMCPt[iCut]->Sumw2();
          fHistoTruePrimaryConvGammaESDPtMCPt[iCut]->Sumw2();
          fHistoTrueSecondaryConvGammaFromXFromK0sMCPtESDPt[iCut]->Sumw2();
          fHistoTrueSecondaryConvGammaFromXFromK0lMCPtESDPt[iCut]->Sumw2();
          fHistoTrueSecondaryConvGammaFromXFromLambdaMCPtESDPt[iCut]->Sumw2();
        }
        fHistoTrueClusGammaPt[iCut]->Sumw2();
        fHistoDoubleCountTrueClusterGammaPt[iCut]->Sumw2();
        fHistoMultipleCountTrueClusterGamma[iCut]->Sumw2();
        if(!fDoLightOutput) {
          fHistoTruePrimaryClusGammaPt[iCut]->Sumw2();
          fHistoTruePrimaryClusGammaESDPtMCPt[iCut]->Sumw2();
          fHistoTrueNLabelsInClusPt[iCut]->Sumw2();
          fHistoTrueClusGammaEM02[iCut]->Sumw2();
          fHistoTrueClusPi0EM02[iCut]->Sumw2();
          fHistoTrueClusElectronPt[iCut]->Sumw2();
          fHistoTrueClusConvGammaPt[iCut]->Sumw2();
          fHistoTruePrimaryClusConvGammaPt[iCut]->Sumw2();
          fHistoTruePrimaryClusConvGammaESDPtMCPt[iCut]->Sumw2();
          fHistoTrueClusEMNonLeadingPt[iCut]->Sumw2();
          fHistoTrueSecondaryClusGammaPt[iCut]->Sumw2();
          fHistoTrueSecondaryClusGammaFromK0sPt[iCut]->Sumw2();
          fHistoTrueSecondaryClusGammaFromK0lPt[iCut]->Sumw2();
          fHistoTrueSecondaryClusGammaFromLambdaPt[iCut]->Sumw2();
        }
      }

      if(!fDoLightOutput){
        fHistoTruePi0InvMassECalib[iCut]          = new TH2F("True_Pi0_InvMass_E_Calib", "True_Pi0_InvMass_E_Calib", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
        fHistoTruePi0InvMassECalib[iCut]->SetXTitle("M_{inv} (GeV/c^{2})");
        fHistoTruePi0InvMassECalib[iCut]->SetYTitle("E_{cluster}(GeV)");
        fESDList[iCut]->Add(fHistoTruePi0InvMassECalib[iCut]);
        fHistoTruePi0PureGammaInvMassECalib[iCut] = new TH2F("True_Pi0PureGamma_InvMass_E_Calib", "True_Pi0PureGamma_InvMass_E_Calib", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
        fHistoTruePi0PureGammaInvMassECalib[iCut]->SetXTitle("M_{inv} (GeV/c^{2})");
        fHistoTruePi0PureGammaInvMassECalib[iCut]->SetYTitle("E_{cluster}(GeV)");
        fESDList[iCut]->Add(fHistoTruePi0PureGammaInvMassECalib[iCut]);
        fHistoTruePi0InvMassECalibPCM[iCut]          = new TH2F("True_Pi0_InvMass_E_Calib_PCM", "True_Pi0_InvMass_E_Calib_PCM", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
        fHistoTruePi0InvMassECalibPCM[iCut]->SetXTitle("M_{inv} (GeV/c^{2})");
        fHistoTruePi0InvMassECalibPCM[iCut]->SetYTitle("E_{#gamma_{PCM}}(GeV)");
        fESDList[iCut]->Add(fHistoTruePi0InvMassECalibPCM[iCut]);
      }

      if (fDoClusterQA > 0 && fDoLightOutput!=2){
        fHistoTrueClusConvGammaFullyPt[iCut]        = new TH1F("TrueClusConvGammaFullyContained_Pt", "TrueClusConvGammaFullyContained_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTrueClusConvGammaFullyPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
        fClusterOutputList[iCut]->Add(fHistoTrueClusConvGammaFullyPt[iCut]);
        fHistoTrueClusMergedGammaPt[iCut]           = new TH1F("TrueClusMergedGamma_Pt", "TrueClusMergedGamma_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTrueClusMergedGammaPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
        fClusterOutputList[iCut]->Add(fHistoTrueClusMergedGammaPt[iCut]);
        fHistoTrueClusMergedPartConvGammaPt[iCut]   = new TH1F("TrueClusMergedPartConvGamma_Pt", "TrueClusMergedPartConvGamma_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTrueClusMergedPartConvGammaPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
        fClusterOutputList[iCut]->Add(fHistoTrueClusMergedPartConvGammaPt[iCut]);
        fHistoTrueClusDalitzPt[iCut]                = new TH1F("TrueClusDalitz_Pt", "TrueClusDalitz_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTrueClusDalitzPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
        fClusterOutputList[iCut]->Add(fHistoTrueClusDalitzPt[iCut]);
        fHistoTrueClusDalitzMergedPt[iCut]          = new TH1F("TrueClusDalitzMerged_Pt", "TrueClusDalitzMerged_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTrueClusDalitzMergedPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
        fClusterOutputList[iCut]->Add(fHistoTrueClusDalitzMergedPt[iCut]);
        fHistoTrueClusPhotonFromElecMotherPt[iCut]  = new TH1F("TrueClusPhotonFromElecMother_Pt", "TrueClusPhotonFromElecMother_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTrueClusPhotonFromElecMotherPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
        fClusterOutputList[iCut]->Add(fHistoTrueClusPhotonFromElecMotherPt[iCut]);
        fHistoTrueClusShowerPt[iCut]                = new TH1F("TrueClusShower_Pt", "TrueClusShower_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTrueClusShowerPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
        fClusterOutputList[iCut]->Add(fHistoTrueClusShowerPt[iCut]);
        fHistoTrueClusSubLeadingPt[iCut]            = new TH1F("TrueClusSubleading_Pt", "TrueClusSubleading_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTrueClusSubLeadingPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
        fClusterOutputList[iCut]->Add(fHistoTrueClusSubLeadingPt[iCut]);
        fHistoTrueClusNMothers[iCut]              = new TH1F("TrueClusNParticles", "TrueClusNParticles", 20, 0, 20);
        fHistoTrueClusNMothers[iCut]->SetXTitle("# particles in cluster");
        fClusterOutputList[iCut]->Add(fHistoTrueClusNMothers[iCut]);
        if (fIsMC > 1){
          fHistoTrueClusConvGammaFullyPt[iCut]->Sumw2();
          fHistoTrueClusMergedGammaPt[iCut]->Sumw2();
          fHistoTrueClusMergedPartConvGammaPt[iCut]->Sumw2();
          fHistoTrueClusDalitzPt[iCut]->Sumw2();
          fHistoTrueClusDalitzMergedPt[iCut]->Sumw2();
          fHistoTrueClusPhotonFromElecMotherPt[iCut]->Sumw2();
          fHistoTrueClusShowerPt[iCut]->Sumw2();
          fHistoTrueClusSubLeadingPt[iCut]->Sumw2();
          fHistoTrueClusNMothers[iCut]->Sumw2();
        }
      }

      if(fDoMesonAnalysis){
        fHistoTruePi0InvMassPt[iCut]                = new TH2F("ESD_TruePi0_InvMass_Pt", "ESD_TruePi0_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
        fHistoTruePi0InvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0}}(GeV/c^{2})");
        fHistoTruePi0InvMassPt[iCut]->SetYTitle("#pi^{0}p_{T}(GeV/c)");
        fTrueList[iCut]->Add(fHistoTruePi0InvMassPt[iCut]);
        fHistoTrueEtaInvMassPt[iCut]                =  new TH2F("ESD_TrueEta_InvMass_Pt", "ESD_TrueEta_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
        fHistoTrueEtaInvMassPt[iCut]->SetXTitle("M_{inv,#eta}(GeV/c^{2})");
        fHistoTrueEtaInvMassPt[iCut]->SetYTitle("#eta p_{T}(GeV/c)");
        fTrueList[iCut]->Add(fHistoTrueEtaInvMassPt[iCut]);
        if (fDoLightOutput!=2 ){
            fHistoTruePi0InvMassPtAdditional[iCut]                = new TH2F("ESD_TruePi0_InvMass_Pt_Additional", "ESD_TruePi0_InvMass_Pt_Additional", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
            fHistoTruePi0InvMassPtAdditional[iCut]->SetXTitle("M_{inv,#pi^{0}}(GeV/c^{2})");
            fHistoTruePi0InvMassPtAdditional[iCut]->SetYTitle("#pi^{0}p_{T}(GeV/c)");
            fTrueList[iCut]->Add(fHistoTruePi0InvMassPtAdditional[iCut]);
            fHistoTrueEtaInvMassPtAdditional[iCut]                =  new TH2F("ESD_TrueEta_InvMass_Pt_Additional", "ESD_TrueEta_InvMass_PtAdditional", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
            fHistoTrueEtaInvMassPtAdditional[iCut]->SetXTitle("M_{inv,#eta}(GeV/c^{2})");
            fHistoTrueEtaInvMassPtAdditional[iCut]->SetYTitle("#eta p_{T}(GeV/c)");
            fTrueList[iCut]->Add(fHistoTrueEtaInvMassPtAdditional[iCut]);
            fHistoTruePi0MatchedInvMassPt[iCut]         = new TH2F("ESD_TruePi0_Matched_InvMass_Pt", "ESD_TruePi0_Matched_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
            fHistoTruePi0MatchedInvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0}}(GeV/c^{2})");
            fHistoTruePi0MatchedInvMassPt[iCut]->SetYTitle("#pi^{0}p_{T}(GeV/c)");
            fTrueList[iCut]->Add(fHistoTruePi0MatchedInvMassPt[iCut]);
            fHistoTrueEtaMatchedInvMassPt[iCut]         = new TH2F("ESD_TrueEta_Matched_InvMass_Pt", "ESD_TrueEta_Matched_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
            fHistoTrueEtaMatchedInvMassPt[iCut]->SetXTitle("M_{inv,#eta}(GeV/c^{2})");
            fHistoTrueEtaMatchedInvMassPt[iCut]->SetYTitle("#eta p_{T}(GeV/c)");
            fTrueList[iCut]->Add(fHistoTrueEtaMatchedInvMassPt[iCut]);

            fHistoDoubleCountTruePi0InvMassPt[iCut]     = new TH2F("ESD_TrueDoubleCountPi0_InvMass_Pt", "ESD_TrueDoubleCountPi0_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
            fHistoDoubleCountTruePi0InvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0}}(GeV/c^{2})");
            fHistoDoubleCountTruePi0InvMassPt[iCut]->SetYTitle("#pi^{0} p_{T}(GeV/c)");
            fTrueList[iCut]->Add(fHistoDoubleCountTruePi0InvMassPt[iCut]);
            fHistoMultipleCountTruePi0[iCut]            = new TH1F("ESD_TrueMultipleCountPi0", "ESD_TrueMultipleCountPi0", 10, 1, 11);
            fHistoMultipleCountTruePi0[iCut]->SetXTitle("# number of multiple counts");
            fTrueList[iCut]->Add(fHistoMultipleCountTruePi0[iCut]);
            fHistoDoubleCountTrueEtaInvMassPt[iCut]     = new TH2F("ESD_TrueDoubleCountEta_InvMass_Pt", "ESD_TrueDoubleCountEta_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
            fHistoDoubleCountTrueEtaInvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0}}(GeV/c^{2})");
            fHistoDoubleCountTrueEtaInvMassPt[iCut]->SetYTitle("#eta p_{T}(GeV/c)");
            fTrueList[iCut]->Add(fHistoDoubleCountTrueEtaInvMassPt[iCut]);
            fHistoMultipleCountTrueEta[iCut]            = new TH1F("ESD_TrueMultipleCountEta", "ESD_TrueMultipleCountEta", 10, 1, 11);
            fHistoMultipleCountTrueEta[iCut]->SetXTitle("# number of multiple counts");
            fTrueList[iCut]->Add(fHistoMultipleCountTrueEta[iCut]);
        }

        fHistoTruePrimaryPi0InvMassPt[iCut]         = new TH2F("ESD_TruePrimaryPi0_InvMass_Pt", "ESD_TruePrimaryPi0_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
        fHistoTruePrimaryPi0InvMassPt[iCut]->SetXTitle("M_{inv,prim #pi^{0}}(GeV/c^{2})");
        fHistoTruePrimaryPi0InvMassPt[iCut]->SetYTitle("#pi^{0}p_{T}(GeV/c)");
        fHistoTruePrimaryPi0InvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTruePrimaryPi0InvMassPt[iCut]);

        fHistoTruePrimaryEtaInvMassPt[iCut]         = new TH2F("ESD_TruePrimaryEta_InvMass_Pt", "ESD_TruePrimaryEta_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
        fHistoTruePrimaryEtaInvMassPt[iCut]->SetXTitle("M_{inv,prim #eta}(GeV/c^{2})");
        fHistoTruePrimaryEtaInvMassPt[iCut]->SetYTitle("#eta p_{T}(GeV/c)");
        fHistoTruePrimaryEtaInvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTruePrimaryEtaInvMassPt[iCut]);

        fHistoTruePrimaryPi0W0WeightingInvMassPt[iCut]  = new TH2F("ESD_TruePrimaryPi0W0Weights_InvMass_Pt", "ESD_TruePrimaryPi0W0Weights_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
        fHistoTruePrimaryPi0W0WeightingInvMassPt[iCut]->SetXTitle("M_{inv,prim #pi^{0}}(GeV/c^{2})");
        fHistoTruePrimaryPi0W0WeightingInvMassPt[iCut]->SetYTitle("#pi^{0}p_{T}(GeV/c)");
        fTrueList[iCut]->Add(fHistoTruePrimaryPi0W0WeightingInvMassPt[iCut]);

        fHistoTruePrimaryEtaW0WeightingInvMassPt[iCut]  = new TH2F("ESD_TruePrimaryEtaW0Weights_InvMass_Pt", "ESD_TruePrimaryEtaW0Weights_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
        fHistoTruePrimaryEtaW0WeightingInvMassPt[iCut]->SetXTitle("M_{inv,prim #eta}(GeV/c^{2})");
        fHistoTruePrimaryEtaW0WeightingInvMassPt[iCut]->SetYTitle("#eta p_{T}(GeV/c)");
        fTrueList[iCut]->Add(fHistoTruePrimaryEtaW0WeightingInvMassPt[iCut]);

        fProfileTruePrimaryPi0WeightsInvMassPt[iCut]    = new TProfile2D("ESD_TruePrimaryPi0Weights_InvMass_Pt", "ESD_TruePrimaryPi0Weights_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
        fProfileTruePrimaryPi0WeightsInvMassPt[iCut]->SetXTitle("M_{inv,prim #pi^{0}}(GeV/c^{2})");
        fProfileTruePrimaryPi0WeightsInvMassPt[iCut]->SetYTitle("#pi^{0}p_{T}(GeV/c)");
        fProfileTruePrimaryPi0WeightsInvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fProfileTruePrimaryPi0WeightsInvMassPt[iCut]);
        fProfileTruePrimaryEtaWeightsInvMassPt[iCut]    = new TProfile2D("ESD_TruePrimaryEtaWeights_InvMass_Pt", "ESD_TruePrimaryEtaWeights_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
        fProfileTruePrimaryEtaWeightsInvMassPt[iCut]->SetXTitle("M_{inv,prim #eta}(GeV/c^{2})");
        fProfileTruePrimaryEtaWeightsInvMassPt[iCut]->SetYTitle("#eta p_{T}(GeV/c)");
        fProfileTruePrimaryEtaWeightsInvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fProfileTruePrimaryEtaWeightsInvMassPt[iCut]);

        fHistoTrueSecondaryPi0InvMassPt[iCut]           = new TH2F("ESD_TrueSecondaryPi0_InvMass_Pt", "ESD_TrueSecondaryPi0_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
        fHistoTrueSecondaryPi0InvMassPt[iCut]->SetXTitle("M_{inv,sec #pi^{0}}(GeV/c^{2})");
        fHistoTrueSecondaryPi0InvMassPt[iCut]->SetYTitle("#pi^{0}p_{T}(GeV/c)");
        fHistoTrueSecondaryPi0InvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTrueSecondaryPi0InvMassPt[iCut]);

        fHistoTrueSecondaryPi0FromK0sInvMassPt[iCut]    = new TH2F("ESD_TrueSecondaryPi0FromK0s_InvMass_Pt", "ESD_TrueSecondaryPi0FromK0s_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
        fHistoTrueSecondaryPi0FromK0sInvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0}from K^{0}_{S}}(GeV/c^{2})");
        fHistoTrueSecondaryPi0FromK0sInvMassPt[iCut]->SetYTitle("#pi^{0}p_{T}(GeV/c)");
        fHistoTrueSecondaryPi0FromK0sInvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTrueSecondaryPi0FromK0sInvMassPt[iCut]);
        fHistoTrueSecondaryPi0FromK0lInvMassPt[iCut]    = new TH2F("ESD_TrueSecondaryPi0FromK0l_InvMass_Pt", "ESD_TrueSecondaryPi0FromK0l_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
        fHistoTrueSecondaryPi0FromK0lInvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0}from K^{0}_{L}}(GeV/c^{2})");
        fHistoTrueSecondaryPi0FromK0lInvMassPt[iCut]->SetYTitle("#pi^{0}p_{T}(GeV/c)");
        fHistoTrueSecondaryPi0FromK0lInvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTrueSecondaryPi0FromK0lInvMassPt[iCut]);
        fHistoTrueSecondaryPi0FromEtaInvMassPt[iCut]    = new TH2F("ESD_TrueSecondaryPi0FromEta_InvMass_Pt", "ESD_TrueSecondaryPi0FromEta_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
        fHistoTrueSecondaryPi0FromEtaInvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0}from #eta}(GeV/c^{2})");
        fHistoTrueSecondaryPi0FromEtaInvMassPt[iCut]->SetYTitle("#pi^{0}p_{T}(GeV/c)");
        fTrueList[iCut]->Add(fHistoTrueSecondaryPi0FromEtaInvMassPt[iCut]);
        fHistoTrueSecondaryPi0FromLambdaInvMassPt[iCut] = new TH2F("ESD_TrueSecondaryPi0FromLambda_InvMass_Pt", "ESD_TrueSecondaryPi0FromLambda_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
        fHistoTrueSecondaryPi0FromLambdaInvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0}from #Lambda}(GeV/c^{2})");
        fHistoTrueSecondaryPi0FromLambdaInvMassPt[iCut]->SetYTitle("#pi^{0}p_{T}(GeV/c)");
        fTrueList[iCut]->Add(fHistoTrueSecondaryPi0FromLambdaInvMassPt[iCut]);

        if(!fDoLightOutput) {
          fHistoTruePrimaryPi0PhotonPairPtconv[iCut]            = new TH2F("ESD_TruePrimaryPi0_InvMass_PtConv", "", 400, 0, 0.4, nBinsPt, arrPtBinning);
          fHistoTruePrimaryPi0PhotonPairPtconv[iCut]->SetXTitle("M_{inv,#pi^{0}}(GeV/c^{2})");
          fHistoTruePrimaryPi0PhotonPairPtconv[iCut]->SetYTitle("#gamma^{conv} p_{T} (GeV/c)");
          fHistoTruePrimaryPi0PhotonPairPtconv[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePrimaryPi0PhotonPairPtconv[iCut]);

          fHistoTruePrimaryPi0W0WeightsPhotonPairPtconv[iCut]   = new TH2F("ESD_TruePrimaryPi0W0Weights_InvMass_Pt", "", 400, 0, 0.4, nBinsPt, arrPtBinning);
          fHistoTruePrimaryPi0W0WeightsPhotonPairPtconv[iCut]->SetXTitle("M_{inv,#pi^{0}}(GeV/c^{2})");
          fHistoTruePrimaryPi0W0WeightsPhotonPairPtconv[iCut]->SetYTitle("#gamma^{conv} p_{T} (GeV/c)");
          fTrueList[iCut]->Add(fHistoTruePrimaryPi0W0WeightsPhotonPairPtconv[iCut]);

          fHistoTrueSecondaryPi0PhotonPairPtconv[iCut]          = new TH2F("ESD_TrueSecondaryPi0_InvMass_PtConv", "", 400, 0, 0.4, nBinsPt, arrPtBinning);
          fHistoTrueSecondaryPi0PhotonPairPtconv[iCut]->SetXTitle("M_{inv,#pi^{0}}(GeV/c^{2})");
          fHistoTrueSecondaryPi0PhotonPairPtconv[iCut]->SetYTitle("#gamma^{conv} p_{T} (GeV/c)");
          fHistoTrueSecondaryPi0PhotonPairPtconv[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTrueSecondaryPi0PhotonPairPtconv[iCut]);

          fHistoTrueSecondaryPi0FromK0sPhotonPairPtconv[iCut]   = new TH2F("ESD_TrueSecondaryPi0FromK0s_InvMass_PtConv", "", 400, 0, 0.4, nBinsPt, arrPtBinning);
          fHistoTrueSecondaryPi0FromK0sPhotonPairPtconv[iCut]->SetXTitle("M_{inv,#pi^{0}}(GeV/c^{2})");
          fHistoTrueSecondaryPi0FromK0sPhotonPairPtconv[iCut]->SetYTitle("#gamma^{conv} p_{T} (GeV/c)");
          fHistoTrueSecondaryPi0FromK0sPhotonPairPtconv[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTrueSecondaryPi0FromK0sPhotonPairPtconv[iCut]);

          fHistoTrueSecondaryPi0FromK0lPhotonPairPtconv[iCut]   = new TH2F("ESD_TrueSecondaryPi0FromK0l_InvMass_PtConv", "", 400, 0, 0.4, nBinsPt, arrPtBinning);
          fHistoTrueSecondaryPi0FromK0lPhotonPairPtconv[iCut]->SetXTitle("M_{inv,#pi^{0}}(GeV/c^{2})");
          fHistoTrueSecondaryPi0FromK0lPhotonPairPtconv[iCut]->SetYTitle("#gamma^{conv} p_{T} (GeV/c)");
          fHistoTrueSecondaryPi0FromK0lPhotonPairPtconv[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTrueSecondaryPi0FromK0lPhotonPairPtconv[iCut]);

          fHistoTrueSecondaryPi0FromLambdaPhotonPairPtconv[iCut]= new TH2F("ESD_TrueSecondaryPi0FromLambda_InvMass_PtConv", "", 400, 0, 0.4, nBinsPt, arrPtBinning);
          fHistoTrueSecondaryPi0FromLambdaPhotonPairPtconv[iCut]->SetXTitle("M_{inv,#pi^{0}}(GeV/c^{2})");
          fHistoTrueSecondaryPi0FromLambdaPhotonPairPtconv[iCut]->SetYTitle("#gamma^{conv} p_{T} (GeV/c)");
          fHistoTrueSecondaryPi0FromLambdaPhotonPairPtconv[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTrueSecondaryPi0FromLambdaPhotonPairPtconv[iCut]);

          fHistoTruePrimaryPi0DCPtconv[iCut]                = new TH1F("ESD_TruePrimaryPi0DC_PtConv", "", nBinsPt, arrPtBinning);
          fHistoTruePrimaryPi0DCPtconv[iCut]->SetXTitle("#gamma^{conv} p_{T} (GeV/c)");
          fHistoTruePrimaryPi0DCPtconv[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePrimaryPi0DCPtconv[iCut]);

          fHistoTrueSecondaryPi0DCPtconvSource[iCut]        = new TH2F("ESD_TrueSecondaryPi0DC_PtConv_Source", "", nBinsPt, arrPtBinning, 4, -0.5, 3.5);
          fHistoTrueSecondaryPi0DCPtconvSource[iCut]->SetXTitle("#gamma^{conv} p_{T} (GeV/c)");
          fHistoTrueSecondaryPi0DCPtconvSource[iCut]->SetYTitle("sec. source");
          fHistoTrueSecondaryPi0DCPtconvSource[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTrueSecondaryPi0DCPtconvSource[iCut]);

          fHistoTruePrimaryPi0MissingPtconv[iCut]           = new TH1F("ESD_TruePrimaryPi0Missing_PtConv", "", nBinsPt, arrPtBinning);
          fHistoTruePrimaryPi0MissingPtconv[iCut]->SetXTitle("#gamma^{conv} p_{T} (GeV/c)");
          fHistoTruePrimaryPi0MissingPtconv[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePrimaryPi0MissingPtconv[iCut]);

          fHistoTrueSecondaryPi0MissingPtconvSource[iCut]   = new TH2F("ESD_TrueSecondaryPi0Missing_PtConv_Source", "", nBinsPt, arrPtBinning, 4, -0.5, 3.5);
          fHistoTrueSecondaryPi0MissingPtconvSource[iCut]->SetXTitle("#gamma^{conv} p_{T} (GeV/c)");
          fHistoTrueSecondaryPi0MissingPtconvSource[iCut]->SetYTitle("sec. source");
          fHistoTrueSecondaryPi0MissingPtconvSource[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTrueSecondaryPi0MissingPtconvSource[iCut]);

          // put proper summing for JJ weights
          if (fIsMC > 1){
            fHistoTruePrimaryPi0W0WeightsPhotonPairPtconv[iCut]->Sumw2();
          }
        }

        if(fIsMC>0 && fDoHBTHistoOutput){
          fHistoTrueGammaInvMassPt[iCut]                = new TH2F("ESD_TrueGamma_InvMass_Pt", "ESD_TrueGamma_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
          fHistoTrueGammaInvMassPt[iCut]->SetXTitle("M_{inv}(GeV/c^{2})");
          fHistoTrueGammaInvMassPt[iCut]->SetYTitle("p_{T}(GeV/c)");
          fHistoTrueGammaInvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTrueGammaInvMassPt[iCut]);

          fHistoBckTrueGammaInvMassPt[iCut]                = new TH2F("ESD_TrueGammaBck_InvMass_Pt", "ESD_TrueGammaBck_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
          fHistoBckTrueGammaInvMassPt[iCut]->SetXTitle("M_{inv}(GeV/c^{2})");
          fHistoBckTrueGammaInvMassPt[iCut]->SetYTitle("p_{T}(GeV/c)");
          fHistoBckTrueGammaInvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoBckTrueGammaInvMassPt[iCut]);

          fHistoMCGammaHBTInvMassPt[iCut]                = new TH2F("ESD_GenGamma_InvMass_Pt", "ESD_GenGamma_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
          fHistoMCGammaHBTInvMassPt[iCut]->SetXTitle("M_{inv}(GeV/c^{2})");
          fHistoMCGammaHBTInvMassPt[iCut]->SetYTitle("p_{T}(GeV/c)");
          fHistoMCGammaHBTInvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoMCGammaHBTInvMassPt[iCut]);

          fHistoBckMCGammaHBTInvMassPt[iCut]                = new TH2F("ESD_GenGammaBck_InvMass_Pt", "ESD_GenGammaBck_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
          fHistoBckMCGammaHBTInvMassPt[iCut]->SetXTitle("M_{inv}(GeV/c^{2})");
          fHistoBckMCGammaHBTInvMassPt[iCut]->SetYTitle("p_{T}(GeV/c)");
          fHistoBckMCGammaHBTInvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoBckMCGammaHBTInvMassPt[iCut]);

          fHistoMCGammaHBTOpeningAnglePt[iCut]                = new TH2F("ESD_GenGamma_OpeningAngle_Pt", "ESD_GenGamma_OpeningAngle_Pt", 700, 0, 3.5, nBinsPt, arrPtBinning);
          fHistoMCGammaHBTOpeningAnglePt[iCut]->SetXTitle("#theta");
          fHistoMCGammaHBTOpeningAnglePt[iCut]->SetYTitle("p_{T}(GeV/c)");
          fHistoMCGammaHBTOpeningAnglePt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoMCGammaHBTOpeningAnglePt[iCut]);

          fHistoBckMCGammaHBTOpeningAnglePt[iCut]                = new TH2F("ESD_GenGammaBck_OpeningAngle_Pt", "ESD_GenGammaBck_OpeningAngle_Pt", 700, 0, 3.5, nBinsPt, arrPtBinning);
          fHistoBckMCGammaHBTOpeningAnglePt[iCut]->SetXTitle("#theta");
          fHistoBckMCGammaHBTOpeningAnglePt[iCut]->SetYTitle("p_{T}(GeV/c)");
          fHistoBckMCGammaHBTOpeningAnglePt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoBckMCGammaHBTOpeningAnglePt[iCut]);

          fHistoMCGammaHBTDeltaEPt[iCut]                = new TH2F("ESD_GenGamma_DeltaE_Pt", "ESD_GenGamma_DeltaE_Pt", 800, 0, 4.0, nBinsPt, arrPtBinning);
          fHistoMCGammaHBTDeltaEPt[iCut]->SetXTitle("#Delta E");
          fHistoMCGammaHBTDeltaEPt[iCut]->SetYTitle("p_{T}(GeV/c)");
          fHistoMCGammaHBTDeltaEPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoMCGammaHBTDeltaEPt[iCut]);

          fHistoBckMCGammaHBTDeltaEPt[iCut]                = new TH2F("ESD_GenGammaBck_DeltaE_Pt", "ESD_GenGammaBck_DeltaE_Pt", 800, 0, 4.0, nBinsPt, arrPtBinning);
          fHistoBckMCGammaHBTDeltaEPt[iCut]->SetXTitle("#Delta E");
          fHistoBckMCGammaHBTDeltaEPt[iCut]->SetYTitle("p_{T}(GeV/c)");
          fHistoBckMCGammaHBTDeltaEPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoBckMCGammaHBTDeltaEPt[iCut]);
        }

        if (fIsMC > 1){
          fHistoTruePi0InvMassPt[iCut]->Sumw2();
          fHistoTrueEtaInvMassPt[iCut]->Sumw2();
          if (fDoLightOutput!=2 ){
            fHistoTruePi0InvMassPtAdditional[iCut]->Sumw2();
            fHistoTrueEtaInvMassPtAdditional[iCut]->Sumw2();
            fHistoTruePi0MatchedInvMassPt[iCut]->Sumw2();
            fHistoTrueEtaMatchedInvMassPt[iCut]->Sumw2();
            fHistoDoubleCountTruePi0InvMassPt[iCut]->Sumw2();
            fHistoMultipleCountTruePi0[iCut]->Sumw2();
            fHistoDoubleCountTrueEtaInvMassPt[iCut]->Sumw2();
            fHistoMultipleCountTrueEta[iCut]->Sumw2();
          }
          fHistoTruePrimaryPi0W0WeightingInvMassPt[iCut]->Sumw2();
          fHistoTruePrimaryEtaW0WeightingInvMassPt[iCut]->Sumw2();
          fHistoTrueSecondaryPi0FromEtaInvMassPt[iCut]->Sumw2();
          fHistoTrueSecondaryPi0FromLambdaInvMassPt[iCut]->Sumw2();
        }

        if (fDoMesonQA > 0){
          if (fIsMC < 2){
            fHistoTruePi0CaloPhotonInvMassPt[iCut]                  = new TH2F("ESD_TruePi0CaloPhoton_InvMass_Pt", "ESD_TruePi0CaloPhoton_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
            fHistoTruePi0CaloPhotonInvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0}}(GeV/c^{2}) #gamma #gamma");
            fHistoTruePi0CaloPhotonInvMassPt[iCut]->SetYTitle("#pi^{0}p_{T}(GeV/c)");
            fTrueList[iCut]->Add(fHistoTruePi0CaloPhotonInvMassPt[iCut]);

            fHistoTrueEtaCaloPhotonInvMassPt[iCut]                  = new TH2F("ESD_TrueEtaCaloPhoton_InvMass_Pt", "ESD_TrueEtaCaloPhoton_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
            fHistoTrueEtaCaloPhotonInvMassPt[iCut]->SetXTitle("M_{inv,#eta}(GeV/c^{2}) #gamma #gamma");
            fHistoTrueEtaCaloPhotonInvMassPt[iCut]->SetYTitle("#eta p_{T}(GeV/c)");
            fTrueList[iCut]->Add(fHistoTrueEtaCaloPhotonInvMassPt[iCut]);

            fHistoTruePi0CaloConvertedPhotonInvMassPt[iCut]         = new TH2F("ESD_TruePi0CaloConvertedPhoton_InvMass_Pt", "ESD_TruePi0CaloConvertedPhoton_InvMass_Pt",
                                                                               nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
            fHistoTruePi0CaloConvertedPhotonInvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0}}(GeV/c^{2}) #gamma #gamma_{conv}");
            fHistoTruePi0CaloConvertedPhotonInvMassPt[iCut]->SetYTitle("#pi^{0}p_{T}(GeV/c)");
            fTrueList[iCut]->Add(fHistoTruePi0CaloConvertedPhotonInvMassPt[iCut]);

            fHistoTruePi0CaloConvertedPhotonMatchedInvMassPt[iCut]  = new TH2F("ESD_TruePi0CaloConvertedPhotonMatched_InvMass_Pt", "ESD_TruePi0CaloConvertedPhotonMatched_InvMass_Pt",
                                                                               nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
            fHistoTruePi0CaloConvertedPhotonMatchedInvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0}}(GeV/c^{2}) #gamma #gamma_{conv,matched}");
            fHistoTruePi0CaloConvertedPhotonMatchedInvMassPt[iCut]->SetYTitle("#pi^{0}p_{T}(GeV/c)");
            fTrueList[iCut]->Add(fHistoTruePi0CaloConvertedPhotonMatchedInvMassPt[iCut]);

            fHistoTruePi0CaloConvPhotonConvRPt[iCut]                = new TH2F("ESD_TruePi0CaloConvPhoton_ConvR_PtE", "ESD_TruePi0CaloConvPhoton_ConvR_PtE", 920, 0, 460, nBinsPt, arrPtBinning);
            fHistoTruePi0CaloConvPhotonConvRPt[iCut]->SetXTitle("R_{conv,e_{calo}}(cm)");
            fHistoTruePi0CaloConvPhotonConvRPt[iCut]->SetYTitle("e^{#pm}p_{T}(GeV/c)");
            fTrueList[iCut]->Add(fHistoTruePi0CaloConvPhotonConvRPt[iCut]);

            fHistoTruePi0CaloConvPhotonConvRAlpha[iCut]             = new TH2F("ESD_TruePi0CaloConvPhoton_ConvR_AlphaE", "ESD_TruePi0CaloConvPhoton_ConvR_AlphaE", 920, 0, 460, 200, -1, 1);
            fHistoTruePi0CaloConvPhotonConvRAlpha[iCut]->SetXTitle("R_{conv,e_{calo}}(cm)");
            fHistoTruePi0CaloConvPhotonConvRAlpha[iCut]->SetYTitle("#alpha converted e-pair");
            fTrueList[iCut]->Add(fHistoTruePi0CaloConvPhotonConvRAlpha[iCut]);

            fHistoTruePi0CaloConvPhotonPtAlpha[iCut]                = new TH2F("ESD_TruePi0CaloConvPhoton_PtE_AlphaE", "ESD_TruePi0CaloConvPhoton_PtE_AlphaE", nBinsPt, arrPtBinning, 200, -1, 1);
            fHistoTruePi0CaloConvPhotonPtAlpha[iCut]->SetXTitle("e^{#pm}p_{T}(GeV/c)");
            fHistoTruePi0CaloConvPhotonPtAlpha[iCut]->SetYTitle("#alpha converted e-pair");
            fTrueList[iCut]->Add(fHistoTruePi0CaloConvPhotonPtAlpha[iCut]);

            fHistoTrueEtaCaloConvertedPhotonInvMassPt[iCut]         = new TH2F("ESD_TrueEtaCaloConvertedPhoton_InvMass_Pt", "ESD_TrueEtaCaloConvertedPhoton_InvMass_Pt",
                                                                               nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
            fHistoTrueEtaCaloConvertedPhotonInvMassPt[iCut]->SetXTitle("M_{inv,#eta}(GeV/c^{2}) #gamma #gamma_{conv}");
            fHistoTrueEtaCaloConvertedPhotonInvMassPt[iCut]->SetYTitle("#eta p_{T}(GeV/c)");
            fTrueList[iCut]->Add(fHistoTrueEtaCaloConvertedPhotonInvMassPt[iCut]);

            fHistoTrueEtaCaloConvertedPhotonMatchedInvMassPt[iCut]  = new TH2F("ESD_TrueEtaCaloConvertedPhotonMatched_InvMass_Pt", "ESD_TrueEtaCaloConvertedPhotonMatched_InvMass_Pt",
                                                                               nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
            fHistoTrueEtaCaloConvertedPhotonMatchedInvMassPt[iCut]->SetXTitle("M_{inv,#eta}(GeV/c^{2}) #gamma #gamma_{conv,matched}");
            fHistoTrueEtaCaloConvertedPhotonMatchedInvMassPt[iCut]->SetYTitle("#eta p_{T}(GeV/c)");
            fTrueList[iCut]->Add(fHistoTrueEtaCaloConvertedPhotonMatchedInvMassPt[iCut]);

            fHistoTrueEtaCaloConvPhotonConvRPt[iCut]                = new TH2F("ESD_TrueEtaCaloConvPhoton_ConvR_PtE", "ESD_TrueEtaCaloConvPhoton_ConvR_PtE", 920, 0, 460, nBinsPt, arrPtBinning);
            fHistoTrueEtaCaloConvPhotonConvRPt[iCut]->SetXTitle("R_{conv,e_{calo}}(cm)");
            fHistoTrueEtaCaloConvPhotonConvRPt[iCut]->SetYTitle("e^{#pm}p_{T}(GeV/c)");
            fTrueList[iCut]->Add(fHistoTrueEtaCaloConvPhotonConvRPt[iCut]);

            fHistoTrueEtaCaloConvPhotonConvRAlpha[iCut]             = new TH2F("ESD_TrueEtaCaloConvPhoton_ConvR_AlphaE", "ESD_TrueEtaCaloConvPhoton_ConvR_AlphaE", 920, 0, 460, 200, -1, 1);
            fHistoTrueEtaCaloConvPhotonConvRAlpha[iCut]->SetXTitle("R_{conv,e_{calo}}(cm)");
            fHistoTrueEtaCaloConvPhotonConvRAlpha[iCut]->SetYTitle("#alpha converted e-pair");
            fTrueList[iCut]->Add(fHistoTrueEtaCaloConvPhotonConvRAlpha[iCut]);

            fHistoTrueEtaCaloConvPhotonPtAlpha[iCut]                = new TH2F("ESD_TrueEtaCaloConvPhoton_PtE_AlphaE", "ESD_TrueEtaCaloConvPhoton_PtE_AlphaE", nBinsPt, arrPtBinning, 200, -1, 1);
            fHistoTrueEtaCaloConvPhotonPtAlpha[iCut]->SetXTitle("e^{#pm}p_{T}(GeV/c)");
            fHistoTrueEtaCaloConvPhotonPtAlpha[iCut]->SetYTitle("#alpha converted e-pair");
            fTrueList[iCut]->Add(fHistoTrueEtaCaloConvPhotonPtAlpha[iCut]);


            fHistoTruePi0CaloElectronInvMassPt[iCut]                = new TH2F("ESD_TruePi0CaloElectron_InvMass_Pt", "ESD_TruePi0CaloElectron_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
            fHistoTruePi0CaloElectronInvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0}}(GeV/c^{2}) #gamma e^{#pm}");
            fHistoTruePi0CaloElectronInvMassPt[iCut]->SetYTitle("#pi^{0}p_{T}(GeV/c)");
            fTrueList[iCut]->Add(fHistoTruePi0CaloElectronInvMassPt[iCut]);
            fHistoTrueEtaCaloElectronInvMassPt[iCut]                = new TH2F("ESD_TrueEtaCaloElectron_InvMass_Pt", "ESD_TrueEtaCaloElectron_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
            fHistoTrueEtaCaloElectronInvMassPt[iCut]->SetXTitle("M_{inv,#eta}(GeV/c^{2}) #gamma e^{#pm}");
            fHistoTrueEtaCaloElectronInvMassPt[iCut]->SetYTitle("#eta p_{T}(GeV/c)");
            fTrueList[iCut]->Add(fHistoTrueEtaCaloElectronInvMassPt[iCut]);

            fHistoTruePi0CaloMergedClusterInvMassPt[iCut]           = new TH2F("ESD_TruePi0CaloMergedCluster_InvMass_Pt", "ESD_TruePi0CaloMergedCluster_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
            fHistoTruePi0CaloMergedClusterInvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0}}(GeV/c^{2}) #gamma merged cluster");
            fHistoTruePi0CaloMergedClusterInvMassPt[iCut]->SetYTitle("#pi^{0}p_{T}(GeV/c)");
            fTrueList[iCut]->Add(fHistoTruePi0CaloMergedClusterInvMassPt[iCut]);
            fHistoTrueEtaCaloMergedClusterInvMassPt[iCut]           = new TH2F("ESD_TrueEtaCaloMergedCluster_InvMass_Pt", "ESD_TrueEtaCaloMergedCluster_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
            fHistoTrueEtaCaloMergedClusterInvMassPt[iCut]->SetXTitle("M_{inv,#eta}(GeV/c^{2}) #gamma merged cluster");
            fHistoTrueEtaCaloMergedClusterInvMassPt[iCut]->SetYTitle("#eta p_{T}(GeV/c)");
            fTrueList[iCut]->Add(fHistoTrueEtaCaloMergedClusterInvMassPt[iCut]);

            fHistoTrueMotherCaloEMNonLeadingInvMassPt[iCut]         = new TH2F("ESD_TrueMotherCaloEMNonLeading_InvMass_Pt", "ESD_TrueMotherCaloEMNonLeading_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
            fHistoTrueMotherCaloEMNonLeadingInvMassPt[iCut]->SetXTitle("M_{inv} (GeV/c^{2}) #gamma cluster no leading EM");
            fHistoTrueMotherCaloEMNonLeadingInvMassPt[iCut]->SetYTitle("#pair p_{T}(GeV/c)");
            fTrueList[iCut]->Add(fHistoTrueMotherCaloEMNonLeadingInvMassPt[iCut]);
            fHistoTruePi0CaloMergedClusterPartConvInvMassPt[iCut]   = new TH2F("ESD_TruePi0CaloMergedClusterPartConv_InvMass_Pt", "ESD_TruePi0CaloMergedClusterPartConv_InvMass_Pt",
                                                                               nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
            fHistoTruePi0CaloMergedClusterPartConvInvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0}}(GeV/c^{2}) #gamma merged cluster, part conv");
            fHistoTruePi0CaloMergedClusterPartConvInvMassPt[iCut]->SetYTitle("#pi^{0}p_{T}(GeV/c)");
            fTrueList[iCut]->Add(fHistoTruePi0CaloMergedClusterPartConvInvMassPt[iCut]);
            fHistoTrueEtaCaloMergedClusterPartConvInvMassPt[iCut]   = new TH2F("ESD_TrueEtaCaloMergedClusterPartConv_InvMass_Pt", "ESD_TrueEtaCaloMergedClusterPartConv_InvMass_Pt",
                                                                               nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
            fHistoTrueEtaCaloMergedClusterPartConvInvMassPt[iCut]->SetXTitle("M_{inv,#eta}(GeV/c^{2}) #gamma merged cluster, part conv");
            fHistoTrueEtaCaloMergedClusterPartConvInvMassPt[iCut]->SetYTitle("#eta p_{T}(GeV/c)");
            fTrueList[iCut]->Add(fHistoTrueEtaCaloMergedClusterPartConvInvMassPt[iCut]);

            fHistoTruePrimaryPi0MCPtResolPt[iCut]                   = new TH2F("ESD_TruePrimaryPi0_MCPt_ResolPt", "ESD_TruePrimaryPi0_ResolPt_MCPt", 500, 0.03, 25, 1000, -1., 1.);
            fHistoTruePrimaryPi0MCPtResolPt[iCut]->SetXTitle("#pi^{0}p_{T,MC}(GeV/c)");
            fHistoTruePrimaryPi0MCPtResolPt[iCut]->SetYTitle("#pi^{0}(p_{T,rec}-p_{T,MC})/p_{T,MC}()");
            fHistoTruePrimaryPi0MCPtResolPt[iCut]->Sumw2();
            SetLogBinningXTH2(fHistoTruePrimaryPi0MCPtResolPt[iCut]);
            fTrueList[iCut]->Add(fHistoTruePrimaryPi0MCPtResolPt[iCut]);

            fHistoTruePrimaryEtaMCPtResolPt[iCut]                   = new TH2F("ESD_TruePrimaryEta_MCPt_ResolPt", "ESD_TruePrimaryEta_ResolPt_MCPt", 500, 0.03, 25, 1000, -1., 1.);
            fHistoTruePrimaryEtaMCPtResolPt[iCut]->SetXTitle("#eta p_{T,MC}(GeV/c)");
            fHistoTruePrimaryEtaMCPtResolPt[iCut]->SetYTitle("#eta (p_{T,rec}-p_{T,MC})/p_{T,MC}()");
            fHistoTruePrimaryEtaMCPtResolPt[iCut]->Sumw2();
            SetLogBinningXTH2(fHistoTruePrimaryEtaMCPtResolPt[iCut]);
            fTrueList[iCut]->Add(fHistoTruePrimaryEtaMCPtResolPt[iCut]);

            fHistoTrueK0sWithPi0DaughterMCPt[iCut]                  = new TH1F("ESD_TrueK0sWithPi0Daughter_MCPt", "ESD_TrueK0sWithPi0Daughter_MCPt", nBinsPt, arrPtBinning);
            fHistoTrueK0sWithPi0DaughterMCPt[iCut]->SetXTitle("K^{0}_{s}p_{MC,T}(GeV/c) for K^{0}_{s}where #pi^{0}rec ");
            fTrueList[iCut]->Add(fHistoTrueK0sWithPi0DaughterMCPt[iCut]);
            fHistoTrueK0lWithPi0DaughterMCPt[iCut]                  = new TH1F("ESD_TrueK0lWithPi0Daughter_MCPt", "ESD_TrueK0lWithPi0Daughter_MCPt", nBinsPt, arrPtBinning);
            fHistoTrueK0lWithPi0DaughterMCPt[iCut]->SetXTitle("K^{0}_{s}p_{MC,T}(GeV/c) for K^{0}_{l}where #pi^{0}rec ");
            fTrueList[iCut]->Add(fHistoTrueK0lWithPi0DaughterMCPt[iCut]);
            fHistoTrueEtaWithPi0DaughterMCPt[iCut]                  = new TH1F("ESD_TrueEtaWithPi0Daughter_MCPt", "ESD_TrueEtaWithPi0Daughter_MCPt", nBinsPt, arrPtBinning);
            fHistoTrueEtaWithPi0DaughterMCPt[iCut]->SetXTitle("#eta p_{MC,T}(GeV/c) for #eta where #pi^{0}rec ");
            fTrueList[iCut]->Add(fHistoTrueEtaWithPi0DaughterMCPt[iCut]);
            fHistoTrueLambdaWithPi0DaughterMCPt[iCut]               = new TH1F("ESD_TrueLambdaWithPi0Daughter_MCPt", "ESD_TrueLambdaWithPi0Daughter_MCPt", nBinsPt, arrPtBinning);
            fHistoTrueLambdaWithPi0DaughterMCPt[iCut]->SetXTitle("#Lambda p_{MC,T}(GeV/c) for #Lambda where #pi^{0}rec ");
            fTrueList[iCut]->Add(fHistoTrueLambdaWithPi0DaughterMCPt[iCut]);
          }
          if(fDoMesonQA>1){
            fHistoTrueBckGGInvMassPt[iCut]                          = new TH2F("ESD_TrueBckGG_InvMass_Pt", "ESD_TrueBckGG_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
            fHistoTrueBckGGInvMassPt[iCut]->SetXTitle("M_{inv} (GeV/c^{2}) #gamma #gamma no signal");
            fHistoTrueBckGGInvMassPt[iCut]->SetYTitle("#pair p_{T}(GeV/c)");
            fTrueList[iCut]->Add(fHistoTrueBckGGInvMassPt[iCut]);
            fHistoTrueBckFullMesonContainedInOneClusterInvMassPt[iCut] = new TH2F("ESD_TrueBckFullMesonContained_InvMass_Pt", "ESD_TrueBckFullMesonContained_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
            fHistoTrueBckFullMesonContainedInOneClusterInvMassPt[iCut]->SetXTitle("M_{inv} (GeV/c^{2}) #gamma #gamma, calo gamma with full pi0");
            fHistoTrueBckFullMesonContainedInOneClusterInvMassPt[iCut]->SetYTitle("#pair p_{T}(GeV/c)");
            fTrueList[iCut]->Add(fHistoTrueBckFullMesonContainedInOneClusterInvMassPt[iCut]);
            fHistoTrueBckAsymEClustersInvMassPt[iCut]   = new TH2F("ESD_TrueBckAsymEClus_InvMass_Pt", "ESD_TrueBckAsymEClus_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
            fHistoTrueBckAsymEClustersInvMassPt[iCut]->SetXTitle("M_{inv} (GeV/c^{2}) #gamma #gamma, calo gamma >70% of pi0 energy");
            fHistoTrueBckAsymEClustersInvMassPt[iCut]->SetYTitle("#pair p_{T}(GeV/c)");
            fTrueList[iCut]->Add(fHistoTrueBckAsymEClustersInvMassPt[iCut]);
            fHistoTrueBckContInvMassPt[iCut]                        = new TH2F("ESD_TrueBckCont_InvMass_Pt", "ESD_TrueBckCont_InvMass_Pt", nBinsMinv, minMinv, maxMinv, nBinsPt, arrPtBinning);
            fHistoTrueBckContInvMassPt[iCut]->SetXTitle("M_{inv} (GeV/c^{2}) contamination");
            fHistoTrueBckContInvMassPt[iCut]->SetYTitle("#pair p_{T}(GeV/c)");
            fTrueList[iCut]->Add(fHistoTrueBckContInvMassPt[iCut]);
          }
          fHistoTruePi0PtY[iCut]                      = new TH2F("ESD_TruePi0_Pt_Y", "ESD_TruePi0_Pt_Y", nBinsQAPt, arrQAPtBinning, 150, -1.5, 1.5);
          fHistoTruePi0PtY[iCut]->SetYTitle("Y_{#pi^{0}}");
          fHistoTruePi0PtY[iCut]->SetXTitle("#pi^{0}p_{T}(GeV/c)");
          fTrueList[iCut]->Add(fHistoTruePi0PtY[iCut]);
          fHistoTrueEtaPtY[iCut]                      = new TH2F("ESD_TrueEta_Pt_Y", "ESD_TrueEta_Pt_Y", nBinsQAPt, arrQAPtBinning, 150, -1.5, 1.5);
          fHistoTrueEtaPtY[iCut]->SetYTitle("Y_{#eta}");
          fHistoTrueEtaPtY[iCut]->SetXTitle("#eta p_{T}(GeV/c)");
          fTrueList[iCut]->Add(fHistoTrueEtaPtY[iCut]);
          fHistoTruePi0PtAlpha[iCut]                  = new TH2F("ESD_TruePi0_Pt_Alpha", "ESD_TruePi0_Pt_Alpha", nBinsQAPt, arrQAPtBinning, 200, -1, 1);
          fHistoTruePi0PtAlpha[iCut]->SetYTitle("#alpha_{#pi^{0}}");
          fHistoTruePi0PtAlpha[iCut]->SetXTitle("#pi^{0}p_{T}(GeV/c)");
          fTrueList[iCut]->Add(fHistoTruePi0PtAlpha[iCut]);
          fHistoTrueEtaPtAlpha[iCut]                  = new TH2F("ESD_TrueEta_Pt_Alpha", "ESD_TrueEta_Pt_Alpha", nBinsQAPt, arrQAPtBinning, 200, -1, 1);
          fHistoTrueEtaPtAlpha[iCut]->SetYTitle("#alpha_{#eta}");
          fHistoTrueEtaPtAlpha[iCut]->SetXTitle("#eta p_{T}(GeV/c)");
          fTrueList[iCut]->Add(fHistoTrueEtaPtAlpha[iCut]);

          fHistoTruePi0PtOpenAngle[iCut]              = new TH2F("ESD_TruePi0_Pt_OpenAngle", "ESD_TruePi0_Pt_OpenAngle", nBinsQAPt, arrQAPtBinning, 100, 0, 1);
          fHistoTruePi0PtOpenAngle[iCut]->SetYTitle("#theta_{#pi^{0}}");
          fHistoTruePi0PtOpenAngle[iCut]->SetXTitle("#pi^{0}p_{T}(GeV/c)");
          fTrueList[iCut]->Add(fHistoTruePi0PtOpenAngle[iCut]);
          fHistoTrueEtaPtOpenAngle[iCut]              = new TH2F("ESD_TrueEta_Pt_OpenAngle", "ESD_TrueEta_Pt_OpenAngle", nBinsQAPt, arrQAPtBinning, 200, 0, TMath::Pi());
          fHistoTrueEtaPtOpenAngle[iCut]->SetYTitle("#theta_{#eta}");
          fHistoTrueEtaPtOpenAngle[iCut]->SetXTitle("#eta p_{T}(GeV/c)");
          fTrueList[iCut]->Add(fHistoTrueEtaPtOpenAngle[iCut]);

          fHistoTrueMotherPi0ConvPhotonEtaPhi[iCut]   = new TH2F("ESD_TrueMotherPi0ConvPhoton_Eta_Phi", "conv photons for true #pi^{0}", 600, 0, 2*TMath::Pi(), 200, -1, 1);
          fHistoTrueMotherPi0ConvPhotonEtaPhi[iCut]->SetXTitle("#phi_{#gamma_{conv}}(rad)");
          fHistoTrueMotherPi0ConvPhotonEtaPhi[iCut]->SetYTitle("#eta_{#gamma_{conv}}");
          fTrueList[iCut]->Add(fHistoTrueMotherPi0ConvPhotonEtaPhi[iCut]);
          fHistoTrueMotherEtaConvPhotonEtaPhi[iCut]   = new TH2F("ESD_TrueMotherEtaConvPhoton_Eta_Phi", "conv photons for true #eta", 600, 0, 2*TMath::Pi(), 200, -1, 1);
          fHistoTrueMotherEtaConvPhotonEtaPhi[iCut]->SetXTitle("#phi_{#gamma_{conv}}(rad)");
          fHistoTrueMotherEtaConvPhotonEtaPhi[iCut]->SetYTitle("#eta_{#gamma_{conv}}");
          fTrueList[iCut]->Add(fHistoTrueMotherEtaConvPhotonEtaPhi[iCut]);

          if (fIsMC > 1){
            if(fDoMesonQA>1){
              fHistoTrueBckGGInvMassPt[iCut]->Sumw2();
              fHistoTrueBckFullMesonContainedInOneClusterInvMassPt[iCut]->Sumw2();
              fHistoTrueBckAsymEClustersInvMassPt[iCut]->Sumw2();
              fHistoTrueBckContInvMassPt[iCut]->Sumw2();
            }
            fHistoTruePi0PtY[iCut]->Sumw2();
            fHistoTrueEtaPtY[iCut]->Sumw2();
            fHistoTruePi0PtAlpha[iCut]->Sumw2();
            fHistoTrueEtaPtAlpha[iCut]->Sumw2();
            fHistoTruePi0PtOpenAngle[iCut]->Sumw2();
            fHistoTrueEtaPtOpenAngle[iCut]->Sumw2();
            fHistoTrueMotherPi0ConvPhotonEtaPhi[iCut]->Sumw2();
            fHistoTrueMotherEtaConvPhotonEtaPhi[iCut]->Sumw2();
          }
        }
      }
    }
  }

  fVectorDoubleCountTruePi0s.clear();
  fVectorDoubleCountTrueEtas.clear();
  fVectorDoubleCountTrueConvGammas.clear();
  fVectorDoubleCountTrueClusterGammas.clear();

  fMapMultipleCountTruePi0s.clear();
  fMapMultipleCountTrueEtas.clear();
  fMapMultipleCountTrueConvGammas.clear();
  fMapMultipleCountTrueClusterGammas.clear();

  fVectorRecTruePi0s.clear();
  fVectorRecTrueEtas.clear();

  if(fV0Reader) {
    if((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())
      if(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms())
        fOutputContainer->Add(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms());
    if((AliConvEventCuts*)fV0Reader->GetEventCuts())
      if(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms())
        fOutputContainer->Add(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms());
    if(fV0Reader->GetProduceV0FindingEfficiency())
      if (fV0Reader->GetV0FindingEfficiencyHistograms())
        fOutputContainer->Add(fV0Reader->GetV0FindingEfficiencyHistograms());
  }

  for(Int_t iMatcherTask = 0; iMatcherTask < 5; iMatcherTask++){
    AliCaloTrackMatcher* temp = 0x0;
    if(!fCorrTaskSetting.CompareTo("")){
      temp = (AliCaloTrackMatcher*) (AliAnalysisManager::GetAnalysisManager()->GetTask(Form("CaloTrackMatcher_%i_%i",iMatcherTask,fTrackMatcherRunningMode)));
    } else {
      temp = (AliCaloTrackMatcher*) (AliAnalysisManager::GetAnalysisManager()->GetTask(Form("CaloTrackMatcher_%i_%i_%s",iMatcherTask,fTrackMatcherRunningMode,fCorrTaskSetting.Data())));
    }
    if(temp) fOutputContainer->Add(temp->GetCaloTrackMatcherHistograms());
  }

  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if(!((AliConvEventCuts*)fEventCutArray->At(iCut))) continue;
    if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms()){
      fCutFolder[iCut]->Add(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms());
    }
    if(!((AliConversionPhotonCuts*)fCutArray->At(iCut))) continue;
    if(((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutHistograms()){
      fCutFolder[iCut]->Add(((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutHistograms());
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
    if(fDoJetAnalysis){
      fCutFolder[iCut]->Add(fJetHistograms[iCut]);
    }
  }

  if (fIsMC > 0 && fDoLightOutput!=2){
    tBrokenFiles = new TTree("BrokenFiles", "BrokenFiles");
    tBrokenFiles->Branch("fileName",&fFileNameBroken);
    fOutputContainer->Add(tBrokenFiles);
  }


  PostData(1, fOutputContainer);
  Int_t nContainerOutput = 2;
  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if(fDoPhotonQA == 2){
        OpenFile(nContainerOutput);
        PostData(nContainerOutput, fTreeConvGammaPtDcazCat[iCut]);
        nContainerOutput++;
    }
    if(fDoMesonQA == 2){
      OpenFile(nContainerOutput);
      PostData(nContainerOutput, tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut]);
      nContainerOutput++;
    }
  }
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskGammaConvCalo::Notify()
{
  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetPeriodEnum() == AliConvEventCuts::kNoPeriod && ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum() != AliConvEventCuts::kNoPeriod){
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->SetPeriodEnumExplicit(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum());
    } else if (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetPeriodEnum() == AliConvEventCuts::kNoPeriod ){
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->SetPeriodEnum(fV0Reader->GetPeriodName());
    }

    if(fIsHeavyIon == 2){
      if(!((AliConvEventCuts*)fEventCutArray->At(iCut))->GetDoEtaShift()){
        fProfileEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
        continue; // No Eta Shift requested, continue
      }
      if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift() == 0.0){ // Eta Shift requested but not set, get shift automatically
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCorrectEtaShiftFromPeriod();
        fProfileEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once
        continue;
      } else{
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
void AliAnalysisTaskGammaConvCalo::UserExec(Option_t *)
{
  //
  // Called for each event
  //
  fInputEvent = InputEvent();

  if(fIsMC > 0) fMCEvent = MCEvent();

  Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();
  if(fInputEvent->IsIncompleteDAQ()==kTRUE) eventQuality = 2;// incomplete event
  // Event Not Accepted due to MC event missing or wrong trigger for V0ReaderV1 or because it is incomplete => abort processing of this event/file
  if(eventQuality == 2 || eventQuality == 3){
    // write out name of broken file for first event
    if (fIsMC > 0 && fDoLightOutput!=2){
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
      if (fIsMC > 1) fHistoNEventsWOWeight[iCut]->Fill(eventQuality);
    }
    return;
  }


//   if(fInputEvent->IsA()==AliAODEvent::Class()){
//     fInputEvent->InitMagneticField();
//   }

  fReaderGammas = fV0Reader->GetReconstructedGammas(); // Gammas from default Cut

  // ------------------- BeginEvent ----------------------------
  AliEventplane *EventPlane = fInputEvent->GetEventplane();
  if(fIsHeavyIon ==1)fEventPlaneAngle = EventPlane->GetEventplane("V0",fInputEvent,2);
  else fEventPlaneAngle=0.0;

  if(fIsMC>0 && fInputEvent->IsA()==AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled())){
    RelabelAODPhotonCandidates(kTRUE);// In case of AODMC relabeling MC
    fV0Reader->RelabelAODs(kTRUE);
  }

  for(Int_t iCut = 0; iCut<fnCuts; iCut++){

    fiCut = iCut;
//     cout << ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber() << "_" <<  ((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutNumber() <<
//             "_" << ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber() <<endl;

    Bool_t isRunningEMCALrelAna = kFALSE;
    if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 1) isRunningEMCALrelAna = kTRUE;

    Int_t eventNotAccepted = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon,isRunningEMCALrelAna);

    if(fIsMC==2){
      Float_t xsection      = -1.;
      Float_t ntrials       = -1.;
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetXSectionAndNTrials(fMCEvent,xsection,ntrials,fInputEvent);
      if((xsection==-1.) || (ntrials==-1.)) AliFatal(Form("ERROR: GetXSectionAndNTrials returned invalid xsection/ntrials, periodName from V0Reader: '%s'",fV0Reader->GetPeriodName().Data()));
      fProfileJetJetXSection[iCut]->Fill(0.,xsection);
      fHistoJetJetNTrials[iCut]->Fill("#sum{NTrials}",ntrials);
    }

    if(fIsMC>0){
      fWeightJetJetMC       = 1;
  //     cout << fMCEvent << endl;
      Float_t maxjetpt      = -1.;
      Float_t pthard = -1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetUseJetFinderForOutliers()) maxjetpt = fOutlierJetReader->GetMaxJetPt();
      Bool_t isMCJet        = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsJetJetMCEventAccepted( fMCEvent, fWeightJetJetMC, pthard,fInputEvent, maxjetpt);
      if (fIsMC == 3){
        Double_t weightMult   = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetWeightForMultiplicity(fV0Reader->GetNumberOfPrimaryTracks());
        fWeightJetJetMC       = fWeightJetJetMC*weightMult;
      }
      if(fIsMC==1) fWeightJetJetMC = 1;
      if (!isMCJet){
        fHistoNEvents[iCut]->Fill(10,fWeightJetJetMC);
        if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(10);
        continue;
      }
    }

    Bool_t triggered = kTRUE;

    if(eventNotAccepted!= 0){
      fHistoNEvents[iCut]->Fill(eventNotAccepted, fWeightJetJetMC); // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1
      if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventNotAccepted);
      // cout << "event rejected due to wrong trigger: " <<eventNotAccepted << endl;
      if (eventNotAccepted==3 && fIsMC > 0){
        triggered = kFALSE;
      }else {
        continue;
      }
    }

    if(eventQuality != 0 && triggered== kTRUE){// Event Not Accepted
      //cout << "event rejected due to: " <<eventQuality << endl;
      fHistoNEvents[iCut]->Fill(eventQuality, fWeightJetJetMC);
      if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventQuality);
      continue;
    }

    if (triggered==kTRUE){
      if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetUseSphericity()!=0){
        if(fV0Reader->GetSphericity() != -1 && fV0Reader->GetSphericity() != 0){
          fHistoEventSphericity[iCut]->Fill(fV0Reader->GetSphericity(), fWeightJetJetMC);
          fHistoNEvents[iCut]->Fill(eventQuality, fWeightJetJetMC); // Should be 0 here
          if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventQuality); // Should be 0 here

          if(fHistoNGoodESDTracks && fHistoNGoodESDTracks[iCut]) fHistoNGoodESDTracks[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks(), fWeightJetJetMC);
          if(fHistoVertexZ && fHistoVertexZ[iCut]) fHistoVertexZ[iCut]->Fill(fInputEvent->GetPrimaryVertex()->GetZ(), fWeightJetJetMC);
        }
      } else {
        fHistoNEvents[iCut]->Fill(eventQuality, fWeightJetJetMC); // Should be 0 here
        if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventQuality); // Should be 0 here
        if(fHistoNGoodESDTracks && fHistoNGoodESDTracks[iCut]) fHistoNGoodESDTracks[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks(), fWeightJetJetMC);
        if(fHistoVertexZ && fHistoVertexZ[iCut]) fHistoVertexZ[iCut]->Fill(fInputEvent->GetPrimaryVertex()->GetZ(), fWeightJetJetMC);
        if(!fDoLightOutput){
          fHistoVertexX[iCut]->Fill(fInputEvent->GetPrimaryVertex()->GetX(), fWeightJetJetMC);
          fHistoVertexY[iCut]->Fill(fInputEvent->GetPrimaryVertex()->GetY(), fWeightJetJetMC);
          fHistoSPDClusterTrackletBackground[iCut]->Fill(fInputEvent->GetMultiplicity()->GetNumberOfTracklets(),(fInputEvent->GetNumberOfITSClusters(0)+fInputEvent->GetNumberOfITSClusters(1)),fWeightJetJetMC);
          ((AliConvEventCuts*)fEventCutArray->At(iCut))->FillTPCPileUpHistograms(fInputEvent);
          if(((AliConvEventCuts*)fEventCutArray->At(iCut))->IsHeavyIon() == 2)  fHistoNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A(), fWeightJetJetMC);
          else fHistoNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A()+fInputEvent->GetVZEROData()->GetMTotV0C(), fWeightJetJetMC);
        }
      }
    }

      if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetDoElecDeDxPostCalibration()){
        if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->LoadElecDeDxPostCalibration(fInputEvent->GetRunNumber())){
          AliFatal(Form("ERROR: LoadElecDeDxPostCalibration returned kFALSE for %d despite being requested!",fInputEvent->GetRunNumber()));
        }
      }

    if(fIsMC>0){
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
    if(fIsMC>0){
      if(fInputEvent->IsA()==AliESDEvent::Class())
        ProcessMCParticles();
      if(fInputEvent->IsA()==AliAODEvent::Class())
        ProcessAODMCParticles();
    }

    if (triggered==kFALSE) continue;

    // it is in the loop to have the same conversion cut string (used also for MC stuff that should be same for V0 and Cluster)
    ProcessClusters();// process calo clusters
    ProcessPhotonCandidates(); // Process this cuts gammas
    if(fDoJetAnalysis)   ProcessJets(); //Process jets
    if(fHistoNGammaCandidates && fHistoNGammaCandidates[iCut]) fHistoNGammaCandidates[iCut]->Fill(fGammaCandidates->GetEntries(),fWeightJetJetMC);
    if(!fDoLightOutput) fHistoNGoodESDTracksVsNGammaCandidates[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks(),fGammaCandidates->GetEntries(),fWeightJetJetMC);
    if(fDoMesonAnalysis){ // Meson Analysis
      if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseMCPSmearing() && fIsMC>0){
        fUnsmearedPx = new Double_t[fGammaCandidates->GetEntries()]; // Store unsmeared Momenta
        fUnsmearedPy = new Double_t[fGammaCandidates->GetEntries()];
        fUnsmearedPz = new Double_t[fGammaCandidates->GetEntries()];
        fUnsmearedE =  new Double_t[fGammaCandidates->GetEntries()];

        for(Int_t gamma=0;gamma<fGammaCandidates->GetEntries();gamma++){ // Smear the AODPhotons in MC
        fUnsmearedPx[gamma] = ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->Px();
        fUnsmearedPy[gamma] = ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->Py();
        fUnsmearedPz[gamma] = ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->Pz();
        fUnsmearedE[gamma] =  ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->E();
        ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->SmearParticle(dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(gamma)));
        }
      }

      CalculatePi0Candidates(); // Combine Gammas from conversion and from calo

      if(fIsMC>0 && fDoHBTHistoOutput && fInputEvent->IsA()==AliAODEvent::Class()){
        ProcessAODMCParticlesForHBT();
      }
      if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->DoBGCalculation()){
        if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->BackgroundHandlerType() == 0){
          CalculateBackground(); // Combinatorial Background
          if(fIsMC>0 && fDoHBTHistoOutput) CalculateHBTBackgroundMC();
          UpdateEventByEventData(); // Store Event for mixed Events
        } else if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->BackgroundHandlerType() == 2){
          CalculateBackgroundSwapp(); // Combinatorial Background
        } else{
          CalculateBackgroundRP(); // Combinatorial Background
          fBGHandlerRP[iCut]->AddEvent(fGammaCandidates,fInputEvent); // Store Event for mixed Events
          fBGClusHandlerRP[iCut]->AddEvent(fClusterCandidates,fInputEvent); // Store Event for mixed Events
        }
      }

      if(fIsMC>0 && fInputEvent->IsA()==AliAODEvent::Class()){
        ProcessConversionPhotonsForMissingTagsAOD(); //Count missing tags
      }else if (fIsMC>0 && fInputEvent->IsA()==AliESDEvent::Class()){
        ProcessConversionPhotonsForMissingTags(); //Count missing tags
      }

      if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseMCPSmearing() && fIsMC>0){
        for(Int_t gamma=0;gamma<fGammaCandidates->GetEntries();gamma++){ // Smear the AODPhotons in MC
        ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->SetPx(fUnsmearedPx[gamma]); // Reset Unsmeared Momenta
        ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->SetPy(fUnsmearedPy[gamma]);
        ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->SetPz(fUnsmearedPz[gamma]);
        ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->SetE(fUnsmearedE[gamma]);
        }
        delete[] fUnsmearedPx; fUnsmearedPx = 0x0;
        delete[] fUnsmearedPy; fUnsmearedPy = 0x0;
        delete[] fUnsmearedPz; fUnsmearedPz = 0x0;
        delete[] fUnsmearedE;fUnsmearedE  = 0x0;
      }

      if(fIsMC>0){
        fVectorRecTruePi0s.clear();
        fVectorRecTrueEtas.clear();
        fVectorDoubleCountTruePi0s.clear();
        fVectorDoubleCountTrueEtas.clear();
        if(fHistoMultipleCountTruePi0 && fHistoMultipleCountTruePi0[iCut]) FillMultipleCountHistoAndClear(fMapMultipleCountTruePi0s,fHistoMultipleCountTruePi0[iCut]);
        if(fHistoMultipleCountTrueEta && fHistoMultipleCountTrueEta[iCut]) FillMultipleCountHistoAndClear(fMapMultipleCountTrueEtas,fHistoMultipleCountTrueEta[iCut]);
      }
    }

    if(fIsMC>0){
      fVectorDoubleCountTrueConvGammas.clear();
      if(!fDoLightOutput) FillMultipleCountHistoAndClear(fMapMultipleCountTrueConvGammas,fHistoMultipleCountTrueConvGamma[iCut]);
      fVectorDoubleCountTrueClusterGammas.clear();
      if(!fDoLightOutput) FillMultipleCountHistoAndClear(fMapMultipleCountTrueClusterGammas,fHistoMultipleCountTrueClusterGamma[iCut]);
    }

    if(fDoInvMassShowerShapeTree) tESDmapIsClusterAcceptedWithoutTrackMatch.clear();

    fGammaCandidates->Clear(); // delete this cuts good gammas
    fClusterCandidates->Clear(); // delete cluster candidates

    if(fIsMC>0 && fDoHBTHistoOutput){
      fTrueGammaCandidatesConv->Clear();
      fTrueGammaCandidatesCalo->Clear();
      fMCGammaCandidates->Clear();
    }
  }

  if(fIsMC>0 && fInputEvent->IsA()==AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled())){
    RelabelAODPhotonCandidates(kFALSE); // Back to ESDMC Label
    fV0Reader->RelabelAODs(kFALSE);
  }

  PostData(1, fOutputContainer);
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::ProcessClusters(){
  Int_t nclus = 0;
  TClonesArray * arrClustersProcess = NULL;
  if(!fCorrTaskSetting.CompareTo("")){
    nclus = fInputEvent->GetNumberOfCaloClusters();
  } else {
    arrClustersProcess = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));
    if(!arrClustersProcess)
      AliFatal(Form("%sClustersBranch was not found in AliAnalysisTaskGammaConvCalo! Check the correction framework settings!",fCorrTaskSetting.Data()));
    nclus = arrClustersProcess->GetEntries();
  }

//   cout << nclus << endl;
  vector<AliAODConversionPhoton*>         vectorCurrentClusters;
  vector<Int_t>                           vectorRejectCluster;
  vector<Double_t>                        vectorPhotonWeight;
  vector<Double_t>                        vectorClusterM02;

  if(nclus == 0)  return;

  // plotting histograms on cell/tower level, only if extendedMatchAndQA > 1
  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FillHistogramsExtendedQA(fInputEvent,fIsMC);

  // match tracks to clusters
  if(fDoPrimaryTrackMatching) ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchTracksToClusters(fInputEvent,fWeightJetJetMC,kFALSE,fMCEvent);

  // vertex
  Double_t vertex[3] = {0,0,0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  // Loop over EMCal clusters
  for(Int_t i = 0; i < nclus; i++){
    AliVCluster* clus = NULL;
    if(fInputEvent->IsA()==AliESDEvent::Class()){
      if(arrClustersProcess)
        clus = new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersProcess->At(i));
      else
        clus = new AliESDCaloCluster(*(AliESDCaloCluster*)fInputEvent->GetCaloCluster(i));
    } else if(fInputEvent->IsA()==AliAODEvent::Class()){
      if(arrClustersProcess)
        clus = new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersProcess->At(i));
      else
        clus = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(i));
    }

    if (!clus) continue;
    if(!((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC,fWeightJetJetMC,i)){
      if(fDoInvMassShowerShapeTree && ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedBeforeTrackMatch() ) tESDmapIsClusterAcceptedWithoutTrackMatch[i] = 1;
      delete clus;
      continue;
    }

    // TLorentzvector with cluster
    TLorentzVector clusterVector;
    clus->GetMomentum(clusterVector,vertex);

    TLorentzVector* tmpvec = new TLorentzVector();
    tmpvec->SetPxPyPzE(clusterVector.Px(),clusterVector.Py(),clusterVector.Pz(),clusterVector.E());

    // convert to AODConversionPhoton
    AliAODConversionPhoton *PhotonCandidate = new AliAODConversionPhoton(tmpvec);
    if(!PhotonCandidate){ delete clus; delete tmpvec; continue;}

    if(fDoInvMassShowerShapeTree ) tESDmapIsClusterAcceptedWithoutTrackMatch[i] = 1;

    // Flag Photon as CaloPhoton
    PhotonCandidate->SetIsCaloPhoton(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType());
    PhotonCandidate->SetCaloClusterRef((Long_t)i);
    // get MC label
    if(fIsMC>0){
      Int_t* mclabelsCluster = clus->GetLabels();
      PhotonCandidate->SetNCaloPhotonMCLabels(clus->GetNLabels());
//       cout << clus->GetNLabels() << endl;
      if (clus->GetNLabels()>0){
        for (Int_t k =0; k<(Int_t)clus->GetNLabels(); k++){
          PhotonCandidate->SetCaloPhotonMCLabel(k,mclabelsCluster[k]);
//           Int_t pdgCode = fMCEvent->Particle(mclabelsCluster[k])->GetPdgCode();
//           cout << "label " << k << "\t" << mclabelsCluster[k] << " pdg code: " << pdgCode << endl;
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
    if(fHistoClusAllHeadersGammaPt && fHistoClusAllHeadersGammaPt[fiCut]) fHistoClusAllHeadersGammaPt[fiCut]->Fill(PhotonCandidate->Pt(),fWeightJetJetMC);
    if (!fIsFromDesiredHeader)
      if(fHistoClusRejectedHeadersGammaPt && fHistoClusRejectedHeadersGammaPt[fiCut]) fHistoClusRejectedHeadersGammaPt[fiCut]->Fill(PhotonCandidate->Pt(),fWeightJetJetMC);
    if (fIsFromDesiredHeader && fIsOverlappingWithOtherHeader)
      if(fHistoClusOverlapHeadersGammaPt && fHistoClusOverlapHeadersGammaPt[fiCut]) fHistoClusOverlapHeadersGammaPt[fiCut]->Fill(PhotonCandidate->Pt(),fWeightJetJetMC);

    if ( (fIsFromDesiredHeader && !fIsOverlappingWithOtherHeader && !fAllowOverlapHeaders) || (fIsFromDesiredHeader && fAllowOverlapHeaders)){
      vectorCurrentClusters.push_back(PhotonCandidate);
      vectorPhotonWeight.push_back(fWeightJetJetMC);
      vectorClusterM02.push_back(clus->GetM02());
    }else{
      delete PhotonCandidate;
    }

    delete clus;
    delete tmpvec;
  }

  //Bool_t rejected = kFALSE;
  // run conversion recovery in addition
  if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetIsConversionRecovery()){
    /*rejected = */((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckForReconstructedConversionPairs(vectorCurrentClusters,vectorRejectCluster);
  }

  //Loop to acquire highest energy clusters, which are not bad from analysis bad map (or trigger bad map)
  Int_t highestClusterE_Iter_AnaBM=-1;
  Int_t highestClusterE_Iter_BothBM=-1;
  Double_t highestClusterE_Value_AnaBM=0;
  Double_t highestClusterE_Value_BothBM=0;
  for (Int_t iter = 0; iter < (Int_t)vectorCurrentClusters.size();iter++){
    if ((vectorCurrentClusters.at(iter)->E())>0.){
        if ((vectorCurrentClusters.at(iter)->E())>highestClusterE_Value_AnaBM){
          highestClusterE_Value_AnaBM=vectorCurrentClusters.at(iter)->E();
          highestClusterE_Iter_AnaBM=iter;
        }

        if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType()==2){
          if (fCaloTriggerMimicHelper[fiCut]){
            if ((fCaloTriggerMimicHelper[fiCut]->IsClusterIDBadMapTrigger(vectorCurrentClusters.at(iter)->GetCaloClusterRef()))>0){
              if ((vectorCurrentClusters.at(iter)->E())>highestClusterE_Value_BothBM){
                if (((((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsSpecialTrigger()==6)&&(fCaloTriggerMimicHelper[fiCut]->IsClusterIDTriggered(vectorCurrentClusters.at(iter)->GetCaloClusterRef())))
                        ||(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsSpecialTrigger()==0)||(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsSpecialTrigger()==1)){
                  highestClusterE_Value_BothBM=vectorCurrentClusters.at(iter)->E();
                  highestClusterE_Iter_BothBM=iter;
                }
              }
            }
          }
        }
    }
  }
  //Fill only highest cluster to Histograms
  //MB and PHI7
  if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType()==2){
    if (fCaloTriggerMimicHelper[fiCut]){
      if (highestClusterE_Iter_BothBM!=-1){
        fHistoClusGammaE_BothBM_highestE[fiCut]->Fill(vectorCurrentClusters.at(highestClusterE_Iter_BothBM)->E(), vectorPhotonWeight.at(highestClusterE_Iter_BothBM));
      }
    }
    //Only MB
    //if ( (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsSpecialTrigger()==1)||(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsSpecialTrigger()==0) ){
      if (highestClusterE_Iter_AnaBM!=-1){
        fHistoClusGammaE_AnaBM_highestE[fiCut]->Fill(vectorCurrentClusters.at(highestClusterE_Iter_AnaBM)->E(), vectorPhotonWeight.at(highestClusterE_Iter_AnaBM));
      }
    //}
  }

  for (Int_t iter = 0; iter < (Int_t)vectorCurrentClusters.size();iter++){

    if (!((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckVectorForIndexAndAdd(vectorRejectCluster, iter,kFALSE)){
      if(fHistoClusGammaPt && fHistoClusGammaPt[fiCut]) fHistoClusGammaPt[fiCut]->Fill(vectorCurrentClusters.at(iter)->Pt(), vectorPhotonWeight.at(iter));
      if(fHistoClusGammaE && fHistoClusGammaE[fiCut]) fHistoClusGammaE[fiCut]->Fill(vectorCurrentClusters.at(iter)->E(), vectorPhotonWeight.at(iter));

      if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType()==2){
        //--------------------------------------------------
        //Only MB
        //if ( (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsSpecialTrigger()==1)||(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsSpecialTrigger()==0) ){
          if (fCaloTriggerMimicHelper[fiCut]){
            if ((vectorCurrentClusters.at(iter)->E())>0){//Analysis bad map protection; Energy of bad clusters set to 0
              if ((fCaloTriggerMimicHelper[fiCut]->IsClusterIDBadMapTrigger(vectorCurrentClusters.at(iter)->GetCaloClusterRef()))>0){//Good cluster by bad trigger map decision
                fHistoClusGammaE_BothBM[fiCut]->Fill(vectorCurrentClusters.at(iter)->E(), vectorPhotonWeight.at(iter));
              }
            }
          }
        //}
        //--------------------------------------------------
        //Only PHI7
        if ( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsSpecialTrigger()==6 ){
          if (fCaloTriggerMimicHelper[fiCut]){
            if (fCaloTriggerMimicHelper[fiCut]->IsClusterIDTriggered(vectorCurrentClusters.at(iter)->GetCaloClusterRef())){
              fHistoClusGammaE_onlyTriggered[fiCut]->Fill(vectorCurrentClusters.at(iter)->E(), vectorPhotonWeight.at(iter));
            }
          }
        }
      }

      if(fIsMC> 0){
        if(fInputEvent->IsA()==AliESDEvent::Class()){
          ProcessTrueClusterCandidates(vectorCurrentClusters.at(iter),vectorClusterM02.at(iter));
        } else {
          ProcessTrueClusterCandidatesAOD(vectorCurrentClusters.at(iter),vectorClusterM02.at(iter));
        }
      }
      fClusterCandidates->Add(vectorCurrentClusters.at(iter));
    }
  }
  vectorClusterM02.clear();
  vectorPhotonWeight.clear();
  vectorRejectCluster.clear();


}

//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::ProcessTrueClusterCandidates(AliAODConversionPhoton *TruePhotonCandidate, Float_t clusM02)
{

  TParticle *Photon = NULL;
  if (TruePhotonCandidate->GetIsCaloPhoton() == 0) AliFatal("CaloPhotonFlag has not been set task will abort");
  if (TruePhotonCandidate->GetCaloPhotonMCLabel(0) < 0) return;
  if(!fDoLightOutput) fHistoTrueNLabelsInClusPt[fiCut]->Fill(TruePhotonCandidate->GetNCaloPhotonMCLabels(),TruePhotonCandidate->Pt(),fWeightJetJetMC);

  if (TruePhotonCandidate->GetNCaloPhotonMCLabels()>0)Photon = fMCEvent->Particle(TruePhotonCandidate->GetCaloPhotonMCLabel(0));
    else return;

  if(Photon == NULL){
  //    cout << "no photon" << endl;
    return;
  }

  TruePhotonCandidate->SetCaloPhotonMCFlags(fMCEvent, fEnableSortForClusMC);

  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX           = primVtxMC->GetX();
  Double_t mcProdVtxY           = primVtxMC->GetY();
  Double_t mcProdVtxZ           = primVtxMC->GetZ();
  Bool_t isPrimary              = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, TruePhotonCandidate->GetCaloPhotonMCLabel(0), mcProdVtxX, mcProdVtxY, mcProdVtxZ);

  // to get primary distrinction right put mother of conversion electron as particle to check
  if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion())
    isPrimary              = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, Photon->GetMother(0), mcProdVtxX, mcProdVtxY, mcProdVtxZ);


  // Fill histograms for inclusive gamma corrections
  // --> a) all clusters with leading real or converted photons
  if (TruePhotonCandidate->IsLargestComponentPhoton() || (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()) ){
    if(fHistoTrueClusGammaPt && fHistoTrueClusGammaPt[fiCut]) fHistoTrueClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC);
    if(!fDoLightOutput){
      // how many of those clusters are from converted photons
      if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
        fHistoTrueClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
      }
      // --> b) which of these are primary
      if(isPrimary){
        fHistoTruePrimaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC);
        fHistoTruePrimaryClusGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt(),fWeightJetJetMC);
        // --> c) which are from conversions? Just additonal information
        if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion() && (Photon->GetMother(0)>-1)){
          fHistoTruePrimaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
          if(Photon->GetMother(0) > -1){
            fHistoTruePrimaryClusConvGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),((TParticle*)fMCEvent->Particle(Photon->GetMother(0)))->Pt(),fWeightJetJetMC);
          }
        }
      // --> d) how do the secondaries look like
      }else {
        Int_t secondaryClass    = -1;
        if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion())
          secondaryClass        = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->SecondaryClassificationPhoton( Photon, fMCEvent, kTRUE);
        else
          secondaryClass        = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->SecondaryClassificationPhoton( Photon, fMCEvent, kFALSE);
        // all secondaries
        if (secondaryClass > 0 ){
          fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC);
        }
        // secondaries from K0s
        if (secondaryClass == 2)
          fHistoTrueSecondaryClusGammaFromK0sPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC);
        // secondaries from Lambda
        else if (secondaryClass == 3)
          fHistoTrueSecondaryClusGammaFromLambdaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC);
        else if (secondaryClass == 5)
          fHistoTrueSecondaryClusGammaFromK0lPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC);
      }
    }

  // How many clusters are from electrons besides conversions
  }else if (TruePhotonCandidate->IsLargestComponentElectron()) {
    if(!fDoLightOutput) fHistoTrueClusElectronPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC);
  // How many clusters don't have an electromagnetic particle as leading particle
  }else {
    if(!fDoLightOutput) fHistoTrueClusEMNonLeadingPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC);
  }

  // Filling histogram for true photons: E vs M02 for paper
  if(!fDoLightOutput) {
    if (  TruePhotonCandidate->IsLargestComponentPhoton() && !TruePhotonCandidate->IsPhotonWithElecMother() &&
          !TruePhotonCandidate->IsMerged() && !TruePhotonCandidate->IsMergedPartConv() && !TruePhotonCandidate->IsDalitzMerged() )
      fHistoTrueClusGammaEM02[fiCut]->Fill(TruePhotonCandidate->E(),clusM02, fWeightJetJetMC);
  }

  // Some additional QA
  if (fDoClusterQA > 0 && fDoLightOutput!=2){
    // how many of the converted photons are fully contained in the cluster
    if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion() && TruePhotonCandidate->IsConversionFullyContained())
      fHistoTrueClusConvGammaFullyPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
    // how often do we have merged pi0/eta...
    if ( (TruePhotonCandidate->IsMerged() || TruePhotonCandidate->IsDalitzMerged()) || TruePhotonCandidate->IsMergedPartConv() )
      fHistoTrueClusMergedGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
    // how often do we have a merged cluster with at least one conversion
    if (TruePhotonCandidate->IsMergedPartConv())
      fHistoTrueClusMergedPartConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
    // how often do we reconstruct Dalitz electrons
    if (TruePhotonCandidate->IsDalitz())
      fHistoTrueClusDalitzPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
    // how often merge Dalitz decays
    if (TruePhotonCandidate->IsDalitzMerged())
      fHistoTrueClusDalitzMergedPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
    // how often do we see Bremstrahlung
    if (TruePhotonCandidate->IsPhotonWithElecMother())
      fHistoTrueClusPhotonFromElecMotherPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
    // how often do we see a shower in the clusters
    if (TruePhotonCandidate->IsShower())
      fHistoTrueClusShowerPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
    // how often is the EM a subleading contributor
    if (TruePhotonCandidate->IsSubLeadingEM())
      fHistoTrueClusSubLeadingPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
    // how many mother particles point to the cluster
    fHistoTrueClusNMothers[fiCut]->Fill(TruePhotonCandidate->GetNCaloPhotonMotherMCLabels(), fWeightJetJetMC);
  }

  // Check if we are double counting photons
  Int_t motherLab = Photon->GetMother(0);
  if (motherLab > -1){
    if (TruePhotonCandidate->IsLargestComponentPhoton()){
      if (CheckVectorForDoubleCount(fVectorDoubleCountTrueClusterGammas,motherLab)){
        if(fHistoDoubleCountTrueClusterGammaPt && fHistoDoubleCountTrueClusterGammaPt[fiCut]) fHistoDoubleCountTrueClusterGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),(Double_t)0,fWeightJetJetMC);
        FillMultipleCountMap(fMapMultipleCountTrueClusterGammas,motherLab);
      }
    }
    if ( TMath::Abs(fMCEvent->Particle(motherLab)->GetPdgCode()) == 111 &&
      TruePhotonCandidate->IsLargestComponentPhoton() &&
      TruePhotonCandidate->IsMerged() &&
      !TruePhotonCandidate->IsDalitzMerged() &&
      !TruePhotonCandidate->IsMergedPartConv()
    ){
      if(!fDoLightOutput) fHistoTrueClusPi0EM02[fiCut]->Fill(TruePhotonCandidate->E(),clusM02,fWeightJetJetMC);
    }
    Int_t grandMotherLab = fMCEvent->Particle(motherLab)->GetMother(0);
    if (grandMotherLab > -1){
      if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
        if (CheckVectorForDoubleCount(fVectorDoubleCountTrueClusterGammas,grandMotherLab)){
          if(fHistoDoubleCountTrueClusterGammaPt && fHistoDoubleCountTrueClusterGammaPt[fiCut]) fHistoDoubleCountTrueClusterGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),(Double_t)1,fWeightJetJetMC);
          FillMultipleCountMap(fMapMultipleCountTrueClusterGammas,grandMotherLab);
        }
      }
    }
  }

  return;
}


//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::ProcessTrueClusterCandidatesAOD(AliAODConversionPhoton *TruePhotonCandidate, Float_t clusM02)
{
  AliAODMCParticle *Photon = NULL;
  if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (fAODMCTrackArray){
    if (TruePhotonCandidate->GetIsCaloPhoton() == 0) AliFatal("CaloPhotonFlag has not been set task will abort");
    if (TruePhotonCandidate->GetNCaloPhotonMCLabels()>0) Photon = (AliAODMCParticle*) fAODMCTrackArray->At(TruePhotonCandidate->GetCaloPhotonMCLabel(0));
      else return;
  }else {
    AliInfo("AODMCTrackArray could not be loaded");
    return;
  }

  if(Photon == NULL){
  //  cout << "no photon" << endl;
    return;
  }

  TruePhotonCandidate->SetCaloPhotonMCFlagsAOD(fAODMCTrackArray, fEnableSortForClusMC);
  if(!fDoLightOutput) fHistoTrueNLabelsInClusPt[fiCut]->Fill(TruePhotonCandidate->GetNCaloPhotonMCLabels(),TruePhotonCandidate->Pt(),fWeightJetJetMC);

  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();
  Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, Photon, mcProdVtxX, mcProdVtxY, mcProdVtxZ);

  // to get primary distrinction right put mother of conversion electron as particle to check
  if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
    if (Photon->GetMother()> -1){
      AliAODMCParticle *Mother  = (AliAODMCParticle*) fAODMCTrackArray->At(Photon->GetMother());
      isPrimary                 = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD( fInputEvent, Mother, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
    }
  }

  // True Photon
  if (TruePhotonCandidate->IsLargestComponentPhoton() || (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()) ){
    if(fHistoTrueClusGammaPt && fHistoTrueClusGammaPt[fiCut]) fHistoTrueClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC);
    if(!fDoLightOutput){
      if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
        fHistoTrueClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
      }
      if(isPrimary){
        fHistoTruePrimaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC);
        fHistoTruePrimaryClusGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt(),fWeightJetJetMC);
        if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
          fHistoTruePrimaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
          if(Photon->GetMother()>-1){
            AliAODMCParticle *Mother = (AliAODMCParticle*) fAODMCTrackArray->At(Photon->GetMother());
            fHistoTruePrimaryClusConvGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Mother->Pt(),fWeightJetJetMC);
          }
        }
      }else {
        Int_t secondaryClass    = -1;
        if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion())
          secondaryClass        = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->SecondaryClassificationPhotonAOD( Photon, fAODMCTrackArray, kTRUE);
        else
          secondaryClass        = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->SecondaryClassificationPhotonAOD( Photon, fAODMCTrackArray, kFALSE);
        // all secondaries
        if (secondaryClass > 0 ){
          fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC);
        }
        // secondaries from K0s
        if (secondaryClass == 2)
          fHistoTrueSecondaryClusGammaFromK0sPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC);
        // secondaries from Lambda
        else if (secondaryClass == 3)
          fHistoTrueSecondaryClusGammaFromLambdaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC);
        else if (secondaryClass == 5)
          fHistoTrueSecondaryClusGammaFromK0lPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC);
      }
    }

  }else if (TruePhotonCandidate->IsLargestComponentElectron()) {
    if(!fDoLightOutput) fHistoTrueClusElectronPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC);
  }else {
    if(!fDoLightOutput) fHistoTrueClusEMNonLeadingPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC);
  }

  if(!fDoLightOutput) {
    if ( TruePhotonCandidate->IsLargestComponentPhoton() && !TruePhotonCandidate->IsPhotonWithElecMother() &&
          !TruePhotonCandidate->IsMerged() && !TruePhotonCandidate->IsMergedPartConv() && !TruePhotonCandidate->IsDalitzMerged() )
      fHistoTrueClusGammaEM02[fiCut]->Fill(TruePhotonCandidate->E(),clusM02, fWeightJetJetMC);
  }
  if (fDoClusterQA > 0 && fDoLightOutput!=2){
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
    fHistoTrueClusNMothers[fiCut]->Fill(TruePhotonCandidate->GetNCaloPhotonMotherMCLabels(), fWeightJetJetMC);

  }
  Int_t motherLab = Photon->GetMother();
  if (motherLab > -1){
    if (TruePhotonCandidate->IsLargestComponentPhoton()){
      if (CheckVectorForDoubleCount(fVectorDoubleCountTrueClusterGammas,motherLab)){
        if(fHistoDoubleCountTrueClusterGammaPt && fHistoDoubleCountTrueClusterGammaPt[fiCut]) fHistoDoubleCountTrueClusterGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),(Double_t)0,fWeightJetJetMC);
        FillMultipleCountMap(fMapMultipleCountTrueClusterGammas,motherLab);
      }
    }
    if(!fDoLightOutput) {
      if ( TMath::Abs(((AliAODMCParticle*) fAODMCTrackArray->At(motherLab))->GetPdgCode()) == 111 && TruePhotonCandidate->IsLargestComponentPhoton() && TruePhotonCandidate->IsMerged()  && !TruePhotonCandidate->IsDalitzMerged() && !TruePhotonCandidate->IsMergedPartConv())
        fHistoTrueClusPi0EM02[fiCut]->Fill(TruePhotonCandidate->E(),clusM02, fWeightJetJetMC);
    }
    Int_t grandMotherLab = ((AliAODMCParticle*) fAODMCTrackArray->At(motherLab))->GetMother();
    if (grandMotherLab > -1){
      if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
        if (CheckVectorForDoubleCount(fVectorDoubleCountTrueClusterGammas,grandMotherLab)){
          if(fHistoDoubleCountTrueClusterGammaPt && fHistoDoubleCountTrueClusterGammaPt[fiCut]) fHistoDoubleCountTrueClusterGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),(Double_t)1,fWeightJetJetMC);
          FillMultipleCountMap(fMapMultipleCountTrueClusterGammas,grandMotherLab);
        }
      }
    }
  }

  return;
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::ProcessPhotonCandidates()
{
  Double_t magField = fInputEvent->GetMagneticField();
  Int_t nV0 = 0;
  TList *GammaCandidatesStepOne = new TList();
  TList *GammaCandidatesStepTwo = new TList();
  // Loop over Photon Candidates allocated by ReaderV1
  for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){
    AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(i);
    if(!PhotonCandidate) continue;
    fIsFromDesiredHeader = kTRUE;

    Double_t weightMatBudgetGamma = 1.;
    if (fDoMaterialBudgetWeightingOfGammasForTrueMesons && ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetMaterialBudgetWeightsInitialized()) {
      weightMatBudgetGamma = ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetMaterialBudgetCorrectingWeightForTrueGamma(PhotonCandidate, magField);
    }


    if(fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
      Int_t isPosFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCEvent, fInputEvent);
      if(isPosFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
      Int_t isNegFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCEvent, fInputEvent);
      if(isNegFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
      if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromDesiredHeader = kFALSE;
    }

    if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelected(PhotonCandidate,fInputEvent)) continue;
    if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(PhotonCandidate->GetPhotonPhi(),fEventPlaneAngle)) continue;
    if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseElecSharingCut() &&
    !((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){
      fGammaCandidates->Add(PhotonCandidate); // if no second loop is required add to events good gammas

      if(fIsFromDesiredHeader){
        if(!fDoLightOutput) fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
      }
      if(fIsMC>0){
        if(fInputEvent->IsA()==AliESDEvent::Class())
        ProcessTruePhotonCandidates(PhotonCandidate);
        if(fInputEvent->IsA()==AliAODEvent::Class())
        ProcessTruePhotonCandidatesAOD(PhotonCandidate);
      }
      if (fIsFromDesiredHeader && fDoPhotonQA == 2){
        if (fIsHeavyIon == 1 && PhotonCandidate->Pt() > 0.399 && PhotonCandidate->Pt() < 30.){
          fPtGamma        = (UShort_t)(PhotonCandidate->Pt()*1000);
          fDCAzPhoton     = (Short_t)(PhotonCandidate->GetDCAzToPrimVtx()*1000);
          if (fIsMC > 0){
            fCharCatPhoton  = fCharPhotonMCInfo*10+PhotonCandidate->GetPhotonQuality();
          } else {
            fCharCatPhoton  = PhotonCandidate->GetPhotonQuality();
          }
          fTreeConvGammaPtDcazCat[fiCut]->Fill();
        }else if ( PhotonCandidate->Pt() > 0.299 && PhotonCandidate->Pt() < 30.){
          fPtGamma        = (UShort_t)(PhotonCandidate->Pt()*1000);
          fDCAzPhoton     = (Short_t)(PhotonCandidate->GetDCAzToPrimVtx()*1000);
          if (fIsMC > 0){
            fCharCatPhoton  = fCharPhotonMCInfo*10+PhotonCandidate->GetPhotonQuality();
          } else {
            fCharCatPhoton  = PhotonCandidate->GetPhotonQuality();
          }
          fTreeConvGammaPtDcazCat[fiCut]->Fill();
        }
      }
    }else if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseElecSharingCut()){ // if Shared Electron cut is enabled, Fill array, add to step one
      ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->FillElectonLabelArray(PhotonCandidate,nV0);
      nV0++;
      GammaCandidatesStepOne->Add(PhotonCandidate);
    }else if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseElecSharingCut() &&
        ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){ // shared electron is disabled, step one not needed -> step two
      GammaCandidatesStepTwo->Add(PhotonCandidate);
    }
  }
  if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseElecSharingCut()){
    for(Int_t i = 0;i<GammaCandidatesStepOne->GetEntries();i++){
      AliAODConversionPhoton *PhotonCandidate= (AliAODConversionPhoton*) GammaCandidatesStepOne->At(i);
      if(!PhotonCandidate) continue;
      fIsFromDesiredHeader = kTRUE;

      Double_t weightMatBudgetGamma = 1.;
      if (fDoMaterialBudgetWeightingOfGammasForTrueMesons && ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetMaterialBudgetWeightsInitialized()) {
	weightMatBudgetGamma = ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetMaterialBudgetCorrectingWeightForTrueGamma(PhotonCandidate, magField);
      }

      if(fMCEvent && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        Int_t isPosFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCEvent, fInputEvent);
        Int_t isNegFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCEvent, fInputEvent);
        if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromDesiredHeader = kFALSE;
      }
      if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->RejectSharedElectronV0s(PhotonCandidate,i,GammaCandidatesStepOne->GetEntries())) continue;
      if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){ // To Colse v0s cut diabled, step two not needed
        fGammaCandidates->Add(PhotonCandidate);
        if(fIsFromDesiredHeader){
          if(!fDoLightOutput) fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
        }
        if(fIsMC>0){
          if(fInputEvent->IsA()==AliESDEvent::Class())
            ProcessTruePhotonCandidates(PhotonCandidate);
          if(fInputEvent->IsA()==AliAODEvent::Class())
            ProcessTruePhotonCandidatesAOD(PhotonCandidate);
        }
        if (fIsFromDesiredHeader && fDoPhotonQA == 2){
          if (fIsHeavyIon ==1 && PhotonCandidate->Pt() > 0.399 && PhotonCandidate->Pt() < 30.){
            fPtGamma        = (UShort_t)(PhotonCandidate->Pt()*1000);
            fDCAzPhoton     = (Short_t)(PhotonCandidate->GetDCAzToPrimVtx()*1000);
            if (fIsMC > 0){
              fCharCatPhoton  = fCharPhotonMCInfo*10+PhotonCandidate->GetPhotonQuality();
            } else {
              fCharCatPhoton  = PhotonCandidate->GetPhotonQuality();
            }
            fTreeConvGammaPtDcazCat[fiCut]->Fill();
          }else if ( PhotonCandidate->Pt() > 0.299 && PhotonCandidate->Pt() < 30.){
            fPtGamma        = (UShort_t)(PhotonCandidate->Pt()*1000);
            fDCAzPhoton     = (Short_t)(PhotonCandidate->GetDCAzToPrimVtx()*1000);
            if (fIsMC > 0){
              fCharCatPhoton  = fCharPhotonMCInfo*10+PhotonCandidate->GetPhotonQuality();
            } else {
              fCharCatPhoton  = PhotonCandidate->GetPhotonQuality();
            }
            fTreeConvGammaPtDcazCat[fiCut]->Fill();
          }
        }
      } else GammaCandidatesStepTwo->Add(PhotonCandidate); // Close v0s cut enabled -> add to list two
    }
  }
  if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){
    for(Int_t i = 0;i<GammaCandidatesStepTwo->GetEntries();i++){
      AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) GammaCandidatesStepTwo->At(i);
      if(!PhotonCandidate) continue;
      fIsFromDesiredHeader = kTRUE;

      Double_t weightMatBudgetGamma = 1.;
      if (fDoMaterialBudgetWeightingOfGammasForTrueMesons && ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetMaterialBudgetWeightsInitialized()) {
	weightMatBudgetGamma = ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetMaterialBudgetCorrectingWeightForTrueGamma(PhotonCandidate, magField);
      }

      if(fMCEvent && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        Int_t isPosFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCEvent, fInputEvent);
        Int_t isNegFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCEvent, fInputEvent);
        if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromDesiredHeader = kFALSE;
      }
      if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->RejectToCloseV0s(PhotonCandidate,GammaCandidatesStepTwo,i)) continue;
      fGammaCandidates->Add(PhotonCandidate); // Add gamma to current cut TList
      if(fIsFromDesiredHeader){
        if(!fDoLightOutput) fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
      }
      if(fIsMC>0){
        if(fInputEvent->IsA()==AliESDEvent::Class())
          ProcessTruePhotonCandidates(PhotonCandidate);
        if(fInputEvent->IsA()==AliAODEvent::Class())
          ProcessTruePhotonCandidatesAOD(PhotonCandidate);
      }
      if (fIsFromDesiredHeader && fDoPhotonQA == 2){
        if (fIsHeavyIon == 1 && PhotonCandidate->Pt() > 0.399 && PhotonCandidate->Pt() < 65.){
          fPtGamma        = (UShort_t)(PhotonCandidate->Pt()*1000);
          fDCAzPhoton     = (Short_t)(PhotonCandidate->GetDCAzToPrimVtx()*1000);
          if (fIsMC > 0){
            fCharCatPhoton  = fCharPhotonMCInfo*10+PhotonCandidate->GetPhotonQuality();
          } else {
            fCharCatPhoton  = PhotonCandidate->GetPhotonQuality();
          }
          fTreeConvGammaPtDcazCat[fiCut]->Fill();
        }else if ( PhotonCandidate->Pt() > 0.299 && PhotonCandidate->Pt() < 65.){
          fPtGamma        = (UShort_t)(PhotonCandidate->Pt()*1000);
          fDCAzPhoton     = (Short_t)(PhotonCandidate->GetDCAzToPrimVtx()*1000);
          if (fIsMC > 0){
            fCharCatPhoton  = fCharPhotonMCInfo*10+PhotonCandidate->GetPhotonQuality();
          } else {
            fCharCatPhoton  = PhotonCandidate->GetPhotonQuality();
          }
          fTreeConvGammaPtDcazCat[fiCut]->Fill();
        }
      }
    }
  }

  delete GammaCandidatesStepOne;
  GammaCandidatesStepOne = 0x0;
  delete GammaCandidatesStepTwo;
  GammaCandidatesStepTwo = 0x0;

}

//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::ProcessJets()
{
  fHistoNJets[fiCut]->Fill(fConvJetReader->GetNJets());
  if(fConvJetReader->GetNJets()>0){
    fVectorJetPt  = fConvJetReader->GetVectorJetPt();
    fVectorJetPx  = fConvJetReader->GetVectorJetPx();
    fVectorJetPy  = fConvJetReader->GetVectorJetPy();
    fVectorJetPz  = fConvJetReader->GetVectorJetPz();
    fVectorJetEta = fConvJetReader->GetVectorJetEta();
    fVectorJetPhi = fConvJetReader->GetVectorJetPhi();
    fVectorJetArea = fConvJetReader->GetVectorJetArea();
    if(fIsMC > 0 && fConvJetReader->GetTrueNJets()>0){
      fTrueVectorJetPt = fConvJetReader->GetTrueVectorJetPt();
      fTrueVectorJetEta = fConvJetReader->GetTrueVectorJetEta();
      fTrueVectorJetPhi = fConvJetReader->GetTrueVectorJetPhi();
    }
    if(fVectorJetPt.size() == fConvJetReader->GetNJets() && fVectorJetEta.size() == fConvJetReader->GetNJets() && fVectorJetPhi.size() == fConvJetReader->GetNJets() && fVectorJetArea.size() == fConvJetReader->GetNJets()){
      for(Int_t i=0; i<fConvJetReader->GetNJets(); i++){
        fHistoPtJet[fiCut]->Fill(fVectorJetPt.at(i));
        fHistoJetEta[fiCut]->Fill(fVectorJetEta.at(i));
        fHistoJetPhi[fiCut]->Fill(fVectorJetPhi.at(i));
        fHistoJetArea[fiCut]->Fill(fVectorJetArea.at(i));
        if(!fDoJetQA) fHistoEventwJets[fiCut]->Fill(0);
        if(fIsMC > 0 && fConvJetReader->GetNJets()>0 && fConvJetReader->GetTrueNJets()>0){
          Double_t min = 100;
          Int_t match = 0;
          for(Int_t j = 0; j<fConvJetReader->GetTrueNJets(); j++){
            Double_t R_jetjet;
            Double_t DeltaEta = fVectorJetEta.at(i)-fTrueVectorJetEta.at(j);
            Double_t DeltaPhi = abs(fVectorJetPhi.at(i)-fTrueVectorJetPhi.at(j));
            if(DeltaPhi > M_PI) {
              DeltaPhi = 2*M_PI - DeltaPhi;
            }
            R_jetjet = TMath::Sqrt(pow((DeltaEta),2)+pow((DeltaPhi),2));
            if(R_jetjet < min){
              min = R_jetjet;
              match = j;
            }
          }
          if(!fDoLightOutput) fHistoTruevsRecJetPt[fiCut]->Fill(fVectorJetPt.at(i), fTrueVectorJetPt.at(match));
          if(fDoJetQA){
            if(fVectorJetPt.at(i) >= 10) fHistoEventwJets[fiCut]->Fill(0);
            if(fVectorJetPt.at(i) < 10 && fTrueVectorJetPt.at(match) >= 10) fHistoEventwJets[fiCut]->Fill(1);
            if(fVectorJetPt.at(i) >= 10 && fTrueVectorJetPt.at(match) < 10) fHistoEventwJets[fiCut]->Fill(2);
          }
        }
      }

      if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoJetMixing()){
        Double_t MaxPt = 0.;
        for(Int_t k = 0; k < fConvJetReader->GetNJets(); k++){
          if(fVectorJetEta.at(k) > (-0.6687 -0.4) && fVectorJetEta.at(k) < (0.66465 + 0.4)){
            if(fVectorJetPhi.at(k) > (1.39626 -0.4) && fVectorJetPhi.at(k) < (3.15 + 0.4)){ //INSIDE EMCAL
              fJetNearEMCal = kTRUE;
              if(fVectorJetPhi.at(k) > MaxPt){
                MaxPt = fVectorJetPt.at(k);
                fMaxPtNearEMCalPlace = k;
              }
            }else fJetNearEMCal = kFALSE;
          }else fJetNearEMCal = kFALSE;
        }
      }
    }
    fVectorJetPt.clear();
    fVectorJetPx.clear();
    fVectorJetPy.clear();
    fVectorJetPz.clear();
    fVectorJetEta.clear();
    fVectorJetPhi.clear();
    fVectorJetArea.clear();
    if(fIsMC > 0 && fConvJetReader->GetTrueNJets()>0){
      fTrueVectorJetPt.clear();
      fTrueVectorJetEta.clear();
      fTrueVectorJetPhi.clear();
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::ProcessTruePhotonCandidatesAOD(AliAODConversionPhoton *TruePhotonCandidate)
{
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  Double_t magField = fInputEvent->GetMagneticField();
  // Double_t magFieldFlip = 1.0;
  // if( magField  < 0.0 ){
  //   magFieldFlip =  1.0;
  // }
  // else {
  //   magFieldFlip =  -1.0;
  // }

  if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (fAODMCTrackArray == NULL) return;
  AliAODMCParticle *posDaughter = (AliAODMCParticle*) fAODMCTrackArray->At(TruePhotonCandidate->GetMCLabelPositive());
  AliAODMCParticle *negDaughter = (AliAODMCParticle*) fAODMCTrackArray->At(TruePhotonCandidate->GetMCLabelNegative());
  fCharPhotonMCInfo = 0;

  if(posDaughter == NULL || negDaughter == NULL) return; // One particle does not exist
  Int_t pdgCode[2] = {TMath::Abs(posDaughter->GetPdgCode()),TMath::Abs(negDaughter->GetPdgCode())};

  if(posDaughter->GetMother() != negDaughter->GetMother()){
    FillPhotonCombinatorialBackgroundHist(TruePhotonCandidate, pdgCode);
    fCharPhotonMCInfo = 1;
    return;
  }
  else if(posDaughter->GetMother() == -1){
    FillPhotonCombinatorialBackgroundHist(TruePhotonCandidate, pdgCode);
    fCharPhotonMCInfo = 1;
    return;
  }

  if(pdgCode[0]!=11 || pdgCode[1]!=11){
    fCharPhotonMCInfo = 1;
    return; //One Particle is not a electron
  }

  if(posDaughter->GetPdgCode()==negDaughter->GetPdgCode()){
    fCharPhotonMCInfo = 1;
    return; // Same Charge
  }

  AliAODMCParticle *Photon = (AliAODMCParticle*) fAODMCTrackArray->At(posDaughter->GetMother());
  if(Photon->GetPdgCode() != 22){
    fCharPhotonMCInfo = 1;
    return; // Mother is no Photon
  }

  if(((posDaughter->GetMCProcessCode())) != 5 || ((negDaughter->GetMCProcessCode())) != 5){
    fCharPhotonMCInfo = 1;
    return;// check if the daughters come from a conversion
  }
  // STILL A BUG IN ALIROOT >>8 HAS TPO BE REMOVED AFTER FIX



  // True Photon

  Double_t weightMatBudgetGamma = 1.;
  if (fDoMaterialBudgetWeightingOfGammasForTrueMesons && ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetMaterialBudgetWeightsInitialized()) {
    weightMatBudgetGamma = ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetMaterialBudgetCorrectingWeightForTrueGamma(TruePhotonCandidate, magField);
  }


  if(fIsFromDesiredHeader){
    if(!fDoLightOutput) fHistoTrueConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
    if (CheckVectorForDoubleCount(fVectorDoubleCountTrueConvGammas,posDaughter->GetMother())){
      if(!fDoLightOutput) fHistoDoubleCountTrueConvGammaRPt[fiCut]->Fill(TruePhotonCandidate->GetConversionRadius(),TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
      FillMultipleCountMap(fMapMultipleCountTrueConvGammas,posDaughter->GetMother());
    }
  }

  Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, Photon, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
  if(isPrimary){
    // Count just primary MC Gammas as true --> For Ratio esdtruegamma / mcconvgamma
    if(fIsFromDesiredHeader){
      fCharPhotonMCInfo = 6;
      if(!fDoLightOutput){
        fHistoTruePrimaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
        fHistoTruePrimaryConvGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt(),fWeightJetJetMC*weightMatBudgetGamma); // Allways Filled
      }
    }
    // (Not Filled for i6, Extra Signal Gamma (parambox) are secondary)
  }else {
    if(fIsFromDesiredHeader){

      fCharPhotonMCInfo = 2;
      // check for secondaries from K0s
      if(((AliAODMCParticle*)fAODMCTrackArray->At(Photon->GetMother()))->GetMother() > -1 &&
          ((AliAODMCParticle*)fAODMCTrackArray->At(((AliAODMCParticle*)fAODMCTrackArray->At(Photon->GetMother()))->GetMother()))->GetPdgCode() == 310){
        fCharPhotonMCInfo = 4;
        if(!fDoLightOutput){
          fHistoTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),0.,fWeightJetJetMC*weightMatBudgetGamma);
          fHistoTrueSecondaryConvGammaMCPt[fiCut]->Fill(Photon->Pt(),0.,fWeightJetJetMC*weightMatBudgetGamma);
          fHistoTrueSecondaryConvGammaFromXFromK0sMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
        }
      // check for secondaries from K0l
      } else if(((AliAODMCParticle*)fAODMCTrackArray->At(Photon->GetMother()))->GetMother() > -1 &&
                  ((AliAODMCParticle*)fAODMCTrackArray->At(((AliAODMCParticle*)fAODMCTrackArray->At(Photon->GetMother()))->GetMother()))->GetPdgCode() == 130) {
        fCharPhotonMCInfo = 7;
        if(!fDoLightOutput){
          fHistoTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),1.,fWeightJetJetMC*weightMatBudgetGamma);
          fHistoTrueSecondaryConvGammaMCPt[fiCut]->Fill(Photon->Pt(),1.,fWeightJetJetMC*weightMatBudgetGamma);
          fHistoTrueSecondaryConvGammaFromXFromK0lMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
        }
      // check for secondaries from Lambda
      } else if(((AliAODMCParticle*)fAODMCTrackArray->At(Photon->GetMother()))->GetMother() > -1 &&
        ((AliAODMCParticle*)fAODMCTrackArray->At(((AliAODMCParticle*)fAODMCTrackArray->At(Photon->GetMother()))->GetMother()))->GetPdgCode() == 3122){
        fCharPhotonMCInfo = 5;
        if(!fDoLightOutput){
          fHistoTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),2.,fWeightJetJetMC*weightMatBudgetGamma);
          fHistoTrueSecondaryConvGammaMCPt[fiCut]->Fill(Photon->Pt(),0.,fWeightJetJetMC*weightMatBudgetGamma);
          fHistoTrueSecondaryConvGammaFromXFromLambdaMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
        }
      } else if (((AliAODMCParticle*)fAODMCTrackArray->At(Photon->GetMother()))->GetMother() > -1 &&
                  ((AliAODMCParticle*)fAODMCTrackArray->At(((AliAODMCParticle*)fAODMCTrackArray->At(Photon->GetMother()))->GetMother()))->GetPdgCode() == 221) {
        fCharPhotonMCInfo = 3;
        if(!fDoLightOutput){
          fHistoTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3.,fWeightJetJetMC*weightMatBudgetGamma);
          fHistoTrueSecondaryConvGammaMCPt[fiCut]->Fill(Photon->Pt(),3.,fWeightJetJetMC*weightMatBudgetGamma);
        }
      } else {
        if(!fDoLightOutput){
          fHistoTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3.,fWeightJetJetMC*weightMatBudgetGamma);
          fHistoTrueSecondaryConvGammaMCPt[fiCut]->Fill(Photon->Pt(),3.,fWeightJetJetMC*weightMatBudgetGamma);
        }
      }
    }
  }
  TruePhotonCandidate->SetIsTrueConvertedPhoton();
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::ProcessTruePhotonCandidates(AliAODConversionPhoton *TruePhotonCandidate)
{
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  Double_t magField = fInputEvent->GetMagneticField();
  // Double_t magFieldFlip = 1.0;

  // if( magField  < 0.0 ){
  //   magFieldFlip =  1.0;
  // }
  // else {
  //   magFieldFlip =  -1.0;
  // }

  // Process True Photons
  TParticle *posDaughter = TruePhotonCandidate->GetPositiveMCDaughter(fMCEvent);
  TParticle *negDaughter = TruePhotonCandidate->GetNegativeMCDaughter(fMCEvent);
  fCharPhotonMCInfo = 0;

  if(posDaughter == NULL || negDaughter == NULL) return; // One particle does not exist
  Int_t pdgCode[2] = {TMath::Abs(posDaughter->GetPdgCode()),TMath::Abs(negDaughter->GetPdgCode())};
  fCharPhotonMCInfo = 1;
  if(posDaughter->GetMother(0) != negDaughter->GetMother(0)){
    FillPhotonCombinatorialBackgroundHist(TruePhotonCandidate, pdgCode);
    return;
  }
  else if(posDaughter->GetMother(0) == -1){
    FillPhotonCombinatorialBackgroundHist(TruePhotonCandidate, pdgCode);
    return;
  }

  if(pdgCode[0]!=11 || pdgCode[1]!=11) return; //One Particle is not a electron

  if(posDaughter->GetPdgCode()==negDaughter->GetPdgCode()) return; // Same Charge

  TParticle *Photon = TruePhotonCandidate->GetMCParticle(fMCEvent);

  if(Photon->GetPdgCode() != 22){
    return; // Mother is no Photon
  }

  if(posDaughter->GetUniqueID() != 5 || negDaughter->GetUniqueID() !=5) return;// check if the daughters come from a conversion

  // True Photon

  Double_t weightMatBudgetGamma = 1.;
  if (fDoMaterialBudgetWeightingOfGammasForTrueMesons && ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetMaterialBudgetWeightsInitialized()) {
    weightMatBudgetGamma = ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetMaterialBudgetCorrectingWeightForTrueGamma(TruePhotonCandidate, magField);
  }

  if(fIsFromDesiredHeader){
    if(!fDoLightOutput) fHistoTrueConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
    if (CheckVectorForDoubleCount(fVectorDoubleCountTrueConvGammas,posDaughter->GetMother(0))){
      if(!fDoLightOutput) fHistoDoubleCountTrueConvGammaRPt[fiCut]->Fill(TruePhotonCandidate->GetConversionRadius(),TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
      FillMultipleCountMap(fMapMultipleCountTrueConvGammas,posDaughter->GetMother(0));
    }
  }
  Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, posDaughter->GetMother(0), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
  if(isPrimary){
    // Count just primary MC Gammas as true --> For Ratio esdtruegamma / mcconvgamma
    if(fIsFromDesiredHeader){
      fCharPhotonMCInfo = 6;
      if(!fDoLightOutput){
        fHistoTruePrimaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
        fHistoTruePrimaryConvGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt(),fWeightJetJetMC*weightMatBudgetGamma); // Allways Filled
      }
    }
    // (Not Filled for i6, Extra Signal Gamma (parambox) are secondary)
  }else {
    // fill secondary histograms
    if(fIsFromDesiredHeader){
      fCharPhotonMCInfo = 2;
      if( Photon->GetMother(0) > -1 && fMCEvent->Particle(Photon->GetMother(0))->GetMother(0) > -1){
        if (fMCEvent->Particle(fMCEvent->Particle(Photon->GetMother(0))->GetMother(0))->GetPdgCode() == 310){
          fCharPhotonMCInfo = 4;
          if (!fDoLightOutput){
            fHistoTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),0.,fWeightJetJetMC*weightMatBudgetGamma);
            fHistoTrueSecondaryConvGammaMCPt[fiCut]->Fill(Photon->Pt(),0.,fWeightJetJetMC*weightMatBudgetGamma);
            fHistoTrueSecondaryConvGammaFromXFromK0sMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
          }
        } else if (fMCEvent->Particle(fMCEvent->Particle(Photon->GetMother(0))->GetMother(0))->GetPdgCode() == 130) {
          fCharPhotonMCInfo = 7;
          if (!fDoLightOutput){
            fHistoTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),1.,fWeightJetJetMC*weightMatBudgetGamma);
            fHistoTrueSecondaryConvGammaMCPt[fiCut]->Fill(Photon->Pt(),1.,fWeightJetJetMC*weightMatBudgetGamma);
            fHistoTrueSecondaryConvGammaFromXFromK0lMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
          }
        } else if (fMCEvent->Particle(fMCEvent->Particle(Photon->GetMother(0))->GetMother(0))->GetPdgCode() == 3122) {
          fCharPhotonMCInfo = 5;
          if (!fDoLightOutput){
            fHistoTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),2.,fWeightJetJetMC*weightMatBudgetGamma);
            fHistoTrueSecondaryConvGammaMCPt[fiCut]->Fill(Photon->Pt(),2.,fWeightJetJetMC*weightMatBudgetGamma);
            fHistoTrueSecondaryConvGammaFromXFromLambdaMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
          }
        } else if (fMCEvent->Particle(fMCEvent->Particle(Photon->GetMother(0))->GetMother(0))->GetPdgCode() == 221) {
          fCharPhotonMCInfo = 3;
          if (!fDoLightOutput){
            fHistoTrueSecondaryConvGammaPt[fiCut]->Fill(Photon->Pt(),3.,fWeightJetJetMC*weightMatBudgetGamma);
            fHistoTrueSecondaryConvGammaMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3.,fWeightJetJetMC*weightMatBudgetGamma);
          }
        } else {
          if (!fDoLightOutput){
            fHistoTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3.,fWeightJetJetMC*weightMatBudgetGamma);
            fHistoTrueSecondaryConvGammaMCPt[fiCut]->Fill(Photon->Pt(),3.,fWeightJetJetMC*weightMatBudgetGamma);
          }
        }
      } else {
        if (!fDoLightOutput){
          fHistoTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3.,fWeightJetJetMC*weightMatBudgetGamma);
          fHistoTrueSecondaryConvGammaMCPt[fiCut]->Fill(Photon->Pt(),3.,fWeightJetJetMC*weightMatBudgetGamma);
        }
      }
    }
  }
  TruePhotonCandidate->SetIsTrueConvertedPhoton();
  return;
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::ProcessAODMCParticles()
{
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (fAODMCTrackArray == NULL) return;

  // Loop over all primary MC particle
  for(Long_t i = 0; i < fAODMCTrackArray->GetEntriesFast(); i++) {

    AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(i));
    if (!particle) continue;

    Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, particle, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
    if (isPrimary) {

      Int_t isMCFromMBHeader = -1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        isMCFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
        if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
      }

      if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(particle->Phi(),fEventPlaneAngle,kFALSE)) continue;
      if(!fDoLightOutput){
        if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(particle,fAODMCTrackArray,kFALSE)){
          fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // All MC Gamma
          if (TMath::Abs(particle->Eta()) < 0.66 ){
            if (particle->Phi() > 1.39626 && particle->Phi() < 3.125) fHistoMCAllGammaEMCALAccPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
          }
          if(particle->GetMother() >-1){ // Meson Decay Gamma
            switch((static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetMother())))->GetPdgCode()){
            case 111: // Pi0
              fHistoMCDecayGammaPi0Pt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
              break;
            case 113: // Rho0
              fHistoMCDecayGammaRhoPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
              break;
            case 221: // Eta
              fHistoMCDecayGammaEtaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
              break;
            case 223: // Omega
              fHistoMCDecayGammaOmegaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
              break;
            case 331: // Eta'
              fHistoMCDecayGammaEtapPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
              break;
            case 333: // Phi
              fHistoMCDecayGammaPhiPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
              break;
            case 3212: // Sigma
              fHistoMCDecayGammaSigmaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
              break;
            }
          }
        }
      }
      if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(particle,fAODMCTrackArray,kTRUE)){
        for(Int_t daughterIndex=particle->GetDaughterLabel(0);daughterIndex<=particle->GetDaughterLabel(1);daughterIndex++){
          AliAODMCParticle *tmpDaughter = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(daughterIndex));
          if(!tmpDaughter) continue;
        }
        if(!fDoLightOutput) fHistoMCConvGammaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
      }
      // Converted MC Gamma
      if(fDoMesonAnalysis){

        Double_t mesonY = 1.e30;
        Double_t ratio  = 0;
        if (particle->E() != TMath::Abs(particle->Pz())){
          ratio         = (particle->E()+particle->Pz()) / (particle->E()-particle->Pz());
        }
        if( !(ratio <= 0) ){
          mesonY = particle->Y()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
        }

        if(fDoLightOutput!=2){
          if ((mesonY > ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetRapidityCutValueMin()) && (mesonY < ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetRapidityCutValueMax())){
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
        }
        // check neutral mesons
        if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
          ->MesonIsSelectedAODMC(particle,fAODMCTrackArray,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
          AliAODMCParticle* daughter0 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetDaughterLabel(0)));
          AliAODMCParticle* daughter1 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetDaughterLabel(1)));
          Float_t weighted= 1;
          if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent)){
            if (particle->Pt()>0.005){
              weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(i, 0x0, fInputEvent);
              //                   if(particle->GetPdgCode() == 221){
              //                      cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
              //                   }
            }
          }
          Double_t mesonY = 1.e30;
          Double_t ratio  = 0;
          if (particle->E() != TMath::Abs(particle->Pz())){
            ratio         = (particle->E()+particle->Pz()) / (particle->E()-particle->Pz());
          }
          if( !(ratio <= 0) ){
            mesonY = particle->Y()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
          }
          Double_t alpha = -10;
          if (particle->GetPdgCode() == 111 || particle->GetPdgCode() == 221){
            alpha = (daughter0->E() - daughter1->E())/(daughter0->E() + daughter1->E());
          }


          if(particle->GetPdgCode() == 111){
            if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoMCPi0Pt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // All MC Pi0
            if(fDoJetAnalysis){
              if(fConvJetReader->GetTrueNJets()>0){
                if(!fDoLightOutput) fHistoMCPi0JetEventGenerated[fiCut]->Fill(particle->Pt(),weighted* fWeightJetJetMC); // MC Pi0 with gamma in acc in jet event
                fTrueVectorJetEta = fConvJetReader->GetTrueVectorJetEta();
                fTrueVectorJetPhi = fConvJetReader->GetTrueVectorJetPhi();
                Double_t RJetPi0Cand;
                for(Int_t i=0; i<fConvJetReader->GetTrueNJets(); i++){
                  Double_t DeltaEta = fTrueVectorJetEta.at(i)-particle->Eta();
                  Double_t DeltaPhi = abs(fTrueVectorJetPhi.at(i)-particle->Phi());
                  if(DeltaPhi > M_PI) {
                    DeltaPhi = 2*M_PI - DeltaPhi;
                  }
                  RJetPi0Cand = TMath::Sqrt(pow((DeltaEta),2)+pow((DeltaPhi),2));
                  if(fConvJetReader->Get_Jet_Radius() > 0 ){
                    if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                      if(!fDoLightOutput) fHistoMCPi0inJetGenerated[fiCut]->Fill(particle->Pt(),weighted* fWeightJetJetMC); // MC Pi0 with gamma in acc in a jet
                      else fHistoMCPi0Pt[fiCut]->Fill(particle->Pt(),weighted* fWeightJetJetMC);
                    }
                  }
                }
                fTrueVectorJetEta.clear();
                fTrueVectorJetPhi.clear();
              }
            }
            fHistoMCPi0WOWeightPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
            if (!fDoLightOutput){
              // fill pi0 hists against gamma pt of first leg
              fHistoMCPi0PtGammaLeg[fiCut]->Fill(daughter0->Pt(),0.,weighted*fWeightJetJetMC);
              fHistoMCPi0WOWeightPtGammaLeg[fiCut]->Fill(daughter0->Pt(),0.,fWeightJetJetMC);
              // fill pi0 hists against gamma pt of second leg
              fHistoMCPi0PtGammaLeg[fiCut]->Fill(daughter1->Pt(),1.,weighted*fWeightJetJetMC);
              fHistoMCPi0WOWeightPtGammaLeg[fiCut]->Fill(daughter1->Pt(),1.,fWeightJetJetMC);
              // fill pi0 hists against gamma pt of both legs
              fHistoMCPi0PtGammaLeg[fiCut]->Fill(daughter0->Pt(),2.,weighted*fWeightJetJetMC);
              fHistoMCPi0WOWeightPtGammaLeg[fiCut]->Fill(daughter0->Pt(),2.,fWeightJetJetMC);
              fHistoMCPi0PtGammaLeg[fiCut]->Fill(daughter1->Pt(),2.,weighted*fWeightJetJetMC);
              fHistoMCPi0WOWeightPtGammaLeg[fiCut]->Fill(daughter1->Pt(),2.,fWeightJetJetMC);
            }

            if (fIsMC > 1) fHistoMCPi0WOEvtWeightPt[fiCut]->Fill(particle->Pt());
            if (fDoMesonQA > 0){
              fHistoMCPi0PtY[fiCut]->Fill(particle->Pt(),mesonY,weighted*fWeightJetJetMC);
              fHistoMCPi0PtAlpha[fiCut]->Fill(particle->Pt(),alpha,fWeightJetJetMC);
              if (fIsMC == 2) fHistoMCPi0PtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),fWeightJetJetMC);
            }
          }else if(particle->GetPdgCode() == 221){
            if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoMCEtaPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // All MC Eta
            if(fDoJetAnalysis){
              if(fConvJetReader->GetTrueNJets()>0){
                if(!fDoLightOutput) fHistoMCEtaJetEventGenerated[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // MC Pi0 with gamma in acc in jet event
                fTrueVectorJetEta = fConvJetReader->GetTrueVectorJetEta();
                fTrueVectorJetPhi = fConvJetReader->GetTrueVectorJetPhi();
                Double_t RJetPi0Cand;
                for(Int_t i=0; i<fConvJetReader->GetTrueNJets(); i++){
                  Double_t DeltaEta = fTrueVectorJetEta.at(i)-particle->Eta();
                  Double_t DeltaPhi = abs(fTrueVectorJetPhi.at(i)-particle->Phi());
                  if(DeltaPhi > M_PI) {
                    DeltaPhi = 2*M_PI - DeltaPhi;
                  }
                  RJetPi0Cand = TMath::Sqrt(pow((DeltaEta),2)+pow((DeltaPhi),2));
                  if(fConvJetReader->Get_Jet_Radius() > 0 ){
                    if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                      if(!fDoLightOutput) fHistoMCEtainJetGenerated[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // MC Pi0 with gamma in acc in a jet
                      else fHistoMCEtaPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC);
                    }
                  }
                }
                fTrueVectorJetEta.clear();
                fTrueVectorJetPhi.clear();
              }
            }
            fHistoMCEtaWOWeightPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
            if (fIsMC > 1) fHistoMCEtaWOEvtWeightPt[fiCut]->Fill(particle->Pt());
            if (fDoMesonQA > 0){
              fHistoMCEtaPtY[fiCut]->Fill(particle->Pt(),mesonY,weighted*fWeightJetJetMC);
              fHistoMCEtaPtAlpha[fiCut]->Fill(particle->Pt(),alpha,fWeightJetJetMC);
              if (fIsMC == 2) fHistoMCEtaPtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),fWeightJetJetMC);
            }
          }

          // Check the acceptance for both gammas
          if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(daughter0,fAODMCTrackArray,kFALSE) &&
          ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(daughter1,fAODMCTrackArray,kFALSE)  &&
          ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter0->Phi(),fEventPlaneAngle,kFALSE) &&
          ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter1->Phi(),fEventPlaneAngle,kFALSE)){
            // check acceptance of clusters as well, true if one of them points into the Calo acceptance
            if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(daughter0,fAODMCTrackArray) ||
              ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(daughter1,fAODMCTrackArray) ){
              if(particle->GetPdgCode() == 111){
                if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoMCPi0InAccPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // MC Pi0 with gamma in acc
                if(fDoJetAnalysis){
                  if(fConvJetReader->GetTrueNJets()>0){
                    if(!fDoLightOutput) fHistoMCPi0JetInAccPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // MC Pi0 with gamma in acc in jet event
                    fTrueVectorJetEta = fConvJetReader->GetTrueVectorJetEta();
                    fTrueVectorJetPhi = fConvJetReader->GetTrueVectorJetPhi();
                    Double_t RJetPi0Cand;
                    for(Int_t i=0; i<fConvJetReader->GetTrueNJets(); i++){
                      Double_t DeltaEta = fTrueVectorJetEta.at(i)-particle->Eta();
                      Double_t DeltaPhi = abs(fTrueVectorJetPhi.at(i)-particle->Phi());
                      if(DeltaPhi > M_PI) {
                        DeltaPhi = 2*M_PI - DeltaPhi;
                      }
                      RJetPi0Cand = TMath::Sqrt(pow((DeltaEta),2)+pow((DeltaPhi),2));
                      if(fConvJetReader->Get_Jet_Radius() > 0 ){
                        if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                          if(!fDoLightOutput) fHistoMCPi0inJetInAccPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // MC Pi0 with gamma in acc in a jet
                          else fHistoMCPi0InAccPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC);
                        }
                      }
                    }
                    fTrueVectorJetEta.clear();
                    fTrueVectorJetPhi.clear();
                  }
                }
                fHistoMCPi0WOWeightInAccPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // MC Pi0 with gamma in acc wo weight
                if(fIsMC > 1) fHistoMCPi0WOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Pi0 with gamma in acc wo any weight
                if (!fDoLightOutput){
                  // fill pi0 hists against gamma pt of first leg
                  fHistoMCPi0InAccPtGammaLeg[fiCut]->Fill(daughter0->Pt(),0.,weighted*fWeightJetJetMC);
                  fHistoMCPi0WOWeightInAccPtGammaLeg[fiCut]->Fill(daughter0->Pt(),0.,fWeightJetJetMC);
                  // fill pi0 hists against gamma pt of second leg
                  fHistoMCPi0InAccPtGammaLeg[fiCut]->Fill(daughter1->Pt(),1.,weighted*fWeightJetJetMC);
                  fHistoMCPi0WOWeightInAccPtGammaLeg[fiCut]->Fill(daughter1->Pt(),1.,fWeightJetJetMC);
                  // fill pi0 hists against gamma pt of both legs
                  fHistoMCPi0InAccPtGammaLeg[fiCut]->Fill(daughter0->Pt(),2.,weighted*fWeightJetJetMC);
                  fHistoMCPi0WOWeightInAccPtGammaLeg[fiCut]->Fill(daughter0->Pt(),2.,fWeightJetJetMC);
                  fHistoMCPi0InAccPtGammaLeg[fiCut]->Fill(daughter1->Pt(),2.,weighted*fWeightJetJetMC);
                  fHistoMCPi0WOWeightInAccPtGammaLeg[fiCut]->Fill(daughter1->Pt(),2.,fWeightJetJetMC);
                }
              }else if(particle->GetPdgCode() == 221){
                if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoMCEtaInAccPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // MC Eta with gamma in acc
                if(fDoJetAnalysis){
                  if(fConvJetReader->GetTrueNJets()>0){
                    if(!fDoLightOutput) fHistoMCEtaJetInAccPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // MC Pi0 with gamma in acc in jet event
                    fTrueVectorJetEta = fConvJetReader->GetTrueVectorJetEta();
                    fTrueVectorJetPhi = fConvJetReader->GetTrueVectorJetPhi();
                    Double_t RJetPi0Cand;
                    for(Int_t i=0; i<fConvJetReader->GetTrueNJets(); i++){
                      Double_t DeltaEta = fTrueVectorJetEta.at(i)-particle->Eta();
                      Double_t DeltaPhi = abs(fTrueVectorJetPhi.at(i)-particle->Phi());
                      if(DeltaPhi > M_PI) {
                        DeltaPhi = 2*M_PI - DeltaPhi;
                      }
                      RJetPi0Cand = TMath::Sqrt(pow((DeltaEta),2)+pow((DeltaPhi),2));
                      if(fConvJetReader->Get_Jet_Radius() > 0 ){
                        if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                          if(!fDoLightOutput) fHistoMCEtainJetInAccPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // MC Pi0 with gamma in acc in a jet
                          else fHistoMCEtaInAccPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC);
                        }
                      }
                    }
                    fTrueVectorJetEta.clear();
                    fTrueVectorJetPhi.clear();
                  }
                }
                fHistoMCEtaWOWeightInAccPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // MC Eta with gamma in acc wo weight
                if(fIsMC > 1) fHistoMCEtaWOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Eta with gamma in acc wo any weight
              }
            }
          }
        }
      }
    } else {
      // fill secondary histograms
      Int_t isMCFromMBHeader = -1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        isMCFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
        if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
      }
      if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(particle->Phi(),fEventPlaneAngle,kFALSE)) {
        if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(particle,fAODMCTrackArray,kFALSE)){
          if (particle->GetMother() > -1) {
            AliAODMCParticle *tmpMother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetMother()));
            if (tmpMother->GetMother() > -1) {
              AliAODMCParticle *tmpGrandMother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(tmpMother->GetMother()));
              if(tmpGrandMother->GetPdgCode() == 310) {
                if (!fDoLightOutput) fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),0.,fWeightJetJetMC);
              } else if (tmpGrandMother->GetPdgCode() == 130) {
                if (!fDoLightOutput) fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),1.,fWeightJetJetMC);
              } else if (tmpGrandMother->GetPdgCode() == 3122) {
                if (!fDoLightOutput) fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),2.,fWeightJetJetMC);
              } else {
                if ( !(TMath::Abs(tmpMother->GetPdgCode()) == 11 && tmpGrandMother->GetPdgCode() == 22) )
                  if (!fDoLightOutput)  fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),3.,fWeightJetJetMC);
              }
            } else {
              if (!fDoLightOutput) fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),3.,fWeightJetJetMC);
            }
          } else {
            if (!fDoLightOutput) fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),3.,fWeightJetJetMC);
          }
        }

        if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(particle,fAODMCTrackArray,kTRUE)){
          if (particle->GetMother() > -1) {
            AliAODMCParticle *tmpMother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetMother()));
            if (tmpMother->GetMother() > -1) {
              AliAODMCParticle *tmpGrandMother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(tmpMother->GetMother()));
              if(tmpGrandMother->GetPdgCode() == 310) {
                if (!fDoLightOutput) fHistoMCSecondaryConvGammaPt[fiCut]->Fill(particle->Pt(),0.,fWeightJetJetMC);
              } else if (tmpGrandMother->GetPdgCode() == 130) {
                if (!fDoLightOutput) fHistoMCSecondaryConvGammaPt[fiCut]->Fill(particle->Pt(),1.,fWeightJetJetMC);
              } else if (tmpGrandMother->GetPdgCode() == 3122) {
                if (!fDoLightOutput) fHistoMCSecondaryConvGammaPt[fiCut]->Fill(particle->Pt(),2.,fWeightJetJetMC);
              } else {
                if ( !(TMath::Abs(tmpMother->GetPdgCode()) == 11 && tmpGrandMother->GetPdgCode() == 22) )
                  if (!fDoLightOutput) fHistoMCSecondaryConvGammaPt[fiCut]->Fill(particle->Pt(),3.,fWeightJetJetMC);
              }
            } else {
              if (!fDoLightOutput) fHistoMCSecondaryConvGammaPt[fiCut]->Fill(particle->Pt(),3.,fWeightJetJetMC);
            }
          } else {
            if (!fDoLightOutput) fHistoMCSecondaryConvGammaPt[fiCut]->Fill(particle->Pt(),3.,fWeightJetJetMC);
          }
        }
      }

      if(fDoMesonAnalysis){
        if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedAODMC(particle,fAODMCTrackArray,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
          AliAODMCParticle* daughter0 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetDaughterLabel(0)));
          AliAODMCParticle* daughter1 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetDaughterLabel(1)));
          AliAODMCParticle* mother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetMother()));
          Int_t pdgCode = mother->GetPdgCode();
          if(particle->GetPdgCode() == 111){
            Int_t source = GetSourceClassification(111,pdgCode);
            if(fHistoMCSecPi0PtvsSource && fHistoMCSecPi0PtvsSource[fiCut]) fHistoMCSecPi0PtvsSource[fiCut]->Fill(particle->Pt(),source,fWeightJetJetMC); // All MC Pi0
            if (!fDoLightOutput){
              fHistoMCSecPi0PtGamma1vsSource[fiCut]->Fill(daughter0->Pt(),source,fWeightJetJetMC); // All MC Pi0
              fHistoMCSecPi0PtGamma2vsSource[fiCut]->Fill(daughter1->Pt(),source,fWeightJetJetMC); // All MC Pi0
            }
            if(fHistoMCSecPi0Source && fHistoMCSecPi0Source[fiCut]) fHistoMCSecPi0Source[fiCut]->Fill(pdgCode);
          }else if(particle->GetPdgCode() == 221){

            if(fHistoMCSecEtaPt && fHistoMCSecEtaPt[fiCut]) fHistoMCSecEtaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // All MC Pi0
            if(fHistoMCSecEtaSource && fHistoMCSecEtaSource[fiCut]) fHistoMCSecEtaSource[fiCut]->Fill(pdgCode);
          }

          // check if conversion where within acceptance
          if( ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(daughter0,fAODMCTrackArray,kFALSE) &&
              ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(daughter1,fAODMCTrackArray,kFALSE)  &&
              ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter0->Phi(),fEventPlaneAngle,kFALSE) &&
              ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter1->Phi(),fEventPlaneAngle,kFALSE)){
              // check acceptance of clusters as well, true if one of them points into the Calo acceptance
            if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(daughter0,fAODMCTrackArray) ||
                ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(daughter1,fAODMCTrackArray) ){
              if (particle->GetPdgCode() == 111){
                Int_t source = GetSourceClassification(111,pdgCode);
                if(fHistoMCSecPi0InAccPtvsSource && fHistoMCSecPi0InAccPtvsSource[fiCut]) fHistoMCSecPi0InAccPtvsSource[fiCut]->Fill(particle->Pt(),source,fWeightJetJetMC); // All MC Pi0
                if (!fDoLightOutput){
                  fHistoMCSecPi0InAccPtGamma1vsSource[fiCut]->Fill(daughter0->Pt(),source,fWeightJetJetMC); // All MC Pi0
                  fHistoMCSecPi0InAccPtGamma2vsSource[fiCut]->Fill(daughter1->Pt(),source,fWeightJetJetMC); // All MC Pi0
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
void AliAnalysisTaskGammaConvCalo::ProcessMCParticles()
{
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();
//   cout << mcProdVtxX <<"\t" << mcProdVtxY << "\t" << mcProdVtxZ << endl;

  // Loop over all primary MC particle
  for(Long_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {
    if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, i, mcProdVtxX, mcProdVtxY, mcProdVtxZ)){
      // fill primary histograms
      TParticle* particle = (TParticle *)fMCEvent->Particle(i);
      if (!particle) continue;

      Int_t isMCFromMBHeader = -1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        isMCFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
        if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
      }


      if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(particle->Phi(),fEventPlaneAngle,kFALSE)) continue;
      if(!fDoLightOutput){
        if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCEvent,kFALSE)){
          fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // All MC Gamma
          if (TMath::Abs(particle->Eta()) < 0.66 ){
            if (particle->Phi() > 1.39626 && particle->Phi() < 3.125) fHistoMCAllGammaEMCALAccPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
          }

          if(particle->GetMother(0) >-1){ // Meson Decay Gamma
            switch(fMCEvent->Particle(particle->GetMother(0))->GetPdgCode()){
            case 111: // Pi0
              fHistoMCDecayGammaPi0Pt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
              break;
            case 113: // Rho0
              fHistoMCDecayGammaRhoPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
              break;
            case 221: // Eta
              fHistoMCDecayGammaEtaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
              break;
            case 223: // Omega
              fHistoMCDecayGammaOmegaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
              break;
            case 331: // Eta'
              fHistoMCDecayGammaEtapPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
              break;
            case 333: // Phi
              fHistoMCDecayGammaPhiPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
              break;
            case 3212: // Sigma
              fHistoMCDecayGammaSigmaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
              break;
            }
          }
        }
      }
      // Converted MC Gamma
      if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCEvent,kTRUE)){
        if(!fDoLightOutput) fHistoMCConvGammaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
      }


      if(fDoMesonAnalysis ){

        // Fill histograms for other particles
        Double_t mesonY = 1.e30;
        Double_t ratio  = 0;
        if (particle->Energy() != TMath::Abs(particle->Pz())){
          ratio         = (particle->Energy()+particle->Pz()) / (particle->Energy()-particle->Pz());
        }
        if( !(ratio <= 0) ){
          mesonY = particle->Y()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
        }

        if(fHistoMCPrimaryPtvsSource[fiCut]){
          if ((mesonY > ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetRapidityCutValueMin()) && (mesonY < ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetRapidityCutValueMax())){
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
        }
        // check neutral mesons
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
          Double_t mesonY = 1.e30;
          Double_t ratio  = 0;
          if (particle->Energy() != TMath::Abs(particle->Pz())){
            ratio         = (particle->Energy()+particle->Pz()) / (particle->Energy()-particle->Pz());
          }
          if( !(ratio <= 0) ){
            mesonY = particle->Y()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
          }

          Double_t alpha = -10;
          if (particle->GetPdgCode() == 111 || particle->GetPdgCode() == 221){
            alpha = (daughter0->Energy() - daughter1->Energy())/(daughter0->Energy() + daughter1->Energy());
          }

          if(particle->GetPdgCode() == 111){
            if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoMCPi0Pt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // All MC Pi0
            if(fDoJetAnalysis){
              if(fConvJetReader->GetTrueNJets()>0){
                if(!fDoLightOutput) fHistoMCPi0JetEventGenerated[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // MC Pi0 with gamma in acc in jet event
                fTrueVectorJetEta = fConvJetReader->GetTrueVectorJetEta();
                fTrueVectorJetPhi = fConvJetReader->GetTrueVectorJetPhi();
                Double_t RJetPi0Cand;
                for(Int_t j=0; j<fConvJetReader->GetTrueNJets(); j++){
                  Double_t DeltaEta = fTrueVectorJetEta.at(j)-particle->Eta();
                  Double_t DeltaPhi = abs(fTrueVectorJetPhi.at(j)-particle->Phi());
                  if(DeltaPhi > M_PI) {
                    DeltaPhi = 2*M_PI - DeltaPhi;
                  }
                  RJetPi0Cand = TMath::Sqrt(pow((DeltaEta),2)+pow((DeltaPhi),2));
                  if(fConvJetReader->Get_Jet_Radius() > 0 ){
                    if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                      if(!fDoLightOutput) fHistoMCPi0inJetGenerated[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // MC Pi0 with gamma in acc in a jet
                      else fHistoMCPi0Pt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC);
                    }
                  }
                }
                fTrueVectorJetEta.clear();
                fTrueVectorJetPhi.clear();
              }
            }
            fHistoMCPi0WOWeightPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
            if (!fDoLightOutput){
              // fill pi0 hists against gamma pt of first leg
              fHistoMCPi0PtGammaLeg[fiCut]->Fill(daughter0->Pt(),0.,weighted*fWeightJetJetMC);
              fHistoMCPi0WOWeightPtGammaLeg[fiCut]->Fill(daughter0->Pt(),0.,fWeightJetJetMC);
              // fill pi0 hists against gamma pt of second leg
              fHistoMCPi0PtGammaLeg[fiCut]->Fill(daughter1->Pt(),1.,weighted*fWeightJetJetMC);
              fHistoMCPi0WOWeightPtGammaLeg[fiCut]->Fill(daughter1->Pt(),1.,fWeightJetJetMC);
              // fill pi0 hists against gamma pt of both legs
              fHistoMCPi0PtGammaLeg[fiCut]->Fill(daughter0->Pt(),2.,weighted*fWeightJetJetMC);
              fHistoMCPi0WOWeightPtGammaLeg[fiCut]->Fill(daughter0->Pt(),2.,fWeightJetJetMC);
              fHistoMCPi0PtGammaLeg[fiCut]->Fill(daughter1->Pt(),2.,weighted*fWeightJetJetMC);
              fHistoMCPi0WOWeightPtGammaLeg[fiCut]->Fill(daughter1->Pt(),2.,fWeightJetJetMC);
            }

            if (fIsMC > 1) fHistoMCPi0WOEvtWeightPt[fiCut]->Fill(particle->Pt());
            if (fDoMesonQA > 0){
              fHistoMCPi0PtY[fiCut]->Fill(particle->Pt(),mesonY,weighted*fWeightJetJetMC); // All MC Pi0
              fHistoMCPi0PtAlpha[fiCut]->Fill(particle->Pt(),alpha,fWeightJetJetMC); // All MC Pi0
              if (fIsMC == 2) fHistoMCPi0PtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),fWeightJetJetMC);
            }
          }else if(particle->GetPdgCode() == 221){
            if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoMCEtaPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // All MC Eta
            if(fDoJetAnalysis){
              if(fConvJetReader->GetTrueNJets()>0){
                if(!fDoLightOutput) fHistoMCEtaJetEventGenerated[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // MC Pi0 with gamma in acc in jet event
                fTrueVectorJetEta = fConvJetReader->GetTrueVectorJetEta();
                fTrueVectorJetPhi = fConvJetReader->GetTrueVectorJetPhi();
                Double_t RJetPi0Cand;
                for(Int_t j=0; j<fConvJetReader->GetTrueNJets(); j++){
                  Double_t DeltaEta = fTrueVectorJetEta.at(j)-particle->Eta();
                  Double_t DeltaPhi = abs(fTrueVectorJetPhi.at(j)-particle->Phi());
                  if(DeltaPhi > M_PI) {
                    DeltaPhi = 2*M_PI - DeltaPhi;
                  }
                  RJetPi0Cand = TMath::Sqrt(pow((DeltaEta),2)+pow((DeltaPhi),2));
                  if(fConvJetReader->Get_Jet_Radius() > 0 ){
                    if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                      if(!fDoLightOutput) fHistoMCEtainJetGenerated[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // MC Pi0 with gamma in acc in a jet
                      else fHistoMCEtaPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC);
                    }
                  }
                }
                fTrueVectorJetEta.clear();
                fTrueVectorJetPhi.clear();
              }
            }
            fHistoMCEtaWOWeightPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
            if (fIsMC > 1) fHistoMCEtaWOEvtWeightPt[fiCut]->Fill(particle->Pt());
            if (fDoMesonQA > 0){
              fHistoMCEtaPtY[fiCut]->Fill(particle->Pt(),mesonY,weighted*fWeightJetJetMC); // All MC Pi0
              fHistoMCEtaPtAlpha[fiCut]->Fill(particle->Pt(),alpha,fWeightJetJetMC); // All MC Pi0
              if (fIsMC == 2) fHistoMCEtaPtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),fWeightJetJetMC);
            }
          }

          // Check the acceptance for both gammas & whether they are counted as primaries as well
          Bool_t kDaughter0IsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, particle->GetFirstDaughter(), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
          Bool_t kDaughter1IsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, particle->GetLastDaughter(), mcProdVtxX, mcProdVtxY, mcProdVtxZ);

          if( kDaughter0IsPrim && kDaughter1IsPrim &&
            ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(daughter0,fMCEvent,kFALSE) &&
            ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(daughter1,fMCEvent,kFALSE)  &&
            ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter0->Phi(),fEventPlaneAngle,kFALSE) &&
            ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter1->Phi(),fEventPlaneAngle,kFALSE)){
            // check acceptance of clusters as well, true if one of them points into the Calo acceptance
            if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter0,fMCEvent) ||
              ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter1,fMCEvent) ){
              if(particle->GetPdgCode() == 111){
                if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoMCPi0InAccPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // MC Pi0 with gamma in acc
                if(fDoJetAnalysis){
                  if(fConvJetReader->GetTrueNJets()>0){
                    if(!fDoLightOutput) fHistoMCPi0JetInAccPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // MC Pi0 with gamma in acc in jet event
                    fTrueVectorJetEta = fConvJetReader->GetTrueVectorJetEta();
                    fTrueVectorJetPhi = fConvJetReader->GetTrueVectorJetPhi();
                    Double_t RJetPi0Cand;
                    for(Int_t j=0; j<fConvJetReader->GetTrueNJets(); j++){
                      Double_t DeltaEta = fTrueVectorJetEta.at(j)-particle->Eta();
                      Double_t DeltaPhi = abs(fTrueVectorJetPhi.at(j)-particle->Phi());
                      if(DeltaPhi > M_PI) {
                        DeltaPhi = 2*M_PI - DeltaPhi;
                      }
                      RJetPi0Cand = TMath::Sqrt(pow((DeltaEta),2)+pow((DeltaPhi),2));
                      if(fConvJetReader->Get_Jet_Radius() > 0 ){
                        if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                          if(!fDoLightOutput) fHistoMCPi0inJetInAccPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // MC Pi0 with gamma in acc in a jet
                          else fHistoMCPi0InAccPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC);
                        }
                      }
                    }
                    fTrueVectorJetEta.clear();
                    fTrueVectorJetPhi.clear();
                  }
                }
                fHistoMCPi0WOWeightInAccPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // MC Pi0 with gamma in acc wo weighting
                if(fIsMC > 1) fHistoMCPi0WOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Pi0 with gamma in acc wo any weight
                if (!fDoLightOutput){
                  // fill pi0 hists against gamma pt of first leg
                  fHistoMCPi0InAccPtGammaLeg[fiCut]->Fill(daughter0->Pt(),0.,weighted*fWeightJetJetMC);
                  fHistoMCPi0WOWeightInAccPtGammaLeg[fiCut]->Fill(daughter0->Pt(),0.,fWeightJetJetMC);
                  // fill pi0 hists against gamma pt of second leg
                  fHistoMCPi0InAccPtGammaLeg[fiCut]->Fill(daughter1->Pt(),1.,weighted*fWeightJetJetMC);
                  fHistoMCPi0WOWeightInAccPtGammaLeg[fiCut]->Fill(daughter1->Pt(),1.,fWeightJetJetMC);
                  // fill pi0 hists against gamma pt of both legs
                  fHistoMCPi0InAccPtGammaLeg[fiCut]->Fill(daughter0->Pt(),2.,weighted*fWeightJetJetMC);
                  fHistoMCPi0WOWeightInAccPtGammaLeg[fiCut]->Fill(daughter0->Pt(),2.,fWeightJetJetMC);
                  fHistoMCPi0InAccPtGammaLeg[fiCut]->Fill(daughter1->Pt(),2.,weighted*fWeightJetJetMC);
                  fHistoMCPi0WOWeightInAccPtGammaLeg[fiCut]->Fill(daughter1->Pt(),2.,fWeightJetJetMC);
                }
              }else if(particle->GetPdgCode() == 221){
                if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoMCEtaInAccPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // MC Eta with gamma in acc
                if(fDoJetAnalysis){
                  if(fConvJetReader->GetTrueNJets()>0){
                    if(!fDoLightOutput) fHistoMCEtaJetInAccPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // MC Pi0 with gamma in acc in jet event
                    fTrueVectorJetEta = fConvJetReader->GetTrueVectorJetEta();
                    fTrueVectorJetPhi = fConvJetReader->GetTrueVectorJetPhi();
                    Double_t RJetPi0Cand;
                    for(Int_t j=0; j<fConvJetReader->GetTrueNJets(); j++){
                      Double_t DeltaEta = fTrueVectorJetEta.at(j)-particle->Eta();
                      Double_t DeltaPhi = abs(fTrueVectorJetPhi.at(j)-particle->Phi());
                      if(DeltaPhi > M_PI) {
                        DeltaPhi = 2*M_PI - DeltaPhi;
                      }
                      RJetPi0Cand = TMath::Sqrt(pow((DeltaEta),2)+pow((DeltaPhi),2));
                      if(fConvJetReader->Get_Jet_Radius() > 0 ){
                        if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                          if(!fDoLightOutput) fHistoMCEtainJetInAccPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // MC Pi0 with gamma in acc in a jet
                          else fHistoMCEtaInAccPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC);
                        }
                      }
                    }
                    fTrueVectorJetEta.clear();
                    fTrueVectorJetPhi.clear();
                  }
                }
                fHistoMCEtaWOWeightInAccPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // MC Eta with gamma in acc wo weighting
                if(fIsMC > 1) fHistoMCEtaWOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Eta with gamma in acc wo any weight
              }
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

      if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(particle->Phi(),fEventPlaneAngle,kFALSE)){
        if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCEvent,kFALSE)){
          if (particle->GetMother(0) > -1 && fMCEvent->Particle(particle->GetMother(0))->GetMother(0) > -1) {
            if (fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 310){
              if (!fDoLightOutput) fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),0.,fWeightJetJetMC);
            } else if (fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 130) {
              if (!fDoLightOutput) fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),1.,fWeightJetJetMC);
            } else if (fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 3122) {
              if (!fDoLightOutput) fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),2.,fWeightJetJetMC);
            } else {
              if ( !(TMath::Abs(fMCEvent->Particle(particle->GetMother(0))->GetPdgCode()) == 11 && fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 22) )
                if (!fDoLightOutput) fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),3.,fWeightJetJetMC);
            }
          } else {
            if (!fDoLightOutput) fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),3.,fWeightJetJetMC);
          }
        }

        if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCEvent,kTRUE)){
          if (particle->GetMother(0) > -1 && fMCEvent->Particle(particle->GetMother(0))->GetMother(0) > -1) {
            if (fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 310){
              if (!fDoLightOutput) fHistoMCSecondaryConvGammaPt[fiCut]->Fill(particle->Pt(),0.,fWeightJetJetMC);
            } else if (fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 130) {
              if (!fDoLightOutput) fHistoMCSecondaryConvGammaPt[fiCut]->Fill(particle->Pt(),1.,fWeightJetJetMC);
            } else if (fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 3122) {
              if (!fDoLightOutput) fHistoMCSecondaryConvGammaPt[fiCut]->Fill(particle->Pt(),2.,fWeightJetJetMC);
            } else {
              if ( !(TMath::Abs(fMCEvent->Particle(particle->GetMother(0))->GetPdgCode()) == 11 && fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 22) )
                if (!fDoLightOutput) fHistoMCSecondaryConvGammaPt[fiCut]->Fill(particle->Pt(),3.,fWeightJetJetMC);
            }
          } else {
            if (!fDoLightOutput) fHistoMCSecondaryConvGammaPt[fiCut]->Fill(particle->Pt(),3.,fWeightJetJetMC);
          }
        }
      }


      if(fDoMesonAnalysis){
        if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedMC(particle,fMCEvent,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
          TParticle* daughter0  = (TParticle*)fMCEvent->Particle(particle->GetFirstDaughter());
          TParticle* daughter1  = (TParticle*)fMCEvent->Particle(particle->GetLastDaughter());
          Int_t pdgCode         = ((TParticle*)fMCEvent->Particle( particle->GetFirstMother() ))->GetPdgCode();

          if(particle->GetPdgCode() == 111){
            Int_t source = GetSourceClassification(111,pdgCode);
            if(fHistoMCSecPi0PtvsSource && fHistoMCSecPi0PtvsSource[fiCut]) fHistoMCSecPi0PtvsSource[fiCut]->Fill(particle->Pt(),source,fWeightJetJetMC); // All MC Pi0
            if (!fDoLightOutput){
              fHistoMCSecPi0PtGamma1vsSource[fiCut]->Fill(daughter0->Pt(),source,fWeightJetJetMC); // All MC Pi0
              fHistoMCSecPi0PtGamma2vsSource[fiCut]->Fill(daughter1->Pt(),source,fWeightJetJetMC); // All MC Pi0
            }
            if(fHistoMCSecPi0Source && fHistoMCSecPi0Source[fiCut]) fHistoMCSecPi0Source[fiCut]->Fill(pdgCode);
          } else if(particle->GetPdgCode() == 221){
            if(fHistoMCSecEtaPt && fHistoMCSecEtaPt[fiCut]) fHistoMCSecEtaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // All MC Pi0
            if(fHistoMCSecEtaSource && fHistoMCSecEtaSource[fiCut]) fHistoMCSecEtaSource[fiCut]->Fill(pdgCode);
          }

          // check if conversion where within acceptance
          if( ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(daughter0,fMCEvent,kFALSE) &&
              ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(daughter1,fMCEvent,kFALSE)  &&
              ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter0->Phi(),fEventPlaneAngle,kFALSE) &&
              ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter1->Phi(),fEventPlaneAngle,kFALSE)){
              // check acceptance of clusters as well, true if one of them points into the Calo acceptance
            if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter0,fMCEvent) ||
                  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter1,fMCEvent) ){
              if ( particle->GetPdgCode() == 111){
                Int_t source = GetSourceClassification(111,pdgCode);
                if(fHistoMCSecPi0InAccPtvsSource && fHistoMCSecPi0InAccPtvsSource[fiCut]) fHistoMCSecPi0InAccPtvsSource[fiCut]->Fill(particle->Pt(),source,fWeightJetJetMC); // All MC Pi0
                if (!fDoLightOutput){
                  fHistoMCSecPi0InAccPtGamma1vsSource[fiCut]->Fill(daughter0->Pt(),source,fWeightJetJetMC); // All MC Pi0
                  fHistoMCSecPi0InAccPtGamma2vsSource[fiCut]->Fill(daughter1->Pt(),source,fWeightJetJetMC); // All MC Pi0
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
void AliAnalysisTaskGammaConvCalo::CalculatePi0Candidates(){
  TClonesArray * arrClustersMesonCand = NULL;
  if(fCorrTaskSetting.CompareTo(""))
    arrClustersMesonCand = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));
  // Conversion Gammas
  if(fGammaCandidates->GetEntries()>0){
    for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries();firstGammaIndex++){
      AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));
      if (gamma0==NULL) continue;

      for(Int_t secondGammaIndex=0;secondGammaIndex<fClusterCandidates->GetEntries();secondGammaIndex++){
        Bool_t matched = kFALSE;
        AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(secondGammaIndex));
        if (gamma1==NULL) continue;

        if (gamma1->GetIsCaloPhoton() > 0){
          AliVCluster* cluster = NULL;
          if(fInputEvent->IsA()==AliESDEvent::Class()){
            if(arrClustersMesonCand)
              cluster = new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersMesonCand->At(gamma1->GetCaloClusterRef()));
            else
              cluster = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef());
          } else if(fInputEvent->IsA()==AliAODEvent::Class()){
            if(arrClustersMesonCand)
              cluster = new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersMesonCand->At(gamma1->GetCaloClusterRef()));
            else
              cluster = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef());
          }

          matched = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gamma0,cluster, fInputEvent, fWeightJetJetMC);
          if(fDoConvGammaShowerShapeTree && matched){
            Float_t clusPos[3]={0,0,0};
            cluster->GetPosition(clusPos);
            TVector3 clusterVector(clusPos[0],clusPos[1],clusPos[2]);
            tESDClusE = cluster->E();
            tESDGammaConvR = gamma0->GetConversionRadius();
            tESDClusterM02 = cluster->GetM02();
            tESDClusterM20 = cluster->GetM20();
            tESDClusterEta = clusterVector.Eta();
            tESDClusterPhi = clusterVector.Phi();
            tESDClusterNCells = cluster->GetNCells();
            tESDClusterMaxECell = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FindLargestCellInCluster(cluster, fInputEvent);
            tESDClusterNLM = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetNumberOfLocalMaxima(cluster, fInputEvent);
            tESDGammaERM02[fiCut]->Fill();
          }
          if(arrClustersMesonCand) delete cluster;
        }

        if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 2){
          if ( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsSpecialTrigger()==6 ){
            if (fCaloTriggerMimicHelper[fiCut]){
              fHistoGoodMesonClusters[fiCut]->Fill(1); //"All Meson Candidates"
              if ( !(fCaloTriggerMimicHelper[fiCut]->IsClusterIDTriggered(gamma1->GetCaloClusterRef())) ){
                fHistoGoodMesonClusters[fiCut]->Fill(3); //"Cluster Not Triggered"
                if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetTriggerMimicking() != 4) {
                    continue;
                }
              }
              fHistoGoodMesonClusters[fiCut]->Fill(2); //"Triggered Meson Candidates"
              Int_t ClusterIDIsInBadDDL;
              Bool_t FlagMaybeBadDDLs=((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetReduceTriggeredPhiDueBadDDLs();
              Int_t DDLIsBadIndex;
              if (FlagMaybeBadDDLs==kFALSE){ //only flag bad DDLs -> ClusterIDIsInBadDDL>=2
                DDLIsBadIndex=2;
              } else { //also flag maybe bad DDLs  -> ClusterIDIsInBadDDL>=1
                DDLIsBadIndex=1;
              }
              ClusterIDIsInBadDDL=fCaloTriggerMimicHelper[fiCut]->IsTriggeredClusterIDInBadDDL(gamma1->GetCaloClusterRef());
              if (ClusterIDIsInBadDDL>=DDLIsBadIndex){ //DDL is bad
                  fHistoGoodMesonClusters[fiCut]->Fill(7); //"DDL not passed"
                  continue;
              } else { //DDL is good
                  fHistoGoodMesonClusters[fiCut]->Fill(6); //"DDL passed"
              }
            }
          }
        }
        if( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoGammaMinEnergyCut() ){
          Int_t minDaughters        = ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetNDaughterEnergyCut();
          Float_t minDaughterEnergy = ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetSingleDaughterMinE();
          if(minDaughters==1){ // calo over threshold
             if( (gamma1->E() < minDaughterEnergy) ) {
                 if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 2){
                   if ( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsSpecialTrigger()==6 ){
                     fHistoGoodMesonClusters[fiCut]->Fill(5); //"Cluster E not passed"
                   }
                 }
                 continue;
             }
          }
        }
        if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 2){
          if ( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsSpecialTrigger()==6 ){
            fHistoGoodMesonClusters[fiCut]->Fill(4); //"Cluster E passed"
          }
        }


        AliAODConversionMother *pi0cand = new AliAODConversionMother(gamma0,gamma1);
        pi0cand->SetLabels(firstGammaIndex,secondGammaIndex);

        if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
          if (matched){
            if(!fDoLightOutput) fHistoMotherMatchedInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(),fWeightJetJetMC);
          }else {
            if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) {
    	         if(!fDoJetAnalysis){
        		      if( fIsMC == 0 ){
        		        fHistoMotherInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(),fWeightJetJetMC);
        		      }  // In case of analysing MC, we fill also this one in the True, for photons with Material Budget Weights
    	          } else{
    		          fHistoMotherInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(),fWeightJetJetMC);
    	          }
    	    }
            if(fDoJetAnalysis){
              if(fConvJetReader->GetNJets()>0){
                fVectorJetPt = fConvJetReader->GetVectorJetPt();
                fVectorJetPx = fConvJetReader->GetVectorJetPx();
                fVectorJetPy = fConvJetReader->GetVectorJetPy();
                fVectorJetPz = fConvJetReader->GetVectorJetPz();
                fVectorJetEta = fConvJetReader->GetVectorJetEta();
                fVectorJetPhi = fConvJetReader->GetVectorJetPhi();
                if(!fDoLightOutput)  fHistoJetMotherInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(),fWeightJetJetMC);
                Double_t RJetPi0Cand = 0;
                if(fVectorJetPt.size() == fConvJetReader->GetNJets() && fVectorJetEta.size() == fConvJetReader->GetNJets() && fVectorJetPhi.size() == fConvJetReader->GetNJets()){
                  Int_t counter = 0;
                  for(Int_t i=0; i<fConvJetReader->GetNJets(); i++){
                    Double_t DeltaEta = fVectorJetEta.at(i)-pi0cand->Eta();
                    Double_t DeltaPhi = abs(fVectorJetPhi.at(i)-pi0cand->Phi());
                    if(DeltaPhi > M_PI) {
                      DeltaPhi = 2*M_PI - DeltaPhi;
                    }
                    RJetPi0Cand = TMath::Sqrt(pow((DeltaEta),2)+pow((DeltaPhi),2));
                    if(!fDoLightOutput){
                      fHistoEtaPhiJetPi0Cand[fiCut]->Fill(DeltaPhi,DeltaEta,fWeightJetJetMC);
                      fHistoRJetPi0Cand[fiCut]->Fill(RJetPi0Cand,pi0cand->Pt(),fWeightJetJetMC);
                    }
                    if(fConvJetReader->Get_Jet_Radius() > 0 ){
                      if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                        counter ++;
                        if(!fDoLightOutput){
                          fHistoPi0InJetMotherInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(),fWeightJetJetMC);
                          fHistoEtaPhiJetWithPi0Cand[fiCut]->Fill(DeltaPhi, DeltaEta,fWeightJetJetMC);
                          Double_t PtRatio = pi0cand->Pt()/(fVectorJetPt.at(i));
                          fHistoJetPi0PtRatio[fiCut]->Fill(PtRatio);
                          Double_t dotproduct = fVectorJetPx.at(i)*pi0cand->Px() + fVectorJetPy.at(i)*pi0cand->Py() + fVectorJetPz.at(i)*pi0cand->Pz();
                          Double_t magn = pow(fVectorJetPx.at(i),2) + pow(fVectorJetPy.at(i),2) + pow(fVectorJetPz.at(i),2);
                          Double_t z = dotproduct/magn;
                          fHistoJetFragmFunc[fiCut]->Fill(z,fVectorJetPt.at(i), fWeightJetJetMC);
                          fHistoJetFragmFuncZInvMass[fiCut]->Fill(pi0cand->M(), z, fWeightJetJetMC);
                        }
                        else fHistoMotherInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(),fWeightJetJetMC);

                        if(fDoJetQA){
                          if(fIsMC > 0 && fConvJetReader->GetTrueNJets()>0){
                            fTrueVectorJetPt = fConvJetReader->GetTrueVectorJetPt();
                            fTrueVectorJetEta = fConvJetReader->GetTrueVectorJetEta();
                            fTrueVectorJetPhi = fConvJetReader->GetTrueVectorJetPhi();
                            Double_t min = 100;
                            Int_t match = 0;
                            for(Int_t j = 0; j<fConvJetReader->GetTrueNJets(); j++){
                              Double_t R_jetjet;
                              DeltaEta = fVectorJetEta.at(i)-fTrueVectorJetEta.at(j);
                              DeltaPhi = abs(fVectorJetPhi.at(i)-fTrueVectorJetPhi.at(j));
                              if(DeltaPhi > M_PI) {
                                DeltaPhi = 2*M_PI - DeltaPhi;
                              }
                              R_jetjet = TMath::Sqrt(pow((DeltaEta),2)+pow((DeltaPhi),2));
                              if(R_jetjet < min){
                                min = R_jetjet;
                                match = j;
                              }
                            }
                            if(min < 0.3){
                              Double_t dotproduct = fVectorJetPx.at(i)*pi0cand->Px() + fVectorJetPy.at(i)*pi0cand->Py() + fVectorJetPz.at(i)*pi0cand->Pz();
                              Double_t magn = pow(fVectorJetPx.at(i),2) + pow(fVectorJetPy.at(i),2) + pow(fVectorJetPz.at(i),2);
                              Double_t z = dotproduct/magn;
                              if(fVectorJetPt.at(i) > 5){
                                fHistoUnfoldingAsData[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(),fWeightJetJetMC);
                                fHistoUnfoldingAsDataInvMassZ[fiCut]->Fill(pi0cand->M(), z, fWeightJetJetMC);
                              }
                              if(fVectorJetPt.at(i) < 5 && fTrueVectorJetPt.at(match) > 5){
                                fHistoUnfoldingMissed[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(),fWeightJetJetMC);
                                fHistoUnfoldingMissedInvMassZ[fiCut]->Fill(pi0cand->M(), z, fWeightJetJetMC);
                              }
                              if(fVectorJetPt.at(i) > 5 && fTrueVectorJetPt.at(match) < 5){
                                fHistoUnfoldingReject[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(),fWeightJetJetMC);
                                fHistoUnfoldingRejectInvMassZ[fiCut]->Fill(pi0cand->M(), z, fWeightJetJetMC);
                              }
                            }
                            fTrueVectorJetPt.clear();
                            fTrueVectorJetEta.clear();
                            fTrueVectorJetPhi.clear();
                          }
                        }
                      }
                    }
                  }
                  if(!fDoLightOutput) fHistoDoubleCounting[fiCut]->Fill(counter);
                }
              }
              fVectorJetPt.clear();
              fVectorJetPx.clear();
              fVectorJetPy.clear();
              fVectorJetPz.clear();
              fVectorJetEta.clear();
              fVectorJetPhi.clear();
            }
            if(fDoHBTHistoOutput){
              fHistoHBTOpeningAnglePt[fiCut]->Fill(pi0cand->GetOpeningAngle(),pi0cand->Pt());
              fHistoHBTDeltaEPt[fiCut]->Fill(abs(gamma0->E()-gamma1->E()),pi0cand->Pt());
            }
          }
          // fill invMass cluster shape tree if requested
          if(fDoInvMassShowerShapeTree){
            Double_t tempIM = pi0cand->M();
            if( (tempIM > 0.05 && tempIM < 0.2) || (tempIM > 0.4 && tempIM < 0.6) ){
              AliVCluster* cluster = NULL;
              if(fInputEvent->IsA()==AliESDEvent::Class()){
                if(arrClustersMesonCand)
                  cluster = new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersMesonCand->At(gamma1->GetCaloClusterRef()));
                else
                  cluster = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef());
              } else if(fInputEvent->IsA()==AliAODEvent::Class()){
                if(arrClustersMesonCand)
                  cluster = new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersMesonCand->At(gamma1->GetCaloClusterRef()));
                else
                  cluster = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef());
              }
              if(cluster->E()>1.){
                tESDIMMesonInvMass = pi0cand->M();
                tESDIMMesonPt = pi0cand->Pt();
                tESDIMClusE = cluster->E();
                tESDIMClusterM02 = cluster->GetM02();
                tESDIMClusterM20 = cluster->GetM20();
                tESDIMClusterLeadCellID = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FindLargestCellInCluster(cluster,fInputEvent);

                Double_t vertex[3] = {0};
                InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

                //determine isolation in cluster Et
                Float_t clsPos[3] = {0.,0.,0.};
                Float_t secondClsPos[3] = {0.,0.,0.};
                TLorentzVector clusterVector;
                cluster->GetPosition(clsPos);
                TVector3 clsPosVec(clsPos);

                Float_t sum_Et = 0;
                Int_t nclus = fInputEvent->GetNumberOfCaloClusters();
                for(Int_t j=0; j<nclus; j++){
                  if( tESDmapIsClusterAcceptedWithoutTrackMatch[j] != 1 ) continue;

                  AliVCluster* secondClus = NULL;
                  if(fInputEvent->IsA()==AliESDEvent::Class()){
                    if(arrClustersMesonCand)
                      secondClus = new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersMesonCand->At(j));
                    else
                      secondClus = new AliESDCaloCluster(*(AliESDCaloCluster*)fInputEvent->GetCaloCluster(j));
                  } else if(fInputEvent->IsA()==AliAODEvent::Class()){
                    if(arrClustersMesonCand)
                      secondClus = new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersMesonCand->At(j));
                    else
                      secondClus = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(j));;
                  }

                  if(!secondClus){
                    if(arrClustersMesonCand) delete secondClus;
                    continue;
                  }
                  if(secondClus->GetID() == cluster->GetID()) {
                    if(arrClustersMesonCand) delete secondClus;
                    continue;
                  }
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
                tESDIMClusterIsoSumClusterEt = sum_Et;

                //get cluster classification
                Bool_t isESD = kTRUE;
                if(fInputEvent->IsA()==AliAODEvent::Class()) isESD = kFALSE;
                if(fIsMC > 0) tESDIMClusterClassification = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClassifyClusterForTMEffi(cluster,fInputEvent,fMCEvent,isESD);

                //determine dEta/dPhi of cluster to closest track
                Int_t labelTrackMatch = -1;
                AliVTrack* currTrack = 0x0;
                if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClosestMatchedTrackToCluster(fInputEvent,cluster,labelTrackMatch)){
                  currTrack  = dynamic_cast<AliVTrack*>(fInputEvent->GetTrack(labelTrackMatch));
                  if(currTrack){
                    Float_t tempEta = -99999;
                    Float_t tempPhi = -99999;
                    ((AliCaloTrackMatcher*)((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetCaloTrackMatcherInstance())->GetTrackClusterMatchingResidual(currTrack->GetID(),cluster->GetID(),tempEta,tempPhi);
                    tESDIMClusMatchedTrackPt = currTrack->Pt();
                    tESDIMClusTrackDeltaEta = tempEta;
                    tESDIMClusTrackDeltaPhi = tempPhi;
                  }
                }else{
                  tESDIMClusMatchedTrackPt = 0.;
                  tESDIMClusTrackDeltaEta = -9999.;
                  tESDIMClusTrackDeltaPhi = -9999.;
                }

                //determine isolation in track Et
                tESDIMClusterIsoSumTrackEt = ((AliCaloTrackMatcher*)((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetCaloTrackMatcherInstance())->SumTrackEtAroundCluster(fInputEvent,cluster->GetID(),0.2);
                //remove Et from matched track
                if(currTrack){
                  TLorentzVector vecTrack;
                  vecTrack.SetPxPyPzE(currTrack->Px(),currTrack->Py(),currTrack->Pz(),currTrack->E());
                  tESDIMClusterIsoSumTrackEt -= vecTrack.Et();
                }

                tESDInvMassShowerShape[fiCut]->Fill();
              }
              if(arrClustersMesonCand) delete cluster;
            }
          }

          // fill new histograms
          if (!matched){
            if(!fDoLightOutput){
              fHistoPhotonPairPtconv[fiCut]->Fill(pi0cand->M(),gamma0->Pt(),fWeightJetJetMC);
              if(TMath::Abs(pi0cand->GetAlpha())<0.1)
                fHistoMotherInvMassPtAlpha[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(),fWeightJetJetMC);
            }

            if (fDoMesonQA > 0){
              if ( pi0cand->M() > 0.05 && pi0cand->M() < 0.17){
                fHistoMotherPi0PtY[fiCut]->Fill(pi0cand->Pt(),pi0cand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
                fHistoMotherPi0PtAlpha[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetAlpha(),fWeightJetJetMC);
                fHistoMotherPi0PtOpenAngle[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetOpeningAngle(),fWeightJetJetMC);
                fHistoMotherPi0ConvPhotonEtaPhi[fiCut]->Fill(gamma0->GetPhotonPhi(), gamma0->GetPhotonEta(),fWeightJetJetMC);
              }
              if ( pi0cand->M() > 0.45 && pi0cand->M() < 0.65){
                fHistoMotherEtaPtY[fiCut]->Fill(pi0cand->Pt(),pi0cand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
                fHistoMotherEtaPtAlpha[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetAlpha(),fWeightJetJetMC);
                fHistoMotherEtaPtOpenAngle[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetOpeningAngle(),fWeightJetJetMC);
                fHistoMotherEtaConvPhotonEtaPhi[fiCut]->Fill(gamma0->GetPhotonPhi(), gamma0->GetPhotonEta(),fWeightJetJetMC);
              }
            }
            if(fDoTHnSparse && ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoBGCalculation()){
              Int_t zbin = 0;
              Int_t mbin = 0;

              if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->BackgroundHandlerType() == 0){
                zbin = fBGHandler[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
                if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
                  mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
                }else {
                  mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fGammaCandidates->GetEntries());
                }
              }else if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->BackgroundHandlerType() !=2){
                zbin = fBGHandlerRP[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
                if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
                  mbin = fBGHandlerRP[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
                }else {
                  mbin = fBGHandlerRP[fiCut]->GetMultiplicityBinIndex(fGammaCandidates->GetEntries());
                }
              }
              Double_t sparesFill[4] = {pi0cand->M(),pi0cand->Pt(),(Double_t)zbin,(Double_t)mbin};
              fSparseMotherInvMassPtZM[fiCut]->Fill(sparesFill,1);
            }
          }

          if(fIsMC>0){
            if(fInputEvent->IsA()==AliESDEvent::Class())
              ProcessTrueMesonCandidates(pi0cand,gamma0,gamma1, matched);
            if(fInputEvent->IsA()==AliAODEvent::Class())
              ProcessTrueMesonCandidatesAOD(pi0cand,gamma0,gamma1, matched);
          }
          if (!matched){
            if (!fDoLightOutput || fDoECalibOutput){
              fHistoMotherInvMassECalib[fiCut]->Fill(pi0cand->M(),gamma1->E(),fWeightJetJetMC);
              fHistoMotherInvMassECalibPCM[fiCut]->Fill(pi0cand->M(),gamma0->E(),fWeightJetJetMC);
            }
          }
          if (fDoMesonQA == 2){
            fInvMass = pi0cand->M();
            fPt  = pi0cand->Pt();
            if (TMath::Abs(gamma0->GetDCAzToPrimVtx()) < TMath::Abs(gamma1->GetDCAzToPrimVtx())){
              fDCAzGammaMin = gamma0->GetDCAzToPrimVtx();
              fDCAzGammaMax = gamma1->GetDCAzToPrimVtx();
            } else {
              fDCAzGammaMin = gamma1->GetDCAzToPrimVtx();
              fDCAzGammaMax = gamma0->GetDCAzToPrimVtx();
            }
            iFlag = gamma0->GetPhotonQuality()+gamma1->GetPhotonQuality();//pi0cand->GetMesonQuality();
            if (fIsHeavyIon == 1 && fPt > 0.399 && fPt < 20. ) {
              if (fInvMass > 0.08 && fInvMass < 0.2) tESDMesonsInvMassPtDcazMinDcazMaxFlag[fiCut]->Fill();
              if ((fInvMass > 0.45 && fInvMass < 0.6) &&  (fPt > 0.999 && fPt < 20.) )tESDMesonsInvMassPtDcazMinDcazMaxFlag[fiCut]->Fill();
            } else if (((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetSingleElectronPtCut() < 0.04 && fPt > 0.099 && fPt < 20. )  {
              if ( (fInvMass > 0.08 && fInvMass < 0.6) ) tESDMesonsInvMassPtDcazMinDcazMaxFlag[fiCut]->Fill();
            } else if (fPt > 0.299 && fPt < 20. )  {
              if ( (fInvMass > 0.08 && fInvMass < 0.6) ) tESDMesonsInvMassPtDcazMinDcazMaxFlag[fiCut]->Fill();
            }
          }
        }
        delete pi0cand;
        pi0cand=0x0;
      }
    }
  }
}
//______________________________________________________________________
void AliAnalysisTaskGammaConvCalo::ProcessTrueMesonCandidates(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1, Bool_t matched)
{

  Double_t magField = fInputEvent->GetMagneticField();
  // obtain MC vertex
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  iMesonMCInfo = 0;
  // Process True Mesons
  if(TrueGammaCandidate0->GetV0Index()<fInputEvent->GetNumberOfV0s()){
    Bool_t isTruePi0 = kFALSE;
    Bool_t isTrueEta = kFALSE;
    Int_t gamma0MCLabel = -1;
    Int_t gamma0MotherLabel = -1;
    if (TrueGammaCandidate0->IsTrueConvertedPhoton()){
      gamma0MCLabel = TrueGammaCandidate0->GetMCParticleLabel(fMCEvent);
      if(gamma0MCLabel>-1){
        TParticle * gammaMC0 = (TParticle*)fMCEvent->Particle(gamma0MCLabel);
        gamma0MotherLabel=gammaMC0->GetFirstMother();
      }

    }
    if (TrueGammaCandidate1->GetIsCaloPhoton() == 0) AliFatal("CaloPhotonFlag has not been set. Aborting");

    Int_t gamma1MCLabel = TrueGammaCandidate1->GetCaloPhotonMCLabel(0); // get most probable MC label
    Int_t gamma1MotherLabel = -1;
    // check if

    TParticle * gammaMC1 = 0x0;
    if(gamma1MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
      // Daughters Gamma 1
      gammaMC1 = (TParticle*)fMCEvent->Particle(gamma1MCLabel);
      if (TrueGammaCandidate1->IsLargestComponentPhoton() || TrueGammaCandidate1->IsLargestComponentElectron()){    // largest component is electro magnetic
        // get mother of interest (pi0 or eta)
        if (TrueGammaCandidate1->IsLargestComponentPhoton()){                            // for photons its the direct mother
          gamma1MotherLabel=gammaMC1->GetMother(0);
        }else if (TrueGammaCandidate1->IsLargestComponentElectron()){                         // for electrons its either the direct mother or for conversions the grandmother
          if (TrueGammaCandidate1->IsConversion() && gammaMC1->GetMother(0)>-1) gamma1MotherLabel=fMCEvent->Particle(gammaMC1->GetMother(0))->GetMother(0);
          else gamma1MotherLabel=gammaMC1->GetMother(0);
        }
      }else {
        if (fDoMesonQA > 0 && fIsMC < 2) fHistoTrueMotherCaloEMNonLeadingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
    }
    Bool_t previouslyNotFoundTrueMesons = kFALSE;
    Int_t tmpGammaMotherlabel = gamma0MotherLabel;
    while (tmpGammaMotherlabel > 0) {
        if(((TParticle*)fMCEvent->Particle(tmpGammaMotherlabel))->GetPdgCode() != 111 && ((TParticle*)fMCEvent->Particle(tmpGammaMotherlabel))->GetPdgCode() != 221) {
            tmpGammaMotherlabel = ((TParticle*)fMCEvent->Particle(tmpGammaMotherlabel))->GetMother(0);
        } else {
            if (tmpGammaMotherlabel != gamma0MotherLabel) {
                previouslyNotFoundTrueMesons = kTRUE;
            }
            gamma0MotherLabel = tmpGammaMotherlabel;
            break;
        }
    }
    tmpGammaMotherlabel = gamma1MotherLabel;
    while (tmpGammaMotherlabel > 0) {
        if(((TParticle*)fMCEvent->Particle(tmpGammaMotherlabel))->GetPdgCode() != 111 && ((TParticle*)fMCEvent->Particle(tmpGammaMotherlabel))->GetPdgCode() != 221) {
            tmpGammaMotherlabel = ((TParticle*)fMCEvent->Particle(tmpGammaMotherlabel))->GetMother(0);
        } else {
            if (tmpGammaMotherlabel != gamma1MotherLabel) {
                previouslyNotFoundTrueMesons = kTRUE;
            }
            gamma1MotherLabel = tmpGammaMotherlabel;
            break;
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
    Double_t weightMatBudgetGamma = 1.;
    if (fDoMaterialBudgetWeightingOfGammasForTrueMesons && ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetMaterialBudgetWeightsInitialized()) {
      weightMatBudgetGamma = ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetMaterialBudgetCorrectingWeightForTrueGamma(TrueGammaCandidate0, magField);
    }

    //For MC we fill now the historgram here with weights (true and reconstructed efficiency will be consistent)
    if(!fDoJetAnalysis && !matched) fHistoMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);

    if(isTruePi0 || isTrueEta){// True Pion or Eta



      if (!matched){
        if (isTruePi0){
          fHistoTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
          if(previouslyNotFoundTrueMesons && fDoLightOutput!=2 && fHistoTruePi0InvMassPtAdditional[fiCut]) fHistoTruePi0InvMassPtAdditional[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
          if(fDoJetAnalysis && !fDoLightOutput){
            if(fConvJetReader->GetTrueNJets()>0){
              fHistoTruePi0JetMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),  fWeightJetJetMC*weightMatBudgetGamma);
              fTrueVectorJetEta = fConvJetReader->GetTrueVectorJetEta();
              fTrueVectorJetPhi = fConvJetReader->GetTrueVectorJetPhi();
              Double_t RJetPi0Cand;
              Int_t counter = 0;
              for(Int_t i=0; i<fConvJetReader->GetTrueNJets(); i++){
                Double_t DeltaEta = fTrueVectorJetEta.at(i)-Pi0Candidate->Eta();
                Double_t DeltaPhi = abs(fTrueVectorJetPhi.at(i)-Pi0Candidate->Phi());
                if(DeltaPhi > M_PI) {
                  DeltaPhi = 2*M_PI - DeltaPhi;
                }
                RJetPi0Cand = TMath::Sqrt(pow((DeltaEta),2)+pow((DeltaPhi),2));
                if(fConvJetReader->Get_Jet_Radius() > 0 ){
                  if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                    counter ++;
                    fHistoTruePi0InJetMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC*weightMatBudgetGamma);
                    fHistoMotherPi0inJetPtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), fWeightJetJetMC*weightMatBudgetGamma);
                    fHistoMotherPi0inJetPtPhi[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Phi(), fWeightJetJetMC*weightMatBudgetGamma);
                  }
                }
                fHistoTrueDoubleCountingPi0Jet[fiCut]->Fill(counter);
              }
              fTrueVectorJetEta.clear();
              fTrueVectorJetPhi.clear();
            }
          }
          if(!fDoLightOutput){
            fHistoTruePi0InvMassECalib[fiCut]->Fill(Pi0Candidate->M(),TrueGammaCandidate1->E(),fWeightJetJetMC*weightMatBudgetGamma);
            fHistoTruePi0InvMassECalibPCM[fiCut]->Fill(Pi0Candidate->M(),TrueGammaCandidate0->E(),fWeightJetJetMC*weightMatBudgetGamma);
            if (TrueGammaCandidate1->IsLargestComponentPhoton()) fHistoTruePi0PureGammaInvMassECalib[fiCut]->Fill(Pi0Candidate->M(),TrueGammaCandidate1->E(),fWeightJetJetMC*weightMatBudgetGamma);
          }
        }
        if (isTrueEta){
            fHistoTrueEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
            if(previouslyNotFoundTrueMesons && fDoLightOutput!=2 && fHistoTrueEtaInvMassPtAdditional[fiCut]) fHistoTrueEtaInvMassPtAdditional[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
            if(fDoJetAnalysis && !fDoLightOutput){
            if(fConvJetReader->GetTrueNJets()>0){
              fHistoTrueEtaJetMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
              Double_t RJetEtaCand;
              fTrueVectorJetEta = fConvJetReader->GetTrueVectorJetEta();
              fTrueVectorJetPhi = fConvJetReader->GetTrueVectorJetPhi();
              Int_t counter = 0;
              for(Int_t i=0; i<fConvJetReader->GetTrueNJets(); i++){
                Double_t DeltaEta = fTrueVectorJetEta.at(i)-Pi0Candidate->Eta();
                Double_t DeltaPhi = abs(fTrueVectorJetPhi.at(i)-Pi0Candidate->Phi());
                if(DeltaPhi > M_PI) {
                  DeltaPhi = 2*M_PI - DeltaPhi;
                }
                RJetEtaCand = TMath::Sqrt(pow((DeltaEta),2)+pow((DeltaPhi),2));
                if(fConvJetReader->Get_Jet_Radius() > 0 ){
                  if(RJetEtaCand < fConvJetReader->Get_Jet_Radius()){
                    counter ++;
                    fHistoTrueEtaInJetMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
                    fHistoMotherEtainJetPtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), fWeightJetJetMC*weightMatBudgetGamma);
                    fHistoMotherEtainJetPtPhi[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Phi(),fWeightJetJetMC*weightMatBudgetGamma);
                  }
                }
                fHistoTrueDoubleCountingEtaJet[fiCut]->Fill(counter);
              }
              fTrueVectorJetEta.clear();
              fTrueVectorJetPhi.clear();
            }
          }
        }
      }else{
        if (isTruePi0 && fHistoTruePi0MatchedInvMassPt && fHistoTruePi0MatchedInvMassPt[fiCut])fHistoTruePi0MatchedInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
        if (isTrueEta && fHistoTrueEtaMatchedInvMassPt && fHistoTrueEtaMatchedInvMassPt[fiCut])fHistoTrueEtaMatchedInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
      }
      if (fDoMesonQA > 0  && fIsMC < 2){
        if (TrueGammaCandidate1->IsLargestComponentPhoton() && !matched){
          if(isTruePi0) fHistoTruePi0CaloPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
          if(isTrueEta) fHistoTrueEtaCaloPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        }
        if (TrueGammaCandidate1->IsLargestComponentElectron() && !matched){
          if(isTruePi0) fHistoTruePi0CaloElectronInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
          if(isTrueEta) fHistoTrueEtaCaloElectronInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        }
        if (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion() ){
          if (isTruePi0 && !matched){
            fHistoTruePi0CaloConvertedPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
            fHistoTruePi0CaloConvPhotonConvRPt[fiCut]->Fill(fMCEvent->Particle(gamma1MCLabel)->R(), fMCEvent->Particle(gamma1MCLabel)->Pt());

            if(gamma1MCLabel>-1 && fMCEvent->Particle(gamma1MCLabel)->GetMother(0)>-1){
              Int_t secondElec = -1;
              if (fMCEvent->Particle(fMCEvent->Particle(gamma1MCLabel)->GetMother(0))->GetDaughter(0) != gamma1MCLabel)
                secondElec =  fMCEvent->Particle(fMCEvent->Particle(gamma1MCLabel)->GetMother(0))->GetDaughter(0);
              else
                secondElec =  fMCEvent->Particle(fMCEvent->Particle(gamma1MCLabel)->GetMother(0))->GetDaughter(1);

              if(secondElec>-1){
                Double_t alphaE = -1;
                if (fMCEvent->Particle(gamma1MCLabel)->Energy()+fMCEvent->Particle(secondElec)->Energy() != 0)
                  alphaE = (fMCEvent->Particle(gamma1MCLabel)->Energy() - fMCEvent->Particle(secondElec)->Energy())/
                      (fMCEvent->Particle(gamma1MCLabel)->Energy() + fMCEvent->Particle(secondElec)->Energy());

                fHistoTruePi0CaloConvPhotonConvRAlpha[fiCut]->Fill(fMCEvent->Particle(gamma1MCLabel)->R(), alphaE );
                fHistoTruePi0CaloConvPhotonPtAlpha[fiCut]->Fill(fMCEvent->Particle(gamma1MCLabel)->Pt(), alphaE );
              }
            }
          }

          if (isTrueEta && !matched){
            fHistoTrueEtaCaloConvertedPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
            fHistoTrueEtaCaloConvPhotonConvRPt[fiCut]->Fill(fMCEvent->Particle(gamma1MCLabel)->R(), fMCEvent->Particle(gamma1MCLabel)->Pt());

            if(gamma1MCLabel>-1 && fMCEvent->Particle(gamma1MCLabel)->GetMother(0)>-1){
              Int_t secondElec = -1;
              if (fMCEvent->Particle(fMCEvent->Particle(gamma1MCLabel)->GetMother(0))->GetDaughter(0) != gamma1MCLabel)
                secondElec =  fMCEvent->Particle(fMCEvent->Particle(gamma1MCLabel)->GetMother(0))->GetDaughter(0);
              else
                secondElec =  fMCEvent->Particle(fMCEvent->Particle(gamma1MCLabel)->GetMother(0))->GetDaughter(1);

              if(secondElec>-1){
                Double_t alphaE = -1;
                if (fMCEvent->Particle(gamma1MCLabel)->Energy()+fMCEvent->Particle(secondElec)->Energy() != 0)
                  alphaE = (fMCEvent->Particle(gamma1MCLabel)->Energy() - fMCEvent->Particle(secondElec)->Energy())/
                      (fMCEvent->Particle(gamma1MCLabel)->Energy() + fMCEvent->Particle(secondElec)->Energy());

                fHistoTrueEtaCaloConvPhotonConvRAlpha[fiCut]->Fill(fMCEvent->Particle(gamma1MCLabel)->R(), alphaE );
                fHistoTrueEtaCaloConvPhotonPtAlpha[fiCut]->Fill(fMCEvent->Particle(gamma1MCLabel)->Pt(), alphaE );
              }
            }
          }
          if ((TrueGammaCandidate0->GetMCLabelPositive() == gamma1MCLabel || TrueGammaCandidate0->GetMCLabelNegative() == gamma1MCLabel) && isTruePi0){
            fHistoTruePi0CaloConvertedPhotonMatchedInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
          }
          if ((TrueGammaCandidate0->GetMCLabelPositive() == gamma1MCLabel || TrueGammaCandidate0->GetMCLabelNegative() == gamma1MCLabel) && isTrueEta){
            fHistoTrueEtaCaloConvertedPhotonMatchedInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
          }
        }
        if ((TrueGammaCandidate1->IsMerged() || TrueGammaCandidate1->IsMergedPartConv() || TrueGammaCandidate1->IsDalitzMerged()) && !matched ){
          if (isTruePi0 )fHistoTruePi0CaloMergedClusterInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
          if (isTrueEta )fHistoTrueEtaCaloMergedClusterInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        }
        if (TrueGammaCandidate1->IsMergedPartConv() && !matched){
          if (isTruePi0) fHistoTruePi0CaloMergedClusterPartConvInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
          if (isTrueEta) fHistoTrueEtaCaloMergedClusterPartConvInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        }
      }
      if (!matched){
//        if(isTruePi0)
//          fHistoTruePi0NonLinearity[fiCut]->Fill(TrueGammaCandidate1->E(),gammaMC1->Energy()/TrueGammaCandidate1->E());
//        if(isTrueEta)
//          fHistoTrueEtaNonLinearity[fiCut]->Fill(TrueGammaCandidate1->E(),gammaMC1->Energy()/TrueGammaCandidate1->E());
        if (fDoMesonQA > 0){
          if (isTruePi0){
            if ( Pi0Candidate->M() > 0.05 && Pi0Candidate->M() < 0.17){
              fHistoTruePi0PtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC*weightMatBudgetGamma);
              fHistoTruePi0PtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetAlpha(),fWeightJetJetMC*weightMatBudgetGamma);
              fHistoTruePi0PtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle(),fWeightJetJetMC*weightMatBudgetGamma);
              fHistoTrueMotherPi0ConvPhotonEtaPhi[fiCut]->Fill(TrueGammaCandidate0->GetPhotonPhi(), TrueGammaCandidate0->GetPhotonEta(),fWeightJetJetMC*weightMatBudgetGamma);
            }
          }else if (isTrueEta){
            if ( Pi0Candidate->M() > 0.45 && Pi0Candidate->M() < 0.65){
              fHistoTrueEtaPtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC*weightMatBudgetGamma);
              fHistoTrueEtaPtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetAlpha(),fWeightJetJetMC*weightMatBudgetGamma);
              fHistoTrueEtaPtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle(),fWeightJetJetMC*weightMatBudgetGamma);
              fHistoTrueMotherEtaConvPhotonEtaPhi[fiCut]->Fill(TrueGammaCandidate0->GetPhotonPhi(), TrueGammaCandidate0->GetPhotonEta(),fWeightJetJetMC*weightMatBudgetGamma);
            }
          }
        }

        Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, gamma0MotherLabel, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
        if(!isPrimary && gamma0MotherLabel>-1){ // Secondary Meson
          Int_t secMotherLabel = ((TParticle*)fMCEvent->Particle(gamma0MotherLabel))->GetMother(0);
          Float_t weightedSec= 1;
          if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(secMotherLabel, fMCEvent, fInputEvent) && fMCEvent->Particle(secMotherLabel)->GetPdgCode()==310){
            weightedSec= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(secMotherLabel, fMCEvent, fInputEvent)/2.; //invariant mass is additive thus the weight for the daughters has to be devide by two for the K0s at a certain pt
            //cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
          }
          if (isTruePi0){
            if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTrueSecondaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
            if(fDoJetAnalysis){
              if(fConvJetReader->GetTrueNJets()>0){
                if(!fDoLightOutput) fHistoTrueSecondaryPi0InvJetMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
                fTrueVectorJetEta = fConvJetReader->GetTrueVectorJetEta();
                fTrueVectorJetPhi = fConvJetReader->GetTrueVectorJetPhi();
                Double_t RJetPi0Cand;
                  for(Int_t i=0; i<fConvJetReader->GetTrueNJets(); i++){
                    Double_t DeltaEta = fTrueVectorJetEta.at(i)-Pi0Candidate->Eta();
                    Double_t DeltaPhi = abs(fTrueVectorJetPhi.at(i)-Pi0Candidate->Phi());
                    if(DeltaPhi > M_PI) {
                      DeltaPhi = 2*M_PI - DeltaPhi;
                    }
                    RJetPi0Cand = TMath::Sqrt(pow((DeltaEta),2)+pow((DeltaPhi),2));
                    if(fConvJetReader->Get_Jet_Radius() > 0 ){
                      if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                        if(!fDoLightOutput) fHistoTrueSecondaryPi0InvinJetMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
                        else fHistoTrueSecondaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
                      }
                    }
                  }
                fTrueVectorJetEta.clear();
                fTrueVectorJetPhi.clear();
              }
            }
            if(!fDoLightOutput) {
              fHistoTrueSecondaryPi0PhotonPairPtconv[fiCut]->Fill(Pi0Candidate->M(),TrueGammaCandidate0->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
              if(CheckVectorOnly(fVectorRecTruePi0s,gamma0MotherLabel)){
                fHistoTrueSecondaryPi0DCPtconvSource[fiCut]->Fill(TrueGammaCandidate0->Pt(),0.,weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
              }
            }
          }
          iMesonMCInfo = 2;
          if (secMotherLabel >-1){
            if(fMCEvent->Particle(secMotherLabel)->GetPdgCode()==310 && isTruePi0 ){
              iMesonMCInfo = 4;
              if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTrueSecondaryPi0FromK0sInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
              if(fDoJetAnalysis){
                if(fConvJetReader->GetTrueNJets()>0){
                  if(!fDoLightOutput) fHistoTrueSecondaryPi0FromK0sJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
                  fTrueVectorJetEta = fConvJetReader->GetTrueVectorJetEta();
                  fTrueVectorJetPhi = fConvJetReader->GetTrueVectorJetPhi();
                  Double_t RJetPi0Cand;
                    for(Int_t i=0; i<fConvJetReader->GetTrueNJets(); i++){
                      Double_t DeltaEta = fTrueVectorJetEta.at(i)-Pi0Candidate->Eta();
                      Double_t DeltaPhi = abs(fTrueVectorJetPhi.at(i)-Pi0Candidate->Phi());
                      if(DeltaPhi > M_PI) {
                        DeltaPhi = 2*M_PI - DeltaPhi;
                      }
                      RJetPi0Cand = TMath::Sqrt(pow((DeltaEta),2)+pow((DeltaPhi),2));
                      if(fConvJetReader->Get_Jet_Radius() > 0 ){
                        if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                          if(!fDoLightOutput) fHistoTrueSecondaryPi0FromK0sinJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
                          else fHistoTrueSecondaryPi0FromK0sInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
                        }
                      }
                    }
                  fTrueVectorJetEta.clear();
                  fTrueVectorJetPhi.clear();
                }
              }
              if(!fDoLightOutput){
                fHistoTrueSecondaryPi0FromK0sPhotonPairPtconv[fiCut]->Fill(Pi0Candidate->M(),TrueGammaCandidate0->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
                if(CheckVectorOnly(fVectorRecTruePi0s,gamma0MotherLabel)){
                  fHistoTrueSecondaryPi0DCPtconvSource[fiCut]->Fill(TrueGammaCandidate0->Pt(),1.,weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
                }
              }
              if (fDoMesonQA > 0 && fIsMC < 2)fHistoTrueK0sWithPi0DaughterMCPt[fiCut]->Fill(fMCEvent->Particle(secMotherLabel)->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
            }
            if(fMCEvent->Particle(secMotherLabel)->GetPdgCode()==130 && isTruePi0 ){
              iMesonMCInfo = 8;
              if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTrueSecondaryPi0FromK0lInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
              if(fDoJetAnalysis){
                if(fConvJetReader->GetTrueNJets()>0){
                if(!fDoLightOutput) fHistoTrueSecondaryPi0FromK0lJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
                fTrueVectorJetEta = fConvJetReader->GetTrueVectorJetEta();
                fTrueVectorJetPhi = fConvJetReader->GetTrueVectorJetPhi();
                Double_t RJetPi0Cand;
                    for(Int_t i=0; i<fConvJetReader->GetTrueNJets(); i++){
                    Double_t DeltaEta = fTrueVectorJetEta.at(i)-Pi0Candidate->Eta();
                    Double_t DeltaPhi = abs(fTrueVectorJetPhi.at(i)-Pi0Candidate->Phi());
                    if(DeltaPhi > M_PI) {
                        DeltaPhi = 2*M_PI - DeltaPhi;
                    }
                    RJetPi0Cand = TMath::Sqrt(pow((DeltaEta),2)+pow((DeltaPhi),2));
                    if(fConvJetReader->Get_Jet_Radius() > 0 ){
                        if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                        if(!fDoLightOutput) fHistoTrueSecondaryPi0FromK0linJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
                        else fHistoTrueSecondaryPi0FromK0lInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
                        }
                    }
                    }
                fTrueVectorJetEta.clear();
                fTrueVectorJetPhi.clear();
                }
              }
              if(!fDoLightOutput){
                fHistoTrueSecondaryPi0FromK0lPhotonPairPtconv[fiCut]->Fill(Pi0Candidate->M(),TrueGammaCandidate0->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
                if(CheckVectorOnly(fVectorRecTruePi0s,gamma0MotherLabel)){
                  fHistoTrueSecondaryPi0DCPtconvSource[fiCut]->Fill(TrueGammaCandidate0->Pt(),2.,weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
                }
              }

              if (fDoMesonQA > 0 && fIsMC < 2)fHistoTrueK0lWithPi0DaughterMCPt[fiCut]->Fill(fMCEvent->Particle(secMotherLabel)->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
            }
            if(fMCEvent->Particle(secMotherLabel)->GetPdgCode()==221 && isTruePi0){
              iMesonMCInfo = 3;
              fHistoTrueSecondaryPi0FromEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
              if (fDoMesonQA > 0 && fIsMC < 2)fHistoTrueEtaWithPi0DaughterMCPt[fiCut]->Fill(fMCEvent->Particle(secMotherLabel)->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
            }
            if(fMCEvent->Particle(secMotherLabel)->GetPdgCode()==3122 && isTruePi0){
              iMesonMCInfo = 7;
              if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTrueSecondaryPi0FromLambdaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
              if(fDoJetAnalysis){
                if(fConvJetReader->GetTrueNJets()>0){
                  if(!fDoLightOutput) fHistoTrueSecondaryPi0FromLambdaJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
                  fTrueVectorJetEta = fConvJetReader->GetTrueVectorJetEta();
                  fTrueVectorJetPhi = fConvJetReader->GetTrueVectorJetPhi();
                  Double_t RJetPi0Cand;
                  for(Int_t i=0; i<fConvJetReader->GetTrueNJets(); i++){
                    Double_t DeltaEta = fTrueVectorJetEta.at(i)-Pi0Candidate->Eta();
                    Double_t DeltaPhi = abs(fTrueVectorJetPhi.at(i)-Pi0Candidate->Phi());
                    if(DeltaPhi > M_PI) {
                      DeltaPhi = 2*M_PI - DeltaPhi;
                    }
                    RJetPi0Cand = TMath::Sqrt(pow((DeltaEta),2)+pow((DeltaPhi),2));
                    if(fConvJetReader->Get_Jet_Radius() > 0 ){
                      if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                        if(!fDoLightOutput) fHistoTrueSecondaryPi0FromLambdainJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
                        else fHistoTrueSecondaryPi0FromLambdaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
                      }
                    }
                  }
                  fTrueVectorJetEta.clear();
                  fTrueVectorJetPhi.clear();
                }
              }
              if(!fDoLightOutput){
                fHistoTrueSecondaryPi0FromLambdaPhotonPairPtconv[fiCut]->Fill(Pi0Candidate->M(),TrueGammaCandidate0->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
                if(CheckVectorOnly(fVectorRecTruePi0s,gamma0MotherLabel)){
                  fHistoTrueSecondaryPi0DCPtconvSource[fiCut]->Fill(TrueGammaCandidate0->Pt(),3.,weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
                }
              }
              if (fDoMesonQA > 0 && fIsMC < 2)fHistoTrueLambdaWithPi0DaughterMCPt[fiCut]->Fill(fMCEvent->Particle(secMotherLabel)->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
            }
          }
        }else { // Only primary pi0 for efficiency calculation
          iMesonMCInfo = 6;
          Float_t weighted= 1;
          if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma1MotherLabel, fMCEvent, fInputEvent)){
            if (((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt()>0.005){
              weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(gamma1MotherLabel, fMCEvent, fInputEvent);
              //                      cout << "rec \t " <<gamma1MotherLabel << "\t" <<  weighted << endl;
            }
          }
          if (isTruePi0){
            if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTruePrimaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudgetGamma);
            if(fDoJetAnalysis){
              if(fConvJetReader->GetTrueNJets()>0){
                if(!fDoLightOutput) fHistoTruePrimaryPi0JetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudgetGamma);
                fTrueVectorJetEta = fConvJetReader->GetTrueVectorJetEta();
                fTrueVectorJetPhi = fConvJetReader->GetTrueVectorJetPhi();
                fTrueVectorJetPt  = fConvJetReader->GetTrueVectorJetPt();
                fTrueVectorJetPx  = fConvJetReader->GetTrueVectorJetPx();
                fTrueVectorJetPy  = fConvJetReader->GetTrueVectorJetPy();
                fTrueVectorJetPz  = fConvJetReader->GetTrueVectorJetPz();
                Double_t RJetPi0Cand;
                for(Int_t i=0; i<fConvJetReader->GetTrueNJets(); i++){
                  Double_t DeltaEta = fTrueVectorJetEta.at(i)-Pi0Candidate->Eta();
                  Double_t DeltaPhi = abs(fTrueVectorJetPhi.at(i)-Pi0Candidate->Phi());
                  if(DeltaPhi > M_PI) {
                    DeltaPhi = 2*M_PI - DeltaPhi;
                  }
                  RJetPi0Cand = TMath::Sqrt(pow((DeltaEta),2)+pow((DeltaPhi),2));
                  if(fConvJetReader->Get_Jet_Radius() > 0 ){
                    if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                      if(!fDoLightOutput){
                        fHistoTruePrimaryPi0inJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudgetGamma);
                        Double_t dotproduct = fTrueVectorJetPx.at(i)*Pi0Candidate->Px() + fTrueVectorJetPy.at(i)*Pi0Candidate->Py() + fTrueVectorJetPz.at(i)*Pi0Candidate->Pz();
                        Double_t magn = pow(fTrueVectorJetPx.at(i),2) + pow(fTrueVectorJetPy.at(i),2) + pow(fTrueVectorJetPz.at(i),2);
                        Double_t z = dotproduct/magn;
                        fHistoTruePi0JetFragmFunc[fiCut]->Fill(z,fTrueVectorJetPt.at(i), weighted*fWeightJetJetMC*weightMatBudgetGamma);
                        fHistoTruePi0JetFragmFuncZInvMass[fiCut]->Fill(Pi0Candidate->M(), z, weighted*fWeightJetJetMC*weightMatBudgetGamma);
                      }
                      else fHistoTruePrimaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudgetGamma);
                    }
                  }
                }
                fTrueVectorJetEta.clear();
                fTrueVectorJetPhi.clear();
                fTrueVectorJetPt.clear();
                fTrueVectorJetPx.clear();
                fTrueVectorJetPy.clear();
                fTrueVectorJetPz.clear();
              }
            }
            if(!fDoLightOutput) {
              fHistoTruePrimaryPi0PhotonPairPtconv[fiCut]->Fill(Pi0Candidate->M(),TrueGammaCandidate0->Pt(),weighted*fWeightJetJetMC*weightMatBudgetGamma);
              fHistoTruePrimaryPi0W0WeightsPhotonPairPtconv[fiCut]->Fill(Pi0Candidate->M(),TrueGammaCandidate0->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
            }
            fHistoTruePrimaryPi0W0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
            fProfileTruePrimaryPi0WeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudgetGamma);
            if (CheckVectorForDoubleCount(fVectorDoubleCountTruePi0s,gamma0MotherLabel)){
              if(fHistoDoubleCountTruePi0InvMassPt && fHistoDoubleCountTruePi0InvMassPt[fiCut]) fHistoDoubleCountTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudgetGamma);
              FillMultipleCountMap(fMapMultipleCountTruePi0s,gamma0MotherLabel);
            }
            if(!fDoLightOutput && CheckVectorOnly(fVectorRecTruePi0s,gamma0MotherLabel)){
              fHistoTruePrimaryPi0DCPtconv[fiCut]->Fill(TrueGammaCandidate0->Pt(),weighted*fWeightJetJetMC*weightMatBudgetGamma);
            }

          }else if (isTrueEta) {
            if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTruePrimaryEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudgetGamma);
            if(fDoJetAnalysis){
              if(fConvJetReader->GetTrueNJets()>0){
                if(!fDoLightOutput) fHistoTruePrimaryEtaJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudgetGamma);
                fTrueVectorJetEta = fConvJetReader->GetTrueVectorJetEta();
                fTrueVectorJetPhi = fConvJetReader->GetTrueVectorJetPhi();
                fTrueVectorJetPt  = fConvJetReader->GetTrueVectorJetPt();
                fTrueVectorJetPx  = fConvJetReader->GetTrueVectorJetPx();
                fTrueVectorJetPy  = fConvJetReader->GetTrueVectorJetPy();
                fTrueVectorJetPz  = fConvJetReader->GetTrueVectorJetPz();
                Double_t RJetPi0Cand;
                for(Int_t i=0; i<fConvJetReader->GetTrueNJets(); i++){
                  Double_t DeltaEta = fTrueVectorJetEta.at(i)-Pi0Candidate->Eta();
                  Double_t DeltaPhi = abs(fTrueVectorJetPhi.at(i)-Pi0Candidate->Phi());
                  if(DeltaPhi > M_PI) {
                    DeltaPhi = 2*M_PI - DeltaPhi;
                  }
                  RJetPi0Cand = TMath::Sqrt(pow((DeltaEta),2)+pow((DeltaPhi),2));
                  if(fConvJetReader->Get_Jet_Radius() > 0 ){
                    if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                      if(!fDoLightOutput){
                        fHistoTruePrimaryEtainJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudgetGamma);
                        Double_t dotproduct = fTrueVectorJetPx.at(i)*Pi0Candidate->Px() + fTrueVectorJetPy.at(i)*Pi0Candidate->Py() + fTrueVectorJetPz.at(i)*Pi0Candidate->Pz();
                        Double_t magn = pow(fTrueVectorJetPx.at(i),2) + pow(fTrueVectorJetPy.at(i),2) + pow(fTrueVectorJetPz.at(i),2);
                        Double_t z = dotproduct/magn;
                        fHistoTrueEtaJetFragmFunc[fiCut]->Fill(z,fTrueVectorJetPt.at(i), weighted*fWeightJetJetMC*weightMatBudgetGamma);
                        fHistoTrueEtaJetFragmFuncZInvMass[fiCut]->Fill(Pi0Candidate->M(), z, weighted*fWeightJetJetMC*weightMatBudgetGamma);
                      }else fHistoTruePrimaryEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudgetGamma);
                    }
                  }
                }
                fTrueVectorJetEta.clear();
                fTrueVectorJetPhi.clear();
                fTrueVectorJetPt.clear();
                fTrueVectorJetPx.clear();
                fTrueVectorJetPy.clear();
                fTrueVectorJetPz.clear();
              }
            }
            fHistoTruePrimaryEtaW0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
            fProfileTruePrimaryEtaWeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudgetGamma);
            if (CheckVectorForDoubleCount(fVectorDoubleCountTrueEtas,gamma0MotherLabel)){
              if(fHistoDoubleCountTrueEtaInvMassPt && fHistoDoubleCountTrueEtaInvMassPt[fiCut]) fHistoDoubleCountTrueEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudgetGamma);
              FillMultipleCountMap(fMapMultipleCountTrueEtas,gamma0MotherLabel);
            }
          }

          if (fDoMesonQA > 0 && fIsMC < 2){
            if(isTruePi0){ // Only primary pi0 for resolution
              fHistoTruePrimaryPi0MCPtResolPt[fiCut]->Fill(((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt(),(Pi0Candidate->Pt()-((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt())/((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt(),weighted*fWeightJetJetMC*weightMatBudgetGamma);
            }
            if (isTrueEta){ // Only primary eta for resolution
              fHistoTruePrimaryEtaMCPtResolPt[fiCut]->Fill(((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt(),(Pi0Candidate->Pt()-((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt())/((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt(),weighted*fWeightJetJetMC*weightMatBudgetGamma);
            }
          }
        }
      }
    }else if(!isTruePi0 && !isTrueEta){ // Background

      if (fDoMesonQA > 1){
        if(gamma0MotherLabel>-1 && gamma1MotherLabel>-1){ // Both Tracks are Photons and have a mother but not Pi0 or Eta
          fHistoTrueBckGGInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
            iMesonMCInfo = 1;

          if(
              ((((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->GetPdgCode() == 111
              || ((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->GetPdgCode() == 221)
              && (TrueGammaCandidate1->IsMerged() || TrueGammaCandidate1->IsMergedPartConv()))
          ){
            fHistoTrueBckFullMesonContainedInOneClusterInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
          }else if( TrueGammaCandidate1->E()/Pi0Candidate->E() > 0.7 ){
            fHistoTrueBckAsymEClustersInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
          }
        }else { // No photon or without mother
          fHistoTrueBckContInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
        }
      }
    }
    if (isTrueEta && !matched){
      fVectorRecTrueEtas.push_back(gamma0MotherLabel);
    }
    if (isTruePi0 && !matched){
      fVectorRecTruePi0s.push_back(gamma0MotherLabel);
    }
  }
}

//______________________________________________________________________
void AliAnalysisTaskGammaConvCalo::ProcessTrueMesonCandidatesAOD(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1, Bool_t matched)
{
  Double_t magField = fInputEvent->GetMagneticField();
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  // Process True Mesons
  if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (fAODMCTrackArray == NULL) return;
  Bool_t isTruePi0 = kFALSE;
  Bool_t isTrueEta = kFALSE;

  AliAODMCParticle *positiveMC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelPositive()));
  AliAODMCParticle *negativeMC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelNegative()));

  iMesonMCInfo = 0;
  Int_t gamma0MCLabel = -1;
  Int_t gamma0MotherLabel = -1;
  if(!positiveMC||!negativeMC)
    return;

  if (TrueGammaCandidate0->IsTrueConvertedPhoton()){
    gamma0MCLabel = positiveMC->GetMother();
    AliAODMCParticle * gammaMC0 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma0MCLabel));
    gamma0MotherLabel=gammaMC0->GetMother();
  }

  if (TrueGammaCandidate1->GetIsCaloPhoton() == 0) AliFatal("CaloPhotonFlag has not been set. Aborting");
  Int_t gamma1MCLabel = TrueGammaCandidate1->GetCaloPhotonMCLabel(0); // get most probable MC label
  Int_t gamma1MotherLabel = -1;
    // check if

  AliAODMCParticle * gammaMC1 = 0x0;
  if(gamma1MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
    // Daughters Gamma 1
    gammaMC1 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MCLabel));
    if (TrueGammaCandidate1->IsLargestComponentPhoton() || TrueGammaCandidate1->IsLargestComponentElectron()){    // largest component is electro magnetic
      // get mother of interest (pi0 or eta)
      if (TrueGammaCandidate1->IsLargestComponentPhoton()){                            // for photons its the direct mother
        gamma1MotherLabel=gammaMC1->GetMother();
      }else if (TrueGammaCandidate1->IsLargestComponentElectron()){                         // for electrons its either the direct mother or for conversions the grandmother
        if (TrueGammaCandidate1->IsConversion()){
          AliAODMCParticle * gammaGrandMotherMC1 =  static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gammaMC1->GetMother()));
          gamma1MotherLabel=gammaGrandMotherMC1->GetMother();
        }else gamma1MotherLabel=gammaMC1->GetMother();
      }
    }else {
      if (fDoMesonQA > 0 && fIsMC < 2) fHistoTrueMotherCaloEMNonLeadingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
    }
  }
  Bool_t previouslyNotFoundTrueMesons = kFALSE;
  Int_t tmpGammaMotherlabel = gamma0MotherLabel;
  while (tmpGammaMotherlabel > 0) {
      if(((AliAODMCParticle*)fAODMCTrackArray->At(tmpGammaMotherlabel))->GetPdgCode() != 111 && ((AliAODMCParticle*)fAODMCTrackArray->At(tmpGammaMotherlabel))->GetPdgCode() != 221) {
          tmpGammaMotherlabel = ((AliAODMCParticle*)fAODMCTrackArray->At(tmpGammaMotherlabel))->GetMother();
      } else {
          if (tmpGammaMotherlabel != gamma0MotherLabel) {
              previouslyNotFoundTrueMesons = kTRUE;
          }
          gamma0MotherLabel = tmpGammaMotherlabel;
          break;
      }
  }
  tmpGammaMotherlabel = gamma1MotherLabel;
  while (tmpGammaMotherlabel > 0) {
      if(((AliAODMCParticle*)fAODMCTrackArray->At(tmpGammaMotherlabel))->GetPdgCode() != 111 && ((AliAODMCParticle*)fAODMCTrackArray->At(tmpGammaMotherlabel))->GetPdgCode() != 221) {
          tmpGammaMotherlabel = ((AliAODMCParticle*)fAODMCTrackArray->At(tmpGammaMotherlabel))->GetMother();
      } else {
          if (tmpGammaMotherlabel != gamma1MotherLabel) {
              previouslyNotFoundTrueMesons = kTRUE;
          }
          gamma1MotherLabel = tmpGammaMotherlabel;
          break;
      }
  }

  if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
    if(((AliAODMCParticle*)fAODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 111){
      isTruePi0=kTRUE;
    }
    if(((AliAODMCParticle*)fAODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 221){
      isTrueEta=kTRUE;
    }
  }

  Double_t weightMatBudgetGamma = 1.;
  if (fDoMaterialBudgetWeightingOfGammasForTrueMesons && ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetMaterialBudgetWeightsInitialized()) {
    weightMatBudgetGamma = ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetMaterialBudgetCorrectingWeightForTrueGamma(TrueGammaCandidate0, magField);
  }

  //For MC we fill now the historgram here with weights (true and reconstructed efficiency will be consistent)
  if(!fDoJetAnalysis && !matched) fHistoMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);

  if(isTruePi0 || isTrueEta){// True Pion or Eta

    if (!matched){
      if (isTruePi0){
        fHistoTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
        if(previouslyNotFoundTrueMesons && fDoLightOutput!=2 && fHistoTruePi0InvMassPtAdditional[fiCut]) fHistoTruePi0InvMassPtAdditional[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
        if(fDoJetAnalysis && !fDoLightOutput){
          if(fConvJetReader->GetTrueNJets()>0){
            fHistoTruePi0JetMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
            fTrueVectorJetEta = fConvJetReader->GetTrueVectorJetEta();
            fTrueVectorJetPhi = fConvJetReader->GetTrueVectorJetPhi();
            Double_t RJetPi0Cand;
            Int_t counter = 0;
            for(Int_t i=0; i<fConvJetReader->GetTrueNJets(); i++){
              Double_t DeltaEta = fTrueVectorJetEta.at(i)-Pi0Candidate->Eta();
              Double_t DeltaPhi = abs(fTrueVectorJetPhi.at(i)-Pi0Candidate->Phi());
              if(DeltaPhi > M_PI) {
                DeltaPhi = 2*M_PI - DeltaPhi;
              }
              RJetPi0Cand = TMath::Sqrt(pow((DeltaEta),2)+pow((DeltaPhi),2));
              if(fConvJetReader->Get_Jet_Radius() > 0 ){
                if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                  counter ++;
                  fHistoTruePi0InJetMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
                  fHistoMotherPi0inJetPtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC*weightMatBudgetGamma);
                  fHistoMotherPi0inJetPtPhi[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Phi(),fWeightJetJetMC*weightMatBudgetGamma);
                }
              }
              fHistoTrueDoubleCountingPi0Jet[fiCut]->Fill(counter);
            }
            fTrueVectorJetEta.clear();
            fTrueVectorJetPhi.clear();
          }
        }
        if(!fDoLightOutput){
          fHistoTruePi0InvMassECalib[fiCut]->Fill(Pi0Candidate->M(),TrueGammaCandidate1->E(),fWeightJetJetMC*weightMatBudgetGamma);
          fHistoTruePi0InvMassECalibPCM[fiCut]->Fill(Pi0Candidate->M(),TrueGammaCandidate0->E(),fWeightJetJetMC*weightMatBudgetGamma);
          if (TrueGammaCandidate1->IsLargestComponentPhoton()) fHistoTruePi0PureGammaInvMassECalib[fiCut]->Fill(Pi0Candidate->M(),TrueGammaCandidate1->E(),fWeightJetJetMC*weightMatBudgetGamma);
        }
      }
      if (isTrueEta){
          fHistoTrueEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
          if(previouslyNotFoundTrueMesons  && fDoLightOutput!=2 && fHistoTrueEtaInvMassPtAdditional ) fHistoTrueEtaInvMassPtAdditional[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
          if(fDoJetAnalysis && !fDoLightOutput){
            if(fConvJetReader->GetTrueNJets()>0){
              fHistoTrueEtaJetMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
              fTrueVectorJetEta = fConvJetReader->GetTrueVectorJetEta();
              fTrueVectorJetPhi = fConvJetReader->GetTrueVectorJetPhi();
              Double_t RJetPi0Cand;
              Int_t counter = 0;
              for(Int_t i=0; i<fConvJetReader->GetTrueNJets(); i++){
                Double_t DeltaEta = fTrueVectorJetEta.at(i)-Pi0Candidate->Eta();
                Double_t DeltaPhi = abs(fTrueVectorJetPhi.at(i)-Pi0Candidate->Phi());
                if(DeltaPhi > M_PI) {
                  DeltaPhi = 2*M_PI - DeltaPhi;
                }
                RJetPi0Cand = TMath::Sqrt(pow((DeltaEta),2)+pow((DeltaPhi),2));
                if(fConvJetReader->Get_Jet_Radius() > 0 ){
                  if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                    counter ++;
                    fHistoTrueEtaInJetMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
                    fHistoMotherEtainJetPtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC*weightMatBudgetGamma);
                    fHistoMotherEtainJetPtPhi[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Phi(),fWeightJetJetMC*weightMatBudgetGamma);
                  }
                }
                fHistoTrueDoubleCountingEtaJet[fiCut]->Fill(counter);
              }
              fTrueVectorJetEta.clear();
              fTrueVectorJetPhi.clear();
            }
          }
      }
    }else{
      if (isTruePi0 && fHistoTruePi0MatchedInvMassPt && fHistoTruePi0MatchedInvMassPt[fiCut])fHistoTruePi0MatchedInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
      if (isTrueEta && fHistoTrueEtaMatchedInvMassPt && fHistoTrueEtaMatchedInvMassPt[fiCut])fHistoTrueEtaMatchedInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
    }
    if (fDoMesonQA > 0 && fIsMC < 2){
      if (TrueGammaCandidate1->IsLargestComponentPhoton() && !matched){
        if (isTruePi0) fHistoTruePi0CaloPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        if (isTrueEta) fHistoTrueEtaCaloPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
      if (TrueGammaCandidate1->IsLargestComponentElectron() && !matched) {
        if (isTruePi0) fHistoTruePi0CaloElectronInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        if (isTrueEta) fHistoTrueEtaCaloElectronInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
      if (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion()){
          Double_t convProdVtx[3]={-999,-999,-999};
          gammaMC1->XvYvZv(convProdVtx);
          Double_t convR = TMath::Sqrt(convProdVtx[0]*convProdVtx[0]+convProdVtx[1]*convProdVtx[1]);
        if (isTruePi0 && !matched){
          fHistoTruePi0CaloConvertedPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
          fHistoTruePi0CaloConvPhotonConvRPt[fiCut]->Fill(convR, gammaMC1->Pt());
        }
        if (isTrueEta && !matched){
          fHistoTrueEtaCaloConvertedPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
          fHistoTrueEtaCaloConvPhotonConvRPt[fiCut]->Fill(convR, gammaMC1->Pt());
        }
        if ((TrueGammaCandidate0->GetMCLabelPositive() == gamma1MCLabel || TrueGammaCandidate0->GetMCLabelNegative() == gamma1MCLabel) && isTruePi0)
          fHistoTruePi0CaloConvertedPhotonMatchedInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        if ((TrueGammaCandidate0->GetMCLabelPositive() == gamma1MCLabel || TrueGammaCandidate0->GetMCLabelNegative() == gamma1MCLabel) && isTrueEta)
          fHistoTrueEtaCaloConvertedPhotonMatchedInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
      if ((TrueGammaCandidate1->IsMerged() || TrueGammaCandidate1->IsMergedPartConv() || TrueGammaCandidate1->IsDalitzMerged()) && !matched ){
        if (isTruePi0) fHistoTruePi0CaloMergedClusterInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        if (isTrueEta) fHistoTrueEtaCaloMergedClusterInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
      if (TrueGammaCandidate1->IsMergedPartConv() && !matched) {
        if (isTruePi0) fHistoTruePi0CaloMergedClusterPartConvInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        if (isTrueEta) fHistoTrueEtaCaloMergedClusterPartConvInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
    }

    if ( !matched){
//      if(isTruePi0)
//        fHistoTruePi0NonLinearity[fiCut]->Fill(TrueGammaCandidate1->E(),gammaMC1->E()/TrueGammaCandidate1->E());
//      if(isTrueEta)
//        fHistoTrueEtaNonLinearity[fiCut]->Fill(TrueGammaCandidate1->E(),gammaMC1->E()/TrueGammaCandidate1->E());
      if (fDoMesonQA > 0){
        if (isTruePi0){
          if ( Pi0Candidate->M() > 0.05 && Pi0Candidate->M() < 0.17){
            fHistoTruePi0PtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC*weightMatBudgetGamma);
            fHistoTruePi0PtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetAlpha(),fWeightJetJetMC*weightMatBudgetGamma);
            fHistoTruePi0PtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle(),fWeightJetJetMC*weightMatBudgetGamma);
            fHistoTrueMotherPi0ConvPhotonEtaPhi[fiCut]->Fill(TrueGammaCandidate0->GetPhotonPhi(), TrueGammaCandidate0->GetPhotonEta(),fWeightJetJetMC*weightMatBudgetGamma);
          }
        }else if (isTrueEta){
          if ( Pi0Candidate->M() > 0.45 && Pi0Candidate->M() < 0.65){
            fHistoTrueEtaPtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC*weightMatBudgetGamma);
            fHistoTrueEtaPtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetAlpha(),fWeightJetJetMC*weightMatBudgetGamma);
            fHistoTrueEtaPtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle(),fWeightJetJetMC*weightMatBudgetGamma);
            fHistoTrueMotherEtaConvPhotonEtaPhi[fiCut]->Fill(TrueGammaCandidate0->GetPhotonPhi(), TrueGammaCandidate0->GetPhotonEta(),fWeightJetJetMC*weightMatBudgetGamma);
          }
        }
      }

      Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma0MotherLabel)), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
      if(!isPrimary){ // Secondary Meson
        Int_t secMotherLabel = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MotherLabel))->GetMother();
        Float_t weightedSec= 1;
        if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(secMotherLabel, 0x0, fInputEvent) && static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->GetPdgCode()==310){
          weightedSec= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(secMotherLabel, 0x0, fInputEvent)/2.; //invariant mass is additive thus the weight for the daughters has to be devide by two for the K0s at a certain pt
          //cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
        }
        if (isTruePi0){
          if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTrueSecondaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
          if(fDoJetAnalysis){
            if(fConvJetReader->GetTrueNJets()>0){
              if(!fDoLightOutput) fHistoTrueSecondaryPi0InvJetMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
              fTrueVectorJetEta = fConvJetReader->GetTrueVectorJetEta();
              fTrueVectorJetPhi = fConvJetReader->GetTrueVectorJetPhi();
              Double_t RJetPi0Cand;
                for(Int_t i=0; i<fConvJetReader->GetTrueNJets(); i++){
                  Double_t DeltaEta = fTrueVectorJetEta.at(i)-Pi0Candidate->Eta();
                  Double_t DeltaPhi = abs(fTrueVectorJetPhi.at(i)-Pi0Candidate->Phi());
                  if(DeltaPhi > M_PI) {
                    DeltaPhi = 2*M_PI - DeltaPhi;
                  }
                  RJetPi0Cand = TMath::Sqrt(pow((DeltaEta),2)+pow((DeltaPhi),2));
                  if(fConvJetReader->Get_Jet_Radius() > 0 ){
                    if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                      if(!fDoLightOutput) fHistoTrueSecondaryPi0InvinJetMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
                      else fHistoTrueSecondaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
                    }
                  }
                }
              fTrueVectorJetEta.clear();
              fTrueVectorJetPhi.clear();
            }
          }
          if(!fDoLightOutput) {
            fHistoTrueSecondaryPi0PhotonPairPtconv[fiCut]->Fill(Pi0Candidate->M(),TrueGammaCandidate0->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
            if(CheckVectorOnly(fVectorRecTruePi0s,gamma0MotherLabel)){
              fHistoTrueSecondaryPi0DCPtconvSource[fiCut]->Fill(TrueGammaCandidate0->Pt(),0.,weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
            }
          }
        }
        iMesonMCInfo = 2;
        if (secMotherLabel >-1){
          if(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->GetPdgCode()==310 && isTruePi0){
            iMesonMCInfo = 4;
            if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTrueSecondaryPi0FromK0sInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
            if(fDoJetAnalysis){
              if(fConvJetReader->GetTrueNJets()>0){
                if(!fDoLightOutput) fHistoTrueSecondaryPi0FromK0sJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
                fTrueVectorJetEta = fConvJetReader->GetTrueVectorJetEta();
                fTrueVectorJetPhi = fConvJetReader->GetTrueVectorJetPhi();
                Double_t RJetPi0Cand;
                  for(Int_t i=0; i<fConvJetReader->GetTrueNJets(); i++){
                    Double_t DeltaEta = fTrueVectorJetEta.at(i)-Pi0Candidate->Eta();
                    Double_t DeltaPhi = abs(fTrueVectorJetPhi.at(i)-Pi0Candidate->Phi());
                    if(DeltaPhi > M_PI) {
                      DeltaPhi = 2*M_PI - DeltaPhi;
                    }
                    RJetPi0Cand = TMath::Sqrt(pow((DeltaEta),2)+pow((DeltaPhi),2));
                    if(fConvJetReader->Get_Jet_Radius() > 0 ){
                      if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                        if(!fDoLightOutput) fHistoTrueSecondaryPi0FromK0sinJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
                        else fHistoTrueSecondaryPi0FromK0sInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
                      }
                    }
                  }
                fTrueVectorJetEta.clear();
                fTrueVectorJetPhi.clear();
              }
            }
            if(!fDoLightOutput){
              fHistoTrueSecondaryPi0FromK0sPhotonPairPtconv[fiCut]->Fill(Pi0Candidate->M(),TrueGammaCandidate0->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
              if(CheckVectorOnly(fVectorRecTruePi0s,gamma0MotherLabel)){
                fHistoTrueSecondaryPi0DCPtconvSource[fiCut]->Fill(TrueGammaCandidate0->Pt(),1.,weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
              }
            }

            if (fDoMesonQA > 0 && fIsMC < 2)fHistoTrueK0sWithPi0DaughterMCPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->Pt());
          }
          if(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->GetPdgCode()==130 && isTruePi0){
            iMesonMCInfo = 8;
            if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTrueSecondaryPi0FromK0lInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
            if(fDoJetAnalysis){
              if(fConvJetReader->GetTrueNJets()>0){
                if(!fDoLightOutput) fHistoTrueSecondaryPi0FromK0lJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
                fTrueVectorJetEta = fConvJetReader->GetTrueVectorJetEta();
                fTrueVectorJetPhi = fConvJetReader->GetTrueVectorJetPhi();
                Double_t RJetPi0Cand;
                for(Int_t i=0; i<fConvJetReader->GetTrueNJets(); i++){
                  Double_t DeltaEta = fTrueVectorJetEta.at(i)-Pi0Candidate->Eta();
                  Double_t DeltaPhi = abs(fTrueVectorJetPhi.at(i)-Pi0Candidate->Phi());
                  if(DeltaPhi > M_PI) {
                    DeltaPhi = 2*M_PI - DeltaPhi;
                  }
                  RJetPi0Cand = TMath::Sqrt(pow((DeltaEta),2)+pow((DeltaPhi),2));
                  if(fConvJetReader->Get_Jet_Radius() > 0 ){
                    if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                      if(!fDoLightOutput) fHistoTrueSecondaryPi0FromK0linJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
                      else fHistoTrueSecondaryPi0FromK0lInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
                    }
                  }
                }
                fTrueVectorJetEta.clear();
                fTrueVectorJetPhi.clear();
              }
            }
            if(!fDoLightOutput){
              fHistoTrueSecondaryPi0FromK0lPhotonPairPtconv[fiCut]->Fill(Pi0Candidate->M(),TrueGammaCandidate0->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
              if(CheckVectorOnly(fVectorRecTruePi0s,gamma0MotherLabel)){
                fHistoTrueSecondaryPi0DCPtconvSource[fiCut]->Fill(TrueGammaCandidate0->Pt(),2.,weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
              }
            }
            if (fDoMesonQA > 0 && fIsMC < 2)fHistoTrueK0lWithPi0DaughterMCPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->Pt());
          }
          if(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->GetPdgCode()==221 && isTruePi0){
            iMesonMCInfo = 3;
            fHistoTrueSecondaryPi0FromEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
            if (fDoMesonQA > 0 && fIsMC < 2)fHistoTrueEtaWithPi0DaughterMCPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->Pt());
          }
          if(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->GetPdgCode()==3122 && isTruePi0){
            iMesonMCInfo = 7;
            if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTrueSecondaryPi0FromLambdaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC);
            if(fDoJetAnalysis){
              if(fConvJetReader->GetTrueNJets()>0){
                if(!fDoLightOutput) fHistoTrueSecondaryPi0FromLambdaJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
                fTrueVectorJetEta = fConvJetReader->GetTrueVectorJetEta();
                fTrueVectorJetPhi = fConvJetReader->GetTrueVectorJetPhi();
                Double_t RJetPi0Cand;
                for(Int_t i=0; i<fConvJetReader->GetTrueNJets(); i++){
                  Double_t DeltaEta = fTrueVectorJetEta.at(i)-Pi0Candidate->Eta();
                  Double_t DeltaPhi = abs(fTrueVectorJetPhi.at(i)-Pi0Candidate->Phi());
                  if(DeltaPhi > M_PI) {
                    DeltaPhi = 2*M_PI - DeltaPhi;
                  }
                  RJetPi0Cand = TMath::Sqrt(pow((DeltaEta),2)+pow((DeltaPhi),2));
                  if(fConvJetReader->Get_Jet_Radius() > 0 ){
                    if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                      if(!fDoLightOutput) fHistoTrueSecondaryPi0FromLambdainJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
                      else fHistoTrueSecondaryPi0FromLambdaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
                    }
                  }
                }
                fTrueVectorJetEta.clear();
                fTrueVectorJetPhi.clear();
              }
            }
            if(!fDoLightOutput){
              fHistoTrueSecondaryPi0FromLambdaPhotonPairPtconv[fiCut]->Fill(Pi0Candidate->M(),TrueGammaCandidate0->Pt(),weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
              if(CheckVectorOnly(fVectorRecTruePi0s,gamma0MotherLabel)){
                fHistoTrueSecondaryPi0DCPtconvSource[fiCut]->Fill(TrueGammaCandidate0->Pt(),3.,weightedSec*fWeightJetJetMC*weightMatBudgetGamma);
              }
            }
            if (fDoMesonQA > 0 && fIsMC < 2)fHistoTrueLambdaWithPi0DaughterMCPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->Pt());
          }
        }
      }else{ // Only primary pi0 for efficiency calculation
        iMesonMCInfo = 6;
        Float_t weighted= 1;
        if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma1MotherLabel, 0x0, fInputEvent)){
          if (static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MotherLabel))->Pt()>0.005){
          weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(gamma1MotherLabel, 0x0, fInputEvent);
          //                      cout << "rec \t " <<gamma1MotherLabel << "\t" <<  weighted << endl;
          }
        }
        if (isTruePi0){
          if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTruePrimaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudgetGamma);
          if(fDoJetAnalysis){
            if(fConvJetReader->GetTrueNJets()>0){
              if(!fDoLightOutput) fHistoTruePrimaryPi0JetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudgetGamma);
              fTrueVectorJetEta = fConvJetReader->GetTrueVectorJetEta();
              fTrueVectorJetPhi = fConvJetReader->GetTrueVectorJetPhi();
              fTrueVectorJetPt  = fConvJetReader->GetTrueVectorJetPt();
              fTrueVectorJetPx  = fConvJetReader->GetTrueVectorJetPx();
              fTrueVectorJetPy  = fConvJetReader->GetTrueVectorJetPy();
              fTrueVectorJetPz  = fConvJetReader->GetTrueVectorJetPz();
              Double_t RJetPi0Cand;
              for(Int_t i=0; i<fConvJetReader->GetTrueNJets(); i++){
                Double_t DeltaEta = fTrueVectorJetEta.at(i)-Pi0Candidate->Eta();
                Double_t DeltaPhi = abs(fTrueVectorJetPhi.at(i)-Pi0Candidate->Phi());
                if(DeltaPhi > M_PI) {
                  DeltaPhi = 2*M_PI - DeltaPhi;
                }
                RJetPi0Cand = TMath::Sqrt(pow((DeltaEta),2)+pow((DeltaPhi),2));
                if(fConvJetReader->Get_Jet_Radius() > 0 ){
                  if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                    if(!fDoLightOutput){
                      fHistoTruePrimaryPi0inJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudgetGamma);
                      Double_t dotproduct = fTrueVectorJetPx.at(i)*Pi0Candidate->Px() + fTrueVectorJetPy.at(i)*Pi0Candidate->Py() + fTrueVectorJetPz.at(i)*Pi0Candidate->Pz();
                      Double_t magn = pow(fTrueVectorJetPx.at(i),2) + pow(fTrueVectorJetPy.at(i),2) + pow(fTrueVectorJetPz.at(i),2);
                      Double_t z = dotproduct/magn;
                      fHistoTruePi0JetFragmFunc[fiCut]->Fill(z,fTrueVectorJetPt.at(i), weighted*fWeightJetJetMC*weightMatBudgetGamma);
                      fHistoTruePi0JetFragmFuncZInvMass[fiCut]->Fill(Pi0Candidate->M(), z, weighted*fWeightJetJetMC*weightMatBudgetGamma);
                    }else fHistoTruePrimaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudgetGamma);
                  }
                }
              }
              fTrueVectorJetEta.clear();
              fTrueVectorJetPhi.clear();
              fTrueVectorJetPt.clear();
              fTrueVectorJetPx.clear();
              fTrueVectorJetPy.clear();
              fTrueVectorJetPz.clear();
            }
          }
          fHistoTruePrimaryPi0W0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC);
          if(!fDoLightOutput){
            fHistoTruePrimaryPi0PhotonPairPtconv[fiCut]->Fill(Pi0Candidate->M(),TrueGammaCandidate0->Pt(),weighted*fWeightJetJetMC*weightMatBudgetGamma);
            fHistoTruePrimaryPi0W0WeightsPhotonPairPtconv[fiCut]->Fill(Pi0Candidate->M(),TrueGammaCandidate0->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
          }
          fProfileTruePrimaryPi0WeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudgetGamma);
          if (CheckVectorForDoubleCount(fVectorDoubleCountTruePi0s,gamma0MotherLabel)){
            if(fHistoDoubleCountTruePi0InvMassPt && fHistoDoubleCountTruePi0InvMassPt[fiCut]) fHistoDoubleCountTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudgetGamma);
            FillMultipleCountMap(fMapMultipleCountTruePi0s,gamma0MotherLabel);
          }
          if(!fDoLightOutput && CheckVectorOnly(fVectorRecTruePi0s,gamma0MotherLabel)){
            fHistoTruePrimaryPi0DCPtconv[fiCut]->Fill(TrueGammaCandidate0->Pt(),weighted*fWeightJetJetMC*weightMatBudgetGamma);
          }
        }else if (isTrueEta){
          if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTruePrimaryEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudgetGamma);
          if(fDoJetAnalysis){
            if(fConvJetReader->GetTrueNJets()>0){
              if(!fDoLightOutput) fHistoTruePrimaryEtaJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudgetGamma);
              fTrueVectorJetEta = fConvJetReader->GetTrueVectorJetEta();
              fTrueVectorJetPhi = fConvJetReader->GetTrueVectorJetPhi();
              fTrueVectorJetPt  = fConvJetReader->GetTrueVectorJetPt();
              fTrueVectorJetPx  = fConvJetReader->GetTrueVectorJetPx();
              fTrueVectorJetPy  = fConvJetReader->GetTrueVectorJetPy();
              fTrueVectorJetPz  = fConvJetReader->GetTrueVectorJetPz();
              Double_t RJetPi0Cand;
              for(Int_t i=0; i<fConvJetReader->GetTrueNJets(); i++){
                Double_t DeltaEta = fTrueVectorJetEta.at(i)-Pi0Candidate->Eta();
                Double_t DeltaPhi = abs(fTrueVectorJetPhi.at(i)-Pi0Candidate->Phi());
                if(DeltaPhi > M_PI) {
                  DeltaPhi = 2*M_PI - DeltaPhi;
                }
                RJetPi0Cand = TMath::Sqrt(pow((DeltaEta),2)+pow((DeltaPhi),2));
                if(fConvJetReader->Get_Jet_Radius() > 0 ){
                  if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                    if(!fDoLightOutput){
                      fHistoTruePrimaryEtainJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudgetGamma);
                      Double_t dotproduct = fTrueVectorJetPx.at(i)*Pi0Candidate->Px() + fTrueVectorJetPy.at(i)*Pi0Candidate->Py() + fTrueVectorJetPz.at(i)*Pi0Candidate->Pz();
                      Double_t magn = pow(fTrueVectorJetPx.at(i),2) + pow(fTrueVectorJetPy.at(i),2) + pow(fTrueVectorJetPz.at(i),2);
                      Double_t z = dotproduct/magn;
                      fHistoTrueEtaJetFragmFunc[fiCut]->Fill(z,fTrueVectorJetPt.at(i), weighted*fWeightJetJetMC*weightMatBudgetGamma);
                      fHistoTrueEtaJetFragmFuncZInvMass[fiCut]->Fill(Pi0Candidate->M(), z, weighted*fWeightJetJetMC*weightMatBudgetGamma);
                    }else fHistoTruePrimaryEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudgetGamma);
                  }
                }
              }
              fTrueVectorJetEta.clear();
              fTrueVectorJetPhi.clear();
              fTrueVectorJetPt.clear();
              fTrueVectorJetPx.clear();
              fTrueVectorJetPy.clear();
              fTrueVectorJetPz.clear();
            }
          }
          fHistoTruePrimaryEtaW0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
          fProfileTruePrimaryEtaWeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudgetGamma);
          if (CheckVectorForDoubleCount(fVectorDoubleCountTrueEtas,gamma0MotherLabel)){
            if(fHistoDoubleCountTrueEtaInvMassPt && fHistoDoubleCountTrueEtaInvMassPt[fiCut]) fHistoDoubleCountTrueEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudgetGamma);
            FillMultipleCountMap(fMapMultipleCountTrueEtas,gamma0MotherLabel);
          }
        }
        if (fDoMesonQA > 0 && fIsMC < 2){
          if(isTruePi0){ // Only primary pi0 for resolution
            fHistoTruePrimaryPi0MCPtResolPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MotherLabel))->Pt(),
                                (Pi0Candidate->Pt()-static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MotherLabel))->Pt())/static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MotherLabel))->Pt(),weighted*fWeightJetJetMC*weightMatBudgetGamma);

          }
          if (isTrueEta){ // Only primary eta for resolution
            fHistoTruePrimaryEtaMCPtResolPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MotherLabel))->Pt(),
                                (Pi0Candidate->Pt()-static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MotherLabel))->Pt())/static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MotherLabel))->Pt(),weighted*fWeightJetJetMC*weightMatBudgetGamma);
          }
        }
      }
    }
  }else if(!isTruePi0 && !isTrueEta) { // Background

    if (fDoMesonQA > 1){
      if(gamma0MotherLabel>-1 && gamma1MotherLabel>-1){ // Both Tracks are Photons and have a mother but not Pi0 or Eta
        fHistoTrueBckGGInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        iMesonMCInfo = 1;
        if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() != 2){
          if(gamma1MotherLabel < fAODMCTrackArray->GetEntries()){
            AliAODMCParticle *trackgamma1 = NULL;
            trackgamma1 = (AliAODMCParticle*)fAODMCTrackArray->At(gamma1MotherLabel);
            if(trackgamma1 != NULL){

              if(
                  ((((AliAODMCParticle*)fAODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 111
                  || ((AliAODMCParticle*)fAODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 221)
                  && (TrueGammaCandidate1->IsMerged() || TrueGammaCandidate1->IsMergedPartConv()))
              ){
                fHistoTrueBckFullMesonContainedInOneClusterInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
              }else if( TrueGammaCandidate1->E()/Pi0Candidate->E() > 0.7 ){
                fHistoTrueBckAsymEClustersInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
              }
            }
          }
        }
      }else { // No photon or without mother
        fHistoTrueBckContInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
    }
  }

  if (isTrueEta && !matched){
    fVectorRecTrueEtas.push_back(gamma0MotherLabel);
  }
  if (isTruePi0 && !matched){
    fVectorRecTruePi0s.push_back(gamma0MotherLabel);
  }
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::CalculateBackground(){

  Int_t zbin = fBGClusHandler[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
  Int_t mbin = 0;

  if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
    mbin = fBGClusHandler[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
  }else {
    mbin = fBGClusHandler[fiCut]->GetMultiplicityBinIndex(fGammaCandidates->GetEntries());
  }


  AliGammaConversionAODBGHandler::GammaConversionVertex *bgEventVertex = NULL;
  if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
    for(Int_t nEventsInBG=0;nEventsInBG<fBGClusHandler[fiCut]->GetNBGEvents();nEventsInBG++){
      AliGammaConversionAODVector *previousEventV0s = fBGClusHandler[fiCut]->GetBGGoodV0s(zbin,mbin,nEventsInBG);
      if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
        bgEventVertex = fBGClusHandler[fiCut]->GetBGEventVertex(zbin,mbin,nEventsInBG);
      }

      for(Int_t iCurrent=0;iCurrent<fGammaCandidates->GetEntries();iCurrent++){
        AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fGammaCandidates->At(iCurrent));
        for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
          AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));
          if(fMoveParticleAccordingToVertex == kTRUE){
            if (bgEventVertex){
              MoveParticleAccordingToVertex(&previousGoodV0,bgEventVertex);
            }
          }
          if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
            if (bgEventVertex){
              RotateParticleAccordingToEP(&previousGoodV0,bgEventVertex->fEP,fEventPlaneAngle);
            }
          }

          AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
          backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
          if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
            ->MesonIsSelected(backgroundCandidate,kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
            fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(),fWeightJetJetMC);
            if(!fDoLightOutput) fHistoPhotonPairMixedEventPtconv[fiCut]->Fill(backgroundCandidate->M(),currentEventGoodV0.Pt());
            if(fDoTHnSparse){
              Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
              fSparseMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
            }
           if(!fDoLightOutput || fDoECalibOutput)  fHistoMotherBackInvMassECalib[fiCut]->Fill(backgroundCandidate->M(),currentEventGoodV0.E(),fWeightJetJetMC);
           if(!fDoLightOutput || fDoECalibOutput)  fHistoMotherBackInvMassECalibPCM[fiCut]->Fill(backgroundCandidate->M(),previousGoodV0.E(),fWeightJetJetMC);
          }
          delete backgroundCandidate;
          backgroundCandidate = 0x0;
        }
      }
    }
  }else if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoJetMixing()){
    if(fDoJetAnalysis){
      if(fConvJetReader->GetNJets()>0){
        fVectorJetEta = fConvJetReader->GetVectorJetEta();
        fVectorJetPhi = fConvJetReader->GetVectorJetPhi();
        fVectorJetPt  = fConvJetReader->GetVectorJetPt();
        Double_t MaxPt = 0.;
        Int_t MaxPtPlace = 0;
        Double_t JetNearEMCal = kFALSE;
        for(Int_t i=0; i<fConvJetReader->GetNJets(); i++){
          if(fVectorJetEta.at(i) > (-0.6687 -0.4) && fVectorJetEta.at(i) < (0.66465 + 0.4)){
            if(fVectorJetPhi.at(i) > (1.39626 -0.4) && fVectorJetPhi.at(i) < (3.15 + 0.4)){ //INSIDE EMCAL
              JetNearEMCal = kTRUE;
              if(fVectorJetPt.at(i) > MaxPt){
                MaxPt = fVectorJetPt.at(i);
                MaxPtPlace = i;
              }
            }
          }
        }
        if(JetNearEMCal){
          Double_t mbinJets = 0;
          if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoJetPtMixing()){
             if(fVectorJetPt.at(MaxPtPlace) > 15) mbinJets = fBGHandler[fiCut]->GetMultiplicityBinIndex(3);
             else mbinJets = fBGHandler[fiCut]->GetMultiplicityBinIndex(2);
          }else{
             mbinJets = fBGHandler[fiCut]->GetMultiplicityBinIndex(2);
          }
          Int_t zbinJets = fBGHandler[fiCut]->GetZBinIndex(2);
          for(Int_t nEventsInBG=0;nEventsInBG <fBGHandler[fiCut]->GetNBGEvents();nEventsInBG++){
            AliGammaConversionAODVector *previousEventV0s = fBGHandler[fiCut]->GetBGGoodV0s(zbinJets,mbinJets,nEventsInBG);
            AliGammaConversionAODBGHandler::GammaConversionVertex* BGVertex = fBGHandler[fiCut]->GetBGEventVertex(zbinJets,mbinJets,nEventsInBG);
            if(previousEventV0s){
              Double_t BGJetEta = BGVertex->fX;
              Double_t BGJetPhi = BGVertex->fY;
              Int_t EtaSwap = 1;
              Bool_t DoPhiSwap = kFALSE;
              if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoJetRotateMixing()){
                if((fVectorJetEta.at(MaxPtPlace) < 0 && BGJetEta > 0) || (fVectorJetEta.at(MaxPtPlace) > 0 && BGJetEta < 0)) EtaSwap = -1;
                if(fVectorJetPhi.at(MaxPtPlace) < 2.27313 && BGJetPhi > 2.27313){
                  Double_t DeltaPhiMid = BGJetPhi - 2.27313;
                  BGJetPhi = 2.27313 - DeltaPhiMid;
                  DoPhiSwap = kTRUE;
                }
                if(fVectorJetPhi.at(MaxPtPlace) > 2.27313 && BGJetPhi < 2.27313){
                  Double_t DeltaPhiMid = 2.27313 - BGJetPhi;
                  BGJetPhi = 2.27313 + DeltaPhiMid;
                  DoPhiSwap = kTRUE;
                }
              }
              Double_t EtaShift = fVectorJetEta.at(MaxPtPlace) - BGJetEta*EtaSwap;
              Double_t PhiShift = fVectorJetPhi.at(MaxPtPlace) - BGJetPhi;
              for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
                AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));
                Double_t EtaBackgroundAdjusted = previousGoodV0.Eta()*EtaSwap + EtaShift;
                Double_t PhiBackgroundAdjusted = 0.;
                if(DoPhiSwap){
                  if(previousGoodV0.Phi() < 2.27313){
                    Double_t DeltaPhiMid = 2.27313 - previousGoodV0.Phi();
                    PhiBackgroundAdjusted = 2.27313 + DeltaPhiMid + PhiShift;
                  }
                  if(previousGoodV0.Phi() > 2.27313){
                    Double_t DeltaPhiMid = previousGoodV0.Phi() - 2.27313;
                    PhiBackgroundAdjusted = 2.27313 - DeltaPhiMid + PhiShift;
                  }
                }else PhiBackgroundAdjusted = previousGoodV0.Phi() + PhiShift;
                previousGoodV0.SetPtEtaPhiE(previousGoodV0.Pt(), EtaBackgroundAdjusted, PhiBackgroundAdjusted, previousGoodV0.E());

                for(Int_t iCurrent=0;iCurrent<fGammaCandidates->GetEntries();iCurrent++){
                  AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fGammaCandidates->At(iCurrent));
                  AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
                  backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
                  if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
                    if(!fDoLightOutput) fHistoMotherBackJetInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(),fWeightJetJetMC);
                    else fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(),fWeightJetJetMC);
                  }
                  delete backgroundCandidate;
                  backgroundCandidate = 0x0;
                }
              }
            }
          }
        }
        fVectorJetEta.clear();
        fVectorJetPhi.clear();
      }
    }
  } else {
    // mixing current conversion photons with previous clusters
    for(Int_t nEventsInBG=0;nEventsInBG <fBGClusHandler[fiCut]->GetNBGEvents();nEventsInBG++){
      AliGammaConversionAODVector *previousEventV0s = fBGClusHandler[fiCut]->GetBGGoodV0s(zbin,mbin,nEventsInBG);
      if(previousEventV0s){
        if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0 ){
          bgEventVertex = fBGClusHandler[fiCut]->GetBGEventVertex(zbin,mbin,nEventsInBG);
        }
        for(Int_t iCurrent=0;iCurrent<fGammaCandidates->GetEntries();iCurrent++){
          AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fGammaCandidates->At(iCurrent));
          for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){

            AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));

            if(fMoveParticleAccordingToVertex == kTRUE){
              if (bgEventVertex){
                MoveParticleAccordingToVertex(&previousGoodV0,bgEventVertex);
              }
            }
            if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
              if (bgEventVertex){
                RotateParticleAccordingToEP(&previousGoodV0,bgEventVertex->fEP,fEventPlaneAngle);
              }
            }

            AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
            backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
            if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate,kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
              fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(),fWeightJetJetMC);
              if(!fDoLightOutput) fHistoPhotonPairMixedEventPtconv[fiCut]->Fill(backgroundCandidate->M(),currentEventGoodV0.Pt());
              if(fDoTHnSparse){
                Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
                fSparseMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
              }
              if(!fDoLightOutput || fDoECalibOutput) fHistoMotherBackInvMassECalib[fiCut]->Fill(backgroundCandidate->M(),currentEventGoodV0.E(),fWeightJetJetMC);
              if(fDoHBTHistoOutput){
                  fHistoBckHBTOpeningAnglePt[fiCut]->Fill(backgroundCandidate->GetOpeningAngle(),backgroundCandidate->Pt());
                  fHistoBckHBTDeltaEPt[fiCut]->Fill(abs(currentEventGoodV0.E()-previousGoodV0.E()),backgroundCandidate->Pt());
              }
            }
            delete backgroundCandidate;
            backgroundCandidate = 0x0;
          }
        }
      }
    }
    if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoConvCaloMixing())){
      // mixing current clusters with previous conversion photons
      for(Int_t nEventsInBG=0;nEventsInBG <fBGHandler[fiCut]->GetNBGEvents();nEventsInBG++){
        AliGammaConversionAODVector *previousEventV0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,nEventsInBG);
        if(previousEventV0s){
          if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0 ){
            bgEventVertex = fBGHandler[fiCut]->GetBGEventVertex(zbin,mbin,nEventsInBG);
          }
          for(Int_t iCurrent=0;iCurrent<fClusterCandidates->GetEntries();iCurrent++){
            AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fClusterCandidates->At(iCurrent));
            for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){

              AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));

              if(fMoveParticleAccordingToVertex == kTRUE){
                if (bgEventVertex){
                  MoveParticleAccordingToVertex(&previousGoodV0,bgEventVertex);
                }
              }
              if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
                if (bgEventVertex){
                  RotateParticleAccordingToEP(&previousGoodV0,bgEventVertex->fEP,fEventPlaneAngle);
                }
              }

              AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
              backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
              if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate,kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
                fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(),fWeightJetJetMC);
                if(!fDoLightOutput) fHistoPhotonPairMixedEventPtconv[fiCut]->Fill(backgroundCandidate->M(),currentEventGoodV0.Pt());
                if(fDoTHnSparse){
                  Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
                  fSparseMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
                }
                if(!fDoLightOutput || fDoECalibOutput) fHistoMotherBackInvMassECalib[fiCut]->Fill(backgroundCandidate->M(),currentEventGoodV0.E(),fWeightJetJetMC);
                if(fDoHBTHistoOutput){
                  fHistoBckHBTOpeningAnglePt[fiCut]->Fill(backgroundCandidate->GetOpeningAngle(),backgroundCandidate->Pt());
                  fHistoBckHBTDeltaEPt[fiCut]->Fill(abs(currentEventGoodV0.E()-previousGoodV0.E()),backgroundCandidate->Pt());
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
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::CalculateBackgroundSwapp(){

  if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoGammaSwappForBg()) {

    Double_t rotationAngle = TMath::Pi()/2.0;

    TLorentzVector lvRotationPhoton1;   // photon candidates which get rotated
    TLorentzVector lvRotationPhoton2;   // photon candidates which get rotated
    TVector3 lvRotationPion;            // reconstructed mother particle from the two photons

    Int_t cellIDRotatedPhoton = -1; // cell ID of the cluster after rotation

    std::vector<std::array<Double_t, 2>> vSwappingInvMassPT;
    std::vector<std::array<Double_t, 2>> vSwappingInvMassPTAlphaCut;
    vSwappingInvMassPT.clear();
    vSwappingInvMassPTAlphaCut.clear();
    vSwappingInvMassPT.resize(0);
    vSwappingInvMassPTAlphaCut.resize(0);

    Double_t tempMultWeightSwapping = 1; // weight taking multiplicity of event into account

    for(Int_t iCurrent1=0;iCurrent1<fClusterCandidates->GetEntries();iCurrent1++){
      AliAODConversionPhoton* currentEventGoodV0Temp1 = (AliAODConversionPhoton*)(fClusterCandidates->At(iCurrent1));

      for(Int_t iCurrent2=0;iCurrent2<fGammaCandidates->GetEntries();iCurrent2++){
        AliAODConversionPhoton *currentEventGoodV0Temp2=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(iCurrent2));
        if (currentEventGoodV0Temp2==NULL) continue;

        for(int iSwapp = 0; iSwapp < ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetNumberOfSwappsForBg(); ++iSwapp){

          lvRotationPhoton1.SetX(currentEventGoodV0Temp1->Px());
          lvRotationPhoton1.SetY(currentEventGoodV0Temp1->Py());
          lvRotationPhoton1.SetZ(currentEventGoodV0Temp1->Pz());
          lvRotationPhoton1.SetE(currentEventGoodV0Temp1->E());

          lvRotationPhoton2.SetX(currentEventGoodV0Temp2->Px());
          lvRotationPhoton2.SetY(currentEventGoodV0Temp2->Py());
          lvRotationPhoton2.SetZ(currentEventGoodV0Temp2->Pz());
          lvRotationPhoton2.SetE(currentEventGoodV0Temp2->E());

          lvRotationPion = (lvRotationPhoton1 + lvRotationPhoton2).Vect();

          // rotate both photons around the momentum vector of their hypothetical mother particle
          if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 0 || ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 1)){
            if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 0){ // rotate by 90 degree
              rotationAngle = TMath::Pi()/2.0;
            } else if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 1){  // rotate by random angle between
               Double_t temp = (fRandom.Rndm() < 0.5) ? 0 : TMath::Pi();
               rotationAngle = temp + TMath::Pi()/3.0 + fRandom.Rndm()*TMath::Pi()/3.0;
            }
            lvRotationPhoton1.Rotate(rotationAngle, lvRotationPion);
            lvRotationPhoton2.Rotate(rotationAngle, lvRotationPion);
          }

          // Fill Eta Phi Map for Calo Photon
          cellIDRotatedPhoton = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetCaloCellIdFromEtaPhi(lvRotationPhoton1.Eta(), lvRotationPhoton1.Phi());
          if(!fDoLightOutput){
            if(!(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedPhoton, lvRotationPhoton1.Phi(), fInputEvent))){
              ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FillEtaPhiMapForClusterInBg(lvRotationPhoton1.Eta(), static_cast<double>((lvRotationPhoton1.Phi()<0) ? lvRotationPhoton1.Phi() + TMath::Pi()*2. : lvRotationPhoton1.Phi()), 1);
            }
          }

          std::unique_ptr<AliAODConversionPhoton> currentEventGoodV0Rotation1 (new AliAODConversionPhoton(&lvRotationPhoton1));
          std::unique_ptr<AliAODConversionPhoton> currentEventGoodV0Rotation2 (new AliAODConversionPhoton(&lvRotationPhoton2));

          for(auto const& kCurrentClusterCandidates  : *fClusterCandidates){
            if(currentEventGoodV0Temp1 == ((AliAODConversionPhoton*) kCurrentClusterCandidates) || currentEventGoodV0Temp2 == ((AliAODConversionPhoton*) kCurrentClusterCandidates)) continue;

            std::unique_ptr<AliAODConversionMother> backgroundCandidate2(new AliAODConversionMother(currentEventGoodV0Rotation2.get(), ((AliAODConversionPhoton*) kCurrentClusterCandidates)));

            if( fabs(currentEventGoodV0Temp2->Eta()) <= ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetEtaCut() )
            {
              if(((AliConversionMesonCuts*) fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate2.get(),kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
              {
                vSwappingInvMassPT.push_back({backgroundCandidate2->M(),backgroundCandidate2->Pt()});
              }
            }
          }
          for(auto const& kCurrentGammaCandidates  : *fGammaCandidates){
            if(currentEventGoodV0Temp1 == ((AliAODConversionPhoton*) kCurrentGammaCandidates) || currentEventGoodV0Temp2 == ((AliAODConversionPhoton*) kCurrentGammaCandidates)) continue;

            std::unique_ptr<AliAODConversionMother> backgroundCandidate1(new AliAODConversionMother(currentEventGoodV0Rotation1.get(), ((AliAODConversionPhoton*) kCurrentGammaCandidates)));
            // std::unique_ptr<AliAODConversionMother> backgroundCandidate2(new AliAODConversionMother(currentEventGoodV0Rotation2.get(), ((AliAODConversionPhoton*) kCurrentGammaCandidates)));

            if(!(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedPhoton, lvRotationPhoton1.Phi(), fInputEvent)))
            {
              if(((AliConversionMesonCuts*) fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate1.get(),kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
              {
                vSwappingInvMassPT.push_back({backgroundCandidate1->M(),backgroundCandidate1->Pt()});
              }
            }
          }
        }
      }
    }
    // Fill the histograms
    if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoWeightingInSwappBg() && vSwappingInvMassPT.size() > 0){
      tempMultWeightSwapping = (fGammaCandidates->GetEntries()*fClusterCandidates->GetEntries())/(vSwappingInvMassPT.size());
    }
    for(Int_t i = 0; i < (Int_t)vSwappingInvMassPT.size(); i++){
      fHistoMotherBackInvMassPt[fiCut]->Fill(vSwappingInvMassPT.at(i)[0], vSwappingInvMassPT.at(i)[1], tempMultWeightSwapping*fWeightJetJetMC);
      if(!fDoLightOutput || fDoECalibOutput) fHistoMotherBackInvMassECalib[fiCut]->Fill(vSwappingInvMassPT.at(i)[0], vSwappingInvMassPT.at(i)[1],tempMultWeightSwapping*fWeightJetJetMC);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::CalculateBackgroundRP(){

  Int_t zbin = 0;
  Int_t mbin = 0;

  if(fDoTHnSparse){
    zbin= fBGHandlerRP[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
    if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
      mbin = fBGHandlerRP[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
    }else {
      mbin = fBGHandlerRP[fiCut]->GetMultiplicityBinIndex(fGammaCandidates->GetEntries());
    }
  }

  //Rotation Method
  if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseRotationMethod()){
    // Correct for the number of rotations
    // BG is for rotation the same, except for factor NRotations
    Double_t weight=1./Double_t(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetNumberOfBGEvents());

    for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries();firstGammaIndex++){

      AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));
      if (gamma0==NULL) continue;
      for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fGammaCandidates->GetEntries();secondGammaIndex++){
        AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(secondGammaIndex));
        if (gamma1 == NULL) continue;
        if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelected(gamma1,fInputEvent))continue;
        for(Int_t nRandom=0;nRandom<((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetNumberOfBGEvents();nRandom++){

          RotateParticle(gamma1);
          AliAODConversionMother backgroundCandidate(gamma0,gamma1);
          backgroundCandidate.CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
          if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
            ->MesonIsSelected(&backgroundCandidate,kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
            fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate.M(),backgroundCandidate.Pt(),fWeightJetJetMC);
            if(!fDoLightOutput) fHistoPhotonPairMixedEventPtconv[fiCut]->Fill(backgroundCandidate.M(),gamma0->Pt());
            if(fDoTHnSparse){
              Double_t sparesFill[4] = {backgroundCandidate.M(),backgroundCandidate.Pt(),(Double_t)zbin,(Double_t)mbin};
              fSparseMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,weight);
            }
          }
        }
      }
    }
  }else {
    // Do Event Mixing
    for(Int_t nEventsInBG=0;nEventsInBG <fBGHandlerRP[fiCut]->GetNBGEvents(fGammaCandidates,fInputEvent);nEventsInBG++){

      AliGammaConversionPhotonVector *previousEventGammas = fBGHandlerRP[fiCut]->GetBGGoodGammas(fGammaCandidates,fInputEvent,nEventsInBG);

      if(previousEventGammas){
        // test weighted background
        Double_t weight=1.0;
        // Correct for the number of eventmixing:
        // N gammas -> (N-1) + (N-2) +(N-3) ...+ (N-(N-1))  using sum formula sum(i)=N*(N-1)/2  -> N*(N-1)/2
        // real combinations (since you cannot combine a photon with its own)
        // but BG leads to N_{a}*N_{b}combinations
        weight*=0.5*(Double_t(fGammaCandidates->GetEntries()-1))/Double_t(previousEventGammas->size());

        for(Int_t iCurrent=0;iCurrent<fGammaCandidates->GetEntries();iCurrent++){
          AliAODConversionPhoton *gamma0 = (AliAODConversionPhoton*)(fGammaCandidates->At(iCurrent));
          for(UInt_t iPrevious=0;iPrevious<previousEventGammas->size();iPrevious++){

            AliAODConversionPhoton *gamma1 = (AliAODConversionPhoton*)(previousEventGammas->at(iPrevious));

            AliAODConversionMother backgroundCandidate(gamma0,gamma1);
            backgroundCandidate.CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
            if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
            ->MesonIsSelected(&backgroundCandidate,kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
              fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate.M(),backgroundCandidate.Pt(),fWeightJetJetMC);
              if(!fDoLightOutput) fHistoPhotonPairMixedEventPtconv[fiCut]->Fill(backgroundCandidate.M(),gamma0->Pt());
              if(fDoTHnSparse){
                Double_t sparesFill[4] = {backgroundCandidate.M(),backgroundCandidate.Pt(),(Double_t)zbin,(Double_t)mbin};
                fSparseMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,weight);
              }
            }
          }
        }
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::RotateParticle(AliAODConversionPhoton *gamma){
  Int_t fNDegreesPMBackground= ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->NDegreesRotation();
  Double_t nRadiansPM = fNDegreesPMBackground*TMath::Pi()/180;
  Double_t rotationValue = fRandom.Rndm()*2*nRadiansPM + TMath::Pi()-nRadiansPM;
  gamma->RotateZ(rotationValue);
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::RotateParticleAccordingToEP(AliAODConversionPhoton *gamma, Double_t previousEventEP, Double_t thisEventEP){
  previousEventEP=previousEventEP+TMath::Pi();
  thisEventEP=thisEventEP+TMath::Pi();
  Double_t rotationValue= thisEventEP-previousEventEP;
  gamma->RotateZ(rotationValue);
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::MoveParticleAccordingToVertex(AliAODConversionPhoton* particle,const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex){
  //see header file for documentation

  Double_t dx = vertex->fX - fInputEvent->GetPrimaryVertex()->GetX();
  Double_t dy = vertex->fY - fInputEvent->GetPrimaryVertex()->GetY();
  Double_t dz = vertex->fZ - fInputEvent->GetPrimaryVertex()->GetZ();

  Double_t movedPlace[3] = {particle->GetConversionX() - dx,particle->GetConversionY() - dy,particle->GetConversionZ() - dz};
  particle->SetConversionPoint(movedPlace);
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::UpdateEventByEventData(){
  //see header file for documentation
  if(fDoJetAnalysis && fConvJetReader->GetNJets() == 0) return;
  if(fGammaCandidates->GetEntries() >0 && fClusterCandidates->GetEntries() >0 ){
    if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
      fBGHandler[fiCut]->AddEvent(fGammaCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fV0Reader->GetNumberOfPrimaryTracks(),fEventPlaneAngle);
      fBGClusHandler[fiCut]->AddEvent(fClusterCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fV0Reader->GetNumberOfPrimaryTracks(),fEventPlaneAngle);
    } else if( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoJetMixing() ){
      if(fJetNearEMCal){
        fVectorJetEta = fConvJetReader->GetVectorJetEta();
        fVectorJetPhi = fConvJetReader->GetVectorJetPhi();
        fVectorJetPt = fConvJetReader->GetVectorJetPt();
        Int_t mBinPt = 0;
        if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoJetPtMixing()){
          if(fVectorJetPt.at(fMaxPtNearEMCalPlace) > 15) mBinPt = 3;
          else mBinPt = 2;
        }else mBinPt = 2;
        fBGHandler[fiCut]->AddEvent(fClusterCandidates,fVectorJetEta.at(fMaxPtNearEMCalPlace),fVectorJetPhi.at(fMaxPtNearEMCalPlace), 2,mBinPt,fVectorJetPt.at(fMaxPtNearEMCalPlace));
      }
    }else { // means we use #V0s for multiplicity
      fBGHandler[fiCut]->AddEvent(fGammaCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fGammaCandidates->GetEntries(),fEventPlaneAngle);
      fBGClusHandler[fiCut]->AddEvent(fClusterCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fGammaCandidates->GetEntries(),fEventPlaneAngle);
    }
  }
  if(fIsMC>0 && fDoHBTHistoOutput){
    if(fTrueGammaCandidatesConv->GetEntries()>0 && fTrueGammaCandidatesCalo->GetEntries() > 0){
        fBGHBTTrueGammaHandler[fiCut]->AddMCParticleEvent(fTrueGammaCandidatesConv,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fTrueGammaCandidatesConv->GetEntries(),fEventPlaneAngle);
    }
    if(fMCGammaCandidates->GetEntries()>1){
        fBGHBTGenGammaHandler[fiCut]->AddMCParticleEvent(fMCGammaCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fMCGammaCandidates->GetEntries(),fEventPlaneAngle);
    }
  }
}


//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::FillPhotonCombinatorialBackgroundHist(AliAODConversionPhoton *TruePhotonCandidate, Int_t pdgCode[])
{
  if(fDoLightOutput) return;
  // Combinatorial Bck = 0 ee, 1 ep,i 2 ek, 3 ep, 4 emu, 5 pipi, 6 pik, 7 pip, 8 pimu, 9 kk, 10 kp, 11 kmu, 12 pp, 13 pmu, 14 mumu, 15 Rest
  if(pdgCode[0]==11   && pdgCode[1]==11){if(fIsFromDesiredHeader)fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),0.,fWeightJetJetMC);}
  else if( (pdgCode[0]==11   && pdgCode[1]==211) || (pdgCode[0]==211  && pdgCode[1]==11) )
  {if(fIsFromDesiredHeader)fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),1.,fWeightJetJetMC);}
  else if( (pdgCode[0]==11   && pdgCode[1]==321) || (pdgCode[0]==321  && pdgCode[1]==11) )
  {if(fIsFromDesiredHeader)fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),2.,fWeightJetJetMC);}
  else if( (pdgCode[0]==11   && pdgCode[1]==2212) || (pdgCode[0]==2212 && pdgCode[1]==11) )
  {if(fIsFromDesiredHeader)fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3.,fWeightJetJetMC);}
  else if( (pdgCode[0]==11   && pdgCode[1]==13) || (pdgCode[0]==13   && pdgCode[1]==11) )
  {if(fIsFromDesiredHeader)fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),4.,fWeightJetJetMC);}
  else if(  pdgCode[0]==211  && pdgCode[1]==211 ){if(fIsFromDesiredHeader)fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),5.,fWeightJetJetMC);}
  else if( (pdgCode[0]==211  && pdgCode[1]==321) || (pdgCode[0]==321  && pdgCode[1]==211) )
  {if(fIsFromDesiredHeader)fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),6.,fWeightJetJetMC);}
  else if( (pdgCode[0]==211  && pdgCode[1]==2212) || (pdgCode[0]==2212 && pdgCode[1]==211) )
  {if(fIsFromDesiredHeader)fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),7.,fWeightJetJetMC);}
  else if( (pdgCode[0]==211  && pdgCode[1]==13) || (pdgCode[0]==13   && pdgCode[1]==211) )
  {if(fIsFromDesiredHeader)fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),8.,fWeightJetJetMC);}
  else if(  pdgCode[0]==321  && pdgCode[1]==321 ){if(fIsFromDesiredHeader)fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),9.,fWeightJetJetMC);}
  else if( (pdgCode[0]==321  && pdgCode[1]==2212) || (pdgCode[0]==2212 && pdgCode[1]==321) )
  {if(fIsFromDesiredHeader)fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),10.,fWeightJetJetMC);}
  else if( (pdgCode[0]==321  && pdgCode[1]==13) || (pdgCode[0]==13   && pdgCode[1]==321) )
  {if(fIsFromDesiredHeader)fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),11.,fWeightJetJetMC);}
  else if(  pdgCode[0]==2212   && pdgCode[1]==2212  ){if(fIsFromDesiredHeader)fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),12.,fWeightJetJetMC
  );}
  else if( (pdgCode[0]==2212  && pdgCode[1]==13) || (pdgCode[0]==13   && pdgCode[1]==2212) )
  {if(fIsFromDesiredHeader)fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),13.,fWeightJetJetMC);}
  else if(  pdgCode[0]==13   && pdgCode[1]==13  ){if(fIsFromDesiredHeader)fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),14.,fWeightJetJetMC);}
  else {if(fIsFromDesiredHeader)fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),15.,fWeightJetJetMC);}
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::RelabelAODPhotonCandidates(Bool_t mode){

  // Relabeling For AOD Event
  // ESDiD -> AODiD
  // MCLabel -> AODMCLabel

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

//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::SetLogBinningXTH2(TH2* histoRebin){
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

//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::Terminate(const Option_t *)
{

  //fOutputContainer->Print(); // Will crash on GRID
}

//________________________________________________________________________
Int_t AliAnalysisTaskGammaConvCalo::GetSourceClassification(Int_t daughter, Int_t pdgCode){

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
Bool_t AliAnalysisTaskGammaConvCalo::CheckVectorOnly(vector<Int_t> &vec, Int_t tobechecked)
{
  if(tobechecked > -1)
  {
    vector<Int_t>::iterator it;
    it = find (vec.begin(), vec.end(), tobechecked);
    if (it != vec.end()) return true;
    else return false;
  }
  return false;
}

//_________________________________________________________________________________
Bool_t AliAnalysisTaskGammaConvCalo::CheckVectorForDoubleCount(vector<Int_t> &vec, Int_t tobechecked)
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

//_________________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::ProcessConversionPhotonsForMissingTags (){

  if (!fMCEvent) return;
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries();firstGammaIndex++){
    AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));
    if (gamma0->IsTrueConvertedPhoton()){
      Int_t gamma0MotherLabel = -1;
      Int_t gamma0MCLabel = gamma0->GetMCParticleLabel(fMCEvent);
      if(gamma0MCLabel > -1){
        TParticle * gammaMC0 = (TParticle*)fMCEvent->Particle(gamma0MCLabel);
        gamma0MotherLabel = gammaMC0->GetFirstMother();
        if (gamma0MotherLabel>-1){
          if(((TParticle*)fMCEvent->Particle(gamma0MotherLabel))->GetPdgCode() == 111){
            if (!CheckVectorForDoubleCount(fVectorRecTruePi0s,gamma0MotherLabel)){
              Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, gamma0MotherLabel, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
              if (!isPrimary){
                Int_t secMotherLabel = ((TParticle*)fMCEvent->Particle(gamma0MotherLabel))->GetMother(0);
                Float_t weightedSec= 1;
                if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(secMotherLabel, fMCEvent, fInputEvent) && fMCEvent->Particle(secMotherLabel)->GetPdgCode()==310){
                  weightedSec= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(secMotherLabel, fMCEvent, fInputEvent)/2.; //invariant mass is additive thus the weight for the daughters has to be devide by two for the K0s at a certain pt
                }
                if(!fDoLightOutput){
                  fHistoTrueSecondaryPi0MissingPtconvSource[fiCut]->Fill(gamma0->Pt(),0.,weightedSec*fWeightJetJetMC);
                  if(fMCEvent->Particle(secMotherLabel)->GetPdgCode()==310  ){
                    fHistoTrueSecondaryPi0MissingPtconvSource[fiCut]->Fill(gamma0->Pt(),1.,weightedSec*fWeightJetJetMC);
                  } else if(fMCEvent->Particle(secMotherLabel)->GetPdgCode()==130 ){
                    fHistoTrueSecondaryPi0MissingPtconvSource[fiCut]->Fill(gamma0->Pt(),2.,weightedSec*fWeightJetJetMC);
                  } else if(fMCEvent->Particle(secMotherLabel)->GetPdgCode()==3122 ){
                    fHistoTrueSecondaryPi0MissingPtconvSource[fiCut]->Fill(gamma0->Pt(),3.,weightedSec*fWeightJetJetMC);
                  }
                }
              }else {
                Float_t weighted= 1;
                if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma0MotherLabel, fMCEvent, fInputEvent)){
                  if (((TParticle*)fMCEvent->Particle(gamma0MotherLabel))->Pt()>0.005){
                    weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(gamma0MotherLabel, fMCEvent, fInputEvent);
                  }
                }
                if(!fDoLightOutput) fHistoTruePrimaryPi0MissingPtconv[fiCut]->Fill(gamma0->Pt(),weighted*fWeightJetJetMC);
              }
            }
          }
        }
      }
    }
  }
  return;
}

//_________________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::ProcessConversionPhotonsForMissingTagsAOD (){

  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));

  for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries();firstGammaIndex++){
    AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));

    if (gamma0->IsTrueConvertedPhoton()){
      AliAODMCParticle *positiveMC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma0->GetMCLabelPositive()));
      AliAODMCParticle *negativeMC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma0->GetMCLabelNegative()));

      Int_t gamma0MCLabel = -1;
      Int_t gamma0MotherLabel = -1;
      if(!positiveMC||!negativeMC)
        return;

      if (gamma0->IsTrueConvertedPhoton()){
        gamma0MCLabel = positiveMC->GetMother();
        AliAODMCParticle * gammaMC0 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma0MCLabel));
        gamma0MotherLabel = gammaMC0->GetMother();

        if (gamma0MotherLabel>-1){
          if(((AliAODMCParticle*)fAODMCTrackArray->At(gamma0MotherLabel))->GetPdgCode() == 111){
            if (!CheckVectorForDoubleCount(fVectorRecTruePi0s,gamma0MotherLabel)){
              Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma0MotherLabel)), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
              if (!isPrimary){
                Int_t secMotherLabel = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma0MotherLabel))->GetMother();
                Float_t weightedSec= 1;
                if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(secMotherLabel, 0x0, fInputEvent) && static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->GetPdgCode()==310){
                  weightedSec= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(secMotherLabel, 0x0, fInputEvent)/2.; //invariant mass is additive thus the weight for the daughters has to be devide by two for the K0s at a certain pt
                }
                if(!fDoLightOutput){
                  fHistoTrueSecondaryPi0MissingPtconvSource[fiCut]->Fill(gamma0->Pt(),0.,weightedSec*fWeightJetJetMC);
                  if(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->GetPdgCode()==310){
                    fHistoTrueSecondaryPi0MissingPtconvSource[fiCut]->Fill(gamma0->Pt(),1.,weightedSec*fWeightJetJetMC);
                  } else if(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->GetPdgCode()==130 ){
                    fHistoTrueSecondaryPi0MissingPtconvSource[fiCut]->Fill(gamma0->Pt(),2.,weightedSec*fWeightJetJetMC);
                  } else if(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->GetPdgCode()==3122 ){
                    fHistoTrueSecondaryPi0MissingPtconvSource[fiCut]->Fill(gamma0->Pt(),3.,weightedSec*fWeightJetJetMC);
                  }
                }
              } else {
                Float_t weighted= 1;
                if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma0MotherLabel, 0x0, fInputEvent)){
                  if (static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma0MotherLabel))->Pt()>0.005){
                    weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(gamma0MotherLabel, 0x0, fInputEvent);
                  }
                }
                if(!fDoLightOutput) fHistoTruePrimaryPi0MissingPtconv[fiCut]->Fill(gamma0->Pt(),weighted*fWeightJetJetMC);
              }
            }
          }
        }
      }
    }
  }
  return;
}

//_________________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::FillMultipleCountMap(map<Int_t,Int_t> &ma, Int_t tobechecked){
  if( ma.find(tobechecked) != ma.end() ) ma[tobechecked] += 1;
  else ma[tobechecked] = 2;
  return;
}

//_________________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::FillMultipleCountHistoAndClear(map<Int_t,Int_t> &ma, TH1F* hist){
  map<Int_t, Int_t>::iterator it;
  for (it = ma.begin(); it != ma.end(); it++){
    hist->Fill(it->second, fWeightJetJetMC);
  }
  ma.clear();
  return;
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::ProcessAODMCParticlesForHBT()
{
  if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (fAODMCTrackArray == NULL) return;
  TLorentzVector p1, p2, pTot;

  // Loop over all primary MC particle
  for(Long_t i = 0; i < fAODMCTrackArray->GetEntriesFast(); i++) {
    AliAODMCParticle* particle1 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(i));
    if (!particle1) continue;
    if (!particle1->IsPhysicalPrimary()) continue;
    if(particle1->E()<0.2) continue;
    if(TMath::Abs(particle1->Eta())>0.8) continue;
    if(particle1->GetPdgCode() == 22){
      p1.SetPxPyPzE(particle1->Px(),particle1->Py(),particle1->Pz(),particle1->E());
      fMCGammaCandidates->Add(particle1);
      for(Long_t j = i+1; j < fAODMCTrackArray->GetEntriesFast(); j++) {
        AliAODMCParticle* particle2 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(j));
        if (!particle2) continue;
        if (!particle2->IsPhysicalPrimary()) continue;
        if(particle2->E()<0.2) continue;
        if(TMath::Abs(particle2->Eta())>0.8) continue;
        if(particle2->GetPdgCode() == 22){
          p2.SetPxPyPzE(particle2->Px(),particle2->Py(),particle2->Pz(),particle2->E());
          if(p1.Angle(p2.Vect())>0.005){
            pTot = p1 + p2;
            //for sure dealing with two generated photons
            fHistoMCGammaHBTInvMassPt[fiCut]->Fill(pTot.M(),pTot.Pt());
            fHistoMCGammaHBTOpeningAnglePt[fiCut]->Fill(p1.Angle(p2.Vect()),pTot.Pt());
            fHistoMCGammaHBTDeltaEPt[fiCut]->Fill(abs(particle1->E()-particle2->E()),pTot.Pt());
          }
        }
      }
    }
  }

  if(fGammaCandidates->GetEntries()>0){
    for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries();firstGammaIndex++){
      AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));
      if (gamma0==NULL) continue;
      AliAODMCParticle *positiveMC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma0->GetMCLabelPositive()));
      AliAODMCParticle *negativeMC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma0->GetMCLabelNegative()));
      Int_t gamma0MCLabel = -1;
      if(!positiveMC||!negativeMC) continue;
      if (gamma0->IsTrueConvertedPhoton()){
        gamma0MCLabel = positiveMC->GetMother();
        AliAODMCParticle * gammaMC0 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma0MCLabel));
        if(gammaMC0->GetPdgCode() == 22) fTrueGammaCandidatesConv->Add(gammaMC0);
      }
    }
  }
  if(fClusterCandidates->GetEntries()>0){
    for(Int_t secondGammaIndex=0;secondGammaIndex<fClusterCandidates->GetEntries();secondGammaIndex++){
          AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(secondGammaIndex));
          if (gamma1==NULL) continue;
          if (gamma1->GetIsCaloPhoton() == 0) AliFatal("CaloPhotonFlag has not been set. Aborting");
          Int_t gamma1MCLabel = gamma1->GetCaloPhotonMCLabel(0); // get most probable MC label
          AliAODMCParticle * gammaMC1 = 0x0;
          if(gamma1MCLabel != -1){
            gammaMC1 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MCLabel));
            if(gammaMC1->GetPdgCode() == 22) fTrueGammaCandidatesCalo->Add(gammaMC1);
          }
    }
  }
  // Fill TrueGamma Inv mass for same event pairs
  if(fTrueGammaCandidatesConv->GetEntries()>0 && fClusterCandidates->GetEntries()>0) {
    for(Int_t firstGammaIndex=0;firstGammaIndex<fTrueGammaCandidatesConv->GetEntries();firstGammaIndex++){
      AliAODMCParticle TrueGammaConv = *(AliAODMCParticle*)(fTrueGammaCandidatesConv->At(firstGammaIndex));
      p1.SetPxPyPzE(TrueGammaConv.Px(),TrueGammaConv.Py(),TrueGammaConv.Pz(),TrueGammaConv.E());
      for(Int_t secondGammaIndex=0;secondGammaIndex<fTrueGammaCandidatesCalo->GetEntries();secondGammaIndex++){
        AliAODMCParticle TrueGammaCalo = *(AliAODMCParticle*)(fTrueGammaCandidatesCalo->At(secondGammaIndex));
        p2.SetPxPyPzE(TrueGammaCalo.Px(),TrueGammaCalo.Py(),TrueGammaCalo.Pz(),TrueGammaCalo.E());
        if(p1.Angle(p2.Vect())>0.005){
          pTot = p1 + p2;
          fHistoTrueGammaInvMassPt[fiCut]->Fill(pTot.M(),pTot.Pt());
        }
      }
    }
  }

}

//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::CalculateHBTBackgroundMC(){

  Int_t zbin;
  Int_t mbin;
  TLorentzVector p1, p2, pTot;

  // mixing MC photons -> fHistoBckTrueGammaInvMassPt
  zbin = fBGHBTTrueGammaHandler[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
  mbin = fBGHBTTrueGammaHandler[fiCut]->GetMultiplicityBinIndex(fTrueGammaCandidatesConv->GetEntries());
  for(Int_t nEventsInBG=0;nEventsInBG <fBGHBTTrueGammaHandler[fiCut]->GetNBGEvents();nEventsInBG++){
    AliAODMCParticleVector *previousEventV0s = fBGHBTTrueGammaHandler[fiCut]->GetBGGoodV0sMC(zbin,mbin,nEventsInBG);
    if(previousEventV0s){
//       cout << "previousEventV0s->size() = " << previousEventV0s->size() << endl;
      for(Int_t iCurrent=0;iCurrent<fTrueGammaCandidatesCalo->GetEntries();iCurrent++){
        AliAODMCParticle currentEventGoodV0 = *(AliAODMCParticle*)(fTrueGammaCandidatesCalo->At(iCurrent));
        p1.SetPxPyPzE(currentEventGoodV0.Px(),currentEventGoodV0.Py(),currentEventGoodV0.Pz(),currentEventGoodV0.E());
        for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){

          AliAODMCParticle previousGoodV0 = (AliAODMCParticle)(*(previousEventV0s->at(iPrevious)));
          p2.SetPxPyPzE(previousGoodV0.Px(),previousGoodV0.Py(),previousGoodV0.Pz(),previousGoodV0.E());

          if(p1.Angle(p2.Vect())>0.005){
            pTot = p1 + p2;
            fHistoBckTrueGammaInvMassPt[fiCut]->Fill(pTot.M(),pTot.Pt());
          }
        }
      }
    }
  }

  // mixing MC photons -> fHistoBckMCGammaHBTInvMassPt
  zbin = fBGHBTGenGammaHandler[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
  mbin = fBGHBTGenGammaHandler[fiCut]->GetMultiplicityBinIndex(fMCGammaCandidates->GetEntries());
  for(Int_t nEventsInBG=0;nEventsInBG <fBGHBTGenGammaHandler[fiCut]->GetNBGEvents();nEventsInBG++){
    AliAODMCParticleVector *previousEventV0sGen = fBGHBTGenGammaHandler[fiCut]->GetBGGoodV0sMC(zbin,mbin,nEventsInBG);
    if(previousEventV0sGen){
//       cout << "previousEventV0sGen->size() = " << previousEventV0sGen->size() << endl;
      for(Int_t iCurrent=0;iCurrent<fMCGammaCandidates->GetEntries();iCurrent++){
        AliAODMCParticle currentEventGoodV0 = *(AliAODMCParticle*)(fMCGammaCandidates->At(iCurrent));
        p1.SetPxPyPzE(currentEventGoodV0.Px(),currentEventGoodV0.Py(),currentEventGoodV0.Pz(),currentEventGoodV0.E());
        for(UInt_t iPrevious=0;iPrevious<previousEventV0sGen->size();iPrevious++){

          AliAODMCParticle previousGoodV0 = (AliAODMCParticle)(*(previousEventV0sGen->at(iPrevious)));
          p2.SetPxPyPzE(previousGoodV0.Px(),previousGoodV0.Py(),previousGoodV0.Pz(),previousGoodV0.E());

          if(p1.Angle(p2.Vect())>0.005){
            pTot = p1 + p2;
            fHistoBckMCGammaHBTInvMassPt[fiCut]->Fill(pTot.M(),pTot.Pt());
            fHistoBckMCGammaHBTOpeningAnglePt[fiCut]->Fill(p1.Angle(p2.Vect()),pTot.Pt());
            fHistoBckMCGammaHBTDeltaEPt[fiCut]->Fill(abs(currentEventGoodV0.E()-previousGoodV0.E()),pTot.Pt());
          }
        }
      }
    }
  }
}
