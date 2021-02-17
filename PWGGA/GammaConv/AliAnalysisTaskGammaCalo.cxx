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
#include "AliCaloTriggerMimicHelper.h"
#include <vector>
#include <map>
#include <fstream>

ClassImp(AliAnalysisTaskGammaCalo)

//________________________________________________________________________
AliAnalysisTaskGammaCalo::AliAnalysisTaskGammaCalo(): AliAnalysisTaskSE(),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fCaloTriggerHelperName(""),
  fCorrTaskSetting(""),
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
  fReaderGammas(NULL),
  fGammaCandidates(NULL),
  fClusterCandidates(NULL),
  fEventCutArray(NULL),
  fEventCuts(NULL),
  fClusterCutArray(NULL),
  fCaloPhotonCuts(NULL),
  fMesonCutArray(NULL),
  fMesonCuts(NULL),
  fConvJetReader(NULL),
  fOutlierJetReader(NULL),
  fConversionCuts(NULL),
  fDoJetAnalysis(kFALSE),
  fDoJetQA(kFALSE),
  fDoTrueSphericity(kFALSE),
  fJetHistograms(NULL),
  fTrueJetHistograms(NULL),
  fJetSector(0),
  fMaxPtNearEMCalPlace(0),
  fJetNearEMCal(kFALSE),
  fDDLRange_HistoClusGamma(NULL),
  fCaloTriggerMimicHelper(NULL),
  fSetEventCutsOutputlist(),
  fHistoMotherInvMassPt(NULL),
  fSparseMotherInvMassPtZM(NULL),
  fHistoMotherBackInvMassPt(NULL),
  fSparseMotherBackInvMassPtZM(NULL),
  fHistoMotherPi0PtY(NULL),
  fHistoMotherEtaPtY(NULL),
  fHistoMotherPi0PtAlpha(NULL),
  fHistoMotherEtaPtAlpha(NULL),
  fHistoMotherPi0PtOpenAngle(NULL),
  fHistoMotherEtaPtOpenAngle(NULL),
  fHistoMotherPtOpenAngle(NULL),
  fHistoMotherPtOpenAngleBck(NULL),
  fHistoMotherPi0NGoodESDTracksPt(NULL),
  fHistoMotherEtaNGoodESDTracksPt(NULL),
  fHistoMotherInvMassECalib(NULL),
  fHistoMotherBackInvMassECalib(NULL),
  fHistoClusGammaPt(NULL),
  fHistoClusGammaE(NULL),
  fHistoClusGammaE_BothBM(NULL),
  fHistoClusGammaE_BothBM_highestE(NULL),
  fHistoClusGammaE_AnaBM_highestE(NULL),
  fHistoClusGammaE_onlyTriggered(NULL),
  fHistoClusGammaE_DDL(NULL),
  fHistoClusGammaE_DDL_TrBM(NULL),
  fHistoClusGammaE_DDL_wL0_TrigEv(NULL),
  fHistoClusGammaE_DDL_woL0_TrigEv(NULL),
  fHistoClusGammaE_DDL_wL0_TrigEv_TrBM(NULL),
  fHistoClusGammaE_DDL_woL0_TrigEv_TrBM(NULL),
  fHistoGoodMesonClusters(NULL),
  fHistoClusOverlapHeadersGammaPt(NULL),
  fHistoClusAllHeadersGammaPt(NULL),
  fHistoClusRejectedHeadersGammaPt(NULL),
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
  fHistoTruePi0InvMassPtAdditional(NULL),
  fHistoTrueEtaInvMassPt(NULL),
  fHistoTrueEtaInvMassPtAdditional(NULL),
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
  fVectorDoubleCountTruePi0s(0),
  fVectorDoubleCountTrueEtas(0),
  fVectorDoubleCountTrueClusterGammas(0),
  fHistoMultipleCountTrueClusterGamma(NULL),
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
  fHistoNGammaCandidatesBasic(NULL),
  fHistoNGoodESDTracksVsNGammaCandidates(NULL),
  fHistoSPDClusterTrackletBackground(NULL),
  fHistoNV0Tracks(NULL),
  fProfileEtaShift(NULL),
  fProfileJetJetXSection(NULL),
  fHistoJetJetNTrials(NULL),
  fHistoPtHardJJWeight(NULL),
  fHistoEventSphericity(NULL),
  fHistoEventSphericityAxis(NULL),
  fHistoEventSphericityvsNtracks(NULL),
  fHistoEventSphericityvsNJets(NULL),
  fHistoEventMultiplicityvsNJets(NULL),
  fHistoTrueSphericityvsRecSphericity(NULL),
  fHistoTrueMultiplicityvsRecMultiplicity(NULL),
  fHistoEventSphericityvsHighpt(NULL),
  fHistoEventSphericityvsTotalpt(NULL),
  fHistoEventSphericityvsMeanpt(NULL),
  fHistoPionSpectrum(NULL),
  fHistoProtonSpectrum(NULL),
  fHistoKaonSpectrum(NULL),
  fHistoNPionSpectrum(NULL),
  fHistoEtaSpectrum(NULL),
  fHistoDMesonSpectrum(NULL),
  tTreeSphericity(NULL),
  fRecSph(0),
  fTrueSph(0),
  fPi0Pt(0),
  fPi0InvMass(0),
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
  fHistoMCPi0GenVsNClus(NULL),
  fHistoMCPi0GenFoundInOneCluster(NULL),
  fHistoMCPi0GenFoundInTwoCluster(NULL),
  fHistoMCEtaGenFoundInOneCluster(NULL),
  fHistoMCEtaGenFoundInTwoCluster(NULL),
  fHistoMCGammaConvRvsPt(NULL),
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
  fNumberOfClusters(NULL),
  fNumberOfClustersinJets(NULL),
  fEnergyRatio(NULL),
  fEnergyRatioinJets(NULL),
  fEnergyRatioGamma1(NULL),
  fEnergyRatioGamma1inJets(NULL),
  fEnergyRatioGammaAnywhere(NULL),
  fEnergyRatioGammaAnywhereinJets(NULL),
  fEnergyDeposit(NULL),
  fEnergyDepositinJets(NULL),
  fEnergyDepGamma1(NULL),
  fEnergyDepGamma1inJets(NULL),
  fEnergyDepGammaAnywhere(NULL),
  fEnergyDepGammaAnywhereinJets(NULL),
  fEnergyRatioGamma1Helped(NULL),
  fEnergyRatioGamma1HelpedinJets(NULL),
  fClusterEtaPhiJets(NULL),
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
  tTrueInvMassROpenABPtFlag(NULL),
  fInvMass(-1),
  fRconv(-1),
  fOpenRPrim(-1),
  fInvMassRTOF(-1),
  fPt(-1),
  iFlag(3),
  tSigInvMassPtAlphaTheta(NULL),
  tBckInvMassPtAlphaTheta(NULL),
  fInvMassTreeInvMass(0),
  fInvMassTreePt(0),
  fInvMassTreeAlpha(0),
  fInvMassTreeTheta(0),
  fInvMassTreeMixPool(0),
  fInvMassTreeZVertex(0),
  fInvMassTreeEta(0),
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
  tClusterTimingEff(NULL),
  fClusterTimeTag(0),
  fClusterTimeProbe(0),
  fClusterETag(0),
  fClusterEProbe(0),
//  fHistoTruePi0NonLinearity(NULL),
//  fHistoTrueEtaNonLinearity(NULL),
  fEventPlaneAngle(-100),
  fRandom(0),
  fnCuts(0),
  fiCut(0),
  fIsHeavyIon(0),
  fDoLightOutput(kFALSE),
  fDoECalibOutput(kFALSE),
  fDoMesonAnalysis(kTRUE),
  fDoMesonQA(0),
  fDoClusterQA(0),
  fIsFromDesiredHeader(kTRUE),
  fIsOverlappingWithOtherHeader(kFALSE),
  fIsMC(0),
  fDoTHnSparse(kTRUE),
  fSetPlotHistsExtQA(kFALSE),
  fDoSoftAnalysis(kFALSE),
  fWeightJetJetMC(1),
  fDoInOutTimingCluster(kFALSE),
  fMinTimingCluster(0),
  fMaxTimingCluster(0),
  fEnableSortForClusMC(kFALSE),
  fProduceCellIDPlots(kFALSE),
  fProduceTreeEOverP(kFALSE),
  tBrokenFiles(NULL),
  fFileNameBroken(NULL),
  tTriggerFiles_wL0(NULL),
  tTriggerFiles_woL0(NULL),
  fFileNameTrigger(NULL),
  tClusterQATree(NULL),
  fCloseHighPtClusters(NULL),
  fGenPhaseSpace(),
  fAODMCTrackArray(NULL),
  fLocalDebugFlag(0),
  fAllowOverlapHeaders(kTRUE),
  fNCurrentClusterBasic(0),
  fTrackMatcherRunningMode(0),
  fDoPi0Only(kFALSE)
{

}

//________________________________________________________________________
AliAnalysisTaskGammaCalo::AliAnalysisTaskGammaCalo(const char *name):
  AliAnalysisTaskSE(name),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fCaloTriggerHelperName(""),
  fCorrTaskSetting(""),
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
  fOutputContainer(0),
  fReaderGammas(NULL),
  fGammaCandidates(NULL),
  fClusterCandidates(NULL),
  fEventCutArray(NULL),
  fEventCuts(NULL),
  fClusterCutArray(NULL),
  fCaloPhotonCuts(NULL),
  fMesonCutArray(NULL),
  fMesonCuts(NULL),
  fConvJetReader(NULL),
  fOutlierJetReader(NULL),
  fConversionCuts(NULL),
  fDoJetAnalysis(kFALSE),
  fDoJetQA(kFALSE),
  fDoTrueSphericity(kFALSE),
  fJetHistograms(NULL),
  fTrueJetHistograms(NULL),
  fJetSector(0),
  fMaxPtNearEMCalPlace(0),
  fJetNearEMCal(kFALSE),
  fDDLRange_HistoClusGamma(NULL),
  fCaloTriggerMimicHelper(NULL),
  fSetEventCutsOutputlist(),
  fHistoMotherInvMassPt(NULL),
  fSparseMotherInvMassPtZM(NULL),
  fHistoMotherBackInvMassPt(NULL),
  fSparseMotherBackInvMassPtZM(NULL),
  fHistoMotherPi0PtY(NULL),
  fHistoMotherEtaPtY(NULL),
  fHistoMotherPi0PtAlpha(NULL),
  fHistoMotherEtaPtAlpha(NULL),
  fHistoMotherPi0PtOpenAngle(NULL),
  fHistoMotherEtaPtOpenAngle(NULL),
  fHistoMotherPtOpenAngle(NULL),
  fHistoMotherPtOpenAngleBck(NULL),
  fHistoMotherPi0NGoodESDTracksPt(NULL),
  fHistoMotherEtaNGoodESDTracksPt(NULL),
  fHistoMotherInvMassECalib(NULL),
  fHistoMotherBackInvMassECalib(NULL),
  fHistoClusGammaPt(NULL),
  fHistoClusGammaE(NULL),
  fHistoClusGammaE_BothBM(NULL),
  fHistoClusGammaE_BothBM_highestE(NULL),
  fHistoClusGammaE_AnaBM_highestE(NULL),
  fHistoClusGammaE_onlyTriggered(NULL),
  fHistoClusGammaE_DDL(NULL),
  fHistoClusGammaE_DDL_TrBM(NULL),
  fHistoClusGammaE_DDL_wL0_TrigEv(NULL),
  fHistoClusGammaE_DDL_woL0_TrigEv(NULL),
  fHistoClusGammaE_DDL_wL0_TrigEv_TrBM(NULL),
  fHistoClusGammaE_DDL_woL0_TrigEv_TrBM(NULL),
  fHistoGoodMesonClusters(NULL),
  fHistoClusOverlapHeadersGammaPt(NULL),
  fHistoClusAllHeadersGammaPt(NULL),
  fHistoClusRejectedHeadersGammaPt(NULL),
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
  fHistoTruePi0InvMassPtAdditional(NULL),
  fHistoTrueEtaInvMassPt(NULL),
  fHistoTrueEtaInvMassPtAdditional(NULL),
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
  fVectorDoubleCountTruePi0s(0),
  fVectorDoubleCountTrueEtas(0),
  fVectorDoubleCountTrueClusterGammas(0),
  fHistoMultipleCountTrueClusterGamma(NULL),
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
  fHistoNGammaCandidatesBasic(NULL),
  fHistoNGoodESDTracksVsNGammaCandidates(NULL),
  fHistoSPDClusterTrackletBackground(NULL),
  fHistoNV0Tracks(NULL),
  fProfileEtaShift(NULL),
  fProfileJetJetXSection(NULL),
  fHistoJetJetNTrials(NULL),
  fHistoPtHardJJWeight(NULL),
  fHistoEventSphericity(NULL),
  fHistoEventSphericityAxis(NULL),
  fHistoEventSphericityvsNtracks(NULL),
  fHistoEventSphericityvsNJets(NULL),
  fHistoEventMultiplicityvsNJets(NULL),
  fHistoTrueSphericityvsRecSphericity(NULL),
  fHistoTrueMultiplicityvsRecMultiplicity(NULL),
  fHistoEventSphericityvsHighpt(NULL),
  fHistoEventSphericityvsTotalpt(NULL),
  fHistoEventSphericityvsMeanpt(NULL),
  fHistoPionSpectrum(NULL),
  fHistoProtonSpectrum(NULL),
  fHistoKaonSpectrum(NULL),
  fHistoNPionSpectrum(NULL),
  fHistoEtaSpectrum(NULL),
  fHistoDMesonSpectrum(NULL),
  tTreeSphericity(NULL),
  fRecSph(0),
  fTrueSph(0),
  fPi0Pt(0),
  fPi0InvMass(0),
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
  fHistoMCPi0GenVsNClus(NULL),
  fHistoMCPi0GenFoundInOneCluster(NULL),
  fHistoMCPi0GenFoundInTwoCluster(NULL),
  fHistoMCEtaGenFoundInOneCluster(NULL),
  fHistoMCEtaGenFoundInTwoCluster(NULL),
  fHistoMCGammaConvRvsPt(NULL),
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
  fNumberOfClusters(NULL),
  fNumberOfClustersinJets(NULL),
  fEnergyRatio(NULL),
  fEnergyRatioinJets(NULL),
  fEnergyRatioGamma1(NULL),
  fEnergyRatioGamma1inJets(NULL),
  fEnergyRatioGammaAnywhere(NULL),
  fEnergyRatioGammaAnywhereinJets(NULL),
  fEnergyDeposit(NULL),
  fEnergyDepositinJets(NULL),
  fEnergyDepGamma1(NULL),
  fEnergyDepGamma1inJets(NULL),
  fEnergyDepGammaAnywhere(NULL),
  fEnergyDepGammaAnywhereinJets(NULL),
  fEnergyRatioGamma1Helped(NULL),
  fEnergyRatioGamma1HelpedinJets(NULL),
  fClusterEtaPhiJets(NULL),
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
  tTrueInvMassROpenABPtFlag(NULL),
  fInvMass(-1),
  fRconv(-1),
  fOpenRPrim(-1),
  fInvMassRTOF(-1),
  fPt(-1),
  iFlag(3),
  tSigInvMassPtAlphaTheta(NULL),
  tBckInvMassPtAlphaTheta(NULL),
  fInvMassTreeInvMass(0),
  fInvMassTreePt(0),
  fInvMassTreeAlpha(0),
  fInvMassTreeTheta(0),
  fInvMassTreeMixPool(0),
  fInvMassTreeZVertex(0),
  fInvMassTreeEta(0),
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
  tClusterTimingEff(NULL),
  fClusterTimeTag(0),
  fClusterTimeProbe(0),
  fClusterETag(0),
  fClusterEProbe(0),
//  fHistoTruePi0NonLinearity(NULL),
//  fHistoTrueEtaNonLinearity(NULL),
  fEventPlaneAngle(-100),
  fRandom(0),
  fnCuts(0),
  fiCut(0),
  fIsHeavyIon(0),
  fDoLightOutput(kFALSE),
  fDoECalibOutput(kFALSE),
  fDoMesonAnalysis(kTRUE),
  fDoMesonQA(0),
  fDoClusterQA(0),
  fIsFromDesiredHeader(kTRUE),
  fIsOverlappingWithOtherHeader(kFALSE),
  fIsMC(0),
  fDoTHnSparse(kTRUE),
  fSetPlotHistsExtQA(kFALSE),
  fDoSoftAnalysis(kFALSE),
  fWeightJetJetMC(1),
  fDoInOutTimingCluster(kFALSE),
  fMinTimingCluster(0),
  fMaxTimingCluster(0),
  fEnableSortForClusMC(kFALSE),
  fProduceCellIDPlots(kFALSE),
  fProduceTreeEOverP(kFALSE),
  tBrokenFiles(NULL),
  fFileNameBroken(NULL),
  tTriggerFiles_wL0(NULL),
  tTriggerFiles_woL0(NULL),
  fFileNameTrigger(NULL),
  tClusterQATree(NULL),
  fCloseHighPtClusters(NULL),
  fGenPhaseSpace(),
  fAODMCTrackArray(NULL),
  fLocalDebugFlag(0),
  fAllowOverlapHeaders(kTRUE),
  fNCurrentClusterBasic(0),
  fTrackMatcherRunningMode(0),
  fDoPi0Only(kFALSE)
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
  DefineOutput(9, TTree::Class());}

AliAnalysisTaskGammaCalo::~AliAnalysisTaskGammaCalo()
{
  if(fClusterCandidates){
    delete fClusterCandidates;
    fClusterCandidates = 0x0;
  }
  if(fGammaCandidates){
    delete fGammaCandidates;
    fGammaCandidates = 0x0;
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

        fSparseMotherBackInvMassPtZM[iCut] = new THnSparseF("Back_Back_InvMass_Pt_z_m", "Back_Back_InvMass_Pt_z_m", nDim,nBins,xMin,xMax);
        fBackList[iCut]->Add(fSparseMotherBackInvMassPtZM[iCut]);

        fMotherList[iCut] = new TList();
        fMotherList[iCut]->SetName(Form("%s_%s_%s Mother histograms", cutstringEvent.Data(), cutstringCalo.Data(), cutstringMeson.Data()));
        fMotherList[iCut]->SetOwner(kTRUE);
        fCutFolder[iCut]->Add(fMotherList[iCut]);

        fSparseMotherInvMassPtZM[iCut] = new THnSparseF("Back_Mother_InvMass_Pt_z_m", "Back_Mother_InvMass_Pt_z_m", nDim,nBins,xMin,xMax);
        fMotherList[iCut]->Add(fSparseMotherInvMassPtZM[iCut]);
      }

      if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->BackgroundHandlerType() == 0){
        if( ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->DoSectorMixing() || ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->DoSectorJetMixing()){
          fBGHandler[iCut] = new AliGammaConversionAODBGHandler(
                                    collisionSystem,centMin,centMax,
                                    ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetNumberOfBGEvents(),
                                    ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseTrackMultiplicity(),
                                    4,9,7);
        } else {
          fBGHandler[iCut] = new AliGammaConversionAODBGHandler(
                                    collisionSystem,centMin,centMax,
                                    ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetNumberOfBGEvents(),
                                    ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseTrackMultiplicity(),
                                    4,8,7);
        }
      }
    }
  }
}
//________________________________________________________________________
void AliAnalysisTaskGammaCalo::UserCreateOutputObjects(){

  fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
  if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader

  if(fDoMesonAnalysis){ //Same Jet Finder MUST be used within same trainconfig
    if( ((AliConversionMesonCuts*)fMesonCutArray->At(0))->DoJetAnalysis())  fDoJetAnalysis = kTRUE;
    if( ((AliConversionMesonCuts*)fMesonCutArray->At(0))->DoJetQA())        fDoJetQA       = kTRUE;
    if(((AliConvEventCuts*)fEventCutArray->At(0))->GetUseSphericityTrue()) fDoTrueSphericity = kTRUE;
  }

  if(fDoJetAnalysis){
    fConvJetReader=(AliAnalysisTaskConvJet*)AliAnalysisManager::GetAnalysisManager()->GetTask("AliAnalysisTaskConvJet");
    if(!fConvJetReader){printf("Error: No AliAnalysisTaskConvJet");return;} // GetV0Reader
  }
  if(((AliConvEventCuts*)fEventCutArray->At(0))->GetUseJetFinderForOutliers()){
    fOutlierJetReader=(AliAnalysisTaskJetOutlierRemoval*)AliAnalysisManager::GetAnalysisManager()->GetTask("AliAnalysisTaskJetOutlierRemoval");
    if(!fOutlierJetReader){AliFatal("Error: No AliAnalysisTaskJetOutlierRemoval");} // GetV0Reader
    else{printf("Found AliAnalysisTaskJetOutlierRemoval used for outlier removal!\n");}
  }
  if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetDoSecondaryTrackMatching()){
    fConversionCuts = new AliConversionPhotonCuts();
    fConversionCuts->SetV0ReaderName(fV0ReaderName.Data());
    fConversionCuts->InitializeCutsFromCutString("00200009327000008250400000"); //Use standard cuts
    fConversionCuts->SetIsHeavyIon(fIsHeavyIon);
    fConversionCuts->SetFillCutHistograms("",kFALSE);
  }

  if (fDoClusterQA == 2) fProduceCellIDPlots = kTRUE;
  if (fIsMC == 2){
    fDoTHnSparse      = kFALSE;
  } else if (fIsMC == 3){
    fDoTHnSparse      = kFALSE;
  }

  // set common binning in pT for mesons and photons
  Float_t binWidthPt          = 0.1;
  Int_t nBinsPt               = 250;
  Float_t minPt               = 0;
  Float_t maxPt               = 25;
  Int_t nBinsQAPt             = 175;
  Float_t maxQAPt             = 25;
  Int_t nBinsClusterPt        = 500;
  Float_t minClusterPt        = 0;
  Float_t maxClusterPt        = 50;
  Int_t nBinsMinv             = 800;
  Float_t maxMinv             = 0.8;
  Double_t *arrPtBinning      = new Double_t[1200];
  Double_t *arrQAPtBinning    = new Double_t[1200];
  Double_t *arrClusPtBinning  = new Double_t[1200];
  if( fDoPi0Only ){
    nBinsMinv                 = 150;
    maxMinv                  = 0.3;
  }

  if (((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::kpPb5TeVR2 ||
      ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::kpPb5TeV   ||
      ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::k5TeV ||
      ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::k8TeV ||
      ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::k7TeV ||
      ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::kpPb8TeV ){
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
  } else if ( ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::k13TeV ||
              ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::k13TeVLowB ){
    if( fDoPi0Only ){
      nBinsPt                   = 218;
      minPt                     = 0;
      maxPt                     = 25;
      binWidthPt                = 0.1;
      for(Int_t i=0; i<nBinsPt+1;i++){
        if (i < 1) arrPtBinning[i]              = 0.3*i;
        else if(i<198) arrPtBinning[i]          = 0.3+0.1*(i-1);
        else if(i<218) arrPtBinning[i]          = 20.+0.25*(i-198);
        else  arrPtBinning[i]                   = maxPt;
      }
      nBinsQAPt                 = 145;
      maxQAPt                   = 25;
      for(Int_t i=0; i<nBinsQAPt+1;i++){
        if(i<100) arrQAPtBinning[i]             = 0.1*i;
        else if(i<140) arrQAPtBinning[i]        = 10.+0.25*(i-100);
        else if(i<145) arrQAPtBinning[i]        = 20.+1.0*(i-140);
        else arrQAPtBinning[i]                  = maxQAPt;
      }
      nBinsClusterPt            = 157;
      minClusterPt              = 0;
      maxClusterPt              = 100;
      for(Int_t i=0; i<nBinsClusterPt+1;i++){
        if (i < 1) arrClusPtBinning[i]          = 0.3*i;
        else if(i<98) arrClusPtBinning[i]       = 0.3+0.1*(i-1);
        else if(i<138) arrClusPtBinning[i]      = 10.+0.25*(i-98);
        else if(i<148) arrClusPtBinning[i]      = 20.+1.0*(i-138);
        else if(i<152) arrClusPtBinning[i]      = 30.+5.0*(i-148);
        else if(i<157) arrClusPtBinning[i]      = 50.+10.*(i-152);
        else  arrClusPtBinning[i]               = maxClusterPt;
      }

    } else {
      nBinsMinv                 = 800;

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
        else  arrClusPtBinning[i]               = maxClusterPt;
      }
    }
  } else if (((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::kXeXe5440GeV  ){
    nBinsPt                   = 88;
    minPt                     = 0;
    maxPt                     = 20;
    for(Int_t i=0; i<nBinsPt+1;i++){
      if (i < 1) arrPtBinning[i]              = 0.5*i;
      else if(i<56) arrPtBinning[i]           = 0.5+0.1*(i-1);
      else if(i<80) arrPtBinning[i]           = 6.+0.25*(i-56);
      else if(i<88) arrPtBinning[i]           = 12.+1.0*(i-80);
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
    nBinsPt                   = 108;
    minPt                     = 0;
    maxPt                     = 40;
    for(Int_t i=0; i<nBinsPt+1;i++){
      if (i < 1) arrPtBinning[i]              = 0.5*i;
      else if(i<56) arrPtBinning[i]           = 0.5+0.1*(i-1);
      else if(i<80) arrPtBinning[i]           = 6.+0.25*(i-56);
      else if(i<108) arrPtBinning[i]          = 12.+1.0*(i-80);
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
      else if(i<175) arrQAPtBinning[i]        = 20.+1.0*(i-170);
      else arrQAPtBinning[i]                  = maxQAPt;
    }
  }
 Double_t* arrLogBinning = new Double_t[51]{1.00e-03, 1.20e-03, 1.45e-03, 1.74e-03, 2.01e-03, 2.51e-03, 3.02e-03, 3.63e-03, 4.37e-03, 5.25e-03, 6.31e-03, 7.60e-03,
                                     9.12e-03, 1.01e-02, 1.32e-02, 1.58e-02, 1.91e-02, 2.29e-02, 2.75e-02, 3.31e-02, 3.98e-02, 4.79e-02, 5.75e-02, 6.91e-02,
                                     8.32e-02, 1.00e-01, 1.20e-01, 1.45e-01, 1.74e-01, 2.09e-01, 2.51e-01, 3.02e-01, 3.63e-01, 4.37e-01, 5.25e-01, 6.31e-01,
                                     7.59e-01, 9.12e-01, 1.10e+00, 1.32e+00, 1.58e+00, 1.91e+00, 2.29e+00, 2.75e+00, 3.31e+00, 3.98e+00, 4.79e+00, 5.75e+00,
                                     6.92e+00, 8.32e+00, 1.00e+01};
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
  if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetDoSecondaryTrackMatching()){
    fGammaCandidates    = new TList();
  }
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
    fHistoPtHardJJWeight     = new TH2F*[fnCuts];
  }

  Bool_t EnableSphericity = kFALSE;
  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetUseSphericity()!=0){
      EnableSphericity = kTRUE;
    }
  }

  fHistoNGoodESDTracks        = new TH1F*[fnCuts];
  fHistoVertexZ               = new TH1F*[fnCuts];
  fHistoNGammaCandidates      = new TH1F*[fnCuts];
  fHistoNGammaCandidatesBasic = new TH1F*[fnCuts];
  if(EnableSphericity){
    fHistoEventSphericity       = new TH1F*[fnCuts];
    fHistoEventSphericityAxis   = new TH1F*[fnCuts];
    fHistoEventSphericityvsNtracks    = new TH2F*[fnCuts];
    fHistoEventSphericityvsNJets      = new TH2F*[fnCuts];
    fHistoEventMultiplicityvsNJets      = new TH2F*[fnCuts];
    if(fIsMC >0){
      fHistoTrueSphericityvsRecSphericity     = new TH2F*[fnCuts];
      fHistoTrueMultiplicityvsRecMultiplicity = new TH2F*[fnCuts];
      if(fDoTrueSphericity){
          tTreeSphericity                         = new TTree*[fnCuts];
      }
      fHistoPionSpectrum                      = new TH1F*[fnCuts];
      fHistoProtonSpectrum                    = new TH1F*[fnCuts];
      fHistoKaonSpectrum                      = new TH1F*[fnCuts];
      fHistoNPionSpectrum                     = new TH1F*[fnCuts];
      fHistoEtaSpectrum                       = new TH1F*[fnCuts];
      fHistoDMesonSpectrum                    = new TH1F*[fnCuts];
    }
    fHistoEventSphericityvsHighpt           = new TH2F*[fnCuts];
    fHistoEventSphericityvsTotalpt          = new TH2F*[fnCuts];
    fHistoEventSphericityvsMeanpt           = new TH2F*[fnCuts];
  }
  if(!fDoLightOutput){
    fHistoNGoodESDTracksVsNGammaCandidates  = new TH2F*[fnCuts];
    fHistoSPDClusterTrackletBackground      = new TH2F*[fnCuts];
    fHistoNV0Tracks                         = new TH1F*[fnCuts];
  }
  if(fIsHeavyIon==2) fProfileEtaShift          = new TProfile*[fnCuts];

  if(fDoMesonAnalysis){
    fHistoMotherInvMassPt           = new TH2F*[fnCuts];
    fHistoMotherBackInvMassPt       = new TH2F*[fnCuts];
    if(!fDoLightOutput || fDoPi0Only || fDoECalibOutput){
      fHistoMotherInvMassECalib         = new TH2F*[fnCuts];
      fHistoMotherBackInvMassECalib     = new TH2F*[fnCuts];
    }
    if (fDoMesonQA > 0 && fDoMesonQA < 3){
      fHistoMotherPi0PtY            = new TH2F*[fnCuts];
      fHistoMotherPi0PtAlpha        = new TH2F*[fnCuts];
      fHistoMotherPi0PtOpenAngle    = new TH2F*[fnCuts];
      fHistoMotherPi0NGoodESDTracksPt    = new TH2F*[fnCuts];
      if( !fDoPi0Only ){
        fHistoMotherEtaPtY            = new TH2F*[fnCuts];
        fHistoMotherEtaPtAlpha        = new TH2F*[fnCuts];
        fHistoMotherEtaPtOpenAngle    = new TH2F*[fnCuts];
        fHistoMotherEtaNGoodESDTracksPt    = new TH2F*[fnCuts];
      }
    }
    if (fDoMesonQA == 2){
      fHistoMotherPtOpenAngle       = new TH2F*[fnCuts];
      fHistoMotherPtOpenAngleBck    = new TH2F*[fnCuts];
    }
  }

  if(fProduceCellIDPlots){
    fHistCellIDvsClusterEnergy      = new TH2F*[fnCuts];
    fHistCellIDvsClusterEnergyMax   = new TH2F*[fnCuts];
  }

  fHistoClusGammaPt                 = new TH1F*[fnCuts];
  fHistoClusGammaE                  = new TH1F*[fnCuts];

  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetClusterType()==2){
      if (fHistoClusGammaE_BothBM_highestE==NULL){
        fHistoClusGammaE_BothBM_highestE = new TH1F*[fnCuts];
      }
      //--------------------------------------------------
      //if ( (((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger()==1)||(((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger()==0) ){ //Only MB
        if (fHistoClusGammaE_BothBM==NULL){
          fHistoClusGammaE_BothBM = new TH1F*[fnCuts];
        }
        if (fHistoClusGammaE_AnaBM_highestE==NULL){
          fHistoClusGammaE_AnaBM_highestE = new TH1F*[fnCuts];
        }
      //}
      //--------------------------------------------------
      //Only PHI7
      if ( ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger()==6 ){//Only PHI7 Trigger
        if (fHistoClusGammaE_onlyTriggered==NULL) {
          fHistoClusGammaE_onlyTriggered = new TH1F*[fnCuts];
        }
      }
    }
  }

  fHistoClusOverlapHeadersGammaPt   = new TH1F*[fnCuts];
  fHistoClusAllHeadersGammaPt       = new TH1F*[fnCuts];
  fHistoClusRejectedHeadersGammaPt  = new TH1F*[fnCuts];
  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType()==2){
        if (fCaloTriggerMimicHelper == NULL){
          fCaloTriggerMimicHelper   = new AliCaloTriggerMimicHelper*[fnCuts];
        }
        if (fHistoGoodMesonClusters==NULL){
          fHistoGoodMesonClusters   = new TH1I*[fnCuts];
        }
    }
  }
  if(!fDoLightOutput && fDoClusterQA > 0){
    for(Int_t iCut = 0; iCut<fnCuts;iCut++){
      if(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetClusterType()==2){
        if (fDDLRange_HistoClusGamma==NULL){
          fDDLRange_HistoClusGamma  = new Int_t [2];
          fDDLRange_HistoClusGamma[0] = 6;
          fDDLRange_HistoClusGamma[1] = 19;
        }
        if (fIsMC==0){
          if (fHistoClusGammaE_DDL==NULL){
            fHistoClusGammaE_DDL      = new TH1F**[fDDLRange_HistoClusGamma[1]-fDDLRange_HistoClusGamma[0]+1];
            for (Int_t DDLRange=0; DDLRange<=fDDLRange_HistoClusGamma[1]-fDDLRange_HistoClusGamma[0]; DDLRange++){
              fHistoClusGammaE_DDL[DDLRange]  = new TH1F*[fnCuts];
            }
          }
          //--------------------------------------------------
          //Only MB
          if ( (((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger()==1)||(((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger()==0) ){
            if (fHistoClusGammaE_DDL_TrBM==NULL){
              fHistoClusGammaE_DDL_TrBM             = new TH1F**[fDDLRange_HistoClusGamma[1]-fDDLRange_HistoClusGamma[0]+1];
              for (Int_t DDLRange=0; DDLRange<=fDDLRange_HistoClusGamma[1]-fDDLRange_HistoClusGamma[0]; DDLRange++){
                fHistoClusGammaE_DDL_TrBM[DDLRange]  = new TH1F*[fnCuts];
              }
            } if (fHistoClusGammaE_DDL_wL0_TrigEv==NULL){
              fHistoClusGammaE_DDL_wL0_TrigEv       = new TH1F**[fDDLRange_HistoClusGamma[1]-fDDLRange_HistoClusGamma[0]+1];
              for (Int_t DDLRange=0; DDLRange<=fDDLRange_HistoClusGamma[1]-fDDLRange_HistoClusGamma[0]; DDLRange++){
                fHistoClusGammaE_DDL_wL0_TrigEv[DDLRange]  = new TH1F*[fnCuts];
              }
            } if (fHistoClusGammaE_DDL_woL0_TrigEv==NULL){
              fHistoClusGammaE_DDL_woL0_TrigEv      = new TH1F**[fDDLRange_HistoClusGamma[1]-fDDLRange_HistoClusGamma[0]+1];
              for (Int_t DDLRange=0; DDLRange<=fDDLRange_HistoClusGamma[1]-fDDLRange_HistoClusGamma[0]; DDLRange++){
                fHistoClusGammaE_DDL_woL0_TrigEv[DDLRange]  = new TH1F*[fnCuts];
              }
            } if (fHistoClusGammaE_DDL_wL0_TrigEv_TrBM==NULL){
              fHistoClusGammaE_DDL_wL0_TrigEv_TrBM  = new TH1F**[fDDLRange_HistoClusGamma[1]-fDDLRange_HistoClusGamma[0]+1];
              for (Int_t DDLRange=0; DDLRange<=fDDLRange_HistoClusGamma[1]-fDDLRange_HistoClusGamma[0]; DDLRange++){
                fHistoClusGammaE_DDL_wL0_TrigEv_TrBM[DDLRange]  = new TH1F*[fnCuts];
              }
            } if (fHistoClusGammaE_DDL_woL0_TrigEv_TrBM==NULL){
              fHistoClusGammaE_DDL_woL0_TrigEv_TrBM = new TH1F**[fDDLRange_HistoClusGamma[1]-fDDLRange_HistoClusGamma[0]+1];
              for (Int_t DDLRange=0; DDLRange<=fDDLRange_HistoClusGamma[1]-fDDLRange_HistoClusGamma[0]; DDLRange++){
                fHistoClusGammaE_DDL_woL0_TrigEv_TrBM[DDLRange]  = new TH1F*[fnCuts];
              }
            }//fHistoClusGammaE_DDL_woL0_TrigEv_TrBM==NULL
          } //Only MB
        }//fIsMC==0
      } //ClusterType==2
    }//Loop over Cuts
    fHistoClusGammaPtM02            = new TH2F*[fnCuts];
  }

  if (fDoMesonQA == 4 && fIsMC == 0){
    fTreeList                       = new TList*[fnCuts];
    tSigInvMassPtAlphaTheta         = new TTree*[fnCuts];
    tBckInvMassPtAlphaTheta         = new TTree*[fnCuts];
  }

  if (fDoMesonQA == 5 && fIsMC == 0){
    tClusterTimingEff               = new TTree*[fnCuts];
  }

  if (fProduceTreeEOverP){
    fClusterTreeList                = new TList*[fnCuts];
    tClusterEOverP                  = new TTree*[fnCuts];
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

  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    TString cutstringEvent    = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
    TString cutstringCalo     = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
    TString cutstringMeson    = "NoMesonCut";
    if(fDoMesonAnalysis)
      cutstringMeson          = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();

    fCutFolder[iCut]        = new TList();
    fCutFolder[iCut]->SetName(Form("Cut Number %s_%s_%s", cutstringEvent.Data(), cutstringCalo.Data(), cutstringMeson.Data()));
    fCutFolder[iCut]->SetOwner(kTRUE);
    fOutputContainer->Add(fCutFolder[iCut]);
    fESDList[iCut]          = new TList();
    fESDList[iCut]->SetName(Form("%s_%s_%s ESD histograms", cutstringEvent.Data(), cutstringCalo.Data(), cutstringMeson.Data()));
    fESDList[iCut]->SetOwner(kTRUE);
    fCutFolder[iCut]->Add(fESDList[iCut]);

    fHistoNEvents[iCut]     = new TH1F("NEvents", "NEvents", 15, -0.5, 14.5);
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
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(10,"EMCAL/TPC problems");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(11,"rejectedForJetJetMC");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(14,"Pileup V0M-TPCout Tracks");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(15,"Sphericity");
    fESDList[iCut]->Add(fHistoNEvents[iCut]);

    if (fIsMC > 1){
      fHistoNEventsWOWeight[iCut] = new TH1F("NEventsWOWeight", "NEventsWOWeight", 15, -0.5, 14.5);
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
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(15,"Sphericity");
      fESDList[iCut]->Add(fHistoNEventsWOWeight[iCut]);
    }
    if (fIsMC == 2){
      fProfileJetJetXSection[iCut]  = new TProfile("XSection", "XSection", 1, -0.5, 0.5);
      fESDList[iCut]->Add(fProfileJetJetXSection[iCut]);
      fHistoJetJetNTrials[iCut]     = new TH1F("NTrials", "#sum{NTrials}", 1, 0, 1);
      fHistoJetJetNTrials[iCut]->GetXaxis()->SetBinLabel(1,"#sum{NTrials}");
      fESDList[iCut]->Add(fHistoJetJetNTrials[iCut]);
      fHistoPtHardJJWeight[iCut]     = new TH2F("fHistoPtHardJJWeight", "fHistoPtHardJJWeight", 400, 0, 200, 60, 0, 30);
      fESDList[iCut]->Add(fHistoPtHardJJWeight[iCut]);
    }

    if(fIsHeavyIon == 1)
      fHistoNGoodESDTracks[iCut]    = new TH1F("GoodESDTracks", "GoodESDTracks", 4000, 0, 4000);
    else if(fIsHeavyIon == 2)
      fHistoNGoodESDTracks[iCut]    = new TH1F("GoodESDTracks", "GoodESDTracks", 400, 0, 400);
    else
      fHistoNGoodESDTracks[iCut]    = new TH1F("GoodESDTracks", "GoodESDTracks", 200, 0, 200);
    fHistoNGoodESDTracks[iCut]->GetXaxis()->SetTitle("#primary tracks");
    fESDList[iCut]->Add(fHistoNGoodESDTracks[iCut]);

    if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetUseSphericity()!=0){
      fHistoEventSphericity[iCut]     = new TH1F("EventSphericity", "EventSphericity", 100, 0, 1);
      fHistoEventSphericity[iCut]->GetXaxis()->SetTitle("S");
      fESDList[iCut]->Add(fHistoEventSphericity[iCut]);
      fV0Reader->SetCalcSphericity(kTRUE);
      fHistoEventSphericityAxis[iCut] = new TH1F("EventSphericityAxisPhi", "EventSphericityAxisPhi", 200, 0, 2*TMath::Pi());
      fHistoEventSphericityAxis[iCut]->GetXaxis()->SetTitle("Phi");
      fESDList[iCut]->Add(fHistoEventSphericityAxis[iCut]);
      fHistoEventSphericityvsNtracks[iCut]  = new TH2F("EventSphericity vs Ntracks", "EventSphericity vs Ntracks", 100, 0, 1, 100, 0, 100);
      fHistoEventSphericityvsNtracks[iCut]->GetXaxis()->SetTitle("S");
      fHistoEventSphericityvsNtracks[iCut]->GetYaxis()->SetTitle("Ntracks");
      fESDList[iCut]->Add(fHistoEventSphericityvsNtracks[iCut]);
      if(fDoJetAnalysis){
      fHistoEventSphericityvsNJets[iCut]  = new TH2F("EventSphericity vs NJets", "EventSphericity vs NJets", 100, 0, 1, 10, 0, 10);
      fHistoEventSphericityvsNJets[iCut]->GetXaxis()->SetTitle("S");
      fHistoEventSphericityvsNJets[iCut]->GetYaxis()->SetTitle("NJets");
      fESDList[iCut]->Add(fHistoEventSphericityvsNJets[iCut]);
      fHistoEventMultiplicityvsNJets[iCut]  = new TH2F("Multiplicity vs NJets", "EventSphericity vs NJets", 100, 0, 100, 10, 0, 10);
      fHistoEventMultiplicityvsNJets[iCut]->GetXaxis()->SetTitle("Mult");
      fHistoEventMultiplicityvsNJets[iCut]->GetYaxis()->SetTitle("NJets");
      fESDList[iCut]->Add(fHistoEventMultiplicityvsNJets[iCut]);
      }
      if(fIsMC>0){
        fHistoTrueSphericityvsRecSphericity[iCut]  = new TH2F("True Sphericity vs rec. Sphericity", "True Sphericity vs rec. Sphericity", 50, 0, 1, 50, 0, 1);
        fHistoTrueSphericityvsRecSphericity[iCut]->GetXaxis()->SetTitle("S_{true}");
        fHistoTrueSphericityvsRecSphericity[iCut]->GetYaxis()->SetTitle("S_{rec.}");
        fESDList[iCut]->Add(fHistoTrueSphericityvsRecSphericity[iCut]);
        fHistoTrueMultiplicityvsRecMultiplicity[iCut]  = new TH2F("True Multiplicity vs rec. Multiplicity", "True Multiplicity vs rec. Multiplicity", 100, 0, 100, 100, 0, 100);
        fHistoTrueMultiplicityvsRecMultiplicity[iCut]->GetXaxis()->SetTitle("N_{true}");
        fHistoTrueMultiplicityvsRecMultiplicity[iCut]->GetYaxis()->SetTitle("N_{rec.}");
        fESDList[iCut]->Add(fHistoTrueMultiplicityvsRecMultiplicity[iCut]);
        fHistoPionSpectrum[iCut]     = new TH1F("Charged pion spectrum", "Charged pion spectrum", 200, 0, 20);
        fHistoPionSpectrum[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoPionSpectrum[iCut]->GetYaxis()->SetTitle("Number of charged pions");
        fESDList[iCut]->Add(fHistoPionSpectrum[iCut]);
        fHistoProtonSpectrum[iCut]     = new TH1F("Proton spectrum", "Proton spectrum", 200, 0, 20);
        fHistoProtonSpectrum[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoProtonSpectrum[iCut]->GetYaxis()->SetTitle("Number of Protons");
        fESDList[iCut]->Add(fHistoProtonSpectrum[iCut]);
        fHistoKaonSpectrum[iCut]     = new TH1F("Charged kaon spectrum", "Charged kaon spectrum", 200, 0, 20);
        fHistoKaonSpectrum[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoKaonSpectrum[iCut]->GetYaxis()->SetTitle("Number of charged kaons");
        fESDList[iCut]->Add(fHistoKaonSpectrum[iCut]);
        fHistoNPionSpectrum[iCut]     = new TH1F("Neutral pion spectrum", "Neutral pion spectrum", 200, 0, 20);
        fHistoNPionSpectrum[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoNPionSpectrum[iCut]->GetYaxis()->SetTitle("Number of neutral pions");
        fESDList[iCut]->Add(fHistoNPionSpectrum[iCut]);
        fHistoEtaSpectrum[iCut]     = new TH1F("Eta spectrum", "Eta spectrum", 200, 0, 20);
        fHistoEtaSpectrum[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoEtaSpectrum[iCut]->GetYaxis()->SetTitle("Number of Etas");
        fESDList[iCut]->Add(fHistoEtaSpectrum[iCut]);
        fHistoDMesonSpectrum[iCut]     = new TH1F("D meson spectrum", "Charged kaon spectrum", 200, 0, 20);
        fHistoDMesonSpectrum[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoDMesonSpectrum[iCut]->GetYaxis()->SetTitle("Number of D mesons");
        fESDList[iCut]->Add(fHistoDMesonSpectrum[iCut]);
        if(fDoTrueSphericity){
            tTreeSphericity[iCut] = new TTree("Sphericity correlations", "Sphericity correlations");
            tTreeSphericity[iCut]->Branch("RecSph",&fRecSph,"fRecSph/F");
            tTreeSphericity[iCut]->Branch("TrueSph",&fTrueSph,"fTrueSph/F");
            tTreeSphericity[iCut]->Branch("Pi0Pt",&fPi0Pt,"fPi0Pt/F");
            tTreeSphericity[iCut]->Branch("Pi0InvMass",&fPi0InvMass,"fPi0InvMass/F");
            fESDList[iCut]->Add(tTreeSphericity[iCut]);
        }
      }
      fHistoEventSphericityvsHighpt[iCut]  = new TH2F("EventSphericity vs Highpt", "EventSphericity vs Highpt", 100, 0, 1, 100, 0, 50);
      fHistoEventSphericityvsHighpt[iCut]->GetXaxis()->SetTitle("S");
      fHistoEventSphericityvsHighpt[iCut]->GetYaxis()->SetTitle("Highest p_{T}");
      fESDList[iCut]->Add(fHistoEventSphericityvsHighpt[iCut]);
      fHistoEventSphericityvsTotalpt[iCut]  = new TH2F("EventSphericity vs Totalpt", "EventSphericity vs Totalpt", 100, 0, 1, 1000, 0, 100);
      fHistoEventSphericityvsTotalpt[iCut]->GetXaxis()->SetTitle("S");
      fHistoEventSphericityvsTotalpt[iCut]->GetYaxis()->SetTitle("Total p_{T}");
      fESDList[iCut]->Add(fHistoEventSphericityvsTotalpt[iCut]);
      fHistoEventSphericityvsMeanpt[iCut]  = new TH2F("EventSphericity vs Meanpt", "EventSphericity vs Meanpt", 100, 0, 1, 100, 0, 10);
      fHistoEventSphericityvsMeanpt[iCut]->GetXaxis()->SetTitle("S");
      fHistoEventSphericityvsMeanpt[iCut]->GetYaxis()->SetTitle("Mean p_{T}");
      fESDList[iCut]->Add(fHistoEventSphericityvsMeanpt[iCut]);
    }

    fHistoVertexZ[iCut]             = new TH1F("VertexZ", "VertexZ", 200, -10, 10);
    fHistoVertexZ[iCut]->GetXaxis()->SetTitle("Z_{vtx} (cm)");
    fESDList[iCut]->Add(fHistoVertexZ[iCut]);

    if(fIsHeavyIon == 1)
      fHistoNGammaCandidatesBasic[iCut]  = new TH1F("GammaCandidatesBasic", "GammaCandidatesBasic", 600, 0, 600);
    else if(fIsHeavyIon == 2)
      fHistoNGammaCandidatesBasic[iCut]  = new TH1F("GammaCandidatesBasic", "GammaCandidatesBasic", 400, 0, 400);
    else
      fHistoNGammaCandidatesBasic[iCut]  = new TH1F("GammaCandidatesBasic", "GammaCandidatesBasic", 100, 0, 100);
    fHistoNGammaCandidatesBasic[iCut]->GetXaxis()->SetTitle("#cluster candidates basic");
    fESDList[iCut]->Add(fHistoNGammaCandidatesBasic[iCut]);


    if(fIsHeavyIon == 1)
      fHistoNGammaCandidates[iCut]  = new TH1F("GammaCandidates", "GammaCandidates", 200, 0, 200);
    else if(fIsHeavyIon == 2)
      fHistoNGammaCandidates[iCut]  = new TH1F("GammaCandidates", "GammaCandidates", 100, 0, 100);
    else
      fHistoNGammaCandidates[iCut]  = new TH1F("GammaCandidates", "GammaCandidates", 50, 0, 50);
    fHistoNGammaCandidates[iCut]->GetXaxis()->SetTitle("#cluster candidates with current cut");
    fESDList[iCut]->Add(fHistoNGammaCandidates[iCut]);

    if(!fDoLightOutput){
      if(fIsHeavyIon == 1)
        fHistoNGoodESDTracksVsNGammaCandidates[iCut]  = new TH2F("GoodESDTracksVsGammaCandidates", "GoodESDTracksVsGammaCandidates", 4000, 0, 4000, 200, 0, 200);
      else if(fIsHeavyIon == 2)
        fHistoNGoodESDTracksVsNGammaCandidates[iCut]  = new TH2F("GoodESDTracksVsGammaCandidates", "GoodESDTracksVsGammaCandidates", 400, 0, 400, 100, 0, 100);
      else
        fHistoNGoodESDTracksVsNGammaCandidates[iCut]  = new TH2F("GoodESDTracksVsGammaCandidates", "GoodESDTracksVsGammaCandidates", 200, 0, 200, 50, 0, 50);
      fHistoNGoodESDTracksVsNGammaCandidates[iCut]->SetXTitle("#good tracks");
      fHistoNGoodESDTracksVsNGammaCandidates[iCut]->SetYTitle("#cluster candidates");
      fESDList[iCut]->Add(fHistoNGoodESDTracksVsNGammaCandidates[iCut]);

      fHistoSPDClusterTrackletBackground[iCut]        = new TH2F("SPD tracklets vs SPD clusters", "SPD tracklets vs SPD clusters", 100, 0, 200, 250, 0, 1000);
      fESDList[iCut]->Add(fHistoSPDClusterTrackletBackground[iCut]);

      if(fIsHeavyIon == 1)
        fHistoNV0Tracks[iCut]       = new TH1F("V0 Multiplicity", "V0 Multiplicity", 30000, 0, 30000);
      else if(fIsHeavyIon == 2)
        fHistoNV0Tracks[iCut]       = new TH1F("V0 Multiplicity", "V0 Multiplicity", 2500, 0, 2500);
      else
        fHistoNV0Tracks[iCut]       = new TH1F("V0 Multiplicity", "V0 Multiplicity", 1500, 0, 1500);
      fHistoNV0Tracks[iCut]->SetXTitle("V0 amplitude");
      fESDList[iCut]->Add(fHistoNV0Tracks[iCut]);
    }

    if(fIsHeavyIon==2) {
      fProfileEtaShift[iCut]        = new TProfile("Eta Shift", "Eta Shift", 1, -0.5, 0.5);
      fProfileEtaShift[iCut]->SetXTitle("#eta shift");
      fESDList[iCut]->Add(fProfileEtaShift[iCut]);
    }

    if (fIsMC > 1){
      fHistoNEvents[iCut]->Sumw2();
      fHistoNGoodESDTracks[iCut]->Sumw2();
      if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetUseSphericity()!=0){
        fHistoEventSphericity[iCut]->Sumw2();
        fHistoEventSphericityAxis[iCut]->Sumw2();
        fHistoEventSphericityvsNtracks[iCut]->Sumw2();
        if(fDoJetAnalysis){
          fHistoEventSphericityvsNJets[iCut]->Sumw2();
          fHistoEventMultiplicityvsNJets[iCut]->Sumw2();
        }
        fHistoTrueSphericityvsRecSphericity[iCut]->Sumw2();
        fHistoTrueMultiplicityvsRecMultiplicity[iCut]->Sumw2();
        fHistoEventSphericityvsHighpt[iCut]->Sumw2();
        fHistoEventSphericityvsTotalpt[iCut]->Sumw2();
        fHistoEventSphericityvsMeanpt[iCut]->Sumw2();
        fHistoPionSpectrum[iCut]->Sumw2();
        fHistoProtonSpectrum[iCut]->Sumw2();
        fHistoKaonSpectrum[iCut]->Sumw2();
        fHistoNPionSpectrum[iCut]->Sumw2();
        fHistoEtaSpectrum[iCut]->Sumw2();
        fHistoDMesonSpectrum[iCut]->Sumw2();
      }
      fHistoVertexZ[iCut]->Sumw2();
      fHistoNGammaCandidates[iCut]->Sumw2();
      fHistoNGammaCandidatesBasic[iCut]->Sumw2();
      if(!fDoLightOutput){
        fHistoNGoodESDTracksVsNGammaCandidates[iCut]->Sumw2();
        fHistoSPDClusterTrackletBackground[iCut]->Sumw2();
        fHistoNV0Tracks[iCut]->Sumw2();
      }
      if(fIsHeavyIon==2) fProfileEtaShift[iCut]->Sumw2();
    }

    fHistoClusGammaPt[iCut]               = new TH1F("ClusGamma_Pt", "ClusGamma_Pt", nBinsClusterPt, arrClusPtBinning);
    fHistoClusGammaPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
    fESDList[iCut]->Add(fHistoClusGammaPt[iCut]);
    fHistoClusGammaE[iCut]               = new TH1F("ClusGamma_E", "ClusGamma_E", nBinsClusterPt, arrClusPtBinning);
    fHistoClusGammaE[iCut]->SetXTitle("E_{clus} (GeV/c)");
    fESDList[iCut]->Add(fHistoClusGammaE[iCut]);
    fHistoClusOverlapHeadersGammaPt[iCut]   = new TH1F("ClusGammaOverlapHeaders_Pt", "ClusGammaOverlapHeaders_Pt", nBinsClusterPt, arrClusPtBinning);
    fHistoClusOverlapHeadersGammaPt[iCut]->SetXTitle("p_{T,clus} (GeV/c), selected header w/ overlap");
    fESDList[iCut]->Add(fHistoClusOverlapHeadersGammaPt[iCut]);
    fHistoClusAllHeadersGammaPt[iCut]       = new TH1F("ClusGammaAllHeaders_Pt", "ClusGammaAllHeaders_Pt", nBinsClusterPt, arrClusPtBinning);
    fHistoClusAllHeadersGammaPt[iCut]->SetXTitle("p_{T,clus} (GeV/c), all headers");
    fESDList[iCut]->Add(fHistoClusAllHeadersGammaPt[iCut]);
    fHistoClusRejectedHeadersGammaPt[iCut]  = new TH1F("ClusGammaRejectedHeaders_Pt", "ClusGammaRejectedHeaders_Pt", nBinsClusterPt, arrClusPtBinning);
    fHistoClusRejectedHeadersGammaPt[iCut]->SetXTitle("p_{T,clus} (GeV/c), rejected headers");
    fESDList[iCut]->Add(fHistoClusRejectedHeadersGammaPt[iCut]);
    if(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetClusterType()==2){
      if ( ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger()==6 ){
        fHistoGoodMesonClusters[iCut]                     = new TH1I( "fHistoGoodMesonClusters", "fHistoGoodMesonClusters", 7, 0.5, 7.5);
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
      //if ( (((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger()==1)||(((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger()==0) ){ //Only MB
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

    if(!fDoLightOutput && fDoClusterQA > 0){
      if(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetClusterType()==2){
        for (Int_t DDLRange=0; DDLRange<=fDDLRange_HistoClusGamma[1]-fDDLRange_HistoClusGamma[0]; DDLRange++){
          if (fIsMC==0){
            fHistoClusGammaE_DDL[DDLRange][iCut]    = new TH1F( Form("ClusGamma_E_DDL%d", (Int_t)(DDLRange+fDDLRange_HistoClusGamma[0])), Form("ClusGamma_E_DDL%d", (Int_t)(DDLRange+fDDLRange_HistoClusGamma[0])), nBinsClusterPt, arrClusPtBinning);
            fHistoClusGammaE_DDL[DDLRange][iCut]->SetXTitle("E_{clus} (GeV/c)");
            fESDList[iCut]->Add(fHistoClusGammaE_DDL[DDLRange][iCut]);
            //--------------------------------------------------
            //Only MB
            if ( (((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger()==1)||(((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger()==0) ){
              fHistoClusGammaE_DDL_TrBM[DDLRange][iCut]       = new TH1F( Form("fHistoClusGammaE_DDL%d_TrBM", (Int_t)(DDLRange+fDDLRange_HistoClusGamma[0])), Form("fHistoClusGammaE_DDL_TrBM%d", (Int_t)(DDLRange+fDDLRange_HistoClusGamma[0])), nBinsClusterPt, arrClusPtBinning);
              fHistoClusGammaE_DDL_TrBM[DDLRange][iCut]->SetXTitle("E_{clus} (GeV/c)");
              fESDList[iCut]->Add(fHistoClusGammaE_DDL_TrBM[DDLRange][iCut]);
              //----------
              fHistoClusGammaE_DDL_wL0_TrigEv[DDLRange][iCut] = new TH1F( Form("fHistoClusGammaE_DDL%d_wL0_TrigEv", (Int_t)(DDLRange+fDDLRange_HistoClusGamma[0])), Form("fHistoClusGammaE_DDL_wL0_TrigEv%d", (Int_t)(DDLRange+fDDLRange_HistoClusGamma[0])), nBinsClusterPt, arrClusPtBinning);
              fHistoClusGammaE_DDL_wL0_TrigEv[DDLRange][iCut]->SetXTitle("E_{clus} (GeV/c)");
              fESDList[iCut]->Add(fHistoClusGammaE_DDL_wL0_TrigEv[DDLRange][iCut]);
              //----------
              fHistoClusGammaE_DDL_woL0_TrigEv[DDLRange][iCut] = new TH1F( Form("fHistoClusGammaE_DDL%d_woL0_TrigEv", (Int_t)(DDLRange+fDDLRange_HistoClusGamma[0])), Form("fHistoClusGammaE_DDL_woL0_TrigEv%d", (Int_t)(DDLRange+fDDLRange_HistoClusGamma[0])), nBinsClusterPt, arrClusPtBinning);
              fHistoClusGammaE_DDL_woL0_TrigEv[DDLRange][iCut]->SetXTitle("E_{clus} (GeV/c)");
              fESDList[iCut]->Add(fHistoClusGammaE_DDL_woL0_TrigEv[DDLRange][iCut]);
              //----------
              fHistoClusGammaE_DDL_wL0_TrigEv_TrBM[DDLRange][iCut] = new TH1F( Form("fHistoClusGammaE_DDL%d_wL0_TrigEv_TrBM", (Int_t)(DDLRange+fDDLRange_HistoClusGamma[0])), Form("fHistoClusGammaE_DDL_wL0_TrigEv_TrBM%d", (Int_t)(DDLRange+fDDLRange_HistoClusGamma[0])), nBinsClusterPt, arrClusPtBinning);
              fHistoClusGammaE_DDL_wL0_TrigEv_TrBM[DDLRange][iCut]->SetXTitle("E_{clus} (GeV/c)");
              fESDList[iCut]->Add(fHistoClusGammaE_DDL_wL0_TrigEv_TrBM[DDLRange][iCut]);
              //----------
              fHistoClusGammaE_DDL_woL0_TrigEv_TrBM[DDLRange][iCut] = new TH1F( Form("fHistoClusGammaE_DDL%d_woL0_TrigEv_TrBM", (Int_t)(DDLRange+fDDLRange_HistoClusGamma[0])), Form("fHistoClusGammaE_DDL_woL0_TrigEv_TrBM%d", (Int_t)(DDLRange+fDDLRange_HistoClusGamma[0])), nBinsClusterPt, arrClusPtBinning);
              fHistoClusGammaE_DDL_woL0_TrigEv_TrBM[DDLRange][iCut]->SetXTitle("E_{clus} (GeV/c)");
              fESDList[iCut]->Add(fHistoClusGammaE_DDL_woL0_TrigEv_TrBM[DDLRange][iCut]);
            }//Only MB
          }
        }
      }
      fHistoClusGammaPtM02[iCut]               = new TH2F("ClusGamma_Pt_M02", "ClusGamma_Pt_M02", nBinsClusterPt, arrClusPtBinning, 100, 0, 1);
      fHistoClusGammaPtM02[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
      fHistoClusGammaPtM02[iCut]->SetYTitle("#sigma_{long}^{2}");
      fESDList[iCut]->Add(fHistoClusGammaPtM02[iCut]);

    }

    if (fIsMC > 1){
      fHistoClusGammaPt[iCut]->Sumw2();
      fHistoClusGammaE[iCut]->Sumw2();
      if(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetClusterType()==2){
        fHistoClusGammaE_BothBM_highestE[iCut]->Sumw2();
        //if ( (((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger()==1)||(((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger()==0) ){ //Only MB
          fHistoClusGammaE_BothBM[iCut]->Sumw2();
          fHistoClusGammaE_AnaBM_highestE[iCut]->Sumw2();
        //}
        if ( ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger()==6 ){
          fHistoClusGammaE_onlyTriggered[iCut]->Sumw2();
        }
      }
      if(!fDoLightOutput && fDoClusterQA > 0){
        if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType()==2){
          //Currently Taken out for MC, maybe add back in later on -> only outcomment
          //for (Int_t DDLRange=0; DDLRange<=fDDLRange_HistoClusGamma[1]-fDDLRange_HistoClusGamma[0]; DDLRange++){
          //  fHistoClusGammaE_DDL[DDLRange][iCut]->Sumw2();
          //}
          if ( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsSpecialTrigger()==6 ){
            fHistoGoodMesonClusters[iCut]->Sumw2();
          }
        }
      }
      fHistoClusOverlapHeadersGammaPt[iCut]->Sumw2();
      fHistoClusAllHeadersGammaPt[iCut]->Sumw2();
      fHistoClusRejectedHeadersGammaPt[iCut]->Sumw2();
      if(!fDoLightOutput && fDoClusterQA > 0)fHistoClusGammaPtM02[iCut]->Sumw2();
    }

    if(fDoMesonAnalysis){
      fHistoMotherInvMassPt[iCut]           = new TH2F("ESD_Mother_InvMass_Pt", "ESD_Mother_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
      fHistoMotherInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
      fHistoMotherInvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
      fESDList[iCut]->Add(fHistoMotherInvMassPt[iCut]);
      fHistoMotherBackInvMassPt[iCut]       = new TH2F("ESD_Background_InvMass_Pt", "ESD_Background_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
      fHistoMotherBackInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
      fHistoMotherBackInvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
      fESDList[iCut]->Add(fHistoMotherBackInvMassPt[iCut]);
      if(!fDoLightOutput || fDoPi0Only || fDoECalibOutput){
        fHistoMotherInvMassECalib[iCut]         = new TH2F("ESD_Mother_InvMass_E_Calib", "ESD_Mother_InvMass_E_Calib", 300, 0, 0.3, nBinsPt, arrPtBinning);
        fHistoMotherInvMassECalib[iCut]->SetXTitle("M_{inv} (GeV/c^{2})");
        fHistoMotherInvMassECalib[iCut]->SetYTitle("E_{cluster}(GeV)");
        fESDList[iCut]->Add(fHistoMotherInvMassECalib[iCut]);
        fHistoMotherBackInvMassECalib[iCut]     = new TH2F("ESD_Background_InvMass_E_Calib", "ESD_Background_InvMass_E_Calib", 300, 0, 0.3, nBinsPt, arrPtBinning);
        fHistoMotherBackInvMassECalib[iCut]->SetXTitle("M_{inv} (GeV/c^{2})");
        fHistoMotherBackInvMassECalib[iCut]->SetYTitle("E_{cluster}(GeV)");
        fESDList[iCut]->Add(fHistoMotherBackInvMassECalib[iCut]);
      }

      if (fIsMC > 1){
        fHistoMotherInvMassPt[iCut]->Sumw2();
        fHistoMotherBackInvMassPt[iCut]->Sumw2();
        if(!fDoLightOutput || fDoPi0Only || fDoECalibOutput){
          fHistoMotherInvMassECalib[iCut]->Sumw2();
          fHistoMotherBackInvMassECalib[iCut]->Sumw2();
        }
      }

      if (fDoMesonQA > 0 && fDoMesonQA < 3 ){
        fHistoMotherPi0PtY[iCut]          = new TH2F("ESD_MotherPi0_Pt_Y", "ESD_MotherPi0_Pt_Y", nBinsQAPt, arrQAPtBinning, 150, -1.5, 1.5);
        fHistoMotherPi0PtY[iCut]->SetXTitle("p_{T} (GeV/c)");
        fHistoMotherPi0PtY[iCut]->SetYTitle("y_{#pi^{0}}");
        fESDList[iCut]->Add(fHistoMotherPi0PtY[iCut]);
        fHistoMotherPi0PtAlpha[iCut]      = new TH2F("ESD_MotherPi0_Pt_Alpha", "ESD_MotherPi0_Pt_Alpha", nBinsQAPt, arrQAPtBinning, 100, 0, 1);
        fHistoMotherPi0PtAlpha[iCut]->SetXTitle("p_{T} (GeV/c)");
        fHistoMotherPi0PtAlpha[iCut]->SetYTitle("#alpha_{#pi^{0}}");
        fESDList[iCut]->Add(fHistoMotherPi0PtAlpha[iCut]);
        fHistoMotherPi0PtOpenAngle[iCut]  = new TH2F("ESD_MotherPi0_Pt_OpenAngle", "ESD_MotherPi0_Pt_OpenAngle", nBinsQAPt, arrQAPtBinning,100,0, 0.5);
        fHistoMotherPi0PtOpenAngle[iCut]->SetXTitle("p_{T} (GeV/c)");
        fHistoMotherPi0PtOpenAngle[iCut]->SetYTitle("#theta_{#pi^{0}}");
        fESDList[iCut]->Add(fHistoMotherPi0PtOpenAngle[iCut]);

        if( !fDoPi0Only ){
          fHistoMotherEtaPtY[iCut]          = new TH2F("ESD_MotherEta_Pt_Y", "ESD_MotherEta_Pt_Y", nBinsQAPt, arrQAPtBinning, 150, -1.5, 1.5);
          fHistoMotherEtaPtY[iCut]->SetXTitle("p_{T} (GeV/c)");
          fHistoMotherEtaPtY[iCut]->SetYTitle("y_{#eta}");
          fESDList[iCut]->Add(fHistoMotherEtaPtY[iCut]);
          fHistoMotherEtaPtAlpha[iCut]      = new TH2F("ESD_MotherEta_Pt_Alpha", "ESD_MotherEta_Pt_Alpha", nBinsQAPt, arrQAPtBinning, 100, 0, 1);
          fHistoMotherEtaPtAlpha[iCut]->SetXTitle("p_{T} (GeV/c)");
          fHistoMotherEtaPtAlpha[iCut]->SetYTitle("#alpha_{#eta}");
          fESDList[iCut]->Add(fHistoMotherEtaPtAlpha[iCut]);
          fHistoMotherEtaPtOpenAngle[iCut]  = new TH2F("ESD_MotherEta_Pt_OpenAngle", "ESD_MotherEta_Pt_OpenAngle", nBinsQAPt, arrQAPtBinning,180,0, 1.8);
          fHistoMotherEtaPtOpenAngle[iCut]->SetXTitle("p_{T} (GeV/c)");
          fHistoMotherEtaPtOpenAngle[iCut]->SetYTitle("#theta_{#eta}");
          fESDList[iCut]->Add(fHistoMotherEtaPtOpenAngle[iCut]);
        }
        if(fIsHeavyIon == 1){
          fHistoMotherPi0NGoodESDTracksPt[iCut]  = new TH2F("ESD_MotherPi0_GoodESDTracks_Pt", "ESD_MotherPi0_GoodESDTracks_Pt", 4000, 0, 4000, nBinsQAPt, arrQAPtBinning);
          fHistoMotherPi0NGoodESDTracksPt[iCut]->SetXTitle("# good tracks");
          fHistoMotherPi0NGoodESDTracksPt[iCut]->SetYTitle("p_{#pi^{0}, T} (GeV/c)");
          if( !fDoPi0Only ){
            fHistoMotherEtaNGoodESDTracksPt[iCut]  = new TH2F("ESD_MotherEta_GoodESDTracks_Pt", "ESD_MotherEta_GoodESDTracks_Pt", 4000, 0, 4000, nBinsQAPt, arrQAPtBinning);
            fHistoMotherEtaNGoodESDTracksPt[iCut]->SetXTitle("# good tracks");
            fHistoMotherEtaNGoodESDTracksPt[iCut]->SetYTitle("p_{#eta, T} (GeV/c)");
          }
        }else if(fIsHeavyIon == 2){
          fHistoMotherPi0NGoodESDTracksPt[iCut]  = new TH2F("ESD_MotherPi0_GoodESDTracks_Pt", "ESD_MotherPi0_GoodESDTracks_Pt", 400, 0, 400, nBinsQAPt, arrQAPtBinning);
          fHistoMotherPi0NGoodESDTracksPt[iCut]->SetXTitle("# good tracks");
          fHistoMotherPi0NGoodESDTracksPt[iCut]->SetYTitle("p_{#pi^{0}, T} (GeV/c)");
          if( !fDoPi0Only ){
            fHistoMotherEtaNGoodESDTracksPt[iCut]  = new TH2F("ESD_MotherEta_GoodESDTracks_Pt", "ESD_MotherEta_GoodESDTracks_Pt", 400, 0, 400, nBinsQAPt, arrQAPtBinning);
            fHistoMotherEtaNGoodESDTracksPt[iCut]->SetXTitle("# good tracks");
            fHistoMotherEtaNGoodESDTracksPt[iCut]->SetYTitle("p_{#eta, T} (GeV/c)");
          }
        }else{
          fHistoMotherPi0NGoodESDTracksPt[iCut]  = new TH2F("ESD_MotherPi0_GoodESDTracks_Pt", "ESD_MotherPi0_GoodESDTracks_Pt", 200, 0, 200, nBinsQAPt, arrQAPtBinning);
          fHistoMotherPi0NGoodESDTracksPt[iCut]->SetXTitle("# good tracks");
          fHistoMotherPi0NGoodESDTracksPt[iCut]->SetYTitle("p_{#pi^{0}, T} (GeV/c)");
          if( !fDoPi0Only ){
            fHistoMotherEtaNGoodESDTracksPt[iCut]  = new TH2F("ESD_MotherEta_GoodESDTracks_Pt", "ESD_MotherEta_GoodESDTracks_Pt", 200, 0, 200, nBinsQAPt, arrQAPtBinning);
            fHistoMotherEtaNGoodESDTracksPt[iCut]->SetXTitle("# good tracks");
            fHistoMotherEtaNGoodESDTracksPt[iCut]->SetYTitle("p_{#eta, T} (GeV/c)");
          }
        }
        fESDList[iCut]->Add(fHistoMotherPi0NGoodESDTracksPt[iCut]);
        if( !fDoPi0Only ) fESDList[iCut]->Add(fHistoMotherEtaNGoodESDTracksPt[iCut]);
      }
      if (fDoMesonQA == 2){
        fHistoMotherPtOpenAngle[iCut]  = new TH2F("ESD_Mother_Pt_OpenAngle", "ESD_Mother_Pt_OpenAngle", nBinsQAPt, arrQAPtBinning,180,0, 1.8);
        fHistoMotherPtOpenAngle[iCut]->SetXTitle("p_{T} (GeV/c)");
        fHistoMotherPtOpenAngle[iCut]->SetYTitle("#theta");
        fESDList[iCut]->Add(fHistoMotherPtOpenAngle[iCut]);
        fHistoMotherPtOpenAngleBck[iCut]  = new TH2F("ESD_MotherBck_Pt_OpenAngle", "ESD_MotherBck_Pt_OpenAngle", nBinsQAPt, arrQAPtBinning,180,0, 1.8);
        fHistoMotherPtOpenAngleBck[iCut]->SetXTitle("p_{T} (GeV/c)");
        fHistoMotherPtOpenAngleBck[iCut]->SetYTitle("#theta");
        fESDList[iCut]->Add(fHistoMotherPtOpenAngleBck[iCut]);
        if (fIsMC == 2){
          fHistoMotherPtOpenAngle[iCut]->Sumw2();
          fHistoMotherPtOpenAngleBck[iCut]->Sumw2();
        }
      }

      if (fIsMC > 1 && fDoMesonQA > 0 && fDoMesonQA < 3){
        fHistoMotherPi0PtY[iCut]->Sumw2();
        fHistoMotherPi0PtAlpha[iCut]->Sumw2();
        fHistoMotherPi0PtOpenAngle[iCut]->Sumw2();
        fHistoMotherPi0NGoodESDTracksPt[iCut]->Sumw2();
        if( !fDoPi0Only ){
          fHistoMotherEtaPtY[iCut]->Sumw2();
          fHistoMotherEtaPtAlpha[iCut]->Sumw2();
          fHistoMotherEtaPtOpenAngle[iCut]->Sumw2();
          fHistoMotherEtaNGoodESDTracksPt[iCut]->Sumw2();
        }
      }

      if (fProduceCellIDPlots){
        Int_t nMaxCells      = 12*48*24;
        if(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetClusterType() == 2) nMaxCells = 5*56*64;
        fHistCellIDvsClusterEnergy[iCut]         = new TH2F("CellIDvsClusterEnergy", "CellIDvsClusterEnergy", 100, 0.5, 100., nMaxCells, 0, nMaxCells);
        fHistCellIDvsClusterEnergy[iCut]->SetXTitle("E_{clus} (GeV)");
        fHistCellIDvsClusterEnergy[iCut]->SetYTitle("Cell ID");
        SetLogBinningXTH2(fHistCellIDvsClusterEnergy[iCut]);
        fESDList[iCut]->Add(fHistCellIDvsClusterEnergy[iCut]);
        fHistCellIDvsClusterEnergyMax[iCut]      = new TH2F("CellIDvsClusterEnergyMax", "CellIDvsClusterEnergyMax", 100, 0.5, 100., nMaxCells, 0, nMaxCells);
        fHistCellIDvsClusterEnergyMax[iCut]->SetXTitle("E_{clus} (GeV)");
        fHistCellIDvsClusterEnergyMax[iCut]->SetYTitle("Cell ID");
        SetLogBinningXTH2(fHistCellIDvsClusterEnergyMax[iCut]);
        fESDList[iCut]->Add(fHistCellIDvsClusterEnergyMax[iCut]);
      }

      if (fDoMesonQA == 4 && fIsMC == 0){
        fTreeList[iCut] = new TList();
        fTreeList[iCut]->SetName(Form("%s_%s_%s InvMass Tree", cutstringEvent.Data(), cutstringCalo.Data(), cutstringMeson.Data()));
        fTreeList[iCut]->SetOwner(kTRUE);
        fCutFolder[iCut]->Add(fTreeList[iCut]);

        tSigInvMassPtAlphaTheta[iCut] = new TTree("Sig_InvMass_Pt_Alpha_Theta_MixPool", "Sig_InvMass_Pt_Alpha_Theta_MixPool");
        tSigInvMassPtAlphaTheta[iCut]->Branch("InvMass",&fInvMassTreeInvMass,"fInvMassTreeInvMass/F");
        tSigInvMassPtAlphaTheta[iCut]->Branch("Pt",&fInvMassTreePt,"fInvMassTreePt/F");
        tSigInvMassPtAlphaTheta[iCut]->Branch("Alpha",&fInvMassTreeAlpha,"fInvMassTreeAlpha/F");
        tSigInvMassPtAlphaTheta[iCut]->Branch("Theta",&fInvMassTreeTheta,"fInvMassTreeTheta/F");
        tSigInvMassPtAlphaTheta[iCut]->Branch("MixPool",&fInvMassTreeMixPool,"fInvMassTreeMixPool/I");
        tSigInvMassPtAlphaTheta[iCut]->Branch("zVtx",&fInvMassTreeZVertex,"fInvMassTreeZVertex/F");
        tSigInvMassPtAlphaTheta[iCut]->Branch("Eta",&fInvMassTreeEta,"fInvMassTreeEta/F");
        fTreeList[iCut]->Add(tSigInvMassPtAlphaTheta[iCut]);

        tBckInvMassPtAlphaTheta[iCut] = new TTree("Bck_InvMass_Pt_Alpha_Theta_MixPool", "Bck_InvMass_Pt_Alpha_Theta_MixPool");
        tBckInvMassPtAlphaTheta[iCut]->Branch("InvMass",&fInvMassTreeInvMass,"fInvMassTreeInvMass/F");
        tBckInvMassPtAlphaTheta[iCut]->Branch("Pt",&fInvMassTreePt,"fInvMassTreePt/F");
        tBckInvMassPtAlphaTheta[iCut]->Branch("Alpha",&fInvMassTreeAlpha,"fInvMassTreeAlpha/F");
        tBckInvMassPtAlphaTheta[iCut]->Branch("Theta",&fInvMassTreeTheta,"fInvMassTreeTheta/F");
        tBckInvMassPtAlphaTheta[iCut]->Branch("MixPool",&fInvMassTreeMixPool,"fInvMassTreeMixPool/I");
        tBckInvMassPtAlphaTheta[iCut]->Branch("zVtx",&fInvMassTreeZVertex,"fInvMassTreeZVertex/F");
        tBckInvMassPtAlphaTheta[iCut]->Branch("Eta",&fInvMassTreeEta,"fInvMassTreeEta/F");
        fTreeList[iCut]->Add(tBckInvMassPtAlphaTheta[iCut]);
      }

      if (fDoMesonQA == 5 && fIsMC == 0){
        tClusterTimingEff[iCut]   = new TTree(Form("%s_%s_%s ClusterTimingEff", cutstringEvent.Data(), cutstringCalo.Data(), cutstringMeson.Data()), Form("%s_%s_%s ClusterTimingEff", cutstringEvent.Data(), cutstringCalo.Data(), cutstringMeson.Data()));
        tClusterTimingEff[iCut]->Branch("InvMass",&fInvMassTreeInvMass,"fInvMassTreeInvMass/F");
        tClusterTimingEff[iCut]->Branch("Pt",&fInvMassTreePt,"fInvMassTreePt/F");
        tClusterTimingEff[iCut]->Branch("ClusterTimeTag",&fClusterTimeTag,"fClusterTimeTag/F");
        tClusterTimingEff[iCut]->Branch("ClusterTimeProbe",&fClusterTimeProbe,"fClusterTimeProbe/F");
        tClusterTimingEff[iCut]->Branch("ClusterETag",&fClusterETag,"fClusterETag/F");
        tClusterTimingEff[iCut]->Branch("ClusterEProbe",&fClusterEProbe,"fClusterEProbe/F");
      }

      if (fProduceTreeEOverP ){
        fClusterTreeList[iCut] = new TList();
        fClusterTreeList[iCut]->SetName(Form("%s_%s EoverP Tree",cutstringEvent.Data(),cutstringCalo.Data()));
        fClusterTreeList[iCut]->SetOwner(kTRUE);
        fCutFolder[iCut]->Add(fClusterTreeList[iCut]);

        tClusterEOverP[iCut] = new TTree("EOverP_ClusE_ClusM02_ClusM20_TrackP_TrackPt", "EOverP_ClusE_ClusM02_ClusM20_TrackP_TrackPt");
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
    if(fDoJetAnalysis){

      fJetHistograms[iCut] = new TList();
      fJetHistograms[iCut]->SetOwner(kTRUE);
      fJetHistograms[iCut]->SetName(Form("%s_%s_%s Jet histograms", cutstringEvent.Data(), cutstringCalo.Data(), cutstringMeson.Data()));

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
      fHistoEventwJets[iCut] = new TH1F("NEvents_with_Jets", "NEvents_with_Jets", 4, 0, 4);
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

    if( ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->DoSectorMixing() ){
      fV0Reader->SetCalcSector(kTRUE);
    }

    if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType()==2) {
        fCaloTriggerMimicHelper[iCut] = NULL;
        if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetTriggerMimicking() > 0){
            fCaloTriggerMimicHelper[iCut] = (AliCaloTriggerMimicHelper*) (AliAnalysisManager::GetAnalysisManager()->GetTask(fCaloTriggerHelperName.Data()));
            if(fCaloTriggerMimicHelper[iCut]) {
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

    if(fDoMesonQA>=10){
      fHistoMCPi0GenVsNClus                   = new TH2F*[fnCuts];
      fHistoMCPi0GenFoundInOneCluster         = new TH2F*[fnCuts];
      fHistoMCPi0GenFoundInTwoCluster         = new TH2F*[fnCuts];
      fHistoMCEtaGenFoundInOneCluster         = new TH2F*[fnCuts];
      fHistoMCEtaGenFoundInTwoCluster         = new TH2F*[fnCuts];
      fHistoMCGammaConvRvsPt                  = new TH2F*[fnCuts];
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

      if(!fDoPi0Only){
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
      if(fIsMC > 0 && fDoClusterQA == kTRUE){
        fNumberOfClusters                              = new TH1F*[fnCuts];
        fNumberOfClustersinJets                        = new TH1F*[fnCuts];
        fEnergyRatio                                   = new TH2F*[fnCuts];
        fEnergyRatioinJets                             = new TH2F*[fnCuts];
        fEnergyRatioGamma1                             = new TH2F*[fnCuts];
        fEnergyRatioGamma1inJets                       = new TH2F*[fnCuts];
        fEnergyRatioGammaAnywhere                      = new TH2F*[fnCuts];
        fEnergyRatioGammaAnywhereinJets                = new TH2F*[fnCuts];
        fEnergyDeposit                                 = new TH2F*[fnCuts];
        fEnergyDepositinJets                           = new TH2F*[fnCuts];
        fEnergyDepGamma1                               = new TH2F*[fnCuts];
        fEnergyDepGamma1inJets                         = new TH2F*[fnCuts];
        fEnergyDepGammaAnywhere                        = new TH2F*[fnCuts];
        fEnergyDepGammaAnywhereinJets                  = new TH2F*[fnCuts];
        fEnergyRatioGamma1Helped                       = new TH1F*[fnCuts];
        fEnergyRatioGamma1HelpedinJets                 = new TH1F*[fnCuts];
        fClusterEtaPhiJets                             = new TH2F*[fnCuts];
      }
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
      if (fDoClusterQA > 1) {
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
      fHistoMCPi0InAccPt        = new TH1F*[fnCuts];

      if (fIsMC > 1){
        fHistoMCPi0WOEvtWeightPt       = new TH1F*[fnCuts];
        fHistoMCPi0WOEvtWeightInAccPt  = new TH1F*[fnCuts];
        if( !fDoPi0Only ){
          fHistoMCEtaWOEvtWeightPt       = new TH1F*[fnCuts];
          fHistoMCEtaWOEvtWeightInAccPt  = new TH1F*[fnCuts];
        }
      }

      fHistoTruePi0InvMassPt                    = new TH2F*[fnCuts];
      if (!fDoLightOutput){
          fHistoTruePi0InvMassPtAdditional            = new TH2F*[fnCuts];
      }
      fHistoDoubleCountTruePi0InvMassPt         = new TH2F*[fnCuts];
      fHistoTruePrimaryPi0InvMassPt             = new TH2F*[fnCuts];
      fHistoTruePrimaryPi0W0WeightingInvMassPt  = new TH2F*[fnCuts];
      fProfileTruePrimaryPi0WeightsInvMassPt    = new TProfile2D*[fnCuts];
      fHistoTrueSecondaryPi0InvMassPt           = new TH2F*[fnCuts];
      fHistoTrueSecondaryPi0FromK0sInvMassPt    = new TH2F*[fnCuts];
      fHistoTrueSecondaryPi0FromK0lInvMassPt    = new TH2F*[fnCuts];
      fHistoTrueSecondaryPi0FromLambdaInvMassPt = new TH2F*[fnCuts];
      fHistoTrueSecondaryPi0FromEtaInvMassPt    = new TH2F*[fnCuts];
      if( !fDoPi0Only ){
        fHistoMCEtaPt                             = new TH1F*[fnCuts];
        fHistoMCEtaWOWeightPt                     = new TH1F*[fnCuts];
        fHistoMCEtaInAccPt                        = new TH1F*[fnCuts];
        fHistoTrueEtaInvMassPt                    = new TH2F*[fnCuts];
        if (!fDoLightOutput){
            fHistoTrueEtaInvMassPtAdditional            = new TH2F*[fnCuts];
        }
        fHistoDoubleCountTrueEtaInvMassPt         = new TH2F*[fnCuts];
        fHistoTruePrimaryEtaInvMassPt             = new TH2F*[fnCuts];
        fHistoTruePrimaryEtaW0WeightingInvMassPt  = new TH2F*[fnCuts];
        fProfileTruePrimaryEtaWeightsInvMassPt    = new TProfile2D*[fnCuts];
        fHistoMCSecEtaPt                          = new TH1F*[fnCuts];
        fHistoMCSecEtaSource                      = new TH1F*[fnCuts];
      }
      if(!fDoLightOutput){
        fHistoTruePi0InvMassPtAlpha               = new TH2F*[fnCuts];
        fHistoTruePi0PureGammaInvMassPtAlpha      = new TH2F*[fnCuts];
      }
      fHistoMCPrimaryPtvsSource                   = new TH2F*[fnCuts];
      fHistoMCSecPi0PtvsSource                    = new TH2F*[fnCuts];
      fHistoMCSecPi0InAccPtvsSource               = new TH2F*[fnCuts];
      fHistoMCSecPi0Source                        = new TH1F*[fnCuts];

      if (fDoMesonQA > 0 && fDoMesonQA < 3 ){
        fHistoMCPi0PtY                            = new TH2F*[fnCuts];
        fHistoMCPi0PtAlpha                        = new TH2F*[fnCuts];
        if( !fDoPi0Only ){
          fHistoMCEtaPtY                          = new TH2F*[fnCuts];
          fHistoMCEtaPtAlpha                      = new TH2F*[fnCuts];
        }
        if (fIsMC == 2){
          fHistoMCPi0PtJetPt                      = new TH2F*[fnCuts];
          if( !fDoPi0Only )fHistoMCEtaPtJetPt     = new TH2F*[fnCuts];
        }

        if (fIsMC < 2){
          fHistoTruePi0CaloPhotonInvMassPt                  = new TH2F*[fnCuts];
          fHistoTruePi0CaloMixedPhotonConvPhotonInvMassPt   = new TH2F*[fnCuts];
          fHistoTruePi0CaloConvertedPhotonInvMassPt         = new TH2F*[fnCuts];
          fHistoTruePi0CaloElectronInvMassPt                = new TH2F*[fnCuts];
          fHistoTruePi0CaloMergedClusterInvMassPt           = new TH2F*[fnCuts];
          fHistoTruePi0CaloMergedClusterPartConvInvMassPt   = new TH2F*[fnCuts];
          fHistoTruePi0NonMergedElectronPhotonInvMassPt     = new TH2F*[fnCuts];
          fHistoTruePi0NonMergedElectronMergedPhotonInvMassPt   = new TH2F*[fnCuts];
          fHistoTruePrimaryPi0MCPtResolPt                   = new TH2F*[fnCuts];
          fHistoTrueK0sWithPi0DaughterMCPt                  = new TH1F*[fnCuts];
          fHistoTrueK0lWithPi0DaughterMCPt                  = new TH1F*[fnCuts];
          fHistoTrueLambdaWithPi0DaughterMCPt               = new TH1F*[fnCuts];
          if( !fDoPi0Only ){
            fHistoTrueEtaCaloPhotonInvMassPt                = new TH2F*[fnCuts];
            fHistoTrueEtaCaloMixedPhotonConvPhotonInvMassPt = new TH2F*[fnCuts];
            fHistoTrueEtaCaloConvertedPhotonInvMassPt       = new TH2F*[fnCuts];
            fHistoTrueEtaCaloElectronInvMassPt              = new TH2F*[fnCuts];
            fHistoTrueEtaCaloMergedClusterInvMassPt         = new TH2F*[fnCuts];
            fHistoTrueEtaCaloMergedClusterPartConvInvMassPt = new TH2F*[fnCuts];
            fHistoTruePrimaryEtaMCPtResolPt                 = new TH2F*[fnCuts];
            fHistoTrueEtaWithPi0DaughterMCPt                = new TH1F*[fnCuts];

          }
        }
        fHistoTruePi0PtY                  = new TH2F*[fnCuts];
        fHistoTruePi0PtAlpha              = new TH2F*[fnCuts];
        fHistoTruePi0PtOpenAngle          = new TH2F*[fnCuts];
        fHistoTrueBckGGInvMassPt          = new TH2F*[fnCuts];
        fHistoTrueBckFullMesonContainedInOneClusterInvMassPt = new TH2F*[fnCuts];
        fHistoTrueBckAsymEClustersInvMassPt                  = new TH2F*[fnCuts];
        fHistoTrueBckContInvMassPt        = new TH2F*[fnCuts];
        if( !fDoPi0Only ){
          fHistoTrueEtaPtOpenAngle          = new TH2F*[fnCuts];
          fHistoTrueEtaPtAlpha              = new TH2F*[fnCuts];
          fHistoTrueEtaPtY                  = new TH2F*[fnCuts];
        }
      }
      if (fDoMesonQA==2){
        fHistoTruePi0Category1            = new TH2F*[fnCuts];
        fHistoTruePi0Category2            = new TH2F*[fnCuts];
        fHistoTruePi0Category3            = new TH2F*[fnCuts];
        fHistoTruePi0Category4_6          = new TH2F*[fnCuts];
        fHistoTruePi0Category5            = new TH2F*[fnCuts];
        fHistoTruePi0Category7            = new TH2F*[fnCuts];
        fHistoTruePi0Category8            = new TH2F*[fnCuts];
        if( !fDoPi0Only ){
          fHistoTrueEtaCategory1            = new TH2F*[fnCuts];
          fHistoTrueEtaCategory2            = new TH2F*[fnCuts];
          fHistoTrueEtaCategory3            = new TH2F*[fnCuts];
          fHistoTrueEtaCategory4_6          = new TH2F*[fnCuts];
          fHistoTrueEtaCategory5            = new TH2F*[fnCuts];
          fHistoTrueEtaCategory7            = new TH2F*[fnCuts];
          fHistoTrueEtaCategory8            = new TH2F*[fnCuts];
        }
      }
    }


    for(Int_t iCut = 0; iCut<fnCuts;iCut++){
      TString cutstringEvent    = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
      TString cutstringCalo     = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
      TString cutstringMeson    = "NoMesonCut";
      if(fDoMesonAnalysis)
        cutstringMeson          = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();

      fMCList[iCut]                   = new TList();
      fMCList[iCut]->SetName(Form("%s_%s_%s MC histograms", cutstringEvent.Data(), cutstringCalo.Data(), cutstringMeson.Data()));
      fMCList[iCut]->SetOwner(kTRUE);
      fCutFolder[iCut]->Add(fMCList[iCut]);

      if(!fDoLightOutput){
        fHistoMCHeaders[iCut]             = new TH1I("MC_Headers", "MC_Headers", 20, 0, 20);
        fHistoMCHeaders[iCut]->SetXTitle("accepted headers");
        fMCList[iCut]->Add(fHistoMCHeaders[iCut]);
        fHistoMCAllGammaPt[iCut]          = new TH1F("MC_AllGamma_Pt", "MC_AllGamma_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoMCAllGammaPt[iCut]->SetXTitle("p_{T} (GeV/c)");
        fMCList[iCut]->Add(fHistoMCAllGammaPt[iCut]);
        fHistoMCAllSecondaryGammaPt[iCut] = new TH2F("MC_AllSecondaryGamma_Pt", "MC_AllSecondaryGamma_Pt", nBinsClusterPt, arrClusPtBinning, 5, -0.5, 4.5);
        fHistoMCAllSecondaryGammaPt[iCut]->GetYaxis()->SetBinLabel(1,"K0s");
        fHistoMCAllSecondaryGammaPt[iCut]->GetYaxis()->SetBinLabel(2,"K0l");
        fHistoMCAllSecondaryGammaPt[iCut]->GetYaxis()->SetBinLabel(3,"Lambda");
        fHistoMCAllSecondaryGammaPt[iCut]->GetYaxis()->SetBinLabel(4,"Eta");
        fHistoMCAllSecondaryGammaPt[iCut]->GetYaxis()->SetBinLabel(5,"rest");
        fHistoMCAllSecondaryGammaPt[iCut]->SetXTitle("p_{T} (GeV/c)");
        fHistoMCAllSecondaryGammaPt[iCut]->SetYTitle("sec. particle");
        fMCList[iCut]->Add(fHistoMCAllSecondaryGammaPt[iCut]);
        fHistoMCDecayGammaPi0Pt[iCut]     = new TH1F("MC_DecayGammaPi0_Pt", "MC_DecayGammaPi0_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoMCDecayGammaPi0Pt[iCut]->SetXTitle("p_{T} (GeV/c)");
        fMCList[iCut]->Add(fHistoMCDecayGammaPi0Pt[iCut]);
        fHistoMCDecayGammaRhoPt[iCut]     = new TH1F("MC_DecayGammaRho_Pt", "MC_DecayGammaRho_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoMCDecayGammaRhoPt[iCut]->SetXTitle("p_{T} (GeV/c)");
        fMCList[iCut]->Add(fHistoMCDecayGammaRhoPt[iCut]);
        fHistoMCDecayGammaEtaPt[iCut]     = new TH1F("MC_DecayGammaEta_Pt", "MC_DecayGammaEta_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoMCDecayGammaEtaPt[iCut]->SetXTitle("p_{T} (GeV/c)");
        fMCList[iCut]->Add(fHistoMCDecayGammaEtaPt[iCut]);
        fHistoMCDecayGammaOmegaPt[iCut]   = new TH1F("MC_DecayGammaOmega_Pt", "MC_DecayGammaOmmega_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoMCDecayGammaOmegaPt[iCut]->SetXTitle("p_{T} (GeV/c)");
        fMCList[iCut]->Add(fHistoMCDecayGammaOmegaPt[iCut]);
        fHistoMCDecayGammaEtapPt[iCut]    = new TH1F("MC_DecayGammaEtap_Pt", "MC_DecayGammaEtap_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoMCDecayGammaEtapPt[iCut]->SetXTitle("p_{T} (GeV/c)");
        fMCList[iCut]->Add(fHistoMCDecayGammaEtapPt[iCut]);
        fHistoMCDecayGammaPhiPt[iCut]     = new TH1F("MC_DecayGammaPhi_Pt", "MC_DecayGammaPhi_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoMCDecayGammaPhiPt[iCut]->SetXTitle("p_{T} (GeV/c)");
        fMCList[iCut]->Add(fHistoMCDecayGammaPhiPt[iCut]);
        fHistoMCDecayGammaSigmaPt[iCut]   = new TH1F("MC_DecayGammaSigma_Pt", "MC_DecayGammaSigma_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoMCDecayGammaSigmaPt[iCut]->SetXTitle("p_{T} (GeV/c)");
        fMCList[iCut]->Add(fHistoMCDecayGammaSigmaPt[iCut]);

        if(fDoMesonQA>=10){
          fHistoMCPi0GenVsNClus[iCut]   = new TH2F("MC_Pi0GenVsNClus", "MC_Pi0GenVsNClus", 400, 0,100,50,0,50);
          fHistoMCPi0GenVsNClus[iCut]->SetXTitle("p_{T} (GeV/c)");
          fHistoMCPi0GenVsNClus[iCut]->SetYTitle("Nclus");
          fMCList[iCut]->Add(fHistoMCPi0GenVsNClus[iCut]);

          fHistoMCPi0GenFoundInOneCluster[iCut]   = new TH2F("MC_Pi0GenFoundInOneCluster_Pt", "MC_Pi0GenFoundInOneCluster_Pt", 400, 0,100,20,-0.9,0.9);
          fHistoMCPi0GenFoundInOneCluster[iCut]->SetXTitle("p_{T} (GeV/c)");
          fHistoMCPi0GenFoundInOneCluster[iCut]->SetYTitle("#eta");
          fMCList[iCut]->Add(fHistoMCPi0GenFoundInOneCluster[iCut]);

          fHistoMCPi0GenFoundInTwoCluster[iCut]   = new TH2F("MC_Pi0GenFoundInTwoCluster_Pt", "MC_Pi0GenFoundInTwoCluster_Pt", 400, 0,100,20,-0.9,0.9);
          fHistoMCPi0GenFoundInTwoCluster[iCut]->SetXTitle("p_{T} (GeV/c)");
          fHistoMCPi0GenFoundInTwoCluster[iCut]->SetYTitle("#eta");
          fMCList[iCut]->Add(fHistoMCPi0GenFoundInTwoCluster[iCut]);

          fHistoMCEtaGenFoundInOneCluster[iCut]   = new TH2F("MC_EtaGenFoundInOneCluster_Pt", "MC_EtaGenFoundInOneCluster_Pt", 400, 0,100,20,-0.9,0.9);
          fHistoMCEtaGenFoundInOneCluster[iCut]->SetXTitle("p_{T} (GeV/c)");
          fHistoMCEtaGenFoundInOneCluster[iCut]->SetYTitle("#eta");
          fMCList[iCut]->Add(fHistoMCEtaGenFoundInOneCluster[iCut]);

          fHistoMCEtaGenFoundInTwoCluster[iCut]   = new TH2F("MC_EtaGenFoundInTwoCluster_Pt", "MC_EtaGenFoundInTwoCluster_Pt", 400, 0,100,20,-0.9,0.9);
          fHistoMCEtaGenFoundInTwoCluster[iCut]->SetXTitle("p_{T} (GeV/c)");
          fHistoMCEtaGenFoundInTwoCluster[iCut]->SetYTitle("#eta");
          fMCList[iCut]->Add(fHistoMCEtaGenFoundInTwoCluster[iCut]);

          fHistoMCGammaConvRvsPt[iCut]   = new TH2F("MC_GammaConvRvsP", "MC_GammaConvRvsP", 400, 0,100,350,0,700);
          fHistoMCGammaConvRvsPt[iCut]->SetXTitle("p_{T} (GeV/c)");
          fHistoMCGammaConvRvsPt[iCut]->SetYTitle("R (cm)");
          fMCList[iCut]->Add(fHistoMCGammaConvRvsPt[iCut]);
        }

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
        fHistoMCPi0Pt[iCut]           = new TH1F("MC_Pi0_Pt", "MC_Pi0_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fHistoMCPi0Pt[iCut]->SetXTitle("p_{T} (GeV/c)");
        fHistoMCPi0Pt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0Pt[iCut]);
        fHistoMCPi0WOWeightPt[iCut]   = new TH1F("MC_Pi0_WOWeights_Pt", "MC_Pi0_WOWeights_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fHistoMCPi0WOWeightPt[iCut]->Sumw2();
        fHistoMCPi0WOWeightPt[iCut]->SetXTitle("p_{T} (GeV/c)");
        fMCList[iCut]->Add(fHistoMCPi0WOWeightPt[iCut]);
        fHistoMCPi0InAccPt[iCut]      = new TH1F("MC_Pi0InAcc_Pt", "MC_Pi0InAcc_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fHistoMCPi0InAccPt[iCut]->SetXTitle("p_{T} (GeV/c)");
        fHistoMCPi0InAccPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0InAccPt[iCut]);
        if( !fDoPi0Only ){
          fHistoMCEtaPt[iCut]           = new TH1F("MC_Eta_Pt", "MC_Eta_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
          fHistoMCEtaPt[iCut]->SetXTitle("p_{T} (GeV/c)");
          fHistoMCEtaPt[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCEtaPt[iCut]);
          fHistoMCEtaWOWeightPt[iCut]   = new TH1F("MC_Eta_WOWeights_Pt", "MC_Eta_WOWeights_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
          fHistoMCEtaWOWeightPt[iCut]->SetXTitle("p_{T} (GeV/c)");
          fHistoMCEtaWOWeightPt[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCEtaWOWeightPt[iCut]);
          fHistoMCEtaInAccPt[iCut]      = new TH1F("MC_EtaInAcc_Pt", "MC_EtaInAcc_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
          fHistoMCEtaInAccPt[iCut]->SetXTitle("p_{T} (GeV/c)");
          fHistoMCEtaInAccPt[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCEtaInAccPt[iCut]);
        }
        if (fIsMC > 1){
          fHistoMCPi0WOWeightPt[iCut]->Sumw2();
          fHistoMCPi0WOEvtWeightPt[iCut] = new TH1F("MC_Pi0_WOEventWeights_Pt", "MC_Pi0_WOEventWeights_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
          fHistoMCPi0WOEvtWeightPt[iCut]->SetXTitle("p_{T} (GeV/c)");
          fMCList[iCut]->Add(fHistoMCPi0WOEvtWeightPt[iCut]);
          fHistoMCPi0WOEvtWeightInAccPt[iCut] = new TH1F("MC_Pi0WOEvtWeightInAcc_Pt", "MC_Pi0WOEvtWeightInAcc_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
          fHistoMCPi0WOEvtWeightInAccPt[iCut]->SetXTitle("p_{T} (GeV/c)");
          fMCList[iCut]->Add(fHistoMCPi0WOEvtWeightInAccPt[iCut]);
          if( !fDoPi0Only ){
            fHistoMCEtaWOWeightPt[iCut]->Sumw2();
            fHistoMCEtaWOEvtWeightPt[iCut] = new TH1F("MC_Eta_WOEventWeights_Pt", "MC_Eta_WOEventWeights_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
            fHistoMCEtaWOEvtWeightPt[iCut]->SetXTitle("p_{T} (GeV/c)");
            fMCList[iCut]->Add(fHistoMCEtaWOEvtWeightPt[iCut]);
            fHistoMCEtaWOEvtWeightInAccPt[iCut] = new TH1F("MC_EtaWOEvtWeightInAcc_Pt", "MC_EtaWOEvtWeightInAcc_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
            fHistoMCEtaWOEvtWeightInAccPt[iCut]->SetXTitle("p_{T} (GeV/c)");
            fMCList[iCut]->Add(fHistoMCEtaWOEvtWeightInAccPt[iCut]);
          }

          if (fDoMesonQA > 0  && fDoMesonQA < 3 && fIsMC == 2){
            fHistoMCPi0PtJetPt[iCut]  = new TH2F("MC_Pi0_Pt_JetPt", "MC_Pi0_Pt_JetPt", nBinsQAPt, arrQAPtBinning, 200, 0, 200);
            fHistoMCPi0PtJetPt[iCut]->SetXTitle("p_{T} (GeV/c)");
            fHistoMCPi0PtJetPt[iCut]->SetYTitle("p_{jet, T} (GeV/c)");
            fHistoMCPi0PtJetPt[iCut]->Sumw2();
            fMCList[iCut]->Add(fHistoMCPi0PtJetPt[iCut]);
            if( !fDoPi0Only ){
              fHistoMCEtaPtJetPt[iCut]  = new TH2F("MC_Eta_Pt_JetPt", "MC_Eta_Pt_JetPt", nBinsQAPt, arrQAPtBinning, 200, 0, 200);
              fHistoMCEtaPtJetPt[iCut]->SetXTitle("p_{T} (GeV/c)");
              fHistoMCEtaPtJetPt[iCut]->SetYTitle("p_{jet, T} (GeV/c)");
              fHistoMCEtaPtJetPt[iCut]->Sumw2();
              fMCList[iCut]->Add(fHistoMCEtaPtJetPt[iCut]);
            }
          }
        }
        fHistoMCPrimaryPtvsSource[iCut]   = new TH2F("MC_Primary_Pt_Source", "MC_Primary_Pt_Source", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt, 7, -0.5, 6.5);
        fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(1,"Pi+");
        fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(2,"Pi-");
        fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(3,"K+");
        fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(4,"K-");
        fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(5,"K0s");
        fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(6,"K0l");
        fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(7,"Lambda");
        fHistoMCPrimaryPtvsSource[iCut]->SetXTitle("p_{T} (GeV/c)");
        fHistoMCPrimaryPtvsSource[iCut]->SetYTitle("particle");
        fMCList[iCut]->Add(fHistoMCPrimaryPtvsSource[iCut]);

        fHistoMCSecPi0Source[iCut]      = new TH1F("MC_SecPi0_Source", "MC_SecPi0_Source", 5000, 0., 5000);
        fHistoMCSecPi0Source[iCut]->SetYTitle("source PDG");
        fMCList[iCut]->Add(fHistoMCSecPi0Source[iCut]);
        if( !fDoPi0Only ){
          fHistoMCSecEtaSource[iCut]      = new TH1F("MC_SecEta_Source", "MC_SecEta_Source", 5000, 0, 5000);
          fHistoMCSecEtaSource[iCut]->SetYTitle("source PDG");
          fMCList[iCut]->Add(fHistoMCSecEtaSource[iCut]);
        }
        fHistoMCSecPi0PtvsSource[iCut]  = new TH2F("MC_SecPi0_Pt_Source", "MC_SecPi0_Pt_Source", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt, 16, -0.5, 15.5);
        fHistoMCSecPi0PtvsSource[iCut]->SetXTitle("p_{T} (GeV/c)");
        fHistoMCSecPi0PtvsSource[iCut]->SetYTitle("source");
        fMCList[iCut]->Add(fHistoMCSecPi0PtvsSource[iCut]);
        fHistoMCSecPi0InAccPtvsSource[iCut]  = new TH2F("MC_SecPi0InAcc_Pt_Source", "MC_SecPi0InAcc_Pt_Source", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt,  16, -0.5, 15.5);
        fHistoMCSecPi0InAccPtvsSource[iCut]->SetXTitle("p_{T} (GeV/c)");
        fHistoMCSecPi0InAccPtvsSource[iCut]->SetYTitle("source");
        fMCList[iCut]->Add(fHistoMCSecPi0InAccPtvsSource[iCut]);
        if( !fDoPi0Only ){
          fHistoMCSecEtaPt[iCut]          = new TH1F("MC_SecEta_Pt", "MC_SecEta_Pt", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
          fHistoMCSecEtaPt[iCut]->SetXTitle("p_{T} (GeV/c)");
          fMCList[iCut]->Add(fHistoMCSecEtaPt[iCut]);
        }
        if (fIsMC == 2) {
          fHistoMCPrimaryPtvsSource[iCut]->Sumw2();
          fHistoMCSecPi0PtvsSource[iCut]->Sumw2();
          fHistoMCSecPi0InAccPtvsSource[iCut]->Sumw2();
          if( !fDoPi0Only ) fHistoMCSecEtaPt[iCut]->Sumw2();
        }


        if (fDoMesonQA > 0 && fDoMesonQA < 3){
          fHistoMCPi0PtY[iCut]        = new TH2F("MC_Pi0_Pt_Y", "MC_Pi0_Pt_Y", nBinsQAPt, arrQAPtBinning, 150, -1.5, 1.5);
          fHistoMCPi0PtY[iCut]->SetXTitle("p_{T} (GeV/c)");
          fHistoMCPi0PtY[iCut]->SetYTitle("y");
          fHistoMCPi0PtY[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCPi0PtY[iCut]);
          if( !fDoPi0Only ){
            fHistoMCEtaPtY[iCut]        = new TH2F("MC_Eta_Pt_Y", "MC_Eta_Pt_Y", nBinsQAPt, arrQAPtBinning, 150, -1.5, 1.5);
            fHistoMCEtaPtY[iCut]->SetXTitle("p_{T} (GeV/c)");
            fHistoMCEtaPtY[iCut]->SetYTitle("y");
            fHistoMCEtaPtY[iCut]->Sumw2();
            fMCList[iCut]->Add(fHistoMCEtaPtY[iCut]);
          }
          fHistoMCPi0PtAlpha[iCut]    = new TH2F("MC_Pi0_Pt_Alpha", "MC_Pi0_Pt_Alpha", nBinsQAPt, arrQAPtBinning, 100, 0, 1);
          fHistoMCPi0PtAlpha[iCut]->SetXTitle("p_{T} (GeV/c)");
          fHistoMCPi0PtAlpha[iCut]->SetYTitle("#alpha");
          fMCList[iCut]->Add(fHistoMCPi0PtAlpha[iCut]);
          if( !fDoPi0Only ){
            fHistoMCEtaPtAlpha[iCut]    = new TH2F("MC_Eta_Pt_Alpha", "MC_Eta_Pt_Alpha", nBinsQAPt, arrQAPtBinning, 100, 0, 1);
            fHistoMCEtaPtAlpha[iCut]->SetXTitle("p_{T} (GeV/c)");
            fHistoMCEtaPtAlpha[iCut]->SetYTitle("#alpha");
            fMCList[iCut]->Add(fHistoMCEtaPtAlpha[iCut]);
          }

          if (fIsMC == 2) {
            fHistoMCPi0PtAlpha[iCut]->Sumw2();
            if( !fDoPi0Only )fHistoMCEtaPtAlpha[iCut]->Sumw2();
          }
        }
      }
      fTrueList[iCut]                 = new TList();
      fTrueList[iCut]->SetName(Form("%s_%s_%s True histograms", cutstringEvent.Data(), cutstringCalo.Data(), cutstringMeson.Data()));
      fTrueList[iCut]->SetOwner(kTRUE);
      fCutFolder[iCut]->Add(fTrueList[iCut]);

      if(fDoJetAnalysis && !fDoLightOutput){
        fTrueJetHistograms[iCut] = new TList();
        fTrueJetHistograms[iCut]->SetName(Form("%s_%s_%s True Jet histograms", cutstringEvent.Data(), cutstringCalo.Data(), cutstringMeson.Data()));
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

        if(!fDoPi0Only){
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
        if(fIsMC > 0 && fDoClusterQA == kTRUE){
          fNumberOfClusters[iCut]      = new TH1F("MC_NumberofClusters", "MC_NumberofClusters", 25, 0 ,25);
          fNumberOfClusters[iCut]->SetXTitle("Number of clusters");
          fTrueJetHistograms[iCut]->Add(fNumberOfClusters[iCut]);
          fNumberOfClustersinJets[iCut]      = new TH1F("MC_NumberofClusters_inJets", "MC_NumberofClusters_inJets", 25, 0 ,25);
          fNumberOfClustersinJets[iCut]->SetXTitle("Number of clusters in Jets");
          fTrueJetHistograms[iCut]->Add(fNumberOfClustersinJets[iCut]);
          fEnergyRatio[iCut]      = new TH2F("MC_EnergyFrac", "MC_EnergyFrac", 25, 0, 1, 30 , 0, 30);
          fEnergyRatio[iCut]->SetXTitle("Energy fraction if NLabels > 1");
          fEnergyRatio[iCut]->SetYTitle("Energy particle");
          fTrueJetHistograms[iCut]->Add(fEnergyRatio[iCut]);
          fEnergyRatioinJets[iCut]      = new TH2F("MC_EnergyFrac_inJets", "MC_EnergyFrac_inJets", 25, 0, 1, 30 , 0, 30);
          fEnergyRatioinJets[iCut]->SetXTitle("Energy fraction if NLabels > 1 in Jets");
          fEnergyRatioinJets[iCut]->SetYTitle("Energy particle");
          fTrueJetHistograms[iCut]->Add(fEnergyRatioinJets[iCut]);
          fEnergyRatioGamma1[iCut]      = new TH2F("MC_EnergyFracGamma", "MC_EnergyFracGamma", 25, 0, 1, 30 , 0, 30);
          fEnergyRatioGamma1[iCut]->SetXTitle("Energy fraction if Label 1 is Gamma (NLabels > 1)");
          fEnergyRatioGamma1[iCut]->SetYTitle("Energy particle");
          fTrueJetHistograms[iCut]->Add(fEnergyRatioGamma1[iCut]);
          fEnergyRatioGamma1inJets[iCut]      = new TH2F("MC_EnergyFracGamma_inJets", "MC_EnergyFracGamma_inJets", 25, 0 ,1, 30 , 0, 30);
          fEnergyRatioGamma1inJets[iCut]->SetXTitle("Energy fraction if Label 1 is Gamma (NLabels > 1) in Jets");
          fEnergyRatioGamma1inJets[iCut]->SetYTitle("Energy particle");
          fTrueJetHistograms[iCut]->Add(fEnergyRatioGamma1inJets[iCut]);
          fEnergyRatioGammaAnywhere[iCut]      = new TH2F("MC_EnergyFracGammaAnywhere", "MC_EnergyFracGammaAnywhere", 25, 0, 1, 30 , 0, 30);
          fEnergyRatioGammaAnywhere[iCut]->SetXTitle("Energy fraction if a Label is Gamma (NLabels > 1)");
          fEnergyRatioGammaAnywhere[iCut]->SetYTitle("Energy particle");
          fTrueJetHistograms[iCut]->Add(fEnergyRatioGammaAnywhere[iCut]);
          fEnergyRatioGammaAnywhereinJets[iCut]      = new TH2F("MC_EnergyFracGammaAnywhere_inJets", "MC_EnergyFracGammaAnywhere_inJets", 25, 0 ,1, 30 , 0, 30);
          fEnergyRatioGammaAnywhereinJets[iCut]->SetXTitle("Energy fraction if a Label is Gamma (NLabels > 1) in Jets");
          fEnergyRatioGammaAnywhereinJets[iCut]->SetYTitle("Energy particle");
          fTrueJetHistograms[iCut]->Add(fEnergyRatioGammaAnywhereinJets[iCut]);
          fEnergyDeposit[iCut]      = new TH2F("MC_EnergyDep", "MC_EnergyDep", 25, 0, 5, 30 , 0, 30);
          fEnergyDeposit[iCut]->SetXTitle("Energy deposit if NLabels > 1");
          fEnergyDeposit[iCut]->SetYTitle("Energy particle");
          fTrueJetHistograms[iCut]->Add(fEnergyDeposit[iCut]);
          fEnergyDepositinJets[iCut]      = new TH2F("MC_EnergyDep_inJets", "MC_EnergyDep_inJets", 25, 0, 5, 30 , 0, 30);
          fEnergyDepositinJets[iCut]->SetXTitle("Energy deposit if NLabels > 1 in Jets");
          fEnergyDepositinJets[iCut]->SetYTitle("Energy particle");
          fTrueJetHistograms[iCut]->Add(fEnergyDepositinJets[iCut]);
          fEnergyDepGamma1[iCut]      = new TH2F("MC_EnergyDepGamma", "MC_EnergyDepGamma", 25, 0, 5, 30 , 0, 30);
          fEnergyDepGamma1[iCut]->SetXTitle("Energy deposit if Label 1 is Gamma (NLabels > 1)");
          fEnergyDepGamma1[iCut]->SetYTitle("Energy particle");
          fTrueJetHistograms[iCut]->Add(fEnergyDepGamma1[iCut]);
          fEnergyDepGamma1inJets[iCut]      = new TH2F("MC_EnergyDepGamma_inJets", "MC_EnergyDepGamma_inJets", 25, 0 ,5, 30 , 0, 30);
          fEnergyDepGamma1inJets[iCut]->SetXTitle("Energy deposit if Label 1 is Gamma (NLabels > 1) in Jets");
          fEnergyDepGamma1inJets[iCut]->SetYTitle("Energy particle");
          fTrueJetHistograms[iCut]->Add(fEnergyDepGamma1inJets[iCut]);
          fEnergyDepGammaAnywhere[iCut]      = new TH2F("MC_EnergyDepGammaAnywhere", "MC_EnergyDepGammaAnywhere", 25, 0, 5, 30 , 0, 30);
          fEnergyDepGammaAnywhere[iCut]->SetXTitle("Energy deposit if a Label is Gamma (NLabels > 1)");
          fEnergyDepGammaAnywhere[iCut]->SetYTitle("Energy particle");
          fTrueJetHistograms[iCut]->Add(fEnergyDepGammaAnywhere[iCut]);
          fEnergyDepGammaAnywhereinJets[iCut]      = new TH2F("MC_EnergyDepGammaAnywhere_inJets", "MC_EnergyDepGammaAnywhere_inJets", 25, 0 ,5, 30 , 0, 30);
          fEnergyDepGammaAnywhereinJets[iCut]->SetXTitle("Energy deposit if a Label is Gamma (NLabels > 1) in Jets");
          fEnergyDepGammaAnywhereinJets[iCut]->SetYTitle("Energy particle");
          fTrueJetHistograms[iCut]->Add(fEnergyDepGammaAnywhereinJets[iCut]);
          fEnergyRatioGamma1Helped[iCut]      = new TH1F("MC_EnergyFracGamma_Helped", "MC_EnergyFracGamma_Helped", 50, 0 ,1);
          fEnergyRatioGamma1Helped[iCut]->SetXTitle("Energy photon if NLabels > 1");
          fTrueJetHistograms[iCut]->Add(fEnergyRatioGamma1Helped[iCut]);
          fEnergyRatioGamma1HelpedinJets[iCut]      = new TH1F("MC_EnergyFracGamma_Helped_inJets", "MC_EnergyFracGamma_Helped_inJets", 50, 0 ,1);
          fEnergyRatioGamma1HelpedinJets[iCut]->SetXTitle("Energy photon if NLabels > 1 in Jets");
          fTrueJetHistograms[iCut]->Add(fEnergyRatioGamma1HelpedinJets[iCut]);
          fClusterEtaPhiJets[iCut]      = new TH2F("MC_Cluster_EtaPhi_Jets", "MC_Cluster_EtaPhi_Jets", 462,0,2*TMath::Pi(),110,-0.7,0.7);
          fClusterEtaPhiJets[iCut]->SetXTitle("#varphi");
          fClusterEtaPhiJets[iCut]->SetYTitle("#eta");
          fTrueJetHistograms[iCut]->Add(fClusterEtaPhiJets[iCut]);
        }
      }
      if(fDoJetQA){
        if(fDoLightOutput){
          fTrueJetHistograms[iCut] = new TList();
          fTrueJetHistograms[iCut]->SetName(Form("%s_%s_%s True Jet histograms", cutstringEvent.Data(), cutstringCalo.Data(), cutstringMeson.Data()));
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
        fHistoClusPhotonBGPt[iCut]          = new TH2F("ESD_TrueClusPhotonBG_Pt", "ESD_TrueClusPhotonBG_Pt", nBinsClusterPt, arrClusPtBinning,10,-0.5,9.5);
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
        fHistoClusPhotonBGPt[iCut]->SetXTitle("p_{T} (GeV/c)");
        fHistoClusPhotonBGPt[iCut]->SetYTitle("source");
        fTrueList[iCut]->Add(fHistoClusPhotonBGPt[iCut]);
        fHistoClusPhotonPlusConvBGPt[iCut]  = new TH2F("ESD_TrueClusPhotonPlusConvBG_Pt", "ESD_TrueClusPhotonPlusConvBG_Pt", nBinsClusterPt, arrClusPtBinning,10,-0.5,9.5);
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
        fHistoClusPhotonPlusConvBGPt[iCut]->SetXTitle("p_{T} (GeV/c)");
        fHistoClusPhotonPlusConvBGPt[iCut]->SetYTitle("source");
        fTrueList[iCut]->Add(fHistoClusPhotonPlusConvBGPt[iCut]);

        if (fDoClusterQA > 1) {
          fHistoClustPhotonElectronBGPtM02[iCut]        = new TH2F("ESD_TrueClusPhotonElectronBG_Pt_M02", "ESD_TrueClusPhotonElectronBG_Pt_M02", nBinsClusterPt, arrClusPtBinning, 100, 0, 1);
          fHistoClustPhotonElectronBGPtM02[iCut]->SetXTitle("p_{T} (GeV/c)");
          fHistoClustPhotonElectronBGPtM02[iCut]->SetYTitle("#sigma_{long}^{2}");
          fTrueList[iCut]->Add(fHistoClustPhotonElectronBGPtM02[iCut]);
          fHistoClustPhotonPionBGPtM02[iCut]            = new TH2F("ESD_TrueClusPhotonPionBG_Pt_M02", "ESD_TrueClusPhotonPionBG_Pt_M02", nBinsClusterPt, arrClusPtBinning, 100, 0, 1);
          fHistoClustPhotonPionBGPtM02[iCut]->SetXTitle("p_{T} (GeV/c)");
          fHistoClustPhotonPionBGPtM02[iCut]->SetYTitle("#sigma_{long}^{2}");
          fTrueList[iCut]->Add(fHistoClustPhotonPionBGPtM02[iCut]);
          fHistoClustPhotonKaonBGPtM02[iCut]            = new TH2F("ESD_TrueClusPhotonKaonBG_Pt_M02", "ESD_TrueClusPhotonKaonBG_Pt_M02", nBinsClusterPt, arrClusPtBinning, 100, 0, 1);
          fHistoClustPhotonKaonBGPtM02[iCut]->SetXTitle("p_{T} (GeV/c)");
          fHistoClustPhotonKaonBGPtM02[iCut]->SetYTitle("#sigma_{long}^{2}");
          fTrueList[iCut]->Add(fHistoClustPhotonKaonBGPtM02[iCut]);
          fHistoClustPhotonK0lBGPtM02[iCut]             = new TH2F("ESD_TrueClusPhotonK0lBG_Pt_M02", "ESD_TrueClusPhotonK0lBG_Pt_M02", nBinsClusterPt, arrClusPtBinning, 100, 0, 1);
          fHistoClustPhotonK0lBGPtM02[iCut]->SetXTitle("p_{T} (GeV/c)");
          fHistoClustPhotonK0lBGPtM02[iCut]->SetYTitle("#sigma_{long}^{2}");
          fTrueList[iCut]->Add(fHistoClustPhotonK0lBGPtM02[iCut]);
          fHistoClustPhotonNeutronBGPtM02[iCut]         = new TH2F("ESD_TrueClusPhotonNeutronBG_Pt_M02", "ESD_TrueClusPhotonNeutronBG_Pt_M02", nBinsClusterPt, arrClusPtBinning, 100, 0, 1);
          fHistoClustPhotonNeutronBGPtM02[iCut]->SetXTitle("p_{T} (GeV/c)");
          fHistoClustPhotonNeutronBGPtM02[iCut]->SetYTitle("#sigma_{long}^{2}");
          fTrueList[iCut]->Add(fHistoClustPhotonNeutronBGPtM02[iCut]);
          fHistoClustPhotonRestBGPtM02[iCut]            = new TH2F("ESD_TrueClusPhotonRestBG_Pt_M02", "ESD_TrueClusPhotonRestBG_Pt_M02", nBinsClusterPt, arrClusPtBinning, 100, 0, 1);
          fHistoClustPhotonRestBGPtM02[iCut]->SetXTitle("p_{T} (GeV/c)");
          fHistoClustPhotonRestBGPtM02[iCut]->SetYTitle("#sigma_{long}^{2}");
          fTrueList[iCut]->Add(fHistoClustPhotonRestBGPtM02[iCut]);
          fHistoClustPhotonPlusConvElectronBGPtM02[iCut]= new TH2F("ESD_TrueClusPhotonPlusConvElectronBG_Pt_M02", "ESD_TrueClusPhotonPlusConvElectronBG_Pt_M02", nBinsClusterPt, arrClusPtBinning, 100, 0, 1);
          fHistoClustPhotonPlusConvElectronBGPtM02[iCut]->SetXTitle("p_{T} (GeV/c)");
          fHistoClustPhotonPlusConvElectronBGPtM02[iCut]->SetYTitle("#sigma_{long}^{2}");
          fTrueList[iCut]->Add(fHistoClustPhotonPlusConvElectronBGPtM02[iCut]);
          fHistoClustPhotonPlusConvPionBGPtM02[iCut]    = new TH2F("ESD_TrueClusPhotonPlusConvPionBG_Pt_M02", "ESD_TrueClusPhotonPlusConvPionBG_Pt_M02", nBinsClusterPt, arrClusPtBinning, 100, 0, 1);
          fHistoClustPhotonPlusConvPionBGPtM02[iCut]->SetXTitle("p_{T} (GeV/c)");
          fHistoClustPhotonPlusConvPionBGPtM02[iCut]->SetYTitle("#sigma_{long}^{2}");
          fTrueList[iCut]->Add(fHistoClustPhotonPlusConvPionBGPtM02[iCut]);
          fHistoClustPhotonPlusConvKaonBGPtM02[iCut]    = new TH2F("ESD_TrueClusPhotonPlusConvKaonBG_Pt_M02", "ESD_TrueClusPhotonPlusConvKaonBG_Pt_M02", nBinsClusterPt, arrClusPtBinning, 100, 0, 1);
          fHistoClustPhotonPlusConvKaonBGPtM02[iCut]->SetXTitle("p_{T} (GeV/c)");
          fHistoClustPhotonPlusConvKaonBGPtM02[iCut]->SetYTitle("#sigma_{long}^{2}");
          fTrueList[iCut]->Add(fHistoClustPhotonPlusConvKaonBGPtM02[iCut]);
          fHistoClustPhotonPlusConvK0lBGPtM02[iCut]     = new TH2F("ESD_TrueClusPhotonPlusConvK0lBG_Pt_M02", "ESD_TrueClusPhotonPlusConvK0lBG_Pt_M02", nBinsClusterPt, arrClusPtBinning, 100, 0, 1);
          fHistoClustPhotonPlusConvK0lBGPtM02[iCut]->SetXTitle("p_{T} (GeV/c)");
          fHistoClustPhotonPlusConvK0lBGPtM02[iCut]->SetYTitle("#sigma_{long}^{2}");
          fTrueList[iCut]->Add(fHistoClustPhotonPlusConvK0lBGPtM02[iCut]);
          fHistoClustPhotonPlusConvNeutronBGPtM02[iCut] = new TH2F("ESD_TrueClusPhotonPlusConvNeutronBG_Pt_M02", "ESD_TrueClusPhotonPlusConvNeutronBG_Pt_M02", nBinsClusterPt, arrClusPtBinning, 100, 0, 1);
          fHistoClustPhotonPlusConvNeutronBGPtM02[iCut]->SetXTitle("p_{T} (GeV/c)");
          fHistoClustPhotonPlusConvNeutronBGPtM02[iCut]->SetYTitle("#sigma_{long}^{2}");
          fTrueList[iCut]->Add(fHistoClustPhotonPlusConvNeutronBGPtM02[iCut]);
          fHistoClustPhotonPlusConvRestBGPtM02[iCut]    = new TH2F("ESD_TrueClusPhotonPlusConvRestBG_Pt_M02", "ESD_TrueClusPhotonPlusConvRestBG_Pt_M02", nBinsClusterPt, arrClusPtBinning, 100, 0, 1);
          fHistoClustPhotonPlusConvRestBGPtM02[iCut]->SetXTitle("p_{T} (GeV/c)");
          fHistoClustPhotonPlusConvRestBGPtM02[iCut]->SetYTitle("#sigma_{long}^{2}");
          fTrueList[iCut]->Add(fHistoClustPhotonPlusConvRestBGPtM02[iCut]);
        }
      }
      fHistoTrueClusGammaPt[iCut]                   = new TH1F("TrueClusGamma_Pt", "ESD_TrueClusGamma_Pt", nBinsClusterPt, arrClusPtBinning);
      fHistoTrueClusGammaPt[iCut]->SetXTitle("p_{T} (GeV/c)");
      fTrueList[iCut]->Add(fHistoTrueClusGammaPt[iCut]);
      if(!fDoLightOutput){
        if (fDoClusterQA > 0) {
          fHistoTrueClusGammaPtM02[iCut]              = new TH2F("TrueClusGamma_Pt_M02", "TrueClusGamma_Pt_M02", nBinsClusterPt, arrClusPtBinning, 100, 0, 1);
          fHistoTrueClusGammaPtM02[iCut]->SetXTitle("p_{T} (GeV/c)");
          fHistoTrueClusGammaPtM02[iCut]->SetYTitle("#sigma_{long}^{2}");
          fTrueList[iCut]->Add(fHistoTrueClusGammaPtM02[iCut]);
        }
        fHistoTruePrimaryClusGammaPt[iCut]            = new TH1F("TruePrimaryClusGamma_Pt", "ESD_TruePrimaryClusGamma_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTruePrimaryClusGammaPt[iCut]->SetXTitle("p_{T} (GeV/c)");
        fTrueList[iCut]->Add(fHistoTruePrimaryClusGammaPt[iCut]);
        fHistoTruePrimaryClusGammaESDPtMCPt[iCut]     = new TH2F("TruePrimaryClusGamma_Pt_MCPt", "ESD_TruePrimaryClusGamma_MCPt", nBinsClusterPt, arrClusPtBinning, nBinsClusterPt, arrClusPtBinning);
        fHistoTruePrimaryClusGammaESDPtMCPt[iCut]->SetXTitle("p_{T} (GeV/c)");
        fHistoTruePrimaryClusGammaESDPtMCPt[iCut]->SetYTitle("p_{T, MC} (GeV/c)");
        fTrueList[iCut]->Add(fHistoTruePrimaryClusGammaESDPtMCPt[iCut]);
        fHistoTruePrimaryClusConvGammaPt[iCut]        = new TH1F("TruePrimaryClusConvGamma_Pt", "ESD_TruePrimaryClusConvGamma_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTruePrimaryClusConvGammaPt[iCut]->SetXTitle("p_{T} (GeV/c)");
        fTrueList[iCut]->Add(fHistoTruePrimaryClusConvGammaPt[iCut]);
        fHistoTruePrimaryClusConvGammaESDPtMCPt[iCut] = new TH2F("TruePrimaryClusConvGamma_Pt_MCPt", "ESD_TruePrimaryClusConvGamma_MCPt", nBinsClusterPt, arrClusPtBinning, nBinsClusterPt, arrClusPtBinning);
        fHistoTruePrimaryClusConvGammaESDPtMCPt[iCut]->SetXTitle("p_{T} (GeV/c)");
        fHistoTruePrimaryClusConvGammaESDPtMCPt[iCut]->SetYTitle("p_{T, MC} (GeV/c)");
        fTrueList[iCut]->Add(fHistoTruePrimaryClusConvGammaESDPtMCPt[iCut]);
        fHistoTrueSecondaryClusGammaPt[iCut]          = new TH2F("ESD_TrueSecondaryClusGamma_Pt", "ESD_TrueSecondaryClusGamma_Pt", nBinsClusterPt, arrClusPtBinning, 5, -0.5, 4.5);
        fHistoTrueSecondaryClusGammaPt[iCut]->SetXTitle("p_{T} (GeV/c)");
        fHistoTrueSecondaryClusGammaPt[iCut]->SetYTitle("source");
        fHistoTrueSecondaryClusGammaPt[iCut]->GetYaxis()->SetBinLabel( 1,"K0s");
        fHistoTrueSecondaryClusGammaPt[iCut]->GetYaxis()->SetBinLabel( 2,"K0l");
        fHistoTrueSecondaryClusGammaPt[iCut]->GetYaxis()->SetBinLabel( 3,"Lambda");
        fHistoTrueSecondaryClusGammaPt[iCut]->GetYaxis()->SetBinLabel( 4,"Eta");
        fHistoTrueSecondaryClusGammaPt[iCut]->GetYaxis()->SetBinLabel( 5,"rest");
        fTrueList[iCut]->Add(fHistoTrueSecondaryClusGammaPt[iCut]);
        fHistoTrueSecondaryClusConvGammaPt[iCut]      = new TH2F("ESD_TrueSecondaryClusConvGamma_Pt", "ESD_TrueSecondaryClusConvGamma_Pt", nBinsClusterPt, arrClusPtBinning, 5, -0.5, 4.5);
        fHistoTrueSecondaryClusConvGammaPt[iCut]->SetXTitle("p_{T} (GeV/c)");
        fHistoTrueSecondaryClusConvGammaPt[iCut]->SetYTitle("source");
        fHistoTrueSecondaryClusConvGammaPt[iCut]->GetYaxis()->SetBinLabel( 1,"K0s");
        fHistoTrueSecondaryClusConvGammaPt[iCut]->GetYaxis()->SetBinLabel( 2,"K0l");
        fHistoTrueSecondaryClusConvGammaPt[iCut]->GetYaxis()->SetBinLabel( 3,"Lambda");
        fHistoTrueSecondaryClusConvGammaPt[iCut]->GetYaxis()->SetBinLabel( 4,"Eta");
        fHistoTrueSecondaryClusConvGammaPt[iCut]->GetYaxis()->SetBinLabel( 5,"rest");
        fTrueList[iCut]->Add(fHistoTrueSecondaryClusConvGammaPt[iCut]);
        fHistoTrueSecondaryClusGammaMCPt[iCut]          = new TH2F("ESD_TrueSecondaryClusGamma_MCPt", "ESD_TrueSecondaryClusGamma_MCPt", nBinsClusterPt, arrClusPtBinning, 5, -0.5, 4.5);
        fHistoTrueSecondaryClusGammaMCPt[iCut]->SetXTitle("p_{T, MC} (GeV/c)");
        fHistoTrueSecondaryClusGammaMCPt[iCut]->SetYTitle("source");
        fHistoTrueSecondaryClusGammaMCPt[iCut]->GetYaxis()->SetBinLabel( 1,"K0s");
        fHistoTrueSecondaryClusGammaMCPt[iCut]->GetYaxis()->SetBinLabel( 2,"K0l");
        fHistoTrueSecondaryClusGammaMCPt[iCut]->GetYaxis()->SetBinLabel( 3,"Lambda");
        fHistoTrueSecondaryClusGammaMCPt[iCut]->GetYaxis()->SetBinLabel( 4,"Eta");
        fHistoTrueSecondaryClusGammaMCPt[iCut]->GetYaxis()->SetBinLabel( 5,"rest");
        fTrueList[iCut]->Add(fHistoTrueSecondaryClusGammaMCPt[iCut]);
        fHistoTrueSecondaryClusConvGammaMCPt[iCut]      = new TH2F("ESD_TrueSecondaryClusConvGamma_MCPt", "ESD_TrueSecondaryClusConvGamma_MCPt", nBinsClusterPt, arrClusPtBinning, 5, -0.5, 4.5);
        fHistoTrueSecondaryClusConvGammaMCPt[iCut]->SetXTitle("p_{T, MC} (GeV/c)");
        fHistoTrueSecondaryClusConvGammaMCPt[iCut]->SetYTitle("source");
        fHistoTrueSecondaryClusConvGammaMCPt[iCut]->GetYaxis()->SetBinLabel( 1,"K0s");
        fHistoTrueSecondaryClusConvGammaMCPt[iCut]->GetYaxis()->SetBinLabel( 2,"K0l");
        fHistoTrueSecondaryClusConvGammaMCPt[iCut]->GetYaxis()->SetBinLabel( 3,"Lambda");
        fHistoTrueSecondaryClusConvGammaMCPt[iCut]->GetYaxis()->SetBinLabel( 4,"Eta");
        fHistoTrueSecondaryClusConvGammaMCPt[iCut]->GetYaxis()->SetBinLabel( 5,"rest");
        fTrueList[iCut]->Add(fHistoTrueSecondaryClusConvGammaMCPt[iCut]);

        fHistoTrueSecondaryClusGammaFromXFromK0sMCPtESDPt[iCut] = new TH2F("ESD_TrueSecondaryClusGammaFromXFromK0s_MCPt_Pt", "ESD_TrueSecondaryClusGammaFromXFromK0s_MCPt_Pt",
                                                                           nBinsClusterPt, arrClusPtBinning, nBinsClusterPt, arrClusPtBinning);
        fHistoTrueSecondaryClusGammaFromXFromK0sMCPtESDPt[iCut]->SetXTitle("p_{T, MC} (GeV/c)");
        fHistoTrueSecondaryClusGammaFromXFromK0sMCPtESDPt[iCut]->SetYTitle("p_{T} (GeV/c)");
        fTrueList[iCut]->Add(fHistoTrueSecondaryClusGammaFromXFromK0sMCPtESDPt[iCut]);
        fHistoTrueSecondaryClusConvGammaFromXFromK0sMCPtESDPt[iCut] = new TH2F("ESD_TrueSecondaryClusConvGammaFromXFromK0s_MCPt_Pt", "ESD_TrueSecondaryClusConvGammaFromXFromK0s_MCPt_Pt",
                                                                               nBinsClusterPt, arrClusPtBinning, nBinsClusterPt, arrClusPtBinning);
        fHistoTrueSecondaryClusConvGammaFromXFromK0sMCPtESDPt[iCut]->SetXTitle("p_{T, MC} (GeV/c)");
        fHistoTrueSecondaryClusConvGammaFromXFromK0sMCPtESDPt[iCut]->SetYTitle("p_{T} (GeV/c)");
        fTrueList[iCut]->Add(fHistoTrueSecondaryClusConvGammaFromXFromK0sMCPtESDPt[iCut]);
        fHistoTrueSecondaryClusGammaFromXFromK0lMCPtESDPt[iCut] = new TH2F("ESD_TrueSecondaryClusGammaFromXFromK0l_MCPt_Pt", "ESD_TrueSecondaryClusGammaFromXFromK0l_MCPt_Pt",
                                                                           nBinsClusterPt, arrClusPtBinning, nBinsClusterPt, arrClusPtBinning);
        fHistoTrueSecondaryClusGammaFromXFromK0lMCPtESDPt[iCut]->SetXTitle("p_{T, MC} (GeV/c)");
        fHistoTrueSecondaryClusGammaFromXFromK0lMCPtESDPt[iCut]->SetYTitle("p_{T} (GeV/c)");
        fTrueList[iCut]->Add(fHistoTrueSecondaryClusGammaFromXFromK0lMCPtESDPt[iCut]);
        fHistoTrueSecondaryClusConvGammaFromXFromK0lMCPtESDPt[iCut] = new TH2F("ESD_TrueSecondaryClusConvGammaFromXFromK0l_MCPt_Pt", "ESD_TrueSecondaryClusConvGammaFromXFromK0l_MCPt_Pt",
                                                                               nBinsClusterPt, arrClusPtBinning, nBinsClusterPt, arrClusPtBinning);
        fHistoTrueSecondaryClusConvGammaFromXFromK0lMCPtESDPt[iCut]->SetXTitle("p_{T, MC} (GeV/c)");
        fHistoTrueSecondaryClusConvGammaFromXFromK0lMCPtESDPt[iCut]->SetYTitle("p_{T} (GeV/c)");
        fTrueList[iCut]->Add(fHistoTrueSecondaryClusConvGammaFromXFromK0lMCPtESDPt[iCut]);
        fHistoTrueSecondaryClusGammaFromXFromLambdaMCPtESDPt[iCut] = new TH2F("ESD_TrueSecondaryClusGammaFromXFromLambda_MCPt_Pt", "ESD_TrueSecondaryClusGammaFromXFromLambda_MCPt_Pt",
                                                                              nBinsClusterPt, arrClusPtBinning, nBinsClusterPt, arrClusPtBinning);
        fHistoTrueSecondaryClusGammaFromXFromLambdaMCPtESDPt[iCut]->SetXTitle("p_{T, MC} (GeV/c)");
        fHistoTrueSecondaryClusGammaFromXFromLambdaMCPtESDPt[iCut]->SetYTitle("p_{T} (GeV/c)");
        fTrueList[iCut]->Add(fHistoTrueSecondaryClusGammaFromXFromLambdaMCPtESDPt[iCut]);
        fHistoTrueSecondaryClusConvGammaFromXFromLambdaMCPtESDPt[iCut] = new TH2F("ESD_TrueSecondaryClusConvGammaFromXFromLambda_MCPt_Pt", "ESD_TrueSecondaryClusConvGammaFromXFromLambda_MCPt_Pt",
                                                                                  nBinsClusterPt, arrClusPtBinning, nBinsClusterPt, arrClusPtBinning);
        fHistoTrueSecondaryClusConvGammaFromXFromLambdaMCPtESDPt[iCut]->SetXTitle("p_{T, MC} (GeV/c)");
        fHistoTrueSecondaryClusConvGammaFromXFromLambdaMCPtESDPt[iCut]->SetYTitle("p_{T} (GeV/c)");
        fTrueList[iCut]->Add(fHistoTrueSecondaryClusConvGammaFromXFromLambdaMCPtESDPt[iCut]);

        fHistoTrueNLabelsInClus[iCut]                 = new TH1F("TrueNLabelsInClus", "TrueNLabelsInClus", 100, -0.5, 99.5);
        fHistoTrueNLabelsInClus[iCut]->SetXTitle("# labels");
        fTrueList[iCut]->Add(fHistoTrueNLabelsInClus[iCut]);
      }
      fHistoDoubleCountTrueClusterGammaPt[iCut]     = new TH2F("TrueDoubleCountClusterGamma_Pt", "TrueDoubleCountClusterGamma_Pt", nBinsClusterPt, arrClusPtBinning, 2, 0, 2);
      fHistoDoubleCountTrueClusterGammaPt[iCut]->SetXTitle("p_{T} (GeV/c)");
      fTrueList[iCut]->Add(fHistoDoubleCountTrueClusterGammaPt[iCut]);
      fHistoMultipleCountTrueClusterGamma[iCut]     = new TH1F("TrueMultipleCountClusterGamma", "TrueMultipleCountClusterGamma", 10, 1, 11);
      fHistoMultipleCountTrueClusterGamma[iCut]->SetXTitle("# multiple");
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

      if (fDoClusterQA > 0){
        fHistoTrueClusUnConvGammaPt[iCut]           = new TH1F("TrueClusUnConvGamma_Pt", "TrueClusUnConvGamma_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTrueClusUnConvGammaPt[iCut]->SetXTitle("p_{T} (GeV/c)");
        fTrueList[iCut]->Add(fHistoTrueClusUnConvGammaPt[iCut]);
        fHistoTrueClusUnConvGammaMCPt[iCut]         = new TH1F("TrueClusUnConvGamma_MCPt", "TrueClusUnConvGamma_MCPt", nBinsClusterPt, arrClusPtBinning);
        fHistoTrueClusUnConvGammaMCPt[iCut]->SetXTitle("p_{T, MC} (GeV/c)");
        fTrueList[iCut]->Add(fHistoTrueClusUnConvGammaMCPt[iCut]);
        if (!fDoLightOutput) {
          fHistoTrueClusUnConvGammaPtM02[iCut]      = new TH2F("TrueClusUnConvGamma_Pt_M02", "TrueClusUnConvGamma_Pt_M02", nBinsClusterPt, arrClusPtBinning, 100, 0, 1);
          fHistoTrueClusUnConvGammaPtM02[iCut]->SetXTitle("p_{T, MC} (GeV/c)");
          fHistoTrueClusUnConvGammaPtM02[iCut]->SetYTitle("#sigma_{long}^{2}");
          fTrueList[iCut]->Add(fHistoTrueClusUnConvGammaPtM02[iCut]);
        }
        fHistoTrueClusElectronPt[iCut]              = new TH1F("TrueClusElectron_Pt", "TrueElectronGamma_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTrueClusElectronPt[iCut]->SetXTitle("p_{T} (GeV/c)");
        fTrueList[iCut]->Add(fHistoTrueClusElectronPt[iCut]);
        fHistoTrueClusConvGammaPt[iCut]             = new TH1F("TrueClusConvGamma_Pt", "TrueClusConvGamma_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTrueClusConvGammaPt[iCut]->SetXTitle("p_{T} (GeV/c)");
        fTrueList[iCut]->Add(fHistoTrueClusConvGammaPt[iCut]);
        fHistoTrueClusConvGammaMCPt[iCut]           = new TH1F("TrueClusConvGamma_MCPt", "TrueClusConvGamma_MCPt", nBinsClusterPt, arrClusPtBinning);
        fHistoTrueClusConvGammaMCPt[iCut]->SetXTitle("p_{T, MC} (GeV/c)");
        fTrueList[iCut]->Add(fHistoTrueClusConvGammaMCPt[iCut]);
        fHistoTrueClusConvGammaFullyPt[iCut]        = new TH1F("TrueClusConvGammaFullyContained_Pt", "TrueClusConvGammaFullyContained_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTrueClusConvGammaFullyPt[iCut]->SetXTitle("p_{T} (GeV/c)");
        fTrueList[iCut]->Add(fHistoTrueClusConvGammaFullyPt[iCut]);
        fHistoTrueClusMergedGammaPt[iCut]           = new TH1F("TrueClusMergedGamma_Pt", "TrueClusMergedGamma_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTrueClusMergedGammaPt[iCut]->SetXTitle("p_{T} (GeV/c)");
        fTrueList[iCut]->Add(fHistoTrueClusMergedGammaPt[iCut]);
        fHistoTrueClusMergedPartConvGammaPt[iCut]   = new TH1F("TrueClusMergedPartConvGamma_Pt", "TrueClusMergedPartConvGamma_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTrueClusMergedPartConvGammaPt[iCut]->SetXTitle("p_{T} (GeV/c)");
        fTrueList[iCut]->Add(fHistoTrueClusMergedPartConvGammaPt[iCut]);
        fHistoTrueClusDalitzPt[iCut]                = new TH1F("TrueClusDalitz_Pt", "TrueClusDalitz_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTrueClusDalitzPt[iCut]->SetXTitle("p_{T} (GeV/c)");
        fTrueList[iCut]->Add(fHistoTrueClusDalitzPt[iCut]);
        fHistoTrueClusDalitzMergedPt[iCut]          = new TH1F("TrueClusDalitzMerged_Pt", "TrueClusDalitzMerged_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTrueClusDalitzMergedPt[iCut]->SetXTitle("p_{T} (GeV/c)");
        fTrueList[iCut]->Add(fHistoTrueClusDalitzMergedPt[iCut]);
        fHistoTrueClusPhotonFromElecMotherPt[iCut]  = new TH1F("TrueClusPhotonFromElecMother_Pt", "TrueClusPhotonFromElecMother_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTrueClusPhotonFromElecMotherPt[iCut]->SetXTitle("p_{T} (GeV/c)");
        fTrueList[iCut]->Add(fHistoTrueClusPhotonFromElecMotherPt[iCut]);
        fHistoTrueClusShowerPt[iCut]                = new TH1F("TrueClusShower_Pt", "TrueClusShower_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTrueClusShowerPt[iCut]->SetXTitle("p_{T} (GeV/c)");
        fTrueList[iCut]->Add(fHistoTrueClusShowerPt[iCut]);
        fHistoTrueClusSubLeadingPt[iCut]            = new TH1F("TrueClusSubleading_Pt", "TrueClusSubleading_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTrueClusSubLeadingPt[iCut]->SetXTitle("p_{T} (GeV/c)");
        fTrueList[iCut]->Add(fHistoTrueClusSubLeadingPt[iCut]);
        fHistoTrueClusNParticles[iCut]              = new TH1F("TrueClusNParticles", "TrueClusNParticles", 20, 0, 20);
        fHistoTrueClusNParticles[iCut]->SetXTitle("#particles");
        fTrueList[iCut]->Add(fHistoTrueClusNParticles[iCut]);
        fHistoTrueClusEMNonLeadingPt[iCut]          = new TH1F("TrueClusEMNonLeading_Pt", "TrueClusEMNonLeading_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTrueClusEMNonLeadingPt[iCut]->SetXTitle("p_{T} (GeV/c)");
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
        fHistoTruePi0InvMassPt[iCut]                    = new TH2F("ESD_TruePi0_InvMass_Pt", "ESD_TruePi0_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fHistoTruePi0InvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
        fHistoTruePi0InvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fTrueList[iCut]->Add(fHistoTruePi0InvMassPt[iCut]);
        if (!fDoLightOutput){
            fHistoTruePi0InvMassPtAdditional[iCut]                = new TH2F("ESD_TruePi0_InvMass_Pt_Additional", "ESD_TruePi0_InvMass_Pt_Additional", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
            fHistoTruePi0InvMassPtAdditional[iCut]->SetXTitle("M_{inv,#pi^{0}}(GeV/c^{2})");
            fHistoTruePi0InvMassPtAdditional[iCut]->SetYTitle("#pi^{0}p_{T}(GeV/c)");
            fTrueList[iCut]->Add(fHistoTruePi0InvMassPtAdditional[iCut]);
        }
        fHistoDoubleCountTruePi0InvMassPt[iCut]         = new TH2F("ESD_TrueDoubleCountPi0_InvMass_Pt", "ESD_TrueDoubleCountPi0_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fHistoDoubleCountTruePi0InvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
        fHistoDoubleCountTruePi0InvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fTrueList[iCut]->Add(fHistoDoubleCountTruePi0InvMassPt[iCut]);
        if( !fDoPi0Only ){
          fHistoTrueEtaInvMassPt[iCut]                    = new TH2F("ESD_TrueEta_InvMass_Pt", "ESD_TrueEta_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
          fHistoTrueEtaInvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
          fHistoTrueEtaInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
          fTrueList[iCut]->Add(fHistoTrueEtaInvMassPt[iCut]);
          if (!fDoLightOutput){
              fHistoTrueEtaInvMassPtAdditional[iCut]                =  new TH2F("ESD_TrueEta_InvMass_Pt_Additional", "ESD_TrueEta_InvMass_PtAdditional", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
              fHistoTrueEtaInvMassPtAdditional[iCut]->SetXTitle("M_{inv,#eta}(GeV/c^{2})");
              fHistoTrueEtaInvMassPtAdditional[iCut]->SetYTitle("#eta p_{T}(GeV/c)");
              fTrueList[iCut]->Add(fHistoTrueEtaInvMassPtAdditional[iCut]);
          }
          fHistoDoubleCountTrueEtaInvMassPt[iCut]         = new TH2F("ESD_TrueDoubleCountEta_InvMass_Pt", "ESD_TrueDoubleCountEta_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
          fHistoDoubleCountTrueEtaInvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
          fHistoDoubleCountTrueEtaInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
          fTrueList[iCut]->Add(fHistoDoubleCountTrueEtaInvMassPt[iCut]);
          fHistoTruePrimaryEtaInvMassPt[iCut]             = new TH2F("ESD_TruePrimaryEta_InvMass_Pt", "ESD_TruePrimaryEta_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
          fHistoTruePrimaryEtaInvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
          fHistoTruePrimaryEtaInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
          fHistoTruePrimaryEtaInvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePrimaryEtaInvMassPt[iCut]);
          fHistoTruePrimaryEtaW0WeightingInvMassPt[iCut]  = new TH2F("ESD_TruePrimaryEtaW0Weights_InvMass_Pt", "ESD_TruePrimaryEtaW0Weights_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
          fHistoTruePrimaryEtaW0WeightingInvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
          fHistoTruePrimaryEtaW0WeightingInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
          fHistoTruePrimaryEtaW0WeightingInvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePrimaryEtaW0WeightingInvMassPt[iCut]);
          fProfileTruePrimaryEtaWeightsInvMassPt[iCut]    = new TProfile2D("ESD_TruePrimaryEtaWeights_InvMass_Pt", "ESD_TruePrimaryEtaWeights_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
          fProfileTruePrimaryEtaWeightsInvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
          fProfileTruePrimaryEtaWeightsInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
          fProfileTruePrimaryEtaWeightsInvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fProfileTruePrimaryEtaWeightsInvMassPt[iCut]);
        }
        fHistoTruePrimaryPi0InvMassPt[iCut]             = new TH2F("ESD_TruePrimaryPi0_InvMass_Pt", "ESD_TruePrimaryPi0_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fHistoTruePrimaryPi0InvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
        fHistoTruePrimaryPi0InvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoTruePrimaryPi0InvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTruePrimaryPi0InvMassPt[iCut]);
        fHistoTruePrimaryPi0W0WeightingInvMassPt[iCut]  = new TH2F("ESD_TruePrimaryPi0W0Weights_InvMass_Pt", "ESD_TruePrimaryPi0W0Weights_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fHistoTruePrimaryPi0W0WeightingInvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
        fHistoTruePrimaryPi0W0WeightingInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoTruePrimaryPi0W0WeightingInvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTruePrimaryPi0W0WeightingInvMassPt[iCut]);
        fProfileTruePrimaryPi0WeightsInvMassPt[iCut]    = new TProfile2D("ESD_TruePrimaryPi0Weights_InvMass_Pt", "ESD_TruePrimaryPi0Weights_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fProfileTruePrimaryPi0WeightsInvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
        fProfileTruePrimaryPi0WeightsInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fProfileTruePrimaryPi0WeightsInvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fProfileTruePrimaryPi0WeightsInvMassPt[iCut]);
        fHistoTrueSecondaryPi0InvMassPt[iCut]           = new TH2F("ESD_TrueSecondaryPi0_InvMass_Pt", "ESD_TrueSecondaryPi0_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fHistoTrueSecondaryPi0InvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
        fHistoTrueSecondaryPi0InvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoTrueSecondaryPi0InvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTrueSecondaryPi0InvMassPt[iCut]);

        if(!fDoLightOutput){
          fHistoTruePi0InvMassPtAlpha[iCut]               = new TH2F("ESD_TruePi0_InvMass_vs_Pt_Alpha", "ESD_TruePi0_InvMass_vs_Pt_Alpha", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
          fHistoTruePi0InvMassPtAlpha[iCut]->SetYTitle("p_{T} (GeV/c)");
          fHistoTruePi0InvMassPtAlpha[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
          fHistoTruePi0InvMassPtAlpha[iCut]->Sumw2();
          fESDList[iCut]->Add(fHistoTruePi0InvMassPtAlpha[iCut]);
          fHistoTruePi0PureGammaInvMassPtAlpha[iCut]      = new TH2F("ESD_TruePi0PureGamma_InvMass_vs_Pt_Alpha", "ESD_TruePi0PureGamma_InvMass_vs_Pt_Alpha", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
          fHistoTruePi0PureGammaInvMassPtAlpha[iCut]->SetYTitle("p_{T} (GeV/c)");
          fHistoTruePi0PureGammaInvMassPtAlpha[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
          fHistoTruePi0PureGammaInvMassPtAlpha[iCut]->Sumw2();
          fESDList[iCut]->Add(fHistoTruePi0PureGammaInvMassPtAlpha[iCut]);
        }
        fHistoTrueSecondaryPi0FromK0sInvMassPt[iCut]    = new TH2F("ESD_TrueSecondaryPi0FromK0s_InvMass_Pt", "ESD_TrueSecondaryPi0FromK0s_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fHistoTrueSecondaryPi0FromK0sInvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
        fHistoTrueSecondaryPi0FromK0sInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoTrueSecondaryPi0FromK0sInvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTrueSecondaryPi0FromK0sInvMassPt[iCut]);
        fHistoTrueSecondaryPi0FromK0lInvMassPt[iCut]    = new TH2F("ESD_TrueSecondaryPi0FromK0l_InvMass_Pt", "ESD_TrueSecondaryPi0FromK0l_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fHistoTrueSecondaryPi0FromK0lInvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
        fHistoTrueSecondaryPi0FromK0lInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoTrueSecondaryPi0FromK0lInvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTrueSecondaryPi0FromK0lInvMassPt[iCut]);
        fHistoTrueSecondaryPi0FromEtaInvMassPt[iCut]    = new TH2F("ESD_TrueSecondaryPi0FromEta_InvMass_Pt", "ESD_TrueSecondaryPi0FromEta_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fHistoTrueSecondaryPi0FromEtaInvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
        fHistoTrueSecondaryPi0FromEtaInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fTrueList[iCut]->Add(fHistoTrueSecondaryPi0FromEtaInvMassPt[iCut]);
        fHistoTrueSecondaryPi0FromLambdaInvMassPt[iCut] = new TH2F("ESD_TrueSecondaryPi0FromLambda_InvMass_Pt", "ESD_TrueSecondaryPi0FromLambda_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fHistoTrueSecondaryPi0FromLambdaInvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
        fHistoTrueSecondaryPi0FromLambdaInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fTrueList[iCut]->Add(fHistoTrueSecondaryPi0FromLambdaInvMassPt[iCut]);
        if (fIsMC > 1){
          fHistoTruePi0InvMassPt[iCut]->Sumw2();
          fHistoDoubleCountTruePi0InvMassPt[iCut]->Sumw2();
          fHistoTruePrimaryPi0W0WeightingInvMassPt[iCut]->Sumw2();
          fHistoTrueSecondaryPi0FromEtaInvMassPt[iCut]->Sumw2();
          fHistoTrueSecondaryPi0FromLambdaInvMassPt[iCut]->Sumw2();
          if (!fDoLightOutput){
              fHistoTruePi0InvMassPtAdditional[iCut]->Sumw2();
          }
          if( !fDoPi0Only ){
            fHistoTrueEtaInvMassPt[iCut]->Sumw2();
            if (!fDoLightOutput){
                fHistoTrueEtaInvMassPtAdditional[iCut]->Sumw2();
            }
            fHistoDoubleCountTrueEtaInvMassPt[iCut]->Sumw2();
            fHistoTruePrimaryEtaW0WeightingInvMassPt[iCut]->Sumw2();
          }
        }

        if (fDoMesonQA > 0 && fDoMesonQA < 3){
          if (fIsMC < 2){
            fHistoTruePi0CaloPhotonInvMassPt[iCut]      = new TH2F("ESD_TruePi0CaloPhoton_InvMass_Pt", "ESD_TruePi0CaloPhoton_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
            fHistoTruePi0CaloPhotonInvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
            fHistoTruePi0CaloPhotonInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
            fTrueList[iCut]->Add(fHistoTruePi0CaloPhotonInvMassPt[iCut]);
            if( !fDoPi0Only ){
              fHistoTrueEtaCaloPhotonInvMassPt[iCut]      = new TH2F("ESD_TrueEtaCaloPhoton_InvMass_Pt", "ESD_TrueEtaCaloPhoton_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
              fHistoTrueEtaCaloPhotonInvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
              fHistoTrueEtaCaloPhotonInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
              fTrueList[iCut]->Add(fHistoTrueEtaCaloPhotonInvMassPt[iCut]);
              fHistoTrueEtaCaloMixedPhotonConvPhotonInvMassPt[iCut]   = new TH2F("ESD_TrueEtaCaloMixedPhotonConvertedPhoton_InvMass_Pt", "ESD_TrueEtaCaloMixedPhotonConvertedPhoton_InvMass_Pt", nBinsMinv, 0, maxMinv,
                                                                                 nBinsPt, arrPtBinning);
              fHistoTrueEtaCaloMixedPhotonConvPhotonInvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
              fHistoTrueEtaCaloMixedPhotonConvPhotonInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
              fTrueList[iCut]->Add(fHistoTrueEtaCaloMixedPhotonConvPhotonInvMassPt[iCut]);
              fHistoTrueEtaCaloConvertedPhotonInvMassPt[iCut]         = new TH2F("ESD_TrueEtaCaloConvertedPhoton_InvMass_Pt", "ESD_TrueEtaCaloConvertedPhoton_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
              fHistoTrueEtaCaloConvertedPhotonInvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
              fHistoTrueEtaCaloConvertedPhotonInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
              fTrueList[iCut]->Add(fHistoTrueEtaCaloConvertedPhotonInvMassPt[iCut]);
              fHistoTrueEtaCaloElectronInvMassPt[iCut]    = new TH2F("ESD_TrueEtaCaloElectron_InvMass_Pt", "ESD_TrueEtaCaloElectron_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
              fHistoTrueEtaCaloElectronInvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
              fHistoTrueEtaCaloElectronInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
              fTrueList[iCut]->Add(fHistoTrueEtaCaloElectronInvMassPt[iCut]);
              fHistoTrueEtaCaloMergedClusterInvMassPt[iCut]           = new TH2F("ESD_TrueEtaCaloMergedCluster_InvMass_Pt", "ESD_TrueEtaCaloMergedCluster_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
              fHistoTrueEtaCaloMergedClusterInvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
              fHistoTrueEtaCaloMergedClusterInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
              fTrueList[iCut]->Add(fHistoTrueEtaCaloMergedClusterInvMassPt[iCut]);
              fHistoTrueEtaCaloMergedClusterPartConvInvMassPt[iCut]   = new TH2F("ESD_TrueEtaCaloMergedClusterPartConv_InvMass_Pt", "ESD_TrueEtaCaloMergedClusterPartConv_InvMass_Pt", nBinsMinv, 0, maxMinv,
                                                                                 nBinsPt, arrPtBinning);
              fHistoTrueEtaCaloMergedClusterPartConvInvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
              fHistoTrueEtaCaloMergedClusterPartConvInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
              fTrueList[iCut]->Add(fHistoTrueEtaCaloMergedClusterPartConvInvMassPt[iCut]);
              fHistoTruePrimaryEtaMCPtResolPt[iCut]       = new TH2F("ESD_TruePrimaryEta_MCPt_ResolPt", "ESD_TruePrimaryEta_ResolPt_MCPt", 500, 0.03, 35, 1000, -1., 1.);
              fHistoTruePrimaryEtaMCPtResolPt[iCut]->SetXTitle("p_{T, MC} (GeV/c)");
              fHistoTruePrimaryEtaMCPtResolPt[iCut]->SetXTitle("#delta (p_{T}) (GeV/c)");
              fHistoTruePrimaryEtaMCPtResolPt[iCut]->Sumw2();
              SetLogBinningXTH2(fHistoTruePrimaryEtaMCPtResolPt[iCut]);
              fTrueList[iCut]->Add(fHistoTruePrimaryEtaMCPtResolPt[iCut]);
            }
            fHistoTruePi0CaloMixedPhotonConvPhotonInvMassPt[iCut]   = new TH2F("ESD_TruePi0CaloMixedPhotonConvertedPhoton_InvMass_Pt", "ESD_TruePi0CaloMixedPhotonConvertedPhoton_InvMass_Pt", nBinsMinv, 0, maxMinv,
                                                                               nBinsPt, arrPtBinning);
            fHistoTruePi0CaloMixedPhotonConvPhotonInvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
            fHistoTruePi0CaloMixedPhotonConvPhotonInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
            fTrueList[iCut]->Add(fHistoTruePi0CaloMixedPhotonConvPhotonInvMassPt[iCut]);
            fHistoTruePi0CaloConvertedPhotonInvMassPt[iCut]         = new TH2F("ESD_TruePi0CaloConvertedPhoton_InvMass_Pt", "ESD_TruePi0CaloConvertedPhoton_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
            fHistoTruePi0CaloConvertedPhotonInvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
            fHistoTruePi0CaloConvertedPhotonInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
            fTrueList[iCut]->Add(fHistoTruePi0CaloConvertedPhotonInvMassPt[iCut]);
            fHistoTruePi0CaloElectronInvMassPt[iCut]    = new TH2F("ESD_TruePi0CaloElectron_InvMass_Pt", "ESD_TruePi0CaloElectron_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
            fHistoTruePi0CaloElectronInvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
            fHistoTruePi0CaloElectronInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
            fTrueList[iCut]->Add(fHistoTruePi0CaloElectronInvMassPt[iCut]);
            fHistoTruePi0CaloMergedClusterInvMassPt[iCut]           = new TH2F("ESD_TruePi0CaloMergedCluster_InvMass_Pt", "ESD_TruePi0CaloMergedCluster_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
            fHistoTruePi0CaloMergedClusterInvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
            fHistoTruePi0CaloMergedClusterInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
            fTrueList[iCut]->Add(fHistoTruePi0CaloMergedClusterInvMassPt[iCut]);
            fHistoTruePi0CaloMergedClusterPartConvInvMassPt[iCut]   = new TH2F("ESD_TruePi0CaloMergedClusterPartConv_InvMass_Pt", "ESD_TruePi0CaloMergedClusterPartConv_InvMass_Pt", nBinsMinv, 0, maxMinv,
                                                                               nBinsPt, arrPtBinning);
            fHistoTruePi0CaloMergedClusterPartConvInvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
            fHistoTruePi0CaloMergedClusterPartConvInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
            fTrueList[iCut]->Add(fHistoTruePi0CaloMergedClusterPartConvInvMassPt[iCut]);
            fHistoTruePi0NonMergedElectronPhotonInvMassPt[iCut]     = new TH2F("ESD_TruePi0NonMergedElectronPhoton_InvMass_Pt", "ESD_TruePi0NonMergedElectronPhoton_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
            fHistoTruePi0NonMergedElectronPhotonInvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
            fHistoTruePi0NonMergedElectronPhotonInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
            fTrueList[iCut]->Add(fHistoTruePi0NonMergedElectronPhotonInvMassPt[iCut]);
            fHistoTruePi0NonMergedElectronMergedPhotonInvMassPt[iCut] = new TH2F("ESD_TruePi0NonMergedElectronMergedPhoton_InvMass_Pt", "ESD_TruePi0NonMergedElectronMergedPhoton_InvMass_Pt", nBinsMinv, 0, maxMinv,
                                                                                 nBinsPt, arrPtBinning);
            fHistoTruePi0NonMergedElectronMergedPhotonInvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
            fHistoTruePi0NonMergedElectronMergedPhotonInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
            fTrueList[iCut]->Add(fHistoTruePi0NonMergedElectronMergedPhotonInvMassPt[iCut]);
            fHistoTruePrimaryPi0MCPtResolPt[iCut]       = new TH2F("ESD_TruePrimaryPi0_MCPt_ResolPt", "ESD_TruePrimaryPi0_ResolPt_MCPt", 500, 0.03, 35, 1000, -1., 1.);
            fHistoTruePrimaryPi0MCPtResolPt[iCut]->SetXTitle("p_{T, MC} (GeV/c)");
            fHistoTruePrimaryPi0MCPtResolPt[iCut]->SetXTitle("#delta (p_{T}) (GeV/c)");
            fHistoTruePrimaryPi0MCPtResolPt[iCut]->Sumw2();
            SetLogBinningXTH2(fHistoTruePrimaryPi0MCPtResolPt[iCut]);
            fTrueList[iCut]->Add(fHistoTruePrimaryPi0MCPtResolPt[iCut]);
            fHistoTrueK0sWithPi0DaughterMCPt[iCut]      = new TH1F("ESD_TrueK0sWithPi0Daughter_MCPt", "ESD_TrueK0sWithPi0Daughter_MCPt", nBinsPt, arrPtBinning);
            fHistoTrueK0sWithPi0DaughterMCPt[iCut]->SetXTitle("p_{T, MC} (GeV/c)");
            fTrueList[iCut]->Add(fHistoTrueK0sWithPi0DaughterMCPt[iCut]);
            fHistoTrueK0lWithPi0DaughterMCPt[iCut]      = new TH1F("ESD_TrueK0lWithPi0Daughter_MCPt", "ESD_TrueK0lWithPi0Daughter_MCPt", nBinsPt, arrPtBinning);
            fHistoTrueK0lWithPi0DaughterMCPt[iCut]->SetXTitle("p_{T, MC} (GeV/c)");
            fTrueList[iCut]->Add(fHistoTrueK0lWithPi0DaughterMCPt[iCut]);
            fHistoTrueEtaWithPi0DaughterMCPt[iCut]      = new TH1F("ESD_TrueEtaWithPi0Daughter_MCPt", "ESD_TrueEtaWithPi0Daughter_MCPt", nBinsPt, arrPtBinning);
            fHistoTrueEtaWithPi0DaughterMCPt[iCut]->SetXTitle("p_{T, MC} (GeV/c)");
            fTrueList[iCut]->Add(fHistoTrueEtaWithPi0DaughterMCPt[iCut]);
            fHistoTrueLambdaWithPi0DaughterMCPt[iCut]   = new TH1F("ESD_TrueLambdaWithPi0Daughter_MCPt", "ESD_TrueLambdaWithPi0Daughter_MCPt", nBinsPt, arrPtBinning);
            fHistoTrueLambdaWithPi0DaughterMCPt[iCut]->SetXTitle("p_{T, MC} (GeV/c)");
            fTrueList[iCut]->Add(fHistoTrueLambdaWithPi0DaughterMCPt[iCut]);
          }

          fHistoTruePi0PtY[iCut]          = new TH2F("ESD_TruePi0_Pt_Y", "ESD_TruePi0_Pt_Y", nBinsQAPt, arrQAPtBinning, 150, -1.5, 1.5);
          fHistoTruePi0PtY[iCut]->SetXTitle("p_{T} (GeV/c)");
          fHistoTruePi0PtY[iCut]->SetYTitle("y");
          fTrueList[iCut]->Add(fHistoTruePi0PtY[iCut]);
          if( !fDoPi0Only ){
            fHistoTrueEtaPtY[iCut]          = new TH2F("ESD_TrueEta_Pt_Y", "ESD_TrueEta_Pt_Y", nBinsQAPt, arrQAPtBinning, 150, -1.5, 1.5);
            fHistoTrueEtaPtY[iCut]->SetXTitle("p_{T} (GeV/c)");
            fHistoTrueEtaPtY[iCut]->SetYTitle("y");
            fTrueList[iCut]->Add(fHistoTrueEtaPtY[iCut]);
            fHistoTrueEtaPtAlpha[iCut]      = new TH2F("ESD_TrueEta_Pt_Alpha", "ESD_TrueEta_Pt_Alpha", nBinsQAPt, arrQAPtBinning, 100, 0, 1);
            fHistoTrueEtaPtAlpha[iCut]->SetXTitle("p_{T} (GeV/c)");
            fHistoTrueEtaPtAlpha[iCut]->SetYTitle("#alpha");
            fTrueList[iCut]->Add(fHistoTrueEtaPtAlpha[iCut]);
            fHistoTrueEtaPtOpenAngle[iCut]  = new TH2F("ESD_TrueEta_Pt_OpenAngle", "ESD_TrueEta_Pt_OpenAngle", nBinsQAPt, arrQAPtBinning, 180, 0, 1.8);
            fHistoTrueEtaPtOpenAngle[iCut]->SetXTitle("p_{T} (GeV/c)");
            fHistoTrueEtaPtOpenAngle[iCut]->SetYTitle("#theta");
            fTrueList[iCut]->Add(fHistoTrueEtaPtOpenAngle[iCut]);
          }
          fHistoTruePi0PtAlpha[iCut]      = new TH2F("ESD_TruePi0_Pt_Alpha", "ESD_TruePi0_Pt_Alpha", nBinsQAPt, arrQAPtBinning, 100, 0, 1);
          fHistoTruePi0PtAlpha[iCut]->SetXTitle("p_{T} (GeV/c)");
          fHistoTruePi0PtAlpha[iCut]->SetYTitle("#alpha");
          fTrueList[iCut]->Add(fHistoTruePi0PtAlpha[iCut]);

          fHistoTruePi0PtOpenAngle[iCut]  = new TH2F("ESD_TruePi0_Pt_OpenAngle", "ESD_TruePi0_Pt_OpenAngle", nBinsQAPt, arrQAPtBinning, 100, 0, 0.5);
          fHistoTruePi0PtOpenAngle[iCut]->SetXTitle("p_{T} (GeV/c)");
          fHistoTruePi0PtOpenAngle[iCut]->SetYTitle("#theta");
          fTrueList[iCut]->Add(fHistoTruePi0PtOpenAngle[iCut]);
          fHistoTrueBckGGInvMassPt[iCut]              = new TH2F("ESD_TrueBckGG_InvMass_Pt", "ESD_TrueBckGG_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
          fHistoTrueBckGGInvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
          fHistoTrueBckGGInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
          fTrueList[iCut]->Add(fHistoTrueBckGGInvMassPt[iCut]);
          fHistoTrueBckFullMesonContainedInOneClusterInvMassPt[iCut] = new TH2F("ESD_TrueBckFullMesonContained_InvMass_Pt", "ESD_TrueBckFullMesonContained_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
          fHistoTrueBckFullMesonContainedInOneClusterInvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
          fHistoTrueBckFullMesonContainedInOneClusterInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
          fTrueList[iCut]->Add(fHistoTrueBckFullMesonContainedInOneClusterInvMassPt[iCut]);
          fHistoTrueBckAsymEClustersInvMassPt[iCut]   = new TH2F("ESD_TrueBckAsymEClus_InvMass_Pt", "ESD_TrueBckAsymEClus_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
          fHistoTrueBckAsymEClustersInvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
          fHistoTrueBckAsymEClustersInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
          fTrueList[iCut]->Add(fHistoTrueBckAsymEClustersInvMassPt[iCut]);
          fHistoTrueBckContInvMassPt[iCut]            = new TH2F("ESD_TrueBckCont_InvMass_Pt", "ESD_TrueBckCont_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
          fHistoTrueBckContInvMassPt[iCut]->SetYTitle("p_{T} (GeV/c)");
          fHistoTrueBckContInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
          fTrueList[iCut]->Add(fHistoTrueBckContInvMassPt[iCut]);

          if (fIsMC > 1){
            fHistoTruePi0PtY[iCut]->Sumw2();
            fHistoTruePi0PtAlpha[iCut]->Sumw2();
            fHistoTruePi0PtOpenAngle[iCut]->Sumw2();
            fHistoTrueBckGGInvMassPt[iCut]->Sumw2();
            fHistoTrueBckFullMesonContainedInOneClusterInvMassPt[iCut]->Sumw2();
            fHistoTrueBckAsymEClustersInvMassPt[iCut]->Sumw2();
            fHistoTrueBckContInvMassPt[iCut]->Sumw2();
            if( !fDoPi0Only ){
              fHistoTrueEtaPtOpenAngle[iCut]->Sumw2();
              fHistoTrueEtaPtAlpha[iCut]->Sumw2();
              fHistoTrueEtaPtY[iCut]->Sumw2();
            }
          }

        }

        if (fDoMesonQA == 2 && fIsMC < 2){
          fHistoTruePi0Category1[iCut]    = new TH2F("ESD_TruePi0Category1_InvMass_Pt", "ESD_TruePi0Category1_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
          fHistoTruePi0Category1[iCut]->SetYTitle("p_{T} (GeV/c)");
          fHistoTruePi0Category1[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
          fTrueList[iCut]->Add(fHistoTruePi0Category1[iCut]);
          fHistoTruePi0Category2[iCut]    = new TH2F("ESD_TruePi0Category2_InvMass_Pt", "ESD_TruePi0Category2_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
          fHistoTruePi0Category2[iCut]->SetYTitle("p_{T} (GeV/c)");
          fHistoTruePi0Category2[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
          fTrueList[iCut]->Add(fHistoTruePi0Category2[iCut]);
          fHistoTruePi0Category3[iCut]    = new TH2F("ESD_TruePi0Category3_InvMass_Pt", "ESD_TruePi0Category3_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
          fHistoTruePi0Category3[iCut]->SetYTitle("p_{T} (GeV/c)");
          fHistoTruePi0Category3[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
          fTrueList[iCut]->Add(fHistoTruePi0Category3[iCut]);
          fHistoTruePi0Category4_6[iCut]  = new TH2F("ESD_TruePi0Category4_6_InvMass_Pt", "ESD_TruePi0Category4_6_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
          fHistoTruePi0Category4_6[iCut]->SetYTitle("p_{T} (GeV/c)");
          fHistoTruePi0Category4_6[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
          fTrueList[iCut]->Add(fHistoTruePi0Category4_6[iCut]);
          fHistoTruePi0Category5[iCut]    = new TH2F("ESD_TruePi0Category5_InvMass_Pt", "ESD_TruePi0Category5_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
          fHistoTruePi0Category5[iCut]->SetYTitle("p_{T} (GeV/c)");
          fHistoTruePi0Category5[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
          fTrueList[iCut]->Add(fHistoTruePi0Category5[iCut]);
          fHistoTruePi0Category7[iCut]    = new TH2F("ESD_TruePi0Category7_InvMass_Pt", "ESD_TruePi0Category7_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
          fHistoTruePi0Category7[iCut]->SetYTitle("p_{T} (GeV/c)");
          fHistoTruePi0Category7[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
          fTrueList[iCut]->Add(fHistoTruePi0Category7[iCut]);
          fHistoTruePi0Category8[iCut]    = new TH2F("ESD_TruePi0Category8_InvMass_Pt", "ESD_TruePi0Category8_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
          fHistoTruePi0Category8[iCut]->SetYTitle("p_{T} (GeV/c)");
          fHistoTruePi0Category8[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
          fTrueList[iCut]->Add(fHistoTruePi0Category8[iCut]);

          if( !fDoPi0Only ){
            fHistoTrueEtaCategory1[iCut]    = new TH2F("ESD_TrueEtaCategory1_InvMass_Pt", "ESD_TrueEtaCategory1_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
            fHistoTrueEtaCategory1[iCut]->SetYTitle("p_{T} (GeV/c)");
            fHistoTrueEtaCategory1[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
            fTrueList[iCut]->Add(fHistoTrueEtaCategory1[iCut]);
            fHistoTrueEtaCategory2[iCut]    = new TH2F("ESD_TrueEtaCategory2_InvMass_Pt", "ESD_TrueEtaCategory2_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
            fHistoTrueEtaCategory2[iCut]->SetYTitle("p_{T} (GeV/c)");
            fHistoTrueEtaCategory2[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
            fTrueList[iCut]->Add(fHistoTrueEtaCategory2[iCut]);
            fHistoTrueEtaCategory3[iCut]    = new TH2F("ESD_TrueEtaCategory3_InvMass_Pt", "ESD_TrueEtaCategory3_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
            fHistoTrueEtaCategory3[iCut]->SetYTitle("p_{T} (GeV/c)");
            fHistoTrueEtaCategory3[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
            fTrueList[iCut]->Add(fHistoTrueEtaCategory3[iCut]);
            fHistoTrueEtaCategory4_6[iCut]  = new TH2F("ESD_TrueEtaCategory4_6_InvMass_Pt", "ESD_TrueEtaCategory4_6_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
            fHistoTrueEtaCategory4_6[iCut]->SetYTitle("p_{T} (GeV/c)");
            fHistoTrueEtaCategory4_6[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
            fTrueList[iCut]->Add(fHistoTrueEtaCategory4_6[iCut]);
            fHistoTrueEtaCategory5[iCut]    = new TH2F("ESD_TrueEtaCategory5_InvMass_Pt", "ESD_TrueEtaCategory5_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
            fHistoTrueEtaCategory5[iCut]->SetYTitle("p_{T} (GeV/c)");
            fHistoTrueEtaCategory5[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
            fTrueList[iCut]->Add(fHistoTrueEtaCategory5[iCut]);
            fHistoTrueEtaCategory7[iCut]    = new TH2F("ESD_TrueEtaCategory7_InvMass_Pt", "ESD_TrueEtaCategory7_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
            fHistoTrueEtaCategory7[iCut]->SetYTitle("p_{T} (GeV/c)");
            fHistoTrueEtaCategory7[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
            fTrueList[iCut]->Add(fHistoTrueEtaCategory7[iCut]);
            fHistoTrueEtaCategory8[iCut]    = new TH2F("ESD_TrueEtaCategory8_InvMass_Pt", "ESD_TrueEtaCategory8_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
            fHistoTrueEtaCategory8[iCut]->SetYTitle("p_{T} (GeV/c)");
            fHistoTrueEtaCategory8[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
            fTrueList[iCut]->Add(fHistoTrueEtaCategory8[iCut]);
          }
        }

        if (fDoMesonQA == 3){
          fTreeList[iCut] = new TList();
          fTreeList[iCut]->SetName(Form("%s_%s_%s True ClusterComb tree", cutstringEvent.Data(), cutstringCalo.Data(), cutstringMeson.Data()));
          fTreeList[iCut]->SetOwner(kTRUE);
          fCutFolder[iCut]->Add(fTreeList[iCut]);

          tTrueInvMassROpenABPtFlag[iCut] = new TTree("True_InvMass_R_OpenA_OpenB_Pt_Flag", "True_InvMass_R_OpenA_OpenB_Pt_Flag");
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

  if (fIsMC > 0 ){
    tBrokenFiles = new TTree("BrokenFiles", "BrokenFiles");
    tBrokenFiles->Branch("fileName",&fFileNameBroken);
    fOutputContainer->Add(tBrokenFiles);
  }

  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
      if ( (fIsMC == 0 ) //Only Data
           &&(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetClusterType()==2) //Only PHOS
           &&(fDoClusterQA > 1) //If QA Flag is set
           &&((((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger()==1)||(((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger()==0)) //Only MB
           ){
          if (tTriggerFiles_wL0==NULL){
            tTriggerFiles_wL0 = new TTree("TriggerFiles_wL0", "TriggerFiles_wL0");
            tTriggerFiles_wL0->Branch("fileName_TriggerFiles_wL0",&fFileNameTrigger);
            fOutputContainer->Add(tTriggerFiles_wL0);
          }
          if (tTriggerFiles_woL0==NULL){
            tTriggerFiles_woL0 = new TTree("TriggerFiles_woL0", "TriggerFiles_woL0");
            tTriggerFiles_woL0->Branch("fileName_TriggerFiles_woL0",&fFileNameTrigger);
            fOutputContainer->Add(tTriggerFiles_woL0);
          }
      }
  }
  if (fDoClusterQA > 1){
    tClusterQATree = new TTree("ClusterQATree", "ClusterQATree");
    tClusterQATree->Branch("closeHighPtClusters",&fCloseHighPtClusters);
    fOutputContainer->Add(tClusterQATree);
  }

  if(fLocalDebugFlag > 0){
    fstream fOutputLocalDebug;
    fOutputLocalDebug.open("debugOutput.txt",ios::out);
    fOutputLocalDebug.close();
  }

  OpenFile(1);
  PostData(1, fOutputContainer);
  Int_t nContainerOutput = 2;
  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
      if(fDoMesonQA == 5 && fIsMC == 0){
          OpenFile(nContainerOutput);
          PostData(nContainerOutput, tClusterTimingEff[iCut]);
          nContainerOutput++;
      }
  }
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
  fInputEvent           = InputEvent();
  fCloseHighPtClusters  = 0x0;
  if(fIsMC> 0) fMCEvent = MCEvent();
  Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();

  if(fInputEvent->IsIncompleteDAQ()==kTRUE) eventQuality = 2;  // incomplete event
  if(eventQuality == 2 || eventQuality == 3){// Event Not Accepted due to MC event missing or wrong trigger for V0ReaderV1 or because it is incomplete
    // write out name of broken file for first event
    if (fIsMC > 0){
      if (fInputEvent->IsA()==AliESDEvent::Class()){
        if (((AliESDEvent*)fInputEvent)->GetEventNumberInFile() == 0){
          fFileNameBroken = new TObjString(Form("%s", ((TString)fV0Reader->GetCurrentFileName()).Data()));
          if (tBrokenFiles) tBrokenFiles->Fill();
          delete fFileNameBroken;
        }
      }
    }

    Bool_t TriggerFiles_Filled=kFALSE;
    Int_t eventHasL0Flag_TrigFiles=(fInputHandler->IsEventSelected() & AliVEvent::kPHI7);
    for(Int_t iCut = 0; iCut<fnCuts; iCut++){
        if (!TriggerFiles_Filled){
            if ((((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetClusterType()==2)){//Only PHOS
                if ( (fIsMC == 0 ) //Only Data
                     && fCaloTriggerMimicHelper[fiCut]
                     &&(fDoClusterQA > 1) //If QA Flag is set
                     &&((((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger()==1)||(((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger()==0)) //Only MB
                     ){
                    if (fCaloTriggerMimicHelper[fiCut]->GetEventChosenByTrigger()){
                        fFileNameTrigger = new TObjString(Form("%s", ((TString)fV0Reader->GetCurrentFileName()).Data()));
                        if (eventHasL0Flag_TrigFiles){
                            if (tTriggerFiles_wL0) tTriggerFiles_wL0->Fill();
                        } else {
                            if (tTriggerFiles_woL0) tTriggerFiles_woL0->Fill();
                        }
                        delete fFileNameTrigger;
                        TriggerFiles_Filled=kTRUE;
                    }
                }
            }
        }
    }

    for(Int_t iCut = 0; iCut<fnCuts; iCut++){
      fHistoNEvents[iCut]->Fill(eventQuality);
      if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventQuality);
    }
    return;
  }
  if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetDoSecondaryTrackMatching())fReaderGammas = fV0Reader->GetReconstructedGammas(); // Gammas from default Cut

  // ------------------- BeginEvent ----------------------------

  AliEventplane *EventPlane = fInputEvent->GetEventplane();
  if(fIsHeavyIon ==1)fEventPlaneAngle = EventPlane->GetEventplane("V0",fInputEvent,2);
  else fEventPlaneAngle=0.0;
  for(Int_t iCut = 0; iCut<fnCuts; iCut++){

    fiCut = iCut;
    fNCurrentClusterBasic       = 0;
    Bool_t isRunningEMCALrelAna = kFALSE;
    if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 1) isRunningEMCALrelAna = kTRUE;

    Int_t eventNotAccepted = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon, isRunningEMCALrelAna);

    if(fIsMC==2){
      Float_t xsection      = -1.;
      Float_t ntrials       = -1.;
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetXSectionAndNTrials(fMCEvent,xsection,ntrials, fInputEvent );
      if((xsection==-1.) || (ntrials==-1.)) AliFatal(Form("ERROR: GetXSectionAndNTrials returned invalid xsection/ntrials, periodName from V0Reader: '%s'",fV0Reader->GetPeriodName().Data()));
      fProfileJetJetXSection[iCut]->Fill(0.,xsection);
      fHistoJetJetNTrials[iCut]->Fill("#sum{NTrials}", ntrials);
    }

    if (fIsMC > 0){
      fWeightJetJetMC       = 1;
      Float_t maxjetpt      = -1.;
      Float_t pthard = -1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetUseJetFinderForOutliers()) maxjetpt = fOutlierJetReader->GetMaxJetPt();
      Bool_t isMCJet        = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsJetJetMCEventAccepted( fMCEvent, fWeightJetJetMC ,pthard, fInputEvent, maxjetpt);
      if(isMCJet && (fIsMC==2))           fHistoPtHardJJWeight[iCut]->Fill(pthard,fWeightJetJetMC);
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
      if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventQuality); // Should be 0 here
      continue;
    }
    if (triggered == kTRUE) {
      if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetUseSphericity()!=0){
        if(fV0Reader->GetSphericity() != -1 && fV0Reader->GetSphericity() != 0){
          if(fDoJetAnalysis && !fDoSoftAnalysis){
            if(fConvJetReader->GetNJets() > 0){
                fHistoNEvents[iCut]->Fill(eventQuality, fWeightJetJetMC); // Should be 0 here
                if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventQuality); // Should be 0 here
                fHistoNGoodESDTracks[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks(), fWeightJetJetMC);
                fHistoVertexZ[iCut]->Fill(fInputEvent->GetPrimaryVertex()->GetZ(), fWeightJetJetMC);
                fHistoEventSphericity[iCut]->Fill(fV0Reader->GetSphericity(), fWeightJetJetMC);
                fHistoEventSphericityAxis[iCut]->Fill(fV0Reader->GetPhiMainSphericityAxis(), fWeightJetJetMC);
                fHistoEventSphericityvsNtracks[iCut]->Fill(fV0Reader->GetSphericity(), fV0Reader->GetNumberOfPrimaryTracks(), fWeightJetJetMC);
            }
            fHistoEventSphericityvsNJets[iCut]->Fill(fV0Reader->GetSphericity(), fConvJetReader->GetNJets(), fWeightJetJetMC);
            fHistoEventMultiplicityvsNJets[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks(), fConvJetReader->GetNJets(), fWeightJetJetMC);
            if(fIsMC>0){
                fHistoTrueSphericityvsRecSphericity[iCut]->Fill(fV0Reader->GetSphericityTrue(), fV0Reader->GetSphericity(), fWeightJetJetMC);
                fHistoTrueMultiplicityvsRecMultiplicity[iCut]->Fill(fV0Reader->GetNumberOfTruePrimaryTracks(), fV0Reader->GetNumberOfRecTracks(), fWeightJetJetMC);
                ProcessAODSphericityParticles();
            }
            fHistoEventSphericityvsHighpt[iCut]->Fill(fV0Reader->GetSphericity(), fV0Reader->GetHighestPt(), fWeightJetJetMC);
            fHistoEventSphericityvsTotalpt[iCut]->Fill(fV0Reader->GetSphericity(), fV0Reader->GetTotalPt(), fWeightJetJetMC);
            fHistoEventSphericityvsMeanpt[iCut]->Fill(fV0Reader->GetSphericity(), fV0Reader->GetMeanPt(), fWeightJetJetMC);
          }else if(fDoJetAnalysis && fDoSoftAnalysis){
            if(fConvJetReader->GetNJets() < 1){
                fHistoNEvents[iCut]->Fill(eventQuality, fWeightJetJetMC); // Should be 0 here
                if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventQuality); // Should be 0 here
                fHistoNGoodESDTracks[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks(), fWeightJetJetMC);
                fHistoVertexZ[iCut]->Fill(fInputEvent->GetPrimaryVertex()->GetZ(), fWeightJetJetMC);
                fHistoEventSphericity[iCut]->Fill(fV0Reader->GetSphericity(), fWeightJetJetMC);
                fHistoEventSphericityAxis[iCut]->Fill(fV0Reader->GetPhiMainSphericityAxis(), fWeightJetJetMC);
                fHistoEventSphericityvsNtracks[iCut]->Fill(fV0Reader->GetSphericity(), fV0Reader->GetNumberOfPrimaryTracks(), fWeightJetJetMC);
            }
            fHistoEventSphericityvsNJets[iCut]->Fill(fV0Reader->GetSphericity(), fConvJetReader->GetNJets(), fWeightJetJetMC);
            fHistoEventMultiplicityvsNJets[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks(), fConvJetReader->GetNJets(), fWeightJetJetMC);
            if(fIsMC>0){
                fHistoTrueSphericityvsRecSphericity[iCut]->Fill(fV0Reader->GetSphericityTrue(), fV0Reader->GetSphericity(), fWeightJetJetMC);
                fHistoTrueMultiplicityvsRecMultiplicity[iCut]->Fill(fV0Reader->GetNumberOfTruePrimaryTracks(), fV0Reader->GetNumberOfRecTracks(), fWeightJetJetMC);
                ProcessAODSphericityParticles();
            }
            fHistoEventSphericityvsHighpt[iCut]->Fill(fV0Reader->GetSphericity(), fV0Reader->GetHighestPt(), fWeightJetJetMC);
            fHistoEventSphericityvsTotalpt[iCut]->Fill(fV0Reader->GetSphericity(), fV0Reader->GetTotalPt(), fWeightJetJetMC);
            fHistoEventSphericityvsMeanpt[iCut]->Fill(fV0Reader->GetSphericity(), fV0Reader->GetMeanPt(), fWeightJetJetMC);
        }else{
            fHistoNEvents[iCut]->Fill(eventQuality, fWeightJetJetMC); // Should be 0 here
            if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventQuality); // Should be 0 here
            fHistoNGoodESDTracks[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks(), fWeightJetJetMC);
            fHistoVertexZ[iCut]->Fill(fInputEvent->GetPrimaryVertex()->GetZ(), fWeightJetJetMC);
            if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetUseSphericityTrue() && fIsMC>0){
                fHistoEventSphericity[iCut]->Fill(fV0Reader->GetSphericityTrue(), fWeightJetJetMC);
                fHistoEventSphericityvsNtracks[iCut]->Fill(fV0Reader->GetSphericityTrue(), fV0Reader->GetNumberOfTruePrimaryTracks(), fWeightJetJetMC);
            }else{
                fHistoEventSphericity[iCut]->Fill(fV0Reader->GetSphericity(), fWeightJetJetMC);
                fHistoEventSphericityAxis[iCut]->Fill(fV0Reader->GetPhiMainSphericityAxis(), fWeightJetJetMC);
                fHistoEventSphericityvsNtracks[iCut]->Fill(fV0Reader->GetSphericity(), fV0Reader->GetNumberOfRecTracks(), fWeightJetJetMC);
            }
            if(fIsMC>0){
                fHistoTrueSphericityvsRecSphericity[iCut]->Fill(fV0Reader->GetSphericityTrue(), fV0Reader->GetSphericity(), fWeightJetJetMC);
                fHistoTrueMultiplicityvsRecMultiplicity[iCut]->Fill(fV0Reader->GetNumberOfTruePrimaryTracks(), fV0Reader->GetNumberOfRecTracks(), fWeightJetJetMC);
                ProcessAODSphericityParticles();
            }
            fHistoEventSphericityvsHighpt[iCut]->Fill(fV0Reader->GetSphericity(), fV0Reader->GetHighestPt(), fWeightJetJetMC);
            fHistoEventSphericityvsTotalpt[iCut]->Fill(fV0Reader->GetSphericity(), fV0Reader->GetTotalPt(), fWeightJetJetMC);
            fHistoEventSphericityvsMeanpt[iCut]->Fill(fV0Reader->GetSphericity(), fV0Reader->GetMeanPt(), fWeightJetJetMC);
          }
        }
      }else{
        fHistoNEvents[iCut]->Fill(eventQuality, fWeightJetJetMC); // Should be 0 here
        if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventQuality); // Should be 0 here
        fHistoNGoodESDTracks[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks(), fWeightJetJetMC);
        fHistoVertexZ[iCut]->Fill(fInputEvent->GetPrimaryVertex()->GetZ(), fWeightJetJetMC);
      }
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
    if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetDoSecondaryTrackMatching()) ProcessConversionCandidates(); // process conversion candidates for secondary track matching
    if(fDoJetAnalysis)   ProcessJets();

    fHistoNGammaCandidatesBasic[iCut]->Fill(fNCurrentClusterBasic, fWeightJetJetMC);
    fHistoNGammaCandidates[iCut]->Fill(fClusterCandidates->GetEntries(), fWeightJetJetMC);
    if(!fDoLightOutput) fHistoNGoodESDTracksVsNGammaCandidates[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks(),fClusterCandidates->GetEntries(), fWeightJetJetMC);
    if(fDoMesonAnalysis){ // Meson Analysis
      CalculatePi0Candidates(); // Combine Gammas from conversion and from calo
      if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->DoBGCalculation()){
        if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->BackgroundHandlerType() == 0){
          CalculateBackground(); // Combinatorial Background
          UpdateEventByEventData(); // Store Event for mixed Events
        } else if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->BackgroundHandlerType() == 2){
          CalculateBackgroundSwapp(); // Combinatorial Background
        }

      }
      fVectorDoubleCountTruePi0s.clear();
      fVectorDoubleCountTrueEtas.clear();
    }
    if(fIsMC> 0){
      fVectorDoubleCountTrueClusterGammas.clear();
      FillMultipleCountHistoAndClear(fMapMultipleCountTrueClusterGammas,fHistoMultipleCountTrueClusterGamma[iCut]);
    }
    if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetDoSecondaryTrackMatching()) fGammaCandidates->Clear();
    fClusterCandidates->Clear(); // delete cluster candidates
  }
  if (fCloseHighPtClusters) delete fCloseHighPtClusters;
  PostData(1, fOutputContainer);
}

//________________________________________________________________________
void AliAnalysisTaskGammaCalo::ProcessClusters()
{

  Int_t nclus                       = 0;
  TClonesArray * arrClustersProcess = NULL;
  fNCurrentClusterBasic             = 0;
  if(!fCorrTaskSetting.CompareTo("")){
    nclus = fInputEvent->GetNumberOfCaloClusters();
  } else {
    arrClustersProcess = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));
    if(!arrClustersProcess)
      AliFatal(Form("%sClustersBranch was not found in AliAnalysisTaskGammaCalo! Check the correction framework settings!",fCorrTaskSetting.Data()));
    nclus = arrClustersProcess->GetEntries();
  }

  vector<AliAODConversionPhoton*>         vectorCurrentClusters;
  vector<Int_t>                           vectorRejectCluster;
  vector<Double_t>                        vectorPhotonWeight;
  vector<Double_t>                        vectorClusterM02;
  vector<Bool_t>                          vectorIsFromDesiredHeader;
  vector<Int_t>                           vectorCurrentClusters_DDL;

  vector<clusterLabel>  fTrueClusterLabels;

  if(nclus == 0)  return;
  // plotting histograms on cell/tower level, only if extendedMatchAndQA > 1
  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FillHistogramsExtendedQA(fInputEvent,fIsMC);

  // match tracks to clusters
  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchTracksToClusters(fInputEvent,fWeightJetJetMC,kTRUE, fMCEvent);

  // vertex
  Double_t vertex[3] = {0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  Int_t totalActiveCells = 0;
  Double_t totalCellsinClusters = 0;
  Double_t totalUnclusteredE = 0;
  if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetDoFlatEnergySubtraction()){
    Double_t totalClusterEnergy = 0;
    totalActiveCells = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetNactiveEmcalCells();
    totalCellsinClusters = 0;
    for(Long_t i = 0; i < nclus; i++){
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
      if(!clus) continue;
      totalClusterEnergy += clus->E();
      totalCellsinClusters += clus->GetNCells();
      delete clus;
    }
    Double_t totalEDeposit = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetTotalEnergyDeposit(fInputEvent);
    totalUnclusteredE = totalEDeposit - totalClusterEnergy;
  }

  Double_t maxClusterEnergy = -1;
  Int_t maxClusterID        = -1;
  map<Long_t,Int_t> mapIsClusterAccepted;
  map<Long_t,Int_t> mapIsClusterAcceptedWithoutTrackMatch;
  // Loop over EMCal clusters
  Int_t NClusinJets = 0;
  fTrueClusterLabels.clear();
  for(Long_t i = 0; i < nclus; i++){
    Double_t tempClusterWeight        = fWeightJetJetMC;
    Double_t tempPhotonWeight         = fWeightJetJetMC;
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
    if(!clus) continue;

    if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetDoFlatEnergySubtraction()){
      if(totalUnclusteredE!=0 && totalActiveCells!=0 && totalCellsinClusters!=0){
        clus->SetE(clus->E()-(clus->GetNCells()*(totalUnclusteredE)/(totalActiveCells-totalCellsinClusters)));
      }
    }

    // Set the jetjet weight to 1 in case the cluster orignated from the minimum bias header
    if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4){
      if( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(clus->GetLabelAt(0), fMCEvent, fInputEvent) == 2)
        tempClusterWeight = 1;
    }

    if(!((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC, tempClusterWeight,i)){
      if(fProduceTreeEOverP && ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedBeforeTrackMatch() ) mapIsClusterAcceptedWithoutTrackMatch[i] = 1;
      if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetIsAcceptedForBasicCounting())fNCurrentClusterBasic++;
      delete clus;
      continue;
    }

    fNCurrentClusterBasic++;

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
    PhotonCandidate->SetIsCaloPhoton(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType());
    PhotonCandidate->SetCaloClusterRef(i);
    PhotonCandidate->SetLeadingCellID(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FindLargestCellInCluster(clus,fInputEvent));
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
        } // end of label loop
      }
      if(fDoMesonQA>=10){
        if(fInputEvent->IsA()==AliESDEvent::Class()){
          DoClusterMergingStudies(clus,fTrueClusterLabels);
        } else if(fInputEvent->IsA()==AliAODEvent::Class()){
          DoClusterMergingStudiesAOD(clus,fTrueClusterLabels);
        }
      }
    }
    if(fDoJetAnalysis == kTRUE && fIsMC > 0 && fDoClusterQA == kTRUE){
      Bool_t InsideJet = kFALSE;
      fNumberOfClusters[fiCut]->Fill(nclus);
      Float_t clusPos[3]={0,0,0};
      clus->GetPosition(clusPos);
      TVector3 clusterVectorJets(clusPos[0],clusPos[1],clusPos[2]);
      Double_t etaCluster = clusterVectorJets.Eta();
      Double_t phiCluster = clusterVectorJets.Phi();
      if(fConvJetReader->GetTrueNJets()>0){
        fTrueVectorJetEta = fConvJetReader->GetTrueVectorJetEta();
        fTrueVectorJetPhi = fConvJetReader->GetTrueVectorJetPhi();
        Double_t RJetPi0Cand;
        for(Int_t i=0; i<fConvJetReader->GetTrueNJets(); i++){
          Double_t DeltaEta = fTrueVectorJetEta.at(i)-etaCluster;
          Double_t DeltaPhi = abs(fTrueVectorJetPhi.at(i)-phiCluster);
          if(DeltaPhi > M_PI) {
            DeltaPhi = 2*M_PI - DeltaPhi;
          }
          RJetPi0Cand = TMath::Sqrt(pow((DeltaEta),2)+pow((DeltaPhi),2));
          if(fConvJetReader->Get_Jet_Radius() > 0 ){
            if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
              InsideJet = kTRUE;
              NClusinJets++;
              fClusterEtaPhiJets[fiCut]->Fill(phiCluster, etaCluster);
            }
          }
        }
        fTrueVectorJetEta.clear();
        fTrueVectorJetPhi.clear();
      }
      if(clus->GetNLabels() > 1){
        Int_t* mclabelsCluster = clus->GetLabels();
        if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
        for(Int_t k =0; k< (Int_t)clus->GetNLabels(); k++){
          AliAODMCParticle* dummy  = (AliAODMCParticle*) fAODMCTrackArray->At(mclabelsCluster[k]);
          Int_t pdgcode = dummy->GetPdgCode();
          if(k == 0){
            fEnergyRatio[fiCut]->Fill(clus->GetClusterMCEdepFraction(k), dummy->E());
            fEnergyDeposit[fiCut]->Fill(clus->GetClusterMCEdepFraction(k)*clus->E(), dummy->E());
            if(InsideJet){
              fEnergyRatioinJets[fiCut]->Fill(clus->GetClusterMCEdepFraction(k), dummy->E());
              fEnergyDepositinJets[fiCut]->Fill(clus->GetClusterMCEdepFraction(k)*clus->E(), dummy->E());
            }
          }
          if(pdgcode == 22 && k == 0){
            fEnergyRatioGamma1[fiCut]->Fill(clus->GetClusterMCEdepFraction(k), dummy->E());
            fEnergyDepGamma1[fiCut]->Fill(clus->GetClusterMCEdepFraction(k)*clus->E(), dummy->E());
            if(InsideJet){
              fEnergyRatioGamma1inJets[fiCut]->Fill(clus->GetClusterMCEdepFraction(k), dummy->E());
              fEnergyDepGamma1inJets[fiCut]->Fill(clus->GetClusterMCEdepFraction(k)*clus->E(), dummy->E());
            }
          }
          if(pdgcode == 22){
            fEnergyRatioGammaAnywhere[fiCut]->Fill(clus->GetClusterMCEdepFraction(k), dummy->E());
            fEnergyDepGammaAnywhere[fiCut]->Fill(clus->GetClusterMCEdepFraction(k)*clus->E(), dummy->E());
            if(InsideJet){
              fEnergyRatioGammaAnywhereinJets[fiCut]->Fill(clus->GetClusterMCEdepFraction(k), dummy->E());
              fEnergyDepGammaAnywhereinJets[fiCut]->Fill(clus->GetClusterMCEdepFraction(k)*clus->E(), dummy->E());
            }
          }
          if(pdgcode == 22){
            Double_t energyGamma = clus->GetClusterMCEdepFraction(k)*clus->E();
            if(energyGamma < 0.7){
              fEnergyRatioGamma1Helped[fiCut]->Fill(energyGamma);
              if(InsideJet){
                fEnergyRatioGamma1HelpedinJets[fiCut]->Fill(energyGamma);
              }
            }
          }
        }
      }
    }

    fIsFromDesiredHeader          = kTRUE;
    fIsOverlappingWithOtherHeader = kFALSE;
    // test whether largest contribution to cluster orginates in added signals
    if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() > 0){
      // Set the jetjet weight to 1 in case the photon candidate orignated from the minimum bias header
      if ( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4) tempPhotonWeight = 1;
      if ( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 0) fIsFromDesiredHeader = kFALSE;
      if (clus->GetNLabels()>1){
        Int_t* mclabelsCluster = clus->GetLabels();
        if (fLocalDebugFlag > 1)   cout << "testing if other labels in cluster belong to different header, need to test " << (Int_t)clus->GetNLabels()-1 << " additional labels" << endl;
        for (Int_t l = 1; l < (Int_t)clus->GetNLabels(); l++ ){
          if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(mclabelsCluster[l], fMCEvent, fInputEvent, fLocalDebugFlag) == 0) fIsOverlappingWithOtherHeader = kTRUE;
        }
        if (fLocalDebugFlag > 1 && fIsOverlappingWithOtherHeader) cout << "found overlapping header: " << endl;
      }
    }

    fHistoClusAllHeadersGammaPt[fiCut]->Fill(PhotonCandidate->Pt(), tempPhotonWeight);
    if (!fIsFromDesiredHeader) fHistoClusRejectedHeadersGammaPt[fiCut]->Fill(PhotonCandidate->Pt(), tempPhotonWeight);
    if (fIsFromDesiredHeader && fIsOverlappingWithOtherHeader) fHistoClusOverlapHeadersGammaPt[fiCut]->Fill(PhotonCandidate->Pt(), tempPhotonWeight);


    if ( (fIsFromDesiredHeader && !fIsOverlappingWithOtherHeader && !fAllowOverlapHeaders) || (fIsFromDesiredHeader && fAllowOverlapHeaders) ){
      if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType()==2){
        Float_t clusPos[3] ;
        Int_t relID[4];
        Int_t PHOS_Module;
        Int_t PHOS_Module_ix;
        Int_t PHOS_DDL;
        clus->GetPosition(clusPos);
        TVector3 clusterVector(clusPos[0],clusPos[1],clusPos[2]);
        ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetGeomPHOS()->GlobalPos2RelId(clusterVector,relID) ;
        PHOS_Module     = relID[0];
        PHOS_Module_ix  = relID[2];
        PHOS_DDL        = WhichDDL(PHOS_Module, PHOS_Module_ix) ;
        vectorCurrentClusters_DDL.push_back(PHOS_DDL);
      }
      vectorCurrentClusters.push_back(PhotonCandidate);
      vectorPhotonWeight.push_back(tempPhotonWeight);
      vectorClusterM02.push_back(clus->GetM02());
      vectorIsFromDesiredHeader.push_back(fIsFromDesiredHeader);
    } else{
      delete PhotonCandidate;
    }
    delete clus;
    delete tmpvec;
  }
  Bool_t rejected = kFALSE;
  // run conversion recovery in addition
  if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetIsConversionRecovery()){
    rejected = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckForReconstructedConversionPairs(vectorCurrentClusters,vectorRejectCluster);
    if (fLocalDebugFlag > 1 && rejected) cout << "found clusters to be rejected" << endl;
  }

  //Loop to acquire highest energy clusters, which are not bad from analysis bad map (or trigger bad map)
  Int_t highestClusterE_Iter_AnaBM=-1;
  Int_t highestClusterE_Iter_BothBM=-1;
  Double_t highestClusterE_Value_AnaBM=0;
  Double_t highestClusterE_Value_BothBM=0;
  Int_t eventHasL0Flag_ClusE=(fInputHandler->IsEventSelected() & AliVEvent::kPHI7);
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
    //if ( (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsSpecialTrigger()==1)||(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsSpecialTrigger()==0) ){//Only MB
      if (highestClusterE_Iter_AnaBM!=-1){
        fHistoClusGammaE_AnaBM_highestE[fiCut]->Fill(vectorCurrentClusters.at(highestClusterE_Iter_AnaBM)->E(), vectorPhotonWeight.at(highestClusterE_Iter_AnaBM));
      }
    //}
  }

  for (Int_t iter = 0; iter < (Int_t)vectorCurrentClusters.size();iter++){

    if (!((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckVectorForIndexAndAdd(vectorRejectCluster, iter,kFALSE)){
      fIsFromDesiredHeader = vectorIsFromDesiredHeader.at(iter);
      fHistoClusGammaPt[fiCut]->Fill(vectorCurrentClusters.at(iter)->Pt(), vectorPhotonWeight.at(iter));
      fHistoClusGammaE[fiCut]->Fill(vectorCurrentClusters.at(iter)->E(), vectorPhotonWeight.at(iter));
      if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType()==2){
        //--------------------------------------------------
        //if ( (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsSpecialTrigger()==1)||(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsSpecialTrigger()==0) ){ //Only MB
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
      if (!fDoLightOutput && fDoClusterQA > 0) {
          if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType()==2){
            if ((vectorCurrentClusters_DDL.at(iter)<=fDDLRange_HistoClusGamma[1])&&(vectorCurrentClusters_DDL.at(iter)>=fDDLRange_HistoClusGamma[0])){
              if (fIsMC==0){
                fHistoClusGammaE_DDL[vectorCurrentClusters_DDL.at(iter)-fDDLRange_HistoClusGamma[0]][fiCut]->Fill(vectorCurrentClusters.at(iter)->E(), vectorPhotonWeight.at(iter));
                if (fCaloTriggerMimicHelper[fiCut]){
                  //--------------------------------------------------
                  //Only MB
                  if ( (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsSpecialTrigger()==1)||(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsSpecialTrigger()==0) ){
                    //On trigger bad map
                    if (fCaloTriggerMimicHelper[fiCut]->IsClusterIDBadMapTrigger(vectorCurrentClusters.at(iter)->GetCaloClusterRef())){
                      fHistoClusGammaE_DDL_TrBM[vectorCurrentClusters_DDL.at(iter)-fDDLRange_HistoClusGamma[0]][fiCut]->Fill(vectorCurrentClusters.at(iter)->E(), vectorPhotonWeight.at(iter));
                    }
                    //Only triggered events
                    if (fCaloTriggerMimicHelper[fiCut]->GetEventChosenByTrigger()){
                      if (eventHasL0Flag_ClusE){
                        fHistoClusGammaE_DDL_wL0_TrigEv[vectorCurrentClusters_DDL.at(iter)-fDDLRange_HistoClusGamma[0]][fiCut]->Fill(vectorCurrentClusters.at(iter)->E(), vectorPhotonWeight.at(iter));
                      } else {
                        fHistoClusGammaE_DDL_woL0_TrigEv[vectorCurrentClusters_DDL.at(iter)-fDDLRange_HistoClusGamma[0]][fiCut]->Fill(vectorCurrentClusters.at(iter)->E(), vectorPhotonWeight.at(iter));
                      }
                    }
                    //Only triggered events & on trigger bad map
                    if ( (fCaloTriggerMimicHelper[fiCut]->GetEventChosenByTrigger())
                        &&(fCaloTriggerMimicHelper[fiCut]->IsClusterIDBadMapTrigger(vectorCurrentClusters.at(iter)->GetCaloClusterRef()))
                    ){
                      if (eventHasL0Flag_ClusE){
                        fHistoClusGammaE_DDL_wL0_TrigEv_TrBM[vectorCurrentClusters_DDL.at(iter)-fDDLRange_HistoClusGamma[0]][fiCut]->Fill(vectorCurrentClusters.at(iter)->E(), vectorPhotonWeight.at(iter));
                      } else {
                        fHistoClusGammaE_DDL_woL0_TrigEv_TrBM[vectorCurrentClusters_DDL.at(iter)-fDDLRange_HistoClusGamma[0]][fiCut]->Fill(vectorCurrentClusters.at(iter)->E(), vectorPhotonWeight.at(iter));
                      }
                    }
                  } //OnlyMB
                } //fCaloTriggerMimicHelper[fiCut]
              } //fIsMC==0
            }
          }
          fHistoClusGammaPtM02[fiCut]->Fill(vectorCurrentClusters.at(iter)->Pt(), vectorClusterM02.at(iter), vectorPhotonWeight.at(iter));
      }
      if(fIsMC> 0){
        if(fInputEvent->IsA()==AliESDEvent::Class()){
          ProcessTrueClusterCandidates(vectorCurrentClusters.at(iter),vectorClusterM02.at(iter));
        } else {
          ProcessTrueClusterCandidatesAOD(vectorCurrentClusters.at(iter),vectorClusterM02.at(iter));
        }
      }
      fClusterCandidates->Add(vectorCurrentClusters.at(iter));
    } else {
      if (fLocalDebugFlag > 1) cout << "removed cluster: " << iter << endl;
    }
  }
  if (fLocalDebugFlag > 1 && rejected) cout << "reduced from: " << vectorCurrentClusters.size() << " to " << fClusterCandidates->GetEntries() << endl;
  vectorRejectCluster.clear();
  vectorPhotonWeight.clear();
  vectorClusterM02.clear();
  vectorIsFromDesiredHeader.clear();

  if(fDoJetAnalysis == kTRUE && fIsMC > 0 && fDoClusterQA == kTRUE){
    if(NClusinJets > 0)  fNumberOfClustersinJets[fiCut]->Fill(NClusinJets);
  }

  if(fProduceCellIDPlots){
    for(Long_t i = 0; i < nclus; i++){
      if( mapIsClusterAccepted[i] != 1 ) continue;

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
        if(fInputEvent->IsA()==AliESDEvent::Class()){
          if(arrClustersProcess)
            secondClus = new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersProcess->At(j));
          else
            secondClus = new AliESDCaloCluster(*(AliESDCaloCluster*)fInputEvent->GetCaloCluster(j));
        } else if(fInputEvent->IsA()==AliAODEvent::Class()){
          if(arrClustersProcess)
            secondClus = new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersProcess->At(j));
          else
            secondClus = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(j));
        }

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

  // Do more merged studies
  if(fDoMesonQA>=10 && fTrueClusterLabels.size()> 1){ // we need two labels for this to be possible

      clusterLabel thisLabel = fTrueClusterLabels.at(0);
      // printf("MesonID=%i clusID=%i daughterID=%i daughterPDG=%i EClus=%f EFrac=%f ETrue=%f \n",
      // thisLabel.mesonID,thisLabel.clusID,thisLabel.daughterID,thisLabel.daughterPDG,thisLabel.EClus,thisLabel.EFrac,thisLabel.ETrue);
      Bool_t isSep = kFALSE;
      Bool_t isMerged    = kFALSE;

      Double_t pi0Pt = thisLabel.PtMeson;
      Double_t pi0Eta = thisLabel.EtaMeson;

      // Double_t angle = thisLabel.OpeningAngle;

      TLorentzVector thisClus = thisLabel.clusVec;
      Int_t pdg = 0;
      // search other labels to find separated
      for (UInt_t b = 1; b < fTrueClusterLabels.size(); b++)
      {
        clusterLabel otherLabel = fTrueClusterLabels.at(b);
        TLorentzVector otherClus = otherLabel.clusVec;
        TLorentzVector sumVec = otherClus + thisClus;
        Double_t m = sumVec.M();
        // printf("MesonID=%i clusID=%i daughterID=%i daughterPDG=%i EClus=%f EFrac=%f ETrue=%f \n",
        // otherLabel.mesonID,otherLabel.clusID,otherLabel.daughterID,otherLabel.daughterPDG,otherLabel.EClus,otherLabel.EFrac,otherLabel.ETrue);
        if(otherLabel.mesonID != thisLabel.mesonID) continue;

        if(fInputEvent->IsA()==AliESDEvent::Class()){
          pdg = ((AliMCParticle *)fMCEvent->GetTrack(thisLabel.mesonID))->PdgCode();
        } else{
          pdg = ((AliAODMCParticle *)fAODMCTrackArray->At(thisLabel.mesonID))->PdgCode();
        }

        if(thisLabel.clusID != otherLabel.clusID){
            if(thisLabel.daughterID != otherLabel.daughterID){
               isSep = kTRUE;
              if(pdg==221){
                  if((m > 0.650) || (m < 0.450) )
                    isSep = kFALSE;
                  // if( angle < 0.0143){
                  //   // angle too small to be not merged
                  //   isSep = kFALSE;
                  // }
              } else if (pdg==111){
                  if((m > (0.134 + (0.134*0.25))) || (m < (0.134 - (0.134*0.25))) )
                    isSep = kFALSE;
              }
            }
        }
      }
      // loop again to find merged
      if(!isSep){
        for (UInt_t b = 1; b < fTrueClusterLabels.size(); b++)
        {
          clusterLabel otherLabel = fTrueClusterLabels.at(b);
          if(otherLabel.mesonID != thisLabel.mesonID) continue;
          if(thisLabel.clusID == otherLabel.clusID){
              if(thisLabel.daughterID != otherLabel.daughterID){
                isMerged = kTRUE;
              }
          }
        }
      }
      // Begin fillig of histos
      if((pdg ==111) && isMerged) fHistoMCPi0GenFoundInOneCluster[fiCut]->Fill(pi0Pt,pi0Eta,fWeightJetJetMC);
      if((pdg ==111) && isSep) fHistoMCPi0GenFoundInTwoCluster[fiCut]->Fill(pi0Pt,pi0Eta,fWeightJetJetMC);
      if((pdg ==221) && isMerged) fHistoMCEtaGenFoundInOneCluster[fiCut]->Fill(pi0Pt,pi0Eta,fWeightJetJetMC);
      if((pdg ==221) && isSep) fHistoMCEtaGenFoundInTwoCluster[fiCut]->Fill(pi0Pt,pi0Eta,fWeightJetJetMC);
  }

  return;
}

//________________________________________________________________________
void AliAnalysisTaskGammaCalo::ProcessConversionCandidates(){

  // Loop over Photon Candidates allocated by ReaderV1
  for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){
    AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(i);
    if(!PhotonCandidate) continue;
    if(!fConversionCuts->PhotonIsSelected(PhotonCandidate,fInputEvent)) continue;
    if(!fConversionCuts->InPlaneOutOfPlaneCut(PhotonCandidate->GetPhotonPhi(),fEventPlaneAngle)) continue;
    if(!fConversionCuts->UseElecSharingCut() && !fConversionCuts->UseToCloseV0sCut()){
      fGammaCandidates->Add(PhotonCandidate);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskGammaCalo::ProcessJets()
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
      fJetSector = 0;
      Double_t tempMaxJetPt = 0.;
      Int_t MaxPtJetID = 0;
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
        if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoSectorJetMixing()){ //Determine jet with max pt in event
          if(fVectorJetPt.at(i) > tempMaxJetPt){
            tempMaxJetPt = fVectorJetPt.at(i);
            MaxPtJetID = i;
          }
        }
      }

      if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoSectorJetMixing()){
        Double_t PtMaxPhi = fVectorJetPhi.at(MaxPtJetID);
        fJetSector = 0;
        if((PtMaxPhi > (TMath::Pi()/4.)) && (PtMaxPhi < (3*TMath::Pi()/4.))) fJetSector = 1;
        if((PtMaxPhi > (3*TMath::Pi()/4.)) && (PtMaxPhi < (5*TMath::Pi()/4.))) fJetSector = 2;
        if((PtMaxPhi > (5*TMath::Pi()/4.)) && (PtMaxPhi < (7*TMath::Pi()/4.))) fJetSector = 3;
        if((PtMaxPhi > (7*TMath::Pi()/4.)) && (PtMaxPhi < (TMath::Pi()/4.))) fJetSector = 4;
        if(fVectorJetEta.at(MaxPtJetID)>0) fJetSector += 4;
      }else if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoJetMixing()){
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
void AliAnalysisTaskGammaCalo::ProcessTrueClusterCandidates(AliAODConversionPhoton *TruePhotonCandidate, Double_t clusterM02)
{
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  Double_t tempPhotonWeight       = fWeightJetJetMC;
  TParticle *Photon = NULL;
  if (TruePhotonCandidate->GetIsCaloPhoton() == 0) AliFatal("CaloPhotonFlag has not been set task will abort");
  if (!fDoLightOutput && fDoClusterQA > 0) fHistoTrueNLabelsInClus[fiCut]->Fill(TruePhotonCandidate->GetNCaloPhotonMCLabels(), tempPhotonWeight);

  if (TruePhotonCandidate->GetCaloPhotonMCLabel(0) < 0) return;
  if (TruePhotonCandidate->GetNCaloPhotonMCLabels() > 0) Photon = fMCEvent->Particle(TruePhotonCandidate->GetCaloPhotonMCLabel(0));
    else return;

  if(Photon == NULL){
  //    cout << "no photon" << endl;
    return;
  }
  // Set the jetjet weight to 1 in case the photon orignated from the minimum bias header
  if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4){
    if ( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(TruePhotonCandidate->GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2) tempPhotonWeight = 1;
  }

  Int_t pdgCodeParticle = Photon->GetPdgCode();
  TruePhotonCandidate->SetCaloPhotonMCFlags(fMCEvent, fEnableSortForClusMC);

  // True Photon

  if (TruePhotonCandidate->IsLargestComponentPhoton() || (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion())) {
    fHistoTrueClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
    if (!fDoLightOutput && fDoClusterQA > 0) fHistoTrueClusGammaPtM02[fiCut]->Fill(TruePhotonCandidate->Pt(), clusterM02, tempPhotonWeight);
  } else if (fDoClusterQA > 0) fHistoTrueClusEMNonLeadingPt[fiCut]->Fill(TruePhotonCandidate->Pt());
  if (fDoClusterQA > 0){
    if (TruePhotonCandidate->IsLargestComponentPhoton()){
      fHistoTrueClusUnConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
      fHistoTrueClusUnConvGammaMCPt[fiCut]->Fill(Photon->Pt(), tempPhotonWeight);
      if (!fDoLightOutput) fHistoTrueClusUnConvGammaPtM02[fiCut]->Fill(TruePhotonCandidate->Pt(), clusterM02, tempPhotonWeight);
    }
    if (TruePhotonCandidate->IsLargestComponentElectron())
      fHistoTrueClusElectronPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
    if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
      fHistoTrueClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
      if(Photon->GetMother(0) > -1){
        TParticle* motherPart = (TParticle*)fMCEvent->Particle(Photon->GetMother(0));
        fHistoTrueClusConvGammaMCPt[fiCut]->Fill(motherPart->Pt(), tempPhotonWeight);
      }
    }
    if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion() && TruePhotonCandidate->IsConversionFullyContained())
      fHistoTrueClusConvGammaFullyPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
    if (TruePhotonCandidate->IsMerged() || TruePhotonCandidate->IsMergedPartConv() || TruePhotonCandidate->IsDalitzMerged())
      fHistoTrueClusMergedGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
    if (TruePhotonCandidate->IsMergedPartConv())
      fHistoTrueClusMergedPartConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
    if (TruePhotonCandidate->IsDalitz())
      fHistoTrueClusDalitzPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
    if (TruePhotonCandidate->IsDalitzMerged())
      fHistoTrueClusDalitzMergedPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
    if (TruePhotonCandidate->IsPhotonWithElecMother())
      fHistoTrueClusPhotonFromElecMotherPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
    if (TruePhotonCandidate->IsShower())
      fHistoTrueClusShowerPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
    if (TruePhotonCandidate->IsSubLeadingEM())
      fHistoTrueClusSubLeadingPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
    fHistoTrueClusNParticles[fiCut]->Fill(TruePhotonCandidate->GetNCaloPhotonMCLabels(), tempPhotonWeight);
    if (!TruePhotonCandidate->IsLargestComponentPhoton()) {
      FillPhotonBackgroundHist(TruePhotonCandidate,pdgCodeParticle);
      if (!fDoLightOutput) FillPhotonBackgroundM02Hist(TruePhotonCandidate,clusterM02,pdgCodeParticle);
    }
    if (!(TruePhotonCandidate->IsLargestComponentPhoton() || (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion())) ) {
      FillPhotonPlusConversionBackgroundHist(TruePhotonCandidate,pdgCodeParticle);
      if (!fDoLightOutput) FillPhotonPlusConversionBackgroundM02Hist(TruePhotonCandidate,clusterM02,pdgCodeParticle);
    }
    Int_t motherLab = Photon->GetMother(0);
    if (motherLab > -1){
      if (TruePhotonCandidate->IsLargestComponentPhoton()){
        if (CheckVectorForDoubleCount(fVectorDoubleCountTrueClusterGammas,motherLab)){
          fHistoDoubleCountTrueClusterGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),(Double_t)0,tempPhotonWeight);
          FillMultipleCountMap(fMapMultipleCountTrueClusterGammas,motherLab);
        }
      }
      Int_t grandMotherLab = fMCEvent->Particle(motherLab)->GetMother(0);
      if (grandMotherLab > -1){
        if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
          if (CheckVectorForDoubleCount(fVectorDoubleCountTrueClusterGammas,grandMotherLab)){
            fHistoDoubleCountTrueClusterGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),(Double_t)1,tempPhotonWeight);
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
        fHistoTruePrimaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
        fHistoTruePrimaryClusGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt(), tempPhotonWeight); // Allways Filled
      }
      if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
        fHistoTruePrimaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
        fHistoTruePrimaryClusConvGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt(), tempPhotonWeight); // Allways Filled
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
          fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),0.,tempPhotonWeight);
          fHistoTrueSecondaryClusGammaMCPt[fiCut]->Fill(Photon->Pt(),0.,tempPhotonWeight);
          fHistoTrueSecondaryClusGammaFromXFromK0sMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),tempPhotonWeight);
        } else if (secondaryClass == 5) {
          fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),1.,tempPhotonWeight);
          fHistoTrueSecondaryClusGammaMCPt[fiCut]->Fill(Photon->Pt(),1.,tempPhotonWeight);
          fHistoTrueSecondaryClusGammaFromXFromK0lMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),tempPhotonWeight);
        } else if (secondaryClass == 3) {
          fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),2.,tempPhotonWeight);
          fHistoTrueSecondaryClusGammaMCPt[fiCut]->Fill(Photon->Pt(),2.,tempPhotonWeight);
          fHistoTrueSecondaryClusGammaFromXFromLambdaMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),tempPhotonWeight);
        } else if (secondaryClass == 4) {
          fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3.,tempPhotonWeight);
          fHistoTrueSecondaryClusGammaMCPt[fiCut]->Fill(Photon->Pt(),3.,tempPhotonWeight);
        } else {
          fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),4.,tempPhotonWeight);
          fHistoTrueSecondaryClusGammaMCPt[fiCut]->Fill(Photon->Pt(),4.,tempPhotonWeight);
        }
      }
      if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
        if (secondaryClass == 2) {
          fHistoTrueSecondaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),0.,tempPhotonWeight);
          fHistoTrueSecondaryClusConvGammaMCPt[fiCut]->Fill(Photon->Pt(),0.,tempPhotonWeight);
          fHistoTrueSecondaryClusConvGammaFromXFromK0sMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),tempPhotonWeight);
        } else if (secondaryClass == 5) {
          fHistoTrueSecondaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),1.,tempPhotonWeight);
          fHistoTrueSecondaryClusConvGammaMCPt[fiCut]->Fill(Photon->Pt(),1.,tempPhotonWeight);
          fHistoTrueSecondaryClusConvGammaFromXFromK0lMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),tempPhotonWeight);
        } else if (secondaryClass == 3) {
          fHistoTrueSecondaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),2.,tempPhotonWeight);
          fHistoTrueSecondaryClusConvGammaMCPt[fiCut]->Fill(Photon->Pt(),2.,tempPhotonWeight);
          fHistoTrueSecondaryClusConvGammaFromXFromLambdaMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),tempPhotonWeight);
        } else if (secondaryClass == 4) {
          fHistoTrueSecondaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3.,tempPhotonWeight);
          fHistoTrueSecondaryClusConvGammaMCPt[fiCut]->Fill(Photon->Pt(),3.,tempPhotonWeight);
        } else {
          fHistoTrueSecondaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),4.,tempPhotonWeight);
          fHistoTrueSecondaryClusConvGammaMCPt[fiCut]->Fill(Photon->Pt(),4.,tempPhotonWeight);
        }
      }
    }
  }
  return;
}


//________________________________________________________________________
void AliAnalysisTaskGammaCalo::ProcessTrueClusterCandidatesAOD(AliAODConversionPhoton *TruePhotonCandidate, Double_t clusterM02)
{
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  Double_t tempPhotonWeight       = fWeightJetJetMC;
  AliAODMCParticle *Photon = NULL;
  if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (!fDoLightOutput && fDoClusterQA > 0) fHistoTrueNLabelsInClus[fiCut]->Fill(TruePhotonCandidate->GetNCaloPhotonMCLabels(), tempPhotonWeight);
  if (fAODMCTrackArray){
    if (TruePhotonCandidate->GetIsCaloPhoton() == 0) AliFatal("CaloPhotonFlag has not been set task will abort");
    if (TruePhotonCandidate->GetNCaloPhotonMCLabels()>0) Photon = (AliAODMCParticle*) fAODMCTrackArray->At(TruePhotonCandidate->GetCaloPhotonMCLabel(0));
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
  TruePhotonCandidate->SetCaloPhotonMCFlagsAOD(fAODMCTrackArray, fEnableSortForClusMC);

  // Set the jetjet weight to 1 in case the cluster orignated from the minimum bias header
  if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4){
    if ( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(TruePhotonCandidate->GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2) tempPhotonWeight = 1;
  }

  // True Photon
  if (TruePhotonCandidate->IsLargestComponentPhoton() || (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion())) {
    fHistoTrueClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
    if (!fDoLightOutput && fDoClusterQA > 0) fHistoTrueClusGammaPtM02[fiCut]->Fill(TruePhotonCandidate->Pt(), clusterM02, tempPhotonWeight);
  } else if (fDoClusterQA > 0) fHistoTrueClusEMNonLeadingPt[fiCut]->Fill(TruePhotonCandidate->Pt());
  if (fDoClusterQA > 0){
    if (TruePhotonCandidate->IsLargestComponentPhoton()) {
      fHistoTrueClusUnConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
      fHistoTrueClusUnConvGammaMCPt[fiCut]->Fill(Photon->Pt(), tempPhotonWeight);
      if (!fDoLightOutput) fHistoTrueClusUnConvGammaPtM02[fiCut]->Fill(TruePhotonCandidate->Pt(), clusterM02, tempPhotonWeight);
    }
    if (TruePhotonCandidate->IsLargestComponentElectron())
      fHistoTrueClusElectronPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
    if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()) {
      fHistoTrueClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
      if(Photon->GetMother()>-1){
        AliAODMCParticle *motherPart = (AliAODMCParticle*) fAODMCTrackArray->At(Photon->GetMother());
        fHistoTrueClusConvGammaMCPt[fiCut]->Fill(motherPart->Pt(), tempPhotonWeight);
      }
    }
    if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion() && TruePhotonCandidate->IsConversionFullyContained())
      fHistoTrueClusConvGammaFullyPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
    if (TruePhotonCandidate->IsMerged() || TruePhotonCandidate->IsMergedPartConv() || TruePhotonCandidate->IsDalitzMerged())
      fHistoTrueClusMergedGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
    if (TruePhotonCandidate->IsMergedPartConv())
      fHistoTrueClusMergedPartConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
    if (TruePhotonCandidate->IsDalitz())
      fHistoTrueClusDalitzPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
    if (TruePhotonCandidate->IsDalitzMerged())
      fHistoTrueClusDalitzMergedPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
    if (TruePhotonCandidate->IsPhotonWithElecMother())
      fHistoTrueClusPhotonFromElecMotherPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
    if (TruePhotonCandidate->IsShower())
      fHistoTrueClusShowerPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
    if (TruePhotonCandidate->IsSubLeadingEM())
      fHistoTrueClusSubLeadingPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
    fHistoTrueClusNParticles[fiCut]->Fill(TruePhotonCandidate->GetNCaloPhotonMCLabels(), tempPhotonWeight);

    if (!TruePhotonCandidate->IsLargestComponentPhoton()) {
      FillPhotonBackgroundHist(TruePhotonCandidate,pdgCodeParticle);
      if (!fDoLightOutput) FillPhotonBackgroundM02Hist(TruePhotonCandidate,clusterM02,pdgCodeParticle);
    }
    if (!(TruePhotonCandidate->IsLargestComponentPhoton() || (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion())) ) {
      FillPhotonPlusConversionBackgroundHist(TruePhotonCandidate,pdgCodeParticle);
      if (!fDoLightOutput) FillPhotonPlusConversionBackgroundM02Hist(TruePhotonCandidate,clusterM02,pdgCodeParticle);
    }
    Int_t motherLab = Photon->GetMother();
    if (motherLab > -1){
      if (TruePhotonCandidate->IsLargestComponentPhoton()){
        if (CheckVectorForDoubleCount(fVectorDoubleCountTrueClusterGammas,motherLab)){
          fHistoDoubleCountTrueClusterGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),(Double_t)0,tempPhotonWeight);
          FillMultipleCountMap(fMapMultipleCountTrueClusterGammas,motherLab);
        }
      }
      Int_t grandMotherLab = ((AliAODMCParticle*) fAODMCTrackArray->At(motherLab))->GetMother();
      if (grandMotherLab > -1){
        if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
          if (CheckVectorForDoubleCount(fVectorDoubleCountTrueClusterGammas,grandMotherLab)){
            fHistoDoubleCountTrueClusterGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),(Double_t)1,tempPhotonWeight);
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
        AliAODMCParticle *Mother  = (AliAODMCParticle*) fAODMCTrackArray->At(Photon->GetMother());
        isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, Mother, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
      }
    }
    if(isPrimary){
      if (TruePhotonCandidate->IsLargestComponentPhoton()){
        fHistoTruePrimaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
        fHistoTruePrimaryClusGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt(), tempPhotonWeight); // Allways Filled
      }
      if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
        fHistoTruePrimaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
        fHistoTruePrimaryClusConvGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt(), tempPhotonWeight); // Allways Filled
      }

    } else {
      // filling secondary histograms
      Int_t secondaryClass    = -1;
      if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion())
        secondaryClass        = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->SecondaryClassificationPhotonAOD( Photon, fAODMCTrackArray, kTRUE);
      else
        secondaryClass        = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->SecondaryClassificationPhotonAOD( Photon, fAODMCTrackArray, kFALSE);

      // all secondaries
      if (TruePhotonCandidate->IsLargestComponentPhoton()){
        if (secondaryClass == 2) {
          fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),0.,tempPhotonWeight);
          fHistoTrueSecondaryClusGammaMCPt[fiCut]->Fill(Photon->Pt(),0.,tempPhotonWeight);
          fHistoTrueSecondaryClusGammaFromXFromK0sMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),tempPhotonWeight);
        } else if (secondaryClass == 5) {
          fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),1.,tempPhotonWeight);
          fHistoTrueSecondaryClusGammaMCPt[fiCut]->Fill(Photon->Pt(),1.,tempPhotonWeight);
          fHistoTrueSecondaryClusGammaFromXFromK0lMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),tempPhotonWeight);
        } else if (secondaryClass == 3) {
          fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),2.,tempPhotonWeight);
          fHistoTrueSecondaryClusGammaMCPt[fiCut]->Fill(Photon->Pt(),2.,tempPhotonWeight);
          fHistoTrueSecondaryClusGammaFromXFromLambdaMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),tempPhotonWeight);
        } else if (secondaryClass == 4) {
          fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3.,tempPhotonWeight);
          fHistoTrueSecondaryClusGammaMCPt[fiCut]->Fill(Photon->Pt(),3.,tempPhotonWeight);
        } else {
          fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),4.,tempPhotonWeight);
          fHistoTrueSecondaryClusGammaMCPt[fiCut]->Fill(Photon->Pt(),4.,tempPhotonWeight);
        }
      }
      if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
        if (secondaryClass == 2) {
          fHistoTrueSecondaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),0.,tempPhotonWeight);
          fHistoTrueSecondaryClusConvGammaMCPt[fiCut]->Fill(Photon->Pt(),0.,tempPhotonWeight);
          fHistoTrueSecondaryClusConvGammaFromXFromK0sMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),tempPhotonWeight);
        } else if (secondaryClass == 5) {
          fHistoTrueSecondaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),1.,tempPhotonWeight);
          fHistoTrueSecondaryClusConvGammaMCPt[fiCut]->Fill(Photon->Pt(),1.,tempPhotonWeight);
          fHistoTrueSecondaryClusConvGammaFromXFromK0lMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),tempPhotonWeight);
        } else if (secondaryClass == 3) {
          fHistoTrueSecondaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),2.,tempPhotonWeight);
          fHistoTrueSecondaryClusConvGammaMCPt[fiCut]->Fill(Photon->Pt(),2.,tempPhotonWeight);
          fHistoTrueSecondaryClusConvGammaFromXFromLambdaMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),tempPhotonWeight);
        } else if (secondaryClass == 4) {
          fHistoTrueSecondaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3.,tempPhotonWeight);
          fHistoTrueSecondaryClusConvGammaMCPt[fiCut]->Fill(Photon->Pt(),3.,tempPhotonWeight);
        } else {
          fHistoTrueSecondaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),4.,tempPhotonWeight);
          fHistoTrueSecondaryClusConvGammaMCPt[fiCut]->Fill(Photon->Pt(),4.,tempPhotonWeight);
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


  if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (fAODMCTrackArray == NULL) return;

  // Loop over all primary MC particle
  for(Long_t i = 0; i < fAODMCTrackArray->GetEntriesFast(); i++) {
    Double_t tempParticleWeight       = fWeightJetJetMC;

    AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(i));
    if (!particle) continue;

    Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, particle, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
    if (isPrimary){

      Int_t isMCFromMBHeader = -1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        isMCFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
        if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
        // Set the jetjet weight to 1 in case the particle orignated from the minimum bias header
        if(isMCFromMBHeader == 2 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4) tempParticleWeight = 1;
      }

      if(!fDoLightOutput){
        if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(particle,fAODMCTrackArray)){
          fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt(), tempParticleWeight); // All MC Gamma
          if(particle->GetMother() >-1){ // Meson Decay Gamma
            switch((static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetMother())))->GetPdgCode()){
            case 111: // Pi0
              fHistoMCDecayGammaPi0Pt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
              break;
            case 113: // Rho0
              fHistoMCDecayGammaRhoPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
              break;
            case 221: // Eta
              fHistoMCDecayGammaEtaPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
              break;
            case 223: // Omega
              fHistoMCDecayGammaOmegaPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
              break;
            case 331: // Eta'
              fHistoMCDecayGammaEtapPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
              break;
            case 333: // Phi
              fHistoMCDecayGammaPhiPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
              break;
            case 3212: // Sigma
              fHistoMCDecayGammaSigmaPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
              break;
            }
          }
        }
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

        if ((mesonY > ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetRapidityCutValueMin()) && (mesonY < ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetRapidityCutValueMax())){
          if ( particle->GetPdgCode() == 211 ){  // positve pions
            fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 0., tempParticleWeight);
          } else if ( particle->GetPdgCode() == -211 ){  // negative pions
            fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 1., tempParticleWeight);
          } else if ( particle->GetPdgCode() == 321 ){  // positve kaons
            fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 2., tempParticleWeight);
          } else if ( particle->GetPdgCode() == -321 ){  // negative kaons
            fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 3., tempParticleWeight);
          } else if ( TMath::Abs(particle->GetPdgCode()) == 310 ){  // K0s
            fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 4., tempParticleWeight);
          } else if ( TMath::Abs(particle->GetPdgCode()) == 130 ){  // K0l
            fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 5., tempParticleWeight);
          } else if ( TMath::Abs(particle->GetPdgCode()) == 3122 ){  // Lambda/ AntiLambda
            fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 6., tempParticleWeight);
          }
        }

        if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
          ->MesonIsSelectedAODMC(particle,fAODMCTrackArray,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
          AliAODMCParticle* daughter0 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetDaughterLabel(0)));
          AliAODMCParticle* daughter1 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetDaughterLabel(1)));
          Float_t weighted= 1;
          if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent)){
            if (particle->Pt()>0.005){
              weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(i, 0x0, fInputEvent);
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
          Double_t alpha = -1;
          if (particle->GetPdgCode() == 111 || particle->GetPdgCode() == 221){
            alpha = TMath::Abs((daughter0->E() - daughter1->E()))/(daughter0->E() + daughter1->E());
          }

          if(particle->GetPdgCode() == 111){
            if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoMCPi0Pt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // All MC Pi0
            if(fDoJetAnalysis){
              if(fConvJetReader->GetTrueNJets()>0){
                if(!fDoLightOutput) fHistoMCPi0JetEventGenerated[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc in jet event
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
                      if(!fDoLightOutput) fHistoMCPi0inJetGenerated[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc in a jet
                      else fHistoMCPi0Pt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight);
                    }
                  }
                }
                fTrueVectorJetEta.clear();
                fTrueVectorJetPhi.clear();
              }
            }
            fHistoMCPi0WOWeightPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
            if (fIsMC > 1) fHistoMCPi0WOEvtWeightPt[fiCut]->Fill(particle->Pt());
            if (fDoMesonQA > 0 && fDoMesonQA < 3){
              fHistoMCPi0PtY[fiCut]->Fill(particle->Pt(),mesonY,weighted* tempParticleWeight);
              fHistoMCPi0PtAlpha[fiCut]->Fill(particle->Pt(),alpha, tempParticleWeight);
              if (fIsMC == 2) fHistoMCPi0PtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),tempParticleWeight);
            }
          } else if(particle->GetPdgCode() == 221 && !fDoPi0Only){
            if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoMCEtaPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // All MC Eta
            if(fDoJetAnalysis){
              if(fConvJetReader->GetTrueNJets()>0){
                if(!fDoLightOutput) fHistoMCEtaJetEventGenerated[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc in jet event
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
                      if(!fDoLightOutput) fHistoMCEtainJetGenerated[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc in a jet
                      else fHistoMCEtaPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight);
                    }
                  }
                }
                fTrueVectorJetEta.clear();
                fTrueVectorJetPhi.clear();
              }
            }
            fHistoMCEtaWOWeightPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
            if (fIsMC > 1) fHistoMCEtaWOEvtWeightPt[fiCut]->Fill(particle->Pt());
            if (fDoMesonQA > 0 && fDoMesonQA < 3){
              fHistoMCEtaPtY[fiCut]->Fill(particle->Pt(),mesonY,weighted* tempParticleWeight);
              fHistoMCEtaPtAlpha[fiCut]->Fill(particle->Pt(),alpha, tempParticleWeight);
              if (fIsMC == 2) fHistoMCEtaPtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),tempParticleWeight);
            }
          }

          // Check the acceptance for both gammas
          if( ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(daughter0,fAODMCTrackArray) &&
              ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(daughter1,fAODMCTrackArray) ){

            if(particle->GetPdgCode() == 111){
              if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoMCPi0InAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc
              if(fIsMC > 1) fHistoMCPi0WOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Pi0 with gamma in acc
              if(fDoJetAnalysis){
                if(fConvJetReader->GetTrueNJets()>0){
                  if(!fDoLightOutput) fHistoMCPi0JetInAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc in jet event
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
                        if(!fDoLightOutput) fHistoMCPi0inJetInAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc in a jet
                        else fHistoMCPi0InAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight);
                      }
                    }
                  }
                  fTrueVectorJetEta.clear();
                  fTrueVectorJetPhi.clear();
                }
              }
            } else if(particle->GetPdgCode() == 221 && !fDoPi0Only){
              if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoMCEtaInAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Eta with gamma in acc
              if(fIsMC > 1) fHistoMCEtaWOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Eta with gamma in acc
              if(fDoJetAnalysis){
                if(fConvJetReader->GetTrueNJets()>0){
                  if(!fDoLightOutput) fHistoMCEtaJetInAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc in jet event
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
                        if(!fDoLightOutput) fHistoMCEtainJetInAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc in a jet
                        else fHistoMCEtaInAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight);
                      }
                    }
                  }
                  fTrueVectorJetEta.clear();
                  fTrueVectorJetPhi.clear();
                }
              }
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
        if(isMCFromMBHeader == 2 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4) tempParticleWeight = 1;
      }

      if(!fDoLightOutput) {
        if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(particle,fAODMCTrackArray)){
          if(particle->GetMother() >-1){
            AliAODMCParticle *tmpMother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetMother()));
            if(tmpMother->GetMother() >-1){
              AliAODMCParticle *tmpGrandMother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(tmpMother->GetMother()));
              if(tmpGrandMother->GetPdgCode() == 310) {
                fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),0.,tempParticleWeight);
              } else if (tmpGrandMother->GetPdgCode() == 130) {
                fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),1.,tempParticleWeight);
              } else if (tmpGrandMother->GetPdgCode() == 3122) {
                fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),2.,tempParticleWeight);
              } else if (tmpGrandMother->GetPdgCode() == 221) {
                fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),3.,tempParticleWeight);
              } else {
                if( !(TMath::Abs(tmpMother->GetPdgCode()) == 11 && tmpGrandMother->GetPdgCode() == 22) )
                  fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),4.,tempParticleWeight);
              }
            } else {
              fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),4.,tempParticleWeight);
            }
          } else {
            fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),4.,tempParticleWeight);
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
            fHistoMCSecPi0PtvsSource[fiCut]->Fill(particle->Pt(),source,tempParticleWeight); // All MC Pi0
            fHistoMCSecPi0Source[fiCut]->Fill(pdgCode);
          }else if(particle->GetPdgCode() == 221 && !fDoPi0Only){

            fHistoMCSecEtaPt[fiCut]->Fill(particle->Pt(),tempParticleWeight); // All MC Pi0
            fHistoMCSecEtaSource[fiCut]->Fill(pdgCode);
          }

          // check if conversion where within acceptance
          if( ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(daughter0,fAODMCTrackArray) &&
              ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(daughter1,fAODMCTrackArray)){
            if(particle->GetPdgCode() == 111){
              Int_t source = GetSourceClassification(111,pdgCode);
              fHistoMCSecPi0InAccPtvsSource[fiCut]->Fill(particle->Pt(),source,tempParticleWeight); // All MC Pi0
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
      Double_t tempParticleWeight       = fWeightJetJetMC;

      TParticle* particle = (TParticle *)fMCEvent->Particle(i);
      if (!particle) continue;

      // This section produces histograms to study cluster merging
      if(fDoMesonQA>=10){
        if(particle->GetPdgCode() == 111){ // only look at pi0
            Int_t nDaughters = particle->GetNDaughters();
            TParticle* particle = (TParticle *)fMCEvent->Particle(i);
            TParticle* daughter[2] = {NULL,NULL};
            if(nDaughters==2){
              Int_t nPhotons = 0;
              for (Int_t d = 0; d < nDaughters; d++)
              {
                Int_t label = particle->GetDaughter(d);
                daughter[d] = (TParticle *)fMCEvent->Particle(label);
                if(daughter[d]->GetPdgCode()==22){
                  nPhotons++;
                  continue;
                } else{
                  daughter[d] = NULL;
                }
              }
              Int_t nclus;

              TClonesArray * arrClustersProcessTmp = NULL;
              if(!fCorrTaskSetting.CompareTo("")){
                nclus = fInputEvent->GetNumberOfCaloClusters();
              } else {
                arrClustersProcessTmp = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));
                  if(!arrClustersProcessTmp)
                  AliFatal(Form("%sClustersBranch was not found in AliAnalysisTaskGammaCalo! Check the correction framework settings!",fCorrTaskSetting.Data()));
                nclus = arrClustersProcessTmp->GetEntries();
              }
              if(nPhotons==2) fHistoMCPi0GenVsNClus[fiCut]->Fill(particle->Pt(),nclus,fWeightJetJetMC);
              // only check decay pi0->gammagamma
            }
         }
         // Do study of conversions
          if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, i, mcProdVtxX, mcProdVtxY, mcProdVtxZ)){
            if (particle->GetPdgCode() == 22){
              // looking for conversion gammas (electron + positron from pairbuilding (= 5) )
              TParticle* ePos = NULL;
              TParticle* eNeg = NULL;
              if(particle->GetNDaughters() >= 2){
                for(Int_t daughterIndex=particle->GetFirstDaughter();daughterIndex<=particle->GetLastDaughter();daughterIndex++){
                  if(daughterIndex<0) continue;
                  TParticle *tmpDaughter = fMCEvent->Particle(daughterIndex);
                  if(tmpDaughter->GetUniqueID() == 5){
                    if(tmpDaughter->GetPdgCode() == 11){
                      eNeg = tmpDaughter;
                    } else if(tmpDaughter->GetPdgCode() == -11){
                      ePos = tmpDaughter;
                    }
                  }
                }
                if(ePos != NULL && eNeg != NULL){
                  // found both daughters
                  Double_t pT = particle->Pt();
                  Double_t R  = ((TParticle*)fMCEvent->Particle(particle->GetFirstDaughter()))->R();
                  fHistoMCGammaConvRvsPt[fiCut]->Fill(pT,R,fWeightJetJetMC);
                }
              }
            }
          }
      } // end of merging studies

      Int_t isMCFromMBHeader = -1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        isMCFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
        if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
        // Set the jetjet weight to 1 in case the particle orignated from the minimum bias header
        if(isMCFromMBHeader == 2 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4) tempParticleWeight = 1;
      }
      if(!fDoLightOutput){
        if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(particle,fMCEvent)){
          fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt(), tempParticleWeight); // All MC Gamma
          if(particle->GetMother(0) >-1){ // Meson Decay Gamma
            switch(fMCEvent->Particle(particle->GetMother(0))->GetPdgCode()){
            case 111: // Pi0
              fHistoMCDecayGammaPi0Pt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
              break;
            case 113: // Rho0
              fHistoMCDecayGammaRhoPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
              break;
            case 221: // Eta
              fHistoMCDecayGammaEtaPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
              break;
            case 223: // Omega
              fHistoMCDecayGammaOmegaPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
              break;
            case 331: // Eta'
              fHistoMCDecayGammaEtapPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
              break;
            case 333: // Phi
              fHistoMCDecayGammaPhiPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
              break;
            case 3212: // Sigma
              fHistoMCDecayGammaSigmaPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
              break;
            }
          }
        }
      }
      if(fDoMesonAnalysis){

        Double_t mesonY = 1.e30;
        Double_t ratio  = 0;
        if (particle->Energy() != TMath::Abs(particle->Pz())){
          ratio         = (particle->Energy()+particle->Pz()) / (particle->Energy()-particle->Pz());
        }
        if( !(ratio <= 0) ){
          mesonY = particle->Y()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
        }

        if ((mesonY > ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetRapidityCutValueMin()) && (mesonY < ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetRapidityCutValueMax())){
          if ( particle->GetPdgCode() == 211 ){  // positve pions
            fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 0., tempParticleWeight);
          } else if ( particle->GetPdgCode() == -211 ){  // negative pions
            fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 1., tempParticleWeight);
          } else if ( particle->GetPdgCode() == 321 ){  // positve kaons
            fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 2., tempParticleWeight);
          } else if ( particle->GetPdgCode() == -321 ){  // negative kaons
            fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 3., tempParticleWeight);
          } else if ( TMath::Abs(particle->GetPdgCode()) == 310 ){  // K0s
            fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 4., tempParticleWeight);
          } else if ( TMath::Abs(particle->GetPdgCode()) == 130 ){  // K0l
            fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 5., tempParticleWeight);
          } else if ( TMath::Abs(particle->GetPdgCode()) == 3122 ){  // Lambda/ AntiLambda
            fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 6., tempParticleWeight);
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
          Double_t mesonY = 1.e30;
          Double_t ratio  = 0;
          if (particle->Energy() != TMath::Abs(particle->Pz())){
            ratio         = (particle->Energy()+particle->Pz()) / (particle->Energy()-particle->Pz());
          }
          if( !(ratio <= 0) ){
            mesonY = particle->Y()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
          }
          Double_t alpha = -1;
          if (particle->GetPdgCode() == 111 || particle->GetPdgCode() == 221){
            alpha = TMath::Abs((daughter0->Energy() - daughter1->Energy()))/(daughter0->Energy() + daughter1->Energy());
          }

          if(particle->GetPdgCode() == 111){
            if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoMCPi0Pt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // All MC Pi0
            if(fDoJetAnalysis){
              if(fConvJetReader->GetTrueNJets()>0){
                if(!fDoLightOutput) fHistoMCPi0JetEventGenerated[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc in jet event
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
                      if(!fDoLightOutput) fHistoMCPi0inJetGenerated[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc in a jet
                      else fHistoMCPi0Pt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight);
                    }
                  }
                }
                fTrueVectorJetEta.clear();
                fTrueVectorJetPhi.clear();
              }
            }
            fHistoMCPi0WOWeightPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
            if (fIsMC > 1)fHistoMCPi0WOEvtWeightPt[fiCut]->Fill(particle->Pt());
            if (fDoMesonQA > 0 && fDoMesonQA < 3){
              fHistoMCPi0PtY[fiCut]->Fill(particle->Pt(),mesonY,weighted* tempParticleWeight);
              fHistoMCPi0PtAlpha[fiCut]->Fill(particle->Pt(),alpha, tempParticleWeight);
              if (fIsMC == 2) fHistoMCPi0PtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),tempParticleWeight);
            }
          } else if(particle->GetPdgCode() == 221 && !fDoPi0Only){
            if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoMCEtaPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // All MC Eta
            if(fDoJetAnalysis){
              if(fConvJetReader->GetTrueNJets()>0){
                if(!fDoLightOutput) fHistoMCEtaJetEventGenerated[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc in jet event
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
                      if(!fDoLightOutput) fHistoMCEtainJetGenerated[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc in a jet
                      else fHistoMCEtaPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight);
                    }
                  }
                }
                fTrueVectorJetEta.clear();
                fTrueVectorJetPhi.clear();
              }
            }
            fHistoMCEtaWOWeightPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
            if (fIsMC > 1)fHistoMCEtaWOEvtWeightPt[fiCut]->Fill(particle->Pt());
            if (fDoMesonQA > 0 && fDoMesonQA < 3){
              fHistoMCEtaPtY[fiCut]->Fill(particle->Pt(),mesonY,weighted* tempParticleWeight);
              fHistoMCEtaPtAlpha[fiCut]->Fill(particle->Pt(),alpha, tempParticleWeight);
              if (fIsMC == 2) fHistoMCEtaPtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),tempParticleWeight);
            }
          }
          // Check the acceptance for both gammas & whether they are counted as primaries as well
          Bool_t kDaughter0IsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, particle->GetFirstDaughter(), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
          Bool_t kDaughter1IsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, particle->GetLastDaughter(), mcProdVtxX, mcProdVtxY, mcProdVtxZ);

          if( kDaughter0IsPrim && kDaughter1IsPrim &&
            ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter0,fMCEvent) &&
            ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter1,fMCEvent) ){
            if(particle->GetPdgCode() == 111){
              if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoMCPi0InAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc
              if(fIsMC > 1) fHistoMCPi0WOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Pi0 with gamma in acc
              if(fDoJetAnalysis){
                if(fConvJetReader->GetTrueNJets()>0){
                  if(!fDoLightOutput) fHistoMCPi0JetInAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc in jet event
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
                        if(!fDoLightOutput) fHistoMCPi0inJetInAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc in a jet
                        else fHistoMCPi0InAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight);
                      }
                    }
                  }
                  fTrueVectorJetEta.clear();
                  fTrueVectorJetPhi.clear();
                }
              }
            } else if(particle->GetPdgCode() == 221 && !fDoPi0Only){
              if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoMCEtaInAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Eta with gamma in acc
              if(fIsMC > 1) fHistoMCEtaWOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Eta with gamma in acc
              if(fDoJetAnalysis){
                if(fConvJetReader->GetTrueNJets()>0){
                  if(!fDoLightOutput) fHistoMCEtaJetInAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc in jet event
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
                        if(!fDoLightOutput) fHistoMCEtainJetInAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc in a jet
                        else fHistoMCEtaInAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight);
                      }
                    }
                  }
                  fTrueVectorJetEta.clear();
                  fTrueVectorJetPhi.clear();
                }
              }
            }
          }
        }
      }
    // fill secondary histograms
    } else {
      TParticle* particle = (TParticle *)fMCEvent->Particle(i);
      if (!particle) continue;
      Double_t tempParticleWeight       = fWeightJetJetMC;
      Int_t isMCFromMBHeader = -1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        isMCFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
        if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
        // Set the jetjet weight to 1 in case the particle orignated from the minimum bias header
        if(isMCFromMBHeader == 2 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4) tempParticleWeight = 1;
      }

      if(!fDoLightOutput) {
        if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(particle,fMCEvent)){
          if (particle->GetMother(0) > -1 && fMCEvent->Particle(particle->GetMother(0))->GetMother(0) > -1) {
            if (fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 310){
              fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),0.,tempParticleWeight);
            } else if (fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 130) {
              fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),1.,tempParticleWeight);
            } else if (fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 3122) {
              fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),2.,tempParticleWeight);
            } else if (fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 221) {
              fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),3.,tempParticleWeight);
            } else {
              if ( !(TMath::Abs(fMCEvent->Particle(particle->GetMother(0))->GetPdgCode()) == 11 && fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 22) )
                fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),4.,tempParticleWeight);
            }
          } else {
            fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),4.,tempParticleWeight);
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
            fHistoMCSecPi0PtvsSource[fiCut]->Fill(particle->Pt(),source,tempParticleWeight);
            fHistoMCSecPi0Source[fiCut]->Fill(pdgCode);
          } else if(particle->GetPdgCode() == 221 && !fDoPi0Only){
            fHistoMCSecEtaPt[fiCut]->Fill(particle->Pt(),tempParticleWeight);
            fHistoMCSecEtaSource[fiCut]->Fill(pdgCode);
          }

          // check if photons where within acceptance
          if( ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter0,fMCEvent) &&
                ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter1,fMCEvent)){
            if(particle->GetPdgCode() == 111){
              Int_t source = GetSourceClassification(111,pdgCode);
              fHistoMCSecPi0InAccPtvsSource[fiCut]->Fill(particle->Pt(),source,tempParticleWeight);
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

        Double_t tempPi0CandWeight       = fWeightJetJetMC;
        // Set the pi0 candidate jetjet weight to 1 in case both photons orignated from the minimum bias header
        if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4){
          if( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma0->GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2 &&
              ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma1->GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2)
            tempPi0CandWeight = 1;
        }


        if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 2){
          if ( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsSpecialTrigger()==6 ){
            if (fCaloTriggerMimicHelper[fiCut]){
              fHistoGoodMesonClusters[fiCut]->Fill(1); //"All Meson Candidates"
              if ( !((fCaloTriggerMimicHelper[fiCut]->IsClusterIDTriggered(gamma0->GetCaloClusterRef()))||(fCaloTriggerMimicHelper[fiCut]->IsClusterIDTriggered(gamma1->GetCaloClusterRef()))) ){
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
              if (fCaloTriggerMimicHelper[fiCut]->IsClusterIDTriggered(gamma0->GetCaloClusterRef())){ //gamma 0 is triggered
                ClusterIDIsInBadDDL=fCaloTriggerMimicHelper[fiCut]->IsTriggeredClusterIDInBadDDL(gamma0->GetCaloClusterRef());
                if (ClusterIDIsInBadDDL>=DDLIsBadIndex){ //DDL is bad
                    fHistoGoodMesonClusters[fiCut]->Fill(7); //"DDL not passed"
                    continue;
                }
              }
              if (fCaloTriggerMimicHelper[fiCut]->IsClusterIDTriggered(gamma1->GetCaloClusterRef())) { //gamma 1 is triggered
                ClusterIDIsInBadDDL=fCaloTriggerMimicHelper[fiCut]->IsTriggeredClusterIDInBadDDL(gamma1->GetCaloClusterRef());
                if (ClusterIDIsInBadDDL>=DDLIsBadIndex){ //DDL is bad
                    fHistoGoodMesonClusters[fiCut]->Fill(7); //"DDL not passed"
                    continue;
                }
              }
              fHistoGoodMesonClusters[fiCut]->Fill(6); //"DDL passed"
            }
          }
        }

        if( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoGammaMinEnergyCut() ){
          Int_t minDaughters        = ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetNDaughterEnergyCut();
          Float_t minDaughterEnergy = ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetSingleDaughterMinE();
          if(minDaughters==1){ // at least one over threshold
             if( (gamma0->E() < minDaughterEnergy)  && (gamma1->E() < minDaughterEnergy)) {
                 if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 2){
                   if ( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsSpecialTrigger()==6 ){
                     fHistoGoodMesonClusters[fiCut]->Fill(5); //"Cluster E not passed"
                   }
                 }
                 continue;
             } else if ( (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 2)&&(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsSpecialTrigger()==6) ){
                 if ( (gamma0->E() < minDaughterEnergy) ){
                     if (!(fCaloTriggerMimicHelper[fiCut]->IsClusterIDTriggered(gamma1->GetCaloClusterRef()))){
                         if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 2){
                           if ( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsSpecialTrigger()==6 ){
                             fHistoGoodMesonClusters[fiCut]->Fill(5); //"Cluster E not passed"
                           }
                         }
                         continue;
                     }
                 } else if ( (gamma1->E() < minDaughterEnergy) ){
                     if (!(fCaloTriggerMimicHelper[fiCut]->IsClusterIDTriggered(gamma0->GetCaloClusterRef()))){
                         if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 2){
                           if ( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsSpecialTrigger()==6 ){
                             fHistoGoodMesonClusters[fiCut]->Fill(5); //"Cluster E not passed"
                           }
                         }
                         continue;
                     }
                 }
             }
          } else if (minDaughters==2){ // both over threshold
             if( (gamma0->E() < minDaughterEnergy)  || (gamma1->E() < minDaughterEnergy)) {
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

        if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetDoSecondaryTrackMatching()){
          TClonesArray * arrClustersProcess = NULL;
          arrClustersProcess = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));
          AliVCluster* Cluster0 = NULL;
          AliVCluster* Cluster1 = NULL;
          if (gamma0->GetIsCaloPhoton() > 0 && gamma1->GetIsCaloPhoton() > 0){
            if(fInputEvent->IsA()==AliESDEvent::Class()){
              if(arrClustersProcess){
                Cluster0 = new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersProcess->At(gamma0->GetCaloClusterRef()));
                Cluster1 = new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersProcess->At(gamma1->GetCaloClusterRef()));
              }else{
                Cluster0 = fInputEvent->GetCaloCluster(gamma0->GetCaloClusterRef());
                Cluster1 = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef());
              }
            } else if(fInputEvent->IsA()==AliAODEvent::Class()){
              if(arrClustersProcess){
                Cluster0 = new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersProcess->At(gamma0->GetCaloClusterRef()));
                Cluster1 = new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersProcess->At(gamma1->GetCaloClusterRef()));
              } else{
                Cluster0 = fInputEvent->GetCaloCluster(gamma0->GetCaloClusterRef());
                Cluster1 = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef());
              }
            }
          }
          Bool_t ClusterMatched = kFALSE;
          for(Int_t ConversionIndex=0;ConversionIndex<fGammaCandidates->GetEntries();ConversionIndex++){
            AliAODConversionPhoton *gammaConversion=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(ConversionIndex));
            if (gammaConversion==NULL) continue;
            Bool_t matchedGamma0 =  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gammaConversion, Cluster0, fInputEvent, tempPi0CandWeight);
            Bool_t matchedGamma1 =  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gammaConversion, Cluster1, fInputEvent, tempPi0CandWeight);
            if(matchedGamma0 || matchedGamma1) {
              ClusterMatched = kTRUE;
            }
          }
          if(arrClustersProcess){
              delete Cluster0;
              delete Cluster1;
          }
          if(ClusterMatched) continue;
        }
        if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),gamma0->GetLeadingCellID(),gamma1->GetLeadingCellID(), gamma0->GetIsCaloPhoton(), gamma1->GetIsCaloPhoton()))){
          if(fLocalDebugFlag == 1) DebugMethodPrint1(pi0cand,gamma0,gamma1);
          if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoMotherInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(), tempPi0CandWeight);
          if(fIsMC && fDoTrueSphericity){
            if(fV0Reader->GetSphericity() != -1 && fV0Reader->GetSphericity() != 0 && fV0Reader->GetSphericityTrue() != -1 && fV0Reader->GetSphericityTrue() != 0){
              fRecSph = fV0Reader->GetSphericity();
              fTrueSph = fV0Reader->GetSphericityTrue();
              fPi0Pt = pi0cand->Pt();
              fPi0InvMass = pi0cand->M();
              tTreeSphericity[fiCut]->Fill();
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
              if(!fDoLightOutput)  fHistoJetMotherInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(), tempPi0CandWeight);
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
                      fHistoEtaPhiJetPi0Cand[fiCut]->Fill(DeltaPhi,DeltaEta, tempPi0CandWeight);
                      fHistoRJetPi0Cand[fiCut]->Fill(RJetPi0Cand,pi0cand->Pt(), tempPi0CandWeight);
                  }
                  if(fConvJetReader->Get_Jet_Radius() > 0 ){
                    if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                      counter ++;
                      if(!fDoLightOutput){
                          fHistoPi0InJetMotherInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(), tempPi0CandWeight);
                          fHistoEtaPhiJetWithPi0Cand[fiCut]->Fill(DeltaPhi, DeltaEta, tempPi0CandWeight);
                          Double_t PtRatio = pi0cand->Pt()/(fVectorJetPt.at(i));
                          fHistoJetPi0PtRatio[fiCut]->Fill(PtRatio);
                          Double_t dotproduct = fVectorJetPx.at(i)*pi0cand->Px() + fVectorJetPy.at(i)*pi0cand->Py() + fVectorJetPz.at(i)*pi0cand->Pz();
                          Double_t magn = pow(fVectorJetPx.at(i),2) + pow(fVectorJetPy.at(i),2) + pow(fVectorJetPz.at(i),2);
                          Double_t z = dotproduct/magn;
                          fHistoJetFragmFunc[fiCut]->Fill(z,fVectorJetPt.at(i), tempPi0CandWeight);
                          fHistoJetFragmFuncZInvMass[fiCut]->Fill(pi0cand->M(), z, tempPi0CandWeight);
                      }
                      else fHistoMotherInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(), tempPi0CandWeight);

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
                            if(fVectorJetPt.at(i) >= 10){
                              fHistoUnfoldingAsData[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(), tempPi0CandWeight);
                              fHistoUnfoldingAsDataInvMassZ[fiCut]->Fill(pi0cand->M(), z, tempPi0CandWeight);
                            }
                            if(fVectorJetPt.at(i) < 10 && fTrueVectorJetPt.at(match) >= 10){
                              fHistoUnfoldingMissed[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(), tempPi0CandWeight);
                              fHistoUnfoldingMissedInvMassZ[fiCut]->Fill(pi0cand->M(), z, tempPi0CandWeight);
                            }
                            if(fVectorJetPt.at(i) >= 10 && fTrueVectorJetPt.at(match) < 10){
                              fHistoUnfoldingReject[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(), tempPi0CandWeight);
                              fHistoUnfoldingRejectInvMassZ[fiCut]->Fill(pi0cand->M(), z, tempPi0CandWeight);
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
              fVectorJetPt.clear();
              fVectorJetPx.clear();
              fVectorJetPy.clear();
              fVectorJetPz.clear();
              fVectorJetEta.clear();
              fVectorJetPhi.clear();
            }
          }
          // fill new histograms
          if((!fDoLightOutput || fDoPi0Only || fDoECalibOutput) && TMath::Abs(pi0cand->GetAlpha())<0.1){
            fHistoMotherInvMassECalib[fiCut]->Fill(pi0cand->M(),pi0cand->E(),tempPi0CandWeight);
          }

          if (fDoMesonQA > 0 && fDoMesonQA < 3){
            if ( pi0cand->M() > 0.05 && pi0cand->M() < 0.17){
              fHistoMotherPi0PtY[fiCut]->Fill(pi0cand->Pt(),pi0cand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), tempPi0CandWeight);
              fHistoMotherPi0PtAlpha[fiCut]->Fill(pi0cand->Pt(),TMath::Abs(pi0cand->GetAlpha()), tempPi0CandWeight);
              fHistoMotherPi0PtOpenAngle[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetOpeningAngle(), tempPi0CandWeight);
              fHistoMotherPi0NGoodESDTracksPt[fiCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks(),pi0cand->Pt(), tempPi0CandWeight);
            }
            if ( pi0cand->M() > 0.45 && pi0cand->M() < 0.65 && !fDoPi0Only){
              fHistoMotherEtaPtY[fiCut]->Fill(pi0cand->Pt(),pi0cand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), tempPi0CandWeight);
              fHistoMotherEtaPtAlpha[fiCut]->Fill(pi0cand->Pt(),TMath::Abs(pi0cand->GetAlpha()), tempPi0CandWeight);
              fHistoMotherEtaPtOpenAngle[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetOpeningAngle(),  tempPi0CandWeight);
              fHistoMotherEtaNGoodESDTracksPt[fiCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks(),pi0cand->Pt(), tempPi0CandWeight);
            }
          }
          if (fDoMesonQA == 2){
            fHistoMotherPtOpenAngle[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetOpeningAngle(), tempPi0CandWeight);
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

          if ( (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 2 && fDoMesonQA == 5 && fIsMC == 0 && ((pi0cand->M() > 0.05 && pi0cand->M() < 0.17) || (pi0cand->Pt() > 15. && pi0cand->M() > 0.45 && pi0cand->M() < 0.65 )))
                || ((((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 4 || ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 1 || ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 3) && fDoMesonQA == 5 && fIsMC == 0 && (pi0cand->M() > 0.0 && pi0cand->M() < 0.3) ) ) {

            TClonesArray * arrClustersProcess = NULL;
            arrClustersProcess = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));
            AliVCluster* Cluster0 = NULL;
            AliVCluster* Cluster1 = NULL;
            if (gamma0->GetIsCaloPhoton() > 0 && gamma1->GetIsCaloPhoton() > 0){
              if(fInputEvent->IsA()==AliESDEvent::Class()){
                if(arrClustersProcess){
                  Cluster0 = new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersProcess->At(gamma0->GetCaloClusterRef()));
                  Cluster1 = new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersProcess->At(gamma1->GetCaloClusterRef()));
                }else{
                  Cluster0 = fInputEvent->GetCaloCluster(gamma0->GetCaloClusterRef());
                  Cluster1 = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef());
                }
              } else if(fInputEvent->IsA()==AliAODEvent::Class()){
                if(arrClustersProcess){
                  Cluster0 = new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersProcess->At(gamma0->GetCaloClusterRef()));
                  Cluster1 = new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersProcess->At(gamma1->GetCaloClusterRef()));
                } else{
                  Cluster0 = fInputEvent->GetCaloCluster(gamma0->GetCaloClusterRef());
                  Cluster1 = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef());
                }
              }
            }
            fInvMassTreeInvMass = pi0cand->M();
            fInvMassTreePt      = pi0cand->Pt();
            fClusterTimeProbe   = Cluster1->GetTOF();
            fClusterTimeTag     = Cluster0->GetTOF();
            fClusterEProbe   = Cluster1->E();
            fClusterETag     = Cluster0->E();
            tClusterTimingEff[fiCut]->Fill();
          }

          if(fIsMC> 0){
            if(fInputEvent->IsA()==AliESDEvent::Class())
              ProcessTrueMesonCandidates(pi0cand,gamma0,gamma1);
            if(fInputEvent->IsA()==AliAODEvent::Class())
              ProcessTrueMesonCandidatesAOD(pi0cand,gamma0,gamma1);
          }

          if((pi0cand->GetOpeningAngle() < 0.017) && (pi0cand->Pt() > 15.) && fDoClusterQA > 1){
            if (fCloseHighPtClusters == NULL){
              fCloseHighPtClusters = new TObjString(Form("%s", ((TString)fV0Reader->GetCurrentFileName()).Data()));
              if (tClusterQATree) tClusterQATree->Fill();
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
void AliAnalysisTaskGammaCalo::ProcessTrueMesonCandidates(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1)
{
  // Process True Mesons
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();
  fRconv = -1;

  Double_t tempTruePi0CandWeight       = fWeightJetJetMC;

  Bool_t isTruePi0         = kFALSE;
  Bool_t isTrueEta         = kFALSE;
  //Bool_t isTrueGamma        = kFALSE;
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
  if (TrueGammaCandidate1->GetIsCaloPhoton() == 0) AliFatal("CaloPhotonFlag has not been set. Aborting");

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

  // Set the pi0 candidate jetjet weight to 1 in case both photons orignated from the minimum bias header
  if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4){
    if( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(TrueGammaCandidate0->GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2 &&
        ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(TrueGammaCandidate1->GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2)
      tempTruePi0CandWeight = 1;
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

    TClonesArray * arrClustersMesonCand = NULL;
    if(fCorrTaskSetting.CompareTo(""))
      arrClustersMesonCand = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));

    AliVCluster* clus1 = NULL;
    if(fInputEvent->IsA()==AliESDEvent::Class()){
      if(arrClustersMesonCand)
        clus1 = new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersMesonCand->At(TrueGammaCandidate0->GetCaloClusterRef()));
      else
        clus1 = fInputEvent->GetCaloCluster(TrueGammaCandidate0->GetCaloClusterRef());
    } else if(fInputEvent->IsA()==AliAODEvent::Class()){
      if(arrClustersMesonCand)
        clus1 = new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersMesonCand->At(TrueGammaCandidate0->GetCaloClusterRef()));
      else
        clus1 = fInputEvent->GetCaloCluster(TrueGammaCandidate0->GetCaloClusterRef());
    }

    if (!clus1) return;
    TLorentzVector clusterVector1;
    clus1->GetMomentum(clusterVector1,vertex);
    if(arrClustersMesonCand) delete clus1;
    TLorentzVector* tmpvec1 = new TLorentzVector();
    tmpvec1->SetPxPyPzE(clusterVector1.Px(),clusterVector1.Py(),clusterVector1.Pz(),clusterVector1.E());
    // convert to AODConversionPhoton
    AliAODConversionPhoton *PhotonCandidate1=new AliAODConversionPhoton(tmpvec1);
    if(!PhotonCandidate1) return;


    AliVCluster* clus2 = NULL;
    if(fInputEvent->IsA()==AliESDEvent::Class()){
      if(arrClustersMesonCand)
        clus2 = new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersMesonCand->At(TrueGammaCandidate1->GetCaloClusterRef()));
      else
        clus2 = fInputEvent->GetCaloCluster(TrueGammaCandidate1->GetCaloClusterRef());
    } else if(fInputEvent->IsA()==AliAODEvent::Class()){
      if(arrClustersMesonCand)
        clus2 = new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersMesonCand->At(TrueGammaCandidate1->GetCaloClusterRef()));
      else
        clus2 = fInputEvent->GetCaloCluster(TrueGammaCandidate1->GetCaloClusterRef());
    }

    if (!clus2) return;
    TLorentzVector clusterVector2;
    clus2->GetMomentum(clusterVector2,vertex);
    if(arrClustersMesonCand) delete clus2;
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
      fHistoTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),  tempTruePi0CandWeight);
      if(previouslyNotFoundTrueMesons && !fDoLightOutput && fHistoTruePi0InvMassPtAdditional[fiCut]) fHistoTruePi0InvMassPtAdditional[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),tempTruePi0CandWeight);
      if(fDoJetAnalysis && !fDoLightOutput){
        if(fConvJetReader->GetTrueNJets()>0){
          fHistoTruePi0JetMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),  tempTruePi0CandWeight);
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
                fHistoTruePi0InJetMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), tempTruePi0CandWeight);
                fHistoMotherPi0inJetPtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), tempTruePi0CandWeight);
                fHistoMotherPi0inJetPtPhi[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Phi(), tempTruePi0CandWeight);
              }
            }
            fHistoTrueDoubleCountingPi0Jet[fiCut]->Fill(counter);
          }
          fTrueVectorJetEta.clear();
          fTrueVectorJetPhi.clear();
        }
      }

//      fHistoTruePi0NonLinearity[fiCut]->Fill(TrueGammaCandidate0->E(),gammaMC0->Energy()/TrueGammaCandidate0->E());
//      fHistoTruePi0NonLinearity[fiCut]->Fill(TrueGammaCandidate1->E(),gammaMC1->Energy()/TrueGammaCandidate1->E());
      if (!fDoLightOutput && TMath::Abs(Pi0Candidate->GetAlpha())< 0.1){
        fHistoTruePi0InvMassPtAlpha[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),  tempTruePi0CandWeight);
        if (TrueGammaCandidate1->IsLargestComponentPhoton() && TrueGammaCandidate0->IsLargestComponentPhoton())
          fHistoTruePi0PureGammaInvMassPtAlpha[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),  tempTruePi0CandWeight);
      }
    }
    if (isTrueEta && !fDoPi0Only){
      fHistoTrueEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),  tempTruePi0CandWeight);
      if(previouslyNotFoundTrueMesons && !fDoLightOutput && fHistoTrueEtaInvMassPtAdditional[fiCut]) fHistoTrueEtaInvMassPtAdditional[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),tempTruePi0CandWeight);
      if(fDoJetAnalysis && !fDoLightOutput){
        if(fConvJetReader->GetTrueNJets()>0){
          fHistoTrueEtaJetMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),  tempTruePi0CandWeight);
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
                fHistoTrueEtaInJetMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), tempTruePi0CandWeight);
                fHistoMotherEtainJetPtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), tempTruePi0CandWeight);
                fHistoMotherEtainJetPtPhi[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Phi(), tempTruePi0CandWeight);
              }
            }
            fHistoTrueDoubleCountingEtaJet[fiCut]->Fill(counter);
          }
          fTrueVectorJetEta.clear();
          fTrueVectorJetPhi.clear();
        }
      }
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
        if (isTrueEta && !fDoPi0Only) fHistoTrueEtaCaloPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
      // both particles are electrons
      if (TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate1->IsLargestComponentElectron() ) {
        if (isTruePi0) fHistoTruePi0CaloElectronInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        if (isTrueEta && !fDoPi0Only) fHistoTrueEtaCaloElectronInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
      // both particles are converted electrons
      if ((TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate0->IsConversion()) && (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion()) ){
        if (isTruePi0 )fHistoTruePi0CaloConvertedPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        if (isTrueEta && !fDoPi0Only)fHistoTrueEtaCaloConvertedPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
      // 1 gamma is converted the other one is normals
      if ( (TrueGammaCandidate0->IsLargestComponentPhoton() && TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion()) ||
        (TrueGammaCandidate1->IsLargestComponentPhoton() && TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate0->IsConversion())
      ) {
        if (isTruePi0) fHistoTruePi0CaloMixedPhotonConvPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        if (isTrueEta && !fDoPi0Only) fHistoTrueEtaCaloMixedPhotonConvPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
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
        if (isTrueEta && !fDoPi0Only) fHistoTruePi0CaloMergedClusterInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
      // at least one of the photon is merged and part conv
      if (TrueGammaCandidate1->IsMergedPartConv() || TrueGammaCandidate0->IsMergedPartConv()) {
        if (isTruePi0) fHistoTruePi0CaloMergedClusterPartConvInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        if (isTrueEta && !fDoPi0Only) fHistoTrueEtaCaloMergedClusterPartConvInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
    }

    if (fDoMesonQA == 2 && fIsMC < 2){
      // category 1: 2 real photons unmerged
      if (TrueGammaCandidate0->IsLargestComponentPhoton() && !TrueGammaCandidate0->IsMerged() && TrueGammaCandidate1->IsLargestComponentPhoton() && !TrueGammaCandidate1->IsMerged()) {
        if (isTruePi0) fHistoTruePi0Category1[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        if (isTrueEta && !fDoPi0Only) fHistoTrueEtaCategory1[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
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
        if (isTrueEta && !fDoPi0Only){
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
        if (isTrueEta && !fDoPi0Only){
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
        if (isTrueEta && !fDoPi0Only){
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
          fHistoTruePi0PtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), tempTruePi0CandWeight);
          fHistoTruePi0PtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),TMath::Abs(Pi0Candidate->GetAlpha()), tempTruePi0CandWeight);
          fHistoTruePi0PtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle(), tempTruePi0CandWeight);
        }
      } else if (isTrueEta && !fDoPi0Only){
        if ( Pi0Candidate->M() > 0.45 && Pi0Candidate->M() < 0.65){
          fHistoTrueEtaPtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), tempTruePi0CandWeight);
          fHistoTrueEtaPtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),TMath::Abs(Pi0Candidate->GetAlpha()), tempTruePi0CandWeight);
          fHistoTrueEtaPtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle(), tempTruePi0CandWeight);
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
      if (isTruePi0) {
        if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTrueSecondaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
        if(fDoJetAnalysis){
          if(fConvJetReader->GetTrueNJets()>0){
            if(!fDoLightOutput) fHistoTrueSecondaryPi0InvJetMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
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
                    if(!fDoLightOutput) fHistoTrueSecondaryPi0InvinJetMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
                    else fHistoTrueSecondaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
                  }
                }
              }
            fTrueVectorJetEta.clear();
            fTrueVectorJetPhi.clear();
          }
        }
      }
      if (secMotherLabel >-1){
        if(fMCEvent->Particle(secMotherLabel)->GetPdgCode()==310 && isTruePi0){
          if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTrueSecondaryPi0FromK0sInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
          if (fDoMesonQA > 0 && fDoMesonQA < 3 && fIsMC < 2)fHistoTrueK0sWithPi0DaughterMCPt[fiCut]->Fill(fMCEvent->Particle(secMotherLabel)->Pt(), tempTruePi0CandWeight);
          if (fDoMesonQA == 3 ){
            iFlag     = 3;
            if (!isSameConvertedGamma)tTrueInvMassROpenABPtFlag[fiCut]->Fill();
          }
          if(fDoJetAnalysis){
            if(fConvJetReader->GetTrueNJets()>0){
              if(!fDoLightOutput) fHistoTrueSecondaryPi0FromK0sJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
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
                      if(!fDoLightOutput) fHistoTrueSecondaryPi0FromK0sinJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
                      else fHistoTrueSecondaryPi0FromK0sInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
                    }
                  }
                }
              fTrueVectorJetEta.clear();
              fTrueVectorJetPhi.clear();
            }
          }
        } else if(fMCEvent->Particle(secMotherLabel)->GetPdgCode()==130 && isTruePi0){
          if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTrueSecondaryPi0FromK0lInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
          if (fDoMesonQA > 0 && fDoMesonQA < 3 && fIsMC < 2)fHistoTrueK0lWithPi0DaughterMCPt[fiCut]->Fill(fMCEvent->Particle(secMotherLabel)->Pt(), tempTruePi0CandWeight);
          if(fDoJetAnalysis){
            if(fConvJetReader->GetTrueNJets()>0){
              if(!fDoLightOutput) fHistoTrueSecondaryPi0FromK0lJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
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
                      if(!fDoLightOutput) fHistoTrueSecondaryPi0FromK0linJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
                      else fHistoTrueSecondaryPi0FromK0lInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
                    }
                  }
                }
              fTrueVectorJetEta.clear();
              fTrueVectorJetPhi.clear();
            }
          }

        } else if(fMCEvent->Particle(secMotherLabel)->GetPdgCode()==221 && isTruePi0){
          fHistoTrueSecondaryPi0FromEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
          if (fDoMesonQA > 0 && fDoMesonQA < 3 && fIsMC < 2)fHistoTrueEtaWithPi0DaughterMCPt[fiCut]->Fill(fMCEvent->Particle(secMotherLabel)->Pt(), tempTruePi0CandWeight);
          if (fDoMesonQA == 3 ){
            iFlag     = 4;
            if (!isSameConvertedGamma)tTrueInvMassROpenABPtFlag[fiCut]->Fill();
          }
        } else if(fMCEvent->Particle(secMotherLabel)->GetPdgCode()==3122 && isTruePi0){
          if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTrueSecondaryPi0FromLambdaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
          if (fDoMesonQA > 0 && fDoMesonQA < 3 && fIsMC < 2)fHistoTrueLambdaWithPi0DaughterMCPt[fiCut]->Fill(fMCEvent->Particle(secMotherLabel)->Pt(), tempTruePi0CandWeight);
          if (fDoMesonQA == 3 ){
            iFlag     = 5;
            if (!isSameConvertedGamma)tTrueInvMassROpenABPtFlag[fiCut]->Fill();
          }
          if(fDoJetAnalysis){
            if(fConvJetReader->GetTrueNJets()>0){
              if(!fDoLightOutput) fHistoTrueSecondaryPi0FromLambdaJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
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
                      if(!fDoLightOutput) fHistoTrueSecondaryPi0FromLambdainJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
                      else fHistoTrueSecondaryPi0FromLambdaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
                    }
                  }
                }
              fTrueVectorJetEta.clear();
              fTrueVectorJetPhi.clear();
            }
          }
        } else if (isTruePi0){
          if (fDoMesonQA == 3 ){
            iFlag     = 6;
            if (!isSameConvertedGamma)tTrueInvMassROpenABPtFlag[fiCut]->Fill();
          }
        } else if (isTrueEta && !fDoPi0Only){
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

        if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTruePrimaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* tempTruePi0CandWeight);
        fHistoTruePrimaryPi0W0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), tempTruePi0CandWeight);
        fProfileTruePrimaryPi0WeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* tempTruePi0CandWeight);
        if (CheckVectorForDoubleCount(fVectorDoubleCountTruePi0s,gamma0MotherLabel)) fHistoDoubleCountTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), weighted*tempTruePi0CandWeight);
        if(fDoJetAnalysis){
          Int_t MCLabelGamma0 = TrueGammaCandidate0->GetCaloPhotonMCLabel(0);
          Int_t MCLabelGamma1 = TrueGammaCandidate1->GetCaloPhotonMCLabel(0);
          TParticle * gamma0 = (TParticle*)fMCEvent->Particle(MCLabelGamma0);
          TParticle * gamma1 = (TParticle*)fMCEvent->Particle(MCLabelGamma1);
          Double_t TruePtPi0 = TMath::Sqrt(pow((gamma0->Px() + gamma1->Px()),2)+pow((gamma0->Py() + gamma1->Py()),2));
          if(fConvJetReader->GetTrueNJets()>0){
            if(!fDoLightOutput) fHistoTruePrimaryPi0JetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* tempTruePi0CandWeight);
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
                      fHistoTruePrimaryPi0inJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* tempTruePi0CandWeight);
                      fHistoTruePrimaryPi0InJetInvMassTruePt[fiCut]->Fill(Pi0Candidate->M(),TruePtPi0, tempTruePi0CandWeight);
                      Double_t dotproduct = fTrueVectorJetPx.at(i)*Pi0Candidate->Px() + fTrueVectorJetPy.at(i)*Pi0Candidate->Py() + fTrueVectorJetPz.at(i)*Pi0Candidate->Pz();
                      Double_t magn = pow(fTrueVectorJetPx.at(i),2) + pow(fTrueVectorJetPy.at(i),2) + pow(fTrueVectorJetPz.at(i),2);
                      Double_t z = dotproduct/magn;
                      fHistoTruePi0JetFragmFunc[fiCut]->Fill(z,fTrueVectorJetPt.at(i),weighted* tempTruePi0CandWeight);
                      fHistoTruePi0JetFragmFuncZInvMass[fiCut]->Fill(Pi0Candidate->M(), z, weighted* tempTruePi0CandWeight);
                    }
                    else fHistoTruePrimaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* tempTruePi0CandWeight);
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
      } else if (isTrueEta && !fDoPi0Only){
        if (fDoMesonQA == 3 ){
          iFlag     = 2;
          if (!isSameConvertedGamma)tTrueInvMassROpenABPtFlag[fiCut]->Fill();
        }

        if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTruePrimaryEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* tempTruePi0CandWeight);
        fHistoTruePrimaryEtaW0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), tempTruePi0CandWeight);
        fProfileTruePrimaryEtaWeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* tempTruePi0CandWeight);
        if (CheckVectorForDoubleCount(fVectorDoubleCountTrueEtas,gamma0MotherLabel)) fHistoDoubleCountTrueEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), weighted*tempTruePi0CandWeight);
        if(fDoJetAnalysis){
          if(fConvJetReader->GetTrueNJets()>0){
            if(!fDoLightOutput && !fDoPi0Only) fHistoTruePrimaryEtaJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* tempTruePi0CandWeight);
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
                if(fConvJetReader->Get_Jet_Radius() > 0 && !fDoPi0Only){
                  if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                    if(!fDoLightOutput){
                      fHistoTruePrimaryEtainJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* tempTruePi0CandWeight);
                      Double_t dotproduct = fTrueVectorJetPx.at(i)*Pi0Candidate->Px() + fTrueVectorJetPy.at(i)*Pi0Candidate->Py() + fTrueVectorJetPz.at(i)*Pi0Candidate->Pz();
                      Double_t magn = pow(fTrueVectorJetPx.at(i),2) + pow(fTrueVectorJetPy.at(i),2) + pow(fTrueVectorJetPz.at(i),2);
                      Double_t z = dotproduct/magn;
                      fHistoTrueEtaJetFragmFunc[fiCut]->Fill(z,fTrueVectorJetPt.at(i), weighted* tempTruePi0CandWeight);
                      fHistoTrueEtaJetFragmFuncZInvMass[fiCut]->Fill(Pi0Candidate->M(), z, weighted* tempTruePi0CandWeight);
                    }else fHistoTruePrimaryEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* tempTruePi0CandWeight);
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
      }

      if (fDoMesonQA > 0 && fDoMesonQA < 3 && fIsMC<2){
        if(isTruePi0){ // Only primary pi0 for resolution
          fHistoTruePrimaryPi0MCPtResolPt[fiCut]->Fill(((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt(),(Pi0Candidate->Pt()-((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt())/((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt(),weighted* tempTruePi0CandWeight);
        }
        if (isTrueEta && !fDoPi0Only){ // Only primary eta for resolution
          fHistoTruePrimaryEtaMCPtResolPt[fiCut]->Fill(((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt(),(Pi0Candidate->Pt()-((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt())/((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt(),weighted* tempTruePi0CandWeight);
        }
      }
    }
  } else if(!isTruePi0 && !isTrueEta){ // Background
    if (fDoMesonQA > 0 && fDoMesonQA < 3){
      if(gamma0MotherLabel>-1 && gamma1MotherLabel>-1){ // Both Tracks are Photons and have a mother but not Pi0 or Eta
        fHistoTrueBckGGInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), tempTruePi0CandWeight);

        if( ((((TParticle*)fMCEvent->Particle(gamma0MotherLabel))->GetPdgCode() == 111
            || ((TParticle*)fMCEvent->Particle(gamma0MotherLabel))->GetPdgCode() == 221)
            && (TrueGammaCandidate0->IsMerged() || TrueGammaCandidate0->IsMergedPartConv()))
            ||
            ((((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->GetPdgCode() == 111
            || ((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->GetPdgCode() == 221)
            && (TrueGammaCandidate1->IsMerged() || TrueGammaCandidate1->IsMergedPartConv()))
        ){
          fHistoTrueBckFullMesonContainedInOneClusterInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), tempTruePi0CandWeight);
        }else if( (TrueGammaCandidate0->E()/Pi0Candidate->E() > 0.7) || (TrueGammaCandidate1->E()/Pi0Candidate->E() > 0.7) ){
          fHistoTrueBckAsymEClustersInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), tempTruePi0CandWeight);
        }
      } else { // No photon or without mother
        fHistoTrueBckContInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), tempTruePi0CandWeight);
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

  Double_t tempTruePi0CandWeight       = fWeightJetJetMC;

  // Process True Mesons
  if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (fAODMCTrackArray == NULL) return;

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
    gammaMC0 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma0MCLabel));
    if (TrueGammaCandidate0->IsLargestComponentPhoton() || TrueGammaCandidate0->IsLargestComponentElectron()){    // largest component is electro magnetic
      // get mother of interest (pi0 or eta)
      if (TrueGammaCandidate0->IsLargestComponentPhoton()){                            // for photons its the direct mother
        gamma0MotherLabel=gammaMC0->GetMother();
      } else if (TrueGammaCandidate0->IsLargestComponentElectron()){                         // for electrons its either the direct mother or for conversions the grandmother
        if (TrueGammaCandidate0->IsConversion()){
          convertedPhotonLabel0 = gammaMC0->GetMother();
          AliAODMCParticle * gammaGrandMotherMC0 =  static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gammaMC0->GetMother()));
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
    gammaMC1 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MCLabel));
    if (TrueGammaCandidate1->IsLargestComponentPhoton() || TrueGammaCandidate1->IsLargestComponentElectron()){    // largest component is electro magnetic
      // get mother of interest (pi0 or eta)
      if (TrueGammaCandidate1->IsLargestComponentPhoton()){                            // for photons its the direct mother
        gamma1MotherLabel=gammaMC1->GetMother();
      } else if (TrueGammaCandidate1->IsLargestComponentElectron()){                         // for electrons its either the direct mother or for conversions the grandmother
        if (TrueGammaCandidate1->IsConversion()){
          convertedPhotonLabel1 = gammaMC1->GetMother();
          AliAODMCParticle * gammaGrandMotherMC1 =  static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gammaMC1->GetMother()));
          gamma1MotherLabel=gammaGrandMotherMC1->GetMother();
        } else gamma1MotherLabel=gammaMC1->GetMother();
      }
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

  // Set the pi0 candidate jetjet weight to 1 in case both photons orignated from the minimum bias header
  if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4){
    if( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(TrueGammaCandidate0->GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2 &&
        ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(TrueGammaCandidate1->GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2)
      tempTruePi0CandWeight = 1;
  }

  if (convertedPhotonLabel0 > -1 && convertedPhotonLabel1 > -1){
    if (convertedPhotonLabel0==convertedPhotonLabel1){
      isSameConvertedGamma = kTRUE;
      if (fDoMesonQA == 3 ){
//         fHistoGammaOpenAngleInvMassPt[fiCut]->Fill(Pi0Candidate->GetOpeningAngle(),Pi0Candidate->M(),Pi0Candidate->Pt(),1);
      }
    }
  }


  if(isTruePi0 || isTrueEta){// True Pion or Eta
    if(isTruePi0){
      fHistoTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), tempTruePi0CandWeight);
      if(previouslyNotFoundTrueMesons && !fDoLightOutput && fHistoTruePi0InvMassPtAdditional[fiCut]) fHistoTruePi0InvMassPtAdditional[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),tempTruePi0CandWeight);
      if(fDoJetAnalysis && !fDoLightOutput){
        if(fConvJetReader->GetTrueNJets()>0){
          fHistoTruePi0JetMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),  tempTruePi0CandWeight);
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
                fHistoTruePi0InJetMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), tempTruePi0CandWeight);
                fHistoMotherPi0inJetPtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), tempTruePi0CandWeight);
                fHistoMotherPi0inJetPtPhi[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Phi(), tempTruePi0CandWeight);
              }
            }
            fHistoTrueDoubleCountingPi0Jet[fiCut]->Fill(counter);
          }
          fTrueVectorJetEta.clear();
          fTrueVectorJetPhi.clear();
        }
      }

//      fHistoTruePi0NonLinearity[fiCut]->Fill(TrueGammaCandidate0->E(),gammaMC0->E()/TrueGammaCandidate0->E());
//      fHistoTruePi0NonLinearity[fiCut]->Fill(TrueGammaCandidate1->E(),gammaMC1->E()/TrueGammaCandidate1->E());
      if (!fDoLightOutput && TMath::Abs(Pi0Candidate->GetAlpha())< 0.1){
        fHistoTruePi0InvMassPtAlpha[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),  tempTruePi0CandWeight);
        if (TrueGammaCandidate1->IsLargestComponentPhoton() && TrueGammaCandidate0->IsLargestComponentPhoton())
          fHistoTruePi0PureGammaInvMassPtAlpha[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),  tempTruePi0CandWeight);
      }
    }
    if(isTrueEta && !fDoPi0Only){
      fHistoTrueEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), tempTruePi0CandWeight);
      if(previouslyNotFoundTrueMesons && !fDoLightOutput && fHistoTrueEtaInvMassPtAdditional[fiCut]) fHistoTrueEtaInvMassPtAdditional[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),tempTruePi0CandWeight);
      if(fDoJetAnalysis && !fDoLightOutput){
        if(fConvJetReader->GetTrueNJets()>0){
          fHistoTrueEtaJetMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),  tempTruePi0CandWeight);
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
                fHistoTrueEtaInJetMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), tempTruePi0CandWeight);
                fHistoMotherEtainJetPtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), tempTruePi0CandWeight);
                fHistoMotherEtainJetPtPhi[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Phi(), tempTruePi0CandWeight);
              }
            }
            fHistoTrueDoubleCountingEtaJet[fiCut]->Fill(counter);
          }
          fTrueVectorJetEta.clear();
          fTrueVectorJetPhi.clear();
        }
      }

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
        if (isTrueEta && !fDoPi0Only)fHistoTrueEtaCaloPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
      // both particles are electrons
      if (TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate1->IsLargestComponentElectron() ) {
        if (isTruePi0) fHistoTruePi0CaloElectronInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        if (isTrueEta && !fDoPi0Only) fHistoTrueEtaCaloElectronInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
      // both particles are converted electrons
      if ((TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate0->IsConversion()) && (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion()) ){
        if (isTruePi0 )fHistoTruePi0CaloConvertedPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        if (isTrueEta && !fDoPi0Only)fHistoTrueEtaCaloConvertedPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
      // 1 gamma is converted the other one is normals
      if ( (TrueGammaCandidate0->IsLargestComponentPhoton() && TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion()) ||
        (TrueGammaCandidate1->IsLargestComponentPhoton() && TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate0->IsConversion())
      ) {
        if (isTruePi0) fHistoTruePi0CaloMixedPhotonConvPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        if (isTrueEta && !fDoPi0Only) fHistoTrueEtaCaloMixedPhotonConvPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
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
        if (isTrueEta && !fDoPi0Only) fHistoTrueEtaCaloMergedClusterInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
      // at least one of the photon is merged and part conv
      if (TrueGammaCandidate1->IsMergedPartConv() || TrueGammaCandidate0->IsMergedPartConv()) {
        if (isTruePi0)fHistoTruePi0CaloMergedClusterPartConvInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        if (isTrueEta && !fDoPi0Only)fHistoTrueEtaCaloMergedClusterPartConvInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
    }

    if (fDoMesonQA == 2 && fIsMC < 2){
      // category 1: 2 real photons unmerged
      if (TrueGammaCandidate0->IsLargestComponentPhoton() && !TrueGammaCandidate0->IsMerged() && TrueGammaCandidate1->IsLargestComponentPhoton() && !TrueGammaCandidate1->IsMerged()) {
        if (isTruePi0) fHistoTruePi0Category1[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        if (isTrueEta && !fDoPi0Only) fHistoTrueEtaCategory1[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
      // category 2, 3: 1 real photons unmerged,  1 electron (category 2 merged, category 3 unmerged )
      // -> photon 0 is unconverted
      if ( (TrueGammaCandidate0->IsLargestComponentPhoton() && !TrueGammaCandidate0->IsMerged()) && (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion())) {
        if (isTruePi0){
          if (TrueGammaCandidate1->IsMergedPartConv())  fHistoTruePi0Category2[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
          if (!TrueGammaCandidate1->IsMergedPartConv())  fHistoTruePi0Category3[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        }
        if (isTrueEta && !fDoPi0Only){
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
        if (isTrueEta && !fDoPi0Only){
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
        if (isTrueEta && !fDoPi0Only){
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
        fHistoTruePi0PtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), tempTruePi0CandWeight);
        fHistoTruePi0PtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),TMath::Abs(Pi0Candidate->GetAlpha()), tempTruePi0CandWeight);
        fHistoTruePi0PtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle(), tempTruePi0CandWeight);
        }
      } else if (isTrueEta && !fDoPi0Only){
        if ( Pi0Candidate->M() > 0.45 && Pi0Candidate->M() < 0.65){
        fHistoTrueEtaPtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), tempTruePi0CandWeight);
        fHistoTrueEtaPtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),TMath::Abs(Pi0Candidate->GetAlpha()), tempTruePi0CandWeight);
        fHistoTrueEtaPtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle(), tempTruePi0CandWeight);
        }
      }
    }

    Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma0MotherLabel)), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
    if(!isPrimary){ // Secondary Meson
      Long_t secMotherLabel = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MotherLabel))->GetMother();
      Float_t weightedSec= 1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(secMotherLabel, 0x0, fInputEvent) && static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->GetPdgCode()==310){
        weightedSec= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(secMotherLabel, 0x0, fInputEvent)/2.; //invariant mass is additive thus the weight for the daughters has to be devide by two for the K0s at a certain pt
        //cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
      }
      if (isTruePi0) {
        if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTrueSecondaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
        if(fDoJetAnalysis){
          if(fConvJetReader->GetTrueNJets()>0){
            if(!fDoLightOutput) fHistoTrueSecondaryPi0InvJetMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
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
                    if(!fDoLightOutput) fHistoTrueSecondaryPi0InvinJetMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
                    else fHistoTrueSecondaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
                  }
                }
              }
            fTrueVectorJetEta.clear();
            fTrueVectorJetPhi.clear();
          }
        }
      }
      if (secMotherLabel >-1){
        if(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->GetPdgCode()==310 && isTruePi0 ){
          if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTrueSecondaryPi0FromK0sInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
          if (fDoMesonQA > 0 && fDoMesonQA < 3 && fIsMC < 2 )fHistoTrueK0sWithPi0DaughterMCPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->Pt(), tempTruePi0CandWeight);
          if(fDoJetAnalysis){
            if(fConvJetReader->GetTrueNJets()>0){
              if(!fDoLightOutput) fHistoTrueSecondaryPi0FromK0sJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
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
                      if(!fDoLightOutput) fHistoTrueSecondaryPi0FromK0sinJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
                      else fHistoTrueSecondaryPi0FromK0sInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
                    }
                  }
                }
              fTrueVectorJetEta.clear();
              fTrueVectorJetPhi.clear();
            }
          }
        }
        if(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->GetPdgCode()==130 && isTruePi0 ){
          if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTrueSecondaryPi0FromK0lInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
          if (fDoMesonQA > 0 && fDoMesonQA < 3 && fIsMC < 2 )fHistoTrueK0lWithPi0DaughterMCPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->Pt(), tempTruePi0CandWeight);
          if(fDoJetAnalysis){
            if(fConvJetReader->GetTrueNJets()>0){
              if(!fDoLightOutput) fHistoTrueSecondaryPi0FromK0lJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
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
                      if(!fDoLightOutput) fHistoTrueSecondaryPi0FromK0linJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
                      else fHistoTrueSecondaryPi0FromK0lInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
                    }
                  }
                }
              fTrueVectorJetEta.clear();
              fTrueVectorJetPhi.clear();
            }
          }
        }
        if(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->GetPdgCode()==221 && isTruePi0){
          fHistoTrueSecondaryPi0FromEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
          if (fDoMesonQA > 0 && fDoMesonQA < 3 && fIsMC < 2)fHistoTrueEtaWithPi0DaughterMCPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->Pt(), tempTruePi0CandWeight);
        }
        if(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->GetPdgCode()==3122 && isTruePi0){
          if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTrueSecondaryPi0FromLambdaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
          if (fDoMesonQA > 0 && fDoMesonQA < 3 && fIsMC < 2)fHistoTrueLambdaWithPi0DaughterMCPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->Pt(), tempTruePi0CandWeight);
          if(fDoJetAnalysis){
            if(fConvJetReader->GetTrueNJets()>0){
              if(!fDoLightOutput) fHistoTrueSecondaryPi0FromLambdaJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
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
                      if(!fDoLightOutput) fHistoTrueSecondaryPi0FromLambdainJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
                      else fHistoTrueSecondaryPi0FromLambdaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
                    }
                  }
                }
              fTrueVectorJetEta.clear();
              fTrueVectorJetPhi.clear();
            }
          }
        }
      }
    } else{ // Only primary pi0 for efficiency calculation
      Float_t weighted= 1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma1MotherLabel, 0x0, fInputEvent)){
        if (static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MotherLabel))->Pt()>0.005){
        weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(gamma1MotherLabel, 0x0, fInputEvent);
        //                      cout << "rec \t " <<gamma1MotherLabel << "\t" <<  weighted << endl;
        }
      }
      if (isTruePi0){
        if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTruePrimaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* tempTruePi0CandWeight);
        fHistoTruePrimaryPi0W0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), tempTruePi0CandWeight);
        fProfileTruePrimaryPi0WeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* tempTruePi0CandWeight);
        if (CheckVectorForDoubleCount(fVectorDoubleCountTruePi0s,gamma0MotherLabel)) fHistoDoubleCountTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), weighted*tempTruePi0CandWeight);
        if(fDoJetAnalysis){
          Int_t MCLabelGamma0 = TrueGammaCandidate0->GetCaloPhotonMCLabel(0);
          Int_t MCLabelGamma1 = TrueGammaCandidate1->GetCaloPhotonMCLabel(0);
          AliAODMCParticle* gamma0 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(MCLabelGamma0));
          AliAODMCParticle* gamma1 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(MCLabelGamma1));
          Double_t TruePtPi0 = TMath::Sqrt(pow((gamma0->Px() + gamma1->Px()),2)+pow((gamma0->Py() + gamma1->Py()),2));
          if(fConvJetReader->GetTrueNJets()>0){
            if(!fDoLightOutput) fHistoTruePrimaryPi0JetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* tempTruePi0CandWeight);
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
                      fHistoTruePrimaryPi0inJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* tempTruePi0CandWeight);
                      fHistoTruePrimaryPi0InJetInvMassTruePt[fiCut]->Fill(Pi0Candidate->M(),TruePtPi0, tempTruePi0CandWeight);
                      Double_t dotproduct = fTrueVectorJetPx.at(i)*Pi0Candidate->Px() + fTrueVectorJetPy.at(i)*Pi0Candidate->Py() + fTrueVectorJetPz.at(i)*Pi0Candidate->Pz();
                      Double_t magn = pow(fTrueVectorJetPx.at(i),2) + pow(fTrueVectorJetPy.at(i),2) + pow(fTrueVectorJetPz.at(i),2);
                      Double_t z = dotproduct/magn;
                      fHistoTruePi0JetFragmFunc[fiCut]->Fill(z,fTrueVectorJetPt.at(i), weighted* tempTruePi0CandWeight);
                      fHistoTruePi0JetFragmFuncZInvMass[fiCut]->Fill(Pi0Candidate->M(), z, weighted* tempTruePi0CandWeight);
                    }else fHistoTruePrimaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* tempTruePi0CandWeight);
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
      } else if (isTrueEta && !fDoPi0Only){
        if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTruePrimaryEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* tempTruePi0CandWeight);
        fHistoTruePrimaryEtaW0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), tempTruePi0CandWeight);
        fProfileTruePrimaryEtaWeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* tempTruePi0CandWeight);
        if (CheckVectorForDoubleCount(fVectorDoubleCountTrueEtas,gamma0MotherLabel)) fHistoDoubleCountTrueEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), weighted*tempTruePi0CandWeight);
        if(fDoJetAnalysis){
          if(fConvJetReader->GetTrueNJets()>0){
            if(!fDoLightOutput) fHistoTruePrimaryEtaJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* tempTruePi0CandWeight);
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
                  if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius() && !fDoPi0Only){
                    if(!fDoLightOutput){
                      fHistoTruePrimaryEtainJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* tempTruePi0CandWeight);
                      Double_t dotproduct = fTrueVectorJetPx.at(i)*Pi0Candidate->Px() + fTrueVectorJetPy.at(i)*Pi0Candidate->Py() + fTrueVectorJetPz.at(i)*Pi0Candidate->Pz();
                      Double_t magn = pow(fTrueVectorJetPx.at(i),2) + pow(fTrueVectorJetPy.at(i),2) + pow(fTrueVectorJetPz.at(i),2);
                      Double_t z = dotproduct/magn;
                      fHistoTrueEtaJetFragmFunc[fiCut]->Fill(z,fTrueVectorJetPt.at(i), weighted* tempTruePi0CandWeight);
                      fHistoTrueEtaJetFragmFuncZInvMass[fiCut]->Fill(Pi0Candidate->M(), z, weighted* tempTruePi0CandWeight);
                    }else fHistoTruePrimaryEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* tempTruePi0CandWeight);
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
      }
      if (fDoMesonQA > 0 && fDoMesonQA < 3 && fIsMC<2){
        if(isTruePi0){ // Only primary pi0 for resolution
          fHistoTruePrimaryPi0MCPtResolPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MotherLabel))->Pt(),
                              (Pi0Candidate->Pt()-static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MotherLabel))->Pt())/static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MotherLabel))->Pt(),weighted* tempTruePi0CandWeight);

        }
        if (isTrueEta && !fDoPi0Only){ // Only primary eta for resolution
          fHistoTruePrimaryEtaMCPtResolPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MotherLabel))->Pt(),
                              (Pi0Candidate->Pt()-static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MotherLabel))->Pt())/static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MotherLabel))->Pt(),weighted* tempTruePi0CandWeight);
        }
      }
    }
  } else if(!isTruePi0 && !isTrueEta) { // Background
    if (fDoMesonQA > 0 && fDoMesonQA < 3){
      if(gamma0MotherLabel>-1 && gamma1MotherLabel>-1){ // Both Tracks are Photons and have a mother but not Pi0 or Eta
        fHistoTrueBckGGInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), tempTruePi0CandWeight);
        if(gamma0MotherLabel < fAODMCTrackArray->GetEntries() && gamma1MotherLabel < fAODMCTrackArray->GetEntries()){
          AliAODMCParticle *trackgamma0 = NULL;
          AliAODMCParticle *trackgamma1 = NULL;
          trackgamma0 = (AliAODMCParticle*)fAODMCTrackArray->At(gamma0MotherLabel);
          trackgamma1 = (AliAODMCParticle*)fAODMCTrackArray->At(gamma1MotherLabel);
          if( trackgamma0 != NULL && trackgamma1 != NULL){

            if( ((((AliAODMCParticle*)fAODMCTrackArray->At(gamma0MotherLabel))->GetPdgCode() == 111
                || ((AliAODMCParticle*)fAODMCTrackArray->At(gamma0MotherLabel))->GetPdgCode() == 221)
                && (TrueGammaCandidate0->IsMerged() || TrueGammaCandidate0->IsMergedPartConv()))
                ||
                ((((AliAODMCParticle*)fAODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 111
                || ((AliAODMCParticle*)fAODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 221)
                && (TrueGammaCandidate1->IsMerged() || TrueGammaCandidate1->IsMergedPartConv()))
            ){
              fHistoTrueBckFullMesonContainedInOneClusterInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), tempTruePi0CandWeight);
            }else if( (TrueGammaCandidate0->E()/Pi0Candidate->E() > 0.7) || (TrueGammaCandidate1->E()/Pi0Candidate->E() > 0.7) ){
              fHistoTrueBckAsymEClustersInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), tempTruePi0CandWeight);
            }
          }
        }
      } else { // No photon or without mother
        fHistoTrueBckContInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), tempTruePi0CandWeight);
      }
    }
  }
}
//________________________________________________________________________
void AliAnalysisTaskGammaCalo::ProcessAODSphericityParticles()
{
  if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (fAODMCTrackArray == NULL) return;
  for(Long_t i = 0; i < fAODMCTrackArray->GetEntriesFast(); i++) {
      AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(i));
      if(!particle) continue;
      if(!particle->IsPhysicalPrimary()) continue;
      if(!particle->Pt()) continue;
      if(TMath::Abs(particle->Eta())>0.8) continue;
      if(particle->GetPdgCode() == 211 || particle->GetPdgCode() == -211){ //pi+-
          fHistoPionSpectrum[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
      }else if(particle->GetPdgCode() == 2212 || particle->GetPdgCode() == -2212){ //proton anti-proton
          fHistoProtonSpectrum[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
      }else if(particle->GetPdgCode() == 321 || particle->GetPdgCode() == -321){ // K+-
          fHistoKaonSpectrum[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
      }else if(particle->GetPdgCode() == 111){ //neutral pion
          fHistoNPionSpectrum[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
      }else if(particle->GetPdgCode() == 221){ // eta
          fHistoEtaSpectrum[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
      }else if(particle->GetPdgCode() == 421){ // D0 meson
          fHistoDMesonSpectrum[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
      }

  }
  return;
}
//________________________________________________________________________
void AliAnalysisTaskGammaCalo::CalculateBackground(){

  Int_t zbin= 0;
  if( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoSectorMixing() ) {
    zbin = fBGHandler[fiCut]->GetZBinIndex(fV0Reader->GetPtMaxSector());
  } else {
    zbin = fBGHandler[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
  }
  Int_t mbin = 0;

  Double_t tempBGCandidateWeight       = fWeightJetJetMC;

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
        for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
          AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));
          AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
          backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());

          // Set the BG candidate jetjet weight to 1 in case both photons orignated from the minimum bias header
          if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4){
            if( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(previousGoodV0.GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2 &&
                ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(currentEventGoodV0.GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2)
              tempBGCandidateWeight = 1;
          }

          if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
            ->MesonIsSelected(backgroundCandidate,kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()),currentEventGoodV0.GetLeadingCellID(),previousGoodV0.GetLeadingCellID(),currentEventGoodV0.GetIsCaloPhoton(), previousGoodV0.GetIsCaloPhoton() )){
            fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(), tempBGCandidateWeight);
            if(fDoTHnSparse){
              Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
              fSparseMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
            }
            if((!fDoLightOutput || fDoPi0Only || fDoECalibOutput) && TMath::Abs(backgroundCandidate->GetAlpha())<0.1){
              fHistoMotherBackInvMassECalib[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->E(),tempBGCandidateWeight);
            }

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

    Double_t currentPtMax   = 0;  Double_t previousPtMax  = 0;
    Double_t currentAvePt   = 0;  Double_t previousAvePt  = 0;
    Double_t currentAveEta  = 0;  Double_t previousAveEta = 0;
    Double_t currentAvePhi  = 0;  Double_t previousAvePhi = 0;
    Bool_t acceptedPtMax    = kFALSE;

    for(Int_t nEventsInBG=0;nEventsInBG <fBGHandler[fiCut]->GetNBGEvents();nEventsInBG++){
      AliGammaConversionAODVector *previousEventV0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,nEventsInBG);
      if(previousEventV0s){
        acceptedPtMax = kFALSE;
        currentPtMax = 0; previousPtMax = 0;
        currentAvePt = 0; previousAvePt = 0; currentAveEta = 0; previousAveEta = 0; currentAvePhi = 0; previousAvePhi = 0;
        for(Int_t iCurrent=0;iCurrent<fClusterCandidates->GetEntries();iCurrent++){
          AliAODConversionPhoton *currentV0 = (AliAODConversionPhoton*)(fClusterCandidates->At(iCurrent));
          currentAvePt += currentV0->GetPhotonPt();
          currentAveEta += currentV0->GetPhotonPt()*currentV0->GetPhotonEta();
          currentAvePhi += currentV0->GetPhotonPt()*currentV0->GetPhotonPhi();
          if(currentV0->GetPhotonPt() > currentPtMax){ currentPtMax = currentV0->GetPhotonPt(); }
        }
        currentAveEta /= currentAvePt;
        currentAvePhi /= currentAvePt;
        currentAvePt /= fClusterCandidates->GetEntries();
        for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
            AliAODConversionPhoton *previousV0 = (AliAODConversionPhoton*)(previousEventV0s->at(iPrevious));
            previousAvePt += previousV0->GetPhotonPt();
            previousAveEta += previousV0->GetPhotonPt()*previousV0->GetPhotonEta();
            previousAvePhi += previousV0->GetPhotonPt()*previousV0->GetPhotonPhi();
            if(previousV0->GetPhotonPt() > previousPtMax){ previousPtMax = previousV0->GetPhotonPt(); }
        }
        previousAveEta /= previousAvePt;
        previousAvePhi /= previousAvePt;
        previousAvePt /= previousEventV0s->size();
        if(currentPtMax > 0. && previousPtMax > 0.){
         //if(TMath::Sqrt(pow((currentEta-previousEta),2)+pow((currentPhi-previousPhi),2)) < 0.2) acceptedPtMax = kTRUE;
         if(TMath::Abs(previousAveEta-currentAveEta)<0.4 && TMath::Abs(previousAvePhi-currentAvePhi)<0.6 && (previousAvePt/currentAvePt)<4. && (previousAvePt/currentAvePt)>0.25) acceptedPtMax = kTRUE;
        }
        if(acceptedPtMax){
          for(Int_t iCurrent=0;iCurrent<fClusterCandidates->GetEntries();iCurrent++){
            AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fClusterCandidates->At(iCurrent));
            for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
              AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));
              AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
              backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());

              // Set the BG candidate jetjet weight to 1 in case both photons orignated from the minimum bias header
              if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4){
                if( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(previousGoodV0.GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2 &&
                    ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(currentEventGoodV0.GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2)
                  tempBGCandidateWeight = 1;
              }

              if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
                ->MesonIsSelected(backgroundCandidate,kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()),currentEventGoodV0.GetLeadingCellID(),previousGoodV0.GetLeadingCellID(), currentEventGoodV0.GetIsCaloPhoton(), previousGoodV0.GetIsCaloPhoton())){
                fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(), tempBGCandidateWeight);
              }
              delete backgroundCandidate;
              backgroundCandidate = 0x0;
            }
          }
        }
      }
    }
  } else if( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoSectorMixing() ) {
    if(fClusterCandidates->GetEntries()>0){
      for(Int_t nEventsInBG=0;nEventsInBG <fBGHandler[fiCut]->GetNBGEvents();nEventsInBG++){
        AliGammaConversionAODVector *previousEventV0s = NULL;
        previousEventV0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,nEventsInBG);
        if(previousEventV0s && previousEventV0s->size()>0){
              for(Int_t iCurrent=0;iCurrent<fClusterCandidates->GetEntries();iCurrent++){
                AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fClusterCandidates->At(iCurrent));
                for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
                  AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));
                  AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
                  backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());

                  if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate,kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),currentEventGoodV0.GetLeadingCellID(),previousGoodV0.GetLeadingCellID(), currentEventGoodV0.GetIsCaloPhoton(), previousGoodV0.GetIsCaloPhoton()))){
                    // Set the BG candidate jetjet weight to 1 in case both photons orignated from the minimum bias header
                    if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4){
                      if( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(previousGoodV0.GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2 &&
                          ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(currentEventGoodV0.GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2)
                        tempBGCandidateWeight = 1;
                    }
                    if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(), tempBGCandidateWeight);
                    if(fDoJetAnalysis){
                      if(fConvJetReader->GetNJets() > 0){
                        if(!fDoLightOutput) fHistoMotherBackJetInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(), tempBGCandidateWeight);
                        else fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(), tempBGCandidateWeight);
                      }
                    }
                    if(fDoTHnSparse){
                      Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
                      fSparseMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
                    }
                    if((!fDoLightOutput || fDoPi0Only || fDoECalibOutput) && TMath::Abs(backgroundCandidate->GetAlpha())<0.1){
                      fHistoMotherBackInvMassECalib[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->E(),tempBGCandidateWeight);
                    }

                    if (fDoMesonQA == 2){
                      fHistoMotherPtOpenAngleBck[fiCut]->Fill(backgroundCandidate->Pt(),backgroundCandidate->GetOpeningAngle(), tempBGCandidateWeight);
                    }
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
  }else if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoJetMixing()){
    if(fDoJetAnalysis){
      for(Int_t iCurrent=0;iCurrent<fClusterCandidates->GetEntries();iCurrent++){
        AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fClusterCandidates->At(iCurrent));
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
                  AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
                  backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());

                  if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),currentEventGoodV0.GetLeadingCellID(),previousGoodV0.GetLeadingCellID(), currentEventGoodV0.GetIsCaloPhoton(), previousGoodV0.GetIsCaloPhoton()  ))){
                    if(!fDoLightOutput) fHistoMotherBackJetInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(), tempBGCandidateWeight);
                    else fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(), tempBGCandidateWeight);
                  }
                  delete backgroundCandidate;
                  backgroundCandidate = 0x0;
                }
              }
            }
          }
          fVectorJetEta.clear();
          fVectorJetPhi.clear();
        }
      }
    }
  } else {
    for(Int_t nEventsInBG=0;nEventsInBG <fBGHandler[fiCut]->GetNBGEvents();nEventsInBG++){
      AliGammaConversionAODVector *previousEventV0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,nEventsInBG);
      if(previousEventV0s){
        for(Int_t iCurrent=0;iCurrent<fClusterCandidates->GetEntries();iCurrent++){
          AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fClusterCandidates->At(iCurrent));
          for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){

            AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));
            std::unique_ptr<AliAODConversionMother> backgroundCandidate (new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0));
            backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());

            if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate.get(),kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),currentEventGoodV0.GetLeadingCellID(),previousGoodV0.GetLeadingCellID(),currentEventGoodV0.GetIsCaloPhoton(),  previousGoodV0.GetIsCaloPhoton() ))){
              // Set the BG candidate jetjet weight to 1 in case both photons orignated from the minimum bias header
              if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4){
                if( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(previousGoodV0.GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2 &&
                    ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(currentEventGoodV0.GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2)
                  tempBGCandidateWeight = 1;
              }
              if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(), tempBGCandidateWeight);
              if(fDoJetAnalysis){
                if(fConvJetReader->GetNJets() > 0){
                  if(!fDoLightOutput) fHistoMotherBackJetInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(), tempBGCandidateWeight);
                  else fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(), tempBGCandidateWeight);
                }
              }
              if(fDoTHnSparse){
                Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
                fSparseMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
              }
              if((!fDoLightOutput || fDoPi0Only || fDoECalibOutput) && TMath::Abs(backgroundCandidate->GetAlpha())<0.1){
                fHistoMotherBackInvMassECalib[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->E(),tempBGCandidateWeight);
              }

              if (fDoMesonQA == 2){
                  fHistoMotherPtOpenAngleBck[fiCut]->Fill(backgroundCandidate->Pt(),backgroundCandidate->GetOpeningAngle(), tempBGCandidateWeight);
              }
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
          }
        }
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskGammaCalo::CalculateBackgroundSwapp(){

  if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoGammaSwappForBg()) {

    Double_t rotationAngle = TMath::Pi()/2.0; //0.78539816339; // rotaion angle 90°

    TLorentzVector lvRotationPhoton1;   // photon candidates which get rotated
    TLorentzVector lvRotationPhoton2;   // photon candidates which get rotated
    TVector3 lvRotationPion;            // reconstructed mother particle from the two photons
    // Needed for TGenPhaseSpace
    TVector3 tvEtaPhigamma1, tvEtaPhigamma2, tvEtaPhigamma1Decay, tvEtaPhigamma2Decay, tvNormBeforeDecay, tvNormAfterDecay;
    Float_t asymBeforeDecay = 0.;
    Float_t asymAfterDecay = 0.;
    Double_t massGamma[2] = {0,0};

    Int_t cellIDRotatedPhoton1 = -1; // cell ID of the cluster after rotation
    Int_t cellIDRotatedPhoton2 = -1; // cell ID of the cluster after rotation

    std::vector<std::array<Double_t, 2>> vSwappingInvMassPT;
    std::vector<std::array<Double_t, 2>> vSwappingInvMassPTAlphaCut;
    vSwappingInvMassPT.clear();
    vSwappingInvMassPTAlphaCut.clear();
    vSwappingInvMassPT.resize(0);
    vSwappingInvMassPTAlphaCut.resize(0);
    Double_t tempMultWeightSwapping = 1; // weight taking multiplicity of event into account

    // curcial requierment is that the event has at least 3 cluster candidates
    if(fClusterCandidates->GetEntries() > 2 ){

      for(Int_t iCurrent1=0;iCurrent1<fClusterCandidates->GetEntries();iCurrent1++){
        AliAODConversionPhoton* currentEventGoodV0Temp1 = (AliAODConversionPhoton*)(fClusterCandidates->At(iCurrent1));

        for(Int_t iCurrent2=iCurrent1+1;iCurrent2<fClusterCandidates->GetEntries();iCurrent2++){
          AliAODConversionPhoton* currentEventGoodV0Temp2 = (AliAODConversionPhoton*)(fClusterCandidates->At(iCurrent2));

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
              if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 0) rotationAngle = TMath::Pi()/2.0; // rotate by 90 degree
              else if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 1){  // rotate by random angle between
                 Double_t temp = (fRandom.Rndm() < 0.5) ? 0 : TMath::Pi();
                 rotationAngle = temp + TMath::Pi()/3.0 + fRandom.Rndm()*TMath::Pi()/3.0;
              }
              lvRotationPhoton1.Rotate(rotationAngle, lvRotationPion);
              lvRotationPhoton2.Rotate(rotationAngle, lvRotationPion);
            } else if (((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() >=10){ // generate new decay with TGenPhaseSpace
              if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 11){
                tvEtaPhigamma1 = lvRotationPhoton1.Vect();
                tvEtaPhigamma2 = lvRotationPhoton2.Vect();
                tvNormBeforeDecay = tvEtaPhigamma1.Cross(tvEtaPhigamma2);
                asymBeforeDecay = fabs((lvRotationPhoton1.E()-lvRotationPhoton2.E())/(lvRotationPhoton1.E()+lvRotationPhoton2.E()));
              }

              TLorentzVector lvRotationMother = lvRotationPhoton1 + lvRotationPhoton2;
              fGenPhaseSpace.SetDecay(lvRotationMother, 2, massGamma);
              fGenPhaseSpace.Generate();
              lvRotationPhoton1 = *fGenPhaseSpace.GetDecay(0);
              lvRotationPhoton2 = *fGenPhaseSpace.GetDecay(1);

              if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 11){
                tvEtaPhigamma1Decay = lvRotationPhoton1.Vect();
                tvEtaPhigamma2Decay = lvRotationPhoton2.Vect();
                tvNormAfterDecay = tvEtaPhigamma1Decay.Cross(tvEtaPhigamma2Decay);  // norm vector to decay plane
                asymAfterDecay = fabs((lvRotationPhoton1.E()-lvRotationPhoton2.E())/(lvRotationPhoton1.E()+lvRotationPhoton2.E()));
                // check if decay is nearly the same as original decay: if yes continue with next decay
                if((tvNormAfterDecay.Angle(tvNormBeforeDecay) < 20*TMath::Pi()/180. || tvNormAfterDecay.Angle(tvNormBeforeDecay) > 340*TMath::Pi()/180.) && ( fabs(asymBeforeDecay - asymAfterDecay) < 0.05 )   ) continue;
              }

            }


            cellIDRotatedPhoton1 = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetCaloCellIdFromEtaPhi(lvRotationPhoton1.Eta(), static_cast<double>((lvRotationPhoton1.Phi()<0) ? lvRotationPhoton1.Phi() + TMath::Pi()*2. : lvRotationPhoton1.Phi()));
            cellIDRotatedPhoton2 = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetCaloCellIdFromEtaPhi(lvRotationPhoton2.Eta(), static_cast<double>((lvRotationPhoton2.Phi()<0) ? lvRotationPhoton2.Phi() + TMath::Pi()*2. : lvRotationPhoton2.Phi()));

            if(!fDoLightOutput){
              if(!(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedPhoton1, lvRotationPhoton1.Phi(), fInputEvent))){
                ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FillEtaPhiMapForClusterInBg(lvRotationPhoton1.Eta(), static_cast<double>((lvRotationPhoton1.Phi()<0) ? lvRotationPhoton1.Phi() + TMath::Pi()*2. : lvRotationPhoton1.Phi()), 1);
              }
              if(!(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedPhoton2, lvRotationPhoton2.Phi(), fInputEvent))){
                ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FillEtaPhiMapForClusterInBg(lvRotationPhoton2.Eta(), static_cast<double>((lvRotationPhoton2.Phi()<0) ? lvRotationPhoton2.Phi() + TMath::Pi()*2. : lvRotationPhoton2.Phi()), 1);
              }
            }

            std::unique_ptr<AliAODConversionPhoton> currentEventGoodV0Rotation1 (new AliAODConversionPhoton(&lvRotationPhoton1));
            std::unique_ptr<AliAODConversionPhoton> currentEventGoodV0Rotation2 (new AliAODConversionPhoton(&lvRotationPhoton2));

            for(auto const& kCurrentClusterCandidates  : *fClusterCandidates){
              if(currentEventGoodV0Temp1 == ((AliAODConversionPhoton*) kCurrentClusterCandidates) || currentEventGoodV0Temp2 == ((AliAODConversionPhoton*) kCurrentClusterCandidates)){ continue;}

              std::unique_ptr<AliAODConversionMother> backgroundCandidate1(new AliAODConversionMother(currentEventGoodV0Rotation1.get(), ((AliAODConversionPhoton*) kCurrentClusterCandidates)));
              std::unique_ptr<AliAODConversionMother> backgroundCandidate2(new AliAODConversionMother(currentEventGoodV0Rotation2.get(), ((AliAODConversionPhoton*) kCurrentClusterCandidates)));

              if(!(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedPhoton1, lvRotationPhoton1.Phi(), fInputEvent)) && lvRotationPhoton1.E() > ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetMinClusterEnergy())
              {
                if(((AliConversionMesonCuts*) fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate1.get(),kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), cellIDRotatedPhoton1, ((AliAODConversionPhoton*) kCurrentClusterCandidates)->GetLeadingCellID()))
                {
                  vSwappingInvMassPT.push_back({backgroundCandidate1->M(),backgroundCandidate1->Pt()});
                  if((!fDoLightOutput || fDoPi0Only) && TMath::Abs(backgroundCandidate1->GetAlpha())<0.1){
                    vSwappingInvMassPTAlphaCut.push_back({backgroundCandidate1->M(),backgroundCandidate1->Pt()});
                  }
                }
              }
              if(!(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedPhoton2, lvRotationPhoton2.Phi(), fInputEvent)) && lvRotationPhoton2.E() > ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetMinClusterEnergy())
              {
                if(((AliConversionMesonCuts*) fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate2.get(),kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), cellIDRotatedPhoton2, ((AliAODConversionPhoton*) kCurrentClusterCandidates)->GetLeadingCellID()))
                {
                  vSwappingInvMassPT.push_back({backgroundCandidate2->M(),backgroundCandidate2->Pt()});
                  if((!fDoLightOutput || fDoPi0Only) && TMath::Abs(backgroundCandidate2->GetAlpha())<0.1){
                    vSwappingInvMassPTAlphaCut.push_back({backgroundCandidate2->M(),backgroundCandidate2->Pt()});
                  }
                }
              }
            }
          }
        }
      }
      // Fill the histograms
      if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoWeightingInSwappBg() && vSwappingInvMassPT.size() > 0){
        tempMultWeightSwapping = (0.5*(fClusterCandidates->GetEntries()*fClusterCandidates->GetEntries() - fClusterCandidates->GetEntries()))/(vSwappingInvMassPT.size());
      }
      for(Int_t i = 0; i < (Int_t)vSwappingInvMassPT.size(); i++){
        fHistoMotherBackInvMassPt[fiCut]->Fill(vSwappingInvMassPT.at(i)[0], vSwappingInvMassPT.at(i)[1], tempMultWeightSwapping*fWeightJetJetMC);
      }
      if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoWeightingInSwappBg() && vSwappingInvMassPTAlphaCut.size() > 0){
        tempMultWeightSwapping = (0.5*(fClusterCandidates->GetEntries()*fClusterCandidates->GetEntries() - fClusterCandidates->GetEntries()))/(vSwappingInvMassPTAlphaCut.size());
      }
      for(Int_t i = 0; i < (Int_t)vSwappingInvMassPTAlphaCut.size(); i++){
        fHistoMotherBackInvMassECalib[fiCut]->Fill(vSwappingInvMassPTAlphaCut.at(i)[0], vSwappingInvMassPTAlphaCut.at(i)[1],tempMultWeightSwapping*fWeightJetJetMC);
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
  if(fDoJetAnalysis && fConvJetReader->GetNJets() == 0) return;
  if(fClusterCandidates->GetEntries() >1 ){
    if( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoSectorMixing() ){
      fBGHandler[fiCut]->AddEvent(fClusterCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fV0Reader->GetPtMaxSector(),fClusterCandidates->GetEntries(),fEventPlaneAngle);
    } else if( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoSectorJetMixing() ){
      fBGHandler[fiCut]->AddEvent(fClusterCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(), fJetSector ,fClusterCandidates->GetEntries(),fEventPlaneAngle);
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
    } else if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
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
  Double_t tempParticleWeight       = fWeightJetJetMC;
  if(!fDoLightOutput){
    if(fIsFromDesiredHeader){
      // Set the jetjet weight to 1 in case the particle orignated from the minimum bias header
      if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4){
        if( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(TruePhotonCandidate->GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2) tempParticleWeight = 1;
      }
      if(TMath::Abs(pdgCode) == 11)        fHistoClusPhotonBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),0., tempParticleWeight);
      else if( TMath::Abs(pdgCode)==211)   fHistoClusPhotonBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),1., tempParticleWeight);
      else if( TMath::Abs(pdgCode)==2212)  fHistoClusPhotonBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),2., tempParticleWeight);
      else if( TMath::Abs(pdgCode)==321)   fHistoClusPhotonBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3., tempParticleWeight);
      else if( TMath::Abs(pdgCode)==2112)  fHistoClusPhotonBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),4., tempParticleWeight);
      else if( TMath::Abs(pdgCode)==310)   fHistoClusPhotonBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),5., tempParticleWeight);
      else if( TMath::Abs(pdgCode)==3122)  fHistoClusPhotonBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),6., tempParticleWeight);
      else if( TMath::Abs(pdgCode)==13)    fHistoClusPhotonBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),7., tempParticleWeight);
      else if( TMath::Abs(pdgCode)==130)   fHistoClusPhotonBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),8., tempParticleWeight);
      else                          fHistoClusPhotonBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),9., tempParticleWeight);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskGammaCalo::FillPhotonPlusConversionBackgroundHist(AliAODConversionPhoton *TruePhotonCandidate, Int_t pdgCode)
{
  // Bck = 0 e+-, 1 pi+-, 2 p+-, 3 K+-, 4 n, 5 K0s, 6 Lambda, 7 mu+-, 8 K0l, 9 rest
  Double_t tempParticleWeight       = fWeightJetJetMC;
  if(!fDoLightOutput && fHistoClusPhotonPlusConvBGPt[fiCut]){
    if(fIsFromDesiredHeader){
      // Set the jetjet weight to 1 in case the particle orignated from the minimum bias header
      if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4){
        if( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(TruePhotonCandidate->GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2) tempParticleWeight = 1;
      }
      if(TMath::Abs(pdgCode) == 11)        fHistoClusPhotonPlusConvBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),0., tempParticleWeight);
      else if( TMath::Abs(pdgCode)==211)   fHistoClusPhotonPlusConvBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),1., tempParticleWeight);
      else if( TMath::Abs(pdgCode)==2212)  fHistoClusPhotonPlusConvBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),2., tempParticleWeight);
      else if( TMath::Abs(pdgCode)==321)   fHistoClusPhotonPlusConvBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3., tempParticleWeight);
      else if( TMath::Abs(pdgCode)==2112)  fHistoClusPhotonPlusConvBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),4., tempParticleWeight);
      else if( TMath::Abs(pdgCode)==310)   fHistoClusPhotonPlusConvBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),5., tempParticleWeight);
      else if( TMath::Abs(pdgCode)==3122)  fHistoClusPhotonPlusConvBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),6., tempParticleWeight);
      else if( TMath::Abs(pdgCode)==13)    fHistoClusPhotonPlusConvBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),7., tempParticleWeight);
      else if( TMath::Abs(pdgCode)==130)   fHistoClusPhotonPlusConvBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),8., tempParticleWeight);
      else                          fHistoClusPhotonPlusConvBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),9., tempParticleWeight);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskGammaCalo::FillPhotonBackgroundM02Hist(AliAODConversionPhoton *TruePhotonCandidate, Double_t clusterM02, Int_t pdgCode)
{
  // Bck = e+-, pi+-, n, K+-, K0l, rest
  Double_t tempParticleWeight       = fWeightJetJetMC;
  if(!fDoLightOutput && fDoClusterQA > 1){
    if(fIsFromDesiredHeader){
      // Set the jetjet weight to 1 in case the particle orignated from the minimum bias header
      if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4){
        if( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(TruePhotonCandidate->GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2) tempParticleWeight = 1;
      }
      if(TMath::Abs(pdgCode) == 11 )          fHistoClustPhotonElectronBGPtM02[fiCut]->Fill(TruePhotonCandidate->Pt(),clusterM02, tempParticleWeight);
      else if( TMath::Abs(pdgCode)==211 )     fHistoClustPhotonPionBGPtM02[fiCut]->Fill(    TruePhotonCandidate->Pt(),clusterM02, tempParticleWeight);
      else if( TMath::Abs(pdgCode)==2212 )    fHistoClustPhotonNeutronBGPtM02[fiCut]->Fill( TruePhotonCandidate->Pt(),clusterM02, tempParticleWeight);
      else if( TMath::Abs(pdgCode)==321 )     fHistoClustPhotonKaonBGPtM02[fiCut]->Fill(    TruePhotonCandidate->Pt(),clusterM02, tempParticleWeight);
      else if( TMath::Abs(pdgCode)==130 )     fHistoClustPhotonK0lBGPtM02[fiCut]->Fill(     TruePhotonCandidate->Pt(),clusterM02, tempParticleWeight);
      else                                    fHistoClustPhotonRestBGPtM02[fiCut]->Fill(    TruePhotonCandidate->Pt(),clusterM02, tempParticleWeight);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskGammaCalo::FillPhotonPlusConversionBackgroundM02Hist(AliAODConversionPhoton *TruePhotonCandidate, Double_t clusterM02, Int_t pdgCode)
{
  // Bck = e+-, pi+-, n, K+-, K0l, rest
  Double_t tempParticleWeight       = fWeightJetJetMC;
  if(!fDoLightOutput  && fDoClusterQA > 1){
    if(fIsFromDesiredHeader){
      // Set the jetjet weight to 1 in case the particle orignated from the minimum bias header
      if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4){
        if( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(TruePhotonCandidate->GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2) tempParticleWeight = 1;
      }
      if(TMath::Abs(pdgCode) == 11 )          fHistoClustPhotonPlusConvElectronBGPtM02[fiCut]->Fill(TruePhotonCandidate->Pt(),clusterM02, tempParticleWeight);
      else if( TMath::Abs(pdgCode)==211 )     fHistoClustPhotonPlusConvPionBGPtM02[fiCut]->Fill(    TruePhotonCandidate->Pt(),clusterM02, tempParticleWeight);
      else if( TMath::Abs(pdgCode)==2212)     fHistoClustPhotonPlusConvNeutronBGPtM02[fiCut]->Fill( TruePhotonCandidate->Pt(),clusterM02, tempParticleWeight);
      else if( TMath::Abs(pdgCode)==321 )     fHistoClustPhotonPlusConvKaonBGPtM02[fiCut]->Fill(    TruePhotonCandidate->Pt(),clusterM02, tempParticleWeight);
      else if( TMath::Abs(pdgCode)==130 )     fHistoClustPhotonPlusConvK0lBGPtM02[fiCut]->Fill(     TruePhotonCandidate->Pt(),clusterM02, tempParticleWeight);
      else                                    fHistoClustPhotonPlusConvRestBGPtM02[fiCut]->Fill(    TruePhotonCandidate->Pt(),clusterM02, tempParticleWeight);
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
Int_t AliAnalysisTaskGammaCalo::WhichDDL(Int_t module, Int_t cellx)
{
  const Int_t Nmod=5;//totally, 5 PHOS modules are designed.
  Int_t ddl = -1;

  if(cellx<1 || 64<cellx) return -1;

  if(module<1 || 4<module){
    return -1;
  }
  else{
    ddl = (Nmod-module) * 4 + (cellx-1)/16;//convert offline module numbering to online.
    return ddl;
  }
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

  Int_t nclus = 0;
  TClonesArray * arrClustersDebug = NULL;
  if(!fCorrTaskSetting.CompareTo("")){
    nclus = fInputEvent->GetNumberOfCaloClusters();
  } else {
    arrClustersDebug = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));
    if(!arrClustersDebug)
      AliFatal(Form("%sClustersBranch was not found in AliAnalysisTaskGammaCalo! Check the correction framework settings!",fCorrTaskSetting.Data()));
    nclus = arrClustersDebug->GetEntries();
  }

  AliVCaloCells* cells = fInputEvent->GetEMCALCells();
  fOutputLocalDebug << "--event--" << endl;
  fOutputLocalDebug << "nclusters " << nclus << endl;

  // vertex
  Double_t vertex[3] = {0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  for(Long_t i = 0; i < nclus; i++){
    AliVCluster* clus = NULL;
    if(fInputEvent->IsA()==AliESDEvent::Class()){
      if(arrClustersDebug)
        clus = new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersDebug->At(i));
      else
        clus = new AliESDCaloCluster(*(AliESDCaloCluster*)fInputEvent->GetCaloCluster(i));
    } else if(fInputEvent->IsA()==AliAODEvent::Class()){
      if(arrClustersDebug)
        clus = new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersDebug->At(i));
      else
        clus = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(i));
    }

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
// Function that searches iteratively for
// contributions of a mcparticle (and its daughters)
// in a cluster and returns the cluster id
Int_t AliAnalysisTaskGammaCalo::CheckClustersForMCContribution(Int_t mclabel, Bool_t leading)
{
    Int_t clusID = -1; // found contribution in ID cluster
    Int_t nclus                       = 0;
    TClonesArray * arrClustersProcess = NULL;
    fNCurrentClusterBasic             = 0;
    TParticle* particle = (TParticle *)fMCEvent->Particle(mclabel);
    if(!fCorrTaskSetting.CompareTo("")){
      nclus = fInputEvent->GetNumberOfCaloClusters();
    } else {
      arrClustersProcess = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));
      if(!arrClustersProcess)
      AliFatal(Form("%sClustersBranch was not found in AliAnalysisTaskGammaCalo! Check the correction framework settings!",fCorrTaskSetting.Data()));
      nclus = arrClustersProcess->GetEntries();
    }

    for(Long_t i = 0; i < nclus; i++){
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
      if(!clus) continue;

      if(!((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC, 1.,i)){
         delete clus;
         continue;
      }

      Int_t *mclabelsCluster = clus->GetLabels();
      if (clus->GetNLabels() > 0)
      {
          if(!leading){
            for (Int_t k = 0; k < (Int_t)clus->GetNLabels(); k++)
            {
              if ((mclabelsCluster[k] == mclabel) ) clusID = i;
            }
          } else{
            if ((mclabelsCluster[0] == mclabel)) clusID = i;
          }
      }

      if(clusID != -1){ // do a check if its a conversion electron
          Bool_t motherIsPhoton = kFALSE;
          Int_t absPdg = TMath::Abs(particle->GetPdgCode());
          if(absPdg == 11){
              Int_t mothLabel = particle->GetMother(0);
              if(mothLabel!= -1){
                 TParticle* mother = (TParticle *)fMCEvent->Particle(mothLabel);
                 if(mother->GetPdgCode() == 22){
                        motherIsPhoton = kTRUE;
                 }
                 if(motherIsPhoton){
                   // dont use conv electrons that have smaller than 0.5 of mother photon
                   if(particle->Energy()/mother->Energy() < 0.5) clusID = -1;
                 }
              }
          }

      }
      delete clus;
    }
    // If i did not find a contribution for this mc label, check iteratively for daughters
    // and break as soon as i found a contribution of a daughter in a cluster
    if(clusID == -1 && (particle->GetNDaughters()>0)){
       for (Int_t daughter = particle->GetFirstDaughter(); daughter <= particle->GetLastDaughter(); daughter++)
       {
         clusID = CheckClustersForMCContribution(daughter, leading);
         if(clusID!=-1){
           break; // break when you found a contribution
         }
       }

    }
    return clusID;
}
Bool_t AliAnalysisTaskGammaCalo::CheckSpecificClusterForMCContribution(Int_t mclabel, Int_t cluslabel)
{
    TClonesArray * arrClustersProcess = NULL;
    fNCurrentClusterBasic             = 0;
    if(fCorrTaskSetting.CompareTo("")){
      arrClustersProcess = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));
      if(!arrClustersProcess)
      AliFatal(Form("%sClustersBranch was not found in AliAnalysisTaskGammaCalo! Check the correction framework settings!",fCorrTaskSetting.Data()));
    }

      AliVCluster* clus = NULL;
      if(fInputEvent->IsA()==AliESDEvent::Class()){
        if(arrClustersProcess)
          clus = new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersProcess->At(cluslabel));
        else
          clus = new AliESDCaloCluster(*(AliESDCaloCluster*)fInputEvent->GetCaloCluster(cluslabel));
      } else if(fInputEvent->IsA()==AliAODEvent::Class()){
        if(arrClustersProcess)
          clus = new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersProcess->At(cluslabel));
        else
          clus = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(cluslabel));
      }
      if(!clus) return kFALSE;
      TParticle* particle = (TParticle *)fMCEvent->Particle(mclabel);

      // check if its an electron and mother is photon
      Int_t absPdg = TMath::Abs(particle->GetPdgCode());
        if(absPdg == 11){
            Int_t mothLabel = particle->GetMother(0);
            if(mothLabel!= -1){
                TParticle* mother = (TParticle *)fMCEvent->Particle(mothLabel);
                if(mother->GetPdgCode() == 22){
                    if(particle->Energy()/mother->Energy() < 0.5) return kFALSE;
                }
            }
        }

      Int_t *mclabelsCluster = clus->GetLabels();

      if (clus->GetNLabels() > 0)
      {
        for (Int_t k = 0; k < (Int_t)clus->GetNLabels(); k++)
        {
          if (mclabelsCluster[k] == mclabel){
             delete clus;
             return kTRUE;
          }
        }
      }
   // Did not find it, try recusively with daughters
   if(particle->GetNDaughters()>0){
      Bool_t foundIt = kFALSE;
       for (Int_t daughter = particle->GetFirstDaughter(); daughter <= particle->GetLastDaughter(); daughter++)
       {
         foundIt = CheckSpecificClusterForMCContribution(daughter, cluslabel);
         if(foundIt){
           return kTRUE; // break when you found a contribution
         }
       }

    }
   delete clus;
   return kFALSE;
}
Int_t AliAnalysisTaskGammaCalo::CountPhotonsInCluster(Int_t cluslabel)
{
    TClonesArray * arrClustersProcess = NULL;
    fNCurrentClusterBasic             = 0;
    Int_t nphotons = 0;
    if(fCorrTaskSetting.CompareTo("")){
      arrClustersProcess = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));
      if(!arrClustersProcess)
      AliFatal(Form("%sClustersBranch was not found in AliAnalysisTaskGammaCalo! Check the correction framework settings!",fCorrTaskSetting.Data()));
    }

    AliVCluster* clus = NULL;
    if(fInputEvent->IsA()==AliESDEvent::Class()){
      if(arrClustersProcess)
        clus = new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersProcess->At(cluslabel));
      else
        clus = new AliESDCaloCluster(*(AliESDCaloCluster*)fInputEvent->GetCaloCluster(cluslabel));
    } else if(fInputEvent->IsA()==AliAODEvent::Class()){
      if(arrClustersProcess)
        clus = new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersProcess->At(cluslabel));
      else
        clus = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(cluslabel));
    }
    if(!clus) return kFALSE;
    Int_t *mclabelsCluster = clus->GetLabels();

    if (clus->GetNLabels() > 0)
    {
      TParticle* particle = NULL;
      for (Int_t k = 0; k < (Int_t)clus->GetNLabels(); k++)
      {
        particle = (TParticle *)fMCEvent->Particle(mclabelsCluster[k]);
        Int_t pdg = particle->GetPdgCode();
        if(pdg==22) nphotons++;
      }
    }
   delete clus;
   return nphotons;
}
void AliAnalysisTaskGammaCalo::DoClusterMergingStudies(AliVCluster* clus,vector<clusterLabel> &labelvect)
{
  // vertex
  Double_t vertex[3] = {0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  Int_t* mclabelsClus = clus->GetLabels();
  if (clus->GetNLabels()>0){
    for (Int_t k =0; k< (Int_t)clus->GetNLabels(); k++){
      // cluster merging studies
      // get label particle
      AliMCParticle* clusParticle = (AliMCParticle *)fMCEvent->GetTrack(mclabelsClus[k]);
      // set to same as particle to use in while loop
      AliMCParticle* clusParticleMother = (AliMCParticle *)fMCEvent->GetTrack(mclabelsClus[k]);
      // Find out if I find a mother as pi0 somewhere
      Int_t safety = 0;
      Int_t pi0Pos = -1;
      Int_t pi0DaughterPos = mclabelsClus[k];
      while(clusParticleMother->GetMother()!=-1){
        Int_t motherID = clusParticleMother->GetMother();
        clusParticleMother = (AliMCParticle *)fMCEvent->GetTrack(motherID);
        if(clusParticleMother->PdgCode() == 111 || clusParticleMother->PdgCode()==221 ){
            // found pi0
            pi0Pos = motherID;
            break;
        }
        pi0DaughterPos = motherID; // previous label
        safety++;
        if(safety>20) break; // safety to avoid infinite loops
      }
      if(pi0Pos==-1) continue; // label does not belong to pi0
      // Check that we only save a label if it carries more than 50% of its true energy
      // which should be only fulfilled for one label per MC particle

      // make sure its decay to two photons
      if(clusParticleMother->GetNDaughters()!=2) continue;


      AliMCParticle* daughter1 = (AliMCParticle *)fMCEvent->GetTrack(clusParticleMother->GetDaughterFirst());
      AliMCParticle* daughter2 = (AliMCParticle *)fMCEvent->GetTrack(clusParticleMother->GetDaughterLast());

      if( (daughter1->PdgCode() != 22) ||  (daughter2->PdgCode() != 22)) continue;
      TVector3 vec1, vec2;
      vec1.SetXYZ(daughter1->Px(),daughter1->Py(),daughter1->Pz());
      vec2.SetXYZ(daughter2->Px(),daughter2->Py(),daughter2->Pz());

      Double_t angle = 9999;
      angle = vec1.Angle(vec2);

      Double_t EFrac = clus->GetClusterMCEdepFraction(k);
      Double_t EClus = clus->E();
      Double_t ETrue = clusParticle->E();
      Double_t FracDepos = (EFrac * EClus)/ ETrue;
      if(FracDepos<=0.5) continue;

      // check for electron with photon mother
      AliMCParticle* tmpMoth = NULL;
      AliMCParticle* pi0Photon = NULL;
      if(TMath::Abs(clusParticle->PdgCode())==11){
          tmpMoth = (AliMCParticle *)fMCEvent->GetTrack(clusParticle->GetMother());
          pi0Photon = (AliMCParticle *)fMCEvent->GetTrack(pi0DaughterPos);
          Int_t pdgMoth = tmpMoth->PdgCode();
          if(pdgMoth!=22) continue; // we only consider labels of electrons if they come from gamma
          // now check if the electron carries at least 50 percent of energy
          // of the mother photon right after pi0 (in case of conv or shower)
          if((clusParticle->E()/pi0Photon->E())<=0.8){
              continue;
          }
      }

      // TLorentzvector with cluster
      TLorentzVector clusterVector;
      clus->GetMomentum(clusterVector,vertex);
      // create cluster label object
      clusterLabel tmpLabel;
      tmpLabel.mesonID = pi0Pos;
      tmpLabel.clusID= clus->GetID();
      tmpLabel.daughterPDG = clusParticle->PdgCode();
      tmpLabel.daughterID = pi0DaughterPos; // always store the id of the particle right after pi0
      tmpLabel.EClus = EClus;
      tmpLabel.EFrac = EFrac;
      tmpLabel.ETrue = ETrue;
      tmpLabel.PtMeson = clusParticleMother->Pt();
      tmpLabel.EtaMeson = clusParticleMother->Eta();
      tmpLabel.OpeningAngle = angle;
      tmpLabel.clusVec = clusterVector;


      labelvect.push_back(tmpLabel);
    } // end of label loop
  }
}
void AliAnalysisTaskGammaCalo::DoClusterMergingStudiesAOD(AliVCluster* clus,vector<clusterLabel> &labelvect)
{
  Int_t* mclabelsClus = clus->GetLabels();
  if (clus->GetNLabels()>0){
    for (Int_t k =0; k< (Int_t)clus->GetNLabels(); k++){
      // cluster merging studies
      // get label particle
      AliAODMCParticle* clusParticle = (AliAODMCParticle *)fAODMCTrackArray->At(mclabelsClus[k]);
      // set to same as particle to use in while loop
      AliAODMCParticle* clusParticleMother = (AliAODMCParticle *)fAODMCTrackArray->At(mclabelsClus[k]);
      // Find out if I find a mother as pi0 somewhere
      Int_t safety = 0;
      Int_t pi0Pos = -1;
      Int_t pi0DaughterPos = mclabelsClus[k];
      while(clusParticleMother->GetMother()!=-1){
        Int_t motherID = clusParticleMother->GetMother();
        clusParticleMother = (AliAODMCParticle *)fAODMCTrackArray->At(motherID);
        if(clusParticleMother->GetPdgCode() == 111 || clusParticleMother->GetPdgCode() == 221){
            // found pi0
            pi0Pos = motherID;
            break;
        }
        pi0DaughterPos = motherID; // previous label
        safety++;
        if(safety>20) break; // safety to avoid infinite loops
      }
      if(pi0Pos==-1) continue; // label does not belong to pi0

      // make sure its decay to two photons
      if(clusParticleMother->GetNDaughters()!=2) continue;


      AliAODMCParticle* daughter1 = (AliAODMCParticle *)fAODMCTrackArray->At(clusParticleMother->GetDaughterFirst());
      AliAODMCParticle* daughter2 = (AliAODMCParticle *)fAODMCTrackArray->At(clusParticleMother->GetDaughterLast());

      if( (daughter1->GetPdgCode() != 22) ||  (daughter2->GetPdgCode() != 22)) continue;
      TVector3 vec1, vec2;
      vec1.SetXYZ(daughter1->Px(),daughter1->Py(),daughter1->Pz());
      vec2.SetXYZ(daughter2->Px(),daughter2->Py(),daughter2->Pz());

      Double_t angle = 9999;
      angle = vec1.Angle(vec2);

      // Check that we only save a label if it carries more than 50% of its true energy
      // which should be only fulfilled for one label per MC particle
      Double_t EFrac = clus->GetClusterMCEdepFraction(k);
      Double_t EClus = clus->E();
      Double_t ETrue = clusParticle->E();
      Double_t FracDepos = (EFrac * EClus)/ ETrue;
      if(FracDepos<=0.5) continue;

      // check for electron with photon mother
      AliAODMCParticle* tmpMoth = NULL;
      AliAODMCParticle* pi0Photon = NULL;
      if(TMath::Abs(clusParticle->GetPdgCode())==11){
          tmpMoth = (AliAODMCParticle *)fAODMCTrackArray->At(clusParticle->GetMother());
          pi0Photon = (AliAODMCParticle *)fAODMCTrackArray->At(pi0DaughterPos);
          Int_t pdgMoth = tmpMoth->GetPdgCode();
          if(pdgMoth!=22) continue; // we only consider labels of electrons if they come from gamma
          // now check if the electron carries at least 50 percent of energy
          // of the mother photon right after pi0 (in case of conv or shower)
          if((clusParticle->E()/pi0Photon->E())<=0.8){
              continue;
          }
      }

      // create cluster label object
      clusterLabel tmpLabel;
      tmpLabel.mesonID = pi0Pos;
      tmpLabel.clusID= clus->GetID();
      tmpLabel.daughterPDG = clusParticle->GetPdgCode();
      tmpLabel.daughterID = pi0DaughterPos; // always store the id of the particle right after pi0
      tmpLabel.EClus = EClus;
      tmpLabel.EFrac = EFrac;
      tmpLabel.ETrue = ETrue;
      tmpLabel.PtMeson = clusParticleMother->Pt();
      tmpLabel.EtaMeson = clusParticleMother->Eta();
      tmpLabel.OpeningAngle = angle;

      labelvect.push_back(tmpLabel);
    } // end of label loop
  }
}
