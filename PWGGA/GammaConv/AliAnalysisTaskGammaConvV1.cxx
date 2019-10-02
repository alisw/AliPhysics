/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                    *
* Author: Martin Wilde, Daniel Lohner, Friederike Bock               *
* Version 1.0                  *
*                    *
* based on: on older version (see aliroot up to v5-04-42-AN)             *
*           AliAnalysisTaskGammaConversion.cxx                           *
*           Authors: Kathrin Koch, Kenneth Aamodt, Ana Marin             *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its    *
* documentation strictly for non-commercial purposes is hereby granted    *
* without fee, provided that the above copyright notice appears in all    *
* copies and that both the copyright notice and this permission notice    *
* appear in the supporting documentation. The authors make no claims    *
* about the suitability of this software for any purpose. It is    *
* provided "as is" without express or implied warranty.             *
**************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------

// Class used to do analysis on conversion pairs
//---------------------------------------------
///////////////////////////////////////////////
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
#include "AliAnalysisTaskGammaConvV1.h"
#include "AliVParticle.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliKFVertex.h"
#include "AliGenCocktailEventHeader.h"
#include "AliConversionAODBGHandlerRP.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliEventplane.h"
#include "AliAODEvent.h"
#include <vector>
#include <map>


ClassImp(AliAnalysisTaskGammaConvV1)

//________________________________________________________________________
AliAnalysisTaskGammaConvV1::AliAnalysisTaskGammaConvV1(): AliAnalysisTaskSE(),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fDoLightOutput(kFALSE),
  fBGHandler(NULL),
  fBGHandlerRP(NULL),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fCutFolder(NULL),
  fESDList(NULL),
  fBackList(NULL),
  fMotherList(NULL),
  fTrueList(NULL),
  fMCList(NULL),
  fHeaderNameList(NULL),
  fOutputContainer(0),
  fReaderGammas(NULL),
  fGammaCandidates(NULL),
  fEventCutArray(NULL),
  fCutArray(NULL),
  fMesonCutArray(NULL),
  fClusterCutArray(NULL),
  fConvJetReader(NULL),
  fDoJetAnalysis(kFALSE),
  fDoJetQA(kFALSE),
  fJetHistograms(NULL),
  fTrueJetHistograms(NULL),
  fMaxPtNearEMCalPlace(0),
  fJetNearEMCal(kFALSE),
  fHistoCaloGammaPt(NULL),
  fHistoCaloGammaE(NULL),
  fHistoConvGammaPt(NULL),
  fHistoConvGammaR(NULL),
  fHistoConvGammaEta(NULL),
  fHistoConvGammaPhi(NULL),
  fHistoConvGammaPsiPairPt(NULL),
  fHistoConvGammaInvMass(NULL),
  fHistoConvGammaInvMassReco(NULL),
  tESDConvGammaPtDcazCat(NULL),
  fPtGamma(0),
  fDCAzPhoton(0),
  fRConvPhoton(0),
  fEtaPhoton(0),
  iCatPhoton(0),
  iPhotonMCInfo(0),
  fHistoMotherInvMassPt(NULL),
  sESDMotherInvMassPtZM(NULL),
  fHistoMotherBackInvMassPt(NULL),
  sESDMotherBackInvMassPtZM(NULL),
  fHistoMotherInvMassEalpha(NULL),
  fHistoMotherPi0PtY(NULL),
  fHistoMotherEtaPtY(NULL),
  fHistoMotherPi0PtAlpha(NULL),
  fHistoMotherEtaPtAlpha(NULL),
  fHistoMotherPi0PtOpenAngle(NULL),
  fHistoMotherEtaPtOpenAngle(NULL),
  sPtRDeltaROpenAngle(NULL),
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
  fHistoMCConvGammaPt(NULL),
  fHistoMCSecondaryConvGammaPt(NULL),
  fHistoMCConvGammaR(NULL),
  fHistoMCConvGammaEta(NULL),
  fHistoMCPi0Pt(NULL),
  fHistoMCPi0WOWeightPt(NULL),
  fHistoMCPi0WOEvtWeightPt(NULL),
  fHistoMCEtaWOEvtWeightPt(NULL),
  fHistoMCEtaPt(NULL),
  fHistoMCEtaWOWeightPt(NULL),
  fHistoMCPi0WOWeightInAccPt(NULL),
  fHistoMCEtaWOWeightInAccPt(NULL),
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
  fHistoMCSecPi0RvsSource(NULL),
  fHistoMCSecPi0Source(NULL),
  fHistoMCSecPi0InAccPtvsSource(NULL),
  fHistoMCSecEtaPt(NULL),
  fHistoMCSecEtaSource(NULL),
  fHistoMCPi0PtJetPt(NULL),
  fHistoMCEtaPtJetPt(NULL),
  fHistoMCPhysicalPrimariesPt(NULL),
  fHistoTrueMotherInvMassPt(NULL),
  fHistoTruePrimaryMotherInvMassPt(NULL),
  fHistoTruePrimaryMotherW0WeightingInvMassPt(NULL),
  pESDTruePrimaryMotherWeightsInvMassPt(NULL),
  fHistoTruePrimaryPi0MCPtResolPt(NULL),
  fHistoTruePrimaryEtaMCPtResolPt(NULL),
  fHistoTrueSecondaryMotherInvMassPt(NULL),
  fHistoTrueSecondaryMotherFromK0sInvMassPt(NULL),
  fHistoTrueK0sWithPi0DaughterMCPt(NULL),
  fHistoTrueSecondaryMotherFromK0lInvMassPt(NULL),
  fHistoTrueK0lWithPi0DaughterMCPt(NULL),
  fHistoTrueSecondaryMotherFromEtaInvMassPt(NULL),
  fHistoTrueEtaWithPi0DaughterMCPt(NULL),
  fHistoTrueSecondaryMotherFromLambdaInvMassPt(NULL),
  fHistoTrueLambdaWithPi0DaughterMCPt(NULL),
  fHistoTrueBckGGInvMassPt(NULL),
  fHistoTrueBckContInvMassPt(NULL),
  fHistoTruePi0PtY(NULL),
  fHistoTrueEtaPtY(NULL),
  fHistoTruePi0PtAlpha(NULL),
  fHistoTrueEtaPtAlpha(NULL),
  fHistoTruePi0PtOpenAngle(NULL),
  fHistoTrueEtaPtOpenAngle(NULL),
  fHistoTrueMotherDalitzInvMassPt(NULL),
  fHistoTrueConvGammaPt(NULL),
  fHistoTrueConvGammaR(NULL),
  fHistoTrueConvGammaPtMC(NULL),
  fHistoTrueConvGammaRMC(NULL),
  fHistoTrueConvGammaEta(NULL),
  fHistoTrueConvGammaPsiPairPt(NULL),
  fHistoTrueConvGammaInvMass(NULL),
  fHistoTrueConvGammaInvMassReco(NULL),
  fHistoCombinatorialPt(NULL),
  fHistoCombinatorialMothersPt(NULL),
  fHistoCombinatorialPtDeltaPhi_ek(NULL),
  fHistoCombinatorialPtDeltaPhi_ep(NULL),
  fHistoCombinatorialPtDeltaPhi_epi(NULL),
  fHistoCombinatorialPtDeltaPhi_pik(NULL),
  fHistoCombinatorialPtDeltaPhi_pip(NULL),
  fHistoTruePrimaryConvGammaPt(NULL),
  fHistoTrueSecondaryConvGammaPt(NULL),
  fHistoTrueSecondaryConvGammaMCPt(NULL),
  fHistoTruePrimaryConvGammaESDPtMCPt(NULL),
  fHistoTrueSecondaryConvGammaFromXFromK0sMCPtESDPt(NULL),
  fHistoTrueSecondaryConvGammaFromXFromK0lMCPtESDPt(NULL),
  fHistoTrueSecondaryConvGammaFromXFromLambdaMCPtESDPt(NULL),
  fHistoTrueDalitzPsiPairDeltaPhi(NULL),
  fHistoTrueGammaPsiPairDeltaPhi(NULL),
  fHistoDoubleCountTruePi0InvMassPt(NULL),
  fHistoDoubleCountTrueEtaInvMassPt(NULL),
  fHistoDoubleCountTrueConvGammaRPt(NULL),
  vecDoubleCountTruePi0s(0),
  vecDoubleCountTrueEtas(0),
  vecDoubleCountTrueConvGammas(0),
  fHistoMultipleCountTruePi0(NULL),
  fHistoMultipleCountTrueEta(NULL),
  fHistoMultipleCountTrueConvGamma(NULL),
  mapMultipleCountTruePi0s(),
  mapMultipleCountTrueEtas(),
  mapMultipleCountTrueConvGammas(),
  fHistoNEvents(NULL),
  fHistoNEventsWOWeight(NULL),
  fHistoNGoodESDTracks(NULL),
  fHistoNEventsWeighted(NULL),
  fHistoNGoodESDTracksWeighted(NULL),
  fHistoVertexZ(NULL),
  fHistoVertexZWeighted(NULL),
  fDoCentralityFlat(0),
  fHistoCentrality(NULL),
  fHistoCentralityFlattened(NULL),
  fHistoCentralityVsPrimaryTracks(NULL),
  fHistoNGammaCandidates(NULL),
  fHistoNGoodESDTracksVsNGammaCandidates(NULL),
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
  fHistoTrueJetMotherInvMassPt(NULL),
  fHistoTrueInJetMotherInvMassPt(NULL),
  fHistoTruePrimaryJetInvMassPt(NULL),
  fHistoTruePrimaryinJetInvMassPt(NULL),
  fHistoTruePrimaryInJetInvMassTruePt(NULL),
  fHistoTrueDoubleCountingJet(NULL),
  fHistoTrueJetFragmFunc(NULL),
  fHistoTrueJetFragmFuncZInvMass(NULL),
  fHistoMCPi0JetInAccPt(NULL),
  fHistoMCPi0inJetInAccPt(NULL),
  fHistoMCEtaJetInAccPt(NULL),
  fHistoMCEtainJetInAccPt(NULL),
  fHistoMCPi0JetEventGenerated(NULL),
  fHistoMCPi0inJetGenerated(NULL),
  fHistoMCEtaJetEventGenerated(NULL),
  fHistoMCEtainJetGenerated(NULL),
  fHistoTrueSecondaryFromK0sJetInvMassPt(NULL),
  fHistoTrueSecondaryFromK0sinJetInvMassPt(NULL),
  fHistoTrueSecondaryFromLambdaJetInvMassPt(NULL),
  fHistoTrueSecondaryFromLambdainJetInvMassPt(NULL),
  fHistoTrueSecondaryFromK0lJetInvMassPt(NULL),
  fHistoTrueSecondaryFromK0linJetInvMassPt(NULL),
  fHistoTrueSecondaryInvJetMassPt(NULL),
  fHistoTrueSecondaryInvinJetMassPt(NULL),
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
  fHistoSPDClusterTrackletBackground(NULL),
  fHistoV0MultVsNumberTPCoutTracks(NULL),
  fHistoNV0Tracks(NULL),
  fHistoBDToutput(NULL),
  fHistoBDToutputPt(NULL),
  fHistoBDToutputMCTrue(NULL),
  fProfileEtaShift(NULL),
  fProfileJetJetXSection(NULL),
  fhJetJetNTrials(NULL),
  fHistoEventSphericity(NULL),
  fHistoEtaShift(NULL),
  tESDMesonsInvMassPtDcazMinDcazMaxFlag(NULL),
  fInvMass(0),
  fBDTvariable(NULL),
  fBDTreader(NULL),
  fFileNameBDT(0),
  fEnableBDT(0),
  fPt(0),
  fDCAzGammaMin(0),
  fDCAzGammaMax(0),
  iFlag(0),
  iMesonMCInfo(0),
  fEventPlaneAngle(-100),
  fRandom(0),
  fnGammaCandidates(0),
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
  fDoMesonAnalysis(kTRUE),
  fDoMesonQA(0),
  fDoPhotonQA(0),
  fDoChargedPrimary(kFALSE),
  fDoPlotVsCentrality(kFALSE),
  fIsFromSelectedHeader(kTRUE),
  fIsMC(0),
  fDoTHnSparse(kFALSE),
  fWeightJetJetMC(1),
  fWeightCentrality(NULL),
  fEnableClusterCutsForTrigger(kFALSE),
  fDoMaterialBudgetWeightingOfGammasForTrueMesons(kFALSE),
  tBrokenFiles(NULL),
  fFileNameBroken(NULL)
{

}

//________________________________________________________________________
AliAnalysisTaskGammaConvV1::AliAnalysisTaskGammaConvV1(const char *name):
  AliAnalysisTaskSE(name),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fDoLightOutput(kFALSE),
  fBGHandler(NULL),
  fBGHandlerRP(NULL),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fCutFolder(NULL),
  fESDList(NULL),
  fBackList(NULL),
  fMotherList(NULL),
  fTrueList(NULL),
  fMCList(NULL),
  fHeaderNameList(NULL),
  fOutputContainer(0),
  fReaderGammas(NULL),
  fGammaCandidates(NULL),
  fEventCutArray(NULL),
  fCutArray(NULL),
  fMesonCutArray(NULL),
  fClusterCutArray(NULL),
  fConvJetReader(NULL),
  fDoJetAnalysis(kFALSE),
  fDoJetQA(kFALSE),
  fJetHistograms(NULL),
  fTrueJetHistograms(NULL),
  fMaxPtNearEMCalPlace(0),
  fJetNearEMCal(kFALSE),
  fHistoCaloGammaPt(NULL),
  fHistoCaloGammaE(NULL),
  fHistoConvGammaPt(NULL),
  fHistoConvGammaR(NULL),
  fHistoConvGammaEta(NULL),
  fHistoConvGammaPhi(NULL),
  fHistoConvGammaPsiPairPt(NULL),
  fHistoConvGammaInvMass(NULL),
  fHistoConvGammaInvMassReco(NULL),
  tESDConvGammaPtDcazCat(NULL),
  fPtGamma(0),
  fDCAzPhoton(0),
  fRConvPhoton(0),
  fEtaPhoton(0),
  iCatPhoton(0),
  iPhotonMCInfo(0),
  fHistoMotherInvMassPt(NULL),
  sESDMotherInvMassPtZM(NULL),
  fHistoMotherBackInvMassPt(NULL),
  sESDMotherBackInvMassPtZM(NULL),
  fHistoMotherInvMassEalpha(NULL),
  fHistoMotherPi0PtY(NULL),
  fHistoMotherEtaPtY(NULL),
  fHistoMotherPi0PtAlpha(NULL),
  fHistoMotherEtaPtAlpha(NULL),
  fHistoMotherPi0PtOpenAngle(NULL),
  fHistoMotherEtaPtOpenAngle(NULL),
  sPtRDeltaROpenAngle(NULL),
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
  fHistoMCConvGammaPt(NULL),
  fHistoMCSecondaryConvGammaPt(NULL),
  fHistoMCConvGammaR(NULL),
  fHistoMCConvGammaEta(NULL),
  fHistoMCPi0Pt(NULL),
  fHistoMCPi0WOWeightPt(NULL),
  fHistoMCPi0WOEvtWeightPt(NULL),
  fHistoMCEtaWOEvtWeightPt(NULL),
  fHistoMCEtaPt(NULL),
  fHistoMCEtaWOWeightPt(NULL),
  fHistoMCPi0WOWeightInAccPt(NULL),
  fHistoMCEtaWOWeightInAccPt(NULL),
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
  fHistoMCSecPi0RvsSource(NULL),
  fHistoMCSecPi0Source(NULL),
  fHistoMCSecPi0InAccPtvsSource(NULL),
  fHistoMCSecEtaPt(NULL),
  fHistoMCSecEtaSource(NULL),
  fHistoMCPi0PtJetPt(NULL),
  fHistoMCEtaPtJetPt(NULL),
  fHistoMCPhysicalPrimariesPt(NULL),
  fHistoTrueMotherInvMassPt(NULL),
  fHistoTruePrimaryMotherInvMassPt(NULL),
  fHistoTruePrimaryMotherW0WeightingInvMassPt(NULL),
  pESDTruePrimaryMotherWeightsInvMassPt(NULL),
  fHistoTruePrimaryPi0MCPtResolPt(NULL),
  fHistoTruePrimaryEtaMCPtResolPt(NULL),
  fHistoTrueSecondaryMotherInvMassPt(NULL),
  fHistoTrueSecondaryMotherFromK0sInvMassPt(NULL),
  fHistoTrueK0sWithPi0DaughterMCPt(NULL),
  fHistoTrueSecondaryMotherFromK0lInvMassPt(NULL),
  fHistoTrueK0lWithPi0DaughterMCPt(NULL),
  fHistoTrueSecondaryMotherFromEtaInvMassPt(NULL),
  fHistoTrueEtaWithPi0DaughterMCPt(NULL),
  fHistoTrueSecondaryMotherFromLambdaInvMassPt(NULL),
  fHistoTrueLambdaWithPi0DaughterMCPt(NULL),
  fHistoTrueBckGGInvMassPt(NULL),
  fHistoTrueBckContInvMassPt(NULL),
  fHistoTruePi0PtY(NULL),
  fHistoTrueEtaPtY(NULL),
  fHistoTruePi0PtAlpha(NULL),
  fHistoTrueEtaPtAlpha(NULL),
  fHistoTruePi0PtOpenAngle(NULL),
  fHistoTrueEtaPtOpenAngle(NULL),
  fHistoTrueMotherDalitzInvMassPt(NULL),
  fHistoTrueConvGammaPt(NULL),
  fHistoTrueConvGammaR(NULL),
  fHistoTrueConvGammaPtMC(NULL),
  fHistoTrueConvGammaRMC(NULL),
  fHistoTrueConvGammaEta(NULL),
  fHistoTrueConvGammaPsiPairPt(NULL),
  fHistoTrueConvGammaInvMass(NULL),
  fHistoTrueConvGammaInvMassReco(NULL),
  fHistoCombinatorialPt(NULL),
  fHistoCombinatorialMothersPt(NULL),
  fHistoCombinatorialPtDeltaPhi_ek(NULL),
  fHistoCombinatorialPtDeltaPhi_ep(NULL),
  fHistoCombinatorialPtDeltaPhi_epi(NULL),
  fHistoCombinatorialPtDeltaPhi_pik(NULL),
  fHistoCombinatorialPtDeltaPhi_pip(NULL),
  fHistoTruePrimaryConvGammaPt(NULL),
  fHistoTrueSecondaryConvGammaPt(NULL),
  fHistoTrueSecondaryConvGammaMCPt(NULL),
  fHistoTruePrimaryConvGammaESDPtMCPt(NULL),
  fHistoTrueSecondaryConvGammaFromXFromK0sMCPtESDPt(NULL),
  fHistoTrueSecondaryConvGammaFromXFromK0lMCPtESDPt(NULL),
  fHistoTrueSecondaryConvGammaFromXFromLambdaMCPtESDPt(NULL),
  fHistoTrueDalitzPsiPairDeltaPhi(NULL),
  fHistoTrueGammaPsiPairDeltaPhi(NULL),
  fHistoDoubleCountTruePi0InvMassPt(NULL),
  fHistoDoubleCountTrueEtaInvMassPt(NULL),
  fHistoDoubleCountTrueConvGammaRPt(NULL),
  vecDoubleCountTruePi0s(0),
  vecDoubleCountTrueEtas(0),
  vecDoubleCountTrueConvGammas(0),
  fHistoMultipleCountTruePi0(NULL),
  fHistoMultipleCountTrueEta(NULL),
  fHistoMultipleCountTrueConvGamma(NULL),
  mapMultipleCountTruePi0s(),
  mapMultipleCountTrueEtas(),
  mapMultipleCountTrueConvGammas(),
  fHistoNEvents(NULL),
  fHistoNEventsWOWeight(NULL),
  fHistoNGoodESDTracks(NULL),
  fHistoNEventsWeighted(NULL),
  fHistoNGoodESDTracksWeighted(NULL),
  fHistoVertexZ(NULL),
  fHistoVertexZWeighted(NULL),
  fDoCentralityFlat(0),
  fHistoCentrality(NULL),
  fHistoCentralityFlattened(NULL),
  fHistoCentralityVsPrimaryTracks(NULL),
  fHistoNGammaCandidates(NULL),
  fHistoNGoodESDTracksVsNGammaCandidates(NULL),
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
  fHistoTrueJetMotherInvMassPt(NULL),
  fHistoTrueInJetMotherInvMassPt(NULL),
  fHistoTruePrimaryJetInvMassPt(NULL),
  fHistoTruePrimaryinJetInvMassPt(NULL),
  fHistoTruePrimaryInJetInvMassTruePt(NULL),
  fHistoTrueDoubleCountingJet(NULL),
  fHistoTrueJetFragmFunc(NULL),
  fHistoTrueJetFragmFuncZInvMass(NULL),
  fHistoMCPi0JetInAccPt(NULL),
  fHistoMCPi0inJetInAccPt(NULL),
  fHistoMCEtaJetInAccPt(NULL),
  fHistoMCEtainJetInAccPt(NULL),
  fHistoMCPi0JetEventGenerated(NULL),
  fHistoMCPi0inJetGenerated(NULL),
  fHistoMCEtaJetEventGenerated(NULL),
  fHistoMCEtainJetGenerated(NULL),
  fHistoTrueSecondaryFromK0sJetInvMassPt(NULL),
  fHistoTrueSecondaryFromK0sinJetInvMassPt(NULL),
  fHistoTrueSecondaryFromLambdaJetInvMassPt(NULL),
  fHistoTrueSecondaryFromLambdainJetInvMassPt(NULL),
  fHistoTrueSecondaryFromK0lJetInvMassPt(NULL),
  fHistoTrueSecondaryFromK0linJetInvMassPt(NULL),
  fHistoTrueSecondaryInvJetMassPt(NULL),
  fHistoTrueSecondaryInvinJetMassPt(NULL),
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
  fHistoSPDClusterTrackletBackground(NULL),
  fHistoV0MultVsNumberTPCoutTracks(NULL),
  fHistoNV0Tracks(NULL),
  fHistoBDToutput(NULL),
  fHistoBDToutputPt(NULL),
  fHistoBDToutputMCTrue(NULL),
  fProfileEtaShift(NULL),
  fProfileJetJetXSection(NULL),
  fhJetJetNTrials(NULL),
  fHistoEventSphericity(NULL),
  fHistoEtaShift(NULL),
  tESDMesonsInvMassPtDcazMinDcazMaxFlag(NULL),
  fInvMass(0),
  fBDTvariable(NULL),
  fBDTreader(NULL),
  fFileNameBDT(0),
  fEnableBDT(0),
  fPt(0),
  fDCAzGammaMin(0),
  fDCAzGammaMax(0),
  iFlag(0),
  iMesonMCInfo(0),
  fEventPlaneAngle(-100),
  fRandom(0),
  fnGammaCandidates(0),
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
  fDoMesonAnalysis(kTRUE),
  fDoMesonQA(0),
  fDoPhotonQA(0),
  fDoChargedPrimary(kFALSE),
  fDoPlotVsCentrality(kFALSE),
  fIsFromSelectedHeader(kTRUE),
  fIsMC(0),
  fDoTHnSparse(kFALSE),
  fWeightJetJetMC(1),
  fWeightCentrality(NULL),
  fEnableClusterCutsForTrigger(kFALSE),
  fDoMaterialBudgetWeightingOfGammasForTrueMesons(kFALSE),
  tBrokenFiles(NULL),
  fFileNameBroken(NULL)
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

AliAnalysisTaskGammaConvV1::~AliAnalysisTaskGammaConvV1()
{
  if(fGammaCandidates){
    delete fGammaCandidates;
    fGammaCandidates = 0x0;
  }
  if(fBGHandler){
    delete[] fBGHandler;
    fBGHandler = 0x0;
  }
  if(fBGHandlerRP){
    delete[] fBGHandlerRP;
    fBGHandlerRP = 0x0;
  }

  if(fWeightCentrality){
    delete[] fWeightCentrality;
    fWeightCentrality = 0x0;
  }

}
//___________________________________________________________
void AliAnalysisTaskGammaConvV1::InitBack(){

  const Int_t nDim = 4;
  Int_t nBins[nDim] = {800,250,7,4};
  Double_t xMin[nDim] = {0,0, 0,0};
  Double_t xMax[nDim] = {0.8,25,7,4};
  Int_t nBinsRP[nDim] = {800,250,7,8};
  Double_t xMinRP[nDim] = {0,0, 0,0};
  Double_t xMaxRP[nDim] = {0.8,25,7,8};

  if(fDoTHnSparse){
    sESDMotherInvMassPtZM = new THnSparseF*[fnCuts];
    sESDMotherBackInvMassPtZM = new THnSparseF*[fnCuts];
  }
  fBGHandler = new AliGammaConversionAODBGHandler*[fnCuts];
  fBGHandlerRP = new AliConversionAODBGHandlerRP*[fnCuts];
  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if (((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->DoBGCalculation()){
      TString cutstringEvent   = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
      TString cutstringPhoton = ((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutNumber();
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
      } else if(collisionSystem == 3 || collisionSystem == 6) {
        centMin = centMin*5;
        centMax = centMax*5;
      } else if(collisionSystem == 4 || collisionSystem == 7) {
        centMin = ((centMin*5)+45);
        centMax = ((centMax*5)+45);
      }

      if(fDoTHnSparse){
        fBackList[iCut] = new TList();
        fBackList[iCut]->SetName(Form("%s_%s_%s Back histograms",cutstringEvent.Data(), cutstringPhoton.Data(),cutstringMeson.Data()));
        fBackList[iCut]->SetOwner(kTRUE);
        fCutFolder[iCut]->Add(fBackList[iCut]);

        if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->BackgroundHandlerType() == 0){
          sESDMotherBackInvMassPtZM[iCut] = new THnSparseF("Back_Back_InvMass_Pt_z_m", "Back_Back_InvMass_Pt_z_m",nDim,nBins,xMin,xMax);
        } else {
          sESDMotherBackInvMassPtZM[iCut] = new THnSparseF("Back_Back_InvMass_Pt_z_psi", "Back_Back_InvMass_Pt_z_psi",nDim,nBinsRP,xMinRP,xMaxRP);
        }
        if(fDoCentralityFlat > 0) sESDMotherBackInvMassPtZM[iCut]->Sumw2();
        fBackList[iCut]->Add(sESDMotherBackInvMassPtZM[iCut]);

        fMotherList[iCut] = new TList();
        fMotherList[iCut]->SetName(Form("%s_%s_%s Mother histograms",cutstringEvent.Data() ,cutstringPhoton.Data(),cutstringMeson.Data()));
        fMotherList[iCut]->SetOwner(kTRUE);
        fCutFolder[iCut]->Add(fMotherList[iCut]);

        if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->BackgroundHandlerType() == 0){
          sESDMotherInvMassPtZM[iCut] = new THnSparseF("Back_Mother_InvMass_Pt_z_m", "Back_Mother_InvMass_Pt_z_m",nDim,nBins,xMin,xMax);
        } else {
          sESDMotherInvMassPtZM[iCut] = new THnSparseF("Back_Mother_InvMass_Pt_z_psi", "Back_Mother_InvMass_Pt_z_psi",nDim,nBinsRP,xMinRP,xMaxRP);
        }
        if(fDoCentralityFlat > 0) sESDMotherInvMassPtZM[iCut]->Sumw2();
        fMotherList[iCut]->Add(sESDMotherInvMassPtZM[iCut]);
      }
      if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->BackgroundHandlerType() == 0){
        fBGHandler[iCut] = new AliGammaConversionAODBGHandler(
                                  collisionSystem,centMin,centMax,
                                  ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetNumberOfBGEvents(),
                                  ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseTrackMultiplicity(),
                                  0,8,5);
        fBGHandlerRP[iCut] = NULL;
      } else {
        fBGHandlerRP[iCut] = new AliConversionAODBGHandlerRP(
                                  ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsHeavyIon(),
                                  ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity(),
                                  ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetNumberOfBGEvents());
        fBGHandler[iCut] = NULL;
      }
    }
  }
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::UserCreateOutputObjects(){

  // get V0 Reader
  fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
  if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader


  if( ((AliConversionPhotonCuts*)fCutArray->At(0))->GetUseBDTPhotonCuts()){
      fEnableBDT  = kTRUE;
      InitializeBDT(); //
  }

  if(fDoMesonAnalysis){ //Same Jet Finder MUST be used within same trainconfig
    if( ((AliConversionMesonCuts*)fMesonCutArray->At(0))->DoJetAnalysis())  fDoJetAnalysis = kTRUE;
    if( ((AliConversionMesonCuts*)fMesonCutArray->At(0))->DoJetQA())        fDoJetQA       = kTRUE;
  }

  if(fDoJetAnalysis){
    fConvJetReader=(AliAnalysisTaskConvJet*)AliAnalysisManager::GetAnalysisManager()->GetTask("AliAnalysisTaskConvJet");
    if(!fConvJetReader){printf("Error: No AliAnalysisTaskConvJet");return;} // GetV0Reader
  }

  // Set dimenstions for THnSparse
  const Int_t nDim2 = 4;
  Int_t nBins2[nDim2]         = {250,180,100,100};
  Double_t xMin2[nDim2]       = {0,0, 0,0};
  Double_t xMax2[nDim2]       = {25,180,10,0.1};

  // Set maximum number of tracks, gammas and tracklets
  Int_t nGammaCand              = 50;
  Int_t nTracks                 = 200;
  Int_t nSPDTracklets           = 200;
  Int_t nSPDClusters            = 400;
  Int_t nBinsTrklCls            = 100;
  if(fIsHeavyIon == 1){
    nGammaCand                  = 200;
    nTracks                     = 4000;
    nSPDTracklets               = 6000;
    nSPDClusters                = 20000;
    nBinsTrklCls                = 200;
  } else if(fIsHeavyIon == 2){
    nGammaCand                  = 50;
    nTracks                     = 400;
    nSPDTracklets               = 300;
    nSPDClusters                = 600;
    nBinsTrklCls                = 150;
  }

//  Float_t binWidthPt          = 0.1;
  Int_t nBinsPt               = 250;
  Float_t minPt               = 0;
  Float_t maxPt               = 25;
  Int_t nBinsQAPt             = 175;
  Float_t maxQAPt             = 25;
  Int_t nBinsClusterPt        = 500;
  Float_t minClusterPt        = 0;
  Float_t maxClusterPt        = 50;
  Double_t *arrPtBinning      = new Double_t[1200];
  Double_t *arrQAPtBinning    = new Double_t[1200];
  Double_t *arrClusPtBinning  = new Double_t[1200];

  // Set special pt binning for pPb 5TeV, pPb 8TeV
  if (((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::kpPb5TeV ||
      ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::kpPb5TeVR2  ){
//    binWidthPt                = 0.05;
    nBinsPt                   = 186;
    minPt                     = 0;
    maxPt                     = 40;
    for(Int_t i=0; i<nBinsPt+1;i++){
      if (i < 1) arrPtBinning[i]              = 0.25*i;
      else if(i<56) arrPtBinning[i]           = 0.25+0.05*(i-1);
      else if(i<126) arrPtBinning[i]          = 3.+0.1*(i-56);
      else if(i<166) arrPtBinning[i]          = 10.+0.25*(i-126);
      else if(i<186) arrPtBinning[i]          = 20.+1.0*(i-166);
      else arrPtBinning[i]                    = maxPt;
    }
    nBinsQAPt                 = 190;
    maxQAPt                   = 40;
    for(Int_t i=0; i<nBinsQAPt+1;i++){
      if(i<60) arrQAPtBinning[i]              = 0.05*i;
      else if(i<130) arrQAPtBinning[i]        = 3.+0.1*(i-60);
      else if(i<170) arrQAPtBinning[i]        = 10.+0.25*(i-130);
      else if(i<190) arrQAPtBinning[i]        = 20.+1.0*(i-170);
      else arrQAPtBinning[i]                  = maxQAPt;
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
  // Set special pt binning for pp 13TeV and 5TeV
  } else if ( ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::k13TeV ||
              ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::k13TeVLowB ||
              ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::k5TeV ){
//    binWidthPt                = 0.05;
    nBinsPt                   = 169;
    minPt                     = 0;
    maxPt                     = 50;
    for(Int_t i=0; i<nBinsPt+1;i++){
      if (i < 100) arrPtBinning[i]            = 0.10*i;
      else if(i<140) arrPtBinning[i]          = 10.+0.25*(i-100);
      else if(i<170) arrPtBinning[i]          = 20.+1.00*(i-140);
      else arrPtBinning[i]                    = maxPt;
    }
    nBinsQAPt                 = 170;
    maxQAPt                   = 50;
    for(Int_t i=0; i<nBinsQAPt+1;i++){
      if(i<100) arrQAPtBinning[i]             = 0.1*i;
      else if(i<140) arrQAPtBinning[i]        = 10.+0.25*(i-100);
      else if(i<170) arrQAPtBinning[i]        = 20.+1.0*(i-140);
      else arrQAPtBinning[i]                  = maxQAPt;
    }
    nBinsClusterPt            = 274;
    minClusterPt              = 0;
    maxClusterPt              = 100;
    for(Int_t i=0; i<nBinsClusterPt+1;i++){
      if (i < 1) arrClusPtBinning[i]          = 0.3*i;
      else if(i<98) arrClusPtBinning[i]       = 0.3+0.1*(i-1);
      else if(i<128) arrClusPtBinning[i]      = 10.+0.2*(i-98);
      else if(i<184) arrClusPtBinning[i]      = 16.+0.25*(i-128);
      else if(i<224) arrClusPtBinning[i]      = 30.+0.5*(i-184);
      else if(i<274) arrClusPtBinning[i]      = 50.+1.0*(i-224);
      else arrClusPtBinning[i]                = maxClusterPt;
    }
  // Set special pt binning for pp 13TeV and 5TeV
  } else if ( ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::k8TeV  ||
      ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::kpPb8TeV  ){
    nBinsPt                   = 139;
    minPt                     = 0;
    maxPt                     = 20;
    for(Int_t i=0; i<nBinsPt+1;i++){
      if (i < 100) arrPtBinning[i]            = 0.10*i;
      else if(i<140) arrPtBinning[i]          = 10.+0.25*(i-100);
      else arrPtBinning[i]                    = maxPt;
    }
    nBinsQAPt                 = 140;
    maxQAPt                   = 20;
    for(Int_t i=0; i<nBinsQAPt+1;i++){
      if(i<100) arrQAPtBinning[i]             = 0.1*i;
      else if(i<140) arrQAPtBinning[i]        = 10.+0.25*(i-100);
      else arrQAPtBinning[i]                  = maxQAPt;
    }
    nBinsClusterPt            = 274;
    minClusterPt              = 0;
    maxClusterPt              = 100;
    for(Int_t i=0; i<nBinsClusterPt+1;i++){
      if (i < 1) arrClusPtBinning[i]          = 0.3*i;
      else if(i<98) arrClusPtBinning[i]       = 0.3+0.1*(i-1);
      else if(i<128) arrClusPtBinning[i]      = 10.+0.2*(i-98);
      else if(i<184) arrClusPtBinning[i]      = 16.+0.25*(i-128);
      else if(i<224) arrClusPtBinning[i]      = 30.+0.5*(i-184);
      else if(i<274) arrClusPtBinning[i]      = 50.+1.0*(i-224);
      else arrClusPtBinning[i]                = maxClusterPt;
    }
  // Set special pt binning for XeXe 5.44TeV
  } else if (((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::kXeXe5440GeV ){
//    binWidthPt                = 0.1;
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
  // Set special pt binning for PbPb 5TeV
  } else if (((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::kPbPb5TeV ){

    nBinsPt     = 38;
    minPt       = 0.0;
    maxPt       = 25.0;
    for(Int_t i=0; i<=nBinsPt;i++){
      if(i<=20)           arrPtBinning[i]  = 0.0   + 0.2*i;             // 0.2 GeV bin width until 4 GeV
      else if(i<=30)      arrPtBinning[i]  = 4.0   + 0.5*(i-20);        // 0.5 GeV                 9 GeV
      else if(i<=36)      arrPtBinning[i]  = 9.0   + 1.0*(i-30);        // 1.0 GeV                 15 GeV
      else /*i<=nBinsPt*/ arrPtBinning[i]  = 15.0  + 5.0*(i-36);        // 5.0 GeV                 25 GeV
    }

    nBinsQAPt = 38;
    maxQAPt   = 25.0;
    for(Int_t i=0; i<=nBinsQAPt;i++){
      if(i<=20)           arrQAPtBinning[i]  = 0.0   + 0.2*i;             // 0.2 GeV bin width until 4 GeV
      else if(i<=30)      arrQAPtBinning[i]  = 4.0   + 0.5*(i-20);        // 0.5 GeV                 9 GeV
      else if(i<=36)      arrQAPtBinning[i]  = 9.0   + 1.0*(i-30);        // 1.0 GeV                 15 GeV
      else /*i<=nBinsPt*/ arrQAPtBinning[i]  = 15.0  + 5.0*(i-36);        // 5.0 GeV                 25 GeV
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
  //----------------------------------------------------------------------------------------------------------------
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
      else if(i<175) arrQAPtBinning[i]        = 20.+1*(i-170);
      else arrQAPtBinning[i]                  = maxQAPt;
    }
  }
  Double_t* arrLogBinning = new Double_t[51]{1.00e-03, 1.20e-03, 1.45e-03, 1.74e-03, 2.01e-03, 2.51e-03, 3.02e-03, 3.63e-03, 4.37e-03, 5.25e-03, 6.31e-03, 7.60e-03,
                                     9.12e-03, 1.01e-02, 1.32e-02, 1.58e-02, 1.91e-02, 2.29e-02, 2.75e-02, 3.31e-02, 3.98e-02, 4.79e-02, 5.75e-02, 6.91e-02,
                                     8.32e-02, 1.00e-01, 1.20e-01, 1.45e-01, 1.74e-01, 2.09e-01, 2.51e-01, 3.02e-01, 3.63e-01, 4.37e-01, 5.25e-01, 6.31e-01,
                                     7.59e-01, 9.12e-01, 1.10e+00, 1.32e+00, 1.58e+00, 1.91e+00, 2.29e+00, 2.75e+00, 3.31e+00, 3.98e+00, 4.79e+00, 5.75e+00,
                                     6.92e+00, 8.32e+00, 1.00e+01};


  if (fIsMC == 2){
    fDoPhotonQA           = 0;
    fDoTHnSparse          = kFALSE;
  } else if (fIsMC == 3){
    fDoTHnSparse          = kFALSE;
  }
  // Create histograms
  if(fOutputContainer != NULL){
    delete fOutputContainer;
    fOutputContainer        = NULL;
  }
  if(fOutputContainer == NULL){
    fOutputContainer        = new TList();
    fOutputContainer->SetOwner(kTRUE);
  }

  // Array of current cut's gammas
  fGammaCandidates          = new TList();

  fCutFolder                = new TList*[fnCuts];
  fESDList                  = new TList*[fnCuts];
  if(fDoTHnSparse){
    fBackList               = new TList*[fnCuts];
    fMotherList             = new TList*[fnCuts];
  }

  // event histos:
  fHistoNEvents             = new TH1F*[fnCuts];
  if (fIsMC > 1){
    fHistoNEventsWOWeight   = new TH1F*[fnCuts];
  }
  if (fIsMC == 2){
    fProfileJetJetXSection  = new TProfile*[fnCuts];
    fhJetJetNTrials         = new TH1F*[fnCuts];
  }
  fHistoNGoodESDTracks      = new TH1F*[fnCuts];
  fHistoVertexZ             = new TH1F*[fnCuts];
  fHistoEventSphericity     = new TH1F*[fnCuts];
  if(fDoPlotVsCentrality){
    fHistoCentrality        = new TH1F*[fnCuts];
    fHistoCentralityVsPrimaryTracks = new TH2F*[fnCuts];
  }
  if(fDoCentralityFlat > 0){
    fWeightCentrality            = new Double_t[fnCuts];
    fHistoNEventsWeighted        = new TH1F*[fnCuts];
    fHistoNGoodESDTracksWeighted = new TH1F*[fnCuts];
    fHistoVertexZWeighted        = new TH1F*[fnCuts];
    fHistoCentralityFlattened    = new TH1F*[fnCuts];
  }

  if(!fDoLightOutput){
    fHistoNGammaCandidates                     = new TH1F*[fnCuts];
    if (fIsMC < 2){
      fHistoNGoodESDTracksVsNGammaCandidates   = new TH2F*[fnCuts];
      fHistoSPDClusterTrackletBackground       = new TH2F*[fnCuts];
    }
    fHistoV0MultVsNumberTPCoutTracks           = new TH2F*[fnCuts];
    fHistoNV0Tracks                            = new TH1F*[fnCuts];
  }

  if(fEnableBDT){
    fHistoBDToutput = new TH1F*[fnCuts];
    fHistoBDToutputPt = new TH1F*[fnCuts];
    if (fIsMC > 0){
      fHistoBDToutputMCTrue = new TH1F*[fnCuts];
    }
  }

  fHistoEtaShift                             = new TProfile*[fnCuts];

  // gamma histos:
  fHistoConvGammaPt                          = new TH1F*[fnCuts];
  if (fDoPhotonQA > 0 && fIsMC < 2 ){
    fHistoConvGammaPsiPairPt    = new TH2F*[fnCuts];
    fHistoConvGammaR            = new TH1F*[fnCuts];
    fHistoConvGammaEta          = new TH1F*[fnCuts];
    fHistoConvGammaPhi          = new TH1F*[fnCuts];
    if ((fDoPhotonQA == 4)||(fDoPhotonQA == 5)){
      fHistoConvGammaInvMass      = new TH1F*[fnCuts];
      fHistoConvGammaInvMassReco  = new TH1F*[fnCuts];
    }
  }
  if ((fDoPhotonQA == 2)||(fDoPhotonQA == 5)){
    tESDConvGammaPtDcazCat      = new TTree*[fnCuts];
  }

  // meson histos:
  if(fDoMesonAnalysis){
    fHistoMotherInvMassPt           = new TH2F*[fnCuts];
    fHistoMotherBackInvMassPt       = new TH2F*[fnCuts];
    fHistoMotherInvMassEalpha       = new TH2F*[fnCuts];
    if(fDoMesonQA > 0){
      fHistoMotherPi0PtY            = new TH2F*[fnCuts];
      fHistoMotherEtaPtY            = new TH2F*[fnCuts];
      fHistoMotherPi0PtAlpha        = new TH2F*[fnCuts];
      fHistoMotherEtaPtAlpha        = new TH2F*[fnCuts];
      fHistoMotherPi0PtOpenAngle    = new TH2F*[fnCuts];
      fHistoMotherEtaPtOpenAngle    = new TH2F*[fnCuts];
    }
    if(fDoMesonQA == 2){
      tESDMesonsInvMassPtDcazMinDcazMaxFlag = new TTree*[fnCuts];
    }
    if(fDoMesonQA == 3){
      sPtRDeltaROpenAngle       = new THnSparseF*[fnCuts];
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
  if (fEnableClusterCutsForTrigger){
    fHistoCaloGammaPt           = new TH1F*[fnCuts];
    fHistoCaloGammaE            = new TH1F*[fnCuts];
  }


  for(Int_t iCut = 0; iCut<fnCuts;iCut++){

    TString cutstringEvent      = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
    TString cutstringPhoton     = ((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutNumber();
    TString cutstringMeson      = "NoMesonCut";
    if(fDoMesonAnalysis)
      cutstringMeson            = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();

    fCutFolder[iCut]            = new TList();
    fCutFolder[iCut]->SetName(Form("Cut Number %s_%s_%s",cutstringEvent.Data() ,cutstringPhoton.Data(),cutstringMeson.Data()));
    fCutFolder[iCut]->SetOwner(kTRUE);
    fOutputContainer->Add(fCutFolder[iCut]);
    fESDList[iCut]              = new TList();
    fESDList[iCut]->SetName(Form("%s_%s_%s ESD histograms",cutstringEvent.Data() ,cutstringPhoton.Data(),cutstringMeson.Data()));
    fESDList[iCut]->SetOwner(kTRUE);
    fCutFolder[iCut]->Add(fESDList[iCut]);

    if(fDoCentralityFlat > 0)
      fHistoNEvents[iCut]            = new TH1F("NEventsUnweighted", "NEventsUnweighted", 14, -0.5, 13.5);
    else
      fHistoNEvents[iCut]            = new TH1F("NEvents", "NEvents", 14, -0.5, 13.5);
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(3,"Miss. MC or inc. ev.");
    if (((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger() > 1 ){
      TString TriggerNames      = "Not Trigger: ";
      TriggerNames              = TriggerNames+ ( (AliConvEventCuts*)fEventCutArray->At(iCut))->GetSpecialTriggerName();
      fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(4,TriggerNames.Data());
    } else {
      fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(4,"Trigger");
    }
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(5,"Vertex Z");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(7,"Pile-Up");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(8,"no SDD");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(9,"no V0AND");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(10,"EMCAL problem");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(11,"rejectedForJetJetMC");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(14,"Pileup V0M-TPCout Tracks");
    fESDList[iCut]->Add(fHistoNEvents[iCut]);
    if (fIsMC > 1){
      fHistoNEventsWOWeight[iCut]    = new TH1F("NEventsWOWeight", "NEventsWOWeight", 14, -0.5, 13.5);
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(3,"Miss. MC or inc. ev.");
      if (((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger() > 1 ){
        TString TriggerNames    = "Not Trigger: ";
        TriggerNames            = TriggerNames+ ( (AliConvEventCuts*)fEventCutArray->At(iCut))->GetSpecialTriggerName();
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
      fProfileJetJetXSection[iCut]  = new TProfile("XSection", "XSection", 1, -0.5, 0.5);
      fESDList[iCut]->Add(fProfileJetJetXSection[iCut]);
      fhJetJetNTrials[iCut]         = new TH1F("NTrials", "#sum{NTrials}", 1, 0, 1);
      fhJetJetNTrials[iCut]->GetXaxis()->SetBinLabel(1,"#sum{NTrials}");
      fESDList[iCut]->Add(fhJetJetNTrials[iCut]);
    }
    if(fDoCentralityFlat > 0){
      fHistoNEventsWeighted[iCut]        = new TH1F("NEvents", "NEvents", 14, -0.5, 13.5);//weighted histogram!!
      fHistoNEventsWeighted[iCut]->Sumw2();
      fHistoNEventsWeighted[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
      fHistoNEventsWeighted[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
      fHistoNEventsWeighted[iCut]->GetXaxis()->SetBinLabel(3,"Miss. MC or inc. ev.");
      if (((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger() > 1 ){
        TString TriggerNames        = "Not Trigger: ";
        TriggerNames                = TriggerNames+ ( (AliConvEventCuts*)fEventCutArray->At(iCut))->GetSpecialTriggerName();
        fHistoNEventsWeighted[iCut]->GetXaxis()->SetBinLabel(4,TriggerNames.Data());
      } else {
        fHistoNEventsWeighted[iCut]->GetXaxis()->SetBinLabel(4,"Trigger");
      }
      fHistoNEventsWeighted[iCut]->GetXaxis()->SetBinLabel(5,"Vertex Z");
      fHistoNEventsWeighted[iCut]->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
      fHistoNEventsWeighted[iCut]->GetXaxis()->SetBinLabel(7,"Pile-Up");
      fHistoNEventsWeighted[iCut]->GetXaxis()->SetBinLabel(8,"no SDD");
      fHistoNEventsWeighted[iCut]->GetXaxis()->SetBinLabel(9,"no V0AND");
      fHistoNEventsWeighted[iCut]->GetXaxis()->SetBinLabel(10,"EMCAL problem");
      fHistoNEventsWeighted[iCut]->GetXaxis()->SetBinLabel(11,"rejectedForJetJetMC");
      fHistoNEventsWeighted[iCut]->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
      fHistoNEventsWeighted[iCut]->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
      fHistoNEventsWeighted[iCut]->GetXaxis()->SetBinLabel(14,"Pileup V0M-TPCout Tracks");
      fESDList[iCut]->Add(fHistoNEventsWeighted[iCut]);
    }

    if(fDoCentralityFlat > 0 && fIsHeavyIon == 1){
      fHistoNGoodESDTracks[iCut]         = new TH1F("GoodESDTracksUnweighted", "GoodESDTracksUnweighted", nTracks, -0.5, nTracks-0.5);
      fESDList[iCut]->Add(fHistoNGoodESDTracks[iCut]);
      fHistoVertexZ[iCut]                = new TH1F("VertexZUnweighted", "VertexZUnweighted", 200, -10, 10);
      fESDList[iCut]->Add(fHistoVertexZ[iCut]);

      fHistoNGoodESDTracksWeighted[iCut] = new TH1F("GoodESDTracks", "GoodESDTracks", nTracks, -0.5, nTracks-0.5); //weighted histogram!!
      fHistoNGoodESDTracksWeighted[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoNGoodESDTracksWeighted[iCut]);
      fHistoVertexZWeighted[iCut]        = new TH1F("VertexZ", "VertexZ", 200, -10, 10);
      fHistoVertexZWeighted[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoVertexZWeighted[iCut]);
    } else {
      fHistoNGoodESDTracks[iCut]       = new TH1F("GoodESDTracks", "GoodESDTracks", nTracks, -0.5, nTracks-0.5);
      fESDList[iCut]->Add(fHistoNGoodESDTracks[iCut]);
      fHistoVertexZ[iCut]                = new TH1F("VertexZ", "VertexZ", 200, -10, 10);
      fESDList[iCut]->Add(fHistoVertexZ[iCut]);
    }

    if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetUseSphericity()!=0){
      fHistoEventSphericity[iCut]     = new TH1F("EventSphericity", "EventSphericity", 100, 0, 1);
      fHistoEventSphericity[iCut]->GetXaxis()->SetTitle("S");
      fESDList[iCut]->Add(fHistoEventSphericity[iCut]);
      fV0Reader->SetCalcSphericity(kTRUE);
    }

    if(fDoPlotVsCentrality){
      fHistoCentrality[iCut]                 = new TH1F("Centrality", "Centrality", 100, -0.5, 99.5);
      fESDList[iCut]->Add(fHistoCentrality[iCut]);
      fHistoCentralityVsPrimaryTracks[iCut]  = new TH2F("Centrality vs Primary Tracks", "Centrality vs Primary Tracks ", 100, -0.5, 99.5, nTracks, -0.5, nTracks-0.5);
      if(fDoCentralityFlat > 0 || fIsMC > 1) fHistoCentralityVsPrimaryTracks[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoCentralityVsPrimaryTracks[iCut]);
    }
    if(fDoCentralityFlat > 0){
      fHistoCentralityFlattened[iCut]    = new TH1F("CentralityFlattened", "CentralityFlattened", 100, -0.5, 99.5);
      fHistoCentralityFlattened[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoCentralityFlattened[iCut]);
    }

    if(!fDoLightOutput){
      fHistoNGammaCandidates[iCut]       = new TH1F("GammaCandidates", "GammaCandidates", nGammaCand, -0.5, nGammaCand-0.5);
      if(fDoCentralityFlat > 0 || fIsMC > 1) fHistoNGammaCandidates[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoNGammaCandidates[iCut]);
      if (fIsMC < 2){
	fHistoNGoodESDTracksVsNGammaCandidates[iCut]   = new TH2F("GoodESDTracksVsGammaCandidates", "GoodESDTracksVsGammaCandidates", nTracks, -0.5, nTracks-0.5, nGammaCand, -0.5, nGammaCand-0.5);
	if(fDoCentralityFlat > 0) fHistoNGoodESDTracksVsNGammaCandidates[iCut]->Sumw2();
	fESDList[iCut]->Add(fHistoNGoodESDTracksVsNGammaCandidates[iCut]);
	fHistoSPDClusterTrackletBackground[iCut]    = new TH2F("SPD tracklets vs SPD clusters", "SPD tracklets vs SPD clusters", nBinsTrklCls, 0, nSPDTracklets, nBinsTrklCls, 0, nSPDClusters);
	if(fDoCentralityFlat > 0) fHistoSPDClusterTrackletBackground[iCut]->Sumw2();
	fESDList[iCut]->Add(fHistoSPDClusterTrackletBackground[iCut]);
      }

      if(fIsHeavyIon == 1)
	fHistoV0MultVsNumberTPCoutTracks[iCut]    = new TH2F("V0Mult vs TPCout Tracks", "V0Mult vs TPCout Tracks", 500, 0, 15000, 500, 0, 40000);
      else if(fIsHeavyIon == 2)
	fHistoV0MultVsNumberTPCoutTracks[iCut]    = new TH2F("V0Mult vs TPCout Tracks", "V0Mult vs TPCout Tracks", 500, 0, 1000, 500, 0, 2500);
      else
	fHistoV0MultVsNumberTPCoutTracks[iCut]    = new TH2F("V0Mult vs TPCout Tracks", "V0Mult vs TPCout Tracks", 200, 0, 400, 500, 0, 1500);
      fESDList[iCut]->Add(fHistoV0MultVsNumberTPCoutTracks[iCut]);

      if(fIsHeavyIon == 1)
	fHistoNV0Tracks[iCut]            = new TH1F("V0 Multiplicity", "V0 Multiplicity", 20000, 0, 40000);
      else if(fIsHeavyIon == 2)
	fHistoNV0Tracks[iCut]            = new TH1F("V0 Multiplicity", "V0 Multiplicity", 2500, 0, 2500);
      else
	fHistoNV0Tracks[iCut]            = new TH1F("V0 Multiplicity", "V0 Multiplicity", 1500, 0, 1500);
      if(fDoCentralityFlat > 0 || fIsMC > 1) fHistoNV0Tracks[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoNV0Tracks[iCut]);
    }

    if(fEnableBDT){
      fHistoBDToutput[iCut] = new TH1F("BDT_output", "BDT_output", 200, -1, 1);
      fESDList[iCut]->Add(fHistoBDToutput[iCut]);

      fHistoBDToutputPt[iCut] = new TH1F("BDT_ConvGamma_Pt", "BDT_ESD_ConvGamma_Pt", nBinsPt, arrPtBinning);
      fESDList[iCut]->Add(fHistoBDToutputPt[iCut]);
    }

    fHistoEtaShift[iCut]               = new TProfile("Eta Shift", "Eta Shift", 1, -0.5, 0.5);
    fESDList[iCut]->Add(fHistoEtaShift[iCut]);

    fHistoConvGammaPt[iCut]         = new TH1F("ESD_ConvGamma_Pt", "ESD_ConvGamma_Pt", nBinsPt, arrPtBinning);
    if(fDoCentralityFlat > 0 ) fHistoConvGammaPt[iCut]->Sumw2();
    fESDList[iCut]->Add(fHistoConvGammaPt[iCut]);

    if (  (fIsMC > 1) || (fIsMC>0 && fDoMaterialBudgetWeightingOfGammasForTrueMesons) ) {
      fHistoConvGammaPt[iCut]->Sumw2();
    }

    if (fIsMC > 1){
      fHistoNEvents[iCut]->Sumw2();
      fHistoNGoodESDTracks[iCut]->Sumw2();
      fHistoVertexZ[iCut]->Sumw2();
      fHistoEtaShift[iCut]->Sumw2();
      if(fDoPlotVsCentrality){
        fHistoCentrality[iCut]->Sumw2();
      }
      if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetUseSphericity()!=0){
        fHistoEventSphericity[iCut]->Sumw2();
      }
    }

    if(fEnableClusterCutsForTrigger){
      fHistoCaloGammaPt[iCut]         = new TH1F("ClusGamma_Pt", "ClusGamma_Pt", nBinsClusterPt, arrClusPtBinning);
      fHistoCaloGammaPt[iCut]->SetXTitle("p_{T,clus}(GeV/c)");
      fESDList[iCut]->Add(fHistoCaloGammaPt[iCut]);
      if (fIsMC > 1) fHistoCaloGammaPt[iCut]->Sumw2();
      fHistoCaloGammaE[iCut]          = new TH1F("ClusGamma_E", "ClusGamma_E", nBinsClusterPt, arrClusPtBinning);
      fHistoCaloGammaE[iCut]->SetXTitle("E_{clus}(GeV)");
      fESDList[iCut]->Add(fHistoCaloGammaE[iCut]);
      if (fIsMC > 1) fHistoCaloGammaE[iCut]->Sumw2();
    }

    if ((fDoPhotonQA == 2)||(fDoPhotonQA == 5)){
      if(fDoMaterialBudgetWeightingOfGammasForTrueMesons && fIsMC > 0){
	tESDConvGammaPtDcazCat[iCut]= new TTree(Form("%s_%s_%s MBW Photon DCA tree", cutstringEvent.Data(), cutstringPhoton.Data(), cutstringMeson.Data()), "ESD_ConvGamma_Pt_Dcaz_R_Eta_Cat");
       }else {
	tESDConvGammaPtDcazCat[iCut]= new TTree(Form("%s_%s_%s Photon DCA tree", cutstringEvent.Data(), cutstringPhoton.Data(), cutstringMeson.Data()), "ESD_ConvGamma_Pt_Dcaz_R_Eta_Cat");
      }
      tESDConvGammaPtDcazCat[iCut]->Branch("Pt",&fPtGamma,"fPtGamma/F");
      tESDConvGammaPtDcazCat[iCut]->Branch("DcaZPhoton",&fDCAzPhoton,"fDCAzPhoton/F");
      tESDConvGammaPtDcazCat[iCut]->Branch("cat",&iCatPhoton,"iCatPhoton/b");
      if(fIsMC>0){
        tESDConvGammaPtDcazCat[iCut]->Branch("photonMCInfo",&iPhotonMCInfo,"iPhotonMCInfo/b");
      }
      if (fIsMC > 1){
        tESDConvGammaPtDcazCat[iCut]->Branch("weightEvent",&fWeightJetJetMC,"fWeightJetJetMC/F");
      }
    }

    if(fDoPhotonQA > 0 && fIsMC < 2){

      fHistoConvGammaPsiPairPt[iCut]= new TH2F("ESD_ConvGamma_PsiPair_Pt", "ESD_ConvGamma_PsiPair_Pt", 500, 0, 5, nBinsQAPt, arrQAPtBinning);
      if(fDoCentralityFlat > 0) fHistoConvGammaPsiPairPt[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoConvGammaPsiPairPt[iCut]);

      fHistoConvGammaR[iCut]        = new TH1F("ESD_ConvGamma_R", "ESD_ConvGamma_R", 800, 0, 200);
      if(fDoCentralityFlat > 0) fHistoConvGammaR[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoConvGammaR[iCut]);
      fHistoConvGammaEta[iCut]      = new TH1F("ESD_ConvGamma_Eta", "ESD_ConvGamma_Eta", 1000, -2, 2);
      if(fDoCentralityFlat > 0) fHistoConvGammaEta[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoConvGammaEta[iCut]);
      fHistoConvGammaPhi[iCut]      = new TH1F("ESD_ConvGamma_Phi", "ESD_ConvGamma_Phi", 360, 0, 2*TMath::Pi());
      if(fDoCentralityFlat > 0) fHistoConvGammaPhi[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoConvGammaPhi[iCut]);
      if ((fDoPhotonQA == 4)||(fDoPhotonQA == 5)){
        fHistoConvGammaInvMass[iCut]      = new TH1F("ESD_ConvGamma_InvMass", "ESD_ConvGamma_InvMass", 500, 0., 5.);
        if(fDoCentralityFlat > 0) fHistoConvGammaInvMass[iCut]->Sumw2();
        fESDList[iCut]->Add(fHistoConvGammaInvMass[iCut]);
        fHistoConvGammaInvMassReco[iCut]      = new TH1F("ESD_ConvGamma_InvMassReco", "ESD_ConvGamma_InvMassReco", 500, 0., 5.);
        if(fDoCentralityFlat > 0) fHistoConvGammaInvMassReco[iCut]->Sumw2();
        fESDList[iCut]->Add(fHistoConvGammaInvMassReco[iCut]);
      }
    }

    if(fDoMesonAnalysis){
      fHistoMotherInvMassPt[iCut]   = new TH2F("ESD_Mother_InvMass_Pt", "ESD_Mother_InvMass_Pt", 800, 0, 0.8, nBinsPt, arrPtBinning);
      fESDList[iCut]->Add(fHistoMotherInvMassPt[iCut]);
      fHistoMotherBackInvMassPt[iCut]   = new TH2F("ESD_Background_InvMass_Pt", "ESD_Background_InvMass_Pt", 800, 0, 0.8, nBinsPt, arrPtBinning);
      fESDList[iCut]->Add(fHistoMotherBackInvMassPt[iCut]);
      fHistoMotherInvMassEalpha[iCut]   = new TH2F("ESD_Mother_InvMass_vs_E_alpha", "ESD_Mother_InvMass_vs_E_alpha", 800, 0, 0.8, nBinsPt, arrPtBinning);
      fESDList[iCut]->Add(fHistoMotherInvMassEalpha[iCut]);

      if (fIsMC > 1 || fDoCentralityFlat > 0){
        fHistoMotherInvMassPt[iCut]->Sumw2();
        fHistoMotherBackInvMassPt[iCut]->Sumw2();
        fHistoMotherInvMassEalpha[iCut]->Sumw2();
      }

      if(fDoMesonQA == 2){
	if(fDoMaterialBudgetWeightingOfGammasForTrueMesons && fIsMC > 0){
	  tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut] = new TTree(Form("%s_%s_%s MBW Meson DCA tree",cutstringEvent.Data() ,cutstringPhoton.Data(),cutstringMeson.Data()), "ESD_Mesons_InvMass_Pt_DcazMin_DcazMax_Flag");
	}else {
	  tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut] = new TTree(Form("%s_%s_%s Meson DCA tree",cutstringEvent.Data() ,cutstringPhoton.Data(),cutstringMeson.Data()), "ESD_Mesons_InvMass_Pt_DcazMin_DcazMax_Flag");
	}

        tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut]->Branch("InvMass",&fInvMass,"fInvMass/F");
        tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut]->Branch("Pt",&fPt,"fPt/F");
        tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut]->Branch("DcaZMin",&fDCAzGammaMin,"fDCAzGammaMin/F");
        tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut]->Branch("DcaZMax",&fDCAzGammaMax,"fDCAzGammaMax/F");
        tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut]->Branch("kind",&iFlag,"iFlag/b");
        if(fIsMC>0){
          tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut]->Branch("mesonMCInfo",&iMesonMCInfo,"iMesonMCInfo/b");
        }
      }
      if(fDoMesonQA > 0 ){
        if (fIsMC < 2){
          fHistoMotherPi0PtY[iCut]  = new TH2F("ESD_MotherPi0_Pt_Y", "ESD_MotherPi0_Pt_Y", nBinsQAPt, arrQAPtBinning, 150, -1.5, 1.5);
          fESDList[iCut]->Add(fHistoMotherPi0PtY[iCut]);
          fHistoMotherEtaPtY[iCut]  = new TH2F("ESD_MotherEta_Pt_Y", "ESD_MotherEta_Pt_Y", nBinsQAPt, arrQAPtBinning, 150, -1.5, 1.5);
          fESDList[iCut]->Add(fHistoMotherEtaPtY[iCut]);
          fHistoMotherPi0PtOpenAngle[iCut]  = new TH2F("ESD_MotherPi0_Pt_OpenAngle", "ESD_MotherPi0_Pt_OpenAngle", nBinsQAPt, arrQAPtBinning, 100, 0, TMath::Pi());
          fESDList[iCut]->Add(fHistoMotherPi0PtOpenAngle[iCut]);
          fHistoMotherEtaPtOpenAngle[iCut]  = new TH2F("ESD_MotherEta_Pt_OpenAngle", "ESD_MotherEta_Pt_OpenAngle", nBinsQAPt, arrQAPtBinning, 100, 0, TMath::Pi());
          fESDList[iCut]->Add(fHistoMotherEtaPtOpenAngle[iCut]);
        }
        fHistoMotherPi0PtAlpha[iCut] = new TH2F("ESD_MotherPi0_Pt_Alpha", "ESD_MotherPi0_Pt_Alpha", nBinsQAPt, arrQAPtBinning, 100, 0, 1);
        fESDList[iCut]->Add(fHistoMotherPi0PtAlpha[iCut]);
        fHistoMotherEtaPtAlpha[iCut] = new TH2F("ESD_MotherEta_Pt_Alpha", "ESD_MotherEta_Pt_Alpha", nBinsQAPt, arrQAPtBinning, 100, 0, 1);
        fESDList[iCut]->Add(fHistoMotherEtaPtAlpha[iCut]);
        if(fIsMC > 1){
          fHistoMotherPi0PtAlpha[iCut]->Sumw2();
          fHistoMotherEtaPtAlpha[iCut]->Sumw2();
        }

      }
      if(fDoMesonQA == 3){
        sPtRDeltaROpenAngle[iCut]   = new THnSparseF("PhotonPair_Pt_R_DeltaR_OpenAngle", "PhotonPair_Pt_R_DeltaR_OpenAngle", nDim2, nBins2, xMin2, xMax2);
        fESDList[iCut]->Add(sPtRDeltaROpenAngle[iCut]);
      }
    }
    if(fDoJetAnalysis){

      fJetHistograms[iCut] = new TList();
      fJetHistograms[iCut]->SetOwner(kTRUE);
      fJetHistograms[iCut]->SetName(Form("%s_%s_%s Jet histograms",cutstringEvent.Data() ,cutstringPhoton.Data(),cutstringMeson.Data()));

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
        fHistoJetMotherInvMassPt[iCut] = new TH2F("ESD_Pi0Jet_Mother_InvMass_Pt", "ESD_Pi0Jet_Mother_InvMass_Pt", 800, 0, 0.8, nBinsPt, arrPtBinning);
        fJetHistograms[iCut]->Add(fHistoJetMotherInvMassPt[iCut]);
        fHistoEtaPhiJetPi0Cand[iCut] = new TH2F("Eta_Phi_Distr_Pi0Jet", "Eta_Phi_Distr_Pi0Jet", 20, 0, M_PI, 20, -1, 1);
        fJetHistograms[iCut]->Add(fHistoEtaPhiJetPi0Cand[iCut]);
        fHistoRJetPi0Cand[iCut] = new TH2F("ESD_RPi0Jet_Pt", "ESD_RPi0Jet_Pt", 35, 0, 3.5, nBinsPt, arrPtBinning);
        fJetHistograms[iCut]->Add(fHistoRJetPi0Cand[iCut]);
        fHistoPi0InJetMotherInvMassPt[iCut] = new TH2F("ESD_Pi0inJet_Mother_InvMass_Pt", "ESD_Pi0inJet_Mother_InvMass_Pt", 800, 0, 0.8, nBinsPt, arrPtBinning);
        fJetHistograms[iCut]->Add(fHistoPi0InJetMotherInvMassPt[iCut]);
        fHistoMotherBackJetInvMassPt[iCut] = new TH2F("ESD_Jet_Background_InvMass_Pt", "ESD_Jet_Background_InvMass_Pt", 800, 0, 0.8, nBinsPt, arrPtBinning);
        fJetHistograms[iCut]->Add(fHistoMotherBackJetInvMassPt[iCut]);
        fHistoEtaPhiJetWithPi0Cand[iCut] = new TH2F("Eta_Phi_Distr_Pi0inJet", "Eta_Phi_Distr_Pi0inJet", 15, 0, 0.4, 15, -0.4, 0.4);
        fJetHistograms[iCut]->Add(fHistoEtaPhiJetWithPi0Cand[iCut]);
        fHistoDoubleCounting[iCut] = new TH1F("Double_Counting_Mesons_Jets", "Double_Counting_Mesons_Jets", 6, 0, 6);
        fJetHistograms[iCut]->Add(fHistoDoubleCounting[iCut]);
        fHistoJetFragmFunc[iCut] = new TH2F("ESD_Pi0inJetPt_FragmentationFunc", "ESD_Pi0inJetPt_FragmentationFunc", 50, arrLogBinning, 150, 0., 150.);
        fJetHistograms[iCut]->Add(fHistoJetFragmFunc[iCut]);
        fHistoJetFragmFuncZInvMass[iCut] = new TH2F("ESD_Pi0inJetPt_Fragm_Z_InvMass", "ESD_Pi0inJetPt_Fragm_Z_InvMass", 800, 0, 0.8, 50, arrLogBinning);
        fJetHistograms[iCut]->Add(fHistoJetFragmFuncZInvMass[iCut]);
      }
    }


  }
  if(fDoMesonAnalysis){
    InitBack(); // Init Background Handler
  }

  if(fIsMC>0){
    // MC Histogramms
    fMCList                       = new TList*[fnCuts];
    // True Histogramms
    fTrueList                     = new TList*[fnCuts];
    // Selected Header List
    fHeaderNameList               = new TList*[fnCuts];
    fHistoMCHeaders                    = new TH1I*[fnCuts];
    fHistoMCAllGammaPt                 = new TH1F*[fnCuts];
    fHistoMCAllSecondaryGammaPt        = new TH2F*[fnCuts];
    fHistoMCDecayGammaPi0Pt            = new TH1F*[fnCuts];
    fHistoMCDecayGammaRhoPt            = new TH1F*[fnCuts];
    fHistoMCDecayGammaEtaPt            = new TH1F*[fnCuts];
    fHistoMCDecayGammaOmegaPt          = new TH1F*[fnCuts];
    fHistoMCDecayGammaEtapPt           = new TH1F*[fnCuts];
    fHistoMCDecayGammaPhiPt            = new TH1F*[fnCuts];
    fHistoMCDecayGammaSigmaPt          = new TH1F*[fnCuts];
    fHistoMCConvGammaPt                = new TH1F*[fnCuts];
    fHistoMCSecondaryConvGammaPt       = new TH2F*[fnCuts];
    fHistoTrueConvGammaPt              = new TH1F*[fnCuts];
    fHistoDoubleCountTrueConvGammaRPt  = new TH2F*[fnCuts];
    fHistoMultipleCountTrueConvGamma   = new TH1F*[fnCuts];

    fHistoCombinatorialPt              = new TH2F*[fnCuts];
    if (fDoPhotonQA == 3){
      if(fIsHeavyIon == 1) fHistoCombinatorialMothersPt = new TH3F*[fnCuts];
      fHistoCombinatorialPtDeltaPhi_ek  = new TH2F*[fnCuts];
      fHistoCombinatorialPtDeltaPhi_ep  = new TH2F*[fnCuts];
      fHistoCombinatorialPtDeltaPhi_epi = new TH2F*[fnCuts];
      fHistoCombinatorialPtDeltaPhi_pik = new TH2F*[fnCuts];
      fHistoCombinatorialPtDeltaPhi_pip = new TH2F*[fnCuts];
    }
    fHistoTruePrimaryConvGammaPt                          = new TH1F*[fnCuts];
    fHistoTruePrimaryConvGammaESDPtMCPt                   = new TH2F*[fnCuts];
    fHistoTrueSecondaryConvGammaFromXFromK0sMCPtESDPt     = new TH2F*[fnCuts];
    fHistoTrueSecondaryConvGammaFromXFromK0lMCPtESDPt     = new TH2F*[fnCuts];
    fHistoTrueSecondaryConvGammaFromXFromLambdaMCPtESDPt  = new TH2F*[fnCuts];
    fHistoTrueSecondaryConvGammaPt                        = new TH2F*[fnCuts];
    fHistoTrueSecondaryConvGammaMCPt                      = new TH2F*[fnCuts];

    fHistoTrueDalitzPsiPairDeltaPhi     = new TH2F*[fnCuts];
    fHistoTrueGammaPsiPairDeltaPhi      = new TH2F*[fnCuts];

    if (fDoPhotonQA > 0 ){
      if (fIsMC < 2 ){
        fHistoMCConvGammaR              = new TH1F*[fnCuts];
        fHistoMCConvGammaEta            = new TH1F*[fnCuts];
        fHistoTrueConvGammaPsiPairPt    = new TH2F*[fnCuts];
        fHistoTrueConvGammaEta          = new TH1F*[fnCuts];
        fHistoTrueConvGammaR            = new TH1F*[fnCuts];
        fHistoTrueConvGammaRMC          = new TH1F*[fnCuts];
        if ((fDoPhotonQA == 4)||(fDoPhotonQA == 5)){
          fHistoTrueConvGammaInvMass      = new TH1F*[fnCuts];
          fHistoTrueConvGammaInvMassReco  = new TH1F*[fnCuts];
        }
      }
      fHistoTrueConvGammaPtMC           = new TH1F*[fnCuts];
    }

    if(fDoMesonAnalysis){
      fHistoMCPi0Pt                    = new TH1F*[fnCuts];
      fHistoMCPi0WOWeightPt            = new TH1F*[fnCuts];
      fHistoMCEtaPt                    = new TH1F*[fnCuts];
      fHistoMCEtaWOWeightPt            = new TH1F*[fnCuts];
      fHistoMCPi0WOWeightInAccPt       = new TH1F*[fnCuts];
      fHistoMCEtaWOWeightInAccPt       = new TH1F*[fnCuts];
      fHistoMCPi0InAccPt               = new TH1F*[fnCuts];
      fHistoMCEtaInAccPt               = new TH1F*[fnCuts];

      if(fIsMC > 1){
        fHistoMCPi0WOEvtWeightPt       = new TH1F*[fnCuts];
        fHistoMCEtaWOEvtWeightPt       = new TH1F*[fnCuts];
        fHistoMCPi0WOEvtWeightInAccPt  = new TH1F*[fnCuts];
        fHistoMCEtaWOEvtWeightInAccPt  = new TH1F*[fnCuts];
      }

      fHistoTrueMotherInvMassPt                      = new TH2F*[fnCuts];
      fHistoDoubleCountTruePi0InvMassPt              = new TH2F*[fnCuts];
      fHistoMultipleCountTruePi0                     = new TH1F*[fnCuts];
      fHistoDoubleCountTrueEtaInvMassPt              = new TH2F*[fnCuts];
      fHistoMultipleCountTrueEta                     = new TH1F*[fnCuts];
      fHistoTruePrimaryMotherInvMassPt               = new TH2F*[fnCuts];
      fHistoTruePrimaryMotherW0WeightingInvMassPt    = new TH2F*[fnCuts];
      pESDTruePrimaryMotherWeightsInvMassPt          = new TProfile2D*[fnCuts];
      fHistoTrueSecondaryMotherInvMassPt             = new TH2F*[fnCuts];
      fHistoTrueSecondaryMotherFromK0sInvMassPt      = new TH2F*[fnCuts];
      fHistoTrueSecondaryMotherFromK0lInvMassPt      = new TH2F*[fnCuts];
      fHistoTrueSecondaryMotherFromEtaInvMassPt      = new TH2F*[fnCuts];
      fHistoTrueSecondaryMotherFromLambdaInvMassPt   = new TH2F*[fnCuts];
      fHistoTrueMotherDalitzInvMassPt                = new TH2F*[fnCuts];

      fHistoMCPrimaryPtvsSource        = new TH2F*[fnCuts];
      fHistoMCSecPi0PtvsSource         = new TH2F*[fnCuts];
      fHistoMCSecPi0InAccPtvsSource    = new TH2F*[fnCuts];
      fHistoMCSecPi0Source             = new TH1F*[fnCuts];
      fHistoMCSecEtaPt                 = new TH1F*[fnCuts];
      fHistoMCSecEtaSource             = new TH1F*[fnCuts];

      if (fDoMesonQA > 0){
        fHistoMCPi0PtAlpha             = new TH2F*[fnCuts];
        fHistoMCEtaPtAlpha             = new TH2F*[fnCuts];
        if (fIsMC == 2){
          fHistoMCPi0PtJetPt           = new TH2F*[fnCuts];
          fHistoMCEtaPtJetPt           = new TH2F*[fnCuts];
        }

        if (fIsMC < 2){
          fHistoMCSecPi0RvsSource      = new TH2F*[fnCuts];
          fHistoMCPi0PtY               = new TH2F*[fnCuts];
          fHistoMCEtaPtY               = new TH2F*[fnCuts];
          fHistoTruePrimaryPi0MCPtResolPt         = new TH2F*[fnCuts];
          fHistoTruePrimaryEtaMCPtResolPt         = new TH2F*[fnCuts];
          fHistoTrueK0sWithPi0DaughterMCPt        = new TH1F*[fnCuts];
          fHistoTrueK0lWithPi0DaughterMCPt        = new TH1F*[fnCuts];
          fHistoTrueEtaWithPi0DaughterMCPt        = new TH1F*[fnCuts];
          fHistoTrueLambdaWithPi0DaughterMCPt     = new TH1F*[fnCuts];
          fHistoTrueBckGGInvMassPt                = new TH2F*[fnCuts];
          fHistoTrueBckContInvMassPt              = new TH2F*[fnCuts];
          fHistoTruePi0PtY          = new TH2F*[fnCuts];
          fHistoTrueEtaPtY          = new TH2F*[fnCuts];
          fHistoTruePi0PtOpenAngle  = new TH2F*[fnCuts];
          fHistoTrueEtaPtOpenAngle  = new TH2F*[fnCuts];
        }
        fHistoTruePi0PtAlpha        = new TH2F*[fnCuts];
        fHistoTrueEtaPtAlpha        = new TH2F*[fnCuts];

      }
    }
    if(fDoJetAnalysis && !fDoLightOutput) {
      fTrueJetHistograms                          = new TList*[fnCuts];
      fHistoTruevsRecJetPt                        = new TH2F*[fnCuts];
      fHistoTrueJetMotherInvMassPt                = new TH2F*[fnCuts];
      fHistoTrueInJetMotherInvMassPt              = new TH2F*[fnCuts];
      fHistoTruePrimaryJetInvMassPt               = new TH2F*[fnCuts];
      fHistoTruePrimaryinJetInvMassPt             = new TH2F*[fnCuts];
      fHistoTruePrimaryInJetInvMassTruePt         = new TH2F*[fnCuts];
      fHistoTrueDoubleCountingJet                 = new TH1F*[fnCuts];
      fHistoTrueJetFragmFunc                      = new TH2F*[fnCuts];
      fHistoTrueJetFragmFuncZInvMass              = new TH2F*[fnCuts];
      fHistoMCPi0JetInAccPt                       = new TH1F*[fnCuts];
      fHistoMCPi0inJetInAccPt                     = new TH1F*[fnCuts];
      fHistoMCPi0JetEventGenerated                = new TH1F*[fnCuts];
      fHistoMCPi0inJetGenerated                   = new TH1F*[fnCuts];
      fHistoTrueSecondaryFromK0sJetInvMassPt      = new TH2F*[fnCuts];
      fHistoTrueSecondaryFromK0sinJetInvMassPt    = new TH2F*[fnCuts];
      fHistoTrueSecondaryFromLambdaJetInvMassPt   = new TH2F*[fnCuts];
      fHistoTrueSecondaryFromLambdainJetInvMassPt = new TH2F*[fnCuts];
      fHistoTrueSecondaryFromK0lJetInvMassPt      = new TH2F*[fnCuts];
      fHistoTrueSecondaryFromK0linJetInvMassPt    = new TH2F*[fnCuts];
      fHistoTrueSecondaryInvJetMassPt             = new TH2F*[fnCuts];
      fHistoTrueSecondaryInvinJetMassPt           = new TH2F*[fnCuts];
      fHistoMCEtaJetInAccPt                       = new TH1F*[fnCuts];
      fHistoMCEtainJetInAccPt                     = new TH1F*[fnCuts];
      fHistoMCEtaJetEventGenerated                = new TH1F*[fnCuts];
      fHistoMCEtainJetGenerated                   = new TH1F*[fnCuts];
    }
     if(fDoJetQA){
        if(fDoLightOutput){
          fTrueJetHistograms                           = new TList*[fnCuts];
          fHistoTruevsRecJetPt                        = new TH2F*[fnCuts];
        }
        fHistoUnfoldingAsData                          = new TH2F*[fnCuts];
        fHistoUnfoldingMissed                          = new TH2F*[fnCuts];
        fHistoUnfoldingReject                          = new TH2F*[fnCuts];
        fHistoUnfoldingAsDataInvMassZ                  = new TH2F*[fnCuts];
        fHistoUnfoldingMissedInvMassZ                  = new TH2F*[fnCuts];
        fHistoUnfoldingRejectInvMassZ                  = new TH2F*[fnCuts];
      }

    for(Int_t iCut = 0; iCut<fnCuts;iCut++){
      TString cutstringEvent          = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
      TString cutstringPhoton         = ((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutNumber();
      TString cutstringMeson          = "NoMesonCut";
      if(fDoMesonAnalysis)
        cutstringMeson                = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();

      fMCList[iCut]                   = new TList();
      fMCList[iCut]->SetName(Form("%s_%s_%s MC histograms",cutstringEvent.Data() ,cutstringPhoton.Data(),cutstringMeson.Data()));
      fMCList[iCut]->SetOwner(kTRUE);
      fCutFolder[iCut]->Add(fMCList[iCut]);
      if (fIsMC < 2){
        fHistoMCHeaders[iCut]              = new TH1I("MC_Headers", "MC_Headers", 20, 0, 20);
        fMCList[iCut]->Add(fHistoMCHeaders[iCut]);
      }
      fHistoMCAllGammaPt[iCut]             = new TH1F("MC_AllGamma_Pt", "MC_AllGamma_Pt", nBinsPt, arrPtBinning);
      fMCList[iCut]->Add(fHistoMCAllGammaPt[iCut]);
      fHistoMCAllSecondaryGammaPt[iCut]    = new TH2F("MC_AllSecondaryGamma_Pt", "MC_AllSecondaryGamma_Pt", nBinsPt, arrPtBinning, 4, -0.5, 3.5);
      fHistoMCAllSecondaryGammaPt[iCut]->GetYaxis()->SetBinLabel(1,"K0s");
      fHistoMCAllSecondaryGammaPt[iCut]->GetYaxis()->SetBinLabel(2,"K0l");
      fHistoMCAllSecondaryGammaPt[iCut]->GetYaxis()->SetBinLabel(3,"Lambda");
      fHistoMCAllSecondaryGammaPt[iCut]->GetYaxis()->SetBinLabel(4,"rest");
      fMCList[iCut]->Add(fHistoMCAllSecondaryGammaPt[iCut]);
      fHistoMCDecayGammaPi0Pt[iCut]        = new TH1F("MC_DecayGammaPi0_Pt", "MC_DecayGammaPi0_Pt", nBinsPt, arrPtBinning);
      fMCList[iCut]->Add(fHistoMCDecayGammaPi0Pt[iCut]);
      fHistoMCDecayGammaRhoPt[iCut]        = new TH1F("MC_DecayGammaRho_Pt", "MC_DecayGammaRho_Pt", nBinsPt, arrPtBinning);
      fMCList[iCut]->Add(fHistoMCDecayGammaRhoPt[iCut]);
      fHistoMCDecayGammaEtaPt[iCut]        = new TH1F("MC_DecayGammaEta_Pt", "MC_DecayGammaEta_Pt", nBinsPt, arrPtBinning);
      fMCList[iCut]->Add(fHistoMCDecayGammaEtaPt[iCut]);
      fHistoMCDecayGammaOmegaPt[iCut]      = new TH1F("MC_DecayGammaOmega_Pt", "MC_DecayGammaOmmega_Pt", nBinsPt, arrPtBinning);
      fMCList[iCut]->Add(fHistoMCDecayGammaOmegaPt[iCut]);
      fHistoMCDecayGammaEtapPt[iCut]       = new TH1F("MC_DecayGammaEtap_Pt", "MC_DecayGammaEtap_Pt", nBinsPt, arrPtBinning);
      fMCList[iCut]->Add(fHistoMCDecayGammaEtapPt[iCut]);
      fHistoMCDecayGammaPhiPt[iCut]        = new TH1F("MC_DecayGammaPhi_Pt", "MC_DecayGammaPhi_Pt", nBinsPt, arrPtBinning);
      fMCList[iCut]->Add(fHistoMCDecayGammaPhiPt[iCut]);
      fHistoMCDecayGammaSigmaPt[iCut]      = new TH1F("MC_DecayGammaSigma_Pt", "MC_DecayGammaSigma_Pt", nBinsPt, arrPtBinning);
      fMCList[iCut]->Add(fHistoMCDecayGammaSigmaPt[iCut]);
      fHistoMCConvGammaPt[iCut]            = new TH1F("MC_ConvGamma_Pt", "MC_ConvGamma_Pt", nBinsPt, arrPtBinning);
      fMCList[iCut]->Add(fHistoMCConvGammaPt[iCut]);
      fHistoMCSecondaryConvGammaPt[iCut]  = new TH2F("MC_SecondaryConvGamma_Pt", "MC_SecondaryConvGamma_Pt", nBinsPt, arrPtBinning, 4, -0.5, 3.5);
      fHistoMCSecondaryConvGammaPt[iCut]->GetYaxis()->SetBinLabel(1,"K0s");
      fHistoMCSecondaryConvGammaPt[iCut]->GetYaxis()->SetBinLabel(2,"K0l");
      fHistoMCSecondaryConvGammaPt[iCut]->GetYaxis()->SetBinLabel(3,"Lambda");
      fHistoMCSecondaryConvGammaPt[iCut]->GetYaxis()->SetBinLabel(4,"rest");
      fMCList[iCut]->Add(fHistoMCSecondaryConvGammaPt[iCut]);

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
        fHistoMCConvGammaPt[iCut]->Sumw2();
        fHistoMCSecondaryConvGammaPt[iCut]->Sumw2();
      }

      if (fDoPhotonQA > 0 && fIsMC < 2){
        fHistoMCConvGammaR[iCut]           = new TH1F("MC_ConvGamma_R", "MC_ConvGamma_R", 800, 0, 200);
        fMCList[iCut]->Add(fHistoMCConvGammaR[iCut]);
        fHistoMCConvGammaEta[iCut]         = new TH1F("MC_ConvGamma_Eta", "MC_ConvGamma_Eta", 1000, -2, 2);
        fMCList[iCut]->Add(fHistoMCConvGammaEta[iCut]);
      }

      if(fDoMesonAnalysis){
        fHistoMCPi0Pt[iCut]                = new TH1F("MC_Pi0_Pt", "MC_Pi0_Pt", nBinsPt, arrPtBinning);
        fHistoMCPi0Pt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0Pt[iCut]);
        fHistoMCPi0WOWeightPt[iCut]        = new TH1F("MC_Pi0_WOWeights_Pt", "MC_Pi0_WOWeights_Pt", nBinsPt, arrPtBinning);
        fHistoMCPi0WOWeightPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0WOWeightPt[iCut]);

        fHistoMCEtaPt[iCut]                = new TH1F("MC_Eta_Pt", "MC_Eta_Pt", nBinsPt, arrPtBinning);
        fHistoMCEtaPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCEtaPt[iCut]);
        fHistoMCEtaWOWeightPt[iCut]        = new TH1F("MC_Eta_WOWeights_Pt", "MC_Eta_WOWeights_Pt", nBinsPt, arrPtBinning);
        fHistoMCEtaWOWeightPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCEtaWOWeightPt[iCut]);

        fHistoMCPi0WOWeightInAccPt[iCut]   = new TH1F("MC_Pi0WOWeightInAcc_Pt", "MC_Pi0WOWeightInAcc_Pt", nBinsPt, arrPtBinning);
        fHistoMCPi0WOWeightInAccPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0WOWeightInAccPt[iCut]);
        fHistoMCEtaWOWeightInAccPt[iCut]   = new TH1F("MC_EtaWOWeightInAcc_Pt", "MC_EtaWOWeightInAcc_Pt", nBinsPt, arrPtBinning);
        fHistoMCEtaWOWeightInAccPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCEtaWOWeightInAccPt[iCut]);
        fHistoMCPi0InAccPt[iCut]           = new TH1F("MC_Pi0InAcc_Pt", "MC_Pi0InAcc_Pt", nBinsPt, arrPtBinning);
        fHistoMCPi0InAccPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0InAccPt[iCut]);
        fHistoMCEtaInAccPt[iCut]           = new TH1F("MC_EtaInAcc_Pt", "MC_EtaInAcc_Pt", nBinsPt, arrPtBinning);
        fHistoMCEtaInAccPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCEtaInAccPt[iCut]);

        if (fIsMC > 1){
          fHistoMCPi0WOEvtWeightPt[iCut]   = new TH1F("MC_Pi0_WOEventWeights_Pt", "MC_Pi0_WOEventWeights_Pt", nBinsPt, arrPtBinning);
          fMCList[iCut]->Add(fHistoMCPi0WOEvtWeightPt[iCut]);
          fHistoMCEtaWOEvtWeightPt[iCut]   = new TH1F("MC_Eta_WOEventWeights_Pt", "MC_Eta_WOEventWeights_Pt", nBinsPt, arrPtBinning);
          fMCList[iCut]->Add(fHistoMCEtaWOEvtWeightPt[iCut]);
          fHistoMCPi0WOEvtWeightInAccPt[iCut]   = new TH1F("MC_Pi0WOEvtWeightInAcc_Pt", "MC_Pi0WOEvtWeightInAcc_Pt", nBinsPt, arrPtBinning);
          fHistoMCPi0WOEvtWeightInAccPt[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCPi0WOEvtWeightInAccPt[iCut]);
          fHistoMCEtaWOEvtWeightInAccPt[iCut]   = new TH1F("MC_EtaWOEvtWeightInAcc_Pt", "MC_EtaWOEvtWeightInAcc_Pt", nBinsPt, arrPtBinning);
          fHistoMCEtaWOEvtWeightInAccPt[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCEtaWOEvtWeightInAccPt[iCut]);

          if (fDoMesonQA > 0 && fIsMC == 2){
            fHistoMCPi0PtJetPt[iCut]       = new TH2F("MC_Pi0_Pt_JetPt", "MC_Pi0_Pt_JetPt", nBinsPt, arrPtBinning, 200, 0, 200);
            fHistoMCPi0PtJetPt[iCut]->Sumw2();
            fMCList[iCut]->Add(fHistoMCPi0PtJetPt[iCut]);
            fHistoMCEtaPtJetPt[iCut]       = new TH2F("MC_Eta_Pt_JetPt", "MC_Eta_Pt_JetPt", nBinsPt, arrPtBinning, 200, 0, 200);
            fHistoMCEtaPtJetPt[iCut]->Sumw2();
            fMCList[iCut]->Add(fHistoMCEtaPtJetPt[iCut]);
          }
        }

        fHistoMCPrimaryPtvsSource[iCut]  = new TH2F("MC_Primary_Pt_Source", "MC_Primary_Pt_Source", nBinsPt, arrPtBinning, 10, -0.5, 9.5);
        fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(1,"Pi+");
        fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(2,"Pi-");
        fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(3,"K+");
        fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(4,"K-");
        fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(5,"K0s");
        fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(6,"K0l");
        fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(7,"Lambda");
        fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(8,"Omega");
        fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(9,"Phi");
        fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(10,"Rho0");
        fMCList[iCut]->Add(fHistoMCPrimaryPtvsSource[iCut]);
        fHistoMCSecPi0Source[iCut]    = new TH1F("MC_SecPi0_Source", "MC_SecPi0_Source", 5000, 0., 5000);
        fMCList[iCut]->Add(fHistoMCSecPi0Source[iCut]);
        fHistoMCSecEtaSource[iCut]    = new TH1F("MC_SecEta_Source", "MC_SecEta_Source", 5000, 0, 5000);
        fMCList[iCut]->Add(fHistoMCSecEtaSource[iCut]);
        fHistoMCSecPi0PtvsSource[iCut]   = new TH2F("MC_SecPi0_Pt_Source", "MC_SecPi0_Pt_Source", nBinsPt, arrPtBinning, 16, -0.5, 15.5);
        fMCList[iCut]->Add(fHistoMCSecPi0PtvsSource[iCut]);
        fHistoMCSecPi0InAccPtvsSource[iCut] = new TH2F("MC_SecPi0InAcc_Pt_Source", "MC_SecPi0InAcc_Pt_Source", nBinsPt, arrPtBinning, 16, -0.5, 15.5);
        fMCList[iCut]->Add(fHistoMCSecPi0InAccPtvsSource[iCut]);
        fHistoMCSecEtaPt[iCut]           = new TH1F("MC_SecEta_Pt", "MC_SecEta_Pt", nBinsPt, arrPtBinning);
        fMCList[iCut]->Add(fHistoMCSecEtaPt[iCut]);

        if (fIsMC > 1){
          fHistoMCPrimaryPtvsSource[iCut]->Sumw2();
          fHistoMCSecPi0PtvsSource[iCut]->Sumw2();
          fHistoMCSecPi0InAccPtvsSource[iCut]->Sumw2();
          fHistoMCSecEtaPt[iCut]->Sumw2();
        }

        if (fDoMesonQA > 0){
          fHistoMCPi0PtAlpha[iCut]         = new TH2F("MC_Pi0_Pt_Alpha", "MC_Pi0_Pt_Alpha", nBinsQAPt, arrQAPtBinning, 100, 0, 1);
          fMCList[iCut]->Add(fHistoMCPi0PtAlpha[iCut]);
          fHistoMCEtaPtAlpha[iCut]         = new TH2F("MC_Eta_Pt_Alpha", "MC_Eta_Pt_Alpha", nBinsQAPt, arrQAPtBinning, 100, 0, 1);
          fMCList[iCut]->Add(fHistoMCEtaPtAlpha[iCut]);

          if (fIsMC < 2){
            fHistoMCPi0PtY[iCut]           = new TH2F("MC_Pi0_Pt_Y", "MC_Pi0_Pt_Y", nBinsQAPt, arrQAPtBinning, 150, -1.5, 1.5);
            fHistoMCPi0PtY[iCut]->Sumw2();
            fMCList[iCut]->Add(fHistoMCPi0PtY[iCut]);
            fHistoMCEtaPtY[iCut]           = new TH2F("MC_Eta_Pt_Y", "MC_Eta_Pt_Y", nBinsQAPt, arrQAPtBinning, 150, -1.5, 1.5);
            fHistoMCEtaPtY[iCut]->Sumw2();
            fMCList[iCut]->Add(fHistoMCEtaPtY[iCut]);
            fHistoMCSecPi0RvsSource[iCut]  = new TH2F("MC_SecPi0_R3D_Source", "MC_SecPi0_R3D_Source", 500, 0.0, 20., 16, -0.5, 15.5);
            fHistoMCSecPi0RvsSource[iCut]->Sumw2();
            fMCList[iCut]->Add(fHistoMCSecPi0RvsSource[iCut]);
          }
        }

      }

      fTrueList[iCut]                 = new TList();
      fTrueList[iCut]->SetName(Form("%s_%s_%s True histograms",cutstringEvent.Data() ,cutstringPhoton.Data(),cutstringMeson.Data()));
      fTrueList[iCut]->SetOwner(kTRUE);
      fCutFolder[iCut]->Add(fTrueList[iCut]);

      fHistoTrueConvGammaPt[iCut]       = new TH1F("ESD_TrueConvGamma_Pt", "ESD_TrueConvGamma_Pt", nBinsPt, arrPtBinning);
      fTrueList[iCut]->Add(fHistoTrueConvGammaPt[iCut]);

      fHistoDoubleCountTrueConvGammaRPt[iCut]  = new TH2F("ESD_TrueDoubleCountConvGamma_R_Pt", "ESD_TrueDoubleCountConvGamma_R_Pt", 800, 0, 200, nBinsQAPt, arrQAPtBinning);
      fTrueList[iCut]->Add(fHistoDoubleCountTrueConvGammaRPt[iCut]);
      fHistoMultipleCountTrueConvGamma[iCut]   = new TH1F("ESD_TrueMultipleCountConvGamma", "ESD_TrueMultipleCountConvGamma", 10, 1, 11);
      fTrueList[iCut]->Add(fHistoMultipleCountTrueConvGamma[iCut]);

      if(fEnableBDT){
        fHistoBDToutputMCTrue[iCut] = new TH1F("BDT_output_MCTrue", "BDT_output_MCTrue", 200, -1, 1);
        fTrueList[iCut]->Add(fHistoBDToutputMCTrue[iCut]);
      }

      fHistoCombinatorialPt[iCut]           = new TH2F("ESD_TrueCombinatorial_Pt", "ESD_TrueCombinatorial_Pt", nBinsPt, arrPtBinning, 16, -0.5, 15.5);
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

      fHistoTruePrimaryConvGammaPt[iCut]    = new TH1F("ESD_TruePrimaryConvGamma_Pt", "ESD_TruePrimaryConvGamma_Pt", nBinsPt, arrPtBinning);
      fTrueList[iCut]->Add(fHistoTruePrimaryConvGammaPt[iCut]);
      fHistoTrueSecondaryConvGammaPt[iCut]  = new TH2F("ESD_TrueSecondaryConvGamma_Pt", "ESD_TrueSecondaryConvGamma_Pt", nBinsPt, arrPtBinning, 4, -0.5, 3.5);
      fHistoTrueSecondaryConvGammaPt[iCut]->GetYaxis()->SetBinLabel(1,"K0s");
      fHistoTrueSecondaryConvGammaPt[iCut]->GetYaxis()->SetBinLabel(2,"K0l");
      fHistoTrueSecondaryConvGammaPt[iCut]->GetYaxis()->SetBinLabel(3,"Lambda");
      fHistoTrueSecondaryConvGammaPt[iCut]->GetYaxis()->SetBinLabel(4,"rest");
      fTrueList[iCut]->Add(fHistoTrueSecondaryConvGammaPt[iCut]);
      fHistoTrueSecondaryConvGammaMCPt[iCut]  = new TH2F("ESD_TrueSecondaryConvGamma_MCPt", "ESD_TrueSecondaryConvGamma_MCPt", nBinsPt, arrPtBinning, 4, -0.5, 3.5);
      fHistoTrueSecondaryConvGammaMCPt[iCut]->GetYaxis()->SetBinLabel(1,"K0s");
      fHistoTrueSecondaryConvGammaMCPt[iCut]->GetYaxis()->SetBinLabel(2,"K0l");
      fHistoTrueSecondaryConvGammaMCPt[iCut]->GetYaxis()->SetBinLabel(3,"Lambda");
      fHistoTrueSecondaryConvGammaMCPt[iCut]->GetYaxis()->SetBinLabel(4,"rest");
      fTrueList[iCut]->Add(fHistoTrueSecondaryConvGammaMCPt[iCut]);

      if(fDoPhotonQA == 3 && fIsMC < 2){
        if(fIsHeavyIon == 1){
          fHistoCombinatorialMothersPt[iCut]           = new TH3F("ESD_TrueCombinatorialMothers_Pt", "ESD_TrueCombinatorialMothers_Pt", 6, 0, 6, 13, 0, 13, 250, 0., 25);
          fHistoCombinatorialMothersPt[iCut]->GetXaxis()->SetBinLabel( 1,"Elec");
          fHistoCombinatorialMothersPt[iCut]->GetXaxis()->SetBinLabel( 2,"Pion");
          fHistoCombinatorialMothersPt[iCut]->GetXaxis()->SetBinLabel( 3,"Kaon");
          fHistoCombinatorialMothersPt[iCut]->GetXaxis()->SetBinLabel( 4,"Proton");
          fHistoCombinatorialMothersPt[iCut]->GetXaxis()->SetBinLabel( 5,"Rest");
          fHistoCombinatorialMothersPt[iCut]->GetXaxis()->SetBinLabel( 6,"Ancestor NP");
          fHistoCombinatorialMothersPt[iCut]->GetYaxis()->SetBinLabel( 1,"Elec");
          fHistoCombinatorialMothersPt[iCut]->GetYaxis()->SetBinLabel( 2,"Pion");
          fHistoCombinatorialMothersPt[iCut]->GetYaxis()->SetBinLabel( 3,"Kaon");
          fHistoCombinatorialMothersPt[iCut]->GetYaxis()->SetBinLabel( 4,"Proton");
          fHistoCombinatorialMothersPt[iCut]->GetYaxis()->SetBinLabel( 5,"Pi0");
          fHistoCombinatorialMothersPt[iCut]->GetYaxis()->SetBinLabel( 6,"#eta");
          fHistoCombinatorialMothersPt[iCut]->GetYaxis()->SetBinLabel( 7,"#omega");
          fHistoCombinatorialMothersPt[iCut]->GetYaxis()->SetBinLabel( 8,"#phi");
          fHistoCombinatorialMothersPt[iCut]->GetYaxis()->SetBinLabel( 9,"#eta'");
          fHistoCombinatorialMothersPt[iCut]->GetYaxis()->SetBinLabel(10,"K0s");
          fHistoCombinatorialMothersPt[iCut]->GetYaxis()->SetBinLabel(11,"#Lambda");
          fHistoCombinatorialMothersPt[iCut]->GetYaxis()->SetBinLabel(12,"#rho^{0}, #rho^{#pm}");
          fHistoCombinatorialMothersPt[iCut]->GetYaxis()->SetBinLabel(13,"Other");
          fTrueList[iCut]->Add(fHistoCombinatorialMothersPt[iCut]);
        }

        fHistoCombinatorialPtDeltaPhi_ek[iCut]  = new TH2F("ESD_TrueCombinatorial_Pt_DeltaPhi_ek", "ESD_TrueCombinatorial_Pt_DeltaPhi_ek", nBinsQAPt, arrQAPtBinning, 90, -0.5*TMath::Pi(), 0.5*TMath::Pi());
        fTrueList[iCut]->Add(fHistoCombinatorialPtDeltaPhi_ek[iCut]);
        fHistoCombinatorialPtDeltaPhi_ep[iCut]  = new TH2F("ESD_TrueCombinatorial_Pt_DeltaPhi_ep", "ESD_TrueCombinatorial_Pt_DeltaPhi_ep", nBinsQAPt, arrQAPtBinning, 90, -0.5*TMath::Pi(), 0.5*TMath::Pi());
        fTrueList[iCut]->Add(fHistoCombinatorialPtDeltaPhi_ep[iCut]);
        fHistoCombinatorialPtDeltaPhi_epi[iCut] = new TH2F("ESD_TrueCombinatorial_Pt_DeltaPhi_epi", "ESD_TrueCombinatorial_Pt_DeltaPhi_epi", nBinsQAPt, arrQAPtBinning, 90, -0.5*TMath::Pi(), 0.5*TMath::Pi());
        fTrueList[iCut]->Add(fHistoCombinatorialPtDeltaPhi_epi[iCut]);
        fHistoCombinatorialPtDeltaPhi_pik[iCut] = new TH2F("ESD_TrueCombinatorial_Pt_DeltaPhi_pik", "ESD_TrueCombinatorial_Pt_DeltaPhi_pik", nBinsQAPt, arrQAPtBinning, 90, -0.5*TMath::Pi(), 0.5*TMath::Pi());
        fTrueList[iCut]->Add(fHistoCombinatorialPtDeltaPhi_pik[iCut]);
        fHistoCombinatorialPtDeltaPhi_pip[iCut] = new TH2F("ESD_TrueCombinatorial_Pt_DeltaPhi_pip", "ESD_TrueCombinatorial_Pt_DeltaPhi_pip", nBinsQAPt, arrQAPtBinning, 90, -0.5*TMath::Pi(), 0.5*TMath::Pi());
        fTrueList[iCut]->Add(fHistoCombinatorialPtDeltaPhi_pip[iCut]);
      }

      fHistoTrueDalitzPsiPairDeltaPhi[iCut]               = new TH2F("ESD_TrueDalitzPsiPairDeltaPhi", "ESD_TrueDalitzPsiPairDeltaPhi", 100, -0.5, 2, 100, -0.5, 0.5);
      fTrueList[iCut]->Add(fHistoTrueDalitzPsiPairDeltaPhi[iCut]);
      fHistoTrueGammaPsiPairDeltaPhi[iCut]                = new TH2F("ESD_TrueGammaPsiPairDeltaPhi", "ESD_TrueGammaPsiPairDeltaPhi", 100, -0.5, 2, 100, -0.5, 0.5);
      fTrueList[iCut]->Add(fHistoTrueGammaPsiPairDeltaPhi[iCut]);

      fHistoTruePrimaryConvGammaESDPtMCPt[iCut]           = new TH2F("ESD_TruePrimaryConvGammaESD_PtMCPt", "ESD_TruePrimaryConvGammaESD_PtMCPt",
                                                                     nBinsPt, arrPtBinning, nBinsPt, arrPtBinning);
      fTrueList[iCut]->Add(fHistoTruePrimaryConvGammaESDPtMCPt[iCut]);


      fHistoTrueSecondaryConvGammaFromXFromK0sMCPtESDPt[iCut]     = new TH2F("ESD_TrueSecondaryConvGammaFromXFromK0sESD_MCPtPt", "ESD_TrueSecondaryConvGammaFromXFromK0sESD_MCPtPt",
                                                                             nBinsPt, arrPtBinning, nBinsPt, arrPtBinning);
      fTrueList[iCut]->Add(fHistoTrueSecondaryConvGammaFromXFromK0sMCPtESDPt[iCut]);
      fHistoTrueSecondaryConvGammaFromXFromK0lMCPtESDPt[iCut]     = new TH2F("ESD_TrueSecondaryConvGammaFromXFromK0lESD_MCPtPt", "ESD_TrueSecondaryConvGammaFromXFromK0lESD_MCPtPt",
                                                                             nBinsPt, arrPtBinning, nBinsPt, arrPtBinning);
      fTrueList[iCut]->Add(fHistoTrueSecondaryConvGammaFromXFromK0lMCPtESDPt[iCut]);
      fHistoTrueSecondaryConvGammaFromXFromLambdaMCPtESDPt[iCut]  = new TH2F("ESD_TrueSecondaryConvGammaFromXFromLambdaESD_MCPtPt", "ESD_TrueSecondaryConvGammaFromXFromLambdaESD_MCPtPt",
                                                                             nBinsPt, arrPtBinning, nBinsPt, arrPtBinning);
      fTrueList[iCut]->Add(fHistoTrueSecondaryConvGammaFromXFromLambdaMCPtESDPt[iCut]);


      if ((fIsMC > 1)  || (fDoMaterialBudgetWeightingOfGammasForTrueMesons && fIsMC > 0) ) {
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
        fHistoTrueDalitzPsiPairDeltaPhi[iCut]->Sumw2();
        fHistoTrueGammaPsiPairDeltaPhi[iCut]->Sumw2();
      }

      if(fDoJetAnalysis && !fDoLightOutput){
        fTrueJetHistograms[iCut] = new TList();
        fTrueJetHistograms[iCut]->SetName(Form("%s_%s_%s True Jet histograms", cutstringEvent.Data() ,cutstringPhoton.Data(),cutstringMeson.Data()));
        fTrueJetHistograms[iCut]->SetOwner(kTRUE);
        fCutFolder[iCut]->Add(fTrueJetHistograms[iCut]);
        fHistoTruevsRecJetPt[iCut] = new TH2F("True_JetPt_vs_Rec_JetPt", "True_JetPt_vs_Rec_JetPt", 150, 0, 150, 150, 0, 150);
        fTrueJetHistograms[iCut]->Add(fHistoTruevsRecJetPt[iCut]);
        fHistoTrueJetMotherInvMassPt[iCut] = new TH2F("ESD_True_Jet_InvMass_Pt", "ESD_True_Jet_InvMass_Pt", 800, 0, 0.8, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoTrueJetMotherInvMassPt[iCut]);
        fHistoTruePrimaryJetInvMassPt[iCut] = new TH2F("ESD_TruePrimaryJet_InvMass_Pt", "ESD_TruePrimaryJet_InvMass_Pt", 800, 0, 0.8, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoTruePrimaryJetInvMassPt[iCut]);
        fHistoTrueInJetMotherInvMassPt[iCut] = new TH2F("ESD_True_inJet_InvMass_Pt", "ESD_True_inJet_InvMass_Pt", 800, 0, 0.8, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoTrueInJetMotherInvMassPt[iCut]);
        fHistoTruePrimaryinJetInvMassPt[iCut] = new TH2F("ESD_TruePrimaryinJet_InvMass_Pt", "ESD_TruePrimaryinJet_InvMass_Pt", 800, 0, 0.8, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoTruePrimaryinJetInvMassPt[iCut]);
        fHistoTrueDoubleCountingJet[iCut] = new TH1F("Double_Counting_True_inJet", "Double_Counting_True_inJet", 6, 0, 6);
        fTrueJetHistograms[iCut]->Add(fHistoTrueDoubleCountingJet[iCut]);
        fHistoTrueJetFragmFunc[iCut] = new TH2F("ESD_TrueinJetPt_FragmentationFunc", "ESD_TrueinJetPt_FragmentationFunc", 50, arrLogBinning, 150, 0., 150.);
        fTrueJetHistograms[iCut]->Add(fHistoTrueJetFragmFunc[iCut]);
        fHistoTrueJetFragmFuncZInvMass[iCut] = new TH2F("ESD_TrueinJetPt_Fragm_Z_InvMass", "ESD_TrueinJetPt_Fragm_Z_InvMass", 800, 0, 0.8, 50, arrLogBinning);
        fTrueJetHistograms[iCut]->Add(fHistoTrueJetFragmFuncZInvMass[iCut]);
        fHistoMCPi0JetInAccPt[iCut]      = new TH1F("MC_Pi0JetInAcc_Pt", "MC_Pi0JetInAcc_Pt", nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoMCPi0JetInAccPt[iCut]);
        fHistoMCPi0inJetInAccPt[iCut]      = new TH1F("MC_Pi0inJetInAcc_Pt", "MC_Pi0inJetInAcc_Pt", nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoMCPi0inJetInAccPt[iCut]);
        fHistoMCPi0JetEventGenerated[iCut]    = new TH1F("MC_Pi0_JetEvent_Generated", "MC_Pi0_JetEvent_Generated", nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoMCPi0JetEventGenerated[iCut]);
        fHistoMCPi0inJetGenerated[iCut]    = new TH1F("MC_Pi0_inJet_Generated", "MC_Pi0_inJet_Generated", nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoMCPi0inJetGenerated[iCut]);
        fHistoTrueSecondaryFromK0sJetInvMassPt[iCut] = new TH2F("ESD_TrueSecondaryFromK0sJet_InvMass_Pt", "ESD_TrueSecondaryFromK0sJet_InvMass_Pt", 800, 0, 0.8, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoTrueSecondaryFromK0sJetInvMassPt[iCut]);
        fHistoTrueSecondaryFromK0sinJetInvMassPt[iCut] = new TH2F("ESD_TrueSecondaryFromK0s_inJet_InvMass_Pt", "ESD_TrueSecondaryFromK0s_inJet_InvMass_Pt", 800, 0, 0.8, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoTrueSecondaryFromK0sinJetInvMassPt[iCut]);
        fHistoTrueSecondaryFromLambdaJetInvMassPt[iCut] = new TH2F("ESD_TrueSecondaryFromLambdaJet_InvMass_Pt", "ESD_TrueSecondaryFromLambdaJet_InvMass_Pt", 800, 0, 0.8, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoTrueSecondaryFromLambdaJetInvMassPt[iCut]);
        fHistoTrueSecondaryFromLambdainJetInvMassPt[iCut] = new TH2F("ESD_TrueSecondaryFromLambda_inJet_InvMass_Pt", "ESD_TrueSecondaryFromLambda_inJet_InvMass_Pt", 800, 0, 0.8, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoTrueSecondaryFromLambdainJetInvMassPt[iCut]);
        fHistoTrueSecondaryFromK0lJetInvMassPt[iCut] = new TH2F("ESD_TrueSecondaryFromK0lJet_InvMass_Pt", "ESD_TrueSecondaryFromK0lJet_InvMass_Pt", 800, 0, 0.8, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoTrueSecondaryFromK0lJetInvMassPt[iCut]);
        fHistoTrueSecondaryFromK0linJetInvMassPt[iCut] = new TH2F("ESD_TrueSecondaryFromK0l_inJet_InvMass_Pt", "ESD_TrueSecondaryFromK0l_inJet_InvMass_Pt", 800, 0, 0.8, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoTrueSecondaryFromK0linJetInvMassPt[iCut]);
        fHistoTrueSecondaryInvJetMassPt[iCut] = new TH2F("ESD_TrueSecondaryJet_InvMass_Pt", "ESD_TrueSecondaryJet_InvMass_Pt", 800, 0, 0.8, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoTrueSecondaryInvJetMassPt[iCut]);
        fHistoTrueSecondaryInvinJetMassPt[iCut] = new TH2F("ESD_TrueSecondary_inJet_InvMass_Pt", "ESD_TrueSecondary_inJet_InvMass_Pt", 800, 0, 0.8, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoTrueSecondaryInvinJetMassPt[iCut]);
        fHistoMCEtaJetInAccPt[iCut]      = new TH1F("MC_EtaJetInAcc_Pt", "MC_EtaJetInAcc_Pt", nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoMCEtaJetInAccPt[iCut]);
        fHistoMCEtainJetInAccPt[iCut]      = new TH1F("MC_EtainJetInAcc_Pt", "MC_EtainJetInAcc_Pt", nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoMCEtainJetInAccPt[iCut]);
        fHistoMCEtaJetEventGenerated[iCut]    = new TH1F("MC_Eta_JetEvent_Generated", "MC_Eta_JetEvent_Generated", nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoMCEtaJetEventGenerated[iCut]);
        fHistoMCEtainJetGenerated[iCut]    = new TH1F("MC_Eta_inJet_Generated", "MC_Eta_inJet_Generated", nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoMCEtainJetGenerated[iCut]);
      }
      if(fDoJetQA){
        if(fDoLightOutput){
          fTrueJetHistograms[iCut] = new TList();
          fTrueJetHistograms[iCut]->SetName(Form("%s_%s_%s True Jet histograms", cutstringEvent.Data() ,cutstringPhoton.Data(),cutstringMeson.Data()));
          fTrueJetHistograms[iCut]->SetOwner(kTRUE);
          fCutFolder[iCut]->Add(fTrueJetHistograms[iCut]);

          fHistoTruevsRecJetPt[iCut] = new TH2F("True_JetPt_vs_Rec_JetPt", "True_JetPt_vs_Rec_JetPt", 150, 0, 150, 150, 0, 150);
          fTrueJetHistograms[iCut]->Add(fHistoTruevsRecJetPt[iCut]);
        }
        fHistoUnfoldingAsData[iCut]      = new TH2F("Unfolding_AsData", "Unfolding_AsData", 800, 0, 0.8, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoUnfoldingAsData[iCut]);
        fHistoUnfoldingMissed[iCut]      = new TH2F("Unfolding_Missed", "Unfolding_Missed", 800, 0, 0.8, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoUnfoldingMissed[iCut]);
        fHistoUnfoldingReject[iCut]      = new TH2F("Unfolding_Reject", "Unfolding_Reject", 800, 0, 0.8, nBinsPt, arrPtBinning);
        fTrueJetHistograms[iCut]->Add(fHistoUnfoldingReject[iCut]);
        fHistoUnfoldingAsDataInvMassZ[iCut]      = new TH2F("Unfolding_AsData_InvMass_Z", "Unfolding_AsData_InvMass_Z", 800, 0, 0.8, 50, arrLogBinning);
        fTrueJetHistograms[iCut]->Add(fHistoUnfoldingAsDataInvMassZ[iCut]);
        fHistoUnfoldingMissedInvMassZ[iCut]      = new TH2F("Unfolding_Missed_InvMass_Z", "Unfolding_Missed_InvMass_Z", 800, 0, 0.8, 50, arrLogBinning);
        fTrueJetHistograms[iCut]->Add(fHistoUnfoldingMissedInvMassZ[iCut]);
        fHistoUnfoldingRejectInvMassZ[iCut]      = new TH2F("Unfolding_Reject_InvMass_Z", "Unfolding_Reject_InvMass_Z", 800, 0, 0.8, 50, arrLogBinning);
        fTrueJetHistograms[iCut]->Add(fHistoUnfoldingRejectInvMassZ[iCut]);
      }

      if (fDoPhotonQA > 0 ){
        if (fIsMC < 2){
          fHistoTrueConvGammaPsiPairPt[iCut]      = new TH2F("ESD_TrueConvGamma_PsiPair_Pt", "ESD_TrueConvGamma_PsiPair_Pt", 500, 0, 5, nBinsPt, arrPtBinning);
          fTrueList[iCut]->Add(fHistoTrueConvGammaPsiPairPt[iCut]);
          fHistoTrueConvGammaEta[iCut]            = new TH1F("ESD_TrueConvGamma_Eta", "ESD_TrueConvGamma_Eta", 1000, -2, 2);
          fTrueList[iCut]->Add(fHistoTrueConvGammaEta[iCut]);
          fHistoTrueConvGammaR[iCut]              = new TH1F("ESD_TrueConvGamma_R", "ESD_TrueConvGamma_R", 800, 0, 200);
          fTrueList[iCut]->Add(fHistoTrueConvGammaR[iCut]);
          fHistoTrueConvGammaRMC[iCut]            = new TH1F("ESD_TrueConvGamma_RMC", "ESD_TrueConvGamma_RMC", 800, 0, 200);
          fTrueList[iCut]->Add(fHistoTrueConvGammaRMC[iCut]);
          if ((fDoPhotonQA == 4)||(fDoPhotonQA == 5)){
            fHistoTrueConvGammaInvMass[iCut]        = new TH1F("ESD_TrueConvGamma_InvMass", "ESD_TrueConvGamma_InvMass", 500, 0., 5.);
            fTrueList[iCut]->Add(fHistoTrueConvGammaInvMass[iCut]);
            fHistoTrueConvGammaInvMassReco[iCut]    = new TH1F("ESD_TrueConvGamma_InvMassReco", "ESD_TrueConvGamma_InvMassReco", 500, 0., 5.);
            fTrueList[iCut]->Add(fHistoTrueConvGammaInvMassReco[iCut]);
          }
        }
        fHistoTrueConvGammaPtMC[iCut]           = new TH1F("ESD_TrueConvGamma_PtMC", "ESD_TrueConvGamma_PtMC", nBinsPt, arrPtBinning);
        fTrueList[iCut]->Add(fHistoTrueConvGammaPtMC[iCut]);
        if ((fIsMC > 1) || (fDoMaterialBudgetWeightingOfGammasForTrueMesons && fIsMC > 0) )
          fHistoTrueConvGammaPtMC[iCut]->Sumw2();

      }

      if(fDoMesonAnalysis){
        fHistoTrueMotherInvMassPt[iCut]         = new TH2F("ESD_TrueMother_InvMass_Pt", "ESD_TrueMother_InvMass_Pt", 800, 0, 0.8, nBinsPt, arrPtBinning);
        fTrueList[iCut]->Add(fHistoTrueMotherInvMassPt[iCut]);
        fHistoDoubleCountTruePi0InvMassPt[iCut]    = new TH2F("ESD_TrueDoubleCountPi0_InvMass_Pt", "ESD_TrueDoubleCountPi0_InvMass_Pt", 800, 0, 0.8, nBinsPt, arrPtBinning);
        fTrueList[iCut]->Add(fHistoDoubleCountTruePi0InvMassPt[iCut]);
        fHistoMultipleCountTruePi0[iCut]           = new TH1F("ESD_TrueMultipleCountPi0", "ESD_TrueMultipleCountPi0", 10, 1, 11);
        fTrueList[iCut]->Add(fHistoMultipleCountTruePi0[iCut]);
        fHistoDoubleCountTrueEtaInvMassPt[iCut]    = new TH2F("ESD_TrueDoubleCountEta_InvMass_Pt", "ESD_TrueDoubleCountEta_InvMass_Pt", 800, 0, 0.8, nBinsPt, arrPtBinning);
        fTrueList[iCut]->Add(fHistoDoubleCountTrueEtaInvMassPt[iCut]);
        fHistoMultipleCountTrueEta[iCut]           = new TH1F("ESD_TrueMultipleCountEta", "ESD_TrueMultipleCountEta", 10, 1, 11);
        fTrueList[iCut]->Add(fHistoMultipleCountTrueEta[iCut]);
        fHistoTruePrimaryMotherInvMassPt[iCut]  = new TH2F("ESD_TruePrimaryMother_InvMass_Pt", "ESD_TruePrimaryMother_InvMass_Pt", 800, 0, 0.8, nBinsPt, arrPtBinning);
        fHistoTruePrimaryMotherInvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTruePrimaryMotherInvMassPt[iCut]);
        fHistoTruePrimaryMotherW0WeightingInvMassPt[iCut]   = new TH2F("ESD_TruePrimaryMotherW0Weights_InvMass_Pt", "ESD_TruePrimaryMotherW0Weights_InvMass_Pt", 800, 0, 0.8, nBinsPt, arrPtBinning);
        fHistoTruePrimaryMotherW0WeightingInvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTruePrimaryMotherW0WeightingInvMassPt[iCut]);
        pESDTruePrimaryMotherWeightsInvMassPt[iCut]       = new TProfile2D("ESD_TruePrimaryMotherWeights_InvMass_Pt", "ESD_TruePrimaryMotherWeights_InvMass_Pt", 800, 0, 0.8, nBinsPt, arrPtBinning);
        pESDTruePrimaryMotherWeightsInvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(pESDTruePrimaryMotherWeightsInvMassPt[iCut]);
        fHistoTrueSecondaryMotherInvMassPt[iCut]            = new TH2F("ESD_TrueSecondaryMother_InvMass_Pt", "ESD_TrueSecondaryMother_InvMass_Pt", 800, 0, 0.8, nBinsPt, arrPtBinning);
        fHistoTrueSecondaryMotherInvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTrueSecondaryMotherInvMassPt[iCut]);
        fHistoTrueSecondaryMotherFromK0sInvMassPt[iCut]     = new TH2F("ESD_TrueSecondaryMotherFromK0s_InvMass_Pt", "ESD_TrueSecondaryMotherFromK0s_InvMass_Pt", 800, 0, 0.8, nBinsPt, arrPtBinning);
        fHistoTrueSecondaryMotherFromK0sInvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTrueSecondaryMotherFromK0sInvMassPt[iCut]);
        fHistoTrueSecondaryMotherFromK0lInvMassPt[iCut]     = new TH2F("ESD_TrueSecondaryMotherFromK0l_InvMass_Pt", "ESD_TrueSecondaryMotherFromK0l_InvMass_Pt", 800, 0, 0.8, nBinsPt, arrPtBinning);
        fHistoTrueSecondaryMotherFromK0lInvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTrueSecondaryMotherFromK0lInvMassPt[iCut]);
        fHistoTrueSecondaryMotherFromEtaInvMassPt[iCut]     = new TH2F("ESD_TrueSecondaryMotherFromEta_InvMass_Pt", "ESD_TrueSecondaryMotherFromEta_InvMass_Pt", 800, 0, 0.8, nBinsPt, arrPtBinning);
        fHistoTrueSecondaryMotherFromEtaInvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTrueSecondaryMotherFromEtaInvMassPt[iCut]);
        fHistoTrueSecondaryMotherFromLambdaInvMassPt[iCut]  = new TH2F("ESD_TrueSecondaryMotherFromLambda_InvMass_Pt", "ESD_TrueSecondaryMotherFromLambda_InvMass_Pt", 800, 0, 0.8, nBinsPt, arrPtBinning);
        fHistoTrueSecondaryMotherFromLambdaInvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTrueSecondaryMotherFromLambdaInvMassPt[iCut]);

        if(fIsMC < 2){
          fHistoTrueMotherDalitzInvMassPt[iCut] = new TH2F("ESD_TrueDalitz_InvMass_Pt", "ESD_TrueDalitz_InvMass_Pt", 800, 0, 0.8, nBinsPt, arrPtBinning);
          fTrueList[iCut]->Add(fHistoTrueMotherDalitzInvMassPt[iCut]);
        }

        if (fIsMC > 1){
          fHistoTrueMotherInvMassPt[iCut]->Sumw2();
          fHistoDoubleCountTruePi0InvMassPt[iCut]->Sumw2();
          fHistoMultipleCountTruePi0[iCut]->Sumw2();

          fHistoDoubleCountTrueEtaInvMassPt[iCut]->Sumw2();
          fHistoMultipleCountTrueEta[iCut]->Sumw2();
        }

        if (fDoMesonQA > 0){
          if (fIsMC < 2){
            fHistoTruePrimaryPi0MCPtResolPt[iCut]   = new TH2F("ESD_TruePrimaryPi0_MCPt_ResolPt", "ESD_TruePrimaryPi0_ResolPt_MCPt", nBinsPt, arrPtBinning, 1000, -1., 1.);
            fHistoTruePrimaryPi0MCPtResolPt[iCut]->Sumw2();
            fTrueList[iCut]->Add(fHistoTruePrimaryPi0MCPtResolPt[iCut]);
            fHistoTruePrimaryEtaMCPtResolPt[iCut]   = new TH2F("ESD_TruePrimaryEta_MCPt_ResolPt", "ESD_TruePrimaryEta_ResolPt_MCPt", nBinsPt, arrPtBinning, 1000, -1., 1.);
            fHistoTruePrimaryEtaMCPtResolPt[iCut]->Sumw2();
            fTrueList[iCut]->Add(fHistoTruePrimaryEtaMCPtResolPt[iCut]);
            fHistoTrueBckGGInvMassPt[iCut]          = new TH2F("ESD_TrueBckGG_InvMass_Pt", "ESD_TrueBckGG_InvMass_Pt", 800, 0, 0.8, nBinsPt, arrPtBinning);
            fTrueList[iCut]->Add(fHistoTrueBckGGInvMassPt[iCut]);
            fHistoTrueBckContInvMassPt[iCut]        = new TH2F("ESD_TrueBckCont_InvMass_Pt", "ESD_TrueBckCont_InvMass_Pt", 800, 0, 0.8, nBinsPt, arrPtBinning);
            fTrueList[iCut]->Add(fHistoTrueBckContInvMassPt[iCut]);
            fHistoTrueK0sWithPi0DaughterMCPt[iCut]  = new TH1F("ESD_TrueK0sWithPi0Daughter_MCPt", "ESD_TrueK0sWithPi0Daughter_MCPt", nBinsPt, arrPtBinning);
            fTrueList[iCut]->Add(fHistoTrueK0sWithPi0DaughterMCPt[iCut]);
            fHistoTrueK0lWithPi0DaughterMCPt[iCut]  = new TH1F("ESD_TrueK0lWithPi0Daughter_MCPt", "ESD_TrueK0lWithPi0Daughter_MCPt", nBinsPt, arrPtBinning);
            fTrueList[iCut]->Add(fHistoTrueK0lWithPi0DaughterMCPt[iCut]);
            fHistoTrueEtaWithPi0DaughterMCPt[iCut]  = new TH1F("ESD_TrueEtaWithPi0Daughter_MCPt", "ESD_TrueEtaWithPi0Daughter_MCPt", nBinsPt, arrPtBinning);
            fTrueList[iCut]->Add(fHistoTrueEtaWithPi0DaughterMCPt[iCut]);
            fHistoTrueLambdaWithPi0DaughterMCPt[iCut] = new TH1F("ESD_TrueLambdaWithPi0Daughter_MCPt", "ESD_TrueLambdaWithPi0Daughter_MCPt", nBinsPt, arrPtBinning);
            fTrueList[iCut]->Add(fHistoTrueLambdaWithPi0DaughterMCPt[iCut]);

            fHistoTruePi0PtY[iCut]            = new TH2F("ESD_TruePi0_Pt_Y", "ESD_TruePi0_Pt_Y", nBinsQAPt, arrQAPtBinning, 150, -1.5, 1.5);
            fTrueList[iCut]->Add(fHistoTruePi0PtY[iCut]);
            fHistoTrueEtaPtY[iCut]            = new TH2F("ESD_TrueEta_Pt_Y", "ESD_TrueEta_Pt_Y", nBinsQAPt, arrQAPtBinning, 150, -1.5, 1.5);
            fTrueList[iCut]->Add(fHistoTrueEtaPtY[iCut]);
            fHistoTruePi0PtOpenAngle[iCut]    = new TH2F("ESD_TruePi0_Pt_OpenAngle", "ESD_TruePi0_Pt_OpenAngle", nBinsQAPt, arrQAPtBinning, 200, 0, 2*TMath::Pi());
            fTrueList[iCut]->Add(fHistoTruePi0PtOpenAngle[iCut]);
            fHistoTrueEtaPtOpenAngle[iCut]    = new TH2F("ESD_TrueEta_Pt_OpenAngle", "ESD_TrueEta_Pt_OpenAngle", nBinsQAPt, arrQAPtBinning, 200, 0, 2*TMath::Pi());
            fTrueList[iCut]->Add(fHistoTrueEtaPtOpenAngle[iCut]);
          }

          fHistoTruePi0PtAlpha[iCut]          = new TH2F("ESD_TruePi0_Pt_Alpha", "ESD_TruePi0_Pt_Alpha", nBinsQAPt, arrQAPtBinning, 100, 0, 1);
          fTrueList[iCut]->Add(fHistoTruePi0PtAlpha[iCut]);
          fHistoTrueEtaPtAlpha[iCut]          = new TH2F("ESD_TrueEta_Pt_Alpha", "ESD_TrueEta_Pt_Alpha", nBinsQAPt, arrQAPtBinning, 100, 0, 1);
          fTrueList[iCut]->Add(fHistoTrueEtaPtAlpha[iCut]);

        }
      }
    }
  }

  vecDoubleCountTruePi0s.clear();
  vecDoubleCountTrueEtas.clear();
  vecDoubleCountTrueConvGammas.clear();

  mapMultipleCountTruePi0s.clear();
  mapMultipleCountTrueEtas.clear();
  mapMultipleCountTrueConvGammas.clear();

  if(fV0Reader)
    if((AliConvEventCuts*)fV0Reader->GetEventCuts())
      if(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms())
        fOutputContainer->Add(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms());

  if(fV0Reader)
    if((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())
      if(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms())
        fOutputContainer->Add(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms());

  if(fV0Reader && fV0Reader->GetProduceV0FindingEfficiency())
    if (fV0Reader->GetV0FindingEfficiencyHistograms())
      fOutputContainer->Add(fV0Reader->GetV0FindingEfficiencyHistograms());

  if(fV0Reader && fV0Reader->GetProduceImpactParamHistograms())fOutputContainer->Add(fV0Reader->GetImpactParamHistograms());

  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if(!((AliConvEventCuts*)fEventCutArray->At(iCut))) continue;
    if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms()){
      fCutFolder[iCut]->Add(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms());
    }
    if(!((AliConversionPhotonCuts*)fCutArray->At(iCut))) continue;
    if(((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutHistograms()){
      fCutFolder[iCut]->Add(((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutHistograms());
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
    if (fEnableClusterCutsForTrigger){
      if(!((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))) continue;
      if(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutHistograms()){
        fCutFolder[iCut]->Add(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutHistograms());
      }
    }

  }
  if (fIsMC > 0){
    tBrokenFiles = new TTree("BrokenFiles", "BrokenFiles");
    tBrokenFiles->Branch("fileName",&fFileNameBroken);
    fOutputContainer->Add(tBrokenFiles);
  }
  OpenFile(1);
  PostData(1, fOutputContainer);
  Int_t nContainerOutput = 2;
  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if((fDoPhotonQA == 2)||(fDoPhotonQA == 5)){
      OpenFile(nContainerOutput);
      PostData(nContainerOutput, tESDConvGammaPtDcazCat[iCut]);
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
Bool_t AliAnalysisTaskGammaConvV1::Notify()
{
  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetPeriodEnum() == AliConvEventCuts::kNoPeriod && ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum() != AliConvEventCuts::kNoPeriod){
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->SetPeriodEnumExplicit(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum());
    } else if (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetPeriodEnum() == AliConvEventCuts::kNoPeriod ){
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->SetPeriodEnum(fV0Reader->GetPeriodName());
    }

    if(!((AliConvEventCuts*)fEventCutArray->At(iCut))->GetDoEtaShift()){
      fHistoEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
      continue; // No Eta Shift requested, continue
    }
    if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift() == 0.0){ // Eta Shift requested but not set, get shift automatically
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCorrectEtaShiftFromPeriod();
      fHistoEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once
      continue;
    }
    else{
      printf(" Gamma Conversion Task %s :: Eta Shift Manually Set to %f \n\n",
          (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber()).Data(),((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift());
      fHistoEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once
    }
  }
  return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskGammaConvV1::UserExec(Option_t *)
{
  //
  // Called for each event
  //
  fInputEvent = InputEvent();

  // Set MC events
  if(fIsMC>0) fMCEvent = MCEvent();

  //calculating the weight for the centrality flattening
  for(Int_t iCut = 0; iCut<fnCuts; iCut++){
    if(fDoCentralityFlat > 0){
      fWeightCentrality[iCut] = 1.;
      fWeightCentrality[iCut] = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetWeightForCentralityFlattening(fInputEvent);
    }
  }

  Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();
  if(fInputEvent->IsIncompleteDAQ()==kTRUE) eventQuality = 2;  // incomplete event
  // Event Not Accepted due to MC event missing or because it is incomplere or  wrong trigger for V0ReaderV1 => skip broken events/files
  if(eventQuality == 2 || eventQuality == 3){
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
      if( fIsMC > 1 ) fHistoNEventsWOWeight[iCut]->Fill(eventQuality);
      if( fDoCentralityFlat > 0) fHistoNEventsWeighted[iCut]->Fill(eventQuality, fWeightCentrality[iCut]);
    }
    return;
  }

  fReaderGammas = fV0Reader->GetReconstructedGammas(); // Gammas from default Cut

  // ------------------- BeginEvent ----------------------------

  AliEventplane *EventPlane = fInputEvent->GetEventplane();
  if(fIsHeavyIon ==1)fEventPlaneAngle = EventPlane->GetEventplane("V0",fInputEvent,2);
  else fEventPlaneAngle=0.0;

  if(fIsMC > 0 && fInputEvent->IsA()==AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled())){
    RelabelAODPhotonCandidates(kTRUE);    // In case of AODMC relabeling MC
    fV0Reader->RelabelAODs(kTRUE);
  }
  for(Int_t iCut = 0; iCut<fnCuts; iCut++){
    fiCut = iCut;

    Int_t eventNotAccepted = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon,kFALSE);

    if( fIsMC == 2 ){
      Float_t xsection      = -1.;
      Float_t ntrials       = -1.;
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetXSectionAndNTrials(fMCEvent,xsection,ntrials,fInputEvent);
      if((xsection==-1.) || (ntrials==-1.)) AliFatal(Form("ERROR: GetXSectionAndNTrials returned invalid xsection/ntrials, periodName from V0Reader: '%s'",fV0Reader->GetPeriodName().Data()));
      fProfileJetJetXSection[iCut]->Fill(0.,xsection);
      fhJetJetNTrials[iCut]->Fill("#sum{NTrials}",ntrials);
    }

    if( fIsMC > 0 ){
      fWeightJetJetMC       = 1;
      Bool_t isMCJet        = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsJetJetMCEventAccepted( fMCEvent, fWeightJetJetMC , fInputEvent);
      if (fIsMC == 3){
        Double_t weightMult   = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetWeightForMultiplicity(fV0Reader->GetNumberOfPrimaryTracks());
        fWeightJetJetMC       = fWeightJetJetMC*weightMult;
      }

      if( fIsMC == 1 ) fWeightJetJetMC = 1;
      if(!isMCJet){
        fHistoNEvents[iCut]->Fill(10,fWeightJetJetMC);
        if( fIsMC > 1 ) fHistoNEventsWOWeight[iCut]->Fill(10);
        continue;
      }
    }

    if(eventNotAccepted){
      // cout << "event rejected due to wrong trigger: " <<eventNotAccepted << endl;
      fHistoNEvents[iCut]->Fill(eventNotAccepted,fWeightJetJetMC); // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1
      if( fIsMC > 1 ) fHistoNEventsWOWeight[iCut]->Fill(eventNotAccepted);
      if(fDoCentralityFlat > 0) fHistoNEventsWeighted[iCut]->Fill(eventNotAccepted, fWeightCentrality[iCut]*fWeightJetJetMC);
      continue;
    }

    if(eventQuality != 0){// Event Not Accepted
      // cout << "event rejected due to: " <<eventQuality << endl;
      fHistoNEvents[iCut]->Fill(eventQuality,fWeightJetJetMC);
      if( fIsMC > 1 ) fHistoNEventsWOWeight[iCut]->Fill(eventQuality);
      if(fDoCentralityFlat > 0) fHistoNEventsWeighted[iCut]->Fill(eventQuality, fWeightCentrality[iCut]*fWeightJetJetMC);
      continue;
    }

    if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetUseSphericity()!=0){
        if(fV0Reader->GetSphericity() != -1 && fV0Reader->GetSphericity() != 0){
          fHistoEventSphericity[iCut]->Fill(fV0Reader->GetSphericity(), fWeightJetJetMC);
          fHistoNEvents[iCut]->Fill(eventQuality,fWeightJetJetMC); // Should be 0 here
          if( fIsMC > 1 ) fHistoNEventsWOWeight[iCut]->Fill(eventQuality);
          fHistoNGoodESDTracks[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks(),fWeightJetJetMC);
          fHistoVertexZ[iCut]->Fill(fInputEvent->GetPrimaryVertex()->GetZ(),fWeightJetJetMC);
        }
    } else {
      fHistoNEvents[iCut]->Fill(eventQuality,fWeightJetJetMC); // Should be 0 here
      if( fIsMC > 1 ) fHistoNEventsWOWeight[iCut]->Fill(eventQuality);
      if(fDoCentralityFlat > 0) fHistoNEventsWeighted[iCut]->Fill(eventQuality, fWeightCentrality[iCut]*fWeightJetJetMC); // Should be 0 here

      fHistoNGoodESDTracks[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks(),fWeightJetJetMC);
      if(fDoCentralityFlat > 0) fHistoNGoodESDTracksWeighted[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks(), fWeightCentrality[iCut]*fWeightJetJetMC);

      fHistoVertexZ[iCut]->Fill(fInputEvent->GetPrimaryVertex()->GetZ(),fWeightJetJetMC);
      if(fDoCentralityFlat > 0) fHistoVertexZWeighted[iCut]->Fill(fInputEvent->GetPrimaryVertex()->GetZ(), fWeightCentrality[iCut]*fWeightJetJetMC);

      if(fDoPlotVsCentrality) fHistoCentrality[iCut]->Fill(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCentrality(fInputEvent),fWeightJetJetMC);
      if(fDoCentralityFlat > 0) fHistoCentralityFlattened[iCut]->Fill(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCentrality(fInputEvent), fWeightCentrality[iCut]*fWeightJetJetMC);

      if(fDoCentralityFlat > 0) fHistoCentralityVsPrimaryTracks[iCut]->Fill(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCentrality(fInputEvent),fV0Reader->GetNumberOfPrimaryTracks(), fWeightCentrality[iCut]*fWeightJetJetMC);
      else if(fDoPlotVsCentrality) fHistoCentralityVsPrimaryTracks[iCut]->Fill(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCentrality(fInputEvent),fV0Reader->GetNumberOfPrimaryTracks(), fWeightJetJetMC);

      if(!fDoLightOutput){
        if( fIsMC < 2 ){
          if(fDoCentralityFlat > 0) fHistoSPDClusterTrackletBackground[iCut]->Fill(fInputEvent->GetMultiplicity()->GetNumberOfTracklets(),(fInputEvent->GetNumberOfITSClusters(0)+fInputEvent->GetNumberOfITSClusters(1)), fWeightCentrality[iCut]);
          else fHistoSPDClusterTrackletBackground[iCut]->Fill(fInputEvent->GetMultiplicity()->GetNumberOfTracklets(),(fInputEvent->GetNumberOfITSClusters(0)+fInputEvent->GetNumberOfITSClusters(1)));
        }

        if(((AliConvEventCuts*)fEventCutArray->At(iCut))->IsHeavyIon() == 2){
          fHistoV0MultVsNumberTPCoutTracks[iCut]->Fill(fV0Reader->GetNumberOfTPCoutTracks(), fInputEvent->GetVZEROData()->GetMTotV0A());
        } else {
          fHistoV0MultVsNumberTPCoutTracks[iCut]->Fill(fV0Reader->GetNumberOfTPCoutTracks(), fInputEvent->GetVZEROData()->GetMTotV0A()+fInputEvent->GetVZEROData()->GetMTotV0C());
        }
        if(((AliConvEventCuts*)fEventCutArray->At(iCut))->IsHeavyIon() == 2) fHistoNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A(),fWeightJetJetMC);
        else if(fDoCentralityFlat > 0){
          fHistoNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A()+fInputEvent->GetVZEROData()->GetMTotV0C(), fWeightCentrality[iCut]*fWeightJetJetMC);
        } else {
          fHistoNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A()+fInputEvent->GetVZEROData()->GetMTotV0C(),fWeightJetJetMC);
        }
      }
    }

    if (fEnableClusterCutsForTrigger ){
      ProcessClusters();
    }

    if(fIsMC > 0){
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
            if(fIsMC < 2){
              TString nameBin= fHistoMCHeaders[iCut]->GetXaxis()->GetBinLabel(i+1);
              if (nameBin.CompareTo("")== 0){
                TString nameHeader = ((TObjString*)((TList*)((AliConvEventCuts*)fEventCutArray->At(iCut))
                                ->GetAcceptedHeader())->At(i))->GetString();
  //               cout << nameHeader << endl;
                fHistoMCHeaders[iCut]->GetXaxis()->SetBinLabel(i+1,nameHeader.Data());
              }
            }
          }
        }
      }
    }

    if( fIsMC > 0 ){
      if(fInputEvent->IsA()==AliESDEvent::Class()) ProcessMCParticles();
      if(fInputEvent->IsA()==AliAODEvent::Class()) ProcessAODMCParticles();
    }

    ProcessPhotonCandidates(); // Process this cuts gammas
    if(fEnableBDT) ProcessPhotonBDT(); //
    if(fDoJetAnalysis)   ProcessJets(); //Process jets

    if(!fDoLightOutput){
      if(fDoCentralityFlat > 0){
	fHistoNGammaCandidates[iCut]->Fill(fGammaCandidates->GetEntries(), fWeightCentrality[iCut]*fWeightJetJetMC);
	if( fIsMC < 2 ) fHistoNGoodESDTracksVsNGammaCandidates[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks(),fGammaCandidates->GetEntries(), fWeightCentrality[iCut]);
      } else {
	fHistoNGammaCandidates[iCut]->Fill(fGammaCandidates->GetEntries(),fWeightJetJetMC);
	if( fIsMC < 2 ) fHistoNGoodESDTracksVsNGammaCandidates[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks(),fGammaCandidates->GetEntries());
      }
    }

    if(fDoMesonAnalysis){ // Meson Analysis
      if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseMCPSmearing() && fIsMC > 0 ){
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

      CalculatePi0Candidates(); // Combine Gammas
      if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->DoBGCalculation()){
        if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->BackgroundHandlerType() == 0){
          CalculateBackground(); // Combinatorial Background
          UpdateEventByEventData(); // Store Event for mixed Events
        } else {
          CalculateBackgroundRP(); // Combinatorial Background
          fBGHandlerRP[iCut]->AddEvent(fGammaCandidates,fInputEvent); // Store Event for mixed Events
        }
      }
      if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseMCPSmearing() && fIsMC > 0 ){
        for(Int_t gamma=0;gamma<fGammaCandidates->GetEntries();gamma++){ // Smear the AODPhotons in MC
          ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->SetPx(fUnsmearedPx[gamma]); // Reset Unsmeared Momenta
          ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->SetPy(fUnsmearedPy[gamma]);
          ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->SetPz(fUnsmearedPz[gamma]);
          ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->SetE(fUnsmearedE[gamma]);
        }
        delete[] fUnsmearedPx; fUnsmearedPx = 0x0;
        delete[] fUnsmearedPy; fUnsmearedPy = 0x0;
        delete[] fUnsmearedPz; fUnsmearedPz = 0x0;
        delete[] fUnsmearedE;  fUnsmearedE  = 0x0;
      }

      if( fIsMC > 0 ){
        vecDoubleCountTruePi0s.clear();
        vecDoubleCountTrueEtas.clear();
        FillMultipleCountHistoAndClear(mapMultipleCountTruePi0s,fHistoMultipleCountTruePi0[iCut]);
        FillMultipleCountHistoAndClear(mapMultipleCountTrueEtas,fHistoMultipleCountTrueEta[iCut]);
      }
    }

    if( fIsMC > 0 ){
      vecDoubleCountTrueConvGammas.clear();
      FillMultipleCountHistoAndClear(mapMultipleCountTrueConvGammas,fHistoMultipleCountTrueConvGamma[iCut]);
    }

    fGammaCandidates->Clear(); // delete this cuts good gammas
  }

  if( fIsMC > 0 && fInputEvent->IsA()==AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled())){
    RelabelAODPhotonCandidates(kFALSE); // Back to ESDMC Label
    fV0Reader->RelabelAODs(kFALSE);
  }

  PostData(1, fOutputContainer);
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::ProcessPhotonCandidates()
{
  Int_t nV0 = 0;
  TList *GammaCandidatesStepOne = new TList();
  TList *GammaCandidatesStepTwo = new TList();
  // Loop over Photon Candidates allocated by ReaderV1
  for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){
    AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(i);
    if(!PhotonCandidate) continue;
    fIsFromSelectedHeader = kTRUE;

    Float_t weightMatBudgetGamma = 1.;
    if (fDoMaterialBudgetWeightingOfGammasForTrueMesons && ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetMaterialBudgetWeightsInitialized()) {
      weightMatBudgetGamma = ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetMaterialBudgetCorrectingWeightForTrueGamma(PhotonCandidate);
    }
    if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetDoElecDeDxPostCalibration()){
      if(!(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->LoadElecDeDxPostCalibration(fInputEvent->GetRunNumber()))){
        AliFatal(Form("ERROR: LoadElecDeDxPostCalibration returned kFALSE for %d despite being requested!",fInputEvent->GetRunNumber()));
      }
    }

    if( fIsMC > 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
      Int_t isPosFromMBHeader
        = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCEvent, fInputEvent);
      if(isPosFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
      Int_t isNegFromMBHeader
        = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCEvent, fInputEvent);
      if(isNegFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;

      if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromSelectedHeader = kFALSE;
    }


    if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelected(PhotonCandidate,fInputEvent)) continue;
    if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(PhotonCandidate->GetPhotonPhi(),fEventPlaneAngle)) continue;
    if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseElecSharingCut() &&
      !((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){
      fGammaCandidates->Add(PhotonCandidate); // if no second loop is required add to events good gammas

      if(fIsFromSelectedHeader){
        if(fDoCentralityFlat > 0) fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt(), fWeightCentrality[fiCut]*fWeightJetJetMC*weightMatBudgetGamma);
        else fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt(), fWeightJetJetMC*weightMatBudgetGamma);
        if (fDoPhotonQA > 0 && fIsMC < 2){
          if(fDoCentralityFlat > 0){
            fHistoConvGammaPsiPairPt[fiCut]->Fill(PhotonCandidate->GetPsiPair(),PhotonCandidate->Pt(), fWeightCentrality[fiCut]*fWeightJetJetMC*weightMatBudgetGamma);
            fHistoConvGammaR[fiCut]->Fill(PhotonCandidate->GetConversionRadius(), fWeightCentrality[fiCut]*fWeightJetJetMC*weightMatBudgetGamma);
            fHistoConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta(), fWeightCentrality[fiCut]*fWeightJetJetMC*weightMatBudgetGamma);
            fHistoConvGammaPhi[fiCut]->Fill(PhotonCandidate->Phi(), fWeightCentrality[fiCut]*fWeightJetJetMC*weightMatBudgetGamma);
            if ((fDoPhotonQA == 4)||(fDoPhotonQA == 5)){
              fHistoConvGammaInvMass[fiCut]->Fill(PhotonCandidate->GetMass(), fWeightCentrality[fiCut]*fWeightJetJetMC*weightMatBudgetGamma);
              fHistoConvGammaInvMassReco[fiCut]->Fill(GetOriginalInvMass(PhotonCandidate,fInputEvent), fWeightCentrality[fiCut]*fWeightJetJetMC*weightMatBudgetGamma);
            }
          } else {
            fHistoConvGammaPsiPairPt[fiCut]->Fill(PhotonCandidate->GetPsiPair(),PhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
            fHistoConvGammaR[fiCut]->Fill(PhotonCandidate->GetConversionRadius(),fWeightJetJetMC*weightMatBudgetGamma);
            fHistoConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta(),fWeightJetJetMC*weightMatBudgetGamma);
            fHistoConvGammaPhi[fiCut]->Fill(PhotonCandidate->Phi(),fWeightJetJetMC*weightMatBudgetGamma);
            if ((fDoPhotonQA == 4)||(fDoPhotonQA == 5)){
              fHistoConvGammaInvMass[fiCut]->Fill(PhotonCandidate->GetMass(),fWeightJetJetMC*weightMatBudgetGamma);
              fHistoConvGammaInvMassReco[fiCut]->Fill(GetOriginalInvMass(PhotonCandidate,fInputEvent),fWeightJetJetMC*weightMatBudgetGamma);
            }
          }
        }
        if( fIsMC > 0 ){
          if(fInputEvent->IsA()==AliESDEvent::Class())
          ProcessTruePhotonCandidates(PhotonCandidate);
          if(fInputEvent->IsA()==AliAODEvent::Class())
          ProcessTruePhotonCandidatesAOD(PhotonCandidate);
        }
        if ((fDoPhotonQA == 2)||(fDoPhotonQA == 5)){
          if (fIsHeavyIon == 1 && PhotonCandidate->Pt() > 0.399 && PhotonCandidate->Pt() < 12.){
            fPtGamma = PhotonCandidate->Pt();
            fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
            fRConvPhoton = PhotonCandidate->GetConversionRadius();
            fEtaPhoton = PhotonCandidate->GetPhotonEta();
            iCatPhoton = PhotonCandidate->GetPhotonQuality();
            tESDConvGammaPtDcazCat[fiCut]->Fill();
            } else if ( ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetSingleElectronPtCut() < 0.04 && PhotonCandidate->Pt() > 0.099 && PhotonCandidate->Pt() < 16.){
            fPtGamma = PhotonCandidate->Pt();
            fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
            fRConvPhoton = PhotonCandidate->GetConversionRadius();
            fEtaPhoton = PhotonCandidate->GetPhotonEta();
            iCatPhoton = PhotonCandidate->GetPhotonQuality();
            tESDConvGammaPtDcazCat[fiCut]->Fill();
          } else if ( PhotonCandidate->Pt() > 0.299 && PhotonCandidate->Pt() < 16.){
            fPtGamma = PhotonCandidate->Pt();
            fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
            fRConvPhoton = PhotonCandidate->GetConversionRadius();
            fEtaPhoton = PhotonCandidate->GetPhotonEta();
            iCatPhoton = PhotonCandidate->GetPhotonQuality();
            tESDConvGammaPtDcazCat[fiCut]->Fill();
          }
        }
      }
    } else if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseElecSharingCut()){ // if Shared Electron cut is enabled, Fill array, add to step one
      ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->FillElectonLabelArray(PhotonCandidate,nV0);
      nV0++;
      GammaCandidatesStepOne->Add(PhotonCandidate);
    } else if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseElecSharingCut() &&
        ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){ // shared electron is disabled, step one not needed -> step two
      GammaCandidatesStepTwo->Add(PhotonCandidate);
    }
  }
  if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseElecSharingCut()){
    for(Int_t i = 0;i<GammaCandidatesStepOne->GetEntries();i++){
      AliAODConversionPhoton *PhotonCandidate= (AliAODConversionPhoton*) GammaCandidatesStepOne->At(i);
      if(!PhotonCandidate) continue;
      fIsFromSelectedHeader = kTRUE;

      Float_t weightMatBudgetGamma = 1.;
      if (fDoMaterialBudgetWeightingOfGammasForTrueMesons && ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetMaterialBudgetWeightsInitialized()) {
	weightMatBudgetGamma = ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetMaterialBudgetCorrectingWeightForTrueGamma(PhotonCandidate);
      }


     if(fMCEvent && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        Int_t isPosFromMBHeader
        = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCEvent, fInputEvent);
        Int_t isNegFromMBHeader
        = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCEvent, fInputEvent);
        if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromSelectedHeader = kFALSE;
      }
      if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->RejectSharedElectronV0s(PhotonCandidate,i,GammaCandidatesStepOne->GetEntries())) continue;
      if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){ // To Colse v0s cut diabled, step two not needed
        fGammaCandidates->Add(PhotonCandidate);

        if(fIsFromSelectedHeader){
          if(fDoCentralityFlat > 0) fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt(), fWeightCentrality[fiCut]*fWeightJetJetMC*weightMatBudgetGamma);
          else fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
          if (fDoPhotonQA > 0 && fIsMC < 2){
            if(fDoCentralityFlat > 0){
              fHistoConvGammaPsiPairPt[fiCut]->Fill(PhotonCandidate->GetPsiPair(),PhotonCandidate->Pt(), fWeightCentrality[fiCut]*fWeightJetJetMC*weightMatBudgetGamma);
              fHistoConvGammaR[fiCut]->Fill(PhotonCandidate->GetConversionRadius(), fWeightCentrality[fiCut]*fWeightJetJetMC*weightMatBudgetGamma);
              fHistoConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta(), fWeightCentrality[fiCut]*fWeightJetJetMC*weightMatBudgetGamma);
              fHistoConvGammaPhi[fiCut]->Fill(PhotonCandidate->Phi(), fWeightCentrality[fiCut]*fWeightJetJetMC*weightMatBudgetGamma);
              if ((fDoPhotonQA == 4)||(fDoPhotonQA == 5)){
                fHistoConvGammaInvMass[fiCut]->Fill(PhotonCandidate->GetMass(), fWeightCentrality[fiCut]*fWeightJetJetMC*weightMatBudgetGamma);
                fHistoConvGammaInvMassReco[fiCut]->Fill(GetOriginalInvMass(PhotonCandidate,fInputEvent), fWeightCentrality[fiCut]*fWeightJetJetMC*weightMatBudgetGamma);
              }
            } else {
              fHistoConvGammaPsiPairPt[fiCut]->Fill(PhotonCandidate->GetPsiPair(),PhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
              fHistoConvGammaR[fiCut]->Fill(PhotonCandidate->GetConversionRadius(),fWeightJetJetMC*weightMatBudgetGamma);
              fHistoConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta(),fWeightJetJetMC*weightMatBudgetGamma);
              fHistoConvGammaPhi[fiCut]->Fill(PhotonCandidate->Phi(),fWeightJetJetMC*weightMatBudgetGamma);
              if ((fDoPhotonQA == 4)||(fDoPhotonQA == 5)){
                fHistoConvGammaInvMass[fiCut]->Fill(PhotonCandidate->GetMass(),fWeightJetJetMC*weightMatBudgetGamma);
                fHistoConvGammaInvMassReco[fiCut]->Fill(GetOriginalInvMass(PhotonCandidate,fInputEvent),fWeightJetJetMC*weightMatBudgetGamma);
              }
            }
          }
          if( fIsMC > 0 ){
            if(fInputEvent->IsA()==AliESDEvent::Class())
              ProcessTruePhotonCandidates(PhotonCandidate);
            if(fInputEvent->IsA()==AliAODEvent::Class())
              ProcessTruePhotonCandidatesAOD(PhotonCandidate);
          }
          if ((fDoPhotonQA == 2)||(fDoPhotonQA == 5)){
            if (fIsHeavyIon ==1 && PhotonCandidate->Pt() > 0.399 && PhotonCandidate->Pt() < 12.){
              fPtGamma = PhotonCandidate->Pt();
              fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
              fRConvPhoton = PhotonCandidate->GetConversionRadius();
              fEtaPhoton = PhotonCandidate->GetPhotonEta();
              iCatPhoton = PhotonCandidate->GetPhotonQuality();
              tESDConvGammaPtDcazCat[fiCut]->Fill();
            } else if ( ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetSingleElectronPtCut() < 0.04 && PhotonCandidate->Pt() > 0.099 && PhotonCandidate->Pt() < 16.){
              fPtGamma = PhotonCandidate->Pt();
              fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
              fRConvPhoton = PhotonCandidate->GetConversionRadius();
              fEtaPhoton = PhotonCandidate->GetPhotonEta();
              iCatPhoton = PhotonCandidate->GetPhotonQuality();
              tESDConvGammaPtDcazCat[fiCut]->Fill();
            } else if ( PhotonCandidate->Pt() > 0.299 && PhotonCandidate->Pt() < 16.){
              fPtGamma = PhotonCandidate->Pt();
              fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
              fRConvPhoton = PhotonCandidate->GetConversionRadius();
              fEtaPhoton = PhotonCandidate->GetPhotonEta();
              iCatPhoton = PhotonCandidate->GetPhotonQuality();
              tESDConvGammaPtDcazCat[fiCut]->Fill();
            }
          }
        }
      } else GammaCandidatesStepTwo->Add(PhotonCandidate); // Close v0s cut enabled -> add to list two
    }
  }
  if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){
    for(Int_t i = 0;i<GammaCandidatesStepTwo->GetEntries();i++){
      AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) GammaCandidatesStepTwo->At(i);
      if(!PhotonCandidate) continue;
      fIsFromSelectedHeader = kTRUE;

     Float_t weightMatBudgetGamma = 1.;
      if (fDoMaterialBudgetWeightingOfGammasForTrueMesons && ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetMaterialBudgetWeightsInitialized()) {
	weightMatBudgetGamma = ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetMaterialBudgetCorrectingWeightForTrueGamma(PhotonCandidate);
      }

      if(fMCEvent && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        Int_t isPosFromMBHeader
        = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCEvent, fInputEvent);
        Int_t isNegFromMBHeader
        = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCEvent, fInputEvent);
        if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromSelectedHeader = kFALSE;
      }
      if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->RejectToCloseV0s(PhotonCandidate,GammaCandidatesStepTwo,i)) continue;
      fGammaCandidates->Add(PhotonCandidate); // Add gamma to current cut TList

      if(fIsFromSelectedHeader){
        if(fDoCentralityFlat > 0) fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt(), fWeightCentrality[fiCut]*fWeightJetJetMC*weightMatBudgetGamma);
        else fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
        if (fDoPhotonQA > 0 && fIsMC < 2 ){
          if(fDoCentralityFlat > 0){
            fHistoConvGammaPsiPairPt[fiCut]->Fill(PhotonCandidate->GetPsiPair(),PhotonCandidate->Pt(), fWeightCentrality[fiCut]*fWeightJetJetMC*weightMatBudgetGamma);
            fHistoConvGammaR[fiCut]->Fill(PhotonCandidate->GetConversionRadius(), fWeightCentrality[fiCut]*fWeightJetJetMC*weightMatBudgetGamma);
            fHistoConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta(), fWeightCentrality[fiCut]*fWeightJetJetMC*weightMatBudgetGamma);
            fHistoConvGammaPhi[fiCut]->Fill(PhotonCandidate->Phi(), fWeightCentrality[fiCut]*fWeightJetJetMC*weightMatBudgetGamma);
            if ((fDoPhotonQA == 4)||(fDoPhotonQA == 5)){
              fHistoConvGammaInvMass[fiCut]->Fill(PhotonCandidate->GetMass(), fWeightCentrality[fiCut]*fWeightJetJetMC*weightMatBudgetGamma);
              fHistoConvGammaInvMassReco[fiCut]->Fill(GetOriginalInvMass(PhotonCandidate,fInputEvent), fWeightCentrality[fiCut]*fWeightJetJetMC*weightMatBudgetGamma);
            }
          } else {
            fHistoConvGammaPsiPairPt[fiCut]->Fill(PhotonCandidate->GetPsiPair(),PhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
            fHistoConvGammaR[fiCut]->Fill(PhotonCandidate->GetConversionRadius(),fWeightJetJetMC*weightMatBudgetGamma);
            fHistoConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta(),fWeightJetJetMC*weightMatBudgetGamma);
            fHistoConvGammaPhi[fiCut]->Fill(PhotonCandidate->Phi(),fWeightJetJetMC*weightMatBudgetGamma);
            if ((fDoPhotonQA == 4)||(fDoPhotonQA == 5)){
              fHistoConvGammaInvMass[fiCut]->Fill(PhotonCandidate->GetMass(),fWeightJetJetMC*weightMatBudgetGamma);
              fHistoConvGammaInvMassReco[fiCut]->Fill(GetOriginalInvMass(PhotonCandidate,fInputEvent),fWeightJetJetMC*weightMatBudgetGamma);
            }
          }
        }
        if( fIsMC > 0 ){
          if(fInputEvent->IsA()==AliESDEvent::Class())
            ProcessTruePhotonCandidates(PhotonCandidate);
          if(fInputEvent->IsA()==AliAODEvent::Class())
            ProcessTruePhotonCandidatesAOD(PhotonCandidate);
        }
        if ((fDoPhotonQA == 2)||(fDoPhotonQA == 5)){
          if (fIsHeavyIon == 1 && PhotonCandidate->Pt() > 0.399 && PhotonCandidate->Pt() < 12.){
            fPtGamma = PhotonCandidate->Pt();
            fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
            fRConvPhoton = PhotonCandidate->GetConversionRadius();
            fEtaPhoton = PhotonCandidate->GetPhotonEta();
            iCatPhoton = PhotonCandidate->GetPhotonQuality();
            tESDConvGammaPtDcazCat[fiCut]->Fill();
          } else if ( ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetSingleElectronPtCut() < 0.04 && PhotonCandidate->Pt() > 0.099 && PhotonCandidate->Pt() < 16.){
            fPtGamma = PhotonCandidate->Pt();
            fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
            fRConvPhoton = PhotonCandidate->GetConversionRadius();
            fEtaPhoton = PhotonCandidate->GetPhotonEta();
            iCatPhoton = PhotonCandidate->GetPhotonQuality();
            tESDConvGammaPtDcazCat[fiCut]->Fill();
          } else if ( PhotonCandidate->Pt() > 0.299 && PhotonCandidate->Pt() < 16.){
            fPtGamma = PhotonCandidate->Pt();
            fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
            fRConvPhoton = PhotonCandidate->GetConversionRadius();
            fEtaPhoton = PhotonCandidate->GetPhotonEta();
            iCatPhoton = PhotonCandidate->GetPhotonQuality();
            tESDConvGammaPtDcazCat[fiCut]->Fill();
          }
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
void AliAnalysisTaskGammaConvV1::ProcessJets()
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
          if(!fDoLightOutput || (fDoLightOutput && fDoJetQA)) fHistoTruevsRecJetPt[fiCut]->Fill(fVectorJetPt.at(i), fTrueVectorJetPt.at(match));
          if(fDoJetQA){
            if(fVectorJetPt.at(i) >= 10) fHistoEventwJets[fiCut]->Fill(0);
            if(fVectorJetPt.at(i) < 10 && fTrueVectorJetPt.at(match) >= 10) fHistoEventwJets[fiCut]->Fill(1);
            if(fVectorJetPt.at(i) >= 10 && fTrueVectorJetPt.at(match) < 10) fHistoEventwJets[fiCut]->Fill(2);
          }
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
void AliAnalysisTaskGammaConvV1::InitializeBDT()
{

    // This loads the library
   TMVA::Tools::Instance();
   fBDTreader = new TMVA::Reader( "!Color:!Silent" );

   fBDTvariable  = new Float_t[18];
   for(Int_t i=0; i<18; i++){
     fBDTvariable[i] = 0;
   }

   // add variables
   fBDTreader->AddVariable( "dEdxElectronITS + dEdxPositronITS" , &fBDTvariable[0] );
   fBDTreader->AddVariable( "photonQt"                          , &fBDTvariable[1] );
   fBDTreader->AddVariable( "photonInvMass"                     , &fBDTvariable[2] );
   fBDTreader->AddVariable( "photonR"                           , &fBDTvariable[3] );
   fBDTreader->AddVariable( "photonAlpha"                       , &fBDTvariable[4] );
   fBDTreader->AddVariable( "photonPsiPair"                     , &fBDTvariable[5] );
   fBDTreader->AddVariable( "photonCosPoint"                    , &fBDTvariable[6] );
   fBDTreader->AddVariable( "fracClsTPCPositron"                , &fBDTvariable[7] );
   fBDTreader->AddVariable( "fracClsTPCElectron"                , &fBDTvariable[8] );
   fBDTreader->AddVariable( "clsITSPositron"                    , &fBDTvariable[9] );
   fBDTreader->AddVariable( "clsITSElectron"                    , &fBDTvariable[10] );
   fBDTreader->AddVariable( "nSigmaTPCElectron"                 , &fBDTvariable[11] );
   fBDTreader->AddVariable( "nSigmaTPCPositron"                 , &fBDTvariable[12] );
   // add spectators like in the training
   fBDTreader->AddSpectator( "photonPt", &fBDTvariable[13] );
   fBDTreader->AddSpectator( "kind", &fBDTvariable[14] );
   fBDTreader->AddSpectator( "photonEta", &fBDTvariable[15] );
   fBDTreader->AddSpectator( "photonP", &fBDTvariable[16] );
   fBDTreader->AddSpectator( "photonPhi", &fBDTvariable[17] );

   fBDTreader->BookMVA( "BDT method", fFileNameBDT.Data());

}
//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::ProcessPhotonBDT()
{
  // Loop over Photon Candidates allocated by ReaderV1
  for(Int_t i = 0; i < fGammaCandidates->GetEntries(); i++){
    AliAODConversionPhoton *PhotonCandidate=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(i));
    if(!PhotonCandidate) continue;

    Float_t SingleTrackBDTVariableValues[7];
    if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetBDTVariableValues(PhotonCandidate,fInputEvent, SingleTrackBDTVariableValues)){
      fBDTvariable[0] = SingleTrackBDTVariableValues[0]; //"dEdxElectronITS + dEdxPositronITS"
      fBDTvariable[1] = PhotonCandidate->GetArmenterosQt(); //"photonQt"
      fBDTvariable[2] = GetOriginalInvMass(PhotonCandidate,fInputEvent); //"photonInvMass"
      fBDTvariable[3] = PhotonCandidate->GetConversionRadius(); //"photonR"
      fBDTvariable[4] = PhotonCandidate->GetArmenterosAlpha(); //"photonAlpha"
      fBDTvariable[5] = PhotonCandidate->GetPsiPair(); //"photonPsiPair"
      fBDTvariable[6] = ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetCosineOfPointingAngle(PhotonCandidate,fInputEvent); //"photonCosPoint"
      fBDTvariable[7] = SingleTrackBDTVariableValues[1]; //"fracClsTPCPositron"
      fBDTvariable[8] = SingleTrackBDTVariableValues[2]; //"fracClsTPCElectron"
      fBDTvariable[9] = SingleTrackBDTVariableValues[3]; //"clsITSPositron"
      fBDTvariable[10] = SingleTrackBDTVariableValues[4]; //"clsITSElectron"
      fBDTvariable[11] = SingleTrackBDTVariableValues[5]; //"nSigmaTPCElectron"
      fBDTvariable[12] = SingleTrackBDTVariableValues[6]; //"nSigmaTPCPositron"

      //Here we want to see the BDT output value
      Float_t BDToutput = fBDTreader->EvaluateMVA("BDT method");
      fHistoBDToutput[fiCut]->Fill( BDToutput );
      if(BDToutput > 0.) fHistoBDToutputPt[fiCut]->Fill( PhotonCandidate->Pt() );

      //going to check the BDT output if the current photon is a true conversion
      if(fIsMC > 0){
        TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
        if (AODMCTrackArray != NULL && PhotonCandidate != NULL){
          AliAODMCParticle *posDaughter = (AliAODMCParticle*) AODMCTrackArray->At(PhotonCandidate->GetMCLabelPositive());
          AliAODMCParticle *negDaughter = (AliAODMCParticle*) AODMCTrackArray->At(PhotonCandidate->GetMCLabelNegative());
          if(posDaughter != NULL && negDaughter != NULL){
              Int_t pdgCode[2] = {TMath::Abs(posDaughter->GetPdgCode()),TMath::Abs(negDaughter->GetPdgCode())};
              if(posDaughter->GetMother() == negDaughter->GetMother()){
                if((pdgCode[0]==11 && pdgCode[1]==11)&&  posDaughter->GetPdgCode()!=negDaughter->GetPdgCode()){
                  AliAODMCParticle *Photon = (AliAODMCParticle*) AODMCTrackArray->At(posDaughter->GetMother());
                  if(Photon->GetPdgCode() == 22) fHistoBDToutputMCTrue[fiCut]->Fill( BDToutput );
                }
              }
          }
        }
      }

    }
  }
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::ProcessTruePhotonCandidatesAOD(AliAODConversionPhoton *TruePhotonCandidate)
{
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  Double_t magField = fInputEvent->GetMagneticField();
  if( magField  < 0.0 ){
    magField =  1.0;
  }
  else {
    magField =  -1.0;
  }

  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (AODMCTrackArray != NULL && TruePhotonCandidate != NULL){

    AliAODMCParticle *posDaughter = (AliAODMCParticle*) AODMCTrackArray->At(TruePhotonCandidate->GetMCLabelPositive());
    AliAODMCParticle *negDaughter = (AliAODMCParticle*) AODMCTrackArray->At(TruePhotonCandidate->GetMCLabelNegative());
    iPhotonMCInfo = 0;

    if(posDaughter == NULL || negDaughter == NULL) return; // One particle does not exist
    Int_t pdgCode[2] = {TMath::Abs(posDaughter->GetPdgCode()),TMath::Abs(negDaughter->GetPdgCode())};

    Double_t PhiParticle[2] = {posDaughter->Phi(),negDaughter->Phi()};

    if(posDaughter->GetMother() != negDaughter->GetMother()){
      FillPhotonCombinatorialBackgroundHist(TruePhotonCandidate, pdgCode, PhiParticle);
      if(fDoPhotonQA == 3 && fIsHeavyIon == 1){
          if(posDaughter->GetMother() > -1){ //contamination is not a primary
              AliAODMCParticle *Mom = (AliAODMCParticle*) AODMCTrackArray->At(posDaughter->GetMother());
              if(Mom->GetMother() == -1){
                  FillPhotonCombinatorialMothersHistAOD(posDaughter,Mom);
              } else {
                  AliAODMCParticle *GranMom = (AliAODMCParticle*) AODMCTrackArray->At(Mom->GetMother());
                  if(GranMom->GetMother() == -1){
                      FillPhotonCombinatorialMothersHistAOD(posDaughter,GranMom);
                  } else {
                      AliAODMCParticle *GranGranMom = (AliAODMCParticle*) AODMCTrackArray->At(GranMom->GetMother());
                      if(GranGranMom->GetMother() == -1){
                          FillPhotonCombinatorialMothersHistAOD(posDaughter,GranGranMom);
                      } else FillPhotonCombinatorialMothersHistAOD(posDaughter,GranGranMom);
                  }
              }
          }
          if(negDaughter->GetMother() > -1){ //contamination is not a primary
              AliAODMCParticle *Mom = (AliAODMCParticle*) AODMCTrackArray->At(negDaughter->GetMother());
              if(Mom->GetMother() == -1){
                  FillPhotonCombinatorialMothersHistAOD(negDaughter,Mom);
              } else {
                  AliAODMCParticle *GranMom = (AliAODMCParticle*) AODMCTrackArray->At(Mom->GetMother());
                  if(GranMom->GetMother() == -1){
                      FillPhotonCombinatorialMothersHistAOD(negDaughter,GranMom);
                  } else {
                      AliAODMCParticle *GranGranMom = (AliAODMCParticle*) AODMCTrackArray->At(GranMom->GetMother());
                      if(GranGranMom->GetMother() == -1){
                          FillPhotonCombinatorialMothersHistAOD(negDaughter,GranGranMom);
                      } else FillPhotonCombinatorialMothersHistAOD(posDaughter,GranGranMom);
                  }
              }
          }
      }
      iPhotonMCInfo = 1;
      return;
    }
    else if(posDaughter->GetMother() == -1){
      FillPhotonCombinatorialBackgroundHist(TruePhotonCandidate, pdgCode, PhiParticle);
      if(fDoPhotonQA == 3 && fIsHeavyIon == 1){
          FillPhotonCombinatorialMothersHistAOD(posDaughter,posDaughter);
          FillPhotonCombinatorialMothersHistAOD(negDaughter,negDaughter);
      }
      iPhotonMCInfo = 1;
      return;
    }

    if(pdgCode[0]!=11 || pdgCode[1]!=11){
      iPhotonMCInfo = 1;
      return; //One Particle is not a electron
    }
    if(posDaughter->GetPdgCode()==negDaughter->GetPdgCode()){
      iPhotonMCInfo = 1;
      return; // Same Charge
    }

    AliAODMCParticle *Photon = (AliAODMCParticle*) AODMCTrackArray->At(posDaughter->GetMother());
    AliVTrack * electronCandidate = ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetTrack(fInputEvent,TruePhotonCandidate->GetTrackLabelNegative() );
    AliVTrack * positronCandidate = ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetTrack(fInputEvent,TruePhotonCandidate->GetTrackLabelPositive() );
    Double_t deltaPhi = magField * TVector2::Phi_mpi_pi( electronCandidate->Phi()-positronCandidate->Phi());

    if(Photon->GetPdgCode() != 22){
      fHistoTrueDalitzPsiPairDeltaPhi[fiCut]->Fill(deltaPhi,TruePhotonCandidate->GetPsiPair(),fWeightJetJetMC);
      iPhotonMCInfo = 1;
      return; // Mother is no Photon
    }

    if(((posDaughter->GetMCProcessCode())) != 5 || ((negDaughter->GetMCProcessCode())) != 5){
      iPhotonMCInfo = 1;
      return;// check if the daughters come from a conversion
    }
    // STILL A BUG IN ALIROOT >>8 HAS TPO BE REMOVED AFTER FIX

    Double_t rConv=0.;
    rConv = sqrt( (posDaughter->Xv()*posDaughter->Xv()) + (posDaughter->Yv()*posDaughter->Yv()) );

    // True Photon

    Float_t weightMatBudgetGamma = 1.;
    if (fDoMaterialBudgetWeightingOfGammasForTrueMesons && ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetMaterialBudgetWeightsInitialized()) {
      weightMatBudgetGamma = ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetMaterialBudgetCorrectingWeightForTrueGamma(TruePhotonCandidate);
    }

    fHistoTrueConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
    if (CheckVectorForDoubleCount(vecDoubleCountTrueConvGammas,posDaughter->GetMother())){
      fHistoDoubleCountTrueConvGammaRPt[fiCut]->Fill(TruePhotonCandidate->GetConversionRadius(),TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
      FillMultipleCountMap(mapMultipleCountTrueConvGammas,posDaughter->GetMother());
    }
    if (fDoPhotonQA > 0 ){
      if (fIsMC < 2){
        fHistoTrueConvGammaPsiPairPt[fiCut]->Fill(TruePhotonCandidate->GetPsiPair(),TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
        fHistoTrueConvGammaEta[fiCut]->Fill(TruePhotonCandidate->Eta(),fWeightJetJetMC*weightMatBudgetGamma);
        fHistoTrueConvGammaR[fiCut]->Fill(TruePhotonCandidate->GetConversionRadius(),fWeightJetJetMC*weightMatBudgetGamma);
        fHistoTrueConvGammaRMC[fiCut]->Fill(rConv,fWeightJetJetMC*weightMatBudgetGamma);
      }
      fHistoTrueConvGammaPtMC[fiCut]->Fill(Photon->Pt(), fWeightJetJetMC*weightMatBudgetGamma);
    }
    fHistoTrueGammaPsiPairDeltaPhi[fiCut]->Fill(deltaPhi,TruePhotonCandidate->GetPsiPair(),fWeightJetJetMC*weightMatBudgetGamma);

    Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, Photon, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
    if(isPrimary){
      // Count just primary MC Gammas as true --> For Ratio esdtruegamma / mcconvgamma
      iPhotonMCInfo = 6;
      fHistoTruePrimaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
      fHistoTruePrimaryConvGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt(),fWeightJetJetMC*weightMatBudgetGamma); // Allways Filled
      // (Not Filled for i6, Extra Signal Gamma (parambox) are secondary)
    } else {
      iPhotonMCInfo = 2;
      if(((AliAODMCParticle*)AODMCTrackArray->At(Photon->GetMother()))->GetMother() > -1 && ((AliAODMCParticle*)AODMCTrackArray->At(((AliAODMCParticle*)AODMCTrackArray->At(Photon->GetMother()))->GetMother()))->GetMother() > -1 ){
        if (((AliAODMCParticle*)AODMCTrackArray->At(((AliAODMCParticle*)AODMCTrackArray->At(Photon->GetMother()))->GetMother()))->GetPdgCode() == 310){
            iPhotonMCInfo = 4;
            fHistoTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),0.,fWeightJetJetMC*weightMatBudgetGamma);
            fHistoTrueSecondaryConvGammaMCPt[fiCut]->Fill(Photon->Pt(),0.,fWeightJetJetMC*weightMatBudgetGamma);
            fHistoTrueSecondaryConvGammaFromXFromK0sMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
        } else if (((AliAODMCParticle*)AODMCTrackArray->At(((AliAODMCParticle*)AODMCTrackArray->At(Photon->GetMother()))->GetMother()))->GetPdgCode() == 130) {
            iPhotonMCInfo = 7;
            fHistoTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),1.,fWeightJetJetMC*weightMatBudgetGamma);
            fHistoTrueSecondaryConvGammaMCPt[fiCut]->Fill(Photon->Pt(),1.,fWeightJetJetMC*weightMatBudgetGamma);
            fHistoTrueSecondaryConvGammaFromXFromK0lMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
        } else if (((AliAODMCParticle*)AODMCTrackArray->At(((AliAODMCParticle*)AODMCTrackArray->At(Photon->GetMother()))->GetMother()))->GetPdgCode() == 3122) {
            iPhotonMCInfo = 5;
            fHistoTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),2.,fWeightJetJetMC*weightMatBudgetGamma);
            fHistoTrueSecondaryConvGammaMCPt[fiCut]->Fill(Photon->Pt(),2.,fWeightJetJetMC*weightMatBudgetGamma);
            fHistoTrueSecondaryConvGammaFromXFromLambdaMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
        } else if (((AliAODMCParticle*)AODMCTrackArray->At(((AliAODMCParticle*)AODMCTrackArray->At(Photon->GetMother()))->GetMother()))->GetPdgCode() == 221) {
            iPhotonMCInfo = 3;
            fHistoTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3.,fWeightJetJetMC*weightMatBudgetGamma);
            fHistoTrueSecondaryConvGammaMCPt[fiCut]->Fill(Photon->Pt(),3.,fWeightJetJetMC*weightMatBudgetGamma);
        } else {
    //         if ( !(TMath::Abs(((AliAODMCParticle*)AODMCTrackArray->At(Photon->GetMother()))->GetPdgCode()) == 11 && ((AliAODMCParticle*)AODMCTrackArray->At(((AliAODMCParticle*)AODMCTrackArray->At(Photon->GetMother()))->GetMother()))->GetPdgCode() == 22) ) {
            fHistoTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3.,fWeightJetJetMC*weightMatBudgetGamma);
            fHistoTrueSecondaryConvGammaMCPt[fiCut]->Fill(Photon->Pt(),3.,fWeightJetJetMC*weightMatBudgetGamma);
    //         }
        }
      }
    }
  }
  return;
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::ProcessTruePhotonCandidates(AliAODConversionPhoton *TruePhotonCandidate)
{
  Double_t magField = fInputEvent->GetMagneticField();
  if( magField  < 0.0 ){
    magField =  1.0;
  }
  else {
    magField =  -1.0;
  }

  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  // Process True Photons
  TParticle *posDaughter = TruePhotonCandidate->GetPositiveMCDaughter(fMCEvent);
  TParticle *negDaughter = TruePhotonCandidate->GetNegativeMCDaughter(fMCEvent);

  iPhotonMCInfo = 0;

  if(posDaughter == NULL || negDaughter == NULL) return; // One particle does not exist
  Int_t pdgCode[2] = {TMath::Abs(posDaughter->GetPdgCode()),TMath::Abs(negDaughter->GetPdgCode())};

  Double_t PhiParticle[2] = {posDaughter->Phi(),negDaughter->Phi()};

  iPhotonMCInfo = 1;
  if(posDaughter->GetMother(0) != negDaughter->GetMother(0)){
    FillPhotonCombinatorialBackgroundHist(TruePhotonCandidate, pdgCode, PhiParticle);
    if(fDoPhotonQA == 3 && fIsHeavyIon == 1){
        if(posDaughter->GetMother(0) > -1){ //contamination is not a primary
            TParticle *Mom = fMCEvent->Particle(posDaughter->GetMother(0));
            if(Mom->GetMother(0) == -1){
                FillPhotonCombinatorialMothersHistESD(posDaughter,Mom);
            } else {
                TParticle *GranMom = fMCEvent->Particle(Mom->GetMother(0));
                if(GranMom->GetMother(0) == -1){
                    FillPhotonCombinatorialMothersHistESD(posDaughter,GranMom);
                } else {
                    TParticle *GranGranMom = fMCEvent->Particle(GranMom->GetMother(0));
                    if(GranGranMom->GetMother(0) == -1){
                        FillPhotonCombinatorialMothersHistESD(posDaughter,GranGranMom);
                    } else FillPhotonCombinatorialMothersHistESD(posDaughter,GranGranMom);
                }
            }
        }
        if(negDaughter->GetMother(0) > -1){ //contamination is not a primary
            TParticle *Mom = fMCEvent->Particle(negDaughter->GetMother(0));
            if(Mom->GetMother(0) == -1){
                FillPhotonCombinatorialMothersHistESD(negDaughter,Mom);
            } else {
                TParticle *GranMom = fMCEvent->Particle(Mom->GetMother(0));
                if(GranMom->GetMother(0) == -1){
                    FillPhotonCombinatorialMothersHistESD(negDaughter,GranMom);
                } else {
                    TParticle *GranGranMom = fMCEvent->Particle(GranMom->GetMother(0));
                    if(GranGranMom->GetMother(0) == -1){
                        FillPhotonCombinatorialMothersHistESD(negDaughter,GranGranMom);
                    } else FillPhotonCombinatorialMothersHistESD(posDaughter,GranGranMom);
                }
            }
        }
    }
    return;
  } else if(posDaughter->GetMother(0) == -1){ //gamma contamination is a primary
    FillPhotonCombinatorialBackgroundHist(TruePhotonCandidate, pdgCode, PhiParticle);
    if(fDoPhotonQA == 3 && fIsHeavyIon == 1){
        FillPhotonCombinatorialMothersHistESD(posDaughter,posDaughter);
        FillPhotonCombinatorialMothersHistESD(negDaughter,negDaughter);
    }
    return;
  }

  if(pdgCode[0]!=11 || pdgCode[1]!=11) return; //One Particle is not a electron

  if(posDaughter->GetPdgCode()==negDaughter->GetPdgCode()) return; // Same Charge

  TParticle *Photon = TruePhotonCandidate->GetMCParticle(fMCEvent);
  AliVTrack * electronCandidate = ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetTrack(fInputEvent,TruePhotonCandidate->GetTrackLabelNegative() );
  AliVTrack * positronCandidate = ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetTrack(fInputEvent,TruePhotonCandidate->GetTrackLabelPositive() );
  Double_t deltaPhi = magField * TVector2::Phi_mpi_pi( electronCandidate->Phi()-positronCandidate->Phi());

  if(Photon->GetPdgCode() != 22){
    fHistoTrueDalitzPsiPairDeltaPhi[fiCut]->Fill(deltaPhi,TruePhotonCandidate->GetPsiPair(),fWeightJetJetMC);
    return; // Mother is no Photon
  }

  if(posDaughter->GetUniqueID() != 5 || negDaughter->GetUniqueID() !=5) return;// check if the daughters come from a conversion

  // True Photon
  Float_t weightMatBudgetGamma = 1.;
  if (fDoMaterialBudgetWeightingOfGammasForTrueMesons && ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetMaterialBudgetWeightsInitialized()) {
    weightMatBudgetGamma = ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetMaterialBudgetCorrectingWeightForTrueGamma(TruePhotonCandidate);
  }

  // cout<< " AM- Material Budget weight Gamma::"<< weightMatBudgetGamma << " "<< TruePhotonCandidate->GetConversionRadius() << endl;

  fHistoTrueConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
  if (CheckVectorForDoubleCount(vecDoubleCountTrueConvGammas,posDaughter->GetMother(0))){
    fHistoDoubleCountTrueConvGammaRPt[fiCut]->Fill(TruePhotonCandidate->GetConversionRadius(),TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
    FillMultipleCountMap(mapMultipleCountTrueConvGammas,posDaughter->GetMother(0));
  }
  if (fDoPhotonQA > 0){
    if (fIsMC < 2){
      fHistoTrueConvGammaPsiPairPt[fiCut]->Fill(TruePhotonCandidate->GetPsiPair(),TruePhotonCandidate->Pt(),weightMatBudgetGamma);
      fHistoTrueConvGammaEta[fiCut]->Fill(TruePhotonCandidate->Eta(),weightMatBudgetGamma);
      fHistoTrueConvGammaR[fiCut]->Fill(TruePhotonCandidate->GetConversionRadius(),weightMatBudgetGamma);
      fHistoTrueConvGammaRMC[fiCut]->Fill(posDaughter->R(),weightMatBudgetGamma);
      if ((fDoPhotonQA == 4)||(fDoPhotonQA == 5)){
        fHistoTrueConvGammaInvMass[fiCut]->Fill(TruePhotonCandidate->GetMass(),weightMatBudgetGamma);
        fHistoTrueConvGammaInvMassReco[fiCut]->Fill(GetOriginalInvMass(TruePhotonCandidate,fInputEvent),weightMatBudgetGamma);
      }
    }
    fHistoTrueConvGammaPtMC[fiCut]->Fill(Photon->Pt(), fWeightJetJetMC*weightMatBudgetGamma);
  }

  fHistoTrueGammaPsiPairDeltaPhi[fiCut]->Fill(deltaPhi,TruePhotonCandidate->GetPsiPair(),fWeightJetJetMC*weightMatBudgetGamma);
  if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, posDaughter->GetMother(0), mcProdVtxX, mcProdVtxY, mcProdVtxZ)){
    // filling primary histograms
    // Count just primary MC Gammas as true --> For Ratio esdtruegamma / mcconvgamma
    iPhotonMCInfo = 6;
    fHistoTruePrimaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
    fHistoTruePrimaryConvGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt(),fWeightJetJetMC*weightMatBudgetGamma); // Allways Filled
    // (Not Filled for i6, Extra Signal Gamma (parambox) are secondary)
  } else {
    // filling secondary photon histograms
    iPhotonMCInfo = 2;
    if( Photon->GetMother(0) > -1 && fMCEvent->Particle(Photon->GetMother(0))->GetMother(0) > -1){
      if (fMCEvent->Particle(fMCEvent->Particle(Photon->GetMother(0))->GetMother(0))->GetPdgCode() == 310){
        iPhotonMCInfo = 4;
        fHistoTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),0.,fWeightJetJetMC*weightMatBudgetGamma);
        fHistoTrueSecondaryConvGammaMCPt[fiCut]->Fill(Photon->Pt(),0.,fWeightJetJetMC*weightMatBudgetGamma);
        fHistoTrueSecondaryConvGammaFromXFromK0sMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
      } else if (fMCEvent->Particle(fMCEvent->Particle(Photon->GetMother(0))->GetMother(0))->GetPdgCode() == 130) {
        iPhotonMCInfo = 7;
        fHistoTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),1.,fWeightJetJetMC*weightMatBudgetGamma);
        fHistoTrueSecondaryConvGammaMCPt[fiCut]->Fill(Photon->Pt(),1.,fWeightJetJetMC*weightMatBudgetGamma);
        fHistoTrueSecondaryConvGammaFromXFromK0lMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
      } else if (fMCEvent->Particle(fMCEvent->Particle(Photon->GetMother(0))->GetMother(0))->GetPdgCode() == 3122) {
        iPhotonMCInfo = 5;
        fHistoTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),2.,fWeightJetJetMC*weightMatBudgetGamma);
        fHistoTrueSecondaryConvGammaMCPt[fiCut]->Fill(Photon->Pt(),2.,fWeightJetJetMC*weightMatBudgetGamma);
        fHistoTrueSecondaryConvGammaFromXFromLambdaMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
      } else if (fMCEvent->Particle(fMCEvent->Particle(Photon->GetMother(0))->GetMother(0))->GetPdgCode() == 221) {
        iPhotonMCInfo = 3;
        fHistoTrueSecondaryConvGammaPt[fiCut]->Fill(Photon->Pt(),3.,fWeightJetJetMC*weightMatBudgetGamma);
        fHistoTrueSecondaryConvGammaMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3.,fWeightJetJetMC*weightMatBudgetGamma);
      } else {
//         if ( !(TMath::Abs(fMCEvent->Particle(Photon->GetMother(0))->GetPdgCode()) == 11 && fMCEvent->Particle(fMCEvent->Particle(Photon->GetMother(0))->GetMother(0))->GetPdgCode() == 22) ) {
          fHistoTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3.,fWeightJetJetMC*weightMatBudgetGamma);
          fHistoTrueSecondaryConvGammaMCPt[fiCut]->Fill(Photon->Pt(),3.,fWeightJetJetMC*weightMatBudgetGamma);
//         }
      }
    } else {
      fHistoTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3.,fWeightJetJetMC*weightMatBudgetGamma);
      fHistoTrueSecondaryConvGammaMCPt[fiCut]->Fill(Photon->Pt(),3.,fWeightJetJetMC*weightMatBudgetGamma);
    }
  }
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::ProcessAODMCParticles()
{
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));

  if (AODMCTrackArray){
    // Loop over all primary MC particle
    for(Int_t i = 0; i < AODMCTrackArray->GetEntriesFast(); i++) {

      AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(i));
      if (!particle) continue;

      Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, particle, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
      if (isPrimary){

        Int_t isMCFromMBHeader = -1;
        if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
          isMCFromMBHeader
            = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
          if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
        }

        if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(particle->Phi(),fEventPlaneAngle,kFALSE)) continue;
        if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(particle,AODMCTrackArray,kFALSE)){
          fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // All MC Gamma
          if(particle->GetMother() >-1){ // Meson Decay Gamma
            switch((static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetMother())))->GetPdgCode()){
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
        if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(particle,AODMCTrackArray,kTRUE)){
          Double_t rConv = 0;
          for(Int_t daughterIndex=particle->GetDaughterLabel(0);daughterIndex<=particle->GetDaughterLabel(1);daughterIndex++){
            AliAODMCParticle *tmpDaughter = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(daughterIndex));
            if(!tmpDaughter) continue;
            if(TMath::Abs(tmpDaughter->GetPdgCode()) == 11){
              rConv = sqrt( (tmpDaughter->Xv()*tmpDaughter->Xv()) + (tmpDaughter->Yv()*tmpDaughter->Yv()) );
            }
          }
          fHistoMCConvGammaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
          if ( fDoPhotonQA > 0 && fIsMC < 2){
            fHistoMCConvGammaR[fiCut]->Fill(rConv,fWeightJetJetMC);
            fHistoMCConvGammaEta[fiCut]->Fill(particle->Eta(),fWeightJetJetMC);
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

          if (TMath::Abs(mesonY) < ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetRapidityCutValue()){
            if ( particle->GetPdgCode() == 211 ){  // positve pions
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(),0.,fWeightJetJetMC);
            } else if ( particle->GetPdgCode() == -211 ){  // negative pions
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(),1.,fWeightJetJetMC);
            } else if ( particle->GetPdgCode() == 321 ){  // positve kaons
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(),2.,fWeightJetJetMC);
            } else if ( particle->GetPdgCode() == -321 ){  // negative kaons
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(),3.,fWeightJetJetMC);
            } else if ( TMath::Abs(particle->GetPdgCode()) == 310 ){  // K0s
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(),4.,fWeightJetJetMC);
            } else if ( TMath::Abs(particle->GetPdgCode()) == 130 ){  // K0l
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(),5.,fWeightJetJetMC);
            } else if ( TMath::Abs(particle->GetPdgCode()) == 3122 ){  // Lambda/ AntiLambda
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(),6.,fWeightJetJetMC);
            } else if ( TMath::Abs(particle->GetPdgCode()) == 223 ){  // Omega
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(),7.,fWeightJetJetMC);
            } else if ( TMath::Abs(particle->GetPdgCode()) == 333 ){  // Phi
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(),8.,fWeightJetJetMC);
            } else if ( TMath::Abs(particle->GetPdgCode()) == 113 ){  // Rho0
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(),9.,fWeightJetJetMC);
            }
          }

          if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedAODMC(particle,AODMCTrackArray,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
            AliAODMCParticle* daughter0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetDaughterLabel(0)));
            AliAODMCParticle* daughter1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetDaughterLabel(1)));
            Float_t weighted= 1;
	    if (particle->Pt()>0.005){
	      weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(i, 0x0, fInputEvent);
	      //                   if(particle->GetPdgCode() == 221){
	      //                      cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
	      //                   }
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
              if ( fIsMC > 1 ) fHistoMCPi0WOEvtWeightPt[fiCut]->Fill(particle->Pt());
              if (fDoMesonQA > 0){
                if ( fIsMC < 2 )fHistoMCPi0PtY[fiCut]->Fill(particle->Pt(),mesonY,weighted); // All MC Pi0
                fHistoMCPi0PtAlpha[fiCut]->Fill(particle->Pt(),alpha,fWeightJetJetMC); // All MC Pi0
                if ( fIsMC == 2 ) fHistoMCPi0PtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),fWeightJetJetMC);
              }
            } else if(particle->GetPdgCode() == 221){
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
              if ( fIsMC > 1 ) fHistoMCEtaWOEvtWeightPt[fiCut]->Fill(particle->Pt());
              if (fDoMesonQA > 0){
                if ( fIsMC < 2 )fHistoMCEtaPtY[fiCut]->Fill(particle->Pt(),mesonY,weighted); // All MC Pi0
                fHistoMCEtaPtAlpha[fiCut]->Fill(particle->Pt(),alpha,fWeightJetJetMC); // All MC Pi0
                if ( fIsMC == 2 ) fHistoMCEtaPtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),fWeightJetJetMC);
              }
            }

            // Check the acceptance for both gammas
            if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(daughter0,AODMCTrackArray,kFALSE) &&
            ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(daughter1,AODMCTrackArray,kFALSE)  &&
            ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter0->Phi(),fEventPlaneAngle,kFALSE) &&
            ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter1->Phi(),fEventPlaneAngle,kFALSE)){

              if(particle->GetPdgCode() == 111){
                if(fIsMC > 1) fHistoMCPi0WOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Pi0 with gamma in acc NOT weighted at all
                fHistoMCPi0WOWeightInAccPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // MC Pi0 with gamma in acc NOT weighted
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
              } else if(particle->GetPdgCode() == 221){
                if(fIsMC > 1) fHistoMCEtaWOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Eta with gamma in acc NOT weighted at all
                fHistoMCEtaWOWeightInAccPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // MC Eta with gamma in acc NOT weighted
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
              }
            }
          }
        }
      // fill secondary histograms
      } else {
        Int_t isMCFromMBHeader = -1;
        if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
          isMCFromMBHeader
            = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
          if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
        }

        if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(particle->Phi(),fEventPlaneAngle,kFALSE)) {
          if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(particle,AODMCTrackArray,kFALSE)){
            if (particle->GetMother() > -1) {
              AliAODMCParticle *tmpMother = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetMother()));
              if (tmpMother->GetMother() > -1) {
                AliAODMCParticle *tmpGrandMother = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(tmpMother->GetMother()));
                if(tmpGrandMother->GetPdgCode() == 310) {
                  fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),0.,fWeightJetJetMC);
                } else if (tmpGrandMother->GetPdgCode() == 130) {
                  fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),1.,fWeightJetJetMC);
                } else if (tmpGrandMother->GetPdgCode() == 3122) {
                  fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),2.,fWeightJetJetMC);
                } else {
                  if ( !(TMath::Abs(tmpMother->GetPdgCode()) == 11 && tmpGrandMother->GetPdgCode() == 22) )
                    fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),3.,fWeightJetJetMC);
                }
              } else {
                fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),3.,fWeightJetJetMC);
              }
            } else {
              fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),3.,fWeightJetJetMC);
            }
          }

          if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(particle,AODMCTrackArray,kTRUE)){
            if (particle->GetMother() > -1) {
              AliAODMCParticle *tmpMother = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetMother()));
              if (tmpMother->GetMother() > -1) {
                AliAODMCParticle *tmpGrandMother = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(tmpMother->GetMother()));
                if(tmpGrandMother->GetPdgCode() == 310) {
                  fHistoMCSecondaryConvGammaPt[fiCut]->Fill(particle->Pt(),0.,fWeightJetJetMC);
                } else if (tmpGrandMother->GetPdgCode() == 130) {
                  fHistoMCSecondaryConvGammaPt[fiCut]->Fill(particle->Pt(),1.,fWeightJetJetMC);
                } else if (tmpGrandMother->GetPdgCode() == 3122) {
                  fHistoMCSecondaryConvGammaPt[fiCut]->Fill(particle->Pt(),2.,fWeightJetJetMC);
                } else {
                  if ( !(TMath::Abs(tmpMother->GetPdgCode()) == 11 && tmpGrandMother->GetPdgCode() == 22) )
                    fHistoMCSecondaryConvGammaPt[fiCut]->Fill(particle->Pt(),3.,fWeightJetJetMC);
                }
              } else {
                fHistoMCSecondaryConvGammaPt[fiCut]->Fill(particle->Pt(),3.,fWeightJetJetMC);
              }
            } else {
              fHistoMCSecondaryConvGammaPt[fiCut]->Fill(particle->Pt(),3.,fWeightJetJetMC);
            }
          }
        }

        if(fDoMesonAnalysis){
          if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedAODMC(particle,AODMCTrackArray,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
            AliAODMCParticle* daughter0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetDaughterLabel(0)));
            AliAODMCParticle* daughter1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetDaughterLabel(1)));
            AliAODMCParticle* mother = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetMother()));
            Int_t pdgCode         = mother->GetPdgCode();
            if(particle->GetPdgCode() == 111){
              Int_t source = GetSourceClassification(111,pdgCode);
              fHistoMCSecPi0PtvsSource[fiCut]->Fill(particle->Pt(),source,fWeightJetJetMC);
              fHistoMCSecPi0Source[fiCut]->Fill(pdgCode);
            } else if(particle->GetPdgCode() == 221){
              fHistoMCSecEtaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
              fHistoMCSecEtaSource[fiCut]->Fill(pdgCode);
            }

            // pi0 really in acceptance/
            if( ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(daughter0,AODMCTrackArray,kFALSE) &&
                ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(daughter1,AODMCTrackArray,kFALSE)  &&
                ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter0->Phi(),fEventPlaneAngle,kFALSE) &&
                ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter1->Phi(),fEventPlaneAngle,kFALSE)){
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
  return;
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::ProcessMCParticles()
{
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();
//   cout << mcProdVtxX <<"\t" << mcProdVtxY << "\t" << mcProdVtxZ << endl;

  // Loop over all primary MC particle
  for(Long_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {

    if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, i, mcProdVtxX, mcProdVtxY, mcProdVtxZ)){
      // fill primary histogram
      TParticle* particle = (TParticle *)fMCEvent->Particle(i);
      if (!particle) continue;

      Int_t isMCFromMBHeader = -1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        isMCFromMBHeader
          = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
        if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
      }

      if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(particle->Phi(),fEventPlaneAngle,kFALSE)) continue;
      if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCEvent,kFALSE)){
        fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // All MC Gamma
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
      if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCEvent,kTRUE)){
        fHistoMCConvGammaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
        if (fDoPhotonQA > 0 && fIsMC < 2){
          fHistoMCConvGammaR[fiCut]->Fill(((TParticle*)fMCEvent->Particle(particle->GetFirstDaughter()))->R(),fWeightJetJetMC);
          fHistoMCConvGammaEta[fiCut]->Fill(particle->Eta(),fWeightJetJetMC);
        }
      } // Converted MC Gamma
      if(fDoMesonAnalysis){

          Double_t mesonY = 1.e30;
          Double_t ratio  = 0;
          if (particle->Energy() != TMath::Abs(particle->Pz())){
            ratio         = (particle->Energy()+particle->Pz()) / (particle->Energy()-particle->Pz());
          }
          if( !(ratio <= 0) ){
            mesonY = particle->Y()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
          }

          if (TMath::Abs(mesonY) < ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetRapidityCutValue()){
            if ( particle->GetPdgCode() == 211 ){  // positve pions
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(),0.,fWeightJetJetMC);
            } else if ( particle->GetPdgCode() == -211 ){  // negative pions
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(),1.,fWeightJetJetMC);
            } else if ( particle->GetPdgCode() == 321 ){  // positve kaons
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(),2.,fWeightJetJetMC);
            } else if ( particle->GetPdgCode() == -321 ){  // negative kaons
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(),3.,fWeightJetJetMC);
            } else if ( TMath::Abs(particle->GetPdgCode()) == 310 ){  // K0s
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(),4.,fWeightJetJetMC);
            } else if ( TMath::Abs(particle->GetPdgCode()) == 130 ){  // K0l
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(),5.,fWeightJetJetMC);
            } else if ( TMath::Abs(particle->GetPdgCode()) == 3122 ){  // Lambda/ AntiLambda
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(),6.,fWeightJetJetMC);
            } else if ( TMath::Abs(particle->GetPdgCode()) == 223 ){  // Omega
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(),7.,fWeightJetJetMC);
            } else if ( TMath::Abs(particle->GetPdgCode()) == 333 ){  // Phi
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(),8.,fWeightJetJetMC);
            } else if ( TMath::Abs(particle->GetPdgCode()) == 113 ){  // Rho0
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(),9.,fWeightJetJetMC);
            }
          }

        if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
          ->MesonIsSelectedMC(particle,fMCEvent,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
          TParticle* daughter0 = (TParticle*)fMCEvent->Particle(particle->GetFirstDaughter());
          TParticle* daughter1 = (TParticle*)fMCEvent->Particle(particle->GetLastDaughter());

          Float_t weighted= 1;
	  if (particle->Pt()>0.005){
	    weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(i, fMCEvent, fInputEvent);
	    //                   if(particle->GetPdgCode() == 221){
	    //                      cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
	    //                   }
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
            if (fIsMC > 1) fHistoMCPi0WOEvtWeightPt[fiCut]->Fill(particle->Pt());
            if (fDoMesonQA > 0){
              if (fIsMC < 2)fHistoMCPi0PtY[fiCut]->Fill(particle->Pt(),mesonY,weighted); // All MC Pi0
              fHistoMCPi0PtAlpha[fiCut]->Fill(particle->Pt(),alpha,fWeightJetJetMC); // All MC Pi0
              if (fIsMC == 2) fHistoMCPi0PtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),fWeightJetJetMC);
            }
          } else if(particle->GetPdgCode() == 221){
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
              if (fIsMC < 2)fHistoMCEtaPtY[fiCut]->Fill(particle->Pt(),mesonY,weighted); // All MC Pi0
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

            if(particle->GetPdgCode() == 111){
              if(fIsMC > 1) fHistoMCPi0WOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Pi0 with gamma in acc NOT weighted at all
              fHistoMCPi0WOWeightInAccPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // MC Pi0 with gamma in acc NOT weighted
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
            } else if(particle->GetPdgCode() == 221){
              if(fIsMC > 1) fHistoMCEtaWOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Eta with gamma in acc NOT weighted at all
              fHistoMCEtaWOWeightInAccPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // MC Eta with gamma in acc NOT weighted
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
            }
          }
        }
      }
    } else {
      // fill secondary histograms
      TParticle* particle = (TParticle *)fMCEvent->Particle(i);
      if (!particle) continue;

      Int_t isMCFromMBHeader = -1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        isMCFromMBHeader
          = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
        if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
      }

      if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(particle->Phi(),fEventPlaneAngle,kFALSE)){
        if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCEvent,kFALSE)){
          if (particle->GetMother(0) > -1 && fMCEvent->Particle(particle->GetMother(0))->GetMother(0) > -1) {
            if (fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 310){
              fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),0.,fWeightJetJetMC);
            } else if (fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 130) {
              fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),1.,fWeightJetJetMC);
            } else if (fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 3122) {
              fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),2.,fWeightJetJetMC);
            } else {
              if ( !(TMath::Abs(fMCEvent->Particle(particle->GetMother(0))->GetPdgCode()) == 11 && fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 22) )
                fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),3.,fWeightJetJetMC);
            }
          } else {
            fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),3.,fWeightJetJetMC);
          }
        }

        if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCEvent,kTRUE)){
          if (particle->GetMother(0) > -1 && fMCEvent->Particle(particle->GetMother(0))->GetMother(0) > -1) {
            if (fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 310){
              fHistoMCSecondaryConvGammaPt[fiCut]->Fill(particle->Pt(),0.,fWeightJetJetMC);
            } else if (fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 130) {
              fHistoMCSecondaryConvGammaPt[fiCut]->Fill(particle->Pt(),1.,fWeightJetJetMC);
            } else if (fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 3122) {
              fHistoMCSecondaryConvGammaPt[fiCut]->Fill(particle->Pt(),2.,fWeightJetJetMC);
            } else {
              if ( !(TMath::Abs(fMCEvent->Particle(particle->GetMother(0))->GetPdgCode()) == 11 && fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 22) )
                fHistoMCSecondaryConvGammaPt[fiCut]->Fill(particle->Pt(),3.,fWeightJetJetMC);
            }
          } else {
            fHistoMCSecondaryConvGammaPt[fiCut]->Fill(particle->Pt(),3.,fWeightJetJetMC);
          }
        }
      }

      if(fDoMesonAnalysis){
        if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedMC(particle,fMCEvent,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
          TParticle* daughter0  = (TParticle*)fMCEvent->Particle(particle->GetFirstDaughter());
          TParticle* daughter1  = (TParticle*)fMCEvent->Particle(particle->GetLastDaughter());
          Int_t pdgCode         = -1;
          if(particle->GetFirstMother()>-1) pdgCode = ((TParticle*)fMCEvent->Particle( particle->GetFirstMother() ))->GetPdgCode();
          if(particle->GetPdgCode() == 111){
            Int_t source = GetSourceClassification(111,pdgCode);
            fHistoMCSecPi0PtvsSource[fiCut]->Fill(particle->Pt(),source,fWeightJetJetMC);

            Double_t deltaX = particle->Vx() - mcProdVtxX;
            Double_t deltaY = particle->Vy() - mcProdVtxY;
            Double_t deltaZ = particle->Vz() - mcProdVtxZ;
            Double_t realRadius3D = TMath::Sqrt(deltaX*deltaX+deltaY*deltaY+deltaZ*deltaZ);
            if(fDoMesonQA > 0 && fIsMC < 2) fHistoMCSecPi0RvsSource[fiCut]->Fill(realRadius3D,source);
            fHistoMCSecPi0Source[fiCut]->Fill(pdgCode);
          } else if(particle->GetPdgCode() == 221){
            fHistoMCSecEtaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
            fHistoMCSecEtaSource[fiCut]->Fill(pdgCode);
          }

          // pi0 really in acceptance/
          if( ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(daughter0,fMCEvent,kFALSE) &&
              ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(daughter1,fMCEvent,kFALSE)  &&
              ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter0->Phi(),fEventPlaneAngle,kFALSE) &&
              ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter1->Phi(),fEventPlaneAngle,kFALSE)){
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
void AliAnalysisTaskGammaConvV1::CalculatePi0Candidates(){
  // Conversion Gammas
  if(fGammaCandidates->GetEntries()>1){
    for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries()-1;firstGammaIndex++){
      AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));
      if (gamma0==NULL) continue;
      for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fGammaCandidates->GetEntries();secondGammaIndex++){
        AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(secondGammaIndex));
        //Check for same Electron ID
        if (gamma1==NULL) continue;
        if(gamma0->GetTrackLabelPositive() == gamma1->GetTrackLabelPositive() ||
        gamma0->GetTrackLabelNegative() == gamma1->GetTrackLabelNegative() ||
        gamma0->GetTrackLabelNegative() == gamma1->GetTrackLabelPositive() ||
        gamma0->GetTrackLabelPositive() == gamma1->GetTrackLabelNegative() ) continue;

        AliAODConversionMother *pi0cand = new AliAODConversionMother(gamma0,gamma1);
        pi0cand->SetLabels(firstGammaIndex,secondGammaIndex);
        pi0cand->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());

        if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
          if(fDoCentralityFlat > 0){
            fHistoMotherInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(), fWeightCentrality[fiCut]*fWeightJetJetMC);
            if(TMath::Abs(pi0cand->GetAlpha())<0.1) fHistoMotherInvMassEalpha[fiCut]->Fill(pi0cand->M(),pi0cand->E(), fWeightCentrality[fiCut]*fWeightJetJetMC);
          } else {
            if(TMath::Abs(pi0cand->GetAlpha())<0.1) fHistoMotherInvMassEalpha[fiCut]->Fill(pi0cand->M(),pi0cand->E(),fWeightJetJetMC);
            if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoMotherInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(),fWeightJetJetMC);
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
          }

          if (fDoMesonQA > 0){

            if(fDoMesonQA == 3 && TMath::Abs(gamma0->GetConversionRadius()-gamma1->GetConversionRadius())<10 && pi0cand->GetOpeningAngle()<0.1){
                    Double_t sparesFill[4] = {gamma0->GetPhotonPt(),gamma0->GetConversionRadius(),TMath::Abs(gamma0->GetConversionRadius()-gamma1->GetConversionRadius()),pi0cand->GetOpeningAngle()};
                    sPtRDeltaROpenAngle[fiCut]->Fill(sparesFill, 1);
            }

            if ( pi0cand->M() > 0.05 && pi0cand->M() < 0.17){
              if (fIsMC < 2){
                fHistoMotherPi0PtY[fiCut]->Fill(pi0cand->Pt(),pi0cand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift());
                fHistoMotherPi0PtOpenAngle[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetOpeningAngle());
              }
              fHistoMotherPi0PtAlpha[fiCut]->Fill(pi0cand->Pt(),TMath::Abs(pi0cand->GetAlpha()),fWeightJetJetMC);

            }
            if ( pi0cand->M() > 0.45 && pi0cand->M() < 0.65){
              if (fIsMC < 2){
                fHistoMotherEtaPtY[fiCut]->Fill(pi0cand->Pt(),pi0cand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift());
                fHistoMotherEtaPtOpenAngle[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetOpeningAngle());
              }
              fHistoMotherEtaPtAlpha[fiCut]->Fill(pi0cand->Pt(),TMath::Abs(pi0cand->GetAlpha()),fWeightJetJetMC);
            }
          }
          if(fDoTHnSparse && ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoBGCalculation()){
            Int_t psibin = 0;
            Int_t zbin = 0;
            Int_t mbin = 0;

            Double_t sparesFill[4];
            if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->BackgroundHandlerType() == 0){
              zbin = fBGHandler[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
              if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
                mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
              } else {
                mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fGammaCandidates->GetEntries());
              }
              sparesFill[0] = pi0cand->M();
              sparesFill[1] = pi0cand->Pt();
              sparesFill[2] = (Double_t)zbin;
              sparesFill[3] = (Double_t)mbin;
            } else {
              psibin = fBGHandlerRP[fiCut]->GetRPBinIndex(TMath::Abs(fEventPlaneAngle));
              zbin = fBGHandlerRP[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
//               if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
//                 mbin = fBGHandlerRP[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
//               } else {
//                 mbin = fBGHandlerRP[fiCut]->GetMultiplicityBinIndex(fGammaCandidates->GetEntries());
//               }
              sparesFill[0] = pi0cand->M();
              sparesFill[1] = pi0cand->Pt();
              sparesFill[2] = (Double_t)zbin;
              sparesFill[3] = (Double_t)psibin;
            }
//             Double_t sparesFill[4] = {pi0cand->M(),pi0cand->Pt(),(Double_t)zbin,(Double_t)mbin};
            if(fDoCentralityFlat > 0) sESDMotherInvMassPtZM[fiCut]->Fill(sparesFill, fWeightCentrality[fiCut]*fWeightJetJetMC); //instead of weight 1
            else  sESDMotherInvMassPtZM[fiCut]->Fill(sparesFill, fWeightJetJetMC);
          }


          if( fIsMC > 0 ){
            if(fInputEvent->IsA()==AliESDEvent::Class())
              ProcessTrueMesonCandidates(pi0cand,gamma0,gamma1);
            if(fInputEvent->IsA()==AliAODEvent::Class())
              ProcessTrueMesonCandidatesAOD(pi0cand,gamma0,gamma1);
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
            iFlag = pi0cand->GetMesonQuality();
    //                   cout << "gamma 0: " << gamma0->GetV0Index()<< "\t" << gamma0->GetPx() << "\t" << gamma0->GetPy() << "\t" <<  gamma0->GetPz() << "\t" << endl;
    //                   cout << "gamma 1: " << gamma1->GetV0Index()<< "\t"<< gamma1->GetPx() << "\t" << gamma1->GetPy() << "\t" <<  gamma1->GetPz() << "\t" << endl;
    //                    cout << "pi0: "<<fInvMass << "\t" << fPt <<"\t" << fDCAzGammaMin << "\t" << fDCAzGammaMax << "\t" << (Int_t)iFlag << "\t" << (Int_t)iMesonMCInfo <<endl;
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
void AliAnalysisTaskGammaConvV1::ProcessTrueMesonCandidates(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1)
{
  // Process True Mesons
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  iMesonMCInfo = 0;
  if(TrueGammaCandidate0->GetV0Index()<fInputEvent->GetNumberOfV0s()){
    Bool_t isTruePi0 = kFALSE;
    Bool_t isTrueEta = kFALSE;
    Bool_t isTruePi0Dalitz = kFALSE;
    Bool_t isTrueEtaDalitz = kFALSE;
    Bool_t gamma0DalitzCand = kFALSE;
    Bool_t gamma1DalitzCand = kFALSE;
    Int_t gamma0MCLabel = TrueGammaCandidate0->GetMCParticleLabel(fMCEvent);
    Int_t gamma0MotherLabel = -1;
    if(gamma0MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
      // Daughters Gamma 0
      TParticle * negativeMC = (TParticle*)TrueGammaCandidate0->GetNegativeMCDaughter(fMCEvent);
      TParticle * positiveMC = (TParticle*)TrueGammaCandidate0->GetPositiveMCDaughter(fMCEvent);
      TParticle * gammaMC0 = (TParticle*)fMCEvent->Particle(gamma0MCLabel);
      if(TMath::Abs(negativeMC->GetPdgCode())==11 && TMath::Abs(positiveMC->GetPdgCode())==11){  // Electrons ...
        if(negativeMC->GetUniqueID() == 5 && positiveMC->GetUniqueID() ==5){ // ... From Conversion ...
          if(gammaMC0->GetPdgCode() == 22){ // ... with Gamma Mother
            gamma0MotherLabel=gammaMC0->GetFirstMother();
          }
        }
        if(gammaMC0->GetPdgCode() ==111){ // Dalitz candidate
          gamma0DalitzCand = kTRUE;
          gamma0MotherLabel=-111;
        }
        if(gammaMC0->GetPdgCode() ==221){ // Dalitz candidate
          gamma0DalitzCand = kTRUE;
          gamma0MotherLabel=-221;
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
          if(gammaMC1->GetPdgCode() ==221){ // Dalitz candidate
            gamma1DalitzCand = kTRUE;
            gamma1MotherLabel=-221;
          }
        }
      }
      if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
        if(((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->GetPdgCode() == 111){
          isTruePi0=kTRUE;
          if (CheckVectorForDoubleCount(vecDoubleCountTruePi0s,gamma0MotherLabel)){
            fHistoDoubleCountTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC);
            FillMultipleCountMap(mapMultipleCountTruePi0s,gamma0MotherLabel);
          }
        }
        if(((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->GetPdgCode() == 221){
          isTrueEta=kTRUE;
          if (CheckVectorForDoubleCount(vecDoubleCountTrueEtas,gamma0MotherLabel)){
            fHistoDoubleCountTrueEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC);
            FillMultipleCountMap(mapMultipleCountTrueEtas,gamma0MotherLabel);
          }
        }
      }

      //Identify Dalitz candidate
      if (gamma1DalitzCand || gamma0DalitzCand){
        if (gamma0DalitzCand && gamma0MCLabel >=0 && gamma0MCLabel==gamma1MotherLabel){
          if (gamma0MotherLabel == -111) isTruePi0Dalitz = kTRUE;
          if (gamma0MotherLabel == -221) isTrueEtaDalitz = kTRUE;
        }
        if (gamma1DalitzCand && gamma1MCLabel >=0 && gamma1MCLabel==gamma0MotherLabel){
          if (gamma1MotherLabel == -111) isTruePi0Dalitz = kTRUE;
          if (gamma1MotherLabel == -221) isTrueEtaDalitz = kTRUE;
        }
      }


      if(isTruePi0 || isTrueEta){// True Pion or Eta

        Float_t weightMatBudget = 1.;
        if (fDoMaterialBudgetWeightingOfGammasForTrueMesons && ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetMaterialBudgetWeightsInitialized()) {
                weightMatBudget = ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetMaterialBudgetCorrectingWeightForTrueGamma(TrueGammaCandidate0) * ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetMaterialBudgetCorrectingWeightForTrueGamma(TrueGammaCandidate1);
        }

        fHistoTrueMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightMatBudget*fWeightJetJetMC);
        if(fDoJetAnalysis && !fDoLightOutput){
          if(fConvJetReader->GetTrueNJets()>0){
            fHistoTrueJetMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),  fWeightJetJetMC*weightMatBudget);
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
                  fHistoTrueInJetMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC*weightMatBudget);
                }
              }
              fHistoTrueDoubleCountingJet[fiCut]->Fill(counter);
            }
            fTrueVectorJetEta.clear();
            fTrueVectorJetPhi.clear();
          }
        }
        if (fDoMesonQA > 0){
          if (isTruePi0){
            if ( Pi0Candidate->M() > 0.05 && Pi0Candidate->M() < 0.17){
              if(fIsMC < 2){
                fHistoTruePi0PtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift());
                fHistoTruePi0PtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle());
              }
              fHistoTruePi0PtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),TMath::Abs(Pi0Candidate->GetAlpha()),fWeightJetJetMC);

            }
          } else if (isTrueEta){
            if ( Pi0Candidate->M() > 0.45 && Pi0Candidate->M() < 0.65){
              if(fIsMC < 2){
                fHistoTrueEtaPtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift());
                fHistoTrueEtaPtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle());
              }
              fHistoTrueEtaPtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),TMath::Abs(Pi0Candidate->GetAlpha()),fWeightJetJetMC);
            }
          }
        }

        Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, gamma0MotherLabel, mcProdVtxX, mcProdVtxY, mcProdVtxZ);

        if(!isPrimary && gamma1MotherLabel>-1){ // Secondary Meson
          Long_t secMotherLabel = ((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->GetMother(0);
          Float_t weightedSec= 1;
          if(fMCEvent->Particle(secMotherLabel)->GetPdgCode()==310){
            weightedSec= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(secMotherLabel, fMCEvent, fInputEvent)/2.; //invariant mass is additive thus the weight for the daughters has to be devide by two for the K0s at a certain pt
            //cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
          }
          if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTrueSecondaryMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
            if(fDoJetAnalysis){
              if(fConvJetReader->GetTrueNJets()>0){
                if(!fDoLightOutput) fHistoTrueSecondaryInvJetMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
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
                        if(!fDoLightOutput) fHistoTrueSecondaryInvinJetMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
                        else fHistoTrueSecondaryMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
                      }
                    }
                  }
                fTrueVectorJetEta.clear();
                fTrueVectorJetPhi.clear();
              }
            }
          iMesonMCInfo = 2;
          if (secMotherLabel >-1){
            if(fMCEvent->Particle(secMotherLabel)->GetPdgCode()==310){
              iMesonMCInfo = 4;
              if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTrueSecondaryMotherFromK0sInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
              if(fDoJetAnalysis){
                if(fConvJetReader->GetTrueNJets()>0){
                  if(!fDoLightOutput) fHistoTrueSecondaryFromK0sJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
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
                          if(!fDoLightOutput) fHistoTrueSecondaryFromK0sinJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
                          else fHistoTrueSecondaryMotherFromK0sInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
                        }
                      }
                    }
                  fTrueVectorJetEta.clear();
                  fTrueVectorJetPhi.clear();
                }
              }
              if (fDoMesonQA > 0 && fIsMC < 2 )fHistoTrueK0sWithPi0DaughterMCPt[fiCut]->Fill(fMCEvent->Particle(secMotherLabel)->Pt());
            }
            if(fMCEvent->Particle(secMotherLabel)->GetPdgCode()==130){
              iMesonMCInfo = 8;
              if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTrueSecondaryMotherFromK0lInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
              if(fDoJetAnalysis){
                if(fConvJetReader->GetTrueNJets()>0){
                if(!fDoLightOutput) fHistoTrueSecondaryFromK0lJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
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
                        if(!fDoLightOutput) fHistoTrueSecondaryFromK0linJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
                        else fHistoTrueSecondaryMotherFromK0lInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
                      }
                    }
                  }
                  fTrueVectorJetEta.clear();
                  fTrueVectorJetPhi.clear();
                }
              }
              if (fDoMesonQA > 0 && fIsMC < 2 )fHistoTrueK0lWithPi0DaughterMCPt[fiCut]->Fill(fMCEvent->Particle(secMotherLabel)->Pt());
            }
            if(fMCEvent->Particle(secMotherLabel)->GetPdgCode()==221){
              iMesonMCInfo = 3;
              fHistoTrueSecondaryMotherFromEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*weightMatBudget*fWeightJetJetMC);
              if (fDoMesonQA > 0 && fIsMC < 2)fHistoTrueEtaWithPi0DaughterMCPt[fiCut]->Fill(fMCEvent->Particle(secMotherLabel)->Pt());
            }
            if(fMCEvent->Particle(secMotherLabel)->GetPdgCode()==3122){
              iMesonMCInfo = 7;
              if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTrueSecondaryMotherFromLambdaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
              if(fDoJetAnalysis){
                if(fConvJetReader->GetTrueNJets()>0){
                  if(!fDoLightOutput) fHistoTrueSecondaryFromLambdaJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
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
                        if(!fDoLightOutput) fHistoTrueSecondaryFromLambdainJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
                        else fHistoTrueSecondaryMotherFromLambdaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
                      }
                    }
                  }
                  fTrueVectorJetEta.clear();
                  fTrueVectorJetPhi.clear();
                }
              }
              if (fDoMesonQA > 0 && fIsMC < 2)fHistoTrueLambdaWithPi0DaughterMCPt[fiCut]->Fill(fMCEvent->Particle(secMotherLabel)->Pt());
            }
          }
        } else { // Only primary pi0 for efficiency calculation
          iMesonMCInfo = 6;
          Float_t weighted= 1;
	  if (((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt()>0.005){
	    weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(gamma1MotherLabel, fMCEvent, fInputEvent);
	    //                      cout << "rec \t " <<gamma1MotherLabel << "\t" <<  weighted << endl;
	  }

          fHistoTruePrimaryMotherW0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC);
          pESDTruePrimaryMotherWeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*weightMatBudget*fWeightJetJetMC);
          if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTruePrimaryMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudget);
            if(fDoJetAnalysis){
              if(fConvJetReader->GetTrueNJets()>0){
                if(!fDoLightOutput) fHistoTruePrimaryJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudget);
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
                        fHistoTruePrimaryinJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudget);
                        Double_t dotproduct = fTrueVectorJetPx.at(i)*Pi0Candidate->Px() + fTrueVectorJetPy.at(i)*Pi0Candidate->Py() + fTrueVectorJetPz.at(i)*Pi0Candidate->Pz();
                        Double_t magn = pow(fTrueVectorJetPx.at(i),2) + pow(fTrueVectorJetPy.at(i),2) + pow(fTrueVectorJetPz.at(i),2);
                        Double_t z = dotproduct/magn;
                        fHistoTrueJetFragmFunc[fiCut]->Fill(z,fTrueVectorJetPt.at(i), weighted*fWeightJetJetMC*weightMatBudget);
                        fHistoTrueJetFragmFuncZInvMass[fiCut]->Fill(Pi0Candidate->M(), z, weighted*fWeightJetJetMC*weightMatBudget);
                      }
                      else fHistoTruePrimaryMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudget);
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


          if (fDoMesonQA > 0 && fIsMC < 2){
            if(isTruePi0){ // Only primary pi0 for resolution
              fHistoTruePrimaryPi0MCPtResolPt[fiCut]->Fill(((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt(),(Pi0Candidate->Pt()-((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt())/((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt(),weighted*weightMatBudget);
            }
            if (isTrueEta){ // Only primary eta for resolution
              fHistoTruePrimaryEtaMCPtResolPt[fiCut]->Fill(((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt(),(Pi0Candidate->Pt()-((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt())/((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt(),weighted*weightMatBudget);
            }
          }
        }
      } else if(!isTruePi0 && !isTrueEta){ // Background
        if (fDoMesonQA > 0 && fIsMC < 2){
          if(gamma0MotherLabel>-1 && gamma1MotherLabel>-1){ // Both Tracks are Photons and have a mother but not Pi0 or Eta
            fHistoTrueBckGGInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
            iMesonMCInfo = 1;
          } else { // No photon or without mother
            fHistoTrueBckContInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
          }
        }
        if( isTruePi0Dalitz || isTrueEtaDalitz ){
        // Dalitz
          iMesonMCInfo = 5;
          if (fIsMC < 2) fHistoTrueMotherDalitzInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        } else if (gamma0DalitzCand || gamma1DalitzCand){
          if (fDoMesonQA > 0 && fIsMC < 2) fHistoTrueBckContInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        }
      }
    }
  }
}

//______________________________________________________________________
void AliAnalysisTaskGammaConvV1::ProcessTrueMesonCandidatesAOD(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1)
{
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  // Process True Mesons
  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  Bool_t isTruePi0 = kFALSE;
  Bool_t isTrueEta = kFALSE;
  Bool_t isTruePi0Dalitz = kFALSE;
  Bool_t isTrueEtaDalitz = kFALSE;
  Bool_t gamma0DalitzCand = kFALSE;
  Bool_t gamma1DalitzCand = kFALSE;

  if (AODMCTrackArray!=NULL && TrueGammaCandidate0 != NULL){
    AliAODMCParticle *positiveMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelPositive()));
    AliAODMCParticle *negativeMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelNegative()));

    iMesonMCInfo = 0;
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
          }
        }
        if(gammaMC0->GetPdgCode() ==111){ // Dalitz candidate
          gamma0DalitzCand = kTRUE;
          gamma0MotherLabel=-111;
        }
        if(gammaMC0->GetPdgCode() ==221){ // Dalitz candidate
          gamma0DalitzCand = kTRUE;
          gamma0MotherLabel=-221;
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
        if(gammaMC1->GetPdgCode() ==221){ // Dalitz candidate
          gamma1DalitzCand = kTRUE;
          gamma1MotherLabel=-221;
        }
      }
    }
    if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
      if(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 111){
        isTruePi0=kTRUE;
        if (CheckVectorForDoubleCount(vecDoubleCountTruePi0s,gamma0MotherLabel)){
          fHistoDoubleCountTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC);
          FillMultipleCountMap(mapMultipleCountTruePi0s,gamma0MotherLabel);
        }
      }
      if(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 221){
        isTrueEta=kTRUE;
        if (CheckVectorForDoubleCount(vecDoubleCountTrueEtas,gamma0MotherLabel)){
          fHistoDoubleCountTrueEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC);
          FillMultipleCountMap(mapMultipleCountTrueEtas,gamma0MotherLabel);
        }
      }
    }

    //Identify Dalitz candidate
    if (gamma1DalitzCand || gamma0DalitzCand){
      if (gamma0DalitzCand && gamma0MCLabel >=0 && gamma0MCLabel==gamma1MotherLabel){
        if (gamma0MotherLabel == -111) isTruePi0Dalitz = kTRUE;
        if (gamma0MotherLabel == -221) isTrueEtaDalitz = kTRUE;
      }
      if (gamma1DalitzCand && gamma1MCLabel >=0 && gamma1MCLabel==gamma0MotherLabel){
        if (gamma1MotherLabel == -111) isTruePi0Dalitz = kTRUE;
        if (gamma1MotherLabel == -221) isTrueEtaDalitz = kTRUE;
      }
    }

    if(isTruePi0 || isTrueEta){// True Pion or Eta

      Float_t weightMatBudget = 1.;
      if (fDoMaterialBudgetWeightingOfGammasForTrueMesons && ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetMaterialBudgetWeightsInitialized()) {
                weightMatBudget = ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetMaterialBudgetCorrectingWeightForTrueGamma(TrueGammaCandidate0) * ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetMaterialBudgetCorrectingWeightForTrueGamma(TrueGammaCandidate1);
      }

      fHistoTrueMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightMatBudget*fWeightJetJetMC);
      if(fDoJetAnalysis && !fDoLightOutput){
        if(fConvJetReader->GetTrueNJets()>0){
          fHistoTrueJetMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),  fWeightJetJetMC*weightMatBudget);
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
                fHistoTrueInJetMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC*weightMatBudget);
              }
            }
            fHistoTrueDoubleCountingJet[fiCut]->Fill(counter);
          }
          fTrueVectorJetEta.clear();
          fTrueVectorJetPhi.clear();
        }
      }
      if (fDoMesonQA > 0){
        if (isTruePi0){
          if ( Pi0Candidate->M() > 0.05 && Pi0Candidate->M() < 0.17){
            if(fIsMC < 2){
              fHistoTruePi0PtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift());
              fHistoTruePi0PtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle());
            }
            fHistoTruePi0PtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),TMath::Abs(Pi0Candidate->GetAlpha()),fWeightJetJetMC);
          }
        } else if (isTrueEta){
          if ( Pi0Candidate->M() > 0.45 && Pi0Candidate->M() < 0.65){
            if(fIsMC < 2){
              fHistoTrueEtaPtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift());
              fHistoTrueEtaPtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle());
            }
            fHistoTrueEtaPtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),TMath::Abs(Pi0Candidate->GetAlpha()),fWeightJetJetMC);
          }
        }
      }
      Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MotherLabel)), mcProdVtxX, mcProdVtxY, mcProdVtxZ);

      if(!isPrimary){ // Secondary Meson
        Long_t secMotherLabel = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->GetMother();
        Float_t weightedSec= 1;
        if(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->GetPdgCode()==310){
          weightedSec= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(secMotherLabel, 0x0, fInputEvent)/2.; //invariant mass is additive thus the weight for the daughters has to be devide by two for the K0s at a certain pt
          //cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
        }
        if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTrueSecondaryMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
        if(fDoJetAnalysis){
          if(fConvJetReader->GetTrueNJets()>0){
            if(!fDoLightOutput) fHistoTrueSecondaryInvJetMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
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
                  if(!fDoLightOutput) fHistoTrueSecondaryInvinJetMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
                  else fHistoTrueSecondaryMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
                }
              }
            }
            fTrueVectorJetEta.clear();
            fTrueVectorJetPhi.clear();
          }
        }
        iMesonMCInfo = 2;
        if (secMotherLabel >-1){
          if(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->GetPdgCode()==310){
            iMesonMCInfo = 4;
            if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTrueSecondaryMotherFromK0sInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
            if(fDoJetAnalysis){
              if(fConvJetReader->GetTrueNJets()>0){
                if(!fDoLightOutput) fHistoTrueSecondaryFromK0sJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
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
                      if(!fDoLightOutput) fHistoTrueSecondaryFromK0sinJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
                      else fHistoTrueSecondaryMotherFromK0sInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
                    }
                  }
                }
                fTrueVectorJetEta.clear();
                fTrueVectorJetPhi.clear();
              }
            }
            if (fDoMesonQA > 0 && fIsMC < 2)fHistoTrueK0sWithPi0DaughterMCPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->Pt());
          }
          if(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->GetPdgCode()==130){
            iMesonMCInfo = 8;
            if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTrueSecondaryMotherFromK0lInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
            if(fDoJetAnalysis){
              if(fConvJetReader->GetTrueNJets()>0){
                if(!fDoLightOutput) fHistoTrueSecondaryFromK0lJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
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
                      if(!fDoLightOutput) fHistoTrueSecondaryFromK0linJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
                      else fHistoTrueSecondaryMotherFromK0lInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
                    }
                  }
                }
                fTrueVectorJetEta.clear();
                fTrueVectorJetPhi.clear();
              }
            }
            if (fDoMesonQA > 0 && fIsMC < 2)fHistoTrueK0lWithPi0DaughterMCPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->Pt());
          }
          if(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->GetPdgCode()==221){
            iMesonMCInfo = 3;
            fHistoTrueSecondaryMotherFromEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*weightMatBudget*fWeightJetJetMC);
            if (fDoMesonQA > 0 && fIsMC < 2)fHistoTrueEtaWithPi0DaughterMCPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->Pt());
          }
          if(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->GetPdgCode()==3122){
            iMesonMCInfo = 7;
            if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTrueSecondaryMotherFromLambdaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
              if(fDoJetAnalysis){
                if(fConvJetReader->GetTrueNJets()>0){
                  if(!fDoLightOutput) fHistoTrueSecondaryFromLambdaJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
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
                        if(!fDoLightOutput) fHistoTrueSecondaryFromLambdainJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
                        else fHistoTrueSecondaryMotherFromLambdaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
                      }
                    }
                  }
                  fTrueVectorJetEta.clear();
                  fTrueVectorJetPhi.clear();
                }
              }
            if (fDoMesonQA > 0 && fIsMC < 2)fHistoTrueLambdaWithPi0DaughterMCPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->Pt());
          }
        }
      } else { // Only primary pi0 for efficiency calculation
        Float_t weighted= 1;
        iMesonMCInfo = 6;
	if (static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt()>0.005){
	  weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(gamma1MotherLabel, 0x0, fInputEvent);
          //                      cout << "rec \t " <<gamma1MotherLabel << "\t" <<  weighted << endl;
	}
        if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoTruePrimaryMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudget);
        if(fDoJetAnalysis){
          if(fConvJetReader->GetTrueNJets()>0){
            if(!fDoLightOutput) fHistoTruePrimaryJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudget);
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
                    fHistoTruePrimaryinJetInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudget);
                    Double_t dotproduct = fTrueVectorJetPx.at(i)*Pi0Candidate->Px() + fTrueVectorJetPy.at(i)*Pi0Candidate->Py() + fTrueVectorJetPz.at(i)*Pi0Candidate->Pz();
                    Double_t magn = pow(fTrueVectorJetPx.at(i),2) + pow(fTrueVectorJetPy.at(i),2) + pow(fTrueVectorJetPz.at(i),2);
                    Double_t z = dotproduct/magn;
                    fHistoTrueJetFragmFunc[fiCut]->Fill(z,fTrueVectorJetPt.at(i), weighted*fWeightJetJetMC*weightMatBudget);
                    fHistoTrueJetFragmFuncZInvMass[fiCut]->Fill(Pi0Candidate->M(), z, weighted*fWeightJetJetMC*weightMatBudget);
                  }
                  else fHistoTruePrimaryMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudget);
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
        fHistoTruePrimaryMotherW0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC);
        pESDTruePrimaryMotherWeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*weightMatBudget*fWeightJetJetMC);

        if (fDoMesonQA > 0 && fIsMC < 2){
          if(isTruePi0){ // Only primary pi0 for resolution
            fHistoTruePrimaryPi0MCPtResolPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt(),
                                (Pi0Candidate->Pt()-static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt())/static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt(),weighted*weightMatBudget);

          }
          if (isTrueEta){ // Only primary eta for resolution
            fHistoTruePrimaryEtaMCPtResolPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt(),
                                (Pi0Candidate->Pt()-static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt())/static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt(),weighted*weightMatBudget);
          }
        }
      }
    } else if(!isTruePi0 && !isTrueEta) { // Background
      if (fDoMesonQA > 0 && fIsMC < 2){
        if(gamma0MotherLabel>-1 && gamma1MotherLabel>-1){ // Both Tracks are Photons and have a mother but not Pi0 or Eta
          fHistoTrueBckGGInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
          iMesonMCInfo = 1;
        } else { // No photon or without mother
          fHistoTrueBckContInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        }
      }
      if( isTruePi0Dalitz || isTrueEtaDalitz ){
        // Dalitz
        iMesonMCInfo = 5;
        if (fIsMC < 2)fHistoTrueMotherDalitzInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      } else if (gamma0DalitzCand || gamma1DalitzCand){
        if (fDoMesonQA > 0 && fIsMC < 2)fHistoTrueBckContInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
    }
  }
  return;
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::CalculateBackground(){
  Int_t zbin = fBGHandler[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
  Int_t mbin = 0;

    if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
        mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
    } else {
        mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fGammaCandidates->GetEntries());
    }

  if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseRotationMethod()){

    for(Int_t iCurrent=0;iCurrent<fGammaCandidates->GetEntries();iCurrent++){
      AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fGammaCandidates->At(iCurrent));
      for(Int_t iCurrent2=iCurrent+1;iCurrent2<fGammaCandidates->GetEntries();iCurrent2++){
        for(Int_t nRandom=0;nRandom<((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetNumberOfBGEvents();nRandom++){
        AliAODConversionPhoton currentEventGoodV02 = *(AliAODConversionPhoton*)(fGammaCandidates->At(iCurrent2));

        if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoBGProbability()){
          AliAODConversionMother *backgroundCandidateProb = new AliAODConversionMother(&currentEventGoodV0,&currentEventGoodV02);
          Double_t massBGprob = backgroundCandidateProb->M();
          if(massBGprob>0.1 && massBGprob<0.14){
            if(fRandom.Rndm()>fBGHandler[fiCut]->GetBGProb(zbin,mbin)){
              delete backgroundCandidateProb;
              continue;
            }
          }
          delete backgroundCandidateProb;
          backgroundCandidateProb = 0x0;
        }

        RotateParticle(&currentEventGoodV02);
        AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&currentEventGoodV02);
        backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
        if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
          ->MesonIsSelected(backgroundCandidate,kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
          if(fDoCentralityFlat > 0) fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(), fWeightCentrality[fiCut]*fWeightJetJetMC);
          else fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(),fWeightJetJetMC);
          if(fDoTHnSparse){
            Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
            if(fDoCentralityFlat > 0) sESDMotherBackInvMassPtZM[fiCut]->Fill(sparesFill, fWeightCentrality[fiCut]*fWeightJetJetMC); //instead of weight 1
            else sESDMotherBackInvMassPtZM[fiCut]->Fill(sparesFill, fWeightJetJetMC);
          }
        }
        delete backgroundCandidate;
        backgroundCandidate = 0x0;
        }
      }
    }
  } else {
    AliGammaConversionAODBGHandler::GammaConversionVertex *bgEventVertex = NULL;

    if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
      for(Int_t nEventsInBG=0;nEventsInBG<fBGHandler[fiCut]->GetNBGEvents();nEventsInBG++){
        AliGammaConversionAODVector *previousEventV0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,nEventsInBG);
        if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
          bgEventVertex = fBGHandler[fiCut]->GetBGEventVertex(zbin,mbin,nEventsInBG);
        }

        for(Int_t iCurrent=0;iCurrent<fGammaCandidates->GetEntries();iCurrent++){
        AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fGammaCandidates->At(iCurrent));
        for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
          AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));
          if(fMoveParticleAccordingToVertex == kTRUE){
            MoveParticleAccordingToVertex(&previousGoodV0,bgEventVertex);
          }
          if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
            RotateParticleAccordingToEP(&previousGoodV0,bgEventVertex->fEP,fEventPlaneAngle);
          }

          AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
          backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
          if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
            ->MesonIsSelected(backgroundCandidate,kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
            if(fDoCentralityFlat > 0) fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(), fWeightCentrality[fiCut]*fWeightJetJetMC);
            else fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(),fWeightJetJetMC);
            if(fDoTHnSparse){
              Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
              if(fDoCentralityFlat > 0) sESDMotherBackInvMassPtZM[fiCut]->Fill(sparesFill, fWeightCentrality[fiCut]*fWeightJetJetMC); //instead of weight 1
              else sESDMotherBackInvMassPtZM[fiCut]->Fill(sparesFill, fWeightJetJetMC);
            }
          }
          delete backgroundCandidate;
          backgroundCandidate = 0x0;
        }
        }
      }
    } else {
      for(Int_t nEventsInBG=0;nEventsInBG <fBGHandler[fiCut]->GetNBGEvents();nEventsInBG++){
        AliGammaConversionAODVector *previousEventV0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,nEventsInBG);
        if(previousEventV0s){
        if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
          bgEventVertex = fBGHandler[fiCut]->GetBGEventVertex(zbin,mbin,nEventsInBG);
        }
        for(Int_t iCurrent=0;iCurrent<fGammaCandidates->GetEntries();iCurrent++){
          AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fGammaCandidates->At(iCurrent));
          for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){

            AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));

            if(fMoveParticleAccordingToVertex == kTRUE){
              MoveParticleAccordingToVertex(&previousGoodV0,bgEventVertex);
            }
            if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
              RotateParticleAccordingToEP(&previousGoodV0,bgEventVertex->fEP,fEventPlaneAngle);
            }


            AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
            backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
            if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
              ->MesonIsSelected(backgroundCandidate,kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
              if(fDoCentralityFlat > 0) fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(), fWeightCentrality[fiCut]*fWeightJetJetMC);
              else{
                  if(!fDoJetAnalysis || (fDoJetAnalysis && !fDoLightOutput)) fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(), fWeightJetJetMC);
                  if(fDoJetAnalysis){
                  if(fConvJetReader->GetNJets() > 0){
                    if(!fDoLightOutput) fHistoMotherBackJetInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(), fWeightJetJetMC);
                    else fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(), fWeightJetJetMC);
                  }
                }
              }
              if(fDoTHnSparse){
                Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
                if(fDoCentralityFlat > 0) sESDMotherBackInvMassPtZM[fiCut]->Fill(sparesFill, fWeightCentrality[fiCut]*fWeightJetJetMC); //instead of weight 1
                else sESDMotherBackInvMassPtZM[fiCut]->Fill(sparesFill, fWeightJetJetMC);
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
void AliAnalysisTaskGammaConvV1::CalculateBackgroundRP(){

  Int_t psibin = 0;
  Int_t zbin = 0;

  if(fDoTHnSparse){
    psibin = fBGHandlerRP[fiCut]->GetRPBinIndex(TMath::Abs(fEventPlaneAngle));
    zbin = fBGHandlerRP[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
//     if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
//       mbin = fBGHandlerRP[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
//     } else {
//       mbin = fBGHandlerRP[fiCut]->GetMultiplicityBinIndex(fGammaCandidates->GetEntries());
//     }
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
          if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(&backgroundCandidate,kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
            if(fDoCentralityFlat > 0) fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate.M(),backgroundCandidate.Pt(), fWeightCentrality[fiCut]*fWeightJetJetMC);
            else fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate.M(),backgroundCandidate.Pt(),fWeightJetJetMC);
            if(fDoTHnSparse){
//               Double_t sparesFill[4] = {backgroundCandidate.M(),backgroundCandidate.Pt(),(Double_t)zbin,(Double_t)mbin};
              Double_t sparesFill[4] = {backgroundCandidate.M(),backgroundCandidate.Pt(),(Double_t)zbin,(Double_t)psibin};
              if(fDoCentralityFlat > 0) sESDMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,weight*fWeightCentrality[fiCut]*fWeightJetJetMC);
              else sESDMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,weight*fWeightJetJetMC);
            }
          }
        }
      }
    }

  } else {
    // Do Event Mixing
    for(Int_t nEventsInBG=0;nEventsInBG <fBGHandlerRP[fiCut]->GetNBGEvents(fGammaCandidates,fInputEvent);nEventsInBG++){

      AliGammaConversionPhotonVector *previousEventGammas = fBGHandlerRP[fiCut]->GetBGGoodGammas(fGammaCandidates,fInputEvent,nEventsInBG);

      if(previousEventGammas){
        // test weighted background
        Double_t weight=1.0;
        // Correct for the number of eventmixing:
        // N gammas -> (N-1) + (N-2) +(N-3) ...+ (N-(N-1))  using sum formula sum(i)=N*(N-1)/2  -> N*(N-1)/2
        // real combinations (since you cannot combine a photon with its own)
        // but BG leads to N_{a}*N_{b} combinations
        weight*=0.5*(Double_t(fGammaCandidates->GetEntries()-1))/Double_t(previousEventGammas->size());

        for(Int_t iCurrent=0;iCurrent<fGammaCandidates->GetEntries();iCurrent++){

                    AliAODConversionPhoton *gamma0 = (AliAODConversionPhoton*)(fGammaCandidates->At(iCurrent));

                    for(UInt_t iPrevious=0;iPrevious<previousEventGammas->size();iPrevious++){

                        AliAODConversionPhoton *gamma1 = (AliAODConversionPhoton*)(previousEventGammas->at(iPrevious));

                        AliAODConversionMother backgroundCandidate(gamma0,gamma1);
                        backgroundCandidate.CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
                        if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
                            ->MesonIsSelected(&backgroundCandidate,kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
                            if(fDoCentralityFlat > 0) fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate.M(),backgroundCandidate.Pt(), fWeightCentrality[fiCut]*fWeightJetJetMC);
              else fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate.M(),backgroundCandidate.Pt(),fWeightJetJetMC);
                            if(fDoTHnSparse){
//                              Double_t sparesFill[4] = {backgroundCandidate.M(),backgroundCandidate.Pt(),(Double_t)zbin,(Double_t)mbin};
                                Double_t sparesFill[4] = {backgroundCandidate.M(),backgroundCandidate.Pt(),(Double_t)zbin,(Double_t)psibin};
                                if(fDoCentralityFlat > 0) sESDMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,weight*fWeightCentrality[fiCut]*fWeightJetJetMC);
                                else sESDMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,weight*fWeightJetJetMC);
                            }
                        }
                    }
        }
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::RotateParticle(AliAODConversionPhoton *gamma){
  Int_t fNDegreesPMBackground= ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->NDegreesRotation();
  Double_t nRadiansPM = fNDegreesPMBackground*TMath::Pi()/180;
  Double_t rotationValue = fRandom.Rndm()*2*nRadiansPM + TMath::Pi()-nRadiansPM;
  gamma->RotateZ(rotationValue);
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::RotateParticleAccordingToEP(AliAODConversionPhoton *gamma, Double_t previousEventEP, Double_t thisEventEP){

  previousEventEP=previousEventEP+TMath::Pi();
  thisEventEP=thisEventEP+TMath::Pi();
  Double_t rotationValue= thisEventEP-previousEventEP;
  gamma->RotateZ(rotationValue);
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::MoveParticleAccordingToVertex(AliAODConversionPhoton* particle,const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex){
  //see header file for documentation

  Double_t dx = vertex->fX - fInputEvent->GetPrimaryVertex()->GetX();
  Double_t dy = vertex->fY - fInputEvent->GetPrimaryVertex()->GetY();
  Double_t dz = vertex->fZ - fInputEvent->GetPrimaryVertex()->GetZ();

  Double_t movedPlace[3] = {particle->GetConversionX() - dx,particle->GetConversionY() - dy,particle->GetConversionZ() - dz};
  particle->SetConversionPoint(movedPlace);
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::UpdateEventByEventData(){
  //see header file for documentation
  if(fDoJetAnalysis && fConvJetReader->GetNJets() == 0) return;
  if(fGammaCandidates->GetEntries() >0 ){
    if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
      fBGHandler[fiCut]->AddEvent(fGammaCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fV0Reader->GetNumberOfPrimaryTracks(),fEventPlaneAngle);
    }
    else{ // means we use #V0s for multiplicity
      fBGHandler[fiCut]->AddEvent(fGammaCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fGammaCandidates->GetEntries(),fEventPlaneAngle);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::FillPhotonCombinatorialBackgroundHist(AliAODConversionPhoton *TruePhotonCandidate, Int_t pdgCode[], Double_t PhiParticle[])
{
  // Combinatorial Bck = 0 ee, 1 epi, 2 ek, 3 ep, 4 emu, 5 pipi, 6 pik, 7 pip, 8 pimu, 9 kk, 10 kp, 11 kmu, 12 pp, 13 pmu, 14 mumu, 15 Rest
  if(pdgCode[0]==11   && pdgCode[1]==11){
    fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),0.,fWeightJetJetMC);
  } else if( (pdgCode[0]==11   && pdgCode[1]==211) || (pdgCode[0]==211  && pdgCode[1]==11) ){
    fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),1.,fWeightJetJetMC);
  } else if( (pdgCode[0]==11   && pdgCode[1]==321) || (pdgCode[0]==321  && pdgCode[1]==11) ){
    fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),2.,fWeightJetJetMC);
  } else if( (pdgCode[0]==11   && pdgCode[1]==2212) || (pdgCode[0]==2212 && pdgCode[1]==11) ){
    fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3.,fWeightJetJetMC);
  } else if( (pdgCode[0]==11   && pdgCode[1]==13) || (pdgCode[0]==13   && pdgCode[1]==11) ){
    fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),4.,fWeightJetJetMC);
  } else if(  pdgCode[0]==211  && pdgCode[1]==211 ){
    fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),5.,fWeightJetJetMC);
  } else if( (pdgCode[0]==211  && pdgCode[1]==321) || (pdgCode[0]==321  && pdgCode[1]==211) ){
    fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),6.,fWeightJetJetMC);
  } else if( (pdgCode[0]==211  && pdgCode[1]==2212) || (pdgCode[0]==2212 && pdgCode[1]==211) ){
    fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),7.,fWeightJetJetMC);
  } else if( (pdgCode[0]==211  && pdgCode[1]==13) || (pdgCode[0]==13   && pdgCode[1]==211) ){
    fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),8.,fWeightJetJetMC);
  } else if(  pdgCode[0]==321  && pdgCode[1]==321 ){
    fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),9.,fWeightJetJetMC);
  } else if( (pdgCode[0]==321  && pdgCode[1]==2212) || (pdgCode[0]==2212 && pdgCode[1]==321) ){
    fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),10.,fWeightJetJetMC);
  } else if( (pdgCode[0]==321  && pdgCode[1]==13) || (pdgCode[0]==13   && pdgCode[1]==321) ){
    fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),11.,fWeightJetJetMC);
  } else if(  pdgCode[0]==2212   && pdgCode[1]==2212  ){
    fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),12.,fWeightJetJetMC);
  } else if( (pdgCode[0]==2212  && pdgCode[1]==13) || (pdgCode[0]==13   && pdgCode[1]==2212) ){
    fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),13.,fWeightJetJetMC);
  } else if(  pdgCode[0]==13   && pdgCode[1]==13  ){
    fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),14.,fWeightJetJetMC);
  } else {
    fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),15.,fWeightJetJetMC);
  }

  if(fDoPhotonQA == 3 && fIsMC < 2){
    if( (pdgCode[0]==11   && pdgCode[1]==211) || (pdgCode[0]==211  && pdgCode[1]==11) ){
      if(pdgCode[0]==11){fHistoCombinatorialPtDeltaPhi_epi[fiCut]->Fill(TruePhotonCandidate->Pt(),TMath::Abs((TruePhotonCandidate->Phi() - PhiParticle[1])));}
      else if(pdgCode[0]==211){fHistoCombinatorialPtDeltaPhi_epi[fiCut]->Fill(TruePhotonCandidate->Pt(),TMath::Abs((TruePhotonCandidate->Phi() - PhiParticle[0])));}
    }
    else if( (pdgCode[0]==11   && pdgCode[1]==321) || (pdgCode[0]==321  && pdgCode[1]==11) ){
      if(pdgCode[0]==11){fHistoCombinatorialPtDeltaPhi_ek[fiCut]->Fill(TruePhotonCandidate->Pt(),TMath::Abs((TruePhotonCandidate->Phi() - PhiParticle[1])));}
      else if(pdgCode[0]==321){fHistoCombinatorialPtDeltaPhi_ek[fiCut]->Fill(TruePhotonCandidate->Pt(),TMath::Abs((TruePhotonCandidate->Phi() - PhiParticle[0])));}
    }
    else if( (pdgCode[0]==11   && pdgCode[1]==2212) || (pdgCode[0]==2212 && pdgCode[1]==11) ){
      if(pdgCode[0]==11){fHistoCombinatorialPtDeltaPhi_ep[fiCut]->Fill(TruePhotonCandidate->Pt(),TMath::Abs((TruePhotonCandidate->Phi() - PhiParticle[1])));}
      else if(pdgCode[0]==2212){fHistoCombinatorialPtDeltaPhi_ep[fiCut]->Fill(TruePhotonCandidate->Pt(),TMath::Abs((TruePhotonCandidate->Phi() - PhiParticle[0])));}
    }
    else if( (pdgCode[0]==211  && pdgCode[1]==321) || (pdgCode[0]==321  && pdgCode[1]==211) ){
      if(pdgCode[0]==211){fHistoCombinatorialPtDeltaPhi_pik[fiCut]->Fill(TruePhotonCandidate->Pt(),TMath::Abs((TruePhotonCandidate->Phi() - PhiParticle[1])));}
      else if(pdgCode[0]==321){fHistoCombinatorialPtDeltaPhi_pik[fiCut]->Fill(TruePhotonCandidate->Pt(),TMath::Abs((TruePhotonCandidate->Phi() - PhiParticle[0])));}
    }
    else if( (pdgCode[0]==211  && pdgCode[1]==2212) || (pdgCode[0]==2212 && pdgCode[1]==211) ){
      if(pdgCode[0]==211){fHistoCombinatorialPtDeltaPhi_pip[fiCut]->Fill(TruePhotonCandidate->Pt(),TMath::Abs((TruePhotonCandidate->Phi() - PhiParticle[1])));}
      else if(pdgCode[0]==2212){fHistoCombinatorialPtDeltaPhi_pip[fiCut]->Fill(TruePhotonCandidate->Pt(),TMath::Abs((TruePhotonCandidate->Phi() - PhiParticle[0])));}
    }

  }
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::FillPhotonCombinatorialMothersHistESD(TParticle *daughter, TParticle *motherCombPart)
{
  Int_t pdgCombPart = TMath::Abs(daughter->GetPdgCode());
  Int_t pdgMotherCombPart = TMath::Abs(motherCombPart->GetPdgCode());

  // Combinatorial Bck mothers: e, pi, k, p, pi0, eta, omega, phi, eta', K0s, Lambda, rhos, other, not prim
  if(pdgCombPart==11){ //comb particle is an electron
    if(pdgMotherCombPart==11)         fHistoCombinatorialMothersPt[fiCut]->Fill(0.,0.,motherCombPart->Pt());
    else if(pdgMotherCombPart==211)   fHistoCombinatorialMothersPt[fiCut]->Fill(0.,1.,motherCombPart->Pt());
    else if(pdgMotherCombPart==321)   fHistoCombinatorialMothersPt[fiCut]->Fill(0.,2.,motherCombPart->Pt());
    else if(pdgMotherCombPart==2212)  fHistoCombinatorialMothersPt[fiCut]->Fill(0.,3.,motherCombPart->Pt());
    else if(pdgMotherCombPart==111)   fHistoCombinatorialMothersPt[fiCut]->Fill(0.,4.,motherCombPart->Pt());
    else if(pdgMotherCombPart==221)   fHistoCombinatorialMothersPt[fiCut]->Fill(0.,5.,motherCombPart->Pt());
    else if(pdgMotherCombPart==223)   fHistoCombinatorialMothersPt[fiCut]->Fill(0.,6.,motherCombPart->Pt());
    else if(pdgMotherCombPart==333)   fHistoCombinatorialMothersPt[fiCut]->Fill(0.,7.,motherCombPart->Pt());
    else if(pdgMotherCombPart==331)   fHistoCombinatorialMothersPt[fiCut]->Fill(0.,8.,motherCombPart->Pt());
    else if(pdgMotherCombPart==310)   fHistoCombinatorialMothersPt[fiCut]->Fill(0.,9.,motherCombPart->Pt());
    else if(pdgMotherCombPart==3122)  fHistoCombinatorialMothersPt[fiCut]->Fill(0.,10.,motherCombPart->Pt());
    else if(pdgMotherCombPart==113 || pdgMotherCombPart==213)fHistoCombinatorialMothersPt[fiCut]->Fill(0.,11.,motherCombPart->Pt());
    else fHistoCombinatorialMothersPt[fiCut]->Fill(0.,12.,motherCombPart->Pt());
  } else if(pdgCombPart==211){ //comb particle is a pion
    if(pdgMotherCombPart==11)         fHistoCombinatorialMothersPt[fiCut]->Fill(1.,0.,motherCombPart->Pt());
    else if(pdgMotherCombPart==211)   fHistoCombinatorialMothersPt[fiCut]->Fill(1.,1.,motherCombPart->Pt());
    else if(pdgMotherCombPart==321)   fHistoCombinatorialMothersPt[fiCut]->Fill(1.,2.,motherCombPart->Pt());
    else if(pdgMotherCombPart==2212)  fHistoCombinatorialMothersPt[fiCut]->Fill(1.,3.,motherCombPart->Pt());
    else if(pdgMotherCombPart==111)   fHistoCombinatorialMothersPt[fiCut]->Fill(1.,4.,motherCombPart->Pt());
    else if(pdgMotherCombPart==221)   fHistoCombinatorialMothersPt[fiCut]->Fill(1.,5.,motherCombPart->Pt());
    else if(pdgMotherCombPart==223)   fHistoCombinatorialMothersPt[fiCut]->Fill(1.,6.,motherCombPart->Pt());
    else if(pdgMotherCombPart==333)   fHistoCombinatorialMothersPt[fiCut]->Fill(1.,7.,motherCombPart->Pt());
    else if(pdgMotherCombPart==331)   fHistoCombinatorialMothersPt[fiCut]->Fill(1.,8.,motherCombPart->Pt());
    else if(pdgMotherCombPart==310)   fHistoCombinatorialMothersPt[fiCut]->Fill(1.,9.,motherCombPart->Pt());
    else if(pdgMotherCombPart==3122)  fHistoCombinatorialMothersPt[fiCut]->Fill(1.,10.,motherCombPart->Pt());
    else if(pdgMotherCombPart==113 || pdgMotherCombPart==213)fHistoCombinatorialMothersPt[fiCut]->Fill(1.,11.,motherCombPart->Pt());
    else fHistoCombinatorialMothersPt[fiCut]->Fill(1.,12.,motherCombPart->Pt());
  } else if(pdgCombPart==321){ //comb particle is a kaon
    if(pdgMotherCombPart==11)         fHistoCombinatorialMothersPt[fiCut]->Fill(2.,0.,motherCombPart->Pt());
    else if(pdgMotherCombPart==211)   fHistoCombinatorialMothersPt[fiCut]->Fill(2.,1.,motherCombPart->Pt());
    else if(pdgMotherCombPart==321)   fHistoCombinatorialMothersPt[fiCut]->Fill(2.,2.,motherCombPart->Pt());
    else if(pdgMotherCombPart==2212)  fHistoCombinatorialMothersPt[fiCut]->Fill(2.,3.,motherCombPart->Pt());
    else if(pdgMotherCombPart==111)   fHistoCombinatorialMothersPt[fiCut]->Fill(2.,4.,motherCombPart->Pt());
    else if(pdgMotherCombPart==221)   fHistoCombinatorialMothersPt[fiCut]->Fill(2.,5.,motherCombPart->Pt());
    else if(pdgMotherCombPart==223)   fHistoCombinatorialMothersPt[fiCut]->Fill(2.,6.,motherCombPart->Pt());
    else if(pdgMotherCombPart==333)   fHistoCombinatorialMothersPt[fiCut]->Fill(2.,7.,motherCombPart->Pt());
    else if(pdgMotherCombPart==331)   fHistoCombinatorialMothersPt[fiCut]->Fill(2.,8.,motherCombPart->Pt());
    else if(pdgMotherCombPart==310)   fHistoCombinatorialMothersPt[fiCut]->Fill(2.,9.,motherCombPart->Pt());
    else if(pdgMotherCombPart==3122)  fHistoCombinatorialMothersPt[fiCut]->Fill(2.,10.,motherCombPart->Pt());
    else if(pdgMotherCombPart==113 || pdgMotherCombPart==213)fHistoCombinatorialMothersPt[fiCut]->Fill(2.,11.,motherCombPart->Pt());
    else fHistoCombinatorialMothersPt[fiCut]->Fill(2.,12.,motherCombPart->Pt());
  } else if(pdgCombPart==2212){ //comb particle is a proton
    if(pdgMotherCombPart==11)         fHistoCombinatorialMothersPt[fiCut]->Fill(3.,0.,motherCombPart->Pt());
    else if(pdgMotherCombPart==211)   fHistoCombinatorialMothersPt[fiCut]->Fill(3.,1.,motherCombPart->Pt());
    else if(pdgMotherCombPart==321)   fHistoCombinatorialMothersPt[fiCut]->Fill(3.,2.,motherCombPart->Pt());
    else if(pdgMotherCombPart==2212)  fHistoCombinatorialMothersPt[fiCut]->Fill(3.,3.,motherCombPart->Pt());
    else if(pdgMotherCombPart==111)   fHistoCombinatorialMothersPt[fiCut]->Fill(3.,4.,motherCombPart->Pt());
    else if(pdgMotherCombPart==221)   fHistoCombinatorialMothersPt[fiCut]->Fill(3.,5.,motherCombPart->Pt());
    else if(pdgMotherCombPart==223)   fHistoCombinatorialMothersPt[fiCut]->Fill(3.,6.,motherCombPart->Pt());
    else if(pdgMotherCombPart==333)   fHistoCombinatorialMothersPt[fiCut]->Fill(3.,7.,motherCombPart->Pt());
    else if(pdgMotherCombPart==331)   fHistoCombinatorialMothersPt[fiCut]->Fill(3.,8.,motherCombPart->Pt());
    else if(pdgMotherCombPart==310)   fHistoCombinatorialMothersPt[fiCut]->Fill(3.,9.,motherCombPart->Pt());
    else if(pdgMotherCombPart==3122)  fHistoCombinatorialMothersPt[fiCut]->Fill(3.,10.,motherCombPart->Pt());
    else if(pdgMotherCombPart==113 || pdgMotherCombPart==213)fHistoCombinatorialMothersPt[fiCut]->Fill(3.,11.,motherCombPart->Pt());
    else fHistoCombinatorialMothersPt[fiCut]->Fill(3.,12.,motherCombPart->Pt());
  } else if(motherCombPart->GetMother(0) >-1){
    if(pdgMotherCombPart==11)         fHistoCombinatorialMothersPt[fiCut]->Fill(5.,0.,motherCombPart->Pt());
    else if(pdgMotherCombPart==211)   fHistoCombinatorialMothersPt[fiCut]->Fill(5.,1.,motherCombPart->Pt());
    else if(pdgMotherCombPart==321)   fHistoCombinatorialMothersPt[fiCut]->Fill(5.,2.,motherCombPart->Pt());
    else if(pdgMotherCombPart==2212)  fHistoCombinatorialMothersPt[fiCut]->Fill(5.,3.,motherCombPart->Pt());
    else if(pdgMotherCombPart==111)   fHistoCombinatorialMothersPt[fiCut]->Fill(5.,4.,motherCombPart->Pt());
    else if(pdgMotherCombPart==221)   fHistoCombinatorialMothersPt[fiCut]->Fill(5.,5.,motherCombPart->Pt());
    else if(pdgMotherCombPart==223)   fHistoCombinatorialMothersPt[fiCut]->Fill(5.,6.,motherCombPart->Pt());
    else if(pdgMotherCombPart==333)   fHistoCombinatorialMothersPt[fiCut]->Fill(5.,7.,motherCombPart->Pt());
    else if(pdgMotherCombPart==331)   fHistoCombinatorialMothersPt[fiCut]->Fill(5.,8.,motherCombPart->Pt());
    else if(pdgMotherCombPart==310)   fHistoCombinatorialMothersPt[fiCut]->Fill(5.,9.,motherCombPart->Pt());
    else if(pdgMotherCombPart==3122)  fHistoCombinatorialMothersPt[fiCut]->Fill(5.,10.,motherCombPart->Pt());
    else if(pdgMotherCombPart==113 || pdgMotherCombPart==213)fHistoCombinatorialMothersPt[fiCut]->Fill(5.,11.,motherCombPart->Pt());
    else fHistoCombinatorialMothersPt[fiCut]->Fill(5.,12.,motherCombPart->Pt());
  } else {
    if(pdgMotherCombPart==11)         fHistoCombinatorialMothersPt[fiCut]->Fill(4.,0.,motherCombPart->Pt());
    else if(pdgMotherCombPart==211)   fHistoCombinatorialMothersPt[fiCut]->Fill(4.,1.,motherCombPart->Pt());
    else if(pdgMotherCombPart==321)   fHistoCombinatorialMothersPt[fiCut]->Fill(4.,2.,motherCombPart->Pt());
    else if(pdgMotherCombPart==2212)  fHistoCombinatorialMothersPt[fiCut]->Fill(4.,3.,motherCombPart->Pt());
    else if(pdgMotherCombPart==111)   fHistoCombinatorialMothersPt[fiCut]->Fill(4.,4.,motherCombPart->Pt());
    else if(pdgMotherCombPart==221)   fHistoCombinatorialMothersPt[fiCut]->Fill(4.,5.,motherCombPart->Pt());
    else if(pdgMotherCombPart==223)   fHistoCombinatorialMothersPt[fiCut]->Fill(4.,6.,motherCombPart->Pt());
    else if(pdgMotherCombPart==333)   fHistoCombinatorialMothersPt[fiCut]->Fill(4.,7.,motherCombPart->Pt());
    else if(pdgMotherCombPart==331)   fHistoCombinatorialMothersPt[fiCut]->Fill(4.,8.,motherCombPart->Pt());
    else if(pdgMotherCombPart==310)   fHistoCombinatorialMothersPt[fiCut]->Fill(4.,9.,motherCombPart->Pt());
    else if(pdgMotherCombPart==3122)  fHistoCombinatorialMothersPt[fiCut]->Fill(4.,10.,motherCombPart->Pt());
    else if(pdgMotherCombPart==113 || pdgMotherCombPart==213)fHistoCombinatorialMothersPt[fiCut]->Fill(4.,11.,motherCombPart->Pt());
    else fHistoCombinatorialMothersPt[fiCut]->Fill(4.,12.,motherCombPart->Pt());
  }
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::FillPhotonCombinatorialMothersHistAOD(AliAODMCParticle *daughter, AliAODMCParticle* motherCombPart)
{
  Int_t pdgCombPart = TMath::Abs(daughter->GetPdgCode());
  Int_t pdgMotherCombPart = TMath::Abs(motherCombPart->GetPdgCode());

  // Combinatorial Bck mothers: e, pi, k, p, pi0, eta, omega, phi, eta', K0s, Lambda, rhos, other, not prim
  if(pdgCombPart==11){ //comb particle is an electron
    if(pdgMotherCombPart==11)         fHistoCombinatorialMothersPt[fiCut]->Fill(0.,0.,motherCombPart->Pt());
    else if(pdgMotherCombPart==211)   fHistoCombinatorialMothersPt[fiCut]->Fill(0.,1.,motherCombPart->Pt());
    else if(pdgMotherCombPart==321)   fHistoCombinatorialMothersPt[fiCut]->Fill(0.,2.,motherCombPart->Pt());
    else if(pdgMotherCombPart==2212)  fHistoCombinatorialMothersPt[fiCut]->Fill(0.,3.,motherCombPart->Pt());
    else if(pdgMotherCombPart==111)   fHistoCombinatorialMothersPt[fiCut]->Fill(0.,4.,motherCombPart->Pt());
    else if(pdgMotherCombPart==221)   fHistoCombinatorialMothersPt[fiCut]->Fill(0.,5.,motherCombPart->Pt());
    else if(pdgMotherCombPart==223)   fHistoCombinatorialMothersPt[fiCut]->Fill(0.,6.,motherCombPart->Pt());
    else if(pdgMotherCombPart==333)   fHistoCombinatorialMothersPt[fiCut]->Fill(0.,7.,motherCombPart->Pt());
    else if(pdgMotherCombPart==331)   fHistoCombinatorialMothersPt[fiCut]->Fill(0.,8.,motherCombPart->Pt());
    else if(pdgMotherCombPart==310)   fHistoCombinatorialMothersPt[fiCut]->Fill(0.,9.,motherCombPart->Pt());
    else if(pdgMotherCombPart==3122)  fHistoCombinatorialMothersPt[fiCut]->Fill(0.,10.,motherCombPart->Pt());
    else if(pdgMotherCombPart==113 || pdgMotherCombPart==213)fHistoCombinatorialMothersPt[fiCut]->Fill(0.,11.,motherCombPart->Pt());
    else fHistoCombinatorialMothersPt[fiCut]->Fill(0.,12.,motherCombPart->Pt());
  } else if(pdgCombPart==211){ //comb particle is a pion
    if(pdgMotherCombPart==11)         fHistoCombinatorialMothersPt[fiCut]->Fill(1.,0.,motherCombPart->Pt());
    else if(pdgMotherCombPart==211)   fHistoCombinatorialMothersPt[fiCut]->Fill(1.,1.,motherCombPart->Pt());
    else if(pdgMotherCombPart==321)   fHistoCombinatorialMothersPt[fiCut]->Fill(1.,2.,motherCombPart->Pt());
    else if(pdgMotherCombPart==2212)  fHistoCombinatorialMothersPt[fiCut]->Fill(1.,3.,motherCombPart->Pt());
    else if(pdgMotherCombPart==111)   fHistoCombinatorialMothersPt[fiCut]->Fill(1.,4.,motherCombPart->Pt());
    else if(pdgMotherCombPart==221)   fHistoCombinatorialMothersPt[fiCut]->Fill(1.,5.,motherCombPart->Pt());
    else if(pdgMotherCombPart==223)   fHistoCombinatorialMothersPt[fiCut]->Fill(1.,6.,motherCombPart->Pt());
    else if(pdgMotherCombPart==333)   fHistoCombinatorialMothersPt[fiCut]->Fill(1.,7.,motherCombPart->Pt());
    else if(pdgMotherCombPart==331)   fHistoCombinatorialMothersPt[fiCut]->Fill(1.,8.,motherCombPart->Pt());
    else if(pdgMotherCombPart==310)   fHistoCombinatorialMothersPt[fiCut]->Fill(1.,9.,motherCombPart->Pt());
    else if(pdgMotherCombPart==3122)  fHistoCombinatorialMothersPt[fiCut]->Fill(1.,10.,motherCombPart->Pt());
    else if(pdgMotherCombPart==113 || pdgMotherCombPart==213)fHistoCombinatorialMothersPt[fiCut]->Fill(1.,11.,motherCombPart->Pt());
    else fHistoCombinatorialMothersPt[fiCut]->Fill(1.,12.,motherCombPart->Pt());
  } else if(pdgCombPart==321){ //comb particle is a kaon
    if(pdgMotherCombPart==11)         fHistoCombinatorialMothersPt[fiCut]->Fill(2.,0.,motherCombPart->Pt());
    else if(pdgMotherCombPart==211)   fHistoCombinatorialMothersPt[fiCut]->Fill(2.,1.,motherCombPart->Pt());
    else if(pdgMotherCombPart==321)   fHistoCombinatorialMothersPt[fiCut]->Fill(2.,2.,motherCombPart->Pt());
    else if(pdgMotherCombPart==2212)  fHistoCombinatorialMothersPt[fiCut]->Fill(2.,3.,motherCombPart->Pt());
    else if(pdgMotherCombPart==111)   fHistoCombinatorialMothersPt[fiCut]->Fill(2.,4.,motherCombPart->Pt());
    else if(pdgMotherCombPart==221)   fHistoCombinatorialMothersPt[fiCut]->Fill(2.,5.,motherCombPart->Pt());
    else if(pdgMotherCombPart==223)   fHistoCombinatorialMothersPt[fiCut]->Fill(2.,6.,motherCombPart->Pt());
    else if(pdgMotherCombPart==333)   fHistoCombinatorialMothersPt[fiCut]->Fill(2.,7.,motherCombPart->Pt());
    else if(pdgMotherCombPart==331)   fHistoCombinatorialMothersPt[fiCut]->Fill(2.,8.,motherCombPart->Pt());
    else if(pdgMotherCombPart==310)   fHistoCombinatorialMothersPt[fiCut]->Fill(2.,9.,motherCombPart->Pt());
    else if(pdgMotherCombPart==3122)  fHistoCombinatorialMothersPt[fiCut]->Fill(2.,10.,motherCombPart->Pt());
    else if(pdgMotherCombPart==113 || pdgMotherCombPart==213)fHistoCombinatorialMothersPt[fiCut]->Fill(2.,11.,motherCombPart->Pt());
    else fHistoCombinatorialMothersPt[fiCut]->Fill(2.,12.,motherCombPart->Pt());
  } else if(pdgCombPart==2212){ //comb particle is a proton
    if(pdgMotherCombPart==11)         fHistoCombinatorialMothersPt[fiCut]->Fill(3.,0.,motherCombPart->Pt());
    else if(pdgMotherCombPart==211)   fHistoCombinatorialMothersPt[fiCut]->Fill(3.,1.,motherCombPart->Pt());
    else if(pdgMotherCombPart==321)   fHistoCombinatorialMothersPt[fiCut]->Fill(3.,2.,motherCombPart->Pt());
    else if(pdgMotherCombPart==2212)  fHistoCombinatorialMothersPt[fiCut]->Fill(3.,3.,motherCombPart->Pt());
    else if(pdgMotherCombPart==111)   fHistoCombinatorialMothersPt[fiCut]->Fill(3.,4.,motherCombPart->Pt());
    else if(pdgMotherCombPart==221)   fHistoCombinatorialMothersPt[fiCut]->Fill(3.,5.,motherCombPart->Pt());
    else if(pdgMotherCombPart==223)   fHistoCombinatorialMothersPt[fiCut]->Fill(3.,6.,motherCombPart->Pt());
    else if(pdgMotherCombPart==333)   fHistoCombinatorialMothersPt[fiCut]->Fill(3.,7.,motherCombPart->Pt());
    else if(pdgMotherCombPart==331)   fHistoCombinatorialMothersPt[fiCut]->Fill(3.,8.,motherCombPart->Pt());
    else if(pdgMotherCombPart==310)   fHistoCombinatorialMothersPt[fiCut]->Fill(3.,9.,motherCombPart->Pt());
    else if(pdgMotherCombPart==3122)  fHistoCombinatorialMothersPt[fiCut]->Fill(3.,10.,motherCombPart->Pt());
    else if(pdgMotherCombPart==113 || pdgMotherCombPart==213)fHistoCombinatorialMothersPt[fiCut]->Fill(3.,11.,motherCombPart->Pt());
    else fHistoCombinatorialMothersPt[fiCut]->Fill(3.,12.,motherCombPart->Pt());
  } else if(motherCombPart->GetMother() >-1){
    if(pdgMotherCombPart==11)         fHistoCombinatorialMothersPt[fiCut]->Fill(5.,0.,motherCombPart->Pt());
    else if(pdgMotherCombPart==211)   fHistoCombinatorialMothersPt[fiCut]->Fill(5.,1.,motherCombPart->Pt());
    else if(pdgMotherCombPart==321)   fHistoCombinatorialMothersPt[fiCut]->Fill(5.,2.,motherCombPart->Pt());
    else if(pdgMotherCombPart==2212)  fHistoCombinatorialMothersPt[fiCut]->Fill(5.,3.,motherCombPart->Pt());
    else if(pdgMotherCombPart==111)   fHistoCombinatorialMothersPt[fiCut]->Fill(5.,4.,motherCombPart->Pt());
    else if(pdgMotherCombPart==221)   fHistoCombinatorialMothersPt[fiCut]->Fill(5.,5.,motherCombPart->Pt());
    else if(pdgMotherCombPart==223)   fHistoCombinatorialMothersPt[fiCut]->Fill(5.,6.,motherCombPart->Pt());
    else if(pdgMotherCombPart==333)   fHistoCombinatorialMothersPt[fiCut]->Fill(5.,7.,motherCombPart->Pt());
    else if(pdgMotherCombPart==331)   fHistoCombinatorialMothersPt[fiCut]->Fill(5.,8.,motherCombPart->Pt());
    else if(pdgMotherCombPart==310)   fHistoCombinatorialMothersPt[fiCut]->Fill(5.,9.,motherCombPart->Pt());
    else if(pdgMotherCombPart==3122)  fHistoCombinatorialMothersPt[fiCut]->Fill(5.,10.,motherCombPart->Pt());
    else if(pdgMotherCombPart==113 || pdgMotherCombPart==213)fHistoCombinatorialMothersPt[fiCut]->Fill(5.,11.,motherCombPart->Pt());
    else fHistoCombinatorialMothersPt[fiCut]->Fill(5.,12.,motherCombPart->Pt());
  } else {
    if(pdgMotherCombPart==11)         fHistoCombinatorialMothersPt[fiCut]->Fill(4.,0.,motherCombPart->Pt());
    else if(pdgMotherCombPart==211)   fHistoCombinatorialMothersPt[fiCut]->Fill(4.,1.,motherCombPart->Pt());
    else if(pdgMotherCombPart==321)   fHistoCombinatorialMothersPt[fiCut]->Fill(4.,2.,motherCombPart->Pt());
    else if(pdgMotherCombPart==2212)  fHistoCombinatorialMothersPt[fiCut]->Fill(4.,3.,motherCombPart->Pt());
    else if(pdgMotherCombPart==111)   fHistoCombinatorialMothersPt[fiCut]->Fill(4.,4.,motherCombPart->Pt());
    else if(pdgMotherCombPart==221)   fHistoCombinatorialMothersPt[fiCut]->Fill(4.,5.,motherCombPart->Pt());
    else if(pdgMotherCombPart==223)   fHistoCombinatorialMothersPt[fiCut]->Fill(4.,6.,motherCombPart->Pt());
    else if(pdgMotherCombPart==333)   fHistoCombinatorialMothersPt[fiCut]->Fill(4.,7.,motherCombPart->Pt());
    else if(pdgMotherCombPart==331)   fHistoCombinatorialMothersPt[fiCut]->Fill(4.,8.,motherCombPart->Pt());
    else if(pdgMotherCombPart==310)   fHistoCombinatorialMothersPt[fiCut]->Fill(4.,9.,motherCombPart->Pt());
    else if(pdgMotherCombPart==3122)  fHistoCombinatorialMothersPt[fiCut]->Fill(4.,10.,motherCombPart->Pt());
    else if(pdgMotherCombPart==113 || pdgMotherCombPart==213)fHistoCombinatorialMothersPt[fiCut]->Fill(4.,11.,motherCombPart->Pt());
    else fHistoCombinatorialMothersPt[fiCut]->Fill(4.,12.,motherCombPart->Pt());
  }
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::RelabelAODPhotonCandidates(Bool_t mode){

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
void AliAnalysisTaskGammaConvV1::SetLogBinningXTH2(TH2* histoRebin){
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
Bool_t AliAnalysisTaskGammaConvV1::CheckVectorForDoubleCount(vector<Int_t> &vec, Int_t tobechecked)
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
void AliAnalysisTaskGammaConvV1::FillMultipleCountMap(map<Int_t,Int_t> &ma, Int_t tobechecked){
  if( ma.find(tobechecked) != ma.end() ) ma[tobechecked] += 1;
  else ma[tobechecked] = 2;
  return;
}

//_________________________________________________________________________________
void AliAnalysisTaskGammaConvV1::FillMultipleCountHistoAndClear(map<Int_t,Int_t> &ma, TH1F* hist){
  map<Int_t, Int_t>::iterator it;
  for (it = ma.begin(); it != ma.end(); it++){
    hist->Fill(it->second,fWeightJetJetMC);
  }
  ma.clear();
  return;
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::Terminate(const Option_t *)
{

  //fOutputContainer->Print(); // Will crash on GRID
}

//________________________________________________________________________
Int_t AliAnalysisTaskGammaConvV1::GetSourceClassification(Int_t daughter, Int_t pdgCode){

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

//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::ProcessClusters(){

  Int_t nclus = 0;
  nclus = fInputEvent->GetNumberOfCaloClusters();

//   cout << nclus << endl;

  if(nclus == 0)  return;

  // vertex
  Double_t vertex[3] = {0,0,0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  // Loop over EMCal clusters
  for(Int_t i = 0; i < nclus; i++){

    AliVCluster* clus = NULL;
    if(fInputEvent->IsA()==AliESDEvent::Class()) clus = new AliESDCaloCluster(*(AliESDCaloCluster*)fInputEvent->GetCaloCluster(i));
    else if(fInputEvent->IsA()==AliAODEvent::Class()) clus = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(i));

    if (!clus) continue;
    if(!((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC,fWeightJetJetMC,i)){ delete clus; continue;}

    // TLorentzvector with cluster
    TLorentzVector clusterVector;
    clus->GetMomentum(clusterVector,vertex);

    TLorentzVector* tmpvec = new TLorentzVector();
    tmpvec->SetPxPyPzE(clusterVector.Px(),clusterVector.Py(),clusterVector.Pz(),clusterVector.E());

    // convert to AODConversionPhoton
    AliAODConversionPhoton *PhotonCandidate = new AliAODConversionPhoton(tmpvec);
    if(!PhotonCandidate){ delete clus; delete tmpvec; continue;}

    // Flag Photon as CaloPhoton
    PhotonCandidate->SetIsCaloPhoton();
    PhotonCandidate->SetCaloClusterRef((Long_t)i);

    fHistoCaloGammaPt[fiCut]->Fill(PhotonCandidate->Pt(),fWeightJetJetMC);
    fHistoCaloGammaE[fiCut]->Fill(PhotonCandidate->E(),fWeightJetJetMC);
    delete PhotonCandidate;
    delete clus;
    delete tmpvec;
  }
}

///________________________________________________________________________
Double_t AliAnalysisTaskGammaConvV1::GetOriginalInvMass( const AliConversionPhotonBase * photon, AliVEvent * event) const{
  // invariant mass of original V0 before recalculation
  Double_t invMass = -999;
  if(event->IsA()==AliESDEvent::Class()){
    AliESDEvent *esdEvent = dynamic_cast<AliESDEvent*>(event);
    if(!esdEvent) return -999;
    AliESDv0 *v0 = esdEvent->GetV0(photon->GetV0Index());   // Get V0 from ESD V0 array using current photons index
    if(!v0) return -999;
    const float emass = 0.55e-3;
    invMass = v0->GetEffMassExplicit(emass,emass);
  }
  return invMass;
}
