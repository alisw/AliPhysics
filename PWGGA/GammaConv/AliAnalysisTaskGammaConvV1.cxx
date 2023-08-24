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
#include "AliGAKFVertex.h"
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
  fDoIsolatedAnalysis(kFALSE),
  fDoHighPtHadronAnalysis(kFALSE),
  fHistoCaloGammaPt(NULL),
  fHistoCaloGammaE(NULL),
  fHistoConvGammaPt(NULL),
  fHistoConvGammaPtwithHighPtHadron(NULL),
  fHistoConvGammaPtwithoutHighPtHadron(NULL),
  fHistoNEventsHighPtHadron(NULL),
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
  fHistoMotherInvMassPtIso(NULL),
  fHistoMotherInvMassPtNonIso(NULL),
  fHistoMotherInvMassPtMCRecIsoTrueNonIso(NULL),
  fHistoMotherInvMassPtMCRecNonIsoTrueIso(NULL),
  fHistoMotherInvMassPtCalib(NULL),
  fHistoMotherEisoPt(NULL),
  fHistoMotherRisoPt(NULL),
  fHistoMotherNtracksIsoPt(NULL),
  sESDMotherInvMassPtZM(NULL),
  fHistoMotherBackInvMassPt(NULL),
  fHistoMotherBackInvMassPtCalib(NULL),
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
  fHistoMCAllGammaPtNotTriggered(NULL),
  fHistoMCAllGammaPtNoVertex(NULL),
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
  fHistoMCPi0PtNotTriggered(NULL),
  fHistoMCPi0PtNoVertex(NULL),
  fHistoMCPi0WOWeightPt(NULL),
  fHistoMCPi0WOEvtWeightPt(NULL),
  fHistoMCEtaWOEvtWeightPt(NULL),
  fHistoMCEtaPt(NULL),
  fHistoMCEtaPtNotTriggered(NULL),
  fHistoMCEtaPtNoVertex(NULL),
  fHistoMCEtaWOWeightPt(NULL),
  fHistoMCPi0WOWeightInAccPt(NULL),
  fHistoMCEtaWOWeightInAccPt(NULL),
  fHistoMCPi0InAccPt(NULL),
  fHistoMCPi0InAccPtNotTriggered(NULL),
  fHistoMCPi0InAccPtNoVertex(NULL),
  fHistoMCEtaInAccPt(NULL),
  fHistoMCEtaInAccPtNotTriggered(NULL),
  fHistoMCEtaInAccPtNoVertex(NULL),
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
  fHistoSPDClusterTrackletBackground(NULL),
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
  fiEventCut(NULL),
  fiPhotonCut(NULL),
  fiMesonCut(NULL),
  fMoveParticleAccordingToVertex(kTRUE),
  fIsHeavyIon(0),
  fDoMesonAnalysis(kTRUE),
  fDoMesonQA(0),
  fDoPhotonQA(0),
  fDoChargedPrimary(kFALSE),
  fDoPlotVsCentrality(kFALSE),
  fIsMC(0),
  fDoTHnSparse(kFALSE),
  fWeightJetJetMC(1),
  fWeightCentrality(NULL),
  fEnableClusterCutsForTrigger(kFALSE),
  fDoMaterialBudgetWeightingOfGammasForTrueMesons(kFALSE),
  tBrokenFiles(NULL),
  fFileNameBroken(NULL),
  fFileWasAlreadyReported(kFALSE),
  fAODMCTrackArray(NULL),
  fAddressChanges(NULL),
  fMapPhotonHeaders()
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
  fDoIsolatedAnalysis(kFALSE),
  fDoHighPtHadronAnalysis(kFALSE),
  fHistoCaloGammaPt(NULL),
  fHistoCaloGammaE(NULL),
  fHistoConvGammaPt(NULL),
  fHistoConvGammaPtwithHighPtHadron(NULL),
  fHistoConvGammaPtwithoutHighPtHadron(NULL),
  fHistoNEventsHighPtHadron(NULL),
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
  fHistoMotherInvMassPtIso(NULL),
  fHistoMotherInvMassPtNonIso(NULL),
  fHistoMotherInvMassPtMCRecIsoTrueNonIso(NULL),
  fHistoMotherInvMassPtMCRecNonIsoTrueIso(NULL),
  fHistoMotherInvMassPtCalib(NULL),
  fHistoMotherEisoPt(NULL),
  fHistoMotherRisoPt(NULL),
  fHistoMotherNtracksIsoPt(NULL),
  sESDMotherInvMassPtZM(NULL),
  fHistoMotherBackInvMassPt(NULL),
  fHistoMotherBackInvMassPtCalib(NULL),
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
  fHistoMCAllGammaPtNotTriggered(NULL),
  fHistoMCAllGammaPtNoVertex(NULL),
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
  fHistoMCPi0PtNotTriggered(NULL),
  fHistoMCPi0PtNoVertex(NULL),
  fHistoMCPi0WOWeightPt(NULL),
  fHistoMCPi0WOEvtWeightPt(NULL),
  fHistoMCEtaWOEvtWeightPt(NULL),
  fHistoMCEtaPt(NULL),
  fHistoMCEtaPtNotTriggered(NULL),
  fHistoMCEtaPtNoVertex(NULL),
  fHistoMCEtaWOWeightPt(NULL),
  fHistoMCPi0WOWeightInAccPt(NULL),
  fHistoMCEtaWOWeightInAccPt(NULL),
  fHistoMCPi0InAccPt(NULL),
  fHistoMCPi0InAccPtNotTriggered(NULL),
  fHistoMCPi0InAccPtNoVertex(NULL),
  fHistoMCEtaInAccPt(NULL),
  fHistoMCEtaInAccPtNotTriggered(NULL),
  fHistoMCEtaInAccPtNoVertex(NULL),
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
  fHistoSPDClusterTrackletBackground(NULL),
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
  fiEventCut(NULL),
  fiPhotonCut(NULL),
  fiMesonCut(NULL),
  fMoveParticleAccordingToVertex(kTRUE),
  fIsHeavyIon(0),
  fDoMesonAnalysis(kTRUE),
  fDoMesonQA(0),
  fDoPhotonQA(0),
  fDoChargedPrimary(kFALSE),
  fDoPlotVsCentrality(kFALSE),
  fIsMC(0),
  fDoTHnSparse(kFALSE),
  fWeightJetJetMC(1),
  fWeightCentrality(NULL),
  fEnableClusterCutsForTrigger(kFALSE),
  fDoMaterialBudgetWeightingOfGammasForTrueMesons(kFALSE),
  tBrokenFiles(NULL),
  fFileNameBroken(NULL),
  fFileWasAlreadyReported(kFALSE),
  fAODMCTrackArray(NULL),
  fAddressChanges(NULL),
  fMapPhotonHeaders()
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
      } else if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->BackgroundHandlerType() != 2){
        fBGHandlerRP[iCut] = new AliConversionAODBGHandlerRP(
                                  ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsHeavyIon(),
                                  ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseTrackMultiplicity(),
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

  if(fDoMesonAnalysis){ 
    if( ((AliConversionMesonCuts*)fMesonCutArray->At(0))->DoIsolatedAnalysis())     fDoIsolatedAnalysis = kTRUE;
  }
  if( ((AliConversionMesonCuts*)fMesonCutArray->At(0))->DoHighPtHadronAnalysis()) fDoHighPtHadronAnalysis = kTRUE;

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


  // Set special pt binning for pp 13TeV and 5TeV
  if ( ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::k13TeV ||
       ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::k13TeVLowB ){
    // binWidthPt                = 0.05;
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
  } else if ( ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::kpPb5TeV ||
              ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::kpPb5TeVR2 ||
              ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::k5TeV ||
              ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::k8TeV  ||
              ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::kpPb8TeV  ){
    nBinsPt                   = 159;
    minPt                     = 0;
    maxPt                     = 40;
    for(Int_t i=0; i<nBinsPt+1;i++){
      if (i < 100) arrPtBinning[i]            = 0.10*i;
      else if(i<140) arrPtBinning[i]          = 10.+0.25*(i-100);
      else if(i<160) arrPtBinning[i]          = 20.+1.0*(i-140);
      else arrPtBinning[i]                    = maxPt;
    }
    nBinsQAPt                 = 159;
    maxQAPt                   = 40;
    for(Int_t i=0; i<nBinsQAPt+1;i++){
      if(i<100) arrQAPtBinning[i]             = 0.1*i;
      else if(i<140) arrQAPtBinning[i]        = 10.+0.25*(i-100);
      else if(i<160) arrQAPtBinning[i]        = 20.+1.0*(i-140);
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
  if(fDoHighPtHadronAnalysis){
    fHistoConvGammaPtwithHighPtHadron        = new TH2F*[fnCuts];
    fHistoConvGammaPtwithoutHighPtHadron     = new TH2F*[fnCuts];
    fHistoNEventsHighPtHadron                = new TH1F*[fnCuts];
  }
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
    if(!fDoLightOutput){
      fHistoMotherInvMassPtCalib    = new TH2F*[fnCuts];
      fHistoMotherBackInvMassPtCalib= new TH2F*[fnCuts];
    }
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
  if(fDoIsolatedAnalysis){
    fHistoMotherInvMassPtIso        = new TH2F*[fnCuts];
    fHistoMotherInvMassPtNonIso     = new TH2F*[fnCuts];
    fHistoMotherEisoPt              = new TH2F*[fnCuts];
    fHistoMotherRisoPt              = new TH2F*[fnCuts];
    fHistoMotherNtracksIsoPt        = new TH2F*[fnCuts];
    if (fIsMC > 0){
      fHistoMotherInvMassPtMCRecIsoTrueNonIso     = new TH2F*[fnCuts];
      fHistoMotherInvMassPtMCRecNonIsoTrueIso     = new TH2F*[fnCuts];
    }
  }
  if (fEnableClusterCutsForTrigger){
    fHistoCaloGammaPt           = new TH1F*[fnCuts];
    fHistoCaloGammaE            = new TH1F*[fnCuts];
  }


  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    AliConvEventCuts* iEventCut = dynamic_cast<AliConvEventCuts*>(fEventCutArray->At(iCut));

    TString cutstringEvent      = iEventCut->GetCutNumber();
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
    if (iEventCut->IsSpecialTrigger() > 1 ){
      TString TriggerNames      = "Not Trigger: ";
      TriggerNames              = TriggerNames+ iEventCut->GetSpecialTriggerName();
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
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(11,"rejectedForJetJetMC");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(14,iEventCut->GetLabelNamePileupCutTPC().Data());
    fESDList[iCut]->Add(fHistoNEvents[iCut]);
    if (fIsMC > 1){
      fHistoNEventsWOWeight[iCut]    = new TH1F("NEventsWOWeight", "NEventsWOWeight", 14, -0.5, 13.5);
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(3,"Miss. MC or inc. ev.");
      if (iEventCut->IsSpecialTrigger() > 1 ){
        TString TriggerNames    = "Not Trigger: ";
        TriggerNames            = TriggerNames+ iEventCut->GetSpecialTriggerName();
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
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(14,iEventCut->GetLabelNamePileupCutTPC().Data());
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
      if (iEventCut->IsSpecialTrigger() > 1 ){
        TString TriggerNames        = "Not Trigger: ";
        TriggerNames                = TriggerNames+ iEventCut->GetSpecialTriggerName();
        fHistoNEventsWeighted[iCut]->GetXaxis()->SetBinLabel(4,TriggerNames.Data());
      } else {
        fHistoNEventsWeighted[iCut]->GetXaxis()->SetBinLabel(4,"Trigger");
      }
      fHistoNEventsWeighted[iCut]->GetXaxis()->SetBinLabel(5,"Vertex Z");
      fHistoNEventsWeighted[iCut]->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
      fHistoNEventsWeighted[iCut]->GetXaxis()->SetBinLabel(7,"Pile-Up");
      fHistoNEventsWeighted[iCut]->GetXaxis()->SetBinLabel(8,"no SDD");
      fHistoNEventsWeighted[iCut]->GetXaxis()->SetBinLabel(9,"no V0AND");
      fHistoNEventsWeighted[iCut]->GetXaxis()->SetBinLabel(10,"EMCAL/TPC problem");
      fHistoNEventsWeighted[iCut]->GetXaxis()->SetBinLabel(11,"rejectedForJetJetMC");
      fHistoNEventsWeighted[iCut]->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
      fHistoNEventsWeighted[iCut]->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
      fHistoNEventsWeighted[iCut]->GetXaxis()->SetBinLabel(14,iEventCut->GetLabelNamePileupCutTPC().Data());
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

    if(iEventCut->GetUseSphericity()!=0){
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
    if(fDoHighPtHadronAnalysis){
      fHistoConvGammaPtwithHighPtHadron[iCut]         = new TH2F("ESD_ConvGamma_Pt_withHighPtHadron", "ESD_ConvGamma_Pt_withHighPtHadron", nBinsPt, arrPtBinning, nTracks, 0, nTracks);
      fESDList[iCut]->Add(fHistoConvGammaPtwithHighPtHadron[iCut]);
      fHistoConvGammaPtwithoutHighPtHadron[iCut]         = new TH2F("ESD_ConvGamma_Pt_withoutHighPtHadron", "ESD_ConvGamma_Pt_withoutHighPtHadron", nBinsPt, arrPtBinning, nTracks, 0, nTracks);
      fESDList[iCut]->Add(fHistoConvGammaPtwithoutHighPtHadron[iCut]);
      fHistoNEventsHighPtHadron[iCut]            = new TH1F("NEventsHighPtHadron", "NEventsHighPtHadron", 2, -0.5, 1.5);
      fHistoNEventsHighPtHadron[iCut]->GetXaxis()->SetBinLabel(1,"With");
      fHistoNEventsHighPtHadron[iCut]->GetXaxis()->SetBinLabel(2,"Without");
      fESDList[iCut]->Add(fHistoNEventsHighPtHadron[iCut]);
    }

    if (  (fIsMC > 1) || (fIsMC>0 && fDoMaterialBudgetWeightingOfGammasForTrueMesons) ) {
      fHistoConvGammaPt[iCut]->Sumw2();
      if(fDoHighPtHadronAnalysis){
        fHistoConvGammaPtwithHighPtHadron[iCut]->Sumw2();
        fHistoConvGammaPtwithoutHighPtHadron[iCut]->Sumw2();
        fHistoNEventsHighPtHadron[iCut]->Sumw2();
      }
    }

    if (fIsMC > 1){
      fHistoNEvents[iCut]->Sumw2();
      fHistoNGoodESDTracks[iCut]->Sumw2();
      fHistoVertexZ[iCut]->Sumw2();
      fHistoEtaShift[iCut]->Sumw2();
      if(fDoPlotVsCentrality){
        fHistoCentrality[iCut]->Sumw2();
      }
      if(iEventCut->GetUseSphericity()!=0){
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

      if(!fDoLightOutput){
        fHistoMotherInvMassPtCalib[iCut]   = new TH2F("ESD_Mother_InvMass_Pt_Calib", "ESD_Mother_InvMass_Pt_Calib", 800, 0, 0.8, nBinsPt, arrPtBinning);
        fESDList[iCut]->Add(fHistoMotherInvMassPtCalib[iCut]);
        fHistoMotherBackInvMassPtCalib[iCut]   = new TH2F("ESD_Background_InvMass_Pt_Calib", "ESD_Background_InvMass_Pt_Calib", 800, 0, 0.8, nBinsPt, arrPtBinning);
        fESDList[iCut]->Add(fHistoMotherBackInvMassPtCalib[iCut]);
      }

      if (fIsMC > 1 || fDoCentralityFlat > 0){
        fHistoMotherInvMassPt[iCut]->Sumw2();
        fHistoMotherBackInvMassPt[iCut]->Sumw2();
        fHistoMotherInvMassEalpha[iCut]->Sumw2();
        if(!fDoLightOutput){
          fHistoMotherInvMassPtCalib[iCut]->Sumw2();
          fHistoMotherBackInvMassPtCalib[iCut]->Sumw2();
        }
      }

      if(fDoIsolatedAnalysis){
        fHistoMotherInvMassPtIso[iCut]   = new TH2F("ESD_Mother_InvMass_Pt_Iso", "ESD_Mother_InvMass_Pt_Iso", 800, 0, 0.8, nBinsPt, arrPtBinning);
        fESDList[iCut]->Add(fHistoMotherInvMassPtIso[iCut]);
        fHistoMotherInvMassPtNonIso[iCut]   = new TH2F("ESD_Mother_InvMass_Pt_NonIso", "ESD_Mother_InvMass_Pt_NonIso", 800, 0, 0.8, nBinsPt, arrPtBinning);
        fESDList[iCut]->Add(fHistoMotherInvMassPtNonIso[iCut]);
        fHistoMotherEisoPt[iCut]   = new TH2F("ESD_Mother_EIso_Pt", "ESD_Mother_EIso_Pt", 400, 0, 40, nBinsPt, arrPtBinning);
        fESDList[iCut]->Add(fHistoMotherEisoPt[iCut]);
        fHistoMotherRisoPt[iCut]   = new TH2F("ESD_Mother_RIso_Pt", "ESD_Mother_RIso_Pt", 800, 0, 8, nBinsPt, arrPtBinning);
        fESDList[iCut]->Add(fHistoMotherRisoPt[iCut]);
        fHistoMotherNtracksIsoPt[iCut]   = new TH2F("ESD_Mother_NTracksIso_Pt", "ESD_Mother_NTracksIso_Pt", 100, 0, 100, nBinsPt, arrPtBinning);
        fESDList[iCut]->Add(fHistoMotherNtracksIsoPt[iCut]);
        if (fIsMC > 0 ){
          fHistoMotherInvMassPtMCRecIsoTrueNonIso[iCut]   = new TH2F("ESD_Mother_InvMass_Pt_MCRecIsoTrueNonIso", "ESD_Mother_InvMass_Pt_MCRecIsoTrueNonIso", 800, 0, 0.8, nBinsPt, arrPtBinning);
          fESDList[iCut]->Add(fHistoMotherInvMassPtMCRecIsoTrueNonIso[iCut]);
          fHistoMotherInvMassPtMCRecNonIsoTrueIso[iCut]   = new TH2F("ESD_Mother_InvMass_Pt_MCRecNonIsoTrueIso", "ESD_Mother_InvMass_Pt_MCRecNonIsoTrueIso", 800, 0, 0.8, nBinsPt, arrPtBinning);
          fESDList[iCut]->Add(fHistoMotherInvMassPtMCRecNonIsoTrueIso[iCut]);
        }
        if (fIsMC > 1 ){
          fHistoMotherInvMassPtIso[iCut]->Sumw2();
          fHistoMotherInvMassPtNonIso[iCut]->Sumw2();
          fHistoMotherEisoPt[iCut]->Sumw2();
          fHistoMotherRisoPt[iCut]->Sumw2();
          fHistoMotherNtracksIsoPt[iCut]->Sumw2();
          fHistoMotherInvMassPtMCRecIsoTrueNonIso[iCut]->Sumw2();
          fHistoMotherInvMassPtMCRecNonIsoTrueIso[iCut]->Sumw2();
        }
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
    fHistoMCAllGammaPtNotTriggered     = new TH1F*[fnCuts];
    fHistoMCAllGammaPtNoVertex         = new TH1F*[fnCuts];
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
      fHistoMCPi0PtNotTriggered        = new TH1F*[fnCuts];
      fHistoMCPi0PtNoVertex            = new TH1F*[fnCuts];
      fHistoMCPi0WOWeightPt            = new TH1F*[fnCuts];
      fHistoMCEtaPt                    = new TH1F*[fnCuts];
      fHistoMCEtaPtNotTriggered        = new TH1F*[fnCuts];
      fHistoMCEtaPtNoVertex            = new TH1F*[fnCuts];
      fHistoMCEtaWOWeightPt            = new TH1F*[fnCuts];
      fHistoMCPi0WOWeightInAccPt       = new TH1F*[fnCuts];
      fHistoMCEtaWOWeightInAccPt       = new TH1F*[fnCuts];
      fHistoMCPi0InAccPt               = new TH1F*[fnCuts];
      fHistoMCPi0InAccPtNotTriggered   = new TH1F*[fnCuts];
      fHistoMCPi0InAccPtNoVertex       = new TH1F*[fnCuts];
      fHistoMCEtaInAccPt               = new TH1F*[fnCuts];
      fHistoMCEtaInAccPtNotTriggered   = new TH1F*[fnCuts];
      fHistoMCEtaInAccPtNoVertex       = new TH1F*[fnCuts];

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
      fHistoMCAllGammaPtNotTriggered[iCut]             = new TH1F("MC_AllGamma_Pt_NotTriggered", "MC_AllGamma_Pt_NotTriggered", nBinsPt, arrPtBinning);
      fMCList[iCut]->Add(fHistoMCAllGammaPtNotTriggered[iCut]);
      fHistoMCAllGammaPtNoVertex[iCut]             = new TH1F("MC_AllGamma_Pt_NoVertex", "MC_AllGamma_Pt_NoVertex", nBinsPt, arrPtBinning);
      fMCList[iCut]->Add(fHistoMCAllGammaPtNoVertex[iCut]);
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

        fHistoMCPi0PtNotTriggered[iCut]     = new TH1F("MC_Pi0_Pt_NotTriggered", "MC_Pi0_Pt_NotTriggered", nBinsPt, arrPtBinning);
        fHistoMCPi0PtNotTriggered[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0PtNotTriggered[iCut]);

        fHistoMCPi0PtNoVertex[iCut]     = new TH1F("MC_Pi0_Pt_NoVertex", "MC_Pi0_Pt_NoVertex", nBinsPt, arrPtBinning);
        fHistoMCPi0PtNoVertex[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0PtNoVertex[iCut]);

        fHistoMCPi0WOWeightPt[iCut]        = new TH1F("MC_Pi0_WOWeights_Pt", "MC_Pi0_WOWeights_Pt", nBinsPt, arrPtBinning);
        fHistoMCPi0WOWeightPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0WOWeightPt[iCut]);

        fHistoMCEtaPtNotTriggered[iCut]     = new TH1F("MC_Eta_Pt_NotTriggered", "MC_Eta_Pt_NotTriggered", nBinsPt, arrPtBinning);
        fHistoMCEtaPtNotTriggered[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCEtaPtNotTriggered[iCut]);

        fHistoMCEtaPtNoVertex[iCut]        = new TH1F("MC_Eta_Pt_NoVertex", "MC_Eta_Pt_NoVertex", nBinsPt, arrPtBinning);
        fHistoMCEtaPtNoVertex[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCEtaPtNoVertex[iCut]);

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
        fHistoMCPi0InAccPtNotTriggered[iCut]           = new TH1F("MC_Pi0InAcc_Pt_NotTriggered", "MC_Pi0InAcc_Pt_NotTriggered", nBinsPt, arrPtBinning);
        fHistoMCPi0InAccPtNotTriggered[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0InAccPtNotTriggered[iCut]);
        fHistoMCPi0InAccPtNoVertex[iCut]           = new TH1F("MC_Pi0InAcc_Pt_NoVertex", "MC_Pi0InAcc_Pt_NoVertex", nBinsPt, arrPtBinning);
        fHistoMCPi0InAccPtNoVertex[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0InAccPtNoVertex[iCut]);
        fHistoMCEtaInAccPt[iCut]           = new TH1F("MC_EtaInAcc_Pt", "MC_EtaInAcc_Pt", nBinsPt, arrPtBinning);
        fHistoMCEtaInAccPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCEtaInAccPt[iCut]);
        fHistoMCEtaInAccPtNotTriggered[iCut]           = new TH1F("MC_EtaInAcc_Pt_NotTriggered", "MC_EtaInAcc_Pt_NotTriggered", nBinsPt, arrPtBinning);
        fHistoMCEtaInAccPtNotTriggered[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCEtaInAccPtNotTriggered[iCut]);
        fHistoMCEtaInAccPtNoVertex[iCut]           = new TH1F("MC_EtaInAcc_Pt_NoVertex", "MC_EtaInAcc_Pt_NoVertex", nBinsPt, arrPtBinning);
        fHistoMCEtaInAccPtNoVertex[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCEtaInAccPtNoVertex[iCut]);

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
    if (fEnableClusterCutsForTrigger){
      if(!((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))) continue;
      if(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutHistograms()){
        fCutFolder[iCut]->Add(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutHistograms());
      }
    }

  }
  // create tree to store broken files
  tBrokenFiles = new TTree("BrokenFiles", "BrokenFiles");
  tBrokenFiles->Branch("fileName",&fFileNameBroken);
  fOutputContainer->Add(tBrokenFiles);
  
  fAddressChanges = new TH1F("fAddressChanges","fAddressChanges", 2, 0., 2.);
  fOutputContainer->Add(fAddressChanges);

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
  fFileWasAlreadyReported = kFALSE;

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
  if(fIsMC>0) {
    fMCEvent = MCEvent();
    
    if(fInputEvent->IsA()==AliAODEvent::Class()){
      TClonesArray* lOldAddress = fAODMCTrackArray;
      fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
      // remove next line  and fAddressChanges after first trainrun 
      if (lOldAddress && fAODMCTrackArray!=lOldAddress) fAddressChanges->Fill(0.5);
      
      if (fAODMCTrackArray == NULL){
        fAddressChanges->Fill(1.5);
        AliInfo("AODMCTrackArray could not be loaded");
        return;
      }
    }
  }

  //calculating the weight for the centrality flattening
  for(Int_t iCut = 0; iCut<fnCuts; iCut++){
    if(fDoCentralityFlat > 0){
      fWeightCentrality[iCut] = 1.;
      fWeightCentrality[iCut] = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetWeightForCentralityFlattening(fInputEvent);
    }
  }

  Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();
  if(fInputEvent->IsIncompleteDAQ()==kTRUE || fV0Reader->GetErrorAODRelabeling()) eventQuality = 2;  // incomplete event or relabeling failed
  // Event Not Accepted due to MC event missing or because it is incomplete or  wrong trigger for V0ReaderV1 => skip broken events/files
  if(eventQuality == 2 || eventQuality == 3){
    // write out name of broken file if it has not been done already
    if(!fFileWasAlreadyReported){
      fFileNameBroken = new TObjString(Form("%s",((TString)fV0Reader->GetCurrentFileName()).Data()));
      AliInfo(Form("Adding %s to tree of broken files.", fFileNameBroken->GetString().Data()));
      tBrokenFiles->Fill();
      delete fFileNameBroken;
      fFileWasAlreadyReported = kTRUE;
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
    fiEventCut = dynamic_cast<AliConvEventCuts*>(fEventCutArray->At(iCut));
    fiPhotonCut = dynamic_cast<AliConversionPhotonCuts*>(fCutArray->At(fiCut));
    fiMesonCut  = dynamic_cast<AliConversionMesonCuts*>(fMesonCutArray->At(fiCut));
    
    // reset the event cuts fAODMCTrackArray to nullptr at the beginning of each event since in AliConvEventCuts it is only obtained from a new event if it's not a nullptr - hence only once at the beginning of a job. 
    if(fIsMC) fiEventCut->ResetAODMCTrackArray(); 
    
    Int_t eventNotAccepted = fiEventCut->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon,kFALSE);
    if( fIsMC == 2 ){
      Float_t xsection      = -1.;
      Float_t ntrials       = -1.;
      fiEventCut->GetXSectionAndNTrials(fMCEvent,xsection,ntrials,fInputEvent);
      if((xsection==-1.) || (ntrials==-1.)) AliFatal(Form("ERROR: GetXSectionAndNTrials returned invalid xsection/ntrials, periodName from V0Reader: '%s'",fV0Reader->GetPeriodName().Data()));
      fProfileJetJetXSection[iCut]->Fill(0.,xsection);
      fhJetJetNTrials[iCut]->Fill("#sum{NTrials}",ntrials);
    }

    if( fIsMC > 0 ){
      fWeightJetJetMC       = 1;
      Float_t pthard = -1;
      Bool_t isMCJet        = fiEventCut->IsJetJetMCEventAccepted( fMCEvent, fWeightJetJetMC , pthard, fInputEvent);
      if (fIsMC == 3){
        Double_t weightMult   = fiEventCut->GetWeightForMultiplicity(fV0Reader->GetNumberOfPrimaryTracks());
        fWeightJetJetMC       = fWeightJetJetMC*weightMult;
      }

      if( fIsMC == 1 ) fWeightJetJetMC = 1;
      if(!isMCJet){
        fHistoNEvents[iCut]->Fill(10,fWeightJetJetMC);
        if( fIsMC > 1 ) fHistoNEventsWOWeight[iCut]->Fill(10);
        continue;
      }
    }

    if(eventNotAccepted!= 0){
      // cout << "event rejected due to wrong trigger: " <<eventNotAccepted << endl;
      fHistoNEvents[iCut]->Fill(eventNotAccepted,fWeightJetJetMC); // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1
      if( fIsMC > 1 ) fHistoNEventsWOWeight[iCut]->Fill(eventNotAccepted);
      if(fDoCentralityFlat > 0) fHistoNEventsWeighted[iCut]->Fill(eventNotAccepted, fWeightCentrality[iCut]*fWeightJetJetMC);

      if(fIsMC > 0){
        if( eventNotAccepted == 3 && (eventQuality == 0 || eventQuality == 3 || eventQuality == 5)){ // wrong trigger selected. However, we still want to count the MC particles fot these events! If MC particles should be rejected in addition, use IsMCTriggerSelected function
          if(fInputEvent->IsA()==AliESDEvent::Class())
            ProcessMCParticles(1);
          if(fInputEvent->IsA()==AliAODEvent::Class())
            ProcessAODMCParticles(1);
        }
      }
      continue;
    }

    if(eventQuality != 0){// Event Not Accepted
      // cout << "event rejected due to: " <<eventQuality << endl;
      fHistoNEvents[iCut]->Fill(eventQuality,fWeightJetJetMC);
      if( fIsMC > 1 ) fHistoNEventsWOWeight[iCut]->Fill(eventQuality);
      if(fDoCentralityFlat > 0) fHistoNEventsWeighted[iCut]->Fill(eventQuality, fWeightCentrality[iCut]*fWeightJetJetMC);

      if(fIsMC > 0){
        if(eventQuality == 3 || eventQuality == 5){  // 3 = wrong trigger, 5 = GetNumberOfContributorsVtx
          if(fInputEvent->IsA()==AliESDEvent::Class())
            ProcessMCParticles(2);
          if(fInputEvent->IsA()==AliAODEvent::Class())
            ProcessAODMCParticles(2);
        }
      }
      continue;
    }

    if(fiEventCut->GetUseSphericity()!=0){
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

      if(fDoPlotVsCentrality) fHistoCentrality[iCut]->Fill(fiEventCut->GetCentrality(fInputEvent),fWeightJetJetMC);
      if(fDoCentralityFlat > 0) fHistoCentralityFlattened[iCut]->Fill(fiEventCut->GetCentrality(fInputEvent), fWeightCentrality[iCut]*fWeightJetJetMC);

      if(fDoCentralityFlat > 0) fHistoCentralityVsPrimaryTracks[iCut]->Fill(fiEventCut->GetCentrality(fInputEvent),fV0Reader->GetNumberOfPrimaryTracks(), fWeightCentrality[iCut]*fWeightJetJetMC);
      else if(fDoPlotVsCentrality) fHistoCentralityVsPrimaryTracks[iCut]->Fill(fiEventCut->GetCentrality(fInputEvent),fV0Reader->GetNumberOfPrimaryTracks(), fWeightJetJetMC);

      if(!fDoLightOutput){

        fiEventCut->FillTPCPileUpHistogram(fInputEvent);

        if( fIsMC < 2 ){
          if(fDoCentralityFlat > 0) fHistoSPDClusterTrackletBackground[iCut]->Fill(fInputEvent->GetMultiplicity()->GetNumberOfTracklets(),(fInputEvent->GetNumberOfITSClusters(0)+fInputEvent->GetNumberOfITSClusters(1)), fWeightCentrality[iCut]);
          else fHistoSPDClusterTrackletBackground[iCut]->Fill(fInputEvent->GetMultiplicity()->GetNumberOfTracklets(),(fInputEvent->GetNumberOfITSClusters(0)+fInputEvent->GetNumberOfITSClusters(1)));
        }

        if(fiEventCut->IsHeavyIon() == 2) fHistoNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A(),fWeightJetJetMC);
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
      if(fiEventCut->GetSignalRejection() != 0){
        if(fInputEvent->IsA()==AliESDEvent::Class()){
          fiEventCut->GetNotRejectedParticles(fiEventCut->GetSignalRejection(),
                                              fiEventCut->GetAcceptedHeader(),
                                              fMCEvent);
        }
        else if(fInputEvent->IsA()==AliAODEvent::Class()){
          if(fiEventCut->GetSignalRejection()==5) fiEventCut->ResetMcHeader();
          fiEventCut->GetNotRejectedParticles(fiEventCut->GetSignalRejection(),
                                              fiEventCut->GetAcceptedHeader(),
                                              fInputEvent);
        }

        /* todo: don't do this for each event. How are differing bin labels even merged?
         * todo: use the actually used headers  (that is the intersection of user specified
         *       cutselection and AliConvEventCuts::fHeaderList and what is contained in the events)*/
        if(fiEventCut->GetAcceptedHeader()){
          for(Int_t i = 0;i< fiEventCut->GetAcceptedHeader()->GetEntries();i++){
            if(fIsMC < 2){
              TString nameBin= fHistoMCHeaders[iCut]->GetXaxis()->GetBinLabel(i+1);
              if (nameBin.CompareTo("")== 0){
                TString nameHeader = ((TObjString*)((TList*)fiEventCut->GetAcceptedHeader())->At(i))->GetString();
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
    if(fDoHighPtHadronAnalysis) ProcessPhotonsHighPtHadronAnalysis();

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
        } else if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->BackgroundHandlerType() == 2){
          CalculateBackgroundSwapp();
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
  // ProcessPhotonCandidates() starts after definition of following lambda function
  auto fillHistosAndTree = [&](AliAODConversionPhoton *thePhoton){

    if( fIsMC > 0 ){
      if(fInputEvent->IsA()==AliESDEvent::Class()) ProcessTruePhotonCandidates(thePhoton);
      if(fInputEvent->IsA()==AliAODEvent::Class()) ProcessTruePhotonCandidatesAOD(thePhoton);
    }

    Double_t lPt = thePhoton->Pt();
    Float_t lWeightMatBudgetGamma = (fDoMaterialBudgetWeightingOfGammasForTrueMesons &&
                                     fiPhotonCut->GetMaterialBudgetWeightsInitialized())
        ? fiPhotonCut->GetMaterialBudgetCorrectingWeightForTrueGamma(thePhoton, fInputEvent->GetMagneticField())
        : 1.;

    Double_t lTotalWeight = (fDoCentralityFlat > 0) ? fWeightJetJetMC*lWeightMatBudgetGamma*fWeightCentrality[fiCut]
                                                    : fWeightJetJetMC*lWeightMatBudgetGamma;
    fHistoConvGammaPt[fiCut]->Fill(lPt, lTotalWeight);

    if (fDoPhotonQA > 0 && fIsMC < 2){
      fHistoConvGammaPsiPairPt[fiCut]->Fill(thePhoton->GetPsiPair(),lPt, lTotalWeight);
      fHistoConvGammaR        [fiCut]->Fill(thePhoton->GetConversionRadius(), lTotalWeight);
      fHistoConvGammaEta      [fiCut]->Fill(thePhoton->Eta(), lTotalWeight);
      fHistoConvGammaPhi      [fiCut]->Fill(thePhoton->Phi(), lTotalWeight);

      if (fDoPhotonQA == 4 || fDoPhotonQA == 5){
        fHistoConvGammaInvMass    [fiCut]->Fill(thePhoton->GetMass(), lTotalWeight);
        fHistoConvGammaInvMassReco[fiCut]->Fill(GetOriginalInvMass(thePhoton, fInputEvent), lTotalWeight);
      }
    }

    if (fDoPhotonQA == 2 || fDoPhotonQA == 5){
      if ((fIsHeavyIon == 1                             && lPt > 0.399 && lPt < 12.) ||
          (fiPhotonCut->GetSingleElectronPtCut() < 0.04 && lPt > 0.099 && lPt < 16.) ||
          (                                                lPt > 0.299 && lPt < 16.))
      {
        fPtGamma     = lPt;
        fDCAzPhoton  = thePhoton->GetDCAzToPrimVtx();
        fRConvPhoton = thePhoton->GetConversionRadius();
        fEtaPhoton   = thePhoton->GetPhotonEta();
        iCatPhoton   = thePhoton->GetPhotonQuality();
        tESDConvGammaPtDcazCat[fiCut]->Fill();
      }
    }
  }; // end of lambda fillHistosAndTree()

  // ProcessPhotonCandidates() starts here
  if(fiPhotonCut->GetDoElecDeDxPostCalibration()){
    if(!(fiPhotonCut->LoadElecDeDxPostCalibration(fInputEvent->GetRunNumber()))){
      AliFatal(Form("ERROR: LoadElecDeDxPostCalibration returned kFALSE for %d despite being requested!",fInputEvent->GetRunNumber()));
    }
  }

  fMapPhotonHeaders.clear();
  Bool_t lUseElecShareCut = fiPhotonCut->UseElecSharingCut();
  Bool_t lUseTooCloseCut  = fiPhotonCut->UseToCloseV0sCut();

  // Loop over Photon Candidates allocated by ReaderV1
  for (TObject *iObj : *fReaderGammas){

    AliAODConversionPhoton *iCandidate = dynamic_cast<AliAODConversionPhoton*>(iObj);
    if (!iCandidate) { AliWarning("Non AliAODConversionPhoton type object in fReaderGammas.\n"); continue; }

    Bool_t lIsFromSelectedHeader = kTRUE;
    if(fIsMC){
      if (!fiEventCut->PhotonPassesAddedParticlesCriterion(fMCEvent, fInputEvent, *iCandidate, lIsFromSelectedHeader)) continue;
    }

    if(!fiPhotonCut->PhotonIsSelected(iCandidate,fInputEvent)) continue;
    if(!fiPhotonCut->InPlaneOutOfPlaneCut(iCandidate->GetPhotonPhi(),fEventPlaneAngle)) continue;

    // if no further cuts, add to fGammaCandidates and we are done. If header criterion is fullfilled, also fill histos and tree
    if (!(lUseElecShareCut || lUseTooCloseCut)){

      fGammaCandidates->Add(iCandidate);
      if (lIsFromSelectedHeader){
        fillHistosAndTree(iCandidate);
      }
    }
    else{
      // we have one of lUseElecShareCut and lUseTooCloseCut -> we cant fill the histos before having looked at all photons
      fMapPhotonHeaders.insert({iCandidate, lIsFromSelectedHeader});
    }
  }

  if (lUseElecShareCut) fiPhotonCut->RemovePhotonsWithSharedTracks(fMapPhotonHeaders);
  if (lUseTooCloseCut) fiPhotonCut->RemoveTooClosePhotons(fMapPhotonHeaders);

  // add remaining candidates to fGammaCandidates. If header criterion is fullfilled, also fill histos and tree
  // note: fMapPhotonHeaders will be empty unless lUseElecShareCut || lUseTooCloseCut
  for (auto &iPhotonHeader : fMapPhotonHeaders){

    fGammaCandidates->Add(iPhotonHeader.first);
    if (iPhotonHeader.second){
      fillHistosAndTree(iPhotonHeader.first);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::ProcessPhotonsHighPtHadronAnalysis()
{
  Bool_t DoesEventContainHighPtHadron = kFALSE;
  Double_t NTotalTracks = fInputEvent->GetNumberOfTracks();
  Int_t NTracks = 0;
  for(Int_t iTracks = 0; iTracks<NTotalTracks; iTracks++){
    AliAODTrack* curTrack = (AliAODTrack*) fInputEvent->GetTrack(iTracks);
    if(curTrack->GetID()<0) continue; // Avoid double counting of tracks
    if(!curTrack->IsHybridGlobalConstrainedGlobal()) continue;
    if(TMath::Abs(curTrack->Eta())>0.8) continue;
    if(curTrack->Pt()<0.50) continue;
    if(curTrack->Pt()>10) DoesEventContainHighPtHadron = kTRUE;
    NTracks++;
  }

  if(fGammaCandidates->GetEntries()>1){
    for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries()-1;firstGammaIndex++){
      AliAODConversionPhoton *gamma=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));
      if(gamma==NULL) continue;
      if(DoesEventContainHighPtHadron){
        fHistoConvGammaPtwithHighPtHadron[fiCut]->Fill(gamma->Pt(),NTracks);
      } else {
        fHistoConvGammaPtwithoutHighPtHadron[fiCut]->Fill(gamma->Pt(),NTracks);
      }
    }
  }

  if(DoesEventContainHighPtHadron){
        fHistoNEventsHighPtHadron[fiCut]->Fill(0);
  } else {
        fHistoNEventsHighPtHadron[fiCut]->Fill(1);
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
    if(fiPhotonCut->GetBDTVariableValues(PhotonCandidate,fInputEvent, SingleTrackBDTVariableValues)){
      fBDTvariable[0] = SingleTrackBDTVariableValues[0]; //"dEdxElectronITS + dEdxPositronITS"
      fBDTvariable[1] = PhotonCandidate->GetArmenterosQt(); //"photonQt"
      fBDTvariable[2] = GetOriginalInvMass(PhotonCandidate,fInputEvent); //"photonInvMass"
      fBDTvariable[3] = PhotonCandidate->GetConversionRadius(); //"photonR"
      fBDTvariable[4] = PhotonCandidate->GetArmenterosAlpha(); //"photonAlpha"
      fBDTvariable[5] = PhotonCandidate->GetPsiPair(); //"photonPsiPair"
      fBDTvariable[6] = fiPhotonCut->GetCosineOfPointingAngle(PhotonCandidate,fInputEvent); //"photonCosPoint"
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
	//        TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
        if (fAODMCTrackArray != NULL && PhotonCandidate != NULL){
          AliAODMCParticle *posDaughter = (AliAODMCParticle*) fAODMCTrackArray->At(PhotonCandidate->GetMCLabelPositive());
          AliAODMCParticle *negDaughter = (AliAODMCParticle*) fAODMCTrackArray->At(PhotonCandidate->GetMCLabelNegative());
          if(posDaughter != NULL && negDaughter != NULL){
              Int_t pdgCode[2] = {TMath::Abs(posDaughter->GetPdgCode()),TMath::Abs(negDaughter->GetPdgCode())};
              if(posDaughter->GetMother() == negDaughter->GetMother()){
                if((pdgCode[0]==11 && pdgCode[1]==11)&&  posDaughter->GetPdgCode()!=negDaughter->GetPdgCode()){
                  AliAODMCParticle *Photon = (AliAODMCParticle*) fAODMCTrackArray->At(posDaughter->GetMother());
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
  Double_t magFieldFlip = 1.0;
  if( magField  < 0.0 ){
    magFieldFlip =  1.0;
  }
  else {
    magFieldFlip =  -1.0;
  }

  if (fAODMCTrackArray != NULL && TruePhotonCandidate != NULL){

    AliAODMCParticle *posDaughter = (AliAODMCParticle*) fAODMCTrackArray->At(TruePhotonCandidate->GetMCLabelPositive());
    AliAODMCParticle *negDaughter = (AliAODMCParticle*) fAODMCTrackArray->At(TruePhotonCandidate->GetMCLabelNegative());
    iPhotonMCInfo = 0;

    if(posDaughter == NULL || negDaughter == NULL) return; // One particle does not exist
    Int_t pdgCode[2] = {TMath::Abs(posDaughter->GetPdgCode()),TMath::Abs(negDaughter->GetPdgCode())};

    Double_t PhiParticle[2] = {posDaughter->Phi(),negDaughter->Phi()};

    if(posDaughter->GetMother() != negDaughter->GetMother()){
      FillPhotonCombinatorialBackgroundHist(TruePhotonCandidate, pdgCode, PhiParticle);
      if(fDoPhotonQA == 3 && fIsHeavyIon == 1){
          if(posDaughter->GetMother() > -1){ //contamination is not a primary
              AliAODMCParticle *Mom = (AliAODMCParticle*) fAODMCTrackArray->At(posDaughter->GetMother());
              if(Mom->GetMother() == -1){
                  FillPhotonCombinatorialMothersHistAOD(posDaughter,Mom);
              } else {
                  AliAODMCParticle *GranMom = (AliAODMCParticle*) fAODMCTrackArray->At(Mom->GetMother());
                  if(GranMom->GetMother() == -1){
                      FillPhotonCombinatorialMothersHistAOD(posDaughter,GranMom);
                  } else {
                      AliAODMCParticle *GranGranMom = (AliAODMCParticle*) fAODMCTrackArray->At(GranMom->GetMother());
                      if(GranGranMom->GetMother() == -1){
                          FillPhotonCombinatorialMothersHistAOD(posDaughter,GranGranMom);
                      } else FillPhotonCombinatorialMothersHistAOD(posDaughter,GranGranMom);
                  }
              }
          }
          if(negDaughter->GetMother() > -1){ //contamination is not a primary
              AliAODMCParticle *Mom = (AliAODMCParticle*) fAODMCTrackArray->At(negDaughter->GetMother());
              if(Mom->GetMother() == -1){
                  FillPhotonCombinatorialMothersHistAOD(negDaughter,Mom);
              } else {
                  AliAODMCParticle *GranMom = (AliAODMCParticle*) fAODMCTrackArray->At(Mom->GetMother());
                  if(GranMom->GetMother() == -1){
                      FillPhotonCombinatorialMothersHistAOD(negDaughter,GranMom);
                  } else {
                      AliAODMCParticle *GranGranMom = (AliAODMCParticle*) fAODMCTrackArray->At(GranMom->GetMother());
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

    AliAODMCParticle *Photon = (AliAODMCParticle*) fAODMCTrackArray->At(posDaughter->GetMother());
    AliVTrack * electronCandidate = fiPhotonCut->GetTrack(fInputEvent,TruePhotonCandidate->GetTrackLabelNegative() );
    AliVTrack * positronCandidate = fiPhotonCut->GetTrack(fInputEvent,TruePhotonCandidate->GetTrackLabelPositive() );
    Double_t deltaPhi = magFieldFlip * TVector2::Phi_mpi_pi( electronCandidate->Phi()-positronCandidate->Phi());

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
    if (fDoMaterialBudgetWeightingOfGammasForTrueMesons && fiPhotonCut->GetMaterialBudgetWeightsInitialized()) {
      weightMatBudgetGamma = fiPhotonCut->GetMaterialBudgetCorrectingWeightForTrueGamma(TruePhotonCandidate,magField);
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

    Bool_t isPrimary = fiEventCut->IsConversionPrimaryAOD(fInputEvent, Photon, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
    if(isPrimary){
      // Count just primary MC Gammas as true --> For Ratio esdtruegamma / mcconvgamma
      iPhotonMCInfo = 6;
      fHistoTruePrimaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
      fHistoTruePrimaryConvGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt(),fWeightJetJetMC*weightMatBudgetGamma); // Allways Filled
      // (Not Filled for i6, Extra Signal Gamma (parambox) are secondary)
    } else {
      iPhotonMCInfo = 2;
      if(((AliAODMCParticle*)fAODMCTrackArray->At(Photon->GetMother()))->GetMother() > -1 && ((AliAODMCParticle*)fAODMCTrackArray->At(((AliAODMCParticle*)fAODMCTrackArray->At(Photon->GetMother()))->GetMother()))->GetMother() > -1 ){
        if (((AliAODMCParticle*)fAODMCTrackArray->At(((AliAODMCParticle*)fAODMCTrackArray->At(Photon->GetMother()))->GetMother()))->GetPdgCode() == 310){
            iPhotonMCInfo = 4;
            fHistoTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),0.,fWeightJetJetMC*weightMatBudgetGamma);
            fHistoTrueSecondaryConvGammaMCPt[fiCut]->Fill(Photon->Pt(),0.,fWeightJetJetMC*weightMatBudgetGamma);
            fHistoTrueSecondaryConvGammaFromXFromK0sMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
        } else if (((AliAODMCParticle*)fAODMCTrackArray->At(((AliAODMCParticle*)fAODMCTrackArray->At(Photon->GetMother()))->GetMother()))->GetPdgCode() == 130) {
            iPhotonMCInfo = 7;
            fHistoTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),1.,fWeightJetJetMC*weightMatBudgetGamma);
            fHistoTrueSecondaryConvGammaMCPt[fiCut]->Fill(Photon->Pt(),1.,fWeightJetJetMC*weightMatBudgetGamma);
            fHistoTrueSecondaryConvGammaFromXFromK0lMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
        } else if (((AliAODMCParticle*)fAODMCTrackArray->At(((AliAODMCParticle*)fAODMCTrackArray->At(Photon->GetMother()))->GetMother()))->GetPdgCode() == 3122) {
            iPhotonMCInfo = 5;
            fHistoTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),2.,fWeightJetJetMC*weightMatBudgetGamma);
            fHistoTrueSecondaryConvGammaMCPt[fiCut]->Fill(Photon->Pt(),2.,fWeightJetJetMC*weightMatBudgetGamma);
            fHistoTrueSecondaryConvGammaFromXFromLambdaMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
        } else if (((AliAODMCParticle*)fAODMCTrackArray->At(((AliAODMCParticle*)fAODMCTrackArray->At(Photon->GetMother()))->GetMother()))->GetPdgCode() == 221) {
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
  Double_t magFieldFlip = 1.0;
  if( magField  < 0.0 ){
    magFieldFlip =  1.0;
  }
  else {
    magFieldFlip =  -1.0;
  }

  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  // Process True Photons
  AliMCParticle *posDaughter = (AliMCParticle*) TruePhotonCandidate->GetPositiveMCDaughter(fMCEvent);
  AliMCParticle *negDaughter = (AliMCParticle*) TruePhotonCandidate->GetNegativeMCDaughter(fMCEvent);

  iPhotonMCInfo = 0;

  if(posDaughter == NULL || negDaughter == NULL) return; // One particle does not exist
  Int_t pdgCode[2] = {TMath::Abs(posDaughter->PdgCode()),TMath::Abs(negDaughter->PdgCode())};

  Double_t PhiParticle[2] = {posDaughter->Phi(),negDaughter->Phi()};

  iPhotonMCInfo = 1;
  if(posDaughter->GetMother() != negDaughter->GetMother()){
    FillPhotonCombinatorialBackgroundHist(TruePhotonCandidate, pdgCode, PhiParticle);
    if(fDoPhotonQA == 3 && fIsHeavyIon == 1){
        if(posDaughter->GetMother() > -1){ //contamination is not a primary
            AliMCParticle *Mom = (AliMCParticle*) fMCEvent->GetTrack(posDaughter->GetMother());
            if(Mom->GetMother() == -1){
                FillPhotonCombinatorialMothersHistESD(posDaughter,Mom);
            } else {
                AliMCParticle *GranMom = (AliMCParticle*) fMCEvent->GetTrack(Mom->GetMother());
                if(GranMom->GetMother() == -1){
                    FillPhotonCombinatorialMothersHistESD(posDaughter,GranMom);
                } else {
                    AliMCParticle *GranGranMom = (AliMCParticle*) fMCEvent->GetTrack(GranMom->GetMother());
                    if(GranGranMom->GetMother() == -1){
                        FillPhotonCombinatorialMothersHistESD(posDaughter,GranGranMom);
                    } else FillPhotonCombinatorialMothersHistESD(posDaughter,GranGranMom);
                }
            }
        }
        if(negDaughter->GetMother() > -1){ //contamination is not a primary
            AliMCParticle *Mom = (AliMCParticle*) fMCEvent->GetTrack(negDaughter->GetMother());
            if(Mom->GetMother() == -1){
                FillPhotonCombinatorialMothersHistESD(negDaughter,Mom);
            } else {
                AliMCParticle *GranMom = (AliMCParticle*) fMCEvent->GetTrack(Mom->GetMother());
                if(GranMom->GetMother() == -1){
                    FillPhotonCombinatorialMothersHistESD(negDaughter,GranMom);
                } else {
                    AliMCParticle *GranGranMom = (AliMCParticle*) fMCEvent->GetTrack(GranMom->GetMother());
                    if(GranGranMom->GetMother() == -1){
                        FillPhotonCombinatorialMothersHistESD(negDaughter,GranGranMom);
                    } else FillPhotonCombinatorialMothersHistESD(posDaughter,GranGranMom);
                }
            }
        }
    }
    return;
  } else if(posDaughter->GetMother() == -1){ //gamma contamination is a primary
    FillPhotonCombinatorialBackgroundHist(TruePhotonCandidate, pdgCode, PhiParticle);
    if(fDoPhotonQA == 3 && fIsHeavyIon == 1){
        FillPhotonCombinatorialMothersHistESD(posDaughter,posDaughter);
        FillPhotonCombinatorialMothersHistESD(negDaughter,negDaughter);
    }
    return;
  }

  if(pdgCode[0]!=11 || pdgCode[1]!=11) return; //One Particle is not a electron

  if(posDaughter->PdgCode()==negDaughter->PdgCode()) return; // Same Charge

  AliVParticle *Photon = TruePhotonCandidate->GetMCParticle(fMCEvent);
  AliVTrack * electronCandidate = fiPhotonCut->GetTrack(fInputEvent,TruePhotonCandidate->GetTrackLabelNegative() );
  AliVTrack * positronCandidate = fiPhotonCut->GetTrack(fInputEvent,TruePhotonCandidate->GetTrackLabelPositive() );
  Double_t deltaPhi = magFieldFlip * TVector2::Phi_mpi_pi( electronCandidate->Phi()-positronCandidate->Phi());

  if(Photon->PdgCode() != 22){
    fHistoTrueDalitzPsiPairDeltaPhi[fiCut]->Fill(deltaPhi,TruePhotonCandidate->GetPsiPair(),fWeightJetJetMC);
    return; // Mother is no Photon
  }

  if(posDaughter->Particle()->GetUniqueID() != 5 || negDaughter->Particle()->GetUniqueID() !=5) return;// check if the daughters come from a conversion

  // True Photon
  Float_t weightMatBudgetGamma = 1.;
  if (fDoMaterialBudgetWeightingOfGammasForTrueMesons && fiPhotonCut->GetMaterialBudgetWeightsInitialized()) {
    weightMatBudgetGamma = fiPhotonCut->GetMaterialBudgetCorrectingWeightForTrueGamma(TruePhotonCandidate,magField);
  }

  // cout<< " AM- Material Budget weight Gamma::"<< weightMatBudgetGamma << " "<< TruePhotonCandidate->GetConversionRadius() << endl;

  fHistoTrueConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
  if (CheckVectorForDoubleCount(vecDoubleCountTrueConvGammas,posDaughter->GetMother())){
    fHistoDoubleCountTrueConvGammaRPt[fiCut]->Fill(TruePhotonCandidate->GetConversionRadius(),TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
    FillMultipleCountMap(mapMultipleCountTrueConvGammas,posDaughter->GetMother());
  }
  if (fDoPhotonQA > 0){
    if (fIsMC < 2){
      fHistoTrueConvGammaPsiPairPt[fiCut]->Fill(TruePhotonCandidate->GetPsiPair(),TruePhotonCandidate->Pt(),weightMatBudgetGamma);
      fHistoTrueConvGammaEta[fiCut]->Fill(TruePhotonCandidate->Eta(),weightMatBudgetGamma);
      fHistoTrueConvGammaR[fiCut]->Fill(TruePhotonCandidate->GetConversionRadius(),weightMatBudgetGamma);
      fHistoTrueConvGammaRMC[fiCut]->Fill(posDaughter->Particle()->R(),weightMatBudgetGamma);
      if ((fDoPhotonQA == 4)||(fDoPhotonQA == 5)){
        fHistoTrueConvGammaInvMass[fiCut]->Fill(TruePhotonCandidate->GetMass(),weightMatBudgetGamma);
        fHistoTrueConvGammaInvMassReco[fiCut]->Fill(GetOriginalInvMass(TruePhotonCandidate,fInputEvent),weightMatBudgetGamma);
      }
    }
    fHistoTrueConvGammaPtMC[fiCut]->Fill(Photon->Pt(), fWeightJetJetMC*weightMatBudgetGamma);
  }

  fHistoTrueGammaPsiPairDeltaPhi[fiCut]->Fill(deltaPhi,TruePhotonCandidate->GetPsiPair(),fWeightJetJetMC*weightMatBudgetGamma);
  if (fiEventCut->IsConversionPrimaryESD( fMCEvent, posDaughter->GetMother(), mcProdVtxX, mcProdVtxY, mcProdVtxZ)){
    // filling primary histograms
    // Count just primary MC Gammas as true --> For Ratio esdtruegamma / mcconvgamma
    iPhotonMCInfo = 6;
    fHistoTruePrimaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
    fHistoTruePrimaryConvGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt(),fWeightJetJetMC*weightMatBudgetGamma); // Allways Filled
    // (Not Filled for i6, Extra Signal Gamma (parambox) are secondary)
  } else {
    // filling secondary photon histograms
    iPhotonMCInfo = 2;
    if( Photon->GetMother() > -1 && fMCEvent->GetTrack(Photon->GetMother())->GetMother() > -1){
      if (fMCEvent->GetTrack(fMCEvent->GetTrack(Photon->GetMother())->GetMother())->PdgCode() == 310){
        iPhotonMCInfo = 4;
        fHistoTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),0.,fWeightJetJetMC*weightMatBudgetGamma);
        fHistoTrueSecondaryConvGammaMCPt[fiCut]->Fill(Photon->Pt(),0.,fWeightJetJetMC*weightMatBudgetGamma);
        fHistoTrueSecondaryConvGammaFromXFromK0sMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
      } else if (fMCEvent->GetTrack(fMCEvent->GetTrack(Photon->GetMother())->GetMother())->PdgCode() == 130) {
        iPhotonMCInfo = 7;
        fHistoTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),1.,fWeightJetJetMC*weightMatBudgetGamma);
        fHistoTrueSecondaryConvGammaMCPt[fiCut]->Fill(Photon->Pt(),1.,fWeightJetJetMC*weightMatBudgetGamma);
        fHistoTrueSecondaryConvGammaFromXFromK0lMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
      } else if (fMCEvent->GetTrack(fMCEvent->GetTrack(Photon->GetMother())->GetMother())->PdgCode() == 3122) {
        iPhotonMCInfo = 5;
        fHistoTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),2.,fWeightJetJetMC*weightMatBudgetGamma);
        fHistoTrueSecondaryConvGammaMCPt[fiCut]->Fill(Photon->Pt(),2.,fWeightJetJetMC*weightMatBudgetGamma);
        fHistoTrueSecondaryConvGammaFromXFromLambdaMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),fWeightJetJetMC*weightMatBudgetGamma);
      } else if (fMCEvent->GetTrack(fMCEvent->GetTrack(Photon->GetMother())->GetMother())->PdgCode() == 221) {
        iPhotonMCInfo = 3;
        fHistoTrueSecondaryConvGammaPt[fiCut]->Fill(Photon->Pt(),3.,fWeightJetJetMC*weightMatBudgetGamma);
        fHistoTrueSecondaryConvGammaMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3.,fWeightJetJetMC*weightMatBudgetGamma);
      } else {
//         if ( !(TMath::Abs(fMCEvent->GetTrack(Photon->GetMother())->PdgCode()) == 11 && fMCEvent->GetTrack(fMCEvent->GetTrack(Photon->GetMother())->GetMother())->PdgCode() == 22) ) {
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
void AliAnalysisTaskGammaConvV1::ProcessAODMCParticles(int isCurrentEventSelected)
{
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  // Check if MC generated particles should be filled for this event using the selected trigger
  if( !((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsMCTriggerSelected(fInputEvent, fMCEvent)){
    return;
  }

  if (fAODMCTrackArray){
    // Loop over all primary MC particle
    for(Int_t i = 0; i < fAODMCTrackArray->GetEntriesFast(); i++) {

      AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(i));
      if (!particle) continue;

      Bool_t isPrimary = fiEventCut->IsConversionPrimaryAOD(fInputEvent, particle, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
      if (isPrimary){

        Int_t isMCFromMBHeader = -1;
        if(fiEventCut->GetSignalRejection() != 0){
          isMCFromMBHeader = fiEventCut->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
          if(isMCFromMBHeader == 0 && fiEventCut->GetSignalRejection() != 3) continue;
        }

        if(!fiPhotonCut->InPlaneOutOfPlaneCut(particle->Phi(),fEventPlaneAngle,kFALSE)) continue;
        if(fiPhotonCut->PhotonIsSelectedAODMC(particle,fAODMCTrackArray,kFALSE)){
          fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // All MC Gamma
          if(isCurrentEventSelected == 1){
            fHistoMCAllGammaPtNotTriggered[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
          } else if(isCurrentEventSelected == 2){
            fHistoMCAllGammaPtNoVertex[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
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
        if(fiPhotonCut->PhotonIsSelectedAODMC(particle,fAODMCTrackArray,kTRUE)){
          Double_t rConv = 0;
          for(Int_t daughterIndex=particle->GetDaughterLabel(0);daughterIndex<=particle->GetDaughterLabel(1);daughterIndex++){
            AliAODMCParticle *tmpDaughter = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(daughterIndex));
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
            mesonY = particle->Y()-fiEventCut->GetEtaShift();
          }

          if ((mesonY > fiMesonCut->GetRapidityCutValueMin()) && (mesonY < fiMesonCut->GetRapidityCutValueMax())){
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

          if(fiMesonCut->MesonIsSelectedAODMC(particle,fAODMCTrackArray,fiEventCut->GetEtaShift())){
            AliAODMCParticle* daughter0 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetDaughterLabel(0)));
            AliAODMCParticle* daughter1 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetDaughterLabel(1)));
            Float_t weighted= 1;
	        if (particle->Pt()>0.005){
	          weighted= fiEventCut->GetWeightForMeson(i, 0x0, fInputEvent);
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
              mesonY = particle->Y()-fiEventCut->GetEtaShift();
            }

            Double_t alpha = -1;
            if (particle->GetPdgCode() == 111 || particle->GetPdgCode() == 221){
              alpha = TMath::Abs((daughter0->E() - daughter1->E()))/(daughter0->E() + daughter1->E());
            }

            if(particle->GetPdgCode() == 111){
              fHistoMCPi0Pt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // All MC Pi0
              if(isCurrentEventSelected == 1){
                fHistoMCPi0PtNotTriggered[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
              } else if(isCurrentEventSelected == 2){
                fHistoMCPi0PtNoVertex[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
              }
              fHistoMCPi0WOWeightPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
              if ( fIsMC > 1 ) fHistoMCPi0WOEvtWeightPt[fiCut]->Fill(particle->Pt());
              if (fDoMesonQA > 0){
                if ( fIsMC < 2 )fHistoMCPi0PtY[fiCut]->Fill(particle->Pt(),mesonY,weighted); // All MC Pi0
                fHistoMCPi0PtAlpha[fiCut]->Fill(particle->Pt(),alpha,fWeightJetJetMC); // All MC Pi0
                if ( fIsMC == 2 ) fHistoMCPi0PtJetPt[fiCut]->Fill(particle->Pt(),fiEventCut->GetMaxPtJet(),fWeightJetJetMC);
              }
            } else if(particle->GetPdgCode() == 221){
              fHistoMCEtaPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // All MC Eta
              if(isCurrentEventSelected == 1){
                fHistoMCEtaPtNotTriggered[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
              } else if(isCurrentEventSelected == 2){
                fHistoMCEtaPtNoVertex[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
              }
            
              fHistoMCEtaWOWeightPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
              if ( fIsMC > 1 ) fHistoMCEtaWOEvtWeightPt[fiCut]->Fill(particle->Pt());
              if (fDoMesonQA > 0){
                if ( fIsMC < 2 )fHistoMCEtaPtY[fiCut]->Fill(particle->Pt(),mesonY,weighted); // All MC Pi0
                fHistoMCEtaPtAlpha[fiCut]->Fill(particle->Pt(),alpha,fWeightJetJetMC); // All MC Pi0
                if ( fIsMC == 2 ) fHistoMCEtaPtJetPt[fiCut]->Fill(particle->Pt(),fiEventCut->GetMaxPtJet(),fWeightJetJetMC);
              }
            }

            // Check the acceptance for both gammas
            if(fiPhotonCut->PhotonIsSelectedAODMC(daughter0,fAODMCTrackArray,kFALSE) &&
               fiPhotonCut->PhotonIsSelectedAODMC(daughter1,fAODMCTrackArray,kFALSE)  &&
               fiPhotonCut->InPlaneOutOfPlaneCut(daughter0->Phi(),fEventPlaneAngle,kFALSE) &&
               fiPhotonCut->InPlaneOutOfPlaneCut(daughter1->Phi(),fEventPlaneAngle,kFALSE)){

              if(particle->GetPdgCode() == 111){
                if(fIsMC > 1) fHistoMCPi0WOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Pi0 with gamma in acc NOT weighted at all
                fHistoMCPi0WOWeightInAccPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // MC Pi0 with gamma in acc NOT weighted
                fHistoMCPi0InAccPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // MC Pi0 with gamma in acc
                if(isCurrentEventSelected == 1){
                  fHistoMCPi0InAccPtNotTriggered[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
                } else if(isCurrentEventSelected == 2){
                  fHistoMCPi0InAccPtNoVertex[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
                }
              } else if(particle->GetPdgCode() == 221){
                if(fIsMC > 1) fHistoMCEtaWOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Eta with gamma in acc NOT weighted at all
                fHistoMCEtaWOWeightInAccPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // MC Eta with gamma in acc NOT weighted
                fHistoMCEtaInAccPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // MC Eta with gamma in acc
                if(isCurrentEventSelected == 1){
                  fHistoMCEtaInAccPtNotTriggered[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
                } else if(isCurrentEventSelected == 2){
                  fHistoMCEtaInAccPtNoVertex[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
                }
              }
            }
          }
        }
      // fill secondary histograms
      } else {
        Int_t isMCFromMBHeader = -1;
        if(fiEventCut->GetSignalRejection() != 0){
          isMCFromMBHeader = fiEventCut->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
          if(isMCFromMBHeader == 0 && fiEventCut->GetSignalRejection() != 3) continue;
        }

        if(fiPhotonCut->InPlaneOutOfPlaneCut(particle->Phi(),fEventPlaneAngle,kFALSE)) {
          if(fiPhotonCut->PhotonIsSelectedAODMC(particle,fAODMCTrackArray,kFALSE)){
            if (particle->GetMother() > -1) {
              AliAODMCParticle *tmpMother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetMother()));
              if (tmpMother->GetMother() > -1) {
                AliAODMCParticle *tmpGrandMother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(tmpMother->GetMother()));
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

          if(fiPhotonCut->PhotonIsSelectedAODMC(particle,fAODMCTrackArray,kTRUE)){
            if (particle->GetMother() > -1) {
              AliAODMCParticle *tmpMother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetMother()));
              if (tmpMother->GetMother() > -1) {
                AliAODMCParticle *tmpGrandMother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(tmpMother->GetMother()));
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
          if(fiMesonCut->MesonIsSelectedAODMC(particle,fAODMCTrackArray,fiEventCut->GetEtaShift())){
            AliAODMCParticle* daughter0 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetDaughterLabel(0)));
            AliAODMCParticle* daughter1 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetDaughterLabel(1)));
            AliAODMCParticle* mother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetMother()));
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
            if( fiPhotonCut->PhotonIsSelectedAODMC(daughter0,fAODMCTrackArray,kFALSE) &&
                fiPhotonCut->PhotonIsSelectedAODMC(daughter1,fAODMCTrackArray,kFALSE)  &&
                fiPhotonCut->InPlaneOutOfPlaneCut(daughter0->Phi(),fEventPlaneAngle,kFALSE) &&
                fiPhotonCut->InPlaneOutOfPlaneCut(daughter1->Phi(),fEventPlaneAngle,kFALSE)){
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
void AliAnalysisTaskGammaConvV1::ProcessMCParticles(int isCurrentEventSelected)
{
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();
//   cout << mcProdVtxX <<"\t" << mcProdVtxY << "\t" << mcProdVtxZ << endl;

  // Check if MC generated particles should be filled for this event using the selected trigger
  if( !((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsMCTriggerSelected(fInputEvent, fMCEvent)){
    return;
  }

  // Loop over all primary MC particle
  for(Long_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {

    if (fiEventCut->IsConversionPrimaryESD( fMCEvent, i, mcProdVtxX, mcProdVtxY, mcProdVtxZ)){
      // fill primary histogram
      AliMCParticle* particle = (AliMCParticle *)fMCEvent->GetTrack(i);
      if (!particle) continue;

      Int_t isMCFromMBHeader = -1;
      if(fiEventCut->GetSignalRejection() != 0){
        isMCFromMBHeader = fiEventCut->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
        if(isMCFromMBHeader == 0 && fiEventCut->GetSignalRejection() != 3) continue;
      }

      if(!fiPhotonCut->InPlaneOutOfPlaneCut(particle->Phi(),fEventPlaneAngle,kFALSE)) continue;
      if(fiPhotonCut->PhotonIsSelectedMC(particle,fMCEvent,kFALSE)){
        fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // All MC Gamma
        if(isCurrentEventSelected == 1){
          fHistoMCAllGammaPtNotTriggered[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
        } else if(isCurrentEventSelected == 2){
          fHistoMCAllGammaPtNoVertex[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
        }

        if(particle->GetMother() >-1){ // Meson Decay Gamma
          switch(fMCEvent->GetTrack(particle->GetMother())->PdgCode()){
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
      if(fiPhotonCut->PhotonIsSelectedMC(particle,fMCEvent,kTRUE)){
        fHistoMCConvGammaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
        if (fDoPhotonQA > 0 && fIsMC < 2){
          fHistoMCConvGammaR[fiCut]->Fill(((AliMCParticle*)fMCEvent->GetTrack(particle->GetDaughterFirst()))->Particle()->R(),fWeightJetJetMC);
          fHistoMCConvGammaEta[fiCut]->Fill(particle->Eta(),fWeightJetJetMC);
        }
      } // Converted MC Gamma
      if(fDoMesonAnalysis){

          Double_t mesonY = 1.e30;
          Double_t ratio  = 0;
          if (particle->E() != TMath::Abs(particle->Pz())){
            ratio         = (particle->E()+particle->Pz()) / (particle->E()-particle->Pz());
          }
          if( !(ratio <= 0) ){
            mesonY = particle->Y()-fiEventCut->GetEtaShift();
          }

          if ((mesonY > fiMesonCut->GetRapidityCutValueMin()) && (mesonY < fiMesonCut->GetRapidityCutValueMax())){
            if ( particle->PdgCode() == 211 ){  // positve pions
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(),0.,fWeightJetJetMC);
            } else if ( particle->PdgCode() == -211 ){  // negative pions
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(),1.,fWeightJetJetMC);
            } else if ( particle->PdgCode() == 321 ){  // positve kaons
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(),2.,fWeightJetJetMC);
            } else if ( particle->PdgCode() == -321 ){  // negative kaons
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(),3.,fWeightJetJetMC);
            } else if ( TMath::Abs(particle->PdgCode()) == 310 ){  // K0s
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(),4.,fWeightJetJetMC);
            } else if ( TMath::Abs(particle->PdgCode()) == 130 ){  // K0l
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(),5.,fWeightJetJetMC);
            } else if ( TMath::Abs(particle->PdgCode()) == 3122 ){  // Lambda/ AntiLambda
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(),6.,fWeightJetJetMC);
            } else if ( TMath::Abs(particle->PdgCode()) == 223 ){  // Omega
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(),7.,fWeightJetJetMC);
            } else if ( TMath::Abs(particle->PdgCode()) == 333 ){  // Phi
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(),8.,fWeightJetJetMC);
            } else if ( TMath::Abs(particle->PdgCode()) == 113 ){  // Rho0
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(),9.,fWeightJetJetMC);
            }
          }

        if(fiMesonCut->MesonIsSelectedMC(particle,fMCEvent,fiEventCut->GetEtaShift())){
          AliMCParticle* daughter0 = (AliMCParticle*)fMCEvent->GetTrack(particle->GetDaughterFirst());
          AliMCParticle* daughter1 = (AliMCParticle*)fMCEvent->GetTrack(particle->GetDaughterLast());

          Float_t weighted= 1;
	      if (particle->Pt()>0.005){
	        weighted= fiEventCut->GetWeightForMeson(i, fMCEvent, fInputEvent);
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
            mesonY = particle->Y()-fiEventCut->GetEtaShift();
          }

          Double_t alpha = -1;
          if (particle->PdgCode() == 111 || particle->PdgCode() == 221){
            alpha = TMath::Abs((daughter0->E() - daughter1->E()))/(daughter0->E() + daughter1->E());
          }

          if(particle->PdgCode() == 111){
            fHistoMCPi0Pt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // All MC Pi0
            if(isCurrentEventSelected == 1){
              fHistoMCPi0PtNotTriggered[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
            } else if(isCurrentEventSelected == 2){
              fHistoMCPi0PtNoVertex[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
            }
            fHistoMCPi0WOWeightPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
            if (fIsMC > 1) fHistoMCPi0WOEvtWeightPt[fiCut]->Fill(particle->Pt());
            if (fDoMesonQA > 0){
              if (fIsMC < 2)fHistoMCPi0PtY[fiCut]->Fill(particle->Pt(),mesonY,weighted); // All MC Pi0
              fHistoMCPi0PtAlpha[fiCut]->Fill(particle->Pt(),alpha,fWeightJetJetMC); // All MC Pi0
              if (fIsMC == 2) fHistoMCPi0PtJetPt[fiCut]->Fill(particle->Pt(),fiEventCut->GetMaxPtJet(),fWeightJetJetMC);
            }
          } else if(particle->PdgCode() == 221){
            fHistoMCEtaPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // All MC Eta
            if(isCurrentEventSelected == 1){
              fHistoMCEtaPtNotTriggered[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
            } else if(isCurrentEventSelected == 2){
              fHistoMCEtaPtNoVertex[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
            }
            fHistoMCEtaWOWeightPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
            if (fIsMC > 1) fHistoMCEtaWOEvtWeightPt[fiCut]->Fill(particle->Pt());
            if (fDoMesonQA > 0){
              if (fIsMC < 2)fHistoMCEtaPtY[fiCut]->Fill(particle->Pt(),mesonY,weighted); // All MC Pi0
              fHistoMCEtaPtAlpha[fiCut]->Fill(particle->Pt(),alpha,fWeightJetJetMC); // All MC Pi0
              if (fIsMC == 2) fHistoMCEtaPtJetPt[fiCut]->Fill(particle->Pt(),fiEventCut->GetMaxPtJet(),fWeightJetJetMC);
            }
          }

          // Check the acceptance for both gammas & whether they are counted as primaries as well
          Bool_t kDaughter0IsPrim = fiEventCut->IsConversionPrimaryESD( fMCEvent, particle->GetDaughterFirst(), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
          Bool_t kDaughter1IsPrim = fiEventCut->IsConversionPrimaryESD( fMCEvent, particle->GetDaughterLast(), mcProdVtxX, mcProdVtxY, mcProdVtxZ);

          if( kDaughter0IsPrim && kDaughter1IsPrim &&
              fiPhotonCut->PhotonIsSelectedMC(daughter0,fMCEvent,kFALSE) &&
              fiPhotonCut->PhotonIsSelectedMC(daughter1,fMCEvent,kFALSE)  &&
              fiPhotonCut->InPlaneOutOfPlaneCut(daughter0->Phi(),fEventPlaneAngle,kFALSE) &&
              fiPhotonCut->InPlaneOutOfPlaneCut(daughter1->Phi(),fEventPlaneAngle,kFALSE)){

            if(particle->PdgCode() == 111){
              if(fIsMC > 1) fHistoMCPi0WOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Pi0 with gamma in acc NOT weighted at all
              fHistoMCPi0WOWeightInAccPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // MC Pi0 with gamma in acc NOT weighted
              fHistoMCPi0InAccPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // MC Pi0 with gamma in acc
              if(isCurrentEventSelected == 1){
                fHistoMCPi0InAccPtNotTriggered[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
              } else if(isCurrentEventSelected == 2){
                fHistoMCPi0InAccPtNoVertex[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
              }
            
            } else if(particle->PdgCode() == 221){
              if(fIsMC > 1) fHistoMCEtaWOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Eta with gamma in acc NOT weighted at all
              fHistoMCEtaWOWeightInAccPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // MC Eta with gamma in acc NOT weighted
              fHistoMCEtaInAccPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // MC Eta with gamma in acc
              if(isCurrentEventSelected == 1){
                fHistoMCEtaInAccPtNotTriggered[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
              } else if(isCurrentEventSelected == 2){
                fHistoMCEtaInAccPtNoVertex[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
              }
            }
          }
        }
      }
    } else {
      // fill secondary histograms
      AliMCParticle* particle = (AliMCParticle *)fMCEvent->GetTrack(i);
      if (!particle) continue;

      Int_t isMCFromMBHeader = -1;
      if(fiEventCut->GetSignalRejection() != 0){
        isMCFromMBHeader = fiEventCut->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
        if(isMCFromMBHeader == 0 && fiEventCut->GetSignalRejection() != 3) continue;
      }

      if(fiPhotonCut->InPlaneOutOfPlaneCut(particle->Phi(),fEventPlaneAngle,kFALSE)){
        if(fiPhotonCut->PhotonIsSelectedMC(particle,fMCEvent,kFALSE)){
          if (particle->GetMother() > -1 && fMCEvent->GetTrack(particle->GetMother())->GetMother() > -1) {
            if (fMCEvent->GetTrack(fMCEvent->GetTrack(particle->GetMother())->GetMother())->PdgCode() == 310){
              fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),0.,fWeightJetJetMC);
            } else if (fMCEvent->GetTrack(fMCEvent->GetTrack(particle->GetMother())->GetMother())->PdgCode() == 130) {
              fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),1.,fWeightJetJetMC);
            } else if (fMCEvent->GetTrack(fMCEvent->GetTrack(particle->GetMother())->GetMother())->PdgCode() == 3122) {
              fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),2.,fWeightJetJetMC);
            } else {
              if ( !(TMath::Abs(fMCEvent->GetTrack(particle->GetMother())->PdgCode()) == 11 && fMCEvent->GetTrack(fMCEvent->GetTrack(particle->GetMother())->GetMother())->PdgCode() == 22) )
                fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),3.,fWeightJetJetMC);
            }
          } else {
            fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),3.,fWeightJetJetMC);
          }
        }

        if(fiPhotonCut->PhotonIsSelectedMC(particle,fMCEvent,kTRUE)){
          if (particle->GetMother() > -1 && fMCEvent->GetTrack(particle->GetMother())->GetMother() > -1) {
            if (fMCEvent->GetTrack(fMCEvent->GetTrack(particle->GetMother())->GetMother())->PdgCode() == 310){
              fHistoMCSecondaryConvGammaPt[fiCut]->Fill(particle->Pt(),0.,fWeightJetJetMC);
            } else if (fMCEvent->GetTrack(fMCEvent->GetTrack(particle->GetMother())->GetMother())->PdgCode() == 130) {
              fHistoMCSecondaryConvGammaPt[fiCut]->Fill(particle->Pt(),1.,fWeightJetJetMC);
            } else if (fMCEvent->GetTrack(fMCEvent->GetTrack(particle->GetMother())->GetMother())->PdgCode() == 3122) {
              fHistoMCSecondaryConvGammaPt[fiCut]->Fill(particle->Pt(),2.,fWeightJetJetMC);
            } else {
              if ( !(TMath::Abs(fMCEvent->GetTrack(particle->GetMother())->PdgCode()) == 11 && fMCEvent->GetTrack(fMCEvent->GetTrack(particle->GetMother())->GetMother())->PdgCode() == 22) )
                fHistoMCSecondaryConvGammaPt[fiCut]->Fill(particle->Pt(),3.,fWeightJetJetMC);
            }
          } else {
            fHistoMCSecondaryConvGammaPt[fiCut]->Fill(particle->Pt(),3.,fWeightJetJetMC);
          }
        }
      }

      if(fDoMesonAnalysis){
        if(fiMesonCut->MesonIsSelectedMC(particle,fMCEvent,fiEventCut->GetEtaShift())){
          AliMCParticle* daughter0  = (AliMCParticle*)fMCEvent->GetTrack(particle->GetDaughterFirst());
          AliMCParticle* daughter1  = (AliMCParticle*)fMCEvent->GetTrack(particle->GetDaughterLast());
          Int_t pdgCode         = -1;
          if(particle->GetMother()>-1) pdgCode = ((AliMCParticle*)fMCEvent->GetTrack( particle->GetMother() ))->PdgCode();
          if(particle->PdgCode() == 111){
            Int_t source = GetSourceClassification(111,pdgCode);
            fHistoMCSecPi0PtvsSource[fiCut]->Fill(particle->Pt(),source,fWeightJetJetMC);

            Double_t deltaX = particle->Xv() - mcProdVtxX;
            Double_t deltaY = particle->Yv() - mcProdVtxY;
            Double_t deltaZ = particle->Zv() - mcProdVtxZ;
            Double_t realRadius3D = TMath::Sqrt(deltaX*deltaX+deltaY*deltaY+deltaZ*deltaZ);
            if(fDoMesonQA > 0 && fIsMC < 2) fHistoMCSecPi0RvsSource[fiCut]->Fill(realRadius3D,source);
            fHistoMCSecPi0Source[fiCut]->Fill(pdgCode);
          } else if(particle->PdgCode() == 221){
            fHistoMCSecEtaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
            fHistoMCSecEtaSource[fiCut]->Fill(pdgCode);
          }

          // pi0 really in acceptance/
          if( fiPhotonCut->PhotonIsSelectedMC(daughter0,fMCEvent,kFALSE) &&
              fiPhotonCut->PhotonIsSelectedMC(daughter1,fMCEvent,kFALSE)  &&
              fiPhotonCut->InPlaneOutOfPlaneCut(daughter0->Phi(),fEventPlaneAngle,kFALSE) &&
              fiPhotonCut->InPlaneOutOfPlaneCut(daughter1->Phi(),fEventPlaneAngle,kFALSE)){
            if(particle->PdgCode() == 111){
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

        if((fiMesonCut->MesonIsSelected(pi0cand,kTRUE,fiEventCut->GetEtaShift()))){
          if(fDoCentralityFlat > 0){
            fHistoMotherInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(), fWeightCentrality[fiCut]*fWeightJetJetMC);
            if(!fDoLightOutput) {
              fHistoMotherInvMassPtCalib[fiCut]->Fill(pi0cand->M(),gamma0->E(), fWeightCentrality[fiCut]*fWeightJetJetMC);
              fHistoMotherInvMassPtCalib[fiCut]->Fill(pi0cand->M(),gamma1->E(), fWeightCentrality[fiCut]*fWeightJetJetMC);
            }
            if(TMath::Abs(pi0cand->GetAlpha())<0.1) fHistoMotherInvMassEalpha[fiCut]->Fill(pi0cand->M(),pi0cand->E(), fWeightCentrality[fiCut]*fWeightJetJetMC);
          } else {
            if(TMath::Abs(pi0cand->GetAlpha())<0.1) fHistoMotherInvMassEalpha[fiCut]->Fill(pi0cand->M(),pi0cand->E(),fWeightJetJetMC);
            fHistoMotherInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(),fWeightJetJetMC);
            if(!fDoLightOutput) {
              fHistoMotherInvMassPtCalib[fiCut]->Fill(pi0cand->M(),gamma0->E(), fWeightJetJetMC);
              fHistoMotherInvMassPtCalib[fiCut]->Fill(pi0cand->M(),gamma1->E(), fWeightJetJetMC);
            }
            if(fDoIsolatedAnalysis){
              //Check if the pi0 is isolated
              Double_t Iso_E_Pi0 = 0;
              Double_t Iso_R_Pi0 = 0;
              Int_t Iso_counter = 0;
              //loop over all charged tracks
              for(Int_t iTracks = 0; iTracks<fInputEvent->GetNumberOfTracks(); iTracks++){
                AliAODTrack* curTrack = (AliAODTrack*) fInputEvent->GetTrack(iTracks);
                if(curTrack->GetID()<0) continue; // Avoid double counting of tracks
                if(!curTrack->IsHybridGlobalConstrainedGlobal()) continue;
                if(TMath::Abs(curTrack->Eta())>0.8) continue;
                if(curTrack->Pt()<0.15) continue;
                Double_t Iso_DeltaEta = curTrack->Eta()-pi0cand->Eta();
                Double_t Iso_DeltaPhi = abs(curTrack->Phi()-pi0cand->Phi());
                if(Iso_DeltaPhi > M_PI) {
                  Iso_DeltaPhi = 2*M_PI - Iso_DeltaPhi;
                }
                Iso_R_Pi0 = TMath::Sqrt(pow((Iso_DeltaEta),2)+pow((Iso_DeltaPhi),2));
                fHistoMotherRisoPt[fiCut]->Fill(Iso_R_Pi0,pi0cand->Pt());
                if(Iso_R_Pi0 < 0.4){
                  Iso_counter++;
                  Iso_E_Pi0 += curTrack->E();
                }
              }
              fHistoMotherEisoPt[fiCut]->Fill(Iso_E_Pi0,pi0cand->Pt());
              if(Iso_E_Pi0 < 2){
                //pi0 is isolated
                fHistoMotherInvMassPtIso[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(),fWeightJetJetMC);
                fHistoMotherNtracksIsoPt[fiCut]->Fill(Iso_counter,pi0cand->Pt());
              }else{
                //pi0 is not isolated
                fHistoMotherInvMassPtNonIso[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(),fWeightJetJetMC);
              }
              if(fIsMC>0){
                //Check for MC true if the pi0 should be isolated or not
                Double_t Iso_E_Pi0_True = 0;
                Iso_R_Pi0 = 0;
                //loop over all charged tracks
		//                TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
		// if (AODMCTrackArray == NULL) continue;
                for(Long_t iTracks = 0; iTracks < fAODMCTrackArray->GetEntriesFast(); iTracks++) {
                  AliAODMCParticle* curTrack = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(iTracks));
                  if (!curTrack) continue;
                  if (!curTrack->IsPhysicalPrimary()) continue;
                  if (!curTrack->Pt()) continue;
                  if(curTrack->Charge()==0) continue;
                  if(curTrack->Pt()<0.15) continue;
                  if(TMath::Abs(curTrack->Eta())>0.8) continue;
                  Double_t Iso_DeltaEta = curTrack->Eta()-pi0cand->Eta();
                  Double_t Iso_DeltaPhi = abs(curTrack->Phi()-pi0cand->Phi());
                  if(Iso_DeltaPhi > M_PI) {
                    Iso_DeltaPhi = 2*M_PI - Iso_DeltaPhi;
                  }
                  Iso_R_Pi0 = TMath::Sqrt(pow((Iso_DeltaEta),2)+pow((Iso_DeltaPhi),2));
                  if(Iso_R_Pi0 < 0.4){
                    Iso_E_Pi0_True += curTrack->E();
                  }
                }
                if(Iso_E_Pi0 < 2 && Iso_E_Pi0_True > 2){
                  //pi0 rec is isolated, but using true particles not
                  fHistoMotherInvMassPtMCRecIsoTrueNonIso[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(),fWeightJetJetMC);
                }else if(Iso_E_Pi0 > 2 && Iso_E_Pi0_True < 2){
                  //pi0 rec is not isolated, but using true particles it is
                  fHistoMotherInvMassPtMCRecNonIsoTrueIso[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(),fWeightJetJetMC);
                }
              }
            }
          }

          if (fDoMesonQA > 0){

            if(fDoMesonQA == 3 && TMath::Abs(gamma0->GetConversionRadius()-gamma1->GetConversionRadius())<10 && pi0cand->GetOpeningAngle()<0.1){
                    Double_t sparesFill[4] = {gamma0->GetPhotonPt(),gamma0->GetConversionRadius(),TMath::Abs(gamma0->GetConversionRadius()-gamma1->GetConversionRadius()),pi0cand->GetOpeningAngle()};
                    sPtRDeltaROpenAngle[fiCut]->Fill(sparesFill, 1);
            }

            if ( pi0cand->M() > 0.05 && pi0cand->M() < 0.17){
              if (fIsMC < 2){
                fHistoMotherPi0PtY[fiCut]->Fill(pi0cand->Pt(),pi0cand->Rapidity()-fiEventCut->GetEtaShift());
                fHistoMotherPi0PtOpenAngle[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetOpeningAngle());
              }
              fHistoMotherPi0PtAlpha[fiCut]->Fill(pi0cand->Pt(),TMath::Abs(pi0cand->GetAlpha()),fWeightJetJetMC);

            }
            if ( pi0cand->M() > 0.45 && pi0cand->M() < 0.65){
              if (fIsMC < 2){
                fHistoMotherEtaPtY[fiCut]->Fill(pi0cand->Pt(),pi0cand->Rapidity()-fiEventCut->GetEtaShift());
                fHistoMotherEtaPtOpenAngle[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetOpeningAngle());
              }
              fHistoMotherEtaPtAlpha[fiCut]->Fill(pi0cand->Pt(),TMath::Abs(pi0cand->GetAlpha()),fWeightJetJetMC);
            }
          }
          if(fDoTHnSparse && fiMesonCut->DoBGCalculation()){
            Int_t psibin = 0;
            Int_t zbin = 0;
            Int_t mbin = 0;

            Double_t sparesFill[4];
            if(fiMesonCut->BackgroundHandlerType() == 0){
              zbin = fBGHandler[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
              if(fiMesonCut->UseTrackMultiplicity()){
                mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
              } else {
                mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fGammaCandidates->GetEntries());
              }
              sparesFill[0] = pi0cand->M();
              sparesFill[1] = pi0cand->Pt();
              sparesFill[2] = (Double_t)zbin;
              sparesFill[3] = (Double_t)mbin;
            } else if(fiMesonCut->BackgroundHandlerType() != 2){
              psibin = fBGHandlerRP[fiCut]->GetRPBinIndex(TMath::Abs(fEventPlaneAngle));
              zbin = fBGHandlerRP[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
//               if(fiMesonCut->UseTrackMultiplicity()){
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
            } else if (fiPhotonCut->GetSingleElectronPtCut() < 0.04 && fPt > 0.099 && fPt < 20. )  {
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
  Double_t magField = fInputEvent->GetMagneticField();
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
      AliMCParticle * negativeMC = (AliMCParticle*)TrueGammaCandidate0->GetNegativeMCDaughter(fMCEvent);
      AliMCParticle * positiveMC = (AliMCParticle*)TrueGammaCandidate0->GetPositiveMCDaughter(fMCEvent);
      AliMCParticle * gammaMC0 = (AliMCParticle*)fMCEvent->GetTrack(gamma0MCLabel);
      if(TMath::Abs(negativeMC->PdgCode())==11 && TMath::Abs(positiveMC->PdgCode())==11){  // Electrons ...
        if(negativeMC->Particle()->GetUniqueID() == 5 && positiveMC->Particle()->GetUniqueID() ==5){ // ... From Conversion ...
          if(gammaMC0->PdgCode() == 22){ // ... with Gamma Mother
            gamma0MotherLabel=gammaMC0->GetMother();
          }
        }
        if(gammaMC0->PdgCode() ==111){ // Dalitz candidate
          gamma0DalitzCand = kTRUE;
          gamma0MotherLabel=-111;
        }
        if(gammaMC0->PdgCode() ==221){ // Dalitz candidate
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
        AliMCParticle * negativeMC = (AliMCParticle*)TrueGammaCandidate1->GetNegativeMCDaughter(fMCEvent);
        AliMCParticle * positiveMC = (AliMCParticle*)TrueGammaCandidate1->GetPositiveMCDaughter(fMCEvent);
        AliMCParticle * gammaMC1 = (AliMCParticle*)fMCEvent->GetTrack(gamma1MCLabel);
        if(TMath::Abs(negativeMC->PdgCode())==11 && TMath::Abs(positiveMC->PdgCode())==11){  // Electrons ...
          if(negativeMC->Particle()->GetUniqueID() == 5 && positiveMC->Particle()->GetUniqueID() ==5){ // ... From Conversion ...
            if(gammaMC1->PdgCode() == 22){ // ... with Gamma Mother
              gamma1MotherLabel=gammaMC1->GetMother();
            }
          }
          if(gammaMC1->PdgCode() ==111 ){ // Dalitz candidate
            gamma1DalitzCand = kTRUE;
            gamma1MotherLabel=-111;
          }
          if(gammaMC1->PdgCode() ==221){ // Dalitz candidate
            gamma1DalitzCand = kTRUE;
            gamma1MotherLabel=-221;
          }
        }
      }
      if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
        if(((AliMCParticle*)fMCEvent->GetTrack(gamma1MotherLabel))->PdgCode() == 111){
          isTruePi0=kTRUE;
          if (CheckVectorForDoubleCount(vecDoubleCountTruePi0s,gamma0MotherLabel)){
            fHistoDoubleCountTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC);
            FillMultipleCountMap(mapMultipleCountTruePi0s,gamma0MotherLabel);
          }
        }
        if(((AliMCParticle*)fMCEvent->GetTrack(gamma1MotherLabel))->PdgCode() == 221){
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
        if (fDoMaterialBudgetWeightingOfGammasForTrueMesons && fiPhotonCut->GetMaterialBudgetWeightsInitialized()) {
	      weightMatBudget = fiPhotonCut->GetMaterialBudgetCorrectingWeightForTrueGamma(TrueGammaCandidate0, magField) * fiPhotonCut->GetMaterialBudgetCorrectingWeightForTrueGamma(TrueGammaCandidate1,magField);
        }

        fHistoTrueMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightMatBudget*fWeightJetJetMC);
        if (fDoMesonQA > 0){
          if (isTruePi0){
            if ( Pi0Candidate->M() > 0.05 && Pi0Candidate->M() < 0.17){
              if(fIsMC < 2){
                fHistoTruePi0PtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-fiEventCut->GetEtaShift());
                fHistoTruePi0PtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle());
              }
              fHistoTruePi0PtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),TMath::Abs(Pi0Candidate->GetAlpha()),fWeightJetJetMC);

            }
          } else if (isTrueEta){
            if ( Pi0Candidate->M() > 0.45 && Pi0Candidate->M() < 0.65){
              if(fIsMC < 2){
                fHistoTrueEtaPtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-fiEventCut->GetEtaShift());
                fHistoTrueEtaPtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle());
              }
              fHistoTrueEtaPtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),TMath::Abs(Pi0Candidate->GetAlpha()),fWeightJetJetMC);
            }
          }
        }

        Bool_t isPrimary = fiEventCut->IsConversionPrimaryESD( fMCEvent, gamma0MotherLabel, mcProdVtxX, mcProdVtxY, mcProdVtxZ);

        if(!isPrimary && gamma1MotherLabel>-1){ // Secondary Meson
          Long_t secMotherLabel = ((AliMCParticle*)fMCEvent->GetTrack(gamma1MotherLabel))->GetMother();
          Float_t weightedSec= 1;
          if(fMCEvent->GetTrack(secMotherLabel)->PdgCode()==310){
            weightedSec= fiEventCut->GetWeightForMeson(secMotherLabel, fMCEvent, fInputEvent)/2.; //invariant mass is additive thus the weight for the daughters has to be devide by two for the K0s at a certain pt
            //cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
          }
          fHistoTrueSecondaryMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
          iMesonMCInfo = 2;
          if (secMotherLabel >-1){
            if(fMCEvent->GetTrack(secMotherLabel)->PdgCode()==310){
              iMesonMCInfo = 4;
              fHistoTrueSecondaryMotherFromK0sInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
              if (fDoMesonQA > 0 && fIsMC < 2 )fHistoTrueK0sWithPi0DaughterMCPt[fiCut]->Fill(fMCEvent->GetTrack(secMotherLabel)->Pt());
            }
            if(fMCEvent->GetTrack(secMotherLabel)->PdgCode()==130){
              iMesonMCInfo = 8;
              fHistoTrueSecondaryMotherFromK0lInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
              if (fDoMesonQA > 0 && fIsMC < 2 )fHistoTrueK0lWithPi0DaughterMCPt[fiCut]->Fill(fMCEvent->GetTrack(secMotherLabel)->Pt());
            }
            if(fMCEvent->GetTrack(secMotherLabel)->PdgCode()==221){
              iMesonMCInfo = 3;
              fHistoTrueSecondaryMotherFromEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*weightMatBudget*fWeightJetJetMC);
              if (fDoMesonQA > 0 && fIsMC < 2)fHistoTrueEtaWithPi0DaughterMCPt[fiCut]->Fill(fMCEvent->GetTrack(secMotherLabel)->Pt());
            }
            if(fMCEvent->GetTrack(secMotherLabel)->PdgCode()==3122){
              iMesonMCInfo = 7;
              fHistoTrueSecondaryMotherFromLambdaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
              if (fDoMesonQA > 0 && fIsMC < 2)fHistoTrueLambdaWithPi0DaughterMCPt[fiCut]->Fill(fMCEvent->GetTrack(secMotherLabel)->Pt());
            }
          }
        } else { // Only primary pi0 for efficiency calculation
          iMesonMCInfo = 6;
          Float_t weighted= 1;
          if (((AliMCParticle*)fMCEvent->GetTrack(gamma1MotherLabel))->Pt()>0.005){
            weighted= fiEventCut->GetWeightForMeson(gamma1MotherLabel, fMCEvent, fInputEvent);
          }
          fHistoTruePrimaryMotherW0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC);
          pESDTruePrimaryMotherWeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*weightMatBudget*fWeightJetJetMC);
          fHistoTruePrimaryMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudget);
          if (fDoMesonQA > 0 && fIsMC < 2){
            if(isTruePi0){ // Only primary pi0 for resolution
              fHistoTruePrimaryPi0MCPtResolPt[fiCut]->Fill(((AliMCParticle*)fMCEvent->GetTrack(gamma1MotherLabel))->Pt(),(Pi0Candidate->Pt()-((AliMCParticle*)fMCEvent->GetTrack(gamma1MotherLabel))->Pt())/((AliMCParticle*)fMCEvent->GetTrack(gamma1MotherLabel))->Pt(),weighted*weightMatBudget);
            }
            if (isTrueEta){ // Only primary eta for resolution
              fHistoTruePrimaryEtaMCPtResolPt[fiCut]->Fill(((AliMCParticle*)fMCEvent->GetTrack(gamma1MotherLabel))->Pt(),(Pi0Candidate->Pt()-((AliMCParticle*)fMCEvent->GetTrack(gamma1MotherLabel))->Pt())/((AliMCParticle*)fMCEvent->GetTrack(gamma1MotherLabel))->Pt(),weighted*weightMatBudget);
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
  Double_t magField = fInputEvent->GetMagneticField();
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  // Process True Mesons
  //  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  Bool_t isTruePi0 = kFALSE;
  Bool_t isTrueEta = kFALSE;
  Bool_t isTruePi0Dalitz = kFALSE;
  Bool_t isTrueEtaDalitz = kFALSE;
  Bool_t gamma0DalitzCand = kFALSE;
  Bool_t gamma1DalitzCand = kFALSE;

  if (fAODMCTrackArray!=NULL && TrueGammaCandidate0 != NULL){
    AliAODMCParticle *positiveMC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelPositive()));
    AliAODMCParticle *negativeMC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelNegative()));

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
      AliAODMCParticle * gammaMC0 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma0MCLabel));
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
    positiveMC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(TrueGammaCandidate1->GetMCLabelPositive()));
    negativeMC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(TrueGammaCandidate1->GetMCLabelNegative()));

    Int_t gamma1MCLabel = -1;
    Int_t gamma1MotherLabel = -1;
    if(!positiveMC||!negativeMC)
      return;

    if(positiveMC->GetMother()>-1&&(negativeMC->GetMother() == positiveMC->GetMother())){
      gamma1MCLabel = positiveMC->GetMother();
    }
    if(gamma1MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
      // Daughters Gamma 1
      AliAODMCParticle * gammaMC1 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MCLabel));
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
      if(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 111){
        isTruePi0=kTRUE;
        if (CheckVectorForDoubleCount(vecDoubleCountTruePi0s,gamma0MotherLabel)){
          fHistoDoubleCountTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC);
          FillMultipleCountMap(mapMultipleCountTruePi0s,gamma0MotherLabel);
        }
      }
      if(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 221){
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
      if (fDoMaterialBudgetWeightingOfGammasForTrueMesons && fiPhotonCut->GetMaterialBudgetWeightsInitialized()) {
	    weightMatBudget = fiPhotonCut->GetMaterialBudgetCorrectingWeightForTrueGamma(TrueGammaCandidate0,magField) * fiPhotonCut->GetMaterialBudgetCorrectingWeightForTrueGamma(TrueGammaCandidate1,magField);
      }

      fHistoTrueMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightMatBudget*fWeightJetJetMC);
      if (fDoMesonQA > 0){
        if (isTruePi0){
          if ( Pi0Candidate->M() > 0.05 && Pi0Candidate->M() < 0.17){
            if(fIsMC < 2){
              fHistoTruePi0PtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-fiEventCut->GetEtaShift());
              fHistoTruePi0PtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle());
            }
            fHistoTruePi0PtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),TMath::Abs(Pi0Candidate->GetAlpha()),fWeightJetJetMC);
          }
        } else if (isTrueEta){
          if ( Pi0Candidate->M() > 0.45 && Pi0Candidate->M() < 0.65){
            if(fIsMC < 2){
              fHistoTrueEtaPtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-fiEventCut->GetEtaShift());
              fHistoTrueEtaPtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle());
            }
            fHistoTrueEtaPtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),TMath::Abs(Pi0Candidate->GetAlpha()),fWeightJetJetMC);
          }
        }
      }
      Bool_t isPrimary = fiEventCut->IsConversionPrimaryAOD(fInputEvent, static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma0MotherLabel)), mcProdVtxX, mcProdVtxY, mcProdVtxZ);

      if(!isPrimary){ // Secondary Meson
        Long_t secMotherLabel = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MotherLabel))->GetMother();
        Float_t weightedSec= 1;
        if(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->GetPdgCode()==310){
          weightedSec= fiEventCut->GetWeightForMeson(secMotherLabel, 0x0, fInputEvent)/2.; //invariant mass is additive thus the weight for the daughters has to be devide by two for the K0s at a certain pt
          //cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
        }
        fHistoTrueSecondaryMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
        iMesonMCInfo = 2;
        if (secMotherLabel >-1){
          if(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->GetPdgCode()==310){
            iMesonMCInfo = 4;
            fHistoTrueSecondaryMotherFromK0sInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
            if (fDoMesonQA > 0 && fIsMC < 2)fHistoTrueK0sWithPi0DaughterMCPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->Pt());
          }
          if(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->GetPdgCode()==130){
            iMesonMCInfo = 8;
            fHistoTrueSecondaryMotherFromK0lInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
            if (fDoMesonQA > 0 && fIsMC < 2)fHistoTrueK0lWithPi0DaughterMCPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->Pt());
          }
          if(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->GetPdgCode()==221){
            iMesonMCInfo = 3;
            fHistoTrueSecondaryMotherFromEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*weightMatBudget*fWeightJetJetMC);
            if (fDoMesonQA > 0 && fIsMC < 2)fHistoTrueEtaWithPi0DaughterMCPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->Pt());
          }
          if(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->GetPdgCode()==3122){
            iMesonMCInfo = 7;
            fHistoTrueSecondaryMotherFromLambdaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*fWeightJetJetMC*weightMatBudget);
            if (fDoMesonQA > 0 && fIsMC < 2)fHistoTrueLambdaWithPi0DaughterMCPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->Pt());
          }
        }
      } else { // Only primary pi0 for efficiency calculation
        Float_t weighted= 1;
        iMesonMCInfo = 6;
        if (static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MotherLabel))->Pt()>0.005){
          weighted= fiEventCut->GetWeightForMeson(gamma1MotherLabel, 0x0, fInputEvent);
        }
        fHistoTruePrimaryMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*fWeightJetJetMC*weightMatBudget);
        fHistoTruePrimaryMotherW0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),fWeightJetJetMC);
        pESDTruePrimaryMotherWeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*weightMatBudget*fWeightJetJetMC);

        if (fDoMesonQA > 0 && fIsMC < 2){
          if(isTruePi0){ // Only primary pi0 for resolution
            fHistoTruePrimaryPi0MCPtResolPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MotherLabel))->Pt(),
                                (Pi0Candidate->Pt()-static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MotherLabel))->Pt())/static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MotherLabel))->Pt(),weighted*weightMatBudget);

          }
          if (isTrueEta){ // Only primary eta for resolution
            fHistoTruePrimaryEtaMCPtResolPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MotherLabel))->Pt(),
                                (Pi0Candidate->Pt()-static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MotherLabel))->Pt())/static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MotherLabel))->Pt(),weighted*weightMatBudget);
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

    if(fiMesonCut->UseTrackMultiplicity()){
        mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
    } else {
        mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fGammaCandidates->GetEntries());
    }

  if(fiMesonCut->UseRotationMethod()){

    for(Int_t iCurrent=0;iCurrent<fGammaCandidates->GetEntries();iCurrent++){
      AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fGammaCandidates->At(iCurrent));
      for(Int_t iCurrent2=iCurrent+1;iCurrent2<fGammaCandidates->GetEntries();iCurrent2++){
        for(Int_t nRandom=0;nRandom<fiMesonCut->GetNumberOfBGEvents();nRandom++){
        AliAODConversionPhoton currentEventGoodV02 = *(AliAODConversionPhoton*)(fGammaCandidates->At(iCurrent2));

        if(fiMesonCut->DoBGProbability()){
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
        if((fiMesonCut->MesonIsSelected(backgroundCandidate,kFALSE,fiEventCut->GetEtaShift()))){
          if(fDoCentralityFlat > 0) fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(), fWeightCentrality[fiCut]*fWeightJetJetMC);
          else fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(),fWeightJetJetMC);
          if(!fDoLightOutput){
            if(fDoCentralityFlat > 0) fHistoMotherBackInvMassPtCalib[fiCut]->Fill(backgroundCandidate->M(),currentEventGoodV0.E(), fWeightCentrality[fiCut]*fWeightJetJetMC);
            else fHistoMotherBackInvMassPtCalib[fiCut]->Fill(backgroundCandidate->M(),currentEventGoodV0.E(),fWeightJetJetMC);
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
  } else {
    AliGammaConversionAODBGHandler::GammaConversionVertex *bgEventVertex = NULL;

    if(fiMesonCut->UseTrackMultiplicity()){
      for(Int_t nEventsInBG=0;nEventsInBG<fBGHandler[fiCut]->GetNBGEvents();nEventsInBG++){
        AliGammaConversionAODVector *previousEventV0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,nEventsInBG);
        if(fMoveParticleAccordingToVertex == kTRUE || fiPhotonCut->GetInPlaneOutOfPlaneCut() != 0){
          bgEventVertex = fBGHandler[fiCut]->GetBGEventVertex(zbin,mbin,nEventsInBG);
        }

        for(Int_t iCurrent=0;iCurrent<fGammaCandidates->GetEntries();iCurrent++){
        AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fGammaCandidates->At(iCurrent));
        for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
          AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));
          if(fMoveParticleAccordingToVertex == kTRUE){
            MoveParticleAccordingToVertex(&previousGoodV0,bgEventVertex);
          }
          if(fiPhotonCut->GetInPlaneOutOfPlaneCut() != 0){
            RotateParticleAccordingToEP(&previousGoodV0,bgEventVertex->fEP,fEventPlaneAngle);
          }

          AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
          backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
          if((fiMesonCut->MesonIsSelected(backgroundCandidate,kFALSE,fiEventCut->GetEtaShift()))){
            if(fDoCentralityFlat > 0) fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(), fWeightCentrality[fiCut]*fWeightJetJetMC);
            else fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(),fWeightJetJetMC);
            if(!fDoLightOutput){
              if(fDoCentralityFlat > 0) fHistoMotherBackInvMassPtCalib[fiCut]->Fill(backgroundCandidate->M(),currentEventGoodV0.E(), fWeightCentrality[fiCut]*fWeightJetJetMC);
              else fHistoMotherBackInvMassPtCalib[fiCut]->Fill(backgroundCandidate->M(),currentEventGoodV0.E(),fWeightJetJetMC);
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
    } else {
      for(Int_t nEventsInBG=0;nEventsInBG <fBGHandler[fiCut]->GetNBGEvents();nEventsInBG++){
        AliGammaConversionAODVector *previousEventV0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,nEventsInBG);
        if(previousEventV0s){
        if(fMoveParticleAccordingToVertex == kTRUE || fiPhotonCut->GetInPlaneOutOfPlaneCut() != 0){
          bgEventVertex = fBGHandler[fiCut]->GetBGEventVertex(zbin,mbin,nEventsInBG);
        }
        for(Int_t iCurrent=0;iCurrent<fGammaCandidates->GetEntries();iCurrent++){
          AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fGammaCandidates->At(iCurrent));
          for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){

            AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));

            if(fMoveParticleAccordingToVertex == kTRUE){
              MoveParticleAccordingToVertex(&previousGoodV0,bgEventVertex);
            }
            if(fiPhotonCut->GetInPlaneOutOfPlaneCut() != 0){
              RotateParticleAccordingToEP(&previousGoodV0,bgEventVertex->fEP,fEventPlaneAngle);
            }


            AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
            backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
            if((fiMesonCut->MesonIsSelected(backgroundCandidate,kFALSE,fiEventCut->GetEtaShift()))){
              if(fDoCentralityFlat > 0) fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(), fWeightCentrality[fiCut]*fWeightJetJetMC);
              else{
                fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(), fWeightJetJetMC);
                if(!fDoLightOutput) fHistoMotherBackInvMassPtCalib[fiCut]->Fill(backgroundCandidate->M(),currentEventGoodV0.E(), fWeightJetJetMC);  
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
void AliAnalysisTaskGammaConvV1::CalculateBackgroundSwapp(){

  if(fiMesonCut->DoGammaSwappForBg()) {

    Double_t rotationAngle = TMath::Pi()/2.0; //0.78539816339; // rotaion angle 90°

    TLorentzVector lvRotationPhoton1;   // photon candidates which get rotated
    TLorentzVector lvRotationPhoton2;   // photon candidates which get rotated
    TVector3 lvRotationPion;            // reconstructed mother particle from the two photons

    std::vector<std::array<Double_t, 2>> vSwappingInvMassPT;
    std::vector<std::array<Double_t, 2>> vSwappingInvMassPTAlphaCut;
    vSwappingInvMassPT.clear();
    vSwappingInvMassPTAlphaCut.clear();
    vSwappingInvMassPT.resize(0);
    vSwappingInvMassPTAlphaCut.resize(0);
    Double_t tempMultWeightSwapping = 1; // weight taking multiplicity of event into account

    // curcial requierment is that the event has at least 3 cluster candidates
    if(fGammaCandidates->GetEntries() > 2 ){
      for(Int_t iCurrent1=0;iCurrent1<fGammaCandidates->GetEntries();iCurrent1++){
        AliAODConversionPhoton* currentEventGoodV0Temp1 = (AliAODConversionPhoton*)(fGammaCandidates->At(iCurrent1));
        if(!currentEventGoodV0Temp1) continue;
        for(Int_t iCurrent2=iCurrent1+1; iCurrent2 < fGammaCandidates->GetEntries();iCurrent2++){
          AliAODConversionPhoton* currentEventGoodV0Temp2 = (AliAODConversionPhoton*)(fGammaCandidates->At(iCurrent2));
          if(!currentEventGoodV0Temp2) continue;
          for(int iSwapp = 0; iSwapp < fiMesonCut->GetNumberOfSwappsForBg(); ++iSwapp){

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
            if((fiMesonCut->GammaSwappMethodBg() == 0 || fiMesonCut->GammaSwappMethodBg() == 1)){
              if(fiMesonCut->GammaSwappMethodBg() == 0) rotationAngle = TMath::Pi()/2.0; // rotate by 90 degree
              else if(fiMesonCut->GammaSwappMethodBg() == 1){  // rotate by random angle between
                 Double_t temp = (fRandom.Rndm() < 0.5) ? 0 : TMath::Pi();
                 rotationAngle = temp + TMath::Pi()/3.0 + fRandom.Rndm()*TMath::Pi()/3.0;
              }
              lvRotationPhoton1.Rotate(rotationAngle, lvRotationPion);
              lvRotationPhoton2.Rotate(rotationAngle, lvRotationPion);
            }
            std::unique_ptr<AliAODConversionPhoton> currentEventGoodV0Rotation1 (new AliAODConversionPhoton(&lvRotationPhoton1));
            std::unique_ptr<AliAODConversionPhoton> currentEventGoodV0Rotation2 (new AliAODConversionPhoton(&lvRotationPhoton2));
            for(auto const& kCurrentGammaCandidates  : *fGammaCandidates){
              if(currentEventGoodV0Temp1 == ((AliAODConversionPhoton*) kCurrentGammaCandidates) || currentEventGoodV0Temp2 == ((AliAODConversionPhoton*) kCurrentGammaCandidates)) continue;

              std::unique_ptr<AliAODConversionMother> backgroundCandidate1(new AliAODConversionMother(currentEventGoodV0Rotation1.get(), ((AliAODConversionPhoton*) kCurrentGammaCandidates)));
              std::unique_ptr<AliAODConversionMother> backgroundCandidate2(new AliAODConversionMother(currentEventGoodV0Rotation2.get(), ((AliAODConversionPhoton*) kCurrentGammaCandidates)));
              if( fabs(currentEventGoodV0Temp1->Eta()) <= fiPhotonCut->GetEtaCut())
              {
                if(((AliConversionMesonCuts*) fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate1.get(),kFALSE,fiEventCut->GetEtaShift()))
                {
                  vSwappingInvMassPT.push_back({backgroundCandidate1->M(),backgroundCandidate1->Pt()});
                }
              }
              if( fabs(currentEventGoodV0Temp2->Eta()) <= fiPhotonCut->GetEtaCut())
              {
                if(((AliConversionMesonCuts*) fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate2.get(),kFALSE,fiEventCut->GetEtaShift()))
                {
                  vSwappingInvMassPT.push_back({backgroundCandidate2->M(),backgroundCandidate2->Pt()});
                }
              }
            }
          }
        }
      }
      // Fill the histograms
      if(fiMesonCut->DoWeightingInSwappBg() && vSwappingInvMassPT.size() > 0){
        tempMultWeightSwapping = (0.5*(fGammaCandidates->GetEntries()*fGammaCandidates->GetEntries() - fGammaCandidates->GetEntries()))/(vSwappingInvMassPT.size());
      }
      for(Int_t i = 0; i < (Int_t)vSwappingInvMassPT.size(); i++){
        if(fDoCentralityFlat > 0) fHistoMotherBackInvMassPt[fiCut]->Fill(vSwappingInvMassPT.at(i)[0], vSwappingInvMassPT.at(i)[1], tempMultWeightSwapping*fWeightCentrality[fiCut]*fWeightJetJetMC);
        else fHistoMotherBackInvMassPt[fiCut]->Fill(vSwappingInvMassPT.at(i)[0], vSwappingInvMassPT.at(i)[1], tempMultWeightSwapping*fWeightJetJetMC);
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
//     if(fiMesonCut->UseTrackMultiplicity()){
//       mbin = fBGHandlerRP[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
//     } else {
//       mbin = fBGHandlerRP[fiCut]->GetMultiplicityBinIndex(fGammaCandidates->GetEntries());
//     }
  }

  //Rotation Method
  if(fiMesonCut->UseRotationMethod()){
    // Correct for the number of rotations
    // BG is for rotation the same, except for factor NRotations
    Double_t weight=1./Double_t(fiMesonCut->GetNumberOfBGEvents());

    for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries();firstGammaIndex++){

      AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));
      if (gamma0==NULL) continue;
      for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fGammaCandidates->GetEntries();secondGammaIndex++){
        AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(secondGammaIndex));
        if (gamma1 == NULL) continue;
        if(!fiPhotonCut->PhotonIsSelected(gamma1,fInputEvent))continue;
        for(Int_t nRandom=0;nRandom<fiMesonCut->GetNumberOfBGEvents();nRandom++){
          RotateParticle(gamma1);
          AliAODConversionMother backgroundCandidate(gamma0,gamma1);
          backgroundCandidate.CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
          if(fiMesonCut->MesonIsSelected(&backgroundCandidate,kFALSE,fiEventCut->GetEtaShift())){
            if(fDoCentralityFlat > 0) fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate.M(),backgroundCandidate.Pt(), fWeightCentrality[fiCut]*fWeightJetJetMC);
            else fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate.M(),backgroundCandidate.Pt(),fWeightJetJetMC);
            if(!fDoLightOutput){
              if(fDoCentralityFlat > 0) fHistoMotherBackInvMassPtCalib[fiCut]->Fill(backgroundCandidate.M(),gamma0->E(), fWeightCentrality[fiCut]*fWeightJetJetMC);
              else fHistoMotherBackInvMassPtCalib[fiCut]->Fill(backgroundCandidate.M(),gamma0->E(),fWeightJetJetMC);
            }
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
            if(fiMesonCut->MesonIsSelected(&backgroundCandidate,kFALSE,fiEventCut->GetEtaShift())){
              if(fDoCentralityFlat > 0) fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate.M(),backgroundCandidate.Pt(), fWeightCentrality[fiCut]*fWeightJetJetMC);
              else fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate.M(),backgroundCandidate.Pt(),fWeightJetJetMC);
              if(!fDoLightOutput){
                if(fDoCentralityFlat > 0) fHistoMotherBackInvMassPtCalib[fiCut]->Fill(backgroundCandidate.M(),gamma0->E(), fWeightCentrality[fiCut]*fWeightJetJetMC);
                else fHistoMotherBackInvMassPtCalib[fiCut]->Fill(backgroundCandidate.M(),gamma0->E(),fWeightJetJetMC);
              }
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
  Int_t fNDegreesPMBackground= fiMesonCut->NDegreesRotation();
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
  if(fGammaCandidates->GetEntries() >1 ){
    if(fiMesonCut->UseTrackMultiplicity()){
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
void AliAnalysisTaskGammaConvV1::FillPhotonCombinatorialMothersHistESD(AliMCParticle *daughter, AliMCParticle *motherCombPart)
{
  Int_t pdgCombPart = TMath::Abs(daughter->PdgCode());
  Int_t pdgMotherCombPart = TMath::Abs(motherCombPart->PdgCode());

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
    PhotonCandidate->SetIsCaloPhoton(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType());
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
