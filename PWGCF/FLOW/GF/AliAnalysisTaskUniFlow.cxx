/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// =================================================================================================
// AliAnalysisTaskUniFlow - ALICE Unified Flow framework
// Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2016-2018
// =================================================================================================
//
// ALICE analysis task for universal study of flow via 2-(or multi-)particle correlations
// using Generic Framework notation for calculations including per-particle weights.
//
// Implemented flow calculation for both reference & pt-differential flow
// of both inclusive charged (done anyway) & identified particles (pi,K,p,K0s,Lambda,phi).
//
// Note: So far implemented only for AOD analysis and tuned on Run2 pp & pPb analyses!
//
// PLEASE READ THE INSTRUCTION BELLOW BEFORE RUNNING !!!
//
// =================================================================================================
// Analysis can run in these modes setup via AliAnalysisTaskUniFlow::SetRunMode(RunMode)
//  -- AnalysisTaskUniFlow::kFull : running mode
//      - full scale analysis
//
//  -- AnalysisTaskUniFlow::kTest : development / testing / debugging mode
//      - only limited number of events is processed (AliAnalysisTaskUniFlow::SetNumEventsAnalyse(Int_t))
//
//  -- AnalysisTaskUniFlow::kSkipFlow : usable for QA and(or) weights estimation before full scale running
//      - events are processed whilst skipping correlation calculations, i.e. event loop ends after particles are filtered
//
//  NOTE: by default, all modes includes :
//  - filling QA plots : can be turned-off by AliAnalysisTaskUniFlow::SetFillQAhistos(kFALSE)
//    (might be heavy while running several tasks at once)
//  - filling GF weights : can be turned-off by AliAnalysisTaskUniFlow::SetFlowFillWeights(kFALSE)
//
// =================================================================================================
// Overview of analysis flow (see implementation of corresonding method for details)
// 1) Event selection : EventSelection()
//
// 2) Particle selection & reconstruction : Filtering()
//        - whether or not are particles processed is driven by 'fProcess?' flag setup by AliAnalysisTaskUniFlow::SetProcess?(kTRUE)
//          (except for incl. charged, which are processed anyway)
//        - filling QA plots
//        - filling GF weights
//
//    !!! here the event loop ends in case of running in 'kSkipFlow' mode
//
// 3) Flow / correlation calculations : DoFlowPOIs(), DoFlowRefs()
//        - to setup harmonics to be calculated, modify 'AliAnalysisTaskUniFlow::fHarmonics[]' (in .cxx) AND 'fNumHarmonics' (in .h)
//        - to setup eta gaps (2part) to be calculated, modify 'AliAnalysisTaskUniFlow::fEtaGap[]' (in .cxx) AND 'fNumEtaGap' (in .h)
//              - for NO gap case, fill in value '-1.0'
//        - 2-particle cumulants are run by default (in 'kFull' mode)
//        - 4-particle cumulants has to be setup by invoking AliAnalysisTaskUniFlow::SetFlowDoFourCorrelations(kTRUE)
//
// =================================================================================================

#include <TDatabasePDG.h>
#include <TPDGCode.h>

#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TList.h"
#include "TClonesArray.h"
#include "TComplex.h"
#include "TRandom3.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisUtils.h"
#include "AliMultSelection.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliLog.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliVTrack.h"
#include "AliPicoTrack.h"
#include "AliAODv0.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"

#include "AliAnalysisTaskUniFlow.h"

class AliAnalysisTaskUniFlow;

ClassImp(AliAnalysisTaskUniFlow);

Int_t AliAnalysisTaskUniFlow::fHarmonics[] = {2};
Double_t AliAnalysisTaskUniFlow::fEtaGap[] = {0.0,0.4,0.8};
Double_t AliAnalysisTaskUniFlow::fMultBins[] = {0.,5.,10.,20.,40.,60.,100.};

AliAnalysisTaskUniFlow::AliAnalysisTaskUniFlow() : AliAnalysisTaskSE(),
  fPDGMassPion(TDatabasePDG::Instance()->GetParticle(211)->Mass()),
  fPDGMassKaon(TDatabasePDG::Instance()->GetParticle(321)->Mass()),
  fPDGMassProton(TDatabasePDG::Instance()->GetParticle(2212)->Mass()),
  fPDGMassPhi(TDatabasePDG::Instance()->GetParticle(333)->Mass()),
  fPDGMassK0s(TDatabasePDG::Instance()->GetParticle(310)->Mass()),
  fPDGMassLambda(TDatabasePDG::Instance()->GetParticle(3122)->Mass()),
  fFlowPOIsPtMin(0.0),
  fFlowPOIsPtMax(15.0),
  fFlowCentMin(0),
  fFlowCentMax(0),

  fEventAOD(0x0),
  fPVz(0.0),
  fPIDResponse(0x0),
  fPIDCombined(0x0),
  fFlowWeightsFile(0x0),
  fInit(kFALSE),
  fMC(kFALSE),
  fArrayMC(0x0),
  fIndexSampling(0),
  fIndexCentrality(-1),
  fEventCounter(0),
  fNumEventsAnalyse(50),
  fRunNumber(-1),

  // FlowPart containers
  fVectorRefs(0x0),
  fVectorCharged(0x0),
  fVectorPion(0x0),
  fVectorKaon(0x0),
  fVectorProton(0x0),
  fVectorK0s(0x0),
  fVectorLambda(0x0),
  fVectorPhi(0x0),

  // analysis selection
  fRunMode(kFull),
  fAnalType(kAOD),
  fSampling(kTRUE),
  fFillQA(kTRUE),
  // fNumSamples(10),
  fProcessPID(kTRUE),
  fProcessV0s(kTRUE),
  fProcessPhi(kTRUE),

  // flow related
  fUseFixedMultBins(kFALSE),
  fCutFlowRFPsPtMin(0.2),
  fCutFlowRFPsPtMax(5.0),
  fCutFlowDoFourCorrelations(kFALSE),
  fFlowFillWeights(kTRUE),
  fFlowUseWeights(kFALSE),
  fFlowUse3Dweights(kFALSE),
  fFlowRunByRunWeights(kTRUE),
  fFlowWeightsPath(),

  // events selection
  fColSystem(kPPb),
  fTrigger(0),
  fMultEstimator("V0A"),
  fUseAliEventCuts(kFALSE),
  fPVtxCutZ(10.0),

  // charged tracks selection
  fCutChargedTrackFilterBit(96),
  fCutChargedNumTPCclsMin(70),
  fCutChargedEtaMax(0.8),
  fCutChargedDCAzMax(0.),
  fCutChargedDCAxyMax(0.),

  // PID tracks selection
  fCutPIDUseAntiProtonOnly(kFALSE),
  fCutPIDnSigmaCombinedTOFrejection(kTRUE),
  fCutPIDnSigmaPionMax(3.),
  fCutPIDnSigmaKaonMax(3.),
  fCutPIDnSigmaProtonMax(3.),
  fCutPIDnSigmaTPCRejectElectron(3.),
  fCutUseBayesPID(kTRUE),
  fCutPIDBayesPionMin(0.8),
  fCutPIDBayesKaonMin(0.8),
  fCutPIDBayesProtonMin(0.8),
  fCutPIDBayesRejectElectron(0.8),
  fCutPIDBayesRejectMuon(0.8),

  // V0s selection
  fCutV0sOnFly(kFALSE),
  fCutV0srefitTPC(kTRUE),
  fCutV0srejectKinks(kTRUE),
  fCutV0sDaughterNumTPCClsMin(0),
  fCutV0sDaughterNumTPCCrossMin(70),
  fCutV0sDaughterNumTPCFindMin(1),
  fCutV0sDaughterNumTPCClsPIDMin(0),
  fCutV0sDaughterRatioCrossFindMin(0.8),
  fCutV0sCrossMassRejection(kTRUE),
  fCutV0sCrossMassCutK0s(0.005),
  fCutV0sCrossMassCutLambda(0.010),
  fCutV0sDCAtoPVMin(0.06),
  fCutV0sDCAtoPVMax(0.),
  fCutV0sDCAtoPVzMax(0.),
  fCutV0sDCADaughtersMin(0.),
  fCutV0sDCADaughtersMax(1.0),
  fCutV0sDecayRadiusMin(0.5),
  fCutV0sDecayRadiusMax(200.0),
  fCutV0sDaughterFilterBit(0),
  fCutV0sDaughterPtMin(0.),
  fCutV0sDaughterPtMax(0.),
  fCutV0sDaughterEtaMax(0.8),
  fCutV0sMotherEtaMax(0.8),
  fCutV0sMotherRapMax(0.),
  fCutV0sCPAK0sMin(0.97),
  fCutV0sCPALambdaMin(0.995),
  fCutV0sNumTauK0sMax(7.46),
  fCutV0sNumTauLambdaMax(3.8),
  fCutV0sInvMassK0sMin(0.4),
  fCutV0sInvMassK0sMax(0.6),
  fCutV0sInvMassLambdaMin(1.08),
  fCutV0sInvMassLambdaMax(1.16),
  fCutV0sArmenterosAlphaK0sMin(0.2),
  fCutV0sArmenterosAlphaLambdaMax(0.),
  fCutV0sK0sPionNumTPCSigmaMax(5.0),
  fCutV0sLambdaPionNumTPCSigmaMax(5.0),
  fCutV0sLambdaProtonNumTPCSigmaMax(5.0),

  // phi selection
  fCutPhiMotherEtaMax(0.8),
  fCutPhiInvMassMin(0.99),
  fCutPhiInvMassMax(1.07),

  // output lists
  fQAEvents(0x0),
  fQACharged(0x0),
  fQAPID(0x0),
  fQAV0s(0x0),
  fQAPhi(0x0),
  fFlowWeights(0x0),
  fFlowRefs(0x0),
  fFlowCharged(0x0),
  fFlowPID(0x0),
  fFlowPhi(0x0),
  fFlowK0s(0x0),
  fFlowLambda(0x0),


  // flow histograms & profiles
  fh3WeightsRefs(0x0),
  fh3WeightsCharged(0x0),
  fh3WeightsPion(0x0),
  fh3WeightsKaon(0x0),
  fh3WeightsProton(0x0),
  fh3WeightsPhi(0x0),
  fh3WeightsK0s(0x0),
  fh3WeightsLambda(0x0),
  fh3AfterWeightsRefs(0x0),
  fh3AfterWeightsCharged(0x0),
  fh3AfterWeightsPion(0x0),
  fh3AfterWeightsKaon(0x0),
  fh3AfterWeightsProton(0x0),
  fh3AfterWeightsPhi(0x0),
  fh3AfterWeightsK0s(0x0),
  fh3AfterWeightsLambda(0x0),
  fh2WeightRefs(0x0),
  fh2WeightCharged(0x0),
  fh2WeightPion(0x0),
  fh2WeightKaon(0x0),
  fh2WeightProton(0x0),
  fh2WeightK0s(0x0),
  fh2WeightLambda(0x0),
  fh2WeightPhi(0x0),
  fh3WeightRefs(0x0),
  fh3WeightCharged(0x0),
  fh3WeightPion(0x0),
  fh3WeightKaon(0x0),
  fh3WeightProton(0x0),
  fh3WeightK0s(0x0),
  fh3WeightLambda(0x0),
  fh3WeightPhi(0x0),

  // event histograms
  fhEventSampling(0x0),
  fhEventCentrality(0x0),
  fh2EventCentralityNumRefs(0x0),
  fhEventCounter(0x0),

  // charged histogram
  fhRefsMult(0x0),
  fhRefsPt(0x0),
  fhRefsEta(0x0),
  fhRefsPhi(0x0),
  fpRefsMult(0x0),
  fhChargedCounter(0x0),

  // PID histogram
  fhPIDPionMult(0x0),
  fhPIDPionPt(0x0),
  fhPIDPionPhi(0x0),
  fhPIDPionEta(0x0),
  fhPIDPionCharge(0x0),
  fhPIDKaonMult(0x0),
  fhPIDKaonPt(0x0),
  fhPIDKaonPhi(0x0),
  fhPIDKaonEta(0x0),
  fhPIDKaonCharge(0x0),
  fhPIDProtonMult(0x0),
  fhPIDProtonPt(0x0),
  fhPIDProtonPhi(0x0),
  fhPIDProtonEta(0x0),
  fhPIDProtonCharge(0x0),

  fh2PIDPionTPCnSigmaPion(0x0),
  fh2PIDPionTOFnSigmaPion(0x0),
  fh2PIDPionTPCnSigmaKaon(0x0),
  fh2PIDPionTOFnSigmaKaon(0x0),
  fh2PIDPionTPCnSigmaProton(0x0),
  fh2PIDPionTOFnSigmaProton(0x0),
  fh2PIDPionBayesPion(0x0),
  fh2PIDPionBayesKaon(0x0),
  fh2PIDPionBayesProton(0x0),
  fh2PIDPionTPCdEdx(0x0),
  fh2PIDPionTOFbeta(0x0),

  fh2PIDKaonTPCdEdx(0x0),
  fh2PIDKaonTOFbeta(0x0),
  fh2PIDProtonTPCdEdx(0x0),
  fh2PIDProtonTOFbeta(0x0),
  fh2PIDKaonTPCnSigmaPion(0x0),
  fh2PIDKaonTOFnSigmaPion(0x0),
  fh2PIDKaonTPCnSigmaKaon(0x0),
  fh2PIDKaonTOFnSigmaKaon(0x0),
  fh2PIDKaonTPCnSigmaProton(0x0),
  fh2PIDKaonTOFnSigmaProton(0x0),
  fh2PIDKaonBayesPion(0x0),
  fh2PIDKaonBayesKaon(0x0),
  fh2PIDKaonBayesProton(0x0),

  fh2PIDProtonTPCnSigmaPion(0x0),
  fh2PIDProtonTOFnSigmaPion(0x0),
  fh2PIDProtonTPCnSigmaKaon(0x0),
  fh2PIDProtonTOFnSigmaKaon(0x0),
  fh2PIDProtonTPCnSigmaProton(0x0),
  fh2PIDProtonTOFnSigmaProton(0x0),
  fh2PIDProtonBayesPion(0x0),
  fh2PIDProtonBayesKaon(0x0),
  fh2PIDProtonBayesProton(0x0),

  fhMCRecoSelectedPionPt(0x0),
  fhMCRecoSelectedTruePionPt(0x0),
  fhMCRecoAllPionPt(0x0),
  fhMCGenAllPionPt(0x0),
  fhMCRecoSelectedKaonPt(0x0),
  fhMCRecoSelectedTrueKaonPt(0x0),
  fhMCRecoAllKaonPt(0x0),
  fhMCGenAllKaonPt(0x0),
  fhMCRecoSelectedProtonPt(0x0),
  fhMCRecoSelectedTrueProtonPt(0x0),
  fhMCRecoAllProtonPt(0x0),
  fhMCGenAllProtonPt(0x0),

  // phi histograms
  fhPhiCounter(0x0),
  fhPhiMult(0x0),
  fhPhiBGMult(0x0),
  fhPhiInvMass(0x0),
  fhPhiBGInvMass(0x0),
  fhPhiCharge(0x0),
  fhPhiBGCharge(0x0),
  fhPhiPt(0x0),
  fhPhiEta(0x0),
  fhPhiPhi(0x0),

  // V0s histogram
  fhV0sCounter(0x0),
  fhV0sCounterK0s(0x0),
  fhV0sCounterLambda(0x0),
  fhV0sInvMassK0s(0x0),
  fhV0sInvMassLambda(0x0),
  fhV0sCompetingInvMassK0s(0x0),
  fhV0sCompetingInvMassLambda(0x0)
{
  // default constructor, don't allocate memory here!
  // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskUniFlow::AliAnalysisTaskUniFlow(const char* name) : AliAnalysisTaskSE(name),
fPDGMassPion(TDatabasePDG::Instance()->GetParticle(211)->Mass()),
fPDGMassKaon(TDatabasePDG::Instance()->GetParticle(321)->Mass()),
fPDGMassProton(TDatabasePDG::Instance()->GetParticle(2212)->Mass()),
fPDGMassPhi(TDatabasePDG::Instance()->GetParticle(333)->Mass()),
fPDGMassK0s(TDatabasePDG::Instance()->GetParticle(310)->Mass()),
fPDGMassLambda(TDatabasePDG::Instance()->GetParticle(3122)->Mass()),
fFlowPOIsPtMin(0.0),
fFlowPOIsPtMax(15.0),
fFlowCentMin(0),
fFlowCentMax(0),

fEventAOD(0x0),
fPVz(0.0),
fPIDResponse(0x0),
fPIDCombined(0x0),
fFlowWeightsFile(0x0),
fInit(kFALSE),
fMC(kFALSE),
fArrayMC(0x0),
fIndexSampling(0),
fIndexCentrality(-1),
fEventCounter(0),
fNumEventsAnalyse(50),
fRunNumber(-1),

// FlowPart containers
fVectorRefs(0x0),
fVectorCharged(0x0),
fVectorPion(0x0),
fVectorKaon(0x0),
fVectorProton(0x0),
fVectorK0s(0x0),
fVectorLambda(0x0),
fVectorPhi(0x0),

// analysis selection
fRunMode(kFull),
fAnalType(kAOD),
fSampling(kTRUE),
fFillQA(kTRUE),
// fNumSamples(10),
fProcessPID(kTRUE),
fProcessV0s(kTRUE),
fProcessPhi(kTRUE),

// flow related
fUseFixedMultBins(kFALSE),
fCutFlowRFPsPtMin(0.2),
fCutFlowRFPsPtMax(5.0),
fCutFlowDoFourCorrelations(kFALSE),
fFlowFillWeights(kTRUE),
fFlowUseWeights(kFALSE),
fFlowUse3Dweights(kFALSE),
fFlowRunByRunWeights(kTRUE),
fFlowWeightsPath(),

// events selection
fColSystem(kPPb),
fTrigger(0),
fMultEstimator("V0A"),
fUseAliEventCuts(kFALSE),
fPVtxCutZ(10.0),

// charged tracks selection
fCutChargedTrackFilterBit(96),
fCutChargedNumTPCclsMin(70),
fCutChargedEtaMax(0.8),
fCutChargedDCAzMax(0.),
fCutChargedDCAxyMax(0.),

// PID tracks selection
fCutPIDUseAntiProtonOnly(kFALSE),
fCutPIDnSigmaCombinedTOFrejection(kTRUE),
fCutPIDnSigmaPionMax(3.),
fCutPIDnSigmaKaonMax(3.),
fCutPIDnSigmaProtonMax(3.),
fCutPIDnSigmaTPCRejectElectron(3.),
fCutUseBayesPID(kTRUE),
fCutPIDBayesPionMin(0.8),
fCutPIDBayesKaonMin(0.8),
fCutPIDBayesProtonMin(0.8),
fCutPIDBayesRejectElectron(0.8),
fCutPIDBayesRejectMuon(0.8),

// V0s selection
fCutV0sOnFly(kFALSE),
fCutV0srefitTPC(kTRUE),
fCutV0srejectKinks(kTRUE),
fCutV0sDaughterNumTPCClsMin(0),
fCutV0sDaughterNumTPCCrossMin(70),
fCutV0sDaughterNumTPCFindMin(1),
fCutV0sDaughterNumTPCClsPIDMin(0),
fCutV0sDaughterRatioCrossFindMin(0.8),
fCutV0sCrossMassRejection(kTRUE),
fCutV0sCrossMassCutK0s(0.005),
fCutV0sCrossMassCutLambda(0.010),
fCutV0sDCAtoPVMin(0.06),
fCutV0sDCAtoPVMax(0.),
fCutV0sDCAtoPVzMax(0.),
fCutV0sDCADaughtersMin(0.),
fCutV0sDCADaughtersMax(1.0),
fCutV0sDecayRadiusMin(0.5),
fCutV0sDecayRadiusMax(200.0),
fCutV0sDaughterFilterBit(0),
fCutV0sDaughterPtMin(0.),
fCutV0sDaughterPtMax(0.),
fCutV0sDaughterEtaMax(0.8),
fCutV0sMotherEtaMax(0.8),
fCutV0sMotherRapMax(0.),
fCutV0sCPAK0sMin(0.97),
fCutV0sCPALambdaMin(0.995),
fCutV0sNumTauK0sMax(7.46),
fCutV0sNumTauLambdaMax(3.8),
fCutV0sInvMassK0sMin(0.4),
fCutV0sInvMassK0sMax(0.6),
fCutV0sInvMassLambdaMin(1.08),
fCutV0sInvMassLambdaMax(1.16),
fCutV0sArmenterosAlphaK0sMin(0.2),
fCutV0sArmenterosAlphaLambdaMax(0.),
fCutV0sK0sPionNumTPCSigmaMax(5.0),
fCutV0sLambdaPionNumTPCSigmaMax(5.0),
fCutV0sLambdaProtonNumTPCSigmaMax(5.0),

// phi selection
fCutPhiMotherEtaMax(0.8),
fCutPhiInvMassMin(0.99),
fCutPhiInvMassMax(1.07),

// output lists
fQAEvents(0x0),
fQACharged(0x0),
fQAPID(0x0),
fQAV0s(0x0),
fQAPhi(0x0),
fFlowWeights(0x0),
fFlowRefs(0x0),
fFlowCharged(0x0),
fFlowPID(0x0),
fFlowPhi(0x0),
fFlowK0s(0x0),
fFlowLambda(0x0),


// flow histograms & profiles
fh3WeightsRefs(0x0),
fh3WeightsCharged(0x0),
fh3WeightsPion(0x0),
fh3WeightsKaon(0x0),
fh3WeightsProton(0x0),
fh3WeightsPhi(0x0),
fh3WeightsK0s(0x0),
fh3WeightsLambda(0x0),
fh3AfterWeightsRefs(0x0),
fh3AfterWeightsCharged(0x0),
fh3AfterWeightsPion(0x0),
fh3AfterWeightsKaon(0x0),
fh3AfterWeightsProton(0x0),
fh3AfterWeightsPhi(0x0),
fh3AfterWeightsK0s(0x0),
fh3AfterWeightsLambda(0x0),
fh2WeightRefs(0x0),
fh2WeightCharged(0x0),
fh2WeightPion(0x0),
fh2WeightKaon(0x0),
fh2WeightProton(0x0),
fh2WeightK0s(0x0),
fh2WeightLambda(0x0),
fh2WeightPhi(0x0),
fh3WeightRefs(0x0),
fh3WeightCharged(0x0),
fh3WeightPion(0x0),
fh3WeightKaon(0x0),
fh3WeightProton(0x0),
fh3WeightK0s(0x0),
fh3WeightLambda(0x0),
fh3WeightPhi(0x0),

// event histograms
fhEventSampling(0x0),
fhEventCentrality(0x0),
fh2EventCentralityNumRefs(0x0),
fhEventCounter(0x0),

// charged histogram
fhRefsMult(0x0),
fhRefsPt(0x0),
fhRefsEta(0x0),
fhRefsPhi(0x0),
fpRefsMult(0x0),
fhChargedCounter(0x0),

// PID histogram
fhPIDPionMult(0x0),
fhPIDPionPt(0x0),
fhPIDPionPhi(0x0),
fhPIDPionEta(0x0),
fhPIDPionCharge(0x0),
fhPIDKaonMult(0x0),
fhPIDKaonPt(0x0),
fhPIDKaonPhi(0x0),
fhPIDKaonEta(0x0),
fhPIDKaonCharge(0x0),
fhPIDProtonMult(0x0),
fhPIDProtonPt(0x0),
fhPIDProtonPhi(0x0),
fhPIDProtonEta(0x0),
fhPIDProtonCharge(0x0),

fh2PIDPionTPCnSigmaPion(0x0),
fh2PIDPionTOFnSigmaPion(0x0),
fh2PIDPionTPCnSigmaKaon(0x0),
fh2PIDPionTOFnSigmaKaon(0x0),
fh2PIDPionTPCnSigmaProton(0x0),
fh2PIDPionTOFnSigmaProton(0x0),
fh2PIDPionBayesPion(0x0),
fh2PIDPionBayesKaon(0x0),
fh2PIDPionBayesProton(0x0),
fh2PIDPionTPCdEdx(0x0),
fh2PIDPionTOFbeta(0x0),

fh2PIDKaonTPCdEdx(0x0),
fh2PIDKaonTOFbeta(0x0),
fh2PIDProtonTPCdEdx(0x0),
fh2PIDProtonTOFbeta(0x0),
fh2PIDKaonTPCnSigmaPion(0x0),
fh2PIDKaonTOFnSigmaPion(0x0),
fh2PIDKaonTPCnSigmaKaon(0x0),
fh2PIDKaonTOFnSigmaKaon(0x0),
fh2PIDKaonTPCnSigmaProton(0x0),
fh2PIDKaonTOFnSigmaProton(0x0),
fh2PIDKaonBayesPion(0x0),
fh2PIDKaonBayesKaon(0x0),
fh2PIDKaonBayesProton(0x0),

fh2PIDProtonTPCnSigmaPion(0x0),
fh2PIDProtonTOFnSigmaPion(0x0),
fh2PIDProtonTPCnSigmaKaon(0x0),
fh2PIDProtonTOFnSigmaKaon(0x0),
fh2PIDProtonTPCnSigmaProton(0x0),
fh2PIDProtonTOFnSigmaProton(0x0),
fh2PIDProtonBayesPion(0x0),
fh2PIDProtonBayesKaon(0x0),
fh2PIDProtonBayesProton(0x0),

fhMCRecoSelectedPionPt(0x0),
fhMCRecoSelectedTruePionPt(0x0),
fhMCRecoAllPionPt(0x0),
fhMCGenAllPionPt(0x0),
fhMCRecoSelectedKaonPt(0x0),
fhMCRecoSelectedTrueKaonPt(0x0),
fhMCRecoAllKaonPt(0x0),
fhMCGenAllKaonPt(0x0),
fhMCRecoSelectedProtonPt(0x0),
fhMCRecoSelectedTrueProtonPt(0x0),
fhMCRecoAllProtonPt(0x0),
fhMCGenAllProtonPt(0x0),

// phi histograms
fhPhiCounter(0x0),
fhPhiMult(0x0),
fhPhiBGMult(0x0),
fhPhiInvMass(0x0),
fhPhiBGInvMass(0x0),
fhPhiCharge(0x0),
fhPhiBGCharge(0x0),
fhPhiPt(0x0),
fhPhiEta(0x0),
fhPhiPhi(0x0),

// V0s histogram
fhV0sCounter(0x0),
fhV0sCounterK0s(0x0),
fhV0sCounterLambda(0x0),
fhV0sInvMassK0s(0x0),
fhV0sInvMassLambda(0x0),
fhV0sCompetingInvMassK0s(0x0),
fhV0sCompetingInvMassLambda(0x0)
{
  // Flow vectors
  for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
  {
    for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
    {
      fFlowVecQpos[iHarm][iPower] = TComplex(0,0,kFALSE);
      fFlowVecQneg[iHarm][iPower] = TComplex(0,0,kFALSE);
      fFlowVecPpos[iHarm][iPower] = TComplex(0,0,kFALSE);
      fFlowVecPneg[iHarm][iPower] = TComplex(0,0,kFALSE);
      fFlowVecS[iHarm][iPower] = TComplex(0,0,kFALSE);
    }
  }

  // Flow profiles & histograms
  for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
  {
    for(Short_t iSample(0); iSample < fNumSamples; iSample++)
    {
      fpRefsCor4[iSample][iHarm] = 0x0;
      fp2ChargedCor4[iSample][iHarm] = 0x0;
      fp2PionCor4[iSample][iHarm] = 0x0;
      fp2KaonCor4[iSample][iHarm] = 0x0;
      fp2ProtonCor4[iSample][iHarm] = 0x0;
    }

    fp3V0sCorrK0sCor4[iHarm] = 0x0;
    fp3V0sCorrLambdaCor4[iHarm] = 0x0;
    fp3PhiCorrCor4[iHarm] = 0x0;

    for(Short_t iGap(0); iGap < fNumEtaGap; iGap++)
    {
      if(iHarm == 0)
      {
        fh3V0sEntriesK0sPos[iGap] = 0x0;
        fh3V0sEntriesK0sNeg[iGap] = 0x0;
        fh3V0sEntriesLambdaPos[iGap] = 0x0;
        fh3V0sEntriesLambdaNeg[iGap] = 0x0;
        fh3PhiEntriesSignalPos[iGap] = 0x0;
        fh3PhiEntriesSignalNeg[iGap] = 0x0;
        fh3PhiEntriesBGPos[iGap] = 0x0;
        fh3PhiEntriesBGNeg[iGap] = 0x0;
      }

      for(Short_t iSample(0); iSample < fNumSamples; iSample++)
      {
        fpRefsCor2[iSample][iGap][iHarm] = 0x0;
        fp2ChargedCor2Pos[iSample][iGap][iHarm] = 0x0;
        fp2ChargedCor2Neg[iSample][iGap][iHarm] = 0x0;
        fp2PionCor2Pos[iSample][iGap][iHarm] = 0x0;
        fp2PionCor2Neg[iSample][iGap][iHarm] = 0x0;
        fp2KaonCor2Pos[iSample][iGap][iHarm] = 0x0;
        fp2KaonCor2Neg[iSample][iGap][iHarm] = 0x0;
        fp2ProtonCor2Pos[iSample][iGap][iHarm] = 0x0;
        fp2ProtonCor2Neg[iSample][iGap][iHarm] = 0x0;
      }

      fp3V0sCorrK0sCor2Pos[iGap][iHarm] = 0x0;
      fp3V0sCorrK0sCor2Neg[iGap][iHarm] = 0x0;
      fp3V0sCorrLambdaCor2Pos[iGap][iHarm] = 0x0;
      fp3V0sCorrLambdaCor2Neg[iGap][iHarm] = 0x0;
      fp3PhiCorrCor2Pos[iGap][iHarm] = 0x0;
      fp3PhiCorrCor2Neg[iGap][iHarm] = 0x0;
    }
  }

  // QA histograms
  for(Short_t iQA(0); iQA < fiNumIndexQA; iQA++)
  {
    // Event histograms
    fhQAEventsPVz[iQA] = 0x0;
    fhQAEventsNumContrPV[iQA] = 0x0;
    fhQAEventsNumSPDContrPV[iQA] = 0x0;
    fhQAEventsDistPVSPD[iQA] = 0x0;
    fhQAEventsSPDresol[iQA] = 0x0;

    // charged
    fhQAChargedMult[iQA] = 0x0;
    fhQAChargedPt[iQA] = 0x0;
    fhQAChargedEta[iQA] = 0x0;
    fhQAChargedPhi[iQA] = 0x0;
    fhQAChargedCharge[iQA] = 0x0;
    fhQAChargedFilterBit[iQA] = 0x0;
    fhQAChargedNumTPCcls[iQA] = 0x0;
    fhQAChargedDCAxy[iQA] = 0x0;
    fhQAChargedDCAz[iQA] = 0x0;

    // PID
    fhQAPIDTPCstatus[iQA] = 0x0;
    fhQAPIDTOFstatus[iQA] = 0x0;
    fhQAPIDTPCdEdx[iQA] = 0x0;
    fhQAPIDTOFbeta[iQA] = 0x0;
    fh3QAPIDnSigmaTPCTOFPtPion[iQA] = 0x0;
    fh3QAPIDnSigmaTPCTOFPtKaon[iQA] = 0x0;
    fh3QAPIDnSigmaTPCTOFPtProton[iQA] = 0x0;
    fh3QAPIDnSigmaBayesElectron[iQA] = 0x0;
    fh3QAPIDnSigmaBayesMuon[iQA] = 0x0;
    fh3QAPIDnSigmaBayesPion[iQA] = 0x0;
    fh3QAPIDnSigmaBayesKaon[iQA] = 0x0;
    fh3QAPIDnSigmaBayesProton[iQA] = 0x0;

    // V0s
    fhQAV0sMultK0s[iQA] = 0x0;
    fhQAV0sMultLambda[iQA] = 0x0;
  	fhQAV0sRecoMethod[iQA] = 0x0;
		fhQAV0sDCAtoPV[iQA] = 0x0;
		fhQAV0sDCADaughters[iQA] = 0x0;
		fhQAV0sDecayRadius[iQA] = 0x0;
    fhQAV0sDaughterTPCRefit[iQA] = 0x0;
    fhQAV0sDaughterKinks[iQA] = 0x0;
    fhQAV0sDaughterNumTPCCls[iQA] = 0x0;
    fhQAV0sDaughterNumTPCFind[iQA] = 0x0;
    fhQAV0sDaughterNumTPCCrossRows[iQA] = 0x0;
    fhQAV0sDaughterTPCCrossFindRatio[iQA] = 0x0;
    fhQAV0sDaughterNumTPCClsPID[iQA] = 0x0;
    fhQAV0sDaughterPt[iQA] = 0x0;
		fhQAV0sDaughterPhi[iQA] = 0x0;
		fhQAV0sDaughterEta[iQA] = 0x0;
    fhQAV0sDaughterCharge[iQA] = 0x0;
    fhQAV0sDaughterTPCdEdxK0s[iQA] = 0x0;
    fhQAV0sDaughterNumSigmaPionK0s[iQA] = 0x0;
    fhQAV0sDaughterTPCstatus[iQA] = 0x0;
    fhQAV0sDaughterTOFstatus[iQA] = 0x0;
    fhQAV0sDaughterTPCdEdxLambda[iQA] = 0x0;
    fhQAV0sDaughterNumSigmaPionLambda[iQA] = 0x0;
    fhQAV0sDaughterNumSigmaProtonLambda[iQA] = 0x0;
    fhQAV0sDaughterNumSigmaPionALambda[iQA] = 0x0;
    fhQAV0sDaughterNumSigmaProtonALambda[iQA] = 0x0;
    fhQAV0sMotherPt[iQA] = 0x0;
		fhQAV0sMotherPhi[iQA] = 0x0;
		fhQAV0sMotherEta[iQA] = 0x0;
    fhQAV0sMotherCharge[iQA] = 0x0;
		fhQAV0sMotherRapK0s[iQA] = 0x0;
		fhQAV0sMotherRapLambda[iQA] = 0x0;
    fhQAV0sInvMassK0s[iQA] = 0x0;
    fhQAV0sInvMassLambda[iQA] = 0x0;
		fhQAV0sCPAK0s[iQA] = 0x0;
		fhQAV0sCPALambda[iQA] = 0x0;
		fhQAV0sNumTauK0s[iQA] = 0x0;
		fhQAV0sNumTauLambda[iQA] = 0x0;
		fhQAV0sArmenterosK0s[iQA] = 0x0;
		fhQAV0sArmenterosLambda[iQA] = 0x0;
		fhQAV0sArmenterosALambda[iQA] = 0x0;
  }

  // defining input/output
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
  DefineOutput(4, TList::Class());
  DefineOutput(5, TList::Class());
  DefineOutput(6, TList::Class());
  DefineOutput(7, TList::Class());
  DefineOutput(8, TList::Class());
  DefineOutput(9, TList::Class());
  DefineOutput(10, TList::Class());
  DefineOutput(11, TList::Class());
  DefineOutput(12, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskUniFlow::~AliAnalysisTaskUniFlow()
{
  // destructor
  // if(fPIDCombined)
  // {
  //   delete fPIDCombined;
  // }

  // clearing vectors before deleting
  ClearVectors();

  // deleting FlowPart vectors (containers)
  if(fVectorRefs) delete fVectorRefs;
  if(fVectorCharged) delete fVectorCharged;
  if(fVectorPion) delete fVectorPion;
  if(fVectorKaon) delete fVectorKaon;
  if(fVectorProton) delete fVectorProton;
  if(fVectorK0s) delete fVectorK0s;
  if(fVectorLambda) delete fVectorLambda;
  if(fVectorPhi) delete fVectorPhi;

  // deleting output lists
  if(fFlowWeights) delete fFlowWeights;
  if(fFlowRefs) delete fFlowRefs;
  if(fFlowCharged) delete fFlowCharged;
  if(fFlowPID) delete fFlowPID;
  if(fFlowPhi) delete fFlowPhi;
  if(fFlowK0s) delete fFlowK0s;
  if(fFlowLambda) delete fFlowLambda;

  if(fQAEvents) delete fQAEvents;
  if(fQACharged) delete fQACharged;
  if(fQAPID) delete fQAPID;
  if(fQAPhi) delete fQAPhi;
  if(fQAV0s) delete fQAV0s;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::UserCreateOutputObjects()
{
  // create output objects
  // this function is called ONCE at the start of your analysis (RUNTIME)
  // *************************************************************

  // list all parameters used in this analysis
  ListParameters();

  // task initialization
  fInit = InitializeTask();
  if(!fInit) return;

  // creating output lists
  fFlowRefs = new TList();
  fFlowRefs->SetOwner(kTRUE);
  fFlowRefs->SetName("fFlowRefs");
  fFlowCharged = new TList();
  fFlowCharged->SetOwner(kTRUE);
  fFlowCharged->SetName("fFlowCharged");
  fFlowPID = new TList();
  fFlowPID->SetOwner(kTRUE);
  fFlowPID->SetName("fFlowPID");
  fFlowPhi = new TList();
  fFlowPhi->SetOwner(kTRUE);
  fFlowPhi->SetName("fFlowPhi");
  fFlowK0s = new TList();
  fFlowK0s->SetOwner(kTRUE);
  fFlowK0s->SetName("fFlowK0s");
  fFlowLambda = new TList();
  fFlowLambda->SetOwner(kTRUE);
  fFlowLambda->SetName("fFlowLambda");
  fFlowWeights = new TList();
  fFlowWeights->SetOwner(kTRUE);
  fFlowWeights->SetName("fFlowWeights");

  fQAEvents = new TList();
  fQAEvents->SetOwner(kTRUE);
  fQACharged = new TList();
  fQACharged->SetOwner(kTRUE);
  fQAPID = new TList();
  fQAPID->SetOwner(kTRUE);
  fQAPhi = new TList();
  fQAPhi->SetOwner(kTRUE);
  fQAV0s = new TList();
  fQAV0s->SetOwner(kTRUE);

  // setting number of bins based on set range with fixed width
  const Int_t iPOIsPtNumBins = (Int_t) (fFlowPOIsPtMax - fFlowPOIsPtMin) / 0.1 + 0.5; // fixed width 0.1 GeV/c
  const Int_t iMultNumBins = fFlowCentMax - fFlowCentMin; // fixed unit percentile or Nch bin width

  // creating histograms
    // event histogram
    fhEventSampling = new TH2D("fhEventSampling","Event sampling; centrality/multiplicity; sample index", iMultNumBins,fFlowCentMin,fFlowCentMax, fNumSamples,0,fNumSamples);
    fQAEvents->Add(fhEventSampling);
    fhEventCentrality = new TH1D("fhEventCentrality",Form("Event centrality (%s); centrality/multiplicity",fMultEstimator.Data()), iMultNumBins,fFlowCentMin,fFlowCentMax);
    fQAEvents->Add(fhEventCentrality);
    fh2EventCentralityNumRefs = new TH2D("fh2EventCentralityNumRefs",Form("Event centrality (%s) vs. N_{RFP}; centrality/multiplicity; N_{RFP}",fMultEstimator.Data()), iMultNumBins,fFlowCentMin,fFlowCentMax, 150,0,150);
    fQAEvents->Add(fh2EventCentralityNumRefs);

    const Short_t iEventCounterBins = 8;
    TString sEventCounterLabel[iEventCounterBins] = {"Input","Physics selection OK","EventCuts OK","PV OK","SPD Vtx OK","Pileup MV OK","PV #it{z} OK","Selected"};
    fhEventCounter = new TH1D("fhEventCounter","Event Counter",iEventCounterBins,0,iEventCounterBins);
    for(Short_t i(0); i < iEventCounterBins; i++) fhEventCounter->GetXaxis()->SetBinLabel(i+1, sEventCounterLabel[i].Data() );
    fQAEvents->Add(fhEventCounter);

    // flow histograms & profiles
    // weights
    if(fFlowFillWeights)
    {
      fh3WeightsRefs = new TH3D("fh3WeightsRefs","Weights: Refs; #varphi; #eta; PV-z (cm)", 100,0,TMath::TwoPi(), 80,-1.0,1.0, 2*fPVtxCutZ,-fPVtxCutZ,fPVtxCutZ);
      fh3WeightsRefs->Sumw2();
      fFlowWeights->Add(fh3WeightsRefs);
      fh3WeightsCharged = new TH3D("fh3WeightsCharged","Weights: Charged; #varphi; #eta; PV-z (cm)", 100,0,TMath::TwoPi(), 80,-1.0,1.0, 2*fPVtxCutZ,-fPVtxCutZ,fPVtxCutZ);
      fh3WeightsCharged->Sumw2();
      fFlowWeights->Add(fh3WeightsCharged);
      fh3WeightsPion = new TH3D("fh3WeightsPion","Weights: #pi; #varphi; #eta; PV-z (cm)", 100,0,TMath::TwoPi(), 80,-1.0,1.0, 2*fPVtxCutZ,-fPVtxCutZ,fPVtxCutZ);
      fh3WeightsPion->Sumw2();
      fFlowWeights->Add(fh3WeightsPion);
      fh3WeightsKaon = new TH3D("fh3WeightsKaon","Weights: K; #varphi; #eta; PV-z (cm)", 100,0,TMath::TwoPi(), 80,-1.0,1.0, 2*fPVtxCutZ,-fPVtxCutZ,fPVtxCutZ);
      fh3WeightsKaon->Sumw2();
      fFlowWeights->Add(fh3WeightsKaon);
      fh3WeightsProton = new TH3D("fh3WeightsProton","Weights: p; #varphi; #eta; PV-z (cm)", 100,0,TMath::TwoPi(), 80,-1.0,1.0, 2*fPVtxCutZ,-fPVtxCutZ,fPVtxCutZ);
      fh3WeightsProton->Sumw2();
      fFlowWeights->Add(fh3WeightsProton);
      fh3WeightsPhi = new TH3D("fh3WeightsPhi","Weights: #phi; #varphi; #eta; PV-z (cm)", 100,-TMath::Pi(),TMath::Pi(), 80,-1.0,1.0, 2*fPVtxCutZ,-fPVtxCutZ,fPVtxCutZ);
      fh3WeightsPhi->Sumw2();
      fFlowWeights->Add(fh3WeightsPhi);
      fh3WeightsK0s = new TH3D("fh3WeightsK0s","Weights: K^{0}_{S}; #varphi; #eta; PV-z (cm)", 100,0,TMath::TwoPi(), 80,-1.0,1.0, 2*fPVtxCutZ,-fPVtxCutZ,fPVtxCutZ);
      fh3WeightsK0s->Sumw2();
      fFlowWeights->Add(fh3WeightsK0s);
      fh3WeightsLambda = new TH3D("fh3WeightsLambda","Weights: #Lambda; #varphi; #eta; PV-z (cm)", 100,0,TMath::TwoPi(), 80,-1.0,1.0, 2*fPVtxCutZ,-fPVtxCutZ,fPVtxCutZ);
      fh3WeightsLambda->Sumw2();
      fFlowWeights->Add(fh3WeightsLambda);
    }

    if(fFlowUseWeights)
    {
      fh3AfterWeightsRefs = new TH3D("fh3AfterWeightsRefs","Weights: Refs; #varphi (After); #eta; PV-z (cm)", 100,0,TMath::TwoPi(), 80,-1.0,1.0, 2*fPVtxCutZ,-fPVtxCutZ,fPVtxCutZ);
      fh3AfterWeightsRefs->Sumw2();
      fFlowWeights->Add(fh3AfterWeightsRefs);
      fh3AfterWeightsCharged = new TH3D("fh3AfterWeightsCharged","Weights: Charged (After); #varphi; #eta; PV-z (cm)", 100,0,TMath::TwoPi(), 80,-1.0,1.0, 2*fPVtxCutZ,-fPVtxCutZ,fPVtxCutZ);
      fh3AfterWeightsCharged->Sumw2();
      fFlowWeights->Add(fh3AfterWeightsCharged);
      fh3AfterWeightsPion = new TH3D("fh3AfterWeightsPion","Weights: #pi (After); #varphi; #eta; PV-z (cm)", 100,0,TMath::TwoPi(), 80,-1.0,1.0, 2*fPVtxCutZ,-fPVtxCutZ,fPVtxCutZ);
      fh3AfterWeightsPion->Sumw2();
      fFlowWeights->Add(fh3AfterWeightsPion);
      fh3AfterWeightsKaon = new TH3D("fh3AfterWeightsKaon","Weights: K (After); #varphi; #eta; PV-z (cm)", 100,0,TMath::TwoPi(), 80,-1.0,1.0, 2*fPVtxCutZ,-fPVtxCutZ,fPVtxCutZ);
      fh3AfterWeightsKaon->Sumw2();
      fFlowWeights->Add(fh3AfterWeightsKaon);
      fh3AfterWeightsProton = new TH3D("fh3AfterWeightsProton","Weights: p (After); #varphi; #eta; PV-z (cm)", 100,0,TMath::TwoPi(), 80,-1.0,1.0, 2*fPVtxCutZ,-fPVtxCutZ,fPVtxCutZ);
      fh3AfterWeightsProton->Sumw2();
      fFlowWeights->Add(fh3AfterWeightsProton);
      fh3AfterWeightsPhi = new TH3D("fh3AfterWeightsPhi","Weights: #phi (After); #varphi; #eta; PV-z (cm)", 100,-TMath::Pi(),TMath::Pi(), 80,-1.0,1.0, 2*fPVtxCutZ,-fPVtxCutZ,fPVtxCutZ);
      fh3AfterWeightsPhi->Sumw2();
      fFlowWeights->Add(fh3AfterWeightsPhi);
      fh3AfterWeightsK0s = new TH3D("fh3AfterWeightsK0s","Weights: K^{0}_{S} (After); #varphi; #eta; PV-z (cm)", 100,0,TMath::TwoPi(), 80,-1.0,1.0, 2*fPVtxCutZ,-fPVtxCutZ,fPVtxCutZ);
      fh3AfterWeightsK0s->Sumw2();
      fFlowWeights->Add(fh3AfterWeightsK0s);
      fh3AfterWeightsLambda = new TH3D("fh3AfterWeightsLambda","Weights: #Lambda (After); #varphi; #eta; PV-z (cm)", 100,0,TMath::TwoPi(), 80,-1.0,1.0, 2*fPVtxCutZ,-fPVtxCutZ,fPVtxCutZ);
      fh3AfterWeightsLambda->Sumw2();
      fFlowWeights->Add(fh3AfterWeightsLambda);
    }

    // candidate distribution for flow-mass method
    for(Short_t iGap(0); iGap < fNumEtaGap; iGap++)
    {
      if(fProcessPhi)
      {
        fh3PhiEntriesSignalPos[iGap] = new TH3D(Form("fh3PhiEntriesSignal_gap%02.2g_Pos",10*fEtaGap[iGap]), Form("#phi: Distribution (Gap %g) | POIs pos); centrality/multiplicity; #it{p}_{T} (GeV/c); #it{m}_{inv} (GeV/#it{c}^{2})",fEtaGap[iGap]), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fPhiNumBinsMass,fCutPhiInvMassMin,fCutPhiInvMassMax);
        fh3PhiEntriesSignalPos[iGap]->Sumw2();
        fFlowPhi->Add(fh3PhiEntriesSignalPos[iGap]);
        fh3PhiEntriesBGPos[iGap] = new TH3D(Form("fh3PhiEntriesBG_gap%02.2g_Pos",10*fEtaGap[iGap]), Form("#phi (BG): Distribution (Gap %g | POIs pos); centrality/multiplicity; #it{p}_{T} (GeV/c); #it{m}_{inv} (GeV/#it{c}^{2})",fEtaGap[iGap]), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fPhiNumBinsMass,fCutPhiInvMassMin,fCutPhiInvMassMax);
        fh3PhiEntriesBGPos[iGap]->Sumw2();
        fFlowPhi->Add(fh3PhiEntriesBGPos[iGap]);

        if(fEtaGap[iGap] != -1.)
        {
          fh3PhiEntriesSignalNeg[iGap] = new TH3D(Form("fh3PhiEntriesSignal_gap%02.2g_Neg",10*fEtaGap[iGap]), Form("#phi: Distribution (Gap %g) | POIs neg); centrality/multiplicity; #it{p}_{T} (GeV/c); #it{m}_{inv} (GeV/#it{c}^{2})",fEtaGap[iGap]), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fPhiNumBinsMass,fCutPhiInvMassMin,fCutPhiInvMassMax);
          fh3PhiEntriesSignalNeg[iGap]->Sumw2();
          fFlowPhi->Add(fh3PhiEntriesSignalNeg[iGap]);
          fh3PhiEntriesBGNeg[iGap] = new TH3D(Form("fh3PhiEntriesBG_gap%02.2g_Neg",10*fEtaGap[iGap]), Form("#phi (BG): Distribution (Gap %g | POIs neg); centrality/multiplicity; #it{p}_{T} (GeV/c); #it{m}_{inv} (GeV/#it{c}^{2})",fEtaGap[iGap]), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fPhiNumBinsMass,fCutPhiInvMassMin,fCutPhiInvMassMax);
          fh3PhiEntriesBGNeg[iGap]->Sumw2();
          fFlowPhi->Add(fh3PhiEntriesBGNeg[iGap]);
        }
      }

      if(fProcessV0s)
      {
        fh3V0sEntriesK0sPos[iGap] = new TH3D(Form("fh3V0sEntriesK0s_gap%02.2g_Pos",10*fEtaGap[iGap]), Form("K_{S}^{0}: Distribution (Gap %g | POIs pos); centrality/multiplicity; #it{p}_{T} (GeV/c); #it{m}_{inv} (GeV/#it{c}^{2})",fEtaGap[iGap]), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fV0sNumBinsMass,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax);
        fh3V0sEntriesK0sPos[iGap]->Sumw2();
        fFlowK0s->Add(fh3V0sEntriesK0sPos[iGap]);
        fh3V0sEntriesLambdaPos[iGap] = new TH3D(Form("fh3V0sEntriesLambda_gap%02.2g_Pos",10*fEtaGap[iGap]), Form("#Lambda/#bar{#Lambda}: Distribution (Gap %g | POIs pos); centrality/multiplicity; #it{p}_{T} (GeV/c); #it{m}_{inv} (GeV/#it{c}^{2})",fEtaGap[iGap]), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fV0sNumBinsMass,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax);
        fh3V0sEntriesLambdaPos[iGap]->Sumw2();
        fFlowLambda->Add(fh3V0sEntriesLambdaPos[iGap]);

        if(fEtaGap[iGap] != -1.)
        {
          fh3V0sEntriesK0sNeg[iGap] = new TH3D(Form("fh3V0sEntriesK0s_gap%02.2g_Neg",10*fEtaGap[iGap]), Form("K_{S}^{0}: Distribution (Gap %g | POIs neg); centrality/multiplicity; #it{p}_{T} (GeV/c); #it{m}_{inv} (GeV/#it{c}^{2})",fEtaGap[iGap]), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fV0sNumBinsMass,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax);
          fh3V0sEntriesK0sNeg[iGap]->Sumw2();
          fFlowK0s->Add(fh3V0sEntriesK0sNeg[iGap]);
          fh3V0sEntriesLambdaNeg[iGap] = new TH3D(Form("fh3V0sEntriesLambda_gap%02.2g_Neg",10*fEtaGap[iGap]), Form("#Lambda/#bar{#Lambda}: Distribution (Gap %g | POIs neg); centrality/multiplicity; #it{p}_{T} (GeV/c); #it{m}_{inv} (GeV/#it{c}^{2})",fEtaGap[iGap]), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fV0sNumBinsMass,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax);
          fh3V0sEntriesLambdaNeg[iGap]->Sumw2();
          fFlowLambda->Add(fh3V0sEntriesLambdaNeg[iGap]);
        }
      }
    }

    // correlations
    if(fRunMode != kSkipFlow)
    {
      for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
      {
        for(Short_t iGap(0); iGap < fNumEtaGap; iGap++)
        {
          for(Short_t iSample(0); iSample < fNumSamples; iSample++)
          {
            if(!fSampling && iSample > 0) break; // define only one sample histogram if sampling is off

            fpRefsCor2[iSample][iGap][iHarm] = new TProfile(Form("fpRefs_<2>_harm%d_gap%02.2g_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("Ref: <<2>> | Gap %g | n=%d | sample %d ; centrality/multiplicity;",fEtaGap[iGap],fHarmonics[iHarm],iSample), iMultNumBins,fFlowCentMin,fFlowCentMax);
            fpRefsCor2[iSample][iGap][iHarm]->Sumw2(kTRUE);
            fFlowRefs->Add(fpRefsCor2[iSample][iGap][iHarm]);

            if(fCutFlowDoFourCorrelations && iGap == 0)
            {
              fpRefsCor4[iSample][iHarm] = new TProfile(Form("fpRefs_<4>_harm%d_gap%02.2g_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("Ref: <<4>> | Gap %g | n=%d | sample %d ; centrality/multiplicity;",fEtaGap[iGap],fHarmonics[iHarm],iSample), iMultNumBins,fFlowCentMin,fFlowCentMax);
              fpRefsCor4[iSample][iHarm]->Sumw2(kTRUE);
              fFlowRefs->Add(fpRefsCor4[iSample][iHarm]);
            }

            fp2ChargedCor2Pos[iSample][iGap][iHarm] = new TProfile2D(Form("fp2Charged_<2>_harm%d_gap%02.2g_Pos_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("Charged: <<2'>> | Gap %g | n=%d | sample %d | POIs pos; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fHarmonics[iHarm],iSample), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
            fp2ChargedCor2Pos[iSample][iGap][iHarm]->Sumw2(kTRUE);
            fFlowCharged->Add(fp2ChargedCor2Pos[iSample][iGap][iHarm]);

            if(fEtaGap[iGap] != -1.)
            {
              fp2ChargedCor2Neg[iSample][iGap][iHarm] = new TProfile2D(Form("fp2Charged_<2>_harm%d_gap%02.2g_Neg_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("Charged: <<2'>> | Gap %g | n=%d | sample %d | POIs neg; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fHarmonics[iHarm],iSample), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
              fp2ChargedCor2Neg[iSample][iGap][iHarm]->Sumw2(kTRUE);
              fFlowCharged->Add(fp2ChargedCor2Neg[iSample][iGap][iHarm]);
            }

            if(fCutFlowDoFourCorrelations && iGap == 0)
            {
              fp2ChargedCor4[iSample][iHarm] = new TProfile2D(Form("fp2Charged_<4>_harm%d_gap%02.2g_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("Charged: <<4'>> | Gap %g | n=%d | sample %d; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fHarmonics[iHarm],iSample), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
              fp2ChargedCor4[iSample][iHarm]->Sumw2(kTRUE);
              fFlowCharged->Add(fp2ChargedCor4[iSample][iHarm]);
            }


            if(fProcessPID)
            {
              fp2PionCor2Pos[iSample][iGap][iHarm] = new TProfile2D(Form("fp2Pion_<2>_harm%d_gap%02.2g_Pos_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("PID #pi: <<2'>> | Gap %g | n=%d | sample %d  | POIs pos; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fHarmonics[iHarm],iSample), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
              fp2PionCor2Pos[iSample][iGap][iHarm]->Sumw2(kTRUE);
              fFlowPID->Add(fp2PionCor2Pos[iSample][iGap][iHarm]);

              fp2KaonCor2Pos[iSample][iGap][iHarm] = new TProfile2D(Form("fp2Kaon_<2>_harm%d_gap%02.2g_Pos_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("PID K: <<2'>> | Gap %g | n=%d | sample %d | POIs pos; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fHarmonics[iHarm],iSample), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
              fp2KaonCor2Pos[iSample][iGap][iHarm]->Sumw2(kTRUE);
              fFlowPID->Add(fp2KaonCor2Pos[iSample][iGap][iHarm]);

              fp2ProtonCor2Pos[iSample][iGap][iHarm] = new TProfile2D(Form("fp2Proton_<2>_harm%d_gap%02.2g_Pos_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("PID p: <<2'>> | Gap %g | n=%d | sample %d | POIs pos; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fHarmonics[iHarm],iSample), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
              fp2ProtonCor2Pos[iSample][iGap][iHarm]->Sumw2(kTRUE);
              fFlowPID->Add(fp2ProtonCor2Pos[iSample][iGap][iHarm]);

              if(fEtaGap[iGap] != -1.)
              {
                fp2PionCor2Neg[iSample][iGap][iHarm] = new TProfile2D(Form("fp2Pion_<2>_harm%d_gap%02.2g_Neg_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("PID #pi: <<2'>> | Gap %g | n=%d | sample %d | POIs neg; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fHarmonics[iHarm],iSample), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
                fp2PionCor2Neg[iSample][iGap][iHarm]->Sumw2(kTRUE);
                fFlowPID->Add(fp2PionCor2Neg[iSample][iGap][iHarm]);

                fp2KaonCor2Neg[iSample][iGap][iHarm] = new TProfile2D(Form("fp2Kaon_<2>_harm%d_gap%02.2g_Neg_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("PID K: <<2'>> | Gap %g | n=%d | sample %d | POIs neg; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fHarmonics[iHarm],iSample), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
                fp2KaonCor2Neg[iSample][iGap][iHarm]->Sumw2(kTRUE);
                fFlowPID->Add(fp2KaonCor2Neg[iSample][iGap][iHarm]);

                fp2ProtonCor2Neg[iSample][iGap][iHarm] = new TProfile2D(Form("fp2Proton_<2>_harm%d_gap%02.2g_Neg_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("PID p: <<2'>> | Gap %g | n=%d | sample %d | POIs neg; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fHarmonics[iHarm],iSample), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
                fp2ProtonCor2Neg[iSample][iGap][iHarm]->Sumw2(kTRUE);
                fFlowPID->Add(fp2ProtonCor2Neg[iSample][iGap][iHarm]);
              }

              if(fCutFlowDoFourCorrelations && iGap == 0)
              {
                fp2PionCor4[iSample][iHarm] = new TProfile2D(Form("fp2Pion_<4>_harm%d_gap%02.2g_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("PID #pi: <<4'>> | Gap %g | n=%d | sample %d; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fHarmonics[iHarm],iSample), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
                fp2PionCor4[iSample][iHarm]->Sumw2(kTRUE);
                fFlowPID->Add(fp2PionCor4[iSample][iHarm]);
                fp2KaonCor4[iSample][iHarm] = new TProfile2D(Form("fp2Kaon_<4>_harm%d_gap%02.2g_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("PID K: <<4'>> | Gap %g | n=%d | sample %d; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fHarmonics[iHarm],iSample), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
                fp2KaonCor4[iSample][iHarm]->Sumw2(kTRUE);
                fFlowPID->Add(fp2KaonCor4[iSample][iHarm]);
                fp2ProtonCor4[iSample][iHarm] = new TProfile2D(Form("fp2Proton_<4>_harm%d_gap%02.2g_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("PID p: <<4'>> | Gap %g | n=%d | sample %d; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fHarmonics[iHarm],iSample), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
                fp2ProtonCor4[iSample][iHarm]->Sumw2(kTRUE);
                fFlowPID->Add(fp2ProtonCor4[iSample][iHarm]);
              }
            }
          }

          if(fProcessPhi)
          {
            fp3PhiCorrCor2Pos[iGap][iHarm] = new TProfile3D(Form("fp3PhiCorr_<2>_harm%d_gap%02.2g_Pos",fHarmonics[iHarm],10*fEtaGap[iGap]), Form("#phi: <<2'>> | Gap %g | n=%d | POIs pos; centrality/multiplicity; #it{p}_{T} (GeV/c); #it{m}_{inv} (GeV/#it{c}^{2})",fEtaGap[iGap],fHarmonics[iHarm]), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fPhiNumBinsMass,fCutPhiInvMassMin,fCutPhiInvMassMax);
            fp3PhiCorrCor2Pos[iGap][iHarm]->Sumw2();
            fFlowPhi->Add(fp3PhiCorrCor2Pos[iGap][iHarm]);

            if(fEtaGap[iGap] != -1.)
            {
              fp3PhiCorrCor2Neg[iGap][iHarm] = new TProfile3D(Form("fp3PhiCorr_<2>_harm%d_gap%02.2g_Neg",fHarmonics[iHarm],10*fEtaGap[iGap]), Form("#phi: <<2'>> | Gap %g | n=%d  | POIs neg; centrality/multiplicity; #it{p}_{T} (GeV/c); #it{m}_{inv} (GeV/#it{c}^{2})",fEtaGap[iGap],fHarmonics[iHarm]), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fPhiNumBinsMass,fCutPhiInvMassMin,fCutPhiInvMassMax);
              fp3PhiCorrCor2Neg[iGap][iHarm]->Sumw2();
              fFlowPhi->Add(fp3PhiCorrCor2Neg[iGap][iHarm]);
            }

            if(fCutFlowDoFourCorrelations && iGap == 0)
            {
              fp3PhiCorrCor4[iHarm] = new TProfile3D(Form("fp3PhiCorr_<4>_harm%d_gap%02.2g",fHarmonics[iHarm],10*fEtaGap[iGap]), Form("#phi: <<4'>> | Gap %g | n=%d; centrality/multiplicity; #it{p}_{T} (GeV/c); #it{m}_{inv} (GeV/#it{c}^{2})",fEtaGap[iGap],fHarmonics[iHarm]), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fPhiNumBinsMass,fCutPhiInvMassMin,fCutPhiInvMassMax);
              fp3PhiCorrCor4[iHarm]->Sumw2();
              fFlowPhi->Add(fp3PhiCorrCor4[iHarm]);
            }
          }

          if(fProcessV0s)
          {
            fp3V0sCorrK0sCor2Pos[iGap][iHarm] = new TProfile3D(Form("fp3V0sCorrK0s_<2>_harm%d_gap%02.2g_Pos",fHarmonics[iHarm],10*fEtaGap[iGap]), Form("K_{S}^{0}: <<2'>> | Gap %g | n=%d | POIs pos; centrality/multiplicity; #it{p}_{T} (GeV/c); #it{m}_{inv} (GeV/#it{c}^{2})",fEtaGap[iGap],fHarmonics[iHarm]), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fV0sNumBinsMass,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax);
            fp3V0sCorrK0sCor2Pos[iGap][iHarm]->Sumw2();
            fFlowK0s->Add(fp3V0sCorrK0sCor2Pos[iGap][iHarm]);
            fp3V0sCorrLambdaCor2Pos[iGap][iHarm] = new TProfile3D(Form("fp3V0sCorrLambda_<2>_harm%d_gap%02.2g_Pos",fHarmonics[iHarm],10*fEtaGap[iGap]), Form("#Lambda/#bar{#Lambda}: <<2'>> | Gap %g | n=%d | POIs pos; centrality/multiplicity; #it{p}_{T} (GeV/c); #it{m}_{inv} (GeV/#it{c}^{2})",fEtaGap[iGap],fHarmonics[iHarm]), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fV0sNumBinsMass,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax);
            fp3V0sCorrLambdaCor2Pos[iGap][iHarm]->Sumw2();
            fFlowLambda->Add(fp3V0sCorrLambdaCor2Pos[iGap][iHarm]);

            if(fEtaGap[iGap] != -1.)
            {
              fp3V0sCorrK0sCor2Neg[iGap][iHarm] = new TProfile3D(Form("fp3V0sCorrK0s_<2>_harm%d_gap%02.2g_Neg",fHarmonics[iHarm],10*fEtaGap[iGap]), Form("K_{S}^{0}: <<2'>> | Gap %g | n=%d | POIs neg; centrality/multiplicity; #it{p}_{T} (GeV/c); #it{m}_{inv} (GeV/#it{c}^{2})",fEtaGap[iGap],fHarmonics[iHarm]), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fV0sNumBinsMass,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax);
              fp3V0sCorrK0sCor2Neg[iGap][iHarm]->Sumw2();
              fFlowK0s->Add(fp3V0sCorrK0sCor2Neg[iGap][iHarm]);

              fp3V0sCorrLambdaCor2Neg[iGap][iHarm] = new TProfile3D(Form("fp3V0sCorrLambda_<2>_harm%d_gap%02.2g_Neg",fHarmonics[iHarm],10*fEtaGap[iGap]), Form("#Lambda/#bar{#Lambda}: <<2'>> | Gap %g | n=%d | POIs neg; centrality/multiplicity; #it{p}_{T} (GeV/c); #it{m}_{inv} (GeV/#it{c}^{2})",fEtaGap[iGap],fHarmonics[iHarm]), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fV0sNumBinsMass,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax);
              fp3V0sCorrLambdaCor2Neg[iGap][iHarm]->Sumw2();
              fFlowLambda->Add(fp3V0sCorrLambdaCor2Neg[iGap][iHarm]);
            }

            if(fCutFlowDoFourCorrelations && iGap == 0)
            {
              fp3V0sCorrK0sCor4[iHarm] = new TProfile3D(Form("fp3V0sCorrK0s_<4>_harm%d_gap%02.2g",fHarmonics[iHarm],10*fEtaGap[iGap]), Form("K_{S}^{0}: <<4'>> | Gap %g | n=%d; centrality/multiplicity; #it{p}_{T} (GeV/c); #it{m}_{inv} (GeV/#it{c}^{2})",fEtaGap[iGap],fHarmonics[iHarm]), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fV0sNumBinsMass,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax);
              fp3V0sCorrK0sCor4[iHarm]->Sumw2();
              fFlowK0s->Add(fp3V0sCorrK0sCor4[iHarm]);
              fp3V0sCorrLambdaCor4[iHarm] = new TProfile3D(Form("fp3V0sCorrLambda_<4>_harm%d_gap%02.2g",fHarmonics[iHarm],10*fEtaGap[iGap]), Form("#Lambda/#bar{#Lambda}: <<4'>> | Gap %g | n=%d; centrality/multiplicity; #it{p}_{T} (GeV/c); #it{m}_{inv} (GeV/#it{c}^{2})",fEtaGap[iGap],fHarmonics[iHarm]), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fV0sNumBinsMass,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax);
              fp3V0sCorrLambdaCor4[iHarm]->Sumw2();
              fFlowLambda->Add(fp3V0sCorrLambdaCor4[iHarm]);
            }
          }
        }
      }
    }

    // charged (tracks) histograms
    fhRefsMult = new TH1D("fhRefsMult","RFPs: Multiplicity; multiplicity", 1000,0,1000);
    fQACharged->Add(fhRefsMult);
    fhRefsPt = new TH1D("fhRefsPt","RFPs: #it{p}_{T};  #it{p}_{T} (GeV/#it{c})", 300,0,30);
    fQACharged->Add(fhRefsPt);
    fhRefsEta = new TH1D("fhRefsEta","RFPs: #eta; #eta", 151,-1.5,1.5);
    fQACharged->Add(fhRefsEta);
    fhRefsPhi = new TH1D("fhRefsPhi","RFPs: #varphi; #varphi", 100,0,TMath::TwoPi());
    fQACharged->Add(fhRefsPhi);
    fpRefsMult = new TProfile("fpRefsMult","Ref mult; centrality/multiplicity", iMultNumBins,fFlowCentMin,fFlowCentMax);
    fpRefsMult->Sumw2();
    fQACharged->Add(fpRefsMult);

    TString sChargedCounterLabel[] = {"Input","Pt","Eta","FB","#TPC-Cls","DCA-z","DCA-xy","Selected"};
    const Short_t iNBinsChargedCounter = sizeof(sChargedCounterLabel)/sizeof(sChargedCounterLabel[0]);
    fhChargedCounter = new TH1D("fhChargedCounter","Charged tracks: Counter",iNBinsChargedCounter,0,iNBinsChargedCounter);
    for(Short_t i(0); i < iNBinsChargedCounter; i++) fhChargedCounter->GetXaxis()->SetBinLabel(i+1, sChargedCounterLabel[i].Data() );
    fQACharged->Add(fhChargedCounter);

    // PID tracks histograms
    if(fProcessPID || fProcessPhi)
    {
      fhPIDPionMult = new TH1D("fhPIDPionMult","PID: #pi: Multiplicity; multiplicity", 200,0,200);
      fQAPID->Add(fhPIDPionMult);
      fhPIDKaonMult = new TH1D("fhPIDKaonMult","PID: K: Multiplicity; multiplicity", 100,0,100);
      fQAPID->Add(fhPIDKaonMult);
      fhPIDProtonMult = new TH1D("fhPIDProtonMult","PID: p: Multiplicity; multiplicity", 100,0,100);
      fQAPID->Add(fhPIDProtonMult);

      if(fMC)
      {
        fhMCRecoSelectedPionPt = new TH1D("fhMCRecoSelectedPionPt","fhMCRecoSelectedPionPt; p_{T} (GeV/c); Counts", iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fQAPID->Add(fhMCRecoSelectedPionPt);
        fhMCRecoSelectedTruePionPt = new TH1D("fhMCRecoSelectedTruePionPt","fhMCRecoSelectedTruePionPt; p_{T} (GeV/c); Counts", iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fQAPID->Add(fhMCRecoSelectedTruePionPt);
        fhMCRecoAllPionPt = new TH1D("fhMCRecoAllPionPt","fhMCRecoAllPionPt; p_{T} (GeV/c); Counts", iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fQAPID->Add(fhMCRecoAllPionPt);
        fhMCGenAllPionPt = new TH1D("fhMCGenAllPionPt","fhMCGenAllPionPt; p_{T} (GeV/c); Counts", iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fQAPID->Add(fhMCGenAllPionPt);

        fhMCRecoSelectedKaonPt = new TH1D("fhMCRecoSelectedKaonPt","fhMCRecoSelectedKaonPt; p_{T} (GeV/c); Counts", iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fQAPID->Add(fhMCRecoSelectedKaonPt);
        fhMCRecoSelectedTrueKaonPt = new TH1D("fhMCRecoSelectedTrueKaonPt","fhMCRecoSelectedTrueKaonPt; p_{T} (GeV/c); Counts", iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fQAPID->Add(fhMCRecoSelectedTrueKaonPt);
        fhMCRecoAllKaonPt = new TH1D("fhMCRecoAllKaonPt","fhMCRecoAllKaonPt; p_{T} (GeV/c); Counts", iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fQAPID->Add(fhMCRecoAllKaonPt);
        fhMCGenAllKaonPt = new TH1D("fhMCGenAllKaonPt","fhMCGenAllKaonPt; p_{T} (GeV/c); Counts", iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fQAPID->Add(fhMCGenAllKaonPt);

        fhMCRecoSelectedProtonPt = new TH1D("fhMCRecoSelectedProtonPt","fhMCRecoSelectedProtonPt; p_{T} (GeV/c); Counts", iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fQAPID->Add(fhMCRecoSelectedProtonPt);
        fhMCRecoSelectedTrueProtonPt = new TH1D("fhMCRecoSelectedTrueProtonPt","fhMCRecoSelectedTrueProtonPt; p_{T} (GeV/c); Counts", iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fQAPID->Add(fhMCRecoSelectedTrueProtonPt);
        fhMCRecoAllProtonPt = new TH1D("fhMCRecoAllProtonPt","fhMCRecoAllProtonPt; p_{T} (GeV/c); Counts", iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fQAPID->Add(fhMCRecoAllProtonPt);
        fhMCGenAllProtonPt = new TH1D("fhMCGenAllProtonPt","fhMCGenAllProtonPt; p_{T} (GeV/c); Counts", iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fQAPID->Add(fhMCGenAllProtonPt);
      }


      if(fFillQA)
      {
        fhPIDPionPt = new TH1D("fhPIDPionPt","PID: #pi: #it{p}_{T}; #it{p}_{T}", 150,0.,30.);
        fQAPID->Add(fhPIDPionPt);
        fhPIDPionPhi = new TH1D("fhPIDPionPhi","PID: #pi: #varphi; #varphi", 100,0,TMath::TwoPi());
        fQAPID->Add(fhPIDPionPhi);
        fhPIDPionEta = new TH1D("fhPIDPionEta","PID: #pi: #eta; #eta", 151,-1.5,1.5);
        fQAPID->Add(fhPIDPionEta);
        fhPIDPionCharge = new TH1D("fhPIDPionCharge","PID: #pi: charge; charge", 3,-1.5,1.5);
        fQAPID->Add(fhPIDPionCharge);
        fhPIDKaonPt = new TH1D("fhPIDKaonPt","PID: K: #it{p}_{T}; #it{p}_{T}", 150,0.,30.);
        fQAPID->Add(fhPIDKaonPt);
        fhPIDKaonPhi = new TH1D("fhPIDKaonPhi","PID: K: #varphi; #varphi", 100,0,TMath::TwoPi());
        fQAPID->Add(fhPIDKaonPhi);
        fhPIDKaonEta = new TH1D("fhPIDKaonEta","PID: K: #eta; #eta", 151,-1.5,1.5);
        fQAPID->Add(fhPIDKaonEta);
        fhPIDKaonCharge = new TH1D("fhPIDKaonCharge","PID: K: charge; charge", 3,-1.5,1.5);
        fQAPID->Add(fhPIDKaonCharge);
        fhPIDProtonPt = new TH1D("fhPIDProtonPt","PID: p: #it{p}_{T}; #it{p}_{T}", 150,0.,30.);
        fQAPID->Add(fhPIDProtonPt);
        fhPIDProtonPhi = new TH1D("fhPIDProtonPhi","PID: p: #varphi; #varphi", 100,0,TMath::TwoPi());
        fQAPID->Add(fhPIDProtonPhi);
        fhPIDProtonEta = new TH1D("fhPIDProtonEta","PID: p: #eta; #eta", 151,-1.5,1.5);
        fQAPID->Add(fhPIDProtonEta);
        fhPIDProtonCharge = new TH1D("fhPIDProtonCharge","PID: p: charge; charge", 3,-1.5,1.5);
        fQAPID->Add(fhPIDProtonCharge);
        fh2PIDPionTPCdEdx = new TH2D("fh2PIDPionTPCdEdx","PID: #pi: TPC dE/dx; #it{p} (GeV/#it{c}); TPC dE/dx", 200,0,20, 131,-10,1000);
        fQAPID->Add(fh2PIDPionTPCdEdx);
        fh2PIDPionTOFbeta = new TH2D("fh2PIDPionTOFbeta","PID: #pi: TOF #beta; #it{p} (GeV/#it{c});TOF #beta", 200,0,20, 101,-0.1,1.5);
        fQAPID->Add(fh2PIDPionTOFbeta);
        fh2PIDKaonTPCdEdx = new TH2D("fh2PIDKaonTPCdEdx","PID: K: TPC dE/dx; #it{p} (GeV/#it{c}); TPC dE/dx", 200,0,20, 131,-10,1000);
        fQAPID->Add(fh2PIDKaonTPCdEdx);
        fh2PIDKaonTOFbeta = new TH2D("fh2PIDKaonTOFbeta","PID: K: TOF #beta; #it{p} (GeV/#it{c});TOF #beta", 200,0,20, 101,-0.1,1.5);
        fQAPID->Add(fh2PIDKaonTOFbeta);
        fh2PIDProtonTPCdEdx = new TH2D("fh2PIDProtonTPCdEdx","PID: p: TPC dE/dx; #it{p} (GeV/#it{c}); TPC dE/dx", 200,0,20, 131,-10,1000);
        fQAPID->Add(fh2PIDProtonTPCdEdx);
        fh2PIDProtonTOFbeta = new TH2D("fh2PIDProtonTOFbeta","PID: p: TOF #beta; #it{p} (GeV/#it{c});TOF #beta", 200,0,20, 101,-0.1,1.5);
        fQAPID->Add(fh2PIDProtonTOFbeta);
        fh2PIDPionTPCnSigmaPion = new TH2D("fh2PIDPionTPCnSigmaPion","PID: #pi: TPC n#sigma (#pi hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDPionTPCnSigmaPion);
        fh2PIDPionTOFnSigmaPion = new TH2D("fh2PIDPionTOFnSigmaPion","PID: #pi: TOF n#sigma (#pi hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDPionTOFnSigmaPion);
        fh2PIDPionTPCnSigmaKaon = new TH2D("fh2PIDPionTPCnSigmaKaon","PID: #pi: TPC n#sigma (K hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDPionTPCnSigmaKaon);
        fh2PIDPionTOFnSigmaKaon = new TH2D("fh2PIDPionTOFnSigmaKaon","PID: #pi: TOF n#sigma (K hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDPionTOFnSigmaKaon);
        fh2PIDPionTPCnSigmaProton = new TH2D("fh2PIDPionTPCnSigmaProton","PID: #pi: TPC n#sigma (p hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDPionTPCnSigmaProton);
        fh2PIDPionTOFnSigmaProton = new TH2D("fh2PIDPionTOFnSigmaProton","PID: #pi: TOF n#sigma (p hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDPionTOFnSigmaProton);
        fh2PIDPionBayesPion = new TH2D("fh2PIDPionBayesPion","PID: #pi: Bayes probability (#pi hyp.); #it{p}_{T} (GeV/#it{c}); Bayes prob.", 200,0,20, 50,0,1);
        fQAPID->Add(fh2PIDPionBayesPion);
        fh2PIDPionBayesKaon = new TH2D("fh2PIDPionBayesKaon","PID: #pi: Bayes probability (K hyp.); #it{p}_{T} (GeV/#it{c}); Bayes prob.", 200,0,20, 50,0,1);
        fQAPID->Add(fh2PIDPionBayesKaon);
        fh2PIDPionBayesProton = new TH2D("fh2PIDPionBayesProton","PID: #pi: Bayes probability (p hyp.); #it{p}_{T} (GeV/#it{c}); Bayes prob.", 200,0,20, 50,0,1);
        fQAPID->Add(fh2PIDPionBayesProton);
        fh2PIDKaonTPCnSigmaPion = new TH2D("fh2PIDKaonTPCnSigmaPion","PID: K: TPC n#sigma (#pi hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDKaonTPCnSigmaPion);
        fh2PIDKaonTOFnSigmaPion = new TH2D("fh2PIDKaonTOFnSigmaPion","PID: K: TOF n#sigma (#pi hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDKaonTOFnSigmaPion);
        fh2PIDKaonTPCnSigmaKaon = new TH2D("fh2PIDKaonTPCnSigmaKaon","PID: K: TPC n#sigma (K hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDKaonTPCnSigmaKaon);
        fh2PIDKaonTOFnSigmaKaon = new TH2D("fh2PIDKaonTOFnSigmaKaon","PID: K: TOF n#sigma (K hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDKaonTOFnSigmaKaon);
        fh2PIDKaonTPCnSigmaProton = new TH2D("fh2PIDKaonTPCnSigmaProton","PID: K: TPC n#sigma (p hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDKaonTPCnSigmaProton);
        fh2PIDKaonTOFnSigmaProton = new TH2D("fh2PIDKaonTOFnSigmaProton","PID: K: TOF n#sigma (p hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDKaonTOFnSigmaProton);
        fh2PIDKaonBayesPion = new TH2D("fh2PIDKaonBayesPion","PID: K: Bayes probability (#pi hyp.); #it{p}_{T} (GeV/#it{c}); Bayes prob.", 200,0,20, 50,0,1);
        fQAPID->Add(fh2PIDKaonBayesPion);
        fh2PIDKaonBayesKaon = new TH2D("fh2PIDKaonBayesKaon","PID: K: Bayes probability (K hyp.); #it{p}_{T} (GeV/#it{c}); Bayes prob.", 200,0,20, 50,0,1);
        fQAPID->Add(fh2PIDKaonBayesKaon);
        fh2PIDKaonBayesProton = new TH2D("fh2PIDKaonBayesProton","PID: K: Bayes probability (p hyp.); #it{p}_{T} (GeV/#it{c}); Bayes prob.", 200,0,20, 50,0,1);
        fQAPID->Add(fh2PIDKaonBayesProton);
        fh2PIDProtonTPCnSigmaPion = new TH2D("fh2PIDProtonTPCnSigmaPion","PID: p: TPC n#sigma (#pi hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDProtonTPCnSigmaPion);
        fh2PIDProtonTOFnSigmaPion = new TH2D("fh2PIDProtonTOFnSigmaPion","PID: p: TOF n#sigma (#pi hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDProtonTOFnSigmaPion);
        fh2PIDProtonTPCnSigmaKaon = new TH2D("fh2PIDProtonTPCnSigmaKaon","PID: p: TPC n#sigma (K hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDProtonTPCnSigmaKaon);
        fh2PIDProtonTOFnSigmaKaon = new TH2D("fh2PIDProtonTOFnSigmaKaon","PID: p: TOF n#sigma (K hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDProtonTOFnSigmaKaon);
        fh2PIDProtonTPCnSigmaProton = new TH2D("fh2PIDProtonTPCnSigmaProton","PID: p: TPC n#sigma (p hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDProtonTPCnSigmaProton);
        fh2PIDProtonTOFnSigmaProton = new TH2D("fh2PIDProtonTOFnSigmaProton","PID: p: TOF n#sigma (p hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDProtonTOFnSigmaProton);
        fh2PIDProtonBayesPion = new TH2D("fh2PIDProtonBayesPion","PID: p: Bayes probability (#pi hyp.); #it{p}_{T} (GeV/#it{c}); Bayes prob.", 200,0,20, 50,0,1);
        fQAPID->Add(fh2PIDProtonBayesPion);
        fh2PIDProtonBayesKaon = new TH2D("fh2PIDProtonBayesKaon","PID: p: Bayes probability (K hyp.); #it{p}_{T} (GeV/#it{c}); Bayes prob.", 200,0,20, 50,0,1);
        fQAPID->Add(fh2PIDProtonBayesKaon);
        fh2PIDProtonBayesProton = new TH2D("fh2PIDProtonBayesProton","PID: p: Bayes probability (p hyp.); #it{p}_{T} (GeV/#it{c}); Bayes prob.", 200,0,20, 50,0,1);
        fQAPID->Add(fh2PIDProtonBayesProton);

      }

    } //endif {fProcessPID}

    if(fProcessPhi)
    {
      TString sPhiCounterLabel[] = {"Input","InvMass","Pt","Eta","Selected","Unlike-sign","BG"};
      const Short_t iNBinsPhiCounter = sizeof(sPhiCounterLabel)/sizeof(sPhiCounterLabel[0]);
      fhPhiCounter = new TH1D("fhPhiCounter","#phi: Counter",iNBinsPhiCounter,0,iNBinsPhiCounter);
      for(Short_t i(0); i < iNBinsPhiCounter; i++) fhPhiCounter->GetXaxis()->SetBinLabel(i+1, sPhiCounterLabel[i].Data() );
      fQAPhi->Add(fhPhiCounter);

      fhPhiMult = new TH1D("fhPhiMult","#phi: Multiplicity; Multiplicity", 150,0,150);
      fQAPhi->Add(fhPhiMult);
      fhPhiBGMult = new TH1D("fhPhiBGMult","#phi (BG): Multiplicity; Multiplicity", 150,0,150);
      fQAPhi->Add(fhPhiBGMult);
      fhPhiInvMass = new TH1D("fhPhiInvMass","#phi: InvMass; #it{m}_{inv} (GeV/#it{c}^{2})", 90,fCutPhiInvMassMin,fCutPhiInvMassMax);
      fQAPhi->Add(fhPhiInvMass);
      fhPhiBGInvMass = new TH1D("fhPhiBGInvMass","#phi (BG): InvMass; #it{m}_{inv} (GeV/#it{c}^{2})", 90,fCutPhiInvMassMin,fCutPhiInvMassMax);
      fQAPhi->Add(fhPhiBGInvMass);
      fhPhiCharge = new TH1D("fhPhiCharge","#phi: charge; charge", 5,-2.5,2.5);
      fQAPhi->Add(fhPhiCharge);
      fhPhiBGCharge = new TH1D("fhPhiBGCharge","#phi (BG): charge; charge", 5,-2.5,2.5);
      fQAPhi->Add(fhPhiBGCharge);
      fhPhiPt = new TH1D("fhPhiPt","#phi: #it{p}_{T}; #it{p}_{T} (GeV/#it{c})", 150,0.,30.);
      fQAPhi->Add(fhPhiPt);
      fhPhiEta = new TH1D("fhPhiEta","#phi: #eta; #eta", 151,-1.5,1.5);
      fQAPhi->Add(fhPhiEta);
      fhPhiPhi = new TH1D("fhPhiPhi","#phi: #varphi; #varphi", 100,-TMath::Pi(),TMath::Pi());
      fQAPhi->Add(fhPhiPhi);
    } //endif {fProcessPhi}

    // V0 candidates histograms
    if(fProcessV0s)
    {
      TString sV0sCounterLabel[] = {"Input","Daughters OK","Mother acceptance","Daughter acceptance","Charge","Reconstruction method","TPC refit","Kinks","Daughters track quality","DCA to PV","Daughters DCA","Decay radius","Common passed","K^{0}_{S}","#Lambda/#bar{#Lambda}","K^{0}_{S} && #Lambda/#bar{#Lambda}"};
      const Short_t iNBinsV0sCounter = sizeof(sV0sCounterLabel)/sizeof(sV0sCounterLabel[0]);
      fhV0sCounter = new TH1D("fhV0sCounter","V^{0}: Counter",iNBinsV0sCounter,0,iNBinsV0sCounter);
      for(Short_t i(0); i < iNBinsV0sCounter; i++) fhV0sCounter->GetXaxis()->SetBinLabel(i+1, sV0sCounterLabel[i].Data() );
      fQAV0s->Add(fhV0sCounter);

      TString sV0sK0sCounterLabel[] = {"Input","#it{y}","InvMass","CPA","Armenteros-Podolanski","c#tau","Daughters PID","Competing InvMass","Selected"};
      const Short_t iNBinsV0sK0sCounter = sizeof(sV0sK0sCounterLabel)/sizeof(sV0sK0sCounterLabel[0]);
      fhV0sCounterK0s = new TH1D("fhV0sCounterK0s","V^{0}: K^{0}_{S} Counter",iNBinsV0sK0sCounter,0,iNBinsV0sK0sCounter);
      for(Short_t i(0); i < iNBinsV0sK0sCounter; i++) fhV0sCounterK0s->GetXaxis()->SetBinLabel(i+1, sV0sK0sCounterLabel[i].Data() );
      fQAV0s->Add(fhV0sCounterK0s);

      TString sV0sLambdaCounterLabel[] = {"Input","#it{y}","InvMass","CPA","Armenteros-Podolanski","c#tau","Daughter PID","Competing InvMass","Selected","only #Lambda","only #bar{#Lambda}","#Lambda && #bar{#Lambda}"};
      const Short_t iNBinsV0sLambdaCounter = sizeof(sV0sLambdaCounterLabel)/sizeof(sV0sLambdaCounterLabel[0]);
      fhV0sCounterLambda = new TH1D("fhV0sCounterLambda","V^{0}: #Lambda/#bar{#Lambda} Counter",iNBinsV0sLambdaCounter,0,iNBinsV0sLambdaCounter);
      for(Short_t i(0); i < iNBinsV0sLambdaCounter; i++) fhV0sCounterLambda->GetXaxis()->SetBinLabel(i+1, sV0sLambdaCounterLabel[i].Data() );
      fQAV0s->Add(fhV0sCounterLambda);

      fhV0sInvMassK0s = new TH2D("fhV0sInvMassK0s","V^{0}: K^{0}_{S}: InvMass (selected); K^{0}_{S} #it{m}_{inv} (GeV/#it{c}^{2}); #Lambda/#bar{#Lambda} #it{m}_{inv} (GeV/#it{c}^{2})", 110,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax, 50,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax);
      fQAV0s->Add(fhV0sInvMassK0s);
      fhV0sInvMassLambda = new TH2D("fhV0sInvMassLambda","V^{0}: #Lambda/#bar{#Lambda}: InvMass (selected); K^{0}_{S} #it{m}_{inv} (GeV/#it{c}^{2}); #Lambda/#bar{#Lambda} #it{m}_{inv} (GeV/#it{c}^{2})", 110,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax, 50,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax);
      fQAV0s->Add(fhV0sInvMassLambda);
      fhV0sCompetingInvMassK0s = new TH2D("fhV0sCompetingInvMassK0s","V^{0}: K^{0}_{S}: Competing InvMass rejection; K^{0}_{S} #it{m}_{inv} (GeV/#it{c}^{2}); #Lambda/#bar{#Lambda} #it{m}_{inv} (GeV/#it{c}^{2})", 110,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax, 50,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax);
      fQAV0s->Add(fhV0sCompetingInvMassK0s);
      fhV0sCompetingInvMassLambda = new TH2D("fhV0sCompetingInvMassLambda","V^{0}: #Lambda/#bar{#Lambda}: Competing InvMass rejection; K^{0}_{S} #it{m}_{inv} (GeV/#it{c}^{2}); #Lambda/#bar{#Lambda} #it{m}_{inv} (GeV/#it{c}^{2})", 110,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax, 50,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax);
      fQAV0s->Add(fhV0sCompetingInvMassLambda);
    } // endif {fProcessV0s}

    const Short_t iNBinsPIDstatus = 4;
    TString sPIDstatus[iNBinsPIDstatus] = {"kDetNoSignal","kDetPidOk","kDetMismatch","kDetNoParams"};
    const Short_t iNFilterMapBinBins = 32;

    // QA histograms
    if(fFillQA)
    {
      TString sQAindex[fiNumIndexQA] = {"Before", "After"};
      for(Short_t iQA(0); iQA < fiNumIndexQA; iQA++)
      {
        // EVENTs QA histograms
        fhQAEventsPVz[iQA] = new TH1D(Form("fhQAEventsPVz_%s",sQAindex[iQA].Data()), "QA Events: PV-#it{z}", 101,-50,50);
        fQAEvents->Add(fhQAEventsPVz[iQA]);
        fhQAEventsNumContrPV[iQA] = new TH1D(Form("fhQAEventsNumContrPV_%s",sQAindex[iQA].Data()), "QA Events: Number of contributors to AOD PV", 20,0,20);
        fQAEvents->Add(fhQAEventsNumContrPV[iQA]);
        fhQAEventsNumSPDContrPV[iQA] = new TH1D(Form("fhQAEventsNumSPDContrPV_%s",sQAindex[iQA].Data()), "QA Events: SPD contributors to PV", 20,0,20);
        fQAEvents->Add(fhQAEventsNumSPDContrPV[iQA]);
        fhQAEventsDistPVSPD[iQA] = new TH1D(Form("fhQAEventsDistPVSPD_%s",sQAindex[iQA].Data()), "QA Events: PV SPD vertex", 50,0,5);
        fQAEvents->Add(fhQAEventsDistPVSPD[iQA]);
        fhQAEventsSPDresol[iQA] = new TH1D(Form("fhQAEventsSPDresol_%s",sQAindex[iQA].Data()), "QA Events: SPD resolution", 150,0,15);
        fQAEvents->Add(fhQAEventsSPDresol[iQA]);

        // Charged tracks QA
        fhQAChargedMult[iQA] = new TH1D(Form("fhQAChargedMult_%s",sQAindex[iQA].Data()),"QA Charged: Number of Charged in selected events; #it{N}^{Charged}", 1500,0,1500);
        fQACharged->Add(fhQAChargedMult[iQA]);
        fhQAChargedCharge[iQA] = new TH1D(Form("fhQAChargedCharge_%s",sQAindex[iQA].Data()),"QA Charged: Track charge; charge;", 3,-1.5,1.5);
        fQACharged->Add(fhQAChargedCharge[iQA]);
        fhQAChargedPt[iQA] = new TH1D(Form("fhQAChargedPt_%s",sQAindex[iQA].Data()),"QA Charged: Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c})", 300,0.,30.);
        fQACharged->Add(fhQAChargedPt[iQA]);
        fhQAChargedEta[iQA] = new TH1D(Form("fhQAChargedEta_%s",sQAindex[iQA].Data()),"QA Charged: Track #it{#eta}; #it{#eta}", 151,-1.5,1.5);
        fQACharged->Add(fhQAChargedEta[iQA]);
        fhQAChargedPhi[iQA] = new TH1D(Form("fhQAChargedPhi_%s",sQAindex[iQA].Data()),"QA Charged: Track #it{#varphi}; #it{#varphi}", 100,0.,TMath::TwoPi());
        fQACharged->Add(fhQAChargedPhi[iQA]);
        fhQAChargedFilterBit[iQA] = new TH1D(Form("fhQAChargedFilterBit_%s",sQAindex[iQA].Data()), "QA Charged: Filter bit",iNFilterMapBinBins,0,iNFilterMapBinBins);
        for(Int_t j = 0x0; j < iNFilterMapBinBins; j++) fhQAChargedFilterBit[iQA]->GetXaxis()->SetBinLabel(j+1, Form("%g",TMath::Power(2,j)));
        fQACharged->Add(fhQAChargedFilterBit[iQA]);
        fhQAChargedNumTPCcls[iQA] = new TH1D(Form("fhQAChargedNumTPCcls_%s",sQAindex[iQA].Data()),"QA Charged: Track number of TPC clusters; #it{N}^{TPC clusters}", 160,0,160);
        fQACharged->Add(fhQAChargedNumTPCcls[iQA]);
        fhQAChargedDCAxy[iQA] = new TH1D(Form("fhQAChargedDCAxy_%s",sQAindex[iQA].Data()),"QA Charged: Track DCA-xy; DCA_{#it{xy}} (cm)", 100,0.,10);
        fQACharged->Add(fhQAChargedDCAxy[iQA]);
        fhQAChargedDCAz[iQA] = new TH1D(Form("fhQAChargedDCAz_%s",sQAindex[iQA].Data()),"QA Charged: Track DCA-z; DCA_{#it{z}} (cm)", 200,-10.,10.);
        fQACharged->Add(fhQAChargedDCAz[iQA]);

        // PID tracks QA
        if(fProcessPID || fProcessPhi)
        {
          fhQAPIDTPCstatus[iQA] = new TH1D(Form("fhQAPIDTPCstatus_%s",sQAindex[iQA].Data()),"QA PID: PID status: TPC;", iNBinsPIDstatus,0,iNBinsPIDstatus);
          fQAPID->Add(fhQAPIDTPCstatus[iQA]);
          fhQAPIDTPCdEdx[iQA] = new TH2D(Form("fhQAPIDTPCdEdx_%s",sQAindex[iQA].Data()),"QA PID: TPC PID information; #it{p} (GeV/#it{c}); TPC dEdx (au)", 100,0,10, 131,-10,1000);
          fQAPID->Add(fhQAPIDTPCdEdx[iQA]);
          fhQAPIDTOFstatus[iQA] = new TH1D(Form("fhQAPIDTOFstatus_%s",sQAindex[iQA].Data()),"QA PID: PID status: TOF;", iNBinsPIDstatus,0,iNBinsPIDstatus);
          fQAPID->Add(fhQAPIDTOFstatus[iQA]);
          fhQAPIDTOFbeta[iQA] = new TH2D(Form("fhQAPIDTOFbeta_%s",sQAindex[iQA].Data()),"QA PID: TOF #beta information; #it{p} (GeV/#it{c}); TOF #beta", 100,0,10, 101,-0.1,1.5);
          fQAPID->Add(fhQAPIDTOFbeta[iQA]);

          fh3QAPIDnSigmaTPCTOFPtPion[iQA] = new TH3D(Form("fh3QAPIDnSigmaTPCTOFPtPion_%s",sQAindex[iQA].Data()), "QA PID: nSigma Pion vs. p_{T}; n#sigma TPC; n#sigma TOF; p_{T} (GeV/c)", 21,-11,10, 21,-11,10, iPOIsPtNumBins, fFlowPOIsPtMin, fFlowPOIsPtMax);
          fQAPID->Add(fh3QAPIDnSigmaTPCTOFPtPion[iQA]);
          fh3QAPIDnSigmaTPCTOFPtKaon[iQA] = new TH3D(Form("fh3QAPIDnSigmaTPCTOFPtKaon_%s",sQAindex[iQA].Data()), "QA PID: nSigma Kaon vs. p_{T}; n#sigma TPC; n#sigma TOF; p_{T} (GeV/c)", 21,-11,10, 21,-11,10, iPOIsPtNumBins, fFlowPOIsPtMin, fFlowPOIsPtMax);
          fQAPID->Add(fh3QAPIDnSigmaTPCTOFPtKaon[iQA]);
          fh3QAPIDnSigmaTPCTOFPtProton[iQA] = new TH3D(Form("fh3QAPIDnSigmaTPCTOFPtProton_%s",sQAindex[iQA].Data()), "QA PID: nSigma Proton vs. p_{T}; n#sigma TPC; n#sigma TOF; p_{T} (GeV/c)", 21,-11,10, 21,-11,10, iPOIsPtNumBins, fFlowPOIsPtMin, fFlowPOIsPtMax);
          fQAPID->Add(fh3QAPIDnSigmaTPCTOFPtProton[iQA]);

          fh3QAPIDnSigmaBayesElectron[iQA] = new TH3D(Form("fh3QAPIDnSigmaBayesElectron_%s",sQAindex[iQA].Data()),"QA PID: e; n#sigma^{TPC}; n#sigma^{TOF}; Bayes prob.", 21,-11,10, 21,-11,10, 22,-0.1,1);
          fQAPID->Add(fh3QAPIDnSigmaBayesElectron[iQA]);
          fh3QAPIDnSigmaBayesMuon[iQA] = new TH3D(Form("fh3QAPIDnSigmaBayesMuon_%s",sQAindex[iQA].Data()),"QA PID: #mu; n#sigma^{TPC}; n#sigma^{TOF}; Bayes prob.", 21,-11,10, 21,-11,10, 22,-0.1,1);
          fQAPID->Add(fh3QAPIDnSigmaBayesMuon[iQA]);
          fh3QAPIDnSigmaBayesPion[iQA] = new TH3D(Form("fh3QAPIDnSigmaBayesPion_%s",sQAindex[iQA].Data()),"QA PID: #pi; n#sigma^{TPC}; n#sigma^{TOF}; Bayes prob.", 21,-11,10, 21,-11,10, 22,-0.1,1);
          fQAPID->Add(fh3QAPIDnSigmaBayesPion[iQA]);
          fh3QAPIDnSigmaBayesKaon[iQA] = new TH3D(Form("fh3QAPIDnSigmaBayesKaon_%s",sQAindex[iQA].Data()),"QA PID: K; n#sigma^{TPC}; n#sigma^{TOF}; Bayes prob.", 21,-11,10, 21,-11,10, 22,-0.1,1);
          fQAPID->Add(fh3QAPIDnSigmaBayesKaon[iQA]);
          fh3QAPIDnSigmaBayesProton[iQA] = new TH3D(Form("fh3QAPIDnSigmaBayesProton_%s",sQAindex[iQA].Data()),"QA PID: p; n#sigma^{TPC}; n#sigma^{TOF}; Bayes prob.", 21,-11,10, 21,-11,10, 22,-0.1,1);
          fQAPID->Add(fh3QAPIDnSigmaBayesProton[iQA]);

          for(Int_t j = 0x0; j < iNBinsPIDstatus; j++)
          {
            fhQAPIDTOFstatus[iQA]->GetXaxis()->SetBinLabel(j+1, sPIDstatus[j].Data() );
            fhQAPIDTPCstatus[iQA]->GetXaxis()->SetBinLabel(j+1, sPIDstatus[j].Data() );
          }
        } // endif {fProcessPID}

        // V0s QA
        if(fProcessV0s)
        {
          fhQAV0sMultK0s[iQA] = new TH1D(Form("fhQAV0sMultK0s_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Number of K^{0}_{S} candidates", 1000,0,1000);
          fQAV0s->Add(fhQAV0sMultK0s[iQA]);
          fhQAV0sMultLambda[iQA] = new TH1D(Form("fhQAV0sMultLambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Number of #Lambda candidates", 1000,0,1000);
          fQAV0s->Add(fhQAV0sMultLambda[iQA]);
          fhQAV0sMultALambda[iQA] = new TH1D(Form("fhQAV0sMultALambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Number of #bar{#Lambda} candidates", 1000,0,1000);
          fQAV0s->Add(fhQAV0sMultALambda[iQA]);
          fhQAV0sRecoMethod[iQA] = new TH1D(Form("fhQAV0sRecoMethod_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Reconstruction method", 2,-0.5,1.5);
          fhQAV0sRecoMethod[iQA]->GetXaxis()->SetBinLabel(1, "offline");
          fhQAV0sRecoMethod[iQA]->GetXaxis()->SetBinLabel(2, "online (on-the-fly)");
          fQAV0s->Add(fhQAV0sRecoMethod[iQA]);
          fhQAV0sDCAtoPV[iQA] = new TH1D(Form("fhQAV0sDCAtoPV_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter DCA to PV; daughter DCA^{PV} (cm)", 250,0.0,5.0);
          fQAV0s->Add(fhQAV0sDCAtoPV[iQA]);
          fhQAV0sDCADaughters[iQA] = new TH1D(Form("fhQAV0sDCADaughters_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: DCA among daughters; DCA^{daughters} (cm)", 200,0.,20.);
          fQAV0s->Add(fhQAV0sDCADaughters[iQA]);
          fhQAV0sDecayRadius[iQA] = new TH1D(Form("fhQAV0sDecayRadius_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Decay radius; #it{r_{xy}}^{decay} (cm)", 400,0.0,200.0);
          fQAV0s->Add(fhQAV0sDecayRadius[iQA]);
          fhQAV0sCPAK0s[iQA] = new TH1D(Form("fhQAV0sCPAK0s_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: K^{0}_{S}: CPA; CPA^{K0s}", 100,0.9,1.);
          fQAV0s->Add(fhQAV0sCPAK0s[iQA]);
          fhQAV0sCPALambda[iQA] = new TH1D(Form("fhQAV0sCPALambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: #Lambda/#bar{#Lambda}: CPA; CPA^{#Lambda}", 100, 0.9,1.);
          fQAV0s->Add(fhQAV0sCPALambda[iQA]);
          fhQAV0sNumTauK0s[iQA] = new TH1D(Form("fhQAV0sNumTauK0s_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}:  K^{0}_{S}: Number of #it{c#tau}; #it{c#tau}^{K0s} (cm)", 100, 0.0,60.0);
          fQAV0s->Add(fhQAV0sNumTauK0s[iQA]);
          fhQAV0sNumTauLambda[iQA] = new TH1D(Form("fhQAV0sNumTauLambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Number of #it{c#tau}; #it{c#tau}^{#Lambda} (cm)", 100, 0.0,60.0);
          fQAV0s->Add(fhQAV0sNumTauLambda[iQA]);
          fhQAV0sArmenterosK0s[iQA] = new TH2D(Form("fhQAV0sArmenterosK0s_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}:  K^{0}_{S}: Armenteros-Podolaski plot; #alpha; #it{p}_{T}^{Arm} (GeV/#it{c});", 100,-1.,1., 100,0.,0.3);
          fQAV0s->Add(fhQAV0sArmenterosK0s[iQA]);
          fhQAV0sArmenterosLambda[iQA] = new TH2D(Form("fhQAV0sArmenterosLambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: #Lambda: Armenteros-Podolaski plot; #alpha; #it{p}_{T}^{Arm} (GeV/#it{c});", 100,-1.,1., 100,0.,0.3);
          fQAV0s->Add(fhQAV0sArmenterosLambda[iQA]);
          fhQAV0sArmenterosALambda[iQA] = new TH2D(Form("fhQAV0sArmenterosALambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: #bar{#Lambda}: Armenteros-Podolaski plot; #alpha; #it{p}_{T}^{Arm} (GeV/#it{c});", 100,-1.,1., 100,0.,0.3);
          fQAV0s->Add(fhQAV0sArmenterosALambda[iQA]);
          fhQAV0sInvMassK0s[iQA] = new TH1D(Form("fhQAV0sInvMassK0s_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: K^{0}_{S}: InvMass; #it{m}_{inv} (GeV/#it{c}^{2});", 200,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax);
          fQAV0s->Add(fhQAV0sInvMassK0s[iQA]);
          fhQAV0sInvMassLambda[iQA] = new TH1D(Form("fhQAV0sInvMassLambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: #Lambda/#bar{#Lambda}: InvMass; #it{m}_{inv} (GeV/#it{c}^{2});", 80,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax);
          fQAV0s->Add(fhQAV0sInvMassLambda[iQA]);
          fhQAV0sMotherPt[iQA] = new TH1D(Form("fhQAV0sMotherPt_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Mother #it{p}_{T}; #it{p}_{T}^{V0} (GeV/#it{c})", 200,0.,20.);
          fQAV0s->Add(fhQAV0sMotherPt[iQA]);
          fhQAV0sMotherPhi[iQA] = new TH1D(Form("fhQAV0sMotherPhi_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Mother #it{#varphi}; #it{#varphi}^{V0} (GeV/#it{c})", 100,0.,TMath::TwoPi());
          fQAV0s->Add(fhQAV0sMotherPhi[iQA]);
          fhQAV0sMotherEta[iQA] = new TH1D(Form("fhQAV0sMotherEta_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Mother #it{#eta}; #it{#eta}^{V0}", 151,-1.5,1.5);
          fQAV0s->Add(fhQAV0sMotherEta[iQA]);
          fhQAV0sMotherCharge[iQA] = new TH1D(Form("fhQAV0sMotherCharge_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Mother charge; V^{0} charge", 3,-1.5,1.5);
          fQAV0s->Add(fhQAV0sMotherCharge[iQA]);
          fhQAV0sMotherRapK0s[iQA] = new TH1D(Form("fhQAV0sMotherRapK0s_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Mother #it{y} (K^{0}_{S} hypo); #it{y}^{V0,K0s}", 151,-1.5,1.5);
          fQAV0s->Add(fhQAV0sMotherRapK0s[iQA]);
          fhQAV0sMotherRapLambda[iQA] = new TH1D(Form("fhQAV0sMotherRapLambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Mother #it{y} (Lambda/#bar{#Lambda} hypo); #it{y}^{V0,#Lambda}", 301,-3.,3.);
          fQAV0s->Add(fhQAV0sMotherRapLambda[iQA]);
          fhQAV0sDaughterTPCRefit[iQA] = new TH1D(Form("fhQAV0sDaughterTPCRefit_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter TPC refit", 2,-0.5,1.5);
          fhQAV0sDaughterTPCRefit[iQA]->GetXaxis()->SetBinLabel(1, "NOT AliAODTrack::kTPCrefit");
          fhQAV0sDaughterTPCRefit[iQA]->GetXaxis()->SetBinLabel(2, "AliAODTrack::kTPCrefit");
          fQAV0s->Add(fhQAV0sDaughterTPCRefit[iQA]);
          fhQAV0sDaughterKinks[iQA] = new TH1D(Form("fhQAV0sDaughterKinks_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter Kinks", 2,-0.5,1.5);
          fhQAV0sDaughterKinks[iQA]->GetXaxis()->SetBinLabel(1, "NOT AliAODVertex::kKink");
          fhQAV0sDaughterKinks[iQA]->GetXaxis()->SetBinLabel(2, "AliAODVertex:kKink");
          fQAV0s->Add(fhQAV0sDaughterKinks[iQA]);
          fhQAV0sDaughterNumTPCCls[iQA] = new TH1D(Form("fhQAV0sDaughterNumTPCCls_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter # of TPC clusters; # cls; Counts;", 180,-10.,170);
          fQAV0s->Add(fhQAV0sDaughterNumTPCCls[iQA]);
          fhQAV0sDaughterNumTPCClsPID[iQA] = new TH1D(Form("fhQAV0sDaughterNumTPCClsPID_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter # of TPC clusters for PID; # cls PID; Counts;", 180,-10.,170);
          fQAV0s->Add(fhQAV0sDaughterNumTPCClsPID[iQA]);
          fhQAV0sDaughterNumTPCFind[iQA] = new TH1D(Form("fhQAV0sDaughterNumTPCFind_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter # of findable TPC clusters; # findable; Counts;", 180,-10.,170);
          fQAV0s->Add(fhQAV0sDaughterNumTPCFind[iQA]);
          fhQAV0sDaughterNumTPCCrossRows[iQA] = new TH1D(Form("fhQAV0sDaughterNumTPCCrossRows_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter # of crossed TPC rows; # crossed; Counts;", 180,-10.,170);
          fQAV0s->Add(fhQAV0sDaughterNumTPCCrossRows[iQA]);
          fhQAV0sDaughterTPCCrossFindRatio[iQA] = new TH1D(Form("fhQAV0sDaughterTPCCrossFindRatio_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter crossed / findable TPC clusters; #crossed/find; Counts;", 50,0,5);
          fQAV0s->Add(fhQAV0sDaughterTPCCrossFindRatio[iQA]);
          fhQAV0sDaughterPt[iQA] = new TH1D(Form("fhQAV0sDaughterPt_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter #it{p}_{T}; #it{p}_{T}^{daughter} (GeV/#it{c})", 200,0.,20.);
          fQAV0s->Add(fhQAV0sDaughterPt[iQA]);
          fhQAV0sDaughterPhi[iQA] = new TH1D(Form("fhQAV0sDaughterPhi_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter #it{#varphi}; #it{#varphi}^{daughter} (GeV/#it{c})", 100,0.,TMath::TwoPi());
          fQAV0s->Add(fhQAV0sDaughterPhi[iQA]);
          fhQAV0sDaughterEta[iQA] = new TH1D(Form("fhQAV0sDaughterEta_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter #it{#eta}; #it{#eta}^{daugter}", 151,-1.5,1.5);
          fQAV0s->Add(fhQAV0sDaughterEta[iQA]);
          fhQAV0sDaughterCharge[iQA] = new TH1D(Form("fhQAV0sDaughterCharge_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter charge; daughter charge", 3,-1.5,1.5);
          fQAV0s->Add(fhQAV0sDaughterCharge[iQA]);
          fhQAV0sDaughterTPCstatus[iQA] = new TH1D(Form("fhQAV0sDaughterTPCstatus_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: PID status: TPC;", iNBinsPIDstatus,0,iNBinsPIDstatus);
          fQAV0s->Add(fhQAV0sDaughterTPCstatus[iQA]);
          fhQAV0sDaughterTOFstatus[iQA] = new TH1D(Form("fhQAV0sDaughterTOFstatus_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: PID status: TOF;", iNBinsPIDstatus,0,iNBinsPIDstatus);
          fQAV0s->Add(fhQAV0sDaughterTOFstatus[iQA]);
          fhQAV0sDaughterTPCdEdxK0s[iQA] = new TH2D(Form("fhQAV0sDaughterTPCdEdxK0s_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: K^{0}_{S}: TPC dEdx daughters; #it{p}^{daughter} (GeV/#it{c}); TPC dEdx (au);", 100,0.,20, 101,-10,1000);
          fQAV0s->Add(fhQAV0sDaughterTPCdEdxK0s[iQA]);
          fhQAV0sDaughterNumSigmaPionK0s[iQA] = new TH2D(Form("fhQAV0sDaughterNumSigmaPionK0s_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: K^{0}_{S}: Daughter PID (#pi); #it{p}_{T}^{daughter} (GeV/#it{c}); #pi PID (#sigma^{TPC});", 200,0.,20, 100,-10.,10.);
          fQAV0s->Add(fhQAV0sDaughterNumSigmaPionK0s[iQA]);
          fhQAV0sDaughterTPCdEdxLambda[iQA] = new TH2D(Form("fhQAV0sDaughterTPCdEdxLambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: #Lambda/#bar{#Lambda}: TPC dEdx daughters; #it{p}^{daughter} (GeV/#it{c}); TPC dEdx (au);", 100,0.,20, 101,-10,1000);
          fQAV0s->Add(fhQAV0sDaughterTPCdEdxLambda[iQA]);
          fhQAV0sDaughterNumSigmaPionLambda[iQA] = new TH2D(Form("fhQAV0sDaughterNumSigmaPionLambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: #Lambda: Daughter PID (#pi); #it{p}_{T}^{pion} (GeV/#it{c}); pion PID (#sigma^{TPC});", 200,0.,20, 100,-10.,10.);
          fQAV0s->Add(fhQAV0sDaughterNumSigmaPionLambda[iQA]);
          fhQAV0sDaughterNumSigmaProtonLambda[iQA] = new TH2D(Form("fhQAV0sDaughterNumSigmaProtonLambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: #Lambda: Daughter PID (p); #it{p}_{T}^{proton} (GeV/#it{c}); proton PID (#sigma^{TPC});", 200,0.,20, 100,-10.,10.);
          fQAV0s->Add(fhQAV0sDaughterNumSigmaProtonLambda[iQA]);
          fhQAV0sDaughterNumSigmaPionALambda[iQA] = new TH2D(Form("fhQAV0sDaughterNumSigmaPionALambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: #bar{#Lambda}: Daughter PID (#pi); #it{p}_{T}^{pion} (GeV/#it{c}); pion PID (#sigma^{TPC});", 200,0.,20, 100,-10.,10.);
          fQAV0s->Add(fhQAV0sDaughterNumSigmaPionALambda[iQA]);
          fhQAV0sDaughterNumSigmaProtonALambda[iQA] = new TH2D(Form("fhQAV0sDaughterNumSigmaProtonALambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: #bar{#Lambda}: Daughter PID (p); #it{p}_{T}^{proton} (GeV/#it{c}); proton PID (#sigma^{TPC});", 200,0.,20, 100,-10.,10.);
          fQAV0s->Add(fhQAV0sDaughterNumSigmaProtonALambda[iQA]);



          for(Int_t j = 0x0; j < iNBinsPIDstatus; j++)
          {
            fhQAV0sDaughterTOFstatus[iQA]->GetXaxis()->SetBinLabel(j+1, sPIDstatus[j].Data() );
            fhQAV0sDaughterTPCstatus[iQA]->GetXaxis()->SetBinLabel(j+1, sPIDstatus[j].Data() );
          }
        } // endif {fProcessV0s}
      }
    }

  // posting data (mandatory)
  PostData(1, fFlowRefs);
  PostData(2, fFlowCharged);
  PostData(3, fFlowPID);
  PostData(4, fFlowPhi);
  PostData(5, fFlowK0s);
  PostData(6, fFlowLambda);
  PostData(7, fQAEvents);
  PostData(8, fQACharged);
  PostData(9, fQAPID);
  PostData(10, fQAPhi);
  PostData(11, fQAV0s);
  PostData(12, fFlowWeights);

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::ListParameters()
{
  // lists all task parameters
  // *************************************************************
  AliInfo("Listing all AliAnalysisTaskUniFlow parameters");
  printf("   -------- Analysis task ---------------------------------------\n");
  printf("      fRunMode: (RunMode) %d\n",    fRunMode);
  printf("      fAnalType: (AnalType) %d\n",    fAnalType);
  printf("      fSampling: (Bool_t) %s\n",    fSampling ? "kTRUE" : "kFALSE");
  printf("      fFillQA: (Bool_t) %s\n",    fFillQA ? "kTRUE" : "kFALSE");
  printf("      fProcessPID: (Bool_t) %s\n",    fProcessPID ? "kTRUE" : "kFALSE");
  printf("      fProcessPhi: (Bool_t) %s\n",    fProcessPhi ? "kTRUE" : "kFALSE");
  printf("      fProcessV0s: (Bool_t) %s\n",    fProcessV0s ? "kTRUE" : "kFALSE");
  printf("   -------- Flow related ----------------------------------------\n");
  printf("      fCutFlowDoFourCorrelations: (Bool_t) %s\n",    fCutFlowDoFourCorrelations ? "kTRUE" : "kFALSE");
  printf("      fCutFlowRFPsPtMin: (Double_t) %g (GeV/c)\n",    fCutFlowRFPsPtMin);
  printf("      fCutFlowRFPsPtMax: (Double_t) %g (GeV/c)\n",    fCutFlowRFPsPtMax);
  printf("      fFlowPOIsPtMin: (Double_t) %g (GeV/c)\n",    fFlowPOIsPtMin);
  printf("      fFlowPOIsPtMax: (Double_t) %g (GeV/c)\n",    fFlowPOIsPtMax);
  printf("      fFlowCentMin: (Int_t) %d (GeV/c)\n",    fFlowCentMin);
  printf("      fFlowCentMax: (Int_t) %d (GeV/c)\n",    fFlowCentMax);
  printf("      fFlowUseWeights: (Bool_t) %s\n",    fFlowUseWeights ? "kTRUE" : "kFALSE");
  printf("      fFlowRunByRunWeights: (Bool_t) %s\n",    fFlowRunByRunWeights ? "kTRUE" : "kFALSE");
  printf("      fFlowWeightsPath: (TString) '%s' \n",    fFlowWeightsPath.Data());
  printf("   -------- Events ----------------------------------------------\n");
  printf("      fColSystem: (ColSystem) %d\n",    fColSystem);
  printf("      fTrigger: (Short_t) %d\n",    fTrigger);
  printf("      fMultEstimator: (TString) '%s'\n",    fMultEstimator.Data());
  printf("      fPVtxCutZ: (Double_t) %g (cm)\n",    fPVtxCutZ);
  printf("      fUseAliEventCuts: (Bool_t) %s\n",    fUseAliEventCuts ? "kTRUE" : "kFALSE");
  printf("   -------- Charged tracks --------------------------------------\n");
  printf("      fCutChargedTrackFilterBit: (UInt) %d\n",    fCutChargedTrackFilterBit);
  printf("      fCutChargedNumTPCclsMin: (UShort_t) %d\n",    fCutChargedNumTPCclsMin);
  printf("      fCutChargedEtaMax: (Double_t) %g\n",    fCutChargedEtaMax);
  printf("      fCutChargedDCAzMax: (Double_t) %g (cm)\n",    fCutChargedDCAzMax);
  printf("      fCutChargedDCAxyMax: (Double_t) %g (cm)\n",    fCutChargedDCAxyMax);
  printf("   -------- PID (pi,K,p) tracks ---------------------------------\n");
  printf("      fCutPIDUseAntiProtonOnly: (Bool_t) %s\n",  fCutPIDUseAntiProtonOnly ? "kTRUE" : "kFALSE");
  printf("      fCutPIDnSigmaCombinedTOFrejection: (Bool_t) %s\n",  fCutPIDnSigmaCombinedTOFrejection ? "kTRUE" : "kFALSE");
  printf("      fCutPIDnSigmaPionMax: (Float_t) %g\n",    fCutPIDnSigmaPionMax);
  printf("      fCutPIDnSigmaKaonMax: (Float_t) %g\n",    fCutPIDnSigmaKaonMax);
  printf("      fCutPIDnSigmaProtonMax: (Float_t) %g\n",    fCutPIDnSigmaProtonMax);
  printf("      fCutPIDnSigmaTPCRejectElectron: (Float_t) %g\n",    fCutPIDnSigmaTPCRejectElectron);
  printf("      fCutUseBayesPID: (Bool_t) %s\n",    fCutUseBayesPID ? "kTRUE" : "kFALSE");
  printf("      fCutPIDBayesPionMin: (Double_t) %g\n",    fCutPIDBayesPionMin);
  printf("      fCutPIDBayesKaonMin: (Double_t) %g\n",    fCutPIDBayesKaonMin);
  printf("      fCutPIDBayesProtonMin: (Double_t) %g\n",    fCutPIDBayesProtonMin);
  printf("      SetPIDBayesRejectElectron: (Double_t) %g\n",    fCutPIDBayesRejectElectron);
  printf("      SetPIDBayesRejectMuon: (Double_t) %g\n",    fCutPIDBayesRejectMuon);
  printf("   -------- Phi candidates --------------------------------------\n");
  printf("      fCutPhiMotherEtaMax: (Double_t) %g\n",    fCutPhiMotherEtaMax);
  printf("      fCutPhiInvMassMin: (Double_t) %g\n",    fCutPhiInvMassMin);
  printf("      fCutPhiInvMassMax: (Double_t) %g\n",    fCutPhiInvMassMax);
  printf("   -------- V0s candidates --------------------------------------\n");
  printf("      fCutV0sOnFly: (Bool_t) %s\n",    fCutV0sOnFly ? "kTRUE" : "kFALSE");
  printf("      fCutV0srefitTPC: (Bool_t) %s\n",     fCutV0srefitTPC ? "kTRUE" : "kFALSE");
  printf("      fCutV0srejectKinks: (Bool_t) %s\n",     fCutV0srejectKinks ? "kTRUE" : "kFALSE");
  printf("      fCutV0sDaughterNumTPCClsMin: (UShort_t) %d\n",     fCutV0sDaughterNumTPCClsMin);
  printf("      fCutV0sDaughterNumTPCClsPIDMin: (UShort_t) %d\n",     fCutV0sDaughterNumTPCClsPIDMin);
  printf("      fCutV0sDaughterNumTPCCrossMin: (UShort_t) %d\n",     fCutV0sDaughterNumTPCCrossMin);
  printf("      fCutV0sDaughterNumTPCFindMin: (UShort_t) %d\n",     fCutV0sDaughterNumTPCFindMin);
  printf("      fCutV0sDaughterRatioCrossFindMin: (Double_t) %g\n",     fCutV0sDaughterRatioCrossFindMin);
  printf("      fCutV0sCrossMassRejection: (Bool_t) %s\n",     fCutV0sCrossMassRejection ? "kTRUE" : "kFALSE");
  printf("      fCutV0sCrossMassCutK0s: (Double_t) %g (GeV/c2)\n",     fCutV0sCrossMassCutK0s);
  printf("      fCutV0sCrossMassCutLambda: (Double_t) %g (GeV/c2)\n",     fCutV0sCrossMassCutLambda);
  printf("      fCutV0sDCAtoPVMin: (Double_t) %g (cm)\n",    fCutV0sDCAtoPVMin);
  printf("      fCutV0sDCAtoPVMax: (Double_t) %g (cm)\n",    fCutV0sDCAtoPVMax);
  printf("      fCutV0sDCAtoPVzMax: (Double_t) %g (cm)\n",    fCutV0sDCAtoPVzMax);
  printf("      fCutV0sDCADaughtersMin: (Double_t) %g (cm)\n",    fCutV0sDCADaughtersMin);
  printf("      fCutV0sDCADaughtersMax: (Double_t) %g (cm)\n",    fCutV0sDCADaughtersMax);
  printf("      fCutV0sDecayRadiusMin: (Double_t) %g (cm)\n",    fCutV0sDecayRadiusMin);
  printf("      fCutV0sDecayRadiusMax: (Double_t) %g (cm)\n",    fCutV0sDecayRadiusMax);
  printf("      fCutV0sDaughterFilterBit: (UInt) %d\n",    fCutV0sDaughterFilterBit);
  printf("      fCutV0sDaughterPtMin: (Double_t) %g (GeV/c)\n",    fCutV0sDaughterPtMin);
  printf("      fCutV0sDaughterPtMax: (Double_t) %g (GeV/c)\n",    fCutV0sDaughterPtMax);
  printf("      fCutV0sDaughterEtaMax: (Double_t) %g ()\n",    fCutV0sDaughterEtaMax);
  printf("      fCutV0sMotherEtaMax: (Double_t) %g ()\n",    fCutV0sMotherEtaMax);
  printf("      fCutV0sMotherRapMax: (Double_t) %g ()\n",    fCutV0sMotherRapMax);
  printf("      fCutV0sCPAK0sMin: (Double_t) %g ()\n",    fCutV0sCPAK0sMin);
  printf("      fCutV0sCPALambdaMin: (Double_t) %g ()\n",    fCutV0sCPALambdaMin);
  printf("      fCutV0sNumTauK0sMax: (Double_t) %g (c*tau)\n",    fCutV0sNumTauK0sMax);
  printf("      fCutV0sNumTauLambdaMax: (Double_t) %g (c*tau)\n",    fCutV0sNumTauLambdaMax);
  printf("      fCutV0sInvMassK0sMin: (Double_t) %g (GeV/c2)\n",    fCutV0sInvMassK0sMin);
  printf("      fCutV0sInvMassK0sMax: (Double_t) %g (GeV/c2)\n",    fCutV0sInvMassK0sMax);
  printf("      fCutV0sInvMassLambdaMin: (Double_t) %g (GeV/c2)\n",    fCutV0sInvMassLambdaMin);
  printf("      fCutV0sInvMassLambdaMax: (Double_t) %g (GeV/c2)\n",    fCutV0sInvMassLambdaMax);
  printf("      fCutV0sArmenterosAlphaK0sMin: (Double_t) %g (alpha)\n",    fCutV0sArmenterosAlphaK0sMin);
  printf("      fCutV0sArmenterosAlphaLambdaMax: (Double_t) %g (alpha)\n",    fCutV0sArmenterosAlphaLambdaMax);
  printf("      fCutV0sK0sPionNumTPCSigmaMax: (Float_t) %g (n*sigma)\n",    fCutV0sK0sPionNumTPCSigmaMax);
  printf("      fCutV0sLambdaPionNumTPCSigmaMax: (Float_t) %g (n*sigma)\n",    fCutV0sLambdaPionNumTPCSigmaMax);
  printf("      fCutV0sLambdaProtonNumTPCSigmaMax: (Float_t) %g (n*sigma)\n",    fCutV0sLambdaProtonNumTPCSigmaMax);
  printf("=====================================================================\n\n");

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::ClearVectors()
{
  // Properly clear all particle vectors (if exists)
  // NOTE: should be called at the end of each event & before vectors deleting
  // *************************************************************

  // pointers owned by Ali*Event containers
  if(fVectorRefs) { fVectorRefs->clear(); }
  if(fVectorCharged) { fVectorCharged->clear(); }
  if(fVectorPion) { fVectorPion->clear(); }
  if(fVectorKaon) { fVectorKaon->clear(); }
  if(fVectorProton) { fVectorProton->clear(); }

  // pointers owned by task
  if(fVectorK0s)
  {
    for(auto part = fVectorK0s->begin(); part != fVectorK0s->end(); ++part) { delete *part; }
    fVectorK0s->clear();
  }

  if(fVectorLambda)
  {
    for(auto part = fVectorLambda->begin(); part != fVectorLambda->end(); ++part) { delete *part; }
    fVectorLambda->clear();
  }

  if(fVectorPhi)
  {
    for(auto part = fVectorPhi->begin(); part != fVectorPhi->end(); ++part) { delete *part; }
    fVectorPhi->clear();
  }

  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::InitializeTask()
{
  // called once on beginning of task (within UserCreateOutputObjects method)
  // check if task parameters are specified and valid
  // returns kTRUE if succesfull
  // *************************************************************
  AliInfo("Checking task setting");

  if(fAnalType != kAOD)
  {
    AliFatal("Analysis type: not kAOD (not implemented for ESDs yet)! Terminating!");
    return kFALSE;
  }

  // checking PID response
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*)mgr->GetInputEventHandler();
  fPIDResponse = inputHandler->GetPIDResponse();
  if(!fPIDResponse)
  {
    AliFatal("AliPIDResponse object not found! Terminating!");
    return kFALSE;
  }

  fPIDCombined = new AliPIDCombined();
  if(!fPIDCombined)
  {
    AliFatal("AliPIDCombined object not found! Terminating!");
    return kFALSE;
  }
  fPIDCombined->SetDefaultTPCPriors();
  fPIDCombined->SetSelectedSpecies(5); // all particle species
  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF); // setting TPC + TOF mask

  // checking the fFlowNumHarmonicsMax, fFlowNumWeightPowersMax dimensions of p,Q,S vectors
  for(Int_t iHarm(0); iHarm < fNumHarmonics; ++iHarm)
  {
    if(fFlowNumWeightPowersMax < 3) { AliFatal("Low range of flow vector weight dimension! Not enought for <2>!"); return kFALSE; }
    if(fCutFlowDoFourCorrelations && fFlowNumWeightPowersMax < 5) { AliFatal("Low range of flow vector weight dimension! Not enought for <4>!"); return kFALSE; }
    if(fFlowNumHarmonicsMax < fHarmonics[iHarm]+1) { AliFatal("Low range of flow vector harmonics dimension!"); return kFALSE; }
  }

  if(fSampling && fNumSamples == 0)
  {
    AliFatal("Sampling used, but number of samples is 0! Terminating!");
    return kFALSE;
  }

  // checking cut setting
  AliInfo("Checking task parameters setting conflicts (ranges, etc)");
  if(fCutFlowRFPsPtMin > 0. && fCutFlowRFPsPtMax > 0. && fCutFlowRFPsPtMin > fCutFlowRFPsPtMax)
  {
    AliFatal("Cut: RFPs Pt range wrong! Terminating!");
    return kFALSE;
  }
  if(fCutV0sInvMassK0sMin > fCutV0sInvMassK0sMax || fCutV0sInvMassK0sMin < 0. || fCutV0sInvMassK0sMax < 0.)
  {
    AliFatal("Cut: InvMass (K0s) range wrong! Terminating! ");
    return kFALSE;
  }
  if(fCutV0sInvMassLambdaMin > fCutV0sInvMassLambdaMax || fCutV0sInvMassLambdaMin < 0. || fCutV0sInvMassLambdaMax < 0.)
  {
    AliFatal("Cut: InvMass (Lambda) range wrong! Terminating!");
    return kFALSE;
  }

  // upper-case for multiplicity estimator
  fMultEstimator.ToUpper();

  // setting fFlowCentMax according to estimator : 0-100 for %'s and 0-150 for CHARGED
  if(fMultEstimator.EqualTo("") || fMultEstimator.EqualTo("CHARGED")) { fFlowCentMin = 0; fFlowCentMax = 200; }
  else { fFlowCentMin = 0; fFlowCentMax = 100; }

  // increasing fFlowCentMax+1 (just to be sure)
  fFlowCentMax++;

  // checking for weights source file
  if(fFlowUseWeights && !fFlowWeightsPath.EqualTo(""))
  {
    fFlowWeightsFile = TFile::Open(Form("%s",fFlowWeightsPath.Data()));
    if(!fFlowWeightsFile)
    {
      AliFatal("Flow weights file not found! Terminating!");
      return kFALSE;
    }

    // if only one (run integrated) set of weights is used, load it here
    if(!fFlowRunByRunWeights)
    {
      TList* listFlowWeights = (TList*) fFlowWeightsFile->Get("weights");
      if(!listFlowWeights) { AliError("TList from flow weights not found."); return kFALSE; }

      if(fFlowUse3Dweights)
      {
        fh3WeightRefs = (TH3D*) listFlowWeights->FindObject("Refs3D"); if(!fh3WeightRefs) { AliFatal("Refs weights not found"); return kFALSE; }
        fh3WeightCharged = (TH3D*) listFlowWeights->FindObject("Charged3D"); if(!fh3WeightCharged) { AliFatal("Charged weights not found"); return kFALSE; }
        fh3WeightPion = (TH3D*) listFlowWeights->FindObject("Pion3D"); if(!fh3WeightPion) { AliFatal("Pion weights not found"); return kFALSE; }
        fh3WeightKaon = (TH3D*) listFlowWeights->FindObject("Kaon3D"); if(!fh3WeightKaon) { AliFatal("Kaon weights not found"); return kFALSE; }
        fh3WeightProton = (TH3D*) listFlowWeights->FindObject("Proton3D"); if(!fh3WeightProton) { AliFatal("Proton weights not found"); return kFALSE; }
        fh3WeightK0s = (TH3D*) listFlowWeights->FindObject("K0s3D"); if(!fh3WeightK0s) { AliFatal("K0s weights not found"); return kFALSE; }
        fh3WeightLambda = (TH3D*) listFlowWeights->FindObject("Lambda3D"); if(!fh3WeightLambda) { AliFatal("Phi weights not found"); return kFALSE; }
        fh3WeightPhi = (TH3D*) listFlowWeights->FindObject("Phi3D"); if(!fh3WeightPhi) { AliFatal("Phi weights not found"); return kFALSE; }
      }
      else
      {
        fh2WeightRefs = (TH2D*) listFlowWeights->FindObject("Refs"); if(!fh2WeightRefs) { AliFatal("Refs weights not found"); return kFALSE; }
        fh2WeightCharged = (TH2D*) listFlowWeights->FindObject("Charged"); if(!fh2WeightCharged) { AliFatal("Charged weights not found"); return kFALSE; }
        fh2WeightPion = (TH2D*) listFlowWeights->FindObject("Pion"); if(!fh2WeightPion) { AliFatal("Pion weights not found"); return kFALSE; }
        fh2WeightKaon = (TH2D*) listFlowWeights->FindObject("Kaon"); if(!fh2WeightKaon) { AliFatal("Kaon weights not found"); return kFALSE; }
        fh2WeightProton = (TH2D*) listFlowWeights->FindObject("Proton"); if(!fh2WeightProton) { AliFatal("Proton weights not found"); return kFALSE; }
        fh2WeightK0s = (TH2D*) listFlowWeights->FindObject("K0s"); if(!fh2WeightK0s) { AliFatal("K0s weights not found"); return kFALSE; }
        fh2WeightLambda = (TH2D*) listFlowWeights->FindObject("Lambda"); if(!fh2WeightLambda) { AliFatal("Phi weights not found"); return kFALSE; }
        fh2WeightPhi = (TH2D*) listFlowWeights->FindObject("Phi"); if(!fh2WeightPhi) { AliFatal("Phi weights not found"); return kFALSE; }
      }
    }
  }

  AliInfo("Preparing particle containers (std::vectors)");
  // creating particle vectors & reserving capacity in order to avoid memory re-allocation
  // when capacity is not enough later during filtering
  // NOTE: system and cuts dependent (should be modified accordingly)

  fVectorRefs = new std::vector<AliVTrack*>;
  fVectorCharged = new std::vector<AliVTrack*>;
  fVectorPion = new std::vector<AliVTrack*>;
  fVectorKaon = new std::vector<AliVTrack*>;
  fVectorProton = new std::vector<AliVTrack*>;
  fVectorK0s = new std::vector<AliVTrack*>;
  fVectorLambda = new std::vector<AliVTrack*>;
  fVectorPhi = new std::vector<AliVTrack*>;

  switch(fColSystem)
  {
    case kPP :
      fVectorRefs->reserve(150);
      fVectorCharged->reserve(150);
      if(fProcessPID) { fVectorPion->reserve(100); fVectorKaon->reserve(20); fVectorProton->reserve(20); }
      if(fProcessV0s) { fVectorK0s->reserve(30); fVectorLambda->reserve(30); }
      if(fProcessPhi) { fVectorPhi->reserve(20); }
      break;

    case kPPb :
      fVectorRefs->reserve(200);
      fVectorCharged->reserve(200);
      if(fProcessPID) { fVectorPion->reserve(100); fVectorKaon->reserve(40); fVectorProton->reserve(30); }
      if(fProcessV0s) { fVectorK0s->reserve(100); fVectorLambda->reserve(150); }
      if(fProcessPhi) { fVectorPhi->reserve(20); }
      break;

    default :
      fVectorRefs->reserve(300);
      fVectorCharged->reserve(300);
      if(fProcessPID) { fVectorPion->reserve(200); fVectorKaon->reserve(100); fVectorProton->reserve(100); }
      if(fProcessV0s) { fVectorK0s->reserve(100); fVectorLambda->reserve(100); }
      if(fProcessPhi) { fVectorPhi->reserve(200); }
  }

  AliInfo("Initialization succesfull!");
  return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::UserExec(Option_t *)
{
  // main method called for each event (event loop)
  // *************************************************************

  // check if initialization succesfull (done within UserCreateOutputObjects())
  if(!fInit) return;

  // local event counter check: if running in test mode, it runs until the 50 events are succesfully processed
  if(fRunMode == kTest && fEventCounter >= fNumEventsAnalyse) return;

  // event selection
  fEventAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!EventSelection()) return;

  // loading array with MC particles
  if(fMC)
  {
    fArrayMC = (TClonesArray*) fEventAOD->FindListObject("mcparticles");
    if(!fArrayMC) { AliError("TClonesArray with MC particle not found!"); return; }
  }

  // processing of selected event
  Bool_t bProcessed = ProcessEvent();

  // need to be done no matter if the event was completely processed or not
  // because processing can return at various steps where tracks can be created/pushed to vectors
  ClearVectors();

  if(!bProcessed) return;

  // posting data (mandatory)
  PostData(1, fFlowRefs);
  PostData(2, fFlowCharged);
  PostData(3, fFlowPID);
  PostData(4, fFlowPhi);
  PostData(5, fFlowK0s);
  PostData(6, fFlowLambda);
  PostData(7, fQAEvents);
  PostData(8, fQACharged);
  PostData(9, fQAPID);
  PostData(10, fQAPhi);
  PostData(11, fQAV0s);
  PostData(12, fFlowWeights);

  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::EventSelection()
{
  // main (envelope) method for event selection
  // Specific event selection methods are called from here.
  // Fill the event QA if event pass selection.
  // returns kTRUE if event pass all selection criteria kFALSE otherwise
  // *************************************************************
  if(!fEventAOD) return kFALSE;
  fhEventCounter->Fill("Input",1);

  // Fill event QA BEFORE cuts
  if(fFillQA) FillEventsQA(0);

  Bool_t eventSelected = kFALSE;

  // event selection for small systems pp, pPb in Run2 (2016)
  if(fColSystem == kPP || fColSystem == kPPb) eventSelected = IsEventSelected_small_2016();

  if(!eventSelected) return kFALSE;

  fhEventCounter->Fill("Selected",1);

  // Fill event QA AFTER cuts
  if(fFillQA) FillEventsQA(1);

  // extract PV-z for weights
  fPVz = fEventAOD->GetPrimaryVertex()->GetZ();

  return eventSelected;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::IsEventSelected_small_2016()
{
  // Event selection for small system collision recorder in Run 2 year 2016
  // pp (LHC16kl...), pPb (LHC16rqts)
  // return kTRUE if event passes all criteria, kFALSE otherwise
  // *************************************************************

  // Physics selection (trigger)
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) mgr->GetInputEventHandler();
  UInt_t fSelectMask = inputHandler->IsEventSelected();

  Bool_t isTriggerSelected = kFALSE;
  switch(fTrigger) // check for high multiplicity trigger
  {
    case 0:
      isTriggerSelected = fSelectMask& AliVEvent::kINT7;
      break;

    case 1:
      isTriggerSelected = fSelectMask& AliVEvent::kHighMultV0;
      break;

    case 2:
      isTriggerSelected = fSelectMask& AliVEvent::kHighMultSPD;
      break;

    default:
      isTriggerSelected = kFALSE;
  }

  if(!isTriggerSelected)
    return kFALSE;

  // events passing physics selection
  fhEventCounter->Fill("Physics selection OK",1);

  // events passing AliEventCuts selection
  if(fUseAliEventCuts && !fEventCuts.AcceptEvent(fEventAOD)) return kFALSE;
  fhEventCounter->Fill("EventCuts OK",1);

  // primary vertex selection
  const AliAODVertex* vtx = dynamic_cast<const AliAODVertex*>(fEventAOD->GetPrimaryVertex());
  if(!vtx || vtx->GetNContributors() < 1)
    return kFALSE;
  fhEventCounter->Fill("PV OK",1);

  // SPD vertex selection
  const AliAODVertex* vtxSPD = dynamic_cast<const AliAODVertex*>(fEventAOD->GetPrimaryVertexSPD());

  Double_t dMaxResol = 0.25; // suggested from DPG
  Double_t cov[6] = {0};
  vtxSPD->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  if ( vtxSPD->IsFromVertexerZ() && (zRes > dMaxResol)) return kFALSE;
  fhEventCounter->Fill("SPD Vtx OK",1);

  // PileUp rejection included in Physics selection
  // but with values for high mult pp (> 5 contrib) => for low ones: do manually (> 3 contrib)

  /*
  if(fTrigger == 0 && fAOD->IsPileupFromSPD(3,0.8) )
  {
    return kFALSE;
  }
  */

  //fhEventCounter->Fill("Pileup SPD OK",1);

  // pileup rejection from multivertexer
  AliAnalysisUtils utils;
  utils.SetMinPlpContribMV(5);
  utils.SetMaxPlpChi2MV(5);
  utils.SetMinWDistMV(15);
  utils.SetCheckPlpFromDifferentBCMV(kFALSE);
  Bool_t isPileupFromMV = utils.IsPileUpMV(fEventAOD);

  if(isPileupFromMV) return kFALSE;
  fhEventCounter->Fill("Pileup MV OK",1);

  // if(fRejectOutOfBunchPU) // out-of-bunch rejection (provided by Christian)
  // {
  //   //out-of-bunch 11 BC
  //   if (utils.IsOutOfBunchPileUp(fEventAOD))
  //   {
  //     return kFALSE;
  //   }
  //   fhEventCounter->Fill("OOBPU OK",1);
  //
  //   if (utils.IsSPDClusterVsTrackletBG(fEventAOD))
  //   {
  //     return kFALSE;
  //   }
  //
  //   fhEventCounter->Fill("SPDClTrBG OK",1);
  //
  //   // SPD pileup
  //   if (utils.IsPileUpSPD(fEventAOD))
  //   {
  //     return kFALSE;
  //   }
  //
  //   fhEventCounter->Fill("SPDPU OK",1);
  // }

  //fhEventCounter->Fill("Utils OK",1);

  // cutting on PV z-distance
  const Double_t aodVtxZ = vtx->GetZ();
  if( TMath::Abs(aodVtxZ) > fPVtxCutZ )
  {
    return kFALSE;
  }
  fhEventCounter->Fill("PV #it{z} OK",1);

  return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FillEventsQA(const Short_t iQAindex)
{
  // Filling various QA plots related with event selection
  // *************************************************************

  const AliAODVertex* aodVtx = fEventAOD->GetPrimaryVertex();
  const Double_t dVtxZ = aodVtx->GetZ();
  const Int_t iNumContr = aodVtx->GetNContributors();
  const AliAODVertex* spdVtx = fEventAOD->GetPrimaryVertexSPD();
  const Int_t iNumContrSPD = spdVtx->GetNContributors();
  const Double_t spdVtxZ = spdVtx->GetZ();

  fhQAEventsPVz[iQAindex]->Fill(dVtxZ);
  fhQAEventsNumContrPV[iQAindex]->Fill(iNumContr);
  fhQAEventsNumSPDContrPV[iQAindex]->Fill(iNumContrSPD);
  fhQAEventsDistPVSPD[iQAindex]->Fill(TMath::Abs(dVtxZ - spdVtxZ));

  // // event / physics selection criteria
  // AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  // AliInputEventHandler* inputHandler = (AliInputEventHandler*) mgr->GetInputEventHandler();
  // UInt_t fSelectMask = inputHandler->IsEventSelected();
  //
  // if( fSelectMask& AliVEvent::kINT7 ) { fQAEventsTriggerSelection[iQAindex]->Fill("kINT7",1); }
  // else if (fSelectMask& AliVEvent::kHighMultV0) { fQAEventsTriggerSelection[iQAindex]->Fill("kHighMultV0",1); }
  // else if (fSelectMask& AliVEvent::kHighMultSPD) { fQAEventsTriggerSelection[iQAindex]->Fill("kHighMultSPD",1); }
  // else { fQAEventsTriggerSelection[iQAindex]->Fill("Other",1); }

  // SPD vertexer resolution
  Double_t cov[6] = {0};
  spdVtx->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  fhQAEventsSPDresol[iQAindex]->Fill(zRes);

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::Filtering()
{
  // main (envelope) method for filtering of all particles used for flow in selected events
  // All particles passing given requirements are pushed to correspodning std::vectors
  //  - in case of Charged & PID particles, only AliAODTrack* is stored (already loaded & stored by AliVEvent)
  //  - otherwise a new AliPicoTrack is constructed with relevant informations (needed to be deleted at destructor)
  // *************************************************************

  FilterCharged();

  fIndexSampling = GetSamplingIndex();
  fhEventSampling->Fill(fIndexCentrality,fIndexSampling);

  // // if neither is ON, filtering is skipped
  // if(!fProcessPID && !fProcessV0s && !fProcessPhi)
  // return;

  if(fProcessPID || fProcessPhi) { FilterPID(); }
  if(fProcessPhi) { FilterPhi(); }
  if(fProcessV0s) { FilterV0s(); }

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FilterCharged()
{
  // Filtering input charged tracks for POIs (stored in fVectorCharged) or RFPs (fVectorRefs)
  // If track passes all requirements its pointer is pushed to relevant vector container
  // *************************************************************

  Int_t iNumTracks = fEventAOD->GetNumberOfTracks();
  if(iNumTracks < 1) return;

  Int_t iNumRefs = 0;
  for(Int_t iTrack(0); iTrack < iNumTracks; iTrack++)
  {
    AliAODTrack* track = static_cast<AliAODTrack*>(fEventAOD->GetTrack(iTrack));
    if(!track) continue;

    if(fFillQA) FillQACharged(0,track); // QA before selection

    if(!IsChargedSelected(track)) continue;

    fVectorCharged->push_back(track);
    if(fFillQA) FillQACharged(1,track); // QA after selection

    if(fFlowFillWeights) { fh3WeightsCharged->Fill(track->Phi(),track->Eta(),fPVz); }
    if(fFlowUseWeights)
    {
      Double_t weight = 1.0;
      if(fFlowUse3Dweights){ weight = fh3WeightCharged->GetBinContent( fh3WeightCharged->FindFixBin(track->Eta(),track->Phi(),fPVz)); }
      else { weight = fh2WeightCharged->GetBinContent( fh2WeightCharged->FindFixBin(track->Eta(),track->Phi())); }
      fh3AfterWeightsCharged->Fill(track->Phi(),track->Eta(),fPVz,weight);
    }

    if(fMC)
    {
      AliAODMCParticle* trackMC = GetMCParticle(track->GetLabel());
      if(trackMC)
      {
        Int_t iPDG = TMath::Abs(trackMC->GetPdgCode());

        // filling info about all (i.e. before selection) reconstructed PID particles
        if(iPDG == 211) { fhMCRecoAllPionPt->Fill(track->Pt()); }
        if(iPDG == 321) { fhMCRecoAllKaonPt->Fill(track->Pt()); }
        if(iPDG == 2212) { fhMCRecoAllProtonPt->Fill(track->Pt()); }
      }
    }

    // Checking if selected track is eligible for Ref. flow
    if(!IsWithinRefs(track)) continue;

    // track is used for Ref. flow
    fVectorRefs->push_back(track);
    iNumRefs++;
    FillQARefs(1,track);
    if(fFlowFillWeights) { fh3WeightsRefs->Fill(track->Phi(),track->Eta(),fPVz); }
    if(fFlowUseWeights)
    {
      Double_t weight = 1.0;
      if(fFlowUse3Dweights) {weight = fh3WeightRefs->GetBinContent( fh3WeightRefs->FindFixBin(track->Eta(),track->Phi(),fPVz)); }
      else { weight = fh2WeightRefs->GetBinContent( fh2WeightRefs->FindFixBin(track->Eta(),track->Phi()) ); }
      fh3AfterWeightsRefs->Fill(track->Phi(),track->Eta(),fPVz,weight);
    }

  }

  // fill QA charged multiplicity
  fhRefsMult->Fill(iNumRefs);

  if(fFillQA)
  {
    fhQAChargedMult[0]->Fill(fEventAOD->GetNumberOfTracks());
    fhQAChargedMult[1]->Fill(fVectorCharged->size());
  }

  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::IsChargedSelected(const AliAODTrack* track)
{
  // Selection of charged track
  // returns kTRUE if track pass all requirements, kFALSE otherwise
  // *************************************************************
  if(!track) return kFALSE;
  fhChargedCounter->Fill("Input",1);

  // pt
  if(track->Pt() >= fFlowPOIsPtMax) return kFALSE;
  fhChargedCounter->Fill("Pt",1);

  // pseudorapidity (eta)
  if(fCutChargedEtaMax > 0. && TMath::Abs(track->Eta()) >= fCutChargedEtaMax) return kFALSE;
  fhChargedCounter->Fill("Eta",1);

  // filter bit
  if( !track->TestFilterBit(fCutChargedTrackFilterBit) ) return kFALSE;
  fhChargedCounter->Fill("FB",1);

  // number of TPC clusters (additional check for not ITS-standalone tracks)
  if( track->GetTPCNcls() < fCutChargedNumTPCclsMin && fCutChargedTrackFilterBit != 2) return kFALSE;
  fhChargedCounter->Fill("#TPC-Cls",1);

  // track DCA coordinates
  // note AliAODTrack::XYZAtDCA() works only for constrained tracks
  Double_t dTrackXYZ[3] = {0.,0.,0.};
  Double_t dVertexXYZ[3] = {0.,0.,0.};
  Double_t dDCAXYZ[3] = {0.,0.,0.};
  if( fCutChargedDCAzMax > 0. || fCutChargedDCAxyMax > 0.)
  {
    const AliAODVertex* vertex = fEventAOD->GetPrimaryVertex();
    if(!vertex) return kFALSE; // event does not have a PV

    track->GetXYZ(dTrackXYZ);
    vertex->GetXYZ(dVertexXYZ);

    for(Short_t i(0); i < 3; i++) { dDCAXYZ[i] = dTrackXYZ[i] - dVertexXYZ[i]; }
  }

  if(fCutChargedDCAzMax > 0. && TMath::Abs(dDCAXYZ[2]) > fCutChargedDCAzMax) return kFALSE;
  fhChargedCounter->Fill("DCA-z",1);

  if(fCutChargedDCAxyMax > 0. && TMath::Sqrt(dDCAXYZ[0]*dDCAXYZ[0] + dDCAXYZ[1]*dDCAXYZ[1]) > fCutChargedDCAxyMax) return kFALSE;
  fhChargedCounter->Fill("DCA-xy",1);


  // track passing all criteria
  fhChargedCounter->Fill("Selected",1);
  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::IsWithinRefs(const AliAODTrack* track)
{
  // Checking if (preselected) track fulfills criteria for RFPs
  // NOTE: This is not a standalone selection, but complementary check for IsChargedSelected()
  // It is used to selecting RFPs out of selected charged tracks
  // OR for estimating autocorrelations for Charged & PID particles
  // *************************************************************
  if(fCutFlowRFPsPtMin > 0.0 && track->Pt() < fCutFlowRFPsPtMin) return kFALSE;
  if(fCutFlowRFPsPtMax > 0.0 && track->Pt() > fCutFlowRFPsPtMax) return kFALSE;

  return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FillQARefs(const Short_t iQAindex, const AliAODTrack* track)
{
  // Filling various QA plots related to RFPs subset of charged track selection
  // *************************************************************

  if(!track) return;
  if(iQAindex == 0) return; // NOTE implemented only for selected RFPs

  fhRefsPt->Fill(track->Pt());
  fhRefsEta->Fill(track->Eta());
  fhRefsPhi->Fill(track->Phi());

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FillQACharged(const Short_t iQAindex, const AliAODTrack* track)
{
  // Filling various QA plots related to charged track selection
  // *************************************************************
  if(!track) return;

  // filter bit testing
  for(Short_t i(0); i < 32; i++)
  {
    if(track->TestFilterBit(TMath::Power(2.,i)))
      fhQAChargedFilterBit[iQAindex]->Fill(i);
  }

  // track charge
  fhQAChargedCharge[iQAindex]->Fill(track->Charge());

  // number of TPC clusters
  fhQAChargedNumTPCcls[iQAindex]->Fill(track->GetTPCNcls());

  // track DCA
  Double_t dDCAXYZ[3] = {-999., -999., -999.};
  const AliAODVertex* vertex = fEventAOD->GetPrimaryVertex();
  if(vertex)
  {
    Double_t dTrackXYZ[3] = {-999., -999., -999.};
    Double_t dVertexXYZ[3] = {-999., -999., -999.};

    track->GetXYZ(dTrackXYZ);
    vertex->GetXYZ(dVertexXYZ);

    for(Short_t i(0); i < 3; i++)
      dDCAXYZ[i] = dTrackXYZ[i] - dVertexXYZ[i];
  }
  fhQAChargedDCAxy[iQAindex]->Fill(TMath::Sqrt(dDCAXYZ[0]*dDCAXYZ[0] + dDCAXYZ[1]*dDCAXYZ[1]));
  fhQAChargedDCAz[iQAindex]->Fill(dDCAXYZ[2]);

  // kinematics
  fhQAChargedPt[iQAindex]->Fill(track->Pt());
  fhQAChargedPhi[iQAindex]->Fill(track->Phi());
  fhQAChargedEta[iQAindex]->Fill(track->Eta());

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FilterV0s()
{
  // Filtering input V0s candidates (K0s, (Anti)Lambda)
  // If track passes all requirements as defined in IsV0sSelected() (and species dependent one)
  // the relevant properties (pT, eta, phi,mass,species) are stored in a new AliPicoTrack
  // and pushed to relevant vector container.
  // *************************************************************

  Int_t iNumK0sSelected = 0;  // counter for selected K0s candidates
  Int_t iNumLambdaSelected = 0; // counter for selected Lambda candidates
  Int_t iNumALambdaSelected = 0; // counter for selected Anti-Lambda candidates

  Int_t iNumV0s = fEventAOD->GetNumberOfV0s();
  if(iNumV0s < 1) return;

  for(Int_t iV0(0); iV0 < iNumV0s; iV0++)
  {
    AliAODv0* v0 = static_cast<AliAODv0*>(fEventAOD->GetV0(iV0));
    if(!v0) continue;

    if(fFillQA) FillQAV0s(0,v0); // QA BEFORE selection

    if(IsV0Selected(v0))
    {
      Bool_t bIsK0s = IsV0aK0s(v0);
      Short_t iIsLambda = IsV0aLambda(v0);

      if(fFillQA && (bIsK0s || iIsLambda != 0))
      FillQAV0s(1,v0,bIsK0s,iIsLambda); // QA AFTER selection

      if(bIsK0s)
      {
        iNumK0sSelected++;
        fhV0sCounter->Fill("K^{0}_{S}",1);
        fhV0sInvMassK0s->Fill(v0->MassK0Short(),v0->MassLambda());

        fVectorK0s->push_back( new AliPicoTrack(v0->Pt(),v0->Eta(),v0->Phi(),v0->Charge(),0,0,0,0,0,0,v0->MassK0Short()) );

        if(fFlowFillWeights) { fh3WeightsK0s->Fill(v0->Phi(),v0->Eta(),fPVz); }
        if(fFlowUseWeights)
        {
          Double_t weight = 1.0;
          if(fFlowUse3Dweights) { weight = fh3WeightK0s->GetBinContent( fh3WeightK0s->FindFixBin(v0->Eta(),v0->Phi(),fPVz) ); }

          else { weight = fh2WeightK0s->GetBinContent( fh2WeightK0s->FindFixBin(v0->Eta(),v0->Phi()) ); }
          fh3AfterWeightsK0s->Fill(v0->Phi(),v0->Eta(),fPVz,weight);
        }
      }

      if(iIsLambda == 1) // lambda
      {
        iNumLambdaSelected++;
        fhV0sCounter->Fill("#Lambda/#bar{#Lambda}",1);
        fhV0sInvMassLambda->Fill(v0->MassK0Short(),v0->MassLambda());

        fVectorLambda->push_back( new AliPicoTrack(v0->Pt(),v0->Eta(),v0->Phi(),v0->Charge(),0,0,0,0,0,0,v0->MassLambda()) );

        if(fFlowFillWeights) { fh3WeightsLambda->Fill(v0->Phi(),v0->Eta(),fPVz); }
        if(fFlowUseWeights)
        {
          Double_t weight = 1.0;
          if(fFlowUse3Dweights) { weight = fh3WeightLambda->GetBinContent( fh3WeightLambda->FindFixBin(v0->Eta(),v0->Phi(),fPVz) ); }
          else { weight = fh2WeightLambda->GetBinContent( fh2WeightLambda->FindFixBin(v0->Eta(),v0->Phi()) ); }
          fh3AfterWeightsLambda->Fill(v0->Phi(),v0->Eta(),fPVz,weight);
        }
      }

      if(iIsLambda == -1) // anti-lambda
      {
        iNumALambdaSelected++;
        fhV0sCounter->Fill("#Lambda/#bar{#Lambda}",1);
        fhV0sInvMassLambda->Fill(v0->MassK0Short(),v0->MassAntiLambda());

        fVectorLambda->push_back( new AliPicoTrack(v0->Pt(),v0->Eta(),v0->Phi(),v0->Charge(),0,0,0,0,0,0,v0->MassAntiLambda()) );

        if(fFlowFillWeights) { fh3WeightsLambda->Fill(v0->Phi(),v0->Eta(),fPVz); }
        if(fFlowUseWeights)
        {
          Double_t weight = 1.0;
          if(fFlowUse3Dweights) { weight = fh3WeightLambda->GetBinContent( fh3WeightLambda->FindFixBin(v0->Eta(),v0->Phi(), fPVz) );}
          else { weight = fh2WeightLambda->GetBinContent( fh2WeightLambda->FindFixBin(v0->Eta(),v0->Phi()) ); }
          fh3AfterWeightsLambda->Fill(v0->Phi(),v0->Eta(),fPVz,weight);
        }
      }

      if(bIsK0s && iIsLambda != 0)
      fhV0sCounter->Fill("K^{0}_{S} && #Lambda/#bar{#Lambda}",1);
    }
  }

  // fill QA charged multiplicity
  if(fFillQA)
  {
    fhQAV0sMultK0s[0]->Fill(fEventAOD->GetNumberOfV0s());
    fhQAV0sMultLambda[0]->Fill(fEventAOD->GetNumberOfV0s());
    fhQAV0sMultALambda[0]->Fill(fEventAOD->GetNumberOfV0s());
    fhQAV0sMultK0s[1]->Fill(iNumK0sSelected);
    fhQAV0sMultLambda[1]->Fill(iNumLambdaSelected);
    fhQAV0sMultALambda[1]->Fill(iNumALambdaSelected);
  }

  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::IsV0aK0s(const AliAODv0* v0)
{
  // Topological reconstruction and selection of V0 candidates
  // specific for K0s candidates
  // return kTRUE if a candidate fulfill all requirements, kFALSE otherwise
  // *************************************************************
  if(!v0) return kFALSE;
  fhV0sCounterK0s->Fill("Input",1);

  // rapidity selection
  if(fCutV0sMotherRapMax > 0. && ( TMath::Abs(v0->RapK0Short()) > fCutV0sMotherRapMax ) ) return kFALSE;
  fhV0sCounterK0s->Fill("#it{y}",1);

  // inv. mass window
  Double_t dMass = v0->MassK0Short();
  if( dMass < fCutV0sInvMassK0sMin || dMass > fCutV0sInvMassK0sMax ) return kFALSE;
  fhV0sCounterK0s->Fill("InvMass",1);

  // cosine of pointing angle (CPA)
  if(fCutV0sCPAK0sMin > 0.)
  {
    Double_t dCPA = v0->CosPointingAngle(fEventAOD->GetPrimaryVertex());
    if( dCPA < fCutV0sCPAK0sMin ) return kFALSE;
  }
  fhV0sCounterK0s->Fill("CPA",1);

  // Armenteros-Podolaski plot
  if(fCutV0sArmenterosAlphaK0sMin > 0.)
  {
    Double_t dPtArm = v0->PtArmV0();
    Double_t dAlpha = v0->AlphaV0();
    if(dPtArm < (fCutV0sArmenterosAlphaK0sMin * TMath::Abs(dAlpha))) return kFALSE;
  }
  fhV0sCounterK0s->Fill("Armenteros-Podolanski",1);

  // proper life-time
  if( fCutV0sNumTauK0sMax > 0. )
  {
    Double_t dPrimVtxCoor[3] = {0.}; // primary vertex position {x,y,z}
    Double_t dSecVtxCoor[3] = {0.}; // secondary vertex position {x,y,z}
    Double_t dDecayCoor[3] = {0.}; // decay vector coor {xyz}
    AliAODVertex* primVtx = fEventAOD->GetPrimaryVertex();
    primVtx->GetXYZ(dPrimVtxCoor);
    v0->GetSecondaryVtx(dSecVtxCoor);

    for(Int_t i(0); i < 3; i++) { dDecayCoor[i] = dSecVtxCoor[i] - dPrimVtxCoor[i]; }

    // implementation in xy plane
    // Double_t dPropLife = ( (fPDGMassK0s / v0->Pt()) * TMath::Sqrt(dDecayCoor[0]*dDecayCoor[0] + dDecayCoor[1]*dDecayCoor[1]) );
    Double_t dPropLife = ( (fPDGMassK0s / (v0->P() + 1e-10) ) * TMath::Sqrt(dDecayCoor[0]*dDecayCoor[0] + dDecayCoor[1]*dDecayCoor[1] + dDecayCoor[2]*dDecayCoor[2]) );
    if( dPropLife > (fCutV0sNumTauK0sMax * 2.68) ) return kFALSE;
  }
  fhV0sCounterK0s->Fill("c#tau",1);

  // Daughter PID
  if(fCutV0sK0sPionNumTPCSigmaMax > 0.)
  {
    const AliAODTrack* daughterPos = (AliAODTrack*) v0->GetDaughter(0);
    const AliAODTrack* daughterNeg = (AliAODTrack*) v0->GetDaughter(1);

    if(!HasTrackPIDTPC(daughterPos) || !HasTrackPIDTPC(daughterNeg)) return kFALSE;

    if (daughterPos->GetTPCsignalN() < fCutV0sDaughterNumTPCClsPIDMin || daughterNeg->GetTPCsignalN() < fCutV0sDaughterNumTPCClsPIDMin) return kFALSE;
    Float_t nSigmaPiPos = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(daughterPos, AliPID::kPion));
    Float_t nSigmaPiNeg = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(daughterNeg, AliPID::kPion));
    if(nSigmaPiPos > fCutV0sK0sPionNumTPCSigmaMax || nSigmaPiNeg > fCutV0sK0sPionNumTPCSigmaMax) return kFALSE;
  }
  fhV0sCounterK0s->Fill("Daughters PID",1);

  // competing V0 rejection based on InvMass
  if(fCutV0sCrossMassRejection)
  {
    Double_t dMassLambda = v0->MassLambda();
    Double_t dMassALambda = v0->MassAntiLambda();

    // K0s candidate is within 10 MeV of (Anti)Lambda InvMass physSelTask
    if(TMath::Abs(dMassLambda - fPDGMassLambda) < fCutV0sCrossMassCutK0s)
    {
      // in Lambda peak
      fhV0sCompetingInvMassK0s->Fill(dMass,dMassLambda);
      return kFALSE;
    }

    if(TMath::Abs(dMassALambda - fPDGMassLambda) < fCutV0sCrossMassCutK0s)
    {
      // in Anti-Lambda peak
      fhV0sCompetingInvMassK0s->Fill(dMass,dMassALambda);
      return kFALSE;
    }
  }
  fhV0sCounterK0s->Fill("Competing InvMass",1);

  // passing all criteria
  fhV0sCounterK0s->Fill("Selected",1);
  return kTRUE;
}
//_____________________________________________________________________________
Short_t AliAnalysisTaskUniFlow::IsV0aLambda(const AliAODv0* v0)
{
  // Topological reconstruction and selection of V0 candidates
  // specific for Lambda candidates
  // return 0 if candidate does not fullfill any Lambda or Anti-Lambda requirements;
  // return 1 if a candidate fulfill all Lambda requirements;
  // return -1 if a candidate fullfill all Anti-Lambda requirements;
  // return 2 if a candidate fulfill all both Lambda & Anti-Lambda requirements
  // *************************************************************
  if(!v0) return 0;
  fhV0sCounterLambda->Fill("Input",1);

  // rapidity selection
  if(fCutV0sMotherRapMax > 0. && ( TMath::Abs(v0->RapLambda()) > fCutV0sMotherRapMax) ) return 0;
  fhV0sCounterLambda->Fill("#it{y}",1);

  // particle species dependent
  Bool_t bIsLambda = kFALSE;
  Bool_t bIsALambda = kFALSE;

  // inv. mass window
  Double_t dMassLambda = v0->MassLambda();
  Double_t dMassALambda = v0->MassAntiLambda();
  if( dMassLambda > fCutV0sInvMassLambdaMin && dMassLambda < fCutV0sInvMassLambdaMax) { bIsLambda = kTRUE; }
  if( dMassALambda > fCutV0sInvMassLambdaMin && dMassALambda < fCutV0sInvMassLambdaMax) { bIsALambda = kTRUE; }

  if(!bIsLambda && !bIsALambda) return 0;
  fhV0sCounterLambda->Fill("InvMass",1);

  // cosine of pointing angle (CPA)
  if( fCutV0sCPALambdaMin > 0. )
  {
    Double_t dCPA = v0->CosPointingAngle(fEventAOD->GetPrimaryVertex());
    if( dCPA < fCutV0sCPALambdaMin ) return 0;
  }
  fhV0sCounterLambda->Fill("CPA",1);

  // Armenteros-Podolaski plot
  if(fCutV0sArmenterosAlphaLambdaMax > 0.)
  {
    Double_t dPtArm = v0->PtArmV0();
    Double_t dAlpha = v0->AlphaV0();
    if(dPtArm > (fCutV0sArmenterosAlphaLambdaMax * TMath::Abs(dAlpha))) return 0;
  }
  fhV0sCounterLambda->Fill("Armenteros-Podolanski",1);

  // // Armenteros-Podolanski for candidates fullfilling both Lambda and Anti-Lambda selection
  // if(bIsLambda && bIsALambda)
  // {
  //   if(v0->AlphaV0() < 0.) bIsLambda = kFALSE;
  //   if(v0->AlphaV0() > 0.) bIsALambda = kFALSE;
  // }

  // proper life-time
  if( fCutV0sNumTauLambdaMax > 0. )
  {
    Double_t dPrimVtxCoor[3] = {0.}; // primary vertex position {x,y,z}
    Double_t dSecVtxCoor[3] = {0.}; // secondary vertex position {x,y,z}
    Double_t dDecayCoor[3] = {0.}; // decay vector coor {xyz}
    AliAODVertex* primVtx = fEventAOD->GetPrimaryVertex();
    primVtx->GetXYZ(dPrimVtxCoor);
    v0->GetSecondaryVtx(dSecVtxCoor);

    for(Int_t i(0); i < 3; i++) { dDecayCoor[i] = dSecVtxCoor[i] - dPrimVtxCoor[i]; }

    // implementation in xy plane
    // Double_t dPropLife = ( (fPDGMassLambda / v0->Pt()) * TMath::Sqrt(dDecayCoor[0]*dDecayCoor[0] + dDecayCoor[1]*dDecayCoor[1]) );
    Double_t dPropLife = ( (fPDGMassLambda / (v0->P() + 1e-10) ) * TMath::Sqrt(dDecayCoor[0]*dDecayCoor[0] + dDecayCoor[1]*dDecayCoor[1] + dDecayCoor[2]*dDecayCoor[2]) );
    if( dPropLife > (fCutV0sNumTauLambdaMax * 7.89) ) return 0;
  }
  fhV0sCounterLambda->Fill("c#tau",1);

  // daughter PID of Lambda Candidates
  if( fCutV0sLambdaProtonNumTPCSigmaMax > 0. || fCutV0sLambdaPionNumTPCSigmaMax > 0. )
  {
    const AliAODTrack* trackDaughterPos = (AliAODTrack*) v0->GetDaughter(0); // positive charge
    const AliAODTrack* trackDaughterNeg = (AliAODTrack*) v0->GetDaughter(1); // negative charge

    Bool_t bIsPosOK = HasTrackPIDTPC(trackDaughterPos);
    Bool_t bIsNegOK = HasTrackPIDTPC(trackDaughterNeg);

    Float_t dSigmaPos = 999.;
    Float_t dSigmaNeg = 999.;

    if(fCutV0sLambdaPionNumTPCSigmaMax > 0.) // check pions
    {
      if(bIsPosOK) dSigmaPos = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackDaughterPos, AliPID::kPion)); else dSigmaPos = 999.;
      if(bIsNegOK) dSigmaNeg = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackDaughterNeg, AliPID::kPion)); else dSigmaNeg = 999.;

      if(bIsLambda && (!bIsNegOK || dSigmaNeg > fCutV0sLambdaPionNumTPCSigmaMax) ) bIsLambda = kFALSE;
      if(bIsALambda && (!bIsPosOK || dSigmaPos > fCutV0sLambdaPionNumTPCSigmaMax) ) bIsALambda = kFALSE;
    }

    if(fCutV0sLambdaProtonNumTPCSigmaMax > 0.) // check protons
    {
      if(bIsPosOK) dSigmaPos = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackDaughterPos, AliPID::kProton)); else dSigmaPos = 999.;
      if(bIsNegOK) dSigmaNeg = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackDaughterNeg, AliPID::kProton)); else dSigmaNeg = 999.;

      if(bIsLambda && (!bIsPosOK || dSigmaPos > fCutV0sLambdaProtonNumTPCSigmaMax) ) bIsLambda = kFALSE;
      if(bIsALambda && (!bIsNegOK || dSigmaNeg > fCutV0sLambdaProtonNumTPCSigmaMax) ) bIsALambda = kFALSE;
    }

    if(!bIsLambda && !bIsALambda) return 0;
  }
  fhV0sCounterLambda->Fill("Daughter PID",1);

  // Lambda(AntiLamda) candidate is within fCutV0sCrossMassCutLambda of K0s InvMass
  if(fCutV0sCrossMassRejection)
  {
    Double_t dMassK0s = v0->MassK0Short();
    if(TMath::Abs(dMassK0s - fPDGMassK0s) < fCutV0sCrossMassCutLambda)
    {
      if(bIsLambda) fhV0sCompetingInvMassLambda->Fill(dMassK0s,dMassLambda);
      if(bIsALambda) fhV0sCompetingInvMassLambda->Fill(dMassK0s,dMassALambda);
      return 0;
    }
  }
  fhV0sCounterLambda->Fill("Competing InvMass",1);

  // passing all criteria
  fhV0sCounterLambda->Fill("Selected",1);

  if(bIsLambda && bIsALambda) { fhV0sCounterLambda->Fill("#Lambda && #bar{#Lambda}",1); return 2; } // both Lambda & Anti-Lambda
  if(bIsLambda) { fhV0sCounterLambda->Fill("only #Lambda",1); return 1; } // only Lambda
  if(bIsALambda) { fhV0sCounterLambda->Fill("only #bar{#Lambda}",1); return -1; } // only Anti-Lambda
  return 0; // neither
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskUniFlow::GetRapidity(Double_t mass, Double_t Pt, Double_t Eta)
{
    Double_t rapid = TMath::Log( (TMath::Sqrt(mass*mass + Pt*Pt*TMath::CosH(Eta)*TMath::CosH(Eta)) + Pt*TMath::SinH(Eta)) / TMath::Sqrt(mass*mass + Pt*Pt) );
    return rapid;
}
//_____________________________________________________________________________
AliAODMCParticle* AliAnalysisTaskUniFlow::GetMCParticle(Int_t label)
{
  if(!fArrayMC) { AliError("fArrayMC not found!"); return 0x0; }
  if(label < 0) { /*AliWarning("MC label negative");*/ return 0x0; }

  AliAODMCParticle* mcTrack = (AliAODMCParticle*) fArrayMC->At(label);
  if(!mcTrack) { AliError("Corresponding MC track not found!"); return 0x0; }
  return mcTrack;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::IsV0Selected(const AliAODv0* v0)
{
  // Topological reconstruction and selection of V0 candidates
  // common for both K0s and (Anti)-Lambdas
  // return kTRUE if a candidate fulfill all requirements, kFALSE otherwise
  // *************************************************************
  if(!v0) return kFALSE;
  fhV0sCounter->Fill("Input",1);

  const AliAODTrack* daughterPos = (AliAODTrack*) v0->GetDaughter(0);
  const AliAODTrack* daughterNeg = (AliAODTrack*) v0->GetDaughter(1);

  // daughter track check
  if(!daughterPos || !daughterNeg) return kFALSE;
  fhV0sCounter->Fill("Daughters OK",1);

  // acceptance checks
  if(fCutV0sMotherEtaMax > 0. && TMath::Abs(v0->Eta()) >= fCutV0sMotherEtaMax ) return kFALSE;
  if(v0->Pt() >= fFlowPOIsPtMax) return kFALSE;
  fhV0sCounter->Fill("Mother acceptance",1);

  if(fCutV0sDaughterPtMin > 0. && (daughterPos->Pt() <= fCutV0sDaughterPtMin  || daughterNeg->Pt() <= fCutV0sDaughterPtMin) ) return kFALSE;
  if(fCutV0sDaughterPtMax > 0. && (daughterPos->Pt() >= fCutV0sDaughterPtMax  || daughterNeg->Pt() >= fCutV0sDaughterPtMax) ) return kFALSE;
  if(fCutV0sDaughterEtaMax > 0. && ( (TMath::Abs(daughterNeg->Eta()) >= fCutV0sDaughterEtaMax) || (TMath::Abs(daughterPos->Eta()) >= fCutV0sDaughterEtaMax) ) ) return kFALSE;
  fhV0sCounter->Fill("Daughter acceptance",1);

  // daughters & mother charge checks
  if(v0->Charge() != 0) return kFALSE;
  if(daughterPos->Charge() == daughterNeg->Charge()) return kFALSE;
  if( (TMath::Abs(daughterPos->Charge()) != 1) || (TMath::Abs(daughterNeg->Charge()) != 1) ) return kFALSE;
  fhV0sCounter->Fill("Charge",1);

  // reconstruction method: online (on-the-fly) OR offline
  if(v0->GetOnFlyStatus() != fCutV0sOnFly) return kFALSE;
  fhV0sCounter->Fill("Reconstruction method",1);

  // TPC refit
  if(fCutV0srefitTPC && ( !daughterPos->IsOn(AliAODTrack::kTPCrefit) || !daughterNeg->IsOn(AliAODTrack::kTPCrefit) ) ) return kFALSE;
  fhV0sCounter->Fill("TPC refit",1);

  // filter bit
  if( fCutV0sDaughterFilterBit > 0 && (!daughterPos->TestFilterBit(fCutV0sDaughterFilterBit) || !daughterNeg->TestFilterBit(fCutV0sDaughterFilterBit) ) ) return kFALSE;
  fhV0sCounter->Fill("Daughter FB",1);

  // Kinks
  const AliAODVertex* prodVtxDaughterPos = (AliAODVertex*) daughterPos->GetProdVertex(); // production vertex of the positive daughter track
  const AliAODVertex* prodVtxDaughterNeg = (AliAODVertex*) daughterNeg->GetProdVertex(); // production vertex of the negative daughter track
  if(fCutV0srejectKinks && ( (prodVtxDaughterPos->GetType() == AliAODVertex::kKink ) || (prodVtxDaughterNeg->GetType() == AliAODVertex::kKink ) ) ) return kFALSE;
  fhV0sCounter->Fill("Kinks",1);

  // Daughter track quality
  if(daughterPos->GetTPCNcls() < fCutV0sDaughterNumTPCClsMin || daughterNeg->GetTPCNcls() < fCutV0sDaughterNumTPCClsMin) return kFALSE;
  if(daughterPos->GetTPCNCrossedRows() < fCutV0sDaughterNumTPCCrossMin || daughterNeg->GetTPCNCrossedRows() < fCutV0sDaughterNumTPCCrossMin) return kFALSE;
  if(daughterPos->GetTPCNclsF() < fCutV0sDaughterNumTPCFindMin || daughterNeg->GetTPCNclsF() < fCutV0sDaughterNumTPCFindMin) return kFALSE;
  if(fCutV0sDaughterRatioCrossFindMin > -1.)
  {
    if(daughterPos->GetTPCNclsF() < 1 || daughterNeg->GetTPCNclsF() < 1) return kFALSE; // at least 1 findable cls for proper division
    Double_t dRatioCrossFindPos = (Double_t) daughterPos->GetTPCNCrossedRows() / (Double_t) daughterPos->GetTPCNclsF();
    Double_t dRatioCrossFindNeg = (Double_t) daughterNeg->GetTPCNCrossedRows() / (Double_t) daughterNeg->GetTPCNclsF();
    if( dRatioCrossFindPos < fCutV0sDaughterRatioCrossFindMin || dRatioCrossFindNeg < fCutV0sDaughterRatioCrossFindMin) return kFALSE;
  }
  fhV0sCounter->Fill("Daughters track quality",1);

  // Daughters DCA to PV
  Double_t dDCAPosToPV = TMath::Abs(v0->DcaPosToPrimVertex());
  Double_t dDCANegToPV = TMath::Abs(v0->DcaNegToPrimVertex());
  if(fCutV0sDCAtoPVMin > 0. && ( dDCAPosToPV < fCutV0sDCAtoPVMin || dDCANegToPV < fCutV0sDCAtoPVMin ) ) return kFALSE;
  if(fCutV0sDCAtoPVMax > 0. && ( dDCAPosToPV > fCutV0sDCAtoPVMax || dDCANegToPV > fCutV0sDCAtoPVMax ) ) return kFALSE;

  // note AliAODTrack::XYZAtDCA() works only for constrained tracks
  if(fCutV0sDCAtoPVzMax > 0.)
  {
    Double_t dVertexXYZ[3] = {0.};
    Double_t dTrackXYZpos[3] = {0.};
    Double_t dTrackXYZneg[3] = {0.};
    Double_t dDCAXYZpos[3] = {0.};
    Double_t dDCAXYZneg[3] = {0.};

    const AliAODVertex* vertex = fEventAOD->GetPrimaryVertex();
    if(!vertex) return kFALSE; // event does not have a PV

    vertex->GetXYZ(dVertexXYZ);
    daughterPos->GetXYZ(dTrackXYZpos);
    daughterNeg->GetXYZ(dTrackXYZneg);

    for(Short_t i(0); i < 3; i++)
    {
      dDCAXYZpos[i] = dTrackXYZpos[i] - dVertexXYZ[i];
      dDCAXYZneg[i] = dTrackXYZneg[i] - dVertexXYZ[i];
    }

    if( TMath::Abs(dDCAXYZpos[2]) > fCutV0sDCAtoPVzMax || TMath::Abs(dDCAXYZneg[2]) > fCutV0sDCAtoPVzMax ) return kFALSE;
  }
  fhV0sCounter->Fill("DCA to PV",1);

  // Daughter DCA among themselves
  Double_t dDCADaughters = v0->DcaV0Daughters();
  if(fCutV0sDCADaughtersMin > 0. && TMath::Abs(dDCADaughters) < fCutV0sDCADaughtersMin) return kFALSE;
  if(fCutV0sDCADaughtersMax > 0. && TMath::Abs(dDCADaughters) > fCutV0sDCADaughtersMax) return kFALSE;
  fhV0sCounter->Fill("Daughters DCA",1);

  // radius of decay vertex in transverse plane
  Double_t dDecayRadius = v0->RadiusV0();
  if( fCutV0sDecayRadiusMin > 0. && (dDecayRadius < fCutV0sDecayRadiusMin) ) return kFALSE;
  if( fCutV0sDecayRadiusMax > 0. && (dDecayRadius > fCutV0sDecayRadiusMax) ) return kFALSE;
  fhV0sCounter->Fill("Decay radius",1);

  // passing all common criteria
  fhV0sCounter->Fill("Common passed",1);
  return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FillQAV0s(const Short_t iQAindex, const AliAODv0* v0, const Bool_t bIsK0s, const Short_t bIsLambda)
{
  // Filling various QA plots related to V0 candidate selection
  // *************************************************************
  // checking mother & daughters
  if(!v0) return;
  AliAODTrack* trackDaughter[2] = {(AliAODTrack*) v0->GetDaughter(0), (AliAODTrack*) v0->GetDaughter(1)};
  if(!trackDaughter[0] || !trackDaughter[1]) return;

  // setting internal flags for Lambdas and Anti-Lambdas
  Bool_t bCandLambda;
  Bool_t bCandAntiLambda;

  switch (bIsLambda)
  {
    case 1:
    {
      bCandLambda = kTRUE;
      bCandAntiLambda = kFALSE;
      break;
    }
    case -1:
    {
      bCandLambda = kFALSE;
      bCandAntiLambda = kTRUE;
      break;
    }
    case 2:
    {
      bCandLambda = kTRUE;
      bCandAntiLambda = kTRUE;
      break;
    }
    default:
    {
      bCandLambda = kFALSE;
      bCandAntiLambda = kFALSE;
    }
  }

  // reconstruction method
  fhQAV0sRecoMethod[iQAindex]->Fill(v0->GetOnFlyStatus());

  // DCA between daughters and PV
  fhQAV0sDCAtoPV[iQAindex]->Fill(v0->DcaPosToPrimVertex());
  fhQAV0sDCAtoPV[iQAindex]->Fill(v0->DcaNegToPrimVertex());

  // Daughter DCA among themselves
  fhQAV0sDCADaughters[iQAindex]->Fill(v0->DcaV0Daughters());

  // charge
  fhQAV0sMotherCharge[iQAindex]->Fill(v0->Charge());

  // radius of decay vertex in transverse plane
  Double_t dSecVtxCoor[3] = {0};
  v0->GetSecondaryVtx(dSecVtxCoor);
  Double_t dDecayRadius = TMath::Sqrt(dSecVtxCoor[0]*dSecVtxCoor[0] + dSecVtxCoor[1]*dSecVtxCoor[1]);
  fhQAV0sDecayRadius[iQAindex]->Fill(dDecayRadius);

  // mother kinematics
  fhQAV0sMotherPt[iQAindex]->Fill(v0->Pt());
  fhQAV0sMotherPhi[iQAindex]->Fill(v0->Phi());
  fhQAV0sMotherEta[iQAindex]->Fill(v0->Eta());

  // proper lifetime preparation (to be filled in particle dependent if scope)
  Double_t dPrimVtxCoor[3] = {0};
  Double_t dDecayCoor[3] = {0};
  AliAODVertex* primVtx2 = fEventAOD->GetPrimaryVertex();
  primVtx2->GetXYZ(dPrimVtxCoor);
  for(Int_t i(0); i < 2; i++)
    dDecayCoor[i] = dSecVtxCoor[i] - dPrimVtxCoor[i];

  // particle dependent
  if(bIsK0s)
  {
    // K0s
    fhQAV0sMotherRapK0s[iQAindex]->Fill(v0->RapK0Short());
    fhQAV0sInvMassK0s[iQAindex]->Fill(v0->MassK0Short());

    // CPA
    AliAODVertex* primVtx = fEventAOD->GetPrimaryVertex();
    fhQAV0sCPAK0s[iQAindex]->Fill(v0->CosPointingAngle(primVtx));

    // Armenteros-Podolanski
    fhQAV0sArmenterosK0s[iQAindex]->Fill(v0->AlphaV0(), v0->PtArmV0());

    // proper lifetime
    Double_t dMassPDGK0s = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();
    Double_t dPropLifeK0s = ( (dMassPDGK0s / v0->Pt()) * TMath::Sqrt(dDecayCoor[0]*dDecayCoor[0] + dDecayCoor[1]*dDecayCoor[1]) );
    fhQAV0sNumTauK0s[iQAindex]->Fill(dPropLifeK0s);
  }
  if(bCandLambda || bCandAntiLambda)
  {
    // (Anti)Lambda
    fhQAV0sMotherRapLambda[iQAindex]->Fill(v0->RapLambda());
    fhQAV0sInvMassLambda[iQAindex]->Fill(v0->MassLambda());
    fhQAV0sInvMassLambda[iQAindex]->Fill(v0->MassAntiLambda());

    // CPA
    AliAODVertex* primVtx = fEventAOD->GetPrimaryVertex();
    fhQAV0sCPALambda[iQAindex]->Fill(v0->CosPointingAngle(primVtx));

    // Armenteros-Podolanski
    if(bCandLambda)
      fhQAV0sArmenterosLambda[iQAindex]->Fill(v0->AlphaV0(), v0->PtArmV0());

    if(bCandAntiLambda)
      fhQAV0sArmenterosALambda[iQAindex]->Fill(v0->AlphaV0(), v0->PtArmV0());

    // proper lifetime
    Double_t dMassPDGLambda = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();
    Double_t dPropLifeLambda = ( (dMassPDGLambda / v0->Pt()) * TMath::Sqrt(dDecayCoor[0]*dDecayCoor[0] + dDecayCoor[1]*dDecayCoor[1]) );
    fhQAV0sNumTauLambda[iQAindex]->Fill(dPropLifeLambda);
  }

  AliPIDResponse::EDetPidStatus pidStatusTPC;
  AliPIDResponse::EDetPidStatus pidStatusTOF;
  UShort_t numTPCfindable = 0;
  Float_t numTPCcrossed = 0;

  // daughters properties
  AliAODVertex* prodVtxDaughter = 0x0;
  for(Short_t i(0); i < 2; i++)
  {
    // TPC refit
    fhQAV0sDaughterTPCRefit[iQAindex]->Fill(trackDaughter[i]->IsOn(AliAODTrack::kTPCrefit));

    // kinks
    prodVtxDaughter = (AliAODVertex*) trackDaughter[i]->GetProdVertex();
    fhQAV0sDaughterKinks[iQAindex]->Fill(prodVtxDaughter->GetType() == AliAODVertex::kKink);

    // track quality
    numTPCcrossed = trackDaughter[i]->GetTPCNCrossedRows();
    numTPCfindable = trackDaughter[i]->GetTPCNclsF();
    fhQAV0sDaughterNumTPCCls[iQAindex]->Fill(trackDaughter[i]->GetTPCNcls());
    fhQAV0sDaughterNumTPCClsPID[iQAindex]->Fill(trackDaughter[i]->GetTPCsignalN());
    fhQAV0sDaughterNumTPCFind[iQAindex]->Fill(numTPCfindable);
    fhQAV0sDaughterNumTPCCrossRows[iQAindex]->Fill(numTPCcrossed);
    if(numTPCfindable > 0.) fhQAV0sDaughterTPCCrossFindRatio[iQAindex]->Fill(numTPCcrossed/numTPCfindable); else fhQAV0sDaughterTPCCrossFindRatio[iQAindex]->Fill(-5.);

    // detector status
    pidStatusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, trackDaughter[i]);
    pidStatusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, trackDaughter[i]);
    fhQAV0sDaughterTPCstatus[iQAindex]->Fill((Int_t) pidStatusTPC );
    fhQAV0sDaughterTOFstatus[iQAindex]->Fill((Int_t) pidStatusTOF );

    // daughter kinematics
    fhQAV0sDaughterPt[iQAindex]->Fill(trackDaughter[i]->Pt());
    fhQAV0sDaughterPhi[iQAindex]->Fill(trackDaughter[i]->Phi());
    fhQAV0sDaughterEta[iQAindex]->Fill(trackDaughter[i]->Eta());

    // daughter charge
    fhQAV0sDaughterCharge[iQAindex]->Fill(trackDaughter[i]->Charge());
  }

  AliPIDResponse::EDetPidStatus pidStatusTPCpos;
  AliPIDResponse::EDetPidStatus pidStatusTPCneg;

  // PID checks
  if(fPIDResponse)
  {
    // checking the detector status
    pidStatusTPCpos = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, trackDaughter[0]);
    pidStatusTPCneg = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, trackDaughter[1]);

    if(pidStatusTPCpos == AliPIDResponse::kDetPidOk && pidStatusTPCneg == AliPIDResponse::kDetPidOk)
    {
      if(bIsK0s)
      {
        // daughter PID
        fhQAV0sDaughterTPCdEdxK0s[iQAindex]->Fill(trackDaughter[0]->P(), trackDaughter[0]->GetTPCsignal());
        fhQAV0sDaughterTPCdEdxK0s[iQAindex]->Fill(trackDaughter[1]->P(), trackDaughter[1]->GetTPCsignal());

        // Pion PID for daughters
        fhQAV0sDaughterNumSigmaPionK0s[iQAindex]->Fill(v0->Pt(), fPIDResponse->NumberOfSigmasTPC(trackDaughter[0], AliPID::kPion));
        fhQAV0sDaughterNumSigmaPionK0s[iQAindex]->Fill(v0->Pt(), fPIDResponse->NumberOfSigmasTPC(trackDaughter[1], AliPID::kPion));
      }

      if(bCandLambda || bCandAntiLambda)
      {
        // daughter PID
        fhQAV0sDaughterTPCdEdxLambda[iQAindex]->Fill(trackDaughter[0]->P(), trackDaughter[0]->GetTPCsignal());
        fhQAV0sDaughterTPCdEdxLambda[iQAindex]->Fill(trackDaughter[1]->P(), trackDaughter[1]->GetTPCsignal());

        if(bCandLambda)
        {
          fhQAV0sDaughterNumSigmaProtonLambda[iQAindex]->Fill(v0->Pt(), fPIDResponse->NumberOfSigmasTPC(trackDaughter[0], AliPID::kProton));
          fhQAV0sDaughterNumSigmaPionLambda[iQAindex]->Fill(v0->Pt(), fPIDResponse->NumberOfSigmasTPC(trackDaughter[1], AliPID::kPion));
        }

        if(bCandAntiLambda)
        {
          fhQAV0sDaughterNumSigmaProtonALambda[iQAindex]->Fill(v0->Pt(), fPIDResponse->NumberOfSigmasTPC(trackDaughter[1], AliPID::kProton));
          fhQAV0sDaughterNumSigmaPionALambda[iQAindex]->Fill(v0->Pt(), fPIDResponse->NumberOfSigmasTPC(trackDaughter[0], AliPID::kPion));
        }
      }
    }
  }

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FilterPhi()
{
  // Reconstruction and filtering of Phi meson candidates out of selected Kaon sample
  // If track passes all requirements, the relevant properties (pT, eta, phi) are stored
  // in FlowPart struct  and pushed to relevant vector container.
  // *************************************************************
  Int_t iNumKaons = (Int_t) fVectorKaon->size();
  // check if there are at least 2 selected kaons in event (minimum for phi reconstruction)
  if(iNumKaons < 2) return;

  // start Phi reconstruction
  Int_t iNumBG = 0;
  for(Int_t iKaon1(0); iKaon1 < iNumKaons; iKaon1++)
  {
    AliAODTrack* kaon1 = dynamic_cast<AliAODTrack*>(fVectorKaon->at(iKaon1));
    if(!kaon1) continue;

    for(Int_t iKaon2(iKaon1+1); iKaon2 < iNumKaons; iKaon2++)
    {
      AliAODTrack* kaon2 = dynamic_cast<AliAODTrack*>(fVectorKaon->at(iKaon2));
      if(!kaon2) continue;

      AliPicoTrack* mother = MakeMother(kaon1,kaon2);
      fhPhiCounter->Fill("Input",1);

      // filling QA BEFORE selection
      if(fFillQA) FillQAPhi(0,mother);

      if(mother->M() < fCutPhiInvMassMin || mother->M() > fCutPhiInvMassMax) { delete mother; continue; }
      fhPhiCounter->Fill("InvMass",1);

      if(mother->Pt() < fFlowPOIsPtMin || mother->Pt() > fFlowPOIsPtMax) { delete mother; continue; }
      fhPhiCounter->Fill("Pt",1);

      if(fCutPhiMotherEtaMax > 0 && TMath::Abs(mother->Eta()) > fCutPhiMotherEtaMax) { delete mother; continue; }
      fhPhiCounter->Fill("Eta",1);

      // mother (phi) candidate passing all criteria
      fhPhiCounter->Fill("Selected",1);

      // filling QA AFTER selection
      if(fFillQA) FillQAPhi(1,mother);

      // filling weights
      if(fFlowFillWeights) { fh3WeightsPhi->Fill(mother->Phi(), mother->Eta(), fPVz); }
      if(fFlowUseWeights)
      {
        Double_t weight = 1.0;
        if(fFlowUse3Dweights) { weight = fh3WeightPhi->GetBinContent( fh3WeightPhi->FindFixBin(mother->Eta(),mother->Phi(), fPVz) ); }
        else { weight = fh2WeightPhi->GetBinContent( fh2WeightPhi->FindFixBin(mother->Eta(),mother->Phi()) ) ; }
        fh3AfterWeightsPhi->Fill(mother->Phi(),mother->Eta(),fPVz,weight);
      }

      if(mother->Charge() == 0)
      {
        // opposite-sign combination (signal+background)
        fhPhiCounter->Fill("Unlike-sign",1);
        fVectorPhi->push_back(mother);
      }

      if(TMath::Abs(mother->Charge()) == 2)
      {
        // like-sign combination (background)
        fhPhiCounter->Fill("BG",1);
        iNumBG++;

        // filing background entries for Phi candidates
        for(Int_t iGap(0); iGap < fNumEtaGap; iGap++)
        {
          if(mother->Eta() > fEtaGap[iGap]/2 ) fh3PhiEntriesBGPos[iGap]->Fill(fIndexCentrality,mother->Pt(),mother->M());
          if(fEtaGap[iGap] > -1. && mother->Eta() < -fEtaGap[iGap]/2 ) fh3PhiEntriesBGNeg[iGap]->Fill(fIndexCentrality,mother->Pt(),mother->M());
        }
      }

    } // endfor {iKaon2} : second kaon
  } // endfor {iKaon1} : first Kaon

  // filling multiplicity distribution
  fhPhiMult->Fill(fVectorPhi->size());
  fhPhiBGMult->Fill(iNumBG);

  return;
}
//____________________________________________________________________
AliPicoTrack* AliAnalysisTaskUniFlow::MakeMother(const AliAODTrack* part1, const AliAODTrack* part2)
{
  // Reconstructing mother particle from two prongs and fill its properties.
  // return ptr to created mother particle
  // *************************************************************

  if(!part1 || !part2) return 0x0;

  // combining momenta
  TVector3 mom1 = TVector3( part1->Px(), part1->Py(), part1->Pz() );
  TVector3 mom2 = TVector3( part2->Px(), part2->Py(), part2->Pz() );
  TVector3 mom = mom1 + mom2;

  Byte_t iCharge = part1->Charge() + part2->Charge();

  // calculating inv. mass
  Double_t dMass = -999.;
  Double_t dE1 = TMath::Sqrt( mom1.Mag2() + TMath::Power(fPDGMassKaon,2) );
  Double_t dE2 = TMath::Sqrt( mom2.Mag2() + TMath::Power(fPDGMassKaon,2) );

  Double_t dMassSq = TMath::Power((dE1+dE2),2) - mom.Mag2();
  if(dMassSq >= 0.) dMass = TMath::Sqrt(dMassSq);

  return new AliPicoTrack(mom.Pt(),mom.Eta(),mom.Phi(),iCharge,0,0,0,0,0,0,dMass);
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FillQAPhi(const Short_t iQAindex, const AliPicoTrack* part)
{
  if(!part) return;

  if(iQAindex == 0) return; // TODO not implemented (do not know what)

  if(part->Charge() == 0)
  {
    // phi candidate (unlike-sign pair)
    fhPhiInvMass->Fill(part->M());
    fhPhiCharge->Fill(part->Charge());
    fhPhiPt->Fill(part->Pt());
    fhPhiEta->Fill(part->Eta());
    fhPhiPhi->Fill(part->Phi());
  }

  if(TMath::Abs(part->Charge()) == 2)
  {
    // phi candidate (unlike-sign pair)
    fhPhiBGInvMass->Fill(part->M());
    fhPhiBGCharge->Fill(part->Charge());
  }

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FilterPID()
{
  // Filtering input PID tracks (pi,K,p)
  // If track passes all requirements as defined in IsPIDSelected() (and species dependent),
  // the relevant properties (pT, eta, phi) are stored in FlowPart struct
  // and pushed to relevant vector container.
  // return kFALSE if any complications occurs
  // *************************************************************

  for(auto part = fVectorCharged->begin(); part != fVectorCharged->end(); ++part)
  {
    AliAODTrack* track = static_cast<AliAODTrack*>(*part);
    if(!track) continue;

    // PID tracks are subset of selected charged tracks (same quality requirements)
    // if(!IsChargedSelected(track)) continue;

    if(fFillQA) FillQAPID(0,track,kUnknown);   // filling QA for tracks before selection (but after charged criteria applied)

    // PID track selection (return most favourable species)
    PartSpecies species = IsPIDSelected(track);
    // check if only protons should be used
    if(fCutPIDUseAntiProtonOnly && species == kProton && track->Charge() == 1) { species = kUnknown; }

    // selection of PID tracks
    switch (species)
    {
      case kPion:
        fVectorPion->push_back(track);
        if(fFlowFillWeights) { fh3WeightsPion->Fill(track->Phi(), track->Eta(), fPVz); }
        if(fFlowUseWeights)
        {
          Double_t weight = 1.0;
          if(fFlowUse3Dweights) { weight = fh3WeightPion->GetBinContent( fh3WeightPion->FindFixBin(track->Eta(),track->Phi(), fPVz)); }
          else { weight = fh2WeightPion->GetBinContent( fh2WeightPion->FindFixBin(track->Eta(),track->Phi()) ); }
          fh3AfterWeightsPion->Fill(track->Phi(),track->Eta(),fPVz,weight);
        }
        break;
      case kKaon:
        fVectorKaon->push_back(track);
        if(fFlowFillWeights) { fh3WeightsKaon->Fill(track->Phi(), track->Eta(), fPVz); }
        if(fFlowUseWeights)
        {
          Double_t weight = 1.0;
          if(fFlowUse3Dweights) { weight = fh3WeightKaon->GetBinContent( fh3WeightKaon->FindFixBin(track->Eta(),track->Phi(),fPVz) ); }
          else { weight = fh2WeightKaon->GetBinContent( fh2WeightKaon->FindFixBin(track->Eta(),track->Phi()) ); }
          fh3AfterWeightsKaon->Fill(track->Phi(),track->Eta(),fPVz,weight);
        }
        break;
      case kProton:
        fVectorProton->push_back(track);
        if(fFlowFillWeights) { fh3WeightsProton->Fill(track->Phi(), track->Eta(), fPVz); }
        if(fFlowUseWeights)
        {
          Double_t weight = 1.0;
          if(fFlowUse3Dweights) { weight = fh3WeightProton->GetBinContent( fh3WeightProton->FindFixBin(track->Eta(),track->Phi(),fPVz) );}
          else { weight = fh2WeightProton->GetBinContent( fh2WeightProton->FindFixBin(track->Eta(),track->Phi()) ); }
          fh3AfterWeightsProton->Fill(track->Phi(),track->Eta(),fPVz,weight);
        }
        break;
      default:
        continue;
    }

    if(fFillQA) FillQAPID(1,track,species); // filling QA for tracks AFTER selection

    // process MC data
    if(fMC)
    {
      AliAODMCParticle* trackMC = GetMCParticle(track->GetLabel());
      if(trackMC)
      {
        Int_t iPDG = TMath::Abs(trackMC->GetPdgCode());

        switch (species)
        {
          case kPion:
            fhMCRecoSelectedPionPt->Fill(track->Pt());
            if(iPDG == 211) { fhMCRecoSelectedTruePionPt->Fill(track->Pt()); }
          break;
          case kKaon:
            fhMCRecoSelectedKaonPt->Fill(track->Pt());
            if(iPDG == 321) { fhMCRecoSelectedTrueKaonPt->Fill(track->Pt()); }
          break;
          case kProton:
            fhMCRecoSelectedProtonPt->Fill(track->Pt());
            if(iPDG == 2212) { fhMCRecoSelectedTrueProtonPt->Fill(track->Pt()); }
          break;

          default :
            continue;
        }

      }
    }

  }

  fhPIDPionMult->Fill(fVectorPion->size());
  fhPIDKaonMult->Fill(fVectorKaon->size());
  fhPIDProtonMult->Fill(fVectorProton->size());

  return;
}
//_____________________________________________________________________________
AliAnalysisTaskUniFlow::PartSpecies AliAnalysisTaskUniFlow::IsPIDSelected(const AliAODTrack* track)
{
  // Selection of PID tracks (pi,K,p) - track identification
  // Based on fCutUseBayesPID flag, either Bayes PID or nSigma cutting is used
  // returns AliAnalysisTaskUniFlow::PartSpecies enum : kPion, kKaon, kProton if any of this passed kUnknown otherwise
  // *************************************************************

  // checking detector statuses
  Bool_t bIsTPCok = HasTrackPIDTPC(track);
  Bool_t bIsTOFok = HasTrackPIDTOF(track);

  if(!bIsTPCok) return kUnknown;
  // TODO: TOF check???

  if(fCutUseBayesPID)
  {
    // use Bayesian PID
    Double_t dProbPID[5] = {0.}; // array for Bayes PID probabilities:  0: electron / 1: muon / 2: pion / 3: kaon / 4: proton
    UInt_t iDetUsed = fPIDCombined->ComputeProbabilities(track, fPIDResponse, dProbPID); // filling probabilities to dPropPID array
		if(iDetUsed > 999) {}
    Double_t dMaxProb = TMath::MaxElement(5,dProbPID);

    // printf("PID Prob: e %g | mu %g | pi %g | K %g | p %g ||| MAX %g \n",dProbPID[0],dProbPID[1],dProbPID[2],dProbPID[3],dProbPID[4],dMaxProb);

    // electron and mion rejection
    if(dProbPID[0] >= fCutPIDBayesRejectElectron || dProbPID[1] >= fCutPIDBayesRejectMuon) return kUnknown;

    // checking the PID probability
    // TODO: think about: if Pion has maximum probibility < fCutBayesPIDPion, track is rejected -> is it good?
    if(dMaxProb == dProbPID[2] && dProbPID[2] >= fCutPIDBayesPionMin) return kPion;
    if(dMaxProb == dProbPID[3] && dProbPID[3] >= fCutPIDBayesKaonMin) return kKaon;
    if(dMaxProb == dProbPID[4] && dProbPID[4] >= fCutPIDBayesProtonMin) return kProton;
  }
  else
  {
    // use nSigma cuts (based on combination of TPC / TOF nSigma cuts
    const Double_t dPt = track->Pt();

    Float_t dNumSigmaTPC[5] = {-99.,-99.,-99.,-99.,-99.}; // TPC nSigma array: 0: electron / 1: muon / 2: pion / 3: kaon / 4: proton
    Float_t dNumSigmaTOF[5] = {-99.,-99.,-99.,-99.,-99.}; // TOF nSigma array: 0: electron / 1: muon / 2: pion / 3: kaon / 4: proton

    // filling nSigma arrays
    if(bIsTPCok) // should be anyway
    {
      dNumSigmaTPC[0] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron));
      dNumSigmaTPC[1] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kMuon));
      dNumSigmaTPC[2] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion));
      dNumSigmaTPC[3] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon));
      dNumSigmaTPC[4] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton));
    }

    if(bIsTOFok) // should be anyway
    {
      dNumSigmaTOF[0] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kElectron));
      dNumSigmaTOF[1] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kMuon));
      dNumSigmaTOF[2] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion));
      dNumSigmaTOF[3] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon));
      dNumSigmaTOF[4] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton));
    }

    // TPC nSigma cuts
    if(dPt <= 0.4)
    {
      Float_t dMinSigmasTPC = TMath::MinElement(5,dNumSigmaTPC);

      // electron rejection
      if(dMinSigmasTPC == dNumSigmaTPC[0] && dNumSigmaTPC[0] <= fCutPIDnSigmaTPCRejectElectron) return kUnknown;
      if(dMinSigmasTPC == dNumSigmaTPC[2] && dNumSigmaTPC[2] <= fCutPIDnSigmaPionMax) return kPion;
      if(dMinSigmasTPC == dNumSigmaTPC[3] && dNumSigmaTPC[3] <= fCutPIDnSigmaKaonMax) return kKaon;
      if(dMinSigmasTPC == dNumSigmaTPC[4] && dNumSigmaTPC[4] <= fCutPIDnSigmaProtonMax) return kProton;
    }

    // combined TPC + TOF nSigma cuts
    // NB: for testing of efficiency removed the upper limmit for PID
    // if(dPt > 0.4 && dPt < 4.0)
    if(dPt > 0.4)
    {
      Float_t dNumSigmaCombined[5] = {-99.,-99.,-99.,-99.,-99.};

      // discard candidates if no TOF is available if cut is on
      if(fCutPIDnSigmaCombinedTOFrejection && !bIsTOFok) return kUnknown;

      // calculating combined nSigmas
      for(Short_t i(0); i < 5; i++)
      {
        if(bIsTOFok) { dNumSigmaCombined[i] = TMath::Sqrt(dNumSigmaTPC[i]*dNumSigmaTPC[i] + dNumSigmaTOF[i]*dNumSigmaTOF[i]); }
        else { dNumSigmaCombined[i] = dNumSigmaTPC[i]; }
      }

      Float_t dMinSigmasCombined = TMath::MinElement(5,dNumSigmaCombined);

      // electron rejection
      if(dMinSigmasCombined == dNumSigmaCombined[0] && dNumSigmaCombined[0] <= fCutPIDnSigmaTPCRejectElectron) return kUnknown;
      if(dMinSigmasCombined == dNumSigmaCombined[2] && dNumSigmaCombined[2] <= fCutPIDnSigmaPionMax) return kPion;
      if(dMinSigmasCombined == dNumSigmaCombined[3] && dNumSigmaCombined[3] <= fCutPIDnSigmaKaonMax) return kKaon;
      if(dMinSigmasCombined == dNumSigmaCombined[4] && dNumSigmaCombined[4] <= fCutPIDnSigmaProtonMax) return kProton;
    }

    // if(dPt >= 4.0)
    // {
      // NOTE: in this pt range, nSigmaTPC is not enought to distinquish well between species
      // all three values are close to each other -> minimum difference is not applied, just nSigma cut

      // TPC dEdx parametrisation (dEdx - <dEdx>)
      // TODO: TPC dEdx parametrisation cuts
      // if(dPt > 3.)
      //
    // }
  }

  return kUnknown;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FillQAPID(const Short_t iQAindex, const AliAODTrack* track, const PartSpecies species)
{
  // Filling various QA plots related to PID (pi,K,p) track selection
  // *************************************************************
  if(!track) return;

  if(!fPIDResponse || !fPIDCombined)
  {
    AliError("AliPIDResponse or AliPIDCombined object not found!");
    return;
  }

  // TPC & TOF statuses & measures
  AliPIDResponse::EDetPidStatus pidStatusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, track);
  AliPIDResponse::EDetPidStatus pidStatusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, track);

  fhQAPIDTOFstatus[iQAindex]->Fill((Int_t) pidStatusTOF );
  fhQAPIDTPCstatus[iQAindex]->Fill((Int_t) pidStatusTPC );

  Bool_t bIsTPCok = HasTrackPIDTPC(track);
  Bool_t bIsTOFok = HasTrackPIDTOF(track);

  Double_t dNumSigmaTPC[5] = {-11}; // array: 0: electron / 1: muon / 2: pion / 3: kaon / 4: proton
  Double_t dNumSigmaTOF[5] = {-11}; // array: 0: electron / 1: muon / 2: pion / 3: kaon / 4: proton
  Double_t dBayesProb[5] = {-0.1}; // Bayesian probability | array: 0: electron / 1: muon / 2: pion / 3: kaon / 4: proton

  Double_t dTPCdEdx = -5; // TPC dEdx for selected particle
  Double_t dTOFbeta = -0.05; //TOF beta for selected particle

  // filling Bayesian PID probabilities to dBayesProb array
  UInt_t iDetUsed = fPIDCombined->ComputeProbabilities(track, fPIDResponse, dBayesProb);
	if(iDetUsed > 999) {}

  Double_t dP = track->P();
  Double_t dPt = track->Pt();

  // detector status dependent
  if(bIsTPCok)
  {
    dNumSigmaTPC[0] = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
    dNumSigmaTPC[1] = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kMuon);
    dNumSigmaTPC[2] = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
    dNumSigmaTPC[3] = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
    dNumSigmaTPC[4] = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);

    dTPCdEdx = track->GetTPCsignal();
    fhQAPIDTPCdEdx[iQAindex]->Fill(track->P(), dTPCdEdx);
  }
  else // TPC status not OK
  {
    dNumSigmaTPC[0] = -11.;
    dNumSigmaTPC[1] = -11.;
    dNumSigmaTPC[2] = -11.;
    dNumSigmaTPC[3] = -11.;
    dNumSigmaTPC[4] = -11.;

    fhQAPIDTPCdEdx[iQAindex]->Fill(track->P(), -5.);
  }

  if(bIsTOFok)
  {
    dNumSigmaTOF[0] = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
    dNumSigmaTOF[1] = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kMuon);
    dNumSigmaTOF[2] = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);
    dNumSigmaTOF[3] = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);
    dNumSigmaTOF[4] = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);

    Double_t dTOF[5];
    track->GetIntegratedTimes(dTOF);
    dTOFbeta = dTOF[0] / track->GetTOFsignal();
    fhQAPIDTOFbeta[iQAindex]->Fill(dP,dTOFbeta);
  }
  else // TOF status not OK
  {
    dNumSigmaTOF[0] = -11.;
    dNumSigmaTOF[1] = -11.;
    dNumSigmaTOF[2] = -11.;
    dNumSigmaTOF[3] = -11.;
    dNumSigmaTOF[4] = -11.;

    fhQAPIDTOFbeta[iQAindex]->Fill(track->P(),-0.05);
  }

  fh3QAPIDnSigmaTPCTOFPtPion[iQAindex]->Fill(dNumSigmaTPC[2],dNumSigmaTOF[2],track->Pt());
  fh3QAPIDnSigmaTPCTOFPtKaon[iQAindex]->Fill(dNumSigmaTPC[3],dNumSigmaTOF[3],track->Pt());
  fh3QAPIDnSigmaTPCTOFPtProton[iQAindex]->Fill(dNumSigmaTPC[4],dNumSigmaTOF[4],track->Pt());

  fh3QAPIDnSigmaBayesElectron[iQAindex]->Fill(dNumSigmaTPC[0],dNumSigmaTOF[0],dBayesProb[0]);
  fh3QAPIDnSigmaBayesMuon[iQAindex]->Fill(dNumSigmaTPC[1],dNumSigmaTOF[1],dBayesProb[1]);
  fh3QAPIDnSigmaBayesPion[iQAindex]->Fill(dNumSigmaTPC[2],dNumSigmaTOF[2],dBayesProb[2]);
  fh3QAPIDnSigmaBayesKaon[iQAindex]->Fill(dNumSigmaTPC[3],dNumSigmaTOF[3],dBayesProb[3]);
  fh3QAPIDnSigmaBayesProton[iQAindex]->Fill(dNumSigmaTPC[4],dNumSigmaTOF[4],dBayesProb[4]);

  // species dependent QA
  switch (species)
  {
    case kPion:
      fhPIDPionPt->Fill(track->Pt());
      fhPIDPionPhi->Fill(track->Phi());
      fhPIDPionEta->Fill(track->Eta());
      fhPIDPionCharge->Fill(track->Charge());
      fh2PIDPionTPCdEdx->Fill(dPt,dTPCdEdx);
      fh2PIDPionTOFbeta->Fill(dPt,dTOFbeta);
      fh2PIDPionTPCnSigmaPion->Fill(dPt,dNumSigmaTPC[2]);
      fh2PIDPionTOFnSigmaPion->Fill(dPt,dNumSigmaTOF[2]);
      fh2PIDPionTPCnSigmaKaon->Fill(dPt,dNumSigmaTPC[3]);
      fh2PIDPionTOFnSigmaKaon->Fill(dPt,dNumSigmaTOF[3]);
      fh2PIDPionTPCnSigmaProton->Fill(dPt,dNumSigmaTPC[4]);
      fh2PIDPionTOFnSigmaProton->Fill(dPt,dNumSigmaTOF[4]);
      fh2PIDPionBayesPion->Fill(dPt,dBayesProb[2]);
      fh2PIDPionBayesKaon->Fill(dPt,dBayesProb[3]);
      fh2PIDPionBayesProton->Fill(dPt,dBayesProb[4]);
      break;

    case kKaon:
      fhPIDKaonPt->Fill(track->Pt());
      fhPIDKaonPhi->Fill(track->Phi());
      fhPIDKaonEta->Fill(track->Eta());
      fhPIDKaonCharge->Fill(track->Charge());
      fh2PIDKaonTPCdEdx->Fill(dP,dTPCdEdx);
      fh2PIDKaonTOFbeta->Fill(dP,dTOFbeta);
      fh2PIDKaonTPCnSigmaPion->Fill(dPt,dNumSigmaTPC[2]);
      fh2PIDKaonTOFnSigmaPion->Fill(dPt,dNumSigmaTOF[2]);
      fh2PIDKaonTPCnSigmaKaon->Fill(dPt,dNumSigmaTPC[3]);
      fh2PIDKaonTOFnSigmaKaon->Fill(dPt,dNumSigmaTOF[3]);
      fh2PIDKaonTPCnSigmaProton->Fill(dPt,dNumSigmaTPC[4]);
      fh2PIDKaonTOFnSigmaProton->Fill(dPt,dNumSigmaTOF[4]);
      fh2PIDKaonBayesPion->Fill(dPt,dBayesProb[2]);
      fh2PIDKaonBayesKaon->Fill(dPt,dBayesProb[3]);
      fh2PIDKaonBayesProton->Fill(dPt,dBayesProb[4]);
      break;

    case kProton:
      fhPIDProtonPt->Fill(track->Pt());
      fhPIDProtonPhi->Fill(track->Phi());
      fhPIDProtonEta->Fill(track->Eta());
      fhPIDProtonCharge->Fill(track->Charge());
      fh2PIDProtonTPCdEdx->Fill(dP,dTPCdEdx);
      fh2PIDProtonTOFbeta->Fill(dP,dTOFbeta);
      fh2PIDProtonTPCnSigmaPion->Fill(dPt,dNumSigmaTPC[2]);
      fh2PIDProtonTOFnSigmaPion->Fill(dPt,dNumSigmaTOF[2]);
      fh2PIDProtonTPCnSigmaKaon->Fill(dPt,dNumSigmaTPC[3]);
      fh2PIDProtonTOFnSigmaKaon->Fill(dPt,dNumSigmaTOF[3]);
      fh2PIDProtonTPCnSigmaProton->Fill(dPt,dNumSigmaTPC[4]);
      fh2PIDProtonTOFnSigmaProton->Fill(dPt,dNumSigmaTOF[4]);
      fh2PIDProtonBayesPion->Fill(dPt,dBayesProb[2]);
      fh2PIDProtonBayesKaon->Fill(dPt,dBayesProb[3]);
      fh2PIDProtonBayesProton->Fill(dPt,dBayesProb[4]);
      break;

    default:
      break;
  }


  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::ProcessEvent()
{
  // main method for processing of selected events:
  // - Filtering of tracks / particles for flow calculations
  // - Phi,eta,pt weights for generic framework are calculated if specified
  // - Flow calculations
  // returns kTRUE if succesfull
  // *************************************************************

  // checking the run number for aplying weights & loading TList with weights
  if(fFlowUseWeights && fFlowRunByRunWeights && (fRunNumber != fEventAOD->GetRunNumber()))
  {
    fRunNumber = fEventAOD->GetRunNumber();
    if(fFlowWeightsFile)
    {
      TList* listFlowWeights = (TList*) fFlowWeightsFile->Get(Form("%d",fRunNumber));
      if(!listFlowWeights) { AliError(Form("TList with flow weights (run %d) not found.",fRunNumber)); return kFALSE; }

      if(fFlowUse3Dweights)
      {
        fh3WeightRefs = (TH3D*) listFlowWeights->FindObject("Refs3D"); if(!fh3WeightRefs) { AliError("Refs weights not found"); return kFALSE; }
        fh3WeightCharged = (TH3D*) listFlowWeights->FindObject("Charged3D"); if(!fh3WeightCharged) { AliError("Charged weights not found"); return kFALSE; }
        fh3WeightPion = (TH3D*) listFlowWeights->FindObject("Pion3D"); if(!fh3WeightPion) { AliError("Pion weights not found"); return kFALSE; }
        fh3WeightKaon = (TH3D*) listFlowWeights->FindObject("Kaon3D"); if(!fh3WeightKaon) { AliError("Kaon weights not found"); return kFALSE; }
        fh3WeightProton = (TH3D*) listFlowWeights->FindObject("Proton3D"); if(!fh3WeightProton) { AliError("Proton weights not found"); return kFALSE; }
        fh3WeightK0s = (TH3D*) listFlowWeights->FindObject("K0s3D"); if(!fh3WeightK0s) { AliError("K0s weights not found"); return kFALSE; }
        fh3WeightLambda = (TH3D*) listFlowWeights->FindObject("Lambda3D"); if(!fh3WeightLambda) { AliError("Phi weights not found"); return kFALSE; }
        fh3WeightPhi = (TH3D*) listFlowWeights->FindObject("Phi3D"); if(!fh3WeightPhi) { AliError("Phi weights not found"); return kFALSE; }
      }
      else
      {
        fh2WeightRefs = (TH2D*) listFlowWeights->FindObject("Refs"); if(!fh2WeightRefs) { AliError("Refs weights not found"); return kFALSE; }
        fh2WeightCharged = (TH2D*) listFlowWeights->FindObject("Charged"); if(!fh2WeightCharged) { AliError("Charged weights not found"); return kFALSE; }
        fh2WeightPion = (TH2D*) listFlowWeights->FindObject("Pion"); if(!fh2WeightPion) { AliError("Pion weights not found"); return kFALSE; }
        fh2WeightKaon = (TH2D*) listFlowWeights->FindObject("Kaon"); if(!fh2WeightKaon) { AliError("Kaon weights not found"); return kFALSE; }
        fh2WeightProton = (TH2D*) listFlowWeights->FindObject("Proton"); if(!fh2WeightProton) { AliError("Proton weights not found"); return kFALSE; }
        fh2WeightK0s = (TH2D*) listFlowWeights->FindObject("K0s"); if(!fh2WeightK0s) { AliError("K0s weights not found"); return kFALSE; }
        fh2WeightLambda = (TH2D*) listFlowWeights->FindObject("Lambda"); if(!fh2WeightLambda) { AliError("Phi weights not found"); return kFALSE; }
        fh2WeightPhi = (TH2D*) listFlowWeights->FindObject("Phi"); if(!fh2WeightPhi) { AliError("Phi weights not found"); return kFALSE; }
      }
    }
  }

  // Selection of relevant particles (pushing into corresponding vectors)
  Filtering();

  // checking if there is at least one charged track selected (min requirement or <<2>>)
  // if not, event is skipped: unable to compute Reference flow (and thus any differential flow)
  if(fVectorRefs->empty()) { return kFALSE; }

  // estimate centrality & assign indexes (centrality/percentile, sampling, ...)
  fIndexCentrality = GetCentralityIndex();
  if(fIndexCentrality < 0) { return kFALSE; }

  fhEventCentrality->Fill(fIndexCentrality);
  fh2EventCentralityNumRefs->Fill(fIndexCentrality,fVectorRefs->size());
  fpRefsMult->Fill(fIndexCentrality,fVectorRefs->size(),1.0);
  // at this point, centrality index (percentile) should be properly estimated, if not, skip event

  // if running in kSkipFlow mode, skip the remaining part
  if(fRunMode == kSkipFlow) { fEventCounter++; return kTRUE; }

  // >>>> flow starts here <<<<
  // >>>> Flow a la General Framework <<<<
  for(Int_t iGap(0); iGap < fNumEtaGap; iGap++)
  {
    // Reference (pT integrated) flow
    DoFlowRefs(iGap);

    // pT differential flow
    if(!fVectorCharged->empty()) { DoFlowPOIs(iGap, kCharged); }

    if(fProcessPID)
    {
      if(!fVectorPion->empty()) { DoFlowPOIs(iGap,kPion); }
      if(!fVectorKaon->empty()) { DoFlowPOIs(iGap,kKaon); }
      if(!fVectorProton->empty()) { DoFlowPOIs(iGap,kProton); }
    }

    if(fProcessPhi)
    {
      if(!fVectorPhi->empty()) { DoFlowPOIs(iGap,kPhi); }
    }

    if(fProcessV0s)
    {
      if(!fVectorK0s->empty()) { DoFlowPOIs(iGap,kK0s); }
      if(!fVectorLambda->empty()) { DoFlowPOIs(iGap,kLambda); }
    }
  } // endfor {iGap} eta gaps

  fEventCounter++; // counter of processed events

  return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::DoFlowRefs(const Int_t iEtaGapIndex)
{
  // Estimate <2> for reference flow for all harmonics based on relevant flow vectors
  // *************************************************************

  Double_t dEtaGap = fEtaGap[iEtaGapIndex];

  FillRefsVectors(iEtaGapIndex); // filling RFPs (Q) flow vectors

  if(dEtaGap == -1.0) // no gap
  {
    // estimating <2>
    Double_t D02 = Two(0,0).Re();
    if(D02 != 0.0)
    {
      for(Int_t iHarm(0); iHarm < fNumHarmonics; ++iHarm)
      {
        Int_t iHarmonics = fHarmonics[iHarm];
        Double_t Nn2 = Two(iHarmonics,-iHarmonics).Re();
        Double_t dValue = Nn2 / D02;
        if( TMath::Abs(dValue) <= 1.0 ) { fpRefsCor2[fIndexSampling][iEtaGapIndex][iHarm]->Fill(fIndexCentrality, dValue, D02); }
      }
    }

    // estimating <4>
    if(fCutFlowDoFourCorrelations)
    {
      Double_t D04 = Four(0,0,0,0).Re();
      if(D04 != 0.0)
      {
        for(Int_t iHarm(0); iHarm < fNumHarmonics; ++iHarm)
        {
          Int_t iHarmonics = fHarmonics[iHarm];
          Double_t Nn4 = Four(iHarmonics,iHarmonics,-iHarmonics,-iHarmonics).Re();
          Double_t dValue = Nn4 / D04;
          if( TMath::Abs(dValue) <= 1.0 ) { fpRefsCor4[fIndexSampling][iHarm]->Fill(fIndexCentrality, dValue, D04); }
        }
      }
    }

  }
  else // with gap
  {
    // estimating <2>
    Double_t D02gap = TwoGap(0,0).Re();
    if(D02gap != 0.0)
    {
      for(Int_t iHarm(0); iHarm < fNumHarmonics; ++iHarm)
      {
        Int_t iHarmonics = fHarmonics[iHarm];
        Double_t Nn2gap = TwoGap(iHarmonics,-iHarmonics).Re();
        Double_t dValue = Nn2gap / D02gap;
        if( TMath::Abs(dValue) <= 1.0 ) { fpRefsCor2[fIndexSampling][iEtaGapIndex][iHarm]->Fill(fIndexCentrality, dValue, D02gap); }
      }
    }
  } // endif {dEtaGap}
  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::DoFlowPOIs(const Int_t iEtaGapIndex, const PartSpecies species)
{
  // Estimate <2> for pT diff flow of POIs particles for all harmonics based on relevant flow vectors
  // *************************************************************
  TProfile2D** prof2TwoPos = 0x0;
  TProfile2D** prof2TwoNeg = 0x0;
  TProfile2D** prof2Four = 0x0;

  TProfile3D** prof3TwoPos = 0x0;
  TProfile3D** prof3TwoNeg = 0x0;
  TProfile3D** prof3Four = 0x0;

  TAxis* axisPt = 0x0;
  TAxis* axisMass = 0x0;

  Bool_t bHasMass = kFALSE;

  // swich based on species
  switch (species)
  {
    case kCharged:
      prof2TwoPos = fp2ChargedCor2Pos[fIndexSampling][iEtaGapIndex];
      prof2TwoNeg = fp2ChargedCor2Neg[fIndexSampling][iEtaGapIndex];
      prof2Four = fp2ChargedCor4[fIndexSampling];
      break;

    case kPion:
      prof2TwoPos = fp2PionCor2Pos[fIndexSampling][iEtaGapIndex];
      prof2TwoNeg = fp2PionCor2Neg[fIndexSampling][iEtaGapIndex];
      prof2Four = fp2PionCor4[fIndexSampling];
      break;

    case kKaon:
      prof2TwoPos = fp2KaonCor2Pos[fIndexSampling][iEtaGapIndex];
      prof2TwoNeg = fp2KaonCor2Neg[fIndexSampling][iEtaGapIndex];
      prof2Four = fp2KaonCor4[fIndexSampling];
      break;

    case kProton:
      prof2TwoPos = fp2ProtonCor2Pos[fIndexSampling][iEtaGapIndex];
      prof2TwoNeg = fp2ProtonCor2Neg[fIndexSampling][iEtaGapIndex];
      prof2Four = fp2ProtonCor4[fIndexSampling];
      break;

    case kK0s:
      prof3TwoPos = fp3V0sCorrK0sCor2Pos[iEtaGapIndex];
      prof3TwoNeg = fp3V0sCorrK0sCor2Neg[iEtaGapIndex];
      prof3Four = fp3V0sCorrK0sCor4;
      bHasMass = kTRUE;
      break;

    case kLambda:
      prof3TwoPos = fp3V0sCorrLambdaCor2Pos[iEtaGapIndex];
      prof3TwoNeg = fp3V0sCorrLambdaCor2Neg[iEtaGapIndex];
      prof3Four = fp3V0sCorrLambdaCor4;
      bHasMass = kTRUE;
      break;

    case kPhi:
      prof3TwoPos = fp3PhiCorrCor2Pos[iEtaGapIndex];
      prof3TwoNeg = fp3PhiCorrCor2Neg[iEtaGapIndex];
      prof3Four = fp3PhiCorrCor4;
      bHasMass = kTRUE;
      break;

    default:
      AliError("Selected species unknown or not supported.");
      return;
  }

  Double_t dEtaGap = fEtaGap[iEtaGapIndex];

  if(bHasMass)
  {
    // particle species with inv. mass (K0s, Lambda, Phi)
    if(!prof3TwoPos) { AliError("'prof3TwoPos' not found!"); return; }
    if(dEtaGap > -1.0 && !prof3TwoNeg) { AliError("'prof3TwoNeg' not found!"); return; }
    if(fCutFlowDoFourCorrelations && !prof3Four) { AliError("'prof3Four' not found!"); return; }
    axisPt = prof3TwoPos[0]->GetYaxis(); if(!axisPt) { AliError("Pt axis object not found!"); return; }
    axisMass = prof3TwoPos[0]->GetZaxis(); if(!axisMass) { AliError("Mass axis object not found!"); return; }
  }
  else
  {
    // particle species without inv. mass (charged, pi, K, p)
    if(!prof2TwoPos) { AliError("'prof2TwoPos' not found!"); return; }
    if(dEtaGap > -1.0 && !prof2TwoNeg) { AliError("'prof2TwoNeg' not found!"); return; }
    if(fCutFlowDoFourCorrelations && !prof2Four) { AliError("'prof2Four' not found!"); return; }
    axisPt = prof2TwoPos[0]->GetYaxis(); if(!axisPt) { AliError("Pt axis object not found!"); return; }
  }

  // iNumMassBins = 1 to run at least once if bHasMass == kFALSE
  Int_t iNumMassBins = 1; if(bHasMass) { iNumMassBins = axisMass->GetNbins(); }
  Int_t iNumPtBins = axisPt->GetNbins();

  for(Int_t iMass(1); iMass < iNumMassBins+1; ++iMass)
  {
    Double_t dMass = 0.0; Double_t dMassLow = 0.0; Double_t dMassHigh = 0.0;
    if(bHasMass)
    {
      dMass = axisMass->GetBinCenter(iMass);
      dMassLow = axisMass->GetBinLowEdge(iMass);
      dMassHigh = axisMass->GetBinUpEdge(iMass);
    }

    for(Int_t iPt(1); iPt < iNumPtBins+1; ++iPt)
    {
      Double_t dPt = axisPt->GetBinCenter(iPt);
      Double_t dPtLow = axisPt->GetBinLowEdge(iPt);
      Double_t dPtHigh = axisPt->GetBinUpEdge(iPt);

      // filling POIs (P,S) flow vectors
      FillPOIsVectors(iEtaGapIndex,species,dPtLow,dPtHigh,dMassLow,dMassHigh);

      if(dEtaGap == -1.0) // no eta gap
      {
        // estimating <2'>
        Double_t D02 = TwoDiff(0,0).Re();
        if(D02 != 0.0)
        {
          for(Int_t iHarm(0); iHarm < fNumHarmonics; ++iHarm)
          {
            Int_t iHarmonics = fHarmonics[iHarm];
            Double_t Nn2 = TwoDiff(iHarmonics,-iHarmonics).Re();
            Double_t dValue = Nn2 / D02;
            if(TMath::Abs(dValue) <= 1.0)
            {
              if(bHasMass) { prof3TwoPos[iHarm]->Fill(fIndexCentrality, dPt, dMass, dValue, D02); }
              else { prof2TwoPos[iHarm]->Fill(fIndexCentrality, dPt, dValue, D02); }
            }
          }
        }

        // estimating <4'>
        if(fCutFlowDoFourCorrelations)
        {
          Double_t D04 = FourDiff(0,0,0,0).Re();
          if(D04 != 0.0)
          {
            for(Int_t iHarm(0); iHarm < fNumHarmonics; ++iHarm)
            {
              Int_t iHarmonics = fHarmonics[iHarm];
              Double_t Nn4 = FourDiff(iHarmonics,iHarmonics,-iHarmonics,-iHarmonics).Re();
              Double_t dValue = Nn4 / D04;
              if( TMath::Abs(dValue) <= 1.0 )
              {
                if(bHasMass) { prof3Four[iHarm]->Fill(fIndexCentrality, dPt, dMass, dValue, D04); }
                else { prof2Four[iHarm]->Fill(fIndexCentrality, dPt, dValue, D04); }
              }
            }
          }
        }

      }
      else // eta gap
      {
        // estimating <2'> : POIs in positive eta
        Double_t D02pos = TwoDiffGapPos(0,0).Re();
        if(D02pos != 0.0)
        {
          for(Int_t iHarm(0); iHarm < fNumHarmonics; ++iHarm)
          {
            Int_t iHarmonics = fHarmonics[iHarm];
            Double_t Nn2pos = TwoDiffGapPos(iHarmonics,-iHarmonics).Re();
            Double_t dValue = Nn2pos / D02pos;
            if( TMath::Abs(dValue) <= 1.0 )
            {
              if(bHasMass) { prof3TwoPos[iHarm]->Fill(fIndexCentrality, dPt, dMass, dValue, D02pos); }
              else { prof2TwoPos[iHarm]->Fill(fIndexCentrality, dPt, dValue, D02pos); }
            }
          }
        }

        // estimating <2'> : POIs in negative eta
        Double_t D02neg = TwoDiffGapNeg(0,0).Re();
        if(D02neg != 0.0)
        {
          for(Int_t iHarm(0); iHarm < fNumHarmonics; ++iHarm)
          {
            Int_t iHarmonics = fHarmonics[iHarm];
            Double_t Nn2neg = TwoDiffGapNeg(iHarmonics,-iHarmonics).Re();
            Double_t dValue = Nn2neg / D02neg;
            if( TMath::Abs(dValue) <= 1.0 )
            {
              if(bHasMass) { prof3TwoNeg[iHarm]->Fill(fIndexCentrality, dPt, dMass, dValue, D02neg); }
              else { prof2TwoNeg[iHarm]->Fill(fIndexCentrality, dPt, dValue, D02neg); }
            }
          }
        }
      } // endif {dEtaGap}

    } // endfor {iPt}
  } // endfor {iMass}

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FillRefsVectors(const Short_t iEtaGapIndex)
{
  // Filling Q flow vector with RFPs
  // return kTRUE if succesfull (i.e. no error occurs), kFALSE otherwise
  // *************************************************************
  Double_t dEtaGap = fEtaGap[iEtaGapIndex];
  Double_t dEtaLimit = dEtaGap / 2.0;
  Bool_t bHasGap = kFALSE;
  if(dEtaGap > -1.0) bHasGap = kTRUE;

  TH2D* h2Weights = fh2WeightRefs;
  TH3D* h3Weights = fh3WeightRefs;
  if(fFlowUseWeights)
  {
    if(fFlowUse3Dweights && !fh3WeightRefs) { AliError("Histogram with Refs weights not found."); return; }
    if(!fFlowUse3Dweights && !fh2WeightRefs) { AliError("Histogram with Refs weights not found."); return; }
  }

  // clearing output (global) flow vectors
  ResetFlowVector(fFlowVecQpos);
  ResetFlowVector(fFlowVecQneg);

  for (auto part = fVectorRefs->begin(); part != fVectorRefs->end(); part++)
  {
    Double_t dPhi = (*part)->Phi();
    Double_t dEta = (*part)->Eta();

    Double_t dWeight = 1.0;

    // loading weights if needed
    if(fFlowUseWeights)
    {
      if(fFlowUse3Dweights) { dWeight = h3Weights->GetBinContent(h3Weights->FindFixBin(dEta,dPhi,fPVz)); }
      else { dWeight = h2Weights->GetBinContent(h2Weights->FindFixBin(dEta,dPhi)); }
      if(dWeight <= 0.0) dWeight = 1.0;
    }

    if(!bHasGap) // no eta gap
    {
      for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
        for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
        {
          Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
          Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
          fFlowVecQpos[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
        }
    }
    else
    {
      if(dEta > dEtaLimit)   // RFP in positive eta acceptance
      {
        for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
          for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
          {
            Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
            Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
            fFlowVecQpos[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
          }
      }
      if(dEta < -dEtaLimit)   // RFP in negative eta acceptance
      {
        for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
          for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
          {
            Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
            Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
            fFlowVecQneg[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
          }
      }
    } // endif {dEtaGap}
  } // endfor {tracks} particle loop

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FillPOIsVectors(const Short_t iEtaGapIndex, const PartSpecies species, const Double_t dPtLow, const Double_t dPtHigh, const Double_t dMassLow, const Double_t dMassHigh)
{
  // Filling p,q and s flow vectors with POIs (given by species) for differential flow calculation
  // *************************************************************
  Bool_t bHasGap = kFALSE;
  Double_t dEtaGap = fEtaGap[iEtaGapIndex];
  if(dEtaGap > -1.0) { bHasGap = kTRUE; }
  Double_t dEtaLimit = dEtaGap / 2.0;
  Bool_t bHasMass = kFALSE;

  std::vector<AliVTrack*>* vector = 0x0;
  TH3D* histPos = 0x0;
  TH3D* histNeg = 0x0;
  TH2D* h2Weights = 0x0;
  TH3D* h3Weights = 0x0;

  // swich based on species
  switch (species)
  {
    case kCharged:
      vector = fVectorCharged;
      h2Weights = fh2WeightCharged;
      h3Weights = fh3WeightCharged;
      break;

    case kPion:
      vector = fVectorPion;
      h2Weights = fh2WeightPion;
      h3Weights = fh3WeightPion;
      break;

    case kKaon:
      vector = fVectorKaon;
      h2Weights = fh2WeightKaon;
      h3Weights = fh3WeightKaon;
      break;

    case kProton:
      vector = fVectorProton;
      h2Weights = fh2WeightProton;
      h3Weights = fh3WeightProton;
      break;

    case kK0s:
      vector = fVectorK0s;
      h2Weights = fh2WeightK0s;
      h3Weights = fh3WeightK0s;
      histPos = fh3V0sEntriesK0sPos[iEtaGapIndex];
      if(bHasGap) { histNeg = fh3V0sEntriesK0sNeg[iEtaGapIndex]; }
      bHasMass = kTRUE;
      break;

    case kLambda:
      vector = fVectorLambda;
      h2Weights = fh2WeightLambda;
      h3Weights = fh3WeightLambda;
      histPos = fh3V0sEntriesLambdaPos[iEtaGapIndex];
      if(bHasGap) { histNeg = fh3V0sEntriesLambdaNeg[iEtaGapIndex]; }
      bHasMass = kTRUE;
      break;

    case kPhi:
      vector = fVectorPhi;
      histPos = fh3PhiEntriesSignalPos[iEtaGapIndex];
      h2Weights = fh2WeightPhi;
      h3Weights = fh3WeightPhi;
      if(bHasGap) { histNeg = fh3PhiEntriesSignalNeg[iEtaGapIndex]; }
      bHasMass = kTRUE;
      break;

    default:
      AliError("Selected species unknown.");
      return;
  }

  if(!vector) { AliError("Vector with selected POIs not found."); return; }
  if(bHasMass && !histPos) { AliError("Histogram for POIs in positive eta not found."); return; }
  if(bHasMass && bHasGap && !histNeg) { AliError("Historgam for POIs in negative eta not found."); return; }
  if(fFlowUseWeights)
  {
    if(!fFlowUse3Dweights && !h2Weights) { AliError("Histogram with weights not found."); return; }
    if(fFlowUse3Dweights && !h3Weights) { AliError("Histogram with weights not found."); return; }
  }

  // clearing output (global) flow vectors
  ResetFlowVector(fFlowVecPpos);
  if(bHasGap) { ResetFlowVector(fFlowVecPneg); } else { ResetFlowVector(fFlowVecS); }

  for(auto part = vector->begin(); part != vector->end(); ++part)
  {
    Double_t dPt = (*part)->Pt();
    Double_t dPhi = (*part)->Phi();
    Double_t dEta = (*part)->Eta();
    Double_t dMass = 0.0;

    // checking if pt is within pt (bin) range
    if(dPt < dPtLow || dPt >= dPtHigh) { continue; }

    // checking if mass is within mass (bin) range
    if(bHasMass) { dMass = (*part)->M(); if(dMass < dMassLow || dMass >= dMassHigh) { continue; } }

    // at this point particles corresponding to this pt (& mass) bin survive

    Double_t dWeight = 1.0;
    if(fFlowUseWeights)
    {
      if(fFlowUse3Dweights) { dWeight = h3Weights->GetBinContent(h3Weights->FindFixBin(dEta,dPhi,fPVz)); }
      else { dWeight = h2Weights->GetBinContent(h2Weights->FindFixBin(dEta,dPhi)); }
      if(dWeight <= 0.0) dWeight = 1.0;
    }

    if(!bHasGap) // no eta gap
    {
      if(bHasMass) { histPos->Fill(fIndexCentrality,dPt,dMass,1.0); }

      for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
        for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
        {
          Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
          Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
          fFlowVecPpos[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);

          // check if track (passing criteria) is overlapping with RFPs pT region; if so, fill S (q) vector
          // in case of charged, pions, kaons or protons (one witout mass)
          if(!bHasMass && IsWithinRefs(static_cast<const AliAODTrack*>(*part)))
          {
            Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
            Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
            fFlowVecS[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
          }
        }

    }
    else // with eta gap
    {
      if(dEta > dEtaLimit) // particle in positive eta acceptance
      {
        if(bHasMass) { histPos->Fill(fIndexCentrality,dPt,dMass,1.0); }

        for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
          for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
          {
            Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
            Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
            fFlowVecPpos[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
          }
       }
       if(dEta < -dEtaLimit) // particle in negative eta acceptance
       {
         if(bHasMass) { histNeg->Fill(fIndexCentrality,dPt,dMass,1.0); }

         for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
           for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
           {
             Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
             Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
             fFlowVecPneg[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
           }
       }
     } // endif {dEtaGap}
   } // endfor {tracks}
   return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::ResetFlowVector(TComplex (&array)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax])
{
  // Reset RFPs (Q) array values to TComplex(0,0,kFALSE) for given array
  // *************************************************************
  for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
    for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
      array[iHarm][iPower] = TComplex(0,0,kFALSE);
  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::ListFlowVector(TComplex (&array)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax])
{
  // List all values of given flow vector TComplex array
  // *************************************************************
  printf(" ### Listing (TComplex) flow vector array ###########################\n");
  for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
  {
    printf("Harm %d (power):",iHarm);
    for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
    {
        printf("|(%d) %g+%g(i)",iPower, array[iHarm][iPower].Re(), array[iHarm][iPower].Im());
    }
    printf("\n");
  }
  return;
}
//_____________________________________________________________________________
Short_t AliAnalysisTaskUniFlow::GetSamplingIndex()
{
  // Assessing sampling index based on generated random number
  // returns centrality index
  // *************************************************************

  Short_t index = 0x0;

  if(fSampling && fNumSamples > 1)
  {
    TRandom3 rr(0);
    Double_t ranNum = rr.Rndm(); // getting random number in (0,1)
    Double_t generated = ranNum * fNumSamples; // getting random number in range (0, fNumSamples)

    // finding right index for sampling based on generated number and total number of samples
    for(Short_t i(0); i < fNumSamples; i++)
    {
      if(generated < (i+1) )
      {
        index = i;
        break;
      }
    }
  }

  return index;
}
//_____________________________________________________________________________
Short_t AliAnalysisTaskUniFlow::GetCentralityIndex()
{
  // Estimating centrality percentile based on selected estimator.
  // (Default) If no multiplicity estimator is specified (fMultEstimator == '' || Charged), percentile is estimated as number of selected / filtered charged tracks.
  // If a valid multiplicity estimator is specified, centrality percentile is estimated via AliMultSelection
  // otherwise -1 is returned (and event is skipped)
  // *************************************************************

  Short_t iCentralityIndex = -1;

  // assigning centrality based on number of selected charged tracks
  if( fMultEstimator.EqualTo("") || fMultEstimator.EqualTo("CHARGED") )
  {
    iCentralityIndex = fVectorRefs->size();
  }
  else if(
    // some of supported AliMultSelection estimators
    fMultEstimator.EqualTo("V0A") || fMultEstimator.EqualTo("V0C") ||
    fMultEstimator.EqualTo("V0M") || fMultEstimator.EqualTo("CL0") ||
    fMultEstimator.EqualTo("CL1") || fMultEstimator.EqualTo("ZNA") ||
    fMultEstimator.EqualTo("ZNC")
  )
  {
    // checking AliMultSelection
    AliMultSelection* multSelection = (AliMultSelection*) fEventAOD->FindListObject("MultSelection");
    if(!multSelection) { AliError("AliMultSelection object not found! Returning -1"); return -1; }

    Float_t dPercentile = multSelection->GetMultiplicityPercentile(fMultEstimator.Data());
    if(dPercentile > 100 || dPercentile < 0)
    { AliWarning("Centrality percentile estimated not within 0-100 range. Returning -1"); return -1; }

    iCentralityIndex = (Short_t) dPercentile;
  }
  else
  {
    AliWarning(Form("Multiplicity estimator '%s' not supported. Returning -1\n",fMultEstimator.Data()));
    return -1;
  }

  // transfering centrality percentile to multiplicity bin (if fixed size bins are used)
  if(fUseFixedMultBins)
  {
    for(Int_t multIndex(0); multIndex < fNumMultBins; multIndex++)
    {
      if(iCentralityIndex >= fMultBins[multIndex] && iCentralityIndex < fMultBins[multIndex+1]) { iCentralityIndex = multIndex; break; }
    }
  }

  return iCentralityIndex;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::HasTrackPIDTPC(const AliAODTrack* track)
{
  // Checks if the track has ok PID information from TPC
  // *************************************************************
  if(!track || !fPIDResponse) return kFALSE;
  AliPIDResponse::EDetPidStatus pidStatusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, track);
  return (pidStatusTPC == AliPIDResponse::kDetPidOk);
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::HasTrackPIDTOF(const AliAODTrack* track)
{
  // Checks if the track has ok PID information from TOF
  // *************************************************************
  if(!track || !fPIDResponse) return kFALSE;
  AliPIDResponse::EDetPidStatus pidStatusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, track);
  return ((pidStatusTOF == AliPIDResponse::kDetPidOk) && (track->GetStatus()& AliVTrack::kTOFout) && (track->GetStatus()& AliVTrack::kTIME));
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::Terminate(Option_t* option)
{
  // called on end of task, after all events are processed
  // *************************************************************
  AliAnalysisTaskSE::Terminate(option);

  return;
}
//_____________________________________________________________________________
// Set of methods returning given complex flow vector based on flow harmonics (n) and weight power indexes (p)
// a la General Framework implementation.
// Q: flow vector of RFPs (with/out eta gap)
// P: flow vector of POIs (with/out eta gap) (in usual notation p)
// S: flow vector of overlaping RFPs and POIs (in usual notation q)

TComplex AliAnalysisTaskUniFlow::Q(const Short_t n, const Short_t p)
{
  if (n < 0) return TComplex::Conjugate(fFlowVecQpos[-n][p]);
  else return fFlowVecQpos[n][p];
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::QGapPos(const Short_t n, const Short_t p)
{
  if (n < 0) return TComplex::Conjugate(fFlowVecQpos[-n][p]);
  else return fFlowVecQpos[n][p];
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::QGapNeg(const Short_t n, const Short_t p)
{
  if(n < 0) return TComplex::Conjugate(fFlowVecQneg[-n][p]);
  else return fFlowVecQneg[n][p];
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::P(const Short_t n, const Short_t p)
{
  if(n < 0) return TComplex::Conjugate(fFlowVecPpos[-n][p]);
  else return fFlowVecPpos[n][p];
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::PGapPos(const Short_t n, const Short_t p)
{
  if(n < 0) return TComplex::Conjugate(fFlowVecPpos[-n][p]);
  else return fFlowVecPpos[n][p];
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::PGapNeg(const Short_t n, const Short_t p)
{
  if(n < 0) return TComplex::Conjugate(fFlowVecPneg[-n][p]);
  else return fFlowVecPneg[n][p];
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::S(const Short_t n, const Short_t p)
{
  if(n < 0) return TComplex::Conjugate(fFlowVecS[-n][p]);
  else return fFlowVecS[n][p];
}
//____________________________________________________________________

// Set of flow calculation methods for cumulants of different orders with/out eta gap

TComplex AliAnalysisTaskUniFlow::Two(const Short_t n1, const Short_t n2)
{
  TComplex formula = Q(n1,1)*Q(n2,1) - Q(n1+n2,2);
  return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::TwoGap(const Short_t n1, const Short_t n2)
{
  TComplex formula = QGapPos(n1,1)*QGapNeg(n2,1);
  return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::TwoDiff(const Short_t n1, const Short_t n2)
{
  TComplex formula = P(n1,1)*Q(n2,1) - S(n1+n2,2);
  return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::TwoDiffGapPos(const Short_t n1, const Short_t n2)
{
  TComplex formula = PGapPos(n1,1)*QGapNeg(n2,1);
  return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::TwoDiffGapNeg(const Short_t n1, const Short_t n2)
{
  TComplex formula = PGapNeg(n1,1)*QGapPos(n2,1);
  return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::Four(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4)
{
  TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*Q(n1+n3,2)*Q(n4,1)
                    - Q(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.0*Q(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)
                    + Q(n2+n3,2)*Q(n1+n4,2)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)+Q(n1+n3,2)*Q(n2+n4,2)
                    + 2.0*Q(n3,1)*Q(n1+n2+n4,3)-Q(n1,1)*Q(n2,1)*Q(n3+n4,2)+Q(n1+n2,2)*Q(n3+n4,2)
                    + 2.0*Q(n2,1)*Q(n1+n3+n4,3)+2.0*Q(n1,1)*Q(n2+n3+n4,3)-6.0*Q(n1+n2+n3+n4,4);
  return formula;
}
// //____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::FourDiff(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4)
{
  TComplex formula = P(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-S(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*S(n1+n3,2)*Q(n4,1)
                    - P(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.0*S(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*S(n1+n4,2)
                    + Q(n2+n3,2)*S(n1+n4,2)-P(n1,1)*Q(n3,1)*Q(n2+n4,2)+S(n1+n3,2)*Q(n2+n4,2)
                    + 2.0*Q(n3,1)*S(n1+n2+n4,3)-P(n1,1)*Q(n2,1)*Q(n3+n4,2)+S(n1+n2,2)*Q(n3+n4,2)
                    + 2.0*Q(n2,1)*S(n1+n3+n4,3)+2.0*P(n1,1)*Q(n2+n3+n4,3)-6.0*S(n1+n2+n3+n4,4);
    return formula;
}
//____________________________________________________________________
// TComplex* AliAnalysisTaskUniFlow::FourGap(int n1, int n2, int n3, int n4)
// {
//
//   TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*Q(n1+n3,2)*Q(n4,1)
//                     - Q(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.0*Q(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)
//                     + Q(n2+n3,2)*Q(n1+n4,2)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)+Q(n1+n3,2)*Q(n2+n4,2)
//                     + 2.0*Q(n3,1)*Q(n1+n2+n4,3)-Q(n1,1)*Q(n2,1)*Q(n3+n4,2)+Q(n1+n2,2)*Q(n3+n4,2)
//                     + 2.0*Q(n2,1)*Q(n1+n3+n4,3)+2.0*Q(n1,1)*Q(n2+n3+n4,3)-6.0*Q(n1+n2+n3+n4,4);
//   TComplex *out = (TComplex*) &formula;
//   return out;
//
// }
//____________________________________________________________________
