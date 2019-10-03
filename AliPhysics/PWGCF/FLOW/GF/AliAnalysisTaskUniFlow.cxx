/**************************************************************************
* Copyright(c) 2016, ALICE Experiment at CERN, All rights reserved. *
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
// Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2016-2019
// Modifications: Zuzana Moravcova (zuzana.moravcova@cern.ch), NBI, 2019-
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
// =================================================================================================
// Overview of analysis flow (see implementation of corresonding method for details)
// 1) Event selection : IsEventSelected()
//        - based on AliEventCuts
//
// 2) Particle selection & reconstruction : Filtering()
//        - whether or not are particles processed is driven by 'fProcess?' flag setup by AliAnalysisTaskUniFlow::SetProcess?(kTRUE)
//          (except for incl. charged, which are processed anyway)
//        - filling QA plots
//        - filling GF weights
//
//    !!! here the event loop ends in case of running in 'kSkipFlow' mode
//
// 3) Flow / correlation calculations : CalculateFlow()
//        - Desired correlations are setup via general AddCorr() method, which acceps arbitrary correlation order, combination of harmonics, and value of eta gaps according to Generic Framework notation
//        - These are stored in std::vector and processed in ProcessCorrTask()
//        - As a "pre-defined" options, see AddTwo(), AddTwoGap(), AddFour(), ..., etc.
//        - Also note, that the desired output profiles are generated "automatically" based on registered CorrTasks
//
// =================================================================================================

#ifndef ALIANALYSISTASKUNIFLOW_CXX
#define ALIANALYSISTASKUNIFLOW_CXX

#include <algorithm>
#include <vector>
#include <TDatabasePDG.h>
#include <TPDGCode.h>

#include "TObjectTable.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TList.h"
#include "TComplex.h"
#include "TRandom3.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliMultSelection.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliLog.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliPicoTrack.h"
#include "AliAODv0.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"

#include "AliAnalysisTaskUniFlow.h"
#include "AliUniFlowCorrTask.h"

ClassImp(AliAnalysisTaskUniFlow);

namespace {
    const Int_t fPDGCode[] = {0,0,211,321,2212,310,3122,333};
    const Double_t fPDGMass[] = {0,0,0.13957,0.493677,0.938272,0.497614,1.11568,1.019455};
}

AliAnalysisTaskUniFlow::AliAnalysisTaskUniFlow() : AliAnalysisTaskSE(),
  fEventCuts{},
  fEventAOD{nullptr},
  fEventMC{nullptr},
  fPVz{},
  fPIDResponse{nullptr},
  fPIDCombined{nullptr},
  fFlowWeightsList{nullptr},
  fMC{kFALSE},
  fInit{kFALSE},
  fIndexSampling{0},
  fIndexCentrality{-1},
  fEventCounter{0},
  fNumEventsAnalyse{50},
  fRunNumber{-1},
  fFlowVecQpos{},
  fFlowVecQneg{},
  fFlowVecQmid{},
  fFlowVecPpos{},
  fFlowVecPneg{},
  fFlowVecSpos{},
  fFlowVecSneg{},
  fVecCorrTask{},
  fVector{},
  fRunMode{kFull},
  fAnalType{kAOD},
  fDumpTObjectTable{kFALSE},
  fSampling{kFALSE},
  fFillQA{kTRUE},
  fProcessSpec{},
  fFlowRFPsPtMin{0.2},
  fFlowRFPsPtMax{5.0},
  fFlowPOIsPtMin{0.0},
  fFlowPOIsPtMax{10.0},
  fFlowPOIsPtBinNum{0},
  fFlowPOIsPtBinEdges{},
  fFlowEtaMax{0.8},
  fFlowEtaBinNum{0},
  fFlowPhiBinNum{100},
  fPhiNumBinsMass{60},
  fV0sNumBinsMass{60},
  fNumSamples{1},
  fEtaCheckRFP{kFALSE},
  fFlowFillWeights{kTRUE},
  fFlowFillAfterWeights{kTRUE},
  fFlowUseWeights{kFALSE},
  fFlowUse3Dweights{kFALSE},
  fFlowRunByRunWeights{kTRUE},
  fFlowWeightsApplyForReco{kTRUE},
  fFlowWeightsTag{},
  fColSystem{kPPb},
  fTrigger{AliVEvent::kINT7},
  fCentEstimator{kV0A},
  fCentMin{0},
  fCentMax{0},
  fCentBinNum{0},
  fCentEstimatorAdd{kRFP},
  fCentMinAdd{0},
  fCentMaxAdd{0},
  fPVtxCutZ{10.0},
  fEventRejectAddPileUp{kFALSE},
  fCutChargedTrackFilterBit{96},
  fCutChargedNumTPCclsMin{70},
  fCutChargedDCAzMax{0.0},
  fCutChargedDCAxyMax{0.0},
  fCutPIDUseAntiProtonOnly{kFALSE},
  fCutPIDnSigmaCombinedTOFrejection{kTRUE},
  fCutUseBayesPID{kFALSE},
  fCutPIDnSigmaTPCRejectElectron{0.0},
  fCutPIDnSigmaMax{},
  fCutPIDBayesMin{},
  fCutV0sOnFly{kFALSE},
  fCutV0srefitTPC{kTRUE},
  fCutV0srejectKinks{kTRUE},
  fCutV0sDaughterNumTPCClsMin{0},
  fCutV0sDaughterNumTPCCrossMin{70},
  fCutV0sDaughterNumTPCFindMin{1},
  fCutV0sDaughterNumTPCClsPIDMin{0},
  fCutV0sDaughterRatioCrossFindMin{0.8},
  fCutV0sCrossMassRejection{kTRUE},
  fCutV0sCrossMassCutK0s{0.005},
  fCutV0sCrossMassCutLambda{0.010},
  fCutV0sDCAtoPVMin{0.06},
  fCutV0sDCAtoPVMax{0.0},
  fCutV0sDCAtoPVzMax{0.0},
  fCutV0sDCADaughtersMin{0.},
  fCutV0sDCADaughtersMax{1.0},
  fCutV0sDecayRadiusMin{0.5},
  fCutV0sDecayRadiusMax{200.0},
  fCutV0sDaughterFilterBit{0},
  fCutV0sDaughterPtMin{0.0},
  fCutV0sDaughterPtMax{0.0},
  fCutV0sDaughterEtaMax{0.8},
  fCutV0sMotherRapMax{0.0},
  fCutV0sCPAK0sMin{0.97},
  fCutV0sCPALambdaMin{0.995},
  fCutV0sNumTauK0sMax{7.46},
  fCutV0sNumTauLambdaMax{3.8},
  fCutV0sInvMassK0sMin{0.4},
  fCutV0sInvMassK0sMax{0.6},
  fCutV0sInvMassLambdaMin{1.08},
  fCutV0sInvMassLambdaMax{1.16},
  fCutV0sArmenterosAlphaK0sMin{0.2},
  fCutV0sArmenterosAlphaLambdaMax{0.0},
  fCutV0sK0sPionNumTPCSigmaMax{5.0},
  fCutV0sLambdaPionNumTPCSigmaMax{5.0},
  fCutV0sLambdaProtonNumTPCSigmaMax{5.0},
  fCutPhiInvMassMin{0.99},
  fCutPhiInvMassMax{1.07},

  fQAEvents{nullptr},
  fQACharged{nullptr},
  fQAPID{nullptr},
  fQAV0s{nullptr},
  fQAPhi{nullptr},
  fFlowWeights{nullptr},
  fListFlow{nullptr},
  fListMC{nullptr},

  fhsCandK0s{nullptr},
  fhsCandLambda{nullptr},
  fhsCandPhi{nullptr},
  fhsCandPhiBg{nullptr},
  fh2Weights{nullptr},
  fh3Weights{nullptr},
  fh2AfterWeights{nullptr},
  fh3AfterWeights{nullptr},
  fhEventSampling{nullptr},
  fhEventCentrality{nullptr},
  fh2EventCentralityNumRefs{nullptr},
  fhEventCounter{nullptr},
  fhRefsMult{nullptr},
  fhRefsPt{nullptr},
  fhRefsEta{nullptr},
  fhRefsPhi{nullptr},
  fpRefsMult{nullptr},
  fhChargedCounter{nullptr},
  fhPIDCounter{nullptr},
  fhPIDMult{nullptr},
  fhPIDPt{nullptr},
  fhPIDPhi{nullptr},
  fhPIDEta{nullptr},
  fhPIDCharge{nullptr},
  fh2PIDTPCdEdx{nullptr},
  fh2PIDTPCdEdxDelta{nullptr},
  fh2PIDTOFbeta{nullptr},
  fh2PIDTOFbetaDelta{nullptr},
  fh2PIDBayesElectron{nullptr},
  fh2PIDBayesMuon{nullptr},
  fh2PIDBayesPion{nullptr},
  fh2PIDBayesKaon{nullptr},
  fh2PIDBayesProton{nullptr},
  fh2PIDTPCnSigmaElectron{nullptr},
  fh2PIDTOFnSigmaElectron{nullptr},
  fh2PIDTPCnSigmaMuon{nullptr},
  fh2PIDTOFnSigmaMuon{nullptr},
  fh2PIDTPCnSigmaPion{nullptr},
  fh2PIDTOFnSigmaPion{nullptr},
  fh2PIDTPCnSigmaKaon{nullptr},
  fh2PIDTOFnSigmaKaon{nullptr},
  fh2PIDTPCnSigmaProton{nullptr},
  fh2PIDTOFnSigmaProton{nullptr},
  fh2MCPtEtaGen{nullptr},
  fh2MCPtEtaReco{nullptr},
  fh2MCPtEtaRecoTrue{nullptr},
  fhPhiCounter{nullptr},
  fhPhiMult{nullptr},
  fhPhiBGMult{nullptr},
  fhPhiInvMass{nullptr},
  fhPhiBGInvMass{nullptr},
  fhPhiCharge{nullptr},
  fhPhiBGCharge{nullptr},
  fhPhiPt{nullptr},
  fhPhiEta{nullptr},
  fhPhiPhi{nullptr},
  fhV0sCounter{nullptr},
  fhV0sCounterK0s{nullptr},
  fhV0sCounterLambda{nullptr},
  fhV0sInvMassK0s{nullptr},
  fhV0sInvMassLambda{nullptr},
  fhV0sCompetingInvMassK0s{nullptr},
  fhV0sCompetingInvMassLambda{nullptr},
  fhQAEventsPVz{nullptr},
  fhQAEventsNumContrPV{nullptr},
  fhQAEventsNumSPDContrPV{nullptr},
  fhQAEventsDistPVSPD{nullptr},
  fhQAEventsSPDresol{nullptr},
  fhQAEventsfMult32vsCentr{nullptr},
  fhQAEventsMult128vsCentr{nullptr},
  fhQAEventsfMultTPCvsTOF{nullptr},
  fhQAEventsfMultTPCvsESD{nullptr},
  fhQAChargedMult{nullptr},
  fhQAChargedPt{nullptr},
  fhQAChargedEta{nullptr},
  fhQAChargedPhi{nullptr},
  fhQAChargedCharge{nullptr},
  fhQAChargedNumTPCcls{nullptr},
  fhQAChargedDCAxy{nullptr},
  fhQAChargedDCAz{nullptr},
  fhQAPIDTPCstatus{nullptr},
  fhQAPIDTOFstatus{nullptr},
  fhQAPIDTPCdEdx{nullptr},
  fhQAPIDTOFbeta{nullptr},
  fh3QAPIDnSigmaTPCTOFPtPion{nullptr},
  fh3QAPIDnSigmaTPCTOFPtKaon{nullptr},
  fh3QAPIDnSigmaTPCTOFPtProton{nullptr},
  fhQAV0sMultK0s{nullptr},
  fhQAV0sMultLambda{nullptr},
  fhQAV0sMultALambda{nullptr},
  fhQAV0sRecoMethod{nullptr},
  fhQAV0sDaughterTPCRefit{nullptr},
  fhQAV0sDaughterKinks{nullptr},
  fhQAV0sDaughterNumTPCCls{nullptr},
  fhQAV0sDaughterNumTPCFind{nullptr},
  fhQAV0sDaughterNumTPCCrossRows{nullptr},
  fhQAV0sDaughterTPCCrossFindRatio{nullptr},
  fhQAV0sDaughterNumTPCClsPID{nullptr},
  fhQAV0sDCAtoPV{nullptr},
  fhQAV0sDCADaughters{nullptr},
  fhQAV0sDecayRadius{nullptr},
  fhQAV0sInvMassK0s{nullptr},
  fhQAV0sInvMassLambda{nullptr},
  fhQAV0sMotherPt{nullptr},
  fhQAV0sMotherPhi{nullptr},
  fhQAV0sMotherEta{nullptr},
  fhQAV0sMotherCharge{nullptr},
  fhQAV0sMotherRapK0s{nullptr},
  fhQAV0sMotherRapLambda{nullptr},
  fhQAV0sDaughterPt{nullptr},
  fhQAV0sDaughterPhi{nullptr},
  fhQAV0sDaughterEta{nullptr},
  fhQAV0sDaughterCharge{nullptr},
  fhQAV0sDaughterTPCstatus{nullptr},
  fhQAV0sDaughterTOFstatus{nullptr},
  fhQAV0sDaughterTPCdEdxK0s{nullptr},
  fhQAV0sDaughterNumSigmaPionK0s{nullptr},
  fhQAV0sDaughterTPCdEdxLambda{nullptr},
  fhQAV0sDaughterNumSigmaPionLambda{nullptr},
  fhQAV0sDaughterNumSigmaProtonLambda{nullptr},
  fhQAV0sDaughterNumSigmaPionALambda{nullptr},
  fhQAV0sDaughterNumSigmaProtonALambda{nullptr},
  fhQAV0sCPAK0s{nullptr},
  fhQAV0sCPALambda{nullptr},
  fhQAV0sNumTauK0s{nullptr},
  fhQAV0sNumTauLambda{nullptr},
  fhQAV0sArmenterosK0s{nullptr},
  fhQAV0sArmenterosLambda{nullptr},
  fhQAV0sArmenterosALambda{nullptr}
{
  // default constructor, don't allocate memory here!
  // this is used by root for IO purposes, it needs to remain empty
}
// ============================================================================
AliAnalysisTaskUniFlow::AliAnalysisTaskUniFlow(const char* name, ColSystem colSys, Bool_t bUseWeights, Bool_t bIsMC) : AliAnalysisTaskSE(name),
  fEventCuts{},
  fEventAOD{nullptr},
  fEventMC{nullptr},
  fPVz{},
  fPIDResponse{nullptr},
  fPIDCombined{nullptr},
  fFlowWeightsList{nullptr},
  fMC{bIsMC},
  fInit{kFALSE},
  fIndexSampling{0},
  fIndexCentrality{-1},
  fEventCounter{0},
  fNumEventsAnalyse{50},
  fRunNumber{-1},
  fFlowVecQpos{},
  fFlowVecQneg{},
  fFlowVecQmid{},
  fFlowVecPpos{},
  fFlowVecPneg{},
  fFlowVecSpos{},
  fFlowVecSneg{},
  fVecCorrTask{},
  fVector{},
  fRunMode{kFull},
  fAnalType{kAOD},
  fDumpTObjectTable{kFALSE},
  fSampling{kFALSE},
  fFillQA{kTRUE},
  fProcessSpec{},
  fFlowRFPsPtMin{0.2},
  fFlowRFPsPtMax{5.0},
  fFlowPOIsPtMin{0.0},
  fFlowPOIsPtMax{10.0},
  fFlowPOIsPtBinNum{0},
  fFlowPOIsPtBinEdges{},
  fFlowEtaMax{0.8},
  fFlowEtaBinNum{0},
  fFlowPhiBinNum{100},
  fPhiNumBinsMass{60},
  fV0sNumBinsMass{60},
  fNumSamples{1},
  fEtaCheckRFP{kFALSE},
  fFlowFillWeights{kTRUE},
  fFlowFillAfterWeights{kTRUE},
  fFlowUseWeights{bUseWeights},
  fFlowUse3Dweights{kFALSE},
  fFlowRunByRunWeights{kTRUE},
  fFlowWeightsApplyForReco{kTRUE},
  fFlowWeightsTag{},
  fColSystem{colSys},
  fTrigger{AliVEvent::kINT7},
  fCentEstimator{kV0A},
  fCentMin{0},
  fCentMax{0},
  fCentBinNum{0},
  fCentEstimatorAdd{kRFP},
  fCentMinAdd{0},
  fCentMaxAdd{0},
  fPVtxCutZ{10.0},
  fEventRejectAddPileUp{kFALSE},
  fCutChargedTrackFilterBit{96},
  fCutChargedNumTPCclsMin{70},
  fCutChargedDCAzMax{0.0},
  fCutChargedDCAxyMax{0.0},
  fCutPIDUseAntiProtonOnly{kFALSE},
  fCutPIDnSigmaCombinedTOFrejection{kTRUE},
  fCutUseBayesPID{kFALSE},
  fCutPIDnSigmaTPCRejectElectron{0.0},
  fCutPIDnSigmaMax{},
  fCutPIDBayesMin{},
  fCutV0sOnFly{kFALSE},
  fCutV0srefitTPC{kTRUE},
  fCutV0srejectKinks{kTRUE},
  fCutV0sDaughterNumTPCClsMin{0},
  fCutV0sDaughterNumTPCCrossMin{70},
  fCutV0sDaughterNumTPCFindMin{1},
  fCutV0sDaughterNumTPCClsPIDMin{0},
  fCutV0sDaughterRatioCrossFindMin{0.8},
  fCutV0sCrossMassRejection{kTRUE},
  fCutV0sCrossMassCutK0s{0.005},
  fCutV0sCrossMassCutLambda{0.010},
  fCutV0sDCAtoPVMin{0.06},
  fCutV0sDCAtoPVMax{0.0},
  fCutV0sDCAtoPVzMax{0.0},
  fCutV0sDCADaughtersMin{0.},
  fCutV0sDCADaughtersMax{1.0},
  fCutV0sDecayRadiusMin{0.5},
  fCutV0sDecayRadiusMax{200.0},
  fCutV0sDaughterFilterBit{0},
  fCutV0sDaughterPtMin{0.0},
  fCutV0sDaughterPtMax{0.0},
  fCutV0sDaughterEtaMax{0.8},
  fCutV0sMotherRapMax{0.0},
  fCutV0sCPAK0sMin{0.97},
  fCutV0sCPALambdaMin{0.995},
  fCutV0sNumTauK0sMax{7.46},
  fCutV0sNumTauLambdaMax{3.8},
  fCutV0sInvMassK0sMin{0.4},
  fCutV0sInvMassK0sMax{0.6},
  fCutV0sInvMassLambdaMin{1.08},
  fCutV0sInvMassLambdaMax{1.16},
  fCutV0sArmenterosAlphaK0sMin{0.2},
  fCutV0sArmenterosAlphaLambdaMax{0.0},
  fCutV0sK0sPionNumTPCSigmaMax{5.0},
  fCutV0sLambdaPionNumTPCSigmaMax{5.0},
  fCutV0sLambdaProtonNumTPCSigmaMax{5.0},
  fCutPhiInvMassMin{0.99},
  fCutPhiInvMassMax{1.07},

  fQAEvents{nullptr},
  fQACharged{nullptr},
  fQAPID{nullptr},
  fQAV0s{nullptr},
  fQAPhi{nullptr},
  fFlowWeights{nullptr},
  fListFlow{nullptr},
  fListMC{nullptr},

  fhsCandK0s{nullptr},
  fhsCandLambda{nullptr},
  fhsCandPhi{nullptr},
  fhsCandPhiBg{nullptr},
  fh2Weights{nullptr},
  fh3Weights{nullptr},
  fh2AfterWeights{nullptr},
  fh3AfterWeights{nullptr},
  fhEventSampling{nullptr},
  fhEventCentrality{nullptr},
  fh2EventCentralityNumRefs{nullptr},
  fhEventCounter{nullptr},
  fhRefsMult{nullptr},
  fhRefsPt{nullptr},
  fhRefsEta{nullptr},
  fhRefsPhi{nullptr},
  fpRefsMult{nullptr},
  fhChargedCounter{nullptr},
  fhPIDCounter{nullptr},
  fhPIDMult{nullptr},
  fhPIDPt{nullptr},
  fhPIDPhi{nullptr},
  fhPIDEta{nullptr},
  fhPIDCharge{nullptr},
  fh2PIDTPCdEdx{nullptr},
  fh2PIDTPCdEdxDelta{nullptr},
  fh2PIDTOFbeta{nullptr},
  fh2PIDTOFbetaDelta{nullptr},
  fh2PIDBayesElectron{nullptr},
  fh2PIDBayesMuon{nullptr},
  fh2PIDBayesPion{nullptr},
  fh2PIDBayesKaon{nullptr},
  fh2PIDBayesProton{nullptr},
  fh2PIDTPCnSigmaElectron{nullptr},
  fh2PIDTOFnSigmaElectron{nullptr},
  fh2PIDTPCnSigmaMuon{nullptr},
  fh2PIDTOFnSigmaMuon{nullptr},
  fh2PIDTPCnSigmaPion{nullptr},
  fh2PIDTOFnSigmaPion{nullptr},
  fh2PIDTPCnSigmaKaon{nullptr},
  fh2PIDTOFnSigmaKaon{nullptr},
  fh2PIDTPCnSigmaProton{nullptr},
  fh2PIDTOFnSigmaProton{nullptr},
  fh2MCPtEtaGen{nullptr},
  fh2MCPtEtaReco{nullptr},
  fh2MCPtEtaRecoTrue{nullptr},
  fhPhiCounter{nullptr},
  fhPhiMult{nullptr},
  fhPhiBGMult{nullptr},
  fhPhiInvMass{nullptr},
  fhPhiBGInvMass{nullptr},
  fhPhiCharge{nullptr},
  fhPhiBGCharge{nullptr},
  fhPhiPt{nullptr},
  fhPhiEta{nullptr},
  fhPhiPhi{nullptr},
  fhV0sCounter{nullptr},
  fhV0sCounterK0s{nullptr},
  fhV0sCounterLambda{nullptr},
  fhV0sInvMassK0s{nullptr},
  fhV0sInvMassLambda{nullptr},
  fhV0sCompetingInvMassK0s{nullptr},
  fhV0sCompetingInvMassLambda{nullptr},
  fhQAEventsPVz{nullptr},
  fhQAEventsNumContrPV{nullptr},
  fhQAEventsNumSPDContrPV{nullptr},
  fhQAEventsDistPVSPD{nullptr},
  fhQAEventsSPDresol{nullptr},
  fhQAEventsfMult32vsCentr{nullptr},
  fhQAEventsMult128vsCentr{nullptr},
  fhQAEventsfMultTPCvsTOF{nullptr},
  fhQAEventsfMultTPCvsESD{nullptr},
  fhQAChargedMult{nullptr},
  fhQAChargedPt{nullptr},
  fhQAChargedEta{nullptr},
  fhQAChargedPhi{nullptr},
  fhQAChargedCharge{nullptr},
  fhQAChargedNumTPCcls{nullptr},
  fhQAChargedDCAxy{nullptr},
  fhQAChargedDCAz{nullptr},
  fhQAPIDTPCstatus{nullptr},
  fhQAPIDTOFstatus{nullptr},
  fhQAPIDTPCdEdx{nullptr},
  fhQAPIDTOFbeta{nullptr},
  fh3QAPIDnSigmaTPCTOFPtPion{nullptr},
  fh3QAPIDnSigmaTPCTOFPtKaon{nullptr},
  fh3QAPIDnSigmaTPCTOFPtProton{nullptr},
  fhQAV0sMultK0s{nullptr},
  fhQAV0sMultLambda{nullptr},
  fhQAV0sMultALambda{nullptr},
  fhQAV0sRecoMethod{nullptr},
  fhQAV0sDaughterTPCRefit{nullptr},
  fhQAV0sDaughterKinks{nullptr},
  fhQAV0sDaughterNumTPCCls{nullptr},
  fhQAV0sDaughterNumTPCFind{nullptr},
  fhQAV0sDaughterNumTPCCrossRows{nullptr},
  fhQAV0sDaughterTPCCrossFindRatio{nullptr},
  fhQAV0sDaughterNumTPCClsPID{nullptr},
  fhQAV0sDCAtoPV{nullptr},
  fhQAV0sDCADaughters{nullptr},
  fhQAV0sDecayRadius{nullptr},
  fhQAV0sInvMassK0s{nullptr},
  fhQAV0sInvMassLambda{nullptr},
  fhQAV0sMotherPt{nullptr},
  fhQAV0sMotherPhi{nullptr},
  fhQAV0sMotherEta{nullptr},
  fhQAV0sMotherCharge{nullptr},
  fhQAV0sMotherRapK0s{nullptr},
  fhQAV0sMotherRapLambda{nullptr},
  fhQAV0sDaughterPt{nullptr},
  fhQAV0sDaughterPhi{nullptr},
  fhQAV0sDaughterEta{nullptr},
  fhQAV0sDaughterCharge{nullptr},
  fhQAV0sDaughterTPCstatus{nullptr},
  fhQAV0sDaughterTOFstatus{nullptr},
  fhQAV0sDaughterTPCdEdxK0s{nullptr},
  fhQAV0sDaughterNumSigmaPionK0s{nullptr},
  fhQAV0sDaughterTPCdEdxLambda{nullptr},
  fhQAV0sDaughterNumSigmaPionLambda{nullptr},
  fhQAV0sDaughterNumSigmaProtonLambda{nullptr},
  fhQAV0sDaughterNumSigmaPionALambda{nullptr},
  fhQAV0sDaughterNumSigmaProtonALambda{nullptr},
  fhQAV0sCPAK0s{nullptr},
  fhQAV0sCPALambda{nullptr},
  fhQAV0sNumTauK0s{nullptr},
  fhQAV0sNumTauLambda{nullptr},
  fhQAV0sArmenterosK0s{nullptr},
  fhQAV0sArmenterosLambda{nullptr},
  fhQAV0sArmenterosALambda{nullptr}
{
  // defining input/output
  DefineInput(0, TChain::Class());
  if(fFlowUseWeights) { DefineInput(1, TList::Class()); }
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
  DefineOutput(13, TList::Class());
  DefineOutput(14, TList::Class());
  if(fMC) DefineOutput(15, TList::Class());
}
// ============================================================================
AliAnalysisTaskUniFlow::~AliAnalysisTaskUniFlow()
{
  // destructor
  // if(fPIDCombined)
  // {
  //   delete fPIDCombined;
  // }


  DumpTObjTable("Destructor: start");

  // clearing vectors before deleting
  ClearVectors();

  // deleting FlowPart vectors (containers)
  for(Int_t iSpec(0); iSpec < kUnknown; ++iSpec)
  {
    if(fVector[iSpec]) delete fVector[iSpec];
    if(fListFlow[iSpec]) delete fListFlow[iSpec];
  }

  for(Size_t i(0); i < fVecCorrTask.size(); ++i) { delete fVecCorrTask.at(i); }

  // deleting output lists
  if(fFlowWeights) delete fFlowWeights;
  if(fQAEvents) delete fQAEvents;
  if(fQACharged) delete fQACharged;
  if(fQAPID) delete fQAPID;
  if(fQAPhi) delete fQAPhi;
  if(fQAV0s) delete fQAV0s;

  DumpTObjTable("Destructor: end");

}
// ============================================================================
const char* AliAnalysisTaskUniFlow::GetSpeciesName(const PartSpecies species) const
{
  const char* name;

  switch(species) {
    case kRefs: name = "Refs"; break;
    case kCharged: name = "Charged"; break;
    case kPion: name = "Pion"; break;
    case kKaon: name = "Kaon"; break;
    case kProton: name = "Proton"; break;
    case kK0s: name = "K0s"; break;
    case kLambda: name = "Lambda"; break;
    case kPhi: name = "Phi"; break;
    default: name = "Unknown";
  }

  return name;
}
// ============================================================================
const char* AliAnalysisTaskUniFlow::GetSpeciesLabel(const PartSpecies species) const
{
  const char* label;

  switch(species) {
    case kRefs: label = "RFP"; break;
    case kCharged: label = "h^{#pm}"; break;
    case kPion: label = "#pi^{#pm}"; break;
    case kKaon: label = "K^{#pm}"; break;
    case kProton: label = "p(#bar{p})"; break;
    case kK0s: label = "K^{0}_{S}"; break;
    case kLambda: label = "#Lambda(#bar{#Lambda})"; break;
    case kPhi: label = "#phi"; break;
    default: label = "NA";
  }

  return label;
}
// ============================================================================
void AliAnalysisTaskUniFlow::ListParameters() const
{
  // lists all task parameters
  // *************************************************************
  AliInfo("Listing all AliAnalysisTaskUniFlow parameters");
  printf("   -------- Analysis task ---------------------------------------\n");
  printf("      fRunMode: (RunMode) %d\n",    fRunMode);
  printf("      fAnalType: (AnalType) %d\n",    fAnalType);
  printf("      fMC: (Bool_t) %s\n",    fMC ? "kTRUE" : "kFALSE");
  printf("      fDumpTObjectTable: (Bool_t) %s\n",    fDumpTObjectTable ? "kTRUE" : "kFALSE");
  printf("      fSampling: (Bool_t) %s\n",    fSampling ? "kTRUE" : "kFALSE");
  printf("      fFillQA: (Bool_t) %s\n",    fFillQA ? "kTRUE" : "kFALSE");
  printf("      fEtaCheckRFP: (Bool_t) %s\n",    fEtaCheckRFP ? "kTRUE" : "kFALSE");
  for(Int_t iSpec(0); iSpec < kUnknown; ++iSpec) { printf("      fProcessSpec[k%s]: (Bool_t) %s\n",   GetSpeciesName(PartSpecies(iSpec)), fProcessSpec[iSpec] ? "kTRUE" : "kFALSE"); }
  printf("   -------- Flow related ----------------------------------------\n");
  printf("      fFlowRFPsPtMin: (Double_t) %g (GeV/c)\n",    fFlowRFPsPtMin);
  printf("      fFlowRFPsPtMax: (Double_t) %g (GeV/c)\n",    fFlowRFPsPtMax);
  printf("      fFlowPOIsPtMin: (Double_t) %g (GeV/c)\n",    fFlowPOIsPtMin);
  printf("      fFlowPOIsPtMax: (Double_t) %g (GeV/c)\n",    fFlowPOIsPtMax);
  printf("      fFlowPOIsPtBinNum: (Int_t) %d\n",    fFlowPOIsPtBinNum);
  printf("      fFlowEtaMax: (Double_t) %g\n",    fFlowEtaMax);
  printf("      fFlowEtaBinNum: (Int_t) %d\n",    fFlowEtaBinNum);
  printf("      fFlowPhiBinNum: (Int_t) %d\n",    fFlowPhiBinNum);
  printf("      fV0sNumBinsMass: (Int_t) %d\n",    fV0sNumBinsMass);
  printf("      fPhiNumBinsMass: (Int_t) %d\n",    fPhiNumBinsMass);
  printf("      fFlowFillWeights: (Bool_t) %s\n",    fFlowFillWeights ? "kTRUE" : "kFALSE");
  printf("      fFlowFillAfterWeights: (Bool_t) %s\n",    fFlowFillAfterWeights ? "kTRUE" : "kFALSE");
  printf("      fFlowUseWeights: (Bool_t) %s\n",    fFlowUseWeights ? "kTRUE" : "kFALSE");
  printf("      fFlowWeightsTag: (TString) '%s'\n",    fFlowWeightsTag.Data());
  printf("      fFlowRunByRunWeights: (Bool_t) %s\n",    fFlowRunByRunWeights ? "kTRUE" : "kFALSE");
  printf("      fFlowUse3Dweights: (Bool_t) %s\n",    fFlowUse3Dweights ? "kTRUE" : "kFALSE");
  printf("      fFlowWeightsApplyForReco: (Bool_t) %s\n",    fFlowWeightsApplyForReco ? "kTRUE" : "kFALSE");
  printf("   -------- Events ----------------------------------------------\n");
  printf("      fColSystem: (ColSystem) %d\n",    fColSystem);
  printf("      fTrigger: (Short_t) %d\n",    fTrigger);
  printf("      fCentEstimator: (CentEst) '%s' (%d)\n",    GetCentEstimatorLabel(fCentEstimator), fCentEstimator);
  printf("      fCentMin: (Int_t) %d\n",    fCentMin);
  printf("      fCentMax: (Int_t) %d\n",    fCentMax);
  printf("      fCentBinNum: (Int_t) %d\n",    fCentBinNum);
  printf("      fCentEstimatorAdd: (CentEst) '%s' (%d)\n",    GetCentEstimatorLabel(fCentEstimatorAdd), fCentEstimatorAdd);
  printf("      fCentMinAdd: (Int_t) %d\n",    fCentMinAdd);
  printf("      fCentMaxAdd: (Int_t) %d\n",    fCentMaxAdd);
  printf("      fPVtxCutZ: (Double_t) %g (cm)\n",    fPVtxCutZ);
  printf("   -------- Charged tracks --------------------------------------\n");
  printf("      fCutChargedTrackFilterBit: (UInt) %d\n",    fCutChargedTrackFilterBit);
  printf("      fCutChargedNumTPCclsMin: (UShort_t) %d\n",    fCutChargedNumTPCclsMin);
  printf("      fCutChargedDCAzMax: (Double_t) %g (cm)\n",    fCutChargedDCAzMax);
  printf("      fCutChargedDCAxyMax: (Double_t) %g (cm)\n",    fCutChargedDCAxyMax);
  printf("   -------- PID (pi,K,p) tracks ---------------------------------\n");
  printf("      fCutPIDUseAntiProtonOnly: (Bool_t) %s\n",  fCutPIDUseAntiProtonOnly ? "kTRUE" : "kFALSE");
  printf("      fCutPIDnSigmaCombinedTOFrejection: (Bool_t) %s\n",  fCutPIDnSigmaCombinedTOFrejection ? "kTRUE" : "kFALSE");
  for(Int_t iSpec(2); iSpec < fPIDNumSpecies; ++iSpec) { printf("      fCutPIDnSigmaPionMax: (Float_t) %g\n",    fCutPIDnSigmaMax[iSpec]); }
  printf("      fCutPIDnSigmaTPCRejectElectron: (Float_t) %g\n",    fCutPIDnSigmaTPCRejectElectron);
  printf("      fCutUseBayesPID: (Bool_t) %s\n",    fCutUseBayesPID ? "kTRUE" : "kFALSE");
  for(Int_t iSpec(2); iSpec < fPIDNumSpecies; ++iSpec) { printf("      fCutPIDBayesMin: (Double_t) %g\n",    fCutPIDBayesMin[iSpec]); }
  printf("   -------- Phi candidates --------------------------------------\n");
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
// ============================================================================
void AliAnalysisTaskUniFlow::ClearVectors()
{
  // Properly clear all particle vectors (if exists)
  // NOTE: should be called at the end of each event & before vectors deleting
  // *************************************************************

  for(Int_t iSpec(0); iSpec < kUnknown; ++iSpec) {
    if(!fProcessSpec[iSpec]) { continue; }
    std::vector<AliVTrack*>* vector = fVector[iSpec];
    if(!vector) { continue; }
    if(HasMass(PartSpecies(iSpec))) { for(auto part = vector->begin(); part != vector->end(); ++part) { delete *part; } }
    vector->clear();
  }

  return;
}
// ============================================================================
void AliAnalysisTaskUniFlow::DumpTObjTable(const char* note, Option_t* opt) const
{
  // Skipping if flag is off
  if(!fDumpTObjectTable) { return; }

  if(note) {
    printf("TObjectTable::%s",note);
  }

  gObjectTable->Print(opt);
}
// ============================================================================
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
  fPIDCombined->SetDefaultTPCPriors();
  fPIDCombined->SetSelectedSpecies(fPIDNumSpecies);
  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF); // setting TPC + TOF mask

  // checking consistency of PartSpecies with AliPID::EParticleType
  if((Int_t) kPion != (Int_t) AliPID::kPion || (Int_t) kKaon != (Int_t) AliPID::kKaon || (Int_t) kProton != (Int_t) AliPID::kProton)
  {
    AliFatal("PartSpecies enum not consistent with AliPID::EParticleType! Needed for PID selection!");
    return kFALSE;
  }

  // checking the fFlowNumHarmonicsMax, fFlowNumWeightPowersMax dimensions of p,Q,S vectors
  if(fFlowNumWeightPowersMax < 5) { AliFatal("Low range of flow vector weight dimension! Not enought for <4>!"); return kFALSE; }
  // TODO fFlowNumHarmonicsMax;

  if(fSampling && fNumSamples < 2)
  {
    AliFatal("Sampling used, but number of samples < 2! Terminating!");
    return kFALSE;
  }

  // checking cut setting
  // Refs
  AliInfo("Checking task parameters setting conflicts (ranges, etc)");
  if(fFlowRFPsPtMin > 0. && fFlowRFPsPtMax > 0. && fFlowRFPsPtMin > fFlowRFPsPtMax)
  {
    AliFatal("Cut: RFPs Pt range wrong! Terminating!");
    return kFALSE;
  }

  // POIs
  if(fFlowPOIsPtMin > 0. && fFlowPOIsPtMax > 0. && fFlowPOIsPtMin > fFlowPOIsPtMax)
  {
    AliFatal("Cut: POIs Pt range wrong! Terminating!");
    return kFALSE;
  }

  // setting POIs Pt binning
  if(fFlowPOIsPtBinNum < 1) { fFlowPOIsPtBinNum = (Int_t) ((fFlowPOIsPtMax - fFlowPOIsPtMin) / 0.1 + 0.5); }

  // Inv. Mass
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

  if(fCutPhiInvMassMin > fCutPhiInvMassMax || fCutPhiInvMassMin < 0. || fCutPhiInvMassMax < 0.)
  {
    AliFatal("Cut: InvMass (Phi) range wrong! Terminating!");
    return kFALSE;
  }

  // eta
  if(fFlowEtaMax < 0.)
  {
    AliFatal("Cut: Eta max. wrong! Terminating!");
    return kFALSE;
  }
  if(fFlowEtaBinNum < 1)
  {
    AliWarning("Cut: fFlowEtaBinNum not set. Setting automatically with 0.05 bin width!");
    fFlowEtaBinNum = 2.0*fFlowEtaMax / 0.05;
  }

  // phi
  if(fFlowPhiBinNum < 1)
  {
    AliFatal("Cut: PhiBinNum wrong! Terminating!");
    return kFALSE;
  }

  // Centraltity // setting default values for centrality estimators
  if(fCentEstimator == kRFP)
  {
    if(fCentMax < 1) { fCentMax = 200; }
    if(fCentBinNum < 1) { fCentBinNum = 200; }
  }
  else
  {
    if(fCentMax < 1) { fCentMax = 100; }
    if(fCentBinNum < 1) { fCentBinNum = 100; }
  }

  if(fCentMin < 0 || fCentMax < 1 || fCentMin > fCentMax)
  {
    AliFatal("Centrality: range wrong! Terminating!");
    return kFALSE;
  }

  if(fCentBinNum < 1)
  {
    AliFatal("Centrality: Number of bins wrong! Terminating!");
    return kFALSE;
  }

  // setting procesing Refs & Charged by default
  fProcessSpec[kRefs] = kTRUE;
  fProcessSpec[kCharged] = kTRUE;

  // setting processing Kaons if Phi is on
  if(fProcessSpec[kPhi] && !fProcessSpec[kKaon])
  {
    fProcessSpec[kKaon] = kTRUE;
    AliInfo("Processing of Phi ON but PID OFF: setting processing of Kaons ON");
  }

  // checking for weights source file
  if(fFlowUseWeights)
  {
    if(fFlowFillWeights) { AliFatal("Cannot fill and run weights at the same time"); return kFALSE; }
    // BUG currently two pointer arrays overlay with each other, to-be-fixed

    fFlowWeightsList = (TList*) GetInputData(1);
    if(!fFlowRunByRunWeights && !LoadWeights()) { AliFatal("Initial flow weights not loaded! Terminating!"); return kFALSE; }
  }

  AliInfo("Preparing particle containers (std::vectors)");
  // creating particle vectors & reserving capacity in order to avoid memory re-allocation
  Int_t iReserve = 0;
  switch(fColSystem)
  {
    case kPP :
      iReserve = 150;
      break;

    case kPPb :
      iReserve = 200;
      break;

    default :
      iReserve = 300;
      break;
  }

  for(Int_t iSpec(0); iSpec < kUnknown; ++iSpec)
  {
    if(!fProcessSpec[iSpec]) { continue; }
    fVector[iSpec] = new std::vector<AliVTrack*>();
    fVector[iSpec]->reserve(iReserve);
  }

  AliInfo("Initialization succesfull!");
  return kTRUE;
}
// ============================================================================
void AliAnalysisTaskUniFlow::UserExec(Option_t *)
{
  // main method called for each event (event loop)
  // *************************************************************
  DumpTObjTable("UserExec: start");

  // check if initialization succesfull (done within UserCreateOutputObjects())
  if(!fInit) { AliFatal("Something went wrong : task not initialized!"); return; }

  // local event counter check: if running in test mode, it runs until the 50 events are succesfully processed
  if(fRunMode == kTest && fEventCounter >= fNumEventsAnalyse) { return; }

  // event selection
  fEventAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fEventAOD) { return; }

  // loading AliPIDResponse
  AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();
  if(!fPIDResponse) { AliFatal("AliPIDResponse not attached!"); return; }

  // loading array with MC particles
  if(fMC) {
    fEventMC = inputHandler->MCEvent();
    if(!fEventMC) { AliFatal("fEventMC with MC particle not found!"); return; }
  }

  // "valid" events before selection
  fhEventCounter->Fill("Input",1);

  // Fill event QA BEFORE cuts
  if(fFillQA) { FillQAEvents(kBefore); }

  Bool_t bEventSelected = IsEventSelected();
  // if(!IsEventSelected()) { return; }

  DumpTObjTable("UserExec: after event selection");
  if(!bEventSelected) { return; }

  fhEventCounter->Fill("Event OK",1);

  // checking the run number for aplying weights & loading TList with weights
  if(fFlowUseWeights && fFlowRunByRunWeights && fRunNumber != fEventAOD->GetRunNumber() && !LoadWeights()) { AliFatal("Weights not loaded!"); return; }

  DumpTObjTable("UserExec: before filtering");

  // Filter charged (& Refs) particles to evaluate event multiplcity / N_RFPs
  // NB: clear charged vectors because it might keep particles from previous event (not happen for other species)
  fVector[kRefs]->clear();
  fVector[kCharged]->clear();
  FilterCharged();

  // checking if there is at least 9 particles: needed to "properly" calculate correlations
  if(fVector[kRefs]->size() < 9) { return; }

  // estimate centrality & assign indexes (centrality/percentile, sampling, ...)
  if(fCentEstimator == kRFP) {
      fIndexCentrality = GetCentralityIndex(fCentEstimator);
      if(fIndexCentrality < 0) { return; }

      if(fCentMin > 0 && fIndexCentrality < fCentMin) { return; }
      if(fCentMax > 0 && fIndexCentrality > fCentMax) { return; }
  }

  if(fCentMaxAdd > 0 && fCentMinAdd > 0 && fCentEstimatorAdd == kRFP) {
    Int_t addCent = GetCentralityIndex(fCentEstimatorAdd);
    if(addCent < fCentMinAdd) { return; }
    if(addCent > fCentMaxAdd) { return; }
  }

  fhEventCounter->Fill("#RPFs OK",1);

  // here events are selected
  fhEventCounter->Fill("Selected",1);

  // event sampling
  fIndexSampling = GetSamplingIndex();

  // extract PV-z for weights
  fPVz = fEventAOD->GetPrimaryVertex()->GetZ();

  // Fill QA AFTER cuts (i.e. only in selected events)
  if(fFillQA) { FillQAEvents(kAfter); }

  // filling Charged QA histos
  // NB: for other species done within Filter*(): expection since # of Refs is part of event selection
  for (auto part = fVector[kRefs]->begin(); part != fVector[kRefs]->end(); part++) {
    fhChargedCounter->Fill("Refs",1);
    if(fFillQA) { FillQARefs(kAfter,static_cast<AliAODTrack*>(*part)); }
    if(!FillFlowWeight(*part, kRefs)) { AliFatal("Flow weight filling failed!"); return; }
  }

  for (auto part = fVector[kCharged]->begin(); part != fVector[kCharged]->end(); part++) {
    fhChargedCounter->Fill("POIs",1);
    if(fFillQA) { FillQACharged(kAfter,static_cast<AliAODTrack*>(*part)); } // QA after selection
    if(!FillFlowWeight(*part, kCharged)) { AliFatal("Flow weight filling failed!"); return; }
  }

  if(fFillQA) {
    // Charged QA before selection
    for(Int_t iTrack(0); iTrack < fEventAOD->GetNumberOfTracks(); iTrack++) {
      AliAODTrack* track = static_cast<AliAODTrack*>(fEventAOD->GetTrack(iTrack));
      if(!track) { continue; }
      FillQACharged(kBefore,track);
    }

    fhQAChargedMult[0]->Fill(fEventAOD->GetNumberOfTracks());
    fhQAChargedMult[1]->Fill(fVector[kCharged]->size());
    fhRefsMult->Fill(fVector[kRefs]->size());
  }

  // sorting charged hadrons
  std::sort(fVector[kCharged]->begin(), fVector[kCharged]->end(), [this](const AliVTrack* a, const AliVTrack* b){ return this->sortPt(a, b); });

  // Filtering other species
  if(fProcessSpec[kPion] || fProcessSpec[kKaon] || fProcessSpec[kProton]) { FilterPID(); }
  if(fProcessSpec[kK0s] || fProcessSpec[kLambda]) {
      FilterV0s();
      std::sort(fVector[kK0s]->begin(), fVector[kK0s]->end(), [this](const AliVTrack* a, const AliVTrack* b){ return this->sortPt(a, b); });
      std::sort(fVector[kLambda]->begin(), fVector[kLambda]->end(), [this](const AliVTrack* a, const AliVTrack* b){ return this->sortPt(a, b); });
  }
  if(fProcessSpec[kPhi]) {
      FilterPhi();
      std::sort(fVector[kPhi]->begin(), fVector[kPhi]->end(), [this](const AliVTrack* a, const AliVTrack* b){ return this->sortPt(a, b); });
  }

  DumpTObjTable("UserExec: after filtering");

  // processing of selected event
  Bool_t bProcessed = CalculateFlow();

  DumpTObjTable("UserExec: after CalculateFlow");

  // should be cleared at the end of processing especially for reconstructed
  // particles (Phi, V0s) because here new AliPicoTracks are created
  ClearVectors();

  if(fMC) { ProcessMC(); }

  // extracting run number here to store run number from previous event (for current run number use info in AliAODEvent)
  fRunNumber = fEventAOD->GetRunNumber();

  DumpTObjTable("UserExec: end");

  if(!bProcessed) return;

  // posting data (mandatory)
  Int_t i = 0;
  for(Int_t iSpec(0); iSpec < kUnknown; ++iSpec) { PostData(++i, fListFlow[iSpec]); }
  PostData(++i, fQAEvents);
  PostData(++i, fQACharged);
  PostData(++i, fQAPID);
  PostData(++i, fQAV0s);
  PostData(++i, fQAPhi);
  PostData(++i, fFlowWeights);
  if(fMC) { PostData(++i, fListMC); }

  return;
}
// ============================================================================
Bool_t AliAnalysisTaskUniFlow::IsEventSelected()
{
  // Event selection for pp & p-Pb collision recorded in Run 2 using AliEventCuts
  // return kTRUE if event passes all criteria, kFALSE otherwise
  // *************************************************************

  // Physics selection (trigger)
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) mgr->GetInputEventHandler();
  UInt_t fSelectMask = inputHandler->IsEventSelected();
  if(!(fSelectMask & fTrigger)) { return kFALSE; }

  // events passing physics && trigger selection
  fhEventCounter->Fill("Physics selection OK",1);

  // events passing AliEventCuts selection
  if(!fEventCuts.AcceptEvent(fEventAOD)) { return kFALSE; }
  fhEventCounter->Fill("EventCuts OK",1);

  // Additional pile-up rejection cuts for LHC15o dataset
  if(fColSystem == kPbPb && fEventRejectAddPileUp && IsEventRejectedAddPileUp()) { return kFALSE; }

  // estimate centrality & assign indexes (only if AliMultEstimator is requested)
  if(fCentEstimator != kRFP) {
      fIndexCentrality = GetCentralityIndex(fCentEstimator);

      if(fIndexCentrality < 0) { return kFALSE; }
      if(fCentMin > 0 && fIndexCentrality < fCentMin) { return kFALSE; }
      if(fCentMax > 0 && fIndexCentrality > fCentMax) { return kFALSE; }
  }

  // additional centrality cut for "double differential" cut
  if(fCentEstimatorAdd != kRFP && fCentMaxAdd > 0) {
    Int_t addCentIndex = GetCentralityIndex(fCentEstimatorAdd);
    if(addCentIndex < fCentMinAdd) { return kFALSE; }
    if(addCentIndex > fCentMaxAdd) { return kFALSE; }
  }

  fhEventCounter->Fill("Cent/Mult OK",1);

  return kTRUE;
}
// ============================================================================
Bool_t AliAnalysisTaskUniFlow::IsEventRejectedAddPileUp() const
{
  // Check for additional pile-up rejection in Run 2 Pb-Pb collisions (15o, 17n)
  // based on multiplicity correlations
  // ***************************************************************************

  Bool_t bIs17n = kFALSE;
  Bool_t bIs15o = kFALSE;

  Int_t iRunNumber = fEventAOD->GetRunNumber();
  if(iRunNumber >= 244824 && iRunNumber <= 246994) { bIs15o = kTRUE; }
  else if(iRunNumber == 280235 || iRunNumber == 20234) { bIs17n = kTRUE; }
  else { return kFALSE; }

  // recounting multiplcities
  const Int_t multESD = ((AliAODHeader*) fEventAOD->GetHeader())->GetNumberOfESDTracks();
  const Int_t nTracks = fEventAOD->GetNumberOfTracks();
  Int_t multTPC32 = 0;
  Int_t multTPC128 = 0;
  Int_t multTOF = 0;
  Int_t multTrk = 0;
  Double_t multESDTPCdif = 0.0;
  Double_t v0Centr = 0.0;

  for(Int_t it(0); it < nTracks; it++)
  {
    AliAODTrack* track = (AliAODTrack*) fEventAOD->GetTrack(it);
    if(!track) { continue; }

    if(track->TestFilterBit(32))
    {
      multTPC32++;
      if(TMath::Abs(track->GetTOFsignalDz()) <= 10.0 && track->GetTOFsignal() >= 12000.0 && track->GetTOFsignal() <= 25000.0) { multTOF++; }
      if((TMath::Abs(track->Eta())) < fFlowEtaMax && (track->GetTPCNcls() >= fCutChargedNumTPCclsMin) && (track->Pt() >= fFlowRFPsPtMin) && (track->Pt() < fFlowRFPsPtMax)) { multTrk++; }
    }

    if(track->TestFilterBit(128)) { multTPC128++; }
  }

  if(bIs17n)
  {
    multESDTPCdif = multESD - (6.6164 + 3.64583*multTPC128 + 0.000126397*multTPC128*multTPC128);
    if(multESDTPCdif > 1000) { return kTRUE; }
    if( ((AliAODHeader*) fEventAOD->GetHeader())->GetRefMultiplicityComb08() < 0) { return kTRUE; }
  }

  if(bIs15o)
  {
    multESDTPCdif = multESD - 3.38*multTPC128;
    if(multESDTPCdif > 500) { return kTRUE; }

    TF1 fMultTOFLowCut = TF1("fMultTOFLowCut", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 4.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x+[9]*x*x*x*x*x)", 0, 10000);
    fMultTOFLowCut.SetParameters(-1.0178, 0.333132, 9.10282e-05, -1.61861e-08, 1.47848, 0.0385923, -5.06153e-05, 4.37641e-08, -1.69082e-11, 2.35085e-15);
    if(Double_t(multTOF) < fMultTOFLowCut.Eval(Double_t (multTPC32))) { return kTRUE; }

    TF1 fMultTOFHighCut = TF1("fMultTOFHighCut", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 4.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x+[9]*x*x*x*x*x)", 0, 10000);
    fMultTOFHighCut.SetParameters(-1.0178, 0.333132, 9.10282e-05, -1.61861e-08, 1.47848, 0.0385923, -5.06153e-05, 4.37641e-08, -1.69082e-11, 2.35085e-15);
    if(Double_t(multTOF) > fMultTOFHighCut.Eval(Double_t (multTPC32))) { return kTRUE; }

    AliMultSelection* multSelection = (AliMultSelection*) fEventAOD->FindListObject("MultSelection");
    if(!multSelection) { AliError("AliMultSelection object not found! Returning -1"); return -1; }
    v0Centr = multSelection->GetMultiplicityPercentile("V0M");

    TF1 fMultCentLowCut = TF1("fMultCentLowCut", "[0]+[1]*x+[2]*exp([3]-[4]*x) - 5.*([5]+[6]*exp([7]-[8]*x))", 0, 100);
    fMultCentLowCut.SetParameters(-6.15980e+02, 4.89828e+00, 4.84776e+03, -5.22988e-01, 3.04363e-02, -1.21144e+01, 2.95321e+02, -9.20062e-01, 2.17372e-02);
    if(Double_t(multTrk) < fMultCentLowCut.Eval(v0Centr)) { return kTRUE; }
  }

  // QA Plots
  if(fFillQA) {
    fhQAEventsfMult32vsCentr->Fill(v0Centr, multTrk);
    fhQAEventsMult128vsCentr->Fill(v0Centr, multTPC128);
    fhQAEventsfMultTPCvsTOF->Fill(multTPC32, multTOF);
    fhQAEventsfMultTPCvsESD->Fill(multTPC128, multESD);
  }

  // fCentralityDis->Fill(centrV0);
  // fV0CentralityDis->Fill(cent);
  //
  // fCentralityV0MCL1->Fill(v0Centr, cl1Centr);
  // fCentralityV0MCL0->Fill(v0Centr, cl0Centr);
  // fCentralityCL0CL1->Fill(cl0Centr, cl1Centr);

  return kFALSE;
}
// ============================================================================
Bool_t AliAnalysisTaskUniFlow::LoadWeights()
{
  // (Re-) Loading of flow vector weights
  // ***************************************************************************
  if(!fFlowWeightsList) { AliError("Flow weights list not found! Terminating!"); return kFALSE; }

  TList* listFlowWeights = nullptr;

  if(!fFlowWeightsTag.IsNull()) {
      // using weights Tag if provided (systematics)
      listFlowWeights = (TList*) fFlowWeightsList->FindObject(fFlowWeightsTag.Data());
      if(!listFlowWeights) { AliError(Form("TList with tag '%s' not found!",fFlowWeightsTag.Data())); fFlowWeightsList->ls(); return kFALSE; }
  } else {
      if(!fFlowRunByRunWeights) {
          // loading run-averaged weights
          listFlowWeights = (TList*) fFlowWeightsList->FindObject("averaged");
          if(!listFlowWeights) { AliError("TList with flow run-averaged weights not found."); fFlowWeightsList->ls(); return kFALSE; }
      } else {
          // loading run-specific weights
          listFlowWeights = (TList*) fFlowWeightsList->FindObject(Form("%d",fEventAOD->GetRunNumber()));

          if(!listFlowWeights) {
              // run-specific weights not found for this run; loading run-averaged instead
              AliWarning(Form("TList with flow weights (run %d) not found. Using run-averaged weights instead (as a back-up)", fEventAOD->GetRunNumber()));
              listFlowWeights = (TList*) fFlowWeightsList->FindObject("averaged");
              if(!listFlowWeights) { AliError("Loading run-averaged weights failed!"); fFlowWeightsList->ls(); return kFALSE; }
          }
      }
  }


  for(Int_t iSpec(0); iSpec < kUnknown; ++iSpec) {
    if(!fProcessSpec[iSpec]) { continue; }
    if(iSpec == kKaon && (!fProcessSpec[kPion] || !fProcessSpec[kProton])) { continue; }

    if(fFlowUse3Dweights) {
      fh3Weights[iSpec] = (TH3D*) listFlowWeights->FindObject(Form("%s3D",GetSpeciesName(PartSpecies(iSpec))));
      if(!fh3Weights[iSpec]) { AliError(Form("Weight 3D (%s) not found",GetSpeciesName(PartSpecies(iSpec)))); return kFALSE; }
    } else {
      fh2Weights[iSpec] = (TH2D*) listFlowWeights->FindObject(GetSpeciesName(PartSpecies(iSpec)));
      if(!fh2Weights[iSpec]) { AliError(Form("Weight 2D (%s) not found",GetSpeciesName(PartSpecies(iSpec)))); return kFALSE; }
    }
  }

  return kTRUE;
}
// ============================================================================
Bool_t AliAnalysisTaskUniFlow::FillFlowWeight(const AliVTrack* track, const PartSpecies species) const
{
  if(!track) { AliError("Track not exists!"); return kFALSE; }
  if(species == kUnknown) { AliError("Invalid species 'Unknown'!"); return kFALSE; }

  if(fFlowFillWeights) {
    if(fFlowUse3Dweights) {
      fh3Weights[species]->Fill(track->Phi(),track->Eta(),fPVz);
    } else {
      fh2Weights[species]->Fill(track->Phi(),track->Eta());
    }
  }

  if(fFlowUseWeights && fFlowFillAfterWeights) {
    Double_t weight = GetFlowWeight(track, species);

    if(fFlowUse3Dweights) {
      fh3AfterWeights[species]->Fill(track->Phi(),track->Eta(),fPVz,weight);
    } else {
      fh2AfterWeights[species]->Fill(track->Phi(),track->Eta(),weight);
    }
  }

  return kTRUE;
}
// ============================================================================
Double_t AliAnalysisTaskUniFlow::GetFlowWeight(const AliVTrack* track, const PartSpecies species) const
{
  // if not applying for reconstructed
  if(!fFlowWeightsApplyForReco && HasMass(species)) { return 1.0; }

  Double_t dWeight = 1.0;
  if(fFlowUse3Dweights) {
    Int_t iBin = fh3Weights[species]->FindFixBin(track->Phi(),track->Eta(),fPVz);
    dWeight = fh3Weights[species]->GetBinContent(iBin);
  } else {
    Int_t iBin = fh2Weights[species]->FindFixBin(track->Phi(),track->Eta());
    dWeight = fh2Weights[species]->GetBinContent(iBin);
  }

  if(dWeight <= 0.0) { dWeight = 1.0; }
  return dWeight;
}
// ============================================================================
void AliAnalysisTaskUniFlow::FillQAEvents(const QAindex iQAindex) const
{
  // Filling various QA plots related with event selection
  // *************************************************************

  if(iQAindex == 1)
  {
    fhEventCentrality->Fill(fIndexCentrality);
    fh2EventCentralityNumRefs->Fill(fIndexCentrality,fVector[kRefs]->size());
    fpRefsMult->Fill(fIndexCentrality,fVector[kRefs]->size(),1.0);
    fhEventSampling->Fill(fIndexCentrality,fIndexSampling);
  }

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
// ============================================================================
void AliAnalysisTaskUniFlow::ProcessMC() const
{
    if(!fEventMC) { AliError("fEventMC with MC particles not found!"); return; }

    const Int_t iNumTracksMC = fEventMC->GetNumberOfTracks();
    for(Int_t iTrackMC(0); iTrackMC < iNumTracksMC; ++iTrackMC) {

        AliAODMCParticle* trackMC = (AliAODMCParticle*) fEventMC->GetTrack(iTrackMC);
        if(!trackMC) { continue; }


        Double_t dEta = trackMC->Eta();
        Double_t dPt = trackMC->Pt();
        Bool_t bCharged = (trackMC->Charge() != 0);
        Bool_t bIsPhysPrimary = trackMC->IsPhysicalPrimary();
        Bool_t bIsWithinPOIs = IsWithinPOIs(trackMC);


        if(bCharged && bIsPhysPrimary && bIsWithinPOIs) {
            // if(IsWithinRefs(trackMC)) { fh2MCPtEtaGen[kRefs]->Fill(dPt, dEta); }
            // if(bIsWithinPOIs) { fh2MCPtEtaGen[kCharged]->Fill(dPt,dEta); }

            fh2MCPtEtaGen[kCharged]->Fill(dPt,dEta); // Fill for both Charged & Refs
        }

        if(!bIsWithinPOIs) { continue; }

        PartSpecies species = kUnknown;
        Int_t iPDG = TMath::Abs(trackMC->GetPdgCode());
        for(Int_t spec(kPion); spec < Int_t(kUnknown); ++spec) {
            if(iPDG == fPDGCode[spec]) { species = PartSpecies(spec); }
        }

        if(species == kUnknown) { continue; }

        if(species == kPhi && (trackMC->GetMother() == -1)) {
            fh2MCPtEtaGen[species]->Fill(dPt, dEta);
        } else {
            if(bIsPhysPrimary) { fh2MCPtEtaGen[species]->Fill(dPt, dEta); }
        }
    }
}
// ============================================================================
void AliAnalysisTaskUniFlow::FilterCharged() const
{
  // Filtering input charged tracks for POIs (stored in fVector[kCharged]) or RFPs (fVector[kRefs])
  // If track passes all requirements its pointer is pushed to relevant vector container
  // *************************************************************

  Int_t iNumTracks = fEventAOD->GetNumberOfTracks();
  if(iNumTracks < 1) { return; }

  for(Int_t iTrack(0); iTrack < iNumTracks; iTrack++) {
    AliAODTrack* track = static_cast<AliAODTrack*>(fEventAOD->GetTrack(iTrack));
    if(!track) { continue; }

    // passing reconstruction criteria
    if(!IsChargedSelected(track)) { continue; }

    // Checking if selected track is eligible for Ref. flow
    if(IsWithinRefs(track)) {
        fVector[kRefs]->push_back(track);
        // if(fMC) { fh2MCPtEtaReco[kRefs]->Fill(track->Pt(), track->Eta()); }
    }

    // Checking if selected track is within POIs pt,eta acceptance
    if(IsWithinPOIs(track)) {
        fVector[kCharged]->push_back(track);
        if(fMC) { fh2MCPtEtaReco[kCharged]->Fill(track->Pt(), track->Eta()); }
    }
  } // end-for {iTrack}

  return;
}
// ============================================================================
Bool_t AliAnalysisTaskUniFlow::IsChargedSelected(const AliAODTrack* track) const
{
  // Selection of charged track
  // returns kTRUE if track pass all requirements, kFALSE otherwise
  // *************************************************************
  if(!track) { return kFALSE; }
  fhChargedCounter->Fill("Input",1);

  // filter bit
  if( !track->TestFilterBit(fCutChargedTrackFilterBit) ) { return kFALSE; }
  fhChargedCounter->Fill("FB",1);

  // number of TPC clusters (additional check for not ITS-standalone tracks)
  if( track->GetTPCNcls() < fCutChargedNumTPCclsMin && fCutChargedTrackFilterBit != 2) { return kFALSE; }
  fhChargedCounter->Fill("#TPC-Cls",1);

  // track DCA coordinates
  // note AliAODTrack::XYZAtDCA() works only for constrained tracks
  Double_t dTrackXYZ[3] = {0.,0.,0.};
  Double_t dVertexXYZ[3] = {0.,0.,0.};
  Double_t dDCAXYZ[3] = {0.,0.,0.};
  if( fCutChargedDCAzMax > 0. || fCutChargedDCAxyMax > 0.)
  {
    const AliAODVertex* vertex = fEventAOD->GetPrimaryVertex();
    if(!vertex) { return kFALSE; } // event does not have a PV

    track->GetXYZ(dTrackXYZ);
    vertex->GetXYZ(dVertexXYZ);

    for(Short_t i(0); i < 3; i++) { dDCAXYZ[i] = dTrackXYZ[i] - dVertexXYZ[i]; }
  }

  if(fCutChargedDCAzMax > 0. && TMath::Abs(dDCAXYZ[2]) > fCutChargedDCAzMax) { return kFALSE; }
  fhChargedCounter->Fill("DCA-z",1);

  if(fCutChargedDCAxyMax > 0. && TMath::Sqrt(dDCAXYZ[0]*dDCAXYZ[0] + dDCAXYZ[1]*dDCAXYZ[1]) > fCutChargedDCAxyMax) { return kFALSE; }
  fhChargedCounter->Fill("DCA-xy",1);

  // track passing all criteria
  fhChargedCounter->Fill("Selected",1);
  return kTRUE;
}
// ============================================================================
Bool_t AliAnalysisTaskUniFlow::IsWithinRefs(const AliVParticle* track) const
{
  // Checking if (preselected) track fulfills acceptance criteria for RFPs
  // NOTE: This is not a standalone selection, but additional check for IsChargedSelected()
  // It is used to selecting RFPs out of selected charged tracks
  // OR for estimating autocorrelations for Charged & PID particles
  // *************************************************************

  if(fFlowEtaMax > 0.0 && TMath::Abs(track->Eta()) > fFlowEtaMax) { return kFALSE; }
  if(fFlowRFPsPtMin > 0.0 && track->Pt() < fFlowRFPsPtMin) { return kFALSE; }
  if(fFlowRFPsPtMax > 0.0 && track->Pt() > fFlowRFPsPtMax) { return kFALSE; }

  return kTRUE;
}
// ============================================================================
Bool_t AliAnalysisTaskUniFlow::IsWithinPOIs(const AliVParticle* track) const
{
  // Checking if (preselected) track fulfills acceptance criteria for POIs
  // *************************************************************

  if(fFlowEtaMax > 0.0 && TMath::Abs(track->Eta()) > fFlowEtaMax) { return kFALSE; }
  if(fFlowPOIsPtMin > 0.0 && track->Pt() < fFlowPOIsPtMin) { return kFALSE; }
  if(fFlowPOIsPtMax > 0.0 && track->Pt() > fFlowPOIsPtMax) { return kFALSE; }

  return kTRUE;
}
// ============================================================================
void AliAnalysisTaskUniFlow::FillSparseCand(THnSparse* sparse, const AliVTrack* track) const
{
  // Fill sparse histogram for inv. mass distribution of candidates (V0s,Phi)
  // *************************************************************
  if(fRunMode == kSkipFlow || fVecCorrTask.size() < 1) { return; } // no sparse required
  if(!sparse) { Error("THnSparse not valid!","FillSparseCand"); return; }
  if(!track) { Error("Track not valid!","FillSparseCand"); return; }

  Double_t dValues[SparseCand::kDim] = {0};
  dValues[SparseCand::kCent] = fIndexCentrality;
  dValues[SparseCand::kInvMass] = track->M();
  dValues[SparseCand::kPt] = track->Pt();
  dValues[SparseCand::kEta] = track->Eta();
  sparse->Fill(dValues);

  return;
}
// ============================================================================
void AliAnalysisTaskUniFlow::FillQARefs(const QAindex iQAindex, const AliAODTrack* track) const
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
// ============================================================================
void AliAnalysisTaskUniFlow::FillQACharged(const QAindex iQAindex, const AliAODTrack* track) const
{
  // Filling various QA plots related to charged track selection
  // *************************************************************
  if(!track) return;

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
// ============================================================================
void AliAnalysisTaskUniFlow::FilterV0s() const
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
  if(iNumV0s < 1) { return; }

  for(Int_t iV0(0); iV0 < iNumV0s; iV0++)
  {
    AliAODv0* v0 = static_cast<AliAODv0*>(fEventAOD->GetV0(iV0));
    if(!v0) { continue; }

    if(fFillQA) { FillQAV0s(kBefore,v0); } // QA BEFORE selection

    if(!IsV0Selected(v0)) { continue; }

    Bool_t bIsK0s = IsV0aK0s(v0);
    Short_t iIsLambda = IsV0aLambda(v0);
    if(!bIsK0s && iIsLambda == 0) { continue; }

    if(fFillQA) { FillQAV0s(kAfter,v0,bIsK0s,iIsLambda); } // QA AFTER selection

    if(bIsK0s)
    {
      iNumK0sSelected++;
      fhV0sCounter->Fill("K^{0}_{S}",1);
      if(fFillQA) { fhV0sInvMassK0s->Fill(v0->MassK0Short(),v0->MassLambda()); }
      AliPicoTrack* pico = new AliPicoTrack(v0->Pt(),v0->Eta(),v0->Phi(),v0->Charge(),0,0,0,0,0,0,v0->MassK0Short());
      fVector[kK0s]->push_back(pico);
      FillSparseCand(fhsCandK0s, pico);

      if(fMC) {
          fh2MCPtEtaReco[kK0s]->Fill(v0->Pt(), v0->Eta());
          if(CheckMCTruthReco(kK0s,v0,(AliAODTrack*)v0->GetDaughter(0),(AliAODTrack*)v0->GetDaughter(1))) { fh2MCPtEtaRecoTrue[kK0s]->Fill(v0->Pt(), v0->Eta()); }
      }

      if(!FillFlowWeight(v0, kK0s)) { AliFatal("Flow weight filling failed!"); return; }
    }

    if(iIsLambda == 1) // lambda
    {
      iNumLambdaSelected++;
      fhV0sCounter->Fill("#Lambda/#bar{#Lambda}",1);
      if(fFillQA) { fhV0sInvMassLambda->Fill(v0->MassK0Short(),v0->MassLambda()); }

      AliPicoTrack* pico = new AliPicoTrack(v0->Pt(),v0->Eta(),v0->Phi(),v0->Charge(),0,0,0,0,0,0,v0->MassLambda());
      fVector[kLambda]->push_back(pico);
      FillSparseCand(fhsCandLambda, pico);

      if(fMC) {
          fh2MCPtEtaReco[kLambda]->Fill(v0->Pt(), v0->Eta());
          if(CheckMCTruthReco(kLambda,v0,(AliAODTrack*)v0->GetDaughter(0),(AliAODTrack*)v0->GetDaughter(1))) { fh2MCPtEtaRecoTrue[kLambda]->Fill(v0->Pt(), v0->Eta()); }
      }

      if(!FillFlowWeight(v0, kLambda)) { AliFatal("Flow weight filling failed!"); return; }
    }

    if(iIsLambda == -1) // anti-lambda
    {
      iNumALambdaSelected++;
      fhV0sCounter->Fill("#Lambda/#bar{#Lambda}",1);
      if(fFillQA) { fhV0sInvMassLambda->Fill(v0->MassK0Short(),v0->MassAntiLambda()); }

      AliPicoTrack* pico = new AliPicoTrack(v0->Pt(),v0->Eta(),v0->Phi(),v0->Charge(),0,0,0,0,0,0,v0->MassAntiLambda());
      fVector[kLambda]->push_back(pico);
      FillSparseCand(fhsCandLambda, pico);

      if(fMC) {
          fh2MCPtEtaReco[kLambda]->Fill(v0->Pt(), v0->Eta());
          if(CheckMCTruthReco(kLambda,v0,(AliAODTrack*)v0->GetDaughter(0),(AliAODTrack*)v0->GetDaughter(1))) { fh2MCPtEtaRecoTrue[kLambda]->Fill(v0->Pt(), v0->Eta()); }
      }

      if(!FillFlowWeight(v0, kLambda)) { AliFatal("Flow weight filling failed!"); return; }
    }

    if(bIsK0s && iIsLambda != 0) { fhV0sCounter->Fill("K^{0}_{S} && #Lambda/#bar{#Lambda}",1); }

  } // end-for {v0}

  // fill QA multiplicity
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
// ============================================================================
Bool_t AliAnalysisTaskUniFlow::IsV0aK0s(const AliAODv0* v0) const
{
  // Topological reconstruction and selection of V0 candidates
  // specific for K0s candidates
  // return kTRUE if a candidate fulfill all requirements, kFALSE otherwise
  // *************************************************************
  if(!v0) { return kFALSE; }
  fhV0sCounterK0s->Fill("Input",1);

  // rapidity selection
  if(fCutV0sMotherRapMax > 0.0 && (TMath::Abs(v0->RapK0Short()) > fCutV0sMotherRapMax) ) { return kFALSE; }
  fhV0sCounterK0s->Fill("#it{y}",1);

  // inv. mass window
  Double_t dMass = v0->MassK0Short();
  if(dMass < fCutV0sInvMassK0sMin || dMass > fCutV0sInvMassK0sMax) { return kFALSE; }
  fhV0sCounterK0s->Fill("InvMass",1);

  // cosine of pointing angle (CPA)
  if(fCutV0sCPAK0sMin > 0.0)
  {
    Double_t dCPA = v0->CosPointingAngle(fEventAOD->GetPrimaryVertex());
    if(dCPA < fCutV0sCPAK0sMin) { return kFALSE; }
  }
  fhV0sCounterK0s->Fill("CPA",1);

  // Armenteros-Podolaski plot
  if(fCutV0sArmenterosAlphaK0sMin > 0.0)
  {
    Double_t dPtArm = v0->PtArmV0();
    Double_t dAlpha = v0->AlphaV0();
    if(dPtArm < (fCutV0sArmenterosAlphaK0sMin * TMath::Abs(dAlpha))) { return kFALSE; }
  }
  fhV0sCounterK0s->Fill("Armenteros-Podolanski",1);

  // proper life-time
  if(fCutV0sNumTauK0sMax > 0.0)
  {
    Double_t dPrimVtxCoor[3] = {0.0,0.0,0.0}; // primary vertex position {x,y,z}
    Double_t dSecVtxCoor[3] = {0.0,0.0,0.0}; // secondary vertex position {x,y,z}
    Double_t dDecayCoor[3] = {0.0,0.0,0.0}; // decay vector coor {xyz}
    AliAODVertex* primVtx = fEventAOD->GetPrimaryVertex();
    primVtx->GetXYZ(dPrimVtxCoor);
    v0->GetSecondaryVtx(dSecVtxCoor);

    for(Int_t i(0); i < 3; i++) { dDecayCoor[i] = dSecVtxCoor[i] - dPrimVtxCoor[i]; }

    // implementation in xy plane
    // Double_t dPropLife = ( (fPDGMass[kK0s] / v0->Pt()) * TMath::Sqrt(dDecayCoor[0]*dDecayCoor[0] + dDecayCoor[1]*dDecayCoor[1]) );
    Double_t dPropLife = ( (fPDGMass[kK0s] / (v0->P() + 1e-10) ) * TMath::Sqrt(dDecayCoor[0]*dDecayCoor[0] + dDecayCoor[1]*dDecayCoor[1] + dDecayCoor[2]*dDecayCoor[2]) );
    if(dPropLife > (fCutV0sNumTauK0sMax * 2.68)) { return kFALSE; }
  }
  fhV0sCounterK0s->Fill("c#tau",1);

  // Daughter PID
  if(fCutV0sK0sPionNumTPCSigmaMax > 0.0)
  {
    const AliAODTrack* daughterPos = (AliAODTrack*) v0->GetDaughter(0);
    const AliAODTrack* daughterNeg = (AliAODTrack*) v0->GetDaughter(1);

    if(!HasTrackPIDTPC(daughterPos) || !HasTrackPIDTPC(daughterNeg)) { return kFALSE; }

    if (daughterPos->GetTPCsignalN() < fCutV0sDaughterNumTPCClsPIDMin || daughterNeg->GetTPCsignalN() < fCutV0sDaughterNumTPCClsPIDMin) { return kFALSE; }
    Float_t nSigmaPiPos = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(daughterPos, AliPID::kPion));
    Float_t nSigmaPiNeg = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(daughterNeg, AliPID::kPion));
    if(nSigmaPiPos > fCutV0sK0sPionNumTPCSigmaMax || nSigmaPiNeg > fCutV0sK0sPionNumTPCSigmaMax) { return kFALSE; }
  }
  fhV0sCounterK0s->Fill("Daughters PID",1);

  // competing V0 rejection based on InvMass
  if(fCutV0sCrossMassRejection)
  {
    Double_t dMassLambda = v0->MassLambda();
    Double_t dMassALambda = v0->MassAntiLambda();

    // K0s candidate is within 10 MeV of (Anti)Lambda InvMass physSelTask
    if(TMath::Abs(dMassLambda - fPDGMass[kLambda]) < fCutV0sCrossMassCutK0s)
    {
      // in Lambda peak
      if(fFillQA) { fhV0sCompetingInvMassK0s->Fill(dMass,dMassLambda); }
      return kFALSE;
    }

    if(TMath::Abs(dMassALambda - fPDGMass[kLambda]) < fCutV0sCrossMassCutK0s)
    {
      // in Anti-Lambda peak
      if(fFillQA) { fhV0sCompetingInvMassK0s->Fill(dMass,dMassALambda); }
      return kFALSE;
    }
  }
  fhV0sCounterK0s->Fill("Competing InvMass",1);

  // passing all criteria
  fhV0sCounterK0s->Fill("Selected",1);
  return kTRUE;
}
// ============================================================================
Int_t AliAnalysisTaskUniFlow::IsV0aLambda(const AliAODv0* v0) const
{
  // Topological reconstruction and selection of V0 candidates
  // specific for Lambda candidates
  // return 0 if candidate does not fullfill any Lambda or Anti-Lambda requirements;
  // return 1 if a candidate fulfill all Lambda requirements;
  // return -1 if a candidate fullfill all Anti-Lambda requirements;
  // return 2 if a candidate fulfill all both Lambda & Anti-Lambda requirements
  // *************************************************************
  if(!v0) { return 0; }
  fhV0sCounterLambda->Fill("Input",1);

  // rapidity selection
  if(fCutV0sMotherRapMax > 0.0 && (TMath::Abs(v0->RapLambda()) > fCutV0sMotherRapMax)) { return 0; }
  fhV0sCounterLambda->Fill("#it{y}",1);

  // particle species dependent
  Bool_t bIsLambda = kFALSE;
  Bool_t bIsALambda = kFALSE;

  // inv. mass window
  Double_t dMassLambda = v0->MassLambda();
  Double_t dMassALambda = v0->MassAntiLambda();
  if( dMassLambda > fCutV0sInvMassLambdaMin && dMassLambda < fCutV0sInvMassLambdaMax) { bIsLambda = kTRUE; }
  if( dMassALambda > fCutV0sInvMassLambdaMin && dMassALambda < fCutV0sInvMassLambdaMax) { bIsALambda = kTRUE; }

  if(!bIsLambda && !bIsALambda) { return 0; }
  fhV0sCounterLambda->Fill("InvMass",1);

  // cosine of pointing angle (CPA)
  if(fCutV0sCPALambdaMin > 0.0)
  {
    Double_t dCPA = v0->CosPointingAngle(fEventAOD->GetPrimaryVertex());
    if( dCPA < fCutV0sCPALambdaMin ) { return 0; }
  }
  fhV0sCounterLambda->Fill("CPA",1);

  // Armenteros-Podolaski plot
  if(fCutV0sArmenterosAlphaLambdaMax > 0.0)
  {
    Double_t dPtArm = v0->PtArmV0();
    Double_t dAlpha = v0->AlphaV0();
    if(dPtArm > (fCutV0sArmenterosAlphaLambdaMax * TMath::Abs(dAlpha))) { return 0; }
  }
  fhV0sCounterLambda->Fill("Armenteros-Podolanski",1);

  // // Armenteros-Podolanski for candidates fullfilling both Lambda and Anti-Lambda selection
  // if(bIsLambda && bIsALambda)
  // {
  //   if(v0->AlphaV0() < 0.) bIsLambda = kFALSE;
  //   if(v0->AlphaV0() > 0.) bIsALambda = kFALSE;
  // }

  // proper life-time
  if(fCutV0sNumTauLambdaMax > 0.0)
  {
    Double_t dPrimVtxCoor[3] = {0.0,0.0,0.0}; // primary vertex position {x,y,z}
    Double_t dSecVtxCoor[3] = {0.0,0.0,0.0}; // secondary vertex position {x,y,z}
    Double_t dDecayCoor[3] = {0.0,0.0,0.0}; // decay vector coor {xyz}
    AliAODVertex* primVtx = fEventAOD->GetPrimaryVertex();
    primVtx->GetXYZ(dPrimVtxCoor);
    v0->GetSecondaryVtx(dSecVtxCoor);

    for(Int_t i(0); i < 3; i++) { dDecayCoor[i] = dSecVtxCoor[i] - dPrimVtxCoor[i]; }

    // implementation in xy plane
    // Double_t dPropLife = ( (fPDGMass[kLambda] / v0->Pt()) * TMath::Sqrt(dDecayCoor[0]*dDecayCoor[0] + dDecayCoor[1]*dDecayCoor[1]) );
    Double_t dPropLife = ((fPDGMass[kLambda] / (v0->P() + 1e-10) ) * TMath::Sqrt(dDecayCoor[0]*dDecayCoor[0] + dDecayCoor[1]*dDecayCoor[1] + dDecayCoor[2]*dDecayCoor[2]));
    if(dPropLife > (fCutV0sNumTauLambdaMax * 7.89) ) { return 0; }
  }
  fhV0sCounterLambda->Fill("c#tau",1);

  // daughter PID of Lambda Candidates
  if(fCutV0sLambdaProtonNumTPCSigmaMax > 0.0 || fCutV0sLambdaPionNumTPCSigmaMax > 0.0)
  {
    const AliAODTrack* trackDaughterPos = (AliAODTrack*) v0->GetDaughter(0); // positive charge
    const AliAODTrack* trackDaughterNeg = (AliAODTrack*) v0->GetDaughter(1); // negative charge

    Bool_t bIsPosOK = HasTrackPIDTPC(trackDaughterPos);
    Bool_t bIsNegOK = HasTrackPIDTPC(trackDaughterNeg);

    Float_t dSigmaPos = 999.9;
    Float_t dSigmaNeg = 999.9;

    if(fCutV0sLambdaPionNumTPCSigmaMax > 0.0) // check pions
    {
      if(bIsPosOK) dSigmaPos = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackDaughterPos, AliPID::kPion)); else dSigmaPos = 999.9;
      if(bIsNegOK) dSigmaNeg = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackDaughterNeg, AliPID::kPion)); else dSigmaNeg = 999.9;

      if(bIsLambda && (!bIsNegOK || dSigmaNeg > fCutV0sLambdaPionNumTPCSigmaMax)) { bIsLambda = kFALSE; }
      if(bIsALambda && (!bIsPosOK || dSigmaPos > fCutV0sLambdaPionNumTPCSigmaMax)) { bIsALambda = kFALSE; }
    }

    if(fCutV0sLambdaProtonNumTPCSigmaMax > 0.0) // check protons
    {
      if(bIsPosOK) dSigmaPos = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackDaughterPos, AliPID::kProton)); else dSigmaPos = 999.9;
      if(bIsNegOK) dSigmaNeg = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackDaughterNeg, AliPID::kProton)); else dSigmaNeg = 999.9;

      if(bIsLambda && (!bIsPosOK || dSigmaPos > fCutV0sLambdaProtonNumTPCSigmaMax)) { bIsLambda = kFALSE; }
      if(bIsALambda && (!bIsNegOK || dSigmaNeg > fCutV0sLambdaProtonNumTPCSigmaMax)) { bIsALambda = kFALSE; }
    }

    if(!bIsLambda && !bIsALambda) { return 0; }
  }
  fhV0sCounterLambda->Fill("Daughter PID",1);

  // Lambda(AntiLamda) candidate is within fCutV0sCrossMassCutLambda of K0s InvMass
  if(fCutV0sCrossMassRejection)
  {
    Double_t dMassK0s = v0->MassK0Short();
    if(TMath::Abs(dMassK0s - fPDGMass[kK0s]) < fCutV0sCrossMassCutLambda)
    {
      if(fFillQA && bIsLambda) { fhV0sCompetingInvMassLambda->Fill(dMassK0s,dMassLambda); }
      if(fFillQA && bIsALambda) { fhV0sCompetingInvMassLambda->Fill(dMassK0s,dMassALambda); }
      return 0;
    }
  }
  fhV0sCounterLambda->Fill("Competing InvMass",1);

  // passing all criteria
  fhV0sCounterLambda->Fill("Selected",1);

  if(bIsLambda && bIsALambda) { fhV0sCounterLambda->Fill("#Lambda && #bar{#Lambda}",1); return 2; } // both Lambda & Anti-Lambda
  if(bIsLambda) { fhV0sCounterLambda->Fill("only #Lambda",1); return 1; } // only Lambda
  if(bIsALambda) { fhV0sCounterLambda->Fill("only #bar{#Lambda}",1); return -1; } // only Anti-Lambda
  return 0;
}
// ============================================================================
Double_t AliAnalysisTaskUniFlow::GetRapidity(const Double_t mass, const Double_t Pt, const Double_t Eta) const
{
    Double_t rapid = TMath::Log( (TMath::Sqrt(mass*mass + Pt*Pt*TMath::CosH(Eta)*TMath::CosH(Eta)) + Pt*TMath::SinH(Eta)) / TMath::Sqrt(mass*mass + Pt*Pt) );
    return rapid;
}
// ============================================================================
AliAODMCParticle* AliAnalysisTaskUniFlow::GetMCParticle(const Int_t label) const
{
  if(!fEventMC) { AliError("fEventMC not found!"); return nullptr; }

  const Int_t labelAbs = TMath::Abs(label);
  // Negative label just indicate track with shared clustes, but otherwise should be used
  // absolute value has to be used
  // if(label < 0) { /*AliWarning("MC label negative");*/ return nullptr; }

  AliAODMCParticle* mcTrack = (AliAODMCParticle*) fEventMC->GetTrack(labelAbs);
  if(!mcTrack) { AliWarning("Corresponding MC track not found!"); return nullptr; }
  return mcTrack;
}
// ============================================================================
Bool_t AliAnalysisTaskUniFlow::CheckMCPDG(const AliVParticle* track, const Int_t iPDGCode) const
{
    if(!track) { AliError("Input track does not exists!"); return kFALSE; }

    AliAODMCParticle* trackMC = GetMCParticle(track->GetLabel());
    if(!trackMC) { return kFALSE; }

    Int_t iPDG = TMath::Abs(trackMC->GetPdgCode());
    return (iPDG == iPDGCode);
}
// ============================================================================
Bool_t AliAnalysisTaskUniFlow::CheckMCPDG(const AliVParticle* track, const PartSpecies species) const
{
    if(!track) { AliError("Input track does not exists!"); return kFALSE; }
    if(species == kUnknown) { AliError("Cannot get PDG code for 'Unknown' species!"); return kFALSE; }
    if(species == kRefs || species == kCharged) { AliError("Cannot get charged hadron code for 'Unknown' species!"); return kFALSE; }

    return CheckMCPDG(track, fPDGCode[species]);
}
// ============================================================================
Bool_t AliAnalysisTaskUniFlow::CheckMCTruthReco(const PartSpecies species, const AliVParticle* track, const AliVParticle* daughterPos, const AliVParticle* daughterNeg) const
{
    if(!track) { AliError("Input track does not exists!"); return kFALSE; }
    if(!HasMass(species)) {
        // pi,K,p
        return CheckMCPDG(track,species);
    } else {
        // reconstructed
        if(!daughterPos || !daughterNeg) { AliError("Input daughter track does not exists!"); return kFALSE; }

        AliAODMCParticle* daughterMCPos = GetMCParticle(daughterPos->GetLabel());
        AliAODMCParticle* daughterMCNeg = GetMCParticle(daughterNeg->GetLabel());
        if(!daughterMCPos || !daughterMCNeg) { return kFALSE; }

        // checking charges (should be opposite)
        if(daughterMCPos->Charge() == daughterMCNeg->Charge()) { return kFALSE; }

        // checking if daughter tracks are secondary
        // NB: Kaons from Phi are considered PhysPrimary !!!
        // if(daughterMCPos->IsPhysicalPrimary() || daughterMCNeg->IsPhysicalPrimary()) { return kFALSE; }

        // Checking mother label
        Int_t iLabelPosMother = daughterMCPos->GetMother();
        Int_t iLabelNegMother = daughterMCNeg->GetMother();
        if(iLabelPosMother != iLabelNegMother) { return kFALSE; }

        AliAODMCParticle* motherMC = GetMCParticle(iLabelPosMother);
        if(!motherMC) { return kFALSE; }

        // checking PDG mother
        Int_t iPDGMother = TMath::Abs(motherMC->GetPdgCode());
        if(iPDGMother != fPDGCode[species]) { return kFALSE; }

        // checking PDG daughters
        Int_t iPDGPos = TMath::Abs(daughterMCPos->GetPdgCode());
        Int_t iPDGNeg = TMath::Abs(daughterMCNeg->GetPdgCode());

        if(species == kK0s) {
            // K0s -> pi+ + pi-
            return (iPDGPos == fPDGCode[kPion] && iPDGNeg == fPDGCode[kPion]);
        }

        if(species == kLambda) {
            // Lambda -> p + pi-
            if(iPDGPos == fPDGCode[kProton] && iPDGNeg == fPDGCode[kPion]) { return kTRUE; }
            // bar{Lambda} -> bar{p} + pi+
            if(iPDGPos == fPDGCode[kPion] && iPDGNeg == fPDGCode[kProton]) { return kTRUE; }

            return kFALSE;
        }

        if(species == kPhi) {
            return (iPDGPos == fPDGCode[kKaon] && iPDGNeg == fPDGCode[kKaon]);
        }
    }

    return kFALSE;
}
// ============================================================================
Bool_t AliAnalysisTaskUniFlow::IsV0Selected(const AliAODv0* v0) const
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
  if(!IsWithinPOIs(v0)) { return kFALSE; }
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
// ============================================================================
void AliAnalysisTaskUniFlow::FillQAV0s(const QAindex iQAindex, const AliAODv0* v0, const Bool_t bIsK0s, const Int_t bIsLambda) const
{
  // Filling various QA plots related to V0 candidate selection
  // *************************************************************
  // checking mother & daughters
  if(!v0) { return; }
  AliAODTrack* trackDaughter[2] = {(AliAODTrack*) v0->GetDaughter(0), (AliAODTrack*) v0->GetDaughter(1)};
  if(!trackDaughter[0] || !trackDaughter[1]) { return; }

  // setting internal flags for Lambdas and Anti-Lambdas
  Bool_t bCandLambda = kTRUE;
  Bool_t bCandAntiLambda = kTRUE;

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
  Double_t dSecVtxCoor[3] = {0.0,0.0,0.0};
  v0->GetSecondaryVtx(dSecVtxCoor);
  Double_t dDecayRadius = TMath::Sqrt(dSecVtxCoor[0]*dSecVtxCoor[0] + dSecVtxCoor[1]*dSecVtxCoor[1]);
  fhQAV0sDecayRadius[iQAindex]->Fill(dDecayRadius);

  // mother kinematics
  fhQAV0sMotherPt[iQAindex]->Fill(v0->Pt());
  fhQAV0sMotherPhi[iQAindex]->Fill(v0->Phi());
  fhQAV0sMotherEta[iQAindex]->Fill(v0->Eta());

  // proper lifetime preparation (to be filled in particle dependent if scope)
  Double_t dPrimVtxCoor[3] = {0.0,0.0,0.0};
  Double_t dDecayCoor[3] = {0.0,0.0,0.0};
  AliAODVertex* primVtx2 = fEventAOD->GetPrimaryVertex();
  primVtx2->GetXYZ(dPrimVtxCoor);
  for(Int_t i(0); i < 2; i++) { dDecayCoor[i] = dSecVtxCoor[i] - dPrimVtxCoor[i]; }

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
    Double_t dMassPDGK0s = fPDGMass[kK0s];
    // Double_t dPropLifeK0s = ( (dMassPDGK0s / v0->Pt() + 1e-10) * TMath::Sqrt(dDecayCoor[0]*dDecayCoor[0] + dDecayCoor[1]*dDecayCoor[1]) );
    Double_t dPropLifeK0s = ( (dMassPDGK0s / (v0->P() + 1e-10) ) * TMath::Sqrt(dDecayCoor[0]*dDecayCoor[0] + dDecayCoor[1]*dDecayCoor[1] + dDecayCoor[2]*dDecayCoor[2]) );
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
    if(bCandLambda) { fhQAV0sArmenterosLambda[iQAindex]->Fill(v0->AlphaV0(), v0->PtArmV0()); }

    if(bCandAntiLambda) { fhQAV0sArmenterosALambda[iQAindex]->Fill(v0->AlphaV0(), v0->PtArmV0()); }

    // proper lifetime
    Double_t dMassPDGLambda = fPDGMass[kLambda];
    // Double_t dPropLifeLambda = ( (dMassPDGLambda / v0->Pt() + 1e-10) * TMath::Sqrt(dDecayCoor[0]*dDecayCoor[0] + dDecayCoor[1]*dDecayCoor[1]) );
    Double_t dPropLifeLambda = ( (dMassPDGLambda / (v0->P() + 1e-10) ) * TMath::Sqrt(dDecayCoor[0]*dDecayCoor[0] + dDecayCoor[1]*dDecayCoor[1] + dDecayCoor[2]*dDecayCoor[2]) );
    fhQAV0sNumTauLambda[iQAindex]->Fill(dPropLifeLambda);
  }

  AliPIDResponse::EDetPidStatus pidStatusTPC;
  AliPIDResponse::EDetPidStatus pidStatusTOF;
  UShort_t numTPCfindable = 0;
  Float_t numTPCcrossed = 0.0;

  // daughters properties
  AliAODVertex* prodVtxDaughter = nullptr;
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
// ============================================================================
void AliAnalysisTaskUniFlow::FilterPhi() const
{
  // Reconstruction and filtering of Phi meson candidates out of selected Kaon sample
  // If track passes all requirements, the relevant properties (pT, eta, phi) are stored
  // in FlowPart struct  and pushed to relevant vector container.
  // *************************************************************
  Int_t iNumKaons = (Int_t) fVector[kKaon]->size();
  // check if there are at least 2 selected kaons in event (minimum for a single phi reconstruction)
  if(iNumKaons < 2) { return; }

  // start Phi reconstruction
  Int_t iNumBG = 0;
  for(Int_t iKaon1(0); iKaon1 < iNumKaons; iKaon1++) {
    AliAODTrack* kaon1 = dynamic_cast<AliAODTrack*>(fVector[kKaon]->at(iKaon1));
    if(!kaon1) { continue; }

    for(Int_t iKaon2(iKaon1+1); iKaon2 < iNumKaons; iKaon2++) {
      AliAODTrack* kaon2 = dynamic_cast<AliAODTrack*>(fVector[kKaon]->at(iKaon2));
      if(!kaon2) { continue; }

      AliPicoTrack* mother = MakeMother(kaon1,kaon2);
      fhPhiCounter->Fill("Input",1);

      // filling QA BEFORE selection
      if(fFillQA) { FillQAPhi(kBefore,mother); }

      if(fCutPhiInvMassMin > 0. && mother->M() < fCutPhiInvMassMin) { delete mother; continue; }
      if(fCutPhiInvMassMax > 0. && mother->M() > fCutPhiInvMassMax) { delete mother; continue; }
      fhPhiCounter->Fill("InvMass",1);

      if(!IsWithinPOIs(mother))  { delete mother; continue; }
      fhPhiCounter->Fill("Acceptance",1);

      // mother (phi) candidate passing all criteria (except for charge)
      fhPhiCounter->Fill("Before charge",1);

      if(TMath::Abs(mother->Charge()) == 2) {
        // like-sign combination (background)
        fhPhiCounter->Fill("BG",1);
        FillSparseCand(fhsCandPhiBg, mother);
        iNumBG++;
      }

      if(mother->Charge() == 0) {
        // opposite-sign combination (signal+background)
        fhPhiCounter->Fill("Unlike-sign",1);
        FillSparseCand(fhsCandPhi, mother);
        fVector[kPhi]->push_back(mother);

        if(fMC) {
          fh2MCPtEtaReco[kPhi]->Fill(mother->Pt(), mother->Eta());
          if(CheckMCTruthReco(kPhi,mother,kaon1,kaon2)) { fh2MCPtEtaRecoTrue[kPhi]->Fill(mother->Pt(), mother->Eta()); }
        }

        if(!FillFlowWeight(mother, kPhi)) { AliFatal("Flow weight filling failed!"); return; }
      }

      // filling QA AFTER selection
      if(fFillQA) { FillQAPhi(kAfter,mother); }

    } // endfor {iKaon2} : second kaon
  } // endfor {iKaon1} : first Kaon

  // filling multiplicity distribution
  if(fFillQA) {
    fhPhiMult->Fill(fVector[kPhi]->size());
    fhPhiBGMult->Fill(iNumBG);
  }

  return;
}
// ============================================================================
AliPicoTrack* AliAnalysisTaskUniFlow::MakeMother(const AliAODTrack* part1, const AliAODTrack* part2) const
{
  // Reconstructing mother particle from two prongs and fill its properties.
  // return ptr to created mother particle
  // *************************************************************

  if(!part1 || !part2) { return nullptr; }

  // combining momenta
  TVector3 mom1 = TVector3( part1->Px(), part1->Py(), part1->Pz() );
  TVector3 mom2 = TVector3( part2->Px(), part2->Py(), part2->Pz() );
  TVector3 mom = mom1 + mom2;

  Byte_t iCharge = part1->Charge() + part2->Charge();

  // calculating inv. mass
  Double_t dMass = -999.;
  Double_t dE1 = TMath::Sqrt( mom1.Mag2() + TMath::Power(fPDGMass[kKaon],2) );
  Double_t dE2 = TMath::Sqrt( mom2.Mag2() + TMath::Power(fPDGMass[kKaon],2) );

  Double_t dMassSq = TMath::Power((dE1+dE2),2) - mom.Mag2();
  if(dMassSq >= 0.) dMass = TMath::Sqrt(dMassSq);

  // maving phi form [-pi,pi] -> [0,2pi] for consistency with other species
  Double_t dPhi = mom.Phi() + TMath::Pi();

  return new AliPicoTrack(mom.Pt(),mom.Eta(),dPhi,iCharge,0,0,0,0,0,0,dMass);
}
// ============================================================================
void AliAnalysisTaskUniFlow::FillQAPhi(const QAindex iQAindex, const AliPicoTrack* part) const
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
// ============================================================================
void AliAnalysisTaskUniFlow::FilterPID() const
{
  // Filtering input PID tracks (pi,K,p)
  // If track passes all requirements as defined in IsPIDSelected() (and species dependent),
  // the relevant properties (pT, eta, phi) are stored in FlowPart struct
  // and pushed to relevant vector container.
  // return kFALSE if any complications occurs
  // *************************************************************

  for(auto part = fVector[kCharged]->begin(); part != fVector[kCharged]->end(); ++part)
  {
    AliAODTrack* track = static_cast<AliAODTrack*>(*part);
    if(!track) { continue; }

    fhPIDCounter->Fill("Input",1);

    if(fFillQA) { FillQAPID(kBefore,track,kUnknown); } // filling QA for tracks before selection (but after charged criteria applied)

    // PID track selection (return most favourable species)
    PartSpecies species = IsPIDSelected(track);
    if(species != kPion && species != kKaon && species != kProton) { continue; }

    // check if only protons should be used
    if(fCutPIDUseAntiProtonOnly && species == kProton && track->Charge() == 1) { continue; }

    if(!fProcessSpec[species]) { continue; }

    fhPIDCounter->Fill("Selected",1);
    fhPIDCounter->Fill(GetSpeciesName(species),1);

    fVector[species]->push_back(track);
    if(fFillQA) { FillQAPID(kAfter,track,species); } // filling QA for tracks AFTER selection }

    if(fProcessSpec[kPion] && fProcessSpec[kKaon] && fProcessSpec[kProton]) { // NB: aka process PID (not just Kaons for Phi)
      if(!FillFlowWeight(track, species)) { AliFatal("Flow weight filling failed!"); return; }
    }

    if(fMC) {
      fh2MCPtEtaReco[species]->Fill(track->Pt(), track->Eta());
      if(CheckMCTruthReco(species,track)) { fh2MCPtEtaRecoTrue[species]->Fill(track->Pt(), track->Eta()); }
    }

  } // end-for {part}

  if(fFillQA)
  {
    if(fProcessSpec[kPion]) { fhPIDMult[0]->Fill(fVector[kPion]->size()); }
    if(fProcessSpec[kKaon]) { fhPIDMult[1]->Fill(fVector[kKaon]->size()); }
    if(fProcessSpec[kProton]) { fhPIDMult[2]->Fill(fVector[kProton]->size()); }
  }

  return;
}
// ============================================================================
AliAnalysisTaskUniFlow::PartSpecies AliAnalysisTaskUniFlow::IsPIDSelected(const AliAODTrack* track) const
{
  // Selection of PID tracks (pi,K,p) - track identification
  // Based on fCutUseBayesPID flag, either Bayes PID or nSigma cutting is used
  // returns AliAnalysisTaskUniFlow::PartSpecies enum : kPion, kKaon, kProton if any of this passed kUnknown otherwise
  // *************************************************************

  // checking detector statuses
  Bool_t bIsTPCok = HasTrackPIDTPC(track);
  Bool_t bIsTOFok = HasTrackPIDTOF(track);

  if(!bIsTPCok) { return kUnknown; }

  // Preparing nSigma/Bayes arrays
  Float_t dNumSigmaTPC[fPIDNumSpecies];
  Float_t dNumSigmaTOF[fPIDNumSpecies];
  Float_t dNumSigmaCombined[fPIDNumSpecies];
  Double_t dProbPID[fPIDNumSpecies];

  for(Int_t iSpec(0); iSpec < fPIDNumSpecies; ++iSpec) {
    dNumSigmaTPC[iSpec] = -99.0;
    dNumSigmaTOF[iSpec] = -99.0;
    dNumSigmaCombined[iSpec] = -99.0;
    dProbPID[iSpec] = -99.0;
  }


  for(Int_t iSpec(0); iSpec < fPIDNumSpecies; ++iSpec) {
    dNumSigmaTPC[iSpec] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::EParticleType(iSpec)));
    dNumSigmaTOF[iSpec] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::EParticleType(iSpec)));
  }

  if(fCutUseBayesPID) {
    UInt_t iDetUsed = fPIDCombined->ComputeProbabilities(track, fPIDResponse, dProbPID); // filling probabilities to dPropPID array

    // check which detector were used
    Bool_t bUsedTOF = iDetUsed & AliPIDResponse::kDetTOF;
    // Bool_t bUsedTPC = iDetUsed & AliPIDResponse::kDetTPC; // Not checked
    // printf("   Selected: TPC:%d && TOF:%d)\n",bUsedTPC,bUsedTOF);

    Double_t dMaxProb = TMath::MaxElement(fPIDNumSpecies,dProbPID);
    // printf("PID Prob: e %g | mu %g | pi %g | K %g | p %g ||| MAX %g \n",dProbPID[0],dProbPID[1],dProbPID[2],dProbPID[3],dProbPID[4],dMaxProb);
    // Double_t dSum = 0.0; for(Int_t i(0); i < fPIDNumSpecies; ++i) { dSum += dProbPID[i]; }
    // printf("dSum = %g\n",dSum);

    // TODO: think about: if Pion has maximum probibility < fCutBayesPIDPion, track is rejected -> is it good?
    for(Int_t iSpec(2); iSpec < fPIDNumSpecies; ++iSpec) {
      if(dMaxProb == dProbPID[iSpec] && dProbPID[iSpec] >= fCutPIDBayesMin[iSpec]) {
        if(fCutPIDnSigmaMax[iSpec] > 0.0) {
          if(dNumSigmaTPC[iSpec] > fCutPIDnSigmaMax[iSpec]) { return kUnknown; }
          if(bUsedTOF && dNumSigmaTOF[iSpec] > fCutPIDnSigmaMax[iSpec]) { return kUnknown; }
        }
        return PartSpecies(iSpec);
      } // end-if {fCutPIDnSigmaMax}
    } // end-for {iSpec}
  } else {

    const Double_t dPt = track->Pt();

    // TPC nSigma cuts
    if(dPt <= 0.4) {
      Float_t dMinSigmasTPC = TMath::MinElement(fPIDNumSpecies,dNumSigmaTPC);

      // electron rejection
      if(dMinSigmasTPC == dNumSigmaTPC[0] && dNumSigmaTPC[0] <= fCutPIDnSigmaTPCRejectElectron) { return kUnknown; }
      for(Int_t iSpec(2); iSpec < fPIDNumSpecies; ++iSpec) {
        if(dMinSigmasTPC == dNumSigmaTPC[iSpec] && dNumSigmaTPC[iSpec] <= fCutPIDnSigmaMax[iSpec]) { return PartSpecies(iSpec); }
      }
    }

    // combined TPC + TOF nSigma cuts
    // NB: for testing of efficiency removed the upper limmit for PID
    // if(dPt > 0.4 && dPt < 4.0)
    if(dPt > 0.4) {
      // discard candidates if no TOF is available if cut is on
      if(fCutPIDnSigmaCombinedTOFrejection && !bIsTOFok) { return kUnknown; }

      // calculating combined nSigmas
      for(Int_t i(0); i < fPIDNumSpecies; ++i) {
        if(bIsTOFok) { dNumSigmaCombined[i] = TMath::Sqrt(dNumSigmaTPC[i]*dNumSigmaTPC[i] + dNumSigmaTOF[i]*dNumSigmaTOF[i]); }
        else { dNumSigmaCombined[i] = dNumSigmaTPC[i]; }
      }

      Float_t dMinSigmasCombined = TMath::MinElement(fPIDNumSpecies,dNumSigmaCombined);

      // electron rejection
      if(dMinSigmasCombined == dNumSigmaCombined[0] && dNumSigmaCombined[0] <= fCutPIDnSigmaTPCRejectElectron) { return kUnknown; }
      for(Int_t iSpec(2); iSpec < fPIDNumSpecies; ++iSpec) {
        if(dMinSigmasCombined == dNumSigmaCombined[iSpec] && dNumSigmaCombined[iSpec] <= fCutPIDnSigmaMax[iSpec]) { return PartSpecies(iSpec); }
      }
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
// ============================================================================
void AliAnalysisTaskUniFlow::FillQAPID(const QAindex iQAindex, const AliAODTrack* track, const PartSpecies species) const
{
  // Filling various QA plots related to PID (pi,K,p) track selection
  // *************************************************************
  if(!track) { return; }

  if(!fPIDResponse || !fPIDCombined) { AliError("AliPIDResponse or AliPIDCombined object not found!"); return; }

  // TPC & TOF statuses & measures
  AliPIDResponse::EDetPidStatus pidStatusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, track);
  AliPIDResponse::EDetPidStatus pidStatusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, track);

  Bool_t bIsTPCok = HasTrackPIDTPC(track);
  Bool_t bIsTOFok = HasTrackPIDTOF(track);

  Double_t dBayesProb[fPIDNumSpecies];
  Float_t dNumSigmaTPC[fPIDNumSpecies];
  Float_t dNumSigmaTOF[fPIDNumSpecies];
  Double_t dTPCdEdxDelta[fPIDNumSpecies];
  Double_t dTOFbetaDelta[fPIDNumSpecies];

  for(Int_t iSpec(0); iSpec < fPIDNumSpecies; ++iSpec) {
    dBayesProb[iSpec] = -0.1;
    dNumSigmaTPC[iSpec] = -11.0;
    dNumSigmaTOF[iSpec] = -11.0;
    dTPCdEdxDelta[iSpec] = -999.9;
    dTOFbetaDelta[iSpec] = -999.9;
  }

  Double_t dP = track->P();
  Double_t dPt = track->Pt();

  // TPC dEdx
  Double_t dTPCdEdx = track->GetTPCsignal();

  // TOF beta
  Double_t dTOF[5];
  track->GetIntegratedTimes(dTOF);
  Double_t dTOFbeta = dTOF[0] / track->GetTOFsignal();

  // filling Bayesian PID probabilities to dBayesProb array
  UInt_t iDetUsed = fPIDCombined->ComputeProbabilities(track, fPIDResponse, dBayesProb);

  // check which detector were used
  Bool_t bUsedTPC = kFALSE;
  Bool_t bUsedTOF = kFALSE;

  if(fCutUseBayesPID) {
    bUsedTPC = (bIsTPCok && (iDetUsed & AliPIDResponse::kDetTPC));
    bUsedTOF = iDetUsed & AliPIDResponse::kDetTOF;
  } else {
    bUsedTPC = bIsTPCok;
    bUsedTOF = (bIsTOFok && dPt > 0.4);
  }


  if(iQAindex == 0 || bUsedTPC) {
      fhQAPIDTPCstatus[iQAindex]->Fill((Int_t) pidStatusTPC );
      fhQAPIDTPCdEdx[iQAindex]->Fill(track->P(), dTPCdEdx);

      for(Int_t iSpec(0); iSpec < fPIDNumSpecies; ++iSpec) {
        dNumSigmaTPC[iSpec] = fPIDResponse->NumberOfSigmasTPC(track, AliPID::EParticleType(iSpec));
        dTPCdEdxDelta[iSpec] = fPIDResponse->GetSignalDelta(AliPIDResponse::kTPC, track, AliPID::EParticleType(iSpec));
      }
  }
  if(iQAindex == 0 || bUsedTOF) {
      fhQAPIDTOFstatus[iQAindex]->Fill((Int_t) pidStatusTOF );
      fhQAPIDTOFbeta[iQAindex]->Fill(dP,dTOFbeta);

      for(Int_t iSpec(0); iSpec < fPIDNumSpecies; ++iSpec) {
        dNumSigmaTOF[iSpec] = fPIDResponse->NumberOfSigmasTOF(track, AliPID::EParticleType(iSpec));
        dTOFbetaDelta[iSpec] = fPIDResponse->GetSignalDelta(AliPIDResponse::kTOF, track, AliPID::EParticleType(iSpec));
      }
  }

  fh3QAPIDnSigmaTPCTOFPtPion[iQAindex]->Fill(dNumSigmaTPC[2],dNumSigmaTOF[2],track->Pt());
  fh3QAPIDnSigmaTPCTOFPtKaon[iQAindex]->Fill(dNumSigmaTPC[3],dNumSigmaTOF[3],track->Pt());
  fh3QAPIDnSigmaTPCTOFPtProton[iQAindex]->Fill(dNumSigmaTPC[4],dNumSigmaTOF[4],track->Pt());

  if(species == kUnknown) { return; }

  // Here only selected particles (and iQAindex == 1 by construction)

  Int_t iPID = species - 2; // NB: translation from PartSpecies to PID QA index

  fhPIDPt[iPID]->Fill(track->Pt());
  fhPIDPhi[iPID]->Fill(track->Phi());
  fhPIDEta[iPID]->Fill(track->Eta());
  fhPIDCharge[iPID]->Fill(track->Charge());

  // TODO: Potentially might be converted to fh2[fPIDNumSpecies][3]

  fh2PIDBayesElectron[iPID]->Fill(dPt,dBayesProb[0]);
  fh2PIDBayesMuon[iPID]->Fill(dPt,dBayesProb[1]);
  fh2PIDBayesPion[iPID]->Fill(dPt,dBayesProb[2]);
  fh2PIDBayesKaon[iPID]->Fill(dPt,dBayesProb[3]);
  fh2PIDBayesProton[iPID]->Fill(dPt,dBayesProb[4]);

  if(bUsedTPC) {
    fh2PIDTPCdEdx[iPID]->Fill(dPt,dTPCdEdx);
    fh2PIDTPCdEdxDelta[iPID]->Fill(dPt,dTPCdEdxDelta[species]);
    fh2PIDTPCnSigmaElectron[iPID]->Fill(dPt,dNumSigmaTPC[0]);
    fh2PIDTPCnSigmaMuon[iPID]->Fill(dPt,dNumSigmaTPC[1]);
    fh2PIDTPCnSigmaPion[iPID]->Fill(dPt,dNumSigmaTPC[2]);
    fh2PIDTPCnSigmaKaon[iPID]->Fill(dPt,dNumSigmaTPC[3]);
    fh2PIDTPCnSigmaProton[iPID]->Fill(dPt,dNumSigmaTPC[4]);
  }

  if(bUsedTOF) {
    fh2PIDTOFbeta[iPID]->Fill(dPt,dTOFbeta);
    fh2PIDTOFbetaDelta[iPID]->Fill(dPt,dTOFbetaDelta[species]);
    fh2PIDTOFnSigmaElectron[iPID]->Fill(dPt,dNumSigmaTOF[0]);
    fh2PIDTOFnSigmaMuon[iPID]->Fill(dPt,dNumSigmaTOF[1]);
    fh2PIDTOFnSigmaPion[iPID]->Fill(dPt,dNumSigmaTOF[2]);
    fh2PIDTOFnSigmaKaon[iPID]->Fill(dPt,dNumSigmaTOF[3]);
    fh2PIDTOFnSigmaProton[iPID]->Fill(dPt,dNumSigmaTOF[4]);
  }

  return;
}
// ============================================================================
Bool_t AliAnalysisTaskUniFlow::ProcessCorrTask(const AliUniFlowCorrTask* task)
{
    if(!task) { AliError("AliUniFlowCorrTask does not exists!"); return kFALSE; }
    // task->Print();

    Int_t iNumHarm = task->fiNumHarm;
    Int_t iNumGaps = task->fiNumGaps;

    if(iNumGaps > 1) { AliError("Too many gaps! Not implemented yet!"); return kFALSE; }
    if(iNumHarm > 8) { AliError("Too many harmonics! Not implemented yet!"); return kFALSE; }


    Double_t dGap = -1.0;
    if(iNumGaps > 0) { dGap = task->fdGaps[0]; }

    // Fill anyway -> needed for any correlations
    FillRefsVectors(task, dGap); // TODO might check if previous task uses different Gap and if so, not fill it

    for(Int_t iSpec(0); iSpec < kUnknown; ++iSpec) {
        AliDebug(2,Form("Processing species '%s'",GetSpeciesName(PartSpecies(iSpec))));

        if(iSpec == kRefs) {
            if(!task->fbDoRefs) { continue; }
            CalculateCorrelations(task, kRefs);
            continue;
        }

        // here-after only POIs survive (Refs are dealt with already)
        if(!task->fbDoPOIs) { continue; }
        if(!fProcessSpec[iSpec]) { continue; }

        // NB: skip flow if Kaons are used only for Phi (flow not needed) not as full PID
        if(iSpec == kKaon && (!fProcessSpec[kPion] || !fProcessSpec[kProton])) { continue; }

        // loading (generic) profile to acess axes and bins
        TH1* genProf = (TH1*) fListFlow[iSpec]->FindObject(Form("%s_Pos_sample0",task->fsName.Data()));
        if(!genProf) { AliError(Form("Generic Profile '%s' not found", task->fsName.Data())); fListFlow[iSpec]->ls(); return kFALSE; }

        TAxis* axisPt = genProf->GetYaxis();
        if(!axisPt) { AliError("Pt axis object not found!"); return kFALSE; }
        Int_t iNumPtBins = axisPt->GetNbins();

        TAxis* axisMass = nullptr;
        Int_t iNumMassBins = 1;

        // check for 'massive' species
        Bool_t bHasMass = HasMass(PartSpecies(iSpec));
        if(bHasMass) {
            axisMass = genProf->GetZaxis();
            if(!axisMass) { AliError("Mass axis object not found!"); return kFALSE; }
            iNumMassBins = axisMass->GetNbins();
        }

        Int_t indexStart = 0;

        Int_t iNumPart = fVector[iSpec]->size();
        Int_t iNumFilled = 0;

        for(Int_t iPt(1); iPt < iNumPtBins+1; ++iPt) {
            Int_t iNumInPtBin = -10;

            Double_t dPt = axisPt->GetBinCenter(iPt);
            Double_t dPtLow = axisPt->GetBinLowEdge(iPt);
            Double_t dPtHigh = axisPt->GetBinUpEdge(iPt);

            for(Int_t iMass(1); iMass < iNumMassBins+1; ++iMass) {

                Double_t dMass = 0.0;
                Double_t dMassLow = 0.0;
                Double_t dMassHigh = 0.0;

                if(bHasMass) {
                    dMass = axisMass->GetBinCenter(iMass);
                    dMassLow = axisMass->GetBinLowEdge(iMass);
                    dMassHigh = axisMass->GetBinUpEdge(iMass);
                }

                Int_t contIndexStart = indexStart;

                // filling POIs (P,S) flow vectors
                Int_t iFilledHere = FillPOIsVectors(task, dGap ,PartSpecies(iSpec), contIndexStart, iNumInPtBin, dPtLow, dPtHigh, dMassLow, dMassHigh);
                CalculateCorrelations(task, PartSpecies(iSpec),dPt,dMass);

                // updating counters with numbers from this step
                iNumFilled += iFilledHere;
                iNumInPtBin -= iFilledHere;

                // switching index when all masses were proccessed in given pt bin (if applicable)
                // so strating point shifts to first particle in next pt bin
                if(iNumInPtBin < 1 || iNumFilled >= iNumPart || iMass == iNumMassBins) {
                    indexStart = contIndexStart;
                    break;
                }

                // if(iNumFilled >= iNumPart) { break; }
            }  // end-for {iMass}

            if(iNumFilled >= iNumPart) { break; }

        } // end-for {iPt}
    } // end-for {iSpecies}

    return kTRUE;
}
// ============================================================================
void AliAnalysisTaskUniFlow::CalculateCorrelations(const AliUniFlowCorrTask* const task, const PartSpecies species, const Double_t dPt, const Double_t dMass) const
{
  if(!task) { AliError("AliUniFlowCorrTask does not exists!"); return; }
  if(species >= kUnknown) { AliError(Form("Invalid species: %s!", GetSpeciesName(species))); return; }

  Bool_t bHasGap = task->HasGap();
  Int_t iNumHarm = task->fiNumHarm;
  Bool_t bDiff = kTRUE;
  Bool_t etaCheck = kFALSE;
  if(species == kRefs) {
    bDiff = kFALSE;
    if(fEtaCheckRFP) etaCheck = kTRUE;
  }


  // results of correlations
  TComplex cNom = TComplex(0.0,0.0,kFALSE);
  TComplex cDenom = TComplex(0.0,0.0,kFALSE);
  TComplex cNomNeg = TComplex(0.0,0.0,kFALSE);
  TComplex cDenomNeg = TComplex(0.0,0.0,kFALSE);

  // calculating correlations
  switch(iNumHarm)
  {
    case 2 : {
      if(!bHasGap) { // no gap
        if(bDiff) {
          cDenom = TwoDiff(0,0);
          cNom = TwoDiff(task->fiHarm[0],task->fiHarm[1]);
        }
        else {
          cDenom = Two(0,0);
          cNom = Two(task->fiHarm[0],task->fiHarm[1]);
        }
      }
      else { // has gap
        if(bDiff) {
          cDenom = TwoDiffGapPos(0,0);
          cDenomNeg = TwoDiffGapNeg(0,0);
          cNom = TwoDiffGapPos(task->fiHarm[0],task->fiHarm[1]);
          cNomNeg = TwoDiffGapNeg(task->fiHarm[0],task->fiHarm[1]);
        }
        else {
          if(!etaCheck) {
            cDenom = TwoGap(0,0);
            cNom = TwoGap(task->fiHarm[0],task->fiHarm[1]);
          }
          else {
            cDenom = TwoPos(0,0);
            cDenomNeg = TwoNeg(0,0);
            cNom = TwoPos(task->fiHarm[0],task->fiHarm[1]);
            cNomNeg = TwoNeg(task->fiHarm[0],task->fiHarm[1]);
            }
          }
        }
        break;
      }

    case 3 : {
      if(!bHasGap) { // no gap
        if(bDiff) {
          cDenom = ThreeDiff(0,0,0);
          cNom = ThreeDiff(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2]);
        }
        else {
          cDenom = Three(0,0,0);
          cNom = Three(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2]);
        }
      }
      else { // has gap
        if(bDiff) {
          cDenom = ThreeDiffGapPos(0,0,0);
          cDenomNeg = ThreeDiffGapNeg(0,0,0);
          cNom = ThreeDiffGapPos(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2]);
          cNomNeg = ThreeDiffGapNeg(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2]);
        }
        else {
          AliWarning("ThreeGap() not implemented!");
          return;
          // cDenom = ThreeGap(0,0,0);
          // cNom = ThreeGap(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2]); }
        }
      }
      break;
    }

    case 4 : {
      if(!bHasGap) { // no gap
        if(bDiff) {
          cDenom = FourDiff(0,0,0,0);
          cNom = FourDiff(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3]);
        }
        else {
          cDenom = Four(0,0,0,0);
          cNom = Four(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3]);
        }
      }
      else { // has gap
        if(bDiff) {
          cDenom = FourDiffGapPos(0,0,0,0);
          cDenomNeg = FourDiffGapNeg(0,0,0,0);
          cNom = FourDiffGapPos(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3]);
          cNomNeg = FourDiffGapNeg(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3]);
        }
        else {
          if(!etaCheck) {
            cDenom = FourGap(0,0,0,0);
            cNom = FourGap(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3]);
          }
          else {
            cDenom = FourPos(0,0,0,0);
            cDenomNeg = FourNeg(0,0,0,0);
            cNom = FourPos(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3]);
            cNomNeg = FourNeg(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3]);
          }
        }
      }
      break;
    }

    case 5: {
      if(!bHasGap) { // no gap
        if(bDiff) {
          AliWarning("Differential 5-particle correlations not implemented!");
          return;
        }
        else {
          cDenom = Five(0,0,0,0,0);
          cNom = Five(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3],task->fiHarm[4]);
        }
      }
      else { // has gap
        AliWarning("5-particle correlations with eta gap not implemented!");
        return;
      }
      break;
    }

    case 6: {
      if(!bHasGap) { // no gap
        if(bDiff) {
          AliWarning("Differential 6-particle correlations with no gap not implemented!");
          return;
        }
        else {
          cDenom = Six(0,0,0,0,0,0);
          cNom = Six(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3],task->fiHarm[4],task->fiHarm[5]);
        }
      }
      else { // has gap
        if(bDiff) {
          cDenom = SixDiffGapPos(0,0,0,0,0,0);
          cDenomNeg = SixDiffGapNeg(0,0,0,0,0,0);
          cNom = SixDiffGapPos(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3],task->fiHarm[4],task->fiHarm[5]);
          cNomNeg = SixDiffGapNeg(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3],task->fiHarm[4],task->fiHarm[5]);
        }
        else {
          if(!etaCheck) {
            cDenom = SixGap(0,0,0,0,0,0);
            cNom = SixGap(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3],task->fiHarm[4],task->fiHarm[5]);
          }
          else {
            cDenom = SixPos(0,0,0,0,0,0);
            cDenomNeg = SixNeg(0,0,0,0,0,0);
            cNom = SixPos(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3],task->fiHarm[4],task->fiHarm[5]);
            cNomNeg = SixNeg(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3],task->fiHarm[4],task->fiHarm[5]);
          }
        }
      }
      break;
    }

    case 7: {
      if(!bHasGap) { // no gap
        if(bDiff) {
          AliWarning("Differential 7-particle correlations not implemented!");
          return;
        }
        else {
          cDenom = Seven(0,0,0,0,0,0,0);
          cNom = Seven(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3],task->fiHarm[4],task->fiHarm[5],task->fiHarm[6]);
        }
      }
      else { // has gap
        AliWarning("7-particle correlations with eta gap not implemented!");
        return;
      }
      break;
    }

    case 8: {
      if(!bHasGap) { // no gap
        if(bDiff) {
          AliWarning("Differential 8-particle correlations with no gap not implemented!");
          return;
        }
        else {
          cDenom = Eight(0,0,0,0,0,0,0,0);
          cNom = Eight(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3],task->fiHarm[4],task->fiHarm[5],task->fiHarm[6],task->fiHarm[7]);
        }
      }
      else { // has gap
        if(bDiff) {
          cDenom = EightDiffGapPos(0,0,0,0,0,0,0,0);
          cDenomNeg = EightDiffGapNeg(0,0,0,0,0,0,0,0);
          cNom = EightDiffGapPos(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3],task->fiHarm[4],task->fiHarm[5],task->fiHarm[6],task->fiHarm[7]);
          cNomNeg = EightDiffGapNeg(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3],task->fiHarm[4],task->fiHarm[5],task->fiHarm[6],task->fiHarm[7]);
        }
        else {
          if(!etaCheck) {
            cDenom = EightGap(0,0,0,0,0,0,0,0);
            cNom = EightGap(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3],task->fiHarm[4],task->fiHarm[5],task->fiHarm[6],task->fiHarm[7]);
          }
          else {
            cDenom = EightPos(0,0,0,0,0,0,0,0);
            cDenomNeg = EightNeg(0,0,0,0,0,0,0,0);
            cNom = EightPos(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3],task->fiHarm[4],task->fiHarm[5],task->fiHarm[6],task->fiHarm[7]);
            cNomNeg = EightNeg(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3],task->fiHarm[4],task->fiHarm[5],task->fiHarm[6],task->fiHarm[7]);
          }
        }
      }
      break;
    }


    default:
      return;
  }

  Double_t dNom = cNom.Re();
  Double_t dDenom = cDenom.Re();
  Double_t dNomNeg = cNomNeg.Re();
  Double_t dDenomNeg = cDenomNeg.Re();

  Double_t dValue = 0.0;
  Double_t dValueNeg = 0.0;

  Bool_t bFillPos = kFALSE;
  Bool_t bFillNeg = kFALSE;

  // check which profiles should be filled (POS/NEG)
  if(dDenom > 0.0) { bFillPos = kTRUE; dValue = dNom / dDenom; }
  if(bFillPos && TMath::Abs(dValue) > 1.0) { bFillPos = kFALSE; }

  if(bHasGap && dDenomNeg > 0.0) { bFillNeg = kTRUE; dValueNeg = dNomNeg / dDenomNeg; }
  if(bFillNeg && TMath::Abs(dValueNeg) > 1.0) { bFillNeg = kFALSE; }

  if(!bFillPos && !bFillNeg) { return; } // To save some CPU time

  // Filling corresponding profiles
  switch(species)
  {
    case kRefs:
    {
      if(bFillPos)
      {
        TProfile* prof = (TProfile*) fListFlow[species]->FindObject(Form("%s_Pos_sample%d",task->fsName.Data(),fIndexSampling));
        if(!prof) { AliError(Form("Profile '%s_Pos_sample%d' not found!", task->fsName.Data(),fIndexSampling)); return; }
        prof->Fill(fIndexCentrality, dValue, dDenom);
      }
      if(bFillNeg)
      {
        TProfile* profNeg = (TProfile*) fListFlow[species]->FindObject(Form("%s_Neg_sample%d",task->fsName.Data(),fIndexSampling));
        if(!profNeg) { AliError(Form("Profile '%s_Neg_sample%d' not found!", task->fsName.Data(),fIndexSampling)); return; }
        profNeg->Fill(fIndexCentrality, dValueNeg, dDenomNeg);
      }
      break;
    }

    case kCharged:
    case kPion:
    case kKaon:
    case kProton:
    {
      if(bFillPos)
      {
        TProfile2D* prof = (TProfile2D*) fListFlow[species]->FindObject(Form("%s_Pos_sample%d",task->fsName.Data(),fIndexSampling));
        if(!prof) { AliError(Form("Profile '%s_Pos_sample%d' not found!", task->fsName.Data(),fIndexSampling)); return; }
        prof->Fill(fIndexCentrality, dPt, dValue, dDenom);
      }

      if(bFillNeg)
      {
        TProfile2D* profNeg = (TProfile2D*) fListFlow[species]->FindObject(Form("%s_Neg_sample%d",task->fsName.Data(),fIndexSampling));
        if(!profNeg) { AliError(Form("Profile '%s_Neg_sample%d' not found!", task->fsName.Data(),fIndexSampling)); return; }
        profNeg->Fill(fIndexCentrality, dPt, dValueNeg, dDenomNeg);
      }
      break;
    }


    case kK0s:
    case kLambda:
    case kPhi:
    {
      if(bFillPos)
      {
        TProfile3D* prof = (TProfile3D*) fListFlow[species]->FindObject(Form("%s_Pos_sample0",task->fsName.Data()));
        if(!prof) { AliError(Form("Profile '%s_Pos_sample0' not found!", task->fsName.Data())); return; }
        prof->Fill(fIndexCentrality, dPt, dMass, dValue, dDenom);
      }

      if(bFillNeg)
      {
        TProfile3D* profNeg = (TProfile3D*) fListFlow[species]->FindObject(Form("%s_Neg_sample0",task->fsName.Data()));
        if(!profNeg) { AliError(Form("Profile '%s_Neg_sample0' not found!", task->fsName.Data())); return; }
        profNeg->Fill(fIndexCentrality, dPt, dMass, dValueNeg, dDenomNeg);
      }
      break;
    }

    case kUnknown:
      return;
  }
  return;
}
// ============================================================================
Bool_t AliAnalysisTaskUniFlow::CalculateFlow()
{
  // main (envelope) method for flow calculations in selected events
  // returns kFALSE if something failes (with error), kTRUE otherwise
  // *************************************************************

  // if running in kSkipFlow mode, skip the remaining part
  if(fRunMode == kSkipFlow) { fEventCounter++; return kTRUE; }

  // >>>> Using AliUniFlowCorrTask <<<<<

  Int_t iNumTasks = fVecCorrTask.size();
  for(Int_t iTask(0); iTask < iNumTasks; ++iTask)
  {
    Bool_t process = ProcessCorrTask(fVecCorrTask.at(iTask));
    if(!process) { AliError("AliUniFlowCorrTask processing failed!\n"); fVecCorrTask.at(iTask)->Print(); return kFALSE; }
  }

  fEventCounter++; // counter of processed events

  return kTRUE;
}
// ============================================================================
void AliAnalysisTaskUniFlow::FillRefsVectors(const AliUniFlowCorrTask* task, const Double_t dGap)
{
  // Filling Q flow vector with RFPs
  // return kTRUE if succesfull (i.e. no error occurs), kFALSE otherwise
  // *************************************************************
  Double_t dEtaGap = dGap;
  Double_t dEtaLimit = dEtaGap / 2.0;
  Bool_t bHasGap = kFALSE;
  Bool_t bHas3sub = kFALSE;
  if(dEtaGap > -1.0) { bHasGap = kTRUE; }

  Int_t maxHarm = task->fMaxHarm;
  Int_t maxWeightPower = task->fMaxWeightPower;

  // clearing output (global) flow vectors
  ResetFlowVector(fFlowVecQpos, maxHarm, maxWeightPower);
  ResetFlowVector(fFlowVecQneg, maxHarm, maxWeightPower);
  if(bHas3sub) { ResetFlowVector(fFlowVecQmid, maxHarm, maxWeightPower); }

  for (auto part = fVector[kRefs]->begin(); part != fVector[kRefs]->end(); part++)
  {
    Double_t dPhi = (*part)->Phi();
    Double_t dEta = (*part)->Eta();

    if(bHasGap && TMath::Abs(dEta) < dEtaLimit) { continue; }

    // loading weights if needed
    Double_t dWeight = 1.0;
    if(fFlowUseWeights) { dWeight = GetFlowWeight(*part, kRefs); }

    if(!bHasGap) // no eta gap
    {
      for(Int_t iHarm(0); iHarm <= maxHarm; iHarm++)
        for(Int_t iPower(0); iPower <= maxWeightPower; iPower++)
        {
          Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
          Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
          fFlowVecQpos[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
        }
    }
    else
    {
      // RFP in positive eta acceptance
      if(dEta > dEtaLimit)
      {
        for(Int_t iHarm(0); iHarm <= maxHarm; iHarm++)
          for(Int_t iPower(0); iPower <= maxWeightPower; iPower++)
          {
            Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
            Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
            fFlowVecQpos[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
          }
      }
      // RFP in negative eta acceptance
      if(dEta < -dEtaLimit)
      {
        for(Int_t iHarm(0); iHarm <= maxHarm; iHarm++)
          for(Int_t iPower(0); iPower <= maxWeightPower; iPower++)
          {
            Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
            Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
            fFlowVecQneg[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
          }
      }

      // RFP in middle (for 3sub) if gap > 0
      if(bHas3sub && (TMath::Abs(dEta) < dEtaLimit) )
      {
        for(Int_t iHarm(0); iHarm <= maxHarm; iHarm++)
          for(Int_t iPower(0); iPower <= maxWeightPower; iPower++)
          {
            Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
            Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
            fFlowVecQmid[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
          }
      }

    } // endif {dEtaGap}
  } // endfor {tracks} particle loop

  return;
}
// ============================================================================
Int_t AliAnalysisTaskUniFlow::FillPOIsVectors(const AliUniFlowCorrTask* task, const Double_t dEtaGap, const PartSpecies species, Int_t& indStart, Int_t& tracksInBin, const Double_t dPtLow, const Double_t dPtHigh, const Double_t dMassLow, const Double_t dMassHigh)
{
  // Filling p,q and s flow vectors with POIs (given by species) for differential flow calculation
  // *************************************************************
  std::vector<AliVTrack*>* vector = fVector[species];
  if(!vector) { AliError("Vector with selected POIs not found."); return 0; }

  Double_t dEtaLimit = dEtaGap / 2.0;
  Bool_t bHasGap = kFALSE; if(dEtaGap > -1.0) { bHasGap = kTRUE; }
  Bool_t bHasMass = HasMass(species);
  if(bHasMass && dMassLow == 0.0 && dMassHigh == 0.0) { AliError("Particle mass low && high limits not specified!"); return 0; }

  Int_t maxHarm = task->fMaxHarm;
  Int_t maxWeightPower = task->fMaxWeightPower;

  // clearing output (global) flow vectors
  ResetFlowVector(fFlowVecPpos, maxHarm, maxWeightPower);
  ResetFlowVector(fFlowVecSpos, maxHarm, maxWeightPower);

  if(bHasGap) {
    ResetFlowVector(fFlowVecPneg, maxHarm, maxWeightPower);
    ResetFlowVector(fFlowVecSneg, maxHarm, maxWeightPower);
  }

  Int_t iTracksFilled = 0; // counter of filled tracks
  Int_t iTracksInPtBin = 0; // counter for all tracks in pt bins

  // for(auto part = vector->begin(); part != vector->end(); ++part)
  for(Int_t index(indStart); index < (Int_t) vector->size(); ++index) {
    AliVTrack* part = vector->at(index);
    if(!part) { AliError("Particle does not exists within given vector"); return -1; }

    Double_t dPhi = part->Phi();
    Double_t dEta = part->Eta();
    Double_t dPt = part->Pt();
    Double_t dMass = (bHasMass ? part->M() : 0.0);

    // checking if pt is within pt (bin) range
    if(dPt < dPtLow) { continue; }
    if(dPt >= dPtHigh) {
        // refresh the starting index value for next pt bin
        indStart = index;
        break;
        // return iTracksFilled;
    }

    iTracksInPtBin++;

    if(bHasMass) {
        if(dMass < dMassLow || dMass >= dMassHigh) { continue; }
    }

    // checking if mass is within mass (bin) range

    if(bHasGap && TMath::Abs(dEta) < dEtaLimit) { continue; }

    // at this point particles corresponding to this pt (& mass) bin and eta acceptance (gap) survives
    iTracksFilled++;

    // loading weights if needed
    Double_t dWeight = 1.0;
    if(fFlowUseWeights) { dWeight = GetFlowWeight(part, species); }

    // check if POI overlaps with RFPs (not for reconstructed)
    Bool_t bIsWithinRefs = (!bHasMass && IsWithinRefs(static_cast<const AliAODTrack*>(part)));

    if(!bHasGap) // no eta gap
    {
      for(Int_t iHarm(0); iHarm <= maxHarm; iHarm++)
        for(Int_t iPower(0); iPower <= maxWeightPower; iPower++)
        {
          Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
          Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
          fFlowVecPpos[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);

          // check if track (passing criteria) is overlapping with RFPs pT region; if so, fill S (q) vector
          // in case of charged, pions, kaons or protons (one witout mass)
          if(bIsWithinRefs)
          {
            Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
            Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
            fFlowVecSpos[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
          }
        }

    }
    else // with eta gap
    {
      if(dEta > dEtaLimit) // particle in positive eta acceptance
      {
        for(Int_t iHarm(0); iHarm <= maxHarm; iHarm++)
          for(Int_t iPower(0); iPower <= maxWeightPower; iPower++)
          {
            Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
            Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
            fFlowVecPpos[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);

            // possible overlap for <<4'>> with single gap (within the same subevent)
            if(bIsWithinRefs)
            {
              Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
              Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
              fFlowVecSpos[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
            }
          }
       }
       if(dEta < -dEtaLimit) // particle in negative eta acceptance
       {
         for(Int_t iHarm(0); iHarm <= maxHarm; iHarm++)
           for(Int_t iPower(0); iPower <= maxWeightPower; iPower++)
           {
             Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
             Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
             fFlowVecPneg[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);

             // possible overlap for <<4'>> with single gap (within the same subevent)
             if(bIsWithinRefs)
             {
               Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
               Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
               fFlowVecSneg[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
             }
           }
       }
     } // endif {dEtaGap}
   } // endfor {tracks}

   // refresh the value only if first go (aka initialized to -10); after that it is used as a counter of remaining particles
   if(tracksInBin < 0) { tracksInBin = iTracksInPtBin; }

   return iTracksFilled;
}
// ============================================================================
void AliAnalysisTaskUniFlow::ResetFlowVector(TComplex (&array)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax], Int_t maxHarm, Int_t maxWeightPower)
{
  // Reset RFPs (Q) array values to TComplex(0,0,kFALSE) for given array
  // *************************************************************
  for(Int_t iHarm(0); iHarm <= maxHarm; ++iHarm) {
    for(Int_t iPower(0); iPower <= maxWeightPower; ++iPower) {
      array[iHarm][iPower](0.0,0.0);
    }
  }
  return;
}
// ============================================================================
void AliAnalysisTaskUniFlow::ListFlowVector(TComplex (&array)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]) const
{
  // List all values of given flow vector TComplex array
  // *************************************************************
  printf(" ### Listing (TComplex) flow vector array ###########################\n");
  for(Int_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
  {
    printf("Harm %d (power):",iHarm);
    for(Int_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
    {
        printf("|(%d) %g+%g(i)",iPower, array[iHarm][iPower].Re(), array[iHarm][iPower].Im());
    }
    printf("\n");
  }
  return;
}
// ============================================================================
Int_t AliAnalysisTaskUniFlow::GetSamplingIndex() const
{
  // Assessing sampling index based on generated random number
  // returns centrality index
  // *************************************************************
  if(!fSampling) { return 0; }

  TRandom3 rr(0);
  Double_t ranNum = rr.Rndm(); // getting random number in (0,1)
  Double_t generated = ranNum * fNumSamples; // getting random number in range (0, fNumSamples)

  // finding right index for sampling based on generated number and total number of samples
  Int_t index = 0;
  for(Int_t i(0); i < fNumSamples; ++i) {
    if(generated < (i+1)) { index = i; break; }
  }

  return index;
}
// ============================================================================
Int_t AliAnalysisTaskUniFlow::GetCentralityIndex(CentEst est) const
{
  // Estimating centrality percentile based on selected estimator.
  // (Default) If no multiplicity estimator is specified, percentile is estimated as number of selected / filtered charged tracks (NRFP).
  // If a valid multiplicity estimator is specified, centrality percentile is estimated via AliMultSelection
  // otherwise -1 is returned (and event is skipped)
  // *************************************************************
  Int_t iCentralityIndex = -1;

  // assigning centrality based on number of selected charged tracks
  if(est == kRFP) {
      iCentralityIndex = fVector[kRefs]->size();
  } else {
    AliMultSelection* multSelection = (AliMultSelection*) fEventAOD->FindListObject("MultSelection");
    if(!multSelection) {
      AliError("AliMultSelection object not found! Returning -1");
      return -1;
    }

    Float_t dPercentile = multSelection->GetMultiplicityPercentile(GetCentEstimatorLabel(est));
    if(dPercentile > 100 || dPercentile < 0) {
      AliWarning("Centrality percentile estimated not within 0-100 range. Returning -1");
      return -1;
    }

    iCentralityIndex = (Int_t) dPercentile;
  }

  return iCentralityIndex;
}
// ============================================================================
const char* AliAnalysisTaskUniFlow::GetCentEstimatorLabel(const CentEst est) const
{
  // Return string with estimator name or 'n/a' if not available
  // *************************************************************
  switch (est) {
    case kRFP: return "N_{RFP}";
    case kV0A: return "V0A";
    case kV0C: return "V0C";
    case kV0M: return "V0M";
    case kCL0: return "CL0";
    case kCL1: return "CL1";
    case kZNA: return "ZNA";
    case kZNC: return "ZNC";
    default: return "n/a";
  }
  return "n/a";
}
// ============================================================================
Bool_t AliAnalysisTaskUniFlow::HasTrackPIDTPC(const AliAODTrack* track) const
{
  // Checks if the track has ok PID information from TPC
  // *************************************************************
  if(!track || !fPIDResponse) return kFALSE;
  AliPIDResponse::EDetPidStatus pidStatusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, track);
  return (pidStatusTPC == AliPIDResponse::kDetPidOk);
}
// ============================================================================
Bool_t AliAnalysisTaskUniFlow::HasTrackPIDTOF(const AliAODTrack* track) const
{
  // Checks if the track has ok PID information from TOF
  // *************************************************************
  if(!track || !fPIDResponse) return kFALSE;
  AliPIDResponse::EDetPidStatus pidStatusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, track);
  return ((pidStatusTOF == AliPIDResponse::kDetPidOk) && (track->GetStatus()& AliVTrack::kTOFout) && (track->GetStatus()& AliVTrack::kTIME));
}
// ============================================================================
void AliAnalysisTaskUniFlow::AddCorr(std::vector<Int_t> harms, std::vector<Double_t> gaps, Bool_t doRFPs, Bool_t doPOIs)
{
    fVecCorrTask.push_back(new AliUniFlowCorrTask(doRFPs, doPOIs, harms, gaps));
}
// ============================================================================
void AliAnalysisTaskUniFlow::Terminate(Option_t* option)
{
  // called on end of task, after all events are processed
  // *************************************************************
  AliAnalysisTaskSE::Terminate(option);
  return;
}
// ============================================================================
// Set of methods returning given complex flow vector based on flow harmonics (n) and weight power indexes (p)
// a la General Framework implementation.
// Q: flow vector of RFPs (with/out eta gap)
// P: flow vector of POIs (with/out eta gap) (in usual notation p)
// S: flow vector of overlaping RFPs and POIs (in usual notation q)

TComplex AliAnalysisTaskUniFlow::Q(const Int_t n, const Int_t p) const
{
  if (n < 0) return TComplex::Conjugate(fFlowVecQpos[-n][p]);
  else return fFlowVecQpos[n][p];
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::QGapPos(const Int_t n, const Int_t p) const
{
  if (n < 0) return TComplex::Conjugate(fFlowVecQpos[-n][p]);
  else return fFlowVecQpos[n][p];
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::QGapNeg(const Int_t n, const Int_t p) const
{
  if(n < 0) return TComplex::Conjugate(fFlowVecQneg[-n][p]);
  else return fFlowVecQneg[n][p];
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::QGapMid(const Int_t n, const Int_t p) const
{
  if(n < 0) return TComplex::Conjugate(fFlowVecQmid[-n][p]);
  else return fFlowVecQmid[n][p];
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::P(const Int_t n, const Int_t p) const
{
  if(n < 0) return TComplex::Conjugate(fFlowVecPpos[-n][p]);
  else return fFlowVecPpos[n][p];
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::PGapPos(const Int_t n, const Int_t p) const
{
  if(n < 0) return TComplex::Conjugate(fFlowVecPpos[-n][p]);
  else return fFlowVecPpos[n][p];
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::PGapNeg(const Int_t n, const Int_t p) const
{
  if(n < 0) return TComplex::Conjugate(fFlowVecPneg[-n][p]);
  else return fFlowVecPneg[n][p];
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::S(const Int_t n, const Int_t p) const
{
  if(n < 0) return TComplex::Conjugate(fFlowVecSpos[-n][p]);
  else return fFlowVecSpos[n][p];
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::SGapPos(const Int_t n, const Int_t p) const
{
  if(n < 0) return TComplex::Conjugate(fFlowVecSpos[-n][p]);
  else return fFlowVecSpos[n][p];
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::SGapNeg(const Int_t n, const Int_t p) const
{
  if(n < 0) return TComplex::Conjugate(fFlowVecSneg[-n][p]);
  else return fFlowVecSneg[n][p];
}
// ============================================================================

// Set of flow calculation methods for cumulants of different orders with/out eta gap

// ============================================================================
TComplex AliAnalysisTaskUniFlow::Two(const Int_t n1, const Int_t n2) const
{
  TComplex formula = Q(n1,1)*Q(n2,1) - Q(n1+n2,2);
  return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::TwoGap(const Int_t n1, const Int_t n2) const
{
  TComplex formula = QGapPos(n1,1)*QGapNeg(n2,1);
  return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::TwoDiff(const Int_t n1, const Int_t n2) const
{
  TComplex formula = P(n1,1)*Q(n2,1) - S(n1+n2,2);
  return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::TwoDiffGapPos(const Int_t n1, const Int_t n2) const
{
  TComplex formula = PGapPos(n1,1)*QGapNeg(n2,1);
  return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::TwoDiffGapNeg(const Int_t n1, const Int_t n2) const
{
  TComplex formula = PGapNeg(n1,1)*QGapPos(n2,1);
  return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::TwoPos(const Int_t n1, const Int_t n2) const
{
  TComplex formula = QGapPos(n1,1)*QGapPos(n2,1) - QGapPos(n1+n2,2);
  return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::TwoNeg(const Int_t n1, const Int_t n2) const
{
  TComplex formula = QGapNeg(n1,1)*QGapNeg(n2,1) - QGapNeg(n1+n2,2);
  return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::TwoDiffPos(const Int_t n1, const Int_t n2) const
{
  TComplex formula = PGapPos(n1,1)*QGapPos(n2,1) - SGapPos(n1+n2,2);
  return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::TwoDiffNeg(const Int_t n1, const Int_t n2) const
{
  TComplex formula = PGapNeg(n1,1)*QGapNeg(n2,1) - SGapNeg(n1+n2,2);
  return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::Three(const Int_t n1, const Int_t n2, const Int_t n3) const
{
  TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)-Q(n1+n2,2)*Q(n3,1)-Q(n2,1)*Q(n1+n3,2)
 		                 - Q(n1,1)*Q(n2+n3,2)+2.*Q(n1+n2+n3,3);
  return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::ThreePos(const Int_t n1, const Int_t n2, const Int_t n3) const
{
  TComplex formula = QGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n3,1)-QGapPos(n1+n2,2)*QGapPos(n3,1)-QGapPos(n2,1)*QGapPos(n1+n3,2)
 		                 - QGapPos(n1,1)*QGapPos(n2+n3,2)+2.*QGapPos(n1+n2+n3,3);
  return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::ThreeNeg(const Int_t n1, const Int_t n2, const Int_t n3) const
{
  TComplex formula = QGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n3,1)-QGapNeg(n1+n2,2)*QGapNeg(n3,1)-QGapNeg(n2,1)*QGapNeg(n1+n3,2)
 		                 - QGapNeg(n1,1)*QGapNeg(n2+n3,2)+2.*QGapNeg(n1+n2+n3,3);
  return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::Four(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4) const
{
  TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*Q(n1+n3,2)*Q(n4,1)
                    - Q(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.0*Q(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)
                    + Q(n2+n3,2)*Q(n1+n4,2)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)+Q(n1+n3,2)*Q(n2+n4,2)
                    + 2.0*Q(n3,1)*Q(n1+n2+n4,3)-Q(n1,1)*Q(n2,1)*Q(n3+n4,2)+Q(n1+n2,2)*Q(n3+n4,2)
                    + 2.0*Q(n2,1)*Q(n1+n3+n4,3)+2.0*Q(n1,1)*Q(n2+n3+n4,3)-6.0*Q(n1+n2+n3+n4,4);
  return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::FourPos(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4) const
{
  TComplex formula = QGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n3,1)*QGapPos(n4,1)-QGapPos(n1+n2,2)*QGapPos(n3,1)*QGapPos(n4,1)-QGapPos(n2,1)*QGapPos(n1+n3,2)*QGapPos(n4,1)
                    - QGapPos(n1,1)*QGapPos(n2+n3,2)*QGapPos(n4,1)+2.0*QGapPos(n1+n2+n3,3)*QGapPos(n4,1)-QGapPos(n2,1)*QGapPos(n3,1)*QGapPos(n1+n4,2)
                    + QGapPos(n2+n3,2)*QGapPos(n1+n4,2)-QGapPos(n1,1)*QGapPos(n3,1)*QGapPos(n2+n4,2)+QGapPos(n1+n3,2)*QGapPos(n2+n4,2)
                    + 2.0*QGapPos(n3,1)*QGapPos(n1+n2+n4,3)-QGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n3+n4,2)+QGapPos(n1+n2,2)*QGapPos(n3+n4,2)
                    + 2.0*QGapPos(n2,1)*QGapPos(n1+n3+n4,3)+2.0*QGapPos(n1,1)*QGapPos(n2+n3+n4,3)-6.0*QGapPos(n1+n2+n3+n4,4);
  return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::FourNeg(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4) const
{
  TComplex formula = QGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n3,1)*QGapNeg(n4,1)-QGapNeg(n1+n2,2)*QGapNeg(n3,1)*QGapNeg(n4,1)-QGapNeg(n2,1)*QGapNeg(n1+n3,2)*QGapNeg(n4,1)
                    - QGapNeg(n1,1)*QGapNeg(n2+n3,2)*QGapNeg(n4,1)+2.0*QGapNeg(n1+n2+n3,3)*QGapNeg(n4,1)-QGapNeg(n2,1)*QGapNeg(n3,1)*QGapNeg(n1+n4,2)
                    + QGapNeg(n2+n3,2)*QGapNeg(n1+n4,2)-QGapNeg(n1,1)*QGapNeg(n3,1)*QGapNeg(n2+n4,2)+QGapNeg(n1+n3,2)*QGapNeg(n2+n4,2)
                    + 2.0*QGapNeg(n3,1)*QGapNeg(n1+n2+n4,3)-QGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n3+n4,2)+QGapNeg(n1+n2,2)*QGapNeg(n3+n4,2)
                    + 2.0*QGapNeg(n2,1)*QGapNeg(n1+n3+n4,3)+2.0*QGapNeg(n1,1)*QGapNeg(n2+n3+n4,3)-6.0*QGapNeg(n1+n2+n3+n4,4);
  return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::Five(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t n5) const
{
    TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)
    - Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)-Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)
    + 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)
    + Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)
    + Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)+2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)
    - Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)+Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)
    + 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)+2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)
    - 6.*Q(n1+n2+n3+n4,4)*Q(n5,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)
    + Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)+Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)
    + Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)-2.*Q(n2+n3+n4,3)*Q(n1+n5,2)
    - Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)+Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)
    + Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)+Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)
    - 2.*Q(n1+n3+n4,3)*Q(n2+n5,2)+2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)
    - 2.*Q(n3+n4,2)*Q(n1+n2+n5,3)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)
    + Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)+Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)
    + Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)-2.*Q(n1+n2+n4,3)*Q(n3+n5,2)
    + 2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)-2.*Q(n2+n4,2)*Q(n1+n3+n5,3)
    + 2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)-2.*Q(n1+n4,2)*Q(n2+n3+n5,3)
    - 6.*Q(n4,1)*Q(n1+n2+n3+n5,4)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)
    + Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)+Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)
    + Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)-2.*Q(n1+n2+n3,3)*Q(n4+n5,2)
    + 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)-2.*Q(n2+n3,2)*Q(n1+n4+n5,3)
    + 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)-2.*Q(n1+n3,2)*Q(n2+n4+n5,3)
    - 6.*Q(n3,1)*Q(n1+n2+n4+n5,4)+2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3)
    - 2.*Q(n1+n2,2)*Q(n3+n4+n5,3)-6.*Q(n2,1)*Q(n1+n3+n4+n5,4)
    - 6.*Q(n1,1)*Q(n2+n3+n4+n5,4)+24.*Q(n1+n2+n3+n4+n5,5);
    return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::FivePos(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t n5) const
{
    TComplex formula = QGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n3,1)*QGapPos(n4,1)*QGapPos(n5,1)-QGapPos(n1+n2,2)*QGapPos(n3,1)*QGapPos(n4,1)*QGapPos(n5,1)
    - QGapPos(n2,1)*QGapPos(n1+n3,2)*QGapPos(n4,1)*QGapPos(n5,1)-QGapPos(n1,1)*QGapPos(n2+n3,2)*QGapPos(n4,1)*QGapPos(n5,1)
    + 2.*QGapPos(n1+n2+n3,3)*QGapPos(n4,1)*QGapPos(n5,1)-QGapPos(n2,1)*QGapPos(n3,1)*QGapPos(n1+n4,2)*QGapPos(n5,1)
    + QGapPos(n2+n3,2)*QGapPos(n1+n4,2)*QGapPos(n5,1)-QGapPos(n1,1)*QGapPos(n3,1)*QGapPos(n2+n4,2)*QGapPos(n5,1)
    + QGapPos(n1+n3,2)*QGapPos(n2+n4,2)*QGapPos(n5,1)+2.*QGapPos(n3,1)*QGapPos(n1+n2+n4,3)*QGapPos(n5,1)
    - QGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n3+n4,2)*QGapPos(n5,1)+QGapPos(n1+n2,2)*QGapPos(n3+n4,2)*QGapPos(n5,1)
    + 2.*QGapPos(n2,1)*QGapPos(n1+n3+n4,3)*QGapPos(n5,1)+2.*QGapPos(n1,1)*QGapPos(n2+n3+n4,3)*QGapPos(n5,1)
    - 6.*QGapPos(n1+n2+n3+n4,4)*QGapPos(n5,1)-QGapPos(n2,1)*QGapPos(n3,1)*QGapPos(n4,1)*QGapPos(n1+n5,2)
    + QGapPos(n2+n3,2)*QGapPos(n4,1)*QGapPos(n1+n5,2)+QGapPos(n3,1)*QGapPos(n2+n4,2)*QGapPos(n1+n5,2)
    + QGapPos(n2,1)*QGapPos(n3+n4,2)*QGapPos(n1+n5,2)-2.*QGapPos(n2+n3+n4,3)*QGapPos(n1+n5,2)
    - QGapPos(n1,1)*QGapPos(n3,1)*QGapPos(n4,1)*QGapPos(n2+n5,2)+QGapPos(n1+n3,2)*QGapPos(n4,1)*QGapPos(n2+n5,2)
    + QGapPos(n3,1)*QGapPos(n1+n4,2)*QGapPos(n2+n5,2)+QGapPos(n1,1)*QGapPos(n3+n4,2)*QGapPos(n2+n5,2)
    - 2.*QGapPos(n1+n3+n4,3)*QGapPos(n2+n5,2)+2.*QGapPos(n3,1)*QGapPos(n4,1)*QGapPos(n1+n2+n5,3)
    - 2.*QGapPos(n3+n4,2)*QGapPos(n1+n2+n5,3)-QGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n4,1)*QGapPos(n3+n5,2)
    + QGapPos(n1+n2,2)*QGapPos(n4,1)*QGapPos(n3+n5,2)+QGapPos(n2,1)*QGapPos(n1+n4,2)*QGapPos(n3+n5,2)
    + QGapPos(n1,1)*QGapPos(n2+n4,2)*QGapPos(n3+n5,2)-2.*QGapPos(n1+n2+n4,3)*QGapPos(n3+n5,2)
    + 2.*QGapPos(n2,1)*QGapPos(n4,1)*QGapPos(n1+n3+n5,3)-2.*QGapPos(n2+n4,2)*QGapPos(n1+n3+n5,3)
    + 2.*QGapPos(n1,1)*QGapPos(n4,1)*QGapPos(n2+n3+n5,3)-2.*QGapPos(n1+n4,2)*QGapPos(n2+n3+n5,3)
    - 6.*QGapPos(n4,1)*QGapPos(n1+n2+n3+n5,4)-QGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n3,1)*QGapPos(n4+n5,2)
    + QGapPos(n1+n2,2)*QGapPos(n3,1)*QGapPos(n4+n5,2)+QGapPos(n2,1)*QGapPos(n1+n3,2)*QGapPos(n4+n5,2)
    + QGapPos(n1,1)*QGapPos(n2+n3,2)*QGapPos(n4+n5,2)-2.*QGapPos(n1+n2+n3,3)*QGapPos(n4+n5,2)
    + 2.*QGapPos(n2,1)*QGapPos(n3,1)*QGapPos(n1+n4+n5,3)-2.*QGapPos(n2+n3,2)*QGapPos(n1+n4+n5,3)
    + 2.*QGapPos(n1,1)*QGapPos(n3,1)*QGapPos(n2+n4+n5,3)-2.*QGapPos(n1+n3,2)*QGapPos(n2+n4+n5,3)
    - 6.*QGapPos(n3,1)*QGapPos(n1+n2+n4+n5,4)+2.*QGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n3+n4+n5,3)
    - 2.*QGapPos(n1+n2,2)*QGapPos(n3+n4+n5,3)-6.*QGapPos(n2,1)*QGapPos(n1+n3+n4+n5,4)
    - 6.*QGapPos(n1,1)*QGapPos(n2+n3+n4+n5,4)+24.*QGapPos(n1+n2+n3+n4+n5,5);
    return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::FiveNeg(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t n5) const
{
    TComplex formula = QGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n3,1)*QGapNeg(n4,1)*QGapNeg(n5,1)-QGapNeg(n1+n2,2)*QGapNeg(n3,1)*QGapNeg(n4,1)*QGapNeg(n5,1)
    - QGapNeg(n2,1)*QGapNeg(n1+n3,2)*QGapNeg(n4,1)*QGapNeg(n5,1)-QGapNeg(n1,1)*QGapNeg(n2+n3,2)*QGapNeg(n4,1)*QGapNeg(n5,1)
    + 2.*QGapNeg(n1+n2+n3,3)*QGapNeg(n4,1)*QGapNeg(n5,1)-QGapNeg(n2,1)*QGapNeg(n3,1)*QGapNeg(n1+n4,2)*QGapNeg(n5,1)
    + QGapNeg(n2+n3,2)*QGapNeg(n1+n4,2)*QGapNeg(n5,1)-QGapNeg(n1,1)*QGapNeg(n3,1)*QGapNeg(n2+n4,2)*QGapNeg(n5,1)
    + QGapNeg(n1+n3,2)*QGapNeg(n2+n4,2)*QGapNeg(n5,1)+2.*QGapNeg(n3,1)*QGapNeg(n1+n2+n4,3)*QGapNeg(n5,1)
    - QGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n3+n4,2)*QGapNeg(n5,1)+QGapNeg(n1+n2,2)*QGapNeg(n3+n4,2)*QGapNeg(n5,1)
    + 2.*QGapNeg(n2,1)*QGapNeg(n1+n3+n4,3)*QGapNeg(n5,1)+2.*QGapNeg(n1,1)*QGapNeg(n2+n3+n4,3)*QGapNeg(n5,1)
    - 6.*QGapNeg(n1+n2+n3+n4,4)*QGapNeg(n5,1)-QGapNeg(n2,1)*QGapNeg(n3,1)*QGapNeg(n4,1)*QGapNeg(n1+n5,2)
    + QGapNeg(n2+n3,2)*QGapNeg(n4,1)*QGapNeg(n1+n5,2)+QGapNeg(n3,1)*QGapNeg(n2+n4,2)*QGapNeg(n1+n5,2)
    + QGapNeg(n2,1)*QGapNeg(n3+n4,2)*QGapNeg(n1+n5,2)-2.*QGapNeg(n2+n3+n4,3)*QGapNeg(n1+n5,2)
    - QGapNeg(n1,1)*QGapNeg(n3,1)*QGapNeg(n4,1)*QGapNeg(n2+n5,2)+QGapNeg(n1+n3,2)*QGapNeg(n4,1)*QGapNeg(n2+n5,2)
    + QGapNeg(n3,1)*QGapNeg(n1+n4,2)*QGapNeg(n2+n5,2)+QGapNeg(n1,1)*QGapNeg(n3+n4,2)*QGapNeg(n2+n5,2)
    - 2.*QGapNeg(n1+n3+n4,3)*QGapNeg(n2+n5,2)+2.*QGapNeg(n3,1)*QGapNeg(n4,1)*QGapNeg(n1+n2+n5,3)
    - 2.*QGapNeg(n3+n4,2)*QGapNeg(n1+n2+n5,3)-QGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n4,1)*QGapNeg(n3+n5,2)
    + QGapNeg(n1+n2,2)*QGapNeg(n4,1)*QGapNeg(n3+n5,2)+QGapNeg(n2,1)*QGapNeg(n1+n4,2)*QGapNeg(n3+n5,2)
    + QGapNeg(n1,1)*QGapNeg(n2+n4,2)*QGapNeg(n3+n5,2)-2.*QGapNeg(n1+n2+n4,3)*QGapNeg(n3+n5,2)
    + 2.*QGapNeg(n2,1)*QGapNeg(n4,1)*QGapNeg(n1+n3+n5,3)-2.*QGapNeg(n2+n4,2)*QGapNeg(n1+n3+n5,3)
    + 2.*QGapNeg(n1,1)*QGapNeg(n4,1)*QGapNeg(n2+n3+n5,3)-2.*QGapNeg(n1+n4,2)*QGapNeg(n2+n3+n5,3)
    - 6.*QGapNeg(n4,1)*QGapNeg(n1+n2+n3+n5,4)-QGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n3,1)*QGapNeg(n4+n5,2)
    + QGapNeg(n1+n2,2)*QGapNeg(n3,1)*QGapNeg(n4+n5,2)+QGapNeg(n2,1)*QGapNeg(n1+n3,2)*QGapNeg(n4+n5,2)
    + QGapNeg(n1,1)*QGapNeg(n2+n3,2)*QGapNeg(n4+n5,2)-2.*QGapNeg(n1+n2+n3,3)*QGapNeg(n4+n5,2)
    + 2.*QGapNeg(n2,1)*QGapNeg(n3,1)*QGapNeg(n1+n4+n5,3)-2.*QGapNeg(n2+n3,2)*QGapNeg(n1+n4+n5,3)
    + 2.*QGapNeg(n1,1)*QGapNeg(n3,1)*QGapNeg(n2+n4+n5,3)-2.*QGapNeg(n1+n3,2)*QGapNeg(n2+n4+n5,3)
    - 6.*QGapNeg(n3,1)*QGapNeg(n1+n2+n4+n5,4)+2.*QGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n3+n4+n5,3)
    - 2.*QGapNeg(n1+n2,2)*QGapNeg(n3+n4+n5,3)-6.*QGapNeg(n2,1)*QGapNeg(n1+n3+n4+n5,4)
    - 6.*QGapNeg(n1,1)*QGapNeg(n2+n3+n4+n5,4)+24.*QGapNeg(n1+n2+n3+n4+n5,5);
    return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::Six(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t n5, const Int_t n6) const
{
  TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)
    - Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)
    + 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)
    + Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)
    + Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)+2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)*Q(n6,1)
    - Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)+Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)
    + 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)*Q(n6,1)+2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)*Q(n6,1)
    - 6.*Q(n1+n2+n3+n4,4)*Q(n5,1)*Q(n6,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)
    + Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)+Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)*Q(n6,1)
    + Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)*Q(n6,1)-2.*Q(n2+n3+n4,3)*Q(n1+n5,2)*Q(n6,1)
    - Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)+Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)
    + Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)*Q(n6,1)+Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)*Q(n6,1)
    - 2.*Q(n1+n3+n4,3)*Q(n2+n5,2)*Q(n6,1)+2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)*Q(n6,1)
    - 2.*Q(n3+n4,2)*Q(n1+n2+n5,3)*Q(n6,1)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)
    + Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)+Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)*Q(n6,1)
    + Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)*Q(n6,1)-2.*Q(n1+n2+n4,3)*Q(n3+n5,2)*Q(n6,1)
    + 2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)*Q(n6,1)-2.*Q(n2+n4,2)*Q(n1+n3+n5,3)*Q(n6,1)
    + 2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)*Q(n6,1)-2.*Q(n1+n4,2)*Q(n2+n3+n5,3)*Q(n6,1)
    - 6.*Q(n4,1)*Q(n1+n2+n3+n5,4)*Q(n6,1)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)
    + Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)+Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)*Q(n6,1)
    + Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)*Q(n6,1)-2.*Q(n1+n2+n3,3)*Q(n4+n5,2)*Q(n6,1)
    + 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)*Q(n6,1)-2.*Q(n2+n3,2)*Q(n1+n4+n5,3)*Q(n6,1)
    + 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)*Q(n6,1)-2.*Q(n1+n3,2)*Q(n2+n4+n5,3)*Q(n6,1)
    - 6.*Q(n3,1)*Q(n1+n2+n4+n5,4)*Q(n6,1)+2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3)*Q(n6,1)
    - 2.*Q(n1+n2,2)*Q(n3+n4+n5,3)*Q(n6,1)-6.*Q(n2,1)*Q(n1+n3+n4+n5,4)*Q(n6,1)
    - 6.*Q(n1,1)*Q(n2+n3+n4+n5,4)*Q(n6,1)+24.*Q(n1+n2+n3+n4+n5,5)*Q(n6,1)
    - Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)+Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)
    + Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n1+n6,2)+Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n1+n6,2)
    - 2.*Q(n2+n3+n4,3)*Q(n5,1)*Q(n1+n6,2)+Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n1+n6,2)
    - Q(n3+n4,2)*Q(n2+n5,2)*Q(n1+n6,2)+Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n1+n6,2)
    - Q(n2+n4,2)*Q(n3+n5,2)*Q(n1+n6,2)-2.*Q(n4,1)*Q(n2+n3+n5,3)*Q(n1+n6,2)
    + Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n1+n6,2)-Q(n2+n3,2)*Q(n4+n5,2)*Q(n1+n6,2)
    - 2.*Q(n3,1)*Q(n2+n4+n5,3)*Q(n1+n6,2)-2.*Q(n2,1)*Q(n3+n4+n5,3)*Q(n1+n6,2)
    + 6.*Q(n2+n3+n4+n5,4)*Q(n1+n6,2)-Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)
    + Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)+Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n2+n6,2)
    + Q(n1,1)*Q(n3+n4,2)*Q(n5,1)*Q(n2+n6,2)-2.*Q(n1+n3+n4,3)*Q(n5,1)*Q(n2+n6,2)
    + Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n2+n6,2)-Q(n3+n4,2)*Q(n1+n5,2)*Q(n2+n6,2)
    + Q(n1,1)*Q(n4,1)*Q(n3+n5,2)*Q(n2+n6,2)-Q(n1+n4,2)*Q(n3+n5,2)*Q(n2+n6,2)
    - 2.*Q(n4,1)*Q(n1+n3+n5,3)*Q(n2+n6,2)+Q(n1,1)*Q(n3,1)*Q(n4+n5,2)*Q(n2+n6,2)
    - Q(n1+n3,2)*Q(n4+n5,2)*Q(n2+n6,2)-2.*Q(n3,1)*Q(n1+n4+n5,3)*Q(n2+n6,2)
    - 2.*Q(n1,1)*Q(n3+n4+n5,3)*Q(n2+n6,2)+6.*Q(n1+n3+n4+n5,4)*Q(n2+n6,2)
    + 2.*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n2+n6,3)-2.*Q(n3+n4,2)*Q(n5,1)*Q(n1+n2+n6,3)
    - 2.*Q(n4,1)*Q(n3+n5,2)*Q(n1+n2+n6,3)-2.*Q(n3,1)*Q(n4+n5,2)*Q(n1+n2+n6,3)
    + 4.*Q(n3+n4+n5,3)*Q(n1+n2+n6,3)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)
    + Q(n1+n2,2)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)+Q(n2,1)*Q(n1+n4,2)*Q(n5,1)*Q(n3+n6,2)
    + Q(n1,1)*Q(n2+n4,2)*Q(n5,1)*Q(n3+n6,2)-2.*Q(n1+n2+n4,3)*Q(n5,1)*Q(n3+n6,2)
    + Q(n2,1)*Q(n4,1)*Q(n1+n5,2)*Q(n3+n6,2)-Q(n2+n4,2)*Q(n1+n5,2)*Q(n3+n6,2)
    + Q(n1,1)*Q(n4,1)*Q(n2+n5,2)*Q(n3+n6,2)-Q(n1+n4,2)*Q(n2+n5,2)*Q(n3+n6,2)
    - 2.*Q(n4,1)*Q(n1+n2+n5,3)*Q(n3+n6,2)+Q(n1,1)*Q(n2,1)*Q(n4+n5,2)*Q(n3+n6,2)
    - Q(n1+n2,2)*Q(n4+n5,2)*Q(n3+n6,2)-2.*Q(n2,1)*Q(n1+n4+n5,3)*Q(n3+n6,2)
    - 2.*Q(n1,1)*Q(n2+n4+n5,3)*Q(n3+n6,2)+6.*Q(n1+n2+n4+n5,4)*Q(n3+n6,2)
    + 2.*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n1+n3+n6,3)-2.*Q(n2+n4,2)*Q(n5,1)*Q(n1+n3+n6,3)
    - 2.*Q(n4,1)*Q(n2+n5,2)*Q(n1+n3+n6,3)-2.*Q(n2,1)*Q(n4+n5,2)*Q(n1+n3+n6,3)
    + 4.*Q(n2+n4+n5,3)*Q(n1+n3+n6,3)+2.*Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n2+n3+n6,3)
    - 2.*Q(n1+n4,2)*Q(n5,1)*Q(n2+n3+n6,3)-2.*Q(n4,1)*Q(n1+n5,2)*Q(n2+n3+n6,3)
    - 2.*Q(n1,1)*Q(n4+n5,2)*Q(n2+n3+n6,3)+4.*Q(n1+n4+n5,3)*Q(n2+n3+n6,3)
    - 6.*Q(n4,1)*Q(n5,1)*Q(n1+n2+n3+n6,4)+6.*Q(n4+n5,2)*Q(n1+n2+n3+n6,4)
    - Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)+Q(n1+n2,2)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)
    + Q(n2,1)*Q(n1+n3,2)*Q(n5,1)*Q(n4+n6,2)+Q(n1,1)*Q(n2+n3,2)*Q(n5,1)*Q(n4+n6,2)
    - 2.*Q(n1+n2+n3,3)*Q(n5,1)*Q(n4+n6,2)+Q(n2,1)*Q(n3,1)*Q(n1+n5,2)*Q(n4+n6,2)
    - Q(n2+n3,2)*Q(n1+n5,2)*Q(n4+n6,2)+Q(n1,1)*Q(n3,1)*Q(n2+n5,2)*Q(n4+n6,2)
    - Q(n1+n3,2)*Q(n2+n5,2)*Q(n4+n6,2)-2.*Q(n3,1)*Q(n1+n2+n5,3)*Q(n4+n6,2)
    + Q(n1,1)*Q(n2,1)*Q(n3+n5,2)*Q(n4+n6,2)-Q(n1+n2,2)*Q(n3+n5,2)*Q(n4+n6,2)
    - 2.*Q(n2,1)*Q(n1+n3+n5,3)*Q(n4+n6,2)-2.*Q(n1,1)*Q(n2+n3+n5,3)*Q(n4+n6,2)
    + 6.*Q(n1+n2+n3+n5,4)*Q(n4+n6,2)+2.*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n1+n4+n6,3)
    - 2.*Q(n2+n3,2)*Q(n5,1)*Q(n1+n4+n6,3)-2.*Q(n3,1)*Q(n2+n5,2)*Q(n1+n4+n6,3)
    - 2.*Q(n2,1)*Q(n3+n5,2)*Q(n1+n4+n6,3)+4.*Q(n2+n3+n5,3)*Q(n1+n4+n6,3)
    + 2.*Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n2+n4+n6,3)-2.*Q(n1+n3,2)*Q(n5,1)*Q(n2+n4+n6,3)
    - 2.*Q(n3,1)*Q(n1+n5,2)*Q(n2+n4+n6,3)-2.*Q(n1,1)*Q(n3+n5,2)*Q(n2+n4+n6,3)
    + 4.*Q(n1+n3+n5,3)*Q(n2+n4+n6,3)-6.*Q(n3,1)*Q(n5,1)*Q(n1+n2+n4+n6,4)
    + 6.*Q(n3+n5,2)*Q(n1+n2+n4+n6,4)+2.*Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n3+n4+n6,3)
    - 2.*Q(n1+n2,2)*Q(n5,1)*Q(n3+n4+n6,3)-2.*Q(n2,1)*Q(n1+n5,2)*Q(n3+n4+n6,3)
    - 2.*Q(n1,1)*Q(n2+n5,2)*Q(n3+n4+n6,3)+4.*Q(n1+n2+n5,3)*Q(n3+n4+n6,3)
    - 6.*Q(n2,1)*Q(n5,1)*Q(n1+n3+n4+n6,4)+6.*Q(n2+n5,2)*Q(n1+n3+n4+n6,4)
    - 6.*Q(n1,1)*Q(n5,1)*Q(n2+n3+n4+n6,4)+6.*Q(n1+n5,2)*Q(n2+n3+n4+n6,4)
    + 24.*Q(n5,1)*Q(n1+n2+n3+n4+n6,5)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)
    + Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)+Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5+n6,2)
    + Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5+n6,2)-2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5+n6,2)
    + Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5+n6,2)-Q(n2+n3,2)*Q(n1+n4,2)*Q(n5+n6,2)
    + Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5+n6,2)-Q(n1+n3,2)*Q(n2+n4,2)*Q(n5+n6,2)
    - 2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5+n6,2)+Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5+n6,2)
    - Q(n1+n2,2)*Q(n3+n4,2)*Q(n5+n6,2)-2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5+n6,2)
    - 2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5+n6,2)+6.*Q(n1+n2+n3+n4,4)*Q(n5+n6,2)
    + 2.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5+n6,3)-2.*Q(n2+n3,2)*Q(n4,1)*Q(n1+n5+n6,3)
    - 2.*Q(n3,1)*Q(n2+n4,2)*Q(n1+n5+n6,3)-2.*Q(n2,1)*Q(n3+n4,2)*Q(n1+n5+n6,3)
    + 4.*Q(n2+n3+n4,3)*Q(n1+n5+n6,3)+2.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5+n6,3)
    - 2.*Q(n1+n3,2)*Q(n4,1)*Q(n2+n5+n6,3)-2.*Q(n3,1)*Q(n1+n4,2)*Q(n2+n5+n6,3)
    - 2.*Q(n1,1)*Q(n3+n4,2)*Q(n2+n5+n6,3)+4.*Q(n1+n3+n4,3)*Q(n2+n5+n6,3)
    - 6.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5+n6,4)+6.*Q(n3+n4,2)*Q(n1+n2+n5+n6,4)
    + 2.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5+n6,3)-2.*Q(n1+n2,2)*Q(n4,1)*Q(n3+n5+n6,3)
    - 2.*Q(n2,1)*Q(n1+n4,2)*Q(n3+n5+n6,3)-2.*Q(n1,1)*Q(n2+n4,2)*Q(n3+n5+n6,3)
    + 4.*Q(n1+n2+n4,3)*Q(n3+n5+n6,3)-6.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5+n6,4)
    + 6.*Q(n2+n4,2)*Q(n1+n3+n5+n6,4)-6.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5+n6,4)
    + 6.*Q(n1+n4,2)*Q(n2+n3+n5+n6,4)+24.*Q(n4,1)*Q(n1+n2+n3+n5+n6,5)
    + 2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5+n6,3)-2.*Q(n1+n2,2)*Q(n3,1)*Q(n4+n5+n6,3)
    - 2.*Q(n2,1)*Q(n1+n3,2)*Q(n4+n5+n6,3)-2.*Q(n1,1)*Q(n2+n3,2)*Q(n4+n5+n6,3)
    + 4.*Q(n1+n2+n3,3)*Q(n4+n5+n6,3)-6.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5+n6,4)
    + 6.*Q(n2+n3,2)*Q(n1+n4+n5+n6,4)-6.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5+n6,4)
    + 6.*Q(n1+n3,2)*Q(n2+n4+n5+n6,4)+24.*Q(n3,1)*Q(n1+n2+n4+n5+n6,5)
    - 6.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5+n6,4)+6.*Q(n1+n2,2)*Q(n3+n4+n5+n6,4)
    + 24.*Q(n2,1)*Q(n1+n3+n4+n5+n6,5)+24.*Q(n1,1)*Q(n2+n3+n4+n5+n6,5)
    - 120.*Q(n1+n2+n3+n4+n5+n6,6);
    return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::SixPos(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t n5, const Int_t n6) const
{
  TComplex formula = QGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n3,1)*QGapPos(n4,1)*QGapPos(n5,1)*QGapPos(n6,1)-QGapPos(n1+n2,2)*QGapPos(n3,1)*QGapPos(n4,1)*QGapPos(n5,1)*QGapPos(n6,1)
    - QGapPos(n2,1)*QGapPos(n1+n3,2)*QGapPos(n4,1)*QGapPos(n5,1)*QGapPos(n6,1)-QGapPos(n1,1)*QGapPos(n2+n3,2)*QGapPos(n4,1)*QGapPos(n5,1)*QGapPos(n6,1)
    + 2.*QGapPos(n1+n2+n3,3)*QGapPos(n4,1)*QGapPos(n5,1)*QGapPos(n6,1)-QGapPos(n2,1)*QGapPos(n3,1)*QGapPos(n1+n4,2)*QGapPos(n5,1)*QGapPos(n6,1)
    + QGapPos(n2+n3,2)*QGapPos(n1+n4,2)*QGapPos(n5,1)*QGapPos(n6,1)-QGapPos(n1,1)*QGapPos(n3,1)*QGapPos(n2+n4,2)*QGapPos(n5,1)*QGapPos(n6,1)
    + QGapPos(n1+n3,2)*QGapPos(n2+n4,2)*QGapPos(n5,1)*QGapPos(n6,1)+2.*QGapPos(n3,1)*QGapPos(n1+n2+n4,3)*QGapPos(n5,1)*QGapPos(n6,1)
    - QGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n3+n4,2)*QGapPos(n5,1)*QGapPos(n6,1)+QGapPos(n1+n2,2)*QGapPos(n3+n4,2)*QGapPos(n5,1)*QGapPos(n6,1)
    + 2.*QGapPos(n2,1)*QGapPos(n1+n3+n4,3)*QGapPos(n5,1)*QGapPos(n6,1)+2.*QGapPos(n1,1)*QGapPos(n2+n3+n4,3)*QGapPos(n5,1)*QGapPos(n6,1)
    - 6.*QGapPos(n1+n2+n3+n4,4)*QGapPos(n5,1)*QGapPos(n6,1)-QGapPos(n2,1)*QGapPos(n3,1)*QGapPos(n4,1)*QGapPos(n1+n5,2)*QGapPos(n6,1)
    + QGapPos(n2+n3,2)*QGapPos(n4,1)*QGapPos(n1+n5,2)*QGapPos(n6,1)+QGapPos(n3,1)*QGapPos(n2+n4,2)*QGapPos(n1+n5,2)*QGapPos(n6,1)
    + QGapPos(n2,1)*QGapPos(n3+n4,2)*QGapPos(n1+n5,2)*QGapPos(n6,1)-2.*QGapPos(n2+n3+n4,3)*QGapPos(n1+n5,2)*QGapPos(n6,1)
    - QGapPos(n1,1)*QGapPos(n3,1)*QGapPos(n4,1)*QGapPos(n2+n5,2)*QGapPos(n6,1)+QGapPos(n1+n3,2)*QGapPos(n4,1)*QGapPos(n2+n5,2)*QGapPos(n6,1)
    + QGapPos(n3,1)*QGapPos(n1+n4,2)*QGapPos(n2+n5,2)*QGapPos(n6,1)+QGapPos(n1,1)*QGapPos(n3+n4,2)*QGapPos(n2+n5,2)*QGapPos(n6,1)
    - 2.*QGapPos(n1+n3+n4,3)*QGapPos(n2+n5,2)*QGapPos(n6,1)+2.*QGapPos(n3,1)*QGapPos(n4,1)*QGapPos(n1+n2+n5,3)*QGapPos(n6,1)
    - 2.*QGapPos(n3+n4,2)*QGapPos(n1+n2+n5,3)*QGapPos(n6,1)-QGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n4,1)*QGapPos(n3+n5,2)*QGapPos(n6,1)
    + QGapPos(n1+n2,2)*QGapPos(n4,1)*QGapPos(n3+n5,2)*QGapPos(n6,1)+QGapPos(n2,1)*QGapPos(n1+n4,2)*QGapPos(n3+n5,2)*QGapPos(n6,1)
    + QGapPos(n1,1)*QGapPos(n2+n4,2)*QGapPos(n3+n5,2)*QGapPos(n6,1)-2.*QGapPos(n1+n2+n4,3)*QGapPos(n3+n5,2)*QGapPos(n6,1)
    + 2.*QGapPos(n2,1)*QGapPos(n4,1)*QGapPos(n1+n3+n5,3)*QGapPos(n6,1)-2.*QGapPos(n2+n4,2)*QGapPos(n1+n3+n5,3)*QGapPos(n6,1)
    + 2.*QGapPos(n1,1)*QGapPos(n4,1)*QGapPos(n2+n3+n5,3)*QGapPos(n6,1)-2.*QGapPos(n1+n4,2)*QGapPos(n2+n3+n5,3)*QGapPos(n6,1)
    - 6.*QGapPos(n4,1)*QGapPos(n1+n2+n3+n5,4)*QGapPos(n6,1)-QGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n3,1)*QGapPos(n4+n5,2)*QGapPos(n6,1)
    + QGapPos(n1+n2,2)*QGapPos(n3,1)*QGapPos(n4+n5,2)*QGapPos(n6,1)+QGapPos(n2,1)*QGapPos(n1+n3,2)*QGapPos(n4+n5,2)*QGapPos(n6,1)
    + QGapPos(n1,1)*QGapPos(n2+n3,2)*QGapPos(n4+n5,2)*QGapPos(n6,1)-2.*QGapPos(n1+n2+n3,3)*QGapPos(n4+n5,2)*QGapPos(n6,1)
    + 2.*QGapPos(n2,1)*QGapPos(n3,1)*QGapPos(n1+n4+n5,3)*QGapPos(n6,1)-2.*QGapPos(n2+n3,2)*QGapPos(n1+n4+n5,3)*QGapPos(n6,1)
    + 2.*QGapPos(n1,1)*QGapPos(n3,1)*QGapPos(n2+n4+n5,3)*QGapPos(n6,1)-2.*QGapPos(n1+n3,2)*QGapPos(n2+n4+n5,3)*QGapPos(n6,1)
    - 6.*QGapPos(n3,1)*QGapPos(n1+n2+n4+n5,4)*QGapPos(n6,1)+2.*QGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n3+n4+n5,3)*QGapPos(n6,1)
    - 2.*QGapPos(n1+n2,2)*QGapPos(n3+n4+n5,3)*QGapPos(n6,1)-6.*QGapPos(n2,1)*QGapPos(n1+n3+n4+n5,4)*QGapPos(n6,1)
    - 6.*QGapPos(n1,1)*QGapPos(n2+n3+n4+n5,4)*QGapPos(n6,1)+24.*QGapPos(n1+n2+n3+n4+n5,5)*QGapPos(n6,1)
    - QGapPos(n2,1)*QGapPos(n3,1)*QGapPos(n4,1)*QGapPos(n5,1)*QGapPos(n1+n6,2)+QGapPos(n2+n3,2)*QGapPos(n4,1)*QGapPos(n5,1)*QGapPos(n1+n6,2)
    + QGapPos(n3,1)*QGapPos(n2+n4,2)*QGapPos(n5,1)*QGapPos(n1+n6,2)+QGapPos(n2,1)*QGapPos(n3+n4,2)*QGapPos(n5,1)*QGapPos(n1+n6,2)
    - 2.*QGapPos(n2+n3+n4,3)*QGapPos(n5,1)*QGapPos(n1+n6,2)+QGapPos(n3,1)*QGapPos(n4,1)*QGapPos(n2+n5,2)*QGapPos(n1+n6,2)
    - QGapPos(n3+n4,2)*QGapPos(n2+n5,2)*QGapPos(n1+n6,2)+QGapPos(n2,1)*QGapPos(n4,1)*QGapPos(n3+n5,2)*QGapPos(n1+n6,2)
    - QGapPos(n2+n4,2)*QGapPos(n3+n5,2)*QGapPos(n1+n6,2)-2.*QGapPos(n4,1)*QGapPos(n2+n3+n5,3)*QGapPos(n1+n6,2)
    + QGapPos(n2,1)*QGapPos(n3,1)*QGapPos(n4+n5,2)*QGapPos(n1+n6,2)-QGapPos(n2+n3,2)*QGapPos(n4+n5,2)*QGapPos(n1+n6,2)
    - 2.*QGapPos(n3,1)*QGapPos(n2+n4+n5,3)*QGapPos(n1+n6,2)-2.*QGapPos(n2,1)*QGapPos(n3+n4+n5,3)*QGapPos(n1+n6,2)
    + 6.*QGapPos(n2+n3+n4+n5,4)*QGapPos(n1+n6,2)-QGapPos(n1,1)*QGapPos(n3,1)*QGapPos(n4,1)*QGapPos(n5,1)*QGapPos(n2+n6,2)
    + QGapPos(n1+n3,2)*QGapPos(n4,1)*QGapPos(n5,1)*QGapPos(n2+n6,2)+QGapPos(n3,1)*QGapPos(n1+n4,2)*QGapPos(n5,1)*QGapPos(n2+n6,2)
    + QGapPos(n1,1)*QGapPos(n3+n4,2)*QGapPos(n5,1)*QGapPos(n2+n6,2)-2.*QGapPos(n1+n3+n4,3)*QGapPos(n5,1)*QGapPos(n2+n6,2)
    + QGapPos(n3,1)*QGapPos(n4,1)*QGapPos(n1+n5,2)*QGapPos(n2+n6,2)-QGapPos(n3+n4,2)*QGapPos(n1+n5,2)*QGapPos(n2+n6,2)
    + QGapPos(n1,1)*QGapPos(n4,1)*QGapPos(n3+n5,2)*QGapPos(n2+n6,2)-QGapPos(n1+n4,2)*QGapPos(n3+n5,2)*QGapPos(n2+n6,2)
    - 2.*QGapPos(n4,1)*QGapPos(n1+n3+n5,3)*QGapPos(n2+n6,2)+QGapPos(n1,1)*QGapPos(n3,1)*QGapPos(n4+n5,2)*QGapPos(n2+n6,2)
    - QGapPos(n1+n3,2)*QGapPos(n4+n5,2)*QGapPos(n2+n6,2)-2.*QGapPos(n3,1)*QGapPos(n1+n4+n5,3)*QGapPos(n2+n6,2)
    - 2.*QGapPos(n1,1)*QGapPos(n3+n4+n5,3)*QGapPos(n2+n6,2)+6.*QGapPos(n1+n3+n4+n5,4)*QGapPos(n2+n6,2)
    + 2.*QGapPos(n3,1)*QGapPos(n4,1)*QGapPos(n5,1)*QGapPos(n1+n2+n6,3)-2.*QGapPos(n3+n4,2)*QGapPos(n5,1)*QGapPos(n1+n2+n6,3)
    - 2.*QGapPos(n4,1)*QGapPos(n3+n5,2)*QGapPos(n1+n2+n6,3)-2.*QGapPos(n3,1)*QGapPos(n4+n5,2)*QGapPos(n1+n2+n6,3)
    + 4.*QGapPos(n3+n4+n5,3)*QGapPos(n1+n2+n6,3)-QGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n4,1)*QGapPos(n5,1)*QGapPos(n3+n6,2)
    + QGapPos(n1+n2,2)*QGapPos(n4,1)*QGapPos(n5,1)*QGapPos(n3+n6,2)+QGapPos(n2,1)*QGapPos(n1+n4,2)*QGapPos(n5,1)*QGapPos(n3+n6,2)
    + QGapPos(n1,1)*QGapPos(n2+n4,2)*QGapPos(n5,1)*QGapPos(n3+n6,2)-2.*QGapPos(n1+n2+n4,3)*QGapPos(n5,1)*QGapPos(n3+n6,2)
    + QGapPos(n2,1)*QGapPos(n4,1)*QGapPos(n1+n5,2)*QGapPos(n3+n6,2)-QGapPos(n2+n4,2)*QGapPos(n1+n5,2)*QGapPos(n3+n6,2)
    + QGapPos(n1,1)*QGapPos(n4,1)*QGapPos(n2+n5,2)*QGapPos(n3+n6,2)-QGapPos(n1+n4,2)*QGapPos(n2+n5,2)*QGapPos(n3+n6,2)
    - 2.*QGapPos(n4,1)*QGapPos(n1+n2+n5,3)*QGapPos(n3+n6,2)+QGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n4+n5,2)*QGapPos(n3+n6,2)
    - QGapPos(n1+n2,2)*QGapPos(n4+n5,2)*QGapPos(n3+n6,2)-2.*QGapPos(n2,1)*QGapPos(n1+n4+n5,3)*QGapPos(n3+n6,2)
    - 2.*QGapPos(n1,1)*QGapPos(n2+n4+n5,3)*QGapPos(n3+n6,2)+6.*QGapPos(n1+n2+n4+n5,4)*QGapPos(n3+n6,2)
    + 2.*QGapPos(n2,1)*QGapPos(n4,1)*QGapPos(n5,1)*QGapPos(n1+n3+n6,3)-2.*QGapPos(n2+n4,2)*QGapPos(n5,1)*QGapPos(n1+n3+n6,3)
    - 2.*QGapPos(n4,1)*QGapPos(n2+n5,2)*QGapPos(n1+n3+n6,3)-2.*QGapPos(n2,1)*QGapPos(n4+n5,2)*QGapPos(n1+n3+n6,3)
    + 4.*QGapPos(n2+n4+n5,3)*QGapPos(n1+n3+n6,3)+2.*QGapPos(n1,1)*QGapPos(n4,1)*QGapPos(n5,1)*QGapPos(n2+n3+n6,3)
    - 2.*QGapPos(n1+n4,2)*QGapPos(n5,1)*QGapPos(n2+n3+n6,3)-2.*QGapPos(n4,1)*QGapPos(n1+n5,2)*QGapPos(n2+n3+n6,3)
    - 2.*QGapPos(n1,1)*QGapPos(n4+n5,2)*QGapPos(n2+n3+n6,3)+4.*QGapPos(n1+n4+n5,3)*QGapPos(n2+n3+n6,3)
    - 6.*QGapPos(n4,1)*QGapPos(n5,1)*QGapPos(n1+n2+n3+n6,4)+6.*QGapPos(n4+n5,2)*QGapPos(n1+n2+n3+n6,4)
    - QGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n3,1)*QGapPos(n5,1)*QGapPos(n4+n6,2)+QGapPos(n1+n2,2)*QGapPos(n3,1)*QGapPos(n5,1)*QGapPos(n4+n6,2)
    + QGapPos(n2,1)*QGapPos(n1+n3,2)*QGapPos(n5,1)*QGapPos(n4+n6,2)+QGapPos(n1,1)*QGapPos(n2+n3,2)*QGapPos(n5,1)*QGapPos(n4+n6,2)
    - 2.*QGapPos(n1+n2+n3,3)*QGapPos(n5,1)*QGapPos(n4+n6,2)+QGapPos(n2,1)*QGapPos(n3,1)*QGapPos(n1+n5,2)*QGapPos(n4+n6,2)
    - QGapPos(n2+n3,2)*QGapPos(n1+n5,2)*QGapPos(n4+n6,2)+QGapPos(n1,1)*QGapPos(n3,1)*QGapPos(n2+n5,2)*QGapPos(n4+n6,2)
    - QGapPos(n1+n3,2)*QGapPos(n2+n5,2)*QGapPos(n4+n6,2)-2.*QGapPos(n3,1)*QGapPos(n1+n2+n5,3)*QGapPos(n4+n6,2)
    + QGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n3+n5,2)*QGapPos(n4+n6,2)-QGapPos(n1+n2,2)*QGapPos(n3+n5,2)*QGapPos(n4+n6,2)
    - 2.*QGapPos(n2,1)*QGapPos(n1+n3+n5,3)*QGapPos(n4+n6,2)-2.*QGapPos(n1,1)*QGapPos(n2+n3+n5,3)*QGapPos(n4+n6,2)
    + 6.*QGapPos(n1+n2+n3+n5,4)*QGapPos(n4+n6,2)+2.*QGapPos(n2,1)*QGapPos(n3,1)*QGapPos(n5,1)*QGapPos(n1+n4+n6,3)
    - 2.*QGapPos(n2+n3,2)*QGapPos(n5,1)*QGapPos(n1+n4+n6,3)-2.*QGapPos(n3,1)*QGapPos(n2+n5,2)*QGapPos(n1+n4+n6,3)
    - 2.*QGapPos(n2,1)*QGapPos(n3+n5,2)*QGapPos(n1+n4+n6,3)+4.*QGapPos(n2+n3+n5,3)*QGapPos(n1+n4+n6,3)
    + 2.*QGapPos(n1,1)*QGapPos(n3,1)*QGapPos(n5,1)*QGapPos(n2+n4+n6,3)-2.*QGapPos(n1+n3,2)*QGapPos(n5,1)*QGapPos(n2+n4+n6,3)
    - 2.*QGapPos(n3,1)*QGapPos(n1+n5,2)*QGapPos(n2+n4+n6,3)-2.*QGapPos(n1,1)*QGapPos(n3+n5,2)*QGapPos(n2+n4+n6,3)
    + 4.*QGapPos(n1+n3+n5,3)*QGapPos(n2+n4+n6,3)-6.*QGapPos(n3,1)*QGapPos(n5,1)*QGapPos(n1+n2+n4+n6,4)
    + 6.*QGapPos(n3+n5,2)*QGapPos(n1+n2+n4+n6,4)+2.*QGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n5,1)*QGapPos(n3+n4+n6,3)
    - 2.*QGapPos(n1+n2,2)*QGapPos(n5,1)*QGapPos(n3+n4+n6,3)-2.*QGapPos(n2,1)*QGapPos(n1+n5,2)*QGapPos(n3+n4+n6,3)
    - 2.*QGapPos(n1,1)*QGapPos(n2+n5,2)*QGapPos(n3+n4+n6,3)+4.*QGapPos(n1+n2+n5,3)*QGapPos(n3+n4+n6,3)
    - 6.*QGapPos(n2,1)*QGapPos(n5,1)*QGapPos(n1+n3+n4+n6,4)+6.*QGapPos(n2+n5,2)*QGapPos(n1+n3+n4+n6,4)
    - 6.*QGapPos(n1,1)*QGapPos(n5,1)*QGapPos(n2+n3+n4+n6,4)+6.*QGapPos(n1+n5,2)*QGapPos(n2+n3+n4+n6,4)
    + 24.*QGapPos(n5,1)*QGapPos(n1+n2+n3+n4+n6,5)-QGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n3,1)*QGapPos(n4,1)*QGapPos(n5+n6,2)
    + QGapPos(n1+n2,2)*QGapPos(n3,1)*QGapPos(n4,1)*QGapPos(n5+n6,2)+QGapPos(n2,1)*QGapPos(n1+n3,2)*QGapPos(n4,1)*QGapPos(n5+n6,2)
    + QGapPos(n1,1)*QGapPos(n2+n3,2)*QGapPos(n4,1)*QGapPos(n5+n6,2)-2.*QGapPos(n1+n2+n3,3)*QGapPos(n4,1)*QGapPos(n5+n6,2)
    + QGapPos(n2,1)*QGapPos(n3,1)*QGapPos(n1+n4,2)*QGapPos(n5+n6,2)-QGapPos(n2+n3,2)*QGapPos(n1+n4,2)*QGapPos(n5+n6,2)
    + QGapPos(n1,1)*QGapPos(n3,1)*QGapPos(n2+n4,2)*QGapPos(n5+n6,2)-QGapPos(n1+n3,2)*QGapPos(n2+n4,2)*QGapPos(n5+n6,2)
    - 2.*QGapPos(n3,1)*QGapPos(n1+n2+n4,3)*QGapPos(n5+n6,2)+QGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n3+n4,2)*QGapPos(n5+n6,2)
    - QGapPos(n1+n2,2)*QGapPos(n3+n4,2)*QGapPos(n5+n6,2)-2.*QGapPos(n2,1)*QGapPos(n1+n3+n4,3)*QGapPos(n5+n6,2)
    - 2.*QGapPos(n1,1)*QGapPos(n2+n3+n4,3)*QGapPos(n5+n6,2)+6.*QGapPos(n1+n2+n3+n4,4)*QGapPos(n5+n6,2)
    + 2.*QGapPos(n2,1)*QGapPos(n3,1)*QGapPos(n4,1)*QGapPos(n1+n5+n6,3)-2.*QGapPos(n2+n3,2)*QGapPos(n4,1)*QGapPos(n1+n5+n6,3)
    - 2.*QGapPos(n3,1)*QGapPos(n2+n4,2)*QGapPos(n1+n5+n6,3)-2.*QGapPos(n2,1)*QGapPos(n3+n4,2)*QGapPos(n1+n5+n6,3)
    + 4.*QGapPos(n2+n3+n4,3)*QGapPos(n1+n5+n6,3)+2.*QGapPos(n1,1)*QGapPos(n3,1)*QGapPos(n4,1)*QGapPos(n2+n5+n6,3)
    - 2.*QGapPos(n1+n3,2)*QGapPos(n4,1)*QGapPos(n2+n5+n6,3)-2.*QGapPos(n3,1)*QGapPos(n1+n4,2)*QGapPos(n2+n5+n6,3)
    - 2.*QGapPos(n1,1)*QGapPos(n3+n4,2)*QGapPos(n2+n5+n6,3)+4.*QGapPos(n1+n3+n4,3)*QGapPos(n2+n5+n6,3)
    - 6.*QGapPos(n3,1)*QGapPos(n4,1)*QGapPos(n1+n2+n5+n6,4)+6.*QGapPos(n3+n4,2)*QGapPos(n1+n2+n5+n6,4)
    + 2.*QGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n4,1)*QGapPos(n3+n5+n6,3)-2.*QGapPos(n1+n2,2)*QGapPos(n4,1)*QGapPos(n3+n5+n6,3)
    - 2.*QGapPos(n2,1)*QGapPos(n1+n4,2)*QGapPos(n3+n5+n6,3)-2.*QGapPos(n1,1)*QGapPos(n2+n4,2)*QGapPos(n3+n5+n6,3)
    + 4.*QGapPos(n1+n2+n4,3)*QGapPos(n3+n5+n6,3)-6.*QGapPos(n2,1)*QGapPos(n4,1)*QGapPos(n1+n3+n5+n6,4)
    + 6.*QGapPos(n2+n4,2)*QGapPos(n1+n3+n5+n6,4)-6.*QGapPos(n1,1)*QGapPos(n4,1)*QGapPos(n2+n3+n5+n6,4)
    + 6.*QGapPos(n1+n4,2)*QGapPos(n2+n3+n5+n6,4)+24.*QGapPos(n4,1)*QGapPos(n1+n2+n3+n5+n6,5)
    + 2.*QGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n3,1)*QGapPos(n4+n5+n6,3)-2.*QGapPos(n1+n2,2)*QGapPos(n3,1)*QGapPos(n4+n5+n6,3)
    - 2.*QGapPos(n2,1)*QGapPos(n1+n3,2)*QGapPos(n4+n5+n6,3)-2.*QGapPos(n1,1)*QGapPos(n2+n3,2)*QGapPos(n4+n5+n6,3)
    + 4.*QGapPos(n1+n2+n3,3)*QGapPos(n4+n5+n6,3)-6.*QGapPos(n2,1)*QGapPos(n3,1)*QGapPos(n1+n4+n5+n6,4)
    + 6.*QGapPos(n2+n3,2)*QGapPos(n1+n4+n5+n6,4)-6.*QGapPos(n1,1)*QGapPos(n3,1)*QGapPos(n2+n4+n5+n6,4)
    + 6.*QGapPos(n1+n3,2)*QGapPos(n2+n4+n5+n6,4)+24.*QGapPos(n3,1)*QGapPos(n1+n2+n4+n5+n6,5)
    - 6.*QGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n3+n4+n5+n6,4)+6.*QGapPos(n1+n2,2)*QGapPos(n3+n4+n5+n6,4)
    + 24.*QGapPos(n2,1)*QGapPos(n1+n3+n4+n5+n6,5)+24.*QGapPos(n1,1)*QGapPos(n2+n3+n4+n5+n6,5)
    - 120.*QGapPos(n1+n2+n3+n4+n5+n6,6);
    return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::SixNeg(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t n5, const Int_t n6) const
{
  TComplex formula = QGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n3,1)*QGapNeg(n4,1)*QGapNeg(n5,1)*QGapNeg(n6,1)-QGapNeg(n1+n2,2)*QGapNeg(n3,1)*QGapNeg(n4,1)*QGapNeg(n5,1)*QGapNeg(n6,1)
    - QGapNeg(n2,1)*QGapNeg(n1+n3,2)*QGapNeg(n4,1)*QGapNeg(n5,1)*QGapNeg(n6,1)-QGapNeg(n1,1)*QGapNeg(n2+n3,2)*QGapNeg(n4,1)*QGapNeg(n5,1)*QGapNeg(n6,1)
    + 2.*QGapNeg(n1+n2+n3,3)*QGapNeg(n4,1)*QGapNeg(n5,1)*QGapNeg(n6,1)-QGapNeg(n2,1)*QGapNeg(n3,1)*QGapNeg(n1+n4,2)*QGapNeg(n5,1)*QGapNeg(n6,1)
    + QGapNeg(n2+n3,2)*QGapNeg(n1+n4,2)*QGapNeg(n5,1)*QGapNeg(n6,1)-QGapNeg(n1,1)*QGapNeg(n3,1)*QGapNeg(n2+n4,2)*QGapNeg(n5,1)*QGapNeg(n6,1)
    + QGapNeg(n1+n3,2)*QGapNeg(n2+n4,2)*QGapNeg(n5,1)*QGapNeg(n6,1)+2.*QGapNeg(n3,1)*QGapNeg(n1+n2+n4,3)*QGapNeg(n5,1)*QGapNeg(n6,1)
    - QGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n3+n4,2)*QGapNeg(n5,1)*QGapNeg(n6,1)+QGapNeg(n1+n2,2)*QGapNeg(n3+n4,2)*QGapNeg(n5,1)*QGapNeg(n6,1)
    + 2.*QGapNeg(n2,1)*QGapNeg(n1+n3+n4,3)*QGapNeg(n5,1)*QGapNeg(n6,1)+2.*QGapNeg(n1,1)*QGapNeg(n2+n3+n4,3)*QGapNeg(n5,1)*QGapNeg(n6,1)
    - 6.*QGapNeg(n1+n2+n3+n4,4)*QGapNeg(n5,1)*QGapNeg(n6,1)-QGapNeg(n2,1)*QGapNeg(n3,1)*QGapNeg(n4,1)*QGapNeg(n1+n5,2)*QGapNeg(n6,1)
    + QGapNeg(n2+n3,2)*QGapNeg(n4,1)*QGapNeg(n1+n5,2)*QGapNeg(n6,1)+QGapNeg(n3,1)*QGapNeg(n2+n4,2)*QGapNeg(n1+n5,2)*QGapNeg(n6,1)
    + QGapNeg(n2,1)*QGapNeg(n3+n4,2)*QGapNeg(n1+n5,2)*QGapNeg(n6,1)-2.*QGapNeg(n2+n3+n4,3)*QGapNeg(n1+n5,2)*QGapNeg(n6,1)
    - QGapNeg(n1,1)*QGapNeg(n3,1)*QGapNeg(n4,1)*QGapNeg(n2+n5,2)*QGapNeg(n6,1)+QGapNeg(n1+n3,2)*QGapNeg(n4,1)*QGapNeg(n2+n5,2)*QGapNeg(n6,1)
    + QGapNeg(n3,1)*QGapNeg(n1+n4,2)*QGapNeg(n2+n5,2)*QGapNeg(n6,1)+QGapNeg(n1,1)*QGapNeg(n3+n4,2)*QGapNeg(n2+n5,2)*QGapNeg(n6,1)
    - 2.*QGapNeg(n1+n3+n4,3)*QGapNeg(n2+n5,2)*QGapNeg(n6,1)+2.*QGapNeg(n3,1)*QGapNeg(n4,1)*QGapNeg(n1+n2+n5,3)*QGapNeg(n6,1)
    - 2.*QGapNeg(n3+n4,2)*QGapNeg(n1+n2+n5,3)*QGapNeg(n6,1)-QGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n4,1)*QGapNeg(n3+n5,2)*QGapNeg(n6,1)
    + QGapNeg(n1+n2,2)*QGapNeg(n4,1)*QGapNeg(n3+n5,2)*QGapNeg(n6,1)+QGapNeg(n2,1)*QGapNeg(n1+n4,2)*QGapNeg(n3+n5,2)*QGapNeg(n6,1)
    + QGapNeg(n1,1)*QGapNeg(n2+n4,2)*QGapNeg(n3+n5,2)*QGapNeg(n6,1)-2.*QGapNeg(n1+n2+n4,3)*QGapNeg(n3+n5,2)*QGapNeg(n6,1)
    + 2.*QGapNeg(n2,1)*QGapNeg(n4,1)*QGapNeg(n1+n3+n5,3)*QGapNeg(n6,1)-2.*QGapNeg(n2+n4,2)*QGapNeg(n1+n3+n5,3)*QGapNeg(n6,1)
    + 2.*QGapNeg(n1,1)*QGapNeg(n4,1)*QGapNeg(n2+n3+n5,3)*QGapNeg(n6,1)-2.*QGapNeg(n1+n4,2)*QGapNeg(n2+n3+n5,3)*QGapNeg(n6,1)
    - 6.*QGapNeg(n4,1)*QGapNeg(n1+n2+n3+n5,4)*QGapNeg(n6,1)-QGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n3,1)*QGapNeg(n4+n5,2)*QGapNeg(n6,1)
    + QGapNeg(n1+n2,2)*QGapNeg(n3,1)*QGapNeg(n4+n5,2)*QGapNeg(n6,1)+QGapNeg(n2,1)*QGapNeg(n1+n3,2)*QGapNeg(n4+n5,2)*QGapNeg(n6,1)
    + QGapNeg(n1,1)*QGapNeg(n2+n3,2)*QGapNeg(n4+n5,2)*QGapNeg(n6,1)-2.*QGapNeg(n1+n2+n3,3)*QGapNeg(n4+n5,2)*QGapNeg(n6,1)
    + 2.*QGapNeg(n2,1)*QGapNeg(n3,1)*QGapNeg(n1+n4+n5,3)*QGapNeg(n6,1)-2.*QGapNeg(n2+n3,2)*QGapNeg(n1+n4+n5,3)*QGapNeg(n6,1)
    + 2.*QGapNeg(n1,1)*QGapNeg(n3,1)*QGapNeg(n2+n4+n5,3)*QGapNeg(n6,1)-2.*QGapNeg(n1+n3,2)*QGapNeg(n2+n4+n5,3)*QGapNeg(n6,1)
    - 6.*QGapNeg(n3,1)*QGapNeg(n1+n2+n4+n5,4)*QGapNeg(n6,1)+2.*QGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n3+n4+n5,3)*QGapNeg(n6,1)
    - 2.*QGapNeg(n1+n2,2)*QGapNeg(n3+n4+n5,3)*QGapNeg(n6,1)-6.*QGapNeg(n2,1)*QGapNeg(n1+n3+n4+n5,4)*QGapNeg(n6,1)
    - 6.*QGapNeg(n1,1)*QGapNeg(n2+n3+n4+n5,4)*QGapNeg(n6,1)+24.*QGapNeg(n1+n2+n3+n4+n5,5)*QGapNeg(n6,1)
    - QGapNeg(n2,1)*QGapNeg(n3,1)*QGapNeg(n4,1)*QGapNeg(n5,1)*QGapNeg(n1+n6,2)+QGapNeg(n2+n3,2)*QGapNeg(n4,1)*QGapNeg(n5,1)*QGapNeg(n1+n6,2)
    + QGapNeg(n3,1)*QGapNeg(n2+n4,2)*QGapNeg(n5,1)*QGapNeg(n1+n6,2)+QGapNeg(n2,1)*QGapNeg(n3+n4,2)*QGapNeg(n5,1)*QGapNeg(n1+n6,2)
    - 2.*QGapNeg(n2+n3+n4,3)*QGapNeg(n5,1)*QGapNeg(n1+n6,2)+QGapNeg(n3,1)*QGapNeg(n4,1)*QGapNeg(n2+n5,2)*QGapNeg(n1+n6,2)
    - QGapNeg(n3+n4,2)*QGapNeg(n2+n5,2)*QGapNeg(n1+n6,2)+QGapNeg(n2,1)*QGapNeg(n4,1)*QGapNeg(n3+n5,2)*QGapNeg(n1+n6,2)
    - QGapNeg(n2+n4,2)*QGapNeg(n3+n5,2)*QGapNeg(n1+n6,2)-2.*QGapNeg(n4,1)*QGapNeg(n2+n3+n5,3)*QGapNeg(n1+n6,2)
    + QGapNeg(n2,1)*QGapNeg(n3,1)*QGapNeg(n4+n5,2)*QGapNeg(n1+n6,2)-QGapNeg(n2+n3,2)*QGapNeg(n4+n5,2)*QGapNeg(n1+n6,2)
    - 2.*QGapNeg(n3,1)*QGapNeg(n2+n4+n5,3)*QGapNeg(n1+n6,2)-2.*QGapNeg(n2,1)*QGapNeg(n3+n4+n5,3)*QGapNeg(n1+n6,2)
    + 6.*QGapNeg(n2+n3+n4+n5,4)*QGapNeg(n1+n6,2)-QGapNeg(n1,1)*QGapNeg(n3,1)*QGapNeg(n4,1)*QGapNeg(n5,1)*QGapNeg(n2+n6,2)
    + QGapNeg(n1+n3,2)*QGapNeg(n4,1)*QGapNeg(n5,1)*QGapNeg(n2+n6,2)+QGapNeg(n3,1)*QGapNeg(n1+n4,2)*QGapNeg(n5,1)*QGapNeg(n2+n6,2)
    + QGapNeg(n1,1)*QGapNeg(n3+n4,2)*QGapNeg(n5,1)*QGapNeg(n2+n6,2)-2.*QGapNeg(n1+n3+n4,3)*QGapNeg(n5,1)*QGapNeg(n2+n6,2)
    + QGapNeg(n3,1)*QGapNeg(n4,1)*QGapNeg(n1+n5,2)*QGapNeg(n2+n6,2)-QGapNeg(n3+n4,2)*QGapNeg(n1+n5,2)*QGapNeg(n2+n6,2)
    + QGapNeg(n1,1)*QGapNeg(n4,1)*QGapNeg(n3+n5,2)*QGapNeg(n2+n6,2)-QGapNeg(n1+n4,2)*QGapNeg(n3+n5,2)*QGapNeg(n2+n6,2)
    - 2.*QGapNeg(n4,1)*QGapNeg(n1+n3+n5,3)*QGapNeg(n2+n6,2)+QGapNeg(n1,1)*QGapNeg(n3,1)*QGapNeg(n4+n5,2)*QGapNeg(n2+n6,2)
    - QGapNeg(n1+n3,2)*QGapNeg(n4+n5,2)*QGapNeg(n2+n6,2)-2.*QGapNeg(n3,1)*QGapNeg(n1+n4+n5,3)*QGapNeg(n2+n6,2)
    - 2.*QGapNeg(n1,1)*QGapNeg(n3+n4+n5,3)*QGapNeg(n2+n6,2)+6.*QGapNeg(n1+n3+n4+n5,4)*QGapNeg(n2+n6,2)
    + 2.*QGapNeg(n3,1)*QGapNeg(n4,1)*QGapNeg(n5,1)*QGapNeg(n1+n2+n6,3)-2.*QGapNeg(n3+n4,2)*QGapNeg(n5,1)*QGapNeg(n1+n2+n6,3)
    - 2.*QGapNeg(n4,1)*QGapNeg(n3+n5,2)*QGapNeg(n1+n2+n6,3)-2.*QGapNeg(n3,1)*QGapNeg(n4+n5,2)*QGapNeg(n1+n2+n6,3)
    + 4.*QGapNeg(n3+n4+n5,3)*QGapNeg(n1+n2+n6,3)-QGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n4,1)*QGapNeg(n5,1)*QGapNeg(n3+n6,2)
    + QGapNeg(n1+n2,2)*QGapNeg(n4,1)*QGapNeg(n5,1)*QGapNeg(n3+n6,2)+QGapNeg(n2,1)*QGapNeg(n1+n4,2)*QGapNeg(n5,1)*QGapNeg(n3+n6,2)
    + QGapNeg(n1,1)*QGapNeg(n2+n4,2)*QGapNeg(n5,1)*QGapNeg(n3+n6,2)-2.*QGapNeg(n1+n2+n4,3)*QGapNeg(n5,1)*QGapNeg(n3+n6,2)
    + QGapNeg(n2,1)*QGapNeg(n4,1)*QGapNeg(n1+n5,2)*QGapNeg(n3+n6,2)-QGapNeg(n2+n4,2)*QGapNeg(n1+n5,2)*QGapNeg(n3+n6,2)
    + QGapNeg(n1,1)*QGapNeg(n4,1)*QGapNeg(n2+n5,2)*QGapNeg(n3+n6,2)-QGapNeg(n1+n4,2)*QGapNeg(n2+n5,2)*QGapNeg(n3+n6,2)
    - 2.*QGapNeg(n4,1)*QGapNeg(n1+n2+n5,3)*QGapNeg(n3+n6,2)+QGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n4+n5,2)*QGapNeg(n3+n6,2)
    - QGapNeg(n1+n2,2)*QGapNeg(n4+n5,2)*QGapNeg(n3+n6,2)-2.*QGapNeg(n2,1)*QGapNeg(n1+n4+n5,3)*QGapNeg(n3+n6,2)
    - 2.*QGapNeg(n1,1)*QGapNeg(n2+n4+n5,3)*QGapNeg(n3+n6,2)+6.*QGapNeg(n1+n2+n4+n5,4)*QGapNeg(n3+n6,2)
    + 2.*QGapNeg(n2,1)*QGapNeg(n4,1)*QGapNeg(n5,1)*QGapNeg(n1+n3+n6,3)-2.*QGapNeg(n2+n4,2)*QGapNeg(n5,1)*QGapNeg(n1+n3+n6,3)
    - 2.*QGapNeg(n4,1)*QGapNeg(n2+n5,2)*QGapNeg(n1+n3+n6,3)-2.*QGapNeg(n2,1)*QGapNeg(n4+n5,2)*QGapNeg(n1+n3+n6,3)
    + 4.*QGapNeg(n2+n4+n5,3)*QGapNeg(n1+n3+n6,3)+2.*QGapNeg(n1,1)*QGapNeg(n4,1)*QGapNeg(n5,1)*QGapNeg(n2+n3+n6,3)
    - 2.*QGapNeg(n1+n4,2)*QGapNeg(n5,1)*QGapNeg(n2+n3+n6,3)-2.*QGapNeg(n4,1)*QGapNeg(n1+n5,2)*QGapNeg(n2+n3+n6,3)
    - 2.*QGapNeg(n1,1)*QGapNeg(n4+n5,2)*QGapNeg(n2+n3+n6,3)+4.*QGapNeg(n1+n4+n5,3)*QGapNeg(n2+n3+n6,3)
    - 6.*QGapNeg(n4,1)*QGapNeg(n5,1)*QGapNeg(n1+n2+n3+n6,4)+6.*QGapNeg(n4+n5,2)*QGapNeg(n1+n2+n3+n6,4)
    - QGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n3,1)*QGapNeg(n5,1)*QGapNeg(n4+n6,2)+QGapNeg(n1+n2,2)*QGapNeg(n3,1)*QGapNeg(n5,1)*QGapNeg(n4+n6,2)
    + QGapNeg(n2,1)*QGapNeg(n1+n3,2)*QGapNeg(n5,1)*QGapNeg(n4+n6,2)+QGapNeg(n1,1)*QGapNeg(n2+n3,2)*QGapNeg(n5,1)*QGapNeg(n4+n6,2)
    - 2.*QGapNeg(n1+n2+n3,3)*QGapNeg(n5,1)*QGapNeg(n4+n6,2)+QGapNeg(n2,1)*QGapNeg(n3,1)*QGapNeg(n1+n5,2)*QGapNeg(n4+n6,2)
    - QGapNeg(n2+n3,2)*QGapNeg(n1+n5,2)*QGapNeg(n4+n6,2)+QGapNeg(n1,1)*QGapNeg(n3,1)*QGapNeg(n2+n5,2)*QGapNeg(n4+n6,2)
    - QGapNeg(n1+n3,2)*QGapNeg(n2+n5,2)*QGapNeg(n4+n6,2)-2.*QGapNeg(n3,1)*QGapNeg(n1+n2+n5,3)*QGapNeg(n4+n6,2)
    + QGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n3+n5,2)*QGapNeg(n4+n6,2)-QGapNeg(n1+n2,2)*QGapNeg(n3+n5,2)*QGapNeg(n4+n6,2)
    - 2.*QGapNeg(n2,1)*QGapNeg(n1+n3+n5,3)*QGapNeg(n4+n6,2)-2.*QGapNeg(n1,1)*QGapNeg(n2+n3+n5,3)*QGapNeg(n4+n6,2)
    + 6.*QGapNeg(n1+n2+n3+n5,4)*QGapNeg(n4+n6,2)+2.*QGapNeg(n2,1)*QGapNeg(n3,1)*QGapNeg(n5,1)*QGapNeg(n1+n4+n6,3)
    - 2.*QGapNeg(n2+n3,2)*QGapNeg(n5,1)*QGapNeg(n1+n4+n6,3)-2.*QGapNeg(n3,1)*QGapNeg(n2+n5,2)*QGapNeg(n1+n4+n6,3)
    - 2.*QGapNeg(n2,1)*QGapNeg(n3+n5,2)*QGapNeg(n1+n4+n6,3)+4.*QGapNeg(n2+n3+n5,3)*QGapNeg(n1+n4+n6,3)
    + 2.*QGapNeg(n1,1)*QGapNeg(n3,1)*QGapNeg(n5,1)*QGapNeg(n2+n4+n6,3)-2.*QGapNeg(n1+n3,2)*QGapNeg(n5,1)*QGapNeg(n2+n4+n6,3)
    - 2.*QGapNeg(n3,1)*QGapNeg(n1+n5,2)*QGapNeg(n2+n4+n6,3)-2.*QGapNeg(n1,1)*QGapNeg(n3+n5,2)*QGapNeg(n2+n4+n6,3)
    + 4.*QGapNeg(n1+n3+n5,3)*QGapNeg(n2+n4+n6,3)-6.*QGapNeg(n3,1)*QGapNeg(n5,1)*QGapNeg(n1+n2+n4+n6,4)
    + 6.*QGapNeg(n3+n5,2)*QGapNeg(n1+n2+n4+n6,4)+2.*QGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n5,1)*QGapNeg(n3+n4+n6,3)
    - 2.*QGapNeg(n1+n2,2)*QGapNeg(n5,1)*QGapNeg(n3+n4+n6,3)-2.*QGapNeg(n2,1)*QGapNeg(n1+n5,2)*QGapNeg(n3+n4+n6,3)
    - 2.*QGapNeg(n1,1)*QGapNeg(n2+n5,2)*QGapNeg(n3+n4+n6,3)+4.*QGapNeg(n1+n2+n5,3)*QGapNeg(n3+n4+n6,3)
    - 6.*QGapNeg(n2,1)*QGapNeg(n5,1)*QGapNeg(n1+n3+n4+n6,4)+6.*QGapNeg(n2+n5,2)*QGapNeg(n1+n3+n4+n6,4)
    - 6.*QGapNeg(n1,1)*QGapNeg(n5,1)*QGapNeg(n2+n3+n4+n6,4)+6.*QGapNeg(n1+n5,2)*QGapNeg(n2+n3+n4+n6,4)
    + 24.*QGapNeg(n5,1)*QGapNeg(n1+n2+n3+n4+n6,5)-QGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n3,1)*QGapNeg(n4,1)*QGapNeg(n5+n6,2)
    + QGapNeg(n1+n2,2)*QGapNeg(n3,1)*QGapNeg(n4,1)*QGapNeg(n5+n6,2)+QGapNeg(n2,1)*QGapNeg(n1+n3,2)*QGapNeg(n4,1)*QGapNeg(n5+n6,2)
    + QGapNeg(n1,1)*QGapNeg(n2+n3,2)*QGapNeg(n4,1)*QGapNeg(n5+n6,2)-2.*QGapNeg(n1+n2+n3,3)*QGapNeg(n4,1)*QGapNeg(n5+n6,2)
    + QGapNeg(n2,1)*QGapNeg(n3,1)*QGapNeg(n1+n4,2)*QGapNeg(n5+n6,2)-QGapNeg(n2+n3,2)*QGapNeg(n1+n4,2)*QGapNeg(n5+n6,2)
    + QGapNeg(n1,1)*QGapNeg(n3,1)*QGapNeg(n2+n4,2)*QGapNeg(n5+n6,2)-QGapNeg(n1+n3,2)*QGapNeg(n2+n4,2)*QGapNeg(n5+n6,2)
    - 2.*QGapNeg(n3,1)*QGapNeg(n1+n2+n4,3)*QGapNeg(n5+n6,2)+QGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n3+n4,2)*QGapNeg(n5+n6,2)
    - QGapNeg(n1+n2,2)*QGapNeg(n3+n4,2)*QGapNeg(n5+n6,2)-2.*QGapNeg(n2,1)*QGapNeg(n1+n3+n4,3)*QGapNeg(n5+n6,2)
    - 2.*QGapNeg(n1,1)*QGapNeg(n2+n3+n4,3)*QGapNeg(n5+n6,2)+6.*QGapNeg(n1+n2+n3+n4,4)*QGapNeg(n5+n6,2)
    + 2.*QGapNeg(n2,1)*QGapNeg(n3,1)*QGapNeg(n4,1)*QGapNeg(n1+n5+n6,3)-2.*QGapNeg(n2+n3,2)*QGapNeg(n4,1)*QGapNeg(n1+n5+n6,3)
    - 2.*QGapNeg(n3,1)*QGapNeg(n2+n4,2)*QGapNeg(n1+n5+n6,3)-2.*QGapNeg(n2,1)*QGapNeg(n3+n4,2)*QGapNeg(n1+n5+n6,3)
    + 4.*QGapNeg(n2+n3+n4,3)*QGapNeg(n1+n5+n6,3)+2.*QGapNeg(n1,1)*QGapNeg(n3,1)*QGapNeg(n4,1)*QGapNeg(n2+n5+n6,3)
    - 2.*QGapNeg(n1+n3,2)*QGapNeg(n4,1)*QGapNeg(n2+n5+n6,3)-2.*QGapNeg(n3,1)*QGapNeg(n1+n4,2)*QGapNeg(n2+n5+n6,3)
    - 2.*QGapNeg(n1,1)*QGapNeg(n3+n4,2)*QGapNeg(n2+n5+n6,3)+4.*QGapNeg(n1+n3+n4,3)*QGapNeg(n2+n5+n6,3)
    - 6.*QGapNeg(n3,1)*QGapNeg(n4,1)*QGapNeg(n1+n2+n5+n6,4)+6.*QGapNeg(n3+n4,2)*QGapNeg(n1+n2+n5+n6,4)
    + 2.*QGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n4,1)*QGapNeg(n3+n5+n6,3)-2.*QGapNeg(n1+n2,2)*QGapNeg(n4,1)*QGapNeg(n3+n5+n6,3)
    - 2.*QGapNeg(n2,1)*QGapNeg(n1+n4,2)*QGapNeg(n3+n5+n6,3)-2.*QGapNeg(n1,1)*QGapNeg(n2+n4,2)*QGapNeg(n3+n5+n6,3)
    + 4.*QGapNeg(n1+n2+n4,3)*QGapNeg(n3+n5+n6,3)-6.*QGapNeg(n2,1)*QGapNeg(n4,1)*QGapNeg(n1+n3+n5+n6,4)
    + 6.*QGapNeg(n2+n4,2)*QGapNeg(n1+n3+n5+n6,4)-6.*QGapNeg(n1,1)*QGapNeg(n4,1)*QGapNeg(n2+n3+n5+n6,4)
    + 6.*QGapNeg(n1+n4,2)*QGapNeg(n2+n3+n5+n6,4)+24.*QGapNeg(n4,1)*QGapNeg(n1+n2+n3+n5+n6,5)
    + 2.*QGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n3,1)*QGapNeg(n4+n5+n6,3)-2.*QGapNeg(n1+n2,2)*QGapNeg(n3,1)*QGapNeg(n4+n5+n6,3)
    - 2.*QGapNeg(n2,1)*QGapNeg(n1+n3,2)*QGapNeg(n4+n5+n6,3)-2.*QGapNeg(n1,1)*QGapNeg(n2+n3,2)*QGapNeg(n4+n5+n6,3)
    + 4.*QGapNeg(n1+n2+n3,3)*QGapNeg(n4+n5+n6,3)-6.*QGapNeg(n2,1)*QGapNeg(n3,1)*QGapNeg(n1+n4+n5+n6,4)
    + 6.*QGapNeg(n2+n3,2)*QGapNeg(n1+n4+n5+n6,4)-6.*QGapNeg(n1,1)*QGapNeg(n3,1)*QGapNeg(n2+n4+n5+n6,4)
    + 6.*QGapNeg(n1+n3,2)*QGapNeg(n2+n4+n5+n6,4)+24.*QGapNeg(n3,1)*QGapNeg(n1+n2+n4+n5+n6,5)
    - 6.*QGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n3+n4+n5+n6,4)+6.*QGapNeg(n1+n2,2)*QGapNeg(n3+n4+n5+n6,4)
    + 24.*QGapNeg(n2,1)*QGapNeg(n1+n3+n4+n5+n6,5)+24.*QGapNeg(n1,1)*QGapNeg(n2+n3+n4+n5+n6,5)
    - 120.*QGapNeg(n1+n2+n3+n4+n5+n6,6);
    return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::Seven(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t n5, const Int_t n6, const Int_t n7) const
{
  TComplex Correlation = {0, 0};
    Int_t Narray[] = {n1, n2, n3, n4, n5, n6};

    for(Int_t k=7; k-- >0; )
    {// backward loop of k from m-1 until 0, where m is the m-particle correlation, in this case m=4

        Int_t array[6] = {0,1,2,3,4,5};
        Int_t iPerm = 0;
        Int_t count = 0;

        // k==6: there is just one combination, we can add it manually
        if(k==6){
            Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
            Six(n1, n2, n3, n4, n5, n6)*Q(n7, 7-k);
        }// k==6

        else if(k==5){
            do{
                iPerm += 1;
                if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4]){
                    count += 1;
                    Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
                    Five(Narray[Int_t(array[0])], Narray[Int_t(array[1])], Narray[Int_t(array[2])],
                         Narray[Int_t(array[3])], Narray[Int_t(array[4])])*
                    Q(Narray[int(array[5])]+n7, 7-k);
                }
            }while(std::next_permutation(array, array+6));
        }// k==5

        else if(k==4){
            do{
                iPerm += 1;
                if(iPerm%2 == 1){
                    if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3]){
                        Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
                        Four(Narray[Int_t(array[0])], Narray[Int_t(array[1])], Narray[Int_t(array[2])],
                             Narray[Int_t(array[3])])*
                        Q(Narray[Int_t(array[4])]+Narray[Int_t(array[5])]+n7, 7-k);
                    }
                }
            }while(std::next_permutation(array, array+6));
        }// k==4

        else if(k==3){
            do{
                iPerm += 1;
                if(iPerm%6 == 1){
                    if(array[0] < array[1] && array[1] < array[2]){
                        Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
                        Three(Narray[Int_t(array[0])], Narray[Int_t(array[1])], Narray[Int_t(array[2])])*
                        Q(Narray[Int_t(array[3])]+Narray[Int_t(array[4])]+Narray[Int_t(array[5])]+n7, 7-k);
                    }
                }
            }while(std::next_permutation(array, array+6));
        }// k==3

        else if(k==2){
            do{
                iPerm += 1;
                if(iPerm%24 == 1){
                    if(array[0] < array[1]){
                        Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
                        Two(Narray[Int_t(array[0])], Narray[Int_t(array[1])])*
                        Q(Narray[Int_t(array[2])]+Narray[Int_t(array[3])]+Narray[Int_t(array[4])]
                          +Narray[Int_t(array[5])]+n7, 7-k);
                    }
                }
            }while(std::next_permutation(array, array+6));
        }// k==2

        else if(k == 1){
            Correlation = Correlation
            + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n1, 1)*Q(n2+n3+n4+n5+n6+n7, 7-k)
            + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n2, 1)*Q(n1+n3+n4+n5+n6+n7, 7-k)
            + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n3, 1)*Q(n1+n2+n4+n5+n6+n7, 7-k)
            + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n4, 1)*Q(n1+n2+n3+n5+n6+n7, 7-k)
            + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n5, 1)*Q(n1+n2+n3+n4+n6+n7, 7-k)
            + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n6, 1)*Q(n1+n2+n3+n4+n5+n7, 7-k);
        }// k==1

        else if(k == 0){
            Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n1+n2+n3+n4+n5+n6+n7, 7-k);
        }// k==0

        else{
            printf("Invalid range of k in Seven()\n");;
            return {0,0};
        }

    }// loop over k
    return Correlation;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::SevenPos(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t n5, const Int_t n6, const Int_t n7) const
{
  TComplex Correlation = {0, 0};
    Int_t Narray[] = {n1, n2, n3, n4, n5, n6};

    for(Int_t k=7; k-- >0; )
    {// backward loop of k from m-1 until 0, where m is the m-particle correlation, in this case m=4

        Int_t array[6] = {0,1,2,3,4,5};
        Int_t iPerm = 0;
        Int_t count = 0;

        // k==6: there is just one combination, we can add it manually
        if(k==6){
            Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
            SixPos(n1, n2, n3, n4, n5, n6)*QGapPos(n7, 7-k);
        }// k==6

        else if(k==5){
            do{
                iPerm += 1;
                if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4]){
                    count += 1;
                    Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
                    FivePos(Narray[Int_t(array[0])], Narray[Int_t(array[1])], Narray[Int_t(array[2])],
                         Narray[Int_t(array[3])], Narray[Int_t(array[4])])*
                    QGapPos(Narray[int(array[5])]+n7, 7-k);
                }
            }while(std::next_permutation(array, array+6));
        }// k==5

        else if(k==4){
            do{
                iPerm += 1;
                if(iPerm%2 == 1){
                    if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3]){
                        Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
                        FourPos(Narray[Int_t(array[0])], Narray[Int_t(array[1])], Narray[Int_t(array[2])],
                             Narray[Int_t(array[3])])*
                        QGapPos(Narray[Int_t(array[4])]+Narray[Int_t(array[5])]+n7, 7-k);
                    }
                }
            }while(std::next_permutation(array, array+6));
        }// k==4

        else if(k==3){
            do{
                iPerm += 1;
                if(iPerm%6 == 1){
                    if(array[0] < array[1] && array[1] < array[2]){
                        Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
                        ThreePos(Narray[Int_t(array[0])], Narray[Int_t(array[1])], Narray[Int_t(array[2])])*
                        QGapPos(Narray[Int_t(array[3])]+Narray[Int_t(array[4])]+Narray[Int_t(array[5])]+n7, 7-k);
                    }
                }
            }while(std::next_permutation(array, array+6));
        }// k==3

        else if(k==2){
            do{
                iPerm += 1;
                if(iPerm%24 == 1){
                    if(array[0] < array[1]){
                        Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
                        TwoPos(Narray[Int_t(array[0])], Narray[Int_t(array[1])])*
                        QGapPos(Narray[Int_t(array[2])]+Narray[Int_t(array[3])]+Narray[Int_t(array[4])]
                          +Narray[Int_t(array[5])]+n7, 7-k);
                    }
                }
            }while(std::next_permutation(array, array+6));
        }// k==2

        else if(k == 1){
            Correlation = Correlation
            + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*QGapPos(n1, 1)*QGapPos(n2+n3+n4+n5+n6+n7, 7-k)
            + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*QGapPos(n2, 1)*QGapPos(n1+n3+n4+n5+n6+n7, 7-k)
            + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*QGapPos(n3, 1)*QGapPos(n1+n2+n4+n5+n6+n7, 7-k)
            + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*QGapPos(n4, 1)*QGapPos(n1+n2+n3+n5+n6+n7, 7-k)
            + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*QGapPos(n5, 1)*QGapPos(n1+n2+n3+n4+n6+n7, 7-k)
            + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*QGapPos(n6, 1)*QGapPos(n1+n2+n3+n4+n5+n7, 7-k);
        }// k==1

        else if(k == 0){
            Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*QGapPos(n1+n2+n3+n4+n5+n6+n7, 7-k);
        }// k==0

        else{
            printf("Invalid range of k in Seven()\n");;
            return {0,0};
        }

    }// loop over k
    return Correlation;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::SevenNeg(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t n5, const Int_t n6, const Int_t n7) const
{
  TComplex Correlation = {0, 0};
    Int_t Narray[] = {n1, n2, n3, n4, n5, n6};

    for(Int_t k=7; k-- >0; )
    {// backward loop of k from m-1 until 0, where m is the m-particle correlation, in this case m=4

        Int_t array[6] = {0,1,2,3,4,5};
        Int_t iPerm = 0;
        Int_t count = 0;

        // k==6: there is just one combination, we can add it manually
        if(k==6){
            Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
            SixNeg(n1, n2, n3, n4, n5, n6)*QGapNeg(n7, 7-k);
        }// k==6

        else if(k==5){
            do{
                iPerm += 1;
                if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4]){
                    count += 1;
                    Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
                    FiveNeg(Narray[Int_t(array[0])], Narray[Int_t(array[1])], Narray[Int_t(array[2])],
                         Narray[Int_t(array[3])], Narray[Int_t(array[4])])*
                    QGapNeg(Narray[int(array[5])]+n7, 7-k);
                }
            }while(std::next_permutation(array, array+6));
        }// k==5

        else if(k==4){
            do{
                iPerm += 1;
                if(iPerm%2 == 1){
                    if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3]){
                        Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
                        FourNeg(Narray[Int_t(array[0])], Narray[Int_t(array[1])], Narray[Int_t(array[2])],
                             Narray[Int_t(array[3])])*
                        QGapNeg(Narray[Int_t(array[4])]+Narray[Int_t(array[5])]+n7, 7-k);
                    }
                }
            }while(std::next_permutation(array, array+6));
        }// k==4

        else if(k==3){
            do{
                iPerm += 1;
                if(iPerm%6 == 1){
                    if(array[0] < array[1] && array[1] < array[2]){
                        Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
                        ThreeNeg(Narray[Int_t(array[0])], Narray[Int_t(array[1])], Narray[Int_t(array[2])])*
                        QGapNeg(Narray[Int_t(array[3])]+Narray[Int_t(array[4])]+Narray[Int_t(array[5])]+n7, 7-k);
                    }
                }
            }while(std::next_permutation(array, array+6));
        }// k==3

        else if(k==2){
            do{
                iPerm += 1;
                if(iPerm%24 == 1){
                    if(array[0] < array[1]){
                        Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
                        TwoNeg(Narray[Int_t(array[0])], Narray[Int_t(array[1])])*
                        QGapNeg(Narray[Int_t(array[2])]+Narray[Int_t(array[3])]+Narray[Int_t(array[4])]
                          +Narray[Int_t(array[5])]+n7, 7-k);
                    }
                }
            }while(std::next_permutation(array, array+6));
        }// k==2

        else if(k == 1){
            Correlation = Correlation
            + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*QGapNeg(n1, 1)*QGapNeg(n2+n3+n4+n5+n6+n7, 7-k)
            + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*QGapNeg(n2, 1)*QGapNeg(n1+n3+n4+n5+n6+n7, 7-k)
            + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*QGapNeg(n3, 1)*QGapNeg(n1+n2+n4+n5+n6+n7, 7-k)
            + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*QGapNeg(n4, 1)*QGapNeg(n1+n2+n3+n5+n6+n7, 7-k)
            + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*QGapNeg(n5, 1)*QGapNeg(n1+n2+n3+n4+n6+n7, 7-k)
            + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*QGapNeg(n6, 1)*QGapNeg(n1+n2+n3+n4+n5+n7, 7-k);
        }// k==1

        else if(k == 0){
            Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*QGapNeg(n1+n2+n3+n4+n5+n6+n7, 7-k);
        }// k==0

        else{
            printf("Invalid range of k in Seven()\n");;
            return {0,0};
        }

    }// loop over k
    return Correlation;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::Eight(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t n5, const Int_t n6, const Int_t n7, const Int_t n8) const
{
  TComplex Correlation = {0, 0};
    Int_t Narray[] = {n1, n2, n3, n4, n5, n6, n7};

    for(Int_t k=8; k-->0; )
    {// backward loop of k from m-1 until 0, where m is the m-particle correlation, in this case m=4

        Int_t array[7] = {0,1,2,3,4,5,6};
        Int_t iPerm = 0;
        Int_t count = 0;

        // k==7: there is just one combination, we can add it manually
        if(k==7){
            Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
            Seven(n1, n2, n3, n4, n5, n6, n7)*Q(n8, 8-k);
        }// k==7

        else if(k==6){
            do{
                iPerm += 1;
                if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4] && array[4] < array[5]){
                    count += 1;
                    Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
                    Six(Narray[Int_t(array[0])], Narray[Int_t(array[1])], Narray[Int_t(array[2])],
                        Narray[Int_t(array[3])], Narray[Int_t(array[4])], Narray[Int_t(array[5])])*
                    Q(Narray[Int_t(array[6])]+n8, 8-k);
                }
            }while(std::next_permutation(array, array+7));
        }// k==6

        else if(k==5){
            do{
                iPerm += 1;
                if(iPerm%2 == 1){
                    if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4]){
                        Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
                        Five(Narray[Int_t(array[0])], Narray[Int_t(array[1])], Narray[Int_t(array[2])],
                             Narray[Int_t(array[3])], Narray[Int_t(array[4])])*
                        Q(Narray[Int_t(array[5])]+Narray[Int_t(array[6])]+n8, 8-k);
                    }
                }
            }while(std::next_permutation(array, array+7));
        }// k==5

        else if(k==4){
            do{
                iPerm += 1;
                if(iPerm%6 == 1){
                    if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3]){
                        Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
                        Four(Narray[Int_t(array[0])], Narray[Int_t(array[1])], Narray[Int_t(array[2])], Narray[Int_t(array[3])])*
                        Q(Narray[Int_t(array[4])]+Narray[Int_t(array[5])]+Narray[Int_t(array[6])]+n8, 8-k);
                    }
                }
            }while(std::next_permutation(array, array+7));
        }// k==4

        else if(k==3){
            do{
                iPerm += 1;
                if(iPerm%24 == 1){
                    if(array[0] < array[1] && array[1] < array[2]){
                        Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
                        Three(Narray[Int_t(array[0])], Narray[Int_t(array[1])], Narray[Int_t(array[2])])*
                        Q(Narray[Int_t(array[3])]+Narray[Int_t(array[4])]+Narray[Int_t(array[5])]+Narray[Int_t(array[6])]+n8, 8-k);
                    }
                }
            }while(std::next_permutation(array, array+7));
        }// k==3

        else if(k==2){
            do{
                iPerm += 1;
                if(iPerm%120 == 1){
                    if(array[0] < array[1]){
                        Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
                        Two(Narray[Int_t(array[0])], Narray[Int_t(array[1])])*
                        Q(Narray[Int_t(array[2])]+Narray[Int_t(array[3])]+Narray[Int_t(array[4])]
                          +Narray[Int_t(array[5])]+Narray[Int_t(array[6])]+n8, 8-k);
                    }
                }
            }while(std::next_permutation(array, array+7));
        }// k==2

        else if(k == 1){
            Correlation = Correlation
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n1, 1)*Q(n2+n3+n4+n5+n6+n7+n8, 8-k)
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n2, 1)*Q(n1+n3+n4+n5+n6+n7+n8, 8-k)
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n3, 1)*Q(n1+n2+n4+n5+n6+n7+n8, 8-k)
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n4, 1)*Q(n1+n2+n3+n5+n6+n7+n8, 8-k)
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n5, 1)*Q(n1+n2+n3+n4+n6+n7+n8, 8-k)
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n6, 1)*Q(n1+n2+n3+n4+n5+n7+n8, 8-k)
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n7, 1)*Q(n1+n2+n3+n4+n5+n6+n8, 8-k);
        }// k==1

        else if(k == 0){
            Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n1+n2+n3+n4+n5+n6+n7+n8, 8-k);
        }// k==0

        else{
            printf("Invalid range of k in Eight() \n");
            return {0,0};
        }

    }// loop over k

    return Correlation;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::EightPos(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t n5, const Int_t n6, const Int_t n7, const Int_t n8) const
{
  TComplex Correlation = {0, 0};
    Int_t Narray[] = {n1, n2, n3, n4, n5, n6, n7};

    for(Int_t k=8; k-->0; )
    {// backward loop of k from m-1 until 0, where m is the m-particle correlation, in this case m=4

        Int_t array[7] = {0,1,2,3,4,5,6};
        Int_t iPerm = 0;
        Int_t count = 0;

        // k==7: there is just one combination, we can add it manually
        if(k==7){
            Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
            SevenPos(n1, n2, n3, n4, n5, n6, n7)*QGapPos(n8, 8-k);
        }// k==7

        else if(k==6){
            do{
                iPerm += 1;
                if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4] && array[4] < array[5]){
                    count += 1;
                    Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
                    SixPos(Narray[Int_t(array[0])], Narray[Int_t(array[1])], Narray[Int_t(array[2])],
                        Narray[Int_t(array[3])], Narray[Int_t(array[4])], Narray[Int_t(array[5])])*
                    QGapPos(Narray[Int_t(array[6])]+n8, 8-k);
                }
            }while(std::next_permutation(array, array+7));
        }// k==6

        else if(k==5){
            do{
                iPerm += 1;
                if(iPerm%2 == 1){
                    if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4]){
                        Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
                        FivePos(Narray[Int_t(array[0])], Narray[Int_t(array[1])], Narray[Int_t(array[2])],
                             Narray[Int_t(array[3])], Narray[Int_t(array[4])])*
                        QGapPos(Narray[Int_t(array[5])]+Narray[Int_t(array[6])]+n8, 8-k);
                    }
                }
            }while(std::next_permutation(array, array+7));
        }// k==5

        else if(k==4){
            do{
                iPerm += 1;
                if(iPerm%6 == 1){
                    if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3]){
                        Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
                        FourPos(Narray[Int_t(array[0])], Narray[Int_t(array[1])], Narray[Int_t(array[2])], Narray[Int_t(array[3])])*
                        QGapPos(Narray[Int_t(array[4])]+Narray[Int_t(array[5])]+Narray[Int_t(array[6])]+n8, 8-k);
                    }
                }
            }while(std::next_permutation(array, array+7));
        }// k==4

        else if(k==3){
            do{
                iPerm += 1;
                if(iPerm%24 == 1){
                    if(array[0] < array[1] && array[1] < array[2]){
                        Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
                        ThreePos(Narray[Int_t(array[0])], Narray[Int_t(array[1])], Narray[Int_t(array[2])])*
                        QGapPos(Narray[Int_t(array[3])]+Narray[Int_t(array[4])]+Narray[Int_t(array[5])]+Narray[Int_t(array[6])]+n8, 8-k);
                    }
                }
            }while(std::next_permutation(array, array+7));
        }// k==3

        else if(k==2){
            do{
                iPerm += 1;
                if(iPerm%120 == 1){
                    if(array[0] < array[1]){
                        Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
                        TwoPos(Narray[Int_t(array[0])], Narray[Int_t(array[1])])*
                        QGapPos(Narray[Int_t(array[2])]+Narray[Int_t(array[3])]+Narray[Int_t(array[4])]
                          +Narray[Int_t(array[5])]+Narray[Int_t(array[6])]+n8, 8-k);
                    }
                }
            }while(std::next_permutation(array, array+7));
        }// k==2

        else if(k == 1){
            Correlation = Correlation
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*QGapPos(n1, 1)*QGapPos(n2+n3+n4+n5+n6+n7+n8, 8-k)
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*QGapPos(n2, 1)*QGapPos(n1+n3+n4+n5+n6+n7+n8, 8-k)
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*QGapPos(n3, 1)*QGapPos(n1+n2+n4+n5+n6+n7+n8, 8-k)
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*QGapPos(n4, 1)*QGapPos(n1+n2+n3+n5+n6+n7+n8, 8-k)
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*QGapPos(n5, 1)*QGapPos(n1+n2+n3+n4+n6+n7+n8, 8-k)
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*QGapPos(n6, 1)*QGapPos(n1+n2+n3+n4+n5+n7+n8, 8-k)
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*QGapPos(n7, 1)*QGapPos(n1+n2+n3+n4+n5+n6+n8, 8-k);
        }// k==1

        else if(k == 0){
            Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*QGapPos(n1+n2+n3+n4+n5+n6+n7+n8, 8-k);
        }// k==0

        else{
            printf("Invalid range of k in Eight() \n");
            return {0,0};
        }

    }// loop over k

    return Correlation;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::EightNeg(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t n5, const Int_t n6, const Int_t n7, const Int_t n8) const
{
  TComplex Correlation = {0, 0};
    Int_t Narray[] = {n1, n2, n3, n4, n5, n6, n7};

    for(Int_t k=8; k-->0; )
    {// backward loop of k from m-1 until 0, where m is the m-particle correlation, in this case m=4

        Int_t array[7] = {0,1,2,3,4,5,6};
        Int_t iPerm = 0;
        Int_t count = 0;

        // k==7: there is just one combination, we can add it manually
        if(k==7){
            Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
            SevenNeg(n1, n2, n3, n4, n5, n6, n7)*QGapNeg(n8, 8-k);
        }// k==7

        else if(k==6){
            do{
                iPerm += 1;
                if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4] && array[4] < array[5]){
                    count += 1;
                    Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
                    SixNeg(Narray[Int_t(array[0])], Narray[Int_t(array[1])], Narray[Int_t(array[2])],
                        Narray[Int_t(array[3])], Narray[Int_t(array[4])], Narray[Int_t(array[5])])*
                    QGapNeg(Narray[Int_t(array[6])]+n8, 8-k);
                }
            }while(std::next_permutation(array, array+7));
        }// k==6

        else if(k==5){
            do{
                iPerm += 1;
                if(iPerm%2 == 1){
                    if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4]){
                        Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
                        FiveNeg(Narray[Int_t(array[0])], Narray[Int_t(array[1])], Narray[Int_t(array[2])],
                             Narray[Int_t(array[3])], Narray[Int_t(array[4])])*
                        QGapNeg(Narray[Int_t(array[5])]+Narray[Int_t(array[6])]+n8, 8-k);
                    }
                }
            }while(std::next_permutation(array, array+7));
        }// k==5

        else if(k==4){
            do{
                iPerm += 1;
                if(iPerm%6 == 1){
                    if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3]){
                        Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
                        FourNeg(Narray[Int_t(array[0])], Narray[Int_t(array[1])], Narray[Int_t(array[2])], Narray[Int_t(array[3])])*
                        QGapNeg(Narray[Int_t(array[4])]+Narray[Int_t(array[5])]+Narray[Int_t(array[6])]+n8, 8-k);
                    }
                }
            }while(std::next_permutation(array, array+7));
        }// k==4

        else if(k==3){
            do{
                iPerm += 1;
                if(iPerm%24 == 1){
                    if(array[0] < array[1] && array[1] < array[2]){
                        Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
                        ThreeNeg(Narray[Int_t(array[0])], Narray[Int_t(array[1])], Narray[Int_t(array[2])])*
                        QGapNeg(Narray[Int_t(array[3])]+Narray[Int_t(array[4])]+Narray[Int_t(array[5])]+Narray[Int_t(array[6])]+n8, 8-k);
                    }
                }
            }while(std::next_permutation(array, array+7));
        }// k==3

        else if(k==2){
            do{
                iPerm += 1;
                if(iPerm%120 == 1){
                    if(array[0] < array[1]){
                        Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
                        TwoNeg(Narray[Int_t(array[0])], Narray[Int_t(array[1])])*
                        QGapNeg(Narray[Int_t(array[2])]+Narray[Int_t(array[3])]+Narray[Int_t(array[4])]
                          +Narray[Int_t(array[5])]+Narray[Int_t(array[6])]+n8, 8-k);
                    }
                }
            }while(std::next_permutation(array, array+7));
        }// k==2

        else if(k == 1){
            Correlation = Correlation
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*QGapNeg(n1, 1)*QGapNeg(n2+n3+n4+n5+n6+n7+n8, 8-k)
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*QGapNeg(n2, 1)*QGapNeg(n1+n3+n4+n5+n6+n7+n8, 8-k)
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*QGapNeg(n3, 1)*QGapNeg(n1+n2+n4+n5+n6+n7+n8, 8-k)
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*QGapNeg(n4, 1)*QGapNeg(n1+n2+n3+n5+n6+n7+n8, 8-k)
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*QGapNeg(n5, 1)*QGapNeg(n1+n2+n3+n4+n6+n7+n8, 8-k)
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*QGapNeg(n6, 1)*QGapNeg(n1+n2+n3+n4+n5+n7+n8, 8-k)
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*QGapNeg(n7, 1)*QGapNeg(n1+n2+n3+n4+n5+n6+n8, 8-k);
        }// k==1

        else if(k == 0){
            Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*QGapNeg(n1+n2+n3+n4+n5+n6+n7+n8, 8-k);
        }// k==0

        else{
            printf("Invalid range of k in Eight() \n");
            return {0,0};
        }

    }// loop over k

    return Correlation;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::FourGap(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4) const
{
  TComplex formula = QGapPos(n1,1)*QGapPos(n2,1)*QGapNeg(n3,1)*QGapNeg(n4,1)-QGapPos(n1+n2,2)*QGapNeg(n3,1)*QGapNeg(n4,1)
                    -QGapPos(n1,1)*QGapPos(n2,1)*QGapNeg(n3+n4,2)+QGapPos(n1+n2,2)*QGapNeg(n3+n4,2);
  //same as
  //TComplex formula = TwoPos(n1,n2)*TwoNeg(n3,n4);
	return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::SixGap(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t n5, const Int_t n6) const
{
  TComplex formula = ThreePos(n1,n2,n3)*ThreeNeg(n4,n5,n6);
	return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::EightGap(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t n5, const Int_t n6, const Int_t n7, const Int_t n8) const
{
  TComplex formula = FourPos(n1,n2,n3,n4)*FourNeg(n5,n6,n7,n8);
	return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::Four3sub(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4) const
{
  // left = neg, middle = mid; rigth = pos
  TComplex formula = QGapMid(n1,1)*QGapMid(n2,1)*QGapNeg(n3,1)*QGapPos(n4,1)-QGapMid(n1+n2,2)*QGapNeg(n3,1)*QGapPos(n4,1);
  return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::ThreeDiff(const Int_t n1, const Int_t n2, const Int_t n3) const
{
  TComplex formula = P(n1,1)*Q(n2,1)*Q(n3,1)-S(n1+n2,2)*Q(n3,1)-S(n1+n3,2)*Q(n2,1)
 		                 - P(n1,1)*Q(n2+n3,2)+2.0*S(n1+n2+n3,3);
  return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::ThreeDiffGapPos(const Int_t n1, const Int_t n2, const Int_t n3) const
{
  TComplex formula = PGapPos(n1,1)*QGapNeg(n2,1)*QGapNeg(n3,1) - PGapPos(n1,1)*QGapNeg(n2+n3,2);
  return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::ThreeDiffGapNeg(const Int_t n1, const Int_t n2, const Int_t n3) const
{
  TComplex formula = PGapNeg(n1,1)*QGapPos(n2,1)*QGapPos(n3,1) - PGapNeg(n1,1)*QGapPos(n2+n3,2);
  return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::ThreeDiffPos(const Int_t n1, const Int_t n2, const Int_t n3) const
{
  TComplex formula = PGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n3,1)-SGapPos(n1+n2,2)*QGapPos(n3,1)-SGapPos(n1+n3,2)*QGapPos(n2,1)
 		                 - PGapPos(n1,1)*QGapPos(n2+n3,2)+2.0*SGapPos(n1+n2+n3,3);
  return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::ThreeDiffNeg(const Int_t n1, const Int_t n2, const Int_t n3) const
{
  TComplex formula = PGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n3,1)-SGapNeg(n1+n2,2)*QGapNeg(n3,1)-SGapNeg(n1+n3,2)*QGapNeg(n2,1)
 		                 - PGapNeg(n1,1)*QGapNeg(n2+n3,2)+2.0*SGapNeg(n1+n2+n3,3);
  return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::FourDiff(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4) const
{
  TComplex formula = P(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-S(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*S(n1+n3,2)*Q(n4,1)
                    - P(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.0*S(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*S(n1+n4,2)
                    + Q(n2+n3,2)*S(n1+n4,2)-P(n1,1)*Q(n3,1)*Q(n2+n4,2)+S(n1+n3,2)*Q(n2+n4,2)
                    + 2.0*Q(n3,1)*S(n1+n2+n4,3)-P(n1,1)*Q(n2,1)*Q(n3+n4,2)+S(n1+n2,2)*Q(n3+n4,2)
                    + 2.0*Q(n2,1)*S(n1+n3+n4,3)+2.0*P(n1,1)*Q(n2+n3+n4,3)-6.0*S(n1+n2+n3+n4,4);
  return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::FourDiffPos(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4) const
{
  TComplex formula = PGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n3,1)*QGapPos(n4,1)-SGapPos(n1+n2,2)*QGapPos(n3,1)*QGapPos(n4,1)-QGapPos(n2,1)*SGapPos(n1+n3,2)*QGapPos(n4,1)
                    - PGapPos(n1,1)*QGapPos(n2+n3,2)*QGapPos(n4,1)+2.0*SGapPos(n1+n2+n3,3)*QGapPos(n4,1)-QGapPos(n2,1)*QGapPos(n3,1)*SGapPos(n1+n4,2)
                    + QGapPos(n2+n3,2)*SGapPos(n1+n4,2)-PGapPos(n1,1)*QGapPos(n3,1)*QGapPos(n2+n4,2)+SGapPos(n1+n3,2)*QGapPos(n2+n4,2)
                    + 2.0*QGapPos(n3,1)*SGapPos(n1+n2+n4,3)-PGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n3+n4,2)+SGapPos(n1+n2,2)*QGapPos(n3+n4,2)
                    + 2.0*QGapPos(n2,1)*SGapPos(n1+n3+n4,3)+2.0*PGapPos(n1,1)*QGapPos(n2+n3+n4,3)-6.0*SGapPos(n1+n2+n3+n4,4);
  return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::FourDiffNeg(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4) const
{
  TComplex formula = PGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n3,1)*QGapNeg(n4,1)-SGapNeg(n1+n2,2)*QGapNeg(n3,1)*QGapNeg(n4,1)-QGapNeg(n2,1)*SGapNeg(n1+n3,2)*QGapNeg(n4,1)
                    - PGapNeg(n1,1)*QGapNeg(n2+n3,2)*QGapNeg(n4,1)+2.0*SGapNeg(n1+n2+n3,3)*QGapNeg(n4,1)-QGapNeg(n2,1)*QGapNeg(n3,1)*SGapNeg(n1+n4,2)
                    + QGapNeg(n2+n3,2)*SGapNeg(n1+n4,2)-PGapNeg(n1,1)*QGapNeg(n3,1)*QGapNeg(n2+n4,2)+SGapNeg(n1+n3,2)*QGapNeg(n2+n4,2)
                    + 2.0*QGapNeg(n3,1)*SGapNeg(n1+n2+n4,3)-PGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n3+n4,2)+SGapNeg(n1+n2,2)*QGapNeg(n3+n4,2)
                    + 2.0*QGapNeg(n2,1)*SGapNeg(n1+n3+n4,3)+2.0*PGapNeg(n1,1)*QGapNeg(n2+n3+n4,3)-6.0*SGapNeg(n1+n2+n3+n4,4);
  return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::FourDiffGapPos(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4) const
{
  TComplex formula = PGapPos(n1,1)*QGapPos(n2,1)*QGapNeg(n3,1)*QGapNeg(n4,1)
                      - SGapPos(n1+n2,2)*QGapNeg(n3,1)*QGapNeg(n4,1)
                      - PGapPos(n1,1)*QGapPos(n2,1)*QGapNeg(n3+n4,2)
                      + SGapPos(n1+n2,2)*QGapNeg(n3+n4,2);
  //same as
  //TComplex formula = TwoDiffPos(n1,n2)*TwoNeg(n3,n4);
  return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::FourDiffGapNeg(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4) const
{
  TComplex formula = PGapNeg(n1,1)*QGapNeg(n2,1)*QGapPos(n3,1)*QGapPos(n4,1)
                      - SGapNeg(n1+n2,2)*QGapPos(n3,1)*QGapPos(n4,1)
                      - PGapNeg(n1,1)*QGapNeg(n2,1)*QGapPos(n3+n4,2)
                      + SGapNeg(n1+n2,2)*QGapPos(n3+n4,2);
  // same as
  // TComplex formula = TwoDiffNeg(n1,n2)*TwoPos(n3,n4);
  return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::SixDiffGapPos(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t n5, const Int_t n6) const
{
  TComplex formula = ThreeDiffPos(n1,n2,n3)*ThreeNeg(n4,n5,n6);
  return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::SixDiffGapNeg(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t n5, const Int_t n6) const
{
  TComplex formula = ThreeDiffNeg(n1,n2,n3)*ThreePos(n4,n5,n6);
  return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::EightDiffGapPos(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t n5, const Int_t n6, const Int_t n7, const Int_t n8) const
{
  TComplex formula = FourDiffPos(n1,n2,n3,n4)*FourNeg(n5,n6,n7,n8);
  return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::EightDiffGapNeg(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t n5, const Int_t n6, const Int_t n7, const Int_t n8) const
{
  TComplex formula = FourDiffNeg(n1,n2,n3,n4)*FourPos(n5,n6,n7,n8);
  return formula;
}
// ============================================================================
std::vector<Double_t> AliAnalysisTaskUniFlow::MakeBinsVector(Int_t num, Double_t min, Double_t max)
{
    std::vector<Double_t> vec = std::vector<Double_t>();
    Double_t step = (max - min) / num;
    Double_t edge = min;

    for(Int_t i(0); i < num+1; ++i) {
        vec.push_back(edge);
        // printf("%d: %f\n",i,edge);
        edge += step;
    }
    // printf("num %d | size %lu\n",num,vec.size() );

    return vec;
}
// ============================================================================
void AliAnalysisTaskUniFlow::UserCreateOutputObjects()
{
  // create output objects
  // this function is called ONCE at the start of your analysis (RUNTIME)
  // *************************************************************

  DumpTObjTable("UserCreateOutputObjects: start");

  // task initialization
  fInit = InitializeTask();
  if(!fInit) { return; }

  DumpTObjTable("UserCreateOutputObjects: after Initialization");


  // list all parameters used in this analysis
  ListParameters();

  // creating output lists
  for(Int_t iSpec(0); iSpec < kUnknown; ++iSpec) {
    fListFlow[iSpec] = new TList();
    fListFlow[iSpec]->SetOwner(kTRUE);
    fListFlow[iSpec]->SetName(Form("fFlow%s",GetSpeciesName(PartSpecies(iSpec))));
  }

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
  if(fMC) {
      fListMC = new TList();
      fListMC->SetOwner(kTRUE);
  }

  // setting number of bins based on set range with fixed width
  const Int_t iFlowRFPsPtBinNum = (Int_t) ((fFlowRFPsPtMax - fFlowRFPsPtMin) / 0.1 + 0.5);

  // creating output correlations profiles based on CorrTasks
  Int_t iNumTasks = fVecCorrTask.size();
  if(fRunMode != kSkipFlow && iNumTasks > 0)
  {
    for(Int_t iTask(0); iTask < iNumTasks; ++iTask)
    {
      AliUniFlowCorrTask* task = fVecCorrTask.at(iTask);
      if(!task) { fInit = kFALSE; AliError(Form("AliUniFlowCorrTask %d does not exists\n",iTask)); return; }

      if(fFlowNumHarmonicsMax < task->fMaxHarm) { fInit = kFALSE; AliError(Form("Max Harm error in task %d\n",iTask)); return; }
      if(fFlowNumWeightPowersMax < task->fMaxWeightPower) { fInit = kFALSE; AliError(Form("Max Weight Power error in task %d\n",iTask));return; }

      Bool_t bHasGap = task->HasGap();
      const char* corName = task->fsName.Data();
      const char* corLabel = task->fsLabel.Data();

      for(Int_t iSpec(0); iSpec < kUnknown; ++iSpec)
      {
        // check if AliUniFlowCorrTask should be done for all flow particles (RFP/POI/Both)
        if(!task->fbDoRefs && iSpec == kRefs) { continue; }
        if(!task->fbDoPOIs && iSpec != kRefs) { continue; }

        if(!fProcessSpec[iSpec]) { continue; }
        if(iSpec == kKaon && (!fProcessSpec[kPion] || !fProcessSpec[kProton])) { continue; }

        Int_t iNumPtFixBins = fFlowPOIsPtBinEdges[iSpec].size() - 1;
        Double_t* dPtFixBinEdges = fFlowPOIsPtBinEdges[iSpec].data();

        for(Int_t iSample(0); iSample < fNumSamples; ++iSample)
        {
          if(iSample > 0 && !fSampling) { break; }
          if(iSample > 0 && HasMass(PartSpecies(iSpec))) { break; } // reconstructed are not sampled

          TH1* profile = nullptr;
          TH1* profileNeg = nullptr;

          switch(iSpec)
          {
            case kRefs :
            {
              profile = new TProfile(Form("%s_Pos_sample%d",corName,iSample), Form("%s: %s; %s",GetSpeciesLabel(PartSpecies(iSpec)),corLabel,GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax);
              if(fEtaCheckRFP) profileNeg = new TProfile(Form("%s_Neg_sample%d",corName,iSample), Form("%s: %s; %s",GetSpeciesLabel(PartSpecies(iSpec)),corLabel,GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax);
              break;
            }

            case kCharged :
            case kPion :
            case kKaon :
            case kProton :
            {
                if(iNumPtFixBins > 0) {
                    profile = new TProfile2D(Form("%s_Pos_sample%d",corName,iSample), Form("%s: %s (Pos); %s; #it{p}_{T} (GeV/#it{c})",GetSpeciesLabel(PartSpecies(iSpec)), corLabel, GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax, iNumPtFixBins,dPtFixBinEdges,"");
                    if(bHasGap) { profileNeg = new TProfile2D(Form("%s_Neg_sample%d",corName,iSample), Form("%s: %s (Neg); %s; #it{p}_{T} (GeV/#it{c})",GetSpeciesLabel(PartSpecies(iSpec)), corLabel,GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax, iNumPtFixBins,dPtFixBinEdges,""); }
                } else {
                    profile = new TProfile2D(Form("%s_Pos_sample%d",corName,iSample), Form("%s: %s (Pos); %s; #it{p}_{T} (GeV/#it{c})",GetSpeciesLabel(PartSpecies(iSpec)), corLabel,GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax, fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax);
                    if(bHasGap) { profileNeg = new TProfile2D(Form("%s_Neg_sample%d",corName,iSample), Form("%s: %s (Neg); %s; #it{p}_{T} (GeV/#it{c})",GetSpeciesLabel(PartSpecies(iSpec)), corLabel,GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax, fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax); }
                }
                break;
            }

            case kK0s:
            {
                if(iNumPtFixBins > 0) {
                    profile = new TProfile3D(Form("%s_Pos_sample%d",corName,iSample), Form("%s: %s (Pos); %s; #it{p}_{T} (GeV/#it{c}); #it{m}_{inv} (GeV/#it{c}^{2})",GetSpeciesLabel(PartSpecies(iSpec)),corLabel,GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,MakeBinsVector(fCentBinNum,fCentMin,fCentMax).data(), iNumPtFixBins,dPtFixBinEdges, fV0sNumBinsMass,MakeBinsVector(fV0sNumBinsMass,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax).data(),"");
                    if(bHasGap) { profileNeg = new TProfile3D(Form("%s_Neg_sample%d",corName,iSample), Form("%s: %s (Neg); %s; #it{p}_{T} (GeV/#it{c}); #it{m}_{inv} (GeV/#it{c}^{2})",GetSpeciesLabel(PartSpecies(iSpec)),corLabel,GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,MakeBinsVector(fCentBinNum,fCentMin,fCentMax).data(), iNumPtFixBins,dPtFixBinEdges, fV0sNumBinsMass,MakeBinsVector(fV0sNumBinsMass,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax).data(),""); }
                } else {
                    profile = new TProfile3D(Form("%s_Pos_sample%d",corName,iSample), Form("%s: %s (Pos); %s; #it{p}_{T} (GeV/#it{c}); #it{m}_{inv} (GeV/#it{c}^{2})",GetSpeciesLabel(PartSpecies(iSpec)),corLabel,GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax, fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax, fV0sNumBinsMass,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax);
                    if(bHasGap) { profileNeg = new TProfile3D(Form("%s_Neg_sample%d",corName,iSample), Form("%s: %s (Neg); %s; #it{p}_{T} (GeV/#it{c}); #it{m}_{inv} (GeV/#it{c}^{2})",GetSpeciesLabel(PartSpecies(iSpec)),corLabel,GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax, fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax, fV0sNumBinsMass,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax); }
                }
                break;
            }

            case kLambda:
            {
                if(iNumPtFixBins > 0) {
                    profile = new TProfile3D(Form("%s_Pos_sample%d",corName,iSample), Form("%s: %s (Pos); %s; #it{p}_{T} (GeV/#it{c}); #it{m}_{inv} (GeV/#it{c}^{2})",GetSpeciesLabel(PartSpecies(iSpec)),corLabel,GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,MakeBinsVector(fCentBinNum,fCentMin,fCentMax).data(), iNumPtFixBins,dPtFixBinEdges, fV0sNumBinsMass,MakeBinsVector(fV0sNumBinsMass,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax).data());
                    if(bHasGap) { profileNeg = new TProfile3D(Form("%s_Neg_sample%d",corName,iSample), Form("%s: %s (Neg); %s; #it{p}_{T} (GeV/#it{c}); #it{m}_{inv} (GeV/#it{c}^{2})",GetSpeciesLabel(PartSpecies(iSpec)),corLabel,GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,MakeBinsVector(fCentBinNum,fCentMin,fCentMax).data(), iNumPtFixBins,dPtFixBinEdges, fV0sNumBinsMass,MakeBinsVector(fV0sNumBinsMass,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax).data()); }
                } else {
                    profile = new TProfile3D(Form("%s_Pos_sample%d",corName,iSample), Form("%s: %s (Pos); %s; #it{p}_{T} (GeV/#it{c}); #it{m}_{inv} (GeV/#it{c}^{2})",GetSpeciesLabel(PartSpecies(iSpec)),corLabel,GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax, fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax, fV0sNumBinsMass,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax);
                    if(bHasGap) { profileNeg = new TProfile3D(Form("%s_Neg_sample%d",corName,iSample), Form("%s: %s (Neg); %s; #it{p}_{T} (GeV/#it{c}); #it{m}_{inv} (GeV/#it{c}^{2})",GetSpeciesLabel(PartSpecies(iSpec)),corLabel,GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax, fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax, fV0sNumBinsMass,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax); }
                }
                break;
            }

            case kPhi:
            {
                if(iNumPtFixBins > 0) {
                    profile = new TProfile3D(Form("%s_Pos_sample%d",corName,iSample), Form("%s: %s (Pos); %s; #it{p}_{T} (GeV/#it{c}); #it{m}_{inv} (GeV/#it{c}^{2})",GetSpeciesLabel(PartSpecies(iSpec)),corLabel,GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,MakeBinsVector(fCentBinNum,fCentMin,fCentMax).data(), iNumPtFixBins,dPtFixBinEdges, fPhiNumBinsMass,MakeBinsVector(fPhiNumBinsMass,fCutPhiInvMassMin,fCutPhiInvMassMax).data());
                    if(bHasGap) { profileNeg = new TProfile3D(Form("%s_Neg_sample%d",corName,iSample), Form("%s: %s (Neg); %s; #it{p}_{T} (GeV/#it{c}); #it{m}_{inv} (GeV/#it{c}^{2})",GetSpeciesLabel(PartSpecies(iSpec)),corLabel,GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,MakeBinsVector(fCentBinNum,fCentMin,fCentMax).data(), iNumPtFixBins,dPtFixBinEdges, fPhiNumBinsMass,MakeBinsVector(fPhiNumBinsMass,fCutPhiInvMassMin,fCutPhiInvMassMax).data()); }
                } else {
                    profile = new TProfile3D(Form("%s_Pos_sample%d",corName,iSample), Form("%s: %s (Pos); %s; #it{p}_{T} (GeV/#it{c}); #it{m}_{inv} (GeV/#it{c}^{2})",GetSpeciesLabel(PartSpecies(iSpec)),corLabel,GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax, fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax, fPhiNumBinsMass,fCutPhiInvMassMin,fCutPhiInvMassMax);
                    if(bHasGap) { profileNeg = new TProfile3D(Form("%s_Neg_sample%d",corName,iSample), Form("%s: %s (Neg); %s; #it{p}_{T} (GeV/#it{c}); #it{m}_{inv} (GeV/#it{c}^{2})",GetSpeciesLabel(PartSpecies(iSpec)),corLabel,GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax, fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax, fPhiNumBinsMass,fCutPhiInvMassMin,fCutPhiInvMassMax); }
                }
                break;
            }
          }

          if(!profile) { fInit = kFALSE; AliError("Profile (Pos) NOT created!"); task->Print(); return; }

          // check if same profile does not exists already
          if(fListFlow[iSpec]->FindObject(profile->GetName())) {
            AliError(Form("AliUniFlowCorrTask %d : Profile '%s' already exists! Please check run macro for AliUniFlowCorrTask duplicates!",iTask,profile->GetName()));
            fInit = kFALSE;
            task->Print();
            delete profile;
            return;
          }

          profile->Sumw2();
          fListFlow[iSpec]->Add(profile);

          if(bHasGap && (iSpec != kRefs || fEtaCheckRFP))
          { // Refs does not distinquish Pos/Neg
            if(!profileNeg) { fInit = kFALSE; AliError("Profile (Neg) NOT created!"); task->Print(); return; }

            // same for Neg
            if(fListFlow[iSpec]->FindObject(profileNeg->GetName())) {
              AliError(Form("AliUniFlowCorrTask %d : Profile '%s' already exists! Please check run macro for AliUniFlowCorrTask duplicates!",iTask,profile->GetName()));
              fInit = kFALSE;
              task->Print();
              delete profileNeg;
              return;
            }

            profileNeg->Sumw2();
            fListFlow[iSpec]->Add(profileNeg);
          }
        } // end-for {iSample}
      } // end-for {iSpec}
    } // end-for {iTask}

    // Making THnSparse distribution of candidates
    // species independent

    TString sLabelCand[SparseCand::kDim];
    sLabelCand[SparseCand::kInvMass] = "#it{m}_{inv} (GeV/#it{c}^{2})";
    sLabelCand[SparseCand::kCent] = GetCentEstimatorLabel(fCentEstimator);
    sLabelCand[SparseCand::kPt] = "#it{p}_{T} (GeV/c)";
    sLabelCand[SparseCand::kEta] = "#eta";
    TString sAxes = TString(); for(Int_t i(0); i < SparseCand::kDim; ++i) { sAxes += Form("%s; ",sLabelCand[i].Data()); }

    Int_t iNumBinsCand[SparseCand::kDim]; Double_t dMinCand[SparseCand::kDim]; Double_t dMaxCand[SparseCand::kDim];
    iNumBinsCand[SparseCand::kCent] = fCentBinNum; dMinCand[SparseCand::kCent] = fCentMin; dMaxCand[SparseCand::kCent] = fCentMax;
    iNumBinsCand[SparseCand::kEta] = fFlowEtaBinNum; dMinCand[SparseCand::kEta] = -fFlowEtaMax; dMaxCand[SparseCand::kEta] = fFlowEtaMax;

    // species dependent
    if(fProcessSpec[kK0s] || fProcessSpec[kLambda])
    {
        if(fFlowPOIsPtBinEdges[kK0s].size() > 0) {
            iNumBinsCand[SparseCand::kPt] = (fFlowPOIsPtBinEdges[kK0s].size() - 1);
            dMinCand[SparseCand::kPt] = fFlowPOIsPtBinEdges[kK0s].front();
            dMaxCand[SparseCand::kPt] = fFlowPOIsPtBinEdges[kK0s].back();
        } else {
            iNumBinsCand[SparseCand::kPt] = fFlowPOIsPtBinNum;
            dMinCand[SparseCand::kPt] = fFlowPOIsPtMin;
            dMaxCand[SparseCand::kPt] = fFlowPOIsPtMax;
        }

        iNumBinsCand[SparseCand::kInvMass] = fV0sNumBinsMass; dMinCand[SparseCand::kInvMass] = fCutV0sInvMassK0sMin; dMaxCand[SparseCand::kInvMass] = fCutV0sInvMassK0sMax;
        fhsCandK0s = new THnSparseD("fhsCandK0s",Form("K_{S}^{0}: Distribution; %s;", sAxes.Data()), SparseCand::kDim, iNumBinsCand, dMinCand, dMaxCand);
        fhsCandK0s->Sumw2();
        fListFlow[kK0s]->Add(fhsCandK0s);

        if(fFlowPOIsPtBinEdges[kLambda].size() > 0) {
            iNumBinsCand[SparseCand::kPt] = (fFlowPOIsPtBinEdges[kLambda].size() - 1);
            dMinCand[SparseCand::kPt] = fFlowPOIsPtBinEdges[kLambda].front();
            dMaxCand[SparseCand::kPt] = fFlowPOIsPtBinEdges[kLambda].back();
        } else {
            iNumBinsCand[SparseCand::kPt] = fFlowPOIsPtBinNum;
            dMinCand[SparseCand::kPt] = fFlowPOIsPtMin;
            dMaxCand[SparseCand::kPt] = fFlowPOIsPtMax;
        }

        iNumBinsCand[SparseCand::kInvMass] = fV0sNumBinsMass; dMinCand[SparseCand::kInvMass] = fCutV0sInvMassLambdaMin; dMaxCand[SparseCand::kInvMass] = fCutV0sInvMassLambdaMax;
        fhsCandLambda = new THnSparseD("fhsCandLambda",Form("#Lambda: Distribution; %s;", sAxes.Data()), SparseCand::kDim, iNumBinsCand, dMinCand, dMaxCand);
        fhsCandLambda->Sumw2();
        fListFlow[kLambda]->Add(fhsCandLambda);
    }

    if(fProcessSpec[kPhi])
    {
        if(fFlowPOIsPtBinEdges[kPhi].size() > 0) {
            iNumBinsCand[SparseCand::kPt] = (fFlowPOIsPtBinEdges[kPhi].size() - 1);
            dMinCand[SparseCand::kPt] = fFlowPOIsPtBinEdges[kPhi].front();
            dMaxCand[SparseCand::kPt] = fFlowPOIsPtBinEdges[kPhi].back();
        } else {
            iNumBinsCand[SparseCand::kPt] = fFlowPOIsPtBinNum;
            dMinCand[SparseCand::kPt] = fFlowPOIsPtMin;
            dMaxCand[SparseCand::kPt] = fFlowPOIsPtMax;
        }

        iNumBinsCand[SparseCand::kInvMass] = fPhiNumBinsMass; dMinCand[SparseCand::kInvMass] = fCutPhiInvMassMin; dMaxCand[SparseCand::kInvMass] = fCutPhiInvMassMax;
        fhsCandPhi = new THnSparseD("fhsCandPhi",Form("#phi (Sig): Distribution; %s;", sAxes.Data()), SparseCand::kDim, iNumBinsCand, dMinCand, dMaxCand);
        fhsCandPhi->Sumw2();
        fListFlow[kPhi]->Add(fhsCandPhi);

        fhsCandPhiBg = new THnSparseD("fhsCandPhiBg",Form("#phi (Bg): Distribution; %s;", sAxes.Data()), SparseCand::kDim, iNumBinsCand, dMinCand, dMaxCand);
        fhsCandPhiBg->Sumw2();
        fListFlow[kPhi]->Add(fhsCandPhiBg);
    }
  } // end-if {fRunMode != fSkipFlow || iNumCorrTask > 0 }

  // creating GF weights
  if(fFlowFillWeights || fFlowUseWeights)
  {
    for(Int_t iSpec(0); iSpec < kUnknown; ++iSpec)
    {
      if(!fProcessSpec[iSpec]) { continue; }
      if(iSpec == kKaon && (!fProcessSpec[kPion] || !fProcessSpec[kProton])) { continue; }

      if(fFlowFillWeights) {
        if(fFlowUse3Dweights) {
          const char* weightName = Form("fh3Weights%s",GetSpeciesName(PartSpecies(iSpec)));
          const char* weightLabel = Form("Weights: %s; #varphi; #eta; PV-z (cm)", GetSpeciesName(PartSpecies(iSpec)));
          fh3Weights[iSpec] = new TH3D(weightName, weightLabel, fFlowPhiBinNum,0.0,TMath::TwoPi(), fFlowEtaBinNum,-fFlowEtaMax,fFlowEtaMax, 2*fPVtxCutZ,-fPVtxCutZ,fPVtxCutZ);
          fh3Weights[iSpec]->Sumw2();
          fFlowWeights->Add(fh3Weights[iSpec]);
        } else {
          const char* weightName = Form("fh2Weights%s",GetSpeciesName(PartSpecies(iSpec)));
          const char* weightLabel = Form("Weights: %s; #varphi; #eta", GetSpeciesName(PartSpecies(iSpec)));
          fh2Weights[iSpec] = new TH2D(weightName, weightLabel, fFlowPhiBinNum,0.0,TMath::TwoPi(), fFlowEtaBinNum,-fFlowEtaMax,fFlowEtaMax);
          fh2Weights[iSpec]->Sumw2();
          fFlowWeights->Add(fh2Weights[iSpec]);
        }
      }

      if(fFlowUseWeights && fFlowFillAfterWeights)
      {
        if(fFlowUse3Dweights) {
          const char* weightName = Form("fh3AfterWeights%s",GetSpeciesName(PartSpecies(iSpec)));
          const char* weightLabel = Form("Weights (after): %s; #varphi; #eta; PV-z (cm)",GetSpeciesLabel(PartSpecies(iSpec)));
          fh3AfterWeights[iSpec] = new TH3D(weightName,weightLabel, fFlowPhiBinNum,0.0,TMath::TwoPi(), fFlowEtaBinNum,-fFlowEtaMax,fFlowEtaMax, 2*fPVtxCutZ,-fPVtxCutZ,fPVtxCutZ);
          fh3AfterWeights[iSpec]->Sumw2();
          fFlowWeights->Add(fh3AfterWeights[iSpec]);
        } else {
          const char* weightName = Form("fh2AfterWeights%s",GetSpeciesName(PartSpecies(iSpec)));
          const char* weightLabel = Form("Weights (after): %s; #varphi; #eta;",GetSpeciesLabel(PartSpecies(iSpec)));
          fh2AfterWeights[iSpec] = new TH2D(weightName,weightLabel, fFlowPhiBinNum,0.0,TMath::TwoPi(), fFlowEtaBinNum,-fFlowEtaMax,fFlowEtaMax);
          fh2AfterWeights[iSpec]->Sumw2();
          fFlowWeights->Add(fh2AfterWeights[iSpec]);
        }
      }
    } // end-for {iSpec}
  }

  // Selection / reconstruction counters : omni-present
  {
    TString sEventCounterLabel[] = {"Input","Physics selection OK","EventCuts OK","Event OK","#RPFs OK","Multiplicity OK","Selected"};
    const Int_t iEventCounterBins = sizeof(sEventCounterLabel)/sizeof(sEventCounterLabel[0]);
    fhEventCounter = new TH1D("fhEventCounter","Event Counter",iEventCounterBins,0,iEventCounterBins);
    for(Int_t i(0); i < iEventCounterBins; ++i) { fhEventCounter->GetXaxis()->SetBinLabel(i+1, sEventCounterLabel[i].Data() ); }
    fQAEvents->Add(fhEventCounter);
  }

  {
    TString sChargedCounterLabel[] = {"Input","FB","#TPC-Cls","DCA-z","DCA-xy","Selected","POIs","Refs"};
    const Int_t iNBinsChargedCounter = sizeof(sChargedCounterLabel)/sizeof(sChargedCounterLabel[0]);
    fhChargedCounter = new TH1D("fhChargedCounter","Charged tracks: Counter",iNBinsChargedCounter,0,iNBinsChargedCounter);
    for(Int_t i(0); i < iNBinsChargedCounter; i++) fhChargedCounter->GetXaxis()->SetBinLabel(i+1, sChargedCounterLabel[i].Data() );
    fQACharged->Add(fhChargedCounter);
  }

  if(fProcessSpec[kPion] || fProcessSpec[kKaon] || fProcessSpec[kProton])
  {
    TString sPIDCounterLabel[] = {"Input","Selected","Pion","Kaon","Proton"};
    const Short_t iNBinsPIDCounter = sizeof(sPIDCounterLabel)/sizeof(sPIDCounterLabel[0]);
    fhPIDCounter = new TH1D("fhPIDCounter","PID: Counter",iNBinsPIDCounter,0,iNBinsPIDCounter);
    for(Int_t i(0); i < iNBinsPIDCounter; i++) fhPIDCounter->GetXaxis()->SetBinLabel(i+1, sPIDCounterLabel[i].Data() );
    fQAPID->Add(fhPIDCounter);
  }

  if(fProcessSpec[kPhi])
  {
    TString sPhiCounterLabel[] = {"Input","InvMass","Acceptance","Before charge","Unlike-sign","BG"};
    const Int_t iNBinsPhiCounter = sizeof(sPhiCounterLabel)/sizeof(sPhiCounterLabel[0]);
    fhPhiCounter = new TH1D("fhPhiCounter","#phi: Counter",iNBinsPhiCounter,0,iNBinsPhiCounter);
    for(Int_t i(0); i < iNBinsPhiCounter; ++i) { fhPhiCounter->GetXaxis()->SetBinLabel(i+1, sPhiCounterLabel[i].Data() ); }
    fQAPhi->Add(fhPhiCounter);
  }

  if(fProcessSpec[kK0s] || fProcessSpec[kLambda])
  {
    TString sV0sCounterLabel[] = {"Input","Daughters OK","Mother acceptance","Daughter acceptance","Charge","Reconstruction method","Daughter FB","TPC refit","Kinks","Daughters track quality","DCA to PV","Daughters DCA","Decay radius","Common passed","K^{0}_{S}","#Lambda/#bar{#Lambda}","K^{0}_{S} && #Lambda/#bar{#Lambda}"};
    const Int_t iNBinsV0sCounter = sizeof(sV0sCounterLabel)/sizeof(sV0sCounterLabel[0]);
    fhV0sCounter = new TH1D("fhV0sCounter","V^{0}: Counter",iNBinsV0sCounter,0,iNBinsV0sCounter);
    for(Int_t i(0); i < iNBinsV0sCounter; ++i) { fhV0sCounter->GetXaxis()->SetBinLabel(i+1, sV0sCounterLabel[i].Data() ); }
    fQAV0s->Add(fhV0sCounter);

    TString sV0sK0sCounterLabel[] = {"Input","#it{y}","InvMass","CPA","Armenteros-Podolanski","c#tau","Daughters PID","Competing InvMass","Selected"};
    const Int_t iNBinsV0sK0sCounter = sizeof(sV0sK0sCounterLabel)/sizeof(sV0sK0sCounterLabel[0]);
    fhV0sCounterK0s = new TH1D("fhV0sCounterK0s","V^{0}: K^{0}_{S} Counter",iNBinsV0sK0sCounter,0,iNBinsV0sK0sCounter);
    for(Int_t i(0); i < iNBinsV0sK0sCounter; ++i) { fhV0sCounterK0s->GetXaxis()->SetBinLabel(i+1, sV0sK0sCounterLabel[i].Data() ); }
    fQAV0s->Add(fhV0sCounterK0s);

    TString sV0sLambdaCounterLabel[] = {"Input","#it{y}","InvMass","CPA","Armenteros-Podolanski","c#tau","Daughter PID","Competing InvMass","Selected","only #Lambda","only #bar{#Lambda}","#Lambda && #bar{#Lambda}"};
    const Int_t iNBinsV0sLambdaCounter = sizeof(sV0sLambdaCounterLabel)/sizeof(sV0sLambdaCounterLabel[0]);
    fhV0sCounterLambda = new TH1D("fhV0sCounterLambda","V^{0}: #Lambda/#bar{#Lambda} Counter",iNBinsV0sLambdaCounter,0,iNBinsV0sLambdaCounter);
    for(Int_t i(0); i < iNBinsV0sLambdaCounter; ++i) { fhV0sCounterLambda->GetXaxis()->SetBinLabel(i+1, sV0sLambdaCounterLabel[i].Data() ); }
    fQAV0s->Add(fhV0sCounterLambda);
  }

  // #### Fill QA[2] plots
  if(fFillQA)
  {
    const Int_t iNBinsPIDstatus = 4;
    TString sPIDstatus[iNBinsPIDstatus] = {"kDetNoSignal","kDetPidOk","kDetMismatch","kDetNoParams"};

    // event histogram
    fEventCuts.AddQAplotsToList(fQAEvents);

    fhEventSampling = new TH2D("fhEventSampling",Form("Event sampling; %s; sample index", GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax, fNumSamples,0,fNumSamples);
    fQAEvents->Add(fhEventSampling);
    fhEventCentrality = new TH1D("fhEventCentrality",Form("Event centrality (%s); %s", GetCentEstimatorLabel(fCentEstimator), GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax);
    fQAEvents->Add(fhEventCentrality);

    Int_t iMinRFPs = 0;
    Int_t iMaxRFPs = 150;
    Int_t iBinsRPFs = 150;
    if(fColSystem == kPbPb) {
        iMinRFPs = 0;
        iMaxRFPs = 5000;
        iBinsRPFs = 100;
    }
    fh2EventCentralityNumRefs = new TH2D("fh2EventCentralityNumRefs",Form("Event centrality (%s) vs. N_{RFP}; %s; N_{RFP}",GetCentEstimatorLabel(fCentEstimator), GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax, iBinsRPFs,iMinRFPs,iMaxRFPs);
    fQAEvents->Add(fh2EventCentralityNumRefs);

    if(fEventRejectAddPileUp)
    {
      fhQAEventsfMult32vsCentr = new TH2D("fhQAEventsfMult32vsCentr", "; centrality V0M; TPC multiplicity (FB32)", 100, 0, 100, 100, 0, 3000);
      fQAEvents->Add(fhQAEventsfMult32vsCentr);
      fhQAEventsMult128vsCentr = new TH2D("fhQAEventsfMult128vsCentr", "; centrality V0M; TPC multiplicity (FB128)", 100, 0, 100, 100, 0, 5000);
      fQAEvents->Add(fhQAEventsMult128vsCentr);
      fhQAEventsfMultTPCvsTOF = new TH2D("fhQAEventsfMultTPCvsTOF", "; TPC FB32 multiplicity; TOF multiplicity", 200, 0, 4000, 200, 0, 2000);
      fQAEvents->Add(fhQAEventsfMultTPCvsTOF);
      fhQAEventsfMultTPCvsESD = new TH2D("fhQAEventsfMultTPCvsESD", "; TPC FB128 multiplicity; ESD multiplicity", 200, 0, 7000, 300, -1000, 35000);
      fQAEvents->Add(fhQAEventsfMultTPCvsESD);
    }

    // charged (tracks) histograms
    fhRefsMult = new TH1D("fhRefsMult","RFPs: Multiplicity; multiplicity", 200,0,1000);
    fQACharged->Add(fhRefsMult);
    fhRefsPt = new TH1D("fhRefsPt","RFPs: #it{p}_{T};  #it{p}_{T} (GeV/#it{c})", iFlowRFPsPtBinNum,fFlowRFPsPtMin,fFlowRFPsPtMax);
    fQACharged->Add(fhRefsPt);
    fhRefsEta = new TH1D("fhRefsEta","RFPs: #eta; #eta", fFlowEtaBinNum,-fFlowEtaMax,fFlowEtaMax);
    fQACharged->Add(fhRefsEta);
    fhRefsPhi = new TH1D("fhRefsPhi","RFPs: #varphi; #varphi", fFlowPhiBinNum,0.0,TMath::TwoPi());
    fQACharged->Add(fhRefsPhi);
    fpRefsMult = new TProfile("fpRefsMult","Ref mult; %s", fCentBinNum,fCentMin,fCentMax);
    fpRefsMult->Sumw2();
    fQACharged->Add(fpRefsMult);

    // PID tracks histograms
    TString sNamePID[3] = {"Pion","Kaon","Proton"};
    TString sLabelPID[3] = {"#pi","K","p"};

    for(Int_t iPID(0); iPID < 3; ++iPID)
    {
      if(!fProcessSpec[iPID+2]) { continue; }

      fhPIDMult[iPID] = new TH1D(Form("fhPID%sMult",sNamePID[iPID].Data()),Form("PID: %s: Multiplicity; multiplicity",sLabelPID[iPID].Data()), 200,0,200);
      fQAPID->Add(fhPIDMult[iPID]);
      fhPIDPt[iPID] = new TH1D(Form("fhPID%sPt",sNamePID[iPID].Data()),Form("PID: %s: #it{p}_{T}; #it{p}_{T}",sLabelPID[iPID].Data()), fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax);
      fQAPID->Add(fhPIDPt[iPID]);
      fhPIDPhi[iPID] = new TH1D(Form("fhPID%sPhi",sNamePID[iPID].Data()),Form("PID: %s: #varphi; #varphi",sLabelPID[iPID].Data()), fFlowPhiBinNum,0,TMath::TwoPi());
      fQAPID->Add(fhPIDPhi[iPID]);
      fhPIDEta[iPID] = new TH1D(Form("fhPID%sEta",sNamePID[iPID].Data()),Form("PID: %s: #eta; #eta",sLabelPID[iPID].Data()), fFlowEtaBinNum,-fFlowEtaMax,fFlowEtaMax);
      fQAPID->Add(fhPIDEta[iPID]);
      fhPIDCharge[iPID] = new TH1D(Form("fhPID%sCharge",sNamePID[iPID].Data()),Form("PID: %s: charge; charge",sLabelPID[iPID].Data()), 3,-1.5,1.5);
      fQAPID->Add(fhPIDCharge[iPID]);
      fh2PIDTPCdEdx[iPID] = new TH2D(Form("fh2PID%sTPCdEdx",sNamePID[iPID].Data()),Form("PID: %s: TPC dE/dx; #it{p} (GeV/#it{c}); TPC dE/dx",sLabelPID[iPID].Data()), fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax, 131,-10,1000);
      fQAPID->Add(fh2PIDTPCdEdx[iPID]);
      fh2PIDTPCdEdxDelta[iPID] = new TH2D(Form("fh2PID%sTPCdEdxDelta",sNamePID[iPID].Data()),Form("PID: %s: TPC #DeltadE/dx; #it{p} (GeV/#it{c}); TPC #DeltadE/dx",sLabelPID[iPID].Data()), fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax, 200,-200,200);
      fQAPID->Add(fh2PIDTPCdEdxDelta[iPID]);
      fh2PIDTOFbeta[iPID] = new TH2D(Form("fh2PID%sTOFbeta",sNamePID[iPID].Data()),Form("PID: %s: TOF #beta; #it{p} (GeV/#it{c});TOF #beta",sLabelPID[iPID].Data()), fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax, 101,-0.1,1.5);
      fQAPID->Add(fh2PIDTOFbeta[iPID]);
      fh2PIDTOFbetaDelta[iPID] = new TH2D(Form("fh2PID%sTOFbetaDelta",sNamePID[iPID].Data()),Form("PID: %s: TOF #Delta#beta; #it{p} (GeV/#it{c});TOF #Delta#beta",sLabelPID[iPID].Data()), fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax, 100,-5000,5000);
      fQAPID->Add(fh2PIDTOFbetaDelta[iPID]);
      fh2PIDTPCnSigmaElectron[iPID] = new TH2D(Form("fh2PID%sTPCnSigmaElectron",sNamePID[iPID].Data()),Form("PID: %s: TPC n#sigma (e hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma",sLabelPID[iPID].Data()), fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax, 21,-11,10);
      fQAPID->Add(fh2PIDTPCnSigmaElectron[iPID]);
      fh2PIDTOFnSigmaElectron[iPID] = new TH2D(Form("fh2PID%sTOFnSigmaElectron",sNamePID[iPID].Data()),Form("PID: %s: TOF n#sigma (e hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma",sLabelPID[iPID].Data()), fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax, 21,-11,10);
      fQAPID->Add(fh2PIDTOFnSigmaElectron[iPID]);
      fh2PIDTPCnSigmaMuon[iPID] = new TH2D(Form("fh2PID%sTPCnSigmaMuon",sNamePID[iPID].Data()),Form("PID: %s: TPC n#sigma (#mu hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma",sLabelPID[iPID].Data()), fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax, 21,-11,10);
      fQAPID->Add(fh2PIDTPCnSigmaMuon[iPID]);
      fh2PIDTOFnSigmaMuon[iPID] = new TH2D(Form("fh2PID%sTOFnSigmaMuon",sNamePID[iPID].Data()),Form("PID: %s: TOF n#sigma (#mu hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma",sLabelPID[iPID].Data()), fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax, 21,-11,10);
      fQAPID->Add(fh2PIDTOFnSigmaMuon[iPID]);
      fh2PIDTPCnSigmaPion[iPID] = new TH2D(Form("fh2PID%sTPCnSigmaPion",sNamePID[iPID].Data()),Form("PID: %s: TPC n#sigma (#pi hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma",sLabelPID[iPID].Data()), fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax, 21,-11,10);
      fQAPID->Add(fh2PIDTPCnSigmaPion[iPID]);
      fh2PIDTOFnSigmaPion[iPID] = new TH2D(Form("fh2PID%sTOFnSigmaPion",sNamePID[iPID].Data()),Form("PID: %s: TOF n#sigma (#pi hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma",sLabelPID[iPID].Data()), fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax, 21,-11,10);
      fQAPID->Add(fh2PIDTOFnSigmaPion[iPID]);
      fh2PIDTPCnSigmaKaon[iPID] = new TH2D(Form("fh2PID%sTPCnSigmaKaon",sNamePID[iPID].Data()),Form("PID: %s: TPC n#sigma (K hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma",sLabelPID[iPID].Data()), fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax, 21,-11,10);
      fQAPID->Add(fh2PIDTPCnSigmaKaon[iPID]);
      fh2PIDTOFnSigmaKaon[iPID] = new TH2D(Form("fh2PID%sTOFnSigmaKaon",sNamePID[iPID].Data()),Form("PID: %s: TOF n#sigma (K hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma",sLabelPID[iPID].Data()), fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax, 21,-11,10);
      fQAPID->Add(fh2PIDTOFnSigmaKaon[iPID]);
      fh2PIDTPCnSigmaProton[iPID] = new TH2D(Form("fh2PID%sTPCnSigmaProton",sNamePID[iPID].Data()),Form("PID: %s: TPC n#sigma (p hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma",sLabelPID[iPID].Data()), fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax, 21,-11,10);
      fQAPID->Add(fh2PIDTPCnSigmaProton[iPID]);
      fh2PIDTOFnSigmaProton[iPID] = new TH2D(Form("fh2PID%sTOFnSigmaProton",sNamePID[iPID].Data()),Form("PID: %s: TOF n#sigma (p hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma",sLabelPID[iPID].Data()), fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax, 21,-11,10);
      fQAPID->Add(fh2PIDTOFnSigmaProton[iPID]);
      fh2PIDBayesElectron[iPID] = new TH2D(Form("fh2PID%sBayesElectron",sNamePID[iPID].Data()),Form("PID: %s: Bayes probability (e hyp.); #it{p}_{T} (GeV/#it{c}); Bayes prob.",sLabelPID[iPID].Data()), fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax, 50,0,1);
      fQAPID->Add(fh2PIDBayesElectron[iPID]);
      fh2PIDBayesMuon[iPID] = new TH2D(Form("fh2PID%sBayesMuon",sNamePID[iPID].Data()),Form("PID: %s: Bayes probability (#mu hyp.); #it{p}_{T} (GeV/#it{c}); Bayes prob.",sLabelPID[iPID].Data()), fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax, 50,0,1);
      fQAPID->Add(fh2PIDBayesMuon[iPID]);
      fh2PIDBayesPion[iPID] = new TH2D(Form("fh2PID%sBayesPion",sNamePID[iPID].Data()),Form("PID: %s: Bayes probability (#pi hyp.); #it{p}_{T} (GeV/#it{c}); Bayes prob.",sLabelPID[iPID].Data()), fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax, 50,0,1);
      fQAPID->Add(fh2PIDBayesPion[iPID]);
      fh2PIDBayesKaon[iPID] = new TH2D(Form("fh2PID%sBayesKaon",sNamePID[iPID].Data()),Form("PID: %s: Bayes probability (K hyp.); #it{p}_{T} (GeV/#it{c}); Bayes prob.",sLabelPID[iPID].Data()), fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax, 50,0,1);
      fQAPID->Add(fh2PIDBayesKaon[iPID]);
      fh2PIDBayesProton[iPID] = new TH2D(Form("fh2PID%sBayesProton",sNamePID[iPID].Data()),Form("PID: %s: Bayes probability (p hyp.); #it{p}_{T} (GeV/#it{c}); Bayes prob.",sLabelPID[iPID].Data()), fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax, 50,0,1);
      fQAPID->Add(fh2PIDBayesProton[iPID]);
    } //end-if {fProcessSpec[kPion]}

    if(fProcessSpec[kPhi])
    {
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
      fhPhiPt = new TH1D("fhPhiPt","#phi: #it{p}_{T}; #it{p}_{T} (GeV/#it{c})", fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax);
      fQAPhi->Add(fhPhiPt);
      fhPhiEta = new TH1D("fhPhiEta","#phi: #eta; #eta", fFlowEtaBinNum,-fFlowEtaMax,fFlowEtaMax);
      fQAPhi->Add(fhPhiEta);
      fhPhiPhi = new TH1D("fhPhiPhi","#phi: #varphi; #varphi", fFlowPhiBinNum,0.,TMath::TwoPi());
      fQAPhi->Add(fhPhiPhi);
    } //end-if {fProcessPhi}

    if(fProcessSpec[kK0s] || fProcessSpec[kLambda])
    {
      fhV0sInvMassK0s = new TH2D("fhV0sInvMassK0s","V^{0}: K^{0}_{S}: InvMass (selected); K^{0}_{S} #it{m}_{inv} (GeV/#it{c}^{2}); #Lambda/#bar{#Lambda} #it{m}_{inv} (GeV/#it{c}^{2})", 110,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax, 50,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax);
      fQAV0s->Add(fhV0sInvMassK0s);
      fhV0sInvMassLambda = new TH2D("fhV0sInvMassLambda","V^{0}: #Lambda/#bar{#Lambda}: InvMass (selected); K^{0}_{S} #it{m}_{inv} (GeV/#it{c}^{2}); #Lambda/#bar{#Lambda} #it{m}_{inv} (GeV/#it{c}^{2})", 110,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax, 50,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax);
      fQAV0s->Add(fhV0sInvMassLambda);
      fhV0sCompetingInvMassK0s = new TH2D("fhV0sCompetingInvMassK0s","V^{0}: K^{0}_{S}: Competing InvMass rejection; K^{0}_{S} #it{m}_{inv} (GeV/#it{c}^{2}); #Lambda/#bar{#Lambda} #it{m}_{inv} (GeV/#it{c}^{2})", 110,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax, 50,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax);
      fQAV0s->Add(fhV0sCompetingInvMassK0s);
      fhV0sCompetingInvMassLambda = new TH2D("fhV0sCompetingInvMassLambda","V^{0}: #Lambda/#bar{#Lambda}: Competing InvMass rejection; K^{0}_{S} #it{m}_{inv} (GeV/#it{c}^{2}); #Lambda/#bar{#Lambda} #it{m}_{inv} (GeV/#it{c}^{2})", 110,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax, 50,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax);
      fQAV0s->Add(fhV0sCompetingInvMassLambda);
    } // end-if {fProcessSpec[kK0s] || fProcessSpec[kLambda]}

    // ####  Selection QA (2-step)
    TString sQAindex[QAindex::kNumQA] = {"Before", "After"};
    for(Int_t iQA(0); iQA < Int_t(QAindex::kNumQA); ++iQA)
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
      fhQAChargedMult[iQA] = new TH1D(Form("fhQAChargedMult_%s",sQAindex[iQA].Data()),"QA Charged: Number of Charged in selected events; #it{N}^{Charged}", 150,0,1500);
      fQACharged->Add(fhQAChargedMult[iQA]);
      fhQAChargedCharge[iQA] = new TH1D(Form("fhQAChargedCharge_%s",sQAindex[iQA].Data()),"QA Charged: Track charge; charge;", 3,-1.5,1.5);
      fQACharged->Add(fhQAChargedCharge[iQA]);
      fhQAChargedPt[iQA] = new TH1D(Form("fhQAChargedPt_%s",sQAindex[iQA].Data()),"QA Charged: Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c})", 200,0.,20.);
      fQACharged->Add(fhQAChargedPt[iQA]);
      fhQAChargedEta[iQA] = new TH1D(Form("fhQAChargedEta_%s",sQAindex[iQA].Data()),"QA Charged: Track #it{#eta}; #it{#eta}", 151,-1.5,1.5);
      fQACharged->Add(fhQAChargedEta[iQA]);
      fhQAChargedPhi[iQA] = new TH1D(Form("fhQAChargedPhi_%s",sQAindex[iQA].Data()),"QA Charged: Track #it{#varphi}; #it{#varphi}", 100,0.,TMath::TwoPi());
      fQACharged->Add(fhQAChargedPhi[iQA]);
      fhQAChargedNumTPCcls[iQA] = new TH1D(Form("fhQAChargedNumTPCcls_%s",sQAindex[iQA].Data()),"QA Charged: Track number of TPC clusters; #it{N}^{TPC clusters}", 160,0,160);
      fQACharged->Add(fhQAChargedNumTPCcls[iQA]);
      fhQAChargedDCAxy[iQA] = new TH1D(Form("fhQAChargedDCAxy_%s",sQAindex[iQA].Data()),"QA Charged: Track DCA-xy; DCA_{#it{xy}} (cm)", 100,0.,10);
      fQACharged->Add(fhQAChargedDCAxy[iQA]);
      fhQAChargedDCAz[iQA] = new TH1D(Form("fhQAChargedDCAz_%s",sQAindex[iQA].Data()),"QA Charged: Track DCA-z; DCA_{#it{z}} (cm)", 200,-10.,10.);
      fQACharged->Add(fhQAChargedDCAz[iQA]);

      if(fProcessSpec[kPion] || fProcessSpec[kKaon] || fProcessSpec[kProton] || fProcessSpec[kPhi])
      {
        fhQAPIDTPCstatus[iQA] = new TH1D(Form("fhQAPIDTPCstatus_%s",sQAindex[iQA].Data()),"QA PID: PID status: TPC;", iNBinsPIDstatus,0,iNBinsPIDstatus);
        fQAPID->Add(fhQAPIDTPCstatus[iQA]);
        fhQAPIDTPCdEdx[iQA] = new TH2D(Form("fhQAPIDTPCdEdx_%s",sQAindex[iQA].Data()),"QA PID: TPC PID information; #it{p} (GeV/#it{c}); TPC dEdx (au)", 100,0,10, 131,-10,1000);
        fQAPID->Add(fhQAPIDTPCdEdx[iQA]);
        fhQAPIDTOFstatus[iQA] = new TH1D(Form("fhQAPIDTOFstatus_%s",sQAindex[iQA].Data()),"QA PID: PID status: TOF;", iNBinsPIDstatus,0,iNBinsPIDstatus);
        fQAPID->Add(fhQAPIDTOFstatus[iQA]);
        fhQAPIDTOFbeta[iQA] = new TH2D(Form("fhQAPIDTOFbeta_%s",sQAindex[iQA].Data()),"QA PID: TOF #beta information; #it{p} (GeV/#it{c}); TOF #beta", 100,0,10, 101,-0.1,1.5);
        fQAPID->Add(fhQAPIDTOFbeta[iQA]);
        fh3QAPIDnSigmaTPCTOFPtPion[iQA] = new TH3D(Form("fh3QAPIDnSigmaTPCTOFPtPion_%s",sQAindex[iQA].Data()), "QA PID: nSigma Pion vs. p_{T}; n#sigma TPC; n#sigma TOF; p_{T} (GeV/c)", 21,-11,10, 21,-11,10, fFlowPOIsPtBinNum, fFlowPOIsPtMin, fFlowPOIsPtMax);
        fQAPID->Add(fh3QAPIDnSigmaTPCTOFPtPion[iQA]);
        fh3QAPIDnSigmaTPCTOFPtKaon[iQA] = new TH3D(Form("fh3QAPIDnSigmaTPCTOFPtKaon_%s",sQAindex[iQA].Data()), "QA PID: nSigma Kaon vs. p_{T}; n#sigma TPC; n#sigma TOF; p_{T} (GeV/c)", 21,-11,10, 21,-11,10, fFlowPOIsPtBinNum, fFlowPOIsPtMin, fFlowPOIsPtMax);
        fQAPID->Add(fh3QAPIDnSigmaTPCTOFPtKaon[iQA]);
        fh3QAPIDnSigmaTPCTOFPtProton[iQA] = new TH3D(Form("fh3QAPIDnSigmaTPCTOFPtProton_%s",sQAindex[iQA].Data()), "QA PID: nSigma Proton vs. p_{T}; n#sigma TPC; n#sigma TOF; p_{T} (GeV/c)", 21,-11,10, 21,-11,10, fFlowPOIsPtBinNum, fFlowPOIsPtMin, fFlowPOIsPtMax);
        fQAPID->Add(fh3QAPIDnSigmaTPCTOFPtProton[iQA]);

        for(Int_t j = 0; j < iNBinsPIDstatus; ++j)
        {
          fhQAPIDTOFstatus[iQA]->GetXaxis()->SetBinLabel(j+1, sPIDstatus[j].Data() );
          fhQAPIDTPCstatus[iQA]->GetXaxis()->SetBinLabel(j+1, sPIDstatus[j].Data() );
        }
      } // end-if {fProcessSpec[kPion] || fProcessSpec[kKaon] || fProcessSpec[kProton]}

      if(fProcessSpec[kK0s] || fProcessSpec[kLambda])
      {
        fhQAV0sMultK0s[iQA] = new TH1D(Form("fhQAV0sMultK0s_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Number of K^{0}_{S} candidates", 100,0,1000);
        fQAV0s->Add(fhQAV0sMultK0s[iQA]);
        fhQAV0sMultLambda[iQA] = new TH1D(Form("fhQAV0sMultLambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Number of #Lambda candidates", 100,0,1000);
        fQAV0s->Add(fhQAV0sMultLambda[iQA]);
        fhQAV0sMultALambda[iQA] = new TH1D(Form("fhQAV0sMultALambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Number of #bar{#Lambda} candidates", 100,0,1000);
        fQAV0s->Add(fhQAV0sMultALambda[iQA]);
        fhQAV0sRecoMethod[iQA] = new TH1D(Form("fhQAV0sRecoMethod_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Reconstruction method", 2,-0.5,1.5);
        fhQAV0sRecoMethod[iQA]->GetXaxis()->SetBinLabel(1, "offline");
        fhQAV0sRecoMethod[iQA]->GetXaxis()->SetBinLabel(2, "online (on-the-fly)");
        fQAV0s->Add(fhQAV0sRecoMethod[iQA]);
        fhQAV0sDCAtoPV[iQA] = new TH1D(Form("fhQAV0sDCAtoPV_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter DCA to PV; daughter DCA^{PV} (cm)", 100,0.0,5.0);
        fQAV0s->Add(fhQAV0sDCAtoPV[iQA]);
        fhQAV0sDCADaughters[iQA] = new TH1D(Form("fhQAV0sDCADaughters_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: DCA among daughters; DCA^{daughters} (cm)", 100,0.,10.);
        fQAV0s->Add(fhQAV0sDCADaughters[iQA]);
        fhQAV0sDecayRadius[iQA] = new TH1D(Form("fhQAV0sDecayRadius_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Decay radius; #it{r_{xy}}^{decay} (cm)", 200,0.0,200.0);
        fQAV0s->Add(fhQAV0sDecayRadius[iQA]);
        fhQAV0sCPAK0s[iQA] = new TH1D(Form("fhQAV0sCPAK0s_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: K^{0}_{S}: CPA; CPA^{K0s}", 100,0.9,1.);
        fQAV0s->Add(fhQAV0sCPAK0s[iQA]);
        fhQAV0sCPALambda[iQA] = new TH1D(Form("fhQAV0sCPALambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: #Lambda/#bar{#Lambda}: CPA; CPA^{#Lambda}", 100, 0.9,1.);
        fQAV0s->Add(fhQAV0sCPALambda[iQA]);
        fhQAV0sNumTauK0s[iQA] = new TH1D(Form("fhQAV0sNumTauK0s_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}:  K^{0}_{S}: Number of #it{c#tau}; #it{c#tau}^{K0s} (cm)", 60, 0.0,60.0);
        fQAV0s->Add(fhQAV0sNumTauK0s[iQA]);
        fhQAV0sNumTauLambda[iQA] = new TH1D(Form("fhQAV0sNumTauLambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Number of #it{c#tau}; #it{c#tau}^{#Lambda} (cm)", 60, 0.0,60.0);
        fQAV0s->Add(fhQAV0sNumTauLambda[iQA]);
        fhQAV0sArmenterosK0s[iQA] = new TH2D(Form("fhQAV0sArmenterosK0s_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}:  K^{0}_{S}: Armenteros-Podolaski plot; #alpha; #it{p}_{T}^{Arm} (GeV/#it{c});", 100,-1.,1., 60,0.,0.3);
        fQAV0s->Add(fhQAV0sArmenterosK0s[iQA]);
        fhQAV0sArmenterosLambda[iQA] = new TH2D(Form("fhQAV0sArmenterosLambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: #Lambda: Armenteros-Podolaski plot; #alpha; #it{p}_{T}^{Arm} (GeV/#it{c});", 100,-1.,1., 60,0.,0.3);
        fQAV0s->Add(fhQAV0sArmenterosLambda[iQA]);
        fhQAV0sArmenterosALambda[iQA] = new TH2D(Form("fhQAV0sArmenterosALambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: #bar{#Lambda}: Armenteros-Podolaski plot; #alpha; #it{p}_{T}^{Arm} (GeV/#it{c});", 100,-1.,1., 60,0.,0.3);
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
        fhQAV0sMotherRapLambda[iQA] = new TH1D(Form("fhQAV0sMotherRapLambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Mother #it{y} (Lambda/#bar{#Lambda} hypo); #it{y}^{V0,#Lambda}", 151,-1.5,1.5);
        fQAV0s->Add(fhQAV0sMotherRapLambda[iQA]);
        fhQAV0sDaughterTPCRefit[iQA] = new TH1D(Form("fhQAV0sDaughterTPCRefit_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter TPC refit", 2,-0.5,1.5);
        fhQAV0sDaughterTPCRefit[iQA]->GetXaxis()->SetBinLabel(1, "NOT AliAODTrack::kTPCrefit");
        fhQAV0sDaughterTPCRefit[iQA]->GetXaxis()->SetBinLabel(2, "AliAODTrack::kTPCrefit");
        fQAV0s->Add(fhQAV0sDaughterTPCRefit[iQA]);
        fhQAV0sDaughterKinks[iQA] = new TH1D(Form("fhQAV0sDaughterKinks_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter Kinks", 2,-0.5,1.5);
        fhQAV0sDaughterKinks[iQA]->GetXaxis()->SetBinLabel(1, "NOT AliAODVertex::kKink");
        fhQAV0sDaughterKinks[iQA]->GetXaxis()->SetBinLabel(2, "AliAODVertex:kKink");
        fQAV0s->Add(fhQAV0sDaughterKinks[iQA]);
        fhQAV0sDaughterNumTPCCls[iQA] = new TH1D(Form("fhQAV0sDaughterNumTPCCls_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter # of TPC clusters; # cls; Counts;", 165,-5,160);
        fQAV0s->Add(fhQAV0sDaughterNumTPCCls[iQA]);
        fhQAV0sDaughterNumTPCClsPID[iQA] = new TH1D(Form("fhQAV0sDaughterNumTPCClsPID_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter # of TPC clusters for PID; # cls PID; Counts;", 165,-5,160);
        fQAV0s->Add(fhQAV0sDaughterNumTPCClsPID[iQA]);
        fhQAV0sDaughterNumTPCFind[iQA] = new TH1D(Form("fhQAV0sDaughterNumTPCFind_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter # of findable TPC clusters; # findable; Counts;", 165,-5,160);
        fQAV0s->Add(fhQAV0sDaughterNumTPCFind[iQA]);
        fhQAV0sDaughterNumTPCCrossRows[iQA] = new TH1D(Form("fhQAV0sDaughterNumTPCCrossRows_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter # of crossed TPC rows; # crossed; Counts;", 165,-5,160);
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

        for(Int_t j = 0; j < iNBinsPIDstatus; ++j)
        {
          fhQAV0sDaughterTOFstatus[iQA]->GetXaxis()->SetBinLabel(j+1, sPIDstatus[j].Data() );
          fhQAV0sDaughterTPCstatus[iQA]->GetXaxis()->SetBinLabel(j+1, sPIDstatus[j].Data() );
        }
      } // end-if {fProcessSpec[kK0s] || fProcessSpec[kLambda]}

    } // end-for {iQA}
  } // end-if {fFillQA}

  if(fMC) {
    // NUE weights
    for(Int_t iSpec(kCharged); iSpec < kUnknown; ++iSpec) {
        // if(!fProcessSpec[iSpec]) { continue; }

        Double_t dPtLow = 0.0;
        Double_t dPtHigh = fFlowPOIsPtMax;
        Int_t iNumBinsPt = floor(fFlowPOIsPtMax / 0.1 + 0.5);

        fh2MCPtEtaGen[iSpec] = new TH2D(Form("fh2MCPtEtaGen%s",GetSpeciesName(iSpec)),Form("MC %s (Gen); #it{p}_{T} (GeV/#it{c}); #it{#eta}", GetSpeciesLabel(iSpec)), iNumBinsPt,dPtLow,dPtHigh, fFlowEtaBinNum,-fFlowEtaMax,fFlowEtaMax);
        fListMC->Add(fh2MCPtEtaGen[iSpec]);
        fh2MCPtEtaReco[iSpec] = new TH2D(Form("fh2MCPtEtaReco%s",GetSpeciesName(iSpec)),Form("MC %s (reco); #it{p}_{T} (GeV/#it{c}); #it{#eta}", GetSpeciesLabel(iSpec)), iNumBinsPt,dPtLow,dPtHigh, fFlowEtaBinNum,-fFlowEtaMax,fFlowEtaMax);
        fListMC->Add(fh2MCPtEtaReco[iSpec]);
        fh2MCPtEtaRecoTrue[iSpec] = new TH2D(Form("fh2MCPtEtaRecoTrue%s",GetSpeciesName(iSpec)),Form("MC %s (reco + true); #it{p}_{T} (GeV/#it{c}); #it{#eta}", GetSpeciesLabel(iSpec)), iNumBinsPt,dPtLow,dPtHigh, fFlowEtaBinNum,-fFlowEtaMax,fFlowEtaMax);
        fListMC->Add(fh2MCPtEtaRecoTrue[iSpec]);
    }

  } // end-if{fMC}

  // posting data (mandatory)
  Int_t i = 0;
  for(Int_t iSpec(0); iSpec < kUnknown; ++iSpec) { PostData(++i, fListFlow[iSpec]); }
  PostData(++i, fQAEvents);
  PostData(++i, fQACharged);
  PostData(++i, fQAPID);
  PostData(++i, fQAV0s);
  PostData(++i, fQAPhi);
  PostData(++i, fFlowWeights);
  if(fMC) { PostData(++i, fListMC); }

  DumpTObjTable("UserCreateOutputObjects: end");

  return;
}

#endif
