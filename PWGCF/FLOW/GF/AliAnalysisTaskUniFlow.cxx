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
#include <array>
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
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliGenEventHeader.h"
#include "AliCollisionGeometry.h"
#include "AliGenHijingEventHeader.h"
#include "AliPDG.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliPicoTrack.h"
#include "AliAODv0.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliEventPoolManager.h"

#include "AliAnalysisTaskUniFlow.h"
#include "AliUniFlowCorrTask.h"

ClassImp(AliAnalysisTaskUniFlow);

namespace {
    const Int_t fPDGCode[] = {0,0,211,321,2212,0,310,3122,333};
    const Double_t fPDGMass[] = {0,0,0.13957,0.493677,0.938272,0,0.497614,1.11568,1.019455};
}

AliAnalysisTaskUniFlow::AliAnalysisTaskUniFlow() : AliAnalysisTaskSE(),
  fEventCuts{},
  fEventAOD{nullptr},
  // fEventMC{nullptr},
  fEvent{nullptr},
  fPVz{},
  fPIDResponse{nullptr},
  fPIDCombined{nullptr},
  fFlowWeightsList{nullptr},
  fMC{kFALSE},
  fNeedPIDCorrection{kFALSE},
  fPIDCorrectionPhi{kFALSE},
  fIs2018data{kFALSE},
  fIsHMpp{kFALSE},
  fInit{kFALSE},
  fUseGeneralFormula{kFALSE},
  fFlowUsePIDWeights{kFALSE},
  fPIDonlyForRefs{kFALSE},
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
  fFlowVecPmid{},
  fFlowVecSpos{},
  fFlowVecSneg{},
  fFlowVecSmid{},
  fFlowVecQ{},
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
  fFlowFillWeightsMultiD{kFALSE},
  fFlowFillAfterWeights{kTRUE},
  fFlowUseWeights{kFALSE},
  fFlowUse3Dweights{kFALSE},
  fFlowRunByRunWeights{kTRUE},
  fFlowPeriodWeights{kFALSE},
  fFlowWeightsApplyForReco{kTRUE},
  fFlowWeightsTag{},
  fEventPoolMgr{nullptr},
  fPool{nullptr},
  fSelectedTracks{nullptr},
  fCorrUsingGF{kFALSE},
  fCorrFill{kFALSE},
  fFillMixed{kTRUE},
  fUsePtBinnedEventPool{kTRUE},
  fPoolSize{-1},
  fMixingTracks{5000},
  fMinEventsToMix{5},
  fCorrDEtaBinNum{32},
  fCorrDPhiBinNum{72},
  fCorrdEtaMin{-1.6},
  fCorrdEtaMax{1.6},
  fCorrdPhiMin{-0.5*TMath::Pi()},
  fCorrdPhiMax{1.5*TMath::Pi()},
  fEtaSlicesArr{},
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
  fVxMax{3.},
  fVyMax{3.},
  fVzMax{10.},
  fImpactParameterMC{0.},
  fEventRejectAddPileUp{kFALSE},
  fPileUpCutESDTPC{500},
  fPileUpCutCentrality{10},
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
  fhWeightsMultiD{nullptr},
  fhEventSampling{nullptr},
  fhEventCentrality{nullptr},
  fh2EventCentralityNumRefs{nullptr},
  fhEventCounter{nullptr},
  fhV0Mamplitude{nullptr},
  fhV0MamplitudeRatio{nullptr},
  fh2V0MnCharged{nullptr},
  fh2MeanMultRFP{nullptr},
  fh2MCip{nullptr},
  fhRefsMult{nullptr},
  fhRefsPt{nullptr},
  fhRefsEta{nullptr},
  fhRefsPhi{nullptr},
  fpRefsMult{nullptr},
  fhChargedCounter{nullptr},
  fh4CorrelationsSE{nullptr},
  fh4CorrelationsME{nullptr},
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
  // fEventMC{nullptr},
  fEvent{nullptr},
  fPVz{},
  fPIDResponse{nullptr},
  fPIDCombined{nullptr},
  fFlowWeightsList{nullptr},
  fMC{bIsMC},
  fNeedPIDCorrection{kFALSE},
  fPIDCorrectionPhi{kFALSE},
  fIs2018data{kFALSE},
  fIsHMpp{kFALSE},
  fInit{kFALSE},
  fUseGeneralFormula{kFALSE},
  fFlowUsePIDWeights{kFALSE},
  fPIDonlyForRefs{kFALSE},
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
  fFlowVecPmid{},
  fFlowVecSpos{},
  fFlowVecSneg{},
  fFlowVecSmid{},
  fFlowVecQ{},
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
  fFlowFillWeightsMultiD{kFALSE},
  fFlowFillAfterWeights{kTRUE},
  fFlowUseWeights{bUseWeights},
  fFlowUse3Dweights{kFALSE},
  fFlowRunByRunWeights{kTRUE},
  fFlowPeriodWeights{kFALSE},
  fFlowWeightsApplyForReco{kTRUE},
  fFlowWeightsTag{},
  fEventPoolMgr{nullptr},
  fPool{nullptr},
  fSelectedTracks{nullptr},
  fCorrUsingGF{kFALSE},
  fCorrFill{kFALSE},
  fFillMixed{kTRUE},
  fUsePtBinnedEventPool{kTRUE},
  fPoolSize{-1},
  fMixingTracks{5000},
  fMinEventsToMix{5},
  fCorrDEtaBinNum{32},
  fCorrDPhiBinNum{72},
  fCorrdEtaMin{-1.6},
  fCorrdEtaMax{1.6},
  fCorrdPhiMin{-0.5*TMath::Pi()},
  fCorrdPhiMax{1.5*TMath::Pi()},
  fEtaSlicesArr{},
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
  fVxMax{3.},
  fVyMax{3.},
  fVzMax{10.},
  fImpactParameterMC{0.},
  fEventRejectAddPileUp{kFALSE},
  fPileUpCutESDTPC{500},
  fPileUpCutCentrality{10},
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
  fhWeightsMultiD{nullptr},
  fhEventSampling{nullptr},
  fhEventCentrality{nullptr},
  fh2EventCentralityNumRefs{nullptr},
  fhEventCounter{nullptr},
  fhV0Mamplitude{nullptr},
  fhV0MamplitudeRatio{nullptr},
  fh2V0MnCharged{nullptr},
  fh2MeanMultRFP{nullptr},
  fh2MCip{nullptr},
  fhRefsMult{nullptr},
  fhRefsPt{nullptr},
  fhRefsEta{nullptr},
  fhRefsPhi{nullptr},
  fpRefsMult{nullptr},
  fhChargedCounter{nullptr},
  fh4CorrelationsSE{nullptr},
  fh4CorrelationsME{nullptr},
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
  DefineOutput(15, TList::Class());
  if(fMC) DefineOutput(16, TList::Class());
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

  //deleting event mixing variables
  if(fEventPoolMgr) delete fEventPoolMgr;
  if(fPool) delete fPool;
  if(fSelectedTracks) delete fSelectedTracks;

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
    case kCharUnidentified: name = "UnidentifiedCharged"; break;
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
    case kCharUnidentified: label = "h^{#pm}_{unID}"; break;
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
  printf("      fUseGeneralFormula: (Bool_t) %s\n",    fUseGeneralFormula ? "kTRUE" : "kFALSE");
  printf("      fIsHMpp: (Bool_t) %s\n",    fIsHMpp ? "kTRUE" : "kFALSE");
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
  printf("      fFlowFillWeightsMultiD: (Bool_t) %s\n",    fFlowFillWeightsMultiD ? "kTRUE" : "kFALSE");
  printf("      fFlowFillAfterWeights: (Bool_t) %s\n",    fFlowFillAfterWeights ? "kTRUE" : "kFALSE");
  printf("      fFlowUseWeights: (Bool_t) %s\n",    fFlowUseWeights ? "kTRUE" : "kFALSE");
  printf("      fFlowWeightsTag: (TString) '%s'\n",    fFlowWeightsTag.Data());
  printf("      fFlowRunByRunWeights: (Bool_t) %s\n",    fFlowRunByRunWeights ? "kTRUE" : "kFALSE");
  printf("      fFlowPeriodWeights: (Bool_t) %s\n",    fFlowPeriodWeights ? "kTRUE" : "kFALSE");
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
    std::vector<AliVParticle*>* vector = fVector[iSpec];
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
  if(fAnalType != kAOD && fAnalType != kMC)
  {
    AliFatal("Analysis type: not kAOD (not implemented for ESDs yet)! Terminating!");
    return kFALSE;
  }

  // checking PID response
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*)mgr->GetInputEventHandler();
  fPIDResponse = inputHandler->GetPIDResponse();
  if(!fPIDResponse && fAnalType != kMC)
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

  if(fProcessSpec[kPion] || fProcessSpec[kKaon] || fProcessSpec[kProton]) fProcessSpec[kCharUnidentified] = kTRUE;

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
    if(!fFlowRunByRunWeights  && !fFlowPeriodWeights && !LoadWeights()) { AliFatal("Initial flow weights not loaded! Terminating!"); return kFALSE; }
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
    fVector[iSpec] = new std::vector<AliVParticle*>();
    fVector[iSpec]->reserve(iReserve);
  }

  if(fCorrUsingGF){
    for(Int_t iBin(0); iBin < fFlowBinNumberEtaSlices+1; iBin++){
      fEtaSlicesArr[iBin] = (2.0*fFlowEtaMax/fFlowBinNumberEtaSlices)*iBin - fFlowEtaMax;
    }
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
  // fEventAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  // if(!fEventAOD) { return; }

  // loading AliPIDResponse
  AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
  AliVEventHandler* inputHandler = (AliVEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();
  if(!fPIDResponse && fAnalType != kMC) { AliFatal("AliPIDResponse not attached!"); return; }

  // loading array with MC particles
  // if(fMC) {
  //   fEventMC = inputHandler->MCEvent();
  //   if(!fEventMC) { AliFatal("fEventMC with MC particle not found!"); return; }
  // }

  //loading VEvent with generated event / data
  if(fAnalType == kMC){ fEvent = dynamic_cast<AliVEvent*>(MCEvent());}
  else{
    fEvent = dynamic_cast<AliVEvent*>(InputEvent());
    fEventAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fEventAOD) { AliFatal("AOD not found!"); return; }
  }
  if(!fEvent) { AliFatal("fEvent not found!"); return; }

  // "valid" events before selection
  fhEventCounter->Fill("Input",1);

  // Fill event QA BEFORE cuts
  if(fFillQA) { FillQAEvents(kBefore); }

  // extract PV-z for weights
  fPVz = fEvent->GetPrimaryVertex()->GetZ();

  Bool_t bEventSelected = kFALSE;
  if(fAnalType != kMC) bEventSelected = IsEventSelected();
  else bEventSelected = IsMCEventSelected();

  DumpTObjTable("UserExec: after event selection");
  if(!bEventSelected) { return; }

  fhEventCounter->Fill("Event OK",1);

  if(fAnalType == kMC && fFlowUseWeights) { AliFatal("Cannot generate events and use weights on in the same time! Terminating!"); return; }

  // checking the run number for aplying weights & loading TList with weights
  if(fAnalType != kMC && fFlowUseWeights && (fFlowRunByRunWeights || fFlowPeriodWeights) && fRunNumber != fEventAOD->GetRunNumber() && !LoadWeights()) { AliFatal("Weights not loaded!"); return; }

  DumpTObjTable("UserExec: before filtering");

  // Filter charged (& Refs) particles to evaluate event multiplcity / N_RFPs
  // NB: clear charged vectors because it might keep particles from previous event (not happen for other species)
  fVector[kRefs]->clear();
  fVector[kCharged]->clear();

  if(!fMC) FilterCharged();
  else FilterChargedMC();

  if(fIsHMpp && fFillQA) {
    AliMultSelection* multSelection = (AliMultSelection*) fEventAOD->FindListObject("MultSelection");
    if(!multSelection) { AliError("AliMultSelection object not found! Returning -1"); return; }
    AliMultEstimator* lEst = multSelection->GetEstimator("V0M");
    fh2V0MnCharged->Fill(lEst->GetValue()/lEst->GetMean(),fVector[kCharged]->size());
  }

  // checking if there is at least 4/6/8 particles: needed to "properly" calculate correlations
  UInt_t minNOfPar = 4;
  if(fColSystem == kPbPb) minNOfPar = 8;
  if(fUseGeneralFormula) minNOfPar = 12;
  if(fVector[kRefs]->size() <= minNOfPar) { return; }

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
    for(Int_t iTrack(0); iTrack < fEvent->GetNumberOfTracks(); iTrack++) {
      AliVParticle* track = static_cast<AliVParticle*>(fEvent->GetTrack(iTrack));
      if(!track) { continue; }
      FillQACharged(kBefore,track);
    }
    fhQAChargedMult[0]->Fill(fEvent->GetNumberOfTracks());
    // fhQAChargedMult[1]->Fill(fVector[kCharged]->size());
    fhRefsMult->Fill(fVector[kRefs]->size());
  }

  // sorting charged hadrons
  std::sort(fVector[kCharged]->begin(), fVector[kCharged]->end(), [this](const AliVParticle* a, const AliVParticle* b){ return this->sortPt(a, b); });

  // Filtering other species
  if(fProcessSpec[kPion] || fProcessSpec[kKaon] || fProcessSpec[kProton]) {
    FilterPID();
  }



  if(fProcessSpec[kK0s] || fProcessSpec[kLambda]) {
      FilterV0s();
      std::sort(fVector[kK0s]->begin(), fVector[kK0s]->end(), [this](const AliVParticle* a, const AliVParticle* b){ return this->sortPt(a, b); });
      std::sort(fVector[kLambda]->begin(), fVector[kLambda]->end(), [this](const AliVParticle* a, const AliVParticle* b){ return this->sortPt(a, b); });
  }
  if(fProcessSpec[kPhi]) {
      FilterPhi();
      std::sort(fVector[kPhi]->begin(), fVector[kPhi]->end(), [this](const AliVParticle* a, const AliVParticle* b){ return this->sortPt(a, b); });
  }

  DumpTObjTable("UserExec: after filtering");

  Bool_t bCorrelations = FillCorrelations();
  DumpTObjTable("UserExec: after FillCorrelations");

  // processing of selected event
  Bool_t bProcessed = CalculateFlow();

  DumpTObjTable("UserExec: after CalculateFlow");

  // should be cleared at the end of processing especially for reconstructed
  // particles (Phi, V0s) because here new AliPicoTracks are created
  ClearVectors();

  // if(fMC) { ProcessMC(); }

  // extracting run number here to store run number from previous event (for current run number use info in AliAODEvent)
  if(fAnalType == kAOD) fRunNumber = fEventAOD->GetRunNumber();

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
  if(fAnalType == kMC) { AliFatal("Simulated event cannot pass 'IsEventSelected'! \n"); return kFALSE; }

  // Physics selection (trigger)
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) mgr->GetInputEventHandler();
  UInt_t fSelectMask = inputHandler->IsEventSelected();
  fhEventCounter->Fill("Loaded OK",1);

  if(!fIs2018data){
    if(!(fSelectMask & fTrigger)) { return kFALSE; }
    // events passing physics && trigger selection
    fhEventCounter->Fill("Triggers OK",1);
  }
  else{
    fEventCuts.SetupPbPb2018();
    fhEventCounter->Fill("2018 OK",1);
  }

  AliMultSelection* multSelection = nullptr;
  if(fIsHMpp){
    if(fColSystem != kPP) {AliWarning("\n\n\n Watch out! Using manual HM pp data for different collision system! \n\n\n"); }

    if(fEventAOD->IsPileupFromSPDInMultBins() ) { return kFALSE; }
    fhEventCounter->Fill("Is not pile up",1);

    multSelection = (AliMultSelection*) fEventAOD->FindListObject("MultSelection");
    if(!multSelection) { AliError("AliMultSelection object not found! Returning -1"); return kFALSE; }
    fhEventCounter->Fill("Multiplicity OK",1);

    if(!multSelection->GetThisEventIsNotPileup() || !multSelection->GetThisEventIsNotPileupInMultBins() || !multSelection->GetThisEventHasNoInconsistentVertices() || !multSelection->GetThisEventPassesTrackletVsCluster()) { return kFALSE; }
    fhEventCounter->Fill("Multiplicity cuts OK",1);

    Int_t nTracksPrim = fEventAOD->GetPrimaryVertex()->GetNContributors();
    if(nTracksPrim < 0.5) { return kFALSE; }
    fhEventCounter->Fill("Contributors OK",1);

    if(TMath::Abs(fPVz) >= fPVtxCutZ) { return kFALSE; }
    fhEventCounter->Fill("PVz OK",1);
  }
  else{
    // events passing AliEventCuts selection
    if(!fEventCuts.AcceptEvent(fEventAOD))  { return kFALSE; }
  }
  fhEventCounter->Fill("EventCuts OK",1);

  // estimate centrality & assign indexes (only if AliMultEstimator is requested)
  if(fCentEstimator != kRFP) {
      fIndexCentrality = GetCentralityIndex(fCentEstimator);

      if(fIndexCentrality < 0) { return kFALSE; }
      if(fCentMin > 0 && fIndexCentrality < fCentMin) { return kFALSE; }
      if(fCentMax > 0 && fIndexCentrality >= fCentMax) { return kFALSE; }
  }
  fhEventCounter->Fill("Centrality OK",1);

  // additional centrality cut for "double differential" cut
  if(fCentEstimatorAdd != kRFP && fCentMaxAdd > 0) {
    Int_t addCentIndex = GetCentralityIndex(fCentEstimatorAdd);
    if(addCentIndex < fCentMinAdd) { return kFALSE; }
    if(addCentIndex > fCentMaxAdd) { return kFALSE; }
  }
  fhEventCounter->Fill("Centrality cuts OK",1);

  if(fIs2018data){
    if((fIndexCentrality<10) || (fIndexCentrality>30 && fIndexCentrality<50)){
      if(!(fSelectMask & (AliVEvent::kCentral|AliVEvent::kSemiCentral|fTrigger))) { return kFALSE; }
    }
    else{
      if(!(fSelectMask & fTrigger)) { return kFALSE;}
    }
    // events passing physics && trigger selection
    fhEventCounter->Fill("Triggers OK",1);
  }

  if(fIsHMpp && fFillQA){
    AliMultEstimator* lEst = multSelection->GetEstimator("V0M");
    fhV0Mamplitude->Fill(lEst->GetValue());
    fhV0MamplitudeRatio->Fill(lEst->GetValue()/lEst->GetMean());
  }

  // Additional pile-up rejection cuts for LHC15o dataset
  if(fColSystem == kPbPb && fEventRejectAddPileUp && fIndexCentrality < fPileUpCutCentrality && IsEventRejectedAddPileUp()) { return kFALSE; }

  fhEventCounter->Fill("PileUp cut OK",1);

  return kTRUE;
}
// ============================================================================
Bool_t AliAnalysisTaskUniFlow::IsMCEventSelected()
{
  // Return kTRUE if event passes all criteria, kFALSE otherwise
  // for simulated events only
  // *************************************************************

  if(fAnalType != kMC) { AliFatal("Real event cannot pass 'IsMCEventSelected'! \n"); return kFALSE; }

  AliMCEvent* ev = dynamic_cast<AliMCEvent*>(fEvent);
  if(!ev) { AliFatal("MC event not found!"); return kFALSE; }

  AliGenEventHeader *header = dynamic_cast<AliGenEventHeader*>(ev->GenEventHeader());
  if(!header) { AliFatal("MC event not generated!"); return kFALSE; }

  fhEventCounter->Fill("Physics selection OK",1);

  const AliVVertex *vertex = ev->GetPrimaryVertex();
  if(!ev) { AliError("Vertex of MC not found!"); }

  if(TMath::Abs(vertex->GetX()) > fVxMax) return kFALSE;
  if(TMath::Abs(vertex->GetY()) > fVyMax) return kFALSE;
  if(TMath::Abs(vertex->GetZ()) > fVzMax) return kFALSE;

  fhEventCounter->Fill("EventCuts OK",1);

  AliCollisionGeometry* headerH;
  TString genName;
  TList *ltgen = (TList*)ev->GetCocktailList();
  if (ltgen) {
  for(auto&& listObject: *ltgen){
    genName = Form("%s",listObject->GetName());
    if (genName.Contains("Hijing")) {
      headerH = dynamic_cast<AliCollisionGeometry*>(listObject);
      break;
      }
    }
  }
  else
    headerH = dynamic_cast<AliCollisionGeometry*>(ev->GenEventHeader());
  if(headerH){
      fImpactParameterMC = headerH->ImpactParameter();
  }

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
    if(multESDTPCdif > fPileUpCutESDTPC) { return kTRUE; }

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
      if(!fFlowRunByRunWeights && !fFlowPeriodWeights) {
          // loading run-averaged weights
          listFlowWeights = (TList*) fFlowWeightsList->FindObject("averaged");
          if(!listFlowWeights) { AliError("TList with flow run-averaged weights not found."); fFlowWeightsList->ls(); return kFALSE; }
      } else if(fFlowPeriodWeights){
        // loading period-specific weights
        listFlowWeights = (TList*) fFlowWeightsList->FindObject(ReturnPPperiod(fEventAOD->GetRunNumber()));
        if(!listFlowWeights) { AliError("Loading period weights failed!"); fFlowWeightsList->ls(); return kFALSE; }
      }
      else {
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
Bool_t AliAnalysisTaskUniFlow::FillFlowWeight(const AliVParticle* track, const PartSpecies species) const
{
  if(!track) { AliError("Track not exists!"); return kFALSE; }
  if(species == kUnknown) { AliError("Invalid species 'Unknown'!"); return kFALSE; }

  if(fFlowFillWeights) {
    if(fFlowUse3Dweights) {
      fh3Weights[species]->Fill(track->Phi(),track->Eta(),fPVz);
    } else {
      fh2Weights[species]->Fill(track->Phi(),track->Eta());
    }

    if(fFlowFillWeightsMultiD){
      if(!fhWeightsMultiD) { Error("THnSparse not valid!","FillFlowWeight"); return kFALSE; }
      Double_t dWeightsValues[SparseWeights::wDim] = {0};
      dWeightsValues[SparseWeights::wPhi] = track->Phi();
      dWeightsValues[SparseWeights::wCent] = fIndexCentrality;
      dWeightsValues[SparseWeights::wPt] = track->Pt();
      dWeightsValues[SparseWeights::wEta] = track->Eta();
      dWeightsValues[SparseWeights::wVz] = fPVz;
      dWeightsValues[SparseWeights::wSpec] = species;
      fhWeightsMultiD->Fill(dWeightsValues);
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
Double_t AliAnalysisTaskUniFlow::GetFlowWeight(const AliVParticle* track, const PartSpecies species) const
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
    fpRefsMult->Fill(fIndexCentrality,fVector[kRefs]->size(),1.0);
    fhEventSampling->Fill(fIndexCentrality,fIndexSampling);
    if(fAnalType == kMC) fh2MCip->Fill(fImpactParameterMC, fIndexCentrality);
  }

  const AliVVertex* aodVtx = fEvent->GetPrimaryVertex();
  const Double_t dVtxZ = aodVtx->GetZ();
  const Int_t iNumContr = aodVtx->GetNContributors();

  fhQAEventsPVz[iQAindex]->Fill(dVtxZ);
  fhQAEventsNumContrPV[iQAindex]->Fill(iNumContr);

  if(fAnalType == kMC) { return; }

  if(iQAindex == 1)
  {
    fhEventCentrality->Fill(fIndexCentrality);
    fh2EventCentralityNumRefs->Fill(fIndexCentrality,fVector[kRefs]->size());
  }

  const AliAODVertex* spdVtx = fEventAOD->GetPrimaryVertexSPD();
  const Int_t iNumContrSPD = spdVtx->GetNContributors();
  const Double_t spdVtxZ = spdVtx->GetZ();
  fhQAEventsNumSPDContrPV[iQAindex]->Fill(iNumContrSPD);
  fhQAEventsDistPVSPD[iQAindex]->Fill(TMath::Abs(dVtxZ - spdVtxZ));

  // SPD vertexer resolution
  Double_t cov[6] = {0};
  spdVtx->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  fhQAEventsSPDresol[iQAindex]->Fill(zRes);

  // // event / physics selection criteria
  // AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  // AliInputEventHandler* inputHandler = (AliInputEventHandler*) mgr->GetInputEventHandler();
  // UInt_t fSelectMask = inputHandler->IsEventSelected();
  //
  // if( fSelectMask& AliVEvent::kINT7 ) { fQAEventsTriggerSelection[iQAindex]->Fill("kINT7",1); }
  // else if (fSelectMask& AliVEvent::kHighMultV0) { fQAEventsTriggerSelection[iQAindex]->Fill("kHighMultV0",1); }
  // else if (fSelectMask& AliVEvent::kHighMultSPD) { fQAEventsTriggerSelection[iQAindex]->Fill("kHighMultSPD",1); }
  // else { fQAEventsTriggerSelection[iQAindex]->Fill("Other",1); }

  return;
}
// ============================================================================
void AliAnalysisTaskUniFlow::ProcessMC() const
{
    AliMCEvent* ev = dynamic_cast<AliMCEvent*>(fEvent);
    if(!ev) { AliFatal("MC event not found!"); return; }

    const Int_t iNumTracksMC = ev->GetNumberOfTracks();
    for(Int_t iTrackMC(0); iTrackMC < iNumTracksMC; ++iTrackMC) {

        AliAODMCParticle* trackMC = (AliAODMCParticle*) ev->GetTrack(iTrackMC);
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

  Int_t iNumTracks = fEvent->GetNumberOfTracks();
  if(iNumTracks < 1) { return; }

  for(Int_t iTrack(0); iTrack < iNumTracks; iTrack++) {
    AliAODTrack* track = static_cast<AliAODTrack*>(fEvent->GetTrack(iTrack));
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
void AliAnalysisTaskUniFlow::FilterChargedMC() const
{
  // Filtering input charged tracks for Monte Carlo POIs (stored in fVector[kCharged]) or RFPs (fVector[kRefs])
  // If track passes all requirements its pointer is pushed to relevant vector container
  // *************************************************************
  AliMCEvent* ev = dynamic_cast<AliMCEvent*>(fEvent);
  if(!ev) { AliFatal("MC event not found!"); return; }

  Int_t iNumTracks = ev->GetNumberOfPrimaries();
  if(iNumTracks < 1) { return; }

  for(Int_t iTrack(0); iTrack < iNumTracks; iTrack++) {
    AliMCParticle* track = dynamic_cast<AliMCParticle*>(ev->GetTrack(iTrack));
    if(!track) { continue; }

    //excluding non stable particles
    if(!(ev->IsPhysicalPrimary(iTrack))) continue;
    if(track->Charge() == 0) continue;

    if(IsWithinRefs(track)) { fVector[kRefs]->push_back(track); }
    if(IsWithinPOIs(track)) { fVector[kCharged]->push_back(track); }
  }

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
  dValues[SparseCand::kSample] = fIndexSampling;
  sparse->Fill(dValues);

  return;
}
// ============================================================================
void AliAnalysisTaskUniFlow::FillQARefs(const QAindex iQAindex, const AliVParticle* track) const
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
void AliAnalysisTaskUniFlow::FillQACharged(const QAindex iQAindex, AliVParticle* track) const
{
  // Filling various QA plots related to charged track selection
  // *************************************************************
  if(!track) return;

  // printf("charge %d, QA idx %d , pdg code %d \n", track->Charge(), iQAindex, track->PdgCode());
  // track charge
  fhQAChargedCharge[iQAindex]->Fill(track->Charge());

  if(fAnalType != kMC){
    AliAODTrack* tr =  dynamic_cast<AliAODTrack*>(track);
    if(!tr) { AliFatal("AOD track not found!"); return; }

    // number of TPC clusters
    fhQAChargedNumTPCcls[iQAindex]->Fill(tr->GetTPCNcls());

    // track DCA
    Double_t dDCAXYZ[3] = {-999., -999., -999.};
    const AliVVertex* vertex = fEventAOD->GetPrimaryVertex();
    if(vertex)
    {
      Double_t dTrackXYZ[3] = {-999., -999., -999.};
      Double_t dVertexXYZ[3] = {-999., -999., -999.};

      tr->GetXYZ(dTrackXYZ);
      vertex->GetXYZ(dVertexXYZ);

      for(Short_t i(0); i < 3; i++)
        dDCAXYZ[i] = dTrackXYZ[i] - dVertexXYZ[i];
    }
    fhQAChargedDCAxy[iQAindex]->Fill(TMath::Sqrt(dDCAXYZ[0]*dDCAXYZ[0] + dDCAXYZ[1]*dDCAXYZ[1]));
    fhQAChargedDCAz[iQAindex]->Fill(dDCAXYZ[2]);
  }

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

  // std::vector<AliAODTrack*> vectorTest[4];

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

    AliAODTrack* daughterPos = (AliAODTrack*) v0->GetDaughter(0);
    AliAODTrack* daughterNeg = (AliAODTrack*) v0->GetDaughter(1);
    if(!daughterPos || !daughterNeg) { AliFatal("Daughters track not found!"); return; }

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

      // vectorTest[0].push_back(daughterPos);
      // vectorTest[0].push_back(daughterNeg);
      //
      // if(IsChargedSelected(daughterPos) && IsWithinPOIs(daughterPos)) vectorTest[1].push_back(daughterPos);
      // if(IsChargedSelected(daughterNeg) && IsWithinPOIs(daughterNeg)) vectorTest[1].push_back(daughterNeg);
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

      // vectorTest[2].push_back(daughterPos);
      // vectorTest[2].push_back(daughterNeg);
      //
      // if(IsChargedSelected(daughterPos) && IsWithinPOIs(daughterPos)) vectorTest[3].push_back(daughterPos);
      // if(IsChargedSelected(daughterNeg) && IsWithinPOIs(daughterNeg)) vectorTest[3].push_back(daughterNeg);

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

    //   vectorTest[2].push_back(daughterPos);
    //   vectorTest[2].push_back(daughterNeg);
    //
    //   if(IsChargedSelected(daughterPos) && IsWithinPOIs(daughterPos)) vectorTest[3].push_back(daughterPos);
    //   if(IsChargedSelected(daughterNeg) && IsWithinPOIs(daughterNeg)) vectorTest[3].push_back(daughterNeg);
    }

    if(bIsK0s && iIsLambda != 0) { fhV0sCounter->Fill("K^{0}_{S} && #Lambda/#bar{#Lambda}",1); }

  } // end-for {v0}

  // printf("\n\n\n\n TESTING V0 daughters \n\n");
  // if(vectorTest[0].size() > 0 && (Double_t) vectorTest[1].size()/vectorTest[0].size() > 0.0) printf("K0s: passed criteria for charged / total: %f ... total: %d, passed criteria: %d \n", (Double_t) vectorTest[1].size()/vectorTest[0].size(), vectorTest[0].size(), vectorTest[1].size());
  // if(vectorTest[2].size() > 2 && (Double_t) vectorTest[3].size()/vectorTest[2].size() > 0.0) printf("Lambda: passed criteria for charged / total: %f ... total: %d, passed criteria: %d \n and that is with # of charged: %d ", (Double_t) vectorTest[3].size()/vectorTest[2].size(), vectorTest[2].size(), vectorTest[3].size(), fVector[kCharged]->size());
  //
  // for(Int_t i(0); i < 4; i++) vectorTest[i].clear();

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
  AliMCEvent* ev = dynamic_cast<AliMCEvent*>(fEvent);
  if(!ev) { AliFatal("MC event not found!"); return nullptr; }

  const Int_t labelAbs = TMath::Abs(label);
  // Negative label just indicate track with shared clustes, but otherwise should be used
  // absolute value has to be used
  // if(label < 0) { /*AliWarning("MC label negative");*/ return nullptr; }

  AliAODMCParticle* mcTrack = (AliAODMCParticle*) ev->GetTrack(labelAbs);
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
    AliVParticle* track = static_cast<AliVParticle*>(*part);
    if(!track) { continue; }

    if((fCentEstimator != kRFP) && (fColSystem == kPP || fColSystem == kPPb)) {
      Int_t counter = fIndexCentrality/10;
      fh2MeanMultRFP[counter]->Fill(track->Pt(), fVector[kCharged]->size());
    }

    fhPIDCounter->Fill("Input",1);

    if(fFillQA) { FillQAPID(kBefore,track,kUnknown); } // filling QA for tracks before selection (but after charged criteria applied)

    // PID track selection (return most favourable species)
    PartSpecies species = kUnknown;
    if(fAnalType != kMC) species = IsPIDSelected(track);
    else species = IsPIDSelectedMC(track);
    if(species != kPion && species != kKaon && species != kProton) {
      // fVector[kCharUnidentified]->push_back(track);
      species = kCharUnidentified; }

    //check pT ranges
    if(fFlowPOIsPtBinEdges[species].size() > 0){
      if(track->Pt() < fFlowPOIsPtBinEdges[species].front() || track->Pt() > fFlowPOIsPtBinEdges[species].back() ) {
        // fVector[kCharUnidentified]->push_back(track);
        species = kCharUnidentified; }
    }

    // check if only protons should be used
    if(fCutPIDUseAntiProtonOnly && species == kProton && track->Charge() == 1) { continue; }

    if(!fProcessSpec[species]) { continue; }

    fhPIDCounter->Fill("Selected",1);
    fhPIDCounter->Fill(GetSpeciesName(species),1);

    fVector[species]->push_back(track);
    if(fFillQA && species != kCharUnidentified) { FillQAPID(kAfter,track,species); } // filling QA for tracks AFTER selection }

    if(fProcessSpec[kPion] && fProcessSpec[kKaon] && fProcessSpec[kProton]) { // NB: aka process PID (not just Kaons for Phi)
      if(fAnalType != kMC && !FillFlowWeight(track, species)) { AliFatal("Flow weight filling failed!"); return; }
    }

      if(fMC) {
      fh2MCPtEtaReco[species]->Fill(track->Pt(), track->Eta());
      if(fAnalType != kMC && CheckMCTruthReco(species,track)) { fh2MCPtEtaRecoTrue[species]->Fill(track->Pt(), track->Eta()); }
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
AliAnalysisTaskUniFlow::PartSpecies AliAnalysisTaskUniFlow::IsPIDSelected(AliVParticle* tr) const
{
  // Selection of PID tracks (pi,K,p) - track identification
  // Based on fCutUseBayesPID flag, either Bayes PID or nSigma cutting is used
  // returns AliAnalysisTaskUniFlow::PartSpecies enum : kPion, kKaon, kProton if any of this passed kUnknown otherwise
  // *************************************************************
  AliAODTrack* track = dynamic_cast<AliAODTrack*>(tr);
  if(!track) {AliError("AOD track not found!"); return kUnknown; }

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
    dNumSigmaTOF[iSpec] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::EParticleType(iSpec)));
    //2018 data PID correction
    if(!fNeedPIDCorrection){
      dNumSigmaTPC[iSpec] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::EParticleType(iSpec)));
    }
    else{
      dNumSigmaTPC[iSpec] = TMath::Abs(PIDCorrection(track, PartSpecies(iSpec)));
    }
  }

  if(fCutUseBayesPID && !fNeedPIDCorrection) {
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
AliAnalysisTaskUniFlow::PartSpecies AliAnalysisTaskUniFlow::IsPIDSelectedMC(AliVParticle* tr) const
{
  // Selection of PID tracks (pi,K,p) - track identification
  // returns AliAnalysisTaskUniFlow::PartSpecies enum : kPion, kKaon, kProton if any of this passed kUnknown otherwise
  // *************************************************************
  if(!tr) {AliError("AliMCParticle not found!"); return kUnknown; }

  Int_t id = TMath::Abs(tr->PdgCode());

  if(id == 211) { return kPion; }
  else if(id == 321) { return kKaon; }
  else if(id == 2212) { return kProton; }
  else { return kUnknown; }
}
// ============================================================================
void AliAnalysisTaskUniFlow::FillQAPID(const QAindex iQAindex, AliVParticle* track, const PartSpecies species) const
{
  // Filling various QA plots related to PID (pi,K,p) track selection
  // *************************************************************
  if(!track) { return; }

  Int_t iPID = species - 2; // NB: translation from PartSpecies to PID QA index

  if(fAnalType == kMC) {
    AliMCParticle* tr = dynamic_cast<AliMCParticle*>(track);
    if(!tr) { AliError("AliMCParticle not found!"); return; }

    fhPIDPt[iPID]->Fill(tr->Pt());
    fhPIDPhi[iPID]->Fill(tr->Phi());
    fhPIDEta[iPID]->Fill(tr->Eta());
    fhPIDCharge[iPID]->Fill(tr->Charge()/3.);

    return;
  }

  if(!fPIDResponse || !fPIDCombined) { AliError("AliPIDResponse or AliPIDCombined object not found!"); return; }

  AliAODTrack* tr = dynamic_cast<AliAODTrack*>(track);
  if(!tr) { AliError("AliAODTrack not found!"); return; }

  // TPC & TOF statuses & measures
  AliPIDResponse::EDetPidStatus pidStatusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, tr);
  AliPIDResponse::EDetPidStatus pidStatusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, tr);

  Bool_t bIsTPCok = HasTrackPIDTPC(tr);
  Bool_t bIsTOFok = HasTrackPIDTOF(tr);

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

  Double_t dP = tr->P();
  Double_t dPt = tr->Pt();

  // TPC dEdx
  Double_t dTPCdEdx = tr->GetTPCsignal();

  // TOF beta
  Double_t dTOF[5];
  tr->GetIntegratedTimes(dTOF);
  Double_t dTOFbeta = dTOF[0] / tr->GetTOFsignal();

  // filling Bayesian PID probabilities to dBayesProb array
  UInt_t iDetUsed = fPIDCombined->ComputeProbabilities(tr, fPIDResponse, dBayesProb);

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
      fhQAPIDTPCdEdx[iQAindex]->Fill(tr->P(), dTPCdEdx);

      for(Int_t iSpec(0); iSpec < fPIDNumSpecies; ++iSpec) {
        if(!fNeedPIDCorrection){
          dNumSigmaTPC[iSpec] = fPIDResponse->NumberOfSigmasTPC(tr, AliPID::EParticleType(iSpec));
        }
        else{
          dNumSigmaTPC[iSpec] = TMath::Abs(PIDCorrection(tr, PartSpecies(iSpec)));
        }
        dTPCdEdxDelta[iSpec] = fPIDResponse->GetSignalDelta(AliPIDResponse::kTPC, tr, AliPID::EParticleType(iSpec));
      }
  }
  if(iQAindex == 0 || bUsedTOF) {
      fhQAPIDTOFstatus[iQAindex]->Fill((Int_t) pidStatusTOF );
      fhQAPIDTOFbeta[iQAindex]->Fill(dP,dTOFbeta);

      for(Int_t iSpec(0); iSpec < fPIDNumSpecies; ++iSpec) {
        dNumSigmaTOF[iSpec] = fPIDResponse->NumberOfSigmasTOF(tr, AliPID::EParticleType(iSpec));
        dTOFbetaDelta[iSpec] = fPIDResponse->GetSignalDelta(AliPIDResponse::kTOF, tr, AliPID::EParticleType(iSpec));
      }
  }

  fh3QAPIDnSigmaTPCTOFPtPion[iQAindex]->Fill(dNumSigmaTPC[2],dNumSigmaTOF[2],tr->Pt());
  fh3QAPIDnSigmaTPCTOFPtKaon[iQAindex]->Fill(dNumSigmaTPC[3],dNumSigmaTOF[3],tr->Pt());
  fh3QAPIDnSigmaTPCTOFPtProton[iQAindex]->Fill(dNumSigmaTPC[4],dNumSigmaTOF[4],tr->Pt());

  if(species == kUnknown) { return; }

  // Here only selected particles (and iQAindex == 1 by construction)

  fhPIDPt[iPID]->Fill(tr->Pt());
  fhPIDPhi[iPID]->Fill(tr->Phi());
  fhPIDEta[iPID]->Fill(tr->Eta());
  fhPIDCharge[iPID]->Fill(tr->Charge());

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
Bool_t AliAnalysisTaskUniFlow::FillCorrelations()
{
  if(!fCorrFill) { return kTRUE; }

  fSelectedTracks = new TObjArray();

  Double_t fillingCorr[4];
  fillingCorr[3] = fIndexCentrality;

  Int_t nRefs = fVector[kRefs]->size();

  for (auto part = fVector[kRefs]->begin(); part != fVector[kRefs]->end(); part++)
  {
    AliAODTrack* track = dynamic_cast<AliAODTrack*>(*part);
    if(!track) AliError("Track was not dynamically recasted.");
    Double_t etaTrig = track->Eta();
    Double_t phiTrig = track->Phi();
    // Double_t ptTrig = track->Pt();
    Int_t dTrigID = track->GetID();

    fSelectedTracks->Add(track);

    for (auto part2 = fVector[kRefs]->begin(); part2 != fVector[kRefs]->end(); part2++)
    {
      AliAODTrack* track2 = dynamic_cast<AliAODTrack*>(*part2);
      if(!track2) AliError("Track was not dynamically recasted.");
      if(track2->GetID() == dTrigID) continue;

      Double_t etaAs = track2->Eta();
      Double_t phiAs = track2->Phi();

      fillingCorr[0] = etaTrig - etaAs;
      fillingCorr[1] = RangePhi(phiTrig - phiAs);
      fillingCorr[2] = 1.0;

      fh4CorrelationsSE[kRefs]->Fill(fillingCorr);
    }
    /*
    for(Int_t iSpec(1); iSpec < kUnknown; iSpec++){
      if(!fProcessSpec[iSpec]) continue;
      if(!fh4CorrelationsSE[iSpec]) {AliError("Sparse (same event) doesn't exist."); return kFALSE; }

      if(iSpec == kCharged || iSpec == kPion || iSpec == kKaon || iSpec == kProton){
        for (auto partAs = fVector[iSpec]->begin(); partAs != fVector[iSpec]->end(); partAs++)
        {
          AliAODTrack* trackAs = dynamic_cast<AliAODTrack*>(*partAs);
          Double_t ptAs = trackAs->Phi();
          if(ptTrig > ptAs) continue;
          if(dTrigID == trackAs->GetID()) continue;

          Double_t etaAs = trackAs->Eta();
          Double_t phiAs = trackAs->Phi();

          fillingCorr[0] = etaTrig - etaAs;
          fillingCorr[1] = RangePhi(phiTrig - phiAs);
          fillingCorr[2] = ptAs;

          fh4CorrelationsSE[kCharged]->Fill(fillingCorr);
        }// end loop particle vector
      } //end direct species
      else{
        AliWarning("Not implemented yet!"); return kFALSE;
      } //end reconstructed species
    }
    */
  } // end loop reference particles (triggers)

  if(!fFillMixed) return kTRUE;

  fPool = fEventPoolMgr->GetEventPool(fIndexCentrality, fPVz);
  if(!fPool) {  AliFatal(Form("No pool found for centrality = %d, zVtx = %f", fIndexCentrality,fPVz)); return kFALSE; }

  //mixed events
  Int_t nMixEvents = fPool->GetCurrentNEvents();
  if(fPool->IsReady() || fPool->NTracksInPool() > 0.1*fMixingTracks ||  nMixEvents >= fMinEventsToMix){
    for(Int_t iMix(0); iMix < nMixEvents; iMix++){
      TObjArray *mixedEvent = fPool->GetEvent(iMix);
      if(!mixedEvent) {  AliFatal("Mixed event not found!"); return kFALSE; }

      for(auto part = fVector[kRefs]->begin(); part != fVector[kRefs]->end(); part++)
      {
        AliAODTrack* track = dynamic_cast<AliAODTrack*>(*part);
        if(!track) AliError("Track was not dynamically recasted.");
        Double_t etaTrig = track->Eta();
        Double_t phiTrig = track->Phi();

        for(Int_t mixTrack(0); mixTrack < mixedEvent->GetEntriesFast(); mixTrack++){
          AliAODTrack* track2 = (AliAODTrack*) mixedEvent->At(mixTrack);
          if(!track2) AliError("Mixed track not here!");

          Double_t etaAs = track2->Eta();
          Double_t phiAs = track2->Phi();

          fillingCorr[0] = etaTrig - etaAs;
          fillingCorr[1] = RangePhi(phiTrig - phiAs);
          fillingCorr[2] = 1.0;

          fh4CorrelationsME[kRefs]->Fill(fillingCorr);
        } //end mixed tracks
      } // end refs (current event)
    } // end loop mixed events
  } // end pool is ready (etc.)

  TObjArray* cloneArray = (TObjArray *)fSelectedTracks->Clone();
  cloneArray->SetOwner(kTRUE);
  fPool->UpdatePool(cloneArray);

  return kTRUE;
}
// ============================================================================
Double_t AliAnalysisTaskUniFlow::RangePhi(Double_t dPhi){
    if (dPhi < -0.5*TMath::Pi()) dPhi += 2 * TMath::Pi();
    if (dPhi > 1.5*TMath::Pi()) dPhi -= 2*TMath::Pi();
    return dPhi;
}
// ============================================================================

Bool_t AliAnalysisTaskUniFlow::ProcessCorrTask(const AliUniFlowCorrTask* task, const Int_t iTask, Bool_t doLowerOrder)
{
    if(!task) { AliError("AliUniFlowCorrTask does not exists!"); return kFALSE; }
    // task->PrintTask();

    Int_t iNumHarm = task->fiNumHarm;
    Int_t iNumGaps = task->fiNumGaps;

    if(iNumGaps > 2) { AliError("Too many gaps! Not implemented yet!"); return kFALSE; }
    if(iNumGaps == 2 && task->fdGaps[0] != task->fdGaps[1]) { AliError("Different position of the border when using 3 subevents! Not implemented yet!"); return kFALSE; }
    if(iNumHarm > 16) { AliError("Too many harmonics! Not implemented yet!"); return kFALSE; }

    Double_t dGap = -1.0;
    if(iNumGaps > 0) { dGap = task->fdGaps[0]; }

    // Fill anyway -> needed for any correlations
    FillRefsVectors(task, dGap);

    for(Int_t iSpec(0); iSpec < kUnknown; ++iSpec) {
        AliDebug(2,Form("Processing species '%s'",GetSpeciesName(PartSpecies(iSpec))));

        if(iSpec == kRefs) {
            if(!task->fbDoRefs) { continue; }
            CalculateCorrelations(task, kRefs);
            if(doLowerOrder){
              if(iNumHarm > 2) CalculateCorrelations(fVecCorrTask.at(iTask-1), kRefs);
              if(iNumHarm > 4) CalculateCorrelations(fVecCorrTask.at(iTask-2), kRefs);
              if(iNumHarm > 6) CalculateCorrelations(fVecCorrTask.at(iTask-3), kRefs);
              if(iNumHarm > 8) CalculateCorrelations(fVecCorrTask.at(iTask-4), kRefs);
              if(iNumHarm > 10) CalculateCorrelations(fVecCorrTask.at(iTask-5), kRefs);
              if(iNumHarm > 12) CalculateCorrelations(fVecCorrTask.at(iTask-6), kRefs);
              if(iNumHarm > 14) CalculateCorrelations(fVecCorrTask.at(iTask-7), kRefs);
            }
            if(fCorrUsingGF) CalculateDihCorr(task);
            continue;
        }

        // here-after only POIs survive (Refs are dealt with already)
        if(!task->fbDoPOIs) { continue; }
        if(!fProcessSpec[iSpec]) { continue; }

        if(iSpec == kCharUnidentified) { continue; }

        if(fPIDonlyForRefs && (iSpec == kPion || iSpec == kKaon || iSpec == kProton)) { continue; }

        if(iNumHarm > 4) { AliError("Too many harmonics! Not implemented yet!"); return kFALSE; }

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
        std::array<Int_t, 4> indexesStart = {0, 0, 0, 0};

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
                Int_t iFilledHere = 0;
                if(iSpec == kCharged && fFlowUsePIDWeights) iFilledHere = FillPOIsVectorsCharged(task, dGap, dPtLow, dPtHigh, indexesStart);
                else iFilledHere = FillPOIsVectors(task, dGap ,PartSpecies(iSpec), contIndexStart, iNumInPtBin, dPtLow, dPtHigh, dMassLow, dMassHigh);
                CalculateCorrelations(task, PartSpecies(iSpec),dPt,dMass);
                if(doLowerOrder)
                {
                  if(iNumHarm > 2) CalculateCorrelations(fVecCorrTask.at(iTask-1), PartSpecies(iSpec),dPt,dMass);
                  // if(iNumHarm > 4) CalculateCorrelations(fVecCorrTask.at(iTask-2), PartSpecies(iSpec),dPt,dMass);
                  // if(iNumHarm > 6) CalculateCorrelations(fVecCorrTask.at(iTask-3), PartSpecies(iSpec),dPt,dMass);
                }

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
void AliAnalysisTaskUniFlow::CalculateDihCorr(const AliUniFlowCorrTask* const task) const
{
  //opening all histos
  if(task->fiHarm[0] != -task->fiHarm[1]) { AliError("Harmonics are not the same! Not possible to calculate dihadron correlations."); return;}
  Int_t harm = (Int_t) task->fiHarm[0];
  Double_t width = 2.0*fFlowEtaMax/fFlowBinNumberEtaSlices;

  TComplex cNom = TComplex(0.0,0.0,kFALSE);
  TComplex cDenom = TComplex(0.0,0.0,kFALSE);
  TComplex hermConj = TComplex(0.0,0.0,kFALSE);

  TProfile* prof = nullptr;
  for(Int_t iBin(0); iBin < fFlowBinNumberEtaSlices; iBin++){
    hermConj = TComplex::Conjugate(fFlowVecQ[iBin][harm][1]);
    cDenom = fFlowVecQ[iBin][0][1]*fFlowVecQ[iBin][0][1] - fFlowVecQ[iBin][0][2];
    cNom = fFlowVecQ[iBin][harm][1]*hermConj - fFlowVecQ[iBin][0][2];

    Double_t dNom = cNom.Re();
    Double_t dDenom = cDenom.Re();
    Double_t dValue = 0.0;
    if(dDenom > 0.0) { dValue = dNom / dDenom; }

    prof = (TProfile*) fListFlow[kRefs]->FindObject(Form("%s_eta_%.3g_%.3g_sample%d",task->fsName.Data(),0.0,width,fIndexSampling));
    if(!prof) { AliError(Form("Profile '%s_eta_%.3g_%.3g_sample%d' not found!", task->fsName.Data(),0.0,width,fIndexSampling)); return; }
    prof->Fill(fIndexCentrality, dValue, dDenom);
  }

  for(Int_t iBin(0); iBin < fFlowBinNumberEtaSlices; iBin++){
    for(Int_t iSecondBin(iBin+1); iSecondBin < fFlowBinNumberEtaSlices; iSecondBin++){
        Int_t diff = iSecondBin - iBin;
        hermConj = TComplex::Conjugate(fFlowVecQ[iSecondBin][harm][1]);
        cDenom = fFlowVecQ[iBin][0][1]*fFlowVecQ[iSecondBin][0][1];
        cNom = fFlowVecQ[iBin][harm][1]*hermConj;

        Double_t dNom = cNom.Re();
        Double_t dDenom = cDenom.Re();
        Double_t dValue = 0.0;
        if(dDenom > 0.0) { dValue = dNom / dDenom; }

        prof = (TProfile*) fListFlow[kRefs]->FindObject(Form("%s_eta_%.3g_%.3g_sample%d",task->fsName.Data(),width*diff,width*(diff+1),fIndexSampling));
        if(!prof) { AliError(Form("Profile '%s_eta_%.3g_%.3g_sample%d' not found!", task->fsName.Data(),width*diff,width*(diff+1),fIndexSampling)); return; }
        prof->Fill(fIndexCentrality, dValue, dDenom);
    } // end second
  } // end iBin

  return;
}
// ============================================================================
void AliAnalysisTaskUniFlow::CalculateCorrelations(const AliUniFlowCorrTask* const task, const PartSpecies species, const Double_t dPt, const Double_t dMass) const
{
  if(!task) { AliError("AliUniFlowCorrTask does not exists!"); return; }
  if(species >= kUnknown) { AliError(Form("Invalid species: %s!", GetSpeciesName(species))); return; }

  Bool_t bHasGap = task->HasGap();
  Bool_t bHas3sub = kFALSE;
  if(task->fiNumGaps > 1) bHas3sub = kTRUE;
  Int_t iNumHarm = task->fiNumHarm;
  Bool_t bDiff = kTRUE;
  Bool_t etaCheck = kFALSE;
  if(species == kRefs) {
    bDiff = kFALSE;
    if(fEtaCheckRFP) etaCheck = kTRUE;
  }
  char sides[] = "LMR";


  // results of correlations
  TComplex cNom = TComplex(0.0,0.0,kFALSE);
  TComplex cDenom = TComplex(0.0,0.0,kFALSE);
  TComplex cNomNeg = TComplex(0.0,0.0,kFALSE);
  TComplex cDenomNeg = TComplex(0.0,0.0,kFALSE);
  TComplex cNom3Sub[3][3];
  TComplex cDenom3Sub[3][3];
  if(bHas3sub && species == kCharged){
    for(Int_t poiPos(0); poiPos < 3; poiPos++)
      for(Int_t twoPos(0); twoPos < 3; twoPos++){
        cNom3Sub[poiPos][twoPos] = TComplex(0.0,0.0,kFALSE);
        cDenom3Sub[poiPos][twoPos] = TComplex(0.0,0.0,kFALSE);
      }
  }

  if(!fUseGeneralFormula){
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
            if(bHas3sub){
              for(Int_t poiPos(0); poiPos < 3; poiPos++)
                for(Int_t rfPos(0); rfPos < 3; rfPos++){
                  if(poiPos == rfPos) continue;
                  cDenom3Sub[poiPos][rfPos] = TwoDiffGap3sub(0,0,poiPos,rfPos);
                  cNom3Sub[poiPos][rfPos] = TwoDiffGap3sub(task->fiHarm[0],task->fiHarm[1],poiPos,rfPos);
                }
            }
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
            if(bHas3sub){
              for(Int_t rf1Pos(0); rf1Pos < 3; rf1Pos++){
                if(rf1Pos > 1) break;
                for(Int_t rf2Pos(0); rf2Pos < 3; rf2Pos++ ){
                    if(rf1Pos >= rf2Pos) continue;
                    cDenom3Sub[rf1Pos][rf2Pos] = TwoGap3sub(0,0,rf1Pos,rf2Pos);
                    cNom3Sub[rf1Pos][rf2Pos] = TwoGap3sub(task->fiHarm[0],task->fiHarm[1],rf1Pos,rf2Pos);
                  }
                }
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
            if(bHas3sub){
              for(Int_t poiPos(0); poiPos < 3; poiPos++)
                for(Int_t twoPos(0); twoPos < 3; twoPos++){
                  cDenom3Sub[poiPos][twoPos] = FourDiff3sub(0,0,0,0,poiPos,twoPos);
                  cNom3Sub[poiPos][twoPos] = FourDiff3sub(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3],poiPos,twoPos);
                }
            }
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
            if(bHas3sub){
              for(Int_t twoPos(0); twoPos < 3; twoPos++){
                cDenom3Sub[0][twoPos] = Four3sub(0,0,0,0,twoPos);
                cNom3Sub[0][twoPos] = Four3sub(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3],twoPos);
              }
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
    } // end switch
  } // end not general formula
  else{
    if(bDiff) {
      AliWarning("Differential correlations using general formula not implemented!");
      return;
    }
    if(bHasGap){
      AliWarning("Correlations with eta gap using general formula not implemented!");
      return;
    }
    std::vector<Int_t> vecHarm = task->fiHarm;
    std::vector<Int_t> vecZeros(iNumHarm,0);
    Int_t* harmAr = vecHarm.data();
    Int_t* harmZerosAr = vecZeros.data();

    cDenom = Correlator(iNumHarm, harmZerosAr);
    cNom = Correlator(iNumHarm, harmAr);
  } // end general formula


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

  Bool_t bFill3sub[3][3] = {kFALSE};
  Double_t dValue3Sub[3][3] = {0.0};
  Double_t dDenom3Sub[3][3] = {0.0};
  Double_t dNom3Sub[3][3] = {0.0};

  if(bHas3sub){
    for(Int_t poiPos(0); poiPos < 3; poiPos++){
      if(species == kRefs && ( (iNumHarm == 2 && poiPos > 1) || (iNumHarm == 4 && poiPos > 0) ) ) break;
      for(Int_t twoPos(0); twoPos < 3; twoPos++){
        if(species == kRefs && iNumHarm == 2 && poiPos >= twoPos ) continue;
        if(species != kRefs && iNumHarm == 2 && poiPos == twoPos ) continue;
        bFill3sub[poiPos][twoPos] = kTRUE;
        dNom3Sub[poiPos][twoPos] = cNom3Sub[poiPos][twoPos].Re();
        dDenom3Sub[poiPos][twoPos] = cDenom3Sub[poiPos][twoPos].Re();
        if(dDenom3Sub[poiPos][twoPos] > 0.0)
          dValue3Sub[poiPos][twoPos] = dNom3Sub[poiPos][twoPos] / dDenom3Sub[poiPos][twoPos];
        else
          bFill3sub[poiPos][twoPos] = kFALSE;
        if(TMath::Abs(dValue3Sub[poiPos][twoPos]) > 1.0)
          bFill3sub[poiPos][twoPos] = kFALSE;
      }
    }
  }

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
      if(bHas3sub)
      {
        if(iNumHarm == 2){
          for(Int_t rf1Pos(0); rf1Pos < 3; rf1Pos++){
            if(rf1Pos > 1) break;
            for(Int_t rf2Pos(0); rf2Pos < 3; rf2Pos++){
              if(rf1Pos >= rf2Pos) continue;
              TProfile* prof = (TProfile*) fListFlow[species]->FindObject(Form("%s_Pos_sample%d_rf1_%c_rf2_%c",task->fsName.Data(),fIndexSampling,sides[rf1Pos],sides[rf2Pos]));
              if(!prof) { AliError(Form("Profile '%s_Pos_sample%d_rf1_%c_rf2_%c' not found!", task->fsName.Data(),fIndexSampling,sides[rf1Pos],sides[rf2Pos])); return; }
              if(bFill3sub[rf1Pos][rf2Pos]) prof->Fill(fIndexCentrality, dValue3Sub[rf1Pos][rf2Pos], dDenom3Sub[rf1Pos][rf2Pos]);
            }
          }
        }
        else if(iNumHarm == 4){
          for(Int_t twoPos(0); twoPos < 3; twoPos++){
            TProfile* prof = (TProfile*) fListFlow[species]->FindObject(Form("%s_Pos_sample%d_two_%c",task->fsName.Data(),fIndexSampling,sides[twoPos]));
            if(!prof) { AliError(Form("Profile '%s_Pos_sample%d_two_%c' not found!", task->fsName.Data(),fIndexSampling,sides[twoPos])); return; }
            if(bFill3sub[0][twoPos]) prof->Fill(fIndexCentrality, dValue3Sub[0][twoPos], dDenom3Sub[0][twoPos]);
          }
        }
        else
          return;
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

      if(bHas3sub)
      {
        if(iNumHarm == 2){
          for(Int_t poiPos(0); poiPos < 3; poiPos++)
            for(Int_t rfPos(0); rfPos < 3; rfPos++){
              if(poiPos == rfPos) continue;
              TProfile2D* prof = (TProfile2D*) fListFlow[species]->FindObject(Form("%s_Pos_sample%d_poi_%c_rfp_%c",task->fsName.Data(),fIndexSampling,sides[poiPos],sides[rfPos]));
              if(!prof) { AliError(Form("Profile '%s_Pos_sample%d_poi_%c_rfp_%c' not found!", task->fsName.Data(),fIndexSampling,sides[poiPos],sides[rfPos])); return; }
              if(bFill3sub[poiPos][rfPos]) prof->Fill(fIndexCentrality, dPt, dValue3Sub[poiPos][rfPos], dDenom3Sub[poiPos][rfPos]);
            }
          }
        else if(iNumHarm == 4){
          for(Int_t poiPos(0); poiPos < 3; poiPos++)
            for(Int_t twoPos(0); twoPos < 3; twoPos++){
              TProfile2D* prof = (TProfile2D*) fListFlow[species]->FindObject(Form("%s_Pos_sample%d_poi_%c_two_%c",task->fsName.Data(),fIndexSampling,sides[poiPos],sides[twoPos]));
              if(!prof) { AliError(Form("Profile '%s_Pos_sample%d_poi_%c_two_%c' not found!", task->fsName.Data(),fIndexSampling,sides[poiPos],sides[twoPos])); return; }
              if(bFill3sub[poiPos][twoPos]) prof->Fill(fIndexCentrality, dPt, dValue3Sub[poiPos][twoPos], dDenom3Sub[poiPos][twoPos]);
            }
        }
        else
          return;
      }
      break;
    }


    case kK0s:
    case kLambda:
    case kPhi:
    {
      if(bFillPos)
      {
        TProfile3D* prof = (TProfile3D*) fListFlow[species]->FindObject(Form("%s_Pos_sample%d",task->fsName.Data(),fIndexSampling));
        if(!prof) { AliError(Form("Profile '%s_Pos_sample%d' not found!", task->fsName.Data(),fIndexSampling)); return; }
        prof->Fill(fIndexCentrality, dPt, dMass, dValue, dDenom);
      }

      if(bFillNeg)
      {
        TProfile3D* profNeg = (TProfile3D*) fListFlow[species]->FindObject(Form("%s_Neg_sample%d",task->fsName.Data(),fIndexSampling));
        if(!profNeg) { AliError(Form("Profile '%s_Neg_sample%d' not found!", task->fsName.Data(),fIndexSampling)); return; }
        profNeg->Fill(fIndexCentrality, dPt, dMass, dValueNeg, dDenomNeg);
      }
      if(bHas3sub)
      {
        if(iNumHarm == 2){
          for(Int_t poiPos(0); poiPos < 3; poiPos++)
            for(Int_t rfPos(0); rfPos < 3; rfPos++){
              if(poiPos == rfPos) continue;
              TProfile3D* prof = (TProfile3D*) fListFlow[species]->FindObject(Form("%s_Pos_sample%d_poi_%c_rfp_%c",task->fsName.Data(),fIndexSampling,sides[poiPos],sides[rfPos]));
              if(!prof) { AliError(Form("Profile '%s_Pos_sample%d_poi_%c_rfp_%c' not found!", task->fsName.Data(),fIndexSampling,sides[poiPos],sides[rfPos])); return; }
              if(bFill3sub[poiPos][rfPos]) prof->Fill(fIndexCentrality, dPt, dMass, dValue3Sub[poiPos][rfPos], dDenom3Sub[poiPos][rfPos]);
            }
          }
        else if(iNumHarm == 4){
          for(Int_t poiPos(0); poiPos < 3; poiPos++)
            for(Int_t twoPos(0); twoPos < 3; twoPos++){
              TProfile3D* prof = (TProfile3D*) fListFlow[species]->FindObject(Form("%s_Pos_sample%d_poi_%c_two_%c",task->fsName.Data(),fIndexSampling,sides[poiPos],sides[twoPos]));
              if(!prof) { AliError(Form("Profile '%s_Pos_sample%d_poi_%c_two_%c' not found!", task->fsName.Data(),fIndexSampling,sides[poiPos],sides[twoPos])); return; }
              if(bFill3sub[poiPos][twoPos]) prof->Fill(fIndexCentrality, dPt, dMass, dValue3Sub[poiPos][twoPos], dDenom3Sub[poiPos][twoPos]);
            }
        }
        else
          return;
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
  Bool_t doLowerOrder = kFALSE;
  for(Int_t iTask(0); iTask < iNumTasks; ++iTask)
  {
    //if you want to save some CPU time
    //you need to add tasks with same gaps one after the other
    //starting with the lowest harmonics
    Bool_t isLowerOrder = kFALSE;
    AliUniFlowCorrTask* thisTask = (AliUniFlowCorrTask*) fVecCorrTask.at(iTask);
    AliUniFlowCorrTask* nextTask = nullptr;
    if(iTask+1 < iNumTasks) {
      nextTask = (AliUniFlowCorrTask*) fVecCorrTask.at(iTask+1);
      if(thisTask->fiNumHarm < nextTask->fiNumHarm && thisTask->fiHarm[0] == nextTask->fiHarm[0] && thisTask->fiNumGaps == nextTask->fiNumGaps){
        isLowerOrder = kTRUE;
        for(Int_t gaps(0); gaps < thisTask->fiNumGaps; gaps++){
          if(thisTask->fdGaps[gaps] != nextTask->fdGaps[gaps]) isLowerOrder = kFALSE;
        }
      }
    }
    if(isLowerOrder) {
      doLowerOrder = kTRUE;
      continue;
    }
    Bool_t process = ProcessCorrTask(fVecCorrTask.at(iTask), iTask, doLowerOrder);
    if(!process) { AliError("AliUniFlowCorrTask processing failed!\n"); fVecCorrTask.at(iTask)->Print(); return kFALSE; }
    doLowerOrder = kFALSE;
  }

  fEventCounter++; // counter of processed events

  return kTRUE;
}
// ============================================================================
void AliAnalysisTaskUniFlow::FillFlowQVectorsForDih(const Double_t dWeight, const Double_t dPhi, const Double_t dEta, const Int_t harm)
{
  // Filling Q flow vector with RFPs for dihadron correlation study
  for(Int_t iBin(0); iBin < fFlowBinNumberEtaSlices; iBin++){
    if(dEta > fEtaSlicesArr[iBin+1]) continue;
    Double_t dCos = dWeight * TMath::Cos(harm * dPhi);
    Double_t dSin = dWeight * TMath::Sin(harm * dPhi);
    fFlowVecQ[iBin][harm][1] += TComplex(dCos,dSin,kFALSE);

    fFlowVecQ[iBin][0][1] += TComplex(dWeight,(Double_t) 0.0,kFALSE);

    dCos = TMath::Power(dWeight,2);
    fFlowVecQ[iBin][0][2] += TComplex(dCos,(Double_t) 0.0,kFALSE);
    break;
  }
  return;
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
  Double_t dEtaLim3sub = dEtaLimit;
  if(dEtaGap > -1.0) { bHasGap = kTRUE; }
  if(task->fiNumGaps > 1) {
    bHas3sub = kTRUE;
    if(task->fdGaps[0] != task->fdGaps[1]) dEtaLim3sub = task->fdGaps[1]/2;
  }

  Int_t maxHarm = task->fMaxHarm;
  Int_t maxWeightPower = task->fMaxWeightPower;

  Bool_t usePowVector = task->fbUsePowerVector;
  std::vector<Int_t> maxPowVec = {0};
  if(usePowVector) maxPowVec = task->fiMaxPow;

  // clearing output (global) flow vectors
  ResetFlowVector(fFlowVecQpos, maxHarm, maxWeightPower, usePowVector, maxPowVec);
  ResetFlowVector(fFlowVecQneg, maxHarm, maxWeightPower, usePowVector, maxPowVec);
  if(bHas3sub) { ResetFlowVector(fFlowVecQmid, maxHarm, maxWeightPower, usePowVector, maxPowVec); }

  if(fCorrUsingGF) { ResetFlowVectorQdih(fFlowVecQ, (Int_t) task->fiHarm[0]); }


  if(!fFlowUsePIDWeights){
    for (auto part = fVector[kRefs]->begin(); part != fVector[kRefs]->end(); part++)
    {
      Double_t dPhi = (*part)->Phi();
      Double_t dEta = (*part)->Eta();

      if(bHasGap && TMath::Abs(dEta) < dEtaLimit && !bHas3sub) { continue; }

      // loading weights if needed
      Double_t dWeight = 1.0;
      if(fFlowUseWeights) { dWeight = GetFlowWeight(*part, kRefs); }

      if(fCorrUsingGF) { FillFlowQVectorsForDih((Double_t) dWeight, (Double_t) dPhi, (Double_t) dEta, (Int_t) task->fiHarm[0]); }

      if(!bHasGap) // no eta gap
      {
        for(Int_t iHarm(0); iHarm <= maxHarm; iHarm++){
          if(usePowVector) maxWeightPower = maxPowVec[iHarm];
          for(Int_t iPower(0); iPower <= maxWeightPower; iPower++)
          {
            Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
            Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
            fFlowVecQpos[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
          }
        }
      }
      else
      {
        // RFP in positive eta acceptance
        if(dEta > dEtaLimit)
        {
          for(Int_t iHarm(0); iHarm <= maxHarm; iHarm++){
            if(usePowVector) maxWeightPower = maxPowVec[iHarm];
            for(Int_t iPower(0); iPower <= maxWeightPower; iPower++)
            {
              Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
              Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
              fFlowVecQpos[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
            }
          }
        }
        // RFP in negative eta acceptance
        if(dEta < -dEtaLimit)
        {
          for(Int_t iHarm(0); iHarm <= maxHarm; iHarm++){
            if(usePowVector) maxWeightPower = maxPowVec[iHarm];
            for(Int_t iPower(0); iPower <= maxWeightPower; iPower++)
            {
              Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
              Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
              fFlowVecQneg[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
            }
          }
        }

        // RFP in middle (for 3sub) if gap > 0
        if(bHas3sub && (TMath::Abs(dEta) < dEtaLim3sub) )
        {
          for(Int_t iHarm(0); iHarm <= maxHarm; iHarm++){
            if(usePowVector) maxWeightPower = maxPowVec[iHarm];
            for(Int_t iPower(0); iPower <= maxWeightPower; iPower++)
            {
              Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
              Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
              fFlowVecQmid[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
            }
          }
        }

      } // endif {dEtaGap}
    } // endfor {tracks} particle loop

    return;
  } //end fFlowUsePIDWeights

  //if use PID weights

  for (auto part = fVector[kCharUnidentified]->begin(); part != fVector[kCharUnidentified]->end(); part++)
  {
    Double_t dPhi = (*part)->Phi();
    Double_t dEta = (*part)->Eta();

    if(!IsWithinRefs(*part)) { continue; }

    if(bHasGap && TMath::Abs(dEta) < dEtaLimit && !bHas3sub) { continue; }

    // loading weights if needed
    Double_t dWeight = 1.0;
    if(fFlowUseWeights) { dWeight = GetFlowWeight(*part, kCharUnidentified); }

    if(!bHasGap) // no eta gap
    {
      for(Int_t iHarm(0); iHarm <= maxHarm; iHarm++){
        if(usePowVector) maxWeightPower = maxPowVec[iHarm];
        for(Int_t iPower(0); iPower <= maxWeightPower; iPower++)
        {
          Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
          Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
          fFlowVecQpos[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
        }
      }
    }
    else
    {
      // RFP in positive eta acceptance
      if(dEta > dEtaLimit)
      {
        for(Int_t iHarm(0); iHarm <= maxHarm; iHarm++){
          if(usePowVector) maxWeightPower = maxPowVec[iHarm];
          for(Int_t iPower(0); iPower <= maxWeightPower; iPower++)
          {
            Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
            Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
            fFlowVecQpos[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
          }
        }
      }
      // RFP in negative eta acceptance
      if(dEta < -dEtaLimit)
      {
        for(Int_t iHarm(0); iHarm <= maxHarm; iHarm++){
          if(usePowVector) maxWeightPower = maxPowVec[iHarm];
          for(Int_t iPower(0); iPower <= maxWeightPower; iPower++)
          {
            Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
            Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
            fFlowVecQneg[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
          }
        }
      }

      // RFP in middle (for 3sub) if gap > 0
      if(bHas3sub && (TMath::Abs(dEta) < dEtaLim3sub) )
      {
        for(Int_t iHarm(0); iHarm <= maxHarm; iHarm++){
          if(usePowVector) maxWeightPower = maxPowVec[iHarm];
          for(Int_t iPower(0); iPower <= maxWeightPower; iPower++)
          {
            Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
            Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
            fFlowVecQmid[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
          }
        }
      }

    } // endif {dEtaGap}
  } // endfor {tracks} particle loop

  for (auto part = fVector[kPion]->begin(); part != fVector[kPion]->end(); part++)
  {
    Double_t dPhi = (*part)->Phi();
    Double_t dEta = (*part)->Eta();

    if(!IsWithinRefs(*part)) { continue; }

    if(bHasGap && TMath::Abs(dEta) < dEtaLimit && !bHas3sub) { continue; }

    // loading weights if needed
    Double_t dWeight = 1.0;
    if(fFlowUseWeights) { dWeight = GetFlowWeight(*part, kPion); }

    if(!bHasGap) // no eta gap
    {
      for(Int_t iHarm(0); iHarm <= maxHarm; iHarm++){
        if(usePowVector) maxWeightPower = maxPowVec[iHarm];
        for(Int_t iPower(0); iPower <= maxWeightPower; iPower++)
        {
          Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
          Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
          fFlowVecQpos[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
        }
      }
    }
    else
    {
      // RFP in positive eta acceptance
      if(dEta > dEtaLimit)
      {
        for(Int_t iHarm(0); iHarm <= maxHarm; iHarm++){
          if(usePowVector) maxWeightPower = maxPowVec[iHarm];
          for(Int_t iPower(0); iPower <= maxWeightPower; iPower++)
          {
            Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
            Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
            fFlowVecQpos[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
          }
        }
      }
      // RFP in negative eta acceptance
      if(dEta < -dEtaLimit)
      {
        for(Int_t iHarm(0); iHarm <= maxHarm; iHarm++){
          if(usePowVector) maxWeightPower = maxPowVec[iHarm];
          for(Int_t iPower(0); iPower <= maxWeightPower; iPower++)
          {
            Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
            Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
            fFlowVecQneg[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
          }
        }
      }

      // RFP in middle (for 3sub) if gap > 0
      if(bHas3sub && (TMath::Abs(dEta) < dEtaLim3sub) )
      {
        for(Int_t iHarm(0); iHarm <= maxHarm; iHarm++){
          if(usePowVector) maxWeightPower = maxPowVec[iHarm];
          for(Int_t iPower(0); iPower <= maxWeightPower; iPower++)
          {
            Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
            Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
            fFlowVecQmid[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
          }
        }
      }

    } // endif {dEtaGap}
  } // endfor {tracks} particle loop

    for (auto part = fVector[kProton]->begin(); part != fVector[kProton]->end(); part++)
  {
    Double_t dPhi = (*part)->Phi();
    Double_t dEta = (*part)->Eta();

    if(!IsWithinRefs(*part)) { continue; }

    if(bHasGap && TMath::Abs(dEta) < dEtaLimit && !bHas3sub) { continue; }

    // loading weights if needed
    Double_t dWeight = 1.0;
    if(fFlowUseWeights) { dWeight = GetFlowWeight(*part, kProton); }

    if(!bHasGap) // no eta gap
    {
      for(Int_t iHarm(0); iHarm <= maxHarm; iHarm++){
        if(usePowVector) maxWeightPower = maxPowVec[iHarm];
        for(Int_t iPower(0); iPower <= maxWeightPower; iPower++)
        {
          Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
          Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
          fFlowVecQpos[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
        }
      }
    }
    else
    {
      // RFP in positive eta acceptance
      if(dEta > dEtaLimit)
      {
        for(Int_t iHarm(0); iHarm <= maxHarm; iHarm++){
          if(usePowVector) maxWeightPower = maxPowVec[iHarm];
          for(Int_t iPower(0); iPower <= maxWeightPower; iPower++)
          {
            Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
            Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
            fFlowVecQpos[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
          }
        }
      }
      // RFP in negative eta acceptance
      if(dEta < -dEtaLimit)
      {
        for(Int_t iHarm(0); iHarm <= maxHarm; iHarm++){
          if(usePowVector) maxWeightPower = maxPowVec[iHarm];
          for(Int_t iPower(0); iPower <= maxWeightPower; iPower++)
          {
            Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
            Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
            fFlowVecQneg[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
          }
        }
      }

      // RFP in middle (for 3sub) if gap > 0
      if(bHas3sub && (TMath::Abs(dEta) < dEtaLim3sub) )
      {
        for(Int_t iHarm(0); iHarm <= maxHarm; iHarm++){
          if(usePowVector) maxWeightPower = maxPowVec[iHarm];
          for(Int_t iPower(0); iPower <= maxWeightPower; iPower++)
          {
            Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
            Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
            fFlowVecQmid[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
          }
        }
      }

    } // endif {dEtaGap}
  } // endfor {tracks} particle loop

  if(fProcessSpec[kPhi] && !fPIDCorrectionPhi) return;

  for (auto part = fVector[kKaon]->begin(); part != fVector[kKaon]->end(); part++)
  {
    Double_t dPhi = (*part)->Phi();
    Double_t dEta = (*part)->Eta();

    if(!IsWithinRefs(*part)) { continue; }

    if(bHasGap && TMath::Abs(dEta) < dEtaLimit && !bHas3sub) { continue; }

    // loading weights if needed
    Double_t dWeight = 1.0;
    if(fFlowUseWeights) { dWeight = GetFlowWeight(*part, kKaon); }

    if(!bHasGap) // no eta gap
    {
      for(Int_t iHarm(0); iHarm <= maxHarm; iHarm++){
        if(usePowVector) maxWeightPower = maxPowVec[iHarm];
        for(Int_t iPower(0); iPower <= maxWeightPower; iPower++)
        {
          Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
          Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
          fFlowVecQpos[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
        }
      }
    }
    else
    {
      // RFP in positive eta acceptance
      if(dEta > dEtaLimit)
      {
        for(Int_t iHarm(0); iHarm <= maxHarm; iHarm++){
          if(usePowVector) maxWeightPower = maxPowVec[iHarm];
          for(Int_t iPower(0); iPower <= maxWeightPower; iPower++)
          {
            Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
            Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
            fFlowVecQpos[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
          }
        }
      }
      // RFP in negative eta acceptance
      if(dEta < -dEtaLimit)
      {
        for(Int_t iHarm(0); iHarm <= maxHarm; iHarm++){
          if(usePowVector) maxWeightPower = maxPowVec[iHarm];
          for(Int_t iPower(0); iPower <= maxWeightPower; iPower++)
          {
            Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
            Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
            fFlowVecQneg[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
          }
        }
      }

      // RFP in middle (for 3sub) if gap > 0
      if(bHas3sub && (TMath::Abs(dEta) < dEtaLim3sub) )
      {
        for(Int_t iHarm(0); iHarm <= maxHarm; iHarm++){
          if(usePowVector) maxWeightPower = maxPowVec[iHarm];
          for(Int_t iPower(0); iPower <= maxWeightPower; iPower++)
          {
            Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
            Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
            fFlowVecQmid[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
          }
        }
      }

    } // endif {dEtaGap}
  } // endfor {tracks} particle loop

  return;
}
// ============================================================================
Int_t AliAnalysisTaskUniFlow::FillPOIsVectorsCharged(const AliUniFlowCorrTask* task, const Double_t dEtaGap, const Double_t dPtLow, const Double_t dPtHigh, std::array<Int_t, 4> &indexStart)
{
  // Filling p,q and s flow vectors with POIs (given by species) for differential flow calculation
  // *************************************************************
  std::vector<AliVParticle*>* pions = fVector[kPion];
  std::vector<AliVParticle*>* kaons = fVector[kKaon];
  std::vector<AliVParticle*>* protons = fVector[kProton];
  std::vector<AliVParticle*>* unidentified = fVector[kCharUnidentified];
  if(!pions || !kaons || !protons || !unidentified) { AliError("Vector with selected POIs not found."); return 0; }

  Double_t dEtaLimit = dEtaGap / 2.0;
  Bool_t bHasGap = kFALSE; if(dEtaGap > -1.0) { bHasGap = kTRUE; }
  Bool_t bHas3sub = kFALSE; if(task->fiNumGaps > 1) { bHas3sub = kTRUE; }

  Int_t maxHarm = task->fMaxHarm;
  Int_t maxWeightPower = task->fMaxWeightPower;

  // clearing output (global) flow vectors
  ResetFlowVector(fFlowVecPpos, maxHarm, maxWeightPower);
  ResetFlowVector(fFlowVecSpos, maxHarm, maxWeightPower);

  if(bHasGap) {
    ResetFlowVector(fFlowVecPneg, maxHarm, maxWeightPower);
    ResetFlowVector(fFlowVecSneg, maxHarm, maxWeightPower);
  }

  if(bHas3sub) {
    ResetFlowVector(fFlowVecPmid, maxHarm, maxWeightPower);
    ResetFlowVector(fFlowVecSmid, maxHarm, maxWeightPower);
  }

  Int_t iTracksFilled = 0; // counter of filled tracks

  //loop over pions
  for(Int_t index(indexStart[0]); index < (Int_t) pions->size(); ++index) {
    AliVParticle* part = pions->at(index);
    if(!part) { AliError("Particle does not exists within given vector"); return -1; }

    Double_t dPhi = part->Phi();
    Double_t dEta = part->Eta();
    Double_t dPt = part->Pt();

    // checking if pt is within pt (bin) range
    if(dPt < dPtLow) { continue; }
    if(dPt >= dPtHigh) {
        // refresh the starting index value for next pt bin
        indexStart[0] = index;
        break;
    }

    if(bHasGap && TMath::Abs(dEta) < dEtaLimit) { continue; }

    iTracksFilled++;

    // loading weights if needed
    Double_t dWeight = 1.0;
    if(fFlowUseWeights) { dWeight = GetFlowWeight(part, kPion); }

    // check if POI overlaps with RFPs (not for reconstructed)
    Bool_t bIsWithinRefs = IsWithinRefs(static_cast<const AliAODTrack*>(part));

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

  //loop over Kaons
  for(Int_t index(indexStart[1]); index < (Int_t) kaons->size(); ++index) {
     AliVParticle* part = kaons->at(index);
     if(!part) { AliError("Particle does not exists within given vector"); return -1; }

     Double_t dPhi = part->Phi();
     Double_t dEta = part->Eta();
     Double_t dPt = part->Pt();

     // checking if pt is within pt (bin) range
     if(dPt < dPtLow) { continue; }
     if(dPt >= dPtHigh) {
         // refresh the starting index value for next pt bin
         indexStart[1] = index;
         break;
     }

     if(bHasGap && TMath::Abs(dEta) < dEtaLimit) { continue; }

     iTracksFilled++;

     // loading weights if needed
     Double_t dWeight = 1.0;
     if(fFlowUseWeights) { dWeight = GetFlowWeight(part, kKaon); }

     // check if POI overlaps with RFPs (not for reconstructed)
     Bool_t bIsWithinRefs = IsWithinRefs(static_cast<const AliAODTrack*>(part));

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

  //loop over protons
  for(Int_t index(indexStart[2]); index < (Int_t) protons->size(); ++index) {
      AliVParticle* part = protons->at(index);
      if(!part) { AliError("Particle does not exists within given vector"); return -1; }

      Double_t dPhi = part->Phi();
      Double_t dEta = part->Eta();
      Double_t dPt = part->Pt();

      // checking if pt is within pt (bin) range
      if(dPt < dPtLow) { continue; }
      if(dPt >= dPtHigh) {
          // refresh the starting index value for next pt bin
          indexStart[2] = index;
          break;
      }

      if(bHasGap && TMath::Abs(dEta) < dEtaLimit) { continue; }

      iTracksFilled++;

      // loading weights if needed
      Double_t dWeight = 1.0;
      if(fFlowUseWeights) { dWeight = GetFlowWeight(part, kProton); }

      // check if POI overlaps with RFPs (not for reconstructed)
      Bool_t bIsWithinRefs = IsWithinRefs(static_cast<const AliAODTrack*>(part));

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

  //loop over unID
  for(Int_t index(indexStart[3]); index < (Int_t) unidentified->size(); ++index) {
      AliVParticle* part = unidentified->at(index);
      if(!part) { AliError("Particle does not exists within given vector"); return -1; }

      Double_t dPhi = part->Phi();
      Double_t dEta = part->Eta();
      Double_t dPt = part->Pt();

      // checking if pt is within pt (bin) range
      if(dPt < dPtLow) { continue; }
      if(dPt >= dPtHigh) {
          // refresh the starting index value for next pt bin
          indexStart[3] = index;
          break;
      }

      if(bHasGap && TMath::Abs(dEta) < dEtaLimit) { continue; }

      iTracksFilled++;

      // loading weights if needed
      Double_t dWeight = 1.0;
      if(fFlowUseWeights) { dWeight = GetFlowWeight(part, kCharUnidentified); }

      // check if POI overlaps with RFPs (not for reconstructed)
      Bool_t bIsWithinRefs = IsWithinRefs(static_cast<const AliAODTrack*>(part));

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

  return iTracksFilled;
}
// ============================================================================
Int_t AliAnalysisTaskUniFlow::FillPOIsVectors(const AliUniFlowCorrTask* task, const Double_t dEtaGap, const PartSpecies species, Int_t& indStart, Int_t& tracksInBin, const Double_t dPtLow, const Double_t dPtHigh, const Double_t dMassLow, const Double_t dMassHigh)
{
  // Filling p,q and s flow vectors with POIs (given by species) for differential flow calculation
  // *************************************************************
  std::vector<AliVParticle*>* vector = fVector[species];
  if(!vector) { AliError("Vector with selected POIs not found."); return 0; }

  Double_t dEtaLimit = dEtaGap / 2.0;
  Bool_t bHasGap = kFALSE; if(dEtaGap > -1.0) { bHasGap = kTRUE; }
  Bool_t bHas3sub = kFALSE; if(task->fiNumGaps > 1) { bHas3sub = kTRUE; }
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

  if(bHas3sub) {
    ResetFlowVector(fFlowVecPmid, maxHarm, maxWeightPower);
    ResetFlowVector(fFlowVecSmid, maxHarm, maxWeightPower);
  }

  Int_t iTracksFilled = 0; // counter of filled tracks
  Int_t iTracksInPtBin = 0; // counter for all tracks in pt bins

  // for(auto part = vector->begin(); part != vector->end(); ++part)
  for(Int_t index(indStart); index < (Int_t) vector->size(); ++index) {
    AliVParticle* part = vector->at(index);
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

    if(bHasGap && TMath::Abs(dEta) < dEtaLimit && !bHas3sub) { continue; }

    // at this point particles corresponding to this pt (& mass) bin and eta acceptance (gap) survives
    iTracksFilled++;

    // loading weights if needed
    Double_t dWeight = 1.0;
    if(fFlowUseWeights) { dWeight = GetFlowWeight(part, species); }

    // Double_t dWeightRef = 1.0;
    // if(fFlowUseWeights) { dWeightRef = GetFlowWeight(part, kRefs); }

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
            Double_t dCos = 0.0;
            Double_t dSin = 0.0;
            // "mode 2" of GF correction
            // if(iPower > 1){
            //   dCos = dWeight * TMath::Power(dWeightRef,iPower-1) * TMath::Cos(iHarm * dPhi);
            //   dSin = dWeight * TMath::Power(dWeightRef,iPower-1) * TMath::Sin(iHarm * dPhi);
            // }
            // else{
            //   dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
            //   dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
            // }
            dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
            dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
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
              Double_t dCos = 0.0;
              Double_t dSin = 0.0;
              // if(iPower > 1){
              //   dCos = dWeight * TMath::Power(dWeightRef,iPower-1) * TMath::Cos(iHarm * dPhi);
              //   dSin = dWeight * TMath::Power(dWeightRef,iPower-1) * TMath::Sin(iHarm * dPhi);
              // }
              // else{
              //   dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
              //   dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
              // }
              dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
              dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
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
               Double_t dCos = 0.0;
               Double_t dSin = 0.0;
               // if(iPower > 1){
               //   dCos = dWeight * TMath::Power(dWeightRef,iPower-1) * TMath::Cos(iHarm * dPhi);
               //   dSin = dWeight * TMath::Power(dWeightRef,iPower-1) * TMath::Sin(iHarm * dPhi);
               // }
               // else{
               //   dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
               //   dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
               // }
               dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
               dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
               fFlowVecSneg[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
             }
           }
       }
       //
       if(bHas3sub && (TMath::Abs(dEta) < dEtaLimit) ) //particle in mid acceptance
       {
         for(Int_t iHarm(0); iHarm <= maxHarm; iHarm++)
           for(Int_t iPower(0); iPower <= maxWeightPower; iPower++)
           {
             Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
             Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
             fFlowVecPmid[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);

             if(bIsWithinRefs)
             {
               Double_t dCos = 0.0;
               Double_t dSin = 0.0;
               // if(iPower > 1){
               //   dCos = dWeight * TMath::Power(dWeightRef,iPower-1) * TMath::Cos(iHarm * dPhi);
               //   dSin = dWeight * TMath::Power(dWeightRef,iPower-1) * TMath::Sin(iHarm * dPhi);
               // }
               // else{
               //   dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
               //   dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
               // }
               dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
               dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
               fFlowVecSmid[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
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
void AliAnalysisTaskUniFlow::ResetFlowVector(TComplex (&array)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax], Int_t maxHarm, Int_t maxWeightPower, Bool_t usePow, std::vector<Int_t> maxPowVec)
{
  // Reset RFPs (Q) array values to TComplex(0,0,kFALSE) for given array
  // *************************************************************
  for(Int_t iHarm(0); iHarm <= maxHarm; ++iHarm) {
    if(usePow) maxWeightPower = maxPowVec[iHarm];
    for(Int_t iPower(0); iPower <= maxWeightPower; ++iPower) {
      array[iHarm][iPower](0.0,0.0);
    }
  }
  return;
}
// ============================================================================
void AliAnalysisTaskUniFlow::ResetFlowVectorQdih(TComplex (&array)[fFlowBinNumberEtaSlices][6][3], Int_t harm)
{
  // Reset RFPs (Q) array values to TComplex(0,0,kFALSE) for given array
  // *************************************************************
  for(Int_t iBin(0); iBin < fFlowBinNumberEtaSlices; iBin++){
    array[iBin][harm][1](0.0,0.0);
    array[iBin][0][1](0.0,0.0);
    array[iBin][0][2](0.0,0.0);
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
    if(fMC) {
      AliWarning("Monte Carlo run with centrality estimator different than Nch. Returning -1");
      return -1;
    }

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
void AliAnalysisTaskUniFlow::AddCorr(std::vector<Int_t> harms, std::vector<Double_t> gaps, Bool_t doRFPs, Bool_t doPOIs, std::vector<Int_t> maxPowVec)
{
    fVecCorrTask.push_back(new AliUniFlowCorrTask(doRFPs, doPOIs, harms, gaps, maxPowVec));
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
TComplex AliAnalysisTaskUniFlow::PGapMid(const Int_t n, const Int_t p) const
{
  if(n < 0) return TComplex::Conjugate(fFlowVecPmid[-n][p]);
  else return fFlowVecPmid[n][p];
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
TComplex AliAnalysisTaskUniFlow::SGapMid(const Int_t n, const Int_t p) const
{
  if(n < 0) return TComplex::Conjugate(fFlowVecSmid[-n][p]);
  else return fFlowVecSmid[n][p];
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::Correlator(Int_t n, Int_t* harmonic, Int_t mult, Int_t skip) const
{
  //general formula
  //for calculation of correlations
  Int_t n_m_1 = n-1;
  TComplex c(Q(harmonic[n_m_1], mult));
  if (n_m_1 == 0) return c;
  c *= Correlator(n_m_1, harmonic);
  if (n_m_1 == skip) return c;

  Int_t mult_p_1 = mult+1;
  Int_t n_m_2 = n-2;
  Int_t counter1 = 0;
  Int_t h_hold = harmonic[counter1];
  harmonic[counter1] = harmonic[n_m_2];
  harmonic[n_m_2] = h_hold + harmonic[n_m_1];
  TComplex c2(Correlator(n_m_1, harmonic, mult_p_1, n_m_2));
  Int_t counter2 = n-3;
  while (counter2 >= skip) {
    harmonic[n_m_2] = harmonic[counter1];
    harmonic[counter1] = h_hold;
    ++counter1;
    h_hold = harmonic[counter1];
    harmonic[counter1] = harmonic[n_m_2];
    harmonic[n_m_2] = h_hold + harmonic[n_m_1];
    c2 += Correlator(n_m_1, harmonic, mult_p_1, counter2);
    --counter2;
  }
  harmonic[n_m_2] = harmonic[counter1];
  harmonic[counter1] = h_hold;

  if (mult == 1) return c-c2;
  return c-Double_t(mult)*c2;;
}
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
TComplex AliAnalysisTaskUniFlow::TwoGap3sub(const Int_t n1, const Int_t n2, const Int_t rf1Pos, const Int_t rf2Pos) const
{
  TComplex formula = TComplex(0.0,0.0,kFALSE);
  if(rf1Pos >= rf2Pos) { AliError("TwoGap3sub: Incorrect position of RFPs."); return 0; }
  if(rf1Pos == 0){
    if(rf2Pos == 1)
      formula = QGapNeg(n1,1)*QGapMid(n2,1);
    else
      formula = QGapNeg(n1,1)*QGapPos(n2,1);
  }
  else
    formula = QGapMid(n1,1)*QGapPos(n2,1);
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
TComplex AliAnalysisTaskUniFlow::TwoMid(const Int_t n1, const Int_t n2) const
{
  TComplex formula = QGapMid(n1,1)*QGapMid(n2,1) - QGapMid(n1+n2,2);
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
TComplex AliAnalysisTaskUniFlow::TwoDiffMid(const Int_t n1, const Int_t n2) const
{
  TComplex formula = PGapMid(n1,1)*QGapMid(n2,1) - SGapMid(n1+n2,2);
  return formula;
}
// ============================================================================
TComplex AliAnalysisTaskUniFlow::TwoDiffGap3sub(const Int_t n1, const Int_t n2, const Int_t poiPos, const Int_t rfPos) const
{
  TComplex formula = TComplex(0.0,0.0,kFALSE);
  if(poiPos == rfPos) { AliError("TwoDiffGap3sub: Incorrect position of POI and RFP."); return 0; }
  switch(poiPos){
    case 0:
    {
      if(rfPos == 1) formula = PGapNeg(n1,1)*QGapMid(n2,1);
      else formula = PGapNeg(n1,1)*QGapPos(n2,1);
      break;
    }
    case 1:
    {
      if(rfPos == 0) formula = PGapMid(n1,1)*QGapNeg(n2,1);
      else formula = PGapMid(n1,1)*QGapPos(n2,1);
      break;
    }
    case 2:
    {
      if(rfPos == 0) formula = PGapPos(n1,1)*QGapNeg(n2,1);
      else formula = PGapPos(n1,1)*QGapMid(n2,1);
      break;
    }
    default:
      return 0;
  }
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
TComplex AliAnalysisTaskUniFlow::Four3sub(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t twoParCorrPosition) const
{
  // left = neg, middle = mid; rigth = pos
  TComplex formula = TComplex(0.0,0.0,kFALSE);
  if(!(n1 == n2 && n1 == -n3 && n3 == n4) ) { AliError("Four par. diff. correlation with different harmonics not implemented!"); return 0; }
  switch(twoParCorrPosition){
    case 0:
    {
      formula = TwoNeg(n1,n2)*QGapMid(n3,1)*QGapPos(n4,1);
      break;
    }
    case 1:
    {
      formula = QGapNeg(n3,1)*TwoMid(n1,n2)*QGapPos(n4,1);
      break;
    }
    case 2:
    {
      formula = QGapNeg(n3,1)*QGapMid(n4,1)*TwoPos(n1,n2);
      break;
    }
    default:
      return 0;
  }
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
TComplex AliAnalysisTaskUniFlow::FourDiff3sub(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t poiPosition, const Int_t twoParCorrPosition) const
{
  /*
  Four particle differential correlations with 3 subevents
  0 = left, 1 = middle, 2 = right subevent
  Important to distinquish the position of POI and the position of 2-pc
  POI can be within 2-pc
  POI always with n1
  2-pc (the same subevent): harmonics have to have the same sign!
  */
  TComplex formula = TComplex(0.0,0.0,kFALSE);
  if(!(n1 == n2 && n1 == -n3 && n3 == n4) ) { AliError("Four par. diff. correlation with different harmonics not implemented!"); return 0; }
  switch (poiPosition) {
    case 0:
    {
      switch (twoParCorrPosition) {
        case 0:
        {
          formula = TwoDiffNeg(n1,n2)*QGapMid(n3,1)*QGapPos(n4,1);
          break;
        }
        case 1:
        {
          formula = PGapNeg(n1,1)*TwoMid(n3,n4)*QGapPos(n2,1);
          break;
        }
        case 2:
        {
          formula = PGapNeg(n1,1)*QGapMid(n2,1)*TwoPos(n3,n4);
          break;
        }
        default:
          return 0;
      }
      break;
    }
    case 1:
    {
      switch (twoParCorrPosition) {
        case 0:
        {
          formula = TwoNeg(n3,n4)*PGapMid(n1,1)*QGapPos(n2,1);
          break;
        }
        case 1:
        {
          formula = QGapNeg(n3,1)*TwoDiffMid(n1,n2)*QGapPos(n4,1);
          break;
        }
        case 2:
        {
          formula = QGapNeg(n2,1)*PGapMid(n1,1)*TwoPos(n3,n4);
          break;
        }
        default:
          return 0;
      }
      break;
    }
    case 2:
    {
      switch (twoParCorrPosition) {
        case 0:
        {
          formula = TwoNeg(n3,n4)*QGapMid(n2,1)*PGapPos(n1,1);
          break;
        }
        case 1:
        {
          formula = QGapNeg(n2,1)*TwoMid(n3,n4)*PGapPos(n1,1);
          break;
        }
        case 2:
        {
          formula = QGapNeg(n3,1)*QGapMid(n4,1)*TwoDiffPos(n1,n2);
          break;
        }
        default:
          return 0;
      }
      break;
    }
    default:
      return 0;
  }
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
  char sides[] = "LMR";

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
      Bool_t bHas3sub = kFALSE;
      if(task->fiNumGaps > 1) bHas3sub = kTRUE;
      Int_t corrOrder = task->fiNumHarm;
      const char* corName = task->fsName.Data();
      const char* corLabel = task->fsLabel.Data();

      for(Int_t iSpec(0); iSpec < kUnknown; ++iSpec)
      {
        // check if AliUniFlowCorrTask should be done for all flow particles (RFP/POI/Both)
        if(!task->fbDoRefs && iSpec == kRefs) { continue; }
        if(!task->fbDoPOIs && iSpec != kRefs) { continue; }

        if(!fProcessSpec[iSpec]) { continue; }
        if(iSpec == kKaon && (!fProcessSpec[kPion] || !fProcessSpec[kProton])) { continue; }
        if(iSpec == kCharUnidentified) continue;
        if(fPIDonlyForRefs && (iSpec == kPion || iSpec == kKaon || iSpec == kProton)) continue;

        Int_t iNumPtFixBins = fFlowPOIsPtBinEdges[iSpec].size() - 1;
        Double_t* dPtFixBinEdges = fFlowPOIsPtBinEdges[iSpec].data();

        for(Int_t iSample(0); iSample < fNumSamples; ++iSample)
        {
          if(iSample > 0 && !fSampling) { break; }
          // if(iSample > 0 && HasMass(PartSpecies(iSpec))) {  } // reconstructed are not sampled

          TH1* profile = nullptr;
          TH1* profileNeg = nullptr;
          TH1* profile3sub[3][3] = {nullptr};
          TProfile* profileEtaSlices[fFlowBinNumberEtaSlices] = {nullptr};

          switch(iSpec)
          {
            case kRefs :
            {
              profile = new TProfile(Form("%s_Pos_sample%d",corName,iSample), Form("%s: %s; %s",GetSpeciesLabel(PartSpecies(iSpec)),corLabel,GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax);
              if(fEtaCheckRFP) profileNeg = new TProfile(Form("%s_Neg_sample%d",corName,iSample), Form("%s: %s; %s",GetSpeciesLabel(PartSpecies(iSpec)),corLabel,GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax);
              if(bHas3sub){
                if(corrOrder == 2) {
                  for(Int_t rf1Pos(0); rf1Pos < 3; rf1Pos++){
                    if(rf1Pos > 1) break;
                    for(Int_t rf2Pos(0); rf2Pos < 3; rf2Pos++ ){
                      if(rf1Pos >= rf2Pos) continue;
                      profile3sub[rf1Pos][rf2Pos] = new TProfile(Form("%s_Pos_sample%d_rf1_%c_rf2_%c",corName,iSample,sides[rf1Pos],sides[rf2Pos]), Form("%s: %s; %s",GetSpeciesLabel(PartSpecies(iSpec)), corLabel,GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax);
                    }
                  }
                }
                else if(corrOrder == 4) {
                  for(Int_t twoPos(0); twoPos < 3; twoPos++){
                    profile3sub[0][twoPos] = new TProfile(Form("%s_Pos_sample%d_two_%c",corName,iSample,sides[twoPos]), Form("%s: %s; %s",GetSpeciesLabel(PartSpecies(iSpec)), corLabel,GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax);
                  }
                }
                else {
                  AliError(Form("AliUniFlowCorrTask %d : 3 subevents implemented only for 2- and 4-particle correlations.\n",iTask));
                  return;
                }
              }
              if(fCorrUsingGF){
                Double_t width = 2.0*fFlowEtaMax/fFlowBinNumberEtaSlices;
                for(Int_t iBin(0); iBin < fFlowBinNumberEtaSlices; iBin++){
                  Double_t lowerLimit = width*iBin;
                  Double_t upperLimit = width*(iBin+1);
                  profileEtaSlices[iBin] = new TProfile(Form("%s_eta_%.3g_%.3g_sample%d",corName,lowerLimit,upperLimit,iSample), Form("%s: %s; %s",GetSpeciesLabel(PartSpecies(iSpec)), corLabel,GetCentEstimatorLabel(fCentEstimator)),fCentBinNum,fCentMin,fCentMax);
                }
              }
              break;
            }

            case kCharged :
            case kPion :
            case kKaon :
            case kProton :
            case kCharUnidentified :
            {
                if(iNumPtFixBins > 0) {
                    profile = new TProfile2D(Form("%s_Pos_sample%d",corName,iSample), Form("%s: %s (Pos); %s; #it{p}_{T} (GeV/#it{c})",GetSpeciesLabel(PartSpecies(iSpec)), corLabel, GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax, iNumPtFixBins,dPtFixBinEdges,"");
                    if(bHasGap) { profileNeg = new TProfile2D(Form("%s_Neg_sample%d",corName,iSample), Form("%s: %s (Neg); %s; #it{p}_{T} (GeV/#it{c})",GetSpeciesLabel(PartSpecies(iSpec)), corLabel,GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax, iNumPtFixBins,dPtFixBinEdges,""); }
                } else {
                    profile = new TProfile2D(Form("%s_Pos_sample%d",corName,iSample), Form("%s: %s (Pos); %s; #it{p}_{T} (GeV/#it{c})",GetSpeciesLabel(PartSpecies(iSpec)), corLabel,GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax, fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax);
                    if(bHasGap) { profileNeg = new TProfile2D(Form("%s_Neg_sample%d",corName,iSample), Form("%s: %s (Neg); %s; #it{p}_{T} (GeV/#it{c})",GetSpeciesLabel(PartSpecies(iSpec)), corLabel,GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax, fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax); }
                }
                if(bHas3sub){
                  if(corrOrder == 2) {
                    for(Int_t poiPos(0); poiPos < 3; poiPos++)
                      for(Int_t rfPos(0); rfPos < 3; rfPos++){
                        if(poiPos == rfPos) continue;
                        profile3sub[poiPos][rfPos] = new TProfile2D(Form("%s_Pos_sample%d_poi_%c_rfp_%c",corName,iSample,sides[poiPos],sides[rfPos]), Form("%s: %s (Pos); %s; #it{p}_{T} (GeV/#it{c})",GetSpeciesLabel(PartSpecies(iSpec)), corLabel,GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax, fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax);
                      }
                  }
                  else if(corrOrder == 4) {
                    for(Int_t poiPos(0); poiPos < 3; poiPos++)
                      for(Int_t twoPos(0); twoPos < 3; twoPos++){
                        profile3sub[poiPos][twoPos] = new TProfile2D(Form("%s_Pos_sample%d_poi_%c_two_%c",corName,iSample,sides[poiPos],sides[twoPos]), Form("%s: %s (Pos); %s; #it{p}_{T} (GeV/#it{c})",GetSpeciesLabel(PartSpecies(iSpec)), corLabel,GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax, fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax);
                      }
                  }
                  else {
                    AliError(Form("AliUniFlowCorrTask %d : 3 subevents implemented only for 2- & 4-particle correlations.\n",iTask));
                    return;
                  }
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
                if(bHas3sub){
                  if(corrOrder == 2) {
                    for(Int_t poiPos(0); poiPos < 3; poiPos++)
                      for(Int_t rfPos(0); rfPos < 3; rfPos++){
                        if(poiPos == rfPos) continue;
                        profile3sub[poiPos][rfPos] = new TProfile3D(Form("%s_Pos_sample%d_poi_%c_rfp_%c",corName,iSample,sides[poiPos],sides[rfPos]), Form("%s: %s (Pos); %s; #it{p}_{T} (GeV/#it{c})",GetSpeciesLabel(PartSpecies(iSpec)), corLabel,GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax, fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax,fV0sNumBinsMass,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax);
                      }
                  }
                  else if(corrOrder == 4) {
                    for(Int_t poiPos(0); poiPos < 3; poiPos++)
                      for(Int_t twoPos(0); twoPos < 3; twoPos++){
                        profile3sub[poiPos][twoPos] = new TProfile3D(Form("%s_Pos_sample%d_poi_%c_two_%c",corName,iSample,sides[poiPos],sides[twoPos]), Form("%s: %s (Pos); %s; #it{p}_{T} (GeV/#it{c})",GetSpeciesLabel(PartSpecies(iSpec)), corLabel,GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax, fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax,fV0sNumBinsMass,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax);
                      }
                  }
                  else {
                    AliError(Form("AliUniFlowCorrTask %d : 3 subevents implemented only for 2- & 4-particle correlations.\n",iTask));
                    return;
                  }
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
                if(bHas3sub){
                  if(corrOrder == 2) {
                    for(Int_t poiPos(0); poiPos < 3; poiPos++)
                      for(Int_t rfPos(0); rfPos < 3; rfPos++){
                        if(poiPos == rfPos) continue;
                        profile3sub[poiPos][rfPos] = new TProfile3D(Form("%s_Pos_sample%d_poi_%c_rfp_%c",corName,iSample,sides[poiPos],sides[rfPos]), Form("%s: %s (Pos); %s; #it{p}_{T} (GeV/#it{c})",GetSpeciesLabel(PartSpecies(iSpec)), corLabel,GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax, fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax,fV0sNumBinsMass,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax);
                      }
                  }
                  else if(corrOrder == 4) {
                    for(Int_t poiPos(0); poiPos < 3; poiPos++)
                      for(Int_t twoPos(0); twoPos < 3; twoPos++){
                        profile3sub[poiPos][twoPos] = new TProfile3D(Form("%s_Pos_sample%d_poi_%c_two_%c",corName,iSample,sides[poiPos],sides[twoPos]), Form("%s: %s (Pos); %s; #it{p}_{T} (GeV/#it{c})",GetSpeciesLabel(PartSpecies(iSpec)), corLabel,GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax, fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax,fV0sNumBinsMass,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax);
                      }
                  }
                  else {
                    AliError(Form("AliUniFlowCorrTask %d : 3 subevents implemented only for 2- & 4-particle correlations.\n",iTask));
                    return;
                  }
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
                if(bHas3sub){
                  if(corrOrder == 2) {
                    for(Int_t poiPos(0); poiPos < 3; poiPos++)
                      for(Int_t rfPos(0); rfPos < 3; rfPos++){
                        if(poiPos == rfPos) continue;
                        profile3sub[poiPos][rfPos] = new TProfile3D(Form("%s_Pos_sample%d_poi_%c_rfp_%c",corName,iSample,sides[poiPos],sides[rfPos]), Form("%s: %s (Pos); %s; #it{p}_{T} (GeV/#it{c})",GetSpeciesLabel(PartSpecies(iSpec)), corLabel,GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax, fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax,fPhiNumBinsMass,fCutPhiInvMassMin,fCutPhiInvMassMax);
                      }
                  }
                  else if(corrOrder == 4) {
                    for(Int_t poiPos(0); poiPos < 3; poiPos++)
                      for(Int_t twoPos(0); twoPos < 3; twoPos++){
                        profile3sub[poiPos][twoPos] = new TProfile3D(Form("%s_Pos_sample%d_poi_%c_two_%c",corName,iSample,sides[poiPos],sides[twoPos]), Form("%s: %s (Pos); %s; #it{p}_{T} (GeV/#it{c})",GetSpeciesLabel(PartSpecies(iSpec)), corLabel,GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax, fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax,fPhiNumBinsMass,fCutPhiInvMassMin,fCutPhiInvMassMax);
                      }
                  }
                  else {
                    AliError(Form("AliUniFlowCorrTask %d : 3 subevents implemented only for 2- & 4-particle correlations.\n",iTask));
                    return;
                  }
                }
                break;
            }
          }

          if(!profile) { fInit = kFALSE; AliError("Profile (Pos) NOT created!"); task->PrintTask(); return; }

          // check if same profile does not exists already
          if(fListFlow[iSpec]->FindObject(profile->GetName())) {
            AliError(Form("AliUniFlowCorrTask %d : Profile '%s' already exists! Please check run macro for AliUniFlowCorrTask duplicates!",iTask,profile->GetName()));
            fInit = kFALSE;
            task->PrintTask();
            delete profile;
            return;
          }

          profile->Sumw2();
          fListFlow[iSpec]->Add(profile);

          if(bHasGap)
          { // Refs does not distinquish Pos/Neg
            if(iSpec != kRefs || fEtaCheckRFP){
              if(!profileNeg) { fInit = kFALSE; AliError("Profile (Neg) NOT created!"); task->PrintTask(); return; }
              // same for Neg
              if(fListFlow[iSpec]->FindObject(profileNeg->GetName())) {
                AliError(Form("AliUniFlowCorrTask %d : Profile '%s' already exists! Please check run macro for AliUniFlowCorrTask duplicates!",iTask,profile->GetName()));
                fInit = kFALSE;
                task->PrintTask();
                delete profileNeg;
                return;
              }
              profileNeg->Sumw2();
              fListFlow[iSpec]->Add(profileNeg);
            }

          if(bHas3sub){
            for(Int_t poiPos(0); poiPos < 3; poiPos++){
              if(iSpec == kRefs && ((corrOrder == 4 && poiPos > 0) || (corrOrder == 2 && poiPos > 1)) ) break;
              for(Int_t twoPos(0); twoPos < 3; twoPos++){
                if(iSpec == kRefs && corrOrder == 2 && poiPos >= twoPos ) continue;
                if(iSpec != kRefs && corrOrder == 2 && poiPos == twoPos) continue;
                if(!profile3sub[poiPos][twoPos]) { fInit = kFALSE; AliError("Profiles combi NOT created!"); task->PrintTask(); return; }
                if(fListFlow[iSpec]->FindObject(profile3sub[poiPos][twoPos]->GetName())) {
                  AliError(Form("AliUniFlowCorrTask %d : Profile '%s' already exists! Please check run macro for AliUniFlowCorrTask duplicates!",iTask,profile3sub[poiPos][twoPos]->GetName()));
                  fInit = kFALSE;
                  task->PrintTask();
                  delete profile3sub[poiPos][twoPos];
                  return;
                }
                profile3sub[poiPos][twoPos]->Sumw2();
                fListFlow[iSpec]->Add(profile3sub[poiPos][twoPos]);
                }
              }
            }
          } // end has gap

          if(fCorrUsingGF){
            for(Int_t iBin(0); iBin < fFlowBinNumberEtaSlices; iBin++){
              if(!profileEtaSlices[iBin]) { fInit = kFALSE; AliError("Profile (eta slices) NOT created!"); task->PrintTask(); return; }
              if(fListFlow[iSpec]->FindObject(profileEtaSlices[iBin]->GetName())) {
                AliError(Form("AliUniFlowCorrTask %d : Profile '%s' already exists! Please check run macro for AliUniFlowCorrTask duplicates!",iTask,profileEtaSlices[iBin]->GetName()));
                fInit = kFALSE;
                task->PrintTask();
                delete profileEtaSlices[iBin];
                return;
              }
              profileEtaSlices[iBin]->Sumw2();
              fListFlow[iSpec]->Add(profileEtaSlices[iBin]);
            }
          } // end di-hadron corr


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
    sLabelCand[SparseCand::kSample] = "iSample";
    TString sAxes = TString(); for(Int_t i(0); i < SparseCand::kDim; ++i) { sAxes += Form("%s; ",sLabelCand[i].Data()); }

    Int_t iNumBinsCand[SparseCand::kDim];
    Double_t dMinCand[SparseCand::kDim];
    Double_t dMaxCand[SparseCand::kDim];
    iNumBinsCand[SparseCand::kCent] = fCentBinNum;
    dMinCand[SparseCand::kCent] = fCentMin;
    dMaxCand[SparseCand::kCent] = fCentMax;
    iNumBinsCand[SparseCand::kEta] = fFlowEtaBinNum;
    dMinCand[SparseCand::kEta] = -fFlowEtaMax;
    dMaxCand[SparseCand::kEta] = fFlowEtaMax;
    iNumBinsCand[SparseCand::kSample] = fNumSamples;
    dMinCand[SparseCand::kSample] = 0;
    dMaxCand[SparseCand::kSample] = fNumSamples;

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

  //creating output 2D histograms for correlations & preparing event pool
  if(fCorrFill){
    std::vector<Double_t> centVec, vertexVec, ptVec, etaVec, phiVec;
    Double_t psiAr[2] = {-999.,999.};
    for(Int_t counter(0); counter < fCentBinNum+1; counter++){ centVec.push_back(fCentMin + counter*(fCentMax - fCentMin)/fCentBinNum); }
    for(Int_t counter(0); counter < 2*fPVtxCutZ+1; counter++){ vertexVec.push_back(-fPVtxCutZ + counter); }
    for(Int_t counter(0); counter < fCorrDEtaBinNum+1; counter++){ etaVec.push_back(fCorrdEtaMin + counter*(fCorrdEtaMax - fCorrdEtaMin)/fCorrDEtaBinNum); }
    for(Int_t counter(0); counter < fCorrDPhiBinNum+1; counter++){ phiVec.push_back(fCorrdPhiMin + counter*(fCorrdPhiMax - fCorrdPhiMin)/fCorrDPhiBinNum); }
    ptVec = {-999., 999.};
    if(fUsePtBinnedEventPool){
      ptVec = fFlowPOIsPtBinEdges[kCharged];
    }
    else { ptVec = {-999., 999.}; }
    Double_t* centAr = centVec.data();
    Double_t* vertAr = vertexVec.data();
    Double_t* ptAr = ptVec.data();
    Double_t* etaAr = etaVec.data();
    Double_t* phiAr = phiVec.data();
    Int_t sizePt = ptVec.size() - 1;

    fEventPoolMgr = new AliEventPoolManager(fPoolSize, fMixingTracks, fCentBinNum, centAr, 2*fPVtxCutZ, vertAr, 1, psiAr, sizePt, ptAr);
    fEventPoolMgr->SetTargetValues(fMixingTracks, 0.1, 5);
    if(!fEventPoolMgr){ AliError("AliEventPoolManager doesn't exist!"); }

    // fEventPoolMgr->Validate();

    const Int_t binsForCor[4] = {fCorrDEtaBinNum, fCorrDPhiBinNum, sizePt, fCentBinNum};
    for(Int_t iSpec(0); iSpec < kUnknown; iSpec++)
    {
      if(!fProcessSpec[iSpec]) continue;

      fh4CorrelationsSE[iSpec] = new THnSparseD(Form("fh4CorrelationsSE_%s",GetSpeciesName(PartSpecies(iSpec))), Form("%s: Distribution; #Delta #eta; #Delta #phi; p_{T} (assoc); centrality",GetSpeciesName(PartSpecies(iSpec))), 4, binsForCor);
      fh4CorrelationsSE[iSpec]->SetBinEdges(0,etaAr);
      fh4CorrelationsSE[iSpec]->SetBinEdges(1,phiAr);
      fh4CorrelationsSE[iSpec]->SetBinEdges(2,ptAr);
      fh4CorrelationsSE[iSpec]->SetBinEdges(3,centAr);
      fListFlow[iSpec]->Add(fh4CorrelationsSE[iSpec]);

      fh4CorrelationsME[iSpec] = new THnSparseD(Form("fh4CorrelationsME_%s",GetSpeciesName(PartSpecies(iSpec))), Form("%s: Distribution; #Delta #eta; #Delta #phi; p_{T} (assoc); centrality",GetSpeciesName(PartSpecies(iSpec))), 4, binsForCor);
      fh4CorrelationsME[iSpec]->SetBinEdges(0,etaAr);
      fh4CorrelationsME[iSpec]->SetBinEdges(1,phiAr);
      fh4CorrelationsME[iSpec]->SetBinEdges(2,ptAr);
      fh4CorrelationsME[iSpec]->SetBinEdges(3,centAr);
      fListFlow[iSpec]->Add(fh4CorrelationsME[iSpec]);
    }
  }

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

    if(fFlowFillWeightsMultiD){
      TString wLabelCand[SparseWeights::wDim];
      wLabelCand[SparseWeights::wPhi] = "#phi";
      wLabelCand[SparseWeights::wCent] = GetCentEstimatorLabel(fCentEstimator);
      wLabelCand[SparseWeights::wPt] = "#it{p}_{T} (GeV/c)";
      wLabelCand[SparseWeights::wEta] = "#eta";
      wLabelCand[SparseWeights::wVz] = "v_{z}";
      wLabelCand[SparseWeights::wSpec] = "Species";
      TString sAxesWeights = TString(); for(Int_t i(0); i < SparseWeights::wDim; ++i) { sAxesWeights += Form("%s; ",wLabelCand[i].Data()); }

      Int_t iNumBinsWeights[SparseWeights::wDim];
      Double_t dMinWeights[SparseWeights::wDim];
      Double_t dMaxWeights[SparseWeights::wDim];
      iNumBinsWeights[SparseWeights::wCent] = fCentBinNum;
      dMinWeights[SparseWeights::wCent] = fCentMin;
      dMaxWeights[SparseWeights::wCent] = fCentMax;
      iNumBinsWeights[SparseWeights::wPhi] = fFlowPhiBinNum;
      dMinWeights[SparseWeights::wPhi] = 0.0;
      dMaxWeights[SparseWeights::wPhi] = TMath::TwoPi();
      iNumBinsWeights[SparseWeights::wEta] = fFlowEtaBinNum;
      dMinWeights[SparseWeights::wEta] = -fFlowEtaMax;
      dMaxWeights[SparseWeights::wEta] = fFlowEtaMax;
      iNumBinsWeights[SparseWeights::wPt] = fFlowPOIsPtBinNum;
      dMinWeights[SparseWeights::wPt] = fFlowPOIsPtMin;
      dMaxWeights[SparseWeights::wPt] = fFlowPOIsPtMax;
      iNumBinsWeights[SparseWeights::wVz] = 2*fPVtxCutZ;
      dMinWeights[SparseWeights::wVz] = -fPVtxCutZ;
      dMaxWeights[SparseWeights::wVz] = fPVtxCutZ;
      iNumBinsWeights[SparseWeights::wSpec] = kUnknown;
      dMinWeights[SparseWeights::wSpec] = 0;
      dMaxWeights[SparseWeights::wSpec] = kUnknown;

      fhWeightsMultiD =
      new THnSparseD("fhWeightsMultiD",Form("Weights distribution; %s;", sAxesWeights.Data()), SparseWeights::wDim, iNumBinsWeights, dMinWeights, dMaxWeights);
      fhWeightsMultiD->Sumw2();
      fFlowWeights->Add(fhWeightsMultiD);
    }
  }

  // Selection / reconstruction counters : omni-present
  {
    // TString sEventCounterLabel[] = {"Input","Physics selection OK","EventCuts OK","Event OK","#RPFs OK","Multiplicity OK","Selected"};
    TString sEventCounterLabel[] = {"Input"};
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
    if( (fCentEstimator != kRFP) && (fColSystem == kPP || fColSystem == kPPb) ){
      Int_t counter = fCentBinNum/10;
      for(Int_t iCen(0); iCen < counter; ++iCen)
      {
        if(iCen > 9) { AliWarning("Incorrect number of centrality bins for pT vs. multiplicity histograms for small systems."); break; }
        fh2MeanMultRFP[iCen] = new TH2D(Form("fh2MeanMultCharged_Cent%d",iCen), "RFPs: pT vs. multiplicity; #it{p}_{T}; multiplicity", fFlowPOIsPtBinNum,fFlowPOIsPtMin,fFlowPOIsPtMax,200,0,200);
        fQACharged->Add(fh2MeanMultRFP[iCen]);
      }
    }
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

    fhEventSampling = new TH2D("fhEventSampling",Form("Event sampling; %s; sample index", GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax, fNumSamples,0,fNumSamples);
    fQAEvents->Add(fhEventSampling);

    if(fIsHMpp){
      fhV0Mamplitude = new TH1D("fhV0Mamplitude","; V0M amplitude; Counts",1000,0,1000);
      fQAEvents->Add(fhV0Mamplitude);

      fhV0MamplitudeRatio = new TH1D("fhV0MamplitudeRatio","; V0M / <V0M>; Counts",150,0,15);
      fQAEvents->Add(fhV0MamplitudeRatio);

      fh2V0MnCharged = new TH2D("fh2V0MnCharged", "; V0M / <V0M>; N_{ch}", 150,0,15,200,0,200);
      fQAEvents->Add(fh2V0MnCharged);
    }

    if(fAnalType == kMC){
      fh2MCip = new TH2D("fh2MCip", "RFPs: impact parameter vs. multiplicity; b; multiplicity", 200,0,20,fCentBinNum,fCentMin,fCentMax);
      fQAEvents->Add(fh2MCip);
    }
    else{
      fhEventCentrality = new TH1D("fhEventCentrality",Form("Event centrality (%s); %s", GetCentEstimatorLabel(fCentEstimator), GetCentEstimatorLabel(fCentEstimator)), fCentBinNum,fCentMin,fCentMax);
      fQAEvents->Add(fhEventCentrality);

      // event histogram
      fEventCuts.AddQAplotsToList(fQAEvents);

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
// ============================================================================
Double_t AliAnalysisTaskUniFlow::PIDCorrection(const AliAODTrack *track, const PartSpecies species) const
{
  // PID correction for 2018 data
  // from PWG-HF

  if(!track) { AliError("Track not exists!"); return -999.9; }
  if(species == kUnknown) { AliError("Invalid species 'Unknown'!"); return -999.9; }

  Int_t iRunNumber = fEventAOD->GetRunNumber();
  Int_t nPbins = 8;
  Float_t pTPClims[] = {0.3,0.5,0.75,1.,1.5,2.,3.,5.,10.};
  Int_t nEtabins = 5;
  Float_t absetalims[] = {0.,0.1,0.2,0.4,0.6,0.8};

  Double_t NumPion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
  Double_t NumKaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
  Double_t NumProton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);

  Double_t eta = TMath::Abs( track->Eta() );
  Double_t pTPC = track->GetTPCmomentum();
  Double_t SigmaValue = -999.9;
  if (species == kPion){
     SigmaValue = NumPion;
  }
  else if (species == kKaon){
     SigmaValue = NumKaon;
  }
  else if (species == kProton){
      SigmaValue = NumProton;
  }

  std::vector<Double_t> meanPion, meanKaon, meanProton, sigmaPion, sigmaKaon, sigmaProton;

  if(pTPC < 0.3) return SigmaValue;
  if(fCentEstimator != kV0M) return SigmaValue;

  if(fIndexCentrality > 50) return SigmaValue;

  if(iRunNumber>=295585 && iRunNumber<=296623){
    //LHC18q
    if(fIndexCentrality >= 0 && fIndexCentrality<= 10){
      meanPion = {
        -0.656082, -0.604754, -0.63195, -0.669819, -0.708323, -0.746162, -0.800557, -0.893548,//eta 0
        -0.711512, -0.686848, -0.711177, -0.752377, -0.781678, -0.838341, -0.80682, -0.917012,//eta 1
        -0.650884, -0.706274, -0.752449, -0.798673, -0.835543, -0.859084, -0.823375, -0.854207,//eta 2
        -0.479705, -0.700093, -0.815494, -0.895986, -0.958185, -0.980633, -0.967819, -0.996364,//eta 3
        -0.264033, -0.637537, -0.828537, -0.952223, -1.06091, -1.10456, -1.09789, -1.15579 }; //eta 4
      meanKaon = {
        -0.48114, -0.672897, -0.625657, -0.776678, -0.786824, -0.708909, -0.822472, -0.491422,//eta 0
        -0.432004, -0.71966, -0.708949, -0.870034, -0.856239, -0.825942, -0.871391, -1.17962,//eta 1
        -0.336167, -0.71892, -0.735327, -0.93073, -0.864011, -0.891611, -0.924604, -0.735026,//eta 2
        -0.391054, -0.716122, -0.796868, -1.0208, -0.984637, -0.998813, -1.01377, -1.06832,//eta 3
        -0.551251, -0.696801, -0.815058, -1.05691, -1.06688, -1.06648, -1.07023, -1.05183};//eta 4
      meanProton = {
        -0.200581, -0.16751, -0.451043, -0.952266, -0.852847, -0.760682, -0.676723, -0.603716,//eta 0
        -0.123522, -0.128086, -0.444873, -0.846087, -0.988114, -1.05446, -0.761678, -0.785548,//eta 1
        -0.100534, -0.1431, -0.448783, -0.8385, -0.843197, -1.05925, -0.878891, -0.69573,//eta 2
        -0.233023, -0.317509, -0.598837, -0.945108, -1.01043, -1.21354, -0.99634, -0.915479,//eta 3
        -0.391233, -0.516667, -0.768414, -1.00696, -1.03589, -1.26272, -1.02806, -0.994112};//eta 4
      sigmaPion = {
        0.986692, 1.01991, 1.00333, 0.986744, 0.981785, 1.04139, 1.0638, 1.09162,//eta 0
        0.968236, 0.999018, 0.984673, 0.963658, 0.963348, 0.991749, 1.00931, 1.07714,//eta 1
        0.948544, 0.971808, 0.957514, 0.938766, 0.936994, 0.987149, 0.957657, 0.994133,//eta 2
        0.92104, 0.931385, 0.916291, 0.896921, 0.890271, 0.926601, 0.891902, 0.905638,//eta 3
        0.909424, 0.89589, 0.881729, 0.860767, 0.842961, 0.873783, 0.83704, 0.758586};//eta 4
      sigmaKaon = {
        0.891168, 1.00326, 1.04053, 0.952367, 0.919632, 0.951424, 0.902434, 1.07759,//eta 0
        0.853805, 0.968212, 1.01839, 0.930692, 0.918841, 0.92428, 0.913715, 0.929855,//eta 1
        0.830364, 0.911894, 0.987625, 0.912264, 0.939659, 0.906548, 0.893721, 1.00634,//eta 2
        0.718803, 0.882484, 0.959253, 0.870302, 0.907273, 0.881795, 0.867757, 0.906373,//eta 3
        0.688955, 0.855596, 0.932222, 0.839553, 0.867639, 0.855212, 0.845177, 0.90448};//eta 4
      sigmaProton = {
        0.771648, 0.841043, 0.917283, 1.12449, 1.0023, 0.952976, 0.963016, 1.01111,//eta 0
        0.752951, 0.825488, 0.883897, 1.02998, 1.07061, 0.866346, 0.952794, 0.916068,//eta 1
        0.72623, 0.799905, 0.860896, 0.996524, 1.02278, 0.866737, 0.834891, 0.988606,//eta 2
        0.70374, 0.759504, 0.811721, 0.944681, 0.96305, 0.81128, 0.85648, 0.900749,//eta 3
        0.702538, 0.723393, 0.781419, 0.83867, 0.940137, 0.785817, 0.841202, 0.859564 };//eta 4
    } //end 0-10% centrality
    else if(fIndexCentrality >= 10 && fIndexCentrality <= 30){
      meanPion = {
        -0.656082, -0.604754, -0.63195, -0.669819, -0.708323, -0.746162, -0.800557, -0.893548,//eta 0
        -0.711512, -0.686848, -0.711177, -0.752377, -0.781678, -0.838341, -0.80682, -0.917012,//eta 1
        -0.650884, -0.706274, -0.752449, -0.798673, -0.835543, -0.859084, -0.823375, -0.854207,//eta 2
        -0.479705, -0.700093, -0.815494, -0.895986, -0.958185, -0.980633, -0.967819, -0.996364,//eta 3
        -0.264033, -0.637537, -0.828537, -0.952223, -1.06091, -1.10456, -1.09789, -1.15579};//eta 4
      meanKaon = {
        -0.48114, -0.672897, -0.625657, -0.776678, -0.786824, -0.708909, -0.822472, -0.491422,//eta 0
        -0.432004, -0.71966, -0.708949, -0.870034, -0.856239, -0.825942, -0.871391, -1.17962,//eta 1
        -0.336167, -0.71892, -0.735327, -0.93073, -0.864011, -0.891611, -0.924604, -0.735026,//eta 2
        -0.391054, -0.716122, -0.796868, -1.0208, -0.984637, -0.998813, -1.01377, -1.06832,//eta 3
        -0.551251, -0.696801, -0.815058, -1.05691, -1.06688, -1.06648, -1.07023, -1.05183};//eta 4
      meanProton = {
        -0.200581, -0.16751, -0.451043, -0.952266, -0.852847, -0.760682, -0.676723, -0.603716,//eta 0
        -0.123522, -0.128086, -0.444873, -0.846087, -0.988114, -1.05446, -0.761678, -0.785548,//eta 1
        -0.100534, -0.1431, -0.448783, -0.8385, -0.843197, -1.05925, -0.878891, -0.69573,//eta 2
        -0.233023, -0.317509, -0.598837, -0.945108, -1.01043, -1.21354, -0.99634, -0.915479,//eta 3
        -0.391233, -0.516667, -0.768414, -1.00696, -1.03589, -1.26272, -1.02806, -0.994112};//eta 4
      sigmaPion = {
        0.986692, 1.01991, 1.00333, 0.986744, 0.981785, 1.04139, 1.0638, 1.09162,//eta 0
        0.968236, 0.999018, 0.984673, 0.963658, 0.963348, 0.991749, 1.00931, 1.07714,//eta 1
        0.948544, 0.971808, 0.957514, 0.938766, 0.936994, 0.987149, 0.957657, 0.994133,//eta 2
        0.92104, 0.931385, 0.916291, 0.896921, 0.890271, 0.926601, 0.891902, 0.905638,//eta 3
        0.909424, 0.89589, 0.881729, 0.860767, 0.842961, 0.873783, 0.83704, 0.758586};//eta 4
      sigmaKaon = {
        0.891168, 1.00326, 1.04053, 0.952367, 0.919632, 0.951424, 0.902434, 1.07759,//eta 0
        0.853805, 0.968212, 1.01839, 0.930692, 0.918841, 0.92428, 0.913715, 0.929855,//eta 1
        0.830364, 0.911894, 0.987625, 0.912264, 0.939659, 0.906548, 0.893721, 1.00634,//eta 2
        0.718803, 0.882484, 0.959253, 0.870302, 0.907273, 0.881795, 0.867757, 0.906373,//eta 3
        0.688955, 0.855596, 0.932222, 0.839553, 0.867639, 0.855212, 0.845177, 0.90448};//eta 4
      sigmaProton = {
        0.771648, 0.841043, 0.917283, 1.12449, 1.0023, 0.952976, 0.963016, 1.01111,//eta 0
        0.752951, 0.825488, 0.883897, 1.02998, 1.07061, 0.866346, 0.952794, 0.916068,//eta 1
        0.72623, 0.799905, 0.860896, 0.996524, 1.02278, 0.866737, 0.834891, 0.988606,//eta 2
        0.70374, 0.759504, 0.811721, 0.944681, 0.96305, 0.81128, 0.85648, 0.900749,//eta 3
        0.702538, 0.723393, 0.781419, 0.83867, 0.940137, 0.785817, 0.841202, 0.859564 };//eta 4
    } //end 10-30% centrality
    else if(fIndexCentrality >= 30 && fIndexCentrality <= 50){
      meanPion = {
        -0.537046, -0.427744, -0.411915, -0.42242, -0.445157, -0.423209, -0.403354, -0.39011,
        -0.451747, -0.358803, -0.347827, -0.360332, -0.379419, -0.373067, -0.309076, -0.201842,
        -0.37701, -0.342516, -0.352131, -0.375476, -0.397523, -0.380363, -0.348125, -0.334354,
        -0.340843, -0.438716, -0.504516, -0.554322, -0.614602, -0.624125, -0.612949, -0.616095,
        -0.273705, -0.510522, -0.643307, -0.739041, -0.830935, -0.860098, -0.885069, -0.956967};
      meanKaon = {
        -0.226216, -0.427422, -0.414774, -0.562412, -0.521969, -0.467143, -0.414713, -0.330372,
        -0.16309, -0.36387, -0.436445, -0.395585, -0.452207, -0.408419, -0.30696, -0.323571,
        -0.106973, -0.382832, -0.361131, -0.402223, -0.452784, -0.400874, -0.342613, -0.185365,
        -0.252347, -0.480454, -0.500724, -0.573321, -0.642967, -0.616602, -0.546648, -0.482116,
        -0.5916, -0.593699, -0.6563, -0.683526, -0.807421, -0.815922, -0.785543, -0.842433};
      meanProton = {
         0.00677222, 0.0347718, -0.211127, -0.466866, -0.323172, -0.53392, -0.504211, -0.334974,
         0.0935506, 0.0970568, -0.138627, -0.392521, -0.399267, -0.43474, -0.200821, -0.23501,
         0.075394, 0.0609517, -0.170246, -0.409987, -0.420188, -0.448851, -0.267424, -0.313302,
         -0.133011, -0.210294, -0.42431, -0.562203, -0.459603, -0.673718, -0.649959, -0.520375,
         -0.324865, -0.495658, -0.69697, -0.814164, -0.710279, -0.778491, -0.80033, -0.76221};
      sigmaPion = {
        0.915632, 0.93365, 0.932587, 0.931425, 0.922551, 0.92571, 0.881836, 0.796746,
        0.906096, 0.925435, 0.924713, 0.919556, 0.906064, 0.911215, 0.866192, 0.773724,
        0.881485, 0.902513, 0.901312, 0.893701, 0.89239, 0.875522, 0.825261, 0.764695,
        0.835562, 0.850935, 0.853431, 0.846889, 0.834667, 0.819543, 0.798447, 0.774576,
        0.809872, 0.807042, 0.80727, 0.799673, 0.785673, 0.757604, 0.74583, 0.726052};
      sigmaKaon = {
        0.814089, 0.877054, 0.944152, 0.90886, 0.92537, 0.931201, 0.926482, 0.722974,
        0.798753, 0.841944, 0.952407, 0.952185, 0.929863, 0.921988, 0.932327, 0.828574,
        0.733393, 0.821445, 0.903447, 0.928461, 0.917579, 0.921567, 0.923788, 0.727546,
        0.694873, 0.772357, 0.864553, 0.886463, 0.866449, 0.86353, 0.872874, 0.990303,
        0.687564, 0.747428, 0.792902, 0.87729, 0.824179, 0.814195, 0.813793, 0.901867};
      sigmaProton = {
         0.758072, 0.796573, 0.838565, 0.942299, 0.990804, 0.914681, 0.900297, 0.872938,
         0.733353, 0.776217, 0.823004, 0.922272, 1.00522, 0.91665, 0.986521, 1.05648,
         0.715624, 0.75819, 0.79931, 0.897033, 0.972885, 0.920213, 0.942288, 0.978577,
         0.691934, 0.724153, 0.7582, 0.793117, 0.879233, 0.832657, 0.8289, 0.932063,
         0.695097, 0.694242, 0.719957, 0.790272, 0.855222, 0.820038, 0.789541, 0.834146};
    } //end 30-50% centrality
  } //end LHC18q
  else if(iRunNumber>=296690 && iRunNumber<=297595){
    //LHC18r
    if(fIndexCentrality >= 0 && fIndexCentrality<= 10){
      meanPion = {
        -0.744242, -0.715713, -0.69801, -0.696472, -0.717912, -0.767909, -0.822175, -0.883157,
        -0.976407, -0.964248, -0.976588, -0.969265, -1.00251, -1.04185, -1.08507, -1.02488,
        -0.938573, -1.02253, -1.0532, -1.06874, -1.09608, -1.11066, -1.07855, -1.06274,
        -0.462091, -0.766549, -0.875959, -0.918783, -0.979887, -0.984493, -0.945828, -0.954307,
         0.154123, -0.361271, -0.568491, -0.667592, -0.782836, -0.751772, -0.732903, -0.749};
      meanKaon = {
        -0.468947, -0.636701, -0.601858, -0.806051, -0.94714, -0.842379, -0.955165, -0.898824,
        -0.588647, -0.883708, -0.894757, -1.09769, -1.11786, -1.08056, -1.15336, -1.44054,
        -0.462369, -0.900829, -0.959231, -1.19209, -1.17182, -1.21788, -1.26831, -1.46315,
        -0.28288, -0.585668, -0.942842, -1.13897, -1.18188, -1.1556, -1.20724, -1.06756,
        -0.0830475, -0.129884, -0.40388, -0.905485, -1.03586, -0.963208, -0.95807, -0.591766};
      meanProton = {
        -0.43448, -0.41261, -0.468653, -0.766399, -0.906529, -0.87423, -0.925983, -0.834281,
        -0.439377, -0.510631, -0.648086, -0.99403, -1.09146, -1.14373, -1.19123, -0.993241,
        -0.400341, -0.54514, -0.649413, -0.979681, -1.22253, -1.27323, -1.27736, -1.12044,
        -0.295184, -0.441872, -0.470702, -0.910766, -1.04581, -1.17824, -1.17277, -0.978326,
        -0.169806, -0.33929, -0.309714, -0.680191, -0.862101, -0.972894, -0.951602, -0.676351 };
      sigmaPion = {
        1.19971, 1.20244, 1.19018, 1.18674, 1.19735, 1.27442, 1.31886, 1.35234,
        1.17433, 1.18271, 1.16982, 1.17218, 1.17712, 1.26327, 1.33523, 1.30084,
        1.16732, 1.17134, 1.16173, 1.15583, 1.16336, 1.24219, 1.23569, 1.24103,
        1.1773, 1.16117, 1.14767, 1.14647, 1.19338, 1.22358, 1.21504, 1.20365,
        1.21022, 1.17586, 1.16182, 1.15354, 1.15374, 1.22138, 1.19871, 1.20838};
      sigmaKaon = {
        1.2537, 1.29365, 1.27439, 1.1386, 1.06835, 1.10702, 1.03569, 1.12542,
        1.20352, 1.25335, 1.24745, 1.12229, 1.12846, 1.11836, 1.0746, 1.33407,
        1.17982, 1.21396, 1.23344, 1.13316, 1.15827, 1.10607, 1.06816, 1.14628,
        1.10113, 1.23381, 1.30511, 1.11268, 1.11029, 1.1049, 1.02285, 1.26782,
        1.09653, 1.23731, 1.28769, 1.10684, 1.09593, 1.11015, 1.04836, 1.12603};
      sigmaProton = {
        1.13405, 1.18163, 1.24085, 1.27282, 1.12543, 1.0912, 1.04366, 1.10697,
        1.09801, 1.16854, 1.22617, 1.26278, 1.25383, 1.07191, 1.04581, 1.11811,
        1.11623, 1.17665, 1.22384, 1.24119, 1.23884, 1.04855, 1.05348, 1.11713,
        1.08417, 1.16621, 1.22578, 1.34172, 1.24733, 1.04892, 1.05745, 1.14111,
        1.07421, 1.15563, 1.23173, 1.33796, 1.18478, 1.07793, 1.05178, 1.18215};
    } //end 0-10% centrality
    else if(fIndexCentrality >= 10 && fIndexCentrality <= 30){
      meanPion = {
        -0.694117,-0.631546,-0.587532,-0.568296,-0.572205,-0.587273,-0.592046,-0.597132,
        -0.845603,-0.794267,-0.782727,-0.761903,-0.770849,-0.773905,-0.763973,-0.764653,
        -0.790168,-0.822974,-0.830126,-0.824832,-0.834853,-0.822694,-0.799294,-0.760148,
        -0.377544,-0.606416,-0.680226,-0.709356,-0.748197,-0.735027,-0.687131,-0.67114,
        0.152892,-0.270204,-0.435937,-0.519352,-0.600842,-0.58243,-0.556651,-0.595698};
      meanKaon = {
        -0.32462,-0.482121,-0.52681,-0.702296,-0.801404,-0.736053,-0.758662,-0.649379,
        -0.418126,-0.666333,-0.708767,-0.90688,-0.927897,-0.882262,-0.865255,-1.01626,
        -0.326047,-0.690534,-0.767469,-0.967351,-0.968556,-0.963625,-0.937979,-0.87335,
        -0.140399,-0.422132,-0.693828,-0.903177,-0.947538,-0.905028,-0.873963,-0.703542,
        -0.0417266,-0.0456749,-0.283011,-0.714638,-0.824359,-0.758533,-0.708575,-0.496457};
      meanProton = {
        -0.292022,-0.284451,-0.320832,-0.571975,-0.719106,-0.745712,-0.750995,-0.728083,
        -0.294409,-0.357698,-0.451651,-0.716477,-0.883707,-0.890497,-0.915806,-0.789553,
        -0.274963,-0.404285,-0.462159,-0.704062,-0.864195,-1.00387,-0.968421,-0.796727,
        -0.213432,-0.360112,-0.346473,-0.658365,-0.824904,-0.895539,-0.902184,-0.814731,
        -0.139315,-0.309929,-0.251431,-0.465551,-0.697841,-0.764322,-0.748536,-0.549089};
      sigmaPion = {
        1.14223,1.15078,1.14811,1.1507,1.16023,1.20415,1.21901,1.16923,
        1.12642,1.13908,1.13519,1.13835,1.14473,1.19325,1.2187,1.22128,
        1.12001,1.13085,1.12843,1.12838,1.13356,1.17258,1.18389,1.15264,
        1.12539,1.12088,1.11835,1.11784,1.14067,1.15473,1.13739,1.10985,
        1.14555,1.12278,1.12099,1.11785,1.11802,1.15036,1.12267,1.11536};
      sigmaKaon = {
        1.17726,1.21852,1.22302,1.11337,1.07194,1.08253,1.01122,1.10989,
        1.13178,1.1722,1.1896,1.10195,1.10604,1.10493,1.0549,1.22236,
        1.11222,1.15166,1.19099,1.10719,1.12581,1.10115,1.0574,1.0145,
        1.06012,1.16507,1.23038,1.09635,1.10469,1.09994,1.03767,1.1339,
        1.08241,1.16377,1.21056,1.09332,1.09307,1.10426,1.03828,1.10098};
      sigmaProton = {
        1.09967,1.14495,1.19153,1.2288,1.11384,1.08263,1.06448,1.12628,
        1.07441,1.13454,1.17493,1.20853,1.151,1.08579,1.07429,1.16777,
        1.09054,1.145,1.17888,1.18472,1.18271,1.06515,1.07872,1.17542,
        1.06934,1.13974,1.18288,1.25817,1.22109,1.07178,1.06567,1.11053,
        1.06432,1.12584,1.19036,1.2397,1.11441,1.08026,1.053,1.15935};
    } //end 10-30% centrality
    else if(fIndexCentrality >= 30 && fIndexCentrality <= 50){
      meanPion = {
        -0.643992, -0.547379, -0.477054, -0.440121, -0.426498, -0.406638, -0.361916, -0.311107,
        -0.714799, -0.624286, -0.588865, -0.554541, -0.539189, -0.50596, -0.442875, -0.504426,
        -0.641763, -0.623418, -0.607051, -0.580924, -0.573627, -0.534728, -0.520039, -0.457556,
        -0.292997, -0.446284, -0.484493, -0.499929, -0.516507, -0.485561, -0.428434, -0.387974,
        0.151661, -0.179137, -0.303383, -0.371113, -0.418849, -0.413089, -0.380399, -0.442395};
      meanKaon = {
        -0.180292, -0.327542, -0.451762, -0.598542, -0.655667, -0.629728, -0.562158, -0.399934,
        -0.247604, -0.448958, -0.522777, -0.71607, -0.737935, -0.683963, -0.577149, -0.59199,
        -0.189726, -0.48024, -0.575707, -0.742612, -0.765291, -0.709371, -0.607647, -0.28355,
        0.00208301, -0.258596, -0.444814, -0.667385, -0.713196, -0.654455, -0.540685, -0.339524,
        -0.000405731, 0.0385342, -0.162143, -0.523791, -0.612858, -0.553858, -0.45908, -0.401148};
      meanProton = {
        -0.149563, -0.156293, -0.17301, -0.377551, -0.531683, -0.617193, -0.576007, -0.621884,
        -0.14944, -0.204765, -0.255216, -0.438923, -0.675954, -0.637264, -0.640382, -0.585865,
        -0.149585, -0.263429, -0.274904, -0.428442, -0.50586, -0.734513, -0.659482, -0.473013,
        -0.13168, -0.278353, -0.222243, -0.405964, -0.603998, -0.612839, -0.631597, -0.651136,
        -0.108823, -0.280568, -0.193149, -0.250911, -0.53358, -0.55575, -0.545469, -0.421828};
      sigmaPion = {
        1.08474, 1.09912, 1.10605, 1.11466, 1.12312, 1.13388, 1.11915, 0.98612,
        1.0785, 1.09545, 1.10056, 1.10451, 1.11234, 1.12323, 1.10217, 1.14173,
        1.0727, 1.09036, 1.09512, 1.10093, 1.10376, 1.10297, 1.1321, 1.06425,
        1.07348, 1.08058, 1.08902, 1.08921, 1.08795, 1.08588, 1.05973, 1.01606,
        1.08087, 1.06969, 1.08016, 1.08215, 1.08229, 1.07935, 1.04663, 1.02234};
      sigmaKaon = {
        1.10082, 1.14338, 1.17166, 1.08815, 1.07553, 1.05804, 0.98676, 1.09437,
        1.06003, 1.09105, 1.13176, 1.08161, 1.08362, 1.0915, 1.0352, 1.11064,
        1.04461, 1.08936, 1.14854, 1.08121, 1.09336, 1.09623, 1.04665, 0.882711,
        1.01912, 1.09634, 1.15565, 1.08002, 1.09909, 1.09498, 1.05249, 0.999984,
        1.0683, 1.09023, 1.13342, 1.07979, 1.0902, 1.09837, 1.0282, 1.07592};
       sigmaProton = {
         1.06529, 1.10828, 1.14221, 1.18478, 1.10225, 1.07407, 1.08529, 1.14558,
         1.05081, 1.10053, 1.1237, 1.15429, 1.04817, 1.09967, 1.10276, 1.21743,
         1.06485, 1.11335, 1.13392, 1.12825, 1.12658, 1.08174, 1.10395, 1.2337,
         1.05451, 1.11328, 1.13998, 1.17462, 1.19486, 1.09464, 1.07388, 1.07995,
         1.05442, 1.09606, 1.14899, 1.14143, 1.04405, 1.08258, 1.05422, 1.13655};
    } //end 30-50% centrality
  }
  else {
    return SigmaValue;
  }

  for (Int_t iEta(0); iEta < nEtabins; iEta++){
     for(Int_t ipTPC(0); ipTPC < nPbins; ipTPC++){
        if(eta > absetalims[iEta+1] || eta < absetalims[iEta]) continue;
        if(pTPC > pTPClims[ipTPC+1] || pTPC < pTPClims[ipTPC]) continue;
        if(species == kPion){
          SigmaValue = (NumPion-meanPion[iEta*8+ipTPC])/sigmaPion[iEta*8+ipTPC];
        }
        else if(species == kKaon){
          SigmaValue = (NumKaon-meanKaon[iEta*8+ipTPC])/sigmaKaon[iEta*8+ipTPC];
        }
        else if(species == kProton){
          SigmaValue = (NumProton-meanProton[iEta*8+ipTPC])/sigmaProton[iEta*8+ipTPC];
        }
      } //end ipTPC
    } //end iEta

    return SigmaValue;
}
// ============================================================================
const char* AliAnalysisTaskUniFlow::ReturnPPperiod(const Int_t runNumber) const
{
  Bool_t isHM = kFALSE;
  if(fTrigger == AliVEvent::kHighMultV0) isHM = kTRUE;

  if(runNumber >= 252235 && runNumber <= 264347){ // LHC16
    if(!isHM && runNumber >= 252235 && runNumber <= 252375) return "LHC16de"; //d
    if(!isHM && runNumber >= 253437 && runNumber <= 253591) return "LHC16de"; //e
    if(runNumber >= 254128 && runNumber <= 254332) return "LHC16ghi"; //g
    if(runNumber >= 254604 && runNumber <= 255467) return "LHC16ghi"; //h
    if(runNumber >= 255539 && runNumber <= 255618) return "LHC16ghi"; //i
    if(runNumber >= 256219 && runNumber <= 256418) return "LHC16j";
    if(runNumber >= 256941 && runNumber <= 258537) return "LHC16k";
    if(runNumber >= 258962 && runNumber <= 259888) return "LHC16l";
    if(runNumber >= 262424 && runNumber <= 264035) return "LHC16o";
    if(runNumber >= 264076 && runNumber <= 264347) return "LHC16p";
  }

  if(runNumber >= 270581 && runNumber <= 282704){ // LHC17
    if(!isHM && runNumber >= 270581 && runNumber <= 270667) return "LHC17ce";
    if(runNumber >= 270822 && runNumber <= 270830){
      if(isHM) return "averaged";
      else return "LHC17ce";
    }
    if(runNumber >= 270854 && runNumber <= 270865){
      if(isHM) return "averaged";
      else return "LHC17f";
    }
    if(runNumber >= 271870 && runNumber <= 273103) return "LHC17h";
    if(runNumber >= 273591 && runNumber <= 274442) return "LHC17i";
    if(!isHM && runNumber >= 274593 && runNumber <= 274671) return "LHC17j";
    if(runNumber >= 274690 && runNumber <= 276508) return "LHC17k";
    if(runNumber >= 276551 && runNumber <= 278216) return "LHC17l";
    if(runNumber >= 278914 && runNumber <= 280140) return "LHC17m";
    if(runNumber >= 280282 && runNumber <= 281961) return "LHC17o";
    if(runNumber >= 282528 && runNumber <= 282704) return "LHC17r";
  }

  if(runNumber >= 285009 && runNumber <= 294925){ // LHC18
    if(runNumber >= 285009 && runNumber <= 285396){
      if(isHM) return "LHC18bd";
      else return "LHC18b";
    }
    if(runNumber >= 285978 && runNumber <= 286350){
      if(isHM) return "LHC18bd";
      else return "LHC18d";
    }
    if(runNumber >= 286380 && runNumber <= 286937) return "LHC18e";
    if(runNumber >= 287000 && runNumber <= 287658) return "LHC18f";
    if(runNumber >= 288804 && runNumber <= 288806){
      if(isHM) return "LHC18hjk";
      else return "LHC18ghijk";
    }
    if(runNumber == 288943){
      if(isHM) return "LHC18hjk";
      else return "LHC18ghijk";
    }
    if(runNumber >= 289165 && runNumber <= 289201){
      if(isHM) return "LHC18hjk";
      else return "LHC18ghijk";
    }
    if(!isHM && runNumber >= 288619 && runNumber <= 288750) return "LHC18ghijk"; //g, no HM event, only MB
    if(!isHM && runNumber >= 288861 && runNumber <= 288909) return "LHC18ghijk"; //i, no HM event, only MB
    if(runNumber >= 289240 && runNumber <= 289971) return "LHC18l";
    if(runNumber >= 290323 && runNumber <= 292839){
      if(isHM) return "LHC18m";
      else return "LHC18mn";
    }
    if(!isHM && runNumber >= 293357 && runNumber <= 293359) return "LHC18mn"; //n, no HM event, only MB
    if(runNumber >= 293475 && runNumber <= 293898) return "LHC18o";
    if(runNumber >= 294009 && runNumber <= 294925) return "LHC18p";
  }

  AliWarning("Unknown period! Returning averaged weights");
  return "averaged";
}

#endif
