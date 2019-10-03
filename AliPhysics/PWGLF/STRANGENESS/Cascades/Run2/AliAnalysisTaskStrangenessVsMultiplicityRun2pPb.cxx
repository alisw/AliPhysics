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

// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// This task is meant to provide a lightweight alternative to the regular
// V0 and cascade analysis tasks that were used in Run 1, which produced full-size
// TTree objects for V0 and cascade candidates. Instead of that, the output
// for this task has been structured as follows:
//
//  Output 1: TList object containing some standard analysis histograms
//            for event counting.
//
//  Output 2: TList object containing all registered output objects for
//            the V0 analysis in AliV0Result format. The output objects will
//            each hold a TH3F with analysis-relevant information and the
//            configuration that was used to get to that.
//
//  Output 3: TList object containing all registered output objects for
//            the cascade analysis in AliCascadeResult format. The output objects will
//            each hold a TH3F with analysis-relevant information and the
//            configuration that was used to get to that.
//
//  Output 4: (optional) TTree object holding event characteristics. At the
//            moment only related to a single centrality estimator (default:
//            V0M). No downscaling option exists.
//
//  Output 5: (optional) TTree object containing V0 candidates to allow for
//            Run 1 legacy code to function and to allow for full information
//            to be saved. To keep output under control, a downscaling factor
//            can be applied (default: 0.001, configurable)
//
//  Output 6: (optional) TTree object containing cascade candidates to allow for
//            Run 1 legacy code to function and to allow for full information
//            to be saved. To keep output under control, a downscaling factor
//            can be applied (default: 0.001, configurable)
//
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class TTree;
class TParticle;
class TVector3;

//class AliMCEventHandler;
//class AliMCEvent;
//class AliStack;

class AliESDVertex;
class AliAODVertex;
class AliESDv0;
class AliAODv0;

#include <numeric>

#include <Riostream.h>
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "THnSparse.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "TProfile.h"
//#include "AliLog.h"

#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"
#include "AliLightV0vertexer.h"
#include "AliLightCascadeVertexer.h"
#include "AliESDpid.h"
#include "AliExternalTrackParam.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliCentrality.h"
#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"

#include "AliAODForwardMult.h"
#include "AliForwardUtil.h"
#include "AliTriggerIR.h"
#include "AliCFContainer.h"
#include "AliMultiplicity.h"
#include "AliAODMCParticle.h"
#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliESDUtils.h"
#include "AliGenEventHeader.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliEventCuts.h"
#include "AliV0Result.h"
#include "AliCascadeResult.h"
#include "AliAnalysisTaskStrangenessVsMultiplicityRun2pPb.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskStrangenessVsMultiplicityRun2pPb)

AliAnalysisTaskStrangenessVsMultiplicityRun2pPb::AliAnalysisTaskStrangenessVsMultiplicityRun2pPb()
: AliAnalysisTaskSE(), fListHist(0), fListV0(0), fListCascade(0), fTreeEvent(0), fTreeV0(0), fTreeCascade(0), fPIDResponse(0), fESDtrackCuts(0), fESDtrackCutsITSsa2010(0), fESDtrackCutsGlobal2015(0), fUtils(0), fRand(0),

//---> Flags controlling Event Tree output
fkSaveEventTree    ( kTRUE ), //no downscaling in this tree so far

//---> Flags controlling V0 TTree output
fkSaveV0Tree       ( kTRUE ),
fkDownScaleV0      ( kTRUE  ),
fDownScaleFactorV0 ( 0.001  ),
fkPreselectDedx ( kFALSE ),
fkUseOnTheFlyV0Cascading( kFALSE ),
fkDebugWrongPIDForTracking ( kFALSE ),
fkDebugBump(kFALSE),
fkDebugOOBPileup(kFALSE),
fkDoExtraEvSels(kTRUE),

//---> Flags controlling Cascade TTree output
fkSaveCascadeTree       ( kTRUE  ),
fkDownScaleCascade      ( kTRUE  ),
fDownScaleFactorCascade ( 0.001  ),

//---> Flags controlling Vertexers
fkRunVertexers    ( kFALSE ),
fkUseLightVertexer ( kTRUE ),
fkDoV0Refit       ( kTRUE ),
fkExtraCleanup    ( kTRUE ),

//---> Flag controlling trigger selection
fTrigType(AliVEvent::kMB),

//---> Min pT to save candidate
fMinPtToSave( 0.55 ),
fMaxPtToSave( 100.00 ),

//---> Variables for fTreeEvent
fCentrality_V0A(0),
fCentrality_V0C(0),
fCentrality_V0M(0),
fMVPileupFlag(kFALSE),
fOOBPileupFlag(kFALSE),
fNTOFClusters(-1),
fNTOFMatches(-1),
fNTracksITSsa2010(-1),
fNTracksGlobal2015(-1),
fNTracksGlobal2015TriggerPP(-1),
fAmplitudeV0A(-1.),
fAmplitudeV0C(-1.),
fClosestNonEmptyBC(-1),

//---> Variables for fTreeV0
fTreeVariableChi2V0(0),
fTreeVariableDcaV0Daughters(0),
fTreeVariableDcaV0ToPrimVertex(0),
fTreeVariableDcaPosToPrimVertex(0),
fTreeVariableDcaNegToPrimVertex(0),
fTreeVariableV0CosineOfPointingAngle(0),
fTreeVariableV0Radius(0),
fTreeVariablePt(0),
fTreeVariableRapK0Short(0),
fTreeVariableRapLambda(0),
fTreeVariableInvMassK0s(0),
fTreeVariableInvMassLambda(0),
fTreeVariableInvMassAntiLambda(0),
fTreeVariableAlphaV0(0),
fTreeVariablePtArmV0(0),
fTreeVariableNegEta(0),
fTreeVariablePosEta(0),
fTreeVariableMaxChi2PerCluster(0),
fTreeVariableMinTrackLength(0),

fTreeVariableNSigmasPosProton(0),
fTreeVariableNSigmasPosPion(0),
fTreeVariableNSigmasNegProton(0),
fTreeVariableNSigmasNegPion(0),

fTreeVariablePosPIDForTracking(-1),
fTreeVariableNegPIDForTracking(-1),
fTreeVariablePosdEdx(-1),
fTreeVariableNegdEdx(-1),
fTreeVariablePosInnerP(-1),
fTreeVariableNegInnerP(-1),
fTreeVariableNegTrackStatus(0),
fTreeVariablePosTrackStatus(0),
fTreeVariableNegDCAz(-1),
fTreeVariablePosDCAz(-1),

fTreeVariablePosITSClusters0(0),
fTreeVariablePosITSClusters1(0),
fTreeVariablePosITSClusters2(0),
fTreeVariablePosITSClusters3(0),
fTreeVariablePosITSClusters4(0),
fTreeVariablePosITSClusters5(0),

fTreeVariableNegITSClusters0(0),
fTreeVariableNegITSClusters1(0),
fTreeVariableNegITSClusters2(0),
fTreeVariableNegITSClusters3(0),
fTreeVariableNegITSClusters4(0),
fTreeVariableNegITSClusters5(0),

fTreeVariablePosITSSharedClusters0(0),
fTreeVariablePosITSSharedClusters1(0),
fTreeVariablePosITSSharedClusters2(0),
fTreeVariablePosITSSharedClusters3(0),
fTreeVariablePosITSSharedClusters4(0),
fTreeVariablePosITSSharedClusters5(0),

fTreeVariableNegITSSharedClusters0(0),
fTreeVariableNegITSSharedClusters1(0),
fTreeVariableNegITSSharedClusters2(0),
fTreeVariableNegITSSharedClusters3(0),
fTreeVariableNegITSSharedClusters4(0),
fTreeVariableNegITSSharedClusters5(0),

fTreeVariableDistOverTotMom(0),
fTreeVariableLeastNbrCrossedRows(0),
fTreeVariableLeastRatioCrossedRowsOverFindable(0),

fTreeVariableNegTOFExpTDiff(99999),
fTreeVariablePosTOFExpTDiff(99999),
fTreeVariableNegTOFSignal(99999),
fTreeVariablePosTOFSignal(99999),
fTreeVariableNegTOFBCid(99999),
fTreeVariablePosTOFBCid(99999),
fTreeVariableAmplitudeV0A(-1.),
fTreeVariableAmplitudeV0C(-1.),
fTreeVariableClosestNonEmptyBC(-1),

fTreeVariableCentrality_V0A(0),
fTreeVariableCentrality_V0C(0),
fTreeVariableCentrality_V0M(0),
fTreeVariableMVPileupFlag(kFALSE),
fTreeVariableOOBPileupFlag(kFALSE),

//---> Variables for fTreeCascade
fTreeCascVarCharge(0),
fTreeCascVarMassAsXi(0),
fTreeCascVarMassAsOmega(0),
fTreeCascVarPt(0),
fTreeCascVarRapXi(0),
fTreeCascVarRapOmega(0),
fTreeCascVarNegEta(0),
fTreeCascVarPosEta(0),
fTreeCascVarBachEta(0),
fTreeCascVarDCACascDaughters(0),
fTreeCascVarDCABachToPrimVtx(0),
fTreeCascVarDCAV0Daughters(0),
fTreeCascVarDCAV0ToPrimVtx(0),
fTreeCascVarDCAPosToPrimVtx(0),
fTreeCascVarDCANegToPrimVtx(0),
fTreeCascVarCascCosPointingAngle(0),
fTreeCascVarCascRadius(0),
fTreeCascVarV0Mass(0),
fTreeCascVarV0MassLambda(0),
fTreeCascVarV0MassAntiLambda(0),
fTreeCascVarV0CosPointingAngle(0),
fTreeCascVarV0CosPointingAngleSpecial(0),
fTreeCascVarV0Radius(0),
fTreeCascVarDCABachToBaryon(0),
fTreeCascVarWrongCosPA(0),
fTreeCascVarLeastNbrClusters(0),
fTreeCascVarDistOverTotMom(0),
fTreeCascVarNegNSigmaPion(0),
fTreeCascVarNegNSigmaProton(0),
fTreeCascVarPosNSigmaPion(0),
fTreeCascVarPosNSigmaProton(0),
fTreeCascVarBachNSigmaPion(0),
fTreeCascVarBachNSigmaKaon(0),

fTreeCascVarNegTOFNSigmaPion(0),
fTreeCascVarNegTOFNSigmaProton(0),
fTreeCascVarPosTOFNSigmaPion(0),
fTreeCascVarPosTOFNSigmaProton(0),
fTreeCascVarBachTOFNSigmaPion(0),
fTreeCascVarBachTOFNSigmaKaon(0),

fTreeCascVarNegITSNSigmaPion(0),
fTreeCascVarNegITSNSigmaProton(0),
fTreeCascVarPosITSNSigmaPion(0),
fTreeCascVarPosITSNSigmaProton(0),
fTreeCascVarBachITSNSigmaPion(0),
fTreeCascVarBachITSNSigmaKaon(0),
fTreeCascVarChiSquareV0(1e+3),
fTreeCascVarChiSquareCascade(1e+3),

fTreeCascVarBachDCAPVSigmaX2(0),
fTreeCascVarBachDCAPVSigmaY2(0),
fTreeCascVarBachDCAPVSigmaZ2(0),
fTreeCascVarPosDCAPVSigmaX2(0),
fTreeCascVarPosDCAPVSigmaY2(0),
fTreeCascVarPosDCAPVSigmaZ2(0),
fTreeCascVarNegDCAPVSigmaX2(0),
fTreeCascVarNegDCAPVSigmaY2(0),
fTreeCascVarNegDCAPVSigmaZ2(0),

fTreeCascVarPosPIDForTracking(-1),
fTreeCascVarNegPIDForTracking(-1),
fTreeCascVarBachPIDForTracking(-1),
fTreeCascVarPosdEdx(-1),
fTreeCascVarNegdEdx(-1),
fTreeCascVarBachdEdx(-1),
fTreeCascVarPosInnerP(-1),
fTreeCascVarNegInnerP(-1),
fTreeCascVarBachInnerP(-1),
fTreeCascVarNegTrackStatus(0), //!
fTreeCascVarPosTrackStatus(0), //!
fTreeCascVarBachTrackStatus(0), //!
fTreeCascVarNegDCAz(-1),
fTreeCascVarPosDCAz(-1),
fTreeCascVarBachDCAz(-1),
//fTreeCascVarPosTotMom(-1),
//fTreeCascVarNegTotMom(-1),
//fTreeCascVarBachTotMom(-1),
fTreeCascVarMaxChi2PerCluster(0),
fTreeCascVarMinTrackLength(0),

fTreeCascVarPosITSClusters0(0),
fTreeCascVarPosITSClusters1(0),
fTreeCascVarPosITSClusters2(0),
fTreeCascVarPosITSClusters3(0),
fTreeCascVarPosITSClusters4(0),
fTreeCascVarPosITSClusters5(0),

fTreeCascVarNegITSClusters0(0),
fTreeCascVarNegITSClusters1(0),
fTreeCascVarNegITSClusters2(0),
fTreeCascVarNegITSClusters3(0),
fTreeCascVarNegITSClusters4(0),
fTreeCascVarNegITSClusters5(0),

fTreeCascVarBachITSClusters0(0),
fTreeCascVarBachITSClusters1(0),
fTreeCascVarBachITSClusters2(0),
fTreeCascVarBachITSClusters3(0),
fTreeCascVarBachITSClusters4(0),
fTreeCascVarBachITSClusters5(0),

fTreeCascVarPosITSSharedClusters0(0),
fTreeCascVarPosITSSharedClusters1(0),
fTreeCascVarPosITSSharedClusters2(0),
fTreeCascVarPosITSSharedClusters3(0),
fTreeCascVarPosITSSharedClusters4(0),
fTreeCascVarPosITSSharedClusters5(0),

fTreeCascVarNegITSSharedClusters0(0),
fTreeCascVarNegITSSharedClusters1(0),
fTreeCascVarNegITSSharedClusters2(0),
fTreeCascVarNegITSSharedClusters3(0),
fTreeCascVarNegITSSharedClusters4(0),
fTreeCascVarNegITSSharedClusters5(0),

fTreeCascVarBachITSSharedClusters0(0),
fTreeCascVarBachITSSharedClusters1(0),
fTreeCascVarBachITSSharedClusters2(0),
fTreeCascVarBachITSSharedClusters3(0),
fTreeCascVarBachITSSharedClusters4(0),
fTreeCascVarBachITSSharedClusters5(0),

//Variables for debugging the invariant mass bump
//Full momentum information
fTreeCascVarNegPx(0),
fTreeCascVarNegPy(0),
fTreeCascVarNegPz(0),
fTreeCascVarPosPx(0),
fTreeCascVarPosPy(0),
fTreeCascVarPosPz(0),
fTreeCascVarBachPx(0),
fTreeCascVarBachPy(0),
fTreeCascVarBachPz(0),
fTreeCascVarV0DecayX(0),
fTreeCascVarV0DecayY(0),
fTreeCascVarV0DecayZ(0),
fTreeCascVarCascadeDecayX(0),
fTreeCascVarCascadeDecayY(0),
fTreeCascVarCascadeDecayZ(0),
fTreeCascVarV0Lifetime(0),
//Track Labels (check for duplicates, etc)
fTreeCascVarNegIndex(0),
fTreeCascVarPosIndex(0),
fTreeCascVarBachIndex(0),
//Event Number (check same-event index mixups)
fTreeCascVarEventNumber(0),

fTreeCascVarNegTOFExpTDiff(99999),
fTreeCascVarPosTOFExpTDiff(99999),
fTreeCascVarBachTOFExpTDiff(99999),
fTreeCascVarNegTOFSignal(99999),
fTreeCascVarPosTOFSignal(99999),
fTreeCascVarBachTOFSignal(99999),
fTreeCascVarNegTOFBCid(99999),
fTreeCascVarPosTOFBCid(99999),
fTreeCascVarBachTOFBCid(99999),
fTreeCascVarAmplitudeV0A(-1.),
fTreeCascVarAmplitudeV0C(-1.),
fTreeCascVarClosestNonEmptyBC(-1),

fTreeCascVarCentrality_V0A(0),
fTreeCascVarCentrality_V0C(0),
fTreeCascVarCentrality_V0M(0),
fTreeCascVarMVPileupFlag(kFALSE),
fTreeCascVarOOBPileupFlag(kFALSE),

//Kink tagging
fTreeCascVarBachIsKink(kFALSE),
fTreeCascVarPosIsKink(kFALSE),
fTreeCascVarNegIsKink(kFALSE),

//Histos
fHistEventCounter(0),
fHistCentrality_V0A(0),
fHistCentrality_V0C(0),
fHistCentrality_V0M(0)

//------------------------------------------------
// Tree Variables
{
    
}

AliAnalysisTaskStrangenessVsMultiplicityRun2pPb::AliAnalysisTaskStrangenessVsMultiplicityRun2pPb(Bool_t lSaveEventTree, Bool_t lSaveV0Tree, Bool_t lSaveCascadeTree, const char *name, TString lExtraOptions)
: AliAnalysisTaskSE(name), fListHist(0), fListV0(0), fListCascade(0), fTreeEvent(0), fTreeV0(0), fTreeCascade(0), fPIDResponse(0), fESDtrackCuts(0), fESDtrackCutsITSsa2010(0), fESDtrackCutsGlobal2015(0), fUtils(0), fRand(0),

//---> Flags controlling Event Tree output
fkSaveEventTree    ( kTRUE ), //no downscaling in this tree so far

//---> Flags controlling V0 TTree output
fkSaveV0Tree       ( kTRUE ),
fkDownScaleV0      ( kFALSE  ),
fDownScaleFactorV0 ( 0.001  ),
fkPreselectDedx ( kFALSE ),
fkUseOnTheFlyV0Cascading( kFALSE ),
fkDebugWrongPIDForTracking ( kFALSE ), //also for cascades...
fkDebugBump( kFALSE ),
fkDebugOOBPileup(kFALSE),
fkDoExtraEvSels(kTRUE),

//---> Flags controlling Cascade TTree output
fkSaveCascadeTree       ( kTRUE  ),
fkDownScaleCascade      ( kFALSE  ),
fDownScaleFactorCascade ( 0.001  ),

//---> Flags controlling Vertexers
fkRunVertexers    ( kFALSE ),
fkUseLightVertexer ( kTRUE ),
fkDoV0Refit       ( kTRUE ),
fkExtraCleanup    ( kTRUE ),

//---> Flag controlling trigger selection
fTrigType(AliVEvent::kMB),

//---> Min pT to save candidate
fMinPtToSave( 0.55 ),
fMaxPtToSave( 100.00 ),

//---> Variables for fTreeEvent
fCentrality_V0A(0),
fCentrality_V0C(0),
fCentrality_V0M(0),
fMVPileupFlag(kFALSE),
fOOBPileupFlag(kFALSE),
fNTOFClusters(-1),
fNTOFMatches(-1),
fNTracksITSsa2010(-1),
fNTracksGlobal2015(-1),
fNTracksGlobal2015TriggerPP(-1),
fAmplitudeV0A(-1.),
fAmplitudeV0C(-1.),
fClosestNonEmptyBC(-1),

//---> Variables for fTreeV0
fTreeVariableChi2V0(0),
fTreeVariableDcaV0Daughters(0),
fTreeVariableDcaV0ToPrimVertex(0),
fTreeVariableDcaPosToPrimVertex(0),
fTreeVariableDcaNegToPrimVertex(0),
fTreeVariableV0CosineOfPointingAngle(0),
fTreeVariableV0Radius(0),
fTreeVariablePt(0),
fTreeVariableRapK0Short(0),
fTreeVariableRapLambda(0),
fTreeVariableInvMassK0s(0),
fTreeVariableInvMassLambda(0),
fTreeVariableInvMassAntiLambda(0),
fTreeVariableAlphaV0(0),
fTreeVariablePtArmV0(0),
fTreeVariableNegEta(0),
fTreeVariablePosEta(0),
fTreeVariableMaxChi2PerCluster(0),
fTreeVariableMinTrackLength(0),

fTreeVariableNSigmasPosProton(0),
fTreeVariableNSigmasPosPion(0),
fTreeVariableNSigmasNegProton(0),
fTreeVariableNSigmasNegPion(0),

fTreeVariablePosPIDForTracking(-1),
fTreeVariableNegPIDForTracking(-1),
fTreeVariablePosdEdx(-1),
fTreeVariableNegdEdx(-1),
fTreeVariablePosInnerP(-1),
fTreeVariableNegInnerP(-1),
fTreeVariableNegTrackStatus(0),
fTreeVariablePosTrackStatus(0),
fTreeVariableNegDCAz(-1),
fTreeVariablePosDCAz(-1),

fTreeVariablePosITSClusters0(0),
fTreeVariablePosITSClusters1(0),
fTreeVariablePosITSClusters2(0),
fTreeVariablePosITSClusters3(0),
fTreeVariablePosITSClusters4(0),
fTreeVariablePosITSClusters5(0),

fTreeVariableNegITSClusters0(0),
fTreeVariableNegITSClusters1(0),
fTreeVariableNegITSClusters2(0),
fTreeVariableNegITSClusters3(0),
fTreeVariableNegITSClusters4(0),
fTreeVariableNegITSClusters5(0),

fTreeVariablePosITSSharedClusters0(0),
fTreeVariablePosITSSharedClusters1(0),
fTreeVariablePosITSSharedClusters2(0),
fTreeVariablePosITSSharedClusters3(0),
fTreeVariablePosITSSharedClusters4(0),
fTreeVariablePosITSSharedClusters5(0),

fTreeVariableNegITSSharedClusters0(0),
fTreeVariableNegITSSharedClusters1(0),
fTreeVariableNegITSSharedClusters2(0),
fTreeVariableNegITSSharedClusters3(0),
fTreeVariableNegITSSharedClusters4(0),
fTreeVariableNegITSSharedClusters5(0),

fTreeVariableDistOverTotMom(0),
fTreeVariableLeastNbrCrossedRows(0),
fTreeVariableLeastRatioCrossedRowsOverFindable(0),

fTreeVariableNegTOFExpTDiff(99999),
fTreeVariablePosTOFExpTDiff(99999),
fTreeVariableNegTOFSignal(99999),
fTreeVariablePosTOFSignal(99999),
fTreeVariableNegTOFBCid(99999),
fTreeVariablePosTOFBCid(99999),
fTreeVariableAmplitudeV0A(-1.),
fTreeVariableAmplitudeV0C(-1.),
fTreeVariableClosestNonEmptyBC(-1),

fTreeVariableCentrality_V0A(0),
fTreeVariableCentrality_V0C(0),
fTreeVariableCentrality_V0M(0),
fTreeVariableMVPileupFlag(kFALSE),
fTreeVariableOOBPileupFlag(kFALSE),

//---> Variables for fTreeCascade
fTreeCascVarCharge(0),
fTreeCascVarMassAsXi(0),
fTreeCascVarMassAsOmega(0),
fTreeCascVarPt(0),
fTreeCascVarRapXi(0),
fTreeCascVarRapOmega(0),
fTreeCascVarNegEta(0),
fTreeCascVarPosEta(0),
fTreeCascVarBachEta(0),
fTreeCascVarDCACascDaughters(0),
fTreeCascVarDCABachToPrimVtx(0),
fTreeCascVarDCAV0Daughters(0),
fTreeCascVarDCAV0ToPrimVtx(0),
fTreeCascVarDCAPosToPrimVtx(0),
fTreeCascVarDCANegToPrimVtx(0),
fTreeCascVarCascCosPointingAngle(0),
fTreeCascVarCascRadius(0),
fTreeCascVarV0Mass(0),
fTreeCascVarV0MassLambda(0),
fTreeCascVarV0MassAntiLambda(0),
fTreeCascVarV0CosPointingAngle(0),
fTreeCascVarV0CosPointingAngleSpecial(0),
fTreeCascVarV0Radius(0),
fTreeCascVarDCABachToBaryon(0),
fTreeCascVarWrongCosPA(0),
fTreeCascVarLeastNbrClusters(0),
fTreeCascVarDistOverTotMom(0),
fTreeCascVarNegNSigmaPion(0),
fTreeCascVarNegNSigmaProton(0),
fTreeCascVarPosNSigmaPion(0),
fTreeCascVarPosNSigmaProton(0),
fTreeCascVarBachNSigmaPion(0),
fTreeCascVarBachNSigmaKaon(0),

fTreeCascVarNegITSNSigmaPion(0),
fTreeCascVarNegITSNSigmaProton(0),
fTreeCascVarPosITSNSigmaPion(0),
fTreeCascVarPosITSNSigmaProton(0),
fTreeCascVarBachITSNSigmaPion(0),
fTreeCascVarBachITSNSigmaKaon(0),

fTreeCascVarNegTOFNSigmaPion(0),
fTreeCascVarNegTOFNSigmaProton(0),
fTreeCascVarPosTOFNSigmaPion(0),
fTreeCascVarPosTOFNSigmaProton(0),
fTreeCascVarBachTOFNSigmaPion(0),
fTreeCascVarBachTOFNSigmaKaon(0),

fTreeCascVarChiSquareV0(1e+3),
fTreeCascVarChiSquareCascade(1e+3),

fTreeCascVarBachDCAPVSigmaX2(0),
fTreeCascVarBachDCAPVSigmaY2(0),
fTreeCascVarBachDCAPVSigmaZ2(0),
fTreeCascVarPosDCAPVSigmaX2(0),
fTreeCascVarPosDCAPVSigmaY2(0),
fTreeCascVarPosDCAPVSigmaZ2(0),
fTreeCascVarNegDCAPVSigmaX2(0),
fTreeCascVarNegDCAPVSigmaY2(0),
fTreeCascVarNegDCAPVSigmaZ2(0),

fTreeCascVarPosPIDForTracking(-1),
fTreeCascVarNegPIDForTracking(-1),
fTreeCascVarBachPIDForTracking(-1),
fTreeCascVarPosdEdx(-1),
fTreeCascVarNegdEdx(-1),
fTreeCascVarBachdEdx(-1),
fTreeCascVarPosInnerP(-1),
fTreeCascVarNegInnerP(-1),
fTreeCascVarBachInnerP(-1),
fTreeCascVarNegTrackStatus(0), //!
fTreeCascVarPosTrackStatus(0), //!
fTreeCascVarBachTrackStatus(0), //!
fTreeCascVarNegDCAz(-1),
fTreeCascVarPosDCAz(-1),
fTreeCascVarBachDCAz(-1),
//fTreeCascVarPosTotMom(-1),
//fTreeCascVarNegTotMom(-1),
//fTreeCascVarBachTotMom(-1),
fTreeCascVarMaxChi2PerCluster(0),
fTreeCascVarMinTrackLength(0),

fTreeCascVarPosITSClusters0(0),
fTreeCascVarPosITSClusters1(0),
fTreeCascVarPosITSClusters2(0),
fTreeCascVarPosITSClusters3(0),
fTreeCascVarPosITSClusters4(0),
fTreeCascVarPosITSClusters5(0),

fTreeCascVarNegITSClusters0(0),
fTreeCascVarNegITSClusters1(0),
fTreeCascVarNegITSClusters2(0),
fTreeCascVarNegITSClusters3(0),
fTreeCascVarNegITSClusters4(0),
fTreeCascVarNegITSClusters5(0),

fTreeCascVarBachITSClusters0(0),
fTreeCascVarBachITSClusters1(0),
fTreeCascVarBachITSClusters2(0),
fTreeCascVarBachITSClusters3(0),
fTreeCascVarBachITSClusters4(0),
fTreeCascVarBachITSClusters5(0),

fTreeCascVarPosITSSharedClusters0(0),
fTreeCascVarPosITSSharedClusters1(0),
fTreeCascVarPosITSSharedClusters2(0),
fTreeCascVarPosITSSharedClusters3(0),
fTreeCascVarPosITSSharedClusters4(0),
fTreeCascVarPosITSSharedClusters5(0),

fTreeCascVarNegITSSharedClusters0(0),
fTreeCascVarNegITSSharedClusters1(0),
fTreeCascVarNegITSSharedClusters2(0),
fTreeCascVarNegITSSharedClusters3(0),
fTreeCascVarNegITSSharedClusters4(0),
fTreeCascVarNegITSSharedClusters5(0),

fTreeCascVarBachITSSharedClusters0(0),
fTreeCascVarBachITSSharedClusters1(0),
fTreeCascVarBachITSSharedClusters2(0),
fTreeCascVarBachITSSharedClusters3(0),
fTreeCascVarBachITSSharedClusters4(0),
fTreeCascVarBachITSSharedClusters5(0),

//Variables for debugging the invariant mass bump
//Full momentum information
fTreeCascVarNegPx(0),
fTreeCascVarNegPy(0),
fTreeCascVarNegPz(0),
fTreeCascVarPosPx(0),
fTreeCascVarPosPy(0),
fTreeCascVarPosPz(0),
fTreeCascVarBachPx(0),
fTreeCascVarBachPy(0),
fTreeCascVarBachPz(0),
fTreeCascVarV0DecayX(0),
fTreeCascVarV0DecayY(0),
fTreeCascVarV0DecayZ(0),
fTreeCascVarCascadeDecayX(0),
fTreeCascVarCascadeDecayY(0),
fTreeCascVarCascadeDecayZ(0),
fTreeCascVarV0Lifetime(0),
//Track Labels (check for duplicates, etc)
fTreeCascVarNegIndex(0),
fTreeCascVarPosIndex(0),
fTreeCascVarBachIndex(0),
//Event Number (check same-event index mixups)
fTreeCascVarEventNumber(0),

fTreeCascVarNegTOFExpTDiff(99999),
fTreeCascVarPosTOFExpTDiff(99999),
fTreeCascVarBachTOFExpTDiff(99999),
fTreeCascVarNegTOFSignal(99999),
fTreeCascVarPosTOFSignal(99999),
fTreeCascVarBachTOFSignal(99999),
fTreeCascVarNegTOFBCid(99999),
fTreeCascVarPosTOFBCid(99999),
fTreeCascVarBachTOFBCid(99999),
fTreeCascVarAmplitudeV0A(-1.),
fTreeCascVarAmplitudeV0C(-1.),
fTreeCascVarClosestNonEmptyBC(-1),

fTreeCascVarCentrality_V0A(0),
fTreeCascVarCentrality_V0C(0),
fTreeCascVarCentrality_V0M(0),
fTreeCascVarMVPileupFlag(kFALSE),
fTreeCascVarOOBPileupFlag(kFALSE),

//Kink tagging
fTreeCascVarBachIsKink(kFALSE),
fTreeCascVarPosIsKink(kFALSE),
fTreeCascVarNegIsKink(kFALSE),

//Histos
fHistEventCounter(0),
fHistCentrality_V0A(0),
fHistCentrality_V0C(0),
fHistCentrality_V0M(0)



{
    
    //Re-vertex: Will only apply for cascade candidates
    
    fV0VertexerSels[0] =  33.  ;  // max allowed chi2
    fV0VertexerSels[1] =   0.02;  // min allowed impact parameter for the 1st daughter (LHC09a4 : 0.05)
    fV0VertexerSels[2] =   0.02;  // min allowed impact parameter for the 2nd daughter (LHC09a4 : 0.05)
    fV0VertexerSels[3] =   2.0 ;  // max allowed DCA between the daughter tracks       (LHC09a4 : 0.5)
    fV0VertexerSels[4] =   0.95;  // min allowed cosine of V0's pointing angle         (LHC09a4 : 0.99)
    fV0VertexerSels[5] =   1.0 ;  // min radius of the fiducial volume                 (LHC09a4 : 0.2)
    fV0VertexerSels[6] = 200.  ;  // max radius of the fiducial volume                 (LHC09a4 : 100.0)
    
    fCascadeVertexerSels[0] =  33.   ;  // max allowed chi2 (same as PDC07)
    fCascadeVertexerSels[1] =   0.05 ;  // min allowed V0 impact parameter                    (PDC07 : 0.05   / LHC09a4 : 0.025 )
    fCascadeVertexerSels[2] =   0.010;  // "window" around the Lambda mass                    (PDC07 : 0.008  / LHC09a4 : 0.010 )
    fCascadeVertexerSels[3] =   0.03 ;  // min allowed bachelor's impact parameter            (PDC07 : 0.035  / LHC09a4 : 0.025 )
    fCascadeVertexerSels[4] =   2.0  ;  // max allowed DCA between the V0 and the bachelor    (PDC07 : 0.1    / LHC09a4 : 0.2   )
    fCascadeVertexerSels[5] =   0.95 ;  // min allowed cosine of the cascade pointing angle   (PDC07 : 0.9985 / LHC09a4 : 0.998 )
    fCascadeVertexerSels[6] =   0.4  ;  // min radius of the fiducial volume                  (PDC07 : 0.9    / LHC09a4 : 0.2   )
    fCascadeVertexerSels[7] = 100.   ;  // max radius of the fiducial volume                  (PDC07 : 100    / LHC09a4 : 100   )
    
    //[0]+[1]*TMath::Exp([2]*x)+[3]*TMath::Exp([4]*x)
    fLambdaMassMean[0]=1.116; //standard fixed
    fLambdaMassMean[1]=0.0;
    fLambdaMassMean[2]=0.0;
    fLambdaMassMean[3]=0.0;
    fLambdaMassMean[4]=0.0;
    
    //[0]+[1]*x+[2]*TMath::Exp([3]*x)
    fLambdaMassSigma[0]=0.002; //standard at roughly the integ val
    fLambdaMassSigma[1]=0.0;
    fLambdaMassSigma[2]=0.0;
    fLambdaMassSigma[3]=0.0;
    
    fkSaveEventTree    = lSaveEventTree;
    fkSaveV0Tree       = lSaveV0Tree;
    fkSaveCascadeTree  = lSaveCascadeTree;
    
    //Standard output
    DefineOutput(1, TList::Class()); // Basic Histograms
    DefineOutput(2, TList::Class()); // V0 Histogram Output
    DefineOutput(3, TList::Class()); // Cascade Histogram Output
    
    //Optional output
    if (fkSaveEventTree)
        DefineOutput(4, TTree::Class()); // Event Tree output
    if (fkSaveV0Tree)
        DefineOutput(5, TTree::Class()); // V0 Tree output
    if (fkSaveCascadeTree)
        DefineOutput(6, TTree::Class()); // Cascade Tree output
    
    //Special Debug Options (more to be added as needed)
    // A - Study Wrong PID for tracking bug
    // B - Study invariant mass *B*ump
    // C - Study OOB pileup in pp 2016 data
    if ( lExtraOptions.Contains("A") ) fkDebugWrongPIDForTracking = kTRUE;
    if ( lExtraOptions.Contains("B") ) fkDebugBump                = kTRUE;
    if ( lExtraOptions.Contains("C") ) fkDebugOOBPileup           = kTRUE;
    
}


AliAnalysisTaskStrangenessVsMultiplicityRun2pPb::~AliAnalysisTaskStrangenessVsMultiplicityRun2pPb()
{
    //------------------------------------------------
    // DESTRUCTOR
    //------------------------------------------------
    
    //Destroy output objects if present
    if (fListHist) {
        delete fListHist;
        fListHist = 0x0;
    }
    if (fListV0) {
        delete fListV0;
        fListV0 = 0x0;
    }
    if (fListCascade) {
        delete fListCascade;
        fListCascade = 0x0;
    }
    if (fTreeEvent) {
        delete fTreeEvent;
        fTreeEvent = 0x0;
    }
    if (fTreeV0) {
        delete fTreeV0;
        fTreeV0 = 0x0;
    }
    if (fTreeCascade) {
        delete fTreeCascade;
        fTreeCascade = 0x0;
    }
    if (fUtils) {
        delete fUtils;
        fUtils = 0x0;
    }
    if (fRand) {
        delete fRand;
        fRand = 0x0;
    }
}

//________________________________________________________________________
void AliAnalysisTaskStrangenessVsMultiplicityRun2pPb::UserCreateOutputObjects()
{
    //------------------------------------------------
    // fTreeEvent: EbyE information
    //------------------------------------------------
    if(fkSaveEventTree){
        fTreeEvent = new TTree("fTreeEvent","Event");
        //Branch Definitions
        fTreeEvent->Branch("fCentrality_V0A",&fCentrality_V0A,"fCentrality_V0A/F");
        fTreeEvent->Branch("fCentrality_V0C",&fCentrality_V0C,"fCentrality_V0C/F");
        fTreeEvent->Branch("fCentrality_V0M",&fCentrality_V0M,"fCentrality_V0M/F");
        fTreeEvent->Branch("fMVPileupFlag",&fMVPileupFlag,"fMVPileupFlag/O");
        //
        if ( fkDebugOOBPileup ){
            fTreeEvent->Branch("fOOBPileupFlag",&fOOBPileupFlag,"fOOBPileupFlag/O");
            fTreeEvent->Branch("fNTOFClusters",&fNTOFClusters,"fNTOFClusters/I");
            fTreeEvent->Branch("fNTOFMatches",&fNTOFMatches,"fNTOFMatches/I");
            fTreeEvent->Branch("fNTracksITSsa2010",&fNTracksITSsa2010,"fNTracksITSsa2010/I");
            fTreeEvent->Branch("fNTracksGlobal2015",&fNTracksGlobal2015,"fNTracksGlobal2015/I");
            fTreeEvent->Branch("fNTracksGlobal2015TriggerPP",&fNTracksGlobal2015TriggerPP,"fNTracksGlobal2015TriggerPP/I");
            fTreeEvent->Branch("fAmplitudeV0A",&fAmplitudeV0A,"fAmplitudeV0A/F");
            fTreeEvent->Branch("fAmplitudeV0C",&fAmplitudeV0C,"fAmplitudeV0C/F");
            fTreeEvent->Branch("fClosestNonEmptyBC",&fClosestNonEmptyBC,"fClosestNonEmptyBC/I");
        }
    }
    
    //------------------------------------------------
    // fTreeV0: V0 Candidate Information
    //------------------------------------------------
    if(fkSaveV0Tree){
        //Create Basic V0 Output Tree
        fTreeV0 = new TTree( "fTreeV0", "V0 Candidates");
        //-----------BASIC-INFO---------------------------
        fTreeV0->Branch("fTreeVariableChi2V0",&fTreeVariableChi2V0,"fTreeVariableChi2V0/F");
        fTreeV0->Branch("fTreeVariableDcaV0Daughters",&fTreeVariableDcaV0Daughters,"fTreeVariableDcaV0Daughters/F");
        fTreeV0->Branch("fTreeVariableDcaV0ToPrimVertex",&fTreeVariableDcaV0ToPrimVertex,"fTreeVariableDcaV0ToPrimVertex/F");
        fTreeV0->Branch("fTreeVariableDcaPosToPrimVertex",&fTreeVariableDcaPosToPrimVertex,"fTreeVariableDcaPosToPrimVertex/F");
        fTreeV0->Branch("fTreeVariableDcaNegToPrimVertex",&fTreeVariableDcaNegToPrimVertex,"fTreeVariableDcaNegToPrimVertex/F");
        fTreeV0->Branch("fTreeVariableV0Radius",&fTreeVariableV0Radius,"fTreeVariableV0Radius/F");
        fTreeV0->Branch("fTreeVariablePt",&fTreeVariablePt,"fTreeVariablePt/F");
        fTreeV0->Branch("fTreeVariableRapK0Short",&fTreeVariableRapK0Short,"fTreeVariableRapK0Short/F");
        fTreeV0->Branch("fTreeVariableRapLambda",&fTreeVariableRapLambda,"fTreeVariableRapLambda/F");
        fTreeV0->Branch("fTreeVariableInvMassK0s",&fTreeVariableInvMassK0s,"fTreeVariableInvMassK0s/F");
        fTreeV0->Branch("fTreeVariableInvMassLambda",&fTreeVariableInvMassLambda,"fTreeVariableInvMassLambda/F");
        fTreeV0->Branch("fTreeVariableInvMassAntiLambda",&fTreeVariableInvMassAntiLambda,"fTreeVariableInvMassAntiLambda/F");
        fTreeV0->Branch("fTreeVariableV0CosineOfPointingAngle",&fTreeVariableV0CosineOfPointingAngle,"fTreeVariableV0CosineOfPointingAngle/F");
        fTreeV0->Branch("fTreeVariableAlphaV0",&fTreeVariableAlphaV0,"fTreeVariableAlphaV0/F");
        fTreeV0->Branch("fTreeVariablePtArmV0",&fTreeVariablePtArmV0,"fTreeVariablePtArmV0/F");
        fTreeV0->Branch("fTreeVariableLeastNbrCrossedRows",&fTreeVariableLeastNbrCrossedRows,"fTreeVariableLeastNbrCrossedRows/I");
        fTreeV0->Branch("fTreeVariableLeastRatioCrossedRowsOverFindable",&fTreeVariableLeastRatioCrossedRowsOverFindable,"fTreeVariableLeastRatioCrossedRowsOverFindable/F");
        fTreeV0->Branch("fTreeVariableMaxChi2PerCluster",&fTreeVariableMaxChi2PerCluster,"fTreeVariableMaxChi2PerCluster/F");
        fTreeV0->Branch("fTreeVariableMinTrackLength",&fTreeVariableMinTrackLength,"fTreeVariableMinTrackLength/F");
        fTreeV0->Branch("fTreeVariableDistOverTotMom",&fTreeVariableDistOverTotMom,"fTreeVariableDistOverTotMom/F");
        fTreeV0->Branch("fTreeVariableNSigmasPosProton",&fTreeVariableNSigmasPosProton,"fTreeVariableNSigmasPosProton/F");
        fTreeV0->Branch("fTreeVariableNSigmasPosPion",&fTreeVariableNSigmasPosPion,"fTreeVariableNSigmasPosPion/F");
        fTreeV0->Branch("fTreeVariableNSigmasNegProton",&fTreeVariableNSigmasNegProton,"fTreeVariableNSigmasNegProton/F");
        fTreeV0->Branch("fTreeVariableNSigmasNegPion",&fTreeVariableNSigmasNegPion,"fTreeVariableNSigmasNegPion/F");
        fTreeV0->Branch("fTreeVariableNegEta",&fTreeVariableNegEta,"fTreeVariableNegEta/F");
        fTreeV0->Branch("fTreeVariablePosEta",&fTreeVariablePosEta,"fTreeVariablePosEta/F");
        //-----------MULTIPLICITY-INFO--------------------
        fTreeV0->Branch("fTreeVariableCentrality_V0A",&fTreeVariableCentrality_V0A,"fTreeVariableCentrality_V0A/F");
        fTreeV0->Branch("fTreeVariableCentrality_V0C",&fTreeVariableCentrality_V0C,"fTreeVariableCentrality_V0C/F");
        fTreeV0->Branch("fTreeVariableCentrality_V0M",&fTreeVariableCentrality_V0M,"fTreeVariableCentrality_V0M/F");
        fTreeV0->Branch("fTreeVariableMVPileupFlag",&fTreeVariableMVPileupFlag,"fTreeVariableMVPileupFlag/O");
        //------------------------------------------------
        if ( fkDebugWrongPIDForTracking ){
            fTreeV0->Branch("fTreeVariablePosPIDForTracking",&fTreeVariablePosPIDForTracking,"fTreeVariablePosPIDForTracking/I");
            fTreeV0->Branch("fTreeVariableNegPIDForTracking",&fTreeVariableNegPIDForTracking,"fTreeVariableNegPIDForTracking/I");
            fTreeV0->Branch("fTreeVariablePosdEdx",&fTreeVariablePosdEdx,"fTreeVariablePosdEdx/F");
            fTreeV0->Branch("fTreeVariableNegdEdx",&fTreeVariableNegdEdx,"fTreeVariableNegdEdx/F");
            fTreeV0->Branch("fTreeVariablePosInnerP",&fTreeVariablePosInnerP,"fTreeVariablePosInnerP/F");
            fTreeV0->Branch("fTreeVariableNegInnerP",&fTreeVariableNegInnerP,"fTreeVariableNegInnerP/F");
            fTreeV0->Branch("fTreeVariableNegTrackStatus",&fTreeVariableNegTrackStatus,"fTreeVariableNegTrackStatus/l");
            fTreeV0->Branch("fTreeVariablePosTrackStatus",&fTreeVariablePosTrackStatus,"fTreeVariablePosTrackStatus/l");
            fTreeV0->Branch("fTreeVariableNegDCAz",&fTreeVariableNegDCAz,"fTreeVariableNegDCAz/F");
            fTreeV0->Branch("fTreeVariablePosDCAz",&fTreeVariablePosDCAz,"fTreeVariablePosDCAz/F");
        }
        if ( fkDebugOOBPileup ) {
            fTreeV0->Branch("fTreeVariablePosITSClusters0",&fTreeVariablePosITSClusters0,"fTreeVariablePosITSClusters0/O");
            fTreeV0->Branch("fTreeVariablePosITSClusters1",&fTreeVariablePosITSClusters1,"fTreeVariablePosITSClusters1/O");
            fTreeV0->Branch("fTreeVariablePosITSClusters2",&fTreeVariablePosITSClusters2,"fTreeVariablePosITSClusters2/O");
            fTreeV0->Branch("fTreeVariablePosITSClusters3",&fTreeVariablePosITSClusters3,"fTreeVariablePosITSClusters3/O");
            fTreeV0->Branch("fTreeVariablePosITSClusters4",&fTreeVariablePosITSClusters4,"fTreeVariablePosITSClusters4/O");
            fTreeV0->Branch("fTreeVariablePosITSClusters5",&fTreeVariablePosITSClusters5,"fTreeVariablePosITSClusters5/O");
            
            fTreeV0->Branch("fTreeVariableNegITSClusters0",&fTreeVariableNegITSClusters0,"fTreeVariableNegITSClusters0/O");
            fTreeV0->Branch("fTreeVariableNegITSClusters1",&fTreeVariableNegITSClusters1,"fTreeVariableNegITSClusters1/O");
            fTreeV0->Branch("fTreeVariableNegITSClusters2",&fTreeVariableNegITSClusters2,"fTreeVariableNegITSClusters2/O");
            fTreeV0->Branch("fTreeVariableNegITSClusters3",&fTreeVariableNegITSClusters3,"fTreeVariableNegITSClusters3/O");
            fTreeV0->Branch("fTreeVariableNegITSClusters4",&fTreeVariableNegITSClusters4,"fTreeVariableNegITSClusters4/O");
            fTreeV0->Branch("fTreeVariableNegITSClusters5",&fTreeVariableNegITSClusters5,"fTreeVariableNegITSClusters5/O");
            
            fTreeV0->Branch("fTreeVariablePosITSSharedClusters0",&fTreeVariablePosITSSharedClusters0,"fTreeVariablePosITSSharedClusters0/O");
            fTreeV0->Branch("fTreeVariablePosITSSharedClusters1",&fTreeVariablePosITSSharedClusters1,"fTreeVariablePosITSSharedClusters1/O");
            fTreeV0->Branch("fTreeVariablePosITSSharedClusters2",&fTreeVariablePosITSSharedClusters2,"fTreeVariablePosITSSharedClusters2/O");
            fTreeV0->Branch("fTreeVariablePosITSSharedClusters3",&fTreeVariablePosITSSharedClusters3,"fTreeVariablePosITSSharedClusters3/O");
            fTreeV0->Branch("fTreeVariablePosITSSharedClusters4",&fTreeVariablePosITSSharedClusters4,"fTreeVariablePosITSSharedClusters4/O");
            fTreeV0->Branch("fTreeVariablePosITSSharedClusters5",&fTreeVariablePosITSSharedClusters5,"fTreeVariablePosITSSharedClusters5/O");
            
            fTreeV0->Branch("fTreeVariableNegITSSharedClusters0",&fTreeVariableNegITSSharedClusters0,"fTreeVariableNegITSSharedClusters0/O");
            fTreeV0->Branch("fTreeVariableNegITSSharedClusters1",&fTreeVariableNegITSSharedClusters1,"fTreeVariableNegITSSharedClusters1/O");
            fTreeV0->Branch("fTreeVariableNegITSSharedClusters2",&fTreeVariableNegITSSharedClusters2,"fTreeVariableNegITSSharedClusters2/O");
            fTreeV0->Branch("fTreeVariableNegITSSharedClusters3",&fTreeVariableNegITSSharedClusters3,"fTreeVariableNegITSSharedClusters3/O");
            fTreeV0->Branch("fTreeVariableNegITSSharedClusters4",&fTreeVariableNegITSSharedClusters4,"fTreeVariableNegITSSharedClusters4/O");
            fTreeV0->Branch("fTreeVariableNegITSSharedClusters5",&fTreeVariableNegITSSharedClusters5,"fTreeVariableNegITSSharedClusters5/O");
            
            fTreeV0->Branch("fTreeVariableNegTOFExpTDiff",&fTreeVariableNegTOFExpTDiff,"fTreeVariableNegTOFExpTDiff/F");
            fTreeV0->Branch("fTreeVariablePosTOFExpTDiff",&fTreeVariablePosTOFExpTDiff,"fTreeVariablePosTOFExpTDiff/F");
            fTreeV0->Branch("fTreeVariableNegTOFSignal",&fTreeVariableNegTOFSignal,"fTreeVariableNegTOFSignal/F");
            fTreeV0->Branch("fTreeVariablePosTOFSignal",&fTreeVariablePosTOFSignal,"fTreeVariablePosTOFSignal/F");
            fTreeV0->Branch("fTreeVariableNegTOFBCid",&fTreeVariableNegTOFBCid,"fTreeVariableNegTOFBCid/I");
            fTreeV0->Branch("fTreeVariablePosTOFBCid",&fTreeVariablePosTOFBCid,"fTreeVariablePosTOFBCid/I");
            // Event info
            fTreeV0->Branch("fTreeVariableOOBPileupFlag",&fTreeVariableOOBPileupFlag,"fTreeVariableOOBPileupFlag/O");
            fTreeV0->Branch("fTreeVariableAmplitudeV0A",&fTreeVariableAmplitudeV0A,"fTreeVariableAmplitudeV0A/F");
            fTreeV0->Branch("fTreeVariableAmplitudeV0C",&fTreeVariableAmplitudeV0C,"fTreeVariableAmplitudeV0C/F");
            fTreeV0->Branch("fTreeVariableClosestNonEmptyBC",&fTreeVariableClosestNonEmptyBC,"fTreeVariableClosestNonEmptyBC/I");
        }
        //------------------------------------------------
    }
    
    //------------------------------------------------
    // fTreeCascade Branch definitions - Cascade Tree
    //------------------------------------------------
    if(fkSaveCascadeTree){
        //Create Cascade output tree
        fTreeCascade = new TTree("fTreeCascade","CascadeCandidates");
        //-----------BASIC-INFO---------------------------
        fTreeCascade->Branch("fTreeCascVarCharge",&fTreeCascVarCharge,"fTreeCascVarCharge/I");
        fTreeCascade->Branch("fTreeCascVarMassAsXi",&fTreeCascVarMassAsXi,"fTreeCascVarMassAsXi/F");
        fTreeCascade->Branch("fTreeCascVarMassAsOmega",&fTreeCascVarMassAsOmega,"fTreeCascVarMassAsOmega/F");
        fTreeCascade->Branch("fTreeCascVarPt",&fTreeCascVarPt,"fTreeCascVarPt/F");
        fTreeCascade->Branch("fTreeCascVarRapXi",&fTreeCascVarRapXi,"fTreeCascVarRapXi/F");
        fTreeCascade->Branch("fTreeCascVarRapOmega",&fTreeCascVarRapOmega,"fTreeCascVarRapOmega/F");
        fTreeCascade->Branch("fTreeCascVarNegEta",&fTreeCascVarNegEta,"fTreeCascVarNegEta/F");
        fTreeCascade->Branch("fTreeCascVarPosEta",&fTreeCascVarPosEta,"fTreeCascVarPosEta/F");
        fTreeCascade->Branch("fTreeCascVarBachEta",&fTreeCascVarBachEta,"fTreeCascVarBachEta/F");
        //-----------INFO-FOR-CUTS------------------------
        fTreeCascade->Branch("fTreeCascVarDCACascDaughters",&fTreeCascVarDCACascDaughters,"fTreeCascVarDCACascDaughters/F");
        fTreeCascade->Branch("fTreeCascVarDCABachToPrimVtx",&fTreeCascVarDCABachToPrimVtx,"fTreeCascVarDCABachToPrimVtx/F");
        fTreeCascade->Branch("fTreeCascVarDCAV0Daughters",&fTreeCascVarDCAV0Daughters,"fTreeCascVarDCAV0Daughters/F");
        fTreeCascade->Branch("fTreeCascVarDCAV0ToPrimVtx",&fTreeCascVarDCAV0ToPrimVtx,"fTreeCascVarDCAV0ToPrimVtx/F");
        fTreeCascade->Branch("fTreeCascVarDCAPosToPrimVtx",&fTreeCascVarDCAPosToPrimVtx,"fTreeCascVarDCAPosToPrimVtx/F");
        fTreeCascade->Branch("fTreeCascVarDCANegToPrimVtx",&fTreeCascVarDCANegToPrimVtx,"fTreeCascVarDCANegToPrimVtx/F");
        fTreeCascade->Branch("fTreeCascVarCascCosPointingAngle",&fTreeCascVarCascCosPointingAngle,"fTreeCascVarCascCosPointingAngle/F");
        fTreeCascade->Branch("fTreeCascVarCascRadius",&fTreeCascVarCascRadius,"fTreeCascVarCascRadius/F");
        fTreeCascade->Branch("fTreeCascVarV0Mass",&fTreeCascVarV0Mass,"fTreeCascVarV0Mass/F");
        fTreeCascade->Branch("fTreeCascVarV0MassLambda",&fTreeCascVarV0MassLambda,"fTreeCascVarV0MassLambda/F");
        fTreeCascade->Branch("fTreeCascVarV0MassAntiLambda",&fTreeCascVarV0MassAntiLambda,"fTreeCascVarV0MassAntiLambda/F");
        fTreeCascade->Branch("fTreeCascVarV0CosPointingAngle",&fTreeCascVarV0CosPointingAngle,"fTreeCascVarV0CosPointingAngle/F");
        fTreeCascade->Branch("fTreeCascVarV0CosPointingAngleSpecial",&fTreeCascVarV0CosPointingAngleSpecial,"fTreeCascVarV0CosPointingAngleSpecial/F");
        fTreeCascade->Branch("fTreeCascVarV0Radius",&fTreeCascVarV0Radius,"fTreeCascVarV0Radius/F");
        fTreeCascade->Branch("fTreeCascVarDCABachToBaryon",&fTreeCascVarDCABachToBaryon,"fTreeCascVarDCABachToBaryon/F");
        fTreeCascade->Branch("fTreeCascVarWrongCosPA",&fTreeCascVarWrongCosPA,"fTreeCascVarWrongCosPA/F");
        fTreeCascade->Branch("fTreeCascVarLeastNbrClusters",&fTreeCascVarLeastNbrClusters,"fTreeCascVarLeastNbrClusters/I");
        fTreeCascade->Branch("fTreeCascVarMaxChi2PerCluster",&fTreeCascVarMaxChi2PerCluster,"fTreeCascVarMaxChi2PerCluster/F");
        fTreeCascade->Branch("fTreeCascVarMinTrackLength",&fTreeCascVarMinTrackLength,"fTreeCascVarMinTrackLength/F");
        //-----------MULTIPLICITY-INFO--------------------
        fTreeCascade->Branch("fTreeCascVarCentrality_V0A",&fTreeCascVarCentrality_V0A,"fTreeCascVarCentrality_V0A/F");
        fTreeCascade->Branch("fTreeCascVarCentrality_V0C",&fTreeCascVarCentrality_V0C,"fTreeCascVarCentrality_V0C/F");
        fTreeCascade->Branch("fTreeCascVarCentrality_V0M",&fTreeCascVarCentrality_V0M,"fTreeCascVarCentrality_V0M/F");
        fTreeCascade->Branch("fTreeCascVarMVPileupFlag",&fTreeCascVarMVPileupFlag,"fTreeCascVarMVPileupFlag/O");
        //-----------DECAY-LENGTH-INFO--------------------
        fTreeCascade->Branch("fTreeCascVarDistOverTotMom",&fTreeCascVarDistOverTotMom,"fTreeCascVarDistOverTotMom/F");
        //------------------------------------------------
        fTreeCascade->Branch("fTreeCascVarNegNSigmaPion",&fTreeCascVarNegNSigmaPion,"fTreeCascVarNegNSigmaPion/F");
        fTreeCascade->Branch("fTreeCascVarNegNSigmaProton",&fTreeCascVarNegNSigmaProton,"fTreeCascVarNegNSigmaProton/F");
        fTreeCascade->Branch("fTreeCascVarPosNSigmaPion",&fTreeCascVarPosNSigmaPion,"fTreeCascVarPosNSigmaPion/F");
        fTreeCascade->Branch("fTreeCascVarPosNSigmaProton",&fTreeCascVarPosNSigmaProton,"fTreeCascVarPosNSigmaProton/F");
        fTreeCascade->Branch("fTreeCascVarBachNSigmaPion",&fTreeCascVarBachNSigmaPion,"fTreeCascVarBachNSigmaPion/F");
        fTreeCascade->Branch("fTreeCascVarBachNSigmaKaon",&fTreeCascVarBachNSigmaKaon,"fTreeCascVarBachNSigmaKaon/F");
        
        fTreeCascade->Branch("fTreeCascVarNegTOFNSigmaPion",&fTreeCascVarNegTOFNSigmaPion,"fTreeCascVarNegTOFNSigmaPion/F");
        fTreeCascade->Branch("fTreeCascVarNegTOFNSigmaProton",&fTreeCascVarNegTOFNSigmaProton,"fTreeCascVarTOFNegNSigmaProton/F");
        fTreeCascade->Branch("fTreeCascVarPosTOFNSigmaPion",&fTreeCascVarPosTOFNSigmaPion,"fTreeCascVarPosTOFNSigmaPion/F");
        fTreeCascade->Branch("fTreeCascVarPosTOFNSigmaProton",&fTreeCascVarPosTOFNSigmaProton,"fTreeCascVarPosTOFNSigmaProton/F");
        fTreeCascade->Branch("fTreeCascVarBachTOFNSigmaPion",&fTreeCascVarBachTOFNSigmaPion,"fTreeCascVarBachTOFNSigmaPion/F");
        fTreeCascade->Branch("fTreeCascVarBachTOFNSigmaKaon",&fTreeCascVarBachTOFNSigmaKaon,"fTreeCascVarBachTOFNSigmaKaon/F");
        
        fTreeCascade->Branch("fTreeCascVarNegITSNSigmaPion",&fTreeCascVarNegITSNSigmaPion,"fTreeCascVarNegITSNSigmaPion/F");
        fTreeCascade->Branch("fTreeCascVarNegITSNSigmaProton",&fTreeCascVarNegITSNSigmaProton,"fTreeCascVarITSNegNSigmaProton/F");
        fTreeCascade->Branch("fTreeCascVarPosITSNSigmaPion",&fTreeCascVarPosITSNSigmaPion,"fTreeCascVarPosITSNSigmaPion/F");
        fTreeCascade->Branch("fTreeCascVarPosITSNSigmaProton",&fTreeCascVarPosITSNSigmaProton,"fTreeCascVarPosITSNSigmaProton/F");
        fTreeCascade->Branch("fTreeCascVarBachITSNSigmaPion",&fTreeCascVarBachITSNSigmaPion,"fTreeCascVarBachITSNSigmaPion/F");
        fTreeCascade->Branch("fTreeCascVarBachITSNSigmaKaon",&fTreeCascVarBachITSNSigmaKaon,"fTreeCascVarBachITSNSigmaKaon/F");
        
        //------------------------------------------------
        fTreeCascade->Branch("fTreeCascVarChiSquareV0",&fTreeCascVarChiSquareV0,"fTreeCascVarChiSquareV0/F");
        fTreeCascade->Branch("fTreeCascVarChiSquareCascade",&fTreeCascVarChiSquareCascade,"fTreeCascVarChiSquareCascade/F");
        //------------------------------------------------
        if ( fkDebugWrongPIDForTracking ){
            fTreeCascade->Branch("fTreeCascVarPosPIDForTracking",&fTreeCascVarPosPIDForTracking,"fTreeCascVarPosPIDForTracking/I");
            fTreeCascade->Branch("fTreeCascVarNegPIDForTracking",&fTreeCascVarNegPIDForTracking,"fTreeCascVarNegPIDForTracking/I");
            fTreeCascade->Branch("fTreeCascVarBachPIDForTracking",&fTreeCascVarBachPIDForTracking,"fTreeCascVarBachPIDForTracking/I");
            fTreeCascade->Branch("fTreeCascVarPosdEdx",&fTreeCascVarPosdEdx,"fTreeCascVarPosdEdx/F");
            fTreeCascade->Branch("fTreeCascVarNegdEdx",&fTreeCascVarNegdEdx,"fTreeCascVarNegdEdx/F");
            fTreeCascade->Branch("fTreeCascVarBachdEdx",&fTreeCascVarBachdEdx,"fTreeCascVarBachdEdx/F");
            fTreeCascade->Branch("fTreeCascVarPosInnerP",&fTreeCascVarPosInnerP,"fTreeCascVarPosInnerP/F");
            fTreeCascade->Branch("fTreeCascVarNegInnerP",&fTreeCascVarNegInnerP,"fTreeCascVarNegInnerP/F");
            fTreeCascade->Branch("fTreeCascVarBachInnerP",&fTreeCascVarBachInnerP,"fTreeCascVarBachInnerP/F");
            fTreeCascade->Branch("fTreeCascVarNegTrackStatus",&fTreeCascVarNegTrackStatus,"fTreeCascVarNegTrackStatus/l");
            fTreeCascade->Branch("fTreeCascVarPosTrackStatus",&fTreeCascVarPosTrackStatus,"fTreeCascVarPosTrackStatus/l");
            fTreeCascade->Branch("fTreeCascVarBachTrackStatus",&fTreeCascVarBachTrackStatus,"fTreeCascVarBachTrackStatus/l");
            fTreeCascade->Branch("fTreeCascVarNegDCAz",&fTreeCascVarNegDCAz,"fTreeCascVarNegDCAz/F");
            fTreeCascade->Branch("fTreeCascVarPosDCAz",&fTreeCascVarPosDCAz,"fTreeCascVarPosDCAz/F");
            fTreeCascade->Branch("fTreeCascVarBachDCAz",&fTreeCascVarBachDCAz,"fTreeCascVarBachDCAz/F");
        }
        //------------------------------------------------
        if ( fkDebugBump ){
            //Variables for debugging the invariant mass bump
            //Full momentum information
            fTreeCascade->Branch("fTreeCascVarPosPx",&fTreeCascVarPosPx,"fTreeCascVarPosPx/F");
            fTreeCascade->Branch("fTreeCascVarPosPy",&fTreeCascVarPosPy,"fTreeCascVarPosPy/F");
            fTreeCascade->Branch("fTreeCascVarPosPz",&fTreeCascVarPosPz,"fTreeCascVarPosPz/F");
            fTreeCascade->Branch("fTreeCascVarNegPx",&fTreeCascVarNegPx,"fTreeCascVarNegPx/F");
            fTreeCascade->Branch("fTreeCascVarNegPy",&fTreeCascVarNegPy,"fTreeCascVarNegPy/F");
            fTreeCascade->Branch("fTreeCascVarNegPz",&fTreeCascVarNegPz,"fTreeCascVarNegPz/F");
            fTreeCascade->Branch("fTreeCascVarBachPx",&fTreeCascVarBachPx,"fTreeCascVarBachPx/F");
            fTreeCascade->Branch("fTreeCascVarBachPy",&fTreeCascVarBachPy,"fTreeCascVarBachPy/F");
            fTreeCascade->Branch("fTreeCascVarBachPz",&fTreeCascVarBachPz,"fTreeCascVarBachPz/F");
            // Decay positions
            fTreeCascade->Branch("fTreeCascVarV0DecayX",&fTreeCascVarV0DecayX,"fTreeCascVarV0DecayX/F");
            fTreeCascade->Branch("fTreeCascVarV0DecayY",&fTreeCascVarV0DecayY,"fTreeCascVarV0DecayY/F");
            fTreeCascade->Branch("fTreeCascVarV0DecayZ",&fTreeCascVarV0DecayZ,"fTreeCascVarV0DecayZ/F");
            fTreeCascade->Branch("fTreeCascVarCascadeDecayX",&fTreeCascVarCascadeDecayX,"fTreeCascVarCascadeDecayX/F");
            fTreeCascade->Branch("fTreeCascVarCascadeDecayY",&fTreeCascVarCascadeDecayY,"fTreeCascVarCascadeDecayY/F");
            fTreeCascade->Branch("fTreeCascVarCascadeDecayZ",&fTreeCascVarCascadeDecayZ,"fTreeCascVarCascadeDecayZ/F");
            
            fTreeCascade->Branch("fTreeCascVarV0Lifetime",&fTreeCascVarV0Lifetime,"fTreeCascVarV0Lifetime/F");
            
            //Track Labels (check for duplicates, etc)
            fTreeCascade->Branch("fTreeCascVarNegIndex",&fTreeCascVarNegIndex,"fTreeCascVarNegIndex/I");
            fTreeCascade->Branch("fTreeCascVarPosIndex",&fTreeCascVarPosIndex,"fTreeCascVarPosIndex/I");
            fTreeCascade->Branch("fTreeCascVarBachIndex",&fTreeCascVarBachIndex,"fTreeCascVarBachIndex/I");
            //Event Number (check same-event index mixups)
            fTreeCascade->Branch("fTreeCascVarEventNumber",&fTreeCascVarEventNumber,"fTreeCascVarEventNumber/l");
        }
        if ( fkDebugOOBPileup ) {
            fTreeCascade->Branch("fTreeCascVarPosITSClusters0",&fTreeCascVarPosITSClusters0,"fTreeCascVarPosITSClusters0/O");
            fTreeCascade->Branch("fTreeCascVarPosITSClusters1",&fTreeCascVarPosITSClusters1,"fTreeCascVarPosITSClusters1/O");
            fTreeCascade->Branch("fTreeCascVarPosITSClusters2",&fTreeCascVarPosITSClusters2,"fTreeCascVarPosITSClusters2/O");
            fTreeCascade->Branch("fTreeCascVarPosITSClusters3",&fTreeCascVarPosITSClusters3,"fTreeCascVarPosITSClusters3/O");
            fTreeCascade->Branch("fTreeCascVarPosITSClusters4",&fTreeCascVarPosITSClusters4,"fTreeCascVarPosITSClusters4/O");
            fTreeCascade->Branch("fTreeCascVarPosITSClusters5",&fTreeCascVarPosITSClusters5,"fTreeCascVarPosITSClusters5/O");
            
            fTreeCascade->Branch("fTreeCascVarNegITSClusters0",&fTreeCascVarNegITSClusters0,"fTreeCascVarNegITSClusters0/O");
            fTreeCascade->Branch("fTreeCascVarNegITSClusters1",&fTreeCascVarNegITSClusters1,"fTreeCascVarNegITSClusters1/O");
            fTreeCascade->Branch("fTreeCascVarNegITSClusters2",&fTreeCascVarNegITSClusters2,"fTreeCascVarNegITSClusters2/O");
            fTreeCascade->Branch("fTreeCascVarNegITSClusters3",&fTreeCascVarNegITSClusters3,"fTreeCascVarNegITSClusters3/O");
            fTreeCascade->Branch("fTreeCascVarNegITSClusters4",&fTreeCascVarNegITSClusters4,"fTreeCascVarNegITSClusters4/O");
            fTreeCascade->Branch("fTreeCascVarNegITSClusters5",&fTreeCascVarNegITSClusters5,"fTreeCascVarNegITSClusters5/O");
            
            fTreeCascade->Branch("fTreeCascVarBachITSClusters0",&fTreeCascVarBachITSClusters0,"fTreeCascVarBachITSClusters0/O");
            fTreeCascade->Branch("fTreeCascVarBachITSClusters1",&fTreeCascVarBachITSClusters1,"fTreeCascVarBachITSClusters1/O");
            fTreeCascade->Branch("fTreeCascVarBachITSClusters2",&fTreeCascVarBachITSClusters2,"fTreeCascVarBachITSClusters2/O");
            fTreeCascade->Branch("fTreeCascVarBachITSClusters3",&fTreeCascVarBachITSClusters3,"fTreeCascVarBachITSClusters3/O");
            fTreeCascade->Branch("fTreeCascVarBachITSClusters4",&fTreeCascVarBachITSClusters4,"fTreeCascVarBachITSClusters4/O");
            fTreeCascade->Branch("fTreeCascVarBachITSClusters5",&fTreeCascVarBachITSClusters5,"fTreeCascVarBachITSClusters5/O");
            
            fTreeCascade->Branch("fTreeCascVarPosITSSharedClusters0",&fTreeCascVarPosITSSharedClusters0,"fTreeCascVarPosITSSharedClusters0/O");
            fTreeCascade->Branch("fTreeCascVarPosITSSharedClusters1",&fTreeCascVarPosITSSharedClusters1,"fTreeCascVarPosITSSharedClusters1/O");
            fTreeCascade->Branch("fTreeCascVarPosITSSharedClusters2",&fTreeCascVarPosITSSharedClusters2,"fTreeCascVarPosITSSharedClusters2/O");
            fTreeCascade->Branch("fTreeCascVarPosITSSharedClusters3",&fTreeCascVarPosITSSharedClusters3,"fTreeCascVarPosITSSharedClusters3/O");
            fTreeCascade->Branch("fTreeCascVarPosITSSharedClusters4",&fTreeCascVarPosITSSharedClusters4,"fTreeCascVarPosITSSharedClusters4/O");
            fTreeCascade->Branch("fTreeCascVarPosITSSharedClusters5",&fTreeCascVarPosITSSharedClusters5,"fTreeCascVarPosITSSharedClusters5/O");
            
            fTreeCascade->Branch("fTreeCascVarNegITSSharedClusters0",&fTreeCascVarNegITSSharedClusters0,"fTreeCascVarNegITSSharedClusters0/O");
            fTreeCascade->Branch("fTreeCascVarNegITSSharedClusters1",&fTreeCascVarNegITSSharedClusters1,"fTreeCascVarNegITSSharedClusters1/O");
            fTreeCascade->Branch("fTreeCascVarNegITSSharedClusters2",&fTreeCascVarNegITSSharedClusters2,"fTreeCascVarNegITSSharedClusters2/O");
            fTreeCascade->Branch("fTreeCascVarNegITSSharedClusters3",&fTreeCascVarNegITSSharedClusters3,"fTreeCascVarNegITSSharedClusters3/O");
            fTreeCascade->Branch("fTreeCascVarNegITSSharedClusters4",&fTreeCascVarNegITSSharedClusters4,"fTreeCascVarNegITSSharedClusters4/O");
            fTreeCascade->Branch("fTreeCascVarNegITSSharedClusters5",&fTreeCascVarNegITSSharedClusters5,"fTreeCascVarNegITSSharedClusters5/O");
            
            fTreeCascade->Branch("fTreeCascVarBachITSSharedClusters0",&fTreeCascVarBachITSSharedClusters0,"fTreeCascVarBachITSSharedClusters0/O");
            fTreeCascade->Branch("fTreeCascVarBachITSSharedClusters1",&fTreeCascVarBachITSSharedClusters1,"fTreeCascVarBachITSSharedClusters1/O");
            fTreeCascade->Branch("fTreeCascVarBachITSSharedClusters2",&fTreeCascVarBachITSSharedClusters2,"fTreeCascVarBachITSSharedClusters2/O");
            fTreeCascade->Branch("fTreeCascVarBachITSSharedClusters3",&fTreeCascVarBachITSSharedClusters3,"fTreeCascVarBachITSSharedClusters3/O");
            fTreeCascade->Branch("fTreeCascVarBachITSSharedClusters4",&fTreeCascVarBachITSSharedClusters4,"fTreeCascVarBachITSSharedClusters4/O");
            fTreeCascade->Branch("fTreeCascVarBachITSSharedClusters5",&fTreeCascVarBachITSSharedClusters5,"fTreeCascVarBachITSSharedClusters5/O");
            fTreeCascade->Branch("fTreeCascVarNegTOFExpTDiff",&fTreeCascVarNegTOFExpTDiff,"fTreeCascVarNegTOFExpTDiff/F");
            fTreeCascade->Branch("fTreeCascVarPosTOFExpTDiff",&fTreeCascVarPosTOFExpTDiff,"fTreeCascVarPosTOFExpTDiff/F");
            fTreeCascade->Branch("fTreeCascVarBachTOFExpTDiff",&fTreeCascVarBachTOFExpTDiff,"fTreeCascVarBachTOFExpTDiff/F");
            fTreeCascade->Branch("fTreeCascVarNegTOFSignal",&fTreeCascVarNegTOFSignal,"fTreeCascVarNegTOFSignal/F");
            fTreeCascade->Branch("fTreeCascVarPosTOFSignal",&fTreeCascVarPosTOFSignal,"fTreeCascVarPosTOFSignal/F");
            fTreeCascade->Branch("fTreeCascVarBachTOFSignal",&fTreeCascVarBachTOFSignal,"fTreeCascVarBachTOFSignal/F");
            fTreeCascade->Branch("fTreeCascVarNegTOFBCid",&fTreeCascVarNegTOFBCid,"fTreeCascVarNegTOFBCid/I");
            fTreeCascade->Branch("fTreeCascVarPosTOFBCid",&fTreeCascVarPosTOFBCid,"fTreeCascVarPosTOFBCid/I");
            fTreeCascade->Branch("fTreeCascVarBachTOFBCid",&fTreeCascVarBachTOFBCid,"fTreeCascVarBachTOFBCid/I");
            // Event info
            fTreeCascade->Branch("fTreeCascVarOOBPileupFlag",&fTreeCascVarOOBPileupFlag,"fTreeCascVarOOBPileupFlag/O");
            fTreeCascade->Branch("fTreeCascVarAmplitudeV0A",&fTreeCascVarAmplitudeV0A,"fTreeCascVarAmplitudeV0A/F");
            fTreeCascade->Branch("fTreeCascVarAmplitudeV0C",&fTreeCascVarAmplitudeV0C,"fTreeCascVarAmplitudeV0C/F");
            fTreeCascade->Branch("fTreeCascVarClosestNonEmptyBC",&fTreeCascVarClosestNonEmptyBC,"fTreeCascVarClosestNonEmptyBC/I");
        }
        //Kink tagging
        fTreeCascade->Branch("fTreeCascVarBachIsKink",&fTreeCascVarBachIsKink,"fTreeCascVarBachIsKink/O");
        fTreeCascade->Branch("fTreeCascVarPosIsKink",&fTreeCascVarPosIsKink,"fTreeCascVarPosIsKink/O");
        fTreeCascade->Branch("fTreeCascVarNegIsKink",&fTreeCascVarNegIsKink,"fTreeCascVarNegIsKink/O");
        //------------------------------------------------
    }
    //------------------------------------------------
    // Particle Identification Setup
    //------------------------------------------------
    
    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    inputHandler->SetNeedField();
    
    // Multiplicity
    if(! fESDtrackCuts ) {
        fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,kFALSE);
        fESDtrackCuts->SetPtRange(0.15);  // adding pt cut
        fESDtrackCuts->SetEtaRange(-1.0, 1.0);
    }
    //Analysis Utils
    if(! fUtils ) {
        fUtils = new AliAnalysisUtils();
    }
    if(! fRand ){
        fRand = new TRandom3();
        // From TRandom3 reference:
        // if seed is 0 (default value) a TUUID is generated and
        // used to fill the first 8 integers of the seed array
        fRand->SetSeed(0);
    }
    
    // OOB Pileup in pp 2016
    if( !fESDtrackCutsGlobal2015 && fkDebugOOBPileup ) {
        fESDtrackCutsGlobal2015 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2015PbPb(kTRUE,kFALSE);
        //Initial set of cuts - to be adjusted
        fESDtrackCutsGlobal2015->SetPtRange(0.15);
        fESDtrackCutsGlobal2015->SetEtaRange(-1.0, 1.0);
    }
    if( !fESDtrackCutsITSsa2010 && fkDebugOOBPileup ) {
        fESDtrackCutsITSsa2010 = AliESDtrackCuts::GetStandardITSSATrackCuts2010();
    }
    
    //------------------------------------------------
    // V0 Multiplicity Histograms
    //------------------------------------------------
    
    // Create histograms
    fListHist = new TList();
    fListHist->SetOwner();  // See http://root.cern.ch/root/html/TCollection.html#TCollection:SetOwner
    
    fEventCuts.AddQAplotsToList(fListHist);
    
    if(! fHistEventCounter ) {
        //Histogram Output: Event-by-Event
        fHistEventCounter = new TH1D( "fHistEventCounter", ";Evt. Sel. Step;Count",2,0,2);
        fHistEventCounter->GetXaxis()->SetBinLabel(1, "Processed");
        fHistEventCounter->GetXaxis()->SetBinLabel(2, "Selected");
        fListHist->Add(fHistEventCounter);
    }
    
    if(! fHistCentrality_V0A ) {
        //Histogram Output: Event-by-Event
        fHistCentrality_V0A = new TH1D( "fHistCentrality_V0A", "WARNING: no pileup rejection applied!;Centrality;Event Count",100,0,100);
        fListHist->Add(fHistCentrality_V0A);
    }
    if(! fHistCentrality_V0C ) {
        //Histogram Output: Event-by-Event
        fHistCentrality_V0C = new TH1D( "fHistCentrality_V0C", "WARNING: no pileup rejection applied!;Centrality;Event Count",100,0,100);
        fListHist->Add(fHistCentrality_V0C);
    } if(! fHistCentrality_V0M ) {
        //Histogram Output: Event-by-Event
        fHistCentrality_V0M = new TH1D( "fHistCentrality_V0M", "WARNING: no pileup rejection applied!;Centrality;Event Count",100,0,100);
        fListHist->Add(fHistCentrality_V0M);
    }
    
    //Superlight mode output
    if ( !fListV0 ){
        fListV0 = new TList();
        fListV0->SetOwner();
    }
    
    if ( !fListCascade ){
        //Superlight mode output
        fListCascade = new TList();
        fListCascade->SetOwner();
    }
    
    //Regular Output: Slots 1, 2, 3
    PostData(1, fListHist    );
    PostData(2, fListV0      );
    PostData(3, fListCascade );
    
    //TTree Objects: Slots 4, 5, 6
    if(fkSaveEventTree)    PostData(4, fTreeEvent   );
    if(fkSaveV0Tree)       PostData(5, fTreeV0      );
    if(fkSaveCascadeTree)  PostData(6, fTreeCascade );
    
}// end UserCreateOutputObjects


//________________________________________________________________________
void AliAnalysisTaskStrangenessVsMultiplicityRun2pPb::UserExec(Option_t *)
{
    // Main loop
    // Called for each event
    AliESDEvent *lESDevent = 0x0;
    
    // Connect to the InputEvent
    // After these lines, we should have an ESD/AOD event + the number of V0s in it.
    
    // Appropriate for ESD analysis!
    
    lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
    if (!lESDevent) {
        AliWarning("ERROR: lESDevent not available \n");
        return;
    }
    
    //Get VZERO Information for multiplicity later
    AliVVZERO* esdV0 = lESDevent->GetVZEROData();
    if (!esdV0) {
        AliError("AliVVZERO not available");
        return;
    }
    
    Double_t lMagneticField = -10;
    lMagneticField = lESDevent->GetMagneticField( );
    
    //------------------------------------------------
    // Retrieving IR info for OOB Pileup rejection
    //------------------------------------------------
    if( fkDebugOOBPileup ) {
        fClosestNonEmptyBC = 10*3564; // start with an isolated event
        AliESDHeader* lESDHeader = (AliESDHeader*)lESDevent->GetHeader();
        Int_t    nIRs       = lESDHeader->GetTriggerIREntries();
        Long64_t lThisOrbit = lESDHeader->GetOrbitNumber();
        Int_t    lThisBC    = lESDHeader->GetBunchCrossNumber();
        
        const AliTriggerIR* lIR;
        for(Int_t i=0; i<nIRs; i++) {
            
            lIR = lESDHeader->GetTriggerIR(i);
            Long64_t lOrbit     = lIR->GetOrbit();
            UInt_t   lNWord     = lIR->GetNWord();
            UShort_t *lBCsForIR = lIR->GetBCs();
            Bool_t   *lInt1     = lIR->GetInt1s();
            Bool_t   *lInt2     = lIR->GetInt2s();
            
            for(UInt_t j=0; j<lNWord; j++) {
                
                if( (lInt1[j]) || (lInt2[j]) ) {
                    
                    Int_t lBC = lBCsForIR[j];
                    
                    if((lOrbit == lThisOrbit) && (lBC != lThisBC)) {
                        Int_t lClosestNonEmptyBC = lBC - lThisBC;
                        if(TMath::Abs(lClosestNonEmptyBC)<TMath::Abs(fClosestNonEmptyBC)) fClosestNonEmptyBC = lClosestNonEmptyBC;
                    }
                    
                    if(lOrbit == (lThisOrbit+1)) {
                        Int_t lClosestNonEmptyBC = (lBC+3564) - lThisBC;
                        if(TMath::Abs(lClosestNonEmptyBC)<TMath::Abs(fClosestNonEmptyBC)) fClosestNonEmptyBC = lClosestNonEmptyBC;
                    }
                    
                    if(lOrbit == (lThisOrbit-1)) {
                        Int_t lClosestNonEmptyBC = (lBC-3564) - lThisBC;
                        if(TMath::Abs(lClosestNonEmptyBC)<TMath::Abs(fClosestNonEmptyBC)) fClosestNonEmptyBC = lClosestNonEmptyBC;
                    }
                }
            }
        }
    }
    // done with IR ----------------------------------
    
    //------------------------------------------------
    // Event Selection ---
    //  --- Performed entirely via AliPPVsMultUtils
    // (except removal of incomplete events and SPDClusterVsTracklets cut)
    //------------------------------------------------
    
    //Copy-paste of steps done in AliAnalysisTaskSkeleton
    
    fHistEventCounter->Fill(0.5);
    
    //------------------------------------------------
    // Primary Vertex Requirements Section:
    //  ---> pp: has vertex, |z|<10cm
    //------------------------------------------------
    
    //classical Proton-proton like selection
    const AliESDVertex *lPrimaryBestESDVtx     = lESDevent->GetPrimaryVertex();
    const AliESDVertex *lPrimaryTrackingESDVtx = lESDevent->GetPrimaryVertexTracks();
    const AliESDVertex *lPrimarySPDVtx         = lESDevent->GetPrimaryVertexSPD();
    
    Double_t lBestPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
    lPrimaryBestESDVtx->GetXYZ( lBestPrimaryVtxPos );
    
    //------------------------------------------------
    // Multiplicity Information Acquistion
    //------------------------------------------------
    
    Float_t lPercentile_V0A = 500;
    Float_t lPercentile_V0C = 500;
    Float_t lPercentile_V0M = 500;
    
    Int_t lEvSelCode = 100;
    AliMultSelection *MultSelection = (AliMultSelection*) lESDevent -> FindListObject("MultSelection");
    if( !MultSelection) {
        //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
        AliWarning("AliMultSelection object not found!");
    } else {
        //V0C Multiplicity Percentile
        lPercentile_V0A = MultSelection->GetMultiplicityPercentile("V0A");
        lPercentile_V0C = MultSelection->GetMultiplicityPercentile("V0C");
        lPercentile_V0M = MultSelection->GetMultiplicityPercentile("V0M");
        //Event Selection Code
        lEvSelCode = MultSelection->GetEvSelCode();
    }
    
    //just ask AliMultSelection. It will know.
    fMVPileupFlag = kFALSE;
    fMVPileupFlag = MultSelection->GetThisEventIsNotPileupMV();
    
    fCentrality_V0A = lPercentile_V0A;
    fCentrality_V0C = lPercentile_V0C;
    fCentrality_V0M = lPercentile_V0M;
    
    
    if( lEvSelCode != 0 ) {
        PostData(1, fListHist    );
        PostData(2, fListV0      );
        PostData(3, fListCascade );
        if( fkSaveEventTree   ) PostData(4, fTreeEvent   );
        if( fkSaveV0Tree      ) PostData(5, fTreeV0      );
        if( fkSaveCascadeTree ) PostData(6, fTreeCascade );
        return;
    }
    
    AliVEvent *ev = InputEvent();
    if( fkDoExtraEvSels ) {
        if( !fEventCuts.AcceptEvent(ev) ) {
            PostData(1, fListHist    );
            PostData(2, fListV0      );
            PostData(3, fListCascade );
            if( fkSaveEventTree   ) PostData(4, fTreeEvent   );
            if( fkSaveV0Tree      ) PostData(5, fTreeV0      );
            if( fkSaveCascadeTree ) PostData(6, fTreeCascade );
            return;
        }
    }
    
    fHistEventCounter->Fill(1.5);
    
    //Bookkeep event number for debugging
    fTreeCascVarEventNumber =
    ( ( ((ULong64_t)lESDevent->GetPeriodNumber() ) << 36 ) |
     ( ((ULong64_t)lESDevent->GetOrbitNumber () ) << 12 ) |
     ((ULong64_t)lESDevent->GetBunchCrossNumber() )  );
    
    //Save info for pileup study (high multiplicity triggers pp 13 TeV - 2016 data)
    if( fkDebugOOBPileup ) {
        fOOBPileupFlag     = !fUtils->IsOutOfBunchPileUp(ev);
        fNTOFClusters      = lESDevent->GetESDTOFClusters()->GetEntriesFast();
        fNTOFMatches       = lESDevent->GetESDTOFMatches()->GetEntriesFast();
        fNTracksITSsa2010  = 0;
        fNTracksGlobal2015 = 0;
        fNTracksGlobal2015TriggerPP = 0;
        //Count tracks with various selections
        for(Long_t itrack = 0; itrack<lESDevent->GetNumberOfTracks(); itrack++) {
            AliVTrack *track = lESDevent -> GetVTrack( itrack );
            if( !track ) continue;
            //Only ITSsa tracks
            if( fESDtrackCutsITSsa2010->AcceptVTrack(track) ) fNTracksITSsa2010++;
            if( !fESDtrackCutsGlobal2015->AcceptVTrack(track) ) continue;
            //Only for accepted tracks
            fNTracksGlobal2015++;
            //Count accepted + TOF time window (info from Alberica)
            //Warning: 12.5 is appropriate for pp (for Pb-Pb use 30)
            if( TMath::Abs( track->GetTOFExpTDiff() ) < 12.5 ) fNTracksGlobal2015TriggerPP++;
        }
        
        //VZERO info
        fAmplitudeV0A = ((AliMultEstimator*)MultSelection->GetEstimator("V0A"))->GetValue();
        fAmplitudeV0C = ((AliMultEstimator*)MultSelection->GetEstimator("V0C"))->GetValue();
        
        
    }
    
    //Fill centrality histogram
    fHistCentrality_V0A->Fill(fCentrality_V0A);
    fHistCentrality_V0C->Fill(fCentrality_V0C);
    fHistCentrality_V0M->Fill(fCentrality_V0M);
    
    
    //Event-level fill
    if ( fkSaveEventTree ) fTreeEvent->Fill() ;
    
    //STOP HERE if skipping event selections (no point in doing the rest...)
    
    //------------------------------------------------
    
    //------------------------------------------------
    // Fill V0 Tree as needed
    //------------------------------------------------
    
    //Variable definition
    Int_t    lOnFlyStatus = 0;// nv0sOn = 0, nv0sOff = 0;
    Double_t lChi2V0 = 0;
    Double_t lDcaV0Daughters = 0, lDcaV0ToPrimVertex = 0;
    Double_t lDcaPosToPrimVertex = 0, lDcaNegToPrimVertex = 0;
    Double_t lV0CosineOfPointingAngle = 0;
    Double_t lV0Radius = 0, lPt = 0;
    Double_t lRapK0Short = 0, lRapLambda = 0;
    Double_t lInvMassK0s = 0, lInvMassLambda = 0, lInvMassAntiLambda = 0;
    Double_t lAlphaV0 = 0, lPtArmV0 = 0;
    
    Double_t fMinV0Pt = 0;
    Double_t fMaxV0Pt = 100;
    
    //------------------------------------------------
    // Rerun V0 Vertexer!
    // WARNING: this will only work if the
    // special "use on the fly cascading" flag
    // is disabled!
    //------------------------------------------------
    
    if( fkRunVertexers && !fkUseOnTheFlyV0Cascading ) {
        //Only reset if not using on-the-fly (or else nothing passes)
        lESDevent->ResetV0s();
        
        //Decide between regular and light vertexer (default: light)
        if ( ! fkUseLightVertexer ){
            //Instantiate vertexer object
            AliV0vertexer lV0vtxer;
            //Set Cuts
            lV0vtxer.SetDefaultCuts(fV0VertexerSels);
            lV0vtxer.SetCuts(fV0VertexerSels);
            //Redo
            lV0vtxer.Tracks2V0vertices(lESDevent);
        } else {
            //Instantiate vertexer object
            AliLightV0vertexer lV0vtxer;
            //Set do or don't do V0 refit for improved precision
            lV0vtxer.SetDoRefit( kFALSE );
            if (fkDoV0Refit) lV0vtxer.SetDoRefit(kTRUE);
            //Set Cuts
            lV0vtxer.SetDefaultCuts(fV0VertexerSels);
            lV0vtxer.SetCuts(fV0VertexerSels);
            //Redo
            lV0vtxer.Tracks2V0vertices(lESDevent);
        }
    }
    
    Int_t nv0s = 0;
    nv0s = lESDevent->GetNumberOfV0s();
    
    for (Int_t iV0 = 0; iV0 < nv0s; iV0++) //extra-crazy test
    {   // This is the begining of the V0 loop
        AliESDv0 *v0 = ((AliESDEvent*)lESDevent)->GetV0(iV0);
        if (!v0) continue;
        
        CheckChargeV0( v0 );
        //Remove like-sign (will not affect offline V0 candidates!)
        if( v0->GetParamN()->Charge() > 0 && v0->GetParamP()->Charge() > 0 ){
            continue;
        }
        if( v0->GetParamN()->Charge() < 0 && v0->GetParamP()->Charge() < 0 ){
            continue;
        }
        
        Double_t tDecayVertexV0[3];
        v0->GetXYZ(tDecayVertexV0[0],tDecayVertexV0[1],tDecayVertexV0[2]);
        
        Double_t tV0mom[3];
        v0->GetPxPyPz( tV0mom[0],tV0mom[1],tV0mom[2] );
        Double_t lV0TotalMomentum = TMath::Sqrt(
                                                tV0mom[0]*tV0mom[0]+tV0mom[1]*tV0mom[1]+tV0mom[2]*tV0mom[2] );
        
        lV0Radius = TMath::Sqrt(tDecayVertexV0[0]*tDecayVertexV0[0]+tDecayVertexV0[1]*tDecayVertexV0[1]);
        
        lPt = v0->Pt();
        lRapK0Short = v0->RapK0Short();
        lRapLambda  = v0->RapLambda();
        if ((lPt<fMinV0Pt)||(fMaxV0Pt<lPt)) continue;
        
        UInt_t lKeyPos = (UInt_t)TMath::Abs(v0->GetPindex());
        UInt_t lKeyNeg = (UInt_t)TMath::Abs(v0->GetNindex());
        
        Double_t lMomPos[3];
        v0->GetPPxPyPz(lMomPos[0],lMomPos[1],lMomPos[2]);
        Double_t lMomNeg[3];
        v0->GetNPxPyPz(lMomNeg[0],lMomNeg[1],lMomNeg[2]);
        
        AliESDtrack *pTrack=((AliESDEvent*)lESDevent)->GetTrack(lKeyPos);
        AliESDtrack *nTrack=((AliESDEvent*)lESDevent)->GetTrack(lKeyNeg);
        
        if (!pTrack || !nTrack) {
            Printf("ERROR: Could not retreive one of the daughter track");
            continue;
        }
        fTreeVariablePosPIDForTracking = pTrack->GetPIDForTracking();
        fTreeVariableNegPIDForTracking = nTrack->GetPIDForTracking();
        
        //Check its clusters
        fTreeVariablePosITSClusters0 = pTrack->HasPointOnITSLayer(0);
        fTreeVariablePosITSClusters1 = pTrack->HasPointOnITSLayer(1);
        fTreeVariablePosITSClusters2 = pTrack->HasPointOnITSLayer(2);
        fTreeVariablePosITSClusters3 = pTrack->HasPointOnITSLayer(3);
        fTreeVariablePosITSClusters4 = pTrack->HasPointOnITSLayer(4);
        fTreeVariablePosITSClusters5 = pTrack->HasPointOnITSLayer(5);
        
        fTreeVariableNegITSClusters0 = nTrack->HasPointOnITSLayer(0);
        fTreeVariableNegITSClusters1 = nTrack->HasPointOnITSLayer(1);
        fTreeVariableNegITSClusters2 = nTrack->HasPointOnITSLayer(2);
        fTreeVariableNegITSClusters3 = nTrack->HasPointOnITSLayer(3);
        fTreeVariableNegITSClusters4 = nTrack->HasPointOnITSLayer(4);
        fTreeVariableNegITSClusters5 = nTrack->HasPointOnITSLayer(5);
        
        //Check its clusters, shared
        fTreeVariablePosITSSharedClusters0 = pTrack->HasSharedPointOnITSLayer(0);
        fTreeVariablePosITSSharedClusters1 = pTrack->HasSharedPointOnITSLayer(1);
        fTreeVariablePosITSSharedClusters2 = pTrack->HasSharedPointOnITSLayer(2);
        fTreeVariablePosITSSharedClusters3 = pTrack->HasSharedPointOnITSLayer(3);
        fTreeVariablePosITSSharedClusters4 = pTrack->HasSharedPointOnITSLayer(4);
        fTreeVariablePosITSSharedClusters5 = pTrack->HasSharedPointOnITSLayer(5);
        
        fTreeVariableNegITSSharedClusters0 = nTrack->HasSharedPointOnITSLayer(0);
        fTreeVariableNegITSSharedClusters1 = nTrack->HasSharedPointOnITSLayer(1);
        fTreeVariableNegITSSharedClusters2 = nTrack->HasSharedPointOnITSLayer(2);
        fTreeVariableNegITSSharedClusters3 = nTrack->HasSharedPointOnITSLayer(3);
        fTreeVariableNegITSSharedClusters4 = nTrack->HasSharedPointOnITSLayer(4);
        fTreeVariableNegITSSharedClusters5 = nTrack->HasSharedPointOnITSLayer(5);
        
        const AliExternalTrackParam *innernegv0=nTrack->GetInnerParam();
        const AliExternalTrackParam *innerposv0=pTrack->GetInnerParam();
        Float_t lThisPosInnerP = -1;
        Float_t lThisNegInnerP = -1;
        Float_t lThisPosInnerPt = -1;
        Float_t lThisNegInnerPt = -1;
        if(innerposv0)  { lThisPosInnerP  = innerposv0 ->GetP(); }
        if(innernegv0)  { lThisNegInnerP  = innernegv0 ->GetP(); }
        if(innerposv0)  { lThisPosInnerPt  = innerposv0 ->Pt(); }
        if(innernegv0)  { lThisNegInnerPt  = innernegv0 ->Pt(); }
        Float_t lThisPosdEdx = pTrack -> GetTPCsignal();
        Float_t lThisNegdEdx = nTrack -> GetTPCsignal();
        
        fTreeVariablePosdEdx = lThisPosdEdx;
        fTreeVariableNegdEdx = lThisNegdEdx;
        
        fTreeVariablePosInnerP = lThisPosInnerP;
        fTreeVariableNegInnerP = lThisNegInnerP;
        
        //Daughter Eta for Eta selection, afterwards
        fTreeVariableNegEta = nTrack->Eta();
        fTreeVariablePosEta = pTrack->Eta();
        
        if ( fkExtraCleanup ){
            if( TMath::Abs(fTreeVariableNegEta)>0.8 || TMath::Abs(fTreeVariablePosEta)>0.8 ) continue;
            if( TMath::Abs(lRapK0Short        )>0.5 && TMath::Abs(lRapLambda         )>0.5 ) continue;
        }
        
        // Filter like-sign V0 (next: add counter and distribution)
        if ( pTrack->GetSign() == nTrack->GetSign()) {
            continue;
        }
        
        //________________________________________________________________________
        // Track quality cuts
        Float_t lPosTrackCrossedRows = pTrack->GetTPCClusterInfo(2,1);
        Float_t lNegTrackCrossedRows = nTrack->GetTPCClusterInfo(2,1);
        fTreeVariableLeastNbrCrossedRows = (Int_t) lPosTrackCrossedRows;
        if( lNegTrackCrossedRows < fTreeVariableLeastNbrCrossedRows )
            fTreeVariableLeastNbrCrossedRows = (Int_t) lNegTrackCrossedRows;
        
        // TPC refit condition (done during reconstruction for Offline but not for On-the-fly)
        if( !(pTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
        if( !(nTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
        
        //Get status flags
        fTreeVariablePosTrackStatus = pTrack->GetStatus();
        fTreeVariableNegTrackStatus = nTrack->GetStatus();
        
        fTreeVariablePosDCAz = GetDCAz(pTrack);
        fTreeVariableNegDCAz = GetDCAz(nTrack);
        
        //GetKinkIndex condition
        if( pTrack->GetKinkIndex(0)>0 || nTrack->GetKinkIndex(0)>0 ) continue;
        
        //Findable clusters > 0 condition
        if( pTrack->GetTPCNclsF()<=0 || nTrack->GetTPCNclsF()<=0 ) continue;
        
        //Compute ratio Crossed Rows / Findable clusters
        //Note: above test avoids division by zero!
        Float_t lPosTrackCrossedRowsOverFindable = lPosTrackCrossedRows / ((double)(pTrack->GetTPCNclsF()));
        Float_t lNegTrackCrossedRowsOverFindable = lNegTrackCrossedRows / ((double)(nTrack->GetTPCNclsF()));
        
        fTreeVariableLeastRatioCrossedRowsOverFindable = lPosTrackCrossedRowsOverFindable;
        if( lNegTrackCrossedRowsOverFindable < fTreeVariableLeastRatioCrossedRowsOverFindable )
            fTreeVariableLeastRatioCrossedRowsOverFindable = lNegTrackCrossedRowsOverFindable;
        
        //Lowest Cut Level for Ratio Crossed Rows / Findable = 0.8, set here
        if ( fTreeVariableLeastRatioCrossedRowsOverFindable < 0.8 ) continue;
        
        //Extra track quality: Chi2/cluster for cross-checks
        Float_t lBiggestChi2PerCluster = -1;
        
        Float_t lPosChi2PerCluster = 1000;
        Float_t lNegChi2PerCluster = 1000;
        
        if( pTrack->GetTPCNcls() > 0 ) lPosChi2PerCluster = pTrack->GetTPCchi2() / ((Float_t)pTrack->GetTPCNcls());
        if( nTrack->GetTPCNcls() > 0 ) lNegChi2PerCluster = nTrack->GetTPCchi2() / ((Float_t)nTrack->GetTPCNcls());
        
        if ( lPosChi2PerCluster  > lBiggestChi2PerCluster ) lBiggestChi2PerCluster = lPosChi2PerCluster;
        if ( lNegChi2PerCluster  > lBiggestChi2PerCluster ) lBiggestChi2PerCluster = lNegChi2PerCluster;
        
        fTreeVariableMaxChi2PerCluster = lBiggestChi2PerCluster;
        
        //Extra track quality: min track length
        Float_t lSmallestTrackLength = 1000;
        Float_t lPosTrackLength = -1;
        Float_t lNegTrackLength = -1;
        
        if (pTrack->GetInnerParam()) lPosTrackLength = pTrack->GetLengthInActiveZone(1, 2.0, 220.0, lESDevent->GetMagneticField());
        if (nTrack->GetInnerParam()) lNegTrackLength = nTrack->GetLengthInActiveZone(1, 2.0, 220.0, lESDevent->GetMagneticField());
        
        if ( lPosTrackLength  < lSmallestTrackLength ) lSmallestTrackLength = lPosTrackLength;
        if ( lNegTrackLength  < lSmallestTrackLength ) lSmallestTrackLength = lNegTrackLength;
        
        fTreeVariableMinTrackLength = lSmallestTrackLength;
        
        if ( ( ( ( pTrack->GetTPCClusterInfo(2,1) ) < 70 ) || ( ( nTrack->GetTPCClusterInfo(2,1) ) < 70 ) ) && lSmallestTrackLength<80 ) continue;
        
        //End track Quality Cuts
        //________________________________________________________________________
        
        lDcaPosToPrimVertex = TMath::Abs(pTrack->GetD(lBestPrimaryVtxPos[0],
                                                      lBestPrimaryVtxPos[1],
                                                      lMagneticField) );
        
        lDcaNegToPrimVertex = TMath::Abs(nTrack->GetD(lBestPrimaryVtxPos[0],
                                                      lBestPrimaryVtxPos[1],
                                                      lMagneticField) );
        
        lOnFlyStatus = v0->GetOnFlyStatus();
        lChi2V0 = v0->GetChi2V0();
        lDcaV0Daughters = v0->GetDcaV0Daughters();
        lDcaV0ToPrimVertex = v0->GetD(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2]);
        lV0CosineOfPointingAngle = v0->GetV0CosineOfPointingAngle(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2]);
        fTreeVariableV0CosineOfPointingAngle=lV0CosineOfPointingAngle;
        
        // Getting invariant mass infos directly from ESD
        v0->ChangeMassHypothesis(310);
        lInvMassK0s = v0->GetEffMass();
        v0->ChangeMassHypothesis(3122);
        lInvMassLambda = v0->GetEffMass();
        v0->ChangeMassHypothesis(-3122);
        lInvMassAntiLambda = v0->GetEffMass();
        lAlphaV0 = v0->AlphaV0();
        lPtArmV0 = v0->PtArmV0();
        
        fTreeVariableMVPileupFlag = fMVPileupFlag;
        
        fTreeVariablePt = v0->Pt();
        fTreeVariableChi2V0 = lChi2V0;
        fTreeVariableDcaV0ToPrimVertex = lDcaV0ToPrimVertex;
        fTreeVariableDcaV0Daughters = lDcaV0Daughters;
        fTreeVariableV0CosineOfPointingAngle = lV0CosineOfPointingAngle;
        fTreeVariableV0Radius = lV0Radius;
        fTreeVariableDcaPosToPrimVertex = lDcaPosToPrimVertex;
        fTreeVariableDcaNegToPrimVertex = lDcaNegToPrimVertex;
        fTreeVariableInvMassK0s = lInvMassK0s;
        fTreeVariableInvMassLambda = lInvMassLambda;
        fTreeVariableInvMassAntiLambda = lInvMassAntiLambda;
        fTreeVariableRapK0Short = lRapK0Short;
        fTreeVariableRapLambda = lRapLambda;
        fTreeVariableAlphaV0 = lAlphaV0;
        fTreeVariablePtArmV0 = lPtArmV0;
        
        //Official means of acquiring N-sigmas
        fTreeVariableNSigmasPosProton = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kProton );
        fTreeVariableNSigmasPosPion   = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kPion );
        fTreeVariableNSigmasNegProton = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kProton );
        fTreeVariableNSigmasNegPion   = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kPion );
        
        //This requires an Invariant Mass Hypothesis afterwards
        fTreeVariableDistOverTotMom = TMath::Sqrt(
                                                  TMath::Power( tDecayVertexV0[0] - lBestPrimaryVtxPos[0] , 2) +
                                                  TMath::Power( tDecayVertexV0[1] - lBestPrimaryVtxPos[1] , 2) +
                                                  TMath::Power( tDecayVertexV0[2] - lBestPrimaryVtxPos[2] , 2)
                                                  );
        fTreeVariableDistOverTotMom /= (lV0TotalMomentum+1e-10); //avoid division by zero, to be sure
        
        //Copy Multiplicity information
        fTreeVariableCentrality_V0A = fCentrality_V0A;
        fTreeVariableCentrality_V0C = fCentrality_V0C;
        fTreeVariableCentrality_V0M = fCentrality_V0M;
        
        
        //Info for pileup studies
        fTreeVariableNegTOFExpTDiff = nTrack->GetTOFExpTDiff( lESDevent->GetMagneticField() );
        fTreeVariablePosTOFExpTDiff = pTrack->GetTOFExpTDiff( lESDevent->GetMagneticField() );
        fTreeVariableNegTOFSignal = nTrack->GetTOFsignal() * 1.e-3; // in ns
        fTreeVariablePosTOFSignal = pTrack->GetTOFsignal() * 1.e-3; // in ns
        fTreeVariableNegTOFBCid = nTrack->GetTOFBunchCrossing( lESDevent->GetMagneticField() );
        fTreeVariablePosTOFBCid = pTrack->GetTOFBunchCrossing( lESDevent->GetMagneticField() );
        //Copy OOB pileup flag for this event
        fTreeVariableOOBPileupFlag = fOOBPileupFlag;
        //Copy VZERO information for this event
        fTreeVariableAmplitudeV0A = fAmplitudeV0A;
        fTreeVariableAmplitudeV0C = fAmplitudeV0C;
        //Copy IR information for this event
        fTreeVariableClosestNonEmptyBC = fClosestNonEmptyBC;
        
        
        //------------------------------------------------
        // Fill Tree!
        //------------------------------------------------
        
        // The conditionals are meant to decrease excessive
        // memory usage!
        
        //First Selection: Reject OnFly
        if( lOnFlyStatus == 0 ) {
            //Second Selection: rough 20-sigma band, parametric.
            //K0Short: Enough to parametrize peak broadening with linear function.
            Double_t lUpperLimitK0Short = (5.63707e-01) + (1.14979e-02)*fTreeVariablePt;
            Double_t lLowerLimitK0Short = (4.30006e-01) - (1.10029e-02)*fTreeVariablePt;
            //Lambda: Linear (for higher pt) plus exponential (for low-pt broadening)
            //[0]+[1]*x+[2]*TMath::Exp(-[3]*x)
            Double_t lUpperLimitLambda = (1.13688e+00) + (5.27838e-03)*fTreeVariablePt + (8.42220e-02)*TMath::Exp(-(3.80595e+00)*fTreeVariablePt);
            Double_t lLowerLimitLambda = (1.09501e+00) - (5.23272e-03)*fTreeVariablePt - (7.52690e-02)*TMath::Exp(-(3.46339e+00)*fTreeVariablePt);
            //Do Selection
            if(
               //Case 1: Lambda Selection
               (fTreeVariableInvMassLambda    < lUpperLimitLambda  && fTreeVariableInvMassLambda     > lLowerLimitLambda &&
                (!fkPreselectDedx || (fkPreselectDedx&&TMath::Abs(fTreeVariableNSigmasPosProton) < 7.0 && TMath::Abs(fTreeVariableNSigmasNegPion) < 7.0) )
                )
               ||
               //Case 2: AntiLambda Selection
               (fTreeVariableInvMassAntiLambda < lUpperLimitLambda  && fTreeVariableInvMassAntiLambda > lLowerLimitLambda &&
                (!fkPreselectDedx || (fkPreselectDedx&&TMath::Abs(fTreeVariableNSigmasNegProton) < 7.0 && TMath::Abs(fTreeVariableNSigmasPosPion) < 7.0) )
                )
               ||
               //Case 3: K0Short Selection
               (fTreeVariableInvMassK0s        < lUpperLimitK0Short && fTreeVariableInvMassK0s        > lLowerLimitK0Short &&
                (!fkPreselectDedx || (fkPreselectDedx&&TMath::Abs(fTreeVariableNSigmasNegPion)   < 7.0 && TMath::Abs(fTreeVariableNSigmasPosPion) < 7.0) )
                ) ) {
                   //Pre-selection in case this is AA...
                   
                   //Random denial
                   Bool_t lKeepV0 = kTRUE;
                   if(fkDownScaleV0 && ( fRand->Uniform() > fDownScaleFactorV0 )) lKeepV0 = kFALSE;
                   
                   //pT window
                   if( fTreeVariablePt < fMinPtToSave ) lKeepV0 = kFALSE;
                   if( fTreeVariablePt > fMaxPtToSave ) lKeepV0 = kFALSE;
                   
                   if ( TMath::Abs(fTreeVariableNegEta)<0.8 && TMath::Abs(fTreeVariablePosEta)<0.8 && fkSaveV0Tree && lKeepV0 ) fTreeV0->Fill();
               }
        }
        
        //------------------------------------------------
        // Fill V0 tree over.
        //------------------------------------------------
        
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        // Superlight adaptive output mode
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        
        //Step 1: Sweep members of the output object TList and fill all of them as appropriate
        Int_t lNumberOfConfigurations = fListV0->GetEntries();
        //AliWarning(Form("[V0 Analyses] Processing different configurations (%i detected)",lNumberOfConfigurations));
        TH3F *histoout_V0A         = 0x0;
        TH3F *histoout_V0C         = 0x0;
        TH3F *histoout_V0M         = 0x0;
        AliV0Result *lV0Result = 0x0;
        for(Int_t lcfg=0; lcfg<lNumberOfConfigurations; lcfg++){
            lV0Result = (AliV0Result*) fListV0->At(lcfg);
            histoout_V0A  = lV0Result->GetHistogram();
            histoout_V0C  = lV0Result->GetHistogram();
            histoout_V0M  = lV0Result->GetHistogram();
            
            
            Float_t lMass = 0;
            Float_t lRap  = 0;
            Float_t lPDGMass = -1;
            Float_t lNegdEdx = 100;
            Float_t lPosdEdx = 100;
            Float_t lBaryonMomentum = -0.5;

            //========================================================================
            //Setting up: Variable V0 CosPA
            Float_t lV0CosPACut = lV0Result -> GetCutV0CosPA();
            Float_t lVarV0CosPApar[5];
            lVarV0CosPApar[0] = lV0Result->GetCutVarV0CosPAExp0Const();
            lVarV0CosPApar[1] = lV0Result->GetCutVarV0CosPAExp0Slope();
            lVarV0CosPApar[2] = lV0Result->GetCutVarV0CosPAExp1Const();
            lVarV0CosPApar[3] = lV0Result->GetCutVarV0CosPAExp1Slope();
            lVarV0CosPApar[4] = lV0Result->GetCutVarV0CosPAConst();
            Float_t lVarV0CosPA = TMath::Cos(
                                             lVarV0CosPApar[0]*TMath::Exp(lVarV0CosPApar[1]*fTreeVariablePt) +
                                             lVarV0CosPApar[2]*TMath::Exp(lVarV0CosPApar[3]*fTreeVariablePt) +
                                             lVarV0CosPApar[4]);
            if( lV0Result->GetCutUseVarV0CosPA() ){
                //Only use if tighter than the non-variable cut
                if( lVarV0CosPA > lV0CosPACut ) lV0CosPACut = lVarV0CosPA;
            }
            //========================================================================
            
            if ( lV0Result->GetMassHypothesis() == AliV0Result::kK0Short     ){
                lMass    = fTreeVariableInvMassK0s;
                lRap     = fTreeVariableRapK0Short;
                lPDGMass = 0.497;
                lNegdEdx = fTreeVariableNSigmasNegPion;
                lPosdEdx = fTreeVariableNSigmasPosPion;
            }
            if ( lV0Result->GetMassHypothesis() == AliV0Result::kLambda      ){
                lMass = fTreeVariableInvMassLambda;
                lRap = fTreeVariableRapLambda;
                lPDGMass = 1.115683;
                lNegdEdx = fTreeVariableNSigmasNegPion;
                lPosdEdx = fTreeVariableNSigmasPosProton;
                lBaryonMomentum = fTreeVariablePosInnerP;
            }
            if ( lV0Result->GetMassHypothesis() == AliV0Result::kAntiLambda  ){
                lMass = fTreeVariableInvMassAntiLambda;
                lRap = fTreeVariableRapLambda;
                lPDGMass = 1.115683;
                lNegdEdx = fTreeVariableNSigmasNegProton;
                lPosdEdx = fTreeVariableNSigmasPosPion;
                lBaryonMomentum = fTreeVariableNegInnerP;
            }
            
            if (
                //Check 1: Offline Vertexer
                lOnFlyStatus == lV0Result->GetUseOnTheFly() &&
                
                //Check 2: Basic Acceptance cuts
                lV0Result->GetCutMinEtaTracks() < fTreeVariableNegEta && fTreeVariableNegEta < lV0Result->GetCutMaxEtaTracks() &&
                lV0Result->GetCutMinEtaTracks() < fTreeVariablePosEta && fTreeVariablePosEta < lV0Result->GetCutMaxEtaTracks() &&
                lRap > lV0Result->GetCutMinRapidity() &&
                lRap < lV0Result->GetCutMaxRapidity() &&
                
                //Check 3: Topological Variables
                fTreeVariableV0Radius > lV0Result->GetCutV0Radius() &&
                fTreeVariableV0Radius < lV0Result->GetCutMaxV0Radius() &&
                fTreeVariableDcaNegToPrimVertex > lV0Result->GetCutDCANegToPV() &&
                fTreeVariableDcaPosToPrimVertex > lV0Result->GetCutDCAPosToPV() &&
                fTreeVariableDcaV0Daughters < lV0Result->GetCutDCAV0Daughters() &&
                fTreeVariableV0CosineOfPointingAngle > lV0CosPACut &&
                fTreeVariableDistOverTotMom*lPDGMass < lV0Result->GetCutProperLifetime() &&
                fTreeVariableLeastNbrCrossedRows > lV0Result->GetCutLeastNumberOfCrossedRows() &&
                fTreeVariableLeastRatioCrossedRowsOverFindable > lV0Result->GetCutLeastNumberOfCrossedRowsOverFindable() &&
                
                //Check 4: Minimum momentum of baryon daughter
                ( lV0Result->GetMassHypothesis() == AliV0Result::kK0Short || lBaryonMomentum > lV0Result->GetCutMinBaryonMomentum() ) &&
                
                //Check 5: TPC dEdx selections
                TMath::Abs(lNegdEdx)<lV0Result->GetCutTPCdEdx() &&
                TMath::Abs(lPosdEdx)<lV0Result->GetCutTPCdEdx() &&
                
                //Check 6: Armenteros-Podolanski space cut (for K0Short analysis)
                ( ( lV0Result->GetCutArmenteros() == kFALSE || lV0Result->GetMassHypothesis() != AliV0Result::kK0Short ) || ( fTreeVariablePtArmV0>lV0Result->GetCutArmenterosParameter()*TMath::Abs(fTreeVariableAlphaV0) ) ) &&
                
                //Check 7: kITSrefit track selection if requested
                (
                 ( (fTreeVariableNegTrackStatus & AliESDtrack::kITSrefit) &&
                  (fTreeVariablePosTrackStatus & AliESDtrack::kITSrefit) )
                 ||
                 !lV0Result->GetCutUseITSRefitTracks()
                 )&&
                
                //Check 8: Max Chi2/Clusters if not absurd
                ( lV0Result->GetCutMaxChi2PerCluster()>1e+3 ||
                 (fTreeVariableMaxChi2PerCluster < lV0Result->GetCutMaxChi2PerCluster())
                 ) &&
                //Check 9: Min Track Length if positive
                ( lV0Result->GetCutMinTrackLength()<0 || //this is a bit paranoid...
                 fTreeVariableMinTrackLength > lV0Result->GetCutMinTrackLength()
                 )
            )
            {
                //This satisfies all my conditionals! Fill histogram
                histoout_V0A -> Fill ( fCentrality_V0A, fTreeVariablePt, lMass );
                histoout_V0C -> Fill ( fCentrality_V0C, fTreeVariablePt, lMass );
                histoout_V0M -> Fill ( fCentrality_V0M, fTreeVariablePt, lMass );
                
            }
        }
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        // End Superlight adaptive output mode
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        
    }// This is the end of the V0 loop
    
    //------------------------------------------------
    // Rerun cascade vertexer!
    //------------------------------------------------
    
    if( fkRunVertexers ) {
        //Remove existing cascades
        lESDevent->ResetCascades();
        
        //Decide between regular and light vertexer (default: light)
        if ( ! fkUseLightVertexer ){
            //Instantiate vertexer object
            AliCascadeVertexer lCascVtxer;
            lCascVtxer.SetDefaultCuts(fCascadeVertexerSels);
            lCascVtxer.SetCuts(fCascadeVertexerSels);
            lCascVtxer.V0sTracks2CascadeVertices(lESDevent);
        } else {
            AliLightCascadeVertexer lCascVtxer;
            lCascVtxer.SetDefaultCuts(fCascadeVertexerSels);
            lCascVtxer.SetCuts(fCascadeVertexerSels);
            if( fkUseOnTheFlyV0Cascading ) lCascVtxer.SetUseOnTheFlyV0(kTRUE);
            lCascVtxer.V0sTracks2CascadeVertices(lESDevent);
        }
    }
    
    //------------------------------------------------
    // MAIN CASCADE LOOP STARTS HERE
    //------------------------------------------------
    // Code Credit: Antonin Maire (thanks^100)
    // ---> This is an adaptation
    
    Long_t ncascades = 0;
    ncascades = lESDevent->GetNumberOfCascades();
                
    for (Int_t iXi = 0; iXi < ncascades; iXi++) {
        //------------------------------------------------
        // Initializations
        //------------------------------------------------
        //Double_t lTrkgPrimaryVtxRadius3D = -500.0;
        //Double_t lBestPrimaryVtxRadius3D = -500.0;
        
        // - 1st part of initialisation : variables needed to store AliESDCascade data members
        Double_t lEffMassXi      = 0. ;
        //Double_t lChi2Xi         = -1. ;
        Double_t lDcaXiDaughters = -1. ;
        Double_t lXiCosineOfPointingAngle = -1. ;
        Double_t lPosXi[3] = { -1000.0, -1000.0, -1000.0 };
        Double_t lXiRadius = -1000. ;
        
        // - 2nd part of initialisation : Nbr of clusters within TPC for the 3 daughter cascade tracks
        Int_t    lPosTPCClusters    = -1; // For ESD only ...//FIXME : wait for availability in AOD
        Int_t    lNegTPCClusters    = -1; // For ESD only ...
        Int_t    lBachTPCClusters   = -1; // For ESD only ...
        
        // - 3rd part of initialisation : about V0 part in cascades
        Double_t lInvMassLambdaAsCascDghter = 0.;
        //Double_t lV0Chi2Xi         = -1. ;
        Double_t lDcaV0DaughtersXi = -1.;
        
        Double_t lDcaBachToPrimVertexXi = -1., lDcaV0ToPrimVertexXi = -1.;
        Double_t lDcaPosToPrimVertexXi  = -1.;
        Double_t lDcaNegToPrimVertexXi  = -1.;
        Double_t lV0CosineOfPointingAngleXi = -1. ;
        Double_t lV0CosineOfPointingAngleXiSpecial = -1. ;
        Double_t lPosV0Xi[3] = { -1000. , -1000., -1000. }; // Position of VO coming from cascade
        Double_t lV0RadiusXi = -1000.0;
        Double_t lV0quality  = 0.;
        
        // - 4th part of initialisation : Effective masses
        Double_t lInvMassXiMinus    = 0.;
        Double_t lInvMassXiPlus     = 0.;
        Double_t lInvMassOmegaMinus = 0.;
        Double_t lInvMassOmegaPlus  = 0.;
        
        fTreeCascVarChiSquareV0      = 1e+3;
        fTreeCascVarChiSquareCascade = 1e+3;
        
        // - 6th part of initialisation : extra info for QA
        Double_t lXiMomX       = 0. , lXiMomY = 0., lXiMomZ = 0.;
        Double_t lXiTransvMom  = 0. ;
        //Double_t lXiTransvMomMC= 0. ;
        Double_t lXiTotMom     = 0. ;
        
        Double_t lBachMomX       = 0., lBachMomY  = 0., lBachMomZ   = 0.;
        //Double_t lBachTransvMom  = 0.;
        //Double_t lBachTotMom     = 0.;
        
        fTreeCascVarNegNSigmaPion   = -100;
        fTreeCascVarNegNSigmaProton = -100;
        fTreeCascVarPosNSigmaPion   = -100;
        fTreeCascVarPosNSigmaProton = -100;
        fTreeCascVarBachNSigmaPion  = -100;
        fTreeCascVarBachNSigmaKaon  = -100;
        
        fTreeCascVarNegTOFNSigmaPion   = -100;
        fTreeCascVarNegTOFNSigmaProton = -100;
        fTreeCascVarPosTOFNSigmaPion   = -100;
        fTreeCascVarPosTOFNSigmaProton = -100;
        fTreeCascVarBachTOFNSigmaPion  = -100;
        fTreeCascVarBachTOFNSigmaKaon  = -100;
        
        fTreeCascVarBachIsKink = kFALSE;
        fTreeCascVarPosIsKink = kFALSE;
        fTreeCascVarNegIsKink = kFALSE;
        
        //fTreeCascVarBachTotMom = -1;
        //fTreeCascVarPosTotMom  = -1;
        //fTreeCascVarNegTotMom  = -1;
        
        Short_t  lChargeXi = -2;
        //Double_t lV0toXiCosineOfPointingAngle = 0. ;
        
        Double_t lRapXi   = -20.0, lRapOmega = -20.0; //  lEta = -20.0, lTheta = 360., lPhi = 720. ;
        //Double_t lAlphaXi = -200., lPtArmXi  = -200.0;
        
        // -------------------------------------
        // II.ESD - Calculation Part dedicated to Xi vertices (ESD)
        
        AliESDcascade *xi = lESDevent->GetCascade(iXi);
        if (!xi) continue;
        
        // - II.Step 2 : Assigning the necessary variables for specific AliESDcascade data members (ESD)
        //-------------
        lV0quality = 0.;
        xi->ChangeMassHypothesis(lV0quality , 3312); // default working hypothesis : cascade = Xi- decay
        
        lEffMassXi  			= xi->GetEffMassXi();
        //ChiSquare implementation
        fTreeCascVarChiSquareV0      = xi->GetChi2V0();
        fTreeCascVarChiSquareCascade = xi->GetChi2Xi();
        
        lDcaXiDaughters 	= xi->GetDcaXiDaughters();
        lXiCosineOfPointingAngle 	            = xi->GetCascadeCosineOfPointingAngle( lBestPrimaryVtxPos[0],
                                                                                      lBestPrimaryVtxPos[1],
                                                                                      lBestPrimaryVtxPos[2] );
        // Take care : the best available vertex should be used (like in AliCascadeVertexer)
        
        xi->GetXYZcascade( lPosXi[0],  lPosXi[1], lPosXi[2] );
        lXiRadius			= TMath::Sqrt( lPosXi[0]*lPosXi[0]  +  lPosXi[1]*lPosXi[1] );
        
        fTreeCascVarCascadeDecayX = lPosXi[0];
        fTreeCascVarCascadeDecayY = lPosXi[1];
        fTreeCascVarCascadeDecayZ = lPosXi[2];
        
        // - II.Step 3 : around the tracks : Bach + V0 (ESD)
        // ~ Necessary variables for ESDcascade data members coming from the ESDv0 part (inheritance)
        //-------------
        
        UInt_t lIdxPosXi 	= (UInt_t) TMath::Abs( xi->GetPindex() );
        UInt_t lIdxNegXi 	= (UInt_t) TMath::Abs( xi->GetNindex() );
        UInt_t lBachIdx 	= (UInt_t) TMath::Abs( xi->GetBindex() );
        // Care track label can be negative in MC production (linked with the track quality)
        // However = normally, not the case for track index ...
        
        // FIXME : rejection of a double use of a daughter track (nothing but just a crosscheck of what is done in the cascade vertexer)
        if(lBachIdx == lIdxNegXi) {
            AliWarning("Pb / Idx(Bach. track) = Idx(Neg. track) ... continue!");
            continue;
        }
        if(lBachIdx == lIdxPosXi) {
            AliWarning("Pb / Idx(Bach. track) = Idx(Pos. track) ... continue!");
            continue;
        }
        
        AliESDtrack *pTrackXi		= lESDevent->GetTrack( lIdxPosXi );
        AliESDtrack *nTrackXi		= lESDevent->GetTrack( lIdxNegXi );
        AliESDtrack *bachTrackXi	= lESDevent->GetTrack( lBachIdx );
        
        fTreeCascVarNegIndex  = lIdxNegXi;
        fTreeCascVarPosIndex  = lIdxPosXi;
        fTreeCascVarBachIndex = lBachIdx;
        
        if (!pTrackXi || !nTrackXi || !bachTrackXi ) {
            AliWarning("ERROR: Could not retrieve one of the 3 ESD daughter tracks of the cascade ...");
            continue;
        }
        
        fTreeCascVarPosEta = pTrackXi->Eta();
        fTreeCascVarNegEta = nTrackXi->Eta();
        fTreeCascVarBachEta = bachTrackXi->Eta();
        
        //Check its clusters
        fTreeCascVarPosITSClusters0 = pTrackXi->HasPointOnITSLayer(0);
        fTreeCascVarPosITSClusters1 = pTrackXi->HasPointOnITSLayer(1);
        fTreeCascVarPosITSClusters2 = pTrackXi->HasPointOnITSLayer(2);
        fTreeCascVarPosITSClusters3 = pTrackXi->HasPointOnITSLayer(3);
        fTreeCascVarPosITSClusters4 = pTrackXi->HasPointOnITSLayer(4);
        fTreeCascVarPosITSClusters5 = pTrackXi->HasPointOnITSLayer(5);
        
        fTreeCascVarNegITSClusters0 = nTrackXi->HasPointOnITSLayer(0);
        fTreeCascVarNegITSClusters1 = nTrackXi->HasPointOnITSLayer(1);
        fTreeCascVarNegITSClusters2 = nTrackXi->HasPointOnITSLayer(2);
        fTreeCascVarNegITSClusters3 = nTrackXi->HasPointOnITSLayer(3);
        fTreeCascVarNegITSClusters4 = nTrackXi->HasPointOnITSLayer(4);
        fTreeCascVarNegITSClusters5 = nTrackXi->HasPointOnITSLayer(5);
        
        fTreeCascVarBachITSClusters0 = bachTrackXi->HasPointOnITSLayer(0);
        fTreeCascVarBachITSClusters1 = bachTrackXi->HasPointOnITSLayer(1);
        fTreeCascVarBachITSClusters2 = bachTrackXi->HasPointOnITSLayer(2);
        fTreeCascVarBachITSClusters3 = bachTrackXi->HasPointOnITSLayer(3);
        fTreeCascVarBachITSClusters4 = bachTrackXi->HasPointOnITSLayer(4);
        fTreeCascVarBachITSClusters5 = bachTrackXi->HasPointOnITSLayer(5);
        
        //Check its clusters, shared
        fTreeCascVarPosITSSharedClusters0 = pTrackXi->HasSharedPointOnITSLayer(0);
        fTreeCascVarPosITSSharedClusters1 = pTrackXi->HasSharedPointOnITSLayer(1);
        fTreeCascVarPosITSSharedClusters2 = pTrackXi->HasSharedPointOnITSLayer(2);
        fTreeCascVarPosITSSharedClusters3 = pTrackXi->HasSharedPointOnITSLayer(3);
        fTreeCascVarPosITSSharedClusters4 = pTrackXi->HasSharedPointOnITSLayer(4);
        fTreeCascVarPosITSSharedClusters5 = pTrackXi->HasSharedPointOnITSLayer(5);
        
        fTreeCascVarNegITSSharedClusters0 = nTrackXi->HasSharedPointOnITSLayer(0);
        fTreeCascVarNegITSSharedClusters1 = nTrackXi->HasSharedPointOnITSLayer(1);
        fTreeCascVarNegITSSharedClusters2 = nTrackXi->HasSharedPointOnITSLayer(2);
        fTreeCascVarNegITSSharedClusters3 = nTrackXi->HasSharedPointOnITSLayer(3);
        fTreeCascVarNegITSSharedClusters4 = nTrackXi->HasSharedPointOnITSLayer(4);
        fTreeCascVarNegITSSharedClusters5 = nTrackXi->HasSharedPointOnITSLayer(5);
        
        fTreeCascVarBachITSSharedClusters0 = bachTrackXi->HasSharedPointOnITSLayer(0);
        fTreeCascVarBachITSSharedClusters1 = bachTrackXi->HasSharedPointOnITSLayer(1);
        fTreeCascVarBachITSSharedClusters2 = bachTrackXi->HasSharedPointOnITSLayer(2);
        fTreeCascVarBachITSSharedClusters3 = bachTrackXi->HasSharedPointOnITSLayer(3);
        fTreeCascVarBachITSSharedClusters4 = bachTrackXi->HasSharedPointOnITSLayer(4);
        fTreeCascVarBachITSSharedClusters5 = bachTrackXi->HasSharedPointOnITSLayer(5);
        
        //GetKinkIndex condition
        if( bachTrackXi->GetKinkIndex(0)>0 ) fTreeCascVarBachIsKink = kTRUE;
        if( pTrackXi->GetKinkIndex(0)>0 ) fTreeCascVarPosIsKink = kTRUE;
        if( nTrackXi->GetKinkIndex(0)>0 ) fTreeCascVarNegIsKink = kTRUE;
        
        //Get track uncertainties
        //WARNING: THIS REFERS TO THE UNCERTAINTIES CLOSEST TO THE PV
        fTreeCascVarNegDCAPVSigmaX2 = TMath::Power(TMath::Sin(nTrackXi->GetAlpha()),2)*nTrackXi->GetSigmaY2();
        fTreeCascVarNegDCAPVSigmaY2 = TMath::Power(TMath::Cos(nTrackXi->GetAlpha()),2)*nTrackXi->GetSigmaY2();
        fTreeCascVarNegDCAPVSigmaZ2 = nTrackXi->GetSigmaZ2();
        
        fTreeCascVarPosDCAPVSigmaX2 = TMath::Power(TMath::Sin(pTrackXi->GetAlpha()),2)*pTrackXi->GetSigmaY2();
        fTreeCascVarPosDCAPVSigmaY2 = TMath::Power(TMath::Cos(pTrackXi->GetAlpha()),2)*pTrackXi->GetSigmaY2();
        fTreeCascVarPosDCAPVSigmaZ2 = pTrackXi->GetSigmaZ2();
        
        fTreeCascVarBachDCAPVSigmaX2 = TMath::Power(TMath::Sin(bachTrackXi->GetAlpha()),2)*bachTrackXi->GetSigmaY2();
        fTreeCascVarBachDCAPVSigmaY2 = TMath::Power(TMath::Cos(bachTrackXi->GetAlpha()),2)*bachTrackXi->GetSigmaY2();
        fTreeCascVarBachDCAPVSigmaZ2 = bachTrackXi->GetSigmaZ2();
        
        Double_t lBMom[3], lNMom[3], lPMom[3];
        xi->GetBPxPyPz( lBMom[0], lBMom[1], lBMom[2] );
        xi->GetPPxPyPz( lPMom[0], lPMom[1], lPMom[2] );
        xi->GetNPxPyPz( lNMom[0], lNMom[1], lNMom[2] );
        
        //fTreeCascVarBachTotMom = TMath::Sqrt( lBMom[0]*lBMom[0] + lBMom[1]*lBMom[1] + lBMom[2]*lBMom[2] );
        //fTreeCascVarPosTotMom  = TMath::Sqrt( lPMom[0]*lPMom[0] + lPMom[1]*lPMom[1] + lPMom[2]*lPMom[2] );
        //fTreeCascVarNegTotMom  = TMath::Sqrt( lNMom[0]*lNMom[0] + lNMom[1]*lNMom[1] + lNMom[2]*lNMom[2] );
        
        fTreeCascVarNegPx = lNMom[0];
        fTreeCascVarNegPy = lNMom[1];
        fTreeCascVarNegPz = lNMom[2];
        fTreeCascVarPosPx = lPMom[0];
        fTreeCascVarPosPy = lPMom[1];
        fTreeCascVarPosPz = lPMom[2];
        fTreeCascVarBachPx = lBMom[0];
        fTreeCascVarBachPy = lBMom[1];
        fTreeCascVarBachPz = lBMom[2];
        
        Int_t lNegTrackSign = 1;
        Int_t lPosTrackSign = 1;
        Int_t lBachTrackSign = 1;
        
        if( nTrackXi->GetSign() < 0 ) lNegTrackSign = -1;
        if( nTrackXi->GetSign() > 0 ) lNegTrackSign = +1;
        
        if( pTrackXi->GetSign() < 0 ) lPosTrackSign = -1;
        if( pTrackXi->GetSign() > 0 ) lPosTrackSign = +1;
        
        if( bachTrackXi->GetSign() < 0 ) lBachTrackSign = -1;
        if( bachTrackXi->GetSign() > 0 ) lBachTrackSign = +1;
        
        //------------------------------------------------
        // TPC dEdx information
        //------------------------------------------------
        fTreeCascVarNegNSigmaPion   = fPIDResponse->NumberOfSigmasTPC( nTrackXi, AliPID::kPion   );
        fTreeCascVarNegNSigmaProton = fPIDResponse->NumberOfSigmasTPC( nTrackXi, AliPID::kProton );
        fTreeCascVarPosNSigmaPion   = fPIDResponse->NumberOfSigmasTPC( pTrackXi, AliPID::kPion );
        fTreeCascVarPosNSigmaProton = fPIDResponse->NumberOfSigmasTPC( pTrackXi, AliPID::kProton );
        fTreeCascVarBachNSigmaPion  = fPIDResponse->NumberOfSigmasTPC( bachTrackXi, AliPID::kPion );
        fTreeCascVarBachNSigmaKaon  = fPIDResponse->NumberOfSigmasTPC( bachTrackXi, AliPID::kKaon );
        
        //------------------------------------------------
        // TOF info (no correction for weak decay traj.)
        //------------------------------------------------
        fTreeCascVarNegTOFNSigmaPion   = fPIDResponse->NumberOfSigmasTOF( nTrackXi, AliPID::kPion   );
        fTreeCascVarNegTOFNSigmaProton = fPIDResponse->NumberOfSigmasTOF( nTrackXi, AliPID::kProton );
        fTreeCascVarPosTOFNSigmaPion   = fPIDResponse->NumberOfSigmasTOF( pTrackXi, AliPID::kPion );
        fTreeCascVarPosTOFNSigmaProton = fPIDResponse->NumberOfSigmasTOF( pTrackXi, AliPID::kProton );
        fTreeCascVarBachTOFNSigmaPion  = fPIDResponse->NumberOfSigmasTOF( bachTrackXi, AliPID::kPion );
        fTreeCascVarBachTOFNSigmaKaon  = fPIDResponse->NumberOfSigmasTOF( bachTrackXi, AliPID::kKaon );
        
        fTreeCascVarNegITSNSigmaPion   = fPIDResponse->NumberOfSigmasITS( nTrackXi, AliPID::kPion   );
        fTreeCascVarNegITSNSigmaProton = fPIDResponse->NumberOfSigmasITS( nTrackXi, AliPID::kProton );
        fTreeCascVarPosITSNSigmaPion   = fPIDResponse->NumberOfSigmasITS( pTrackXi, AliPID::kPion );
        fTreeCascVarPosITSNSigmaProton = fPIDResponse->NumberOfSigmasITS( pTrackXi, AliPID::kProton );
        fTreeCascVarBachITSNSigmaPion  = fPIDResponse->NumberOfSigmasITS( bachTrackXi, AliPID::kPion );
        fTreeCascVarBachITSNSigmaKaon  = fPIDResponse->NumberOfSigmasITS( bachTrackXi, AliPID::kKaon );
        
        //------------------------------------------------
        // Raw TPC dEdx + PIDForTracking information
        //------------------------------------------------
        
        //Step 1: Acquire TPC inner wall total momentum
        const AliExternalTrackParam *innerneg=nTrackXi->GetInnerParam();
        const AliExternalTrackParam *innerpos=pTrackXi->GetInnerParam();
        const AliExternalTrackParam *innerbach=bachTrackXi->GetInnerParam();
        fTreeCascVarPosInnerP = -1;
        fTreeCascVarNegInnerP = -1;
        fTreeCascVarBachInnerP = -1;
        if(innerpos)  { fTreeCascVarPosInnerP  = innerpos ->GetP(); }
        if(innerneg)  { fTreeCascVarNegInnerP  = innerneg ->GetP(); }
        if(innerbach) { fTreeCascVarBachInnerP = innerbach->GetP(); }
        
        //Step 2: Acquire TPC Signals
        fTreeCascVarPosdEdx = pTrackXi->GetTPCsignal();
        fTreeCascVarNegdEdx = nTrackXi->GetTPCsignal();
        fTreeCascVarBachdEdx = bachTrackXi->GetTPCsignal();
        
        //Step 3: Acquire PID For Tracking
        fTreeCascVarPosPIDForTracking = pTrackXi->GetPIDForTracking();
        fTreeCascVarNegPIDForTracking = nTrackXi->GetPIDForTracking();
        fTreeCascVarBachPIDForTracking = bachTrackXi->GetPIDForTracking();
        
        //------------------------------------------------
        // TPC Number of clusters info
        // --- modified to save the smallest number
        // --- of TPC clusters for the 3 tracks
        //------------------------------------------------
        
        lPosTPCClusters   = pTrackXi->GetTPCNcls();
        lNegTPCClusters   = nTrackXi->GetTPCNcls();
        lBachTPCClusters  = bachTrackXi->GetTPCNcls();
        
        // 1 - Poor quality related to TPCrefit
        ULong_t pStatus    = pTrackXi->GetStatus();
        ULong_t nStatus    = nTrackXi->GetStatus();
        ULong_t bachStatus = bachTrackXi->GetStatus();
        
        //fTreeCascVarkITSRefitBachelor = kTRUE;
        //fTreeCascVarkITSRefitNegative = kTRUE;
        //fTreeCascVarkITSRefitPositive = kTRUE;
        
        if ((pStatus&AliESDtrack::kTPCrefit)    == 0) {
            AliDebug(1, "Pb / V0 Pos. track has no TPCrefit ... continue!");
            continue;
        }
        if ((nStatus&AliESDtrack::kTPCrefit)    == 0) {
            AliDebug(1, "Pb / V0 Neg. track has no TPCrefit ... continue!");
            continue;
        }
        if ((bachStatus&AliESDtrack::kTPCrefit) == 0) {
            AliDebug(1, "Pb / Bach.   track has no TPCrefit ... continue!");
            continue;
        }
        
        //Get status flags
        fTreeCascVarPosTrackStatus = pTrackXi->GetStatus();
        fTreeCascVarNegTrackStatus = nTrackXi->GetStatus();
        fTreeCascVarBachTrackStatus = bachTrackXi->GetStatus();
        
        fTreeCascVarPosDCAz = GetDCAz(pTrackXi);
        fTreeCascVarNegDCAz = GetDCAz(nTrackXi);
        fTreeCascVarBachDCAz = GetDCAz(bachTrackXi);
        
        Float_t lPosChi2PerCluster = pTrackXi->GetTPCchi2() / ((Float_t) lPosTPCClusters);
        Float_t lNegChi2PerCluster = nTrackXi->GetTPCchi2() / ((Float_t) lNegTPCClusters);
        Float_t lBachChi2PerCluster = bachTrackXi->GetTPCchi2() / ((Float_t) lBachTPCClusters);
        
        Int_t leastnumberofclusters = 1000;
        Float_t lBiggestChi2PerCluster = -1;
        
        //Pick minimum
        if( lPosTPCClusters < leastnumberofclusters ) leastnumberofclusters = lPosTPCClusters;
        if( lNegTPCClusters < leastnumberofclusters ) leastnumberofclusters = lNegTPCClusters;
        if( lBachTPCClusters < leastnumberofclusters ) leastnumberofclusters = lBachTPCClusters;
        
        //Pick maximum
        if( lPosChi2PerCluster > lBiggestChi2PerCluster ) lBiggestChi2PerCluster = lPosChi2PerCluster;
        if( lNegChi2PerCluster > lBiggestChi2PerCluster ) lBiggestChi2PerCluster = lNegChi2PerCluster;
        if( lBachChi2PerCluster > lBiggestChi2PerCluster ) lBiggestChi2PerCluster = lBachChi2PerCluster;
        
        //Extra track quality: min track length
        Float_t lSmallestTrackLength = 1000;
        Float_t lPosTrackLength = -1;
        Float_t lNegTrackLength = -1;
        Float_t lBachTrackLength = -1;
        
        if (pTrackXi->GetInnerParam()) lPosTrackLength = pTrackXi->GetLengthInActiveZone(1, 2.0, 220.0, lESDevent->GetMagneticField());
        if (nTrackXi->GetInnerParam()) lNegTrackLength = nTrackXi->GetLengthInActiveZone(1, 2.0, 220.0, lESDevent->GetMagneticField());
        if (bachTrackXi->GetInnerParam()) lBachTrackLength = bachTrackXi->GetLengthInActiveZone(1, 2.0, 220.0, lESDevent->GetMagneticField());
        
        if ( lPosTrackLength  < lSmallestTrackLength ) lSmallestTrackLength = lPosTrackLength;
        if ( lNegTrackLength  < lSmallestTrackLength ) lSmallestTrackLength = lNegTrackLength;
        if ( lBachTrackLength  < lSmallestTrackLength ) lSmallestTrackLength = lBachTrackLength;
        
        fTreeCascVarMinTrackLength = lSmallestTrackLength;
        
        // 2 - Poor quality related to TPC clusters: lowest cut of 70 clusters
        if(lPosTPCClusters  < 70 && lSmallestTrackLength < 80) {
            AliDebug(1, "Pb / V0 Pos. track has less than 70 TPC clusters ... continue!");
            continue;
        }
        if(lNegTPCClusters  < 70 && lSmallestTrackLength < 80) {
            AliDebug(1, "Pb / V0 Neg. track has less than 70 TPC clusters ... continue!");
            continue;
        }
        if(lBachTPCClusters < 70 && lSmallestTrackLength < 80) {
            AliDebug(1, "Pb / Bach.   track has less than 70 TPC clusters ... continue!");
            continue;
        }
        
        lInvMassLambdaAsCascDghter	= xi->GetEffMass();
        // This value shouldn't change, whatever the working hyp. is : Xi-, Xi+, Omega-, Omega+
        lDcaV0DaughtersXi 		= xi->GetDcaV0Daughters();
        //lV0Chi2Xi 			= xi->GetChi2V0();
        
        
        lV0CosineOfPointingAngleXi 	= xi->GetV0CosineOfPointingAngle( lBestPrimaryVtxPos[0],
                                                                     lBestPrimaryVtxPos[1],
                                                                     lBestPrimaryVtxPos[2] );
        //Modification: V0 CosPA wrt to Cascade decay vertex
        lV0CosineOfPointingAngleXiSpecial 	= xi->GetV0CosineOfPointingAngle( lPosXi[0],
                                                                             lPosXi[1],
                                                                             lPosXi[2] );
        
        lDcaV0ToPrimVertexXi 		= xi->GetD( lBestPrimaryVtxPos[0],
                                               lBestPrimaryVtxPos[1],
                                               lBestPrimaryVtxPos[2] );
        
        lDcaBachToPrimVertexXi = TMath::Abs( bachTrackXi->GetD(	lBestPrimaryVtxPos[0],
                                                               lBestPrimaryVtxPos[1],
                                                               lMagneticField  ) );
        // Note : AliExternalTrackParam::GetD returns an algebraic value ...
        
        xi->GetXYZ( lPosV0Xi[0],  lPosV0Xi[1], lPosV0Xi[2] );
        lV0RadiusXi		= TMath::Sqrt( lPosV0Xi[0]*lPosV0Xi[0]  +  lPosV0Xi[1]*lPosV0Xi[1] );
        
        fTreeCascVarV0DecayX = lPosV0Xi[0];
        fTreeCascVarV0DecayY = lPosV0Xi[1];
        fTreeCascVarV0DecayZ = lPosV0Xi[2];
        
        //========================================================================================
        //Calculate V0 lifetime for adaptive decay radius cut
        //3D Distance travelled by the V0 in the cascade
        Float_t lV0DistanceTrav =  TMath::Sqrt(  TMath::Power( lPosV0Xi[0]-lPosXi[0] , 2)
                                               + TMath::Power( lPosV0Xi[1]-lPosXi[1] , 2)
                                               + TMath::Power( lPosV0Xi[2]-lPosXi[2] , 2) );
        
        //Total V0 momentum
        Float_t lV0TotMomentum = TMath::Sqrt(  TMath::Power( lNMom[0]+lPMom[0] , 2)
                                             + TMath::Power( lNMom[1]+lPMom[1] , 2)
                                             + TMath::Power( lNMom[2]+lPMom[2] , 2) );
        
        //V0 transverse momentum
        Float_t lV0Pt = TMath::Sqrt(  TMath::Power( lNMom[0]+lPMom[0] , 2)
                                    + TMath::Power( lNMom[1]+lPMom[1] , 2) );
        
        //Calculate V0 lifetime: mL/p
        if( TMath::Abs(lV0TotMomentum)>1e-5 ){
            fTreeCascVarV0Lifetime = 1.115683*lV0DistanceTrav / lV0TotMomentum;
        }else{
            fTreeCascVarV0Lifetime = -1;
        }
        //========================================================================================
        
        lDcaPosToPrimVertexXi 	= TMath::Abs( pTrackXi	->GetD(	lBestPrimaryVtxPos[0],
                                                               lBestPrimaryVtxPos[1],
                                                               lMagneticField  )     );
        
        lDcaNegToPrimVertexXi 	= TMath::Abs( nTrackXi	->GetD(	lBestPrimaryVtxPos[0],
                                                               lBestPrimaryVtxPos[1],
                                                               lMagneticField  )     );
        
        // - II.Step 4 : around effective masses (ESD)
        // ~ change mass hypotheses to cover all the possibilities :  Xi-/+, Omega -/+
        
        if( bachTrackXi->Charge() < 0 )	{
            lV0quality = 0.;
            xi->ChangeMassHypothesis(lV0quality , 3312);
            // Calculate the effective mass of the Xi- candidate.
            // pdg code 3312 = Xi-
            lInvMassXiMinus = xi->GetEffMassXi();
            
            lV0quality = 0.;
            xi->ChangeMassHypothesis(lV0quality , 3334);
            // Calculate the effective mass of the Xi- candidate.
            // pdg code 3334 = Omega-
            lInvMassOmegaMinus = xi->GetEffMassXi();
            
            lV0quality = 0.;
            xi->ChangeMassHypothesis(lV0quality , 3312); 	// Back to default hyp.
        }// end if negative bachelor
        
        
        if( bachTrackXi->Charge() >  0 ) {
            lV0quality = 0.;
            xi->ChangeMassHypothesis(lV0quality , -3312);
            // Calculate the effective mass of the Xi+ candidate.
            // pdg code -3312 = Xi+
            lInvMassXiPlus = xi->GetEffMassXi();
            
            lV0quality = 0.;
            xi->ChangeMassHypothesis(lV0quality , -3334);
            // Calculate the effective mass of the Xi+ candidate.
            // pdg code -3334  = Omega+
            lInvMassOmegaPlus = xi->GetEffMassXi();
            
            lV0quality = 0.;
            xi->ChangeMassHypothesis(lV0quality , -3312); 	// Back to "default" hyp.
        }// end if positive bachelor
        
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        //Recalculate from scratch, please
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        
        //WARNING: This will not be checked for correctness (charge-wise, etc)
        //         It will be up to the user to use the correct variable whenever needed!
        
        //+-+ Recalculate Xi Masses from scratch: will not change with lambda mass as
        //the perfect lambda mass is always assumed
        
        //+-+ Recalculate Lambda mass from scratch
        //Under Lambda hypothesis, the positive daughter is the proton, negative pion
        Double_t m1 = TDatabasePDG::Instance()->GetParticle(kProton)->Mass();
        Double_t m2 = TDatabasePDG::Instance()->GetParticle(kPiPlus)->Mass();
        Double_t e12   = m1*m1+lPMom[0]*lPMom[0]+lPMom[1]*lPMom[1]+lPMom[2]*lPMom[2];
        Double_t e22   = m2*m2+lNMom[0]*lNMom[0]+lNMom[1]*lNMom[1]+lNMom[2]*lNMom[2];
        fTreeCascVarV0MassLambda = TMath::Sqrt(TMath::Max(m1*m1+m2*m2
                                                          +2.*(TMath::Sqrt(e12*e22)-lPMom[0]*lNMom[0]-lPMom[1]*lNMom[1]-lPMom[2]*lNMom[2]),0.));
        
        //+-+ Recalculate AntiLambda mass from scratch
        //Under Lambda hypothesis, the positive daughter is the pion, negative antiproton
        m1 = TDatabasePDG::Instance()->GetParticle(kPiPlus)->Mass();
        m2 = TDatabasePDG::Instance()->GetParticle(kProton)->Mass();
        e12   = m1*m1+lPMom[0]*lPMom[0]+lPMom[1]*lPMom[1]+lPMom[2]*lPMom[2];
        e22   = m2*m2+lNMom[0]*lNMom[0]+lNMom[1]*lNMom[1]+lNMom[2]*lNMom[2];
        fTreeCascVarV0MassAntiLambda = TMath::Sqrt(TMath::Max(m1*m1+m2*m2
                                                              +2.*(TMath::Sqrt(e12*e22)-lPMom[0]*lNMom[0]-lPMom[1]*lNMom[1]-lPMom[2]*lNMom[2]),0.));
        
        
        
        // - II.Step 6 : extra info for QA (ESD)
        // miscellaneous pieces of info that may help regarding data quality assessment.
        //-------------
        xi->GetPxPyPz( lXiMomX, lXiMomY, lXiMomZ );
        lXiTransvMom  	= TMath::Sqrt( lXiMomX*lXiMomX   + lXiMomY*lXiMomY );
        lXiTotMom  	= TMath::Sqrt( lXiMomX*lXiMomX   + lXiMomY*lXiMomY   + lXiMomZ*lXiMomZ );
        
        xi->GetBPxPyPz(  lBachMomX,  lBachMomY,  lBachMomZ );
        //lBachTransvMom  = TMath::Sqrt( lBachMomX*lBachMomX   + lBachMomY*lBachMomY );
        //lBachTotMom  	= TMath::Sqrt( lBachMomX*lBachMomX   + lBachMomY*lBachMomY  +  lBachMomZ*lBachMomZ  );
        
        lChargeXi = xi->Charge();
        
        //lV0toXiCosineOfPointingAngle = xi->GetV0CosineOfPointingAngle( lPosXi[0], lPosXi[1], lPosXi[2] );
        
        lRapXi    = xi->RapXi();
        lRapOmega = xi->RapOmega();
        //lEta      = xi->Eta();
        //lTheta    = xi->Theta() *180.0/TMath::Pi();
        //lPhi      = xi->Phi()   *180.0/TMath::Pi();
        //lAlphaXi  = xi->AlphaXi();
        //lPtArmXi  = xi->PtArmXi();
        
        //----------------------------------------
        // Bump studies: perform propagation
        //----------------------------------------
        
        AliESDtrack *lBaryonTrack = 0x0;
        AliESDtrack *lBachelorTrack = 0x0;
        if ( lChargeXi == -1 ){
            lBaryonTrack = pTrackXi;
            lBachelorTrack = bachTrackXi;
        }
        if ( lChargeXi == +1 ){
            lBaryonTrack = nTrackXi;
            lBachelorTrack = bachTrackXi;
        }
        
        fTreeCascVarDCABachToBaryon = -100;
        
        Double_t bMag = lESDevent->GetMagneticField();
        Double_t xn, xp;
        
        //Care has to be taken here
        if ( lBaryonTrack && lBachelorTrack ){
            //Attempt zero: Calculate DCA between bachelor and baryon daughter
            fTreeCascVarDCABachToBaryon = lBaryonTrack->GetDCA(lBachelorTrack, bMag, xn, xp);
        }
        
        fTreeCascVarWrongCosPA = -1;
        if( bachTrackXi->Charge() < 0 )
            fTreeCascVarWrongCosPA = GetCosPA( bachTrackXi , pTrackXi, lESDevent );
        if( bachTrackXi->Charge() > 0 )
            fTreeCascVarWrongCosPA = GetCosPA( bachTrackXi , nTrackXi, lESDevent );
        
        
        //------------------------------------------------
        // Set Variables for adding to tree
        //------------------------------------------------
        
        fTreeCascVarMVPileupFlag = fMVPileupFlag;
        
        fTreeCascVarCharge	= lChargeXi;
        if (lChargeXi < 0 ){
            fTreeCascVarMassAsXi    = lInvMassXiMinus;
            fTreeCascVarMassAsOmega = lInvMassOmegaMinus;
        }
        if (lChargeXi > 0 ){
            fTreeCascVarMassAsXi    = lInvMassXiPlus;
            fTreeCascVarMassAsOmega = lInvMassOmegaPlus;
        }
        fTreeCascVarPt = lXiTransvMom;
        fTreeCascVarRapXi = lRapXi ;
        fTreeCascVarRapOmega = lRapOmega ;
        fTreeCascVarDCACascDaughters = lDcaXiDaughters;
        fTreeCascVarDCABachToPrimVtx = lDcaBachToPrimVertexXi;
        fTreeCascVarDCAV0Daughters = lDcaV0DaughtersXi;
        fTreeCascVarDCAV0ToPrimVtx = lDcaV0ToPrimVertexXi;
        fTreeCascVarDCAPosToPrimVtx = lDcaPosToPrimVertexXi;
        fTreeCascVarDCANegToPrimVtx = lDcaNegToPrimVertexXi;
        fTreeCascVarCascCosPointingAngle = lXiCosineOfPointingAngle;
        fTreeCascVarCascRadius = lXiRadius;
        fTreeCascVarV0Mass = lInvMassLambdaAsCascDghter;
        fTreeCascVarV0CosPointingAngle = lV0CosineOfPointingAngleXi;
        fTreeCascVarV0CosPointingAngleSpecial = lV0CosineOfPointingAngleXiSpecial;
        fTreeCascVarV0Radius = lV0RadiusXi;
        fTreeCascVarLeastNbrClusters = leastnumberofclusters;
        fTreeCascVarMaxChi2PerCluster = lBiggestChi2PerCluster;
        
        //Copy Multiplicity information
        fTreeCascVarCentrality_V0A = fCentrality_V0A;
        fTreeCascVarCentrality_V0C = fCentrality_V0C;
        fTreeCascVarCentrality_V0M = fCentrality_V0M;
        
        
        fTreeCascVarDistOverTotMom = TMath::Sqrt(
                                                 TMath::Power( lPosXi[0] - lBestPrimaryVtxPos[0] , 2) +
                                                 TMath::Power( lPosXi[1] - lBestPrimaryVtxPos[1] , 2) +
                                                 TMath::Power( lPosXi[2] - lBestPrimaryVtxPos[2] , 2)
                                                 );
        fTreeCascVarDistOverTotMom /= (lXiTotMom+1e-13);
        
        //Info for pileup studies
        fTreeCascVarBachTOFExpTDiff = bachTrackXi->GetTOFExpTDiff( bMag );
        fTreeCascVarNegTOFExpTDiff = nTrackXi->GetTOFExpTDiff( bMag );
        fTreeCascVarPosTOFExpTDiff = pTrackXi->GetTOFExpTDiff( bMag );
        
        fTreeCascVarBachTOFSignal = bachTrackXi->GetTOFsignal() * 1.e-3; // in ns
        fTreeCascVarNegTOFSignal = nTrackXi->GetTOFsignal() * 1.e-3; // in ns
        fTreeCascVarPosTOFSignal = pTrackXi->GetTOFsignal() * 1.e-3; // in ns
        fTreeCascVarBachTOFBCid = bachTrackXi->GetTOFBunchCrossing( bMag );
        fTreeCascVarNegTOFBCid = nTrackXi->GetTOFBunchCrossing( bMag );
        fTreeCascVarPosTOFBCid = pTrackXi->GetTOFBunchCrossing( bMag );
        //Copy OOB pileup flag for this event
        fTreeCascVarOOBPileupFlag = fOOBPileupFlag;
        //Copy VZERO information for this event
        fTreeCascVarAmplitudeV0A = fAmplitudeV0A;
        fTreeCascVarAmplitudeV0C = fAmplitudeV0C;
        //Copy IR information for this event
        fTreeCascVarClosestNonEmptyBC = fClosestNonEmptyBC;
        
        if ( fkExtraCleanup ){
            //Meant to provide extra level of cleanup
            if( TMath::Abs(fTreeCascVarPosEta)>0.8 || TMath::Abs(fTreeCascVarNegEta)>0.8 || TMath::Abs(fTreeCascVarBachEta)>0.8 ) continue;
            if( TMath::Abs(fTreeCascVarRapXi)>0.5 && TMath::Abs(fTreeCascVarRapOmega)>0.5 ) continue;
            if ( fkPreselectDedx ){
                Bool_t lPassesPreFilterdEdx = kFALSE;
                //XiMinus Pre-selection
                if( fTreeCascVarMassAsXi<1.32+0.250&&fTreeCascVarMassAsXi>1.32-0.250 && TMath::Abs(fTreeCascVarPosNSigmaProton) < 5.0 && TMath::Abs(fTreeCascVarNegNSigmaPion) < 5.0 && TMath::Abs(fTreeCascVarBachNSigmaPion) < 5.0 && fTreeCascVarCharge == -1 ) lPassesPreFilterdEdx = kTRUE;
                if( fTreeCascVarMassAsXi<1.32+0.250&&fTreeCascVarMassAsXi>1.32-0.250 && TMath::Abs(fTreeCascVarPosNSigmaPion) < 5.0 && TMath::Abs(fTreeCascVarNegNSigmaProton) < 5.0 && TMath::Abs(fTreeCascVarBachNSigmaPion) < 5.0 && fTreeCascVarCharge == +1 ) lPassesPreFilterdEdx = kTRUE;
                if(fTreeCascVarMassAsOmega<1.68+0.250&&fTreeCascVarMassAsOmega>1.68-0.250 && TMath::Abs(fTreeCascVarPosNSigmaProton) < 5.0 && TMath::Abs(fTreeCascVarNegNSigmaPion) < 5.0 && TMath::Abs(fTreeCascVarBachNSigmaKaon) < 5.0 && fTreeCascVarCharge == -1  ) lPassesPreFilterdEdx = kTRUE;
                if(fTreeCascVarMassAsOmega<1.68+0.250&&fTreeCascVarMassAsOmega>1.68-0.250 && TMath::Abs(fTreeCascVarPosNSigmaPion) < 5.0 && TMath::Abs(fTreeCascVarNegNSigmaProton) < 5.0 && TMath::Abs(fTreeCascVarBachNSigmaKaon) < 5.0 && fTreeCascVarCharge == +1) lPassesPreFilterdEdx = kTRUE;
                if( !lPassesPreFilterdEdx ) continue;
            }
        }
        
        //All vars not specified here: specified elsewhere!
        
        //------------------------------------------------
        // Fill Tree!
        //------------------------------------------------
        
        // The conditional is meant to decrease excessive
        // memory usage! Be careful when loosening the
        // cut!
        
        //Xi    Mass window: 150MeV wide
        //Omega mass window: 150MeV wide
        
        //Random denial
        Bool_t lKeepCascade = kTRUE;
        if(fkDownScaleCascade && ( fRand->Uniform() > fDownScaleFactorCascade )) lKeepCascade = kFALSE;
        
        //Lowest pT cutoff (this is all background anyways)
        if( fTreeCascVarPt < fMinPtToSave ) lKeepCascade = kFALSE;
        if( fTreeCascVarPt > fMaxPtToSave ) lKeepCascade = kFALSE;
        
        if( fkSaveCascadeTree && lKeepCascade &&
           (
            (//START XI SELECTIONS
             (fTreeCascVarMassAsXi<1.32+0.075&&fTreeCascVarMassAsXi>1.32-0.075)
             )//end Xi Selections
            ||
            (//START OMEGA SELECTIONS
             (fTreeCascVarMassAsOmega<1.68+0.075&&fTreeCascVarMassAsOmega>1.68-0.075)
             )//end Xi Selections
            )
           )
        {
            fTreeCascade->Fill();
        }
        //------------------------------------------------
        // Fill tree over.
        //------------------------------------------------
        
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        // Superlight adaptive output mode
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        
        //Step 1: Sweep members of the output object TList and fill all of them as appropriate
        Int_t lNumberOfConfigurationsCascade = fListCascade->GetEntries();
        //AliWarning(Form("[Cascade Analyses] Processing different configurations (%i detected)",lNumberOfConfigurationsCascade));
        TH3F *histoout_V0A         = 0x0;
        TH3F *histoout_V0C         = 0x0;
        TH3F *histoout_V0M         = 0x0;
        
        AliCascadeResult *lCascadeResult = 0x0;
        for(Int_t lcfg=0; lcfg<lNumberOfConfigurationsCascade; lcfg++){
            lCascadeResult = (AliCascadeResult*) fListCascade->At(lcfg);
            histoout_V0A  = lCascadeResult->GetHistogram();
            histoout_V0C  = lCascadeResult->GetHistogram();
            histoout_V0M  = lCascadeResult->GetHistogram();
            
            
            Float_t lMass = 0;
            Float_t lV0Mass = 0;
            Float_t lRap  = 0;
            Float_t lPDGMass = -1;
            Float_t lNegdEdx = 100;
            Float_t lPosdEdx = 100;
            Float_t lBachdEdx = 100;
            Float_t lNegTOFsigma = 100;
            Float_t lPosTOFsigma = 100;
            Float_t lBachTOFsigma = 100;
            Short_t  lCharge = -2;
            Int_t lChargePos =  1;
            Int_t lChargeNeg = -1;
            Float_t lprpx, lprpy, lprpz, lpipx, lpipy, lpipz;
            lpipx = fTreeCascVarBachPx;
            lpipy = fTreeCascVarBachPy;
            lpipz = fTreeCascVarBachPz;
            
            //For parametric V0 Mass selection
            Float_t lExpV0Mass =
            fLambdaMassMean[0]+
            fLambdaMassMean[1]*TMath::Exp(fLambdaMassMean[2]*lV0Pt)+
            fLambdaMassMean[3]*TMath::Exp(fLambdaMassMean[4]*lV0Pt);
            
            Float_t lExpV0Sigma =
            fLambdaMassSigma[0]+fLambdaMassSigma[1]*lV0Pt+
            fLambdaMassSigma[2]*TMath::Exp(fLambdaMassSigma[3]*lV0Pt);
            
            //========================================================================
            //For 2.76TeV-like parametric V0 CosPA
            Float_t l276TeVV0CosPA = 0.998;
            Float_t pThr=1.5;
            if (lV0TotMomentum<pThr) {
                //Below the threshold "pThr", try a momentum dependent cos(PA) cut
                const Double_t bend=0.03; // approximate Xi bending angle
                const Double_t qt=0.211;  // max Lambda pT in Omega decay
                const Double_t cpaThr=TMath::Cos(TMath::ATan(qt/pThr) + bend);
                Double_t
                cpaCut=(0.998/cpaThr)*TMath::Cos(TMath::ATan(qt/lV0TotMomentum) + bend);
                l276TeVV0CosPA = cpaCut;
            }
            //========================================================================
            
            //========================================================================
            //Setting up: Variable Cascade CosPA
            Float_t lCascCosPACut = lCascadeResult -> GetCutCascCosPA();
            Float_t lVarCascCosPApar[5];
            lVarCascCosPApar[0] = lCascadeResult->GetCutVarCascCosPAExp0Const();
            lVarCascCosPApar[1] = lCascadeResult->GetCutVarCascCosPAExp0Slope();
            lVarCascCosPApar[2] = lCascadeResult->GetCutVarCascCosPAExp1Const();
            lVarCascCosPApar[3] = lCascadeResult->GetCutVarCascCosPAExp1Slope();
            lVarCascCosPApar[4] = lCascadeResult->GetCutVarCascCosPAConst();
            Float_t lVarCascCosPA = TMath::Cos(
                                               lVarCascCosPApar[0]*TMath::Exp(lVarCascCosPApar[1]*fTreeCascVarPt) +
                                               lVarCascCosPApar[2]*TMath::Exp(lVarCascCosPApar[3]*fTreeCascVarPt) +
                                               lVarCascCosPApar[4]);
            if( lCascadeResult->GetCutUseVarCascCosPA() ){
                //Only use if tighter than the non-variable cut
                if( lVarCascCosPA > lCascCosPACut ) lCascCosPACut = lVarCascCosPA;
            }
            //========================================================================
            
            //========================================================================
            //Setting up: Variable V0 CosPA
            Float_t lV0CosPACut = lCascadeResult -> GetCutV0CosPA();
            Float_t lVarV0CosPApar[5];
            lVarV0CosPApar[0] = lCascadeResult->GetCutVarV0CosPAExp0Const();
            lVarV0CosPApar[1] = lCascadeResult->GetCutVarV0CosPAExp0Slope();
            lVarV0CosPApar[2] = lCascadeResult->GetCutVarV0CosPAExp1Const();
            lVarV0CosPApar[3] = lCascadeResult->GetCutVarV0CosPAExp1Slope();
            lVarV0CosPApar[4] = lCascadeResult->GetCutVarV0CosPAConst();
            Float_t lVarV0CosPA = TMath::Cos(
                                             lVarV0CosPApar[0]*TMath::Exp(lVarV0CosPApar[1]*fTreeCascVarPt) +
                                             lVarV0CosPApar[2]*TMath::Exp(lVarV0CosPApar[3]*fTreeCascVarPt) +
                                             lVarV0CosPApar[4]);
            if( lCascadeResult->GetCutUseVarV0CosPA() ){
                //Only use if tighter than the non-variable cut
                if( lVarV0CosPA > lV0CosPACut ) lV0CosPACut = lVarV0CosPA;
            }
            //========================================================================
            
            //========================================================================
            //Setting up: Variable BB CosPA
            Float_t lBBCosPACut = lCascadeResult -> GetCutBachBaryonCosPA();
            Float_t lVarBBCosPApar[5];
            lVarBBCosPApar[0] = lCascadeResult->GetCutVarBBCosPAExp0Const();
            lVarBBCosPApar[1] = lCascadeResult->GetCutVarBBCosPAExp0Slope();
            lVarBBCosPApar[2] = lCascadeResult->GetCutVarBBCosPAExp1Const();
            lVarBBCosPApar[3] = lCascadeResult->GetCutVarBBCosPAExp1Slope();
            lVarBBCosPApar[4] = lCascadeResult->GetCutVarBBCosPAConst();
            Float_t lVarBBCosPA = TMath::Cos(
                                             lVarBBCosPApar[0]*TMath::Exp(lVarBBCosPApar[1]*fTreeCascVarPt) +
                                             lVarBBCosPApar[2]*TMath::Exp(lVarBBCosPApar[3]*fTreeCascVarPt) +
                                             lVarBBCosPApar[4]);
            if( lCascadeResult->GetCutUseVarBBCosPA() ){
                //Only use if looser than the non-variable cut (WARNING: BEWARE INVERSE LOGIC)
                if( lVarBBCosPA > lBBCosPACut ) lBBCosPACut = lVarBBCosPA;
            }
            //========================================================================
            
            if ( lCascadeResult->GetMassHypothesis() == AliCascadeResult::kXiMinus     ){
                lCharge  = -1;
                if ( lCascadeResult->GetSwapBachelorCharge() ) lCharge *= -1;
                lMass    = fTreeCascVarMassAsXi;
                lV0Mass  = fTreeCascVarV0MassLambda;
                lRap     = fTreeCascVarRapXi;
                lPDGMass = 1.32171;
                lNegdEdx = fTreeCascVarNegNSigmaPion;
                lPosdEdx = fTreeCascVarPosNSigmaProton;
                lBachdEdx= fTreeCascVarBachNSigmaPion;
                lNegTOFsigma = fTreeCascVarNegTOFNSigmaPion;
                lPosTOFsigma = fTreeCascVarPosTOFNSigmaProton;
                lBachTOFsigma = fTreeCascVarBachTOFNSigmaPion;
                lprpx = fTreeCascVarPosPx;
                lprpy = fTreeCascVarPosPy;
                lprpz = fTreeCascVarPosPz;
            }
            if ( lCascadeResult->GetMassHypothesis() == AliCascadeResult::kXiPlus      ){
                lCharge  = +1;
                if ( lCascadeResult->GetSwapBachelorCharge() ) lCharge *= -1;
                lMass    = fTreeCascVarMassAsXi;
                lV0Mass  = fTreeCascVarV0MassAntiLambda;
                lRap     = fTreeCascVarRapXi;
                lPDGMass = 1.32171;
                lNegdEdx = fTreeCascVarNegNSigmaProton;
                lPosdEdx = fTreeCascVarPosNSigmaPion;
                lBachdEdx= fTreeCascVarBachNSigmaPion;
                lNegTOFsigma = fTreeCascVarNegTOFNSigmaProton;
                lPosTOFsigma = fTreeCascVarPosTOFNSigmaPion;
                lBachTOFsigma = fTreeCascVarBachTOFNSigmaPion;
                lprpx = fTreeCascVarNegPx;
                lprpy = fTreeCascVarNegPy;
                lprpz = fTreeCascVarNegPz;
            }
            if ( lCascadeResult->GetMassHypothesis() == AliCascadeResult::kOmegaMinus     ){
                lCharge  = -1;
                if ( lCascadeResult->GetSwapBachelorCharge() ) lCharge *= -1;
                lMass    = fTreeCascVarMassAsOmega;
                lV0Mass  = fTreeCascVarV0MassLambda;
                lRap     = fTreeCascVarRapOmega;
                lPDGMass = 1.67245;
                lNegdEdx = fTreeCascVarNegNSigmaPion;
                lPosdEdx = fTreeCascVarPosNSigmaProton;
                lBachdEdx= fTreeCascVarBachNSigmaKaon;
                lNegTOFsigma = fTreeCascVarNegTOFNSigmaPion;
                lPosTOFsigma = fTreeCascVarPosTOFNSigmaProton;
                lBachTOFsigma = fTreeCascVarBachTOFNSigmaKaon;
                lprpx = fTreeCascVarPosPx;
                lprpy = fTreeCascVarPosPy;
                lprpz = fTreeCascVarPosPz;
            }
            if ( lCascadeResult->GetMassHypothesis() == AliCascadeResult::kOmegaPlus      ){
                lCharge  = +1;
                if ( lCascadeResult->GetSwapBachelorCharge() ) lCharge *= -1;
                lMass    = fTreeCascVarMassAsOmega;
                lV0Mass  = fTreeCascVarV0MassAntiLambda;
                lRap     = fTreeCascVarRapOmega;
                lPDGMass = 1.67245;
                lNegdEdx = fTreeCascVarNegNSigmaProton;
                lPosdEdx = fTreeCascVarPosNSigmaPion;
                lBachdEdx= fTreeCascVarBachNSigmaKaon;
                lNegTOFsigma = fTreeCascVarNegTOFNSigmaProton;
                lPosTOFsigma = fTreeCascVarPosTOFNSigmaPion;
                lBachTOFsigma = fTreeCascVarBachTOFNSigmaKaon;
                lprpx = fTreeCascVarNegPx;
                lprpy = fTreeCascVarNegPy;
                lprpz = fTreeCascVarNegPz;
            }
            if (lCascadeResult->GetCutUseTOFUnchecked() == kFALSE ){
                //Always-pass values
                lNegTOFsigma = 0;
                lPosTOFsigma = 0;
                lBachTOFsigma = 0;
            }
            
            if (
                //Check 1: Charge consistent with expectations
                fTreeCascVarCharge == lCharge &&
                
                //Check 2: Basic Acceptance cuts
                lCascadeResult->GetCutMinEtaTracks() < fTreeCascVarPosEta && fTreeCascVarPosEta < lCascadeResult->GetCutMaxEtaTracks() &&
                lCascadeResult->GetCutMinEtaTracks() < fTreeCascVarNegEta && fTreeCascVarNegEta < lCascadeResult->GetCutMaxEtaTracks() &&
                lCascadeResult->GetCutMinEtaTracks() < fTreeCascVarBachEta && fTreeCascVarBachEta < lCascadeResult->GetCutMaxEtaTracks() &&
                lRap > lCascadeResult->GetCutMinRapidity() &&
                lRap < lCascadeResult->GetCutMaxRapidity() &&
                
                //Check 3: Topological Variables
                // - V0 Selections
                fTreeCascVarDCANegToPrimVtx > lCascadeResult->GetCutDCANegToPV() &&
                fTreeCascVarDCAPosToPrimVtx > lCascadeResult->GetCutDCAPosToPV() &&
                fTreeCascVarDCAV0Daughters < lCascadeResult->GetCutDCAV0Daughters() &&
                fTreeCascVarV0CosPointingAngle > lV0CosPACut &&
                fTreeCascVarV0Radius > lCascadeResult->GetCutV0Radius() &&
                // - Cascade Selections
                fTreeCascVarDCAV0ToPrimVtx > lCascadeResult->GetCutDCAV0ToPV() &&
                TMath::Abs(lV0Mass-1.116) < lCascadeResult->GetCutV0Mass() &&
                fTreeCascVarDCABachToPrimVtx > lCascadeResult->GetCutDCABachToPV() &&
                fTreeCascVarDCACascDaughters < lCascadeResult->GetCutDCACascDaughters() &&
                fTreeCascVarCascCosPointingAngle > lCascCosPACut &&
                fTreeCascVarCascRadius > lCascadeResult->GetCutCascRadius() &&
                
                // - Implementation of a parametric V0 Mass cut if requested
                (
                 ( lCascadeResult->GetCutV0MassSigma() > 50 ) || //anything goes
                 (TMath::Abs( (lV0Mass-lExpV0Mass) / lExpV0Sigma ) < lCascadeResult->GetCutV0MassSigma() )
                 ) &&
                
                // - Miscellaneous
                fTreeCascVarDistOverTotMom*lPDGMass < lCascadeResult->GetCutProperLifetime() &&
                fTreeCascVarLeastNbrClusters > lCascadeResult->GetCutLeastNumberOfClusters() &&
                
                //Check 4: TPC dEdx selections
                TMath::Abs(lNegdEdx )<lCascadeResult->GetCutTPCdEdx() &&
                TMath::Abs(lPosdEdx )<lCascadeResult->GetCutTPCdEdx() &&
                TMath::Abs(lBachdEdx)<lCascadeResult->GetCutTPCdEdx() &&
                
                //Check 5: Xi rejection for Omega analysis
                ( ( lCascadeResult->GetMassHypothesis() != AliCascadeResult::kOmegaMinus && lCascadeResult->GetMassHypothesis() != AliCascadeResult::kOmegaPlus  ) || ( TMath::Abs( fTreeCascVarMassAsXi - 1.32171 ) > lCascadeResult->GetCutXiRejection() ) ) &&
                
                //Check 6: Experimental DCA Bachelor to Baryon cut
                ( fTreeCascVarDCABachToBaryon > lCascadeResult->GetCutDCABachToBaryon() ) &&
                
                //Check 7: Experimental Bach Baryon CosPA
                ( fTreeCascVarWrongCosPA < lBBCosPACut  ) &&
                
                //Check 8: Min/Max V0 Lifetime cut
                ( ( fTreeCascVarV0Lifetime > lCascadeResult->GetCutMinV0Lifetime() ) &&
                 ( fTreeCascVarV0Lifetime < lCascadeResult->GetCutMaxV0Lifetime() ||
                  lCascadeResult->GetCutMaxV0Lifetime() > 1e+3 ) ) &&
                
                //Check 9: kITSrefit track selection if requested
                (
                 ( (fTreeCascVarPosTrackStatus & AliESDtrack::kITSrefit) &&
                  (fTreeCascVarNegTrackStatus & AliESDtrack::kITSrefit) &&
                  (fTreeCascVarBachTrackStatus & AliESDtrack::kITSrefit)
                  )
                 ||
                 !lCascadeResult->GetCutUseITSRefitTracks()
                 ) &&
                
                //Check 10: Max Chi2/Clusters if not absurd
                ( lCascadeResult->GetCutMaxChi2PerCluster()>1e+3 ||
                 (fTreeCascVarMaxChi2PerCluster < lCascadeResult->GetCutMaxChi2PerCluster())
                 )&&
                
                //Check 11: Min Track Length if positive
                ( lCascadeResult->GetCutMinTrackLength()<0 || //this is a bit paranoid...
                 fTreeCascVarMinTrackLength > lCascadeResult->GetCutMinTrackLength()
                 )&&
                
                //Check 12: Check if special V0 CosPA cut used
                //either don't use the cut at all, or make sure it's above threshold
                ( lCascadeResult->GetCutUse276TeVV0CosPA()==kFALSE ||
                 fTreeCascVarV0CosPointingAngle>l276TeVV0CosPA
                 )
                
                ) {
                //This satisfies all my conditionals! Fill histogram
                histoout_V0A -> Fill ( fCentrality_V0A, fTreeCascVarPt, lMass );
                histoout_V0C -> Fill ( fCentrality_V0C, fTreeCascVarPt, lMass );
                histoout_V0M -> Fill ( fCentrality_V0M, fTreeCascVarPt, lMass );
                
            }
        }
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        // End Superlight adaptive output mode
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        
    }// end of the Cascade loop (ESD or AOD)
    
    // Post output data.
    PostData(1, fListHist    );
    PostData(2, fListV0      );
    PostData(3, fListCascade );
    if( fkSaveEventTree   ) PostData(4, fTreeEvent   );
    if( fkSaveV0Tree      ) PostData(5, fTreeV0      );
    if( fkSaveCascadeTree ) PostData(6, fTreeCascade );
}

//________________________________________________________________________
void AliAnalysisTaskStrangenessVsMultiplicityRun2pPb::Terminate(Option_t *)
{
    // Draw result to the screen
    // Called once at the end of the query
    
    TList *cRetrievedList = 0x0;
    cRetrievedList = (TList*)GetOutputData(1);
    if(!cRetrievedList) {
        Printf("ERROR - AliAnalysisTaskStrangenessVsMultiplicityRun2pPb : ouput data container list not available\n");
        return;
    }
    
    fHistEventCounter = dynamic_cast<TH1D*> (  cRetrievedList->FindObject("fHistEventCounter")  );
    if (!fHistEventCounter) {
        Printf("ERROR - AliAnalysisTaskStrangenessVsMultiplicityRun2pPb : fHistEventCounter not available");
        return;
    }
    
    TCanvas *canCheck = new TCanvas("AliAnalysisTaskStrangenessVsMultiplicityRun2pPb","V0 Multiplicity",10,10,510,510);
    canCheck->cd(1)->SetLogy();
    
    fHistEventCounter->SetMarkerStyle(22);
    fHistEventCounter->DrawCopy("E");
}

//________________________________________________________________________
Double_t AliAnalysisTaskStrangenessVsMultiplicityRun2pPb::MyRapidity(Double_t rE, Double_t rPz) const
{
    // Local calculation for rapidity
    Double_t ReturnValue = -100;
    if( (rE-rPz+1.e-13) != 0 && (rE+rPz) != 0 ) {
        ReturnValue =  0.5*TMath::Log((rE+rPz)/(rE-rPz+1.e-13));
    }
    return ReturnValue;
}

//________________________________________________________________________
void AliAnalysisTaskStrangenessVsMultiplicityRun2pPb::AddConfiguration( AliV0Result *lV0Result )
{
    if (!fListV0){
        Printf("fListV0 does not exist. Creating...");
        fListV0 = new TList();
        fListV0->SetOwner();
        
    }
    fListV0->Add(lV0Result);
}

//________________________________________________________________________
void AliAnalysisTaskStrangenessVsMultiplicityRun2pPb::AddConfiguration( AliCascadeResult *lCascadeResult )
{
    if (!fListCascade){
        Printf("fListCascade does not exist. Creating...");
        fListCascade = new TList();
        fListCascade->SetOwner();
        
    }
    fListCascade->Add(lCascadeResult);
}

//________________________________________________________________________
void AliAnalysisTaskStrangenessVsMultiplicityRun2pPb::SetupStandardVertexing()
//Meant to store standard re-vertexing configuration
{
    //Tell the task to re-run vertexers
    SetRunVertexers(kTRUE);
    SetDoV0Refit(kTRUE);
    
    //V0-Related topological selections
    SetV0VertexerDCAFirstToPV(0.05);
    SetV0VertexerDCASecondtoPV(0.05);
    SetV0VertexerDCAV0Daughters(1.20);
    SetV0VertexerCosinePA(0.98);
    SetV0VertexerMinRadius(0.9);
    SetV0VertexerMaxRadius(200);
    
    //Cascade-Related topological selections
    SetCascVertexerMinV0ImpactParameter(0.05);
    SetCascVertexerV0MassWindow(0.006);
    SetCascVertexerDCABachToPV(0.02);
    SetCascVertexerDCACascadeDaughters(1.2);
    SetCascVertexerCascadeMinRadius(.8);
    SetCascVertexerCascadeCosinePA(.98);
}

//________________________________________________________________________
void AliAnalysisTaskStrangenessVsMultiplicityRun2pPb::SetupLooseVertexing()
//Meant to store standard re-vertexing configuration
{
    //Tell the task to re-run vertexers
    SetRunVertexers(kTRUE);
    SetDoV0Refit(kTRUE);
    
    //V0-Related topological selections
    SetV0VertexerDCAFirstToPV(0.1);
    SetV0VertexerDCASecondtoPV(0.1);
    SetV0VertexerDCAV0Daughters(1.40);
    SetV0VertexerCosinePA(0.95);
    SetV0VertexerMinRadius(0.9);
    SetV0VertexerMaxRadius(200);
    
    //Cascade-Related topological selections
    SetCascVertexerMinV0ImpactParameter(0.05);
    SetCascVertexerV0MassWindow(0.006);
    SetCascVertexerDCABachToPV(0.02);
    SetCascVertexerDCACascadeDaughters(1.4);
    SetCascVertexerCascadeMinRadius(.5);
    SetCascVertexerCascadeCosinePA(.95);
}

//________________________________________________________________________
void AliAnalysisTaskStrangenessVsMultiplicityRun2pPb::AddTopologicalQAV0(Int_t lRecNumberOfSteps)
//Add all configurations to do QA of topological variables for the V0 analysis
{
    // STEP 1: Decide on binning (needed to improve on memory consumption)
    
    // pT binning
    Double_t lPtbinlimits[] ={0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
        2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.5, 8.0, 10, 12, 15};
    Long_t lPtbinnumb = sizeof(lPtbinlimits)/sizeof(Double_t) - 1;
    
    Double_t lPtbinlimitsCascade[] ={0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
        2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.5, 8.0, 10, 12, 14, 17, 20};
    Long_t lPtbinnumbCascade = sizeof(lPtbinlimitsCascade)/sizeof(Double_t) - 1;
    
    // centrality binning
    Double_t lCentbinlimits[] = {0, 10};
    Long_t lCentbinnumb = sizeof(lCentbinlimits)/sizeof(Double_t) - 1;
    
    // TStrings for output names
    TString lParticleName[] = {"K0Short", "Lambda",  "AntiLambda"};
    
    //STEP 3: Creation of output objects
    
    //Map to mass hypothesis
    AliV0Result::EMassHypo lMassHypoV0[3];
    lMassHypoV0[0] = AliV0Result::kK0Short;
    lMassHypoV0[1] = AliV0Result::kLambda;
    lMassHypoV0[2] = AliV0Result::kAntiLambda;
    
    Float_t lLifetimeCut[3];
    lLifetimeCut[0] = 20.0;
    lLifetimeCut[1] = 30.0;
    lLifetimeCut[2] = 30.0;
    
    Float_t lMass[3];
    lMass[0] = 0.497;
    lMass[1] = 1.116;
    lMass[2] = 1.116;
    
    Float_t lMWindow[3];
    lMWindow[0] = 0.075;
    lMWindow[1] = 0.050;
    lMWindow[2] = 0.050;
    
    //Array of results
    AliV0Result *lV0Result[5000];
    Long_t lNV0 = 0;
    
    //Central results: Stored in indices 0, 1, 2 (careful!)
    for(Int_t i = 0 ; i < 3 ; i ++){
        //Central result, customized binning: the one to use, usually
        lV0Result[lNV0] = new AliV0Result( Form("%s_Central",lParticleName[i].Data() ),lMassHypoV0[i],"",lCentbinnumb,lCentbinlimits, lPtbinnumb,lPtbinlimits, 100,lMass[i]-lMWindow[i],lMass[i]+lMWindow[i]);
        //if ( i>0 ) not neeed for real data
        //    lV0Result[lNV0]->InitializeFeeddownMatrix( lPtbinnumb, lPtbinlimits, lPtbinnumbCascade, lPtbinlimitsCascade, lCentbinnumb, lCentbinlimits );
        
        //Setters for V0 Cuts
        lV0Result[lNV0]->SetCutDCANegToPV            ( 0.05 ) ;
        lV0Result[lNV0]->SetCutDCAPosToPV            ( 0.05 ) ;
        lV0Result[lNV0]->SetCutDCAV0Daughters        ( 1.2 ) ;
        lV0Result[lNV0]->SetCutV0CosPA               ( 0.98 ) ;
        lV0Result[lNV0]->SetCutV0Radius              ( 0.9 ) ;
        
        //Miscellaneous
        lV0Result[lNV0]->SetCutProperLifetime        ( lLifetimeCut[i] ) ;
        lV0Result[lNV0]->SetCutLeastNumberOfCrossedRows ( 70 ) ;
        lV0Result[lNV0]->SetCutLeastNumberOfCrossedRowsOverFindable               ( 0.8 ) ;
        lV0Result[lNV0]->SetCutTPCdEdx               ( 4 ) ;
        
        //Add result to pool
        lNV0++;
    }
    
    //Will now proceed to sweep individual variables
    //Number of steps used for the variable sweeps
    const Int_t lNumberOfSteps = lRecNumberOfSteps;
    
    //________________________________________________________
    // Variable 1: DCA Neg to PV
    Float_t lMaxDCANegToPV = 20.00;
    
    for(Int_t i = 0 ; i < 3 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lV0Result[lNV0] = new AliV0Result( lV0Result[i], Form("%s_%s_%i",lParticleName[i].Data(),"DCANegToPVSweep",icut) );
            //Add result to pool
            Float_t lThisCut = ((Float_t)icut+1)*lMaxDCANegToPV / ((Float_t) lNumberOfSteps) ;
            lV0Result[lNV0] -> SetCutDCANegToPV ( lThisCut );
            lNV0++;
        }
    }
    //________________________________________________________
    // Variable 2: DCA Pos to PV
    Float_t lMaxDCAPosToPV = 20.00;
    
    for(Int_t i = 0 ; i < 3 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lV0Result[lNV0] = new AliV0Result( lV0Result[i], Form("%s_%s_%i",lParticleName[i].Data(),"DCAPosToPVSweep",icut) );
            //Add result to pool
            Float_t lThisCut = ((Float_t)icut+1)*lMaxDCAPosToPV / ((Float_t) lNumberOfSteps) ;
            lV0Result[lNV0] -> SetCutDCAPosToPV ( lThisCut );
            lNV0++;
        }
    }
    //________________________________________________________
    // Variable 3: DCA V0 Daughters
    Float_t lMaxDCAV0Daughters = 1.20;
    
    for(Int_t i = 0 ; i < 3 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lV0Result[lNV0] = new AliV0Result( lV0Result[i], Form("%s_%s_%i",lParticleName[i].Data(),"DCAV0DaughtersSweep",icut) );
            //Add result to pool
            Float_t lThisCut = ((Float_t)icut+1)*lMaxDCAV0Daughters / ((Float_t) lNumberOfSteps) ;
            lV0Result[lNV0] -> SetCutDCAV0Daughters ( lThisCut );
            lNV0++;
        }
    }
    //________________________________________________________
    // Variable 4: V0 CosPA
    Float_t lMinV0CosPA = 0.98;
    Float_t lMaxV0CosPA = 1.00;
    Double_t lV0CosPAVals[lNumberOfSteps];
    Double_t lMinV0PA = 0.0;
    Double_t lMaxV0PA = TMath::ACos(lMinV0CosPA);
    Double_t lDeltaV0PA = lMaxV0PA / ((Double_t)(lNumberOfSteps));
    for(Int_t iStep = 0; iStep<lNumberOfSteps; iStep++){
        lV0CosPAVals[iStep] = TMath::Cos( ((Float_t)(iStep+1))*lDeltaV0PA );
    }
    for(Int_t i = 0 ; i < 3 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lV0Result[lNV0] = new AliV0Result( lV0Result[i], Form("%s_%s_%i",lParticleName[i].Data(),"V0CosPASweep",icut) );
            //Add result to pool
            lV0Result[lNV0] -> SetCutV0CosPA ( lV0CosPAVals[icut] );
            lNV0++;
        }
    }
    //________________________________________________________
    // Variable 5: V0 Radius
    Float_t lMinV0Radius = 2.0;
    Float_t lMaxV0Radius = 20.00;
    for(Int_t i = 0 ; i < 3 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lV0Result[lNV0] = new AliV0Result( lV0Result[i], Form("%s_%s_%i",lParticleName[i].Data(),"V0RadiusSweep",icut) );
            //Add result to pool
            Float_t lThisCut = lMinV0Radius + (lMaxV0Radius-lMinV0Radius)*(((Float_t)icut)+1)/((Float_t)lNumberOfSteps);
            lV0Result[lNV0] -> SetCutV0Radius ( lThisCut );
            lNV0++;
        }
    }
    for (Int_t iconf = 0; iconf<lNV0; iconf++)
        AddConfiguration(lV0Result[iconf]);
    
    cout<<"Added "<<lNV0<<" V0 configurations to output."<<endl;
}

//________________________________________________________________________
void AliAnalysisTaskStrangenessVsMultiplicityRun2pPb::AddTopologicalQACascade(Int_t lRecNumberOfSteps)
//Add all configurations to do QA of topological variables for the V0 analysis
{
    // STEP 1: Decide on binning (needed to improve on memory consumption)
    
    // pT binning
    Double_t lPtbinlimits[] = {0.4, 0.5, 0.6,
        0.7,0.8,.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,
        2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,
        4.4,4.6,4.8,5.0,5.5,6.0,6.5,7.0,8.0,9.0,10.,11.,12.};
    //Double_t lPtbinlimits[] = {0.2,0.3, 0.4, 0.5, 0.6,
    //    0.7,0.8,.9,1.0,1.2, 1.4, 1.6, 1.8 ,2.0,
    //    2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0,
    //    4.4,4.8,5.0,6.0,7.0,8.0,9.0,10.,11.,12.};
    
    //Double_t lPtbinlimits[] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2,
    //3.6, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.5, 8.5, 10, 12};
    
    Long_t lPtbinnumb = sizeof(lPtbinlimits)/sizeof(Double_t) - 1;
    
    // centrality binning
    Double_t lCentbinlimits[] = {0, 10}; //optimize in 0-10%
    Long_t lCentbinnumb = sizeof(lCentbinlimits)/sizeof(Double_t) - 1;
    
    //Just a counter and one array, please. Nothing else needed
    AliCascadeResult *lCascadeResult[5000];
    Long_t lN = 0;
    
    //Map to mass hypothesis
    AliCascadeResult::EMassHypo lMassHypo[4];
    lMassHypo[0] = AliCascadeResult::kXiMinus;
    lMassHypo[1] = AliCascadeResult::kXiPlus;
    lMassHypo[2] = AliCascadeResult::kOmegaMinus;
    lMassHypo[3] = AliCascadeResult::kOmegaPlus;
    
    Float_t lLifetimeCut[4];
    lLifetimeCut[0] = 15.0;
    lLifetimeCut[1] = 15.0;
    lLifetimeCut[2] = 12.0;
    lLifetimeCut[3] = 12.0;
    
    Float_t lMass[4];
    lMass[0] = 1.322;
    lMass[1] = 1.322;
    lMass[2] = 1.672;
    lMass[3] = 1.672;
    
    TString lParticleName[] = {"XiMinus", "XiPlus",  "OmegaMinus", "OmegaPlus"};
    
    //Number of steps used for the variable sweeps
    const Int_t lNumberOfSteps = lRecNumberOfSteps;
    
    //Central results: Stored in indices 0, 1, 2, 3 (careful!)
    for(Int_t i = 0 ; i < 4 ; i ++){
        //Central result, customized binning: the one to use, usually
        lCascadeResult[lN] = new AliCascadeResult( Form("%s_VertexerLevel",lParticleName[i].Data() ),lMassHypo[i],"",lCentbinnumb,lCentbinlimits, lPtbinnumb,lPtbinlimits,100,lMass[i]-0.050,lMass[i]+0.050);
        
        //Default cuts: use vertexer level ones
        //Setters for V0 Cuts
        lCascadeResult[lN]->SetCutDCANegToPV            ( 0.2 ) ;
        lCascadeResult[lN]->SetCutDCAPosToPV            ( 0.2 ) ;
        lCascadeResult[lN]->SetCutDCAV0Daughters        (  1. ) ;
        lCascadeResult[lN]->SetCutV0CosPA               ( 0.95 ) ; //+variable
        lCascadeResult[lN]->SetCutVarV0CosPA            (TMath::Exp(10.853),
                                                         -25.0322,
                                                         TMath::Exp(-0.843948),
                                                         -0.890794,
                                                         0.057553);
        lCascadeResult[lN]->SetCutV0Radius              (  3 ) ;
        //Setters for Cascade Cuts
        lCascadeResult[lN]->SetCutDCAV0ToPV             ( 0.1 ) ;
        lCascadeResult[lN]->SetCutV0Mass                ( 0.006 ) ;
        lCascadeResult[lN]->SetCutDCABachToPV           ( 0.1 ) ;
        lCascadeResult[lN]->SetCutDCACascDaughters      ( 1.0) ;
        lCascadeResult[lN]->SetCutCascRadius            ( 1.2 ) ;
        if(i==2||i==3)
            lCascadeResult[lN]->SetCutCascRadius            ( 1.0 ) ; //omega case
        lCascadeResult[lN]->SetCutCascCosPA             ( 0.95 ) ; //+variable
        lCascadeResult[lN]->SetCutVarCascCosPA          (TMath::Exp(4.86664),
                                                         -10.786,
                                                         TMath::Exp(-1.33411),
                                                         -0.729825,
                                                         0.0695724);
        //Miscellaneous
        lCascadeResult[lN]->SetCutProperLifetime        ( lLifetimeCut[i] ) ;
        lCascadeResult[lN]->SetCutLeastNumberOfClusters ( 70.0 ) ;
        lCascadeResult[lN]->SetCutTPCdEdx               ( 4.0 ) ;
        lCascadeResult[lN]->SetCutXiRejection           ( 0.008 ) ;
        lCascadeResult[lN]->SetCutBachBaryonCosPA       ( TMath::Cos(0.04) ) ; //+variable
        lCascadeResult[lN]->SetCutVarBBCosPA            (TMath::Exp(-2.29048),
                                                         -20.2016,
                                                         TMath::Exp(-2.9581),
                                                         -0.649153,
                                                         0.00526455);
        //Add result to pool
        lN++;
    }
    
    //Will now proceed to sweep individual variables
    
    //________________________________________________________
    // Variable 1: DCA Neg to PV
    Float_t lMaxDCANegToPV = 1.5;
    
    for(Int_t i = 0 ; i < 4 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%i",lParticleName[i].Data(),"DCANegToPVSweep",icut) );
            //Add result to pool
            Float_t lThisCut = ((Float_t)icut+1)*lMaxDCANegToPV / ((Float_t) lNumberOfSteps) ;
            lCascadeResult[lN] -> SetCutDCANegToPV ( lThisCut );
            lN++;
        }
    }
    //________________________________________________________
    // Variable 2: DCA Pos to PV
    Float_t lMaxDCAPosToPV = 1.5;
    
    for(Int_t i = 0 ; i < 4 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%i",lParticleName[i].Data(),"DCAPosToPVSweep",icut) );
            //Add result to pool
            Float_t lThisCut = ((Float_t)icut+1)*lMaxDCAPosToPV / ((Float_t) lNumberOfSteps) ;
            lCascadeResult[lN] -> SetCutDCAPosToPV ( lThisCut );
            lN++;
        }
    }
    //________________________________________________________
    // Variable 3: DCA V0 Daughters
    Float_t lMaxDCAV0Daughters = 1.40;
    
    for(Int_t i = 0 ; i < 4 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%i",lParticleName[i].Data(),"DCAV0DaughtersSweep",icut) );
            //Add result to pool
            Float_t lThisCut = ((Float_t)icut+1)*lMaxDCAV0Daughters / ((Float_t) lNumberOfSteps) ;
            lCascadeResult[lN] -> SetCutDCAV0Daughters ( lThisCut );
            lN++;
        }
    }
    //________________________________________________________
    // Variable 4: V0 CosPA
    Float_t lMinV0CosPA = 0.95;
    Float_t lMaxV0CosPA = 1.00;
    Double_t lV0CosPAVals[lNumberOfSteps];
    Double_t lMinV0PA = 0.0;
    Double_t lMaxV0PA = TMath::ACos(lMinV0CosPA);
    Double_t lDeltaV0PA = lMaxV0PA / ((Double_t)(lNumberOfSteps));
    for(Int_t iStep = 0; iStep<lNumberOfSteps; iStep++){
        lV0CosPAVals[iStep] = TMath::Cos( ((Float_t)(iStep+1))*lDeltaV0PA );
    }
    for(Int_t i = 0 ; i < 4 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%i",lParticleName[i].Data(),"V0CosPASweep",icut) );
            //Add result to pool
            lCascadeResult[lN] -> SetCutUseVarV0CosPA( kFALSE );
            lCascadeResult[lN] -> SetCutV0CosPA ( lV0CosPAVals[icut] );
            lN++;
        }
    }
    //________________________________________________________
    // Variable 5: V0 Radius
    Float_t lMinV0Radius = 0.0;
    Float_t lMaxV0Radius = 20.00;
    for(Int_t i = 0 ; i < 4 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%i",lParticleName[i].Data(),"V0RadiusSweep",icut) );
            //Add result to pool
            Float_t lThisCut = lMinV0Radius + (lMaxV0Radius-lMinV0Radius)*(((Float_t)icut)+1)/((Float_t)lNumberOfSteps);
            lCascadeResult[lN] -> SetCutV0Radius ( lThisCut );
            lN++;
        }
    }
    //________________________________________________________
    // Variable 6:
    Float_t lMaxDCAV0ToPV = 0.5;
    for(Int_t i = 0 ; i < 4 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%i",lParticleName[i].Data(),"DCAV0ToPVSweep",icut) );
            //Add result to pool
            Float_t lThisCut = ((Float_t)icut+1)*lMaxDCAV0ToPV / ((Float_t) lNumberOfSteps) ;
            lCascadeResult[lN] -> SetCutDCAV0ToPV ( lThisCut );
            lN++;
        }
    }
    //________________________________________________________
    // Variable 7: DCA Bach To PV
    Float_t lMaxDCABachToPV = 0.5;
    for(Int_t i = 0 ; i < 4 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%i",lParticleName[i].Data(),"DCABachToPVSweep",icut) );
            //Add result to pool
            Float_t lThisCut = ((Float_t)icut+1)*lMaxDCABachToPV / ((Float_t) lNumberOfSteps) ;
            lCascadeResult[lN] -> SetCutDCABachToPV ( lThisCut );
            lN++;
        }
    }
    //________________________________________________________
    // Variable 8: DCA Casc Daughters
    Float_t lMaxDCACascDaughters = 1.40;
    
    for(Int_t i = 0 ; i < 4 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%i",lParticleName[i].Data(),"DCACascDaughtersSweep",icut) );
            //Add result to pool
            Float_t lThisCut = ((Float_t)icut+1)*lMaxDCACascDaughters / ((Float_t) lNumberOfSteps) ;
            lCascadeResult[lN] -> SetCutDCACascDaughters ( lThisCut );
            lN++;
        }
    }
    //________________________________________________________
    // Variable 9: Cascade Radius
    Float_t lMinCascRadius = 0.5;
    Float_t lMaxCascRadius = 7.0;
    for(Int_t i = 0 ; i < 4 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%i",lParticleName[i].Data(),"CascRadiusSweep",icut) );
            //Add result to pool
            Float_t lThisCut = lMinCascRadius + (lMaxCascRadius-lMinCascRadius)*(((Float_t)icut)+1)/((Float_t)lNumberOfSteps);
            lCascadeResult[lN] -> SetCutCascRadius ( lThisCut );
            lN++;
        }
    }
    //________________________________________________________
    // Variable 10: Cascade CosPA
    Float_t lMinCascCosPA = 0.95;
    Float_t lMaxCascCosPA = 1.00;
    Double_t lCascCosPAVals[lNumberOfSteps];
    Double_t lMinCascPA = 0.0;
    Double_t lMaxCascPA = TMath::ACos(lMinCascCosPA);
    Double_t lDeltaCascPA = lMaxCascPA / ((Double_t)(lNumberOfSteps));
    for(Int_t iStep = 0; iStep<lNumberOfSteps; iStep++){
        lCascCosPAVals[iStep] = TMath::Cos( ((Float_t)(iStep+1))*lDeltaCascPA );
    }
    for(Int_t i = 0 ; i < 4 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%i",lParticleName[i].Data(),"CascCosPASweep",icut) );
            //Add result to pool
            lCascadeResult[lN] -> SetCutUseVarCascCosPA( kFALSE );
            lCascadeResult[lN] -> SetCutCascCosPA ( lCascCosPAVals[icut] );
            lN++;
        }
    }
    //________________________________________________________
    // Variable 11: Bach-Baryon CosPA
    Float_t lMinBBCosPA = TMath::Cos(0.1);
    Float_t lMaxBBCosPA = 1.000;
    Double_t lBBCosPAVals[lNumberOfSteps];
    Double_t lMinBBPA = 0.0;
    Double_t lMaxBBPA = TMath::ACos(lMinBBCosPA);
    Double_t lDeltaBBPA = lMaxBBPA / ((Double_t)(lNumberOfSteps));
    for(Int_t iStep = 0; iStep<lNumberOfSteps; iStep++){
        lBBCosPAVals[iStep] = TMath::Cos( ((Float_t)(iStep+1))*lDeltaBBPA );
    }
    for(Int_t i = 0 ; i < 4 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%i",lParticleName[i].Data(),"BBCosPASweep",icut) );
            //Add result to pool
            lCascadeResult[lN] -> SetCutUseVarBBCosPA( kFALSE );
            lCascadeResult[lN] -> SetCutBachBaryonCosPA ( lBBCosPAVals[icut] );
            lN++;
        }
    }
    //________________________________________________________
    // Variable 12: Cascade Lifetime Sweep
    
    Int_t lLifetimeSteps = 15;
    for(Int_t i = 0 ; i < 4 ; i ++){
        Float_t lMinLifetime = 5.00;
        Float_t lMaxLifetime = 20.00;
        for(Int_t icut = 0; icut<lLifetimeSteps; icut++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%i",lParticleName[i].Data(),"CascLifetimeSweep",icut) );
            Float_t lThisCut = lMinLifetime + (lMaxLifetime-lMinLifetime)*(((Float_t)icut)+1)/((Float_t)lLifetimeSteps);
            //Add result to pool
            lCascadeResult[lN] -> SetCutProperLifetime ( lThisCut );
            lN++;
        }
    }
    //________________________________________________________
    // Variable 13: V0 Lifetime Sweep
    Float_t lMinV0Lifetime = 8.00;
    Float_t lMaxV0Lifetime = 40.00;
    Int_t lV0LifetimeSteps = 32;
    for(Int_t i = 0 ; i < 4 ; i ++){
        for(Int_t icut = 0; icut<lV0LifetimeSteps; icut++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%i",lParticleName[i].Data(),"MaxV0LifetimeSweep",icut) );
            Float_t lThisCut = lMinV0Lifetime + (lMaxV0Lifetime-lMinV0Lifetime)*(((Float_t)icut)+1)/((Float_t)lV0LifetimeSteps);
            //Add result to pool
            lCascadeResult[lN] -> SetCutMaxV0Lifetime ( lThisCut );
            lN++;
        }
    }
    
    for (Int_t iconf = 0; iconf<lN; iconf++)
        AddConfiguration(lCascadeResult[iconf]);
    
    cout<<"Added "<<lN<<" Cascade configurations to output."<<endl;
}

//________________________________________________________________________
void AliAnalysisTaskStrangenessVsMultiplicityRun2pPb::AddStandardV0Configuration()
//Meant to add some standard V0 analysis Configuration + its corresponding systematics
{
    // STEP 1: Decide on binning (needed to improve on memory consumption)
    
    // pT binning
    Double_t lPtbinlimits[] ={0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
        2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.5, 8.0, 10, 12, 15};
    Long_t lPtbinnumb = sizeof(lPtbinlimits)/sizeof(Double_t) - 1;
    
    // centrality binning
    Double_t lCentbinlimits[] = {0, 1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 85, 90};
    Long_t lCentbinnumb = sizeof(lCentbinlimits)/sizeof(Double_t) - 1;
    
    // TStrings for output names
    TString lParticleName[] = {"K0Short", "Lambda",  "AntiLambda"};
    TString lConfName[]     = {"Loose",   "Central", "Tight"     };
    TString lCutName[]      = {"DCANegToPV","DCAPosToPV","DCAV0Daughters","V0CosPA","V0Radius",
        "ProperLifetime","LeastNbrCrs","LeastNbrCrsOvFind","TPCdEdx"};
    
    // STEP 2: Decide on a set of selections
    
    //1st index: Particle Species
    //2nd index: Loose / Central / Tight
    //3rd index: Number of selection (as ordered above)
    Double_t lcutsV0[3][3][9];
    
    //N.B.: These are mostly symmetric, except for the proper lifetime, which is different
    //      for the two particle species. Another potential improvement could be asymmetric
    //      DCA selections for the Neg / Pos tracks for the (anti)Lambda decay, as decay
    //      kinematics would prefer having the pion daughter with a larger DCA.
    
    /***
     //Information from Michal, 09th April
     === Lambda ===
     --- 0-10% ---
     Signal loss: 5%
     Parameters from real data:
     par 0   0.18945
     par 1  -0.57882
     par 2   0.01302
     Signal loss: 5%
     Parameters from MC:
     par 0   0.26081
     par 1  -1.16285
     par 2   0.02692
     
     === AntiLambda ===
     --- 0-10% ---
     Signal loss: 5%
     Parameters from real data:
     par 0   0.21861
     par 1  -0.67273
     par 2   0.01200
     Signal loss: 5%
     Parameters from MC:
     par 0   0.24144
     par 1  -1.04444
     par 2   0.02684
     
     === K0Short ===
     --- 0-10% ---
     Signal loss: 5%
     Parameters from real data:
     par 0   0.21320
     par 1  -0.91380
     par 2   0.02483
     Signal loss: 5%
     Parameters from MC:
     par 0   0.17816
     par 1  -0.79000
     par 2   0.02184
     */
    
    Double_t parExpConst[3] = { 0.26081, 0.24144, 0.17816 };
    Double_t parExpSlope[3] = { -1.16285, -1.04444, -0.79000 };
    Double_t parConst[3]    = { 0.02692, 0.02684, 0.02184 };
    
    //================================================================================
    // K0SHORT SELECTIONS
    //--------------------------------------------------------------------------------
    //                  LOOSE                        CENTRAL                   TIGHT
    lcutsV0[0][0][ 0] = 0.06;    lcutsV0[0][1][ 0] =   0.1; lcutsV0[0][2][0] = 0.17; //DCANegToPV
    lcutsV0[0][0][ 1] = 0.06;    lcutsV0[0][1][ 1] =   0.1; lcutsV0[0][2][1] = 0.17; //DCAPosToPV
    lcutsV0[0][0][ 2] = 0.95;    lcutsV0[0][1][ 2] =   0.8; lcutsV0[0][2][2] =  0.7; //DCAV0Daughters
    lcutsV0[0][0][ 3] = 0.95;    lcutsV0[0][1][ 3] =  0.95; lcutsV0[0][2][3] = 0.95; //V0CosPA
    lcutsV0[0][0][ 4] = 4.50;    lcutsV0[0][1][ 4] =  5.00; lcutsV0[0][2][4] = 5.50; //V0Radius
    lcutsV0[0][0][ 5] =   12;    lcutsV0[0][1][ 5] =    10; lcutsV0[0][2][5] =    8; //Proper Lifetime (in cm)
    lcutsV0[0][0][ 6] =   70;    lcutsV0[0][1][ 6] =    70; lcutsV0[0][2][6] =   80; //Least Nbr Crossed Rows
    lcutsV0[0][0][ 7] =  0.7;    lcutsV0[0][1][ 7] =   0.8; lcutsV0[0][2][7] = 0.85; //Least Ratio CrdRows/Findable
    lcutsV0[0][0][ 8] =  4.0;    lcutsV0[0][1][ 8] =   3.0; lcutsV0[0][2][8] =  2.5; //TPC dE/dx
    //================================================================================
    
    //================================================================================
    // LAMBDA SELECTIONS
    //--------------------------------------------------------------------------------
    //                  LOOSE                        CENTRAL                   TIGHT
    lcutsV0[1][0][ 0] =  0.1;    lcutsV0[1][1][ 0] =   0.2; lcutsV0[1][2][0] = 0.30; //DCANegToPV
    lcutsV0[1][0][ 1] = 0.08;    lcutsV0[1][1][ 1] =   0.1; lcutsV0[1][2][1] = 0.13; //DCAPosToPV
    lcutsV0[1][0][ 2] =  1.0;    lcutsV0[1][1][ 2] =   0.8; lcutsV0[1][2][2] = 0.65; //DCAV0Daughters
    lcutsV0[1][0][ 3] = 0.95;    lcutsV0[1][1][ 3] =  0.95; lcutsV0[1][2][3] = 0.95; //V0CosPA
    lcutsV0[1][0][ 4] = 4.00;    lcutsV0[1][1][ 4] =  5.00; lcutsV0[1][2][4] = 6.00; //V0Radius
    lcutsV0[1][0][ 5] =   24;    lcutsV0[1][1][ 5] =    20; lcutsV0[1][2][5] =   17; //Proper Lifetime (in cm)
    lcutsV0[1][0][ 6] =   70;    lcutsV0[1][1][ 6] =    70; lcutsV0[1][2][6] =   80; //Least Nbr Crossed Rows
    lcutsV0[1][0][ 7] =  0.7;    lcutsV0[1][1][ 7] =   0.8; lcutsV0[1][2][7] = 0.85; //Least Ratio CrdRows/Findable
    lcutsV0[1][0][ 8] =  4.0;    lcutsV0[1][1][ 8] =   3.0; lcutsV0[1][2][8] =  2.5; //TPC dE/dx
    //================================================================================
    
    //================================================================================
    // ANTILAMBDA SELECTIONS
    //--------------------------------------------------------------------------------
    //                  LOOSE                        CENTRAL                   TIGHT
    lcutsV0[2][0][ 0] = 0.08;    lcutsV0[2][1][ 0] =   0.1; lcutsV0[2][2][0] = 0.13; //DCANegToPV
    lcutsV0[2][0][ 1] =  0.1;    lcutsV0[2][1][ 1] =   0.2; lcutsV0[2][2][1] = 0.30; //DCAPosToPV
    lcutsV0[2][0][ 2] =  1.0;    lcutsV0[2][1][ 2] =   0.8; lcutsV0[2][2][2] = 0.65; //DCAV0Daughters
    lcutsV0[2][0][ 3] = 0.95;    lcutsV0[2][1][ 3] =  0.95; lcutsV0[2][2][3] = 0.95; //V0CosPA
    lcutsV0[2][0][ 4] = 4.00;    lcutsV0[2][1][ 4] =  5.00; lcutsV0[2][2][4] = 6.00; //V0Radius
    lcutsV0[2][0][ 5] =   24;    lcutsV0[2][1][ 5] =    20; lcutsV0[2][2][5] =   17; //Proper Lifetime (in cm)
    lcutsV0[2][0][ 6] =   70;    lcutsV0[2][1][ 6] =    70; lcutsV0[2][2][6] =   80; //Least Nbr Crossed Rows
    lcutsV0[2][0][ 7] =  0.7;    lcutsV0[2][1][ 7] =   0.8; lcutsV0[2][2][7] = 0.85; //Least Ratio CrdRows/Findable
    lcutsV0[2][0][ 8] =  4.0;    lcutsV0[2][1][ 8] =   3.0; lcutsV0[2][2][8] =  2.5; //TPC dE/dx
    //================================================================================
    
    
    //STEP 3: Creation of output objects
    
    //Map to mass hypothesis
    AliV0Result::EMassHypo lMassHypoV0[3];
    lMassHypoV0[0] = AliV0Result::kK0Short;
    lMassHypoV0[1] = AliV0Result::kLambda;
    lMassHypoV0[2] = AliV0Result::kAntiLambda;
    
    //Array of results
    AliV0Result *lV0Result[500];
    Long_t lNV0 = 0;
    
    //Central results: Stored in indices 0, 1, 2 (careful!)
    for(Int_t i = 0 ; i < 3 ; i ++){
        //Central result, customized binning: the one to use, usually
        lV0Result[lNV0] = new AliV0Result( Form("%s_Central",lParticleName[i].Data() ),lMassHypoV0[i],"",lCentbinnumb,lCentbinlimits, lPtbinnumb,lPtbinlimits);
        
        //Setters for V0 Cuts
        lV0Result[lNV0]->SetCutDCANegToPV            ( lcutsV0[i][1][ 0] ) ;
        lV0Result[lNV0]->SetCutDCAPosToPV            ( lcutsV0[i][1][ 1] ) ;
        lV0Result[lNV0]->SetCutDCAV0Daughters        ( lcutsV0[i][1][ 2] ) ;
        lV0Result[lNV0]->SetCutV0CosPA               ( lcutsV0[i][1][ 3] ) ;
        //Set Variable cut
        lV0Result[lNV0]->SetCutVarV0CosPA               ( parExpConst[i], parExpSlope[i], 0.0, 1.0, parConst[i] ) ;
        
        lV0Result[lNV0]->SetCutV0Radius              ( lcutsV0[i][1][ 4] ) ;
        
        //Miscellaneous
        lV0Result[lNV0]->SetCutProperLifetime        ( lcutsV0[i][1][ 5] ) ;
        lV0Result[lNV0]->SetCutLeastNumberOfCrossedRows ( lcutsV0[i][1][ 6] ) ;
        lV0Result[lNV0]->SetCutLeastNumberOfCrossedRowsOverFindable               ( lcutsV0[i][1][ 7] ) ;
        lV0Result[lNV0]->SetCutTPCdEdx               ( lcutsV0[i][1][ 8] ) ;
        
        //Add result to pool
        lNV0++;
    }
    
    //Central full results: indices 2, 3, 4
    for(Int_t i = 0 ; i < 3 ; i ++){
        //Central Result, Full: No rebinning. Will use a significant amount of memory,
        //not to be replicated several times for systematics!
        lV0Result[lNV0] = new AliV0Result( Form("%s_Central_Full",lParticleName[i].Data() ),lMassHypoV0[i]);
        
        //Setters for V0 Cuts
        lV0Result[lNV0]->SetCutDCANegToPV            ( lcutsV0[i][1][ 0] ) ;
        lV0Result[lNV0]->SetCutDCAPosToPV            ( lcutsV0[i][1][ 1] ) ;
        lV0Result[lNV0]->SetCutDCAV0Daughters        ( lcutsV0[i][1][ 2] ) ;
        lV0Result[lNV0]->SetCutV0CosPA               ( lcutsV0[i][1][ 3] ) ;
        //Set Variable cut
        lV0Result[lNV0]->SetCutVarV0CosPA               ( parExpConst[i], parExpSlope[i], 0.0, 1.0, parConst[i] ) ;
        lV0Result[lNV0]->SetCutV0Radius              ( lcutsV0[i][1][ 4] ) ;
        
        //Miscellaneous
        lV0Result[lNV0]->SetCutProperLifetime        ( lcutsV0[i][1][ 5] ) ;
        lV0Result[lNV0]->SetCutLeastNumberOfCrossedRows ( lcutsV0[i][1][ 6] ) ;
        lV0Result[lNV0]->SetCutLeastNumberOfCrossedRowsOverFindable               ( lcutsV0[i][1][ 7] ) ;
        lV0Result[lNV0]->SetCutTPCdEdx               ( lcutsV0[i][1][ 8] ) ;
        
        //Add result to pool
        lNV0++;
    }
    
    
    // STEP 4: Creation of objects to be used in systematics
    // Optimized via use of copy constructors
    for(Int_t i = 0 ; i < 3 ; i ++){
        for(Int_t iCut = 0 ; iCut < 9 ; iCut ++){
            
            //LOOSE VARIATIONS
            //Create a new object from default
            lV0Result[lNV0] = new AliV0Result( lV0Result[i], Form("%s_%s_%s",lParticleName[i].Data(),lCutName[iCut].Data(),lConfName[0].Data()) );
            
            if(iCut ==  0 ) lV0Result[lNV0]->SetCutDCANegToPV            ( lcutsV0[i][0][iCut] ) ;
            if(iCut ==  1 ) lV0Result[lNV0]->SetCutDCAPosToPV            ( lcutsV0[i][0][iCut] ) ;
            if(iCut ==  2 ) lV0Result[lNV0]->SetCutDCAV0Daughters        ( lcutsV0[i][0][iCut] ) ;
            if(iCut ==  3 ) lV0Result[lNV0]->SetCutV0CosPA               ( lcutsV0[i][0][iCut] ) ;
            
            
            if(iCut ==  4 ) lV0Result[lNV0]->SetCutV0Radius              ( lcutsV0[i][0][iCut] ) ;
            
            //Miscellaneous
            if(iCut ==  5 ) lV0Result[lNV0]->SetCutProperLifetime        ( lcutsV0[i][0][iCut] ) ;
            if(iCut ==  6 ) lV0Result[lNV0]->SetCutLeastNumberOfCrossedRows ( lcutsV0[i][0][iCut] ) ;
            if(iCut ==  7 ) lV0Result[lNV0]->SetCutLeastNumberOfCrossedRowsOverFindable               ( lcutsV0[i][0][iCut] ) ;
            if(iCut ==  8 ) lV0Result[lNV0]->SetCutTPCdEdx               ( lcutsV0[i][0][iCut] ) ;
            
            //Print this variation, add to pool
            lV0Result[lNV0]->Print();
            lNV0++;
            
            //TIGHT VARIATIONS
            //Create a new object from default
            lV0Result[lNV0] = new AliV0Result( lV0Result[i], Form("%s_%s_%s",lParticleName[i].Data(),lCutName[iCut].Data(),lConfName[2].Data()) );
            
            if(iCut ==  0 ) lV0Result[lNV0]->SetCutDCANegToPV            ( lcutsV0[i][2][iCut] ) ;
            if(iCut ==  1 ) lV0Result[lNV0]->SetCutDCAPosToPV            ( lcutsV0[i][2][iCut] ) ;
            if(iCut ==  2 ) lV0Result[lNV0]->SetCutDCAV0Daughters        ( lcutsV0[i][2][iCut] ) ;
            if(iCut ==  3 ) lV0Result[lNV0]->SetCutV0CosPA               ( lcutsV0[i][2][iCut] ) ;
            if(iCut ==  4 ) lV0Result[lNV0]->SetCutV0Radius              ( lcutsV0[i][2][iCut] ) ;
            
            //Miscellaneous
            if(iCut ==  5 ) lV0Result[lNV0]->SetCutProperLifetime        ( lcutsV0[i][2][iCut] ) ;
            if(iCut ==  6 ) lV0Result[lNV0]->SetCutLeastNumberOfCrossedRows ( lcutsV0[i][2][iCut] ) ;
            if(iCut ==  7 ) lV0Result[lNV0]->SetCutLeastNumberOfCrossedRowsOverFindable               ( lcutsV0[i][2][iCut] ) ;
            if(iCut ==  8 ) lV0Result[lNV0]->SetCutTPCdEdx               ( lcutsV0[i][2][iCut] ) ;
            
            //Print this variation, add to pool
            lV0Result[lNV0]->Print();
            lNV0++;
        }
    }
    for (Int_t iconf = 0; iconf<lNV0; iconf++)
        AddConfiguration(lV0Result[iconf]);
    
    cout<<"Added "<<lNV0<<" V0 configurations to output."<<endl;
}

//________________________________________________________________________
void AliAnalysisTaskStrangenessVsMultiplicityRun2pPb::AddStandardCascadeConfiguration()
//Meant to add some standard cascade analysis Configuration + its corresponding systematics
{
    // STEP 1: Decide on binning (needed to improve on memory consumption)
    
    // pT binning
    Double_t lPtbinlimits[] = {0.4, 0.5, 0.6,
        0.7,0.8,.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,
        2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,
        4.4,4.6,4.8,5.0,5.5,6.0,6.5,7.0,8.0,9.0,10.,11.,12.};
    Long_t lPtbinnumb = sizeof(lPtbinlimits)/sizeof(Double_t) - 1;
    
    // centrality binning
    Double_t lCentbinlimits[] = {0, 1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 85, 90};
    Long_t lCentbinnumb = sizeof(lCentbinlimits)/sizeof(Double_t) - 1;
    
    // TStrings for output names
    TString lParticleName[] = {"XiMinus", "XiPlus",  "OmegaMinus", "OmegaPlus"};
    TString lConfName[]     = {"Loose",   "Central", "Tight"     };
    TString lCutName[]      = {
        "DCANegToPV", //1
        "DCAPosToPV", //2
        "DCAV0Daughters", //3
        "V0Radius", //4
        "DCAV0ToPV", //5
        "V0Mass", //6
        "DCABachToPV", //7
        "DCACascDaughters", //8
        "CascRadius", //9
        "ProperLifetime", //10
        "ProperLifetimeV0", //11
        "LeastNumberOfClusters", //12
        "TPCdEdx", //13
        "Competing" //14
    };
    
    // STEP 2: Decide on a set of selections
    
    //1st index: Particle Species
    //2nd index: Loose / Central / Tight
    //3rd index: Number of selection (as ordered above)
    Double_t lcuts[4][3][14];
    
    //N.B.: These are mostly symmetric, except for the proper lifetime, which is different
    //      for the two particle species. Another potential improvement could be asymmetric
    //      DCA selections for the Neg / Pos tracks for the (anti)Lambda decay, as decay
    //      kinematics would prefer having the pion daughter with a larger DCA.
    
    Int_t lIdx = 0;
    
    //================================================================================
    // XIMINUS SELECTIONS
    lIdx = 0; //Master XiMinus Index
    //--------------------------------------------------------------------------------
    //                  LOOSE                        CENTRAL                   TIGHT
    lcuts[lIdx][0][ 0] = 0.10;    lcuts[lIdx][1][ 0] =  0.20; lcuts[lIdx][2][ 0] =  0.30; //DCANegToPV 1
    lcuts[lIdx][0][ 1] = 0.10;    lcuts[lIdx][1][ 1] =  0.20; lcuts[lIdx][2][ 1] =  0.30; //DCAPostToPV 2
    lcuts[lIdx][0][ 2] =  1.2;    lcuts[lIdx][1][ 2] =   1.0; lcuts[lIdx][2][ 2] =   0.8; //DCAV0Daughters 3
    lcuts[lIdx][0][ 3] = 2.00;    lcuts[lIdx][1][ 3] =  3.00; lcuts[lIdx][2][ 3] =   4.0; //V0Radius 4
    lcuts[lIdx][0][ 4] = 0.05;    lcuts[lIdx][1][ 4] =   0.1; lcuts[lIdx][2][ 4] =  0.15; //DCAV0ToPV 5
    lcuts[lIdx][0][ 5] =0.006;    lcuts[lIdx][1][ 5] = 0.005; lcuts[lIdx][2][ 5] = 0.004; //V0Mass 6
    lcuts[lIdx][0][ 6] = 0.05;    lcuts[lIdx][1][ 6] =  0.10; lcuts[lIdx][2][ 6] =  0.15; //DCABachToPV 7
    lcuts[lIdx][0][ 7] = 1.20;    lcuts[lIdx][1][ 7] =   1.0; lcuts[lIdx][2][ 7] =   0.8; //DCACascDaughters 8
    lcuts[lIdx][0][ 8] =  0.8;    lcuts[lIdx][1][ 8] =   1.2; lcuts[lIdx][2][ 8] =  3.00; //CascRadius 9
    lcuts[lIdx][0][ 9] = 17.5;    lcuts[lIdx][1][ 9] =  15.0; lcuts[lIdx][2][ 9] =  12.5; //ProperLifetime 10
    lcuts[lIdx][0][10] = 40.0;    lcuts[lIdx][1][10] =  30.0; lcuts[lIdx][2][10] =  20.0; //ProperLifetimeV0 11
    lcuts[lIdx][0][11] =   -1;    lcuts[lIdx][1][11] =    70; lcuts[lIdx][2][11] =    80; //LeastNumberOfClusters 12
    lcuts[lIdx][0][12] =    5;    lcuts[lIdx][1][12] =     4; lcuts[lIdx][2][12] =     3; //TPCdEdx 13
    lcuts[lIdx][0][13] =  0.0;    lcuts[lIdx][1][13] = 0.008; lcuts[lIdx][2][13] = 0.010; //Competing 14
    //================================================================================
    
    //================================================================================
    // XIPLUS SELECTIONS
    lIdx = 1; //Master XiPlus Index
    //--------------------------------------------------------------------------------
    //                  LOOSE                        CENTRAL                   TIGHT
    lcuts[lIdx][0][ 0] = 0.10;    lcuts[lIdx][1][ 0] =  0.20; lcuts[lIdx][2][ 0] =  0.30; //DCANegToPV 1
    lcuts[lIdx][0][ 1] = 0.10;    lcuts[lIdx][1][ 1] =  0.20; lcuts[lIdx][2][ 1] =  0.30; //DCAPostToPV 2
    lcuts[lIdx][0][ 2] =  1.2;    lcuts[lIdx][1][ 2] =   1.0; lcuts[lIdx][2][ 2] =   0.8; //DCAV0Daughters 3
    lcuts[lIdx][0][ 3] = 2.00;    lcuts[lIdx][1][ 3] =  3.00; lcuts[lIdx][2][ 3] =   4.0; //V0Radius 4
    lcuts[lIdx][0][ 4] = 0.05;    lcuts[lIdx][1][ 4] =   0.1; lcuts[lIdx][2][ 4] =  0.15; //DCAV0ToPV 5
    lcuts[lIdx][0][ 5] =0.006;    lcuts[lIdx][1][ 5] = 0.005; lcuts[lIdx][2][ 5] = 0.004; //V0Mass 6
    lcuts[lIdx][0][ 6] = 0.05;    lcuts[lIdx][1][ 6] =  0.10; lcuts[lIdx][2][ 6] =  0.15; //DCABachToPV 7
    lcuts[lIdx][0][ 7] = 1.20;    lcuts[lIdx][1][ 7] =   1.0; lcuts[lIdx][2][ 7] =   0.8; //DCACascDaughters 8
    lcuts[lIdx][0][ 8] =  0.8;    lcuts[lIdx][1][ 8] =   1.2; lcuts[lIdx][2][ 8] =  3.00; //CascRadius 9
    lcuts[lIdx][0][ 9] = 17.5;    lcuts[lIdx][1][ 9] =  15.0; lcuts[lIdx][2][ 9] =  12.5; //ProperLifetime 10
    lcuts[lIdx][0][10] = 40.0;    lcuts[lIdx][1][10] =  30.0; lcuts[lIdx][2][10] =  20.0; //ProperLifetimeV0 11
    lcuts[lIdx][0][11] =   -1;    lcuts[lIdx][1][11] =    70; lcuts[lIdx][2][11] =    80; //LeastNumberOfClusters 12
    lcuts[lIdx][0][12] =    5;    lcuts[lIdx][1][12] =     4; lcuts[lIdx][2][12] =     3; //TPCdEdx 13
    lcuts[lIdx][0][13] =  0.0;    lcuts[lIdx][1][13] = 0.008; lcuts[lIdx][2][13] = 0.010; //Competing 14
    //================================================================================
    
    //================================================================================
    // OMEGAMINUS SELECTIONS
    lIdx = 2; //Master OmegaMinus Index
    //--------------------------------------------------------------------------------
    //                  LOOSE                        CENTRAL                   TIGHT
    lcuts[lIdx][0][ 0] = 0.10;    lcuts[lIdx][1][ 0] =  0.20; lcuts[lIdx][2][ 0] =  0.30; //DCANegToPV 1
    lcuts[lIdx][0][ 1] = 0.10;    lcuts[lIdx][1][ 1] =  0.20; lcuts[lIdx][2][ 1] =  0.30; //DCAPostToPV 2
    lcuts[lIdx][0][ 2] =  1.2;    lcuts[lIdx][1][ 2] =   1.0; lcuts[lIdx][2][ 2] =   0.8; //DCAV0Daughters 3
    lcuts[lIdx][0][ 3] = 2.00;    lcuts[lIdx][1][ 3] =  3.00; lcuts[lIdx][2][ 3] =   4.0; //V0Radius 4
    lcuts[lIdx][0][ 4] = 0.05;    lcuts[lIdx][1][ 4] =   0.1; lcuts[lIdx][2][ 4] =  0.15; //DCAV0ToPV 5
    lcuts[lIdx][0][ 5] =0.006;    lcuts[lIdx][1][ 5] = 0.005; lcuts[lIdx][2][ 5] = 0.004; //V0Mass 6
    lcuts[lIdx][0][ 6] = 0.05;    lcuts[lIdx][1][ 6] =  0.10; lcuts[lIdx][2][ 6] =  0.15; //DCABachToPV 7
    lcuts[lIdx][0][ 7] = 1.00;    lcuts[lIdx][1][ 7] =   0.7; lcuts[lIdx][2][ 7] =   0.5; //DCACascDaughters 8
    lcuts[lIdx][0][ 8] =  0.6;    lcuts[lIdx][1][ 8] =   1.0; lcuts[lIdx][2][ 8] =  2.50; //CascRadius 9
    lcuts[lIdx][0][ 9] = 14.0;    lcuts[lIdx][1][ 9] =  12.0; lcuts[lIdx][2][ 9] =  10.0; //ProperLifetime 10
    lcuts[lIdx][0][10] = 40.0;    lcuts[lIdx][1][10] =  30.0; lcuts[lIdx][2][10] =  20.0; //ProperLifetimeV0 11
    lcuts[lIdx][0][11] =   -1;    lcuts[lIdx][1][11] =    70; lcuts[lIdx][2][11] =    80; //LeastNumberOfClusters 12
    lcuts[lIdx][0][12] =    5;    lcuts[lIdx][1][12] =     4; lcuts[lIdx][2][12] =     3; //TPCdEdx 13
    lcuts[lIdx][0][13] =  0.0;    lcuts[lIdx][1][13] = 0.008; lcuts[lIdx][2][13] = 0.010; //Competing 14
    //================================================================================
    
    //================================================================================
    // OMEGAPLUS SELECTIONS
    lIdx = 3; //Master OmegaPlus Index
    //--------------------------------------------------------------------------------
    //                  LOOSE                        CENTRAL                   TIGHT
    lcuts[lIdx][0][ 0] = 0.10;    lcuts[lIdx][1][ 0] =  0.20; lcuts[lIdx][2][ 0] =  0.30; //DCANegToPV 1
    lcuts[lIdx][0][ 1] = 0.10;    lcuts[lIdx][1][ 1] =  0.20; lcuts[lIdx][2][ 1] =  0.30; //DCAPostToPV 2
    lcuts[lIdx][0][ 2] =  1.2;    lcuts[lIdx][1][ 2] =   1.0; lcuts[lIdx][2][ 2] =   0.8; //DCAV0Daughters 3
    lcuts[lIdx][0][ 3] = 2.00;    lcuts[lIdx][1][ 3] =  3.00; lcuts[lIdx][2][ 3] =   4.0; //V0Radius 4
    lcuts[lIdx][0][ 4] = 0.05;    lcuts[lIdx][1][ 4] =   0.1; lcuts[lIdx][2][ 4] =  0.15; //DCAV0ToPV 5
    lcuts[lIdx][0][ 5] =0.006;    lcuts[lIdx][1][ 5] = 0.005; lcuts[lIdx][2][ 5] = 0.004; //V0Mass 6
    lcuts[lIdx][0][ 6] = 0.05;    lcuts[lIdx][1][ 6] =  0.10; lcuts[lIdx][2][ 6] =  0.15; //DCABachToPV 7
    lcuts[lIdx][0][ 7] = 1.00;    lcuts[lIdx][1][ 7] =   0.7; lcuts[lIdx][2][ 7] =   0.5; //DCACascDaughters 8
    lcuts[lIdx][0][ 8] =  0.6;    lcuts[lIdx][1][ 8] =   1.0; lcuts[lIdx][2][ 8] =  2.50; //CascRadius 9
    lcuts[lIdx][0][ 9] = 14.0;    lcuts[lIdx][1][ 9] =  12.0; lcuts[lIdx][2][ 9] =  10.0; //ProperLifetime 10
    lcuts[lIdx][0][10] = 40.0;    lcuts[lIdx][1][10] =  30.0; lcuts[lIdx][2][10] =  20.0; //ProperLifetimeV0 11
    lcuts[lIdx][0][11] =   -1;    lcuts[lIdx][1][11] =    70; lcuts[lIdx][2][11] =    80; //LeastNumberOfClusters 12
    lcuts[lIdx][0][12] =    5;    lcuts[lIdx][1][12] =     4; lcuts[lIdx][2][12] =     3; //TPCdEdx 13
    lcuts[lIdx][0][13] =  0.0;    lcuts[lIdx][1][13] = 0.008; lcuts[lIdx][2][13] = 0.010; //Competing 14
    //================================================================================
    
    //STEP 3: Creation of output objects
    
    //Just a counter and one array, please. Nothing else needed
    AliCascadeResult *lCascadeResult[600];
    Long_t lN = 0;
    
    //Map to mass hypothesis
    AliCascadeResult::EMassHypo lMassHypo[4];
    lMassHypo[0] = AliCascadeResult::kXiMinus;
    lMassHypo[1] = AliCascadeResult::kXiPlus;
    lMassHypo[2] = AliCascadeResult::kOmegaMinus;
    lMassHypo[3] = AliCascadeResult::kOmegaPlus;
    
    
    //Central results: Stored in indices 0, 1, 2, 3 (careful!)
    for(Int_t i = 0 ; i < 4 ; i ++){
        //Central result, customized binning: the one to use, usually
        lCascadeResult[lN] = new AliCascadeResult( Form("%s_Central",lParticleName[i].Data() ),lMassHypo[i],"",lCentbinnumb,lCentbinlimits, lPtbinnumb,lPtbinlimits);
        
        //Setters for V0 Cuts
        lCascadeResult[lN]->SetCutDCANegToPV            ( lcuts[i][1][ 0] ) ;
        lCascadeResult[lN]->SetCutDCAPosToPV            ( lcuts[i][1][ 1] ) ;
        lCascadeResult[lN]->SetCutDCAV0Daughters        ( lcuts[i][1][ 2] ) ;
        lCascadeResult[lN]->SetCutV0Radius              ( lcuts[i][1][ 3] ) ;
        //Setters for Cascade Cuts
        lCascadeResult[lN]->SetCutDCAV0ToPV             ( lcuts[i][1][ 4] ) ;
        lCascadeResult[lN]->SetCutV0Mass                ( lcuts[i][1][ 5] ) ;
        lCascadeResult[lN]->SetCutDCABachToPV           ( lcuts[i][1][ 6] ) ;
        lCascadeResult[lN]->SetCutDCACascDaughters      ( lcuts[i][1][ 7] ) ;
        lCascadeResult[lN]->SetCutCascRadius            ( lcuts[i][1][ 8] ) ;
        //Miscellaneous
        lCascadeResult[lN]->SetCutProperLifetime        ( lcuts[i][1][ 9] ) ;
        lCascadeResult[lN]->SetCutMaxV0Lifetime         ( lcuts[i][1][10] ) ;
        lCascadeResult[lN]->SetCutLeastNumberOfClusters ( lcuts[i][1][11] ) ;
        lCascadeResult[lN]->SetCutTPCdEdx               ( lcuts[i][1][12] ) ;
        lCascadeResult[lN]->SetCutXiRejection           ( lcuts[i][1][13] ) ;
        
        //Parametric angle cut initializations
        //V0 cosine of pointing angle
        lCascadeResult[lN]->SetCutV0CosPA               ( 0.95 ) ; //+variable
        lCascadeResult[lN]->SetCutVarV0CosPA            (TMath::Exp(10.853),
                                                         -25.0322,
                                                         TMath::Exp(-0.843948),
                                                         -0.890794,
                                                         0.057553);
        
        //Cascade cosine of pointing angle
        lCascadeResult[lN]->SetCutCascCosPA             ( 0.95 ) ; //+variable
        lCascadeResult[lN]->SetCutVarCascCosPA          (TMath::Exp(4.86664),
                                                         -10.786,
                                                         TMath::Exp(-1.33411),
                                                         -0.729825,
                                                         0.0695724);
        
        //BB cosine of pointing angle
        lCascadeResult[lN]->SetCutBachBaryonCosPA       ( TMath::Cos(0.04) ) ; //+variable
        lCascadeResult[lN]->SetCutVarBBCosPA            (TMath::Exp(-2.29048),
                                                         -20.2016,
                                                         TMath::Exp(-2.9581),
                                                         -0.649153,
                                                         0.00526455);
        
        //Add result to pool
        lN++;
    }
    
    //Central Full results: Stored in indices 4, 5, 6, 7 (careful!)
    for(Int_t i = 0 ; i < 4 ; i ++){
        lCascadeResult[lN] = new AliCascadeResult( Form("%s_Central_Full",lParticleName[i].Data() ),lMassHypo[i]);
        
        //Setters for V0 Cuts
        lCascadeResult[lN]->SetCutDCANegToPV            ( lcuts[i][1][ 0] ) ;
        lCascadeResult[lN]->SetCutDCAPosToPV            ( lcuts[i][1][ 1] ) ;
        lCascadeResult[lN]->SetCutDCAV0Daughters        ( lcuts[i][1][ 2] ) ;
        lCascadeResult[lN]->SetCutV0Radius              ( lcuts[i][1][ 3] ) ;
        //Setters for Cascade Cuts
        lCascadeResult[lN]->SetCutDCAV0ToPV             ( lcuts[i][1][ 4] ) ;
        lCascadeResult[lN]->SetCutV0Mass                ( lcuts[i][1][ 5] ) ;
        lCascadeResult[lN]->SetCutDCABachToPV           ( lcuts[i][1][ 6] ) ;
        lCascadeResult[lN]->SetCutDCACascDaughters      ( lcuts[i][1][ 7] ) ;
        lCascadeResult[lN]->SetCutCascRadius            ( lcuts[i][1][ 8] ) ;
        //Miscellaneous
        lCascadeResult[lN]->SetCutProperLifetime        ( lcuts[i][1][ 9] ) ;
        lCascadeResult[lN]->SetCutMaxV0Lifetime         ( lcuts[i][1][10] ) ;
        lCascadeResult[lN]->SetCutLeastNumberOfClusters ( lcuts[i][1][11] ) ;
        lCascadeResult[lN]->SetCutTPCdEdx               ( lcuts[i][1][12] ) ;
        lCascadeResult[lN]->SetCutXiRejection           ( lcuts[i][1][13] ) ;
        
        //Parametric angle cut initializations
        //V0 cosine of pointing angle
        lCascadeResult[lN]->SetCutV0CosPA               ( 0.95 ) ; //+variable
        lCascadeResult[lN]->SetCutVarV0CosPA            (TMath::Exp(10.853),
                                                         -25.0322,
                                                         TMath::Exp(-0.843948),
                                                         -0.890794,
                                                         0.057553);
        
        //Cascade cosine of pointing angle
        lCascadeResult[lN]->SetCutCascCosPA             ( 0.95 ) ; //+variable
        lCascadeResult[lN]->SetCutVarCascCosPA          (TMath::Exp(4.86664),
                                                         -10.786,
                                                         TMath::Exp(-1.33411),
                                                         -0.729825,
                                                         0.0695724);
        
        //BB cosine of pointing angle
        lCascadeResult[lN]->SetCutBachBaryonCosPA       ( TMath::Cos(0.04) ) ; //+variable
        lCascadeResult[lN]->SetCutVarBBCosPA            (TMath::Exp(-2.29048),
                                                         -20.2016,
                                                         TMath::Exp(-2.9581),
                                                         -0.649153,
                                                         0.00526455);
        
        //Add result to pool
        lN++;
    }
    
    //Explore restricted rapidity range check
    for(Int_t i = 0 ; i < 4 ; i ++){
        lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_Central_y03",lParticleName[i].Data() ) );
        
        lCascadeResult[lN] -> SetCutMinRapidity(-0.3);
        lCascadeResult[lN] -> SetCutMaxRapidity(+0.3);
        
        //Add result to pool
        lN++;
    }
    
    Float_t lLowRap = -0.6;
    Float_t lHighRap = -0.5;
    for(Int_t i=0;i<4;i++){
        lLowRap = -0.6;
        lHighRap = -0.5;
        for(Int_t irapbin=0;irapbin<12;irapbin++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%f_%f",lParticleName[i].Data(),"DefaultRapiditySweep",lLowRap,lHighRap ) );
            lCascadeResult[lN]->SetCutMinRapidity(lLowRap);
            lCascadeResult[lN]->SetCutMaxRapidity(lHighRap);
            lN++;
            lLowRap+=0.1;
            lHighRap+=0.1;
        }
    }
    
    // STEP 4: Creation of objects to be used in systematics
    // Optimized via use of copy constructors
    for(Int_t i = 0 ; i < 4 ; i ++){
        for(Int_t iCut = 0 ; iCut < 14 ; iCut ++){
            
            //LOOSE VARIATIONS
            //Create a new object from default
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%s",lParticleName[i].Data(),lCutName[iCut].Data(),lConfName[0].Data()) );
            
            if(iCut ==  0 ) lCascadeResult[lN]->SetCutDCANegToPV            ( lcuts[i][0][iCut] ) ;
            if(iCut ==  1 ) lCascadeResult[lN]->SetCutDCAPosToPV            ( lcuts[i][0][iCut] ) ;
            if(iCut ==  2 ) lCascadeResult[lN]->SetCutDCAV0Daughters        ( lcuts[i][0][iCut] ) ;
            if(iCut ==  3 ) lCascadeResult[lN]->SetCutV0Radius              ( lcuts[i][0][iCut] ) ;
            
            //Setters for Cascade Cuts
            if(iCut ==  4 ) lCascadeResult[lN]->SetCutDCAV0ToPV             ( lcuts[i][0][iCut] ) ;
            if(iCut ==  5 ) lCascadeResult[lN]->SetCutV0Mass                ( lcuts[i][0][iCut] ) ;
            if(iCut ==  6 ) lCascadeResult[lN]->SetCutDCABachToPV           ( lcuts[i][0][iCut] ) ;
            if(iCut ==  7 ) lCascadeResult[lN]->SetCutDCACascDaughters      ( lcuts[i][0][iCut] ) ;
            if(iCut ==  8 ) lCascadeResult[lN]->SetCutCascRadius            ( lcuts[i][0][iCut] ) ;
            
            //Miscellaneous
            if(iCut ==  9 ) lCascadeResult[lN]->SetCutProperLifetime        ( lcuts[i][0][iCut] ) ;
            if(iCut == 10 ) lCascadeResult[lN]->SetCutMaxV0Lifetime         ( lcuts[i][0][iCut] ) ;
            if(iCut == 11 ) lCascadeResult[lN]->SetCutLeastNumberOfClusters ( lcuts[i][0][iCut] ) ;
            if(iCut == 12 ) lCascadeResult[lN]->SetCutTPCdEdx               ( lcuts[i][0][iCut] ) ;
            if(iCut == 13 ) lCascadeResult[lN]->SetCutXiRejection           ( lcuts[i][0][iCut] ) ;
            
            //Print this variation, add to pool
            //lCascadeResult[lN]->Print();
            lN++;
            
            //TIGHT VARIATIONS
            //Create a new object from default
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%s",lParticleName[i].Data(),lCutName[iCut].Data(),lConfName[2].Data()) );
            
            if(iCut ==  0 ) lCascadeResult[lN]->SetCutDCANegToPV            ( lcuts[i][2][iCut] ) ;
            if(iCut ==  1 ) lCascadeResult[lN]->SetCutDCAPosToPV            ( lcuts[i][2][iCut] ) ;
            if(iCut ==  2 ) lCascadeResult[lN]->SetCutDCAV0Daughters        ( lcuts[i][2][iCut] ) ;
            if(iCut ==  3 ) lCascadeResult[lN]->SetCutV0Radius              ( lcuts[i][2][iCut] ) ;
            
            //Setters for Cascade Cuts
            if(iCut ==  4 ) lCascadeResult[lN]->SetCutDCAV0ToPV             ( lcuts[i][2][iCut] ) ;
            if(iCut ==  5 ) lCascadeResult[lN]->SetCutV0Mass                ( lcuts[i][2][iCut] ) ;
            if(iCut ==  6 ) lCascadeResult[lN]->SetCutDCABachToPV           ( lcuts[i][2][iCut] ) ;
            if(iCut ==  7 ) lCascadeResult[lN]->SetCutDCACascDaughters      ( lcuts[i][2][iCut] ) ;
            if(iCut ==  8 ) lCascadeResult[lN]->SetCutCascRadius            ( lcuts[i][2][iCut] ) ;
            
            //Miscellaneous
            if(iCut ==  9 ) lCascadeResult[lN]->SetCutProperLifetime        ( lcuts[i][2][iCut] ) ;
            if(iCut == 10 ) lCascadeResult[lN]->SetCutMaxV0Lifetime         ( lcuts[i][2][iCut] ) ;
            if(iCut == 11 ) lCascadeResult[lN]->SetCutLeastNumberOfClusters ( lcuts[i][2][iCut] ) ;
            if(iCut == 12 ) lCascadeResult[lN]->SetCutTPCdEdx               ( lcuts[i][2][iCut] ) ;
            if(iCut == 13 ) lCascadeResult[lN]->SetCutXiRejection           ( lcuts[i][2][iCut] ) ;
            
            //Print this variation, add to pool
            //lCascadeResult[lN]->Print();
            lN++;
        }
    }
    
    //STEP 5: re-parametrization of cosines for tight and loose variations (done manually)
    for(Int_t i = 0 ; i < 4 ; i ++){
        //======================================================
        //V0CosPA Variations
        //======================================================
        lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%s",lParticleName[i].Data(),"V0CosPA","Loose") );
        lCascadeResult[lN]->SetCutVarV0CosPA(TMath::Exp(  -1.77429),
                                             -0.692453,
                                             TMath::Exp( -2.01938),
                                             -0.201574,
                                             0.0776465);
        lN++;
        lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%s",lParticleName[i].Data(),"V0CosPA","Tight") );
        lCascadeResult[lN]->SetCutVarV0CosPA(TMath::Exp(  -1.21892),
                                             -41.8521,
                                             TMath::Exp(   -1.278),
                                             -0.894064,
                                             0.0303932);
        lN++;
        lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%s",lParticleName[i].Data(),"V0CosPA","VeryTight") );
        lCascadeResult[lN]->SetCutVarV0CosPA(TMath::Exp(   12.8077),
                                             -21.2944,
                                             TMath::Exp( -1.53357),
                                             -0.920017,
                                             0.0262315);
        lN++;
        //======================================================
        //CascCosPA Variations
        //======================================================
        lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%s",lParticleName[i].Data(),"CascCosPA","Loose") );
        lCascadeResult[lN]->SetCutVarCascCosPA(TMath::Exp(  -1.77429),
                                               -0.692453,
                                               TMath::Exp( -2.01938),
                                               -0.201574,
                                               0.0776465);
        lN++;
        lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%s",lParticleName[i].Data(),"CascCosPA","Tight") );
        lCascadeResult[lN]->SetCutVarCascCosPA(TMath::Exp(   12.8752),
                                               -21.522,
                                               TMath::Exp( -1.49906),
                                               -0.813472,
                                               0.0480962);
        lN++;
        lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%s",lParticleName[i].Data(),"CascCosPA","VeryTight") );
        lCascadeResult[lN]->SetCutVarCascCosPA(TMath::Exp(    12.801),
                                               -21.6157,
                                               TMath::Exp( -1.66297),
                                               -0.889246,
                                               0.0346838);
        lN++;
        //======================================================
        //BBCosPA Variations
        //======================================================
        lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%s",lParticleName[i].Data(),"BBCosPA","Loose") );
        lCascadeResult[lN]->SetCutBachBaryonCosPA        ( TMath::Cos(0.03) ) ;
        lCascadeResult[lN]->SetCutVarBBCosPA(TMath::Exp(    -2.8798),
                                             -20.9876,
                                             TMath::Exp(  -3.10847),
                                             -0.73045,
                                             0.00235147);
        lN++;
        lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%s",lParticleName[i].Data(),"BBCosPA","Tight") );
        lCascadeResult[lN]->SetCutBachBaryonCosPA        ( TMath::Cos(0.05) ) ;
        lCascadeResult[lN]->SetCutVarBBCosPA(TMath::Exp(   12.4606),
                                             -20.578,
                                             TMath::Exp( -2.41442),
                                             -0.709588,
                                             0.01079);
        lN++;
    }
    
    //STEP 6: V0 Mass sweep
    //for(Int_t i = 0 ; i < 4 ; i ++){
    //    for(Int_t isweep=0; isweep<20;isweep++){
    //        lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_V0MassSweep_%i",lParticleName[i].Data(),isweep) );
    //        lCascadeResult[lN]->SetCutV0MassSigma( ((Double_t)(isweep)/4000.0)); //in GeV/c^2
    //        lN++;
    //    }
    //}
    
    Float_t lLifetimeCut[4];
    lLifetimeCut[0] = 15.0;
    lLifetimeCut[1] = 15.0;
    lLifetimeCut[2] = 12.0;
    lLifetimeCut[3] = 12.0;
    
    Float_t lMass[4];
    lMass[0] = 1.322;
    lMass[1] = 1.322;
    lMass[2] = 1.672;
    lMass[3] = 1.672;
    
    //Old vertexer-level configuration for cross-checks
    for(Int_t i = 0 ; i < 4 ; i ++){
        //Central result, customized binning: the one to use, usually
        lCascadeResult[lN] = new AliCascadeResult( Form("%s_VertexerLevel",lParticleName[i].Data() ),lMassHypo[i],"",lCentbinnumb,lCentbinlimits, lPtbinnumb,lPtbinlimits,100,lMass[i]-0.050,lMass[i]+0.050);
        
        //Default cuts: use vertexer level ones
        //Setters for V0 Cuts
        lCascadeResult[lN]->SetCutDCANegToPV            ( 0.2 ) ;
        lCascadeResult[lN]->SetCutDCAPosToPV            ( 0.2 ) ;
        lCascadeResult[lN]->SetCutDCAV0Daughters        (  1. ) ;
        lCascadeResult[lN]->SetCutV0CosPA               ( 0.98 ) ;
        lCascadeResult[lN]->SetCutV0Radius              (  3 ) ;
        //Setters for Cascade Cuts
        lCascadeResult[lN]->SetCutDCAV0ToPV             ( 0.1 ) ;
        lCascadeResult[lN]->SetCutV0Mass                ( 0.006 ) ;
        lCascadeResult[lN]->SetCutDCABachToPV           ( 0.03 ) ;
        lCascadeResult[lN]->SetCutDCACascDaughters      ( 1. ) ;
        lCascadeResult[lN]->SetCutCascRadius            ( 1.2 ) ;
        if(i==2||i==3)
            lCascadeResult[lN]->SetCutCascRadius            ( 1.0 ) ; //omega case
        lCascadeResult[lN]->SetCutCascCosPA             ( 0.98 ) ;
        //Miscellaneous
        lCascadeResult[lN]->SetCutProperLifetime        ( lLifetimeCut[i] ) ;
        lCascadeResult[lN]->SetCutLeastNumberOfClusters ( 70.0 ) ;
        lCascadeResult[lN]->SetCutTPCdEdx               ( 4.0 ) ;
        lCascadeResult[lN]->SetCutXiRejection           ( 0.008 ) ;
        lCascadeResult[lN]->SetCutBachBaryonCosPA        ( TMath::Cos(0.006) ) ;
        //Add result to pool
        lN++;
    }
    
    cout<<"Added "<<lN<<" Cascade configurations to output."<<endl;
    for (Int_t iconf = 0; iconf<lN; iconf++)
        AddConfiguration(lCascadeResult[iconf]);
    
    //Add V0 configuration for the determination of mass and width of the Lambda peak
    
    //Map to mass hypothesis
    AliV0Result::EMassHypo lMassHypoV0[3];
    lMassHypoV0[0] = AliV0Result::kK0Short;
    lMassHypoV0[1] = AliV0Result::kLambda;
    lMassHypoV0[2] = AliV0Result::kAntiLambda;
    
    //Array of results
    AliV0Result *lV0Result[500];
    Long_t lNV0 = 0;
    
    TString lParticleNameV0[] = {"K0Short", "Lambda",  "AntiLambda"};
    
    //Central results: Stored in indices 0, 1, 2 (careful!)
    for(Int_t i = 0 ; i < 3 ; i ++){
        //Central result, customized binning: the one to use, usually
        lV0Result[lNV0] = new AliV0Result( Form("%s_Central",lParticleNameV0[i].Data() ),lMassHypoV0[i],"",lCentbinnumb,lCentbinlimits, lPtbinnumb,lPtbinlimits);
        
        //Setters for V0 Cuts
        lV0Result[lNV0]->SetCutDCANegToPV            ( lcuts[0][1][ 0] ) ;
        lV0Result[lNV0]->SetCutDCAPosToPV            ( lcuts[0][1][ 1] ) ;
        lV0Result[lNV0]->SetCutDCAV0Daughters        ( lcuts[0][1][ 2] ) ;
        lV0Result[lNV0]->SetCutV0CosPA               ( 0.95 ) ; //+variable
        lV0Result[lNV0]->SetCutVarV0CosPA            (TMath::Exp(10.853),
                                                      -25.0322,
                                                      TMath::Exp(-0.843948),
                                                      -0.890794,
                                                      0.057553);
        
        lV0Result[lNV0]->SetCutV0Radius              ( lcuts[0][1][ 3] ) ;
        
        //Add result to pool
        lNV0++;
    }
    
    for (Int_t iconf = 0; iconf<lNV0; iconf++)
        AddConfiguration(lV0Result[iconf]);
    
    cout<<"Added "<<lN<<" Cascade configurations to output."<<endl;
    cout<<"Added "<<lNV0<<" V0 configurations to output."<<endl;
}


//________________________________________________________________________
void AliAnalysisTaskStrangenessVsMultiplicityRun2pPb::AddCascadeConfiguration276TeV()
//Adds 2.76 TeV cascade analysis configuration
{
    // STEP 1: Decide on binning (needed to improve on memory consumption)
    
    // pT binning
    Double_t lPtbinlimits[] = {0.4, 0.5, 0.6,
        0.7,0.8,.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,
        2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,
        4.4,4.6,4.8,5.0,5.5,6.0,6.5,7.0,8.0,9.0,10.,11.,12.};
    Long_t lPtbinnumb = sizeof(lPtbinlimits)/sizeof(Double_t) - 1;
    
    // centrality binning
    Double_t lCentbinlimits[] = {0, 1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 85, 90};
    Long_t lCentbinnumb = sizeof(lCentbinlimits)/sizeof(Double_t) - 1;
    
    // TStrings for output names
    TString lParticleName[] = {"XiMinus", "XiPlus",  "OmegaMinus", "OmegaPlus"};
    
    //Just a counter and one array, please. Nothing else needed
    AliCascadeResult *lCascadeResult[100];
    Long_t lN = 0;
    
    //Map to mass hypothesis
    AliCascadeResult::EMassHypo lMassHypo[4];
    lMassHypo[0] = AliCascadeResult::kXiMinus;
    lMassHypo[1] = AliCascadeResult::kXiPlus;
    lMassHypo[2] = AliCascadeResult::kOmegaMinus;
    lMassHypo[3] = AliCascadeResult::kOmegaPlus;
    
    for(Int_t i = 0 ; i < 4 ; i ++){
        //2.76 Config result, customized binning: the one to use, usually
        lCascadeResult[lN] = new AliCascadeResult( Form("%s_276TeV",lParticleName[i].Data() ),lMassHypo[i],"",lCentbinnumb,lCentbinlimits, lPtbinnumb,lPtbinlimits);
        
        //Setters for V0 Cuts
        lCascadeResult[lN]->SetCutDCANegToPV            ( 0.1    ) ;
        lCascadeResult[lN]->SetCutDCAPosToPV            ( 0.1    ) ;
        lCascadeResult[lN]->SetCutDCAV0Daughters        ( 0.8    ) ;
        lCascadeResult[lN]->SetCutV0CosPA               ( 0.95   ) ; // + variable
        lCascadeResult[lN]->SetCutUse276TeVV0CosPA      ( kTRUE  ) ;
        lCascadeResult[lN]->SetCutV0Radius              ( 3.0    ) ;
        //Setters for Cascade Cuts
        lCascadeResult[lN]->SetCutDCAV0ToPV             ( 0.1    ) ;
        lCascadeResult[lN]->SetCutV0Mass                ( 0.005  ) ;
        lCascadeResult[lN]->SetCutDCABachToPV           ( 0.03   ) ;
        lCascadeResult[lN]->SetCutDCACascDaughters      ( 0.3    ) ;
        lCascadeResult[lN]->SetCutCascRadius            ( 1.5    ) ;
        lCascadeResult[lN]->SetCutCascCosPA             ( 0.9992 ) ;
        //Miscellaneous
        lCascadeResult[lN]->SetCutProperLifetime        ( 15.0   ) ;
        lCascadeResult[lN]->SetCutLeastNumberOfClusters ( 70     ) ;
        lCascadeResult[lN]->SetCutTPCdEdx               ( 4      ) ;
        lCascadeResult[lN]->SetCutXiRejection           ( 0.008  ) ;
        lCascadeResult[lN]->SetCutDCABachToBaryon       ( 0      ) ;
        
        if(i > 1){
            lCascadeResult[lN]->SetCutCascRadius            ( 1.0 ) ;
            lCascadeResult[lN]->SetCutProperLifetime        ( 8.0 ) ;
        }
        
        //Add result to pool
        lN++;
    }
    
    //Explore restricted rapidity range check
    for(Int_t i = 0 ; i < 4 ; i ++){
        lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_276TeV_y03",lParticleName[i].Data() ) );
        
        lCascadeResult[lN] -> SetCutMinRapidity(-0.3);
        lCascadeResult[lN] -> SetCutMaxRapidity(+0.3);
        
        //Add result to pool
        lN++;
    }
    
    Float_t lLowRap = -0.6;
    Float_t lHighRap = -0.5;
    for(Int_t i=0;i<4;i++){
        lLowRap = -0.6;
        lHighRap = -0.5;
        for(Int_t irapbin=0;irapbin<12;irapbin++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%f_%f",lParticleName[i].Data(),"276TeVRapiditySweep",lLowRap,lHighRap ) );
            lCascadeResult[lN]->SetCutMinRapidity(lLowRap);
            lCascadeResult[lN]->SetCutMaxRapidity(lHighRap);
            lN++;
            lLowRap+=0.1;
            lHighRap+=0.1;
        }
    }
    
    for (Int_t iconf = 0; iconf<lN; iconf++)
        AddConfiguration(lCascadeResult[iconf]);
    
    cout<<"Added "<<lN<<" cascade configurations to output (corresponding to 2.76 TeV analysis cuts)"<<endl;
}


//________________________________________________________________________
Float_t AliAnalysisTaskStrangenessVsMultiplicityRun2pPb::GetDCAz(AliESDtrack *lTrack)
//Encapsulation of DCAz calculation
{
    Float_t b[2];
    Float_t bCov[3];
    lTrack->GetImpactParameters(b,bCov);
    if (bCov[0]<=0 || bCov[2]<=0) {
        AliDebug(1, "Estimated b resolution lower or equal to zero!");
        bCov[0]=0; bCov[2]=0;
    }
    Float_t dcaToVertexXY = b[0];
    Float_t dcaToVertexZ = b[1];
    
    return dcaToVertexZ;
}


//________________________________________________________________________
Float_t AliAnalysisTaskStrangenessVsMultiplicityRun2pPb::GetCosPA(AliESDtrack *lPosTrack, AliESDtrack *lNegTrack, AliESDEvent *lEvent)
//Encapsulation of CosPA calculation (warning: considers AliESDtrack clones)
{
    Float_t lCosPA = -1;
    AliESDtrack* lNegClone = (AliESDtrack*) lNegTrack->Clone("lNegClone"); //need clone, in order not to change track parameters
    AliESDtrack* lPosClone = (AliESDtrack*) lPosTrack->Clone("lPosClone"); //need clone, in order not to change track parameters
    
    //Get Magnetic field and primary vertex
    Double_t b=lEvent->GetMagneticField();
    const AliESDVertex *vtxT3D=lEvent->GetPrimaryVertex();
    Double_t xPrimaryVertex=vtxT3D->GetX();
    Double_t yPrimaryVertex=vtxT3D->GetY();
    Double_t zPrimaryVertex=vtxT3D->GetZ();
    
    //Get ExternalTrackParam
    AliExternalTrackParam nt(*lNegClone), pt(*lPosClone);
    
    //Find DCA
    Double_t xn, xp, dca=lNegClone->GetDCA(lPosClone,b,xn,xp);
    
    //Propagate to it
    nt.PropagateTo(xn,b); pt.PropagateTo(xp,b);
    
    //Create V0 object to do propagation
    AliESDv0 vertex(nt,1,pt,2); //Never mind indices, won't use
    
    //Get CosPA
    lCosPA=vertex.GetV0CosineOfPointingAngle(xPrimaryVertex,yPrimaryVertex,zPrimaryVertex);
    
    //Cleanup
    delete lNegClone;
    delete lPosClone;
    
    //Return value
    return lCosPA;
}

//________________________________________________________________________
void AliAnalysisTaskStrangenessVsMultiplicityRun2pPb::CheckChargeV0(AliESDv0 *v0)
{
    // This function checks charge of negative and positive daughter tracks.
    // If incorrectly defined (onfly vertexer), swaps out.
    if( v0->GetParamN()->Charge() > 0 && v0->GetParamP()->Charge() < 0 ) {
        //V0 daughter track swapping is required! Note: everything is swapped here... P->N, N->P
        Long_t lCorrectNidx = v0->GetPindex();
        Long_t lCorrectPidx = v0->GetNindex();
        Double32_t	lCorrectNmom[3];
        Double32_t	lCorrectPmom[3];
        v0->GetPPxPyPz( lCorrectNmom[0], lCorrectNmom[1], lCorrectNmom[2] );
        v0->GetNPxPyPz( lCorrectPmom[0], lCorrectPmom[1], lCorrectPmom[2] );
        
        AliExternalTrackParam	lCorrectParamN(
                                               v0->GetParamP()->GetX() ,
                                               v0->GetParamP()->GetAlpha() ,
                                               v0->GetParamP()->GetParameter() ,
                                               v0->GetParamP()->GetCovariance()
                                               );
        AliExternalTrackParam	lCorrectParamP(
                                               v0->GetParamN()->GetX() ,
                                               v0->GetParamN()->GetAlpha() ,
                                               v0->GetParamN()->GetParameter() ,
                                               v0->GetParamN()->GetCovariance()
                                               );
        lCorrectParamN.SetMostProbablePt( v0->GetParamP()->GetMostProbablePt() );
        lCorrectParamP.SetMostProbablePt( v0->GetParamN()->GetMostProbablePt() );
        
        //Get Variables___________________________________________________
        Double_t lDcaV0Daughters = v0 -> GetDcaV0Daughters();
        Double_t lCosPALocal     = v0 -> GetV0CosineOfPointingAngle();
        Bool_t lOnFlyStatusLocal = v0 -> GetOnFlyStatus();
        
        //Create Replacement Object_______________________________________
        AliESDv0 *v0correct = new AliESDv0(lCorrectParamN,lCorrectNidx,lCorrectParamP,lCorrectPidx);
        v0correct->SetDcaV0Daughters          ( lDcaV0Daughters   );
        v0correct->SetV0CosineOfPointingAngle ( lCosPALocal       );
        v0correct->ChangeMassHypothesis       ( kK0Short          );
        v0correct->SetOnFlyStatus             ( lOnFlyStatusLocal );
        
        //Reverse Cluster info..._________________________________________
        v0correct->SetClusters( v0->GetClusters( 1 ), v0->GetClusters ( 0 ) );
        
        *v0 = *v0correct;
        //Proper cleanup..._______________________________________________
        v0correct->Delete();
        v0correct = 0x0;
        
        //Just another cross-check and output_____________________________
        if( v0->GetParamN()->Charge() > 0 && v0->GetParamP()->Charge() < 0 ) {
            AliWarning("Found Swapped Charges, tried to correct but something FAILED!");
        } else {
            //AliWarning("Found Swapped Charges and fixed.");
        }
        //________________________________________________________________
    } else {
        //Don't touch it! ---
        //Printf("Ah, nice. Charges are already ordered...");
    }
    return;
}
