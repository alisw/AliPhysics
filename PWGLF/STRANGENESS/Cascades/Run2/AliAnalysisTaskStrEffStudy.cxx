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

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// This task is meant to study V0 and cascade detection efficiencies
// decomposing tracking and pos/neg, V0/bach pairing efficiencies.
//
// The main intent is to compare 'findable' V0s and cascades, i.e. those
// whose daughter particles have been successfully tracked, and correctly
// reconstructed candidates.
//
// The output is composed of two TTree objects storing all 'findable'
// candidates and all reconstruction info (if available) and intermediate
// vertexing information to determine why that particular findable V0 or
// cascade was not found by ALICE.
//
//   --- Questions? Bugs? Please write to:
//       david.dobrigkeit.chinellato@cern.ch
//
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

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

#include <array>
#include <numeric>
#include <utility>

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
#include "AliAnalysisTaskStrEffStudy.h"

using std::cout;
using std::endl;

namespace {
    struct TrackMC {
        AliESDtrack* track;
        TParticle* mother;
        TParticle* particle;
        int motherId;
    };

    struct CandidateMC {
        AliESDtrack *track_deu;
        AliESDtrack *track_p;
        AliESDtrack *track_pi;
        TParticle *part1;
        TParticle *part2;
        TParticle *part3;
        TParticle *mother;
        int motherId;
    };
}

ClassImp(AliAnalysisTaskStrEffStudy)

AliAnalysisTaskStrEffStudy::AliAnalysisTaskStrEffStudy()
: AliAnalysisTaskSE(), fListHist(0), fListV0(0), fListCascade(0), fTreeEvent(0), fTreeV0(0), fTreeCascade(0), fTreeHyperTriton3Body(nullptr), fPIDResponse(0), fESDtrackCuts(0), fESDtrackCutsITSsa2010(0), fESDtrackCutsGlobal2015(0), fUtils(0), fEventCuts(), fRand(0),

//---> Flags controlling Event Tree output
fkSaveEventTree    ( kTRUE ), //no downscaling in this tree so far

//---> Flags controlling V0 TTree output
fkSaveV0Tree       ( kTRUE ),
fkDownScaleV0      ( kTRUE  ),
fDownScaleFactorV0 ( 0.001  ),
fkPreselectDedx ( kFALSE ),
fkPreselectPID  ( kTRUE  ),
fkUseOnTheFlyV0Cascading( kFALSE ),
fkDoImprovedCascadeVertexFinding(kFALSE),
fkDoImprovedDCAV0DauPropagation(kFALSE),
fkDoImprovedDCACascDauPropagation(kFALSE),
fkIfImprovedPerformInitialLinearPropag( kFALSE ),
fkIfImprovedExtraPrecisionFactor ( 1.0 ),
fkDebugWrongPIDForTracking ( kFALSE ),
fkDebugBump( kFALSE ),
fkDebugOOBPileup(kFALSE),
fkDoExtraEvSels( kTRUE ),

//---> Flags controlling Cascade TTree output
fkSaveCascadeTree       ( kTRUE  ),
fkDownScaleCascade      ( kTRUE  ),
fDownScaleFactorCascade ( 0.001  ),
fMinPtToSave( 0.00   ) ,
fMaxPtToSave( 100.00 ) ,

//---> Flags controlling HyperTriton3Body TTree output
fkSaveHyperTriton3BodyTree(false),
fkDownScaleHyperTriton3Body(false),
fDownScaleFactorHyperTriton3Body(1.),

//---> Flags controlling Vertexers
fkRunVertexers    ( kFALSE ),
fkUseLightVertexer ( kTRUE ),
fkDoV0Refit ( kTRUE ),
fkExtraCleanup    ( kTRUE ),
fkSaveGoodTracks( kTRUE ),
fkSandboxV0( kTRUE ),
fkSandboxCascade( kFALSE ),

//---> Flag controlling trigger selection
fTrigType(AliVEvent::kMB),

fPrecisionCutoffCascadeDCA(1e-4),
fMaxIterationsCascadeDCA(27),

//---> Variables for fTreeEvent
fCentrality(0),
fMVPileupFlag(kFALSE),
fOOBPileupFlag(kFALSE),
fNTOFClusters(-1),
fNTOFMatches(-1),
fNTracksITSsa2010(-1),
fNTracksGlobal2015(-1),
fNTracksGlobal2015TriggerPP(-1),
fAmplitudeV0A(-1.),
fAmplitudeV0C(-1.),
fNHitsFMDA(-1.),
fNHitsFMDC(-1.),

//---> Variables for fTreeV0
fTreeVariableGoodV0(kFALSE),
fTreeVariableCentrality(0),
fTreeVariablePosLength(0),
fTreeVariableNegLength(0),
fTreeVariablePosCrossedRows(0),
fTreeVariableNegCrossedRows(0),
fTreeVariablePosTrackStatus(0),
fTreeVariableNegTrackStatus(0),
fTreeVariablePosDxy(0),
fTreeVariableNegDxy(0),
fTreeVariablePosDz(0),
fTreeVariableNegDz(0),
fTreeVariableDcaV0Daughters(0),
fTreeVariableDcaV0DaughtersGeometric(0),
fTreeVariablePosPropagStatus(0),
fTreeVariableNegPropagStatus(0),
fTreeVariableV0Radius(0),
fTreeVariableV0CosineOfPointingAngle(0),
fTreeVariableDecayX(0),
fTreeVariableDecayY(0),
fTreeVariableDecayZ(0),
fTreeVariableDecayXMC(0),
fTreeVariableDecayYMC(0),
fTreeVariableDecayZMC(0),
fTreeVariableNegPxMC(0),
fTreeVariableNegPyMC(0),
fTreeVariableNegPzMC(0),
fTreeVariablePosPxMC(0),
fTreeVariablePosPyMC(0),
fTreeVariablePosPzMC(0),
fTreeVariableInvMassK0s(0),
fTreeVariableInvMassLambda(0),
fTreeVariableInvMassAntiLambda(0),
fTreeVariablePID(0),
fTreeVariablePIDPositive(0),
fTreeVariablePIDNegative(0),
fTreeVariablePtMC(0),
fTreeVariableRapMC(0),
fTreeVariableNegTOFSignal(99999), 
fTreeVariablePosTOFSignal(99999),
fTreeVariablePosAlpha(0),
fTreeVariablePosSigmaY2(0),
fTreeVariablePosSigmaZ2(0),
fTreeVariableNegAlpha(0),
fTreeVariableNegSigmaY2(0),
fTreeVariableNegSigmaZ2(0),

fTreeVariablePosTrack(0x0),
fTreeVariableNegTrack(0x0),
fTreeVariableOTFV0(0x0),
fTreeVariableFoundOTFV0(kFALSE),
fTreeVariableMagneticField(0),
fTreeVariablePosOriginalX(0),
fTreeVariableNegOriginalX(0),
fTreeVariablePVx(0),
fTreeVariablePVy(0),
fTreeVariablePVz(0),
fTreeVariableAliESDvertex(0),
fTreeVariableRun(0),

//---> Variables for fTreeCascade
fTreeCascVarCentrality(0),
fTreeCascVarPosSign(0),
fTreeCascVarNegSign(0),
fTreeCascVarBachSign(0),
fTreeCascVarPosLength(0),
fTreeCascVarNegLength(0),
fTreeCascVarBachLength(0),
fTreeCascVarPosCrossedRows(0),
fTreeCascVarNegCrossedRows(0),
fTreeCascVarBachCrossedRows(0),
//Tracking flags
fTreeCascVarPosTrackStatus(0),
fTreeCascVarNegTrackStatus(0),
fTreeCascVarBachTrackStatus(0),
//DCAxy to PV
fTreeCascVarPosDxy(0),
fTreeCascVarNegDxy(0),
fTreeCascVarBachDxy(0),
//DCAz
fTreeCascVarPosDz(0),
fTreeCascVarNegDz(0),
fTreeCascVarBachDz(0),
fTreeCascVarDcaV0Daughters(0),
fTreeCascVarNegPropagStatus(0),
fTreeCascVarPosPropagStatus(0),
fTreeCascVarV0Radius(0),
fTreeCascVarV0DecayX(0),
fTreeCascVarV0DecayY(0),
fTreeCascVarV0DecayZ(0),
fTreeCascVarV0DecayXMC(0),
fTreeCascVarV0DecayYMC(0),
fTreeCascVarV0DecayZMC(0),
fTreeCascVarV0CosineOfPointingAngle(0),
fTreeCascVarDCAV0ToPrimVtx(0),
fTreeCascVarDCAxyV0ToPrimVtx(0),
fTreeCascVarInvMassLambda(0),
fTreeCascVarInvMassAntiLambda(0),

fTreeCascVarDCACascDaughters(0),
fTreeCascVarCascPropagation(0),

fTreeCascVarDecayX(0),
fTreeCascVarDecayY(0),
fTreeCascVarDecayZ(0),
fTreeCascVarDecayXMC(0),
fTreeCascVarDecayYMC(0),
fTreeCascVarDecayZMC(0),
fTreeCascVarCascCosPointingAngle(0),

fTreeCascVarInvMassXiMinus(0),
fTreeCascVarInvMassXiPlus(0),
fTreeCascVarInvMassOmegaMinus(0),
fTreeCascVarInvMassOmegaPlus(0),

fTreeCascVarPIDPositive(0),
fTreeCascVarPIDNegative(0),
fTreeCascVarPIDBachelor(0),
fTreeCascVarPID(0),
fTreeCascVarPtMC(0),
fTreeCascVarRapMC(0),

fTreeCascVarNegTOFSignal(99999),
fTreeCascVarPosTOFSignal(99999),
fTreeCascVarBachTOFSignal(99999),

fTreeCascVarPosDistanceToTrueDecayPt(0),
fTreeCascVarNegDistanceToTrueDecayPt(0),
fTreeCascVarBachDistanceToTrueDecayPt(0),
fTreeCascVarV0DistanceToTrueDecayPt(0),

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
fTreeCascVarNegPxMC(0),
fTreeCascVarNegPyMC(0),
fTreeCascVarNegPzMC(0),
fTreeCascVarPosPxMC(0),
fTreeCascVarPosPyMC(0),
fTreeCascVarPosPzMC(0),
fTreeCascVarBachPxMC(0),
fTreeCascVarBachPyMC(0),
fTreeCascVarBachPzMC(0),

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//Save full info for full re-vertex offline replay ('sandbox mode')
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
fTreeCascVarBachTrack(0),
fTreeCascVarPosTrack(0),
fTreeCascVarNegTrack(0),
fTreeCascVarOTFV0(0x0),
fTreeCascVarOTFV0NegBach(0x0),
fTreeCascVarOTFV0PosBach(0x0),
fTreeCascVarV0AsOTF(0),
fTreeCascVarNegBachAsOTF(0),
fTreeCascVarPosBachAsOTF(0),
fTreeCascVarMagneticField(0),
fTreeCascVarBachOriginalX(0),
fTreeCascVarPosOriginalX(0),
fTreeCascVarNegOriginalX(0),
fTreeCascVarPVx(0),
fTreeCascVarPVy(0),
fTreeCascVarPVz(0),
fTreeCascVarAliESDvertex(0),

/// Hypertriton 3 Body sandbox
fTreeHyp3BodyVarTracks(),
fTreeHyp3BodyVarPDGcodes(),
fTreeHyp3BodyVarEventId(0u),
fTreeHyp3BodyVarMotherId(0),
fTreeHyp3BodyVarTruePx(0.f),
fTreeHyp3BodyVarTruePy(0.f),
fTreeHyp3BodyVarTruePz(0.f),
fTreeHyp3BodyVarDecayVx(0.f),
fTreeHyp3BodyVarDecayVy(0.f),
fTreeHyp3BodyVarDecayVz(0.f),
fTreeHyp3BodyVarDecayT(0.f),
fTreeHyp3BodyVarPVx(0.f),
fTreeHyp3BodyVarPVy(0.f),
fTreeHyp3BodyVarPVz(0.f),
fTreeHyp3BodyVarPVt(0.f),
fTreeHyp3BodyVarMagneticField(0.f),

//Histos
fHistEventCounter(0),
fHistCentrality(0),
//V0s
fHistGeneratedPtVsYVsCentralityK0Short(0),
fHistGeneratedPtVsYVsCentralityLambda(0),
fHistGeneratedPtVsYVsCentralityAntiLambda(0),
//Cascades
fHistGeneratedPtVsYVsCentralityXiMinus(0),
fHistGeneratedPtVsYVsCentralityXiPlus(0),
fHistGeneratedPtVsYVsCentralityOmegaMinus(0),
fHistGeneratedPtVsYVsCentralityOmegaPlus(0),
//Hypertriton
fHistGeneratedPtVsYVsCentralityHypTrit(0),
fHistGeneratedPtVsYVsCentralityAntiHypTrit(0)
//------------------------------------------------
// Tree Variables
{

}

AliAnalysisTaskStrEffStudy::AliAnalysisTaskStrEffStudy(Bool_t lSaveEventTree, Bool_t lSaveV0Tree, Bool_t lSaveCascadeTree, Bool_t lSaveHyperTriton, const char *name, TString lExtraOptions)
: AliAnalysisTaskSE(name), fListHist(0), fListV0(0), fListCascade(0), fTreeEvent(0), fTreeV0(0), fTreeCascade(0), fTreeHyperTriton3Body(nullptr),fPIDResponse(0), fESDtrackCuts(0), fESDtrackCutsITSsa2010(0), fESDtrackCutsGlobal2015(0), fUtils(0), fEventCuts(), fRand(0),

//---> Flags controlling Event Tree output
fkSaveEventTree    ( kFALSE ), //no downscaling in this tree so far

//---> Flags controlling V0 TTree output
fkSaveV0Tree       ( kTRUE ),
fkDownScaleV0      ( kTRUE  ),
fDownScaleFactorV0 ( 0.001  ),
fkPreselectDedx ( kFALSE ),
fkPreselectPID  ( kTRUE  ),
fkUseOnTheFlyV0Cascading( kFALSE ),
fkDoImprovedCascadeVertexFinding(kFALSE),
fkDoImprovedDCAV0DauPropagation(kFALSE),
fkDoImprovedDCACascDauPropagation(kFALSE),
fkIfImprovedPerformInitialLinearPropag( kFALSE ),
fkIfImprovedExtraPrecisionFactor ( 1.0 ),
fkDebugWrongPIDForTracking ( kFALSE ), //also for cascades...
fkDebugBump( kFALSE ),
fkDebugOOBPileup(kFALSE),
fkDoExtraEvSels( kTRUE ),

//---> Flags controlling Cascade TTree output
fkSaveCascadeTree       ( kTRUE  ),
fkDownScaleCascade      ( kTRUE  ),
fDownScaleFactorCascade ( 0.001  ),
fMinPtToSave( 0.00   ) ,
fMaxPtToSave( 100.00 ) ,

//---> Flags controlling HyperTriton3Body TTree output
fkSaveHyperTriton3BodyTree(lSaveHyperTriton),
fkDownScaleHyperTriton3Body(false),
fDownScaleFactorHyperTriton3Body(1.),

//---> Flags controlling Vertexers
fkRunVertexers    ( kFALSE ),
fkUseLightVertexer ( kTRUE ),
fkDoV0Refit ( kTRUE ),
fkExtraCleanup    ( kTRUE ),
fkSaveGoodTracks( kTRUE ),
fkSandboxV0( kTRUE ),
fkSandboxCascade( kFALSE ),

//---> Flag controlling trigger selection
fTrigType(AliVEvent::kMB),

fPrecisionCutoffCascadeDCA(1e-4),
fMaxIterationsCascadeDCA(27),

//---> Variables for fTreeEvent
fCentrality(0),
fMVPileupFlag(kFALSE),
fOOBPileupFlag(kFALSE),
fNTOFClusters(-1),
fNTOFMatches(-1),
fNTracksITSsa2010(-1),
fNTracksGlobal2015(-1),
fNTracksGlobal2015TriggerPP(-1),
fAmplitudeV0A(-1.),
fAmplitudeV0C(-1.),
fNHitsFMDA(-1.),
fNHitsFMDC(-1.),

//---> Variables for fTreeV0
fTreeVariableGoodV0(kFALSE),
fTreeVariableCentrality(0),
fTreeVariablePosLength(0),
fTreeVariableNegLength(0),
fTreeVariablePosCrossedRows(0),
fTreeVariableNegCrossedRows(0),
fTreeVariablePosTrackStatus(0),
fTreeVariableNegTrackStatus(0),
fTreeVariablePosDxy(0),
fTreeVariableNegDxy(0),
fTreeVariablePosDz(0),
fTreeVariableNegDz(0),
fTreeVariableDcaV0Daughters(0),
fTreeVariableDcaV0DaughtersGeometric(0),
fTreeVariablePosPropagStatus(0),
fTreeVariableNegPropagStatus(0),
fTreeVariableV0Radius(0),
fTreeVariableV0CosineOfPointingAngle(0),
fTreeVariableDecayX(0),
fTreeVariableDecayY(0),
fTreeVariableDecayZ(0),
fTreeVariableDecayXMC(0),
fTreeVariableDecayYMC(0),
fTreeVariableDecayZMC(0),
fTreeVariableNegPxMC(0),
fTreeVariableNegPyMC(0),
fTreeVariableNegPzMC(0),
fTreeVariablePosPxMC(0),
fTreeVariablePosPyMC(0),
fTreeVariablePosPzMC(0),
fTreeVariableInvMassK0s(0),
fTreeVariableInvMassLambda(0),
fTreeVariableInvMassAntiLambda(0),
fTreeVariablePID(0),
fTreeVariablePIDPositive(0),
fTreeVariablePIDNegative(0),
fTreeVariablePtMC(0),
fTreeVariableRapMC(0),
fTreeVariableNegTOFSignal(99999), 
fTreeVariablePosTOFSignal(99999),
fTreeVariablePosAlpha(0),
fTreeVariablePosSigmaY2(0),
fTreeVariablePosSigmaZ2(0),
fTreeVariableNegAlpha(0),
fTreeVariableNegSigmaY2(0),
fTreeVariableNegSigmaZ2(0),

fTreeVariablePosTrack(0x0),
fTreeVariableNegTrack(0x0),
fTreeVariableOTFV0(0x0),
fTreeVariableFoundOTFV0(kFALSE),
fTreeVariableMagneticField(0),
fTreeVariablePosOriginalX(0),
fTreeVariableNegOriginalX(0),
fTreeVariablePVx(0),
fTreeVariablePVy(0),
fTreeVariablePVz(0),
fTreeVariableAliESDvertex(0),
fTreeVariableRun(0),

//---> Variables for fTreeCascade
fTreeCascVarCentrality(0),
fTreeCascVarPosSign(0),
fTreeCascVarNegSign(0),
fTreeCascVarBachSign(0),
fTreeCascVarPosLength(0),
fTreeCascVarNegLength(0),
fTreeCascVarBachLength(0),
fTreeCascVarPosCrossedRows(0),
fTreeCascVarNegCrossedRows(0),
fTreeCascVarBachCrossedRows(0),
//Tracking flags
fTreeCascVarPosTrackStatus(0),
fTreeCascVarNegTrackStatus(0),
fTreeCascVarBachTrackStatus(0),
//DCAxy to PV
fTreeCascVarPosDxy(0),
fTreeCascVarNegDxy(0),
fTreeCascVarBachDxy(0),
//DCAz
fTreeCascVarPosDz(0),
fTreeCascVarNegDz(0),
fTreeCascVarBachDz(0),
fTreeCascVarDcaV0Daughters(0),
fTreeCascVarNegPropagStatus(0),
fTreeCascVarPosPropagStatus(0),
fTreeCascVarV0Radius(0),
fTreeCascVarV0DecayX(0),
fTreeCascVarV0DecayY(0),
fTreeCascVarV0DecayZ(0),
fTreeCascVarV0DecayXMC(0),
fTreeCascVarV0DecayYMC(0),
fTreeCascVarV0DecayZMC(0),
fTreeCascVarV0CosineOfPointingAngle(0),
fTreeCascVarDCAV0ToPrimVtx(0),
fTreeCascVarDCAxyV0ToPrimVtx(0),
fTreeCascVarInvMassLambda(0),
fTreeCascVarInvMassAntiLambda(0),

fTreeCascVarDCACascDaughters(0),
fTreeCascVarCascPropagation(0),

fTreeCascVarDecayX(0),
fTreeCascVarDecayY(0),
fTreeCascVarDecayZ(0),
fTreeCascVarDecayXMC(0),
fTreeCascVarDecayYMC(0),
fTreeCascVarDecayZMC(0),
fTreeCascVarCascCosPointingAngle(0),

fTreeCascVarInvMassXiMinus(0),
fTreeCascVarInvMassXiPlus(0),
fTreeCascVarInvMassOmegaMinus(0),
fTreeCascVarInvMassOmegaPlus(0),

fTreeCascVarPIDPositive(0),
fTreeCascVarPIDNegative(0),
fTreeCascVarPIDBachelor(0),
fTreeCascVarPID(0),
fTreeCascVarPtMC(0),
fTreeCascVarRapMC(0),

fTreeCascVarNegTOFSignal(99999),
fTreeCascVarPosTOFSignal(99999),
fTreeCascVarBachTOFSignal(99999),

fTreeCascVarPosDistanceToTrueDecayPt(0),
fTreeCascVarNegDistanceToTrueDecayPt(0),
fTreeCascVarBachDistanceToTrueDecayPt(0),
fTreeCascVarV0DistanceToTrueDecayPt(0),

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
fTreeCascVarNegPxMC(0),
fTreeCascVarNegPyMC(0),
fTreeCascVarNegPzMC(0),
fTreeCascVarPosPxMC(0),
fTreeCascVarPosPyMC(0),
fTreeCascVarPosPzMC(0),
fTreeCascVarBachPxMC(0),
fTreeCascVarBachPyMC(0),
fTreeCascVarBachPzMC(0),

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//Save full info for full re-vertex offline replay ('sandbox mode')
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
fTreeCascVarBachTrack(0),
fTreeCascVarPosTrack(0),
fTreeCascVarNegTrack(0),
fTreeCascVarOTFV0(0),
fTreeCascVarOTFV0NegBach(0x0),
fTreeCascVarOTFV0PosBach(0x0),
fTreeCascVarV0AsOTF(0),
fTreeCascVarNegBachAsOTF(0),
fTreeCascVarPosBachAsOTF(0),
fTreeCascVarMagneticField(0),
fTreeCascVarBachOriginalX(0),
fTreeCascVarPosOriginalX(0),
fTreeCascVarNegOriginalX(0),
fTreeCascVarPVx(0),
fTreeCascVarPVy(0),
fTreeCascVarPVz(0),
fTreeCascVarAliESDvertex(0),

/// Hypertriton 3 Body sandbox
fTreeHyp3BodyVarTracks(),
fTreeHyp3BodyVarPDGcodes(),
fTreeHyp3BodyVarEventId(0u),
fTreeHyp3BodyVarMotherId(0),
fTreeHyp3BodyVarTruePx(0.f),
fTreeHyp3BodyVarTruePy(0.f),
fTreeHyp3BodyVarTruePz(0.f),
fTreeHyp3BodyVarDecayVx(0.f),
fTreeHyp3BodyVarDecayVy(0.f),
fTreeHyp3BodyVarDecayVz(0.f),
fTreeHyp3BodyVarDecayT(0.f),
fTreeHyp3BodyVarPVx(0.f),
fTreeHyp3BodyVarPVy(0.f),
fTreeHyp3BodyVarPVz(0.f),
fTreeHyp3BodyVarPVt(0.f),
fTreeHyp3BodyVarMagneticField(0.f),

//Histos
fHistEventCounter(0),
fHistCentrality(0),
//V0s
fHistGeneratedPtVsYVsCentralityK0Short(0),
fHistGeneratedPtVsYVsCentralityLambda(0),
fHistGeneratedPtVsYVsCentralityAntiLambda(0),
//Cascades
fHistGeneratedPtVsYVsCentralityXiMinus(0),
fHistGeneratedPtVsYVsCentralityXiPlus(0),
fHistGeneratedPtVsYVsCentralityOmegaMinus(0),
fHistGeneratedPtVsYVsCentralityOmegaPlus(0),
//Hypertriton
fHistGeneratedPtVsYVsCentralityHypTrit(0),
fHistGeneratedPtVsYVsCentralityAntiHypTrit(0)

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
    DefineOutput(4, TTree::Class()); // Event Tree output
    DefineOutput(5, TTree::Class()); // V0 Tree output
    DefineOutput(6, TTree::Class()); // Cascade Tree output
    DefineOutput(7, TTree::Class()); // HyperTriton Tree output
    
    //Special Debug Options (more to be added as needed)
    // A - Study Wrong PID for tracking bug
    // B - Study invariant mass *B*ump
    // C - Study OOB pileup in pp 2016 data
    if ( lExtraOptions.Contains("A") ) fkDebugWrongPIDForTracking = kTRUE;
    if ( lExtraOptions.Contains("B") ) fkDebugBump                = kTRUE;
    if ( lExtraOptions.Contains("C") ) fkDebugOOBPileup           = kTRUE;
    
}

AliAnalysisTaskStrEffStudy::~AliAnalysisTaskStrEffStudy()
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
    if (fTreeHyperTriton3Body) {
        delete fTreeHyperTriton3Body;
        fTreeHyperTriton3Body = nullptr;
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
void AliAnalysisTaskStrEffStudy::UserCreateOutputObjects()
{
    //------------------------------------------------
    // fTreeEvent: EbyE information
    //------------------------------------------------
    fTreeEvent = new TTree("fTreeEvent","Event");
    //Branch Definitions
    fTreeEvent->Branch("fCentrality",&fCentrality,"fCentrality/F");
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
        fTreeEvent->Branch("fNHitsFMDA",&fNHitsFMDA,"fNHitsFMDA/F");
        fTreeEvent->Branch("fNHitsFMDC",&fNHitsFMDC,"fNHitsFMDC/F");
    }
    //------------------------------------------------
    // fTreeV0: V0 Candidate Information
    //------------------------------------------------
    
    //Create Basic V0 Output Tree
    fTreeV0 = new TTree( "fTreeV0", "Findable V0 Candidates");
    //-----------BASIC-INFO---------------------------
    fTreeV0->Branch("fTreeVariableGoodV0",&fTreeVariableGoodV0,"fTreeVariableGoodV0/O");
    fTreeV0->Branch("fTreeVariableCentrality",&fTreeVariableCentrality,"fTreeVariableCentrality/F");
    fTreeV0->Branch("fTreeVariablePosLength",&fTreeVariablePosLength,"fTreeVariablePosLength/F");
    fTreeV0->Branch("fTreeVariableNegLength",&fTreeVariableNegLength,"fTreeVariableNegLength/F");
    fTreeV0->Branch("fTreeVariablePosCrossedRows",&fTreeVariablePosCrossedRows,"fTreeVariablePosCrossedRows/F");
    fTreeV0->Branch("fTreeVariableNegCrossedRows",&fTreeVariableNegCrossedRows,"fTreeVariableNegCrossedRows/F");
    fTreeV0->Branch("fTreeVariableNegTrackStatus",&fTreeVariableNegTrackStatus,"fTreeVariableNegTrackStatus/l");
    fTreeV0->Branch("fTreeVariablePosTrackStatus",&fTreeVariablePosTrackStatus,"fTreeVariablePosTrackStatus/l");
    fTreeV0->Branch("fTreeVariablePosDxy",&fTreeVariablePosDxy,"fTreeVariablePosDxy/F");
    fTreeV0->Branch("fTreeVariableNegDxy",&fTreeVariableNegDxy,"fTreeVariableNegDxy/F");
    fTreeV0->Branch("fTreeVariablePosDz",&fTreeVariablePosDz,"fTreeVariablePosDz/F");
    fTreeV0->Branch("fTreeVariableNegDz",&fTreeVariableNegDz,"fTreeVariableNegDz/F");
    fTreeV0->Branch("fTreeVariableDcaV0Daughters",&fTreeVariableDcaV0Daughters,"fTreeVariableDcaV0Daughters/F");
    fTreeV0->Branch("fTreeVariableDcaV0DaughtersGeometric",&fTreeVariableDcaV0DaughtersGeometric,"fTreeVariableDcaV0DaughtersGeometric/F");
    fTreeV0->Branch("fTreeVariablePosPropagStatus",&fTreeVariablePosPropagStatus,"fTreeVariablePosPropagStatus/O");
    fTreeV0->Branch("fTreeVariableNegPropagStatus",&fTreeVariableNegPropagStatus,"fTreeVariableNegPropagStatus/O");
    fTreeV0->Branch("fTreeVariableV0Radius",&fTreeVariableV0Radius,"fTreeVariableV0Radius/F");
    fTreeV0->Branch("fTreeVariableV0CosineOfPointingAngle",&fTreeVariableV0CosineOfPointingAngle,"fTreeVariableV0CosineOfPointingAngle/F");
    fTreeV0->Branch("fTreeVariableDecayX",&fTreeVariableDecayX,"fTreeVariableDecayX/F");
    fTreeV0->Branch("fTreeVariableDecayY",&fTreeVariableDecayY,"fTreeVariableDecayY/F");
    fTreeV0->Branch("fTreeVariableDecayZ",&fTreeVariableDecayZ,"fTreeVariableDecayZ/F");
    fTreeV0->Branch("fTreeVariableDecayXMC",&fTreeVariableDecayXMC,"fTreeVariableDecayXMC/F");
    fTreeV0->Branch("fTreeVariableDecayYMC",&fTreeVariableDecayYMC,"fTreeVariableDecayYMC/F");
    fTreeV0->Branch("fTreeVariableDecayZMC",&fTreeVariableDecayZMC,"fTreeVariableDecayZMC/F");
    
    fTreeV0->Branch("fTreeVariableNegPxMC",&fTreeVariableNegPxMC,"fTreeVariableNegPxMC/F");
    fTreeV0->Branch("fTreeVariableNegPyMC",&fTreeVariableNegPyMC,"fTreeVariableNegPyMC/F");
    fTreeV0->Branch("fTreeVariableNegPzMC",&fTreeVariableNegPzMC,"fTreeVariableNegPzMC/F");
    fTreeV0->Branch("fTreeVariablePosPxMC",&fTreeVariablePosPxMC,"fTreeVariablePosPxMC/F");
    fTreeV0->Branch("fTreeVariablePosPyMC",&fTreeVariablePosPyMC,"fTreeVariablePosPyMC/F");
    fTreeV0->Branch("fTreeVariablePosPzMC",&fTreeVariablePosPzMC,"fTreeVariablePosPzMC/F");
    
    fTreeV0->Branch("fTreeVariableInvMassK0s",       &fTreeVariableInvMassK0s,       "fTreeVariableInvMassK0s/F");
    fTreeV0->Branch("fTreeVariableInvMassLambda",    &fTreeVariableInvMassLambda,    "fTreeVariableInvMassLambda/F");
    fTreeV0->Branch("fTreeVariableInvMassAntiLambda",&fTreeVariableInvMassAntiLambda,"fTreeVariableInvMassAntiLambda/F");
    fTreeV0->Branch("fTreeVariablePID",&fTreeVariablePID,"fTreeVariablePID/I");
    fTreeV0->Branch("fTreeVariablePIDPositive",&fTreeVariablePIDPositive,"fTreeVariablePIDPositive/I");
    fTreeV0->Branch("fTreeVariablePIDNegative",&fTreeVariablePIDNegative,"fTreeVariablePIDNegative/I");
    fTreeV0->Branch("fTreeVariablePtMC",&fTreeVariablePtMC,"fTreeVariablePtMC/F");
    fTreeV0->Branch("fTreeVariableRapMC",&fTreeVariableRapMC,"fTreeVariableRapMC/F");
    fTreeV0->Branch("fTreeVariableNegTOFSignal",&fTreeVariableNegTOFSignal,"fTreeVariableNegTOFSignal/F");
    fTreeV0->Branch("fTreeVariablePosTOFSignal",&fTreeVariablePosTOFSignal,"fTreeVariablePosTOFSignal/F");
    
    //Uncertainties
    fTreeV0->Branch("fTreeVariablePosAlpha",&fTreeVariablePosAlpha,"fTreeVariablePosAlpha/F");
    fTreeV0->Branch("fTreeVariablePosSigmaY2",&fTreeVariablePosSigmaY2,"fTreeVariablePosSigmaY2/F");
    fTreeV0->Branch("fTreeVariablePosSigmaZ2",&fTreeVariablePosSigmaZ2,"fTreeVariablePosSigmaZ2/F");
    fTreeV0->Branch("fTreeVariableNegAlpha",&fTreeVariableNegAlpha,"fTreeVariableNegAlpha/F");
    fTreeV0->Branch("fTreeVariableNegSigmaY2",&fTreeVariableNegSigmaY2,"fTreeVariableNegSigmaY2/F");
    fTreeV0->Branch("fTreeVariableNegSigmaZ2",&fTreeVariableNegSigmaZ2,"fTreeVariableNegSigmaZ2/F");
    
    //Sandbox mode
    if(fkSandboxV0){
        fTreeV0->Branch("fTreeVariablePosTrack", &fTreeVariablePosTrack,16000,99);
        fTreeV0->Branch("fTreeVariableNegTrack", &fTreeVariableNegTrack,16000,99);
        fTreeV0->Branch("fTreeVariableOTFV0", &fTreeVariableOTFV0,16000,99);
        fTreeV0->Branch("fTreeVariableFoundOTFV0", &fTreeVariableFoundOTFV0,"fTreeVariableFoundOTFV0/O");
        fTreeV0->Branch("fTreeVariableMagneticField",&fTreeVariableMagneticField,"fTreeVariableMagneticField/F");
        fTreeV0->Branch("fTreeVariablePosOriginalX",&fTreeVariablePosOriginalX,"fTreeVariablePosOriginalX/F");
        fTreeV0->Branch("fTreeVariableNegOriginalX",&fTreeVariableNegOriginalX,"fTreeVariableNegOriginalX/F");
        fTreeV0->Branch("fTreeVariablePVx",&fTreeVariablePVx,"fTreeVariablePVx/F");
        fTreeV0->Branch("fTreeVariablePVy",&fTreeVariablePVy,"fTreeVariablePVy/F");
        fTreeV0->Branch("fTreeVariablePVz",&fTreeVariablePVz,"fTreeVariablePVz/F");
        fTreeV0->Branch("fTreeVariableAliESDvertex", &fTreeVariableAliESDvertex,16000,99);
    }
    
    fTreeV0->Branch("fTreeVariableRun",&fTreeVariableRun,"fTreeVariableRun/I");
    //------------------------------------------------
    
    //------------------------------------------------
    // fTreeCascade Branch definitions - Cascade Tree
    //------------------------------------------------
    //Create Cascade output tree
    fTreeCascade = new TTree("fTreeCascade","CascadeCandidates");
    fTreeCascade->Branch("fTreeCascVarCentrality",&fTreeCascVarCentrality,"fTreeCascVarCentrality/F");
    //-----------BASIC-INFO---------------------------
    fTreeCascade->Branch("fTreeCascVarPosSign",&fTreeCascVarPosSign,"fTreeCascVarPosSign/I");
    fTreeCascade->Branch("fTreeCascVarNegSign",&fTreeCascVarNegSign,"fTreeCascVarNegSign/I");
    fTreeCascade->Branch("fTreeCascVarBachSign",&fTreeCascVarBachSign,"fTreeCascVarBachSign/I");
    fTreeCascade->Branch("fTreeCascVarPosLength",&fTreeCascVarPosLength,"fTreeCascVarPosLength/F");
    fTreeCascade->Branch("fTreeCascVarNegLength",&fTreeCascVarNegLength,"fTreeCascVarNegLength/F");
    fTreeCascade->Branch("fTreeCascVarBachLength",&fTreeCascVarBachLength,"fTreeCascVarBachLength/F");
    fTreeCascade->Branch("fTreeCascVarPosCrossedRows",&fTreeCascVarPosCrossedRows,"fTreeCascVarPosCrossedRows/F");
    fTreeCascade->Branch("fTreeCascVarNegCrossedRows",&fTreeCascVarNegCrossedRows,"fTreeCascVarNegCrossedRows/F");
    fTreeCascade->Branch("fTreeCascVarBachCrossedRows",&fTreeCascVarBachCrossedRows,"fTreeCascVarBachCrossedRows/F");
    fTreeCascade->Branch("fTreeCascVarPosTrackStatus",&fTreeCascVarPosTrackStatus,"fTreeCascVarPosTrackStatus/l");
    fTreeCascade->Branch("fTreeCascVarNegTrackStatus",&fTreeCascVarNegTrackStatus,"fTreeCascVarNegTrackStatus/l");
    fTreeCascade->Branch("fTreeCascVarBachTrackStatus",&fTreeCascVarBachTrackStatus,"fTreeCascVarBachTrackStatus/l");
    //DCAxy to PV
    fTreeCascade->Branch("fTreeCascVarPosDxy",&fTreeCascVarPosDxy,"fTreeCascVarPosDxy/F");
    fTreeCascade->Branch("fTreeCascVarNegDxy",&fTreeCascVarNegDxy,"fTreeCascVarNegDxy/F");
    fTreeCascade->Branch("fTreeCascVarBachDxy",&fTreeCascVarBachDxy,"fTreeCascVarBachDxy/F");
    //DCAz
    fTreeCascade->Branch("fTreeCascVarPosDz",&fTreeCascVarPosDz,"fTreeCascVarPosDz/F");
    fTreeCascade->Branch("fTreeCascVarNegDz",&fTreeCascVarNegDz,"fTreeCascVarNegDz/F");
    fTreeCascade->Branch("fTreeCascVarBachDz",&fTreeCascVarBachDz,"fTreeCascVarBachDz/F");
    
    fTreeCascade->Branch("fTreeCascVarDcaV0Daughters",&fTreeCascVarDcaV0Daughters,"fTreeCascVarDcaV0Daughters/F");
    fTreeCascade->Branch("fTreeCascVarNegPropagStatus",&fTreeCascVarNegPropagStatus,"fTreeCascVarNegPropagStatus/O");
    fTreeCascade->Branch("fTreeCascVarPosPropagStatus",&fTreeCascVarPosPropagStatus,"fTreeCascVarPosPropagStatus/O");
    fTreeCascade->Branch("fTreeCascVarDcaV0Daughters",&fTreeCascVarDcaV0Daughters,"fTreeCascVarDcaV0Daughters/F");
    fTreeCascade->Branch("fTreeCascVarV0DecayX",&fTreeCascVarV0DecayX,"fTreeCascVarV0DecayX/F");
    fTreeCascade->Branch("fTreeCascVarV0DecayY",&fTreeCascVarV0DecayY,"fTreeCascVarV0DecayY/F");
    fTreeCascade->Branch("fTreeCascVarV0DecayZ",&fTreeCascVarV0DecayZ,"fTreeCascVarV0DecayZ/F");
    fTreeCascade->Branch("fTreeCascVarV0DecayXMC",&fTreeCascVarV0DecayXMC,"fTreeCascVarV0DecayXMC/F");
    fTreeCascade->Branch("fTreeCascVarV0DecayYMC",&fTreeCascVarV0DecayYMC,"fTreeCascVarV0DecayYMC/F");
    fTreeCascade->Branch("fTreeCascVarV0DecayZMC",&fTreeCascVarV0DecayZMC,"fTreeCascVarV0DecayZMC/F");
    fTreeCascade->Branch("fTreeCascVarV0CosineOfPointingAngle",&fTreeCascVarV0CosineOfPointingAngle,"fTreeCascVarV0CosineOfPointingAngle/F");
    fTreeCascade->Branch("fTreeCascVarDCAV0ToPrimVtx",&fTreeCascVarDCAV0ToPrimVtx,"fTreeCascVarDCAV0ToPrimVtx/F");
    fTreeCascade->Branch("fTreeCascVarDCAxyV0ToPrimVtx",&fTreeCascVarDCAxyV0ToPrimVtx,"fTreeCascVarDCAxyV0ToPrimVtx/F");
    fTreeCascade->Branch("fTreeCascVarInvMassLambda",&fTreeCascVarInvMassLambda,"fTreeCascVarInvMassLambda/F");
    fTreeCascade->Branch("fTreeCascVarInvMassAntiLambda",&fTreeCascVarInvMassAntiLambda,"fTreeCascVarInvMassAntiLambda/F");
    
    //CLASSICAL BACK-PROPAGATION VARIABLES
    fTreeCascade->Branch("fTreeCascVarDCACascDaughters",&fTreeCascVarDCACascDaughters,"fTreeCascVarDCACascDaughters/F");
    fTreeCascade->Branch("fTreeCascVarCascPropagation",&fTreeCascVarCascPropagation,"fTreeCascVarCascPropagation/O");
    
    fTreeCascade->Branch("fTreeCascVarDecayX",&fTreeCascVarDecayX,"fTreeCascVarDecayX/F");
    fTreeCascade->Branch("fTreeCascVarDecayY",&fTreeCascVarDecayY,"fTreeCascVarDecayY/F");
    fTreeCascade->Branch("fTreeCascVarDecayZ",&fTreeCascVarDecayZ,"fTreeCascVarDecayZ/F");
    fTreeCascade->Branch("fTreeCascVarDecayXMC",&fTreeCascVarDecayXMC,"fTreeCascVarDecayXMC/F");
    fTreeCascade->Branch("fTreeCascVarDecayYMC",&fTreeCascVarDecayYMC,"fTreeCascVarDecayYMC/F");
    fTreeCascade->Branch("fTreeCascVarDecayZMC",&fTreeCascVarDecayZMC,"fTreeCascVarDecayZMC/F");
    fTreeCascade->Branch("fTreeCascVarCascCosPointingAngle",&fTreeCascVarCascCosPointingAngle,"fTreeCascVarCascCosPointingAngle/F");
    
    fTreeCascade->Branch("fTreeCascVarInvMassXiMinus",&fTreeCascVarInvMassXiMinus,"fTreeCascVarInvMassXiMinus/F");
    fTreeCascade->Branch("fTreeCascVarInvMassXiPlus",&fTreeCascVarInvMassXiPlus,"fTreeCascVarInvMassXiPlus/F");
    fTreeCascade->Branch("fTreeCascVarInvMassOmegaMinus",&fTreeCascVarInvMassOmegaMinus,"fTreeCascVarInvMassOmegaMinus/F");
    fTreeCascade->Branch("fTreeCascVarInvMassOmegaPlus",&fTreeCascVarInvMassOmegaPlus,"fTreeCascVarInvMassOmegaPlus/F");
    
    //MC VARIABLES
    fTreeCascade->Branch("fTreeCascVarPIDPositive",&fTreeCascVarPIDPositive,"fTreeCascVarPIDPositive/I");
    fTreeCascade->Branch("fTreeCascVarPIDNegative",&fTreeCascVarPIDNegative,"fTreeCascVarPIDNegative/I");
    fTreeCascade->Branch("fTreeCascVarPIDBachelor",&fTreeCascVarPIDBachelor,"fTreeCascVarPIDBachelor/I");
    fTreeCascade->Branch("fTreeCascVarPID",&fTreeCascVarPID,"fTreeCascVarPID/I");
    fTreeCascade->Branch("fTreeCascVarPtMC",&fTreeCascVarPtMC,"fTreeCascVarPtMC/F");
    fTreeCascade->Branch("fTreeCascVarRapMC",&fTreeCascVarRapMC,"fTreeCascVarRapMC/F");

    //TOF SIGNAL
    fTreeCascade->Branch("fTreeCascVarNegTOFSignal",&fTreeCascVarNegTOFSignal,"fTreeCascVarNegTOFSignal/F");
    fTreeCascade->Branch("fTreeCascVarPosTOFSignal",&fTreeCascVarPosTOFSignal,"fTreeCascVarPosTOFSignal/F");
    fTreeCascade->Branch("fTreeCascVarBachTOFSignal",&fTreeCascVarBachTOFSignal,"fTreeCascVarBachTOFSignal/F");
    
    fTreeCascade->Branch("fTreeCascVarPosDistanceToTrueDecayPt",&fTreeCascVarPosDistanceToTrueDecayPt,"fTreeCascVarPosDistanceToTrueDecayPt/F");
    fTreeCascade->Branch("fTreeCascVarNegDistanceToTrueDecayPt",&fTreeCascVarNegDistanceToTrueDecayPt,"fTreeCascVarNegDistanceToTrueDecayPt/F");
    fTreeCascade->Branch("fTreeCascVarBachDistanceToTrueDecayPt",&fTreeCascVarBachDistanceToTrueDecayPt,"fTreeCascVarBachDistanceToTrueDecayPt/F");
    fTreeCascade->Branch("fTreeCascVarV0DistanceToTrueDecayPt",&fTreeCascVarV0DistanceToTrueDecayPt,"fTreeCascVarV0DistanceToTrueDecayPt/F");
    
    //full momentum info
    fTreeCascade->Branch("fTreeCascVarPosPx",&fTreeCascVarPosPx,"fTreeCascVarPosPx/F");
    fTreeCascade->Branch("fTreeCascVarPosPy",&fTreeCascVarPosPy,"fTreeCascVarPosPy/F");
    fTreeCascade->Branch("fTreeCascVarPosPz",&fTreeCascVarPosPz,"fTreeCascVarPosPz/F");
    fTreeCascade->Branch("fTreeCascVarNegPx",&fTreeCascVarNegPx,"fTreeCascVarNegPx/F");
    fTreeCascade->Branch("fTreeCascVarNegPy",&fTreeCascVarNegPy,"fTreeCascVarNegPy/F");
    fTreeCascade->Branch("fTreeCascVarNegPz",&fTreeCascVarNegPz,"fTreeCascVarNegPz/F");
    fTreeCascade->Branch("fTreeCascVarBachPx",&fTreeCascVarBachPx,"fTreeCascVarBachPx/F");
    fTreeCascade->Branch("fTreeCascVarBachPy",&fTreeCascVarBachPy,"fTreeCascVarBachPy/F");
    fTreeCascade->Branch("fTreeCascVarBachPz",&fTreeCascVarBachPz,"fTreeCascVarBachPz/F");
    fTreeCascade->Branch("fTreeCascVarPosPxMC",&fTreeCascVarPosPxMC,"fTreeCascVarPosPxMC/F");
    fTreeCascade->Branch("fTreeCascVarPosPyMC",&fTreeCascVarPosPyMC,"fTreeCascVarPosPyMC/F");
    fTreeCascade->Branch("fTreeCascVarPosPzMC",&fTreeCascVarPosPzMC,"fTreeCascVarPosPzMC/F");
    fTreeCascade->Branch("fTreeCascVarNegPxMC",&fTreeCascVarNegPxMC,"fTreeCascVarNegPxMC/F");
    fTreeCascade->Branch("fTreeCascVarNegPyMC",&fTreeCascVarNegPyMC,"fTreeCascVarNegPyMC/F");
    fTreeCascade->Branch("fTreeCascVarNegPzMC",&fTreeCascVarNegPzMC,"fTreeCascVarNegPzMC/F");
    fTreeCascade->Branch("fTreeCascVarBachPxMC",&fTreeCascVarBachPxMC,"fTreeCascVarBachPxMC/F");
    fTreeCascade->Branch("fTreeCascVarBachPyMC",&fTreeCascVarBachPyMC,"fTreeCascVarBachPyMC/F");
    fTreeCascade->Branch("fTreeCascVarBachPzMC",&fTreeCascVarBachPzMC,"fTreeCascVarBachPzMC/F");
    
    if( fkSandboxCascade ){
        //Full track info for DCA minim optimization
        fTreeCascade->Branch("fTreeCascVarBachTrack", &fTreeCascVarBachTrack,16000,99);
        fTreeCascade->Branch("fTreeCascVarPosTrack", &fTreeCascVarPosTrack,16000,99);
        fTreeCascade->Branch("fTreeCascVarNegTrack", &fTreeCascVarNegTrack,16000,99);
        fTreeCascade->Branch("fTreeCascVarOTFV0", &fTreeCascVarOTFV0,16000,99);
        fTreeCascade->Branch("fTreeCascVarOTFV0NegBach", &fTreeCascVarOTFV0NegBach,16000,99);
        fTreeCascade->Branch("fTreeCascVarOTFV0PosBach", &fTreeCascVarOTFV0PosBach,16000,99);
        
        //for sandbox mode
        fTreeCascade->Branch("fTreeCascVarMagneticField",&fTreeCascVarMagneticField,"fTreeCascVarMagneticField/F");
        
        fTreeCascade->Branch("fTreeCascVarBachOriginalX",&fTreeCascVarBachOriginalX,"fTreeCascVarBachOriginalX/F");
        fTreeCascade->Branch("fTreeCascVarPosOriginalX",&fTreeCascVarPosOriginalX,"fTreeCascVarPosOriginalX/F");
        fTreeCascade->Branch("fTreeCascVarNegOriginalX",&fTreeCascVarNegOriginalX,"fTreeCascVarNegOriginalX/F");
        
        fTreeCascade->Branch("fTreeCascVarPVx",&fTreeCascVarPVx,"fTreeCascVarPVx/F");
        fTreeCascade->Branch("fTreeCascVarPVy",&fTreeCascVarPVy,"fTreeCascVarPVy/F");
        fTreeCascade->Branch("fTreeCascVarPVz",&fTreeCascVarPVz,"fTreeCascVarPVz/F");
        
        fTreeCascade->Branch("fTreeCascVarAliESDvertex", &fTreeCascVarAliESDvertex,16000,99);
    }
    fTreeCascade->Branch("fTreeCascVarV0AsOTF",&fTreeCascVarV0AsOTF,"fTreeCascVarV0AsOTF/O");
    fTreeCascade->Branch("fTreeCascVarNegBachAsOTF",&fTreeCascVarNegBachAsOTF,"fTreeCascVarNegBachAsOTF/O");
    fTreeCascade->Branch("fTreeCascVarPosBachAsOTF",&fTreeCascVarPosBachAsOTF,"fTreeCascVarPosBachAsOTF/O");

    fTreeHyperTriton3Body = new TTree("fTreeHyperTriton3Body","HyperTriton3BodyCandidates");
    if (fkSaveHyperTriton3BodyTree) {
        fTreeHyperTriton3Body->Branch("fTreeHyp3BodyVarTrack0", &fTreeHyp3BodyVarTracks[0],16000,99);
        fTreeHyperTriton3Body->Branch("fTreeHyp3BodyVarTrack1", &fTreeHyp3BodyVarTracks[1],16000,99);
        fTreeHyperTriton3Body->Branch("fTreeHyp3BodyVarTrack2", &fTreeHyp3BodyVarTracks[2],16000,99);

        fTreeHyperTriton3Body->Branch("fTreeHyp3BodyVarPDGcode0", &fTreeHyp3BodyVarPDGcodes[0],"fTreeHyp3BodyVarPDGcode0/I");
        fTreeHyperTriton3Body->Branch("fTreeHyp3BodyVarPDGcode1", &fTreeHyp3BodyVarPDGcodes[1],"fTreeHyp3BodyVarPDGcode1/I");
        fTreeHyperTriton3Body->Branch("fTreeHyp3BodyVarPDGcode2", &fTreeHyp3BodyVarPDGcodes[2],"fTreeHyp3BodyVarPDGcode2/I");

        fTreeHyperTriton3Body->Branch("fTreeHyp3BodyVarEventId",&fTreeHyp3BodyVarEventId,"fTreeHyp3BodyVarEventId/l");
        fTreeHyperTriton3Body->Branch("fTreeHyp3BodyVarMotherId",&fTreeHyp3BodyVarMotherId,"fTreeHyp3BodyVarMotherId/I");

        fTreeHyperTriton3Body->Branch("fTreeHyp3BodyVarTruePx",&fTreeHyp3BodyVarTruePx,"fTreeHyp3BodyVarTruePx/F");
        fTreeHyperTriton3Body->Branch("fTreeHyp3BodyVarTruePy",&fTreeHyp3BodyVarTruePy,"fTreeHyp3BodyVarTruePy/F");
        fTreeHyperTriton3Body->Branch("fTreeHyp3BodyVarTruePz",&fTreeHyp3BodyVarTruePz,"fTreeHyp3BodyVarTruePz/F");

        fTreeHyperTriton3Body->Branch("fTreeHyp3BodyVarDecayVx",&fTreeHyp3BodyVarDecayVx,"fTreeHyp3BodyVarDecayVx/F");
        fTreeHyperTriton3Body->Branch("fTreeHyp3BodyVarDecayVy",&fTreeHyp3BodyVarDecayVy,"fTreeHyp3BodyVarDecayVy/F");
        fTreeHyperTriton3Body->Branch("fTreeHyp3BodyVarDecayVz",&fTreeHyp3BodyVarDecayVz,"fTreeHyp3BodyVarDecayVz/F");
        fTreeHyperTriton3Body->Branch("fTreeHyp3BodyVarDecayT",&fTreeHyp3BodyVarDecayT,"fTreeHyp3BodyVarDecayT/F");

        fTreeHyperTriton3Body->Branch("fTreeHyp3BodyVarPVx",&fTreeHyp3BodyVarPVx,"fTreeHyp3BodyVarPVx/F");
        fTreeHyperTriton3Body->Branch("fTreeHyp3BodyVarPVy",&fTreeHyp3BodyVarPVy,"fTreeHyp3BodyVarPVy/F");
        fTreeHyperTriton3Body->Branch("fTreeHyp3BodyVarPVz",&fTreeHyp3BodyVarPVz,"fTreeHyp3BodyVarPVz/F");
        fTreeHyperTriton3Body->Branch("fTreeHyp3BodyVarPVt",&fTreeHyp3BodyVarPVt,"fTreeHyp3BodyVarPVt/F");

        fTreeHyperTriton3Body->Branch("fTreeHyp3BodyVarMagneticField",&fTreeHyp3BodyVarMagneticField,"fTreeHyp3BodyVarMagneticField/F");
    }
    //------------------------------------------------

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
        // used to fill the first 8 integers of the seed array.
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
    
    if(! fHistCentrality ) {
        //Histogram Output: Event-by-Event
        fHistCentrality = new TH1D( "fHistCentrality", "WARNING: no pileup rejection applied!;Centrality;Event Count",100,0,100);
        fListHist->Add(fHistCentrality);
    }
    
    if(! fHistGeneratedPtVsYVsCentralityK0Short ) {
        //Histogram Output: Efficiency Denominator
        fHistGeneratedPtVsYVsCentralityK0Short = new TH3D( "fHistGeneratedPtVsYVsCentralityK0Short", ";pT;y;centrality",500,0,25,40,-1.0,1.0,100,0,100);
        fListHist->Add(fHistGeneratedPtVsYVsCentralityK0Short);
    }
    if(! fHistGeneratedPtVsYVsCentralityLambda ) {
        //Histogram Output: Efficiency Denominator
        fHistGeneratedPtVsYVsCentralityLambda = new TH3D( "fHistGeneratedPtVsYVsCentralityLambda", ";pT;y;centrality",500,0,25,40,-1.0,1.0,100,0,100);
        fListHist->Add(fHistGeneratedPtVsYVsCentralityLambda);
    }
    if(! fHistGeneratedPtVsYVsCentralityAntiLambda ) {
        //Histogram Output: Efficiency Denominator
        fHistGeneratedPtVsYVsCentralityAntiLambda = new TH3D( "fHistGeneratedPtVsYVsCentralityAntiLambda", ";pT;y;centrality",500,0,25,40,-1.0,1.0,100,0,100);
        fListHist->Add(fHistGeneratedPtVsYVsCentralityAntiLambda);
    }
    if(! fHistGeneratedPtVsYVsCentralityXiMinus ) {
        //Histogram Output: Efficiency Denominator
        fHistGeneratedPtVsYVsCentralityXiMinus = new TH3D( "fHistGeneratedPtVsYVsCentralityXiMinus", ";pT;y;centrality",500,0,25,40,-1.0,1.0,100,0,100);
        fListHist->Add(fHistGeneratedPtVsYVsCentralityXiMinus);
    }
    if(! fHistGeneratedPtVsYVsCentralityXiPlus ) {
        //Histogram Output: Efficiency Denominator
        fHistGeneratedPtVsYVsCentralityXiPlus = new TH3D( "fHistGeneratedPtVsYVsCentralityXiPlus", ";pT;y;centrality",500,0,25,40,-1.0,1.0,100,0,100);
        fListHist->Add(fHistGeneratedPtVsYVsCentralityXiPlus);
    }
    if(! fHistGeneratedPtVsYVsCentralityOmegaMinus ) {
        //Histogram Output: Efficiency Denominator
        fHistGeneratedPtVsYVsCentralityOmegaMinus = new TH3D( "fHistGeneratedPtVsYVsCentralityOmegaMinus", ";pT;y;centrality",500,0,25,40,-1.0,1.0,100,0,100);
        fListHist->Add(fHistGeneratedPtVsYVsCentralityOmegaMinus);
    }
    if(! fHistGeneratedPtVsYVsCentralityOmegaPlus ) {
        //Histogram Output: Efficiency Denominator
        fHistGeneratedPtVsYVsCentralityOmegaPlus = new TH3D( "fHistGeneratedPtVsYVsCentralityOmegaPlus", ";pT;y;centrality",500,0,25,40,-1.0,1.0,100,0,100);
        fListHist->Add(fHistGeneratedPtVsYVsCentralityOmegaPlus);
    }
    if(! fHistGeneratedPtVsYVsCentralityHypTrit ) {
        //Histogram Output: Efficiency Denominator
        fHistGeneratedPtVsYVsCentralityHypTrit = new TH3D( "fHistGeneratedPtVsYVsCentralityHypTrit", ";pT;y;centrality",500,0,25,40,-1.0,1.0,100,0,100);
        fListHist->Add(fHistGeneratedPtVsYVsCentralityHypTrit);
    }
    if(! fHistGeneratedPtVsYVsCentralityAntiHypTrit ) {
        //Histogram Output: Efficiency Denominator
        fHistGeneratedPtVsYVsCentralityAntiHypTrit = new TH3D( "fHistGeneratedPtVsYVsCentralityAntiHypTrit", ";pT;y;centrality",500,0,25,40,-1.0,1.0,100,0,100);
        fListHist->Add(fHistGeneratedPtVsYVsCentralityAntiHypTrit);
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
    PostData(4, fTreeEvent   );
    PostData(5, fTreeV0      );
    PostData(6, fTreeCascade );
    PostData(7, fTreeHyperTriton3Body);
    
}// end UserCreateOutputObjects


//________________________________________________________________________
void AliAnalysisTaskStrEffStudy::UserExec(Option_t *)
{
    // Main loop
    // Called for each event
    AliESDEvent *lESDevent = 0x0;
    AliMCEvent  *lMCevent  = 0x0;
    AliStack    *lMCstack  = 0x0;
    
    fTreeCascVarPVx = -100;
    fTreeCascVarPVy = -100;
    fTreeCascVarPVz = -100;
    
    // Connect to the InputEvent
    // After these lines, we should have an ESD/AOD event + the number of V0s in it.
    
    // Appropriate for ESD analysis!
    
    lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
    if (!lESDevent) {
        AliWarning("ERROR: lESDevent not available \n");
        return;
    }
    
    //=================================================
    // Monte Carlo-related information
    lMCevent = MCEvent();
    if (!lMCevent) {
        Printf("ERROR: Could not retrieve MC event \n");
        cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
        return;
    }
    lMCstack = lMCevent->Stack();
    if (!lMCstack) {
        Printf("ERROR: Could not retrieve MC stack \n");
        cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
        return;
    }
    //=================================================
    
    //Get VZERO Information for multiplicity later
    AliVVZERO* esdV0 = lESDevent->GetVZEROData();
    if (!esdV0) {
        AliError("AliVVZERO not available");
        return;
    }
    
    Double_t lMagneticField = -10;
    lMagneticField = lESDevent->GetMagneticField( );
    
    //sandbox mode
    fTreeCascVarMagneticField = lMagneticField;
    fTreeVariableMagneticField = lMagneticField;
    
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

    Double_t lBestPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
    lPrimaryBestESDVtx->GetXYZ( lBestPrimaryVtxPos );

    AliESDVertex lPVobject(*lPrimaryBestESDVtx), *lPVpointer=&lPVobject;
    fTreeVariableAliESDvertex = lPVpointer;

    AliESDVertex lPVobject2(*lPrimaryBestESDVtx), *lPVpointer2=&lPVobject2;
    fTreeCascVarAliESDvertex = lPVpointer2;

    fTreeVariableRun = lESDevent->GetRunNumber();

    //sandbox info
    fTreeCascVarPVx = lBestPrimaryVtxPos[0];
    fTreeCascVarPVy = lBestPrimaryVtxPos[1];
    fTreeCascVarPVz = lBestPrimaryVtxPos[2];
    fTreeVariablePVx = lBestPrimaryVtxPos[0];
    fTreeVariablePVy = lBestPrimaryVtxPos[1];
    fTreeVariablePVz = lBestPrimaryVtxPos[2];

    //------------------------------------------------
    // Multiplicity Information Acquistion
    //------------------------------------------------

    Float_t lPercentile = 500;
    Float_t lPercentileEmbeddedSelection = 500;
    Int_t lEvSelCode = 100;
    AliMultSelection *MultSelection = (AliMultSelection*) lESDevent -> FindListObject("MultSelection");
    if( !MultSelection) {
        //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
        //AliWarning("AliMultSelection object not found! Trying to resort to AliCentrality now...");
        AliCentrality* centrality = 0x0;
        centrality = lESDevent->GetCentrality();
        if( centrality ){
            lPercentile = centrality->GetCentralityPercentileUnchecked("V0M");
            lPercentileEmbeddedSelection = lPercentile;
            lEvSelCode = 0;
            if(centrality->GetQuality()>1){
                //Not good!
                lEvSelCode = 999;
            }
        }
    } else {
        //V0M Multiplicity Percentile
        lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
        lPercentileEmbeddedSelection = MultSelection->GetMultiplicityPercentile("V0M", kTRUE );
        //Event Selection Code
        lEvSelCode = MultSelection->GetEvSelCode();
    }
    
    //just ask AliMultSelection. It will know.
    //fMVPileupFlag = kFALSE;
    //fMVPileupFlag = MultSelection->GetThisEventIsNotPileupMV();
    
    fCentrality = lPercentile;
    
    //Let's find out why efficiency is so centrality dependent, please!
    fTreeCascVarCentrality  = lPercentile;
    fTreeVariableCentrality = lPercentile;
    
    if( lEvSelCode != 0 ) {
        PostData(1, fListHist    );
        PostData(2, fListV0      );
        PostData(3, fListCascade );
        PostData(4, fTreeEvent   );
        PostData(5, fTreeV0      );
        PostData(6, fTreeCascade );
        PostData(7, fTreeHyperTriton3Body );
        return;
    }
    
    AliVEvent *ev = InputEvent();
    if( fkDoExtraEvSels ) {
        if( !fEventCuts.AcceptEvent(ev) ) {
            PostData(1, fListHist    );
            PostData(2, fListV0      );
            PostData(3, fListCascade );
            PostData(4, fTreeEvent   );
            PostData(5, fTreeV0      );
            PostData(6, fTreeCascade );
            PostData(7, fTreeHyperTriton3Body );
            return;
        }
    }
    
    fHistEventCounter->Fill(1.5);
    
    //Fill centrality histogram
    fHistCentrality->Fill(fCentrality);
    
    //Event-level fill
    if ( fkSaveEventTree ) fTreeEvent->Fill() ;
    
    //STOP HERE if skipping event selections (no point in doing the rest...)
    
    //------------------------------------------------
    
    //----- Loop on Generated Particles --------------
    Int_t    lThisPDG  = 0;
    Double_t lThisRap  = 0;
    Double_t lThisPt   = 0;
    for (Int_t ilab = 0;  ilab < (lMCstack->GetNtrack()); ilab++)
    {   // This is the begining of the loop on tracks
        
        TParticle* lPart = 0x0;
        lPart = lMCstack->Particle( ilab );
        if(!lPart) {
            Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", ilab );
            continue;
        }
        
        lThisPDG = lPart->GetPdgCode();
        
        //This if is necessary in some situations (rapidity calculation and PYTHIA junctions, etc)
        if ( (TMath::Abs(lThisPDG) == 3312) || (TMath::Abs(lThisPDG) == 3334) || (TMath::Abs(lThisPDG) == 3122) || lThisPDG == 310 || TMath::Abs(lThisPDG)==1010010030 )
        {
            lThisRap   = MyRapidity(lPart->Energy(),lPart->Pz());
            lThisPt    = lPart->Pt();
            
            //Use Physical Primaries only for filling These Histos
            if ( lMCstack->IsPhysicalPrimary(ilab)!=kTRUE ) continue;
            
            if( lThisPDG ==   310 ) {
                fHistGeneratedPtVsYVsCentralityK0Short       -> Fill (lThisPt, lThisRap, lPercentileEmbeddedSelection);
            }
            if( lThisPDG ==  3122 ) {
                fHistGeneratedPtVsYVsCentralityLambda       -> Fill (lThisPt, lThisRap, lPercentileEmbeddedSelection);
            }
            if( lThisPDG == -3122 ) {
                fHistGeneratedPtVsYVsCentralityAntiLambda       -> Fill (lThisPt, lThisRap, lPercentileEmbeddedSelection);
            }
            if( lThisPDG ==  3312 ) {
                fHistGeneratedPtVsYVsCentralityXiMinus       -> Fill (lThisPt, lThisRap, lPercentileEmbeddedSelection);
            }
            if( lThisPDG == -3312 ) {
                fHistGeneratedPtVsYVsCentralityXiPlus       -> Fill (lThisPt, lThisRap, lPercentileEmbeddedSelection);
            }
            if( lThisPDG ==  3334 ) {
                fHistGeneratedPtVsYVsCentralityOmegaMinus       -> Fill (lThisPt, lThisRap, lPercentileEmbeddedSelection);
            }
            if( lThisPDG == -3334 ) {
                fHistGeneratedPtVsYVsCentralityOmegaPlus       -> Fill (lThisPt, lThisRap, lPercentileEmbeddedSelection);
            }
            if( lThisPDG ==  1010010030 ) {
                fHistGeneratedPtVsYVsCentralityHypTrit       -> Fill (lThisPt, lThisRap, lPercentileEmbeddedSelection);
            }
            if( lThisPDG == -1010010030 ) {
                fHistGeneratedPtVsYVsCentralityAntiHypTrit       -> Fill (lThisPt, lThisRap, lPercentileEmbeddedSelection);
            }
        }
    }//End of loop on tracks
    //----- End Loop on Cascades ------------------------------------------------------------
    
    //------------------------------------------------
    // Fill V0 Tree as needed
    //------------------------------------------------
    
    //-------------------------------------------------
    // V0s from scratch: locate findable V0 candidates
    //-------------------------------------------------
    
    //Particles of interest
    constexpr Int_t lNV0Types = 5;
    Int_t lV0Types[lNV0Types]          = { 310, 3122, -3122,  1010010030, -1010010030};
    
    //Number of tracks
    Long_t lNTracks = lESDevent->GetNumberOfTracks();
    Double_t b      = lESDevent->GetMagneticField();
    
    //pos/neg daughters
    TArrayI lTrackArray      (lNTracks);
    TArrayI lTrackMotherArray(lNTracks);
    
    Long_t nTracksOfInterest = 0;
    
    auto lRemoveDeltaRayFromDaughters = [](const AliMCEvent* lMCev, const TParticle* lMum) {
        int lNDaughters = 0;
        for (int iPart = lMum->GetFirstDaughter(); iPart <= lMum->GetLastDaughter(); ++iPart) {
            TParticle* lParticle = lMCev->Particle(iPart);
            if (lParticle->GetPdgCode() != 11)
                lNDaughters++;
        }
        return lNDaughters;
    };

    //____________________________________________________________________________
    //Step 1: establish list of tracks coming from desired type
    for(Long_t iTrack = 0; iTrack<lNTracks; iTrack++){
        AliESDtrack *esdTrack = lESDevent->GetTrack(iTrack);
        if(!esdTrack) continue;
        
        Int_t lLabel = (Int_t) TMath::Abs( esdTrack->GetLabel() );
        TParticle* lParticle = lMCstack->Particle( lLabel );
        
        Int_t lLabelMother = lParticle->GetFirstMother();
        if( lLabelMother < 0 ) continue;
        
        //Do not select on primaries so that this list can be used for cascades too
        //if( lMCstack->IsPhysicalPrimary(lLabelMother) ) continue;
        
        TParticle *lParticleMother = lMCstack->Particle( lLabelMother );
        Int_t lParticleMotherPDG = lParticleMother->GetPdgCode();
        
        //Skip three-body decays and the like
        int lMotherType = -1;
        for(Int_t iType=0; iType<lNV0Types; iType++){
            if( lParticleMotherPDG == lV0Types[iType] ) lMotherType = iType;
        }
        if( lMotherType < 0 ) continue;

        int lNDaughters = lParticleMother->GetNDaughters();
        if ( lMotherType > 2 ) { ///Count real hypertriton daughters
            lNDaughters = lRemoveDeltaRayFromDaughters(lMCevent, lParticleMother);
        }

        if ( lNDaughters!=2 ) continue;
        
        //If here: this is a daughter of a mother particle of desired type, add
        lTrackArray        [nTracksOfInterest] = iTrack;
        lTrackMotherArray  [nTracksOfInterest] = lLabelMother;
        nTracksOfInterest++;
    }
    
    TArrayI lPosTrackArray      (lNTracks);
    TArrayI lNegTrackArray      (lNTracks);
    Long_t lFindableV0s = 0;
    
    //____________________________________________________________________________
    //Step 2: determine findable V0s: look for pairs having shared mother label
    for(Long_t iTrack = 0; iTrack<nTracksOfInterest; iTrack++){
        //Start nested loop from iTrack+1: avoid permutations + combination with self
        for(Long_t jTrack = iTrack+1; jTrack<nTracksOfInterest; jTrack++){
            if( lTrackMotherArray[iTrack]==lTrackMotherArray[jTrack]){
                //This is a findable V0! Yay! Check daughters before indexing
                AliESDtrack *esdTrack1 = 0x0;
                AliESDtrack *esdTrack2 = 0x0;
                esdTrack1 = lESDevent->GetTrack( lTrackArray[iTrack] );
                esdTrack2 = lESDevent->GetTrack( lTrackArray[jTrack] );
                
                //Check for non-existing
                if ( !esdTrack1 || !esdTrack2 ) continue;
                
                if( esdTrack1->GetSign() < 0 && esdTrack2->GetSign() < 0 ) continue;
                if( esdTrack1->GetSign() > 0 && esdTrack2->GetSign() > 0 ) continue;
                
                //1 = Positive, 2 = Negative case
                if( esdTrack1->GetSign() > 0 && esdTrack2->GetSign() < 0 ){
                    lPosTrackArray[lFindableV0s] = lTrackArray[iTrack];
                    lNegTrackArray[lFindableV0s] = lTrackArray[jTrack];
                    lFindableV0s++; //add this to findable
                }else{
                    lPosTrackArray[lFindableV0s] = lTrackArray[jTrack];
                    lNegTrackArray[lFindableV0s] = lTrackArray[iTrack];
                    lFindableV0s++; //add this to findable
                }
                
            }
        }
    }
    
    //____________________________________________________________________________
    //Step 3: Loop over findable V0s and save their relevant characteristics
    for(Long_t iV0 = 0; iV0<lFindableV0s; iV0++){
        //Get the two tracks we're talking about
        AliESDtrack *esdTrackPos = lESDevent->GetTrack( lPosTrackArray[iV0] );
        AliESDtrack *esdTrackNeg = lESDevent->GetTrack( lNegTrackArray[iV0] );
        
        fTreeVariableNegOriginalX = esdTrackNeg->GetX();
        fTreeVariablePosOriginalX = esdTrackPos->GetX();
        
        //store original esd tracks
        fTreeVariableNegTrack = esdTrackNeg;
        fTreeVariablePosTrack = esdTrackPos;
        
        //-----------------------------------------------------------------
        //3a: get basic track characteristics
        fTreeVariablePosLength = -1;
        fTreeVariableNegLength = -1;
        
        if( esdTrackPos->GetInnerParam() )
            fTreeVariablePosLength = esdTrackPos->GetLengthInActiveZone(1, 2.0, 220.0, lESDevent->GetMagneticField());
        if( esdTrackNeg->GetInnerParam() )
            fTreeVariableNegLength = esdTrackNeg->GetLengthInActiveZone(1, 2.0, 220.0, lESDevent->GetMagneticField());
        
        fTreeVariablePosCrossedRows = esdTrackPos ->GetTPCClusterInfo(2,1);
        fTreeVariableNegCrossedRows = esdTrackNeg ->GetTPCClusterInfo(2,1);
        //Tracking flags
        fTreeVariablePosTrackStatus = esdTrackPos->GetStatus();
        fTreeVariableNegTrackStatus = esdTrackNeg->GetStatus();
        //DCAxy to PV
        fTreeVariablePosDxy = TMath::Abs(esdTrackPos->GetD(lBestPrimaryVtxPos[0],
                                                           lBestPrimaryVtxPos[1],
                                                           lMagneticField) );
        fTreeVariableNegDxy = TMath::Abs(esdTrackNeg->GetD(lBestPrimaryVtxPos[0],
                                                           lBestPrimaryVtxPos[1],
                                                           lMagneticField) );
        //DCAz
        Float_t dz[2];
        esdTrackPos->GetDZ(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2], lMagneticField, dz );
        fTreeVariablePosDz = dz[1];
        esdTrackNeg->GetDZ(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2], lMagneticField, dz );
        fTreeVariableNegDz = dz[1];
        
        //-----------------------------------------------------------------
        //3b: Do V0 combination and see if it works, please
        //Step 1: propagate to DCA
        fTreeVariableDcaV0Daughters = -1;
        fTreeVariableDcaV0DaughtersGeometric = -1;
        Double_t xn, xp, dca; //=esdTrackNeg->GetDCA(esdTrackPos,lMagneticField,xn,xp);
        
        AliExternalTrackParam nt(*esdTrackNeg), pt(*esdTrackPos);//, *pointnt=&nt, *pointpt=&pt;
        dca=GetDCAV0Dau(&pt, &nt, xp, xn, lMagneticField);
        
        //Correct for beam pipe material
        //Warning: this is an unfinished implementation and should not do much at this stage
        
        Bool_t corrected=kFALSE;
        if ((nt.GetX() > 3.) && (xn < 3.)) {
            //correct for the beam pipe material
            corrected=kTRUE;
        }
        if ((pt.GetX() > 3.) && (xp < 3.)) {
            //correct for the beam pipe material
            corrected=kTRUE;
        }
        if (corrected) {
            //dca=nt.GetDCA(&pt,lMagneticField,xn,xp);
        }
        
        fTreeVariableDcaV0Daughters = dca; //Pass to TTree object, please
        
        //Actual propagation
        fTreeVariableNegPropagStatus = nt.PropagateTo(xn,lMagneticField);
        fTreeVariablePosPropagStatus = pt.PropagateTo(xp,lMagneticField);
        

        
        //=================================================================================
        //OTF loop: try to find equivalent OTF V0, store empty object if not found
        fTreeVariableOTFV0 = 0x0;
        
        //lNegTrackArray[iV0], lPosTrackArray[iV0]
        Int_t nv0s = lESDevent->GetNumberOfV0s();
        
        fTreeVariableFoundOTFV0 = kFALSE;
        
        for(Long_t iOTFv0=0; iOTFv0<nv0s; iOTFv0++){
            AliESDv0 *v0 = ((AliESDEvent*)lESDevent)->GetV0(iOTFv0);
            if (!v0) continue;
            
            if (v0->GetOnFlyStatus()  &&
                (
                 (v0->GetPindex() == lPosTrackArray[iV0] && v0->GetNindex() == lNegTrackArray[iV0] ) ||
                 (v0->GetNindex() == lPosTrackArray[iV0] && v0->GetPindex() == lNegTrackArray[iV0] )
                 )
                ){
                //Found corresponding OTF V0! Save it to TTree, please
                AliESDv0 lV0ToStore(*v0), *lPointerToV0ToStore=&lV0ToStore;
                fTreeVariableOTFV0 = lPointerToV0ToStore;
                fTreeVariableFoundOTFV0 = kTRUE;
                break; //stop looking
            }
        }
        if( !fTreeVariableFoundOTFV0 ) {
            AliESDv0 lV0ToStore, *lPointerToV0ToStore=&lV0ToStore;
            fTreeVariableOTFV0 = lPointerToV0ToStore;
        }
        //=================================================================================
        
        //Tag OK V0s (will probably tag >99%? will still have to be studied!)
        if ( fTreeVariableNegPropagStatus == kTRUE &&
            fTreeVariablePosPropagStatus == kTRUE )
            fTreeVariableGoodV0 = kTRUE;
        
        //Acquire the DCA that's not strictly computed with uncertainties (geometric only) for comparison
        Double_t lx1, ly1, lz1, lx2, ly2, lz2, tmp[3];
        nt.GetXYZ(tmp);
        lx1 = tmp[0]; ly1 = tmp[1]; lz1 = tmp[2];
        pt.GetXYZ(tmp);
        lx2 = tmp[0]; ly2 = tmp[1]; lz2 = tmp[2];
        fTreeVariableDcaV0DaughtersGeometric = TMath::Sqrt(
                                                           TMath::Power(lx1-lx2,2)+
                                                           TMath::Power(ly1-ly2,2)+
                                                           TMath::Power(lz1-lz2,2)
                                                           );
        
        //Get track uncertainties
        fTreeVariableNegAlpha   = nt.GetAlpha();
        fTreeVariableNegSigmaY2 = nt.GetSigmaY2();
        fTreeVariableNegSigmaZ2 = nt.GetSigmaZ2();
        fTreeVariablePosAlpha   = pt.GetAlpha();
        fTreeVariablePosSigmaY2 = pt.GetSigmaY2();
        fTreeVariablePosSigmaZ2 = pt.GetSigmaZ2();
        
        //Step 2: Attempt creation of a V0 vertex in these conditions
        AliESDv0 vertex(nt,lNegTrackArray[iV0],pt,lPosTrackArray[iV0]);
        
        //Get 2D decay radius from V0 vertex
        Double_t x=vertex.Xv(), y=vertex.Yv();
        Double_t r2D = TMath::Sqrt(x*x + y*y);
        fTreeVariableV0Radius = r2D;
        
        //Get Estimated decay position
        fTreeVariableDecayX = x;
        fTreeVariableDecayY = y;
        fTreeVariableDecayZ = vertex.Zv();
        
        //Get Cosine of pointing angle
        Float_t cpa=vertex.GetV0CosineOfPointingAngle(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2]);
        fTreeVariableV0CosineOfPointingAngle = cpa;

        //Get TOF signal
        fTreeVariableNegTOFSignal = esdTrackNeg->GetTOFsignal() * 1.e-3; // in ns
        fTreeVariablePosTOFSignal = esdTrackPos->GetTOFsignal() * 1.e-3; // in ns
        
        //Final step: get estimated masses under different mass hypotheses
        vertex.ChangeMassHypothesis(310);
        fTreeVariableInvMassK0s = vertex.GetEffMass();
        vertex.ChangeMassHypothesis(3122);
        fTreeVariableInvMassLambda = vertex.GetEffMass();
        vertex.ChangeMassHypothesis(-3122);
        fTreeVariableInvMassAntiLambda = vertex.GetEffMass();
        
        //-----------------------------------------------------------------
        //3c: Get perfect MC information for bookkeeping
        Int_t lblPosV0Dghter = (Int_t) TMath::Abs( esdTrackPos->GetLabel() );
        Int_t lblNegV0Dghter = (Int_t) TMath::Abs( esdTrackNeg->GetLabel() );
        
        TParticle* mcPosV0Dghter = lMCstack->Particle( lblPosV0Dghter );
        TParticle* mcNegV0Dghter = lMCstack->Particle( lblNegV0Dghter );
        
        //-----------------------------------------------------------------
        //3c: Get perfect MC information for bookkeeping
        fTreeVariableDecayXMC = mcPosV0Dghter->Vx();
        fTreeVariableDecayYMC = mcPosV0Dghter->Vy();
        fTreeVariableDecayZMC = mcPosV0Dghter->Vz();

        fTreeVariableNegPxMC = mcNegV0Dghter->Px();
        fTreeVariableNegPyMC = mcNegV0Dghter->Py();
        fTreeVariableNegPzMC = mcNegV0Dghter->Pz();
        fTreeVariablePosPxMC = mcPosV0Dghter->Px();
        fTreeVariablePosPyMC = mcPosV0Dghter->Py();
        fTreeVariablePosPzMC = mcPosV0Dghter->Pz();

        fTreeVariablePIDPositive = mcPosV0Dghter -> GetPdgCode();
        fTreeVariablePIDNegative = mcNegV0Dghter -> GetPdgCode();
        
        Int_t lblMotherV0 = mcPosV0Dghter->GetFirstMother();
        TParticle* pThisV0 = lMCstack->Particle( lblMotherV0 );
        //Set tree variables
        fTreeVariablePID   = pThisV0->GetPdgCode(); //PDG Code
        fTreeVariablePtMC  = pThisV0->Pt(); //Perfect Pt
        fTreeVariableRapMC = pThisV0->Y();
        
        //IMPORTANT: select only physical primaries, please
        if( ! lMCstack->IsPhysicalPrimary( lblMotherV0 ) ) continue; //won't fill TTree
        
        if( fkSaveGoodTracks ){
            //...where good -> kTPCrefit, at least length zero (more still needed!)
            if ((fTreeVariablePosTrackStatus&AliESDtrack::kTPCrefit)==0) continue;
            if ((fTreeVariableNegTrackStatus&AliESDtrack::kTPCrefit)==0) continue;
            if(fTreeVariablePosLength<0) continue;
            if(fTreeVariableNegLength<0) continue;
        }
        
        //End step 3: fill findable ttree
        if( fkSaveV0Tree ) fTreeV0->Fill();
    }
    
    //--] END V0 PART [-------------------------------
    
    //------------------------------------------------
    // Cascades from scratch: locate findable cascades
    //------------------------------------------------
    
    //Particles of interest
    Int_t lCascadeTypes[4] = {3312, -3312, 3334, -3334};
    
    //pos/neg daughters
    TArrayI lBachelorArray      (lNTracks);
    TArrayI lBachelorMotherArray(lNTracks);
    
    Long_t nBachelorsOfInterest = 0;
    
    //________________________________________________________________________fTreeCascade____
    //Step 1: establish list of bachelors from cascades
    for(Long_t iTrack = 0; iTrack<lNTracks; iTrack++){
        AliESDtrack *esdTrack = lESDevent->GetTrack(iTrack);
        if(!esdTrack) continue;
        
        Int_t lLabel = (Int_t) TMath::Abs( esdTrack->GetLabel() );
        TParticle* lParticle = lMCstack->Particle( lLabel );
        
        Int_t lLabelMother = lParticle->GetFirstMother();
        if( lLabelMother < 0 ) continue;
        
        //Only interested in tracks whose mother was a primary (cascade)
        if( !lMCstack->IsPhysicalPrimary(lLabelMother) ) continue;
        
        TParticle *lParticleMother = lMCstack->Particle( lLabelMother );
        Int_t lParticleMotherPDG = lParticleMother->GetPdgCode();
        
        //Skip three-body decays and the like (has to be bach+V0)
        if ( lParticleMother->GetNDaughters()!=2 ) continue;

        Bool_t lOfDesiredType = kFALSE;
        for(Int_t iType=0; iType<4; iType++){
            if( lParticleMotherPDG == lCascadeTypes[iType] ) lOfDesiredType = kTRUE;
        }
        if( !lOfDesiredType ) continue;

        //If here: this is a daughter of a mother particle of desired type, add
        lBachelorArray        [nBachelorsOfInterest] = iTrack;
        lBachelorMotherArray  [nBachelorsOfInterest] = lLabelMother;
        nBachelorsOfInterest++;
    }
    cout<<"Findable bachelors: "<<nBachelorsOfInterest<<endl;

    TArrayI lCascPosTrackArray      (lNTracks);
    TArrayI lCascNegTrackArray      (lNTracks);
    TArrayI lCascBachTrackArray     (lNTracks);
    Long_t lFindableCascades = 0;

    //____________________________________________________________________________
    //Step 2: Loop over findable V0s and check if they share a mother with bach
    for(Long_t iV0 = 0; iV0<lFindableV0s; iV0++){
        //Get the two tracks we're talking about
        AliESDtrack *esdTrackPos = lESDevent->GetTrack( lPosTrackArray[iV0] );
        
        //3c: Get perfect MC information for bookkeeping
        Int_t lblPosV0Dghter = (Int_t) TMath::Abs( esdTrackPos->GetLabel() );
        
        TParticle* mcPosV0Dghter = lMCstack->Particle( lblPosV0Dghter );
        
        Int_t lblMotherV0 = mcPosV0Dghter->GetFirstMother();
        TParticle* pThisV0 = lMCstack->Particle( lblMotherV0 );
        
        //Only interested in lambdas
        if( TMath::Abs( pThisV0 -> GetPdgCode() ) != 3122 ) continue;
        
        Int_t lblGrandMotherV0 = pThisV0->GetFirstMother();
        if( lblGrandMotherV0 < 0 ) continue; //Not a cascade
        
        TParticle* pV0Origin = lMCstack->Particle( lblGrandMotherV0 );
        Int_t pV0OriginPDG = pV0Origin->GetPdgCode();
        //Pre-filter on primaries
        if( ! lMCstack->IsPhysicalPrimary( lblGrandMotherV0 ) ) continue;
        
        //Check if this is actually from a real cascade
        Bool_t lOfDesiredType = kFALSE;
        for(Int_t iType=0; iType<4; iType++){
            if( pV0OriginPDG == lCascadeTypes[iType] ) lOfDesiredType = kTRUE;
        }
        if( !lOfDesiredType ) continue;
        
        //Nested loop (all filters pre-applied already), combine with bach
        for(Long_t iBach=0; iBach<nBachelorsOfInterest; iBach++){
            if( lblGrandMotherV0 == lBachelorMotherArray[iBach] ){
                //Interesting! This is a valid cascade. Please save the track indices
                lCascNegTrackArray [lFindableCascades] = lNegTrackArray[iV0];
                lCascPosTrackArray [lFindableCascades] = lPosTrackArray[iV0];
                lCascBachTrackArray[lFindableCascades] = lBachelorArray[iBach];
                lFindableCascades++;
            }
        }
    }
    cout<<"Findable Cascades: "<<lFindableCascades<<endl;
    
    //____________________________________________________________________________
    //Step 3: Loop over findable cascades and determine their relevant characteristics
    for(Long_t iCasc = 0; iCasc<lFindableCascades; iCasc++){
        //Get the three tracks we're talking about
        AliESDtrack *esdTrackPos  = lESDevent->GetTrack( lCascPosTrackArray[iCasc] );
        AliESDtrack *esdTrackNeg  = lESDevent->GetTrack( lCascNegTrackArray[iCasc] );
        AliESDtrack *esdTrackBach = lESDevent->GetTrack( lCascBachTrackArray[iCasc] );
        
        //Sandbox information: always, regardless of status
        fTreeCascVarBachTrack = esdTrackBach;
        fTreeCascVarPosTrack = esdTrackPos;
        fTreeCascVarNegTrack = esdTrackNeg;
        
        //=================================================================================
        //OTF loop: try to find equivalent OTF V0, store empty object if not found
        fTreeCascVarOTFV0NegBach = 0x0;
        fTreeCascVarOTFV0PosBach = 0x0;
        fTreeCascVarPosBachAsOTF = kFALSE;
        fTreeCascVarNegBachAsOTF = kFALSE;
        
        //lNegTrackArray[iV0], lPosTrackArray[iV0]
        Int_t nv0s = lESDevent->GetNumberOfV0s();
        
        AliESDv0 lV0ToStorePosBach, *lPointerToV0ToStorePosBach=&lV0ToStorePosBach;
        AliESDv0 lV0ToStoreNegBach, *lPointerToV0ToStoreNegBach=&lV0ToStoreNegBach;
        
        Bool_t lFoundOTF = kFALSE;
        for(Long_t iOTFv0=0; iOTFv0<nv0s; iOTFv0++){
            AliESDv0 *v0 = ((AliESDEvent*)lESDevent)->GetV0(iOTFv0);
            if (!v0) continue;
            if ( v0->GetOnFlyStatus()  &&
                (
                 (v0->GetPindex() == lCascPosTrackArray[iCasc] && v0->GetNindex() == lCascBachTrackArray[iCasc] ) ||
                 (v0->GetNindex() == lCascPosTrackArray[iCasc] && v0->GetPindex() == lCascBachTrackArray[iCasc] )
                 )
                ){
                //Found corresponding OTF V0! Save it to TTree, please
                fTreeCascVarOTFV0PosBach = v0;
                lFoundOTF = kTRUE;
                fTreeCascVarPosBachAsOTF = kTRUE;
                break; //stop looking
            }
        }
        if( !lFoundOTF ) {
            fTreeCascVarOTFV0PosBach = lPointerToV0ToStorePosBach;
        }
        //---->
        lFoundOTF = kFALSE;
        for(Long_t iOTFv0=0; iOTFv0<nv0s; iOTFv0++){
            AliESDv0 *v0 = ((AliESDEvent*)lESDevent)->GetV0(iOTFv0);
            if (!v0) continue;
            if ( v0->GetOnFlyStatus()  &&
                (
                 (v0->GetPindex() == lCascNegTrackArray[iCasc] && v0->GetNindex() == lCascBachTrackArray[iCasc] ) ||
                 (v0->GetNindex() == lCascNegTrackArray[iCasc] && v0->GetPindex() == lCascBachTrackArray[iCasc] )
                 )
                ){
                //Found corresponding OTF V0! Save it to TTree, please
                fTreeCascVarOTFV0NegBach = v0;
                lFoundOTF = kTRUE;
                fTreeCascVarNegBachAsOTF = kTRUE;
                break; //stop looking
            }
        }
        if( !lFoundOTF ) {
            fTreeCascVarOTFV0NegBach = lPointerToV0ToStoreNegBach;
        }
        
        //=================================================================================
        //OTF loop: try to find equivalent OTF V0, store empty object if not found
        fTreeCascVarOTFV0 = 0x0;
        
        lFoundOTF = kFALSE;
        AliESDv0 lV0ToStore, *lPointerToV0ToStore=&lV0ToStore;
        fTreeCascVarV0AsOTF = kFALSE;
        for(Long_t iOTFv0=0; iOTFv0<nv0s; iOTFv0++){
            AliESDv0 *v0 = ((AliESDEvent*)lESDevent)->GetV0(iOTFv0);
            if (!v0) continue;
            
            if ( v0->GetOnFlyStatus()  &&
                (
                 //Check if this is the bump -> case switch for positive hyperons
                 (v0->GetPindex() == lCascPosTrackArray[iCasc] && v0->GetNindex() == lCascNegTrackArray[iCasc] ) ||
                 (v0->GetNindex() == lCascPosTrackArray[iCasc] && v0->GetPindex() == lCascNegTrackArray[iCasc] )
                 )
                ){
                //Found corresponding OTF V0! Save it to TTree, please
                fTreeCascVarOTFV0 = v0;
                lFoundOTF = kTRUE;
                fTreeCascVarV0AsOTF = kTRUE;
                break; //stop looking
            }
        }
        if( !lFoundOTF ) {
            fTreeCascVarOTFV0 = lPointerToV0ToStore;
        }
        //=================================================================================
        
        
        //get original X values (sandbox mode)
        fTreeCascVarBachOriginalX = esdTrackBach->GetX();
        fTreeCascVarPosOriginalX  = esdTrackPos->GetX();
        fTreeCascVarNegOriginalX  = esdTrackNeg->GetX();
        
        fTreeCascVarNegPx = 0.0;
        fTreeCascVarNegPy = 0.0;
        fTreeCascVarNegPz = 0.0;
        fTreeCascVarPosPx = 0.0;
        fTreeCascVarPosPy = 0.0;
        fTreeCascVarPosPz = 0.0;
        fTreeCascVarBachPx = 0.0;
        fTreeCascVarBachPy = 0.0;
        fTreeCascVarBachPz = 0.0;
        fTreeCascVarNegPxMC = 0.0;
        fTreeCascVarNegPyMC = 0.0;
        fTreeCascVarNegPzMC = 0.0;
        fTreeCascVarPosPxMC = 0.0;
        fTreeCascVarPosPyMC = 0.0;
        fTreeCascVarPosPzMC = 0.0;
        fTreeCascVarBachPxMC = 0.0;
        fTreeCascVarBachPyMC = 0.0;
        fTreeCascVarBachPzMC = 0.0;
        
        //-----------------------------------------------------------------
        //3a: get basic track characteristics
        fTreeCascVarPosLength  = -1;
        fTreeCascVarNegLength  = -1;
        fTreeCascVarBachLength = -1;
        
        if( esdTrackPos->GetInnerParam() )
            fTreeCascVarPosLength = esdTrackPos->GetLengthInActiveZone(1, 2.0, 220.0, lESDevent->GetMagneticField());
        if( esdTrackNeg->GetInnerParam() )
            fTreeCascVarNegLength = esdTrackNeg->GetLengthInActiveZone(1, 2.0, 220.0, lESDevent->GetMagneticField());
        if( esdTrackBach->GetInnerParam() )
            fTreeCascVarBachLength = esdTrackBach->GetLengthInActiveZone(1, 2.0, 220.0, lESDevent->GetMagneticField());
        
        fTreeCascVarPosCrossedRows  = esdTrackPos  ->GetTPCClusterInfo(2,1);
        fTreeCascVarNegCrossedRows  = esdTrackNeg  ->GetTPCClusterInfo(2,1);
        fTreeCascVarBachCrossedRows = esdTrackBach ->GetTPCClusterInfo(2,1);
        //Tracking flags
        fTreeCascVarPosTrackStatus = esdTrackPos->GetStatus();
        fTreeCascVarNegTrackStatus = esdTrackNeg->GetStatus();
        fTreeCascVarBachTrackStatus = esdTrackBach->GetStatus();
        //Charge
        fTreeCascVarPosSign = esdTrackPos -> GetSign();
        fTreeCascVarNegSign = esdTrackNeg -> GetSign();
        fTreeCascVarBachSign = esdTrackBach -> GetSign();
        //DCAxy to PV
        fTreeCascVarPosDxy = TMath::Abs(esdTrackPos->GetD(lBestPrimaryVtxPos[0],
                                                          lBestPrimaryVtxPos[1],
                                                          lMagneticField) );
        fTreeCascVarNegDxy = TMath::Abs(esdTrackNeg->GetD(lBestPrimaryVtxPos[0],
                                                          lBestPrimaryVtxPos[1],
                                                          lMagneticField) );
        fTreeCascVarBachDxy = TMath::Abs(esdTrackBach->GetD(lBestPrimaryVtxPos[0],
                                                           lBestPrimaryVtxPos[1],
                                                           lMagneticField) );
        //DCAz
        Float_t dztrack[2];
        esdTrackPos->GetDZ(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2], lMagneticField, dztrack );
        fTreeCascVarPosDz = dztrack[1];
        esdTrackNeg->GetDZ(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2], lMagneticField, dztrack );
        fTreeCascVarNegDz = dztrack[1];
        esdTrackBach->GetDZ(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2], lMagneticField, dztrack );
        fTreeCascVarBachDz = dztrack[1];
        
        
        //-----------------------------------------------------------------
        //3b: Do V0+cascade combination and see if it works, please
        
        //---] V0 PART [---------------------------------------------------
        //Step 1: propagate to DCA
        fTreeCascVarDcaV0Daughters = -1;
        
        //Correct for beam pipe material
        //Warning: this is an unfinished implementation and should not do much at this stage
        AliExternalTrackParam nt(*esdTrackNeg), pt(*esdTrackPos), *ntp=&nt, *ptp=&pt;
        Double_t xn, xp, dca;
        
        if( fkDoImprovedDCAV0DauPropagation ){
            //Improved: use own call
            dca=GetDCAV0Dau(ptp, ntp, xp, xn, b);
        }else{
            //Old: use old call
            dca=nt.GetDCA(&pt,b,xn,xp);
        }
        
        fTreeCascVarDcaV0Daughters = dca; //Pass to TTree object, please
        
        //Actual propagation
        fTreeCascVarNegPropagStatus = nt.PropagateTo(xn,lMagneticField);
        fTreeCascVarPosPropagStatus = pt.PropagateTo(xp,lMagneticField);
        
        //Step 2: Attempt creation of a V0 vertex in these conditions
        AliESDv0 vertex(nt,lCascNegTrackArray[iCasc],pt,lCascPosTrackArray[iCasc]);
        
        //Get 2D decay radius from V0 vertex
        Double_t x=vertex.Xv(), y=vertex.Yv();
        Double_t r2D = TMath::Sqrt(x*x + y*y);
        fTreeCascVarV0Radius = r2D;
        
        //Get Estimated decay position
        fTreeCascVarV0DecayX = x;
        fTreeCascVarV0DecayY = y;
        fTreeCascVarV0DecayZ = vertex.Zv();
        
        //Get Cosine of pointing angle
        Float_t cpa=vertex.GetV0CosineOfPointingAngle(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2]);
        fTreeCascVarV0CosineOfPointingAngle = cpa;
        
        //DCA to PV
        fTreeCascVarDCAV0ToPrimVtx = vertex.GetD(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2]);
        fTreeCascVarDCAxyV0ToPrimVtx = vertex.GetD(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1]);

        //Get TOF signal
        fTreeCascVarBachTOFSignal = esdTrackBach->GetTOFsignal() * 1.e-3; // in ns
        fTreeCascVarNegTOFSignal  = esdTrackNeg->GetTOFsignal() * 1.e-3; // in ns
        fTreeCascVarPosTOFSignal  = esdTrackPos->GetTOFsignal() * 1.e-3; // in ns
        
        //Final step: get estimated masses under different mass hypotheses
        vertex.ChangeMassHypothesis(3122);
        fTreeCascVarInvMassLambda = vertex.GetEffMass();
        vertex.ChangeMassHypothesis(-3122);
        fTreeCascVarInvMassAntiLambda = vertex.GetEffMass();
        
        //---] Cascade PART [----------------------------------------------
        
        //Step 1: propagation (encapsulated in PropagateToDCA, see appropriate options)
        AliESDv0 v0(vertex);
        AliESDv0 *pv0=&v0;
        AliExternalTrackParam bt(*esdTrackBach), *pbt=&bt;
        Double_t cascdca = PropagateToDCA(pv0,pbt,lESDevent,lMagneticField);
        
        fTreeCascVarDCACascDaughters = 1e+10;
        fTreeCascVarCascPropagation = kFALSE;
        
        fTreeCascVarNegPx = -100;
        fTreeCascVarNegPy = -100;
        fTreeCascVarNegPz = -100;
        fTreeCascVarPosPx = -100;
        fTreeCascVarPosPy = -100;
        fTreeCascVarPosPz = -100;
        fTreeCascVarBachPx = -100;
        fTreeCascVarBachPy = -100;
        fTreeCascVarBachPz = -100;
        
        fTreeCascVarDecayX = -100;
        fTreeCascVarDecayY = -100;
        fTreeCascVarDecayZ = -100;
        fTreeCascVarCascCosPointingAngle = -100;
        
        fTreeCascVarInvMassXiMinus = -100;
        fTreeCascVarInvMassXiPlus = -100;
        fTreeCascVarInvMassOmegaMinus = -100;
        fTreeCascVarInvMassOmegaPlus = -100;
        
        //Check if propagation successful
        if (cascdca < 1e+4){
            fTreeCascVarDCACascDaughters = cascdca;
            fTreeCascVarCascPropagation = kTRUE;
            
            //Construct cascade
            AliESDcascade cascade(*pv0,*pbt,lCascBachTrackArray[iCasc]);
            
            //Decay Position
            Double_t xcasc,ycasc,zcasc; cascade.GetXYZcascade(xcasc,ycasc,zcasc);
            fTreeCascVarDecayX = xcasc;
            fTreeCascVarDecayY = ycasc;
            fTreeCascVarDecayZ = zcasc;
            
            fTreeCascVarCascCosPointingAngle = cascade.GetCascadeCosineOfPointingAngle(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2]);
            
            Double_t lV0quality  = 0.;
            cascade.ChangeMassHypothesis(lV0quality , 3312);
            fTreeCascVarInvMassXiMinus = cascade.GetEffMassXi();
            cascade.ChangeMassHypothesis(lV0quality ,-3312);
            fTreeCascVarInvMassXiPlus = cascade.GetEffMassXi();
            cascade.ChangeMassHypothesis(lV0quality , 3334);
            fTreeCascVarInvMassOmegaMinus = cascade.GetEffMassXi();
            cascade.ChangeMassHypothesis(lV0quality ,-3334);
            fTreeCascVarInvMassOmegaPlus = cascade.GetEffMassXi();
            
            Double_t lBMom[3], lNMom[3], lPMom[3];
            cascade.GetBPxPyPz( lBMom[0], lBMom[1], lBMom[2] );
            cascade.GetPPxPyPz( lPMom[0], lPMom[1], lPMom[2] );
            cascade.GetNPxPyPz( lNMom[0], lNMom[1], lNMom[2] );
            
            fTreeCascVarNegPx = lNMom[0];
            fTreeCascVarNegPy = lNMom[1];
            fTreeCascVarNegPz = lNMom[2];
            fTreeCascVarPosPx = lPMom[0];
            fTreeCascVarPosPy = lPMom[1];
            fTreeCascVarPosPz = lPMom[2];
            fTreeCascVarBachPx = lBMom[0];
            fTreeCascVarBachPy = lBMom[1];		
            fTreeCascVarBachPz = lBMom[2];
        }
        
        //-----------------------------------------------------------------
        //3c: Get perfect MC information for bookkeeping
        Int_t lblPosCascDghter = (Int_t) TMath::Abs( esdTrackPos->GetLabel() );
        Int_t lblNegCascDghter = (Int_t) TMath::Abs( esdTrackNeg->GetLabel() );
        Int_t lblBachCascDghter = (Int_t) TMath::Abs( esdTrackBach->GetLabel() );
        
        TParticle* mcPosCascDghter  = lMCstack->Particle( lblPosCascDghter );
        TParticle* mcNegCascDghter  = lMCstack->Particle( lblNegCascDghter );
        TParticle* mcBachCascDghter = lMCstack->Particle( lblBachCascDghter );
        
        //Get V0/bachelor decay position
        //Be careful: Vx, Vy, Vz: Creation vertex. So decay position is the
        //Creation vertex of any one of the daughters!
        fTreeCascVarV0DecayXMC = mcPosCascDghter->Vx();
        fTreeCascVarV0DecayYMC = mcPosCascDghter->Vy();
        fTreeCascVarV0DecayZMC = mcPosCascDghter->Vz();
        fTreeCascVarDecayXMC = mcBachCascDghter->Vx();
        fTreeCascVarDecayYMC = mcBachCascDghter->Vy();
        fTreeCascVarDecayZMC = mcBachCascDghter->Vz();
        
        //Get MC information
        fTreeCascVarNegPxMC = mcNegCascDghter->Px();
        fTreeCascVarNegPyMC = mcNegCascDghter->Py();
        fTreeCascVarNegPzMC = mcNegCascDghter->Pz();
        fTreeCascVarPosPxMC = mcPosCascDghter->Px();
        fTreeCascVarPosPyMC = mcPosCascDghter->Py();
        fTreeCascVarPosPzMC = mcPosCascDghter->Pz();
        fTreeCascVarBachPxMC = mcBachCascDghter->Px();
        fTreeCascVarBachPyMC = mcBachCascDghter->Py();
        fTreeCascVarBachPzMC = mcBachCascDghter->Pz();
        
        fTreeCascVarPIDPositive = mcPosCascDghter -> GetPdgCode();
        fTreeCascVarPIDNegative = mcNegCascDghter -> GetPdgCode();
        fTreeCascVarPIDBachelor = mcBachCascDghter -> GetPdgCode();
        
        Int_t lLabelMother = mcBachCascDghter->GetFirstMother();
        
        TParticle *lParticleMother = lMCstack->Particle( lLabelMother );
        //Int_t lParticleMotherPDG = lParticleMother->GetPdgCode();
        
        //Set tree variables
        fTreeCascVarPID   = lParticleMother->GetPdgCode(); //PDG Code
        fTreeCascVarPtMC  = lParticleMother->Pt(); //Perfect Pt
        fTreeCascVarRapMC = lParticleMother->Y();
        
        if( fkSaveGoodTracks ){
            //...where good -> kTPCrefit, at least length zero (more still needed!)
            if ((fTreeCascVarPosTrackStatus&AliESDtrack::kTPCrefit)==0) continue;
            if ((fTreeCascVarNegTrackStatus&AliESDtrack::kTPCrefit)==0) continue;
            if ((fTreeCascVarBachTrackStatus&AliESDtrack::kTPCrefit)==0) continue;
            if(fTreeCascVarPosLength<0) continue;
            if(fTreeCascVarNegLength<0) continue;
            if(fTreeCascVarBachLength<0) continue;
        }
        
        //Check how close the daughter tracks passed to the relevant decay points
        Float_t dzspec[2];
        esdTrackPos->GetDZ( fTreeCascVarV0DecayXMC, fTreeCascVarV0DecayYMC, fTreeCascVarV0DecayZMC, lMagneticField, dzspec );
        fTreeCascVarPosDistanceToTrueDecayPt = TMath::Sqrt(dzspec[0]*dzspec[0]+dzspec[1]*dzspec[1]);
        esdTrackNeg->GetDZ( fTreeCascVarV0DecayXMC, fTreeCascVarV0DecayYMC, fTreeCascVarV0DecayZMC, lMagneticField, dzspec );
        fTreeCascVarNegDistanceToTrueDecayPt = TMath::Sqrt(dzspec[0]*dzspec[0]+dzspec[1]*dzspec[1]);
        esdTrackBach->GetDZ( fTreeCascVarDecayXMC, fTreeCascVarDecayYMC, fTreeCascVarDecayZMC, lMagneticField, dzspec );
        fTreeCascVarBachDistanceToTrueDecayPt = TMath::Sqrt(dzspec[0]*dzspec[0]+dzspec[1]*dzspec[1]);
        
        //Check how close the reconstructed V0 passed to the actual cascade decay point
        fTreeCascVarV0DistanceToTrueDecayPt = vertex.GetD(fTreeCascVarDecayXMC,fTreeCascVarDecayYMC,fTreeCascVarDecayZMC);
        
        //Fill Findable cascade tree
        if( fkSaveCascadeTree ) fTreeCascade->Fill();
    }
    
    //--] END CASCADE PART [--------------------------

    //------------------------------------------------
    // HyperTriton in 3 body from scratch: locate findable HyperTriton
    //------------------------------------------------

    //pos/neg daughters
    std::vector<TrackMC> lTrackOfInterest;
    lTrackOfInterest.reserve(lNTracks);

    //_________________________________________________________
    //Step 1: establish list of tracks coming from des
    for(Long_t iTrack = 0; iTrack < lNTracks; iTrack++){
        AliESDtrack *esdTrack = lESDevent->GetTrack(iTrack);
        if (!esdTrack) continue;
        /// The minimal TPC/ITS reconstruction criteria must be statisfied
        if (((esdTrack->GetStatus() & AliVTrack::kTPCrefit) == 0 &&
             (esdTrack->GetStatus() & AliVTrack::kITSrefit) == 0) ||
            esdTrack->GetKinkIndex(0) > 0)
            continue;
        Int_t lLabel = (Int_t) TMath::Abs( esdTrack->GetLabel() );
        TParticle* lParticle = lMCevent->Particle( lLabel );
        const int pdgAbs = std::abs(lParticle->GetPdgCode());

        Int_t lLabelMother = lParticle->GetFirstMother();
        if (lLabelMother < 0) continue;

        if (!lMCevent->IsPhysicalPrimary(lLabelMother)) continue;

        TParticle *lParticleMother = lMCevent->Particle( lLabelMother );

        Int_t lParticleMotherPDG = lParticleMother->GetPdgCode();
        if (std::abs(lParticleMotherPDG) != 1010010030) continue;
        int lNDaughters = lRemoveDeltaRayFromDaughters(lMCevent, lParticleMother);
        if (lNDaughters!=3) continue;
        //If here: this is a daughter of a mother particle of desired type, add
        lTrackOfInterest.push_back({esdTrack, lParticleMother, lParticle, lLabelMother});
    }

    if (!lTrackOfInterest.empty()) {
        bool lNewEvent = true;
        fTreeHyp3BodyVarMagneticField = b;
        fTreeHyp3BodyVarEventId++;
        fTreeHyp3BodyVarPVt = lTrackOfInterest.back().mother->T();
        fTreeHyp3BodyVarPVx = lTrackOfInterest.back().mother->Vx();
        fTreeHyp3BodyVarPVy = lTrackOfInterest.back().mother->Vy();
        fTreeHyp3BodyVarPVz = lTrackOfInterest.back().mother->Vz();

        /// This makes the output tree sorted, having possible clones close to each other.
        /// At the same time this quick sort will speed up the following loops
        std::sort(lTrackOfInterest.begin(), lTrackOfInterest.end(), [](const TrackMC & a, const TrackMC & b)
        {
            return a.motherId > b.motherId;
        });
        //____________________________________________________________________________
        //Step 2: determine findable hypertritons
        std::vector<CandidateMC> candidate;
        for (size_t iTrack = 0; iTrack < lTrackOfInterest.size(); iTrack++) {
            std::array<std::pair<int,int>,3> index;
            Int_t pdg1 = lTrackOfInterest[iTrack].particle->GetPdgCode();
            index[0] = {pdg1, iTrack};
            //Start nested loop from iTrack+1: avoid permutations + combination with self
            for (size_t jTrack = iTrack+1; jTrack < lTrackOfInterest.size(); jTrack++) {
                if (lTrackOfInterest[iTrack].motherId != lTrackOfInterest[jTrack].motherId) continue;
                Int_t pdg2 = lTrackOfInterest[jTrack].particle->GetPdgCode();
                index[1] = {pdg2, jTrack};
                for (size_t zTrack = jTrack+1; zTrack < lTrackOfInterest.size(); zTrack++) {
                    if(lTrackOfInterest[iTrack].motherId != lTrackOfInterest[zTrack].motherId) continue;
                    /// Reject all the triplets with +++ and ---
                    if (lTrackOfInterest[iTrack].track->GetSign() == lTrackOfInterest[jTrack].track->GetSign() &&
                        lTrackOfInterest[iTrack].track->GetSign() == lTrackOfInterest[zTrack].track->GetSign())
                        continue;
                    Int_t pdg3 = lTrackOfInterest[zTrack].particle->GetPdgCode();
                    index[2] = {pdg3, zTrack};
                    std::sort(index.begin(),index.end(),[](const std::pair<int,int> & a, const std::pair<int,int> & b)
                    {
                        return std::abs(a.first) > std::abs(b.first);
                    });
                    CandidateMC c;
                    c.track_deu = lTrackOfInterest[index[0].second].track;
                    c.track_p   = lTrackOfInterest[index[1].second].track;
                    c.track_pi  = lTrackOfInterest[index[2].second].track;
                    c.part1     = lTrackOfInterest[index[0].second].particle;
                    c.part2     = lTrackOfInterest[index[1].second].particle;
                    c.part3     = lTrackOfInterest[index[2].second].particle;
                    c.mother    = lTrackOfInterest[index[0].second].mother;
                    c.motherId  = lTrackOfInterest[index[0].second].motherId;
                    candidate.push_back(c);
                }
            }
        }
        /// sorting hypertriton candidates respect the motherId
        std::sort(candidate.begin(), candidate.end(), [](const CandidateMC &a, const CandidateMC &b)
        { 
            return a.motherId > b.motherId; 
        });
        //____________________________________________________________________________
        // Step 3: checks on the candidate vector
        for (size_t iCand = 0; iCand < candidate.size(); iCand++) {

            fTreeHyp3BodyVarTracks[0] = candidate[iCand].track_deu;
            fTreeHyp3BodyVarTracks[1] = candidate[iCand].track_p;
            fTreeHyp3BodyVarTracks[2] = candidate[iCand].track_pi;
            fTreeHyp3BodyVarPDGcodes[0] = candidate[iCand].part1->GetPdgCode();
            fTreeHyp3BodyVarPDGcodes[1] = candidate[iCand].part2->GetPdgCode();
            fTreeHyp3BodyVarPDGcodes[2] = candidate[iCand].part3->GetPdgCode();

            TParticle* lHyperTriton = candidate[iCand].mother;
            fTreeHyp3BodyVarTruePx = lHyperTriton->Px();
            fTreeHyp3BodyVarTruePy = lHyperTriton->Py();
            fTreeHyp3BodyVarTruePz = lHyperTriton->Pz();

            TParticle* prong = candidate[iCand].part1;
            fTreeHyp3BodyVarDecayVx = prong->Vx();
            fTreeHyp3BodyVarDecayVy = prong->Vy();
            fTreeHyp3BodyVarDecayVz = prong->Vz();
            fTreeHyp3BodyVarDecayT =  prong->T();

            fTreeHyp3BodyVarMotherId = candidate[iCand].motherId;
            if (lNewEvent) {
                fTreeHyp3BodyVarMotherId *= -1;
                lNewEvent = false;
            }
            fTreeHyperTriton3Body->Fill();
        }
    }

    //--] END HYPERTRITON3BODY PART [--------------------------

    // Post output data.
    PostData(1, fListHist    );
    PostData(2, fListV0      );
    PostData(3, fListCascade );
    PostData(4, fTreeEvent   );
    PostData(5, fTreeV0      );
    PostData(6, fTreeCascade );
    PostData(7, fTreeHyperTriton3Body );
}

//________________________________________________________________________
void AliAnalysisTaskStrEffStudy::Terminate(Option_t *)
{
    // Draw result to the screen
    // Called once at the end of the query
    
    TList *cRetrievedList = 0x0;
    cRetrievedList = (TList*)GetOutputData(1);
    if(!cRetrievedList) {
        Printf("ERROR - AliAnalysisTaskStrEffStudy : ouput data container list not available\n");
        return;
    }
    
    fHistEventCounter = dynamic_cast<TH1D*> (  cRetrievedList->FindObject("fHistEventCounter")  );
    if (!fHistEventCounter) {
        Printf("ERROR - AliAnalysisTaskStrEffStudy : fHistEventCounter not available");
        return;
    }
    
    TCanvas *canCheck = new TCanvas("AliAnalysisTaskStrEffStudy","V0 Multiplicity",10,10,510,510);
    canCheck->cd(1)->SetLogy();
    
    fHistEventCounter->SetMarkerStyle(22);
    fHistEventCounter->DrawCopy("E");
}

//________________________________________________________________________
Double_t AliAnalysisTaskStrEffStudy::MyRapidity(Double_t rE, Double_t rPz) const
{
    // Local calculation for rapidity
    Double_t ReturnValue = -100;
    if( (rE-rPz+1.e-13) != 0 && (rE+rPz) != 0 ) {
        ReturnValue =  0.5*TMath::Log((rE+rPz)/(rE-rPz+1.e-13));
    }
    return ReturnValue;
}

//________________________________________________________________________
void AliAnalysisTaskStrEffStudy::AddConfiguration( AliV0Result *lV0Result )
{
    if (!fListV0){
        Printf("fListV0 does not exist. Creating...");
        fListV0 = new TList();
        fListV0->SetOwner();
        
    }
    fListV0->Add(lV0Result);
}

//________________________________________________________________________
void AliAnalysisTaskStrEffStudy::AddConfiguration( AliCascadeResult *lCascadeResult )
{
    if (!fListCascade){
        Printf("fListCascade does not exist. Creating...");
        fListCascade = new TList();
        fListCascade->SetOwner();
        
    }
    fListCascade->Add(lCascadeResult);
}

//________________________________________________________________________
void AliAnalysisTaskStrEffStudy::SetupStandardVertexing()
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
void AliAnalysisTaskStrEffStudy::SetupLooseVertexing()
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
void AliAnalysisTaskStrEffStudy::AddTopologicalQAV0(Int_t lRecNumberOfSteps)
//Add all configurations to do QA of topological variables for the V0 analysis
{
    //Deprecated! Use main analysis task!
}

//________________________________________________________________________
void AliAnalysisTaskStrEffStudy::AddTopologicalQACascade(Int_t lRecNumberOfSteps)
//Add all configurations to do QA of topological variables for the V0 analysis
{
    //Deprecated! Use main analysis task!
}

//________________________________________________________________________
void AliAnalysisTaskStrEffStudy::AddStandardV0Configuration()
//Meant to add some standard V0 analysis Configuration + its corresponding systematics
{
    //Deprecated! Use main analysis task!
}

//________________________________________________________________________
void AliAnalysisTaskStrEffStudy::AddStandardCascadeConfiguration()
//Meant to add some standard cascade analysis Configuration + its corresponding systematics
{
    //Deprecated! Use main analysis task!
}

//________________________________________________________________________
void AliAnalysisTaskStrEffStudy::AddCascadeConfiguration276TeV()
//Adds 2.76 TeV cascade analysis configuration
{
    //Deprecated! Use main analysis task!
}

//________________________________________________________________________
Float_t AliAnalysisTaskStrEffStudy::GetDCAz(AliESDtrack *lTrack)
//Encapsulation of DCAz calculation
{
    Float_t b[2];
    Float_t bCov[3];
    lTrack->GetImpactParameters(b,bCov);
    if (bCov[0]<=0 || bCov[2]<=0) {
        AliDebug(1, "Estimated b resolution lower or equal to zero!");
        bCov[0]=0; bCov[2]=0;
    }
    //Float_t dcaToVertexXY = b[0];
    Float_t dcaToVertexZ = b[1];
    
    return dcaToVertexZ;
}

//________________________________________________________________________
Float_t AliAnalysisTaskStrEffStudy::GetCosPA(AliESDtrack *lPosTrack, AliESDtrack *lNegTrack, AliESDEvent *lEvent)
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
    Double_t xn, xp;
    lNegClone->GetDCA(lPosClone,b,xn,xp);
    
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
void AliAnalysisTaskStrEffStudy::CheckChargeV0(AliESDv0 *v0)
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

//______________________________________________________________________
AliAnalysisTaskStrEffStudy::FMDhits AliAnalysisTaskStrEffStudy::GetFMDhits(AliAODEvent* aodEvent) const
// Relies on the event being vaild (no extra checks if object exists done here)
{
    AliAODForwardMult* aodForward = static_cast<AliAODForwardMult*>(aodEvent->FindListObject("Forward"));
    // Shape of d2Ndetadphi: 200, -4, 6, 20, 0, 2pi
    const TH2D& d2Ndetadphi = aodForward->GetHistogram();
    Int_t nEta = d2Ndetadphi.GetXaxis()->GetNbins();
    Int_t nPhi = d2Ndetadphi.GetYaxis()->GetNbins();
    FMDhits ret_vector;
    for (Int_t iEta = 1; iEta <= nEta; iEta++) {
        Int_t valid = Int_t(d2Ndetadphi.GetBinContent(iEta, 0));
        if (!valid) {
            // No data expected for this eta
            continue;
        }
        Float_t eta = d2Ndetadphi.GetXaxis()->GetBinCenter(iEta);
        for (Int_t iPhi = 1; iPhi <= nPhi; iPhi++) {
            // Bin content is most likely number of particles!
            Float_t mostProbableN = d2Ndetadphi.GetBinContent(iEta, iPhi);
            if (mostProbableN > 0) {
                Float_t phi = d2Ndetadphi.GetYaxis()->GetBinCenter(iPhi);
                ret_vector.push_back(AliAnalysisTaskStrEffStudy::FMDhit(eta, phi, mostProbableN));
            }
        }
    }
    return ret_vector;
}

//________________________________________________________________________
Double_t AliAnalysisTaskStrEffStudy::Det(Double_t a00, Double_t a01, Double_t a10, Double_t a11) const {
    //--------------------------------------------------------------------
    // This function calculates locally a 2x2 determinant
    //--------------------------------------------------------------------
    return a00*a11 - a01*a10;
}

//________________________________________________________________________
Double_t AliAnalysisTaskStrEffStudy::Det(Double_t a00,Double_t a01,Double_t a02,
                                         Double_t a10,Double_t a11,Double_t a12,
                                         Double_t a20,Double_t a21,Double_t a22) const {
    //--------------------------------------------------------------------
    // This function calculates locally a 3x3 determinant
    //--------------------------------------------------------------------
    return  a00*Det(a11,a12,a21,a22)-a01*Det(a10,a12,a20,a22)+a02*Det(a10,a11,a20,a21);
}

//________________________________________________________________________
Double_t AliAnalysisTaskStrEffStudy::PropagateToDCA(AliESDv0 *v, AliExternalTrackParam *t, AliESDEvent *event, Double_t b) {
    //--------------------------------------------------------------------
    // This function returns the DCA between the V0 and the track
    //--------------------------------------------------------------------
    
    //Count received
    //fHistV0ToBachelorPropagationStatus->Fill(0.5);
    
    Double_t alpha=t->GetAlpha(), cs1=TMath::Cos(alpha), sn1=TMath::Sin(alpha);
    Double_t r[3]; t->GetXYZ(r);
    Double_t x1=r[0], y1=r[1], z1=r[2];
    Double_t p[3]; t->GetPxPyPz(p);
    Double_t px1=p[0], py1=p[1], pz1=p[2];
    
    Double_t x2,y2,z2;     // position and momentum of V0
    Double_t px2,py2,pz2;
    
    v->GetXYZ(x2,y2,z2);
    v->GetPxPyPz(px2,py2,pz2);
    
    Double_t dca = 1e+33;
    if ( !fkDoImprovedCascadeVertexFinding || fkIfImprovedPerformInitialLinearPropag ){
        // calculation dca
        Double_t dd= Det(x2-x1,y2-y1,z2-z1,px1,py1,pz1,px2,py2,pz2);
        Double_t ax= Det(py1,pz1,py2,pz2);
        Double_t ay=-Det(px1,pz1,px2,pz2);
        Double_t az= Det(px1,py1,px2,py2);
        
        dca=TMath::Abs(dd)/TMath::Sqrt(ax*ax + ay*ay + az*az);
        
        //points of the DCA
        Double_t t1 = Det(x2-x1,y2-y1,z2-z1,px2,py2,pz2,ax,ay,az)/
        Det(px1,py1,pz1,px2,py2,pz2,ax,ay,az);
        
        x1 += px1*t1; y1 += py1*t1; //z1 += pz1*t1;
        
        //propagate track to the points of DCA
        
        x1=x1*cs1 + y1*sn1;
        if (!t->PropagateTo(x1,b)) {
            //Count linear propagation failures
            //fHistV0ToBachelorPropagationStatus->Fill(1.5);
            Error("PropagateToDCA","Propagation failed !");
            return 1.e+33;
        }
        //Count linear propagation successes
        //fHistV0ToBachelorPropagationStatus->Fill(2.5);
    }
    
    if( fkDoImprovedCascadeVertexFinding ){
        //Count Improved Cascade propagation received
        //fHistV0ToBachelorPropagationStatus->Fill(3.5); //bin 4
        
        //DCA Calculation improved -> non-linear propagation
        //Preparatory step 1: get two tracks corresponding to V0
        UInt_t lKeyPos = (UInt_t)TMath::Abs(v->GetPindex());
        UInt_t lKeyNeg = (UInt_t)TMath::Abs(v->GetNindex());
        AliESDtrack *pTrack=((AliESDEvent*)event)->GetTrack(lKeyPos);
        AliESDtrack *nTrack=((AliESDEvent*)event)->GetTrack(lKeyNeg);
        
        //Uncertainties: bachelor track as well as V0
        Double_t dy2=t->GetSigmaY2() + pTrack->GetSigmaY2() + nTrack->GetSigmaY2();
        Double_t dz2=t->GetSigmaZ2() + pTrack->GetSigmaZ2() + nTrack->GetSigmaZ2();
        Double_t dx2=dy2;
        
        if( TMath::Abs(fkIfImprovedExtraPrecisionFactor-1.0)>1e-4 ){
            //For testing purposes: override uncertainties, please
            dx2 = fkIfImprovedExtraPrecisionFactor;
            dy2 = fkIfImprovedExtraPrecisionFactor;
            dz2 = fkIfImprovedExtraPrecisionFactor;
        }
        
        //Create dummy V0 track
        //V0 properties to get started
        Double_t xyz[3], pxpypz[3], cv[21];
        for(Int_t ii=0;ii<21;ii++) cv[ii]=0.0; //something small
        
        v->GetXYZ(xyz[0],xyz[1],xyz[2]);
        v->GetPxPyPz( pxpypz[0],pxpypz[1],pxpypz[2] );
        
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        //EXPERIMENTAL: Improve initial position guess based on (neutral!) cowboy/sailor
        //Check bachelor trajectory properties
        Double_t p1[8]; t->GetHelixParameters(p1,b);
        p1[6]=TMath::Sin(p1[2]); p1[7]=TMath::Cos(p1[2]);
        
        if ( fkDoImprovedDCACascDauPropagation ) {
            //Look for XY plane characteristics: determine relevant helix properties
            Double_t lBachRadius = TMath::Abs(1./p1[4]);
            Double_t lBachCenter[2];
            GetHelixCenter( t, lBachCenter, b);
            
            //Algebra: define V0 momentum unit vector
            Double_t ux = pxpypz[0];
            Double_t uy = pxpypz[1];
            Double_t uz = pxpypz[2]; //needed to propagate in 3D (but will be norm to 2D modulus!)
            Double_t umod = TMath::Sqrt(ux*ux+uy*uy);
            ux /= umod; uy /= umod; uz /= umod;
            //perpendicular vector (for projection)
            Double_t vx = -uy;
            Double_t vy = +ux;
            
            //Step 1: calculate distance between line and helix center
            Double_t lDist = (xyz[0]-lBachCenter[0])*vx + (xyz[1]-lBachCenter[1])*vy;
            Double_t lDistSign = lDist / TMath::Abs(lDist);
            
            //Step 2: two cases
            if( TMath::Abs(lDist) > lBachRadius ){
                //only one starting point would sound reasonable
                //check necessary distance to travel for V0
                Double_t lV0travel = (lBachCenter[0]-xyz[0])*ux + (lBachCenter[1]-xyz[1])*uy;
                
                //move V0 forward, please: I already know where to!
                xyz[0] += lV0travel*ux;
                xyz[1] += lV0travel*uy;
                xyz[2] += lV0travel*uz;
                
                //find helix intersection point
                Double_t bX = lBachCenter[0]+lDistSign*lBachRadius*vx;
                Double_t bY = lBachCenter[1]+lDistSign*lBachRadius*vy;
                
                Double_t cs=TMath::Cos(t->GetAlpha());
                Double_t sn=TMath::Sin(t->GetAlpha());
                Double_t lPreprocessX = bX*cs + bY*sn;
                
                //Propagate bachelor track: already know where to!
                t->PropagateTo(lPreprocessX,b);
                
            }else{
                //test two points in which DCAxy=0 for their DCA3D, pick smallest
                //Step 1: find V0-to-center DCA
                Double_t aX = lBachCenter[0]+lDistSign*lDist*vx;
                Double_t aY = lBachCenter[1]+lDistSign*lDist*vy;
                
                //Step 2: find half-axis distance
                Double_t lh = TMath::Sqrt(lBachRadius*lBachRadius - lDist*lDist); //always positive
                
                //Step 3: find 2 points in which XY intersection happens
                Double_t lptAx = aX + lh*ux;
                Double_t lptAy = aY + lh*uy;
                Double_t lptBx = aX - lh*ux;
                Double_t lptBy = aY - lh*uy;
                
                //Step 4: calculate 3D DCA in each point: bachelor
                Double_t xyzptA[3], xyzptB[3];
                Double_t csBach=TMath::Cos(t->GetAlpha());
                Double_t snBach=TMath::Sin(t->GetAlpha());
                Double_t xBachA = lptAx*csBach + lptAy*snBach;
                Double_t xBachB = lptBx*csBach + lptBy*snBach;
                t->GetXYZAt(xBachA,b, xyzptA);
                t->GetXYZAt(xBachB,b, xyzptB);
                
                //Propagate V0 to relevant points
                Double_t lV0travelA = (lptAx-xyz[0])*ux + (lptAy-xyz[1])*uy;
                Double_t lV0travelB = (lptBx-xyz[0])*ux + (lptBy-xyz[1])*uy;
                Double_t lV0xyzptA[3], lV0xyzptB[3];
                lV0xyzptA[0] = xyz[0] + lV0travelA*ux;
                lV0xyzptA[1] = xyz[1] + lV0travelA*uy;
                lV0xyzptA[2] = xyz[2] + lV0travelA*uz;
                lV0xyzptB[0] = xyz[0] + lV0travelB*ux;
                lV0xyzptB[1] = xyz[1] + lV0travelB*uy;
                lV0xyzptB[2] = xyz[2] + lV0travelB*uz;
                
                //Enough info now available to decide on 3D distance
                Double_t l3DdistA = TMath::Sqrt(
                                                TMath::Power(lV0xyzptA[0] - xyzptA[0], 2) +
                                                TMath::Power(lV0xyzptA[1] - xyzptA[1], 2) +
                                                TMath::Power(lV0xyzptA[2] - xyzptA[2], 2)
                                                );
                Double_t l3DdistB = TMath::Sqrt(
                                                TMath::Power(lV0xyzptB[0] - xyzptB[0], 2) +
                                                TMath::Power(lV0xyzptB[1] - xyzptB[1], 2) +
                                                TMath::Power(lV0xyzptB[2] - xyzptB[2], 2)
                                                );
                
                
                if( l3DdistA + 1e-6 < l3DdistB ){
                    //A is the better point! move there, if DCA isn't crazy + x is OK
                    if( l3DdistA < 999 && xBachA > 0.0 && xBachA < 200.0 ) {
                        for(Int_t icoord = 0; icoord<3; icoord++) {
                            xyz[icoord] = lV0xyzptA[icoord];
                        }
                        t->PropagateTo( xBachA , b );
                    }
                }else{
                    //B is the better point! move there, if DCA isn't crazy + x is OK
                    if( l3DdistB < 999 && xBachB > 0.0 && xBachB < 200.0 ) {
                        for(Int_t icoord = 0; icoord<3; icoord++) {
                            xyz[icoord] = lV0xyzptB[icoord];
                        }
                        t->PropagateTo( xBachB , b );
                    }
                }
            }
        }
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        
        //Mockup track for V0 trajectory (no covariance)
        //AliExternalTrackParam *hV0Traj = new AliExternalTrackParam(xyz,pxpypz,cv,+1);
        AliExternalTrackParam lV0TrajObject(xyz,pxpypz,cv,+1), *hV0Traj = &lV0TrajObject;
        hV0Traj->ResetCovariance(1); //won't use
        
        //Re-acquire helix parameters for bachelor (necessary!)
        t->GetHelixParameters(p1,b);
        p1[6]=TMath::Sin(p1[2]);
        p1[7]=TMath::Cos(p1[2]);
        
        Double_t p2[8]; hV0Traj->GetHelixParameters(p2,0.0); //p2[4]=0 -> no curvature (fine, predicted in Evaluate)
        p2[6]=TMath::Sin(p2[2]); p2[7]=TMath::Cos(p2[2]);
        
        Double_t r1[3],g1[3],gg1[3]; Double_t t1=0.;
        Evaluate(p1,t1,r1,g1,gg1);
        Double_t r2[3],g2[3],gg2[3]; Double_t t2=0.;
        Evaluate(p2,t2,r2,g2,gg2);
        
        Double_t dx=r2[0]-r1[0], dy=r2[1]-r1[1], dz=r2[2]-r1[2];
        Double_t dm=dx*dx/dx2 + dy*dy/dy2 + dz*dz/dz2;
        
        Int_t max=27;
        while (max--) {
            Double_t gt1=-(dx*g1[0]/dx2 + dy*g1[1]/dy2 + dz*g1[2]/dz2);
            Double_t gt2=+(dx*g2[0]/dx2 + dy*g2[1]/dy2 + dz*g2[2]/dz2);
            Double_t h11=(g1[0]*g1[0] - dx*gg1[0])/dx2 +
            (g1[1]*g1[1] - dy*gg1[1])/dy2 +
            (g1[2]*g1[2] - dz*gg1[2])/dz2;
            Double_t h22=(g2[0]*g2[0] + dx*gg2[0])/dx2 +
            (g2[1]*g2[1] + dy*gg2[1])/dy2 +
            (g2[2]*g2[2] + dz*gg2[2])/dz2;
            Double_t h12=-(g1[0]*g2[0]/dx2 + g1[1]*g2[1]/dy2 + g1[2]*g2[2]/dz2);
            
            Double_t det=h11*h22-h12*h12;
            
            Double_t dt1,dt2;
            if (TMath::Abs(det)<1.e-33) {
                //(quasi)singular Hessian
                dt1=-gt1; dt2=-gt2;
            } else {
                dt1=-(gt1*h22 - gt2*h12)/det;
                dt2=-(h11*gt2 - h12*gt1)/det;
            }
            
            if ((dt1*gt1+dt2*gt2)>0) {dt1=-dt1; dt2=-dt2;}
            
            //check delta(phase1) ?
            //check delta(phase2) ?
            
            if (TMath::Abs(dt1)/(TMath::Abs(t1)+1.e-3) < 1.e-4)
                if (TMath::Abs(dt2)/(TMath::Abs(t2)+1.e-3) < 1.e-4) {
                    if ((gt1*gt1+gt2*gt2) > 1.e-4/dy2/dy2){
                        AliDebug(1," stopped at not a stationary point !");
                        //Count not stationary point
                        //fHistV0ToBachelorPropagationStatus->Fill(4.5); //bin 5
                    }
                    Double_t lmb=h11+h22; lmb=lmb-TMath::Sqrt(lmb*lmb-4*det);
                    if (lmb < 0.){
                        //Count stopped at not a minimum
                        //fHistV0ToBachelorPropagationStatus->Fill(5.5);
                        AliDebug(1," stopped at not a minimum !");
                    }
                    break;
                }
            
            Double_t dd=dm;
            for (Int_t div=1 ; ; div*=2) {
                Evaluate(p1,t1+dt1,r1,g1,gg1);
                Evaluate(p2,t2+dt2,r2,g2,gg2);
                dx=r2[0]-r1[0]; dy=r2[1]-r1[1]; dz=r2[2]-r1[2];
                dd=dx*dx/dx2 + dy*dy/dy2 + dz*dz/dz2;
                if (dd<dm) break;
                dt1*=0.5; dt2*=0.5;
                if (div>512) {
                    AliDebug(1," overshoot !"); break;
                    //Count overshoots
                    //fHistV0ToBachelorPropagationStatus->Fill(6.5);
                }
            }
            dm=dd;
            
            t1+=dt1;
            t2+=dt2;
            
        }
        
        if (max<=0){
            AliDebug(1," too many iterations !");
            //Count excessive iterations
            //fHistV0ToBachelorPropagationStatus->Fill(7.5);
        }
        
        Double_t cs=TMath::Cos(t->GetAlpha());
        Double_t sn=TMath::Sin(t->GetAlpha());
        Double_t xthis=r1[0]*cs + r1[1]*sn;
        
        //Propagate bachelor to the point of DCA
        if (!t->PropagateTo(xthis,b)) {
            //AliWarning(" propagation failed !";
            //Count curved propagation failures
            //fHistV0ToBachelorPropagationStatus->Fill(8.5);
            return 1e+33;
        }
    
        //V0 distance to bachelor: the desired distance
        Double_t rBachDCAPt[3]; t->GetXYZ(rBachDCAPt);
        dca = v->GetD(rBachDCAPt[0],rBachDCAPt[1],rBachDCAPt[2]);
        //fHistV0ToBachelorPropagationStatus->Fill(9.5);
    }
    
    return dca;
}

//________________________________________________________________________
void AliAnalysisTaskStrEffStudy::Evaluate(const Double_t *h, Double_t t,
                                          Double_t r[3],  //radius vector
                                          Double_t g[3],  //first defivatives
                                          Double_t gg[3]) //second derivatives
{
    //--------------------------------------------------------------------
    // Calculate position of a point on a track and some derivatives
    //--------------------------------------------------------------------
    Double_t phase=h[4]*t+h[2];
    Double_t sn=TMath::Sin(phase), cs=TMath::Cos(phase);
    
    r[0] = h[5];
    r[1] = h[0];
    if (TMath::Abs(h[4])>kAlmost0) {
        r[0] += (sn - h[6])/h[4];
        r[1] -= (cs - h[7])/h[4];
    } else {
        r[0] += t*cs;
        r[1] -= -t*sn;
    }
    r[2] = h[1] + h[3]*t;
    
    g[0] = cs; g[1]=sn; g[2]=h[3];
    
    gg[0]=-h[4]*sn; gg[1]=h[4]*cs; gg[2]=0.;
}

//________________________________________________________________________
Double_t AliAnalysisTaskStrEffStudy::GetErrorInPosition(AliExternalTrackParam *t1) const {
    Double_t alpha=t1->GetAlpha(), /*cs=TMath::Cos(alpha),*/ sn=TMath::Sin(alpha);
    Double_t tmp[3];
    t1->GetPxPyPz(tmp);
    //Double_t px1=tmp[0], py1=tmp[1], pz1=tmp[2];
    t1->GetXYZ(tmp);
    //Double_t  x1=tmp[0],  y1=tmp[1],  z1=tmp[2];
    const Double_t ss=0.0005*0.0005;//a kind of a residual misalignment precision
    Double_t sx1=sn*sn*t1->GetSigmaY2()+ss;// sy1=cs*cs*t1->GetSigmaY2()+ss;
    return sx1;
}


Double_t AliAnalysisTaskStrEffStudy::GetDCAV0Dau( AliExternalTrackParam *pt, AliExternalTrackParam *nt, Double_t &xp, Double_t &xn, Double_t b) {
    //--------------------------------------------------------------
    // Propagates this track and the argument track to the position of the
    // distance of closest approach.
    // Returns the (weighed !) distance of closest approach.
    //--------------------------------------------------------------
    Double_t dy2=nt -> GetSigmaY2() + pt->GetSigmaY2();
    Double_t dz2=nt -> GetSigmaZ2() + pt->GetSigmaZ2();
    Double_t dx2=dy2;
    
    Double_t p1[8]; nt->GetHelixParameters(p1,b);
    p1[6]=TMath::Sin(p1[2]); p1[7]=TMath::Cos(p1[2]);
    Double_t p2[8]; pt->GetHelixParameters(p2,b);
    p2[6]=TMath::Sin(p2[2]); p2[7]=TMath::Cos(p2[2]);
    
    if( fkDoImprovedDCAV0DauPropagation){
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        // V0 preprocessing: analytical estimate of DCAxy position
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        Double_t nhelix[6], phelix[6];
        nt->GetHelixParameters(nhelix,b);
        pt->GetHelixParameters(phelix,b);
        Double_t lNegCenterR[2], lPosCenterR[2];
        
        //Negative track parameters in XY
        GetHelixCenter( nt , lNegCenterR, b);
        Double_t xNegCenter = lNegCenterR[0];
        Double_t yNegCenter = lNegCenterR[1];
        Double_t NegRadius = TMath::Abs(1./nhelix[4]);
        
        //Positive track parameters in XY
        GetHelixCenter( pt , lPosCenterR, b );
        Double_t xPosCenter = lPosCenterR[0];
        Double_t yPosCenter = lPosCenterR[1];
        Double_t PosRadius = TMath::Abs(1./phelix[4]);
        
        //Define convenient coordinate system
        //Logical zero: position of negative center
        Double_t ux = xPosCenter - xNegCenter;
        Double_t uy = yPosCenter - yNegCenter;
        
        //Check center-to-center distance
        Double_t lDist = TMath::Sqrt(
                                     TMath::Power( xNegCenter - xPosCenter , 2) +
                                     TMath::Power( yNegCenter - yPosCenter , 2)
                                     );
        //Normalize ux, uz to unit vector
        ux /= lDist; uy /= lDist;
        
        //Calculate perpendicular vector (normalized)
        Double_t vx = -uy;
        Double_t vy = +ux;
        
        Double_t lPreprocessDCAxy = 1e+3; //define outside scope
        Double_t lPreprocessxp = pt->GetX(); //start at current location
        Double_t lPreprocessxn = nt->GetX(); //start at current location
        
        if( lDist > NegRadius + PosRadius ){
            //================================================================
            //Case 1: distance bigger than sum of radii ("gamma-like")
            //        re-position tracks along the center-to-center axis
            //Re-position negative track
            Double_t xNegOptPosition = xNegCenter + NegRadius*ux;
            Double_t yNegOptPosition = yNegCenter + NegRadius*uy;
            Double_t csNeg=TMath::Cos(nt->GetAlpha());
            Double_t snNeg=TMath::Sin(nt->GetAlpha());
            Double_t xThisNeg=xNegOptPosition*csNeg + yNegOptPosition*snNeg;
            
            //Re-position positive track
            Double_t xPosOptPosition = xPosCenter - PosRadius*ux;
            Double_t yPosOptPosition = yPosCenter - PosRadius*uy;
            Double_t csPos=TMath::Cos(pt->GetAlpha());
            Double_t snPos=TMath::Sin(pt->GetAlpha());
            Double_t xThisPos=xPosOptPosition*csPos + yPosOptPosition*snPos;
            
            if( xThisNeg < fV0VertexerSels[6] && xThisPos < fV0VertexerSels[6] && xThisNeg > 0.0 && xThisPos > 0.0){
                Double_t lCase1NegR[3]; nt->GetXYZAt(xThisNeg,b, lCase1NegR);
                Double_t lCase1PosR[3]; pt->GetXYZAt(xThisPos,b, lCase1PosR);
                lPreprocessDCAxy = TMath::Sqrt(
                                               TMath::Power(lCase1NegR[0]-lCase1PosR[0],2)+
                                               TMath::Power(lCase1NegR[1]-lCase1PosR[1],2)+
                                               TMath::Power(lCase1NegR[2]-lCase1PosR[2],2)
                                               );
                //Pass coordinates
                if( lPreprocessDCAxy<999){
                    lPreprocessxp = xThisPos;
                    lPreprocessxn = xThisNeg;
                }
            }
            //================================================================
        } else {
            if( lDist > TMath::Abs(NegRadius-PosRadius) ){ //otherwise this algorithm will fail!
                //================================================================
                //Case 2: distance smaller than sum of radii (cowboy/sailor configs)
                
                //Calculate coordinate for radical line
                Double_t lRadical = (lDist*lDist - PosRadius*PosRadius + NegRadius*NegRadius) / (2*lDist);
                
                //Calculate absolute displacement from center-to-center axis
                Double_t lDisplace = (0.5/lDist) * TMath::Sqrt(
                                                               (-lDist + PosRadius - NegRadius) *
                                                               (-lDist - PosRadius + NegRadius) *
                                                               (-lDist + PosRadius + NegRadius) *
                                                               ( lDist + PosRadius + NegRadius)
                                                               );
                
                Double_t lCase2aDCA = 1e+3;
                Double_t lCase2bDCA = 1e+3;
                
                //2 cases: positive and negative displacement
                Double_t xNegOptPosition[2], yNegOptPosition[2], xPosOptPosition[2], yPosOptPosition[2];
                Double_t csNeg, snNeg, csPos, snPos;
                Double_t xThisNeg[2], xThisPos[2];
                
                csNeg=TMath::Cos(nt->GetAlpha());
                snNeg=TMath::Sin(nt->GetAlpha());
                csPos=TMath::Cos(pt->GetAlpha());
                snPos=TMath::Sin(pt->GetAlpha());
                
                //Case 2a: Positive displacement along v vector
                //Re-position negative track
                xNegOptPosition[0] = xNegCenter + lRadical*ux + lDisplace*vx;
                yNegOptPosition[0] = yNegCenter + lRadical*uy + lDisplace*vy;
                xThisNeg[0] = xNegOptPosition[0]*csNeg + yNegOptPosition[0]*snNeg;
                //Re-position positive track
                xPosOptPosition[0] = xNegCenter + lRadical*ux + lDisplace*vx;
                yPosOptPosition[0] = yNegCenter + lRadical*uy + lDisplace*vy;
                xThisPos[0] = xPosOptPosition[0]*csPos + yPosOptPosition[0]*snPos;
                
                //Case 2b: Negative displacement along v vector
                //Re-position negative track
                xNegOptPosition[1] = xNegCenter + lRadical*ux - lDisplace*vx;
                yNegOptPosition[1] = yNegCenter + lRadical*uy - lDisplace*vy;
                xThisNeg[1] = xNegOptPosition[1]*csNeg + yNegOptPosition[1]*snNeg;
                //Re-position positive track
                xPosOptPosition[1] = xNegCenter + lRadical*ux - lDisplace*vx;
                yPosOptPosition[1] = yNegCenter + lRadical*uy - lDisplace*vy;
                xThisPos[1] = xPosOptPosition[1]*csPos + yPosOptPosition[1]*snPos;
                
                //Test the two cases, please
                
                //Case 2a
                if( xThisNeg[0] < 200 && xThisPos[0] < 200 && xThisNeg[0] > 0.0 && xThisPos[0] > 0.0 ){
                    Double_t lCase2aNegR[3]; nt->GetXYZAt(xThisNeg[0],b, lCase2aNegR);
                    Double_t lCase2aPosR[3]; pt->GetXYZAt(xThisPos[0],b, lCase2aPosR);
                    lCase2aDCA = TMath::Sqrt(
                                             TMath::Power(lCase2aNegR[0]-lCase2aPosR[0],2)+
                                             TMath::Power(lCase2aNegR[1]-lCase2aPosR[1],2)+
                                             TMath::Power(lCase2aNegR[2]-lCase2aPosR[2],2)
                                             );
                }
                
                //Case 2b
                if( xThisNeg[1] < 200 && xThisPos[1] < 200 && xThisNeg[1] > 0.0 && xThisPos[1] > 0.0 ){
                    Double_t lCase2bNegR[3]; nt->GetXYZAt(xThisNeg[1],b, lCase2bNegR);
                    Double_t lCase2bPosR[3]; pt->GetXYZAt(xThisPos[1],b, lCase2bPosR);
                    lCase2bDCA = TMath::Sqrt(
                                             TMath::Power(lCase2bNegR[0]-lCase2bPosR[0],2)+
                                             TMath::Power(lCase2bNegR[1]-lCase2bPosR[1],2)+
                                             TMath::Power(lCase2bNegR[2]-lCase2bPosR[2],2)
                                             );
                }
                
                //Minor detail: all things being equal, prefer closest X
                Double_t lCase2aSumX = xThisPos[0]+xThisNeg[0];
                Double_t lCase2bSumX = xThisPos[1]+xThisNeg[1];
                
                Double_t lDCAxySmallestR = lCase2aDCA;
                Double_t lxpSmallestR = xThisPos[0];
                Double_t lxnSmallestR = xThisNeg[0];
                
                Double_t lDCAxyLargestR = lCase2bDCA;
                Double_t lxpLargestR = xThisPos[1];
                Double_t lxnLargestR = xThisNeg[1];
                
                if( lCase2bSumX+1e-6 < lCase2aSumX ){
                    lDCAxySmallestR = lCase2bDCA;
                    lxpSmallestR = xThisPos[1];
                    lxnSmallestR = xThisNeg[1];
                    lDCAxyLargestR = lCase2aDCA;
                    lxpLargestR = xThisPos[0];
                    lxnLargestR = xThisNeg[0];
                }
                
                //Pass conclusion to lPreprocess variables, please
                lPreprocessDCAxy = lDCAxySmallestR;
                lPreprocessxp = lxpSmallestR;
                lPreprocessxn = lxnSmallestR;
                if( lDCAxyLargestR+1e-6 < lDCAxySmallestR ){ //beware epsilon: numerical calculations are unstable here
                    lPreprocessDCAxy = lDCAxyLargestR;
                    lPreprocessxp = lxpLargestR;
                    lPreprocessxn = lxnLargestR;
                }
                //Protection against something too crazy, please
                if( lPreprocessDCAxy>999){
                    lPreprocessxp = pt->GetX(); //start at current location
                    lPreprocessxn = nt->GetX(); //start at current location
                }
            }
        }
        //End of preprocessing stage!
        //at this point lPreprocessxp, lPreprocessxn are already good starting points: update helixparams
        if( lPreprocessDCAxy < 999 ) { //some improvement... otherwise discard in all cases, please
            nt->PropagateTo(lPreprocessxn, b);
            pt->PropagateTo(lPreprocessxp, b);
        }
        
        //don't redefine!
        nt->GetHelixParameters(p1,b);
        p1[6]=TMath::Sin(p1[2]);
        p1[7]=TMath::Cos(p1[2]);
        pt->GetHelixParameters(p2,b);
        p2[6]=TMath::Sin(p2[2]);
        p2[7]=TMath::Cos(p2[2]);
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    }
    
    Double_t r1[3],g1[3],gg1[3]; Double_t t1=0.;
    Evaluate(p1,t1,r1,g1,gg1);
    Double_t r2[3],g2[3],gg2[3]; Double_t t2=0.;
    Evaluate(p2,t2,r2,g2,gg2);
    
    Double_t dx=r2[0]-r1[0], dy=r2[1]-r1[1], dz=r2[2]-r1[2];
    Double_t dm=dx*dx/dx2 + dy*dy/dy2 + dz*dz/dz2;
    
    Int_t max=27;
    while (max--) {
        Double_t gt1=-(dx*g1[0]/dx2 + dy*g1[1]/dy2 + dz*g1[2]/dz2);
        Double_t gt2=+(dx*g2[0]/dx2 + dy*g2[1]/dy2 + dz*g2[2]/dz2);
        Double_t h11=(g1[0]*g1[0] - dx*gg1[0])/dx2 +
        (g1[1]*g1[1] - dy*gg1[1])/dy2 +
        (g1[2]*g1[2] - dz*gg1[2])/dz2;
        Double_t h22=(g2[0]*g2[0] + dx*gg2[0])/dx2 +
        (g2[1]*g2[1] + dy*gg2[1])/dy2 +
        (g2[2]*g2[2] + dz*gg2[2])/dz2;
        Double_t h12=-(g1[0]*g2[0]/dx2 + g1[1]*g2[1]/dy2 + g1[2]*g2[2]/dz2);
        
        Double_t det=h11*h22-h12*h12;
        
        Double_t dt1,dt2;
        if (TMath::Abs(det)<1.e-33) {
            //(quasi)singular Hessian
            dt1=-gt1; dt2=-gt2;
        } else {
            dt1=-(gt1*h22 - gt2*h12)/det;
            dt2=-(h11*gt2 - h12*gt1)/det;
        }
        
        if ((dt1*gt1+dt2*gt2)>0) {dt1=-dt1; dt2=-dt2;}
        
        //check delta(phase1) ?
        //check delta(phase2) ?
        
        if (TMath::Abs(dt1)/(TMath::Abs(t1)+1.e-3) < 1.e-4)
            if (TMath::Abs(dt2)/(TMath::Abs(t2)+1.e-3) < 1.e-4) {
                if ((gt1*gt1+gt2*gt2) > 1.e-4/dy2/dy2)
                    AliDebug(1," stopped at not a stationary point !");
                Double_t lmb=h11+h22; lmb=lmb-TMath::Sqrt(lmb*lmb-4*det);
                if (lmb < 0.)
                    AliDebug(1," stopped at not a minimum !");
                break;
            }
        
        Double_t dd=dm;
        for (Int_t div=1 ; ; div*=2) {
            Evaluate(p1,t1+dt1,r1,g1,gg1);
            Evaluate(p2,t2+dt2,r2,g2,gg2);
            dx=r2[0]-r1[0]; dy=r2[1]-r1[1]; dz=r2[2]-r1[2];
            dd=dx*dx/dx2 + dy*dy/dy2 + dz*dz/dz2;
            if (dd<dm) break;
            dt1*=0.5; dt2*=0.5;
            if (div>512) {
                AliDebug(1," overshoot !"); break;
            }
        }
        dm=dd;
        
        t1+=dt1;
        t2+=dt2;

    }

    if (max<=0) AliDebug(1," too many iterations !");

    Double_t cs=TMath::Cos(nt->GetAlpha());
    Double_t sn=TMath::Sin(nt->GetAlpha());
    xn=r1[0]*cs + r1[1]*sn;

    cs=TMath::Cos(pt->GetAlpha());
    sn=TMath::Sin(pt->GetAlpha());
    xp=r2[0]*cs + r2[1]*sn;

    return TMath::Sqrt(dm*TMath::Sqrt(dy2*dz2));
}

///________________________________________________________________________
void AliAnalysisTaskStrEffStudy::GetHelixCenter(const AliExternalTrackParam *track,Double_t center[2], Double_t b){
    // Copied from AliV0ReaderV1::GetHelixCenter
    // Get Center of the helix track parametrization

    Int_t charge=track->Charge();

    Double_t	helix[6];
    track->GetHelixParameters(helix,b);

    Double_t xpos =	helix[5];
    Double_t ypos =	helix[0];
    Double_t radius = TMath::Abs(1./helix[4]);
    Double_t phi = helix[2];
    if(phi < 0){
        phi = phi + 2*TMath::Pi();
    }
    phi -= TMath::Pi()/2.;
    Double_t xpoint =	radius * TMath::Cos(phi);
    Double_t ypoint =	radius * TMath::Sin(phi);
    if(b<0){
        if(charge > 0){
            xpoint = - xpoint;
            ypoint = - ypoint;
        }
        /*if(charge < 0){
            xpoint =	xpoint;
            ypoint =	ypoint;
        }*/
    }
    if(b>0){
        /*if(charge > 0){
            xpoint =	xpoint;
            ypoint =	ypoint;
        }*/
        if(charge < 0){
            xpoint = - xpoint;
            ypoint = - ypoint;
        }
    }
    center[0] =	xpos + xpoint;
    center[1] =	ypos + ypoint;
    return;
}
