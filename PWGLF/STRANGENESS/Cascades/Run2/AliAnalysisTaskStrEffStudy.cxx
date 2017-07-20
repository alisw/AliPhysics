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

ClassImp(AliAnalysisTaskStrEffStudy)

AliAnalysisTaskStrEffStudy::AliAnalysisTaskStrEffStudy()
: AliAnalysisTaskSE(), fListHist(0), fListV0(0), fListCascade(0), fTreeEvent(0), fTreeV0(0), fTreeCascade(0), fPIDResponse(0), fESDtrackCuts(0), fESDtrackCutsITSsa2010(0), fESDtrackCutsGlobal2015(0), fUtils(0), fRand(0),

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
fkIfImprovedPerformInitialLinearPropag( kFALSE ),
fkIfImprovedExtraPrecisionFactor ( 1.0 ),
fkDebugWrongPIDForTracking ( kFALSE ),
fkDebugBump( kFALSE ),
fkDebugOOBPileup(kFALSE),
fkDoExtraEvSels( kTRUE ),
fPrecisionCutoffCascadeDCA(1e-4),
fMaxIterationsCascadeDCA(27),

//---> Flags controlling Cascade TTree output
fkSaveCascadeTree       ( kTRUE  ),
fkDownScaleCascade      ( kTRUE  ),
fDownScaleFactorCascade ( 0.001  ),
fMinPtToSave( 0.00   ) ,
fMaxPtToSave( 100.00 ) ,

//---> Flags controlling Vertexers
fkRunVertexers    ( kFALSE ),
fkUseLightVertexer ( kTRUE ),
fkDoV0Refit ( kTRUE ),
fkExtraCleanup    ( kTRUE ),
fkSaveGoodTracks( kTRUE ),

//---> Flag controlling trigger selection
fTrigType(AliVEvent::kMB),

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
fTreeVariableInvMassK0s(0),
fTreeVariableInvMassLambda(0),
fTreeVariableInvMassAntiLambda(0),
fTreeVariablePID(0),
fTreeVariablePIDPositive(0),
fTreeVariablePIDNegative(0),
fTreeVariablePtMC(0),
fTreeVariableRapMC(0),
fTreeVariablePosAlpha(0),
fTreeVariablePosSigmaY2(0),
fTreeVariablePosSigmaZ2(0),
fTreeVariableNegAlpha(0),
fTreeVariableNegSigmaY2(0),
fTreeVariableNegSigmaZ2(0),

//---> Variables for fTreeCascade
fTreeCascVarCentrality(0),
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

fTreeCascVarDCACascDaughtersClassical(0),
fTreeCascVarCascPropagationClassical(0),

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

fTreeCascVarCascPropagationImprovedIterations(0),
fTreeCascVarCascPropagationImprovedStatus(0),
fTreeCascVarDCACascDaughtersImproved(0),

fTreeCascVarImprovedDecayX(0),
fTreeCascVarImprovedDecayY(0),
fTreeCascVarImprovedDecayZ(0),
fTreeCascVarImprovedCascCosPointingAngle(0),
fTreeCascVarImprovedCascDCAxyToPV(0),
fTreeCascVarImprovedCascDCAzToPV(0),

fTreeCascVarImprovedInvMassXiMinus(0),
fTreeCascVarImprovedInvMassXiPlus(0),
fTreeCascVarImprovedInvMassOmegaMinus(0),
fTreeCascVarImprovedInvMassOmegaPlus(0),

fTreeCascVarPIDPositive(0),
fTreeCascVarPIDNegative(0),
fTreeCascVarPIDBachelor(0),
fTreeCascVarPID(0),
fTreeCascVarPtMC(0),
fTreeCascVarRapMC(0),

fTreeCascVarPosDistanceToTrueDecayPt(0),
fTreeCascVarNegDistanceToTrueDecayPt(0),
fTreeCascVarBachDistanceToTrueDecayPt(0),
fTreeCascVarV0DistanceToTrueDecayPt(0),

fTreeCascVarBachPropagationParameterClassical(0),
fTreeCascVarBachPropagationParameterImproved(0),

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
fHistGeneratedPtVsYVsCentralityOmegaPlus(0)
//------------------------------------------------
// Tree Variables
{
    
}

AliAnalysisTaskStrEffStudy::AliAnalysisTaskStrEffStudy(Bool_t lSaveEventTree, Bool_t lSaveV0Tree, Bool_t lSaveCascadeTree, const char *name, TString lExtraOptions)
: AliAnalysisTaskSE(name), fListHist(0), fListV0(0), fListCascade(0), fTreeEvent(0), fTreeV0(0), fTreeCascade(0), fPIDResponse(0), fESDtrackCuts(0), fESDtrackCutsITSsa2010(0), fESDtrackCutsGlobal2015(0), fUtils(0), fRand(0),

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
fkIfImprovedPerformInitialLinearPropag( kFALSE ),
fkIfImprovedExtraPrecisionFactor ( 1.0 ),
fkDebugWrongPIDForTracking ( kFALSE ), //also for cascades...
fkDebugBump( kFALSE ),
fkDebugOOBPileup(kFALSE),
fkDoExtraEvSels( kTRUE ),
fPrecisionCutoffCascadeDCA(1e-4),
fMaxIterationsCascadeDCA(27),

//---> Flags controlling Cascade TTree output
fkSaveCascadeTree       ( kTRUE  ),
fkDownScaleCascade      ( kTRUE  ),
fDownScaleFactorCascade ( 0.001  ),
fMinPtToSave( 0.00   ) ,
fMaxPtToSave( 100.00 ) ,

//---> Flags controlling Vertexers
fkRunVertexers    ( kFALSE ),
fkUseLightVertexer ( kTRUE ),
fkDoV0Refit ( kTRUE ),
fkExtraCleanup    ( kTRUE ),
fkSaveGoodTracks( kTRUE ),

//---> Flag controlling trigger selection
fTrigType(AliVEvent::kMB),

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
fTreeVariableInvMassK0s(0),
fTreeVariableInvMassLambda(0),
fTreeVariableInvMassAntiLambda(0),
fTreeVariablePID(0),
fTreeVariablePIDPositive(0),
fTreeVariablePIDNegative(0),
fTreeVariablePtMC(0),
fTreeVariableRapMC(0),
fTreeVariablePosAlpha(0),
fTreeVariablePosSigmaY2(0),
fTreeVariablePosSigmaZ2(0),
fTreeVariableNegAlpha(0),
fTreeVariableNegSigmaY2(0),
fTreeVariableNegSigmaZ2(0),

//---> Variables for fTreeCascade
fTreeCascVarCentrality(0),
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
fTreeCascVarV0CosineOfPointingAngle(0),
fTreeCascVarDCAV0ToPrimVtx(0),
fTreeCascVarDCAxyV0ToPrimVtx(0),
fTreeCascVarInvMassLambda(0),
fTreeCascVarInvMassAntiLambda(0),

fTreeCascVarDCACascDaughtersClassical(0),
fTreeCascVarCascPropagationClassical(0),

fTreeCascVarDecayX(0),
fTreeCascVarDecayY(0),
fTreeCascVarDecayZ(0),
fTreeCascVarCascCosPointingAngle(0),

fTreeCascVarInvMassXiMinus(0),
fTreeCascVarInvMassXiPlus(0),
fTreeCascVarInvMassOmegaMinus(0),
fTreeCascVarInvMassOmegaPlus(0),

fTreeCascVarCascPropagationImprovedIterations(0),
fTreeCascVarCascPropagationImprovedStatus(0),
fTreeCascVarDCACascDaughtersImproved(0),

fTreeCascVarImprovedDecayX(0),
fTreeCascVarImprovedDecayY(0),
fTreeCascVarImprovedDecayZ(0),
fTreeCascVarImprovedCascCosPointingAngle(0),
fTreeCascVarImprovedCascDCAxyToPV(0),
fTreeCascVarImprovedCascDCAzToPV(0),

fTreeCascVarImprovedInvMassXiMinus(0),
fTreeCascVarImprovedInvMassXiPlus(0),
fTreeCascVarImprovedInvMassOmegaMinus(0),
fTreeCascVarImprovedInvMassOmegaPlus(0),

fTreeCascVarPIDPositive(0),
fTreeCascVarPIDNegative(0),
fTreeCascVarPIDBachelor(0),
fTreeCascVarPID(0),
fTreeCascVarPtMC(0),
fTreeCascVarRapMC(0),

fTreeCascVarPosDistanceToTrueDecayPt(0),
fTreeCascVarNegDistanceToTrueDecayPt(0),
fTreeCascVarBachDistanceToTrueDecayPt(0),
fTreeCascVarV0DistanceToTrueDecayPt(0),

fTreeCascVarBachPropagationParameterClassical(0),
fTreeCascVarBachPropagationParameterImproved(0),

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
fHistGeneratedPtVsYVsCentralityOmegaPlus(0)
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
    DefineOutput(5, TTree::Class()); // V0 Tree output
    DefineOutput(6, TTree::Class()); // Cascade Tree output
    
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
    fTreeV0->Branch("fTreeVariableInvMassK0s",       &fTreeVariableInvMassK0s,       "fTreeVariableInvMassK0s/F");
    fTreeV0->Branch("fTreeVariableInvMassLambda",    &fTreeVariableInvMassLambda,    "fTreeVariableInvMassLambda/F");
    fTreeV0->Branch("fTreeVariableInvMassAntiLambda",&fTreeVariableInvMassAntiLambda,"fTreeVariableInvMassAntiLambda/F");
    fTreeV0->Branch("fTreeVariablePID",&fTreeVariablePID,"fTreeVariablePID/I");
    fTreeV0->Branch("fTreeVariablePIDPositive",&fTreeVariablePIDPositive,"fTreeVariablePIDPositive/I");
    fTreeV0->Branch("fTreeVariablePIDNegative",&fTreeVariablePIDNegative,"fTreeVariablePIDNegative/I");
    fTreeV0->Branch("fTreeVariablePtMC",&fTreeVariablePtMC,"fTreeVariablePtMC/F");
    fTreeV0->Branch("fTreeVariableRapMC",&fTreeVariableRapMC,"fTreeVariableRapMC/F");
    
    //Uncertainties
    fTreeV0->Branch("fTreeVariablePosAlpha",&fTreeVariablePosAlpha,"fTreeVariablePosAlpha/F");
    fTreeV0->Branch("fTreeVariablePosSigmaY2",&fTreeVariablePosSigmaY2,"fTreeVariablePosSigmaY2/F");
    fTreeV0->Branch("fTreeVariablePosSigmaZ2",&fTreeVariablePosSigmaZ2,"fTreeVariablePosSigmaZ2/F");
    fTreeV0->Branch("fTreeVariableNegAlpha",&fTreeVariableNegAlpha,"fTreeVariableNegAlpha/F");
    fTreeV0->Branch("fTreeVariableNegSigmaY2",&fTreeVariableNegSigmaY2,"fTreeVariableNegSigmaY2/F");
    fTreeV0->Branch("fTreeVariableNegSigmaZ2",&fTreeVariableNegSigmaZ2,"fTreeVariableNegSigmaZ2/F");
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
    fTreeCascade->Branch("fTreeCascVarDCACascDaughtersClassical",&fTreeCascVarDCACascDaughtersClassical,"fTreeCascVarDCACascDaughtersClassical/F");
    fTreeCascade->Branch("fTreeCascVarCascPropagationClassical",&fTreeCascVarCascPropagationClassical,"fTreeCascVarCascPropagationClassical/O");
    
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
    
    //IMPROVED BACK-PROPAGATION VARIABLES
    fTreeCascade->Branch("fTreeCascVarCascPropagationImprovedIterations",&fTreeCascVarCascPropagationImprovedIterations,"fTreeCascVarCascPropagationImprovedIterations/I");
    fTreeCascade->Branch("fTreeCascVarCascPropagationImprovedStatus",&fTreeCascVarCascPropagationImprovedStatus,"fTreeCascVarCascPropagationImprovedStatus/I");
    fTreeCascade->Branch("fTreeCascVarDCACascDaughtersImproved",&fTreeCascVarDCACascDaughtersImproved,"fTreeCascVarDCACascDaughtersImproved/F");
    
    fTreeCascade->Branch("fTreeCascVarImprovedDecayX",&fTreeCascVarImprovedDecayX,"fTreeCascVarImprovedDecayX/F");
    fTreeCascade->Branch("fTreeCascVarImprovedDecayY",&fTreeCascVarImprovedDecayY,"fTreeCascVarImprovedDecayY/F");
    fTreeCascade->Branch("fTreeCascVarImprovedDecayZ",&fTreeCascVarImprovedDecayZ,"fTreeCascVarImprovedDecayZ/F");
    fTreeCascade->Branch("fTreeCascVarImprovedCascCosPointingAngle",&fTreeCascVarImprovedCascCosPointingAngle,"fTreeCascVarImprovedCascCosPointingAngle/F");
    fTreeCascade->Branch("fTreeCascVarImprovedCascDCAxyToPV",&fTreeCascVarImprovedCascDCAxyToPV,"fTreeCascVarImprovedCascDCAxyToPV/F");
    fTreeCascade->Branch("fTreeCascVarImprovedCascDCAzToPV",&fTreeCascVarImprovedCascDCAzToPV,"fTreeCascVarImprovedCascDCAzToPV/F");
    
    fTreeCascade->Branch("fTreeCascVarImprovedInvMassXiMinus",&fTreeCascVarImprovedInvMassXiMinus,"fTreeCascVarImprovedInvMassXiMinus/F");
    fTreeCascade->Branch("fTreeCascVarImprovedInvMassXiPlus",&fTreeCascVarImprovedInvMassXiPlus,"fTreeCascVarImprovedInvMassXiPlus/F");
    fTreeCascade->Branch("fTreeCascVarImprovedInvMassOmegaMinus",&fTreeCascVarImprovedInvMassOmegaMinus,"fTreeCascVarImprovedInvMassOmegaMinus/F");
    fTreeCascade->Branch("fTreeCascVarImprovedInvMassOmegaPlus",&fTreeCascVarImprovedInvMassOmegaPlus,"fTreeCascVarImprovedInvMassOmegaPlus/F");
    
    
    //MC VARIABLES
    fTreeCascade->Branch("fTreeCascVarPIDPositive",&fTreeCascVarPIDPositive,"fTreeCascVarPIDPositive/I");
    fTreeCascade->Branch("fTreeCascVarPIDNegative",&fTreeCascVarPIDNegative,"fTreeCascVarPIDNegative/I");
    fTreeCascade->Branch("fTreeCascVarPIDBachelor",&fTreeCascVarPIDBachelor,"fTreeCascVarPIDBachelor/I");
    fTreeCascade->Branch("fTreeCascVarPID",&fTreeCascVarPID,"fTreeCascVarPID/I");
    fTreeCascade->Branch("fTreeCascVarPtMC",&fTreeCascVarPtMC,"fTreeCascVarPtMC/F");
    fTreeCascade->Branch("fTreeCascVarRapMC",&fTreeCascVarRapMC,"fTreeCascVarRapMC/F");
    
    fTreeCascade->Branch("fTreeCascVarPosDistanceToTrueDecayPt",&fTreeCascVarPosDistanceToTrueDecayPt,"fTreeCascVarPosDistanceToTrueDecayPt/F");
    fTreeCascade->Branch("fTreeCascVarNegDistanceToTrueDecayPt",&fTreeCascVarNegDistanceToTrueDecayPt,"fTreeCascVarNegDistanceToTrueDecayPt/F");
    fTreeCascade->Branch("fTreeCascVarBachDistanceToTrueDecayPt",&fTreeCascVarBachDistanceToTrueDecayPt,"fTreeCascVarBachDistanceToTrueDecayPt/F");
    fTreeCascade->Branch("fTreeCascVarV0DistanceToTrueDecayPt",&fTreeCascVarV0DistanceToTrueDecayPt,"fTreeCascVarV0DistanceToTrueDecayPt/F");
    
    fTreeCascade->Branch("fTreeCascVarBachPropagationParameterClassical",&fTreeCascVarBachPropagationParameterClassical,"fTreeCascVarBachPropagationParameterClassical/F");
    fTreeCascade->Branch("fTreeCascVarBachPropagationParameterImproved",&fTreeCascVarBachPropagationParameterImproved,"fTreeCascVarBachPropagationParameterImproved/F");
    
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
    
}// end UserCreateOutputObjects


//________________________________________________________________________
void AliAnalysisTaskStrEffStudy::UserExec(Option_t *)
{
    // Main loop
    // Called for each event
    AliESDEvent *lESDevent = 0x0;
    AliMCEvent  *lMCevent  = 0x0;
    AliStack    *lMCstack  = 0x0;
    
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
    
    Float_t lPercentile = 500;
    Float_t lPercentileEmbeddedSelection = 500;
    Int_t lEvSelCode = 100;
    AliMultSelection *MultSelection = (AliMultSelection*) lESDevent -> FindListObject("MultSelection");
    if( !MultSelection) {
        //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
        AliWarning("AliMultSelection object not found!");
    } else {
        //V0M Multiplicity Percentile
        lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
        lPercentileEmbeddedSelection = MultSelection->GetMultiplicityPercentile("V0M", kTRUE );
        //Event Selection Code
        lEvSelCode = MultSelection->GetEvSelCode();
    }
    
    //just ask AliMultSelection. It will know.
    fMVPileupFlag = kFALSE;
    fMVPileupFlag = MultSelection->GetThisEventIsNotPileupMV();
    
    fCentrality = lPercentile;
    
    //Let's find out why efficiency is so centrality dependent, please 
    fTreeCascVarCentrality = lPercentile;
    
    if( lEvSelCode != 0 ) {
        PostData(1, fListHist    );
        PostData(2, fListV0      );
        PostData(3, fListCascade );
        PostData(4, fTreeEvent   );
        PostData(5, fTreeV0      );
        PostData(6, fTreeCascade );
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
        if ( (TMath::Abs(lThisPDG) == 3312) || (TMath::Abs(lThisPDG) == 3334) || (TMath::Abs(lThisPDG) == 3122) || lThisPDG == 310 )
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
        }
    }//End of loop on tracks
    //----- End Loop on Cascades ------------------------------------------------------------
    
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
    
    //-------------------------------------------------
    // V0s from scratch: locate findable V0 candidates
    //-------------------------------------------------
    
    //Particles of interest
    Int_t lV0Types[3]          = { 310, 3122, -3122};
    Int_t lV0TypesPDau[3]      = { 211, 2212,   211};
    Int_t lV0TypesNDau[3]      = {-211, -211, -2212};
    
    //Number of tracks
    Long_t lNTracks = lESDevent->GetNumberOfTracks();
    Double_t b      = lESDevent->GetMagneticField();
    
    //pos/neg daughters
    TArrayI lTrackArray      (lNTracks);
    TArrayI lTrackMotherArray(lNTracks);
    
    Long_t nTracksOfInterest = 0;
    
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
        if ( lParticleMother->GetNDaughters()!=2 ) continue;
        
        Bool_t lOfDesiredType = kFALSE;
        for(Int_t iType=0; iType<3; iType++){
            if( lParticleMotherPDG == lV0Types[iType] ) lOfDesiredType = kTRUE;
        }
        if( !lOfDesiredType ) continue;
        
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
                AliESDtrack *esdTrack1 = lESDevent->GetTrack( lTrackArray[iTrack] );
                AliESDtrack *esdTrack2 = lESDevent->GetTrack( lTrackArray[jTrack] );
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
        Double_t xn, xp, dca=esdTrackNeg->GetDCA(esdTrackPos,lMagneticField,xn,xp);
        
        //Correct for beam pipe material
        //Warning: this is an unfinished implementation and should not do much at this stage
        AliExternalTrackParam nt(*esdTrackNeg), pt(*esdTrackPos);
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
            dca=nt.GetDCA(&pt,lMagneticField,xn,xp);
        }
        
        fTreeVariableDcaV0Daughters = dca; //Pass to TTree object, please
        
        //Actual propagation
        fTreeVariableNegPropagStatus = nt.PropagateTo(xn,lMagneticField);
        fTreeVariablePosPropagStatus = pt.PropagateTo(xp,lMagneticField);
        
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
    
    //____________________________________________________________________________
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
        AliESDtrack *esdTrackNeg = lESDevent->GetTrack( lNegTrackArray[iV0] );
        
        //3c: Get perfect MC information for bookkeeping
        Int_t lblPosV0Dghter = (Int_t) TMath::Abs( esdTrackPos->GetLabel() );
        Int_t lblNegV0Dghter = (Int_t) TMath::Abs( esdTrackNeg->GetLabel() );
        
        TParticle* mcPosV0Dghter = lMCstack->Particle( lblPosV0Dghter );
        TParticle* mcNegV0Dghter = lMCstack->Particle( lblNegV0Dghter );
        
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
        
        Double_t xn, xp, dca=esdTrackNeg->GetDCA(esdTrackPos,lMagneticField,xn,xp);
        
        //Correct for beam pipe material
        //Warning: this is an unfinished implementation and should not do much at this stage
        AliExternalTrackParam nt(*esdTrackNeg), pt(*esdTrackPos);
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
            dca=nt.GetDCA(&pt,lMagneticField,xn,xp);
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
        
        //Final step: get estimated masses under different mass hypotheses
        vertex.ChangeMassHypothesis(3122);
        fTreeCascVarInvMassLambda = vertex.GetEffMass();
        vertex.ChangeMassHypothesis(-3122);
        fTreeCascVarInvMassAntiLambda = vertex.GetEffMass();
        
        //---] Cascade PART [----------------------------------------------
        
        //Step 1: propagation (CLASSICAL)
        AliESDv0 v0(vertex);
        AliESDv0 *pv0=&v0;
        AliExternalTrackParam bt(*esdTrackBach), *pbt=&bt;
        Double_t cascdca = PropagateToDCA(pv0,pbt,lESDevent,lMagneticField);
        
        fTreeCascVarDCACascDaughtersClassical = 1e+10;
        fTreeCascVarCascPropagationClassical = kFALSE;
        fTreeCascVarBachPropagationParameterClassical = -100;
        
        fTreeCascVarDecayX = -100;
        fTreeCascVarDecayY = -100;
        fTreeCascVarDecayZ = -100;
        fTreeCascVarCascCosPointingAngle = -100;
        
        fTreeCascVarInvMassXiMinus = -100;
        fTreeCascVarInvMassXiPlus = -100;
        fTreeCascVarInvMassOmegaMinus = -100;
        fTreeCascVarInvMassOmegaPlus = -100;
        
        //Check if propagation successful
        if (cascdca < 1e+32){
            fTreeCascVarDCACascDaughtersClassical = cascdca;
            fTreeCascVarCascPropagationClassical = kTRUE;
            fTreeCascVarBachPropagationParameterClassical = pbt->GetX();
            
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
        
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        //Step 1 IMPROVED: Propagate one to each other using improved techniques
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        //Status code
        fTreeCascVarCascPropagationImprovedStatus = 0; //started
        fTreeCascVarImprovedDecayX = -100;
        fTreeCascVarImprovedDecayY = -100;
        fTreeCascVarImprovedDecayZ = -100;
        fTreeCascVarImprovedCascCosPointingAngle = -100;
        fTreeCascVarImprovedCascDCAxyToPV = -100;
        fTreeCascVarImprovedCascDCAzToPV = -100;
        fTreeCascVarBachPropagationParameterImproved = -100;
        fTreeCascVarDCACascDaughtersImproved = -100;
        
        fTreeCascVarImprovedInvMassXiMinus = -100;
        fTreeCascVarImprovedInvMassXiPlus = -100;
        fTreeCascVarImprovedInvMassOmegaMinus = -100;
        fTreeCascVarImprovedInvMassOmegaPlus = -100;
        
        //Use copy
        AliESDv0 v0imp(vertex);
        AliESDv0 *pv0imp=&v0imp;
        AliExternalTrackParam btimp(*esdTrackBach), *pbtimp=&btimp;
        
        //Uncertainties: bachelor track as well as V0
        //Harmonic mean-like: errors as weights to combine the errors (approx)
        Double_t dy2V0 = TMath::Power(1/(esdTrackPos->GetSigmaY2()+1e-8)+1/(esdTrackNeg->GetSigmaY2()+1e-8)+1e-8,-1);
        Double_t dz2V0 = TMath::Power(1/(esdTrackPos->GetSigmaZ2()+1e-8)+1/(esdTrackNeg->GetSigmaZ2()+1e-8)+1e-8,-1);
        
        Double_t dy2=pbtimp->GetSigmaY2() + dy2V0;
        Double_t dz2=pbtimp->GetSigmaZ2() + dz2V0;
        Double_t dx2=dy2;
        
        //For testing purposes: try geometric criterion only, disregard uncertainties
        if( fkIfImprovedExtraPrecisionFactor<0.9){
            dx2 = fkIfImprovedExtraPrecisionFactor;
            dy2 = fkIfImprovedExtraPrecisionFactor;
            dz2 = fkIfImprovedExtraPrecisionFactor;
        }
        
        //Create dummy V0 track
        //V0 properties to get started
        Double_t xyz[3], pxpypz[3], cv[21];
        for(Int_t ii=0;ii<21;ii++) cv[ii]=0.0; //something small
        
        pv0imp->GetXYZ(xyz[0],xyz[1],xyz[2]);
        pv0imp->GetPxPyPz( pxpypz[0],pxpypz[1],pxpypz[2] );
        
        //Mockup track for V0 trajectory (no covariance)
        AliExternalTrackParam lV0TrajObject(xyz,pxpypz,cv,+1), *hV0Traj = &lV0TrajObject;
        hV0Traj->ResetCovariance(1); //won't use
        
        Double_t p1[8]; pbtimp->GetHelixParameters(p1,lMagneticField);
        p1[6]=TMath::Sin(p1[2]); p1[7]=TMath::Cos(p1[2]);
        Double_t p2[8]; hV0Traj->GetHelixParameters(p2,0.0); //p2[4]=0 -> no curvature (fine, predicted in Evaluate)
        p2[6]=TMath::Sin(p2[2]); p2[7]=TMath::Cos(p2[2]);
        
        Double_t r1[3],g1[3],gg1[3]; Double_t t1=0.;
        Evaluate(p1,t1,r1,g1,gg1);
        Double_t r2[3],g2[3],gg2[3]; Double_t t2=0.;
        Evaluate(p2,t2,r2,g2,gg2);
        
        Double_t dx=r2[0]-r1[0], dy=r2[1]-r1[1], dz=r2[2]-r1[2];
        Double_t dm=dx*dx/dx2 + dy*dy/dy2 + dz*dz/dz2;
        
        Int_t max=fMaxIterationsCascadeDCA;
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
            
            if (TMath::Abs(dt1)/(TMath::Abs(t1)+1.e-3) < fPrecisionCutoffCascadeDCA)
                if (TMath::Abs(dt2)/(TMath::Abs(t2)+1.e-3) < fPrecisionCutoffCascadeDCA) {
                    if ((gt1*gt1+gt2*gt2) > 1.e-4/dy2/dy2){
                        AliDebug(1," stopped at not a stationary point !");
                        //Count not stationary point
                        fTreeCascVarCascPropagationImprovedStatus += 1; //not stationary
                    }
                    Double_t lmb=h11+h22; lmb=lmb-TMath::Sqrt(lmb*lmb-4*det);
                    if (lmb < 0.){
                        
                        AliDebug(1," stopped at not a minimum !");
                        //Count stopped at not a minimum
                        fTreeCascVarCascPropagationImprovedStatus += 10; //not minimum
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
                    fTreeCascVarCascPropagationImprovedStatus += 100; //overshoot!
                }
            }
            dm=dd;
            
            t1+=dt1;
            t2+=dt2;
            
        }
        
        if (max<=0){
            AliDebug(1," too many iterations !");
        }
        
        fTreeCascVarCascPropagationImprovedIterations = max;
        
        Double_t cs=TMath::Cos(pbtimp->GetAlpha());
        Double_t sn=TMath::Sin(pbtimp->GetAlpha());
        Double_t xthis=r1[0]*cs + r1[1]*sn;
        fTreeCascVarBachPropagationParameterImproved = xthis;
        
        //Memory cleanup
        hV0Traj->Delete();
        hV0Traj=0x0;
        
        //Propagate bachelor to the point of DCA
        if (!pbtimp->PropagateTo(xthis,lMagneticField)) {
            //AliWarning(" propagation failed !";
            //HOLY CRAP, propagation failed!
            fTreeCascVarCascPropagationImprovedStatus += 1000; //failure: add 1000
        }
        
        if( fTreeCascVarCascPropagationImprovedStatus < 999 ){
            
            //V0 distance to bachelor: the desired distance
            Double_t rBachDCAPt[3]; pbtimp->GetXYZ(rBachDCAPt);
            fTreeCascVarDCACascDaughtersImproved = pv0imp->GetD(rBachDCAPt[0],rBachDCAPt[1],rBachDCAPt[2]);
            
            //Construct cascade
            AliESDcascade cascadeimproved(*pv0imp,*pbtimp,lCascBachTrackArray[iCasc]);
            
            //Decay Position
            Double_t xcascimp,ycascimp,zcascimp; cascadeimproved.GetXYZcascade(xcascimp,ycascimp,zcascimp);
            fTreeCascVarImprovedDecayX = xcascimp;
            fTreeCascVarImprovedDecayY = ycascimp;
            fTreeCascVarImprovedDecayZ = zcascimp;
            
            fTreeCascVarImprovedCascCosPointingAngle = cascadeimproved.GetCascadeCosineOfPointingAngle(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2]);
            
            
            //V0 properties to get started
            Double_t xyzCascade[3], pxpypzCascade[3], cvCascade[21];
            for(Int_t ii=0;ii<21;ii++) cvCascade[ii]=0.0; //something small
            
            Int_t lChargeCascade = -1;
            if ( fTreeCascVarBachSign > 0 ) lChargeCascade = +1;
            
            //Don't forget: set properties, please...
            cascadeimproved.GetXYZcascade( xyzCascade[0],  xyzCascade[1], xyzCascade[2] );
            cascadeimproved.GetPxPyPz( pxpypzCascade[0], pxpypzCascade[1], pxpypzCascade[2] );
            
            AliExternalTrackParam lCascTrajObject(xyzCascade,pxpypzCascade,cvCascade,lChargeCascade), *hCascTraj = &lCascTrajObject;
            
            fTreeCascVarImprovedCascDCAxyToPV = TMath::Abs(hCascTraj->GetD(lBestPrimaryVtxPos[0],
                                                                            lBestPrimaryVtxPos[1],
                                                                            lMagneticField) );
            Float_t dzcascade[2];
            hCascTraj->GetDZ(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2], lMagneticField, dzcascade );
            fTreeCascVarImprovedCascDCAzToPV = dzcascade[1];
            
            Double_t lV0quality  = 0.;
            cascadeimproved.ChangeMassHypothesis(lV0quality , 3312);
            fTreeCascVarImprovedInvMassXiMinus = cascadeimproved.GetEffMassXi();
            cascadeimproved.ChangeMassHypothesis(lV0quality ,-3312);
            fTreeCascVarImprovedInvMassXiPlus = cascadeimproved.GetEffMassXi();
            cascadeimproved.ChangeMassHypothesis(lV0quality , 3334);
            fTreeCascVarImprovedInvMassOmegaMinus = cascadeimproved.GetEffMassXi();
            cascadeimproved.ChangeMassHypothesis(lV0quality ,-3334);
            fTreeCascVarImprovedInvMassOmegaPlus = cascadeimproved.GetEffMassXi();
        }
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        
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
        Int_t lParticleMotherPDG = lParticleMother->GetPdgCode();
        
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
    
    // Post output data.
    PostData(1, fListHist    );
    PostData(2, fListV0      );
    PostData(3, fListCascade );
    PostData(4, fTreeEvent   );
    PostData(5, fTreeV0      );
    PostData(6, fTreeCascade );
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
        if ( i>0 )
            lV0Result[lNV0]->InitializeFeeddownMatrix( lPtbinnumb, lPtbinlimits, lPtbinnumbCascade, lPtbinlimitsCascade, lCentbinnumb, lCentbinlimits );
        
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
void AliAnalysisTaskStrEffStudy::AddTopologicalQACascade(Int_t lRecNumberOfSteps)
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
void AliAnalysisTaskStrEffStudy::AddStandardV0Configuration()
//Meant to add some standard V0 analysis Configuration + its corresponding systematics
{
    // STEP 1: Decide on binning (needed to improve on memory consumption)
    
    // pT binning
    Double_t lPtbinlimitsV0[] ={0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.5, 8.0, 10, 12, 14, 15, 17, 20};
    Long_t lPtbinnumbV0 = sizeof(lPtbinlimitsV0)/sizeof(Double_t) - 1;
    Double_t lPtbinlimitsXi[] ={0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.5, 8.0, 10, 12, 14, 17, 20};
    Long_t lPtbinnumbXi = sizeof(lPtbinlimitsXi)/sizeof(Double_t) - 1;
    
    // centrality binning
    Double_t lCentbinlimitsV0[] = {0, 1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 85, 90};
    Long_t lCentbinnumbV0 = sizeof(lCentbinlimitsV0)/sizeof(Double_t) - 1;
    
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
        lV0Result[lNV0] = new AliV0Result( Form("%s_Central",lParticleName[i].Data() ),lMassHypoV0[i],"",lCentbinnumbV0,lCentbinlimitsV0, lPtbinnumbV0,lPtbinlimitsV0);
        if ( i!=0 ) lV0Result[lNV0] -> InitializeFeeddownMatrix( lPtbinnumbV0,lPtbinlimitsV0, lPtbinnumbXi,lPtbinlimitsXi, lCentbinnumbV0,lCentbinlimitsV0);
        
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
        if ( i!=0 ) lV0Result[lNV0] -> InitializeFeeddownMatrix( lPtbinnumbV0,lPtbinlimitsV0, lPtbinnumbXi,lPtbinlimitsXi, lCentbinnumbV0,lCentbinlimitsV0);
        
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
void AliAnalysisTaskStrEffStudy::AddStandardCascadeConfiguration()
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
        
        //This is MC: generate profile for G3/F (if ever needed)
        lCascadeResult[lN] -> InitializeProtonProfile(lPtbinnumb,lPtbinlimits);
        
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
        
        //This is MC: generate profile for G3/F (if ever needed)
        lCascadeResult[lN] -> InitializeProtonProfile(lPtbinnumb,lPtbinlimits);
        
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
        
        
        //This is MC: generate profile for G3/F (if ever needed)
        lCascadeResult[lN] -> InitializeProtonProfile(lPtbinnumb,lPtbinlimits);
        
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
void AliAnalysisTaskStrEffStudy::AddCascadeConfiguration276TeV()
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
        lCascadeResult[lN] -> InitializeProtonProfile(lPtbinnumb,lPtbinlimits);
        
        //Setters for V0 Cuts
        lCascadeResult[lN]->SetCutDCANegToPV            ( 0.1    ) ;
        lCascadeResult[lN]->SetCutDCAPosToPV            ( 0.1    ) ;
        lCascadeResult[lN]->SetCutDCAV0Daughters        ( 0.8    ) ;
        lCascadeResult[lN]->SetCutV0CosPA               ( 0.998  ) ;
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
    Float_t dcaToVertexXY = b[0];
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
            Error("PropagateToDCA","Propagation failed !");
            return 1.e+33;
        }
    }
    
    if( fkDoImprovedCascadeVertexFinding ){
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
        
        //For testing purposes
        dx2 = dx2/fkIfImprovedExtraPrecisionFactor;
        dy2 = dy2/fkIfImprovedExtraPrecisionFactor;
        dz2 = dz2/fkIfImprovedExtraPrecisionFactor;
        
        //Create dummy V0 track
        //V0 properties to get started
        Double_t xyz[3], pxpypz[3], cv[21];
        for(Int_t ii=0;ii<21;ii++) cv[ii]=0.0; //something small
        
        v->GetXYZ(xyz[0],xyz[1],xyz[2]);
        v->GetPxPyPz( pxpypz[0],pxpypz[1],pxpypz[2] );
        
        //Mockup track for V0 trajectory (no covariance)
        AliExternalTrackParam *hV0Traj = new AliExternalTrackParam(xyz,pxpypz,cv,+1);
        hV0Traj->ResetCovariance(1); //won't use
        
        Double_t p1[8]; t->GetHelixParameters(p1,b);
        p1[6]=TMath::Sin(p1[2]); p1[7]=TMath::Cos(p1[2]);
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
        
        Double_t cs=TMath::Cos(t->GetAlpha());
        Double_t sn=TMath::Sin(t->GetAlpha());
        Double_t xthis=r1[0]*cs + r1[1]*sn;
        
        //Memory cleanup
        hV0Traj->Delete();
        hV0Traj=0x0;
        
        //Propagate bachelor to the point of DCA
        if (!t->PropagateTo(xthis,b)) {
            //AliWarning(" propagation failed !";
            return 1e+33;
        }
        
        
        //V0 distance to bachelor: the desired distance
        Double_t rBachDCAPt[3]; t->GetXYZ(rBachDCAPt);
        dca = v->GetD(rBachDCAPt[0],rBachDCAPt[1],rBachDCAPt[2]);
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
    Double_t alpha=t1->GetAlpha(), cs=TMath::Cos(alpha), sn=TMath::Sin(alpha);
    Double_t tmp[3];
    t1->GetPxPyPz(tmp);
    Double_t px1=tmp[0], py1=tmp[1], pz1=tmp[2];
    t1->GetXYZ(tmp);
    Double_t  x1=tmp[0],  y1=tmp[1],  z1=tmp[2];
    const Double_t ss=0.0005*0.0005;//a kind of a residual misalignment precision
    Double_t sx1=sn*sn*t1->GetSigmaY2()+ss, sy1=cs*cs*t1->GetSigmaY2()+ss;
    return sx1;
}