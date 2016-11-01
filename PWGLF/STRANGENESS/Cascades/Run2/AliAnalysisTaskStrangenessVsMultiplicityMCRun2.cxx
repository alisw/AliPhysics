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
//#include "AliLog.h"

#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"
#include "AliLightV0vertexer.h"
#include "AliLightCascadeVertexer.h"
#include "AliESDpid.h"
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

#include "AliCFContainer.h"
#include "AliMultiplicity.h"
#include "AliAODMCParticle.h"
#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliESDUtils.h"
#include "AliGenEventHeader.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliV0Result.h"
#include "AliCascadeResult.h"
#include "AliAnalysisTaskStrangenessVsMultiplicityMCRun2.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskStrangenessVsMultiplicityMCRun2)

AliAnalysisTaskStrangenessVsMultiplicityMCRun2::AliAnalysisTaskStrangenessVsMultiplicityMCRun2()
    : AliAnalysisTaskSE(), fListHist(0), fListV0(0), fListCascade(0), fTreeEvent(0), fTreeV0(0), fTreeCascade(0), fPIDResponse(0), fESDtrackCuts(0), fUtils(0), fRand(0),

//---> Flags controlling Event Tree output
        fkSaveEventTree    ( kTRUE ), //no downscaling in this tree so far

//---> Flags controlling V0 TTree output
        fkSaveV0Tree       ( kTRUE ),
        fkDownScaleV0      ( kTRUE  ),
        fDownScaleFactorV0 ( 0.001  ),
        fkPreselectDedx ( kFALSE ),
        fkPreselectPID  ( kTRUE  ),
        fkDebugWrongPIDForTracking ( kFALSE ),

//---> Flags controlling Cascade TTree output
        fkSaveCascadeTree       ( kTRUE  ),
        fkDownScaleCascade      ( kTRUE  ),
        fDownScaleFactorCascade ( 0.001  ),

//---> Flags controlling Vertexers
        fkRunVertexers    ( kFALSE ),
        fkUseLightVertexer ( kTRUE ),

//---> Flag controlling trigger selection
        fTrigType(AliVEvent::kMB),

//---> Variables for fTreeEvent
      fCentrality(0),

//---> Variables for fTreeV0
      fTreeVariableChi2V0(0),
      fTreeVariableDcaV0Daughters(0),
      fTreeVariableDcaV0ToPrimVertex(0),
      fTreeVariableDcaPosToPrimVertex(0),
      fTreeVariableDcaNegToPrimVertex(0),
      fTreeVariableV0CosineOfPointingAngle(0),
      fTreeVariableV0Radius(0),
      fTreeVariablePt(0),
      fTreeVariablePtMC(0),
      fTreeVariableRapK0Short(0),
      fTreeVariableRapLambda(0),
      fTreeVariableRapMC(0),
      fTreeVariableInvMassK0s(0),
      fTreeVariableInvMassLambda(0),
      fTreeVariableInvMassAntiLambda(0),
      fTreeVariableAlphaV0(0),
      fTreeVariablePtArmV0(0),
      fTreeVariableNegEta(0),
      fTreeVariablePosEta(0),

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

      fTreeVariableDistOverTotMom(0),
      fTreeVariableLeastNbrCrossedRows(0),
      fTreeVariableLeastRatioCrossedRowsOverFindable(0),

      fTreeVariableCentrality(0),
//MC Variables
      fTreeVariablePtMother(0),
      fTreeVariableRapMother(0),
      fTreeVariablePID(0),
      fTreeVariablePIDPositive(0),
      fTreeVariablePIDNegative(0),
      fTreeVariablePIDMother(0),
      fTreeVariablePrimaryStatus(0),
      fTreeVariablePrimaryStatusMother(0),

//---> Variables for fTreeCascade
      fTreeCascVarCharge(0),
      fTreeCascVarMassAsXi(0),
      fTreeCascVarMassAsOmega(0),
      fTreeCascVarPt(0),
      fTreeCascVarPtMC(0),
      fTreeCascVarRapXi(0),
      fTreeCascVarRapOmega(0),
      fTreeCascVarRapMC(0),
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
      fTreeCascVarV0CosPointingAngle(0),
      fTreeCascVarV0CosPointingAngleSpecial(0),
      fTreeCascVarV0Radius(0),
      fTreeCascVarLeastNbrClusters(0),
      fTreeCascVarDistOverTotMom(0),
      fTreeCascVarNegNSigmaPion(0),
      fTreeCascVarNegNSigmaProton(0),
      fTreeCascVarPosNSigmaPion(0),
      fTreeCascVarPosNSigmaProton(0),
      fTreeCascVarBachNSigmaPion(0),
      fTreeCascVarBachNSigmaKaon(0),
fTreeCascVarPosPIDForTracking(-1),
fTreeCascVarNegPIDForTracking(-1),
fTreeCascVarBachPIDForTracking(-1),
fTreeCascVarPosdEdx(-1),
fTreeCascVarNegdEdx(-1),
fTreeCascVarBachdEdx(-1),
fTreeCascVarPosInnerP(-1),
fTreeCascVarNegInnerP(-1),
fTreeCascVarBachInnerP(-1),
fTreeCascVarCentrality(0),
fTreeCascVarPID(0),
fTreeCascVarSwappedPID(0),
fTreeCascVarIsPhysicalPrimary(0),

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

AliAnalysisTaskStrangenessVsMultiplicityMCRun2::AliAnalysisTaskStrangenessVsMultiplicityMCRun2(Bool_t lSaveEventTree, Bool_t lSaveV0Tree, Bool_t lSaveCascadeTree, const char *name, TString lExtraOptions)
    : AliAnalysisTaskSE(name), fListHist(0), fListV0(0), fListCascade(0), fTreeEvent(0), fTreeV0(0), fTreeCascade(0), fPIDResponse(0), fESDtrackCuts(0), fUtils(0), fRand(0),

//---> Flags controlling Event Tree output
fkSaveEventTree    ( kFALSE ), //no downscaling in this tree so far

//---> Flags controlling V0 TTree output
fkSaveV0Tree       ( kTRUE ),
fkDownScaleV0      ( kTRUE  ),
fDownScaleFactorV0 ( 0.001  ),
fkPreselectDedx ( kFALSE ),
fkPreselectPID  ( kTRUE  ),
fkDebugWrongPIDForTracking ( kFALSE ), //also for cascades...

//---> Flags controlling Cascade TTree output
fkSaveCascadeTree       ( kTRUE  ),
fkDownScaleCascade      ( kTRUE  ),
fDownScaleFactorCascade ( 0.001  ),

//---> Flags controlling Vertexers
fkRunVertexers    ( kFALSE ),
fkUseLightVertexer ( kTRUE ),

//---> Flag controlling trigger selection
fTrigType(AliVEvent::kMB),

//---> Variables for fTreeEvent
fCentrality(0),

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

      fTreeVariableDistOverTotMom(0),
      fTreeVariableLeastNbrCrossedRows(0),
      fTreeVariableLeastRatioCrossedRowsOverFindable(0),

      fTreeVariableCentrality(0),

//MC Variables
fTreeVariablePtMother(0),
fTreeVariableRapMother(0),
fTreeVariablePID(0),
fTreeVariablePIDPositive(0),
fTreeVariablePIDNegative(0),
fTreeVariablePIDMother(0),
fTreeVariablePrimaryStatus(0),
fTreeVariablePrimaryStatusMother(0),

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
      fTreeCascVarV0CosPointingAngle(0),
      fTreeCascVarV0CosPointingAngleSpecial(0),
      fTreeCascVarV0Radius(0),
      fTreeCascVarLeastNbrClusters(0),
      fTreeCascVarDistOverTotMom(0),
      fTreeCascVarNegNSigmaPion(0),
      fTreeCascVarNegNSigmaProton(0),
      fTreeCascVarPosNSigmaPion(0),
      fTreeCascVarPosNSigmaProton(0),
      fTreeCascVarBachNSigmaPion(0),
      fTreeCascVarBachNSigmaKaon(0),
fTreeCascVarPosPIDForTracking(-1),
fTreeCascVarNegPIDForTracking(-1),
fTreeCascVarBachPIDForTracking(-1),
fTreeCascVarPosdEdx(-1),
fTreeCascVarNegdEdx(-1),
fTreeCascVarBachdEdx(-1),
fTreeCascVarPosInnerP(-1),
fTreeCascVarNegInnerP(-1),
fTreeCascVarBachInnerP(-1),
fTreeCascVarCentrality(0),
fTreeCascVarPID(0),
fTreeCascVarSwappedPID(0),
fTreeCascVarIsPhysicalPrimary(0),
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
    if ( lExtraOptions.Contains("A") ) fkDebugWrongPIDForTracking = kTRUE;

}


AliAnalysisTaskStrangenessVsMultiplicityMCRun2::~AliAnalysisTaskStrangenessVsMultiplicityMCRun2()
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
void AliAnalysisTaskStrangenessVsMultiplicityMCRun2::UserCreateOutputObjects()
{
    //------------------------------------------------
    // fTreeEvent: EbyE information
    //------------------------------------------------
    if(fkSaveEventTree){
        fTreeEvent = new TTree("fTreeEvent","Event");
        //Branch Definitions
        fTreeEvent->Branch("fCentrality",&fCentrality,"fCentrality/F");
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
        fTreeV0->Branch("fTreeVariablePtMC",&fTreeVariablePtMC,"fTreeVariablePtMC/F");
        fTreeV0->Branch("fTreeVariableRapK0Short",&fTreeVariableRapK0Short,"fTreeVariableRapK0Short/F");
        fTreeV0->Branch("fTreeVariableRapLambda",&fTreeVariableRapLambda,"fTreeVariableRapLambda/F");
        fTreeV0->Branch("fTreeVariableRapMC",&fTreeVariableRapMC,"fTreeVariableRapMC/F");
        fTreeV0->Branch("fTreeVariableInvMassK0s",&fTreeVariableInvMassK0s,"fTreeVariableInvMassK0s/F");
        fTreeV0->Branch("fTreeVariableInvMassLambda",&fTreeVariableInvMassLambda,"fTreeVariableInvMassLambda/F");
        fTreeV0->Branch("fTreeVariableInvMassAntiLambda",&fTreeVariableInvMassAntiLambda,"fTreeVariableInvMassAntiLambda/F");
        fTreeV0->Branch("fTreeVariableV0CosineOfPointingAngle",&fTreeVariableV0CosineOfPointingAngle,"fTreeVariableV0CosineOfPointingAngle/F");
        fTreeV0->Branch("fTreeVariableAlphaV0",&fTreeVariableAlphaV0,"fTreeVariableAlphaV0/F");
        fTreeV0->Branch("fTreeVariablePtArmV0",&fTreeVariablePtArmV0,"fTreeVariablePtArmV0/F");
        fTreeV0->Branch("fTreeVariableLeastNbrCrossedRows",&fTreeVariableLeastNbrCrossedRows,"fTreeVariableLeastNbrCrossedRows/I");
        fTreeV0->Branch("fTreeVariableLeastRatioCrossedRowsOverFindable",&fTreeVariableLeastRatioCrossedRowsOverFindable,"fTreeVariableLeastRatioCrossedRowsOverFindable/F");
        fTreeV0->Branch("fTreeVariableDistOverTotMom",&fTreeVariableDistOverTotMom,"fTreeVariableDistOverTotMom/F");
        fTreeV0->Branch("fTreeVariableNSigmasPosProton",&fTreeVariableNSigmasPosProton,"fTreeVariableNSigmasPosProton/F");
        fTreeV0->Branch("fTreeVariableNSigmasPosPion",&fTreeVariableNSigmasPosPion,"fTreeVariableNSigmasPosPion/F");
        fTreeV0->Branch("fTreeVariableNSigmasNegProton",&fTreeVariableNSigmasNegProton,"fTreeVariableNSigmasNegProton/F");
        fTreeV0->Branch("fTreeVariableNSigmasNegPion",&fTreeVariableNSigmasNegPion,"fTreeVariableNSigmasNegPion/F");
        fTreeV0->Branch("fTreeVariableNegEta",&fTreeVariableNegEta,"fTreeVariableNegEta/F");
        fTreeV0->Branch("fTreeVariablePosEta",&fTreeVariablePosEta,"fTreeVariablePosEta/F");
        //-----------MULTIPLICITY-INFO--------------------
        fTreeV0->Branch("fTreeVariableCentrality",&fTreeVariableCentrality,"fTreeVariableCentrality/F");
        //------------------------------------------------
        if ( fkDebugWrongPIDForTracking ){
            fTreeV0->Branch("fTreeVariablePosPIDForTracking",&fTreeVariablePosPIDForTracking,"fTreeVariablePosPIDForTracking/I");
            fTreeV0->Branch("fTreeVariableNegPIDForTracking",&fTreeVariableNegPIDForTracking,"fTreeVariableNegPIDForTracking/I");
            fTreeV0->Branch("fTreeVariablePosdEdx",&fTreeVariablePosdEdx,"fTreeVariablePosdEdx/F");
            fTreeV0->Branch("fTreeVariableNegdEdx",&fTreeVariableNegdEdx,"fTreeVariableNegdEdx/F");
            fTreeV0->Branch("fTreeVariablePosInnerP",&fTreeVariablePosInnerP,"fTreeVariablePosInnerP/F");
            fTreeV0->Branch("fTreeVariableNegInnerP",&fTreeVariableNegInnerP,"fTreeVariableNegInnerP/F");
        }
        //-----------MC Exclusive info--------------------
        fTreeV0->Branch("fTreeVariablePtMother",&fTreeVariablePtMother,"fTreeVariablePtMother/F");
        fTreeV0->Branch("fTreeVariableRapMother",&fTreeVariableRapMother,"fTreeVariableRapMother/F");
        fTreeV0->Branch("fTreeVariablePID",&fTreeVariablePID,"fTreeVariablePID/I");
        fTreeV0->Branch("fTreeVariablePIDPositive",&fTreeVariablePIDPositive,"fTreeVariablePIDPositive/I");
        fTreeV0->Branch("fTreeVariablePIDNegative",&fTreeVariablePIDNegative,"fTreeVariablePIDNegative/I");
        fTreeV0->Branch("fTreeVariablePIDMother",&fTreeVariablePIDMother,"fTreeVariablePIDMother/I");
        fTreeV0->Branch("fTreeVariablePrimaryStatus",&fTreeVariablePrimaryStatus,"fTreeVariablePrimaryStatus/I");
        fTreeV0->Branch("fTreeVariablePrimaryStatusMother",&fTreeVariablePrimaryStatusMother,"fTreeVariablePrimaryStatusMother/I");
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
        fTreeCascade->Branch("fTreeCascVarPtMC",&fTreeCascVarPtMC,"fTreeCascVarPtMC/F");
        fTreeCascade->Branch("fTreeCascVarRapXi",&fTreeCascVarRapXi,"fTreeCascVarRapXi/F");
        fTreeCascade->Branch("fTreeCascVarRapOmega",&fTreeCascVarRapOmega,"fTreeCascVarRapOmega/F");
        fTreeCascade->Branch("fTreeCascVarRapMC",&fTreeCascVarRapMC,"fTreeCascVarRapMC/F");
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
        fTreeCascade->Branch("fTreeCascVarV0CosPointingAngle",&fTreeCascVarV0CosPointingAngle,"fTreeCascVarV0CosPointingAngle/F");
        fTreeCascade->Branch("fTreeCascVarV0CosPointingAngleSpecial",&fTreeCascVarV0CosPointingAngleSpecial,"fTreeCascVarV0CosPointingAngleSpecial/F");
        fTreeCascade->Branch("fTreeCascVarV0Radius",&fTreeCascVarV0Radius,"fTreeCascVarV0Radius/F");
        fTreeCascade->Branch("fTreeCascVarLeastNbrClusters",&fTreeCascVarLeastNbrClusters,"fTreeCascVarLeastNbrClusters/I");
        //-----------MULTIPLICITY-INFO--------------------
        fTreeCascade->Branch("fTreeCascVarCentrality",&fTreeCascVarCentrality,"fTreeCascVarCentrality/F");
        //-----------DECAY-LENGTH-INFO--------------------
        fTreeCascade->Branch("fTreeCascVarDistOverTotMom",&fTreeCascVarDistOverTotMom,"fTreeCascVarDistOverTotMom/F");
        //------------------------------------------------
        fTreeCascade->Branch("fTreeCascVarNegNSigmaPion",&fTreeCascVarNegNSigmaPion,"fTreeCascVarNegNSigmaPion/F");
        fTreeCascade->Branch("fTreeCascVarNegNSigmaProton",&fTreeCascVarNegNSigmaProton,"fTreeCascVarNegNSigmaProton/F");
        fTreeCascade->Branch("fTreeCascVarPosNSigmaPion",&fTreeCascVarPosNSigmaPion,"fTreeCascVarPosNSigmaPion/F");
        fTreeCascade->Branch("fTreeCascVarPosNSigmaProton",&fTreeCascVarPosNSigmaProton,"fTreeCascVarPosNSigmaProton/F");
        fTreeCascade->Branch("fTreeCascVarBachNSigmaPion",&fTreeCascVarBachNSigmaPion,"fTreeCascVarBachNSigmaPion/F");
        fTreeCascade->Branch("fTreeCascVarBachNSigmaKaon",&fTreeCascVarBachNSigmaKaon,"fTreeCascVarBachNSigmaKaon/F");
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
        }
        //-----------MC Exclusive info--------------------
        fTreeCascade->Branch("fTreeCascVarIsPhysicalPrimary",&fTreeCascVarIsPhysicalPrimary,"fTreeCascVarIsPhysicalPrimary/I");
        fTreeCascade->Branch("fTreeCascVarPID",&fTreeCascVarPID,"fTreeCascVarPID/I");
        fTreeCascade->Branch("fTreeCascVarSwappedPID",&fTreeCascVarSwappedPID,"fTreeCascVarSwappedPID/I");
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
        // used to fill the first 8 integers of the seed array.
        fRand->SetSeed(0);
    }

    //------------------------------------------------
    // V0 Multiplicity Histograms
    //------------------------------------------------

    // Create histograms
    fListHist = new TList();
    fListHist->SetOwner();  // See http://root.cern.ch/root/html/TCollection.html#TCollection:SetOwner

    if(! fHistEventCounter ) {
        //Histogram Output: Event-by-Event
        fHistEventCounter = new TH1D( "fHistEventCounter", ";Evt. Sel. Step;Count",2,0,2);
        fHistEventCounter->GetXaxis()->SetBinLabel(1, "Processed");
        fHistEventCounter->GetXaxis()->SetBinLabel(2, "Selected");
        fListHist->Add(fHistEventCounter);
    }
    
    
    if(! fHistCentrality ) {
        //Histogram Output: Event-by-Event
        fHistCentrality = new TH1D( "fHistCentrality", ";Centrality;Event Count",100,0,100);
        fListHist->Add(fHistCentrality);
    }
    
    if(! fHistGeneratedPtVsYVsCentralityK0Short ) {
        //Histogram Output: Efficiency Denominator
        fHistGeneratedPtVsYVsCentralityK0Short = new TH3D( "fHistGeneratedPtVsYVsCentralityK0Short", ";pT;y;centrality",200,0,20,20,-1.0,1.0,100,0,100);
        fListHist->Add(fHistGeneratedPtVsYVsCentralityK0Short);
    }
    if(! fHistGeneratedPtVsYVsCentralityLambda ) {
        //Histogram Output: Efficiency Denominator
        fHistGeneratedPtVsYVsCentralityLambda = new TH3D( "fHistGeneratedPtVsYVsCentralityLambda", ";pT;y;centrality",200,0,20,20,-1.0,1.0,100,0,100);
        fListHist->Add(fHistGeneratedPtVsYVsCentralityLambda);
    }
    if(! fHistGeneratedPtVsYVsCentralityAntiLambda ) {
        //Histogram Output: Efficiency Denominator
        fHistGeneratedPtVsYVsCentralityAntiLambda = new TH3D( "fHistGeneratedPtVsYVsCentralityAntiLambda", ";pT;y;centrality",200,0,20,20,-1.0,1.0,100,0,100);
        fListHist->Add(fHistGeneratedPtVsYVsCentralityAntiLambda);
    }
    if(! fHistGeneratedPtVsYVsCentralityXiMinus ) {
        //Histogram Output: Efficiency Denominator
        fHistGeneratedPtVsYVsCentralityXiMinus = new TH3D( "fHistGeneratedPtVsYVsCentralityXiMinus", ";pT;y;centrality",200,0,20,20,-1.0,1.0,100,0,100);
        fListHist->Add(fHistGeneratedPtVsYVsCentralityXiMinus);
    }
    if(! fHistGeneratedPtVsYVsCentralityXiPlus ) {
        //Histogram Output: Efficiency Denominator
        fHistGeneratedPtVsYVsCentralityXiPlus = new TH3D( "fHistGeneratedPtVsYVsCentralityXiPlus", ";pT;y;centrality",200,0,20,20,-1.0,1.0,100,0,100);
        fListHist->Add(fHistGeneratedPtVsYVsCentralityXiPlus);
    }
    if(! fHistGeneratedPtVsYVsCentralityOmegaMinus ) {
        //Histogram Output: Efficiency Denominator
        fHistGeneratedPtVsYVsCentralityOmegaMinus = new TH3D( "fHistGeneratedPtVsYVsCentralityOmegaMinus", ";pT;y;centrality",200,0,20,20,-1.0,1.0,100,0,100);
        fListHist->Add(fHistGeneratedPtVsYVsCentralityOmegaMinus);
    }
    if(! fHistGeneratedPtVsYVsCentralityOmegaPlus ) {
        //Histogram Output: Efficiency Denominator
        fHistGeneratedPtVsYVsCentralityOmegaPlus = new TH3D( "fHistGeneratedPtVsYVsCentralityOmegaPlus", ";pT;y;centrality",200,0,20,20,-1.0,1.0,100,0,100);
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
    if(fkSaveEventTree)    PostData(4, fTreeEvent   );
    if(fkSaveV0Tree)       PostData(5, fTreeV0      );
    if(fkSaveCascadeTree)  PostData(6, fTreeCascade );

}// end UserCreateOutputObjects


//________________________________________________________________________
void AliAnalysisTaskStrangenessVsMultiplicityMCRun2::UserExec(Option_t *)
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
    
    fCentrality = lPercentile;
    
    if( lEvSelCode != 0 ) {
        PostData(1, fListHist    );
        PostData(2, fListV0      );
        PostData(3, fListCascade );
        if( fkSaveEventTree   ) PostData(4, fTreeEvent   );
        if( fkSaveV0Tree      ) PostData(5, fTreeV0      );
        if( fkSaveCascadeTree ) PostData(6, fTreeCascade );
        return;
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

    Int_t nv0s = 0;
    nv0s = lESDevent->GetNumberOfV0s();

    for (Int_t iV0 = 0; iV0 < nv0s; iV0++) //extra-crazy test
    {   // This is the begining of the V0 loop
        AliESDv0 *v0 = ((AliESDEvent*)lESDevent)->GetV0(iV0);
        if (!v0) continue;

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
        
        const AliExternalTrackParam *innernegv0=nTrack->GetInnerParam();
        const AliExternalTrackParam *innerposv0=pTrack->GetInnerParam();
        Float_t lThisPosInnerP = -1;
        Float_t lThisNegInnerP = -1;
        if(innerposv0)  { lThisPosInnerP  = innerposv0 ->GetP(); }
        if(innernegv0)  { lThisNegInnerP  = innernegv0 ->GetP(); }
        Float_t lThisPosdEdx = pTrack -> GetTPCsignal();
        Float_t lThisNegdEdx = nTrack -> GetTPCsignal();
        
        fTreeVariablePosdEdx = lThisPosdEdx;
        fTreeVariableNegdEdx = lThisNegdEdx;
        
        fTreeVariablePosInnerP = lThisPosInnerP;
        fTreeVariableNegInnerP = lThisNegInnerP;
        
        //Daughter Eta for Eta selection, afterwards
        fTreeVariableNegEta = nTrack->Eta();
        fTreeVariablePosEta = pTrack->Eta();

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


        if ( ( ( pTrack->GetTPCClusterInfo(2,1) ) < 70 ) || ( ( nTrack->GetTPCClusterInfo(2,1) ) < 70 ) ) continue;

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
        fTreeVariableCentrality = fCentrality;

        //===============================================
        // V0 Monte Carlo Association starts here
        //===============================================
        
        //---> Set Everything to "I don't know" before starting
        //---> This is strictly necessary!
        
        fTreeVariablePIDPositive = 0;
        fTreeVariablePIDNegative = 0;
        
        fTreeVariablePtMother = -1;
        fTreeVariableRapMother = -100;
        fTreeVariablePtMC = -1;
        fTreeVariableRapMC = -100;
        
        fTreeVariablePID = -1;
        fTreeVariablePIDMother = -1;
        
        fTreeVariablePrimaryStatus = 0;
        fTreeVariablePrimaryStatusMother = 0;
        
        //fTreeVariablePosTransvMomentumMC = -1;
        //fTreeVariableNegTransvMomentumMC = -1;
        
        Int_t lblPosV0Dghter = (Int_t) TMath::Abs( pTrack->GetLabel() );
        Int_t lblNegV0Dghter = (Int_t) TMath::Abs( nTrack->GetLabel() );
        
        TParticle* mcPosV0Dghter = lMCstack->Particle( lblPosV0Dghter );
        TParticle* mcNegV0Dghter = lMCstack->Particle( lblNegV0Dghter );
        
        Int_t lPIDPositive = mcPosV0Dghter -> GetPdgCode();
        Int_t lPIDNegative = mcNegV0Dghter -> GetPdgCode();
        
        //fTreeVariablePosTransvMomentumMC = mcPosV0Dghter->Pt();
        //fTreeVariableNegTransvMomentumMC = mcNegV0Dghter->Pt();
        
        fTreeVariablePIDPositive = lPIDPositive;
        fTreeVariablePIDNegative = lPIDNegative;
        
        Int_t lblMotherPosV0Dghter = mcPosV0Dghter->GetFirstMother() ;
        Int_t lblMotherNegV0Dghter = mcNegV0Dghter->GetFirstMother();
        
        if( lblMotherPosV0Dghter == lblMotherNegV0Dghter && lblMotherPosV0Dghter > -1 ) {
            //either label is fine, they're equal at this stage
            TParticle* pThisV0 = lMCstack->Particle( lblMotherPosV0Dghter );
            //Set tree variables
            fTreeVariablePID   = pThisV0->GetPdgCode(); //PDG Code
            fTreeVariablePtMC  = pThisV0->Pt(); //Perfect Pt
            
            //Only Interested if it's a Lambda, AntiLambda or K0s
            //Avoid the Junction Bug! PYTHIA has particles with Px=Py=Pz=E=0 occasionally,
            //having particle code 88 (unrecognized by PDG), for documentation purposes.
            //Even ROOT's TParticle::Y() is not prepared to deal with that exception!
            //Note that TParticle::Pt() is immune (that would just return 0)...
            //Though granted that that should be extremely rare in this precise condition...
            if( TMath::Abs(fTreeVariablePID) == 3122 || fTreeVariablePID==310 ) {
                fTreeVariableRapMC = pThisV0->Y(); //Perfect Y
            }
            if( lMCstack->IsPhysicalPrimary       (lblMotherPosV0Dghter) ) fTreeVariablePrimaryStatus = 1; //Is Primary!
            if( lMCstack->IsSecondaryFromWeakDecay(lblMotherPosV0Dghter) ) fTreeVariablePrimaryStatus = 2; //Weak Decay!
            if( lMCstack->IsSecondaryFromMaterial (lblMotherPosV0Dghter) ) fTreeVariablePrimaryStatus = 3; //Material Int!
            
            //Now we try to acquire the V0 parent particle, if possible
            Int_t lblThisV0Parent = pThisV0->GetFirstMother();
            if ( lblThisV0Parent > -1 ) { //if it has a parent, get it and store specs
                TParticle* pThisV0Parent = lMCstack->Particle( lblThisV0Parent );
                fTreeVariablePIDMother   = pThisV0Parent->GetPdgCode(); //V0 Mother PDG
                fTreeVariablePtMother    = pThisV0Parent->Pt();         //V0 Mother Pt
                //NOTE: Fill only for charged xi
                if ( TMath::Abs(fTreeVariablePIDMother)==3312) fTreeVariableRapMother   = pThisV0Parent->Y();         //V0 Mother Pt
                //Primary Status for the V0 Mother particle
                if( lMCstack->IsPhysicalPrimary       (lblThisV0Parent) ) fTreeVariablePrimaryStatusMother = 1; //Is Primary!
                if( lMCstack->IsSecondaryFromWeakDecay(lblThisV0Parent) ) fTreeVariablePrimaryStatusMother = 2; //Weak Decay!
                if( lMCstack->IsSecondaryFromMaterial (lblThisV0Parent) ) fTreeVariablePrimaryStatusMother = 3; //Material Int!
            }
        }
        
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
                 (!fkPreselectDedx || (fkPreselectDedx&&TMath::Abs(fTreeVariableNSigmasPosProton) < 7.0 && TMath::Abs(fTreeVariableNSigmasNegPion) < 7.0) ) &&
                 (!fkPreselectPID || fTreeVariablePID == 3122 )
                )
                ||
                //Case 2: AntiLambda Selection
                (fTreeVariableInvMassAntiLambda < lUpperLimitLambda  && fTreeVariableInvMassAntiLambda > lLowerLimitLambda &&
                 (!fkPreselectDedx || (fkPreselectDedx&&TMath::Abs(fTreeVariableNSigmasNegProton) < 7.0 && TMath::Abs(fTreeVariableNSigmasPosPion) < 7.0) ) &&
                 (!fkPreselectPID || fTreeVariablePID == -3122 )
                )
                ||
                //Case 3: K0Short Selection
                (fTreeVariableInvMassK0s        < lUpperLimitK0Short && fTreeVariableInvMassK0s        > lLowerLimitK0Short &&
                 (!fkPreselectDedx || (fkPreselectDedx&&TMath::Abs(fTreeVariableNSigmasNegPion)   < 7.0 && TMath::Abs(fTreeVariableNSigmasPosPion) < 7.0) ) &&
                 (!fkPreselectPID || fTreeVariablePID == 310 )
                ) ) {
                //Pre-selection in case this is AA...
                    
                    //Random denial
                    Bool_t lKeepV0 = kTRUE;
                    if(fkDownScaleV0 && ( fRand->Uniform() > fDownScaleFactorV0 )) lKeepV0 = kFALSE;
                    
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
        TH3F *histoout         = 0x0;
        AliV0Result *lV0Result = 0x0;
        for(Int_t lcfg=0; lcfg<lNumberOfConfigurations; lcfg++){
            lV0Result = (AliV0Result*) fListV0->At(lcfg);
            histoout  = lV0Result->GetHistogram();
            
            Float_t lMass = 0;
            Float_t lRap  = 0;
            Float_t lPDGMass = -1;
            Float_t lNegdEdx = 100;
            Float_t lPosdEdx = 100;
            Int_t lPDGCode = 0;
            Int_t lPDGCodeXiMother = 0;
            
            if ( lV0Result->GetMassHypothesis() == AliV0Result::kK0Short     ){
                lMass    = fTreeVariableInvMassK0s;
                lRap     = fTreeVariableRapK0Short;
                lPDGMass = 0.497;
                lNegdEdx = fTreeVariableNSigmasNegPion;
                lPosdEdx = fTreeVariableNSigmasPosPion;
                lPDGCode = 310;
            }
            if ( lV0Result->GetMassHypothesis() == AliV0Result::kLambda      ){
                lMass = fTreeVariableInvMassLambda;
                lRap = fTreeVariableRapLambda;
                lPDGMass = 1.115683;
                lNegdEdx = fTreeVariableNSigmasNegPion;
                lPosdEdx = fTreeVariableNSigmasPosProton;
                lPDGCode = 3122;
                lPDGCodeXiMother = 3312;
            }
            if ( lV0Result->GetMassHypothesis() == AliV0Result::kAntiLambda  ){
                lMass = fTreeVariableInvMassAntiLambda;
                lRap = fTreeVariableRapLambda;
                lPDGMass = 1.115683;
                lNegdEdx = fTreeVariableNSigmasNegProton;
                lPosdEdx = fTreeVariableNSigmasPosPion;
                lPDGCode = -3122;
                lPDGCodeXiMother = -3312;
            }
            
            //Override rapidity for true rapidity if requested to do so
            if ( lV0Result -> GetCutMCUseMCProperties() ) {
                lRap = fTreeVariableRapMC;
            }
            
            if (
                //Check 1: Offline Vertexer
                lOnFlyStatus == 0 &&
                
                //Check 2: Basic Acceptance cuts
                TMath::Abs(fTreeVariableNegEta)<0.8 &&
                TMath::Abs(fTreeVariablePosEta)<0.8 &&
                TMath::Abs(lRap) < 0.5 &&
                
                //Check 3: Topological Variables
                fTreeVariableV0Radius > lV0Result->GetCutV0Radius() &&
                fTreeVariableDcaNegToPrimVertex > lV0Result->GetCutDCANegToPV() &&
                fTreeVariableDcaPosToPrimVertex > lV0Result->GetCutDCAPosToPV() &&
                fTreeVariableDcaV0Daughters < lV0Result->GetCutDCAV0Daughters() &&
                fTreeVariableV0CosineOfPointingAngle > lV0Result->GetCutV0CosPA() &&
                fTreeVariableDistOverTotMom*lPDGMass < lV0Result->GetCutProperLifetime() &&
                fTreeVariableLeastNbrCrossedRows > lV0Result->GetCutLeastNumberOfCrossedRows() &&
                fTreeVariableLeastRatioCrossedRowsOverFindable > lV0Result->GetCutLeastNumberOfCrossedRowsOverFindable() &&
                
                // - MC specific: either don't associate (if not requested) or associate
                // - FIXME: the 2D efficiency matrix needs to be filled according to Xi momentum as well.
                // - This problem still has to be addressed in this case.
                // - Different Options:
                //    - Use TTree (will eventually become a problem as well, though with association might be OK)
                //    - Create another out histogram in AliV0Result to be able to do the momentum correlations?
                ( ! (lV0Result->GetCutMCPhysicalPrimary())    || fTreeVariablePrimaryStatus == 1 ) &&
                ( ! (lV0Result->GetCutMCLambdaFromPrimaryXi())|| (fTreeVariablePrimaryStatusMother == 1 && fTreeVariablePIDMother == lPDGCodeXiMother) ) &&
                ( ! (lV0Result->GetCutMCPDGCodeAssociation()) || fTreeVariablePID == lPDGCode     ) &&
                
                //Check 4: TPC dEdx selections
                TMath::Abs(lNegdEdx)<lV0Result->GetCutTPCdEdx() &&
                TMath::Abs(lPosdEdx)<lV0Result->GetCutTPCdEdx() &&
            
                //Check 5: Armenteros-Podolanski space cut (for K0Short analysis)
                ( ( lV0Result->GetCutArmenteros() == kFALSE || lV0Result->GetMassHypothesis() != AliV0Result::kK0Short ) || ( fTreeVariablePtArmV0*5>TMath::Abs(fTreeVariableAlphaV0) ) )
                )
            {
                //This satisfies all my conditionals! Fill histogram
                if( !lV0Result -> GetCutMCUseMCProperties() ){
                    histoout -> Fill ( fCentrality, fTreeVariablePt, lMass );
                }else{
                    histoout -> Fill ( fCentrality, fTreeVariablePtMC, lMass );
                }
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
        lESDevent->ResetCascades();
        lESDevent->ResetV0s();
        
        //Decide between regular and light vertexer (default: light)
        if ( ! fkUseLightVertexer ){
            AliV0vertexer lV0vtxer;
            AliCascadeVertexer lCascVtxer;
            
            lV0vtxer.SetDefaultCuts(fV0VertexerSels);
            lCascVtxer.SetDefaultCuts(fCascadeVertexerSels);

            lV0vtxer.SetCuts(fV0VertexerSels);
            lCascVtxer.SetCuts(fCascadeVertexerSels);
            
            lV0vtxer.Tracks2V0vertices(lESDevent);
            lCascVtxer.V0sTracks2CascadeVertices(lESDevent);
        } else {
            AliLightV0vertexer lV0vtxer;
            AliLightCascadeVertexer lCascVtxer;
            
            lV0vtxer.SetDefaultCuts(fV0VertexerSels);
            lCascVtxer.SetDefaultCuts(fCascadeVertexerSels);
            
            lV0vtxer.SetCuts(fV0VertexerSels);
            lCascVtxer.SetCuts(fCascadeVertexerSels);
            
            lV0vtxer.Tracks2V0vertices(lESDevent);
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

        Short_t  lChargeXi = -2;
        //Double_t lV0toXiCosineOfPointingAngle = 0. ;

        Double_t lRapXi   = -20.0, lRapOmega = -20.0, lRapMC = -20.0; //  lEta = -20.0, lTheta = 360., lPhi = 720. ;
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
        //lChi2Xi 			    = xi->GetChi2Xi();
        lDcaXiDaughters 	= xi->GetDcaXiDaughters();
        lXiCosineOfPointingAngle 	            = xi->GetCascadeCosineOfPointingAngle( lBestPrimaryVtxPos[0],
                lBestPrimaryVtxPos[1],
                lBestPrimaryVtxPos[2] );
        // Take care : the best available vertex should be used (like in AliCascadeVertexer)

        xi->GetXYZcascade( lPosXi[0],  lPosXi[1], lPosXi[2] );
        lXiRadius			= TMath::Sqrt( lPosXi[0]*lPosXi[0]  +  lPosXi[1]*lPosXi[1] );

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

        if (!pTrackXi || !nTrackXi || !bachTrackXi ) {
            AliWarning("ERROR: Could not retrieve one of the 3 ESD daughter tracks of the cascade ...");
            continue;
        }

        fTreeCascVarPosEta = pTrackXi->Eta();
        fTreeCascVarNegEta = nTrackXi->Eta();
        fTreeCascVarBachEta = bachTrackXi->Eta();

        Double_t lBMom[3], lNMom[3], lPMom[3];
        xi->GetBPxPyPz( lBMom[0], lBMom[1], lBMom[2] );
        xi->GetPPxPyPz( lPMom[0], lPMom[1], lPMom[2] );
        xi->GetNPxPyPz( lNMom[0], lNMom[1], lNMom[2] );

        //fTreeCascVarBachTransMom = TMath::Sqrt( lBMom[0]*lBMom[0] + lBMom[1]*lBMom[1] );
        //fTreeCascVarPosTransMom  = TMath::Sqrt( lPMom[0]*lPMom[0] + lPMom[1]*lPMom[1] );
        //fTreeCascVarNegTransMom  = TMath::Sqrt( lNMom[0]*lNMom[0] + lNMom[1]*lNMom[1] );
        
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

        // 2 - Poor quality related to TPC clusters: lowest cut of 70 clusters
        if(lPosTPCClusters  < 70) {
            AliDebug(1, "Pb / V0 Pos. track has less than 70 TPC clusters ... continue!");
            continue;
        }
        if(lNegTPCClusters  < 70) {
            AliDebug(1, "Pb / V0 Neg. track has less than 70 TPC clusters ... continue!");
            continue;
        }
        if(lBachTPCClusters < 70) {
            AliDebug(1, "Pb / Bach.   track has less than 70 TPC clusters ... continue!");
            continue;
        }
        Int_t leastnumberofclusters = 1000;
        if( lPosTPCClusters < leastnumberofclusters ) leastnumberofclusters = lPosTPCClusters;
        if( lNegTPCClusters < leastnumberofclusters ) leastnumberofclusters = lNegTPCClusters;
        if( lBachTPCClusters < leastnumberofclusters ) leastnumberofclusters = lBachTPCClusters;

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
        
        
        //------------------------------------------------
        // Associate Cascade Candidates to Monte Carlo!
        //------------------------------------------------
        
        //Warning: Not using Continues... Need to fill tree later!
        
        Double_t lXiTransvMomMC= -100. ;
        Int_t lPDGCodeCascade = 0;
        Int_t lPID_BachMother = 0;
        Int_t lPID_NegMother = 0;
        Int_t lPID_PosMother = 0;
        fTreeCascVarIsPhysicalPrimary = 0; // 0: not defined, any candidate may have this
        //fTreeCascVarPosTransvMomentumMC = -1;
        //fTreeCascVarNegTransvMomentumMC = -1;
        
        if(fDebug > 5)
            cout 	<< "MC EventNumber : " << lMCevent->Header()->GetEvent()
            << " / MC event Number in Run : " << lMCevent->Header()->GetEventNrInRun() << endl;
        
        //----------------------------------------
        // Regular MC ASSOCIATION STARTS HERE
        //----------------------------------------
        
        Int_t lblPosV0Dghter = (Int_t) TMath::Abs( pTrackXi->GetLabel() );
        // Abs value = needed ! question of quality track association ...
        Int_t lblNegV0Dghter = (Int_t) TMath::Abs( nTrackXi->GetLabel() );
        Int_t lblBach        = (Int_t) TMath::Abs( bachTrackXi->GetLabel() );
        
        TParticle* mcPosV0Dghter = lMCstack->Particle( lblPosV0Dghter );
        TParticle* mcNegV0Dghter = lMCstack->Particle( lblNegV0Dghter );
        TParticle* mcBach        = lMCstack->Particle( lblBach );
        
        //fTreeCascVarPosTransvMomentumMC = mcPosV0Dghter->Pt();
        //fTreeCascVarNegTransvMomentumMC = mcNegV0Dghter->Pt();
        
        //fTreeCascVarPIDPositive = mcPosV0Dghter -> GetPdgCode();
        //fTreeCascVarPIDNegative = mcNegV0Dghter -> GetPdgCode();
        //fTreeCascVarPIDBachelor = mcBach->GetPdgCode();
        
        // - Step 4.2 : level of the Xi daughters
        
        Int_t lblMotherPosV0Dghter = mcPosV0Dghter->GetFirstMother() ;
        Int_t lblMotherNegV0Dghter = mcNegV0Dghter->GetFirstMother();
        
        //Rather uncivilized: Open brackets for each 'continue'
        if(! (lblMotherPosV0Dghter != lblMotherNegV0Dghter) ) { // same mother
            if(! (lblMotherPosV0Dghter < 0) ) { // mother != primary (!= -1)
                if(! (lblMotherNegV0Dghter < 0) ) {
                    
                    // mothers = Lambda candidate ... a priori
                    
                    TParticle* mcMotherPosV0Dghter = lMCstack->Particle( lblMotherPosV0Dghter );
                    TParticle* mcMotherNegV0Dghter = lMCstack->Particle( lblMotherNegV0Dghter );
                    
                    // - Step 4.3 : level of Xi candidate
                    
                    Int_t lblGdMotherPosV0Dghter =   mcMotherPosV0Dghter->GetFirstMother() ;
                    Int_t lblGdMotherNegV0Dghter =   mcMotherNegV0Dghter->GetFirstMother() ;
                    
                    if(! (lblGdMotherPosV0Dghter != lblGdMotherNegV0Dghter) ) {
                        if(! (lblGdMotherPosV0Dghter < 0) ) { // primary lambda ...
                            if(! (lblGdMotherNegV0Dghter < 0) ) { // primary lambda ...
                                
                                // Gd mothers = Xi candidate ... a priori
                                
                                TParticle* mcGdMotherPosV0Dghter = lMCstack->Particle( lblGdMotherPosV0Dghter );
                                TParticle* mcGdMotherNegV0Dghter = lMCstack->Particle( lblGdMotherNegV0Dghter );
                                
                                Int_t lblMotherBach = (Int_t) TMath::Abs( mcBach->GetFirstMother()  );
                                
                                //		if( lblMotherBach != lblGdMotherPosV0Dghter ) continue; //same mother for bach and V0 daughters
                                if(!(lblMotherBach != lblGdMotherPosV0Dghter)) { //same mother for bach and V0 daughters
                                    
                                    TParticle* mcMotherBach = lMCstack->Particle( lblMotherBach );
                                    
                                    lPID_BachMother = mcMotherBach->GetPdgCode();
                                    lPID_NegMother = mcGdMotherPosV0Dghter->GetPdgCode();
                                    lPID_PosMother = mcGdMotherNegV0Dghter->GetPdgCode();
                                    
                                    if(lPID_BachMother==lPID_NegMother && lPID_BachMother==lPID_PosMother) {
                                        lPDGCodeCascade = lPID_BachMother;
                                        lXiTransvMomMC = mcMotherBach->Pt();
                                        if( lMCstack->IsPhysicalPrimary       (lblMotherBach) ) fTreeCascVarIsPhysicalPrimary = 1; //Is Primary!
                                        if( lMCstack->IsSecondaryFromWeakDecay(lblMotherBach) ) fTreeCascVarIsPhysicalPrimary = 2; //Weak Decay!
                                        if( lMCstack->IsSecondaryFromMaterial (lblMotherBach) ) fTreeCascVarIsPhysicalPrimary = 3; //From Material!
                                        if ( (mcMotherBach->Energy() + mcMotherBach->Pz()) / (mcMotherBach->Energy() - mcMotherBach->Pz() +1.e-13) !=0 ) {
                                            lRapMC = 0.5*TMath::Log( (mcMotherBach->Energy() + mcMotherBach->Pz()) / (mcMotherBach->Energy() - mcMotherBach->Pz() +1.e-13) );
                                        }
                                    }
                                    
                                }
                            }
                        }
                    }
                }
            }
        } //Ends all conditionals above...
        
        //----------------------------------------
        // Regular MC ASSOCIATION ENDS HERE
        //----------------------------------------

        //----------------------------------------
        // Swapped MC ASSOCIATION STARTS HERE
        // WARNING: THIS IS EXPERIMENTAL!
        //----------------------------------------
        
        // Abs value = needed ! question of quality track association ...
        Int_t lblPosV0DghterSwapped = -1;
        Int_t lblNegV0DghterSwapped = -1;
        Int_t lblBachSwapped = -1;
        
        //Here's where it all happens: please swap stuff around
        //...and "Keep Calm and Carry On"
        if (lChargeXi < 0){
            lblPosV0DghterSwapped = lblPosV0Dghter;
            lblNegV0DghterSwapped = lblBach;
            lblBachSwapped        = lblNegV0Dghter;
        }else{
            lblPosV0DghterSwapped = lblBach;
            lblNegV0DghterSwapped = lblNegV0Dghter;
            lblBachSwapped        = lblPosV0Dghter;
        }
        
        TParticle* mcPosV0DghterSwapped = lMCstack->Particle( lblPosV0DghterSwapped );
        TParticle* mcNegV0DghterSwapped = lMCstack->Particle( lblNegV0DghterSwapped );
        TParticle* mcBachSwapped        = lMCstack->Particle( lblBachSwapped );
        
        //Needed, just for checks
        Int_t lPID_BachMotherSwapped = 0;
        Int_t lPID_NegMotherSwapped = 0;
        Int_t lPID_PosMotherSwapped = 0;
        Int_t lPDGCodeCascadeSwapped = 0;
        
        // - Step 4.2 : level of the Xi daughters
        
        Int_t lblMotherPosV0DghterSwapped = mcPosV0DghterSwapped->GetFirstMother() ;
        Int_t lblMotherNegV0DghterSwapped = mcNegV0DghterSwapped->GetFirstMother();
        
        //Rather uncivilized: Open brackets for each 'continue'
        if(! (lblMotherPosV0DghterSwapped != lblMotherNegV0DghterSwapped) ) { // same mother
            if(! (lblMotherPosV0DghterSwapped < 0) ) { // mother != primary (!= -1)
                if(! (lblMotherNegV0DghterSwapped < 0) ) {
                    
                    // mothers = Lambda candidate ... a priori
                    
                    TParticle* mcMotherPosV0DghterSwapped = lMCstack->Particle( lblMotherPosV0DghterSwapped );
                    TParticle* mcMotherNegV0DghterSwapped = lMCstack->Particle( lblMotherNegV0DghterSwapped );
                    
                    // - Step 4.3 : level of Xi candidate
                    
                    Int_t lblGdMotherPosV0DghterSwapped =   mcMotherPosV0DghterSwapped->GetFirstMother() ;
                    Int_t lblGdMotherNegV0DghterSwapped =   mcMotherNegV0DghterSwapped->GetFirstMother() ;
                    
                    if(! (lblGdMotherPosV0DghterSwapped != lblGdMotherNegV0DghterSwapped) ) {
                        if(! (lblGdMotherPosV0DghterSwapped < 0) ) { // primary lambda ...
                            if(! (lblGdMotherNegV0DghterSwapped < 0) ) { // primary lambda ...
                                
                                // Gd mothers = Xi candidate ... a priori
                                
                                TParticle* mcGdMotherPosV0DghterSwapped = lMCstack->Particle( lblGdMotherPosV0DghterSwapped );
                                TParticle* mcGdMotherNegV0DghterSwapped = lMCstack->Particle( lblGdMotherNegV0DghterSwapped );
                                
                                Int_t lblMotherBachSwapped = (Int_t) TMath::Abs( mcBachSwapped->GetFirstMother()  );
                                
                                //		if( lblMotherBach != lblGdMotherPosV0Dghter ) continue; //same mother for bach and V0 daughters
                                if(!(lblMotherBachSwapped != lblGdMotherPosV0DghterSwapped)) { //same mother for bach and V0 daughters
                                    
                                    TParticle* mcMotherBachSwapped = lMCstack->Particle( lblMotherBachSwapped );
                                    
                                    lPID_BachMotherSwapped = mcMotherBachSwapped->GetPdgCode();
                                    lPID_NegMotherSwapped = mcGdMotherPosV0DghterSwapped->GetPdgCode();
                                    lPID_PosMotherSwapped = mcGdMotherNegV0DghterSwapped->GetPdgCode();
                                    
                                    if(lPID_BachMotherSwapped==lPID_NegMotherSwapped && lPID_BachMotherSwapped==lPID_PosMotherSwapped) {
                                        lPDGCodeCascadeSwapped = lPID_BachMotherSwapped;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } //Ends all conditionals above...
        
        //----------------------------------------
        // Swapped MC ASSOCIATION ENDS HERE
        //----------------------------------------

        
        //------------------------------------------------
        // Set Variables for adding to tree
        //------------------------------------------------
        
        fTreeCascVarCharge	= lChargeXi;
        if (lChargeXi < 0 ){
            fTreeCascVarMassAsXi    = lInvMassXiMinus;
            fTreeCascVarMassAsOmega = lInvMassOmegaMinus;
        }
        if (lChargeXi > 0 ){
            fTreeCascVarMassAsXi    = lInvMassXiPlus;
            fTreeCascVarMassAsOmega = lInvMassOmegaPlus;
        }
        fTreeCascVarPID = lPDGCodeCascade; 
        fTreeCascVarSwappedPID = lPDGCodeCascadeSwapped;
        fTreeCascVarPt = lXiTransvMom;
        fTreeCascVarPtMC = lXiTransvMomMC;
        fTreeCascVarRapXi = lRapXi ;
        fTreeCascVarRapOmega = lRapOmega ;
        fTreeCascVarRapMC = lRapMC ;
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

        //Copy Multiplicity information
        fTreeCascVarCentrality = fCentrality;

        fTreeCascVarDistOverTotMom = TMath::Sqrt(
                                         TMath::Power( lPosXi[0] - lBestPrimaryVtxPos[0] , 2) +
                                         TMath::Power( lPosXi[1] - lBestPrimaryVtxPos[1] , 2) +
                                         TMath::Power( lPosXi[2] - lBestPrimaryVtxPos[2] , 2)
                                     );
        fTreeCascVarDistOverTotMom /= (lXiTotMom+1e-13);

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

        if( fkSaveCascadeTree && lKeepCascade && (
                                                  //Xi Conditionals
                                                  (fTreeCascVarMassAsXi<1.32+0.075&&fTreeCascVarMassAsXi>1.32-0.075 &&
                                                   (!fkPreselectPID || TMath::Abs(fTreeCascVarPID) == 3312 )) ||
                                                  //Omega Conditionals
                                                  (fTreeCascVarMassAsOmega<1.68+0.075&&fTreeCascVarMassAsOmega>1.68-0.075 &&
                                                   (!fkPreselectPID || TMath::Abs(fTreeCascVarPID) == 3334 ))
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
        TH3F *histoout         = 0x0;
        AliCascadeResult *lCascadeResult = 0x0;
        for(Int_t lcfg=0; lcfg<lNumberOfConfigurationsCascade; lcfg++){
            lCascadeResult = (AliCascadeResult*) fListCascade->At(lcfg);
            histoout  = lCascadeResult->GetHistogram();
            
            Float_t lMass = 0;
            Float_t lRap  = 0;
            Float_t lPDGMass = -1;
            Float_t lNegdEdx = 100;
            Float_t lPosdEdx = 100;
            Float_t lBachdEdx = 100;
            Short_t  lCharge = -2;
            Int_t lPDGCode = 0;
            
            if ( lCascadeResult->GetMassHypothesis() == AliCascadeResult::kXiMinus     ){
                lCharge  = -1;
                lMass    = fTreeCascVarMassAsXi;
                lRap     = fTreeCascVarRapXi;
                lPDGMass = 1.32171;
                lNegdEdx = fTreeCascVarNegNSigmaPion;
                lPosdEdx = fTreeCascVarPosNSigmaProton;
                lBachdEdx= fTreeCascVarBachNSigmaPion;
                lPDGCode = 3312;
            }
            if ( lCascadeResult->GetMassHypothesis() == AliCascadeResult::kXiPlus      ){
                lCharge  = +1;
                lMass    = fTreeCascVarMassAsXi;
                lRap     = fTreeCascVarRapXi;
                lPDGMass = 1.32171;
                lNegdEdx = fTreeCascVarNegNSigmaProton;
                lPosdEdx = fTreeCascVarPosNSigmaPion;
                lBachdEdx= fTreeCascVarBachNSigmaPion;
                lPDGCode = -3312;
            }
            if ( lCascadeResult->GetMassHypothesis() == AliCascadeResult::kOmegaMinus     ){
                lCharge  = -1;
                lMass    = fTreeCascVarMassAsOmega;
                lRap     = fTreeCascVarRapOmega;
                lPDGMass = 1.67245;
                lNegdEdx = fTreeCascVarNegNSigmaPion;
                lPosdEdx = fTreeCascVarPosNSigmaProton;
                lBachdEdx= fTreeCascVarBachNSigmaKaon;
                lPDGCode = 3334;
            }
            if ( lCascadeResult->GetMassHypothesis() == AliCascadeResult::kOmegaPlus      ){
                lCharge  = +1;
                lMass    = fTreeCascVarMassAsOmega;
                lRap     = fTreeCascVarRapOmega;
                lPDGMass = 1.67245;
                lNegdEdx = fTreeCascVarNegNSigmaProton;
                lPosdEdx = fTreeCascVarPosNSigmaPion;
                lBachdEdx= fTreeCascVarBachNSigmaKaon;
                lPDGCode = -3334;
            }
            
            //Override rapidity for true rapidity if requested to do so
            if ( lCascadeResult -> GetCutMCUseMCProperties() ) {
                lRap = fTreeVariableRapMC;
            }
            
            
            if (
                //Check 1: Charge consistent with expectations
                fTreeCascVarCharge == lCharge &&
                
                //Check 2: Basic Acceptance cuts
                TMath::Abs(fTreeCascVarPosEta)<0.8 &&
                TMath::Abs(fTreeCascVarNegEta)<0.8 &&
                TMath::Abs(fTreeCascVarBachEta)<0.8 &&
                TMath::Abs(lRap) < 0.5 &&
                
                //Check 3: Topological Variables
                // - V0 Selections
                fTreeCascVarDCANegToPrimVtx > lCascadeResult->GetCutDCANegToPV() &&
                fTreeCascVarDCAPosToPrimVtx > lCascadeResult->GetCutDCAPosToPV() &&
                fTreeCascVarDCAV0Daughters < lCascadeResult->GetCutDCAV0Daughters() &&
                fTreeCascVarV0CosPointingAngle > lCascadeResult->GetCutV0CosPA() &&
                fTreeCascVarV0Radius > lCascadeResult->GetCutV0Radius() &&
                // - Cascade Selections
                fTreeCascVarDCAV0ToPrimVtx > lCascadeResult->GetCutDCAV0ToPV() &&
                TMath::Abs(fTreeCascVarV0Mass-1.116) < lCascadeResult->GetCutV0Mass() &&
                fTreeCascVarDCABachToPrimVtx > lCascadeResult->GetCutDCABachToPV() &&
                fTreeCascVarDCACascDaughters < lCascadeResult->GetCutDCACascDaughters() &&
                fTreeCascVarCascCosPointingAngle > lCascadeResult->GetCutCascCosPA() &&
                fTreeCascVarCascRadius > lCascadeResult->GetCutCascRadius() &&
                
                // - Miscellaneous
                fTreeCascVarDistOverTotMom*lPDGMass < lCascadeResult->GetCutProperLifetime() &&
                fTreeCascVarLeastNbrClusters > lCascadeResult->GetCutLeastNumberOfClusters() &&
                
                // - MC specific: either don't associate (if not requested) or associate
                ( ! (lCascadeResult->GetCutMCPhysicalPrimary())    || fTreeCascVarIsPhysicalPrimary == kTRUE ) &&
                ( ! (lCascadeResult->GetCutMCPDGCodeAssociation()) || fTreeCascVarPID == lPDGCode            ) &&
                
                //Check 4: TPC dEdx selections
                TMath::Abs(lNegdEdx )<lCascadeResult->GetCutTPCdEdx() &&
                TMath::Abs(lPosdEdx )<lCascadeResult->GetCutTPCdEdx() &&
                TMath::Abs(lBachdEdx)<lCascadeResult->GetCutTPCdEdx() &&
                
                //Check 5: Xi rejection for Omega analysis
                ( ( lCascadeResult->GetMassHypothesis() != AliCascadeResult::kOmegaMinus || lCascadeResult->GetMassHypothesis() != AliCascadeResult::kOmegaPlus  ) || ( TMath::Abs( fTreeCascVarMassAsXi - 1.32171 ) > lCascadeResult->GetCutXiRejection() ) )
                ){
                
                //This satisfies all my conditionals! Fill histogram
                histoout -> Fill ( fCentrality, fTreeCascVarPt, lMass );
                
                //This satisfies all my conditionals! Fill histogram
                if( !lCascadeResult -> GetCutMCUseMCProperties() ){
                    histoout -> Fill ( fCentrality, fTreeCascVarPt, lMass );
                }else{
                    histoout -> Fill ( fCentrality, fTreeCascVarPtMC, lMass );
                }
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
void AliAnalysisTaskStrangenessVsMultiplicityMCRun2::Terminate(Option_t *)
{
    // Draw result to the screen
    // Called once at the end of the query

    TList *cRetrievedList = 0x0;
    cRetrievedList = (TList*)GetOutputData(1);
    if(!cRetrievedList) {
        Printf("ERROR - AliAnalysisTaskStrangenessVsMultiplicityMCRun2 : ouput data container list not available\n");
        return;
    }

    fHistEventCounter = dynamic_cast<TH1D*> (  cRetrievedList->FindObject("fHistEventCounter")  );
    if (!fHistEventCounter) {
        Printf("ERROR - AliAnalysisTaskStrangenessVsMultiplicityMCRun2 : fHistEventCounter not available");
        return;
    }

    TCanvas *canCheck = new TCanvas("AliAnalysisTaskStrangenessVsMultiplicityMCRun2","V0 Multiplicity",10,10,510,510);
    canCheck->cd(1)->SetLogy();

    fHistEventCounter->SetMarkerStyle(22);
    fHistEventCounter->DrawCopy("E");
}

//________________________________________________________________________
Double_t AliAnalysisTaskStrangenessVsMultiplicityMCRun2::MyRapidity(Double_t rE, Double_t rPz) const
{
    // Local calculation for rapidity
    Double_t ReturnValue = -100;
    if( (rE-rPz+1.e-13) != 0 && (rE+rPz) != 0 ) {
        ReturnValue =  0.5*TMath::Log((rE+rPz)/(rE-rPz+1.e-13));
    }
    return ReturnValue;
}

//________________________________________________________________________
void AliAnalysisTaskStrangenessVsMultiplicityMCRun2::AddConfiguration( AliV0Result *lV0Result )
{
    if (!fListV0){
        Printf("fListV0 does not exist. Creating...");
        fListV0 = new TList();
        fListV0->SetOwner();
        
    }
    fListV0->Add(lV0Result);
}

//________________________________________________________________________
void AliAnalysisTaskStrangenessVsMultiplicityMCRun2::AddConfiguration( AliCascadeResult *lCascadeResult )
{
    if (!fListCascade){
        Printf("fListCascade does not exist. Creating...");
        fListCascade = new TList();
        fListCascade->SetOwner();
        
    }
    fListCascade->Add(lCascadeResult);
}
