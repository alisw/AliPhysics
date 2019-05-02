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

#ifndef AliAnalysisTaskStrEffStudy_H
#define AliAnalysisTaskStrEffStudy_H

class TList;
class TH1D;
class TH3D;
class TVector3;
class TRandom3;
class TTree;

class AliESDtrackCuts;
class AliAnalysisUtils;
class AliESDEvent;
class AliESDtrack;
class AliPIDResponse;
class AliV0Result;
class AliCascadeResult;
class AliExternalTrackParam;

#include "TString.h"
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"

class AliAnalysisTaskStrEffStudy : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskStrEffStudy();
    AliAnalysisTaskStrEffStudy(Bool_t lSaveEventTree, Bool_t lSaveV0Tree, Bool_t lSaveCascadeTree, Bool_t lSaveHyperTriton, const char *name, TString lExtraOptions = "");
    virtual ~AliAnalysisTaskStrEffStudy();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    Double_t MyRapidity(Double_t rE, Double_t rPz) const;

    //Fix on-the-fly v0s
    void CheckChargeV0(AliESDv0 *v0);

    void SetSaveV0s                (Bool_t lSaveV0s        = kTRUE ) {
        fkSaveV0Tree        = lSaveV0s;
    }
    void SetSaveCascades           (Bool_t lSaveCascades   = kTRUE ) {
        fkSaveCascadeTree   = lSaveCascades;
    }
    void SetSaveHyperTriton3Body (Bool_t lSaveHyp3Body   = kTRUE ) {
        fkSaveHyperTriton3BodyTree   = lSaveHyp3Body;
    }
    void SetPreselectDedx (Bool_t lPreselectDedx= kTRUE ) {
        fkPreselectDedx   = lPreselectDedx;
    }
    void SetPreselectPID (Bool_t lPreselectPID = kTRUE ) {
        fkPreselectPID   = lPreselectPID;
    }
    void SetUseOnTheFlyV0Cascading( Bool_t lUseOnTheFlyV0Cascading = kTRUE ){
        //Highly experimental, use with care!
        fkUseOnTheFlyV0Cascading = lUseOnTheFlyV0Cascading;
    }
    void SetDoImprovedCascadeVertexFinding( Bool_t lOpt = kTRUE ){
        //Highly experimental, use with care!
        fkDoImprovedCascadeVertexFinding = lOpt;
    }
    void SetDoImprovedDCAV0DauPropagation( Bool_t lOpt = kTRUE ){
        //Highly experimental, use with care!
        fkDoImprovedDCAV0DauPropagation = lOpt;
    }
    void SetDoImprovedDCACascDauPropagation( Bool_t lOpt = kTRUE ){
        //Highly experimental, use with care!
        fkDoImprovedDCACascDauPropagation = lOpt;
    }
    void SetIfImprovedPerformInitialLinearPropag( Bool_t lOpt = kTRUE ){
        //Highly experimental, use with care!
        fkIfImprovedPerformInitialLinearPropag = lOpt;
    }
    void SetIfImprovedExtraPrecisionFactor( Double_t lOpt ){
        //Highly experimental, use with care!
        fkIfImprovedExtraPrecisionFactor = lOpt;
    }

//---------------------------------------------------------------------------------------
    //Task Configuration: trigger selection
    void SetSelectedTriggerClass(AliVEvent::EOfflineTriggerTypes trigType) { fTrigType = trigType;}
//---------------------------------------------------------------------------------------
    //Task Configuration: Meant to enable quick re-execution of vertexer if needed
    void SetRunVertexers ( Bool_t lRunVertexers = kTRUE) {
        fkRunVertexers = lRunVertexers;
    }
    void SetUseLightVertexers ( Bool_t lUseLightVertexers = kTRUE) {
        fkUseLightVertexer = lUseLightVertexers;
    }
    void SetDoV0Refit ( Bool_t lDoV0Refit = kTRUE) {
        fkDoV0Refit = lDoV0Refit;
    }
    void SetExtraCleanup ( Bool_t lExtraCleanup = kTRUE) {
        fkExtraCleanup = lExtraCleanup;
    }
    void SetSaveGoodTracks ( Bool_t lOpt = kTRUE) {
        fkSaveGoodTracks = lOpt;
    }
    void SetSandboxV0 ( Bool_t lOpt = kTRUE) {
        fkSandboxV0 = lOpt;
    }
    void SetSandboxCascade ( Bool_t lOpt = kTRUE) {
        fkSandboxCascade = lOpt;
    }

//---------------------------------------------------------------------------------------
    void SetUseExtraEvSels ( Bool_t lUseExtraEvSels = kTRUE) {
        fkDoExtraEvSels = lUseExtraEvSels;
    }
//---------------------------------------------------------------------------------------
    //Task Configuration: Skip Event Selections after trigger (VZERO test)
    void SetDownScaleV0 ( Bool_t lOpt = kTRUE, Float_t lVal = 0.001) {
        fkDownScaleV0 = lOpt;
        fDownScaleFactorV0 = lVal;
    }
    void SetDownScaleCascade ( Bool_t lOpt = kTRUE, Float_t lVal = 0.001 ) {
        fkDownScaleCascade = lOpt;
        fDownScaleFactorCascade = lVal;
    }
    void SetDownScaleHyperTriton3Body ( Bool_t lOpt = kTRUE, Float_t lVal = 0.001 ) {
        fkDownScaleHyperTriton3Body = lOpt;
        fDownScaleFactorHyperTriton3Body = lVal;
    }
//---------------------------------------------------------------------------------------
//Setters for the V0 Vertexer Parameters
    void SetV0VertexerMaxChisquare   ( Double_t lParameter ) {
        fV0VertexerSels[0] = lParameter;
    }
    void SetV0VertexerDCAFirstToPV   ( Double_t lParameter ) {
        fV0VertexerSels[1] = lParameter;
    }
    void SetV0VertexerDCASecondtoPV  ( Double_t lParameter ) {
        fV0VertexerSels[2] = lParameter;
    }
    void SetV0VertexerDCAV0Daughters ( Double_t lParameter ) {
        fV0VertexerSels[3] = lParameter;
    }
    void SetV0VertexerCosinePA       ( Double_t lParameter ) {
        fV0VertexerSels[4] = lParameter;
    }
    void SetV0VertexerMinRadius      ( Double_t lParameter ) {
        fV0VertexerSels[5] = lParameter;
    }
    void SetV0VertexerMaxRadius      ( Double_t lParameter ) {
        fV0VertexerSels[6] = lParameter;
    }
//---------------------------------------------------------------------------------------
//Setters for the Cascade Vertexer Parameters
    void SetCascVertexerMaxChisquare         ( Double_t lParameter ) {
        fCascadeVertexerSels[0] = lParameter;
    }
    void SetCascVertexerMinV0ImpactParameter ( Double_t lParameter ) {
        fCascadeVertexerSels[1] = lParameter;
    }
    void SetCascVertexerV0MassWindow         ( Double_t lParameter ) {
        fCascadeVertexerSels[2] = lParameter;
    }
    void SetCascVertexerDCABachToPV          ( Double_t lParameter ) {
        fCascadeVertexerSels[3] = lParameter;
    }
    void SetCascVertexerDCACascadeDaughters  ( Double_t lParameter ) {
        fCascadeVertexerSels[4] = lParameter;
    }
    void SetCascVertexerCascadeCosinePA      ( Double_t lParameter ) {
        fCascadeVertexerSels[5] = lParameter;
    }
    void SetCascVertexerCascadeMinRadius     ( Double_t lParameter ) {
        fCascadeVertexerSels[6] = lParameter;
    }
    void SetCascVertexerCascadeMaxRadius     ( Double_t lParameter ) {
        fCascadeVertexerSels[7] = lParameter;
    }
    //---------------------------------------------------------------------------------------
    void SetMinPt     ( Float_t lMinPt ) {
        fMinPtToSave = lMinPt;
    }
    void SetMaxPt     ( Float_t lMaxPt ) {
        fMaxPtToSave = lMaxPt;
    }
    void SetLambdaWindowParameters     ( Double_t *fMeanPars, Double_t *fSigmaPars ) {
        for(Int_t ipar=0; ipar<5; ipar++) fLambdaMassMean[ipar]  = fMeanPars[ipar];
        for(Int_t ipar=0; ipar<4; ipar++) fLambdaMassSigma[ipar] = fSigmaPars[ipar];
    }
    void SetLambdaWindowParametersStandard (){
        fLambdaMassMean[0] =  1.15768e+00;
        fLambdaMassMean[1] = -4.15945e-02;
        fLambdaMassMean[2] = -7.14294e-04;
        fLambdaMassMean[3] = -1.62793e-02;
        fLambdaMassMean[4] = -7.84067e+00;
        fLambdaMassSigma[0] = 1.30345e-03;
        fLambdaMassSigma[1] = 2.89679e-04;
        fLambdaMassSigma[2] = 1.52661e-03;
        fLambdaMassSigma[3] =-2.58251e+00;
    }
    //---------------------------------------------------------------------------------------
    void SetPrecisionCutoffCascadeDCA     ( Double_t lPrecision ) {
        fPrecisionCutoffCascadeDCA = lPrecision;
    }
    void SetMaxIterationsCascadeDCA     ( Int_t lMaxIter ) {
        fMaxIterationsCascadeDCA = lMaxIter;
    }
    //---------------------------------------------------------------------------------------
    //Superlight mode: add another configuration, please
    void AddConfiguration( AliV0Result      *lV0Result      );
    void AddConfiguration( AliCascadeResult *lCascadeResult );
    //---------------------------------------------------------------------------------------
    //Functions for analysis Bookkeepinp
    // 1- Configure standard vertexing
    void SetupStandardVertexing();
    void SetupLooseVertexing();
    // 2- Standard Topological Selection QA Sweeps
    void AddTopologicalQAV0(Int_t lRecNumberOfSteps = 100);
    void AddTopologicalQACascade(Int_t lRecNumberOfSteps = 100);
    // 3 - Standard analysis configurations + systematics
    void AddStandardV0Configuration();
    void AddStandardCascadeConfiguration();
    void AddCascadeConfiguration276TeV();
    //---------------------------------------------------------------------------------------
    Float_t GetDCAz(AliESDtrack *lTrack);
    Float_t GetCosPA(AliESDtrack *lPosTrack, AliESDtrack *lNegTrack, AliESDEvent *lEvent);
    //---------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------
    // A simple struct to handle FMD hits information
    // Nothe: this struct is based on what is implemented in AliAnalysisTaskValidation
    //        defined as 'Track' (thanks to C. Bourjau). It was slightly changed here and
    //        renamed to 'FMDhit' in order to avoid any confusion.
    struct FMDhit {
        Float_t eta;
        Float_t phi;
        Float_t weight;
        //Constructor
        FMDhit(Float_t _eta, Float_t _phi, Float_t _weight)
            :eta(_eta), phi(_phi), weight(_weight) {};
    };
    typedef std::vector<AliAnalysisTaskStrEffStudy::FMDhit> FMDhits;
//---------------------------------------------------------------------------------------
    AliAnalysisTaskStrEffStudy::FMDhits GetFMDhits(AliAODEvent* aodEvent) const;
//---------------------------------------------------------------------------------------
    Double_t PropagateToDCA(AliESDv0 *v, AliExternalTrackParam *t, AliESDEvent *event, Double_t b);
    //Helper functions
    Double_t Det(Double_t a00, Double_t a01, Double_t a10, Double_t a11) const;
    Double_t Det(Double_t a00,Double_t a01,Double_t a02,
                 Double_t a10,Double_t a11,Double_t a12,
                 Double_t a20,Double_t a21,Double_t a22) const;
    void Evaluate(const Double_t *h, Double_t t,
                  Double_t r[3],  //radius vector
                  Double_t g[3],  //first defivatives
                  Double_t gg[3]); //second derivatives
    Double_t GetErrorInPosition(AliExternalTrackParam *t1) const;
    //---------------------------------------------------------------------------------------
    //Improved DCA V0 Dau
    Double_t GetDCAV0Dau ( AliExternalTrackParam *pt, AliExternalTrackParam *nt, Double_t &xp, Double_t &xn, Double_t b);
    void GetHelixCenter(const AliExternalTrackParam *track,Double_t center[2], Double_t b);
    //---------------------------------------------------------------------------------------

private:
    // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
    // your data member object is created on the worker nodes and streaming is not needed.
    // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14
    TList  *fListHist;      //! List of Cascade histograms
    TList  *fListV0;        // List of Cascade histograms
    TList  *fListCascade;   // List of Cascade histograms
    TTree  *fTreeEvent;              //! Output Tree, Events
    TTree  *fTreeV0;              //! Output Tree, V0s
    TTree  *fTreeCascade;              //! Output Tree, Cascades
    TTree  *fTreeHyperTriton3Body;              //! Output Tree, Cascades

    AliPIDResponse *fPIDResponse;     // PID response object
    AliESDtrackCuts *fESDtrackCuts;   // ESD track cuts used for primary track definition
    AliESDtrackCuts *fESDtrackCutsITSsa2010;  // ESD track cuts used for ITSsa track definition
    AliESDtrackCuts *fESDtrackCutsGlobal2015; // ESD track cuts used for global track definition
    AliAnalysisUtils *fUtils;         // analysis utils (for MV pileup selection)

    //Implementation of event selection utility
    AliEventCuts fEventCuts; /// Event cuts class

    TRandom3 *fRand;

    //Objects Controlling Task Behaviour
    Bool_t fkSaveEventTree;           //if true, save Event TTree
    Bool_t fkSaveV0Tree;              //if true, save TTree
    Bool_t fkDownScaleV0;
    Double_t fDownScaleFactorV0;
    Bool_t fkPreselectDedx;
    Bool_t fkPreselectPID;
    Bool_t fkUseOnTheFlyV0Cascading;
    Bool_t fkDoImprovedCascadeVertexFinding;
    Bool_t fkDoImprovedDCAV0DauPropagation;
    Bool_t fkDoImprovedDCACascDauPropagation;
    Bool_t fkIfImprovedPerformInitialLinearPropag;
    Double_t fkIfImprovedExtraPrecisionFactor;
    Bool_t fkDebugWrongPIDForTracking; //if true, add extra information to TTrees for debugging
    Bool_t fkDebugBump; //if true, add extra information to TTrees for debugging
    Bool_t fkDebugOOBPileup; // if true, add extra information to TTrees for pileup study
    Bool_t fkDoExtraEvSels; //use AliEventCuts for event selection

    Bool_t fkSaveCascadeTree;         //if true, save TTree
    Bool_t fkDownScaleCascade;
    Double_t fDownScaleFactorCascade;

    Float_t fMinPtToSave; //minimum pt above which we keep candidates in TTree output
    Float_t fMaxPtToSave; //maximum pt below which we keep candidates in TTree output

    Bool_t fkSaveHyperTriton3BodyTree;         //if true, save TTree
    Bool_t fkDownScaleHyperTriton3Body;
    Double_t fDownScaleFactorHyperTriton3Body;

    //Objects Controlling Task Behaviour: has to be streamed!
    Bool_t    fkRunVertexers;           // if true, re-run vertexer with loose cuts *** only for CASCADES! ***
    Bool_t    fkUseLightVertexer;       // if true, use AliLightVertexers instead of regular ones
    Bool_t    fkDoV0Refit;              // if true, will invoke AliESDv0::Refit() to improve precision
    Bool_t    fkExtraCleanup;           //if true, perform pre-rejection of useless candidates before going through configs

    //Save only decent tracks
    Bool_t fkSaveGoodTracks;

    //Sandbox mode
    Bool_t fkSandboxV0;
    Bool_t fkSandboxCascade;

    AliVEvent::EOfflineTriggerTypes fTrigType; // trigger type

    Double_t  fV0VertexerSels[7];        // Array to store the 7 values for the different selections V0 related
    Double_t  fCascadeVertexerSels[8];   // Array to store the 8 values for the different selections Casc. related

    Double_t fLambdaMassMean[5]; //Array to store the lambda mass mean parametrization
    //[0]+[1]*TMath::Exp([2]*x)+[3]*TMath::Exp([4]*x)

    Double_t fLambdaMassSigma[4]; //Array to store the lambda mass sigma parametrization
    //[0]+[1]*x+[2]*TMath::Exp([3]*x)

    Double_t fPrecisionCutoffCascadeDCA; //Precision cutoff for GetDCA numerical recipe
    Int_t fMaxIterationsCascadeDCA; //Max N Iter for cascade DCA calculation

//===========================================================================================
//   Variables for Event Tree
//===========================================================================================
    Float_t fCentrality; //!
    Bool_t fMVPileupFlag; //!
    Bool_t fOOBPileupFlag; //!

    //TOF info for OOB pileuo study
    Int_t  fNTOFClusters;  //!
    Int_t  fNTOFMatches;   //!
    Int_t  fNTracksITSsa2010; //!
    Int_t  fNTracksGlobal2015; //!
    Int_t  fNTracksGlobal2015TriggerPP; //!

    //V0 info for OOB pileup study
    Float_t fAmplitudeV0A; //!
    Float_t fAmplitudeV0C; //!

    //FMD info for OOB pileup study
    Float_t fNHitsFMDA; //!
    Float_t fNHitsFMDC; //!



//===========================================================================================
//   Variables for V0 Tree
//===========================================================================================
    
    Bool_t fTreeVariableGoodV0;
    Float_t fTreeVariableCentrality;
    Float_t fTreeVariablePosLength;
    Float_t fTreeVariableNegLength;
    Float_t fTreeVariablePosCrossedRows;
    Float_t fTreeVariableNegCrossedRows;
    ULong64_t fTreeVariablePosTrackStatus;
    ULong64_t fTreeVariableNegTrackStatus;
    Float_t fTreeVariablePosDxy;
    Float_t fTreeVariableNegDxy;
    Float_t fTreeVariablePosDz;
    Float_t fTreeVariableNegDz;
    Float_t fTreeVariableDcaV0Daughters; //!;
    Float_t fTreeVariableDcaV0DaughtersGeometric; //purely geometric, no use of uncertainties
    Bool_t fTreeVariablePosPropagStatus;
    Bool_t fTreeVariableNegPropagStatus;
    Float_t fTreeVariableV0Radius; //!
    Float_t fTreeVariableV0CosineOfPointingAngle; //!
    Float_t fTreeVariableDecayX;
    Float_t fTreeVariableDecayY;
    Float_t fTreeVariableDecayZ;
    Float_t fTreeVariableDecayXMC;
    Float_t fTreeVariableDecayYMC;
    Float_t fTreeVariableDecayZMC;
    //p of daughters
    Float_t fTreeVariableNegPxMC; //!
    Float_t fTreeVariableNegPyMC; //!
    Float_t fTreeVariableNegPzMC; //!
    Float_t fTreeVariablePosPxMC; //!
    Float_t fTreeVariablePosPyMC; //!
    Float_t fTreeVariablePosPzMC; //!
    Float_t fTreeVariableInvMassK0s; //!
    Float_t fTreeVariableInvMassLambda; //!
    Float_t fTreeVariableInvMassAntiLambda; //!
    Int_t fTreeVariablePID;
    Int_t fTreeVariablePIDPositive;
    Int_t fTreeVariablePIDNegative;
    Float_t fTreeVariablePtMC;
    Float_t fTreeVariableRapMC;
    //TOF info
    Float_t fTreeVariableNegTOFSignal; //!
    Float_t fTreeVariablePosTOFSignal; //!
    //Uncertainties
    Float_t fTreeVariablePosAlpha;
    Float_t fTreeVariablePosSigmaY2;
    Float_t fTreeVariablePosSigmaZ2;
    Float_t fTreeVariableNegAlpha;
    Float_t fTreeVariableNegSigmaY2;
    Float_t fTreeVariableNegSigmaZ2;
    
    //Sandbox mode
    AliESDtrack *fTreeVariablePosTrack;
    AliESDtrack *fTreeVariableNegTrack;
    
    AliESDv0 *fTreeVariableOTFV0;
    Bool_t fTreeVariableFoundOTFV0; 
    
    Float_t fTreeVariableMagneticField;
    
    Float_t fTreeVariablePosOriginalX;
    Float_t fTreeVariableNegOriginalX;
    
    Float_t fTreeVariablePVx;
    Float_t fTreeVariablePVy;
    Float_t fTreeVariablePVz;
    
    AliESDVertex *fTreeVariableAliESDvertex;
    
    Int_t fTreeVariableRun;
    
//===========================================================================================
//   Variables for Cascade Candidate Tree
//===========================================================================================
    Float_t fTreeCascVarCentrality;
    Int_t fTreeCascVarPosSign;
    Int_t fTreeCascVarNegSign;
    Int_t fTreeCascVarBachSign;
    Float_t fTreeCascVarPosLength;
    Float_t fTreeCascVarNegLength;
    Float_t fTreeCascVarBachLength;
    Float_t fTreeCascVarPosCrossedRows;
    Float_t fTreeCascVarNegCrossedRows;
    Float_t fTreeCascVarBachCrossedRows;
    //Tracking flags
    ULong64_t fTreeCascVarPosTrackStatus;
    ULong64_t fTreeCascVarNegTrackStatus;
    ULong64_t fTreeCascVarBachTrackStatus;
    //DCAxy to PV
    Float_t fTreeCascVarPosDxy;
    Float_t fTreeCascVarNegDxy;
    Float_t fTreeCascVarBachDxy;
    //DCAz
    Float_t fTreeCascVarPosDz;
    Float_t fTreeCascVarNegDz;
    Float_t fTreeCascVarBachDz;
    Float_t fTreeCascVarDcaV0Daughters;
    Bool_t fTreeCascVarNegPropagStatus;
    Bool_t fTreeCascVarPosPropagStatus;
    Float_t fTreeCascVarV0Radius;
    Float_t fTreeCascVarV0DecayX;
    Float_t fTreeCascVarV0DecayY;
    Float_t fTreeCascVarV0DecayZ;
    Float_t fTreeCascVarV0DecayXMC;
    Float_t fTreeCascVarV0DecayYMC;
    Float_t fTreeCascVarV0DecayZMC;
    Float_t fTreeCascVarV0CosineOfPointingAngle;
    Float_t fTreeCascVarDCAV0ToPrimVtx;
    Float_t fTreeCascVarDCAxyV0ToPrimVtx;
    Float_t fTreeCascVarInvMassLambda;
    Float_t fTreeCascVarInvMassAntiLambda;
    
    Float_t fTreeCascVarDCACascDaughters;
    Bool_t fTreeCascVarCascPropagation;
    
    Float_t fTreeCascVarDecayX;
    Float_t fTreeCascVarDecayY;
    Float_t fTreeCascVarDecayZ;
    Float_t fTreeCascVarDecayXMC;
    Float_t fTreeCascVarDecayYMC;
    Float_t fTreeCascVarDecayZMC;
    Float_t fTreeCascVarCascCosPointingAngle;
    
    Float_t fTreeCascVarInvMassXiMinus;
    Float_t fTreeCascVarInvMassXiPlus;
    Float_t fTreeCascVarInvMassOmegaMinus;
    Float_t fTreeCascVarInvMassOmegaPlus;
    
    Int_t fTreeCascVarPIDPositive;
    Int_t fTreeCascVarPIDNegative;
    Int_t fTreeCascVarPIDBachelor;

    //Set tree variables
    Int_t fTreeCascVarPID;
    Float_t fTreeCascVarPtMC;
    Float_t fTreeCascVarRapMC;
    
    //TOF info
    Float_t fTreeCascVarNegTOFSignal; //!
    Float_t fTreeCascVarPosTOFSignal; //!
    Float_t fTreeCascVarBachTOFSignal; //!

    //Super-control vars
    Float_t fTreeCascVarPosDistanceToTrueDecayPt;
    Float_t fTreeCascVarNegDistanceToTrueDecayPt;
    Float_t fTreeCascVarBachDistanceToTrueDecayPt;
    Float_t fTreeCascVarV0DistanceToTrueDecayPt;
    
    Float_t fTreeCascVarNegPx; //!
    Float_t fTreeCascVarNegPy; //!
    Float_t fTreeCascVarNegPz; //!
    Float_t fTreeCascVarPosPx; //!
    Float_t fTreeCascVarPosPy; //!
    Float_t fTreeCascVarPosPz; //!
    Float_t fTreeCascVarBachPx; //!
    Float_t fTreeCascVarBachPy; //!
    Float_t fTreeCascVarBachPz; //!
    
    Float_t fTreeCascVarNegPxMC; //!
    Float_t fTreeCascVarNegPyMC; //!
    Float_t fTreeCascVarNegPzMC; //!
    Float_t fTreeCascVarPosPxMC; //!
    Float_t fTreeCascVarPosPyMC; //!
    Float_t fTreeCascVarPosPzMC; //!
    Float_t fTreeCascVarBachPxMC; //!
    Float_t fTreeCascVarBachPyMC; //!
    Float_t fTreeCascVarBachPzMC; //!
    
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    //Save full info for full re-vertex offline replay ('sandbox mode')
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    AliESDtrack *fTreeCascVarBachTrack;
    AliESDtrack *fTreeCascVarPosTrack;
    AliESDtrack *fTreeCascVarNegTrack;
    
    //Sandbox on-the-fly V0 for comparison, please
    AliESDv0 *fTreeCascVarOTFV0;
    AliESDv0 *fTreeCascVarOTFV0NegBach;
    AliESDv0 *fTreeCascVarOTFV0PosBach;
    Bool_t fTreeCascVarV0AsOTF; 
    Bool_t fTreeCascVarNegBachAsOTF;
    Bool_t fTreeCascVarPosBachAsOTF;
    
    Float_t fTreeCascVarMagneticField;
    
    Float_t fTreeCascVarBachOriginalX;
    Float_t fTreeCascVarPosOriginalX;
    Float_t fTreeCascVarNegOriginalX;
    
    Float_t fTreeCascVarPVx;
    Float_t fTreeCascVarPVy;
    Float_t fTreeCascVarPVz;
    AliESDVertex *fTreeCascVarAliESDvertex;
    
    
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    
    AliESDtrack *fTreeHyp3BodyVarTracks[3];
    Int_t        fTreeHyp3BodyVarPDGcodes[3];
    
    ULong64_t    fTreeHyp3BodyVarEventId;
    Int_t        fTreeHyp3BodyVarMotherId;

    Float_t      fTreeHyp3BodyVarTruePx;
    Float_t      fTreeHyp3BodyVarTruePy;
    Float_t      fTreeHyp3BodyVarTruePz;

    Float_t      fTreeHyp3BodyVarDecayVx;
    Float_t      fTreeHyp3BodyVarDecayVy;
    Float_t      fTreeHyp3BodyVarDecayVz;
    Float_t      fTreeHyp3BodyVarDecayT;

    Float_t      fTreeHyp3BodyVarPVx;
    Float_t      fTreeHyp3BodyVarPVy;
    Float_t      fTreeHyp3BodyVarPVz;
    Float_t      fTreeHyp3BodyVarPVt;

    Float_t      fTreeHyp3BodyVarMagneticField;

//===========================================================================================
//   Histograms
//===========================================================================================

    TH1D *fHistEventCounter; //!
    TH1D *fHistCentrality; //!

    //Histograms for efficiency denominators (at final analysis level selections)
    //V0s
    TH3D *fHistGeneratedPtVsYVsCentralityK0Short;
    TH3D *fHistGeneratedPtVsYVsCentralityLambda;
    TH3D *fHistGeneratedPtVsYVsCentralityAntiLambda;

    //Cascades
    TH3D *fHistGeneratedPtVsYVsCentralityXiMinus;
    TH3D *fHistGeneratedPtVsYVsCentralityXiPlus;
    TH3D *fHistGeneratedPtVsYVsCentralityOmegaMinus;
    TH3D *fHistGeneratedPtVsYVsCentralityOmegaPlus;
    
    //Hypertriton
    TH3D *fHistGeneratedPtVsYVsCentralityHypTrit;
    TH3D *fHistGeneratedPtVsYVsCentralityAntiHypTrit;
    
    AliAnalysisTaskStrEffStudy(const AliAnalysisTaskStrEffStudy&);            // not implemented
    AliAnalysisTaskStrEffStudy& operator=(const AliAnalysisTaskStrEffStudy&); // not implemented

    ClassDef(AliAnalysisTaskStrEffStudy, 1);
    //1: first implementation
};

#endif
