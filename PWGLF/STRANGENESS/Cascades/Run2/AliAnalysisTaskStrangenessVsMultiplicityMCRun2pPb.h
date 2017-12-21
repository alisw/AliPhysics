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
// Modified version of AliAnalysisTaskCheckCascade.h
// Used bits of code from AliAnalysisTaskCheckPerformanceStrange
//
// --- David Dobrigkeit Chinellato
//
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#ifndef AliAnalysisTaskStrangenessVsMultiplicityMCRun2pPb_H
#define AliAnalysisTaskStrangenessVsMultiplicityMCRun2pPb_H

class TList;
class TH1F;
class TH2F;
class TH3F;
class TVector3;
class THnSparse;
class TRandom3;
class TProfile; 

class AliESDpid;
class AliESDtrackCuts;
class AliAnalysisUtils;
class AliESDEvent;
class AliPhysicsSelection;
class AliESDFMD;
class AliCFContainer;
class AliV0Result;
class AliCascadeResult;
class AliExternalTrackParam; 

//#include "TString.h"
//#include "AliESDtrackCuts.h"
//#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"

class AliAnalysisTaskStrangenessVsMultiplicityMCRun2pPb : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskStrangenessVsMultiplicityMCRun2pPb();
    AliAnalysisTaskStrangenessVsMultiplicityMCRun2pPb(Bool_t lSaveEventTree, Bool_t lSaveV0Tree, Bool_t lSaveCascadeTree, const char *name, TString lExtraOptions = "");
    virtual ~AliAnalysisTaskStrangenessVsMultiplicityMCRun2pPb();

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
    typedef std::vector<AliAnalysisTaskStrangenessVsMultiplicityMCRun2pPb::FMDhit> FMDhits;
//---------------------------------------------------------------------------------------
    AliAnalysisTaskStrangenessVsMultiplicityMCRun2pPb::FMDhits GetFMDhits(AliAODEvent* aodEvent) const;
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

    //Objects Controlling Task Behaviour: has to be streamed!
    Bool_t    fkRunVertexers;           // if true, re-run vertexer with loose cuts *** only for CASCADES! ***
    Bool_t    fkUseLightVertexer;       // if true, use AliLightVertexers instead of regular ones
    Bool_t    fkDoV0Refit;              // if true, will invoke AliESDv0::Refit() to improve precision
    Bool_t    fkExtraCleanup;           //if true, perform pre-rejection of useless candidates before going through configs

    AliVEvent::EOfflineTriggerTypes fTrigType; // trigger type

    Double_t  fV0VertexerSels[7];        // Array to store the 7 values for the different selections V0 related
    Double_t  fCascadeVertexerSels[8];   // Array to store the 8 values for the different selections Casc. related

    Double_t fLambdaMassMean[5]; //Array to store the lambda mass mean parametrization
    //[0]+[1]*TMath::Exp([2]*x)+[3]*TMath::Exp([4]*x)

    Double_t fLambdaMassSigma[4]; //Array to store the lambda mass sigma parametrization
    //[0]+[1]*x+[2]*TMath::Exp([3]*x)

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
    Float_t fTreeVariableChi2V0;         //!
    Float_t fTreeVariableDcaV0Daughters; //!
    Float_t fTreeVariableDcaV0ToPrimVertex; //!
    Float_t fTreeVariableDcaPosToPrimVertex; //!
    Float_t fTreeVariableDcaNegToPrimVertex; //!
    Float_t fTreeVariableV0CosineOfPointingAngle; //!
    Float_t fTreeVariableV0Radius; //!
    Float_t fTreeVariablePt; //!
    Float_t fTreeVariablePtMC; //!
    Float_t fTreeVariableRapK0Short; //!
    Float_t fTreeVariableRapLambda; //!
    Float_t fTreeVariableRapMC; //!
    Float_t fTreeVariableInvMassK0s; //!
    Float_t fTreeVariableInvMassLambda; //!
    Float_t fTreeVariableInvMassAntiLambda; //!
    Float_t fTreeVariableAlphaV0; //!
    Float_t fTreeVariablePtArmV0;//!
    Float_t fTreeVariableNegEta; //!
    Float_t fTreeVariablePosEta; //!

    Float_t fTreeVariableNSigmasPosProton; //!
    Float_t fTreeVariableNSigmasPosPion; //!
    Float_t fTreeVariableNSigmasNegProton; //!
    Float_t fTreeVariableNSigmasNegPion; //!

    Float_t fTreeVariableDistOverTotMom;//!
    Int_t   fTreeVariableLeastNbrCrossedRows;//!
    Float_t fTreeVariableLeastRatioCrossedRowsOverFindable;//!
    Float_t fTreeVariableMaxChi2PerCluster; //!
    Float_t fTreeVariableMinTrackLength; //!

    //Variables for debugging Wrong PID hypothesis in tracking bug
    // more info at: https://alice.its.cern.ch/jira/browse/PWGPP-218
    Int_t fTreeVariablePosPIDForTracking; //! uses AliPID EParticleType code (0=electron, 1=muon, 2=pion, etc)
    Int_t fTreeVariableNegPIDForTracking; //!
    Float_t fTreeVariablePosdEdx; //!
    Float_t fTreeVariableNegdEdx; //!
    Float_t fTreeVariablePosInnerP; //!
    Float_t fTreeVariableNegInnerP; //!
    //Decay Length issue debugging: ULong_t with track status
    ULong64_t fTreeVariableNegTrackStatus; //!
    ULong64_t fTreeVariablePosTrackStatus; //!
    Float_t fTreeVariableNegDCAz; //!
    Float_t fTreeVariablePosDCAz; //!

    //Variables for OOB pileup study (high-multiplicity triggers pp 13 TeV - 2016 data)
    Float_t fTreeVariableNegTOFExpTDiff;      //!
    Float_t fTreeVariablePosTOFExpTDiff;      //!
    //Event info
    Float_t fTreeVariableAmplitudeV0A; //!
    Float_t fTreeVariableAmplitudeV0C; //!
    Float_t fTreeVariableNHitsFMDA; //!
    Float_t fTreeVariableNHitsFMDC; //!

    //Event Multiplicity Variables
    Float_t fTreeVariableCentrality; //!
    Bool_t fTreeVariableMVPileupFlag; //!
    Bool_t fTreeVariableOOBPileupFlag; //!

    //MC exclusive Characteristics: 7, also required for feeddown tests
    Float_t fTreeVariablePtMother; //!
    Float_t fTreeVariableRapMother; //!
    Int_t fTreeVariablePID; //!
    Int_t fTreeVariablePIDPositive; //!
    Int_t fTreeVariablePIDNegative; //!
    Int_t fTreeVariablePIDMother; //!
    Int_t fTreeVariablePrimaryStatus; //!
    Int_t fTreeVariablePrimaryStatusMother; //!

//===========================================================================================
//   Variables for Cascade Candidate Tree
//===========================================================================================
    Int_t fTreeCascVarCharge;         //!
    Float_t fTreeCascVarMassAsXi;     //!
    Float_t fTreeCascVarMassAsOmega;  //!
    Float_t fTreeCascVarPt;           //!
    Float_t fTreeCascVarPtMC;         //!
    Float_t fTreeCascVarRapXi;        //!
    Float_t fTreeCascVarRapOmega;     //!
    Float_t fTreeCascVarRapMC;        //!
    Float_t fTreeCascVarNegEta;       //!
    Float_t fTreeCascVarPosEta;       //!
    Float_t fTreeCascVarBachEta;      //!
    Float_t fTreeCascVarDCACascDaughters; //!
    Float_t fTreeCascVarDCABachToPrimVtx; //!
    Float_t fTreeCascVarDCAV0Daughters;   //!
    Float_t fTreeCascVarDCAV0ToPrimVtx;   //!
    Float_t fTreeCascVarDCAPosToPrimVtx;  //!
    Float_t fTreeCascVarDCANegToPrimVtx;  //!
    Float_t fTreeCascVarCascCosPointingAngle;         //!
    Float_t fTreeCascVarCascRadius;                   //!
    Float_t fTreeCascVarV0Mass;                       //!
    Float_t fTreeCascVarV0CosPointingAngle;           //!
    Float_t fTreeCascVarV0CosPointingAngleSpecial;    //!
    Float_t fTreeCascVarV0Radius;                     //!
    Float_t fTreeCascVarDCABachToBaryon;              //!
    Float_t fTreeCascVarWrongCosPA;                   //!
    Int_t   fTreeCascVarLeastNbrClusters;             //!
    Float_t fTreeCascVarDistOverTotMom;               //!
    Float_t fTreeCascVarMaxChi2PerCluster; //!
    Float_t fTreeCascVarMinTrackLength; //!

    //TPC dEdx
    Float_t fTreeCascVarNegNSigmaPion;   //!
    Float_t fTreeCascVarNegNSigmaProton; //!
    Float_t fTreeCascVarPosNSigmaPion;   //!
    Float_t fTreeCascVarPosNSigmaProton; //!
    Float_t fTreeCascVarBachNSigmaPion;  //!
    Float_t fTreeCascVarBachNSigmaKaon;  //!

    //Variables for debugging Wrong PID hypothesis in tracking bug
    // more info at: https://alice.its.cern.ch/jira/browse/PWGPP-218
    Int_t fTreeCascVarPosPIDForTracking; //! uses AliPID EParticleType code (0=electron, 1=muon, 2=pion, etc)
    Int_t fTreeCascVarNegPIDForTracking; //!
    Int_t fTreeCascVarBachPIDForTracking; //!
    Float_t fTreeCascVarNegInnerP; //!
    Float_t fTreeCascVarPosInnerP; //!
    Float_t fTreeCascVarBachInnerP; //!
    Float_t fTreeCascVarNegdEdx; //!
    Float_t fTreeCascVarPosdEdx; //!
    Float_t fTreeCascVarBachdEdx; //!

    //Decay Length issue debugging: ULong_t with track status
    ULong64_t fTreeCascVarNegTrackStatus; //!
    ULong64_t fTreeCascVarPosTrackStatus; //!
    ULong64_t fTreeCascVarBachTrackStatus; //!
    Float_t fTreeCascVarNegDCAz; //!
    Float_t fTreeCascVarPosDCAz; //!
    Float_t fTreeCascVarBachDCAz; //!

    //Variables for debugging the invariant mass bump
    //Full momentum information
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

    Float_t fTreeCascVarV0DecayX; //!
    Float_t fTreeCascVarV0DecayY; //!
    Float_t fTreeCascVarV0DecayZ; //!
    Float_t fTreeCascVarCascadeDecayX; //!
    Float_t fTreeCascVarCascadeDecayY; //!
    Float_t fTreeCascVarCascadeDecayZ; //!

    Float_t fTreeCascVarV0DecayXMC; //!
    Float_t fTreeCascVarV0DecayYMC; //!
    Float_t fTreeCascVarV0DecayZMC; //!
    Float_t fTreeCascVarCascadeDecayXMC; //!
    Float_t fTreeCascVarCascadeDecayYMC; //!
    Float_t fTreeCascVarCascadeDecayZMC; //!

    //Vars for studying cascade decay point calculations
    Float_t fTreeCascVarBachelorDCAptX; //!
    Float_t fTreeCascVarBachelorDCAptY; //!
    Float_t fTreeCascVarBachelorDCAptZ; //!
    Float_t fTreeCascVarV0DCAptX; //!
    Float_t fTreeCascVarV0DCAptY; //!
    Float_t fTreeCascVarV0DCAptZ; //!
    Float_t fTreeCascVarDCADaughters_Test; //!
    Float_t fTreeCascVarBachelorDCAptSigmaX2; //
    Float_t fTreeCascVarBachelorDCAptSigmaY2; //
    Float_t fTreeCascVarBachelorDCAptSigmaZ2; //
    Float_t fTreeCascVarV0DCAptUncertainty_V0Pos; //
    Float_t fTreeCascVarV0DCAptUncertainty_V0Ang; //
    
    Float_t fTreeCascVarV0DCAptPosSigmaX2; //
    Float_t fTreeCascVarV0DCAptPosSigmaY2; //
    Float_t fTreeCascVarV0DCAptPosSigmaZ2; //
    Float_t fTreeCascVarV0DCAptPosSigmaSnp2; //
    Float_t fTreeCascVarV0DCAptPosSigmaTgl2; //
    
    Float_t fTreeCascVarV0DCAptNegSigmaX2; //
    Float_t fTreeCascVarV0DCAptNegSigmaY2; //
    Float_t fTreeCascVarV0DCAptNegSigmaZ2; //
    Float_t fTreeCascVarV0DCAptNegSigmaSnp2; //
    Float_t fTreeCascVarV0DCAptNegSigmaTgl2; //

    Float_t fTreeCascVarPrimVertexX;
    Float_t fTreeCascVarPrimVertexY;
    Float_t fTreeCascVarPrimVertexZ;
    
    Float_t fTreeCascVarMagField; // for X-checks
    
    Float_t fTreeCascVarV0Lifetime; //! //V0 lifetime (actually, mL/p)
    //Track Labels (check for duplicates, etc)
    Int_t fTreeCascVarNegIndex; //!
    Int_t fTreeCascVarPosIndex; //!
    Int_t fTreeCascVarBachIndex; //!
    Int_t fTreeCascVarNegLabel; //!
    Int_t fTreeCascVarPosLabel; //!
    Int_t fTreeCascVarBachLabel; //!
    Int_t fTreeCascVarNegLabelMother; //!
    Int_t fTreeCascVarPosLabelMother; //!
    Int_t fTreeCascVarBachLabelMother; //!
    Int_t fTreeCascVarNegLabelGrandMother; //!
    Int_t fTreeCascVarPosLabelGrandMother; //!
    Int_t fTreeCascVarBachLabelGrandMother; //!
    //Event Number (check same-event index mixups)
    ULong64_t fTreeCascVarEventNumber; //!

    //Variables for OOB pileup study (high-multiplicity triggers pp 13 TeV - 2016 data)
    Float_t fTreeCascVarNegTOFExpTDiff; //!
    Float_t fTreeCascVarPosTOFExpTDiff; //!
    Float_t fTreeCascVarBachTOFExpTDiff; //!
    //Event info
    Float_t fTreeCascVarAmplitudeV0A; //!
    Float_t fTreeCascVarAmplitudeV0C; //!
    Float_t fTreeCascVarNHitsFMDA; //!
    Float_t fTreeCascVarNHitsFMDC; //!

    //Event Multiplicity Variables
    Float_t fTreeCascVarCentrality; //!
    Bool_t fTreeCascVarMVPileupFlag; //!
    Bool_t fTreeCascVarOOBPileupFlag; //!

    //MC-only Variabless
    Int_t   fTreeCascVarIsPhysicalPrimary; //!
    Int_t   fTreeCascVarPID;         //!
    Int_t   fTreeCascVarPIDNegative;         //!
    Int_t   fTreeCascVarPIDPositive;         //!
    Int_t   fTreeCascVarPIDBachelor;         //!
    Int_t   fTreeCascVarPIDNegativeMother;         //!
    Int_t   fTreeCascVarPIDPositiveMother;         //!
    Int_t   fTreeCascVarPIDBachelorMother;         //!
    Int_t   fTreeCascVarPIDNegativeGrandMother;         //!
    Int_t   fTreeCascVarPIDPositiveGrandMother;         //!
    Int_t   fTreeCascVarPIDBachelorGrandMother;         //!

    Bool_t fTreeCascVarIsPhysicalPrimaryNegative;
    Bool_t fTreeCascVarIsPhysicalPrimaryPositive;
    Bool_t fTreeCascVarIsPhysicalPrimaryBachelor;
    Bool_t fTreeCascVarIsPhysicalPrimaryNegativeMother;
    Bool_t fTreeCascVarIsPhysicalPrimaryPositiveMother;
    Bool_t fTreeCascVarIsPhysicalPrimaryBachelorMother;
    Bool_t fTreeCascVarIsPhysicalPrimaryNegativeGrandMother;
    Bool_t fTreeCascVarIsPhysicalPrimaryPositiveGrandMother;
    Bool_t fTreeCascVarIsPhysicalPrimaryBachelorGrandMother;

    //Well, why not? Let's give it a shot
    Int_t   fTreeCascVarSwappedPID;         //!

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

    AliAnalysisTaskStrangenessVsMultiplicityMCRun2pPb(const AliAnalysisTaskStrangenessVsMultiplicityMCRun2pPb&);            // not implemented
    AliAnalysisTaskStrangenessVsMultiplicityMCRun2pPb& operator=(const AliAnalysisTaskStrangenessVsMultiplicityMCRun2pPb&); // not implemented

    ClassDef(AliAnalysisTaskStrangenessVsMultiplicityMCRun2pPb, 1);
    //1: first implementation
};

#endif
