
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

#ifndef AliAnalysisTaskStrangenessVsMultiplicityEERun2_H
#define AliAnalysisTaskStrangenessVsMultiplicityEERun2_H

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
class AliESDtrack;
class AliPhysicsSelection;
class AliPIDResponse;
class AliCFContainer;
class AliV0Result;
class AliCascadeResult;
class AliExternalTrackParam;

//#include "TString.h"
//#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"


class AliAnalysisTaskStrangenessVsMultiplicityEERun2 : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskStrangenessVsMultiplicityEERun2();
    AliAnalysisTaskStrangenessVsMultiplicityEERun2(Bool_t lSaveEventTree, Bool_t lSaveV0Tree, Bool_t lSaveCascadeTree, const char *name, TString lExtraOptions = "");
    virtual ~AliAnalysisTaskStrangenessVsMultiplicityEERun2();

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
    void SetUseOnTheFlyV0Cascading( Bool_t lUseOnTheFlyV0Cascading = kTRUE ){
        //Highly experimental, use with care!
        fkUseOnTheFlyV0Cascading = lUseOnTheFlyV0Cascading;
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
    void SetPileupRejectionMode ( Int_t lMode = 1 ){
        //mode switch
        // 0 -> no rejection
        // 1 -> Ionut
        // 2 -> Anti-Ionut
        fkPileupRejectionMode = lMode;
    }
    void SetUseOldCentrality ( Bool_t lUseOldCent = kTRUE) {
        fkUseOldCentrality = lUseOldCent;
    }
//---------------------------------------------------------------------------------------
    void SetSelectCharge ( Int_t lCharge = -1) {
        fkSelectCharge = lCharge;
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
    void AddStandardV0Configuration(Bool_t lUseFull = kFALSE, Bool_t lDoSweepLooseTight = kFALSE, Int_t lSweepFullNumb = 0);
    void AddStandardV0RadiusSweep();
    void AddStandardCascadeConfiguration(Bool_t lUseFull = kFALSE, Bool_t lDoSystematics = kTRUE);
    void AddCascadeConfiguration276TeV(); //Adds old 2.76 PbPb cut level analyses
    void AddCascadeConfigurationPreliminaryCrosscheck(); //
//---------------------------------------------------------------------------------------
    Float_t GetDCAz(AliESDtrack *lTrack);
    Float_t GetCosPA(AliESDtrack *lPosTrack, AliESDtrack *lNegTrack, AliESDEvent *lEvent);
//---------------------------------------------------------------------------------------


private:
    // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
    // your data member object is created on the worker nodes and streaming is not needed.
    // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14
    TList  *fListHist;      //! List of Cascade histograms
    TList  *fListK0Short;        // List of Cascade histograms
    TList  *fListLambda;        // List of Cascade histograms
    TList  *fListAntiLambda;        // List of Cascade histograms
    TList  *fListXiMinus;   // List of XiMinus outputs
    TList  *fListXiPlus;   // List of XiPlus outputs
    TList  *fListOmegaMinus;   // List of XiMinus outputs
    TList  *fListOmegaPlus;   // List of XiPlus outputs
    TTree  *fTreeEvent;              //! Output Tree, Events
    TTree  *fTreeV0;              //! Output Tree, V0s
    TTree  *fTreeCascade;              //! Output Tree, Cascades

    AliPIDResponse *fPIDResponse;     //! PID response object
    AliESDtrackCuts *fESDtrackCuts;   //! ESD track cuts used for primary track definition
    AliESDtrackCuts *fESDtrackCutsITSsa2010;  //! ESD track cuts used for ITSsa track definition
    AliESDtrackCuts *fESDtrackCutsGlobal2015; //! ESD track cuts used for global track definition
    AliAnalysisUtils *fUtils;         //! analysis utils (for MV pileup selection)
    
    AliEventCuts fEventCuts;                 /// Event cuts class
    AliEventCuts fEventCutsStrictAntipileup; /// Event cuts class

    TRandom3 *fRand; //!

    //Objects Controlling Task Behaviour
    Bool_t fkSaveEventTree;           //if true, save Event TTree
    Bool_t fkSaveV0Tree;              //if true, save TTree
    Bool_t fkDownScaleV0;
    Double_t fDownScaleFactorV0;
    Bool_t fkPreselectDedx;
    Bool_t fkUseOnTheFlyV0Cascading;
    Bool_t fkDebugWrongPIDForTracking; //if true, add extra information to TTrees for debugging
    Bool_t fkDebugBump; //if true, add extra information to TTrees for debugging
    Bool_t fkDebugOOBPileup; // if true, add extra information to TTrees for pileup study
    Bool_t fkDoExtraEvSels; //if true, rely on AliEventCuts
    Int_t fkPileupRejectionMode; //pileup rejection mode (0=none, 1=ionut, 2=anti-ionut)
    Bool_t fkUseOldCentrality; //if true, use AliCentrality instead of AliMultSelection 
    Bool_t fkDebugZDCInfo; //if true, add extra information from ZDC
    Bool_t fkSwitchOffITSClusters; //if true not interested in ITS Cluster info

    Bool_t fkSaveCascadeTree;         //if true, save TTree
    Bool_t fkDownScaleCascade;
    Double_t fDownScaleFactorCascade;

    //Objects Controlling Task Behaviour: has to be streamed!
    Bool_t    fkRunVertexers;           // if true, re-run vertexer with loose cuts *** only for CASCADES! ***
    Bool_t    fkUseLightVertexer;       // if true, use AliLightVertexers instead of regular ones
    Bool_t    fkDoV0Refit;              // if true, will invoke AliESDv0::Refit in the vertexing procedure
    Bool_t    fkExtraCleanup;           //if true, perform pre-rejection of useless candidates before going through configs

    AliVEvent::EOfflineTriggerTypes fTrigType; // trigger type

    Double_t  fV0VertexerSels[7];        // Array to store the 7 values for the different selections V0 related
    Double_t  fCascadeVertexerSels[8];   // Array to store the 8 values for the different selections Casc. related

    Double_t fLambdaMassMean[5]; //Array to store the lambda mass mean parametrization
    //[0]+[1]*TMath::Exp([2]*x)+[3]*TMath::Exp([4]*x)

    Double_t fLambdaMassSigma[4]; //Array to store the lambda mass sigma parametrization
    //[0]+[1]*x+[2]*TMath::Exp([3]*x)

    Float_t fMinPtToSave; //minimum pt above which we keep candidates in TTree output
    Float_t fMaxPtToSave; //maximum pt below which we keep candidates in TTree output

    //if true, save sandbox mode info (beware large files!)
    Bool_t fkSandboxMode; 
    
//===========================================================================================
//   Variables for Event Tree
//===========================================================================================
   	Float_t fZNApp;//!
    Float_t fZNCpp;//!
    Float_t fZPApp;//!
    Float_t fZPCpp;//!
    Float_t fCentrality; //!
    Bool_t fMVPileupFlag; //!
    Bool_t fOOBPileupFlag; //!
    Float_t fTestVariable; //!
    Int_t fRun;//!

    //TOF info for OOB pileuo study
    Int_t  fNTOFClusters;  //!
    Int_t  fNTOFMatches;   //!
    Int_t  fNTracksITSsa2010; //!
    Int_t  fNTracksGlobal2015; //!
    Int_t  fNTracksGlobal2015TriggerPP; //!

    //V0 info for OOB pileup study
    Float_t fAmplitudeV0A; //!
    Float_t fAmplitudeV0C; //!

    //IR info for OOB pileup study
    Int_t fClosestNonEmptyBC; //!

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

    Float_t fTreeVariableRapK0Short; //!
    Float_t fTreeVariableRapLambda; //!
    Float_t fTreeVariableInvMassK0s; //!
    Float_t fTreeVariableInvMassLambda; //!
    Float_t fTreeVariableInvMassAntiLambda; //!
    Float_t fTreeVariableAlphaV0; //!
    Float_t fTreeVariablePtArmV0;//!
    Float_t fTreeVariableNegEta; //!
    Float_t fTreeVariablePosEta; //!
    Float_t fTreeVariablePosPx; //!
	Float_t fTreeVariablePosPy; //!
	Float_t fTreeVariablePosPz; //!
	Float_t fTreeVariableNegPx; //!
	Float_t fTreeVariableNegPy; //!
	Float_t fTreeVariableNegPz; //!

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

    //Cluster information for all daughter tracks
    Bool_t fTreeVariablePosITSClusters0;
    Bool_t fTreeVariablePosITSClusters1;
    Bool_t fTreeVariablePosITSClusters2;
    Bool_t fTreeVariablePosITSClusters3;
    Bool_t fTreeVariablePosITSClusters4;
    Bool_t fTreeVariablePosITSClusters5;
    
    Bool_t fTreeVariableNegITSClusters0;
    Bool_t fTreeVariableNegITSClusters1;
    Bool_t fTreeVariableNegITSClusters2;
    Bool_t fTreeVariableNegITSClusters3;
    Bool_t fTreeVariableNegITSClusters4;
    Bool_t fTreeVariableNegITSClusters5;
    
    //Cluster information for all daughter tracks
    Bool_t fTreeVariablePosITSSharedClusters0;
    Bool_t fTreeVariablePosITSSharedClusters1;
    Bool_t fTreeVariablePosITSSharedClusters2;
    Bool_t fTreeVariablePosITSSharedClusters3;
    Bool_t fTreeVariablePosITSSharedClusters4;
    Bool_t fTreeVariablePosITSSharedClusters5;
    
    Bool_t fTreeVariableNegITSSharedClusters0;
    Bool_t fTreeVariableNegITSSharedClusters1;
    Bool_t fTreeVariableNegITSSharedClusters2;
    Bool_t fTreeVariableNegITSSharedClusters3;
    Bool_t fTreeVariableNegITSSharedClusters4;
    Bool_t fTreeVariableNegITSSharedClusters5;
    
    Bool_t fTreeVariableIsCowboy; //store if V0 is cowboy-like or sailor-like in XY plane

    //Variables for OOB pileup study (high-multiplicity triggers pp 13 TeV - 2016 data)
    Float_t fTreeVariableNegTOFExpTDiff; //!
    Float_t fTreeVariablePosTOFExpTDiff; //!
    Float_t fTreeVariableNegTOFSignal; //!
    Float_t fTreeVariablePosTOFSignal; //!
    Int_t   fTreeVariableNegTOFBCid; //!
    Int_t   fTreeVariablePosTOFBCid; //! 
    //Event info
    Float_t fTreeVariableAmplitudeV0A; //!
    Float_t fTreeVariableAmplitudeV0C; //!
    Int_t   fTreeVariableClosestNonEmptyBC; //!

    //Event Multiplicity Variables
    Float_t fTreeVariableCentrality; //!
    Float_t fTreeVariableZNApp;//!
    Float_t fTreeVariableZNCpp;//!
    Float_t fTreeVariableZPApp;//!
    Float_t fTreeVariableZPCpp;//!
    Bool_t fTreeVariableMVPileupFlag; //!
    Bool_t fTreeVariableOOBPileupFlag; //!
    
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    //Sandbox V0
    Float_t fTreeVariablePrimVertexX;
    Float_t fTreeVariablePrimVertexY;
    Float_t fTreeVariablePrimVertexZ;
    
    AliExternalTrackParam *fTreeVariablePosTrack; //!
    AliExternalTrackParam *fTreeVariableNegTrack; //!
    
    Float_t fTreeVariableMagneticField;
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

//===========================================================================================
//   Variables for Cascade Candidate Tree
//===========================================================================================
    Int_t fTreeCascVarCharge;         //!
    Float_t fTreeCascVarMassAsXi;     //!
    Float_t fTreeCascVarMassAsOmega;  //!
    Float_t fTreeCascVarPt;           //!
    Float_t fTreeCascVarRapXi;        //!
    Float_t fTreeCascVarRapOmega;     //!
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
    Float_t fTreeCascVarCascDCAtoPVxy;         //!
    Float_t fTreeCascVarCascDCAtoPVz;         //!
    Float_t fTreeCascVarCascRadius;                   //!
    Float_t fTreeCascVarV0Mass;                       //!
    Float_t fTreeCascVarV0MassLambda;                       //!
    Float_t fTreeCascVarV0MassAntiLambda;                       //!
    Float_t fTreeCascVarV0CosPointingAngle;           //!
    Float_t fTreeCascVarV0CosPointingAngleSpecial;    //!
    Float_t fTreeCascVarV0Radius;                     //!
    Float_t fTreeCascVarDCABachToBaryon;                     //!
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
    
    //TOF (experimental, not corrected for weak decay traj)
    Float_t fTreeCascVarNegTOFNSigmaPion;   //!
    Float_t fTreeCascVarNegTOFNSigmaProton; //!
    Float_t fTreeCascVarPosTOFNSigmaPion;   //!
    Float_t fTreeCascVarPosTOFNSigmaProton; //!
    Float_t fTreeCascVarBachTOFNSigmaPion;  //!
    Float_t fTreeCascVarBachTOFNSigmaKaon;  //!
    
    Float_t fTreeCascVarNegITSNSigmaPion;   //!
    Float_t fTreeCascVarNegITSNSigmaProton; //!
    Float_t fTreeCascVarPosITSNSigmaPion;   //!
    Float_t fTreeCascVarPosITSNSigmaProton; //!
    Float_t fTreeCascVarBachITSNSigmaPion;  //!
    Float_t fTreeCascVarBachITSNSigmaKaon;  //!
    
    //ChiSquares
    Float_t fTreeCascVarChiSquareV0;
    Float_t fTreeCascVarChiSquareCascade;
    
    //Extended information: uncertainties at point closest to pV
    Float_t fTreeCascVarBachDCAPVSigmaX2; //
    Float_t fTreeCascVarBachDCAPVSigmaY2; //
    Float_t fTreeCascVarBachDCAPVSigmaZ2; //
    Float_t fTreeCascVarPosDCAPVSigmaX2; //
    Float_t fTreeCascVarPosDCAPVSigmaY2; //
    Float_t fTreeCascVarPosDCAPVSigmaZ2; //
    Float_t fTreeCascVarNegDCAPVSigmaX2; //
    Float_t fTreeCascVarNegDCAPVSigmaY2; //
    Float_t fTreeCascVarNegDCAPVSigmaZ2; //

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

    Float_t fTreeCascVarV0DecayX; //!
    Float_t fTreeCascVarV0DecayY; //!
    Float_t fTreeCascVarV0DecayZ; //!
    Float_t fTreeCascVarCascadeDecayX; //!
    Float_t fTreeCascVarCascadeDecayY; //!
    Float_t fTreeCascVarCascadeDecayZ; //!
    Float_t fTreeCascVarV0Lifetime; //! //V0 lifetime (actually, mL/p)
    //Track Labels (check for duplicates, etc)
    Int_t fTreeCascVarNegIndex; //!
    Int_t fTreeCascVarPosIndex; //!
    Int_t fTreeCascVarBachIndex; //!
    //Event Number (check same-event index mixups)
    ULong64_t fTreeCascVarEventNumber; //!

    //Cluster information for all daughter tracks
    Bool_t fTreeCascVarPosITSClusters0;
    Bool_t fTreeCascVarPosITSClusters1;
    Bool_t fTreeCascVarPosITSClusters2;
    Bool_t fTreeCascVarPosITSClusters3;
    Bool_t fTreeCascVarPosITSClusters4;
    Bool_t fTreeCascVarPosITSClusters5;
    
    Bool_t fTreeCascVarNegITSClusters0;
    Bool_t fTreeCascVarNegITSClusters1;
    Bool_t fTreeCascVarNegITSClusters2;
    Bool_t fTreeCascVarNegITSClusters3;
    Bool_t fTreeCascVarNegITSClusters4;
    Bool_t fTreeCascVarNegITSClusters5;
    
    Bool_t fTreeCascVarBachITSClusters0;
    Bool_t fTreeCascVarBachITSClusters1;
    Bool_t fTreeCascVarBachITSClusters2;
    Bool_t fTreeCascVarBachITSClusters3;
    Bool_t fTreeCascVarBachITSClusters4;
    Bool_t fTreeCascVarBachITSClusters5;
    
    //Cluster information for all daughter tracks
    Bool_t fTreeCascVarPosITSSharedClusters0;
    Bool_t fTreeCascVarPosITSSharedClusters1;
    Bool_t fTreeCascVarPosITSSharedClusters2;
    Bool_t fTreeCascVarPosITSSharedClusters3;
    Bool_t fTreeCascVarPosITSSharedClusters4;
    Bool_t fTreeCascVarPosITSSharedClusters5;
    
    Bool_t fTreeCascVarNegITSSharedClusters0;
    Bool_t fTreeCascVarNegITSSharedClusters1;
    Bool_t fTreeCascVarNegITSSharedClusters2;
    Bool_t fTreeCascVarNegITSSharedClusters3;
    Bool_t fTreeCascVarNegITSSharedClusters4;
    Bool_t fTreeCascVarNegITSSharedClusters5;
    
    Bool_t fTreeCascVarBachITSSharedClusters0;
    Bool_t fTreeCascVarBachITSSharedClusters1;
    Bool_t fTreeCascVarBachITSSharedClusters2;
    Bool_t fTreeCascVarBachITSSharedClusters3;
    Bool_t fTreeCascVarBachITSSharedClusters4;
    Bool_t fTreeCascVarBachITSSharedClusters5;

    //Variables for OOB pileup study (high-multiplicity triggers pp 13 TeV - 2016 data)
    Float_t fTreeCascVarNegTOFExpTDiff; //!
    Float_t fTreeCascVarPosTOFExpTDiff; //!
    Float_t fTreeCascVarBachTOFExpTDiff; //!
    Float_t fTreeCascVarNegTOFSignal; //!
    Float_t fTreeCascVarPosTOFSignal; //!
    Float_t fTreeCascVarBachTOFSignal; //!
    Int_t   fTreeCascVarNegTOFBCid; //!
    Int_t   fTreeCascVarPosTOFBCid; //!
    Int_t   fTreeCascVarBachTOFBCid; //!
    //Event info
    Float_t fTreeCascVarAmplitudeV0A; //!
    Float_t fTreeCascVarAmplitudeV0C; //!
    Int_t   fTreeCascVarClosestNonEmptyBC; //!

    //Event Multiplicity Variables
    Float_t fTreeCascVarCentrality; //!
    Float_t fTreeCascVarZNApp;//!
    Float_t fTreeCascVarZNCpp;//!
    Float_t fTreeCascVarZPApp;//!
    Float_t fTreeCascVarZPCpp;//!
    Bool_t fTreeCascVarMVPileupFlag; //!
    Bool_t fTreeCascVarOOBPileupFlag; //!
    
    //Kink tagging
    Bool_t fTreeCascVarBachIsKink;
    Bool_t fTreeCascVarPosIsKink;
    Bool_t fTreeCascVarNegIsKink;
    
    //Cowboy/sailor studies
    Bool_t  fTreeCascVarIsCowboy;   //store if V0 is cowboy-like or sailor-like in XY plane
    Float_t fTreeCascVarCowboyness; //negative -> cowboy, positive -> sailor
    Bool_t  fTreeCascVarIsCascadeCowboy;   //store if V0 is cowboy-like or sailor-like in XY plane
    Float_t fTreeCascVarCascadeCowboyness; //negative -> cowboy, positive -> sailor
    
    //Select charge (testing / checks)
    Int_t fkSelectCharge;
    
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    //Sandbox mode for cascades
    Float_t fTreeCascVarPrimVertexX;
    Float_t fTreeCascVarPrimVertexY;
    Float_t fTreeCascVarPrimVertexZ;
    
    AliExternalTrackParam *fTreeCascVarBachTrack;//!
    AliExternalTrackParam *fTreeCascVarPosTrack; //!
    AliExternalTrackParam *fTreeCascVarNegTrack; //!
    
    Float_t fTreeCascVarMagneticField;
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

//===========================================================================================
//   Histograms
//===========================================================================================

    TH1D *fHistEventCounter; //!
    TH1D *fHistEventCounterDifferential; //!
    TH1D *fHistCentrality; //!

    AliAnalysisTaskStrangenessVsMultiplicityEERun2(const AliAnalysisTaskStrangenessVsMultiplicityEERun2&);            // not implemented
    AliAnalysisTaskStrangenessVsMultiplicityEERun2& operator=(const AliAnalysisTaskStrangenessVsMultiplicityEERun2&); // not implemented

    ClassDef(AliAnalysisTaskStrangenessVsMultiplicityEERun2, 4);
    //1: first implementation
};

#endif

 
