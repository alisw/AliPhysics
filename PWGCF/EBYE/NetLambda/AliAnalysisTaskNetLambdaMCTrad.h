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

#ifndef AliAnalysisTaskNetLambdaMCTrad_H
#define AliAnalysisTaskNetLambdaMCTrad_H

class TList;
class TH1F;
class TH2F;
class TH3F;
class TVector3;
class THnSparse;
class TRandom3;

class AliESDpid;
class AliESDtrackCuts;
class AliAnalysisUtils;
class AliESDEvent;
class AliPhysicsSelection;
class AliCFContainer;
class AliV0Result; 
class AliCascadeResult;

//#include "TString.h"
//#include "AliESDtrackCuts.h"
//#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskNetLambdaMCTrad : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskNetLambdaMCTrad();
    AliAnalysisTaskNetLambdaMCTrad(Bool_t lSaveEventTree, Bool_t lSaveV0Tree, const char *name, TString lExtraOptions = "");
    virtual ~AliAnalysisTaskNetLambdaMCTrad();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    Double_t MyRapidity(Double_t rE, Double_t rPz) const;
 
    void SetSaveV0s                (Bool_t lSaveV0s        = kTRUE ) {
        fkSaveV0Tree        = lSaveV0s;
    }
    void SetPreselectDedx (Bool_t lPreselectDedx= kTRUE ) {
        fkPreselectDedx   = lPreselectDedx;
    }
    void SetPreselectPID (Bool_t lPreselectPID = kTRUE ) {
        fkPreselectPID   = lPreselectPID;
    }
    
    
//---------------------------------------------------------------------------------------
    //Task Configuration: trigger selection 
    void SetSelectedTriggerClass(AliVEvent::EOfflineTriggerTypes trigType) { fTrigType = trigType;}

    //Task Configuration: Skip Event Selections after trigger (VZERO test)
    void SetDownScaleV0 ( Bool_t lOpt = kTRUE, Float_t lVal = 0) {
        fkDownScaleV0 = lOpt;
        fDownScaleFactorV0 = lVal;
    }
    

    void SetMinPt     ( Float_t lMinPt ) {
        fMinPtToSave = lMinPt;
    }
    void SetMaxPt     ( Float_t lMaxPt ) {
        fMaxPtToSave = lMaxPt;
    }
    
    void SetNumberOfPtBins(Int_t nPtBins){ fNptBins = nPtBins;}

 
    Float_t GetDCAz(AliESDtrack *lTrack);
    Float_t GetCosPA(AliESDtrack *lPosTrack, AliESDtrack *lNegTrack, AliESDEvent *lEvent);
    //---------------------------------------------------------------------------------------

private:
    // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
    // your data member object is created on the worker nodes and streaming is not needed.
    // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14
    TList  *fListHist;      //! List of Cascade histograms
    TList  *fListV0;        // List of Cascade histograms
    TTree  *fTreeEvent;              //! Output Tree, Events
    TTree  *fTreeV0;              //! Output Tree, V0s

    AliPIDResponse *fPIDResponse;     // PID response object
    AliESDtrackCuts *fESDtrackCuts;   // ESD track cuts used for primary track definition
    AliAnalysisUtils *fUtils;         // analysis utils (for MV pileup selection)
    TRandom3 *fRand; 

    //Objects Controlling Task Behaviour
    Bool_t fkSaveEventTree;           //if true, save Event TTree
    Bool_t fkSaveV0Tree;              //if true, save TTree
    Bool_t fkDownScaleV0;
    Double_t fDownScaleFactorV0;
    Bool_t fkPreselectDedx;
    Bool_t fkPreselectPID;
    Bool_t fkUseOnTheFlyV0Cascading; 
    Bool_t fkDebugWrongPIDForTracking; //if true, add extra information to TTrees for debugging
    Bool_t fkDebugBump; //if true, add extra information to TTrees for debugging
    
    
    Float_t fMinPtToSave; //minimum pt above which we keep candidates in TTree output
    Float_t fMaxPtToSave; //maximum pt below which we keep candidates in TTree output
    
    //Objects Controlling Task Behaviour: has to be streamed!
    
    AliVEvent::EOfflineTriggerTypes fTrigType; // trigger type
    

//===========================================================================================
//   Variables for Event Tree
//===========================================================================================
    Float_t fCentrality; //!
    Bool_t fMVPileupFlag; //!

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
    
    //Event Multiplicity Variables
    Float_t fTreeVariableCentrality; //!
    Bool_t fTreeVariableMVPileupFlag;         //!
    
    //MC exclusive Characteristics: 7, also required for feeddown tests
    Float_t fTreeVariablePtMother; //!
    Float_t fTreeVariableRapMother; //!
    Int_t fTreeVariablePID; //!
    Int_t fTreeVariablePIDPositive; //!
    Int_t fTreeVariablePIDNegative; //!
    Int_t fTreeVariablePIDMother; //!
    Int_t fTreeVariablePrimaryStatus; //!
    Int_t fTreeVariablePrimaryStatusMother; //!

    Bool_t     fIsMC;                         // Is MC event - Auto set by Add Task
    Int_t      fNptBins;                     //pt bin numbers
    Double_t        *fPtArray;       //! pt array

    Int_t            fPidType;                  // 0=charge, 1=pion, 2=kaon, 3=proton

    Int_t    GetPtBin(Double_t pt); //to get the bin number of the pt
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
    
    //addded hists
    
    TH3F* f3dHistInvMassVsPtVsCentLambda;
    TH3F* f3dHistInvMassVsPtVsCentAntiLambda;
    TH2F* fhistLambdaRecPT;
    TH2F* fhistAntiLambdaRecPT;
    TH2F* fhistLambdaRecpri;
    TH2F* fhistAntiLambdaRecpri;
    TH2F* fHistGeneratedCentralityVsLambda;
    TH2F* fHistGeneratedCentralityVsAntiLambda;  
    
    THnSparse *fPtBinNplusNminusCh;
    THnSparse *fPtBinNplusNminusChTruth;

    
   
    
    AliAnalysisTaskNetLambdaMCTrad(const AliAnalysisTaskNetLambdaMCTrad&);            // not implemented
    AliAnalysisTaskNetLambdaMCTrad& operator=(const AliAnalysisTaskNetLambdaMCTrad&); // not implemented

    ClassDef(AliAnalysisTaskNetLambdaMCTrad, 1);
    //1: first implementation
};

#endif
