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

#ifndef AliAnalysisTaskStrangenessVsMultiplicity_H
#define AliAnalysisTaskStrangenessVsMultiplicity_H

class TList;
class TH1F;
class TH2F;
class TH3F;
class TVector3;
class THnSparse;

class AliESDpid;
class AliESDtrackCuts;
class AliAnalysisUtils;
class AliPPVsMultUtils;
class AliESDEvent;
class AliPhysicsSelection;
class AliCFContainer;

//#include "TString.h"
//#include "AliESDtrackCuts.h"
//#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskStrangenessVsMultiplicity : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskStrangenessVsMultiplicity();
    AliAnalysisTaskStrangenessVsMultiplicity(const char *name);
    virtual ~AliAnalysisTaskStrangenessVsMultiplicity();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    Double_t MyRapidity(Double_t rE, Double_t rPz) const;
    Bool_t   IsINELgtZERO(AliESDEvent *lESDevent, TString lType) const;
    Int_t GetNumberTPCtracks(AliESDEvent *event) const;
 
    void SetSaveV0s                (Bool_t lSaveV0s        = kTRUE ) {
        fkSaveV0Tree        = lSaveV0s;
    }
    void SetSaveCascades           (Bool_t lSaveCascades   = kTRUE ) {
        fkSaveCascadeTree   = lSaveCascades;
    }
    //Output reduction measures
    void SetSaveLambda (Bool_t lSaveLambda = kTRUE ) {
        fkSaveLambda   = lSaveLambda;
    }
    void SetSaveAntiLambda (Bool_t lSaveAntiLambda = kTRUE ) {
        fkSaveAntiLambda   = lSaveAntiLambda;
    }
    void SetSaveK0Short (Bool_t lSaveK0Short = kTRUE ) {
        fkSaveK0Short   = lSaveK0Short;
    }
    void SetPreselectDedx (Bool_t lPreselectDedx= kTRUE ) {
        fkPreselectDedx   = lPreselectDedx;
    }

    void SetApplySPDClsVsTrackletsCut(Bool_t lSPDClsVsTrk = kTRUE) {
        fkApplyTrackletsVsClustersCut = lSPDClsVsTrk;
    }

    Bool_t SelectEventsMultiplicityDependentAnalysis();
    Bool_t SelectEventsMinimumBiasAnalysis();

//---------------------------------------------------------------------------------------
    //Task Configuration: trigger selection 
    void SetSelectedTriggerClass(AliVEvent::EOfflineTriggerTypes trigType) { fTrigType = trigType;}
    void SetSelectedTriggerClass(TString trigName) { fkSelectTriggerByName = kTRUE; fTrigName = trigName;}
    void SetMininumBiasAnalysis(Bool_t lMinbAn) { fMinbAnalysis = lMinbAn;} 
 
//---------------------------------------------------------------------------------------
    //Task Configuration: Meant to enable quick re-execution of vertexer if needed
    void SetRunVertexers ( Bool_t lRunVertexers = kTRUE) {
        fkRunVertexers = lRunVertexers;
    }
//---------------------------------------------------------------------------------------
    //Task Configuration: Skip Event Selections after trigger (VZERO test)
    void SetSkipEventSelection ( Bool_t lSkipEventSelection = kTRUE) {
        fkSkipEventSelection = lSkipEventSelection;
    }
//---------------------------------------------------------------------------------------
    //Task Configuration: Skip Event Selections after trigger (VZERO test)
    void SetUseMultSelection ( Bool_t lUseMultSelection = kTRUE) {
        fkMultSelection = lUseMultSelection;
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

private:
    // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
    // your data member object is created on the worker nodes and streaming is not needed.
    // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14
    TList  *fListHist;  //! List of Cascade histograms
    TTree  *fTreeEvent;              //! Output Tree, Events
    TTree  *fTreeV0;              //! Output Tree, V0s
    TTree  *fTreeCascade;              //! Output Tree, Cascades

    AliPIDResponse *fPIDResponse;     // PID response object
    AliESDtrackCuts *fESDtrackCuts;   // ESD track cuts used for primary track definition
    AliPPVsMultUtils *fPPVsMultUtils; //
    AliAnalysisUtils *fUtils;         // analysis utils (for MV pileup selection)

    //Objects Controlling Task Behaviour
    Bool_t fkSaveV0Tree;              //if true, save TTree
    Bool_t fkSaveCascadeTree;         //if true, save TTree
    Bool_t fkSaveLambda;            //if true, removes tree fill confirmation for lambda mass window
    Bool_t fkSaveAntiLambda;        //if true, removes tree fill confirmation for antilambda mass window
    Bool_t fkSaveK0Short;           //if true, removes tree fill confirmation for k0short mass window
    Bool_t fkPreselectDedx;         //if true, applies a 7-sigma TPC dEdx band 
    Bool_t fkSaveExtendedRefMultInfo; //if true, save 20 integers per event extra for eta-dependence of refmult
    Bool_t fkSelectTriggerByName;  // to select trigger by name (if it's not availble in AliVEvent)
    Bool_t fkMultSelection; //if true, will use new multiplicity selection framework

    //Objects Controlling Task Behaviour: has to be streamed!
    Bool_t    fkRunVertexers;           // if true, re-run vertexer with loose cuts *** only for CASCADES! ***
    Bool_t    fkSkipEventSelection;     // if true, will only perform TRIGGER selection (currently kMB, to change)
    Bool_t    fkApplyTrackletsVsClustersCut; //if true, applies Tracklet vs clusters cut together with PS
    Bool_t    fMinbAnalysis;                // if true event selection is done according to minb analysis 
    AliVEvent::EOfflineTriggerTypes fTrigType; // trigger type
    TString   fTrigName; // trigger name (if it's not available in AliVEvent)
    Double_t  fV0VertexerSels[7];        // Array to store the 7 values for the different selections V0 related
    Double_t  fCascadeVertexerSels[8];   // Array to store the 8 values for the different selections Casc. related

//===========================================================================================
//   Variables for Event Tree
//===========================================================================================
    Float_t fAmplitude_V0A;   //!
    Float_t fAmplitude_V0C;   //!
    Float_t fAmplitude_V0M;   //!
    Float_t fAmplitude_V0Apartial;   //!
    Float_t fAmplitude_V0Cpartial;   //!
    Float_t fAmplitude_V0AEq; //!
    Float_t fAmplitude_V0CEq; //!
    Float_t fAmplitude_V0MEq; //!
    Float_t fCentrality_V0A;         //!
    Float_t fCentrality_V0C;         //!
    Float_t fCentrality_V0M;         //!
    Float_t fCentrality_OnlineV0A;         //!
    Float_t fCentrality_OnlineV0C;         //!
    Float_t fCentrality_OnlineV0M;         //!
    Float_t fCentrality_ADA;         //!
    Float_t fCentrality_ADC;         //!
    Float_t fCentrality_ADM;         //!
    Float_t fCentrality_V0AEq;       //!
    Float_t fCentrality_V0CEq;       //!
    Float_t fCentrality_V0MEq;       //!
    Float_t fCentrality_V0B;         //!
    Float_t fCentrality_V0Apartial;  //!
    Float_t fCentrality_V0Cpartial;  //!
    Float_t fCentrality_V0S;         //!
    Float_t fCentrality_V0SB;        //!
    Int_t fRefMultEta5;              //!
    Int_t fRefMultEta8;              //!
    Int_t fRunNumber;                //!
    Int_t fBunchCrossNumber;         //!

    //Differential reference multiplicity
    Int_t  fRefMultDiffEta[20]; //!

    //Not to be filled completely: tracklet eta
    Float_t fEvent_TrackletEta[500]; //!



    //Event Characterization Variables - optional
    Bool_t fEvSel_HasAtLeastSPDVertex;      //!
    Bool_t fEvSel_VtxZCut;                  //!
    Bool_t fEvSel_IsNotPileup;              //!
    Bool_t fEvSel_IsNotPileupMV;            //!
    Bool_t fEvSel_IsNotPileupInMultBins;    //!
    Bool_t fEvSel_Triggered;                //!
    Bool_t fEvSel_INELgtZERO;               //!
    Bool_t fEvSel_INELgtZEROtracklets;      //!

    Bool_t fEvSel_INELgtZERORefMult;           //!
    Bool_t fEvSel_INELgtZERORefMultTracklets;  //!
 
    //Other Selections: more dedicated filtering to be studied!
    Int_t   fEvSel_nTracklets;              //!
    Int_t   fEvSel_nTrackletsEta10; //!
    Int_t   fEvSel_nSPDClusters;            //!
    Int_t   fEvSel_nTPCtracks;            //!
    Float_t fEvSel_VtxZ; //! the actual value
    Int_t   fEvSel_nSPDPrimVertices; //! pileup vertices
    Float_t fEvSel_distZ; //! distance between largest vertices
    Int_t   fEvSel_nContributors; //!
    Int_t   fEvSel_nContributorsPileup; //!

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

    Float_t fTreeVariableNSigmasPosProton; //!
    Float_t fTreeVariableNSigmasPosPion; //!
    Float_t fTreeVariableNSigmasNegProton; //!
    Float_t fTreeVariableNSigmasNegPion; //!

    Float_t fTreeVariableDistOverTotMom;//!
    Int_t   fTreeVariableLeastNbrCrossedRows;//!
    Float_t fTreeVariableLeastRatioCrossedRowsOverFindable;//!

    //Event Multiplicity Variables
    Float_t fTreeVariableCentV0M;    //!
    Float_t fTreeVariableCentV0A;    //!
    Float_t fTreeVariableCentV0C;    //!
    Float_t fTreeVariableCentOnlineV0M;    //!
    Float_t fTreeVariableCentOnlineV0A;    //!
    Float_t fTreeVariableCentOnlineV0C;    //!
    Float_t fTreeVariableCentADM;    //!
    Float_t fTreeVariableCentADA;    //!
    Float_t fTreeVariableCentADC;    //!
    
    Float_t fTreeVariableCentV0MEq;  //!
    Float_t fTreeVariableCentV0AEq;  //!
    Float_t fTreeVariableCentV0CEq;  //!
    Float_t fTreeVariableCentV0B;  //!
    Float_t fTreeVariableCentV0Apartial;  //!
    Float_t fTreeVariableCentV0Cpartial;  //!
    Float_t fTreeVariableCentV0S;  //!
    Float_t fTreeVariableCentV0SB;  //!
    Int_t   fTreeVariableRefMultEta8;  //!
    Int_t   fTreeVariableRefMultEta5;  //!
    Int_t   fTreeVariableRunNumber; //! //want to re-quantile per run? here's your ticket
    Int_t   fTreeVariableBunchCrossNumber; //! bunch cross number to select candidates from trains / isolated BC

    //Differential reference multiplicity
    Int_t  fTreeVariableRefMultDiffEta[20]; //!

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
    Float_t fTreeCascVarCascRadius;                   //!
    Float_t fTreeCascVarV0Mass;                       //!
    Float_t fTreeCascVarV0CosPointingAngle;           //!
    Float_t fTreeCascVarV0CosPointingAngleSpecial;    //!
    Float_t fTreeCascVarV0Radius;                     //!
    Int_t   fTreeCascVarLeastNbrClusters;             //!
    Float_t fTreeCascVarDistOverTotMom;               //!

    //TPC dEdx
    Float_t fTreeCascVarNegNSigmaPion;   //!
    Float_t fTreeCascVarNegNSigmaProton; //!
    Float_t fTreeCascVarPosNSigmaPion;   //!
    Float_t fTreeCascVarPosNSigmaProton; //!
    Float_t fTreeCascVarBachNSigmaPion;  //!
    Float_t fTreeCascVarBachNSigmaKaon;  //!

    //Event Multiplicity Variables
    Float_t fTreeCascVarCentV0M;    //!
    Float_t fTreeCascVarCentV0A;    //!
    Float_t fTreeCascVarCentV0C;    //!
    Float_t fTreeCascVarCentOnlineV0M;    //!
    Float_t fTreeCascVarCentOnlineV0A;    //!
    Float_t fTreeCascVarCentOnlineV0C;    //!
    Float_t fTreeCascVarCentADM;    //!
    Float_t fTreeCascVarCentADA;    //!
    Float_t fTreeCascVarCentADC;    //!
    Float_t fTreeCascVarCentV0MEq;  //!
    Float_t fTreeCascVarCentV0AEq;  //!
    Float_t fTreeCascVarCentV0CEq;  //!
    Float_t fTreeCascVarCentV0B;  //!
    Float_t fTreeCascVarCentV0Apartial;  //!
    Float_t fTreeCascVarCentV0Cpartial;  //!
    Float_t fTreeCascVarCentV0S;  //!
    Float_t fTreeCascVarCentV0SB;  //!
    Int_t fTreeCascVarRefMultEta8;  //!
    Int_t fTreeCascVarRefMultEta5;  //!
    Int_t fTreeCascVarRunNumber;    //! //want to re-quantile per run? here's your ticket
    Int_t fTreeCascVarBunchCrossNumber; //! bunch cross number to select candidates from trains / isolated BC

    //Differential reference multiplicity
    Int_t  fTreeCascVarRefMultDiffEta[20]; //!



//===========================================================================================
//   Histograms
//===========================================================================================

    TH1D *fHistEventCounter; //!

    AliAnalysisTaskStrangenessVsMultiplicity(const AliAnalysisTaskStrangenessVsMultiplicity&);            // not implemented
    AliAnalysisTaskStrangenessVsMultiplicity& operator=(const AliAnalysisTaskStrangenessVsMultiplicity&); // not implemented

    ClassDef(AliAnalysisTaskStrangenessVsMultiplicity, 11);
};

#endif
