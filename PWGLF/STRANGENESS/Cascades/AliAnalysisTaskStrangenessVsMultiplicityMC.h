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

#ifndef AliAnalysisTaskStrangenessVsMultiplicityMC_H
#define AliAnalysisTaskStrangenessVsMultiplicityMC_H

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

class AliAnalysisTaskStrangenessVsMultiplicityMC : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskStrangenessVsMultiplicityMC();
    AliAnalysisTaskStrangenessVsMultiplicityMC(const char *name);
    virtual ~AliAnalysisTaskStrangenessVsMultiplicityMC();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    Double_t MyRapidity(Double_t rE, Double_t rPz) const;
    Bool_t   IsINELgtZERO(AliESDEvent *lESDevent, TString lType) const;
    //Bool_t   IsINELgtZEROPrimOnly  (AliESDEvent *lESDevent, AliStack *lStack) const;
    //Bool_t   IsINELgtZEROGenerator (AliStack *lStack) const;

    void SetSaveV0s                (Bool_t lSaveV0s        = kTRUE ) {
        fkSaveV0Tree        = lSaveV0s;
    }
    void SetSaveCascades           (Bool_t lSaveCascades   = kTRUE ) {
        fkSaveCascadeTree   = lSaveCascades;
    }
    void SetMonteCarloAssociation  (Bool_t lMCAssociation   = kTRUE ) {
        fkMCAssociation   = lMCAssociation;
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
    AliAnalysisUtils *fUtils; //

    //Objects Controlling Task Behaviour
    Bool_t fkSaveV0Tree;              //if true, save TTree
    Bool_t fkSaveCascadeTree;         //if true, save TTree
    Bool_t fkMCAssociation;           //if true, save only MC associated particles
    Bool_t fkSaveLambda;            //if true, removes tree fill confirmation for lambda mass window
    Bool_t fkSaveAntiLambda;        //if true, removes tree fill confirmation for antilambda mass window
    Bool_t fkSaveK0Short;           //if true, removes tree fill confirmation for k0short mass window 
    Bool_t fkSaveExtendedRefMultInfo; //if true, save 20 integers per event extra for eta-dependence of refmult
    Bool_t fkSelectTriggerByName;  // to select trigger by name (if it's not availble in AliVEvent)
    Bool_t fkMultSelection; //if true, use new framework (under tests) 
 
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
    Int_t fTrueMultEta5;             //!
    Int_t fTrueMultEta8;             //!
    Int_t fTrueMultEta10;             //!
    Int_t fTrueMultVZEROA;           //!
    Int_t fTrueMultVZEROC;           //!
    Int_t fRunNumber;                //!
    
    Int_t fEvSel_nTrackletsEta10; //!

    //Differential reference multiplicity
    Int_t  fRefMultDiffEta[20]; //!

    //Event Characterization Variables - optional
    Bool_t fEvSel_HasAtLeastSPDVertex;      //!
    Bool_t fEvSel_VtxZCut;                  //!
    Bool_t fEvSel_IsNotPileup;              //!
    Bool_t fEvSel_IsNotPileupMV;            //!
    Bool_t fEvSel_IsNotPileupInMultBins;    //!
    Bool_t fEvSel_HasVtxContributor;        //!
    Bool_t fEvSel_Triggered;                //!

    Bool_t fEvSel_INELgtZERO;               //!
    Bool_t fEvSel_INELgtZEROStackPrimaries; //!
    Bool_t fEvSel_INELgtZEROtracklets;      //!

    Bool_t fEvSel_INELgtZERORefMult;           //!
    Bool_t fEvSel_INELgtZERORefMultTracklets;  //!

    Float_t fEvSel_VtxZ; //! pv z position (cm)
    Float_t fEvSel_VtxZMC; //! pv z pos from mc record
    Int_t fEvSel_MCType; //! type of event (to be used in PYTHIA, specifically)

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
    Float_t fTreeVariableNegTransvMomentumMC; //!
    Float_t fTreeVariablePosTransvMomentumMC; //!

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

    //Differential reference multiplicity
    Int_t  fTreeVariableRefMultDiffEta[20]; //!

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
    Int_t   fTreeCascVarLeastNbrClusters;             //!
    Float_t fTreeCascVarDistOverTotMom;               //!

    //TPC dEdx
    Float_t fTreeCascVarNegNSigmaPion;   //!
    Float_t fTreeCascVarNegNSigmaProton; //!
    Float_t fTreeCascVarPosNSigmaPion;   //!
    Float_t fTreeCascVarPosNSigmaProton; //!
    Float_t fTreeCascVarBachNSigmaPion;  //!
    Float_t fTreeCascVarBachNSigmaKaon;  //!
    Float_t fTreeCascVarNegTransvMomentumMC; //!
    Float_t fTreeCascVarPosTransvMomentumMC; //!

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
    Int_t fTreeCascVarTrueMultEta5; //!
    Int_t fTreeCascVarTrueMultEta8; //!
    Int_t fTreeCascVarTrueMultVZEROA; //!
    Int_t fTreeCascVarTrueMultVZEROC; //!

    //Differential reference multiplicity
    Int_t  fTreeCascVarRefMultDiffEta[20]; //!

    //MC-only Variabless
    Int_t   fTreeCascVarIsPhysicalPrimary; //!
    Int_t   fTreeCascVarPID;         //!

    Int_t   fTreeCascVarRunNumber;         //!

//===========================================================================================
//   Histograms
//===========================================================================================

    TH1D *fHistEventCounter; //!

//===========================================================================================
//  Histos needed for efficiency computation || Xis
//===========================================================================================

    //These are all analysis-level

    //Let's not care too much about rapidity at the moment!
    TH1D *fHistPt_GenK0Short;      //!
    TH1D *fHistPt_GenLambda;       //!
    TH1D *fHistPt_GenAntiLambda;   //!
    TH1D *fHistPt_GenXiMinus;      //!
    TH1D *fHistPt_GenXiPlus;       //!
    TH1D *fHistPt_GenOmegaMinus;   //!
    TH1D *fHistPt_GenOmegaPlus;    //!

    //VsRefMult
    TH2D *fHistPtVsRefMultEta5_GenK0Short;      //!
    TH2D *fHistPtVsRefMultEta5_GenLambda;       //!
    TH2D *fHistPtVsRefMultEta5_GenAntiLambda;   //!
    TH2D *fHistPtVsRefMultEta5_GenXiMinus;      //!
    TH2D *fHistPtVsRefMultEta5_GenXiPlus;       //!
    TH2D *fHistPtVsRefMultEta5_GenOmegaMinus;   //!
    TH2D *fHistPtVsRefMultEta5_GenOmegaPlus;    //!
    TH2D *fHistPtVsRefMultEta8_GenK0Short;      //!
    TH2D *fHistPtVsRefMultEta8_GenLambda;       //!
    TH2D *fHistPtVsRefMultEta8_GenAntiLambda;   //!
    TH2D *fHistPtVsRefMultEta8_GenXiMinus;      //!
    TH2D *fHistPtVsRefMultEta8_GenXiPlus;       //!
    TH2D *fHistPtVsRefMultEta8_GenOmegaMinus;   //!
    TH2D *fHistPtVsRefMultEta8_GenOmegaPlus;    //!

    //VsCentralities
    TH2D *fHistPtVsCentV0A_GenK0Short;      //!
    TH2D *fHistPtVsCentV0A_GenLambda;       //!
    TH2D *fHistPtVsCentV0A_GenAntiLambda;   //!
    TH2D *fHistPtVsCentV0A_GenXiMinus;      //!
    TH2D *fHistPtVsCentV0A_GenXiPlus;       //!
    TH2D *fHistPtVsCentV0A_GenOmegaMinus;   //!
    TH2D *fHistPtVsCentV0A_GenOmegaPlus;    //!
    TH2D *fHistPtVsCentV0C_GenK0Short;      //!
    TH2D *fHistPtVsCentV0C_GenLambda;       //!
    TH2D *fHistPtVsCentV0C_GenAntiLambda;   //!
    TH2D *fHistPtVsCentV0C_GenXiMinus;      //!
    TH2D *fHistPtVsCentV0C_GenXiPlus;       //!
    TH2D *fHistPtVsCentV0C_GenOmegaMinus;   //!
    TH2D *fHistPtVsCentV0C_GenOmegaPlus;    //!
    TH2D *fHistPtVsCentV0M_GenK0Short;      //!
    TH2D *fHistPtVsCentV0M_GenLambda;       //!
    TH2D *fHistPtVsCentV0M_GenAntiLambda;   //!
    TH2D *fHistPtVsCentV0M_GenXiMinus;      //!
    TH2D *fHistPtVsCentV0M_GenXiPlus;       //!
    TH2D *fHistPtVsCentV0M_GenOmegaMinus;   //!
    TH2D *fHistPtVsCentV0M_GenOmegaPlus;    //!
    
    //Pre-selection stage
    TH2D *fHistPtVsCentV0MUnselected_GenK0Short;      //!
    TH2D *fHistPtVsCentV0MUnselected_GenLambda;       //!
    TH2D *fHistPtVsCentV0MUnselected_GenAntiLambda;   //!
    TH2D *fHistPtVsCentV0MUnselected_GenXiMinus;      //!
    TH2D *fHistPtVsCentV0MUnselected_GenXiPlus;       //!
    TH2D *fHistPtVsCentV0MUnselected_GenOmegaMinus;   //!
    TH2D *fHistPtVsCentV0MUnselected_GenOmegaPlus;    //!  

    //Equalized
    TH2D *fHistPtVsCentV0AEq_GenXiMinus;      //!
    TH2D *fHistPtVsCentV0AEq_GenXiPlus;       //!
    TH2D *fHistPtVsCentV0AEq_GenOmegaMinus;   //!
    TH2D *fHistPtVsCentV0AEq_GenOmegaPlus;    //!
    TH2D *fHistPtVsCentV0CEq_GenXiMinus;      //!
    TH2D *fHistPtVsCentV0CEq_GenXiPlus;       //!
    TH2D *fHistPtVsCentV0CEq_GenOmegaMinus;   //!
    TH2D *fHistPtVsCentV0CEq_GenOmegaPlus;    //!
    TH2D *fHistPtVsCentV0MEq_GenXiMinus;      //!
    TH2D *fHistPtVsCentV0MEq_GenXiPlus;       //!
    TH2D *fHistPtVsCentV0MEq_GenOmegaMinus;   //!
    TH2D *fHistPtVsCentV0MEq_GenOmegaPlus;    //!

    //VsAmp
    TH2D *fHistPtVsAmpV0A_GenXiMinus;      //!
    TH2D *fHistPtVsAmpV0A_GenXiPlus;       //!
    TH2D *fHistPtVsAmpV0A_GenOmegaMinus;   //!
    TH2D *fHistPtVsAmpV0A_GenOmegaPlus;    //!
    TH2D *fHistPtVsAmpV0C_GenXiMinus;      //!
    TH2D *fHistPtVsAmpV0C_GenXiPlus;       //!
    TH2D *fHistPtVsAmpV0C_GenOmegaMinus;   //!
    TH2D *fHistPtVsAmpV0C_GenOmegaPlus;    //!
    TH2D *fHistPtVsAmpV0M_GenXiMinus;      //!
    TH2D *fHistPtVsAmpV0M_GenXiPlus;       //!
    TH2D *fHistPtVsAmpV0M_GenOmegaMinus;   //!
    TH2D *fHistPtVsAmpV0M_GenOmegaPlus;    //!
    //Equalized Amps
    TH2D *fHistPtVsAmpV0AEq_GenXiMinus;      //!
    TH2D *fHistPtVsAmpV0AEq_GenXiPlus;       //!
    TH2D *fHistPtVsAmpV0AEq_GenOmegaMinus;   //!
    TH2D *fHistPtVsAmpV0AEq_GenOmegaPlus;    //!
    TH2D *fHistPtVsAmpV0CEq_GenXiMinus;      //!
    TH2D *fHistPtVsAmpV0CEq_GenXiPlus;       //!
    TH2D *fHistPtVsAmpV0CEq_GenOmegaMinus;   //!
    TH2D *fHistPtVsAmpV0CEq_GenOmegaPlus;    //!
    TH2D *fHistPtVsAmpV0MEq_GenXiMinus;      //!
    TH2D *fHistPtVsAmpV0MEq_GenXiPlus;       //!
    TH2D *fHistPtVsAmpV0MEq_GenOmegaMinus;   //!
    TH2D *fHistPtVsAmpV0MEq_GenOmegaPlus;    //!

    TH2D *fHistVZEROResponseStudy; //!
    TH2D *fHistVZEROResponseStudyTotal; //!


    AliAnalysisTaskStrangenessVsMultiplicityMC(const AliAnalysisTaskStrangenessVsMultiplicityMC&);            // not implemented
    AliAnalysisTaskStrangenessVsMultiplicityMC& operator=(const AliAnalysisTaskStrangenessVsMultiplicityMC&); // not implemented

    ClassDef(AliAnalysisTaskStrangenessVsMultiplicityMC, 11);
};

#endif
