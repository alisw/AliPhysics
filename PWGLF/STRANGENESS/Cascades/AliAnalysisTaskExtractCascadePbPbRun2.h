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

#ifndef ALIANALYSISTASKEXTRACTCASCADEPBPBRUN2_H
#define ALIANALYSISTASKEXTRACTCASCADEPBPBRUN2_H

class TList;
class TH1F;
class TH2F;
class TH3F;
class TVector3;
class THnSparse;

class AliESDpid;
class AliESDtrackCuts;
class AliESDEvent;
class AliPhysicsSelection;
class AliCFContainer;
class AliAnalysisUtils;

//#include "TString.h"
//#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskExtractCascadePbPbRun2 : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskExtractCascadePbPbRun2();
    AliAnalysisTaskExtractCascadePbPbRun2(const char *name, TString lExtraOptions="dEdx");
    virtual ~AliAnalysisTaskExtractCascadePbPbRun2();
    
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    Double_t MyRapidity(Double_t rE, Double_t rPz) const;
    
    //---------------------------------------------------------------------------------------
    //Task Configuration: Meant to enable quick re-execution of vertexer if needed
    void SetRunVertexers ( Bool_t lRunVertexers = kTRUE) { fkRunVertexers = lRunVertexers; }
    //---------------------------------------------------------------------------------------
    //Set Peripheral event debugging mode (Pb-Pb X-check)
    void SetSelectCentrality ( Bool_t lSelectCentrality = kTRUE, Double_t lCentSelLow = 0.0, Double_t lCentSelHigh = 10.0) {
        fkSelectCentrality = lSelectCentrality;
        fCentSel_Low = lCentSelLow;
        fCentSel_High = lCentSelHigh;
    }
    //---------------------------------------------------------------------------------------
    void SetLowPtCutoff ( Double_t lLowPtCutoff = 1.0) {
        fLowPtCutoff = lLowPtCutoff;
    }
    //---------------------------------------------------------------------------------------
    void SetCascadeMassWindow ( Double_t lCascadeMassWindow = 0.060) {
        fCascadeMassWindow = lCascadeMassWindow;
    }
    //---------------------------------------------------------------------------------------
    //Setters for the V0 Vertexer Parameters
    void SetV0VertexerMaxChisquare   ( Double_t lParameter ){ fV0VertexerSels[0] = lParameter; }
    void SetV0VertexerDCAFirstToPV   ( Double_t lParameter ){ fV0VertexerSels[1] = lParameter; }
    void SetV0VertexerDCASecondtoPV  ( Double_t lParameter ){ fV0VertexerSels[2] = lParameter; }
    void SetV0VertexerDCAV0Daughters ( Double_t lParameter ){ fV0VertexerSels[3] = lParameter; }
    void SetV0VertexerCosinePA       ( Double_t lParameter ){ fV0VertexerSels[4] = lParameter; }
    void SetV0VertexerMinRadius      ( Double_t lParameter ){ fV0VertexerSels[5] = lParameter; }
    void SetV0VertexerMaxRadius      ( Double_t lParameter ){ fV0VertexerSels[6] = lParameter; }
    //---------------------------------------------------------------------------------------
    //Setters for the Cascade Vertexer Parameters
    void SetCascVertexerMaxChisquare         ( Double_t lParameter ){ fCascadeVertexerSels[0] = lParameter; }
    void SetCascVertexerMinV0ImpactParameter ( Double_t lParameter ){ fCascadeVertexerSels[1] = lParameter; }
    void SetCascVertexerV0MassWindow         ( Double_t lParameter ){ fCascadeVertexerSels[2] = lParameter; }
    void SetCascVertexerDCABachToPV          ( Double_t lParameter ){ fCascadeVertexerSels[3] = lParameter; }
    void SetCascVertexerDCACascadeDaughters  ( Double_t lParameter ){ fCascadeVertexerSels[4] = lParameter; }
    void SetCascVertexerCascadeCosinePA      ( Double_t lParameter ){ fCascadeVertexerSels[5] = lParameter; }
    void SetCascVertexerCascadeMinRadius     ( Double_t lParameter ){ fCascadeVertexerSels[6] = lParameter; }
    void SetCascVertexerCascadeMaxRadius     ( Double_t lParameter ){ fCascadeVertexerSels[7] = lParameter; }
    //---------------------------------------------------------------------------------------
    
    //---------------------------------------------------------------------------------------
private:
    // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
    // your data member object is created on the worker nodes and streaming is not needed.
    // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14
    TList  *fListHist;  //! List of Cascade histograms
    TTree  *fTreeCascade;              //! Output Tree, Cascades
    
    //Objects that have to be streamed:
    AliPIDResponse *fPIDResponse;     // PID response object
    AliESDtrackCuts *fESDtrackCuts;   // ESD track cuts used for primary track definition
    AliAnalysisUtils *fUtils;         // analysis utils (for pA vertex selection)
    
    //Objects Controlling Task Behaviour
    // (have to be streamed too or configuration is lost)
    //Objects Controlling Task Behaviour: has to be streamed!
    Bool_t    fkRunVertexers;           // if true, re-run vertexer with loose cuts. CARE MUST BE TAKEN in PbPb!
    Double_t  fV0VertexerSels[7];        // Array to store the 7 values for the different selections V0 related
    Double_t  fCascadeVertexerSels[8];   // Array to store the 8 values for the different selections Casc. related
    
    Bool_t fkSaveTree; // if true, saves TTree object
    Bool_t fkSaveRawdEdxSignals; // if true, will save raw dEdx signals for later use
    
    Bool_t fkSelectCentrality; // if true, perform cut on centrality 
    Double_t fCentSel_Low;
    Double_t fCentSel_High;
    
    Double_t fLowPtCutoff;          //Reduction of data volume
    Double_t fCascadeMassWindow;    //Reduction of data volume: mass window selection
    
    //===========================================================================================
    //   Variables for tree, cascades
    //===========================================================================================
    
    Int_t fTreeCascVarCharge;         //!
    Float_t fTreeCascVarMassAsXi;     //!
    Float_t fTreeCascVarMassAsOmega;  //!
    Float_t fTreeCascVarPt;           //!
    Float_t fTreeCascVarPtMC;         //!
    Float_t fTreeCascVarRapMC;        //!
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
    Float_t fTreeCascVarCascCosPointingAngle; //!
    Float_t fTreeCascVarCascRadius;           //!
    Float_t fTreeCascVarV0Mass;               //!
    Float_t fTreeCascVarV0CosPointingAngle;   //!
    Float_t fTreeCascVarV0CosPointingAngleSpecial;   //!
    Float_t fTreeCascVarV0Radius;             //!
    Int_t   fTreeCascVarLeastNbrClusters;     //!
    Float_t fTreeCascVarCentrality;         //!
    Float_t fTreeCascVarDistOverTotMom;       //!

    Float_t fTreeCascVarNegNSigmaPion;   //!
    Float_t fTreeCascVarNegNSigmaProton; //!
    Float_t fTreeCascVarPosNSigmaPion;   //!
    Float_t fTreeCascVarPosNSigmaProton; //!
    Float_t fTreeCascVarBachNSigmaPion;  //!
    Float_t fTreeCascVarBachNSigmaKaon;  //!

	Bool_t fTreeCascVarNegInDistortedRegion; //!
	Bool_t fTreeCascVarPosInDistortedRegion; //!
	Bool_t fTreeCascVarBachInDistortedRegion; //!
	
	Float_t fTreeCascVarNegInnerP; //!
	Float_t fTreeCascVarPosInnerP; //!
    Float_t fTreeCascVarBachInnerP; //!

	Float_t fTreeCascVarNegdEdx; //!
	Float_t fTreeCascVarPosdEdx; //!
    Float_t fTreeCascVarBachdEdx; //!
    
    //Debugging information, if requested
    //Part A: EbyE info, Run number
    Int_t     fTreeCascVarRunNumber; //!
    
    //===========================================================================================
    //   Histograms
    //===========================================================================================
    
    TH1F      *fHistEventCounter; //! Histogram with basic event counting 
    TH1F      *fHistCentrality;   //! Histogram with Centrality Distribution 
    
    //For TPC dEdx reparametrization, if desired 
    TH2F *fHistdEdx;                   //! inclusive dEdx plot (all tracks passing some set of standard selections) 
    TH2F *fHistdEdxPionsFromLambda;    //! Histogram filled only with what looks like pions from Lambda decays
    TH2F *fHistdEdxProtonsFromLambda;  //! Histogram filled only with what looks like protons from Lambda decays
    TH2F *fHistdEdxPionsFromK0s;       //! Histogram filled only with what looks like pions from K0s decays
    
        
    //=======================================================
    //     --- Superlight Output Mode - Experimental ---
    //=======================================================
        
    AliAnalysisTaskExtractCascadePbPbRun2(const AliAnalysisTaskExtractCascadePbPbRun2&);            // not implemented
    AliAnalysisTaskExtractCascadePbPbRun2& operator=(const AliAnalysisTaskExtractCascadePbPbRun2&); // not implemented
    
    ClassDef(AliAnalysisTaskExtractCascadePbPbRun2, 11);
};

#endif

