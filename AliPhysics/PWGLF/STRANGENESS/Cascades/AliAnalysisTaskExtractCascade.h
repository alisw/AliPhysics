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

#ifndef ALIANALYSISTASKEXTRACTCASCADE_H
#define ALIANALYSISTASKEXTRACTCASCADE_H

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

class AliAnalysisTaskExtractCascade : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskExtractCascade();
    AliAnalysisTaskExtractCascade(const char *name);
    virtual ~AliAnalysisTaskExtractCascade();
    
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    Double_t MyRapidity(Double_t rE, Double_t rPz) const;
    
    void SetIsNuclear           (Bool_t lIsNuclear   = kTRUE ) { fkIsNuclear   = lIsNuclear;   }
    void SetINT7Trigger         (Bool_t lSwitchINT7  = kTRUE ) { fkSwitchINT7   = lSwitchINT7; }
    void SetCentralityEstimator (TString lCentralityEstimator = "V0M" ) { fCentralityEstimator = lCentralityEstimator; }
    void SetpAVertexSelection   (Bool_t lpAVertexSelection = kTRUE) {fkpAVertexSelection = lpAVertexSelection;  }
    void SetEtaRefMult ( Double_t lEtaRefMult = 0.5 ) { fEtaRefMult = lEtaRefMult; }
    
    //---------------------------------------------------------------------------------------
    //Task Configuration: Meant to enable quick re-execution of vertexer if needed
    void SetRunVertexers ( Bool_t lRunVertexers = kTRUE) { fkRunVertexers = lRunVertexers; }
    //---------------------------------------------------------------------------------------
    //Set Debug Mode
    void SetDebugMode ( Bool_t lDebugMode = kTRUE) { fkDebugMode = lDebugMode; }
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
    void SetSuperLightMode ( Bool_t lSuperLight = kTRUE) {
        fkSuperLight = lSuperLight;
    }
    //---------------------------------------------------------------------------------------
    void SetLightMode ( Bool_t lLight = kTRUE) {
        fkLight = lLight;
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
    //Setters for SuperLight analysis
    void SetSuperLightV0Radius               ( Double_t lParameter ){ fCut_V0Radius   = lParameter;         }
    void SetSuperLightCascRadius             ( Double_t lParameter ){ fCut_CascRadius = lParameter;         }
    void SetSuperLightV0Mass                 ( Double_t lParameter ){ fCut_V0Mass     = lParameter;         }
    void SetSuperLightV0CosPA                ( Double_t lParameter ){ fCut_V0CosPA    = lParameter;         }
    void SetSuperLightCascCosPA              ( Double_t lParameter ){ fCut_CascCosPA  = lParameter;         }
    void SetSuperLightDCANegToPV             ( Double_t lParameter ){ fCut_DCANegToPV = lParameter;         }
    void SetSuperLightDCAPosToPV             ( Double_t lParameter ){ fCut_DCAPosToPV = lParameter;         }
    void SetSuperLightDCABachToPV            ( Double_t lParameter ){ fCut_DCABachToPV = lParameter;        }
    void SetSuperLightDCAV0Daughters         ( Double_t lParameter ){ fCut_DCAV0Daughters = lParameter;     }
    void SetSuperLightDCACascDaughters       ( Double_t lParameter ){ fCut_DCACascDaughters = lParameter;   }
    void SetSuperLightDCAV0ToPV              ( Double_t lParameter ){ fCut_DCAV0ToPV = lParameter;          }
    void SetSuperLightCTau                   ( Double_t lParameter ){ fCut_CTau = lParameter;               }
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
    Bool_t fkIsNuclear;   //if true, replace multiplicity est. by centrality (default FALSE)
    Bool_t fkSwitchINT7;  //if true, skip FASTOnly (default FALSE)
    TString fCentralityEstimator; //Centrality Estimator String value (default V0M)
    Bool_t fkpAVertexSelection; //if true, select vertex with pPb Methods
    Double_t fEtaRefMult; //Reference multiplicity eta
    //Objects Controlling Task Behaviour: has to be streamed!
    Bool_t    fkRunVertexers;           // if true, re-run vertexer with loose cuts. CARE MUST BE TAKEN in PbPb!
    Bool_t    fkDebugMode; //Store extra debugging information in tree
    Bool_t    fkSelectCentrality; //Switch to skip anything other than 60-80% V0M
    Double_t  fV0VertexerSels[7];        // Array to store the 7 values for the different selections V0 related
    Double_t  fCascadeVertexerSels[8];   // Array to store the 8 values for the different selections Casc. related
    
    Double_t fCentSel_Low;
    Double_t fCentSel_High;
    
    Double_t fLowPtCutoff;          //Reduction of data volume
    Double_t fCascadeMassWindow;    //Reduction of data volume: mass window selection
    
	//Double_t        fV0Sels[7];                     // Array to store the 7 values for the different selections V0 related
	//Double_t        fCascSels[8];                   // Array to store the 8 values for the different selections Casc. related
    
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
    Int_t   fTreeCascVarMultiplicity;         //!
    Int_t   fTreeCascVarMultiplicityV0A;         //!
    Int_t   fTreeCascVarMultiplicityZNA;         //!
    Int_t   fTreeCascVarMultiplicityTRK;         //!
    Int_t   fTreeCascVarMultiplicitySPD;         //!
    Float_t fTreeCascVarDistOverTotMom;       //!
    Int_t   fTreeCascVarPID;         //!
    Int_t   fTreeCascVarPIDBachelor; //!
    Int_t   fTreeCascVarPIDNegative; //!
    Int_t   fTreeCascVarPIDPositive; //!
    Float_t fTreeCascVarBachTransMom;   //!
    Float_t fTreeCascVarPosTransMom;   //!
    Float_t fTreeCascVarNegTransMom;   //!
    Float_t fTreeCascVarPosTransMomMC; //!
    Float_t fTreeCascVarNegTransMomMC; //!
    
    Float_t fTreeCascVarNegNSigmaPion;   //!
    Float_t fTreeCascVarNegNSigmaProton; //!
    Float_t fTreeCascVarPosNSigmaPion;   //!
    Float_t fTreeCascVarPosNSigmaProton; //!
    Float_t fTreeCascVarBachNSigmaPion;  //!
    Float_t fTreeCascVarBachNSigmaKaon;  //!
    
    Bool_t fTreeCascVarkITSRefitBachelor; //!
    Bool_t fTreeCascVarkITSRefitNegative; //!
    Bool_t fTreeCascVarkITSRefitPositive; //!
    
    //Debugging information, if requested
    //Part A: EbyE info, Run number
    Int_t     fTreeCascVarRunNumber; //!
    ULong64_t fTreeCascVarEventNumber; //!
    
    //Part B: Shared Clusters
    Int_t fTreeCascVarNegClusters; //!
    Int_t fTreeCascVarPosClusters; //!
    Int_t fTreeCascVarBachClusters; //!
    Int_t fTreeCascVarNegSharedClusters; //!
    Int_t fTreeCascVarPosSharedClusters; //!
    Int_t fTreeCascVarBachSharedClusters; //!
    
    //Part C: All momenta
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
    
    Bool_t fTreeCascVarBadCascadeJai; //!
    Float_t fTreeCascVarDeltaDCA; //!
    
    //===========================================================================================
    //   Histograms
    //===========================================================================================
    
    TH1F      *fHistV0MultiplicityBeforeTrigSel;              //! V0 multiplicity distribution
    TH1F      *fHistV0MultiplicityForTrigEvt;                 //! V0 multiplicity distribution
    TH1F      *fHistV0MultiplicityForSelEvt;                  //! V0 multiplicity distribution
    TH1F      *fHistV0MultiplicityForSelEvtNoTPCOnly;         //! V0 multiplicity distribution
    TH1F      *fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup; //! V0 multiplicity distribution
    
    //V0M Centrality (default)
    TH1F      *fHistMultiplicityBeforeTrigSel;     //! multiplicity distribution
    TH1F      *fHistMultiplicityForTrigEvt;        //! multiplicity distribution
    TH1F      *fHistMultiplicity;                  //! multiplicity distribution
    TH1F      *fHistMultiplicityNoTPCOnly;         //! multiplicity distribution
    TH1F      *fHistMultiplicityNoTPCOnlyNoPileup; //! multiplicity distribution
    
    //V0A Centrality
    TH1F    *fHistMultiplicityV0ABeforeTrigSel; 	        //! multiplicity distribution
    TH1F    *fHistMultiplicityV0AForTrigEvt;  		        //! multiplicity distribution
    TH1F    *fHistMultiplicityV0AAfterPVSelection;     //! multiplicity distribution
    TH1F    *fHistMultiplicityV0A;     					        //! multiplicity distribution
    TH1F    *fHistMultiplicityV0ANoTPCOnly;			        //! multiplicity distribution
    TH1F    *fHistMultiplicityV0ANoTPCOnlyNoPileup;			//! multiplicity distribution
    
    //ZNA Centrality
    TH1F    *fHistMultiplicityZNABeforeTrigSel; 	        //! multiplicity distribution
	TH1F    *fHistMultiplicityZNAForTrigEvt;  		        //! multiplicity distribution
	TH1F    *fHistMultiplicityZNA;     					        //! multiplicity distribution
	TH1F    *fHistMultiplicityZNANoTPCOnly;			        //! multiplicity distribution
	TH1F    *fHistMultiplicityZNANoTPCOnlyNoPileup;			//! multiplicity distribution
    
    //TRK Centrality
    TH1F    *fHistMultiplicityTRKBeforeTrigSel; 	        //! multiplicity distribution
	TH1F    *fHistMultiplicityTRKForTrigEvt;  		        //! multiplicity distribution
	TH1F    *fHistMultiplicityTRK;     					        //! multiplicity distribution
	TH1F    *fHistMultiplicityTRKNoTPCOnly;			        //! multiplicity distribution
	TH1F    *fHistMultiplicityTRKNoTPCOnlyNoPileup;			//! multiplicity distribution
    
    //SPD Centrality
    TH1F    *fHistMultiplicitySPDBeforeTrigSel; 	        //! multiplicity distribution
	TH1F    *fHistMultiplicitySPDForTrigEvt;  		        //! multiplicity distribution
	TH1F    *fHistMultiplicitySPD;     					        //! multiplicity distribution
	TH1F    *fHistMultiplicitySPDNoTPCOnly;			        //! multiplicity distribution
	TH1F    *fHistMultiplicitySPDNoTPCOnlyNoPileup;			//! multiplicity distribution

  //
	TH1F*  fHistPVZDistribution;  // PV z-distribution
    
    //---> Generated Histograms
    
    TH1F      *fHistPVx;                      //! PVx distrib
    TH1F      *fHistPVy;                      //! PVy distrib
    TH1F      *fHistPVz;                      //! PVz distrib
    TH1F      *fHistPVxAnalysis;                      //! PVx distrib
    TH1F      *fHistPVyAnalysis;                      //! PVy distrib
    TH1F      *fHistPVzAnalysis;                      //! PVz distrib
    
    //=======================================================
    //     --- Superlight Output Mode - Experimental ---
    //=======================================================
    
    Bool_t fkLight; //Switch for intermediate mode
    Bool_t fkSuperLight; //Switch for super light output mode
    
    //Selection criteria for superlight analysis
    Double_t fCut_V0Radius;
    Double_t fCut_CascRadius;
    Double_t fCut_V0Mass;
    Double_t fCut_V0CosPA;
    Double_t fCut_CascCosPA;
    Double_t fCut_DCANegToPV;
    Double_t fCut_DCAPosToPV;
    Double_t fCut_DCABachToPV;
    Double_t fCut_DCAV0Daughters;
    Double_t fCut_DCACascDaughters;
    Double_t fCut_DCAV0ToPV;
    Double_t fCut_CTau;
    
    //---> Super-lightweight histogram output (rebinnable)
    
    TH2F *f2dHist_MassVsPt_XiMinus;
    TH2F *f2dHist_MassVsPt_XiPlus;
    TH2F *f2dHist_MassVsPt_OmegaMinus;
    TH2F *f2dHist_MassVsPt_OmegaPlus;
    
    AliAnalysisTaskExtractCascade(const AliAnalysisTaskExtractCascade&);            // not implemented
    AliAnalysisTaskExtractCascade& operator=(const AliAnalysisTaskExtractCascade&); // not implemented
    
    ClassDef(AliAnalysisTaskExtractCascade, 11);
};

#endif
