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

#ifndef ALIANALYSISTASKEXTRACTPERFORMANCECASCADE_H
#define ALIANALYSISTASKEXTRACTPERFORMANCECASCADE_H

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

class AliAnalysisTaskExtractPerformanceCascade : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskExtractPerformanceCascade();
  AliAnalysisTaskExtractPerformanceCascade(const char *name);
  virtual ~AliAnalysisTaskExtractPerformanceCascade();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  Double_t MyRapidity(Double_t rE, Double_t rPz) const;

  void SetIsNuclear           (Bool_t lIsNuclear   = kTRUE ) { fkIsNuclear   = lIsNuclear;   }
  void SetINT7Trigger         (Bool_t lSwitchINT7  = kTRUE ) { fkSwitchINT7   = lSwitchINT7; }
  void SetpARapidityShift     (Double_t lRapShift = 0.465 ) { fpArapidityShift = lRapShift; }
  void SetCentralityEstimator (TString lCentralityEstimator = "V0M" ) { fCentralityEstimator = lCentralityEstimator; }
  void SetpAVertexSelection   (Bool_t lpAVertexSelection = kTRUE) {fkpAVertexSelection = lpAVertexSelection;  }
  void SetEtaRefMult ( Double_t lEtaRefMult = 0.5 ) { fEtaRefMult = lEtaRefMult; }
  
//---------------------------------------------------------------------------------------
  //Task Configuration: Meant to enable quick re-execution of vertexer if needed
  void SetRunVertexers ( Bool_t lRunVertexers = kTRUE) { fkRunVertexers = lRunVertexers; }
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
  //Bachelor and Pion Swapping Check
  void SetCheckSwapping ( Bool_t lCheckSwapping = kTRUE) { fkCheckSwapping = lCheckSwapping; }
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
  Bool_t fkSwitchINT7 ; //if true, skip FASTOnly (default FALSE)
  Double_t fpArapidityShift; //pA rapidity shift (should be 0.465, usually)
  TString fCentralityEstimator; //Centrality Estimator String value (default V0M)
  Bool_t fkpAVertexSelection; //if true, select vertex with pPb Methods
  Double_t fEtaRefMult; //Reference multiplicity eta
  //Objects Controlling Task Behaviour: has to be streamed! 
  Bool_t    fkRunVertexers;           // if true, re-run vertexer with loose cuts. CARE MUST BE TAKEN in PbPb!
  Double_t  fV0VertexerSels[7];        // Array to store the 7 values for the different selections V0 related
  Double_t  fCascadeVertexerSels[8];   // Array to store the 8 values for the different selections Casc. related
  //Meson Swapping Check Switch
  Bool_t fkCheckSwapping; // if true, will perform association with mesons switched (in ADDITION to reg. association)
    Bool_t    fkSelectCentrality; //Switch to skip anything other than 60-80% V0M
    Double_t fCentSel_Low;
    Double_t fCentSel_High;
    Double_t fLowPtCutoff; //Reduction of data volume
    
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
  
  Int_t   fTreeCascVarMultiplicityMC;         //!
  Float_t fTreeCascVarDistOverTotMom;       //!
  Int_t   fTreeCascVarIsPhysicalPrimary; //!
  Int_t   fTreeCascVarPID;         //!
  Int_t   fTreeCascVarPIDSwapped;  //!
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
  
    //Part B: Shared Clusters
    Int_t fTreeCascVarNegClusters; //!
    Int_t fTreeCascVarPosClusters; //!
    Int_t fTreeCascVarBachClusters; //!
    Int_t fTreeCascVarNegSharedClusters; //!
    Int_t fTreeCascVarPosSharedClusters; //!
    Int_t fTreeCascVarBachSharedClusters; //!

  Bool_t fTreeCascVarEvHasXiMinus;    //!
  Bool_t fTreeCascVarEvHasXiPlus;     //!
  Bool_t fTreeCascVarEvHasOmegaMinus; //!
  Bool_t fTreeCascVarEvHasOmegaPlus;  //!
  Bool_t fTreeCascVarEvHasLambda;     //!
  Bool_t fTreeCascVarEvHasAntiLambda; //!

  Bool_t fTreeCascVarEvHasLowPtXiMinus;    //!
  Bool_t fTreeCascVarEvHasLowPtXiPlus;     //!
  Bool_t fTreeCascVarEvHasLowPtOmegaMinus; //!
  Bool_t fTreeCascVarEvHasLowPtOmegaPlus;  //!
  Bool_t fTreeCascVarEvHasLowPtLambda;     //!
  Bool_t fTreeCascVarEvHasLowPtAntiLambda; //!

  Bool_t fTreeCascVarEvHasVeryLowPtXiMinus;    //!
  Bool_t fTreeCascVarEvHasVeryLowPtXiPlus;     //!
  Bool_t fTreeCascVarEvHasVeryLowPtOmegaMinus; //!
  Bool_t fTreeCascVarEvHasVeryLowPtOmegaPlus;  //!
  Bool_t fTreeCascVarEvHasVeryLowPtLambda;     //!
  Bool_t fTreeCascVarEvHasVeryLowPtAntiLambda; //!
  
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

   TH3F      *f3dHistGenPtVsYVsMultXiMinus;      //! Generated Xi- Distrib
   TH3F      *f3dHistGenPtVsYVsMultXiPlus;       //! Generated Xi+ Distrib
   TH3F      *f3dHistGenPtVsYVsMultOmegaMinus;      //! Generated Omega- Distrib
   TH3F      *f3dHistGenPtVsYVsMultOmegaPlus;       //! Generated Omega+ Distrib

   TH3F      *f3dHistGenSelectedPtVsYVsMultXiMinus;      //! Generated Xi- Distrib, at event selection level
   TH3F      *f3dHistGenSelectedPtVsYVsMultXiPlus;       //! Generated Xi+ Distrib, at event selection level
   TH3F      *f3dHistGenSelectedPtVsYVsMultOmegaMinus;      //! Generated Omega- Distrib, at event selection level
   TH3F      *f3dHistGenSelectedPtVsYVsMultOmegaPlus;       //! Generated Omega+ Distrib, at event selection level

   TH3F      *f3dHistGenPtVsYCMSVsMultXiMinus;      //! Generated Xi- Distrib
   TH3F      *f3dHistGenPtVsYCMSVsMultXiPlus;       //! Generated Xi+ Distrib
   TH3F      *f3dHistGenPtVsYCMSVsMultOmegaMinus;      //! Generated Omega- Distrib
   TH3F      *f3dHistGenPtVsYCMSVsMultOmegaPlus;       //! Generated Omega+ Distrib

   TH3F      *f3dHistGenSelectedPtVsYCMSVsMultXiMinus;      //! Generated Xi- Distrib, at event selection level
   TH3F      *f3dHistGenSelectedPtVsYCMSVsMultXiPlus;       //! Generated Xi+ Distrib, at event selection level
   TH3F      *f3dHistGenSelectedPtVsYCMSVsMultOmegaMinus;      //! Generated Omega- Distrib, at event selection level
   TH3F      *f3dHistGenSelectedPtVsYCMSVsMultOmegaPlus;       //! Generated Omega+ Distrib, at event selection level

//---> Multiplicity -> MC multiplicity
  
  TH3F      *f3dHistGenPtVsYVsMultMCXiMinus;      //! Generated Xi- Distrib
  TH3F      *f3dHistGenPtVsYVsMultMCXiPlus;       //! Generated Xi+ Distrib
  TH3F      *f3dHistGenPtVsYVsMultMCOmegaMinus;      //! Generated Omega- Distrib
  TH3F      *f3dHistGenPtVsYVsMultMCOmegaPlus;       //! Generated Omega+ Distrib
  
  TH3F      *f3dHistGenSelectedPtVsYVsMultMCXiMinus;      //! Generated Xi- Distrib, at event selection level
  TH3F      *f3dHistGenSelectedPtVsYVsMultMCXiPlus;       //! Generated Xi+ Distrib, at event selection level
  TH3F      *f3dHistGenSelectedPtVsYVsMultMCOmegaMinus;      //! Generated Omega- Distrib, at event selection level
  TH3F      *f3dHistGenSelectedPtVsYVsMultMCOmegaPlus;       //! Generated Omega+ Distrib, at event selection level
  
  TH3F      *f3dHistGenPtVsYCMSVsMultMCXiMinus;      //! Generated Xi- Distrib
  TH3F      *f3dHistGenPtVsYCMSVsMultMCXiPlus;       //! Generated Xi+ Distrib
  TH3F      *f3dHistGenPtVsYCMSVsMultMCOmegaMinus;      //! Generated Omega- Distrib
  TH3F      *f3dHistGenPtVsYCMSVsMultMCOmegaPlus;       //! Generated Omega+ Distrib
  
  TH3F      *f3dHistGenSelectedPtVsYCMSVsMultMCXiMinus;      //! Generated Xi- Distrib, at event selection level
  TH3F      *f3dHistGenSelectedPtVsYCMSVsMultMCXiPlus;       //! Generated Xi+ Distrib, at event selection level
  TH3F      *f3dHistGenSelectedPtVsYCMSVsMultMCOmegaMinus;      //! Generated Omega- Distrib, at event selection level
  TH3F      *f3dHistGenSelectedPtVsYCMSVsMultMCOmegaPlus;       //! Generated Omega+ Distrib, at event selection level
  
    //---> Multiplicity -> V0A centrality
    
    TH3F      *f3dHistGenTriggPtVsYVsMultV0APiMinus;      //! Generated Xi- Distrib, at event selection level
    TH3F      *f3dHistGenTriggPtVsYVsMultV0APiPlus;       //! Generated Xi+ Distrib, at event selection level
    
    TH3F      *f3dHistGenPtVsYVsMultV0AXiMinus;      //! Generated Xi- Distrib
    TH3F      *f3dHistGenPtVsYVsMultV0AXiPlus;       //! Generated Xi+ Distrib
    TH3F      *f3dHistGenPtVsYVsMultV0AOmegaMinus;      //! Generated Omega- Distrib
    TH3F      *f3dHistGenPtVsYVsMultV0AOmegaPlus;       //! Generated Omega+ Distrib
    
    TH3F      *f3dHistGenTriggPtVsYCMSVsMultV0APiMinus;      //! Generated Xi- Distrib, at event selection level
    TH3F      *f3dHistGenTriggPtVsYCMSVsMultV0APiPlus;       //! Generated Xi+ Distrib, at event selection level
    
  TH3F      *f3dHistGenTriggPtVsYCMSVsMultV0AXiMinus;      //! Generated Xi- Distrib, at triggered events level
  TH3F      *f3dHistGenTriggPtVsYCMSVsMultV0AXiPlus;       //! Generated Xi+ Distrib, at triggered events level
  TH3F      *f3dHistGenTriggPtVsYCMSVsMultV0AOmegaMinus;      //! Generated Omega- Distrib, at triggered events level
  TH3F      *f3dHistGenTriggPtVsYCMSVsMultV0AOmegaPlus;       //! Generated Omega+ Distrib, at triggered events level

  TH3F      *f3dHistGenVtxPtVsYCMSVsMultV0AXiMinus;      //! Generated Xi- Distrib, at triggered events level
  TH3F      *f3dHistGenVtxPtVsYCMSVsMultV0AXiPlus;       //! Generated Xi+ Distrib, at triggered events level
  TH3F      *f3dHistGenVtxPtVsYCMSVsMultV0AOmegaMinus;      //! Generated Omega- Distrib, at triggered events level
  TH3F      *f3dHistGenVtxPtVsYCMSVsMultV0AOmegaPlus;       //! Generated Omega+ Distrib, at triggered events level

  TH3F      *f3dHistGenSelectedPtVsYVsMultV0APiMinus;      //! Generated Xi- Distrib, at event selection level
  TH3F      *f3dHistGenSelectedPtVsYVsMultV0APiPlus;       //! Generated Xi+ Distrib, at event selection level

  TH3F      *f3dHistGenSelectedPtVsYVsMultV0AXiMinus;      //! Generated Xi- Distrib, at event selection level
  TH3F      *f3dHistGenSelectedPtVsYVsMultV0AXiPlus;       //! Generated Xi+ Distrib, at event selection level
  TH3F      *f3dHistGenSelectedPtVsYVsMultV0AOmegaMinus;      //! Generated Omega- Distrib, at event selection level
  TH3F      *f3dHistGenSelectedPtVsYVsMultV0AOmegaPlus;       //! Generated Omega+ Distrib, at event selection level
  
  TH3F      *f3dHistGenPtVsYCMSVsMultV0AXiMinus;      //! Generated Xi- Distrib
  TH3F      *f3dHistGenPtVsYCMSVsMultV0AXiPlus;       //! Generated Xi+ Distrib
  TH3F      *f3dHistGenPtVsYCMSVsMultV0AOmegaMinus;      //! Generated Omega- Distrib
  TH3F      *f3dHistGenPtVsYCMSVsMultV0AOmegaPlus;       //! Generated Omega+ Distrib
  
  TH3F      *f3dHistGenSelectedPtVsYCMSVsMultV0AXiMinus;      //! Generated Xi- Distrib, at event selection level
  TH3F      *f3dHistGenSelectedPtVsYCMSVsMultV0AXiPlus;       //! Generated Xi+ Distrib, at event selection level
  TH3F      *f3dHistGenSelectedPtVsYCMSVsMultV0AOmegaMinus;      //! Generated Omega- Distrib, at event selection level
  TH3F      *f3dHistGenSelectedPtVsYCMSVsMultV0AOmegaPlus;       //! Generated Omega+ Distrib, at event selection level

    TH3F      *f3dHistGenSelectedPrimPtVsYCMSVsMultV0AXiMinus;      //! Generated Xi- Distrib, at event selection level
    TH3F      *f3dHistGenSelectedPrimPtVsYCMSVsMultV0AXiPlus;       //! Generated Xi+ Distrib, at event selection level
    TH3F      *f3dHistGenSelectedPrimPtVsYCMSVsMultV0AOmegaMinus;      //! Generated Omega- Distrib, at event selection level
    TH3F      *f3dHistGenSelectedPrimPtVsYCMSVsMultV0AOmegaPlus;       //! Generated Omega+ Distrib, at event selection level
    
//---> Multiplicity -> ZNA centrality
  
  TH3F      *f3dHistGenPtVsYVsMultZNAXiMinus;      //! Generated Xi- Distrib
  TH3F      *f3dHistGenPtVsYVsMultZNAXiPlus;       //! Generated Xi+ Distrib
  TH3F      *f3dHistGenPtVsYVsMultZNAOmegaMinus;      //! Generated Omega- Distrib
  TH3F      *f3dHistGenPtVsYVsMultZNAOmegaPlus;       //! Generated Omega+ Distrib
  
  TH3F      *f3dHistGenSelectedPtVsYVsMultZNAXiMinus;      //! Generated Xi- Distrib, at event selection level
  TH3F      *f3dHistGenSelectedPtVsYVsMultZNAXiPlus;       //! Generated Xi+ Distrib, at event selection level
  TH3F      *f3dHistGenSelectedPtVsYVsMultZNAOmegaMinus;      //! Generated Omega- Distrib, at event selection level
  TH3F      *f3dHistGenSelectedPtVsYVsMultZNAOmegaPlus;       //! Generated Omega+ Distrib, at event selection level
  
  TH3F      *f3dHistGenPtVsYCMSVsMultZNAXiMinus;      //! Generated Xi- Distrib
  TH3F      *f3dHistGenPtVsYCMSVsMultZNAXiPlus;       //! Generated Xi+ Distrib
  TH3F      *f3dHistGenPtVsYCMSVsMultZNAOmegaMinus;      //! Generated Omega- Distrib
  TH3F      *f3dHistGenPtVsYCMSVsMultZNAOmegaPlus;       //! Generated Omega+ Distrib
  
  TH3F      *f3dHistGenSelectedPtVsYCMSVsMultZNAXiMinus;      //! Generated Xi- Distrib, at event selection level
  TH3F      *f3dHistGenSelectedPtVsYCMSVsMultZNAXiPlus;       //! Generated Xi+ Distrib, at event selection level
  TH3F      *f3dHistGenSelectedPtVsYCMSVsMultZNAOmegaMinus;      //! Generated Omega- Distrib, at event selection level
  TH3F      *f3dHistGenSelectedPtVsYCMSVsMultZNAOmegaPlus;       //! Generated Omega+ Distrib, at event selection level
  
//---> Multiplicity -> TRK centrality
  
  TH3F      *f3dHistGenPtVsYVsMultTRKXiMinus;      //! Generated Xi- Distrib
  TH3F      *f3dHistGenPtVsYVsMultTRKXiPlus;       //! Generated Xi+ Distrib
  TH3F      *f3dHistGenPtVsYVsMultTRKOmegaMinus;      //! Generated Omega- Distrib
  TH3F      *f3dHistGenPtVsYVsMultTRKOmegaPlus;       //! Generated Omega+ Distrib
  
  TH3F      *f3dHistGenSelectedPtVsYVsMultTRKXiMinus;      //! Generated Xi- Distrib, at event selection level
  TH3F      *f3dHistGenSelectedPtVsYVsMultTRKXiPlus;       //! Generated Xi+ Distrib, at event selection level
  TH3F      *f3dHistGenSelectedPtVsYVsMultTRKOmegaMinus;      //! Generated Omega- Distrib, at event selection level
  TH3F      *f3dHistGenSelectedPtVsYVsMultTRKOmegaPlus;       //! Generated Omega+ Distrib, at event selection level
  
  TH3F      *f3dHistGenPtVsYCMSVsMultTRKXiMinus;      //! Generated Xi- Distrib
  TH3F      *f3dHistGenPtVsYCMSVsMultTRKXiPlus;       //! Generated Xi+ Distrib
  TH3F      *f3dHistGenPtVsYCMSVsMultTRKOmegaMinus;      //! Generated Omega- Distrib
  TH3F      *f3dHistGenPtVsYCMSVsMultTRKOmegaPlus;       //! Generated Omega+ Distrib
  
  TH3F      *f3dHistGenSelectedPtVsYCMSVsMultTRKXiMinus;      //! Generated Xi- Distrib, at event selection level
  TH3F      *f3dHistGenSelectedPtVsYCMSVsMultTRKXiPlus;       //! Generated Xi+ Distrib, at event selection level
  TH3F      *f3dHistGenSelectedPtVsYCMSVsMultTRKOmegaMinus;      //! Generated Omega- Distrib, at event selection level
  TH3F      *f3dHistGenSelectedPtVsYCMSVsMultTRKOmegaPlus;       //! Generated Omega+ Distrib, at event selection level
  
  //---> Multiplicity -> SPD centrality
  
  TH3F      *f3dHistGenPtVsYVsMultSPDXiMinus;      //! Generated Xi- Distrib
  TH3F      *f3dHistGenPtVsYVsMultSPDXiPlus;       //! Generated Xi+ Distrib
  TH3F      *f3dHistGenPtVsYVsMultSPDOmegaMinus;      //! Generated Omega- Distrib
  TH3F      *f3dHistGenPtVsYVsMultSPDOmegaPlus;       //! Generated Omega+ Distrib
  
  TH3F      *f3dHistGenSelectedPtVsYVsMultSPDXiMinus;      //! Generated Xi- Distrib, at event selection level
  TH3F      *f3dHistGenSelectedPtVsYVsMultSPDXiPlus;       //! Generated Xi+ Distrib, at event selection level
  TH3F      *f3dHistGenSelectedPtVsYVsMultSPDOmegaMinus;      //! Generated Omega- Distrib, at event selection level
  TH3F      *f3dHistGenSelectedPtVsYVsMultSPDOmegaPlus;       //! Generated Omega+ Distrib, at event selection level
  
  TH3F      *f3dHistGenPtVsYCMSVsMultSPDXiMinus;      //! Generated Xi- Distrib
  TH3F      *f3dHistGenPtVsYCMSVsMultSPDXiPlus;       //! Generated Xi+ Distrib
  TH3F      *f3dHistGenPtVsYCMSVsMultSPDOmegaMinus;      //! Generated Omega- Distrib
  TH3F      *f3dHistGenPtVsYCMSVsMultSPDOmegaPlus;       //! Generated Omega+ Distrib
  
  TH3F      *f3dHistGenSelectedPtVsYCMSVsMultSPDXiMinus;      //! Generated Xi- Distrib, at event selection level
  TH3F      *f3dHistGenSelectedPtVsYCMSVsMultSPDXiPlus;       //! Generated Xi+ Distrib, at event selection level
  TH3F      *f3dHistGenSelectedPtVsYCMSVsMultSPDOmegaMinus;      //! Generated Omega- Distrib, at event selection level
  TH3F      *f3dHistGenSelectedPtVsYCMSVsMultSPDOmegaPlus;       //! Generated Omega+ Distrib, at event selection level
  
  //---------------------
  
   TH1F      *fHistPVx;                      //! PVx distrib
   TH1F      *fHistPVy;                      //! PVy distrib
   TH1F      *fHistPVz;                      //! PVz distrib
   TH1F      *fHistPVxAnalysis;                      //! PVx distrib
   TH1F      *fHistPVyAnalysis;                      //! PVy distrib
   TH1F      *fHistPVzAnalysis;                      //! PVz distrib

   AliAnalysisTaskExtractPerformanceCascade(const AliAnalysisTaskExtractPerformanceCascade&);            // not implemented
   AliAnalysisTaskExtractPerformanceCascade& operator=(const AliAnalysisTaskExtractPerformanceCascade&); // not implemented
   
   ClassDef(AliAnalysisTaskExtractPerformanceCascade, 11);
};

#endif
