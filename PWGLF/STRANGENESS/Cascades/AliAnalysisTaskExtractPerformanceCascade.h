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
  
 private:
        // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
        // your data member object is created on the worker nodes and streaming is not needed.
        // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14
  TList  *fListHist;  //! List of Cascade histograms
  TTree  *fTreeCascade;              //! Output Tree, Cascades

  //Objects that have to be streamed:
  AliPIDResponse *fPIDResponse;     // PID response object
  AliESDtrackCuts *fESDtrackCuts;   // ESD track cuts used for primary track definition

  //Objects Controlling Task Behaviour 
  // (have to be streamed too or configuration is lost)
  
  Bool_t fkIsNuclear;   //if true, replace multiplicity est. by centrality (default FALSE) 
  Bool_t fkSwitchINT7 ; //if true, skip FASTOnly (default FALSE)
  Double_t fpArapidityShift; //pA rapidity shift (should be 0.465, usually)
  TString fCentralityEstimator; //Centrality Estimator String value (default V0M)

	Double_t        fV0Sels[7];                     // Array to store the 7 values for the different selections V0 related
	Double_t        fCascSels[8];                   // Array to store the 8 values for the different selections Casc. related


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
