#ifndef AliAnalysisTaskExtractV0AODRun2_H
#define AliAnalysisTaskExtractV0AODRun2_H

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

//-----------------------------------------------------------------
//      AliAnalysisTaskExtractV0AODRun2 class
//      ---------------------------------
//
//    Please see cxx file for more details.   
//             
//-----------------------------------------------------------------

// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// --- This version: 23rd March 2012 
//
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class TList;
class TH1F;
class TH2F;
class TH3F;
class TVector3;
class THnSparse;

class AliESDpid;
class AliESDtrackCuts;
class AliESDEvent;
class AliAODEvent;
class AliPhysicsSelection;
class AliCFContainer;

//#include "TString.h"
//#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskExtractV0AODRun2 : public AliAnalysisTaskSE {
 public:
	AliAnalysisTaskExtractV0AODRun2();
	AliAnalysisTaskExtractV0AODRun2(const char *name);
	virtual ~AliAnalysisTaskExtractV0AODRun2();
	
	virtual void	 UserCreateOutputObjects();
	virtual void	 UserExec(Option_t *option);
	virtual void	 Terminate(Option_t *);
  void CheckChargeV0(AliESDv0 *thisv0);

  void SetIsNuclear           (Bool_t lIsNuclear   = kTRUE  ) { fkIsNuclear   = lIsNuclear;   }
  void SetIsLowEnergyPP       (Bool_t lLowEnergyPP = kTRUE  ) { fkLowEnergyPP = lLowEnergyPP; }
  void SetUseOnTheFly         (Bool_t lUseOnTheFly = kTRUE  ) { fkUseOnTheFly = lUseOnTheFly; }
    void SetTriggerMask         (TString lTriggerMask = "kMB" ) { fTriggerMask  = lTriggerMask; }
  void SetLowE         (Bool_t lLowE = kTRUE ) { fkLowE  = lLowE; }
    void SetSaveAllInvMasses  (Bool_t lSaveAllInvMasses = kTRUE ) { fkSaveAllInvMasses  = lSaveAllInvMasses; }
    void SetPreSelect  (Bool_t lPreSelect = kTRUE ) { fkPreSelect  = lPreSelect; }

 private:
				// Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
				// your data member object is created on the worker nodes and streaming is not needed.
				// http://root.cern.ch/download/doc/11InputOutput.pdf, page 14
	TList	*fListHistV0;	//! List of output objects
	TTree	*fTree;							//! Output Tree

	AliPIDResponse *fPIDResponse;     // PID response object

  //Objects Controlling Task Behaviour 
  
  Bool_t fkIsNuclear;   // if true, replace multiplicity est. by centrality (default FALSE)
  Bool_t fkLowEnergyPP; // if true, skip FASTOnly (default FALSE)
  Bool_t fkUseOnTheFly; // if true, will use On-the-fly V0s instead of Offline V0s (default FALSE)

  TString fTriggerMask; // Selected trigger mask (kMB, kINT7, kINT8, kAnyINT)
    
    Bool_t fkLowE;   // if true, use old centrality framework (default FALSE)
    Bool_t fkSaveAllInvMasses; //Do not cut on invariant masses for the tree (default FALSE)
    Bool_t fkPreSelect; //Pre-select the V0s to reduce the size of the tree (default FALSE)
    
  //Variables for Tree
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
    Float_t fTreeVariableAlphaV0BoostAsK0; //!
	Float_t fTreeVariablePtArmV0;//!
    
	Float_t fTreeVariableNegEta; //!
	Float_t fTreeVariablePosEta; //!

	Float_t fTreeVariableNSigmasPosProton; //!
	Float_t fTreeVariableNSigmasPosPion; //! 
	Float_t fTreeVariableNSigmasNegProton; //!
	Float_t fTreeVariableNSigmasNegPion; //!
    
    Float_t fTreeVariableNegTPCSignal; //!
    Float_t fTreeVariablePosTPCSignal; //!
    Float_t fTreeVariableNegInnerP; //!
    Float_t fTreeVariablePosInnerP; //!
	
	Float_t fTreeVariableDistOverTotMom;//!
	Int_t   fTreeVariableLeastNbrCrossedRows;//!
	Float_t fTreeVariableLeastRatioCrossedRowsOverFindable;//!
	Int_t   fTreeVariableCentrality ;//!

  Int_t   fTreeVariableRunNumber; //! 
  Int_t fTreeVariablePeriodNumber; //!
    Int_t fTreeVariableOrbitNumber; //!
    Int_t fTreeVariableBunchCrossNumber; //!

    Float_t fTreeVariableV0Px; //!
    Float_t fTreeVariableV0Py; //!
    Float_t fTreeVariableV0Pz; //!
    
    Float_t fTreeVariableDecayMomLambda; //!
    Float_t fTreeVariableDecayMomAntiLambda; //!
    Float_t fTreeVariableDecayMomK0Short; //!

    Float_t fTreeVariablePrimVX; //!
    Float_t fTreeVariablePrimVY; //!
    Float_t fTreeVariablePrimVZ; //!
    
    Bool_t fTreeVariableCowboy; //!
    

//Note: TDistOverTotMom needs a mass hypothesis to be converted to proper decaylength.

	TH1F    *fHistV0MultiplicityBeforeTrigSel;             //! V0 multiplicity distribution
	TH1F    *fHistV0MultiplicityForTrigEvt;                //! V0 multiplicity distribution
	TH1F    *fHistV0MultiplicityForSelEvt;                 //! V0 multiplicity distribution
	TH1F    *fHistV0MultiplicityForSelEvtNoTPCOnly;        //! V0 multiplicity distribution
	TH1F    *fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup;//! V0 multiplicity distribution

	TH1F    *fHistMultiplicityBeforeTrigSel; 	        //! multiplicity distribution    
	TH1F    *fHistMultiplicityForTrigEvt;  		        //! multiplicity distribution
	TH1F    *fHistMultiplicity;     					        //! multiplicity distribution
	TH1F    *fHistMultiplicityNoTPCOnly;			        //! multiplicity distribution
	TH1F    *fHistMultiplicityNoTPCOnlyNoPileup;			//! multiplicity distribution

    TH1F    *fHistCentralityBeforeTrigSel; 	        //! Centrality distribution
    TH1F    *fHistCentralityForTrigEvt;  		        //! Centrality distribution
    TH1F    *fHistCentrality;     					        //! Centrality distribution
    TH1F    *fHistCentralityNoTPCOnly;			        //! Centrality distribution
    TH1F    *fHistCentralityNoTPCOnlyNoPileup;			//! Centrality distribution
    
	TH1F    *fHistPVx;     					        //! Primary Vertex X
	TH1F    *fHistPVy;     					        //! Primary Vertex Y
	TH1F    *fHistPVz;     					        //! Primary Vertex Z
	TH1F    *fHistPVxAnalysis;     					        //! Primary Vertex X
	TH1F    *fHistPVyAnalysis;     					        //! Primary Vertex Y
	TH1F    *fHistPVzAnalysis;     					        //! Primary Vertex Z
	TH1F    *fHistSwappedV0Counter;     					        //! Swapped V0 Counter

    
    TH1D    *fHistRunNumber;                        //! run number histogram
    TH1D    *fHistPeriodNumber;                     //! period number histogram
    TH1D    *fHistOrbitNumber;                      //! orbit number histogram
    TH1D    *fHistBunchCrossNumber;                 //! bunch cross number histogram
    
    // histograms for TPC signal
    TH2D			   *f2dHistTPCSignal;   //! TPC Signal vs Inner Wall Momentum - tracks
    TH2D			   *f2dHistTPCSignalPionsFromK0;   //! TPC Signal vs Inner Wall Momentum - Pions from K0
    TH2D			   *f2dHistTPCSignalPionsFromLambda;   //! TPC Signal vs Inner Wall Momentum - Pions from Lambda and AntiLambda
    TH2D			   *f2dHistTPCSignalProtonsFromLambda;   //! TPC Signal vs Inner Wall Momentum - Protons from Lambda and AntiLambda
    
    //3D histos
    TH3F* f3dHistInvMassVsPtVsCentLambda;      //! Lambda
    TH3F* f3dHistInvMassVsPtVsCentAntiLambda;  //! AntiLambda
    TH3F* f3dHistInvMassVsPtVsCentK0Short;     //! K0Short
    
    TH3F* f3dHistInvMassVsPtVsCentLambdaCowboy;      //! Lambda Cowboy
    TH3F* f3dHistInvMassVsPtVsCentAntiLambdaCowboy;  //! AntiLambda Cowboy
    TH3F* f3dHistInvMassVsPtVsCentK0ShortCowboy;     //! K0Short Cowboy

    TH3F* f3dHistInvMassVsPtVsCentLambdaSailor;      //! Lambda Sailor
    TH3F* f3dHistInvMassVsPtVsCentAntiLambdaSailor;  //! AntiLambda Sailor
    TH3F* f3dHistInvMassVsPtVsCentK0ShortSailor;     //! K0Short Sailor

    
   AliAnalysisTaskExtractV0AODRun2(const AliAnalysisTaskExtractV0AODRun2&);            // not implemented
   AliAnalysisTaskExtractV0AODRun2& operator=(const AliAnalysisTaskExtractV0AODRun2&); // not implemented
  
   ClassDef(AliAnalysisTaskExtractV0AODRun2, 11);
};

#endif
