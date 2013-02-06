#ifndef ALIANALYSISTASKEXTRACTV0_H
#define ALIANALYSISTASKEXTRACTV0_H

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
//      AliAnalysisTaskExtractV0 class
//      ------------------------------
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
class AliPhysicsSelection;
class AliCFContainer;

//#include "TString.h"
//#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskExtractV0 : public AliAnalysisTaskSE {
 public:
	AliAnalysisTaskExtractV0();
	AliAnalysisTaskExtractV0(const char *name);
	virtual ~AliAnalysisTaskExtractV0();
	
	virtual void	 UserCreateOutputObjects();
	virtual void	 UserExec(Option_t *option);
	virtual void	 Terminate(Option_t *);
  void CheckChargeV0(AliESDv0 *thisv0);

  void SetIsNuclear           (Bool_t lIsNuclear   = kTRUE ) { fkIsNuclear   = lIsNuclear;   }
  void SetINT7Trigger         (Bool_t lSwitchINT7  = kTRUE ) { fkSwitchINT7  = lSwitchINT7; }
  void SetUseOnTheFly         (Bool_t lUseOnTheFly = kTRUE ) { fkUseOnTheFly = lUseOnTheFly; }
  void SetTakeAllTracks       (Bool_t lTakeAllTracks = kTRUE ) { fkTakeAllTracks = lTakeAllTracks; }
  void SetCentralityEstimator (TString lCentralityEstimator = "V0M" ) { fCentralityEstimator = lCentralityEstimator; }
  void SetLightWeightAnalysis (Bool_t lLightWeight = kTRUE) {fkLightWeight = lLightWeight;  }
  
 private:
				// Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
				// your data member object is created on the worker nodes and streaming is not needed.
				// http://root.cern.ch/download/doc/11InputOutput.pdf, page 14
	TList	*fListHistV0;	//! List of output objects
	TTree	*fTree;							//! Output Tree

  AliPIDResponse *fPIDResponse;     // PID response object
  AliESDtrackCuts *fESDtrackCuts;   // ESD track cuts used for primary track definition

  //Objects Controlling Task Behaviour 
  
  Bool_t fkIsNuclear;   // if true, replace multiplicity est. by centrality (default FALSE) 
  Bool_t fkSwitchINT7; // if true, skip FASTOnly (default FALSE)
  Bool_t fkUseOnTheFly; // if true, will use On-the-fly V0s instead of Offline V0s (default FALSE)
  Bool_t fkTakeAllTracks; // if true, no TPC crossed rows and ratio cut
  TString fCentralityEstimator; //Centrality Estimator String value (default V0M)
  
  Bool_t fkLightWeight; // if true, analysis output will exclude some non-fundamental
                        // debugging information. This creates smaller output. 

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
	Float_t fTreeVariablePtArmV0;//!
	Float_t fTreeVariableNegTotMomentum; //!	
	Float_t fTreeVariablePosTotMomentum; //!
	Float_t fTreeVariableNegdEdxSig; //!
	Float_t fTreeVariablePosdEdxSig; //!
	Float_t fTreeVariableNegEta; //!
	Float_t fTreeVariablePosEta; //!

	Float_t fTreeVariableNSigmasPosProton; //!
	Float_t fTreeVariableNSigmasPosPion; //! 
	Float_t fTreeVariableNSigmasNegProton; //!
	Float_t fTreeVariableNSigmasNegPion; //! 
	
	Float_t fTreeVariableDistOverTotMom;//!
	Int_t   fTreeVariableLeastNbrCrossedRows;//!
	Float_t fTreeVariableLeastRatioCrossedRowsOverFindable;//!
	Int_t   fTreeVariableMultiplicity ;//!

  Int_t   fTreeVariableRunNumber; //! 
  ULong64_t fTreeVariableEventNumber; //!

  Float_t fTreeVariableV0x; //!
  Float_t fTreeVariableV0y; //!
  Float_t fTreeVariableV0z; //!

  Float_t fTreeVariableV0Px; //!
  Float_t fTreeVariableV0Py; //!
  Float_t fTreeVariableV0Pz; //!

  Float_t fTreeVariablePVx; //!
  Float_t fTreeVariablePVy; //!
  Float_t fTreeVariablePVz; //!
  
  //Decay Length issue debugging: ULong_t with track status
  ULong64_t fTreeVariableNegTrackStatus;
  ULong64_t fTreeVariablePosTrackStatus;
  
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

  //Raw Data for Vertex Z position estimator change
	TH2F    *f2dHistMultiplicityVsVertexZBeforeTrigSel; 	        //! multiplicity distribution    
	TH2F    *f2dHistMultiplicityVsVertexZForTrigEvt;  		        //! multiplicity distribution
	TH2F    *f2dHistMultiplicityVsVertexZ;     					        //! multiplicity distribution
	TH2F    *f2dHistMultiplicityVsVertexZNoTPCOnly;			        //! multiplicity distribution
	TH2F    *f2dHistMultiplicityVsVertexZNoTPCOnlyNoPileup;			//! multiplicity distribution


	TH1F    *fHistPVx;     					        //! multiplicity distribution
	TH1F    *fHistPVy;     					        //! multiplicity distribution
	TH1F    *fHistPVz;     					        //! multiplicity distribution
	TH1F    *fHistPVxAnalysis;     					        //! multiplicity distribution
	TH1F    *fHistPVyAnalysis;     					        //! multiplicity distribution
	TH1F    *fHistPVzAnalysis;     					        //! multiplicity distribution
	TH1F    *fHistSwappedV0Counter;     					        //! Swapped V0 Counter

   AliAnalysisTaskExtractV0(const AliAnalysisTaskExtractV0&);            // not implemented
   AliAnalysisTaskExtractV0& operator=(const AliAnalysisTaskExtractV0&); // not implemented
  
   ClassDef(AliAnalysisTaskExtractV0, 11);
};

#endif
