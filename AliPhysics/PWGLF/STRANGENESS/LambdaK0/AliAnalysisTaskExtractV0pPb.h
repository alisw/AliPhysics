#ifndef AliAnalysisTaskExtractV0pPb_H
#define AliAnalysisTaskExtractV0pPb_H

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
//      AliAnalysisTaskExtractV0pPb class
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
class AliAnalysisUtils;
class AliESDEvent;
class AliPhysicsSelection;
class AliCFContainer;

//#include "TString.h"
//#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskExtractV0pPb : public AliAnalysisTaskSE {
 public:
	AliAnalysisTaskExtractV0pPb();
	AliAnalysisTaskExtractV0pPb(const char *name);
	virtual ~AliAnalysisTaskExtractV0pPb();
	
	virtual void	 UserCreateOutputObjects();
	virtual void	 UserExec(Option_t *option);
	virtual void	 Terminate(Option_t *);

//---------------------------------------------------------------------------------------
//   Minimalistic Setters for configuration of task: pPb version, V0A
//---------------------------------------------------------------------------------------
  
  void SetTPCdEdxSelection ( Bool_t lTPCdEdxSelection = kTRUE ) { fTPCdEdxSelection = lTPCdEdxSelection; }
  
 private:
				// Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
				// your data member object is created on the worker nodes and streaming is not needed.
				// http://root.cern.ch/download/doc/11InputOutput.pdf, page 14
	TList	*fListHistV0;	//! List of output objects
	TTree	*fTree;				//! Output Tree filled with V0 candidates
	TTree	*fTreeEvents;	//! Output Tree filled with Events

  AliPIDResponse *fPIDResponse;     // PID response object
  
  //Objects Controlling Task Behaviour 
  
  //Extra selections
  Bool_t fTPCdEdxSelection; //Configuration to apply extra TPC dE/dx selection for better filling of tree

  //Variables for Tree
	Float_t fTreeVariableChi2V0;         //!
	Float_t fTreeVariableDcaV0Daughters; //!
	Float_t fTreeVariableDcaV0ToPrimVertex; //!
	Float_t fTreeVariableDcaPosToPrimVertex; //!
	Float_t fTreeVariableDcaNegToPrimVertex; //!
	Float_t fTreeVariableDCAV0ToPrimVertex; //!
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
  //This is only V0A centrality
	Float_t fTreeVariableCentrality ;//! 

  //Variables for Event Tree
  Float_t fTreeEventsCentrality; //!
 
//Note: TDistOverTotMom needs a mass hypothesis to be converted to proper decaylength.

//---------------------------------------------------------------------------------------  

  //V0A Centrality Distributions
	TH1F    *fHistCentralityProcessed;  //! All processed
	TH1F    *fHistCentralityTrigEvt;    //! Those selected by PS / Trigger
	TH1F    *fHistCentralityHasVtx;     //! Those that have a well-established vertex
	TH1F    *fHistCentralityVtxZ;			  //! Those whose vertex falls within |z|<10cm 


   AliAnalysisTaskExtractV0pPb(const AliAnalysisTaskExtractV0pPb&);            // not implemented
   AliAnalysisTaskExtractV0pPb& operator=(const AliAnalysisTaskExtractV0pPb&); // not implemented
  
   ClassDef(AliAnalysisTaskExtractV0pPb, 11);
};

#endif

