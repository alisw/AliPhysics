#ifndef AliAnalysisTaskQAV0_H
#define AliAnalysisTaskQAV0_H

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
//      AliAnalysisTaskQAV0 class
//      -------------------------
//
//    Please see cxx file for more details.   
//             
//-----------------------------------------------------------------

// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// --- This version: 23 Oct 2013
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
#include "AliPIDResponse.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskQAV0 : public AliAnalysisTaskSE {
 public:
	AliAnalysisTaskQAV0();
	AliAnalysisTaskQAV0(const char *name);
	virtual ~AliAnalysisTaskQAV0();
	
	virtual void	 UserCreateOutputObjects();
	virtual void	 UserExec(Option_t *option);
	virtual void	 Terminate(Option_t *);

//---------------------------------------------------------------------------------------
  //Task Configuration: Meant to enable quick re-execution of vertexer if needed
  void SetRunV0Vertexer ( Bool_t lRunV0Vertexer = kTRUE) { fkRunV0Vertexer = lRunV0Vertexer; }
  void SetUseLightV0Vertexer ( Bool_t lGoLightV0Vertexer = kTRUE) { fkUseLightV0Vertexer = lGoLightV0Vertexer; }
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
//Setters for the V0 Extraction
  void SetV0SelectionMaxChisquare   ( Double_t lParameter ){ fV0Sels[0] = lParameter; }
  void SetV0SelectionDCAFirstToPV   ( Double_t lParameter ){ fV0Sels[1] = lParameter; }
  void SetV0SelectionDCASecondtoPV  ( Double_t lParameter ){ fV0Sels[2] = lParameter; }
  void SetV0SelectionDCAV0Daughters ( Double_t lParameter ){ fV0Sels[3] = lParameter; }
  void SetV0SelectionCosinePA       ( Double_t lParameter ){ fV0Sels[4] = lParameter; }
  void SetV0SelectionMinRadius      ( Double_t lParameter ){ fV0Sels[5] = lParameter; }
  void SetV0SelectionMaxRadius      ( Double_t lParameter ){ fV0Sels[6] = lParameter; }
//---------------------------------------------------------------------------------------
//Setters for dE/dx selection
  void SetTPCdEdxSelection ( Double_t lParameter ) { fdEdxCut = lParameter; }
//---------------------------------------------------------------------------------------
  void SetSelectedTriggerClass(AliVEvent::EOfflineTriggerTypes trigType) { fTrigType = trigType;}
//---------------------------------------------------------------------------------------
    
 private:
				// Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
				// your data member object is created on the worker nodes and streaming is not needed.
				// http://root.cern.ch/download/doc/11InputOutput.pdf, page 14

  //Top Structure
	TList	*fOutput;	              //! List of output objects

  //Per-Event Histograms
  TH1D* fHistEvent;             //! Event Selection Histogram (no further info)

  //Per-Candidate Histograms
  
  //Base Topological Variables for All Reconstructed V0 Vertices
  TH1D *fHistTopDCANegToPV;       //!
  TH1D *fHistTopDCAPosToPV;       //!
  TH1D *fHistTopDCAV0Daughters;   //!
  TH1D *fHistTopCosinePA;         //!
  TH1D *fHistTopV0Radius;         //!

  //Zoomed into Selection criteria (fV0Sels)
  TH1D *fHistSelectedTopDCANegToPV;       //!
  TH1D *fHistSelectedTopDCAPosToPV;       //!
  TH1D *fHistSelectedTopDCAV0Daughters;   //!
  TH1D *fHistSelectedTopCosinePA;         //!
  TH1D *fHistSelectedTopV0Radius;         //!

  //Histograms for Storing Invariant Mass: 2D 
  //Stored only if the 5 topological selections in fV0Sels are satisfied 
  TH2D *f2dHistInvMassK0Short;      //!
  TH2D *f2dHistInvMassLambda;       //!
  TH2D *f2dHistInvMassAntiLambda;   //!

  //With dE/dx Selection (extra)
  TH2D *f2dHistInvMassWithdEdxK0Short;      //!
  TH2D *f2dHistInvMassWithdEdxLambda;       //!
  TH2D *f2dHistInvMassWithdEdxAntiLambda;   //!

    //With dE/dx Selection + are cowboys
    TH2D *f2dHistInvMassWithdEdxCowboysK0Short;      //!
    TH2D *f2dHistInvMassWithdEdxCowboysLambda;       //!
    TH2D *f2dHistInvMassWithdEdxCowboysAntiLambda;   //!
    
    //With dE/dx Selection + are sailors
    TH2D *f2dHistInvMassWithdEdxSailorsK0Short;      //!
    TH2D *f2dHistInvMassWithdEdxSailorsLambda;       //!
    TH2D *f2dHistInvMassWithdEdxSailorsAntiLambda;   //!
    
  //dEdx QA Histograms (extra) 
  //PIDFrameWork
  TH2D *f2dHistResponseNegativeAsPion;   //! 
  TH2D *f2dHistResponseNegativeAsProton; //! 
  TH2D *f2dHistResponsePositiveAsPion;   //! 
  TH2D *f2dHistResponsePositiveAsProton; //! 

  //Raw Stuff, Clean Version: Proton and Pion From Lambdas
  //Potentially useful for checking calibration of dE/dx 
  TH2D *f2dHistdEdxSignalPionFromLambda;   //! 
  TH2D *f2dHistdEdxSignalProtonFromLambda; //! 
  TH2D *f2dHistResponsePionFromLambda;     //! 
  TH2D *f2dHistResponseProtonFromLambda;   //!
    
    AliVEvent::EOfflineTriggerTypes fTrigType; // trigger type

  AliPIDResponse *fPIDResponse;     // PID response object

  //Objects Controlling Task Behaviour: has to be streamed! 
  Bool_t    fkRunV0Vertexer;      //if true, re-run vertexer with loose cuts. CARE MUST BE TAKEN in PbPb!
  Bool_t    fkUseLightV0Vertexer;  //if true, use light vertexers
  Double_t  fV0VertexerSels[7];     // Array to store the 7 values for the different selections V0 related
  Double_t  fV0Sels[7];           // Array to store the 7 values for the different selections V0 related

  //Variables controlling task behaviour (don't "//!" them!)
  //For setting a dEdx cut: can be strict... (will anyhow depend on analysis!) 
  Double_t fdEdxCut; 

   AliAnalysisTaskQAV0(const AliAnalysisTaskQAV0&);            // not implemented
   AliAnalysisTaskQAV0& operator=(const AliAnalysisTaskQAV0&); // not implemented
  
   ClassDef(AliAnalysisTaskQAV0, 11);
};

#endif
