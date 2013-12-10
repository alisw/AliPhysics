#ifndef ALIANALYSISTASKQAV0AOD_H
#define ALIANALYSISTASKQAV0AOD_H

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
//      AliAnalysisTaskQAV0AOD class
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

class AliAnalysisTaskQAV0AOD : public AliAnalysisTaskSE {
 public:
	AliAnalysisTaskQAV0AOD();
	AliAnalysisTaskQAV0AOD(const char *name);
	virtual ~AliAnalysisTaskQAV0AOD();
	
	virtual void	 UserCreateOutputObjects();
	virtual void	 UserExec(Option_t *option);
	virtual void	 Terminate(Option_t *);

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


  //Objects Controlling Task Behaviour 
  AliPIDResponse *fPIDResponse;     // PID response object

  //Objects Controlling Task Behaviour: has to be streamed! 
  Double_t  fV0Sels[7];           // Array to store the 7 values for the different selections V0 related

   AliAnalysisTaskQAV0AOD(const AliAnalysisTaskQAV0AOD&);            // not implemented
   AliAnalysisTaskQAV0AOD& operator=(const AliAnalysisTaskQAV0AOD&); // not implemented
  
   ClassDef(AliAnalysisTaskQAV0AOD, 11);
};

#endif
