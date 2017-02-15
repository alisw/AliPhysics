/*************************************************************************
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
//
// Select events according to gap conditions, analyze two track events in pp
// collisions
//
// Author:
//  Paul Buehler <paul.buehler@oeaw.ac.at>

#ifndef AliAnalysisTaskCEP_H
#define AliAnalysisTaskCEP_H

#ifndef ALIANALYSISTASK_H
#include "AliAnalysisTaskSE.h"
#endif

class AliESDEvent;
class AliPIDResponse;
class AliPIDCombined;
class AliPhysicsSelection;

#include "AliEventCuts.h"
#include "AliCEPUtils.h"
#include "CEPEventBuffer.h"


class AliAnalysisTaskCEP : public AliAnalysisTaskSE
{
public:
  const static Int_t kTrackInfo = 6;

  // two  class  constructors
	AliAnalysisTaskCEP();
	AliAnalysisTaskCEP(const char* name, Long_t state = 0x0, Int_t numTracksMax = 7);
  // class  destructor
	virtual ~AliAnalysisTaskCEP();
  // called  once  at  beginning  or  runtime
	virtual void UserCreateOutputObjects();
  // called  for  each  event
	virtual void UserExec(Option_t *);
  // called  at  end  of  analysis
  // virtual void Terminate(Option_t* option);
  
private:

	AliAnalysisTaskCEP(const AliAnalysisTaskCEP  &p);
	AliAnalysisTaskCEP& operator=(const AliAnalysisTaskCEP  &p);

	// functions called by the UserExec(...), not to be called elsewhere!
	//-------------------------------------------------------------------
	Bool_t CheckInput();
	void PostOutputs();

  // events are save if the number of tracks is <= fnumTracksMax
  Int_t fnumTracksMax;

	// event information
	Int_t fRun;                     //  run number
  AliESDRun      *fESDRun;        //! esd run object
  AliESDEvent    *fESDEvent;      //! esd event object
  CEPEventBuffer *fCEPEvent;      //! event buffer
  TObjArray      *fTracks;        //! array of tracks
  TArrayI        *fTrackStatus;   //! array of TrackStatus
  TString         fLHCPeriod;     //! LHC period
    
  TVector3        fVtxPos;        // vertex position

	// analysis task status
	Long_t fAnalysisStatus; // stores the status bits used to specify
                          // processing steps to be performed
  // MC CEP system
  TLorentzVector fMCCEPSystem;
  
	// PID information
  AliPIDResponse *fPIDResponse;   //! esd pid object
  AliPIDCombined *fPIDCombined1;  //! PID Combined object with priors
  AliPIDCombined *fPIDCombined2;  //! PID Combined object without priors
  AliPIDCombined *fPIDCombined3;  //! PID Combined object without priors

	// some Ali tools
	AliTriggerAnalysis *fTrigger;           //! trigger object
	AliPhysicsSelection *fPhysicsSelection; //! physics selection object
  AliEventCuts *fEventCuts;               //! standard event cuts
  AliESDtrackCuts *fTrackCuts;            //! standard track cuts
  AliMultiplicitySelectionCP *fMartinSel; //! Martin's selection
  AliCEPUtils *fCEPUtil;                  //! object of type AliCEPUtils

  TH1F *fhStatsFlow;    //! histogram with event selection statistics
  
	// output objects
	TList *fHist;         //! output list (contains all histograms)
	TTree *fCEPtree;      //! Tree containing the detailed information

	ClassDef(AliAnalysisTaskCEP, 1);
  
};


#endif
