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
	AliAnalysisTaskCEP(const char* name,
    Long_t state=AliCEPBase::kBitConfigurationSet,
    Int_t rnummin=100000, Int_t rnummax=300000,
    Int_t numTracksMax=6,
    Double_t fracDG=1.0, Double_t fracNDG=0.0,
    UInt_t ETmaskDG=AliCEPBase::kETBaseLine, UInt_t ETpatternDG=AliCEPBase::kETBaseLine,
    UInt_t ETmaskNDG=AliCEPBase::kETBaseLine, UInt_t ETpatternNDG=AliCEPBase::kETBaseLine,
    UInt_t TTmask=AliCEPBase::kTTBaseLine,  UInt_t TTpattern=AliCEPBase::kTTBaseLine);
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
  Bool_t IsSTGFired(TBits* fFOmap,Int_t dphiMin=0,Int_t dphiMax=10);

  // events are saved if (ET=Event test, TT=Track test)
  // . ET conditions are met (conditions for DG and NDG)
  // . 0 < nTrack='number of tracks meeting the TT conditions' <= fnumTracksMax
  Int_t frnummin, frnummax;
  Int_t fnumTracksMax;
  Double_t ffracDG, ffracNDG;
  UInt_t fETmaskDG,  fETpatternDG;
  UInt_t fETmaskNDG, fETpatternNDG;
  UInt_t fTTmask, fTTpattern;

	// event information
	Int_t fRun;                     //  run number
  AliESDRun      *fESDRun;        //! esd run object
  Bool_t          fisESD;         //!
  Bool_t          fisAOD;         //!
  AliVEvent      *fEvent;         //! event object
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

  TList *flQArnum;      //! list of QA histograms for QA vs rnum study
  TList *flBBFlag;      //! list of QA histograms for BBFlag study
  TList *flSPDpileup;   //! list of QA histograms for SPD pile-up study
  TList *flnClunTra;    //! list of QA histograms for nClunTra BG rejection
  TList *flVtx     ;    //! list of QA histograms for vertex selection
  TList *flV0      ;    //! list of QA histograms for V0 study
  TList *flFMD     ;    //! list of QA histograms for FMD study
  TH1F *fhStatsFlow;    //! histogram with event selection statistics
  
	// output objects
	TList *fHist;         //! output list (contains all histograms)
	TTree *fCEPtree;      //! Tree containing the detailed information

	ClassDef(AliAnalysisTaskCEP, 1);
  
};


#endif
