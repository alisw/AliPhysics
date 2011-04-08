#ifndef ALIEPSELECTIONTASK_H
#define ALIEPSELECTIONTASK_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*****************************************************
//   Class AliEPSelectionTask
//   author: Alberica Toia, Johanna Gramling
//*****************************************************

#include "AliAnalysisTaskSE.h"

class TFile;
class TH1F;
class TH2F;
class TList;
class TString;
class TVector2;
class TObjArray;

class AliESDEvent;
class AliESDtrackCuts;
class AliESDtrack;
class AliEventplane;

class AliEPSelectionTask : public AliAnalysisTaskSE {

 public:

  AliEPSelectionTask();
  AliEPSelectionTask(const char *name);
  virtual ~AliEPSelectionTask();

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  
  TVector2 GetQ(AliEventplane* EP, TObjArray* event);
  void GetQsub(TVector2& Qsub1, TVector2& Qsub2, TObjArray* event);
  Double_t GetWeight(AliESDtrack* track);
  Double_t GetPhiWeight(AliESDtrack* track);

  virtual void  SetDebugLevel(Int_t level) {fDebug = level;}
  void SetInput(const char* input)         {fAnalysisInput = input;}
  void SetUseMCRP()			   {fUseMCRP = kTRUE;}
  void SetUsePhiWeight()		   {fUsePhiWeight = kTRUE;}
  void SetUsePtWeight()			   {fUsePtWeight = kTRUE;}
  void SetSaveTrackContribution()	   {fSaveTrackContribution = kTRUE;}
  void SetESDtrackCuts(TString status);
  void SetPhiDistribution(char* filename, char* listname);
  
 private:
   
  AliEPSelectionTask(const AliEPSelectionTask& ep);
  AliEPSelectionTask& operator= (const AliEPSelectionTask& ep); 

  Int_t    fDebug;	   		// Debug flag
  TString  fAnalysisInput; 		// "ESD", "AOD"
  TString  fStatus;			// "GLOBAL", "TPC"
  Bool_t   fUseMCRP;			// i.e. usable for Therminator, when MC RP is provided
  Bool_t   fUsePhiWeight;		// use of phi weights
  Bool_t   fUsePtWeight;		// use of pT weights
  Bool_t   fSaveTrackContribution;	// storage of contribution of each track to Q-Vector
  
  AliESDtrackCuts* fESDtrackCuts;
  
  TObjArray* ftracklist;		// list of accepted tracks for Q-Vector
  TH1F*	 fPhiDist;			// Phi distribution used to calculate phi weights

  TVector2* fQVector;			// Q-Vector of the event  
  Double_t* fQContributionX;		// array of the tracks' contributions to X component of Q-Vector - index = track ID
  Double_t* fQContributionY;		// array of the tracks' contributions to Y component of Q-Vector - index = track ID
  Double_t  fEventplaneQ; 		// Event plane angle from Q-Vector
  TVector2* fQsub1;			// Q-Vector of sub-event 1
  TVector2* fQsub2;			// Q-Vector of sub-event 2
  Double_t  fQsubRes;			// Difference of EP angles of subevents
  
  TList* fOutputList;  
  TH1F*  fHOutEventplaneQ;    		// control histogram: Event Plane angle
  TH1F*  fHOutPhi;			// control histogram: original phi distribution
  TH1F*	 fHOutPhiCorr;			// control histogram: corrected phi distribution
  TH2F*  fHOutsub1sub2;			// control histogram: correlation of EP from subevents
  TH2F*  fHOutNTEPRes;			// control histogram: Difference of EP angles of subevents vs Nch
  TH2F*  fHOutPTPsi;			// control histogram: Difference of EP angle and emission angle of track vs track pT
  TH2F*	 fHOutDiff;			// control histogram: Difference of MC RP and EP - only filled if fUseMCRP is true!
  TH2F*  fHOutleadPTPsi;		// control histogram: emission angle of leading pT track vs EP angle

  ClassDef(AliEPSelectionTask,2); 
};

#endif

