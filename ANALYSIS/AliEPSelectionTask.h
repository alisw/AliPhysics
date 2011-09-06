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

class AliESDEvent;
class AliESDtrackCuts;
class AliESDtrack;
class AliEventplane;
class AliOADBContainer;

class AliEPSelectionTask : public AliAnalysisTaskSE {

 public:
  
  enum ResoMethod{kRandom,kEta};

  AliEPSelectionTask();
  AliEPSelectionTask(const char *name);
  virtual ~AliEPSelectionTask();

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  
  TVector2 GetQ(AliEventplane* EP, TObjArray* event);
  void GetQsub(TVector2& Qsub1, TVector2& Qsub2, TObjArray* event,AliEventplane* EP);
  Double_t GetWeight(TObject* track1);
  Double_t GetPhiWeight(TObject* track1);

  virtual void  SetDebugLevel(Int_t level)   {fDebug = level;}
  void SetInput(const char* input)           {fAnalysisInput = input;}
  void SetUseMCRP()			     {fUseMCRP = kTRUE;}
  void SetUsePhiWeight()		     {fUsePhiWeight = kTRUE;}
  void SetUsePtWeight()			     {fUsePtWeight = kTRUE;}
  void SetSaveTrackContribution()	     {fSaveTrackContribution = kTRUE;}
  void SetTrackType(TString tracktype);
  void SetPhiDist();
  void SetPersonalESDtrackCuts(AliESDtrackCuts* trackcuts);
  void SetPersonalAODtrackCuts(UInt_t filterbit = 1, Float_t etalow = -0.8, Float_t etaup = 0.8, Float_t ptlow = 0.5, Float_t ptup = 20.);
  void SetPersonalPhiDistribution(const char* filename, char* listname);
  void SetEtaGap(Float_t etagap)             {fEtaGap = etagap;}
  void SetSubeventsSplitMethod(Int_t method) {fSplitMethod = method;}
  
 private:
   
  AliEPSelectionTask(const AliEPSelectionTask& ep);
  AliEPSelectionTask& operator= (const AliEPSelectionTask& ep); 

  TObjArray* GetAODTracksAndMaxID(AliAODEvent* aod, Int_t& maxid);

  TString  fAnalysisInput; 		// "ESD", "AOD"
  TString  fTrackType;			// "GLOBAL", "TPC"
  Bool_t   fUseMCRP;			// i.e. usable for Therminator, when MC RP is provided
  Bool_t   fUsePhiWeight;		// use of phi weights
  Bool_t   fUsePtWeight;		// use of pT weights
  Bool_t   fSaveTrackContribution;	// storage of contribution of each track to Q-Vector
  Bool_t   fUserphidist;		// bool, if personal phi distribution should be used
  Bool_t   fUsercuts;			// bool, if personal cuts should be used
  Int_t    fRunNumber;			// runnumber
  UInt_t   fAODfilterbit;               // AOD filter bit for AOD track selection  
  Float_t  fEtaGap;                     // Eta Gap between Subevent A and B
  Int_t    fSplitMethod;                // Splitting Method for subevents
  

  AliESDtrackCuts* fESDtrackCuts;       // track cuts
  
  AliOADBContainer* fEPContainer;	//! OADB Container
  TH1F*	 fPhiDist;			// Phi distribution used to calculate phi weights

  TVector2* fQVector;			//! Q-Vector of the event  
  Double_t* fQContributionX;		//! array of the tracks' contributions to X component of Q-Vector - index = track ID
  Double_t* fQContributionY;		//! array of the tracks' contributions to Y component of Q-Vector - index = track ID
  Double_t  fEventplaneQ; 		//! Event plane angle from Q-Vector
  TVector2* fQsub1;			//! Q-Vector of sub-event 1
  TVector2* fQsub2;			//! Q-Vector of sub-event 2
  Double_t  fQsubRes;			//! Difference of EP angles of subevents
  
  TList* fOutputList;                   // Output histograms
  TH1F*  fHOutEventplaneQ;    		//! control histogram: Event Plane angle
  TH1F*  fHOutPhi;			//! control histogram: original phi distribution
  TH1F*	 fHOutPhiCorr;			//! control histogram: corrected phi distribution
  TH2F*  fHOutsub1sub2;			//! control histogram: correlation of EP from subevents
  TH2F*  fHOutNTEPRes;			//! control histogram: Difference of EP angles of subevents vs Nch
  TH2F*  fHOutPTPsi;			//! control histogram: Difference of EP angle and emission angle of track vs track pT
  TH2F*	 fHOutDiff;			//! control histogram: Difference of MC RP and EP - only filled if fUseMCRP is true!
  TH2F*  fHOutleadPTPsi;		//! control histogram: emission angle of leading pT track vs EP angle

  ClassDef(AliEPSelectionTask,4); 
};

#endif

