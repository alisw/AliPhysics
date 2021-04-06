#ifndef AliAnalysisTaskPbPbTree_SingleMuons_H
#define AliAnalysisTaskPbPbTree_SingleMuons_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"
#include "TTreeStream.h"

class TObjArray;
class AliVParticle;
class AliAODEvent;
class TLorentzVector;
class AliMuonTrackCuts;
class AliTriggerAnalysis;
class AliVMultiplicity;
class AliAODTracklets;

//class TList;

class AliAnalysisTaskPbPbTree_SingleMuons: public AliAnalysisTaskSE {
  public:

  AliAnalysisTaskPbPbTree_SingleMuons();
  AliAnalysisTaskPbPbTree_SingleMuons(const char *name);
  virtual ~AliAnalysisTaskPbPbTree_SingleMuons();

  void UserCreateOutputObjects();
  void UserExec(Option_t *option);
  void Terminate(Option_t *);
  virtual void NotifyRun();

  void SetBeamEnergy(Double_t en) {fBeamEnergy=en;}
  void SetAnalysisType(const char* type) {fkAnalysisType=type;}
  void SetPeriod(TString period) {fPeriod=period;}

 private:
  AliAnalysisTaskPbPbTree_SingleMuons(const AliAnalysisTaskPbPbTree_SingleMuons&);
  AliAnalysisTaskPbPbTree_SingleMuons& operator=(const AliAnalysisTaskPbPbTree_SingleMuons&);

 //protected:

  TTree     *fOutputTree;      //! tree output
  Double_t   fNevt;             // event counter
  TH1D      *fhNEv;            //! histo

  Double_t fBeamEnergy;            // Energy of the beam (required for the CS angle)

  const char* fkAnalysisType;      //ESD or AOD based analysis
  TString fPeriod;                 //period
  Int_t fRun;                      // run number - for calibration

  Int_t fCountTotEv;               // counter
  Int_t fCountTrigger;             // counter
  Int_t fCountCINT7;               // counter
  Int_t fCountCMUL7;               //counter
  Int_t fCountCMLL7;               //counter
  Int_t fCountCMSL7;               //counter
  Int_t fCountCMSH7;               //counter

  Int_t		fNMuons;		              // muon tracks in the event
  Int_t		fNTracks;		              // tracks in the event
  Int_t		fNTracklets;		          // spd tracklets
  Int_t		fNContributors;	        	// n contributors
  Bool_t  fIsPhysSelected;    //check on phys selection
  AliAODEvent*  fAODEvent;          //! AOD event
  char   	fTrigClass[800];	        // fired trigger classes

  AliMuonTrackCuts* fMuonTrackCuts;
  Double_t	fVertex[3];		          // x,y,z vertex
  Float_t   fPercentV0M;             //! percentile V0
  Double_t	fPt[1500];		          // single mu pT
  Double_t	fE[1500];		 	          // single mu E
  Double_t	fPx[1500];		          // single mu px
  Double_t	fPy[1500];		          // single mu py
  Double_t	fPz[1500];		          // single mu pz
  Double_t	fY[1500];		            // single mu y
  Double_t	fEta[1500];		          // single mu eta
  Double_t	fPhi[1500];		          // single mu phi
  Int_t		fMatchTrig[1500];		      // single mu match trigger
  Double_t	fTrackChi2[1500];		    // single mu chi2 track
  Double_t	fMatchTrigChi2[1500];	  // single mu chi2 of match trigger
  Double_t	fDCA[1500];		          // single mu DCA
  Int_t	        fCharge[1500];	  	// single mu charge
  Double_t	fRAtAbsEnd[1500];		    // single mu distance from beam center at end abs
  Int_t         fpDCA[1500];               //pDCA

 ClassDef(AliAnalysisTaskPbPbTree_SingleMuons,1);
};

#endif
