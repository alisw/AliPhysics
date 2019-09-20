#ifndef AliAnalysisTaskTree_MCut_H
#define AliAnalysisTaskTree_MCut_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"
#include "TTreeStream.h"//why?

class TObjArray;
class AliVParticle;
class AliAODEvent;
class TLorentzVector;
class AliMuonTrackCuts;

//class TList;

class AliAnalysisTaskTree_MCut : public AliAnalysisTaskSE {
  public:

  AliAnalysisTaskTree_MCut();
  AliAnalysisTaskTree_MCut(const char *name);
  virtual ~AliAnalysisTaskTree_MCut();

  void UserCreateOutputObjects();
  void UserExec(Option_t *option);
  void Terminate(Option_t *);
  virtual void NotifyRun(); 
  
  void SetBeamEnergy(Double_t en) {fBeamEnergy=en;}
  void SetAnalysisType(const char* type) {fkAnalysisType=type;}
  void SetPeriod(TString period) {fPeriod=period;}
    
 private:
  AliAnalysisTaskTree_MCut(const AliAnalysisTaskTree_MCut&);
  AliAnalysisTaskTree_MCut& operator=(const AliAnalysisTaskTree_MCut&);
     
 //protected:
     
  TTree         *fOutputTree	;      //! tree output
  Double_t      fNevt		;      // event counter

  Double_t fBeamEnergy;   // Energy of the beam (required for the CS angle)    
  const char* fkAnalysisType; //ESD or AOD based analysis
  TString fPeriod; //period
  Int_t fCountTotEv;   // counter
  Int_t fCountTrigger;   // counter
  Int_t fCountCINT7;   // counter
  Int_t fCountCMUL7; //counter
  Int_t fCountCMLL7; //counter
  Int_t fCountCMSL7; //counter
  Int_t fCountCMSH7; //counter
  
  Int_t		fNMuons;		// muon tracks in the event
  Int_t		fNTracklets;		// spd tracklets
  Int_t		fNContributors;		// n contributors
  Int_t		fNDimu;			// dimuons in the event
  Bool_t        fIsPhysSelected;   //check on phys selection 
  AliAODEvent*  fAODEvent;      //! AOD event 
  char   	fTrigClass[800];	// fired trigger classes
  UInt_t        finpmask;               // trigger input mask
  
  AliMuonTrackCuts* fMuonTrackCuts;
  Double_t	fVertex[3];		// x,y,z vertex
  Double_t	fPt[200];		 // single mu pT
  Double_t	fE[200];			// single mu E
  Double_t	fPx[200];		// single mu px
  Double_t	fPy[200];		// single mu py
  Double_t	fPz[200];		// single mu pz
  Double_t	fY[200];		// single mu y
  Double_t	fEta[200];		// single mu eta
  Int_t		fMatchTrig[200];		// single mu match trigger
  Double_t	fTrackChi2[200];		// single mu chi2 track
  Double_t	fMatchTrigChi2[200];	// single mu chi2 of match trigger
  Double_t	fDCA[200];		// single mu DCA
  Int_t	        fCharge[200];		// single mu charge
  Double_t	fRAtAbsEnd[200];		// single mu distance from beam center at end abs
  Int_t	        fpDCA[200];             //pDCA
  
  Double_t	fDimuPt[100];			    // dimuon pT
  Double_t	fDimuPx[100]; 		    // dimuon px
  Double_t	fDimuPy[100]; 		    // dimuon py
  Double_t	fDimuPz[100]; 		    // dimuon pz
  Double_t	fDimuY[100];			// dimuon y
  Double_t	fDimuMass[100];		// dimuon invariant mass
  Int_t	        fDimuCharge[100];		// dimuon charge
  Int_t	        fDimuMatch[100];		// dimuon match

  Double_t      fDimuCostHE[100]; // cost Helicty frame
  Double_t      fDimuPhiHE[100]; // phi Helicty fame
  Double_t      fDimuCostCS[100]; // cost Collins-Soper
  Double_t      fDimuPhiCS[100]; // phi Collins-Soper
  Int_t		fDimuMu[100][2];	// reference to single mus
  
  
//  TList *fOutput;  //!< List of histograms for data
  
 ClassDef(AliAnalysisTaskTree_MCut,1);
};

#endif
