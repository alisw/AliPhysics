#ifndef AliAnalysisTaskPbPbTree_MCut_H
#define AliAnalysisTaskPbPbTree_MCut_H

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

class AliAnalysisTaskPbPbTree_MCut: public AliAnalysisTaskSE {
  public:

  AliAnalysisTaskPbPbTree_MCut();
  AliAnalysisTaskPbPbTree_MCut(const char *name);
  virtual ~AliAnalysisTaskPbPbTree_MCut();

  void UserCreateOutputObjects();
  void UserExec(Option_t *option);
  void Terminate(Option_t *);
  virtual void NotifyRun();

  void SetBeamEnergy(Double_t en) {fBeamEnergy=en;}
  void SetMassCut(Double_t MassCut) {fMassCut=MassCut;}
  void SetAnalysisType(const char* type) {fkAnalysisType=type;}
  void SetPeriod(TString period) {fPeriod=period;}

 private:
  AliAnalysisTaskPbPbTree_MCut(const AliAnalysisTaskPbPbTree_MCut&);
  AliAnalysisTaskPbPbTree_MCut& operator=(const AliAnalysisTaskPbPbTree_MCut&);

  virtual Float_t GetVertex(AliAODEvent* aod) const;
  virtual Short_t GetVertexZ(Float_t vtxZ) const;
  Double_t CalcCorPhi(Double_t phi, Double_t dPhi) const;
  void OpenInfoCalbration(Int_t run);

 //protected:

  TTree         *fOutputTree;      //! tree output
  Double_t      fNevt;             // event counter
  TH1D          *fhNEv;            //! histo

  Double_t fBeamEnergy;            // Energy of the beam (required for the CS angle)
  Double_t fMassCut;               //  Mass Cut

  const char* fkAnalysisType;      //ESD or AOD based analysis
  TString fPeriod;                 //period
  Int_t fRun;                      // run number - for calibration
  Double_t  fVtxCut;               // vertex cut = 14 cm

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
  Int_t		fNDimu;			              // dimuons in the event
  Bool_t        fIsPhysSelected;    //check on phys selection
  AliAODEvent*  fAODEvent;          //! AOD event
  char   	fTrigClass[800];	        // fired trigger classes
  UInt_t        finpmask;           // trigger input mask

  AliMuonTrackCuts* fMuonTrackCuts;
  Double_t	fVertex[3];		          // x,y,z vertex
  Float_t  fPercentV0M;             //! percentile V0
  Float_t  fPercentCL0;             //! percentile n. cluster layer 0
  Float_t  fPercentCL1;             //! percentile n. cluster layer 1
  Float_t  fPercentV0A;             //! percentile V0A
  Float_t  fPercentV0C;             //! percentile V0C
  Float_t  fPercentZNA;             //! percentile ZNA
  Float_t  fPercentZNC;             //! percentile ZNC
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
  Int_t  fpDCA[1500];               //pDCA

  Double_t	fDimuPt[400];			      // dimuon pT
  Double_t	fDimuPx[400]; 		      // dimuon px
  Double_t	fDimuPy[400]; 		      // dimuon py
  Double_t	fDimuPz[400]; 		      // dimuon pz
  Double_t	fDimuY[400];			      // dimuon y
  Double_t	fDimuMass[400];		      // dimuon invariant mass
  Int_t	        fDimuCharge[400];		// dimuon charge
  Int_t	        fDimuMatch[400];		// dimuon match

  Double_t      fPsi2Trkl;          //! Event plane Psi2
  Double_t      fPsi3Trkl;          //! Event plane Psi3
  Double_t      fPsi2RP;            //! Random plane Psi2
  Double_t      fDimuPhi[400];      //! JPsi azimuthal angle in ALICE reference frame

  Double_t      fDimuCostHE[400];   // cost Helicty frame
  Double_t      fDimuPhiHE[400];    // phi Helicty fame
  Double_t      fDimuCostCS[400];   // cost Collins-Soper
  Double_t      fDimuPhiCS[400];    // phi Collins-Soper
  Double_t      fDimuCostEPnB[400]; // cost Event-Plane (not boosted)
  Double_t      fDimuCostRPnB[400]; // cost Random-Plane (not boosted)
  Int_t		fDimuMu[400][2];	        // reference to single mus

 ClassDef(AliAnalysisTaskPbPbTree_MCut,3);
};

#endif
