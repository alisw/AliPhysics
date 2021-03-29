#ifndef AliAnalysisTaskQuarkoniumTreeEmbedding_H
#define AliAnalysisTaskQuarkoniumTreeEmbedding_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"

class TObjArray;
class AliVParticle;
class AliAODEvent;
class TLorentzVector;
class AliMuonTrackCuts;

class AliAnalysisTaskQuarkoniumTreeEmbedding : public AliAnalysisTaskSE {
  public:

  AliAnalysisTaskQuarkoniumTreeEmbedding();
  AliAnalysisTaskQuarkoniumTreeEmbedding(const char *name);
  virtual ~AliAnalysisTaskQuarkoniumTreeEmbedding();

  void UserCreateOutputObjects();
  void UserExec(Option_t *option);
  void Terminate(Option_t *);
  virtual void NotifyRun();

  void SetBeamEnergy(Double_t en) {fBeamEnergy=en;}
  void SetResonance(TString resonance) {fResonance=resonance;}
  void SetAnalysisType(const char* type) {fkAnalysisType=type;}
  void SetPeriod(TString period) {fPeriod=period;}
  AliMuonTrackCuts* fMuonTrackCuts;

 private:
  AliAnalysisTaskQuarkoniumTreeEmbedding(const AliAnalysisTaskQuarkoniumTreeEmbedding&);
  AliAnalysisTaskQuarkoniumTreeEmbedding& operator=(const AliAnalysisTaskQuarkoniumTreeEmbedding&);

 //protected:

  TTree         *fOutputTree	;      //! tree output
  Double_t fBeamEnergy;   // Energy of the beam (required for the CS angle)
  const char* fkAnalysisType; //ESD or AOD based analysis
  TString fPeriod; //period
  TString fResonance; //resonance

  Int_t		fNMuons_gen;		// muon tracks in the event
  Int_t		fNDimu_gen;			// dimuons in the event
  Int_t		fNMuons_rec;		// muon tracks in the event
  Int_t		fNDimu_rec;			// dimuons in the event
  Double_t      fPercentV0M;

  Double_t	fPt_rec[100];		 // single mu pT
  Double_t	fE_rec[100];			// single mu E
  Double_t	fPx_rec[100];		// single mu px
  Double_t	fPy_rec[100];		// single mu py
  Double_t	fPz_rec[100];		// single mu pz
  Double_t	fY_rec[100];		// single mu y
  Double_t	fEta_rec[100];		// single mu eta
  Int_t		  fMatchTrig_rec[100];		// single mu match trigger
  Double_t	fTrackChi2_rec[100];		// single mu chi2 track
  Double_t	fMatchTrigChi2_rec[100];	// single mu chi2 of match trigger
  Int_t	    fCharge_rec[100];		// single mu charge
  Double_t	fRAtAbsEnd_rec[100];		// single mu distance from beam center at end abs
  Int_t	    fpDCA_rec[100];		// single mu charge

  Double_t	fDimuPt_gen[1000];			    // dimuon pT
  Double_t	fDimuPx_gen[1000]; 		    // dimuon px
  Double_t	fDimuPy_gen[1000]; 		    // dimuon py
  Double_t	fDimuPz_gen[1000]; 		    // dimuon pz
  Double_t	fDimuY_gen[1000];			// dimuon y
  Double_t	fDimuMass_gen[1000];		// dimuon invariant mass
  Int_t	    fDimuCharge_gen[1000];		// dimuon charge

  Int_t		  fDimuMu_rec[1000][2];	// reference to single mus
  Double_t	fDimuPt_rec[1000];			    // dimuon pT
  Double_t	fDimuPx_rec[1000]; 		    // dimuon px
  Double_t	fDimuPy_rec[1000]; 		    // dimuon py
  Double_t	fDimuPz_rec[1000]; 		    // dimuon pz
  Double_t	fDimuY_rec[1000];			// dimuon y
  Double_t	fDimuMass_rec[1000];		// dimuon invariant mass
  Int_t	        fDimuCharge_rec[1000];		// dimuon charge
  Int_t	        fDimuMatch_rec[1000];		// dimuon match

  AliAODEvent* fAODEvent;      //! AOD event  //tolgo !

 ClassDef(AliAnalysisTaskQuarkoniumTreeEmbedding,1);
};

#endif
