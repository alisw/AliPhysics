#ifndef ALIANALYSISTASKMUONTREEBUILDER_H
#define ALIANALYSISTASKMUONTREEBUILDER_H

/* $Id$ */ 

#include "AliAnalysisTaskSE.h"
#include "TMath.h"

//	Analysis task for muon-dimuon analysis
//	Works for real and MC events
//	author: L. Bianchi - Universita' & INFN Torino

class TH1I;
class TParticle ;
class TLorentzVector ;
class TFile ;
class AliStack ;
class AliESDtrack;
class AliVParticle;
class AliMCParticle;
class AliMCEvent;


class AliAnalysisTaskMuonTreeBuilder : public AliAnalysisTaskSE {
  public:

  AliAnalysisTaskMuonTreeBuilder();
  AliAnalysisTaskMuonTreeBuilder(const Char_t* name);
  AliAnalysisTaskMuonTreeBuilder& operator= (const AliAnalysisTaskMuonTreeBuilder& c);
  AliAnalysisTaskMuonTreeBuilder(const AliAnalysisTaskMuonTreeBuilder& c);
  virtual ~AliAnalysisTaskMuonTreeBuilder();

  // ANALYSIS FRAMEWORK STUFF to loop on data and fill output objects
  void     UserExec(Option_t *option);
  void     Terminate(Option_t *);
  void     UserCreateOutputObjects();
  

  // Data types
  void   SetIsMC     (Bool_t flagMC)  	{fIsMC=flagMC;}
  void	 SetBeamEnergy   (Double_t en)		{fBeamEnergy=en;}
  
 protected:
  
  Double_t           	fNevt            ;    	// event counter
  Double_t		fBeamEnergy      ;    	// Energy of the beam (required for the CS angle)
  TList 		*fOutput         ;    	// output
  TTree 		*fOutputTree     ;    	//! tree output
  
  Bool_t	fIsMC;				// if MC truth has to be read
  
  Bool_t	fIsSelected;			// physics selection flag
  char		fTrigClass[100];		// fired trigger classes

  Int_t		fNumMuonTracks	;		// muon tracks in the event
  Int_t		fNumSPDTracklets;		// spd tracklets
  Int_t		fNumContributors;		// n contributors
  Double_t	fVertex[3];			// x,y,z vertex
  Double_t	fpT[10];			// single mu pT
  Double_t	fE[10];				// single mu E
  Double_t	fpx[10];			// single mu px
  Double_t	fpy[10];			// single mu py
  Double_t	fpz[10];			// single mu pz
  Double_t	fpxUncorr[10];			// single mu px uncorrected
  Double_t	fpyUncorr[10];			// single mu py uncorrected
  Double_t	fpzUncorr[10];			// single mu pz uncorrected
  Double_t	fy[10];				// single mu y
  Double_t	feta[10];			// single mu eta
  Double_t	fphi[10];			// single mu phi
  Int_t		fMatchTrig[10];			// single mu match trigger
  Double_t	fTrackChi2[10];			// single mu chi2 track
  Double_t	fMatchTrigChi2[10];		// single mu chi2 of match trigger
  Double_t	fDCA[10];			// single mu DCA
  Short_t	fCharge[10];			// single mu charge
  Int_t		fMuFamily[10];			// single mu provenience
  Double_t	fRAtAbsEnd[10];			// single mu distance from beam center at end abs

  Int_t		fNumDimuons;			// dimuons in the event
  Int_t		fDimuonConstituent[45][2];	// reference to single mus
  Double_t	fpTdimuon[45];			// dimuon pT
  Double_t	fpxdimuon[45];			// dimuon px
  Double_t	fpydimuon[45];			// dimuon py
  Double_t	fpzdimuon[45];			// dimuon pz
  Double_t	fydimuon[45];			// dimuon y
  Double_t	fiMassdimuon[45];		// dimuon invariant mass
  Double_t	fcostCS[45];			// dimuon cos theta Collins-Soper
  Double_t	fcostHE[45];			// dimuon cos theta Helicity
  Double_t	fphiCS[45];			// dimuon phi Collins-Soper
  Double_t	fphiHE[45];			// dimuon phi Helicity
  
  //Int_t		fIsPrimary[10];
  Int_t		fPDG[10];			// PDG single mu
  Int_t		fPDGmother[10];			// PDG mother single mu
  Int_t		fPDGdimu[45];			// PDG dimu



  
  
  Int_t   FindDimuFamily(AliMCParticle* mcTrack1,AliMCParticle* mcTrack2, AliMCEvent* mcEvent) const;
  Int_t	  FindMuFamily(AliMCParticle* mcTrack, AliMCEvent* mcEvent) const;
  Double_t Imass  (Double_t e1, Double_t px1, Double_t py1, Double_t pz1, Double_t e2, Double_t px2, Double_t py2, Double_t pz2) const;
  Double_t Rap	  (Double_t e, Double_t pz) const;
  Double_t Phideg(Double_t phi) const;
  
  Double_t CostCS (Double_t px1, Double_t py1, Double_t pz1, Double_t e1, Double_t charge1, Double_t px2, Double_t py2, Double_t pz2, Double_t e2);
  Double_t CostHE (Double_t px1, Double_t py1, Double_t pz1, Double_t e1, Double_t charge1, Double_t px2, Double_t py2, Double_t pz2, Double_t e2);
  Double_t PhiCS  (Double_t px1, Double_t py1, Double_t pz1, Double_t e1, Double_t charge1, Double_t px2, Double_t py2, Double_t pz2, Double_t e2);
  Double_t PhiHE  (Double_t px1, Double_t py1, Double_t pz1, Double_t e1, Double_t charge1, Double_t px2, Double_t py2, Double_t pz2, Double_t e2);
  
  ClassDef(AliAnalysisTaskMuonTreeBuilder,1);
};

#endif
