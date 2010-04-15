#ifndef ALIANALYSISTASKMUONTREEBUILDER_H
#define ALIANALYSISTASKMUONTREEBUILDER_H

#include "AliAnalysisTaskSE.h"
#include "TMath.h"

//	Analysis task for muon-dimuon analysis
//	Works for real and MC events
//	author: L. Bianchi - Universita' & INFN Torino

class TH1I;
class TParticle ;
class TLorentzVector ;
//class TVector3 ;
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
  
  Double_t           	fNevt            ;    // event counter
  Double_t		fBeamEnergy      ;    // Energy of the beam (required for the CS angle)
  TList 		*fOutput         ;    // output
  TTree 		*fOutputTree     ;    //! tree output
  
  Bool_t	fIsMC;
  
  Bool_t	fIsSelected;
  char		fTrigClass[100];

  Int_t		fNumMuonTracks	;		// variables for single mu
  Int_t		fNumSPDTracklets;
  Int_t		fNumContributors;
  Double_t	fVertex[3];
  Double_t	fpT[10];
  Double_t	fE[10];
  Double_t	fpx[10];
  Double_t	fpy[10];
  Double_t	fpz[10];
  Double_t	fpxUncorr[10];
  Double_t	fpyUncorr[10];
  Double_t	fpzUncorr[10];
  Double_t	fy[10];
  Double_t	feta[10];
  Double_t	fphi[10];
  Int_t		fMatchTrig[10];
  Double_t	fTrackChi2[10];
  Double_t	fMatchTrigChi2[10];
  Double_t	fDCA[10];
  Short_t	fCharge[10];
  Int_t		fMuFamily[10];
  Double_t	fRAtAbsEnd[10];

  Int_t		fNumDimuons;			// variables for dimuons
  Int_t		fDimuonConstituent[45][2];
  Double_t	fpTdimuon[45];
  Double_t	fpxdimuon[45];
  Double_t	fpydimuon[45];
  Double_t	fpzdimuon[45];
  Double_t	fydimuon[45];
  Double_t	fiMassdimuon[45];
  Double_t	fcostCS[45];
  Double_t	fcostHE[45];
  Double_t	fphiCS[45];
  Double_t	fphiHE[45];
  
  //Int_t		fIsPrimary[10];
  Int_t		fPDG[10];
  Int_t		fPDGmother[10];
  Int_t		fPDGdimu[45];



  
  
  Int_t   FindDimuFamily(AliMCParticle* mcTrack1,AliMCParticle* mcTrack2, AliMCEvent* mcEvent) const;
  Int_t	  FindMuFamily(AliMCParticle* mcTrack, AliMCEvent* mcEvent) const;
  Double_t Imass  (Double_t e1, Double_t px1, Double_t py1, Double_t pz1, Double_t e2, Double_t px2, Double_t py2, Double_t pz2) const;
  Double_t Rap	  (Double_t e, Double_t pz) const;
//   Double_t Imass(Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t) const;
//   Double_t Rap(Double_t, Double_t) const;
  Double_t Phideg(Double_t phi) const;
  
//   Double_t CostCS (Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t);
//   Double_t CostHE (Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t);
//   Double_t PhiCS  (Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t);
//   Double_t PhiHE  (Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t);
  Double_t CostCS (Double_t px1, Double_t py1, Double_t pz1, Double_t e1, Double_t charge1, Double_t px2, Double_t py2, Double_t pz2, Double_t e2);
  Double_t CostHE (Double_t px1, Double_t py1, Double_t pz1, Double_t e1, Double_t charge1, Double_t px2, Double_t py2, Double_t pz2, Double_t e2);
  Double_t PhiCS  (Double_t px1, Double_t py1, Double_t pz1, Double_t e1, Double_t charge1, Double_t px2, Double_t py2, Double_t pz2, Double_t e2);
  Double_t PhiHE  (Double_t px1, Double_t py1, Double_t pz1, Double_t e1, Double_t charge1, Double_t px2, Double_t py2, Double_t pz2, Double_t e2);
  
  ClassDef(AliAnalysisTaskMuonTreeBuilder,1);
};

#endif
