#ifndef ALIANAGAMMAJETLEADCONE_H
#define ALIANAGAMMAJETLEADCONE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.2  2007/08/17 12:40:04  schutz
 * New analysis classes by Gustavo Conesa
 *
 * Revision 1.1.2.1  2007/07/26 10:32:09  schutz
 * new analysis classes in the the new analysis framework
 *
 *
 */

//_________________________________________________________________________
// Class that contains the algorithm for the reconstruction of jet, cone around leading particle
// 1)Take the prompt photon found with AliAnaGammaDirect,
// 2) Search for the highest pt leading particle opposite to the photon within a phi, pt window
// 3) Take all particles around leading in a cone R with pt larger than threshold and construct the jet
//
//  Class created from old AliPHOSGammaJet
//  (see AliRoot versions previous Release 4-09)
//-- Author: Gustavo Conesa (INFN-LNF)

#include "AliAnaGammaCorrelation.h"

class AliAnaGammaJetLeadCone : public AliAnaGammaCorrelation {

public: 
  
  AliAnaGammaJetLeadCone() ; // default ctor
  AliAnaGammaJetLeadCone(const AliAnaGammaJetLeadCone & g) ; // cpy ctor
  AliAnaGammaJetLeadCone & operator = (const AliAnaGammaJetLeadCone & g) ;//cpy assignment
  virtual ~AliAnaGammaJetLeadCone() ; //virtual dtor
  
  TList * GetCreateOutputObjects();

  void InitParameters();

  void Print(const Option_t * opt) const;
 
  Bool_t   AreSeveralConeAndPtCuts() const {return fSeveralConeAndPtCuts ; }
  void SetSeveralConeAndPtCuts(Bool_t several){fSeveralConeAndPtCuts = several ;}

  Bool_t   IsPbPb() const {return fPbPb ; }
  void SetPbPb(Bool_t opt){fPbPb = opt; }

  Double_t GetEtaEMCALCut() const {return fEtaEMCALCut;}
  Double_t GetPhiEMCALCut(Int_t i) const {return fPhiEMCALCut[i];}
  void SetEtaEMCALCut(Double_t etacut) {fEtaEMCALCut = etacut ; }
  void SetPhiEMCALCut(Double_t phi, Int_t i){fPhiEMCALCut[i] = phi; }

  Double_t GetPtJetSelectionCut() const {return fPtJetSelectionCut ; }
  Double_t GetJetRatioMaxCut() const {return fJetRatioMaxCut ; }
  Double_t GetJetRatioMinCut() const {return fJetRatioMinCut ; }
 
  void SetPtJetSelectionCut(Double_t cut){fPtJetSelectionCut = cut; }
  void SetJetSelection(UInt_t select){ fSelect= select ; }

  Int_t       GetJetNCones() const {return fJetNCone ; }
  Int_t       GetJetNPtThres() const {return fJetNPt ; }
  Float_t    GetJetCone() const {return fJetCone ; }
  Float_t    GetJetPtThreshold() const {return fJetPtThreshold ; }
  Float_t    GetJetPtThresPbPb() const {return fJetPtThresPbPb ; }
  Float_t    GetJetCones(Int_t i) const {return fJetCones[i] ; }
  Float_t    GetJetPtThreshold(Int_t i) const {return fJetPtThres[i] ; }
  TString   GetJetConeName(Int_t i) const {return fJetNameCones[i] ; }
  TString   GetJetPtThresName(Int_t i) const {return fJetNamePtThres[i] ; }

  void SetJetNCones(Int_t n){fJetNCone = n ; }
  void SetJetNPtThresholds(Int_t n){fJetNPt = n ; }
  void SetJetCones(Int_t i, Float_t cone, TString sc)
    {fJetCones[i] = cone ; fJetNameCones[i] = sc; };
  void SetCone(Float_t cone)
    {fJetCone = cone; }
  void SetJetPtThreshold(Float_t pt){fJetPtThreshold = pt; };
  void SetJetPtThresPbPb(Float_t pt){fJetPtThresPbPb = pt; };
  void SetJetPtThresholds(Int_t i,Float_t pt, TString spt){fJetPtThres[i] = pt ; 
  fJetNamePtThres[i] = spt; };

  void SetJetRatioCutRange(Double_t ratiomin, Double_t ratiomax)
    {fJetRatioMaxCut =ratiomax;  fJetRatioMinCut = ratiomin ; }
  void SetJetCTSRatioCutRange(Double_t ratiomin, Double_t ratiomax)
    {fJetCTSRatioMaxCut =ratiomax;  fJetCTSRatioMinCut = ratiomin ; }
  
  void MakeGammaCorrelation(TParticle * pGamma, TClonesArray * plCTS,   TClonesArray * plNe) ;
  
  private:

  Double_t CalculateJetRatioLimit(const Double_t ptg, const Double_t *param, 
				  const Double_t *x);
  void FillJetHistos(TClonesArray * pl, Double_t ptg, Double_t ptl,TString type, TString lastname);

  Bool_t IsJetSelected(const Double_t ptg, const Double_t ptjet);

  void MakeJet(TClonesArray * plCTS, TClonesArray * plNe, 
	       TParticle * pGamma, TParticle* pLeading, TString lastname); 
  
  void GetLeadingCharge(TClonesArray * pl, TParticle *pGamma, TParticle * pLeading) const ;
  void GetLeadingPi0   (TClonesArray * pl, TParticle *pGamma, TParticle * pLeading)  ;
  Bool_t GetLeadingParticle(TClonesArray * plCTS, TClonesArray * plNe, 
			    TParticle *pGamma, TParticle * pLeading) ;
  
  void SetJet(TParticle * part, Bool_t & b, Float_t cone, Double_t eta, 
	      Double_t phi);

  private:

  Bool_t       fPbPb;          // PbPb event
  Bool_t       fSeveralConeAndPtCuts;     //  To play with the jet cone size and pt th.

  Double_t   fEtaEMCALCut ;         // Eta EMCAL acceptance
  Double_t   fPhiEMCALCut[2] ; // Phi cut maximum

  //Jet selection parameters
  //Fixed cuts (old)
  Double_t   fJetCTSRatioMaxCut ; // Leading particle/gamma Ratio cut maximum
  Double_t   fJetCTSRatioMinCut ; // Leading particle/gamma Ratio cut minimum
  Double_t   fJetRatioMaxCut ; // Jet/gamma Ratio cut maximum
  Double_t   fJetRatioMinCut ; // Jet/gamma Ratio cut minimum

  //Cuts depending on jet pt
  Double_t fJetE1[2];    //Rec. jet energy parameters
  Double_t fJetE2[2];    //Rec. jet energy parameters
  Double_t fJetSigma1[2];//Rec. sigma of jet energy  parameters
  Double_t fJetSigma2[2];//Rec. sigma of jet energy  parameters
  Double_t fBkgMean[6];  //Background mean energy 
  Double_t fBkgRMS[6];   //Background RMS
  Double_t fJetXMin1[6]; //X Factor to set jet min limit for pp
  Double_t fJetXMin2[6]; //X Factor to set jet min limit for PbPb
  Double_t fJetXMax1[6]; //X Factor to set jet max limit for pp
  Double_t fJetXMax2[6]; //X Factor to set jet max limit for PbPb

  Int_t         fJetNCone ;            // Number of jet cones sizes
  Int_t         fJetNPt   ;            // Number of jet particle pT threshold
  Double_t   fJetCone  ;            // Jet cone sizes under study (!fSeveralConeAndPtCuts)
  Double_t   fJetCones[10];         // Jet cone sizes under study (fSeveralConeAndPtCuts)
  TString     fJetNameCones[10];     // String name of cone to append to histos
  Double_t   fJetPtThreshold;       // Jet pT threshold under study(!fSeveralConeAndPtCuts)
  Double_t   fJetPtThresPbPb;       // Jet pT threshold under study(!fSeveralConeAndPtCuts)
  Double_t   fJetPtThres[10];       // Jet pT threshold under study(fSeveralConeAndPtCuts)
  TString     fJetNamePtThres[10];   // String name of pt th to append to histos
  Double_t   fPtJetSelectionCut; // Jet pt to change to low pt jets analysis
  UInt_t       fSelect  ;   //kTRUE: Selects all jets, no limits.
 
  //Histograms
  //Particle distributions
  TH2F * fhPhiCharged  ; //Phi distribution of charged particles
  TH2F * fhPhiNeutral   ;  //Phi distribution of neutral particles
  TH2F * fhEtaCharged  ;  //Eta distribution of charged particles
  TH2F * fhEtaNeutral   ;  //Eta distribution of neutral particles
  //Leading particle distributions
  TH2F * fhDeltaPhiGammaCharged  ;   //Difference of charged particle phi and prompt gamma phi as function of gamma pT
  TH2F * fhDeltaPhiGammaNeutral   ;  //Difference of neutral particle phi and prompt gamma phi as function of gamma pT
  TH2F * fhDeltaEtaGammaCharged  ;  //Difference of charged particle eta and prompt gamma eta as function of gamma pT
  TH2F * fhDeltaEtaGammaNeutral  ;   //Difference of charged particle eta and prompt gamma eta as function of charged pT

  TH2F * fhAnglePairLeading  ; //Aperture angle of decay photons of leading pi0
  TH2F * fhInvMassPairLeading  ; //Invariant mass of decay photons of leading pi0
  TH2F * fhChargedRatio  ; //Ratio of leading charge and prompt gamma
  TH2F * fhNeutralRatio   ;  //Ratio of leading neutral and prompt gamma
  TH1F * fhNBkg   ; //Bakground multiplicity
  TH2F * fhNLeading  ; //Accepted leading particle pt distribution

  //Jet distributions
  //Fixed cone and pt threshold
  TH1F * fhNJet  ; //Accepted reconstructed Jet pt distribution 
  TH2F * fhJetRatio  ; //Ratio of pt jet and pt gamma
  TH2F * fhJetPt   ; //reconstructed pt jet vs prompt pt gamma
  TH2F * fhBkgRatio   ; //leading pt bakground / pt gamma 
  TH2F * fhBkgPt  ; //leading pt bakground vs pt gamma 
  TH2F * fhJetFragment  ; //Accepted reconstructed jet fragmentation function
  TH2F * fhBkgFragment  ;  //Background "fragmentation function"
  TH2F * fhJetPtDist  ; //Jet particle pt distribution
  TH2F * fhBkgPtDist  ; //Background jet particle pt distribution

  //Variable cone and pt threshold
  TH2F * fhJetRatios[5][5]; //Ratio of pt jet and pt gamma
  TH2F * fhJetPts[5][5]; //reconstructed pt jet vs prompt pt gamma
  TH2F * fhBkgRatios[5][5];  //leading pt bakground / pt gamma 
  TH2F * fhBkgPts[5][5]; //leading pt bakground vs pt gamma 
  
  TH2F * fhNLeadings[5][5];  //Accepted leading particle pt distribution
  TH1F * fhNJets[5][5]; //Accepted reconstructed Jet pt distribution 
  TH1F * fhNBkgs[5][5]; //Bakground multiplicity
  
  TH2F * fhJetFragments[5][5];//Accepted reconstructed jet fragmentation function
  TH2F * fhBkgFragments[5][5];  //Background "fragmentation function"
  TH2F * fhJetPtDists[5][5]; //Jet particle pt distribution
  TH2F * fhBkgPtDists[5][5]; //Background jet particle pt distribution
  

  ClassDef(AliAnaGammaJetLeadCone,1)
} ;
 

#endif //ALIANAGAMMAJETLEADCONE_H



