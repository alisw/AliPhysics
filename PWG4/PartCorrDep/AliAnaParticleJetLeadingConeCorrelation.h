#ifndef ALIANAPARTICLEJETLEADINGCONECORRELATION_H
#define ALIANAPARTICLEJETLEADINGCONECORRELATION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id: $ */

//_________________________________________________________________________
// Class that contains the algorithm for the reconstruction of jet, cone around leading particle
// The seed is a backward particle (direct photon)
// 1) Take the a trigger particle found stored in AliAODPWG4ParticleCorrelation,
// 2) Search for the highest pt leading particle opposite to the trigger within a phi, pt window
// 3) Take all particles around leading in a cone R with pt larger than threshold and construct the jet
//
//  Class created from old AliPHOSGammaJet
//  (see AliRoot versions previous Release 4-09)
//
//-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
class TH2F;

//---- Analysis system ----
class AliAODTrack;
class AliVCluster;
class AliCaloTrackReader;
class AliNeutralMesonSelection;

#include "AliAnaPartCorrBaseClass.h"

class AliAnaParticleJetLeadingConeCorrelation : public AliAnaPartCorrBaseClass {

public: 
  AliAnaParticleJetLeadingConeCorrelation() ; // default ctor
  virtual ~AliAnaParticleJetLeadingConeCorrelation() ; //virtual dtor

  TList * GetCreateOutputObjects();

  void InitParameters();
  
  void Print(const Option_t * opt) const;
  
  Bool_t AreJetsRecalculated() const {return fReMakeJet ; } 
  void SwitchOnJetsRecalculation(){fReMakeJet = kTRUE; }
  void SwitchOffJetsRecalculation(){fReMakeJet = kFALSE; }
  
  Bool_t AreJetsOnlyInCTS() const {return fJetsOnlyInCTS ; } 
  void SwitchOnJetsOnlyInCTS(){fJetsOnlyInCTS = kTRUE; }
  void SwitchOffJetsOnlyInCTS(){fJetsOnlyInCTS = kFALSE; }
  
  Bool_t AreSeveralConeAndPtCuts() const {return fSeveralConeAndPtCuts ; }
  void SwitchOnSeveralConeAndPtCuts(){fSeveralConeAndPtCuts = kTRUE ;}
  void SwitchOffSeveralConeAndPtCuts(){fSeveralConeAndPtCuts = kFALSE ;}
  
  Bool_t IsPbPb() const {return fPbPb ; }
  void SetppCollisions(){fPbPb = kFALSE; }
  void SetPbPbCollisions(){fPbPb = kTRUE; }
  
  Double_t GetDeltaPhiMaxCut() const {return fDeltaPhiMaxCut ; }
  Double_t GetDeltaPhiMinCut() const {return fDeltaPhiMinCut ; }
  Double_t GetLeadingRatioMaxCut() const {return fLeadingRatioMaxCut ; }
  Double_t GetLeadingRatioMinCut() const {return fLeadingRatioMinCut ; }

  Double_t GetPtTriggerSelectionCut() const {return fPtTriggerSelectionCut ; }
  Double_t GetJetRatioMaxCut() const {return fJetRatioMaxCut ; }
  Double_t GetJetRatioMinCut() const {return fJetRatioMinCut ; }  

  void SetPtTriggerSelectionCut(Double_t cut){fPtTriggerSelectionCut = cut; }
  void SetJetSelectionMode(UInt_t select){ fSelect= select ; }
  
  Int_t     GetJetNCones()             const {return fJetNCone ; }
  Int_t     GetJetNPtThres()           const {return fJetNPt ; }
  Float_t   GetJetCone()               const {return fJetCone ; }
  Float_t   GetJetPtThreshold()        const {return fJetPtThreshold ; }
  Float_t   GetJetPtThresPbPb()        const {return fJetPtThresPbPb ; }
  Float_t   GetJetCones(Int_t i)       const {return fJetCones[i] ; }
  Float_t   GetJetPtThreshold(Int_t i) const {return fJetPtThres[i] ; }
  TString   GetJetConeName(Int_t i)    const {return fJetNameCones[i] ; }
  TString   GetJetPtThresName(Int_t i) const {return fJetNamePtThres[i] ; }
  

  void SetDeltaPhiCutRange(Double_t phimin, Double_t phimax)
  {fDeltaPhiMaxCut =phimax;  fDeltaPhiMinCut =phimin;}
  void SetLeadingRatioCutRange(Double_t ratiomin, Double_t ratiomax)
  {fLeadingRatioMaxCut =ratiomax;  fLeadingRatioMinCut = ratiomin ; }

  void SetJetNCones(Int_t n){fJetNCone = n ; }
  void SetJetNPtThresholds(Int_t n){fJetNPt = n ; }
  void SetJetCones(Int_t i, Float_t cone, TString sc) {fJetCones[i] = cone ; fJetNameCones[i] = sc; };
  void SetCone(Float_t cone) {fJetCone = cone; }
  void SetJetPtThreshold(Float_t pt){fJetPtThreshold = pt; };
  void SetJetPtThresPbPb(Float_t pt){fJetPtThresPbPb = pt; };
  void SetJetPtThresholds(Int_t i,Float_t pt, TString spt){fJetPtThres[i] = pt ; fJetNamePtThres[i] = spt; };
  
  void SetJetRatioCutRange(Double_t ratiomin, Double_t ratiomax)
  {fJetRatioMaxCut =ratiomax;  fJetRatioMinCut = ratiomin ; }
  void SetJetCTSRatioCutRange(Double_t ratiomin, Double_t ratiomax)
  {fJetCTSRatioMaxCut =ratiomax;  fJetCTSRatioMinCut = ratiomin ; }
  
  Bool_t OnlyIsolated() const {return fSelectIsolated ; }
  void SelectIsolated(Bool_t select) {fSelectIsolated = select ; }
    
 private:
  
  Double_t CalculateJetRatioLimit(const Double_t ptTrig, const Double_t *param, const Double_t *x) const ;
  
  void FillJetHistos(AliAODPWG4ParticleCorrelation * particle, const TLorentzVector  leading, const TLorentzVector jet, const TString type, const TString lastname);
  
  TList * GetOutputContainer() const {return fOutCont; }
  
  Bool_t IsJetSelected(const Double_t ptTrig, const Double_t ptjet) const ;
  Bool_t IsParticleInJetCone(const Double_t eta, Double_t phi, const Double_t etal, Double_t phil) const ;
  
  void GetLeadingCharge(AliAODPWG4ParticleCorrelation* const particle, TLorentzVector & pLeading) const ;
  void GetLeadingPi0   (AliAODPWG4ParticleCorrelation* const particle, TLorentzVector & pLeading) ;
  Bool_t GetLeadingParticle(AliAODPWG4ParticleCorrelation *particle, TLorentzVector &  pLeading)  ;
  
  void MakeAnalysisFillAOD();
  void MakeAnalysisFillHistograms();   
  void MakeAODJet(AliAODPWG4ParticleCorrelation * particle, const TLorentzVector pLeading) const ; 
  void MakeJetFromAOD(AliAODPWG4ParticleCorrelation * particle, const TLorentzVector pLeading, 
		      TLorentzVector & jet, TLorentzVector & bkg) const ; 
  
  Bool_t  SelectCluster(AliVCluster * calo, Double_t *vertex, TLorentzVector & mom, Int_t & pdg) ;
  
 private:
  
  Bool_t     fJetsOnlyInCTS ;    // Jets measured only in TPC+ITS.
  Bool_t     fPbPb;          // PbPb event
  Bool_t     fSeveralConeAndPtCuts;     //  To play with the jet cone size and pt th.
  Bool_t     fReMakeJet ; //Re make the jet reconstruction from AODParticleCorrelation input

  //Leading particle selection parameters  
  Double_t   fDeltaPhiMaxCut ;      // Minimum Delta Phi Gamma-Leading
  Double_t   fDeltaPhiMinCut ;      //  Maximum Delta Phi Gamma-Leading
  Double_t   fLeadingRatioMaxCut ; // Leading /gamma Ratio cut maximum
  Double_t   fLeadingRatioMinCut ; // Leading/gamma Ratio cut minimum

  //Jet selection parameters
  //Fixed cuts (old)
  Double_t   fJetCTSRatioMaxCut ; // Jet(CTS) /gamma Ratio cut maximum
  Double_t   fJetCTSRatioMinCut ; // Jet(CTS) /gamma Ratio cut maximum
  Double_t   fJetRatioMaxCut ; // Jet(EMCAL+CTS)/gamma Ratio cut maximum
  Double_t   fJetRatioMinCut ; // Jet(EMCAL+CTS)/gamma Ratio cut minimum
  
  //Cuts depending on jet pt
  Double_t   fJetE1[2];    //Rec. jet energy parameters
  Double_t   fJetE2[2];    //Rec. jet energy parameters
  Double_t   fJetSigma1[2];//Rec. sigma of jet energy  parameters
  Double_t   fJetSigma2[2];//Rec. sigma of jet energy  parameters
  Double_t   fBkgMean[6];  //Background mean energy 
  Double_t   fBkgRMS[6];   //Background RMS
  Double_t   fJetXMin1[6]; //X Factor to set jet min limit for pp
  Double_t   fJetXMin2[6]; //X Factor to set jet min limit for PbPb
  Double_t   fJetXMax1[6]; //X Factor to set jet max limit for pp
  Double_t   fJetXMax2[6]; //X Factor to set jet max limit for PbPb
  
  Int_t      fJetNCone ;            // Number of jet cones sizes, maximum 5
  Int_t      fJetNPt   ;            // Number of jet particle pT threshold, maximum 5
  Double_t   fJetCone  ;            // Jet cone sizes under study (!fSeveralConeAndPtCuts)
  Double_t   fJetCones[5];         // Jet cone sizes under study (fSeveralConeAndPtCuts)
  TString    fJetNameCones[5];     // String name of cone to append to histos
  Double_t   fJetPtThreshold;       // Jet pT threshold under study(!fSeveralConeAndPtCuts)
  Double_t   fJetPtThresPbPb;       // Jet pT threshold under study(!fSeveralConeAndPtCuts)
  Double_t   fJetPtThres[5];       // Jet pT threshold under study(fSeveralConeAndPtCuts)
  TString    fJetNamePtThres[5];   // String name of pt th to append to histos
  Double_t   fPtTriggerSelectionCut; // Jet pt to change to low pt jets analysis
  UInt_t     fSelect  ;   //kTRUE: Selects all jets, no limits.
  Bool_t     fSelectIsolated ;      // Select only trigger particles isolated

  //Histograms
  //Leading particle distributions
  TList *  fOutCont ; //! Container for histograms

  TH2F * fhChargedLeadingPt  ;    //! Pt(Pt trigger) distribution of charged hadrons
  TH2F * fhChargedLeadingPhi  ;   //! Phi(Pt trigger) distribution of charged hadrons
  TH2F * fhChargedLeadingEta  ;   //! Eta(Pt trigger) distribution of charged hadrons
  TH2F * fhChargedLeadingDeltaPt  ;   //! Difference of charged hadron and trigger  pT as function of trigger p
  TH2F * fhChargedLeadingDeltaPhi  ;  //! Difference of charged hadron and trigger  phi as function of trigger pT
  TH2F * fhChargedLeadingDeltaEta ;   //! Difference of charged particle and trigger eta as function of trigger pT
  TH2F * fhChargedLeadingRatioPt  ; //! Ratio of Pt leading charge and trigger

  TH2F * fhNeutralLeadingPt   ;   //! Pt(Pt trigger) distribution of neutral hadrons
  TH2F * fhNeutralLeadingPhi   ;  //! Phi(Pt trigger) distribution of neutral hadrons
  TH2F * fhNeutralLeadingEta   ;  //! Eta(Pt trigger) distribution of neutral hadrons
  TH2F * fhNeutralLeadingDeltaPt   ;  //! Difference of neutral hadron and trigger pT as function of trigger pT
  TH2F * fhNeutralLeadingDeltaPhi  ;  //! Difference of neutral hadron and trigger phi as function of trigger pT
  TH2F * fhNeutralLeadingDeltaEta ;   //! Difference of charged particle and trigger eta as function of trigger pT
  TH2F * fhNeutralLeadingRatioPt   ;  //! Ratio of Pt leading neutral and trigger
	
  TH2F * fhChargedLeadingXi  ;   //! Ln (pt leading charge / pt trigger)
  TH2F * fhNeutralLeadingXi  ;   //! Ln (pt leading neutral / pt trigger)
	
  TH2F * fhChargedLeadingDeltaPhiRatioPt30  ;  //! Difference of charged hadron and trigger  phi as function of pT leading / trigger pT, pT Trigger > 30 GeV
  TH2F * fhNeutralLeadingDeltaPhiRatioPt30  ;  //! Difference of neutral hadron and trigger  phi as function of pT leading / trigger pT, pT Trigger > 30 GeV
  TH2F * fhChargedLeadingDeltaPhiRatioPt50  ;  //! Difference of charged hadron and trigger  phi as function of pT leading / trigger pT, pT Trigger > 50 GeV
  TH2F * fhNeutralLeadingDeltaPhiRatioPt50  ;  //! Difference of neutral hadron and trigger  phi as function of pT leading / trigger pT, pT Trigger > 50 GeV
	
  // Jet distributions
  // Fixed cone and pt threshold
  TH2F * fhJetPt  ; //! leading pt jet vs pt trigger
  TH2F * fhJetRatioPt  ; //! Ratio of pt jet and pt trigger
  TH2F * fhJetDeltaPhi  ; //! Delta phi jet-trigger
  TH2F * fhJetDeltaEta  ; //! Delta eta jet-trigger
  TH2F * fhJetLeadingRatioPt  ; //! Ratio of pt leading and pt jet
  TH2F * fhJetLeadingDeltaPhi  ; //! Delta phi jet-leading
  TH2F * fhJetLeadingDeltaEta  ; //! Delta eta jet-leading
  TH2F * fhJetFFz; //! Accepted reconstructed jet fragmentation function, z=ptjet/pttrig
  TH2F * fhJetFFxi; //! Accepted reconstructed jet fragmentation function, xsi = ln(pttrig/ptjet)
  TH2F * fhJetFFpt; //! Jet particle pt distribution in cone
  TH2F * fhJetNTracksInCone   ; //! jet multiplicity in cone

  TH2F * fhBkgPt  ; //! leading pt bakground vs pt trigger
  TH2F * fhBkgRatioPt  ; //! Ratio of pt background and pt trigger
  TH2F * fhBkgDeltaPhi  ; //! Delta phi background-trigger
  TH2F * fhBkgDeltaEta  ; //! Delta eta background-trigger
  TH2F * fhBkgLeadingRatioPt  ; //! Ratio of pt leading and pt background
  TH2F * fhBkgLeadingDeltaPhi  ; //! Delta phi background-leading
  TH2F * fhBkgLeadingDeltaEta  ; //! Delta eta background-leading
  TH2F * fhBkgFFz; //! Accepted reconstructed background fragmentation function, z=ptjet/pttrig
  TH2F * fhBkgFFxi; //! Accepted reconstructed background fragmentation function, xsi = ln(pttrig/ptjet)
  TH2F * fhBkgFFpt; //! Background particle pt distribution in cone
  TH2F * fhBkgNTracksInCone   ; //! Background multiplicity in cone

  // Variable cone and pt threshold

  TH2F * fhJetPts[5][5]; //! leading pt jet vs pt trigger
  TH2F * fhJetRatioPts[5][5]; //! Ratio of pt jet and pt trigger
  TH2F * fhJetDeltaPhis[5][5]; //! Delta phi jet-trigger
  TH2F * fhJetDeltaEtas[5][5]; //! Delta eta jet-trigger
  TH2F * fhJetLeadingRatioPts[5][5]; //! Ratio of pt leading and pt jet
  TH2F * fhJetLeadingDeltaPhis[5][5]; //! Delta phi jet-leading
  TH2F * fhJetLeadingDeltaEtas[5][5]; //! Delta eta jet-leading
  TH2F * fhJetFFzs[5][5]; //! Accepted reconstructed jet fragmentation function, z=ptjet/pttrig
  TH2F * fhJetFFxis[5][5]; //! Accepted reconstructed jet fragmentation function, xsi = ln(pttrig/ptjet)
  TH2F * fhJetFFpts[5][5]; //! Jet particle pt distribution in cone
  TH2F * fhJetNTracksInCones[5][5]; //! jet multiplicity in cone

  TH2F * fhBkgPts[5][5]; //! leading pt bakground vs pt trigger
  TH2F * fhBkgRatioPts[5][5]; //! Ratio of pt background and pt trigger
  TH2F * fhBkgDeltaPhis[5][5]; //! Delta phi background-trigger
  TH2F * fhBkgDeltaEtas[5][5]; //! Delta eta background-trigger
  TH2F * fhBkgLeadingRatioPts[5][5]; //! Ratio of pt leading and pt background
  TH2F * fhBkgLeadingDeltaPhis[5][5]; //! Delta phi background-leading
  TH2F * fhBkgLeadingDeltaEtas[5][5]; //! Delta eta background-leading
  TH2F * fhBkgFFzs[5][5]; //! Accepted reconstructed background fragmentation function, z=ptjet/pttrig
  TH2F * fhBkgFFxis[5][5]; //! Accepted reconstructed background fragmentation function, xsi = ln(pttrig/ptjet)
  TH2F * fhBkgFFpts[5][5]; //! Background particle pt distribution in cone
  TH2F * fhBkgNTracksInCones[5][5]; //! Background multiplicity in cone
  
  AliAnaParticleJetLeadingConeCorrelation(const AliAnaParticleJetLeadingConeCorrelation & g) ; // cpy ctor
  AliAnaParticleJetLeadingConeCorrelation & operator = (const AliAnaParticleJetLeadingConeCorrelation & g) ;//cpy assignment
  
  ClassDef(AliAnaParticleJetLeadingConeCorrelation,1)
 } ;
 

#endif //ALIANAPARTICLEJETLEADINGCONECORRELATION_H



