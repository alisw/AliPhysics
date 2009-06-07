#ifndef ALIANACALORIMETERQA_H
#define ALIANACALORIMETERQA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id: $ */

//_________________________________________________________________________
// Class to check results from simulations or reconstructed real data. 
// Fill few histograms and do some checking plots
//
//-- Author: Gustavo Conesa (INFN-LNF)

// --- Root system ---
class TH2F;
class TH1F;
//class TH3F;

// --- Analysis system --- 
#include "AliAnaPartCorrBaseClass.h"
 
class AliAnaCalorimeterQA : public AliAnaPartCorrBaseClass {
  
 public: 
  
  AliAnaCalorimeterQA() ; // default ctor
  AliAnaCalorimeterQA(const AliAnaCalorimeterQA & g) ; // cpy ctor
  AliAnaCalorimeterQA & operator = (const AliAnaCalorimeterQA & g) ;//cpy assignment
  virtual ~AliAnaCalorimeterQA() {;} //virtual dtor
  
  TList * GetCreateOutputObjects();
  
  void Init();
  void InitParameters();
  
  void Print(const Option_t * opt) const;
    
  void MakeAnalysisFillHistograms() ; 
	
  TString GetCalorimeter() const {return fCalorimeter ;}
  void SetCalorimeter( TString calo ) {fCalorimeter = calo; }
  TString GetStyleMacro() const {return fStyleMacro ;}
  void SetStyleMacro( TString macro ) {fStyleMacro = macro; }
  
  void Terminate(TList * outputList);
  void ReadHistograms(TList * outputList); //Fill histograms with histograms in ouput list, needed in Terminate.

 private:
  
  TString fCalorimeter ; //Calorimeter selection
  TString fStyleMacro  ; //Location of macro for plots style
  
  //Histograms
  //CaloClusters 
  TH1F * fhE  ; //! E distribution, Reco
  TH1F * fhPt ; //! pT distribution, Reco
  TH1F * fhPhi; //! phi distribution, Reco 
  TH1F * fhEta; //! eta distribution, Reco 
  TH2F * fhEtaPhi  ; //! eta vs phi, Reco 
  TH1F * fhECharged  ; //! E distribution, Reco, matched with track
  TH1F * fhPtCharged ; //! pT distribution, Reco, matched with track
  TH1F * fhPhiCharged; //! phi distribution, Reco, matched with track 
  TH1F * fhEtaCharged; //! eta distribution, Reco, matched with track 
  TH2F * fhEtaPhiCharged  ; //! eta vs phi, Reco, matched with track 
  TH1F * fhEChargedNoOut  ; //! E distribution, Reco, matched with track, no outer param
  TH1F * fhPtChargedNoOut ; //! pT distribution, Reco, matched with track, no outer param
  TH1F * fhPhiChargedNoOut; //! phi distribution, Reco, matched with track, no outer param
  TH1F * fhEtaChargedNoOut; //! eta distribution, Reco, matched with track, no outer param
  TH2F * fhEtaPhiChargedNoOut  ; //! eta vs phi, Reco, matched with track, no outer param
  TH1F * fhDeltaE  ; //! MC-Reco E distribution	
  TH1F * fhDeltaPt ; //! MC-Reco pT distribution
  TH1F * fhDeltaPhi; //! MC-Reco phi distribution
  TH1F * fhDeltaEta; //! MC-Reco eta distribution
  TH1F * fhRatioE  ; //! Reco/MC E distribution	
  TH1F * fhRatioPt ; //! Reco/MC pT distribution
  TH1F * fhRatioPhi; //! Reco/MC phi distribution
  TH1F * fhRatioEta; //! Reco/MC eta distribution
  TH2F * fh2E  ; //! E distribution, Reco vs MC
  TH2F * fh2Pt ; //! pT distribution, Reco vs MC
  TH2F * fh2Phi; //! phi distribution, Reco vs MC
  TH2F * fh2Eta; //! eta distribution, Reco vs MC
  TH2F * fhIM; //! cluster pairs invariant mass
  TH2F * fhAsym; //! cluster pairs invariant mass	
  TH2F * fhNCellsPerCluster; //! N cells per cluster	
  TH1F * fhNClusters; //! Number of clusters
	
  //Calo Cells
  TH1F * fhNCells; //! Number of towers/crystals with signal
  TH1F * fhAmplitude; //! Amplitude measured in towers/crystals

  //MC  
  TH1F *fhGenGamPt  ; // pt of primary gamma
  TH1F *fhGenGamEta ; // eta of primart gamma
  TH1F *fhGenGamPhi ; // phi of primary gamma	
  TH1F *fhGenPi0Pt  ; // pt of primary pi0
  TH1F *fhGenPi0Eta ; // eta of primart pi0
  TH1F *fhGenPi0Phi ; // phi of primary pi0	
  TH1F *fhGenEtaPt  ; // pt of primary eta
  TH1F *fhGenEtaEta ; // eta of primart eta
  TH1F *fhGenEtaPhi ; // phi of primary eta
  TH1F *fhGenOmegaPt  ; // pt of primary omega
  TH1F *fhGenOmegaEta ; // eta of primart omega
  TH1F *fhGenOmegaPhi ; // phi of primary omega
  TH1F *fhGenElePt  ; // pt of primary electron
  TH1F *fhGenEleEta ; // eta of primart electron
  TH1F *fhGenElePhi ; // phi of primary electron
	
  //TH3F * fhEMVxyz    ; // Electromagnetic particle production vertex
  TH2F * fhEMVxyz    ; // Electromagnetic particle production vertex
  TH2F * fhEMR       ; // Electromagnetic distance to vertex vs rec energy  
  //TH3F * fhHaVxyz    ; // Hadron production vertex
  TH2F * fhHaVxyz    ; // Hadron production vertex
  TH2F * fhHaR       ; // Hadron distance to vertex vs rec energy  

  TH2F * fhGamE  ; //! E distribution of generated photons, Reco
  TH2F * fhGamPt ; //! pT distribution of generated photons, Reco
  TH2F * fhGamPhi; //! phi distribution of generated photon, Reco 
  TH2F * fhGamEta; //! eta distribution of generated photons, Reco 
  TH1F * fhGamDeltaE  ; //! MC-Reco E distribution of generated photons	
  TH1F * fhGamDeltaPt ; //! MC-Reco pT distribution of generated photons
  TH1F * fhGamDeltaPhi; //! MC-Reco phi distribution of generated photons
  TH1F * fhGamDeltaEta; //! MC-Reco eta distribution of generated photons
  TH1F * fhGamRatioE  ; //! Reco/MC E distribution of generated photons	
  TH1F * fhGamRatioPt ; //! Reco/MC pT distribution of generated photons
  TH1F * fhGamRatioPhi; //! Reco/MC phi distribution of generated photons
  TH1F * fhGamRatioEta; //! Reco/MC eta distribution of generated photons
  TH2F * fhGam2E  ; //! E distribution, Reco vs MC of generated photons
  TH2F * fhGam2Pt ; //! pT distribution, Reco vs MC of generated photons
  TH2F * fhGam2Phi; //! phi distribution, Reco vs MC of generated photons
  TH2F * fhGam2Eta; //! eta distribution, Reco vs MC of generated photons	
  TH2F * fhEleE  ; //! E distribution of generated electrons, Reco
  TH2F * fhElePt ; //! pT distribution of generated electrons, Reco
  TH2F * fhElePhi; //! phi distribution of generated electron, Reco 
  TH2F * fhEleEta; //! eta distribution of generated electrons, Reco 		
  TH2F * fhPi0E  ; //! E distribution of generated pi0, Reco, gamma decay overlapped
  TH2F * fhPi0Pt ; //! pT distribution of generated pi0, Reco, gamma decay overlapped
  TH2F * fhPi0Phi; //! phi distribution of generated pi0, Reco, gamma decay overlapped
  TH2F * fhPi0Eta; //! eta distribution of generated pi0, Reco, gamma decay overlapped
  TH2F * fhNeHadE  ; //! E distribution of generated neutral hadron, Reco
  TH2F * fhNeHadPt ; //! pT distribution of generated neutral hadron, Reco
  TH2F * fhNeHadPhi; //! phi distribution of generated neutral hadron, Reco 
  TH2F * fhNeHadEta; //! eta distribution of generated neutral hadron, Reco 	
  TH2F * fhChHadE  ; //! E distribution of generated charged hadron, Reco
  TH2F * fhChHadPt ; //! pT distribution of generated charged hadron, Reco
  TH2F * fhChHadPhi; //! phi distribution of generated charged hadron, Reco 
  TH2F * fhChHadEta; //! eta distribution of generated charged hadron, Reco 
	
  TH1F *fhGenGamAccE   ; // E of primary gamma
  TH1F *fhGenGamAccPt  ; // pt of primary gamma
  TH1F *fhGenGamAccEta ; // eta of primart gamma
  TH1F *fhGenGamAccPhi ; // phi of primary gamma	
  TH1F *fhGenPi0AccE   ; // E of primary pi0
  TH1F *fhGenPi0AccPt  ; // pt of primary pi0
  TH1F *fhGenPi0AccEta ; // eta of primart pi0
  TH1F *fhGenPi0AccPhi ; // phi of primary pi0		
	
  //Histograms for track-matching
  TH1F *fh1pOverE;     //! p/E for track-cluster matches
  TH1F *fh1dR;         //! distance between projected track and cluster
  TH2F *fh2EledEdx;    //! dE/dx vs. momentum for electron candidates
  TH2F *fh2MatchdEdx;  //! dE/dx vs. momentum for all matches
	
	ClassDef(AliAnaCalorimeterQA,1)
} ;


#endif //ALIANACALORIMETERQA_H



