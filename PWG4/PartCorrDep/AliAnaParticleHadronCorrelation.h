#ifndef ALIANAPARTICLEHADRONCORRELATION_H
#define ALIANAPARTICLEHADRONCORRELATION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id:  $ */

//_________________________________________________________________________
// Class that contains the algorithm for the analysis of particle - hadron correlations
// Particle (for example direct gamma) must be found in a previous analysis 
//-- Author: Gustavo Conesa (INFN-LNF)

//  Modified by Yaxian Mao:
// 1. add the UE subtraction for corrlation study
// 2. change the correlation variable
// 3. Only use leading particle(cluster/track) as trigger for correlation (2010/07/02)

// --- ROOT system ---
class TH2F;

// --- Analysis system ---
#include "AliAnaPartCorrBaseClass.h"
class AliAODPWG4ParticleCorrelation ;

class AliAnaParticleHadronCorrelation : public AliAnaPartCorrBaseClass {
  
 public: 
  AliAnaParticleHadronCorrelation() ; // default ctor
  virtual ~AliAnaParticleHadronCorrelation() {;} //virtual dtor
 private:  
  AliAnaParticleHadronCorrelation(const AliAnaParticleHadronCorrelation & ph) ; // cpy ctor
  AliAnaParticleHadronCorrelation & operator = (const AliAnaParticleHadronCorrelation & ph) ;//cpy assignment

 public:

  TList * GetCreateOutputObjects();
  
  Double_t GetDeltaPhiMaxCut() const {return fDeltaPhiMaxCut ; }
  Double_t GetDeltaPhiMinCut() const {return fDeltaPhiMinCut ; }
  void SetDeltaPhiCutRange(Double_t phimin, Double_t phimax)
  {fDeltaPhiMaxCut =phimax;  fDeltaPhiMinCut =phimin;}

  Double_t GetUeDeltaPhiMaxCut() const {return fUeDeltaPhiMaxCut ; }
  Double_t GetUeDeltaPhiMinCut() const {return fUeDeltaPhiMinCut ; }
  void SetUeDeltaPhiCutRange(Double_t uephimin, Double_t uephimax)
  {fUeDeltaPhiMaxCut =uephimax;  fUeDeltaPhiMinCut =uephimin;}
  Bool_t IsSeveralUEOn() const {return fMakeSeveralUE ; }
  void SwitchOnSeveralUECalculation()  { fMakeSeveralUE = kTRUE;}
  void SwitchOffSeveralUECalculation() { fMakeSeveralUE = kFALSE;}

  Bool_t OnlyIsolated() const {return fSelectIsolated ; }
  void SelectIsolated(Bool_t select) {fSelectIsolated = select ; }
  
  void InitParameters();
  
  void Print(const Option_t * opt) const;
  
  void MakeChargedCorrelation(AliAODPWG4ParticleCorrelation * aodParticle,TObjArray* const pl, const Bool_t bFillHisto) ;
  void MakeNeutralCorrelationFillAOD(AliAODPWG4ParticleCorrelation* const aodParticle, TObjArray* const pl, TString detector)  ;
  void MakeNeutralCorrelationFillHistograms(AliAODPWG4ParticleCorrelation* const aodParticle)  ;
	
  void MakeAnalysisFillAOD()  ;
  void MakeAnalysisFillHistograms() ; 
  
  Bool_t SelectCluster(AliAODCaloCluster * calo, Double_t *vertex, TLorentzVector & mom, Int_t & pdg) const ;
  
 private:
  
  Double_t   fDeltaPhiMaxCut ;      // Minimum Delta Phi Gamma-Hadron
  Double_t   fDeltaPhiMinCut ;      // Maximum Delta Phi Gamma-Hadron
  Bool_t     fSelectIsolated ;      // Select only trigger particles isolated
  Bool_t   fMakeSeveralUE ; // Do analysis for several underlying events contribution
  Double_t   fUeDeltaPhiMaxCut ;      // Minimum Delta Phi Gamma-Underlying Hadron
  Double_t   fUeDeltaPhiMinCut ;      // Maximum Delta Phi Gamma-Underlying Hadron

  
  //Histograms
  //leading particles 
  TH1F * fhPtLeading;         //! pT distribution of leading particles
  TH2F * fhPhiLeading;        //! phi distribution vs pT of leading particles
  TH2F * fhEtaLeading;        //! eta distribution vs pT of leading particles
  TH2F * fhDeltaPhiDeltaEtaCharged ; //! differences of eta and phi between trigger and charged hadrons
  TH2F * fhDeltaPhiDeltaEtaNeutral ; //! differences of eta and phi between trigger and neutral hadrons (pi0)
	
  TH2F * fhPhiCharged  ; //! Phi distribution of charged particles
  TH2F * fhPhiNeutral   ;  //! Phi distribution of neutral particles
  TH2F * fhEtaCharged  ; //! Eta distribution of charged particles
  TH2F * fhEtaNeutral   ; //! Eta distribution of neutral particles
  TH2F * fhDeltaPhiCharged  ;  //! Difference of charged particle phi and trigger particle  phi as function of  trigger particle pT
  TH2F * fhDeltaPhiNeutral   ;  //! Difference of neutral particle phi and trigger particle  phi as function of  trigger particle pT
  TH2F * fhDeltaEtaCharged  ;  //! Difference of charged particle eta and trigger particle  eta as function of  trigger particle pT
  TH2F * fhDeltaEtaNeutral  ;  //! Difference of neutral particle eta and trigger particle  eta as function of  trigger particle pT
  TH2F * fhDeltaPhiChargedPt  ;  //! Difference of charged particle phi and trigger particle  phi as function of charged particle pT
  TH2F * fhDeltaPhiNeutralPt  ;  //! Difference of neutral particle phi and trigger particle  phi as function of neutral particle particle pT
  TH2F * fhDeltaPhiUeChargedPt  ;  //! Difference of charged particle from underlying events phi and trigger particle  phi as function of charged particle pT
  TH2F * fhDeltaPhiUeNeutralPt  ;  //! Difference of neutral particle phi and trigger particle  phi as function of neutral particle particle pT

  TH2F * fhPtImbalanceNeutral  ; //! Trigger particle - neutral hadron momentum imbalance histogram 
  TH2F * fhPtImbalanceCharged  ; //! Trigger particle -charged hadron momentim imbalance histogram
  TH2F * fhPtImbalanceUeCharged  ; //! Trigger particle -underlying charged hadron momentim imbalance histogram  
  TH2F * fhPtImbalanceUeNeutral  ; //! Trigger particle - neutral hadron momentum imbalance histogram 

  //with different imblance varible defination HBP distribution
  TH2F * fhPtHbpCharged  ; //! Trigger particle -charged hadron momentim HBP histogram
  TH2F * fhPtHbpUeCharged  ; //! Trigger particle -underlying charged hadron momentim HBP histogram  
  TH2F * fhPtHbpNeutral  ; //! Trigger particle -neutral particle momentim HBP histogram
  TH2F * fhPtHbpUeNeutral  ; //! Trigger particle -underlying neutral hadron momentim HBP histogram  

  //if several UE calculation is on, most useful for jet-jet events contribution
  TH2F * fhDeltaPhiUeLeftCharged  ;  //! Difference of charged particle from underlying events phi and trigger particle  phi as function of charged particle pT
  TH2F * fhDeltaPhiUeRightCharged  ;  //! Difference of charged particle from underlying events phi and trigger particle  phi 
  TH2F * fhDeltaPhiUeLeftNeutral  ;  //! Difference of charged particle from underlying events phi and trigger particle  phi as function of neutral particle pT
  TH2F * fhDeltaPhiUeRightNeutral  ;  //! Difference of charged particle from underlying events phi and trigger particle  phi 
  TH2F * fhPtImbalanceUeLeftCharged  ; //! Trigger particle -underlying charged hadron momentim imbalance histogram 
  TH2F * fhPtImbalanceUeRightCharged  ; //! Trigger particle -underlying charged hadron momentim imbalance histogram  
  TH2F * fhPtImbalanceUeLeftNeutral  ; //! Trigger particle -underlying neutral hadron momentim imbalance histogram 
  TH2F * fhPtImbalanceUeRightNeutral  ; //! Trigger particle -underlying neutral hadron momentim imbalance histogram 
  TH2F * fhPtHbpUeLeftCharged  ; //! Trigger particle -underlying charged hadron momentim HBP histogram 
  TH2F * fhPtHbpUeRightCharged  ; //! Trigger particle -underlying charged hadron momentim HBP histogram  
  TH2F * fhPtHbpUeLeftNeutral  ; //! Trigger particle -underlying neutral hadron momentim HBP histogram 
  TH2F * fhPtHbpUeRightNeutral  ; //! Trigger particle -underlying neutral hadron momentim HBP histogram  
	
	
  ClassDef(AliAnaParticleHadronCorrelation,3)
} ;
 

#endif //ALIANAPARTICLEHADRONCORRELATION_H



