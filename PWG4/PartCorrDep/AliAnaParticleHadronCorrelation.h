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
// 4. Make decay photon-hadron correlations where decay contribute pi0 mass (2010/09/09)
// 5. fill the pout to extract kt at the end, also to study charge asymmetry(2010/10/06) 
// 6. Add the possibility for event selection analysis based on vertex and multiplicity bins (10/10/2010)
// 7. change the way of delta phi cut for UE study due to memory issue (reduce histograms)
// 8. Add the possibility to request the absolute leading particle at the near side or not, set trigger bins, general clean-up (08/2011)

// --- ROOT system ---
//class TH3D;

// --- Analysis system ---
#include "AliAnaPartCorrBaseClass.h"
class AliAODPWG4ParticleCorrelation ;

class AliAnaParticleHadronCorrelation : public AliAnaPartCorrBaseClass {
  
 public: 
  
  AliAnaParticleHadronCorrelation() ; // default ctor
  virtual ~AliAnaParticleHadronCorrelation() {;} //virtual dtor
  
  // General methods
  
  TObjString * GetAnalysisCuts();
  
  TList      * GetCreateOutputObjects();
  
  void         InitParameters();
  
  void         MakeAnalysisFillAOD()  ;
  
  void         MakeAnalysisFillHistograms() ; 
  
  void         Print(const Option_t * opt) const;
  
  // Main analysis methods
  
  Bool_t       MakeChargedCorrelation  (AliAODPWG4ParticleCorrelation * aodParticle, const TObjArray* pl, const Bool_t bFillHisto) ;
  
  Bool_t       MakeNeutralCorrelation  (AliAODPWG4ParticleCorrelation * aodParticle, const TObjArray* pl, const Bool_t bFillHisto) ;
  
  void         MakeMCChargedCorrelation(AliAODPWG4ParticleCorrelation *aodParticle);
  
  
  // Parameter setter and getter
  
  Double_t     GetDeltaPhiMaxCut()       const { return fDeltaPhiMaxCut   ; }
  Double_t     GetDeltaPhiMinCut()       const { return fDeltaPhiMinCut   ; }
  Double_t     GetUeDeltaPhiMaxCut()     const { return fUeDeltaPhiMaxCut ; }
  Double_t     GetUeDeltaPhiMinCut()     const { return fUeDeltaPhiMinCut ; }
  
  void         SetDeltaPhiCutRange(Double_t phimin, Double_t phimax)
    { fDeltaPhiMaxCut   = phimax ;   fDeltaPhiMinCut   = phimin   ; }
  void         SetUeDeltaPhiCutRange(Double_t uephimin, Double_t uephimax)
    { fUeDeltaPhiMaxCut = uephimax;  fUeDeltaPhiMinCut = uephimin ; }
  
  Bool_t       IsSeveralUEOn()           const { return fMakeSeveralUE    ; }
  void         SwitchOnSeveralUECalculation()  { fMakeSeveralUE = kTRUE   ; }
  void         SwitchOffSeveralUECalculation() { fMakeSeveralUE = kFALSE  ; }

  // Do trigger-neutral correlation
  Bool_t       DoNeutralCorr()           const { return fNeutralCorr      ; }
  void         SwitchOnNeutralCorr()           { fNeutralCorr = kTRUE     ; }
  void         SwitchOffNeutralCorr()          { fNeutralCorr = kFALSE    ; }  
  
  // Taking the absolute leading as the trigger or not
  Bool_t       DoAbsoluteLeading()       const { return fMakeAbsoluteLeading   ; }
  void         SwitchOnAbsoluteLeading()       { fMakeAbsoluteLeading = kTRUE  ; }
  void         SwitchOffAbsoluteLeading()      { fMakeAbsoluteLeading = kFALSE ; }
  
  // Do decay-hadron correlation if it is pi0 trigger
  Bool_t       IsPi0Trigger()            const { return fPi0Trigger       ; }
  void         SwitchOnDecayCorr()             { fPi0Trigger = kTRUE      ; }
  void         SwitchOffDecayCorr()            { fPi0Trigger = kFALSE     ; }  
  
  Bool_t       OnlyIsolated()            const { return fSelectIsolated   ; }
  void         SelectIsolated(Bool_t s)        { fSelectIsolated   = s    ; }

  void         SetPi0AODBranchName(TString n)  { fPi0AODBranchName = n    ; }
  
  void         SetNAssocPtBins(Int_t n) ;     
  void         SetAssocPtBinLimit(Int_t ibin, Float_t pt) ;
                
 private:
  
  Double_t     fDeltaPhiMaxCut ;               // Minimum Delta Phi Gamma-Hadron
  Double_t     fDeltaPhiMinCut ;               // Maximum Delta Phi Gamma-Hadron
  Bool_t       fSelectIsolated ;               // Select only trigger particles isolated
  Bool_t       fMakeSeveralUE ;                // Do analysis for several underlying events contribution
  Double_t     fUeDeltaPhiMaxCut ;             // Minimum Delta Phi Gamma-Underlying Hadron
  Double_t     fUeDeltaPhiMinCut ;             // Maximum Delta Phi Gamma-Underlying Hadron
  TString      fPi0AODBranchName;              // Name of AOD branch with pi0, not trigger
  Bool_t       fNeutralCorr ;                  // switch the analysis with neutral particles
  Bool_t       fPi0Trigger ;                   // switch the analysis with decay photon from pi0 trigger
  Bool_t       fMakeAbsoluteLeading ;          // requesting absolute leading while it is cluster triggers
  Int_t        fLeadingTriggerIndex ;          // Store here per event the trigger index, to avoid too many loops
  
  Int_t        fNAssocPtBins ;                 // Number of associated pT bins under study
  Float_t      fAssocPtBinLimit[10] ;          // Associated pT under study
  
  //Histograms

  //leading particles 
  TH1F *       fhPtLeading;                    //! pT distribution of leading particles
  TH2F *       fhPhiLeading;                   //! phi distribution vs pT of leading particles
  TH2F *       fhEtaLeading;                   //! eta distribution vs pT of leading particles
  
  //trigger-charged histograms
  TH2F *       fhDeltaPhiDeltaEtaCharged ;     //! differences of eta and phi between trigger and charged hadrons
  TH2F *       fhPhiCharged  ;                 //! Phi distribution of charged particles
  TH2F *       fhEtaCharged  ;                 //! Eta distribution of charged particles
  TH2F *       fhDeltaPhiCharged  ;            //! Difference of charged particle phi and trigger particle  phi as function of  trigger particle pT
  TH2F *       fhDeltaEtaCharged  ;            //! Difference of charged particle eta and trigger particle  eta as function of  trigger particle pT
  TH2F *       fhDeltaPhiChargedPt  ;          //! Difference of charged particle phi and trigger particle  phi as function of charged particle pT
  TH2F *       fhDeltaPhiUeChargedPt ;         //! Difference of charged particle from underlying events phi and trigger particle  phi as function of charged particle pT
  TH2F *       fhPtImbalanceCharged  ;         //! Trigger particle -charged hadron momentim imbalance histogram
  TH2F *       fhPtImbalanceUeCharged  ;       //! Trigger particle -underlying charged hadron momentum imbalance histogram  
  TH2F *       fhPtImbalancePosCharged  ;      //! Trigger particle -positive charged hadron momentum imbalance histogram
  TH2F *       fhPtImbalanceNegCharged  ;      //! Trigger particle -negative charged hadron momentum imbalance histogram 
  
  //with different imblance varible defination HBP distribution
  TH2F *       fhPtHbpCharged  ;               //! Trigger particle -charged hadron momentim HBP histogram
  TH2F *       fhPtHbpUeCharged  ;             //! Trigger particle -underlying charged hadron momentim HBP histogram  
  
  //if several UE calculation is on, most useful for jet-jet events contribution
  TH2F *       fhDeltaPhiUeLeftCharged  ;      //! Difference of charged particle from underlying events phi and trigger particle  phi as function of charged particle pT
  TH2F *       fhDeltaPhiUeRightCharged  ;     //! Difference of charged particle from underlying events phi and trigger particle  phi 
  TH2F *       fhPtImbalanceUeLeftCharged  ;   //! Trigger particle -underlying charged hadron momentim imbalance histogram 
  TH2F *       fhPtImbalanceUeRightCharged ;   //! Trigger particle -underlying charged hadron momentim imbalance histogram  
  TH2F *       fhPtHbpUeLeftCharged  ;         //! Trigger particle -underlying charged hadron momentim HBP histogram 
  TH2F *       fhPtHbpUeRightCharged  ;        //! Trigger particle -underlying charged hadron momentim HBP histogram  

  //for pout and kt extraction
  TH2F *       fhPtTrigPout  ;                 //! Pout =associated pt*sin(delta phi) distribution vs trigger pt 
  TH2F *       fhPtTrigCharged ;               //! trigger and correlated particl pt, to be used for mean value for kt	
  
  //if different multiplicity analysis asked
  TH2F **      fhTrigDeltaPhiCharged ;         //![GetMultiBin()] differences of phi between trigger and charged hadrons
  TH2F **      fhTrigDeltaEtaCharged ;         //![GetMultiBin()] differences of eta between trigger and charged hadrons
  TH2F **      fhTrigCorr  ;                   //![GetMultiBin()] Trigger particle -charged hadron momentim imbalance histogram
  TH2F **      fhTrigUeCorr  ;                 //![GetMultiBin()] Trigger particle -UE charged hadron momentim imbalance histogram
    
  TH2F *       fhAssocPt ;                     //! Trigger pT vs associated pT 
  TH2F *       fhAssocPtBkg;                   //! Trigger pT vs associated pT for background
  TH2F **      fhDeltaPhiAssocPtBin;           //![fNAssocPtBins] Trigger pT vs dPhi for different associated pt bins
  TH2F **      fhDeltaPhiBradAssocPtBin;       //![fNAssocPtBins] Trigger pT vs dPhi Brad (?) for different associated pt bins
  TH2F **      fhXEAssocPtBin ;                //![fNAssocPtBins] Trigger pT vs xE for different associated pt bins
  
  //trigger-neutral histograms
  TH2F *       fhDeltaPhiDeltaEtaNeutral ;     //! differences of eta and phi between trigger and neutral hadrons (pi0)
  TH2F *       fhPhiNeutral   ;                //! Phi distribution of neutral particles  
  TH2F *       fhEtaNeutral   ;                //! Eta distribution of neutral particles
  TH2F *       fhDeltaPhiNeutral   ;           //! Difference of neutral particle phi and trigger particle  phi as function of  trigger particle pT
  TH2F *       fhDeltaEtaNeutral  ;            //! Difference of neutral particle eta and trigger particle  eta as function of  trigger particle pT
  TH2F *       fhDeltaPhiNeutralPt  ;          //! Difference of neutral particle phi and trigger particle  phi as function of neutral particle particle pT
  TH2F *       fhDeltaPhiUeNeutralPt ;         //! Difference of neutral particle phi and trigger particle  phi as function of neutral particle particle pT  
  TH2F *       fhPtImbalanceNeutral  ;         //! Trigger particle - neutral hadron momentum imbalance histogram 
  TH2F *       fhPtImbalanceUeNeutral  ;       //! Trigger particle - neutral hadron momentum imbalance histogram 
  
  //with different imblance varible defination HBP distribution  
  TH2F *       fhPtHbpNeutral  ;               //! Trigger particle -neutral particle momentim HBP histogram
  TH2F *       fhPtHbpUeNeutral  ;             //! Trigger particle -underlying neutral hadron momentim HBP histogram  

  //if several UE calculation is on, most useful for jet-jet events contribution
  TH2F *       fhDeltaPhiUeLeftNeutral  ;      //! Difference of charged particle from underlying events phi and trigger particle  phi as function of neutral particle pT
  TH2F *       fhDeltaPhiUeRightNeutral  ;     //! Difference of charged particle from underlying events phi and trigger particle  phi 
  TH2F *       fhPtImbalanceUeLeftNeutral  ;   //! Trigger particle -underlying neutral hadron momentim imbalance histogram 
  TH2F *       fhPtImbalanceUeRightNeutral ;   //! Trigger particle -underlying neutral hadron momentim imbalance histogram 
  TH2F *       fhPtHbpUeLeftNeutral  ;         //! Trigger particle -underlying neutral hadron momentim HBP histogram 
  TH2F *       fhPtHbpUeRightNeutral  ;        //! Trigger particle -underlying neutral hadron momentim HBP histogram 
  
  //for decay photon trigger correlation
  TH2F *       fhPtPi0DecayRatio ;             //! for pi0 pt and ratio of decay photon pt
  TH2F *       fhDeltaPhiDecayCharged  ;       //! Difference of charged particle phi and decay trigger
  TH2F *       fhPtImbalanceDecayCharged ;     //! Trigger particle (decay from pi0)-charged hadron momentim imbalance histogram    
  TH2F *       fhDeltaPhiDecayNeutral  ;       //! Difference of neutral particle phi and decay trigger
  TH2F *       fhPtImbalanceDecayNeutral ;     //! Trigger particle (decay from pi0)-neutral hadron momentim imbalance histogram  
  
  //if the data is MC, fill MC information
  TH2F *       fh2phiLeadingParticle;          //! #phi resolution for triggers
  TH1F *       fhMCLeadingCount;               //! add explanation
  TH2F *       fhMCEtaCharged;                 //! add explanation
  TH2F *       fhMCPhiCharged;                 //! add explanation
  TH2F *       fhMCDeltaEtaCharged;            //! add explanation
  TH2F *       fhMCDeltaPhiCharged;            //! add explanation
  TH2F *       fhMCDeltaPhiDeltaEtaCharged;    //! add explanation
  TH2F *       fhMCDeltaPhiChargedPt;          //! add explanation
  TH2F *       fhMCPtImbalanceCharged;         //! add explanation
  TH2F *       fhMCPtHbpCharged;               //! add explanation
  TH2F *       fhMCPtTrigPout ;                //! add explanation
  TH2F *       fhMCPtAssocDeltaPhi  ;          //! Pout =associated pt*sin(delta phi) distribution

  AliAnaParticleHadronCorrelation(const AliAnaParticleHadronCorrelation & ph) ; // cpy ctor
  AliAnaParticleHadronCorrelation & operator = (const AliAnaParticleHadronCorrelation & ph) ;//cpy assignment
	
  ClassDef(AliAnaParticleHadronCorrelation,8)
} ;
 

#endif //ALIANAPARTICLEHADRONCORRELATION_H



