#ifndef ALIANAPARTICLEHADRONCORRELATION_H
#define ALIANAPARTICLEHADRONCORRELATION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

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

// --- Analysis system ---
#include "AliAnaCaloTrackCorrBaseClass.h"
class AliAODPWG4ParticleCorrelation ;

class AliAnaParticleHadronCorrelation : public AliAnaCaloTrackCorrBaseClass {
  
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
  
  Bool_t       GetDecayPhotonMomentum(const AliAODPWG4Particle* trigger, TLorentzVector & mom1,TLorentzVector & mom2);
  
  Bool_t       MakeChargedCorrelation  (AliAODPWG4ParticleCorrelation * aodParticle, const TObjArray* pl, const Bool_t bFillHisto) ;
  
  Bool_t       MakeNeutralCorrelation  (AliAODPWG4ParticleCorrelation * aodParticle, const TObjArray* pl, const Bool_t bFillHisto) ;
  
  void         MakeMCChargedCorrelation(AliAODPWG4ParticleCorrelation * aodParticle);
  
  // Filling histogram methods
  
  void         FillChargedAngularCorrelationHistograms  (const Float_t ptAssoc,  const Float_t ptTrig,   const Int_t   assocBin,
                                                         const Float_t phiAssoc, const Float_t phiTrig,  Float_t &     deltaPhi,
                                                         const Float_t etaAssoc, const Float_t etaTrig,  
                                                         const Bool_t  decay,    const Float_t hmpidSignal,const Int_t nTracks);
  
  Bool_t       FillChargedMCCorrelationHistograms       (const Float_t mcAssocPt,      Float_t mcAssocPhi, const Float_t mcAssocEta,
                                                         const Float_t mcTrigPt, const Float_t mcTrigPhi,  const Float_t mcTrigEta  );

  
  void         FillChargedMomentumImbalanceHistograms   (const Float_t ptTrig,   const Float_t ptAssoc, 
                                                         const Float_t xE,       const Float_t hbpXE, 
                                                         const Float_t zT,       const Float_t hbpZT, 
                                                         const Float_t pout,     const Float_t deltaPhi,
                                                         const Int_t   nTracks,  const Int_t   charge,
                                                         const Int_t   assocBin, const Bool_t  decay );
  
  void         FillChargedUnderlyingEventHistograms     (const Float_t ptTrig,   const Float_t ptAssoc, 
                                                         const Float_t deltaPhi, const Int_t nTracks);
  
  void         FillChargedUnderlyingEventSidesHistograms(const Float_t ptTrig,   const Float_t ptAssoc, 
                                                         const Float_t xE,       const Float_t hbpXE, 
                                                         const Float_t zT,       const Float_t hbpZT, 
                                                         const Float_t deltaPhi);
  
  void         FillDecayPhotonCorrelationHistograms     (const Float_t ptAssoc,     const Float_t phiAssoc, 
                                                         const TLorentzVector mom1, const TLorentzVector mom2, 
                                                         const Bool_t bChargedOrNeutral); 
  
  
  void         FillNeutralAngularCorrelationHistograms  (const Float_t ptAssoc,  const Float_t ptTrig,
                                                         const Float_t phiAssoc, const Float_t phiTrig,  Float_t &     deltaPhi,
                                                         const Float_t etaAssoc, const Float_t etaTrig);
  
  void         FillNeutralUnderlyingEventSidesHistograms(const Float_t ptTrig,   const Float_t ptAssoc, 
                                                         const Float_t xE,       const Float_t hbpXE, 
                                                         const Float_t zT,       const Float_t hbpZT, 
                                                         const Float_t deltaPhi);  
  
  // Parameter setter and getter
  
  Float_t      GetMinimumTriggerPt()       const { return fMinTriggerPt          ; }
  
  Float_t      GetMaximumAssociatedPt()    const { return fMaxAssocPt            ; }
  Float_t      GetMinimumAssociatedPt()    const { return fMinAssocPt            ; }
  
  Double_t     GetDeltaPhiMaxCut()         const { return fDeltaPhiMaxCut        ; }
  Double_t     GetDeltaPhiMinCut()         const { return fDeltaPhiMinCut        ; }
  
  Double_t     GetUeDeltaPhiMaxCut()       const { return fUeDeltaPhiMaxCut      ; }
  Double_t     GetUeDeltaPhiMinCut()       const { return fUeDeltaPhiMinCut      ; }
  
  void         SetMinimumTriggerPt(Float_t min)  { fMinTriggerPt = min           ; }
  
  void         SetAssociatedPtRange(Float_t min, Float_t max)
                  { fMaxAssocPt   = max ;           fMinAssocPt  = min           ; }
  
  void         SetDeltaPhiCutRange(Double_t phimin, Double_t phimax)
                  { fDeltaPhiMaxCut   = phimax ;    fDeltaPhiMinCut   = phimin   ; }
  
  void         SetUeDeltaPhiCutRange(Double_t uephimin, Double_t uephimax)
                  { fUeDeltaPhiMaxCut = uephimax ;  fUeDeltaPhiMinCut = uephimin ; }
  
  Bool_t       IsSeveralUEOn()             const { return fMakeSeveralUE         ; }
  void         SwitchOnSeveralUECalculation()    { fMakeSeveralUE      = kTRUE   ; }
  void         SwitchOffSeveralUECalculation()   { fMakeSeveralUE      = kFALSE  ; }

  // Do trigger-neutral correlation
  Bool_t       DoNeutralCorr()             const { return fNeutralCorr           ; }
  void         SwitchOnNeutralCorr()             { fNeutralCorr      = kTRUE     ; }
  void         SwitchOffNeutralCorr()            { fNeutralCorr      = kFALSE    ; }  
  
  // Taking the absolute leading as the trigger or not
  Bool_t       DoAbsoluteLeading()         const { return fMakeAbsoluteLeading   ; }
  void         SwitchOnAbsoluteLeading()         { fMakeAbsoluteLeading = kTRUE  ; }
  void         SwitchOffAbsoluteLeading()        { fMakeAbsoluteLeading = kFALSE ; }
  
  // Taking the near side leading as the trigger or not
  Bool_t       DoNearSideLeading()         const { return fMakeNearSideLeading   ; }
  void         SwitchOnNearSideLeading()         { fMakeNearSideLeading = kTRUE  ; }
  void         SwitchOffNearSideLeading()        { fMakeNearSideLeading = kFALSE ; }
  
  // Do decay-hadron correlation if it is pi0 trigger
  Bool_t       IsPi0Trigger()              const { return fPi0Trigger            ; }
  void         SwitchOnPi0TriggerDecayCorr()     { fPi0Trigger          = kTRUE  ; }
  void         SwitchOffPi0TriggerDecayCorr()    { fPi0Trigger          = kFALSE ; }  

  Bool_t       IsDecayTrigger()            const { return fDecayTrigger          ; }
  void         SwitchOnDecayTriggerDecayCorr()   { fDecayTrigger        = kTRUE  ; }
  void         SwitchOffDecayTriggerDecayCorr()  { fDecayTrigger        = kFALSE ; }  

  Bool_t       IsHMPIDCorrelation()        const { return fHMPIDCorrelation      ; }
  void         SwitchOnHMPIDCorrelation()        { fHMPIDCorrelation    = kTRUE  ; }
  void         SwitchOffHMPIDCorrelation()       { fHMPIDCorrelation    = kFALSE ; }  
  
  void         SwitchOnFillBradHistograms()      { fFillBradHisto       = kTRUE  ; }
  void         SwitchOffFillBradHistograms()     { fFillBradHisto       = kFALSE ; }  
    
  Bool_t       OnlyIsolated()              const { return fSelectIsolated        ; }
  void         SelectIsolated(Bool_t s)          { fSelectIsolated   = s         ; }

  void         SetPi0AODBranchName(TString n)    { fPi0AODBranchName = n         ; }
  
  void         SetNAssocPtBins(Int_t n) ;     
  void         SetAssocPtBinLimit(Int_t ibin, Float_t pt) ;
                
 private:
  Float_t      fMinTriggerPt ;                 // Minimum trigger hadron pt
  Float_t      fMaxAssocPt ;                   // Maximum associated hadron pt
  Float_t      fMinAssocPt ;                   // Minimum associated hadron pt
  Double_t     fDeltaPhiMaxCut ;               // Minimum Delta Phi Gamma-Hadron
  Double_t     fDeltaPhiMinCut ;               // Maximum Delta Phi Gamma-Hadron
  Bool_t       fSelectIsolated ;               // Select only trigger particles isolated
  Bool_t       fMakeSeveralUE ;                // Do analysis for several underlying events contribution
  Double_t     fUeDeltaPhiMaxCut ;             // Minimum Delta Phi Gamma-Underlying Hadron
  Double_t     fUeDeltaPhiMinCut ;             // Maximum Delta Phi Gamma-Underlying Hadron
  TString      fPi0AODBranchName;              // Name of AOD branch with pi0, not trigger
  Bool_t       fNeutralCorr ;                  // switch the analysis with neutral particles
  Bool_t       fPi0Trigger ;                   // switch the analysis with decay photon from pi0 trigger
  Bool_t       fDecayTrigger ;                 // switch the analysis with decay photon from photon trigger
  Bool_t       fMakeAbsoluteLeading ;          // requesting absolute leading triggers
  Bool_t       fMakeNearSideLeading ;          // requesting near side leading (+-90º from trigger particle) triggers
  Int_t        fLeadingTriggerIndex ;          // Store here per event the trigger index, to avoid too many loops
  Bool_t       fHMPIDCorrelation    ;          // Correlate with particles on HMPID or its acceptance
  Bool_t       fFillBradHisto ;                // DPhi histograms calculated differently
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
  TH2F *       fhXECharged  ;                  //! Trigger particle -charged hadron momentum imbalance histogram
  TH2F *       fhXEUeCharged  ;                //! Trigger particle -underlying charged hadron momentum imbalance histogram  
  TH2F *       fhXEPosCharged  ;               //! Trigger particle -positive charged hadron momentum imbalance histogram
  TH2F *       fhXENegCharged  ;               //! Trigger particle -negative charged hadron momentum imbalance histogram 
  TH2F *       fhPtHbpXECharged  ;             //! Trigger particle -charged hadron momentum HBP histogram
  TH2F *       fhPtHbpXEUeCharged  ;           //! Trigger particle -underlying charged hadron momentum HBP histogram  
  TH2F *       fhZTCharged  ;                  //! Trigger particle -charged hadron momentum imbalance histogram
  TH2F *       fhZTUeCharged  ;                //! Trigger particle -underlying charged hadron momentum imbalance histogram  
  TH2F *       fhZTPosCharged  ;               //! Trigger particle -positive charged hadron momentum imbalance histogram
  TH2F *       fhZTNegCharged  ;               //! Trigger particle -negative charged hadron momentum imbalance histogram 
  TH2F *       fhPtHbpZTCharged  ;             //! Trigger particle -charged hadron momentum HBP histogram
  TH2F *       fhPtHbpZTUeCharged  ;           //! Trigger particle -underlying charged hadron momentum HBP histogram  
 
  //if several UE calculation is on, most useful for jet-jet events contribution
  TH2F *       fhDeltaPhiUeLeftCharged  ;      //! Difference of charged particle from underlying events phi and trigger particle  phi as function of charged particle pT
  TH2F *       fhDeltaPhiUeRightCharged  ;     //! Difference of charged particle from underlying events phi and trigger particle  phi 
  TH2F *       fhXEUeLeftCharged  ;            //! Trigger particle -underlying charged hadron momentum imbalance histogram 
  TH2F *       fhXEUeRightCharged ;            //! Trigger particle -underlying charged hadron momentum imbalance histogram  
  TH2F *       fhPtHbpXEUeLeftCharged  ;       //! Trigger particle -underlying charged hadron momentum HBP histogram 
  TH2F *       fhPtHbpXEUeRightCharged  ;      //! Trigger particle -underlying charged hadron momentum HBP histogram  
  TH2F *       fhZTUeLeftCharged  ;            //! Trigger particle -underlying charged hadron momentum imbalance histogram 
  TH2F *       fhZTUeRightCharged ;            //! Trigger particle -underlying charged hadron momentum imbalance histogram  
  TH2F *       fhPtHbpZTUeLeftCharged  ;       //! Trigger particle -underlying charged hadron momentum HBP histogram 
  TH2F *       fhPtHbpZTUeRightCharged  ;      //! Trigger particle -underlying charged hadron momentum HBP histogram 
  
  //for pout and kt extraction
  TH2F *       fhPtTrigPout  ;                 //! Pout =associated pt*sin(delta phi) distribution vs trigger pt 
  TH2F *       fhPtTrigCharged ;               //! trigger and correlated particl pt, to be used for mean value for kt	
  
  //if different multiplicity analysis asked
  TH2F **      fhTrigDeltaPhiCharged ;         //![GetMultiBin()] differences of phi between trigger and charged hadrons
  TH2F **      fhTrigDeltaEtaCharged ;         //![GetMultiBin()] differences of eta between trigger and charged hadrons
  TH2F **      fhTrigXECorr  ;                 //![GetMultiBin()] Trigger particle -charged hadron momentum imbalance histogram
  TH2F **      fhTrigXEUeCorr  ;               //![GetMultiBin()] Trigger particle -UE charged hadron momentum imbalance histogram
  TH2F **      fhTrigZTCorr  ;                 //![GetMultiBin()] Trigger particle -charged hadron momentum imbalance histogram
  TH2F **      fhTrigZTUeCorr  ;               //![GetMultiBin()] Trigger particle -UE charged hadron momentum imbalance histogram
  
  TH2F *       fhAssocPtBkg;                   //! Trigger pT vs associated pT for background
  TH2F **      fhDeltaPhiAssocPtBin;           //![fNAssocPtBins] Trigger pT vs dPhi for different associated pt bins
  TH2F **      fhDeltaPhiAssocPtBinHMPID;      //![fNAssocPtBins] Trigger pT vs dPhi for different associated pt bins, track with HMPID  
  TH2F **      fhDeltaPhiAssocPtBinHMPIDAcc;   //![fNAssocPtBins] Trigger pT vs dPhi for different associated pt bins, track with HMPIDAcc  
  TH2F **      fhDeltaPhiBradAssocPtBin;       //![fNAssocPtBins] Trigger pT vs dPhi Brad (?) for different associated pt bins
  TH2F *       fhDeltaPhiBrad;                 //! Trigger pT vs dPhi Brad (?) for different associated pt bins
  TH2F **      fhXEAssocPtBin ;                //![fNAssocPtBins] Trigger pT vs xE for different associated pt bins
  TH2F **      fhZTAssocPtBin ;                //![fNAssocPtBins] Trigger pT vs zT for different associated pt bins

  //trigger-neutral histograms
  TH2F *       fhDeltaPhiDeltaEtaNeutral ;     //! differences of eta and phi between trigger and neutral hadrons (pi0)
  TH2F *       fhPhiNeutral   ;                //! Phi distribution of neutral particles  
  TH2F *       fhEtaNeutral   ;                //! Eta distribution of neutral particles
  TH2F *       fhDeltaPhiNeutral   ;           //! Difference of neutral particle phi and trigger particle  phi as function of  trigger particle pT
  TH2F *       fhDeltaEtaNeutral  ;            //! Difference of neutral particle eta and trigger particle  eta as function of  trigger particle pT
  TH2F *       fhDeltaPhiNeutralPt  ;          //! Difference of neutral particle phi and trigger particle  phi as function of neutral particle particle pT
  TH2F *       fhDeltaPhiUeNeutralPt ;         //! Difference of neutral particle phi and trigger particle  phi as function of neutral particle particle pT  
  TH2F *       fhXENeutral  ;                  //! Trigger particle - neutral hadron momentum imbalance histogram 
  TH2F *       fhXEUeNeutral  ;                //! Trigger particle - neutral hadron momentum imbalance histogram 
  TH2F *       fhPtHbpXENeutral  ;             //! Trigger particle - neutral particle momentum HBP histogram
  TH2F *       fhPtHbpXEUeNeutral  ;           //! Trigger particle - underlying neutral hadron momentum HBP histogram  
  TH2F *       fhZTNeutral  ;                  //! Trigger particle - neutral hadron momentum imbalance histogram 
  TH2F *       fhZTUeNeutral  ;                //! Trigger particle - neutral hadron momentum imbalance histogram 
  TH2F *       fhPtHbpZTNeutral  ;             //! Trigger particle - neutral particle momentum HBP histogram
  TH2F *       fhPtHbpZTUeNeutral  ;           //! Trigger particle - underlying neutral hadron momentum HBP histogram  
  
  //if several UE calculation is on, most useful for jet-jet events contribution
  TH2F *       fhDeltaPhiUeLeftNeutral  ;      //! Difference of charged particle from underlying events phi and trigger particle  phi as function of neutral particle pT
  TH2F *       fhDeltaPhiUeRightNeutral  ;     //! Difference of charged particle from underlying events phi and trigger particle  phi 
  TH2F *       fhXEUeLeftNeutral  ;            //! Trigger particle -underlying neutral hadron momentum imbalance histogram 
  TH2F *       fhXEUeRightNeutral ;            //! Trigger particle -underlying neutral hadron momentum imbalance histogram 
  TH2F *       fhPtHbpXEUeLeftNeutral  ;       //! Trigger particle -underlying neutral hadron momentum HBP histogram 
  TH2F *       fhPtHbpXEUeRightNeutral  ;      //! Trigger particle -underlying neutral hadron momentum HBP histogram 
  TH2F *       fhZTUeLeftNeutral  ;            //! Trigger particle -underlying neutral hadron momentum imbalance histogram 
  TH2F *       fhZTUeRightNeutral ;            //! Trigger particle -underlying neutral hadron momentum imbalance histogram 
  TH2F *       fhPtHbpZTUeLeftNeutral  ;       //! Trigger particle -underlying neutral hadron momentum HBP histogram 
  TH2F *       fhPtHbpZTUeRightNeutral  ;      //! Trigger particle -underlying neutral hadron momentum HBP histogram 
  
  //for decay photon trigger correlation
  TH2F *       fhPtPi0DecayRatio ;             //! for pi0 pt and ratio of decay photon pt
  TH2F *       fhDeltaPhiDecayCharged  ;       //! Difference of charged particle phi and decay trigger
  TH2F *       fhXEDecayCharged ;              //! Trigger particle (decay from pi0)-charged hadron momentum imbalance histogram    
  TH2F *       fhZTDecayCharged ;              //! Trigger particle (decay from pi0)-charged hadron momentum imbalance histogram   

  TH2F *       fhDeltaPhiDecayNeutral  ;       //! Difference of neutral particle phi and decay trigger
  TH2F *       fhXEDecayNeutral ;              //! Trigger particle (decay from pi0)-neutral hadron momentum imbalance histogram  
  TH2F *       fhZTDecayNeutral ;              //! Trigger particle (decay from pi0)-neutral hadron momentum imbalance histogram  

  TH2F **      fhDeltaPhiDecayChargedAssocPtBin;//![fNAssocPtBins] Tagged as decay Trigger pT vs dPhi for different associated pt bins
  TH2F **      fhXEDecayChargedAssocPtBin ;    //![fNAssocPtBins] Tagged as decay Trigger pT vs xE for different associated pt bins
  TH2F **      fhZTDecayChargedAssocPtBin ;    //![fNAssocPtBins] Tagged as decay Trigger pT vs xE for different associated pt bins  
  
  //if the data is MC, fill MC information
  TH2F *       fh2phiLeadingParticle;          //! #phi resolution for triggers
  TH2F *       fhMCEtaCharged;                 //! MC pure particles charged primary pt vs eta (both associated) 
  TH2F *       fhMCPhiCharged;                 //! MC pure particles charged primary pt vs phi (both associated) 
  TH2F *       fhMCDeltaEtaCharged;            //! MC pure particles charged trigger primary pt vs delta eta (associated-trigger) 
  TH2F *       fhMCDeltaPhiCharged;            //! MC pure particles charged trigger primary pt vs delta phi (associated-trigger) 
  TH2F *       fhMCDeltaPhiDeltaEtaCharged;    //! MC pure particles charged associated primary pt vs delta phi (associated-trigger), in away side 
  TH2F *       fhMCDeltaPhiChargedPt;          //! MC pure particles charged delta phi vs delta eta (associated-trigger) 
  TH2F *       fhMCPtXECharged;                //! MC pure particles charged trigger primary pt vs xE
  TH2F *       fhMCPtHbpXECharged;             //! MC pure particles charged trigger primary pt vs ln(1/xE)
  TH2F *       fhMCPtZTCharged;                //! MC pure particles charged trigger primary pt vs zT
  TH2F *       fhMCPtHbpZTCharged;             //! MC pure particles charged trigger primary pt vs ln(1/zT)
  TH2F *       fhMCPtTrigPout ;                //! MC pure particles charged trigger primary pt vs pOut
  TH2F *       fhMCPtAssocDeltaPhi  ;          //! MC pure particles charged associated primary pt vs delta phi (associated-trigger) 

  AliAnaParticleHadronCorrelation(              const AliAnaParticleHadronCorrelation & ph) ; // cpy ctor
  AliAnaParticleHadronCorrelation & operator = (const AliAnaParticleHadronCorrelation & ph) ; // cpy assignment
	
  ClassDef(AliAnaParticleHadronCorrelation,12)
} ;
 

#endif //ALIANAPARTICLEHADRONCORRELATION_H



