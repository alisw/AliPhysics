#ifndef ALIANAPARTICLEHADRONCORRELATION_H
#define ALIANAPARTICLEHADRONCORRELATION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
// Class that contains the algorithm for the analysis of
// particle - hadron correlations
// Particle (for example direct gamma) must be found in a previous analysis
//
//-- Author: Gustavo Conesa (LNF-INFN) (LPSC-IN2P3-CNRS)
//           Yaxian Mao (LPSC-IN2P3-CNRS) and (CNWU) first usable implementation.
//           Xiangrong Zhu (CNWU), implementtion of own mixing.
//

// --- Analysis system ---

#include "AliAnaCaloTrackCorrBaseClass.h"
class AliAODPWG4ParticleCorrelation ;

class AliAnaParticleHadronCorrelation : public AliAnaCaloTrackCorrBaseClass {
  
 public: 
  
  AliAnaParticleHadronCorrelation() ;          // default ctor
  virtual ~AliAnaParticleHadronCorrelation() ; // virtual dtor
  
  // General methods
      
  TObjString * GetAnalysisCuts();
  
  TList      * GetCreateOutputObjects();
  
  void         Init();
  
  void         InitParameters();

  void         FillEventMixPool() ;
    
  void         MakeAnalysisFillHistograms() ; 
  
  void         Print(const Option_t * opt) const;
  
  // Main analysis methods
  
  Bool_t       FindLeadingOppositeHadronInWindow(AliAODPWG4ParticleCorrelation * particle);
  
  Bool_t       GetDecayPhotonMomentum   (AliAODPWG4Particle* trigger, TLorentzVector & mom1, TLorentzVector & mom2);
  
  void         MakeChargedCorrelation   (AliAODPWG4ParticleCorrelation * particle) ;
  
  void         MakeNeutralCorrelation   (AliAODPWG4ParticleCorrelation * particle) ;
  
  void         MakeMCChargedCorrelation (Int_t triggerMCLable, Int_t histoIndex) ;
  
  void         MakeChargedMixCorrelation(AliAODPWG4ParticleCorrelation * particle) ;
  
  // Filling histogram methods
  
  void         FillChargedAngularCorrelationHistograms  (Float_t ptAssoc,  Float_t ptTrig,      Int_t   assocBin,
                                                         Float_t phiAssoc, Float_t phiTrig,     Float_t deltaPhi,
                                                         Float_t etaAssoc, Float_t etaTrig,  
                                                         Bool_t  decay,    Float_t hmpidSignal, Int_t outTOF,
                                                         Int_t   cenbin,   Int_t   mcTag);
  
  void         FillChargedEventMixPool();
  
  Bool_t       FillChargedMCCorrelationHistograms       (Float_t mcAssocPt, Float_t mcAssocPhi, Float_t mcAssocEta,
                                                         Float_t mcTrigPt,  Float_t mcTrigPhi,  Float_t mcTrigEta, Int_t histoIndex);

  
  void         FillChargedMomentumImbalanceHistograms   (Float_t ptTrig,   Float_t ptAssoc, 
                                                         Float_t deltaPhi, Int_t   cenbin, Int_t charge,
                                                         Int_t   assocBin, Bool_t  decay,
                                                         Int_t outTOF,     Int_t mcTag );
  
  void         FillChargedUnderlyingEventHistograms     (Float_t ptTrig,   Float_t ptAssoc, 
                                                         Float_t deltaPhi, Int_t cenbin, Int_t outTOF);
  
  void         FillChargedUnderlyingEventSidesHistograms(Float_t ptTrig,   Float_t ptAssoc, 
                                                         Float_t deltaPhi);
  
  void         FillDecayPhotonCorrelationHistograms     (Float_t ptAssoc,     Float_t phiAssoc, 
                                                         TLorentzVector mom1, TLorentzVector mom2, 
                                                         Bool_t bChargedOrNeutral); 
  
  void         FillNeutralEventMixPool();
  
  
  void         FillNeutralUnderlyingEventSidesHistograms(Float_t ptTrig,   Float_t ptAssoc, 
                                                         Float_t zT,       Float_t hbpZT, 
                                                         Float_t deltaPhi);  
    
  Int_t        GetMCTagHistogramIndex(Int_t tag);
  
  Bool_t       IsTriggerTheEventLeadingParticle();
  
  // Parameter setter and getter
  
  Float_t      GetMinimumTriggerPt()       const { return GetMinPt()             ; }
  Float_t      GetMaximumTriggerPt()       const { return GetMaxPt()             ; }
  void         SetTriggerPtRange(Float_t min, Float_t max)
               { SetMinPt(min), SetMaxPt(max)                                    ; }
  

  Float_t      GetMaximumAssociatedPt()    const { return fMaxAssocPt            ; }
  Float_t      GetMinimumAssociatedPt()    const { return fMinAssocPt            ; }
  void         SetAssociatedPtRange(Float_t min, Float_t max)
               { fMaxAssocPt   = max ;           fMinAssocPt  = min              ; }

  Double_t     GetDeltaPhiMaxCut()         const { return fDeltaPhiMaxCut        ; }
  Double_t     GetDeltaPhiMinCut()         const { return fDeltaPhiMinCut        ; }
  void         SetDeltaPhiCutRange(Double_t phimin, Double_t phimax)
               { fDeltaPhiMaxCut   = phimax ;    fDeltaPhiMinCut   = phimin      ; }
  
  // Leading Hadron
  Double_t     GetLeadHadronPhiMaxCut()    const { return fMaxLeadHadPhi         ; }
  Double_t     GetLeadHadronPhiMinCut()    const { return fMinLeadHadPhi         ; }
  void         SetLeadHadronPhiCut(Float_t min, Float_t max)
               { fMaxLeadHadPhi = max ;          fMinLeadHadPhi  = min           ; }

  Double_t     GetLeadHadronPtMinCut()     const { return fMinLeadHadPt          ; }
  Double_t     GetLeadHadronPtMaxCut()     const { return fMaxLeadHadPt          ; }
  void         SetLeadHadronPtCut(Float_t min, Float_t max)
               { fMaxLeadHadPt  = max ;           fMinLeadHadPt  = min           ; }
  
  Bool_t       IsLeadHadronCutOn()        const { return fSelectLeadingHadronAngle ; }
  void         SwitchOnLeadHadronSelection()    { fSelectLeadingHadronAngle = kTRUE  ; }
  void         SwitchOffLeadHadronSelection()   { fSelectLeadingHadronAngle = kFALSE ; }
  
  // UE
  
  Double_t     GetUeDeltaPhiMaxCut()       const { return fUeDeltaPhiMaxCut      ; }
  Double_t     GetUeDeltaPhiMinCut()       const { return fUeDeltaPhiMinCut      ; }
  
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
  
  Bool_t       IsMixStoredInReaderOn()     const { return fUseMixStoredInReader  ; }
  void         SwitchOnUseMixStoredInReader()    { fUseMixStoredInReader = kTRUE ; }
  void         SwitchOffUseMixStoredInReader()   { fUseMixStoredInReader = kFALSE; }
  
  void         SwitchOnFillNeutralInMixedEvent() { fFillNeutralEventMixPool = kTRUE  ; }
  void         SwitchOffFillNeutralInMixedEvent(){ fFillNeutralEventMixPool = kFALSE ; }
  
  void         SetM02Cut(Float_t min=0, Float_t max=10)  { fM02MinCut   = min ; fM02MaxCut  = max ; }

  void         SwitchOnCorrelationVzBin()        { fCorrelVzBin          = kTRUE  ; }
  void         SwitchOffCorrelationVzBin()       { fCorrelVzBin          = kFALSE ; }  
  
  void         SwitchOnFillPileUpHistograms()    { fFillPileUpHistograms = kTRUE  ; }
  void         SwitchOffFillPileUpHistograms()   { fFillPileUpHistograms = kFALSE ; }

  void         SwitchOnFillHighMultiplicityHistograms() { fFillHighMultHistograms = kTRUE  ; }
  void         SwitchOffFillHighMultiplicityHistograms(){ fFillHighMultHistograms = kFALSE ; }
  
  void         SwitchOnFillTriggerAODWithReferences()   { fFillAODWithReferences = kTRUE  ; }
  void         SwitchOffFillTriggerAODWithReferences()  { fFillAODWithReferences = kFALSE ; }

  void         SwitchOnCheckNeutralClustersForLeading() { fCheckLeadingWithNeutralClusters = kTRUE  ; }
  void         SwitchOffCheckNeutralClustersForLeading(){ fCheckLeadingWithNeutralClusters = kFALSE ; }
  
  void         SwitchOnFillEtaGapHistograms()    { fFillEtaGapsHisto    = kTRUE  ; }
  void         SwitchOffFillEtaGapHistograms()   { fFillEtaGapsHisto    = kFALSE ; }
  
  void         SwitchOnFillPtImbalancePerPtABinHistograms()  { fFillMomImbalancePtAssocBinsHisto = kTRUE  ; }
  void         SwitchOffFillPtImbalancePerPtABinHistograms() { fFillMomImbalancePtAssocBinsHisto = kFALSE ; }
  
  void         SetMCGenType(Int_t min = 0, Int_t max = 6) { if(min >= 0 && min < 7) fMCGenTypeMin = min ;
                                                            if(max >= 0 && max < 7) fMCGenTypeMax = max ; }
  
 private:

  Bool_t       fFillAODWithReferences;         // Add to the trigger particle AOD the reference to the tracks or neutrals in correlation.
  Bool_t       fCheckLeadingWithNeutralClusters;// Compare the trigger candidate to Leading pT with the clusters pT, by default only charged
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
  Float_t      fAssocPtBinLimit[20] ;          // Associated pT under study
  Bool_t       fCorrelVzBin ;                  // Fill one histogram per vz bin
  
  TList **     fListMixTrackEvents ;           //![GetNCentrBin()*GetNZvertBin()*GetNRPBin()] Containers for tracks in stored events for mixing
  TList **     fListMixCaloEvents ;            //![GetNCentrBin()*GetNZvertBin()*GetNRPBin()] Containers for calo clusters in stored events for mixing

  Bool_t       fUseMixStoredInReader;          // Signal if in the current event the pool was filled
  Bool_t       fFillNeutralEventMixPool;       // Add clusters to pool if requested
  
  Float_t      fM02MaxCut   ;                  // Study photon clusters with l0 smaller than cut
  Float_t      fM02MinCut   ;                  // Study photon clusters with l0 larger than cut
  
  Bool_t       fFillPileUpHistograms;          // Fill pile-up related histograms
  Bool_t       fFillHighMultHistograms;        // Histograms with centrality and event plane for triggers pT
  
  Bool_t       fSelectLeadingHadronAngle;      // Select events with leading particle within a range
  Float_t      fMinLeadHadPhi;                 // Minimum angle between the trigger and leading hadron
  Float_t      fMaxLeadHadPhi;                 // Maximum ange between the trigger and leading hadron
  Float_t      fMinLeadHadPt;                  // Minimum pT of leading hadron
  Float_t      fMaxLeadHadPt;                  // Maximum pT of leading hadron

  Bool_t       fFillEtaGapsHisto;              // Fill azimuthal correlation histograms in 2 eta gaps, |eta|>0.8 and |eta|<0.01
  Bool_t       fFillMomImbalancePtAssocBinsHisto; // momentum imbalance histograms in bins of pT associated
  
  Int_t        fMCGenTypeMin;                  // Of the 7 possible types, select those between fMCGenTypeMin and fMCGenTypeMax
  Int_t        fMCGenTypeMax;                  // Of the 7 possible types, select those between fMCGenTypeMin and fMCGenTypeMax
  
  //Histograms

  //trigger particles
  TH1F *       fhPtTriggerInput;               //! pT distribution of trigger particles before selection
  TH1F *       fhPtTriggerSSCut;               //! pT distribution of trigger particles after shower shape selection
  TH1F *       fhPtTriggerIsoCut;              //! pT distribution of trigger particles after isolation cut selection
  TH1F *       fhPtTriggerFidCut;              //! pT distribution of trigger particles after fiducial selection
  TH1F *       fhPtTrigger;                    //! pT distribution of trigger particles
  TH1F *       fhPtTriggerVtxBC0;              //! pT distribution of trigger particles
  TH1F *       fhPtTriggerPileUp[7];           //! pT distribution of trigger particles
  TH2F *       fhPtTriggerVzBin;               //! pT distribution of trigger particles vs vz bin
  TH2F *       fhPtTriggerBin;                 //! pT distribution of trigger particles, vs mixing bin
  TH2F *       fhPhiTrigger;                   //! phi distribution vs pT of trigger particles
  TH2F *       fhEtaTrigger;                   //! eta distribution vs pT of trigger particles
  
  TH1F *       fhPtTriggerMC[7];               //! pT distribution of trigger particles, check the origin of the cluster : "Photon","Pi0","Pi0Decay","EtaDecay","OtherDecay","Electron","Hadron"

  TH2F *       fhPtTriggerCentrality;          //! pT distribution of trigger particles vs centrality
  TH2F *       fhPtTriggerEventPlane;          //! pT distribution of trigger particles vs centrality
  TH2F *       fhTriggerEventPlaneCentrality;  //! event plane vs centrality for trigger particles
  
  TH1F *       fhPtTriggerMixed;               //! pT distribution of trigger particles, used in mixing
  TH2F *       fhPtTriggerMixedVzBin;          //! pT distribution of trigger particles, used in mixing, vs vz bin
  TH2F *       fhPtTriggerMixedBin;            //! pT distribution of trigger particles vs mixing bin
  TH2F *       fhPhiTriggerMixed;              //! phi distribution vs pT of trigger particles, used in mixing
  TH2F *       fhEtaTriggerMixed;              //! eta distribution vs pT of trigger particles, used in mixing  

  // Leading hadron in the opposite side of the trigger
  TH2F *       fhPtLeadingOppositeHadron;        //! pT trigger : pT distribution of leading hadron oposite to trigger
  TH2F *       fhPtDiffPhiLeadingOppositeHadron; //! pT trigger : difference phi distribution of leading hadron oposite and trigger
  TH2F *       fhPtDiffEtaLeadingOppositeHadron; //! pT trigger: difference eta distribution of leading hadron oposite and trigger

  //trigger-charged histograms
  TH2F *       fhDeltaPhiDeltaEtaCharged ;     //! differences of eta and phi between trigger and charged hadrons
  TH2F *       fhPhiCharged  ;                 //! Phi distribution of charged particles
  TH2F *       fhEtaCharged  ;                 //! Eta distribution of charged particles
  TH2F *       fhDeltaPhiCharged  ;            //! Difference of charged particle phi and trigger particle  phi as function of  trigger particle pT
  TH2F *       fhDeltaEtaCharged  ;            //! Difference of charged particle eta and trigger particle  eta as function of  trigger particle pT
  TH2F *       fhDeltaPhiChargedPt  ;          //! Difference of charged particle phi and trigger particle  phi as function of charged particle pT
  TH2F *       fhDeltaPhiUeChargedPt ;         //! Difference of charged particle from underlying events phi and trigger particle  phi as function of charged particle pT
  TH1F *       fhUePart;                       //! UE particles distribution vs pt trig
  TH2F *       fhXECharged  ;                  //! Trigger particle -charged hadron momentum imbalance histogram
  TH2F *       fhXECharged_Cone2  ;            //! Trigger particle -charged hadron momentum imbalance histogram in cone2 (5pi/6-7pi/6)
  TH2F *       fhXEUeCharged  ;                //! Trigger particle -underlying charged hadron momentum imbalance histogram  
  TH2F *       fhXEPosCharged  ;               //! Trigger particle -positive charged hadron momentum imbalance histogram
  TH2F *       fhXENegCharged  ;               //! Trigger particle -negative charged hadron momentum imbalance histogram 
  TH2F *       fhPtHbpXECharged  ;             //! Trigger particle -charged hadron momentum HBP histogram
  TH2F *       fhPtHbpXECharged_Cone2  ;       //! Trigger particle -charged hadron momentum HBP histogram in cone2 (5pi/6-7pi/6)
  TH2F *       fhPtHbpXEUeCharged  ;           //! Trigger particle -underlying charged hadron momentum HBP histogram  
  TH2F *       fhZTCharged  ;                  //! Trigger particle -charged hadron momentum imbalance histogram
  TH2F *       fhZTUeCharged  ;                //! Trigger particle -underlying charged hadron momentum imbalance histogram  
  TH2F *       fhZTPosCharged  ;               //! Trigger particle -positive charged hadron momentum imbalance histogram
  TH2F *       fhZTNegCharged  ;               //! Trigger particle -negative charged hadron momentum imbalance histogram 
  TH2F *       fhPtHbpZTCharged  ;             //! Trigger particle -charged hadron momentum HBP histogram
  TH2F *       fhPtHbpZTUeCharged  ;           //! Trigger particle -underlying charged hadron momentum HBP histogram  
  
  TH2F *       fhXEChargedMC[7]  ;             //! Trigger particle -charged hadron momentum imbalance histogram, check the origin of the cluster : decay photon (pi0, eta, other), merged photon (pi0), hadron, rest of photons (prompt, FSR, ISR)
  TH2F *       fhDeltaPhiChargedMC[7]  ;       //! Trigger particle -charged hadron delta phi histogram, check the origin of the cluster : decay photon (pi0, eta, other), merged photon (pi0), hadron, rest of photons (prompt, FSR, ISR)

  TH2F *       fhDeltaPhiDeltaEtaChargedPtA3GeV;//! differences of eta and phi between trigger and charged hadrons, pTa > 3 GeV
  TH2F *       fhDeltaPhiChargedPtA3GeV  ;      //! Difference of charged particle phi and trigger particle  phi as function of  trigger particle pT, pTa > 3 GeV
  TH2F *       fhDeltaEtaChargedPtA3GeV  ;      //! Difference of charged particle eta and trigger particle  eta as function of  trigger particle pT, pTa > 3 GeV
  
  // Events tagged as pileup by SDD,EMCal, or combination
  TH2F *       fhDeltaPhiChargedPileUp[7]  ;    //! Difference of charged particle phi and trigger particle  phi as function of  trigger particle pT
  TH2F *       fhDeltaEtaChargedPileUp[7]  ;    //! Difference of charged particle eta and trigger particle  eta as function of  trigger particle pT
  TH2F *       fhDeltaPhiChargedPtA3GeVPileUp[7] ; //! Difference of charged particle phi and trigger particle  phi as function of  trigger particle pT, pTa > 3 GeV
  TH2F *       fhDeltaEtaChargedPtA3GeVPileUp[7] ; //! Difference of charged particle eta and trigger particle  eta as function of  trigger particle pT, pTa > 3 GeV
  TH2F *       fhXEChargedPileUp[7]  ;          //! Trigger particle -charged hadron momentum imbalance histogram
  TH2F *       fhXEUeChargedPileUp[7]  ;        //! Trigger particle -charged hadron momentum imbalance histogram
  TH2F *       fhZTChargedPileUp[7]  ;          //! Trigger particle -charged hadron momentum imbalance histogram
  TH2F *       fhZTUeChargedPileUp[7]  ;        //! Trigger particle -charged hadron momentum imbalance histogram
  TH2F *       fhPtTrigChargedPileUp[7] ;       //! trigger and correlated particl pt, to be used for mean value for kt
  
  TH2F *       fhDeltaPhiChargedOtherBC  ;       //! Difference of charged particle phi and trigger particle  phi as function of  trigger particle pT
  TH2F *       fhDeltaPhiChargedPtA3GeVOtherBC ; //! Difference of charged particle phi and trigger particle  phi as function of  trigger particle pT, pTa > 3 GeV
  TH2F *       fhXEChargedOtherBC  ;             //! Trigger particle -charged hadron momentum imbalance histogram
  TH2F *       fhXEUeChargedOtherBC  ;           //! Trigger particle -charged hadron momentum imbalance histogram
  TH2F *       fhZTChargedOtherBC  ;             //! Trigger particle -charged hadron momentum imbalance histogram
  TH2F *       fhZTUeChargedOtherBC  ;           //! Trigger particle -charged hadron momentum imbalance histogram
  TH2F *       fhPtTrigChargedOtherBC ;          //! trigger and correlated particl pt, to be used for mean value for kt

  TH2F *       fhDeltaPhiChargedBC0  ;           //! Difference of charged particle phi and trigger particle  phi as function of  trigger particle pT
  TH2F *       fhDeltaPhiChargedPtA3GeVBC0 ;     //! Difference of charged particle phi and trigger particle  phi as function of  trigger particle pT, pTa > 3 GeV
  TH2F *       fhXEChargedBC0  ;                 //! Trigger particle -charged hadron momentum imbalance histogram
  TH2F *       fhXEUeChargedBC0  ;               //! Trigger particle -charged hadron momentum imbalance histogram
  TH2F *       fhZTChargedBC0  ;                 //! Trigger particle -charged hadron momentum imbalance histogram
  TH2F *       fhZTUeChargedBC0  ;               //! Trigger particle -charged hadron momentum imbalance histogram
  TH2F *       fhPtTrigChargedBC0 ;              //! trigger and correlated particl pt, to be used for mean value for kt

  TH2F *       fhDeltaPhiChargedVtxBC0  ;        //! Difference of charged particle phi and trigger particle  phi as function of  trigger particle pT
  TH2F *       fhDeltaPhiChargedPtA3GeVVtxBC0 ;  //! Difference of charged particle phi and trigger particle  phi as function of  trigger particle pT, pTa > 3 GeV
  TH2F *       fhXEChargedVtxBC0  ;              //! Trigger particle -charged hadron momentum imbalance histogram
  TH2F *       fhXEUeChargedVtxBC0  ;            //! Trigger particle -charged hadron momentum imbalance histogram
  TH2F *       fhZTChargedVtxBC0  ;              //! Trigger particle -charged hadron momentum imbalance histogram
  TH2F *       fhZTUeChargedVtxBC0  ;            //! Trigger particle -charged hadron momentum imbalance histogram
  TH2F *       fhPtTrigChargedVtxBC0 ;           //! trigger and correlated particl pt, to be used for mean value for kt
  
  //if several UE calculation is on, most useful for jet-jet events contribution
  TH2F *       fhDeltaPhiUeLeftCharged  ;      //! Difference of charged particle from underlying events phi and trigger particle  phi as function of charged particle pT
  TH2F *       fhDeltaPhiUeLeftUpCharged;      //! Difference of charged particle from underlying events phi and trigger particle  phi
  TH2F *       fhDeltaPhiUeRightUpCharged;     //! Difference of charged particle from underlying events phi and trigger particle  phi 
  TH2F *       fhDeltaPhiUeLeftDownCharged;    //! Difference of charged particle from underlying events phi and trigger particle  phi 
  TH2F *       fhDeltaPhiUeRightDownCharged;   //! Difference of charged particle from underlying events phi and trigger particle  phi 
  TH2F *       fhXEUeLeftCharged  ;            //! Trigger particle -underlying charged hadron momentum imbalance histogram 
  TH2F *       fhXEUeLeftUpCharged  ;          //! Trigger particle -underlying charged hadron momentum imbalance histogram
  TH2F *       fhXEUeRightUpCharged ;          //! Trigger particle -underlying charged hadron momentum imbalance histogram  
  TH2F *       fhXEUeLeftDownCharged  ;        //! Trigger particle -underlying charged hadron momentum imbalance histogram 
  TH2F *       fhXEUeRightDownCharged ;        //! Trigger particle -underlying charged hadron momentum imbalance histogram  
  TH2F *       fhPtHbpXEUeLeftCharged  ;       //! Trigger particle -underlying charged hadron momentum HBP histogram 
  TH2F *       fhZTUeLeftCharged  ;            //! Trigger particle -underlying charged hadron momentum imbalance histogram
  TH2F *       fhPtHbpZTUeLeftCharged  ;       //! Trigger particle -underlying charged hadron momentum HBP histogram
  
  //for pout and kt extraction
  TH2F *       fhPtTrigPout  ;                 //! Pout =associated pt*sin(delta phi) distribution vs trigger pt 
  TH2F *       fhPtTrigCharged ;               //! trigger and correlated particl pt, to be used for mean value for kt	
  
  //if different multiplicity analysis asked
  TH2F **      fhDeltaPhiChargedMult ;         //![GetNCentrBin()] differences of phi between trigger and charged hadrons: multiplicity bin
  TH2F **      fhDeltaEtaChargedMult ;         //![GetNCentrBin()] differences of eta between trigger and charged hadrons: multiplicity bin
  TH2F **      fhXEMult  ;                     //![GetNCentrBin()] Trigger particle -charged hadron momentum imbalance histogram: multiplicity bin
  TH2F **      fhXEUeMult  ;                   //![GetNCentrBin()] Trigger particle -UE charged hadron momentum imbalance histogram: multiplicity bin
  TH2F **      fhZTMult  ;                     //![GetNCentrBin()] Trigger particle -charged hadron momentum imbalance histogram: multiplicity bin
  TH2F **      fhZTUeMult  ;                   //![GetNCentrBin()] Trigger particle -UE charged hadron momentum imbalance histogram: multiplicity bin
  
  TH2F *       fhAssocPtBkg;                   //! Trigger pT vs associated pT for background
  TH2F **      fhDeltaPhiDeltaEtaAssocPtBin;   //![fNAssocPtBins*GetNZvertBin()] Difference of charged particle phi and trigger particle  phi as function eta difference, for different associated bins
  TH2F **      fhDeltaPhiAssocPtBin;           //![fNAssocPtBins*GetNZvertBin()] Trigger pT vs dPhi for different associated pt and vz bins
  TH2F **      fhDeltaPhiAssocPtBinDEta08;     //![fNAssocPtBins*GetNZvertBin()] Trigger pT vs dPhi for different associated pt and vz bins for Delta eta > 0.8
  TH2F **      fhDeltaPhiAssocPtBinDEta0 ;     //![fNAssocPtBins*GetNZvertBin()] Trigger pT vs dPhi for different associated pt and vz bins for Delta eta = 0
  TH2F **      fhDeltaPhiAssocPtBinHMPID;      //![fNAssocPtBins*GetNZvertBin()] Trigger pT vs dPhi for different associated pt and vz bins, track with HMPID
  TH2F **      fhDeltaPhiAssocPtBinHMPIDAcc;   //![fNAssocPtBins*GetNZvertBin()] Trigger pT vs dPhi for different associated pt and vz bins, track with HMPIDAcc
  TH2F **      fhDeltaPhiBradAssocPtBin;       //![fNAssocPtBins*GetNZvertBin()] Trigger pT vs dPhi Brad (?) for different associated pt bins
  TH2F *       fhDeltaPhiBrad;                 //! Trigger pT vs dPhi Brad (?) for different associated pt bins
  TH2F **      fhXEAssocPtBin ;                //![fNAssocPtBins] Trigger pT vs xE for different associated pt bins
  TH2F **      fhZTAssocPtBin ;                //![fNAssocPtBins] Trigger pT vs zT for different associated pt bins
  TH2F **      fhXEVZ ;                        //![GetNZvertBin()] Trigger pT vs xE for different vz bins
  TH2F **      fhZTVZ ;                        //![GetNZvertBin()] Trigger pT vs zT for different vz bins

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
  
  // If several UE calculation is on,
  TH2F *       fhDeltaPhiUeLeftNeutral  ;      //! Difference of charged particle from underlying events phi and trigger particle  phi as function of neutral particle pT
  TH2F *       fhXEUeLeftNeutral  ;            //! Trigger particle -underlying neutral hadron momentum imbalance histogram
  TH2F *       fhPtHbpXEUeLeftNeutral  ;       //! Trigger particle -underlying neutral hadron momentum HBP histogram
  TH2F *       fhZTUeLeftNeutral  ;            //! Trigger particle -underlying neutral hadron momentum imbalance histogram
  TH2F *       fhPtHbpZTUeLeftNeutral  ;       //! Trigger particle -underlying neutral hadron momentum HBP histogram
  
  //for decay photon trigger correlation
  TH2F *       fhPtPi0DecayRatio ;             //! for pi0 pt and ratio of decay photon pt
  TH2F *       fhDeltaPhiDecayCharged  ;       //! Difference of charged particle phi and decay trigger
  TH2F *       fhXEDecayCharged ;              //! Trigger particle (decay from pi0)-charged hadron momentum imbalance histogram    
  TH2F *       fhZTDecayCharged ;              //! Trigger particle (decay from pi0)-charged hadron momentum imbalance histogram   

  TH2F *       fhDeltaPhiDecayNeutral  ;       //! Difference of neutral particle phi and decay trigger
  TH2F *       fhXEDecayNeutral ;              //! Trigger particle (decay from pi0)-neutral hadron momentum imbalance histogram  
  TH2F *       fhZTDecayNeutral ;              //! Trigger particle (decay from pi0)-neutral hadron momentum imbalance histogram  

  TH2F **      fhDeltaPhiDecayChargedAssocPtBin;//![fNAssocPtBins*GetNZvertBin()] Tagged as decay Trigger pT vs dPhi for different associated pt bins
  
  // If the data is MC, correlation with generated particles
  // check the origin of the cluster : decay photon (pi0, eta, other), merged photon (pi0),
  // hadron, rest of photons (prompt, FSR, ISR)
  TH1F *       fhMCPtTrigger[7];               //! MC pure pT distribution of trigger particles
  TH2F *       fhMCPhiTrigger[7];              //! MC pure Phi distribution of trigger particles
  TH2F *       fhMCEtaTrigger[7];              //! MC pure Eta distribution of trigger particles
  TH1F *       fhMCPtTriggerNotLeading[7];     //! MC pure pT distribution of trigger not leading particles
  TH2F *       fhMCPhiTriggerNotLeading[7];    //! MC pure Phi distribution of trigger not leading particles
  TH2F *       fhMCEtaTriggerNotLeading[7];    //! MC pure Eta distribution of trigger not leading particles
  TH2F *       fhMCEtaCharged[7];              //! MC pure particles charged primary pt vs eta (both associated)
  TH2F *       fhMCPhiCharged[7];              //! MC pure particles charged primary pt vs phi (both associated)
  TH2F *       fhMCDeltaEtaCharged[7];         //! MC pure particles charged trigger primary pt vs delta eta (associated-trigger)
  TH2F *       fhMCDeltaPhiCharged[7];         //! MC pure particles charged trigger primary pt vs delta phi (associated-trigger)
  TH2F *       fhMCDeltaPhiDeltaEtaCharged[7]; //! MC pure particles charged associated primary pt vs delta phi (associated-trigger), in away side
  TH2F *       fhMCDeltaPhiChargedPt[7];       //! MC pure particles charged delta phi vs delta eta (associated-trigger)
  TH2F *       fhMCPtXECharged[7];             //! MC pure particles charged trigger primary pt vs xE
  TH2F *       fhMCPtXEUeCharged[7];           //! MC pure particles charged trigger primary pt vs xE (underlying event)
  TH2F *       fhMCPtXEUeLeftCharged[7];       //! MC pure particles charged trigger primary pt vs xE (underlying event,left cone)
  TH2F *       fhMCPtHbpXECharged[7];          //! MC pure particles charged trigger primary pt vs ln(1/xE)
  TH2F *       fhMCPtHbpXEUeCharged[7];        //! MC pure particles charged trigger primary pt vs ln(1/xE) (underlying event)
  TH2F *       fhMCPtHbpXEUeLeftCharged[7];    //! MC pure particles charged trigger primary pt vs ln(1/xE) (underlying event, left cone)
  TH1F *       fhMCUePart[7];                  //! MC pure UE particles distribution vs pt trig
  TH2F *       fhMCPtZTCharged[7];             //! MC pure particles charged trigger primary pt vs zT
  TH2F *       fhMCPtZTUeCharged[7];           //! MC pure particles charged trigger primary pt vs zT (underlying event)
  TH2F *       fhMCPtZTUeLeftCharged[7];       //! MC pure particles charged trigger primary pt vs zT (underlying event, left cone)
  TH2F *       fhMCPtHbpZTCharged[7];          //! MC pure particles charged trigger primary pt vs ln(1/zT)
  TH2F *       fhMCPtHbpZTUeCharged[7];        //! MC pure particles charged trigger primary pt vs ln(1/zT) (underlying event)
  TH2F *       fhMCPtHbpZTUeLeftCharged[7];    //! MC pure particles charged trigger primary pt vs ln(1/zT) (underlying event, left cone)
  TH2F *       fhMCPtTrigPout[7];              //! MC pure particles charged trigger primary pt vs pOut
  TH2F *       fhMCPtAssocDeltaPhi[7];         //! MC pure particles charged associated primary pt vs delta phi (associated-trigger)

  // Mixing
  TH1I *       fhNEventsTrigger;               //! number of analyzed triggered events
  TH2F *       fhNtracksMB;                    //! total number of tracks in MB events
  TH2F *       fhNclustersMB;                  //! total number of clusters in MB events
  TH2F *       fhMixDeltaPhiCharged  ;         //! Difference of charged particle phi and trigger particle  phi as function of  trigger particle pT
  TH2F *       fhMixDeltaPhiDeltaEtaCharged  ; //! Difference of charged particle phi and trigger particle  phi as function eta difference
  TH2F *       fhMixXECharged;                 //! xE for mixed event
  TH2F *       fhMixXEUeCharged;               //! xE for mixed event in Ue region
  TH2F *       fhMixHbpXECharged;              //! ln(1/xE) for mixed event
  TH2F **      fhMixDeltaPhiChargedAssocPtBin; //![fNAssocPtBins*GetNZvertBin()] Difference of charged particle phi and trigger particle  phi as function of  trigger particle pT, for different associated bins
  TH2F **      fhMixDeltaPhiChargedAssocPtBinDEta08;   //![fNAssocPtBins*GetNZvertBin()] Difference of charged particle phi and trigger particle  phi as function of  trigger particle pT, for different associated bins, delta eta > 0.8
  TH2F **      fhMixDeltaPhiChargedAssocPtBinDEta0;    //![fNAssocPtBins*GetNZvertBin()] Difference of charged particle phi and trigger particle  phi as function of  trigger particle pT, for different associated bins, delta eta = 0
  TH2F **      fhMixDeltaPhiDeltaEtaChargedAssocPtBin; //![fNAssocPtBins*GetNZvertBin()] Difference of charged particle phi and trigger particle  phi as function eta difference, for different associated bins

  TH1I *       fhEventBin;                     //! Number of triggers in a particular event bin (cen,vz,rp)
  TH1I *       fhEventMixBin;                  //! Number of triggers mixed in a particular bin (cen,vz,rp)
  TH1I *       fhEventMBBin;                   //! Number of MB events in a particular bin (cen,vz,rp)
  
  AliAnaParticleHadronCorrelation(              const AliAnaParticleHadronCorrelation & ph) ; // cpy ctor
  AliAnaParticleHadronCorrelation & operator = (const AliAnaParticleHadronCorrelation & ph) ; // cpy assignment
	
  ClassDef(AliAnaParticleHadronCorrelation,34)
} ;
 

#endif //ALIANAPARTICLEHADRONCORRELATION_H



