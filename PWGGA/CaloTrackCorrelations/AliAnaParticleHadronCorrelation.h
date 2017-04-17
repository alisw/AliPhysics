#ifndef ALIANAPARTICLEHADRONCORRELATION_H
#define ALIANAPARTICLEHADRONCORRELATION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliAnaParticleHadronCorrelation
/// \ingroup CaloTrackCorrelationsAnalysis 
/// \brief Correlate trigger particles (photon, pi0, tracks) and charged tracks: Azimuthal correlations, xE distributions.
///
/// Class that contains the algorithm for the analysis of
/// trigger particle - hadron correlations
/// The trigger particle must be found in a previous analysis like Isolated photon
/// (AliAnaPhoton+AliAnaParticleIsolation), High pT pi0 (AliAnaPi0EbE), charged tracks (AliAnaChargedParticle).
///
/// More information can be found in this [twiki](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations)
/// and particularly in this [section](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations#AliAnaParticleHadronCorrelation).
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
/// \author Yaxian Mao <Yaxian.Mao@cern.ch>, CCNU, first implementation.
/// \author Xiangrong Zhu <Xiangrong.Zhu@cern.ch>, CCNU, mixing implementation.
//_________________________________________________________________________

#include "AliAnaCaloTrackCorrBaseClass.h"
class AliAODPWG4ParticleCorrelation ;

class AliAnaParticleHadronCorrelation : public AliAnaCaloTrackCorrBaseClass {
  
public:
  
  AliAnaParticleHadronCorrelation() ;
  
  virtual ~AliAnaParticleHadronCorrelation() ;
  
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
  
  Bool_t       GetDecayPhotonMomentum   (Int_t indexPhoton1, Int_t indexPhoton2, Int_t idetector);
  
  void         MakeChargedCorrelation   (AliAODPWG4ParticleCorrelation * particle) ;
  
  void         MakeNeutralCorrelation   (AliAODPWG4ParticleCorrelation * particle) ;
  
  void         MakeMCChargedCorrelation (Int_t triggerMCLable, Int_t histoIndex, Bool_t lostDecayPair) ;
  
  void         MakeChargedMixCorrelation(AliAODPWG4ParticleCorrelation * particle) ;
  
  // Filling histogram methods
  
  void         FillChargedAngularCorrelationHistograms  (Float_t ptAssoc,  Float_t ptTrig,      Int_t   assocBin,
                                                         Float_t phiAssoc, Float_t phiTrig,     Float_t deltaPhi,
                                                         Float_t etaAssoc, Float_t etaTrig,
                                                         Int_t   decayTag, Float_t hmpidSignal, Int_t outTOF,
                                                         Int_t   cenbin,   Int_t   mcTag);
  
  void         FillChargedEventMixPool();
  
  Bool_t       FillChargedMCCorrelationHistograms       (Float_t mcAssocPt, Float_t mcAssocPhi, Float_t mcAssocEta,
                                                         Float_t mcTrigPt,  Float_t mcTrigPhi,  Float_t mcTrigEta,
                                                         Int_t histoIndex,  Bool_t  lostDecayPair);
  
  void         FillChargedMomentumImbalanceHistograms   (Float_t ptTrig,   Float_t ptAssoc,
                                                         Float_t deltaPhi, Int_t cenbin, Int_t charge,
                                                         Int_t   assocBin, Int_t decayTag,
                                                         Int_t   outTOF,   Int_t mcTag );
  
  void         FillChargedUnderlyingEventHistograms     (Float_t ptTrig,   Float_t ptAssoc,
                                                         Float_t deltaPhi, Int_t cenbin, Int_t outTOF, Int_t   mcTag);
  
  void         FillChargedUnderlyingEventSidesHistograms(Float_t ptTrig,   Float_t ptAssoc,
                                                         Float_t deltaPhi, Int_t   mcTag);
  
  void         FillDecayPhotonCorrelationHistograms     (Float_t ptAssoc,  Float_t phiAssoc, Bool_t bChargedOrNeutral);
  
  void         FillNeutralEventMixPool();
  
  
  void         FillNeutralUnderlyingEventSidesHistograms(Float_t ptTrig,   Float_t ptAssoc,
                                                         Float_t zT,       Float_t hbpZT,
                                                         Float_t deltaPhi);
  
  void         InvMassHisto(AliAODPWG4ParticleCorrelation * trigger, Int_t mcIndex);
  
  Int_t        GetMCTagHistogramIndex(Int_t tag);
  
  /// For histograms in arrays, index in the array, corresponding to any particle origin.
  enum mcTypes     { kmcPhoton     = 0, kmcPi0      = 1, kmcPi0Decay = 2, kmcEta              = 3, kmcEtaDecay         = 4,
    kmcOtherDecay = 5, kmcElectron = 6, kmcHadron   = 7, kmcPi0DecayLostPair = 8, kmcEtaDecayLostPair = 9} ;
  
  static const Int_t fgkNmcTypes = 10;    ///< Number of MC trigger particles checked when filling MC histograms
  
  
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
  
  Bool_t       IsLeadHadronCutOn()        const { return fSelectLeadingHadronAngle   ; }
  void         SwitchOnLeadHadronSelection()    { fSelectLeadingHadronAngle = kTRUE  ; }
  void         SwitchOffLeadHadronSelection()   { fSelectLeadingHadronAngle = kFALSE ; }
  
  void         SwitchOnFillLeadHadronHistograms() { fFillLeadHadOppositeHisto = kTRUE  ; }
  void         SwitchOffFillLeadHadronHistograms(){ fFillLeadHadOppositeHisto = kFALSE ; }
  
  void         SwitchOnInvariantMassHistograms() { fFillInvMassHisto = kTRUE      ; }
  void         SwitchOffInvariantMassHistograms(){ fFillInvMassHisto = kFALSE     ; }
  
  void         SwitchOnBackgroundBinsPtInConeHistograms() { fFillBkgBinsHisto = kTRUE  ; }
  void         SwitchOffBackgroundBinsPtInConeHistograms(){ fFillBkgBinsHisto = kFALSE ; }
  
  void         SwitchOnBackgroundBinsTaggedDecayPtInConeHistograms() { fFillTaggedDecayHistograms = kTRUE  ; }
  void         SwitchOffBackgroundBinsTaggedDecayPtInConeHistograms(){ fFillTaggedDecayHistograms = kFALSE ; }
  
  // Background bins
  void         SetNBackgroundBins(Int_t n)           { if(n < 19) fNBkgBin = n ; }
  void         SetBackgroundLimits(Int_t i,Float_t l){ if(i <= fNBkgBin) fBkgBinLimit[i] = l; }
  
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
  void         SetNDecayBits(Int_t n)            { fNDecayBits = n               ; }
  void         SetDecayBits(Int_t i, UInt_t bit)
  { if(i < AliNeutralMesonSelection::fgkMaxNDecayBits) fDecayBits[i] = bit     ; }
  
  Bool_t       IsHMPIDCorrelation()        const { return fHMPIDCorrelation      ; }
  void         SwitchOnHMPIDCorrelation()        { fHMPIDCorrelation    = kTRUE  ; }
  void         SwitchOffHMPIDCorrelation()       { fHMPIDCorrelation    = kFALSE ; }
  
  void         SwitchOnFillBradHistograms()      { fFillBradHisto       = kTRUE  ; }
  void         SwitchOffFillBradHistograms()     { fFillBradHisto       = kFALSE ; }
  
  Bool_t       OnlyIsolated()              const { return fSelectIsolated        ; }
  void         SelectIsolated(Bool_t s)          { fSelectIsolated   = s         ; }
  
  void         SetPi0AODBranchName(TString n)    { fPi0AODBranchName = n         ; }
  
  void         SetAODNamepTInConeHisto(TString m){ fAODNamepTInConeHisto = m         ; }
  
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
  
  void         SwitchOnFillTriggerAODWithReferences()   { fFillAODWithReferences = kTRUE  ; }
  void         SwitchOffFillTriggerAODWithReferences()  { fFillAODWithReferences = kFALSE ; }
  
  void         SwitchOnCheckNeutralClustersForLeading() { fCheckLeadingWithNeutralClusters = kTRUE  ; }
  void         SwitchOffCheckNeutralClustersForLeading(){ fCheckLeadingWithNeutralClusters = kFALSE ; }
  
  void         SwitchOnFillEtaGapHistograms()    { fFillEtaGapsHisto    = kTRUE  ; }
  void         SwitchOffFillEtaGapHistograms()   { fFillEtaGapsHisto    = kFALSE ; }
  
  void         SwitchOnFillPtImbalancePerPtABinHistograms()  { fFillMomImbalancePtAssocBinsHisto = kTRUE  ; }
  void         SwitchOffFillPtImbalancePerPtABinHistograms() { fFillMomImbalancePtAssocBinsHisto = kFALSE ; }
  
  void         SetMCGenType(Int_t min = 0, Int_t max = 6) { if(min >= 0 && min < fgkNmcTypes) fMCGenTypeMin = min ;
    if(max >= 0 && max < fgkNmcTypes) fMCGenTypeMax = max ; }
  
private:
  
  Bool_t       fFillAODWithReferences;                   ///<  Add to the trigger particle AOD the reference to the tracks or neutrals in correlation.
  
  Bool_t       fCheckLeadingWithNeutralClusters;         ///<  Compare the trigger candidate to Leading pT with the clusters pT, by default only charged.
  
  Float_t      fMaxAssocPt ;                             ///<  Maximum associated hadron pt.
  Float_t      fMinAssocPt ;                             ///<  Minimum associated hadron pt.
  
  Double_t     fDeltaPhiMaxCut ;                         ///<  Minimum Delta Phi Gamma-Hadron.
  Double_t     fDeltaPhiMinCut ;                         ///<  Maximum Delta Phi Gamma-Hadron.
  
  Bool_t       fSelectIsolated ;                         ///<  Select only trigger particles isolated.
  
  Bool_t       fMakeSeveralUE ;                          ///<  Do analysis for several underlying events contribution.
  
  Double_t     fUeDeltaPhiMaxCut ;                       ///<  Minimum Delta Phi Gamma-Underlying Hadron.
  Double_t     fUeDeltaPhiMinCut ;                       ///<  Maximum Delta Phi Gamma-Underlying Hadron.
  
  TString      fPi0AODBranchName;                        ///<  Name of AOD branch with pi0, not trigger.
  
  TString      fAODNamepTInConeHisto;                    ///<  Name of AOD array to fill pT in cone histograms.
  
  Bool_t       fNeutralCorr ;                            ///<  switch the analysis with neutral particles.
  
  Bool_t       fPi0Trigger ;                             ///<  switch the analysis with decay photon from pi0 trigger.
  
  Bool_t       fDecayTrigger ;                           ///<  switch the analysis with decay photon from photon trigger.
  
  Int_t        fNDecayBits ;                             ///<  in case of study of decay triggers, select the decay bit.
  
  ///  in case of study of decay triggers, select the decay bit
  UInt_t       fDecayBits[AliNeutralMesonSelection::fgkMaxNDecayBits] ;
  
  Int_t        fNBkgBin;                                 ///<  Number of bins on pt content in cone.
  Float_t      fBkgBinLimit[20];                         ///<  Pt bin limits on pt content in the cone.
  
  Bool_t       fMakeAbsoluteLeading ;                    ///<  Requesting absolute leading triggers.
  Bool_t       fMakeNearSideLeading ;                    ///<  Requesting near side leading (+-90º from trigger particle) triggers.
  
  Int_t        fLeadingTriggerIndex ;                    ///<  Store here per event the trigger index, to avoid too many loops.
  
  Bool_t       fHMPIDCorrelation    ;                    ///<  Correlate with particles on HMPID or its acceptance.
  
  Bool_t       fFillBradHisto ;                          ///<  DPhi histograms calculated differently.
  
  Int_t        fNAssocPtBins ;                           ///<  Number of associated pT bins under study.
  
  Float_t      fAssocPtBinLimit[20] ;                    ///<  Associated pT under study.
  
  Bool_t       fCorrelVzBin ;                            ///<  Fill one histogram per vz bin.
  
  /// Containers for tracks in stored events for mixing.
  TList **     fListMixTrackEvents ;                     //![GetNCentrBin()*GetNZvertBin()*GetNRPBin()]
  
  /// Containers for calo clusters in stored events for mixing.
  TList **     fListMixCaloEvents ;                      //![GetNCentrBin()*GetNZvertBin()*GetNRPBin()]
  
  Bool_t       fUseMixStoredInReader;                    ///<  Signal if in the current event the pool was filled.
  
  Bool_t       fFillNeutralEventMixPool;                 ///<  Add clusters to pool if requested.
  
  Float_t      fM02MaxCut   ;                            ///<  Study photon clusters with l0 smaller than cut.
  Float_t      fM02MinCut   ;                            ///<  Study photon clusters with l0 larger than cut.
  
  Bool_t       fSelectLeadingHadronAngle;                ///<  Select events with leading particle within a range.
  
  Bool_t       fFillLeadHadOppositeHisto;                ///<  Fill histograms for leading hadrons in opposite side of trigger.
  
  Float_t      fMinLeadHadPhi;                           ///<  Minimum angle between the trigger and leading hadron.
  Float_t      fMaxLeadHadPhi;                           ///<  Maximum ange between the trigger and leading hadron.
  
  Float_t      fMinLeadHadPt;                            ///<  Minimum pT of leading hadron.
  Float_t      fMaxLeadHadPt;                            ///<  Maximum pT of leading hadron.
  
  Bool_t       fFillEtaGapsHisto;                        ///<  Fill azimuthal correlation histograms in 2 eta gaps, |eta|>0.8 and |eta|<0.01.
  
  Bool_t       fFillMomImbalancePtAssocBinsHisto;        ///<  Momentum imbalance histograms in bins of pT associated.
  
  Bool_t       fFillInvMassHisto;                        ///<  Fill invariant mass histograms for trigger.
  
  Bool_t       fFillBkgBinsHisto;                        ///<  Fill pT in cone in background bins distributions.
  
  Bool_t       fFillTaggedDecayHistograms;               ///<  Fill pT in cone distributions in background bins for decay particles.
  
  Float_t      fDecayTagsM02Cut;                         ///<  Lambda0 cut for decay particles.
  
  Int_t        fMCGenTypeMin;                            ///<  Of the fgkNmcTypes possible types, select those between fMCGenTypeMin and fMCGenTypeMax.
  Int_t        fMCGenTypeMax;                            ///<  Of the fgkNmcTypes possible types, select those between fMCGenTypeMin and fMCGenTypeMax.
  
  TVector3       fTrackVector;                           //!<! Track momentum vector.
  TLorentzVector fMomentum;                              //!<! Trigger momentum.
  TLorentzVector fMomentumIM;                            //!<! Cluster momentum from Invariant mass.
  TLorentzVector fDecayMom1;                             //!<! Decay particle momentum.
  TLorentzVector fDecayMom2;                             //!<! Decay particle momentum.
  
  // Histograms
  
  // Trigger particles
  TH1F *       fhPtTriggerInput;                         //!<! pT distribution of trigger particles before selection.
  TH1F *       fhPtTriggerSSCut;                         //!<! pT distribution of trigger particles after shower shape selection.
  TH1F *       fhPtTriggerIsoCut;                        //!<! pT distribution of trigger particles after isolation cut selection.
  TH1F *       fhPtTriggerFidCut;                        //!<! pT distribution of trigger particles after fiducial selection.
  TH1F *       fhPtTrigger;                              //!<! pT distribution of trigger particles.
  TH1F *       fhPtTriggerVtxBC0;                        //!<! pT distribution of trigger particles when vertex is BC0.
  TH1F *       fhPtTriggerPileUp[7];                     //!<! pT distribution of trigger particles for different pile-up definition.
  TH2F *       fhPtTriggerVzBin;                         //!<! pT distribution of trigger particles vs vz bin.
  TH2F *       fhPtTriggerBin;                           //!<! pT distribution of trigger particles vs mixing bin.
  TH2F *       fhPhiTrigger;                             //!<! phi distribution vs pT of trigger particles.
  TH2F *       fhEtaTrigger;                             //!<! eta distribution vs pT of trigger particles.
  
  TH1F *       fhPtTriggerMC[fgkNmcTypes];               //!<! pT distribution of trigger particles, check the origin of the cluster : "Photon","Pi0","Pi0Decay","EtaDecay","OtherDecay","Electron","Hadron".
  
  TH1F *       fhPtDecayTrigger  [AliNeutralMesonSelection::fgkMaxNDecayBits];              //!<! pT distribution of trigger particles, tagged as decay.
  TH1F *       fhPtDecayTriggerMC[AliNeutralMesonSelection::fgkMaxNDecayBits][fgkNmcTypes]; //!<! pT distribution of trigger particles, tagged as decay, check the origin of the cluster.
  
  TH2F *       fhPtTriggerCentrality;                    //!<! pT distribution of trigger particles vs centrality.
  TH2F *       fhPtTriggerEventPlane;                    //!<! pT distribution of trigger particles vs centrality.
  TH2F *       fhTriggerEventPlaneCentrality;            //!<! Event plane vs centrality for trigger particles.
  
  TH1F *       fhPtTriggerMixed;                         //!<! pT distribution of trigger particles, used in mixing.
  TH2F *       fhPtTriggerMixedVzBin;                    //!<! pT distribution of trigger particles, used in mixing, vs vz bin.
  TH2F *       fhPtTriggerMixedBin;                      //!<! pT distribution of trigger particles vs mixing bin.
  TH2F *       fhPhiTriggerMixed;                        //!<! phi distribution vs pT of trigger particles, used in mixing.
  TH2F *       fhEtaTriggerMixed;                        //!<! eta distribution vs pT of trigger particles, used in mixing.
  
  // Leading hadron in the opposite side of the trigger
  TH2F *       fhPtLeadingOppositeHadron;                //!<! pT trigger : pT distribution of leading hadron oposite to trigger.
  TH2F *       fhPtDiffPhiLeadingOppositeHadron;         //!<! pT trigger : difference phi distribution of leading hadron oposite and trigger.
  TH2F *       fhPtDiffEtaLeadingOppositeHadron;         //!<! pT trigger: difference eta distribution of leading hadron oposite and trigger.
  TH1F *       fhPtNoLeadingOppositeHadron;              //!<! pT trigger for events without opposite hadrons.
  TH2F *       fhEtaPhiNoLeadingOppositeHadron;          //!<! Location of trigger when no hadron is found on the opposite side.
  
  //trigger-charged histograms
  TH2F *       fhDeltaPhiDeltaEtaCharged ;               //!<! Differences of eta and phi between trigger and charged hadrons.
  TH2F *       fhPhiCharged  ;                           //!<! Phi distribution of charged particles.
  TH2F *       fhEtaCharged  ;                           //!<! Eta distribution of charged particles.
  TH2F *       fhDeltaPhiCharged  ;                      //!<! Difference of charged particle phi and trigger particle  phi as function of  trigger. particle pT
  TH2F *       fhDeltaEtaCharged  ;                      //!<! Difference of charged particle eta and trigger particle  eta as function of  trigger. particle pT
  TH2F *       fhDeltaPhiChargedPt  ;                    //!<! Difference of charged particle phi and trigger particle  phi as function of charged. particle pT
  TH2F *       fhDeltaPhiUeChargedPt ;                   //!<! Difference of charged particle from underlying events phi and trigger particle  phi as function of charged particle pT.
  TH1F *       fhUePart;                                 //!<! UE particles distribution vs pt trigger.
  TH2F *       fhXECharged  ;                            //!<! Trigger particle -charged hadron momentum imbalance histogram.
  TH2F *       fhXECharged_Cone2  ;                      //!<! Trigger particle -charged hadron momentum imbalance histogram in cone2 (5pi/6-7pi/6).
  TH2F *       fhXEUeCharged  ;                          //!<! Trigger particle -underlying charged hadron momentum imbalance histogram.
  TH2F *       fhXEUeChargedSmallCone  ;                 //!<! Trigger particle -underlying charged hadron momentum imbalance histogram for small cone [80,100].
  TH2F *       fhXEUeChargedMediumCone  ;                //!<! Trigger particle -underlying charged hadron momentum imbalance histogram for medium cone [70,110].
  TH2F *       fhXEUeChargedLargeCone  ;                 //!<! Trigger particle -underlying charged hadron momentum imbalance histogram for large cone [60,120].
  TH2F *       fhXEPosCharged  ;                         //!<! Trigger particle -positive charged hadron momentum imbalance histogram.
  TH2F *       fhXENegCharged  ;                         //!<! Trigger particle -negative charged hadron momentum imbalance histogram.
  TH2F *       fhPtHbpXECharged  ;                       //!<! Trigger particle -charged hadron momentum HBP histogram.
  TH2F *       fhPtHbpXECharged_Cone2  ;                 //!<! Trigger particle -charged hadron momentum HBP histogram in cone2 (5pi/6-7pi/6).
  TH2F *       fhPtHbpXEUeCharged  ;                     //!<! Trigger particle -underlying charged hadron momentum HBP histogram.
  TH2F *       fhZTCharged  ;                            //!<! Trigger particle -charged hadron momentum imbalance histogram.
  TH2F *       fhZTUeCharged  ;                          //!<! Trigger particle -underlying charged hadron momentum imbalance histogram.
  TH2F *       fhZTPosCharged  ;                         //!<! Trigger particle -positive charged hadron momentum imbalance histogram.
  TH2F *       fhZTNegCharged  ;                         //!<! Trigger particle -negative charged hadron momentum imbalance histogram.
  TH2F *       fhPtHbpZTCharged  ;                       //!<! Trigger particle -charged hadron momentum HBP histogram.
  TH2F *       fhPtHbpZTUeCharged  ;                     //!<! Trigger particle -underlying charged hadron momentum HBP histogram.
  
  TH2F *       fhXEChargedMC[fgkNmcTypes]  ;             //!<! Trigger particle -charged hadron momentum imbalance histogram, check the origin of the cluster : decay photon (pi0, eta, other), merged photon (pi0), hadron, rest of photons (prompt, FSR, ISR).
  TH2F *       fhDeltaPhiChargedMC[fgkNmcTypes];         //!<! Trigger particle -charged hadron delta phi histogram, check the origin of the cluster : decay photon (pi0, eta, other), merged photon (pi0), hadron, rest of photons (prompt, FSR, ISR).
  
  TH2F *       fhXEUeChargedRightMC[fgkNmcTypes];        //!<! Trigger particle -underlying hadron momentum imbalance histogram for UE in right cone, check the origin of the cluster : decay photon (pi0, eta, other), merged photon (pi0), hadron, rest of photons (prompt, FSR, ISR).
  TH2F *       fhXEUeChargedLeftMC[fgkNmcTypes];         //!<! Trigger particle -underlying hadron momentum imbalance histogram for UE in left cone, check the origin of the cluster : decay photon (pi0, eta, other), merged photon (pi0), hadron, rest of photons (prompt, FSR, ISR).
  
  TH2F *       fhDeltaPhiDeltaEtaChargedPtA3GeV;         //!<! differences of eta and phi between trigger and charged hadrons, pTa > 3 GeV/c.
  TH2F *       fhDeltaPhiChargedPtA3GeV  ;               //!<! Difference of charged particle phi and trigger particle  phi as function of  trigger particle pT, pTa > 3 GeV/c.
  TH2F *       fhDeltaEtaChargedPtA3GeV  ;               //!<! Difference of charged particle eta and trigger particle  eta as function of  trigger particle pT, pTa > 3 GeV/c.
  
  // Events tagged as pileup by SDD,EMCal, or combination
  TH2F *       fhDeltaPhiChargedPileUp[7]  ;             //!<! Difference of charged particle phi and trigger particle  phi as function of  trigger particle pT.
  TH2F *       fhDeltaEtaChargedPileUp[7]  ;             //!<! Difference of charged particle eta and trigger particle  eta as function of  trigger particle pT.
  TH2F *       fhDeltaPhiChargedPtA3GeVPileUp[7] ;       //!<! Difference of charged particle phi and trigger particle  phi as function of  trigger particle pT, pTa > 3 GeV/c.
  TH2F *       fhDeltaEtaChargedPtA3GeVPileUp[7] ;       //!<! Difference of charged particle eta and trigger particle  eta as function of  trigger particle pT, pTa > 3 GeV/c.
  TH2F *       fhXEChargedPileUp[7]  ;                   //!<! Trigger particle -charged hadron momentum imbalance histogram.
  TH2F *       fhXEUeChargedPileUp[7]  ;                 //!<! Trigger particle -charged hadron momentum imbalance histogram.
  TH2F *       fhZTChargedPileUp[7]  ;                   //!<! Trigger particle -charged hadron momentum imbalance histogram.
  TH2F *       fhZTUeChargedPileUp[7]  ;                 //!<! Trigger particle -charged hadron momentum imbalance histogram.
  TH2F *       fhPtTrigChargedPileUp[7] ;                //!<! trigger and correlated particl pt, to be used for mean value for kT/c.
  
  TH2F *       fhDeltaPhiChargedOtherBC  ;               //!<! Difference of charged particle phi and trigger particle  phi as function of  trigger particle pT.
  TH2F *       fhDeltaPhiChargedPtA3GeVOtherBC ;         //!<! Difference of charged particle phi and trigger particle  phi as function of  trigger particle pT, pTa > 3 GeV/c.
  TH2F *       fhXEChargedOtherBC  ;                     //!<! Trigger particle -charged hadron momentum imbalance histogram.
  TH2F *       fhXEUeChargedOtherBC  ;                   //!<! Trigger particle -charged hadron momentum imbalance histogram.
  TH2F *       fhZTChargedOtherBC  ;                     //!<! Trigger particle -charged hadron momentum imbalance histogram.
  TH2F *       fhZTUeChargedOtherBC  ;                   //!<! Trigger particle -charged hadron momentum imbalance histogram.
  TH2F *       fhPtTrigChargedOtherBC ;                  //!<! trigger and correlated particl pt, to be used for mean value for kT.
  
  TH2F *       fhDeltaPhiChargedBC0  ;                   //!<! Difference of charged particle phi and trigger particle  phi as function of  trigger particle pT.
  TH2F *       fhDeltaPhiChargedPtA3GeVBC0 ;             //!<! Difference of charged particle phi and trigger particle  phi as function of  trigger particle pT, pTa > 3 GeV/c.
  TH2F *       fhXEChargedBC0  ;                         //!<! Trigger particle -charged hadron momentum imbalance histogram.
  TH2F *       fhXEUeChargedBC0  ;                       //!<! Trigger particle -charged hadron momentum imbalance histogram.
  TH2F *       fhZTChargedBC0  ;                         //!<! Trigger particle -charged hadron momentum imbalance histogram.
  TH2F *       fhZTUeChargedBC0  ;                       //!<! Trigger particle -charged hadron momentum imbalance histogram.
  TH2F *       fhPtTrigChargedBC0 ;                      //!<! trigger and correlated particl pt, to be used for mean value for kT.
  
  TH2F *       fhDeltaPhiChargedVtxBC0  ;                //!<! Difference of charged particle phi and trigger particle  phi as function of  trigger particle pT.
  TH2F *       fhDeltaPhiChargedPtA3GeVVtxBC0 ;          //!<! Difference of charged particle phi and trigger particle  phi as function of  trigger particle pT, pTa > 3 GeV/c.
  TH2F *       fhXEChargedVtxBC0  ;                      //!<! Trigger particle -charged hadron momentum imbalance histogram.
  TH2F *       fhXEUeChargedVtxBC0  ;                    //!<! Trigger particle -charged hadron momentum imbalance histogram.
  TH2F *       fhZTChargedVtxBC0  ;                      //!<! Trigger particle -charged hadron momentum imbalance histogram.
  TH2F *       fhZTUeChargedVtxBC0  ;                    //!<! Trigger particle -charged hadron momentum imbalance histogram.
  TH2F *       fhPtTrigChargedVtxBC0 ;                   //!<! trigger and correlated particl pt, to be used for mean value for kT.
  
  // If several UE calculation is on, most useful for jet-jet events contribution
  TH2F *       fhDeltaPhiUeLeftCharged  ;                //!<! Difference of charged particle from underlying events phi and trigger particle  phi as function of charged particle pT
  TH2F *       fhDeltaPhiUeLeftUpCharged;                //!<! Difference of charged particle from underlying events phi and trigger particle  phi
  TH2F *       fhDeltaPhiUeRightUpCharged;               //!<! Difference of charged particle from underlying events phi and trigger particle  phi
  TH2F *       fhDeltaPhiUeLeftDownCharged;              //!<! Difference of charged particle from underlying events phi and trigger particle  phi
  TH2F *       fhDeltaPhiUeRightDownCharged;             //!<! Difference of charged particle from underlying events phi and trigger particle  phi
  TH2F *       fhXEUeLeftCharged  ;                      //!<! Trigger particle -underlying charged hadron momentum imbalance histogram
  TH2F *       fhXEUeLeftUpCharged  ;                    //!<! Trigger particle -underlying charged hadron momentum imbalance histogram
  TH2F *       fhXEUeRightUpCharged ;                    //!<! Trigger particle -underlying charged hadron momentum imbalance histogram
  TH2F *       fhXEUeLeftDownCharged  ;                  //!<! Trigger particle -underlying charged hadron momentum imbalance histogram
  TH2F *       fhXEUeRightDownCharged ;                  //!<! Trigger particle -underlying charged hadron momentum imbalance histogram
  TH2F *       fhPtHbpXEUeLeftCharged  ;                 //!<! Trigger particle -underlying charged hadron momentum HBP histogram
  TH2F *       fhZTUeLeftCharged  ;                      //!<! Trigger particle -underlying charged hadron momentum imbalance histogram
  TH2F *       fhPtHbpZTUeLeftCharged  ;                 //!<! Trigger particle -underlying charged hadron momentum HBP histogram
  
  // For pout and kt extraction
  TH2F *       fhPtTrigPout  ;                           //!<! Pout =associated pt*sin(delta phi) distribution vs trigger pt
  TH2F *       fhPtTrigCharged ;                         //!<! trigger and correlated particl pt, to be used for mean value for kt
  
  /// Differences of phi between trigger and charged hadrons: multiplicity bin
  TH2F **      fhDeltaPhiChargedMult ;                   //![GetNCentrBin()]
  
  /// Differences of eta between trigger and charged hadrons: multiplicity bin.
  TH2F **      fhDeltaEtaChargedMult ;                   //![GetNCentrBin()]
  
  /// Trigger particle -charged hadron momentum imbalance histogram: multiplicity bin.
  TH2F **      fhXEMult  ;                               //![GetNCentrBin()]
  
  /// Trigger particle -UE charged hadron momentum imbalance histogram: multiplicity bin.
  TH2F **      fhXEUeMult  ;                             //![GetNCentrBin()]
  
  /// Trigger particle -charged hadron momentum imbalance histogram: multiplicity bin.
  TH2F **      fhZTMult  ;                               //![GetNCentrBin()]
  
  /// Trigger particle -UE charged hadron momentum imbalance histogram: multiplicity bin.
  TH2F **      fhZTUeMult  ;                             //![GetNCentrBin()]
  
  TH2F *       fhAssocPtBkg;                             //!<! Trigger pT vs associated pT for background.
  
  /// Difference of charged particle phi and trigger particle  phi as function eta difference, for different associated bins.
  TH2F **      fhDeltaPhiDeltaEtaAssocPtBin;             //![fNAssocPtBins*GetNZvertBin()]
  
  /// Trigger pT vs dPhi for different associated pt and vz bins.
  TH2F **      fhDeltaPhiAssocPtBin;                     //![fNAssocPtBins*GetNZvertBin()]
  
  /// Trigger pT vs dPhi for different associated pt and vz bins for Delta eta > 0.8.
  TH2F **      fhDeltaPhiAssocPtBinDEta08;               //![fNAssocPtBins*GetNZvertBin()]
  
  /// Trigger pT vs dPhi for different associated pt and vz bins for Delta eta = 0.
  TH2F **      fhDeltaPhiAssocPtBinDEta0 ;               //![fNAssocPtBins*GetNZvertBin()]
  
  /// Trigger pT vs dPhi for different associated pt and vz bins, track with HMPID.
  TH2F **      fhDeltaPhiAssocPtBinHMPID;                //![fNAssocPtBins*GetNZvertBin()]
  
  /// Trigger pT vs dPhi for different associated pt and vz bins, track with HMPIDAcc.
  TH2F **      fhDeltaPhiAssocPtBinHMPIDAcc;             //![fNAssocPtBins*GetNZvertBin()]
  
  /// Trigger pT vs dPhi Brad (?) for different associated pt bins.
  TH2F **      fhDeltaPhiBradAssocPtBin;                 //![fNAssocPtBins*GetNZvertBin()]
  
  TH2F *       fhDeltaPhiBrad;                           //!<! Trigger pT vs dPhi Brad (?) for different associated pt bins.
  
  /// Trigger pT vs xE for different associated pt bins.
  TH2F **      fhXEAssocPtBin ;                          //![fNAssocPtBins]
  
  /// Trigger pT vs zT for different associated pt bins.
  TH2F **      fhZTAssocPtBin ;                          //![fNAssocPtBins]
  
  /// Trigger pT vs xE for different vz bins.
  TH2F **      fhXEVZ ;                                  //![GetNZvertBin()]
  
  /// Trigger pT vs zT for different vz bins.
  TH2F **      fhZTVZ ;                                  //![GetNZvertBin()]
  
  // Trigger-neutral histograms
  TH2F *       fhDeltaPhiDeltaEtaNeutral ;               //!<! Differences of eta and phi between trigger and neutral hadrons (pi0)
  TH2F *       fhPhiNeutral   ;                          //!<! Phi distribution of neutral particles
  TH2F *       fhEtaNeutral   ;                          //!<! Eta distribution of neutral particles
  TH2F *       fhDeltaPhiNeutral   ;                     //!<! Difference of neutral particle phi and trigger particle  phi as function of  trigger particle pT
  TH2F *       fhDeltaEtaNeutral  ;                      //!<! Difference of neutral particle eta and trigger particle  eta as function of  trigger particle pT
  TH2F *       fhDeltaPhiNeutralPt  ;                    //!<! Difference of neutral particle phi and trigger particle  phi as function of neutral particle particle pT
  TH2F *       fhDeltaPhiUeNeutralPt ;                   //!<! Difference of neutral particle phi and trigger particle  phi as function of neutral particle particle pT
  TH2F *       fhXENeutral  ;                            //!<! Trigger particle - neutral hadron momentum imbalance histogram
  TH2F *       fhXEUeNeutral  ;                          //!<! Trigger particle - neutral hadron momentum imbalance histogram
  TH2F *       fhPtHbpXENeutral  ;                       //!<! Trigger particle - neutral particle momentum HBP histogram
  TH2F *       fhPtHbpXEUeNeutral  ;                     //!<! Trigger particle - underlying neutral hadron momentum HBP histogram
  TH2F *       fhZTNeutral  ;                            //!<! Trigger particle - neutral hadron momentum imbalance histogram
  TH2F *       fhZTUeNeutral  ;                          //!<! Trigger particle - neutral hadron momentum imbalance histogram
  TH2F *       fhPtHbpZTNeutral  ;                       //!<! Trigger particle - neutral particle momentum HBP histogram
  TH2F *       fhPtHbpZTUeNeutral  ;                     //!<! Trigger particle - underlying neutral hadron momentum HBP histogram
  
  // If several UE calculation is on,
  TH2F *       fhDeltaPhiUeLeftNeutral  ;                //!<! Difference of charged particle from underlying events phi and trigger particle  phi as function of neutral particle pT
  TH2F *       fhXEUeLeftNeutral  ;                      //!<! Trigger particle -underlying neutral hadron momentum imbalance histogram
  TH2F *       fhPtHbpXEUeLeftNeutral  ;                 //!<! Trigger particle -underlying neutral hadron momentum HBP histogram
  TH2F *       fhZTUeLeftNeutral  ;                      //!<! Trigger particle -underlying neutral hadron momentum imbalance histogram
  TH2F *       fhPtHbpZTUeLeftNeutral  ;                 //!<! Trigger particle -underlying neutral hadron momentum HBP histogram
  
  // Pi0/Eta trigger correlation, recover input photons
  TH2F *       fhPtPi0DecayRatio ;                       //!<! for pi0 trigger pt and ratio of decay photon pt
  TH2F *       fhDeltaPhiPi0DecayCharged  ;              //!<! Difference of charged particle phi and decay photon from pi0/eta trigger
  TH2F *       fhXEPi0DecayCharged ;                     //!<! Trigger particle (decay from pi0/eta trigger)-charged hadron momentum imbalance histogram
  TH2F *       fhZTPi0DecayCharged ;                     //!<! Trigger particle (decay from pi0/eta trigger)-charged hadron momentum imbalance histogram
  
  TH2F *       fhDeltaPhiPi0DecayNeutral  ;              //!<! Difference of neutral particle phi and decay photon from pi0/eta trigger
  TH2F *       fhXEPi0DecayNeutral ;                     //!<! Trigger particle (decay from pi0/eta trigger)-neutral hadron momentum imbalance histogram
  TH2F *       fhZTPi0DecayNeutral ;                     //!<! Trigger particle (decay from pi0/eta trigger)-neutral hadron momentum imbalance histogram
  
  // Decay photon trigger correlation
  TH2F *       fhDeltaPhiDecayCharged[AliNeutralMesonSelection::fgkMaxNDecayBits]; //!<! Difference of charged particle phi and photon decay trigger.
  TH2F *       fhXEDecayCharged[AliNeutralMesonSelection::fgkMaxNDecayBits];       //!<! Trigger particle (decay from pi0)-charged hadron momentum imbalance histogram.
  TH2F *       fhZTDecayCharged[AliNeutralMesonSelection::fgkMaxNDecayBits];       //!<! Trigger particle (decay from pi0)-charged hadron momentum imbalance histogram.
  
  /// Tagged as decay (fDecayBits[0]) Trigger pT vs dPhi for different associated pt bins
  TH2F **      fhDeltaPhiDecayChargedAssocPtBin;         //![fNAssocPtBins*GetNZvertBin()]
  
  // If the data is MC, correlation with generated particles
  // check the origin of the cluster : decay photon (pi0, eta, other), merged photon (pi0),
  // hadron, rest of photons (prompt, FSR, ISR)
  TH1F *       fhMCPtTrigger[fgkNmcTypes];               //!<! MC pure pT distribution of trigger particles
  TH2F *       fhMCPhiTrigger[fgkNmcTypes];              //!<! MC pure Phi distribution of trigger particles
  TH2F *       fhMCEtaTrigger[fgkNmcTypes];              //!<! MC pure Eta distribution of trigger particles
  TH1F *       fhMCPtTriggerNotLeading[fgkNmcTypes];     //!<! MC pure pT distribution of trigger not leading particles
  TH2F *       fhMCPhiTriggerNotLeading[fgkNmcTypes];    //!<! MC pure Phi distribution of trigger not leading particles
  TH2F *       fhMCEtaTriggerNotLeading[fgkNmcTypes];    //!<! MC pure Eta distribution of trigger not leading particles
  TH2F *       fhMCEtaCharged[fgkNmcTypes];              //!<! MC pure particles charged primary pt vs eta (both associated)
  TH2F *       fhMCPhiCharged[fgkNmcTypes];              //!<! MC pure particles charged primary pt vs phi (both associated)
  TH2F *       fhMCDeltaEtaCharged[fgkNmcTypes];         //!<! MC pure particles charged trigger primary pt vs delta eta (associated-trigger)
  TH2F *       fhMCDeltaPhiCharged[fgkNmcTypes];         //!<! MC pure particles charged trigger primary pt vs delta phi (associated-trigger)
  TH2F *       fhMCDeltaPhiDeltaEtaCharged[fgkNmcTypes]; //!<! MC pure particles charged associated primary pt vs delta phi (associated-trigger), in away side
  TH2F *       fhMCDeltaPhiChargedPt[fgkNmcTypes];       //!<! MC pure particles charged delta phi vs delta eta (associated-trigger)
  TH2F *       fhMCPtXECharged[fgkNmcTypes];             //!<! MC pure particles charged trigger primary pt vs xE
  TH2F *       fhMCPtXEUeCharged[fgkNmcTypes];           //!<! MC pure particles charged trigger primary pt vs xE (underlying event)
  TH2F *       fhMCPtXEUeLeftCharged[fgkNmcTypes];       //!<! MC pure particles charged trigger primary pt vs xE (underlying event,left cone)
  TH2F *       fhMCPtHbpXECharged[fgkNmcTypes];          //!<! MC pure particles charged trigger primary pt vs ln(1/xE)
  TH2F *       fhMCPtHbpXEUeCharged[fgkNmcTypes];        //!<! MC pure particles charged trigger primary pt vs ln(1/xE) (underlying event)
  TH2F *       fhMCPtHbpXEUeLeftCharged[fgkNmcTypes];    //!<! MC pure particles charged trigger primary pt vs ln(1/xE) (underlying event, left cone)
  TH1F *       fhMCUePart[fgkNmcTypes];                  //!<! MC pure UE particles distribution vs pt trig
  TH2F *       fhMCPtZTCharged[fgkNmcTypes];             //!<! MC pure particles charged trigger primary pt vs zT
  TH2F *       fhMCPtZTUeCharged[fgkNmcTypes];           //!<! MC pure particles charged trigger primary pt vs zT (underlying event)
  TH2F *       fhMCPtZTUeLeftCharged[fgkNmcTypes];       //!<! MC pure particles charged trigger primary pt vs zT (underlying event, left cone)
  TH2F *       fhMCPtHbpZTCharged[fgkNmcTypes];          //!<! MC pure particles charged trigger primary pt vs ln(1/zT)
  TH2F *       fhMCPtHbpZTUeCharged[fgkNmcTypes];        //!<! MC pure particles charged trigger primary pt vs ln(1/zT) (underlying event)
  TH2F *       fhMCPtHbpZTUeLeftCharged[fgkNmcTypes];    //!<! MC pure particles charged trigger primary pt vs ln(1/zT) (underlying event, left cone)
  TH2F *       fhMCPtTrigPout[fgkNmcTypes];              //!<! MC pure particles charged trigger primary pt vs pOut
  TH2F *       fhMCPtAssocDeltaPhi[fgkNmcTypes];         //!<! MC pure particles charged associated primary pt vs delta phi (associated-trigger)
  
  // Mixing
  TH1I *       fhNEventsTrigger;                         //!<! Number of analyzed triggered events.
  TH2F *       fhNtracksMB;                              //!<! Total number of tracks in MB events.
  TH2F *       fhNclustersMB;                            //!<! Total number of clusters in MB events.
  TH2F *       fhMixDeltaPhiCharged;                     //!<! Difference of charged particle phi and trigger particle  phi as function of  trigger particle pT.
  TH2F *       fhMixDeltaPhiDeltaEtaCharged;             //!<! Difference of charged particle phi and trigger particle  phi as function eta difference
  TH2F *       fhMixXECharged;                           //!<! xE for mixed event.
  TH2F *       fhMixXEUeCharged;                         //!<! xE for mixed event in Ue region.
  TH2F *       fhMixHbpXECharged;                        //!<! ln(1/xE) for mixed event.
  
  /// Difference of charged particle phi and trigger particle  phi as function of  trigger particle pT, for different associated bins.
  TH2F **      fhMixDeltaPhiChargedAssocPtBin;           //![fNAssocPtBins*GetNZvertBin()]
  
  /// Difference of charged particle phi and trigger particle  phi as function of  trigger particle pT, for different associated bins, delta eta > 0.8.
  TH2F **      fhMixDeltaPhiChargedAssocPtBinDEta08;     //![fNAssocPtBins*GetNZvertBin()]
  
  /// Difference of charged particle phi and trigger particle  phi as function of  trigger particle pT, for different associated bins, delta eta = 0.
  TH2F **      fhMixDeltaPhiChargedAssocPtBinDEta0;      //![fNAssocPtBins*GetNZvertBin()]
  
  /// Difference of charged particle phi and trigger particle  phi as function eta difference, for different associated bins.
  TH2F **      fhMixDeltaPhiDeltaEtaChargedAssocPtBin;   //![fNAssocPtBins*GetNZvertBin()]
  
  TH1I *       fhEventBin;                               //!<! Number of triggers in a particular event bin (cen,vz,rp).
  TH1I *       fhEventMixBin;                            //!<! Number of triggers mixed in a particular bin (cen,vz,rp).
  TH1I *       fhEventMBBin;                             //!<! Number of MB events in a particular bin (cen,vz,rp).
  
  // Check invariant mass
  TH2F *       fhMassPtTrigger;                          //!<! Invariant mass of the trigger.
  TH2F *       fhMCMassPtTrigger[fgkNmcTypes];           //!<! Invariant mass of the trigger vs MC origin.
  
  // pT in isolation cone bins histograms
  
  ///  pT trig distribution for each pT lead in cone bin
  TH1F **       fhPtLeadInConeBin;                       //![fNBkgBin]
  
  ///  pT trig distribution for each pT sum in cone bin
  TH1F **       fhPtSumInConeBin ;                       //![fNBkgBin]
  
  /// pT trig distribution for each pT lead in cone bin for decay particles
  TH1F **       fhPtLeadConeBinDecay;                    //![fNBkgBin*fNDecayBits]
  
  /// pT trig distribution for each pT sum in cone bin for decay particles
  TH1F **       fhSumPtConeBinDecay;                     //![fNBkgBin*fNDecayBits]
  
  /// pT trig distribution for each pT lead in cone bin for MC tag
  TH1F **       fhPtLeadConeBinMC;                       //![fNBkgBin*fgkNmcTypes]
  
  /// pT trig distribution for each pT sum in cone bin for MC tag
  TH1F **       fhSumPtConeBinMC;                        //![fNBkgBin*fgkNmcTypes]
  
  TH2F *        fhTrackResolution;                       //!<! track resolution sigma pT vs pT, away side, ESDs.
  
  TH2F *        fhTrackResolutionUE;                     //!<! track resolution sigma pT vs pT, UE side, ESDs.
  
  /// Copy constructor not implemented.
  AliAnaParticleHadronCorrelation(              const AliAnaParticleHadronCorrelation & ph) ;
  
  /// Assignment operator not implemented.
  AliAnaParticleHadronCorrelation & operator = (const AliAnaParticleHadronCorrelation & ph) ;
  
  /// \cond CLASSIMP
  ClassDef(AliAnaParticleHadronCorrelation,36) ;
  /// \endcond
  
} ;

#endif //ALIANAPARTICLEHADRONCORRELATION_H



