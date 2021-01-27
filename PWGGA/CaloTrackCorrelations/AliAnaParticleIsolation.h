#ifndef ALIANAPARTICLEISOLATION_H
#define ALIANAPARTICLEISOLATION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliAnaParticleIsolation
/// \ingroup CaloTrackCorrelationsAnalysis 
/// \brief Select clusters/tracks with low particle environment in their vecinity,
/// isolated within a cone.
///
/// This class takes a particle AOD object with format AliCaloTrackParticle
/// produced by any of the identified particle classes (AliAnaPhoton, AliAnaElectron,
/// AliAnaPi0EbE, AliAnaChargedParticle) and checks if there is low particle environment
/// around it with the utils of AliIsolationCut, declaring the particle AOD object as isolated or not.
///
/// Class created from old AliPHOSGammaJet (see AliRoot versions previous Release 4-09).
///
/// From june 2019, old methods dealing with UE subtraction moved to AliIsolationCut. 
/// Also removal of studies to use cells as cone input.
///
/// More information can be found in this [twiki](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations)
/// and particularly in this [section](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations#AliAnaParticleIsolation).
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
///_________________________________________________________________________

// --- ROOT system ---
class TH2F;
class TH3F;
class TList ;
class TObjString;

// --- ANALYSIS system ---
#include "AliAnaCaloTrackCorrBaseClass.h"
class AliCaloTrackParticle;
class AliCaloTrackParticleCorrelation ;

class AliAnaParticleIsolation : public AliAnaCaloTrackCorrBaseClass {

 public:
    
  AliAnaParticleIsolation() ;
    
  /// Virtual destructor.
  virtual ~AliAnaParticleIsolation() { ; }

  // Main general methods
  
  TObjString * GetAnalysisCuts() ;
  
  TList      * GetCreateOutputObjects() ;
  
  void         Init() ;

  void         InitParameters() ;
  
  void         MakeAnalysisFillAOD()  ;
  
  void         MakeAnalysisFillHistograms() ;
  
  void         Print( const Option_t * opt ) const ;
 
  // Analysis specific methods
  
  void         FillPileUpHistograms(AliCaloTrackParticleCorrelation* pCandidate) ; //Int_t clusterID) ;
  
  void         FillAcceptanceHistograms();
  
  void         FillShowerShapeControlHistograms(AliCaloTrackParticleCorrelation  * pCandidate,
                                                Int_t mcIndex, Int_t noverlaps) ;
  
  void         FillTrackMatchingControlHistograms(AliCaloTrackParticleCorrelation  * pCandidate,
                                                  Int_t mcIndex) ;
  
  Bool_t       IsTriggerTheNearSideEventLeadingParticle(Int_t & idLeading);
    
  void         StudyEMCALRegions(Float_t pt, Float_t phi, Float_t eta, Float_t m02, 
                                 Float_t coneptsumTrack, Float_t coneptsumCluster, 
                                 Bool_t isolated, Int_t iSM) ;
  
  void         StudyMCConversionRadius(Float_t  pt, Bool_t isolated, Int_t iSM, 
                                       Float_t m02, Int_t     mcTag, Int_t label) ;
  
  void         StudyClustersInCone  (AliCaloTrackParticleCorrelation * aodParticle) ;
  
  void         StudyTracksInCone    (AliCaloTrackParticleCorrelation * aodParticle) ;
  
  void         StudyClustersUEInCone(AliCaloTrackParticleCorrelation * aodParticle) ;
  
  void         StudyTracksUEInCone  (AliCaloTrackParticleCorrelation * aodParticle) ;
  
  // Analysis Setters and Getters
  
  TString      GetTriggerDetectorString()      const { return fIsoDetectorString ; }
  TString      GetTriggerDetector()            const { return fIsoDetector       ; }
  
  Int_t        GetMCIndex(Int_t mcTag);
  
  void         SetTriggerDetector(TString det)     ;
  void         SetTriggerDetector(Int_t   det)     ;

  void         SetMinCellsAngleOverlap(Float_t n)    { fMinCellsAngleOverlap = n ; }
  
  void         SwitchOnTMHistoFill()                 { fFillTMHisto   = kTRUE    ; }
  void         SwitchOffTMHistoFill()                { fFillTMHisto   = kFALSE   ; }
  
  void         SwitchOnSSHistoFill()                 { fFillSSHisto   = kTRUE    ; }
  void         SwitchOffSSHistoFill()                { fFillSSHisto   = kFALSE   ; }

  void         SwitchOnFillHistogramsPerSM()         { fFillPerSMHistograms = kTRUE  ; }
  void         SwitchOffFillHistogramsPerSM()        { fFillPerSMHistograms = kFALSE ; }  
  
  void         SwitchOnFillHistogramsPerTCardIndex()  { fFillPerTCardIndexHistograms = kTRUE  ; }
  void         SwitchOffFillHistogramsPerTCardIndex() { fFillPerTCardIndexHistograms = kFALSE ; }  
  
  void         SwitchOnFillEMCALRegionHistograms()   { fFillEMCALRegionHistograms = kTRUE  ; }
  void         SwitchOffFillEMCALRegionHistograms()  { fFillEMCALRegionHistograms = kFALSE ; }  
  
  Bool_t       IsLeadingOnlyOn()               const { return fLeadingOnly       ; }
  void         SwitchOnLeadingOnly()                 { fLeadingOnly    = kTRUE   ; }
  void         SwitchOffLeadingOnly()                { fLeadingOnly    = kFALSE  ; }
  
  void         SwitchOnCheckNeutralClustersForLeading() { fCheckLeadingWithNeutralClusters = kTRUE  ; }
  void         SwitchOffCheckNeutralClustersForLeading(){ fCheckLeadingWithNeutralClusters = kFALSE ; }

  void         SwitchOnNLMHistoFill()                { fFillNLMHistograms = kTRUE ; }
  void         SwitchOffNLMHistoFill()               { fFillNLMHistograms = kFALSE; }

  void         SwitchOnOnlyTH3HistoFill()            { fFillOnlyTH3Histo = kTRUE ; }
  void         SwitchOffOnlyTH3HistoFill()           { fFillOnlyTH3Histo = kFALSE; }
  
  void         SwitchOnIsolationControlHistoFill()   { fFillIsolationControlHistograms = kTRUE ; }
  void         SwitchOffIsolationControlHistoFill()  { fFillIsolationControlHistograms = kFALSE; }
  
  void         SwitchOnDecayTaggedHistoFill()        { fFillTaggedDecayHistograms = kTRUE ; }
  void         SwitchOffDecayTaggedHistoFill()       { fFillTaggedDecayHistograms = kFALSE; }
  void         SetNDecayBits(Int_t n)                { fNDecayBits = n               ; }
  void         SetDecayBits(Int_t i, UInt_t bit)     { if(i < AliNeutralMesonSelection::fgkMaxNDecayBits)
                                                       fDecayBits[i] = bit           ; }

  void         SwitchOnPrimariesInConeSelection()    { fSelectPrimariesInCone = kTRUE ; }
  void         SwitchOffPrimariesInConeSelection()   { fSelectPrimariesInCone = kFALSE; }

  void         SwitchOnPrimariesPi0DecayStudy()      { fMakePrimaryPi0DecayStudy = kTRUE ; }
  void         SwitchOffPrimariesPi0DecayStudy()     { fMakePrimaryPi0DecayStudy = kFALSE; }
  
  void         SwitchOnOverlapHistograms()           { fFillOverlapHistograms = kTRUE ; }
  void         SwitchOffOverlapHistograms()          { fFillOverlapHistograms = kFALSE; }
  
  void         SwitchOnStudyTracksInCone()           { fStudyTracksInCone = kTRUE  ; }
  void         SwitchOffStudyTracksInCone()          { fStudyTracksInCone = kFALSE ; }

  void         SwitchOnStudyMCConversionRadius()     { fStudyMCConversionRadius = kTRUE  ; }
  void         SwitchOffStudyMCConversionRadius()    { fStudyMCConversionRadius = kFALSE ; }

  void         SwitchOnStudyExoticTrigger()          { fStudyExoticTrigger = kTRUE  ; }
  void         SwitchOffStudyExoticTrigger()         { fStudyExoticTrigger = kFALSE ; }

  void         SwitchOnStudyNCellsCut()              { fStudyNCellsCut     = kTRUE  ; }
  void         SwitchOffStudyNCellsCut()             { fStudyNCellsCut     = kFALSE ; }
  
  void         SwitchOnFillTrackOriginHistograms()   { fFillTrackOriginHistograms = kTRUE  ; }
  void         SwitchOffFillTrackOriginHistograms()  { fFillTrackOriginHistograms = kFALSE ; }
  
  // Study of pT cut in cone
  void         SwitchOnStudyPtCutInCone()            { fStudyPtCutInCone = kTRUE ; }
  void         SwitchOffStudyPtCutInCone()           { fStudyPtCutInCone = kFALSE; }

  void         SwitchOnStudyEtaCutInCone()           { fStudyEtaCutInCone = kTRUE ; }
  void         SwitchOffStudyEtaCutInCone()          { fStudyEtaCutInCone = kFALSE; }

  void         SwitchOnStudyRCutInCone()             { fStudyRCutInCone = kTRUE ; }
  void         SwitchOffStudyRCutInCone()            { fStudyRCutInCone = kFALSE; }
  
  void         SetNPtCutInCone(Int_t n)              { if(n < 19) fNPtCutsInCone = n ; }
  void         SetMinPtCutInConeAt(Int_t i,Float_t l){ if(i <= fNPtCutsInCone) fMinPtCutInCone[i] = l; }  
  void         SetMaxPtCutInConeAt(Int_t i,Float_t l){ if(i <= fNPtCutsInCone) fMaxPtCutInCone[i] = l; }

  void         SetNEtaCutInCone(Int_t n)             { if(n < 10) fNEtaCutsInCone = n ; }
  void         SetEtaCutInConeAt(Int_t i,Float_t l)  { if(i <= fNEtaCutsInCone) fEtaCutInCone[i] = l; }

  void         SetNRCutInCone(Int_t n)               { if(n < 10) fNRCutsInCone = n ; }
  void         SetRCutInConeAt(Int_t i,Float_t l)    { if(i <= fNRCutsInCone) fRCutInCone[i] = l; }

  void         SetNNCellsInCandidate(Int_t n)        { if(n < 19) fNNCellsInCandidate= n ; }
  void         SetNCellsInCandidateAt(Int_t i,Int_t l){ if(i <= fNNCellsInCandidate) fNCellsInCandidate[i] = l; }
  
  void         SetNExoCutInCandidate(Int_t n)        { if(n < 19) fNExoCutInCandidate= n ; }
  void         SetExoCutInCandidateAt(Int_t i,Float_t l) { if(i <= fNExoCutInCandidate) fExoCutInCandidate[i] = l; }
  
  void         SetM02CutForSignal    (Float_t min, Float_t max ) { fM02Narrow[0] = min ; fM02Narrow[1] = max; }
  void         SetM02CutForBackground(Float_t min, Float_t max ) { fM02Wide  [0] = min ; fM02Wide  [1] = max; }

  /// For primary histograms in arrays, index in the array, corresponding to a photon origin.
  enum mcPrimTypes { kmcPrimPhoton = 0, kmcPrimPi0Decay = 1, kmcPrimEtaDecay  = 2, kmcPrimOtherDecay  = 3,
                     kmcPrimPrompt = 4, kmcPrimFrag     = 5, kmcPrimISR       = 6,
                     kmcPrimPi0    = 7, kmcPrimEta      = 8                                               } ;
  
  static const Int_t fgkNmcPrimTypes = 9; ///< Number of MC primary particle types used in the analysis in the histogram arrays.
  
  /// For histograms in arrays, index in the array, corresponding to any particle origin.
  enum mcTypes     { kmcPhoton   = 0, kmcPrompt     = 1, kmcFragment         = 2,
                     kmcPi0      = 3, kmcPi0Decay   = 4, kmcPi0DecayLostPair = 10,
                     kmcEta      = 5, kmcEtaDecay   = 6, kmcEtaDecayLostPair = 11,
                     kmcOtherDecay=7, kmcElectron   = 8, kmcHadron           = 9  } ;
  
  static const Int_t fgkNmcTypes = 12; ///< Number of MC type particles originating the clusters used in the analysis in the histogram arrays.

  void      SetMaximumNumberOfMCParticleCases(Int_t n) { 
    if(n > 0 && n < fgkNmcTypes ) fNumberMCParticleCases = n; else fNumberMCParticleCases = fgkNmcTypes; }
  
 private:
  
  Int_t    fIsoDetector ;                             ///<  Candidate particle for isolation detector.
  TString  fIsoDetectorString ;                       ///<  Candidate particle for isolation detector.
  Bool_t   fFillTMHisto;                              ///<  Fill track matching plots.
  Bool_t   fFillSSHisto;                              ///<  Fill Shower shape plots. Activate it only on photon analysis, enables filling of wide/narrow shape histograms.
  Bool_t   fFillPerSMHistograms ;                     ///<  Fill histograms per SM
  Bool_t   fFillPerTCardIndexHistograms ;             ///<  Fill histograms per T-Card index.
  Int_t    fTCardIndex;                               ///<  Store here the T-Card index per trigger cluster.
  Bool_t   fFillEMCALRegionHistograms ;               ///<  Fill histograms in EMCal slices
  Bool_t   fFillOverlapHistograms;                    ///<  Fill histograms that depend on number of overlaps
  Bool_t   fStudyTracksInCone;                        ///<  Study tracks depending on different track info
  Bool_t   fStudyMCConversionRadius;                  ///<  Study shower shape depending the conversion radius
  Bool_t   fFillTrackOriginHistograms;                ///< Fill histograms checking the MC origin of the tracks in cone
  
  Bool_t   fFillTaggedDecayHistograms;                ///<  Fill histograms for clusters tagged as decay.
  Int_t    fNDecayBits ;                              ///<  In case of study of decay triggers, select the decay bit.
  UInt_t   fDecayBits[AliNeutralMesonSelection::fgkMaxNDecayBits] ; ///< In case of study of decay triggers, select the decay. bit
  
  Bool_t   fFillNLMHistograms;                        ///<  Fill NLM histograms.
  Bool_t   fFillOnlyTH3Histo;                         ///< Fill only TH3 histograms when duplication
  Bool_t   fFillIsolationControlHistograms;           ///< Fill control histograms from AliIsolationCut
  
  Bool_t   fLeadingOnly;                              ///<  Do isolation with leading particle.
  Bool_t   fCheckLeadingWithNeutralClusters;          ///<  Compare the trigger candidate to Leading pT with the clusters pT, by default only charged.
  Bool_t   fSelectPrimariesInCone;                    ///<  In primary particle isolation studies, select only particles in isolation cone within detector acceptance and E cut.
  Bool_t   fMakePrimaryPi0DecayStudy;                 ///<  Fill dedicated histograms for primary decay photons.
  
  Float_t  fMinCellsAngleOverlap;                     ///<  Number of cells that define the cluster overlap.
  
  Float_t  fM02Narrow[2];                             ///<  Long axis signal region
  Float_t  fM02Wide  [2];                             ///<  Long axis background region
  
  Int_t    fNumberMCParticleCases;                    ///< Number of histograms per MC particle type, maximum is fgkNmcPrimTypes
  
  Bool_t   fStudyPtCutInCone;                         ///<  Activate study of track/cluster min pT on sum of pT in cone
  Int_t    fNPtCutsInCone;                            ///<  Number of track/cluster min pT cut to test in cone for sum pT calculation.
  Float_t  fMinPtCutInCone[20];                       ///<  List of track/cluster min pT cut to test in cone for sum pT calculation.
  Float_t  fMaxPtCutInCone[20];                       ///<  List of track/cluster max pT cut to test in cone for sum pT calculation.
  
  Float_t  fConeNClusterPerMinCut    [20];            ///< Temporal container of n clusters per pT min cut
  Float_t  fConeptsumClusterPerMinCut[20];            ///< Temporal container of sum pT clusters per pT min cut
  
  Float_t  fConeptsumEtaBandClusterPerMinCut   [20];  ///< Temporal container of n clusterss in eta band per pT min cut
  Float_t  fConeNEtaBandClusterPerMinCut       [20];  ///< Temporal container of sum pT clusters in eta band per pT min cut
  Float_t  fConeptsumClusterSubEtaBandPerMinCut[20];  ///< Temporal container of n clusters in cone minus  eta band per pT min cut
  Float_t  fConeNClusterSubEtaBandPerMinCut    [20];  ///< Temporal container of sum pT clusters in cone minus eta band per pT min cut

  Float_t  fConeptsumPhiBandClusterPerMinCut   [20];  ///< Temporal container of n clusters in phi band per pT min cut
  Float_t  fConeNPhiBandClusterPerMinCut       [20];  ///< Temporal container of sum pT clusters in phi band per pT min cut
  Float_t  fConeptsumClusterSubPhiBandPerMinCut[20];  ///< Temporal container of n clusters in cone minus  phi band per pT min cut
  Float_t  fConeNClusterSubPhiBandPerMinCut    [20];  ///< Temporal container of sum pT clusters in cone minus phi band per pT min cut
  
  Float_t  fConeNTrackPerMinCut      [20];            ///< Temporal container of n tracks per pT min cut
  Float_t  fConeptsumTrackPerMinCut  [20];            ///< Temporal container of sum pT tracks per pT min cut
 
  Float_t  fConeptsumPerpTrackPerMinCut   [20];       ///< Temporal container of n tracks in perpendicular cone per pT min cut
  Float_t  fConeNPerpTrackPerMinCut       [20];       ///< Temporal container of sum pT tracks in perpendicular cone per pT min cut
  Float_t  fConeptsumTrackSubPerpPerMinCut[20];       ///< Temporal container of n tracks in cone minus  perpendicular cone per pT min cut
  Float_t  fConeNTrackSubPerpPerMinCut    [20];       ///< Temporal container of sum pT tracks in cone minus perpendicular cone per pT min cut
  
  Float_t  fConeptsumEtaBandTrackPerMinCut   [20];    ///< Temporal container of n tracks in eta band per pT min cut
  Float_t  fConeNEtaBandTrackPerMinCut       [20];    ///< Temporal container of sum pT tracks in eta band per pT min cut
  Float_t  fConeptsumTrackSubEtaBandPerMinCut[20];    ///< Temporal container of n tracks in cone minus  eta band per pT min cut
  Float_t  fConeNTrackSubEtaBandPerMinCut    [20];    ///< Temporal container of sum pT tracks in cone minus eta band per pT min cut

  Float_t  fConeptsumPhiBandTrackPerMinCut   [20];    ///< Temporal container of n tracks in phi band per pT min cut
  Float_t  fConeNPhiBandTrackPerMinCut       [20];    ///< Temporal container of sum pT tracks in phi band per pT min cut
  Float_t  fConeptsumTrackSubPhiBandPerMinCut[20];    ///< Temporal container of n tracks in cone minus  phi band per pT min cut
  Float_t  fConeNTrackSubPhiBandPerMinCut    [20];    ///< Temporal container of sum pT tracks in cone minus phi band per pT min cut
  
  Bool_t   fStudyEtaCutInCone;                        ///<  Activate study of track/cluster max eta on sum of pT in cone
  Int_t    fNEtaCutsInCone;                           ///<  Number of track/cluster max eta cut to test in cone for sum pT calculation.
  Float_t  fEtaCutInCone[10];                         ///<  List of track/cluster max eta cut to test in cone for sum pT calculation.

  Bool_t   fStudyRCutInCone;                          ///<  Activate study of track/cluster sum of pT in cone with variable size
  Int_t    fNRCutsInCone;                             ///<  Number of track/cluster max R cut to test in cone for sum pT calculation.
  Float_t  fRCutInCone[10];                           ///<  List of track/cluster max R cut to test in cone for sum pT calculation.

  Bool_t   fStudyNCellsCut;                           ///<  Fill histograms with track and cluster pT depending n cells in cluster
  Int_t    fNNCellsInCandidate;                       ///<  Number of cells in cluster selection to test in cone for sum pT calculation.
  Int_t    fNCellsInCandidate[20];                    ///<  List of Number of cells in cluster selection to test in cone for sum pT calculation.
  Int_t    fNCellsWithWeight;                         ///<  number of cells in cluster with enough energy for shower shape, internal
  Int_t    fTrigSupMod;                               ///<  super module number of trigger cluster
 
  Bool_t   fStudyExoticTrigger;                       ///<  Fill histograms with track and cluster pT when the trigger is exotic
  Int_t    fNExoCutInCandidate;                       ///<  Number of exoticity cuts in cluster selection to test in cone for sum pT calculation.
  Float_t  fExoCutInCandidate[20];                    ///<  List of exoticity cuts in cluster selection to test in cone for sum pT calculation.
  
  TLorentzVector fMomentum;                           //!<! Temporary vector, avoid creation per event.
  TLorentzVector fMomIso;                             //!<! Temporary vector, avoid creation per event.
  TLorentzVector fMomDaugh1;                          //!<! Temporary vector, avoid creation per event.
  TLorentzVector fMomDaugh2;                          //!<! Temporary vector, avoid creation per event.
  TVector3       fTrackVector;                        //!<! Temporary vector, avoid creation per event.
  TVector3       fProdVertex;                         //!<! Temporary vector, avoid creation per event.
 
  AliVCluster*   fCluster;                            //!<! Temporary vcluster, avoid creation per event.
  TObjArray  *   fClustersArr;                        //!<! Temporary ClustersArray, avoid creation per event.
  AliVCaloCells* fCaloCells;                          //!<! Temporary AliVCaloCells pointer for selected calorimeter candidate, avoid creation per event.
  Bool_t         fIsExoticTrigger;                    //!<! Trigger cluster considered as exotic
  Float_t        fClusterExoticity;                   //!<! Temporary container or currently analyzed cluster exoticity

  // Histograms  
  
  TH1F *   fhPt[2][2] ;                                //!<! Number of non/isolated narrow/wide particles vs pT.
  TH2F *   fhPtCentrality[2][2] ;                      //!<! Number of non/isolated narrow/wide particles centrality vs pT.
  TH2F *   fhPtEventPlane[2][2] ;                      //!<! Number of non/isolated narrow/wide particles event plane angle vs pT.
  TH2F *   fhPtNLocMax[2][2] ;                         //!<! Number of non/isolated narrow/wide particles vs NLM in cluster.
  TH3F *   fhPtEtaPhi[2][2] ;                          //!<! cluster pt vs eta vs phi of non/isolated narraw/wide particles.
  TH1F *   fhPtExoTrigger[2];                          //!<! Number of non/isolated exotic cluster vs pT.
  
  TH1F *   fhPtDecay       [2][AliNeutralMesonSelection::fgkMaxNDecayBits]; //!<! Number of (non) isolated Pi0 decay particles (invariant mass tag).
  TH2F *   fhEtaPhiDecay   [2][AliNeutralMesonSelection::fgkMaxNDecayBits]; //!<! eta vs phi of (not) isolated leading Pi0 decay particles.
  TH2F *   fhPtLambda0Decay[2][AliNeutralMesonSelection::fgkMaxNDecayBits]; //!<! Shower shape of (non) isolated leading Pi0 decay particles (do not apply SS cut previously).

  TH2F *   fhPtTrackInConeMCPrimary  [4] ;             //!<! Track Pt in the cone for tracks originating from primary charged pions, kaons, protons and else, reconstructed pT.
  TH2F *   fhPtTrackInConeMCSecondary[4] ;             //!<! Track Pt in the cone for tracks originating from secondary charged pions, kaons, protons and else, reconstructed pT.
  TH2F *   fhPtTrackInConeMCPrimaryGener  [4] ;        //!<! Track Pt in the cone for tracks originating from primary charged pions, kaons, protons and else, generated pT.
  TH2F *   fhPtTrackInConeMCSecondaryGener[4] ;        //!<! Track Pt in the cone for tracks originating from secondary charged pions, kaons, protons and else, generated pT.

  TH2F *   fhPtInConeExoTrigger ;                      //!<! Cluster and tracks  Pt in the cone. Trigger is exotic
  TH2F *   fhPtClusterInConeExoTrigger ;               //!<! Clusters Pt in the cone. Trigger is exotic
  TH2F *   fhPtTrackInConeExoTrigger ;                 //!<! Tracks Pt in the cone. Trigger considered exotic
  
  TH2F *   fhPtTrackInConeOtherBCPileUpSPD ;           //!<! Track Pt in the cone, tracks out of main BC Time window.
  TH2F *   fhPtTrackInConeVtxBC0 ;                     //!<! Track Pt in the cone, tracks in BC=0.
  TH2F *   fhPtTrackInConeBC0PileUpSPD ;               //!<! Track Pt in the cone, tracks in BC=0.
  TH2F *   fhPtInConePileUp[7] ;                       //!<! Particle Pt in the cone, if event is from pile-up (SPD method).
  TH2F *   fhPerpConeSumPtTOFBC0 ;                     //!<! Sum Pt in cone at the perpendicular phi region to trigger axis  (phi +90), TOF BC=0
  TH2F *   fhPtInPerpConeTOFBC0 ;                      //!<! Particle Pt  in cone at the perpendicular phi region to trigger axis  (phi +90), TOF BC=0
  TH2F *   fhEtaPhiInPerpConeTOFBC0 ;                  //!<! Eta vs. phi of tracks in perpendicular cone, with TOF BC=0.
  
  TH3F *   fhPtM02SumPtCone;                           //!<! ABCD TH3F histogram Pt, Shower Shape and sum(ET)+sum(pT) cone
  TH3F *   fhPtM02SumPtConeMC[fgkNmcTypes];            //!<! ABCD TH3F histogram Pt, Shower Shape and sum(ET)+sum(pT) cone, per MC particle
 
  TH3F *   fhPtM02SumPtConeCharged;                    //!<! ABCD TH3F histogram Pt, Shower Shape and sum(ET)+sum(pT) cone, charged in cone
  TH3F *   fhPtM02SumPtConeChargedMC[fgkNmcTypes];     //!<! ABCD TH3F histogram Pt, Shower Shape and sum(ET)+sum(pT) cone, per MC particle, charged in cone
  
  /// ABCD TH3F histogram Pt, Shower Shape and sum(ET)+sum(pT) cone vs centrality
  TH3F **  fhPtM02SumPtConeCent;                       //![GetNCentrBin()] 
  
  /// ABCD TH3F histogram Pt, Shower Shape and sum(ET)+sum(pT) cone vs centrality, charged particles in cone
  TH3F **  fhPtM02SumPtConeChargedCent;                //![GetNCentrBin()] 

  /// ABCD TH3F histogram Pt, Shower Shape and sum(ET)+sum(pT) cone vs centrality
  /// Different centrality bins and MC particle type origin
  TH3F **  fhPtM02SumPtConeCentMC;                     //![GetNCentrBin()*fNumberMCParticleCases] 
  
  /// ABCD TH3F histogram Pt, Shower Shape and sum(ET)+sum(pT) cone vs centrality, charged particles in cone
  /// Different centrality bins and MC particle type origin
  TH3F **  fhPtM02SumPtConeChargedCentMC;              //![GetNCentrBin()*fNumberMCParticleCases] 
  

  TH2F *   fhConeSumPtM02Cut[2] ;                      //!<! Cluster and tracks Sum Pt in the cone for wide or narrow clusters
  TH2F *   fhConeSumPtM02CutMC[fgkNmcTypes][2] ;       //!<! Cluster and tracks Sum Pt in the cone for wide or narrow clusters, per MC particle
  TH3F *   fhConeSumPtCentM02Cut[2] ;                  //!<! Cluster and tracks Sum Pt in the cone for wide or narrow clusters vs centrality

  TH2F *   fhConeSumPtExoTrigger ;                     //!<! Cluster and tracks Sum Pt in the cone. Trigger is exotic
  TH2F *   fhConeSumPtClusterExoTrigger ;              //!<! Clusters Sum Pt  in the cone. Trigger is exotic
  TH2F *   fhConeSumPtTrackExoTrigger ;                //!<! Tracks Sum Pt  in the cone. Trigger considered exotic
  
  // MC
  
  TH2F *   fhEtaPrimMC  [fgkNmcPrimTypes];             //!<! Pt vs Eta of generated photon.
  TH2F *   fhPhiPrimMC  [fgkNmcPrimTypes];             //!<! Pt vs Phi of generated photon.
  TH1F *   fhPtPrimMC   [fgkNmcPrimTypes];             //!<! Number of generated photon vs pT.
  TH1F *   fhPtPrimMCiso[fgkNmcPrimTypes];             //!<! Number of generated isolated photon vs pT.

  TH2F *   fhConeSumPtPrimMC          [fgkNmcPrimTypes]; //!<! Number of generated isolated photon vs photon pT vs sum of primaries pT in cone.
  TH2F *   fhConeSumPtChargedPrimMC   [fgkNmcPrimTypes]; //!<! Number of generated isolated photon vs photon pT vs sum of charged primaries pT in cone.
  TH3F *   fhConeSumPtCenPrimMC       [fgkNmcPrimTypes]; //!<! Number of generated isolated photon vs photon pT vs sum of primaries pT in cone vs centrality.
  TH3F *   fhConeSumPtCenChargedPrimMC[fgkNmcPrimTypes]; //!<! Number of generated isolated photon vs photon pT vs sum of charged primaries pT in cone vs centrality.

  TH2F *   fhConeSumPtUESubPrimMC          [fgkNmcPrimTypes]; //!<! Number of generated isolated photon vs photon pT vs sum of primaries pT in cone UE subtracted.
  TH2F *   fhConeSumPtUESubChargedPrimMC   [fgkNmcPrimTypes]; //!<! Number of generated isolated photon vs photon pT vs sum of charged primaries pT in cone  UE subtracted.
  TH3F *   fhConeSumPtUESubCenPrimMC       [fgkNmcPrimTypes]; //!<! Number of generated isolated photon vs photon pT vs sum of primaries pT in cone  UE subtracted vs centrality.
  TH3F *   fhConeSumPtUESubCenChargedPrimMC[fgkNmcPrimTypes]; //!<! Number of generated isolated photon vs photon pT vs sum of charged primaries pT in cone  UE subtracted vs centrality.

  TH2F *   fhConeSumPtNeutralChargedRatioPrimMC        [fgkNmcPrimTypes]; //!<!  photon pT vs ratio of sum of primaries pT in cone neutral over charged.
  TH2F *   fhConeSumPtPerpConeNeutralChargedRatioPrimMC[fgkNmcPrimTypes]; //!<!  photon pT vs ratio of sum of primaries pT in perpendicular cone neutral over charged.
  
  TH2F *   fhConeSumPtPrimMCEmb          [fgkNmcPrimTypes]; //!<! Number of generated isolated photon vs photon pT vs sum of primaries and embedded data pT in cone.
  TH2F *   fhConeSumPtChargedPrimMCEmb   [fgkNmcPrimTypes]; //!<! Number of generated isolated photon vs photon pT vs sum of charged primaries  and embedded data tracks pT in cone.
  TH3F *   fhConeSumPtCenPrimMCEmb       [fgkNmcPrimTypes]; //!<! Number of generated isolated photon vs photon pT vs sum of primaries  and embedded data pT in cone vs centrality.
  TH3F *   fhConeSumPtCenChargedPrimMCEmb[fgkNmcPrimTypes]; //!<! Number of generated isolated photon vs photon pT vs sum of charged primaries  and embedded data  pT in cone vs centrality.

  TH2F *   fhConeSumPtUESubPrimMCEmb          [fgkNmcPrimTypes]; //!<! Number of generated isolated photon vs photon pT vs sum of primaries  and embedded data pT in cone UE subtracted.
  TH2F *   fhConeSumPtUESubChargedPrimMCEmb   [fgkNmcPrimTypes]; //!<! Number of generated isolated photon vs photon pT vs sum of charged  and embedded data primaries pT in cone  UE subtracted.
  TH3F *   fhConeSumPtUESubCenPrimMCEmb       [fgkNmcPrimTypes]; //!<! Number of generated isolated photon vs photon pT vs sum of primaries  and embedded data pT in cone  UE subtracted vs centrality.
  TH3F *   fhConeSumPtUESubCenChargedPrimMCEmb[fgkNmcPrimTypes]; //!<! Number of generated isolated photon vs photon pT vs sum of charged  and embedded data primaries pT in cone  UE subtracted vs centrality.

  TH1F *   fhPtPrimMCPi0DecayPairOutOfCone;            //!<! Pi0 decay photons, with decay pair out of isolation cone.
  TH1F *   fhPtPrimMCPi0DecayPairOutOfAcceptance;      //!<! Pi0 decay photons, with decay pair out of detector acceptance.
  TH1F *   fhPtPrimMCPi0DecayPairOutOfAcceptanceNoOverlap;   //!<! Pi0 decay photons, with decay pair out of detector acceptance.
  TH1F *   fhPtPrimMCPi0DecayPairAcceptInConeLowPt;    //!<! Pi0 decay photons, with decay pair in cone and acceptance and lower pT than threshold.
  TH1F *   fhPtPrimMCPi0DecayPairAcceptInConeLowPtNoOverlap; //!<! Pi0 decay photons, with decay pair in cone and acceptance and lower pT than threshold, and do not overlap.
  TH1F *   fhPtPrimMCPi0DecayPairAcceptInConeLowPtNoOverlapCaloE; //!<! Pi0 decay photons, with decay pair in cone and acceptance and lower pT than threshold, and larger than detector threshold, and do not overlap.
  TH1F *   fhPtPrimMCPi0DecayPairNoOverlap;            //!<! Pi0 decay photons, not overlapped decay.

  TH1F *   fhPtPrimMCPi0DecayIsoPairOutOfCone;         //!<! Pi0 decay photons, with decay pair out of isolation cone, isolated.
  TH1F *   fhPtPrimMCPi0DecayIsoPairOutOfAcceptance;   //!<! Pi0 decay photons, with decay pair out of detector acceptance, isolated.
  TH1F *   fhPtPrimMCPi0DecayIsoPairOutOfAcceptanceNoOverlap; //!<! Pi0 decay photons, with decay pair out of detector acceptance, isolated.
  TH1F *   fhPtPrimMCPi0DecayIsoPairAcceptInConeLowPt; //!<! Pi0 decay photons, with decay pair in cone and acceptance and lower pT than threshold, isolated.
  TH1F *   fhPtPrimMCPi0DecayIsoPairAcceptInConeLowPtNoOverlap;      //!<! Pi0 decay photons, with decay pair in cone and acceptance and lower pT than threshold, and do not overlap, isolated.
  TH1F *   fhPtPrimMCPi0DecayIsoPairAcceptInConeLowPtNoOverlapCaloE; //!<! Pi0 decay photons, with decay pair in cone and acceptance and lower pT than threshold, and larger than detector threshold, and do not overlap, isolated.
  TH1F *   fhPtPrimMCPi0DecayIsoPairNoOverlap;         //!<! Pi0 decay photons isolated, not overlapped decay.

  TH1F *   fhPtPrimMCPi0Overlap;                       //!<! Pi0 with overlapped decay photons.
  TH1F *   fhPtPrimMCPi0IsoOverlap;                    //!<! Pi0 isolated with overlapped decay photons.

  TH1F *   fhPtPrimMCEtaDecayPairOutOfCone;            //!<! Eta decay photons, with decay pair out of isolation cone.
  TH1F *   fhPtPrimMCEtaDecayPairOutOfAcceptance;      //!<! Eta decay photons, with decay pair out of detector acceptance.
  TH1F *   fhPtPrimMCEtaDecayPairOutOfAcceptanceNoOverlap;   //!<! Eta decay photons, with decay pair out of detector acceptance.
  TH1F *   fhPtPrimMCEtaDecayPairAcceptInConeLowPt;    //!<! Eta decay photons, with decay pair in cone and acceptance and lower pT than threshold.
  TH1F *   fhPtPrimMCEtaDecayPairAcceptInConeLowPtNoOverlap; //!<! Eta decay photons, with decay pair in cone and acceptance and lower pT than threshold, and do not overlap.
  TH1F *   fhPtPrimMCEtaDecayPairAcceptInConeLowPtNoOverlapCaloE; //!<! Eta decay photons, with decay pair in cone and acceptance and lower pT than threshold, and larger than detector threshold, and do not overlap.
  TH1F *   fhPtPrimMCEtaDecayPairNoOverlap;            //!<! Eta decay photons, not overlapped decay.
  
  TH1F *   fhPtPrimMCEtaDecayIsoPairOutOfCone;         //!<! Eta decay photons, with decay pair out of isolation cone, isolated.
  TH1F *   fhPtPrimMCEtaDecayIsoPairOutOfAcceptance;   //!<! Eta decay photons, with decay pair out of detector acceptance, isolated.
  TH1F *   fhPtPrimMCEtaDecayIsoPairOutOfAcceptanceNoOverlap; //!<! Eta decay photons, with decay pair out of detector acceptance, isolated.
  TH1F *   fhPtPrimMCEtaDecayIsoPairAcceptInConeLowPt; //!<! Eta decay photons, with decay pair in cone and acceptance and lower pT than threshold, isolated.
  TH1F *   fhPtPrimMCEtaDecayIsoPairAcceptInConeLowPtNoOverlap; //!<! Eta decay photons, with decay pair in cone and acceptance and lower pT than threshold, and do not overlap, isolated.
  TH1F *   fhPtPrimMCEtaDecayIsoPairAcceptInConeLowPtNoOverlapCaloE; //!<! Eta decay photons, with decay pair in cone and acceptance and lower pT than threshold, and larger than detector threshold, and do not overlap, isolated.
  TH1F *   fhPtPrimMCEtaDecayIsoPairNoOverlap;         //!<! Eta decay photons isolated, not overlapped decay.
  
  TH1F *   fhPtPrimMCEtaOverlap;                       //!<! Eta with overlapped decay photons.
  TH1F *   fhPtPrimMCEtaIsoOverlap;                    //!<! Eta isolated with overlapped decay photons.

  TH1F *   fhPtMC      [fgkNmcTypes][2][2];            //!<! Number of not/isolated narrow/wide mcTypes particle.
  TH3F *   fhPtEtaPhiMC[fgkNmcTypes][2];               //!<! pT vs eta vs phi of not/isolated mcTypes particle.
  
  TH1F *   fhPtDecayMC  [2][AliNeutralMesonSelection::fgkMaxNDecayBits][fgkNmcTypes] ; //!<! Number of (not) isolated Pi0 decay particles (invariant mass tag) for a mcTypes particle.
  
  TH2F *   fhPtLambda0MC    [fgkNmcTypes][2];           //!<! Shower shape of (non) isolated candidates originated by mcTypes particle (do not apply SS cut previously).
  TH2F *   fhPtLambda0MCConv[fgkNmcTypes][2];           //!<! Shower shape of (non) isolated candidates originated by mcTypes particle that converted (do not apply SS cut previously).
  TH2F *   fhPtLambda0MCNCellCut[fgkNmcTypes][2];       //!<! Shower shape of (non) isolated candidates originated by mcTypes particle with n cell_w > 4.

  TH2F *   fhPtLambda0MCWith1Overlap    [fgkNmcTypes][2];           //!<! Shower shape of (non) isolated candidates originated by mcTypes particle (do not apply SS cut previously). At least one overlap from other particles.
  TH2F *   fhPtLambda0MCConvWith1Overlap[fgkNmcTypes][2];           //!<! Shower shape of (non) isolated candidates originated by mcTypes particle that converted (do not apply SS cut previously). At least one overlap from other particles.

  TH2F *   fhPtLambda0MCWithNoOverlap    [fgkNmcTypes][2];          //!<! Shower shape of (non) isolated candidates originated by mcTypes particle (do not apply SS cut previously). More tha one overlap from other particles.
  TH2F *   fhPtLambda0MCConvWithNoOverlap[fgkNmcTypes][2];          //!<! Shower shape of (non) isolated candidates originated by mcTypes particle that converted (do not apply SS cut previously). More tha one overlap from other particles.
  
  TH2F *   fhPtNOverlap    [fgkNmcTypes][2];                        //!<! Number of overlaps of (non) isolated candidates originated by mcTypes (do not apply SS cut previously). More tha one overlap from other particles.
  TH2F *   fhPtNOverlapConv[fgkNmcTypes][2];                        //!<! Number of overlaps of (non) isolated candidates originated by mcTypes particle that converted (do not apply SS cut previously). More tha one overlap from other particles.

  // Track matching studies
  TH2F *   fhTrackMatchedDEta[2]     ;                 //!<! Eta distance between track and cluster vs cluster E.
  TH2F *   fhTrackMatchedDPhi[2]     ;                 //!<! Phi distance between track and cluster vs cluster E.
  TH2F *   fhTrackMatchedDEtaDPhi[2] ;                 //!<! Eta vs Phi distance between track and cluster, E cluster > 0.5 GeV.
  TH2F *   fhTrackMatchedDEtaMC[fgkNmcTypes][2]     ;  //!<! Eta distance between track and cluster vs cluster E for mcTypes particle.
  TH2F *   fhTrackMatchedDPhiMC[fgkNmcTypes][2]     ;  //!<! Phi distance between track and cluster vs cluster E for mcTypes particle.
  TH2F *   fhTrackMatchedDEtaDPhiMC[fgkNmcTypes][2] ;  //!<! Eta vs Phi distance between track and cluster, E cluster > 0.5 GeV for mcTypes particle.
  TH2F *   fhdEdx[2]  ;                                //!<! matched track dEdx vs cluster E.
  TH2F *   fhEOverP[2];                                //!<! matched track E cluster over P track vs cluster E, after dEdx cut.
  TH2F *   fhTrackMatchedMCParticle[2];                //!<! Trace origin of matched particle.

  // Shower Shape histograms
  TH2F *   fhPtLambda0[2];                             //!<! Shower shape of (non) isolated photons (do not apply SS cut previously).
  TH2F *   fhPtLambda0TRD[2];                          //!<! Shower shape of (non) isolated photons, SM behind TRD (do not apply SS cut previously).
  TH3F *   fhPtLambda0Cent[2];                         //!<! Shower shape of (non) isolated photons (do not apply SS cut previously) vs centrality.

  // Selection parameters per supermodule number
  TH2F *   fhPtPerSM[2];                               //!<! Input particle pT distribution per SM
  TH2F *   fhPtLambda0PerSM[2][20];                    //!<! Shower shape of (non) isolated photons per supermodule (do not apply shower shape cut previously).
  TH2F *   fhPtLambda0PerSMNCellCut[2][20];            //!<! Shower shape of (non) isolated photons per supermodule (do not apply shower shape cut previously). N cell with weight > 4
  TH2F *   fhPtNCellPerSM       [2][20];               //!<! N cells with weight in cluster per cluster pT, per SM
  TH2F *   fhPtNCellLowM02PerSM [2][20];               //!<! N cells with weight in cluster per cluster pT for 0.1 < M02 < 0.3, per SM
  TH2F *   fhPtNCellHighM02PerSM[2][20];               //!<! N cells with weight in cluster per cluster pT for 0.5 < M02 < 2, per SM
 
  TH2F *   fhConeSumPtPerSM[20] ;                      //!<! Cluster and tracks Sum Pt in the cone, per supermodule.
  TH2F *   fhConeSumPtClusterPerSM[20] ;               //!<! Clusters Sum Pt in the cone, per supermodule.
  TH2F *   fhConeSumPtTrackPerSM[20] ;                 //!<! Tracks Sum Pt in the cone, per supermodule.
  
  TH2F *   fhPtInConePerSM[20] ;                       //!<! Cluster and tracks Pt in the cone, per supermodule.
  TH2F *   fhPtClusterInConePerSM[20] ;                //!<! Clusters Pt in the cone, per supermodule.
  TH2F *   fhPtTrackInConePerSM[20] ;                  //!<! Tracks Pt in the cone, per supermodule.

  // Selection parameters per T-Card index
  TH2F *   fhPtPerTCardIndex[2];                       //!<! Input particle pT distribution per T-Card index.
  TH2F *   fhPtLambda0PerTCardIndex[2][16];            //!<! Shower shape of (non) isolated photons per T-Card index (do not apply shower shape cut previously).
  
  TH2F *   fhConeSumPtPerTCardIndex[16] ;              //!<! Cluster and tracks Sum Pt in the cone, per T-Card index.
  TH2F *   fhConeSumPtClusterPerTCardIndex[16] ;       //!<! Clusters Sum Pt in the cone, per T-Card index.
  TH2F *   fhConeSumPtTrackPerTCardIndex[16] ;         //!<! Tracks Sum Pt in the cone, per T-Card index.
  
  TH2F *   fhPtInConePerTCardIndex[16] ;               //!<! Cluster and tracks Pt in the cone, per T-Card index.
  TH2F *   fhPtClusterInConePerTCardIndex[16] ;        //!<! Clusters Pt in the cone, per T-Card index.
  TH2F *   fhPtTrackInConePerTCardIndex[16] ;          //!<! Tracks Pt in the cone, per T-Card index.
  
  // Local maxima
  TH2F *   fhPtLambda0LocMax1[2] ;                     //!<! Pt vs lambda0 of selected cluster, 1 local maxima in cluster.
  TH2F *   fhPtLambda1LocMax1[2] ;                     //!<! Pt vs lambda1 of selected cluster, 1 local maxima in cluster.
  TH2F *   fhPtLambda0LocMax2[2] ;                     //!<! Pt vs lambda0 of selected cluster, 2 local maxima in cluster.
  TH2F *   fhPtLambda1LocMax2[2] ;                     //!<! Pt vs lambda1 of selected cluster, 2 local maxima in cluster.
  TH2F *   fhPtLambda0LocMaxN[2] ;                     //!<! Pt vs lambda0 of selected cluster, N>2 local maxima in cluster.
  TH2F *   fhPtLambda1LocMaxN[2] ;                     //!<! Pt vs lambda1 of selected cluster, N>2 local maxima in cluster.
  
  // Pile-up
  TH1F *   fhPtPileUp[7][2] ;                          //!<! Number of isolated particles.
  
  TH2F *   fhTimeENoCut;                               //!<! Time of cluster vs E, no cut.
  TH2F *   fhTimeESPD;                                 //!<! Time of cluster vs E, IsSPDPileUp.
  TH2F *   fhTimeESPDMulti;                            //!<! Time of cluster vs E, IsSPDPileUpMulti.
  TH2F *   fhTimeNPileUpVertSPD;                       //!<! Time of cluster vs n pile-up vertices from SPD.
  TH2F *   fhTimeNPileUpVertTrack;                     //!<! Time of cluster vs n pile-up vertices from Tracks.
  TH2F *   fhTimeNPileUpVertContributors;              //!<! Time of cluster vs n pile-up vertex from SPD contributors.
  TH2F *   fhTimePileUpMainVertexZDistance;            //!<! Time of cluster vs difference of z main vertex and pile-up vertex.
  TH2F *   fhTimePileUpMainVertexZDiamond;             //!<! Time of cluster vs difference of z diamond and pile-up vertex.
  
  TH2F *   fhMCConversionVertex[2];                    //!<! Conversion distance for photon clusters that have at least a contributor from the conversion. Iso and not iso
  TH2F *   fhMCConversionVertexTRD[2];                 //!<! Conversion distance for photon clusters that have at least a contributor from the conversion. Iso and not iso, SM covered by TRD
  TH2F *   fhMCConversionLambda0Rcut[6][2];            //!<! Shower shape of photon conversions, depending on conversion vertex.
  TH2F *   fhMCConversionLambda0RcutTRD[6][2];         //!<! Shower shape of photon conversions, depending on conversion vertex. SM covered by TRD
  
//TH2F *   fhLam0EMCALRegion   [2][4][3];                //!<! Cluster lambda0 vs  E, in different EMCal regions
//TH2F *   fhLam0EMCALRegionTRD[2][4][3];                //!<! Cluster lambda0 vs  E, in different EMCal regions, SM covered by TRD
//TH2F *   fhLam0EMCALRegionMCConvRcut   [2][4][3][6];   //!<! Cluster lambda0 vs  E, in different EMCal regions, MC photon conversions, depending on conversion vertex
//TH2F *   fhLam0EMCALRegionTRDMCConvRcut[2][4][3][6];   //!<! Cluster lambda0 vs  E, in different EMCal regions, SM covered by TRD,  MC photon conversions, depending on conversion vertex

  TH2F *   fhLam0EMCALRegionPerSM         [2][4][3][20]; //!<! Cluster lambda0 vs  E, in different EMCal regions
  TH2F *   fhConeSumPtTrackEMCALRegionPerSM  [4][3][20]; //!<! Track pT sum in cone vs  trigger pT, in different EMCal regions
  TH2F *   fhConeSumPtClusterEMCALRegionPerSM[4][3][20]; //!<! Cluster pT sum in cone vs  trigger pT, in different EMCal regions
  TH2F *   fhEtaPhiLam0BinPtBin[2][7];                   //!<! Cluster eta/phi for a given l0 bin (0.3-0.4) and different E bins 2-3,3-4,4-5,5-6,6-8,8-10,10-12

  TH2F *   fhPtClusterInConePerRCut;                     //!<! Clusters Pt in the cone for different cone sizes, x axis.
  TH2F *   fhPtTrackInConePerRCut;                       //!<! Tracks Pt in the cone for different cone sizes, x axis.
  TH2F *   fhConeSumPtClusterPerRCut;                    //!<! Clusters Sum Pt in the cone for different cone sizes, x axis.
  TH2F *   fhConeSumPtTrackPerRCut;                      //!<! Tracks Sum Pt in the cone for different cone sizes, x axis.

  // Variation of min pt in cone constituents
  //
  TH2F *   fhConeNClusterPerMinPtCut;                    //!<! N Clusters in the cone for different min pT cuts, x axis.
  TH2F *   fhEtaBandConeNClusterPerMinPtCut;             //!<! N Clusters in the eta band for different min pT cuts, x axis.
  TH2F *   fhConeNClusterSubEtaBandPerMinPtCut;          //!<! N Clusters in cone minus in the eta band for different min pT cuts, x axis.
  TH2F *   fhPhiBandConeNClusterPerMinPtCut;             //!<! N Clusters in the phi band for different min pT cuts, x axis.
  TH2F *   fhConeNClusterSubPhiBandPerMinPtCut;          //!<! N Clusters in cone minus in the phi band for different min pT cuts, x axis.
  
  TH2F *   fhConeNTrackPerMinPtCut;                      //!<! N Tracks in the cone for different min pT cuts, x axis.
  TH2F *   fhPerpConeNTrackPerMinPtCut;                  //!<! N Tracks in the perpendicular cone for different min pT cuts, x axis.
  TH2F *   fhConeNTrackSubPerpPerMinPtCut;               //!<! N Tracks in cone minus in the perpendicular cone for different min pT cuts, x axis.
  TH2F *   fhEtaBandConeNTrackPerMinPtCut;               //!<! N Tracks in the eta band for different min pT cuts, x axis.
  TH2F *   fhConeNTrackSubEtaBandPerMinPtCut;            //!<! N Tracks in cone minus in the eta band for different min pT cuts, x axis.
  TH2F *   fhPhiBandConeNTrackPerMinPtCut;               //!<! N Tracks in the phi band for different min pT cuts, x axis.
  TH2F *   fhConeNTrackSubPhiBandPerMinPtCut;            //!<! N Tracks in cone minus in the phi band for different min pT cuts, x axis.
  
  TH3F *   fhConeNClusterPerMinPtCutCent;                //!<! N Clusters in the cone for different min pT cuts, x axis vs centrality.
  TH3F *   fhEtaBandConeNClusterPerMinPtCutCent;         //!<! N Clusters in the eta band for different min pT cuts, x axis vs centrality.
  TH3F *   fhConeNClusterSubEtaBandPerMinPtCutCent;      //!<! N Clusters in cone minus in the eta band for different min pT cuts, x axis vs centrality.
  TH3F *   fhPhiBandConeNClusterPerMinPtCutCent;         //!<! N Clusters in the phi band for different min pT cuts, x axis vs centrality.
  TH3F *   fhConeNClusterSubPhiBandPerMinPtCutCent;      //!<! N Clusters in cone minus in the phi band for different min pT cuts, x axis vs centrality.

  TH3F *   fhConeNTrackPerMinPtCutCent;                  //!<! N Tracks in the cone for different min pT cuts, x axis vs centrality.
  TH3F *   fhPerpConeNTrackPerMinPtCutCent;              //!<! N Tracks in the perpendicular cone for different min pT cuts, x axis vs centrality.
  TH3F *   fhConeNTrackSubPerpPerMinPtCutCent;           //!<! N Tracks in cone minux the perpendicular cone for different min pT cuts, x axis vs centrality.
  TH3F *   fhEtaBandConeNTrackPerMinPtCutCent;           //!<! N Tracks in the eta band for different min pT cuts, x axisvs centrality.
  TH3F *   fhConeNTrackSubEtaBandPerMinPtCutCent;        //!<! N Tracks in cone minus in the eta band for different min pT cuts, x axis vs centrality.
  TH3F *   fhPhiBandConeNTrackPerMinPtCutCent;           //!<! N Tracks in the phi band for different min pT cuts, x axis vs centrality.
  TH3F *   fhConeNTrackSubPhiBandPerMinPtCutCent;        //!<! N Tracks in cone minus in the phi band for different min pT cuts, x axis vs centrality.
  
  TH2F *   fhConeSumPtClusterPerMinPtCut;                //!<! Clusters Sum Pt in the cone for different min pT cuts, x axis.
  TH2F *   fhEtaBandConeSumPtClusterPerMinPtCut;         //!<! Clusters Sum Pt in the eta band for different min pT cuts, x axis.
  TH2F *   fhConeSumPtClusterSubEtaBandPerMinPtCut;      //!<! Clusters Sum Pt in cone minus the eta band for different min pT cuts, x axis.
  TH2F *   fhPhiBandConeSumPtClusterPerMinPtCut;         //!<! Clusters Sum Pt in the phi band for different min pT cuts, x axis.
  TH2F *   fhConeSumPtClusterSubPhiBandPerMinPtCut;      //!<! Clusters Sum Pt in cone minus the phi band for different min pT cuts, x axis.
  
  TH2F *   fhConeSumPtTrackPerMinPtCut;                  //!<! Tracks Sum Pt in the cone for different min pT cuts, x axis.
  TH2F *   fhPerpConeSumPtTrackPerMinPtCut;              //!<! Tracks Sum Pt in the perpendicular cone for different min pT cuts, x axis.
  TH2F *   fhConeSumPtTrackSubPerpPerMinPtCut;           //!<! Tracks Sum Pt in cone minus the perpendicular cone for different min pT cuts, x axis.
  TH2F *   fhEtaBandConeSumPtTrackPerMinPtCut;           //!<! Tracks Sum Pt in the eta band for different min pT cuts, x axis.
  TH2F *   fhConeSumPtTrackSubEtaBandPerMinPtCut;        //!<! Tracks Sum Pt in cone minus the eta band for different min pT cuts, x axis.
  TH2F *   fhPhiBandConeSumPtTrackPerMinPtCut;           //!<! Tracks Sum Pt in the phi band for different min pT cuts, x axis.
  TH2F *   fhConeSumPtTrackSubPhiBandPerMinPtCut;        //!<! Tracks Sum Pt in cone minus the phi band for different min pT cuts, x axis.
  
  TH3F *   fhConeSumPtClusterPerMinPtCutCent;            //!<! Clusters Sum Pt in the cone for different min pT cuts, x axis vs centrality.
  TH3F *   fhEtaBandConeSumPtClusterPerMinPtCutCent;     //!<! Clusters Sum Pt in the eta band for different min pT cuts, x axis vs centrality.
  TH3F *   fhConeSumPtClusterSubEtaBandPerMinPtCutCent;  //!<! Clusters Sum Pt in cone minus the eta band for different min pT cuts, x axis vs centrality.
  TH3F *   fhPhiBandConeSumPtClusterPerMinPtCutCent;     //!<! Clusters Sum Pt in the phi band for different min pT cuts, x axis vs centrality.
  TH3F *   fhConeSumPtClusterSubPhiBandPerMinPtCutCent;  //!<! Clusters Sum Pt in cone minus the phi band for different min pT cuts, x axis vs centrality.
  
  TH3F *   fhConeSumPtTrackPerMinPtCutCent;              //!<! Tracks Sum Pt in the cone for different min pT cuts, x axis vs centrality.
  TH3F *   fhPerpConeSumPtTrackPerMinPtCutCent;          //!<! Tracks Sum Pt in the perpendicular cone for different min pT cuts, x axis vs centrality.
  TH3F *   fhConeSumPtTrackSubPerpPerMinPtCutCent;       //!<! Tracks Sum Pt in cone minus the perpendicular cone for different min pT cuts, x axis vs centrality.
  TH3F *   fhEtaBandConeSumPtTrackPerMinPtCutCent;       //!<! Tracks Sum Pt in the eta band for different min pT cuts, x axis vs centrality.
  TH3F *   fhConeSumPtTrackSubEtaBandPerMinPtCutCent;    //!<! Tracks Sum Pt in cone minus the eta band for different min pT cuts, x axis vs centrality.
  TH3F *   fhPhiBandConeSumPtTrackPerMinPtCutCent;       //!<! Tracks Sum Pt in the phi band for different min pT cuts, x axis vs centrality.
  TH3F *   fhConeSumPtTrackSubPhiBandPerMinPtCutCent;    //!<! Tracks Sum Pt in cone minus the phi band for different min pT cuts, x axis vs centrality.
  
  TH2F *   fhConeNSubEtaBandPerMinPtCut ;                //!<! N Clusters+Tracks in cone minus in the eta band for different min pT cuts, x axis.
  TH2F *   fhConeNSubPhiBandPerMinPtCut ;                //!<! N Clusters+Tracks in cone minus in the eta band for different min pT cuts, x axis.
  TH2F *   fhConeSumPtSubEtaBandPerMinPtCut;             //!<! Clusters+Tracks Sum Pt in cone minus the eta band for different min pT cuts, x axis.
  TH2F *   fhConeSumPtSubPhiBandPerMinPtCut;             //!<! Clusters+Tracks Sum Pt in cone minus the phi band for different min pT cuts, x axis.

  TH3F *   fhConeNSubEtaBandPerMinPtCutCent ;            //!<! N Clusters+Tracks in cone minus in the eta band for different min pT cuts, x axis vs centrality.
  TH3F *   fhConeNSubPhiBandPerMinPtCutCent ;            //!<! N Clusters+Tracks in cone minus in the eta band for different min pT cuts, x axis vs centrality.
  TH3F *   fhConeSumPtSubEtaBandPerMinPtCutCent;         //!<! Clusters+Tracks Sum Pt in cone minus the eta band for different min pT cuts, x axis vs centrality.
  TH3F *   fhConeSumPtSubPhiBandPerMinPtCutCent;         //!<! Clusters+Tracks Sum Pt in cone minus the phi band for different min pT cuts, x axis vs centrality.
  
  TH2F *   fhConeSumPtClusterPerMaxPtCut;                //!<! Clusters Sum Pt in the cone for different max pT cuts, x axis.
  TH2F *   fhConeSumPtTrackPerMaxPtCut;                  //!<! Tracks Sum Pt in the cone for different max pT cuts, x axis.
  TH2F *   fhConeSumPtTrackPerEtaCut;                    //!<! Tracks Sum Pt in the cone for different min eta cuts, x axis.
  
  TH2F *   fhPtClusterInConePerNCellCut;                 //!<! Clusters Pt in the cone for different min cluster n cell cut, x axis.
  TH2F *   fhPtTrackInConePerNCellCut;                   //!<! Tracks Pt in the cone for different min cluster n cell cut, x axis.

  TH2F *   fhConeSumPtClusterPerNCellCut;                //!<! Clusters Sum Pt in the cone for different min cluster n cell cut, x axis.
  TH2F *   fhConeSumPtTrackPerNCellCut;                  //!<! Tracks Sum Pt in the cone for different min cluster n cell cut, x axis.
  
  TH3F *   fhPtClusterInConePerNCellPerSM [4];           //!<! Clusters Pt in the cone for different min cluster n cell cut, x axis vs SM number, 8<E<12 GeV, 3 shower bins
  TH3F *   fhPtTrackInConePerNCellPerSM   [4];           //!<! Tracks Pt in the cone for different min cluster n cell cut, x axis, vs SM number, 8<E<12 GeV, 3 shower bins
  TH3F *   fhConeSumPtClusterPerNCellPerSM[4];           //!<! Clusters Sum Pt in the cone for different min cluster n cell cut, x axis, vs SM number, 8<E<12 GeV, 3 shower bins
  TH3F *   fhConeSumPtTrackPerNCellPerSM  [4];           //!<! Tracks Sum Pt in the cone for different min cluster n cell cut, x axis, vs SM number, 8<E<12 GeV, 3 shower bins

  TH2F *   fhPtClusterInConePerExoCut;                   //!<! Clusters Pt in the cone for different exoticity cut, x axis.
  TH2F *   fhPtTrackInConePerExoCut;                     //!<! Tracks Pt in the cone for different exoticity cut, x axis.

  TH2F *   fhConeSumPtClusterPerExoCut;                  //!<! Clusters Sum Pt in the cone for different exoticity cut, x axis.
  TH2F *   fhConeSumPtTrackPerExoCut;                    //!<! Tracks Sum Pt in the cone for different exoticity cut, x axis.

  TH2F *   fhConeSumPtTrackTOFBC0;                       //!<! track with TOF hit sum pt, tof in BC0 
  TH2F *   fhConeSumPtTrackTOFBCN;                       //!<! track with TOF hit sum pt, tof not in BC0 
  TH2F *   fhConeSumPtTrackTOFNo ;                       //!<! track without TOF hit sum pt 
  TH2F *   fhPtTrackInConeTOFBC0;                        //!<! track with TOF hit, pt, tof in BC0 
  TH2F *   fhPtTrackInConeTOFBCN;                        //!<! track with TOF hit, pt, tof not in BC0 
  TH2F *   fhPtTrackInConeTOFNo ;                        //!<! track without TOF hit, pt 

  TH2F *   fhPhiTrackInCone;                             //!<! track azhimuthal angle
  TH2F *   fhEtaTrackInCone;                             //!<! track pseudo-rapidity
  TH2F *   fhEtaPhiTrackInCone;                          //!<! track azhimuthal angle vs pseudo-rapidity
  
  TH2F *   fhPhiTrackInConeTOFBC0;                       //!<! track with TOF hit, phi, tof in BC0 
  TH2F *   fhPhiTrackInConeTOFBCN;                       //!<! track with TOF hit, phi, tof not in BC0 
  TH2F *   fhPhiTrackInConeTOFNo ;                       //!<! track without TOF hit, phi 
  TH2F *   fhEtaTrackInConeTOFBC0;                       //!<! track with TOF hit, eta, tof in BC0 
  TH2F *   fhEtaTrackInConeTOFBCN;                       //!<! track with TOF hit, eta, tof not in BC0 
  TH2F *   fhEtaTrackInConeTOFNo ;                       //!<! track without TOF hit, eta 
  TH2F *   fhEtaPhiTrackInConeTOFBC0;                    //!<! track with TOF hit, eta-phi, tof in BC0 
  TH2F *   fhEtaPhiTrackInConeTOFBCN;                    //!<! track with TOF hit, eta-phi, tof not in BC0 
  TH2F *   fhEtaPhiTrackInConeTOFNo ;                    //!<! track without TOF hit, eta-phi 
  TH2F *   fhTrackTOFInCone ;                            //!<! track TOF in cone
  TH2F *   fhTrackTOFInConeBC0 ;                         //!<! track TOF in cone and BC0
  TH2F *   fhTrackTOFInConeExoTrigger ;                  //!<! track TOF in cone, trigger is exotic

  TH2F *   fhConeSumPtTrackITSRefitOnSPDOn;              //!<! track with ITS Refit On SPD On 
  TH2F *   fhConeSumPtTrackITSRefitOnSPDOff;             //!<! track with ITS Refit On SPD Off 
  TH2F *   fhConeSumPtTrackITSRefitOffSPDOff ;           //!<! track with ITS Refit Off SPD Off 
  TH2F *   fhPtTrackInConeITSRefitOnSPDOn;               //!<! track with ITS Refit On SPD On
  TH2F *   fhPtTrackInConeITSRefitOnSPDOff;              //!<! track with ITS Refit On SPD Off 
  TH2F *   fhPtTrackInConeITSRefitOffSPDOff ;            //!<! track with ITS Refit Off SPD Off
  TH2F *   fhPhiTrackInConeITSRefitOnSPDOn;              //!<! track with ITS Refit On SPD On 
  TH2F *   fhPhiTrackInConeITSRefitOnSPDOff;             //!<! track with ITS Refit On SPD Off
  TH2F *   fhPhiTrackInConeITSRefitOffSPDOff ;           //!<! track with ITS Refit Off SPD Off 
  TH2F *   fhEtaTrackInConeITSRefitOnSPDOn;              //!<! track with ITS Refit On SPD On
  TH2F *   fhEtaTrackInConeITSRefitOnSPDOff;             //!<! track with ITS Refit On SPD Off 
  TH2F *   fhEtaTrackInConeITSRefitOffSPDOff ;           //!<! track with ITS Refit Off SPD Off 
  TH2F *   fhEtaPhiTrackInConeITSRefitOnSPDOn;           //!<! track with ITS Refit On SPD On
  TH2F *   fhEtaPhiTrackInConeITSRefitOnSPDOff;          //!<! track with ITS Refit On SPD Off
  TH2F *   fhEtaPhiTrackInConeITSRefitOffSPDOff ;        //!<! track with ITS Refit Off SPD Off

  TH2F *   fhConeSumPtTrackTOFBC0ITSRefitOnSPDOn;        //!<! track with ITS Refit On SPD On, TOF BC=0 
  TH2F *   fhPtTrackInConeTOFBC0ITSRefitOnSPDOn;         //!<! track with ITS Refit On SPD On, TOF BC=0
  TH2F *   fhPhiTrackInConeTOFBC0ITSRefitOnSPDOn;        //!<! track with ITS Refit On SPD On, TOF BC=0 
  TH2F *   fhEtaTrackInConeTOFBC0ITSRefitOnSPDOn;        //!<! track with ITS Refit On SPD On, TOF BC=0 
  TH2F *   fhEtaPhiTrackInConeTOFBC0ITSRefitOnSPDOn;     //!<! track with ITS Refit On SPD On, TOF BC=0
  
  TH2F *   fhPerpConeSumPtITSRefitOnSPDOn ;              //!<! Sum track Pt in cone at the perpendicular phi region to trigger axis  (phi +90), ITS Refit On, SPD On
  TH2F *   fhPtInPerpConeITSRefitOnSPDOn ;               //!<! track Pt  in cone at the perpendicular phi region to trigger axis  (phi +90), ITS Refit On, SPD On
  TH2F *   fhEtaPhiInPerpConeITSRefitOnSPDOn ;           //!<! tracl eta vs phi in cone at the perpendicular phi region to trigger axis  (phi +90), ITS Refit On, SPD On

  TH2F *   fhPerpConeSumPtTOFBC0ITSRefitOnSPDOn ;        //!<! Sum track Pt in cone at the perpendicular phi region to trigger axis  (phi +90), ITS Refit On, SPD On, TOF BC=0
  TH2F *   fhPtInPerpConeTOFBC0ITSRefitOnSPDOn ;         //!<! track Pt  in cone at the perpendicular phi region to trigger axis  (phi +90), ITS Refit On, SPD On, TOF BC=0
  TH2F *   fhEtaPhiInPerpConeTOFBC0ITSRefitOnSPDOn ;     //!<! tracl eta vs phi in cone at the perpendicular phi region to trigger axis  (phi +90), ITS Refit On, SPD On, TOF BC=0
  
  TH2F *   fhPtTrackInConeDCA[3];                        //!<! track DCAxy,z,constrained vs track pT, in cone with trigger pT > 10 GeV
  TH2F *   fhPtTrackInPerpConeDCA[3];                    //!<! track DCAxy,z,constrained vs track pT, in perpendicular cone trigger pT > 10 GeV
   
  /// Copy constructor not implemented.
  AliAnaParticleIsolation(              const AliAnaParticleIsolation & iso) ;
    
  /// Assignment operator not implemented.
  AliAnaParticleIsolation & operator = (const AliAnaParticleIsolation & iso) ;
  
  /// \cond CLASSIMP
  ClassDef(AliAnaParticleIsolation,51) ;
  /// \endcond

} ;


#endif //ALIANAPARTICLEISOLATION_H



