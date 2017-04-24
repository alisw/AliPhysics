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
/// This class takes a particle AOD object with format AliAODPWG4Particle
/// produced by any of the identified particle classes (AliAnaPhoton, AliAnaElectron,
/// AliAnaPi0EbE, AliAnaChargedParticle) and checks if there is low particle environment
/// around it with the utils of AliIsolationCut, declaring the particle AOD object as isolated or not.
///
///  Class created from old AliPHOSGammaJet (see AliRoot versions previous Release 4-09).
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
class AliAODPWG4Particle;
class AliAODPWG4ParticleCorrelation ;

class AliAnaParticleIsolation : public AliAnaCaloTrackCorrBaseClass {

 public:
    
  AliAnaParticleIsolation() ;
    
  /// Virtual destructor.
  virtual ~AliAnaParticleIsolation() { ; }

  // Main general methods
    
  void         CalculateCaloUEBand    (AliAODPWG4ParticleCorrelation * pCandidate,
                                       Float_t & etaBand, Float_t & phiBand) ;
    
  void         CalculateCaloCellUEBand(AliAODPWG4ParticleCorrelation * pCandidate,
                                       Float_t & etaBand, Float_t & phiBand) ;
    
  void         CalculateTrackUEBand   (AliAODPWG4ParticleCorrelation * pCandidate,
                                       Float_t & etaBand, Float_t & phiBand) ;
  
  void         CalculateCaloSignalInCone    (AliAODPWG4ParticleCorrelation * aodParticle, Float_t & coneptsumCluster, Float_t & coneptLeadCluster) ;
    
  void         CalculateCaloCellSignalInCone(AliAODPWG4ParticleCorrelation * aodParticle, Float_t & coneptsumCell) ;

  void         CalculateTrackSignalInCone   (AliAODPWG4ParticleCorrelation * aodParticle, Float_t & coneptsumTrack  , Float_t & coneptLeadTrack  ) ;


  void         CalculateNormalizeUEBandPerUnitArea(AliAODPWG4ParticleCorrelation * pCandidate, Float_t coneptsumCluster,
                                                   Float_t coneptsumCell,  Float_t coneptsumTrack,
                                                   Float_t &coneptsumSubEtaBand, Float_t &coneptsumSubPhiBand ) ;
  
  TObjString * GetAnalysisCuts() ;
  
  TList      * GetCreateOutputObjects() ;
  
  void         Init() ;

  void         InitParameters() ;
  
  void         MakeAnalysisFillAOD()  ;
  
  void         MakeAnalysisFillHistograms() ;
  
  void         Print( const Option_t * opt ) const ;
 
  // Analysis specific methods
  
  void         FillPileUpHistograms(Float_t energy, Float_t time) ; //Int_t clusterID) ;
  
  void         FillAcceptanceHistograms();
 
  void         FillTrackMatchingShowerShapeControlHistograms(AliAODPWG4ParticleCorrelation  * pCandidate,
                                                             Float_t coneptsum,  Float_t coneptsumTrack, Float_t coneptsumCluster, 
                                                             Float_t coneleadpt, Int_t mcIndex) ;
  
  Bool_t       IsTriggerTheNearSideEventLeadingParticle(Int_t & idLeading);
  
  void         MakeSeveralICAnalysis( AliAODPWG4ParticleCorrelation * ph, Int_t mcIndex ) ;
  
  void         StudyEMCALRegions(Float_t pt, Float_t phi, Float_t eta, Float_t m02, 
                                 Float_t coneptsumTrack, Float_t coneptsumCluster, 
                                 Bool_t isolated, Int_t iSM) ;
  
  void         StudyMCConversionRadius(Float_t  pt, Bool_t isolated, Int_t iSM, 
                                       Float_t m02, Int_t     mcTag, Int_t label) ;
  
  // Analysis Setters and Getters
  
  TString      GetTriggerDetectorString()      const { return fIsoDetectorString ; }
  TString      GetTriggerDetector()            const { return fIsoDetector       ; }
  Int_t        GetNCones()                     const { return fNCones            ; }
  Int_t        GetNPtThresFrac()               const { return fNPtThresFrac      ; }
  Float_t      GetConeSizes(Int_t i)           const { return fConeSizes[i]      ; }
  Float_t      GetPtThresholds(Int_t i)        const { return fPtThresholds[i]   ; }
  Float_t      GetSumPtThresholds(Int_t i)     const { return fSumPtThresholds[i]; }
  Float_t      GetPtFractions(Int_t i)         const { return fPtFractions[i]    ; }
  
  Int_t        GetMCIndex(Int_t mcTag);
  
  void         SetTriggerDetector(TString & det)     ;
  void         SetTriggerDetector(Int_t  det)        ;
  void         SetNCones(Int_t ncs)                  { fNCones          = ncs    ; }
  void         SetNPtThresFrac(Int_t npt)            { fNPtThresFrac    = npt    ; }
  void         SetConeSizes(Int_t i, Float_t r)      { fConeSizes[i]    = r      ; }
  void         SetPtThresholds(Int_t i, Float_t pt)  { fPtThresholds[i] = pt     ; }
  void         SetPtFractions(Int_t i, Float_t pt)   { fPtFractions[i]  = pt     ; } 
  void 	       SetSumPtThresholds(Int_t i, Float_t pt){ fSumPtThresholds[i] = pt ; }

  void         SetMinCellsAngleOverlap(Float_t n)    { fMinCellsAngleOverlap = n ; }
  
  Bool_t       IsReIsolationOn()               const { return fReMakeIC          ; }
  void         SwitchOnReIsolation()                 { fReMakeIC      = kTRUE    ; }
  void         SwitchOffReIsolation()                { fReMakeIC      = kFALSE   ; }
  
  Bool_t       IsSeveralIsolationOn()          const { return fMakeSeveralIC     ; }
  void         SwitchOnSeveralIsolation()            { fMakeSeveralIC = kTRUE    ; }
  void         SwitchOffSeveralIsolation()           { fMakeSeveralIC = kFALSE   ; }
  
  void         SwitchOnTMHistoFill()                 { fFillTMHisto   = kTRUE    ; }
  void         SwitchOffTMHistoFill()                { fFillTMHisto   = kFALSE   ; }
  
  void         SwitchOnSSHistoFill()                 { fFillSSHisto   = kTRUE    ; }
  void         SwitchOffSSHistoFill()                { fFillSSHisto   = kFALSE   ; }

  void         SwitchOnFillEMCALRegionHistograms()   { fFillEMCALRegionHistograms = kTRUE  ; }
  void         SwitchOffFillEMCALRegionHistograms()  { fFillEMCALRegionHistograms = kFALSE ; }  
  
  Bool_t       IsLeadingOnlyOn()               const { return fLeadingOnly       ; }
  void         SwitchOnLeadingOnly()                 { fLeadingOnly    = kTRUE   ; }
  void         SwitchOffLeadingOnly()                { fLeadingOnly    = kFALSE  ; }
  
  void         SwitchOnCheckNeutralClustersForLeading() { fCheckLeadingWithNeutralClusters = kTRUE  ; }
  void         SwitchOffCheckNeutralClustersForLeading(){ fCheckLeadingWithNeutralClusters = kFALSE ; }
  
  void         SwitchOnUEBandSubtractionHistoFill(Int_t lev = 1)  { fFillUEBandSubtractHistograms = lev ; }
  void         SwitchOffUEBandSubtractionHistoFill() { fFillUEBandSubtractHistograms = 0   ; }

  void         SwitchOnCellHistoFill()               { fFillCellHistograms = kTRUE ; }
  void         SwitchOffCellHistoFill()              { fFillCellHistograms = kFALSE; }

  void         SwitchOnNLMHistoFill()                { fFillNLMHistograms = kTRUE ; }
  void         SwitchOffNLMHistoFill()               { fFillNLMHistograms = kFALSE; }
  
  void         SwitchOnDecayTaggedHistoFill()        { fFillTaggedDecayHistograms = kTRUE ; }
  void         SwitchOffDecayTaggedHistoFill()       { fFillTaggedDecayHistograms = kFALSE; }
  void         SetNDecayBits(Int_t n)                { fNDecayBits = n               ; }
  void         SetDecayBits(Int_t i, UInt_t bit)     { if(i < AliNeutralMesonSelection::fgkMaxNDecayBits)
                                                       fDecayBits[i] = bit           ; }
  void         SetM02CutForTaggedDecays(Float_t m02) { fDecayTagsM02Cut        = m02 ; }
  
  void         SwitchOnBackgroundBinHistoFill()      { fFillBackgroundBinHistograms = kTRUE ; }
  void         SwitchOffBackgroundBinHistoFill()     { fFillBackgroundBinHistograms = kFALSE; }
  void         SetNBackgroundBins(Int_t n)           { if(n < 19) fNBkgBin = n ; }
  void         SetBackgroundLimits(Int_t i,Float_t l){ if(i <= fNBkgBin) fBkgBinLimit[i] = l; }

  void         SwitchOnPtTrigBinHistoFill()          { fFillPtTrigBinHistograms = kTRUE ; }
  void         SwitchOffPtTrigBinHistoFill()         { fFillPtTrigBinHistograms = kFALSE; }
  void         SetNPtTrigBins(Int_t n)               { if(n < 19) fNPtTrigBin = n ; }
  void         SetPtTrigLimits(Int_t i,Float_t l)    { if(i <= fNPtTrigBin) fPtTrigBinLimit[i] = l; }

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
  
  /// For primary histograms in arrays, index in the array, corresponding to a photon origin.
  enum mcPrimTypes { kmcPrimPhoton = 0, kmcPrimPi0Decay = 1, kmcPrimEtaDecay  = 2, kmcPrimOtherDecay  = 3,
                     kmcPrimPrompt = 4, kmcPrimFrag     = 5, kmcPrimISR       = 6,
                     kmcPrimPi0    = 7, kmcPrimEta      = 8                                               } ;
  
  static const Int_t fgkNmcPrimTypes = 9; ///< Number of MC primary particle types used in the analysis in the histogram arrays.
  
  /// For histograms in arrays, index in the array, corresponding to any particle origin.
  enum mcTypes     { kmcPhoton   = 0, kmcPrompt     = 1, kmcFragment         = 2,
                     kmcPi0      = 3, kmcPi0Decay   = 4, kmcPi0DecayLostPair = 5,
                     kmcEta      = 6, kmcEtaDecay   = 7, kmcEtaDecayLostPair = 8,
                     kmcOtherDecay=9, kmcElectron   =10, kmcHadron           =11  } ;
  
  static const Int_t fgkNmcTypes = 12; ///< Number of MC type particles originating the clusters used in the analysis in the histogram arrays.

 private:
  
  Int_t    fIsoDetector ;                             ///<  Candidate particle for isolation detector.
  TString  fIsoDetectorString ;                       ///<  Candidate particle for isolation detector.
  Bool_t   fReMakeIC ;                                ///<  Do isolation analysis.
  Bool_t   fMakeSeveralIC ;                           ///<  Do analysis for different IC.
  Bool_t   fFillTMHisto;                              ///<  Fill track matching plots.
  Bool_t   fFillSSHisto;                              ///<  Fill Shower shape plots.
  Bool_t   fFillEMCALRegionHistograms ;               ///<  Fill histograms in EMCal slices
  Int_t    fFillUEBandSubtractHistograms;             ///<  Fill histograms working on the UE subtraction. 1 fill basic histograms, 2 more detailed
  Bool_t   fFillCellHistograms;                       ///<  Fill cell histograms.
  Bool_t   fFillOverlapHistograms;                    ///<  Fill histograms that depend on number of overlaps
  Bool_t   fStudyTracksInCone;                        ///<  Study tracks depending on different track info
  Bool_t   fStudyMCConversionRadius;                  ///<  Study shower shape depending the conversion radius
  Bool_t   fStudyExoticTrigger;                       ///<  Fill histograms with track and cluster pT when the trigger is exotic
  
  Bool_t   fFillTaggedDecayHistograms;                ///<  Fill histograms for clusters tagged as decay.
  Int_t    fNDecayBits ;                              ///<  In case of study of decay triggers, select the decay bit.
  UInt_t   fDecayBits[AliNeutralMesonSelection::fgkMaxNDecayBits] ; ///< In case of study of decay triggers, select the decay. bit
  Float_t  fDecayTagsM02Cut ;                         ///<  Apply a m02 cut to clusters tagged as decay.
  
  Bool_t   fFillNLMHistograms;                        ///<  Fill NLM histograms.
  Bool_t   fLeadingOnly;                              ///<  Do isolation with leading particle.
  Bool_t   fCheckLeadingWithNeutralClusters;          ///<  Compare the trigger candidate to Leading pT with the clusters pT, by default only charged.
  Bool_t   fSelectPrimariesInCone;                    ///<  In primary particle isolation studies, select only particles in isolation cone within detector acceptance and E cut.
  Bool_t   fMakePrimaryPi0DecayStudy;                 ///<  Fill dedicated histograms for primary decay photons.
  
  Bool_t   fFillBackgroundBinHistograms;              ///<  Fill histograms for different bins in pt content of the cone.
  Int_t    fNBkgBin;                                  ///<  Number of bins on pt content in cone.
  Float_t  fBkgBinLimit[20];                          ///<  Pt bin limits on pt content in the cone.

  Bool_t   fFillPtTrigBinHistograms;                  ///<  Fill histograms for different bins in pt trigger.
  Int_t    fNPtTrigBin;                               ///<  Number of bins on pt trigger.
  Float_t  fPtTrigBinLimit[20];                       ///<  Pt bin limits on pt trigger.
  
  Float_t  fMinCellsAngleOverlap;                     ///<  Number of cells that define the cluster overlap.
  
  //  Analysis data members for multiple cones and pt thresholds
  Int_t    fNCones ;                                  ///<  Number of cone sizes to test. Multiple cones and pt thresholds analysis.
  Int_t    fNPtThresFrac ;                            ///<  Number of ptThres and ptFrac to test. Multiple cones and pt thresholds analysis.
  
  Float_t  fConeSizes[5] ;                            ///<  Array with cones to test. Multiple cones and pt thresholds analysis.
  Float_t  fPtThresholds[5] ;                         ///<  Array with pt thresholds to test. Multiple cones and pt thresholds analysis.
  Float_t  fPtFractions[5] ;                          ///<  Array with pt thresholds to test frac. Multiple cones and pt thresholds analysis.
  Float_t  fSumPtThresholds[5] ;                      ///<  Array with pt thresholds to test frac. Multiple cones and pt thresholds analysis.
  
  Bool_t   fStudyPtCutInCone;                         ///<  Activate study of track/cluster min pT on sum of pT in cone
  Int_t    fNPtCutsInCone;                            ///<  Number of track/cluster min pT cut to test in cone for sum pT calculation.
  Float_t  fMinPtCutInCone[20];                       ///<  List of track/cluster min pT cut to test in cone for sum pT calculation.
  Float_t  fMaxPtCutInCone[20];                       ///<  List of track/cluster max pT cut to test in cone for sum pT calculation.

  Bool_t   fStudyEtaCutInCone;                        ///<  Activate study of track/cluster max eta on sum of pT in cone
  Int_t    fNEtaCutsInCone;                           ///<  Number of track/cluster max eta cut to test in cone for sum pT calculation.
  Float_t  fEtaCutInCone[10];                         ///<  List of track/cluster max eta cut to test in cone for sum pT calculation.

  Bool_t   fStudyRCutInCone;                          ///<  Activate study of track/cluster sum of pT in cone with variable size
  Int_t    fNRCutsInCone;                             ///<  Number of track/cluster max R cut to test in cone for sum pT calculation.
  Float_t  fRCutInCone[10];                           ///<  List of track/cluster max R cut to test in cone for sum pT calculation.

  Int_t    fNNCellsInCandidate;                       ///<  Number of cells in cluster selection to test in cone for sum pT calculation.
  Int_t    fNCellsInCandidate[20];                    ///<  List of Number of cells in cluster selection to test in cone for sum pT calculation.

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
  Bool_t         fIsExoticTrigger;                    //!<! Trigger cluster considered as exotic
  Float_t        fClusterExoticity;                   //!<! Temporary container or currently analyzed cluster exoticity

  //Histograms  
  
  TH1F *   fhEIso ;                                    //!<! Number of isolated particles vs energy.
  TH1F *   fhPtIso ;                                   //!<! Number of isolated particles vs pT.
  TH2F *   fhPtCentralityIso ;                         //!<! Centrality vs pT.
  TH2F *   fhPtEventPlaneIso ;                         //!<! Event plane angle vs pT.
  TH2F *   fhPtNLocMaxIso ;                            //!<! Number of isolated particles vs NLM in cluster.
  TH2F *   fhPhiIso ;                                  //!<! phi of isolated particles.
  TH2F *   fhEtaIso ;                                  //!<! eta of isolated particles.
  TH2F *   fhEtaPhiIso ;                               //!<! eta vs phi of isolated particles.
  TH2F *   fhEtaPhiNoIso ;                             //!<! eta vs phi of not isolated leading particles.
  TH1F *   fhENoIso ;                                  //!<! Number of not isolated leading particles vs Energy.
  TH1F *   fhPtNoIso ;                                 //!<! Number of not isolated leading particles vs pT.
  TH2F *   fhPtNLocMaxNoIso ;                          //!<! Number of not isolated particles vs NLM in cluster.
  TH1F *   fhEIsoExoTrigger;                           //!<! Number of isolated exotic cluster vs E.
  TH1F *   fhENoIsoExoTrigger;                         //!<! Number of not isolated exotic cluster vs E.
  TH1F *   fhPtIsoExoTrigger;                          //!<! Number of isolated exotic cluster vs pT.
  TH1F *   fhPtNoIsoExoTrigger;                        //!<! Number of not isolated exotic cluster vs pT.

  
  TH1F *   fhPtDecay       [2][AliNeutralMesonSelection::fgkMaxNDecayBits]; //!<! Number of (non) isolated Pi0 decay particles (invariant mass tag).
  TH2F *   fhEtaPhiDecay   [2][AliNeutralMesonSelection::fgkMaxNDecayBits]; //!<! eta vs phi of (not) isolated leading Pi0 decay particles.
  TH2F *   fhPtLambda0Decay[2][AliNeutralMesonSelection::fgkMaxNDecayBits]; //!<! Shower shape of (non) isolated leading Pi0 decay particles (do not apply SS cut previously).

  TH2F *   fhPtInCone ;                                //!<! Cluster/track Pt in the cone.
  TH2F *   fhPtClusterInCone ;                         //!<! Cluster Pt in the cone.
  TH2F *   fhPtCellInCone ;                            //!<! Cell amplitude in the cone.
  TH2F *   fhPtTrackInCone ;                           //!<! Track Pt in the cone.
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
  TH2F *   fhPtInConeCent ;                            //!<! Particle Pt in the cone versus centrality.
  TH2F *   fhPerpConeSumPt ;                           //!<! Sum Pt in cone at the perpendicular phi region to trigger axis  (phi +90).
  TH2F *   fhPerpConeSumPtTOFBC0 ;                     //!<! Sum Pt in cone at the perpendicular phi region to trigger axis  (phi +90), TOF BC=0
  TH2F *   fhPtInPerpCone ;                            //!<! Particle Pt  in cone at the perpendicular phi region to trigger axis  (phi +90).
  TH2F *   fhPtInPerpConeTOFBC0 ;                      //!<! Particle Pt  in cone at the perpendicular phi region to trigger axis  (phi +90), TOF BC=0
  
  TH2F *   fhEtaPhiInConeCluster ;                     //!<! Eta vs. phi of clusters in cone.
  TH2F *   fhEtaPhiCluster ;                           //!<! Eta vs. phi of all clusters.
  TH2F *   fhEtaPhiInConeTrack ;                       //!<! Eta vs. phi of tracks in cone.
  TH2F *   fhEtaPhiTrack ;                             //!<! Eta vs. phi of all tracks.
  TH2F *   fhEtaPhiInPerpCone ;                        //!<! Eta vs. phi of tracks in perpendicular cone
  TH2F *   fhEtaPhiInPerpConeTOFBC0 ;                  //!<! Eta vs. phi of tracks in perpendicular cone, with TOF BC=0.
  
  TH2F *   fhEtaBandClusterEtaPhi ;                    //!<! Eta vs Phi in Eta band to estimate UE in cone, only clusters. 
  TH2F *   fhPhiBandClusterEtaPhi ;                    //!<! Eta vs Phi in Phi band to estimate UE in cone, only clusters.
  TH2F *   fhEtaBandTrackEtaPhi   ;                    //!<! Eta vs Phi in Eta band to estimate UE in cone, only tracks.
  TH2F *   fhPhiBandTrackEtaPhi   ;                    //!<! Eta vs Phi in Phi band to estimate UE in cone, only tracks. 
  TH2F *   fhEtaBandClusterPt ;                        //!<! pT in Eta band to estimate UE in cone, only clusters.
  TH2F *   fhPhiBandClusterPt ;                        //!<! pT in Phi band to estimate UE in cone, only clusters.
  TH2F *   fhEtaBandTrackPt   ;                        //!<! pT in Eta band to estimate UE in cone, only tracks.
  TH2F *   fhPhiBandTrackPt   ;                        //!<! pT in Phi band to estimate UE in cone, only tracks.
  TH2F *   fhEtaBandCell ;                             //!<! Row vs Column in Eta band to estimate UE in cone, only cells.
  TH2F *   fhPhiBandCell ;                             //!<! Row vs Column in Phi band to estimate UE in cone, only cells.

  TH2F *   fhConePtLead ;                              //!<! Cluster and tracks leading pt in the cone.
  TH2F *   fhConePtLeadCluster ;                       //!<! Clusters leading pt in the cone.
  TH2F *   fhConePtLeadTrack ;                         //!<! Tracks leading pt in the cone.
  TH2F *   fhConePtLeadClustervsTrack;                 //!<! Tracks vs Clusters leading pt.
  TH2F *   fhConePtLeadClusterTrackFrac;               //!<! Trigger pt vs cluster/track leading pt.
  
  TH2F *   fhConeSumPt ;                               //!<! Cluster and tracks Sum Pt in the cone.
//TH3F *   fhPtLambda0Eiso;                            //!<! ABCD TH3F histogram Pt, ShowerShape and sum(ET)+sum(pT) cone
  TH2F *   fhConeSumPtCellTrack ;                      //!<! Cells and tracks Sum Pt in the cone.
  TH2F *   fhConeSumPtCell ;                           //!<! Cells Sum Pt in the cone.
  TH2F *   fhConeSumPtCluster ;                        //!<! Clusters Sum Pt in the cone.
  TH2F *   fhConeSumPtTrack ;                          //!<! Tracks Sum Pt in the cone.
 
  TH2F *   fhConeSumPtExoTrigger ;                     //!<! Cluster and tracks Sum Pt in the cone. Trigger is exotic
  TH2F *   fhConeSumPtClusterExoTrigger ;              //!<! Clusters Sum Pt  in the cone. Trigger is exotic
  TH2F *   fhConeSumPtTrackExoTrigger ;                //!<! Tracks Sum Pt  in the cone. Trigger considered exotic

  TH2F *   fhConeSumPtEtaBandUECluster;                //!<! Cluster Sum Pt in the eta band for clusters, before normalization.
  TH2F *   fhConeSumPtPhiBandUECluster;                //!<! Cluster Sum Pt in the phi band for clusters, before normalization.
  TH2F *   fhConeSumPtEtaBandUETrack;                  //!<! Track Sum Pt in the eta band for tracks, before normalization.
  TH2F *   fhConeSumPtPhiBandUETrack;                  //!<! Track Sum Pt in the phi badn for tracks, before normalization.
  TH2F *   fhConeSumPtEtaBandUECell;                   //!<! Cell Sum amplitude in the eta band for cells, before normalization.
  TH2F *   fhConeSumPtPhiBandUECell;                   //!<! Cell Sum amplitude in the phi band for cells, before normalization.

  TH2F *   fhConeSumPtTrigEtaPhi ;                     //!<! Cluster and tracks Sum Pt Sum Pt in the cone, per eta-phi bin of trigger.
  TH2F *   fhConeSumPtCellTrackTrigEtaPhi ;            //!<! Cell and tracks Sum Pt Sum Pt in the cone, per eta-phi bin of trigger.
  TH2F *   fhConeSumPtEtaBandUEClusterTrigEtaPhi;      //!<! Cluster Sum Pt in the eta band for clusters, per eta-phi bin of trigger,before normalization.
  TH2F *   fhConeSumPtPhiBandUEClusterTrigEtaPhi;      //!<! Cluster Sum Pt in the phi band for clusters, per eta-phi bin of trigger, before normalization.
  TH2F *   fhConeSumPtEtaBandUETrackTrigEtaPhi;        //!<! Track Sum Pt in the eta band for tracks, per eta-phi bin of trigger, before normalization.
  TH2F *   fhConeSumPtPhiBandUETrackTrigEtaPhi;        //!<! Track Sum Pt in the phi badn for tracks, per eta-phi bin of trigger, before normalization.
  TH2F *   fhConeSumPtEtaBandUECellTrigEtaPhi;         //!<! Cluster Sum amplitude in the eta band for cells, per eta-phi bin of trigger, before normalization.
  TH2F *   fhConeSumPtPhiBandUECellTrigEtaPhi;         //!<! Cluster Sum amplitude in the phi band for cells, per eta-phi bin of trigger, before normalization.
  
  TH2F *   fhConeSumPtEtaUESub;                        //!<! Cluster and tracks Sum Pt in the cone after bkg subtraction, vs pT trigger.
  TH2F *   fhConeSumPtPhiUESub;                        //!<! Cluster and tracks Sum Pt in the cone after bkg subtraction, vs pT trigger.
  TH2F *   fhConeSumPtEtaUESubTrigEtaPhi;              //!<! Cluster and tracks Sum Pt in the cone after bkg subtraction, vs eta-phi trigger.
  TH2F *   fhConeSumPtPhiUESubTrigEtaPhi;              //!<! Cluster and tracks Sum Pt in the cone after bkg subtraction, vs eta-phi trigger.
  
  TH2F *   fhConeSumPtEtaUESubTrackCell;               //!<! Cluster and tracks Sum Pt in the cone after bkg subtraction, vs pT trigger.
  TH2F *   fhConeSumPtPhiUESubTrackCell;               //!<! Cluster and tracks Sum Pt in the cone after bkg subtraction, vs pT trigger.
  TH2F *   fhConeSumPtEtaUESubTrackCellTrigEtaPhi;     //!<! Cluster and tracks Sum Pt in the cone after bkg subtraction, vs eta-phi trigger.
  TH2F *   fhConeSumPtPhiUESubTrackCellTrigEtaPhi;     //!<! Cluster and tracks Sum Pt in the cone after bkg subtraction, vs eta-phi trigger.
  
  TH2F *   fhConeSumPtEtaUENormCluster;                //!<! Cluster Sum Pt in the normalized eta UE cone vs pT trigger.
  TH2F *   fhConeSumPtPhiUENormCluster;                //!<! Cluster Sum Pt in the normalized phi UE cone vs pT trigger.
  TH2F *   fhConeSumPtEtaUESubCluster;                 //!<! Cluster Sum Pt in the cone after bkg subtraction, vs pT trigger.
  TH2F *   fhConeSumPtPhiUESubCluster;                 //!<! Cluster Sum Pt in the cone after bkg subtraction, vs pT trigger.
  TH2F *   fhConeSumPtEtaUESubClusterTrigEtaPhi;       //!<! Cluster Sum Pt in the cone after bkg subtraction, vs eta-phi trigger.
  TH2F *   fhConeSumPtPhiUESubClusterTrigEtaPhi;       //!<! Cluster Sum Pt in the cone after bkg subtraction, vs eta-phi trigger.

  TH2F *   fhConeSumPtEtaUESubCell;                    //!<! Cell Sum amplitude in the cone after bkg subtraction, vs pT trigger.
  TH2F *   fhConeSumPtPhiUESubCell;                    //!<! Cell Sum amplitude in the cone after bkg subtraction, vs pT trigger.
  TH2F *   fhConeSumPtEtaUESubCellTrigEtaPhi;          //!<! Cell Sum amplitude in the cone after bkg subtraction, vs eta-phi trigger.
  TH2F *   fhConeSumPtPhiUESubCellTrigEtaPhi;          //!<! Cell Sum amplitude in the cone after bkg subtraction, vs eta-phi trigger.
  
  TH2F *   fhConeSumPtEtaUENormTrack;                  //!<! Track Sum Pt in the normalized eta UE cone vs pT trigger.
  TH2F *   fhConeSumPtPhiUENormTrack;                  //!<! Track Sum Pt in the normalized phi UE cone vs pT trigger.
  TH2F *   fhConeSumPtEtaUESubTrack;                   //!<! Track Sum Pt in the cone after bkg subtraction, vs pT trigger.
  TH2F *   fhConeSumPtPhiUESubTrack;                   //!<! Track Sum Pt in the cone after bkg subtraction, vs pT trigger.
  TH2F *   fhConeSumPtEtaUESubTrackTrigEtaPhi;         //!<! Track Sum Pt in the cone after bkg subtraction, vs eta-phi trigger.
  TH2F *   fhConeSumPtPhiUESubTrackTrigEtaPhi;         //!<! Track Sum Pt in the cone after bkg subtraction, vs eta-phi trigger.
  
  TH2F *   fhFractionTrackOutConeEta;                  //!<! Fraction of cone out of tracks acceptance in eta.
  TH2F *   fhFractionTrackOutConeEtaTrigEtaPhi;        //!<! Fraction of cone out of tracks acceptance in eta, vs trigger eta-phi.
  TH2F *   fhFractionClusterOutConeEta;                //!<! Fraction of cone out of clusters acceptance in eta.
  TH2F *   fhFractionClusterOutConeEtaTrigEtaPhi;      //!<! Fraction of cone out of clusters acceptance in eta, vs trigger eta-phi.
  TH2F *   fhFractionClusterOutConePhi;                //!<! Fraction of cone out of clusters acceptance in phi.
  TH2F *   fhFractionClusterOutConePhiTrigEtaPhi;      //!<! Fraction of cone out of clusters acceptance in phi, vs trigger eta-phi.
  
  TH2F *   fhFractionCellOutConeEta;                   //!<! Fraction of cone out of cells acceptance in eta.
  TH2F *   fhFractionCellOutConeEtaTrigEtaPhi;         //!<! Fraction of cone out of cells acceptance in eta, vs trigger eta-phi.
  TH2F *   fhFractionCellOutConePhi;                   //!<! Fraction of cone out of cells acceptance in phi.
  TH2F *   fhFractionCellOutConePhiTrigEtaPhi;         //!<! Fraction of cone out of cells acceptance in phi, vs trigger eta-phi.
  
  TH2F *   fhConeSumPtClustervsTrack ;                 //!<! Cluster vs tracks Sum Pt Sum Pt in the cone.
  TH2F *   fhConeSumPtClusterTrackFrac ;               //!<! Cluster / tracks Sum Pt Sum Pt in the cone.
  TH2F *   fhConeSumPtEtaUESubClustervsTrack ;         //!<! Cluster vs tracks Sum Pt Sum Pt in the cone, after subtraction in eta band.
  TH2F *   fhConeSumPtPhiUESubClustervsTrack ;         //!<! Cluster vs tracks Sum Pt Sum Pt in the cone, after subtraction in phi band.
  TH2F *   fhConeSumPtCellvsTrack;                     //!<! Cell vs tracks Sum Pt Sum Pt in the cone.
  TH2F *   fhConeSumPtEtaUESubCellvsTrack ;            //!<! Cell vs tracks Sum Pt Sum Pt in the cone, after subtraction in eta band.
  TH2F *   fhConeSumPtPhiUESubCellvsTrack ;            //!<! Cell vs tracks Sum Pt Sum Pt in the cone, after subtraction in phi band.

  TH2F *   fhEtaBandClustervsTrack ;                   //!<! Accumulated pT in Eta band to estimate UE in cone, clusters vs tracks.
  TH2F *   fhPhiBandClustervsTrack ;                   //!<! Accumulated pT in Phi band to estimate UE in cone, clusters vs tracks.
  TH2F *   fhEtaBandNormClustervsTrack ;               //!<! Accumulated pT in Eta band to estimate UE in cone, normalized to cone size, clusters vs tracks.
  TH2F *   fhPhiBandNormClustervsTrack ;               //!<! Accumulated pT in Phi band to estimate UE in cone, normalized to cone size, clusters vs tracks.
  TH2F *   fhEtaBandCellvsTrack ;                      //!<! Accumulated pT in Eta band to estimate UE in cone, cells vs tracks.
  TH2F *   fhPhiBandCellvsTrack ;                      //!<! Accumulated pT in Phi band to estimate UE in cone, cells vs tracks.
  TH2F *   fhEtaBandNormCellvsTrack ;                  //!<! Accumulated pT cell in Eta band to estimate UE in cone, normalized to cone size, clusters vs tracks.
  TH2F *   fhPhiBandNormCellvsTrack ;                  //!<! Accumulated pT cell in Phi band to estimate UE in cone, normalized to cone.

  TH2F *   fhConeSumPtSubvsConeSumPtTotPhiTrack;       //!<! Tracks, phi band: sum pT in cone after bkg sub vs sum pT in cone before bkg sub.
  TH2F *   fhConeSumPtSubNormvsConeSumPtTotPhiTrack;   //!<! Tracks, phi band: sum pT in cone after bkg sub normalized by sum pT in cone before bkg sub vs sum pT in cone before bkg sub.
  TH2F *   fhConeSumPtSubvsConeSumPtTotEtaTrack;       //!<! Tracks, eta band: sum pT in cone after bkg sub vs sum pT in cone before bkg sub
  TH2F *   fhConeSumPtSubNormvsConeSumPtTotEtaTrack;   //!<! Tracks, eta band: sum pT in cone after bkg sub normalized by sum pT in cone before bkg sub vs sum pT in cone before bkg sub.
  TH2F *   fhConeSumPtSubvsConeSumPtTotPhiCluster;     //!<! Clusters, phi band: sum pT in cone after bkg sub vs sum pT in cone before bkg sub.
  TH2F *   fhConeSumPtSubNormvsConeSumPtTotPhiCluster; //!<! Clusters, phi band: sum pT in cone after bkg sub normalized by sum pT in cone before bkg sub vs sum pT in cone before bkg sub.
  TH2F *   fhConeSumPtSubvsConeSumPtTotEtaCluster;     //!<! Clusters, eta band: sum pT in cone after bkg sub vs sum pT in cone before bkg sub
  TH2F *   fhConeSumPtSubNormvsConeSumPtTotEtaCluster; //!<! Clusters, eta band: sum pT in cone after bkg sub normalized by sum pT in cone before bkg sub vs sum pT in cone before bkg sub.
  TH2F *   fhConeSumPtSubvsConeSumPtTotPhiCell;        //!<! Cells, phi band: sum pT in cone after bkg sub vs sum pT in cone before bkg sub.
  TH2F *   fhConeSumPtSubNormvsConeSumPtTotPhiCell;    //!<! Cells, phi band: sum pT in cone after bkg sub normalized by sum pT in cone before bkg sub vs sum pT in cone before bkg sub.
  TH2F *   fhConeSumPtSubvsConeSumPtTotEtaCell;        //!<! Cells, eta band: sum pT in cone after bkg sub vs sum pT in cone before bkg sub.
  TH2F *   fhConeSumPtSubNormvsConeSumPtTotEtaCell;    //!<! Cells, eta band: sum pT in cone after bkg sub normalized by sum pT in cone before bkg sub vs sum pT in cone before bkg sub.
  TH2F *   fhConeSumPtVSUETracksEtaBand;               //!<! Tracks, eta band: sum pT in cone vs bkg to subtract.
  TH2F *   fhConeSumPtVSUETracksPhiBand;               //!<! Tracks, phi band:  sum pT in cone vs bkg to subtract.
  TH2F *   fhConeSumPtVSUEClusterEtaBand;              //!<! Clusters, eta band: sum pT in cone vs bkg to subtract.
  TH2F *   fhConeSumPtVSUEClusterPhiBand;              //!<! Clusters, phi band:  sum pT in cone vs bkg to subtract.
  
  // MC
  
  TH2F *   fhEtaPrimMC  [fgkNmcPrimTypes];             //!<! Pt vs Eta of generated photon.
  TH2F *   fhPhiPrimMC  [fgkNmcPrimTypes];             //!<! Pt vs Phi of generated photon.
  TH1F *   fhEPrimMC    [fgkNmcPrimTypes];             //!<! Number of generated photon vs E.
  TH1F *   fhPtPrimMC   [fgkNmcPrimTypes];             //!<! Number of generated photon vs pT.
  TH1F *   fhPtPrimMCiso[fgkNmcPrimTypes];             //!<! Number of generated isolated photon vs pT.
  
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

  TH1F *   fhPtNoIsoMC  [fgkNmcTypes];                 //!<! Number of not isolated mcTypes particle.
  TH1F *   fhPtIsoMC    [fgkNmcTypes];                 //!<! Number of isolated mcTypes particle.
  TH2F *   fhPhiIsoMC   [fgkNmcTypes];                 //!<! phi of isolated mcTypes particle.
  TH2F *   fhEtaIsoMC   [fgkNmcTypes];                 //!<! eta of isolated mcTypes particle.
  
  TH1F *   fhPtDecayMC  [2][AliNeutralMesonSelection::fgkMaxNDecayBits][fgkNmcTypes] ; //!<! Number of (not) isolated Pi0 decay particles (invariant mass tag) for a mcTypes particle.
  
  TH2F *   fhPtLambda0MC    [fgkNmcTypes][2];           //!<! Shower shape of (non) isolated candidates originated by mcTypes particle (do not apply SS cut previously).
  TH2F *   fhPtLambda0MCConv[fgkNmcTypes][2];           //!<! Shower shape of (non) isolated candidates originated by mcTypes particle that converted (do not apply SS cut previously).

  TH2F *   fhPtLambda0MCWith1Overlap    [fgkNmcTypes][2];           //!<! Shower shape of (non) isolated candidates originated by mcTypes particle (do not apply SS cut previously). At least one overlap from other particles.
  TH2F *   fhPtLambda0MCConvWith1Overlap[fgkNmcTypes][2];           //!<! Shower shape of (non) isolated candidates originated by mcTypes particle that converted (do not apply SS cut previously). At least one overlap from other particles.

  TH2F *   fhPtLambda0MCWithNoOverlap    [fgkNmcTypes][2];          //!<! Shower shape of (non) isolated candidates originated by mcTypes particle (do not apply SS cut previously). More tha one overlap from other particles.
  TH2F *   fhPtLambda0MCConvWithNoOverlap[fgkNmcTypes][2];          //!<! Shower shape of (non) isolated candidates originated by mcTypes particle that converted (do not apply SS cut previously). More tha one overlap from other particles.
  
  TH2F *   fhPtNOverlap    [fgkNmcTypes][2];                        //!<! Number of overlaps of (non) isolated candidates originated by mcTypes (do not apply SS cut previously). More tha one overlap from other particles.
  TH2F *   fhPtNOverlapConv[fgkNmcTypes][2];                        //!<! Number of overlaps of (non) isolated candidates originated by mcTypes particle that converted (do not apply SS cut previously). More tha one overlap from other particles.
  
  // Multiple cut analysis
  TH2F *   fhSumPtLeadingPt[5] ;                       //!<! Sum Pt in the cone.
  TH2F *   fhPtLeadingPt[5] ;                          //!<! Particle Pt in the cone.
  TH2F *   fhPerpSumPtLeadingPt[5] ;                   //!<! Sum Pt in the cone at the perpendicular phi region to trigger axis  (phi +90).
  TH2F *   fhPerpPtLeadingPt[5];                       //!<! Sum Pt in the cone at the perpendicular phi region to trigger axis  (phi +90).

  TH1F *   fhPtThresIsolated[5][5] ;                   //!<! Isolated particle with pt threshold.
  TH1F *   fhPtFracIsolated[5][5] ;                    //!<! Isolated particle with pt threshold frac.
  TH1F *   fhSumPtIsolated[5][5] ;                     //!<! Isolated particle with threshold on cone pt sum.
  
  TH2F *   fhEtaPhiPtThresIso[5][5] ;                  //!<! eta vs phi of isolated particles with pt threshold.
  TH2F *   fhEtaPhiPtThresDecayIso[5][5] ;             //!<! eta vs phi of isolated particles with pt threshold, only for decay bit fDecayBits[0].
  TH1F *   fhPtPtThresDecayIso[5][5] ;                 //!<! Number of isolated Pi0 decay particles (invariant mass tag) with pt threshold, only for decay bit fDecayBits[0].
  
  TH2F *   fhEtaPhiPtFracIso[5][5] ;                   //!<! eta vs phi of isolated particles with pt frac.
  TH2F *   fhEtaPhiPtFracDecayIso[5][5] ;              //!<! eta vs phi of isolated particles with pt frac, only for decay bit fDecayBits[0].
  TH1F *   fhPtPtFracDecayIso[5][5] ;                  //!<! Number of isolated Pi0 decay particles (invariant mass tag) with pt fra, only for decay bit fDecayBits[0].

  TH2F *   fhEtaPhiPtSumIso[5][5] ;                    //!<! eta vs phi of isolated particles with pt sum.
  TH2F *   fhEtaPhiPtSumDecayIso[5][5] ;               //!<! eta vs phi of isolated particles with pt sum, only for decay bit fDecayBits[0].
  TH1F *   fhPtPtSumDecayIso[5][5] ;                   //!<! Number of isolated Pi0 decay particles (invariant mass tag) with pt sum, only for decay bit fDecayBits[0].
  
  TH2F *   fhEtaPhiSumDensityIso[5][5];                //!<! Isolated particle with threshold on cone sum density.
  TH2F *   fhEtaPhiSumDensityDecayIso[5][5];           //!<! Isolated particle with threshold on cone sum density, only for decay bit fDecayBits[0].
  TH1F *   fhPtSumDensityIso[5][5];                    //!<! Isolated particle with threshold on cone sum density.
  TH1F *   fhPtSumDensityDecayIso[5][5];               //!<! Isolated decay particle with threshold on cone sum density, only for decay bit fDecayBits[0].
  
  TH1F *   fhPtFracPtSumIso[5][5] ;                    //!<! Number of isolated Pi0 decay particles (invariant mass tag) with pt sum.
  TH1F *   fhPtFracPtSumDecayIso[5][5] ;               //!<! Number of isolated Pi0 decay particles (invariant mass tag) with pt sum, only for decay bit fDecayBits[0].
  TH2F *   fhEtaPhiFracPtSumIso[5][5];                 //!<! Isolated particle with threshold on cone sum density.
  TH2F *   fhEtaPhiFracPtSumDecayIso[5][5];            //!<! Isolated particle with threshold on cone sum density, only for decay bit fDecayBits[0].
 
  // Multiple cut MC
  TH1F *   fhPtThresIsolatedMC[fgkNmcTypes][5][5];     //!<! Isolated mcTypes particle with pt threshold.
  TH1F *   fhPtFracIsolatedMC [fgkNmcTypes][5][5];     //!<! Isolated mcTypes particle with pt frac.
  TH1F *   fhSumPtIsolatedMC  [fgkNmcTypes][5][5];     //!<! Isolated mcTypes particle with threshold on cone pt sum.
  TH2F *   fhSumPtLeadingPtMC [fgkNmcTypes][5];        //!<! mcTypes particle for sum Pt, different cone.

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
  TH2F *   fhELambda0[2];                              //!<! Shower shape of (non) isolated photons (do not apply SS cut previously).
  TH2F *   fhPtLambda0[2];                             //!<! Shower shape of (non) isolated photons (do not apply SS cut previously).
//TH2F *   fhELambda1[2];                              //!<! Shower shape of (non) isolated photons (do not apply SS cut previously).
  TH2F *   fhELambda0TRD[2];                           //!<! Shower shape of (non) isolated photons, SM behind TRD (do not apply SS cut previously).
  TH2F *   fhPtLambda0TRD[2];                          //!<! Shower shape of (non) isolated photons, SM behind TRD (do not apply SS cut previously).
//TH2F *   fhELambda1TRD[2];                           //!<! Shower shape of (non) isolated photons, SM behind TRD (do not apply SS cut previously).

  /// Candidate Pt distribution depending on bin of cone leading particle.
  TH1F **  fhPtLeadConeBin ;                           //![fNBkgBin]
    
  /// Candidate Pt distribution depending on bin of cone sum pt.
  TH1F **  fhSumPtConeBin  ;                           //![fNBkgBin]
    
  /// Candidate Pt distribution depending on bin of cone leading particle, per MC particle.
  TH1F **  fhPtLeadConeBinMC ;                         //![fNBkgBin*fgkNmcTypes]
    
  /// Candidate Pt distribution depending on bin of cone sum pt, per MC particle.
  TH1F **  fhSumPtConeBinMC  ;                         //![fNBkgBin*fgkNmcTypes]

  /// Candidate Pt distribution depending on bin of cone leading particle, tagged as decay.
  TH1F **  fhPtLeadConeBinDecay ;                      //![fNBkgBin*fNDecayBits]
    
  /// Candidate Pt distribution depending on bin of cone sum pt, tagged as decay.
  TH1F **  fhSumPtConeBinDecay  ;                      //![fNBkgBin*fNDecayBits]
  
  /// Candidate shower shape distribution depending on bin of cone leading particle.
  TH2F **  fhPtLeadConeBinLambda0 ;                    //![fNBkgBin]
    
  /// Candidate shower shape distribution depending on bin of cone sum pt.
  TH2F **  fhSumPtConeBinLambda0  ;                    //![fNBkgBin]
    
  /// Candidate shower shape distribution depending on bin of cone leading particle, per MC particle.
  TH2F **  fhPtLeadConeBinLambda0MC ;                  //![fNBkgBin*fgkNmcTypes]
    
  /// Candidate shower shape distribution depending on bin of cone sum pt, per MC particle.
  TH2F **  fhSumPtConeBinLambda0MC  ;                  //![fNBkgBin*fgkNmcTypes]

  /// Candidate pt bin, distribution of cone leading particle pt.
  TH1F **  fhPtTrigBinPtLeadCone ;                     //![fNPtTrigBin]
    
  /// Candidate pt bin, distribution of cone sum particle pt.
  TH1F **  fhPtTrigBinSumPtCone  ;                     //![fNPtTrigBin]
  /// Candidate pt bin, distribution of cone sum track pt.
  TH1F **  fhPtTrigBinSumPtTrackCone  ;                //![fNPtTrigBin]
  /// Candidate pt bin, distribution of cone sum cluster pt.
  TH1F **  fhPtTrigBinSumPtClusterCone  ;              //![fNPtTrigBin]

  /// Candidate pt bin, distribution of cone leading particle pt, per MC particle.
  TH1F **  fhPtTrigBinPtLeadConeMC ;                   //![fNPtTrigBin*fgkNmcTypes]
    
  /// Candidate pt bin, distribution of cone sum particle pt, per MC particle.
  TH1F **  fhPtTrigBinSumPtConeMC  ;                   //![fNPtTrigBin*fgkNmcTypes]
  /// Candidate pt bin, distribution of cone sum track pt, per MC particle.
  TH1F **  fhPtTrigBinSumPtTrackConeMC  ;              //![fNPtTrigBin*fgkNmcTypes]
  /// Candidate pt bin, distribution of cone sum cluster pt, per MC particle.
  TH1F **  fhPtTrigBinSumPtClusterConeMC  ;            //![fNPtTrigBin*fgkNmcTypes]
  
  /// Candidate pt bin, distribution of cone leading particle pt, tagged as decay.
  TH1F **  fhPtTrigBinPtLeadConeDecay ;                //![fNBkgBin*fNDecayBits]
  
  /// Candidate pt bin, distribution of cone sum particle pt, tagged as decay.
  TH1F **  fhPtTrigBinSumPtConeDecay  ;                //![fNBkgBin*fNDecayBits]
  /// Candidate pt bin, distribution of cone sum track pt, tagged as decay.
  TH1F **  fhPtTrigBinSumPtTrackConeDecay  ;           //![fNBkgBin*fNDecayBits]
  /// Candidate pt bin, distribution of cone sum cluster pt, tagged as decay.
  TH1F **  fhPtTrigBinSumPtClusterConeDecay  ;         //![fNBkgBin*fNDecayBits]

  /// Candidate shower shape distribution depending vs cone leading particle in pT trigger bins.
  TH2F **  fhPtTrigBinLambda0vsPtLeadCone ;            //![fNPtTrigBin]
    
  /// Candidate shower shape distribution depending vs of cone sum pt in pT trigger bins.
  TH2F **  fhPtTrigBinLambda0vsSumPtCone  ;            //![fNPtTrigBin]
  /// Candidate shower shape distribution depending vs of cone sum track pt in pT trigger bins.
  TH2F **  fhPtTrigBinLambda0vsSumPtTrackCone  ;       //![fNPtTrigBin]
  /// Candidate shower shape distribution depending vs of cone sum cluster pt in pT trigger bins.
  TH2F **  fhPtTrigBinLambda0vsSumPtClusterCone  ;     //![fNPtTrigBin]

  
  /// Candidate shower shape distribution depending vs cone leading particle in pT trigger bins, per MC particle.
  TH2F **  fhPtTrigBinLambda0vsPtLeadConeMC ;          //![fNPtTrigBin*fgkNmcTypes]
    
  /// Candidate shower shape distribution depending vs cone sum pt in pT trigger bins, per MC particle.
  TH2F **  fhPtTrigBinLambda0vsSumPtConeMC  ;          //![fNPtTrigBin*fgkNmcTypes]
  /// Candidate shower shape distribution depending vs cone sum pt track in pT trigger bins, per MC particle.
  TH2F **  fhPtTrigBinLambda0vsSumPtTrackConeMC  ;     //![fNPtTrigBin*fgkNmcTypes]
  /// Candidate shower shape distribution depending vs cone sum pt cluster in pT trigger bins, per MC particle.
  TH2F **  fhPtTrigBinLambda0vsSumPtClusterConeMC  ;   //![fNPtTrigBin*fgkNmcTypes]

  /// Candidate shower shape distribution depending vs cone sum pt in pT trigger bins, per MC particle, candidate without overlaps
  TH2F **  fhPtTrigBinLambda0vsSumPtConeMCNoOverlap  ;          //![fNPtTrigBin*fgkNmcTypes]
  /// Candidate shower shape distribution depending vs cone sum pt track in pT trigger bins, per MC particle,, candidate without overlaps
  TH2F **  fhPtTrigBinLambda0vsSumPtTrackConeMCNoOverlap  ;     //![fNPtTrigBin*fgkNmcTypes]
  /// Candidate shower shape distribution depending vs cone sum pt cluster in pT trigger bins, per MC particle, candidate without overlaps
  TH2F **  fhPtTrigBinLambda0vsSumPtClusterConeMCNoOverlap  ;   //![fNPtTrigBin*fgkNmcTypes]

  /// Candidate shower shape distribution depending vs cone sum pt in pT trigger bins, per MC particle, candidate with at least 1 overlap
  TH2F **  fhPtTrigBinLambda0vsSumPtConeMC1Overlap  ;          //![fNPtTrigBin*fgkNmcTypes]
  /// Candidate shower shape distribution depending vs cone sum pt track in pT trigger bins, per MC particle, candidate with at least 1 overlap
  TH2F **  fhPtTrigBinLambda0vsSumPtTrackConeMC1Overlap  ;     //![fNPtTrigBin*fgkNmcTypes]
  /// Candidate shower shape distribution depending vs cone sum pt cluster in pT trigger bins, per MC particle, candidate with at least 1 overlap
  TH2F **  fhPtTrigBinLambda0vsSumPtClusterConeMC1Overlap  ;   //![fNPtTrigBin*fgkNmcTypes]
  
  // Local maxima
  TH2F *   fhNLocMax[2];                               //!<! Number of maxima in selected clusters.
  TH2F *   fhELambda0LocMax1[2] ;                      //!<! E vs lambda0 of selected cluster, 1 local maxima in cluster.
  TH2F *   fhELambda1LocMax1[2] ;                      //!<! E vs lambda1 of selected cluster, 1 local maxima in cluster.
  TH2F *   fhELambda0LocMax2[2] ;                      //!<! E vs lambda0 of selected cluster, 2 local maxima in cluster.
  TH2F *   fhELambda1LocMax2[2] ;                      //!<! E vs lambda1 of selected cluster, 2 local maxima in cluster.
  TH2F *   fhELambda0LocMaxN[2] ;                      //!<! E vs lambda0 of selected cluster, N>2 local maxima in cluster.
  TH2F *   fhELambda1LocMaxN[2] ;                      //!<! E vs lambda1 of selected cluster, N>2 local maxima in cluster.
  
  // Pile-up
  TH1F *   fhEIsoPileUp[7] ;                           //!<! Number of isolated particles.
  TH1F *   fhPtIsoPileUp[7] ;                          //!<! Number of isolated particles.
  TH1F *   fhENoIsoPileUp[7] ;                         //!<! Number of not isolated particles.
  TH1F *   fhPtNoIsoPileUp[7] ;                        //!<! Number of not isolated particles.
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
  TH2F *   fhPtClusterInConePerRCutLargePtTrig;          //!<! Clusters Pt in the cone for different cone sizes, x axis. Trigger pT > 10 GeV fixed
  TH2F *   fhPtTrackInConePerRCut;                       //!<! Tracks Pt in the cone for different cone sizes, x axis.
  TH2F *   fhPtTrackInConePerRCutLargePtTrig;            //!<! Tracks Pt in the cone for different cone sizes, x axis. Trigger pT > 10 GeV fixed
  TH2F *   fhConeSumPtClusterPerRCut;                    //!<! Clusters Sum Pt in the cone for different cone sizes, x axis.
  TH2F *   fhConeSumPtClusterPerRCutLargePtTrig;         //!<! Clusters Sum Pt in the cone for different cone sizes, x axis. Trigger pT > 10 GeV fixed
  TH2F *   fhConeSumPtTrackPerRCut;                      //!<! Tracks Sum Pt in the cone for different cone sizes, x axis.
  TH2F *   fhConeSumPtTrackPerRCutLargePtTrig;           //!<! Tracks Sum Pt in the cone for different cone sizes, x axis. Trigger pT > 10 GeV fixed

  TH2F *   fhConeNClusterPerMinPtCut;                    //!<! N Clusters in the cone for different min pT cuts, x axis.
  TH2F *   fhConeNClusterPerMinPtCutLargePtTrig;         //!<! N Clusters in the cone for different min pT cuts, x axis. Trigger pT > 10 GeV fixed
  TH2F *   fhConeNTrackPerMinPtCut;                      //!<! N Tracks in the cone for different min pT cuts, x axis.
  TH2F *   fhConeNTrackPerMinPtCutLargePtTrig;           //!<! N Tracks in the cone for different min pT cuts, x axis. Trigger pT > 10 GeV fixed
  TH2F *   fhPerpConeNTrackPerMinPtCut;                  //!<! N Tracks in the perpendicular cone for different min pT cuts, x axis.
  TH2F *   fhPerpConeNTrackPerMinPtCutLargePtTrig;       //!<! N Tracks in the perpendicular cone for different min pT cuts, x axis. Trigger pT > 10 GeV fixed

  TH2F *   fhConeSumPtClusterPerMinPtCut;                //!<! Clusters Sum Pt in the cone for different min pT cuts, x axis.
  TH2F *   fhConeSumPtClusterPerMinPtCutLargePtTrig;     //!<! Clusters Sum Pt in the cone for different min pT cuts, x axis. Trigger pT > 10 GeV fixed
  TH2F *   fhConeSumPtTrackPerMinPtCut;                  //!<! Tracks Sum Pt in the cone for different min pT cuts, x axis.
  TH2F *   fhConeSumPtTrackPerMinPtCutLargePtTrig;       //!<! Tracks Sum Pt in the cone for different min pT cuts, x axis. Trigger pT > 10 GeV fixed
  TH2F *   fhPerpConeSumPtTrackPerMinPtCut;              //!<! Tracks Sum Pt in the perpendicular cone for different min pT cuts, x axis.
  TH2F *   fhPerpConeSumPtTrackPerMinPtCutLargePtTrig;   //!<! Tracks Sum Pt in the perpendicular cone for different min pT cuts, x axis. Trigger pT > 10 GeV fixed
  
  TH2F *   fhConeSumPtClusterPerMaxPtCut;                //!<! Clusters Sum Pt in the cone for different max pT cuts, x axis.
  TH2F *   fhConeSumPtClusterPerMaxPtCutLargePtTrig;     //!<! Clusters Sum Pt in the cone for different max pT cuts, x axis. Trigger pT > 10 GeV fixed
  TH2F *   fhConeSumPtTrackPerMaxPtCut;                  //!<! Tracks Sum Pt in the cone for different max pT cuts, x axis.
  TH2F *   fhConeSumPtTrackPerMaxPtCutLargePtTrig;       //!<! Tracks Sum Pt in the cone for different max pT cuts, x axis. Trigger pT > 10 GeV fixed
  TH2F *   fhConeSumPtTrackPerEtaCut;                    //!<! Tracks Sum Pt in the cone for different min eta cuts, x axis.
  TH2F *   fhConeSumPtTrackPerEtaCutLargePtTrig;         //!<! Tracks Sum Pt in the cone for different min eta cuts, x axis. Trigger pT > 10 GeV fixed
  
  TH2F *   fhPtClusterInConePerNCellCut;                 //!<! Clusters Pt in the cone for different min cluster n cell cut, x axis.
  TH2F *   fhPtClusterInConePerNCellCutLargePtTrig;      //!<! Clusters Pt in the cone for different min cluster n cell cut, x axis. Trigger pT > 10 GeV fixed
  TH2F *   fhPtTrackInConePerNCellCut;                   //!<! Tracks Pt in the cone for different min cluster n cell cut, x axis.
  TH2F *   fhPtTrackInConePerNCellCutLargePtTrig;        //!<! Tracks Pt in the cone for different min cluster n cell cut, x axis. Trigger pT > 10 GeV fixed

  TH2F *   fhConeSumPtClusterPerNCellCut;                //!<! Clusters Sum Pt in the cone for different min cluster n cell cut, x axis.
  TH2F *   fhConeSumPtClusterPerNCellCutLargePtTrig;     //!<! Clusters Sum Pt in the cone for different min cluster n cell cut, x axis. Trigger pT > 10 GeV fixed
  TH2F *   fhConeSumPtTrackPerNCellCut;                  //!<! Tracks Sum Pt in the cone for different min cluster n cell cut, x axis.
  TH2F *   fhConeSumPtTrackPerNCellCutLargePtTrig;       //!<! Tracks Sum Pt in the cone for different min cluster n cell cut, x axis. Trigger pT > 10 GeV fixed

  TH2F *   fhPtClusterInConePerExoCut;                   //!<! Clusters Pt in the cone for different exoticity cut, x axis.
  TH2F *   fhPtClusterInConePerExoCutLargePtTrig;        //!<! Clusters Pt in the cone for different exoticity cut, x axis. Trigger pT > 10 GeV fixed
  TH2F *   fhPtTrackInConePerExoCut;                     //!<! Tracks Pt in the cone for different exoticity cut, x axis.
  TH2F *   fhPtTrackInConePerExoCutLargePtTrig;          //!<! Tracks Pt in the cone for different exoticity cut, x axis. Trigger pT > 10 GeV fixed

  TH2F *   fhConeSumPtClusterPerExoCut;                  //!<! Clusters Sum Pt in the cone for different exoticity cut, x axis.
  TH2F *   fhConeSumPtClusterPerExoCutLargePtTrig;       //!<! Clusters Sum Pt in the cone for different exoticity cut, x axis. Trigger pT > 10 GeV fixed
  TH2F *   fhConeSumPtTrackPerExoCut;                    //!<! Tracks Sum Pt in the cone for different exoticity cut, x axis.
  TH2F *   fhConeSumPtTrackPerExoCutLargePtTrig;         //!<! Tracks Sum Pt in the cone for different exoticity cut, x axis. Trigger pT > 10 GeV fixed

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
  
  /// Copy constructor not implemented.
  AliAnaParticleIsolation(              const AliAnaParticleIsolation & iso) ;
    
  /// Assignment operator not implemented.
  AliAnaParticleIsolation & operator = (const AliAnaParticleIsolation & iso) ;
  
  /// \cond CLASSIMP
  ClassDef(AliAnaParticleIsolation,39) ;
  /// \endcond

} ;


#endif //ALIANAPARTICLEISOLATION_H



