#ifndef ALIANAPARTICLEISOLATION_H
#define ALIANAPARTICLEISOLATION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________

// Class for the analysis of particle isolation
// Input is selected particles put in AOD branch (AliAODPWG4ParticleCorrelation)
//
//  Class created from old AliPHOSGammaJet
//  (see AliRoot versions previous Release 4-09)

//-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
class TH2F;
class TList ;
class TObjString;

// --- ANALYSIS system ---
#include "AliAnaCaloTrackCorrBaseClass.h"
class AliAODPWG4Particle;
class AliAODPWG4ParticleCorrelation ;


class AliAnaParticleIsolation : public AliAnaCaloTrackCorrBaseClass {

 public:   
  AliAnaParticleIsolation() ; // default ctor
  virtual ~AliAnaParticleIsolation() { ; } //virtual dtor

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


  void         CalculateNormalizeUEBandPerUnitArea(AliAODPWG4ParticleCorrelation * pCandidate,
                                                   Float_t coneptsumCluster,       Float_t coneptsumCell,     Float_t coneptsumTrack,
                                                   Float_t &etaBandptsumTrackNorm, Float_t &etaBandptsumClusterNorm ) ;
  
  TObjString * GetAnalysisCuts() ;
  
  TList      * GetCreateOutputObjects() ;
  
  void         Init() ;

  void         InitParameters() ;
  
  void         MakeAnalysisFillAOD()  ;
  
  void         MakeAnalysisFillHistograms() ; 
  
  
  void         Print( const Option_t * opt ) const ;
 
  //Analysis specific methods 
  
  void         FillPileUpHistograms(Int_t clusterID) ;
  
  void         FillAcceptanceHistograms();
 
  void         FillTrackMatchingShowerShapeControlHistograms(AliAODPWG4ParticleCorrelation  * pCandidate,
                                                             Float_t coneptsum, Float_t coneleadpt, Int_t mcIndex) ;
  
  Bool_t       IsTriggerTheNearSideEventLeadingParticle(Int_t & idLeading);
  
  void         MakeSeveralICAnalysis( AliAODPWG4ParticleCorrelation * ph, Int_t mcIndex ) ;
  
  // Analysis Setters and Getters
  
  TString      GetCalorimeter()                const { return fCalorimeter       ; }
  TString      GetTriggerDetector()            const { return fIsoDetector       ; }
  Int_t        GetNCones()                     const { return fNCones            ; }
  Int_t        GetNPtThresFrac()               const { return fNPtThresFrac      ; }
  Float_t      GetConeSizes(Int_t i)           const { return fConeSizes[i]      ; }
  Float_t      GetPtThresholds(Int_t i)        const { return fPtThresholds[i]   ; }
  Float_t      GetSumPtThresholds(Int_t i)     const { return fSumPtThresholds[i]; }
  Float_t      GetPtFractions(Int_t i)         const { return fPtFractions[i]    ; }
  
  Int_t        GetMCIndex(Int_t mcTag);
  
  void         SetCalorimeter(TString & det)         { fCalorimeter     = det    ; }
  void         SetTriggerDetector(TString & det)     { fIsoDetector     = det    ; }
  void         SetNCones(Int_t ncs)                  { fNCones          = ncs    ; }
  void         SetNPtThresFrac(Int_t npt)            { fNPtThresFrac    = npt    ; }
  void         SetConeSizes(Int_t i, Float_t r)      { fConeSizes[i]    = r      ; }
  void         SetPtThresholds(Int_t i, Float_t pt)  { fPtThresholds[i] = pt     ; }
  void         SetPtFractions(Int_t i, Float_t pt)   { fPtFractions[i]  = pt     ; } 
  void 	       SetSumPtThresholds(Int_t i, Float_t pt){ fSumPtThresholds[i] = pt ; }

  Bool_t       IsReIsolationOn()               const { return fReMakeIC          ; }
  void         SwitchOnReIsolation()                 { fReMakeIC      = kTRUE    ; }
  void         SwitchOffReIsolation()                { fReMakeIC      = kFALSE   ; }
  
  Bool_t       IsSeveralIsolationOn()          const { return fMakeSeveralIC     ; }
  void         SwitchOnSeveralIsolation()            { fMakeSeveralIC = kTRUE    ; }
  void         SwitchOffSeveralIsolation()           { fMakeSeveralIC = kFALSE   ; }

  void         SwitchOnFillPileUpHistograms()        { fFillPileUpHistograms = kTRUE  ; }
  void         SwitchOffFillPileUpHistograms()       { fFillPileUpHistograms = kFALSE ; }    
  
  void         SwitchOnTMHistoFill()                 { fFillTMHisto   = kTRUE    ; }
  void         SwitchOffTMHistoFill()                { fFillTMHisto   = kFALSE   ; }
  
  void         SwitchOnSSHistoFill()                 { fFillSSHisto   = kTRUE    ; }
  void         SwitchOffSSHistoFill()                { fFillSSHisto   = kFALSE   ; }

  Bool_t       IsLeadingOnlyOn()               const { return fLeadingOnly       ; }
  void         SwitchOnLeadingOnly()                 { fLeadingOnly    = kTRUE   ; }
  void         SwitchOffLeadingOnly()                { fLeadingOnly    = kFALSE  ; }
  
  void         SwitchOnCheckNeutralClustersForLeading() { fCheckLeadingWithNeutralClusters = kTRUE  ; }
  void         SwitchOffCheckNeutralClustersForLeading(){ fCheckLeadingWithNeutralClusters = kFALSE ; }
  
  void         SwitchOnUEBandSubtractionHistoFill()  { fFillUEBandSubtractHistograms   = kTRUE    ; }
  void         SwitchOffUEBandSubtractionHistoFill() { fFillUEBandSubtractHistograms   = kFALSE   ; }

  void         SwitchOnCellHistoFill()               { fFillCellHistograms = kTRUE ; }
  void         SwitchOffCellHistoFill()              { fFillCellHistograms = kFALSE; }

  void         SwitchOnHighMultiplicityHistoFill()   { fFillHighMultHistograms = kTRUE ; }
  void         SwitchOffHighMultiplicityHistoFill()  { fFillHighMultHistograms = kFALSE; }

  void         SwitchOnNLMHistoFill()                { fFillNLMHistograms = kTRUE ; }
  void         SwitchOffNLMHistoFill()               { fFillNLMHistograms = kFALSE; }
  
  void         SwitchOnDecayTaggedHistoFill()        { fFillTaggedDecayHistograms = kTRUE ; }
  void         SwitchOffDecayTaggedHistoFill()       { fFillTaggedDecayHistograms = kFALSE; }
  void         SetNDecayBits(Int_t n)                { fNDecayBits = n               ; }
  void         SetDecayBits(Int_t i, UInt_t bit)     { if(i < 4) fDecayBits[i] = bit ; }
  
  void         SwitchOnBackgroundBinHistoFill()      { fFillBackgroundBinHistograms = kTRUE ; }
  void         SwitchOffBackgroundBinHistoFill()     { fFillBackgroundBinHistograms = kFALSE; }
  void         SetNBackgroundBins(Int_t n)           { if(n < 19) fNBkgBin = n ; }
  void         SetBackgroundLimits(Int_t i,Float_t l){ if(i <= fNBkgBin) fBkgBinLimit[i] = l; }

  void         SwitchOnPrimariesInConeSelection()    { fSelectPrimariesInCone = kTRUE ; }
  void         SwitchOffPrimariesInConeSelection()   { fSelectPrimariesInCone = kFALSE; }

  void         SwitchOnPrimariesPi0DecayStudy()      { fMakePrimaryPi0DecayStudy = kTRUE ; }
  void         SwitchOffPrimariesPi0DecayStudy()     { fMakePrimaryPi0DecayStudy = kFALSE; }
  
  // For primary histograms in arrays, index in the array, corresponding to a photon origin
  enum mcPrimTypes { kmcPrimPhoton = 0, kmcPrimPi0Decay = 1, kmcPrimOtherDecay  = 2,
                     kmcPrimPrompt = 3, kmcPrimFrag     = 4, kmcPrimISR         = 5, kmcPrimPi0 = 6 } ;
  static const Int_t fgkNmcPrimTypes = 7;
  
  // For histograms in arrays, index in the array, corresponding to any particle origin
  enum mcTypes     { kmcPhoton   = 0, kmcPrompt   = 1, kmcFragment = 2,
                     kmcPi0      = 3, kmcPi0Decay = 4, kmcEtaDecay = 5, kmcOtherDecay = 6,
                     kmcElectron = 7, kmcHadron   = 8                                     } ;
  static const Int_t fgkNmcTypes = 9;

 private:
  
  TString  fCalorimeter ;                         // Calorimeter where neutral particles in cone for isolation are;
  TString  fIsoDetector ;                         // Candidate particle for isolation detector ;
  Bool_t   fReMakeIC ;                            // Do isolation analysis
  Bool_t   fMakeSeveralIC ;                       // Do analysis for different IC
  Bool_t   fFillPileUpHistograms;                 // Fill pile-up related histograms
  Bool_t   fFillTMHisto;                          // Fill track matching plots
  Bool_t   fFillSSHisto;                          // Fill Shower shape plots
  Bool_t   fFillUEBandSubtractHistograms;         // Fill histograms working on the UE subtraction
  Bool_t   fFillCellHistograms;                   // Fill cell histograms
  Bool_t   fFillHighMultHistograms;               // Fill high multiplicity histograms
  Bool_t   fFillTaggedDecayHistograms;            // Fill histograms for clusters tagged as decay
  Int_t    fNDecayBits ;                          // in case of study of decay triggers, select the decay bit
  UInt_t   fDecayBits[4] ;                        // in case of study of decay triggers, select the decay bit
  Bool_t   fFillNLMHistograms;                    // Fill NLM histograms
  Bool_t   fLeadingOnly;                          // Do isolation with leading particle
  Bool_t   fCheckLeadingWithNeutralClusters;      // Compare the trigger candidate to Leading pT with the clusters pT, by default only charged
  Bool_t   fSelectPrimariesInCone;                // In primary particle isolation studies, select only particles in isolation cone within detector acceptance and E cut.
  Bool_t   fMakePrimaryPi0DecayStudy;             // Fill dedicated histograms for primary decay photons
  
  Bool_t   fFillBackgroundBinHistograms;          // Fill histograms for different bins in pt content of the cone
  Int_t    fNBkgBin;                              // Number of bins on pt content in cone
  Float_t  fBkgBinLimit[20];                      // Pt bin limits on pt content in the cone

  // Analysis data members for multiple cones and pt thresholds
  Int_t    fNCones ;                              //! Number of cone sizes to test
  Int_t    fNPtThresFrac ;                        //! Number of ptThres and ptFrac to test
  
  Float_t  fConeSizes[5] ;                        //! Array with cones to test
  Float_t  fPtThresholds[5] ;                     //! Array with pt thresholds to test
  Float_t  fPtFractions[5] ;                      //! Array with pt thresholds to test frac
  Float_t  fSumPtThresholds[5] ;                  //! Array with pt thresholds to test frac
  
  //Histograms  
  
  TH1F *   fhEIso ;                               //! Number of isolated particles vs energy
  TH1F *   fhPtIso ;                              //! Number of isolated particles vs pT
  TH2F *   fhPtCentralityIso ;                    //! centrality vs pT
  TH2F *   fhPtEventPlaneIso ;                    //! event plane angle vs pT
  TH2F *   fhPtNLocMaxIso ;                       //! Number of isolated particles vs NLM in cluster
  TH2F *   fhPhiIso ;                             //! Phi of isolated particles
  TH2F *   fhEtaIso ;                             //! eta of isolated particles
  TH2F *   fhEtaPhiIso ;                          //! eta vs phi of isolated particles
  TH2F *   fhEtaPhiNoIso ;                        //! eta vs phi of not isolated leading particles
  TH1F *   fhENoIso ;                             //! Number of not isolated leading particles vs Energy
  TH1F *   fhPtNoIso ;                            //! Number of not isolated leading particles vs pT
  TH2F *   fhPtNLocMaxNoIso ;                     //! Number of not isolated particles vs NLM in cluster
  TH1F *   fhPtDecayIso[4] ;                      //! Number of isolated Pi0 decay particles (invariant mass tag)
  TH1F *   fhPtDecayNoIso[4] ;                    //! Number of not isolated Pi0 decay leading particles (invariant mass tag)
  TH2F *   fhEtaPhiDecayIso[4] ;                  //! eta vs phi of isolated Pi0 decay particles
  TH2F *   fhEtaPhiDecayNoIso[4] ;                //! eta vs phi of not isolated leading Pi0 decay particles
  TH2F *   fhPtLambda0Decay[2][4];                //! Shower shape of (non) isolated leading Pi0 decay particles (do not apply SS cut previously)

  TH2F *   fhPtInCone ;                           //! Cluster/track Pt in the cone
  TH2F *   fhPtClusterInCone ;                    //! Cluster Pt in the cone
  TH2F *   fhPtCellInCone ;                       //! Cell amplitude in the cone
  TH2F *   fhPtTrackInCone ;                      //! Track Pt in the cone
  TH2F *   fhPtTrackInConeOtherBC ;               //! Track Pt in the cone, tracks out of main BC Time window
  TH2F *   fhPtTrackInConeOtherBCPileUpSPD ;      //! Track Pt in the cone, tracks out of main BC Time window
  TH2F *   fhPtTrackInConeBC0 ;                   //! Track Pt in the cone, tracks in BC=0
  TH2F *   fhPtTrackInConeVtxBC0 ;                //! Track Pt in the cone, tracks in BC=0
  TH2F *   fhPtTrackInConeBC0PileUpSPD ;          //! Track Pt in the cone, tracks in BC=0
  TH2F *   fhPtInConePileUp[7] ;                  //! Particle Pt in the cone, if event is from pile-up (SPD method)
  TH2F *   fhPtInConeCent ;                       //! Particle Pt in the cone versus centrality
  TH2F *   fhPerpConeSumPt ;                      //! Sum Pt in cone at the perpendicular phi region to trigger axis  (phi +90)
  TH2F *   fhPtInPerpCone ;                       //! Particle Pt  in cone at the perpendicular phi region to trigger axis  (phi +90)
  
  TH2F *   fhEtaPhiInConeCluster ;                //! Eta vs. phi of clusters in cone
  TH2F *   fhEtaPhiCluster ;                      //! Eta vs. phi of all clusters
  TH2F *   fhEtaPhiInConeTrack ;                  //! Eta vs. phi of tracks in cone
  TH2F *   fhEtaPhiTrack ;                        //! Eta vs. phi of all tracks
  
  TH2F *   fhEtaBandCluster ;                     //! Accumulated pT in Eta band to estimate UE in cone, only clusters
  TH2F *   fhPhiBandCluster ;                     //! Accumulated pT in Phi band to estimate UE in cone, only clusters
  TH2F *   fhEtaBandTrack   ;                     //! Accumulated pT in Eta band to estimate UE in cone, only tracks
  TH2F *   fhPhiBandTrack   ;                     //! Accumulated pT in Phi band to estimate UE in cone, only tracks
  TH2F *   fhEtaBandCell ;                        //! Accumulated pT in Eta band to estimate UE in cone, only cells
  TH2F *   fhPhiBandCell ;                        //! Accumulated pT in Phi band to estimate UE in cone, only cells

  TH2F *   fhConePtLead ;                         //! Cluster and tracks leading pt in the cone
  TH2F *   fhConePtLeadCluster ;                  //! Clusters leading pt in the cone
  TH2F *   fhConePtLeadTrack ;                    //! Tracks leading pt in the cone
  TH2F *   fhConePtLeadClustervsTrack;            //! Tracks vs Clusters leading pt
  TH2F *   fhConePtLeadClusterTrackFrac;          //! Trigger pt vs cluster/track leading pt
  
  TH2F *   fhConeSumPt ;                          //! Cluster and tracks Sum Pt Sum Pt in the cone
  TH2F *   fhConeSumPtCellTrack ;                 //! Cells and tracks Sum Pt Sum Pt in the cone
  TH2F *   fhConeSumPtCell ;                      //! Cells Sum Pt Sum Pt in the cone
  TH2F *   fhConeSumPtCluster ;                   //! Clusters Sum Pt Sum Pt in the cone
  TH2F *   fhConeSumPtTrack ;                     //! Tracks Sum Pt Sum Pt in the cone
  TH2F *   fhConeSumPtEtaBandUECluster;           //! Cluster Sum Pt in the eta band for clusters, before normalization
  TH2F *   fhConeSumPtPhiBandUECluster;           //! Cluster Sum Pt in the phi band for clusters, before normalization
  TH2F *   fhConeSumPtEtaBandUETrack;             //! Track Sum Pt in the eta band for tracks  , before normalization
  TH2F *   fhConeSumPtPhiBandUETrack;             //! Track Sum Pt in the phi badn for tracks  , before normalization
  TH2F *   fhConeSumPtEtaBandUECell;              //! Cell Sum amplitude in the eta band for cells, before normalization
  TH2F *   fhConeSumPtPhiBandUECell;              //! Cell Sum amplitude in the phi band for cells, before normalization

  TH2F *   fhConeSumPtTrigEtaPhi ;                //! Cluster and tracks Sum Pt Sum Pt in the cone, per eta-phi bin of trigger,
  TH2F *   fhConeSumPtCellTrackTrigEtaPhi ;       //! Cell and tracks Sum Pt Sum Pt in the cone, per eta-phi bin of trigger,
  TH2F *   fhConeSumPtEtaBandUEClusterTrigEtaPhi; //! Cluster Sum Pt in the eta band for clusters, per eta-phi bin of trigger,before normalization
  TH2F *   fhConeSumPtPhiBandUEClusterTrigEtaPhi; //! Cluster Sum Pt in the phi band for clusters, per eta-phi bin of trigger, before normalization
  TH2F *   fhConeSumPtEtaBandUETrackTrigEtaPhi;   //! Track Sum Pt in the eta band for tracks  , per eta-phi bin of trigger, before normalization
  TH2F *   fhConeSumPtPhiBandUETrackTrigEtaPhi;   //! Track Sum Pt in the phi badn for tracks  , per eta-phi bin of trigger, before normalization
  TH2F *   fhConeSumPtEtaBandUECellTrigEtaPhi;    //! Cluster Sum amplitude in the eta band for cells, per eta-phi bin of trigger, before normalization
  TH2F *   fhConeSumPtPhiBandUECellTrigEtaPhi;    //! Cluster Sum amplitude in the phi band for cells, per eta-phi bin of trigger, before normalization
  
  TH2F *   fhConeSumPtEtaUESub;                   //! Cluster and tracks Sum Pt in the cone after bkg subtraction, vs pT trigger
  TH2F *   fhConeSumPtPhiUESub;                   //! Cluster and tracks Sum Pt in the cone after bkg subtraction, vs pT trigger
  TH2F *   fhConeSumPtEtaUESubTrigEtaPhi;         //! Cluster and tracks Sum Pt in the cone after bkg subtraction, vs eta-phi trigger
  TH2F *   fhConeSumPtPhiUESubTrigEtaPhi;         //! Cluster and tracks Sum Pt in the cone after bkg subtraction, vs eta-phi trigger
  
  TH2F *   fhConeSumPtEtaUESubTrackCell;          //! Cluster and tracks Sum Pt in the cone after bkg subtraction, vs pT trigger
  TH2F *   fhConeSumPtPhiUESubTrackCell;          //! Cluster and tracks Sum Pt in the cone after bkg subtraction, vs pT trigger
  TH2F *   fhConeSumPtEtaUESubTrackCellTrigEtaPhi;//! Cluster and tracks Sum Pt in the cone after bkg subtraction, vs eta-phi trigger
  TH2F *   fhConeSumPtPhiUESubTrackCellTrigEtaPhi;//! Cluster and tracks Sum Pt in the cone after bkg subtraction, vs eta-phi trigger
  
  TH2F *   fhConeSumPtEtaUESubCluster;            //! Cluster Sum Pt in the cone after bkg subtraction, vs pT trigger
  TH2F *   fhConeSumPtPhiUESubCluster;            //! Cluster Sum Pt in the cone after bkg subtraction, vs pT trigger
  TH2F *   fhConeSumPtEtaUESubClusterTrigEtaPhi;  //! Cluster Sum Pt in the cone after bkg subtraction, vs eta-phi trigger
  TH2F *   fhConeSumPtPhiUESubClusterTrigEtaPhi;  //! Cluster Sum Pt in the cone after bkg subtraction, vs eta-phi trigger

  TH2F *   fhConeSumPtEtaUESubCell;               //! Cell Sum amplitude in the cone after bkg subtraction, vs pT trigger
  TH2F *   fhConeSumPtPhiUESubCell;               //! Cell Sum amplitude in the cone after bkg subtraction, vs pT trigger
  TH2F *   fhConeSumPtEtaUESubCellTrigEtaPhi;     //! Cell Sum amplitude in the cone after bkg subtraction, vs eta-phi trigger
  TH2F *   fhConeSumPtPhiUESubCellTrigEtaPhi;     //! Cell Sum amplitude in the cone after bkg subtraction, vs eta-phi trigger
  
  TH2F *   fhConeSumPtEtaUESubTrack;              //! Track Sum Pt in the cone after bkg subtraction, vs pT trigger
  TH2F *   fhConeSumPtPhiUESubTrack;              //! Track Sum Pt in the cone after bkg subtraction, vs pT trigger
  TH2F *   fhConeSumPtEtaUESubTrackTrigEtaPhi;    //! Track Sum Pt in the cone after bkg subtraction, vs eta-phi trigger
  TH2F *   fhConeSumPtPhiUESubTrackTrigEtaPhi;    //! Track Sum Pt in the cone after bkg subtraction, vs eta-phi trigger
  
  TH2F *   fhFractionTrackOutConeEta;             //! Fraction of cone out of tracks acceptance in eta
  TH2F *   fhFractionTrackOutConeEtaTrigEtaPhi;   //! Fraction of cone out of tracks acceptance in eta, vs trigger eta-phi
  TH2F *   fhFractionClusterOutConeEta;           //! Fraction of cone out of clusters acceptance in eta
  TH2F *   fhFractionClusterOutConeEtaTrigEtaPhi; //! Fraction of cone out of clusters acceptance in eta, vs trigger eta-phi
  TH2F *   fhFractionClusterOutConePhi;           //! Fraction of cone out of clusters acceptance in phi
  TH2F *   fhFractionClusterOutConePhiTrigEtaPhi; //! Fraction of cone out of clusters acceptance in phi, vs trigger eta-phi
  
  TH2F *   fhFractionCellOutConeEta;              //! Fraction of cone out of cells acceptance in eta
  TH2F *   fhFractionCellOutConeEtaTrigEtaPhi;    //! Fraction of cone out of cells acceptance in eta, vs trigger eta-phi
  TH2F *   fhFractionCellOutConePhi;              //! Fraction of cone out of cells acceptance in phi
  TH2F *   fhFractionCellOutConePhiTrigEtaPhi;    //! Fraction of cone out of cells acceptance in phi, vs trigger eta-phi
  
  TH2F *   fhConeSumPtClustervsTrack ;            //! Cluster vs tracks Sum Pt Sum Pt in the cone
  TH2F *   fhConeSumPtClusterTrackFrac ;          //! Cluster / tracks Sum Pt Sum Pt in the cone
  TH2F *   fhConeSumPtEtaUESubClustervsTrack ;    //! Cluster vs tracks Sum Pt Sum Pt in the cone, after subtraction in eta band
  TH2F *   fhConeSumPtPhiUESubClustervsTrack ;    //! Cluster vs tracks Sum Pt Sum Pt in the cone, after subtraction in phi band
  TH2F *   fhConeSumPtCellvsTrack;                //! Cell vs tracks Sum Pt Sum Pt in the cone
  TH2F *   fhConeSumPtEtaUESubCellvsTrack ;       //! Cell vs tracks Sum Pt Sum Pt in the cone, after subtraction in eta band
  TH2F *   fhConeSumPtPhiUESubCellvsTrack ;       //! Cell vs tracks Sum Pt Sum Pt in the cone, after subtraction in phi band

  TH2F *   fhEtaBandClustervsTrack ;              //! Accumulated pT in Eta band to estimate UE in cone, clusters vs tracks
  TH2F *   fhPhiBandClustervsTrack ;              //! Accumulated pT in Phi band to estimate UE in cone, clusters vs tracks
  TH2F *   fhEtaBandNormClustervsTrack ;          //! Accumulated pT in Eta band to estimate UE in cone, normalized to cone size, clusters vs tracks
  TH2F *   fhPhiBandNormClustervsTrack ;          //! Accumulated pT in Phi band to estimate UE in cone, normalized to cone size, clusters vs tracks
  TH2F *   fhEtaBandCellvsTrack ;                 //! Accumulated pT in Eta band to estimate UE in cone, cells vs tracks
  TH2F *   fhPhiBandCellvsTrack ;                 //! Accumulated pT in Phi band to estimate UE in cone, cells vs tracks
  TH2F *   fhEtaBandNormCellvsTrack ;             //! Accumulated pT cell in Eta band to estimate UE in cone, normalized to cone size, clusters vs tracks
  TH2F *   fhPhiBandNormCellvsTrack ;             //! Accumulated pT cell in Phi band to estimate UE in cone, normalized to cone

  TH2F *   fhConeSumPtSubvsConeSumPtTotPhiTrack;       //! Tracks, phi band: sum pT in cone after bkg sub vs sum pT in cone before bkg sub
  TH2F *   fhConeSumPtSubNormvsConeSumPtTotPhiTrack;   //! Tracks, phi band: sum pT in cone after bkg sub normalized by sum pT in cone before bkg sub vs sum pT in cone before bkg sub
  TH2F *   fhConeSumPtSubvsConeSumPtTotEtaTrack;       //! Tracks, eta band: sum pT in cone after bkg sub vs sum pT in cone before bkg sub
  TH2F *   fhConeSumPtSubNormvsConeSumPtTotEtaTrack;   //! Tracks, eta band: sum pT in cone after bkg sub normalized by sum pT in cone before bkg sub vs sum pT in cone before bkg sub
  TH2F *   fhConeSumPtSubvsConeSumPtTotPhiCluster;     //! Clusters, phi band: sum pT in cone after bkg sub vs sum pT in cone before bkg sub
  TH2F *   fhConeSumPtSubNormvsConeSumPtTotPhiCluster; //! Clusters, phi band: sum pT in cone after bkg sub normalized by sum pT in cone before bkg sub vs sum pT in cone before bkg sub
  TH2F *   fhConeSumPtSubvsConeSumPtTotEtaCluster;     //! Clusters, eta band: sum pT in cone after bkg sub vs sum pT in cone before bkg sub
  TH2F *   fhConeSumPtSubNormvsConeSumPtTotEtaCluster; //! Clusters, eta band: sum pT in cone after bkg sub normalized by sum pT in cone before bkg sub vs sum pT in cone before bkg sub
  TH2F *   fhConeSumPtSubvsConeSumPtTotPhiCell;        //! Cells, phi band: sum pT in cone after bkg sub vs sum pT in cone before bkg sub
  TH2F *   fhConeSumPtSubNormvsConeSumPtTotPhiCell;    //! Cells, phi band: sum pT in cone after bkg sub normalized by sum pT in cone before bkg sub vs sum pT in cone before bkg sub
  TH2F *   fhConeSumPtSubvsConeSumPtTotEtaCell;        //! Cells, eta band: sum pT in cone after bkg sub vs sum pT in cone before bkg sub
  TH2F *   fhConeSumPtSubNormvsConeSumPtTotEtaCell;    //! Cells, eta band: sum pT in cone after bkg sub normalized by sum pT in cone before bkg sub vs sum pT in cone before bkg sub
  TH2F *   fhConeSumPtVSUETracksEtaBand;               //! fhConeSumPtVSUETracksEtaBand
  TH2F *   fhConeSumPtVSUETracksPhiBand;               //! fhConeSumPtVSUETracksPhiBand
  TH2F *   fhConeSumPtVSUEClusterEtaBand;              //! fhConeSumPtVSUEClusterEtaBand
  TH2F *   fhConeSumPtVSUEClusterPhiBand;              //! fhConeSumPtVSUEClusterPhiBand
  
  //MC
  //
  TH2F *   fhEtaPrimMC  [fgkNmcPrimTypes];        //! Pt vs Eta of generated photon
  TH2F *   fhPhiPrimMC  [fgkNmcPrimTypes];        //! Pt vs Phi of generated photon
  TH1F *   fhEPrimMC    [fgkNmcPrimTypes];        //! Number of generated photon vs E
  TH1F *   fhPtPrimMC   [fgkNmcPrimTypes];        //! Number of generated photon vs pT
  TH1F *   fhPtPrimMCiso[fgkNmcPrimTypes];        //! Number of generated isolated photon vs pT
  
  TH1F *   fhPtPrimMCPi0DecayPairOutOfCone;       //! Pi0 decay photons, with decay pair out of isolation cone
  TH1F *   fhPtPrimMCPi0DecayPairOutOfAcceptance; //! Pi0 decay photons, with decay pair out of detector acceptance
  TH1F *   fhPtPrimMCPi0DecayPairAcceptInConeLowPt;//! Pi0 decay photons, with decay pair in cone and acceptance and lower pT than threshold
  TH1F *   fhPtPrimMCPi0DecayPairAcceptInConeLowPtNoOverlap; //! Pi0 decay photons, with decay pair in cone and acceptance and lower pT than threshold, and do not overlap
  TH1F *   fhPtPrimMCPi0DecayPairNoOverlap;        //! Pi0 decay photons, not overlapped decay

  TH1F *   fhPtPrimMCPi0DecayIsoPairOutOfCone;       //! Pi0 decay photons, with decay pair out of isolation cone, isolated
  TH1F *   fhPtPrimMCPi0DecayIsoPairOutOfAcceptance; //! Pi0 decay photons, with decay pair out of detector acceptance, isolated
  TH1F *   fhPtPrimMCPi0DecayIsoPairAcceptInConeLowPt;//! Pi0 decay photons, with decay pair in cone and acceptance and lower pT than threshold, isolated
  TH1F *   fhPtPrimMCPi0DecayIsoPairAcceptInConeLowPtNoOverlap; //! Pi0 decay photons, with decay pair in cone and acceptance and lower pT than threshold, and do not overlap, isolated
  TH1F *   fhPtPrimMCPi0DecayIsoPairNoOverlap;    //! Pi0 decay photons isolated, not overlapped decay

  TH1F *   fhPtPrimMCPi0Overlap;                  //! Pi0 with overlapped decay photons
  TH1F *   fhPtPrimMCPi0IsoOverlap;               //! Pi0 isolated with overlapped decay photons

  
  TH1F *   fhPtNoIsoMC  [fgkNmcTypes];            //! Number of not isolated mcTypes particle
  TH1F *   fhPtIsoMC    [fgkNmcTypes];            //! Number of isolated mcTypes particle
  TH2F *   fhPhiIsoMC   [fgkNmcTypes];            //! Phi of isolated mcTypes particle
  TH2F *   fhEtaIsoMC   [fgkNmcTypes];            //! eta of isolated mcTypes particle
  
  TH1F *   fhPtDecayIsoMC  [4][fgkNmcTypes] ;     //! Number of isolated Pi0 decay particles (invariant mass tag) for a mcTypes particle
  TH1F *   fhPtDecayNoIsoMC[4][fgkNmcTypes] ;     //! Number of not isolated Pi0 decay particles (invariant mass tag) for a mcTypes particle

  TH2F *   fhPtLambda0MC   [fgkNmcTypes][2];      //! Shower shape of (non) isolated candidates originated by mcTypes particle (do not apply SS cut previously)
 
  // Multiple cut analysis
  TH2F *   fhSumPtLeadingPt[5] ;                  //! Sum Pt in the cone
  TH2F *   fhPtLeadingPt[5] ;                     //! Particle Pt in the cone
  TH2F *   fhPerpSumPtLeadingPt[5] ;              //! Sum Pt in the cone at the perpendicular phi region to trigger axis  (phi +90)
  TH2F *   fhPerpPtLeadingPt[5];                  //! Sum Pt in the cone at the perpendicular phi region to trigger axis  (phi +90)

  TH1F *   fhPtThresIsolated[5][5] ;              //! Isolated particle with pt threshold 
  TH1F *   fhPtFracIsolated[5][5] ;               //! Isolated particle with pt threshold frac
  TH1F *   fhSumPtIsolated[5][5] ;                //! Isolated particle with threshold on cone pt sum
  
  TH2F *   fhEtaPhiPtThresIso[5][5] ;             //! eta vs phi of isolated particles with pt threshold
  TH2F *   fhEtaPhiPtThresDecayIso[5][5] ;        //! eta vs phi of isolated particles with pt threshold, only for decay bit fDecayBits[0]
  TH1F *   fhPtPtThresDecayIso[5][5] ;            //! Number of isolated Pi0 decay particles (invariant mass tag) with pt threshold,, only for decay bit fDecayBits[0]
  
  TH2F *   fhEtaPhiPtFracIso[5][5] ;              //! eta vs phi of isolated particles with pt frac
  TH2F *   fhEtaPhiPtFracDecayIso[5][5] ;         //! eta vs phi of isolated particles with pt frac,, only for decay bit fDecayBits[0]
  TH1F *   fhPtPtFracDecayIso[5][5] ;             //! Number of isolated Pi0 decay particles (invariant mass tag) with pt fra, only for decay bit fDecayBits[0]

  TH2F *   fhEtaPhiPtSumIso[5][5] ;               //! eta vs phi of isolated particles with pt sum
  TH2F *   fhEtaPhiPtSumDecayIso[5][5] ;          //! eta vs phi of isolated particles with pt sum,, only for decay bit fDecayBits[0]
  TH1F *   fhPtPtSumDecayIso[5][5] ;              //! Number of isolated Pi0 decay particles (invariant mass tag) with pt sum, only for decay bit fDecayBits[0]
  
  TH2F *   fhEtaPhiSumDensityIso[5][5];           //! Isolated particle with threshold on cone sum density
  TH2F *   fhEtaPhiSumDensityDecayIso[5][5];      //! Isolated particle with threshold on cone sum density, only for decay bit fDecayBits[0]
  TH1F *   fhPtSumDensityIso[5][5];               //! Isolated particle with threshold on cone sum density
  TH1F *   fhPtSumDensityDecayIso[5][5];          //! Isolated decay particle with threshold on cone sum density, only for decay bit fDecayBits[0]
  
  TH1F *   fhPtFracPtSumIso[5][5] ;               //! Number of isolated Pi0 decay particles (invariant mass tag) with pt sum
  TH1F *   fhPtFracPtSumDecayIso[5][5] ;          //! Number of isolated Pi0 decay particles (invariant mass tag) with pt sum, only for decay bit fDecayBits[0]
  TH2F *   fhEtaPhiFracPtSumIso[5][5];            //! Isolated particle with threshold on cone sum density
  TH2F *   fhEtaPhiFracPtSumDecayIso[5][5];       //! Isolated particle with threshold on cone sum density, only for decay bit fDecayBits[0]
 
  // Multiple cut MC
  TH1F *   fhPtThresIsolatedMC[fgkNmcTypes][5][5];//! Isolated mcTypes particle with pt threshold
  TH1F *   fhPtFracIsolatedMC [fgkNmcTypes][5][5];//! Isolated mcTypes particle with pt frac
  TH1F *   fhSumPtIsolatedMC  [fgkNmcTypes][5][5];//! Isolated mcTypes particle with threshold on cone pt sum
  TH2F *   fhSumPtLeadingPtMC [fgkNmcTypes][5];   //! mcTypes particle for sum Pt, different cone

  
  // Track matching studies
  TH2F *   fhTrackMatchedDEta[2]     ;            //! Eta distance between track and cluster vs cluster E
  TH2F *   fhTrackMatchedDPhi[2]     ;            //! Phi distance between track and cluster vs cluster E
  TH2F *   fhTrackMatchedDEtaDPhi[2] ;            //! Eta vs Phi distance between track and cluster, E cluster > 0.5 GeV
  TH2F *   fhdEdx[2]  ;                           //! matched track dEdx vs cluster E 
  TH2F *   fhEOverP[2];                           //! matched track E cluster over P track vs cluster E, after dEdx cut 
  TH2F *   fhTrackMatchedMCParticle[2];           //! Trace origin of matched particle

  // Shower Shape histograms
  TH2F *   fhELambda0[2];                         //! Shower shape of (non) isolated photons (do not apply SS cut previously)  
  TH2F *   fhPtLambda0[2];                        //! Shower shape of (non) isolated photons (do not apply SS cut previously)
  TH2F *   fhELambda1[2];                         //! Shower shape of (non) isolated photons (do not apply SS cut previously)
  TH2F *   fhELambda0TRD[2];                      //! Shower shape of (non) isolated photons, SM behind TRD (do not apply SS cut previously)
  TH2F *   fhPtLambda0TRD[2];                     //! Shower shape of (non) isolated photons, SM behind TRD (do not apply SS cut previously)
  TH2F *   fhELambda1TRD[2];                      //! Shower shape of (non) isolated photons, SM behind TRD (do not apply SS cut previously)

  TH2F **  fhPtLeadConeBinLambda0 ;               //![fNBkgBin] Candidate shower shape distribution depending on bin of cone leading particle
  TH2F **  fhSumPtConeBinLambda0  ;               //![fNBkgBin] Candidate shower shape distribution depending on bin of cone sum pt
  TH2F **  fhPtLeadConeBinLambda0MC ;             //![fNBkgBin*fgkNmcTypes] Candidate shower shape distribution depending on bin of cone leading particle, per MC particle
  TH2F **  fhSumPtConeBinLambda0MC  ;             //![fNBkgBin*fgkNmcTypes] Candidate shower shape distribution depending on bin of cone sum pt, per MC particle
  
  // Local maxima
  TH2F *   fhNLocMax[2];                          //! number of maxima in selected clusters
  TH2F *   fhELambda0LocMax1[2] ;                 //! E vs lambda0 of selected cluster, 1 local maxima in cluster 
  TH2F *   fhELambda1LocMax1[2] ;                 //! E vs lambda1 of selected cluster, 1 local maxima in cluster 
  TH2F *   fhELambda0LocMax2[2] ;                 //! E vs lambda0 of selected cluster, 2 local maxima in cluster 
  TH2F *   fhELambda1LocMax2[2] ;                 //! E vs lambda1 of selected cluster, 2 local maxima in cluster
  TH2F *   fhELambda0LocMaxN[2] ;                 //! E vs lambda0 of selected cluster, N>2 local maxima in cluster 
  TH2F *   fhELambda1LocMaxN[2] ;                 //! E vs lambda1 of selected cluster, N>2 local maxima in cluster 
  
  // Pile-up
  TH1F *   fhEIsoPileUp[7] ;                      //! Number of isolated particles
  TH1F *   fhPtIsoPileUp[7] ;                     //! Number of isolated particles
  TH1F *   fhENoIsoPileUp[7] ;                    //! Number of not isolated particles
  TH1F *   fhPtNoIsoPileUp[7] ;                   //! Number of not isolated particles
  TH2F *   fhTimeENoCut;                          //! time of cluster vs E, no cut 
  TH2F *   fhTimeESPD;                            //! time of cluster vs E, IsSPDPileUp
  TH2F *   fhTimeESPDMulti;                       //! time of cluster vs E, IsSPDPileUpMulti
  TH2F *   fhTimeNPileUpVertSPD;                  //! time of cluster vs n pile-up vertices from SPD
  TH2F *   fhTimeNPileUpVertTrack;                //! time of cluster vs n pile-up vertices from Tracks
  TH2F *   fhTimeNPileUpVertContributors;         //! time of cluster vs n pile-up vertex from SPD contributors
  TH2F *   fhTimePileUpMainVertexZDistance;       //! time of cluster vs difference of z main vertex and pile-up vertex 
  TH2F *   fhTimePileUpMainVertexZDiamond;        //! time of cluster vs difference of z diamond and pile-up vertex 
  
  AliAnaParticleIsolation(              const AliAnaParticleIsolation & iso) ; // cpy ctor
  AliAnaParticleIsolation & operator = (const AliAnaParticleIsolation & iso) ; // cpy assignment
  
  ClassDef(AliAnaParticleIsolation,29)
} ;


#endif //ALIANAPARTICLEISOLATION_H



