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
    
  void         CalculateCaloUEBand (AliAODPWG4ParticleCorrelation * pCandidate,
                                    Float_t & etaBand, Float_t & phiBand) ;
  void         CalculateTrackUEBand(AliAODPWG4ParticleCorrelation * pCandidate,
                                    Float_t & etaBand, Float_t & phiBand) ;
  
  void         CalculateCaloSignalInCone  (AliAODPWG4ParticleCorrelation * aodParticle, Float_t & coneptsumCluster) ;
  void         CalculateTrackSignalInCone (AliAODPWG4ParticleCorrelation * aodParticle, Float_t & coneptsumTrack  ) ;
  
  Float_t      CalculateExcessAreaFraction(const Float_t excess, const Float_t conesize);
  
  void         CalculateNormalizeUEBandPerUnitArea(AliAODPWG4ParticleCorrelation * pCandidate,
                                                   const Float_t coneptsumCluster, const Float_t coneptsumTrack ) ;
  
  TObjString * GetAnalysisCuts() ;
  
  TList      * GetCreateOutputObjects() ;
  
  void         Init() ;

  void         InitParameters() ;
  
  void         MakeAnalysisFillAOD()  ;
  
  void         MakeAnalysisFillHistograms() ; 
  
  
  void         Print( const Option_t * opt ) const ;
 
  //Analysis specific methods 
  
  void         FillPileUpHistograms(Int_t clusterID) ;
  
  void         FillTrackMatchingShowerShapeControlHistograms(AliAODPWG4ParticleCorrelation  * pCandidate,
                                                             const AliCaloTrackReader * reader,
                                                             const AliCaloPID * pid) ;
  
  void         MakeSeveralICAnalysis( AliAODPWG4ParticleCorrelation * ph ) ; 
  
  // Analysis Setters and Getters
  
  TString      GetCalorimeter()                const { return fCalorimeter       ; }
  Int_t        GetNCones()                     const { return fNCones            ; }
  Int_t        GetNPtThresFrac()               const { return fNPtThresFrac      ; }
  Float_t      GetConeSizes(Int_t i)           const { return fConeSizes[i]      ; }
  Float_t      GetPtThresholds(Int_t i)        const { return fPtThresholds[i]   ; }
  Float_t      GetSumPtThresholds(Int_t i)     const { return fSumPtThresholds[i]; }
  Float_t      GetPtFractions(Int_t i)         const { return fPtFractions[i]    ; }
  
  void         SetCalorimeter(TString & det)         { fCalorimeter     = det    ; }
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

  //Histogrammes setters and getters
  
  void         SetHistoPtSumRangeAndNBins(Float_t min, Float_t max, Int_t n)       {
    fHistoNPtSumBins = n ;    fHistoPtSumMax = max ;    fHistoPtSumMin = min ;     }
  
  Int_t        GetHistoNPtSumBins()            const { return fHistoNPtSumBins   ; }
  Float_t      GetHistoPtSumMin()              const { return fHistoPtSumMin     ; }
  Float_t      GetHistoPtSumMax()              const { return fHistoPtSumMax     ; }
  
  void         SetHistoPtInConeRangeAndNBins(Float_t min, Float_t max, Int_t n)    {
    fHistoNPtInConeBins = n ; fHistoPtInConeMax = max ; fHistoPtInConeMin = min  ; }
  
  Int_t        GetHistoNPtInConeBins()         const { return fHistoNPtInConeBins; }
  Float_t      GetHistoPtInConeMin()           const { return fHistoPtInConeMin  ; }
  Float_t      GetHistoPtInConeMax()           const { return fHistoPtInConeMax  ; }
  
 private:
  
  TString  fCalorimeter ;                         // Calorimeter where neutral particles in cone for isolation are;
  Bool_t   fReMakeIC ;                            // Do isolation analysis
  Bool_t   fMakeSeveralIC ;                       // Do analysis for different IC
  Bool_t   fFillPileUpHistograms;                 // Fill pile-up related histograms
  Bool_t   fFillTMHisto;                          // Fill track matching plots
  Bool_t   fFillSSHisto;                          // Fill Shower shape plots

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
  TH1F *   fhPtDecayIso ;                         //! Number of isolated Pi0 decay particles (invariant mass tag)
  TH1F *   fhPtDecayNoIso ;                       //! Number of not isolated Pi0 decay leading particles (invariant mass tag)
  TH2F *   fhEtaPhiDecayIso ;                     //! eta vs phi of isolated Pi0 decay particles
  TH2F *   fhEtaPhiDecayNoIso ;                   //! eta vs phi of not isolated leading Pi0 decay particles

  TH2F *   fhPtInCone ;                           //! Cluster/track Pt in the cone
  TH2F *   fhPtClusterInCone ;                    //! Cluster Pt in the cone
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
  
  TH2F *   fhEtaBandCluster ;                     //! Accumulated pT in Eta band to estimate UE in cone, only clusters
  TH2F *   fhPhiBandCluster ;                     //! Accumulated pT in Phi band to estimate UE in cone, only clusters
  TH2F *   fhEtaBandTrack   ;                     //! Accumulated pT in Eta band to estimate UE in cone, only tracks
  TH2F *   fhPhiBandTrack   ;                     //! Accumulated pT in Phi band to estimate UE in cone, only tracks
  
  TH2F *   fhConeSumPt ;                          //! Cluster and tracks Sum Pt Sum Pt in the cone
  TH2F *   fhConeSumPtCluster ;                   //! Clusters Sum Pt Sum Pt in the cone
  TH2F *   fhConeSumPtTrack ;                     //! Tracks Sum Pt Sum Pt in the cone
  TH2F *   fhConeSumPtEtaBandUECluster;           //! Cluster Sum Pt in the eta band for clusters, before normalization
  TH2F *   fhConeSumPtPhiBandUECluster;           //! Cluster Sum Pt in the phi band for clusters, before normalization
  TH2F *   fhConeSumPtEtaBandUETrack;             //! Track Sum Pt in the eta band for tracks  , before normalization
  TH2F *   fhConeSumPtPhiBandUETrack;             //! Track Sum Pt in the phi badn for tracks  , before normalization
  
  TH2F *   fhConeSumPtTrigEtaPhi ;                //! Cluster and tracks Sum Pt Sum Pt in the cone, per eta-phi bin of trigger,
  TH2F *   fhConeSumPtEtaBandUEClusterTrigEtaPhi; //! Cluster Sum Pt in the eta band for clusters, per eta-phi bin of trigger,before normalization
  TH2F *   fhConeSumPtPhiBandUEClusterTrigEtaPhi; //! Cluster Sum Pt in the phi band for clusters, per eta-phi bin of trigger, before normalization
  TH2F *   fhConeSumPtEtaBandUETrackTrigEtaPhi;   //! Track Sum Pt in the eta band for tracks  , per eta-phi bin of trigger, before normalization
  TH2F *   fhConeSumPtPhiBandUETrackTrigEtaPhi;   //! Track Sum Pt in the phi badn for tracks  , per eta-phi bin of trigger, before normalization

  TH2F *   fhConeSumPtEtaUESub;                   //! Cluster and tracks Sum Pt in the cone after bkg subtraction, vs pT trigger
  TH2F *   fhConeSumPtPhiUESub;                   //! Cluster and tracks Sum Pt in the cone after bkg subtraction, vs pT trigger
  TH2F *   fhConeSumPtEtaUESubTrigEtaPhi;         //! Cluster and tracks Sum Pt in the cone after bkg subtraction, vs eta-phi trigger
  TH2F *   fhConeSumPtPhiUESubTrigEtaPhi;         //! Cluster and tracks Sum Pt in the cone after bkg subtraction, vs eta-phi trigger
    
  TH2F *   fhConeSumPtEtaUESubCluster;            //! Cluster Sum Pt in the cone after bkg subtraction, vs pT trigger
  TH2F *   fhConeSumPtPhiUESubCluster;            //! Cluster Sum Pt in the cone after bkg subtraction, vs pT trigger
  TH2F *   fhConeSumPtEtaUESubClusterTrigEtaPhi;  //! Cluster Sum Pt in the cone after bkg subtraction, vs eta-phi trigger
  TH2F *   fhConeSumPtPhiUESubClusterTrigEtaPhi;  //! Cluster Sum Pt in the cone after bkg subtraction, vs eta-phi trigger

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
  
  TH2F *   fhConeSumPtClustervsTrack ;            //! Cluster vs tracks Sum Pt Sum Pt in the cone
  TH2F *   fhConeSumPtEtaUESubClustervsTrack ;    //! Cluster vs tracks Sum Pt Sum Pt in the cone, after subtraction in eta band
  TH2F *   fhConeSumPtPhiUESubClustervsTrack ;    //! Cluster vs tracks Sum Pt Sum Pt in the cone, after subtraction in phi band
  TH2F *   fhEtaBandClustervsTrack ;              //! Accumulated pT in Eta band to estimate UE in cone, clusters vs tracks
  TH2F *   fhPhiBandClustervsTrack ;              //! Accumulated pT in Phi band to estimate UE in cone, clusters vs tracks
  TH2F *   fhEtaBandNormClustervsTrack ;          //! Accumulated pT in Eta band to estimate UE in cone, normalized to cone size, clusters vs tracks
  TH2F *   fhPhiBandNormClustervsTrack ;          //! Accumulated pT in Phi band to estimate UE in cone, normalized to cone size, clusters vs tracks
  
  //MC
  TH1F *   fhPtIsoPrompt;                         //! Number of isolated prompt gamma 
  TH2F *   fhPhiIsoPrompt;                        //! Phi of isolated prompt gamma
  TH2F *   fhEtaIsoPrompt;                        //! eta of isolated prompt gamma
  TH1F *   fhPtThresIsolatedPrompt[5][5];         //! Isolated prompt gamma with pt threshold 
  TH1F *   fhPtFracIsolatedPrompt[5][5];          //! Isolated prompt gamma with pt frac
  TH2F *   fhPtSumIsolatedPrompt[5];              //! Isolated prompt gamma with threshold on cone pt sume
  TH1F *   fhPtIsoFragmentation;                  //! Number of isolated fragmentation gamma 
  TH2F *   fhPhiIsoFragmentation;                 //! Phi of isolated fragmentation gamma
  TH2F *   fhEtaIsoFragmentation;                 //! eta of isolated fragmentation gamma
  TH1F *   fhPtThresIsolatedFragmentation[5][5];  //! Isolated fragmentation gamma with pt threshold 
  TH1F *   fhPtFracIsolatedFragmentation[5][5];   //! Isolated fragmentation gamma with pt frac
  TH2F *   fhPtSumIsolatedFragmentation[5];       //! Isolated fragmentation gamma with threshold on cone pt sume
  TH1F *   fhPtIsoPi0;                            //! Number of isolated pi0 (2 gamma)
  TH2F *   fhPhiIsoPi0;                           //! Phi of isolated pi0 (2 gamma)
  TH2F *   fhEtaIsoPi0;                           //! eta of isolated pi0 (2 gamma)
  TH1F *   fhPtThresIsolatedPi0[5][5];            //! Isolated pi0 (2 gamma) with pt threshold 
  TH1F *   fhPtFracIsolatedPi0[5][5];             //! Isolated pi0 (2 gamma) with pt frac
  TH2F *   fhPtSumIsolatedPi0[5];                 //! Isolated pi0 (2 gamma) with threshold on cone pt sum
  TH1F *   fhPtIsoPi0Decay;                       //! Number of isolated pi0 decay gamma 
  TH2F *   fhPhiIsoPi0Decay;                      //! Phi of isolated pi0 decay gamma
  TH2F *   fhEtaIsoPi0Decay;                      //! eta of isolated pi0 decay gamma
  TH1F *   fhPtThresIsolatedPi0Decay[5][5];       //! Isolated pi0 decay gamma with pt threshold 
  TH1F *   fhPtFracIsolatedPi0Decay[5][5];        //! Isolated pi0 decay gamma with pt frac
  TH2F *   fhPtSumIsolatedPi0Decay[5];            //! Isolated pi0 decay gamma with threshold on cone pt sume
  TH1F *   fhPtIsoEtaDecay;                       //! Number of isolated eta decay gamma 
  TH2F *   fhPhiIsoEtaDecay;                      //! Phi of isolated eta decay gamma
  TH2F *   fhEtaIsoEtaDecay;                      //! eta of isolated eta decay gamma
  TH1F *   fhPtThresIsolatedEtaDecay[5][5];       //! Isolated eta decay gamma with pt threshold 
  TH1F *   fhPtFracIsolatedEtaDecay[5][5];        //! Isolated eta decay gamma with pt frac
  TH2F *   fhPtSumIsolatedEtaDecay[5];            //! Isolated eta fecay gamma with threshold on cone pt sume  
  TH1F *   fhPtIsoOtherDecay;                     //! Number of isolated other decay gamma 
  TH2F *   fhPhiIsoOtherDecay;                    //! Phi of isolated other decay gamma
  TH2F *   fhEtaIsoOtherDecay;                    //! eta of isolated other decay gamma
  TH1F *   fhPtThresIsolatedOtherDecay[5][5];     //! Isolated OtherDecay gamma with pt threshold 
  TH1F *   fhPtFracIsolatedOtherDecay[5][5];      //! Isolated OtherDecay gamma with pt frac
  TH2F *   fhPtSumIsolatedOtherDecay[5];          //! Isolated OtherDecay gamma with threshold on cone pt sume	
//  TH1F *   fhPtIsoConversion;                     //! Number of isolated Conversion gamma 
//  TH2F *   fhPhiIsoConversion;                    //! Phi of isolated Conversion gamma
//  TH2F *   fhEtaIsoConversion;                    //! eta of isolated Conversion gamma
//  TH1F *   fhPtThresIsolatedConversion[5][5];     //! Isolated Conversion gamma with pt threshold 
//  TH1F *   fhPtFracIsolatedConversion[5][5];      //! Isolated Conversion gamma with pt frac
//  TH2F *   fhPtSumIsolatedConversion[5];          //! Isolated Conversion gamma with threshold on cone pt sume
  TH1F *   fhPtIsoHadron;                        //! Number of isolated Hadron 
  TH2F *   fhPhiIsoHadron;                       //! Phi of isolated Hadron
  TH2F *   fhEtaIsoHadron;                       //! eta of isolated Hadron
  TH1F *   fhPtThresIsolatedHadron[5][5];        //! Isolated Hadron gamma with pt threshold 
  TH1F *   fhPtFracIsolatedHadron[5][5];         //! Isolated Hadron gamma with pt frac
  TH2F *   fhPtSumIsolatedHadron[5];             //! Isolated Hadron gamma with threshold on cone pt sume

  // Multi Cut analysis Several IC
  TH1F *   fhPtNoIsoPi0;                          //! Number of not isolated leading pi0 (2 gamma)
  TH1F *   fhPtNoIsoPi0Decay;                     //! Number of not isolated leading pi0 decay gamma 
  TH1F *   fhPtNoIsoEtaDecay;                     //! Number of not isolated leading eta decay gamma 
  TH1F *   fhPtNoIsoOtherDecay;                   //! Number of not isolated leading other decay gamma 
  TH1F *   fhPtNoIsoPrompt;                       //! Number of not isolated leading prompt gamma 
  TH1F *   fhPtIsoMCPhoton;                       //! Number of isolated leading gamma 
  TH1F *   fhPtNoIsoMCPhoton;                     //! Number of not isolated leading gamma 
//  TH1F *   fhPtNoIsoConversion;                   //! Number of not isolated leading conversion gamma 
  TH1F *   fhPtNoIsoFragmentation;                //! Number of not isolated leading fragmentation gamma 
  TH1F *   fhPtNoIsoHadron;                      //! Number of not isolated leading hadrons 
  
  TH2F *   fhSumPtLeadingPt[5] ;                  //! Sum Pt in the cone
  TH2F *   fhPtLeadingPt[5] ;                     //! Particle Pt in the cone
  TH2F *   fhPerpSumPtLeadingPt[5] ;              //! Sum Pt in the cone at the perpendicular phi region to trigger axis  (phi +90)
  TH2F *   fhPerpPtLeadingPt[5];                  //! Sum Pt in the cone at the perpendicular phi region to trigger axis  (phi +90)

  TH1F *   fhPtThresIsolated[5][5] ;              //! Isolated particle with pt threshold 
  TH1F *   fhPtFracIsolated[5][5] ;               //! Isolated particle with pt threshold frac
  TH1F *   fhPtSumIsolated[5][5] ;                //! Isolated particle with threshold on cone pt sum
  
  TH2F *   fhEtaPhiPtThresIso[5][5] ;             //! eta vs phi of isolated particles with pt threshold
  TH2F *   fhEtaPhiPtThresDecayIso[5][5] ;        //! eta vs phi of isolated particles with pt threshold
  TH1F *   fhPtPtThresDecayIso[5][5] ;            //! Number of isolated Pi0 decay particles (invariant mass tag) with pt threshold   
  
  TH2F *   fhEtaPhiPtFracIso[5][5] ;              //! eta vs phi of isolated particles with pt frac
  TH2F *   fhEtaPhiPtFracDecayIso[5][5] ;         //! eta vs phi of isolated particles with pt frac
  TH1F *   fhPtPtFracDecayIso[5][5] ;             //! Number of isolated Pi0 decay particles (invariant mass tag) with pt fra

  TH2F *   fhEtaPhiPtSumIso[5][5] ;               //! eta vs phi of isolated particles with pt sum
  TH2F *   fhEtaPhiPtSumDecayIso[5][5] ;          //! eta vs phi of isolated particles with pt sum
  TH1F *   fhPtPtSumDecayIso[5][5] ;              //! Number of isolated Pi0 decay particles (invariant mass tag) with pt sum
  
  TH2F *   fhEtaPhiSumDensityIso[5][5];           //! Isolated particle with threshold on cone sum density
  TH2F *   fhEtaPhiSumDensityDecayIso[5][5];      //! Isolated particle with threshold on cone sum density
  TH1F *   fhPtSumDensityIso[5][5];               //! Isolated particle with threshold on cone sum density
  TH1F *   fhPtSumDensityDecayIso[5][5];          //! Isolated decay particle with threshold on cone sum density
  
  TH1F *   fhPtFracPtSumIso[5][5] ;               //! Number of isolated Pi0 decay particles (invariant mass tag) with pt sum
  TH1F *   fhPtFracPtSumDecayIso[5][5] ;          //! Number of isolated Pi0 decay particles (invariant mass tag) with pt sum
  TH2F *   fhEtaPhiFracPtSumIso[5][5];            //! Isolated particle with threshold on cone sum density
  TH2F *   fhEtaPhiFracPtSumDecayIso[5][5];       //! Isolated particle with threshold on cone sum density
 
  // Track matching studies
  TH2F *   fhTrackMatchedDEta[2]     ;            //! Eta distance between track and cluster vs cluster E
  TH2F *   fhTrackMatchedDPhi[2]     ;            //! Phi distance between track and cluster vs cluster E
  TH2F *   fhTrackMatchedDEtaDPhi[2] ;            //! Eta vs Phi distance between track and cluster, E cluster > 0.5 GeV
  TH2F *   fhdEdx[2]  ;                           //! matched track dEdx vs cluster E 
  TH2F *   fhEOverP[2];                           //! matched track E cluster over P track vs cluster E, after dEdx cut 
  TH2F *   fhTrackMatchedMCParticle[2];           //! Trace origin of matched particle

  // Shower Shape histograms
  TH2F *   fhELambda0[2];                         //! Shower shape of (non) isolated photons (do not apply SS cut previously)
  TH2F *   fhELambda1[2];                         //! Shower shape of (non) isolated photons (do not apply SS cut previously)
  TH2F *   fhELambda0SSBkg;                       //! Shower shape of non isolated photons close to isolation threshold (do not apply SS cut previously)
  TH2F *   fhELambda0TRD[2];                      //! Shower shape of (non) isolated photons, SM behind TRD (do not apply SS cut previously)
  TH2F *   fhELambda1TRD[2];                      //! Shower shape of (non) isolated photons, SM behind TRD (do not apply SS cut previously)
  TH2F *   fhELambda0MCPhoton[2];                 //! Shower shape of (non) isolated photon candidates originated by particle photon not decay (do not apply SS cut previously)
  TH2F *   fhELambda0MCPi0[2];                    //! Shower shape of (non) isolated photon candidates originated by particle 2 merged pi0 photons (do not apply SS cut previously)
  TH2F *   fhELambda0MCPi0Decay[2];               //! Shower shape of (non) isolated photon candidates originated by particle pi0 decay photon (do not apply SS cut previously)
  TH2F *   fhELambda0MCEtaDecay[2];               //! Shower shape of (non) isolated photon candidates originated by particle eta decay photon (do not apply SS cut previously)
  TH2F *   fhELambda0MCOtherDecay[2];             //! Shower shape of (non) isolated photon candidates originated by particle other decay photon (do not apply SS cut previously)
  TH2F *   fhELambda0MCHadron[2];                 //! Shower shape of (non) isolated photon candidates originated by particle other hadrons (do not apply SS cut previously)

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
  
  //Histograms settings
  Int_t    fHistoNPtSumBins;                      // Number of bins in PtSum histograms
  Float_t  fHistoPtSumMax;                        // PtSum maximum in histogram
  Float_t  fHistoPtSumMin;	                  // PtSum minimum in histogram
  Int_t    fHistoNPtInConeBins;                   // Number of bins in PtInCone histogram
  Float_t  fHistoPtInConeMax;                     // PtInCone maximum in histogram
  Float_t  fHistoPtInConeMin;                     // PtInCone maximum in histogram 

  AliAnaParticleIsolation(              const AliAnaParticleIsolation & iso) ; // cpy ctor
  AliAnaParticleIsolation & operator = (const AliAnaParticleIsolation & iso) ; // cpy assignment
  
  ClassDef(AliAnaParticleIsolation,23)
} ;


#endif //ALIANAPARTICLEISOLATION_H



