#ifndef ALIANAPARTICLEJETFINDERCORRELATION_H
#define ALIANAPARTICLEJETFINDERCORRELATION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliAnaParticleJetLeadingConeCorrelation
/// \ingroup CaloTrackCorrelationsAnalysis 
/// \brief Correlate trigger particle and reconstructed jet
///
/// Class that contains the algorithm for the analysis of particle (direct gamma) - jet
/// (standard jet found with JETAN) correlation
/// Particle and jet for correlation found by independent algorithms.
/// For Example direct isolated photon found in AliAnaGammaDirect and the jet with JETAN
///
/// More information can be found in this [twiki](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations)
/// and particularly in this [section](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations#AliAnaParticleJetFinderCorrelation).
///
/// \author Adam Matyja <Adam.Matyja@cern.ch>, INP-PAN-Krakow, main developper
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS, first implementation
///_________________________________________________________________________

// --- ROOT system ---
class TH2F;
class TTree;
class TRandom2;

//---- Analysis system ----
#include "AliAnaCaloTrackCorrBaseClass.h"

class AliAnaParticleJetFinderCorrelation : public AliAnaCaloTrackCorrBaseClass {
       
 public:
    
           AliAnaParticleJetFinderCorrelation() ;
    
  virtual ~AliAnaParticleJetFinderCorrelation() ;

  // General methods
  
  void     InitParameters();
  
  TList *  GetCreateOutputObjects();

  void     MakeAnalysisFillAOD() ;
  
  void     MakeAnalysisFillHistograms() ;
  
  // To access non standard branch
  Int_t    SelectJet(AliAODPWG4Particle * particle, TClonesArray * aodRecJets) ;
    
  void     Print(const Option_t * opt)           const;
  
  // Settings
  
  Bool_t   OnlyIsolated()                        const { return fSelectIsolated              ; }
  void     SelectIsolated(Bool_t select)               { fSelectIsolated = select            ; }
  
  Float_t  GetConeSize()                         const { return fConeSize                    ; }
  Float_t  GetPtThresholdInCone()                const { return fPtThresholdInCone           ; }	   
  Double_t GetDeltaPhiMaxCut()                   const { return fDeltaPhiMaxCut              ; }
  Double_t GetDeltaPhiMinCut()                   const { return fDeltaPhiMinCut              ; }
  Double_t GetRatioMaxCut()                      const { return fRatioMaxCut                 ; }
  Double_t GetRatioMinCut()                      const { return fRatioMinCut                 ; }  	   
  Bool_t   AreJetRefTracks()                     const { return fUseJetRefTracks             ; }
  Bool_t   IsCorrelationMadeInHistoMaker()       const { return fMakeCorrelationInHistoMaker ; } 
  Double_t GetJetConeSize()                      const { return fJetConeSize                 ; } 
  Double_t GetJetMinPt()                         const { return fJetMinPt                    ; }
  Double_t GetJetMinPtBkgSub()                   const { return fJetMinPtBkgSub              ; }
  Double_t GetJetAreaFraction()                  const { return fJetAreaFraction             ; }
  Double_t GetGammaConeSize()                    const { return fGammaConeSize               ; }

  void     SetConeSize(Float_t cone)                   { fConeSize = cone                    ; }
  void     SetPtThresholdInCone(Float_t pt)            { fPtThresholdInCone = pt             ; }
    
  void     SetDeltaPhiCutRange(Double_t phimin, Double_t phimax)
            { fDeltaPhiMaxCut = phimax  ;  fDeltaPhiMinCut = phimin                          ; }
    
  void     SetRatioCutRange(Double_t ratiomin, Double_t ratiomax)
            { fRatioMaxCut    = ratiomax;  fRatioMinCut    = ratiomin                        ; }
    
  void     UseJetRefTracks(Bool_t use)                 { fUseJetRefTracks = use              ; }	
  void     SetMakeCorrelationInHistoMaker(Bool_t make) { fMakeCorrelationInHistoMaker = make ; }
  void     SetJetConeSize(Double_t cone)               { fJetConeSize     = cone             ; }
  void     SetJetMinPt(Double_t minpt)                 { fJetMinPt        = minpt            ; }
  void     SetJetMinPtBkgSub(Double_t minpt)           { fJetMinPtBkgSub  = minpt            ; }
  void     SetJetAreaFraction(Double_t areafr)         { fJetAreaFraction = areafr           ; }
  void     SetGammaConeSize(Float_t cone)              { fGammaConeSize   = cone             ; }

  // Settings for non standard jet branch
  TString  GetJetBranchName()                    const { return fJetBranchName               ; }
  void     SetJetBranchName(const char *name)          { fJetBranchName = name               ; }
//void     SwitchOnNonStandardJetFromReader()          { fNonStandardJetFromReader = kTRUE   ; }
//void     SwitchOffNonStandardJetFromReader()         { fNonStandardJetFromReader = kFALSE  ; }
//Bool_t   IsNonStandardJetFromReader()                { return fNonStandardJetFromReader    ; }

  TString  GetBkgJetBranchName()                 const { return fBkgJetBranchName            ; }
  void     SetBkgJetBranchName(const char *name)       { fBkgJetBranchName = name            ; }
  void     SwitchOnBackgroundJetFromReader()           { fBackgroundJetFromReader = kTRUE    ; }
  void     SwitchOffBackgroundJetFromReader()          { fBackgroundJetFromReader = kFALSE   ; }
  Bool_t   IsBackgroundJetFromReader()                 { return fBackgroundJetFromReader     ; }

  //switches for photons
  void     SwitchOnBackgroundSubtractionGamma()        { fUseBackgroundSubtractionGamma = kTRUE ; }
  void     SwitchOffBackgroundSubtractionGamma()       { fUseBackgroundSubtractionGamma = kFALSE; }
  Bool_t   IsBackgroundSubtractionGamma()        const { return fUseBackgroundSubtractionGamma  ; }

  void     CalculateBkg(TVector3 gamma, TVector3 jet,Double_t *vector,Int_t type);

  void     FindMCgenInfo();//gives information on generated level

  void     SwitchOnSaveGJTree()                         { fSaveGJTree = kTRUE  ; }
  void     SwitchOffSaveGJTree()                        { fSaveGJTree = kFALSE ; }
  Bool_t   IsSaveGJTree()                         const { return fSaveGJTree   ; }

  void     SwitchOnMostEnergetic()                      { fMostEnergetic = kTRUE  ; fMostOpposite = kFALSE ; }
  void     SwitchOffMostEnergetic()                     { fMostEnergetic = kFALSE ; fMostOpposite = kTRUE  ; }
  void     SwitchOffMostOpposite()                      { fMostEnergetic = kTRUE  ; fMostOpposite = kFALSE ; }
  void     SwitchOnMostOpposite()                       { fMostEnergetic = kFALSE ; fMostOpposite = kTRUE  ; }
  Bool_t   IsMostEnergetic()                      const { return fMostEnergetic ; }
  Bool_t   IsMostOpposite()                       const { return fMostOpposite  ; }

  //switches for histograms
  void     SwitchOnHistogramJetBkg()                    { fUseHistogramJetBkg = kTRUE  ; }
  void     SwitchOffHistogramJetBkg()                   { fUseHistogramJetBkg = kFALSE ; }
  Bool_t   IsHistogramJetBkg()                    const { return fUseHistogramJetBkg   ; }

  void     SwitchOnHistogramTracks()                    { fUseHistogramTracks = kTRUE  ; }
  void     SwitchOffHistogramTracks()                   { fUseHistogramTracks = kFALSE ; }
  Bool_t   IsHistogramTracks()                    const { return fUseHistogramTracks   ; }

  void     SwitchOnHistogramJetTracks()                 { fUseHistogramJetTracks = kTRUE  ; }
  void     SwitchOffHistogramJetTracks()                { fUseHistogramJetTracks = kFALSE ; }
  Bool_t   IsHistogramJetTracks()                 const { return fUseHistogramJetTracks   ; }

  void     SwitchOnMCStudies()                          { fMCStudies = kTRUE  ; }
  void     SwitchOffMCStudies()                         { fMCStudies = kFALSE ; }
  Bool_t   IsMCStudies()                          const { return fMCStudies   ; }

private:

  // selection parameters
  Double_t   fDeltaPhiMaxCut ;                    ///<  Minimum Delta Phi Gamma-Leading
  Double_t   fDeltaPhiMinCut ;                    ///<  Maximum Delta Phi Gamma-Leading
  Double_t   fRatioMaxCut ;                       ///<  Jet/particle Ratio cut maximum
  Double_t   fRatioMinCut ;                       ///<  Jet/particle Ratio cut minimum
  
  Double_t   fConeSize  ;                         ///<  Jet cone size to calculate fragmentation function
  Double_t   fPtThresholdInCone ;                 ///<  Jet pT threshold in jet cone
  Bool_t     fUseJetRefTracks ;                   ///<  Use track references from JETAN not the AOD tracks to calculate fragmentation function
  Bool_t     fMakeCorrelationInHistoMaker ;       ///< Make particle-jet correlation in histogram maker
  Bool_t     fSelectIsolated ;                    ///<  Select only trigger particles isolated

  Double_t   fJetConeSize ;                       ///<  Reconstructed jet cone size
  Double_t   fJetMinPt ;                          ///<  Minumum jet pt, default 5GeV/c
  Double_t   fJetMinPtBkgSub ;                    ///<  Minumum jet pt after bkg subtraction, default -100 GeV/c
  Double_t   fJetAreaFraction ;                   ///<  Jet area fraction X in X*pi*R^2, default 0.6
//Bool_t     fNonStandardJetFromReader;           ///<  use non standard jet from reader //new
  TString    fJetBranchName ;                     ///<  name of jet branch not set in reader part //new
  Bool_t     fBackgroundJetFromReader;            ///<  use background jet from reader //new
  TString    fBkgJetBranchName ;                  ///<  name of background jet branch not set in reader part //new

  Double_t   fGammaConeSize ;                     ///<  Isolation cone radius
  Bool_t     fUseBackgroundSubtractionGamma;      ///<  flag to use backgrouind subtraction for photons or not
  Bool_t     fSaveGJTree;                         ///<  flag to save gamma-jet tree
  Bool_t     fMostEnergetic;                      ///<  flag to choose gamma-jet pairs most energetic
  Bool_t     fMostOpposite;                       ///<  flag to choose gamma-jet pairs most opposite

  Bool_t     fUseHistogramJetBkg;                 ///<  flag to save bkg jet histograms
  Bool_t     fUseHistogramTracks;                 ///<  flag to save CTS tracks features
  Bool_t     fUseHistogramJetTracks;              ///<  flag to save jet tracks features

  Bool_t     fMCStudies;                          ///<  flag to use MC methods

  TRandom2 * fGenerator;                          //!<! pointer to random generator object

  TLorentzVector fMomentum;                       //!<! momentum
  
  // Histograms
  TH2F *     fhDeltaEta;                          //!<!  Difference of jet eta and trigger particle eta as function of trigger particle pT
//TH2F *     fhDeltaPhi;                          //!<!  Difference of jet phi and trigger particle phi as function of trigger particle pT
  TH2F *     fhDeltaPhiCorrect;                   //!<!  Difference of jet phi and trigger particle phi as function of trigger particle pT
  TH2F *     fhDeltaPhi0PiCorrect;                //!<!  Difference of jet phi and trigger particle phi as function of trigger particle pT

  TH2F *     fhDeltaPt;                           //!<!  Difference of jet pT and trigger particle pT as function of trigger particle pT
  TH2F *     fhPtRatio;                           //!<!  Ratio of jet pT and trigger particle pT as function of trigger particle pT
  TH2F *     fhPt;                                //!<!  jet pT vs trigger particle pT
  
  TH2F *     fhFFz ;                              //!<!  Accepted reconstructed jet fragmentation function, z=pt^particle,jet/pttrig
  TH2F *     fhFFxi;                              //!<!  Accepted reconstructed jet fragmentation function, xsi = ln(pttrig/pt^particle,jet)
  TH2F *     fhFFpt;                              //!<!  Jet particle pt distribution in cone
  TH2F *     fhNTracksInCone;                     //!<!  jet multiplicity in cone
  
  TH2F *     fhJetFFz ;                           //!<!  Accepted reconstructed jet fragmentation function, z=pt^particle,jet/ptjet
  TH2F *     fhJetFFxi;                           //!<!  Accepted reconstructed jet fragmentation function, xsi = ln(ptjet/pt^particle,jet)
  TH2F *     fhJetFFpt;                           //!<!  Jet particle pt distribution in jet cone
  TH2F *     fhJetFFzCor ;                        //!<!  Accepted reconstructed jet fragmentation function, z=pt^particle,jet*-cos(jet,trig)/ptjet
  TH2F *     fhJetFFxiCor;                        //!<!  Accepted reconstructed jet fragmentation function, xsi = ln(ptjet/pt^particle*-cos(jet,trig),jet)

  TH1F * fhGamPtPerTrig ;                         //!<!  per trigger normalisation
  TH2F * fhPtGamPtJet ;                           //!<!  gamma jet correlation filling

  // background from RC
  TH2F *     fhBkgFFz[5] ;                        //!<!  Background fragmentation function, z=ptjet/pttrig
  TH2F *     fhBkgFFxi[5];                        //!<!  Background fragmentation function, xsi = ln(pttrig/ptjet)
  TH2F *     fhBkgFFpt[5];                        //!<!  Background particle pt distribution in cone
  TH2F *     fhBkgNTracksInCone[5];               //!<!  Background multiplicity in cone
  TH2F *     fhBkgSumPtInCone[5];                 //!<!  Background sum pt in cone
  TH2F *     fhBkgSumPtnTracksInCone[5];          //!<!  Background sum pt over ntracks in cone

  // temporary histograms
  TH2F * fhNjetsNgammas;                          //!<!  Number of jets vs number of photons in the event
  TH1F * fhCuts;                                  //!<!  Number of events after cuts

  TH2F * fhDeltaEtaBefore;                        //!<!  Difference of jet eta and trigger particle eta as function of trigger particle pT
  TH2F * fhDeltaPhiBefore;                        //!<!  Difference of jet phi and trigger particle phi as function of trigger particle pT
  TH2F * fhDeltaPtBefore;                         //!<!  Difference of jet pT and trigger particle pT as function of trigger particle pT
  TH2F * fhPtRatioBefore;                         //!<!  Ratio of jet pT and trigger particle pT as function of trigger particle pT
  TH2F * fhPtBefore;                              //!<!  jet pT vs trigger particle pT
  TH2F * fhDeltaPhi0PiCorrectBefore;              //!<!  Difference of jet phi and trigger particle phi (0,pi) as function of trigger particle pT

  // temporary jet histograms
  TH1F * fhJetPtBefore;                           //!<!  Pt of all jets
  TH1F * fhJetPtBeforeCut;                        //!<!  Pt of all jets after bkg correction, raw jet pt>fJetMinPt
  TH1F * fhJetPt;                                 //!<!  Pt of all jets after bkg correction
  TH1F * fhJetPtMostEne;                          //!<!  Pt of the most energetic jet
  TH1F * fhJetPhi;	                              //!<!  Phi of all jets
  TH1F * fhJetEta;	                              //!<!  Eta of all jets
  TH2F * fhJetEtaVsPt;	                          //!<!  Eta of all jets vs pt
  TH2F * fhJetPhiVsEta;	                          //!<!  Phi vs eta of all jets
  TH2F * fhJetEtaVsNpartInJet;                    //!<!  Eta vs number of particles in jet for all jets
  TH2F * fhJetEtaVsNpartInJetBkg;                 //!<!  Eta vs number of particles in jet for background subtracted jets
  TH2F * fhJetChBkgEnergyVsPt;                    //!<!  background energy of each charged jet vs jet pt
  TH2F * fhJetChAreaVsPt;                         //!<!  area of each charged jet vs jet pt
  TH2F * fhTrackPhiVsEta;                         //!<!  Phi vs eta of all chosen tracks in all events
  TH1F * fhTrackAveTrackPt;                       //!<!  average track pt in event
  TH1F * fhJetNjetOverPtCut[10];                  //!<!  number of reconstructed jets in event over pT threshold
  TH2F * fhJetChBkgEnergyVsArea;                  //!<!  area of each charged jet vs jet background
  TH2F * fhJetRhoVsPt;                            //!<!  jet energy density vs jet pt
  TH2F * fhJetRhoVsCentrality;                    //!<!  jet energy density vs centrality
  TH1F * fhJetNparticlesInJet;                    //!<!  number of particles in jets
  TH2F * fhJetDeltaEtaDeltaPhi;                   //!<!  delta eta vs delta phi for (jet-track) <-0.8,0.8>
  TH2F * fhJetDeltaEtaDeltaPhiAllTracks;          //!<!  delta eta vs delta phi for (jet-track) <-pi,pi>

  TH1F * fhJetAveTrackPt;                         //!<!  average track from jets pt in event
  TH2F * fhJetNtracksInJetAboveThr[6];            //!<!  number of tracks in jet with pt above 0,1,2,3,4,5GeV
  TH2F * fhJetRatioNTrkAboveToNTrk[5];            //!<!  ratio tracks in jet with pt above 1,2,3,4,5GeV to ntracks
  TH2F * fhJetNtrackRatioMostEne  [5];            //!<!  the same for most energetic jet
  TH2F * fhJetNtrackRatioJet5GeV  [5];            //!<!  the same for pt jet above 5 GeV
  TH2F * fhJetNtrackRatioLead5GeV [5];            //!<!  the same for jet with leading particle pt>5GeV

  // temporary background jet histograms
  TH1F * fhBkgJetBackground[4];                   //!<!  background from jet bkg branch
  TH1F * fhBkgJetSigma[4];                        //!<!  sigma of jet in backgroud branch
  TH1F * fhBkgJetArea[4];                         //!<!  area of jet in bkg branch

  // temporary photon histograms
  TH1F * fhPhotonPtMostEne;                       //!<!  most pt photon
  TH1F * fhPhotonAverageEnergy;                   //!<!  average energy of photon
  TH1F * fhPhotonRatioAveEneToMostEne;            //!<!  ratio average energy to most energetic photon
  TH1F * fhPhotonAverageEnergyMinus1;             //!<!  average energy of photon w/o most ene photon
  TH1F * fhPhotonRatioAveEneMinus1ToMostEne;      //!<!  ratio average energy of photon w/o most ene photon to most energetic photon
  TH1F * fhPhotonNgammaMoreAverageToNgamma;       //!<!  number of gammas with ene. more than average ene divided by no. of gammas
  TH1F * fhPhotonNgammaMoreAverageMinus1ToNgamma; //!<!  number of gammas with ene. more than average ene (w/o most ene gamma) divided by no. of gammas
  TH1F * fhPhotonNgammaOverPtCut[10];             //!<!  number of photons in event over pT threshold
  TH2F * fhPhotonBkgRhoVsNtracks;                 //!<!  average energy in one cell vs n tracks
  TH2F * fhPhotonBkgRhoVsNclusters;               //!<!  average energy in one cell vs n clusters
  TH2F * fhPhotonBkgRhoVsCentrality;              //!<!  average energy in one cell vs centrality
  TH2F * fhPhotonBkgRhoVsNcells;                  //!<!  average energy in one cell vs n cells
  TH1F * fhPhotonPt;                              //!<!  pt of gamma before bkg correction
  TH1F * fhPhotonPtCorrected;                     //!<!  pt of gamma after background correction
  TH1F * fhPhotonPtCorrectedZoom;                 //!<!  pt of gamma after background correction in +-5 GeV/c
  TH1F * fhPhotonPtDiff;                          //!<!  bkg correction = n_cells * median_rho
  TH2F * fhPhotonPtDiffVsCentrality;              //!<!  correction vs centrality
  TH2F * fhPhotonPtDiffVsNcells;                  //!<!  correction vs Ncells
  TH2F * fhPhotonPtDiffVsNtracks;                 //!<!  correction vs Ntracks
  TH2F * fhPhotonPtDiffVsNclusters;               //!<!  correction vs Nclustres

  TH1F * fhPhotonSumPtInCone;                     //!<!  sum pt in cone before correction
  TH1F * fhPhotonSumPtCorrectInCone;              //!<!  sum pt in cone afrer correction
  TH1F * fhPhotonSumPtChargedInCone;              //!<!  sum pt of charged tracks in the cone before correction


  // temporary jet histograms after selection
  TH2F * fhSelectedJetPhiVsEta;	                  //!<!  phi vs eta of selected jet
  TH2F * fhSelectedJetChBkgEnergyVsPtJet;         //!<!  background energy of selected charged jet vs jet pt
  TH2F * fhSelectedJetChAreaVsPtJet;              //!<!  area of selected charged jet vs jet pt
  TH1F * fhSelectedJetNjet;                       //!<!  number of jets in selected event
  TH1F * fhSelectedNtracks;                       //!<!  number of tracks in selected event
  TH2F * fhSelectedTrackPhiVsEta;                 //!<!  Phi vs eta of all chosen tracks in selected events

  TH1F * fhCuts2;                                 //!<!  efficienct cuts

  // temporary photon histogram after selection
  TH2F * fhSelectedPhotonNLMVsPt;                 //!<!  nlm vs pt for selected photons
  TH2F * fhSelectedPhotonLambda0VsPt;             //!<!  lambda0 vs pt for selected photons
  TH2F * fhRandomPhiEta[5];                       //!<!  eta and phi from random generator

  // MC generated histograms
  TH1F * fhMCPhotonCuts;                          //!<!  generated photon cuts
  TH1F * fhMCPhotonPt;                            //!<!  generated direct photon pt
  TH2F * fhMCPhotonEtaPhi;                        //!<!  generated direct photon eta vs phi
  TH1F * fhMCJetOrigin;                           //!<!  generated origin of jet
  TH2F * fhMCJetNPartVsPt;                        //!<!  generated N parts vs pt full jet
  TH2F * fhMCJetChNPartVsPt;                      //!<!  generated N parts vs pt charged jet
  TH2F * fhMCJetNPart150VsPt;                     //!<!  generated N parts (pt>150 MeV/c) vs pt full jet
  TH2F * fhMCJetChNPart150VsPt;                   //!<!  generated N parts (pt>150 MeV/c) vs pt charged jet
  TH2F * fhMCJetChNPart150ConeVsPt;               //!<!  generated N parts (pt>150 MeV/c) vs pt charged jet R=0.4
  TH1F * fhMCJetRatioChFull;                      //!<!  generated ratio pt charged/full jet
  TH1F * fhMCJetRatioCh150Ch;                     //!<!  generated ratio pt charged(pt>150MeV/c)/charged jet
  TH2F * fhMCJetEtaPhi;                           //!<!  generated jet eta vs phi for full jet
  TH2F * fhMCJetChEtaPhi;                         //!<!  generated jet eta vs phi for charged jet
  TH2F * fhMCJet150EtaPhi;                        //!<!  generated jet eta vs phi full jet (pt>150 MeV/c)
  TH2F * fhMCJetCh150EtaPhi;                      //!<!  generated jet eta vs phi charged jet (pt>150 MeV/c)
  TH2F * fhMCJetCh150ConeEtaPhi;                  //!<!  generated jet eta vs phi charged jet (pt>150 MeV/c) R=0.4
 
  // tree with data gamma and jet
  TTree *  fTreeGJ;                               //!<!  gamma-jet tree
  Double_t fGamPt;                                ///<  pt
  Double_t fGamLambda0;                           ///<  lambda 0
  Int_t    fGamNLM;                               ///<  NLM
  Double_t fGamSumPtCh;                           ///<  energy in isolation cone charged
  Double_t fGamTime;                              ///<  time
  Int_t    fGamNcells;                            ///<  ncells
  Double_t fGamEta;                               ///<  eta photon
  Double_t fGamPhi;                               ///<  phi photon
  Double_t fGamSumPtNeu;                          ///<  energy in isolation cone neutral
  Int_t    fGamNtracks;                           ///<  number of tracks in iso cone
  Int_t    fGamNclusters;                         ///<  number of clusters in iso cone
  Double_t fGamAvEne;                             ///<  average energy of photons (without most ene)
  Double_t fJetPhi;                               ///<  jet phi
  Double_t fJetEta;                               ///<  eta phi
  Double_t fJetPt;                                ///<  jet pt
  Double_t fJetBkgChEne;                          ///<  bkg energy of jet
  Double_t fJetArea;                              ///<  jet area
  Int_t    fJetNtracks;                           ///<  number of jet tracks
  Int_t    fJetNtracks1;                          ///<  number of jet tracks with pt>1 GeV/c
  Int_t    fJetNtracks2;                          ///<  number of jet tracks with pt>2 GeV/c
  Double_t fJetRho;                               ///<  jet rho in event
  Int_t    fEventNumber;                          ///<  event number
  Int_t    fNtracks;                              ///<  n tracks in event
  Double_t fZvertex;                              ///<  z vertex
  Double_t fCentrality;                           ///<  centrality
  Bool_t   fIso;                                  ///<  flag isolated or not
  Double_t fGamRho;                               ///<  background energy for photons per cell in EMCal

  Double_t fMCGamPt;                              ///<  MC gen pt photon
  Double_t fMCGamEta;                             ///<  MC gen eta photon
  Double_t fMCGamPhi;                             ///<  MC gen phi photon
  Int_t    fMCPartonType;                         ///<  MC gen parton type origin of jet
  Double_t fMCJetPt;                              ///<  MC gen full jet pt
  Double_t fMCJetChPt;                            ///<  MC gen charged jet pt
  Double_t fMCJet150Pt;                           ///<  MC gen full jet (pt^particles>150MeV/c) pt
  Double_t fMCJetCh150Pt;                         ///<  MC gen charged jet (pt^particles>150MeV/c) pt
  Int_t    fMCJetNPart;                           ///<  MC gen number of full jet particles
  Int_t    fMCJetChNPart;                         ///<  MC gen number of charged jet particles
  Int_t    fMCJet150NPart;                        ///<  MC gen number of full jet particles (pt>150MeV/c)
  Int_t    fMCJetCh150NPart;                      ///<  MC gen number of charged jet particles (pt>150MeV/c)
  Double_t fMCJetEta;                             ///<  MC gen full jet eta
  Double_t fMCJetPhi;                             ///<  MC gen full jet phi
  Double_t fMCJetChEta;                           ///<  MC gen charged jet eta
  Double_t fMCJetChPhi;                           ///<  MC gen charged jet phi
  Double_t fMCJet150Eta;                          ///<  MC gen full jet eta (pt>150MeV/c)
  Double_t fMCJet150Phi;                          ///<  MC gen full jet phi (pt>150MeV/c)
  Double_t fMCJetCh150Eta;                        ///<  MC gen charged jet eta (pt>150MeV/c)
  Double_t fMCJetCh150Phi;                        ///<  MC gen charged jet phi (pt>150MeV/c)

  Double_t fMCJetCh150ConePt;                     ///<  MC gen charged jet (pt^particles>150MeV/c),R=0.4 pt
  Int_t    fMCJetCh150ConeNPart;                  ///<  MC gen number of charged jet particles (pt>150MeV/c),R=0.4
  Double_t fMCJetCh150ConeEta;                    ///<  MC gen charged jet eta (pt>150MeV/c),R=0.4
  Double_t fMCJetCh150ConePhi;                    ///<  MC gen charged jet phi (pt>150MeV/c),R=0.4

  /// Copy constructor not implemented.
  AliAnaParticleJetFinderCorrelation(              const AliAnaParticleJetFinderCorrelation & g) ;
    
  /// Assignment operator not implemented.
  AliAnaParticleJetFinderCorrelation & operator = (const AliAnaParticleJetFinderCorrelation & g) ;
  
  /// \cond CLASSIMP
  ClassDef(AliAnaParticleJetFinderCorrelation,3) ;
  /// \endcond

 } ;

#endif //ALIANAPARTICLEJETFINDERCORRELATION_H



