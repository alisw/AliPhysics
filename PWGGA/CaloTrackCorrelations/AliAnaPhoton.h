#ifndef ALIANAPHOTON_H
#define ALIANAPHOTON_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliAnaPhoton
/// \ingroup CaloTrackCorrelationsAnalysis 
/// \brief Filter EMCal/PHOS clusters for photon analysis.
///
/// Class for the photon identification.
/// Clusters from calorimeters are identified/selected as photons (candidates)
/// and kept in AOD format (AliAODPWG4Particle). Few histograms produced.
/// Produces input for other analysis classes like AliAnaPi0,
/// AliAnaParticleHadronCorrelation ...
///
/// More information can be found in this [twiki](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations)
/// and particularly in this [section](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations#AliAnaPhoton).
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
///_________________________________________________________________________

// --- ROOT system ---
class TH2F ;
class TH1F;
class TObjString;
class TList ;

// --- ANALYSIS system ---
#include "AliAnaCaloTrackCorrBaseClass.h"

class AliAnaPhoton : public AliAnaCaloTrackCorrBaseClass {

 public:
    
  AliAnaPhoton() ;
    
  /// Virtual destructor, not implemented.
  virtual ~AliAnaPhoton() { ; }
	
  //---------------------------------------
  // General analysis frame methods
  //---------------------------------------
  
  TObjString * GetAnalysisCuts();
  
  TList      * GetCreateOutputObjects();
  
  void         Init();

  void         InitParameters();

  void         MakeAnalysisFillAOD()  ;

  void         MakeAnalysisFillHistograms() ; 
  
  void         Print(const Option_t * opt)const;
    
  
  // Analysis methods
  
  Bool_t       ClusterSelected(AliVCluster* cl, Int_t nlm) ;
  
  void         FillAcceptanceHistograms();
  
//  void         DistanceToAddedSignalAtGeneratorLevel(Int_t label, Int_t nprim, 
//                                     Float_t photonE, Float_t photonEta, Float_t photonPhi);
  
  void         FillShowerShapeHistograms( AliVCluster* cluster, Int_t mcTag, Int_t nlm,
                                         Float_t maxCellEFraction, Int_t & largeTimeInside) ;
  
  void         SwitchOnFillShowerShapeHistograms()        { fFillSSHistograms      = kTRUE  ; }
  void         SwitchOffFillShowerShapeHistograms()       { fFillSSHistograms      = kFALSE ; }  

  void         SwitchOnFillShowerShapeHistogramsPerNLM()  { fFillSSNLocMaxHisto    = kTRUE  ; }
  void         SwitchOffFillShowerShapeHistogramsPerNLM() { fFillSSNLocMaxHisto    = kFALSE ; }  
  
  void         SwitchOnFillEMCALRegionSSHistograms()      { fFillEMCALRegionSSHistograms = kTRUE  ; }
  void         SwitchOffFillEMCALRegionSSHistograms()     { fFillEMCALRegionSSHistograms = kFALSE ; }  

  void         SwitchOnFillConversionVertexHistograms()   { fFillConversionVertexHisto = kTRUE  ; }
  void         SwitchOffFillConversionVertexHistograms()  { fFillConversionVertexHisto = kFALSE ; }  
  
  void         SwitchOnOnlySimpleSSHistoFill()            { fFillOnlySimpleSSHisto = kTRUE  ; }
  void         SwitchOffOnlySimpleHistoFill()             { fFillOnlySimpleSSHisto = kFALSE ; }
  
  void         SwitchOnFillTrackMultiplicityHistograms()  { fFillTrackMultHistograms = kTRUE  ; }
  void         SwitchOffFillTrackMultiplicityHistograms() { fFillTrackMultHistograms = kFALSE ; }
  
  void         FillTrackMatchingResidualHistograms(AliVCluster* calo, Int_t cut);
  
  void         SwitchOnTMHistoFill()                      { fFillTMHisto           = kTRUE  ; }
  void         SwitchOffTMHistoFill()                     { fFillTMHisto           = kFALSE ; }

  void         SwitchOnTMTrackPtHistoFill()               { fFillTMHistoTrackPt    = kTRUE  ; }
  void         SwitchOffTMTrackPtHistoFill()              { fFillTMHistoTrackPt    = kFALSE ; }
  
  void         FillPileUpHistograms(AliVCluster* cluster, AliVCaloCells *cells, Int_t absIdMax) ;
 
  void         SetConstantTimeShift(Float_t shift)        { fConstantTimeShift     = shift  ; }
  
  // Cocktail generator studies
  void         CocktailGeneratorsClusterOverlaps(AliVCluster* calo, Int_t mctag);
  
  void         ActivityNearCluster(Int_t icalo, Float_t en, Float_t eta, Float_t phi, Int_t mctag, TObjArray *clusterList) ;
  void         SwitchOnStudyClusterLocalActivity()        { fStudyActivityNearCluster = kTRUE  ; }
  void         SwitchOffStudyClusterLocalActivity()       { fStudyActivityNearCluster = kFALSE ; }  
    
  // Analysis parameters setters getters
    
  // ** Cluster selection methods **
  
  void         SetMinDistanceToBadChannel(Float_t m1, Float_t m2, Float_t m3) {
                fMinDist = m1; fMinDist2 = m2; fMinDist3 = m3; }

  void         SetTimeCut(Double_t min, Double_t max) { fTimeCutMin = min; 
                                                        fTimeCutMax = max          ; }
  Double_t     GetTimeCutMin()                  const { return fTimeCutMin         ; }
  Double_t     GetTimeCutMax()                  const { return fTimeCutMax         ; }	
	
  void         SetNCellCut(Int_t n)                   { fNCellsCut = n             ; }
  Double_t     GetNCellCut()                    const { return fNCellsCut          ; }
  
  void         SetNLMCut(Int_t min, Int_t max)        { fNLMCutMin = min; 
    fNLMCutMax = max                ; }
  Int_t        GetNLMCutMin()                   const { return fNLMCutMin          ; }
  Int_t        GetNLMCutMax()                   const { return fNLMCutMax          ; }	
  
  Bool_t       IsTrackMatchRejectionOn()        const { return fRejectTrackMatch   ; }
  void         SwitchOnTrackMatchRejection()          { fRejectTrackMatch = kTRUE  ; }
  void         SwitchOffTrackMatchRejection()         { fRejectTrackMatch = kFALSE ; }  

  void         SwitchOnAcceptanceHistoPerEBin()       { fFillEBinAcceptanceHisto = kTRUE  ; }
  void         SwitchOffAcceptanceHistoPerEBin()      { fFillEBinAcceptanceHisto = kFALSE ; }
  
  void         SetNEBinCuts(Int_t nb)           { fNEBinCuts = nb            ; }
  void         SetEBinCutsAt(Int_t i, Float_t va) { if(i < 15) fEBinCuts[i] = va ; }
  
  void         FillNOriginHistograms(Int_t n)         { fNOriginHistograms = n ; 
    if(n > fgkNmcTypes    ) fNOriginHistograms  = fgkNmcTypes     ; }
  void         FillNPrimaryHistograms(Int_t n)        { fNPrimaryHistograms= n ;
    if(n > fgkNmcPrimTypes) fNPrimaryHistograms = fgkNmcPrimTypes ; }

  /// For MC histograms in arrays, index in the array corresponds to a MC originating particle type
  enum mcTypes    { kmcPhoton     =  0,    kmcPi0Decay = 1,       kmcEtaDecay      = 2,  kmcOtherDecay = 3,
                    kmcPi0        =  4,    kmcEta      = 5,       kmcElectron      = 6,
                    kmcConversion =  7,    kmcOther    = 8,       kmcAntiNeutron   = 9,
                    kmcAntiProton = 10,    kmcNeutron  = 11,      kmcProton        = 12, kmcChPion = 13, 
                    kmcPrompt     = 14,    kmcFragmentation = 15,
                    kmcISR        = 16,    kmcString   = 17  };

  /// Total number of cluster MC origin histograms
  static const Int_t fgkNmcTypes = 18;

  /// For MC histograms in arrays, index in the array corresponds to a MC generated primary particle type
  enum mcPTypes   { kmcPPhoton = 0,       kmcPPi0Decay      = 1,  kmcPEtaDecay = 2,     kmcPOtherDecay = 3,
                    kmcPPrompt = 4,       kmcPFragmentation = 5,  kmcPISR      = 6 };
  
  /// Total number of MC primary histograms
  static const Int_t fgkNmcPrimTypes = 7;
  
  /// For MC histograms with shower shape in arrays, index in the array corresponds to a MC originating particle type
  enum mcssTypes  { kmcssPhoton = 0,      kmcssOther = 1,       kmcssPi0 = 2,
                    kmcssEta = 3,         kmcssConversion = 4,  kmcssElectron = 5  };  
  
  /// Total number of MC histograms for shower shape studies.
  static const Int_t fgkNssTypes = 6 ;

  /// For MC histograms with cocktail generator checks in arrays, index in the array corresponds to a MC originating particle type
  enum mcGenTypes { kmcGenPi0Merged = 1,  kmcGenPi0Decay = 2,   kmcGenEtaDecay = 3,
                    kmcGenPhoton    = 4,  kmcGenElectron = 5,   kmcGenOther    = 6 };  
  
  /// Total number of MC histograms for cocktail generator checks.
  static const Int_t fgkNGenTypes = 7 ;

  
  private:
 
  Float_t  fMinDist ;                               ///<  Minimal distance to bad channel to accept cluster
  Float_t  fMinDist2;                               ///<  Cuts on Minimal distance to study acceptance evaluation
  Float_t  fMinDist3;                               ///<  One more cut on distance used for acceptance-efficiency study
    
  Bool_t   fRejectTrackMatch ;                      ///<  If PID on, reject clusters which have an associated TPC track
    
  Bool_t   fFillTMHisto;                            ///<  Fill track matching plots
  Bool_t   fFillTMHistoTrackPt;                     ///<  Fill track matching plots depending on Track pT
    
  Double_t fTimeCutMin  ;                           ///<  Remove clusters/cells with time smaller than this value, in ns
  Double_t fTimeCutMax  ;                           ///<  Remove clusters/cells with time larger than this value, in ns
    
  Int_t    fNCellsCut ;                             ///<  Accept for the analysis clusters with more than fNCellsCut cells
    
  Int_t    fNLMCutMin  ;                            ///<  Remove clusters/cells with number of local maxima smaller than this value
  Int_t    fNLMCutMax  ;                            ///<  Remove clusters/cells with number of local maxima larger than this value
    
  Bool_t   fFillSSHistograms ;                      ///<  Fill shower shape histograms

  Bool_t   fFillEMCALRegionSSHistograms ;           ///<  Fill shower shape histograms in EMCal slices
    
  Bool_t   fFillConversionVertexHisto   ;           ///<  Fill histograms depending on the conversion vertex
  
  Bool_t   fFillOnlySimpleSSHisto;                  ///<  Fill selected cluster histograms, selected SS histograms
    
  Bool_t   fFillSSNLocMaxHisto;                     ///<  Fill shower shape histograms for different NLM
  
  Bool_t   fFillTrackMultHistograms;             ///<  Fill cluster/photon pT spectrum histograms vs track multiplicity or track sum pt

  Int_t    fNOriginHistograms;                      ///<  Fill only NOriginHistograms of the 14 defined types
  Int_t    fNPrimaryHistograms;                     ///<  Fill only NPrimaryHistograms of the 7 defined types
  
  TLorentzVector fMomentum;                         //!<! Cluster momentum, temporary container
  TLorentzVector fMomentum2;                        //!<! Cluster momentum, temporary container
  TLorentzVector fPrimaryMom;                       //!<! Primary MC momentum, temporary container
  TVector3       fProdVertex;                       //!<! Primary MC production vertex, temporary container
  
  Float_t  fConstantTimeShift;                      ///<  Apply a 600 ns time shift in case of simulation, shift in ns.
  
  Bool_t   fFillEBinAcceptanceHisto;                ///<  Fill histograms with cluster eta-phi distribution and column-row cell, for different energy bins
  Float_t  fEBinCuts[15] ;                          ///<  Energy bins cut 
  Int_t    fNEBinCuts;                              ///<  Number of energy bin cuts

  Bool_t   fStudyActivityNearCluster;               ///<  Activate analysis of multiplicity and energy deposit near the cluster
  
  //
  // Histograms
  //
    
  TH1F * fhClusterCutsE [10];                       //!<! control histogram on the different photon selection cuts, E
  TH1F * fhClusterCutsPt[10];                       //!<! control histogram on the different photon selection cuts, pT
  TH2F * fhNCellsE;                                 //!<! number of cells in cluster vs E
  TH2F * fhCellsE;                                  //!<! energy of cells in cluster vs E of cluster
  TH2F * fhMaxCellDiffClusterE;                     //!<! Fraction of energy carried by cell with maximum energy
  TH2F * fhTimePt;                                  //!<! Time of photon cluster vs pt
  TH2F * fhEtaPhi  ;                                //!<! Pseudorapidity vs Phi of clusters for E > 0.5
  
  TH1F * fhEPhoton    ;                             //!<! Number of identified photon vs energy
  TH1F * fhPtPhoton   ;                             //!<! Number of identified photon vs transerse momentum
  TH2F * fhPhiPhoton  ;                             //!<! Azimuthal angle of identified  photon vs transerse momentum
  TH2F * fhEtaPhoton  ;                             //!<! Pseudorapidity of identified  photon vs transerse momentum
  TH2F * fhEtaPhiPhoton  ;                          //!<! Pseudorapidity vs Phi of identified  photon for E > 0.5
  TH2F * fhEtaPhi05Photon  ;                        //!<! Pseudorapidity vs Phi of identified  photon for E < 0.5

  TH2F * fhPtCentralityPhoton    ;                  //!<! centrality  vs photon pT
  TH2F * fhPtEventPlanePhoton    ;                  //!<! event plane vs photon pT
  
  TH2F * fhPtPhotonNTracks    [10];                 //!<! Track multiplicity distribution per event vs track pT, different pT cuts
  TH2F * fhPtPhotonSumPtTracks[10];                 //!<! Track sum pT distribution per event vs track pT, different pT cuts  

  // Shower shape
  TH2F * fhNLocMax;                                 //!<! number of maxima in selected clusters

  TH2F * fhDispE;                                   //!<! Cluster dispersion vs E
  TH2F * fhDispPt;                                  //!<! Cluster dispersion vs pT
  TH2F * fhLam0E;                                   //!<! Cluster lambda0 vs  E
  TH2F * fhLam0Pt;                                  //!<! Cluster lambda0 vs  pT
  TH2F * fhLam1E;                                   //!<! Cluster lambda1 vs  E
  TH2F * fhLam1Pt;                                  //!<! Cluster lambda1 vs  pT

  TH2F * fhLam0PtNLM1;                              //!<! Cluster lambda0 vs  pT, for clusters with NLM=1 
  TH2F * fhLam0PtNLM2;                              //!<! Cluster lambda0 vs  pT, for clusters with NLM=2 
  TH2F * fhLam1PtNLM1;                              //!<! Cluster lambda0 vs  pT, for clusters with NLM=1 
  TH2F * fhLam1PtNLM2;                              //!<! Cluster lambda0 vs  pT, for clusters with NLM=2 
  
  TH2F * fhDispETRD;                                //!<! Cluster dispersion vs E, SM covered by TRD
  TH2F * fhLam0ETRD;                                //!<! Cluster lambda0 vs  E, SM covered by TRD
  TH2F * fhLam0PtTRD;                               //!<! Cluster lambda0 vs  pT, SM covered by TRD
  TH2F * fhLam1ETRD;                                //!<! Cluster lambda1 vs  E, SM covered by TRD

  TH2F * fhDispETM;                                 //!<! Cluster dispersion vs E, cut on Track Matching residual
  TH2F * fhLam0ETM;                                 //!<! Cluster lambda0 vs  E, cut on Track Matching residual
  TH2F * fhLam0PtTM;                                //!<! Cluster lambda0 vs  pT, cut on Track Matching residual
  TH2F * fhLam1ETM;                                 //!<! Cluster lambda1 vs  E, cut on Track Matching residual
  
  TH2F * fhDispETMTRD;                              //!<! Cluster dispersion vs E, SM covered by TRD, cut on Track Matching residual
  TH2F * fhLam0ETMTRD;                              //!<! Cluster lambda0 vs  E, SM covered by TRD, cut on Track Matching residual
  TH2F * fhLam0PtTMTRD;                             //!<! Cluster lambda0 vs  pT, SM covered by TRD, cut on Track Matching residual
  TH2F * fhLam1ETMTRD;                              //!<! Cluster lambda1 vs  E, SM covered by TRD, cut on Track Matching residual
  
  TH2F * fhNCellsLam0LowE;                          //!<! number of cells in cluster vs lambda0
  TH2F * fhNCellsLam1LowE;                          //!<! number of cells in cluster vs lambda1
  TH2F * fhNCellsDispLowE;                          //!<! number of cells in cluster vs dispersion
  TH2F * fhNCellsLam0HighE;                         //!<! number of cells in cluster vs lambda0, E>2
  TH2F * fhNCellsLam1HighE;                         //!<! number of cells in cluster vs lambda1, E>2
  TH2F * fhNCellsDispHighE;                         //!<! number of cells in cluster vs dispersion, E>2
  
  TH2F * fhEtaLam0LowE;                             //!<! Cluster eta vs lambda0, E<2
  TH2F * fhPhiLam0LowE;                             //!<! Cluster phi vs lambda0, E<2
  TH2F * fhEtaLam0HighE;                            //!<! Cluster eta vs lambda0, E>2
  TH2F * fhPhiLam0HighE;                            //!<! Cluster phi vs lambda0, E>2
  TH2F * fhLam0DispLowE;                            //!<! Cluster lambda0 vs dispersion, E<2
  TH2F * fhLam0DispHighE;                           //!<! Cluster lambda0 vs dispersion, E>2
  TH2F * fhLam1Lam0LowE;                            //!<! Cluster lambda1 vs lambda0, E<2
  TH2F * fhLam1Lam0HighE;                           //!<! Cluster lambda1 vs lambda0, E>2
  TH2F * fhDispLam1LowE;                            //!<! Cluster disp vs lambda1, E<2
  TH2F * fhDispLam1HighE;                           //!<! Cluster disp vs lambda1, E>2
    
  TH2F * fhDispEtaE ;                               //!<! shower dispersion in eta direction
  TH2F * fhDispPhiE ;                               //!<! shower dispersion in phi direction
  TH2F * fhSumEtaE ;                                //!<! shower dispersion in eta direction
  TH2F * fhSumPhiE ;                                //!<! shower dispersion in phi direction
  TH2F * fhSumEtaPhiE ;                             //!<! shower dispersion in eta and phi direction
  TH2F * fhDispEtaPhiDiffE ;                        //!<! shower dispersion eta - phi
  TH2F * fhSphericityE ;                            //!<! shower sphericity in eta vs phi
  TH2F * fhDispSumEtaDiffE ;                        //!<! difference of 2 eta dispersions
  TH2F * fhDispSumPhiDiffE ;                        //!<! difference of 2 phi dispersions
  TH2F * fhDispEtaDispPhi[7] ;                      //!<! shower dispersion in eta direction vs phi direction for 5 E bins [0-2],[2-4],[4-6],[6-10],[> 10]
  TH2F * fhLambda0DispEta[7] ;                      //!<! shower shape correlation l0 vs disp eta
  TH2F * fhLambda0DispPhi[7] ;                      //!<! shower shape correlation l0 vs disp phi
  
  // Fill MC dependent histograms, Origin of this cluster is ...

  TH2F * fhMCDeltaE [fgkNmcTypes];                  //!<! MC-Reco E distribution coming from MC particle
  TH2F * fhMCDeltaPt[fgkNmcTypes];                  //!<! MC-Reco pT distribution coming from MC particle
  TH2F * fhMC2E     [fgkNmcTypes];                  //!<! E distribution, Reco vs MC coming from MC particle
  TH2F * fhMC2Pt    [fgkNmcTypes];                  //!<! pT distribution, Reco vs MC coming from MC particle
  
  TH1F * fhMCE  [fgkNmcTypes];                      //!<! Number of identified photon vs cluster energy coming from MC particle
  TH1F * fhMCPt [fgkNmcTypes];                      //!<! Number of identified photon vs cluster pT     coming from MC particle
  TH2F * fhMCPhi[fgkNmcTypes];                      //!<! Phi of identified photon coming from MC particle
  TH2F * fhMCEta[fgkNmcTypes];                      //!<! eta of identified photon coming from MC particle

  TH1F * fhEPrimMC  [fgkNmcPrimTypes];              //!<! Number of generated photon vs energy
  TH1F * fhPtPrimMC [fgkNmcPrimTypes];              //!<! Number of generated photon vs pT
  TH2F * fhPhiPrimMC[fgkNmcPrimTypes];              //!<! Phi of generted photon
  TH2F * fhYPrimMC  [fgkNmcPrimTypes];              //!<! Rapidity of generated photon
  TH2F * fhEtaPrimMC[fgkNmcPrimTypes];              //!<! Eta of generated photon
  
  TH1F * fhEPrimMCAcc  [fgkNmcPrimTypes];           //!<! Number of generated photon vs energy, in calorimeter acceptance
  TH1F * fhPtPrimMCAcc [fgkNmcPrimTypes];           //!<! Number of generated photon vs pT, in calorimeter acceptance
  TH2F * fhPhiPrimMCAcc[fgkNmcPrimTypes];           //!<! Phi of generted photon, in calorimeter acceptance
  TH2F * fhEtaPrimMCAcc[fgkNmcPrimTypes];           //!<! Phi of generted photon, in calorimeter acceptance
  TH2F * fhYPrimMCAcc  [fgkNmcPrimTypes];           //!<! Rapidity of generated photon, in calorimeter acceptance
  
  // Shower Shape MC
    
  TH2F * fhMCELambda0   [fgkNssTypes] ;             //!<! E vs Lambda0     from MC particle
  TH2F * fhMCPtLambda0  [fgkNssTypes] ;             //!<! pT vs Lambda0     from MC particle
  TH2F * fhMCELambda1   [fgkNssTypes] ;             //!<! E vs Lambda1     from MC particle
  TH2F * fhMCEDispersion[fgkNssTypes] ;             //!<! E vs Dispersion  from MC particle
  
  TH2F * fhMCPtLambda0Overlaps[fgkNssTypes][3] ;          //!<! pT vs Lambda0 from MC particle, check if there were 0, 1 or +1 overlaps
  
//  TH2F * fhMCPhotonELambda0NoOverlap ;              //!<! E vs Lambda0     from MC photons, no overlap
//  TH2F * fhMCPhotonELambda0TwoOverlap ;             //!<! E vs Lambda0     from MC photons, 2 particles overlap
//  TH2F * fhMCPhotonELambda0NOverlap ;               //!<! E vs Lambda0     from MC photons, N particles overlap
  
  TH2F * fhMCLambda0vsClusterMaxCellDiffE0[fgkNssTypes]; //!<! Lambda0 vs fraction of energy of max cell for E < 2 GeV
  TH2F * fhMCLambda0vsClusterMaxCellDiffE2[fgkNssTypes]; //!<! Lambda0 vs fraction of energy of max cell for 2< E < 6 GeV
  TH2F * fhMCLambda0vsClusterMaxCellDiffE6[fgkNssTypes]; //!<! Lambda0 vs fraction of energy of max cell for E > 6 GeV
  TH2F * fhMCNCellsvsClusterMaxCellDiffE0 [fgkNssTypes]; //!<! NCells  vs fraction of energy of max cell for E < 2
  TH2F * fhMCNCellsvsClusterMaxCellDiffE2 [fgkNssTypes]; //!<! NCells  vs fraction of energy of max cell for 2 < E < 6 GeV
  TH2F * fhMCNCellsvsClusterMaxCellDiffE6 [fgkNssTypes]; //!<! NCells  vs fraction of energy of max cell for E > 6
  TH2F * fhMCNCellsE            [fgkNssTypes];           //!<! NCells per cluster vs energy
  TH2F * fhMCMaxCellDiffClusterE[fgkNssTypes];           //!<! Fraction of energy carried by cell with maximum energy

  TH2F * fhMCEDispEta       [fgkNssTypes] ;              //!<! shower dispersion in eta direction
  TH2F * fhMCEDispPhi       [fgkNssTypes] ;              //!<! shower dispersion in phi direction
  TH2F * fhMCESumEtaPhi     [fgkNssTypes] ;              //!<! shower dispersion in eta vs phi direction
  TH2F * fhMCEDispEtaPhiDiff[fgkNssTypes] ;              //!<! shower dispersion in eta -phi direction
  TH2F * fhMCESphericity    [fgkNssTypes] ;              //!<! shower sphericity, eta vs phi
  TH2F * fhMCDispEtaDispPhi [7][fgkNssTypes] ;           //!<! shower dispersion in eta direction vs phi direction for 5 E bins [0-2],[2-4],[4-6],[6-10],[> 10]
  TH2F * fhMCLambda0DispEta [7][fgkNssTypes] ;           //!<! shower shape correlation l0 vs disp eta
  TH2F * fhMCLambda0DispPhi [7][fgkNssTypes] ;           //!<! shower shape correlation l0 vs disp phi

  //Embedding
  TH2F * fhEmbeddedSignalFractionEnergy ;           //!<! Fraction of photon energy of embedded signal vs cluster energy
  
  TH2F * fhEmbedPhotonELambda0FullSignal ;          //!<!  Lambda0 vs E for embedded photons with more than 90% of the cluster energy
  TH2F * fhEmbedPhotonELambda0MostlySignal ;        //!<!  Lambda0 vs E for embedded photons with 90%<fraction<50%
  TH2F * fhEmbedPhotonELambda0MostlyBkg ;           //!<!  Lambda0 vs E for embedded photons with 50%<fraction<10%
  TH2F * fhEmbedPhotonELambda0FullBkg ;             //!<!  Lambda0 vs E for embedded photons with less than 10% of the cluster energy
  
  TH2F * fhEmbedPi0ELambda0FullSignal ;             //!<!  Lambda0 vs E for embedded photons with more than 90% of the cluster energy
  TH2F * fhEmbedPi0ELambda0MostlySignal ;           //!<!  Lambda0 vs E for embedded photons with 90%<fraction<50%
  TH2F * fhEmbedPi0ELambda0MostlyBkg ;              //!<!  Lambda0 vs E for embedded photons with 50%<fraction<10%
  TH2F * fhEmbedPi0ELambda0FullBkg ;                //!<!  Lambda0 vs E for embedded photons with less than 10% of the cluster energy
  
  // Track Matching
    
//TH2F * fhTrackMatchedDEta[2]           ;          //!<! Eta distance between track and cluster vs cluster E, after and before photon cuts
//TH2F * fhTrackMatchedDPhi[2]           ;          //!<! Phi distance between track and cluster vs cluster E, after and before photon cuts
//TH2F * fhTrackMatchedDEtaDPhi[2]       ;          //!<! Eta vs Phi distance between track and cluster, E cluster > 0.5 GeV, after and before
  
  TH2F * fhTrackMatchedDEtaPos[2]        ;          //!<! Eta distance between track and cluster vs cluster E, after and before photon cuts
  TH2F * fhTrackMatchedDPhiPos[2]        ;          //!<! Phi distance between track and cluster vs cluster E, after and before photon cuts
  TH2F * fhTrackMatchedDEtaDPhiPos[2]    ;          //!<! Eta vs Phi distance between track and cluster, E cluster > 0.5 GeV, after and before
  
  TH2F * fhTrackMatchedDEtaNeg[2]        ;          //!<! Eta distance between track and cluster vs cluster E, after and before photon cuts
  TH2F * fhTrackMatchedDPhiNeg[2]        ;          //!<! Phi distance between track and cluster vs cluster E, after and before photon cuts
  TH2F * fhTrackMatchedDEtaDPhiNeg[2]    ;          //!<! Eta vs Phi distance between track and cluster, E cluster > 0.5 GeV, after and before photon cuts

//TH2F * fhTrackMatchedDEtaTrackPt[2]    ;          //!<! Eta distance between track and cluster vs track pT, after and before photon cuts
//TH2F * fhTrackMatchedDPhiTrackPt[2]    ;          //!<! Phi distance between track and cluster vs track pT, after and before photon cuts
//TH2F * fhTrackMatchedDEtaDPhiTrackPt[2];          //!<! Eta vs Phi distance between track and cluster, track pT > 0.5 GeV, after and before
  
  TH2F * fhTrackMatchedDEtaPosTrackPt[2] ;          //!<! Eta distance between track and cluster vs track pT, after and before photon cuts
  TH2F * fhTrackMatchedDPhiPosTrackPt[2] ;          //!<! Phi distance between track and cluster vs track pT, after and before photon cuts
  TH2F * fhTrackMatchedDEtaDPhiPosTrackPt[2];       //!<! Eta vs Phi distance between track and cluster, track pT > 0.5 GeV, after and before
  
  TH2F * fhTrackMatchedDEtaNegTrackPt[2] ;          //!<! Eta distance between track and cluster vs track pT, after and before photon cuts
  TH2F * fhTrackMatchedDPhiNegTrackPt[2] ;          //!<! Phi distance between track and cluster vs track pT, after and before photon cuts
  TH2F * fhTrackMatchedDEtaDPhiNegTrackPt[2];       //!<! Eta vs Phi distance between track and cluster, track pT > 0.5 GeV, after and before photon cuts

  TH2F * fhTrackMatchedDEtaTRD[2]        ;          //!<! Eta distance between track and cluster vs cluster E, after and before photon cuts, behind TRD
  TH2F * fhTrackMatchedDPhiTRD[2]        ;          //!<! Phi distance between track and cluster vs cluster E, after and before photon cuts, behind TRD
  
  TH2F * fhTrackMatchedDEtaMCOverlap[2]  ;          //!<! Eta distance between track and cluster vs cluster E, several particle overlap, after and before photon cuts 
  TH2F * fhTrackMatchedDPhiMCOverlap[2]  ;          //!<! Phi distance between track and cluster vs cluster E, several particle overlap, after and before photon cuts
  TH2F * fhTrackMatchedDEtaMCNoOverlap[2];          //!<! Eta distance between track and cluster vs cluster E, not other particle overlap, after and before photon cuts
  TH2F * fhTrackMatchedDPhiMCNoOverlap[2];          //!<! Phi distance between track and cluster vs cluster E, not other particle overlap, after and before photon cuts
  TH2F * fhTrackMatchedDEtaMCConversion[2];         //!<! Eta distance between track and cluster vs cluster E, originated in conversion, after and before photon cuts 
  TH2F * fhTrackMatchedDPhiMCConversion[2];         //!<! Phi distance between track and cluster vs cluster E, originated in conversion, after and before photon cuts 
  
  TH2F * fhTrackMatchedMCParticle[2];               //!<! Trace origin of matched particle
  TH2F * fhdEdx[2];                                 //!<! Matched track dEdx vs cluster E, after and before photon cuts
  TH2F * fhEOverP[2];                               //!<! Matched track E cluster over P track vs cluster E, after dEdx cut, after and before photon cuts
  TH2F * fhEOverPTRD[2];                            //!<! Matched track E cluster over P track vs cluster E, after dEdx cut, after and before photon cuts, behind TRD

  TH2F * fhdEdxTrackPt[2];                          //!<! Matched track dEdx vs track pT, after and before photon cuts
  TH2F * fhEOverPTrackPt[2];                        //!<! Matched track E cluster over P track vs track pT, after dEdx cut, after and before photon cuts

  
  // Pile-up
    
  TH1F * fhPtPhotonPileUp[7];                       //!<! pT distribution of selected photons
  TH2F * fhClusterTimeDiffPhotonPileUp[7];          //!<! E vs Time difference inside cluster for selected photons
  TH2F * fhTimePtPhotonNoCut;                       //!<! Time of photon cluster vs Pt, no cut
  TH2F * fhTimePtPhotonSPD;                         //!<! Time of photon cluster vs Pt, IsSPDPileUp
  TH2F * fhTimeNPileUpVertSPD;                      //!<! Time of cluster vs n pile-up vertices from SPD
  TH2F * fhTimeNPileUpVertTrack;                    //!<! Time of cluster vs n pile-up vertices from Tracks

  TH2F * fhPtPhotonNPileUpSPDVtx;	                  //!<! photon pt vs number of spd pile-up vertices
  TH2F * fhPtPhotonNPileUpTrkVtx;                   //!<! photon pt vs number of track pile-up vertices
  TH2F * fhPtPhotonNPileUpSPDVtxTimeCut;            //!<! photon pt vs number of spd pile-up vertices, time cut +-25 ns
  TH2F * fhPtPhotonNPileUpTrkVtxTimeCut;            //!<! photon pt vs number of track pile-up vertices, time cut +- 25 ns
  TH2F * fhPtPhotonNPileUpSPDVtxTimeCut2;           //!<! photon pt vs number of spd pile-up vertices, time cut +-75 ns
  TH2F * fhPtPhotonNPileUpTrkVtxTimeCut2;           //!<! photon pt vs number of track pile-up vertices, time cut +- 75 ns
	
  TH2F * fhEClusterSM ;                             //!<! Cluster E distribution per SM, before any selection, after reader
  TH2F * fhEPhotonSM  ;                             //!<! photon-like cluster E distribution per SM
  TH2F * fhPtClusterSM;                             //!<! Cluster E distribution per SM, before any selection, after reader
  TH2F * fhPtPhotonSM ;                             //!<! photon-like cluster E distribution per SM
  
  TH2F * fhMCConversionVertex;                      //!<! Conversion distance for photon clusters that have at least a contributor from the conversion. 
  TH2F * fhMCConversionVertexTRD;                   //!<! Conversion distance for photon clusters that have at least a contributor from the conversion, SM covered by TRD.
  TH2F * fhMCConversionLambda0Rcut[6];              //!<! Shower shape M02 of photon conversions, depending on conversion vertex.
  TH2F * fhMCConversionLambda0RcutTRD[6];           //!<! Shower shape M02 of photon conversions, depending on conversion vertex, SM covered by TRD
  TH2F * fhMCConversionLambda1Rcut[6];              //!<! Shower shape M20 of photon conversions, depending on conversion vertex.
  TH2F * fhMCConversionLambda1RcutTRD[6];           //!<! Shower shape M20 of photon conversions, depending on conversion vertex, SM covered by TRD

//TH2F * fhLam0EMCALRegion   [4][3];                //!<! Cluster lambda0 vs  E, in different EMCal regions
//TH2F * fhLam0EMCALRegionTRD[4][3];                //!<! Cluster lambda0 vs  E, in different EMCal regions, SM covered by TRD
//TH2F * fhLam0EMCALRegionMCConvRcut   [4][3][6];   //!<! Cluster lambda0 vs  E, in different EMCal regions, MC photon conversions, depending on conversion vertex
//TH2F * fhLam0EMCALRegionTRDMCConvRcut[4][3][6];   //!<! Cluster lambda0 vs  E, in different EMCal regions, SM covered by TRD,  MC photon conversions, depending on conversion vertex
  
  TH2F * fhLam0EMCALRegionPerSM[4][3][20];          //!<! Cluster lambda0 vs  Pt, in different EMCal regions
  TH2F * fhLam1EMCALRegionPerSM[4][3][20];          //!<! Cluster lambda1 vs  Pt, in different EMCal regions
  
  //
  // Clusters within a shower shape bin, 2 bins [0.23,0.26] (photon) and [0.3,0.4] (tail), if cell weight > 0
  //
  TH2F * fhLam1Lam0BinPerSM                [2][20]; //!<! Cluster lambda1, in a l0 bin per SM 
  TH2F * fhTimeLam0BinPerSM                [2][20]; //!<! Cell time, not maximum cluster cell, in a l0 bin per SM 
  TH2F * fhTimeLam0BinPerSMWeighted        [2][20]; //!<! Cell time, not maximum cluster cell, in a l0 bin per SM, log weight Cell E / Cluster E
  TH2F * fhDTimeLam0BinPerSM               [2][20]; //!<! t_max-t_cell, not maximum cluster cell, in a l0 bin per SM 
  TH2F * fhDTimeLam0BinPerSMWeighted       [2][20]; //!<! t_max-t_cell, not maximum cluster cell, in a l0 bin per SM, log weight Cell E / Cluster E
  TH2F * fhCellClusterEFracLam0BinPerSM    [2][20]; //!<! Cell E / Cluster E vs cluster pT, not maximum cluster cell, in a l0 bin per SM  
//TH2F * fhCellClusterEFracLam0BinPerSMWeighted[2][20]; //!<! Cell E / Cluster E vs cluster pT, not maximum cluster cell, in a l0 bin per SM, log weight  
  TH2F * fhCellClusterELam0BinPerSM        [2][20]; //!<! Cell E vs cluster pT, not maximum cluster cell, in a l0 bin per SM  
  TH2F * fhCellClusterELam0BinPerSMWeighted[2][20]; //!<! Cell E vs cluster pT, not maximum cluster cell, in a l0 bin per SM, log weight Cell E / Cluster E

  // plus different Pt bins 2-3,3-4,4-5,5-6,6-8,8-10,10-12
  TH2F * fhEtaPhiLam0BinPtBin              [2][7] ; //!<! Cluster eta/phi in a l0 bin, different Pt bins
  TH2F * fhColRowLam0BinPtBin              [2][7] ; //!<! Cell hits, not maximum cluster cell, in a l0 bin, different Pt bins 
  TH2F * fhColRowLam0BinPtBinWeighted      [2][7] ; //!<! Cell hits, not maximum cluster cell, in a l0 bin, different Pt bins and log weight Cell E / Cluster E 
  TH2F * fhCellClusterIndexEAndTime        [2][7] ; //!<! Cell in Cluster index (low index high energy, high index low energy) vs cell Time  
  TH2F * fhCellClusterEAndTime             [2][7] ; //!<! Cell in Cluster E cell vs cell Time,  in a l0 bin, different Pt bins
  TH2F * fhCellClusterEFracAndTime         [2][7] ; //!<! Cell in Cluster E cell/ E cluster vs cell Time,  in a l0 bin, different Pt bins

//  // Shared clusters
//  TH2F * fhLam0PerSMShared                   [12] ; //!<! Cluster lambda0 vs  Pt, for shared clusters, EMCal
//  TH2F * fhLam1PerSMShared                   [12] ; //!<! Cluster lambda1 vs  Pt, for shared clusters, EMCal
//  TH2F * fhTimeLam0BinPerSMShared         [2][12] ; //!<! Cell time, not maximum cluster cell, in a l0 bin per SM, shared SM
//  
//  TH2F * fhEtaPhiLam0BinPtBinSMShared      [2][7] ; //!<! Cluster eta/phi in a l0 bin, different Pt bins, SM shared clusters
//  TH2F * fhColRowLam0BinPtBinSMShared      [2][7] ; //!<! Cell hits, not maximum cluster cell, in a l0 bin, different Pt bins 
  
  // Cells with large time
  TH2F * fhEtaPhiLargeTimeInClusterCell    [2][7] ; //!<! Cluster eta/phi, with at least one significant cell with large time
  TH2F * fhColRowLam0BinPtBinLargeTime     [2][7] ; //!<! Cell hits, not maximum cluster cell, in a l0 bin, different Pt bins, cell with large time  
  TH2F * fhCellClusterEFracLam0BinPerSMLargeTime [2][20]; //!<! Cell in Cluster E cell / Ecluster vs cluster pT, with large time
  TH2F * fhCellClusterEFracLam0BinPerSMLargeTimeTotal[2][20]; //!<! Sum of all Cell in Cluster , with large time E cell / Ecluster vs cluster pT
  TH2F * fhCellClusterELam0BinPerSMLargeTime     [2][20]; //!<! Cell in Cluster E cell vs cluster pT, with large time
  TH2F * fhCellClusterIndexELam0BinPerSMLargeTime[2][20]; //!<! Cell in Cluster index (low index high energy, high index low energy) vs cluster pT, with large time  
  TH2F * fhNCellsWithLargeTimeInCluster    [2][20]; //!<! Number of cells in cluster with large time  
  
  TH2F * fhLam0PerSM                         [20] ; //!<! Cluster lambda0 vs  Pt, in different SM
  TH2F * fhLam1PerSM                         [20] ; //!<! Cluster lambda0 vs  Pt, in different SM
  TH2F * fhLam0PerSMLargeTimeInClusterCell   [20] ; //!<! Cluster lambda0 vs  Pt, when any secondary cell has t > 50 ns, in different SM
  TH2F * fhLam1PerSMLargeTimeInClusterCell   [20] ; //!<! Cluster lambda1 vs  Pt, when any secondary cell has t > 50 ns, in different SM  
  TH2F * fhLam0PerNLargeTimeInClusterCell     [5] ; //!<! Cluster lambda0 vs  Pt, when any secondary cell has t > 50 ns, per number of large time secondary cells
  TH2F * fhLam1PerNLargeTimeInClusterCell     [5] ; //!<! Cluster lambda1 vs  Pt, when any secondary cell has t > 50 ns, per number of large time secondary cells 
//TH2F * fhLam0PerSMSPDPileUp                [20] ; //!<! Cluster lambda0 vs  Pt, when event tagged as pile-up by SPD, in different SM
//TH2F * fhLam1PerSMSPDPileUp                [20] ; //!<! Cluster lambda0 vs  Pt, when event tagged as pile-up by SPD, in different SM  
  
  TH2F *  fhEBinClusterEtaPhi[14] ;                 //!<! Eta-Phi location of cluster in different energy bins.
  TH2F *  fhEBinClusterColRow[14] ;                 //!<! Column and row location of cluster max E cell in different energy bins.
  TH2F *  fhEBinClusterEtaPhiPID[14] ;              //!<! Eta-Phi location of cluster in different energy bins, after PID cut
  TH2F *  fhEBinClusterColRowPID[14] ;              //!<! Column and row location of cluster max E cell in different energy bins, after PID cut

  TH2F *  fhLocalRegionClusterEtaPhi[6]  ;                       //!<! Pseudorapidity vs Phi of clusters with cone R within the EMCal, for different cocktail merging cases 
  TH2F *  fhLocalRegionClusterEnergySum[6] ;                     //!<! Sum of energy near the cluster, R<0.2, vs cluster E, for different cocktail merging cases
  TH2F *  fhLocalRegionClusterMultiplicity[6];                   //!<! Cluster multiplicity near cluster, R<0.2, vs cluster E, for different cocktail merging cases
  TH2F *  fhLocalRegionClusterEnergySumPerCentrality[6] ;        //!<! Sum of energy near the cluster, R<0.2, vs centrality percentile, for different cocktail merging cases
  TH2F *  fhLocalRegionClusterMultiplicityPerCentrality[6];      //!<! Cluster multiplicity near cluster, R<0.2, vs centrality percentile, for different cocktail merging cases

  TH2F *  fhLocalRegionClusterEnergySumHijing[6] ;               //!<! Sum of energy near the cluster, R<0.2, vs cluster E, hijing tagged mc clusters, for different cocktail merging cases
  TH2F *  fhLocalRegionClusterMultiplicityHijing[6];             //!<! Cluster multiplicity near cluster, R<0.2, vs cluster E, hijing tagged mc clusters, for different cocktail merging cases
  TH2F *  fhLocalRegionClusterEnergySumPerCentralityHijing[6] ;  //!<! Sum of energy near the cluster, R<0.2, vs centrality percentile, hijing tagged mc clusters, for different cocktail merging cases
  TH2F *  fhLocalRegionClusterMultiplicityPerCentralityHijing[6];//!<! Cluster multiplicity near cluster, R<0.2, vs centrality percentile, hijing tagged mc clusters, for different cocktail merging cases

  TH2F *  fhLocalRegionClusterEnergySumHijing2 ;               //!<! Sum of energy near the cluster, R<0.2, vs cluster E, hijing tagged mc clusters, 
  TH2F *  fhLocalRegionClusterMultiplicityHijing2;             //!<! Cluster multiplicity near cluster, R<0.2, vs cluster E, hijing tagged mc clusters, 
  TH2F *  fhLocalRegionClusterEnergySumPerCentralityHijing2 ;  //!<! Sum of energy near the cluster, R<0.2, vs centrality percentile, hijing tagged mc clusters, 
  TH2F *  fhLocalRegionClusterMultiplicityPerCentralityHijing2;//!<! Cluster multiplicity near cluster, R<0.2, vs centrality percentile, hijing tagged mc clusters, for different cocktail merging cases
  
  TH2F *  fhLocalRegionClusterEnergySumAdded[6] ;                //!<! Sum of energy near the cluster, R<0.2, vs cluster E, not hijing (added signal) tagged mc clusters, for different cocktail merging cases
  TH2F *  fhLocalRegionClusterMultiplicityAdded[6];             //!<! Cluster multiplicity near cluster, R<0.2, vs cluster E, not hijing (added signal) tagged mc clusters, for different cocktail merging cases
  TH2F *  fhLocalRegionClusterEnergySumPerCentralityAdded[6] ;   //!<! Sum of energy near the cluster, R<0.2, vs centrality percentile, not hijing (added signal) tagged mc clusters, for different cocktail merging cases
  TH2F *  fhLocalRegionClusterMultiplicityPerCentralityAdded[6]; //!<! Cluster multiplicity near cluster, R<0.2, vs centrality percentile, not hijing (added signal) tagged mc clusters, for different cocktail merging cases
  
  TH2F *  fhLocalRegionClusterEnergySumMCPi0Decay[6] ;                     //!<! Sum of energy near the cluster, R<0.2, vs cluster E, for different cocktail merging cases
  TH2F *  fhLocalRegionClusterMultiplicityMCPi0Decay[6];                   //!<! Cluster multiplicity near cluster, R<0.2, vs cluster E, for different cocktail merging cases
  TH2F *  fhLocalRegionClusterEnergySumPerCentralityMCPi0Decay[6] ;        //!<! Sum of energy near the cluster, R<0.2, vs centrality percentile, for different cocktail merging cases
  TH2F *  fhLocalRegionClusterMultiplicityPerCentralityMCPi0Decay[6];      //!<! Cluster multiplicity near cluster, R<0.2, vs centrality percentile, for different cocktail merging cases
  
  TH2F *  fhLocalRegionClusterEnergySumHijingMCPi0Decay[6] ;               //!<! Sum of energy near the cluster, R<0.2, vs cluster E, hijing tagged mc clusters, for different cocktail merging cases
  TH2F *  fhLocalRegionClusterMultiplicityHijingMCPi0Decay[6];             //!<! Cluster multiplicity near cluster, R<0.2, vs cluster E, hijing tagged mc clusters, for different cocktail merging cases
  TH2F *  fhLocalRegionClusterEnergySumPerCentralityHijingMCPi0Decay[6] ;  //!<! Sum of energy near the cluster, R<0.2, vs centrality percentile, hijing tagged mc clusters, for different cocktail merging cases
  TH2F *  fhLocalRegionClusterMultiplicityPerCentralityHijingMCPi0Decay[6];//!<! Cluster multiplicity near cluster, R<0.2, vs centrality percentile, hijing tagged mc clusters, for different cocktail merging cases
  
  TH2F *  fhLocalRegionClusterEnergySumAddedMCPi0Decay[6] ;                //!<! Sum of energy near the cluster, R<0.2, vs cluster E, not hijing (added signal) tagged mc clusters, for different cocktail merging cases
  TH2F *  fhLocalRegionClusterMultiplicityAddedMCPi0Decay[6];             //!<! Cluster multiplicity near cluster, R<0.2, vs cluster E, not hijing (added signal) tagged mc clusters, for different cocktail merging cases
  TH2F *  fhLocalRegionClusterEnergySumPerCentralityAddedMCPi0Decay[6] ;   //!<! Sum of energy near the cluster, R<0.2, vs centrality percentile, not hijing (added signal) tagged mc clusters, for different cocktail merging cases
  TH2F *  fhLocalRegionClusterMultiplicityPerCentralityAddedMCPi0Decay[6]; //!<! Cluster multiplicity near cluster, R<0.2, vs centrality percentile, not hijing (added signal) tagged mc clusters, for different cocktail merging cases
  

  
  TH1F *  fhMergeGeneratorCluster                 [10][fgkNGenTypes]; //!<! Cluster energy, at least 2 generators contributions, for different generator origins and different particles.  
  TH1F *  fhMergeGeneratorClusterNotHijingBkg     [10][fgkNGenTypes]; //!<! Cluster energy, at least 2 generators contributions, none is HIJING, for different generator origins and different particles.
  TH1F *  fhMergeGeneratorClusterHijingAndOtherBkg[10][fgkNGenTypes]; //!<! Cluster energy, at least 3 generators contributions, one is HIJING, for different generator origins and different particles.
  TH1F *  fhMergeGeneratorClusterHijingBkg        [10][fgkNGenTypes]; //!<! Cluster energy, at least 2 generators contributions, one is HIJING, for different generator origins and different particles.
  TH1F *  fhCleanGeneratorCluster                 [10][fgkNGenTypes]; //!<! Cluster energy, only one generator is the contributor, for different generator origins and different particles.

  
  TH2F *  fhMergeGeneratorClusterEPrimRecoRatio                 [10][fgkNGenTypes]; //!<! Cluster energy, at least 2 generators contributions, for different generator origins and different particles.  
  TH2F *  fhMergeGeneratorClusterNotHijingBkgEPrimRecoRatio     [10][fgkNGenTypes]; //!<! Cluster energy, at least 2 generators contributions, none is HIJING, for different generator origins and different particles.
  TH2F *  fhMergeGeneratorClusterHijingAndOtherBkgEPrimRecoRatio[10][fgkNGenTypes]; //!<! Cluster energy, at least 3 generators contributions, one is HIJING, for different generator origins and different particles.
  TH2F *  fhMergeGeneratorClusterHijingBkgEPrimRecoRatio        [10][fgkNGenTypes]; //!<! Cluster energy, at least 2 generators contributions, one is HIJING, for different generator origins and different particles.
  TH2F *  fhCleanGeneratorClusterEPrimRecoRatio                 [10][fgkNGenTypes]; //!<! Cluster energy, only one generator is the contributor, for different generator origins and different particles.
  
  TH2F *  fhMergeGeneratorClusterEPrimRecoDiff                 [10][fgkNGenTypes]; //!<! Cluster energy, at least 2 generators contributions, for different generatororigins and different particles.  
  TH2F *  fhMergeGeneratorClusterNotHijingBkgEPrimRecoDiff     [10][fgkNGenTypes]; //!<! Cluster energy, at least 2 generators contributions, none is HIJING, for different generator origins and different particles.
  TH2F *  fhMergeGeneratorClusterHijingAndOtherBkgEPrimRecoDiff[10][fgkNGenTypes]; //!<! Cluster energy, at least 3 generators contributions, one is HIJING, for different generator origins and different particles.
  TH2F *  fhMergeGeneratorClusterHijingBkgEPrimRecoDiff        [10][fgkNGenTypes]; //!<! Cluster energy, at least 2 generators contributions, one is HIJING, for different generator origins and different particles.
  TH2F *  fhCleanGeneratorClusterEPrimRecoDiff                 [10][fgkNGenTypes]; //!<! Cluster energy, only one generator is the contributor, for different generator origins and different particles.
  
//  
//  TH2F *  fhDistanceAddedPhotonAddedPrimarySignal   ; //!<! Generated added cocktail photon vs distance to other generator primary particle.
//  TH2F *  fhDistanceHijingPhotonAddedPrimarySignal  ; //!<! Generated hijing photon vs distance to other generator primary particle.
//  TH2F *  fhDistanceAddedPhotonAddedSecondarySignal ; //!<! Generated added cocktail photon vs distance to other generator secondary particle.
//  TH2F *  fhDistanceHijingPhotonAddedSecondarySignal; //!<! Generated hijing photon vs distance to other generator secondary particle.
//  TH2F *  fhDistanceAddedPhotonHijingSecondary      ; //!<! Generated added cocktail photon vs distance to hijing secondary particle.

  TH2F *  fhDistance2AddedSignals     ;   //!<! Distance between 2 clusters generated by added cocktails 
  TH2F *  fhDistanceAddedSignalsHijing;   //!<! Distance between 2 clusters generated by added cocktails and hijing
  TH2F *  fhDistance2Hijing           ;   //!<! Distance between 2 clusters generated by hijing

  /// Copy constructor not implemented.
  AliAnaPhoton(              const AliAnaPhoton & g) ;
    
  /// Assignment operator not implemented.
  AliAnaPhoton & operator = (const AliAnaPhoton & g) ;
  
  /// \cond CLASSIMP
  ClassDef(AliAnaPhoton,46) ;
  /// \endcond

} ;
 
#endif//ALIANAPHOTON_H



