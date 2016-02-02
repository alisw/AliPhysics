#ifndef ALIANAPHOTON_H
#define ALIANAPHOTON_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliAnaPhoton
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
  
  void         FillShowerShapeHistograms( AliVCluster* cluster, Int_t mcTag, Float_t maxCellEFraction) ;
  
  void         SwitchOnFillShowerShapeHistograms()    { fFillSSHistograms      = kTRUE  ; }
  void         SwitchOffFillShowerShapeHistograms()   { fFillSSHistograms      = kFALSE ; }  

  void         SwitchOnFillEMCALRegionSSHistograms()  { fFillEMCALRegionSSHistograms = kTRUE  ; }
  void         SwitchOffFillEMCALRegionSSHistograms() { fFillEMCALRegionSSHistograms = kFALSE ; }  
  
  void         SwitchOnOnlySimpleSSHistoFill()        { fFillOnlySimpleSSHisto = kTRUE  ; }
  void         SwitchOffOnlySimpleHistoFill()         { fFillOnlySimpleSSHisto = kFALSE ; }
  
  void         FillTrackMatchingResidualHistograms(AliVCluster* calo, Int_t cut);
  
  void         SwitchOnTMHistoFill()                  { fFillTMHisto           = kTRUE  ; }
  void         SwitchOffTMHistoFill()                 { fFillTMHisto           = kFALSE ; }

  void         FillPileUpHistograms(AliVCluster* cluster, AliVCaloCells *cells, Int_t absIdMax) ;
 
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
  
  void         FillNOriginHistograms(Int_t n)         { fNOriginHistograms = n ; 
    if(n > fgkNmcTypes    ) fNOriginHistograms  = fgkNmcTypes     ; }
  void         FillNPrimaryHistograms(Int_t n)        { fNPrimaryHistograms= n ;
    if(n > fgkNmcPrimTypes) fNPrimaryHistograms = fgkNmcPrimTypes ; }

  /// For MC histograms in arrays, index in the array corresponds to a MC originating particle type
  enum mcTypes    { kmcPhoton     =  0,    kmcPi0Decay = 1,       kmcEtaDecay      = 2,  kmcOtherDecay = 3,
                    kmcPi0        =  4,    kmcEta      = 5,       kmcElectron      = 6,
                    kmcConversion =  7,    kmcOther    = 8,       kmcAntiNeutron   = 9,
                    kmcAntiProton = 10,    kmcPrompt   = 11,      kmcFragmentation = 12,
                    kmcISR        = 13,    kmcString   = 14  };

  /// Total number of cluster MC origin histograms
  static const Int_t fgkNmcTypes = 15;

  /// For MC histograms in arrays, index in the array corresponds to a MC generated primary particle type
  enum mcPTypes   { kmcPPhoton = 0,       kmcPPi0Decay      = 1,  kmcPEtaDecay = 2,     kmcPOtherDecay = 3,
                    kmcPPrompt = 4,       kmcPFragmentation = 5,  kmcPISR      = 6 };
  
  /// Total number of MC primary histograms
  static const Int_t fgkNmcPrimTypes = 7;
  
  /// For MC histograms with shower shape in arrays, index in the array corresponds to a MC originating particle type
  enum mcssTypes  { kmcssPhoton = 0,      kmcssOther = 1,       kmcssPi0 = 2,
                    kmcssEta = 3,         kmcssConversion = 4,  kmcssElectron = 5       };  
  
  /// Total number of MC histograms for shower shape studies.
  static const Int_t fgkNssTypes = 6 ;
  
  private:
 
  Float_t  fMinDist ;                               ///<  Minimal distance to bad channel to accept cluster
  Float_t  fMinDist2;                               ///<  Cuts on Minimal distance to study acceptance evaluation
  Float_t  fMinDist3;                               ///<  One more cut on distance used for acceptance-efficiency study
    
  Bool_t   fRejectTrackMatch ;                      ///<  If PID on, reject clusters which have an associated TPC track
    
  Bool_t   fFillTMHisto;                            ///<  Fill track matching plots
    
  Double_t fTimeCutMin  ;                           ///<  Remove clusters/cells with time smaller than this value, in ns
  Double_t fTimeCutMax  ;                           ///<  Remove clusters/cells with time larger than this value, in ns
    
  Int_t    fNCellsCut ;                             ///<  Accept for the analysis clusters with more than fNCellsCut cells
    
  Int_t    fNLMCutMin  ;                            ///<  Remove clusters/cells with number of local maxima smaller than this value
  Int_t    fNLMCutMax  ;                            ///<  Remove clusters/cells with number of local maxima larger than this value
    
  Bool_t   fFillSSHistograms ;                      ///<  Fill shower shape histograms

  Bool_t   fFillEMCALRegionSSHistograms ;           ///<  Fill shower shape histograms in EMCal slices
    
  Bool_t   fFillOnlySimpleSSHisto;                  ///<  Fill selected cluster histograms, selected SS histograms
    
  Int_t    fNOriginHistograms;                      ///<  Fill only NOriginHistograms of the 14 defined types
  Int_t    fNPrimaryHistograms;                     ///<  Fill only NPrimaryHistograms of the 7 defined types
  
  TLorentzVector fMomentum;                         //!<! Cluster momentum, temporary container
  TLorentzVector fPrimaryMom;                       //!<! Primary MC momentum, temporary container
  TVector3       fProdVertex;                       //!<! Primary MC production vertex, temporary container
  
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
  
  // Shower shape
  TH2F * fhNLocMax;                                 //!<! number of maxima in selected clusters

  TH2F * fhDispE;                                   //!<! Cluster dispersion vs E
  TH2F * fhLam0E;                                   //!<! Cluster lambda0 vs  E
  TH2F * fhLam0Pt;                                  //!<! Cluster lambda0 vs  pT
  TH2F * fhLam1E;                                   //!<! Cluster lambda1 vs  E

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
  
  TH2F * fhMCPhotonELambda0NoOverlap ;              //!<! E vs Lambda0     from MC photons, no overlap
  TH2F * fhMCPhotonELambda0TwoOverlap ;             //!<! E vs Lambda0     from MC photons, 2 particles overlap
  TH2F * fhMCPhotonELambda0NOverlap ;               //!<! E vs Lambda0     from MC photons, N particles overlap
  
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
    
  TH2F * fhTrackMatchedDEta[2]           ;          //!<! Eta distance between track and cluster vs cluster E, after and before photon cuts
  TH2F * fhTrackMatchedDPhi[2]           ;          //!<! Phi distance between track and cluster vs cluster E, after and before photon cuts
  TH2F * fhTrackMatchedDEtaDPhi[2]       ;          //!<! Eta vs Phi distance between track and cluster, E cluster > 0.5 GeV, after and before
  
  TH2F * fhTrackMatchedDEtaPos[2]        ;          //!<! Eta distance between track and cluster vs cluster E, after and before photon cuts
  TH2F * fhTrackMatchedDPhiPos[2]        ;          //!<! Phi distance between track and cluster vs cluster E, after and before photon cuts
  TH2F * fhTrackMatchedDEtaDPhiPos[2]    ;          //!<! Eta vs Phi distance between track and cluster, E cluster > 0.5 GeV, after and before
  
  TH2F * fhTrackMatchedDEtaNeg[2]        ;          //!<! Eta distance between track and cluster vs cluster E, after and before photon cuts
  TH2F * fhTrackMatchedDPhiNeg[2]        ;          //!<! Phi distance between track and cluster vs cluster E, after and before photon cuts
  TH2F * fhTrackMatchedDEtaDPhiNeg[2]    ;          //!<! Eta vs Phi distance between track and cluster, E cluster > 0.5 GeV, after and before photon cuts
  
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
  TH2F * fhMCConversionLambda0Rcut[6];              //!<! Shower shape of photon conversions, depending on conversion vertex.
  TH2F * fhMCConversionLambda0RcutTRD[6];           //!<! Shower shape of photon conversions, depending on conversion vertex, SM covered by TRD

//TH2F * fhLam0EMCALRegion   [4][3];                //!<! Cluster lambda0 vs  E, in different EMCal regions
//TH2F * fhLam0EMCALRegionTRD[4][3];                //!<! Cluster lambda0 vs  E, in different EMCal regions, SM covered by TRD
//TH2F * fhLam0EMCALRegionMCConvRcut   [4][3][6];   //!<! Cluster lambda0 vs  E, in different EMCal regions, MC photon conversions, depending on conversion vertex
//TH2F * fhLam0EMCALRegionTRDMCConvRcut[4][3][6];   //!<! Cluster lambda0 vs  E, in different EMCal regions, SM covered by TRD,  MC photon conversions, depending on conversion vertex
  
  TH2F * fhLam0EMCALRegionPerSM[4][3][20];          //!<! Cluster lambda0 vs  E, in different EMCal regions
  TH2F * fhEtaPhiLam0BinPtBin[7];                   //!<! Cluster eta/phi for a given l0 bin (0.3-0.4) and different E bins 2-3,3-4,4-5,5-6,6-8,8-10,10-12
  
  /// Copy constructor not implemented.
  AliAnaPhoton(              const AliAnaPhoton & g) ;
    
  /// Assignment operator not implemented.
  AliAnaPhoton & operator = (const AliAnaPhoton & g) ;
  
  /// \cond CLASSIMP
  ClassDef(AliAnaPhoton,42) ;
  /// \endcond

} ;
 
#endif//ALIANAPHOTON_H



