#ifndef ALIANAPHOTON_H
#define ALIANAPHOTON_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
//
// Class for the photon identification.
// Clusters from calorimeters are identified as photons
// and kept in the AOD. Few histograms produced.
// Produces input for other analysis classes like AliAnaPi0, 
// AliAnaParticleHadronCorrelation ... 
//

//-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
class TH2F ;
class TH1F;
class TString ;
class TObjString;
class TList ;

// --- ANALYSIS system ---
#include "AliAnaCaloTrackCorrBaseClass.h"

class AliAnaPhoton : public AliAnaCaloTrackCorrBaseClass {

 public: 
  AliAnaPhoton() ;              // default ctor
  virtual ~AliAnaPhoton() { ; } // virtual dtor
	
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
  
  Bool_t       ClusterSelected(AliVCluster* cl, TLorentzVector mom, Int_t nlm) ;
  
  void         FillAcceptanceHistograms();

  void         FillShowerShapeHistograms( AliVCluster* cluster, Int_t mcTag) ;
  
  void         SwitchOnFillShowerShapeHistograms()    { fFillSSHistograms = kTRUE  ; }
  void         SwitchOffFillShowerShapeHistograms()   { fFillSSHistograms = kFALSE ; }  
  
  void         SwitchOnOnlySimpleSSHistoFill()        { fFillOnlySimpleSSHisto = kTRUE  ; }
  void         SwitchOffOnlySimpleHistoFill()         { fFillOnlySimpleSSHisto = kFALSE ; }
  
  void         FillTrackMatchingResidualHistograms(AliVCluster* calo, Int_t cut);
  
  void         SwitchOnTMHistoFill()                  { fFillTMHisto      = kTRUE  ; }
  void         SwitchOffTMHistoFill()                 { fFillTMHisto      = kFALSE ; }

  void         FillPileUpHistograms(Float_t energy, Float_t pt, Float_t time) ;
  void         FillPileUpHistogramsPerEvent(TObjArray * clusters) ;

  void         SwitchOnFillPileUpHistograms()         { fFillPileUpHistograms = kTRUE  ; }
  void         SwitchOffFillPileUpHistograms()        { fFillPileUpHistograms = kFALSE ; }    
  
  // Analysis parameters setters getters
  
  TString      GetCalorimeter()                 const { return fCalorimeter        ; }
  void         SetCalorimeter(TString  & det)         { fCalorimeter = det         ; }
    
  // ** Cluster selection methods **
  
  void         SetMinDistanceToBadChannel(Float_t m1, Float_t m2, Float_t m3) {
                fMinDist = m1; fMinDist2 = m2; fMinDist3 = m3; }

  void         SetTimeCut(Double_t min, Double_t max) { fTimeCutMin = min; 
                                                        fTimeCutMax = max          ; }
  Double_t     GetTimeCutMin()                  const { return fTimeCutMin         ; }
  Double_t     GetTimeCutMax()                  const { return fTimeCutMax         ; }	
	
  void         SetNCellCut(Int_t n)                   { fNCellsCut = n             ; }
  Double_t     GetNCellCut()                    const { return fNCellsCut          ; }
  
  void           SetNLMCut(Int_t min, Int_t max)             { fNLMCutMin = min; 
    fNLMCutMax = max                ; }
  Int_t          GetNLMCutMin()                        const { return fNLMCutMin               ; }
  Int_t          GetNLMCutMax()                        const { return fNLMCutMax               ; }	
  
  
  Bool_t       IsTrackMatchRejectionOn()        const { return fRejectTrackMatch   ; }
  void         SwitchOnTrackMatchRejection()          { fRejectTrackMatch = kTRUE  ; }
  void         SwitchOffTrackMatchRejection()         { fRejectTrackMatch = kFALSE ; }  
	  
  void         FillNOriginHistograms(Int_t n)         { fNOriginHistograms = n ; 
    if(n > 14) fNOriginHistograms = 14; }
  void         FillNPrimaryHistograms(Int_t n)        { fNPrimaryHistograms= n ;
    if(n > 7)  fNPrimaryHistograms = 7; }

  // For histograms in arrays, index in the array, corresponding to a particle
  enum mcTypes    { kmcPhoton = 0,        kmcPi0Decay = 1,       kmcOtherDecay = 2,  
                    kmcPi0 = 3,           kmcEta = 4,            kmcElectron = 5,       
                    kmcConversion = 6,    kmcOther = 7,          kmcAntiNeutron = 8,    
                    kmcAntiProton = 9,    kmcPrompt = 10,        kmcFragmentation = 11, 
                    kmcISR = 12,          kmcString = 13                               };  

  enum mcPTypes   { kmcPPhoton = 0,       kmcPPi0Decay = 1,       kmcPOtherDecay = 2,  kmcPOther = 3,
                    kmcPPrompt = 4,       kmcPFragmentation = 5,  kmcPISR = 6           };  
  
  enum mcssTypes  { kmcssPhoton = 0,      kmcssOther = 1,       kmcssPi0 = 2,         
                    kmcssEta = 3,         kmcssConversion = 4,  kmcssElectron = 5       };  
  
  private:
 
  TString  fCalorimeter ;                // Calorimeter where the gamma is searched;
  Float_t  fMinDist ;                    // Minimal distance to bad channel to accept cluster
  Float_t  fMinDist2;                    // Cuts on Minimal distance to study acceptance evaluation
  Float_t  fMinDist3;                    // One more cut on distance used for acceptance-efficiency study
  Bool_t   fRejectTrackMatch ;           // If PID on, reject clusters which have an associated TPC track
  Bool_t   fFillTMHisto;                 // Fill track matching plots
  Double_t fTimeCutMin  ;                // Remove clusters/cells with time smaller than this value, in ns
  Double_t fTimeCutMax  ;                // Remove clusters/cells with time larger than this value, in ns
  Int_t    fNCellsCut ;                  // Accept for the analysis clusters with more than fNCellsCut cells
  Int_t    fNLMCutMin  ;                 // Remove clusters/cells with number of local maxima smaller than this value
  Int_t    fNLMCutMax  ;                 // Remove clusters/cells with number of local maxima larger than this value
  Bool_t   fFillSSHistograms ;           // Fill shower shape histograms
  Bool_t   fFillOnlySimpleSSHisto;       // Fill selected cluster histograms, selected SS histograms
  Int_t    fNOriginHistograms;           // Fill only NOriginHistograms of the 14 defined types
  Int_t    fNPrimaryHistograms;          // Fill only NPrimaryHistograms of the 7 defined types
  Bool_t   fFillPileUpHistograms;        // Fill pile-up related histograms
  
  //Histograms 
  TH1F * fhClusterCuts[10];              //! control histogram on the different photon selection cuts
  TH2F * fhNCellsE;                      //! number of cells in cluster vs E 
  TH2F * fhCellsE;                       //! energy of cells in cluster vs E of cluster
  TH2F * fhMaxCellDiffClusterE;          //! Fraction of energy carried by cell with maximum energy
  TH2F * fhTimeE;                        //! time of cluster vs E 

  TH1F * fhEPhoton    ;                  //! Number of identified photon vs energy
  TH1F * fhPtPhoton   ;                  //! Number of identified photon vs transerse momentum 
  TH2F * fhPhiPhoton  ;                  //! Azimuthal angle of identified  photon vs transerse momentum 
  TH2F * fhEtaPhoton  ;                  //! Pseudorapidity of identified  photon vs transerse momentum 
  TH2F * fhEtaPhiPhoton  ;               //! Pseudorapidity vs Phi of identified  photon for transerse momentum > 0.5
  TH2F * fhEtaPhi05Photon  ;             //! Pseudorapidity vs Phi of identified  photon for transerse momentum < 0.5
  TH2F * fhPtCentralityPhoton    ;       //! centrality  vs photon pT
  TH2F * fhPtEventPlanePhoton    ;       //! event plane vs photon pT
  
  //Shower shape
  TH2F * fhNLocMax;                       //! number of maxima in selected clusters

  TH2F * fhDispE;                         //! cluster dispersion vs E
  TH2F * fhLam0E;                         //! cluster lambda0 vs  E
  TH2F * fhLam1E;                         //! cluster lambda1 vs  E  

  TH2F * fhDispETRD;                      //! cluster dispersion vs E, SM covered by TRD
  TH2F * fhLam0ETRD;                      //! cluster lambda0 vs  E, SM covered by TRD
  TH2F * fhLam1ETRD;                      //! cluster lambda1 vs  E, SM covered by TRD 

  TH2F * fhDispETM;                       //! cluster dispersion vs E, cut on Track Matching residual
  TH2F * fhLam0ETM;                       //! cluster lambda0 vs  E, cut on Track Matching residual
  TH2F * fhLam1ETM;                       //! cluster lambda1 vs  E, cut on Track Matching residual  
  
  TH2F * fhDispETMTRD;                    //! cluster dispersion vs E, SM covered by TRD, cut on Track Matching residual
  TH2F * fhLam0ETMTRD;                    //! cluster lambda0 vs  E, SM covered by TRD, cut on Track Matching residual
  TH2F * fhLam1ETMTRD;                    //! cluster lambda1 vs  E, SM covered by TRD, cut on Track Matching residual 
  
  TH2F * fhNCellsLam0LowE;                //! number of cells in cluster vs lambda0
  TH2F * fhNCellsLam1LowE;                //! number of cells in cluster vs lambda1
  TH2F * fhNCellsDispLowE;                //! number of cells in cluster vs dispersion
  TH2F * fhNCellsLam0HighE;               //! number of cells in cluster vs lambda0, E>2
  TH2F * fhNCellsLam1HighE;               //! number of cells in cluster vs lambda1, E>2
  TH2F * fhNCellsDispHighE;               //! number of cells in cluster vs dispersion, E>2
  
  TH2F * fhEtaLam0LowE;                   //! cluster eta vs lambda0, E<2
  TH2F * fhPhiLam0LowE;                   //! cluster phi vs lambda0, E<2
  TH2F * fhEtaLam0HighE;                  //! cluster eta vs lambda0, E>2
  TH2F * fhPhiLam0HighE;                  //! cluster phi vs lambda0, E>2
  TH2F * fhLam0DispLowE;                  //! cluster lambda0 vs dispersion, E<2
  TH2F * fhLam0DispHighE;                 //! cluster lambda0 vs dispersion, E>2
  TH2F * fhLam1Lam0LowE;                  //! cluster lambda1 vs lambda0, E<2
  TH2F * fhLam1Lam0HighE;                 //! cluster lambda1 vs lambda0, E>2
  TH2F * fhDispLam1LowE;                  //! cluster disp vs lambda1, E<2
  TH2F * fhDispLam1HighE;                 //! cluster disp vs lambda1, E>2
    
  TH2F * fhDispEtaE ;                     //! shower dispersion in eta direction
  TH2F * fhDispPhiE ;                     //! shower dispersion in phi direction
  TH2F * fhSumEtaE ;                      //! shower dispersion in eta direction
  TH2F * fhSumPhiE ;                      //! shower dispersion in phi direction
  TH2F * fhSumEtaPhiE ;                   //! shower dispersion in eta and phi direction
  TH2F * fhDispEtaPhiDiffE ;              //! shower dispersion eta - phi
  TH2F * fhSphericityE ;                  //! shower sphericity in eta vs phi
  TH2F * fhDispSumEtaDiffE ;              //! difference of 2 eta dispersions
  TH2F * fhDispSumPhiDiffE ;              //! difference of 2 phi dispersions
  TH2F * fhDispEtaDispPhi[7] ;            //! shower dispersion in eta direction vs phi direction for 5 E bins [0-2],[2-4],[4-6],[6-10],[> 10]
  TH2F * fhLambda0DispEta[7] ;            //! shower shape correlation l0 vs disp eta
  TH2F * fhLambda0DispPhi[7] ;            //! shower shape correlation l0 vs disp phi
  
  //Fill MC dependent histograms, Origin of this cluster is ...

  TH2F * fhMCDeltaE[14]  ;                      //! MC-Reco E distribution coming from MC particle     
  TH2F * fhMCDeltaPt[14] ;                      //! MC-Reco pT distribution coming from MC particle
  TH2F * fhMC2E[14]  ;                          //! E distribution, Reco vs MC coming from MC particle
  TH2F * fhMC2Pt[14] ;                          //! pT distribution, Reco vs MC coming from MC particle
  
  TH1F * fhMCE[14];                             //! Number of identified photon vs cluster energy coming from MC particle
  TH1F * fhMCPt[14];                            //! Number of identified photon vs cluster pT     coming from MC particle
  TH2F * fhMCPhi[14];                           //! Phi of identified photon coming from MC particle
  TH2F * fhMCEta[14];                           //! eta of identified photon coming from MC particle

  TH1F * fhEPrimMC[7];                          //! Number of generated photon vs energy
  TH1F * fhPtPrimMC[7];                         //! Number of generated photon vs pT   
  TH2F * fhPhiPrimMC[7];                        //! Phi of generted photon
  TH2F * fhYPrimMC[7];                          //! Rapidity of generated photon 
  
  TH1F * fhEPrimMCAcc[7];                       //! Number of generated photon vs energy, in calorimeter acceptance
  TH1F * fhPtPrimMCAcc[7];                      //! Number of generated photon vs pT, in calorimeter acceptance   
  TH2F * fhPhiPrimMCAcc[7];                     //! Phi of generted photon, in calorimeter acceptance
  TH2F * fhYPrimMCAcc[7];                       //! Rapidity of generated photon, in calorimeter acceptance   
  
  // Shower Shape MC

  TH2F * fhMCELambda0[6] ;                      //! E vs Lambda0     from MC particle
  TH2F * fhMCELambda1[6] ;                      //! E vs Lambda1     from MC particle
  TH2F * fhMCEDispersion[6] ;                   //! E vs Dispersion  from MC particle
  
  TH2F * fhMCPhotonELambda0NoOverlap ;          //! E vs Lambda0     from MC photons, no overlap
  TH2F * fhMCPhotonELambda0TwoOverlap ;         //! E vs Lambda0     from MC photons, 2 particles overlap
  TH2F * fhMCPhotonELambda0NOverlap ;           //! E vs Lambda0     from MC photons, N particles overlap
  
  TH2F * fhMCLambda0vsClusterMaxCellDiffE0[6];  //! Lambda0 vs fraction of energy of max cell for E < 2 GeV
  TH2F * fhMCLambda0vsClusterMaxCellDiffE2[6];  //! Lambda0 vs fraction of energy of max cell for 2< E < 6 GeV
  TH2F * fhMCLambda0vsClusterMaxCellDiffE6[6];  //! Lambda0 vs fraction of energy of max cell for E > 6 GeV
  TH2F * fhMCNCellsvsClusterMaxCellDiffE0[6];   //! NCells  vs fraction of energy of max cell for E < 2
  TH2F * fhMCNCellsvsClusterMaxCellDiffE2[6];   //! NCells  vs fraction of energy of max cell for 2 < E < 6 GeV
  TH2F * fhMCNCellsvsClusterMaxCellDiffE6[6];   //! NCells  vs fraction of energy of max cell for E > 6
  TH2F * fhMCNCellsE[6];                        //! NCells per cluster vs energy
  TH2F * fhMCMaxCellDiffClusterE[6];            //! Fraction of energy carried by cell with maximum energy

  TH2F * fhMCEDispEta[6] ;                      //! shower dispersion in eta direction
  TH2F * fhMCEDispPhi[6] ;                      //! shower dispersion in phi direction
  TH2F * fhMCESumEtaPhi[6] ;                    //! shower dispersion in eta vs phi direction
  TH2F * fhMCEDispEtaPhiDiff[6] ;               //! shower dispersion in eta -phi direction
  TH2F * fhMCESphericity[6] ;                   //! shower sphericity, eta vs phi
  TH2F * fhMCDispEtaDispPhi[7][6] ;             //! shower dispersion in eta direction vs phi direction for 5 E bins [0-2],[2-4],[4-6],[6-10],[> 10]
  TH2F * fhMCLambda0DispEta[7][6] ;             //! shower shape correlation l0 vs disp eta
  TH2F * fhMCLambda0DispPhi[7][6] ;             //! shower shape correlation l0 vs disp phi

  //Embedding
  TH2F * fhEmbeddedSignalFractionEnergy ;       //! Fraction of photon energy of embedded signal vs cluster energy
  
  TH2F * fhEmbedPhotonELambda0FullSignal ;      //!  Lambda0 vs E for embedded photons with more than 90% of the cluster energy
  TH2F * fhEmbedPhotonELambda0MostlySignal ;    //!  Lambda0 vs E for embedded photons with 90%<fraction<50% 
  TH2F * fhEmbedPhotonELambda0MostlyBkg ;       //!  Lambda0 vs E for embedded photons with 50%<fraction<10% 
  TH2F * fhEmbedPhotonELambda0FullBkg ;         //!  Lambda0 vs E for embedded photons with less than 10% of the cluster energy
  
  TH2F * fhEmbedPi0ELambda0FullSignal ;         //!  Lambda0 vs E for embedded photons with more than 90% of the cluster energy
  TH2F * fhEmbedPi0ELambda0MostlySignal ;       //!  Lambda0 vs E for embedded photons with 90%<fraction<50% 
  TH2F * fhEmbedPi0ELambda0MostlyBkg ;          //!  Lambda0 vs E for embedded photons with 50%<fraction<10% 
  TH2F * fhEmbedPi0ELambda0FullBkg ;            //!  Lambda0 vs E for embedded photons with less than 10% of the cluster energy
  
  // Track Matching
  TH2F * fhTrackMatchedDEta[2]           ;      //! Eta distance between track and cluster vs cluster E, after and before photon cuts
  TH2F * fhTrackMatchedDPhi[2]           ;      //! Phi distance between track and cluster vs cluster E, after and before photon cuts
  TH2F * fhTrackMatchedDEtaDPhi[2]       ;      //! Eta vs Phi distance between track and cluster, E cluster > 0.5 GeV, after and before photon cuts
  
  TH2F * fhTrackMatchedDEtaTRD[2]        ;      //! Eta distance between track and cluster vs cluster E, after and before photon cuts, behind TRD
  TH2F * fhTrackMatchedDPhiTRD[2]        ;      //! Phi distance between track and cluster vs cluster E, after and before photon cuts, behind TRD
  
  TH2F * fhTrackMatchedDEtaMCOverlap[2]  ;      //! Eta distance between track and cluster vs cluster E, several particle overlap, after and before photon cuts 
  TH2F * fhTrackMatchedDPhiMCOverlap[2]  ;      //! Phi distance between track and cluster vs cluster E, several particle overlap, after and before photon cuts 
  TH2F * fhTrackMatchedDEtaMCNoOverlap[2];      //! Eta distance between track and cluster vs cluster E, not other particle overlap, after and before photon cuts 
  TH2F * fhTrackMatchedDPhiMCNoOverlap[2];      //! Phi distance between track and cluster vs cluster E, not other particle overlap, after and before photon cuts 
  TH2F * fhTrackMatchedDEtaMCConversion[2];     //! Eta distance between track and cluster vs cluster E, originated in conversion, after and before photon cuts 
  TH2F * fhTrackMatchedDPhiMCConversion[2];     //! Phi distance between track and cluster vs cluster E, originated in conversion, after and before photon cuts 
  
  TH2F * fhTrackMatchedMCParticle[2];           //! Trace origin of matched particle
  TH2F * fhdEdx[2];                             //! matched track dEdx vs cluster E, after and before photon cuts 
  TH2F * fhEOverP[2];                           //! matched track E cluster over P track vs cluster E, after dEdx cut, after and before photon cuts 
  TH2F * fhEOverPTRD[2];                        //! matched track E cluster over P track vs cluster E, after dEdx cut, after and before photon cuts, behind TRD 

  // Pile-up
  TH1F * fhPtPileUp[7];                         //! pT distribution of clusters before any selection
  TH1F * fhPtChargedPileUp[7];                  //! pT distribution of track matched clusters
  TH1F * fhPtPhotonPileUp[7];                   //! pT distribution of selected photons
  TH2F * fhLambda0PileUp[7];                    //! E vs M02 distribution of clusters, before any selection
  TH2F * fhLambda0ChargedPileUp[7];             //! E vs M02 distribution of clusters, track matched clusters
  TH2F * fhClusterTimeDiffPileUp[7];            //! E vs Time difference inside cluster, before any selection
  TH2F * fhClusterTimeDiffChargedPileUp[7];     //! E vs Time difference inside cluster for track matched clusters
  TH2F * fhClusterTimeDiffPhotonPileUp[7];      //! E vs Time difference inside cluster for selected photons
  TH2F * fhClusterEFracLongTimePileUp[7];       //! E vs fraction of cluster energy from cells with large time
  TH2F * fhTimeENoCut;                          //! time of cluster vs E, no cut
  TH2F * fhTimeESPD;                            //! time of cluster vs E, IsSPDPileUp
  TH2F * fhTimeESPDMulti;                       //! time of cluster vs E, IsSPDPileUpMulti
  TH2F * fhTimeNPileUpVertSPD;                  //! time of cluster vs n pile-up vertices from SPD
  TH2F * fhTimeNPileUpVertTrack;                //! time of cluster vs n pile-up vertices from Tracks
  TH2F * fhTimeNPileUpVertContributors;         //! time of cluster vs n pile-up vertex from SPD contributors
  TH2F * fhTimePileUpMainVertexZDistance;       //! time of cluster vs difference of z main vertex and pile-up vertex 
  TH2F * fhTimePileUpMainVertexZDiamond;        //! time of cluster vs difference of z diamond and pile-up vertex 
  TH2F * fhClusterMultSPDPileUp[4];             //! E max cluster vs event cluster multiplicity, for tmax-tdiff cuts, pile up event
  TH2F * fhClusterMultNoPileUp[4];              //! E max cluster vs event cluster multiplicity, for tmax-tdiff cuts, not pile up event
  TH2F * fhEtaPhiBC0;                           //! eta/phi of clusters in BC=0
  TH2F * fhEtaPhiBCPlus;                        //! eta/phi of clusters in BC>0
  TH2F * fhEtaPhiBCMinus;                       //! eta/phi of clusters in BC<0
  TH2F * fhEtaPhiBC0PileUpSPD;                  //! eta/phi of clusters in BC=0, SPD pile-up
  TH2F * fhEtaPhiBCPlusPileUpSPD;               //! eta/phi of clusters in BC>0, SPD pile-up
  TH2F * fhEtaPhiBCMinusPileUpSPD;              //! eta/phi of clusters in BC<0, SPD pile-up

  AliAnaPhoton(              const AliAnaPhoton & g) ; // cpy ctor
  AliAnaPhoton & operator = (const AliAnaPhoton & g) ; // cpy assignment
  
  ClassDef(AliAnaPhoton,29)

} ;
 
#endif//ALIANAPHOTON_H



