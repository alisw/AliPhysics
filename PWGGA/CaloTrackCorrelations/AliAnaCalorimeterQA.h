#ifndef ALIANACALORIMETERQA_H
#define ALIANACALORIMETERQA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
// Class to check results from simulations or reconstructed real data. 
// Fill few histograms and do some checking plots
//
//-- Author: Gustavo Conesa (INFN-LNF)

// --- Root system ---
class TH3F;
class TH2F;
class TH1F;
class TObjString;
class TObjArray;

// --- Analysis system --- 
class AliVCaloCells;
class AliVCaloCluster;
class AliVTrack;

#include "AliAnaCaloTrackCorrBaseClass.h"
 
class AliAnaCalorimeterQA : public AliAnaCaloTrackCorrBaseClass {
  
public: 
  AliAnaCalorimeterQA() ; // default ctor	
  virtual ~AliAnaCalorimeterQA() {;} //virtual dtor
    
  // General methods
  
  TObjString * GetAnalysisCuts();

  TList      * GetCreateOutputObjects();
  
  void         Init();
  
  void         InitParameters();
    
  void         MakeAnalysisFillHistograms() ;
  
  void         Print(const Option_t * opt) const;
    
  // Main methods
  
  void         BadClusterHistograms(AliVCluster* clus, const TObjArray *caloClusters,  AliVCaloCells * cells,
                                    Int_t absIdMax,    Double_t maxCellFraction, Float_t eCrossFrac, Double_t tmax);
    
  void         CalculateAverageTime(AliVCluster *clus, AliVCaloCells *cells, Double_t timeAverages[2]);
  
  void         CellHistograms(AliVCaloCells * cells);

  void         CellInClusterPositionHistograms(AliVCluster* cluster);
    
  void         ClusterAsymmetryHistograms(AliVCluster* clus, Int_t absIdMax, const Bool_t goodCluster );
  
  void         ClusterHistograms(AliVCluster* cluster, const TObjArray *caloClusters,  AliVCaloCells * cells, 
                                 Int_t absIdMax, Double_t maxCellFraction, Float_t eCrossFrac, Double_t tmax);
  
  void         ClusterLoopHistograms(const TObjArray * clusters, AliVCaloCells * cells);
  
  Bool_t       ClusterMCHistograms(TLorentzVector mom, Bool_t matched,
                                   const Int_t * labels, Int_t nLabels, Int_t & pdg );

  void         ClusterMatchedWithTrackHistograms(AliVCluster* clus, TLorentzVector mom, 
                                                 Bool_t mcOK, Int_t pdg);

  void         Correlate();
  
  void         ExoticHistograms(Int_t absIdMax, Float_t ampMax,
                                AliVCluster *clus, AliVCaloCells* cells);
  
  Float_t      GetECross(Int_t absId, AliVCaloCells* cells,Float_t dtcut = 10000);
  
  void         InvariantMassHistograms(Int_t iclus, TLorentzVector mom, Int_t nModule,
                                       const TObjArray* caloClusters, AliVCaloCells * cells);

  Bool_t       IsGoodCluster(Int_t absIdMax, AliVCaloCells *cells);
  
  void         MCHistograms();  
  
  void         MCHistograms(const TLorentzVector mom, Int_t pdg);
    
  void         WeightHistograms(AliVCluster *clus, AliVCaloCells* cells);

  // Setters and Getters

  
  Float_t      GetEMCALCellAmpMin()      const  { return fEMCALCellAmpMin    ; }
  void         SetEMCALCellAmpMin(Float_t amp)  { fEMCALCellAmpMin = amp     ; }
  
  Float_t      GetPHOSCellAmpMin()       const  { return fPHOSCellAmpMin     ; }
  void         SetPHOSCellAmpMin (Float_t amp)  { fPHOSCellAmpMin  = amp     ; }
    
  TString      GetCalorimeter()          const  { return fCalorimeter        ; }
  void         SetCalorimeter(TString calo)     { fCalorimeter = calo        ; }
  
  Double_t     GetTimeCutMin()           const  { return fTimeCutMin         ; }
  Double_t     GetTimeCutMax()           const  { return fTimeCutMax         ; }
  void         SetTimeCut(Double_t min, Double_t max) {
                          fTimeCutMin = min ; fTimeCutMax = max              ; }
    
  // Histogram Switchs
  
  void SwitchOnFillAllCellTimeHisto()           { fFillAllCellTimeHisto = kTRUE  ; }
  void SwitchOffFillAllCellTimeHisto()          { fFillAllCellTimeHisto = kFALSE ; }
  
  void SwitchOnFillAllPositionHistogram()       { fFillAllPosHisto  = kTRUE  ; }
  void SwitchOffFillAllPositionHistogram()      { fFillAllPosHisto  = kFALSE ; }
  
  void SwitchOnFillAllPositionHistogram2()      { fFillAllPosHisto2 = kTRUE  ; }
  void SwitchOffFillAllPositionHistogram2()     { fFillAllPosHisto2 = kFALSE ; }
  
  void SwitchOnFillAllTH3Histogram()            { fFillAllTH3       = kTRUE  ; }
  void SwitchOffFillAllTH3Histogram()           { fFillAllTH3       = kFALSE ; }
  
  void SwitchOnFillAllTrackMatchingHistogram()  { fFillAllTMHisto   = kTRUE  ; }
  void SwitchOffFillAllTrackMatchingHistogram() { fFillAllTMHisto   = kFALSE ; }
  
  void SwitchOnFillAllPi0Histogram()            { fFillAllPi0Histo  = kTRUE  ; }
  void SwitchOffFillAllPi0Histogram()           { fFillAllPi0Histo  = kFALSE ; }

  void SwitchOnCorrelation()                    { fCorrelate        = kTRUE  ; }
  void SwitchOffCorrelation()                   { fCorrelate        = kFALSE ; }

  void SwitchOnStudyBadClusters()               { fStudyBadClusters = kTRUE  ; }
  void SwitchOffStudyBadClusters()              { fStudyBadClusters = kFALSE ; }
  
  void SwitchOnStudyClustersAsymmetry()         { fStudyClustersAsymmetry = kTRUE  ; }
  void SwitchOffStudyClustersAsymmetry()        { fStudyClustersAsymmetry = kFALSE ; }

  void SwitchOnStudyWeight()                    { fStudyWeight      = kTRUE  ; }
  void SwitchOffStudyWeight()                   { fStudyWeight      = kFALSE ; }
  
  void SwitchOnStudyExotic()                    { fStudyExotic      = kTRUE  ; }
  void SwitchOffStudyExotic()                   { fStudyExotic      = kFALSE ; }
  
  void SetNECrossCuts(Int_t n)                  { fExoNECrossCuts   = n      ; }
  void SetNDTimeCuts (Int_t n)                  { fExoNDTimeCuts    = n      ; }
  
  void SetExoECrossCuts (Int_t i, Float_t c)    { if (i<fExoNECrossCuts) fExoECrossCuts[i] = c ; }
  void SetExoDTimeCuts  (Int_t i, Float_t c)    { if (i<fExoNDTimeCuts ) fExoDTimeCuts [i] = c ; }
  
 private:
  
  TString  fCalorimeter ;                     // Calorimeter selection
  
  //Switches
  Bool_t   fFillAllCellTimeHisto;             // Fill all cell time histo
  Bool_t   fFillAllPosHisto;                  // Fill all the position related histograms 
  Bool_t   fFillAllPosHisto2;                 // Fill all the position related histograms 2
  Bool_t   fFillAllTH3 ;                      // Fill TH3 histograms
  Bool_t   fFillAllTMHisto ;                  // Fill track matching histograms
  Bool_t   fFillAllPi0Histo ;                 // Fill invariant mass histograms
  Bool_t   fCorrelate   ;                     // Correlate PHOS/EMCAL cells/clusters, also with V0 and track multiplicity
  Bool_t   fStudyBadClusters;                 // Study bad clusters
  Bool_t   fStudyClustersAsymmetry;           // Study asymmetry of clusters
  Bool_t   fStudyExotic;                      // Study the exotic cluster for different cuts
  Bool_t   fStudyWeight;                      // Study the energy weight used in different cluster calculations
  
  // Parameters
  Int_t    fNModules    ;                     // Number of EMCAL/PHOS modules
  Int_t    fNRCU        ;                     // Number of EMCAL/PHOS RCU 
  Int_t    fNMaxCols    ;                     // Number of EMCAL/PHOS rows 
  Int_t    fNMaxRows    ;                     // Number of EMCAL/PHOS columns
  
  //Cuts
  Double_t fTimeCutMin  ;                     // Remove clusters/cells with time smaller than this value, in ns
  Double_t fTimeCutMax  ;                     // Remove clusters/cells with time larger than this value, in ns
  Float_t  fCellAmpMin;                       // amplitude Threshold on calorimeter cells, set at execution time
  Float_t  fEMCALCellAmpMin;                  // amplitude Threshold on emcal cells
  Float_t  fPHOSCellAmpMin ;                  // amplitude Threshold on phos cells
  
  // Exotic studies
  Float_t  fExoNECrossCuts   ;                // Number of ecross cuts
  Float_t  fExoECrossCuts[10];                // List of ecross cuts
  Float_t  fExoNDTimeCuts    ;                // Number of time cuts
  Float_t  fExoDTimeCuts[5]  ;                // List of time cuts
  
  //CaloClusters 
  TH1F *   fhE  ;                             //! E distribution, Reco
  TH1F *   fhPt ;                             //! pT distribution, Reco
  TH1F *   fhPhi;                             //! phi distribution, Reco 
  TH1F *   fhEta;                             //! eta distribution, Reco 
  TH3F *   fhEtaPhiE  ;                       //! eta vs phi vs E, Reco
  TH1F *   fhECharged  ;                      //! E distribution, Reco, matched with track
  TH1F *   fhPtCharged ;                      //! pT distribution, Reco, matched with track
  TH1F *   fhPhiCharged;                      //! phi distribution, Reco, matched with track 
  TH1F *   fhEtaCharged;                      //! eta distribution, Reco, matched with track 
  TH3F *   fhEtaPhiECharged;                  //! eta vs phi vs E, Reco, matched with track 
    
  TH2F *   fhIM;                              //! cluster pairs invariant mass
  TH2F *   fhAsym;                            //! cluster pairs invariant mass	
  
  TH2F *   fhNCellsPerCluster;                //! N cells per cluster vs cluster energy vs eta of cluster	
  TH2F *   fhNCellsPerClusterNoCut;           //! N cells per cluster vs cluster energy vs eta of cluster	

  TH1F *   fhNClusters;                       //! Number of clusters

  TH2F *   fhClusterTimeEnergy;               //! Cluster Time vs Energy 
  TH2F *   fhCellTimeSpreadRespectToCellMax;  //! Difference of the time of cell with maximum dep energy and the rest of cells
  TH1F *   fhCellIdCellLargeTimeSpread;       //! Cells with large time respect to max (diff > 100 ns)
  TH2F *   fhClusterPairDiffTimeE;            //! Pair of clusters time difference vs E
  
  TH2F *   fhClusterMaxCellCloseCellRatio;    //! Ratio between max cell energy and cell energy of the same cluster 
  TH2F *   fhClusterMaxCellCloseCellDiff;     //! Difference between max cell energy and cell energy of the same cluster   
  TH2F *   fhClusterMaxCellDiff;              //! Difference between cluster energy and energy of cell with more energy, good clusters only
  TH2F *   fhClusterMaxCellDiffNoCut;         //! Difference between cluster energy and energy of cell with more energy, no bad cluster rejection
  
  TH2F *   fhClusterMaxCellDiffAverageTime;      //! Difference between cluster average time and time of cell with more energy
  TH2F *   fhClusterMaxCellDiffWeightedTime;     //! Difference between cluster weighted time and time of cell with more energy
  TH2F *   fhClusterMaxCellECross;               //! 1 - Energy in cross around max energy cell / max energy cell vs cluster energy, good clusters
  
  TH2F *   fhLambda0;                         //! cluster Lambda0    vs Energy
  TH2F *   fhLambda1;                         //! cluster Lambda1    vs Energy
  TH2F *   fhDispersion;                      //! cluster Dispersion vs Energy
  
  // Bad clusters histograms
  TH1F *   fhBadClusterEnergy;                //! energy of bad cluster
  TH2F *   fhBadClusterTimeEnergy;            //! Time Max cell of bad cluster
  TH2F *   fhBadClusterPairDiffTimeE;         //! Pair of clusters time difference vs E, bad cluster
  TH2F *   fhBadCellTimeSpreadRespectToCellMax; //! Difference of the time of cell with maximum dep energy and the rest of cells for bad clusters
  
  TH2F *   fhBadClusterMaxCellCloseCellRatio; //! Ratio between max cell energy and cell energy of the same cluster for bad clusters 
  TH2F *   fhBadClusterMaxCellCloseCellDiff ; //! Difference between max cell energy and cell energy of the same cluster for bad clusters 
  TH2F *   fhBadClusterMaxCellDiff;           //! Difference between cluster energy and energy of cell with more energy
  
  TH2F *   fhBadClusterMaxCellDiffAverageTime;      //! Difference between cluster average time and time of cell with more energy
  TH2F *   fhBadClusterMaxCellDiffWeightedTime;     //! Difference between cluster weighted time and time of cell with more energy
  TH2F *   fhBadClusterMaxCellECross;               //! 1 - Energy in cross around max energy cell / max energy cell vs cluster energy, bad clusters

  // Cluster cell size
  TH2F *   fhDeltaIEtaDeltaIPhiE0[2];         //! Difference between max cell index and farthest cell, eta vs phi, E < 2 GeV, with and without matching; 
  TH2F *   fhDeltaIEtaDeltaIPhiE2[2];         //! Difference between max cell index and farthest cell, eta vs phi, 2 < E < 6 GeV, with and without matching; 
  TH2F *   fhDeltaIEtaDeltaIPhiE6[2];         //! Difference between max cell index and farthest cell, eta vs phi, E > 6 GeV, with and without matching; 
  TH2F *   fhDeltaIA[2];                      //! Cluster "asymmetry" in cell terms vs E, with and without matching
  TH2F *   fhDeltaIAL0[2];                    //! Cluster "asymmetry" in cell units vs Lambda0    for E > 0.5 GeV, n cells in cluster > 3, with and without matching
  TH2F *   fhDeltaIAL1[2];                    //! Cluster "asymmetry" in cell units vs Lambda1    for E > 0.5 GeV, n cells in cluster > 3, with and without matching
  TH2F *   fhDeltaIANCells[2] ;               //! Cluster "asymmetry" in cell units vs number of cells in cluster for E > 0.5, with and without matching
  TH2F *   fhDeltaIAMC[4];                    //! Cluster "asymmetry" in cell terms vs E, from MC photon, electron, conversion or hadron.
  TH2F *   fhBadClusterDeltaIEtaDeltaIPhiE0;  //! Difference between max cell index and farthest cell, eta vs phi, E < 2 GeV, with and without matching; bad clusters. 
  TH2F *   fhBadClusterDeltaIEtaDeltaIPhiE2;  //! Difference between max cell index and farthest cell, eta vs phi, 2 < E < 6 GeV, with and without matching; bad clusters.
  TH2F *   fhBadClusterDeltaIEtaDeltaIPhiE6;  //! Difference between max cell index and farthest cell, eta vs phi, E > 6 GeV, with and without matching; bad clusters.
  TH2F *   fhBadClusterDeltaIA;               //! Cluster "asymmetry" in cell terms vs E, with and without matching; bad clusters.
  
  //Cluster/cell Position
  TH2F *   fhRNCells ;                        //! R=sqrt(x^2+y^2) (cm) cluster distribution vs N cells in cluster
  TH2F *   fhXNCells ;                        //! X (cm) cluster distribution vs N cells in cluster
  TH2F *   fhYNCells ;                        //! Y (cm) cluster distribution vs N cells in cluster
  TH2F *   fhZNCells ;                        //! Z (cm) cluster distribution vs N cells in cluster
	
  TH2F *   fhRE ;                             //! R=sqrt(x^2+y^2) (cm) cluster distribution vs cluster energy
  TH2F *   fhXE ;                             //! X (cm) cluster distribution vs cluster energy
  TH2F *   fhYE ;                             //! Y (cm) cluster distribution vs cluster energy
  TH2F *   fhZE ;                             //! Z (cm) cluster distribution vs cluster energy
  TH3F *   fhXYZ;                             //! cluster X vs Y vs Z (cm)
	
  TH2F *   fhRCellE ;                         //! R=sqrt(x^2+y^2) (cm) cell distribution vs cell energy
  TH2F *   fhXCellE ;                         //! X (cm) cell distribution vs cell energy
  TH2F *   fhYCellE ;                         //! Y (cm) cell distribution vs cell energy
  TH2F *   fhZCellE ;                         //! Z (cm) cell distribution vs cell energy
  TH3F *   fhXYZCell;                         //! cell X vs Y vs Z (cm)
  
  TH2F *   fhDeltaCellClusterRNCells ;        //! R cluster - R cell distribution (cm) vs N cells in cluster
  TH2F *   fhDeltaCellClusterXNCells ;        //! X cluster - X cell distribution (cm) vs N cells in cluster
  TH2F *   fhDeltaCellClusterYNCells ;        //! Y cluster - Y cell distribution (cm) vs N cells in cluster
  TH2F *   fhDeltaCellClusterZNCells ;        //! Z cluster - Z cell distribution (cm) vs N cells in cluster
	
  TH2F *   fhDeltaCellClusterRE ;             //! R cluster - R cell distribution (cm) vs cluster energy
  TH2F *   fhDeltaCellClusterXE ;             //! X cluster - X cell distribution (cm) vs cluster energy
  TH2F *   fhDeltaCellClusterYE ;             //! Y cluster - Y cell distribution (cm) vs cluster energy
  TH2F *   fhDeltaCellClusterZE ;             //! Z cluster - Z cell distribution (cm) vs cluster energy
	
  //Calo Cells
  TH1F *   fhNCells;                          //! Number of towers/crystals with signal
  TH1F *   fhNCellsCutAmpMin;                 //! Number of towers/crystals with signal, with min amplitude
  TH1F *   fhAmplitude;                       //! Amplitude measured in towers/crystals
  TH2F *   fhAmpId;                           //! Amplitude measured in towers/crystals vs id of tower.	
  TH3F *   fhEtaPhiAmp;                       //! eta vs phi vs amplitude, cells
   
  TH1F *   fhTime;                            //! Time measured in towers/crystals
  TH2F *   fhTimeVz;                          //! Time measured in towers/crystals vs vertex z component, for E > 0.5
  TH2F *   fhTimeId;                          //! Time vs Absolute cell Id
  TH2F *   fhTimeAmp;                         //! Time vs Amplitude 
  
  TH2F *   fhCellECross;                      //! 1 - Energy in cross around cell /  cell energy 
  
  //Calorimeters Correlation
  TH2F *   fhCaloCorrNClusters;               //! EMCAL vs PHOS, number of clusters	
  TH2F *   fhCaloCorrEClusters;               //! EMCAL vs PHOS, total measured cluster energy
  TH2F *   fhCaloCorrNCells;                  //! EMCAL vs PHOS, number of cells
  TH2F *   fhCaloCorrECells;                  //! EMCAL vs PHOS,  total measured cell energy
	
  //V0 Correlation
  TH2F *   fhCaloV0SCorrNClusters;            //! Calo vs V0 signal , number of clusters	
  TH2F *   fhCaloV0SCorrEClusters;            //! Calo vs V0 signal, total measured cluster energy
  TH2F *   fhCaloV0SCorrNCells;               //! Calo vs V0 signal, number of cells
  TH2F *   fhCaloV0SCorrECells;               //! Calo vs V0 signal,  total measured cell energy
  TH2F *   fhCaloV0MCorrNClusters;            //! Calo vs V0 multiplicity , number of clusters	
  TH2F *   fhCaloV0MCorrEClusters;            //! Calo vs V0 multiplicity, total measured cluster energy
  TH2F *   fhCaloV0MCorrNCells;               //! Calo vs V0 multiplicity, number of cells
  TH2F *   fhCaloV0MCorrECells;               //! Calo vs V0 multiplicity,  total measured cell energy
  
  //Track Correlation
  TH2F *   fhCaloTrackMCorrNClusters;         //! Calo vs Track Multiplicity, number of clusters	
  TH2F *   fhCaloTrackMCorrEClusters;         //! Calo vs Track Multiplicity, total measured cluster energy
  TH2F *   fhCaloTrackMCorrNCells;            //! Calo vs V0 Track Multiplicity, number of cells
  TH2F *   fhCaloTrackMCorrECells;            //! Calo vs V0 Track Multipliticy,  total measured cell energy
  
  // Centrality
  TH2F *   fhCaloCenNClusters;               //! Calo vs centrality, number of clusters
  TH2F *   fhCaloCenEClusters;               //! Calo vs centrality, total measured cluster energy
  TH2F *   fhCaloCenNCells;                  //! Calo vs centrality, number of cells
  TH2F *   fhCaloCenECells;                  //! Calo vs centrality,  total measured cell energy

  // Event plane
  TH2F *   fhCaloEvPNClusters;               //! Calo vs event plane angle, number of clusters
  TH2F *   fhCaloEvPEClusters;               //! Calo vs event plane angle, total measured cluster energy
  TH2F *   fhCaloEvPNCells;                  //! Calo vs event plane angle, number of cells
  TH2F *   fhCaloEvPECells;                  //! Calo vs event plane angle,  total measured cell energy
  
  //Module histograms
  TH2F *   fhEMod  ;                          //! cluster E distribution for different module, Reco
  TH2F *   fhAmpMod ;                         //! cell amplitude distribution for different module, Reco
  TH2F *   fhTimeMod ;                        //! cell time distribution for different module, Reco
  TH2F *   fhNClustersMod ;                   //! Number of clusters for different module, Reco
  TH2F *   fhNCellsMod ;                      //! Number of towers/crystals with signal different module, Reco
  TH2F **  fhNCellsPerClusterMod ;            //! N cells per clusters different module, Reco
  TH2F **  fhNCellsPerClusterModNoCut ;       //! N cells per clusters different module, Reco, No cut
  TH2F *   fhGridCells ;                      //! Cells ordered in column/row for different module, Reco
  TH2F *   fhGridCellsE ;                     //! Cells ordered in column/row for different module, weighted with energy, Reco
  TH2F *   fhGridCellsTime ;                  //! Cells ordered in column/row for different module, weighted with time, Reco
  TH2F **  fhTimeAmpPerRCU;                   //! Time vs Amplitude measured in towers/crystals different RCU
  TH2F **  fhIMMod;                           //! cluster pairs invariant mass, different module,
	
  // Weight studies
  
  TH2F *   fhECellClusterRatio;               //! e cell / e cluster vs e cluster
  TH2F *   fhECellClusterLogRatio;            //! log (e cell / e cluster)  vs e cluster
  TH2F *   fhEMaxCellClusterRatio;            //! e max cell / e cluster vs e cluster
  TH2F *   fhEMaxCellClusterLogRatio;         //! log (e max cell / e cluster) vs e cluster
  
  TH2F *   fhLambda0ForW0[14];                //! L0 for 7 defined w0= 3, 3.5 ... 6
  //TH2F * fhLambda1ForW0[7];                  //! L1 for 7 defined w0= 3, 3.5 ... 6

  TH2F *   fhLambda0ForW0MC[14][5];            //! L0 for 7 defined w0= 3, 3.5 ... 6, depending on the particle of origin
  //TH2F * fhLambda1ForW0MC[7][5];              //! L1 for 7 defined w0= 3, 3.5 ... 6, depending on the particle of origin
  
  TH2F *   fhECellTotalRatio;                 //! e cell / e total vs e total
  TH2F *   fhECellTotalLogRatio;              //! log (e cell / e total)  vs e total
  TH2F **  fhECellTotalRatioMod;              //! e cell / e total vs e total, per SM
  TH2F **  fhECellTotalLogRatioMod;           //! log (e cell / e total)  vs e total, per SM

  // Exotic studies
  
  TH2F *   fhExoNCell  [10][5] ;               //! Number of cells per cluster for different cuts
  TH2F *   fhExoL0     [10][5] ;               //! Long shower shape axis for exotic
  TH2F *   fhExoL1     [10][5] ;               //! Short shower shape axis for exotic
  TH2F *   fhExoECross [10][5] ;               //! E cross for max cell in cluster, for different cuts
  TH2F *   fhExoTime   [10][5] ;               //! Time of exotic cluster, for different cuts
  TH2F *   fhExoDTime  [10]    ;               //! Difference in time between cell with max energy and rest of cells for exotic
  TH2F *   fhExoL0NCell[10][5] ;               //! Lambda0 vs n cells in cluster for several E cross cuts and cluster with E > 5
  TH2F *   fhExoL0ECross      ;                //! Lambda0 vs E cross fraction for clusters with E > 5 GeV
  TH2F *   fhExoL1NCell[10][5] ;               //! Lambda1 vs n cells in cluster for several E cross cuts and cluster with E > 5
  TH2F *   fhExoL1ECross      ;                //! Lambda1 vs E cross fraction for clusters with E > 5 GeV
  
  //Pure MC

  enum mcTypes {kmcPhoton = 0, kmcPi0 = 1, kmcEta = 2, kmcElectron = 3, kmcNeHadron = 4, kmcChHadron = 5 };
  
  TH2F *   fhRecoMCE[6][2]  ;                 //! E   generated particle vs reconstructed E
  TH2F *   fhRecoMCPhi[6][2] ;                //! phi generated particle vs reconstructed phi
  TH2F *   fhRecoMCEta[6][2] ;                //! eta generated particle vs reconstructed Eta
  TH2F *   fhRecoMCDeltaE[6][2]  ;            //! Gen-Reco E    generated particle vs reconstructed E
  TH2F *   fhRecoMCRatioE[6][2]  ;            //! Reco/Gen E    generated particle vs reconstructed E
  TH2F *   fhRecoMCDeltaPhi[6][2];            //! Gen-Reco phi  generated particle vs reconstructed E
  TH2F *   fhRecoMCDeltaEta[6][2];            //! Gen-Reco eta  generated particle vs reconstructed E
  
  TH1F *   fhGenMCE[4]     ;                  //! pt of primary particle
  TH2F *   fhGenMCEtaPhi[4] ;                 //! eta vs phi of primary particle
  TH1F *   fhGenMCAccE[4]     ;               //! pt of primary particle, in acceptance
  TH2F *   fhGenMCAccEtaPhi[4] ;              //! eta vs phi of primary particle, in acceptance
  
  TH2F *   fhEMVxyz    ;                      //! Electromagnetic particle production vertex
  TH2F *   fhEMR       ;                      //! Electromagnetic distance to vertex vs rec energy  
  TH2F *   fhHaVxyz    ;                      //! Hadron production vertex
  TH2F *   fhHaR       ;                      //! Hadron distance to vertex vs rec energy  
	
  //Histograms for MC track-matching
  TH2F *   fh1EOverP;                         //! p/E for track-cluster matches
  TH2F *   fh2dR;                             //! distance between projected track and cluster (eta-phi units)
  TH2F *   fh2EledEdx;                        //! dE/dx vs. momentum for electron candidates
  TH2F *   fh2MatchdEdx;                      //! dE/dx vs. momentum for all matches
	
  TH2F *   fhMCEle1EOverP;                    //! p/E for track-cluster matches, MC electrons
  TH1F *   fhMCEle1dR;                        //! distance between projected track and cluster, MC electrons
  TH2F *   fhMCEle2MatchdEdx;                 //! dE/dx vs. momentum for all matches, MC electrons	
	
  TH2F *   fhMCChHad1EOverP;                  //! p/E for track-cluster matches, MC charged hadrons
  TH1F *   fhMCChHad1dR;                      //! distance between projected track and cluster, MC charged hadrons
  TH2F *   fhMCChHad2MatchdEdx;               //! dE/dx vs. momentum for all matches, MC charged
	
  TH2F *   fhMCNeutral1EOverP;                //! p/E for track-cluster matches, MC neutral
  TH1F *   fhMCNeutral1dR;                    //! distance between projected track and cluster, MC neutral
  TH2F *   fhMCNeutral2MatchdEdx;             //! dE/dx vs. momentum for all matches, MC neutral	
	
  TH2F *   fh1EOverPR02;                      //! p/E for track-cluster matches, dR < 0.2
  TH2F *   fhMCEle1EOverPR02;                 //! p/E for track-cluster matches, dR < 0.2, MC electrons
  TH2F *   fhMCChHad1EOverPR02;               //! p/E for track-cluster matches, dR < 0.2, MC charged hadrons
  TH2F *   fhMCNeutral1EOverPR02;             //! p/E for track-cluster matches, dR < 0.2, MC neutral

  TH2F *   fh1EleEOverP;                      //! p/E for track-cluster matches, dR < 0.2, 60 < dEdx < 100
  TH2F *   fhMCEle1EleEOverP;                 //! p/E for track-cluster matches, dR < 0.2, 60 < dEdx < 100, MC electrons
  TH2F *   fhMCChHad1EleEOverP;               //! p/E for track-cluster matches, dR < 0.2, 60 < dEdx < 100, MC charged hadrons
  TH2F *   fhMCNeutral1EleEOverP;             //! p/E for track-cluster matches, dR < 0.2, 60 < dEdx < 100, MC neutral

  TH2F *   fhTrackMatchedDEta    ;            //! Eta distance between track and cluster vs cluster E, after and before photon cuts
  TH2F *   fhTrackMatchedDPhi    ;            //! Phi distance between track and cluster vs cluster E, after and before photon cuts
  TH2F *   fhTrackMatchedDEtaDPhi;            //! Eta vs Phi distance between track and cluster, E cluster > 0.5 GeV, after and before
  
  TH2F *   fhTrackMatchedDEtaPos;             //! Eta distance between track and cluster vs cluster E, after and before photon cuts
  TH2F *   fhTrackMatchedDPhiPos;             //! Phi distance between track and cluster vs cluster E, after and before photon cuts
  TH2F *   fhTrackMatchedDEtaDPhiPos;         //! Eta vs Phi distance between track and cluster, E cluster > 0.5 GeV, after and before

  AliAnaCalorimeterQA & operator = (const AliAnaCalorimeterQA & qa) ;//cpy assignment
  AliAnaCalorimeterQA(              const AliAnaCalorimeterQA & qa) ; // cpy ctor
  
  ClassDef(AliAnaCalorimeterQA,28)
} ;


#endif //ALIANACALORIMETERQA_H



